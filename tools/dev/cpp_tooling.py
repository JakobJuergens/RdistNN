#!/usr/bin/env python3
"""Isolated native builds, real compile databases and focused native checks.

Python 3.9+. Linux developer tooling; never edits production build configuration.
Every build is fresh, but compiler objects use a bounded private ccache. Test
results are never cached. See CONTRIBUTING.md for host sanitizer limitations.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
IS_R = (ROOT / 'DESCRIPTION').is_file()
OUTPUT = ROOT / ('archives' if IS_R else 'build') / 'cpp-tooling'
FILTERS = ('native-index-boundaries', 'native-threads', 'native-plqr-threshold-recovery',
           'errors-fastpaths', 'score-plqr-regressions', 'plqr-comparison-range',
           'seeded-rng-scope', 'compiled-fallback-parity', 'compensated-accumulation',
           'distance-scaled-arithmetic', 'distance-scaled-underflow',
           'distance-range-stability', 'block-numerical-regression')
TESTS = ('core', 'distance', 'distance_plugin', 'distance_rounding', 'blocks', 'loo_plugin')
SAN_FLAGS = ['-fsanitize=address,undefined', '-fno-sanitize-recover=all',
             '-fno-omit-frame-pointer', '-g', '-O1']


class ToolError(RuntimeError):
    """A requested check could not complete successfully."""


def require_tool(name):
    path = shutil.which(name)
    if not path:
        raise ToolError('Required tool {} not found on PATH; install it before retrying.'.format(name))
    return str(Path(path).absolute())


def run(command, cwd, env, logfile):
    reject_symlink_tree(cwd)
    reject_symlink_path(logfile)
    merged = os.environ.copy()
    merged.update(env)
    with logfile.open('a', encoding='utf-8') as log:
        log.write('$ ' + shlex.join([str(x) for x in command]) + '\n')
        log.flush()
        result = subprocess.run(command, cwd=cwd, env=merged, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise ToolError('Command failed (exit {}); see {}\n{}'.format(
            result.returncode, logfile, '\n'.join(safe_read_text(logfile).splitlines()[-12:])))


def capture(command):
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise ToolError('{} failed: {}'.format(shlex.join(command), result.stderr.strip()))
    return result.stdout.strip()


def reject_symlink_path(path):
    for component in (path, *path.parents):
        if component.is_symlink():
            raise ToolError('Unsafe symlink path: {}'.format(component))


def reject_symlink_tree(path):
    """Reject redirected output paths and every existing descendant before use."""
    reject_symlink_path(path)
    if path.is_dir():
        for directory, folders, files in os.walk(path, followlinks=False):
            for name in folders + files:
                candidate = Path(directory) / name
                if candidate.is_symlink():
                    raise ToolError('Unsafe symlink path: {}'.format(candidate))


def safe_mkdir(path):
    reject_symlink_tree(path)
    path.mkdir(parents=True, exist_ok=True)


def safe_write_text(path, contents):
    reject_symlink_path(path)
    path.write_text(contents, encoding='utf-8')


def safe_read_bytes(path):
    reject_symlink_path(path)
    return path.read_bytes()


def safe_read_text(path):
    reject_symlink_path(path)
    return path.read_text(encoding='utf-8', errors='replace')


def snapshot(root, target):
    """Copy current working inputs, including uncommitted changes, without outputs."""
    reject_symlink_path(root)
    reject_symlink_tree(target)
    root, target = root.absolute(), target.absolute()
    if target == root or target in root.parents:
        raise ToolError('Unsafe snapshot destination would replace source: {}'.format(target))
    if root in target.parents and target.relative_to(root).parts[0] not in ('archives', 'build'):
        raise ToolError('Snapshot inside source must use the excluded archives/build tree.')
    if target.exists():
        shutil.rmtree(target)
    excluded = {'.git', 'archives', 'build', '_dist', '.Rproj.user', '__pycache__', '.pytest_cache',
                '.ruff_cache', 'node_modules', 'revdep', 'doc', 'Meta',
                '.superpowers', '.codex', '.agents', '.claude', 'GH_logs', '.venv', 'venv'}
    if IS_R:
        excluded.add('docs')  # Generated pkgdown output; Stata docs are maintained source.
    suffixes = {'.o', '.so', '.dll', '.dylib', '.plugin', '.gcda', '.gcno', '.pyc'}

    def ignore(directory, names):
        ignored = [name for name in names if name in excluded or Path(name).suffix in suffixes]
        for name in set(names) - set(ignored):
            if (Path(directory) / name).is_symlink():
                raise ToolError('Source symlink cannot enter an isolated snapshot: {}'.format(Path(directory) / name))
        return ignored

    shutil.copytree(root, target, ignore=ignore)


def validate_database(entries, owned, compiler=None):
    seen = set()
    for entry in entries:
        directory = Path(entry['directory'])
        source = Path(entry['file'])
        if not source.is_absolute():
            source = directory / source
        if not directory.is_dir() or not source.is_file():
            raise ToolError('Compilation database path does not exist: {}'.format(source))
        args = entry.get('arguments') or shlex.split(entry.get('command', ''))
        if not args or not any('++' in token or 'clang' in token for token in args[:2]):
            raise ToolError('Compilation database lacks a C++ compiler: {}'.format(source))
        if compiler is not None and Path(args[0]).name != Path(compiler).name:
            raise ToolError('Compilation database compiler differs from selected profile: {}'.format(source))
        if not any(token in ('-std=c++17', '-std=gnu++17') for token in args):
            raise ToolError('Compilation database lacks the C++17 contract: {}'.format(source))
        if not any(token.startswith('-I') or token.startswith('-isystem') for token in args):
            raise ToolError('Compilation database lacks include flags: {}'.format(source))
        seen.add(source.resolve())
    missing = sorted(str(path) for path in owned if path.resolve() not in seen)
    if missing:
        raise ToolError('Compilation database missing owned translation units: ' + ', '.join(missing))
    if not owned:
        raise ToolError('No owned translation units found; refresh the source with compdb.')


def config(profile, sanitized=False):
    compiler = require_tool('g++' if profile == 'gcc' else 'clang++-18')
    tools = {'compiler': compiler, 'ccache': require_tool('ccache'),
             'make': require_tool('make')}
    if IS_R:
        tools.update(R=require_tool('R'), Rscript=require_tool('Rscript'))
    if sanitized:
        tools['symbolizer'] = require_tool('llvm-symbolizer-18')
        runtime = capture([compiler, '-print-file-name=libclang_rt.asan-x86_64.so'])
        if not Path(runtime).is_file():
            raise ToolError('Clang 18 ASan shared runtime is unavailable: ' + runtime)
        tools['asan_runtime'] = runtime
    return tools


def environment(work, tools, profile):
    version = capture([tools['compiler'], '--version']).splitlines()[0]
    r_version = capture([tools['R'], '--version']).splitlines()[0] if IS_R else 'standalone'
    cache_key = hashlib.sha256((version + r_version + profile).encode()).hexdigest()[:16]
    cache = OUTPUT / 'cache' / cache_key
    safe_mkdir(cache)
    cache_config = cache / 'ccache.conf'
    safe_write_text(cache_config, 'max_size = 500M\ncompiler_check = content\n')
    return {'CCACHE_DIR': str(cache), 'CCACHE_CONFIGPATH': str(cache_config),
            'CCACHE_MAXSIZE': '500M', 'CCACHE_COMPILERCHECK': 'content', 'CCACHE_SLOPPINESS': '',
            'MAKEFLAGS': '-j{}'.format(work['jobs']), 'OMP_THREAD_LIMIT': '2',
            'OMP_NUM_THREADS': '2', 'OPENBLAS_NUM_THREADS': '1',
            'R_MAKEVARS_SITE': '/dev/null', 'R_ENVIRON_USER': '/dev/null',
            'R_PROFILE_USER': '/dev/null', 'DEBUGINFOD_URLS': ''}, version, r_version


def openmp_probe(work, tools, profile, env):
    if profile == 'clang-serial':
        return ''
    source = work['base'] / 'openmp-probe.cpp'
    safe_write_text(source, '#include <omp.h>\nint main(){return omp_get_max_threads()<1;}\n')
    flags = '-fopenmp' if profile == 'gcc' else '-fopenmp=libomp'
    run([tools['compiler'], flags, str(source), '-o', str(work['base'] / 'openmp-probe')],
        work['base'], env, work['log'])
    run([str(work['base'] / 'openmp-probe')], work['base'], env, work['log'])
    return flags


def r_makevars(work, tools, omp, sanitized):
    # Override both compile and link owners; strip only incompatible GNU LTO flags.
    def flags(name):
        values = shlex.split(capture([tools['R'], 'CMD', 'config', name]))
        if 'clang' in Path(tools['compiler']).name:
            values = [value for value in values if not value.startswith('-flto') and value != '-ffat-lto-objects']
        return values
    common = flags('CXX17FLAGS')
    if sanitized:
        common = [value for value in common if not value.startswith('-O')] + SAN_FLAGS
    linker = flags('LDFLAGS') + (SAN_FLAGS if sanitized else [])
    compiler = shlex.join([tools['ccache'], tools['compiler']])
    text = [
        'CXX = ' + compiler, 'CXX17 = ' + compiler,
        'CXX17STD = -std=gnu++17', 'CXXFLAGS = ' + shlex.join(common),
        'CXX17FLAGS = ' + shlex.join(common),
        'SHLIB_CXXLD = ' + compiler, 'SHLIB_CXX17LD = ' + compiler,
        'SHLIB_CXXLDFLAGS = -shared', 'SHLIB_CXX17LDFLAGS = -shared',
        'LDFLAGS = ' + shlex.join(linker), 'SHLIB_OPENMP_CXXFLAGS = ' + omp,
    ]
    path = work['base'] / 'Makevars'
    safe_write_text(path, '\n'.join(text) + '\n')
    return path


def test_r(work, tools, env, sanitized):
    script = work['base'] / 'installed-tests.R'
    safe_write_text(script, '''args <- commandArgs(TRUE)
.libPaths(c(args[[1]], .libPaths()))
library(RdistNN)
stopifnot(normalizePath(find.package("RdistNN")) == normalizePath(file.path(args[[1]], "RdistNN")))
options(RdistNN.use_compiled = TRUE, RdistNN.num_threads = 2L)
required <- c("rustat_dist_l1", "rustat_dist_l2", "rustat_thread_info",
  "rustat_subset_row_sums_g", "rustat_exact_rank_table", "rustat_compensated_sum",
  "rustat_compensated_crossprod", "rustat_block_rank_loo", "rustat_jackknife_fast_rank",
  "rustat_jackknife_rank_full", "rustat_plqr_ctx_create", "rustat_plqr_g_all")
for (name in required) stopifnot(is.function(RdistNN:::.rdistnn_compiled_fun(name)))
cat("Installed package:", find.package("RdistNN"), "\\n")
cat("Loaded native DLL:", getLoadedDLLs()[["RdistNN"]][["path"]], "\\n")
stopifnot(normalizePath(getLoadedDLLs()[["RdistNN"]][["path"]]) ==
  normalizePath(file.path(args[[1]], "RdistNN", "libs", paste0("RdistNN", .Platform$dynlib.ext))))
print(RdistNN:::rustat_thread_info())
stopifnot(length(RdistNN:::rustat_dist_l2(matrix(1:6, 3), c(0,0))) == 3L)
result <- testthat::test_package("RdistNN", filter=args[[2]], reporter="summary", stop_on_failure=TRUE)
if (!length(result)) stop("No installed native tests ran")
flat <- unlist(lapply(result, function(x) x$results), recursive=FALSE)
if (any(vapply(flat, inherits, logical(1), "expectation_skip"))) stop("Native profile may not skip tests")
cat("Installed native proof and focused tests passed\\n")
''')
    test_env = dict(env)
    if sanitized:
        test_env.update(LD_PRELOAD=tools['asan_runtime'],
                        ASAN_OPTIONS='detect_leaks=0:halt_on_error=1:strict_string_checks=1',
                        UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1',
                        ASAN_SYMBOLIZER_PATH=tools['symbolizer'])
    run([tools['Rscript'], '--vanilla', str(script), str(work['base'] / 'library'),
         '^(' + '|'.join(FILTERS) + ')$'], work['source'], test_env, work['log'])


def build_r(work, tools, env, omp, sanitized, bear):
    env['R_MAKEVARS_USER'] = str(r_makevars(work, tools, omp, sanitized))
    library = work['base'] / 'library'
    safe_mkdir(library)
    command = [tools['R'], 'CMD', 'INSTALL', '--preclean', '--install-tests', '--no-multiarch',
               '--no-test-load', '--library=' + str(library), str(work['source'])]
    if bear:
        env['CCACHE_DISABLE'] = '1'
        # Bear's local gRPC collector must bypass inherited network proxies.
        env['grpc_proxy'] = ''
        env['no_proxy'] = 'localhost,127.0.0.1'
        command = [require_tool('bear'), '--output', str(work['database']), '--'] + command
    run(command, work['source'], env, work['log'])
    test_r(work, tools, env, sanitized)


def build_stata(work, tools, env, sanitized, bear):
    src = work['source'] / 'src/plugin'
    objects = work['base'] / 'objects'
    safe_mkdir(objects)
    flags = ['-std=c++17', '-O3', '-Wall', '-Wextra', '-fPIC', '-DSYSTEM=OPUNIX', '-I' + str(src)]
    if sanitized:
        flags = [flag for flag in flags if not flag.startswith('-O')] + SAN_FLAGS
    compiler = shlex.join([tools['ccache'], tools['compiler']])
    lines = ['CXX := ' + compiler, 'FLAGS := ' + shlex.join(flags), '.PHONY: all',
             'all: ' + ' '.join(['rdistnn.plugin'] + ['test_' + name for name in TESTS])]
    sources = ['stplugin.c', 'rdistnn_plugin.cpp'] + ['test_rdistnn_' + name + '.cpp' for name in TESTS]
    for name in sources:
        obj = name + '.o'
        # Release g++ compiles stplugin.c as C++; preserve that explicitly for Clang.
        lines += [obj + ':', '\t$(CXX) $(FLAGS) -x c++ -c ' + shlex.quote(str(src / name)) + ' -o $@']
    lines += ['rdistnn.plugin: stplugin.c.o rdistnn_plugin.cpp.o',
              '\t$(CXX) $(FLAGS) -shared $^ -o $@']
    for name in TESTS:
        lines += ['test_' + name + ': test_rdistnn_' + name + '.cpp.o', '\t$(CXX) $(FLAGS) $^ -o $@']
    makefile = work['base'] / 'developer.mk'
    safe_write_text(makefile, '\n'.join(lines) + '\n')
    command = [tools['make'], '-B', '-j' + str(work['jobs']), '-f', str(makefile)]
    if bear:
        env['CCACHE_DISABLE'] = '1'
        # Bear's local gRPC collector must bypass inherited network proxies.
        env['grpc_proxy'] = ''
        env['no_proxy'] = 'localhost,127.0.0.1'
        command = [require_tool('bear'), '--output', str(work['database']), '--'] + command
    run(command, objects, env, work['log'])
    test_env = dict(env)
    if sanitized:
        test_env.update(ASAN_OPTIONS='detect_leaks=1:halt_on_error=1',
                        UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1',
                        ASAN_SYMBOLIZER_PATH=tools['symbolizer'])
    for name in TESTS:
        command = [str(objects / ('test_' + name))]
        if name in ('distance', 'distance_rounding'):
            command += ['2', '2', '0.5790338377234229', '1.4758492921057114']
        run(command, objects, test_env, work['log'])


def owned_sources(source):
    if IS_R:
        return sorted(path for path in (source / 'src').glob('*.cpp') if path.name != 'RcppExports.cpp')
    return sorted((source / 'src/plugin').glob('*.cpp'))


def captured_sources(source):
    if IS_R:
        return sorted((source / 'src').glob('*.cpp'))
    return owned_sources(source) + [source / 'src/plugin/stplugin.c']


def sanitizer_smoke(work, tools, env):
    source = work['base'] / 'sanitizer-negative.cpp'
    safe_write_text(source, 'int main(int argc,char**){ int* p=new int[1]; p[argc]=42; volatile int x=p[argc]; delete[] p; return x; }\n')
    binary = work['base'] / 'sanitizer-negative'
    run([tools['compiler'], *SAN_FLAGS, str(source), '-o', str(binary)], work['base'], env, work['log'])
    # Unique output prevents old diagnostics from certifying this invocation.
    log = work['base'] / ('sanitizer-negative-{}.log'.format(time.time_ns()))
    reject_symlink_path(log)
    # Exclusive creation also refuses an unexpected pre-existing filename.
    with log.open('x', encoding='utf-8'):
        pass
    try:
        run([str(binary)], work['base'], dict(env, ASAN_OPTIONS='detect_leaks=1:halt_on_error=1', ASAN_SYMBOLIZER_PATH=tools['symbolizer']), log)
    except ToolError:
        current_diagnostics = safe_read_text(log)
        if 'ERROR: AddressSanitizer: heap-buffer-overflow' in current_diagnostics:
            print('Negative sanitizer fixture detected; nonzero status propagated.')
            return
        raise
    raise ToolError('Negative sanitizer fixture incorrectly succeeded')


def file_sha256(path):
    return hashlib.sha256(safe_read_bytes(path)).hexdigest()


def input_fingerprint(root):
    paths = [root / '.clang-tidy', root / 'tools/dev/cpp_tooling.py']
    if IS_R:
        paths.append(root / 'DESCRIPTION')
    native_roots = ['src', 'inst/include'] if IS_R else ['src/plugin']
    for directory in native_roots:
        paths.extend(path for path in (root / directory).rglob('*')
                     if path.is_file() and (path.suffix in ('.c', '.cpp', '.h', '.hpp')
                                            or path.name.startswith(('Makevars', 'Makefile'))))
    fingerprints = {}
    for path in sorted(set(paths)):
        reject_symlink_path(path)
        if path.is_file():
            fingerprints[str(path.relative_to(root))] = file_sha256(path)
    return fingerprints


def verify_inputs(fingerprint, source, command):
    """Use the same native/config identity for live and retained inputs."""
    if (not isinstance(fingerprint, dict) or not fingerprint
            or fingerprint != input_fingerprint(ROOT)
            or fingerprint != input_fingerprint(source)):
        raise ToolError('Live/retained native or configuration inputs changed; rerun {} with the same --profile after source edits settle.'.format(command))


def verify_capture(work, profile, compiler_version, r_version, compiler, receipt_sha256=None):
    """Only a complete successful capture can authorize the retained analysis inputs."""
    try:
        receipt = work['base'] / 'provenance.json'
        receipt_contents = safe_read_bytes(receipt)
        receipt_digest = hashlib.sha256(receipt_contents).hexdigest()
        previous = json.loads(receipt_contents)
        if (not isinstance(previous, dict) or previous.get('command') != 'compdb'
                or previous.get('status') != 'passed' or previous.get('profile') != profile):
            raise ToolError('Missing or incomplete successful compdb capture proof.')
        if previous.get('compiler') != compiler_version or previous.get('R') != r_version:
            raise ToolError('Captured compiler/R identity changed.')
        if receipt_sha256 is not None and receipt_sha256 != receipt_digest:
            raise ToolError('Capture proof changed during analysis.')
        fingerprint = previous.get('fingerprint')
        verify_inputs(fingerprint, work['source'], 'compdb')
        sources = previous.get('sources')
        if not isinstance(sources, dict) or any(sources.get(name) != digest for name, digest in fingerprint.items()):
            raise ToolError('Captured source hashes are missing or inconsistent.')
        database_contents = safe_read_bytes(work['database'])
        if previous.get('database_sha256') != hashlib.sha256(database_contents).hexdigest():
            raise ToolError('Compilation database differs from the successful capture.')
        entries = json.loads(database_contents)
        validate_database(entries, captured_sources(work['source']), compiler)
        return receipt_digest
    except (ToolError, OSError, ValueError, TypeError, KeyError) as error:
        raise ToolError('Cannot verify retained analysis: {}; run compdb with the same --profile before tidy.'.format(error)) from error


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['doctor', 'build', 'compdb', 'tidy', 'sanitize'])
    parser.add_argument('--profile', choices=['clang', 'clang-serial', 'gcc'], default='clang')
    parser.add_argument('--jobs', type=int, choices=[1, 2], default=2)
    parser.add_argument('--sanitizer-smoke', action='store_true', help='Prove sanitizer detection with a failing private fixture')
    args = parser.parse_args(argv)
    if args.command == 'sanitize' and args.profile == 'gcc':
        parser.error('sanitize requires --profile clang or clang-serial')
    if args.sanitizer_smoke and args.command != 'sanitize':
        parser.error('--sanitizer-smoke requires sanitize')
    sanitized = args.command == 'sanitize'
    profile = args.profile + ('-sanitize' if sanitized else '')
    base = OUTPUT / profile
    try:
        reject_symlink_tree(OUTPUT)
        reject_symlink_tree(base)
    except ToolError as error:
        print('cpp-tooling: ' + str(error), file=sys.stderr)
        return 1
    safe_mkdir(base)
    work = {'base': base, 'source': base / 'source', 'jobs': args.jobs,
            'database': base / 'compile_commands.json',
            'log': base / ('{}-{}.log'.format(args.command, time.time_ns()))}
    receipt = None
    provenance = None
    try:
        tools = config(args.profile, sanitized)
        env, compiler_version, r_version = environment(work, tools, profile)
        if args.command == 'tidy':
            receipt_digest = verify_capture(work, profile, compiler_version, r_version, tools['compiler'])
            owned = owned_sources(work['source'])
            if not owned:
                raise ToolError('No owned translation units to analyze; refresh the source with compdb.')
            tidy = require_tool('clang-tidy-18')
            for source in owned:
                run([tidy, '-p', str(base), '--config-file=' + str(work['source'] / '.clang-tidy'), str(source)],
                    work['source'], env, work['log'])
            verify_capture(work, profile, compiler_version, r_version, tools['compiler'], receipt_digest)
            print('Owned clang-tidy checks passed: {} translation units; {}'.format(len(owned), work['log']))
            return 0
        omp = openmp_probe(work, tools, args.profile, env) if IS_R else ''
        if args.command == 'doctor':
            runtime_tools = config('clang', True)
            tools['asan_runtime'] = runtime_tools['asan_runtime']
            for name in ('bear', 'clang-tidy-18', 'llvm-symbolizer-18'):
                tools[name] = require_tool(name)
            if IS_R:
                run([tools['Rscript'], '--vanilla', '-e', 'for(p in c("Rcpp","testthat")) if(!requireNamespace(p,quietly=TRUE)) stop(p," is required")'], base, env, work['log'])
            print(json.dumps({'tools': tools, 'compiler': compiler_version, 'R': r_version,
                              'openmp': omp or 'serial', 'output': str(base)}, indent=2))
            return 0
        if args.sanitizer_smoke:
            smoke_tools = config('clang', True)
            sanitizer_smoke(work, smoke_tools, env)
        # Every fresh attempt invalidates the old proof/database before replacing source.
        receipt = base / 'provenance.json'
        fingerprint = input_fingerprint(ROOT)
        provenance = {'command': args.command, 'profile': profile, 'tools': tools, 'fingerprint': fingerprint,
                      'compiler': compiler_version, 'R': r_version, 'jobs': args.jobs,
                      'environment': env, 'openmp': omp or 'serial', 'sources': {},
                      'log': str(work['log']), 'status': 'running',
                      'sanitizer_note': 'R host is unsanitized; package address/UB checks use matching Clang runtime. R-host leak detection disabled; standalone Stata tests check leaks.'}
        safe_write_text(receipt, json.dumps(provenance, indent=2) + '\n')
        if work['database'].exists():
            reject_symlink_path(work['database'])
            work['database'].unlink()
        snapshot(ROOT, work['source'])
        verify_inputs(fingerprint, work['source'], args.command)
        provenance['sources'] = {str(path.relative_to(work['source'])): file_sha256(path)
                                 for path in work['source'].rglob('*') if path.is_file()}
        safe_write_text(receipt, json.dumps(provenance, indent=2) + '\n')
        before = subprocess.run(
            [tools['ccache'], '--print-stats'], env=dict(os.environ, **env), text=True, capture_output=True, check=True).stdout
        builder = build_r if IS_R else build_stata
        if IS_R:
            builder(work, tools, env, omp, sanitized, args.command == 'compdb')
        else:
            builder(work, tools, env, sanitized, args.command == 'compdb')
        verify_inputs(fingerprint, work['source'], args.command)
        if args.command == 'compdb':
            database_contents = safe_read_bytes(work['database'])
            validate_database(json.loads(database_contents), captured_sources(work['source']), tools['compiler'])
            provenance['database_sha256'] = hashlib.sha256(database_contents).hexdigest()
        after = subprocess.run([tools['ccache'], '--print-stats'], env=dict(os.environ, **env), text=True, capture_output=True, check=True).stdout
        safe_write_text(base / 'cache-before.txt', before)
        safe_write_text(base / 'cache-after.txt', after)
        def hits(stats):
            parsed = dict(line.split() for line in stats.splitlines() if len(line.split()) == 2)
            return int(parsed.get('direct_cache_hit', 0)) + int(parsed.get('preprocessed_cache_hit', 0))
        cache_hits = hits(after) - hits(before)
        verify_inputs(fingerprint, work['source'], args.command)
        if args.command == 'compdb' and provenance['database_sha256'] != file_sha256(work['database']):
            raise ToolError('Compilation database changed during capture; run compdb again.')
        provenance.update(status='passed', environment=env, cache_hits=cache_hits)
        print('Compiler cache hits this fresh build: {}'.format(cache_hits))
        safe_write_text(receipt, json.dumps(provenance, indent=2) + '\n')
        print('{} passed ({}); log: {}'.format(args.command, profile, work['log']))
        return 0
    except (ToolError, OSError, ValueError, subprocess.SubprocessError) as error:
        if provenance is not None and receipt is not None:
            provenance.update(status='failed', error=str(error))
            safe_write_text(receipt, json.dumps(provenance, indent=2) + '\n')
        print('cpp-tooling: ' + str(error), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
