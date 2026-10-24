"""Behavioral tests for the isolated native developer command (no R/Stata needed)."""
import importlib.util
import contextlib
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'tools/dev/cpp_tooling.py'


class NativeToolingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('cpp_tooling', SCRIPT)
        cls.mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.mod)

    @contextlib.contextmanager
    def capture_fixture(self, capture=True):
        """Run the real capture/proof lifecycle with compilation replaced by a database writer."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'repo'
            native = root / ('src' if self.mod.IS_R else 'src/plugin')
            native.mkdir(parents=True)
            (root / 'tools/dev').mkdir(parents=True)
            (root / 'tools/dev/cpp_tooling.py').write_text('# fixture driver\n')
            (root / '.clang-tidy').write_text("Checks: 'clang-analyzer-core.*'\n")
            if self.mod.IS_R:
                (root / 'DESCRIPTION').write_text('Package: Fixture\n')
                (root / 'inst/include').mkdir(parents=True)
                (root / 'inst/include/public.hpp').write_text('int public_header;\n')
            (native / 'owned.cpp').write_text('int current_source;\n')
            (native / 'owned.hpp').write_text('int current_header;\n')
            (native / 'Makevars').write_text('CXX_STD = CXX17\n')
            sdk = 'RcppExports.cpp' if self.mod.IS_R else 'stplugin.c'
            (native / sdk).write_text('int generated_or_sdk;\n')
            output = Path(directory) / 'output'
            marker = Path(directory) / 'analyzed.txt'
            analyzer = Path(directory) / 'fixture-tidy'
            analyzer.write_text('#!' + sys.executable + '\n'
                                'import os, pathlib, sys\n'
                                'source = pathlib.Path(sys.argv[-1])\n'
                                'with pathlib.Path(os.environ["ANALYZER_MARKER"]).open("a") as out:\n'
                                '    out.write(str(source) + "\\n" + source.read_text())\n')
            analyzer.chmod(0o755)
            cache = Path(directory) / 'fixture-ccache'
            cache.write_text('#!' + sys.executable + '\nprint("direct_cache_hit 0")\n')
            cache.chmod(0o755)
            tools = {'compiler': '/usr/bin/clang++-18', 'ccache': str(cache)}
            def builder(work, *args):
                entries = [
                    {'directory': str(path.parent), 'file': str(path),
                     'arguments': [tools['compiler'], '-std=c++17', '-I' + str(path.parent), '-c', str(path)]}
                    for path in self.mod.captured_sources(work['source'])
                ]
                work['database'].write_text(json.dumps(entries))
            with contextlib.ExitStack() as stack:
                for name, value in (('ROOT', root), ('OUTPUT', output)):
                    stack.enter_context(mock.patch.object(self.mod, name, value))
                stack.enter_context(mock.patch.object(self.mod, 'config', return_value=tools))
                stack.enter_context(mock.patch.object(self.mod, 'environment', return_value=(
                    {'ANALYZER_MARKER': str(marker)}, 'fixture compiler', 'fixture R')))
                stack.enter_context(mock.patch.object(self.mod, 'openmp_probe', return_value=''))
                stack.enter_context(mock.patch.object(self.mod, 'build_r', side_effect=builder))
                stack.enter_context(mock.patch.object(self.mod, 'build_stata', side_effect=builder))
                stack.enter_context(mock.patch.object(self.mod, 'require_tool', return_value=str(analyzer)))
                fixture = {'root': root, 'base': output / 'clang', 'marker': marker,
                           'relative_native': native.relative_to(root), 'builder': builder}
                if capture:
                    with contextlib.redirect_stdout(io.StringIO()):
                        self.assertEqual(self.mod.main(['compdb']), 0)
                yield fixture

    def assert_tidy_rejected(self, fixture):
        error = io.StringIO()
        with contextlib.redirect_stderr(error), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(self.mod.main(['tidy']), 1, error.getvalue())
        self.assertIn('compdb', error.getvalue())
        self.assertFalse(fixture['marker'].exists(), 'Analyzer must not run on unverified capture')

    def test_unchanged_capture_analyzes_the_retained_source(self):
        with self.capture_fixture() as fixture:
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(self.mod.main(['tidy']), 0)
            analyzed = fixture['marker'].read_text()
            self.assertIn(str(fixture['base'] / 'source'), analyzed)
            self.assertIn('int current_source;', analyzed)

    def test_tidy_rejects_changed_or_missing_retained_native_and_config(self):
        names = ('owned.cpp', 'owned.hpp', 'Makevars',
                 'RcppExports.cpp' if self.mod.IS_R else 'stplugin.c')
        relatives = [Path('src' if self.mod.IS_R else 'src/plugin') / name for name in names]
        relatives += [Path('.clang-tidy'), Path('tools/dev/cpp_tooling.py')]
        if self.mod.IS_R:
            relatives += [Path('DESCRIPTION'), Path('inst/include/public.hpp')]
        for relative in relatives:
            for action in ('edit', 'delete'):
                with self.subTest(path=relative, action=action), self.capture_fixture() as fixture:
                    retained = fixture['base'] / 'source' / relative
                    if action == 'edit':
                        retained.write_text('altered retained input\n')
                    else:
                        retained.unlink()
                    self.assert_tidy_rejected(fixture)

    def test_tidy_rejects_live_native_and_config_changes_before_analysis(self):
        for relative in (Path('.clang-tidy'), Path('tools/dev/cpp_tooling.py'),
                         Path('src' if self.mod.IS_R else 'src/plugin') / 'owned.hpp'):
            for action in ('edit', 'delete'):
                with self.subTest(path=relative, action=action), self.capture_fixture() as fixture:
                    live = fixture['root'] / relative
                    if action == 'edit':
                        live.write_text('changed live input\n')
                    else:
                        live.unlink()
                    self.assert_tidy_rejected(fixture)

    def test_tidy_rejects_edited_database_even_when_command_contract_still_valid(self):
        for action in ('flags', 'whitespace', 'missing', 'invalid-json'):
            with self.subTest(action=action), self.capture_fixture() as fixture:
                database = fixture['base'] / 'compile_commands.json'
                if action == 'flags':
                    entries = json.loads(database.read_text())
                    entries[0]['arguments'].insert(1, '-DCHANGED_INPUT=1')
                    database.write_text(json.dumps(entries))
                elif action == 'whitespace':
                    database.write_text(database.read_text() + '\n')
                elif action == 'missing':
                    database.unlink()
                else:
                    database.write_text('invalid JSON')
                self.assert_tidy_rejected(fixture)

    def test_tidy_rejects_failed_missing_or_edited_capture_proof(self):
        for action in ('missing', 'failed', 'running', 'build', 'profile', 'sources',
                       'fingerprint', 'database_sha256', 'invalid-json', 'non-object'):
            with self.subTest(action=action), self.capture_fixture() as fixture:
                receipt = fixture['base'] / 'provenance.json'
                previous = json.loads(receipt.read_text())
                if action == 'missing':
                    receipt.unlink()
                elif action == 'invalid-json':
                    receipt.write_text('invalid JSON')
                elif action == 'non-object':
                    receipt.write_text('[]')
                else:
                    if action in ('failed', 'running'):
                        previous['status'] = action
                    elif action == 'build':
                        previous['command'] = 'build'
                    elif action == 'profile':
                        previous['profile'] = 'gcc'
                    elif action == 'sources':
                        name = str(fixture['relative_native'] / 'owned.cpp')
                        previous['sources'][name] = '0' * 64
                    elif action == 'fingerprint':
                        previous.pop('fingerprint')
                    else:
                        previous.pop('database_sha256', None)
                    receipt.write_text(json.dumps(previous))
                self.assert_tidy_rejected(fixture)

    def test_snapshot_live_edit_cannot_certify_a_different_retained_copy(self):
        with self.capture_fixture(capture=False) as fixture:
            original = self.mod.snapshot
            def racing_snapshot(root, target):
                original(root, target)
                (root / fixture['relative_native'] / 'owned.cpp').write_text('int newer_live_source;\n')
            error = io.StringIO()
            with mock.patch.object(self.mod, 'snapshot', side_effect=racing_snapshot), contextlib.redirect_stderr(error), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(self.mod.main(['compdb']), 1)
            self.assertIn('compdb', error.getvalue())
            self.assertFalse((fixture['base'] / 'compile_commands.json').exists())
            self.assert_tidy_rejected(fixture)

    def test_capture_mutation_and_failed_recapture_never_reuse_prior_success(self):
        for action in ('retained', 'live', 'failure', 'no-database'):
            with self.subTest(action=action), self.capture_fixture() as fixture:
                def unsuccessful_capture(work, *args):
                    if action == 'failure':
                        raise self.mod.ToolError('fixture build failure')
                    if action != 'no-database':
                        fixture['builder'](work, *args)
                        parent = work['source'] if action == 'retained' else fixture['root']
                        (parent / fixture['relative_native'] / 'owned.hpp').write_text('int edited_during_build;\n')
                name = 'build_r' if self.mod.IS_R else 'build_stata'
                with mock.patch.object(self.mod, name, side_effect=unsuccessful_capture), contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(self.mod.main(['compdb']), 1)
                self.assert_tidy_rejected(fixture)

    def test_tidy_does_not_accept_inputs_changed_during_analysis(self):
        for action in ('retained', 'live', 'database', 'receipt'):
            with self.subTest(action=action), self.capture_fixture() as fixture:
                original = self.mod.run
                def racing_analysis(command, cwd, env, logfile):
                    original(command, cwd, env, logfile)
                    if action in ('retained', 'live'):
                        parent = fixture['base'] / 'source' if action == 'retained' else fixture['root']
                        changed = parent / fixture['relative_native'] / 'owned.hpp'
                    else:
                        changed = fixture['base'] / ('compile_commands.json' if action == 'database' else 'provenance.json')
                    changed.write_text(changed.read_text() + '\n')
                error = io.StringIO()
                with mock.patch.object(self.mod, 'run', side_effect=racing_analysis), contextlib.redirect_stderr(error), contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(self.mod.main(['tidy']), 1)
                self.assertIn('compdb', error.getvalue())
                self.assertTrue(fixture['marker'].exists(), 'Fixture must reach analysis before the race')

    @unittest.skipUnless(shutil.which('clang-tidy-18') and shutil.which('clang++-18'),
                         'Clang 18 tools are required for the real analyzer control')
    def test_real_analyzer_detects_current_bug_and_rejects_altered_retained_copy(self):
        with self.capture_fixture(capture=False) as fixture:
            live = fixture['root'] / fixture['relative_native'] / 'owned.cpp'
            live.write_text('int f() { int *p = nullptr; return *p; }\n')
            (fixture['root'] / '.clang-tidy').write_text(
                "Checks: '-*,clang-analyzer-core.*'\nWarningsAsErrors: '*'\n")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(self.mod.main(['compdb']), 0)
            error = io.StringIO()
            with mock.patch.object(self.mod, 'require_tool', return_value=shutil.which('clang-tidy-18')):
                with contextlib.redirect_stderr(error), contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(self.mod.main(['tidy']), 1)
                logs = list(fixture['base'].glob('tidy-*.log'))
                self.assertEqual(len(logs), 1)
                self.assertIn('clang-analyzer-core.NullDereference', logs[0].read_text())
                retained = fixture['base'] / 'source' / fixture['relative_native'] / 'owned.cpp'
                retained.write_text('int f() { return 0; }\n')
                self.assert_tidy_rejected(fixture)
                self.assertEqual(list(fixture['base'].glob('tidy-*.log')), logs)
                live.write_text('int f() { return 0; }\n')
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(self.mod.main(['compdb']), 0)
                    self.assertEqual(self.mod.main(['tidy']), 0)
                self.assertEqual(len(list(fixture['base'].glob('tidy-*.log'))), 2)

    def test_installed_r_proof_requires_current_production_primitives_and_real_dll(self):
        # This inspects the runner contract only. Actual native execution is
        # established by the isolated R install/build, never by this mock.
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            with mock.patch.object(self.mod, 'run') as execute:
                self.mod.test_r({'base': base, 'source': base, 'log': base / 'run.log'},
                                {'Rscript': 'Rscript'}, {}, False)
            script = (base / 'installed-tests.R').read_text()
            for name in ('rustat_subset_row_sums_g', 'rustat_exact_rank_table',
                         'rustat_compensated_sum', 'rustat_compensated_crossprod',
                         'rustat_block_rank_loo', 'rustat_jackknife_fast_rank',
                         'rustat_jackknife_rank_full', 'rustat_plqr_ctx_create', 'rustat_plqr_g_all'):
                self.assertIn('"' + name + '"', script)
            for retired in ('rustat_cate_ctx_create', 'rustat_cate_g_all',
                            'rustat_cate_g_indexed', 'rustat_plqr_ctx_set_beta', 'rustat_plqr_g_indexed'):
                self.assertNotIn(retired, script)
            self.assertIn('getLoadedDLLs()[["RdistNN"]][["path"]]', script)
            self.assertIn('normalizePath(file.path(args[[1]], "RdistNN", "libs"', script)
            self.assertIn('expectation_skip', script)
            self.assertIn('No installed native tests ran', script)
            selection = execute.call_args.args[0][-1]
            for suite in ('native-index-boundaries', 'native-plqr-threshold-recovery',
                          'errors-fastpaths', 'score-plqr-regressions', 'seeded-rng-scope'):
                self.assertIn(suite, selection)

    def test_openmp_exception_checker_remains_enabled(self):
        config = (ROOT / '.clang-tidy').read_text()
        self.assertIn('openmp-exception-escape', config)
        for checker in ('clang-analyzer-core.*', 'clang-analyzer-cplusplus.*',
                        'clang-analyzer-deadcode.*', 'bugprone-use-after-move'):
            self.assertIn(checker, config)
        self.assertIn("WarningsAsErrors: '*'", config)

    def test_help_outside_working_directory(self):
        result = subprocess.run([sys.executable, str(SCRIPT), '--help'], cwd='/tmp', capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('sanitize', result.stdout)

    def test_invalid_jobs_and_profile_fail_before_building(self):
        for args in [('build', '--jobs', '3'), ('build', '--profile', 'unknown')]:
            result = subprocess.run([sys.executable, str(SCRIPT), *args], cwd='/tmp', capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn('invalid choice', result.stderr)

    def test_missing_tool_is_actionable(self):
        with self.assertRaisesRegex(self.mod.ToolError, 'Required tool.*not found'):
            self.mod.require_tool('rdistnn-tool-that-does-not-exist')

    def test_failed_command_propagates_and_retains_log(self):
        with tempfile.TemporaryDirectory() as directory:
            logfile = Path(directory) / 'failure.log'
            with self.assertRaisesRegex(self.mod.ToolError, 'exit 7'):
                self.mod.run([sys.executable, '-c', 'print("fixture failed"); raise SystemExit(7)'], Path(directory), {}, logfile)
            self.assertIn('fixture failed', logfile.read_text())

    def test_snapshot_excludes_outputs_preserves_current_source(self):
        for is_r in (True, False):
            with self.subTest(is_r=is_r), tempfile.TemporaryDirectory() as directory, mock.patch.object(self.mod, 'IS_R', is_r):
                root = Path(directory) / 'repo'
                root.mkdir()
                (root / 'src').mkdir()
                (root / 'src/current.cpp').write_text('current uncommitted source')
                (root / 'src/stale.o').write_text('object')
                (root / 'archives').mkdir()
                (root / 'archives/private').write_text('private')
                (root / '.superpowers').mkdir()
                (root / '.superpowers/private').write_text('private')
                (root / '_dist').mkdir()
                (root / '_dist/source.zip').write_bytes(b'generated archive')
                (root / 'docs').mkdir()
                (root / 'docs/index.html').write_text('generated site or source documentation')
                (root / 'docs/article.md').write_text('current Stata article')
                (root / 'inst/include').mkdir(parents=True)
                (root / 'inst/include/current.hpp').write_text('current native header')
                target = Path(directory) / 'snapshot'
                self.mod.snapshot(root, target)
                self.assertEqual((target / 'src/current.cpp').read_text(), 'current uncommitted source')
                self.assertFalse((target / 'src/stale.o').exists())
                self.assertFalse((target / 'archives').exists())
                self.assertFalse((target / '.superpowers').exists())
                self.assertFalse((target / '_dist').exists())
                self.assertEqual((target / 'inst/include/current.hpp').read_text(), 'current native header')
                if is_r:
                    self.assertFalse((target / 'docs').exists())
                else:
                    self.assertEqual((target / 'docs/article.md').read_text(), 'current Stata article')
                    self.assertEqual((target / 'docs/index.html').read_text(), 'generated site or source documentation')
                self.assertTrue((root / '_dist/source.zip').exists())
                self.assertTrue((root / 'docs/article.md').exists())
                self.assertTrue((root / 'src/stale.o').exists())

    def test_snapshot_rejects_destructive_destination_and_symlink_escape(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'repo'
            root.mkdir()
            source = root / 'current.cpp'
            source.write_text('keep source')
            with self.assertRaisesRegex(self.mod.ToolError, 'replace source'):
                self.mod.snapshot(root, root)
            outside = Path(directory) / 'outside'
            outside.mkdir()
            (root / 'linked').symlink_to(outside, target_is_directory=True)
            with self.assertRaisesRegex(self.mod.ToolError, 'Source symlink'):
                self.mod.snapshot(root, Path(directory) / 'copy')
            self.assertEqual(source.read_text(), 'keep source')
            self.assertTrue(outside.is_dir())

    def test_compiler_invocation_preserves_cplusplus_driver_name(self):
        # clang++ often resolves to clang; resolving that symlink loses stdlib linking.
        with mock.patch.object(self.mod.shutil, 'which', return_value='/usr/bin/clang++-18'):
            self.assertEqual(self.mod.require_tool('clang++-18'), '/usr/bin/clang++-18')

    def test_fingerprint_detects_source_and_config_changes_ignores_objects(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            native = root / ('src' if self.mod.IS_R else 'src/plugin')
            native.mkdir(parents=True)
            source = native / 'owned.cpp'
            source.write_text('int first;')
            config = root / '.clang-tidy'
            config.write_text('Checks: clang-analyzer-core.*')
            initial = self.mod.input_fingerprint(root)
            (native / 'old.o').write_text('cached object')
            self.assertEqual(initial, self.mod.input_fingerprint(root))
            source.write_text('int second;')
            changed = self.mod.input_fingerprint(root)
            self.assertNotEqual(initial, changed)
            config.write_text('Checks: bugprone-use-after-move')
            self.assertNotEqual(changed, self.mod.input_fingerprint(root))

    def test_database_rejects_missing_owned_translation_unit(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'owned.cpp'
            source.write_text('int x;')
            with self.assertRaisesRegex(self.mod.ToolError, 'missing.*owned.cpp'):
                self.mod.validate_database([], [source])

    def test_cache_symlinks_never_overwrite_external_configuration(self):
        for kind in ('directory', 'config', 'descendant'):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as directory:
                output = Path(directory) / 'output'
                outside = Path(directory) / 'outside'
                outside.mkdir()
                sentinel = outside / 'ccache.conf'
                sentinel.write_text('external configuration must survive')
                version = 'fixture clang'
                r_version = 'fixture R' if self.mod.IS_R else 'standalone'
                key = hashlib.sha256((version + r_version + 'clang').encode()).hexdigest()[:16]
                cache = output / 'cache' / key
                cache.parent.mkdir(parents=True)
                if kind == 'directory':
                    cache.symlink_to(outside, target_is_directory=True)
                else:
                    cache.mkdir()
                    if kind == 'config':
                        (cache / 'ccache.conf').symlink_to(sentinel)
                    else:
                        (cache / 'nested').symlink_to(outside, target_is_directory=True)
                versions = [version, r_version] if self.mod.IS_R else [version]
                tools = {'compiler': 'fixture-cxx', 'R': 'fixture-R'}
                with mock.patch.object(self.mod, 'OUTPUT', output), mock.patch.object(self.mod, 'capture', side_effect=versions):
                    with self.assertRaisesRegex(self.mod.ToolError, 'symlink'):
                        self.mod.environment({'jobs': 2}, tools, 'clang')
                self.assertEqual(sentinel.read_text(), 'external configuration must survive')

    def test_profile_descendant_links_fail_before_tools_or_output_changes(self):
        outputs = ('doctor-1234.log', 'Makevars', 'objects', 'objects/stale.o',
                   'library', 'library/RdistNN/DESCRIPTION', 'compile_commands.json',
                   'compile_commands.events.json', 'provenance.json', 'cache-before.txt',
                   'installed-tests.R', 'developer.mk', 'openmp-probe.cpp', 'openmp-probe',
                   'source', 'source/src/owned.cpp')
        for name in outputs:
            with self.subTest(output=name), tempfile.TemporaryDirectory() as directory:
                output = Path(directory) / 'output'
                base = output / 'clang'
                base.mkdir(parents=True)
                safe = base / 'safe-snapshot.txt'
                safe.write_text('preserve safe snapshot')
                external = Path(directory) / 'external'
                external.write_text('do not read or overwrite')
                link = base / name
                link.parent.mkdir(parents=True, exist_ok=True)
                link.symlink_to(external)
                error = io.StringIO()
                with mock.patch.object(self.mod, 'OUTPUT', output), mock.patch.object(self.mod.time, 'time_ns', return_value=1234), mock.patch.object(self.mod, 'config', side_effect=AssertionError('must reject before tools')), contextlib.redirect_stderr(error):
                    self.assertEqual(self.mod.main(['doctor']), 1)
                self.assertIn('symlink', error.getvalue())
                self.assertEqual(external.read_text(), 'do not read or overwrite')
                self.assertEqual(safe.read_text(), 'preserve safe snapshot')

    def test_sanitizer_smoke_rejects_unrelated_failure_after_old_detection(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            (base / 'sanitizer-negative.log').write_text('OLD AddressSanitizer detection')
            tools = {'compiler': 'fixture-cxx', 'symbolizer': 'fixture-symbolizer'}
            def unrelated_failure(command, cwd, env, logfile):
                if command[0] == 'fixture-cxx':
                    return
                with logfile.open('a') as log:
                    log.write('current unrelated launch failure: AddressSanitizer:DEADLYSIGNAL')
                raise self.mod.ToolError('exit 127 unrelated launch failure')
            with mock.patch.object(self.mod, 'run', side_effect=unrelated_failure):
                with self.assertRaisesRegex(self.mod.ToolError, 'unrelated launch failure'):
                    self.mod.sanitizer_smoke({'base': base, 'log': base / 'compile.log'}, tools, {})

    def test_database_checks_real_command_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'owned.cpp'
            source.write_text('int x;')
            entry = {'directory': directory, 'file': str(source),
                     'arguments': ['/usr/bin/clang++-18', '-std=c++17', '-I' + directory, '-c', str(source)]}
            self.mod.validate_database([entry], [source], '/usr/bin/clang++-18')
            with self.assertRaisesRegex(self.mod.ToolError, 'compiler differs'):
                self.mod.validate_database([entry], [source], '/usr/bin/g++')
            entry['arguments'].remove('-std=c++17')
            with self.assertRaisesRegex(self.mod.ToolError, 'C[+][+]17'):
                self.mod.validate_database([entry], [source])
            entry['arguments'].insert(1, '-std=c++17')
            entry['arguments'].remove('-I' + directory)
            with self.assertRaisesRegex(self.mod.ToolError, 'include flags'):
                self.mod.validate_database([entry], [source])

    def test_database_rejects_stale_paths(self):
        entry = {'directory': '/tmp', 'file': '/tmp/rdistnn-nonexistent.cpp', 'arguments': ['clang++-18', '-I/tmp', '-c', '/tmp/rdistnn-nonexistent.cpp']}
        with self.assertRaisesRegex(self.mod.ToolError, 'does not exist'):
            self.mod.validate_database([entry], [])


if __name__ == '__main__':
    unittest.main()
