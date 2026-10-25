import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


class NativeThreadPolicyTests(unittest.TestCase):
    def test_all_parallel_regions_use_the_package_limit(self):
        regions = []
        for path in (ROOT / "src").glob("*.cpp"):
            for number, line in enumerate(path.read_text().splitlines(), 1):
                if re.match(r"\s*#pragma omp parallel\b", line):
                    regions.append((path.name, number))
                    self.assertIn("num_threads(rustat::thread_limit())", line,
                                  f"Unbounded native team at {path.name}:{number}")
            self.assertNotRegex(path.read_text(), r"\bomp_set_num_threads\s*\(")
        self.assertTrue(regions)

    def test_exact_rank_table_host_callback_preserves_values_and_propagates_cancellation(self):
        compiler = shutil.which("c++")
        if compiler is None:
            self.skipTest("C++ compiler unavailable for portable callback integration")
        source = r'''
#include "rank_weight_tables.hpp"
#include <stdexcept>
static int calls = 0;
static void observed() { ++calls; }
static void cancel() {
  if (++calls == 8) throw std::runtime_error("host cancellation");
}
int main() {
  for (int k : {1, 3, 6}) {
    const auto ordinary = rdistnn_weights::exact_table(24, k, 6);
    calls = 0;
    const auto checked = rdistnn_weights::exact_table(24, k, 6, observed);
    if (calls == 0 || ordinary != checked) return 1;
  }
  // These long phases span multiple chunks. Counts require interior checkpoints;
  // retaining only the first recurrence/phase check must fail this test.
  for (int k : {1, 2}) {
    calls = 0;
    rdistnn_weights::exact_table(262144, k, k + 1, observed);
    if (calls != (k == 1 ? 30 : 47)) return 4;
  }
  calls = 0;
  rdistnn_weights::exact_table(262144, 2, 2, observed);
  if (calls != 5) return 5;
  calls = 0;
  try {
    rdistnn_weights::exact_table(262144, 1, 2, cancel);
  } catch (const std::runtime_error& error) {
    return calls == 8 && std::string(error.what()) == "host cancellation" ? 0 : 2;
  }
  return 3;
}
'''
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "callback.cpp").write_text(source)
            subprocess.run([compiler, "-std=c++17", "-O2", "-I", str(ROOT / "inst/include"),
                            str(path / "callback.cpp"), "-o", str(path / "callback")], check=True,
                           capture_output=True, text=True)
            subprocess.run([str(path / "callback")], check=True, capture_output=True, text=True)
