import re
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
