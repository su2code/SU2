#!/usr/bin/env python3

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

for module in (
    "PyQt5",
    "PyQt5.QtCore",
    "PyQt5.QtGui",
    "PyQt5.QtWidgets",
    "matplotlib.backends.backend_qt5agg",
):
    sys.modules.setdefault(module, unittest.mock.MagicMock())

MODULE_PATH = Path(__file__).parents[2] / "bin" / "monitor.py"
SPEC = importlib.util.spec_from_file_location("su2_monitor", MODULE_PATH)
MONITOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MONITOR)

NAMES = '"Inner_Iter","rms[Rho]","rms[RhoE]","CD"'
LEGACY = 'TITLE = "SU2 history"\nVARIABLES = ' + NAMES + '\nZONE T= "Convergence"\n'
TECPLOT = 'TITLE = "SU2 history"\nVARIABLES = \\\n' + NAMES + '\nZONE T= "Convergence"\n'
CSV = NAMES + "\n"


def row(i):
    return f"{i}, {-1.0 - 0.1 * i}, {-2.0 - 0.1 * i}, 0.5\n"


class MonitorTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.path = Path(self.dir.name) / "history.csv"

    def tearDown(self):
        self.dir.cleanup()

    def reader(self, text):
        self.path.write_text(text)
        reader = MONITOR.HistoryReader(str(self.path))
        reader.update()
        return reader

    def test_header_formats(self):
        for header in (LEGACY, TECPLOT, CSV):
            reader = self.reader(header + row(0) + row(1))
            self.assertEqual(reader.names, ["rms[Rho]", "rms[RhoE]"], header)
            self.assertEqual(reader.n, 2, header)

    def test_time_iterations(self):
        names = '"Time_Iter","Inner_Iter","rms[Rho]"\n'
        body = "".join(f"{t}, {i}, {-t - 0.1 * i}\n" for t in range(3) for i in range(4))
        reader = self.reader(names + body)
        self.assertTrue(reader.hasTime())
        np.testing.assert_array_equal(reader.iters[: reader.n][:5], [0, 1, 2, 3, 0])
        rows = reader.rowsFrom(1, byTime=True)
        np.testing.assert_array_equal(reader.times[rows], [1, 2])
        np.testing.assert_array_equal(reader.iters[rows], [3, 3])
        self.assertFalse(self.reader(CSV + row(0) + row(1)).hasTime())

    def test_incomplete_header_is_deferred(self):
        self.path.write_text('TITLE = "x"\nVARIABLES = \\\n"Inner_Iter","rms[Rho]')
        reader = MONITOR.HistoryReader(str(self.path))
        self.assertFalse(reader.update())
        self.assertIsNone(reader.cols)
        self.path.write_text(TECPLOT + row(0))
        self.assertTrue(reader.update())
        self.assertEqual(reader.n, 1)

    def test_sparse_iterations_and_start_filter(self):
        reader = self.reader(CSV + row(0) + row(10) + row(20) + row(30))
        np.testing.assert_array_equal(reader.iters[: reader.n], [0, 10, 20, 30])
        rows = reader.rowsFrom(15)
        np.testing.assert_array_equal(reader.iters[rows], [20, 30])

    def test_incremental_append_and_incomplete_line(self):
        reader = self.reader(CSV + row(0) + row(1))
        with open(self.path, "a") as f:
            f.write(row(2)[:5])
        self.assertFalse(reader.update())
        self.assertEqual(reader.n, 2)
        with open(self.path, "a") as f:
            f.write(row(2)[5:] + row(3))
        version = reader.version
        self.assertTrue(reader.update())
        self.assertEqual(reader.n, 4)
        self.assertGreater(reader.version, version)

    def test_reset_invalidates_cache_key(self):
        reader = self.reader(CSV + "".join(row(i) for i in range(50)))
        version = reader.version
        self.path.write_text(CSV + row(0) + row(1))
        reader.update()
        self.assertEqual(reader.n, 2)
        self.assertGreater(reader.version, version)

    def test_reset_without_complete_rows_changes_version(self):
        reader = self.reader(CSV + "".join(row(i) for i in range(50)))
        version = reader.version
        self.path.write_text("")
        reader.update()
        self.assertGreater(reader.version, version)


if __name__ == "__main__":
    unittest.main()
