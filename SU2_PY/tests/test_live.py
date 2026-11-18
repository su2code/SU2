#!/usr/bin/env python3

## \file test_live.py
#  \brief Unit tests for the SU2_LIVE terminal monitor.
#  \version 8.6.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import importlib.machinery
import importlib.util
import math
import sys
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = Path(__file__).parents[2] / "bin" / "SU2_LIVE"
LOADER = importlib.machinery.SourceFileLoader("su2_live", str(MODULE_PATH))
SPEC = importlib.util.spec_from_loader("su2_live", LOADER)
LIVE = importlib.util.module_from_spec(SPEC)
sys.modules["su2_live"] = LIVE
LOADER.exec_module(LIVE)

STEADY = list(range(3000))

HEADER = (
    '"Time_Iter","Outer_Iter","Inner_Iter",    "rms[Rho]"    ,'
    '    "rms[RhoU]"   ,    "rms[RhoV]"   ,    "rms[RhoE]"    '
)


def exp_decay(n=3000, tau=1800, floor=-12.0, start=-2.0):
    return [start + (floor - start) * (1 - math.exp(-i / tau)) for i in range(n)]


def rows(n, unsteady=False):
    out = []
    for i in range(n):
        t, k = (i // 15, i % 15) if unsteady else (0, i)
        b = -2 - 7 * (1 - math.exp(-i / 500))
        out.append(
            "%12d,           0,%12d,%18.9f,%18.9f,%18.9f,%18.9f\n"
            % (t, k, b, b + 2, b + 1.5, b + 3.5)
        )
    return "".join(out)


class AnalysisTests(unittest.TestCase):
    """Classification of residual histories seen on SU2 test cases."""

    def status(self, ys, xs=STEADY):
        return LIVE.analyse(xs, ys, -12).status

    def test_converging(self):
        self.assertEqual(self.status(exp_decay()), LIVE.CONVERGING)

    def test_slow_but_converging(self):
        d = LIVE.analyse(
            STEADY, [-6 - 4 * (1 - math.exp(-i / 28000)) for i in STEADY], -12
        )
        self.assertEqual(d.status, LIVE.CONVERGING)
        self.assertGreater(d.eta, 10000)

    def test_frozen(self):
        ys = [
            -7.580046571 if i > 800 else -2 - 6 * (1 - math.exp(-i / 300))
            for i in STEADY
        ]
        self.assertEqual(self.status(ys), LIVE.FROZEN)

    def test_slow_limit_cycle(self):
        self.assertEqual(
            self.status([-9 + 1.4 * math.sin(i / 260) for i in STEADY]), LIVE.CYCLING
        )

    def test_fast_limit_cycle(self):
        self.assertEqual(
            self.status([-9 + 0.5 * math.sin(i / 40) for i in STEADY]), LIVE.CYCLING
        )

    def test_noisy_converging(self):
        # explicit schemes alternate from one iteration to the next on top of
        # a slow decrease, which must not be read as a limit cycle
        ys = [-1.5 - 0.003 * i + 0.03 * (-1) ** i for i in range(300)]
        self.assertEqual(self.status(ys, list(range(300))), LIVE.CONVERGING)

    def test_cycle_after_transient(self):
        # an initial drop followed by a slow cycle, as seen with ROE and the
        # Venkatakrishnan limiter on the supersonic biparabolic case
        ys = [-1.6 - 2.2 * (1 - math.exp(-i / 80)) for i in range(300)]
        ys += [-3.82 + 0.18 * math.sin(2 * math.pi * i / 150) for i in range(700)]
        self.assertEqual(self.status(ys, list(range(1000))), LIVE.CYCLING)

    def test_diverging(self):
        self.assertEqual(
            self.status([-5 + 0.004 * i if i > 500 else -5 for i in STEADY]),
            LIVE.DIVERGING,
        )

    def test_stalled(self):
        self.assertEqual(
            self.status([-8.0 + 0.001 * math.sin(i / 50) for i in STEADY]), LIVE.STALLED
        )

    def test_nan_is_not_converging(self):
        ys = exp_decay(1200, tau=300, floor=-5) + [float("nan")] * 1800
        self.assertEqual(self.status(ys), LIVE.DIVERGING)

    def test_nan_early_then_finite(self):
        self.assertEqual(
            self.status([float("nan")] * 50 + exp_decay(2950)), LIVE.CONVERGING
        )

    def test_constant_iteration_column(self):
        ys = [-2 - 6 * (1 - math.exp(-i / 400)) for i in range(2000)]
        self.assertEqual(self.status(ys, [0] * 2000), LIVE.CONVERGING)

    def test_sawtooth_iteration_column(self):
        xs, ys = [], []
        for t in range(200):
            for k in range(20):
                xs.append(k)
                ys.append(-2 - 6 * (1 - math.exp(-(t * 20 + k) / 900)))
        d = LIVE.analyse(xs, ys, -12)
        self.assertEqual(d.status, LIVE.CONVERGING)
        self.assertLess(d.rate, 0)

    def test_slope_handles_degenerate_x(self):
        self.assertEqual(LIVE.slope([1, 1, 1], [1, 2, 3]), 0.0)


class HistoryTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.path = Path(self.dir.name) / "history.csv"

    def tearDown(self):
        self.dir.cleanup()

    def read(self, text):
        self.path.write_text(text)
        hist = LIVE.History(str(self.path))
        hist.poll()
        return hist

    def test_reads_columns(self):
        hist = self.read(HEADER + "\n" + rows(100))
        self.assertEqual(hist.available(), ["Rho", "RhoU", "RhoV", "RhoE"])
        self.assertEqual(len(hist.iters), 100)

    def test_tecplot_header(self):
        hist = self.read("VARIABLES = \\\n" + HEADER + "\n" + rows(100, unsteady=True))
        self.assertEqual(hist.available(), ["Rho", "RhoU", "RhoV", "RhoE"])
        self.assertEqual(hist.iter_name, "Time_Iter")

    def test_legacy_tecplot_header(self):
        text = (
            'TITLE = "SU2"\nVARIABLES = '
            + HEADER
            + '\nZONE T= "Convergence history"\n'
            + rows(100, True)
        )
        hist = self.read(text)
        self.assertEqual(hist.header[0], "Time_Iter")
        self.assertEqual(hist.iter_name, "Time_Iter")
        self.assertEqual(len(hist.iters), 100)

    def test_unsteady_axis_is_monotonic(self):
        hist = self.read(HEADER + "\n" + rows(300, unsteady=True))
        self.assertTrue(all(b >= a for a, b in zip(hist.iters, hist.iters[1:])))
        self.assertGreater(len(set(hist.iters)), 1)

    def test_multizone_summary(self):
        header = '"Time_Iter","Outer_Iter","Inner_Iter",   "bgs[Rho][0]"  ,  "bgs[RhoE][0]"  ,    "bgs[T][1]"   \n'
        body = "".join(
            f"0, {i}, 0, {-1 - 0.01 * i}, {-2 - 0.01 * i}, {-3 - 0.01 * i}\n"
            for i in range(100)
        )
        hist = self.read(header + body)
        self.assertEqual(hist.available(), ["Rho[0]", "RhoE[0]", "T[1]"])
        self.assertEqual(hist.iter_name, "Outer_Iter")
        self.assertEqual(hist.iters[-1], 99)

    def test_rms_and_max_do_not_collide(self):
        header = '"Inner_Iter","rms[Rho]","max[Rho]"\n'
        hist = self.read(header + "".join(f"{i}, -1, -2\n" for i in range(10)))
        self.assertEqual(hist.available(), ["Rho", "max[Rho]"])

    def test_incremental_read(self):
        hist = self.read(HEADER + "\n" + rows(50))
        first = len(hist.iters)
        with open(self.path, "a") as f:
            f.write("0, 0, 99, -9.0, -7.0, -7.5, -5.5\n0, 0, 100, -9.1")
        hist.poll()
        self.assertEqual(len(hist.iters), first + 1)
        with open(self.path, "a") as f:
            f.write(", -7.1, -7.6, -5.6\n")
        hist.poll()
        self.assertEqual(len(hist.iters), first + 2)
        hist.poll()
        self.assertEqual(len(hist.iters), first + 2)

    def test_truncation_resets(self):
        hist = self.read(HEADER + "\n" + rows(50))
        self.path.write_text(HEADER + "\n" + rows(5))
        hist.poll()
        self.assertEqual(len(hist.iters), 5)

    def test_selection(self):
        hist = self.read(HEADER + "\n" + rows(50))
        self.assertEqual(
            [lab for lab, _ in hist.select(["Rho", "RhoE"])], ["Rho", "RhoE"]
        )
        self.assertEqual(
            [lab for lab, _ in hist.select(["nonsense"])], hist.available()
        )


class ConfigTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.dir.cleanup()

    def path_for(self, text, name="case.cfg"):
        cfg = Path(self.dir.name) / name
        cfg.write_text(text)
        return Path(LIVE.history_path(str(cfg))).name

    def test_default(self):
        self.assertEqual(self.path_for("SOLVER= EULER\n"), "history.csv")

    def test_tecplot(self):
        self.assertEqual(
            self.path_for("CONV_FILENAME= conv.csv\nTABULAR_FORMAT= TECPLOT\n"),
            "conv.dat",
        )

    def test_unsteady_restart(self):
        text = "TIME_DOMAIN= YES\nRESTART_SOL= YES\nRESTART_ITER= 20\n"
        self.assertEqual(self.path_for(text), "history_00020.csv")

    def test_multizone(self):
        text = "SOLVER= MULTIPHYSICS\nCONFIG_LIST= (a.cfg, b.cfg)\n"
        self.assertEqual(self.path_for(text, "cht_2d.cfg"), "cht_2d.csv")

    def test_commented_options_are_ignored(self):
        self.assertEqual(self.path_for("% TABULAR_FORMAT= TECPLOT\n"), "history.csv")

    def test_target(self):
        cfg = Path(self.dir.name) / "case.cfg"
        cfg.write_text("CONV_RESIDUAL_MINVAL= -10 % comment\n")
        self.assertEqual(LIVE.convergence_target(str(cfg)), -10.0)


class PlotTests(unittest.TestCase):
    def test_frame_has_requested_height(self):
        series = [("Rho", exp_decay(200)), ("RhoE", exp_decay(200, start=-1))]
        lines = LIVE.frame(
            series, list(range(200)), 60, 10, "title", annotations=["note"]
        )
        # title, 10 plot rows, axis line, x labels, legend, one annotation
        self.assertEqual(len(lines), 15)

    def test_frame_survives_nan(self):
        lines = LIVE.frame([("Rho", [float("nan")] * 10)], list(range(10)), 40, 5)
        self.assertTrue(lines)


if __name__ == "__main__":
    unittest.main()
