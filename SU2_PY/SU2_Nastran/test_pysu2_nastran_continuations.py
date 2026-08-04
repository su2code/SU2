#!/usr/bin/env python

## \file test_pysu2_nastran_continuations.py
#  \brief Tests for continuation and page break handling in the Nastran mesh reader.
#  \version 8.5.0 "Harrier"
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

import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pysu2_nastran import Solver

# Page header block interrupting the sorted bulk data echo, as printed by
# Nastran at every page break (structure taken from the f06 in issue #2313).
PAGE_HEADER = [
    "1    NX NASTRAN MODES ANALYSIS SET                                         "
    "NOVEMBER  15, 2022  MSC Nastran 11/19/16   PAGE  4257",
    " " * 132,
    "0" + " " * 131,
    " " * 50 + "S O R T E D   B U L K   D A T A   E C H O",
    " " * 17 + "ENTRY",
    " " * 17 + "COUNT        .   1  ..   2  ..   3  ..   4  ..   5  ..   6  "
    "..   7  ..   8  ..   9  ..  10  .",
]


def echo(count, *fields):
    """Builds one numbered line of the sorted bulk data echo."""

    line = ("%6d-" % count).rjust(22) + " " * 8
    for field in fields:
        line += str(field).ljust(8)
    return line.rstrip()


def grid_lines(ids, start=1):
    return [
        echo(start + i, "GRID", gid, 0, "1.0", "2.0", "3.0", 0)
        for i, gid in enumerate(ids)
    ]


def set1_lines(ids, start, markers):
    """Echoes SET1 1 with the given continuation markers between chunks."""

    chunks = [ids[:7]] + [ids[i : i + 8] for i in range(7, len(ids), 8)]
    lines = []
    for i, chunk in enumerate(chunks):
        head = ["SET1", 1] if i == 0 else [markers[i - 1]]
        tail = [markers[i]] if i < len(chunks) - 1 else []
        lines.append(echo(start + i, *(head + list(chunk) + tail)))
    return lines


class TestContinuations(unittest.TestCase):
    def read(self, lines, marker="1"):
        handle, path = tempfile.mkstemp(suffix=".f06")
        with os.fdopen(handle, "w") as mesh_file:
            mesh_file.write("\n".join(lines) + "\n")
        self.addCleanup(os.remove, path)
        solver = Solver.__new__(Solver)
        solver.Mesh_file = path
        solver.FSI_marker = marker
        solver.node = []
        solver.markers = {}
        solver.refsystems = []
        solver._Solver__readNastranMesh()
        return solver

    def marker_ids(self, solver, marker="1"):
        return sorted(solver.node[i].GetID() for i in solver.markers[marker])

    def test_page_break_inside_set1_continuations(self):
        # The scenario of issue #2313: headers interrupt the continuations.
        set1 = set1_lines(list(range(1, 24)), 25, ["+", "+", "+"])
        lines = grid_lines(range(1, 24)) + set1[:2] + PAGE_HEADER + set1[2:]
        solver = self.read(lines)
        self.assertEqual(self.marker_ids(solver), list(range(1, 24)))

    def test_page_break_between_set1_parent_and_continuation(self):
        set1 = set1_lines(list(range(1, 16)), 25, ["+"])
        lines = grid_lines(range(1, 16)) + set1[:1] + PAGE_HEADER + set1[1:]
        solver = self.read(lines)
        self.assertEqual(self.marker_ids(solver), list(range(1, 16)))

    def test_page_break_between_cord2r_lines(self):
        lines = (
            [echo(1, "CORD2R", 1, 0, "0.", "0.", "0.", "0.", "0.", "1.", "+")]
            + PAGE_HEADER
            + [echo(2, "+", "1.", "0.", "0.")]
            + grid_lines([1], start=3)
            + [echo(4, "SET1", 1, 1)]
        )
        solver = self.read(lines)
        self.assertEqual(solver.nRefSys, 1)
        np.testing.assert_allclose(
            solver.refsystems[0].GetRotMatrix(), np.eye(3), atol=1e-12
        )
        self.assertEqual(self.marker_ids(solver), [1])

    def test_generated_continuation_markers(self):
        # Nastran generates +000001 style markers for blank continuations,
        # they must not be read as the grid point 1.
        markers = ["+000001", "+000002"]
        lines = grid_lines(range(1, 17)) + set1_lines(list(range(1, 17)), 20, markers)
        solver = self.read(lines)
        self.assertEqual(self.marker_ids(solver), list(range(1, 17)))

    def test_user_named_continuation_markers(self):
        lines = grid_lines(range(1, 17)) + set1_lines(
            list(range(1, 17)), 20, ["+PB1", "+PB2"]
        )
        solver = self.read(lines)
        self.assertEqual(self.marker_ids(solver), list(range(1, 17)))

    def test_thru(self):
        # SET1 6 29 32 THRU 50 61 THRU 70, the second example of the QRG entry.
        ids = [29] + list(range(32, 51)) + list(range(61, 71))
        lines = grid_lines(ids) + [
            echo(40, "SET1", 6, 29, 32, "THRU", 50, 61, "THRU", 70)
        ]
        solver = self.read(lines, marker="6")
        self.assertEqual(self.marker_ids(solver, "6"), sorted(ids))

    def test_missing_continuation_is_reported(self):
        lines = grid_lines(range(1, 9)) + [
            echo(10, "SET1", 1, 1, 2, 3, 4, 5, 6, 7, "+"),
            echo(11, "SPC1", 2, 123456, 8),
        ]
        with self.assertRaisesRegex(Exception, "continuation of SET1 1.*SPC1"):
            self.read(lines)

    def test_truncated_file_is_reported(self):
        lines = grid_lines(range(1, 9)) + [
            echo(10, "SET1", 1, 1, 2, 3, 4, 5, 6, 7, "+")
        ]
        with self.assertRaisesRegex(Exception, "end of .*continuation of SET1 1"):
            self.read(lines)

    def test_point_missing_from_mesh_is_reported(self):
        lines = grid_lines(range(1, 4)) + [echo(5, "SET1", 1, 1, 2, 99)]
        with self.assertRaisesRegex(Exception, "Point 99 in the set 1 was not found"):
            self.read(lines)

    def test_unpaginated_echo_still_reads(self):
        # The tutorial shaped file must parse exactly as before.
        lines = grid_lines(range(1, 24)) + set1_lines(
            list(range(1, 24)), 25, ["+", "+", "+"]
        )
        solver = self.read(lines)
        self.assertEqual(solver.nPoint, 23)
        self.assertEqual(solver.nMarker, 1)
        self.assertEqual(self.marker_ids(solver), list(range(1, 24)))


if __name__ == "__main__":
    unittest.main()
