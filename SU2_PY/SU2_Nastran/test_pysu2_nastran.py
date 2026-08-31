#!/usr/bin/env python

## \file test_pysu2_nastran.py
#  \brief Tests for the bulk data parsing of the Nastran structural solver.
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

from pysu2_nastran import Solver, nastran_float

# All the spellings of the real number seven listed in the Quick Reference
# Guide, "Format of Bulk Data Entries", plus the "D" exponent and the lower
# case forms accepted by the bulk data readers.
SEVEN = [
    "7.0",
    ".7E1",
    "0.7+1",
    ".70+1",
    "7.E+0",
    "70.-1",
    ".7e1",
    "700.e-2",
    "7.0D0",
    ".7d1",
    "70.-01",
    "7000.-3",
]


def justify(value, style):
    """Places a value in an eight character field, as Nastran allows."""

    if style == "left":
        return value.ljust(8)
    if style == "right":
        return value.rjust(8)
    return value.center(8)


def card(fields, style, prefix=30):
    """Builds one line of a small field bulk data echo."""

    line = fields[0].ljust(8)
    for field in fields[1:]:
        line += justify(field, style)
    return " " * prefix + line.rstrip() + "\n"


class TestNastranFloat(unittest.TestCase):
    def test_spellings_of_seven(self):
        for value in SEVEN:
            for style in ("left", "right", "centre"):
                field = justify(value, style)
                self.assertAlmostEqual(nastran_float(field), 7.0, msg=repr(field))

    def test_negative_spellings_of_seven(self):
        for value in SEVEN:
            for style in ("left", "right", "centre"):
                field = justify("-" + value, style)
                self.assertAlmostEqual(nastran_float(field), -7.0, msg=repr(field))

    def test_leading_plus_is_not_an_exponent(self):
        for style in ("left", "right", "centre"):
            self.assertAlmostEqual(nastran_float(justify("+7.0", style)), 7.0)

    def test_omitted_exponent_letter(self):
        self.assertAlmostEqual(nastran_float("1.23-5"), 1.23e-5)
        self.assertAlmostEqual(nastran_float("  -1.23-5"), -1.23e-5)
        self.assertAlmostEqual(nastran_float("-1.23+5"), -1.23e5)
        self.assertAlmostEqual(nastran_float("     1+5"), 1.0e5)

    def test_invalid_field_is_rejected(self):
        for field in ("", "        ", "abc", "1.2.3"):
            with self.assertRaises(ValueError):
                nastran_float(field)


class TestReadNastranMesh(unittest.TestCase):
    """
    Reads the same model written with different, equally valid, spellings.
    The parsed geometry must not depend on how the fields were written.
    """

    def mesh(self, style, omit_optional_fields=False):
        cd = [] if omit_optional_fields else ["0"]
        lines = [
            card(
                ["CORD2R", "1", "0", "-1.5", "-2.5", "-3.5", "-1.5", "-2.5", "-0.5"],
                style,
            ),
            card(["+", "0.5", "-2.5", "-3.5"], style),
            card(["GRID", "1", "0", "-1.25", "-2.5", "3.75"] + cd, style),
            # A blank CP field means the basic coordinate system.
            card(["GRID", "2", "", "-1.5-2", "2.5-2", "-3.5+1"] + cd, style),
            card(["GRID", "3", "1", "1.0", "2.0", "-3.0"] + cd, style),
            card(["SET1", "1", "1", "2", "3"], style),
        ]
        handle, path = tempfile.mkstemp(suffix=".f06")
        with os.fdopen(handle, "w") as mesh_file:
            mesh_file.writelines(lines)
        self.addCleanup(os.remove, path)
        return path

    def read(self, style, omit_optional_fields=False):
        solver = Solver.__new__(Solver)
        solver.Mesh_file = self.mesh(style, omit_optional_fields)
        solver.FSI_marker = "1"
        solver.node = []
        solver.markers = {}
        solver.refsystems = []
        solver._Solver__readNastranMesh()
        return solver

    def coordinates(self, solver):
        return np.array([point.GetCoord0().ravel() for point in solver.node])

    def test_justification_does_not_change_the_model(self):
        reference = self.coordinates(self.read("left"))
        self.assertEqual(reference.shape, (3, 3))
        for style in ("right", "centre"):
            np.testing.assert_allclose(self.coordinates(self.read(style)), reference)

    def test_coordinates_are_read_correctly(self):
        for style in ("left", "right", "centre"):
            coordinates = self.coordinates(self.read(style))
            # Point 1 is given in the basic system.
            np.testing.assert_allclose(coordinates[0], [-1.25, -2.5, 3.75])
            # Point 2 uses the omitted exponent letter.
            np.testing.assert_allclose(coordinates[1], [-0.015, 0.025, -35.0])
            # Point 3 is given in the reference system defined by the CORD2R,
            # whose origin is (-1.5, -2.5, -3.5) and whose x axis points to
            # (0.5, -2.5, -3.5), i.e. the basic x axis.
            np.testing.assert_allclose(coordinates[2], [-0.5, -0.5, -6.5])

    def test_blank_reference_system_field_is_the_basic_system(self):
        for style in ("left", "right", "centre"):
            solver = self.read(style)
            self.assertEqual([point.GetCP() for point in solver.node], [0, 0, 1])

    def test_optional_trailing_field_may_be_missing(self):
        for style in ("left", "right", "centre"):
            reference = self.coordinates(self.read(style))
            without_cd = self.coordinates(self.read(style, omit_optional_fields=True))
            np.testing.assert_allclose(without_cd, reference)

    def test_reference_system_is_read_correctly(self):
        for style in ("left", "right", "centre"):
            solver = self.read(style)
            self.assertEqual(len(solver.refsystems), 1)
            system = solver.refsystems[0]
            self.assertEqual(system.GetCID(), 1)
            np.testing.assert_allclose(system.GetOrigin().ravel(), [-1.5, -2.5, -3.5])
            np.testing.assert_allclose(system.GetRotMatrix(), np.eye(3), atol=1e-12)


if __name__ == "__main__":
    unittest.main()
