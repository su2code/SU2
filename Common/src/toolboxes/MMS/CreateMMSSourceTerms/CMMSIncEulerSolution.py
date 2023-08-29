#!/usr/bin/env python

## \file CMMSIncEulerSolution.py
#  \brief Python script that generates the source terms for a
#         manufactured solution for the incompressible Euler eqns.
#  \author T. Economon
#  \version 8.0.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

from __future__ import division
from sympy import *
from sympy.printing import print_ccode

init_printing()

x = Symbol("x", real=True)
y = Symbol("y", real=True)

u_0 = Symbol("u_0", real=True)
v_0 = Symbol("v_0", real=True)
P_0 = Symbol("P_0", real=True)

Density = Symbol("Density", real=True)
Viscosity = Symbol("Viscosity", real=True)
epsilon = Symbol("epsilon", real=True)

# Set the manufactured solution. The solution
# is from Salari K, and Knupp P, "Code verification
# by the method of manufactured solutions,"
# SAND 2000-1444, Sandia National Laboratories,
# Albuquerque, NM, 2000.

u = u_0 * (sin(x * x + y * y) + epsilon)
v = v_0 * (cos(x * x + y * y) + epsilon)
P = P_0 * (sin(x * x + y * y) + 2.0)

# Put the manufactured solution through the Euler
# equations to get the source term for each eqn.

print("\nSource-Rho: ")
Sp = diff(Density * u, x) + diff(Density * v, y)
print_ccode(Sp.simplify())

print("\nSource-U: ")
Su = diff(Density * u * u + P, x) + diff(Density * u * v, y)
print_ccode(Su.simplify())

print("\nSource-V: ")
Sv = diff(Density * u * v, x) + diff(Density * v * v + P, y)
print_ccode(Sv.simplify())
