#!/usr/bin/env python

## \file run.py
#  \brief Steady solid heat conduction with a volumetric heat source in an ellipsoid.
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

from mpi4py import MPI
import math
import sys
import pysu2  # imports the SU2 wrapped module

# Slab of width L with both x faces at T_WALL and insulated y faces.
L = 0.016
H = 0.004
N_X = 17
T_WALL = 300.0
DENSITY = 5000.0
CP = 500.0
CONDUCTIVITY = 15.0

# Volumetric heat source (W/m^3) inside an ellipsoid that can be rotated about the z axis.
# The ellipsoid encloses the whole slab here, so the source is uniform.
SOURCE_VAL = 1e6
SOURCE_CENTER = (0.5 * L, 0.5 * H, 0.0)
SOURCE_AXES = (0.02, 0.02, 1.0)
SOURCE_ROTATION_Z = 0.0  # degrees

settings = f"""
SOLVER= HEAT_EQUATION
INC_NONDIM= DIMENSIONAL

FREESTREAM_TEMPERATURE= {T_WALL}
MATERIAL_DENSITY= {DENSITY}
SPECIFIC_HEAT_CP= {CP}
THERMAL_CONDUCTIVITY_CONSTANT= {CONDUCTIVITY}

MARKER_ISOTHERMAL= ( x_minus, {T_WALL}, x_plus, {T_WALL} )
MARKER_HEATFLUX= ( y_minus, 0, y_plus, 0 )
MARKER_MONITORING= ( x_minus, x_plus )

PYTHON_CUSTOM_SOURCE= YES

NUM_METHOD_GRAD= GREEN_GAUSS
CFL_NUMBER= 1e8
TIME_DISCRE_HEAT= EULER_IMPLICIT

LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-15
LINEAR_SOLVER_ITER= 10

MESH_FORMAT= RECTANGLE
MESH_BOX_SIZE= ( {N_X}, 5, 0 )
MESH_BOX_LENGTH= ( {L}, {H}, 0 )

ITER= 100
CONV_FIELD= RMS_TEMPERATURE
CONV_RESIDUAL_MINVAL= -14
CONV_STARTITER= 5

SCREEN_OUTPUT= INNER_ITER, RMS_TEMPERATURE, MAX_TEMPERATURE, TOTAL_HEATFLUX, LINSOL_ITER
HISTORY_OUTPUT= ITER, RMS_RES, HEAT, LINSOL
OUTPUT_FILES= RESTART, PARAVIEW
"""


def HeatSource(coord):
    """Returns the source (W/m^3) at a point, SOURCE_VAL inside the ellipsoid and zero outside."""
    alpha = math.radians(SOURCE_ROTATION_Z)
    dx = coord[0] - SOURCE_CENTER[0]
    dy = coord[1] - SOURCE_CENTER[1]
    x = dx * math.cos(alpha) + dy * math.sin(alpha)
    y = -dx * math.sin(alpha) + dy * math.cos(alpha)
    check = (x / SOURCE_AXES[0])**2 + (y / SOURCE_AXES[1])**2
    if len(coord) == 3:
        check += ((coord[2] - SOURCE_CENTER[2]) / SOURCE_AXES[2])**2
    return SOURCE_VAL if check <= 1 else 0.0


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        with open('config.cfg', 'w') as f:
            f.write(settings)
    comm.Barrier()

    # Initialize the corresponding driver of SU2, this includes solver preprocessing.
    try:
        driver = pysu2.CSinglezoneDriver('config.cfg', 1, comm)
    except TypeError as exception:
        print('A TypeError occured in pysu2.CDriver : ', exception)
        raise

    iHEATSOLVER = driver.GetSolverIndices()['HEAT']
    nNode = driver.GetNumberNodes() - driver.GetNumberHaloNodes()
    Coords = driver.Coordinates()

    # The heat solver works with the temperature equation divided by rho*cp,
    # so a source in W/m^3 has to be divided by rho*cp as well.
    Source = driver.UserDefinedSource(iHEATSOLVER)
    for iPoint in range(nNode):
        Source.Set(iPoint, 0, HeatSource(Coords.Get(iPoint)) / (DENSITY * CP))

    if rank == 0:
        print("\n------------------------------ Begin Solver -----------------------------\n")
    sys.stdout.flush()

    driver.StartSolver()

    # Compare with the exact solution T(x) = T_WALL + q*x*(L-x)/(2*k).
    Solution = driver.Solution(iHEATSOLVER)
    max_error = 0.0
    for iPoint in range(nNode):
        x = Coords.Get(iPoint)[0]
        exact = T_WALL + SOURCE_VAL * x * (L - x) / (2 * CONDUCTIVITY)
        max_error = max(max_error, abs(Solution.Get(iPoint)[0] - exact))
    max_error = comm.allreduce(max_error, op=MPI.MAX)

    # All the heat added by the source has to leave through the isothermal walls.
    heat_flux = driver.GetOutputValue('TOTAL_HEATFLUX')

    driver.Finalize()

    # The isothermal wall is imposed weakly, which shifts the solution by q*h^2/(2*k).
    h = L / (N_X - 1)
    tolerance = 1.01 * SOURCE_VAL * h**2 / (2 * CONDUCTIVITY)
    if rank == 0:
        print(f"Maximum temperature error: {max_error:.6f} K (tolerance {tolerance:.6f} K)")
        print(f"Total heat flux: {heat_flux:.6f} W (exact {-SOURCE_VAL * L * H:.6f} W)")

    assert max_error < tolerance, f"Test FAILED, maximum temperature error = {max_error}"
    assert abs(heat_flux / (-SOURCE_VAL * L * H) - 1) < 1e-6, f"Test FAILED, total heat flux = {heat_flux}"


if __name__ == '__main__':
    main()
