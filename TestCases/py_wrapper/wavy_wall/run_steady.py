#!/usr/bin/env python

## \file run.py
#  \brief Channel with wave-like motion on walls (steady state version).
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

import pysu2ad as pysu2
from mpi4py import MPI
import numpy as np

common_settings = """
SOLVER= NAVIER_STOKES

FLUID_MODEL= IDEAL_GAS
GAMMA_VALUE= 1.4
GAS_CONSTANT= 287.87

MACH_NUMBER= 0.034
INIT_OPTION= TD_CONDITIONS
FREESTREAM_OPTION= TEMPERATURE_FS
FREESTREAM_TEMPERATURE= 299.99
FREESTREAM_PRESSURE= 98750
REYNOLDS_NUMBER= 2000
REYNOLDS_LENGTH= 0.055

VISCOSITY_MODEL= CONSTANT_VISCOSITY
MU_CONSTANT= 0.0005

MARKER_HEATFLUX= ( y_minus, 0, y_plus, 0 )
MARKER_INLET= ( x_minus, 300, 100000, 1, 0, 0 )
MARKER_FAR= ( x_plus )

MARKER_MONITORING= ( y_minus, y_plus )
MARKER_PLOTTING= ( x_minus, x_plus, y_minus, y_plus )

DEFORM_MESH= YES
MARKER_DEFORM_MESH= ( y_minus, y_plus )
DEFORM_CONSOLE_OUTPUT= YES
DEFORM_LINEAR_SOLVER= CONJUGATE_GRADIENT
DEFORM_LINEAR_SOLVER_PREC= ILU
DEFORM_LINEAR_SOLVER_ERROR= 1e-12
DEFORM_LINEAR_SOLVER_ITER= 1000
DEFORM_STIFFNESS_TYPE= CONSTANT_STIFFNESS
DEFORM_POISSONS_RATIO= 0.35

CONV_NUM_METHOD_FLOW= ROE
ENTROPY_FIX_COEFF= 1e-5
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= NONE
NUM_METHOD_GRAD= GREEN_GAUSS

MESH_FORMAT= RECTANGLE
MESH_BOX_SIZE= ( 193, 33, 0 )
MESH_BOX_LENGTH= ( 0.75, __WIDTH__, 0 )

CUSTOM_OUTPUTS= 'p_tot : Macro{PRESSURE + 0.5 * DENSITY * (pow(VELOCITY_X, 2) + pow(VELOCITY_Y, 2))};\
                 p_in : MassFlowAvg{$p_tot}[x_minus];\
                 p_out : MassFlowAvg{$p_tot}[x_plus];\
                 p_drop : Function{p_out - p_in}'

CUSTOM_OBJFUNC= 'p_drop'
OBJECTIVE_FUNCTION= CUSTOM_OBJFUNC

SOLUTION_FILENAME= restart.dat
OUTPUT_FILES= RESTART, PARAVIEW_MULTIBLOCK
OUTPUT_WRT_FREQ= 9999
SCREEN_WRT_FREQ_INNER= 10

INNER_ITER= 3000
CONV_STARTITER= 10
"""

primal_settings = """
MATH_PROBLEM= DIRECT
RESTART_SOL= YES

CFL_NUMBER= 1000
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 0.1
LINEAR_SOLVER_ITER= 20

SCREEN_OUTPUT= INNER_ITER, RMS_RES, LINSOL, p_drop
HISTORY_OUTPUT= ITER, RMS_RES, CUSTOM
VOLUME_OUTPUT= PRIMITIVE, VORTEX_IDENTIFICATION

CONV_RESIDUAL_MINVAL= -13
CONV_FIELD= RMS_DENSITY
"""

adjoint_settings = """
MATH_PROBLEM= DISCRETE_ADJOINT

CFL_NUMBER= 1000
DISCADJ_LIN_SOLVER= FGMRES
DISCADJ_LIN_PREC= ILU
LINEAR_SOLVER_ERROR= 0.001
LINEAR_SOLVER_ITER= 20

SCREEN_OUTPUT= INNER_ITER, RMS_RES, LINSOL
HISTORY_OUTPUT= ITER, RMS_RES, SENSITIVITY

QUASI_NEWTON_NUM_SAMPLES= 20

CONV_RESIDUAL_MINVAL= -8
CONV_FIELD= REL_RMS_ADJ_DENSITY
"""


def GetMarkerIds(driver):
    """Get the ID of the markers where the deformation is applied."""
    all_marker_ids = driver.GetMarkerIndices()
    marker_names = ['y_minus', 'y_plus']
    marker_ids = []
    for name in marker_names:
        marker_ids.append(all_marker_ids[name] if name in all_marker_ids else -1)
    return marker_ids


def RunPrimal(channel_width, deform_amplitude):
    """
    Runs the primal solver for a given amplitude of the deformation.
    Returns the objective function.
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        with open('config.cfg', 'w') as f:
            f.write(common_settings.replace('__WIDTH__', str(channel_width)) + primal_settings)
    comm.Barrier()

    # Initialize the corresponding driver of SU2, this includes solver preprocessing.
    try:
        driver = pysu2.CSinglezoneDriver('config.cfg', 1, comm)
    except TypeError as exception:
        print('A TypeError occured in pysu2.CSinglezoneDriver : ', exception)
        raise

    marker_ids = GetMarkerIds(driver)

    for marker_id in marker_ids:
        Deformation(0, deform_amplitude, marker_id, driver)

    driver.StartSolver()
    avg_p_drop = driver.GetOutputValue('p_drop')

    # Finalize the solver and exit cleanly
    driver.Finalize()

    return avg_p_drop


def RunAdjoint(channel_width, deform_amplitude):
    """
    Runs the adjoint solver for a given amplitude of the deformation, and channel width.
    Returns the sensitivity of the objective function to the amplitude and width.
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        with open('config_ad.cfg', 'w') as f:
            f.write(common_settings.replace('__WIDTH__', str(channel_width)) + adjoint_settings)
    comm.Barrier()

    # Initialize the corresponding driver of SU2, this includes solver preprocessing.
    try:
        driver = pysu2.CDiscAdjSinglezoneDriver('config_ad.cfg', 1, comm)
    except TypeError as exception:
        print('A TypeError occured in pysu2.CDiscAdjSinglezoneDriver : ', exception)
        raise

    marker_ids = GetMarkerIds(driver)

    # Preprocess adjoint to load the primal solution.
    driver.Preprocess(0)

    # Check the surface deformation is the expected from loading the primal results.
    for marker_id in marker_ids:
        Deformation(0, deform_amplitude, marker_id, driver, True)

    driver.Run()
    driver.Postprocess()
    driver.Update()

    # Accumulate the sensitivity to the deformation amplitude.
    sens = 0
    for marker_id in marker_ids:
        sens += SumSensitivity(0, marker_id, driver)

    driver.Monitor(0)
    driver.Output(0)

    size_sens = 0.0
    sensitivity = driver.Sensitivity(driver.GetSolverIndices()['ADJ.FLOW'])
    coords = driver.InitialCoordinates()

    for i_node in range(driver.GetNumberNodes() - driver.GetNumberHaloNodes()):
        y = coords.Get(i_node, 1)
        dy_dw = y / channel_width
        size_sens += dy_dw * sensitivity(i_node, 1)

    # Finalize the solver and exit cleanly
    driver.Finalize()

    return comm.allreduce(size_sens), comm.allreduce(sens)


def Phase(t):
    # 10 cycles per second
    return 10 * 2 * np.pi * t


def DeformationBottom(x, phase):
    # Deformation between 0.1 and 0.4.
    if x < 0.1 or x > 0.4:
        return 0
    t = (x - 0.1) / 0.3 * 2 * np.pi
    ampl = (1 - np.cos(t)) / 2

    # Shape parameter, makes the peaks narrower or wider.
    N = 1

    p = N * (t - phase)
    denom = N * 2 * np.pi
    n = round(p / denom)
    p -= n * denom
    wave = (1 - np.cos(min(max(0, p + np.pi), 2 * np.pi))) / 2

    return ampl * wave


def DeformationTop(x, phase):
    # Top deforms downward with a phase shift.
    return -DeformationBottom(x, phase + np.pi)


def Deformation(t, amplitude, marker_id, driver, check = False):
    """
    Applies the deformation to the marker or checks that the current mesh
    coordinates (loaded by the adjoint solver) match the expected for the time.
    """
    if marker_id < 0:
        return
    ini_coords = driver.MarkerInitialCoordinates(marker_id)
    coords = driver.MarkerCoordinates(marker_id)
    phase = Phase(t)

    for i_vertex in range(ini_coords.Shape()[0]):
        x, y = ini_coords.Get(i_vertex)
        if (y > 1e-3):
            dy = DeformationTop(x, phase)
        else:
            dy = DeformationBottom(x, phase)
        if not check:
            driver.SetMarkerCustomDisplacement(marker_id, i_vertex, (0.0, amplitude * dy))
        else:
            assert abs((coords(i_vertex, 1) - y) / amplitude - dy) < 1e-9


def SumSensitivity(t, marker_id, driver):
    """Integrates the sensitivity with respect to the amplitude of the deformation."""
    if marker_id < 0:
        return 0.0
    ini_coords = driver.MarkerInitialCoordinates(marker_id)
    phase = Phase(t)
    sens = 0

    for i_vertex in range(ini_coords.Shape()[0]):
        i_point = driver.GetMarkerNode(marker_id, i_vertex)
        if not driver.GetNodeDomain(i_point):
            continue

        x, y = ini_coords.Get(i_vertex)
        if (y > 1e-3):
            dy = DeformationTop(x, phase)
        else:
            dy = DeformationBottom(x, phase)
        sens_x, sens_y = driver.GetMarkerDisplacementSensitivity(marker_id, i_vertex)
        sens += sens_y * dy

    return sens


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    p_drop = RunPrimal(0.055, 0.045)
    w_sens, a_sens = RunAdjoint(0.055, 0.045)

    # Print results for the regression script to check.
    # The derivates were verified with finite difference steps of 5e-5.
    if rank == 0:
        print("\n------------------------------ Begin Solver -----------------------------\n")
        print(100, 100, p_drop / 1000, w_sens / 100000, a_sens / 100000, '\n')


if __name__ == '__main__':
  main()
