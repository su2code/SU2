#!/usr/bin/env python

## \file run.py
#  \brief Unsteady adjoint heat transfer case with custom heat flux.
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

import pysu2ad
from mpi4py import MPI

common_settings = """
SOLVER= HEAT_EQUATION
INC_NONDIM= DIMENSIONAL

TIME_DOMAIN= YES
TIME_STEP= 0.05
TIME_MARCHING= DUAL_TIME_STEPPING-2ND_ORDER

FREESTREAM_TEMPERATURE= 345
MATERIAL_DENSITY= 2700
SPECIFIC_HEAT_CP = 870
% Very high value to make the case more interesting.
THERMAL_CONDUCTIVITY_CONSTANT= 1000

MARKER_HEATFLUX= ( x_minus, 0, x_plus, 0, y_minus, 0, y_plus, 0 )
MARKER_PYTHON_CUSTOM= ( x_minus, x_plus, y_minus, y_plus )
MARKER_MONITORING= ( x_minus, x_plus, y_minus, y_plus )

NUM_METHOD_GRAD= GREEN_GAUSS
CFL_NUMBER= 1e8

LINEAR_SOLVER= CONJUGATE_GRADIENT
DISCADJ_LIN_SOLVER= CONJUGATE_GRADIENT
LINEAR_SOLVER_PREC= ILU
DISCADJ_LIN_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-5
LINEAR_SOLVER_ITER= 100

MESH_FORMAT= RECTANGLE
MESH_BOX_SIZE= ( 33, 33, 0 )
MESH_BOX_LENGTH= ( __SIZE__, __SIZE__, 0 )

OUTPUT_FILES= RESTART, PARAVIEW
OUTPUT_WRT_FREQ= 1
OBJECTIVE_FUNCTION= AVG_TEMPERATURE

INNER_ITER= 20
CONV_RESIDUAL_MINVAL= -4
CONV_STARTITER= 2
MAX_TIME= 2.0
TIME_ITER= 41
"""

primal_settings = """
MATH_PROBLEM= DIRECT
SCREEN_OUTPUT= TIME_ITER, CUR_TIME, INNER_ITER, RMS_RES, LINSOL, AVG_TEMPERATURE, TAVG_AVG_TEMPERATURE
HISTORY_OUTPUT= ITER, RMS_RES, HEAT, TAVG_HEAT
"""

adjoint_settings = """
MATH_PROBLEM= DISCRETE_ADJOINT
SCREEN_OUTPUT= TIME_ITER, CUR_TIME, INNER_ITER, RMS_RES, LINSOL
UNST_ADJOINT_ITER= 41
% Only interested in the final value of the objective.
ITER_AVERAGE_OBJ= 1
SOLUTION_FILENAME= restart.dat
VOLUME_OUTPUT= SOLUTION, SENSITIVITY, SENSITIVITY_N
"""


def ApplyHeatFlux(time, driver, marker_ids):
  """Applies a heat flux on all boundaries for a period of time and then removes it."""
  for marker_id in marker_ids:
    if marker_id < 0:
      continue
    hf = (-1 * 2700 * 870, 0)[time > 0.2]
    for i_vertex in range(driver.GetNumberMarkerNodes(marker_id)):
      driver.SetMarkerCustomNormalHeatFlux(marker_id, i_vertex, hf)


def RunPrimal(size):
  """
  Run the heat solver with a custom heat (function of time) flux on all boundaries.
  Returns the final average boundary temperature.
  """
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  if rank == 0:
    with open('config_unsteady.cfg', 'w') as f:
      f.write(common_settings.replace('__SIZE__', str(size)) + primal_settings)
  comm.Barrier()

  # Initialize the primal driver of SU2, this includes solver preprocessing.
  try:
    driver = pysu2ad.CSinglezoneDriver('config_unsteady.cfg', 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2ad.CSinglezoneDriver : ', exception)
    raise

  # Get the ID of the markers where the heat flux is applied.
  all_marker_ids = driver.GetMarkerIndices()
  marker_names = ['x_minus', 'x_plus', 'y_minus', 'y_plus']
  marker_ids = []
  for name in marker_names:
    marker_ids.append(all_marker_ids[name] if name in all_marker_ids else -1)

  # Run the time loop in python to vary the heat flux.
  dt = driver.GetUnsteadyTimeStep()

  for time_iter in range(driver.GetNumberTimeIter()):
    # Custom heat flux.
    ApplyHeatFlux(time_iter * dt, driver, marker_ids)

    driver.Preprocess(time_iter)

    # Run one time iteration.
    driver.Run()
    driver.Postprocess()
    driver.Update()

    # Monitor the solver and output solution to file if required.
    driver.Monitor(time_iter)
    driver.Output(time_iter)

  # Get the final average temperature.
  avg_temperature = driver.GetOutputValue('AVG_TEMPERATURE')

  # Finalize the solver and exit cleanly.
  driver.Finalize()

  return avg_temperature


def RunAdjoint(size):
  """
  Runs the adjoint heat solver and returns the sensitivity of the objective function to
  size of the domain and to the initial temperature.
  """
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  if rank == 0:
    with open('config_unsteady_ad.cfg', 'w') as f:
      f.write(common_settings.replace('__SIZE__', str(size)) + adjoint_settings)
  comm.Barrier()

  # Initialize the adjoint driver of SU2, this includes solver preprocessing.
  try:
    driver = pysu2ad.CDiscAdjSinglezoneDriver('config_unsteady_ad.cfg', 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2ad.CDiscAdjSinglezoneDriver : ', exception)
    raise

  # Get the ID of the markers where the heat flux is applied.
  all_marker_ids = driver.GetMarkerIndices()
  marker_names = ['x_minus', 'x_plus', 'y_minus', 'y_plus']
  marker_ids = []
  for name in marker_names:
    marker_ids.append(all_marker_ids[name] if name in all_marker_ids else -1)

  # Run the time loop in python to vary the heat flux.
  dt = driver.GetUnsteadyTimeStep()

  # Run the time loop in python to extract sensitivities at each step.
  for time_iter in range(driver.GetNumberTimeIter()):
    # Note that time runs in reverse for the adjoint solver.
    ApplyHeatFlux((driver.GetNumberTimeIter() - time_iter - 1) * dt, driver, marker_ids)

    # Preprocess adjoint iteration (AD recording).
    driver.Preprocess(time_iter)

    # Run one time iteration.
    driver.Run()
    driver.Postprocess()
    driver.Update()

    # Monitor the solver and output solution to file if required.
    driver.Monitor(time_iter)
    driver.Output(time_iter)

  size_sens = 0.0
  temp_sens = 0.0
  sensitivity = driver.Sensitivity(driver.GetSolverIndices()['ADJ.HEAT'])
  sens_n = driver.SolutionTimeN(driver.GetSolverIndices()['ADJ.HEAT'])
  sens_n1 = driver.SolutionTimeN1(driver.GetSolverIndices()['ADJ.HEAT'])
  coords = driver.Coordinates()

  for i_node in range(driver.GetNumberNodes() - driver.GetNumberHaloNodes()):
    x, y = coords.Get(i_node)
    dx_ds = x / size
    dy_ds = y / size
    size_sens += dx_ds * sensitivity(i_node, 0) + dy_ds * sensitivity(i_node, 1)
    # For 2nd order BDF the initial conditions apply to time n and n-1,
    # hence we combine the two sensitivities. n-1 is multiplied by 2 because
    # the initial conditions are used as the n-1 solution twice (in the 1st and
    # second iterations). However, note that this is an approximation, the
    # correct value would be sens_n_t0 + sens_n1_t0 + sens_n1_t1 (we approximate
    # sens_n1 at t1 as sens_n1 at t0).
    temp_sens += sens_n(i_node, 0) + 2 * sens_n1(i_node, 0)

  # Finalize the solver and exit cleanly.
  driver.Finalize()

  return comm.allreduce(size_sens), comm.allreduce(temp_sens)


def main():
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  obj_pert_size = RunPrimal(0.100005)
  obj = RunPrimal(0.1)
  sens_size_fd = (obj_pert_size - obj) / 0.000005

  sens_size, sens_temp = RunAdjoint(0.1)

  if rank == 0:
    print("     Finite Differences\tDiscrete Adjoint")
    print(f"Size {sens_size_fd}\t{sens_size}")

  assert abs(sens_size / sens_size_fd - 1) < 1e-4, "Error in geometric derivatives."
  # We expect the final average temperature to be directly proportional to the initial
  # temperature since the applied heat flux is not a function of temperature.
  assert abs(sens_temp - 1) < 1e-5, "Error in initial condition derivatives."

  # Print results for the regression script to check.
  if rank == 0:
    print("\n------------------------------ Begin Solver -----------------------------")
    print(100, 100, sens_size_fd / 100, sens_size / 100, sens_temp)

if __name__ == '__main__':
  main()
