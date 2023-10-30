#!/usr/bin/env python

## \file run.py
#  \brief Unsteady adjoint FEA case with custom load.
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
SOLVER= ELASTICITY
GEOMETRIC_CONDITIONS= LARGE_DEFORMATIONS
FORMULATION_ELASTICITY_2D= PLANE_STRESS

TIME_DOMAIN= YES
TIME_STEP=0.01
TIME_DISCRE_FEA= NEWMARK_IMPLICIT
NEWMARK_GAMMA= 0.5
NEWMARK_BETA= 0.25

MATERIAL_MODEL= NEO_HOOKEAN
ELASTICITY_MODULUS= 10000
POISSON_RATIO= 0.3
MATERIAL_DENSITY= __DENSITY__

MARKER_CLAMPED= ( x_minus )
MARKER_FLUID_LOAD= ( x_plus, y_minus, y_plus )

LINEAR_SOLVER= CONJUGATE_GRADIENT
DISCADJ_LIN_SOLVER= CONJUGATE_GRADIENT
LINEAR_SOLVER_PREC= ILU
DISCADJ_LIN_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-5
LINEAR_SOLVER_ITER= 100

MESH_FORMAT= RECTANGLE
MESH_BOX_SIZE= ( 17, 5, 0 )
MESH_BOX_LENGTH= ( 0.5, __HEIGHT__, 0 )

OUTPUT_FILES= RESTART, PARAVIEW
OUTPUT_WRT_FREQ= 1
OBJECTIVE_FUNCTION= STRESS_PENALTY
STRESS_PENALTY_PARAM= ( 500, 20 )

INNER_ITER= 20
CONV_RESIDUAL_MINVAL= -4
CONV_STARTITER= 5
MAX_TIME= 0.2
TIME_ITER= 21
"""

primal_settings = """
MATH_PROBLEM= DIRECT
SCREEN_OUTPUT= TIME_ITER, CUR_TIME, INNER_ITER, RMS_RES, LINSOL, VMS, STRESS_PENALTY, TAVG_STRESS_PENALTY
HISTORY_OUTPUT= ITER, RMS_RES, STRUCT_COEFF, TAVG_STRUCT_COEFF
CONV_FIELD= REL_RMS_RTOL
"""

adjoint_settings = """
MATH_PROBLEM= DISCRETE_ADJOINT
SCREEN_OUTPUT= TIME_ITER, CUR_TIME, INNER_ITER, ADJOINT_DISP_X, ADJOINT_DISP_Y, LINSOL, SENSITIVITY
CONV_FIELD= ADJOINT_DISP_X, ADJOINT_DISP_Y
FEA_ADVANCED_MODE= YES
UNST_ADJOINT_ITER= 21
ITER_AVERAGE_OBJ= 0
SOLUTION_FILENAME= restart.dat
VOLUME_OUTPUT= SOLUTION, SENSITIVITY, SENSITIVITY_N
"""


def ApplyLoad(driver, marker_id, peak_load):
  """
  Apply a load based on the coordinates and return the derivatives
  of the nodal forces with respect to the peak load.
  """
  derivatives = []
  if marker_id < 0: return derivatives

  marker_coords = driver.MarkerCoordinates(marker_id)
  l = 0.5
  dx = l / 16 # known from mesh settings in this case.
  for i_vertex in range(driver.GetNumberMarkerNodes(marker_id)):
    x = marker_coords(i_vertex, 0)
    nodal_force = (peak_load * x / l) * dx
    # Half load due to half dx on first and last node.
    if abs(x) < 1e-6 or abs(x - l) < 1e-6:
      nodal_force = nodal_force / 2
    driver.SetMarkerCustomFEALoad(marker_id, i_vertex, (0, nodal_force))
    derivatives.append(nodal_force / peak_load)

  return derivatives


def RunPrimal(density, peak_load, height):
  """
  Runs the primal solver for a given density, peak load, and beam height.
  Returns the time average objective function.
  """
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  if rank == 0:
    with open('config_unsteady.cfg', 'w') as f:
      f.write(common_settings.replace('__DENSITY__', str(density)).replace('__HEIGHT__', str(height)) +
              primal_settings)
  comm.Barrier()

  # Initialize the primal driver of SU2, this includes solver preprocessing.
  try:
    driver = pysu2ad.CSinglezoneDriver('config_unsteady.cfg', 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2ad.CSinglezoneDriver : ', exception)
    raise

  # Get the ID of the marker where the load is applied.
  all_marker_ids = driver.GetMarkerIndices()
  marker_name = 'y_minus'
  marker_id = all_marker_ids[marker_name] if marker_name in all_marker_ids else -1

  # Apply a load based on the coordinates.
  ApplyLoad(driver, marker_id, peak_load)

  # Solve.
  driver.StartSolver()

  # Get the time average
  tavg_stress_penalty = driver.GetOutputValue('TAVG_STRESS_PENALTY')

  # Finalize the solver and exit cleanly.
  driver.Finalize()

  return tavg_stress_penalty


def RunAdjoint(density, peak_load, height):
  """
  Runs the adjoint solver and returns the sensitivity of the objective function to the peak
  load, to the material density, and to the beam height.
  """
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  if rank == 0:
    with open('config_unsteady_ad.cfg', 'w') as f:
      f.write(common_settings.replace('__DENSITY__', str(density)).replace('__HEIGHT__', str(height)) +
              adjoint_settings)
  comm.Barrier()

  # Initialize the adjoint driver of SU2, this includes solver preprocessing.
  try:
    driver = pysu2ad.CDiscAdjSinglezoneDriver('config_unsteady_ad.cfg', 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2ad.CDiscAdjSinglezoneDriver : ', exception)
    raise

  # Get the ID of the marker where the load is applied.
  all_marker_ids = driver.GetMarkerIndices()
  marker_name = 'y_minus'
  marker_id = all_marker_ids[marker_name] if marker_name in all_marker_ids else -1

  # Apply the same load that was used in the primal problem.
  derivatives = ApplyLoad(driver, marker_id, peak_load)

  n_vertex = driver.GetNumberMarkerNodes(marker_id) if marker_id >= 0 else 0
  load_sens = 0.0

  # Run the time loop in python to extract sensitivities at each step.
  for time_iter in range(driver.GetNumberTimeIter()):
    # Preprocess adjoint iteration (AD recording).
    driver.Preprocess(time_iter)

    # Run one time iteration.
    driver.Run()
    driver.Postprocess()
    driver.Update()

    # Accumulate load sensitivies (the solver doesn't accumulate
    # these for when they are used for FSI adjoints).
    for i_vertex in range(n_vertex):
      load_sens += derivatives[i_vertex] * driver.GetMarkerFEALoadSensitivity(marker_id, i_vertex)[1]

    # Monitor the solver and output solution to file if required.
    driver.Monitor(time_iter)
    driver.Output(time_iter)

  rho_sens = driver.GetOutputValue('SENS_RHO_0')

  height_sens = 0.0
  sensitivity = driver.Sensitivity(driver.GetSolverIndices()['ADJ.FEA'])
  coords = driver.Coordinates()

  for i_node in range(driver.GetNumberNodes() - driver.GetNumberHaloNodes()):
    y = coords(i_node, 1)
    dy_dh = y / height
    height_sens += dy_dh * sensitivity(i_node, 1)

  # Finalize the solver and exit cleanly.
  driver.Finalize()

  return comm.allreduce(load_sens), rho_sens, comm.allreduce(height_sens)


def main():
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  # Run the primal with 2 loads to compute the sensitivity via finite differences.
  obj_pert_height = RunPrimal(1, 2, 0.05001)
  obj_pert_load = RunPrimal(1, 2.002, 0.05)
  obj_pert_rho = RunPrimal(1.0001, 2, 0.05)
  # Run the un-perturbed last to use the restarts for the adjoint.
  obj = RunPrimal(1, 2, 0.05)
  sens_height_fd = (obj_pert_height - obj) / 0.00001
  sens_load_fd = (obj_pert_load - obj) / 0.002
  sens_rho_fd = (obj_pert_rho - obj) / 0.0001

  sens_load, sens_rho, sens_height = RunAdjoint(1, 2, 0.05)

  if rank == 0:
    print("       Finite Differences\tDiscrete Adjoint")
    print(f"Height {sens_height_fd}\t{sens_height}")
    print(f"Load   {sens_load_fd}\t{sens_load}")
    print(f"Rho    {sens_rho_fd}\t{sens_rho}")

  assert abs(sens_height / sens_height_fd - 1) < 1e-4, "Error in geometric derivatives."
  assert abs(sens_load / sens_load_fd - 1) < 1e-4, "Error in load derivative."
  assert abs(sens_rho / sens_rho_fd - 1) < 1e-3, "Error in material derivative."

  # Print results for the regression script to check.
  if rank == 0:
    print("\n")
    print(100, 100, sens_load_fd, sens_load, sens_rho_fd * 10, sens_rho * 10, sens_height_fd / 100, sens_height / 100)


if __name__ == '__main__':
  main()
