#!/usr/bin/env python

## \file run.py
#  \brief Unsteady adjoint heat transfer case with custom heat flux.
#  \version 8.3.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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
import sys
import pysu2			            # imports the SU2 wrapped module

common_settings = """
SOLVER= HEAT_EQUATION
INC_NONDIM= DIMENSIONAL

TIME_DOMAIN= YES
TIME_STEP= 0.002
TIME_MARCHING= DUAL_TIME_STEPPING-2ND_ORDER

FREESTREAM_TEMPERATURE= 300
MATERIAL_DENSITY= 1000
SPECIFIC_HEAT_CP = 500
THERMAL_CONDUCTIVITY_CONSTANT= 10.0

MARKER_HEATFLUX= ( x_minus, 0, x_plus, 0, y_minus, 0, y_plus, 0 )
MARKER_PYTHON_CUSTOM= ( x_minus, x_plus, y_minus, y_plus )
MARKER_MONITORING= ( x_minus, x_plus, y_minus, y_plus )

PYTHON_CUSTOM_SOURCE= YES

NUM_METHOD_GRAD= GREEN_GAUSS
CFL_NUMBER= 100

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
OUTPUT_WRT_FREQ= 100,100
OBJECTIVE_FUNCTION= AVG_TEMPERATURE

INNER_ITER= 20
CONV_RESIDUAL_MINVAL= -4
CONV_STARTITER= 2
MAX_TIME= 10.0
TIME_ITER= 101
"""

primal_settings = """
MATH_PROBLEM= DIRECT
SCREEN_OUTPUT= TIME_ITER, CUR_TIME, INNER_ITER, RMS_RES, LINSOL, AVG_TEMPERATURE, TAVG_AVG_TEMPERATURE
HISTORY_OUTPUT= ITER, RMS_RES, HEAT, TAVG_HEAT
"""

def HeatSource(time, driver, iPoint):
  """Applies a source in the lower left quadrant for a period of time."""
  allCoords = driver.Coordinates()
  coord = allCoords.Get(iPoint)
  x = coord[0]
  y = coord[1]
  xp = 0.0125
  yp = 0.0125
  R = 0.0125
  Source = 0.0
  if (time > 0.0):
      if ((x-xp)*(x-xp) + (y-yp)*(y-yp) < R*R):
          Source = 100.0

  return Source

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
    print("\n------------------------------ Begin Solver -----------------------------\n")
  sys.stdout.flush()

  comm.Barrier()

 # Initialize the corresponding driver of SU2, this includes solver preprocessing
  try:
    driver = pysu2.CSinglezoneDriver("config_unsteady.cfg", 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ',exception)
    raise

  # Run the time loop in python to vary the heat flux.
  dt = driver.GetUnsteadyTimeStep()

  iHEATSOLVER = driver.GetSolverIndices()['HEAT']
  Source = driver.UserDefinedSource(iHEATSOLVER)
  

  for time_iter in range(driver.GetNumberTimeIter()):
    # Custom heat flux.
    #ApplyHeatFlux(time_iter * dt, driver, marker_ids)

    # set the source term, per point
    for i_node in range(driver.GetNumberNodes() - driver.GetNumberHaloNodes()):
      # add source term:
      S = HeatSource(time_iter * dt, driver, i_node)
      Source.Set(i_node, 0, S)


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


def main():

  RunPrimal(0.1)

if __name__ == '__main__':
  main()
