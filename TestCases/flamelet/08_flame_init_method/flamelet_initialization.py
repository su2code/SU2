#!/usr/bin/env python

## \file run.py
#  \brief Flame simulation initialization test case.
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

import pysu2 
from mpi4py import MPI 

# Generatl config options
config_options = """
FILENAMES_INTERPOLATOR= (MLP_Group1.mlp,MLP_Group2.mlp,MLP_Group3.mlp,MLP_Group4.mlp,MLP_Group5.mlp,MLP_Group6.mlp)

INC_VELOCITY_INIT= (1e-5, 0.0, 0.0)
MARKER_SYM=(x_minus, x_plus, y_plus, y_minus)

SOLVER = INC_NAVIER_STOKES
INC_DENSITY_MODEL= VARIABLE
INC_ENERGY_EQUATION = YES
INC_TEMPERATURE_INIT= 300
INC_NONDIM= DIMENSIONAL


FLUID_MODEL= FLUID_FLAMELET
INTERPOLATION_METHOD= MLP

CONTROLLING_VARIABLE_NAMES=(ProgressVariable,EnthalpyTot,MixtureFraction)
CONTROLLING_VARIABLE_SOURCE_NAMES=(ProdRateTot_PV,NULL,NULL)

KIND_SCALAR_MODEL= FLAMELET      
DIFFUSIVITY_MODEL= FLAMELET
VISCOSITY_MODEL= FLAMELET
CONDUCTIVITY_MODEL= FLAMELET

PREFERENTIAL_DIFFUSION=YES

FLAME_INIT_METHOD= __FLAME_INIT_METHOD__
SPARK_INIT=(5e-4, 2.5e-4, 0.0, 1e-4, 1, 50)
SPARK_REACTION_RATES=(1.5e4, 0,0)
FLAME_INIT= (4e-4, 0.00, 0.00, 1.0, 0.0, 0.0, 5e-5, 1.0)
FLAME_INIT_IGNITION= __IGNITION_TEMP__

CONV_NUM_METHOD_SPECIES= BOUNDED_SCALAR

MUSCL_SPECIES= YES
TIME_DISCRE_SPECIES= EULER_IMPLICIT
SPECIES_CLIPPING= NO
SPECIES_INIT=(-2.251135e-01,+2.262913e+03,+1.446752e-02)

CFL_REDUCTION_SPECIES= 1.0

INC_OUTLET_TYPE= PRESSURE_OUTLET

NUM_METHOD_GRAD= GREEN_GAUSS
CFL_NUMBER= 50
CFL_ADAPT= NO
ITER= __N_ITER__

LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= LU_SGS
LINEAR_SOLVER_ERROR= 1E-03
LINEAR_SOLVER_ITER= 15

CONV_NUM_METHOD_FLOW= FDS
MUSCL_FLOW= YES
TIME_DISCRE_FLOW= EULER_IMPLICIT

SCREEN_OUTPUT = CUSTOM
HISTORY_OUTPUT = CUSTOM
VOLUME_OUTPUT = SOLUTION

MESH_FORMAT= RECTANGLE
MESH_BOX_SIZE= ( 200, 50, 0 )
MESH_BOX_LENGTH= ( 1e-3, 5e-4, 0 )

OUTPUT_FILES=(RESTART)
CUSTOM_OUTPUTS= 'Tprobe : Probe{TEMPERATURE}[5e-4, 2.5e-4];'

"""

comm = MPI.COMM_WORLD

def main():
  
    comm = MPI.COMM_WORLD
    
    T_ignition = 1600
    if comm.Get_rank() == 0:
        
        # Prepare configuration files for flame front and spark ignition methods.
        config_flame_front = config_options.replace("__FLAME_INIT_METHOD__", "FLAME_FRONT")
        config_flame_front = config_flame_front.replace("__N_ITER__", "1")
        config_flame_front = config_flame_front.replace("__IGNITION_TEMP__", "%f" % T_ignition)

        config_spark = config_options.replace("__FLAME_INIT_METHOD__", "SPARK")
        config_spark = config_spark.replace("__N_ITER__", "20")
        config_spark = config_spark.replace("__IGNITION_TEMP__", "%f" % T_ignition)

        with open("config_flame_front.cfg", "w") as f:
            f.write(config_flame_front)
        with open("config_spark.cfg", "w") as f:
            f.write(config_spark)
    comm.Barrier()

    # Initialize simulation with flame front ignition
    driver_flame_front = pysu2.CSinglezoneDriver("config_flame_front.cfg", 1, comm)
    driver_flame_front.StartSolver()

    # Retrieve temperature from the burnt zone of the initial flame front.
    T_flame_front = driver_flame_front.GetOutputValue("Tprobe")
    driver_flame_front.Finalize()

    # Initialize simulation with spark ignition.
    driver_spark = pysu2.CSinglezoneDriver("config_spark.cfg", 1, comm)
    driver_spark.StartSolver()

    # Retrieve temperature from center of the spark.
    T_spark = driver_spark.GetOutputValue("Tprobe")
    driver_spark.Finalize()
    if comm.Get_rank() == 0:
        print("| %i | %.2f | %.2f|" % (1, min(T_ignition, T_flame_front), min(T_ignition, T_spark)))

if __name__ == '__main__':
  main()
