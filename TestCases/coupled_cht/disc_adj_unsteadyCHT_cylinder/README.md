# Unsteady CHT Adjoint Testcase

## Short Description
This is a 2D cylinder in freestream testcase. The flow is incompressible and laminar at Re=200.
A uniform vortex shedding forms behind the cylinder and each shedding cycle is resolved by 54 timesteps.
The pin is heated on the inner surface.

## Mesh
The mesh is for testing purposes only and contains ~4000 elements for the flow and ~1000 for the heat domain.
A gmsh .geo file is added such that the mesh can be recreated and modified.

## Recreating full primal
The primal for a full cycle can be restarted with the `solution_*_00000.dat` and `solution_*_00001.dat`.
The primal solution is necessary for the Discrete Adjoint sweep and for the gradient of the full
shedding cycle the full primal is necessary. The necessary changes to `chtMaster.cfg` are documented
in the config itself

## Discrete Adjoint
In the regression testcase of SU2 only 2 reverse steps are taken.
For that, the solution files 52-55 for the adjoint are added.
The objective Function is the average temperature on the inner pin surface, averaged over the full time.

## Finite Differences using FADO
In order to validate the Discrete Adjoint gradient a Finite Differences python script `findiff.py`
using [FADO](www.github.com/su2code/FADO) is added.
The script starts a baseline simulation and an additional for each of the 18 Design Variables.
Using the `postprocess.py` to print the absolute difference and relative difference in percent to screen.
Make sure to adapt the `chtMaster.cfg` as for the primal.
