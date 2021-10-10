# Periodic Pins

This a solid heat conduction testcase in which the `HEAT_EQUATION` solver runs standalone (i.e. not as CHT).
The simulation domain is the solid domain of the `incomp_navierstokes/streamwise_periodic/chtPinArray_2d`-Testcase.
Therefore the provided gmsh `.geo` file contains the full CHT mesh but only writes out the solid zone when called.

Note that using periodic boundary conditions for the solid zone made the solution take ~10x more iterations to converge , compared to the same setup using adiabatic walls.
This was found for solid only as well as CHT cases.

Expected results are perfectly matched Temperatures at the periodic interface. Compare e.g. using Paraview's `Transform`-Filter with domain length 0.0111544m.
