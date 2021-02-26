# Streamwise Periodicity testcases

This folder contains the additional Testcases for streamwise periodic flow.
A Tutorial can be found on the SU2 website.
For all Testcases a gmsh .geo file is provided which allows to recreate/modify the mesh.

## `pipe_slice_3d`

Hagen Poiseuille flow through a 1-primal-cell thick pipe slice in 3D.

Analytical solution of the velocity magnitude for steady laminar pipe flow in a round pipe `v_mag (r) = -1/(4*mu) * (Delta p / Delta x) * (R**2 - r**2)` therefore a pressure drop Delta p is prescribed.

`Re = rho * v * L / mu = 1.0 * 0.6 * 5e-3 / 1.8e-5` makes Re=167, with the critical Reynolds number being Re~=2300.

This testcase is a regression test.

## `chtPinArray_2d`

Extension of the tutorial case to a CHT problem with 1 additional solid zone.
A gradient validation between discrete and finite differences for this setup is described in the README of that folder.

This gradient validation is also part of the regression tests.

## `chtPinArray_3d`

Extension of the `chtPinArray_2d` to the 3rd dimension with again one solid zone.
The mesh provided is coarse to keep the filesize and computation time low, but using the gmsh .geo script much higher mesh resolutions can be created.

This primal simulation is part of the regression tests.
