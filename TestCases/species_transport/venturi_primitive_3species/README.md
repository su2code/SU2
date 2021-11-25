# Primitive Venturi with 2 and 3 species

All .cfg's in this folder make use of the same mesh. The geometry is channel with an orthogonal duct.
The 'Venturi Primitive' is not exactly correct as there is no chokepoint behind the mixing region.
The walls are straight and the diameters are equal.
With just ~3k elements the mesh is excellent for debugging.

Note that all simulations in this folder use passive species transport, i.e. no mixture fluid materials.

- `species2_venturiPrimitive.cfg` Just 1 additional transport equation is solved, i.e. it is a 2 species flow.
In the regression test values the `SURFACE_SPECIES_VARIANCE` output is tested which is a scale for mixing quality.
Zero indicates a uniform species distribution.

- `species3_venturiPrimitive.cfg` Here, 2 additional transport equations are solved.
In the output the `SURFACE_SPECIES_0` is tested. This is simply the surface averaged species on all `MARKER_ANALYZE` and can be weighted with `AREA` or `MASSFLOW` just like the other quantities.

- `species3_venturiPrimitive_inletFile.cfg` With the `test_inlet_files.sh` a simple sanity check for inlet files is performed.
SU2 writes an `example_inlet_file.dat` when the specified inlet file is not available, with the values of the specified `MARKER_INLET` content.
Therefore comparing a simulation with this example inlet file and without inlet files should result in exactly the same results.
This is done by this script and additionally compares the history files of those 2 simulations which is supposed to be empty.
At the time of writing there is just one last digit in these history files that differs.
The case alone is unique in that it is using inlet files for 2 separate inlets and is non-dimensional and uses Turbulence.
Altogether the case is quite complex on the inlet file mechanism.
