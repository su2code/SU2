# Primitive Venturi with 2 and 3 species

Additional explanattion goes here.

For now this is just a proof-of-concept for 3 species (i.e. 2 transport equations).

The inlet file test should be made non-dimensional

A good sanity check for the surface quantities is to set inlets of the same size to the same vel and just containing 1 species -> then at the outlet the avg should be 0.5 
(Check via Paraview using the Integrate filter)
