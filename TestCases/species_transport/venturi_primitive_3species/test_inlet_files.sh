#! /bin/bash
#
# T. Kattmann, 2021.11.13, Script to test inlet files

# Specify location of SU2-binaries, and other quantities
cfg="species3_primitiveVenturi.cfg"
inlet_cfg="species3_primitiveVenturi_inletFile.cfg"
inlet_file="inlet_venturi.dat"

# Create fresh inlet-file that acts the same as the non-inlet file boundary condition
rm $inlet_file
SU2_CFD $inlet_cfg | tee CFD_inlet.log
mv example_${inlet_file} $inlet_file

# Run the two simulations that should now produce identical results
SU2_CFD $inlet_cfg | tee CFD_inlet.log
SU2_CFD $cfg | tee CFD.log

# Therefore the diff on the history should be empty.
# Note: for SST non-dim some error accumulation was already seen
echo "--> diff hist*"
diff hist*
