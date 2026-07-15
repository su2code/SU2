#!/bin/bash
# Quick test script for FIML forward mode
# Usage: ./test_fiml_forward.sh

echo "========================================="
echo "FIML Forward Mode Test"
echo "========================================="
echo ""

# Check if SU2_CFD exists
if [ ! -f "../../bin/SU2_CFD" ]; then
    echo "ERROR: SU2_CFD not found in ../../bin/"
    echo "Please build SU2 first:"
    echo "  cd ../../build && ninja"
    exit 1
fi

# Check mesh file
if [ ! -f "S809_struct_coarse.su2" ]; then
    echo "ERROR: Mesh file not found!"
    exit 1
fi

echo "Running FIML forward simulation..."
echo "Config: config_fiml_forward.cfg"
echo ""

# Run simulation
../../bin/SU2_CFD config_fiml_forward.cfg

echo ""
echo "========================================="
echo "Simulation complete!"
echo "Check restart_flow.dat and flow.vtu for results"
echo "========================================="
