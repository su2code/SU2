#!/usr/bin/env python
"""Test script for probe performance fix following SU2's testing pattern."""

from __future__ import print_function
import sys
import os

# Add parent directory to path to import TestCase
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from TestCase import TestCase

def main():
    """Run probe performance tests."""
    
    test_list = []
    
    # Test 1: 1 probe (linear search path, ≤10)
    probe_1 = TestCase('probe_1_linear')
    probe_1.command.exec = "SU2_CFD"
    probe_1.cfg_dir = "user_defined_functions"
    probe_1.cfg_file = "test_1_probe.cfg"
    probe_1.test_iter = 5
    probe_1.test_vals = []  # Will extract from run
    probe_1.tol = 1e-6
    test_list.append(probe_1)
    
    # Test 2: 5 probes (linear search path, ≤10)
    probe_5 = TestCase('probe_5_linear')
    probe_5.command.exec = "SU2_CFD"
    probe_5.cfg_dir = "user_defined_functions"
    probe_5.cfg_file = "test_5_probes.cfg"
    probe_5.test_iter = 5
    probe_5.test_vals = []  # Will extract from run
    probe_5.tol = 1e-6
    test_list.append(probe_5)
    
    # Test 3: 11 probes (ADT path, >10)
    probe_11 = TestCase('probe_11_adt')
    probe_11.command.exec = "SU2_CFD"
    probe_11.cfg_dir = "user_defined_functions"
    probe_11.cfg_file = "test_11_probes.cfg"
    probe_11.test_iter = 5
    probe_11.test_vals = []  # Will extract from run
    probe_11.tol = 1e-6
    test_list.append(probe_11)
    
    # Test 4: 100 probes (ADT path, >10, large count)
    probe_100 = TestCase('probe_100_adt')
    probe_100.command.exec = "SU2_CFD"
    probe_100.cfg_dir = "user_defined_functions"
    probe_100.cfg_file = "test_100_probes.cfg"
    probe_100.test_iter = 5
    probe_100.test_vals = []  # Will extract from run
    probe_100.tol = 1e-6
    test_list.append(probe_100)
    
    # Run tests
    passed = True
    for test in test_list:
        test.timeout = 600
        result = test.run_test()
        
        # Check logs for probe output and print it
        log_name = os.path.splitext(test.cfg_file)[0] + "_check.log"
        log_path = os.path.join(test.cfg_dir, log_name)
        if os.path.exists(log_path):
            print(f"--- Probe Check for {test.tag} ---")
            expected_map = {
                'probe_1_linear': {'probe1': 101325},
                'probe_5_linear': {'probe1': 101343, 'probe3': 101326, 'probe5': 101323},
                'probe_11_adt': {'probe1': 101406, 'probe6': 101325, 'probe11': 100928},
                'probe_100_adt': {'probe1': 101406, 'probe51': 101436, 'probe100': 101231}
            }
            
            try:
                with open(log_path, 'r') as f:
                    content = f.read()
                    import re
                    # Find all "Probe <name>: <value>" lines
                    matches = re.findall(r"Probe (\w+): ([\d\.\-\+eE]+)", content)
                    
                    found_probes = {}
                    for name, val_str in matches:
                        found_probes[name] = float(val_str)
                        print(f"Found {name}: {found_probes[name]}")

                    # Verify against expected
                    if test.tag in expected_map:
                        for pname, pval in expected_map[test.tag].items():
                            if pname not in found_probes:
                                print(f"ERROR: Expected probe {pname} not found in output!")
                                passed = False
                            else:
                                diff = abs(found_probes[pname] - pval)
                                if diff > 1.0: # Tolerance of 1.0 Pa is reasonable for 1e5 Pa
                                    print(f"ERROR: Probe {pname} value mismatch! Expected {pval}, got {found_probes[pname]}")
                                    passed = False
                                else:
                                    print(f"Probe {pname} OK (Diff: {diff})")
                    else:
                        print(f"No expected values defined for {test.tag}")

            except Exception as e:
                print(f"Error reading/parsing log: {e}")
                passed = False
            print("-----------------------------------")
        else:
            print(f"Log file {log_path} not found!")
            passed = False
            
        if not result:
            passed = False
    
    if passed:
        print("\n✅ All probe performance tests PASSED")
        return 0
    else:
        print("\n❌ Some probe performance tests FAILED")
        return 1

if __name__ == '__main__':
    sys.exit(main())
