# Fix for Issue #1604 - Tutorial Results Mismatch

## Problem Analysis
1. Unsteady pitching NACA0012 test case did not exist
2. Tutorial on website referenced non-existent case
3. Config parameters caused convergence issues

## Solution Implemented

### Created Missing Test Case
- Location: `TestCases/unsteady/pitching_naca0012/`
- Based on similar pitching airfoil case structure
- Updated parameters per forum solutions

### Fixed Parameters
- TIME_STEP = 0.005 (matches website)
- CFL_NUMBER = 1e12
- LINEAR_SOLVER_ERROR = 0.1
- LINEAR_SOLVER_ITER = 10

### Updated Laminar Cylinder
- Corrected non-dimensionalization settings
- Ensures reproducible tutorial results

## Verification
Test cases now match tutorial intent and forum-verified solutions.

## References
- Forum: https://www.cfd-online.com/Forums/su2/242253-tutorial-cases-naca0012-up-date.html
- Issue: #1604
