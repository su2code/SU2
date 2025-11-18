# Test case for Equivalent area calculation 

## Introduction
The purpose of this test case is to confirm equivalent area distribution and its adjoint can be calculated by SU2. Equivalent area is a distribution calculated from pressure along a line parallel to the aircraft axis. It can be used to supersonic aircraft design. See the paper below for an example of use of SU2 for supersonic aircraft design using equivalent area distribution.

Palacios, F., Alonso, J. J., Colonno, M., Hicken, J., and Lukaczyk, T., "Adjoint-based method for supersonic aircraft design using equivalent area distribution," AIAA 2012-0269, 2012. DOI: [10.2514/6.2012-269](https://arc.aiaa.org/doi/10.2514/6.2012-269)

## Expected outcome
- Direct solver
  - Equivalent area is calculated without errors.
- Discrete adjoint solver
  - Adjoint solution is calculated without errors.

## Files for this test case
Below is a list of files for this test case and explanation for each file.
- SU2 repository
  - NACA64206.cfg
    - Configuration file.
- TestCases repository
  - NACA64206_FFD.su2
    - Mesh file. NACA64-206 was used as an airfoil. Arbitrary swept angle and taper ratio were used to create 3D geometry. It has a circumferential nearfield boundary at 3 aircraft lengths.
  - solution_flow.csv
    - Solution file of direct solver after 4000 iterations. This is required for running adjoint solver.
  - TargetEA.dat
    - Target equivalent area distribution. See the file for format. In this file, equivalent area is increased by 5% from direct solution to use as an example target.

## How to use equivalent area calculation function

### Config file
See SUPERSONIC SIMULATION section of config_template.cfg.

### Mesh
A mesh has to have a circumferential boundary around the aircraft axis within the calculation domain and has to be labeled as MARKER_NEARFIELD. This boundary has to have a structured grid with the same number of nodes along each azimuthal angle.

### Target equivalent area
TargetEA.dat is required for shape optimization using equivalent area. After running direct solver, Equivalent_Area.dat is output, which can be used to define a target since it has the same format as TargetEA.dat. Equivalent area is calculated at each node on a surface with MARKER_NEARFIELD. Difference from target is calculated by sum of squared difference at each node.

## Implementation on SU2
Equivalent area calculation function was temporarily unavailable after introduction of ver. 7.0.0. See [PR1329](https://github.com/su2code/SU2/pull/1329) on Github for recovery.