# Unsteady Pitching NACA0012

Test case for unsteady compressible flow over pitching NACA0012 airfoil.

## Configuration
- **Solver:** Euler
- **Motion:** Pitching
- **Mach:** 0.8
- **AOA:** 1.25°
- **Time Step:** 0.005s

## Fixed Parameters (Issue #1604)
Updated settings based on forum feedback:
- CFL_NUMBER = 1e12
- LINEAR_SOLVER_ERROR = 0.1
- LINEAR_SOLVER_ITER = 10

These changes enable proper unsteady behavior development.

## References
- Website Tutorial: https://su2code.github.io/tutorials/Unsteady_NACA0012/
- Forum Discussion: https://www.cfd-online.com/Forums/su2/242253-tutorial-cases-naca0012-up-date.html
- GitHub Issue: #1604
