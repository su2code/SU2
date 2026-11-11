# FIML field inversion (classic), NACA0012

Field inversion of the SA production correction beta at every grid point, with gradients from the
discrete adjoint and L-BFGS-B from SciPy.

## Case

- **Mesh**: `n0012_225-65.su2` from the SU2 TestCases repository
  (`disc_adj_incomp_rans/naca0012`), 14,576 points.
- **Flow**: incompressible SA, same conditions as the `disc_adj_incomp_rans` NACA0012 test
  (AoA 10 deg); baseline CL = 1.0722, CD = 0.01738.
- **Target**: CL = 0.95, a synthetic value that requires beta < 1.

## Problem

    OBJECTIVE_FUNCTION= CUSTOM_OBJFUNC
    CUSTOM_OBJFUNC= 'pow(LIFT - 0.95, 2) + 1e-6 * BETA_REG'

- `KIND_TURB_MODEL= SA_FIML` with `KIND_TRAIN_NN= NONE`: beta is a field, no neural network.
- beta is read from and written to the column `Beta_FIML` of the flow restart.
- `BETA_REG` is the history output sum_i (beta_i - 1)^2 (group `FIML`). Any history output can enter the
  expression, for instance `DRAG` or `INVERSE_DESIGN_PRESSURE` (Cp difference to `TargetCp.dat`, group `CP_DIFF`).
- The discrete adjoint writes dJ/dbeta as `Sensitivity_Beta_FIML` to the adjoint volume output and to
  `restart_adj_custom.csv`.

## Running

1. Baseline from free stream, which writes `restart_flow.csv` with beta = 1:

       sed 's/^RESTART_SOL=.*/RESTART_SOL= NO/' config_inversion.cfg > baseline.cfg
       mpirun -n 8 SU2_CFD baseline.cfg

2. Inversion (needs `SU2_CFD` and `SU2_CFD_AD`, built with `-Denable-autodiff=true`):

       python3 ../fiml_inversion.py -n 8 --maxiter 10 --adjoint_iter 3000

   Each evaluation runs the direct solver from the previous restart, then the adjoint. The script
   logs J and the history fields given with `--fields` per evaluation to `inversion_history.csv`, keeps beta
   of every evaluation as `x_eval_NNN.npy`, and writes the final field to `beta_inverted.dat` (ordered by
   global point index).

## Gradient verification

Central differences along three directions d, at a non-uniform beta and a weight of 1e-4 on `BETA_REG`,
against g.d from the adjoint (direct runs converged to rms[P] = -12):

| Direction | FD, eps = 1e-2 | Adjoint | Ratio |
|---|---|---|---|
| normalized gradient | 6.842e-3 | 6.890e-3 | 0.993 |
| patch 0.3 < x < 0.9, \|y\| < 0.05 | -5.159e-3 | -5.104e-3 | 1.011 |
| random, within one chord | -7.716e-4 | -7.108e-4 | 1.086 |

The adjoint gradient changes by 0.24% between 4,000 and 16,000 adjoint iterations; the
remaining differences are finite-difference error (they grow at eps = 1e-3).

## Result

`../fiml_inversion.py -n 8 --maxiter 10 --adjoint_iter 3000`, starting from beta = 1, weight 1e-6 on `BETA_REG`:

| Evaluation | J | CL | CD |
|---|---|---|---|
| 0 | 1.492e-2 | 1.0721 | 0.01738 |
| 4 | 8.96e-5 | 0.9444 | 0.02033 |
| 9 | 5.94e-5 | 0.9500 | 0.01991 |
| 13 | 5.73e-5 | 0.9505 | 0.01994 |

CL reaches the target within 9 evaluations (14 in total, about 2.4 min each on 8 ranks); J is then
the regularization term. beta differs from 1 by more than 0.1 at 368 points, all in the boundary
layer (|y| < 0.07); at 18 points beta is below 0.1. A larger weight on `BETA_REG` gives a smoother, less
intrusive correction.
