# S809 FIML: field inversion, training and prediction

Field inversion and machine learning for the S809 airfoil at AoA = 14.24 deg, Re = 2e6, against the
measured upper-surface pressure and lift. Mesh, flow conditions and numerics are those of
`../S809_FIML_Forward` (SA, 113,747 points).

## Data

- `s809_cp_upper_exp.dat`: upper-surface Cp from the Delft low-turbulence tunnel (Somers,
  NREL/SR-440-6918, 1997, fig. 8), as digitized by J. Holland (`Cp_Grabit_Upper.mat` in
  github.com/jholland1/SU2, branch `feature_fiml`). The file does not state the angle of attack; the
  separated plateau from midchord and the suction peak of -7 fit the 14.24 deg case. With the SA
  lower surface it integrates to CL ~ 1.07.
- Lift: CL = 1.0546, the target Holland used. Wolfe & Ochs (AIAA-97-0973) list 1.1104 for the same Delft
  data at 14.24 deg.

## Objective

    CUSTOM_OBJFUNC= '2500 * pow(LIFT - 1.0546, 2) + 667 * INVERSE_DESIGN_PRESSURE + 0.05 * BETA_REG'

- `INVERSE_DESIGN_PRESSURE` = sum over the wall nodes in `TargetCp.dat` of (Cp - Cp_target)^2 times the
  vertex length (m). Wall nodes without a target do not contribute.
- `BETA_REG` = sum over all points of (beta - 1)^2.
- The weights are 1/sigma^2 with sigma_CL = 0.02 and sigma_Cp = 0.05 (per unit of the 0.6 m upper surface).

## Steps

Every step reads the files of the previous one and prints the quantities to check. All steps take
`-n` (MPI ranks), `--su2_cfd` and, where needed, `--su2_cfd_ad`.

| Step | Script | Writes | Check |
|---|---|---|---|
| 1 | `step1_target.py` | `TargetCp.dat` | number of wall nodes with a target, x/c range |
| 2 | `step2_baseline.py` | `restart_baseline.csv`, `surface_baseline.csv`, `history_baseline.csv` | CL, CD and J of SA; converged by the Cauchy test |
| 3 | `step3_invert.py` | `restart_inverted.csv`, `inversion_history.csv` | J and its terms per evaluation; beta range |
| 4 | `step4_train.py` | `features.csv`, `nn_weights.dat` | training and validation error of the network |
| 5 | `step5_predict.py` | `restart_predicted.csv` | CL and J with the network; predicted vs inverted beta |
| 6 | `step6_direct.py` | `nn_weights_direct.dat`, `direct_history.csv` | J when the weights are optimized with the adjoint |
| - | `compare.py` | `cp_upper.png` | table of all steps, Cp plot |

- Step 3 runs `../fiml_inversion.py --mode field` (L-BFGS-B on beta at every point, `KIND_TRAIN_NN= NONE`).
- Step 4 trains offline, in NumPy, on the points that `FILTER_SHIELD` keeps. It uses the scaling, network and
  file layout of `CTurbSASolver`, so SU2 reproduces the trained beta.
- Step 6 runs `../fiml_inversion.py --mode weights`: `KIND_TRAIN_NN= WEIGHTS`, the adjoint gives dJ/dw
  through the network and its features (`nn_weights_grad.dat`).

## Results

`compare.py` prints the table and writes `cp_upper.png` (upper-surface Cp of every step against the experiment).

| Step | CL | CD | Cp_Diff | Beta_Reg | J | CL term | Cp term | beta term |
|---|---|---|---|---|---|---|---|---|
| 2 baseline (SA, beta = 1) | 1.4957 | 0.0403 | 0.4121 | 0 | 761.4 | 486 | 275 | 0 |
| 3 field inversion | 1.0385 | 0.0998 | 0.0640 | 751 | 80.9 | 0.65 | 42.7 | 37.6 |
| 5 network (trained on step 3) | 1.0872 | 0.0904 | 0.0834 | 681 | 92.3 | 2.7 | 55.7 | 34.0 |
| 6 network, weights optimized | 1.0269 | 0.1066 | 0.0549 | 862 | 81.6 | 1.9 | 36.6 | 43.1 |

Target CL = 1.0546. Upper-surface Cp at x/c = 0.1, 0.3, 0.5, 0.7, 0.9:

| | 0.1 | 0.3 | 0.5 | 0.7 | 0.9 |
|---|---|---|---|---|---|
| baseline | -2.49 | -1.55 | -0.97 | -0.25 | -0.14 |
| field inversion | -1.75 | -0.83 | -0.50 | -0.51 | -0.47 |
| network | -1.85 | -0.96 | -0.47 | -0.48 | -0.45 |
| network, weights optimized | -1.70 | -0.73 | -0.54 | -0.53 | -0.48 |
| experiment | -1.2 | -0.58 | -0.61 | -0.63 | -0.64 |

- Step 3: J falls from 761 to 81 in the first 7 evaluations, then little progress; beta in [0, 1.23], differing from 1 by
  more than 0.1 at 2,982 of 113,747 points. The reported runs restarted each adjoint from the previous one for 2,000
  iterations (about 11 min on 8 ranks, `--warm_adjoint --adjoint_iter 2000`); a cold adjoint needs about 5,000
  iterations (28 min). A direct solve takes 2.5 to 12 min.
- Step 4: 17,323 points inside the filter; validation RMS(beta) 0.021, R^2 = 0.989 for the 961 weights.
  The 1,058 points outside the filter with |beta - 1| > 0.05 stay at beta = 1 in SU2.
- Step 5: stops at ITER with CL still falling by 0.001 per 1000 iterations.
- Step 6: 5 L-BFGS-B iterations (10 evaluations), first step of length 0.01 in weight space (`--x_scale`),
  started with `--resume` from weights that an earlier run had brought to J = 86.4.
- The remaining misfit is ahead of x/c = 0.4, where the plateau of the experiment starts at x/c ~ 0.15 and
  the computed suction stays higher. The beta penalty keeps the correction from going further.
- The adjoint runs at a fixed CFL (`--adjoint_cfl`, default 100) and the driver stops if its residual does not
  fall. `--warm_adjoint` starts each adjoint from the previous one: faster, but a diverged adjoint then
  spoils the following ones.
