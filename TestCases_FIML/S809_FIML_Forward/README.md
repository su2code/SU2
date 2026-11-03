# S809 Airfoil FIML Test Case

S809 airfoil with a Field Inversion and Machine Learning (FIML) correction of the
Spalart-Allmaras production term.

## Case

- **Mesh**: `S809_struct_coarse.su2`, 63,400 points, 62,766 quadrilaterals, chord c = 0.6 m,
  far field at about 83 chords, first cell height 7.5e-7 m (y+ ~ 0.1).
- **Flow**: incompressible RANS, rho = 1.2886 kg/m^3, mu = 1.853e-5 Pa s,
  |V| = 47.9332 m/s at AoA = 14.24 deg, Re = 2.0e6 based on chord.
- **Target**: experimental CL = 1.0546 at this angle.

## Configuration files

| File | Purpose |
|---|---|
| `config_sa_baseline.cfg` | Baseline SA from free stream; writes `restart_sa.dat`. |
| `config_fiml_forward.cfg` | SA-FIML from `restart_sa.dat` with the weights in `nn_weights.dat`. |
| `config_fiml_training.cfg` | SA-FIML from `restart_sa.dat`, trains the network on `BETA_TARGET_FILE`. |

All three share the same numerics. The cells behind the sharp trailing edge are strongly
skewed (minimum angle 16 deg), which requires:

- `NUM_METHOD_GRAD_RECON= LEAST_SQUARES` with the Venkatakrishnan limiter;
- CFL ramped from 1 and capped at 10 (a cap of 20 diverges at 14.24 deg).

## Workflow

1. Run the baseline: `mpirun -n 8 SU2_CFD config_sa_baseline.cfg`.
2. Provide `beta_target.dat`: one beta value per line, ordered by global point index
   (63,400 lines). Points without a value get beta = 1.
3. Train: `mpirun -n 8 SU2_CFD config_fiml_training.cfg`. One epoch runs per solver
   iteration after `ITER_START_NN`; `nn_weights.dat` is written after the last epoch.
4. Run forward: `mpirun -n 8 SU2_CFD config_fiml_forward.cfg`.

Without `nn_weights.dat`, hidden layers start from random weights and the output layer from
zero, so beta = 1 and SA-FIML reproduces SA exactly.

SA-FIML writes the field `Beta_FIML` (output group `FIML`) to the volume output and to the
restart file. On restart it is read back, so beta is continuous before `ITER_START_NN`.

`nn_weights.dat` holds a header line, then a line with the 4 input means, 4 stds and
4 Box-Cox lambdas, then all weights layer by layer (bias last in each layer).

## Network

- Inputs: production/destruction ratio, chi, normalized Q-criterion, strain/vorticity ratio.
- Scaling: Box-Cox with lambda fitted per input by maximum likelihood, then z-score.
- 3 hidden tanh layers of 20 neurons; output beta = 1 + tanh(.).

## Baseline results (SA, second order)

| AoA | CL | CD | Note |
|---|---|---|---|
| 4.10 deg | 0.599 | 0.0142 | converged to rms[P] = -9.9 in 4,000 iterations |
| 14.24 deg | 1.383 | 0.0524 | after 29,000 iterations, still drifting; extrapolated CL 1.370, CD 0.0536 |

SA over-predicts lift at 14.24 deg by about 30% (experiment 1.0546), so the correction has
to reduce production (beta < 1) where the flow should separate.

## Not yet implemented

- Field inversion: the FIML objectives (`INVERSE_DESIGN_*_FIML`) and the `FIML` design
  variable are defined but not connected to the adjoint.
- `KIND_TRAIN_NN= WEIGHTS`, `SA_FIML_CORRECTION` other than `PRODUCTION`, `FILTER_SHIELD`.

## References

- Parish & Duraisamy, "A paradigm for data-driven predictive modeling using field inversion
  and machine learning", J. Comput. Phys. 305 (2016).
- Singh, Medida & Duraisamy, "Machine-learning-augmented predictive modeling of turbulent
  separated flows over airfoils", AIAA J. 55 (2017).
