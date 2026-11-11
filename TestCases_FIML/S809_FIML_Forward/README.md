# S809 Airfoil FIML Test Case

S809 airfoil with a Field Inversion and Machine Learning (FIML) correction of the
Spalart-Allmaras production term.

## Case

- **Mesh**: `s809_cmesh.su2`, structured C-mesh with 113,747 points and 112,896 quadrilaterals,
  chord c = 0.6 m, far field and wake at 100 chords, first cell height 1.8e-6 m (y+ ~ 0.3).
- **Flow**: incompressible RANS, rho = 1.2886 kg/m^3, mu = 1.853e-5 Pa s,
  |V| = 47.9332 m/s at AoA = 14.24 deg, Re = 2.0e6 based on chord.
- **Target**: experimental CL = 1.0546 at this angle.

## Mesh generation

`s809_cmesh.py` builds the mesh from `s809.dat` (NREL S809 coordinates from airfoiltools.com,
`s809-nr`) with the gmsh Python API, with the markers `airfoil` and `farfield`:

    python3 s809_cmesh.py --output s809_cmesh.su2              # the mesh used by the configurations
    python3 s809_cmesh.py --output s809_fine.su2 --scale 2     # all spacings halved, for a grid study

Without `--output` it writes `s809_cmesh.msh`.

- Grid lines leave the wall along the normal, then turn to their far-field point; the turn of the
  start directions at the sharp trailing edge is smoothed so that grid lines do not converge.
- The far-field end of each grid line depends on its position along the C-line, not on the wall or
  wake spacings, so refining one spacing leaves the other grid lines in place.
- Defaults: first cell 3e-6 c (y+ ~ 0.3 at Re 2e6), growth 1.1 up to 0.04 c, 129 points normal
  to the wall, 564 points on the airfoil with spacing 1e-4 c at the leading edge and 6e-4 c at the
  trailing edge, wake and far field at 100 c.
- Quality: no cell has a corner angle below 45 deg; cells within 0.05 c of the wall stay above 76 deg.

## Configuration files

| File | Purpose |
|---|---|
| `config_sa_baseline.cfg` | Baseline SA from free stream; writes `restart_sa.dat`. |
| `config_fiml_forward.cfg` | SA-FIML from `restart_sa.dat` with the weights in `nn_weights.dat`. |
| `config_fiml_training.cfg` | SA-FIML from `restart_sa.dat`, trains the network on `BETA_TARGET_FILE`. |

All three share the same numerics:

- `NUM_METHOD_GRAD_RECON= LEAST_SQUARES` with the Venkatakrishnan limiter, `MUSCL_TURB= YES`;
- CFL ramped from 1 and capped at 100;
- the Cauchy test on `LIFT` (`CONV_CAUCHY_EPS= 1E-7` over 500 iterations) stops the run, `ITER= 15000` is an upper limit.

## Workflow

1. Run the baseline: `mpirun -n 8 SU2_CFD config_sa_baseline.cfg`.
2. Provide `beta_target.dat`: one beta value per line, ordered by global point index
   (113,747 lines). Points without a value get beta = 1.
3. Train: `mpirun -n 8 SU2_CFD config_fiml_training.cfg`. One Adam step runs per solver
   iteration after `ITER_START_NN`; `nn_weights.dat` is written after the last epoch.
4. Run forward: `mpirun -n 8 SU2_CFD config_fiml_forward.cfg`.

Without `nn_weights.dat`, hidden layers start from random weights and the output layer from
zero, so beta = 1 and SA-FIML reproduces SA exactly.

SA-FIML writes the field `Beta_FIML` (output group `FIML`) to the volume output and to the
restart file. On restart it is read back, so beta is continuous before `ITER_START_NN`.

`nn_weights.dat` holds a header line, then a line with the 4 input means, 4 stds and
4 Box-Cox lambdas, then all weights layer by layer (bias last in each layer).

## Network

- Inputs (Holland et al.): the SA production/destruction ratio taken from the source term (before beta),
  chi, Delta = mu_t S / (1.5 tau_w) with tau_w = 0.027 Re_x^(-1/7) rho U^2 at x = REF_LENGTH/2, and S/Omega.
- Scaling: Box-Cox with lambda fitted per input by maximum likelihood, then z-score.
- `FILTER_SHIELD= YES` keeps beta = 1 where Delta < 0.25 or S/Omega > 2, and leaves those points out of
  the scaling and the training, so that only the boundary layer and the wake are trained.
- The network starts at `max(ITER_START_NN, 1)`, since the features use the previous source evaluation.
- 3 hidden tanh layers of 20 neurons; output beta = 1 + tanh(.).

## Baseline results (SA, second order)

| AoA | CL | CD | Note |
|---|---|---|---|
| 4.10 deg | 0.6016 | 0.0135 | stopped by the Cauchy test after 2,407 iterations |
| 14.24 deg | 1.4957 | 0.0403 | stopped by the Cauchy test after 10,422 iterations; separation at x/c = 0.62 |

SA over-predicts lift at 14.24 deg by about 42% (experiment 1.0546), so the correction has
to reduce production (beta < 1) where the flow should separate.

### Grid dependence at 14.24 deg

The separation point, and with it CL, depends mainly on the wall spacing at the leading edge, where
the suction peak is. With `--ds_te 8e-4 --ds_upper 5e-3 --blend 2` and the configuration above:

| `--ds_le` | CL | separation x/c |
|---|---|---|
| 8e-4 | 1.4006 | 0.565 |
| 4e-4 | 1.4635 | 0.600 |
| 2e-4 | 1.4867 | 0.614 |
| 1e-4 | 1.4949 | 0.624 |

Starting from the first row, refining the wall-normal direction by a factor 1.4 changes CL by -0.3%.
Refining all spacings along the wall and the wake by factors 1.4 and 2 gives CL 1.4567 and 1.4853,
which extrapolates to a grid-converged CL of about 1.515, 1.3% above the default mesh.

## Field inversion

`../S809_FIML_Inversion` inverts beta against the measured upper-surface Cp and lift on this mesh, trains the
network on the inverted field, and optimizes the network weights with the adjoint (`KIND_TRAIN_NN= WEIGHTS`).

## References

- Parish & Duraisamy, "A paradigm for data-driven predictive modeling using field inversion
  and machine learning", J. Comput. Phys. 305 (2016).
- Singh, Medida & Duraisamy, "Machine-learning-augmented predictive modeling of turbulent
  separated flows over airfoils", AIAA J. 55 (2017).
