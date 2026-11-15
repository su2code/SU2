# Tutorial: setting up and running a FIML test case

Field inversion and machine learning (FIML) corrects the Spalart-Allmaras (SA) production term with a factor
beta, `P -> beta * P`, so that the RANS solution matches reference data (a lift coefficient, a surface
pressure distribution). SU2 finds beta with the discrete adjoint and, optionally, teaches a neural network to
predict it from local flow features, so that the correction carries over to new flows.

This tutorial sets up a case from scratch and runs both methods that SU2 provides:

- **FIML-Classic**: beta is a free field, one value per mesh point.
- **FIML-Direct**: beta is the output of a small neural network, and the network weights are the unknowns.

Worked examples: `FIML_Classic_NACA0012` (small, minutes) and `S809_FIML_Inversion` (production size, hours).
`S809_FIML_Forward` holds the S809 mesh and the baseline configurations.

## 1. The two methods

| | FIML-Classic | FIML-Direct |
|---|---|---|
| Unknowns | beta at every point (about 10^5) | network weights (961 for 3 x 20 neurons) |
| `KIND_TRAIN_NN` | `NONE` | `WEIGHTS` |
| Where beta lives | `Beta_FIML` column of the flow restart | computed from the flow features every iteration |
| Gradient file | `Sensitivity_Beta_FIML` in `restart_adj_custom.csv` | `nn_weights_grad.dat` |
| Result | a beta field for this flow only | a network that predicts beta for any flow |
| Generalizes to other flows | no | yes, if trained on several flows |
| Cost of the gradient | one adjoint | one adjoint (the network is part of the adjoint) |
| Extra input | none | `nn_weights.dat` with weights and feature scaling |
| Typical use | find out what correction the data asks for | produce a usable model |

Field inversion alone fits the data but cannot predict. The usual route is therefore Classic first, then
train a network on the inverted beta, then (optionally) refine the weights with Direct mode:

    baseline SA -> Classic inversion -> train network on beta -> forward run with network -> Direct refinement

A third mode, `KIND_TRAIN_NN= BACKPROP` with `TRAIN_NN= YES`, trains the network inside the flow solver
against a given beta field (`BETA_TARGET_FILE`). It writes `nn_weights.dat` and needs no adjoint, so it is a
shortcut for the training step. The S809 case does that step offline in NumPy (`step4_train.py`) instead.

Both methods share the same SA modification and the same objective. They differ in what the adjoint
differentiates with respect to.

## 2. Requirements

- SU2 built normally (`SU2_CFD`) and with algorithmic differentiation (`SU2_CFD_AD`):

      ./meson.py setup build
      ./meson.py setup build_ad -Denable-autodiff=true
      ninja -C build -j4; ninja -C build_ad -j4

- Python 3 with NumPy and SciPy for the driver `fiml_inversion.py` (L-BFGS-B).
- MPI. Examples below use `mpirun -n 8`.
- Tested with `SOLVER= INC_RANS`. The SA-FIML source terms are also compiled for the compressible solver,
  but the Cp objective output (`INVERSE_DESIGN_PRESSURE`) is only hooked into the incompressible output.

## 3. Common setup

### 3.1 Configuration

Start from a converged SA configuration and add:

    KIND_TURB_MODEL= SA_FIML
    FIML_SHIELD= LEGACY              % beta = 1 where Delta < 0.25 or S/Omega > 2 (boundary layer and wake only)

Output and restart settings that the driver relies on:

    READ_BINARY_RESTART= NO
    OUTPUT_FILES= RESTART_ASCII, SURFACE_CSV
    WRT_RESTART_COMPACT= NO          % otherwise the surface CSV lacks Pressure_Coefficient
    HISTORY_OUTPUT= ITER, RMS_RES, AERO_COEFF, CP_DIFF, FIML, COMBO
    CONV_FILENAME= history
    RESTART_FILENAME= restart_flow
    SOLUTION_FILENAME= restart_flow

With `SA_FIML`, the flow restart carries the extra column `Beta_FIML`, and `BETA_REG` appears in the history
(group `FIML`).

### 3.2 Baseline

Run SA from free stream with `RESTART_SOL= NO`. With `KIND_TRAIN_NN= NONE` and beta = 1 everywhere,
SA-FIML reproduces SA exactly, so this is the baseline to compare against. Check that it converges: all later
steps start from its restart.

    mpirun -n 8 SU2_CFD baseline.cfg

### 3.3 Objective

The objective is a `CUSTOM_OBJFUNC` expression over history fields. Fields useful for FIML:

| Field | Meaning |
|---|---|
| `LIFT`, `DRAG`, `MOMENT_Z`, ... | integrated coefficients |
| `INVERSE_DESIGN_PRESSURE` | sum over wall nodes with a target of vertex length x (Cp - Cp_target)^2 |
| `BETA_REG` | sum over points of (beta - 1)^2, the regularization |

    OBJECTIVE_FUNCTION= CUSTOM_OBJFUNC
    CUSTOM_OBJFUNC= '2500 * pow(LIFT - 1.0546, 2) + 667 * INVERSE_DESIGN_PRESSURE + 0.05 * BETA_REG'

Choosing the weights: use 1/sigma^2 of each data term for the accuracy you trust in it, and set the weight on
`BETA_REG` so that the correction stays moderate. The reference data cannot fix beta uniquely, so this
weight decides how far the field moves. The expression is evaluated on the AD tape, so any combination of
history fields works.

**Cp target file.** `INVERSE_DESIGN_PRESSURE` reads `TargetCp.dat` in the working directory. It is a CSV
file with the header `"PointID","x","y","Pressure","Pressure_Coefficient"`, one row per wall node that has a
target (the global point ID, the coordinates, any pressure value, and the target Cp). Wall nodes not listed do
not contribute, so a measurement that covers only the upper surface can be used as is. Without the file the
term is zero. `S809_FIML_Inversion/step1_target.py` writes this file by interpolating measured data onto the
wall nodes of the mesh.

## 4. FIML-Classic: field inversion

### 4.1 Settings

    KIND_TRAIN_NN= NONE
    TRAIN_NN= NO

Beta is read from the `Beta_FIML` column of `restart_flow.csv` and is initialized to 1 by the baseline run.
The driver edits this column for every evaluation.

The adjoint does not adapt the CFL, so it runs at the fixed CFL that you give with `--adjoint_cfl`.

### 4.2 Running

    python3 ../fiml_inversion.py -n 8 --mode field \
        --su2_cfd  /path/to/build/SU2_CFD/src/SU2_CFD \
        --su2_cfd_ad /path/to/build_ad/SU2_CFD/src/SU2_CFD_AD \
        --maxiter 10 --adjoint_iter 3000 --adjoint_cfl 25 \
        --fields CL CD Beta_Reg

Each evaluation (L-BFGS-B asks for J and the gradient at a new beta):

1. writes beta to `restart_flow.csv` and runs the direct solver from that restart;
2. runs the discrete adjoint on the converged flow and reads `Sensitivity_Beta_FIML`;
3. logs J and the requested history fields to `inversion_history.csv` and keeps beta as `x_eval_NNN.npy`.

At the end the restart holds the flow of the best beta, and `beta_inverted.dat` holds beta in global point
order (one value per line). The direct solves start from the previous restart, so after the first one they
are short.

Useful options:

| Option | Use |
|---|---|
| `--maxiter` | L-BFGS-B iterations; a run usually needs a few more evaluations than iterations (line search) |
| `--adjoint_iter`, `--adjoint_cfl` | adjoint iteration limit and fixed CFL |
| `--warm_adjoint` | start each adjoint from the previous one; faster, but a diverged adjoint spoils the next ones |
| `--scale` | scales the objective seen by the optimizer, to get a first step of a useful size |
| `--beta_bounds` | bounds on beta, default 0 to 2 |
| `--restart`, `--log`, `--output` | file names |

The driver stops if the adjoint diverges (non-finite gradient or no residual reduction). Lower
`--adjoint_cfl` and check that the direct flow is converged.

### 4.3 What to check

- J drops quickly in the first iterations and then slowly; stop when the drop per iteration is small.
- `Beta_Reg` and the beta range (`beta_min`, `beta_max` in `inversion_history.csv`): beta at the bound of 0
  means the optimizer wants to switch production off there, a sign that the target is not reachable with
  this regularization.
- The surface fit: compare the Cp of the inverted flow (`surface_flow.csv`) with the data.
- The direct flow must be converged at each evaluation. A rising `rms[P]` at the end of a direct solve means
  that the gradient will be inaccurate.

NACA0012 example (`FIML_Classic_NACA0012`, CL target 0.95, 4 ranks): J falls from 1.49e-2 to 5.6e-5 in
three iterations (CL 1.072 to 0.949) in about 12 minutes.

## 5. FIML-Direct: network weights as unknowns

### 5.1 The network

- Inputs: the SA production/destruction ratio, chi, Delta = mu_t S / (1.5 tau_w) with a flat-plate wall shear
  stress at x = `REF_LENGTH`/2, and S/Omega.
- Scaling: Box-Cox with one lambda per input (fitted by maximum likelihood), then z-score.
- Hidden layers: `N_HIDDEN_LAYERS` (default 3) of `N_NEURONS` (default 20) tanh neurons. Output
  beta = 1 + tanh(output).
- Points removed by `FIML_SHIELD` keep beta = 1 and are excluded from scaling and training.

### 5.2 The weight file

`nn_weights.dat` must exist for the adjoint, because the feature scaling is stored in it:

    4 Inputs,3 Hidden Layers,20 Hidden Nodes
    <4 means> <4 standard deviations> <4 Box-Cox lambdas>      (one line)
    <weights, one per line, layer by layer, the bias of each layer last>

Three ways to produce it:

1. Offline, from an inverted beta (`S809_FIML_Inversion/step4_train.py`: runs one iteration to dump the
   features, fits Box-Cox, trains in NumPy, writes the file).
2. In SU2 with `KIND_TRAIN_NN= BACKPROP`, `TRAIN_NN= YES`, `BETA_TARGET_FILE= beta_inverted.dat`,
   `LEARNING_RATE`, `NUM_EPOCH`: one Adam step per solver iteration, the file is written after the last epoch.
3. From scratch: without a file, hidden weights start random, the output weights zero (beta = 1), and the
   scaling is computed from the flow at `ITER_START_NN`. A forward run starting this way does not write the
   file, so use 1 or 2 to get a usable one.

### 5.3 Settings

    KIND_TRAIN_NN= WEIGHTS
    TRAIN_NN= NO
    ITER_START_NN= 1                 % the network acts from iteration 1 (the features need a previous iteration)
    N_NEURONS= 20
    N_HIDDEN_LAYERS= 3
    FIML_SHIELD= LEGACY

The flow restart for the first evaluation should come from a forward run with the same network, so that the
direct solves are short (`S809_FIML_Inversion/step5_predict.py`).

### 5.4 Running

    python3 ../fiml_inversion.py -n 8 --mode weights \
        --su2_cfd ... --su2_cfd_ad ... \
        --weights nn_weights.dat --x_scale 0.01 \
        --maxiter 5 --adjoint_iter 5000 --adjoint_cfl 100 \
        --fields CL CD Cp_Diff Beta_Reg --log direct_history.csv

- `--x_scale` sets the length of the first L-BFGS-B step in weight space. The default of 1 is far too large for
  a network (beta jumps to its bounds and J rises by a factor of 20), so choose 0.01 to start and increase if
  the first steps are tiny.
- The adjoint writes `nn_weights_grad.dat` (summed over ranks); the driver reads it, writes the new weights to
  `nn_weights.dat` and runs the next direct solve.
- The direct solves are not warm-started from a separate restart: they continue from `restart_flow.csv`, which
  the previous evaluation wrote. Trial points far from the best can need many more iterations.

### 5.5 Notes on the gradient

- The network and its features are recorded on the AD tape, so the gradient includes how the features respond
  to the flow. The production/destruction input is treated as constant, because the source terms are lagged.
- `FIML_SHIELD= LEGACY` switches points to beta = 1 when a feature crosses its threshold, which makes J slightly
  non-smooth in the weights. For a finite-difference check, set `FIML_SHIELD= NONE`.

## 6. The full pipeline with the S809 case

`S809_FIML_Inversion` runs every step as a separate script. Each reads the files of the previous one and prints
what to check; every script takes `-n`, `--su2_cfd` and, where needed, `--su2_cfd_ad`.

| Step | Script | What it does |
|---|---|---|
| 1 | `step1_target.py` | writes `TargetCp.dat` from measured Cp |
| 2 | `step2_baseline.py` | SA from free stream |
| 3 | `step3_invert.py` | Classic inversion |
| 4 | `step4_train.py` | extracts features, trains the network on the inverted beta |
| 5 | `step5_predict.py` | forward run with the network |
| 6 | `step6_direct.py` | Direct refinement of the weights |
| - | `compare.py` | table of CL, CD, Cp misfit and J per step, and the Cp plot |

Result on the S809 at 14.24 degrees, objective J as above:

| Step | CL | Cp misfit | J |
|---|---|---|---|
| baseline SA | 1.496 | 0.412 | 761.4 |
| Classic | 1.039 | 0.064 | 80.9 |
| network (trained on Classic beta) | 1.087 | 0.083 | 92.3 |
| network, weights refined | 1.027 | 0.055 | 81.6 |

The inversion needs several hours on 8 ranks: an adjoint of 5,000 iterations takes about 28 minutes, 11
minutes with `--warm_adjoint --adjoint_iter 2000`, and a direct solve 2.5 to 12 minutes.

## 7. Setting up your own case

1. Make a converged SA baseline and check that SA-FIML with beta = 1 gives the same coefficients.
2. Prepare the data: the target lift, and a `TargetCp.dat` if you have pressure data.
3. Write `CUSTOM_OBJFUNC` and evaluate it on the baseline with a 1-iteration run from the baseline restart;
   all terms should be finite and J should be the number you expect.
4. Choose the adjoint CFL: run the adjoint once on the baseline (`MATH_PROBLEM= DISCRETE_ADJOINT`,
   `RESTART_SOL= NO`, fixed `CFL_NUMBER`, `CFL_ADAPT= NO`; the flow restart is always read) and make sure the
   residual falls. On a small mesh CFL 25 worked, on the S809 CFL 100. The history field `SENS_FIML`, the norm
   of dJ/dbeta or dJ/dw, shows when the gradient has settled.
5. Check the gradient once, before an expensive run: central differences of J along the gradient direction
   of beta or of the weights (`FIML_SHIELD= NONE` or `FD`; `LEGACY` makes J jump) against the adjoint gradient.
   Converge the perturbed flows until J stops changing, not only until the residual is small: with
   `1e3 * pow(LIFT - 0.95, 2)` in J, a residual of -12 still leaves a drift of about 1e-4 in J, which makes the
   difference quotient wrong by 50% for a step of 1e-4. Converged to -14, a step of 5e-3 in the weights of a
   network with bounded inputs and `FIML_SHIELD= FD` agrees with the adjoint to 0.1% (NACA0012).
6. Run Classic (section 4). Look at beta and the fit.
7. If you want a model, train (section 5.2), run forward, and refine with Direct mode.

## 8. Troubleshooting

| Symptom | Likely cause |
|---|---|
| "The adjoint of the neural network needs nn_weights.dat." | Direct mode adjoint without a weight file; create it (section 5.2) |
| Adjoint residual rises, or the driver stops with "the adjoint diverged" | CFL of the adjoint too high, or the direct flow not converged |
| Surface CSV has no `Pressure_Coefficient` | `WRT_RESTART_COMPACT= YES`; set it to `NO` |
| `INVERSE_DESIGN_PRESSURE` is zero | `TargetCp.dat` missing from the working directory, or not listing wall nodes |
| J does not fall in the first Direct iterations | `--x_scale` too small, or beta of most points frozen by `FIML_SHIELD` |
| J jumps up by a large factor at the first Direct evaluation | `--x_scale` too large |
| Network gives a worse J than the inverted field | beta outside the filter cannot be reproduced; refine with Direct mode |

## 9. Reference: FIML options

| Option | Default | Meaning |
|---|---|---|
| `KIND_TURB_MODEL= SA_FIML` | - | SA with the production factor beta |
| `KIND_TRAIN_NN` | `WEIGHTS` | `NONE` (Classic), `BACKPROP` (train inside the solver), `WEIGHTS` (Direct) |
| `TRAIN_NN` | `NO` | train the network inside the solver (`BACKPROP`) |
| `N_NEURONS`, `N_HIDDEN_LAYERS` | 20, 3 | network size |
| `ITER_START_NN` | 1000 | iteration from which the network acts (at least 1 is used) |
| `NUM_EPOCH`, `LEARNING_RATE` | 1, 0.01 | in-solver training |
| `BETA_TARGET_FILE` | `beta_target.dat` | one beta per line, global point order, for in-solver training |
| `FIML_INPUTS` | `(PROD_DEST, CHI, DELTA, STRAIN_VORT)` | network inputs, in order (section 9.1) |
| `FIML_SHIELD` | `NONE` | where the correction is switched off (section 9.2) |
| `FIML_SHIELD_COEFF` | 8 | C in f_d = 1 - tanh((C r_d)^3) for `FIML_SHIELD= FD` |
| `SA_OPTIONS` | - | the usual SA options also apply to `SA_FIML` |

History fields: `BETA_REG` (group `FIML`), `INVERSE_DESIGN_PRESSURE` (group `CP_DIFF`, column `Cp_Diff`), and in
the adjoint `SENS_FIML` (group `SENSITIVITY`), the L2 norm of dJ/dbeta (`NONE`) or dJ/dw (`WEIGHTS`).
Volume fields: `Beta_FIML` in the group `FIML`, and all inputs of section 9.1 in the group `FIML_INPUTS`
(columns `FIML_Prod_Dest`, ..., `FIML_Fd`). The group `FIML` is always written, so beta is part of every flow restart.

### 9.1 Network inputs

`FIML_INPUTS` lists the inputs in the order the network uses them. All are local, dimensionless and Galilean
invariant, except `DELTA`, which divides by a wall stress built from the free stream.

| Name | Definition | Range | Scaling |
|---|---|---|---|
| `PROD_DEST` | SA production / destruction | >= 0 | Box-Cox + z-score |
| `CHI` | nu_tilde / nu | >= 0 | Box-Cox + z-score |
| `DELTA` | mu_t S / (1.5 tau_w), tau_w of a flat plate at half `REF_LENGTH` in the free stream | >= 0 | Box-Cox + z-score |
| `STRAIN_VORT` | S / Omega | >= 0 | Box-Cox + z-score |
| `VISC_RATIO` | nu_t / (nu_t + 100 nu) | [0, 1) | none |
| `PROD_DEST_RATIO` | P / (P + D) of SA | [0, 1] | none |
| `Q_CRIT` | (Omega^2 - S^2) / (Omega^2 + S^2) | [-1, 1] | none |
| `BL_RATIO` | 1 / (1 + r_d), r_d = (nu_t + nu) / (kappa^2 d^2 sqrt(u_i,j u_i,j)), r_d <= 10 | [1/11, 1] | none |
| `FD` | 1 - tanh((8 r_d)^3) | [0, 1] | none |

- The Box-Cox constants are fitted on the training case and stored in `nn_weights.dat`; the bounded inputs need
  no constants, so a new case cannot push them outside the range the network was trained on.
- `PROD_DEST` and `PROD_DEST_RATIO` use the SA source terms of the previous iteration; in the adjoint they are
  constants. The other inputs are differentiated.
- The network file names its inputs on the line after the header, `Inputs: VISC_RATIO Q_CRIT BL_RATIO`, and SU2
  stops if they differ from `FIML_INPUTS`. A file without that line holds `PROD_DEST CHI DELTA STRAIN_VORT`.
- `step4_train.py --inputs VISC_RATIO,Q_CRIT,BL_RATIO --shield FD` trains on the columns of `FIML_INPUTS` and
  writes this file.

### 9.2 Shield

| `FIML_SHIELD` | Effect |
|---|---|
| `NONE` | the correction acts everywhere |
| `LEGACY` | the network sets beta = 1 where `DELTA` < 0.25 or `STRAIN_VORT` > 2; no effect in Classic mode |
| `FD` | the source term uses 1 + (beta - 1) f_d, f_d = 1 - tanh((C r_d)^3), in Classic and network modes |

f_d is the shielding function of delayed DES: 0 in attached boundary layers, 1 outside them. With `FD` the
correction cannot change attached boundary layers: on the turbulent flat plate (`rans/flatplate`) beta = 2
everywhere changes the skin friction by 13-17% without a shield, by 0.7% with C = 8, by 3-4% with C = 4 and by
7-10% with C = 2. A smaller C leaves the outer part of the boundary layer to the correction, which matters when
the error lies in an adverse-pressure-gradient boundary layer, as at trailing-edge stall.

Regression tests:

| Test | Script | What it checks |
|---|---|---|
| `inc_turb_s809_sa`, `inc_turb_s809_fiml` | `parallel_regression.py` | S809 against the measured Cp: SA, and SA with the trained network |
| `inc_turb_s809_fiml` | `hybrid_regression.py` | the network with OpenMP threads |
| `discadj_incomp_turb_NACA0012_fiml_beta` | `parallel_regression_AD.py` | adjoint with beta as the design variable (Classic) |
| `discadj_incomp_turb_NACA0012_fiml_weights` | `parallel_regression_AD.py` | adjoint with the network weights as design variables (Direct) |

## 10. Background

- Parish & Duraisamy, "A paradigm for data-driven predictive modeling using field inversion and machine
  learning", J. Comput. Phys. 305 (2016).
- Singh, Medida & Duraisamy, "Machine-learning-augmented predictive modeling of turbulent separated flows
  over airfoils", AIAA J. 55 (2017).
- J. Holland, Field inversion and machine learning in SU2 v5 (github.com/jholland1/SU2, branch `feature_fiml`),
  which has FIML-Classic, -Embedded (`BACKPROP`) and -Direct. The objectives and options here follow it,
  except that the objective is a `CUSTOM_OBJFUNC` and the regularization is an ordinary history field.
