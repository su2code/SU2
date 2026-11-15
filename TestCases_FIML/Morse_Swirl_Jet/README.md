# Swirling jet in still air (Morse 1980)

Swirling free jet at swirl number S = 0.36 and Re = 56,000, ERCOFTAC Classic Collection case 026
("Single-stream swirling jet in still air"). Swirl increases the mixing of a jet: it spreads faster and its
centreline velocity and swirl decay faster. Eddy-viscosity models underpredict this, which makes the case a test
for a FIML correction of SA.

## Data

The ERCOFTAC files (not in this repository) hold the swirling case only:

- `morse0.dat`: centreline velocity U_cl(x) up to x/D = 15, maximum swirl velocity W_max and half-width r_1/2 at
  x/D = 0.5, 1, 2, 4, 6, 10.
- `morse1.dat` ... `morse6.dat`: profiles of U/U_cl, V/U_cl, W/W_max and the six Reynolds stresses (divided by
  U_cl^2) against r/x at x/D = 0.5, 1, 2, 4, 6, 10.

The files do not state the nozzle diameter. With U_exit = 34.46 m/s, Re = 56,000 and air, D = 0.0246 m.

## Files

| File | Purpose |
|---|---|
| `morse_wedge.geo` | gmsh script: 5 deg wedge, x/D from 0.5 to 50, r/D up to 30, rotationally periodic sides |
| `morse_inlet.py` | inlet profile at x/D = 0.5 from morse1.dat; nu_t = −u'v'/(dU/dr) from the measured shear stress, or uniform with --nut_ratio| 
| `morse_sa.cfg` | SA, incompressible |
| `morse_data.py` | readers for the measurements and the SU2 restart file, used by the plot scripts |
| `plot_centreline.py` | U_cl, W_max and r_1/2 against x/D, `morse_centreline.png`; prints the errors |
| `plot_profiles.py` | U/U_cl, V/U_cl, W/W_max against r/D at the stations, `morse_profiles.png` |
| `plot_stresses.py` | u'v' and v'w' (Boussinesq, from the eddy viscosity) at the stations, `morse_stresses.png` |
| `morse_target.py` | FIML objective: probes with the measured velocities as targets; writes `config_inversion.cfg`, `config_predict.cfg`, `morse_targets.csv` |

## Running

    gmsh morse_wedge.geo -3 -format su2 -o morse_wedge.su2
    SU2_CFD morse_sa.cfg                      # writes example_inlet_jet.dat and stops
    python3 morse_inlet.py --data <dir>/morse1.dat 
    mpirun -n 8 SU2_CFD morse_sa.cfg
    python3 plot_centreline.py --data_dir <dir> --case SA=restart_flow.csv
    python3 plot_profiles.py --data_dir <dir> --case SA=restart_flow.csv
    python3 plot_stresses.py --data_dir <dir> --case SA=restart_flow.csv

`--case label=file` may be repeated to compare several solutions in one figure.

## Setup

- **Mesh.** 165,100 points, one cell in the azimuthal direction. Uniform radial spacing 0.02 D up to r = 4 D,
  axial spacing 0.02 D at the inlet growing to 0.25 D.
- **Axis.** A slip wall at r = 0.02 D replaces the axis. At r = 0.002 D the cells next to it are 1.7e-4 D wide and
  give an artificial velocity deficit on the axis, also without swirl.
- **Far field and outlet.** `MARKER_FAR` with a co-flow of 1% of U_exit on the lateral boundary, the inlet plane
  outside the jet and the outlet. The incompressible far field fixes the velocity, so all entrained air enters
  through the inlet plane; with r_max = 30 D this suffices up to x/D = 100. At the outlet the jet is forced to the
  co-flow velocity, which decelerates it over the last 15 D. A pressure outlet (`MARKER_OUTLET`) lets a reversed flow
  enter along the axis; openings at the ambient pressure on the far field were unstable.
- **Wall distance.** The domain has no viscous wall. SU2 then sets the wall distance to 1e15 and SA takes its
  free-shear form (production cb1 Omega nu_tilde, no destruction).

## SA baseline

12,000 iterations; rms[P] stalls at -7.7. Of the remaining residual, 90% lies at the lateral far field (r > 20 D)
and less than 5% at x/D < 15. Errors against the measurements:

| Inflow nu_t / nu | U_cl at x/D = 4, 10, 15 | W_max at x/D = 2, 4, 6 [m/s] | r_1/2 at x/D = 2, 4, 6, 10 |
|---|---|---|---|
| 100 | +41%, +20%, +5% | +1.8, +3.1, +1.4 | -14%, -32%, -32%, -23% |
| 300 | +16%, -6%, -13% | -2.6, 0.0, +0.1 | -5%, -18%, -13%, -3% |

- SA mixes too little in the swirling near field: the jet spreads too slowly up to x/D = 6, and with
  nu_t / nu = 100 also the centreline velocity and the swirl decay too slowly.
- The inflow eddy viscosity is not measured. At x/D = 1 the shear stress u'v' matches with nu_t / nu = 100 and is
  60% too large with 300; downstream, 300 is closer. The difference between 100 and 300 is as large as the
  model error, so the inflow has to be fixed before a FIML correction near the inlet can be told apart from an
  inflow error.
- The far-field outlet decelerates the jet to the co-flow velocity over the last 15 D. Moving the outlet from
  x/D = 30 to 50 changes U_cl at x/D = 15 by less than 0.5%.

## FIML

The field inversion finds a factor beta(x) on the SA production that brings the solution closest to the
measurements; a network then learns beta from local flow quantities. The steps are those of
`../S809_FIML_Inversion`, with velocity probes instead of the surface pressure as data.

### Requirements

- `SU2_CFD` and `SU2_CFD_AD` of branch `feature_fiml`. The AD build uses `-Denable-autodiff=true
  -Denable-normal=false`; below, both executables are on the `PATH`.
- A converged SA solution `restart_flow.csv` (ASCII) of `morse_sa.cfg`, with `inlet_jet.dat` and `morse_wedge.su2`
  in the case directory.

### Objective

    J = Misfit_Ucl + Misfit_U + Misfit_W + 0.05 * BETA_REG,    Misfit = sum ((u - u_exp) / sigma)^2

- `Misfit_Ucl`: centreline velocity at the 37 measured x/D from 1 to 15.
- `Misfit_U`, `Misfit_W`: axial and swirl velocity at the 79 measured radii of the stations x/D = 1, 2, 4, 6, 10.
- sigma = 0.02 U_exit = 0.69 m/s. `BETA_REG` = sum over all points of (beta - 1)^2.
- Each target is a `Probe` at a mesh node of the first azimuthal plane. The measurement is interpolated to the radius
  of that node, so the probe needs no interpolation. The swirl velocity is the tangential component in that plane.
- An SU2 expression holds at most 255 nodes, so each misfit is a sum of partial sums of 30 terms (`part_*`).

### Steps

1. **Objective and configurations.**

       python3 morse_target.py --data_dir <dir>

   Writes `config_inversion.cfg` (SA_FIML, `KIND_TRAIN_NN= NONE`, probes, objective), `config_predict.cfg` (the same
   with the network, `KIND_TRAIN_NN= WEIGHTS` and `FIML_INPUTS` from `--inputs`) and `morse_targets.csv` (name,
   position and target of every probe). `--stations`, `--centreline`, `--sigma`, `--beta_reg` and `--no_swirl`
   change the objective.

2. **Baseline, beta = 1.** SA_FIML restarts from the SA solution; the restart has no `Beta_FIML`, so beta starts at 1.

       cp restart_flow.csv restart_sa.csv
       mpirun -n 8 SU2_CFD config_inversion.cfg
       cp restart_flow.csv restart_baseline.csv

   Check: the run stops after a few iterations with J = 2224 (`Misfit_Ucl` 1036, `Misfit_U` 987, `Misfit_W` 201).

3. **Field inversion.** L-BFGS-B on beta at every point; each evaluation is a direct run from the previous solution
   and an adjoint run.

       python3 ../fiml_inversion.py -f config_inversion.cfg -n 8 --mode field \
           --adjoint_iter 300 --adjoint_cfl 100 --maxiter 15 \
           --fields Misfit_Ucl Misfit_U Misfit_W Beta_Reg
       cp restart_flow.csv restart_inverted.csv
       python3 plot_centreline.py --data_dir <dir> --case SA=restart_baseline.csv --case inverted=restart_inverted.csv

   `inversion_history.csv` holds J and its terms per evaluation, `beta_inverted.dat` the final beta. Run it in a copy
   of the case directory: the driver overwrites `restart_flow.csv`.

4. **Training.** One iteration from `restart_inverted.csv` writes the network inputs (`features.csv`); the network is
   trained offline on beta and written to `nn_weights.dat`.

       python3 ../S809_FIML_Inversion/step4_train.py -n 8 --inputs CHI,VISC_RATIO,STRAIN_VORT,Q_CRIT --shield NONE

   `--inputs` must equal `--inputs` of step 1. The other inputs need a wall: without one, PROD_DEST is infinite and
   PROD_DEST_RATIO, BL_RATIO and FD are constant.

5. **Prediction.** SA_FIML with beta from the network, from the baseline.

       cp restart_baseline.csv restart_flow.csv
       mpirun -n 8 SU2_CFD config_predict.cfg
       cp restart_flow.csv restart_predicted.csv

   Compare J and the figures with steps 2 and 3.

6. **Weight optimization (optional).** The adjoint gives dJ/dw through the network; it starts from the trained
   weights and the predicted solution.

       cp nn_weights.dat nn_weights_trained.dat
       python3 ../fiml_inversion.py -f config_predict.cfg -n 8 --mode weights --x_scale 0.01 \
           --adjoint_iter 300 --adjoint_cfl 100 --maxiter 10 \
           --fields Misfit_Ucl Misfit_U Misfit_W Beta_Reg --log direct_history.csv

### Notes

- The adjoint at CFL 100 reaches its lowest rms[A_P] (-4.3) after 50 iterations and then rises slowly, -3.9 after
  300; rms[A_nu] keeps falling. Do not use `--warm_adjoint`: the driver stops when rms[A_P] ends above its first
  value, which a warm start near the minimum gives. 300 adjoint iterations take 7.5 min on 4 ranks.
- Gradient check at the baseline: beta = 1 -/+ 0.05 g / max|g|, g the adjoint sensitivity after 300 iterations, gives
  J = 2351.76 and 2099.77; the central difference -126.0 against -124.2 from the adjoint (ratio 1.015).
- The sensitivity dJ/dbeta is negative everywhere: more production gives more mixing and a lower J. It lies in the
  shear layer up to x/D = 6 and is zero beyond x/D = 15, where there are no data.
- The probes and the misfit functions need the custom outputs evaluated in the order probes, then functions
  (`CFlowOutput::SetCustomOutputs` on `feature_fiml`). With the order of develop, the functions use the probe values
  of the previous evaluation, so the recorded objective does not depend on the current solution.
