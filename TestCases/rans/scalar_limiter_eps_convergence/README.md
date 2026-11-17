# Turbulence limiter epsilon: convergence comparisons

These configurations reproduce the full convergence comparisons for PR #2978.
They complement the short `turb_naca0012_{sa,sst}_limiter_eps` MPI regressions
in `parallel_regression.py`, which retain `LIMITER_ITER=3` to exercise frozen
point limiting within ten iterations. These longer cases are manual benchmarks;
they do not add 10,000-iteration jobs to the regression suite.

The [recorded configurations, compressed full histories and convergence plot](https://github.com/rois1995/SU2/tree/pr-images/scalarLimiterEps)
include the three compared heads and the SST NACA0012 controls.

## Inputs and reproduction

Use the existing inputs from the [SU2 TestCases repository](https://github.com/su2code/TestCases):

| Configuration | Input directory | Required files | Initial condition |
|---|---|---|---|
| `sa_naca.cfg` | `rans/naca0012` | `n0012_225-65.su2`, `solution_flow_sa.dat` | Existing SA restart |
| `sst_flatplate.cfg` | `rans/flatplate` | `mesh_flatplate_turb_137x97.su2` | Free stream |
| `sst_naca.cfg` | `rans/naca0012` | `n0012_225-65.su2` | Free stream; limiter updated throughout |
| `sst_naca_frozen500.cfg` | `rans/naca0012` | `n0012_225-65.su2` | Free stream; limiter frozen at iteration 500 |

SHA-256 checksums of the inputs used:

```text
20dfc96875b1d1a8377c2e64248b6ee42a1ecf525c50d822fb5f539d06d68d15  n0012_225-65.su2
3e665e2d239f71c80c8625235be7e03a8fd5bd9440ef2d930e997291ed412fbc  solution_flow_sa.dat
8c720c21c0fc5a84f33565df9c0b0b361e4d3ac67d97fa96ab928783a52c6707  mesh_flatplate_turb_137x97.su2
```

Run each configuration in a separate directory, with its input files linked or
copied there. For example:

```bash
mkdir -p /tmp/su2-sst-frozen500
cd /tmp/su2-sst-frozen500
cp /path/to/SU2/TestCases/rans/scalar_limiter_eps_convergence/sst_naca_frozen500.cfg run.cfg
ln -s /path/to/TestCases/rans/naca0012/n0012_225-65.su2 .
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpirun -n 2 /path/to/SU2_CFD run.cfg
```

For the baseline comparisons, run the same configuration and inputs with each
head's executable. Recorded builds use GCC 9.4.0, release/skylake, MPI and OpenMP
enabled, two MPI ranks, and one OpenMP thread per rank. CGNS and TecIO are
disabled. Compare iteration counts; concurrent runs make their wall times
unsuitable for speedup claims.

## Results

Compared numerical source revisions:

- `develop`: `3aa1a95558f0af7ff43ccc3ec2fe403aff655811`
- PR #2956: `a3968c94950c18314a61e76703022883071cf655`
- PR #2978: `f1b5967a57cba9d1aa1eea9d33d33ca07531d6f6`

The PR #2956 executable was built at `eef97f4ae5`; its `Common`/`SU2_CFD`
sources are identical to `a3968c9495`. Subsequent changes updated test references
and unit-test attribution. The configurations added here preserve the numerical
settings of the recorded runs.

Counts include iteration zero. A stop below the **density residual** target is
not a claim that every flow and turbulence residual has converged.

| Configuration | Density target (log10 RMS) | develop | PR #2956 | PR #2978 |
|---|---:|---:|---:|---:|
| SA NACA0012 | -12 | 272 | 272 | 268 |
| SST flat plate | -12 | Not reached in 10,000 | Not reached in 10,000 | 758 |
| SST NACA0012, limiter updated throughout | -11 | Not reached in 10,000 | 467 | Not reached in 10,000 |
| SST NACA0012, limiter frozen at 500 | -11 | Not run | Not run | 620 |

The improvement is case dependent: the flat plate reaches its target, SA changes
little, and continuously updating the SST NACA0012 limiter stalls near a density
residual of -10 on the epsilon branch.

## SST NACA0012 and the current workaround

Both SST NACA0012 configurations have `MUSCL_TURB=YES`,
`SLOPE_LIMITER_TURB=VENKATAKRISHNAN`, `VENKAT_LIMITER_COEFF=0.03` and fixed
`CFL_NUMBER=1000`. Their only numerical difference is `LIMITER_ITER`:
99999 versus 500.

Diagnostics point to nonlinear cycling of the turbulence reconstruction in the
wake. In two consecutive stalled steps, the largest change of the k point
limiter is about 0.103 at x=2.284, y=0.063, beside the largest density residuals.
Nine reconstructed wake faces cross the lower k bound. On edge 7122--7123, the
reconstructed k changes from 2.57e-8 to -8.67e-9; the existing positivity guard
then switches that face to its first-order nodal value (about 4.50e-7).

This switching is an observed contributor to a discontinuous reconstruction,
not a proven complete explanation. The implicit convective Jacobian omits the
gradient/limiter dependence of the reconstructed states, which is a plausible
contributor to the cycle. Lowering CFL to 100 did not reach the target in 10,000
iterations. A more accurate linear solve and isolated continuous face-clipping
experiments also failed to remove the stalled cycle within their tested budgets;
no general production correction has been validated.

For this case, the current workaround is to **freeze the stored point limiter**:

```text
MUSCL_TURB= YES
SLOPE_LIMITER_TURB= VENKATAKRISHNAN
LIMITER_ITER= 500
```

The stored limiter continues to be applied after freezing, and turbulence MUSCL
remains enabled. This run reaches the density target in 620 iterations
(log10 RMS density = -11.0072, CL = 1.070484, CD = 0.0157461). The freeze
iteration is a demonstrated setting for this configuration, rather than a
validated universal choice. `VAN_ALBADA_EDGE` is evaluated at edges and is not
frozen by `LIMITER_ITER`.
