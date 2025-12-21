# Quick Fix Summary: Residual Stagnation Solution

## The Problem
```
Species equation → converges ✓
Update T from species → T = (1-c)*Tu + c*Tf
Flow equations → STAGNATE ✗
```

**Why?** Density wasn't updated when temperature changed!

## The Fix

### BEFORE (Broken):
```python
def update_temperature(SU2Driver, iPoint):
    C = SU2Driver.Solution(iSPECIESSOLVER)(iPoint,0)
    T = Tu*(1-C) + Tf*C
    Cp = SU2Driver.Primitives()(iPoint,iCp)
    SU2Driver.Solution(iFLOWSOLVER).Set(iPoint,iTEMP,Cp*T)
    # ⚠️ Density not updated! ρ still has old value
```

### AFTER (Fixed):
```python
def update_temperature(SU2Driver, iPoint):
    C = SU2Driver.Solution(iSPECIESSOLVER)(iPoint,0)
    T_target = Tu*(1-C) + Tf*C
    T_old = SU2Driver.Primitives()(iPoint,iTEMPERATURE)

    # ✓ Under-relaxation for stability
    T_new = RELAX_FACTOR * T_target + (1.0 - RELAX_FACTOR) * T_old

    # ✓ Compute gas constant from old state
    R_gas = P / (rho_old * T_old)

    # ✓ Update density consistently: ρ = P/(RT)
    rho_new = P / (R_gas * T_new)

    # ✓ Update both temperature and density
    SU2Driver.Solution(iFLOWSOLVER).Set(iPoint, iTEMP, T_new)
    SU2Driver.Primitives().Set(iPoint, iDENSITY, rho_new)
```

## Main Loop Sequence

### BEFORE:
```python
driver.Preprocess()
driver.Run()
# Set source terms (WRONG: after Run!)
for i in nodes: Source.Set(i, zimont(i))
# Update temperature (no density update)
for i in nodes: update_temperature(i)
driver.Postprocess()
```

### AFTER:
```python
driver.Preprocess()
# ✓ Set source terms BEFORE Run
for i in nodes: Source.Set(i, zimont(i))
# ✓ Run solvers
driver.Run()
# ✓ Update temperature AND density AFTER Run
for i in nodes: update_temperature(i)
# ✓ Monitor coupling
print(f"Max ΔT: {max_delta_T}, Max Δρ/ρ: {max_delta_rho}")
driver.Postprocess()
```

## Key Parameters

```python
# Add at top of file:
RELAX_FACTOR = 0.7  # Adjust 0.5-0.8 for stability vs speed
```

## Expected Behavior

| Before Fix | After Fix |
|------------|-----------|
| Species: ✓ converging | Species: ✓ converging |
| Pressure: ✗ stagnating | Pressure: ✓ converging |
| Velocity: ✗ stagnating | Velocity: ✓ converging |
| ρ inconsistent with T | ρ = P/(RT) ✓ |

## Physics

**In combustion:**
- c: 0 → 1 (progress)
- T: 673K → 1777K ↑
- ρ: 2.52 → ~0.95 ↓ (MUST decrease!)

**If ρ not updated:**
- Continuity: wrong mass fluxes
- Momentum: wrong inertia
- → Residuals can't converge!

## Monitoring

Look for in output:
```
Source term - Max: 1.23e+02, Min: 0.00e+00, Avg: 3.45e+01
Max temperature change: 5.23e+01 K
Max relative density change: 2.15e-02
```

## Troubleshooting

**Still not converging?**
1. Reduce `RELAX_FACTOR` to 0.5
2. Reduce CFL number in config
3. Check source term isn't too strong
4. Verify mass conservation: Σ(source terms) ≈ 0

**Oscillations?**
1. Reduce `RELAX_FACTOR`
2. Check grid quality in flame region
3. Ensure boundary conditions are correct

## One-Line Summary

**Always update density when temperature changes: `ρ_new = P/(R*T_new)`**
