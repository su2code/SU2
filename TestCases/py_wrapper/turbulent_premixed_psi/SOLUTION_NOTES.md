# Solution for Residual Stagnation in Turbulent Premixed Combustion

## Problem Description

When using the Python wrapper to add species transport with strong source terms and algebraically updating temperature/enthalpy from species composition, the pressure and velocity residuals stagnate in regions with high source terms while the species residuals converge properly.

## Root Cause

The issue occurs because:

1. **Species solver** converges and updates the progress variable `c`
2. **Temperature is updated algebraically**: `T = (1-c)*Tu + c*Tf`
3. **Density is NOT updated** to match the new temperature
4. **Thermodynamic inconsistency**: The ideal gas law `ρ = P/(RT)` is violated
5. **Flow equations see inconsistent state**: Continuity and momentum equations operate with outdated density

This breaks the coupling between species, energy, and momentum equations, causing residual stagnation.

## Solution Implemented (Solution #1)

### Key Changes in `run.py`:

#### 1. **Updated `update_temperature()` function** (Lines ~90-130)
   - Now computes new density when temperature changes: `ρ_new = P/(R*T_new)`
   - Updates both temperature AND density in the primitive variables
   - Added under-relaxation factor for stability: `T_new = ω*T_target + (1-ω)*T_old`
   - Maintains thermodynamic consistency with ideal gas law

#### 2. **Improved iteration sequence** (Lines ~250-290)
   - **STEP 1**: Set species source terms
   - **STEP 2**: Run all solvers (species equation solved)
   - **STEP 3**: Update temperature AND density from new species field
   - **STEP 4**: Postprocess and update
   - Ensures proper coupling between species and flow equations

#### 3. **Added under-relaxation** (Line ~58)
   - `RELAX_FACTOR = 0.7` prevents oscillations
   - Can be tuned between 0.5-0.8 depending on source term strength
   - Lower values = more stable but slower convergence

#### 4. **Added coupling diagnostics** (Lines ~270-280)
   - Monitors `max_delta_T` and `max_delta_rho`
   - Helps identify coupling issues
   - Shows when thermodynamic updates are stabilizing

#### 5. **Added source term monitoring** (Lines ~185-215, ~260-265)
   - Tracks max/min/avg source term values
   - Helps diagnose convergence problems
   - Important for checking mass conservation

## Physical Interpretation

For a premixed flame:
- Progress variable goes from `c=0` (unburnt) to `c=1` (burnt)
- Temperature increases: `T: 673K → 1777K`
- Density MUST decrease (isobaric): `ρ ∝ 1/T`
- Without density update: continuity equation sees wrong mass fluxes
- Result: momentum residuals can't converge

## Expected Results

After this fix:
- ✅ Species residuals converge (as before)
- ✅ Temperature correctly updates from species
- ✅ **Density correctly updates from temperature**
- ✅ **Pressure residuals converge** (no longer stagnate)
- ✅ **Velocity residuals converge** (no longer stagnate)
- ✅ Thermodynamic consistency maintained

## Tuning Parameters

If convergence is still slow, adjust:

### 1. **Under-relaxation factor** (Line 58)
```python
RELAX_FACTOR = 0.5-0.8  # Lower = more stable, slower
```

### 2. **CFL number in config file**
```
CFL_NUMBER= 10.0  # Reduce if needed
CFL_ADAPT= YES
CFL_ADAPT_PARAM= ( 0.3, 0.5, 10.0, 100.0, 1e-6 )
```

### 3. **Number of inner iterations**
```python
for inner_iter in range(10):  # Increase for more coupling iterations
```

## Additional Recommendations

1. **Monitor residual ratios**: Watch for regions where flow residuals are high relative to species residuals

2. **Check mass conservation**: The source term should not create/destroy mass:
   ```
   ∂(ρc)/∂t + ∇·(ρcv) = S_c
   ```
   where `S_c` is per unit volume

3. **Verify ideal gas law**: At any point, check if `P = ρRT` holds after updates

4. **Use implicit time stepping**: Helps with stiff source terms (already in config)

5. **Consider local time stepping**: Can help in regions with varying source terms

## Technical Details

### Incompressible Flow Considerations
This case uses `INC.FLOW` (incompressible solver), which:
- Solves for pressure directly (not density-based)
- Assumes low Mach number
- Still requires density updates for variable density flows (like combustion!)
- Energy equation is decoupled but affects momentum through density

### Variable Density Incompressible Flow
Even though it's "incompressible", combustion creates variable density:
- ∇·(ρv) = 0 (mass conservation, not ∇·v = 0)
- ρ varies due to temperature, not pressure waves
- Density MUST be updated when temperature changes

## References
- See main SU2 documentation for species transport
- Zimont turbulent flame speed closure: Zimont, V.L., 2000. "Gas premixed combustion at high turbulence."
- Variable density incompressible flow formulation

## Contact
For issues or questions about this fix, refer to the SU2 forums or GitHub issues.
