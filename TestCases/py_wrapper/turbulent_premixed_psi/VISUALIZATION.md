# Visualization of the Coupling Problem and Solution

## The Coupling Problem

```
┌─────────────────────────────────────────────────────────────┐
│                    ITERATION N                               │
└─────────────────────────────────────────────────────────────┘

WITHOUT FIX (BROKEN):
━━━━━━━━━━━━━━━━━━━━━
┌─────────────┐
│ Species Eq  │  ∂(ρc)/∂t + ∇·(ρcv) = S_source
│             │  → Solves for c with source term
│ c: 0.3→0.5  │  → Converges well ✓
└─────────────┘
       ↓
┌─────────────┐
│ Update T    │  T = (1-c)*Tu + c*Tf
│ T: 900→1100 │  → Temperature updates ✓
└─────────────┘
       ↓
       ✗ ← NO DENSITY UPDATE!
       ↓
┌─────────────┐
│ Flow Eqs    │  ρ still has OLD value (ρ_old = 2.1)
│ (next iter) │  But T is NEW (T_new = 1100 K)
│             │  → Violates ρ = P/(RT) ✗
│ Residuals:  │  → Momentum sees wrong inertia
│  Pressure ✗ │  → Continuity sees wrong mass flux
│  Velocity ✗ │  → STAGNATION!
└─────────────┘


WITH FIX (CORRECT):
━━━━━━━━━━━━━━━━━━━
┌─────────────┐
│ Species Eq  │  ∂(ρc)/∂t + ∇·(ρcv) = S_source
│             │  → Solves for c with source term
│ c: 0.3→0.5  │  → Converges well ✓
└─────────────┘
       ↓
┌─────────────┐
│ Update T    │  T = (1-c)*Tu + c*Tf
│ T: 900→1100 │  → Temperature updates ✓
│             │
│ Update ρ    │  ρ = P/(R*T_new)
│ ρ: 2.1→1.6  │  → Density updates ✓
└─────────────┘
       ↓
       ✓ ← CONSISTENT STATE!
       ↓
┌─────────────┐
│ Flow Eqs    │  ρ and T are consistent
│ (next iter) │  ρ*R*T = P ✓
│             │  → Momentum sees correct inertia
│ Residuals:  │  → Continuity sees correct mass flux
│  Pressure ✓ │  → CONVERGES!
│  Velocity ✓ │
└─────────────┘
```

## Physical Picture

```
UNBURNT REGION (c=0):          FLAME REGION (0<c<1):        BURNT REGION (c=1):
━━━━━━━━━━━━━━━━━              ━━━━━━━━━━━━━━━━━━━          ━━━━━━━━━━━━━━━━
T = 673 K                      T = 900-1500 K               T = 1777 K
ρ = 2.52 kg/m³                 ρ = 1.0-2.0 kg/m³            ρ = 0.96 kg/m³
P = 5 bar                      P = 5 bar                    P = 5 bar

Heavy, cold gas                STRONG SOURCE TERM HERE!     Light, hot gas
                               c changes rapidly →
                               T changes rapidly →
                               ρ MUST change rapidly →

                               ⚠️ If ρ not updated:
                               Flow equations blow up here!
```

## The Thermodynamic State

```
Ideal Gas Law:  P = ρ * R * T

┌──────────────────────────────────────────────────────────┐
│  WHAT HAPPENS WITHOUT THE FIX:                           │
│                                                           │
│  Iteration N:   P = ρ_old * R * T_old  ✓ (consistent)   │
│                 5 = 2.1   * R * 900    ✓                 │
│                                                           │
│  Update c → T changes, but ρ doesn't:                    │
│                                                           │
│  Iteration N+1: P ≠ ρ_old * R * T_new  ✗ (BROKEN!)      │
│                 5 ≠ 2.1   * R * 1100   ✗                 │
│                     ↑         ↑                           │
│                   still    updated                        │
│                    old!                                   │
│                                                           │
│  → Momentum equation uses ρ_old                          │
│  → Energy equation uses T_new                            │
│  → INCONSISTENT! Residuals can't converge                │
└──────────────────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────────┐
│  WHAT HAPPENS WITH THE FIX:                              │
│                                                           │
│  Iteration N:   P = ρ_old * R * T_old  ✓ (consistent)   │
│                 5 = 2.1   * R * 900    ✓                 │
│                                                           │
│  Update c → T changes → ρ also changes:                  │
│                                                           │
│  Iteration N+1: P = ρ_new * R * T_new  ✓ (CONSISTENT!)  │
│                 5 = 1.6   * R * 1100   ✓                 │
│                     ↑         ↑                           │
│                 updated   updated                         │
│                together!                                  │
│                                                           │
│  → All equations see consistent state                    │
│  → Residuals converge properly ✓                         │
└──────────────────────────────────────────────────────────┘
```

## Residual Behavior

```
WITHOUT FIX:                    WITH FIX:
━━━━━━━━━━━                    ━━━━━━━━━━

Species residual:               Species residual:
  │                               │
  │\                              │\
  │ \                             │ \
  │  \___                         │  \___
  │      ───___                   │      ───___
  └────────────→ iter            └────────────→ iter
  Converges ✓                     Converges ✓

Pressure residual:              Pressure residual:
  │                               │
  │\ /\/\/\/\                     │\
  │ ────────                      │ \
  │  STAGNATES!                   │  \___
  │      ✗                        │      ───___
  └────────────→ iter            └────────────→ iter
                                  Converges ✓

Velocity residual:              Velocity residual:
  │                               │
  │  /\/\/\/\                     │\
  │ ─────────                     │ \
  │  STAGNATES!                   │  \___
  │      ✗                        │      ───___
  └────────────→ iter            └────────────→ iter
                                  Converges ✓
```

## The Fix in One Picture

```
         c changes
            ↓
    ┌───────────────┐
    │  T = f(c)     │  ← Algebraic relation
    └───────────────┘
            ↓
    ┌───────────────┐
    │ ρ = P/(RT)    │  ← MUST ADD THIS!
    └───────────────┘
            ↓
    Consistent state for flow solver
```

## Key Insight

**In variable-density flows (like combustion), whenever you update
temperature, you MUST update density to maintain P = ρRT.**

This is true even for "incompressible" solvers when density varies
due to temperature (not pressure waves).

```
Temperature ↑  →  Density ↓  (at constant pressure)
   +63%               -62%    (for this test case)
```

**Forgetting this → Coupling breaks → Residuals stagnate → No convergence!**
