---
title: Adding a Transition Model
permalink: /docs_v7/Adding-a-Transition-Model/
---

This page shows how the transition model of SU2 is built, how it is coupled with the turbulence models, and what to change to add a correlation, a variant or a new transition model. The code examples are copied from the SU2 source.

SU2 has one transport-equation transition model: the Langtry-Menter (LM) model. It works with SA and with SST. It has two equations:

- the intermittency, which is near 0 in laminar regions and 1 in turbulent regions;
- the transition momentum thickness Reynolds number (`Re_theta_t`), the Reynolds number at which transition starts, based on the momentum thickness of the boundary layer.

The SA option `SA_OPTIONS= BCM` (Bas-Cakmakcioclu) is different: it is an algebraic transition term inside the SA source term. It is described in [Adding a Turbulence Model](../Adding-a-Turbulence-Model/), not here.

Read [Adding a Turbulence Model](../Adding-a-Turbulence-Model/) first. The transition solver uses the same base classes, edge kernels and source term structure.

## Files you will touch

| File | What goes there |
|---|---|
| `Common/include/option_structure.hpp` | `TURB_TRANS_MODEL` and `Trans_Model_Map`, `LM_OPTIONS` and `LM_Options_Map`, `TURB_TRANS_CORRELATION`, `LM_ParsedOptions`, `ParseLMOptions` |
| `Common/include/CConfig.hpp`, `Common/src/CConfig.cpp` | `KIND_TRANS_MODEL`, `LM_OPTIONS`, `HROUGHNESS`, the parsing and checks in `CConfig::SetPostprocessing`, the summary in `CConfig::SetOutput` |
| `SU2_CFD/include/solvers/CTransLMSolver.hpp`, `SU2_CFD/src/solvers/CTransLMSolver.cpp` | The solver |
| `SU2_CFD/include/variables/CTransLMVariable.hpp`, `SU2_CFD/src/variables/CTransLMVariable.cpp` | The variables, with the separation and effective intermittency |
| `SU2_CFD/include/numerics/turbulent/transition/trans_sources.hpp` | The source terms |
| `SU2_CFD/include/numerics/turbulent/transition/trans_correlations.hpp` | The correlations (`TransLMCorrelations`) |
| `SU2_CFD/include/numerics/turbulent/transition/trans_edge_flux.hpp` | Convection and diffusion (`CScalarFlux_TransLM`) |
| `SU2_CFD/src/solvers/CTurbSASolver.cpp`, `CTurbSSTSolver.cpp`, `SU2_CFD/src/variables/CTurbSSTVariable.cpp`, `SU2_CFD/include/numerics/turbulent/turb_sources.hpp` | The coupling with the turbulence model |
| `SU2_CFD/src/solvers/CSolverFactory.cpp`, `SU2_CFD/src/drivers/CDriver.cpp`, `SU2_CFD/src/iteration/CFluidIteration.cpp` | Creation of the solver and numerics, order of the solvers |
| `SU2_CFD/src/output/CFlowOutput.cpp` | Output fields |

---

- [Overview](#overview)
- [Options](#options)
- [Solver and variable classes](#solver-and-variable-classes)
- [Source terms](#source-terms)
- [Correlations](#correlations)
- [Convection and diffusion](#convection-and-diffusion)
- [Coupling with the turbulence model](#coupling-with-the-turbulence-model)
- [Order of the solvers](#order-of-the-solvers)
- [Boundary conditions](#boundary-conditions)
- [Free-stream and initial values](#free-stream-and-initial-values)
- [Output and restart](#output-and-restart)
- [Discrete adjoint](#discrete-adjoint)
- [Tests](#tests)
- [How to add a correlation, a variant or a new model](#how-to-add-a-correlation-a-variant-or-a-new-model)
- [Checklist](#checklist)

---

## Overview

| Item | LM in SU2 |
|---|---|
| Solver | [`CTransLMSolver`](../../doxygen/classCTransLMSolver.html), at position `TRANS_SOL` of the solver container, equation system `RUNTIME_TRANS_SYS` |
| Variables | [`CTransLMVariable`](../../doxygen/classCTransLMVariable.html): intermittency (index 0) and transition momentum thickness Reynolds number (index 1) |
| Source term | [`CSourcePieceWise_TransLM`](../../doxygen/classCSourcePieceWise__TransLM.html) |
| Correlations | [`TransLMCorrelations`](../../doxygen/classTransLMCorrelations.html) |
| Convection and diffusion | [`CScalarFlux_TransLM`](../../doxygen/classCScalarFlux__TransLM.html) |
| Numerical settings | the turbulence ones: `CONV_NUM_METHOD_TURB`, `MUSCL_TURB`, `SLOPE_LIMITER_TURB`, `TIME_DISCRE_TURB`, `CFL_REDUCTION_TURB` |

The transition model has no numerical options of its own. `CConfig::SetGlobalParam` applies the turbulence settings for `RUNTIME_TRANS_SYS`.

## Options

Two options control the model:

- `KIND_TRANS_MODEL` chooses the model. Its enum is `TURB_TRANS_MODEL`, its map is `Trans_Model_Map`, and its values are `NONE` and `LM`.
- `LM_OPTIONS` chooses the variants and the correlation. It is an enum list option, so the user can give several keywords:

```cpp
enum class LM_OPTIONS {
  NONE,         /*!< \brief No option / default. */
  LM2015,       /*!< \brief Cross-flow corrections. */
  MALAN,        /*!< \brief Kind of transition correlation model (Malan). */
  SULUKSNA,     /*!< \brief Kind of transition correlation model (Suluksna). */
  KRAUSE,       /*!< \brief Kind of transition correlation model (Krause). */
  KRAUSE_HYPER, /*!< \brief Kind of transition correlation model (Krause hypersonic). */
  MEDIDA_BAEDER,/*!< \brief Kind of transition correlation model (Medida-Baeder). */
  MEDIDA,       /*!< \brief Kind of transition correlation model (Medida). */
  MENTER_LANGTRY,   /*!< \brief Kind of transition correlation model (Menter-Langtry). */
  DEFAULT       /*!< \brief Kind of transition correlation model (Menter-Langtry if SST, MALAN if SA). */
};
```

`LM_Options_Map` has the keywords `NONE`, `LM2015`, `MALAN`, `SULUKSNA`, `KRAUSE`, `KRAUSE_HYPER`, `MEDIDA_BAEDER`, `MENTER_LANGTRY` and `DEFAULT`. `MEDIDA` is in the enum but not in the map, so users cannot select it.

The correlation has its own enum, `TURB_TRANS_CORRELATION`, with the same names as the correlation keywords. The parsed options struct holds the result of reading the list:

```cpp
struct LM_ParsedOptions {
  LM_OPTIONS version = LM_OPTIONS::NONE;  /*!< \brief LM base model. */
  bool LM2015 = false;                    /*!< \brief Use cross-flow corrections. */
  TURB_TRANS_CORRELATION Correlation = TURB_TRANS_CORRELATION::DEFAULT;
};
```

`ParseLMOptions(LM_Options, nLM_Options, rank, Kind_Turb_Model)` does four things:

1. It sets `LM2015`.
2. It sets `Correlation` from the correlation keyword.
3. It stops if more than one correlation is given.
4. It replaces `DEFAULT` by the default of the turbulence model:

```cpp
  if (LMParsedOptions.Correlation == TURB_TRANS_CORRELATION::DEFAULT){
    if (Kind_Turb_Model == TURB_MODEL::SST) {
      LMParsedOptions.Correlation = TURB_TRANS_CORRELATION::MENTER_LANGTRY;
    } else if (Kind_Turb_Model == TURB_MODEL::SA) {
      LMParsedOptions.Correlation = TURB_TRANS_CORRELATION::MALAN;
    }
  }
```

`ParseLMOptions` does not set `version`. It stays `LM_OPTIONS::NONE`.

`CConfig::SetPostprocessing` calls the parser and checks the options:

```cpp
  if (Kind_Turb_Model == TURB_MODEL::NONE && Kind_Trans_Model != TURB_TRANS_MODEL::NONE) {
    SU2_MPI::Error("KIND_TURB_MODEL cannot be NONE to use a transition model", CURRENT_FUNCTION);
  }
  // ...
  /*--- Postprocess LM_OPTIONS into structure. ---*/
  if (Kind_Trans_Model == TURB_TRANS_MODEL::LM) {
    lmParsedOptions = ParseLMOptions(LM_Options, nLM_Options, rank, Kind_Turb_Model);

    /*--- Check if problem is 2D and LM2015 has been selected ---*/
    if (lmParsedOptions.LM2015 && val_nDim == 2) {
      SU2_MPI::Error("LM2015 is available only for 3D problems", CURRENT_FUNCTION);
    }
  }
```

The code reads the result with `config->GetLMParsedOptions()`. `CConfig::SetOutput` prints the model (2009 or 2015) and the correlation. `HROUGHNESS` (`GethRoughness()`) is the wall roughness used by the `LM2015` cross-flow term.

New variants of the LM model should use the same mechanism: a keyword in `LM_OPTIONS` and a field in `LM_ParsedOptions`.

## Solver and variable classes

The transition model has a solver class and a variable class, built like those of the turbulence models.

`CTransLMSolver` derives from `CTurbSolver` with `conservative = true`. The transported quantities are the density times the intermittency and the density times the Reynolds number. Its constructor (`SU2_CFD/src/solvers/CTransLMSolver.cpp`) follows the turbulence solvers: `nVar = 2`, the linear system, the limits, `Solution_Inf`, `nodes` and `ghostNodes`, `Inlet_TurbVars`, the local CFL. It also stores the options, the correlations and the family of the turbulence model:

```cpp
  /*--- Define variables needed for transition from config file */
  options = config->GetLMParsedOptions();
  TransCorrelations.SetOptions(options);
  TurbFamily = TurbModelFamily(config->GetKind_Turb_Model());
```

The implicit update limits the intermittency to `[1e-4, 5.0]` and the Reynolds number to `[1e-4, 1e15]`.

`CTransLMVariable` derives from `CTurbVariable`. It adds two per-point fields: the separation intermittency `Intermittency_Sep` and the effective intermittency `Intermittency_Eff`. The effective value is computed when it is set:

```cpp
void CTransLMVariable::SetIntermittencyEff(unsigned long iPoint, su2double val_Intermittency_sep) {

  /*--- Effective intermittency ---*/
  Intermittency_Eff(iPoint) = max(Solution(iPoint,0), val_Intermittency_sep);

}
```

`CTransLMSolver::Postprocessing` first computes the gradients. Then, at each point, it computes the separation intermittency from the flow, the turbulence and the correlations (`ReThetaC_Correlations`). It stores both values:

```cpp
    su2double Intermittency_Sep = 2.0*max(0.0, Re_v/(3.235*Corr_Rec)-1.0)*f_reattach;
    Intermittency_Sep = min(Intermittency_Sep,2.0)*f_theta;
    Intermittency_Sep = min(max(0.0, Intermittency_Sep), 2.0);
    nodes->SetIntermittencySep(iPoint, Intermittency_Sep);
    nodes->SetIntermittencyEff(iPoint, Intermittency_Sep);
```

The solver also implements `Preprocessing` (which calls `CommonPreprocessing`), `Upwind_Residual`, `Source_Residual`, the boundary conditions and `LoadRestart`.

`CSolverFactory::CreateTransSolver` creates the solver for `RANS` and `INC_RANS`. It then runs the flow `Preprocessing`, the transition `Postprocessing` (which sets the effective intermittency) and the flow `Preprocessing` again:

```cpp
      case TURB_TRANS_MODEL::LM :
        transSolver = new CTransLMSolver(geometry, config, solver[FLOW_SOL], iMGLevel);
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
        transSolver->Postprocessing(geometry, solver, config, iMGLevel);
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
        break;
```

## Source terms

`CSourcePieceWise_TransLM<FlowIndices>` (`trans_sources.hpp`) computes the sources of both equations at one point. `CDriver::InstantiateTransitionNumerics` creates it for every multigrid level:

```cpp
    if (LM) trans_source_first_term = new CSourcePieceWise_TransLM<Indices>(nDim, nVar_Trans, config);
```

`CTransLMSolver::Source_Residual` sets its inputs. The turbulence solution goes in the "scalar" slots, and the transition solution goes in the "trans" slots:

```cpp
    numerics->SetScalarVar(turbNodes->GetSolution(iPoint), nullptr);
    numerics->SetScalarVarGradient(turbNodes->GetGradient(iPoint), nullptr);

    /*--- Transition variables w/o reconstruction, and its gradient ---*/

    numerics->SetTransVar(nodes->GetSolution(iPoint), nullptr);
    numerics->SetTransVarGradient(nodes->GetGradient(iPoint), nullptr);
```

It also sets:

- the flow primitives and their gradients;
- the volume, the wall distance and the coordinates;
- the vorticity and the strain rate magnitude;
- for `LM2015`, the local grid length (`SetLocalGridLength(geometry->nodes->GetMaxLength(iPoint))`).

The result is subtracted from the residual, as for the turbulence models.

`ComputeResidual` follows the LM model of Langtry and Menter (2009), with the `LM2015` cross-flow terms. Some terms depend on the turbulence model. They branch on `TurbFamily`:

| Quantity | SST (`TURB_FAMILY::KW`) | SA (`TURB_FAMILY::SA`) |
|---|---|---|
| Turbulence intensity `Tu` (in percent) | from the local k and velocity, at least 0.027 | `FREESTREAM_TURBULENCEINTENSITY` |
| `R_t` | `rho*k/(mu*omega)` | `mu_t/mu` |
| `F_length` | correlation blended with 40 near the wall (`f_sub`) | correlation |
| `F_onset2`, `F_onset3` | limits 2.0 and `1 - (R_t/2.5)^3` | limits 4.0 and `2 - (R_t/2.5)^3` |
| `f_wake` | from `re_omega` | 1 |

The source of the intermittency is production `Pg` minus destruction `Dg`. The source of the Reynolds number is `PRethetat` minus `DRethetat`. `DRethetat` is non-zero only with `LM2015`:

```cpp
      /*--- Source ---*/
      Residual[0] += (Pg - Dg) * Volume;
      Residual[1] += (PRethetat - DRethetat) * Volume;
```

`ComputeResidual` also computes the equilibrium Reynolds number `Corr_Ret`. It uses a fixed-point loop on the pressure gradient parameter `lambda`, with at most 100 iterations. The Jacobian is diagonal: `Jacobian_i[0][0]` and `Jacobian_i[1][1]`.

## Correlations

The correlations give two quantities from the turbulence intensity and the transported Reynolds number:

- `ReThetaC_Correlations` gives the critical Reynolds number;
- `FLength_Correlations` gives the transition length function.

Both are functions of `TransLMCorrelations` (`trans_correlations.hpp`). This class holds a copy of `LM_ParsedOptions` (`SetOptions`), and both functions switch on `options.Correlation`:

```cpp
  su2double ReThetaC_Correlations(const su2double Tu, const su2double Re_theta_t) const {
  // ...
  su2double FLength_Correlations(const su2double Tu, const su2double Re_theta_t) const {
```

`CSourcePieceWise_TransLM::ComputeResidual` uses them. So does `CTransLMSolver::Postprocessing`, for the separation intermittency. The `DEFAULT` case stops with an error, because `ParseLMOptions` always replaces `DEFAULT`.

## Convection and diffusion

The edge kernel of the transition model is the smallest one of the scalar models. It inherits the conservative upwind convection of `CUpwScalarFlux`. It only gives the diffusion coefficients: `mu + mu_t` for the intermittency, and `2(mu + mu_t)` for the Reynolds number:

```cpp
template <class Double, class FlowIndices, int nDim, size_t nVar = 2>
class CScalarFlux_TransLM
    : public CUpwScalarBase<Double, CScalarFlux_TransLM<Double, FlowIndices, nDim, nVar>, FlowIndices, nDim, nVar> {
 public:
  static constexpr bool Conservative = true;
  static constexpr bool DiagonalDiffusion = true;
  // ...
    const Double diff_i_gamma = mu_i + muT_i;
    const Double diff_j_gamma = mu_j + muT_j;
    const Double diff_i_ReThetaT = 2.0 * (mu_i + muT_i);
    const Double diff_j_ReThetaT = 2.0 * (mu_j + muT_j);

    Vector<Double, nVar> D;
    D(0) = 0.5 * (diff_i_gamma + diff_j_gamma);
    D(1) = 0.5 * (diff_i_ReThetaT + diff_j_ReThetaT);
    return {D, D};
```

`CTransLMSolver::Upwind_Residual` runs it with `DispatchScheme<CScalarFlux_TransLM, 2>` and `ScalarFluxOptions::Interior(*config, config->GetBounded_Turb())`. [Adding a Turbulence Model](../Adding-a-Turbulence-Model/#5-convection-and-diffusion-the-edge-kernel) explains the kernel structure.

## Coupling with the turbulence model

The transition and turbulence models exchange data in both directions.

**From the turbulence model to the transition model:**

- the turbulence solution: k and omega for SST, through `SetScalarVar` in `Source_Residual` and `turbNodes->GetSolution` in `Postprocessing`;
- the eddy viscosity. `Postprocessing` reads it with `turbNodes->GetmuT`. The source term reads it from the flow primitives (`V_i[idx.EddyViscosity()]`), and so does the edge kernel.

**From the transition model to the turbulence model:** the effective intermittency. SA also gets the transported intermittency.

SST, in `CTurbSSTSolver::Source_Residual`:

```cpp
    /*--- Effective Intermittency ---*/
    if (config->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      numerics->SetIntermittencyEff(solver_container[TRANS_SOL]->GetNodes()->GetIntermittencyEff(iPoint));
    }
```

Then `CSourcePieceWise_TurbSST::ComputeResidual` scales the production and destruction of k:

```cpp
      /*--- LM model coupling with production and dissipation term for k transport equation---*/
      if (config->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
        pk = pk * eff_intermittency;
        dk = min(max(eff_intermittency, 0.1), 1.0) * dk;
      }
```

SST also changes the blending function `F1` in `CTurbSSTVariable::SetBlendingFunc`. `CTurbSSTSolver::Postprocessing` calls this function with the transition model:

```cpp
  /*--- LM model for F1 ---*/
  if (trans_model == TURB_TRANS_MODEL::LM) {
    su2double Ry = val_density*val_dist*sqrt(Solution(iPoint,0))/val_viscosity;
    su2double F3 = exp(-pow(Ry/120.0, 8.0));
    F1(iPoint) = max(F1(iPoint), F3);
  }
```

SA, in `CTurbSASolver::Source_Residual`:

```cpp
    if (config->GetKind_Trans_Model() != TURB_TRANS_MODEL::NONE) {
      numerics->SetIntermittencyEff(solver_container[TRANS_SOL]->GetNodes()->GetIntermittencyEff(iPoint));
      numerics->SetIntermittency(solver_container[TRANS_SOL]->GetNodes()->GetSolution(iPoint, 0));
    }
```

Then, in `CSourceBase_TurbSA::ComputeResidual`, the effective intermittency multiplies the production. The transported intermittency limits the destruction:

```cpp
      } else if (transition_LM){

        var.intermittency = intermittency_eff_i;
        // Is wrong the reference from NASA?
        // Original max(min(gamma, 0.5), 1.0) always gives 1 as result.
        var.interDestrFactor = min(max(intermittency_i, 0.5), 1.0);
```

Both turbulence solvers store the effective intermittency they used, with `nodes->SetIntermittency`. With a transition model, their `Postprocessing` also computes the turbulence index at the walls (`SetTurbIndex`). It is written as the output field `TURB_INDEX`.

The flow solver does not read the transition variables. The transition model acts on the flow only through the eddy viscosity of the turbulence model.

## Order of the solvers

In each inner iteration, the transition model runs after the flow and **before** the turbulence model. `CFluidIteration::Iterate` solves the flow. Then `CFluidIteration::CommonAuxiliarySolvers` solves the transition model and the turbulence model:

```cpp
    /*--- Solve transition model ---*/

    if (config[val_iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_TRANS_SYS);
      integration[val_iZone][val_iInst][TRANS_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                         RUNTIME_TRANS_SYS, val_iZone, val_iInst);
    }

    /*--- Solve the turbulence model ---*/
```

So:

1. The flow uses the eddy viscosity of the previous iteration.
2. The transition model uses the new flow and the turbulence of the previous iteration. At the end of its iteration, `CTransLMSolver::Postprocessing` updates the effective intermittency.
3. The turbulence model uses the new effective intermittency. It computes the eddy viscosity for the next flow iteration.

The transition model is solved on the finest grid only. Unlike the turbulence solution, its solution is not restricted (transferred) to the coarse grids: `CSingleGridIntegration::SingleGrid_Iteration` restricts only `RUNTIME_TURB_SYS`. For unsteady problems, `CFluidIteration::Update` calls `SetDualTime_Solver` for `TRANS_SOL` after `TURB_SOL`.

## Boundary conditions

`CTransLMSolver` implements the conditions in the table. All of them write the outside state into `ghostNodes`. Then they call the shared boundary flux (`BoundaryFlux`) with the convective term only.

| Function | Outside state |
|---|---|
| `BC_HeatFlux_Wall` (and `BC_Isothermal_Wall`, which calls it) | the free-stream state; the comment in the code calls it "wall-normal zero flux, convective only" |
| `BC_Inlet` | `Inlet_TurbVars` (initialized to the free-stream state) |
| `BC_Outlet` | calls `BC_Far_Field` |
| `BC_Far_Field` | the free-stream state |
| `BC_Fluid_Interface` | the donor state (`FluidInterfaceFluxResidual`) |

The other conditions come from the base classes:

- symmetry planes, Euler walls, periodic and Riemann boundaries from `CScalarSolver`;
- the Giles and turbomachinery Riemann boundaries from `CTurbSolver`.

The transition solver does not override `SetInletAtVertex` or `SetUniformInlet`. For new boundary conditions, see [Adding a Boundary Condition](../Adding-a-Boundary-Condition/).

## Free-stream and initial values

`CTransLMSolver` sets the free-stream state in its constructor and initializes all points with it. The intermittency is fixed at 1. The Reynolds number comes from the free-stream turbulence intensity (`FREESTREAM_TURBULENCEINTENSITY`, in percent here):

```cpp
  const su2double Intensity = config->GetTurbulenceIntensity_FreeStream()*100.0;

  const su2double Intermittency_Inf  = 1.0;
  su2double ReThetaT_Inf = 100.0;

  /*--- Momentum thickness Reynolds number, initialized from freestream turbulent intensity*/
  if (Intensity <= 1.3) {
    if(Intensity >=0.027) {
      ReThetaT_Inf = (1173.51-589.428*Intensity+0.2196/(Intensity*Intensity));
    }
    else {
      ReThetaT_Inf = (1173.51-589.428*Intensity+0.2196/(0.27*0.27));
    }
  }
  else if(Intensity>1.3) {
    ReThetaT_Inf = 331.5*pow(Intensity-0.5658,-0.671);
  }

  Solution_Inf[0] = Intermittency_Inf;
  Solution_Inf[1] = ReThetaT_Inf;
```

`CEulerSolver::SetNondimensionalization` computes the same Reynolds number and stores it with `config->SetReThetaT_FreeStream`. This value and `FREESTREAM_INTERMITTENCY` (`GetIntermittency_FreeStream()`) are only printed in the start-up table. The solver does not read them.

The turbulence model sets its own free-stream values (see [Adding a Turbulence Model](../Adding-a-Turbulence-Model/#free-stream-initial-and-inlet-values)).

## Output and restart

The output fields are in `CFlowOutput` (`SU2_CFD/src/output/CFlowOutput.cpp`). They are in the `switch (config->GetKind_Trans_Model())` blocks of the scalar output functions:

| Kind | Fields |
|---|---|
| History | `RMS_INTERMITTENCY`, `RMS_RE_THETA_T`, `MAX_INTERMITTENCY`, `MAX_RE_THETA_T`, `BGS_INTERMITTENCY`, `BGS_RE_THETA_T`, `LINSOL_ITER_TRANS`, `LINSOL_RESIDUAL_TRANS` |
| Volume, `SOLUTION` group | `INTERMITTENCY`, `RE_THETA_T` |
| Volume, `PRIMITIVE` group | `INTERMITTENCY_SEP`, `INTERMITTENCY_EFF`, `TURB_INDEX` |
| Volume, `RESIDUAL` group | `RES_INTERMITTENCY`, `RES_RE_THETA_T` |

To add more, see [Adding Output Fields](../Adding-Output-Fields/).

**Restart.** In the `SOLUTION` group, the transition fields come after the turbulence fields (`CFlowOutput::SetVolumeOutputFieldsScalarSolution`). `CTransLMSolver::LoadRestart` skips the coordinates, the flow variables and the turbulence variables:

```cpp
    unsigned short skipVars = nDim + solver[MESH_0][FLOW_SOL]->GetnVar() + solver[MESH_0][TURB_SOL] ->GetnVar();
```

Then, when no species model is active, it runs the flow `Preprocessing`, the turbulence `Postprocessing` and the transition `Postprocessing`. `CTurbSolver::LoadRestart` skips this step when a transition model is active. This way the step runs once, after both solutions are read.

## Discrete adjoint

The transition model has no discrete adjoint. For `DISC_ADJ_RANS`, `CSolverFactory::CreateSolverContainer` creates no `TRANS_SOL`. `CDiscAdjFluidIteration` handles only the flow, turbulence and species solvers.

Adding the adjoint would need:

- the transition solver in those two places;
- a `CDiscAdjSolver` for it;
- all preaccumulation inputs in `CSourcePieceWise_TransLM` and in the coupling terms of the turbulence sources. Preaccumulation is explained in [Adding a Turbulence Model](../Adding-a-Turbulence-Model/#discrete-adjoint).

## Tests

These tests cover the LM model:

| Case | Where | Script |
|---|---|---|
| T3A and T3A- flat plates (`transitional_LM_model_ConfigFile.cfg`), E387 airfoil with SA-LM and SST-LM | `Tutorials` repository, `compressible_flow/Transitional_Flat_Plate/Langtry_and_Menter` and `compressible_flow/Transitional_Airfoil/Langtry_and_Menter` | `tutorials.py` |
| E387 SST-LM, 2 iterations, sanitizer runs only (`enabled_with_regular = False`) | same configuration as the tutorial | `serial_regression.py` (`tutorial_trans_e387_sst_asan`) |

The folder `TestCases/transition` holds SA-BCM cases (`transitional_BC_model_ConfigFile.cfg`), not LM. `TestCases/rans/s809/trans_s809.cfg` uses SA-LM, but no test script runs it.

To add a test, follow [Running Regression Tests](../Running-Regression-Tests/#add-a-new-test). For a new correlation or variant, add an entry to `parallel_regression.py`. Use a small LM case with the new keyword in `LM_OPTIONS`, and one configuration file per test. Before you open the pull request, check the transition location against experiments, for example the T3A flat plate.

## How to add a correlation, a variant or a new model

### A new correlation

1. Add the value to `LM_OPTIONS`, to `LM_Options_Map` and to `TURB_TRANS_CORRELATION` in `option_structure.hpp`.
2. In `ParseLMOptions`, add an `if (IsPresent(...))` block that sets `Correlation` and increments `NFoundCorrelations`.
3. Add a `case` to `TransLMCorrelations::ReThetaC_Correlations` and `TransLMCorrelations::FLength_Correlations`.
4. Add the name to the `switch (lmParsedOptions.Correlation)` in `CConfig::SetOutput`.
5. Add the keyword to the comment of `LM_OPTIONS` in `config_template.cfg`.

### A variant of the LM model

1. Add the keyword to `LM_OPTIONS` and `LM_Options_Map`. Add a field to `LM_ParsedOptions`, and set it in `ParseLMOptions` (as for `LM2015`).
2. Add the checks in `CConfig::SetPostprocessing` (as the 3D check of `LM2015`).
3. Use the field in `CSourcePieceWise_TransLM::ComputeResidual`, which keeps a copy in `options`. If needed, also use it in `CTransLMSolver`, which also has `options`. If the variant changes the number of equations, these parts must follow: `nVar` in the solver, the variable class, the edge kernel dispatch (`DispatchScheme<CScalarFlux_TransLM, 2>`), the boundary conditions, the restart and the output.
4. If the variant changes how the turbulence model uses the intermittency, change the coupling code listed in [Coupling with the turbulence model](#coupling-with-the-turbulence-model).
5. Print it in `CConfig::SetOutput`.

### A new transition model

1. Add a value to `TURB_TRANS_MODEL` and `Trans_Model_Map`. Every `switch (config->GetKind_Trans_Model())` in `CFlowOutput` must handle it.
2. Write four classes, as for a turbulence model:
   - the variable class, derived from `CTurbVariable`;
   - the solver, derived from `CTurbSolver`;
   - the source class, in `SU2_CFD/include/numerics/turbulent/transition/`;
   - the edge kernel.
3. Create the solver in `CSolverFactory::CreateTransSolver`, and the numerics in `CDriver::InstantiateTransitionNumerics`.
4. Many places test `TURB_TRANS_MODEL::LM` directly: `CDriver::InitializeNumerics` (the `transition` flag), `CFluidIteration::CommonAuxiliarySolvers` and `CFluidIteration::Update`, the turbulence solvers and sources, and `CTurbSSTVariable::SetBlendingFunc`. Find them with `git grep -n "TURB_TRANS_MODEL::LM"`.
5. Add the new `.cpp` files to `SU2_CFD/src/meson.build`.

## Checklist

1. Keyword in `LM_OPTIONS` (and `TURB_TRANS_CORRELATION`) or in `TURB_TRANS_MODEL`, with its map entry.
2. Parsing in `ParseLMOptions`, checks in `CConfig::SetPostprocessing`, summary in `CConfig::SetOutput`, keyword in `config_template.cfg`.
3. Source terms in `CSourcePieceWise_TransLM`, correlations in `TransLMCorrelations`.
4. Coupling terms in the SA and SST sources and in `CTurbSSTVariable::SetBlendingFunc` reviewed.
5. Output fields, and the order of the `SOLUTION` fields for the restart.
6. A regression test with the new keyword.

Last checked against SU2 develop `11c177b904`.
