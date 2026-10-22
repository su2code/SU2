---
title: Adding a Turbulence Model
permalink: /docs_v7/Adding-a-Turbulence-Model/
---

This page is for developers. It shows how the RANS (Reynolds-averaged Navier-Stokes) turbulence models of SU2 are built. It then explains what to change in two cases: a variant of an existing model, and a new model with its own equations. The code examples are copied from the SU2 source.

Read [Code Structure](../Code-Structure/) first. The parts that are not specific to turbulence have their own pages: [Adding a Config Option](../Adding-a-Config-Option/), [Adding Output Fields](../Adding-Output-Fields/), [Adding a Solver](../Adding-a-Solver/), [Adding a Boundary Condition](../Adding-a-Boundary-Condition/) and [Running Regression Tests](../Running-Regression-Tests/). Transition models are covered in [Adding a Transition Model](../Adding-a-Transition-Model/).

## Files you will touch

| File | What goes there |
|---|---|
| `Common/include/option_structure.hpp` | `TURB_MODEL` and `Turb_Model_Map`, `TURB_FAMILY` and `TurbModelFamily`, the model options (`SA_OPTIONS`, `SST_OPTIONS`), their maps, parsed structs and parse functions |
| `Common/include/CConfig.hpp`, `Common/src/CConfig.cpp` | The options, the parsing in `CConfig::SetPostprocessing`, the number of turbulence variables, the start-up summary in `CConfig::SetOutput` |
| `SU2_CFD/include/numerics/turbulent/turb_sources.hpp` | Source terms and their Jacobians |
| `SU2_CFD/include/numerics/turbulent/turb_sa_edge_flux.hpp`, `turb_sst_edge_flux.hpp` | Convection and diffusion, in the edge kernels. An edge kernel is a class that computes the fluxes of one mesh edge. |
| `SU2_CFD/include/numerics/scalar/scalar_edge_flux.hpp` | The shared base of the edge kernels. You use it, but you rarely change it. |
| `SU2_CFD/include/solvers/CTurb*Solver.hpp`, `SU2_CFD/src/solvers/CTurb*Solver.cpp` | The solver: residuals, boundary conditions, eddy viscosity |
| `SU2_CFD/include/variables/CTurb*Variable.hpp`, `SU2_CFD/src/variables/CTurb*Variable.cpp` | The data stored at each point |
| `SU2_CFD/src/solvers/CSolverFactory.cpp` | Creates the turbulence solver (`CSolverFactory::CreateTurbSolver`) |
| `SU2_CFD/src/drivers/CDriver.cpp` | Creates the source numerics (`CDriver::InstantiateTurbulentNumerics`) |
| `SU2_CFD/src/output/CFlowOutput.cpp`, `SU2_CFD/src/output/CAdjFlowOutput.cpp` | Screen, history and volume fields |
| `SU2_CFD/src/meson.build` | New source files |
| `config_template.cfg`, `TestCases/` | Option documentation and tests |

---

- [How the turbulence equations are solved](#how-the-turbulence-equations-are-solved)
- [Class hierarchy](#class-hierarchy)
- [Case A: a variant of an existing model](#case-a-a-variant-of-an-existing-model)
- [Case B: a new model with its own equations](#case-b-a-new-model-with-its-own-equations)
- [Coupling with the flow solver](#coupling-with-the-flow-solver)
- [Walls, wall functions and wall distance](#walls-wall-functions-and-wall-distance)
- [Free-stream, initial and inlet values](#free-stream-initial-and-inlet-values)
- [Output and restart](#output-and-restart)
- [Discrete adjoint](#discrete-adjoint)
- [Hybrid RANS/LES](#hybrid-ransles)
- [Tests](#tests)
- [Checklist](#checklist)

---

## How the turbulence equations are solved

The turbulence model has its own solver object. It sits at position `TURB_SOL` of the solver container, on every multigrid level. The coupling with the flow is segregated: the flow and the turbulence model are solved one after the other, not as one system. In each inner iteration, `CFluidIteration::Iterate` (`SU2_CFD/src/iteration/CFluidIteration.cpp`) first runs the flow (`MultiGrid_Iteration`). Then `CFluidIteration::CommonAuxiliarySolvers` runs the transition model and the turbulence model:

```cpp
  if (config[val_iZone]->GetKind_Turb_Model() != TURB_MODEL::NONE && !frozen_visc) {

    /*--- Solve transition model ---*/

    if (config[val_iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_TRANS_SYS);
      integration[val_iZone][val_iInst][TRANS_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                         RUNTIME_TRANS_SYS, val_iZone, val_iInst);
    }

    /*--- Solve the turbulence model ---*/

    config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_TURB_SYS);
    integration[val_iZone][val_iInst][TURB_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                      RUNTIME_TURB_SYS, val_iZone, val_iInst);
  }
```

`SetGlobalParam(..., RUNTIME_TURB_SYS)` selects the turbulence settings: `CONV_NUM_METHOD_TURB`, `MUSCL_TURB`, `SLOPE_LIMITER_TURB` and `TIME_DISCRE_TURB`. Then `CSingleGridIntegration::SingleGrid_Iteration` (`SU2_CFD/src/integration/CSingleGridIntegration.cpp`) runs these steps on the finest grid:

1. `Preprocessing`: gradients and limiters of the turbulence variables (`CScalarSolver::CommonPreprocessing`).
2. `SetTime_Step`: the local time step. It is the flow time step, scaled by the ratio of the turbulence and flow CFL numbers.
3. `CIntegration::Space_Integration`: `Upwind_Residual` (convection **and** diffusion, see below), `Viscous_Residual` (empty for the turbulence solvers), `Source_Residual`, then the boundary conditions.
4. `CIntegration::Time_Integration`: implicit or explicit update.
5. `Postprocessing`: gradients and the eddy viscosity.
6. For `RUNTIME_TURB_SYS` only: restriction of the solution and of the eddy viscosity to the coarse grids (`SetRestricted_Solution`, `SetRestricted_EddyVisc`). Restriction transfers values from the fine grid to the coarse grids. The flow multigrid uses these values.

The flow solver reads the new eddy viscosity in its next `Preprocessing` (see [Coupling with the flow solver](#coupling-with-the-flow-solver)).

## Class hierarchy

Each model is a set of classes: a solver, a variable class, a source term and an edge kernel. The table lists the base classes and the SA and SST classes for each role.

| Role | Base classes | SA | SST |
|---|---|---|---|
| Solver | [`CSolver`](../../doxygen/classCSolver.html) → [`CScalarSolver`](../../doxygen/classCScalarSolver.html)`<CTurbVariable>` → [`CTurbSolver`](../../doxygen/classCTurbSolver.html) | [`CTurbSASolver`](../../doxygen/classCTurbSASolver.html) | [`CTurbSSTSolver`](../../doxygen/classCTurbSSTSolver.html) |
| Variables | [`CVariable`](../../doxygen/classCVariable.html) → [`CScalarVariable`](../../doxygen/classCScalarVariable.html) → [`CTurbVariable`](../../doxygen/classCTurbVariable.html) | [`CTurbSAVariable`](../../doxygen/classCTurbSAVariable.html) | [`CTurbSSTVariable`](../../doxygen/classCTurbSSTVariable.html) |
| Source term | [`CNumerics`](../../doxygen/classCNumerics.html) | [`CSourceBase_TurbSA`](../../doxygen/classCSourceBase__TurbSA.html), made by `SAFactory` | [`CSourcePieceWise_TurbSST`](../../doxygen/classCSourcePieceWise__TurbSST.html) |
| Convection and diffusion | [`CAvgGradScalarBase`](../../doxygen/classCAvgGradScalarBase.html) → [`CUpwScalarFlux`](../../doxygen/classCUpwScalarFlux.html) → [`CUpwScalarBase`](../../doxygen/classCUpwScalarBase.html) | [`CScalarFlux_SA`](../../doxygen/classCScalarFlux__SA.html) | [`CScalarFlux_SST`](../../doxygen/classCScalarFlux__SST.html) |

The transition solver [`CTransLMSolver`](../../doxygen/classCTransLMSolver.html) also derives from `CTurbSolver`.

The base classes already do most of the work:

- `CScalarSolver` (`SU2_CFD/include/solvers/CScalarSolver.hpp` and `.inl`) provides:
  - gradients and limiters;
  - the generic edge loop `EdgeFluxResidual` and the generic boundary pass `BoundaryFluxResidual`;
  - the local time step, the implicit and explicit updates, and dual time stepping;
  - periodic boundaries;
  - the clipping of the solution between `lowerlimit` and `upperlimit`, in `CompleteImplicitIteration`.

  Its constructor argument `conservative` says if the transported quantity is the solution times the density. It is true for SST and false for SA.
- `CTurbSolver` (`SU2_CFD/src/solvers/CTurbSolver.cpp`) provides `LoadRestart`, `Impose_Fixed_Values` (option `TURB_FIXED_VALUES`), `RegisterSolutionExtra`, `BC_Giles`, `BC_TurboRiemann`, the inlet storage `Inlet_TurbVars` and `ComputeUnderRelaxationFactorHelper`. `RegisterSolutionExtra` registers the eddy viscosity for the discrete adjoint.
- `CTurbVariable` (`SU2_CFD/include/variables/CTurbVariable.hpp`) holds the eddy viscosity `muT` (`GetmuT`, `SetmuT`), `turb_index` and `intermittency`. It sets `MAXNVAR = 4`, the size of the static arrays of the scalar solver.

## Case A: a variant of an existing model

A variant is a correction or another version of SA or SST. It does not need a new solver. It needs three things:

1. A keyword in `SA_OPTIONS` or `SST_OPTIONS`.
2. A field in the parsed options struct. This struct holds the result of reading the option list.
3. A change in the source term. Some variants change the diffusion, the eddy viscosity or the wall values instead.

### How the model options are parsed

`SA_OPTIONS` and `SST_OPTIONS` are enum list options. So is `LM_OPTIONS`, for transition. The user can give several keywords, for example `SA_OPTIONS= NEGATIVE, WITHFT2`. `Common/include/option_structure.hpp` defines four things for each option: an enum, a map from keywords to enum values, a struct with the parsed result, and a parse function:

```cpp
enum class SA_OPTIONS {
  NONE,     /*!< \brief No option / default. */
  NEG,      /*!< \brief Negative SA. */
  EDW,      /*!< \brief Edwards version. */
  FT2,      /*!< \brief Use FT2 term. */
  // ...
  EXP,      /*!< \brief Allow experimental combinations of options (according to TMR). */
};
static const MapType<std::string, SA_OPTIONS> SA_Options_Map = {
  MakePair("NONE", SA_OPTIONS::NONE)
  MakePair("NEGATIVE", SA_OPTIONS::NEG)
  MakePair("EDWARDS", SA_OPTIONS::EDW)
  MakePair("WITHFT2", SA_OPTIONS::FT2)
  // ...
};

struct SA_ParsedOptions {
  SA_OPTIONS version = SA_OPTIONS::NONE;  /*!< \brief SA base model. */
  bool ft2 = false;                       /*!< \brief Use ft2 term. */
  bool qcr2000 = false;                   /*!< \brief Use QCR-2000. */
  bool comp = false;                      /*!< \brief Use compressibility correction. */
  bool rot = false;                       /*!< \brief Use rotation correction. */
  bool bc = false;                        /*!< \brief BC transition. */
};
```

`ParseSAOptions` turns the list into the struct. It checks for conflicts with a small lambda:

```cpp
  auto IsPresent = [&](SA_OPTIONS option) {
    const auto sa_options_end = SA_Options + nSA_Options;
    return std::find(SA_Options, sa_options_end, option) != sa_options_end;
  };
```

Then it sets `version` (`NONE`, `NEG` or `EDW`) and the booleans. Some combinations are not standard versions of SA. The parser stops on them, unless `EXPERIMENTAL` is in the list. Its error message refers to the NASA Turbulence Modeling Resource.

`SST_ParsedOptions` and `ParseSSTOptions` work the same way. The struct has these fields:

- `version`: `V1994` or `V2003`.
- `production`: one production modifier, `V`, `KL`, `UQ`, `COMP_Wilcox` or `COMP_Sarkar`.
- the booleans `sust`, `uq`, `modified`, `compWilcox`, `compSarkar` and `dll`.

With no version keyword, the parser selects `V2003` with `modified = true`. The keywords `V1994` and `V2003` (the unmodified versions) are commented out in `SST_Options_Map`, so users cannot select them.

`CConfig::SetConfig_Options` registers the options:

```cpp
  addEnumListOption("SST_OPTIONS", nSST_Options, SST_Options, SST_Options_Map);
  // ...
  addEnumListOption("SA_OPTIONS", nSA_Options, SA_Options, SA_Options_Map);
```

`CConfig::SetPostprocessing` parses them, only for the selected model:

```cpp
  /*--- Postprocess SST_OPTIONS into structure. ---*/
  if (Kind_Turb_Model == TURB_MODEL::SST) {
    sstParsedOptions = ParseSSTOptions(SST_Options, nSST_Options, rank);
  } else if (Kind_Turb_Model == TURB_MODEL::SA) {
    saParsedOptions = ParseSAOptions(SA_Options, nSA_Options, rank);
  }
```

The rest of the code reads the result with `config->GetSAParsedOptions()` and `config->GetSSTParsedOptions()`.

### Steps for a new option

1. Add the value to `enum class SA_OPTIONS` (or `SST_OPTIONS`). Add a `MakePair` line to `SA_Options_Map` (or `SST_Options_Map`).
2. Add a field to `SA_ParsedOptions` (or `SST_ParsedOptions`). Set it in `ParseSAOptions` (or `ParseSSTOptions`) with `IsPresent`. Add the conflict checks. For SA, decide if the new option is allowed without `EXPERIMENTAL`.
3. Use the field where the model changes. This is usually the source term (see the next sections). Sometimes it is the diffusion, the eddy viscosity, the limits or the wall values.
4. Add the checks that involve other options in `CConfig::SetPostprocessing`. Example from the code:

   ```cpp
   if (Kind_Solver == MAIN_SOLVER::INC_RANS && sstParsedOptions.compSarkar){
     SU2_MPI::Error("COMPRESSIBILITY-SARKAR only supported for SOLVER= RANS", CURRENT_FUNCTION);
   }
   ```

5. Print the variant in the start-up summary in `CConfig::SetOutput`. This is the `switch (Kind_Turb_Model)` that prints `Spalart-Allmaras-noft2-R...` or `Menter's k-omega SST-2003m ...`.
6. Add the keyword to the comment of `SA_OPTIONS` or `SST_OPTIONS` in `config_template.cfg`.
7. Add a regression test (see [Tests](#tests)).

### Where an SA variant changes the source term

The SA source term is a class template in `SU2_CFD/include/numerics/turbulent/turb_sources.hpp`. Each piece of the model that changes between versions is a template parameter:

```cpp
template <class FlowIndices, class Omega, class ft2, class ModVort, class rFunc, class SourceTerms>
class CSourceBase_TurbSA : public CNumerics {
```

| Parameter | Computes | Implementations in `namespace detail` |
|---|---|---|
| `Omega` | the vorticity (or strain) magnitude `var.Omega` | `Omega::Bsl`, `Omega::Edw` |
| `ft2` | `var.ft2` and `var.d_ft2` | `ft2::Zero`, `ft2::Nonzero` |
| `ModVort` | the modified vorticity `var.Shat` and `var.d_Shat` | `ModVort::Bsl`, `ModVort::Edw`, `ModVort::Neg` |
| `rFunc` | `var.r` and `var.d_r` | `r::Bsl`, `r::Edw` |
| `SourceTerms` | production, destruction and the Jacobian | `SourceTerms::Bsl`, `SourceTerms::Neg` |

Each implementation is a struct with a static `get` function. It reads and writes the shared struct `CSAVariables`. This struct holds constants such as `cb1`, `cw2` and `cv1_3`. It also holds the auxiliary values and their derivatives with respect to the SA variable, named `d_...`. For example:

```cpp
/*! \brief Non-zero ft2 term according to the literature. */
struct Nonzero {
  static void get(CSAVariables& var) {
    const su2double xsi2 = pow(var.Ji, 2);
    var.ft2 = var.ct3 * exp(-var.ct4 * xsi2);
    var.d_ft2 = -2.0 * var.ct4 * var.Ji * var.ft2 * var.d_Ji;
  }
};
```

`CSourceBase_TurbSA::ComputeResidual` calls them in this order: `Omega::get`, `ModVort::get`, `ft2::get`, `rFunc::get`, the `fw` function, then `SourceTerms::get`. It returns the residual and the 1x1 Jacobian:

```cpp
      /*--- Compute production, destruction and jacobian ---*/
      su2double Production = 0.0, Destruction = 0.0;
      SourceTerms::get(ScalarVar_i[0], var, Production, Destruction, Jacobian_i[0][0]);
      // ...
      Residual[0] = (Production - Destruction) * Volume;
      // ...
      Jacobian_i[0][0] *= Volume;
```

Helper structs and `SAFactory` instantiate the template. The driver calls `SAFactory`:

```cpp
template <class FlowIndices>
struct BaselineSA {
  template <class Ft2>
  using type = CSourceBase_TurbSA<FlowIndices, Omega::Bsl, Ft2, ModVort::Bsl, r::Bsl, SourceTerms::Bsl>;
};
// ... NegativeSA, EdwardsSA

template <class FlowIndices>
CNumerics* SAFactory(unsigned short nDim, const CConfig* config) {
  const auto options = config->GetSAParsedOptions();

  switch (options.version) {
    case SA_OPTIONS::NONE:
      return detail::SAFactoryImpl<detail::BaselineSA<FlowIndices>>(options.ft2, options.comp, nDim, config);
    case SA_OPTIONS::NEG:
      return detail::SAFactoryImpl<detail::NegativeSA<FlowIndices>>(options.ft2, options.comp, nDim, config);
    case SA_OPTIONS::EDW:
      return detail::SAFactoryImpl<detail::EdwardsSA<FlowIndices>>(options.ft2, options.comp, nDim, config);
    default:
      break;
  }
  return nullptr;
}
```

`SAFactoryImpl` picks `ft2::Nonzero` or `ft2::Zero`. When `options.comp` is set, `AddCompressibilityCorrection` wraps the class in the decorator `CCompressibilityCorrection`. A decorator is a class that derives from the source class and adds to its result.

There are three ways to add an SA variant. Choose the one that matches the change:

1. **A different version of one piece**, for example another modified vorticity:
   - write a new struct with a static `get` in the matching group (`ModVort`, `r`, ...);
   - set the value **and** its derivative;
   - add a helper like `BaselineSA` that uses it;
   - add a `case` to `SAFactory`.
2. **An extra term added to the residual**, like SA-comp. Write a decorator like `CCompressibilityCorrection<ParentClass>`. Its `ComputeResidual` calls `ParentClass::ComputeResidual(config)`. It then changes `this->Residual[0]` and `this->Jacobian_i[0][0]`. Select the decorator in a factory helper like `AddCompressibilityCorrection`.
3. **A small switch inside the common code.** Test the option in `CSourceBase_TurbSA::ComputeResidual`, which keeps its own copy of the options (`options`). The rotation correction does this:

   ```cpp
      /*--- Dacles-Mariani et. al. rotation correction ("-R"). ---*/
      if (options.rot) {
        var.Prod += var.CRot * min(0.0, StrainMag_i - var.Omega);
        /*--- Do not allow negative production for SA-neg. ---*/
        if (ScalarVar_i[0] < 0) var.Prod = abs(var.Prod);
      }
   ```

Some SA options also act outside the source term:

- `NEGATIVE` changes three things:
  - the diffusion coefficient (`fn` in `CScalarFlux_SA`, `turb_sa_edge_flux.hpp`);
  - the eddy viscosity, which `CTurbSASolver::Postprocessing` sets to zero for a negative SA variable;
  - the lower limit. The `CTurbSASolver` constructor sets `lowerlimit[0] = EPS` only for the other versions.
- `QCR2000` changes the viscous stresses of the flow, in `SU2_CFD/src/numerics/flow/flow_diffusion.cpp` and `SU2_CFD/include/numerics_simd/flow/diffusion/viscous_fluxes.hpp`.
- `BCM` changes the source term (`options.bc` in `ComputeResidual`). It also changes the free-stream and inlet values, with the factor `0.005` in the `CTurbSASolver` constructor and in `CTurbSASolver::BC_Inlet`.

### Where an SST variant changes the source term

The SST source term is `CSourcePieceWise_TurbSST<FlowIndices>`, in the same file. It reads the parsed options from `sstParsedOptions`. This is a member of `CNumerics`, set by the `CNumerics` constructor with `config->GetSSTParsedOptions()`. The production modifiers are a `switch`:

```cpp
      /*--- Apply production term modifications ---*/
      switch (sstParsedOptions.production) {
        case SST_OPTIONS::UQ:
          ComputePerturbedRSM(nDim, Eig_Val_Comp, uq_permute, uq_delta_b, uq_urlx, PrimVar_Grad_i + idx.Velocity(),
                            Density_i, Eddy_Viscosity_i, ScalarVar_i[0], MeanPerturbedRSM);
          P_Base = PerturbedStrainMag(ScalarVar_i[0]);
          break;

        case SST_OPTIONS::V:
          P_Base = VorticityMag;
          break;
        // ... KL, COMP_Wilcox, COMP_Sarkar

        default:
          /*--- Base production term for SST-1994 and SST-2003 ---*/
          P_Base = StrainMag_i;
          break;
      }
```

`production` holds one value. `ParseSSTOptions` first sets it from `VORTICITY`, `KATO-LAUNDER` or `UQ`. Then it sets it from the compressibility keywords, which overwrite the first value.

The other SST options act in these places:

| Option field | Where |
|---|---|
| `version` | The constants `constants[8]`, `constants[9]`, `constants[10]` in the `CTurbSSTSolver` constructor. The `pw` term in `CSourcePieceWise_TurbSST::ComputeResidual`. The eddy viscosity, which uses the vorticity for `V1994` and the strain rate otherwise, both in the source term (`zeta`) and in `CTurbSSTSolver::Postprocessing`. |
| `sust` | The sustaining terms (`kAmb`, `omegaAmb`) in `ComputeResidual` |
| `dll` | The lower limits `lowerlimit[0]`, `lowerlimit[1]` in the `CTurbSSTSolver` constructor (options `LOWER_LIMIT_K_FACTOR`, `LOWER_LIMIT_OMEGA_FACTOR`) |
| `uq` | The Reynolds stress of the flow viscous flux (`viscous_fluxes.hpp`) |

The source Jacobian of SST has only the dissipation terms: `Jacobian_i[0][0]`, `Jacobian_i[0][1]` and `Jacobian_i[1][1]`. The production is treated explicitly. `CTurbSSTSolver` is conservative, so the Jacobian is with respect to the conservative variables: density times k, and density times omega.

To add an SST production modifier:

1. Add the keyword.
2. Set `production` in `ParseSSTOptions`, and update the error message that lists the modifiers.
3. Add a `case` to the `switch`.
4. Print it in `CConfig::SetOutput`.

## Case B: a new model with its own equations

A new model needs its own solver, variable, source and edge kernel classes. Start from a copy of `CTurbSSTSolver` (2 equations, conservative) or `CTurbSASolver` (1 equation, not conservative). For the general steps of a new solver, see [Adding a Solver](../Adding-a-Solver/).

### 1. Model keyword and family

In `Common/include/option_structure.hpp`, add a value to `TURB_MODEL` and to `Turb_Model_Map`. Then map it to a family in `TurbModelFamily`:

```cpp
inline TURB_FAMILY TurbModelFamily(TURB_MODEL model) {
  switch (model) {
    case TURB_MODEL::NONE:
      return TURB_FAMILY::NONE;
    case TURB_MODEL::SA:
      return TURB_FAMILY::SA;
    case TURB_MODEL::SST:
      return TURB_FAMILY::KW;
  }
  return TURB_FAMILY::NONE;
}
```

`TURB_FAMILY` has three values: `NONE`, `SA` (SA variants) and `KW` (k-omega models). Most of the code switches on the family, not on the model. If the new model does not fit an existing family, add a value to `TURB_FAMILY`. The compiler then warns at each `switch (TurbModelFamily(...))` that does not handle it. These switches are in:

- `CConfig::SetPostprocessing`: the number of turbulence variables `nTurbVar`;
- `CSolverFactory::CreateTurbSolver`: which solver class is created;
- `CFlowOutput` and `CAdjFlowOutput`: the output fields;
- `CSolver::LoadInletProfile`: the column names of the inlet profile template file.

Some code tests `TURB_MODEL::SST` or `TURB_MODEL::SA` directly. Check each place with `git grep -n "TURB_MODEL::SST"` and `git grep -n "TURB_MODEL::SA\b"`. The main ones are:

- `tkeNeeded`. The flow solvers add the turbulent kinetic energy to the energy and to the stresses only for `TURB_MODEL::SST`. This happens in `CNSSolver::SetPrimitive_Variables`, `CIncNSSolver`, `CFVMFlowSolverBase::Viscous_Residual_impl`, and many boundary conditions in `CEulerSolver.cpp` and `CIncEulerSolver.cpp`. A new model that transports k must be added there.
- `CConfig::GetInlet_TurbVal`, and the check of the number of values of `MARKER_INLET_TURBULENT` in `CConfig::SetPostprocessing`.
- `CDriver::InstantiateTurbulentNumerics`.
- The turbomachinery averaging in `CEulerSolver::PreprocessAverage` and `CEulerSolver::TurboAverageProcess`.

### 2. Configuration

Add the checks of the new model in `CConfig::SetPostprocessing`, and its name in `CConfig::SetOutput`. If users can change options or constants of the model, add them as described in [Adding a Config Option](../Adding-a-Config-Option/). For variants, follow the `SA_OPTIONS` pattern above.

### 3. Variable class

Derive from `CTurbVariable`. The constructor sets the initial solution and allocates `muT`. `CTurbVariable` does not allocate it. From `CTurbSSTVariable`:

```cpp
CTurbSSTVariable::CTurbSSTVariable(su2double kine, su2double omega, su2double mut, unsigned long npoint, unsigned long ndim, unsigned long nvar, const su2double* constants, CConfig *config)
  : CTurbVariable(npoint, ndim, nvar, config) {
  // ...
  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
  {
    Solution(iPoint,0) = kine;
    Solution(iPoint,1) = omega;
  }

  Solution_Old = Solution;
  // ...
  F1.resize(nPoint) = su2double(1.0);
  F2.resize(nPoint) = su2double(0.0);
  CDkw.resize(nPoint) = su2double(0.0);

  muT.resize(nPoint) = mut;
}
```

Add the per-point fields that the model needs. For example, SST stores the blending functions `F1` and `F2` and the cross diffusion `CDkw`, computed in `CTurbSSTVariable::SetBlendingFunc`. Expose them with virtual getters of `CVariable` (for example `GetF1blending`, `GetCrossDiff`), or add new ones. The model can have at most `CTurbVariable::MAXNVAR` (4) equations, unless you increase it.

### 4. Solver class

Derive from `CTurbSolver`. Its constructor is:

```cpp
  CTurbSolver(CGeometry* geometry, CConfig *config, const CSolver* flow_solver, bool conservative);
```

The constructor of the new solver does what `CTurbSSTSolver::CTurbSSTSolver` does:

1. Set `nVar`, `nPrimVar`, `nVarGrad`, `nPoint`, `nPointDomain`.
2. On `MESH_0`, or with a full multigrid cycle, allocate the residual vectors and the linear system: `Jacobian.Initialize(...)`, `LinSysSol.Initialize(...)`, `LinSysRes.Initialize(...)`. Also allocate `EdgeFluxes` and `EdgeFluxesDiff` when `ReducerStrategy` is set. With this strategy, used when edge coloring is not efficient, the fluxes are stored per edge and summed per point afterwards.
3. Set the model constants, `Solution_Inf`, `lowerlimit` and `upperlimit`.
4. Create `nodes` and call `SetBaseClassPointerToNodes()`. Create `ghostNodes`, sized to the largest marker. The boundary conditions use `ghostNodes` as the outside state.
5. Exchange the halo values with `InitiateComms` / `CompleteComms` and `MPI_QUANTITIES::SOLUTION_EDDY`. Halo points are local copies of points owned by other MPI ranks.
6. Allocate `SlidingState` for `FLUID_INTERFACE` markers, and `Inlet_TurbVars` for all markers.
7. Set the local CFL (`config->GetCFL(MGLevel)*config->GetCFLRedCoeff_Turb()`) and `SolverName`.

Then implement these functions:

| Function | What it does in SST |
|---|---|
| `Preprocessing` | `CommonPreprocessing(geometry, config, Output)` |
| `Postprocessing` | gradients, blending functions, eddy viscosity (see [Coupling with the flow solver](#coupling-with-the-flow-solver)) |
| `Upwind_Residual` | the edge loop with the model kernel (step 5) |
| `Source_Residual` | the source term at each point (step 6) |
| `BC_HeatFlux_Wall`, `BC_Isothermal_Wall` | wall values (see [Walls](#walls-wall-functions-and-wall-distance)) |
| `BC_Inlet`, `BC_Outlet`, `BC_Far_Field`, `BC_Inlet_MixingPlane`, `BC_Inlet_Turbo`, `BC_Fluid_Interface` | ghost state, then the shared boundary flux |
| `SetInletAtVertex`, `GetInletAtVertex`, `SetUniformInlet` | inlet profile files |
| `ComputeUnderRelaxationFactor` | calls `ComputeUnderRelaxationFactorHelper` with `MAX_UPDATE_SST` (SA uses `MAX_UPDATE_SA`) |
| `GetConstants`, `GetTke_Inf`, `GetOmega_Inf` | used by the driver to build the source numerics |

The solver `.cpp` file includes `SU2_CFD/include/solvers/CScalarSolver.inl`. This is needed because `EdgeFluxResidual` and `BoundaryFluxResidual` are member templates, instantiated with the kernel of the model.

### 5. Convection and diffusion: the edge kernel

Each model has an edge kernel. It computes the convection and the diffusion of one edge, with their Jacobians, in one pass. The turbulence solvers do not use `CNumerics` objects for these terms. The interior loop is one call:

```cpp
void CTurbSSTSolver::Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics**,
                                     CConfig* config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  const auto opt = ScalarFluxOptions::Interior(*config, config->GetBounded_Turb(),
                                               config->GetUse_Accurate_Turb_Jacobians());

  DispatchScheme<CScalarFlux_SST, 2>(config, [&](auto tag) {
    EdgeFluxResidual<typename decltype(tag)::type>(geometry, solver_container, config, opt);
  });
}
```

`DispatchScheme` turns run-time settings into a compile-time kernel type. It picks the flow indices (compressible or incompressible), the dimension and the number of equations at run time. It then calls the lambda with the matching kernel type. The template arguments after the kernel are the allowed numbers of equations. SA uses `DispatchScheme<CScalarFlux_SA, 1, 4>`. `ScalarFluxOptions` has named constructors for the interior (`Interior`) and for the boundaries (`BoundaryConvective`, `BoundaryFull`, `BoundaryDiffusive`).

A kernel derives from `CUpwScalarBase<Double, Derived, FlowIndices, nDim, nVar>` (`SU2_CFD/include/numerics/scalar/scalar_edge_flux.hpp`). It provides:

- `static constexpr bool Conservative`: the convected quantity is density times the solution;
- `static constexpr bool DiagonalDiffusion`: `coefficients` returns a vector (one coefficient per equation) instead of a matrix;
- `coefficients(...)`: the diffusion coefficients for both orientations of the edge.

It can also provide:

- `coefficientJacobians(...)`: extra Jacobian terms when the coefficients depend on the solution, used with `USE_ACCURATE_TURB_JACOBIANS`;
- `finalizeFlux(...)`: a convective flux different from the default upwind one, as SA does;
- `extraDiffusionTerms(...)`: diffusion terms beyond the standard one, as the flamelet kernel does;
- `DiffusionReadsDensity`: set it to true when the diffusion coefficients need the density, as SA does. The base class then reads the density.

The SST kernel declares:

```cpp
template <class Double, class FlowIndices, int nDim, size_t nVar = 2>
class CScalarFlux_SST
    : public CUpwScalarBase<Double, CScalarFlux_SST<Double, FlowIndices, nDim, nVar>, FlowIndices, nDim, nVar> {
 public:
  static constexpr bool Conservative = true;
  static constexpr bool DiagonalDiffusion = false;
```

Its `coefficients` returns a 2x2 matrix per side. The matrix holds the blended diffusion of k and omega and the cross-diffusion term. The coefficients read, through `gatherVariables`:

- the flow primitives (`GetPrimitive`);
- the turbulence containers (`GetSolution`, `GetF1blending`).

So any variable container that a kernel reads must be stored at every point, like `CTurbVariable::GetF1blending()`. `CTransLMSolver` has the smallest kernel, `CScalarFlux_TransLM`. See [Adding a Transition Model](../Adding-a-Transition-Model/).

### 6. Source term

Write a class in `SU2_CFD/include/numerics/turbulent/turb_sources.hpp`. It derives from `CNumerics` and is a template on `FlowIndices`. `FlowIndices` is the index helper of the flow primitives: `idx.Density()`, `idx.Velocity()`, `idx.LaminarViscosity()`, ... Follow `CSourcePieceWise_TurbSST`:

- Call `CNumerics(nDim, nVar, config)` in the constructor. Keep the residual and the Jacobian in fixed arrays: `Residual[2]`, and `Jacobian_i[2]` pointing into `Jacobian_Buffer[4]`.
- In `ComputeResidual(const CConfig* config)`, fill `Residual` and `Jacobian_i`, and return `ResidualType<>(Residual, Jacobian_i, nullptr)`. `Residual` is the source integrated over the volume, positive for production. `Jacobian_i` is its derivative with respect to the conservative variables of the model.
- Register every input with `AD::SetPreaccIn`, between `AD::StartPreacc()` and `AD::SetPreaccOut(Residual, nVar); AD::EndPreacc();` (see [Discrete adjoint](#discrete-adjoint)).

The solver sets the inputs point by point in `Source_Residual`. It adds the result with a minus sign, because the residual is "fluxes minus sources":

```cpp
    /*--- Compute the source term ---*/

    auto residual = numerics->ComputeResidual(config);
    // ...
    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, residual);
    if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
```

The inputs are set with the setters of `CNumerics`: `SetPrimitive`, `SetPrimVarGradient`, `SetScalarVar`, `SetScalarVarGradient`, `SetVolume`, `SetDistance`, `SetVorticity`, `SetStrainMag`, `SetCoord` and `SetIntermittencyEff`. The SST setters are virtual: `SetF1blending`, `SetF2blending`, `SetCrossDiff`. If the model needs a new input, add a virtual setter to `CNumerics` (`SU2_CFD/include/numerics/CNumerics.hpp`).

### 7. Create the source numerics

`CDriver::InstantiateTurbulentNumerics` (`SU2_CFD/src/drivers/CDriver.cpp`) creates one source object per multigrid level. It also creates one per OpenMP thread, through `offset`. Add the model to the `switch` on `GetKind_Turb_Model()` at the top, and to the loop:

```cpp
  for (auto iMGlevel = 0u; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    auto& turb_source_first_term = numerics[iMGlevel][TURB_SOL][source_first_term];

    if (spalart_allmaras)
      turb_source_first_term = SAFactory<Indices>(nDim, config);
    else if (menter_sst)
      turb_source_first_term = new CSourcePieceWise_TurbSST<Indices>(nDim, nVar_Turb, constants, kine_Inf, omega_Inf,
                                                                     config);

    numerics[iMGlevel][TURB_SOL][source_second_term] = new CSourceNothing(nDim, nVar_Turb, config);
  }
```

The function is a template on the flow indices. At the end of its definition, it is instantiated for `CEulerVariable::CIndices<unsigned short>` and `CIncEulerVariable::CIndices<unsigned short>`. It also checks that `CONV_NUM_METHOD_TURB` is an upwind scheme. It creates no convective or viscous numerics, because the edge kernels replace them.

### 8. Create the solver

`CSolverFactory::CreateTurbSolver` (`SU2_CFD/src/solvers/CSolverFactory.cpp`) switches on the family:

```cpp
      case TURB_FAMILY::KW:
        turbSolver = new CTurbSSTSolver(geometry, config, solver[FLOW_SOL], iMGLevel);
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
        turbSolver->Postprocessing(geometry, solver, config, iMGLevel);
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
        break;
```

The calls after the constructor compute the initial eddy viscosity and pass it to the flow. When `adjoint` is true, the same function also creates the adjoint solver (`CDiscAdjSolver`). This part needs no change.

### 9. Build

Add the new `.cpp` files to `SU2_CFD/src/meson.build`, next to `solvers/CTurbSSTSolver.cpp` and `variables/CTurbSSTVariable.cpp`.

## Coupling with the flow solver

The flow reads two things from the turbulence solver: the eddy viscosity and, for SST, the turbulent kinetic energy. `CNSSolver::SetPrimitive_Variables` (`SU2_CFD/src/solvers/CNSSolver.cpp`) copies them into the primitive variables of each point:

```cpp
    if (turb_model != TURB_MODEL::NONE && solver_container[TURB_SOL] != nullptr) {
      eddy_visc = solver_container[TURB_SOL]->GetNodes()->GetmuT(iPoint);
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);
      // ...
    }
    // ...
    bool physical = ns_nodes->SetPrimVar(iPoint, eddy_visc, turb_ke, GetFluidModel());
```

`CIncNSSolver::SetPrimitive_Variables` does the same. `CNEMONSSolver` reads the eddy viscosity in the same way. `CFVMFlowSolverBase::Viscous_Residual_impl` passes k to the viscous numerics with `SetTurbKineticEnergy`.

The turbulence solver computes the eddy viscosity in `Postprocessing` and stores it with `nodes->SetmuT`. For SST:

```cpp
    const auto& eddy_visc_var = sstParsedOptions.version == SST_OPTIONS::V1994 ? VorticityMag : StrainMag;
    const su2double muT = max(0.0, rho * a1 * kine / max(a1 * omega, eddy_visc_var * F2));

    nodes->SetmuT(iPoint, muT);
```

SA (`CTurbSASolver::Postprocessing`) computes `muT = rho*fv1*nu_hat`, with the roughness correction of `Ji`. With `SA_OPTIONS= NEGATIVE`, it sets `muT = 0` for a negative SA variable.

In the other direction, the turbulence solver reads the flow through `solver_container[FLOW_SOL]`:

- density, laminar viscosity and velocity gradients;
- vorticity and strain rate magnitude (`GetVorticity`, `GetStrainMag`);
- in the edge kernels, the velocity and the eddy viscosity from the primitive variables.

## Walls, wall functions and wall distance

**Wall distance.** `CGeometry::ComputeWallDistance` (`Common/src/geometry/CGeometry.cpp`) computes the distance to the nearest viscous wall. It does this for the `RANS` and `INC_RANS` solvers and their adjoints, and for `FEM_RANS` and `FEM_LES`. The drivers call it at start-up. `CSinglezoneDriver` calls it again when the mesh moves or deforms. Read it with `geometry->nodes->GetWall_Distance(iPoint)`. SA adds `0.03` times the roughness height (`GetRoughnessHeight`) to it.

**Wall values.** `BC_HeatFlux_Wall` imposes the wall values directly, as Dirichlet conditions. `BC_Isothermal_Wall` calls it. The pattern has three steps:

1. Set the value in the solution.
2. Zero the residual row.
3. Replace the Jacobian row by the identity.

From `CTurbSSTSolver::BC_HeatFlux_Wall` (smooth wall):

```cpp
      solution[0] = 0.0;
      solution[1] = 60.0*laminar_viscosity/(density*beta_1*pow(wall_dist,2));
    }

    /*--- Set the solution values and zero the residual ---*/
    nodes->SetSolution_Old(iPoint, solution);
    nodes->SetSolution(iPoint, solution);
    LinSysRes.SetBlock_Zero(iPoint);

    if (implicit) {
      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      Jacobian.DeleteValsRowi(iPoint, 0);
      Jacobian.DeleteValsRowi(iPoint, 1);
    }
```

Rough walls are read with `config->GetWallRoughnessProperties(Marker_Tag)` (option `WALL_ROUGHNESS`). SST chooses the rough wall values with `KIND_ROUGHSST_MODEL` (`ROUGHSST_MODEL` in `option_structure.hpp`). On rough walls, SA imposes a flux instead of a value. The wall-normal derivative of the SA variable equals the SA variable divided by 0.03 times the roughness height.

**Wall functions.** With `MARKER_WALL_FUNCTIONS`, the flow solver computes the wall shear stress and the eddy viscosity at the wall (`CNSSolver::SetTau_Wall_WF`). `BC_HeatFlux_Wall` then calls `SetTurbVars_WF`. This function fixes the turbulence variables at the first point off the wall (`GetNormal_Neighbor()`). It uses `GetYPlus`, `GetUTau` and `GetEddyViscWall`, and acts only when y+ is above `WALLMODEL_MINYPLUS`. A new model needs its own `SetTurbVars_WF`.

The other boundary conditions (inlet, outlet, far field, ...) write the outside state into `ghostNodes` and call the shared boundary flux. See [Adding a Boundary Condition](../Adding-a-Boundary-Condition/).

## Free-stream, initial and inlet values

Three options set the free-stream turbulence. They are registered in `CConfig::SetConfig_Options`:

| Option | Getter | Used by |
|---|---|---|
| `FREESTREAM_TURBULENCEINTENSITY` | `GetTurbulenceIntensity_FreeStream()` | SST (k), LM, SA-BCM |
| `FREESTREAM_TURB2LAMVISCRATIO` | `GetTurb2LamViscRatio_FreeStream()` | SST (omega) |
| `FREESTREAM_NU_FACTOR` | `GetNuFactor_FreeStream()` | SA |

Each turbulence solver computes its free-stream state `Solution_Inf` in its constructor. It uses the non-dimensional free-stream flow: `GetDensity_FreeStreamND()`, `GetVelocity_FreeStreamND()`, `GetViscosity_FreeStreamND()`. It then initializes all points with this state. For SST:

```cpp
  su2double kine_Inf  = 3.0/2.0*(VelMag2*Intensity*Intensity);
  su2double omega_Inf = rhoInf*kine_Inf/(muLamInf*viscRatio);

  Solution_Inf[0] = kine_Inf;
  Solution_Inf[1] = omega_Inf;
```

`CEulerSolver::SetNondimensionalization` also computes the dimensional and non-dimensional free-stream k and omega. The incompressible and NEMO solvers do the same in their own `SetNondimensionalization`. The values are stored in `CConfig` (`SetTke_FreeStream`, `SetOmega_FreeStream`, ...). They are printed in the start-up table and used to write the inlet profile template.

`Solution_Inf` has two uses: it is the far-field state (`BC_Far_Field`), and it is the value imposed by `TURB_FIXED_VALUES`. `Inlet_TurbVars` holds the inlet values for each vertex. Where they come from depends on the inlet type:

- Without an inlet file, the values come from `MARKER_INLET_TURBULENT` or from the free-stream options, through `CConfig::GetInlet_TurbVal`.
- With an inlet file (`SPECIFIED_INLET_PROFILE= YES`), `SetInletAtVertex` reads the turbulence columns after the flow columns. `SetUniformInlet` fills the default values.

## Output and restart

The turbulence fields are declared and set in `CFlowOutput` (`SU2_CFD/src/output/CFlowOutput.cpp`). The functions that handle them switch on `TurbModelFamily(config->GetKind_Turb_Model())`:

- `AddHistoryOutputFields_ScalarRMS_RES`, `AddHistoryOutputFields_ScalarMAX_RES`, `AddHistoryOutputFields_ScalarBGS_RES`;
- `LoadHistoryDataScalar`;
- `SetVolumeOutputFieldsScalarSolution`, `SetVolumeOutputFieldsScalarResidual`, `SetVolumeOutputFieldsScalarLimiter`;
- `LoadVolumeDataScalar`.

The SST fields are `RMS_TKE`, `RMS_DISSIPATION`, `TKE`, `DISSIPATION`, `RES_TKE`, ... The SA fields are `RMS_NU_TILDE`, `NU_TILDE`, ... The fields `EDDY_VISCOSITY`, `TURB_DELTA_TIME` and `TURB_CFL` exist for every model. The adjoint fields are in `CAdjFlowOutput` (`RMS_ADJ_TKE`, `ADJ_TKE`, ...). To add fields, see [Adding Output Fields](../Adding-Output-Fields/).

**Restart.** The restart file is the `SOLUTION` group of the volume output. `CTurbSolver::LoadRestart` reads the turbulence variables right after the coordinates and the flow variables:

```cpp
    /*--- Skip flow variables ---*/

    unsigned short skipVars = nDim + solver[MESH_0][FLOW_SOL]->GetnVar();
```

So declare the `SOLUTION` fields of a new model in `SetVolumeOutputFieldsScalarSolution`, in the order of the solution vector. `CFlowCompOutput::SetVolumeOutputFields` calls it right after the flow `SOLUTION` fields.

After reading, `LoadRestart` restricts the solution to the coarse grids. When no species or transition model is active, it also computes the eddy viscosity: flow `Preprocessing`, then turbulence `Postprocessing`. Otherwise, the restart of those solvers does it.

## Discrete adjoint

A new model needs no adjoint equations. The discrete adjoint differentiates the primal code with CoDiPack, an automatic differentiation (AD) tool. Check these points:

- `CSolverFactory::CreateTurbSolver` creates a `CDiscAdjSolver` for `TURB_SOL`, unless `FROZEN_VISC_DISC= YES`. `CDiscAdjFluidIteration` registers, records and extracts the turbulence solution together with the flow.
- `CTurbSolver::RegisterSolutionExtra` registers `muT` as an input and output of the recording. Keep it if the model computes the eddy viscosity in `Postprocessing`.
- The numerics use preaccumulation. Preaccumulation stores only the local derivatives of a code section's outputs with respect to its inputs on the tape (the record of operations that AD replays backwards). So every value that the result depends on must be declared as an input with `AD::SetPreaccIn`. A missing input gives wrong gradients without any error. For example, `CSourcePieceWise_TurbSST::ComputeResidual` registers the intermittency only when the LM model is on:

  ```cpp
    if (config->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      AD::SetPreaccIn(intermittency_eff_i);
      eff_intermittency = intermittency_eff_i;
    }
  ```

- Loops over points use `AD::StartNoSharedReading()` / `AD::EndNoSharedReading()`. These mark a section where no two OpenMP threads read the same value, which makes the reverse sweep faster.
- Add a test to `serial_regression_AD.py`. See `discadj_rans_naca0012_sa` and `discadj_rans_naca0012_sst`.

## Hybrid RANS/LES

Detached eddy simulation (DES) is implemented for SA only. `HYBRID_RANSLES` takes `SA_DES`, `SA_DDES`, `SA_ZDES` or `SA_EDDES` (`HybridRANSLES_Map`). `DES_CONST` sets the constant. The hooks are in `CTurbSASolver`:

1. `CTurbSASolver::Preprocessing` calls `SetDES_LengthScale`. This stores the length scale and the LES sensor in the variable (`SetDES_LengthScale`, `SetLES_Mode`).
2. `CTurbSASolver::Source_Residual` passes the DES length scale to the source term instead of the wall distance:

   ```cpp
      /*--- Set DES length scale ---*/

      numerics->SetDistance(nodes->GetDES_LengthScale(iPoint), 0.0);
   ```

3. `CNSSolver::SetPrimitive_Variables` copies the length scale and the LES sensor into the flow variables. The output reads them there (`DES_LENGTHSCALE`, `LES_SENSOR`).

The stochastic backscatter model (`STOCHASTIC_BACKSCATTER`) adds three Langevin (stochastic) equations to the SA solver. The solver then has 4 variables. This is why `CScalarFlux_SA` is dispatched with 1 or 4 equations. A new model with a DES variant needs the same hooks in its own solver.

## Tests

These regression tests cover the turbulence models. The scripts are in `TestCases/`.

| Model or feature | Test cases | Script |
|---|---|---|
| SA versions (baseline, negative, Edwards, compressibility, QCR) | `turbmod_sa_*_rae2822` in `TestCases/turbulence_models/sa/rae2822` | `parallel_regression.py` |
| SA, SST on airfoils and wings | `rae2822_sa`, `rae2822_sst`, `rae2822_sst_sust`, `turb_naca0012_sa`, `turb_naca0012_sst*`, `turb_oneram6` | `parallel_regression.py`, `serial_regression.py` |
| Flat plate, compressibility corrections, roughness | `turb_flatplate`, `turb_flatplate_CC_Wilcox`, `turb_flatplate_CC_Sarkar`, `turb_sst_flatplate_rough*` in `TestCases/rans/flatplate` | `parallel_regression.py`, `serial_regression.py` |
| Wall functions | `turb_sst_wallfunction_flatplate`, `turb_sa_wallfunction_flatplate`, `inc_turb_sst_wallfunction_flatplate`, ... in `TestCases/wallfunctions` | `serial_regression.py`, `parallel_regression.py` |
| SST-UQ | `turb_naca0012_1c`, ... in `TestCases/rans_uq` | `serial_regression.py`, `parallel_regression.py`, `hybrid_regression.py` |
| DDES, stochastic backscatter | `ddes_flatplate`, `sbs_backward_step` | `parallel_regression.py`, `hybrid_regression.py` |
| SA-BCM transition | `Schubauer_Klebanoff` in `TestCases/transition` | `serial_regression.py`, `parallel_regression.py` |
| Discrete adjoint | `discadj_rans_naca0012_sa`, `discadj_rans_naca0012_sst`, `discadj_incomp_turb_NACA0012_sa`, `discadj_incomp_turb_NACA0012_sst` | `serial_regression_AD.py`, `parallel_regression_AD.py` |

**Verification and validation.** `TestCases/vandv.py` runs cases from `TestCases/vandv/rans`:

- the flat plate (`flatplate_sst1994m`);
- the bump in channel (`bump_sst1994m`);
- a shock/boundary layer interaction (`swbli_sa`, `swbli_sst`);
- `dsma661_sa` and `dsma661_sst`;
- the 30P30N airfoil.

The script runs each case for a few iterations only, like a regression test. For a new model, run the flat plate and the bump on several grid levels. Compare the skin friction and the eddy viscosity profiles with published results before you open the pull request.

To add tests, follow [Running Regression Tests](../Running-Regression-Tests/#add-a-new-test). For a new model or option, add at least one entry to `parallel_regression.py`, for example a copy of a `turbulence_models/sa/rae2822` case with the new option. If the adjoint must work, also add one entry to `serial_regression_AD.py`. Use one configuration file per test.

## Checklist

1. Keyword: `TURB_MODEL` + `Turb_Model_Map` (+ `TurbModelFamily`, `TURB_FAMILY`), or `SA_OPTIONS`/`SST_OPTIONS` + map + parsed struct + parse function.
2. Checks in `CConfig::SetPostprocessing`, summary in `CConfig::SetOutput`, keyword in `config_template.cfg`.
3. Variable class derived from `CTurbVariable`, with `muT` allocated.
4. Solver class derived from `CTurbSolver`: constructor, `Preprocessing`, `Postprocessing` (eddy viscosity), `Upwind_Residual`, `Source_Residual`, wall and other boundary conditions, inlet functions.
5. Edge kernel derived from `CUpwScalarBase`. Source class derived from `CNumerics`, with the Jacobian and all preaccumulation inputs.
6. `CDriver::InstantiateTurbulentNumerics` and `CSolverFactory::CreateTurbSolver`.
7. Every `switch (TurbModelFamily(...))` and every test on `TURB_MODEL::SST` / `TURB_MODEL::SA` reviewed (`tkeNeeded`, inlet values, turbomachinery).
8. Output fields in `CFlowOutput` and `CAdjFlowOutput`. `SOLUTION` fields in the order the restart expects.
9. New files in `SU2_CFD/src/meson.build`.
10. Regression tests (primal and adjoint) and a verification case.

Last checked against SU2 develop `11c177b904`.
