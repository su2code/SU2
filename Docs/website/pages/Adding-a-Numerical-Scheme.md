---
title: Adding a Numerical Scheme
permalink: /docs_v7/Adding-a-Numerical-Scheme/
---

This page shows how to add a new convective scheme, upwind or centered, to the finite volume flow solvers of SU2. It also covers, more briefly, how to add a viscous flux or a source term. It follows the code of the compressible solver (`CEulerSolver`, `CNSSolver`) and points to the incompressible and NEMO versions.

For the overall layout of the code, see [Code Structure](../Code-Structure/). Turbulence, transition and species schemes are covered in other pages: [Adding a Turbulence Model](../Adding-a-Turbulence-Model/), [Adding a Transition Model](../Adding-a-Transition-Model/) and [Adding a Solver](../Adding-a-Solver/).

---

- [Files you will touch](#files-you-will-touch)
- [How the configuration picks the scheme](#how-the-configuration-picks-the-scheme)
- [Two code paths: classic and vectorized](#two-code-paths-classic-and-vectorized)
- [The CNumerics interface](#the-cnumerics-interface)
- [Add a classic upwind scheme](#add-a-classic-upwind-scheme)
- [MUSCL reconstruction and limiters](#muscl-reconstruction-and-limiters)
- [Add a vectorized scheme](#add-a-vectorized-scheme)
- [Centered schemes](#centered-schemes)
- [Jacobians](#jacobians)
- [Low-Mach and entropy fix options](#low-mach-and-entropy-fix-options)
- [Incompressible and NEMO schemes](#incompressible-and-nemo-schemes)
- [Viscous fluxes and source terms](#viscous-fluxes-and-source-terms)
- [Discrete adjoint](#discrete-adjoint)
- [Tests](#tests)
- [Checklist](#checklist)

---

## Files you will touch

| File | What goes there |
|---|---|
| `Common/include/option_structure.hpp` | A value in `enum class UPWIND` or `enum class CENTERED`, and its keyword in `Upwind_Map` or `Centered_Map` |
| `Common/src/CConfig.cpp` | Optional: checks in `CConfig::SetPostprocessing`, the screen summary in `CConfig::SetOutput` |
| `config_template.cfg` | The new keyword in the comment above `CONV_NUM_METHOD_FLOW` |
| `SU2_CFD/include/numerics/flow/convection/*.hpp`, `SU2_CFD/src/numerics/flow/convection/*.cpp` | The classic numerics class (one edge at a time) |
| `SU2_CFD/src/meson.build` | A new `.cpp` file in the `numerics/...` list |
| `SU2_CFD/src/drivers/CDriver.cpp` | The `new` of the numerics class in `CDriver::InitializeNumerics` |
| `SU2_CFD/include/numerics_simd/flow/convection/upwind.hpp` or `centered.hpp`, `SU2_CFD/include/numerics_simd/CNumericsSIMD.cpp` | Optional: the vectorized version and its factory entry |
| `SU2_CFD/src/solvers/CEulerSolver.cpp` | Only if the scheme needs new inputs, new sensors, or the vectorized path |
| `TestCases/parallel_regression.py`, `TestCases/serial_regression.py` (and the `_AD` lists) | A regression test |

## How the configuration picks the scheme

This section shows where to register the keyword of your scheme, and how the solver finds it.

The option `CONV_NUM_METHOD_FLOW` is registered in `CConfig::SetConfig_Options` (`Common/src/CConfig.cpp`) with `addConvectOption`:

```cpp
  addConvectOption("CONV_NUM_METHOD_FLOW", Kind_ConvNumScheme_Flow, Kind_Centered_Flow, Kind_Upwind_Flow);
```

The parser is `COptionConvect` in `Common/include/option_structure.inl`. It looks up the keyword in `Centered_Map` first, then in `Upwind_Map`. It sets three values:

- `Kind_ConvNumScheme_Flow`: `SPACE_CENTERED` or `SPACE_UPWIND` (`enum ENUM_SPACE`). Read it with `config->GetKind_ConvNumScheme_Flow()`.
- `Kind_Centered_Flow`: a `CENTERED` value. It is `CENTERED::NONE` for upwind schemes. Read it with `config->GetKind_Centered_Flow()`.
- `Kind_Upwind_Flow`: an `UPWIND` value. It is `UPWIND::NONE` for centered schemes. Read it with `config->GetKind_Upwind_Flow()`.

The same kind of option is used for `CONV_NUM_METHOD_ADJFLOW`, `CONV_NUM_METHOD_TURB`, `CONV_NUM_METHOD_ADJTURB`, `CONV_NUM_METHOD_SPECIES` and `CONV_NUM_METHOD_HEAT`.

The enums and maps are in `Common/include/option_structure.hpp`. For an upwind scheme, add a value to `enum class UPWIND` and a line to `Upwind_Map`:

```cpp
static const MapType<std::string, UPWIND> Upwind_Map = {
  MakePair("NONE", UPWIND::NONE)
  MakePair("ROE", UPWIND::ROE)
  MakePair("TURKEL_PREC", UPWIND::TURKEL)
  // ...
  MakePair("LAX-FRIEDRICH", UPWIND::LAX_FRIEDRICH)
  MakePair("UDS", UPWIND::UDS)
};
```

For a centered scheme, do the same with `enum class CENTERED` and `Centered_Map`.

Notes:

- Give a new scheme a keyword that is in one map only. `Centered_Map` is searched first, so `CONV_NUM_METHOD_FLOW` reads a keyword that is in both maps as centered. `LAX-FRIEDRICH` is such a keyword. `UPWIND::LAX_FRIEDRICH` is reached only through `RIEMANN_SOLVER_FEM`, an `addEnumOption` with `Upwind_Map` that the DG solver uses.
- The comment above the two enums mentions an `ENUM_CONVECTIVE`. It no longer exists. The two maps are all you need.
- Add the keyword to the list in `config_template.cfg` above `CONV_NUM_METHOD_FLOW`. For the general rules on enums and maps, see [Adding a Config Option](../Adding-a-Config-Option/#4-add-an-enum).

At run time, the solver picks the residual function as follows:

1. Before each equation system is solved, `CConfig::SetGlobalParam` copies the values of that system into the generic ones (`SetKind_ConvNumScheme`).
2. `CIntegration::Space_Integration` (`SU2_CFD/src/integration/CIntegration.cpp`) then calls `Centered_Residual` if `config->GetKind_ConvNumScheme()` is `SPACE_CENTERED`, and `Upwind_Residual` if it is `SPACE_UPWIND`.

## Two code paths: classic and vectorized

The compressible solver `CEulerSolver` and its child `CNSSolver` can compute the convective fluxes on the edges in two ways. This section tells you which one your scheme needs.

1. **Classic numerics.** A [`CNumerics`](../../doxygen/classCNumerics.html) object computes the flux of one edge. The solver does the MUSCL reconstruction and calls `ComputeResidual`. It then adds the result to the residual and the Jacobian. Files: `SU2_CFD/include/numerics/flow/convection/`.
2. **Vectorized numerics.** A [`CNumericsSIMD`](../../doxygen/classCNumericsSIMD.html) object computes the fluxes of several edges at once. It uses SIMD: one CPU instruction acts on several values together. The object also does the reconstruction and the viscous flux, and updates the residual and the Jacobian itself. Files: `SU2_CFD/include/numerics_simd/`.

The table shows which path is used. See `CEulerSolver::Centered_Residual` and `CEulerSolver::Upwind_Residual` in `SU2_CFD/src/solvers/CEulerSolver.cpp`.

| Case | Path |
|---|---|
| Centered schemes `JST`, `JST_KE`, `JST_MAT`, `LAX-FRIEDRICH` | Vectorized, always. Works with any fluid model. `LOW_MACH_CORR= YES` is rejected. Coarse multigrid levels use the Lax-Friedrich scheme. |
| Upwind `ROE` or `MSW`, `FLUID_MODEL` = `STANDARD_AIR` or `IDEAL_GAS`, `LOW_MACH_CORR= NO` | Vectorized |
| Any other upwind case (other schemes, other fluid models, `LOW_MACH_CORR= YES`) | Classic |
| Boundary conditions (`conv_bound_term`) | Classic, always |

This is the start of `CEulerSolver::Upwind_Residual`:

```cpp
  const bool ideal_gas = (config->GetKind_FluidModel() == STANDARD_AIR) ||
                         (config->GetKind_FluidModel() == IDEAL_GAS);
  const bool low_mach_corr = config->Low_Mach_Correction();

  /*--- Use vectorization if the scheme supports it. ---*/
  if (ideal_gas && !low_mach_corr) {
    if (config->GetKind_Upwind_Flow() == UPWIND::ROE || config->GetKind_Upwind_Flow() == UPWIND::MSW) {
      EdgeFluxResidual(geometry, solver_container, config);
      return;
    }
  }
```

A new upwind scheme needs at least a classic version, because the boundary conditions use the `CONV_BOUND_TERM` numerics. A vectorized version is optional. Write one only for schemes that will be used a lot with an ideal gas.

## The CNumerics interface

Every classic flux class derives from [`CNumerics`](../../doxygen/classCNumerics.html) (`SU2_CFD/include/numerics/CNumerics.hpp`). A flux class:

- takes `(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)` in its constructor, and passes them to `CNumerics(val_nDim, val_nVar, config)`;
- allocates its flux vector and its two Jacobian matrices once, in the constructor, and frees them in the destructor;
- gets the data of the edge through setters. The solver calls them before each call: `SetNormal`, `SetPrimitive`, `SetSecondary`, `SetGridVel`, `SetCoord`, `SetSensor`, `SetDissipation`, `SetVelocity2_Inf`, ...;
- overrides `ComputeResidual(const CConfig* config)`, which returns a `ResidualType<>`.

`ResidualType` is a read-only view of the flux and of the two Jacobians. The data stays inside the numerics object.

```cpp
  template<class Vector_t = const su2double*,
           class Matrix_t = const Vector_t*>
  struct ResidualType {
    const Vector_t residual;
    const Matrix_t jacobian_i;
    const Matrix_t jacobian_j;
    // ...
  };
```

The flow solvers use these conventions:

- The normal (`Normal`) points from point i to point j. Its norm is the area of the face. The flux is the flux that leaves the control volume of i, integrated over the face.
- `jacobian_i` and `jacobian_j` are the derivatives of the flux with respect to the **conservative** variables of i and j.
- The primitive vector of the compressible solver is `(T, vx, vy, [vz,] P, rho, h, c, ...)`. `CEulerVariable::CIndices` gives the positions (`Pressure()`, `Density()`, `Enthalpy()`, ...). Older numerics use `V_i[nDim+1]` for the pressure, `V_i[nDim+2]` for the density and `V_i[nDim+3]` for the enthalpy.
- For a non-ideal fluid model, `S_i[0]` is dP/drho at constant e, and `S_i[1]` is dP/de at constant rho. See `CEulerVariable::SetSecondaryVar`.
- Each OpenMP thread has its own numerics object (`numerics[...][CONV_TERM + offset]`). So the class can keep work arrays as members.

`CNumerics` also has older `ComputeResidual` overloads. They write into arrays passed by the caller, for example `ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig* config)`. The template classes in `SU2_CFD/include/numerics/template.hpp` still use them. The flow solvers call the `ResidualType<>` version, so use that one.

`CNumerics` has helpers for the Euler flux and its eigensystem: `GetInviscidProjFlux`, `GetInviscidProjJac`, `GetPMatrix`, `GetPMatrix_inv`. `GetRoe_Dissipation` serves the low-dissipation Roe variants.

### Example: the Roe family

The Roe schemes in `SU2_CFD/include/numerics/flow/convection/roe.hpp` share an intermediate base class, [`CUpwRoeBase_Flow`](../../doxygen/classCUpwRoeBase__Flow.html). Its `ComputeResidual` works in two steps:

1. It computes the central part of the flux.
2. It calls the pure virtual `FinalizeResidual`. Each child adds its own dissipation there: [`CUpwRoe_Flow`](../../doxygen/classCUpwRoe__Flow.html), `CUpwL2Roe_Flow`, `CUpwLMRoe_Flow`.

The `AD::...Preacc` calls belong to the discrete adjoint. See [Discrete adjoint](#discrete-adjoint).

Trimmed from `SU2_CFD/src/numerics/flow/convection/roe.cpp`:

```cpp
CNumerics::ResidualType<> CUpwRoeBase_Flow::ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  unsigned short iVar, jVar, iDim;
  su2double ProjGridVel = 0.0, Energy_i, Energy_j;

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+4); AD::SetPreaccIn(V_j, nDim+4); AD::SetPreaccIn(Normal, nDim);
  // ...

  /*--- Face area (norm of the normal vector) and unit normal ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Primitive variables at point i ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_i[iDim] = V_i[iDim+1];
  Pressure_i = V_i[nDim+1];
  Density_i  = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;

  // ...

  /*--- Negative RoeSoundSpeed^2, the jump variables is too large, clear fluxes and exit. ---*/

  if (RoeSoundSpeed2 <= 0.0) {
    // ...
    AD::SetPreaccOut(Flux, nVar);
    AD::EndPreacc();

    return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
  }

  // ...

  /*--- Apply Mavriplis' entropy correction to eigenvalues ---*/

  su2double MaxLambda = fabs(ProjVelocity) + RoeSoundSpeed;

  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = max(fabs(Lambda[iVar]), config->GetEntropyFix_Coeff()*MaxLambda);

  // ...

  /*--- Initialize residual (flux) and Jacobians ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);

  if (implicit) {
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, Jacobian_j);
  }

  /*--- Finalize in children class ---*/

  FinalizeResidual(Flux, Jacobian_i, Jacobian_j, config);

  /*--- Correct for grid motion ---*/
  // ...

  AD::SetPreaccOut(Flux, nVar);
  AD::SetPreaccOut(Jacobian_i, nVar, nVar);
  AD::SetPreaccOut(Jacobian_j, nVar, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}
```

The AUSM+ and SLAU family in `ausm_slau.hpp` uses the same pattern. `CUpwAUSMPLUS_SLAU_Base_Flow::ComputeResidual` is `final`, and the children implement `ComputeMassAndPressureFluxes`. If your scheme belongs to one of these families, derive from the base class and implement only the part that is specific to your scheme.

## Add a classic upwind scheme

Follow these steps to add an upwind scheme that works with every fluid model and on the boundaries.

1. Add the enum value and the keyword (see [above](#how-the-configuration-picks-the-scheme)).
2. Declare the class in a header in `SU2_CFD/include/numerics/flow/convection/`. Use an existing header or a new one. Derive the class from `CNumerics` or from a family base class. Override `ResidualType<> ComputeResidual(const CConfig* config)`.
3. Implement it in the matching `.cpp` in `SU2_CFD/src/numerics/flow/convection/`. If you add a new `.cpp` file, add it to the `numerics/...` list in `SU2_CFD/src/meson.build`.
4. Create the object in `SU2_CFD/src/drivers/CDriver.cpp`:
   - Include the header. The convection headers are included at the top of the file.
   - Add a `case` in the `switch (config->GetKind_Upwind_Flow())` of the compressible block of `CDriver::InitializeNumerics`.
   - Create the object for every multigrid level, both as `conv_term` and as `conv_bound_term`.

   For example, the MSW case:

   ```cpp
            case UPWIND::MSW:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
              }
              break;
   ```

   `InitializeNumerics` runs this code once per OpenMP thread (`conv_term = CONV_TERM + offset`), so each thread gets its own objects. `CDriver::FinalizeNumerics` deletes them.

   If the scheme works only with an ideal gas, follow the `UPWIND::ROE` and `UPWIND::HLLC` cases. They test `ideal_gas`, and for other fluids they create `CUpwGeneralRoe_Flow` or `CUpwGeneralHLLC_Flow` instead.
5. The scheme may need data that `CEulerSolver::Upwind_Residual` does not set yet. In that case, add the `Set...` call there, in the edge loop, before `numerics->ComputeResidual(config)`. If it needs a new quantity per point (a sensor, for example), compute it in `CEulerSolver::CommonPreprocessing`, next to `SetCentered_Dissipation_Sensor` and `SetRoe_Dissipation`.
6. Optional: print the scheme name in `CConfig::SetOutput`, next to the other `Kind_Upwind_Flow` lines. Add checks in `CConfig::SetPostprocessing`. For example, non-ideal fluids accept only `ROE`, `HLLC` and `JST` (search for "Non-Ideal Compressible Fluids" in `CConfig.cpp`).

Nothing else is needed. The solver gets the object from `numerics_container[CONV_TERM + omp_get_thread_num()*MAX_TERMS]`, calls it on every edge, and adds the result:

```cpp
    auto conv_residual = numerics->ComputeResidual(config);
    // ...
    if (ReducerStrategy) {
      EdgeFluxes.SetBlock(iEdge, conv_residual);
    } else {
      LinSysRes.AddBlock(iPoint, conv_residual);
      LinSysRes.SubtractBlock(jPoint, conv_residual);
    }

    /*--- Viscous contribution, returns its Jacobians so that the matrix is updated once. ---*/

    const auto visc_residual = Viscous_Residual(
        iEdge, geometry, solver_container, numerics_container[VISC_TERM + omp_get_thread_num()*MAX_TERMS], config);

    if (implicit) UpdateJacobian(iEdge, iPoint, jPoint, conv_residual, visc_residual);
```

`UpdateJacobian` (`CFVMFlowSolverBase.hpp`) updates four blocks of the matrix, together with the viscous Jacobians:

- `jacobian_i` is added to block (i,i) and subtracted from block (j,i);
- `jacobian_j` is added to block (i,j) and subtracted from block (j,j).

## MUSCL reconstruction and limiters

Your scheme does not reconstruct. The solver reconstructs the left and right states and passes them through `SetPrimitive` (and `SetSecondary`). This section shows where that happens.

- **Gradients and limiters.** `CEulerSolver::Preprocessing` computes them when `MUSCL_FLOW= YES`, on the finest grid. It calls `SetPrimitive_Gradient_GG` or `SetPrimitive_Gradient_LS` (option `NUM_METHOD_GRAD_RECON`), then `SetPrimitive_Limiter` (option `SLOPE_LIMITER_FLOW`). The limiter functions are in `SU2_CFD/include/limiters/` (`computeLimiters.hpp`, `CLimiterDetails.hpp`).
- **Reconstruction.** The edge loop of `CEulerSolver::Upwind_Residual` calls `CSolver::MUSCL_Reconstruction` (`SU2_CFD/include/solvers/CSolver.hpp`). It supports `MUSCL_KAPPA_FLOW` and the MUSCL ramp (`RAMP_MUSCL`). With `VAN_ALBADA_EDGE`, the limiter is computed per edge (`LimiterHelpers<>::vanAlbadaFunction`).
- **Non-ideal fluids and `LOW_MACH_CORR= YES`.** `CEulerSolver::ComputeConsistentExtrapolation` (and `LowMachPrimitiveCorrection`) makes the reconstructed states thermodynamically consistent.
- **Non-physical states.** If the reconstruction gives a negative pressure, density or Roe-averaged speed of sound, the edge falls back to first order (`UpdateNonPhysicalEdgeCounter`).

In the vectorized path, `reconstructPrimitives` in `SU2_CFD/include/numerics_simd/flow/convection/common.hpp` does the same steps. `CUpwindBase::ComputeFlux` calls it.

The incompressible solver (`CIncEulerSolver::Upwind_Residual`) and the NEMO solver (`CNEMOEulerSolver::Upwind_Residual`) do their own reconstruction in the same way.

## Add a vectorized scheme

This section shows how to add the vectorized (SIMD) version of a scheme. The framework is in `SU2_CFD/include/numerics_simd/`:

| File | Content |
|---|---|
| `CNumericsSIMD.hpp` | The interface: the pure virtual `ComputeFlux(Int iEdge, ...)` and the factory `CNumericsSIMD::CreateNumerics` |
| `CNumericsSIMD.cpp` | The factory functions `createUpwindIdealNumerics`, `createUpwindGeneralNumerics` (returns `nullptr`, no real-gas schemes yet), `createCenteredNumerics`, `createNumerics<nDim>`. It is compiled through `SU2_CFD/src/meson.build`. |
| `flow/convection/upwind.hpp` | `CUpwindBase`, `CRoeScheme`, `CMSWScheme` |
| `flow/convection/centered.hpp` | `CCenteredBase`, `CJSTScheme`, `CJSTmatScheme`, `CJSTkeScheme`, `CLaxScheme` |
| `flow/convection/common.hpp` | `reconstructPrimitives`, `pMatrix`, `pMatrixInv`, `inviscidProjFlux`, `inviscidProjJac`, `roeDissipation`, ... |
| `flow/diffusion/viscous_fluxes.hpp` | The viscous "decorators" `CNoViscousFlux`, `CCompressibleViscousFlux`, `CGeneralCompressibleViscousFlux` |
| `flow/variables.hpp` | `CCompressiblePrimitives`, `CCompressibleConservatives`, `roeAveragedVariables` |

The types `Double` and `Int` (`CNumericsSIMD.hpp`) are `simd::Array`s. Each operation on them acts on several edges at once. Helpers such as `gatherVariables`, `updateLinearSystem` and `stopPreacc` are in `SU2_CFD/include/numerics/util.hpp`.

Each scheme is a class template. Its template parameter is a viscous decorator: a class that adds the viscous flux, or no flux for `CNoViscousFlux`. The work is split this way:

- The base class (`CUpwindBase`) does the common work in `ComputeFlux`: geometry, reconstruction, conservative variables, viscous terms, and the update of the linear system.
- It then calls `finalizeFlux` of the derived class. The call is resolved at compile time through the template parameter, without virtual functions (static polymorphism).

This is the start of `CRoeScheme`:

```cpp
template<class Decorator>
class CRoeScheme : public CUpwindBase<CRoeScheme<Decorator>, Decorator> {
private:
  using Base = CUpwindBase<CRoeScheme<Decorator>, Decorator>;
  // ...
public:
  template<class... Ts>
  CRoeScheme(const CConfig& config, Ts&... args) : Base(config, args...),
    kappa(config.GetRoe_Kappa()),
    entropyFix(config.GetEntropyFix_Coeff()),
    typeDissip(static_cast<ENUM_ROELOWDISS>(config.GetKind_RoeLowDiss())) {
  }

  template<class PrimVarType, class ConsVarType, class... Ts>
  FORCEINLINE void finalizeFlux(VectorDbl<nVar>& flux,
                                MatrixDbl<nVar>& jac_i,
                                MatrixDbl<nVar>& jac_j,
                                const bool implicit,
                                const Double& area,
                                const VectorDbl<nDim>& unitNormal,
                                const VectorDbl<nDim>& normal,
                                const CPair<PrimVarType>& V,
                                const CPair<ConsVarType>& U,
                                const Int& iPoint,
                                const Int& jPoint,
                                const Double& nonPhysical,
                                const CEulerVariable& solution,
                                const CGeometry& geometry,
                                Ts&...) const {
    // ...
  }
};
```

To add a vectorized upwind scheme for an ideal gas:

1. In `upwind.hpp`, write `template<class Decorator> class CMyScheme : public CUpwindBase<CMyScheme<Decorator>, Decorator>`. Give it a constructor that forwards its arguments to the base, and a `finalizeFlux` with the signature above. Write all arithmetic with `Double`. Do not use `if` on values. Use `fmax`, `fmin`, `abs`, `sign` and masks instead, as `CMSWScheme` does.
2. Add a `case` in `createUpwindIdealNumerics` (`CNumericsSIMD.cpp`):

   ```cpp
     switch (config.GetKind_Upwind_Flow()) {
       case UPWIND::ROE:
         obj = new CRoeScheme<ViscousDecorator>(config, iMesh, turbVars);
         break;
       case UPWIND::MSW:
         obj = new CMSWScheme<ViscousDecorator>(config, iMesh, turbVars);
         break;
       default:
         break;
     }
   ```

3. Add the scheme to the test at the start of `CEulerSolver::Upwind_Residual` (shown [above](#two-code-paths-classic-and-vectorized)). Otherwise the classic path is used.
4. Keep the classic version in `CDriver::InitializeNumerics` for `conv_bound_term`.

How the solver calls it:

1. `CFVMFlowSolverBase::EdgeFluxResidual` (`SU2_CFD/include/solvers/CFVMFlowSolverBase.inl`) creates the object on first use (`CEulerSolver::InstantiateEdgeNumerics`).
2. It loops over the edge colors in chunks of `Double::Size` edges, and calls `edgeNumerics->ComputeFlux(...)`. Edge colors are groups of edges that threads can process in parallel.
3. If the factory returns `nullptr`, the run stops with "The numerical scheme + gas model in use do not support vectorization.".

With OpenMP, `EDGE_COLORING_GROUP_SIZE` must be a multiple of the SIMD length. The SIMD length depends on the build option `cpu-arch`.

The decorator adds the viscous flux in the same edge loop. So a viscous run with a vectorized scheme does not use `CAvgGrad_Flow` for the interior edges.

## Centered schemes

All compressible centered schemes are vectorized, so a new one goes into the vectorized framework. `CEulerSolver::Centered_Residual` only calls `EdgeFluxResidual`. `CDriver::InitializeNumerics` does not create a `conv_term` for these schemes. It creates only `CUpwRoe_Flow` as `conv_bound_term`.

To add a centered scheme:

1. Add the value to `enum class CENTERED` and `Centered_Map`.
2. Write the class in `SU2_CFD/include/numerics_simd/flow/convection/centered.hpp`, derived from `CCenteredBase<CMyScheme<Decorator>, Decorator>`. Follow `CJSTScheme` or `CLaxScheme`.
3. Add a `case` in `createCenteredNumerics` (`CNumericsSIMD.cpp`). On coarse grids (`iMesh != MESH_0`) the factory always uses `CENTERED::LAX`. The `switch` lists every `CENTERED` value and has no `default`. So add the new value there even if it is not vectorized.
4. The scheme may need the pressure sensor, the undivided Laplacian or the spectral radius. In that case, add it to the booleans `center_jst`, `center_jst_ke`, `center_jst_mat` in `CEulerSolver::CommonPreprocessing`. These booleans trigger the calls to `SetCentered_Dissipation_Sensor`, `SetUndivided_Laplacian` and `SetMax_Eigenvalue`.

The existing centered schemes use these options: `JST_SENSOR_COEFF`, `LAX_SENSOR_COEFF`, `CENTRAL_JACOBIAN_FIX_FACTOR`.

`CENTERED::LD2` and `CENTERED::CDS` exist only for the incompressible solver. For the compressible solver, they make the factory return `nullptr`.

## Jacobians

With `TIME_DISCRE_FLOW= EULER_IMPLICIT`, the solver needs the Jacobians of the flux. They go into the block matrix `Jacobian` (`CSysMatrix<su2mixedfloat>`). `ImplicitEuler_Iteration` solves this matrix together with the time-step term.

- The Jacobians do not need to be exact. Most schemes use an approximation:
  - the Roe schemes use `GetInviscidProjJac` and the frozen dissipation matrix;
  - the AUSM+ and SLAU family uses the Roe Jacobian (`ApproximateJacobian`). With `USE_ACCURATE_FLUX_JACOBIANS= YES` it uses `AccurateJacobian` instead, which mixes finite differences and manual differentiation.

  A more accurate Jacobian usually lets you use a higher CFL number.
- The matrix stores passive values (`PassiveAssign` in `CSysMatrix::UpdateBlocks`). Passive values are plain numbers that AD does not track. So AD does not differentiate the Jacobians.
- Fill the Jacobians only when `implicit` is true. This saves time in explicit runs.
- With `NEWTON_KRYLOV= YES`, `CNewtonIntegration` uses matrix-free products with the true Jacobian, computed by finite differences of the residual. It uses this matrix only to build the preconditioner. Expect less accurate finite differences if your flux is not smooth.

## Low-Mach and entropy fix options

These options change the convective fluxes. All are in `Common/src/CConfig.cpp`.

| Option | Getter | Used by |
|---|---|---|
| `ENTROPY_FIX_COEFF` | `GetEntropyFix_Coeff()` | Roe schemes (Mavriplis' correction of the eigenvalues) |
| `ROE_KAPPA` | `GetRoe_Kappa()` | Roe schemes (weight of the dissipation) |
| `ROE_LOW_DISSIPATION` | `GetKind_RoeLowDiss()` | `CUpwRoe_Flow`, `CUpwSLAU_Flow`, `CUpwSLAU2_Flow`, `CRoeScheme`; values in `RoeLowDiss_Map` |
| `LOW_MACH_CORR` | `Low_Mach_Correction()` | Correction of the reconstructed velocity in `CEulerSolver::LowMachPrimitiveCorrection`; not supported by the vectorized path |
| `LOW_MACH_PREC` | `Low_Mach_Preconditioning()` | Low-Mach preconditioning of the implicit time term in `CEulerSolver::PrepareImplicitIteration` (always on with `TURKEL_PREC`) |
| `MSW_ALPHA` | `GetMSW_Alpha()` | MSW scheme |
| `USE_ACCURATE_FLUX_JACOBIANS` | `GetUse_Accurate_Jacobians()` | AUSM+-up, AUSM+-up2, SLAU, SLAU2 |

How `ROE_LOW_DISSIPATION` works:

1. The solver computes the blending in `CEulerSolver::CommonPreprocessing` (`SetRoe_Dissipation`, `SetUpwind_Ducros_Sensor`).
2. It passes the result to the numerics with `SetDissipation` and `SetSensor`.
3. The numerics combine them with `CNumerics::GetRoe_Dissipation`.

A new scheme that supports this option takes a `bool val_low_dissipation` in its constructor, as `CUpwRoe_Flow` does. `CDriver` passes `roe_low_dissipation` for `conv_term` and `false` for `conv_bound_term`.

## Incompressible and NEMO schemes

The other finite volume solvers have their own scheme classes.

**Incompressible, density based** (`CIncEulerSolver`):

- Schemes: upwind `FDS` (`CUpwFDSInc_Flow`, `SU2_CFD/include/numerics/flow/convection/fds.hpp`); centered `JST`/`LD2` (`CCentJSTInc_Flow`) and `LAX-FRIEDRICH` (`CCentLaxInc_Flow`) in `centered.hpp`.
- They are created in the `incompressible` part of the convective `switch` in `CDriver::InitializeNumerics`.
- `CIncEulerSolver::Upwind_Residual` and `CIncEulerSolver::Centered_Residual` call them.
- There is no vectorized path.

**Incompressible, pressure based** (`INCOMP_SYSTEM::PRESSURE_BASED`): `UDS` (`CPBConvection_Upwind`) and `CDS` (`CPBConvection_Central`) in `pressure_based.hpp`.

**NEMO** (`CNEMOEulerSolver`):

- The schemes derive from [`CNEMONumerics`](../../doxygen/classCNEMONumerics.html) and are in `SU2_CFD/include/numerics/NEMO/convection/`: `CUpwRoe_NEMO`, `CUpwAUSM_NEMO`, `CUpwAUSMPLUSUP2_NEMO`, `CUpwAUSMPLUSM_NEMO`, `CUpwMSW_NEMO`, `CCentLax_NEMO`.
- They are created in the `NEMO_euler || NEMO_ns` block of `CDriver::InitializeNumerics`.
- Their constructors also take `nPrimVar` and `nPrimVarGrad`.

**FEM (DG)**: `RIEMANN_SOLVER_FEM` chooses the Riemann solver. `ROE` and `LAX-FRIEDRICH` are hard coded in `CFEM_DG_EulerSolver`. The others reuse the classic finite volume classes.

## Viscous fluxes and source terms

This section shows where to change the viscous flux and how to add a source term.

**Viscous flux.**

- The classic classes are in `SU2_CFD/include/numerics/flow/flow_diffusion.hpp`: `CAvgGrad_Base`, `CAvgGrad_Flow` (ideal gas), `CGeneralAvgGrad_Flow` (other fluid models), `CAvgGradInc_Flow`.
- They are created as `visc_term` and `visc_bound_term`, after the convective `switch` in `CDriver::InitializeNumerics`.
- The flow solvers call them per edge, from the convective loop, through `Viscous_Residual(iEdge, ...)` (`CFVMFlowSolverBase::Viscous_Residual_impl`). The residual is subtracted from i and added to j.
- In the compressible solver, the vectorized path uses the decorators in `numerics_simd/flow/diffusion/viscous_fluxes.hpp` instead. So make any change to the viscous flux in both places.

**Source term.** The flow source terms are in `SU2_CFD/include/numerics/flow/flow_sources.hpp`. They derive from `CSourceBase_Flow`, which allocates `residual` and `jacobian`. To add one:

1. Write the class with `ResidualType<> ComputeResidual(const CConfig* config) override`, returning `ResidualType<>(residual, jacobian, nullptr)`. The residual is minus the source integrated over the control volume (see `CSourceBodyForce::ComputeResidual`).
2. Create it in `CDriver::InitializeNumerics`, in the `if`/`else if` chain for `source_first_term` or `source_second_term`. Each solver has only these two source slots. `CSourceNothing` (an alias of `CNumerics`) fills an unused slot.
3. Call it in the `Source_Residual` of the solver (`CEulerSolver::Source_Residual`, `CIncEulerSolver::Source_Residual`):
   - loop over the points;
   - set the inputs (`SetConservative`, `SetPrimitive`, `SetVolume`, `SetCoord`, ...);
   - call `ComputeResidual`;
   - call `LinSysRes.AddBlock(iPoint, residual)` and, if implicit, `Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i)`.

## Discrete adjoint

A new scheme does not need adjoint code, but you must write it so that automatic differentiation (AD) works. The discrete adjoint drivers record one primal iteration, with all its flux evaluations, on the AD tape. The tape is the record of every operation. Evaluating the tape gives the adjoint. `CDiscAdjSolver` holds the adjoint variables.

Rules for the scheme:

- Use `su2double` for every value that depends on the solution or on the mesh. `passivedouble` and `SU2_TYPE::GetValue` cut the derivative.
- Use preaccumulation correctly. Preaccumulation condenses the operations of one function into a small Jacobian from its inputs to its outputs, so the tape stays smaller.
  - Start with `AD::StartPreacc()`.
  - Declare **all** inputs with `AD::SetPreaccIn`: primitives, normal, grid velocity, sensors, secondary variables, ...
  - End **every** return path with `AD::SetPreaccOut` and `AD::EndPreacc()`, as in the early return of `CUpwRoeBase_Flow::ComputeResidual` above.

  A missing input gives wrong gradients without any error. In the vectorized framework, `gatherVariables` registers the inputs and `stopPreacc` ends the section.
- Branches on solution values (`max`, `fabs`, `if`) are recorded as taken at the converged state. The derivative is the one of the active branch. Avoid branches whose outcome changes between iterations at convergence.
- The Jacobians are not differentiated, because the matrix is passive.
- A tape debug build helps to find preaccumulation errors. Build with `-Dcodi-tape=Tag` and use the options `CHECK_TAPE_TYPE` and `CHECK_TAPE_VARIABLES`. See [Advanced AD Techniques](../Advanced-AD-Techniques/).

## Tests

Add these tests for a new scheme:

- **Regression test.** Add a short case to `TestCases/parallel_regression.py` and `TestCases/serial_regression.py`. Existing scheme tests you can copy: `wedge` (`euler/wedge/inv_wedge_HLLC.cfg`) and `ramp_msw` (`euler/ramp/inv_ramp_msw.cfg`). If the scheme is meant for the discrete adjoint, also add a case to `parallel_regression_AD.py` / `serial_regression_AD.py` (for example, the `disc_adj_euler` cases). See [Running Regression Tests](../Running-Regression-Tests/).
- **Order of accuracy.** Use the method of manufactured solutions (MMS) cases in `TestCases/mms/`. The `mms_fvm_*` tests in `parallel_regression.py` use them:
  - `mms/fvm_euler/inv_mms_vortex.cfg` (`ROE`, `KIND_VERIFICATION_SOLUTION= INVISCID_VORTEX`)
  - `mms/fvm_navierstokes/lam_mms_roe.cfg` (`ROE`, `MMS_NS_UNIT_QUAD`)
  - `mms/fvm_incomp_euler/inv_mms_jst.cfg` (`JST`, `MMS_INC_EULER`)
  - `mms/fvm_incomp_navierstokes/lam_mms_fds.cfg` (`FDS`, `MMS_INC_NS`)

  Change `CONV_NUM_METHOD_FLOW` and run them on a sequence of meshes. The exact solutions are in `Common/include/toolboxes/MMS/`. `CFVMFlowSolverBase::ComputeVerificationError` computes the error norms.
- **Unit test.** You can test numerics alone with Catch2. See `UnitTests/SU2_CFD/numerics/CNumerics_tests.cpp` and [Writing Unit Tests](../Writing-Unit-Tests/).
- **OpenMP.** If you changed the solver loops, run at least one case with OpenMP threads (`hybrid_regression.py`).

## Checklist

- [ ] Enum value and keyword in `option_structure.hpp`, keyword in only one of `Centered_Map` / `Upwind_Map`.
- [ ] Keyword listed in `config_template.cfg`.
- [ ] Classic class with `ResidualType<> ComputeResidual(const CConfig*)`, flux and Jacobians w.r.t. conservative variables.
- [ ] New `.cpp` in `SU2_CFD/src/meson.build`.
- [ ] `conv_term` and `conv_bound_term` created for every multigrid level in `CDriver::InitializeNumerics`.
- [ ] Vectorized version (optional): class in `numerics_simd`, factory `case`, test in `CEulerSolver::Upwind_Residual`.
- [ ] AD: all inputs in preaccumulation, every return path closes it.
- [ ] Regression test, and MMS runs for the order of accuracy.

Last checked against SU2 develop `11c177b904`.
