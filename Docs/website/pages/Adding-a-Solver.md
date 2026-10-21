---
title: Adding a Solver
permalink: /docs_v7/Adding-a-Solver/
---

This page shows how to add a new solver to `SU2_CFD`. It covers three cases:

- a new set of equations with its own `SOLVER=` keyword;
- a new equation that runs next to the flow solver, like the species or turbulence solvers;
- a new `CSolver` class for physics that already exist.

It follows the code path from the configuration file to the restart file. Read [Code Structure](../Code-Structure/) first; this page does not repeat the class overview.

## Files you will touch

| File | What goes there |
|---|---|
| `Common/include/option_structure.hpp` | `enum class MAIN_SOLVER` and `Solver_Map`, the position in the solver container (`enum SOLVER_TYPE`, `MAX_SOLS`), the equation system (`enum RUNTIME_TYPE`), `enum class MPI_QUANTITIES` if you need new communications |
| `Common/include/CConfig.hpp`, `Common/src/CConfig.cpp` | new options, `CConfig::GetContainerPosition`, `CConfig::SetGlobalParam`, the discrete adjoint mapping in `CConfig::SetPostprocessing`, problem helpers such as `GetHeatProblem` |
| `SU2_CFD/include/solvers/CSolverFactory.hpp`, `SU2_CFD/src/solvers/CSolverFactory.cpp` | `enum class SUB_SOLVER_TYPE`, the solver creation and its `INTEGRATION_TYPE` |
| `SU2_CFD/include/solvers/CMySolver.hpp`, `SU2_CFD/src/solvers/CMySolver.cpp` | the new `CSolver` child |
| `SU2_CFD/include/variables/CMyVariable.hpp`, `SU2_CFD/src/variables/CMyVariable.cpp` | the new `CVariable` child |
| `SU2_CFD/src/drivers/CDriver.cpp` | numerics objects in `CDriver::InitializeNumerics` (if the solver uses `CNumerics`) |
| `SU2_CFD/src/iteration/CIterationFactory.cpp` and an iteration class | which `CIteration` runs the solver |
| `SU2_CFD/src/drivers/CSinglezoneDriver.cpp`, `SU2_CFD/src/drivers/CMultizoneDriver.cpp` | the call to `SetInitialCondition` |
| `SU2_CFD/src/output/COutputFactory.cpp` and an output class | screen, history and volume output |
| `SU2_CFD/src/meson.build` | the new `.cpp` files |
| `config_template.cfg`, `TestCases/` | documentation of new options, a regression test |

---

- [Choose the kind of solver and a template](#choose-the-kind-of-solver-and-a-template)
- [1. The SOLVER keyword](#1-the-solver-keyword)
- [2. Position in the solver container and equation system](#2-position-in-the-solver-container-and-equation-system)
- [3. Create the solver in CSolverFactory](#3-create-the-solver-in-csolverfactory)
- [4. The CVariable class](#4-the-cvariable-class)
- [5. The CSolver class](#5-the-csolver-class)
- [6. Numerics](#6-numerics)
- [7. Integration, iteration and initial condition](#7-integration-iteration-and-initial-condition)
- [8. Output](#8-output)
- [9. Restart files](#9-restart-files)
- [10. MPI communication](#10-mpi-communication)
- [11. OpenMP](#11-openmp)
- [12. Discrete adjoint](#12-discrete-adjoint)
- [13. Build and test](#13-build-and-test)
- [Non-conservative systems of equations](#non-conservative-systems-of-equations)
- [Checklist](#checklist)

---

## Choose the kind of solver and a template

First find which of three cases you are in. Each case needs a different subset of the steps below.

| Case | Example in the code | Steps |
|---|---|---|
| New main solver with its own `SOLVER=` keyword | `HEAT_EQUATION` (`CHeatSolver`) | all |
| New equation that runs after the flow solver in the same zone | species transport (`CSpeciesSolver`, option `KIND_SCALAR_MODEL`) | 2 to 13. No new keyword: a model option instead. |
| New class for existing physics | a turbulence model chosen by `KIND_TURB_MODEL` in `CSolverFactory::CreateTurbSolver` | 3 to 6, 8, 13 |

Then pick a class to copy:

- **`CTemplateSolver`** (`SU2_CFD/include/solvers/CTemplateSolver.hpp`) lists the main virtual functions with empty bodies. It is a list, not a working solver. Its constructor does not create `nodes`. `SOLVER= TEMPLATE_SOLVER` stops during set-up, at the latest in `CIterationFactory::CreateIteration` ("No iteration found for specified solver."). Use it only to see which functions exist.
- **[`CHeatSolver`](../../doxygen/classCHeatSolver.html)** is the smallest complete main solver. It has its own keyword (`HEAT_EQUATION`), iteration (`CHeatIteration`), output (`CHeatOutput`), restart and discrete adjoint (`DISC_ADJ_HEAT`). It also runs as a sub-solver of the incompressible solver (`WEAKLY_COUPLED_HEAT_EQUATION= YES`). Copy it for a new set of scalar equations.
- **[`CSpeciesSolver`](../../doxygen/classCSpeciesSolver.html)** is the template for an equation that reads the flow solution.
- Both derive from [`CScalarSolver`](../../doxygen/classCScalarSolver.html)`<VariableType>`. This base class already implements the implicit and explicit updates, dual time stepping, periodic boundaries, OpenMP and the scalar edge loop. Derive from it if your equations are transport equations for scalars. Derive from `CSolver` directly only for something different. Examples are `CRadSolver` and `CRadP1Solver`, or the flow solvers through [`CFVMFlowSolverBase`](../../doxygen/classCFVMFlowSolverBase.html).

## 1. The SOLVER keyword

Add a keyword so that users can select your solver with `SOLVER=`. The option `SOLVER` is read into `CConfig::Kind_Solver`:

```cpp
  addEnumOption("SOLVER", Kind_Solver, Solver_Map, MAIN_SOLVER::NONE);
```

Add a value to `enum class MAIN_SOLVER` and a keyword to `Solver_Map` in `Common/include/option_structure.hpp`:

```cpp
enum class MAIN_SOLVER {
  NONE,                        /*!< \brief Definition of no solver. */
  EULER,                       /*!< \brief Definition of the Euler's solver. */
  // ...
  HEAT_EQUATION,               /*!< \brief Definition of the finite volume heat solver. */
  // ...
  TEMPLATE_SOLVER,             /*!< \brief Definition of template solver. */
  // ...
  DISC_ADJ_HEAT,               /*!< \brief Definition of the discrete adjoint heat solver. */
  // ...
};
static const MapType<std::string, MAIN_SOLVER> Solver_Map = {
  // ...
  MakePair("HEAT_EQUATION", MAIN_SOLVER::HEAT_EQUATION)
  MakePair("ELASTICITY", MAIN_SOLVER::FEM_ELASTICITY)
  MakePair("TEMPLATE_SOLVER", MAIN_SOLVER::TEMPLATE_SOLVER)
  // ...
};
```

**Discrete adjoint.** The `DISC_ADJ_*` and `ADJ_*` values are not in the map. Users write the primal keyword and `MATH_PROBLEM= DISCRETE_ADJOINT`. `CConfig::SetPostprocessing` then replaces the value at its end:

```cpp
    /*--- Note that this is deliberately done at the end of this routine! ---*/
    switch(Kind_Solver) {
      // ...
      case MAIN_SOLVER::HEAT_EQUATION:
        Kind_Solver = MAIN_SOLVER::DISC_ADJ_HEAT;
        break;
      default:
        break;
    }
```

If your solver gets a discrete adjoint, add a `DISC_ADJ_` value and a case here.

**Problem helpers.** Many places ask which kind of problem runs. They use helpers in `Common/include/CConfig.hpp`: `GetFluidProblem()`, `GetStructuralProblem()`, `GetHeatProblem()`, `GetFEMSolver()`. Search for them and decide whether the new value belongs in one of them.

**Sub-models.** A new sub-model of an existing solver does not get a `SOLVER=` value. It gets its own enum option instead, for example:

- `KIND_TURB_MODEL` (`enum class TURB_MODEL`, `Turb_Model_Map`, grouped by `TurbModelFamily()`);
- `KIND_TRANS_MODEL` (`enum class TURB_TRANS_MODEL`);
- `KIND_SCALAR_MODEL` (`enum class SPECIES_MODEL`).

See [Adding a Config Option](../Adding-a-Config-Option/) for options and enums.

## 2. Position in the solver container and equation system

Give the solver a slot in the solver container and an equation system value.

The solver container is the array that holds the solvers of a zone. There is one array per multigrid level: `solver_container[iZone][iInst][iMesh][iSol]`. The positions `iSol` are in `enum SOLVER_TYPE` in `Common/include/option_structure.hpp`:

```cpp
 enum SOLVER_TYPE : const int {
   FLOW_SOL=0,       /*!< \brief Position of the mean flow solution in the solver container array. */
   ADJFLOW_SOL=1,    /*!< \brief Position of the continuous adjoint flow solution in the solver container array. */
   TURB_SOL=2,       /*!< \brief Position of the turbulence model solution in the solver container array. */
   // ...
   SPECIES_SOL=11,    /*!< \brief Position of the species solver. */
   ADJSPECIES_SOL=12, /*!< \brief Position of the adjoint of the species solver. */
   FEA_SOL=0,        /*!< \brief Position of the Finite Element flow solution in the solver container array. */
   ADJFEA_SOL=1,     /*!< \brief Position of the continuous adjoint Finite Element flow solution in the solver container array. */
   POISSON_SOL=13,   /*!< \brief Position of the poisson solution in the solver container array */
   TEMPLATE_SOL=0,   /*!< \brief Position of the template solution. */
 };
```

The array size is `MAX_SOLS`:

```cpp
enum: unsigned int{MAX_SOLS = 14};            /*!< \brief Maximum number of solutions at the same time (dimension of solution container array). */
```

Positions 0 to 13 are taken. You have two options:

- Reuse a slot, if your solver is never created together with the solver that holds it. `FEA_SOL` and `TEMPLATE_SOL` reuse 0.
- Otherwise, add a new position (and one for its adjoint) and increase `MAX_SOLS`.

The drivers, integrations and outputs loop over `MAX_SOLS` and skip null entries.

The solver also needs an **equation system** value in `enum RUNTIME_TYPE`, for example `RUNTIME_HEAT_SYS = 21`. The integrations receive this value. They turn it into the container position with `CConfig::GetContainerPosition` (`Common/src/CConfig.cpp`):

```cpp
unsigned short CConfig::GetContainerPosition(unsigned short val_eqsystem) {

  switch (val_eqsystem) {
    case RUNTIME_FLOW_SYS:      return FLOW_SOL;
    // ...
    case RUNTIME_HEAT_SYS:      return HEAT_SOL;
    // ...
  }
  return 0;
}
```

An unknown value returns 0, the flow solver, without an error. Add your case.

Then set the schemes of your system in `CConfig::SetGlobalParam`. `CIntegration::Space_Integration` and `CIntegration::Time_Integration` read the convective scheme (`GetKind_ConvNumScheme()`) and the time scheme (`GetKind_TimeIntScheme()`) of the *current* system. `CConfig::SetGlobalParam(MAIN_SOLVER, RUNTIME_...)` sets them. The iteration and the solver's `Preprocessing` call it before each system. Add a case for your solver and system:

```cpp
    case MAIN_SOLVER::HEAT_EQUATION:
    case MAIN_SOLVER::DISC_ADJ_HEAT:
      if (val_system == RUNTIME_HEAT_SYS) {
        SetKind_ConvNumScheme(NONE, CENTERED::NONE, UPWIND::NONE, LIMITER::NONE, NONE, 0.0, NONE);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Heat);
      }
      break;
```

With `NONE` as convective scheme, `Space_Integration` calls neither `Centered_Residual` nor `Upwind_Residual`. The heat solver does all its work in `Viscous_Residual`.

A sub-solver of the flow solver needs a case in the branch of each flow `MAIN_SOLVER` value it runs with. See the `SetSpeciesParam` and `SetHeatParam` lambdas in `SetGlobalParam`.

## 3. Create the solver in CSolverFactory

Make `CSolverFactory` create your solver and choose its integration class. `CDriver::InitializeSolver` calls `CSolverFactory::CreateSolverContainer` once per multigrid level, so each level has its own solver object:

```cpp
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++){
    solver[iMesh] = CSolverFactory::CreateSolverContainer(kindSolver, config, geometry[iMesh], iMesh);
  }
```

In `SU2_CFD/include/solvers/CSolverFactory.hpp`, add a value to `enum class SUB_SOLVER_TYPE`. Then, in `SU2_CFD/src/solvers/CSolverFactory.cpp`:

1. In `CSolverFactory::CreateSolverContainer`, add a case for the `MAIN_SOLVER` value. The case fills the positions of your solvers:

   ```cpp
       case MAIN_SOLVER::HEAT_EQUATION:
         solver[HEAT_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
         break;
       // ...
       case MAIN_SOLVER::DISC_ADJ_HEAT:
         solver[HEAT_SOL]    = CreateSubSolver(SUB_SOLVER_TYPE::HEAT, solver, geometry, config, iMGLevel);
         solver[ADJHEAT_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::DISC_ADJ_HEAT, solver, geometry, config, iMGLevel);
         break;
   ```

   For a sub-solver, add the line to the case of each flow solver it can run with (look at `SUB_SOLVER_TYPE::SPECIES`). After the switch, the mesh solvers are always added (`MESH_SOL`, `ADJMESH_SOL`).

2. In `CSolverFactory::CreateSubSolver`, create the object and set its `INTEGRATION_TYPE`:

   ```cpp
       case SUB_SOLVER_TYPE::HEAT:
         genericSolver = CreateHeatSolver(solver, geometry, config, iMGLevel, false);
         metaData.integrationType = INTEGRATION_TYPE::SINGLEGRID;
         break;
       case SUB_SOLVER_TYPE::DISC_ADJ_HEAT:
         genericSolver = CreateHeatSolver(solver, geometry, config, iMGLevel, true);
         metaData.integrationType = INTEGRATION_TYPE::DEFAULT;
         break;
   ```

   The sub-solver may return `nullptr`, for example when the model is off. Only non-null solvers are stored in the metadata map, and only they get an integration.

3. Add the new `#include` at the top of `SU2_CFD/src/solvers/CSolverFactory.cpp`.

`INTEGRATION_TYPE` (same header) picks the integration class in `CIntegrationFactory::CreateIntegration`:

| `INTEGRATION_TYPE` | Class | Used by |
|---|---|---|
| `MULTIGRID` | `CMultiGridIntegration` | finite volume flow solvers, continuous adjoint flow |
| `NEWTON` | `CNewtonIntegration` | flow solvers with `NEWTON_KRYLOV= YES` (primal only) |
| `SINGLEGRID` | `CSingleGridIntegration` | turbulence, transition, species, heat, radiation, Poisson, template |
| `STRUCTURAL` | `CStructuralIntegration` | `CFEASolver` |
| `FEM_DG` | `CFEM_DG_Integration` | DG flow solvers |
| `DEFAULT` | `CIntegration` | discrete adjoint solvers, baseline solvers |
| `NONE` | none | mesh solver |

## 4. The CVariable class

Write a `CVariable` child to hold the per-point data of your solver. [`CVariable`](../../doxygen/classCVariable.html) stores the data of all points of one solver (see [Code Structure](../Code-Structure/#variables)).

For scalar equations, derive from [`CScalarVariable`](../../doxygen/classCScalarVariable.html). It allocates gradients, limiters, the local time step, the under-relaxation and the local CFL. `CHeatVariable` is the shortest example:

```cpp
class CHeatVariable final : public CScalarVariable {
public:
  static constexpr size_t MAXNVAR = 1; /*!< \brief Max number of variables, for static arrays. */
  // ...
  CHeatVariable(su2double heat, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);
  // ...
  inline su2double GetTemperature(unsigned long iPoint) const final { return Solution(iPoint, 0); }
};
```

`CScalarSolver<VariableType>` needs `MAXNVAR`: it sizes its static arrays with `VariableType::MAXNVAR`.

The constructor (`SU2_CFD/src/variables/CHeatVariable.cpp`) sets `Solution` and `Solution_Old`. It allocates only what the solver uses:

- `Res_TruncError`;
- `Solution_time_n` and `Solution_time_n1`, for dual time stepping;
- `Max_Lambda_Visc`;
- `Set_BGSSolution_k()`, for multizone problems.

Store extra per-point quantities as members of the variable class (`su2activevector`, `su2activematrix`). Give them getters and setters that take the point index.

## 5. The CSolver class

The solver class holds the residual loops, the solution update and the restart. This section covers its constructor, the virtual functions that the integration calls, and the implicit update.

### Constructor

Follow the order of `CHeatSolver::CHeatSolver` (`SU2_CFD/src/solvers/CHeatSolver.cpp`):

1. Call the base constructor. For `CScalarSolver` it is `CScalarSolver(geometry, config, flow_solver, conservative, bounded_scalar)`. `flow_solver` may be `nullptr` when there is no flow.
2. Set `nVar`, `nPrimVar`, `nPoint`, `nPointDomain`, `nVarGrad`, and size `Residual_RMS`, `Residual_Max`, `Point_Max`, `Point_Max_Coord`.
3. Initialize the linear system. `ReducerStrategy` is one of the two OpenMP strategies for edge loops (see [step 11](#11-openmp)); it needs per-edge flux storage:

   ```cpp
     Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
     LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
     LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
     if (ReducerStrategy) {
       EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
       EdgeFluxesDiff.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
     }
   ```

4. Create the variables and set the base class pointer. `GetBaseClassPointerToNodes()` is the only pure virtual function of `CSolver`, and `CScalarSolver` implements it. `CSolver::GetNodes()` asserts that the pointer was set.

   ```cpp
     nodes = new CHeatVariable(Solution_Inf[0], nPoint, nDim, nVar, config);

     SetBaseClassPointerToNodes();
   ```

5. Set up the periodic and MPI communication:
   - `InitiatePeriodicComms` / `CompletePeriodicComms` with `PERIODIC_VOLUME` and `PERIODIC_NEIGHBORS`;
   - `SetImplicitPeriodic(...)`;
   - `InitiateComms` / `CompleteComms` with `MPI_QUANTITIES::SOLUTION`.
6. Set the local CFL of every point, and `SolverName` (used in the output).

### Virtual functions

All are declared in `SU2_CFD/include/solvers/CSolver.hpp` with an empty default body, except `GetBaseClassPointerToNodes`. The integration calls them in the order of this table (see `CSingleGridIntegration::SingleGrid_Iteration` and `CIntegration::Space_Integration`):

| Function | Called by | What it does | In `CScalarSolver`? |
|---|---|---|---|
| `Preprocessing(geometry, solver_container, config, iMesh, iRKStep, RunTime_EqSystem, Output)` | integration, restart, adjoint recording | call `config->SetGlobalParam`, compute gradients and limiters, clear residual and Jacobian (`CScalarSolver::CommonPreprocessing`) | helper only |
| `Set_OldSolution()` | `SingleGrid_Iteration` | copy `Solution` to `Solution_Old` | yes (`CSolver`) |
| `SetTime_Step(geometry, solver_container, config, iMesh, Iteration)` | integration | local time step | copies the flow time step; override for a standalone solver (`CHeatSolver` does) |
| `Centered_Residual` / `Upwind_Residual` | `Space_Integration`, by `GetKind_ConvNumScheme()` | convective edge loop | no, implement it |
| `Viscous_Residual(geometry, solver_container, numerics_container, config, iMesh, iRKStep)` | `Space_Integration`, always | diffusive edge loop | no |
| `Source_Residual(geometry, solver_container, numerics_container, config, iMesh)` | `Space_Integration`, always | point source terms | no |
| `SetResidual_DualTime(...)` | `Space_Integration`, dual time only | time derivative terms | yes |
| `BC_...` | `Space_Integration` | boundary conditions, see [Adding a Boundary Condition](../Adding-a-Boundary-Condition/) | some |
| `ImplicitEuler_Iteration` / `ExplicitEuler_Iteration` / `ExplicitRK_Iteration` / `ClassicalRK4_Iteration` | `Time_Integration`, by `GetKind_TimeIntScheme()` | update the solution | implicit and explicit Euler |
| `Postprocessing(geometry, solver_container, config, iMesh)` | integration | quantities derived from the new solution | no |
| `SetInitialCondition(geometry, solver_container, config, TimeIter)` | driver, each time step | initial or analytic condition | no |
| `LoadRestart(geometry, solver, config, val_iter, val_update_geo)` | `CDriver::RestartSolver` | read the restart file | pure virtual in `CScalarSolver` |
| `GetHasHybridParallel()` | integration, restart | whether the solver supports OpenMP | `true` |

`SingleGrid_Iteration` also contains code for specific systems. Examples are the heat fluxes for `RUNTIME_HEAT_SYS`, and the restriction of the turbulence solution to the coarse levels for `RUNTIME_TURB_SYS`. Check it if your solver needs something similar.

### Implicit iteration

A solver derived from `CSolver` directly must implement `ImplicitEuler_Iteration`. The function must:

1. add the time term to the diagonal of `Jacobian`;
2. copy `-LinSysRes` into the right-hand side;
3. call `System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config)` (see the existing implementations);
4. update the solution;
5. communicate the solution.

`CScalarSolver` splits this into `PrepareImplicitIteration` and `CompleteImplicitIteration` (both `final`). Use them as the reference.

Check `config->GetKind_TimeIntScheme() == EULER_IMPLICIT` inside each function, not once in the constructor. The reason: `CNewtonIntegration::ComputeResiduals` temporarily switches the scheme to `EULER_EXPLICIT`, to evaluate residuals without the Jacobian.

## 6. Numerics

The numerics compute the fluxes and source terms. There is no numerics factory. The code uses two styles:

- `CNumerics` objects, for the flow solvers, source terms and radiation;
- scalar edge kernels, for all `CScalarSolver` children.

**CNumerics objects.** `CDriver::InitializeNumerics` (`SU2_CFD/src/drivers/CDriver.cpp`) creates them, one per thread and per multigrid level. It stores them in `numerics[iMGlevel][iSol][iTerm + offset]`, with `offset = thread * MAX_TERMS`. The term indices are `CONV_TERM`, `VISC_TERM`, `SOURCE_FIRST_TERM`, `SOURCE_SECOND_TERM`, `CONV_BOUND_TERM`, `VISC_BOUND_TERM` and `GRAD_TERM`. The radiation solver is a short example:

```cpp
  if (config->AddRadiation()) {
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    numerics[MESH_0][RAD_SOL][VISC_TERM] = new CAvgGradCorrected_P1(nDim, nVar_Rad, config);

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    numerics[MESH_0][RAD_SOL][SOURCE_FIRST_TERM] = new CSourceP1(nDim, nVar_Rad, config);

    /*--- Definition of the boundary condition method ---*/
    numerics[MESH_0][RAD_SOL][VISC_BOUND_TERM] = new CAvgGradCorrected_P1(nDim, nVar_Rad, config);
  }
```

To add your solver:

- `InitializeNumerics` first sets booleans from `config->GetKind_Solver()` (`template_solver`, `heat`, `turbulent`, ...) and reads `nVar` from the solvers. Add your solver there.
- The helpers `InstantiateTurbulentNumerics`, `InstantiateTransitionNumerics` and `InstantiateSpeciesNumerics` create the turbulence, transition and species source terms.
- The solver receives the array `numerics_container` for its position. It picks its object with `numerics_container[CONV_TERM + omp_get_thread_num()*MAX_TERMS]`.
- Template numerics classes are in `SU2_CFD/include/numerics/template.hpp` (`CConvective_Template`, `CViscous_Template`, `CSource_Template`).

**Scalar edge kernels.** An edge kernel is a small class that computes the flux of one edge. In all `CScalarSolver` children, a kernel computes the convection, diffusion and boundary fluxes, not `CNumerics`. The kernel derives from `CUpwScalarBase` (`SU2_CFD/include/numerics/scalar/scalar_edge_flux.hpp`). `InitializeNumerics` creates only the source terms of these solvers. `CONV_TERM`, `VISC_TERM`, `CONV_BOUND_TERM` and `VISC_BOUND_TERM` stay `nullptr`. The solver calls its kernel from `Upwind_Residual` or `Viscous_Residual`:

```cpp
void CHeatSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                                  CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {
  // ...
  const auto opt = ScalarFluxOptions::Interior(*config, false);

  DispatchScheme<CScalarFlux_Heat, 1>(config, [&](auto tag) {
    EdgeFluxResidual<typename decltype(tag)::type>(geometry, solver_container, config, opt);
  });
}
```

`DispatchScheme` picks, at compile time, the compressible or incompressible flow indices, the dimension and the number of equations.

Kernels exist in `SU2_CFD/include/numerics` for:

- SA (`turbulent/turb_sa_edge_flux.hpp`);
- SST (`turbulent/turb_sst_edge_flux.hpp`);
- LM (`turbulent/transition/trans_edge_flux.hpp`);
- species and flamelet (`species/`);
- heat (`heat_edge_flux.hpp`);
- Poisson (`poisson_edge_flux.hpp`).

See [Non-conservative systems of equations](#non-conservative-systems-of-equations) for how a kernel fills the residual.

## 7. Integration, iteration and initial condition

Choose the iteration class that runs your solver, and set up its initial condition.

The integration is set by `INTEGRATION_TYPE` (step 3). `CIterationFactory::CreateIteration` (`SU2_CFD/src/iteration/CIterationFactory.cpp`) creates the iteration; `CDriver::PreprocessIteration` calls it. Add a case for the new `MAIN_SOLVER` value:

```cpp
    case MAIN_SOLVER::HEAT_EQUATION:
      if (rank == MASTER_NODE)
        cout << "Heat iteration (finite volume method)." << endl;
      iteration = new CHeatIteration(config);
      break;
```

`NONE`, `TEMPLATE_SOLVER`, `MULTIPHYSICS` and `POISSON_EQUATION` end in `SU2_MPI::Error("No iteration found for specified solver.", ...)`. The switch has no `default`. So a new value without a case returns a null iteration, without an error. Add the case.

**New iteration class.** Copy [`CHeatIteration`](../../doxygen/classCHeatIteration.html), a child of `CFluidIteration`. Override `Iterate` and `Update`:

```cpp
void CHeatIteration::Iterate(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                             CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                             CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                             unsigned short val_iInst) {
  SU2_ZONE_SCOPED

  /*--- Update global parameters ---*/

  config[val_iZone]->SetGlobalParam(MAIN_SOLVER::HEAT_EQUATION, RUNTIME_HEAT_SYS);

  integration[val_iZone][val_iInst][HEAT_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                    RUNTIME_HEAT_SYS, val_iZone, val_iInst);
}
```

For dual time stepping, `CHeatIteration::Update` calls `SetDualTime_Solver` and `SetDualTime_Geometry` of the integration on every level.

**Sub-solver.** A sub-solver of the flow solver needs no new iteration. Add it to `CFluidIteration::CommonAuxiliarySolvers` (`SU2_CFD/src/iteration/CFluidIteration.cpp`), which runs after the flow solver:

```cpp
  if (config[val_iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
    config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_SPECIES_SYS);
    integration[val_iZone][val_iInst][SPECIES_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                         RUNTIME_SPECIES_SYS, val_iZone, val_iInst);
    // ...
  }
```

Also look at `CFluidIteration::Update` for the dual time update of each solver.

**Initial condition.** The drivers call `SetInitialCondition` of each solver at every time step. The calls are in `CSinglezoneDriver::Preprocess` and `CMultizoneDriver::Preprocess`, and for harmonic balance in `CFluidDriver::Preprocess` (in `CDriver.cpp`). They test `GetFluidProblem()`, `GetHeatProblem()` and the species model. Add a call for your solver there if its initial condition goes beyond the free-stream values set in its constructor.

## 8. Output

Add the output class that writes the screen, history and volume output. `CDriver::PreprocessOutput` calls `COutputFactory::CreateOutput(kindSolver, config, nDim)` (`SU2_CFD/src/output/COutputFactory.cpp`). Add a case:

```cpp
    case MAIN_SOLVER::HEAT_EQUATION:
      output = new CHeatOutput(config, nDim);
      break;
```

A value without a case gets a plain `COutput(config, nDim, false)`, which writes no solver fields.

- For a main solver, write an output class (copy `CHeatOutput`).
- For a sub-solver of the flow, add its fields to `CFlowOutput`.

[Adding Output Fields](../Adding-Output-Fields/) explains the hooks, the groups, and which groups go into the restart file.

## 9. Restart files

Implement `LoadRestart` so that the solver can read its solution back.

`CDriver::RestartSolver` calls `LoadRestart` for every non-null, non-adjoint solver when:

- `RESTART_SOL= YES`, or
- `config->GetRestart_Flow()` is true. Adjoint runs set it, because they read the primal solution.

The call runs inside an OpenMP region if `GetHasHybridParallel()` is true. A second loop loads the adjoint solvers.

The output class writes the restart file. The file holds the coordinates, then the volume fields of the restart groups, in the order of the `AddVolumeOutput` calls (see [Adding Output Fields](../Adding-Output-Fields/#groups-and-restart-files)). SU2 reads the file **by position**. So `LoadRestart` must skip the columns in front of its own variables. `CHeatSolver::LoadRestart` shows the pattern:

```cpp
  string restart_filename = config->GetSolution_FileName();

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
  /*--- Skip coordinates ---*/

  unsigned short skipVars = nDim;

  if (flow) {
    // P, vx, vy (,vz)
    skipVars += 1 + nDim;
  }
  // ...
  if (config->GetRead_Binary_Restart()) {
    restart_filename = config->GetFilename(restart_filename, ".dat", val_iter);
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    restart_filename = config->GetFilename(restart_filename, ".csv", val_iter);
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }
  // ... copy Restart_Data[counter*Restart_Vars[1] + skipVars + iVar] into nodes, check counter == nPointDomain
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
```

After reading the file, `CHeatSolver::LoadRestart`:

1. communicates the solution;
2. calls `Preprocessing`;
3. restricts the solution to the coarse levels with `MultigridRestriction`;
4. clears `Restart_Vars` and `Restart_Data`.

The file name comes from `SOLUTION_FILENAME`. The format comes from `READ_BINARY_RESTART`: SU2 binary (`.dat`) or ASCII (`.csv`). If the layout is simple, `CSolver::BasicLoadRestart` copies `nVar` values per point after `skipVars`.

Your solver may add fields to the restart file of another solver, for example a sub-solver that writes after the flow fields. Then the `LoadRestart` of the solvers that come after it must skip more columns.

## 10. MPI communication

Update the halo points after you change per-point data. Each rank stores halo points: copies of points that a neighbouring rank owns. After you change a per-point quantity that other ranks read, update the halos with this pair:

```cpp
  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);
```

`enum class MPI_QUANTITIES` (`Common/include/option_structure.hpp`) lists what can be sent (`SOLUTION`, `SOLUTION_OLD`, `SOLUTION_GRADIENT`, `SOLUTION_LIMITER`, `PRIMITIVE_GRADIENT`, `SOLUTION_EDDY`, ...). To send a new quantity:

1. Add a value to `MPI_QUANTITIES`.
2. Add a case to `CSolver::GetCommCountAndType` (`SU2_CFD/src/solvers/CSolver.cpp`) with the count per point. Without it, the `default` case stops with "Unrecognized quantity for point-to-point MPI comms."
3. Pack it in the `switch (commType)` of `CSolver::InitiateComms`, and unpack it in `CSolver::CompleteComms`.

Periodic boundaries use `InitiatePeriodicComms` / `CompletePeriodicComms`, with their own enum values (`PERIODIC_VOLUME`, `PERIODIC_NEIGHBORS`, ...). For any other communication, use the `SU2_MPI` wrapper (see [Code Structure](../Code-Structure/#parallelism-mpi-openmp-and-vectorization)).

## 11. OpenMP

Make every solver function thread safe before you turn on OpenMP.

`CSolver::GetHasHybridParallel()` returns `false` by default. The integration then opens its parallel region with a single thread, and the solver works as a pure MPI code. Start like this. Return `true` only when every function is thread safe. `CScalarSolver`, `CEulerSolver` and `CIncEulerSolver` return `true`.

When it is `true`, the **whole** iteration runs inside one parallel region (`SU2_OMP_PARALLEL_(if(...GetHasHybridParallel()))` in `SingleGrid_Iteration` and `MultiGrid_Iteration`). Every thread calls every solver function. So:

- Split loops with work-sharing macros, and close them. Use `SU2_OMP_FOR_STAT(omp_chunk_size)` for point loops and `SU2_OMP_FOR_DYN(OMP_MIN_SIZE)` for vertex loops. Close each one with `END_SU2_OMP_FOR`.
- Function-local variables are private to each thread; members are shared. Change shared state (members, `config`) inside `SU2_OMP_SAFE_GLOBAL_ACCESS(...)` or `BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS { ... } END_SU2_OMP_SAFE_GLOBAL_ACCESS`. `CHeatSolver::Preprocessing` does this for `config->SetGlobalParam`.
- Edge loops write to both points of an edge. There are two ways to avoid conflicts between threads:
  - Loop over `EdgeColoring`. Edges of one color share no point.
  - Support the reducer strategy (`ReducerStrategy`, `EdgeFluxes`, `SumEdgeFluxes`). The solver switches to it by itself when the coloring is poor (`parallelEff < COLORING_EFF_THRESH`). The edge loop then stores one flux per edge in `EdgeFluxes`. `SumEdgeFluxes` then adds them to the points.

  `CScalarSolver::EdgeFluxResidual` does both.

The macro list is in [Code Structure](../Code-Structure/#parallelism-mpi-openmp-and-vectorization).

## 12. Discrete adjoint

The discrete adjoint of a solver is [`CDiscAdjSolver`](../../doxygen/classCDiscAdjSolver.html), created on top of the primal solver:

```cpp
        heatSolver = new CDiscAdjSolver(geometry, config, solver[HEAT_SOL], RUNTIME_HEAT_SYS, iMGLevel);
```

The adjoint uses the AD tape: the record of the primal operations that automatic differentiation runs backwards. `CDiscAdjSolver` registers the primal solution on the tape:

- `direct_solver->GetNodes()->RegisterSolution(true)` for the inputs;
- `RegisterSolution(false)` for the outputs.

These call `AD::RegisterInput` / `AD::RegisterOutput` for every entry of `Solution`. Extra state that must be part of the fixed point goes through `RegisterSolutionExtra`, `ExtractAdjoint_SolutionExtra` and `SetAdjoint_SolutionExtra`.

To add the adjoint:

1. Add a `DISC_ADJ_` value of `MAIN_SOLVER` and its mapping in `CConfig::SetPostprocessing` (step 1).
2. Add an adjoint position (like `ADJHEAT_SOL`), and the `CDiscAdjSolver` in `CreateSolverContainer` with `INTEGRATION_TYPE::DEFAULT` (step 3).
3. Set up the recording.
   - Standalone solver: the constructor of `CDiscAdjSinglezoneDriver` picks `direct_iteration`, `direct_output`, `MainVariables`, `SecondaryVariables` and `MainSolver` for each `DISC_ADJ_` value. `CIterationFactory::CreateIteration` creates the adjoint iteration (`CDiscAdjHeatIteration` for `DISC_ADJ_HEAT`).
   - Sub-solver of the flow: follow `SPECIES_SOL` / `ADJSPECIES_SOL` through `CDiscAdjFluidIteration` (`Preprocess`, `LoadUnsteady_Solution`, `IterateDiscAdj`, `InitializeAdjoint`, `RegisterInput`, `SetDependencies`, `RegisterOutput`).
4. Add an adjoint output class (`CAdjHeatOutput` for heat) in `COutputFactory::CreateOutput`.

Follow these rules in the primal code, so that the tape is correct:

- Use `su2double` for everything that depends on the solution or on design variables. `passivedouble` and `SU2_TYPE::GetValue` cut the derivative.
- The Jacobian assembly is not recorded; the scalar kernels wrap it in `AD::BeginPassive` / `AD::EndPassive`. Only the residual matters for the adjoint.
- Kernels may use preaccumulation (`AD::StartPreacc`, `AD::SetPreaccIn`, `AD::SetPreaccOut`, `AD::EndPreacc`). Preaccumulation stores the local derivatives of a short code section on the tape, instead of all its operations. In edge loops with the reducer strategy, preaccumulation is paused with `AD::PausePreaccumulation`. With coloring, `AD::StartNoSharedReading` is used. See `CScalarSolver::EdgeFluxResidual`.

See [Advanced AD Techniques](../Advanced-AD-Techniques/).

## 13. Build and test

1. Add the `.cpp` files to the source list in `SU2_CFD/src/meson.build` (next to `'solvers/CHeatSolver.cpp'`, `'variables/CHeatVariable.cpp'`, ...).
2. Add the new options to `config_template.cfg` (see [Adding a Config Option](../Adding-a-Config-Option/)).
3. Do a dry run (`SU2_CFD -d case.cfg`). It builds the driver on a dummy geometry and lists the output fields.
4. Run serial, with MPI and, if `GetHasHybridParallel()` is true, with OpenMP (`SU2_CFD -t <threads>`). The results must agree.
5. Add a regression test (and an `_AD` one for the adjoint) as in [Running Regression Tests](../Running-Regression-Tests/#add-a-new-test).

## Non-conservative systems of equations

A term is conservative when an edge adds a flux to the row of one point and the opposite flux to the row of the other point. It is non-conservative when the two row contributions differ.

In short:

- **Scalar solvers already support non-conservative edge terms.** Every `CScalarSolver` child can have them. The SA and SST kernels use them today (see the [table below](#scalar-solvers-two-row-contributions-per-edge)).
- **To use this path**, derive your solver from `CScalarSolver` and write a kernel derived from `CUpwScalarBase`. In the kernel, make the two row contributions differ. The steps are below.
- **The flow solvers assume one conservative flux per edge.** Their edge loops and helpers are not written for two different row contributions. If your system does not fit the scalar model, write your own edge loop (see [Advice](#advice)).

### Scalar solvers: two row contributions per edge

The scalar edge kernels return the contributions to the two rows separately, in an `EdgeResidual` (`SU2_CFD/include/numerics/util.hpp`):

```cpp
/*!
 * \brief Residual of one edge, accumulated by the convective and the diffusive terms.
 * \note flux_i and flux_j are the contributions to the rows of i and j. They are opposite
 *       for a conservative term and independent for a non-conservative one. The Jacobians
 *       map onto the ii, ij, ji and jj blocks of the edge. A dynamic model sizes the storage
 *       with the maximum and iterates to nVar, which the scheme sets from the solver.
 */
template <class Double, size_t nVar_>
struct EdgeResidual {
  // ...
  Vector<Double, Size> flux_i, flux_j;
  Matrix<Double, Size, Size> jac_ii, jac_ij, jac_ji, jac_jj;
  const size_t nVar;
  // ...
};
```

The second `updateLinearSystem` overload in the same file writes it into the linear system. It handles both OpenMP strategies of [step 11](#11-openmp):

- **Coloring:** `flux_i` is added to the row of `iPoint`, and `flux_j` to the row of `jPoint`. `CSysMatrix::SetBlocks(iEdge, iPoint, jPoint, jac_ii, jac_ij, jac_ji, jac_jj)` adds the two diagonal blocks and **sets** the two off-diagonal blocks.
- **Reducer:** `flux_i` goes to `EdgeFluxes`, and `flux_j` to `EdgeFluxesDiff`. Only the off-diagonal blocks are set (`SetOffDiagBlocks`). `CScalarSolver::SumEdgeFluxes` then adds `EdgeFluxes` to the first point of each edge, and `EdgeFluxesDiff` to the second. `Jacobian.SetDiagonalAsColumnSum()` builds the diagonal.

  The note on that overload says this diagonal equals `jac_ii` and `jac_jj` only for a conservative flux. With a non-conservative model, the solver then converges along a slightly different path, to the same solution.

To use this path in a new scalar solver:

1. In its constructor, allocate both `EdgeFluxes` and `EdgeFluxesDiff` when `ReducerStrategy` is true. The `CScalarSolver` constructor does not do it. All existing scalar solvers do, for example `CHeatSolver`, `CTurbSASolver`, `CTurbSSTSolver`, `CTransLMSolver`, `CSpeciesSolver`, `CPoissonSolver`.
2. Call `EdgeFluxResidual<Scheme>` with a kernel derived from `CUpwScalarBase` (`SU2_CFD/include/numerics/scalar/scalar_edge_flux.hpp`), through `DispatchScheme` as in [step 6](#6-numerics).
3. In the kernel, make the two rows differ where the term is non-conservative.

How the base kernel fills the two rows:

- The convective term is conservative: `finalizeFlux` does `flux_i += flux` and `flux_j -= flux`.
- The diffusive term is written per row:

```cpp
        res.flux_i(iVar) -= D.i(iVar) * projGrad(iVar);
        // ...
          res.flux_j(iVar) += D.j(iVar) * projGrad(iVar);
```

`D` comes from the model's `coefficients()` function, which returns a pair `{D_i, D_j}`:

- if the two are equal, the diffusion is conservative;
- if they differ, it is not.

The Jacobian blocks follow the same split: `jac_ii` and `jac_ij` come from `D.i`, `jac_ji` and `jac_jj` from `D.j`.

The kernel has more hooks:

- `finalizeFlux`: your own convective term;
- `coefficientJacobians`: the derivatives of the coefficients;
- `extraDiffusionTerms`: extra fluxes (the flamelet kernel uses it).

Models that use this today:

| Kernel | `coefficients()` | Discretization |
|---|---|---|
| `CScalarFlux_SA` (`SU2_CFD/include/numerics/turbulent/turb_sa_edge_flux.hpp`) | `D_i` and `D_j` use the `nu_tilde` of the row's own point (the `cb2` term) | non-conservative |
| `CScalarFlux_SST` (`SU2_CFD/include/numerics/turbulent/turb_sst_edge_flux.hpp`) | the cross-diffusion part reads `omega` of the row's own point | non-conservative |
| `CScalarFlux_Heat` (`SU2_CFD/include/numerics/heat_edge_flux.hpp`), LM (`SU2_CFD/include/numerics/turbulent/transition/trans_edge_flux.hpp`), species (`SU2_CFD/include/numerics/species/species_edge_flux.hpp`) | `return {D, D};` | conservative |

The kernel flag `Conservative` means something else. It says whether the transported variable is weighted by the density. SST, LM and species set it to `true`; SA and heat set it to `false`.

### Flow solvers: one flux per edge

The flow solvers assemble each edge flux conservatively: the same flux is added to the row of `iPoint` and subtracted from the row of `jPoint`. In the non-vectorized loop of `CEulerSolver::Upwind_Residual`:

```cpp
    if (ReducerStrategy) {
      EdgeFluxes.SetBlock(iEdge, conv_residual);
    } else {
      LinSysRes.AddBlock(iPoint, conv_residual);
      LinSysRes.SubtractBlock(jPoint, conv_residual);
    }
```

The other parts follow the same rule:

- **Jacobian:** `CFVMFlowSolverBase::UpdateJacobian` updates it with `CSysMatrix::UpdateBlocks`, which takes two blocks. `block_i` is added to `ii` and subtracted from `ji`. `block_j` is added to `ij` and subtracted from `jj`. With the reducer strategy it calls `CSysMatrix::SetBlocks` instead. That writes only the off-diagonal blocks `ij` and `ji`; the diagonal blocks are computed afterwards from them (`SetDiagonalAsColumnSum`).
- **Vectorized path:** `CFVMFlowSolverBase::EdgeFluxResidual` processes several edges at once with SIMD instructions. It calls the first `updateLinearSystem` overload in `SU2_CFD/include/numerics/util.hpp` ("one conservative flux"). That overload uses `CSysVector::UpdateBlocks`: add to `iPoint`, subtract from `jPoint`.
- **Reducer strategy:** `CFVMFlowSolverBase::SumEdgeFluxes` adds `EdgeFluxes` for the first point of the edge and subtracts it for the second. `Jacobian.SetDiagonalAsColumnSum()` builds the diagonal blocks from the off-diagonal ones.

All of this assumes a conservative flux.

### Advice

- **Your system fits the scalar model** (transport of scalars, upwind convection, gradient diffusion, point sources). Derive the solver from `CScalarSolver`. Copy the SA kernel for asymmetric coefficients, or the heat kernel as the simplest one.
- **Your system does not fit.** Write your own edge loop that adds separate `i` and `j` contributions to `LinSysRes`, and the four Jacobian blocks.
  - Do not reuse the conservative helpers (`CSysVector::UpdateBlocks`, `CSysMatrix::UpdateBlocks`, the flow `SumEdgeFluxes`).
  - Write all four blocks of an edge in one call. In quantized Jacobian mode, the matrix stores the off-diagonal blocks in a compressed 8-bit form. These blocks can only be overwritten, not accumulated (see the note on `CFVMFlowSolverBase::UpdateJacobian`).
- **Check every part that relies on fluxes cancelling between neighbours:**
  - the reducer-strategy diagonal (above);
  - the multigrid forcing term: `CMultiGridIntegration::SetForcing_Term` sums the fine residuals of the children of a coarse cell (a coarse-grid cell made of several fine cells);
  - integrals or outputs built from boundary fluxes;
  - the periodic and MPI handling of edges and residuals.
- **A point source is an alternative.** You can also write a non-conservative term as a point source in `Source_Residual`, using the gradients stored at the point. SA used this form before. An older `CSourceBase_TurbSA` had a `ComputeCrossProduction` function for the `cb2` term. Commit `bcb0f94f5a` removed it. Today the term is in `CScalarFlux_SA::coefficients`.

## Checklist

1. `MAIN_SOLVER` value and `Solver_Map` keyword (or a sub-model enum and option); `DISC_ADJ_` mapping in `CConfig::SetPostprocessing`; problem helpers in `Common/include/CConfig.hpp`.
2. Container position (and `MAX_SOLS`), `RUNTIME_TYPE` value, `CConfig::GetContainerPosition`, `CConfig::SetGlobalParam`.
3. `SUB_SOLVER_TYPE` value, `CreateSolverContainer` and `CreateSubSolver` cases with the right `INTEGRATION_TYPE`.
4. `CVariable` child with `MAXNVAR` (for `CScalarSolver`).
5. `CSolver` child: constructor (linear system, `nodes`, `SetBaseClassPointerToNodes`, comms), residual loops, time step, BCs, update, `LoadRestart`.
6. Numerics in `CDriver::InitializeNumerics`, or a scalar edge kernel.
7. Iteration (`CIterationFactory::CreateIteration` or `CFluidIteration::CommonAuxiliarySolvers`), `SetInitialCondition` calls in the drivers.
8. Output class in `COutputFactory::CreateOutput`, restart columns in `LoadRestart`.
9. MPI and OpenMP correct; `GetHasHybridParallel()` only `true` when thread safe.
10. Adjoint position, recording and output, if needed.
11. `SU2_CFD/src/meson.build`, `config_template.cfg`, regression tests.

Last checked against SU2 develop `11c177b904`.
