---
title: Code Structure
permalink: /docs_v7/Code-Structure/
---

This page is a map of the SU2 source code for new developers. It shows where things live, how a simulation runs through the classes of `SU2_CFD`, and which files to change for common tasks.

For the full class documentation, see the [Code Reference](../Code-Reference/) (Doxygen). Class names on this page link to their Doxygen page.

---

- [Top-level folders](#top-level-folders)
- [Inside SU2_CFD](#inside-su2_cfd)
- [How SU2_CFD runs](#how-su2_cfd-runs)
- [Main classes](#main-classes)
- [Linear algebra](#linear-algebra)
- [Containers and toolboxes](#containers-and-toolboxes)
- [Parallelism: MPI, OpenMP and vectorization](#parallelism-mpi-openmp-and-vectorization)
- [Algorithmic differentiation](#algorithmic-differentiation)
- [Where do I change...](#where-do-i-change)

---

## Top-level folders

| Folder or file | What it holds |
|---|---|
| `Common` | Code shared by all programs: configuration (`CConfig`), mesh and dual grid (`CGeometry`), linear algebra, containers, MPI/OpenMP wrappers, AD types, grid movement, interface interpolation. Headers in `Common/include`, sources in `Common/src`. |
| `SU2_CFD` | The solver program `SU2_CFD`: drivers, iterations, integrations, solvers, variables, numerics, fluid models and output. |
| `SU2_DEF` | The mesh deformation program `SU2_DEF` (`SU2_DEF/src/SU2_DEF.cpp`). Its drivers `CDeformationDriver` and `CDiscAdjDeformationDriver` are in `SU2_DEF/include/drivers`. |
| `SU2_DOT` | The gradient projection program `SU2_DOT` (`SU2_DOT/src/SU2_DOT.cpp`). It uses `CDiscAdjDeformationDriver` from `SU2_DEF`. |
| `SU2_GEO` | The geometry evaluation program `SU2_GEO` (sections, thickness, volume, used as constraints in shape design). |
| `SU2_SOL` | The program `SU2_SOL`, which writes solution files from restart files (uses `CBaselineSolver` and `CBaselineOutput`). |
| `SU2_PY` | Python scripts (`parallel_computation.py`, `shape_optimization.py`, `discrete_adjoint.py`, ...), the `SU2` Python package (`SU2_PY/SU2`) and the SWIG interface of the Python wrapper (`SU2_PY/pySU2`). |
| `UnitTests` | Unit tests with Catch2 (`externals/catch2`), split into `UnitTests/Common` and `UnitTests/SU2_CFD`. Built with `-Denable-tests=true`. |
| `TestCases` | Regression tests. The test lists are `serial_regression.py`, `parallel_regression.py`, `hybrid_regression.py` and their `_AD` versions, plus `tutorials.py` and `vandv.py`. The shared runner is `TestCase.py`. |
| `externals` | Third-party code: CGNS, METIS, ParMETIS, TecIO, CLI11, Catch2, and the git submodules CoDiPack (`codi`), MeDiPack (`medi`), OpDiLib (`opdi`), MEL (`mel`), Eigen (`eigen`), FADO, meson and ninja. |
| `subprojects` | Optional libraries built with meson: CoolProp, Mutation++ (`Mutationpp`), MLPCpp. |
| `meson.build`, `meson_options.txt` | The build system. `meson_options.txt` lists the build options (`with-mpi`, `with-omp`, `enable-autodiff`, `enable-directdiff`, `enable-pywrapper`, `enable-tests`, ...). Each program folder has its own `meson.build` in `src`. |
| `meson.py`, `preconfigure.py`, `meson_scripts` | `meson.py` runs the meson version bundled in `externals/meson`. Before that, `preconfigure.py` and `meson_scripts/init.py` download the submodules and build ninja. |
| `SU2_IDE`, `docker`, `QuickStart` | IDE project files, container files, and the quick start case. |
| `Docs` | The Doxygen main page and the website pages that are kept with the code (`Docs/website`). |
| `config_template.cfg` | Every configuration option with a short comment. |

## Inside SU2_CFD

Headers are in `SU2_CFD/include/<folder>` and sources in `SU2_CFD/src/<folder>`. The folders `gradients` and `numerics_simd` only exist under `include`.

| Folder | Content |
|---|---|
| `drivers` | `CDriver` and its children. They create all objects and run the time loop. |
| `iteration` | `CIteration` and its children. One iteration of one physics (fluid, heat, structure, adjoint). |
| `integration` | `CIntegration` and its children. Multigrid, single grid, Newton-Krylov, structural and DG time integration. |
| `solvers` | `CSolver` and its children. Residuals, boundary conditions, linear system of one set of equations. Also `CSolverFactory`. |
| `variables` | `CVariable` and its children. The data stored at the mesh points. |
| `numerics` | `CNumerics` and its children. Fluxes and source terms for one edge or one point, in the subfolders `flow` (`convection/`, `flow_diffusion`, `flow_sources`), `turbulent`, `species`, `scalar`, `NEMO`, `elasticity`, `continuous_adjoint`. |
| `numerics_simd` | Vectorized flux schemes (`CNumericsSIMD`), used by `CEulerSolver`. |
| `fluid` | Fluid models (`CFluidModel`, `CIdealGas`, `CPengRobinson`, `CCoolProp`, ...), viscosity and conductivity models. |
| `gradients`, `limiters` | Green-Gauss and least squares gradients, slope limiters. |
| `interfaces` | Data transfer between zones (`CInterface`; `cfd`, `cht`, `fsi`). |
| `output` | `COutput` and its children, the file writers (`output/filewriter`) and the data sorters. See [Adding Output Fields](../Adding-Output-Fields/). |

`SU2_CFD/src/SU2_CFD.cpp` holds `main()`. `SU2_CFD/src/python_wrapper_structure.cpp` holds driver functions used by the Python wrapper.

## How SU2_CFD runs

```
main()                                         SU2_CFD/src/SU2_CFD.cpp
 |-- CConfig: read the config file (number of zones, solver kind)
 |-- create the driver:
 |     -d (dry run)            -> CDummyDriver
 |     single zone             -> CSinglezoneDriver / CDiscAdjSinglezoneDriver
 |     multizone               -> CMultizoneDriver  / CDiscAdjMultizoneDriver
 |     harmonic balance        -> CHBDriver
 |   The CDriver constructor creates everything:
 |     PreprocessInput, PreprocessOutput, InitializeGeometry, InitializeSolver,
 |     InitializeNumerics, InitializeIntegration, PreprocessIteration, ...
 |
 |-- driver->StartSolver()
 |     loop over time steps (CSinglezoneDriver):
 |       Preprocess -> Run -> Postprocess -> Update -> Monitor -> Output
 |       Run() calls CIteration::Solve()
 |         loop over inner iterations:
 |           CIteration::Iterate()                  e.g. CFluidIteration
 |             CIntegration::MultiGrid_Iteration()  flow (CMultiGridIntegration)
 |             CIntegration::SingleGrid_Iteration() turbulence, species, ... (CSingleGridIntegration)
 |               CSolver::Preprocessing()           primitives, gradients, limiters
 |               CSolver::SetTime_Step()
 |               CIntegration::Space_Integration()
 |                 CSolver::Upwind_Residual() or Centered_Residual()
 |                 CSolver::Viscous_Residual(), Source_Residual()
 |                   -> CNumerics::ComputeResidual() on each edge or point
 |                 CSolver::BC_...() on each marker
 |               CIntegration::Time_Integration()
 |                 CSolver::ImplicitEuler_Iteration() -> CSysSolve with CSysMatrix, CSysVector
 |               CSolver::Postprocessing()
 |           CIteration::Monitor()  -> COutput::SetHistoryOutput()
 |           CIteration::Output()   -> COutput::SetResultFiles()
 |
 |-- driver->Finalize()
```

Notes:

- `CMultiGridIntegration::MultiGrid_Iteration` runs the multigrid cycle. On each grid level it does the same steps as above (preprocessing, space and time integration, postprocessing).
- In `CFluidIteration::Iterate`, the flow is solved first. Then `CommonAuxiliarySolvers` runs `SingleGrid_Iteration` for the transition, turbulence and species solvers. The coupling is segregated.
- `CIntegration::Space_Integration` and `CIntegration::Time_Integration` are in `SU2_CFD/src/integration/CIntegration.cpp`. They show exactly which solver functions are called, and in which order.
- The discrete adjoint drivers (`CDiscAdjSinglezoneDriver`, `CDiscAdjMultizoneDriver`) record one primal iteration on the AD tape, then evaluate the tape many times (`AD::ComputeAdjoint`) to solve the adjoint fixed-point problem.

### The containers in the driver

The driver holds arrays of pointers. They are passed down to the iterations and integrations.

| Member | Indices |
|---|---|
| `config_container` | `[iZone]` |
| `geometry_container` | `[iZone][iInst][iMesh]` (iMesh = multigrid level, `MESH_0` is the finest) |
| `solver_container` | `[iZone][iInst][iMesh][iSol]` (iSol = `FLOW_SOL`, `TURB_SOL`, `TRANS_SOL`, `SPECIES_SOL`, `HEAT_SOL`, ...) |
| `numerics_container` | `[iZone][iInst][iMesh][iSol][iTerm]` (iTerm = `CONV_TERM`, `VISC_TERM`, `SOURCE_FIRST_TERM`, `CONV_BOUND_TERM`, ..., one set per OpenMP thread) |
| `integration_container` | `[iZone][iInst][iSol]` |
| `iteration_container` | `[iZone][iInst]` |
| `output_container` | `[iZone]` |

The solver indices and term indices are defined in `Common/include/option_structure.hpp`.

## Main classes

### Drivers

- [`CDriverBase`](../../doxygen/classCDriverBase.html): common base of the `SU2_CFD` and `SU2_DEF` drivers. It also holds the functions used by the Python wrapper.
- [`CDriver`](../../doxygen/classCDriver.html): creates geometry, solvers, numerics, integrations, iterations and outputs.
- [`CSinglezoneDriver`](../../doxygen/classCSinglezoneDriver.html), [`CMultizoneDriver`](../../doxygen/classCMultizoneDriver.html): primal problems. The multizone driver couples zones with block Gauss-Seidel or Jacobi iterations and transfers data through `CInterface` objects.
- [`CDiscAdjSinglezoneDriver`](../../doxygen/classCDiscAdjSinglezoneDriver.html), [`CDiscAdjMultizoneDriver`](../../doxygen/classCDiscAdjMultizoneDriver.html): discrete adjoint problems.
- [`CHBDriver`](../../doxygen/classCHBDriver.html) (child of [`CFluidDriver`](../../doxygen/classCFluidDriver.html)): harmonic balance.
- [`CDummyDriver`](../../doxygen/classCDummyDriver.html): dry run (`SU2_CFD -d`). It prints the available output fields.

### Iterations

[`CIteration`](../../doxygen/classCIteration.html) runs one iteration of one physics in one zone (`Preprocess`, `Iterate`, `Postprocess`, `Monitor`, `Output`, `Solve`). The children are created by `CIterationFactory`:

- [`CFluidIteration`](../../doxygen/classCFluidIteration.html), with children [`CTurboIteration`](../../doxygen/classCTurboIteration.html), [`CPBFluidIteration`](../../doxygen/classCPBFluidIteration.html) (pressure-based incompressible), [`CFEMFluidIteration`](../../doxygen/classCFEMFluidIteration.html), [`CHeatIteration`](../../doxygen/classCHeatIteration.html) and [`CAdjFluidIteration`](../../doxygen/classCAdjFluidIteration.html) (continuous adjoint).
- [`CFEAIteration`](../../doxygen/classCFEAIteration.html): structural analysis.
- [`CDiscAdjFluidIteration`](../../doxygen/classCDiscAdjFluidIteration.html), [`CDiscAdjHeatIteration`](../../doxygen/classCDiscAdjHeatIteration.html), [`CDiscAdjFEAIteration`](../../doxygen/classCDiscAdjFEAIteration.html): discrete adjoint.

### Integrations

[`CIntegration`](../../doxygen/classCIntegration.html) holds `Space_Integration` and `Time_Integration`. `CIntegrationFactory` picks the child from the solver type (`INTEGRATION_TYPE` in `SU2_CFD/include/solvers/CSolverFactory.hpp`):

- [`CMultiGridIntegration`](../../doxygen/classCMultiGridIntegration.html): finite volume flow solvers (also with `MGLEVEL= 0`).
- [`CSingleGridIntegration`](../../doxygen/classCSingleGridIntegration.html): turbulence, transition, species, heat, radiation solvers.
- [`CNewtonIntegration`](../../doxygen/classCNewtonIntegration.html): flow solvers with `NEWTON_KRYLOV= YES` (primal only).
- [`CStructuralIntegration`](../../doxygen/classCStructuralIntegration.html): structural (FEA) solver.
- [`CFEM_DG_Integration`](../../doxygen/classCFEM__DG__Integration.html): discontinuous Galerkin flow solver.
- The base `CIntegration` is used for the discrete adjoint solvers.

### Solvers

[`CSolver`](../../doxygen/classCSolver.html) solves one set of equations on one grid level. It computes residuals, applies boundary conditions (`BC_...` functions), and owns the linear system (`Jacobian`, `LinSysRes`, `LinSysSol`, `System`). `CSolverFactory` (`SU2_CFD/src/solvers/CSolverFactory.cpp`) decides which solvers are created for each `SOLVER` option.

Main families:

- Finite volume flow: [`CFVMFlowSolverBase`](../../doxygen/classCFVMFlowSolverBase.html) (template with the shared code) and its children [`CEulerSolver`](../../doxygen/classCEulerSolver.html) / [`CNSSolver`](../../doxygen/classCNSSolver.html) (compressible), [`CIncEulerSolver`](../../doxygen/classCIncEulerSolver.html) / [`CIncNSSolver`](../../doxygen/classCIncNSSolver.html) (incompressible), [`CNEMOEulerSolver`](../../doxygen/classCNEMOEulerSolver.html) / [`CNEMONSSolver`](../../doxygen/classCNEMONSSolver.html) (non-equilibrium).
- Scalar transport: [`CScalarSolver`](../../doxygen/classCScalarSolver.html) (template), with [`CTurbSolver`](../../doxygen/classCTurbSolver.html) → [`CTurbSASolver`](../../doxygen/classCTurbSASolver.html), [`CTurbSSTSolver`](../../doxygen/classCTurbSSTSolver.html), [`CTransLMSolver`](../../doxygen/classCTransLMSolver.html); [`CSpeciesSolver`](../../doxygen/classCSpeciesSolver.html) → [`CSpeciesFlameletSolver`](../../doxygen/classCSpeciesFlameletSolver.html); [`CHeatSolver`](../../doxygen/classCHeatSolver.html); [`CPoissonSolver`](../../doxygen/classCPoissonSolver.html).
- Others: [`CRadP1Solver`](../../doxygen/classCRadP1Solver.html), [`CFEASolver`](../../doxygen/classCFEASolver.html), [`CMeshSolver`](../../doxygen/classCMeshSolver.html) (mesh deformation), [`CFEM_DG_EulerSolver`](../../doxygen/classCFEM__DG__EulerSolver.html) / [`CFEM_DG_NSSolver`](../../doxygen/classCFEM__DG__NSSolver.html), [`CBaselineSolver`](../../doxygen/classCBaselineSolver.html).
- Adjoint: [`CDiscAdjSolver`](../../doxygen/classCDiscAdjSolver.html) (discrete adjoint of any solver), [`CDiscAdjFEASolver`](../../doxygen/classCDiscAdjFEASolver.html), [`CDiscAdjMeshSolver`](../../doxygen/classCDiscAdjMeshSolver.html), and the continuous adjoint [`CAdjEulerSolver`](../../doxygen/classCAdjEulerSolver.html), [`CAdjNSSolver`](../../doxygen/classCAdjNSSolver.html), [`CAdjTurbSolver`](../../doxygen/classCAdjTurbSolver.html).
- [`CTemplateSolver`](../../doxygen/classCTemplateSolver.html): an empty starting point for a new solver.

### Variables

[`CVariable`](../../doxygen/classCVariable.html) stores the data of **all** points of one solver in one object (`solver->GetNodes()`). Each quantity is its own container indexed by point: `Solution`, `Solution_Old`, `Gradient`, `Limiter`, `Delta_Time`, ... Matrices are row-major (`su2activematrix`), so the values of one point are next to each other in memory. Access is by point index, for example `nodes->GetSolution(iPoint, iVar)`.

Children follow the solvers: [`CFlowVariable`](../../doxygen/classCFlowVariable.html) → [`CEulerVariable`](../../doxygen/classCEulerVariable.html) → [`CNSVariable`](../../doxygen/classCNSVariable.html); [`CIncEulerVariable`](../../doxygen/classCIncEulerVariable.html) → [`CIncNSVariable`](../../doxygen/classCIncNSVariable.html); [`CScalarVariable`](../../doxygen/classCScalarVariable.html) → [`CTurbVariable`](../../doxygen/classCTurbVariable.html) → [`CTurbSAVariable`](../../doxygen/classCTurbSAVariable.html), [`CTurbSSTVariable`](../../doxygen/classCTurbSSTVariable.html); and so on.

### Numerics

[`CNumerics`](../../doxygen/classCNumerics.html) computes a flux across one edge or a source term at one point. `ComputeResidual(config)` returns a `ResidualType` with the residual and the Jacobians for the two points of the edge. The solver loops over the edges or points, sets the inputs of the numerics object, calls `ComputeResidual`, and adds the result to `LinSysRes` and `Jacobian`.

The numerics objects are created in `CDriver::InitializeNumerics` (`SU2_CFD/src/drivers/CDriver.cpp`), with the helpers `InstantiateTurbulentNumerics`, `InstantiateTransitionNumerics` and `InstantiateSpeciesNumerics`.

### Geometry

[`CGeometry`](../../doxygen/classCGeometry.html) holds the mesh and the dual grid used by the edge-based finite volume method:

- `nodes` ([`CPoint`](../../doxygen/classCPoint.html)): coordinates, volumes, neighbors, of all points.
- `edges` ([`CEdge`](../../doxygen/classCEdge.html)): the two points and the normal of each edge.
- `vertex[iMarker][iVertex]` ([`CVertex`](../../doxygen/classCVertex.html)): boundary points and their normals.
- `InitiateComms` / `CompleteComms`: halo exchange between MPI ranks.

Children:

- [`CPhysicalGeometry`](../../doxygen/classCPhysicalGeometry.html): the mesh read from file (SU2 or CGNS format, readers in `Common/include/geometry/meshreader`), partitioned with ParMETIS.
- [`CMultiGridGeometry`](../../doxygen/classCMultiGridGeometry.html): coarse levels made by agglomeration of control volumes.
- [`CDummyGeometry`](../../doxygen/classCDummyGeometry.html): empty geometry for the dry run.
- [`CMeshFEM_DG`](../../doxygen/classCMeshFEM__DG.html): geometry of the DG solver (`Common/include/fem`).

Grid movement (`CSurfaceMovement`, `CVolumetricMovement`, `CFreeFormDefBox`) is in `Common/include/grid_movement`. Interpolation between non-matching zone interfaces (`CInterpolator` and children) is in `Common/include/interface_interpolation`.

### Configuration

[`CConfig`](../../doxygen/classCConfig.html) (`Common/include/CConfig.hpp`, `Common/src/CConfig.cpp`) reads the configuration file and stores every option. There is one `CConfig` per zone, plus a `driver_config` read from the main configuration file (it holds the settings of the driver, for example for multizone problems).

- Options are declared in `CConfig::SetConfig_Options` with `addDoubleOption`, `addBoolOption`, `addEnumOption`, `addStringListOption`, ...
- Enums and the maps from option strings to enum values are in `Common/include/option_structure.hpp`. The option parser classes are in `Common/include/option_structure.inl`.
- Derived settings are computed in `CConfig::SetPostprocessing` and `CConfig::SetMarkers`.

See [Adding a Config Option](../Adding-a-Config-Option/) and the [Option Reference](../Option-Reference/).

### Output

[`COutput`](../../doxygen/classCOutput.html) writes the screen output, the history file and the volume and surface files. Each solver kind has a child class, created by `COutputFactory`, for example [`CFlowCompOutput`](../../doxygen/classCFlowCompOutput.html), [`CFlowIncOutput`](../../doxygen/classCFlowIncOutput.html), [`CHeatOutput`](../../doxygen/classCHeatOutput.html), [`CElasticityOutput`](../../doxygen/classCElasticityOutput.html). [`CMultizoneOutput`](../../doxygen/classCMultizoneOutput.html) combines the zones of a multizone run. File writers (Paraview, Tecplot, CGNS, CSV, SU2 restart, ...) are in `SU2_CFD/include/output/filewriter`.

See [Adding Output Fields](../Adding-Output-Fields/) for the details.

## Linear algebra

Files: `Common/include/linear_algebra` and `Common/src/linear_algebra`.

- [`CSysMatrix`](../../doxygen/classCSysMatrix.html): block sparse matrix (one block per pair of neighbor points). Used for the Jacobian of implicit solvers.
- [`CSysVector`](../../doxygen/classCSysVector.html): block vector, for residual and solution update.
- [`CSysSolve`](../../doxygen/classCSysSolve.html): Krylov solvers (`FGMRES_LinSolver`, `RFGMRES_LinSolver`, `FGCRODR_LinSolver`, `BCGSTAB_LinSolver`, `CG_LinSolver`, `Smoother_LinSolver`). `Solve` reads the method and the preconditioner from the configuration (for the primal flow: `LINEAR_SOLVER` and `LINEAR_SOLVER_PREC`; the mesh, adjoint and other problems have their own options).
- [`CPreconditioner`](../../doxygen/classCPreconditioner.html) and its children [`CJacobiPreconditioner`](../../doxygen/classCJacobiPreconditioner.html), [`CILUPreconditioner`](../../doxygen/classCILUPreconditioner.html), [`CLU_SGSPreconditioner`](../../doxygen/classCLU__SGSPreconditioner.html), [`CLineletPreconditioner`](../../doxygen/classCLineletPreconditioner.html), [`CPastixPreconditioner`](../../doxygen/classCPastixPreconditioner.html). `CPreconditioner::Create` builds the one that was chosen.
- [`CMatrixVectorProduct`](../../doxygen/classCMatrixVectorProduct.html): the product used inside the Krylov solvers.

Each solver has the members `Jacobian` (`CSysMatrix<su2mixedfloat>`), `LinSysRes`, `LinSysSol` and `System` (`CSysSolve<su2mixedfloat>`). `su2mixedfloat` is `float` when SU2 is built with `-Denable-mixedprec=true`. See also [Linear Solvers and Preconditioners](../Linear-Solvers-and-Preconditioners/).

## Containers and toolboxes

Containers (`Common/include/containers`):

- [`C2DContainer`](../../doxygen/classC2DContainer.html) (`C2DContainer.hpp`): aligned 2D storage. The common type aliases are defined in the same file:
  - `su2vector<T>`, `su2matrix<T>` (row-major), `ColMajorMatrix<T>`;
  - `su2activevector`, `su2activematrix` (of `su2double`);
  - `su2passivevector`, `su2passivematrix` (of `passivedouble`).
- `container_decorators.hpp`: `C3DDoubleMatrix` (also called `CVectorOfMatrix`, used for gradients) and `CMatrixView`.
- `CLookUpTable`, `CFileReaderLUT`, `CTrapezoidalMap`: look-up tables, used by the flamelet fluid model (`CFluidFlamelet`) and the data-driven fluid model (`CDataDrivenFluid`).
- `CVertexMap`, `CFastFindAndEraseQueue`, `CPyWrapperMatrixView`.

Toolboxes (`Common/include/toolboxes`): `geometry_toolbox.hpp` (vector operations), `graph_toolbox.hpp` (sparse patterns, coloring), `printing_toolbox.hpp` (`PrintingToolbox::CTablePrinter`), `ndflattener.hpp`, `CSquareMatrixCM`, `CSymmetricMatrix`, `C1DInterpolation`, `CLinearPartitioner`, `CQuasiNewtonInvLeastSquares`, and the verification solutions in `toolboxes/MMS`.

## Parallelism: MPI, OpenMP and vectorization

Files: `Common/include/parallelization`.

**MPI.** Use the `SU2_MPI` wrapper from `mpi_structure.hpp`, never MPI directly. `SU2_MPI` is `CBaseMPIWrapper`, or `CMediMPIWrapper` (MeDiPack) in AD builds, so that MPI calls work with `su2double`. The mesh is partitioned at start-up. Each rank also stores halo points. Halo values are updated with `InitiateComms` / `CompleteComms` of `CGeometry` and `CSolver`, with an `MPI_QUANTITIES` value that says what to send.

**OpenMP.** Do not include `omp.h`; include `omp_structure.hpp` and use its macros. They compile to nothing when OpenMP is off (`-Dwith-omp=false`), and to OpDiLib calls in AD builds:

- `SU2_OMP_PARALLEL` ... `END_SU2_OMP_PARALLEL`, `SU2_OMP_PARALLEL_(ARGS)`
- `SU2_OMP_FOR_STAT(CHUNK)`, `SU2_OMP_FOR_DYN(CHUNK)`, `SU2_OMP_FOR_(ARGS)` ... `END_SU2_OMP_FOR`
- `SU2_OMP_MASTER` ... `END_SU2_OMP_MASTER`, `SU2_OMP_CRITICAL` ... `END_SU2_OMP_CRITICAL`, `SU2_OMP_BARRIER`, `SU2_OMP_ATOMIC`
- `SU2_OMP_SAFE_GLOBAL_ACCESS(...)`: run code on one thread with barriers around it
- `SU2_OMP_SIMD`
- helpers: `omp_get_num_threads()`, `omp_get_thread_num()`, `roundUpDiv`, `atomicAdd`

The integrations open one parallel region for a whole iteration, when the solver supports it (`CSolver::GetHasHybridParallel()`). Edge loops write to both points of an edge. To avoid race conditions, edges are grouped by color (`CGeometry::GetEdgeColoring`, option `EDGE_COLORING_GROUP_SIZE`). If the coloring is poor, the flow and scalar solvers store the edge fluxes first and sum them per point afterwards (`ReducerStrategy`, `EdgeFluxes`). The number of threads per rank is set with `SU2_CFD -t <threads>`.

**Vectorization.** `vectorization.hpp` defines [`simd::Array`](../../doxygen/classsimd_1_1Array.html), a short fixed-size array with arithmetic operators. The compiler turns its loops into SIMD instructions, and `special_vectorization.hpp` has explicit versions for some instruction sets. [`CNumericsSIMD`](../../doxygen/classCNumericsSIMD.html) (`SU2_CFD/include/numerics_simd`) uses it to compute the fluxes of several edges at once. `CEulerSolver` uses these schemes for all centered schemes, and for the `ROE` and `MSW` upwind schemes with an ideal gas and no low-Mach correction (see `CEulerSolver::Centered_Residual` and `CEulerSolver::Upwind_Residual`). The build option `cpu-arch` (default `native`) sets the target instruction set.

## Algorithmic differentiation

SU2 uses operator overloading AD with [CoDiPack](https://github.com/scicompkl/CoDiPack) (`externals/codi`).

- `su2double` (`Common/include/code_config.hpp`) is the type for all values that may need derivatives. It is `passivedouble` (normally `double`) in a normal build, a CoDiPack reverse type with `-Denable-autodiff=true` (builds `SU2_CFD_AD`, `SU2_DOT_AD`), and a forward type with `-Denable-directdiff=true` (builds `SU2_CFD_DIRECTDIFF`).
- `passivedouble` is never differentiated. Use `SU2_TYPE::GetValue(x)` to get the plain value of a `su2double`.
- The [`AD`](../../doxygen/namespaceAD.html) namespace (`Common/include/basic_types/ad_structure.hpp`) wraps the tape: `StartRecording`, `StopRecording`, `RegisterInput`, `RegisterOutput`, `ComputeAdjoint`, `ClearAdjoints`, and preaccumulation (`StartPreacc`, `SetPreaccIn`, `SetPreaccOut`, `EndPreacc`). The functions do nothing in non-AD builds, so the code can call them always.
- MeDiPack (`externals/medi`) differentiates MPI calls. OpDiLib (`externals/opdi`) differentiates OpenMP regions.

See [Advanced AD Techniques](../Advanced-AD-Techniques/) for more.

## Where do I change...

| Task | Where |
|---|---|
| **New boundary condition** | 1. Add the kind to `enum BC_TYPE` in `Common/include/option_structure.hpp`. 2. Add the marker option in `CConfig::SetConfig_Options` and set `Marker_CfgFile_KindBC` in `CConfig::SetMarkers` (`Common/src/CConfig.cpp`). 3. Add a virtual `BC_...` function to `CSolver` (`SU2_CFD/include/solvers/CSolver.hpp`). 4. Implement it in the solvers that need it (for example `CEulerSolver`, `CIncEulerSolver`, and the scalar solvers through `CScalarSolver`). 5. Call it from the `switch` on `GetMarker_All_KindBC` in `CIntegration::Space_Integration`. |
| **New convective scheme** | 1. Add the value to `enum class UPWIND` (or `CENTERED`) and to `Upwind_Map` in `option_structure.hpp`. 2. Write the numerics class in `SU2_CFD/include/numerics/flow/convection/` and `SU2_CFD/src/numerics/flow/convection/`. 3. Create it in `CDriver::InitializeNumerics` (switch on `GetKind_Upwind_Flow()`). 4. Optional, compressible ideal gas only: a vectorized version in `SU2_CFD/include/numerics_simd/flow/convection/upwind.hpp`, created in `createUpwindIdealNumerics` (`CNumericsSIMD.cpp`) and selected in `CEulerSolver::Upwind_Residual`. |
| **New turbulence model source term** | Source terms are in `SU2_CFD/include/numerics/turbulent/turb_sources.hpp` (`CSourceBase_TurbSA` and `SAFactory` for SA, `CSourcePieceWise_TurbSST` for SST). They are created in `CDriver::InstantiateTurbulentNumerics` and used in `CTurbSASolver::Source_Residual` / `CTurbSSTSolver::Source_Residual`. Model variants are chosen with `SA_OPTIONS` / `SST_OPTIONS` (`enum class SA_OPTIONS`, `SST_OPTIONS` in `option_structure.hpp`). |
| **New output field** | The output class of the solver in `SU2_CFD/src/output/`, or `CFlowOutput` for turbulence and species fields. See [Adding Output Fields](../Adding-Output-Fields/). |
| **New config option** | `CConfig::SetConfig_Options`, a member and a getter in `CConfig.hpp`, an entry in `config_template.cfg`. See [Adding a Config Option](../Adding-a-Config-Option/). |
| **New solver** | Start from `CTemplateSolver` (`SU2_CFD/include/solvers/CTemplateSolver.hpp`) and the template numerics (`CConvective_Template`, `CViscous_Template`, `CSource_Template` in `SU2_CFD/include/numerics/template.hpp`). Then: a value in `enum class MAIN_SOLVER` and `Solver_Map` (`option_structure.hpp`); the solver in `CSolverFactory::CreateSolverContainer` and `CreateSubSolver` (with its `INTEGRATION_TYPE`); its numerics in `CDriver::InitializeNumerics`; an iteration in `CIterationFactory::CreateIteration`; an output class in `COutputFactory::CreateOutput`. |

When you change code, add or update a regression test (see [Running Regression Tests](../Running-Regression-Tests/)) or a unit test (see [Writing Unit Tests](../Writing-Unit-Tests/)), and follow the [Style Guide](../Style-Guide/).

Last checked against SU2 develop `11c177b904`.
