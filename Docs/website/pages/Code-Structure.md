---
title: Code Structure
permalink: /docs_v7/Code-Structure/
---

This page is a map of the SU2 source code for new developers. It shows:

- where things live;
- how a simulation runs through the classes of `SU2_CFD`;
- which files to change for common tasks.

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
| `Common` | Code shared by all programs: configuration (`CConfig`), mesh and dual grid (`CGeometry`), linear algebra, containers, MPI/OpenMP wrappers, AD types, grid movement, interface interpolation. Headers are in `Common/include`, sources in `Common/src`. |
| `SU2_CFD` | The solver program `SU2_CFD`: drivers, iterations, integrations, solvers, variables, numerics, fluid models and output. |
| `SU2_DEF` | The mesh deformation program `SU2_DEF` (`SU2_DEF/src/SU2_DEF.cpp`). Its drivers `CDeformationDriver` and `CDiscAdjDeformationDriver` are in `SU2_DEF/include/drivers`. |
| `SU2_DOT` | The gradient projection program `SU2_DOT` (`SU2_DOT/src/SU2_DOT.cpp`). It uses `CDiscAdjDeformationDriver` from `SU2_DEF`. |
| `SU2_GEO` | The geometry evaluation program `SU2_GEO`. It computes sections, thickness and volume, which shape design uses as constraints. |
| `SU2_SOL` | The program `SU2_SOL`. It writes solution files from restart files, with `CBaselineSolver` and `CBaselineOutput`. |
| `SU2_PY` | Python scripts (`parallel_computation.py`, `shape_optimization.py`, `discrete_adjoint.py`, ...), the `SU2` Python package (`SU2_PY/SU2`) and the SWIG interface of the Python wrapper (`SU2_PY/pySU2`). |
| `UnitTests` | Unit tests with Catch2 (`externals/catch2`), split into `UnitTests/Common` and `UnitTests/SU2_CFD`. Build them with `-Denable-tests=true`. |
| `TestCases` | Regression tests. The test lists are `serial_regression.py`, `parallel_regression.py`, `hybrid_regression.py` and their `_AD` versions, plus `tutorials.py` and `vandv.py`. All of them use the shared runner `TestCase.py`. |
| `externals` | Third-party code: CGNS, METIS, ParMETIS, TecIO, CLI11, Catch2. Also the git submodules CoDiPack (`codi`), MeDiPack (`medi`), OpDiLib (`opdi`), MEL (`mel`), Eigen (`eigen`), FADO, meson and ninja. |
| `subprojects` | Optional libraries built with meson: CoolProp, Mutation++ (`Mutationpp`), MLPCpp. |
| `meson.build`, `meson_options.txt` | The build system. `meson_options.txt` lists the build options (`with-mpi`, `with-omp`, `enable-autodiff`, `enable-directdiff`, `enable-pywrapper`, `enable-tests`, ...). Each program folder has its own `meson.build` in `src`. |
| `meson.py`, `preconfigure.py`, `meson_scripts` | `meson.py` runs the meson version bundled in `externals/meson`. Before that, `preconfigure.py` and `meson_scripts/init.py` download the submodules and build ninja. |
| `SU2_IDE`, `docker`, `QuickStart` | IDE project files, container files, and the quick start case. |
| `Docs` | The Doxygen main page, and the website pages that are kept with the code (`Docs/website`). |
| `config_template.cfg` | Every configuration option, with a short comment. |

## Inside SU2_CFD

Headers are in `SU2_CFD/include/<folder>` and sources are in `SU2_CFD/src/<folder>`. The folders `gradients` and `numerics_simd` exist only under `include`.

| Folder | Content |
|---|---|
| `drivers` | `CDriver` and its children. They create all objects and run the time loop. |
| `iteration` | `CIteration` and its children. Each one runs one iteration of one physics (fluid, heat, structure, adjoint). |
| `integration` | `CIntegration` and its children: multigrid, single grid, Newton-Krylov, structural and DG time integration. |
| `solvers` | `CSolver` and its children. A solver computes the residuals, applies the boundary conditions and holds the linear system of one set of equations. Also `CSolverFactory`. |
| `variables` | `CVariable` and its children. They store the data at the mesh points. |
| `numerics` | `CNumerics` and its children. They compute the fluxes for one edge or the source terms for one point. Subfolders: `flow` (`convection/`, `flow_diffusion`, `flow_sources`), `turbulent`, `species`, `scalar`, `NEMO`, `elasticity`, `continuous_adjoint`. |
| `numerics_simd` | Vectorized flux schemes (`CNumericsSIMD`), used by `CEulerSolver`. |
| `fluid` | Fluid models (`CFluidModel`, `CIdealGas`, `CPengRobinson`, `CCoolProp`, ...), viscosity and conductivity models. |
| `gradients`, `limiters` | Green-Gauss and least squares gradients, slope limiters. |
| `interfaces` | Data transfer between zones (`CInterface`; `cfd`, `cht`, `fsi`). |
| `output` | `COutput` and its children, the file writers (`output/filewriter`) and the data sorters. See [Adding Output Fields](../Adding-Output-Fields/). |

`SU2_CFD/src/SU2_CFD.cpp` holds `main()`. `SU2_CFD/src/python_wrapper_structure.cpp` holds the driver functions that the Python wrapper uses.

## How SU2_CFD runs

The tree below shows the calls of one `SU2_CFD` run, from `main()` down to the numerics.

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

- `CMultiGridIntegration::MultiGrid_Iteration` runs the multigrid cycle. On each grid level it does the same steps as above: preprocessing, space and time integration, postprocessing.
- `CFluidIteration::Iterate` solves the flow first. Then `CommonAuxiliarySolvers` runs `SingleGrid_Iteration` for the transition, turbulence and species solvers. The coupling is segregated: each solver is solved on its own, one after the other.
- `CIntegration::Space_Integration` and `CIntegration::Time_Integration` are in `SU2_CFD/src/integration/CIntegration.cpp`. Read them to see exactly which solver functions are called, and in which order.
- The discrete adjoint drivers (`CDiscAdjSinglezoneDriver`, `CDiscAdjMultizoneDriver`) record one primal iteration on the AD tape. The tape is the stored list of operations that AD later runs backwards to get derivatives. The drivers then evaluate the tape many times (`AD::ComputeAdjoint`) to solve the adjoint fixed-point problem.

### The containers in the driver

The driver holds arrays of pointers to all objects. It passes them down to the iterations and integrations.

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

- [`CDriverBase`](../../doxygen/classCDriverBase.html): common base of the `SU2_CFD` and `SU2_DEF` drivers. It also holds the functions that the Python wrapper uses.
- [`CDriver`](../../doxygen/classCDriver.html): creates geometry, solvers, numerics, integrations, iterations and outputs.
- [`CSinglezoneDriver`](../../doxygen/classCSinglezoneDriver.html), [`CMultizoneDriver`](../../doxygen/classCMultizoneDriver.html): primal problems. The multizone driver couples the zones with block Gauss-Seidel or Jacobi iterations. It transfers data between zones through `CInterface` objects.
- [`CDiscAdjSinglezoneDriver`](../../doxygen/classCDiscAdjSinglezoneDriver.html), [`CDiscAdjMultizoneDriver`](../../doxygen/classCDiscAdjMultizoneDriver.html): discrete adjoint problems.
- [`CHBDriver`](../../doxygen/classCHBDriver.html) (child of [`CFluidDriver`](../../doxygen/classCFluidDriver.html)): harmonic balance.
- [`CDummyDriver`](../../doxygen/classCDummyDriver.html): dry run (`SU2_CFD -d`). It prints the available output fields.

### Iterations

[`CIteration`](../../doxygen/classCIteration.html) runs one iteration of one physics in one zone (`Preprocess`, `Iterate`, `Postprocess`, `Monitor`, `Output`, `Solve`). `CIterationFactory` creates the children:

- [`CFluidIteration`](../../doxygen/classCFluidIteration.html), with the children [`CTurboIteration`](../../doxygen/classCTurboIteration.html), [`CPBFluidIteration`](../../doxygen/classCPBFluidIteration.html) (pressure-based incompressible), [`CFEMFluidIteration`](../../doxygen/classCFEMFluidIteration.html), [`CHeatIteration`](../../doxygen/classCHeatIteration.html) and [`CAdjFluidIteration`](../../doxygen/classCAdjFluidIteration.html) (continuous adjoint).
- [`CFEAIteration`](../../doxygen/classCFEAIteration.html): structural analysis.
- [`CDiscAdjFluidIteration`](../../doxygen/classCDiscAdjFluidIteration.html), [`CDiscAdjHeatIteration`](../../doxygen/classCDiscAdjHeatIteration.html), [`CDiscAdjFEAIteration`](../../doxygen/classCDiscAdjFEAIteration.html): discrete adjoint.

### Integrations

[`CIntegration`](../../doxygen/classCIntegration.html) holds `Space_Integration` and `Time_Integration`. `CIntegrationFactory` picks the child from the solver type (`INTEGRATION_TYPE` in `SU2_CFD/include/solvers/CSolverFactory.hpp`):

- [`CMultiGridIntegration`](../../doxygen/classCMultiGridIntegration.html): finite volume flow solvers (also with `MGLEVEL= 0`).
- [`CSingleGridIntegration`](../../doxygen/classCSingleGridIntegration.html): turbulence, transition, species, heat and radiation solvers.
- [`CNewtonIntegration`](../../doxygen/classCNewtonIntegration.html): flow solvers with `NEWTON_KRYLOV= YES` (primal only).
- [`CStructuralIntegration`](../../doxygen/classCStructuralIntegration.html): structural (FEA) solver.
- [`CFEM_DG_Integration`](../../doxygen/classCFEM__DG__Integration.html): discontinuous Galerkin flow solver.
- The discrete adjoint solvers use the base `CIntegration`.

### Solvers

[`CSolver`](../../doxygen/classCSolver.html) solves one set of equations on one grid level. It:

- computes the residuals;
- applies the boundary conditions (`BC_...` functions);
- owns the linear system (`Jacobian`, `LinSysRes`, `LinSysSol`, `System`).

`CSolverFactory` (`SU2_CFD/src/solvers/CSolverFactory.cpp`) decides which solvers to create for each `SOLVER` option.

Main families:

- Finite volume flow: [`CFVMFlowSolverBase`](../../doxygen/classCFVMFlowSolverBase.html) is a template with the shared code. Its children are [`CEulerSolver`](../../doxygen/classCEulerSolver.html) / [`CNSSolver`](../../doxygen/classCNSSolver.html) (compressible), [`CIncEulerSolver`](../../doxygen/classCIncEulerSolver.html) / [`CIncNSSolver`](../../doxygen/classCIncNSSolver.html) (incompressible) and [`CNEMOEulerSolver`](../../doxygen/classCNEMOEulerSolver.html) / [`CNEMONSSolver`](../../doxygen/classCNEMONSSolver.html) (non-equilibrium).
- Scalar transport: [`CScalarSolver`](../../doxygen/classCScalarSolver.html) (template), with:
  - [`CTurbSolver`](../../doxygen/classCTurbSolver.html) → [`CTurbSASolver`](../../doxygen/classCTurbSASolver.html), [`CTurbSSTSolver`](../../doxygen/classCTurbSSTSolver.html), [`CTransLMSolver`](../../doxygen/classCTransLMSolver.html);
  - [`CSpeciesSolver`](../../doxygen/classCSpeciesSolver.html) → [`CSpeciesFlameletSolver`](../../doxygen/classCSpeciesFlameletSolver.html);
  - [`CHeatSolver`](../../doxygen/classCHeatSolver.html);
  - [`CPoissonSolver`](../../doxygen/classCPoissonSolver.html).
- Others: [`CRadP1Solver`](../../doxygen/classCRadP1Solver.html), [`CFEASolver`](../../doxygen/classCFEASolver.html), [`CMeshSolver`](../../doxygen/classCMeshSolver.html) (mesh deformation), [`CFEM_DG_EulerSolver`](../../doxygen/classCFEM__DG__EulerSolver.html) / [`CFEM_DG_NSSolver`](../../doxygen/classCFEM__DG__NSSolver.html), [`CBaselineSolver`](../../doxygen/classCBaselineSolver.html).
- Adjoint: [`CDiscAdjSolver`](../../doxygen/classCDiscAdjSolver.html) (discrete adjoint of any solver), [`CDiscAdjFEASolver`](../../doxygen/classCDiscAdjFEASolver.html), [`CDiscAdjMeshSolver`](../../doxygen/classCDiscAdjMeshSolver.html), and the continuous adjoint solvers [`CAdjEulerSolver`](../../doxygen/classCAdjEulerSolver.html), [`CAdjNSSolver`](../../doxygen/classCAdjNSSolver.html), [`CAdjTurbSolver`](../../doxygen/classCAdjTurbSolver.html).
- [`CTemplateSolver`](../../doxygen/classCTemplateSolver.html): an empty starting point for a new solver.

### Variables

[`CVariable`](../../doxygen/classCVariable.html) stores the data of **all** points of one solver in one object (`solver->GetNodes()`). Each quantity is its own container, indexed by point: `Solution`, `Solution_Old`, `Gradient`, `Limiter`, `Delta_Time`, ... The matrices are row-major (`su2activematrix`), so the values of one point sit next to each other in memory. You access a value by point index, for example `nodes->GetSolution(iPoint, iVar)`.

The children follow the solvers:

- [`CFlowVariable`](../../doxygen/classCFlowVariable.html) → [`CEulerVariable`](../../doxygen/classCEulerVariable.html) → [`CNSVariable`](../../doxygen/classCNSVariable.html);
- [`CIncEulerVariable`](../../doxygen/classCIncEulerVariable.html) → [`CIncNSVariable`](../../doxygen/classCIncNSVariable.html);
- [`CScalarVariable`](../../doxygen/classCScalarVariable.html) → [`CTurbVariable`](../../doxygen/classCTurbVariable.html) → [`CTurbSAVariable`](../../doxygen/classCTurbSAVariable.html), [`CTurbSSTVariable`](../../doxygen/classCTurbSSTVariable.html);
- and so on.

### Numerics

[`CNumerics`](../../doxygen/classCNumerics.html) computes a flux across one edge, or a source term at one point. `ComputeResidual(config)` returns a `ResidualType` with the residual and the Jacobians for the two points of the edge.

The solver uses it in this way:

1. It loops over the edges or points.
2. It sets the inputs of the numerics object.
3. It calls `ComputeResidual`.
4. It adds the result to `LinSysRes` and `Jacobian`.

`CDriver::InitializeNumerics` (`SU2_CFD/src/drivers/CDriver.cpp`) creates the numerics objects. It uses the helpers `InstantiateTurbulentNumerics`, `InstantiateTransitionNumerics` and `InstantiateSpeciesNumerics`.

### Geometry

[`CGeometry`](../../doxygen/classCGeometry.html) holds the mesh and the dual grid. The dual grid is the set of control volumes, one around each mesh point, that the edge-based finite volume method uses. The main members are:

- `nodes` ([`CPoint`](../../doxygen/classCPoint.html)): coordinates, volumes and neighbors of all points.
- `edges` ([`CEdge`](../../doxygen/classCEdge.html)): the two points and the normal of each edge.
- `vertex[iMarker][iVertex]` ([`CVertex`](../../doxygen/classCVertex.html)): boundary points and their normals.
- `InitiateComms` / `CompleteComms`: halo exchange between MPI ranks (see [MPI](#parallelism-mpi-openmp-and-vectorization) below).

Children:

- [`CPhysicalGeometry`](../../doxygen/classCPhysicalGeometry.html): the mesh read from file, partitioned with ParMETIS. The file is in SU2 or CGNS format; the readers are in `Common/include/geometry/meshreader`.
- [`CMultiGridGeometry`](../../doxygen/classCMultiGridGeometry.html): the coarse levels. Each one is made by agglomeration, which merges neighboring control volumes of the finer level into larger ones.
- [`CDummyGeometry`](../../doxygen/classCDummyGeometry.html): empty geometry for the dry run.
- [`CMeshFEM_DG`](../../doxygen/classCMeshFEM__DG.html): geometry of the DG solver (`Common/include/fem`).

Grid movement (`CSurfaceMovement`, `CVolumetricMovement`, `CFreeFormDefBox`) is in `Common/include/grid_movement`. Interpolation between non-matching zone interfaces (`CInterpolator` and its children) is in `Common/include/interface_interpolation`.

### Configuration

[`CConfig`](../../doxygen/classCConfig.html) (`Common/include/CConfig.hpp`, `Common/src/CConfig.cpp`) reads the configuration file and stores every option. There is one `CConfig` per zone. There is also a `driver_config`, read from the main configuration file. It holds the settings of the driver, for example for multizone problems.

- `CConfig::SetConfig_Options` declares the options with `addDoubleOption`, `addBoolOption`, `addEnumOption`, `addStringListOption`, ...
- The enums, and the maps from option strings to enum values, are in `Common/include/option_structure.hpp`. The option parser classes are in `Common/include/option_structure.inl`.
- `CConfig::SetPostprocessing` and `CConfig::SetMarkers` compute the derived settings.

See [Adding a Config Option](../Adding-a-Config-Option/) and the [Option Reference](../Option-Reference/).

### Output

[`COutput`](../../doxygen/classCOutput.html) writes the screen output, the history file, and the volume and surface files. Each solver kind has a child class, created by `COutputFactory`. Examples are [`CFlowCompOutput`](../../doxygen/classCFlowCompOutput.html), [`CFlowIncOutput`](../../doxygen/classCFlowIncOutput.html), [`CHeatOutput`](../../doxygen/classCHeatOutput.html) and [`CElasticityOutput`](../../doxygen/classCElasticityOutput.html). [`CMultizoneOutput`](../../doxygen/classCMultizoneOutput.html) combines the zones of a multizone run. The file writers (Paraview, Tecplot, CGNS, CSV, SU2 restart, ...) are in `SU2_CFD/include/output/filewriter`.

See [Adding Output Fields](../Adding-Output-Fields/) for the details.

## Linear algebra

Files: `Common/include/linear_algebra` and `Common/src/linear_algebra`.

- [`CSysMatrix`](../../doxygen/classCSysMatrix.html): block sparse matrix, with one block per pair of neighbor points. Implicit solvers use it for the Jacobian.
- [`CSysVector`](../../doxygen/classCSysVector.html): block vector, for the residual and the solution update.
- [`CSysSolve`](../../doxygen/classCSysSolve.html): Krylov solvers (`FGMRES_LinSolver`, `RFGMRES_LinSolver`, `FGCRODR_LinSolver`, `BCGSTAB_LinSolver`, `CG_LinSolver`, `Smoother_LinSolver`). `Solve` reads the method and the preconditioner from the configuration. For the primal flow these are `LINEAR_SOLVER` and `LINEAR_SOLVER_PREC`. The mesh, adjoint and other problems have their own options.
- [`CPreconditioner`](../../doxygen/classCPreconditioner.html) and its children [`CJacobiPreconditioner`](../../doxygen/classCJacobiPreconditioner.html), [`CILUPreconditioner`](../../doxygen/classCILUPreconditioner.html), [`CLU_SGSPreconditioner`](../../doxygen/classCLU__SGSPreconditioner.html), [`CLineletPreconditioner`](../../doxygen/classCLineletPreconditioner.html), [`CPastixPreconditioner`](../../doxygen/classCPastixPreconditioner.html). `CPreconditioner::Create` builds the chosen one.
- [`CMatrixVectorProduct`](../../doxygen/classCMatrixVectorProduct.html): the product used inside the Krylov solvers.

Each solver has the members `Jacobian` (`CSysMatrix<su2mixedfloat>`), `LinSysRes`, `LinSysSol` and `System` (`CSysSolve<su2mixedfloat>`). `su2mixedfloat` is `float` when SU2 is built with `-Denable-mixedprec=true`. See also [Linear Solvers and Preconditioners](../Linear-Solvers-and-Preconditioners/).

## Containers and toolboxes

Containers (`Common/include/containers`):

- [`C2DContainer`](../../doxygen/classC2DContainer.html) (`C2DContainer.hpp`): aligned 2D storage. The same file defines the common type aliases:
  - `su2vector<T>`, `su2matrix<T>` (row-major), `ColMajorMatrix<T>`;
  - `su2activevector`, `su2activematrix` (of `su2double`);
  - `su2passivevector`, `su2passivematrix` (of `passivedouble`).
- `container_decorators.hpp`: `C3DDoubleMatrix` (also called `CVectorOfMatrix`, used for gradients) and `CMatrixView`.
- `CLookUpTable`, `CFileReaderLUT`, `CTrapezoidalMap`: look-up tables. The flamelet fluid model (`CFluidFlamelet`) and the data-driven fluid model (`CDataDrivenFluid`) use them.
- `CVertexMap`, `CFastFindAndEraseQueue`, `CPyWrapperMatrixView`.

Toolboxes (`Common/include/toolboxes`):

- `geometry_toolbox.hpp` (vector operations);
- `graph_toolbox.hpp` (sparse patterns, coloring);
- `printing_toolbox.hpp` (`PrintingToolbox::CTablePrinter`);
- `ndflattener.hpp`, `CSquareMatrixCM`, `CSymmetricMatrix`, `C1DInterpolation`, `CLinearPartitioner`, `CQuasiNewtonInvLeastSquares`;
- the verification solutions in `toolboxes/MMS`.

## Parallelism: MPI, OpenMP and vectorization

Files: `Common/include/parallelization`.

**MPI.** Use the `SU2_MPI` wrapper from `mpi_structure.hpp`, never MPI directly. `SU2_MPI` is `CBaseMPIWrapper`, or `CMediMPIWrapper` (MeDiPack) in AD builds. This way MPI calls work with `su2double`.

SU2 partitions the mesh at start-up. Each rank also stores halo points: copies of the points next to its partition that belong to other ranks. To update the halo values, call `InitiateComms` / `CompleteComms` of `CGeometry` and `CSolver`. An `MPI_QUANTITIES` value says what to send.

**OpenMP.** Do not include `omp.h`. Include `omp_structure.hpp` and use its macros. They compile to nothing when OpenMP is off (`-Dwith-omp=false`), and to OpDiLib calls in AD builds:

- `SU2_OMP_PARALLEL` ... `END_SU2_OMP_PARALLEL`, `SU2_OMP_PARALLEL_(ARGS)`
- `SU2_OMP_FOR_STAT(CHUNK)`, `SU2_OMP_FOR_DYN(CHUNK)`, `SU2_OMP_FOR_(ARGS)` ... `END_SU2_OMP_FOR`
- `SU2_OMP_MASTER` ... `END_SU2_OMP_MASTER`, `SU2_OMP_CRITICAL` ... `END_SU2_OMP_CRITICAL`, `SU2_OMP_BARRIER`, `SU2_OMP_ATOMIC`
- `SU2_OMP_SAFE_GLOBAL_ACCESS(...)`: run code on one thread, with barriers around it
- `SU2_OMP_SIMD`
- helpers: `omp_get_num_threads()`, `omp_get_thread_num()`, `roundUpDiv`, `atomicAdd`

When the solver supports it (`CSolver::GetHasHybridParallel()`), the integrations open one parallel region for a whole iteration.

Edge loops write to both points of an edge. Two threads could then write to the same point at the same time (a race condition). To avoid this, SU2 groups the edges by color (`CGeometry::GetEdgeColoring`, option `EDGE_COLORING_GROUP_SIZE`). If the coloring is poor, the flow and scalar solvers first store the edge fluxes, then sum them for each point (`ReducerStrategy`, `EdgeFluxes`).

Set the number of threads per rank with `SU2_CFD -t <threads>`.

**Vectorization.** SIMD (single instruction, multiple data) instructions let the CPU apply one operation to several values at once.

- `vectorization.hpp` defines [`simd::Array`](../../doxygen/classsimd_1_1Array.html), a short fixed-size array with arithmetic operators. The compiler turns its loops into SIMD instructions. `special_vectorization.hpp` has explicit versions for some instruction sets.
- [`CNumericsSIMD`](../../doxygen/classCNumericsSIMD.html) (`SU2_CFD/include/numerics_simd`) uses it to compute the fluxes of several edges at once.
- `CEulerSolver` uses these schemes for all centered schemes. It also uses them for the `ROE` and `MSW` upwind schemes, with an ideal gas and no low-Mach correction. See `CEulerSolver::Centered_Residual` and `CEulerSolver::Upwind_Residual`.
- The build option `cpu-arch` (default `native`) sets the target instruction set.

## Algorithmic differentiation

SU2 uses operator overloading AD with [CoDiPack](https://github.com/scicompkl/CoDiPack) (`externals/codi`).

- `su2double` (`Common/include/code_config.hpp`) is the type for all values that may need derivatives. What it is depends on the build:
  - in a normal build, it is `passivedouble` (normally `double`);
  - with `-Denable-autodiff=true` (builds `SU2_CFD_AD`, `SU2_DOT_AD`), it is a CoDiPack reverse type;
  - with `-Denable-directdiff=true` (builds `SU2_CFD_DIRECTDIFF`), it is a forward type.
- `passivedouble` is never differentiated. Use `SU2_TYPE::GetValue(x)` to get the plain value of a `su2double`.
- The [`AD`](../../doxygen/namespaceAD.html) namespace (`Common/include/basic_types/ad_structure.hpp`) wraps the tape:
  - `StartRecording`, `StopRecording`, `RegisterInput`, `RegisterOutput`, `ComputeAdjoint`, `ClearAdjoints`;
  - preaccumulation (`StartPreacc`, `SetPreaccIn`, `SetPreaccOut`, `EndPreacc`). Preaccumulation stores only the Jacobian of a code section on the tape, instead of all its operations. This makes the tape smaller.

  The functions do nothing in non-AD builds, so the code can always call them.
- MeDiPack (`externals/medi`) differentiates MPI calls. OpDiLib (`externals/opdi`) differentiates OpenMP regions.

See [Advanced AD Techniques](../Advanced-AD-Techniques/) for more.

## Where do I change...

This section lists the files and functions to change for common tasks.

### New boundary condition

1. Add the kind to `enum BC_TYPE` in `Common/include/option_structure.hpp`.
2. In `Common/src/CConfig.cpp`, add the marker option in `CConfig::SetConfig_Options`, and set `Marker_CfgFile_KindBC` in `CConfig::SetMarkers`.
3. Add a virtual `BC_...` function to `CSolver` (`SU2_CFD/include/solvers/CSolver.hpp`).
4. Implement it in the solvers that need it, for example `CEulerSolver`, `CIncEulerSolver`, and the scalar solvers through `CScalarSolver`.
5. Call it from the `switch` on `GetMarker_All_KindBC` in `CIntegration::Space_Integration`.

### New convective scheme

1. Add the value to `enum class UPWIND` (or `CENTERED`) and to `Upwind_Map` in `option_structure.hpp`.
2. Write the numerics class in `SU2_CFD/include/numerics/flow/convection/` and `SU2_CFD/src/numerics/flow/convection/`.
3. Create it in `CDriver::InitializeNumerics` (switch on `GetKind_Upwind_Flow()`).
4. Optional, for a compressible ideal gas only: write a vectorized version in `SU2_CFD/include/numerics_simd/flow/convection/upwind.hpp`. Create it in `createUpwindIdealNumerics` (`CNumericsSIMD.cpp`) and select it in `CEulerSolver::Upwind_Residual`.

### New turbulence model source term

- The source terms are in `SU2_CFD/include/numerics/turbulent/turb_sources.hpp`: `CSourceBase_TurbSA` and `SAFactory` for SA, `CSourcePieceWise_TurbSST` for SST.
- `CDriver::InstantiateTurbulentNumerics` creates them. `CTurbSASolver::Source_Residual` / `CTurbSSTSolver::Source_Residual` use them.
- Users choose the model variants with `SA_OPTIONS` / `SST_OPTIONS` (`enum class SA_OPTIONS`, `SST_OPTIONS` in `option_structure.hpp`).

### New output field

Change the output class of the solver in `SU2_CFD/src/output/`, or `CFlowOutput` for turbulence and species fields. See [Adding Output Fields](../Adding-Output-Fields/).

### New config option

Add the option in `CConfig::SetConfig_Options`, a member and a getter in `CConfig.hpp`, and an entry in `config_template.cfg`. See [Adding a Config Option](../Adding-a-Config-Option/).

### New solver

Start from `CTemplateSolver` (`SU2_CFD/include/solvers/CTemplateSolver.hpp`) and the template numerics (`CConvective_Template`, `CViscous_Template`, `CSource_Template` in `SU2_CFD/include/numerics/template.hpp`). Then add:

1. a value in `enum class MAIN_SOLVER` and `Solver_Map` (`option_structure.hpp`);
2. the solver in `CSolverFactory::CreateSolverContainer` and `CreateSubSolver`, with its `INTEGRATION_TYPE`;
3. its numerics in `CDriver::InitializeNumerics`;
4. an iteration in `CIterationFactory::CreateIteration`;
5. an output class in `COutputFactory::CreateOutput`.

### Tests and style

When you change code, add or update a regression test (see [Running Regression Tests](../Running-Regression-Tests/)) or a unit test (see [Writing Unit Tests](../Writing-Unit-Tests/)). Follow the [Style Guide](../Style-Guide/).

Last checked against SU2 develop `11c177b904`.
