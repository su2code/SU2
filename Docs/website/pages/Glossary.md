---
title: Glossary
permalink: /docs_v7/Glossary/
---

This page explains the words that SU2 users and developers meet in the code, the configuration file and the
discussions on GitHub. Where it helps, an entry names the class or option, and links to more.

For a map of the source code, see [Code Structure](../Code-Structure/). For all options, see the
[Option Reference](../Option-Reference/). Class names link to the [Code Reference](../Code-Reference/) (Doxygen).

---

- [Mesh and geometry](#mesh-and-geometry)
- [Solvers and numerics](#solvers-and-numerics)
- [Multizone and coupling](#multizone-and-coupling)
- [Output](#output)
- [Parallelism](#parallelism)
- [Algorithmic differentiation](#algorithmic-differentiation)
- [Configuration](#configuration)

---

## Mesh and geometry

**Agglomeration, multigrid level.**
Multigrid solves the flow on a sequence of coarser grids to speed up convergence. SU2 builds each coarse grid by
agglomeration: it merges neighboring control volumes of the finer grid into larger ones
([`CMultiGridGeometry`](../../doxygen/classCMultiGridGeometry.html)). The index `iMesh` selects the level.
`MESH_0` is the finest grid, the one read from the file. The option `MGLEVEL` sets the number of coarse levels
(default 0, no multigrid).

**Dual grid, control volume.**
SU2 solves the equations on control volumes, one around each mesh point. Together they form the dual grid.
[`CGeometry`](../../doxygen/classCGeometry.html) builds it from the mesh elements. The volume of each control
volume is stored in [`CPoint`](../../doxygen/classCPoint.html) (`nodes->GetVolume(iPoint)`).

**Edge.**
An edge joins two neighboring mesh points. Each edge carries the normal of the face that the two control volumes
share. Flux loops in the finite volume solvers run over edges.
See [`CEdge`](../../doxygen/classCEdge.html) (`geometry->edges`, `geometry->GetnEdge()`).

**Halo point, ghost point.**
After the mesh is split between MPI ranks, each rank also stores copies of the points next to its part that
belong to other ranks. These are halo (or ghost) points. The rank does not compute them; it receives their values
from the owner. `geometry->GetnPointDomain()` counts the points the rank owns, `geometry->GetnPoint()` counts all
points, halos included. The owned points come first, so a loop up to `nPointDomain` visits only owned points.
`nodes->GetDomain(iPoint)` is `false` for a halo point. See
[Parallelism](../Code-Structure/#parallelism-mpi-openmp-and-vectorization).

**Interface (marker).**
A marker where two zones meet. `MARKER_ZONE_INTERFACE` lists the interface markers in pairs, one marker from each
zone. `MARKER_FLUID_INTERFACE` and `MARKER_CHT_INTERFACE` set the kind of boundary condition on such markers.
When the meshes do not match at the interface, an interpolator
([`CInterpolator`](../../doxygen/classCInterpolator.html), option `KIND_INTERPOLATION`) maps one side to the
other. See [Multizone](../Multizone/).

**Marker.**
A named group of boundary faces, for example a wall or a far field. The mesh file defines the markers. The
`MARKER_*` options of the configuration file set their boundary conditions. In the code, `iMarker` indexes them,
and `config->GetMarker_All_TagBound(iMarker)` gives the name. See [Markers and BC](../Markers-and-BC/).

**Periodic pair.**
Two markers that SU2 treats as the same surface, moved by a rotation or a translation. `MARKER_PERIODIC` lists,
for each pair, the periodic marker, its donor marker, the rotation center, the rotation angles and the
translation. SU2 adds the points of the other side as extra halo points. The solvers exchange data across the pair
with `InitiatePeriodicComms` / `CompletePeriodicComms` of [`CSolver`](../../doxygen/classCSolver.html).

**Vertex.**
A boundary point seen from one marker. `geometry->vertex[iMarker][iVertex]` is a
[`CVertex`](../../doxygen/classCVertex.html). It holds the index of the mesh point (`GetNode()`) and the boundary
normal (`GetNormal()`). A mesh point at a corner belongs to several markers, so it has one vertex on each.

**Zone.**
A part of the domain with its own mesh, its own configuration and its own solvers, for example the fluid and the
solid in a conjugate heat transfer problem. `iZone` indexes the zones. The driver holds one
[`CConfig`](../../doxygen/classCConfig.html) per zone. Set `MULTIZONE= YES` and `CONFIG_LIST` to use several zones.
See [Multizone](../Multizone/).

## Solvers and numerics

**CFL adaptation.**
The CFL (Courant-Friedrichs-Lewy) number sets the size of the local pseudo time step. With `CFL_ADAPT= YES`, SU2 raises the CFL number when
the iterations go well, and lowers it when they do not (`CSolver::AdaptCFLNumber`). `CFL_ADAPT_PARAM` sets the
factors and the limits. The history fields `MIN_CFL`, `MAX_CFL` and `AVG_CFL` show the result.

**Dual time stepping.**
The method that SU2 uses for most unsteady problems. Each physical time step is solved as a steady problem in a
pseudo time, with inner iterations. Set it with `TIME_MARCHING= DUAL_TIME_STEPPING-1ST_ORDER` or
`DUAL_TIME_STEPPING-2ND_ORDER`. `CSolver::SetResidual_DualTime` adds the time terms to the residual.

**Edge kernel.**
The code that computes the flux across one edge from the states of its two points. Examples are
`CNumerics::ComputeResidual`, the vectorized `CNumericsSIMD::ComputeFlux`, and the scalar flux classes such as
`CScalarFlux_SA` and `CScalarFlux_SST` (`SU2_CFD/include/numerics/**/*_edge_flux.hpp`). The solver calls the
kernel inside its edge loop.

**Inner and outer iterations.**
Inner iterations (`INNER_ITER`) solve one zone, or one time step, towards convergence. Outer iterations
(`OUTER_ITER`) repeat the exchange of data between zones in a multizone problem. Time iterations (`TIME_ITER`)
count the physical time steps. In a single-zone steady problem `ITER` is the same as `INNER_ITER`. SU2 stops with
an error if you use `ITER` in a multizone or unsteady problem.

**Jacobian.**
The matrix of the derivatives of the residual with respect to the solution. Implicit solvers use it to compute
the update of each iteration. In SU2 it is a block sparse matrix, with one block per pair of neighbor points:
`Jacobian` of [`CSolver`](../../doxygen/classCSolver.html), of type
[`CSysMatrix`](../../doxygen/classCSysMatrix.html). It is often approximate: it only has to make the iterations
converge.

**Limiter.**
A factor between 0 and 1 that reduces the reconstructed gradient near shocks and strong changes. It prevents new
extrema and oscillations. Set it with `SLOPE_LIMITER_FLOW`, `SLOPE_LIMITER_TURB`, ... (for example
`VENKATAKRISHNAN`). The code is in `SU2_CFD/include/limiters` (`computeLimiters`). The volume output group
`LIMITER` writes the values. See [Slope Limiters](../Slope-Limiters-and-Shock-Resolution/).

**MUSCL.**
Monotonic Upstream-centered Scheme for Conservation Laws: a way to reach second order accuracy with an upwind scheme. Before the flux is computed, the states at the two
ends of an edge are extrapolated to the middle of the edge with the gradients. Set it with `MUSCL_FLOW`,
`MUSCL_TURB`, ... The gradient method is `NUM_METHOD_GRAD_RECON` (or `NUM_METHOD_GRAD` if that is not set). See
[Convective Schemes](../Convective-Schemes/).

**Numerics object.**
An object of a [`CNumerics`](../../doxygen/classCNumerics.html) child class. It computes one term, such as the
convective flux of an edge or the source term of a point, with its Jacobians. The driver creates them in
`CDriver::InitializeNumerics`. They are stored in `numerics_container[iZone][iInst][iMesh][iSol][iTerm]`, where
`iTerm` is `CONV_TERM`, `VISC_TERM`, `SOURCE_FIRST_TERM`, ... See
[Adding a Numerical Scheme](../Adding-a-Numerical-Scheme/).

**Residual.**
The sum of the fluxes and source terms of a control volume. It is zero when the steady solution is reached. The
solvers store it in `LinSysRes`. The history fields in the groups `RMS_RES` and `MAX_RES` show the base 10
logarithm of its root mean square and maximum over all points, one value per equation.

**Runtime system (`RUNTIME_*_SYS`).**
A value of `enum RUNTIME_TYPE` (`Common/include/option_structure.hpp`), such as `RUNTIME_FLOW_SYS`,
`RUNTIME_TURB_SYS` or `RUNTIME_HEAT_SYS`. The iterations pass it to the integrations to say which set of equations
to advance. `CConfig::GetContainerPosition` turns it into the solver position, for example
`RUNTIME_TURB_SYS` into `TURB_SOL`.

**Solver container.**
The array of all solver objects of a run: `solver_container[iZone][iInst][iMesh][iSol]` in the driver. Inside the
solvers and the integrations you mostly see the last level, `CSolver** solver_container`, indexed by the solver
position. See [Code Structure](../Code-Structure/#the-containers-in-the-driver).

**Solver position (`FLOW_SOL`, ...).**
The index of one solver in the solver container: `FLOW_SOL`, `TURB_SOL`, `TRANS_SOL`, `SPECIES_SOL`, `HEAT_SOL`,
`RAD_SOL`, `MESH_SOL`, `ADJFLOW_SOL`, ... (`enum SOLVER_TYPE` in `Common/include/option_structure.hpp`). Some
values are shared: `FEA_SOL` is also 0, like `FLOW_SOL`, because a zone never has both.

## Multizone and coupling

**Driver.**
The object that runs the whole simulation. It creates the geometry, solvers, numerics and outputs, and runs the
time loop. [`CSinglezoneDriver`](../../doxygen/classCSinglezoneDriver.html) and
[`CMultizoneDriver`](../../doxygen/classCMultizoneDriver.html) are the main ones; all derive from
[`CDriver`](../../doxygen/classCDriver.html). See [Code Structure](../Code-Structure/#drivers).

**Instance.**
A copy of the solvers of one zone. Only harmonic balance uses more than one: one instance per time instance (a sampled time of the period).
`iInst` indexes them. In all other cases there is one instance.

**Interface (object).**
An object of a [`CInterface`](../../doxygen/classCInterface.html) child class. It knows what to send from one
zone to another, for example the tractions from a fluid to a structure
([`CFlowTractionInterface`](../../doxygen/classCFlowTractionInterface.html)), or the temperature and heat flux in
conjugate heat transfer. The driver stores them in `interface_container[donorZone][targetZone]`.

**Iteration.**
An object of a [`CIteration`](../../doxygen/classCIteration.html) child class. It runs one iteration of one
physics in one zone: `Preprocess`, `Iterate`, `Postprocess`, `Monitor`, `Output`. For example,
[`CFluidIteration`](../../doxygen/classCFluidIteration.html) advances the flow and then the turbulence, transition
and species solvers.

**Transfer.**
The exchange of boundary data between two zones, through an interface object. `CDriver::TransferData` (and
`CMultizoneDriver::TransferData`) call it. `CInterface::BroadcastData` reads the donor values
(`GetDonor_Variable`), maps them with the interpolator and sets them in the target solver
(`SetTarget_Variable`). The kind of transfer is a value of `enum ENUM_TRANSFER`, such as `FLOW_TRACTION` or
`CONJUGATE_HEAT_FS`.

## Output

**History field.**
A single number per iteration, such as a residual, a force coefficient or the time. History fields go to the
history file (`history.csv` by default). Choose them with `HISTORY_OUTPUT`. In the code they are added with
`COutput::AddHistoryOutput`. See [Custom Output](../Custom-Output/) and the
[Output Fields Reference](../Output-Fields-Reference/).

**Output group.**
A name that stands for several fields, for example `RMS_RES`, `AERO_COEFF`, `SOLUTION` or `PRIMITIVE`. Each field
belongs to one group (an argument of `AddHistoryOutput` and `AddVolumeOutput`). Writing the group name in
`SCREEN_OUTPUT`, `HISTORY_OUTPUT` or `VOLUME_OUTPUT` selects all its fields.

**Screen field.**
A history field that is also printed on the screen. Choose them with `SCREEN_OUTPUT`.

**Volume field.**
A value at each mesh point, such as the pressure or the residual. Volume fields go to the volume and surface files
(Paraview, Tecplot, CGNS, ...). Choose them with `VOLUME_OUTPUT`. In the code they are added with
`COutput::AddVolumeOutput`. See [Adding Output Fields](../Adding-Output-Fields/).

## Parallelism

**Edge coloring (OpenMP).**
An edge loop writes to both points of each edge. If two threads handle edges that share a point, they write to
the same memory at the same time. To avoid this, SU2 groups the edges by color: edges of one color share no
point, so the threads can work on them safely (`CGeometry::GetEdgeColoring`). The option
`EDGE_COLORING_GROUP_SIZE` sets the size of the groups of edges.

**MPI rank.**
One of the processes of a parallel run started with `mpirun`. Each rank works on one part of the mesh. Get its
number with `SU2_MPI::GetRank()`. Rank 0 is `MASTER_NODE`, which prints the screen output. Use the `SU2_MPI`
wrapper, never MPI directly.

**Partition.**
The part of the mesh that one rank owns. SU2 splits the mesh at start-up with ParMETIS
(`CPhysicalGeometry::SetColorGrid_Parallel`). The volume field `RANK` (group `MPI`) shows the partitions.

**Reducer strategy.**
The fallback when the edge coloring is poor, that is, when its expected parallel efficiency is below
`COLORING_EFF_THRESH` (0.875). The solver first stores the flux of every edge (`EdgeFluxes`), then sums them for
each point in a point loop. The member `ReducerStrategy` of the flow and scalar solvers says if it is in use.
`EDGE_COLORING_GROUP_SIZE= 0` forces it.

**SIMD.**
Single instruction, multiple data: CPU instructions that apply one operation to several values at once. SU2 uses
them through [`simd::Array`](../../doxygen/classsimd_1_1Array.html) and
[`CNumericsSIMD`](../../doxygen/classCNumericsSIMD.html), which compute the fluxes of several edges together. The
build option `cpu-arch` sets which instructions the compiler may use.

## Algorithmic differentiation

See [Automatic Differentiation](../Automatic-Differentiation/) for the full guide.

**Preaccumulation.**
A way to make the tape smaller. For a code section with few inputs and outputs, SU2 stores only the Jacobian of
the section on the tape, not all its operations (`AD::StartPreacc`, `AD::SetPreaccIn`, `AD::SetPreaccOut`,
`AD::EndPreacc`).

**`su2double`.**
The type of every value that may need derivatives. In a normal build it is a plain `double`. In an AD build it is
a CoDiPack type that records its operations. `passivedouble` is never differentiated.

**Tape.**
The record of all operations done on `su2double` values while recording is on. The discrete adjoint runs the tape
backwards to get derivatives (`AD::StartRecording`, `AD::StopRecording`, `AD::ComputeAdjoint`).

## Configuration

**Config file.**
The text file with the options of a run, one `NAME= value` per line. `%` starts a comment. A line that ends with
`\` continues on the next line. [`CConfig`](../../doxygen/classCConfig.html) reads it. `config_template.cfg` in
the SU2 folder lists every option. See [Configuration File](../Configuration-File/).

**Keyword map.**
The table that turns the text value of an option into a value of a C++ enum, for example `Solver_Map` turns
`RANS` into `MAIN_SOLVER::RANS`. The maps are in `Common/include/option_structure.hpp`, built with `MapType` and
`MakePair`. They exist only where `ENABLE_MAPS` is defined, such as in `CConfig.cpp`; elsewhere they compile to
nothing, which keeps the build fast.

**Option.**
One setting of the config file. `CConfig::SetConfig_Options` declares each option with `addDoubleOption`,
`addEnumOption`, `addStringListOption`, ... and stores it in `option_map`. Each option is an object of a
`COptionBase` child class (`Common/include/option_structure.inl`). See
[Adding a Config Option](../Adding-a-Config-Option/).

**Sub-config.**
A config file for one zone of a multizone problem, listed in `CONFIG_LIST` of the main config file. It holds only
the options that differ in that zone. The zone takes the other options from the main file
(`CConfig::SetDefaultFromConfig`), except `SCREEN_OUTPUT` and `HISTORY_OUTPUT`. See
[Multizone](../Multizone/#sub-config-files).

Last checked against SU2 develop `11c177b904`.
