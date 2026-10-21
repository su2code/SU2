---
title: Adding a Boundary Condition
permalink: /docs_v7/Adding-a-Boundary-Condition/
---

This page shows how to add a boundary condition (BC) to the finite volume solvers of `SU2_CFD`. It covers:

- the marker option;
- how the marker reaches the solvers;
- how the BC is called in each iteration;
- how to write the BC for the flow solver, and for the turbulence, species and heat solvers.

The examples are copied from the code. See [Code Structure](../Code-Structure/) for the classes, and [Adding a Config Option](../Adding-a-Config-Option/) for options in general.

## Files you will touch

| File | What goes there |
|---|---|
| `Common/include/option_structure.hpp` | a value in `enum BC_TYPE` |
| `Common/include/CConfig.hpp` | the marker count and names, the BC data, their getters |
| `Common/src/CConfig.cpp` | the `MARKER_...` option in `CConfig::SetConfig_Options`, the marker list in `CConfig::SetMarkers`, the getters, the table printed by `CConfig::SetOutput`, `GetViscous_Wall` / `GetSolid_Wall` for walls |
| `config_template.cfg` | the option with a short explanation |
| `SU2_CFD/src/integration/CIntegration.cpp` | the call in `CIntegration::Space_Integration` |
| `SU2_CFD/include/solvers/CSolver.hpp` | a virtual `BC_...` function with an empty body |
| the solver files, e.g. `SU2_CFD/src/solvers/CEulerSolver.cpp`, `SU2_CFD/src/solvers/CIncEulerSolver.cpp`, `SU2_CFD/src/solvers/CTurbSASolver.cpp`, `SU2_CFD/src/solvers/CTurbSSTSolver.cpp`, `SU2_CFD/src/solvers/CSpeciesSolver.cpp`, `SU2_CFD/src/solvers/CHeatSolver.cpp` | the implementations |
| `TestCases/` | a regression test |

---

- [1. Add the BC kind](#1-add-the-bc-kind)
- [2. Add the marker option](#2-add-the-marker-option)
- [3. How the marker reaches the solvers](#3-how-the-marker-reaches-the-solvers)
- [4. Call the BC in Space_Integration](#4-call-the-bc-in-space_integration)
- [5. Declare the BC in CSolver](#5-declare-the-bc-in-csolver)
- [6. Implement it in the flow solvers](#6-implement-it-in-the-flow-solvers)
- [7. Strong (Dirichlet) and weak BCs](#7-strong-dirichlet-and-weak-bcs)
- [8. Turbulence, species and heat solvers](#8-turbulence-species-and-heat-solvers)
- [9. Other places that test the BC kind](#9-other-places-that-test-the-bc-kind)
- [10. Multigrid](#10-multigrid)
- [11. Adjoint](#11-adjoint)
- [12. Output and Python](#12-output-and-python)
- [13. Test](#13-test)
- [Checklist](#checklist)

---

## 1. Add the BC kind

Add a value for the new kind to `enum BC_TYPE`. The kinds are an unscoped enum with fixed integer values in `Common/include/option_structure.hpp`:

```cpp
enum BC_TYPE {
  EULER_WALL = 1,             /*!< \brief Boundary Euler wall definition. */
  FAR_FIELD = 2,              /*!< \brief Boundary far-field definition. */
  SYMMETRY_PLANE = 3,         /*!< \brief Boundary symmetry plane definition. */
  INLET_FLOW = 4,             /*!< \brief Boundary inlet flow definition. */
  OUTLET_FLOW = 5,            /*!< \brief Boundary outlet flow definition. */
  // ...
  SUPERSONIC_INLET = 19,      /*!< \brief Boundary supersonic inlet definition. */
  SUPERSONIC_OUTLET = 20,     /*!< \brief Boundary supersonic inlet definition. */
  // ...
  CHT_WALL_INTERFACE = 50,    /*!< \brief Domain interface definition. */
  SMOLUCHOWSKI_MAXWELL = 55,  /*!< \brief Smoluchoski/Maxwell wall boundary condition. */
  SEND_RECEIVE = 99,          /*!< \brief Boundary send-receive definition. */
};
```

Pick a value that is not used yet. `SEND_RECEIVE` is reserved for the markers that hold the MPI halo points (local copies of points that another rank owns). `CConfig` stores the kinds as `unsigned short`.

## 2. Add the marker option

Add the `MARKER_...` option that users write in the configuration file. The rest of this page follows two examples: `MARKER_SUPERSONIC_OUTLET` (no data) and `MARKER_SUPERSONIC_INLET` (with data).

1. **Members.** In `Common/include/CConfig.hpp`, add the count and the names to the private members, next to the others:

   ```cpp
     nMarker_Supersonic_Outlet,      /*!< \brief Number of supersonic outlet flow markers. */
     // ...
     *Marker_Supersonic_Outlet,      /*!< \brief Supersonic outlet flow markers. */
   ```

2. **Option.** Add the option in `CConfig::SetConfig_Options` (`Common/src/CConfig.cpp`):
   - a list of names uses `addStringListOption`;
   - a BC with values uses one of the special parsers, for example `addStringDoubleListOption` (one value per marker, as in `MARKER_OUTLET` and `MARKER_ISOTHERMAL`) or `addInletOption`.

   ```cpp
     /*!\brief MARKER_SUPERSONIC_INLET  \n DESCRIPTION: Supersonic inlet boundary marker(s)
      * \n   Format: (inlet marker, temperature, static pressure, velocity_x,   velocity_y, velocity_z, ... ), i.e. primitive variables specified. \ingroup Config*/
     addInletOption("MARKER_SUPERSONIC_INLET", nMarker_Supersonic_Inlet, Marker_Supersonic_Inlet, Inlet_Temperature, Inlet_Pressure, Inlet_Velocity);
     /*!\brief MARKER_SUPERSONIC_OUTLET \n DESCRIPTION: Supersonic outlet boundary marker(s) \ingroup Config*/
     addStringListOption("MARKER_SUPERSONIC_OUTLET", nMarker_Supersonic_Outlet, Marker_Supersonic_Outlet);
   ```

   The parsers are the `COption...` classes in `Common/include/option_structure.inl`. See [Adding a Config Option](../Adding-a-Config-Option/) for the comment format and `config_template.cfg`.

3. **Marker list.** In `CConfig::SetMarkers`, add the count to the sum that sizes the arrays:

   ```cpp
     nMarker_CfgFile = nMarker_Euler + nMarker_FarField + nMarker_SymWall +
     // ...
     nMarker_Supersonic_Inlet + nMarker_Supersonic_Outlet + nMarker_Displacement + nMarker_Load +
     // ...
   ```

   Then add a loop that writes the names and the kind:

   ```cpp
     for (iMarker_Supersonic_Outlet = 0; iMarker_Supersonic_Outlet < nMarker_Supersonic_Outlet; iMarker_Supersonic_Outlet++) {
       Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Supersonic_Outlet[iMarker_Supersonic_Outlet];
       Marker_CfgFile_KindBC[iMarker_CfgFile] = SUPERSONIC_OUTLET;
       iMarker_CfgFile++;
     }
   ```

   Declare the loop counter (`iMarker_Supersonic_Outlet`) at the top of the function.

4. **Getters for the data.** The data is stored per marker of the option, in the order of the configuration file. The getters search the marker by name, for example:

   ```cpp
   su2double CConfig::GetOutlet_Pressure(const string& val_marker) const {
     unsigned short iMarker_Outlet;
     for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++)
       if (Marker_Outlet[iMarker_Outlet] == val_marker) break;
     return Outlet_Pressure[iMarker_Outlet];
   }
   ```

   The solvers call them with `config->GetMarker_All_TagBound(val_marker)`. The search is a loop over strings. So call the getter once per marker, outside the vertex loop.

5. **Screen table.** At start-up, `CConfig::SetOutput` prints one line per BC kind (`BoundaryTable << "Supersonic outlet boundary";`). Add one for the new kind.

## 3. How the marker reaches the solvers

You need no new code in this step for a new kind. This section explains how the marker data reaches the solvers, and which functions the solvers use to read it.

`CConfig` keeps two marker lists:

| List | Indexed by | Filled in |
|---|---|---|
| `Marker_CfgFile_TagBound`, `Marker_CfgFile_KindBC`, ... | the markers of the configuration file | `CConfig::SetMarkers` |
| `Marker_All_TagBound`, `Marker_All_KindBC`, `Marker_All_SendRecv`, ... | the markers of the local mesh partition (`iMarker` in the solvers) | the geometry, after partitioning |

`CPhysicalGeometry` copies the kind of each physical marker of the mesh. It sets `SEND_RECEIVE` for the halo markers:

```cpp
    if (Marker_Tag != "SEND_RECEIVE") {
      /*--- Update config information storing the boundary information in the right place ---*/

      Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
      config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
      // ...
    }
```

If a mesh marker is in no `MARKER_...` option, the run stops in `CConfig::GetMarker_CfgFile_TagBound` ("The configuration file doesn't have any definition for marker ...").

In the solvers, use:

- `config->GetnMarker_All()`, `config->GetMarker_All_KindBC(iMarker)`, `config->GetMarker_All_TagBound(iMarker)`;
- `geometry->nVertex[iMarker]` and `geometry->vertex[iMarker][iVertex]` (`GetNode()`, `GetNormal(...)`, `GetNormal_Neighbor()`);
- `geometry->nodes->GetDomain(iPoint)`, which is `false` for halo points. A marker of a partition can contain halo points. Their values come from the owning rank through MPI, so the BC skips them.

## 4. Call the BC in Space_Integration

Add a `case` for the new kind in `CIntegration::Space_Integration` (`SU2_CFD/src/integration/CIntegration.cpp`). All finite volume solvers call their BCs from this function. It serves:

- the flow solver, from `CMultiGridIntegration` on every grid level, and from `CNewtonIntegration`;
- the turbulence, transition, species, heat, radiation and Poisson solvers, from `CSingleGridIntegration`.

The equation system selects the solver:

```cpp
  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);
```

So one `case` serves every solver. Each solver decides what the BC means for its equations, through its own `BC_...` function.

The order in `Space_Integration` is:

1. `Centered_Residual` or `Upwind_Residual`, `Viscous_Residual`, `Source_Residual`, and `SetResidual_DualTime` for dual time stepping.
2. `BC_Fluid_Interface` (sliding interfaces) and the Giles preprocessing.
3. The **weak** BCs, which add a boundary flux (see [step 7](#7-strong-dirichlet-and-weak-bcs)). They are one `switch` over all markers:

   ```cpp
     for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
       KindBC = config->GetMarker_All_KindBC(iMarker);
       switch (KindBC) {
         // ...
         case SUPERSONIC_INLET:
           solver_container[MainSolver]->BC_Supersonic_Inlet(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
           break;
         case OUTLET_FLOW:
           solver_container[MainSolver]->BC_Outlet(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
           break;
         case SUPERSONIC_OUTLET:
           solver_container[MainSolver]->BC_Supersonic_Outlet(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
           break;
         // ...
       }
     }
   ```

4. `Impose_Fixed_Values`.
5. The **strong** BCs, which impose values at the points. They are a second `switch` (`ISOTHERMAL`, `HEAT_FLUX`, `HEAT_TRANSFER`, `CUSTOM_BOUNDARY`, `CHT_WALL_INTERFACE`, `SMOLUCHOWSKI_MAXWELL`).
6. `BC_Periodic`, if there are periodic markers.
7. `BC_Sym_Plane` and `BC_Euler_Wall`, last.

Add your `case` to the weak or the strong `switch`. A BC that overwrites the solution or the residual of its points belongs in the strong one. It then runs after the weak BCs of neighbouring markers that share corner points.

The numerics objects passed in are those of the current thread:

```cpp
  CNumerics* conv_bound_numerics = numerics[CONV_BOUND_TERM + omp_get_thread_num()*MAX_TERMS];
  CNumerics* visc_bound_numerics = numerics[VISC_BOUND_TERM + omp_get_thread_num()*MAX_TERMS];
```

A kind with no `case` is not applied at all. Some cases depend on the grid level: `GILES_BOUNDARY` and the turbomachinery form of `RIEMANN_BOUNDARY` only run on `MESH_0`.

Other solvers have their own dispatch: the structural solver in `CStructuralIntegration` (`CLAMPED_BOUNDARY`, `LOAD_BOUNDARY`, ...) and the DG solver in `CFEM_DG_EulerSolver`.

## 5. Declare the BC in CSolver

Declare the new BC function in the base class, so that every solver has it. Add a virtual function with an empty body to `SU2_CFD/include/solvers/CSolver.hpp`, with the usual signature:

```cpp
  inline virtual void BC_Supersonic_Outlet(CGeometry *geometry,
                                           CSolver **solver_container,
                                           CNumerics *conv_numerics,
                                           CNumerics *visc_numerics,
                                           CConfig *config,
                                           unsigned short val_marker) { }
```

The empty default means: **a solver that does not override the function adds nothing at this boundary**. In an edge-based finite volume scheme, a boundary face that gets no flux acts as a wall with zero convective and zero diffusive flux for that equation. Decide for every solver whether that is right (see [step 8](#8-turbulence-species-and-heat-solvers)).

Declare the override in each solver header that implements it (for example `SU2_CFD/include/solvers/CEulerSolver.hpp`), with `override` or `final`.

## 6. Implement it in the flow solvers

Write the BC for the flow solvers. The flow BCs are in these classes:

- compressible: `CEulerSolver` and `CNSSolver`;
- incompressible: `CIncEulerSolver` and `CIncNSSolver`;
- shared by both: `CFVMFlowSolverBase` (`BC_Sym_Plane`, `BC_Euler_Wall`, `BC_Custom`, `BC_Periodic`, `BC_Fluid_Interface`, in `SU2_CFD/include/solvers/CFVMFlowSolverBase.inl`);
- NEMO: the NEMO solvers have their own.

A weak BC builds a ghost state at each vertex: an outer state, as if there were a point outside the boundary. It then computes the boundary flux with the convective numerics, as for an edge between the boundary point and its ghost. `CEulerSolver::BC_Supersonic_Inlet` is a short example (some lines removed):

```cpp
void CEulerSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                                       CNumerics *conv_numerics, CNumerics *visc_numerics,
                                       CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED
  const su2double Gas_Constant = config->GetGas_ConstantND();
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  // ...

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- Retrieve the inlet profile, note that total conditions are reused as static. ---*/

    const su2double Temperature = Inlet_Ttotal[val_marker][iVertex] / config->GetTemperature_Ref();
    const su2double Pressure = Inlet_Ptotal[val_marker][iVertex] / config->GetPressure_Ref();
    // ... velocity, density, energy

    /*--- Primitive variables, using the derived quantities. ---*/

    auto* V_inlet = GetCharacPrimVar(val_marker, iVertex);
    V_inlet[prim_idx.Temperature()] = Temperature;
    V_inlet[prim_idx.Pressure()] = Pressure;
    V_inlet[prim_idx.Density()] = Density;
    V_inlet[prim_idx.Enthalpy()] = Energy + Pressure / Density;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      V_inlet[iDim+prim_idx.Velocity()] = Velocity[iDim];

    /*--- Current solution at this boundary node. ---*/

    const auto* V_domain = nodes->GetPrimitive(iPoint);

    /*--- Normal vector for this vertex (negate for outward convention). ---*/

    su2double Normal[MAXNDIM] = {0.0};
    geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
    for (unsigned short iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

    /*--- Set various quantities in the solver class. ---*/

    conv_numerics->SetNormal(Normal);
    conv_numerics->SetPrimitive(V_domain, V_inlet);

    if (dynamic_grid)
      conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                geometry->nodes->GetGridVel(iPoint));

    /*--- Compute the residual using an upwind scheme. ---*/

    auto residual = conv_numerics->ComputeResidual(config);

    LinSysRes.AddBlock(iPoint, residual);

    /*--- Jacobian contribution for implicit integration. ---*/

    if (implicit)
      Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    /*--- Viscous contribution, omited to improve convergence. ---*/

  }
  END_SU2_OMP_FOR

}
```

The pattern:

- **`implicit`** is read from `config` at each call. `CNewtonIntegration` switches the time scheme to explicit for some residual evaluations. The BC must then not touch the Jacobian.
- **Loop**: `SU2_OMP_FOR_DYN(OMP_MIN_SIZE)` ... `END_SU2_OMP_FOR`. The integration opens a parallel region, and every thread of it calls the BC. The macro splits the vertices among the threads. Local variables are private to each thread.
- **Halo check**: `if (!geometry->nodes->GetDomain(iPoint)) continue;`.
- **Ghost state**: write it to `GetCharacPrimVar(val_marker, iVertex)`, with the primitive indices `prim_idx` (`Temperature()`, `Pressure()`, `Density()`, `Velocity()`, `Enthalpy()`, ...). `CFVMFlowSolverBase::Allocate` allocates this array for every vertex of every marker. The scalar solvers read it later in the same iteration to build their own ghost state (see [step 8](#8-turbulence-species-and-heat-solvers)). So fill all the primitives they need, including the viscosities if they use the diffusive boundary flux.
- **Normal**: the vertex normal points into the domain. Negate it for the flux.
- **Numerics**: call `SetNormal`, then `SetPrimitive` (domain state, ghost state). Add `SetSecondary` for non-ideal gases (see `BC_Sym_Plane`), and `SetGridVel` for moving grids. Then call `ComputeResidual(config)`.
- **Residual**: add the flux to the row of the boundary point only, with `LinSysRes.AddBlock(iPoint, residual)`. The ghost has no row.
- **Jacobian**: add only the diagonal block, with `Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i)`.

A **viscous** contribution uses `visc_numerics` in the same way (see `CEulerSolver::BC_Supersonic_Outlet`):

- the ghost coordinate is the coordinate of the normal neighbour, reflected;
- the gradients are those of the boundary point.

```cpp
       visc_numerics->SetNormal(Normal);
       su2double Coord_Reflected[MAXNDIM];
       unsigned long Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
       GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                                geometry->nodes->GetCoord(iPoint), Coord_Reflected);
       visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);

       /*--- Primitive variables, and gradient ---*/

       visc_numerics->SetPrimitive(V_domain, V_domain);
       visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));
```

For the **incompressible** solver, `CIncEulerSolver::BC_Far_Field` follows the same pattern. It also sets the entries that only the incompressible primitive vector has (`prim_idx.Beta()`, `prim_idx.CpTotal()`). It handles the pressure-based formulation separately.

**Per-vertex data.** Some BCs have data that varies along the marker, such as an inlet profile or values set from Python. Store such data per vertex, like `Inlet_Ttotal[val_marker][iVertex]`, allocated with `AllocVectorOfVectors(nVertex, ...)` in the solver. Each grid level has its own solver object, so the data exists on each level.

`CDriver::PreprocessInlet` fills the inlet data:

- from a profile file with `LoadInletProfile` (only `INLET_FLOW` and `SUPERSONIC_INLET`);
- otherwise with `SetUniformInlet` of each solver.

`CFVMFlowSolverBase::SetUniformInlet` copies the values of `MARKER_SUPERSONIC_INLET` into the inlet arrays. This is why `BC_Supersonic_Inlet` reads `Inlet_Ttotal`.

## 7. Strong (Dirichlet) and weak BCs

Choose how the BC acts on the equations:

- A **weak** BC adds a boundary flux to the residual (step 6).
- A **strong** (Dirichlet) BC imposes the solution values at the boundary points.

To impose a value strongly:

1. Set the value in `Solution_Old`, because the update adds the correction to it. Also set it in `Solution` if it is read before the update.
2. Zero the residual of the imposed equations.
3. Zero the multigrid truncation error of those equations (`SetRes_TruncErrorZero`, `SetVel_ResTruncError_Zero`).
4. If implicit, replace the rows of the Jacobian by identity rows with `Jacobian.DeleteValsRowi(iPoint, iVar)`. This function sets row `iVar` of the diagonal block to zero, with 1 on the diagonal. It also zeroes the same row in all off-diagonal blocks. The linear solver then returns a zero correction for that equation.

`CHeatSolver::BC_ConjugateHeat_Interface` (fluid side) is a complete strong BC for one equation:

```cpp
    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

      const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

      if (geometry->nodes->GetDomain(iPoint)) {
        const su2double T_Conjugate = GetConjugateHeatVariable(val_marker, iVertex, 0) / Temperature_Ref;

        nodes->SetSolution_Old(iPoint, &T_Conjugate);
        LinSysRes(iPoint, 0) = 0.0;
        nodes->SetRes_TruncErrorZero(iPoint);

        if (implicit) Jacobian.DeleteValsRowi(iPoint, 0);
      }
    }
    END_SU2_OMP_FOR
```

For a system, impose only some rows. The no-slip wall of the incompressible solver (`CIncNSSolver`) fixes the velocity rows and keeps the others:

```cpp
    if (dynamic_grid) {
      nodes->SetVelocity_Old(iPoint, geometry->nodes->GetGridVel(iPoint));
    } else {
      su2double zero[MAXNDIM] = {0.0};
      nodes->SetVelocity_Old(iPoint, zero);
    }

    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      LinSysRes(iPoint, iDim+1) = 0.0;
    nodes->SetVel_ResTruncError_Zero(iPoint);

    if (pressure_based) nodes->SetStrongBC(iPoint);

    /*--- Enforce the no-slip boundary condition in a strong way by
     modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/

    if (implicit) {
      for (unsigned short iVar = 1; iVar <= nDim; iVar++)
        Jacobian.DeleteValsRowi(iPoint, iVar);
    }
```

A third form modifies the residual and the Jacobian rows instead of replacing them. `CFVMFlowSolverBase::BC_Sym_Plane` removes the normal component of the momentum residual. It multiplies the Jacobian blocks in the row of the point by `(I - n n^T)`.

The pressure-based incompressible solver also needs `nodes->SetStrongBC(iPoint)` at strongly imposed points (see `CIncEulerSolver::BC_Far_Field` and `CIncNSSolver`).

## 8. Turbulence, species and heat solvers

Decide what each scalar solver does at the new boundary, and implement it.

These solvers derive from `CScalarSolver`. They get the same `switch` in `Space_Integration`, with `MainSolver` = `TURB_SOL`, `TRANS_SOL`, `SPECIES_SOL`, `HEAT_SOL`, ... For them, `conv_numerics` and `visc_numerics` are `nullptr`: `CDriver::InitializeNumerics` does not create `CONV_BOUND_TERM` or `VISC_BOUND_TERM` for scalar solvers. Their edge kernel computes their boundary fluxes. The edge kernel is the class that computes the flux of one edge (see [Adding a Solver](../Adding-a-Solver/#6-numerics)).

What each one does today if it does not override a BC:

| Function | Default for scalar solvers |
|---|---|
| any `BC_...` not overridden | `CSolver`'s empty body: no boundary flux |
| `BC_Sym_Plane` | `CScalarSolver`: empty, "Convective and viscous fluxes across symmetry plane are equal to zero." |
| `BC_Euler_Wall` | `CScalarSolver`: empty, `final` |
| `BC_Supersonic_Inlet`, `BC_Supersonic_Outlet` | `CScalarSolver`: call `BC_Inlet` and `BC_Outlet` (`final`) |
| `BC_Riemann`, `BC_Periodic` | `CScalarSolver`, `final` |
| `BC_Far_Field`, `BC_Fluid_Interface` | `CScalarSolver`, overridden by most models |

Overrides in the scalar solvers:

| Solver | `BC_...` it implements |
|---|---|
| `CTurbSASolver` | `Far_Field`, `Inlet`, `Outlet`, `HeatFlux_Wall`, `Isothermal_Wall`, `Engine_Inflow`, `Engine_Exhaust`, `ActDisk_Inlet`, `ActDisk_Outlet`, `Inlet_Turbo`, `Inlet_MixingPlane`, `Fluid_Interface` |
| `CTurbSSTSolver` | `Far_Field`, `Inlet`, `Outlet`, `HeatFlux_Wall`, `Isothermal_Wall`, `Inlet_Turbo`, `Inlet_MixingPlane`, `Fluid_Interface` |
| `CTransLMSolver` | `Far_Field`, `Inlet`, `Outlet`, `HeatFlux_Wall`, `Isothermal_Wall`, `Fluid_Interface` |
| `CSpeciesSolver` | `Far_Field`, `Inlet`, `Outlet`, `HeatFlux_Wall`, `Isothermal_Wall`, `Fluid_Interface` |
| `CHeatSolver` | `Far_Field`, `Inlet`, `Outlet`, `HeatFlux_Wall`, `Isothermal_Wall`, `ConjugateHeat_Interface`, `Fluid_Interface` |

So for a new flow BC, decide for each scalar solver:

- Is a zero flux right (walls, symmetry)?
- Or does the scalar need an inflow value or an outflow extrapolation?

Then override the function, or forward to an existing one as `CScalarSolver` does for the supersonic BCs.

A weak scalar BC has two passes:

1. fill the ghost state of every vertex;
2. call the boundary flux of the model.

`CTurbSSTSolver::BC_Outlet`:

```cpp
void CTurbSSTSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics*, CNumerics*,
                               CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  auto* flowSolver = solver_container[FLOW_SOL];

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Neumann: the turbulent variable is copied from the interior before computing the flux. ---*/
    for (auto iVar = 0u; iVar < nVar; iVar++) ghostNodes->SetSolution(iVertex, iVar, nodes->GetSolution(iPoint, iVar));

    SetGhostPrimitives(iVertex, flowSolver->GetCharacPrimVar(val_marker, iVertex));

    SetGhostGeometry(geometry, val_marker, iVertex);
  }
  END_SU2_OMP_FOR

  /*--- The diffusive term at this boundary causes serious convergence problems. ---*/
  BoundaryFlux(geometry, solver_container, config,
               ScalarFluxOptions::BoundaryConvective(*config, config->GetBounded_Turb()), val_marker);
}
```

- `ghostNodes` holds the scalar ghost state.
- `SetGhostPrimitives` copies the flow ghost state that the flow BC wrote to `GetCharacPrimVar` (density, velocity, viscosities).
- `SetGhostGeometry` stores the outward normal and marks the vertex as active (`ghostSkip`).
- `SetGhostDiffusionState` adds the reflected coordinate and gradient when the diffusive term is included.
- `BoundaryFlux` calls `CScalarSolver::BoundaryFluxResidual` with the model's kernel. It skips halo points and skipped vertices. It adds the flux to the row of the boundary point, and adds the diagonal Jacobian block.
- `ScalarFluxOptions::BoundaryConvective`, `BoundaryFull` and `BoundaryDiffusive` choose the terms.

Walls of the turbulence models are strong BCs. For example, `CTurbSSTSolver::BC_HeatFlux_Wall` sets `k` and `omega`, zeroes the residual and calls `Jacobian.DeleteValsRowi` for both rows. The heat solver writes its own wall fluxes (`CHeatSolver::BC_Isothermal_Wall`, `BC_HeatFlux_Wall`).

## 9. Other places that test the BC kind

Find and update every other place that tests the BC kind. The kind is tested in many places besides the dispatch. Before you finish, search the code for an existing kind that behaves like yours, for example `grep -rn "SUPERSONIC_INLET" Common SU2_CFD`. The main places:

| Where | What |
|---|---|
| `CConfig::GetViscous_Wall`, `CConfig::GetSolid_Wall` (`Common/src/CConfig.cpp`) | the list of no-slip walls (`HEAT_FLUX`, `ISOTHERMAL`, `HEAT_TRANSFER`, `SMOLUCHOWSKI_MAXWELL`, `CHT_WALL_INTERFACE`) and of solid walls (those plus `EULER_WALL`) |
| `CGeometry::ComputeWallDistance` | wall distance for the turbulence models, from `GetViscous_Wall` markers |
| `CPhysicalGeometry` | `SetSolidBoundary` / `SetViscousBoundary` of the points, from the two functions above |
| `CFVMFlowSolverBase::Pressure_Forces`, `Friction_Forces`; `CFlowOutput` | forces, skin friction, y+ and heat flux only on solid / viscous walls |
| `CFVMFlowSolverBase::Allocate` | `VertexTraction` (fluid-structure interaction) only on solid walls |
| `CMultiGridIntegration` | corrections on viscous walls (see below) |
| `CDriver::PreprocessInlet`, `CFVMFlowSolverBase::SetUniformInlet`, `CSpeciesSolver::SetUniformInlet`, `CSolver::LoadInletProfile` | inlet data for `INLET_FLOW` and `SUPERSONIC_INLET` |
| `CDriverBase::GetMarkerTypes` | the kind names given to the Python wrapper |

Add a new **wall** kind to `GetViscous_Wall` (no-slip) or `GetSolid_Wall` (slip). Otherwise the turbulence models do not see it in the wall distance, and the forces and wall outputs skip it.

## 10. Multigrid

Check that the flow BC works on the coarse multigrid levels.

The coarse grids (`CMultiGridGeometry`) are built by agglomeration: groups of fine-grid control volumes are merged into one coarse control volume. The coarse grids keep the markers of the fine grid. Boundary control volumes are agglomerated first. Points where two different BC kinds meet are kept apart (`FindMixedBoundaryNodes`). The flow BC is then called on every level with `iMesh`, through `CMultiGridIntegration::MultiGrid_Cycle`. Check that:

- the per-vertex data of the BC exists on each level (each level has its own solver);
- level-dependent parts are guarded with `iMesh == MESH_0`, as for `GILES_BOUNDARY`;
- a strong BC zeroes the truncation error (step 7). `CMultiGridIntegration` also treats some walls:
  - on `GetViscous_Wall` markers, it zeroes the velocity corrections (`GetProlongated_Correction`) and the velocity truncation error (`SetForcing_Term`). It also sets the wall velocity in the restricted solution (`SetRestricted_Solution`);
  - on `EULER_WALL` markers, it projects the corrections on the wall (`MultigridProjectEulerWall`);
- `CMultiGridIntegration::SmoothProlongated_Correction` restores the residuals at all boundary markers except `INTERNAL_BOUNDARY`, `NEARFIELD_BOUNDARY`, `PERIODIC_BOUNDARY` and `SEND_RECEIVE`.

The scalar solvers run on one grid only, `config->GetFinestMesh()` (`CSingleGridIntegration`).

## 11. Adjoint

- **Discrete adjoint.** The primal iteration, BCs included, is recorded on the AD tape and differentiated automatically. The tape is the record of operations that automatic differentiation runs backwards. Use `su2double` for all values that depend on the solution, and do not apply `SU2_TYPE::GetValue` to them. Values read from `config` are constants for the adjoint, unless they are registered as inputs. `CDiscAdjSolver::RegisterVariables` registers the free-stream Mach number, angle of attack, temperature and pressure for compressible flow. The data of a new BC is not registered.
- **Continuous adjoint.** `CAdjEulerSolver` and `CAdjNSSolver` have their own `BC_...` functions. They get the same dispatch, with `MainSolver = ADJFLOW_SOL`. Without an override, the empty default applies.

## 12. Output and Python

Any marker can be listed in `MARKER_MONITORING`, `MARKER_PLOTTING` and `MARKER_ANALYZE`. `CConfig::SetMarkers` and the geometry handle that for every kind. Values that only your BC computes (for example a mass flow) need output fields: see [Adding Output Fields](../Adding-Output-Fields/).

For Python:

- Add a name in `CDriverBase::GetMarkerTypes`. Otherwise the wrapper reports `"UNKNOWN_TYPE"`.
- Add the option to `read_config` in `SU2_PY/SU2/io/config.py`, but only if a script must read the marker list.

## 13. Test

1. Run a small case with the new marker: serial, with MPI, and with OpenMP if the solver supports it (`SU2_CFD -t <threads>`). Check a partition where the marker is split between ranks.
2. Check implicit and explicit time integration, and multigrid (`MGLEVEL`) for the flow solver.
3. Run it with a turbulence model, to see what the turbulence solver does at the new marker.
4. Add a regression test. For example, `MARKER_SUPERSONIC_INLET` is tested by `TestCases/euler/ramp/inv_ramp.cfg` (the test `ramp` in `TestCases/parallel_regression.py`). See [Running Regression Tests](../Running-Regression-Tests/#add-a-new-test).

## Checklist

1. `BC_TYPE` value in `Common/include/option_structure.hpp`.
2. Members, `MARKER_...` option, `nMarker_CfgFile` sum and loop in `CConfig::SetMarkers`, getters, `CConfig::SetOutput`, `config_template.cfg`.
3. `case` in the weak or strong `switch` of `CIntegration::Space_Integration`.
4. Virtual `BC_...` in `SU2_CFD/include/solvers/CSolver.hpp`.
5. Implementation in the flow solvers (compressible, incompressible, NEMO as needed), with halo check, OpenMP loop, `implicit` check and ghost state in `CharacPrimVar`.
6. Decision and implementation for each scalar solver (turbulence, transition, species, heat).
7. `GetViscous_Wall` / `GetSolid_Wall` for walls; inlet data if it is an inlet.
8. Multigrid, adjoint and Python checks.
9. Regression test.

Last checked against SU2 develop `11c177b904`.
