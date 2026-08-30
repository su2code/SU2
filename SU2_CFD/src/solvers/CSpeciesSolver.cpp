/*!
 * \file CSpeciesSolver.cpp
 * \brief Main subroutines of CSpeciesSolver class
 * \author T. Kattmann
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/solvers/CSpeciesSolver.hpp"

#include "../../../Common/include/option_structure.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/solvers/CScalarSolver.inl"
#include "../../include/numerics/species/species_edge_flux.hpp"

/*--- Explicit instantiation of the parent class of CSpeciesSolver. ---*/
template class CScalarSolver<CSpeciesVariable>;

CSpeciesSolver::CSpeciesSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh)
    : CScalarSolver<CSpeciesVariable>(geometry, config, true, config->GetBounded_Species()) {
  SU2_ZONE_SCOPED

  /*--- Dimension of the problem. ---*/

  nVar = config->GetnSpecies();

  Initialize(geometry, config, iMesh, nVar);

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CSpeciesVariable(Solution_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- Ghost states for boundary conditions, sized to the largest marker (see BoundaryFluxResidual). ---*/
  unsigned long maxMarkerVertices = 0;
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++)
    maxMarkerVertices = max(maxMarkerVertices, nVertex[iMarker]);
  ghostNodes = make_unique<CSpeciesVariable>(Solution_Inf, maxMarkerVertices, nDim, nVar, config);

  /*--- Initialize the mass diffusivity. Nondimensionalization done in the flow solver. ---*/
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    for (auto iVar = 0u; iVar <= nVar; iVar++) {
      const auto MassDiffusivity = config->GetDiffusivity_ConstantND();
      nodes->SetDiffusivity(iPoint, MassDiffusivity, iVar);
    }
  }
  END_SU2_OMP_FOR

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

  SlidingState.resize(nMarker);
  SlidingStateNodes.resize(nMarker);

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
      SlidingState[iMarker].resize(nVertex[iMarker], nPrimVar+1) = nullptr;
      SlidingStateNodes[iMarker].resize(nVertex[iMarker],0);
    }
  }

  /*--- Store the initial CFL number for all grid points. ---*/

  const su2double CFL = config->GetCFL(MGLevel) * config->GetCFLRedCoeff_Species();
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0u; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  END_SU2_OMP_FOR
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name. ---*/
  SolverName = "SPECIES";
}


void CSpeciesSolver::Initialize(CGeometry* geometry, CConfig* config, unsigned short iMesh, unsigned short nVar) {
  SU2_ZONE_SCOPED
  /*--- Store if an implicit scheme is used, for use during periodic boundary conditions. ---*/
  SetImplicitPeriodic(config->GetKind_TimeIntScheme_Species() == EULER_IMPLICIT);

  nPrimVar = nVar;

  if (nVar > MAXNVAR)
    SU2_MPI::Error("Increase static array size MAXNVAR for CSpeciesVariable and proceed.", CURRENT_FUNCTION);

  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  AllocVectorOfMatrices( nVertex, nVar,CustomBoundaryScalar);

  if (iMesh == MESH_0 || config->GetMGCycle() == MG_CYCLE::FULL) {

    /*--- Define some auxiliary vector related with the residual ---*/

    Residual_RMS.resize(nVar, 0.0);
    Residual_Max.resize(nVar, 0.0);
    Point_Max.resize(nVar, 0);
    Point_Max_Coord.resize(nVar, nDim) = su2double(0.0);

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (config->GetMultizone_Problem()) {
      Residual_BGS.resize(nVar, 0.0);
      Residual_Max_BGS.resize(nVar, 0.0);
      Point_Max_BGS.resize(nVar, 0);
      Point_Max_Coord_BGS.resize(nVar, nDim) = su2double(0.0);
    }

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (species transport model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    System.SetxIsZero(true);

    if (ReducerStrategy) {
      EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
      EdgeFluxesDiff.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
    }
  }

  /*--- Initialize lower and upper limits---*/

  if (config->GetSpecies_Clipping()) {
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      lowerlimit[iVar] = config->GetSpecies_Clipping_Min(iVar);
      upperlimit[iVar] = config->GetSpecies_Clipping_Max(iVar);
    }
  } else {
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      lowerlimit[iVar] = -1.0e15;
      upperlimit[iVar] = 1.0e15;
    }
  }

  /*--- Scalar variable state at the far-field. ---*/

  for (auto iVar = 0u; iVar < nVar; iVar++) {
    Solution_Inf[iVar] = config->GetSpecies_Init()[iVar];
  }

  /*--- Set the column number for species in inlet-files.
   * e.g. Coords(nDim), Temp(1), VelMag(1), Normal(nDim), Turb(1 or 2), Species(arbitrary) ---*/
  Inlet_Position = nDim + 2 + nDim + config->GetnTurbVar();

  /*-- Allocation of inlet-values. Will be filled either by an inlet file,
   * or uniformly by a uniform boundary condition. ---*/

  Inlet_SpeciesVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_SpeciesVars[iMarker].resize(nVertex[iMarker], nVar);
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        Inlet_SpeciesVars[iMarker](iVertex, iVar) = Solution_Inf[iVar];
      }
    }
  }

  Wall_SpeciesVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Wall_SpeciesVars[iMarker].resize(nVertex[iMarker], nVar);
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        Wall_SpeciesVars[iMarker](iVertex, iVar) = Solution_Inf[iVar];
      }
    }
  }
}


void CSpeciesSolver::LoadRestart(CGeometry** geometry, CSolver*** solver, CConfig* config, int val_iter,
                                 bool val_update_geo) {
  SU2_ZONE_SCOPED
  /*--- Restart the solution from file information ---*/

  string restart_filename = config->GetSolution_FileName();

  /*--- To make this routine safe to call in parallel most of it can only be executed by one thread. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

    if (config->GetRead_Binary_Restart()) {
      restart_filename = config->GetFilename(restart_filename, ".dat", val_iter);
      Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
    } else {
      restart_filename = config->GetFilename(restart_filename, ".csv", val_iter);
      Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
    }

    /*--- Skip flow variables and turbulence variables. ---*/

    unsigned short skipVars = nDim + solver[MESH_0][FLOW_SOL]->GetnVar() + config->GetnTurbVar();

    /*--- Adjust the number of solution variables in the incompressible
     restart. We always carry a space in nVar for the energy equation in the
     mean flow solver, but we only write it to the restart if it is active.
     Therefore, we must reduce skipVars here if energy is inactive so that
     the turbulent variables are read correctly. ---*/

    const bool incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
    const bool energy = config->GetEnergy_Equation();
    const bool weakly_coupled_heat = config->GetWeakly_Coupled_Heat();

    /*--- for the flamelet model, the temperature is saved to file, but the energy equation is off ---*/

   if (incompressible && ((!energy) && (!weakly_coupled_heat))) skipVars--;

    /*--- Load data from the restart into correct containers. ---*/

    unsigned long counter = 0;
    for (auto iPoint_Global = 0ul; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {
      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      const auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {
        /*--- We need to store this point's data, so jump to the correct
         offset in the buffer of data from the restart file and load it. ---*/

        const auto index = counter * Restart_Vars[1] + skipVars;
        for (auto iVar = 0u; iVar < nVar; iVar++)
          nodes->SetSolution(iPoint_Local, iVar, Restart_Data[index + iVar]);

        /*--- Increment the overall counter for how many points have been loaded. ---*/
        counter++;
      }
    }

    /*--- Detect a wrong solution file ---*/

    if (counter != nPointDomain) {
      SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                         string("It could be empty lines at the end of the file."),
                     CURRENT_FUNCTION);
    }

  }  // end safe global access, pre and postprocessing are thread-safe.
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- MPI solution and compute the eddy viscosity ---*/

  solver[MESH_0][SPECIES_SOL]->InitiateComms(geometry[MESH_0], config, MPI_QUANTITIES::SOLUTION);
  solver[MESH_0][SPECIES_SOL]->CompleteComms(geometry[MESH_0], config, MPI_QUANTITIES::SOLUTION);

  // Flow-Pre computes/sets mixture properties
  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER,
                                          RUNTIME_FLOW_SYS, true);
  // Update eddy-visc which needs correct mixture density and mixture lam-visc. Note that after this, another Flow-Pre
  // at the start of the Iteration sets the updated eddy-visc into the Flow-Solvers Primitives.
  if (config->GetKind_Turb_Model() != TURB_MODEL::NONE)
    solver[MESH_0][TURB_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
  // For feature_multicomp this Scalar-Pre only computes the laminar contribution to mass diffusivity
  solver[MESH_0][SPECIES_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER,
                                             RUNTIME_SPECIES_SYS, false);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (auto iMesh = 1u; iMesh <= config->GetnMGLevels(); iMesh++) {
    MultigridRestriction(*geometry[iMesh - 1], solver[iMesh - 1][SPECIES_SOL]->GetNodes()->GetSolution(),
                         *geometry[iMesh], solver[iMesh][SPECIES_SOL]->GetNodes()->GetSolution());
    solver[iMesh][SPECIES_SOL]->InitiateComms(geometry[iMesh], config, MPI_QUANTITIES::SOLUTION);
    solver[iMesh][SPECIES_SOL]->CompleteComms(geometry[iMesh], config, MPI_QUANTITIES::SOLUTION);

    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS,
                                           true);

    if (config->GetKind_Turb_Model() != TURB_MODEL::NONE)
      solver[iMesh][TURB_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);

    solver[iMesh][SPECIES_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER,
                                              RUNTIME_SPECIES_SYS, false);
  }

  /*--- Go back to single threaded execution. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Delete the class memory that is used to load the restart. ---*/

    Restart_Vars = decltype(Restart_Vars){};
    Restart_Data = decltype(Restart_Data){};
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

void CSpeciesSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                   unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem,
                                   bool Output) {
  SU2_ZONE_SCOPED
  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0u; iPoint < nPoint; iPoint++) {
    const su2double temperature = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
    const su2double* scalar = solver_container[SPECIES_SOL]->GetNodes()->GetSolution(iPoint);
    solver_container[FLOW_SOL]->GetFluidModel()->SetMassDiffusivityModel(config);
    solver_container[FLOW_SOL]->GetFluidModel()->SetTDState_T(temperature, scalar);
    /*--- Recompute viscosity, important  to get diffusivity correct across MPI ranks. ---*/
    nodes->SetLaminarViscosity(iPoint, solver_container[FLOW_SOL]->GetFluidModel()->GetLaminarViscosity());
    /*--- Set the laminar mass Diffusivity for the species solver. ---*/
    for (auto iVar = 0u; iVar <= nVar; iVar++) {
      const su2double mass_diffusivity = solver_container[FLOW_SOL]->GetFluidModel()->GetMassDiffusivity(iVar);
      nodes->SetDiffusivity(iPoint, mass_diffusivity, iVar);
    }
  }  // iPoint
  END_SU2_OMP_FOR

  /*--- Clear Residual and Jacobian. Upwind second order reconstruction and gradients. ---*/
  CommonPreprocessing(geometry, config, Output);

  EnsureGhostFlowContainers(solver_container, config);
}

void CSpeciesSolver::Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics**,
                                     CConfig* config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  /*--- The diffusion coefficient does not depend on the species mass fractions, so there is no
   * accurate-Jacobian correction to apply. ---*/
  const auto opt = ScalarFluxOptions::Interior(*config, config->GetBounded_Species());

  DispatchScheme<CScalarFlux_Species, Dynamic>(config, [&](auto tag) {
    EdgeFluxResidual<typename decltype(tag)::type>(geometry, solver_container, config, opt);
  });
}

void CSpeciesSolver::BoundaryFlux(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                  const ScalarFluxOptions& opt, unsigned short val_marker) {
  DispatchScheme<CScalarFlux_Species, Dynamic>(config, [&](auto tag) {
    BoundaryFluxResidual<typename decltype(tag)::type>(geometry, solver_container, config, opt, val_marker);
  });
}

void CSpeciesSolver::BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics*, CNumerics*, CConfig* config,
                              unsigned short val_marker) {
  SU2_ZONE_SCOPED

  const string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  if (config->GetMarker_StrongBC(Marker_Tag)) {
    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
      const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

      if (geometry->nodes->GetDomain(iPoint)) {
        nodes->SetSolution_Old(iPoint, Inlet_SpeciesVars[val_marker][iVertex]);

        LinSysRes.SetBlock_Zero(iPoint);

        /*--- Includes 1 in the diagonal ---*/
        for (auto iVar = 0u; iVar < nVar; iVar++) {
          Jacobian.DeleteValsRowi(iPoint, iVar);
        }
      }
    }
    END_SU2_OMP_FOR
    return;
  }

  /*--- Weak BC: fill the ghost row from the inlet species state, then let the edge kernel
   * compute the purely convective flux. ---*/

  auto* flowSolver = solver_container[FLOW_SOL];

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    for (auto iVar = 0u; iVar < nVar; iVar++)
      ghostNodes->SetSolution(iVertex, iVar, Inlet_SpeciesVars[val_marker][iVertex][iVar]);

    SetGhostPrimitives(iVertex, flowSolver->GetCharacPrimVar(val_marker, iVertex));

    SetGhostGeometry(geometry, val_marker, iVertex);
  }
  END_SU2_OMP_FOR

  /*--- This site imposes a convective flux alone; a diffusive one has not been validated. ---*/
  BoundaryFlux(geometry, solver_container, config,
               ScalarFluxOptions::BoundaryConvective(*config, config->GetBounded_Species()), val_marker);
}



void CSpeciesSolver::BC_Isothermal_Wall(CGeometry* geometry, CSolver** solver_container,
                                                CNumerics* conv_numerics, CNumerics* visc_numerics,
                                                CConfig* config, unsigned short val_marker) {
  SU2_ZONE_SCOPED
  BC_Wall_Generic(geometry, solver_container, config, val_marker);
}

void CSpeciesSolver::BC_HeatFlux_Wall(CGeometry* geometry, CSolver** solver_container,
                                                CNumerics* conv_numerics, CNumerics* visc_numerics,
                                                CConfig* config, unsigned short val_marker) {
  SU2_ZONE_SCOPED
  BC_Wall_Generic(geometry, solver_container, config, val_marker);
}

void CSpeciesSolver::BC_Wall_Generic(CGeometry* geometry, CSolver** solver_container,
                                                        CConfig* config, unsigned short val_marker) {
  SU2_ZONE_SCOPED
  const bool implicit = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  const bool py_custom = config->GetMarker_All_PyCustom(val_marker);

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  for (auto iVar = 0u; iVar < nVar; iVar++) {

    // Get wall species boundary condition type and value for this marker and species
    const su2double WallSpeciesValue = config->GetWall_SpeciesVal(Marker_Tag, iVar);
    const WALL_SPECIES_TYPE wallspeciestype = config->GetWall_SpeciesType(Marker_Tag, iVar);

    SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

      const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

      if (!geometry->nodes->GetDomain(iPoint)) continue;

      const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

      su2double Area = GeometryToolbox::Norm(nDim, Normal);

      su2double WallSpecies = WallSpeciesValue;

      /*--- Get the scalar values from the python wrapper. ---*/
      if (py_custom) {
        WallSpecies = CustomBoundaryScalar[val_marker](iVertex,iVar);
      }

      switch(wallspeciestype) {
        case WALL_SPECIES_TYPE::FLUX:
          //Flux Boundary condition
          LinSysRes(iPoint, iVar) -= WallSpecies * Area;
        break;
        case WALL_SPECIES_TYPE::VALUE:
          //Dirichlet Strong Boundary Condition
          nodes->SetSolution(iPoint, iVar, WallSpecies);
          nodes->SetSolution_Old(iPoint, iVar, WallSpecies);
          LinSysRes(iPoint, iVar) = 0.0;
          if (implicit) {
            Jacobian.DeleteValsRowi(iPoint, iVar);
          }
        break;
      }
    }
    END_SU2_OMP_FOR
  }
}



void CSpeciesSolver::SetInletAtVertex(const su2double *val_inlet,
                                      unsigned short iMarker,
                                      unsigned long iVertex) {
  SU2_ZONE_SCOPED

  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Inlet_SpeciesVars[iMarker][iVertex][iVar] = val_inlet[Inlet_Position+iVar];

}

su2double CSpeciesSolver::GetInletAtVertex(unsigned short iMarker, unsigned long iVertex,
                                           const CGeometry* geometry, su2double* val_inlet) const {
  SU2_ZONE_SCOPED
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    val_inlet[Inlet_Position + iVar] = Inlet_SpeciesVars[iMarker][iVertex][iVar];

  /*--- Compute boundary face area for this vertex. ---*/

  su2double Normal[MAXNDIM] = {0.0};
  geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
  return GeometryToolbox::Norm(nDim, Normal);
}

void CSpeciesSolver::SetUniformInlet(const CConfig* config, unsigned short iMarker) {
  SU2_ZONE_SCOPED
  bool riemann_inlet = false;

  const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
  if (config->GetMarker_All_KindBC(iMarker) == RIEMANN_BOUNDARY) {
    switch (config->GetKind_Data_Riemann(Marker_Tag)) {
      case TOTAL_CONDITIONS_PT: case STATIC_SUPERSONIC_INFLOW_PT: case STATIC_SUPERSONIC_INFLOW_PD: case DENSITY_VELOCITY:
      riemann_inlet = true;
      break;
    }
  }
  /*--- Find BC string to the numeric-identifier. ---*/
  if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW || config->GetMarker_All_KindBC(iMarker) == SUPERSONIC_INLET || riemann_inlet) {
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        Inlet_SpeciesVars[iMarker][iVertex][iVar] = config->GetInlet_SpeciesVal(Marker_Tag)[iVar];
      }
    }
  }
}

void CSpeciesSolver::BC_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics*, CNumerics*,
                               CConfig* config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  const string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  if (config->GetMarker_StrongBC(Marker_Tag)) {
    /*--- Strong zero flux Neumann boundary condition at the outlet ---*/
    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
      const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

      if (geometry->nodes->GetDomain(iPoint)) {
        const auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        nodes->SetSolution_Old(iPoint, nodes->GetSolution(Point_Normal));

        LinSysRes.SetBlock_Zero(iPoint);

        /*--- Includes 1 on the diagonal ---*/
        for (auto iVar = 0u; iVar < nVar; iVar++) {
          Jacobian.DeleteValsRowi(iPoint, iVar);
        }
      }
    }
    END_SU2_OMP_FOR
    return;
  }

  /*--- Weak BC: Neumann, the species variable is copied from the interior of the domain to the
   * ghost row before the edge kernel computes the purely convective flux. ---*/

  auto* flowSolver = solver_container[FLOW_SOL];

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    for (auto iVar = 0u; iVar < nVar; iVar++) ghostNodes->SetSolution(iVertex, iVar, nodes->GetSolution(iPoint, iVar));

    SetGhostPrimitives(iVertex, flowSolver->GetCharacPrimVar(val_marker, iVertex));

    SetGhostGeometry(geometry, val_marker, iVertex);
  }
  END_SU2_OMP_FOR

  /*--- This site imposes a convective flux alone; a diffusive one has not been validated. ---*/
  BoundaryFlux(geometry, solver_container, config,
               ScalarFluxOptions::BoundaryConvective(*config, config->GetBounded_Species()), val_marker);
}

void CSpeciesSolver::BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics*, CNumerics*,
                                  CConfig* config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  auto* flowSolver = solver_container[FLOW_SOL];

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    for (auto iVar = 0u; iVar < nVar; iVar++) ghostNodes->SetSolution(iVertex, iVar, Solution_Inf[iVar]);

    SetGhostPrimitives(iVertex, flowSolver->GetCharacPrimVar(val_marker, iVertex));

    SetGhostGeometry(geometry, val_marker, iVertex);
  }
  END_SU2_OMP_FOR

  BoundaryFlux(geometry, solver_container, config,
               ScalarFluxOptions::BoundaryConvective(*config, config->GetBounded_Species()), val_marker);
}

void CSpeciesSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics*,
                                        CNumerics*, CConfig *config) {
  SU2_ZONE_SCOPED

  if (solver_container[FLOW_SOL] == nullptr) return;

  const auto optConv = ScalarFluxOptions::BoundaryConvective(*config, config->GetBounded_Species());
  const auto optVisc = ScalarFluxOptions::BoundaryDiffusive(*config, false);

  /*--- The diffusion coefficients read the interior point's own diffusivity on both sides of the
   * edge, so the ghost row takes it too. ---*/
  const auto fillGhostExtras = [&](unsigned long iVertex, unsigned long iPoint) {
    for (auto iVar = 0u; iVar < nVar; iVar++)
      ghostNodes->SetDiffusivity(iVertex, nodes->GetDiffusivity(iPoint, iVar), iVar);
  };

  DispatchScheme<CScalarFlux_Species, Dynamic>(config, [&](auto tag) {
    FluidInterfaceFluxResidual<typename decltype(tag)::type>(geometry, solver_container, config, optConv, optVisc,
                                                             fillGhostExtras);
  });
}

void CSpeciesSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                      CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool axisymmetric = config->GetAxisymmetric();

  if (axisymmetric) {
    CNumerics *numerics  = numerics_container[SOURCE_FIRST_TERM  + omp_get_thread_num()*MAX_TERMS];

    SU2_OMP_FOR_DYN(omp_chunk_size)
    for (auto iPoint = 0u; iPoint < nPointDomain; iPoint++) {
      /*--- Set primitive variables w/o reconstruction ---*/

      numerics->SetPrimitive(solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint), nullptr);

      /*--- Set scalar variables w/o reconstruction ---*/

      numerics->SetScalarVar(nodes->GetSolution(iPoint), nullptr);

      numerics->SetDiffusionCoeff(nodes->GetDiffusivity(iPoint), nullptr);

      /*--- Set volume of the dual cell. ---*/

      numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Update scalar sources in the fluidmodel ---*/

      /*--- Axisymmetry source term for the scalar equation. ---*/
      /*--- Set y coordinate ---*/

      numerics->SetCoord(geometry->nodes->GetCoord(iPoint), nullptr);

      /*--- Set gradients ---*/

      numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nullptr);

      auto residual = numerics->ComputeResidual(config);

      /*--- Add Residual ---*/

      LinSysRes.SubtractBlock(iPoint, residual);

      /*--- Implicit part ---*/

      if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

  /*--- Custom user defined source term (from the python wrapper) ---*/
  if (config->GetPyCustomSource()) {
    CustomSourceResidual(geometry, solver_container, numerics_container, config, iMesh);
  }

}

void CSpeciesSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {
  SU2_ZONE_SCOPED

  const bool restart = config->GetRestart() || config->GetRestart_Flow();

  PushSolutionBackInTime(TimeIter, restart, solver_container, geometry, config);
}
