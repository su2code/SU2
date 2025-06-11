/*!
 * \file CSpeciesSolver.cpp
 * \brief Main subroutines of CSpeciesSolver class
 * \author T. Kattmann
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/solvers/CScalarSolver.inl"

/*--- Explicit instantiation of the parent class of CSpeciesSolver. ---*/
template class CScalarSolver<CSpeciesVariable>;

CSpeciesSolver::CSpeciesSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh)
    : CScalarSolver<CSpeciesVariable>(geometry, config, true) {

  /*--- Dimension of the problem. ---*/

  nVar = config->GetnSpecies();

  Initialize(geometry, config, iMesh, nVar);

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CSpeciesVariable(Solution_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

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

  SpeciesPointSource.resize(nPointDomain,nVar);
  SpeciesPointSource.setConstant(0.0);


  if (iMesh == MESH_0 || config->GetMGCycle() == FULLMG_CYCLE) {

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

    if (ReducerStrategy) EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
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
        cout << "Solution=" << Solution_Inf[iVar] << endl;
      }
    }
  }
}


void CSpeciesSolver::LoadRestart(CGeometry** geometry, CSolver*** solver, CConfig* config, int val_iter,
                                 bool val_update_geo) {
  /*--- Restart the solution from file information ---*/
  const string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  /*--- To make this routine safe to call in parallel most of it can only be executed by one thread. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

    if (config->GetRead_Binary_Restart()) {
      Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
    } else {
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
    const bool flamelet = (config->GetKind_FluidModel() == FLUID_FLAMELET);
    const bool weakly_coupled_heat = config->GetWeakly_Coupled_Heat();

    /*--- for the flamelet model, the temperature is saved to file, but the energy equation is off ---*/

   if (incompressible && ((!energy) && (!weakly_coupled_heat) && (!flamelet))) skipVars--;

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
                                          RUNTIME_FLOW_SYS, false);
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
                                           false);

    if (config->GetKind_Turb_Model() != TURB_MODEL::NONE)
      solver[iMesh][TURB_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);

    solver[iMesh][SPECIES_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER,
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
  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)

  /*--- Set the laminar mass Diffusivity for the species solver. ---*/
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0u; iPoint < nPoint; iPoint++) {
    const su2double temperature = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
    const su2double* scalar = solver_container[SPECIES_SOL]->GetNodes()->GetSolution(iPoint);
    solver_container[FLOW_SOL]->GetFluidModel()->SetMassDiffusivityModel(config);
    solver_container[FLOW_SOL]->GetFluidModel()->SetTDState_T(temperature, scalar);
    for (auto iVar = 0u; iVar <= nVar; iVar++) {
      const su2double mass_diffusivity = solver_container[FLOW_SOL]->GetFluidModel()->GetMassDiffusivity(iVar);
      nodes->SetDiffusivity(iPoint, mass_diffusivity, iVar);
    }

  }  // iPoint
  END_SU2_OMP_FOR

  /*--- Clear Residual and Jacobian. Upwind second order reconstruction and gradients ---*/
  CommonPreprocessing(geometry, config, Output);
}

void CSpeciesSolver::Viscous_Residual(const unsigned long iEdge, const CGeometry* geometry, CSolver** solver_container,
                                      CNumerics* numerics, const CConfig* config) {
  /*--- Define an object to set solver specific numerics contribution. ---*/
  auto SolverSpecificNumerics = [&](unsigned long iPoint, unsigned long jPoint) {
    /*--- Mass diffusivity coefficients. ---*/

    numerics->SetDiffusionCoeff(nodes->GetDiffusivity(iPoint), nodes->GetDiffusivity(jPoint));
  };

  /*--- Now instantiate the generic implementation with the functor above. ---*/

  Viscous_Residual_impl(SolverSpecificNumerics, iEdge, geometry, solver_container, numerics, config);
}

void CSpeciesSolver::BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                              CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- Identify the boundary by string name ---*/
    string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

    if (config->GetMarker_StrongBC(Marker_Tag)==true) {
      nodes->SetSolution_Old(iPoint, Inlet_SpeciesVars[val_marker][iVertex]);

      LinSysRes.SetBlock_Zero(iPoint);

      /*--- Includes 1 in the diagonal ---*/
      for (auto iVar = 0u; iVar < nVar; iVar++) {
        auto total_index = iPoint * nVar + iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    } else {  // weak BC
      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++) Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/

      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the species variable state at the inlet. ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Inlet_SpeciesVars[val_marker][iVertex]);

      /*--- Set various other quantities in the solver class ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      if (conv_numerics->GetBoundedScalar()) {
        const su2double* velocity = &V_inlet[prim_idx.Velocity()];
        const su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        conv_numerics->SetMassFlux(BoundedScalarBCFlux(iPoint, implicit, density, velocity, Normal));
      }

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      // Unfinished viscous contribution removed before right after d8a0da9a00. Further testing required.

    }
  }
  END_SU2_OMP_FOR
}

void CSpeciesSolver::SetInletAtVertex(const su2double *val_inlet,
                                      unsigned short iMarker,
                                      unsigned long iVertex) {

  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Inlet_SpeciesVars[iMarker][iVertex][iVar] = val_inlet[Inlet_Position+iVar];

}

su2double CSpeciesSolver::GetInletAtVertex(unsigned short iMarker, unsigned long iVertex,
                                           const CGeometry* geometry, su2double* val_inlet) const {
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    val_inlet[Inlet_Position + iVar] = Inlet_SpeciesVars[iMarker][iVertex][iVar];

  /*--- Compute boundary face area for this vertex. ---*/

  su2double Normal[MAXNDIM] = {0.0};
  geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
  return GeometryToolbox::Norm(nDim, Normal);
}

void CSpeciesSolver::SetUniformInlet(const CConfig* config, unsigned short iMarker) {
  /*--- Find BC string to the numeric-identifier. ---*/
  if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
    const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        Inlet_SpeciesVars[iMarker][iVertex][iVar] = config->GetInlet_SpeciesVal(Marker_Tag)[iVar];
      }
    }
  }
}

void CSpeciesSolver::BC_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                               CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    /*--- Strong zero flux Neumann boundary condition at the outlet ---*/
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- Identify the boundary by string name ---*/
    string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

    if (config->GetMarker_StrongBC(Marker_Tag)==true) {
      /*--- Allocate the value at the outlet ---*/
      auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      nodes->SetSolution_Old(iPoint, nodes->GetSolution(Point_Normal));

      LinSysRes.SetBlock_Zero(iPoint);

      /*--- Includes 1 on the diagonal ---*/
      for (auto iVar = 0u; iVar < nVar; iVar++) {
        auto total_index = iPoint * nVar + iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    } else {  // weak BC

      /*--- Allocate the value at the outlet ---*/
      auto V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_outlet);

      /*--- Set the species variables. Here we use a Neumann BC such
      that the species variable is copied from the interior of the
      domain to the outlet before computing the residual. ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), nodes->GetSolution(iPoint));

      /*--- Set Normal (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++) Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      if (conv_numerics->GetBoundedScalar()) {
        const su2double* velocity = &V_outlet[prim_idx.Velocity()];
        const su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        conv_numerics->SetMassFlux(BoundedScalarBCFlux(iPoint, implicit, density, velocity, Normal));
      }

      /*--- Compute the residual using an upwind scheme ---*/
      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      // Unfinished viscous contribution removed before right after d8a0da9a00. Further testing required.

    }
  }
  END_SU2_OMP_FOR
}

void CSpeciesSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                      CConfig *config, unsigned short iMesh) {

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
}


void CSpeciesSolver::Custom_Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                      CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[SOURCE_SECOND_TERM + omp_get_thread_num()*MAX_TERMS];

  unsigned short iVar;
  unsigned long iPoint;
  AD::StartNoSharedReading();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Load the volume of the dual mesh cell ---*/
    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Get control volume size. ---*/
      su2double Volume = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the residual for this control volume and subtract. ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        LinSysRes[iPoint*nVar+iVar] -= SpeciesPointSource[iPoint][iVar] * Volume;
      }
  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

}


void CSpeciesSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container,
  CNumerics *conv_numerics, CNumerics *visc_numerics,
  CConfig *config, unsigned short val_marker) {

const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
const su2double Temperature_Ref = config->GetTemperature_Ref();
const su2double Prandtl_Lam = config->GetPrandtl_Lam();
const su2double Prandtl_Turb = config->GetPrandtl_Turb();
const su2double Gas_Constant = config->GetGas_ConstantND();
const su2double Gamma = config->GetGamma();
const su2double Cp = (Gamma / (Gamma-1)) * Gas_Constant;

/*--- Identify the boundary and retrieve the specified wall temperature from
the config (for non-CHT problems) as well as the wall function treatment. ---*/

const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);
su2double Twall = 0.0;
Twall = config->GetIsothermal_Temperature(Marker_Tag) / Temperature_Ref;


//  Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
//  if (Wall_Function != WALL_FUNCTION::NONE) {
//    SU2_MPI::Error("Wall function treatment not implemented yet", CURRENT_FUNCTION);
//  }

su2double **Jacobian_i = nullptr;
if (implicit) {
Jacobian_i = new su2double* [nVar];
for (auto iVar = 0u; iVar < nVar; iVar++)
Jacobian_i[iVar] = new su2double [nVar] ();
}

/*--- Loop over boundary points ---*/

SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

if (!geometry->nodes->GetDomain(iPoint)) continue;

/*--- Compute dual-grid area and boundary normal ---*/

const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

su2double Area = GeometryToolbox::Norm(nDim, Normal);

su2double UnitNormal[MAXNDIM] = {0.0};
for (auto iDim = 0u; iDim < nDim; iDim++)
UnitNormal[iDim] = -Normal[iDim]/Area;

/*--- Compute closest normal neighbor ---*/

const auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

/*--- Get coordinates of i & nearest normal and compute distance ---*/

const auto Coord_i = geometry->nodes->GetCoord(iPoint);
const auto Coord_j = geometry->nodes->GetCoord(Point_Normal);
const su2double dist_ij = GeometryToolbox::NormalDistance(nDim, UnitNormal, Coord_i, Coord_j);

su2double zero[MAXNDIM] = {0.0};

LinSysRes(iPoint, 1) = 0.0;


/*--- Get transport coefficients ---*/

su2double laminar_viscosity    = nodes->GetLaminarViscosity(iPoint);
su2double eddy_viscosity       = nodes->GetEddyViscosity(iPoint);
su2double thermal_conductivity = nodes->GetThermalConductivity(iPoint);  // * (laminar_viscosity/Prandtl_Lam + eddy_viscosity/Prandtl_Turb);

// work in progress on real-gases...
//thermal_conductivity = nodes->GetThermalConductivity(iPoint);
//Cp = nodes->GetSpecificHeatCp(iPoint);
//thermal_conductivity += Cp*eddy_viscosity/Prandtl_Turb;

/*--- If it is a customizable or CHT patch, retrieve the specified wall temperature. ---*/

const su2double There = nodes->GetTemperature(Point_Normal);



/*--- Compute the normal gradient in temperature using Twall ---*/

su2double dTdn = -(There - Twall)/dist_ij;

/*--- Apply a weak boundary condition for the energy equation.
Compute the residual due to the prescribed heat flux. ---*/

su2double Res_Conv = 0.0;
su2double Res_Visc = thermal_conductivity * dTdn * Area;

/*--- Calculate Jacobian for implicit time stepping ---*/

if (implicit) {

/*--- Add contributions to the Jacobian from the weak enforcement of the energy equations. ---*/

su2double Density = nodes->GetDensity(iPoint);
//su2double Vel2 = GeometryToolbox::SquaredNorm(nDim, &nodes->GetPrimitive(iPoint)[prim_idx.Velocity()]);
//su2double dTdrho = 1.0/Density * ( -Twall + (Gamma-1.0)/Gas_Constant*(Vel2/2.0) );

// Jacobian_i[1][0] = thermal_conductivity/dist_ij * dTdrho * Area;

// for (auto jDim = 0u; jDim < nDim; jDim++)
// Jacobian_i[1][jDim+1] = 0.0;

Jacobian_i[1][1] = thermal_conductivity/dist_ij * (Gamma-1.0)/(Gas_Constant*Density) * Area;
}



/*--- Convective and viscous contributions to the residual at the wall ---*/

LinSysRes(iPoint, 1) += Res_Conv - Res_Visc;

/*--- Enforce the no-slip boundary condition in a strong way by
modifying the velocity-rows of the Jacobian (1 on the diagonal).
And add the contributions to the Jacobian due to energy. ---*/

if (implicit) {
Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

for (auto iVar = 1u; iVar <= nDim; iVar++) {
auto total_index = iPoint*nVar+iVar;
Jacobian.DeleteValsRowi(total_index);
}
}
}
END_SU2_OMP_FOR

if (Jacobian_i)
for (auto iVar = 0u; iVar < nVar; iVar++)
delete [] Jacobian_i[iVar];
delete [] Jacobian_i;

}


void CSpeciesSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver**, CNumerics*,
  CNumerics*, CConfig *config, unsigned short val_marker) {
/*--- Identify the boundary by string name and get the specified wall
heat flux from config as well as the wall function treatment. ---*/

const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);

/*--- Get the specified wall heat flux, temperature or heat transfer coefficient from config ---*/

su2double Wall_HeatFlux = 0.0, Tinfinity = 0.0, Transfer_Coefficient = 0.0;

Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag) / config->GetHeat_Flux_Ref();


//  Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
//  if (Wall_Function != WALL_FUNCTION::NONE) {
//    SU2_MPI::Error("Wall function treatment not implemented yet", CURRENT_FUNCTION);
//  }

/*--- Jacobian, initialized to zero if needed. ---*/
su2double **Jacobian_i = nullptr;
if (implicit) {
Jacobian_i = new su2double* [nVar];
for (auto iVar = 0u; iVar < nVar; iVar++)
Jacobian_i[iVar] = new su2double [nVar] ();
}

/*--- Loop over all of the vertices on this boundary marker ---*/

SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

if (!geometry->nodes->GetDomain(iPoint)) continue;


/*--- Compute dual-grid area and boundary normal ---*/

const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

su2double Area = GeometryToolbox::Norm(nDim, Normal);

su2double UnitNormal[MAXNDIM] = {0.0};
for (auto iDim = 0u; iDim < nDim; iDim++)
  UnitNormal[iDim] = -Normal[iDim]/Area;

/*--- Apply a weak boundary condition for the energy equation.
Compute the residual due to the prescribed heat flux.
The convective part will be zero if the grid is not moving. ---*/

su2double Res_Conv = 0.0;
su2double Res_Visc = Wall_HeatFlux * Area;

/*--- Impose the value of the velocity as a strong boundary
condition (Dirichlet). Fix the velocity and remove any
contribution to the residual at this node. ---*/

if (dynamic_grid) {
nodes->SetVelocity_Old(iPoint, geometry->nodes->GetGridVel(iPoint));
}
else {
su2double zero[MAXNDIM] = {0.0};
nodes->SetVelocity_Old(iPoint, zero);
}

for (auto iDim = 0u; iDim < nDim; iDim++)
LinSysRes(iPoint, iDim+1) = 0.0;
nodes->SetVel_ResTruncError_Zero(iPoint);

/*--- If the wall is moving, there are additional residual contributions
due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/

/*--- Convective and viscous contributions to the residual at the wall ---*/

LinSysRes(iPoint, 1) += Res_Conv - Res_Visc;

/*--- Enforce the no-slip boundary condition in a strong way by
modifying the velocity-rows of the Jacobian (1 on the diagonal).
And add the contributions to the Jacobian due to energy. ---*/

if (implicit) {

//Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

for (auto iVar = 1u; iVar <= nDim; iVar++) {
auto total_index = iPoint*nVar+iVar;
Jacobian.DeleteValsRowi(total_index);
}
}
}
END_SU2_OMP_FOR

if (Jacobian_i)
for (auto iVar = 0u; iVar < nVar; iVar++)
delete [] Jacobian_i[iVar];
delete [] Jacobian_i;

}