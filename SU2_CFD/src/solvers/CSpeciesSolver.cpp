/*!
 * \file CSpeciesSolver.cpp
 * \brief Main subrotuines of CSpeciesSolver class
 * \author T. Kattmann
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

CSpeciesSolver::CSpeciesSolver() : CScalarSolver<CSpeciesVariable>(true) {}

CSpeciesSolver::CSpeciesSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh)
    : CScalarSolver<CSpeciesVariable>(geometry, config, true) {
  unsigned short nLineLets;

  bool multizone = config->GetMultizone_Problem();

  /*--- Dimension of the problem --> passive scalar will only ever
   have a single equation. Other child classes of CScalarLegacySolver
   can have variable numbers of equations. ---*/

  nVar = config->GetNSpeciesInit();
  nPrimVar = nVar;

  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Single grid simulation ---*/
  /*--- do not see a reason to use only single grid ---*/
  // if (iMesh == MESH_0) {

  /*--- Define some auxiliary vector related with the residual ---*/

  Residual_RMS.resize(nVar, 0.0);
  Residual_Max.resize(nVar, 0.0);
  Point_Max.resize(nVar, 0);
  Point_Max_Coord.resize(nVar, nDim) = su2double(0.0);

  /*--- Initialization of the structure of the whole Jacobian ---*/

  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (passive scalar model)." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);

  if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
    if (rank == MASTER_NODE)
      cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
  }

  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  System.SetxIsZero(true);

  if (ReducerStrategy) EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);

  /*--- Initialize the BGS residuals in multizone problems. ---*/
  if (multizone) {
    Residual_BGS.resize(nVar, 0.0);
    Residual_Max_BGS.resize(nVar, 0.0);
    Point_Max_BGS.resize(nVar, 0);
    Point_Max_Coord_BGS.resize(nVar, nDim) = su2double(0.0);
  }

  //} //iMESH_0

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
  /*--- Far-field flow state quantities and initialization. ---*/
  // su2double Density_Inf, Viscosity_Inf;
  // Density_Inf   = config->GetDensity_FreeStreamND();
  // Viscosity_Inf = config->GetViscosity_FreeStreamND();

  /*--- Set up fluid model for the diffusivity ---*/

  su2double Diffusivity_Ref = 1.0;
  su2double DiffusivityND = config->GetDiffusivity_Constant() / Diffusivity_Ref;
  config->SetDiffusivity_ConstantND(DiffusivityND);

  /*--- Scalar variable state at the far-field. ---*/

  for (auto iVar = 0u; iVar < nVar; iVar++) {
    Species_Inf[iVar] = config->GetSpecies_Init(iVar);
  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CSpeciesVariable(Species_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /// NOTE TK:: This should use a MassDiff model!
  /*--- initialize the mass diffusivity ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      const auto massDiffusivity = config->GetDiffusivity_Constant();  // TK this should be ND
      nodes->SetDiffusivity(iPoint, massDiffusivity, iVar);
    }

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /// NOTE TK:: Hinder sliding mesh in cfg postprocessing

  /*-- Allocation of inlets has to happen in derived classes
   (not CScalarLegacySolver), due to arbitrary number of scalar variables.
   First, we also set the column index for any inlet profiles. ---*/

  const bool turbulent = (config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS);
  const bool turb_SST = turbulent && (config->GetKind_Turb_Model() == TURB_MODEL::SST);
  const bool turb_SA = turbulent && (config->GetKind_Turb_Model() == TURB_MODEL::SA);

  /// NOTE TK:: Make inlet files possible. Note that we have to count for turb as well! See also CDriver.cpp for the
  /// loading

  Inlet_Position = nDim * 2 + 2;
  if (turbulent) {
    if (turb_SA)
      Inlet_Position += 1;
    else if (turb_SST)
      Inlet_Position += 2;
  }

  /*-- Allocation of inlets has to happen in derived classes
    (not CScalarLegacySolver), due to arbitrary number of scalar variables.
    First, we also set the column index for any inlet profiles. ---*/
  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
   * due to arbitrary number of turbulence variables ---*/

  Inlet_SpeciesVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_SpeciesVars[iMarker].resize(nVertex[iMarker], nVar);
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        Inlet_SpeciesVars[iMarker](iVertex, iVar) = Species_Inf[iVar];
      }
    }
  }

  /*--- The turbulence models are always solved implicitly, so set the
  implicit flag in case we have periodic BCs. ---*/

  SetImplicitPeriodic(true);

  /*--- Store the initial CFL number for all grid points. ---*/

  const su2double CFL = config->GetCFL(MGLevel) * config->GetCFLRedCoeff_Species();
  for (auto iPoint = 0u; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "SPECIES";
}

void CSpeciesSolver::LoadRestart(CGeometry** geometry, CSolver*** solver, CConfig* config, int val_iter,
                                 bool val_update_geo) {
  /*--- Restart the solution from file information ---*/

  unsigned short iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine;
  su2double Area_Children, Area_Parent;
  const su2double* Solution_Fine = nullptr;

  string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  /*--- To make this routine safe to call in parallel most of it can only be executed by one thread. ---*/
  SU2_OMP_MASTER {
    /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

    if (config->GetRead_Binary_Restart()) {
      Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
    } else {
      Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
    }

    /*--- Skip flow variables ---*/

    unsigned short skipVars = nDim + solver[MESH_0][FLOW_SOL]->GetnVar();

    /*--- Skip turbulence variables ---*/
    const bool turbulent =
        ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == INC_RANS) ||
         (config->GetKind_Solver() == DISC_ADJ_INC_RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS));

    if (turbulent) skipVars += solver[MESH_0][TURB_SOL]->GetnVar();

    /*--- Adjust the number of solution variables in the incompressible
     restart. We always carry a space in nVar for the energy equation in the
     mean flow solver, but we only write it to the restart if it is active.
     Therefore, we must reduce skipVars here if energy is inactive so that
     the turbulent variables are read correctly. ---*/

    bool incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
    bool energy = config->GetEnergy_Equation();
    bool weakly_coupled_heat = config->GetWeakly_Coupled_Heat();

    if (incompressible && ((!energy) && (!weakly_coupled_heat))) skipVars--;

    /*--- Load data from the restart into correct containers. ---*/

    unsigned long counter = 0, iPoint_Global = 0;
    for (; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {
      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {
        /*--- We need to store this point's data, so jump to the correct
         offset in the buffer of data from the restart file and load it. ---*/

        index = counter * Restart_Vars[1] + skipVars;
        for (iVar = 0; iVar < nVar; ++iVar) nodes->SetSolution(iPoint_Local, iVar, Restart_Data[index + iVar]);

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

  }  // end SU2_OMP_MASTER, pre and postprocessing are thread-safe.
  END_SU2_OMP_MASTER
  SU2_OMP_BARRIER

  /*--- MPI solution and compute the eddy viscosity ---*/

  solver[MESH_0][SPECIES_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][SPECIES_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);

  // Flow-Pre computes/sets mixture properties
  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER,
                                          RUNTIME_FLOW_SYS, false);
  // Update eddy-visc which needs correct mixture density and mixture lam-visc. Note that after this, another Flow-Pre
  // at the start of the Iteration sets the updated eddy-visc into the Flow-Solvers Primitives.
  if (config->GetKind_Turb_Model() != TURB_MODEL::NONE)
    solver[MESH_0][TURB_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
  // For feature_multicomp this Scalar-Pre only computes the laminar contribution to mass diffusivity
  solver[MESH_0][SPECIES_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER,
                                             RUNTIME_FLOW_SYS, false);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    SU2_MPI::Error("Untested code: Remove error and proceed at own risk!", CURRENT_FUNCTION);
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->nodes->GetVolume(iPoint);
      su2double Solution_Coarse[MAXNVAR] = {0.0};
      for (iChildren = 0; iChildren < geometry[iMesh]->nodes->GetnChildren_CV(iPoint); iChildren++) {
        Point_Fine = geometry[iMesh]->nodes->GetChildren_CV(iPoint, iChildren);
        Area_Children = geometry[iMesh - 1]->nodes->GetVolume(Point_Fine);
        Solution_Fine = solver[iMesh - 1][SPECIES_SOL]->GetNodes()->GetSolution(Point_Fine);
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution_Coarse[iVar] += Solution_Fine[iVar] * Area_Children / Area_Parent;
        }
      }
      solver[iMesh][SPECIES_SOL]->GetNodes()->SetSolution(iPoint, Solution_Coarse);
    }
    END_SU2_OMP_FOR

    solver[iMesh][SPECIES_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][SPECIES_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);

    /// NOTE TK:: This probably need to be adapted like above
    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS,
                                           false);
    solver[iMesh][SPECIES_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
  }

  /*--- Go back to single threaded execution. ---*/
  SU2_OMP_MASTER {
    /*--- Delete the class memory that is used to load the restart. ---*/

    delete[] Restart_Vars;
    Restart_Vars = nullptr;
    delete[] Restart_Data;
    Restart_Data = nullptr;
  }
  END_SU2_OMP_MASTER
  SU2_OMP_BARRIER
}

void CSpeciesSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                   unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem,
                                   bool Output) {
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool muscl = config->GetMUSCL_Species();
  const bool limiter =
      (config->GetKind_SlopeLimit_Species() != NO_LIMITER) && (config->GetInnerIter() <= config->GetLimiterIter());

  for (auto iPoint = 0u; iPoint < nPoint; iPoint++) {
    /// NOTE TK:: needs to be value by a model
    const su2double mass_diffusivity = config->GetDiffusivity_Constant();

    for (auto iVar = 0u; iVar < nVar; iVar++) {
      nodes->SetDiffusivity(iPoint, mass_diffusivity, iVar);
    }

  }  // iPoint

  /*--- Clear residual and system matrix, not needed for
   * reducer strategy as we write over the entire matrix. ---*/
  if (!ReducerStrategy && !Output) {
    LinSysRes.SetValZero();
    if (implicit)
      Jacobian.SetValZero();
    else {
      SU2_OMP_BARRIER
    }
  }

  /*--- Upwind second order reconstruction and gradients ---*/

  if (config->GetReconstructionGradientRequired()) {
    switch(config->GetKind_Gradient_Method_Recon()) {
      case GREEN_GAUSS: SetSolution_Gradient_GG(geometry, config, true); break;
      case LEAST_SQUARES: SetSolution_Gradient_LS(geometry, config, true); break;
      case WEIGHTED_LEAST_SQUARES: SetSolution_Gradient_LS(geometry, config, true); break;
    }
  }

  switch(config->GetKind_Gradient_Method()) {
    case GREEN_GAUSS: SetSolution_Gradient_GG(geometry, config); break;
    case WEIGHTED_LEAST_SQUARES: SetSolution_Gradient_LS(geometry, config); break;
  }

  if (limiter && muscl) SetSolution_Limiter(geometry, config);
}

void CSpeciesSolver::SetInitialCondition(CGeometry** geometry, CSolver*** solver_container, CConfig* config,
                                         unsigned long ExtIter) {
  const bool Restart = (config->GetRestart() || config->GetRestart_Flow());

  /// NOTE TK:: Check whether this is done in the constructor for example and whether this is mnecessary. The TurbSolver
  /// does not do this!
  /*--- For a restart do nothing here. Otherwise initialize the solution with the Init value. ---*/
  if ((!Restart) && ExtIter == 0) {
    if (rank == MASTER_NODE) {
      cout << "Initializing passive scalar (initial condition)." << endl;
      cout << "initialization = " << nVar << " " << config->GetSpecies_Init(0) << endl;
    }

    // TK:: arbitrary number which would be easy to recognize during debugging
    su2double scalar_init[MAXNVAR] = {1234.0};

    for (unsigned long iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (unsigned long iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {

        for (int iVar = 0; iVar < nVar; iVar++) {
          scalar_init[iVar] = config->GetSpecies_Init(iVar);
        }

        solver_container[iMesh][SPECIES_SOL]->GetNodes()->SetSolution(iPoint, scalar_init);
      }

      solver_container[iMesh][SPECIES_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
      solver_container[iMesh][SPECIES_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);

      solver_container[iMesh][FLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
      solver_container[iMesh][FLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);

      solver_container[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver_container[iMesh], config, iMesh,
                                                       NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    }
  }
}

void CSpeciesSolver::Viscous_Residual(unsigned long iEdge, CGeometry* geometry, CSolver** solver_container,
                                      CNumerics* numerics, CConfig* config) {
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
  const bool use_strong_BC = !config->GetMUSCL_AdjFlow();  // hook to cfg

  /// NOTE TK:: This is a strong impl whereas TurbSA and inceuler implement a weak version. Testing required.
  // bool grid_movement  = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  // su2double   temp_inlet    = config->GetInlet_Ttotal       (Marker_Tag);
  const su2double* inlet_species = config->GetInlet_SpeciesVal(Marker_Tag);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    if (use_strong_BC) {
      nodes->SetSolution_Old(iPoint, inlet_species);

      LinSysRes.SetBlock_Zero(iPoint);

      for (auto iVar = 0u; iVar < nVar; iVar++) {
        nodes->SetVal_ResTruncError_Zero(iPoint, iVar);
      }

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

      /*--- Set the turbulent variable states. Use free-stream SST
      values for the turbulent state at the inflow. ---*/
      /*--- Load the inlet turbulence variables (uniform by default). ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), inlet_species);

      /*--- Set various other quantities in the solver class ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/
      const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      //      /*--- Viscous contribution, commented out because serious convergence problems ---*/
      //
      //      su2double Coord_Reflected[MAXNDIM];
      //      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
      //                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      //      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      //      visc_numerics->SetNormal(Normal);
      //
      //      /*--- Conservative variables w/o reconstruction ---*/
      //
      //      visc_numerics->SetPrimitive(V_domain, V_inlet);
      //
      //      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      //
      //     visc_numerics->SetScalarVar(Solution_i, Solution_j);
      //     visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      //
      //      /*--- Menter's first blending function ---*/
      //
      //      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());
      //
      //      /*--- Compute residual, and Jacobians ---*/
      //
      //      auto residual = visc_numerics->ComputeResidual(config);
      //
      //      /*--- Subtract residual, and update Jacobians ---*/
      //
      //      LinSysRes.SubtractBlock(iPoint, residual);
      //      Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
    }
  }
  END_SU2_OMP_FOR
}

void CSpeciesSolver::BC_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                               CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  const bool use_strong_BC = !config->GetMUSCL_AdjFlow();  // hook to cfg
  /// NOTE TK:: This is a strong impl whereas TurbSA and inceuler implement a weak version. Testing required.
  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    /* strong zero flux Neumann boundary condition at the outlet */
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    if (use_strong_BC) {
      /*--- Allocate the value at the outlet ---*/
      auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      nodes->SetSolution_Old(iPoint, nodes->GetSolution(Point_Normal));

      LinSysRes.SetBlock_Zero(iPoint);

      for (auto iVar = 0u; iVar < nVar; iVar++) {
        nodes->SetVal_ResTruncError_Zero(iPoint, iVar);
      }

      /*--- Includes 1 in the diagonal ---*/
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

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/
      const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      //      /*--- Viscous contribution, commented out because serious convergence problems ---*/
      //
      //      su2double Coord_Reflected[MAXNDIM];
      //      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
      //                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      //      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      //      visc_numerics->SetNormal(Normal);
      //
      //      /*--- Conservative variables w/o reconstruction ---*/
      //
      //      visc_numerics->SetPrimitive(V_domain, V_outlet);
      //
      //      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      //
      //      visc_numerics->SetScalarVar(Solution_i, Solution_j);
      //      visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      //
      //      /*--- Menter's first blending function ---*/
      //
      //      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());
      //
      //      /*--- Compute residual, and Jacobians ---*/
      //
      //      auto residual = visc_numerics->ComputeResidual(config);
      //
      //      /*--- Subtract residual, and update Jacobians ---*/
      //
      //      LinSysRes.SubtractBlock(iPoint, residual);
      //      Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
    }
  }
  END_SU2_OMP_FOR
}

void CSpeciesSolver::BC_HeatFlux_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                      CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {}

void CSpeciesSolver::BC_Isothermal_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                        CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CSpeciesSolver::SetUniformInlet(const CConfig* config, unsigned short iMarker) {
  for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Inlet_SpeciesVars[iMarker][iVertex][iVar] = Species_Inf[iVar];
  }
}