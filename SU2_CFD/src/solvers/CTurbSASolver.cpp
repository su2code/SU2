/*!
 * \file CTurbSASolver.cpp
 * \brief Main subrotuines of CTurbSASolver class
 * \author F. Palacios, A. Bueno
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

#include "../../include/solvers/CTurbSASolver.hpp"
#include "../../include/variables/CTurbSAVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"


CTurbSASolver::CTurbSASolver(void) : CTurbSolver(false) { }

CTurbSASolver::CTurbSASolver(CGeometry *geometry, CConfig *config, unsigned short iMesh, CFluidModel* FluidModel)
             : CTurbSolver(geometry, config, false) {

  unsigned short nLineLets;
  unsigned long iPoint;
  su2double Density_Inf, Viscosity_Inf, Factor_nu_Inf, Factor_nu_Engine, Factor_nu_ActDisk;

  bool multizone = config->GetMultizone_Problem();

  /*--- Dimension of the problem --> dependent of the turbulent model ---*/

  nVar = 1;
  nPrimVar = 1;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Single grid simulation ---*/

  if (iMesh == MESH_0 || config->GetMGCycle() == FULLMG_CYCLE) {

    /*--- Define some auxiliar vector related with the residual ---*/

    Residual_RMS.resize(nVar,0.0);
    Residual_Max.resize(nVar,0.0);
    Point_Max.resize(nVar,0);
    Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (SA model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);

    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    System.SetxIsZero(true);

    if (ReducerStrategy)
      EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);

    if (config->GetExtraOutput()) {
      if (nDim == 2) { nOutputVariables = 13; }
      else if (nDim == 3) { nOutputVariables = 19; }
      OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
      OutputHeadingNames = new string[nOutputVariables];
    }

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (multizone){
      Residual_BGS.resize(nVar,0.0);
      Residual_Max_BGS.resize(nVar,0.0);
      Point_Max_BGS.resize(nVar,0);
      Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
    }

  }

  /*--- Read farfield conditions from config ---*/

  Density_Inf   = config->GetDensity_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();

  /*--- Factor_nu_Inf in [3.0, 5.0] ---*/

  Factor_nu_Inf = config->GetNuFactor_FreeStream();
  su2double nu_tilde_Inf  = Factor_nu_Inf*Viscosity_Inf/Density_Inf;
  if (config->GetKind_Trans_Model() == BC) {
    nu_tilde_Inf  = 0.005*Factor_nu_Inf*Viscosity_Inf/Density_Inf;
  }

  Solution_Inf[0] = nu_tilde_Inf;

  /*--- Factor_nu_Engine ---*/
  Factor_nu_Engine = config->GetNuFactor_Engine();
  nu_tilde_Engine  = Factor_nu_Engine*Viscosity_Inf/Density_Inf;
  if (config->GetKind_Trans_Model() == BC) {
    nu_tilde_Engine  = 0.005*Factor_nu_Engine*Viscosity_Inf/Density_Inf;
  }

  /*--- Factor_nu_ActDisk ---*/
  Factor_nu_ActDisk = config->GetNuFactor_Engine();
  nu_tilde_ActDisk  = Factor_nu_ActDisk*Viscosity_Inf/Density_Inf;

  /*--- Eddy viscosity at infinity ---*/
  su2double Ji, Ji_3, fv1, cv1_3 = 7.1*7.1*7.1;
  su2double muT_Inf;
  Ji = nu_tilde_Inf/Viscosity_Inf*Density_Inf;
  Ji_3 = Ji*Ji*Ji;
  fv1 = Ji_3/(Ji_3+cv1_3);
  muT_Inf = Density_Inf*fv1*nu_tilde_Inf;

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CTurbSAVariable(nu_tilde_Inf, muT_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION_EDDY);
  CompleteComms(geometry, config, SOLUTION_EDDY);

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  SlidingState.resize(nMarker);
  SlidingStateNodes.resize(nMarker);

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
      SlidingState[iMarker].resize(nVertex[iMarker], nPrimVar+1) = nullptr;
      SlidingStateNodes[iMarker].resize(nVertex[iMarker],0);
    }
  }

  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
   * due to arbitrary number of turbulence variables ---*/

  Inlet_TurbVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++)
    Inlet_TurbVars[iMarker].resize(nVertex[iMarker],nVar) = nu_tilde_Inf;

  /*--- The turbulence models are always solved implicitly, so set the
   implicit flag in case we have periodic BCs. ---*/

  SetImplicitPeriodic(true);

  /*--- Store the initial CFL number for all grid points. ---*/

  const su2double CFL = config->GetCFL(MGLevel)*config->GetCFLRedCoeff_Turb();
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "SA";

}

void CTurbSASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
        unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool muscl = config->GetMUSCL_Turb();
  const bool limiter = (config->GetKind_SlopeLimit_Turb() != NO_LIMITER) &&
                       (config->GetInnerIter() <= config->GetLimiterIter());
  const auto kind_hybridRANSLES = config->GetKind_HybridRANSLES();

  /*--- Clear residual and system matrix, not needed for
   * reducer strategy as we write over the entire matrix. ---*/
  if (!ReducerStrategy) {
    LinSysRes.SetValZero();
    if (implicit) Jacobian.SetValZero();
    else {SU2_OMP_BARRIER}
  }

  /*--- Upwind second order reconstruction and gradients ---*/

  if (config->GetReconstructionGradientRequired()) {
    if (config->GetKind_Gradient_Method_Recon() == GREEN_GAUSS)
      SetSolution_Gradient_GG(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == LEAST_SQUARES)
      SetSolution_Gradient_LS(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES)
      SetSolution_Gradient_LS(geometry, config, true);
  }

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
    SetSolution_Gradient_GG(geometry, config);

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
    SetSolution_Gradient_LS(geometry, config);

  if (limiter && muscl) SetSolution_Limiter(geometry, config);

  if (kind_hybridRANSLES != NO_HYBRIDRANSLES) {

    /*--- Set the vortex tilting coefficient at every node if required ---*/

    if (kind_hybridRANSLES == SA_EDDES){
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        auto Vorticity = solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint);
        auto PrimGrad_Flow = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);
        auto Laminar_Viscosity  = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
        nodes->SetVortex_Tilting(iPoint, PrimGrad_Flow, Vorticity, Laminar_Viscosity);
      }
      END_SU2_OMP_FOR
    }

    /*--- Compute the DES length scale ---*/

    SetDES_LengthScale(solver_container, geometry, config);

  }

}

void CTurbSASolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

  const su2double cv1_3 = 7.1*7.1*7.1, cR1 = 0.5, rough_const = 0.03;

  const bool neg_spalart_allmaras = (config->GetKind_Turb_Model() == SA_NEG);

  /*--- Compute eddy viscosity ---*/

  AD::StartNoSharedReading();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    su2double rho = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    su2double mu  = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);

    su2double nu  = mu/rho;
    su2double nu_hat = nodes->GetSolution(iPoint,0);
    su2double roughness = geometry->nodes->GetRoughnessHeight(iPoint);
    su2double dist = geometry->nodes->GetWall_Distance(iPoint);

    dist += rough_const*roughness;

    su2double Ji   = nu_hat/nu ;
    if (roughness > 1.0e-10)
      Ji+= cR1*roughness/(dist+EPS);

    su2double Ji_3 = Ji*Ji*Ji;
    su2double fv1  = Ji_3/(Ji_3+cv1_3);

    su2double muT = rho*fv1*nu_hat;

    if (neg_spalart_allmaras) muT = max(muT,0.0);

    nodes->SetmuT(iPoint,muT);

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();
}

void CTurbSASolver::Viscous_Residual(unsigned long iEdge, CGeometry* geometry, CSolver** solver_container,
                                     CNumerics* numerics, CConfig* config) {

  /*--- Define an object to set solver specific numerics contribution. ---*/
  struct SolverSpecificNumerics {

    FORCEINLINE void operator() (const CGeometry& geometry, CNumerics& numerics, CConfig& config, const CScalarVariable& nodes, unsigned long iPoint, unsigned long jPoint) const {
      /*--- Roughness heights. ---*/
      if (config.GetKind_Turb_Model() == SA)
        numerics.SetRoughness(geometry.nodes->GetRoughnessHeight(iPoint), geometry.nodes->GetRoughnessHeight(jPoint));
    }

  } SolverSpecificNumerics;

  /*--- Now instantiate the generic implementation with the functor above. ---*/

  Viscous_Residual_impl(SolverSpecificNumerics, iEdge, geometry, solver_container, numerics, config);
}

void CTurbSASolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                    CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool harmonic_balance = (config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE);
  const bool transition    = (config->GetKind_Trans_Model() == LM);
  const bool transition_BC = (config->GetKind_Trans_Model() == BC);

  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  AD::StartNoSharedReading();

  /*--- Loop over all points. ---*/

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), nullptr);

    /*--- Gradient of the primitive and conservative variables ---*/

    numerics->SetPrimVarGradient(flowNodes->GetGradient_Primitive(iPoint), nullptr);

    /*--- Set vorticity and strain rate magnitude ---*/

    numerics->SetVorticity(flowNodes->GetVorticity(iPoint), nullptr);

    numerics->SetStrainMag(flowNodes->GetStrainMag(iPoint), 0.0);

    /*--- Set intermittency ---*/

    if (transition) {
      numerics->SetIntermittency(solver_container[TRANS_SOL]->GetNodes()->GetIntermittency(iPoint));
    }

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/

    numerics->SetScalarVar(nodes->GetSolution(iPoint), nullptr);
    numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nullptr);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Get Hybrid RANS/LES Type and set the appropriate wall distance ---*/

    if (config->GetKind_HybridRANSLES() == NO_HYBRIDRANSLES) {

    /*--- For the SA model, wall roughness is accounted by modifying the computed wall distance
       *                              d_new = d + 0.03 k_s
       *    where k_s is the equivalent sand grain roughness height that is specified in cfg file.
       *    For smooth walls, wall roughness is zero and computed wall distance remains the same. */

      su2double modifiedWallDistance = geometry->nodes->GetWall_Distance(iPoint);

      modifiedWallDistance += 0.03*geometry->nodes->GetRoughnessHeight(iPoint);

      /*--- Set distance to the surface ---*/

      numerics->SetDistance(modifiedWallDistance, 0.0);

      /*--- Set the roughness of the closest wall. ---*/

      numerics->SetRoughness(geometry->nodes->GetRoughnessHeight(iPoint), 0.0 );

    } else {

      /*--- Set DES length scale ---*/

      numerics->SetDistance(nodes->GetDES_LengthScale(iPoint), 0.0);

    }

    /*--- Compute the source term ---*/

    auto residual = numerics->ComputeResidual(config);

    /*--- Store the intermittency ---*/

    if (transition_BC) {
      nodes->SetGammaBC(iPoint,numerics->GetGammaBC());
    }

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, residual);

    if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

  }
  END_SU2_OMP_FOR

  if (harmonic_balance) {

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

      su2double Volume = geometry->nodes->GetVolume(iPoint);

      /*--- Access stored harmonic balance source term ---*/

      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        su2double Source = nodes->GetHarmonicBalance_Source(iPoint,iVar);
        LinSysRes(iPoint,iVar) += Source*Volume;
      }
    }
    END_SU2_OMP_FOR
  }

  AD::EndNoSharedReading();

}

void CTurbSASolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                    CConfig *config, unsigned short iMesh) {
}

void CTurbSASolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Evaluate nu tilde at the closest point to the surface using the wall functions. ---*/

  if (config->GetWall_Functions()) {
    SU2_OMP_MASTER
    SetTurbVars_WF(geometry, solver_container, config, val_marker);
    END_SU2_OMP_MASTER
    SU2_OMP_BARRIER
    return;
  }

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool rough_wall = false;
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  WALL_TYPE WallType;
  su2double Roughness_Height;
  tie(WallType, Roughness_Height) = config->GetWallRoughnessProperties(Marker_Tag);
  if (WallType == WALL_TYPE::ROUGH) rough_wall = true;

  /*--- The dirichlet condition is used only without wall function, otherwise the
   convergence is compromised as we are providing nu tilde values for the
   first point of the wall  ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {
      if (!rough_wall) {
        for (auto iVar = 0u; iVar < nVar; iVar++)
          nodes->SetSolution_Old(iPoint,iVar,0.0);

        LinSysRes.SetBlock_Zero(iPoint);

        /*--- Includes 1 in the diagonal ---*/

        if (implicit) Jacobian.DeleteValsRowi(iPoint);
       } else {
         /*--- For rough walls, the boundary condition is given by
          * (\frac{\partial \nu}{\partial n})_wall = \frac{\nu}{0.03*k_s}
          * where \nu is the solution variable, $n$ is the wall normal direction
          * and k_s is the equivalent sand grain roughness specified. ---*/

         /*--- Compute dual-grid area and boundary normal ---*/
         su2double Normal[MAXNDIM] = {0.0};
         for (auto iDim = 0u; iDim < nDim; iDim++)
           Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);

         su2double Area = GeometryToolbox::Norm(nDim, Normal);

         /*--- Get laminar_viscosity and density ---*/
         su2double sigma = 2.0/3.0;
         su2double laminar_viscosity = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
         su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);

         su2double nu_total = (laminar_viscosity/density + nodes->GetSolution(iPoint,0));

         su2double coeff = (nu_total/sigma);
         su2double RoughWallBC = nodes->GetSolution(iPoint,0)/(0.03*Roughness_Height);

         su2double Res_Wall;// = new su2double [nVar];
         Res_Wall = coeff*RoughWallBC*Area;
         LinSysRes.SubtractBlock(iPoint, &Res_Wall);

         su2double Jacobian_i = (laminar_viscosity*Area)/(0.03*Roughness_Height*sigma);
         Jacobian_i += 2.0*RoughWallBC*Area/sigma;
         if (implicit) Jacobian.AddVal2Diag(iPoint, -Jacobian_i);
      }
    }
  }
  END_SU2_OMP_FOR
}

void CTurbSASolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                       CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}

void CTurbSASolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                 CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Allocate the value at the infinity ---*/

      auto V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Grid Movement ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Set turbulent variable at the wall, and at infinity ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Solution_Inf);

      /*--- Set Normal (it is necessary to change the sign) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Compute residuals and Jacobians ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Add residuals and Jacobians ---*/

      LinSysRes.AddBlock(iPoint, residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR

}

void CTurbSASolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                             CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);

      /*--- Allocate the value at the inlet ---*/

      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Load the inlet turbulence variable (uniform by default). ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Inlet_TurbVars[val_marker][iVertex]);

      /*--- Set various other quantities in the conv_numerics class ---*/

      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

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
//      visc_numerics->SetScalarVar(Solution_i, Solution_j);
//      visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
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

void CTurbSASolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                              CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Allocate the value at the outlet ---*/

      auto V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_outlet);

      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual. ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), nodes->GetSolution(iPoint));

      /*--- Set Normal (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

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

void CTurbSASolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Allocate the value at the infinity ---*/

      auto V_inflow = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inflow);

      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual. ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), nodes->GetSolution(iPoint));

      /*--- Set Normal (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Set grid movement ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

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
//      visc_numerics->SetPrimitive(V_domain, V_inflow);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetScalarVar(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
//      visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
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

void CTurbSASolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);

      /*--- Allocate the value at the infinity ---*/

      auto V_exhaust = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_exhaust);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), &nu_tilde_Engine);

      /*--- Set various other quantities in the conv_numerics class ---*/

      conv_numerics->SetNormal(Normal);

      /*--- Set grid movement ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

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
//      visc_numerics->SetPrimitive(V_domain, V_exhaust);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetScalarVar(Solution_i, Solution_j);
//      visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
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

void CTurbSASolver::BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics, config,  val_marker, true);
}

void CTurbSASolver::BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics, config,  val_marker, false);
}

void CTurbSASolver::BC_ActDisk(CGeometry *geometry, CSolver **solver_container,
                               CNumerics *conv_numerics, CNumerics *visc_numerics,
                               CConfig *config, unsigned short val_marker, bool val_inlet_surface) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    const auto GlobalIndex_donor = solver_container[FLOW_SOL]->GetDonorGlobalIndex(val_marker, iVertex);
    const auto GlobalIndex = geometry->nodes->GetGlobalIndex(iPoint);

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (!geometry->nodes->GetDomain(iPoint) || (GlobalIndex == GlobalIndex_donor)) {
      continue;
    }

    /*--- Normal vector for this vertex (negate for outward convention) ---*/

    su2double Normal[MAXNDIM] = {0.0};
    for (auto iDim = 0u; iDim < nDim; iDim++)
      Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
    conv_numerics->SetNormal(Normal);

    su2double Area = GeometryToolbox::Norm(nDim, Normal);

    su2double UnitNormal[MAXNDIM] = {0.0};
    for (auto iDim = 0u; iDim < nDim; iDim++)
      UnitNormal[iDim] = Normal[iDim]/Area;

    /*--- Retrieve solution at the farfield boundary node ---*/

    auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

    /*--- Check the flow direction. Project the flow into the normal to the inlet face ---*/

    su2double Vn = GeometryToolbox::DotProduct(nDim, &V_domain[1], UnitNormal);

    bool ReverseFlow = false;
    if ((val_inlet_surface) && (Vn < 0.0)) { ReverseFlow = true; }
    if ((!val_inlet_surface) && (Vn > 0.0)) { ReverseFlow = true; }

    /*--- Do not anything if there is a
     reverse flow, Euler b.c. for the direct problem ---*/

    if (ReverseFlow) continue;

    /*--- Allocate the value at the infinity ---*/

    if (val_inlet_surface) {
      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      conv_numerics->SetPrimitive(V_domain, V_inlet);
    }
    else {
      auto V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      conv_numerics->SetPrimitive(V_domain, V_outlet);
    }

    /*--- Set the turb. variable solution
     set  the turbulent variables. Here we use a Neumann BC such
     that the turbulent variable is copied from the interior of the
     domain to the outlet before computing the residual.
     or set the turbulent variable states (prescribed for an inflow)  ----*/

    //      if (val_inlet_surface) Solution_j[0] = 0.5*(nodes->GetSolution(iPoint,0)+V_outlet [nDim+9]);
    //      else Solution_j[0] = 0.5*(nodes->GetSolution(iPoint,0)+V_inlet [nDim+9]);

    //      /*--- Inflow analysis (interior extrapolation) ---*/
    //      if (((val_inlet_surface) && (!ReverseFlow)) || ((!val_inlet_surface) && (ReverseFlow))) {
    //        Solution_j[0] = 2.0*node[iPoint]->GetSolution(0) - node[iPoint_Normal]->GetSolution(0);
    //      }

    //      /*--- Outflow analysis ---*/
    //      else {
    //        if (val_inlet_surface) Solution_j[0] = Factor_nu_ActDisk*V_outlet [nDim+9];
    //        else { Solution_j[0] = Factor_nu_ActDisk*V_inlet [nDim+9]; }
    //      }

    if (((val_inlet_surface) && (!ReverseFlow)) || ((!val_inlet_surface) && (ReverseFlow))) {
      /*--- Inflow analysis (interior extrapolation) ---*/
      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), nodes->GetSolution(iPoint));
    }
    else {
      /*--- Outflow analysis ---*/
      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), &nu_tilde_ActDisk);
    }

    /*--- Grid Movement ---*/

    if (dynamic_grid)
      conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

    /*--- Compute the residual using an upwind scheme ---*/

    auto residual = conv_numerics->ComputeResidual(config);
    LinSysRes.AddBlock(iPoint, residual);

    /*--- Jacobian contribution for implicit integration ---*/

    if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

//        /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//        visc_numerics->SetNormal(Normal);
//        su2double Coord_Reflected[MAXNDIM];
//        GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
//                                                 geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//
//        /*--- Conservative variables w/o reconstruction ---*/
//
//        if (val_inlet_surface) visc_numerics->SetPrimitive(V_domain, V_inlet);
//        else visc_numerics->SetPrimitive(V_domain, V_outlet);
//
//        /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//        visc_numerics->SetScalarVar(Solution_i, Solution_j);
//
//        visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//        /*--- Compute residual, and Jacobians ---*/
//
//        auto residual = visc_numerics->ComputeResidual(config);
//
//        /*--- Subtract residual, and update Jacobians ---*/
//
//        LinSysRes.SubtractBlock(iPoint, residual);
//        Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

  }
  END_SU2_OMP_FOR

}

void CTurbSASolver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                         CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const auto nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (auto iSpan = 0u; iSpan < nSpanWiseSections ; iSpan++){

    su2double extAverageNu = solver_container[FLOW_SOL]->GetExtAverageNu(val_marker, iSpan);

    /*--- Loop over all the vertices on this boundary marker ---*/

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {

      /*--- find the node related to the vertex ---*/
      const auto iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      const auto oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      const auto Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][oldVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/
      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), &extAverageNu);

      /*--- Set various other quantities in the conv_numerics class ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto conv_residual = conv_numerics->ComputeResidual(config);

      /*--- Jacobian contribution for implicit integration ---*/

      LinSysRes.AddBlock(iPoint, conv_residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, conv_residual.jacobian_i);

      /*--- Viscous contribution ---*/
      su2double Coord_Reflected[MAXNDIM];
      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/

      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/

      visc_numerics->SetScalarVar(nodes->GetSolution(iPoint), &extAverageNu);

      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint),
                                        nodes->GetGradient(iPoint));

      /*--- Compute residual, and Jacobians ---*/

      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/

      LinSysRes.SubtractBlock(iPoint, visc_residual);
      if (implicit) Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

}

void CTurbSASolver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                   CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const auto nSpanWiseSections = config->GetnSpanWiseSections();

  CFluidModel *FluidModel = solver_container[FLOW_SOL]->GetFluidModel();

  su2double Factor_nu_Inf = config->GetNuFactor_FreeStream();

  /*--- Loop over all the spans on this boundary marker ---*/
  for (auto iSpan = 0; iSpan < nSpanWiseSections; iSpan++) {

    su2double rho       = solver_container[FLOW_SOL]->GetAverageDensity(val_marker, iSpan);
    su2double pressure  = solver_container[FLOW_SOL]->GetAveragePressure(val_marker, iSpan);

    FluidModel->SetTDState_Prho(pressure, rho);
    su2double muLam = FluidModel->GetLaminarViscosity();

    su2double nu_tilde  = Factor_nu_Inf*muLam/rho;

    /*--- Loop over all the vertices on this boundary marker ---*/

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {

      /*--- find the node related to the vertex ---*/
      const auto iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      const auto oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      const auto Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][oldVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/
      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), &nu_tilde);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto conv_residual = conv_numerics->ComputeResidual(config);

      /*--- Jacobian contribution for implicit integration ---*/

      LinSysRes.AddBlock(iPoint, conv_residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, conv_residual.jacobian_i);

      /*--- Viscous contribution ---*/

      su2double Coord_Reflected[MAXNDIM];
      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);

      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/

      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/

      visc_numerics->SetScalarVar(nodes->GetSolution(iPoint), &nu_tilde);

      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint),
                                        nodes->GetGradient(iPoint));

      /*--- Compute residual, and Jacobians ---*/

      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/

      LinSysRes.SubtractBlock(iPoint, visc_residual);
      if (implicit) Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

}

void CTurbSASolver::SetTurbVars_WF(CGeometry *geometry, CSolver **solver_container,
                                  const CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- We use a very high max nr of iterations, but we only need this the first couple of iterations ---*/
  constexpr unsigned short max_iter = 200;

  /* --- tolerance has LARGE impact on convergence, do not increase this value! --- */
  const su2double tol = 1e-12;
  su2double relax = 0.5;            /*--- relaxation factor for the Newton solver ---*/

  /*--- Typical constants from boundary layer theory ---*/

  const su2double cv1_3 = 7.1*7.1*7.1;

  CVariable* flow_nodes = solver_container[FLOW_SOL]->GetNodes();

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    const auto iPoint_Neighbor = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint_Neighbor)) {

      su2double Y_Plus = solver_container[FLOW_SOL]->GetYPlus(val_marker, iVertex);

      /*--- note that we do not do anything for y+ < 5, meaning that we have a zero flux (Neumann) boundary condition ---*/

      if (Y_Plus < 5.0) continue;

      su2double Lam_Visc_Normal = flow_nodes->GetLaminarViscosity(iPoint_Neighbor);
      su2double Density_Normal = flow_nodes->GetDensity(iPoint_Neighbor);
      su2double Kin_Visc_Normal = Lam_Visc_Normal/Density_Normal;

      su2double Eddy_Visc = solver_container[FLOW_SOL]->GetEddyViscWall(val_marker, iVertex);

      /*--- Solve for the new value of nu_tilde given the eddy viscosity and using a Newton method ---*/

      // start with positive value of nu_til_old
      su2double nu_til = 0.0;
      su2double nu_til_old = nodes->GetSolution(iPoint,0);

      unsigned short counter = 0;
      su2double diff = 1.0;
      relax = 0.5;
      while (diff > tol) {
        // note the error in Nichols and Nelson
        su2double func = pow(nu_til_old,4) - (Eddy_Visc/Density_Normal)*(pow(nu_til_old,3) + pow(Kin_Visc_Normal,3)*cv1_3);
        su2double func_prim = 4.0 * pow(nu_til_old,3) - 3.0*(Eddy_Visc/Density_Normal)*pow(nu_til_old,2);

        // damped Newton method
        nu_til = nu_til_old - relax*(func/func_prim);

        diff = fabs(nu_til-nu_til_old);
        nu_til_old = nu_til;

        // sometimes we get negative values when the solution has not converged yet, we just reset the nu_tilde in that case.
        if (nu_til_old<tol) {
          relax /= 2.0;
          nu_til_old = nodes->GetSolution(iPoint,0)/relax;
        }

        counter++;
        if (counter > max_iter) break;
      }

      nodes->SetSolution_Old(iPoint_Neighbor, &nu_til);
      LinSysRes.SetBlock_Zero(iPoint_Neighbor);

      /*--- includes 1 in the diagonal ---*/

      if (implicit) Jacobian.DeleteValsRowi(iPoint_Neighbor);
    }
  }
}

void CTurbSASolver::SetDES_LengthScale(CSolver **solver, CGeometry *geometry, CConfig *config){

  unsigned short kindHybridRANSLES = config->GetKind_HybridRANSLES();
  unsigned long iPoint = 0, jPoint = 0;
  unsigned short iDim = 0, jDim = 0, iNeigh = 0, nNeigh = 0;

  su2double constDES = config->GetConst_DES();

  su2double density = 0.0, laminarViscosity = 0.0, kinematicViscosity = 0.0,
      eddyViscosity = 0.0, kinematicViscosityTurb = 0.0, wallDistance = 0.0, lengthScale = 0.0;

  su2double maxDelta = 0.0, deltaAux = 0.0, distDES = 0.0, uijuij = 0.0, k2 = 0.0, r_d = 0.0, f_d = 0.0;
  su2double deltaDDES = 0.0, omega = 0.0, ln_max = 0.0, ln[3] = {0.0}, aux_ln = 0.0, f_kh = 0.0;

  su2double nu_hat, fw_star = 0.424, cv1_3 = pow(7.1, 3.0); k2 = pow(0.41, 2.0);
  su2double cb1   = 0.1355, ct3 = 1.2, ct4   = 0.5;
  su2double sigma = 2./3., cb2 = 0.622, f_max=1.0, f_min=0.1, a1=0.15, a2=0.3;
  su2double cw1 = 0.0, Ji = 0.0, Ji_2 = 0.0, Ji_3 = 0.0, fv1 = 0.0, fv2 = 0.0, ft2 = 0.0, psi_2 = 0.0;
  const su2double *coord_i = nullptr, *coord_j = nullptr, *vorticity = nullptr;
  su2double delta[3] = {0.0}, ratioOmega[3] = {0.0}, vortexTiltingMeasure = 0.0;

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){

    coord_i                 = geometry->nodes->GetCoord(iPoint);
    nNeigh                  = geometry->nodes->GetnPoint(iPoint);
    wallDistance            = geometry->nodes->GetWall_Distance(iPoint);
    const auto primVarGrad  = solver[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);
    vorticity               = solver[FLOW_SOL]->GetNodes()->GetVorticity(iPoint);
    density                 = solver[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    laminarViscosity        = solver[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    eddyViscosity           = solver[TURB_SOL]->GetNodes()->GetmuT(iPoint);
    kinematicViscosity      = laminarViscosity/density;
    kinematicViscosityTurb  = eddyViscosity/density;

    uijuij = 0.0;
    for(iDim = 0; iDim < nDim; iDim++){
      for(jDim = 0; jDim < nDim; jDim++){
        uijuij += primVarGrad[1+iDim][jDim]*primVarGrad[1+iDim][jDim];
      }
    }
    uijuij = sqrt(fabs(uijuij));
    uijuij = max(uijuij,1e-10);

    /*--- Low Reynolds number correction term ---*/

    nu_hat = nodes->GetSolution(iPoint,0);
    Ji   = nu_hat/kinematicViscosity;
    Ji_2 = Ji * Ji;
    Ji_3 = Ji*Ji*Ji;
    fv1  = Ji_3/(Ji_3+cv1_3);
    fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    ft2 = ct3*exp(-ct4*Ji_2);
    cw1 = cb1/k2+(1.0+cb2)/sigma;

    psi_2 = (1.0 - (cb1/(cw1*k2*fw_star))*(ft2 + (1.0 - ft2)*fv2))/(fv1 * max(1.0e-10,1.0-ft2));
    psi_2 = min(100.0,psi_2);

    switch(kindHybridRANSLES){
      case SA_DES:
        /*--- Original Detached Eddy Simulation (DES97)
        Spalart
        1997
        ---*/

        maxDelta = geometry->nodes->GetMaxLength(iPoint);
        distDES = constDES * maxDelta;
        lengthScale = min(distDES,wallDistance);

        break;

      case SA_DDES:
        /*--- A New Version of Detached-eddy Simulation, Resistant to Ambiguous Grid Densities.
         Spalart et al.
         Theoretical and Computational Fluid Dynamics - 2006
         ---*/

        maxDelta = geometry->nodes->GetMaxLength(iPoint);

        r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
        f_d = 1.0-tanh(pow(8.0*r_d,3.0));

        distDES = constDES * maxDelta;
        lengthScale = wallDistance-f_d*max(0.0,(wallDistance-distDES));

        break;
      case SA_ZDES:
        /*--- Recent improvements in the Zonal Detached Eddy Simulation (ZDES) formulation.
         Deck
         Theoretical and Computational Fluid Dynamics - 2012
         ---*/

        for (iNeigh = 0; iNeigh < nNeigh; iNeigh++){
            jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
            coord_j = geometry->nodes->GetCoord(jPoint);
            for ( iDim = 0; iDim < nDim; iDim++){
              deltaAux = abs(coord_j[iDim] - coord_i[iDim]);
              delta[iDim] = max(delta[iDim], deltaAux);
            }
            deltaDDES = geometry->nodes->GetMaxLength(iPoint);
        }

        omega = sqrt(vorticity[0]*vorticity[0] +
                     vorticity[1]*vorticity[1] +
                     vorticity[2]*vorticity[2]);

        for (iDim = 0; iDim < 3; iDim++){
          ratioOmega[iDim] = vorticity[iDim]/omega;
        }

        maxDelta = sqrt(pow(ratioOmega[0],2.0)*delta[1]*delta[2] +
                        pow(ratioOmega[1],2.0)*delta[0]*delta[2] +
                        pow(ratioOmega[2],2.0)*delta[0]*delta[1]);

        r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
        f_d = 1.0-tanh(pow(8.0*r_d,3.0));

        if (f_d < 0.99){
          maxDelta = deltaDDES;
        }

        distDES = constDES * maxDelta;
        lengthScale = wallDistance-f_d*max(0.0,(wallDistance-distDES));

        break;

      case SA_EDDES:

        /*--- An Enhanced Version of DES with Rapid Transition from RANS to LES in Separated Flows.
         Shur et al.
         Flow Turbulence Combust - 2015
         ---*/

        vortexTiltingMeasure = nodes->GetVortex_Tilting(iPoint);

        omega = sqrt(vorticity[0]*vorticity[0] +
                     vorticity[1]*vorticity[1] +
                     vorticity[2]*vorticity[2]);

        for (iDim = 0; iDim < 3; iDim++){
          ratioOmega[iDim] = vorticity[iDim]/omega;
        }

        ln_max = 0.0;
        deltaDDES = 0.0;
        for (iNeigh = 0;iNeigh < nNeigh; iNeigh++){
          jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
          coord_j = geometry->nodes->GetCoord(jPoint);
          for (iDim = 0; iDim < nDim; iDim++){
            delta[iDim] = fabs(coord_j[iDim] - coord_i[iDim]);
          }
          deltaDDES = geometry->nodes->GetMaxLength(iPoint);
          ln[0] = delta[1]*ratioOmega[2] - delta[2]*ratioOmega[1];
          ln[1] = delta[2]*ratioOmega[0] - delta[0]*ratioOmega[2];
          ln[2] = delta[0]*ratioOmega[1] - delta[1]*ratioOmega[0];
          aux_ln = sqrt(ln[0]*ln[0] + ln[1]*ln[1] + ln[2]*ln[2]);
          ln_max = max(ln_max,aux_ln);
          vortexTiltingMeasure += nodes->GetVortex_Tilting(jPoint);
        }

        vortexTiltingMeasure = (vortexTiltingMeasure/fabs(nNeigh + 1.0));

        f_kh = max(f_min, min(f_max, f_min + ((f_max - f_min)/(a2 - a1)) * (vortexTiltingMeasure - a1)));

        r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
        f_d = 1.0-tanh(pow(8.0*r_d,3.0));

        maxDelta = (ln_max/sqrt(3.0)) * f_kh;
        if (f_d < 0.999){
          maxDelta = deltaDDES;
        }

        distDES = constDES * maxDelta;
        lengthScale=wallDistance-f_d*max(0.0,(wallDistance-distDES));

        break;

    }

    nodes->SetDES_LengthScale(iPoint, lengthScale);

  }
  END_SU2_OMP_FOR
}

void CTurbSASolver::SetInletAtVertex(const su2double *val_inlet,
                                    unsigned short iMarker,
                                    unsigned long iVertex) {

  Inlet_TurbVars[iMarker][iVertex][0] = val_inlet[nDim+2+nDim];

}

su2double CTurbSASolver::GetInletAtVertex(su2double *val_inlet,
                                          unsigned long val_inlet_point,
                                          unsigned short val_kind_marker,
                                          string val_marker,
                                          const CGeometry *geometry,
                                          const CConfig *config) const {
  /*--- Local variables ---*/

  unsigned short iMarker;
  unsigned long iPoint, iVertex;
  su2double Area = 0.0;
  su2double Normal[3] = {0.0,0.0,0.0};

  /*--- Alias positions within inlet file for readability ---*/

  if (val_kind_marker == INLET_FLOW) {

    unsigned short position = nDim+2+nDim;

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) &&
          (config->GetMarker_All_TagBound(iMarker) == val_marker)) {

        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++){

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (iPoint == val_inlet_point) {

            /*-- Compute boundary face area for this vertex. ---*/

            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            Area = GeometryToolbox::Norm(nDim, Normal);

            /*--- Access and store the inlet variables for this vertex. ---*/

            val_inlet[position] = Inlet_TurbVars[iMarker][iVertex][0];

            /*--- Exit once we find the point. ---*/

            return Area;

          }
        }
      }
    }

  }

  /*--- If we don't find a match, then the child point is not on the
   current inlet boundary marker. Return zero area so this point does
   not contribute to the restriction operator and continue. ---*/

  return Area;

}

void CTurbSASolver::SetUniformInlet(const CConfig* config, unsigned short iMarker) {

  for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
    Inlet_TurbVars[iMarker][iVertex][0] = GetNuTilde_Inf();
  }

}

void CTurbSASolver::ComputeUnderRelaxationFactor(const CConfig *config) {

  /* Apply the turbulent under-relaxation to the SA variants. The
   SA_NEG model is more robust due to allowing for negative nu_tilde,
   so the under-relaxation is not applied to that variant. */

  if (config->GetKind_Turb_Model() == SA_NEG) return;

  /* Loop over the solution update given by relaxing the linear
   system for this nonlinear iteration. */

  su2double localUnderRelaxation =  1.00;
  const su2double allowableRatio =  0.99;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    localUnderRelaxation = 1.0;
    su2double ratio = fabs(LinSysSol[iPoint]) / (fabs(nodes->GetSolution(iPoint, 0)) + EPS);
    /* We impose a limit on the maximum percentage that the
      turbulence variables can change over a nonlinear iteration. */
    if (ratio > allowableRatio) {
      localUnderRelaxation = min(allowableRatio / ratio, localUnderRelaxation);
    }

    /* Threshold the relaxation factor in the event that there is
     a very small value. This helps avoid catastrophic crashes due
     to non-realizable states by canceling the update. */

    if (localUnderRelaxation < 1e-10) localUnderRelaxation = 0.0;

    /* Store the under-relaxation factor for this point. */

    nodes->SetUnderRelaxation(iPoint, localUnderRelaxation);

  }
  END_SU2_OMP_FOR

}