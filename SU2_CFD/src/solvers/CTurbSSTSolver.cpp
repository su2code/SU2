/*!
 * \file CTurbSSTSolver.cpp
 * \brief Main subrotuines of CTurbSSTSolver class
 * \author F. Palacios, A. Bueno
 * \version 7.3.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/solvers/CTurbSSTSolver.hpp"
#include "../../include/variables/CTurbSSTVariable.hpp"
#include "../../include/variables/CFlowVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"


CTurbSSTSolver::CTurbSSTSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
    : CTurbSolver(geometry, config, true) {
  unsigned short nLineLets;
  unsigned long iPoint;
  ifstream restart_file;
  string text_line;

  bool multizone = config->GetMultizone_Problem();

  /*--- Dimension of the problem --> dependent on the turbulence model. ---*/

  nVar = 2;
  nPrimVar = 2;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Single grid simulation ---*/

  if (iMesh == MESH_0) {

    /*--- Define some auxiliary vector related with the residual ---*/

    Residual_RMS.resize(nVar,0.0);
    Residual_Max.resize(nVar,0.0);
    Point_Max.resize(nVar,0);
    Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (SST model)." << endl;
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

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (multizone){
      Residual_BGS.resize(nVar,0.0);
      Residual_Max_BGS.resize(nVar,0.0);
      Point_Max_BGS.resize(nVar,0);
      Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
    }

  }

  /*--- Initialize value for model constants ---*/
  constants[0] = 0.85;   //sigma_k1
  constants[1] = 1.0;    //sigma_k2
  constants[2] = 0.5;    //sigma_om1
  constants[3] = 0.856;  //sigma_om2
  constants[4] = 0.075;  //beta_1
  constants[5] = 0.0828; //beta_2
  constants[6] = 0.09;   //betaStar
  constants[7] = 0.31;   //a1
  /*--- SST constants ---*/
  // constants[8] = constants[4]/constants[6] - constants[2]*0.41*0.41/sqrt(constants[6]);  //alfa_1
  // constants[9] = constants[5]/constants[6] - constants[3]*0.41*0.41/sqrt(constants[6]);  //alfa_2
  /*--- SST2003m constants---*/
  constants[8] = 5.0 / 9.0; 
  constants[9] = 0.44;

  /*--- Initialize lower and upper limits---*/
  // lowerlimit[0] = 1.0e-10;
  // upperlimit[0] = 1.0e10;

  // lowerlimit[1] = 1.0e-4;
  // upperlimit[1] = 1.0e15;
  lowerlimit[0] = 1.0e-10;
  upperlimit[0] = 1.0e15;

  lowerlimit[1] = 1.0e-10;
  upperlimit[1] = 1.0e15;

  /*--- Far-field flow state quantities and initialization. ---*/
  su2double rhoInf, *VelInf, muLamInf, Intensity, viscRatio, muT_Inf;

  rhoInf    = config->GetDensity_FreeStreamND();
  VelInf    = config->GetVelocity_FreeStreamND();
  muLamInf  = config->GetViscosity_FreeStreamND();
  Intensity = config->GetTurbulenceIntensity_FreeStream();
  viscRatio = config->GetTurb2LamViscRatio_FreeStream();

  su2double VelMag2 = GeometryToolbox::SquaredNorm(nDim, VelInf);

  su2double kine_Inf  = 3.0/2.0*(VelMag2*Intensity*Intensity);
  su2double omega_Inf = rhoInf*kine_Inf/(muLamInf*viscRatio);

  Solution_Inf[0] = kine_Inf;
  Solution_Inf[1] = omega_Inf;

  /*--- Eddy viscosity, initialized without stress limiter at the infinity ---*/
  muT_Inf = rhoInf*kine_Inf/omega_Inf;

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CTurbSSTVariable(kine_Inf, omega_Inf, muT_Inf, nPoint, nDim, nVar, constants, config);
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
    due to arbitrary number of turbulence variables ---*/

  Inlet_TurbVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_TurbVars[iMarker].resize(nVertex[iMarker],nVar);
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex) {
      Inlet_TurbVars[iMarker](iVertex,0) = kine_Inf;
      Inlet_TurbVars[iMarker](iVertex,1) = omega_Inf;
    }
  }

  /*--- Store the initial CFL number for all grid points. ---*/

  const su2double CFL = config->GetCFL(MGLevel)*config->GetCFLRedCoeff_Turb();
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "K-W SST";

}

void CTurbSSTSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
         unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);

  /*--- Upwind second order reconstruction and gradients ---*/
  CommonPreprocessing(geometry, config, Output);

  /*--- Compute primitives and gradients ---*/
  Postprocessing(geometry, solver_container, config, iMesh);
}

void CTurbSSTSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container,
                                    CConfig *config, unsigned short iMesh) {

  const su2double a1 = constants[7];

  /*--- Compute turbulence gradients. ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetSolution_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetSolution_Gradient_LS(geometry, config);
  }

  AD::StartNoSharedReading();

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

  SetPrimitive_Variables(solver_container);

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Compute blending functions and cross diffusion ---*/

    su2double rho = flowNodes->GetDensity(iPoint);
    su2double mu  = flowNodes->GetLaminarViscosity(iPoint);

    su2double dist = geometry->nodes->GetWall_Distance(iPoint);

    //su2double VorticityMag = GeometryToolbox::Norm(3, flowNodes->GetVorticity(iPoint));
    //VorticityMag = max(VorticityMag, 1e-12); // safety against division by zero

    su2double StrainMag = nodes->GetStrainMag(iPoint);

    nodes->SetBlendingFunc(iPoint, mu, dist, rho);

    su2double F2 = nodes->GetF2blending(iPoint);

    /*--- Compute the eddy viscosity ---*/

    su2double kine  = nodes->GetPrimitive(iPoint,0);
    su2double omega = nodes->GetPrimitive(iPoint,1);
    //su2double zeta  = min(1.0/omega, a1/(VorticityMag*F2));
    su2double zeta  = max(omega, (StrainMag*F2/a1));
    /*--- compute eddy viscosity according to SST2003m ---*/
    su2double muT   = max(rho*kine/zeta,0.0);

    nodes->SetmuT(iPoint,muT);

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();
}

void CTurbSSTSolver::SetPrimitive_Variables(CSolver **solver_container) {
  
  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++)
    for (auto iVar = 0u; iVar < nVar; iVar++)
      nodes->SetPrimitive(iPoint, iVar, nodes->GetSolution(iPoint, iVar)/flowNodes->GetDensity(iPoint));

}

void CTurbSSTSolver::Viscous_Residual(unsigned long iEdge, CGeometry* geometry, CSolver** solver_container,
                                     CNumerics* numerics, CConfig* config) {

  /*--- Define an object to set solver specific numerics contribution. ---*/
  auto SolverSpecificNumerics = [&](unsigned long iPoint, unsigned long jPoint) {
    /*--- Menter's first blending function (only SST)---*/
    numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(jPoint));
  };

  /*--- Now instantiate the generic implementation with the functor above. ---*/

  Viscous_Residual_impl(SolverSpecificNumerics, iEdge, geometry, solver_container, numerics, config);
}

void CTurbSSTSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  bool axisymmetric = config->GetAxisymmetric();

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

  /*--- Pick one numerics object per thread. ---*/
  auto* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Loop over all points. ---*/

  AD::StartNoSharedReading();

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), nullptr);

    /*--- Gradient of the primitive and conservative variables ---*/

    numerics->SetPrimVarGradient(flowNodes->GetGradient_Primitive(iPoint), nullptr);

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/

    numerics->SetScalarVar(nodes->GetPrimitive(iPoint), nullptr);
    numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nullptr);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Set distance to the surface ---*/

    numerics->SetDistance(geometry->nodes->GetWall_Distance(iPoint), 0.0);

    /*--- Menter's first blending function ---*/

    numerics->SetF1blending(nodes->GetF1blending(iPoint),0.0);

    /*--- Menter's second blending function ---*/

    numerics->SetF2blending(nodes->GetF2blending(iPoint),0.0);

    /*--- Set vorticity and strain rate magnitude ---*/

    numerics->SetVorticity(flowNodes->GetVorticity(iPoint), nullptr);

    numerics->SetStrainMag(flowNodes->GetStrainMag(iPoint), 0.0);

    /*--- Cross diffusion ---*/

    numerics->SetCrossDiff(nodes->GetCrossDiff(iPoint),0.0);

    if (axisymmetric){
      /*--- Set y coordinate ---*/
      numerics->SetCoord(geometry->nodes->GetCoord(iPoint), geometry->nodes->GetCoord(iPoint));
    }

    /*--- Compute the source term ---*/

    auto residual = numerics->ComputeResidual(config);

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, residual);
    if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

}

void CTurbSSTSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {
}

void CTurbSSTSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  bool rough_wall = false;
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  WALL_TYPE WallType; su2double Roughness_Height;
  tie(WallType, Roughness_Height) = config->GetWallRoughnessProperties(Marker_Tag);
  if (WallType == WALL_TYPE::ROUGH) rough_wall = true;

  /*--- Evaluate nu tilde at the closest point to the surface using the wall functions. ---*/

  if (config->GetWall_Functions()) {
    SU2_OMP_MASTER
    SetTurbVars_WF(geometry, solver_container, config, val_marker);
    END_SU2_OMP_MASTER
    SU2_OMP_BARRIER
    return;
  }

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      if (rough_wall) {

        /*--- Set wall values ---*/
        su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        su2double laminar_viscosity = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
        su2double WallShearStress = solver_container[FLOW_SOL]->GetWallShearStress(val_marker, iVertex);

        /*--- Compute non-dimensional velocity ---*/
        su2double FrictionVel = sqrt(fabs(WallShearStress)/density);

        /*--- Compute roughness in wall units. ---*/
        //su2double Roughness_Height = config->GetWall_RoughnessHeight(Marker_Tag);
        su2double kPlus = FrictionVel*Roughness_Height*density/laminar_viscosity;

        su2double S_R= 0.0;
        /*--- Reference 1 original Wilcox (1998) ---*/
        /*if (kPlus <= 25)
            S_R = (50/(kPlus+EPS))*(50/(kPlus+EPS));
          else
            S_R = 100/(kPlus+EPS);*/

        /*--- Reference 2 from D.C. Wilcox Turbulence Modeling for CFD (2006) ---*/
        if (kPlus <= 5)
          S_R = (200/(kPlus+EPS))*(200/(kPlus+EPS));
        else
          S_R = 100/(kPlus+EPS) + ((200/(kPlus+EPS))*(200/(kPlus+EPS)) - 100/(kPlus+EPS))*exp(5-kPlus);

        /*--- Modify the omega to account for a rough wall. ---*/
        su2double solution[2];
        solution[0] = 0.0;
        solution[1] = FrictionVel*FrictionVel*S_R/(laminar_viscosity/density);

        /*--- Set the solution values and zero the residual ---*/
        nodes->SetSolution_Old(iPoint,solution);
        nodes->SetSolution(iPoint,solution);
        LinSysRes.SetBlock_Zero(iPoint);

      } else { // smooth wall

        /*--- distance to closest neighbor ---*/
        const auto jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        // su2double distance2 = GeometryToolbox::SquaredDistance(nDim,
        //                                                      geometry->nodes->GetCoord(iPoint),
        //                                                      geometry->nodes->GetCoord(jPoint));

        su2double Vector[MAXNDIM] = {0.0};
        GeometryToolbox::Distance(nDim,
                                  geometry->nodes->GetCoord(iPoint),
                                  geometry->nodes->GetCoord(jPoint),
                                  Vector);

        const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

        su2double Area = GeometryToolbox::Norm(nDim, Normal);

        su2double UnitNormal[MAXNDIM] = {0.0};
        for (auto iDim = 0u; iDim < nDim; iDim++)
          UnitNormal[iDim] = -Normal[iDim]/Area;

        su2double distance  = GeometryToolbox::DotProduct(nDim, Vector, UnitNormal);
        su2double distance2 = pow(distance, 2.0);

        /*--- Set wall values ---*/

        su2double density_wall = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(jPoint);
        su2double laminar_viscosity = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(jPoint);

        su2double beta_1 = constants[4];
        su2double solution[MAXNVAR];
        solution[0] = 0.0;
        solution[1] = 60.0*laminar_viscosity*density_wall/(density*beta_1*distance2);

        /*--- Set the solution values and zero the residual ---*/
        nodes->SetSolution_Old(iPoint,solution);
        nodes->SetSolution(iPoint,solution);
        LinSysRes.SetBlock_Zero(iPoint);
      }

      if (implicit) {
        /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
        Jacobian.DeleteValsRowi(iPoint*nVar);
        Jacobian.DeleteValsRowi(iPoint*nVar+1);
      }
    }
  }
  END_SU2_OMP_FOR
}


void CTurbSSTSolver::SetTurbVars_WF(CGeometry *geometry, CSolver **solver_container,
                                    const CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- von Karman constant from boundary layer theory ---*/
  const su2double kappa = config->GetwallModel_Kappa();
  const su2double minYPlus = config->GetwallModel_MinYPlus();
  /*--- relaxation factor for k-omega values ---*/
  const su2double relax = config->GetwallModel_RelFac();

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    const auto iPoint_Neighbor = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

    su2double Y_Plus = solver_container[FLOW_SOL]->GetYPlus(val_marker, iVertex);
    su2double Lam_Visc_Wall = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);

    /*--- Do not use wall model at the ipoint when y+ < "limit", use zero flux (Neumann) conditions. ---*/

    if (Y_Plus < minYPlus) {
      /* --- Use zero flux (Neumann) conditions, i.e. nothing has to be done. --- */
      continue;
    }

    su2double Eddy_Visc = solver_container[FLOW_SOL]->GetEddyViscWall(val_marker, iVertex);
    su2double k = nodes->GetPrimitive(iPoint_Neighbor,0);
    su2double omega = nodes->GetPrimitive(iPoint_Neighbor,1);
    su2double Density_Wall = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    su2double U_Tau = solver_container[FLOW_SOL]->GetUTau(val_marker, iVertex);
    su2double y = Y_Plus*Lam_Visc_Wall/(Density_Wall*U_Tau);

    su2double omega1 = 6.0*Lam_Visc_Wall/(0.075*Density_Wall*y*y);  // eq. 19
    su2double omega0 = U_Tau/(sqrt(0.09)*kappa*y);                  // eq. 20
    su2double omega_new = sqrt(omega0*omega0 + omega1*omega1);      // eq. 21 Nichols & Nelson
    su2double k_new = omega_new * Eddy_Visc/Density_Wall;           // eq. 22 Nichols & Nelson
                                           // (is this the correct density? paper says rho and not rho_w)

    /*--- put some relaxation factor on the k-omega values ---*/
    k += relax*(k_new - k);
    omega += relax*(omega_new - omega);

    su2double solution[MAXNVAR] = {k, omega};

    nodes->SetSolution_Old(iPoint_Neighbor,solution);
    nodes->SetSolution(iPoint,solution);

    LinSysRes.SetBlock_Zero(iPoint_Neighbor);

    if (implicit) {
      /*--- includes 1 in the diagonal ---*/
      Jacobian.DeleteValsRowi(iPoint_Neighbor*nVar);
      Jacobian.DeleteValsRowi(iPoint_Neighbor*nVar+1);
    }
  }
}

void CTurbSSTSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}

void CTurbSSTSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
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

      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Set turbulent variable at the wall, and at infinity ---*/

      conv_numerics->SetScalarVar(nodes->GetPrimitive(iPoint), Solution_Inf);

      /*--- Set Normal (it is necessary to change the sign) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Grid Movement ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute residuals and Jacobians ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Add residuals and Jacobians ---*/

      LinSysRes.AddBlock(iPoint, residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
    }
  }
  END_SU2_OMP_FOR

}

void CTurbSSTSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
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
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/

      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Non-dimensionalize Inlet_TurbVars if Inlet-Files are used. ---*/
      su2double Inlet_Vars[MAXNVAR];
      Inlet_Vars[0] = Inlet_TurbVars[val_marker][iVertex][0];
      Inlet_Vars[1] = Inlet_TurbVars[val_marker][iVertex][1];
      if (config->GetInlet_Profile_From_File()) {
        Inlet_Vars[0] /= pow(config->GetVelocity_Ref(), 2);
        Inlet_Vars[1] *= config->GetViscosity_Ref() / (config->GetDensity_Ref() * pow(config->GetVelocity_Ref(), 2));
      }

      /*--- Set the turbulent variable states. Use free-stream SST
       values for the turbulent state at the inflow. ---*/
      /*--- Load the inlet turbulence variables (uniform by default). ---*/

      conv_numerics->SetScalarVar(nodes->GetPrimitive(iPoint), Inlet_Vars);

      /*--- Set various other quantities in the solver class ---*/

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

void CTurbSSTSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
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

      conv_numerics->SetScalarVar(nodes->GetPrimitive(iPoint),
                                nodes->GetPrimitive(iPoint));

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


void CTurbSSTSolver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                          CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  const auto nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (auto iSpan = 0u; iSpan < nSpanWiseSections ; iSpan++){

    su2double extAverageKine = solver_container[FLOW_SOL]->GetExtAverageKine(val_marker, iSpan);
    su2double extAverageOmega = solver_container[FLOW_SOL]->GetExtAverageOmega(val_marker, iSpan);
    su2double solution_j[] = {extAverageKine, extAverageOmega};

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

      conv_numerics->SetScalarVar(nodes->GetPrimitive(iPoint), solution_j);

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
      visc_numerics->SetScalarVar(nodes->GetPrimitive(iPoint), solution_j);
      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(iPoint));

      /*--- Compute residual, and Jacobians ---*/
      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, visc_residual);
      if (implicit) Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

}

void CTurbSSTSolver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  const auto nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- Quantities for computing the  kine and omega to impose at the inlet boundary. ---*/
  CFluidModel *FluidModel = solver_container[FLOW_SOL]->GetFluidModel();

  su2double Intensity = config->GetTurbulenceIntensity_FreeStream();
  su2double viscRatio = config->GetTurb2LamViscRatio_FreeStream();

  for (auto iSpan = 0u; iSpan < nSpanWiseSections ; iSpan++){

    /*--- Compute the inflow kine and omega using the span wise averge quntities---*/

    su2double rho       = solver_container[FLOW_SOL]->GetAverageDensity(val_marker, iSpan);
    su2double pressure  = solver_container[FLOW_SOL]->GetAveragePressure(val_marker, iSpan);
    su2double kine      = solver_container[FLOW_SOL]->GetAverageKine(val_marker, iSpan);

    FluidModel->SetTDState_Prho(pressure, rho);
    su2double muLam = FluidModel->GetLaminarViscosity();

    su2double VelMag2 = GeometryToolbox::SquaredNorm(nDim,
                          solver_container[FLOW_SOL]->GetAverageTurboVelocity(val_marker, iSpan));

    su2double kine_b  = 3.0/2.0*(VelMag2*Intensity*Intensity);
    su2double omega_b = rho*kine/(muLam*viscRatio);

    su2double solution_j[] = {kine_b, omega_b};

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

      /*--- Set the turbulent variable states. Use average span-wise values
             values for the turbulent state at the inflow. ---*/

      conv_numerics->SetScalarVar(nodes->GetPrimitive(iPoint), solution_j);

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
      visc_numerics->SetScalarVar(nodes->GetPrimitive(iPoint), solution_j);

      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(iPoint));

      /*--- Compute residual, and Jacobians ---*/
      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, visc_residual);
      if (implicit) Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

}

void CTurbSSTSolver::SetInletAtVertex(const su2double *val_inlet,
                                     unsigned short iMarker,
                                     unsigned long iVertex) {

  Inlet_TurbVars[iMarker][iVertex][0] = val_inlet[nDim+2+nDim];
  Inlet_TurbVars[iMarker][iVertex][1] = val_inlet[nDim+2+nDim+1];

}

su2double CTurbSSTSolver::GetInletAtVertex(su2double *val_inlet,
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

    unsigned short tke_position   = nDim+2+nDim;
    unsigned short omega_position = nDim+2+nDim+1;

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

            val_inlet[tke_position]   = Inlet_TurbVars[iMarker][iVertex][0];
            val_inlet[omega_position] = Inlet_TurbVars[iMarker][iVertex][1];

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

void CTurbSSTSolver::SetUniformInlet(const CConfig* config, unsigned short iMarker) {

  for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
    Inlet_TurbVars[iMarker][iVertex][0] = GetTke_Inf();
    Inlet_TurbVars[iMarker][iVertex][1] = GetOmega_Inf();
  }

}

void CTurbSSTSolver::ComputeUnderRelaxationFactor(const CConfig *config) {

  /* Loop over the solution update given by relaxing the linear
   system for this nonlinear iteration. */

  su2double localUnderRelaxation =  1.00;
  const su2double allowableRatio =  0.99;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    localUnderRelaxation = 1.0;
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      const unsigned long index = iPoint * nVar + iVar;
      su2double ratio = fabs(LinSysSol[index]) / (fabs(nodes->GetSolution(iPoint, iVar)) + EPS);
      /* We impose a limit on the maximum percentage that the
        turbulence variables can change over a nonlinear iteration. */
      if (ratio > allowableRatio) {
        localUnderRelaxation = min(allowableRatio / ratio, localUnderRelaxation);
      }
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

// TODO: Move convective SST terms here
void CTurbSSTSolver::ConvectiveError(CSolver **solver, const CGeometry *geometry, const CConfig *config,
                                     unsigned long iPoint, vector<vector<su2double> > &weights) { }

// TODO: Move viscous SST terms here
void CTurbSSTSolver::ViscousError(CSolver **solver, const CGeometry *geometry, const CConfig *config,
                                  unsigned long iPoint, vector<vector<su2double> > &weights) { }

void CTurbSSTSolver::TurbulentError(CSolver **solver, const CGeometry *geometry, const CConfig *config,
                                    unsigned long iPoint, vector<vector<su2double> > &weights) {

  CVariable *varFlo    = solver[FLOW_SOL]->GetNodes(),
            *varTur    = solver[TURB_SOL]->GetNodes(),
            *varAdjFlo = solver[ADJFLOW_SOL]->GetNodes(),
            *varAdjTur = solver[ADJTURB_SOL]->GetNodes();

  unsigned short iDim, jDim, iVar;
  const unsigned short nVarFlo = solver[FLOW_SOL]->GetnVar();
  const unsigned short nVarTur = solver[TURB_SOL]->GetnVar();
  
  //--- First-order terms (error due to viscosity)
  su2double r, u[3], k, omega,
            mu, mut,
            R, cp, g, Prt;

  r = varFlo->GetDensity(iPoint);
  u[0] = varFlo->GetVelocity(iPoint, 0);
  u[1] = varFlo->GetVelocity(iPoint, 1);
  if (nDim == 3) u[2] = varFlo->GetVelocity(iPoint, 2);
  k = varTur->GetPrimitive(iPoint, 0);
  omega = varTur->GetPrimitive(iPoint, 1);
  
  mu  = varFlo->GetLaminarViscosity(iPoint);
  mut = nodes->GetmuT(iPoint);

  g    = config->GetGamma();
  R    = config->GetGas_ConstantND();
  cp   = (g/(g-1.))*R;
  Prt  = config->GetPrandtl_Turb();

  su2double gradu[3][3], gradT[3], gradk[3], grado[3], divu, taut[3][3], tautomut[3][3],
            delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}},
            pk = 0, pw = 0;

  const su2double F1 = varTur->GetF1blending(iPoint);
  const su2double F2 = varTur->GetF2blending(iPoint);

  const su2double alfa     = F1*constants[8] + (1.0 - F1)*constants[9];
  const su2double sigmak   = F1*constants[0] + (1.0 - F1)*constants[1];
  const su2double sigmao   = F1*constants[2] + (1.0 - F1)*constants[3];
  const su2double sigmao2  = constants[3];
  const su2double beta     = F1*constants[4] + (1.0 - F1)*constants[5];
  const su2double betastar = constants[6];
  const su2double a1       = constants[7];
  const su2double CDkw     = varTur->GetCrossDiff(iPoint);

  // const su2double VorticityMag = max(GeometryToolbox::Norm(3, varFlo->GetVorticity(iPoint)), 1.0e-12);
  const su2double StrainMag = nodes->GetStrainMag(iPoint);

  // const bool stress_limited = (omega < VorticityMag*F2/a1);
  // const su2double zeta = max(omega,VorticityMag*F2/a1);
  const bool stress_limited = (omega < StrainMag*F2/a1);
  const su2double zeta = max(omega,StrainMag*F2/a1);

  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0 ; jDim < nDim; jDim++) {
      gradu[iDim][jDim] = varFlo->GetGradient_Primitive(iPoint, iDim+1, jDim);
    }
    gradT[iDim] = varFlo->GetGradient_Primitive(iPoint, 0, iDim);
    gradk[iDim] = varTur->GetGradient(iPoint, 0, iDim);
    grado[iDim] = varTur->GetGradient(iPoint, 1, iDim);
  }
  
  //--- Account for wall functions
  // su2double wf = varFlo->GetTauWallFactor(iPoint);
  su2double wf = 1.0;

  divu = 0.0; for (iDim = 0 ; iDim < nDim; ++iDim) divu += gradu[iDim][iDim];

  for (iDim = 0; iDim < nDim; ++iDim) {
    for (jDim = 0; jDim < nDim; ++jDim) {
      taut[iDim][jDim] = wf*(mut*( gradu[jDim][iDim] + gradu[iDim][jDim] )
                       - TWO3*r*k*delta[iDim][jDim]); // SST2003
                      //  - TWO3*mut*divu*delta[iDim][jDim]); // SST2003m
      tautomut[iDim][jDim] = wf*( gradu[jDim][iDim] + gradu[iDim][jDim] 
                           - TWO3*zeta*delta[iDim][jDim]); // SST2003
                          //  - TWO3*divu*delta[iDim][jDim]); // SST2003m
      pk += 1./wf*taut[iDim][jDim]*gradu[iDim][jDim];
      pw += 1./wf*tautomut[iDim][jDim]*gradu[iDim][jDim];
    }
  }

  const bool pk_limited    = (pk > 10.*betastar*r*omega*k);
  const bool pk_positive   = (pk >= 0);
  const bool pw_positive   = (pw >= 0);
  const bool cdkw_positive = (!varTur->GetCrossDiffLimited(iPoint));

  //--- Momentum weights
  vector<su2double> TmpWeights(weights[0].size(), 0.0);
  su2double factor = 0.0;
  for (iDim = 0; iDim < nDim; ++iDim) {
    factor = -TWO3*divu*alfa*varAdjTur->GetGradient_Adaptation(iPoint, 1, iDim)*pw_positive;
    // if (!pk_limited) {   
    factor += -TWO3*divu*mut/r*varAdjTur->GetGradient_Adaptation(iPoint, 0, iDim)*pk_positive;
    // }
    for (jDim = 0; jDim < nDim; ++jDim) {
      factor += (tautomut[iDim][jDim]+(gradu[iDim][jDim]+gradu[jDim][iDim]))*alfa*varAdjTur->GetGradient_Adaptation(iPoint, 1, jDim)*pw_positive;
      // if (!pk_limited) {
      factor += (taut[iDim][jDim]+mut*(gradu[iDim][jDim]+gradu[jDim][iDim]))/r*varAdjTur->GetGradient_Adaptation(iPoint, 0, jDim)*pk_positive;
      // }
    }
    TmpWeights[iDim+1] += factor;
  }

  //--- k and omega weights
  factor = 0.0;
  for (iDim = 0; iDim < nDim; ++iDim) {
    for (jDim = 0; jDim < nDim; ++jDim) {
      iVar = iDim+1;
      factor += (tautomut[iDim][jDim]+wf*TWO3*zeta*delta[iDim][jDim]) // SST2003
      // factor += (tautomut[iDim][jDim]) // SST2003m
              * (varAdjFlo->GetGradient_Adaptation(iPoint, iVar, jDim)
              + u[jDim]*varAdjFlo->GetGradient_Adaptation(iPoint, (nVarFlo-1), iDim));
    }
    factor += cp/Prt*gradT[iDim]*varAdjFlo->GetGradient_Adaptation(iPoint, (nVarFlo-1), iDim);
    factor += sigmak*gradk[iDim]*varAdjTur->GetGradient_Adaptation(iPoint, 0, iDim)
            // + sigmak*gradk[iDim]*varAdjFlo->GetGradient_Adaptation(iPoint, (nVarFlo-1), iDim)
            + sigmao*grado[iDim]*varAdjTur->GetGradient_Adaptation(iPoint, 1, iDim);
  }

  TmpWeights[nVarFlo+0] += factor/zeta;
  // TmpWeights[nVarFlo+1] += -k*factor/pow(zeta,2.)*(!stress_limited);
  TmpWeights[nVarFlo+1] += -k*factor/pow(zeta,2.);
  if (cdkw_positive) {
    for (iDim = 0; iDim < nDim; ++iDim) {
      TmpWeights[nVarFlo+0] += 2.*(1.-F1)*sigmao2/omega*grado[iDim]*varAdjTur->GetGradient_Adaptation(iPoint, 1, iDim);
      TmpWeights[nVarFlo+1] += 2.*(1.-F1)*sigmao2/omega*gradk[iDim]*varAdjTur->GetGradient_Adaptation(iPoint, 1, iDim);
    }
  }

  //--- Density weight
  for (iDim = 0; iDim < nDim; ++iDim) TmpWeights[0] += -u[iDim]*TmpWeights[iDim+1];
  TmpWeights[0] += -k*TmpWeights[nVarFlo+0] - omega*TmpWeights[nVarFlo+1]
                 // + k/zeta*factor*(!stress_limited);
                 + k/zeta*factor;

  //--- Add TmpWeights to weights, then reset for second-order terms
  for (iVar = 0; iVar < nVarFlo+nVarTur; ++iVar) weights[1][iVar] += TmpWeights[iVar];
  std::fill(TmpWeights.begin(), TmpWeights.end(), 0.0);

  //--- Second-order terms (error due to gradients)
  if(nDim == 3) {
    const unsigned short rki = 0, romegai = 1, rei = (nVarFlo - 1), xxi = 0, yyi = 3, zzi = 5;
    TmpWeights[nVarFlo+0] += -(mu+sigmak*mut)/r*(varAdjTur->GetHessian(iPoint, rki, xxi)
                                                +varAdjTur->GetHessian(iPoint, rki, yyi)
                                                +varAdjTur->GetHessian(iPoint, rki, zzi)); // Hk
                             // -(mu+mut*sigmak)/r*(varAdjFlo->GetHessian(iPoint, rei, xxi)
                             //                    +varAdjFlo->GetHessian(iPoint, rei, yyi)
                             //                    +varAdjFlo->GetHessian(iPoint, rei, zzi)); // Hk
    TmpWeights[nVarFlo+1] += -(mu+sigmao*mut)/r*(varAdjTur->GetHessian(iPoint, romegai, xxi)
                                                +varAdjTur->GetHessian(iPoint, romegai, yyi)
                                                +varAdjTur->GetHessian(iPoint, romegai, zzi)); // Homega

  }
  else {
    const unsigned short rki = 0, romegai = 1, rei = (nVarFlo - 1), xxi = 0, yyi = 2;
    TmpWeights[nVarFlo+0] += -(mu+sigmak*mut)/r*(varAdjTur->GetHessian(iPoint, rki, xxi)
                                                +varAdjTur->GetHessian(iPoint, rki, yyi)); // Hk
                             // -(mu+mut*sigmak)/r*(varAdjFlo->GetHessian(iPoint, rei, xxi)
                             //                    +varAdjFlo->GetHessian(iPoint, rei, yyi)); // Hk
    TmpWeights[nVarFlo+1] += -(mu+sigmao*mut)/r*(varAdjTur->GetHessian(iPoint, romegai, xxi)
                                                +varAdjTur->GetHessian(iPoint, romegai, yyi)); // Homega
  }
  TmpWeights[0] += -k*TmpWeights[nVarFlo+0]-omega*TmpWeights[nVarFlo+1];

  //--- Add TmpWeights to weights
  weights[2][0]         += TmpWeights[0];
  weights[2][nVarFlo+0] += TmpWeights[nVarFlo+0];
  weights[2][nVarFlo+1] += TmpWeights[nVarFlo+1];

  //--- Zeroth-order terms due to production
  // if (!pk_limited){
  weights[0][nVarFlo+0] += TWO3*divu*varAdjTur->GetSolution(iPoint,0)*pk_positive;
  // }
  // else {
  //   weights[0][0]         += 20.*betastar*k*omega*varAdjTur->GetSolution(iPoint,0);
  //   weights[0][nVarFlo+0] += -20.*betastar*omega*varAdjTur->GetSolution(iPoint,0);
  //   weights[0][nVarFlo+1] += -20.*betastar*k*varAdjTur->GetSolution(iPoint,0);
  // }
  // weights[0][nVarFlo+1] += TWO3*alfa*divu*varAdjTur->GetSolution(iPoint,1)*(!stress_limited)*pw_positive;
  weights[0][nVarFlo+1] += TWO3*alfa*divu*varAdjTur->GetSolution(iPoint,1)*pw_positive;

  //--- Zeroth-order terms due to dissipation
  weights[0][0]         += -betastar*k*omega*varAdjTur->GetSolution(iPoint,0)
                         - beta*pow(omega,2.)*varAdjTur->GetSolution(iPoint,1);
  weights[0][nVarFlo+0] += betastar*omega*varAdjTur->GetSolution(iPoint,0);
  weights[0][nVarFlo+1] += betastar*k*varAdjTur->GetSolution(iPoint,0)
                         + 2.*beta*omega*varAdjTur->GetSolution(iPoint,1);
  
  //--- Zeroth-order terms due to cross-diffusion
  weights[0][nVarFlo+1] += (1. - F1)*CDkw/(r*omega)*varAdjTur->GetSolution(iPoint,1)*cdkw_positive;

}