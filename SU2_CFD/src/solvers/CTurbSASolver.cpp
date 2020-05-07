/*!
 * \file CTurbSASolver.cpp
 * \brief Main subrotuines of CTurbSASolver class
 * \author F. Palacios, A. Bueno
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../../Common/include/omp_structure.hpp"


CTurbSASolver::CTurbSASolver(void) : CTurbSolver() { }

CTurbSASolver::CTurbSASolver(CGeometry *geometry, CConfig *config, unsigned short iMesh, CFluidModel* FluidModel)
             : CTurbSolver(geometry, config) {

  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  su2double Density_Inf, Viscosity_Inf, Factor_nu_Inf, Factor_nu_Engine, Factor_nu_ActDisk;

  bool multizone = config->GetMultizone_Problem();

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

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

    Residual = new su2double[nVar]();
    Residual_RMS = new su2double[nVar]();
    Residual_i = new su2double[nVar]();
    Residual_j = new su2double[nVar]();
    Residual_Max = new su2double[nVar]();

    /*--- Define some structures for locating max residuals ---*/

    Point_Max = new unsigned long[nVar]();
    Point_Max_Coord = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord[iVar] = new su2double[nDim]();
    }

    /*--- Define some auxiliar vector related with the solution ---*/

    Solution = new su2double[nVar];
    Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];

    /*--- Jacobians and vector structures for implicit computations ---*/

    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (SA model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);

    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

    if (ReducerStrategy)
      EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);

    if (config->GetExtraOutput()) {
      if (nDim == 2) { nOutputVariables = 13; }
      else if (nDim == 3) { nOutputVariables = 19; }
      OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
      OutputHeadingNames = new string[nOutputVariables];
    }

    /*--- Computation of gradients by least squares ---*/

    if (config->GetLeastSquaresRequired()) {
      /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
      Smatrix = new su2double* [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Smatrix[iDim] = new su2double [nDim];

      /*--- c vector := transpose(WA)*(Wb) ---*/
      Cvector = new su2double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++)
        Cvector[iVar] = new su2double [nDim];
    }

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (multizone){
      Residual_BGS      = new su2double[nVar]();
      Residual_Max_BGS  = new su2double[nVar]();

      /*--- Define some structures for locating max residuals ---*/

      Point_Max_BGS       = new unsigned long[nVar]();
      Point_Max_Coord_BGS = new su2double*[nVar];
      for (iVar = 0; iVar < nVar; iVar++) {
        Point_Max_Coord_BGS[iVar] = new su2double[nDim] ();
      }
    }

  }

  /*--- Read farfield conditions from config ---*/

  Density_Inf   = config->GetDensity_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();

  /*--- Factor_nu_Inf in [3.0, 5.0] ---*/

  Factor_nu_Inf = config->GetNuFactor_FreeStream();
  nu_tilde_Inf  = Factor_nu_Inf*Viscosity_Inf/Density_Inf;
  if (config->GetKind_Trans_Model() == BC) {
    nu_tilde_Inf  = 0.005*Factor_nu_Inf*Viscosity_Inf/Density_Inf;
  }

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

  unsigned long iMarker;

  SlidingState       = new su2double*** [nMarker] ();
  SlidingStateNodes  = new int*         [nMarker] ();

  for (iMarker = 0; iMarker < nMarker; iMarker++){

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){

      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)] ();
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)] ();

      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++)
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1] ();
    }

  }

  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
   * due to arbitrary number of turbulence variables ---*/

  Inlet_TurbVars = new su2double**[nMarker];
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_TurbVars[iMarker] = new su2double*[nVertex[iMarker]];
    for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
      Inlet_TurbVars[iMarker][iVertex] = new su2double[nVar] ();
      Inlet_TurbVars[iMarker][iVertex][0] = nu_tilde_Inf;
    }
  }

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

CTurbSASolver::~CTurbSASolver(void) {

  unsigned long iMarker, iVertex;
  unsigned short iVar;

  if ( SlidingState != nullptr ) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if ( SlidingState[iMarker] != nullptr ) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          if ( SlidingState[iMarker][iVertex] != nullptr ){
            for (iVar = 0; iVar < nPrimVar+1; iVar++)
              delete [] SlidingState[iMarker][iVertex][iVar];
            delete [] SlidingState[iMarker][iVertex];
          }
        delete [] SlidingState[iMarker];
      }
    }
    delete [] SlidingState;
  }

  if ( SlidingStateNodes != nullptr ){
    for (iMarker = 0; iMarker < nMarker; iMarker++){
      if (SlidingStateNodes[iMarker] != nullptr)
        delete [] SlidingStateNodes[iMarker];
    }
    delete [] SlidingStateNodes;
  }

}

void CTurbSASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
        unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  bool limiter_turb = (config->GetKind_SlopeLimit_Turb() != NO_LIMITER) &&
                      (config->GetInnerIter() <= config->GetLimiterIter());
  unsigned short kind_hybridRANSLES = config->GetKind_HybridRANSLES();
  const su2double* const* PrimGrad_Flow = nullptr;
  const su2double* Vorticity = nullptr;
  su2double Laminar_Viscosity = 0.0;

  /*--- Clear residual and system matrix, not needed for
   * reducer strategy as we write over the entire matrix. ---*/
  if (!ReducerStrategy) {
    LinSysRes.SetValZero();
    Jacobian.SetValZero();
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

  if (limiter_turb) SetSolution_Limiter(geometry, config);

  if (kind_hybridRANSLES != NO_HYBRIDRANSLES) {

    /*--- Set the vortex tilting coefficient at every node if required ---*/

    if (kind_hybridRANSLES == SA_EDDES){
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        Vorticity = solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint);
        PrimGrad_Flow = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);
        Laminar_Viscosity  = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
        nodes->SetVortex_Tilting(iPoint, PrimGrad_Flow, Vorticity, Laminar_Viscosity);
      }
    }

    /*--- Compute the DES length scale ---*/

    SetDES_LengthScale(solver_container, geometry, config);

  }

}

void CTurbSASolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

  const su2double cv1_3 = 7.1*7.1*7.1;

  const bool neg_spalart_allmaras = (config->GetKind_Turb_Model() == SA_NEG);

  /*--- Compute eddy viscosity ---*/

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    su2double rho = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    su2double mu  = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);

    su2double nu  = mu/rho;
    su2double nu_hat = nodes->GetSolution(iPoint,0);

    su2double Ji   = nu_hat/nu;
    su2double Ji_3 = Ji*Ji*Ji;
    su2double fv1  = Ji_3/(Ji_3+cv1_3);

    su2double muT = rho*fv1*nu_hat;

    if (neg_spalart_allmaras) muT = max(muT,0.0);

    nodes->SetmuT(iPoint,muT);

  }

}

void CTurbSASolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                    CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  const bool harmonic_balance = (config->GetTime_Marching() == HARMONIC_BALANCE);
  const bool transition    = (config->GetKind_Trans_Model() == LM);
  const bool transition_BC = (config->GetKind_Trans_Model() == BC);

  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

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

    numerics->SetTurbVar(nodes->GetSolution(iPoint), nullptr);
    numerics->SetTurbVarGradient(nodes->GetGradient(iPoint), nullptr);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->node[iPoint]->GetVolume());

    /*--- Get Hybrid RANS/LES Type and set the appropriate wall distance ---*/

    if (config->GetKind_HybridRANSLES() == NO_HYBRIDRANSLES) {

      /*--- Set distance to the surface ---*/

      numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);

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

    Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

  }

  if (harmonic_balance) {

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

      su2double Volume = geometry->node[iPoint]->GetVolume();

      /*--- Access stored harmonic balance source term ---*/

      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        su2double Source = nodes->GetHarmonicBalance_Source(iPoint,iVar);
        LinSysRes(iPoint,iVar) += Source*Volume;
      }
    }
  }

}

void CTurbSASolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                    CConfig *config, unsigned short iMesh) {
}

void CTurbSASolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iVar;

  /*--- The dirichlet condition is used only without wall function, otherwise the
   convergence is compromised as we are providing nu tilde values for the
   first point of the wall  ---*/

  if (!config->GetWall_Functions()) {

    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
      iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

      /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

      if (geometry->node[iPoint]->GetDomain()) {

        /*--- Get the velocity vector ---*/

        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = 0.0;

        nodes->SetSolution_Old(iPoint,Solution);
        LinSysRes.SetBlock_Zero(iPoint);

        /*--- Includes 1 in the diagonal ---*/

        Jacobian.DeleteValsRowi(iPoint);
      }
    }
  }
  else {

    /*--- Evaluate nu tilde at the closest point to the surface using the wall functions ---*/

    SetNuTilde_WF(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

  }

}

void CTurbSASolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                       CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iVar;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Get the velocity vector ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] = 0.0;

      nodes->SetSolution_Old(iPoint,Solution);
      LinSysRes.SetBlock_Zero(iPoint);

      /*--- Includes 1 in the diagonal ---*/

      Jacobian.DeleteValsRowi(iPoint);
    }
  }

}

void CTurbSASolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                 CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  su2double *Normal, *V_infty, *V_domain;

  Normal = new su2double[nDim];

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Allocate the value at the infinity ---*/

      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Grid Movement ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Set turbulent variable at the wall, and at infinity ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);
      Solution_j[0] = nu_tilde_Inf;
      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set Normal (it is necessary to change the sign) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Compute residuals and Jacobians ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Add residuals and Jacobians ---*/

      LinSysRes.AddBlock(iPoint, residual);
      Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }

  delete [] Normal;

}

void CTurbSASolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                             CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double *V_inlet, *V_domain, *Normal;

  Normal = new su2double[nDim];

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Allocate the value at the inlet ---*/

      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      Solution_i[0] = nodes->GetSolution(iPoint,0);

      /*--- Load the inlet turbulence variable (uniform by default). ---*/

      Solution_j[0] = Inlet_TurbVars[val_marker][iVertex][0];

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the conv_numerics class ---*/

      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_inlet);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetTurbVar(Solution_i, Solution_j);
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
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

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}

void CTurbSASolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                              CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  su2double *V_outlet, *V_domain, *Normal;

  Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Allocate the value at the outlet ---*/

      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_outlet);

      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> TurbVar_internal,
       Solution_j --> TurbVar_outlet ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);
        Solution_j[iVar] = nodes->GetSolution(iPoint,iVar);
      }
      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set Normal (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_outlet);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetTurbVar(Solution_i, Solution_j);
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
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

  /*--- Free locally allocated memory ---*/

  delete[] Normal;

}

void CTurbSASolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iDim;
  su2double *V_inflow, *V_domain, *Normal;

  Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Allocate the value at the infinity ---*/

      V_inflow = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inflow);

      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual. ---*/

      conv_numerics->SetTurbVar(nodes->GetSolution(iPoint), nodes->GetSolution(iPoint));

      /*--- Set Normal (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Set grid movement ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_inflow);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetTurbVar(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
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

  /*--- Free locally allocated memory ---*/

  delete[] Normal;

}

void CTurbSASolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double *V_exhaust, *V_domain, *Normal;

  Normal = new su2double[nDim];

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Allocate the value at the infinity ---*/

      V_exhaust = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_exhaust);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      Solution_i[0] = nodes->GetSolution(iPoint,0);
      Solution_j[0] = nu_tilde_Engine;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the conv_numerics class ---*/

      conv_numerics->SetNormal(Normal);

      /*--- Set grid movement ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_exhaust);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetTurbVar(Solution_i, Solution_j);
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
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

  /*--- Free locally allocated memory ---*/

  delete[] Normal;

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

  unsigned long iPoint, iVertex, GlobalIndex_donor, GlobalIndex;
  su2double *V_outlet, *V_inlet, *V_domain, *Normal, *UnitNormal, Area, Vn;
  bool ReverseFlow;
  unsigned short iDim;

  Normal = new su2double[nDim];
  UnitNormal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex_donor = solver_container[FLOW_SOL]->GetDonorGlobalIndex(val_marker, iVertex);
    GlobalIndex = geometry->node[iPoint]->GetGlobalIndex();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex != GlobalIndex_donor)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Check the flow direction. Project the flow into the normal to the inlet face ---*/

      Vn = 0.0; ReverseFlow = false;
      for (iDim = 0; iDim < nDim; iDim++) {  Vn += V_domain[iDim+1]*UnitNormal[iDim]; }

      if ((val_inlet_surface) && (Vn < 0.0)) { ReverseFlow = true; }
      if ((!val_inlet_surface) && (Vn > 0.0)) { ReverseFlow = true; }

      /*--- Do not anything if there is a
       reverse flow, Euler b.c. for the direct problem ---*/

      if (!ReverseFlow) {

        /*--- Allocate the value at the infinity ---*/

        if (val_inlet_surface) {
          V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
          V_outlet = solver_container[FLOW_SOL]->GetDonorPrimVar(val_marker, iVertex);
          conv_numerics->SetPrimitive(V_domain, V_inlet);
        }
        else {
          V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
          V_inlet = solver_container[FLOW_SOL]->GetDonorPrimVar(val_marker, iVertex);
          conv_numerics->SetPrimitive(V_domain, V_outlet);
        }

        /*--- Set the turb. variable solution
         set  the turbulent variables. Here we use a Neumann BC such
         that the turbulent variable is copied from the interior of the
         domain to the outlet before computing the residual.
         or set the turbulent variable states (prescribed for an inflow)  ----*/

        Solution_i[0] = nodes->GetSolution(iPoint,0);

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

        /*--- Inflow analysis (interior extrapolation) ---*/
        if (((val_inlet_surface) && (!ReverseFlow)) || ((!val_inlet_surface) && (ReverseFlow))) {
          Solution_j[0] = nodes->GetSolution(iPoint,0);
        }

        /*--- Outflow analysis ---*/
        else {
          Solution_j[0] = nu_tilde_ActDisk;
        }

        conv_numerics->SetTurbVar(Solution_i, Solution_j);

        /*--- Grid Movement ---*/

        if (dynamic_grid)
          conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

        /*--- Compute the residual using an upwind scheme ---*/

        auto residual = conv_numerics->ComputeResidual(config);
        LinSysRes.AddBlock(iPoint, residual);

        /*--- Jacobian contribution for implicit integration ---*/

        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

//        /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//        visc_numerics->SetNormal(Normal);
//        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint_Normal]->GetCoord());
//
//        /*--- Conservative variables w/o reconstruction ---*/
//
//        if (val_inlet_surface) visc_numerics->SetPrimitive(V_domain, V_inlet);
//        else visc_numerics->SetPrimitive(V_domain, V_outlet);
//
//        /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//        visc_numerics->SetTurbVar(Solution_i, Solution_j);
//
//        visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
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
    }
  }

  /*--- Free locally allocated memory ---*/

  delete[] Normal;
  delete[] UnitNormal;

}

void CTurbSASolver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                         CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim, iSpan;
  unsigned long  oldVertex, iPoint, Point_Normal, iVertex;
  su2double *V_inlet, *V_domain, *Normal;
  su2double extAverageNu;
  Normal = new su2double[nDim];

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  unsigned short nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iSpan= 0; iSpan < nSpanWiseSections ; iSpan++){
    extAverageNu = solver_container[FLOW_SOL]->GetExtAverageNu(val_marker, iSpan);

    /*--- Loop over all the vertices on this boundary marker ---*/

    for (iVertex = 0; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {

      /*--- find the node related to the vertex ---*/
      iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][oldVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Allocate the value at the inlet ---*/
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      Solution_i[0] = nodes->GetSolution(iPoint,0);
      Solution_j[0] = extAverageNu;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the conv_numerics class ---*/

      conv_numerics->SetNormal(Normal);

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the conv_numerics class ---*/

      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
            geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      auto conv_residual = conv_numerics->ComputeResidual(config);

      /*--- Jacobian contribution for implicit integration ---*/

      LinSysRes.AddBlock(iPoint, conv_residual);
      Jacobian.AddBlock2Diag(iPoint, conv_residual.jacobian_i);

      /*--- Viscous contribution ---*/

      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/

      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/

      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Compute residual, and Jacobians ---*/

      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/

      LinSysRes.SubtractBlock(iPoint, visc_residual);
      Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}

void CTurbSASolver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                   CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iSpan;
  unsigned long  oldVertex, iPoint, Point_Normal, iVertex;
  su2double *V_inlet, *V_domain, *Normal;

  su2double rho, pressure, muLam, Factor_nu_Inf, nu_tilde;
  Normal = new su2double[nDim];

  unsigned short nSpanWiseSections = config->GetnSpanWiseSections();
  CFluidModel *FluidModel;

  FluidModel = solver_container[FLOW_SOL]->GetFluidModel();
  Factor_nu_Inf = config->GetNuFactor_FreeStream();


  /*--- Loop over all the spans on this boundary marker ---*/
  for (iSpan= 0; iSpan < nSpanWiseSections ; iSpan++){
    rho       = solver_container[FLOW_SOL]->GetAverageDensity(val_marker, iSpan);
    pressure  = solver_container[FLOW_SOL]->GetAveragePressure(val_marker, iSpan);

    FluidModel->SetTDState_Prho(pressure, rho);
    muLam = FluidModel->GetLaminarViscosity();

    nu_tilde  = Factor_nu_Inf*muLam/rho;


    /*--- Loop over all the vertices on this boundary marker ---*/
    for (iVertex = 0; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {

      /*--- find the node related to the vertex ---*/
      iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][oldVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Allocate the value at the inlet ---*/
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      Solution_i[0] = nodes->GetSolution(iPoint,0);
      Solution_j[0] =  nu_tilde;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the conv_numerics class ---*/

      conv_numerics->SetNormal(Normal);

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the conv_numerics class ---*/

      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
            geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      auto conv_residual = conv_numerics->ComputeResidual(config);

      /*--- Jacobian contribution for implicit integration ---*/

      LinSysRes.AddBlock(iPoint, conv_residual);
      Jacobian.AddBlock2Diag(iPoint, conv_residual.jacobian_i);

      /*--- Viscous contribution ---*/

      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/

      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/

      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Compute residual, and Jacobians ---*/

      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/

      LinSysRes.SubtractBlock(iPoint, visc_residual);
      Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}

void CTurbSASolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                          CConfig *config, unsigned short val_marker) {

  //  unsigned long iVertex, iPoint, jPoint;
  //  unsigned short iVar, iDim;
  //
  //  su2double *Vector = new su2double[nDim];
  //
  //#ifndef HAVE_MPI
  //
  //  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //
  //    if (geometry->node[iPoint]->GetDomain()) {
  //
  //      /*--- Find the associate pair to the original node ---*/
  //      jPoint = geometry->vertex[val_marker][iVertex]->GetDonorPoint();
  //
  //      if (iPoint != jPoint) {
  //
  //        /*--- Store the solution for both points ---*/
  //        for (iVar = 0; iVar < nVar; iVar++) {
  //          Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);
  //          Solution_j[iVar] = nodes->GetSolution(jPoint,iVar);
  //        }
  //
  //        /*--- Set Conservative Variables ---*/
  //        numerics->SetTurbVar(Solution_i, Solution_j);
  //
  //        /*--- Retrieve flow solution for both points ---*/
  //        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
  //          FlowPrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
  //          FlowPrimVar_j[iVar] = solver_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
  //        }
  //
  //        /*--- Set Flow Variables ---*/
  //        numerics->SetConservative(FlowPrimVar_i, FlowPrimVar_j);
  //
  //        /*--- Set the normal vector ---*/
  //        geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
  //        for (iDim = 0; iDim < nDim; iDim++)
  //          Vector[iDim] = -Vector[iDim];
  //        numerics->SetNormal(Vector);
  //
  //        /*--- Add Residuals and Jacobians ---*/
  //        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
  //        LinSysRes.AddBlock(iPoint, Residual);
  //        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
  //
  //      }
  //    }
  //  }
  //
  //#else
  //
  //  int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
  //  su2double *Conserv_Var, *Flow_Var;
  //  bool compute;
  //
  //  unsigned short Buffer_Size = nVar+solver_container[FLOW_SOL]->GetnVar();
  //  su2double *Buffer_Send_U = new su2double [Buffer_Size];
  //  su2double *Buffer_Receive_U = new su2double [Buffer_Size];
  //
  //  /*--- Do the send process, by the moment we are sending each
  //   node individually, this must be changed ---*/
  //  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //    if (geometry->node[iPoint]->GetDomain()) {
  //
  //      /*--- Find the associate pair to the original node ---*/
  //      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
  //      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
  //
  //      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
  //      else compute = true;
  //
  //      /*--- We only send the information that belong to other boundary ---*/
  //      if ((jProcessor != rank) && compute) {
  //
  //        Conserv_Var = nodes->GetSolution(iPoint);
  //        Flow_Var = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
  //
  //        for (iVar = 0; iVar < nVar; iVar++)
  //          Buffer_Send_U[iVar] = Conserv_Var[iVar];
  //
  //        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
  //          Buffer_Send_U[nVar+iVar] = Flow_Var[iVar];
  //
  //        MPI::COMM_WORLD.Bsend(Buffer_Send_U, Buffer_Size, MPI::DOUBLE, jProcessor, iPoint);
  //
  //      }
  //    }
  //  }
  //
  //  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //
  //    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //
  //    if (geometry->node[iPoint]->GetDomain()) {
  //
  //      /*--- Find the associate pair to the original node ---*/
  //      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
  //      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
  //
  //      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
  //      else compute = true;
  //
  //      if (compute) {
  //
  //        /*--- We only receive the information that belong to other boundary ---*/
  //        if (jProcessor != rank) {
  //          MPI::COMM_WORLD.Recv(Buffer_Receive_U, Buffer_Size, MPI::DOUBLE, jProcessor, jPoint);
  //        }
  //        else {
  //
  //          for (iVar = 0; iVar < nVar; iVar++)
  //            Buffer_Receive_U[iVar] = nodes->GetSolution(jPoint,iVar);
  //
  //          for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
  //            Buffer_Send_U[nVar+iVar] = solver_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
  //
  //        }
  //
  //        /*--- Store the solution for both points ---*/
  //        for (iVar = 0; iVar < nVar; iVar++) {
  //          Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);
  //          Solution_j[iVar] = Buffer_Receive_U[iVar];
  //        }
  //
  //        /*--- Set Turbulent Variables ---*/
  //        numerics->SetTurbVar(Solution_i, Solution_j);
  //
  //        /*--- Retrieve flow solution for both points ---*/
  //        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
  //          FlowPrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
  //          FlowPrimVar_j[iVar] = Buffer_Receive_U[nVar + iVar];
  //        }
  //
  //        /*--- Set Flow Variables ---*/
  //        numerics->SetConservative(FlowPrimVar_i, FlowPrimVar_j);
  //
  //        geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
  //        for (iDim = 0; iDim < nDim; iDim++)
  //          Vector[iDim] = -Vector[iDim];
  //        numerics->SetNormal(Vector);
  //
  //        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
  //        LinSysRes.AddBlock(iPoint, Residual);
  //        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
  //
  //      }
  //    }
  //  }
  //
  //  delete[] Buffer_Send_U;
  //  delete[] Buffer_Receive_U;
  //
  //#endif
  //
  //  delete[] Vector;
  //
}

void CTurbSASolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                          CConfig *config, unsigned short val_marker) {

  //  unsigned long iVertex, iPoint, jPoint;
  //  unsigned short iVar, iDim;
  //
  //  su2double *Vector = new su2double[nDim];
  //
  //#ifndef HAVE_MPI
  //
  //  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //
  //    if (geometry->node[iPoint]->GetDomain()) {
  //
  //      /*--- Find the associate pair to the original node ---*/
  //      jPoint = geometry->vertex[val_marker][iVertex]->GetDonorPoint();
  //
  //      if (iPoint != jPoint) {
  //
  //        /*--- Store the solution for both points ---*/
  //        for (iVar = 0; iVar < nVar; iVar++) {
  //          Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);
  //          Solution_j[iVar] = nodes->GetSolution(jPoint,iVar);
  //        }
  //
  //        /*--- Set Conservative Variables ---*/
  //        numerics->SetTurbVar(Solution_i, Solution_j);
  //
  //        /*--- Retrieve flow solution for both points ---*/
  //        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
  //          FlowPrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
  //          FlowPrimVar_j[iVar] = solver_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
  //        }
  //
  //        /*--- Set Flow Variables ---*/
  //        numerics->SetConservative(FlowPrimVar_i, FlowPrimVar_j);
  //
  //        /*--- Set the normal vector ---*/
  //        geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
  //        for (iDim = 0; iDim < nDim; iDim++)
  //          Vector[iDim] = -Vector[iDim];
  //        numerics->SetNormal(Vector);
  //
  //        /*--- Add Residuals and Jacobians ---*/
  //        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
  //        LinSysRes.AddBlock(iPoint, Residual);
  //        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
  //
  //      }
  //    }
  //  }
  //
  //#else
  //
  //  int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
  //  su2double *Conserv_Var, *Flow_Var;
  //  bool compute;
  //
  //  unsigned short Buffer_Size = nVar+solver_container[FLOW_SOL]->GetnVar();
  //  su2double *Buffer_Send_U = new su2double [Buffer_Size];
  //  su2double *Buffer_Receive_U = new su2double [Buffer_Size];
  //
  //  /*--- Do the send process, by the moment we are sending each
  //   node individually, this must be changed ---*/
  //  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //    if (geometry->node[iPoint]->GetDomain()) {
  //
  //      /*--- Find the associate pair to the original node ---*/
  //      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
  //      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
  //
  //      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
  //      else compute = true;
  //
  //      /*--- We only send the information that belong to other boundary ---*/
  //      if ((jProcessor != rank) && compute) {
  //
  //        Conserv_Var = nodes->GetSolution(iPoint);
  //        Flow_Var = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
  //
  //        for (iVar = 0; iVar < nVar; iVar++)
  //          Buffer_Send_U[iVar] = Conserv_Var[iVar];
  //
  //        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
  //          Buffer_Send_U[nVar+iVar] = Flow_Var[iVar];
  //
  //        MPI::COMM_WORLD.Bsend(Buffer_Send_U, Buffer_Size, MPI::DOUBLE, jProcessor, iPoint);
  //
  //      }
  //    }
  //  }
  //
  //  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //
  //    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //
  //    if (geometry->node[iPoint]->GetDomain()) {
  //
  //      /*--- Find the associate pair to the original node ---*/
  //      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
  //      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
  //
  //      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
  //      else compute = true;
  //
  //      if (compute) {
  //
  //        /*--- We only receive the information that belong to other boundary ---*/
  //        if (jProcessor != rank) {
  //          MPI::COMM_WORLD.Recv(Buffer_Receive_U, Buffer_Size, MPI::DOUBLE, jProcessor, jPoint);
  //        }
  //        else {
  //
  //          for (iVar = 0; iVar < nVar; iVar++)
  //            Buffer_Receive_U[iVar] = nodes->GetSolution(jPoint,iVar);
  //
  //          for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
  //            Buffer_Send_U[nVar+iVar] = solver_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
  //
  //        }
  //
  //        /*--- Store the solution for both points ---*/
  //        for (iVar = 0; iVar < nVar; iVar++) {
  //          Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);
  //          Solution_j[iVar] = Buffer_Receive_U[iVar];
  //        }
  //
  //        /*--- Set Turbulent Variables ---*/
  //        numerics->SetTurbVar(Solution_i, Solution_j);
  //
  //        /*--- Retrieve flow solution for both points ---*/
  //        for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
  //          FlowPrimVar_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
  //          FlowPrimVar_j[iVar] = Buffer_Receive_U[nVar + iVar];
  //        }
  //
  //        /*--- Set Flow Variables ---*/
  //        numerics->SetConservative(FlowPrimVar_i, FlowPrimVar_j);
  //
  //        geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
  //        for (iDim = 0; iDim < nDim; iDim++)
  //          Vector[iDim] = -Vector[iDim];
  //        numerics->SetNormal(Vector);
  //
  //        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
  //        LinSysRes.AddBlock(iPoint, Residual);
  //        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
  //
  //      }
  //    }
  //  }
  //
  //  delete[] Buffer_Send_U;
  //  delete[] Buffer_Receive_U;
  //
  //#endif
  //
  //  delete[] Vector;
  //
}

void CTurbSASolver::SetNuTilde_WF(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                  CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Local variables ---*/

  unsigned short iDim, jDim, iVar, iNode;
  unsigned long iVertex, iPoint, iPoint_Neighbor, counter;

  su2double func, func_prim;
  su2double *Normal, Area;
  su2double div_vel, UnitNormal[3];
  su2double **grad_primvar, tau[3][3];
  su2double Vel[3], VelNormal, VelTang[3], VelTangMod, VelInfMod, WallDist[3], WallDistMod;
  su2double Lam_Visc_Normal, Kin_Visc_Normal, dypw_dyp, Eddy_Visc, nu_til_old, nu_til, cv1_3;
  su2double T_Normal, P_Normal, Density_Normal;
  su2double Density_Wall, T_Wall, P_Wall, Lam_Visc_Wall, Tau_Wall, Tau_Wall_Old;
  su2double *Coord, *Coord_Normal;
  su2double diff, Delta;
  su2double U_Tau, U_Plus = 0.0, Gam = 0.0, Beta = 0.0, Phi, Q = 0.0, Y_Plus_White = 0.0, Y_Plus;
  su2double TauElem[3], TauNormal, TauTangent[3], WallShearStress;
  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

  unsigned short max_iter = 100;
  su2double tol = 1e-10;

  /*--- Get the freestream velocity magnitude for non-dim. purposes ---*/

  su2double *VelInf = config->GetVelocity_FreeStreamND();
  VelInfMod = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    VelInfMod += VelInf[iDim];
  VelInfMod = sqrt(VelInfMod);

  /*--- Compute the recovery factor ---*/
  // su2double-check: laminar or turbulent Pr for this?
  su2double Recovery = pow(config->GetPrandtl_Lam(),(1.0/3.0));

  /*--- Typical constants from boundary layer theory ---*/

  su2double kappa = 0.4;
  su2double B = 5.5;

  /*--- Identify the boundary by string name ---*/

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall heat flux from config ---*/

  // Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- We can use also GetNormal_Neighbor, and eliminate the following loop ---*/

    iPoint_Neighbor = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

    for(iNode = 0; iNode < geometry->node[iPoint]->GetnPoint(); iNode++) {
      iPoint_Neighbor = geometry->node[iPoint]->GetPoint(iNode);

      /*--- Check if the node belongs to the domain (i.e, not a halo node)
       and the neighbor is not part of the physical boundary ---*/

      if (geometry->node[iPoint]->GetDomain() && (!geometry->node[iPoint_Neighbor]->GetBoundary())) {

        /*--- Get coordinates of the current vertex and nearest normal point ---*/

        Coord = geometry->node[iPoint]->GetCoord();
        Coord_Normal = geometry->node[iPoint_Neighbor]->GetCoord();

        /*--- Compute dual-grid area and boundary normal ---*/

        Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Area += Normal[iDim]*Normal[iDim];
        Area = sqrt (Area);

        for (iDim = 0; iDim < nDim; iDim++)
          UnitNormal[iDim] = -Normal[iDim]/Area;

        /*--- Get the velocity, pressure, and temperature at the nearest
         (normal) interior point. ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          Vel[iDim]    = solver_container[FLOW_SOL]->GetNodes()->GetVelocity(iPoint_Neighbor,iDim);
        P_Normal       = solver_container[FLOW_SOL]->GetNodes()->GetPressure(iPoint_Neighbor);
        T_Normal       = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint_Neighbor);

        /*--- Compute the wall-parallel velocity at first point off the wall ---*/

        VelNormal = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          VelNormal += Vel[iDim] * UnitNormal[iDim];
        for (iDim = 0; iDim < nDim; iDim++)
          VelTang[iDim] = Vel[iDim] - VelNormal*UnitNormal[iDim];

        VelTangMod = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          VelTangMod += VelTang[iDim]*VelTang[iDim];
        VelTangMod = sqrt(VelTangMod);

        /*--- Compute normal distance of the interior point from the wall ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);

        WallDistMod = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          WallDistMod += WallDist[iDim]*WallDist[iDim];
        WallDistMod = sqrt(WallDistMod);

        /*--- Compute mach number ---*/

        // M_Normal = VelTangMod / sqrt(Gamma * Gas_Constant * T_Normal);

        /*--- Compute the wall temperature using the Crocco-Buseman equation ---*/

        //T_Wall = T_Normal * (1.0 + 0.5*Gamma_Minus_One*Recovery*M_Normal*M_Normal);
        T_Wall = T_Normal + Recovery*pow(VelTangMod,2.0)/(2.0*Cp);

        /*--- Extrapolate the pressure from the interior & compute the
         wall density using the equation of state ---*/

        P_Wall = P_Normal;
        Density_Wall = P_Wall/(Gas_Constant*T_Wall);

        /*--- Compute the shear stress at the wall in the regular fashion
         by using the stress tensor on the surface ---*/

        Lam_Visc_Wall = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
        grad_primvar  = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);

        div_vel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          div_vel += grad_primvar[iDim+1][iDim];

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Delta = 0.0; if (iDim == jDim) Delta = 1.0;
            tau[iDim][jDim] = Lam_Visc_Wall*(  grad_primvar[jDim+1][iDim]
                                             + grad_primvar[iDim+1][jDim]) -
            TWO3*Lam_Visc_Wall*div_vel*Delta;
          }
          TauElem[iDim] = 0.0;
          for (jDim = 0; jDim < nDim; jDim++)
            TauElem[iDim] += tau[iDim][jDim]*UnitNormal[jDim];
        }

        /*--- Compute wall shear stress as the magnitude of the wall-tangential
         component of the shear stress tensor---*/

        TauNormal = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          TauNormal += TauElem[iDim] * UnitNormal[iDim];

        for (iDim = 0; iDim < nDim; iDim++)
          TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];

        WallShearStress = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          WallShearStress += TauTangent[iDim]*TauTangent[iDim];
        WallShearStress = sqrt(WallShearStress);

        /*--- Calculate the quantities from boundary layer theory and
         iteratively solve for a new wall shear stress. Use the current wall
         shear stress as a starting guess for the wall function. ---*/

        Tau_Wall_Old = WallShearStress;
        counter = 0; diff = 1.0;

        while (diff > tol) {

          /*--- Friction velocity and u+ ---*/

          U_Tau = sqrt(Tau_Wall_Old/Density_Wall);
          U_Plus = VelTangMod/U_Tau;

          /*--- Gamma, Beta, Q, and Phi, defined by Nichols & Nelson (2004) ---*/

          Gam  = Recovery*U_Tau*U_Tau/(2.0*Cp*T_Wall);
          Beta = 0.0; // For adiabatic flows only
          Q    = sqrt(Beta*Beta + 4.0*Gam);
          Phi  = asin(-1.0*Beta/Q);

          /*--- Y+ defined by White & Christoph (compressibility and heat transfer) ---*/

          Y_Plus_White = exp((kappa/sqrt(Gam))*(asin((2.0*Gam*U_Plus - Beta)/Q) - Phi))*exp(-1.0*kappa*B);

          /*--- Spalding's universal form for the BL velocity with the
           outer velocity form of White & Christoph above. ---*/

          Y_Plus = U_Plus + Y_Plus_White - (exp(-1.0*kappa*B)*
                                            (1.0 + kappa*U_Plus + kappa*kappa*U_Plus*U_Plus/2.0 +
                                             kappa*kappa*kappa*U_Plus*U_Plus*U_Plus/6.0));

          /*--- Calculate an updated value for the wall shear stress
           using the y+ value, the definition of y+, and the definition of
           the friction velocity. ---*/

          Tau_Wall = (1.0/Density_Wall)*pow(Y_Plus*Lam_Visc_Wall/WallDistMod,2.0);

          /*--- Difference between the old and new Tau. Update old value. ---*/

          diff = fabs(Tau_Wall-Tau_Wall_Old);
          Tau_Wall_Old += 0.25*(Tau_Wall-Tau_Wall_Old);

          counter++;
          if (counter > max_iter) {
            cout << "WARNING: Tau_Wall evaluation has not converged in solver_direct_turbulent" << endl;
            break;
          }

        }

        /*--- Now compute the Eddy viscosity at the first point off of the wall ---*/

        Lam_Visc_Normal = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint_Neighbor);
        Density_Normal = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint_Neighbor);
        Kin_Visc_Normal = Lam_Visc_Normal/Density_Normal;

        dypw_dyp = 2.0*Y_Plus_White*(kappa*sqrt(Gam)/Q)*sqrt(1.0 - pow(2.0*Gam*U_Plus - Beta,2.0)/(Q*Q));
        Eddy_Visc = Lam_Visc_Wall*(1.0 + dypw_dyp - kappa*exp(-1.0*kappa*B)*
                                             (1.0 + kappa*U_Plus
                                              + kappa*kappa*U_Plus*U_Plus/2.0)
                                             - Lam_Visc_Normal/Lam_Visc_Wall);

        /*--- Eddy viscosity should be always a positive number ---*/

        Eddy_Visc = max(0.0, Eddy_Visc);

        /*--- Solve for the new value of nu_tilde given the eddy viscosity and using a Newton method ---*/

        nu_til_old = 0.0; nu_til = 0.0; cv1_3 = 7.1*7.1*7.1;
        nu_til_old = nodes->GetSolution(iPoint,0);
        counter = 0; diff = 1.0;

        while (diff > tol) {

          func = nu_til_old*nu_til_old*nu_til_old*nu_til_old - (Eddy_Visc/Density_Normal)*(nu_til_old*nu_til_old*nu_til_old + Kin_Visc_Normal*Kin_Visc_Normal*Kin_Visc_Normal*cv1_3);
          func_prim = 4.0 * nu_til_old*nu_til_old*nu_til_old - 3.0*(Eddy_Visc/Density_Normal)*(nu_til_old*nu_til_old);
          nu_til = nu_til_old - func/func_prim;

          diff = fabs(nu_til-nu_til_old);
          nu_til_old = nu_til;

          counter++;
          if (counter > max_iter) {
            cout << "WARNING: Nu_tilde evaluation has not converged." << endl;
            break;
          }

        }

        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = nu_til;

        nodes->SetSolution_Old(iPoint_Neighbor,Solution);
        LinSysRes.SetBlock_Zero(iPoint_Neighbor);

        /*--- includes 1 in the diagonal ---*/

        Jacobian.DeleteValsRowi(iPoint_Neighbor);

      }

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
  const su2double *coord_i = nullptr, *coord_j = nullptr, *const *primVarGrad = nullptr, *vorticity = nullptr;
  su2double delta[3] = {0.0}, ratioOmega[3] = {0.0}, vortexTiltingMeasure = 0.0;

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){

    coord_i                 = geometry->node[iPoint]->GetCoord();
    nNeigh                  = geometry->node[iPoint]->GetnPoint();
    wallDistance            = geometry->node[iPoint]->GetWall_Distance();
    primVarGrad             = solver[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);
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

        maxDelta = geometry->node[iPoint]->GetMaxLength();
        distDES = constDES * maxDelta;
        lengthScale = min(distDES,wallDistance);

        break;

      case SA_DDES:
        /*--- A New Version of Detached-eddy Simulation, Resistant to Ambiguous Grid Densities.
         Spalart et al.
         Theoretical and Computational Fluid Dynamics - 2006
         ---*/

        maxDelta = geometry->node[iPoint]->GetMaxLength();

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
            jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
            coord_j = geometry->node[jPoint]->GetCoord();
            for ( iDim = 0; iDim < nDim; iDim++){
              deltaAux = abs(coord_j[iDim] - coord_i[iDim]);
              delta[iDim] = max(delta[iDim], deltaAux);
            }
            deltaDDES = geometry->node[iPoint]->GetMaxLength();
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
          jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
          coord_j = geometry->node[jPoint]->GetCoord();
          for (iDim = 0; iDim < nDim; iDim++){
            delta[iDim] = fabs(coord_j[iDim] - coord_i[iDim]);
          }
          deltaDDES = geometry->node[iPoint]->GetMaxLength();
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
}

void CTurbSASolver::SetInletAtVertex(su2double *val_inlet,
                                    unsigned short iMarker,
                                    unsigned long iVertex) {

  Inlet_TurbVars[iMarker][iVertex][0] = val_inlet[nDim+2+nDim];

}

su2double CTurbSASolver::GetInletAtVertex(su2double *val_inlet,
                                          unsigned long val_inlet_point,
                                          unsigned short val_kind_marker,
                                          string val_marker,
                                          CGeometry *geometry,
                                          CConfig *config) const {

  /*--- Local variables ---*/

  unsigned short iMarker, iDim;
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
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            Area = sqrt(Area);

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

void CTurbSASolver::SetUniformInlet(CConfig* config, unsigned short iMarker) {

  for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
    Inlet_TurbVars[iMarker][iVertex][0] = nu_tilde_Inf;
  }

}
