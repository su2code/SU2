/*!
 * \file CTurbSSTSolver.cpp
 * \brief Main subrotuines of CTurbSSTSolver class
 * \author F. Palacios, A. Bueno
 * \version 7.0.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

CTurbSSTSolver::CTurbSSTSolver(void) : CTurbSolver() {

  /*--- Array initialization ---*/
  constants = NULL;
  Inlet_TurbVars = NULL;

}

CTurbSSTSolver::CTurbSSTSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
    : CTurbSolver(geometry, config) {
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  ifstream restart_file;
  string text_line;

  bool multizone = config->GetMultizone_Problem();

  /*--- Array initialization ---*/

  constants = NULL;

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

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

    Residual = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]  = 0.0;
    Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
    Residual_i = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]  = 0.0;
    Residual_j = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]  = 0.0;
    Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max = new unsigned long[nVar];
    for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
    Point_Max_Coord = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
    }

    /*--- Define some auxiliary vector related with the solution ---*/

    Solution = new su2double[nVar];
    Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];

    /*--- Define some auxiliary vector related with the geometry ---*/

    Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];

    /*--- Define some auxiliary vector related with the flow solution ---*/

    FlowPrimVar_i = new su2double [nDim+9]; FlowPrimVar_j = new su2double [nDim+9];

    /*--- Jacobians and vector structures for implicit computations ---*/

    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (SST model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (multizone){
      Residual_BGS      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]  = 0.0;
      Residual_Max_BGS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 0.0;

      /*--- Define some structures for locating max residuals ---*/

      Point_Max_BGS       = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar]  = 0;
      Point_Max_Coord_BGS = new su2double*[nVar];
      for (iVar = 0; iVar < nVar; iVar++) {
        Point_Max_Coord_BGS[iVar] = new su2double[nDim];
        for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
      }
    }

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

  /*--- Initialize value for model constants ---*/
  constants = new su2double[10];
  constants[0] = 0.85;   //sigma_k1
  constants[1] = 1.0;    //sigma_k2
  constants[2] = 0.5;    //sigma_om1
  constants[3] = 0.856;  //sigma_om2
  constants[4] = 0.075;  //beta_1
  constants[5] = 0.0828; //beta_2
  constants[6] = 0.09;   //betaStar
  constants[7] = 0.31;   //a1
  // constants[8] = constants[4]/constants[6] - constants[2]*0.41*0.41/sqrt(constants[6]);  //alfa_1
  // constants[9] = constants[5]/constants[6] - constants[3]*0.41*0.41/sqrt(constants[6]);  //alfa_2
  constants[8] = 5.0/9.0;  //alfa_1
  constants[9] = 0.44;  //alfa_2

  /*--- Initialize lower and upper limits---*/
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];

  lowerlimit[0] = 1.0e-10;
  upperlimit[0] = 1.0e10;

  lowerlimit[1] = 1.0e-4;
  upperlimit[1] = 1.0e15;

  /*--- Far-field flow state quantities and initialization. ---*/
  su2double rhoInf, *VelInf, muLamInf, Intensity, viscRatio, muT_Inf;

  rhoInf    = config->GetDensity_FreeStreamND();
  VelInf    = config->GetVelocity_FreeStreamND();
  muLamInf  = config->GetViscosity_FreeStreamND();
  Intensity = config->GetTurbulenceIntensity_FreeStream();
  viscRatio = config->GetTurb2LamViscRatio_FreeStream();

  su2double VelMag = 0;
  for (iDim = 0; iDim < nDim; iDim++)
    VelMag += VelInf[iDim]*VelInf[iDim];
  VelMag = sqrt(VelMag);

  kine_Inf  = 3.0/2.0*(VelMag*VelMag*Intensity*Intensity);
  omega_Inf = rhoInf*kine_Inf/(muLamInf*viscRatio);

  /*--- Eddy viscosity, initialized without stress limiter at the infinity ---*/
  muT_Inf = rhoInf*kine_Inf/omega_Inf;

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CTurbSSTVariable(kine_Inf, omega_Inf, muT_Inf, nPoint, nDim, nVar, constants, config);
  SetBaseClassPointerToNodes();

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION_EDDY);
  CompleteComms(geometry, config, SOLUTION_EDDY);

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  unsigned long iMarker;

  SlidingState       = new su2double*** [nMarker];
  SlidingStateNodes  = new int*         [nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++){

    SlidingState[iMarker]      = NULL;
    SlidingStateNodes[iMarker] = NULL;

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){

      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)];
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)];

      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++){
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1];

        SlidingStateNodes[iMarker][iPoint] = 0;
        for (iVar = 0; iVar < nPrimVar+1; iVar++)
          SlidingState[iMarker][iPoint][iVar] = NULL;
      }

    }
  }

  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
    due to arbitrary number of turbulence variables ---*/

      Inlet_TurbVars = new su2double**[nMarker];
      for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
        Inlet_TurbVars[iMarker] = new su2double*[nVertex[iMarker]];
        for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
          Inlet_TurbVars[iMarker][iVertex] = new su2double[nVar];
          Inlet_TurbVars[iMarker][iVertex][0] = kine_Inf;
          Inlet_TurbVars[iMarker][iVertex][1] = omega_Inf;
        }
      }

  /*--- The turbulence models are always solved implicitly, so set the
  implicit flag in case we have periodic BCs. ---*/

  SetImplicitPeriodic(true);

  /* Store the initial CFL number for all grid points. */

  const su2double CFL = config->GetCFL(MGLevel);
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "K-W SST";

}

CTurbSSTSolver::~CTurbSSTSolver(void) {

  if (constants != NULL) delete [] constants;

  unsigned long iMarker, iVertex;
  unsigned short iVar;

  if ( SlidingState != NULL ) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if ( SlidingState[iMarker] != NULL ) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          if ( SlidingState[iMarker][iVertex] != NULL ){
            for (iVar = 0; iVar < nPrimVar+1; iVar++)
              delete [] SlidingState[iMarker][iVertex][iVar];
            delete [] SlidingState[iMarker][iVertex];
          }
        delete [] SlidingState[iMarker];
      }
    }
    delete [] SlidingState;
  }

  if ( SlidingStateNodes != NULL ){
    for (iMarker = 0; iMarker < nMarker; iMarker++){
        if (SlidingStateNodes[iMarker] != NULL)
            delete [] SlidingStateNodes[iMarker];
    }
    delete [] SlidingStateNodes;
  }

}

void CTurbSSTSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint;

  bool limiter_turb = (config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (config->GetInnerIter() <= config->GetLimiterIter());


  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Initialize the residual vector ---*/

    LinSysRes.SetBlock_Zero(iPoint);

  }

  /*--- Initialize the Jacobian matrices ---*/

  Jacobian.SetValZero();

  /*--- Upwind second order reconstruction and gradients ---*/

  if (config->GetReconstructionGradientRequired()) {
    if (config->GetKind_Gradient_Method_Recon() == GREEN_GAUSS)
      SetSolution_Gradient_GG(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == LEAST_SQUARES)
      SetSolution_Gradient_LS(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES)
      SetSolution_Gradient_LS(geometry, config, true);
  }
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

  if (limiter_turb) SetSolution_Limiter(geometry, config);

}

void CTurbSSTSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  su2double rho = 0.0, mu = 0.0, dist, omega, kine, F2, muT, zeta;
  su2double a1 = constants[7];
  unsigned long iPoint;

  /*--- Compute mean flow and turbulence gradients ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetSolution_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetSolution_Gradient_LS(geometry, config);
  }

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Compute blending functions and cross diffusion ---*/

    rho = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    mu  = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);

    dist = geometry->node[iPoint]->GetWall_Distance();

    // su2double *Vorticity = solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint);
    // su2double VorticityMag = sqrt(Vorticity[0]*Vorticity[0] +
    //                               Vorticity[1]*Vorticity[1] +
    //                               Vorticity[2]*Vorticity[2]);
    // su2double StrainMag = solver_container[FLOW_SOL]->GetNodes()->GetStrainMag(iPoint);
    su2double StrainMag = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++) {
      StrainMag += pow(solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint,iDim+1,iDim), 2.0);
    }
    StrainMag += 2.0*pow(0.5*(solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint,1,1) 
                            + solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint,2,0)), 2);

    if (nDim == 3) {
      StrainMag += 2.0*pow(0.5*(solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint,1,2) 
                              + solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint,3,0)), 2);
      StrainMag += 2.0*pow(0.5*(solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint,2,2) 
                              + solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint,3,1)), 2);
    }

    StrainMag = sqrt(2.0*StrainMag);

    nodes->SetBlendingFunc(iPoint,mu, dist, rho);

    F2 = nodes->GetF2blending(iPoint);

    /*--- Compute the eddy viscosity ---*/

    kine  = nodes->GetSolution(iPoint,0);
    omega = nodes->GetSolution(iPoint,1);
    // zeta  = min(1.0/omega, a1/(VorticityMag*F2));
    zeta  = min(1.0/omega, a1/(StrainMag*F2));
    muT   = max(rho*kine*zeta,0.0);
    nodes->SetmuT(iPoint,muT);

  }

}

void CTurbSSTSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) {

  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint), NULL);

    /*--- Gradient of the primitive and conservative variables ---*/

    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint), NULL);

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/

    numerics->SetTurbVar(nodes->GetSolution(iPoint), NULL);
    numerics->SetTurbVarGradient(nodes->GetGradient(iPoint), NULL);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->node[iPoint]->GetVolume());

    /*--- Set distance to the surface ---*/

    numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);

    /*--- Menter's first blending function ---*/

    numerics->SetF1blending(nodes->GetF1blending(iPoint),0.0);

    /*--- Menter's second blending function ---*/

    numerics->SetF2blending(nodes->GetF2blending(iPoint),0.0);

    /*--- Set vorticity and strain rate magnitude ---*/

    numerics->SetVorticity(solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint), NULL);

    numerics->SetStrainMag(solver_container[FLOW_SOL]->GetNodes()->GetStrainMag(iPoint), 0.0);

    /*--- Cross diffusion ---*/

    numerics->SetCrossDiff(nodes->GetCrossDiff(iPoint),0.0);

    /*--- Compute the source term ---*/

    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);

    if(config->GetError_Estimate()) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) nodes->SetSource(iPoint, iVar, Residual[iVar]/geometry->node[iPoint]->GetVolume());
    }

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

  }

}

void CTurbSSTSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {

}

void CTurbSSTSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar;
  su2double distance, density = 0.0, laminar_viscosity = 0.0, beta_1;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- distance to closest neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      distance = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        distance += (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim))*
        (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      distance = sqrt(distance);

      /*--- Set wall values ---*/

      density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(jPoint);
      laminar_viscosity = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(jPoint);

      beta_1 = constants[4];

      Solution[0] = 0.0;
      Solution[1] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);

      /*--- Set the solution values and zero the residual ---*/
      nodes->SetSolution_Old(iPoint,Solution);
      nodes->SetSolution(iPoint,Solution);
      LinSysRes.SetBlock_Zero(iPoint);

      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }

    }
  }

}

void CTurbSSTSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                        unsigned short val_marker) {

  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar;
  su2double distance, density = 0.0, laminar_viscosity = 0.0, beta_1;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- distance to closest neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      distance = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        distance += (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim))*
        (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      distance = sqrt(distance);

      /*--- Set wall values ---*/

      density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(jPoint);
      laminar_viscosity = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(jPoint);

      beta_1 = constants[4];

      Solution[0] = 0.0;
      Solution[1] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);

      /*--- Set the solution values and zero the residual ---*/
      nodes->SetSolution_Old(iPoint,Solution);
      nodes->SetSolution(iPoint,Solution);
      LinSysRes.SetBlock_Zero(iPoint);

      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }

    }
  }

}

void CTurbSSTSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;
  su2double *Normal, *V_infty, *V_domain;
  unsigned short iVar, iDim;

  Normal = new su2double[nDim];

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Allocate the value at the infinity ---*/

      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Set turbulent variable at the wall, and at infinity ---*/

      for (iVar = 0; iVar < nVar; iVar++)
      Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);

      Solution_j[0] = kine_Inf;
      Solution_j[1] = omega_Inf;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set Normal (it is necessary to change the sign) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
      Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Grid Movement ---*/

      if (dynamic_grid)
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute residuals and Jacobians ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Add residuals and Jacobians ---*/

      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  delete [] Normal;

}

void CTurbSSTSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                              unsigned short val_marker) {

  unsigned short iVar, iDim;
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

      /*--- Set the turbulent variable states. Use free-stream SST
       values for the turbulent state at the inflow. ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);

      /*--- Load the inlet turbulence variables (uniform by default). ---*/

      Solution_j[0] = Inlet_TurbVars[val_marker][iVertex][0];
      Solution_j[1] = Inlet_TurbVars[val_marker][iVertex][1];

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the solver class ---*/

      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/

      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

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
      //     visc_numerics->SetTurbVar(Solution_i, Solution_j);
      //      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      //
      //      /*--- Menter's first blending function ---*/
      //
      //      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());
      //
      //      /*--- Compute residual, and Jacobians ---*/
      //
      //      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      //
      //      /*--- Subtract residual, and update Jacobians ---*/
      //
      //      LinSysRes.SubtractBlock(iPoint, Residual);
      //      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

    }

  }

  /*--- Free locally allocated memory ---*/

  delete [] Normal;

}

void CTurbSSTSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

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

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/

      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

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
//      /*--- Menter's first blending function ---*/
//
//      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());
//
//      /*--- Compute residual, and Jacobians ---*/
//
//      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//
//      LinSysRes.SubtractBlock(iPoint, Residual);
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}


void CTurbSSTSolver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                              unsigned short val_marker) {

  unsigned short iVar, iSpan, iDim;
  unsigned long  oldVertex, iPoint, Point_Normal, iVertex;
  su2double *V_inlet, *V_domain, *Normal;
  su2double extAverageKine, extAverageOmega;
  unsigned short nSpanWiseSections = config->GetnSpanWiseSections();

  Normal = new su2double[nDim];

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iSpan= 0; iSpan < nSpanWiseSections ; iSpan++){
    extAverageKine = solver_container[FLOW_SOL]->GetExtAverageKine(val_marker, iSpan);
    extAverageOmega = solver_container[FLOW_SOL]->GetExtAverageOmega(val_marker, iSpan);


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

      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);

      Solution_j[0]= extAverageKine;
      Solution_j[1]= extAverageOmega;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
            geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(iPoint));

      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}

void CTurbSSTSolver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                              unsigned short val_marker) {

  unsigned short iVar, iSpan, iDim;
  unsigned long  oldVertex, iPoint, Point_Normal, iVertex;
  su2double *V_inlet, *V_domain, *Normal;
  unsigned short nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- Quantities for computing the  kine and omega to impose at the inlet boundary. ---*/
  su2double rho, pressure, *Vel, VelMag, muLam, Intensity, viscRatio, kine_b, omega_b, kine;
  CFluidModel *FluidModel;

  FluidModel = solver_container[FLOW_SOL]->GetFluidModel();
  Intensity = config->GetTurbulenceIntensity_FreeStream();
  viscRatio = config->GetTurb2LamViscRatio_FreeStream();

  Normal = new su2double[nDim];
  Vel = new su2double[nDim];

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);


  for (iSpan= 0; iSpan < nSpanWiseSections ; iSpan++){

    /*--- Compute the inflow kine and omega using the span wise averge quntities---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Vel[iDim] = solver_container[FLOW_SOL]->GetAverageTurboVelocity(val_marker, iSpan)[iDim];

    rho       = solver_container[FLOW_SOL]->GetAverageDensity(val_marker, iSpan);
    pressure  = solver_container[FLOW_SOL]->GetAveragePressure(val_marker, iSpan);
    kine      = solver_container[FLOW_SOL]->GetAverageKine(val_marker, iSpan);

    FluidModel->SetTDState_Prho(pressure, rho);
    muLam = FluidModel->GetLaminarViscosity();

    VelMag = 0;
    for (iDim = 0; iDim < nDim; iDim++)
      VelMag += Vel[iDim]*Vel[iDim];
    VelMag = sqrt(VelMag);

    kine_b  = 3.0/2.0*(VelMag*VelMag*Intensity*Intensity);
    omega_b = rho*kine/(muLam*viscRatio);

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

      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);

      /*--- Set the turbulent variable states. Use average span-wise values
             values for the turbulent state at the inflow. ---*/

      Solution_j[0]= kine_b;
      Solution_j[1]= omega_b;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
            geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(iPoint));

      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  delete[] Vel;

}


void CTurbSSTSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
    CNumerics *visc_numerics, CConfig *config){

  unsigned long iVertex, jVertex, iPoint, Point_Normal = 0;
  unsigned short iDim, iVar, iMarker;

  unsigned short nPrimVar = solver_container[FLOW_SOL]->GetnPrimVar();
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  su2double *tmp_residual = new su2double[nVar];

  unsigned long nDonorVertex;
  su2double weight;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

        if (geometry->node[iPoint]->GetDomain()) {

          nDonorVertex = GetnSlidingStates(iMarker, iVertex);

          /*--- Initialize Residual, this will serve to accumulate the average ---*/

          for (iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] = 0.0;

          /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/

          for (jVertex = 0; jVertex < nDonorVertex; jVertex++){

            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

            for (iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_i[iVar] = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint,iVar);
              PrimVar_j[iVar] = solver_container[FLOW_SOL]->GetSlidingState(iMarker, iVertex, iVar, jVertex);
            }

            /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/

            weight = solver_container[FLOW_SOL]->GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);

            /*--- Set primitive variables ---*/

            conv_numerics->SetPrimitive( PrimVar_i, PrimVar_j );

            /*--- Set the turbulent variable states ---*/
            Solution_i[0] = nodes->GetSolution(iPoint,0);
            Solution_i[1] = nodes->GetSolution(iPoint,1);

            Solution_j[0] = GetSlidingState(iMarker, iVertex, 0, jVertex);
            Solution_j[1] = GetSlidingState(iMarker, iVertex, 1, jVertex);

            conv_numerics->SetTurbVar(Solution_i, Solution_j);

            /*--- Set the normal vector ---*/

            conv_numerics->SetNormal(Normal);

            if (dynamic_grid)
              conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

            conv_numerics->ComputeResidual(tmp_residual, Jacobian_i, Jacobian_j, config);

            /*--- Accumulate the residuals to compute the average ---*/

            for (iVar = 0; iVar < nVar; iVar++)
              Residual[iVar] += weight*tmp_residual[iVar];
          }

          /*--- Add Residuals and Jacobians ---*/

          LinSysRes.AddBlock(iPoint, Residual);

          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

          /*--- Set the normal vector and the coordinates ---*/

          visc_numerics->SetNormal(Normal);
          visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

          /*--- Primitive variables, and gradient ---*/

          visc_numerics->SetPrimitive(PrimVar_i, PrimVar_j);
          //          visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

          /*--- Turbulent variables and its gradients  ---*/

          visc_numerics->SetTurbVar(Solution_i, Solution_j);
          visc_numerics->SetTurbVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

          /*--- Compute and update residual ---*/

          visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

          LinSysRes.SubtractBlock(iPoint, Residual);

          /*--- Jacobian contribution for implicit integration ---*/

          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

        }
      }
    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] tmp_residual;
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;

}

void CTurbSSTSolver::SetInletAtVertex(su2double *val_inlet,
                                     unsigned short iMarker,
                                     unsigned long iVertex) {

  Inlet_TurbVars[iMarker][iVertex][0] = val_inlet[nDim+2+nDim];
  Inlet_TurbVars[iMarker][iVertex][1] = val_inlet[nDim+2+nDim+1];

}

su2double CTurbSSTSolver::GetInletAtVertex(su2double *val_inlet,
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
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            Area = sqrt(Area);

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

void CTurbSSTSolver::SetUniformInlet(CConfig* config, unsigned short iMarker) {

  for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
    Inlet_TurbVars[iMarker][iVertex][0] = kine_Inf;
    Inlet_TurbVars[iMarker][iVertex][1] = omega_Inf;
  }

}

void CTurbSSTSolver::SetTurbGradient_L2Proj2(CGeometry *geometry, CConfig *config, CSolver *solver_flow) {

  unsigned long iPoint, nPoint = geometry->GetnPoint(), iElem, nElem = geometry->GetnElem();
  unsigned short iVar, iFlux;
  unsigned short nVarMetr = 2, nFluxMetr = 2;  //--- TODO: adjust size of grad vector later for goal vs. feature
  su2double density, velocity[2];
  su2double vnx[3], vny[3];
  su2double graTri[2], graTriVisc[2], graTriSource[2];
  su2double Crd[3][2], Sens[3][nVarMetr][nFluxMetr], SensVisc[3][nVarMetr][nFluxMetr], SensSource[3][nVarMetr];
  su2double laminar_viscosity, eddy_viscosity;
  su2double k, omega, dk[2], domega[2];
  bool dummy_bool;

  const su2double sigmak1 = constants[0], sigmaom1 = constants[2],
                  sigmak2 = constants[1], sigmaom2 = constants[3];

  //--- note: currently only implemented for Tri

  for (iElem=0; iElem<nElem; ++iElem) {
    for (unsigned short iNode=0; iNode<3; ++iNode) {
      const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
      //--- store coordinates
      for (unsigned short iDim = 0; iDim<2; ++iDim) {
        Crd[iNode][iDim] = geometry->node[kNode]->GetCoord(iDim);
      }
      //--- store sensors (goal-oriented)
      density     = solver_flow->GetNodes()->GetDensity(kNode);
      velocity[0] = solver_flow->GetNodes()->GetVelocity(kNode, 0);
      velocity[1] = solver_flow->GetNodes()->GetVelocity(kNode, 1);

      laminar_viscosity = solver_flow->GetNodes()->GetLaminarViscosity(kNode);
      eddy_viscosity    = solver_flow->GetNodes()->GetEddyViscosity(kNode);

      k     = nodes->GetSolution(kNode, 0);
      dk[0] = nodes->GetGradient(kNode, 0,0);
      dk[1] = nodes->GetGradient(kNode, 0,1);

      omega     = nodes->GetSolution(kNode, 1);
      domega[0] = nodes->GetGradient(kNode, 1,0);
      domega[1] = nodes->GetGradient(kNode, 1,1);

      const su2double F1 = nodes->GetF1blending(kNode);
      const su2double sigmak = F1*sigmak1 + (1.0 - F1)*sigmak2;
      const su2double sigmaomega = F1*sigmaom1 + (1.0 - F1)*sigmaom2;

      Sens[iNode][0][0] = density*velocity[0]*k;
      Sens[iNode][0][1] = density*velocity[1]*k;

      SensVisc[iNode][0][0] = (laminar_viscosity+sigmak*eddy_viscosity)*dk[0];
      SensVisc[iNode][0][1] = (laminar_viscosity+sigmak*eddy_viscosity)*dk[1];

      Sens[iNode][1][0] = density*velocity[0]*omega;
      Sens[iNode][1][1] = density*velocity[1]*omega;

      SensVisc[iNode][1][0] = (laminar_viscosity+sigmaomega*eddy_viscosity)*domega[0];
      SensVisc[iNode][1][1] = (laminar_viscosity+sigmaomega*eddy_viscosity)*domega[1];

      SensSource[iNode][0] = nodes->GetSource(kNode)[0];
      SensSource[iNode][1] = nodes->GetSource(kNode)[1];
    }// iNode

    //--- inward edge's normals : edg[0]=P1P2, edg[1]=P2P3P00, edg[2]=P0P1
    vnx[0] = Crd[1][1]-Crd[2][1];
    vny[0] = Crd[2][0]-Crd[1][0];

    vnx[1] = Crd[2][1]-Crd[0][1];
    vny[1] = Crd[0][0]-Crd[2][0];

    vnx[2] = Crd[0][1]-Crd[1][1];
    vny[2] = Crd[1][0]-Crd[0][0];

    //--- check if inward normal
    for(unsigned short iNode = 0; iNode < 3; ++iNode) {
      su2double CrdAvg[2] = {0.0, 0.0};
      for(unsigned short jNode = 0; jNode < 3; ++jNode) {
        if(iNode != jNode) {
          CrdAvg[0] += Crd[jNode][0];
          CrdAvg[1] += Crd[jNode][1];
        }
      }// jNode
      CrdAvg[0] /= 2.;
      CrdAvg[1] /= 2.;
      su2double u[2] = {CrdAvg[0]-Crd[iNode][0],
                        CrdAvg[1]-Crd[iNode][1]};
      if((vnx[iNode]*u[0] + vny[iNode]*u[1]) > 0.) {
        vnx[iNode] *= -1.0;
        vny[iNode] *= -1.0;
      }
    }// iNode

    //--- loop over conservative variables
    for(iVar = 0; iVar < nVarMetr; iVar++){

      //--- loop over directions
      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){

        //--- gradient at the element ( graTri = 2*|T|*gradT ) 
        graTri[0] = Sens[0][iVar][iFlux]*vnx[0] + Sens[1][iVar][iFlux]*vnx[1] + Sens[2][iVar][iFlux]*vnx[2];
        graTri[1] = Sens[0][iVar][iFlux]*vny[0] + Sens[1][iVar][iFlux]*vny[1] + Sens[2][iVar][iFlux]*vny[2];

        graTriVisc[0] = SensVisc[0][iVar][iFlux]*vnx[0] + SensVisc[1][iVar][iFlux]*vnx[1] + SensVisc[2][iVar][iFlux]*vnx[2];
        graTriVisc[1] = SensVisc[0][iVar][iFlux]*vny[0] + SensVisc[1][iVar][iFlux]*vny[1] + SensVisc[2][iVar][iFlux]*vny[2];
    
        //--- assembling
        const unsigned short i = iFlux*nVarMetr*nDim + iVar*nDim;
        for (unsigned short iNode=0; iNode<3; ++iNode) {
          const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
          const su2double Area = geometry->node[kNode]->GetVolume();
          const su2double rap = 1./(Area*6.);
          nodes->AddAnisoGrad(kNode, i+0, graTri[0] * rap);
          nodes->AddAnisoGrad(kNode, i+1, graTri[1] * rap);

          nodes->AddAnisoViscGrad(kNode, i+0, graTriVisc[0] * rap);
          nodes->AddAnisoViscGrad(kNode, i+1, graTriVisc[1] * rap);
        }// iNode
      }// iFlux

      graTriSource[0] = SensSource[0][iVar]*vnx[0] + SensSource[1][iVar]*vnx[1] + SensSource[2][iVar]*vnx[2];
      graTriSource[1] = SensSource[0][iVar]*vny[0] + SensSource[1][iVar]*vny[1] + SensSource[2][iVar]*vny[2];
      //--- assembling
      const unsigned short i = iVar*nDim;
      for (unsigned short iNode=0; iNode<3; ++iNode) {
        const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
        const su2double Area = geometry->node[kNode]->GetVolume();
        const su2double rap = 1./(Area*6.);
        nodes->AddAnisoSourceGrad(kNode, i+0, graTriSource[0] * rap);
        nodes->AddAnisoSourceGrad(kNode, i+1, graTriSource[1] * rap);
      }// iNode
    }// iVar
  }// iElem

  /*--- Communicate the gradient values via MPI. ---*/
  
  InitiateComms(geometry, config, ANISO_GRADIENT);
  CompleteComms(geometry, config, ANISO_GRADIENT);

  InitiateComms(geometry, config, ANISO_GRADIENT_VISC);
  CompleteComms(geometry, config, ANISO_GRADIENT_VISC);

  InitiateComms(geometry, config, ANISO_GRADIENT_SOURCE);
  CompleteComms(geometry, config, ANISO_GRADIENT_SOURCE);
}

void CTurbSSTSolver::SetHessian_L2Proj2(CGeometry *geometry, CConfig *config){

  unsigned long iPoint, nPoint = geometry->GetnPoint(), nPointDomain = geometry->GetnPointDomain(), iElem, nElem = geometry->GetnElem();
  unsigned short iVar, iFlux;
  unsigned short nVarMetr = 2, nFluxMetr = 2;  //--- TODO: adjust size of grad vector later for goal vs. feature
  unsigned short nMetr = 3;
  su2double vnx[3], vny[3];
  su2double hesTri[3], hesTriVisc[3], hesTriSource[3];
  su2double Crd[3][2], Grad[3][2][nVarMetr][nFluxMetr], GradVisc[3][2][nVarMetr][nFluxMetr], GradSource[3][2][nVarMetr];

  su2double **A      = new su2double*[nDim],
            **EigVec = new su2double*[nDim], 
            *EigVal  = new su2double[nDim];

  for(unsigned short iDim = 0; iDim < nDim; ++iDim){
    A[iDim]      = new su2double[nDim];
    EigVec[iDim] = new su2double[nDim];
  }

  //--- note: currently only implemented for Tri

  for (iElem=0; iElem<nElem; ++iElem) {
    for (unsigned short iNode=0; iNode<3; ++iNode) {
      const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
      //--- store coordinates
      for (unsigned short iDim = 0; iDim<2; ++iDim) {
        Crd[iNode][iDim] = geometry->node[kNode]->GetCoord(iDim);
      }
      //--- store gradient
      for(iVar = 0; iVar < nVarMetr; iVar++){
        for(iFlux = 0; iFlux < nFluxMetr; iFlux++){
          const unsigned short i = iFlux*nVarMetr*nDim + iVar*nDim;
          Grad[iNode][0][iVar][iFlux] = nodes->GetAnisoGrad(kNode, i+0);
          Grad[iNode][1][iVar][iFlux] = nodes->GetAnisoGrad(kNode, i+1);

          GradVisc[iNode][0][iVar][iFlux] = nodes->GetAnisoViscGrad(kNode, i+0);
          GradVisc[iNode][1][iVar][iFlux] = nodes->GetAnisoViscGrad(kNode, i+1);
        }// iFlux

        const unsigned short i = iVar*nDim;
        GradSource[iNode][0][iVar] = nodes->GetAnisoSourceGrad(kNode, i+0);
        GradSource[iNode][1][iVar] = nodes->GetAnisoSourceGrad(kNode, i+1);
      }// iVar
    }// iNode

    //--- inward edge's normals : edg[0]=P1P2, edg[1]=P2P0, edg[2]=P0P1
    vnx[0] = Crd[1][1]-Crd[2][1];
    vny[0] = Crd[2][0]-Crd[1][0];

    vnx[1] = Crd[2][1]-Crd[0][1];
    vny[1] = Crd[0][0]-Crd[2][0];

    vnx[2] = Crd[0][1]-Crd[1][1];
    vny[2] = Crd[1][0]-Crd[0][0];

    //--- check if inward normal
    for(unsigned short iNode = 0; iNode < 3; ++iNode) {
      su2double CrdAvg[2] = {0.0, 0.0};
      for(unsigned short jNode = 0; jNode < 3; ++jNode) {
        if(iNode != jNode) {
          CrdAvg[0] += Crd[jNode][0];
          CrdAvg[1] += Crd[jNode][1];
        }
      }
      CrdAvg[0] /= 2.;
      CrdAvg[1] /= 2.;
      su2double u[2] = {CrdAvg[0]-Crd[iNode][0],
                        CrdAvg[1]-Crd[iNode][1]};
      if((vnx[iNode]*u[0] + vny[iNode]*u[1]) > 0.) {
        vnx[iNode] *= -1.0;
        vny[iNode] *= -1.0;
      }// jNode
    }// iNode

    //--- loop over conservative variables
    for(iVar = 0; iVar < nVarMetr; iVar++){

      //--- loop over directions
      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){

        //--- hessian at the element ( hesTri = 2*|T|*hessienT ) 
        hesTri[0] =         Grad[0][0][iVar][iFlux]*vnx[0] 
                          + Grad[1][0][iVar][iFlux]*vnx[1] 
                          + Grad[2][0][iVar][iFlux]*vnx[2];
        hesTri[1] = 0.5 * ( Grad[0][0][iVar][iFlux]*vny[0] 
                          + Grad[1][0][iVar][iFlux]*vny[1] 
                          + Grad[2][0][iVar][iFlux]*vny[2]
                          + Grad[0][1][iVar][iFlux]*vnx[0] 
                          + Grad[1][1][iVar][iFlux]*vnx[1] 
                          + Grad[2][1][iVar][iFlux]*vnx[2] );
        hesTri[2] =         Grad[0][1][iVar][iFlux]*vny[0] 
                          + Grad[1][1][iVar][iFlux]*vny[1] 
                          + Grad[2][1][iVar][iFlux]*vny[2];

        hesTriVisc[0] =         GradVisc[0][0][iVar][iFlux]*vnx[0] 
                              + GradVisc[1][0][iVar][iFlux]*vnx[1] 
                              + GradVisc[2][0][iVar][iFlux]*vnx[2];
        hesTriVisc[1] = 0.5 * ( GradVisc[0][0][iVar][iFlux]*vny[0] 
                              + GradVisc[1][0][iVar][iFlux]*vny[1] 
                              + GradVisc[2][0][iVar][iFlux]*vny[2]
                              + GradVisc[0][1][iVar][iFlux]*vnx[0] 
                              + GradVisc[1][1][iVar][iFlux]*vnx[1] 
                              + GradVisc[2][1][iVar][iFlux]*vnx[2] );
        hesTriVisc[2] =         GradVisc[0][1][iVar][iFlux]*vny[0] 
                              + GradVisc[1][1][iVar][iFlux]*vny[1] 
                              + GradVisc[2][1][iVar][iFlux]*vny[2];
        
        //--- assembling
        const unsigned short i = iFlux*nVarMetr*nMetr + iVar*nMetr;
        for (unsigned short iNode=0; iNode<3; ++iNode) {
          const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
          const su2double Area = geometry->node[kNode]->GetVolume();
          const su2double rap = 1./(Area*6.);
          nodes->AddAnisoHess(kNode, i+0, hesTri[0] * rap);
          nodes->AddAnisoHess(kNode, i+1, hesTri[1] * rap);
          nodes->AddAnisoHess(kNode, i+2, hesTri[2] * rap);

          nodes->AddAnisoViscHess(kNode, i+0, hesTriVisc[0] * rap);
          nodes->AddAnisoViscHess(kNode, i+1, hesTriVisc[1] * rap);
          nodes->AddAnisoViscHess(kNode, i+2, hesTriVisc[2] * rap);
        }// iNode
      }// iFlux

      hesTriSource[0] =         GradSource[0][0][iVar]*vnx[0] 
                              + GradSource[1][0][iVar]*vnx[1] 
                              + GradSource[2][0][iVar]*vnx[2];
      hesTriSource[1] = 0.5 * ( GradSource[0][0][iVar]*vny[0] 
                              + GradSource[1][0][iVar]*vny[1] 
                              + GradSource[2][0][iVar]*vny[2]
                              + GradSource[0][1][iVar]*vnx[0] 
                              + GradSource[1][1][iVar]*vnx[1] 
                              + GradSource[2][1][iVar]*vnx[2] );
      hesTriSource[2] =         GradSource[0][1][iVar]*vny[0] 
                              + GradSource[1][1][iVar]*vny[1] 
                              + GradSource[2][1][iVar]*vny[2];
    
      //--- assembling
      const unsigned short i = iVar*nMetr;
      for (unsigned short iNode=0; iNode<3; ++iNode) {
        const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
        const su2double Area = geometry->node[kNode]->GetVolume();
        const su2double rap = 1./(Area*6.);
        nodes->AddAnisoSourceHess(kNode, i+0, hesTriSource[0] * rap);
        nodes->AddAnisoSourceHess(kNode, i+1, hesTriSource[1] * rap);
        nodes->AddAnisoSourceHess(kNode, i+2, hesTriSource[2] * rap);
      }// iNode
    }// iVar
  }// iElem

  /*--- Communicate the Hessian values via MPI. ---*/
  
  InitiateComms(geometry, config, ANISO_HESSIAN);
  CompleteComms(geometry, config, ANISO_HESSIAN);
  InitiateComms(geometry, config, ANISO_HESSIAN_VISC);
  CompleteComms(geometry, config, ANISO_HESSIAN_VISC);
  InitiateComms(geometry, config, ANISO_HESSIAN_SOURCE);
  CompleteComms(geometry, config, ANISO_HESSIAN_SOURCE);

  CorrectBoundAnisoHess(geometry, config);
  CorrectBoundAnisoSourceHess(geometry, config);

  //--- Make positive definite matrix
  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    for(iVar = 0; iVar < nVarMetr; iVar++){
      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){
        const unsigned short i = iFlux*nVarMetr*nMetr + iVar*nMetr;

        const su2double a = nodes->GetAnisoHess(iPoint, i+0);
        const su2double b = nodes->GetAnisoHess(iPoint, i+1);
        const su2double c = nodes->GetAnisoHess(iPoint, i+2);
        
        A[0][0] = a; A[0][1] = b;
        A[1][0] = b; A[1][1] = c;

        CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

        for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = abs(EigVal[iDim]);

        CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

        nodes->SetAnisoHess(iPoint, i+0, A[0][0]);
        nodes->SetAnisoHess(iPoint, i+1, A[0][1]);
        nodes->SetAnisoHess(iPoint, i+2, A[1][1]);
      }// iFlux

      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){
        const unsigned short i = iFlux*nVarMetr*nMetr + iVar*nMetr;

        const su2double a = nodes->GetAnisoViscHess(iPoint, i+0);
        const su2double b = nodes->GetAnisoViscHess(iPoint, i+1);
        const su2double c = nodes->GetAnisoViscHess(iPoint, i+2);
        
        A[0][0] = a; A[0][1] = b;
        A[1][0] = b; A[1][1] = c;

        CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

        for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = abs(EigVal[iDim]);

        CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

        nodes->SetAnisoViscHess(iPoint, i+0, A[0][0]);
        nodes->SetAnisoViscHess(iPoint, i+1, A[0][1]);
        nodes->SetAnisoViscHess(iPoint, i+2, A[1][1]);
      }// iFlux

      const unsigned short i = iVar*nMetr;

      const su2double a = nodes->GetAnisoSourceHess(iPoint, i+0);
      const su2double b = nodes->GetAnisoSourceHess(iPoint, i+1);
      const su2double c = nodes->GetAnisoSourceHess(iPoint, i+2);

      A[0][0] = a; A[0][1] = b;
      A[1][0] = b; A[1][1] = c;

      CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

      for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = abs(EigVal[iDim]);

      CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

      nodes->SetAnisoSourceHess(iPoint, i+0, A[0][0]);
      nodes->SetAnisoSourceHess(iPoint, i+1, A[0][1]);
      nodes->SetAnisoSourceHess(iPoint, i+2, A[1][1]);
    }// iVar
  }// iPoint
}

void CTurbSSTSolver::SetTurbGradient_L2Proj3(CGeometry *geometry, CConfig *config, CSolver *solver_flow) {

  unsigned long iPoint, nPoint = geometry->GetnPoint(), iElem, nElem = geometry->GetnElem();
  unsigned short iVar, iDim, iFlux;
  unsigned short nVarMetr = 2, nFluxMetr = 3;
  su2double density, velocity[3];
  su2double vnx[4], vny[4], vnz[4];
  su2double graTet[3], graTetVisc[3], graTetSource[3];
  su2double Crd[4][3], Sens[4][nVarMetr][nFluxMetr], SensVisc[4][nVarMetr][nFluxMetr], SensSource[4][nVarMetr];
  su2double laminar_viscosity, eddy_viscosity;
  su2double k, omega, dk[3], domega[3];
  bool dummy_bool;

  const su2double sigmak1 = constants[0], sigmaom1 = constants[2],
                  sigmak2 = constants[1], sigmaom2 = constants[3];

  //--- note: currently only implemented for Tet

  for (iElem=0; iElem<nElem; ++iElem) {
    for (unsigned short iNode=0; iNode<4; ++iNode) {
      const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
      //--- store coordinates
      for (iDim = 0; iDim<3; ++iDim) {
        Crd[iNode][iDim] = geometry->node[kNode]->GetCoord(iDim);
      }// iDim
      //--- store sensors (goal-oriented)
      density     = solver_flow->GetNodes()->GetDensity(kNode);
      velocity[0] = solver_flow->GetNodes()->GetVelocity(kNode, 0);
      velocity[1] = solver_flow->GetNodes()->GetVelocity(kNode, 1);
      velocity[2] = solver_flow->GetNodes()->GetVelocity(kNode, 2);

      laminar_viscosity = solver_flow->GetNodes()->GetLaminarViscosity(kNode);
      eddy_viscosity    = solver_flow->GetNodes()->GetEddyViscosity(kNode);

      k     = nodes->GetSolution(kNode, 0);
      dk[0] = nodes->GetGradient(kNode, 0,0);
      dk[1] = nodes->GetGradient(kNode, 0,1);
      dk[2] = nodes->GetGradient(kNode, 0,2);

      omega     = nodes->GetSolution(kNode, 1);
      domega[0] = nodes->GetGradient(kNode, 1,0);
      domega[1] = nodes->GetGradient(kNode, 1,1);
      domega[2] = nodes->GetGradient(kNode, 1,2);

      const su2double F1 = nodes->GetF1blending(kNode);
      const su2double sigmak = F1*sigmak1 + (1.0 - F1)*sigmak2;
      const su2double sigmaomega = F1*sigmaom1 + (1.0 - F1)*sigmaom2;

      Sens[iNode][0][0] = density*velocity[0]*k;
      Sens[iNode][0][1] = density*velocity[1]*k;
      Sens[iNode][0][2] = density*velocity[2]*k;

      SensVisc[iNode][0][0] = (laminar_viscosity+sigmak*eddy_viscosity)*dk[0];
      SensVisc[iNode][0][1] = (laminar_viscosity+sigmak*eddy_viscosity)*dk[1];
      SensVisc[iNode][0][2] = (laminar_viscosity+sigmak*eddy_viscosity)*dk[2];

      Sens[iNode][1][0] = density*velocity[0]*omega;
      Sens[iNode][1][1] = density*velocity[1]*omega;
      Sens[iNode][1][2] = density*velocity[2]*omega;

      SensVisc[iNode][1][0] = (laminar_viscosity+sigmaomega*eddy_viscosity)*domega[0];
      SensVisc[iNode][1][1] = (laminar_viscosity+sigmaomega*eddy_viscosity)*domega[1];
      SensVisc[iNode][1][2] = (laminar_viscosity+sigmaomega*eddy_viscosity)*domega[2];

      SensSource[iNode][0] = nodes->GetSource(kNode)[0];
      SensSource[iNode][1] = nodes->GetSource(kNode)[1];
    }// iNode

    //--- inward face's normals : fac[0]=P1P2P3, fac[1]=P2P3P0, fac[2]=P3P0P1, fac[3]=P0P1P2
    vnx[0] = (Crd[2][1]-Crd[1][1])*(Crd[3][2]-Crd[1][2]) - (Crd[2][2]-Crd[1][2])*(Crd[3][1]-Crd[1][1]);
    vny[0] = (Crd[2][2]-Crd[1][2])*(Crd[3][0]-Crd[1][0]) - (Crd[2][0]-Crd[1][0])*(Crd[3][2]-Crd[1][2]);
    vnz[0] = (Crd[2][0]-Crd[1][0])*(Crd[3][1]-Crd[1][1]) - (Crd[2][1]-Crd[1][1])*(Crd[3][0]-Crd[1][0]);

    vnx[1] = (Crd[3][1]-Crd[2][1])*(Crd[0][2]-Crd[2][2]) - (Crd[3][2]-Crd[2][2])*(Crd[0][1]-Crd[2][1]);
    vny[1] = (Crd[3][2]-Crd[2][2])*(Crd[0][0]-Crd[2][0]) - (Crd[3][0]-Crd[2][0])*(Crd[0][2]-Crd[2][2]);
    vnz[1] = (Crd[3][0]-Crd[2][0])*(Crd[0][1]-Crd[2][1]) - (Crd[3][1]-Crd[2][1])*(Crd[0][0]-Crd[2][0]);

    vnx[2] = (Crd[0][1]-Crd[3][1])*(Crd[1][2]-Crd[3][2]) - (Crd[0][2]-Crd[3][2])*(Crd[1][1]-Crd[3][1]);
    vny[2] = (Crd[0][2]-Crd[3][2])*(Crd[1][0]-Crd[3][0]) - (Crd[0][0]-Crd[3][0])*(Crd[1][2]-Crd[3][2]);
    vnz[2] = (Crd[0][0]-Crd[3][0])*(Crd[1][1]-Crd[3][1]) - (Crd[0][1]-Crd[3][1])*(Crd[1][0]-Crd[3][0]);

    vnx[3] = (Crd[1][1]-Crd[0][1])*(Crd[2][2]-Crd[0][2]) - (Crd[1][2]-Crd[0][2])*(Crd[2][1]-Crd[0][1]);
    vny[3] = (Crd[1][2]-Crd[0][2])*(Crd[2][0]-Crd[0][0]) - (Crd[1][0]-Crd[0][0])*(Crd[2][2]-Crd[0][2]);
    vnz[3] = (Crd[1][0]-Crd[0][0])*(Crd[2][1]-Crd[0][1]) - (Crd[1][1]-Crd[0][1])*(Crd[2][0]-Crd[0][0]);

    //--- check if inward normal
    for(unsigned short iNode = 0; iNode < 4; ++iNode) {
      su2double CrdAvg[3] = {0.0, 0.0, 0.0};
      for(unsigned short jNode = 0; jNode < 4; ++jNode) {
        if(iNode != jNode) {
          CrdAvg[0] += Crd[jNode][0];
          CrdAvg[1] += Crd[jNode][1];
          CrdAvg[2] += Crd[jNode][2];
        }
      }// jNode
      CrdAvg[0] /= 3.;
      CrdAvg[1] /= 3.;
      CrdAvg[2] /= 3.;
      su2double u[3] = {CrdAvg[0]-Crd[iNode][0],
                        CrdAvg[1]-Crd[iNode][1],
                        CrdAvg[2]-Crd[iNode][2]};
      if((vnx[iNode]*u[0] + vny[iNode]*u[1] + vnz[iNode]*u[2]) > 0.) {
        vnx[iNode] *= -1.0;
        vny[iNode] *= -1.0;
        vnz[iNode] *= -1.0;
      }
    }// iNode

    //--- loop over conservative variables
    for(iVar = 0; iVar < nVarMetr; iVar++){

      //--- loop over directions
      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){

        //--- gradient at the element ( graTet = 6*|T|*gradT ) 
        graTet[0] = Sens[0][iVar][iFlux]*vnx[0] + Sens[1][iVar][iFlux]*vnx[1] + Sens[2][iVar][iFlux]*vnx[2] + Sens[3][iVar][iFlux]*vnx[3];
        graTet[1] = Sens[0][iVar][iFlux]*vny[0] + Sens[1][iVar][iFlux]*vny[1] + Sens[2][iVar][iFlux]*vny[2] + Sens[3][iVar][iFlux]*vny[3];
        graTet[2] = Sens[0][iVar][iFlux]*vnz[0] + Sens[1][iVar][iFlux]*vnz[1] + Sens[2][iVar][iFlux]*vnz[2] + Sens[3][iVar][iFlux]*vnz[3];

        graTetVisc[0] = SensVisc[0][iVar][iFlux]*vnx[0] + SensVisc[1][iVar][iFlux]*vnx[1] + SensVisc[2][iVar][iFlux]*vnx[2] + SensVisc[3][iVar][iFlux]*vnx[3];
        graTetVisc[1] = SensVisc[0][iVar][iFlux]*vny[0] + SensVisc[1][iVar][iFlux]*vny[1] + SensVisc[2][iVar][iFlux]*vny[2] + SensVisc[3][iVar][iFlux]*vny[3];
        graTetVisc[2] = SensVisc[0][iVar][iFlux]*vnz[0] + SensVisc[1][iVar][iFlux]*vnz[1] + SensVisc[2][iVar][iFlux]*vnz[2] + SensVisc[3][iVar][iFlux]*vnz[3];
    
        //--- assembling
        const unsigned short i = iFlux*nVarMetr*nDim + iVar*nDim;
        for (unsigned short iNode=0; iNode<4; ++iNode) {
          const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
          const su2double Vol = geometry->node[kNode]->GetVolume();
          const su2double rap = 1./(Vol*24.);
          nodes->AddAnisoGrad(kNode, i+0, graTet[0] * rap);
          nodes->AddAnisoGrad(kNode, i+1, graTet[1] * rap);
          nodes->AddAnisoGrad(kNode, i+2, graTet[2] * rap);

          nodes->AddAnisoViscGrad(kNode, i+0, graTetVisc[0] * rap);
          nodes->AddAnisoViscGrad(kNode, i+1, graTetVisc[1] * rap);
          nodes->AddAnisoViscGrad(kNode, i+2, graTetVisc[2] * rap);
        }// iNode
      }// iFlux

      graTetSource[0] = SensSource[0][iVar]*vnx[0] + SensSource[1][iVar]*vnx[1] + SensSource[2][iVar]*vnx[2] + SensSource[3][iVar]*vnx[3];
      graTetSource[1] = SensSource[0][iVar]*vny[0] + SensSource[1][iVar]*vny[1] + SensSource[2][iVar]*vny[2] + SensSource[3][iVar]*vny[3];
      graTetSource[2] = SensSource[0][iVar]*vnz[0] + SensSource[1][iVar]*vnz[1] + SensSource[2][iVar]*vnz[2] + SensSource[3][iVar]*vnz[3];
      //--- assembling
      const unsigned short i = iVar*nDim;
      for (unsigned short iNode=0; iNode<4; ++iNode) {
        const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
        const su2double Vol = geometry->node[kNode]->GetVolume();
        const su2double rap = 1./(Vol*24.);
        nodes->AddAnisoSourceGrad(kNode, i+0, graTetSource[0] * rap);
        nodes->AddAnisoSourceGrad(kNode, i+1, graTetSource[1] * rap);
        nodes->AddAnisoSourceGrad(kNode, i+2, graTetSource[2] * rap);
      }// iNode
    }// iVar
  }// iElem

  /*--- Communicate the gradient values via MPI. ---*/
  
  InitiateComms(geometry, config, ANISO_GRADIENT);
  CompleteComms(geometry, config, ANISO_GRADIENT);

  InitiateComms(geometry, config, ANISO_GRADIENT_VISC);
  CompleteComms(geometry, config, ANISO_GRADIENT_VISC);

  InitiateComms(geometry, config, ANISO_GRADIENT_SOURCE);
  CompleteComms(geometry, config, ANISO_GRADIENT_SOURCE);
}

void CTurbSSTSolver::SetHessian_L2Proj3(CGeometry *geometry, CConfig *config){

  unsigned long iPoint, nPoint = geometry->GetnPoint(), nPointDomain = geometry->GetnPointDomain(), iElem, nElem = geometry->GetnElem();
  unsigned short iVar, iFlux;
  unsigned short nVarMetr = 2, nFluxMetr = 3;  //--- TODO: adjust size of grad vector later for goal vs. feature
  unsigned short nMetr = 6;
  su2double vnx[4], vny[4], vnz[4];
  su2double hesTet[6], hesTetVisc[6], hesTetSource[6];
  su2double Crd[4][3], Grad[4][3][nVarMetr][nFluxMetr], GradVisc[4][3][nVarMetr][nFluxMetr], GradSource[4][3][nVarMetr];

  //--- note: currently only implemented for Tri

  for (iElem=0; iElem<nElem; ++iElem) {
    for (unsigned short iNode=0; iNode<4; ++iNode) {
      const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
      //--- store coordinates
      for (unsigned short iDim = 0; iDim<3; ++iDim) {
        Crd[iNode][iDim] = geometry->node[kNode]->GetCoord(iDim);
      }
      //--- store gradient
      for(iVar = 0; iVar < nVarMetr; iVar++){
        for(iFlux = 0; iFlux < nFluxMetr; iFlux++){
          const unsigned short i = iFlux*nVarMetr*nDim + iVar*nDim;
          Grad[iNode][0][iVar][iFlux] = nodes->GetAnisoGrad(kNode, i+0);
          Grad[iNode][1][iVar][iFlux] = nodes->GetAnisoGrad(kNode, i+1);
          Grad[iNode][2][iVar][iFlux] = nodes->GetAnisoGrad(kNode, i+2);

          GradVisc[iNode][0][iVar][iFlux] = nodes->GetAnisoViscGrad(kNode, i+0);
          GradVisc[iNode][1][iVar][iFlux] = nodes->GetAnisoViscGrad(kNode, i+1);
          GradVisc[iNode][2][iVar][iFlux] = nodes->GetAnisoViscGrad(kNode, i+2);
        }// iFlux

        const unsigned short i = iVar*nDim;
        GradSource[iNode][0][iVar] = nodes->GetAnisoSourceGrad(kNode, i+0);
        GradSource[iNode][1][iVar] = nodes->GetAnisoSourceGrad(kNode, i+1);
        GradSource[iNode][2][iVar] = nodes->GetAnisoSourceGrad(kNode, i+2);
      }// iVar
    }// iNode

    //--- inward face's normals : fac[0]=P1P2P3, fac[1]=P2P3P0, fac[2]=P3P0P1, fac[3]=P0P1P2
    vnx[0] = (Crd[2][1]-Crd[1][1])*(Crd[3][2]-Crd[1][2]) - (Crd[2][2]-Crd[1][2])*(Crd[3][1]-Crd[1][1]);
    vny[0] = (Crd[2][2]-Crd[1][2])*(Crd[3][0]-Crd[1][0]) - (Crd[2][0]-Crd[1][0])*(Crd[3][2]-Crd[1][2]);
    vnz[0] = (Crd[2][0]-Crd[1][0])*(Crd[3][1]-Crd[1][1]) - (Crd[2][1]-Crd[1][1])*(Crd[3][0]-Crd[1][0]);

    vnx[1] = (Crd[3][1]-Crd[2][1])*(Crd[0][2]-Crd[2][2]) - (Crd[3][2]-Crd[2][2])*(Crd[0][1]-Crd[2][1]);
    vny[1] = (Crd[3][2]-Crd[2][2])*(Crd[0][0]-Crd[2][0]) - (Crd[3][0]-Crd[2][0])*(Crd[0][2]-Crd[2][2]);
    vnz[1] = (Crd[3][0]-Crd[2][0])*(Crd[0][1]-Crd[2][1]) - (Crd[3][1]-Crd[2][1])*(Crd[0][0]-Crd[2][0]);

    vnx[2] = (Crd[0][1]-Crd[3][1])*(Crd[1][2]-Crd[3][2]) - (Crd[0][2]-Crd[3][2])*(Crd[1][1]-Crd[3][1]);
    vny[2] = (Crd[0][2]-Crd[3][2])*(Crd[1][0]-Crd[3][0]) - (Crd[0][0]-Crd[3][0])*(Crd[1][2]-Crd[3][2]);
    vnz[2] = (Crd[0][0]-Crd[3][0])*(Crd[1][1]-Crd[3][1]) - (Crd[0][1]-Crd[3][1])*(Crd[1][0]-Crd[3][0]);

    vnx[3] = (Crd[1][1]-Crd[0][1])*(Crd[2][2]-Crd[0][2]) - (Crd[1][2]-Crd[0][2])*(Crd[2][1]-Crd[0][1]);
    vny[3] = (Crd[1][2]-Crd[0][2])*(Crd[2][0]-Crd[0][0]) - (Crd[1][0]-Crd[0][0])*(Crd[2][2]-Crd[0][2]);
    vnz[3] = (Crd[1][0]-Crd[0][0])*(Crd[2][1]-Crd[0][1]) - (Crd[1][1]-Crd[0][1])*(Crd[2][0]-Crd[0][0]);

    //--- check if inward normal
    for(unsigned short iNode = 0; iNode < 4; ++iNode) {
      su2double CrdAvg[3] = {0.0, 0.0, 0.0};
      for(unsigned short jNode = 0; jNode < 4; ++jNode) {
        if(iNode != jNode) {
          CrdAvg[0] += Crd[jNode][0];
          CrdAvg[1] += Crd[jNode][1];
          CrdAvg[2] += Crd[jNode][2];
        }
      }// jNode
      CrdAvg[0] /= 3.;
      CrdAvg[1] /= 3.;
      CrdAvg[2] /= 3.;
      su2double u[3] = {CrdAvg[0]-Crd[iNode][0],
                        CrdAvg[1]-Crd[iNode][1],
                        CrdAvg[2]-Crd[iNode][2]};
      if((vnx[iNode]*u[0] + vny[iNode]*u[1] + vnz[iNode]*u[2]) > 0.) {
        vnx[iNode] *= -1.0;
        vny[iNode] *= -1.0;
        vnz[iNode] *= -1.0;
      }
    }// iNode

    //--- loop over conservative variables
    for(iVar = 0; iVar < nVarMetr; iVar++){

      //--- loop over directions
      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){

        //--- hessian at the element ( hesTet = 6*|T|*hessienT ) 
        hesTet[0] =         Grad[0][0][iVar][iFlux]*vnx[0] 
                          + Grad[1][0][iVar][iFlux]*vnx[1] 
                          + Grad[2][0][iVar][iFlux]*vnx[2]
                          + Grad[3][0][iVar][iFlux]*vnx[3];

        hesTet[1] = 0.5 * ( Grad[0][0][iVar][iFlux]*vny[0] 
                          + Grad[1][0][iVar][iFlux]*vny[1] 
                          + Grad[2][0][iVar][iFlux]*vny[2]
                          + Grad[3][0][iVar][iFlux]*vny[3]
                          + Grad[0][1][iVar][iFlux]*vnx[0] 
                          + Grad[1][1][iVar][iFlux]*vnx[1] 
                          + Grad[2][1][iVar][iFlux]*vnx[2]
                          + Grad[3][1][iVar][iFlux]*vnx[3] );

        hesTet[2] = 0.5 * ( Grad[0][0][iVar][iFlux]*vnz[0] 
                          + Grad[1][0][iVar][iFlux]*vnz[1] 
                          + Grad[2][0][iVar][iFlux]*vnz[2]
                          + Grad[3][0][iVar][iFlux]*vnz[3]
                          + Grad[0][2][iVar][iFlux]*vnx[0] 
                          + Grad[1][2][iVar][iFlux]*vnx[1] 
                          + Grad[2][2][iVar][iFlux]*vnx[2]
                          + Grad[3][2][iVar][iFlux]*vnx[3] );

        hesTet[3] =         Grad[0][1][iVar][iFlux]*vny[0] 
                          + Grad[1][1][iVar][iFlux]*vny[1] 
                          + Grad[2][1][iVar][iFlux]*vny[2]
                          + Grad[3][1][iVar][iFlux]*vny[3];

        hesTet[4] = 0.5 * ( Grad[0][1][iVar][iFlux]*vnz[0] 
                          + Grad[1][1][iVar][iFlux]*vnz[1] 
                          + Grad[2][1][iVar][iFlux]*vnz[2]
                          + Grad[3][1][iVar][iFlux]*vnz[3]
                          + Grad[0][2][iVar][iFlux]*vny[0] 
                          + Grad[1][2][iVar][iFlux]*vny[1] 
                          + Grad[2][2][iVar][iFlux]*vny[2]
                          + Grad[3][2][iVar][iFlux]*vny[3] );

        hesTet[5] =         Grad[0][2][iVar][iFlux]*vnz[0] 
                          + Grad[1][2][iVar][iFlux]*vnz[1] 
                          + Grad[2][2][iVar][iFlux]*vnz[2]
                          + Grad[3][2][iVar][iFlux]*vnz[3];

        hesTetVisc[0] =         GradVisc[0][0][iVar][iFlux]*vnx[0] 
                              + GradVisc[1][0][iVar][iFlux]*vnx[1] 
                              + GradVisc[2][0][iVar][iFlux]*vnx[2]
                              + GradVisc[3][0][iVar][iFlux]*vnx[3];

        hesTetVisc[1] = 0.5 * ( GradVisc[0][0][iVar][iFlux]*vny[0] 
                              + GradVisc[1][0][iVar][iFlux]*vny[1] 
                              + GradVisc[2][0][iVar][iFlux]*vny[2]
                              + GradVisc[3][0][iVar][iFlux]*vny[3]
                              + GradVisc[0][1][iVar][iFlux]*vnx[0] 
                              + GradVisc[1][1][iVar][iFlux]*vnx[1] 
                              + GradVisc[2][1][iVar][iFlux]*vnx[2]
                              + GradVisc[3][1][iVar][iFlux]*vnx[3] );

        hesTetVisc[2] = 0.5 * ( GradVisc[0][0][iVar][iFlux]*vnz[0] 
                              + GradVisc[1][0][iVar][iFlux]*vnz[1] 
                              + GradVisc[2][0][iVar][iFlux]*vnz[2]
                              + GradVisc[3][0][iVar][iFlux]*vnz[3]
                              + GradVisc[0][2][iVar][iFlux]*vnx[0] 
                              + GradVisc[1][2][iVar][iFlux]*vnx[1] 
                              + GradVisc[2][2][iVar][iFlux]*vnx[2]
                              + GradVisc[3][2][iVar][iFlux]*vnx[3] );

        hesTetVisc[3] =         GradVisc[0][1][iVar][iFlux]*vny[0] 
                              + GradVisc[1][1][iVar][iFlux]*vny[1] 
                              + GradVisc[2][1][iVar][iFlux]*vny[2]
                              + GradVisc[3][1][iVar][iFlux]*vny[3];

        hesTetVisc[4] = 0.5 * ( GradVisc[0][1][iVar][iFlux]*vnz[0] 
                              + GradVisc[1][1][iVar][iFlux]*vnz[1] 
                              + GradVisc[2][1][iVar][iFlux]*vnz[2]
                              + GradVisc[3][1][iVar][iFlux]*vnz[3]
                              + GradVisc[0][2][iVar][iFlux]*vny[0] 
                              + GradVisc[1][2][iVar][iFlux]*vny[1] 
                              + GradVisc[2][2][iVar][iFlux]*vny[2]
                              + GradVisc[3][2][iVar][iFlux]*vny[3] );

        hesTetVisc[5] =         GradVisc[0][2][iVar][iFlux]*vnz[0] 
                              + GradVisc[1][2][iVar][iFlux]*vnz[1] 
                              + GradVisc[2][2][iVar][iFlux]*vnz[2]
                              + GradVisc[3][2][iVar][iFlux]*vnz[3];
        
        //--- assembling
        const unsigned short i = iFlux*nVarMetr*nMetr + iVar*nMetr;
        for (unsigned short iNode=0; iNode<4; ++iNode) {
          const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
          const su2double Vol = geometry->node[kNode]->GetVolume();
          const su2double rap = 1./(Vol*24.);
          nodes->AddAnisoHess(kNode, i+0, hesTet[0] * rap);
          nodes->AddAnisoHess(kNode, i+1, hesTet[1] * rap);
          nodes->AddAnisoHess(kNode, i+2, hesTet[2] * rap);
          nodes->AddAnisoHess(kNode, i+3, hesTet[3] * rap);
          nodes->AddAnisoHess(kNode, i+4, hesTet[4] * rap);
          nodes->AddAnisoHess(kNode, i+5, hesTet[5] * rap);

          nodes->AddAnisoViscHess(kNode, i+0, hesTetVisc[0] * rap);
          nodes->AddAnisoViscHess(kNode, i+1, hesTetVisc[1] * rap);
          nodes->AddAnisoViscHess(kNode, i+2, hesTetVisc[2] * rap);
          nodes->AddAnisoViscHess(kNode, i+3, hesTetVisc[3] * rap);
          nodes->AddAnisoViscHess(kNode, i+4, hesTetVisc[4] * rap);
          nodes->AddAnisoViscHess(kNode, i+5, hesTetVisc[5] * rap);
        }// iNode
      }// iFlux

      hesTetSource[0] =         GradSource[0][0][iVar]*vnx[0] 
                              + GradSource[1][0][iVar]*vnx[1] 
                              + GradSource[2][0][iVar]*vnx[2]
                              + GradSource[3][0][iVar]*vnx[3];

      hesTetSource[1] = 0.5 * ( GradSource[0][0][iVar]*vny[0] 
                              + GradSource[1][0][iVar]*vny[1] 
                              + GradSource[2][0][iVar]*vny[2]
                              + GradSource[3][0][iVar]*vny[3]
                              + GradSource[0][1][iVar]*vnx[0] 
                              + GradSource[1][1][iVar]*vnx[1] 
                              + GradSource[2][1][iVar]*vnx[2]
                              + GradSource[3][1][iVar]*vnx[3] );

      hesTetSource[2] = 0.5 * ( GradSource[0][0][iVar]*vnz[0] 
                              + GradSource[1][0][iVar]*vnz[1] 
                              + GradSource[2][0][iVar]*vnz[2]
                              + GradSource[3][0][iVar]*vnz[3]
                              + GradSource[0][2][iVar]*vnx[0] 
                              + GradSource[1][2][iVar]*vnx[1] 
                              + GradSource[2][2][iVar]*vnx[2]
                              + GradSource[3][2][iVar]*vnx[3] );

      hesTetSource[3] =         GradSource[0][1][iVar]*vny[0] 
                              + GradSource[1][1][iVar]*vny[1] 
                              + GradSource[2][1][iVar]*vny[2]
                              + GradSource[3][1][iVar]*vny[3];

      hesTetSource[4] = 0.5 * ( GradSource[0][1][iVar]*vnz[0] 
                              + GradSource[1][1][iVar]*vnz[1] 
                              + GradSource[2][1][iVar]*vnz[2]
                              + GradSource[3][1][iVar]*vnz[3]
                              + GradSource[0][2][iVar]*vny[0] 
                              + GradSource[1][2][iVar]*vny[1] 
                              + GradSource[2][2][iVar]*vny[2]
                              + GradSource[3][2][iVar]*vny[3] );

      hesTetSource[5] =         GradSource[0][2][iVar]*vnz[0] 
                              + GradSource[1][2][iVar]*vnz[1] 
                              + GradSource[2][2][iVar]*vnz[2]
                              + GradSource[3][2][iVar]*vnz[3];
    
      //--- assembling
      const unsigned short i = iVar*nMetr;
      for (unsigned short iNode=0; iNode<3; ++iNode) {
        const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
        const su2double Vol = geometry->node[kNode]->GetVolume();
        const su2double rap = 1./(Vol*24.);
        nodes->AddAnisoSourceHess(kNode, i+0, hesTetSource[0] * rap);
        nodes->AddAnisoSourceHess(kNode, i+1, hesTetSource[1] * rap);
        nodes->AddAnisoSourceHess(kNode, i+2, hesTetSource[2] * rap);
        nodes->AddAnisoSourceHess(kNode, i+3, hesTetSource[3] * rap);
        nodes->AddAnisoSourceHess(kNode, i+4, hesTetSource[4] * rap);
        nodes->AddAnisoSourceHess(kNode, i+5, hesTetSource[5] * rap);
      }// iNode
    }// iVar
  }// iElem

  /*--- Communicate the Hessian values via MPI. ---*/
  
  InitiateComms(geometry, config, ANISO_HESSIAN);
  CompleteComms(geometry, config, ANISO_HESSIAN);
  InitiateComms(geometry, config, ANISO_HESSIAN_VISC);
  CompleteComms(geometry, config, ANISO_HESSIAN_VISC);
  InitiateComms(geometry, config, ANISO_HESSIAN_SOURCE);
  CompleteComms(geometry, config, ANISO_HESSIAN_SOURCE);

  CorrectBoundAnisoHess(geometry, config);
  CorrectBoundAnisoSourceHess(geometry, config);

  //--- Make positive definite matrix
  su2double **A      = new su2double*[nDim],
            **EigVec = new su2double*[nDim], 
            *EigVal  = new su2double[nDim];

  for(unsigned short iDim = 0; iDim < nDim; ++iDim){
    A[iDim]      = new su2double[nDim];
    EigVec[iDim] = new su2double[nDim];
  }

  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    for(iVar = 0; iVar < nVarMetr; iVar++){
      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){
        const unsigned short i = iFlux*nVarMetr*nMetr + iVar*nMetr;

        const su2double a = nodes->GetAnisoHess(iPoint, i+0);
        const su2double b = nodes->GetAnisoHess(iPoint, i+1);
        const su2double c = nodes->GetAnisoHess(iPoint, i+2);
        const su2double d = nodes->GetAnisoHess(iPoint, i+3);
        const su2double e = nodes->GetAnisoHess(iPoint, i+4);
        const su2double f = nodes->GetAnisoHess(iPoint, i+5);

        A[0][0] = a; A[0][1] = b; A[0][2] = c;
        A[1][0] = b; A[1][1] = d; A[1][2] = e;
        A[2][0] = c; A[2][1] = e; A[2][2] = f;

        CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

        for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = abs(EigVal[iDim]);

        CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

        nodes->SetAnisoHess(iPoint, i+0, A[0][0]);
        nodes->SetAnisoHess(iPoint, i+1, A[0][1]);
        nodes->SetAnisoHess(iPoint, i+2, A[0][2]);
        nodes->SetAnisoHess(iPoint, i+3, A[1][1]);
        nodes->SetAnisoHess(iPoint, i+4, A[1][2]);
        nodes->SetAnisoHess(iPoint, i+5, A[2][2]);
      }// iFlux

      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){
        const unsigned short i = iFlux*nVarMetr*nMetr + iVar*nMetr;

        const su2double a = nodes->GetAnisoViscHess(iPoint, i+0);
        const su2double b = nodes->GetAnisoViscHess(iPoint, i+1);
        const su2double c = nodes->GetAnisoViscHess(iPoint, i+2);
        const su2double d = nodes->GetAnisoViscHess(iPoint, i+3);
        const su2double e = nodes->GetAnisoViscHess(iPoint, i+4);
        const su2double f = nodes->GetAnisoViscHess(iPoint, i+5);

        A[0][0] = a; A[0][1] = b; A[0][2] = c;
        A[1][0] = b; A[1][1] = d; A[1][2] = e;
        A[2][0] = c; A[2][1] = e; A[2][2] = f;

        CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

        for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = abs(EigVal[iDim]);

        CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

        nodes->SetAnisoViscHess(iPoint, i+0, A[0][0]);
        nodes->SetAnisoViscHess(iPoint, i+1, A[0][1]);
        nodes->SetAnisoViscHess(iPoint, i+2, A[0][2]);
        nodes->SetAnisoViscHess(iPoint, i+3, A[1][1]);
        nodes->SetAnisoViscHess(iPoint, i+4, A[1][2]);
        nodes->SetAnisoViscHess(iPoint, i+5, A[2][2]);
      }// iFlux

      const unsigned short i = iVar*nMetr;

      const su2double a = nodes->GetAnisoSourceHess(iPoint, i+0);
      const su2double b = nodes->GetAnisoSourceHess(iPoint, i+1);
      const su2double c = nodes->GetAnisoSourceHess(iPoint, i+2);
      const su2double d = nodes->GetAnisoSourceHess(iPoint, i+3);
      const su2double e = nodes->GetAnisoSourceHess(iPoint, i+4);
      const su2double f = nodes->GetAnisoSourceHess(iPoint, i+5);

      A[0][0] = a; A[0][1] = b; A[0][2] = c;
      A[1][0] = b; A[1][1] = d; A[1][2] = e;
      A[2][0] = c; A[2][1] = e; A[2][2] = f;

      CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

      for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = abs(EigVal[iDim]);

      CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

      nodes->SetAnisoSourceHess(iPoint, i+0, A[0][0]);
      nodes->SetAnisoSourceHess(iPoint, i+1, A[0][1]);
      nodes->SetAnisoSourceHess(iPoint, i+2, A[0][2]);
      nodes->SetAnisoSourceHess(iPoint, i+3, A[1][1]);
      nodes->SetAnisoSourceHess(iPoint, i+4, A[1][2]);
      nodes->SetAnisoSourceHess(iPoint, i+5, A[2][2]);
    }// iVar
  }// iPoint
}
