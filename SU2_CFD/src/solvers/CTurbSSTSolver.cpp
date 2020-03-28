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
#include "../../include/gradients/computeGradientsGreenGauss.hpp"
#include "../../include/gradients/computeGradientsLeastSquares.hpp"
#include "../../include/limiters/computeLimiters.hpp"

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
    Primitive = new su2double[nVar];
    Primitive_i = new su2double[nVar]; Primitive_j = new su2double[nVar];

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
  constants[8] = constants[4]/constants[6] - constants[2]*0.41*0.41/sqrt(constants[6]);  //alfa_1
  constants[9] = constants[5]/constants[6] - constants[3]*0.41*0.41/sqrt(constants[6]);  //alfa_2
  
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

  unsigned long iPoint, ErrorCounter = 0;
  unsigned short iVar;

  bool limiter_turb = (config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (config->GetInnerIter() <= config->GetLimiterIter());

  /*--- Set the primitive variables ---*/

  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Initialize the residual vector ---*/

    LinSysRes.SetBlock_Zero(iPoint);

  }

  /*--- Initialize the Jacobian matrices ---*/

  Jacobian.SetValZero();

  /*--- Upwind second order reconstruction and gradients ---*/

  if (config->GetReconstructionGradientRequired()) {
    if (config->GetKind_Gradient_Method_Recon() == GREEN_GAUSS)
      SetPrimitive_Gradient_GG(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == LEAST_SQUARES)
      SetPrimitive_Gradient_LS(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES)
      SetPrimitive_Gradient_LS(geometry, config, true);
  }
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetPrimitive_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetPrimitive_Gradient_LS(geometry, config);

  if (limiter_turb) SetPrimitive_Limiter(geometry, config);

}

void CTurbSSTSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  su2double rho = 0.0, mu = 0.0, dist, omega, kine, F2, muT, zeta;
  su2double a1 = constants[7];
  unsigned long iPoint;

  /*--- Set the primitive variables ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    su2double rhokine  = nodes->GetSolution(iPoint, 0);
    su2double rhoomega = nodes->GetSolution(iPoint, 1);

    if((rhokine < 0.) || (rhoomega < 0.)) {
      rhokine  = nodes->GetSolution_Old(iPoint, 0);
      rhoomega = nodes->GetSolution_Old(iPoint, 1);
      nodes->SetSolution(iPoint, 0, rhokine);
      nodes->SetSolution(iPoint, 1, rhoomega);

    }

    rho = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    nodes->SetPrimitive(iPoint, 0, rhokine/rho);
    nodes->SetPrimitive(iPoint, 1, rhoomega/rho);

  }

  /*--- Compute mean flow and turbulence gradients ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
  }

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Compute blending functions and cross diffusion ---*/

    rho = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    mu  = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);

    dist = geometry->node[iPoint]->GetWall_Distance();

    su2double *Vorticity = solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint);
    su2double VorticityMag = sqrt(Vorticity[0]*Vorticity[0] +
                                  Vorticity[1]*Vorticity[1] +
                                  Vorticity[2]*Vorticity[2]);
    
//    su2double StrainMag = solver_container[FLOW_SOL]->GetNodes()->GetStrainMag(iPoint);

    nodes->SetBlendingFunc(iPoint,mu, dist, rho);

    F2 = nodes->GetF2blending(iPoint);

    /*--- Compute the eddy viscosity ---*/

    kine  = nodes->GetPrimitive(iPoint,0);
    omega = nodes->GetPrimitive(iPoint,1);
    zeta  = min(1.0/omega, a1/(VorticityMag*F2));
//  zeta  = min(1.0/omega, a1/(StrainMag*F2));
    muT   = max(rho*kine*zeta,0.0);
    nodes->SetmuT(iPoint,muT);

  }

}

unsigned long CTurbSSTSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {

  unsigned long iPoint, nonPhysicalPoints = 0;
  bool physical = true;

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    su2double rhokine  = nodes->GetSolution(iPoint, 0);
    su2double rhoomega = nodes->GetSolution(iPoint, 1);

    if((rhokine < 0.) || (rhoomega < 0.)) {
      rhokine  = nodes->GetSolution_Old(iPoint, 0);
      rhoomega = nodes->GetSolution_Old(iPoint, 1);
      nodes->SetSolution(iPoint, 0, rhokine);
      nodes->SetSolution(iPoint, 1, rhoomega);

      nonPhysicalPoints++;
    }

    su2double rho = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    nodes->SetPrimitive(iPoint, 0, rhokine/rho);
    nodes->SetPrimitive(iPoint, 1, rhoomega/rho);

  }

  return nonPhysicalPoints;
}

//--- This is hacky, fix later
void CTurbSSTSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config, bool reconstruction) {

  const auto& primitives = nodes->GetPrimitive();
  auto& gradient = reconstruction? nodes->GetGradient_Reconstruction() : nodes->GetGradient();

  computeGradientsGreenGauss(this, SOLUTION_GRADIENT, PERIODIC_SOL_GG, *geometry,
                             *config, primitives, 0, nVar, gradient);
}

void CTurbSSTSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config, bool reconstruction) {

  /*--- Set a flag for unweighted or weighted least-squares. ---*/
  bool weighted;

  if (reconstruction)
    weighted = (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES);
  else
    weighted = (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES);

  const auto& primitives = nodes->GetPrimitive();
  auto& rmatrix = nodes->GetRmatrix();
  auto& gradient = reconstruction? nodes->GetGradient_Reconstruction() : nodes->GetGradient();
  PERIODIC_QUANTITIES kindPeriodicComm = weighted? PERIODIC_SOL_LS : PERIODIC_SOL_ULS;

  computeGradientsLeastSquares(this, SOLUTION_GRADIENT, kindPeriodicComm, *geometry, *config,
                               weighted, primitives, 0, nVar, gradient, rmatrix);
}

void CTurbSSTSolver::SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) {

  auto kindLimiter = static_cast<ENUM_LIMITER>(config->GetKind_SlopeLimit());
  const auto& primitives = nodes->GetPrimitive();
  const auto& gradient = nodes->GetGradient_Reconstruction();
  auto& primMin = nodes->GetSolution_Min();
  auto& primMax = nodes->GetSolution_Max();
  auto& limiter = nodes->GetLimiter();

  computeLimiters(kindLimiter, this, SOLUTION_LIMITER, PERIODIC_LIM_SOL_1, PERIODIC_LIM_SOL_2,
            *geometry, *config, 0, nVar, primitives, gradient, primMin, primMax, limiter);
}
//--- End hacky bit

void CTurbSSTSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) {

  unsigned long iPoint;
  unsigned long iEdge, jPoint;
  unsigned short iNeigh, iDim;
  su2double *DivTurbVarGrad = new su2double[2], *Normal;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint), NULL);

    /*--- Gradient of the primitive and conservative variables ---*/

    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint), NULL);

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/

    numerics->SetTurbVar(nodes->GetPrimitive(iPoint), NULL);
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
    
    /*--- Divergence of turbulent variable gradients ---*/
    
//    DivTurbVarGrad[0] = DivTurbVarGrad[1] = 0.;
//    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
//      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
//      iEdge = geometry->FindEdge(iPoint,jPoint);
//      Normal = geometry->edge[iEdge]->GetNormal();
//      for (iDim = 0; iDim < nDim; iDim++) {
//        if (iPoint < jPoint) {
//          DivTurbVarGrad[0] += nodes->GetGradient(iPoint, 0, iDim)*Normal[iDim];
//          DivTurbVarGrad[1] += nodes->GetGradient(iPoint, 1, iDim)*Normal[iDim];
//        }
//        else {
//          DivTurbVarGrad[0] -= nodes->GetGradient(iPoint, 0, iDim)*Normal[iDim];
//          DivTurbVarGrad[1] -= nodes->GetGradient(iPoint, 1, iDim)*Normal[iDim];
//        }
//      }
//    }
//    numerics->SetDivTurbVarGrad(DivTurbVarGrad, NULL);

    /*--- Compute the source term ---*/

    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

  }
  
  delete [] DivTurbVarGrad;

}

void CTurbSSTSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {

}

void CTurbSSTSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar;
  su2double distance, density = 0.0, laminar_viscosity = 0.0, beta_1 = constants[4];

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

      Solution[0] = 0.0;
      Solution[1] = 60.0*laminar_viscosity/(beta_1*distance*distance+EPS*EPS);
      Solution[1] = min(max(Solution[1]/density, lowerlimit[1]), upperlimit[1])*density;

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

      // Solution[0] = 0.0;
      // Solution[1] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);
      Solution[0] = 0.0;
      Solution[1] = 60.0*laminar_viscosity/(beta_1*distance*distance);

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
      Primitive_i[iVar] = nodes->GetPrimitive(iPoint,iVar);

      Primitive_j[0] = kine_Inf;
      Primitive_j[1] = omega_Inf;

      conv_numerics->SetTurbVar(Primitive_i, Primitive_j);

      /*--- Set Normal (it is necessary to change the sign) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
      Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Set partial control volume ---*/
//      const unsigned long iEdge = geometry->FindEdge(iPoint,geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor());
//      conv_numerics->SetPartialVolume(geometry->edge[iEdge]->GetPartialVolume(0),
//                                      geometry->edge[iEdge]->GetPartialVolume(1));

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
        Primitive_i[iVar] = nodes->GetPrimitive(iPoint,iVar);

      /*--- Load the inlet turbulence variables (uniform by default). ---*/

      Primitive_j[0] = Inlet_TurbVars[val_marker][iVertex][0];
      Primitive_j[1] = Inlet_TurbVars[val_marker][iVertex][1];

      conv_numerics->SetTurbVar(Primitive_i, Primitive_j);

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
      //     visc_numerics->SetTurbVar(Primitive_i, Primitive_j);
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
       Primitive_i --> TurbVar_internal,
       Primitive_j --> TurbVar_outlet ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        Primitive_i[iVar] = nodes->GetPrimitive(iPoint,iVar);
        Primitive_j[iVar] = nodes->GetPrimitive(iPoint,iVar);
      }
      conv_numerics->SetTurbVar(Primitive_i, Primitive_j);

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
//      visc_numerics->SetTurbVar(Primitive_i, Primitive_j);
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
        Primitive_i[iVar] = nodes->GetPrimitive(iPoint,iVar);

      Primitive_j[0]= extAverageKine;
      Primitive_j[1]= extAverageOmega;

      conv_numerics->SetTurbVar(Primitive_i, Primitive_j);

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
      visc_numerics->SetTurbVar(Primitive_i, Primitive_j);
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
        Primitive_i[iVar] = nodes->GetPrimitive(iPoint,iVar);

      /*--- Set the turbulent variable states. Use average span-wise values
             values for the turbulent state at the inflow. ---*/

      Primitive_j[0]= kine_b;
      Primitive_j[1]= omega_b;

      conv_numerics->SetTurbVar(Primitive_i, Primitive_j);

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
      visc_numerics->SetTurbVar(Primitive_i, Primitive_j);
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
            Primitive_i[0] = nodes->GetPrimitive(iPoint,0);
            Primitive_i[1] = nodes->GetPrimitive(iPoint,1);

            Primitive_j[0] = GetSlidingState(iMarker, iVertex, 0, jVertex);
            Primitive_j[1] = GetSlidingState(iMarker, iVertex, 1, jVertex);

            conv_numerics->SetTurbVar(Primitive_i, Primitive_j);

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

          visc_numerics->SetTurbVar(Primitive_i, Primitive_j);
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

void CTurbSSTSolver::Correct_Omega_WF(CGeometry      *geometry,
                                      CSolver        **solver,
                                      CConfig        *config) {
  
  unsigned long jPoint, total_index;
  unsigned short iVar, kNode;
  long iElem;
  su2double distance, density = 0.0, laminar_viscosity = 0.0, eddy_viscosity = 0.0, k = 0.0, beta_1 = constants[4];
  su2double *weights;
  
  /*--- Communicate values needed for WF ---*/
  WF_Comms(geometry, solver, config);

  for (jPoint = 0; jPoint < nPointDomain; jPoint++) {
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[jPoint]->GetBool_Wall_Neighbor()) {
      
      iElem = geometry->node[jPoint]->GetWall_Element();
      
      density = solver[FLOW_SOL]->GetNodes()->GetDensity(jPoint);
      laminar_viscosity = solver[FLOW_SOL]->GetNodes()->GetLaminarViscosity(jPoint);
      eddy_viscosity = solver[FLOW_SOL]->GetNodes()->GetEddyViscosity(jPoint);
      k = nodes->GetPrimitive(jPoint, 0);
      
      weights = geometry->node[jPoint]->GetWall_Interpolation_Weights();
      
      distance = geometry->node[jPoint]->GetWall_Distance();
      
      su2double Omega_0 = sqrt(k) / (pow(0.09,0.25) * 0.41 * distance + EPS);
      su2double Omega = 0.0;
      
      for (kNode = 0; kNode < 4; kNode++) {
        
        const su2double DensityWall = nodes->GetWallDensity(jPoint, kNode);
        const su2double LamViscWall = nodes->GetWallLamVisc(jPoint, kNode);
        
        const su2double Omega_i = 6. * LamViscWall / (beta_1 * DensityWall * pow(distance, 2.0) + EPS*EPS);
        Omega += sqrt(pow(Omega_0, 2.) + pow(Omega_i, 2.))*weights[kNode];
      }
      Solution[1] = density*Omega;

      for (iVar = 1; iVar < nVar; iVar++) {
        Solution[iVar] = min(max(Solution[iVar]/density, lowerlimit[iVar]), upperlimit[iVar])*density;
        nodes->SetSolution_Old(jPoint,iVar,Solution[iVar]);
        nodes->SetSolution(jPoint,iVar,Solution[iVar]);
        LinSysRes.SetBlock_Zero(jPoint,iVar);
        
        /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
        total_index = jPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    }
  }
  
}

void CTurbSSTSolver::WF_Comms(CGeometry *geometry,
                              CSolver   **solver,
                              CConfig   *config) {
  
  unsigned long iPoint;
  unsigned long counter = 0;
  
  /*--- Initialize the wall index map and solution matrices on the first iteration. ---*/
  if (config->GetInnerIter() == 0) {
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (geometry->node[iPoint]->GetBool_Wall_Neighbor()) {
        nodes->SetWallMap(iPoint,counter);
        counter++;
      }
    }
    if (counter > 0) nodes->InitializeWallSolution(counter);
    counter = 0;
  }
  
#ifdef HAVE_MPI
  
  /*--- Begin the communication. ---*/
  
  /*--- Local variables. ---*/

  int iRank, iSend, iRecv;
  
  unsigned short countPerElem = 8;
  
  int *nElemSend = new int[size+1]; nElemSend[0] = 0;
  int *nElemRecv = new int[size+1]; nElemRecv[0] = 0;

  for (iRank = 0; iRank < size; iRank++) {
    nElemSend[iRank] = 0; nElemRecv[iRank] = 0;
  }
  nElemSend[size] = 0; nElemRecv[size] = 0;

  /*--- Loop through all of our nodes and track our sends with each rank. ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    if (geometry->node[iPoint]->GetBool_Wall_Neighbor()) {

      /*--- Get the source rank and number of points to send. ---*/

      iRank = geometry->node[iPoint]->GetWall_Rank();
      nElemSend[iRank+1]++;

    }
  }
  
  /*--- Communicate the number of points to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will recv from each other processor. ---*/
  
  SU2_MPI::Alltoall(&(nElemSend[1]), 1, MPI_INT,
                    &(nElemRecv[1]), 1, MPI_INT, MPI_COMM_WORLD);
  
  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
    
  unsigned short nSend = 0, nRecv = 0;

  for (iRank = 0; iRank < size; iRank++) {
    if ((iRank != rank) && (nElemSend[iRank+1] > 0)) nSend++;
    if ((iRank != rank) && (nElemRecv[iRank+1] > 0)) nRecv++;

    nElemSend[iRank+1] += nElemSend[iRank];
    nElemRecv[iRank+1] += nElemRecv[iRank];
  }
    
  /*-- Allocate our communication memory. Note that the "ElemSend" arrays
   correspond to the wall element indices we sent out. The variables
   being sent back thus correspond to "ElemRecv" arrays. ---*/
  
  su2double *bufDSend = new su2double[countPerElem*nElemRecv[size]];
  for (iRecv = 0; iRecv < countPerElem*nElemRecv[size]; iRecv++)
    bufDSend[iRecv] = 0.0;

  su2double *bufDRecv = new su2double[countPerElem*nElemSend[size]];
  for (iSend = 0; iSend < countPerElem*nElemSend[size]; iSend++)
    bufDRecv[iSend] = 0.0;
  
  unsigned long *bufLSend = new unsigned long[nElemSend[size]];
  for (iSend = 0; iSend < nElemSend[size]; iSend++)
    bufLSend[iSend] = 0;

  unsigned long *bufLRecv = new unsigned long[nElemRecv[size]];
  for (iRecv = 0; iRecv < nElemRecv[size]; iRecv++)
    bufLRecv[iRecv] = 0;
  
  unsigned short *bufSSend = new unsigned short[nElemSend[size]];
  for (iSend = 0; iSend < nElemSend[size]; iSend++)
    bufSSend[iSend] = 0;

  unsigned short *bufSRecv = new unsigned short[nElemRecv[size]];
  for (iRecv = 0; iRecv < nElemRecv[size]; iRecv++)
    bufSRecv[iRecv] = 0;
  
  /*--- Allocate memory for the MPI requests if we need to communicate. ---*/
    
  SU2_MPI::Request *sendReq = NULL, *recvReq = NULL;

  if (nSend > 0) {
    sendReq = new SU2_MPI::Request[nSend];
  }
  if (nRecv > 0) {
    recvReq = new SU2_MPI::Request[nRecv];
  }
    
  /*--- Local variables ---*/
  
  int iMessage, iProc, offset, nElem, count, source, dest, tag;
  unsigned short commType;
  
  for (unsigned short i = 0; i < 3; i++) {
    
    /*--- Elements ---*/
    if (i == 0) {
      commType = COMM_TYPE_UNSIGNED_LONG;
      countPerElem = 1;
      unsigned long ProcCounter[size];
      for (iProc = 0; iProc < size; iProc++) ProcCounter[iProc] = 0;
      
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
        if (geometry->node[iPoint]->GetBool_Wall_Neighbor()) {
          const unsigned short RankID = geometry->node[iPoint]->GetWall_Rank();
          if (RankID != rank) {
            offset = nElemSend[RankID]+ProcCounter[RankID];
            bufLSend[offset] = geometry->node[iPoint]->GetWall_Element();
            bufSSend[offset] = geometry->node[iPoint]->GetWall_Marker();
            ProcCounter[RankID]++;
          }
        }
      }
    }
  
    /*--- Markers ---*/
    else if (i == 1) {
      commType = COMM_TYPE_UNSIGNED_SHORT;
      countPerElem = 1;
    }
    
    /*--- Variables ---*/
    else if (i == 2) {
      commType = COMM_TYPE_DOUBLE;
      countPerElem = 8;
      
      /*--- Fill send buffers with variables based on received
            markers and elements. ---*/
      for (iProc = 0; iProc < size; iProc++) {
        if ((nElemRecv[iProc+1] > nElemRecv[iProc]) && (iProc != rank)) {
          nElem = nElemRecv[iProc+1] - nElemRecv[iProc];
          for (unsigned long iElem = 0; iElem < nElem; iElem++) {
            offset = countPerElem*(nElemRecv[iProc]+iElem);
            
            const unsigned long ElemID = bufLRecv[nElemRecv[iProc]+iElem];
            const unsigned short MarkerID = bufSRecv[nElemRecv[iProc]+iElem];
            const unsigned short nNodeElem = geometry->bound[MarkerID][ElemID]->GetnNodes();
            for (unsigned short kNode = 0; kNode < nNodeElem; kNode++) {
              const unsigned long kPoint = geometry->bound[MarkerID][ElemID]->GetNode(kNode);
              bufDSend[offset+kNode]   = solver[FLOW_SOL]->GetNodes()->GetDensity(kPoint);
              bufDSend[offset+kNode+4] = solver[FLOW_SOL]->GetNodes()->GetLaminarViscosity(kPoint);
            }
          }
        }
      }
    }
  
    /*--- Launch the non-blocking recv's first. ---*/
    
    iMessage = 0;
    
    for (iProc = 0; iProc < size; iProc++) {
      
      /*--- Post recv's only if another proc is sending us data. We do
       not communicate with ourselves or post recv's for zero length
       messages to keep overhead down. ---*/
      
      if ((((nElemRecv[iProc+1] > nElemRecv[iProc]) && (commType != COMM_TYPE_DOUBLE)) ||
           ((nElemSend[iProc+1] > nElemSend[iProc]) && (commType == COMM_TYPE_DOUBLE))) &&
          (iProc != rank)) {
        
        /*--- Compute our location in the recv buffer. ---*/
        
        if (commType == COMM_TYPE_DOUBLE) offset = countPerElem*nElemSend[iProc];
        else                              offset = countPerElem*nElemRecv[iProc];
        
        /*--- Take advantage of cumulative storage format to get the number
         of elems that we need to recv. ---*/
        
        if (commType == COMM_TYPE_DOUBLE) nElem = nElemSend[iProc+1] - nElemSend[iProc];
        else                              nElem = nElemRecv[iProc+1] - nElemRecv[iProc];
        
        /*--- Total count can include multiple pieces of data per element. ---*/
        
        count = countPerElem*nElem;
        
        /*--- Post non-blocking recv for this proc. ---*/
        
        source = iProc; tag = iProc + 1;
        
        switch (commType) {
          case COMM_TYPE_DOUBLE:
            SU2_MPI::Irecv(&(static_cast<su2double*>(bufDRecv)[offset]),
                           count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD,
                           &(sendReq[iMessage]));
            break;
          case COMM_TYPE_UNSIGNED_LONG:
            SU2_MPI::Irecv(&(static_cast<unsigned long*>(bufLRecv)[offset]),
                           count, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD,
                           &(recvReq[iMessage]));
            break;
          case COMM_TYPE_UNSIGNED_SHORT:
            SU2_MPI::Irecv(&(static_cast<unsigned short*>(bufSRecv)[offset]),
                           count, MPI_UNSIGNED_SHORT, source, tag, MPI_COMM_WORLD,
                           &(recvReq[iMessage]));
            break;
          default:
            break;
        }
        
        /*--- Increment message counter. ---*/
        
        iMessage++;
        
      }
    }
    
    /*--- Launch the non-blocking sends next. ---*/
    
    iMessage = 0;
    for (iProc = 0; iProc < size; iProc++) {
      
      /*--- Post sends only if we are sending another proc data. We do
       not communicate with ourselves or post sends for zero length
       messages to keep overhead down. ---*/
      
      if ((((nElemSend[iProc+1] > nElemSend[iProc]) && (commType != COMM_TYPE_DOUBLE)) ||
           ((nElemRecv[iProc+1] > nElemRecv[iProc]) && (commType == COMM_TYPE_DOUBLE))) &&
          (iProc != rank)) {
        
        /*--- Compute our location in the send buffer. ---*/
        
        if (commType == COMM_TYPE_DOUBLE) offset = countPerElem*nElemRecv[iProc];
        else                              offset = countPerElem*nElemSend[iProc];
        
        /*--- Take advantage of cumulative storage format to get the number
         of elems that we need to send. ---*/
        
        if (commType == COMM_TYPE_DOUBLE) nElem = nElemRecv[iProc+1] - nElemRecv[iProc];
        else                              nElem = nElemSend[iProc+1] - nElemSend[iProc];
        
        /*--- Total count can include multiple pieces of data per element. ---*/
        
        count = countPerElem*nElem;
        
        /*--- Post non-blocking send for this proc. ---*/
        
        dest = iProc; tag = rank + 1;
        
        switch (commType) {
          case COMM_TYPE_DOUBLE:
            SU2_MPI::Isend(&(static_cast<su2double*>(bufDSend)[offset]),
                           count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD,
                           &(recvReq[iMessage]));
            break;
          case COMM_TYPE_UNSIGNED_LONG:
            SU2_MPI::Isend(&(static_cast<unsigned long*>(bufLSend)[offset]),
                           count, MPI_UNSIGNED_LONG, dest, tag, MPI_COMM_WORLD,
                           &(sendReq[iMessage]));
            break;
          case COMM_TYPE_UNSIGNED_SHORT:
            SU2_MPI::Isend(&(static_cast<unsigned short*>(bufSSend)[offset]),
                           count, MPI_UNSIGNED_SHORT, dest, tag, MPI_COMM_WORLD,
                           &(sendReq[iMessage]));
            break;
          default:
            break;
        }
        
        /*--- Increment message counter. ---*/
        
        iMessage++;
        
      }
    }
      
    int ind;
    SU2_MPI::Status status;

    /*--- Wait for the non-blocking sends to complete. ---*/

    for (iSend = 0; iSend < nSend; iSend++)
      SU2_MPI::Waitany(nSend, sendReq, &ind, &status);

    /*--- Wait for the non-blocking recvs to complete. ---*/

    for (iRecv = 0; iRecv < nRecv; iRecv++)
      SU2_MPI::Waitany(nRecv, recvReq, &ind, &status);
    
  }
  
  if (sendReq != NULL) delete [] sendReq;
  if (recvReq != NULL) delete [] recvReq;
  
  /*--- Now that the wall elements have been communicated, store them. ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    if (geometry->node[iPoint]->GetBool_Wall_Neighbor()) {
      const unsigned short RankID = geometry->node[iPoint]->GetWall_Rank();
      if ((nElemSend[RankID+1] > nElemSend[RankID]) && (RankID != rank)) {
        nElem = nElemSend[RankID+1] - nElemSend[RankID];
        for (unsigned long iElem = 0; iElem < nElem; iElem++) {
          offset = nElemSend[RankID]+iElem;
          const unsigned long ElemID = bufLSend[offset];
          const unsigned short MarkerID = bufSSend[offset];
          if ((geometry->node[iPoint]->GetWall_Marker() == MarkerID) &&
              (geometry->node[iPoint]->GetWall_Element() == ElemID)) {
            for (unsigned short kNode = 0; kNode < 4; kNode++) {
              nodes->SetWallDensity(iPoint, kNode, bufDRecv[countPerElem*offset+kNode]);
              nodes->SetWallLamVisc(iPoint, kNode, bufDRecv[countPerElem*offset+kNode+4]);
            }
          }
        }
      }
      else if (RankID == rank) {
        const unsigned short MarkerID = geometry->node[iPoint]->GetWall_Marker();
        const unsigned long ElemID = geometry->node[iPoint]->GetWall_Element();
        for (unsigned short kNode = 0; kNode < geometry->bound[MarkerID][ElemID]->GetnNodes(); kNode++) {
          const unsigned long kPoint = geometry->bound[MarkerID][ElemID]->GetNode(kNode);
          nodes->SetWallDensity(iPoint, kNode, solver[FLOW_SOL]->GetNodes()->GetDensity(kPoint));
          nodes->SetWallLamVisc(iPoint, kNode, solver[FLOW_SOL]->GetNodes()->GetLaminarViscosity(kPoint));
        }
      }
    }
  }
  
  /*--- Free some memory. ---*/
  
  delete [] nElemSend;
  delete [] nElemRecv;
  delete [] bufLSend;
  delete [] bufLRecv;
  delete [] bufSSend;
  delete [] bufSRecv;
  delete [] bufDSend;
  delete [] bufDRecv;
  
#else

  /*--- Serial mode. Just store the values. Easy game. ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    if (geometry->node[iPoint]->GetBool_Wall_Neighbor()) {
      const unsigned short MarkerID = geometry->node[iPoint]->GetWall_Marker();
      const unsigned long ElemID = geometry->node[iPoint]->GetWall_Element();
      for (unsigned short kNode = 0; kNode < geometry->bound[MarkerID][ElemID]->GetnNodes(); kNode++) {
        const unsigned long kPoint = geometry->bound[MarkerID][ElemID]->GetNode(kNode);
        nodes->SetWallDensity(iPoint, kNode, solver[FLOW_SOL]->GetNodes()->GetDensity(kPoint));
        nodes->SetWallLamVisc(iPoint, kNode, solver[FLOW_SOL]->GetNodes()->GetLaminarViscosity(kPoint));
      }
    }
  }

#endif
  
}

void CTurbSSTSolver::TurbulentMetric(CSolver           **solver,
                                     CGeometry         *geometry,
                                     CConfig           *config,
                                     unsigned long     iPoint,
                                     vector<su2double> &weights) {
  CVariable *varFlo    = solver[FLOW_SOL]->GetNodes(),
            *varTur    = solver[TURB_SOL]->GetNodes(),
            *varAdjFlo = solver[ADJFLOW_SOL]->GetNodes(),
            *varAdjTur = solver[ADJTURB_SOL]->GetNodes();

  unsigned short iDim, jDim, iVar;
  const unsigned short nMetr = 3*(nDim-1);
  const unsigned short nVarFlo = solver[FLOW_SOL]->GetnVar();
  const unsigned short nVarTur = solver[TURB_SOL]->GetnVar();

  //--- First-order terms (error due to viscosity)
  su2double r, u[3], e, k, omega,
            T, mu, mut, lam, lamt, dmudT, 
            Tref, S, muref, 
            R, cv, cp, g, Pr, Prt;

  r = varFlo->GetDensity(iPoint);
  u[0] = varFlo->GetVelocity(iPoint, 0);
  u[1] = varFlo->GetVelocity(iPoint, 1);
  if (nDim == 3) u[2] = varFlo->GetVelocity(iPoint, 2);
  e = varFlo->GetEnergy(iPoint);
  k = varTur->GetPrimitive(iPoint, 0);
  omega = varTur->GetPrimitive(iPoint, 1);
  
  T   = varFlo->GetTemperature(iPoint);
  mu  = varFlo->GetLaminarViscosity(iPoint);
  mut = varFlo->GetEddyViscosity(iPoint);

  Tref  = config->GetMu_Temperature_RefND();
  S     = config->GetMu_SND();
  muref = config->GetMu_RefND();
  dmudT = muref*(Tref+S)/pow(Tref,1.5) * (3.*S*sqrt(T) + pow(T,1.5))/(2.*pow((T+S),2.));

  g    = config->GetGamma();
  R    = config->GetGas_ConstantND();
  cp   = (g/(g-1.))*R;
  cv   = cp/g;
  Pr   = config->GetPrandtl_Lam();
  Prt  = config->GetPrandtl_Turb();
  lam  = cp*mu/Pr;
  lamt = cp*mut/Prt;

  su2double gradu[3][3], gradT[3], gradk[3], gradomega[3], divu, taut[3][3],
            delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0 ; jDim < nDim; jDim++) {
      gradu[iDim][jDim] = varFlo->GetGradient_Primitive(iPoint, iDim+1, jDim);
    }
    gradT[iDim]     = varFlo->GetGradient_Primitive(iPoint, 0, iDim);
    gradk[iDim]     = varTur->GetGradient(iPoint, 0, iDim);
    gradomega[iDim] = varTur->GetGradient(iPoint, 1, iDim);
  }

  divu = 0.0; for (iDim = 0 ; iDim < nDim; ++iDim) divu += gradu[iDim][iDim];

  for (iDim = 0; iDim < nDim; ++iDim) {
    for (jDim = 0; jDim < nDim; ++jDim) {
      taut[iDim][jDim] = mut*( gradu[jDim][iDim] + gradu[iDim][jDim] ) 
                       - (2./3.)*mut*divu*delta[iDim][jDim]
                       - (2./3.)*r*k*delta[iDim][jDim];
    }
  }

  const su2double F1 = varTur->GetF1blending(iPoint);

  const su2double alfa        = F1*constants[8] + (1.0 - F1)*constants[9];
  const su2double sigmak      = F1*constants[0] + (1.0 - F1)*constants[1];
  const su2double sigmaomega  = F1*constants[2] + (1.0 - F1)*constants[3];
  const su2double sigmaomega2 = constants[3];
  const su2double beta        = F1*constants[4] + (1.0 - F1)*constants[5];
  const su2double betastar    = constants[6];

  //--- Momentum weights
  vector<su2double> TmpWeights(weights.size(), 0.0);
  su2double factor = 0.0;
  for (iDim = 0; iDim < nDim; ++iDim) {
    factor = -(2./3.)*(mut/r*varAdjTur->GetAnisoGrad(iPoint, 0*nDim+iDim)
                                +alfa*varAdjTur->GetAnisoGrad(iPoint, 1*nDim+iDim));
    for (jDim = 0; jDim < nDim; ++jDim) {
      factor += (taut[iDim][jDim]+mut*gradu[iDim][jDim])*(1./r*varAdjTur->GetAnisoGrad(iPoint, 0*nDim+jDim)
                                                         +alfa/mut*varAdjTur->GetAnisoGrad(iPoint, 1*nDim+jDim));
    }
    TmpWeights[iDim+1] += factor;
  }

  //--- k and omega weights
  factor = 0.0;
  for (iDim = 0; iDim < nDim; ++iDim) {
    for (jDim = 0; jDim < nDim; ++jDim) {
      iVar = iDim+1;
      factor += (taut[iDim][jDim]+(2./3.)*r*k*delta[iDim][jDim])/mut
              * (varAdjFlo->GetAnisoGrad(iPoint, iVar*nDim+jDim)
              + u[jDim]*varAdjFlo->GetAnisoGrad(iPoint, (nVarFlo-1)*nDim+iDim));
    }
    factor += cp/Prt*gradT[iDim]*varAdjFlo->GetAnisoGrad(iPoint, (nVarFlo-1)*nDim+iDim);
    factor += sigmak*gradk[iDim]*varAdjTur->GetAnisoGrad(iPoint, 0*nDim+iDim)
            + sigmaomega*gradomega[iDim]*varAdjTur->GetAnisoGrad(iPoint, 1*nDim+iDim);
  }

  TmpWeights[nVarFlo+0] += 1./omega*factor;
  TmpWeights[nVarFlo+1] += -k/pow(omega,2.)*factor;
  for (iDim = 0; iDim < nDim; ++iDim) {
    TmpWeights[nVarFlo+0] += 2.*(1.-F1)*sigmaomega2/omega*gradomega[iDim]*varAdjTur->GetAnisoGrad(iPoint, 1*nDim+iDim);
    TmpWeights[nVarFlo+1] += 2.*(1.-F1)*sigmaomega2/omega*gradk[iDim]*varAdjTur->GetAnisoGrad(iPoint, 1*nDim+iDim);
  }

  //--- Density weight
  for (iDim = 0; iDim < nDim; ++iDim) TmpWeights[0] -= u[iDim]*TmpWeights[iDim+1];
  TmpWeights[0] -= k*TmpWeights[nVarFlo+0] + omega*TmpWeights[nVarFlo+1]
                 + k/omega*factor;

  //--- Add TmpWeights to weights, then reset for second-order terms
  for (iVar = 0; iVar < nVarFlo+nVarTur; ++iVar) weights[iVar] += TmpWeights[iVar];
  fill(TmpWeights.begin(), TmpWeights.end(), 0.0);

  //--- Second-order terms (error due to gradients)
  if(nDim == 3) {
    const unsigned short rki = 0*nMetr, romegai = 1*nMetr, xxi = 0, yyi = 3, zzi = 5;
    TmpWeights[nVarFlo+0] -= (mu+sigmak*mut)*(varAdjTur->GetAnisoHess(iPoint, rki+xxi)
                                             +varAdjTur->GetAnisoHess(iPoint, rki+yyi)
                                             +varAdjTur->GetAnisoHess(iPoint, rki+zzi)); // Hk
    TmpWeights[nVarFlo+1] -= (mu+sigmaomega*mut)*(varAdjTur->GetAnisoHess(iPoint, romegai+xxi)
                                                 +varAdjTur->GetAnisoHess(iPoint, romegai+yyi)
                                                 +varAdjTur->GetAnisoHess(iPoint, romegai+zzi)); // Homega

  }
  else {
    const unsigned short rki = 0*nMetr, romegai = 1*nMetr, xxi = 0, yyi = 2;
    TmpWeights[nVarFlo+0] -= (mu+sigmak*mut)*(varAdjTur->GetAnisoHess(iPoint, rki+xxi)
                                             +varAdjTur->GetAnisoHess(iPoint, rki+yyi)); // Hk
    TmpWeights[nVarFlo+1] -= (mu+sigmaomega*mut)*(varAdjTur->GetAnisoHess(iPoint, romegai+xxi)
                                                 +varAdjTur->GetAnisoHess(iPoint, romegai+yyi)); // Homega
  }
  TmpWeights[0] -= k*TmpWeights[nVarFlo+0]+omega*TmpWeights[nVarFlo+1];

  //--- Add TmpWeights to weights
  weights[0]         += TmpWeights[0];
  weights[nVarFlo+0] += TmpWeights[nVarFlo+0];
  weights[nVarFlo+1] += TmpWeights[nVarFlo+1];

  //--- Zeroth-order terms (due to production and dissipation)
  weights[0]         += -betastar*k*omega*varAdjTur->GetSolution(iPoint,0)
                      - beta*pow(omega,2.)*varAdjTur->GetSolution(iPoint,1);
  weights[nVarFlo+0] += betastar*omega*varAdjTur->GetSolution(iPoint,0)
                      + (2./3.)*divu*varAdjTur->GetSolution(iPoint,0);
  weights[nVarFlo+1] += betastar*k*varAdjTur->GetSolution(iPoint,0)
                      + 2.*beta*omega*varAdjTur->GetSolution(iPoint,1)
                      + (2./3.)*alfa*divu*varAdjTur->GetSolution(iPoint,1);

}
