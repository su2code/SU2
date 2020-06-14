/*!
 * \file CTurbSSTSolver.cpp
 * \brief Main subrotuines of CTurbSSTSolver class
 * \author F. Palacios, A. Bueno
 * \version 7.0.3 "Blackbird"
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

#include "../../include/solvers/CTurbSSTSolver.hpp"
#include "../../include/variables/CTurbSSTVariable.hpp"
#include "../../include/gradients/computeGradientsGreenGauss.hpp"
#include "../../include/gradients/computeGradientsLeastSquares.hpp"
#include "../../include/limiters/computeLimiters.hpp"
#include "../../../Common/include/omp_structure.hpp"

CTurbSSTSolver::CTurbSSTSolver(void) : CTurbSolver() { }

CTurbSSTSolver::CTurbSSTSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
    : CTurbSolver(geometry, config) {
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  ifstream restart_file;
  string text_line;

  bool multizone = config->GetMultizone_Problem();

  /*--- Array initialization ---*/

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

    /*--- Define some auxiliary vector related with the solution ---*/

    Solution = new su2double[nVar];
    Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
    Primitive = new su2double[nVar];
    Primitive_i = new su2double[nVar]; Primitive_j = new su2double[nVar];

    /*--- Jacobians and vector structures for implicit computations ---*/

    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (SST model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);

    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

    if (ReducerStrategy)
      EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (multizone){
      Residual_BGS = new su2double[nVar]();
      Residual_Max_BGS = new su2double[nVar]();

      /*--- Define some structures for locating max residuals ---*/

      Point_Max_BGS = new unsigned long[nVar]();
      Point_Max_Coord_BGS = new su2double*[nVar];
      for (iVar = 0; iVar < nVar; iVar++) {
        Point_Max_Coord_BGS[iVar] = new su2double[nDim]();
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
//  constants[8] = 5./9.;  //alfa_1
//  constants[9] = 0.44;  //alfa_2
      
  /*--- Initialize lower and upper limits---*/
  lowerlimit[0] = numeric_limits<passivedouble>::epsilon();
  upperlimit[0] = 1.0e15;

  lowerlimit[1] = numeric_limits<passivedouble>::epsilon();
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

  SlidingState       = new su2double*** [nMarker]();
  SlidingStateNodes  = new int*         [nMarker]();

  for (iMarker = 0; iMarker < nMarker; iMarker++){

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){

      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)]();
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)]();

      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++)
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1]();
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

CTurbSSTSolver::~CTurbSSTSolver(void) {

  unsigned long iMarker, iVertex;
  unsigned short iVar;
  
  if (Primitive != NULL)   delete [] Primitive;
  if (Primitive_i != NULL) delete [] Primitive_i;
  if (Primitive_j != NULL) delete [] Primitive_j;

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


void CTurbSSTSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
         unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  const bool limiter_turb = (config->GetKind_SlopeLimit_Turb() != NO_LIMITER) &&
                            (config->GetInnerIter() <= config->GetLimiterIter());
  
  /*--- Clear residual and system matrix, not needed for
   * reducer strategy as we write over the entire matrix. ---*/
  if (!ReducerStrategy) {
    LinSysRes.SetValZero();
    Jacobian.SetValZero();
  }
  
  /*--- Set flow solver primitives to values stored in turb solver ---*/

//  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();
//
//  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
//    for (unsigned short iVar = 0; iVar < nDim+7; iVar++) {
//      flowNodes->SetPrimitive(iPoint,iVar,nodes->GetFlowPrimitive(iPoint,iVar));
//    }
//  }
  
  /*--- Set primitives and gradients since flow primitives have updated ---*/
  
  Postprocessing(geometry, solver_container, config, iMesh);
    
  if (config->GetReconstructionGradientRequired()) {
    if (config->GetKind_Gradient_Method_Recon() == GREEN_GAUSS)
      SetPrimitive_Gradient_GG(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == LEAST_SQUARES)
      SetPrimitive_Gradient_LS(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES)
      SetPrimitive_Gradient_LS(geometry, config, true);
  }

  if (limiter_turb) SetPrimitive_Limiter(geometry, config);

}

void CTurbSSTSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  
  SetPrimitive_Variables(solver_container);
  
  /*--- Compute mean flow and turbulence gradients ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetPrimitive_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetPrimitive_Gradient_LS(geometry, config);
  
  /*--- Compute eddy viscosity ---*/

//  solver_container[FLOW_SOL]->Preprocessing(geometry, solver_container, config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
  SetEddyViscosity(geometry, solver_container);
  
  /*--- Store variables from the mean flow solver ---*/

//  SetFlowPrimitive(solver_container);
//  SetFlowGradient(solver_container);

}

void CTurbSSTSolver::SetPrimitive_Variables(CSolver **solver_container) {
  
  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    const su2double rho = flowNodes->GetDensity(iPoint);
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      const su2double cons = nodes->GetSolution(iPoint,iVar);
      nodes->SetPrimitive(iPoint,iVar,cons/rho);
    }
  }

}

void CTurbSSTSolver::SetFlowPrimitive(CSolver **solver_container) {
  
  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    for (unsigned short iVar = 0; iVar < nDim+7; iVar++) {
      if (iVar == nDim+2) nodes->SetFlowPrimitive(iPoint,iVar,flowNodes->GetDensity(iPoint));
      else if (iVar == nDim+6) nodes->SetFlowPrimitive(iPoint,iVar,nodes->GetmuT(iPoint));
      else nodes->SetFlowPrimitive(iPoint,iVar,flowNodes->GetPrimitive(iPoint,iVar));
    }
  }
}

void CTurbSSTSolver::SetFlowGradient(CSolver **solver_container) {
  
  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    for (unsigned short iVar = 0; iVar < nDim+1; iVar++) {
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        nodes->SetFlowGradient(iPoint,iVar,iDim,flowNodes->GetGradient_Primitive(iPoint,iVar,iDim));
      }
    }
  }
}

void CTurbSSTSolver::SetEddyViscosity(CGeometry *geometry, CSolver **solver_container) {
  
  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();
  
  const su2double a1 = constants[7];
  
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Compute blending functions and cross diffusion ---*/
    
    const su2double rho = flowNodes->GetDensity(iPoint);
    const su2double mu  = flowNodes->GetLaminarViscosity(iPoint);

    const su2double dist = geometry->node[iPoint]->GetWall_Distance();

    const su2double *Vorticity = flowNodes->GetVorticity(iPoint);
    const su2double VorticityMag = sqrt(Vorticity[0]*Vorticity[0] +
                                        Vorticity[1]*Vorticity[1] +
                                        Vorticity[2]*Vorticity[2]);
    
    nodes->SetVorticityMag(iPoint, VorticityMag);
    
    nodes->SetBlendingFunc(iPoint, mu, dist, rho);

    const su2double F2 = nodes->GetF2blending(iPoint);

    /*--- Compute the eddy viscosity ---*/

    const su2double kine  = nodes->GetPrimitive(iPoint,0);
    const su2double omega = nodes->GetPrimitive(iPoint,1);
    const su2double zeta  = max(omega, VorticityMag*F2/a1);
    const su2double muT   = rho*kine/zeta;

    nodes->SetmuT(iPoint,muT);
        
  }

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

void CTurbSSTSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Loop over all points. ---*/

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    if (geometry->node[iPoint]->GetWall_Distance() > 1.0e-10) {

      /*--- Conservative variables w/o reconstruction ---*/

      numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), nullptr);

      /*--- Gradient of the primitive and conservative variables ---*/

      numerics->SetPrimVarGradient(flowNodes->GetGradient_Primitive(iPoint), nullptr);
//      numerics->SetPrimVarGradient(nodes->GetFlowGradient(iPoint), nullptr);

      /*--- Turbulent variables w/o reconstruction, and its gradient ---*/

      numerics->SetTurbVar(nodes->GetPrimitive(iPoint), nullptr);
      numerics->SetTurbVarGradient(nodes->GetGradient(iPoint), nullptr);

      /*--- Set volume ---*/

      numerics->SetVolume(geometry->node[iPoint]->GetVolume());

      /*--- Set distance to the surface ---*/

      numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);

      /*--- Menter's first blending function ---*/

      numerics->SetF1blending(nodes->GetF1blending(iPoint),0.0);

      /*--- Menter's second blending function ---*/

      numerics->SetF2blending(nodes->GetF2blending(iPoint),0.0);

      /*--- Set vorticity and strain rate magnitude ---*/

      numerics->SetVorticity(flowNodes->GetVorticity(iPoint), nullptr);
      
      numerics->SetVorticityMag(nodes->GetVorticityMag(iPoint), 0.0);

      numerics->SetStrainMag(flowNodes->GetStrainMag(iPoint), 0.0);

      /*--- Cross diffusion ---*/

      numerics->SetCrossDiff(nodes->GetCrossDiff(iPoint),0.0);

      /*--- Compute the source term ---*/

      auto residual = numerics->ComputeResidual(config);

      /*--- Subtract residual and the Jacobian ---*/

      LinSysRes.SubtractBlock(iPoint, residual);
      Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
      
      /*--- Compute Jacobian for gradient terms in cross-diffusion ---*/
      Cross_Diffusion_Jacobian(geometry, solver_container, config, iPoint);
      
    }

  }
  
}

void CTurbSSTSolver::Cross_Diffusion_Jacobian(CGeometry *geometry,
                                              CSolver **solver_container,
                                              CConfig *config,
                                              unsigned long iPoint) {
  
  AD_BEGIN_PASSIVE
  
  const CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();
  const su2double sigma_om2 = constants[3];
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    
    if (geometry->node[iPoint]->GetWall_Distance() > 1.0e-10) {
      const su2double F1_i     = nodes->GetF1blending(iPoint);
      const su2double r_i      = flowNodes->GetPrimitive(iPoint, nDim+2);
      const su2double om_i     = nodes->GetPrimitive(iPoint,1);
      
      Jacobian_i[0][0] = 0.; Jacobian_i[0][1] = 0.;
      Jacobian_i[1][0] = 0.; Jacobian_i[1][1] = 0.;
      Jacobian_j[0][0] = 0.; Jacobian_j[0][1] = 0.;
      Jacobian_j[1][0] = 0.; Jacobian_j[1][1] = 0.;
      
      /*--- Contribution of TurbVar_{i,j} to cross diffusion gradient Jacobian at i ---*/
      for (unsigned short iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
        const unsigned long jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
        const unsigned long iEdge = geometry->FindEdge(iPoint,jPoint);
        const su2double *Normal = geometry->edge[iEdge]->GetNormal();
        const su2double r_j  = flowNodes->GetPrimitive(jPoint, nDim+2);
        const su2double sign = (iPoint < jPoint) ? 1.0 : -1.0;

        Jacobian_i[1][0] = 0.; Jacobian_i[1][1] = 0.;
        Jacobian_j[1][0] = 0.; Jacobian_j[1][1] = 0.;

        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          Jacobian_i[1][0] += sign*(1. - F1_i)*sigma_om2
                            * Normal[iDim]*nodes->GetGradient(iPoint,1,iDim)/(om_i);
          Jacobian_i[1][1] += sign*(1. - F1_i)*sigma_om2
                            * Normal[iDim]*nodes->GetGradient(iPoint,0,iDim)/(om_i);
          Jacobian_j[1][0] += sign*(1. - F1_i)*sigma_om2*r_i
                            * Normal[iDim]*nodes->GetGradient(iPoint,1,iDim)/(r_j*om_i);
          Jacobian_j[1][1] += sign*(1. - F1_i)*sigma_om2*r_i
                            * Normal[iDim]*nodes->GetGradient(iPoint,0,iDim)/(r_j*om_i);
        }
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
        Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
      }
      
      /*--- Boundary contribution to cross diffusion gradient Jacobian at i ---*/
      if (geometry->node[iPoint]->GetPhysicalBoundary()) {
        Jacobian_i[1][0] = 0.; Jacobian_i[1][1] = 0.;
        for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
          const long iVertex = geometry->node[iPoint]->GetVertex(iMarker);
          if (iVertex > -1) {
            const su2double *Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            for (unsigned short iDim = 0; iDim < nDim; iDim++) {
              Jacobian_i[1][0] -= 2.*(1. - F1_i)*sigma_om2
                                * Normal[iDim]*nodes->GetGradient(iPoint,1,iDim)/(om_i);
              Jacobian_i[1][1] -= 2.*(1. - F1_i)*sigma_om2
                                * Normal[iDim]*nodes->GetGradient(iPoint,0,iDim)/(om_i);
            }
          }
        }
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
      }
    }
  }
  
  AD_END_PASSIVE
}

void CTurbSSTSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {
}

void CTurbSSTSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iVar;
  unsigned short iDim;
  su2double distance, density_s = 0.0, density_v = 0.0, laminar_viscosity_v = 0.0, beta_1 = constants[4];
  su2double energy_v = 0.0, vel2_v = 0.0, staticenergy_v, k_v;
  
  const CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();
  CFluidModel *fluidModel = solver_container[FLOW_SOL]->GetFluidModel();
  
  /*--- The dirichlet condition is used only without wall function, otherwise the
  convergence is compromised ---*/
  
//  if (!config->GetWall_Functions()) {

    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
      iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

      /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
      if (geometry->node[iPoint]->GetDomain()) {
        
        jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
        
        distance = geometry->node[jPoint]->GetWall_Distance();
        
        density_s = flowNodes->GetDensity(iPoint);
        density_v = flowNodes->GetDensity(jPoint);
        
        energy_v = flowNodes->GetEnergy(jPoint);
        vel2_v = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) vel2_v += pow(flowNodes->GetSolution(jPoint,iDim+1)/density_v, 2.);
        k_v = nodes->GetSolution(jPoint,0)/density_v;
        staticenergy_v = energy_v - 0.5*vel2_v - k_v;

        fluidModel->SetTDState_rhoe(density_v, staticenergy_v);
        laminar_viscosity_v = fluidModel->GetLaminarViscosity();
        
        Solution[0] = 0.0;
        Solution[1] = 60.0*density_s*laminar_viscosity_v/(density_v*beta_1*distance*distance);

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
//  }
}

void CTurbSSTSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Call the equivalent heat flux wall boundary condition. ---*/
  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}

void CTurbSSTSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                  CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex, Point_Normal;
  su2double *Normal, *V_infty, *V_domain;
  const su2double *Vel_Infty = config->GetVelocity_FreeStreamND();
  su2double Vn_Infty = 0., Velocity2 = 0.;
  const su2double Intensity = config->GetTurbulenceIntensity_FreeStream();
  su2double Kine_Infty, Omega_Infty;
  unsigned short iVar, iDim;

  Normal = new su2double[nDim];

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Allocate the value at the infinity ---*/

      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Set turbulent variable at the wall, and at infinity ---*/

      for (iVar = 0; iVar < nVar; iVar++) Primitive_i[iVar] = nodes->GetPrimitive(iPoint,iVar);

      /*--- Set Normal (it is necessary to change the sign) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Set primitive state based on flow direction ---*/
      
      Vn_Infty = 0;
      for (iDim = 0; iDim < nDim; iDim++) Vn_Infty += Vel_Infty[iDim+1]*Normal[iDim];

      if (Vn_Infty > 0.0) {
        /*--- Outflow conditions ---*/
        Primitive_j[0] = Primitive_i[0];
        Primitive_j[1] = Primitive_i[1];
      }
      else {
        /*--- Inflow conditions ---*/
        Velocity2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) Velocity2 += pow(V_infty[iDim+1],2.);
        const su2double Rho_Infty = V_infty[nDim+2];
        const su2double muT_Infty = V_infty[nDim+6];
        Kine_Infty  = 3.0/2.0*(Velocity2*Intensity*Intensity);
        Omega_Infty = Rho_Infty*Kine_Infty/muT_Infty;

        Primitive_j[0] = Kine_Infty;
        Primitive_j[1] = Omega_Infty;
      }
      
      conv_numerics->SetTurbVar(Primitive_i, Primitive_j);

      /*--- Grid Movement ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute residuals and Jacobians ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Add residuals and Jacobians ---*/

      LinSysRes.AddBlock(iPoint, residual);
      Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      /*--- Viscous contribution ---*/
      
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                              geometry->node[iPoint]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_domain);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Primitive_i, Primitive_i);
      visc_numerics->SetTurbVarGradient(nodes->GetGradient(iPoint),
                                        nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint),
                                   nodes->GetF1blending(iPoint));
      
      /*--- Menter's second blending function ---*/
      visc_numerics->SetF2blending(nodes->GetF2blending(iPoint),
                                   nodes->GetF2blending(iPoint));
      
      /*--- Vorticity ---*/
      visc_numerics->SetVorticity(solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint),
                                  solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint));
      
      visc_numerics->SetVorticityMag(nodes->GetVorticityMag(iPoint),
                                     nodes->GetVorticityMag(iPoint));

      /*--- Set values for gradient Jacobian ---*/
      visc_numerics->SetVolume(geometry->node[iPoint]->GetVolume(),
                               geometry->node[iPoint]->GetVolume());

      /*--- Compute residual, and Jacobians ---*/
      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, visc_residual);
      Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);
      Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_j);
      
      /*--- Compute Jacobian correction for influence from all neighbors ---*/
      CorrectJacobian(geometry, solver_container, config, iPoint, iPoint, visc_residual.jacobian_ic, nullptr);
      CorrectJacobian(geometry, solver_container, config, iPoint, iPoint, visc_residual.jacobian_jc, nullptr);
        
    }
  }

  delete [] Normal;

}

void CTurbSSTSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                              CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

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

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      /*--- Viscous contribution, commented out because serious convergence problems ---*/

      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/

      visc_numerics->SetPrimitive(V_domain, V_domain);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/

      visc_numerics->SetTurbVar(Primitive_i, Primitive_i);
      visc_numerics->SetTurbVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint),
                                   nodes->GetF1blending(iPoint));
      
      /*--- Menter's second blending function ---*/
      visc_numerics->SetF2blending(nodes->GetF2blending(iPoint),
                                   nodes->GetF2blending(iPoint));
      
      /*--- Vorticity ---*/
      visc_numerics->SetVorticity(solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint),
                                  solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint));
      
      visc_numerics->SetVorticityMag(nodes->GetVorticityMag(iPoint),
                                     nodes->GetVorticityMag(iPoint));

      /*--- Set values for gradient Jacobian ---*/
      visc_numerics->SetVolume(geometry->node[iPoint]->GetVolume(),
                               geometry->node[iPoint]->GetVolume());

      /*--- Compute residual, and Jacobians ---*/
      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, visc_residual);
      Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);
      Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_j);
      
      /*--- Compute Jacobian correction for influence from all neighbors ---*/
      CorrectJacobian(geometry, solver_container, config, iPoint, iPoint, visc_residual.jacobian_ic, nullptr);
      CorrectJacobian(geometry, solver_container, config, iPoint, iPoint, visc_residual.jacobian_jc, nullptr);

    }

  }

  /*--- Free locally allocated memory ---*/

  delete [] Normal;

}

void CTurbSSTSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
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
//      visc_numerics->SetPrimitive(V_domain, V_domain);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetTurbVar(Primitive_i, Primitive_i);
//      visc_numerics->SetTurbVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));
//
//      /*--- Menter's first blending function ---*/
//      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint),
//                                   nodes->GetF1blending(iPoint));
//
//      /*--- Menter's second blending function ---*/
//      visc_numerics->SetF2blending(nodes->GetF2blending(iPoint),
//                                   nodes->GetF2blending(iPoint));
//
//      /*--- Vorticity ---*/
//      visc_numerics->SetVorticity(solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint),
//                                  solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint));
//
//      visc_numerics->SetVorticityMag(nodes->GetVorticityMag(iPoint),
//                                     nodes->GetVorticityMag(iPoint));
//
//      /*--- Set values for gradient Jacobian ---*/
//      visc_numerics->SetVolume(geometry->node[iPoint]->GetVolume(),
//                               geometry->node[iPoint]->GetVolume());
//
//      /*--- Compute residual, and Jacobians ---*/
//      auto visc_residual = visc_numerics->ComputeResidual(config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//      LinSysRes.SubtractBlock(iPoint, visc_residual);
//      Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);
//      Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_j);
//
//      /*--- Compute Jacobian correction for influence from all neighbors ---*/
//      CorrectJacobian(geometry, solver_container, config, iPoint, iPoint, visc_residual.jacobian_ic, nullptr);
//      CorrectJacobian(geometry, solver_container, config, iPoint, iPoint, visc_residual.jacobian_jc, nullptr);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}


void CTurbSSTSolver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                          CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

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
      visc_numerics->SetTurbVar(Primitive_i, Primitive_j);
      visc_numerics->SetTurbVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(iPoint));

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

void CTurbSSTSolver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

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
      visc_numerics->SetTurbVar(Primitive_i, Primitive_j);
      visc_numerics->SetTurbVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(iPoint));

      /*--- Compute residual, and Jacobians ---*/
      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, visc_residual);
      Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  delete[] Vel;

}

void CTurbSSTSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                        CNumerics *visc_numerics, CConfig *config){

  unsigned long iVertex, jVertex, iPoint, Point_Normal = 0;
  unsigned short iDim, iVar, jVar, iMarker;

  unsigned short nPrimVar = solver_container[FLOW_SOL]->GetnPrimVar();
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];

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

          for (iVar = 0; iVar < nVar; iVar++) {
            Residual[iVar] = 0.0;
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          }

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

            auto residual = conv_numerics->ComputeResidual(config);

            /*--- Accumulate the residuals to compute the average ---*/

            for (iVar = 0; iVar < nVar; iVar++) {
              Residual[iVar] += weight*residual.residual[iVar];
              for (jVar = 0; jVar < nVar; jVar++)
                Jacobian_i[iVar][jVar] += weight*residual.jacobian_i[iVar][jVar];
            }
          }

          /*--- Add Residuals and Jacobians ---*/

          LinSysRes.AddBlock(iPoint, Residual);

          Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

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

          auto residual = visc_numerics->ComputeResidual(config);

          LinSysRes.SubtractBlock(iPoint, residual);

          /*--- Jacobian contribution for implicit integration ---*/

          Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

        }
      }
    }
  }

  /*--- Free locally allocated memory ---*/

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

//void CTurbSSTSolver::ComputeWallFunction(CGeometry *geometry, CSolver **solver, CConfig *config) {
//
//  unsigned long jPoint, total_index;
//  unsigned short kNode;
//  long iElem;
//  su2double distance, density = 0.0, laminar_viscosity = 0.0, k = 0.0, beta_1 = constants[4];
//  su2double *weights;
//
//  /*--- Communicate values needed for WF ---*/
//  WallFunctionComms(geometry, solver, config);
//
//  for (jPoint = 0; jPoint < nPointDomain; jPoint++) {
//    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
//    if (geometry->node[jPoint]->GetBool_Wall_Neighbor()) {
//
//      iElem = geometry->node[jPoint]->GetWall_Element();
//
//      density = solver[FLOW_SOL]->GetNodes()->GetDensity(jPoint);
//      laminar_viscosity = solver[FLOW_SOL]->GetNodes()->GetLaminarViscosity(jPoint);
//      k = nodes->GetPrimitive(jPoint, 0);
//
//      weights = geometry->node[jPoint]->GetWall_Interpolation_Weights();
//
//      distance = geometry->node[jPoint]->GetWall_Distance();
//
//      su2double Omega_0 = sqrt(k) / (pow(0.09,0.25) * 0.41 * distance);
//      su2double Omega = 0.0;
//
//      unsigned short nWall = geometry->node[jPoint]->GetWall_nNode();
//
//      for (kNode = 0; kNode < nWall; kNode++) {
//
//        const su2double DensityWall = nodes->GetWallDensity(jPoint, kNode);
//        const su2double LamViscWall = nodes->GetWallLamVisc(jPoint, kNode);
//
//        if (DensityWall != DensityWall) cout << "Nan detected at rank " << rank << ", node "<< jPoint << endl;
//        if (DensityWall < 1.0E-16) cout << "Small value detected at " << rank << ", node " << jPoint << endl;
//
//        const su2double Omega_i = 6. * LamViscWall / (beta_1 * DensityWall * pow(distance, 2.0));
//        Omega += sqrt(pow(Omega_0, 2.) + pow(Omega_i, 2.))*weights[kNode];
//      }
//
//      nodes->SetSolution_Old(jPoint,1,density*Omega);
//      nodes->SetSolution(jPoint,1,density*Omega);
//      LinSysRes.SetBlock_Zero(jPoint,1);
//
//      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
//      total_index = jPoint*nVar+1;
//      Jacobian.DeleteValsRowi(total_index);
//    }
//  }
//
//}

void CTurbSSTSolver::ComputeWallFunction(CGeometry *geometry, CSolver **solver, CConfig *config) {

  /*--- Local variables ---*/
  
  CVariable *flowNodes = solver[FLOW_SOL]->GetNodes();

  unsigned short iDim, jDim, iNode, iVar, iMarker;
  unsigned long iVertex, iPoint, jPoint, counter;

  su2double Area;
  su2double div_vel, Normal[3], UnitNormal[3];
  su2double **grad_primvar, tau[3][3];
  su2double Vel[3], VelNormal, VelTang[3], VelTangMod, VelInfMod, WallDist[3], WallDistMod;
  su2double Lam_Visc_Normal, dypw_dyp, Eddy_Visc;
  su2double T_Normal, P_Normal, Density_Normal;
  su2double Density_Wall, T_Wall, P_Wall, Lam_Visc_Wall, Tau_Wall, Tau_Wall_Old;
  su2double *Coord, *Coord_Normal;
  su2double diff, Delta;
  su2double U_Tau, U_Plus = 0.0, Gam = 0.0, Beta = 0.0, Phi, Q = 0.0, Y_Plus_White = 0.0, Y_Plus;
  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double k, Omega, Omega_vis, Omega_log;
  
  const su2double beta_1 = constants[4];

  unsigned short max_iter = 100;
  su2double tol = 1e-10;

  /*--- Compute the recovery factor ---*/
  // su2double-check: laminar or turbulent Pr for this?
  su2double Recovery = pow(config->GetPrandtl_Lam(),(1.0/3.0));

  /*--- Typical constants from boundary layer theory ---*/

  su2double kappa = 0.4;
  su2double B = 5.5;
  
  for (jPoint = 0; jPoint < nPointDomain; jPoint++) {

    if (geometry->node[jPoint]->GetBool_Wall_Neighbor()) {

      Omega = 0.;
      Eddy_Visc = 0.;
      
      for(iNode = 0; iNode < geometry->node[jPoint]->GetWall_nNode(); iNode++) {

        /*--- Compute dual-grid area and boundary normal ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          Normal[iDim] = flowNodes->GetWallNormal(jPoint, iNode, iDim);

        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Area += Normal[iDim]*Normal[iDim];
        Area = sqrt (Area);

        for (iDim = 0; iDim < nDim; iDim++)
          UnitNormal[iDim] = -Normal[iDim]/Area;

        /*--- Get the velocity, pressure, and temperature at the nearest
         (normal) interior point. ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          Vel[iDim]    = flowNodes->GetVelocity(jPoint,iDim);
        P_Normal       = flowNodes->GetPressure(jPoint);
        T_Normal       = flowNodes->GetTemperature(jPoint);

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
        
        WallDistMod = geometry->node[jPoint]->GetWall_Distance();

        /*--- Compute the wall temperature using the Crocco-Buseman equation ---*/

//        T_Wall = T_Normal + Recovery*pow(VelTangMod,2.0)/(2.0*Cp);
        T_Wall = T_Normal/(1.+Recovery*Gamma_Minus_One/2.*pow(VelTangMod,2.));

        /*--- Extrapolate the pressure from the interior & compute the
         wall density using the equation of state ---*/

        P_Wall = P_Normal;
        Density_Wall = P_Wall/(Gas_Constant*T_Wall);

        /*--- Compute the shear stress at the wall in the regular fashion
         by using the stress tensor on the surface ---*/

        Lam_Visc_Wall = flowNodes->GetWallLamVisc(jPoint, iNode);
        
        Tau_Wall_Old = flowNodes->GetWallTau(jPoint, iNode);
        
        U_Tau = sqrt(Tau_Wall_Old/Density_Wall);
        U_Plus = VelTangMod/U_Tau;
        Gam  = Recovery*U_Tau*U_Tau/(2.0*Cp*T_Wall);
        Beta = 0.0; // For adiabatic flows only
        Q    = sqrt(Beta*Beta + 4.0*Gam);
        Phi  = asin(-1.0*Beta/Q);
        
        /*--- Y+ defined by White & Christoph (compressibility and heat transfer) ---*/

        Y_Plus_White = exp((kappa/sqrt(Gam))*(asin((2.0*Gam*U_Plus - Beta)/Q) - Phi))*exp(-1.0*kappa*B);

        Lam_Visc_Normal = flowNodes->GetLaminarViscosity(jPoint);
        
        const su2double w = geometry->node[jPoint]->GetWall_Interpolation_Weights()[iNode];

        dypw_dyp = 2.0*Y_Plus_White*(kappa*sqrt(Gam)/Q)*sqrt(1.0 - pow(2.0*Gam*U_Plus - Beta,2.0)/(Q*Q));
        Eddy_Visc += w*max(Lam_Visc_Wall*(1.0 + dypw_dyp - kappa*exp(-1.0*kappa*B)*
                          (1.0 + kappa*U_Plus + kappa*kappa*U_Plus*U_Plus/2.0)
                         - Lam_Visc_Normal/Lam_Visc_Wall), 0.0);
        
        Omega_vis = 6. * Lam_Visc_Wall/ (beta_1 * Density_Wall* pow(WallDistMod, 2.0));
        Omega_log = U_Tau/( 0.3 * 0.41 * WallDistMod);
        
        Omega += w*sqrt(pow(Omega_vis, 2.) + pow(Omega_log, 2.));
        
      }
      
      /*--- Eddy viscosity should be always a positive number ---*/
      
      Density_Normal  = flowNodes->GetDensity(jPoint);
      
      /*--- Now compute the TKE and dissipation at the first point off of the wall ---*/
      
      
      k = Eddy_Visc*Omega/Density_Normal;

      Solution[0] = Density_Normal*k;
      Solution[1] = Density_Normal*Omega;

      /*--- Set the solution values and zero the residual ---*/
      nodes->SetSolution_Old(jPoint,Solution);
      nodes->SetSolution(jPoint,Solution);
      LinSysRes.SetBlock_Zero(jPoint);

      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        const unsigned long total_index = jPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }

    }
  }
}

void CTurbSSTSolver::TurbulentMetric(CSolver                    **solver,
                                     CGeometry                  *geometry,
                                     CConfig                    *config,
                                     unsigned long              iPoint,
                                     vector<vector<su2double> > &weights) {
  CVariable *varFlo    = solver[FLOW_SOL]->GetNodes(),
            *varTur    = solver[TURB_SOL]->GetNodes(),
            *varAdjFlo = solver[ADJFLOW_SOL]->GetNodes(),
            *varAdjTur = solver[ADJTURB_SOL]->GetNodes();

  unsigned short iDim, jDim, iVar;
  const unsigned short nVarFlo = solver[FLOW_SOL]->GetnVar();
  const unsigned short nVarTur = solver[TURB_SOL]->GetnVar();
  
  const su2double eps = numeric_limits<passivedouble>::epsilon();

  //--- First-order terms (error due to viscosity)
  su2double r, u[3], k, omega,
            mu, mut, lam, lamt,
            R, cv, cp, g, Pr, Prt,
            walldist;

  r = varFlo->GetDensity(iPoint);
  u[0] = varFlo->GetVelocity(iPoint, 0);
  u[1] = varFlo->GetVelocity(iPoint, 1);
  if (nDim == 3) u[2] = varFlo->GetVelocity(iPoint, 2);
  k = varTur->GetPrimitive(iPoint, 0);
  omega = varTur->GetPrimitive(iPoint, 1);
  
  mu  = varFlo->GetLaminarViscosity(iPoint);
  mut = nodes->GetmuT(iPoint);
//  mut = r*k/omega;

  g    = config->GetGamma();
  R    = config->GetGas_ConstantND();
  cp   = (g/(g-1.))*R;
  cv   = cp/g;
  Pr   = config->GetPrandtl_Lam();
  Prt  = config->GetPrandtl_Turb();
  lam  = cp*mu/Pr;
  lamt = cp*mut/Prt;
  
  walldist = geometry->node[iPoint]->GetWall_Distance();

  su2double gradu[3][3], gradT[3], gradk[3], gradomega[3], divu, taut[3][3],
            delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}},
            pk = 0.;

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
      pk += taut[iDim][jDim]*gradu[iDim][jDim];
    }
  }

  const su2double F1   = varTur->GetF1blending(iPoint);
  const su2double CDkw = varTur->GetCrossDiff(iPoint);

  const su2double alfa        = F1*constants[8] + (1.0 - F1)*constants[9];
  const su2double sigmak      = F1*constants[0] + (1.0 - F1)*constants[1];
  const su2double sigmaomega  = F1*constants[2] + (1.0 - F1)*constants[3];
  const su2double sigmaomega2 = constants[3];
  const su2double beta        = F1*constants[4] + (1.0 - F1)*constants[5];
  const su2double betastar    = constants[6];
  const su2double a1          = constants[7];

  const su2double* Vorticity = varFlo->GetVorticity(iPoint);
  const su2double VorticityMag = sqrt(Vorticity[0]*Vorticity[0] +
                                      Vorticity[1]*Vorticity[1] +
                                      Vorticity[2]*Vorticity[2]);
  
//  const su2double lim = (omega > VorticityMag*F1/a1) ? 1.0 : 0.0;
//  const su2double zeta = min(1./omega,a1/(VorticityMag*F1));
  
  const su2double lim = 1.0;
  const su2double zeta = 1./omega;

  //--- Momentum weights
  vector<su2double> TmpWeights(weights[0].size(), 0.0);
  su2double factor = 0.0;
  if (pk > 0.) {
    for (iDim = 0; iDim < nDim; ++iDim) {
      factor = 0.0;
//      if (pk <= 20.*betastar*r*omega*k) {
        factor += -(2./3.)*divu*alfa*varAdjTur->GetGradient_Adaptation(iPoint, 1, iDim);
        factor += -(2./3.)*divu*mut/r*varAdjTur->GetGradient_Adaptation(iPoint, 0, iDim);
//      }
      for (jDim = 0; jDim < nDim; ++jDim) {
//        if (pk <= 20.*betastar*r*omega*k) {
          factor += (taut[iDim][jDim]+mut*(gradu[iDim][jDim]+gradu[jDim][iDim]))*alfa/max(mut,eps)*varAdjTur->GetGradient_Adaptation(iPoint, 1, jDim);
          factor += (taut[iDim][jDim]+mut*(gradu[iDim][jDim]+gradu[jDim][iDim]))/r*varAdjTur->GetGradient_Adaptation(iPoint, 0, jDim);
//        }
      }
      TmpWeights[iDim+1] += factor;
    }
  }

  //--- k and omega weights
  factor = 0.0;
  for (iDim = 0; iDim < nDim; ++iDim) {
    for (jDim = 0; jDim < nDim; ++jDim) {
      iVar = iDim+1;
      factor += (taut[iDim][jDim]+(2./3.)*r*k*delta[iDim][jDim])/max(mut,eps)
              * (varAdjFlo->GetGradient_Adaptation(iPoint, iVar, jDim)
              + u[jDim]*varAdjFlo->GetGradient_Adaptation(iPoint, (nVarFlo-1), iDim));
    }
    factor += cp/Prt*gradT[iDim]*varAdjFlo->GetGradient_Adaptation(iPoint, (nVarFlo-1), iDim);
    factor += sigmak*gradk[iDim]*varAdjTur->GetGradient_Adaptation(iPoint, 0, iDim)
            + sigmaomega*gradomega[iDim]*varAdjTur->GetGradient_Adaptation(iPoint, 1, iDim);
  }

  TmpWeights[nVarFlo+0] += zeta*factor;
  TmpWeights[nVarFlo+1] += -lim*k/pow(omega,2.)*factor;
  for (iDim = 0; iDim < nDim; ++iDim) {
    TmpWeights[nVarFlo+0] += 2.*(1.-F1)*sigmaomega2/omega*gradomega[iDim]*varAdjTur->GetGradient_Adaptation(iPoint, 1, iDim);
    TmpWeights[nVarFlo+1] += 2.*(1.-F1)*sigmaomega2/omega*gradk[iDim]*varAdjTur->GetGradient_Adaptation(iPoint, 1, iDim);
  }

  //--- Density weight
  for (iDim = 0; iDim < nDim; ++iDim) TmpWeights[0] += -u[iDim]*TmpWeights[iDim+1];
  TmpWeights[0] += -k*TmpWeights[nVarFlo+0] - omega*TmpWeights[nVarFlo+1]
                 + lim*k/omega*factor;

  //--- Add TmpWeights to weights, then reset for second-order terms
  for (iVar = 0; iVar < nVarFlo+nVarTur; ++iVar) weights[1][iVar] += TmpWeights[iVar];
  fill(TmpWeights.begin(), TmpWeights.end(), 0.0);

  //--- Second-order terms (error due to gradients)
  if(nDim == 3) {
    const unsigned short rki = 0, romegai = 1, rei = (nVarFlo - 1), xxi = 0, yyi = 3, zzi = 5;
    TmpWeights[nVarFlo+0] += -(mu+sigmak*mut)/r*(varAdjTur->GetHessian(iPoint, rki, xxi)
                                                +varAdjTur->GetHessian(iPoint, rki, yyi)
                                                +varAdjTur->GetHessian(iPoint, rki, zzi))
                             +(lam+lamt)/(r*cv)*(varAdjFlo->GetHessian(iPoint, rei, xxi)
                                                 +varAdjFlo->GetHessian(iPoint, rei, yyi)
                                                 +varAdjFlo->GetHessian(iPoint, rei, zzi)); // Hk
    TmpWeights[nVarFlo+1] += -(mu+sigmaomega*mut)/r*(varAdjTur->GetHessian(iPoint, romegai, xxi)
                                                    +varAdjTur->GetHessian(iPoint, romegai, yyi)
                                                    +varAdjTur->GetHessian(iPoint, romegai, zzi)); // Homega

  }
  else {
    const unsigned short rki = 0, romegai = 1, rei = (nVarFlo - 1), xxi = 0, yyi = 2;
    TmpWeights[nVarFlo+0] += -(mu+sigmak*mut)/r*(varAdjTur->GetHessian(iPoint, rki, xxi)
                                                +varAdjTur->GetHessian(iPoint, rki, yyi))
                             +(lam+lamt)/(r*cv)*(varAdjFlo->GetHessian(iPoint, rei, xxi)
                                                 +varAdjFlo->GetHessian(iPoint, rei, yyi)); // Hk
    TmpWeights[nVarFlo+1] += -(mu+sigmaomega*mut)/r*(varAdjTur->GetHessian(iPoint, romegai, xxi)
                                                    +varAdjTur->GetHessian(iPoint, romegai, yyi)); // Homega
  }
  TmpWeights[0] += -k*TmpWeights[nVarFlo+0]-omega*TmpWeights[nVarFlo+1];

  //--- Add TmpWeights to weights
  weights[2][0]         += TmpWeights[0];
  weights[2][nVarFlo+0] += TmpWeights[nVarFlo+0];
  weights[2][nVarFlo+1] += TmpWeights[nVarFlo+1];

  //--- Zeroth-order terms due to production
  if (pk > 0.) {
//    if (pk <= 20.*betastar*r*omega*k){
      weights[0][nVarFlo+0] += (2./3.)*divu*varAdjTur->GetSolution(iPoint,0);
      weights[0][nVarFlo+1] += (2./3.)*lim*alfa*divu*varAdjTur->GetSolution(iPoint,1);
//    }
//    else {
//      weights[0][0]         += 20.0*betastar*k*omega*varAdjTur->GetSolution(iPoint,0)
//                             + 20.0*lim*alfa*betastar*pow(omega,2.)*varAdjTur->GetSolution(iPoint,1);;
//      weights[0][nVarFlo+0] += -20.0*betastar*omega*varAdjTur->GetSolution(iPoint,0);
//      weights[0][nVarFlo+1] += -20.0*betastar*k*varAdjTur->GetSolution(iPoint,0)
//                             - 40.0*lim*alfa*betastar*omega*varAdjTur->GetSolution(iPoint,1)
//                             - 20.0*(1.-lim)*alfa*betastar*VorticityMag*F1/a1*varAdjTur->GetSolution(iPoint,1);
//    }
  }
  
  //--- Zeroth-order terms due to dissipation
  weights[0][0]         += -betastar*k*omega*varAdjTur->GetSolution(iPoint,0)
                         - beta*pow(omega,2.)*varAdjTur->GetSolution(iPoint,1);
  weights[0][nVarFlo+0] += betastar*omega*varAdjTur->GetSolution(iPoint,0);
  weights[0][nVarFlo+1] += betastar*k*varAdjTur->GetSolution(iPoint,0)
                         + 2.*beta*omega*varAdjTur->GetSolution(iPoint,1);
  
  //--- Zeroth-order terms due to cross-diffusion
  weights[0][nVarFlo+1] += (1. - F1)*CDkw/(r*omega)*varAdjTur->GetSolution(iPoint,1);

}
