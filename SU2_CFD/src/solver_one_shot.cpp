/*!
 * \file solver_one_shot.cpp
 * \brief Main subroutines for solving the one-shot problem.
 * \author B. Munguía
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../include/solver_structure.hpp"
#include "../include/variables/CDiscAdjVariable.hpp"

COneShotSolver::COneShotSolver(void) : CDiscAdjSolver () {

}

COneShotSolver::COneShotSolver(CGeometry *geometry, CConfig *config)  : CDiscAdjSolver(geometry, config) {

}

COneShotSolver::COneShotSolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver, unsigned short Kind_Solver, unsigned short iMesh)  : CDiscAdjSolver(geometry, config, direct_solver, Kind_Solver, iMesh) {

  rho       = 0.0;
  rho_old   = 1.0;
  theta     = 0.0;
  theta_old = 1.0;
  nConstr   = config->GetnConstr();
  nActiveDV = 0;
  grad_norm = 0.0;

  DConsVec = new su2double** [nConstr];
  for (unsigned short iConstr=0; iConstr<nConstr;iConstr++){
    DConsVec[iConstr] = new su2double* [nPointDomain];
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++){
      DConsVec[iConstr][iPoint] = new su2double [nVar];
      for (unsigned short iVar = 0; iVar < nVar; iVar++){
        DConsVec[iConstr][iPoint][iVar]=0.0;
      }
    }
  }

  geometry->InitializeSensitivity();

}

COneShotSolver::~COneShotSolver(void) {
  for (unsigned short iConstr=0; iConstr < nConstr; iConstr++){
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++){
      delete [] DConsVec[iConstr][iPoint];
    }
    delete [] DConsVec[iConstr];
  }
  delete [] DConsVec;
}

void COneShotSolver::SetRecording(CGeometry* geometry, CConfig *config){


  // bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
  //     (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)),
  // time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

  // unsigned long iPoint;
  // unsigned short iVar;

  // /*--- For the one-shot solver the solution is not reset in each iteration step to the initial solution ---*/

  // if (time_n_needed) {
  //   for (iPoint = 0; iPoint < nPoint; iPoint++) {
  //     for (iVar = 0; iVar < nVar; iVar++) {
  //       AD::ResetInput(direct_solver->GetNodes()->GetSolution_time_n(iPoint)[iVar]);
  //     }
  //   }
  // }
  // if (time_n1_needed) {
  //   for (iPoint = 0; iPoint < nPoint; iPoint++) {
  //     for (iVar = 0; iVar < nVar; iVar++) {
  //       AD::ResetInput(direct_solver->GetNodes()->GetSolution_time_n1(iPoint)[iVar]);
  //     }
  //   }
  // }

  /*--- Set the Jacobian to zero since this is not done inside the fluid iteration
   * when running the discrete adjoint solver. ---*/

  direct_solver->Jacobian.SetValZero();

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}

void COneShotSolver::SetGeometrySensitivityLagrangian(CGeometry *geometry){

    unsigned short iDim;
    unsigned long iPoint;

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        geometry->SetSensitivity(iPoint, iDim, nodes->GetSensitivity_AugmentedLagrangian(iPoint,iDim));
      }
    }
}

void COneShotSolver::SetGeometrySensitivityGradient(CGeometry *geometry){

    unsigned short iDim;
    unsigned long iPoint;

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        geometry->SetSensitivity(iPoint, iDim, nodes->GetSensitivity_ShiftedLagrangian(iPoint,iDim));
      }
    }
}

void COneShotSolver::SetSensitivityShiftedLagrangian(CGeometry *geometry){
    unsigned short iDim;
    unsigned long iPoint;

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        nodes->SetSensitivity_ShiftedLagrangian(iPoint, iDim, nodes->GetSensitivity(iPoint,iDim));
      }
    }
}

void COneShotSolver::ResetSensitivityLagrangian(CGeometry *geometry){
    unsigned short iDim;
    unsigned long iPoint;

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        nodes->SetSensitivity_AugmentedLagrangian(iPoint, iDim, 0.0);
      }
    }
}

void COneShotSolver::UpdateSensitivityLagrangian(CGeometry *geometry, su2double factor){
    unsigned short iDim;
    unsigned long iPoint;

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        nodes->SetSensitivity_AugmentedLagrangian(iPoint,iDim, nodes->GetSensitivity_AugmentedLagrangian(iPoint,iDim)+factor*nodes->GetSensitivity(iPoint,iDim));
      }
    }
}

void COneShotSolver::SetMeshPointsOld(CConfig *config, CGeometry *geometry){
    unsigned long iVertex, jPoint;
    unsigned short iMarker;
    for (jPoint=0; jPoint < nPoint; jPoint++){
        geometry->node[jPoint]->SetCoord_Old(geometry->node[jPoint]->GetCoord());
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          geometry->vertex[iMarker][iVertex]->SetNormal_Old(geometry->vertex[iMarker][iVertex]->GetNormal());
        }
    }
}

void COneShotSolver::LoadMeshPointsOld(CConfig *config, CGeometry *geometry){
    unsigned long iVertex, jPoint;
    unsigned short iMarker;
    for (jPoint=0; jPoint < nPoint; jPoint++){
        geometry->node[jPoint]->SetCoord(geometry->node[jPoint]->GetCoord_Old());
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          geometry->vertex[iMarker][iVertex]->SetNormal(geometry->vertex[iMarker][iVertex]->GetNormal_Old());
        }
    }
}

void COneShotSolver::SetStoreSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution_Store(iPoint);
    nodes->SetSolution_Store(iPoint);
  }
}

void COneShotSolver::LoadSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution(iPoint,direct_solver->GetNodes()->GetSolution_Store(iPoint));
    nodes->SetSolution(iPoint,nodes->GetSolution_Store(iPoint));
  }
}

void COneShotSolver::SetSaveSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution_Save(iPoint);
    nodes->SetSolution_Save(iPoint);
  }
}

void COneShotSolver::LoadSaveSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution(iPoint,direct_solver->GetNodes()->GetSolution_Save(iPoint));
    nodes->SetSolution(iPoint,nodes->GetSolution_Save(iPoint));
  }
}

void COneShotSolver::LoadStepSolution(su2double stepsize){
  unsigned long iPoint;
  unsigned short iVar;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for(iVar = 0; iVar < nVar; iVar++){
      const su2double dy    = direct_solver->GetNodes()->GetSolution_Save(iPoint,iVar) - direct_solver->GetNodes()->GetSolution_Store(iPoint,iVar);
      const su2double dbary = nodes->GetSolution_Save(iPoint,iVar) - nodes->GetSolution_Store(iPoint,iVar);
      direct_solver->GetNodes()->SetSolution(iPoint,iVar,direct_solver->GetNodes()->GetSolution_Store(iPoint,iVar) + stepsize*dy);
      nodes->SetSolution(iPoint,iVar,nodes->GetSolution_Store(iPoint,iVar) + stepsize*dbary);
    }
  }
}

void COneShotSolver::CalculateRhoTheta(CConfig *config){
  unsigned short iVar;
  unsigned long iPoint;
  su2double normDelta=0.0,    myNormDelta=0.0;
  su2double normDeltaNew=0.0, myNormDeltaNew=0.0;
  su2double helper=0.0,       myHelper=0.0;

  /* --- Estimate rho and theta values --- */
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      myNormDelta    += direct_solver->GetNodes()->GetSolution_DeltaStore(iPoint,iVar)
                      * direct_solver->GetNodes()->GetSolution_DeltaStore(iPoint,iVar);
      myNormDeltaNew += (direct_solver->GetNodes()->GetSolution(iPoint,iVar)-direct_solver->GetNodes()->GetSolution_Store(iPoint,iVar))
                      * (direct_solver->GetNodes()->GetSolution(iPoint,iVar)-direct_solver->GetNodes()->GetSolution_Store(iPoint,iVar));
      myHelper       += direct_solver->GetNodes()->GetSolution_DeltaStore(iPoint,iVar)*(nodes->GetSolution(iPoint,iVar)-nodes->GetSolution_Store(iPoint,iVar))
                      - nodes->GetSolution_DeltaStore(iPoint,iVar)*(direct_solver->GetNodes()->GetSolution(iPoint,iVar)-direct_solver->GetNodes()->GetSolution_Store(iPoint,iVar));
    }
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&myNormDelta, &normDelta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&myNormDeltaNew, &normDeltaNew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&myHelper, &helper, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  normDelta    = myNormDelta;
  normDeltaNew = myNormDeltaNew;
  helper       = myHelper;
#endif

  rho   = min(max(sqrt(normDeltaNew)/sqrt(normDelta), 0.9*rho_old), 1.0-1E-6); // Saturate contractivity
  theta = max(sqrt(fabs(helper)/normDelta*theta_old), 0.9*theta_old);
}

void COneShotSolver::CalculateAlphaBetaGamma(CConfig *config, su2double val_bcheck_norm){

  /* --- Estimate alpha, beta, and gamma values --- */
  // su2double alpha = 2./((1.-rho)*(1.-rho));
  // su2double beta  = 2.;
  // su2double gamma = 2./(beta*val_bcheck_norm);
  su2double alpha = 2.*theta/((1.-rho)*(1.-rho));
  su2double beta  = 2./theta;
  su2double gamma = 2./(beta*val_bcheck_norm);

  config->SetOneShotAlpha(alpha);
  config->SetOneShotBeta(beta);
  config->SetOneShotGamma(gamma);

  /* --- Store rho and theta values for this iteration --- */
  rho_old   = rho;
  theta_old = theta;
}

su2double COneShotSolver::CalculateLagrangian(CConfig *config){
  unsigned short iVar;
  unsigned long iPoint;
  su2double Lagrangian=0.0, myLagrangian=0.0;
  su2double helper=0.0;

  /* --- Calculate augmented Lagrangian terms (alpha and beta) --- */
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      helper+=direct_solver->GetNodes()->GetSolution_Delta(iPoint,iVar)*direct_solver->GetNodes()->GetSolution_Delta(iPoint,iVar);
    }
  }
  myLagrangian+=helper*(config->GetOneShotAlpha()/2);
  helper=0.0;
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      helper+=nodes->GetSolution_Delta(iPoint,iVar)*nodes->GetSolution_Delta(iPoint,iVar);
    }
  }
  myLagrangian+=helper*(config->GetOneShotBeta()/2);

  helper=0.0;
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      helper+=direct_solver->GetNodes()->GetSolution_Delta(iPoint,iVar)*nodes->GetSolution_Store(iPoint,iVar);
    }
  }
  myLagrangian+=helper;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&myLagrangian, &Lagrangian, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  Lagrangian = myLagrangian;
#endif

  return Lagrangian;
}

void COneShotSolver::SetAdjoint_OutputUpdate(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->GetNodes()->SetAdjointSolution(iPoint,direct_solver->GetNodes()->GetSolution_Delta(iPoint));
  }
}

void COneShotSolver::SetAdjoint_OutputZero(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;
  unsigned short iVar;
  su2double *ZeroSolution = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++){
      ZeroSolution[iVar] = 0.0;
  }

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->GetNodes()->SetAdjointSolution(iPoint,ZeroSolution);
  }

  delete [] ZeroSolution;
}

void COneShotSolver::ExtractAdjoint_Solution_Clean(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Extract the adjoint solution ---*/

    direct_solver->GetNodes()->GetAdjointSolution(iPoint,Solution);

    /*--- Store the adjoint solution ---*/

    nodes->SetSolution(iPoint,Solution);
  }

}

void COneShotSolver::UpdateStateVariable(CConfig *config){
  unsigned long iPoint;
  unsigned short iVar;
  su2double fd_step=config->GetFDStep();
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      Solution[iVar] = direct_solver->GetNodes()->GetSolution_Store(iPoint,iVar)+fd_step*nodes->GetSolution_Delta(iPoint,iVar);
    }
    direct_solver->GetNodes()->SetSolution(iPoint,Solution);
  }
}

void COneShotSolver::SetFiniteDifferenceSens(CGeometry *geometry, CConfig* config){
  unsigned short iDim;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      nodes->SetSensitivity(iPoint,iDim, (nodes->GetSensitivity(iPoint,iDim)-nodes->GetSensitivity_ShiftedLagrangian(iPoint,iDim))*(1./config->GetFDStep()));
    }
  }
}

void COneShotSolver::SetSolutionDelta(CGeometry *geometry){
  unsigned short iVar;
  unsigned long iPoint;

  // for(iVar = 0; iVar < nVar; iVar++) {
  //   SetRes_RMS(iVar,0.0);
  //   SetRes_Max(iVar,0.0,0);
  //   direct_solver->SetRes_RMS(iVar,0.0);
  //   direct_solver->SetRes_Max(iVar,0.0,0);
  // }

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      const su2double res_y    = direct_solver->GetNodes()->GetSolution(iPoint,iVar)-direct_solver->GetNodes()->GetSolution_Store(iPoint,iVar);
      const su2double res_bary = nodes->GetSolution(iPoint,iVar)-nodes->GetSolution_Store(iPoint,iVar);

      direct_solver->GetNodes()->SetSolution_Delta(iPoint, iVar, res_y);
      nodes->SetSolution_Delta(iPoint, iVar, res_bary);

      // direct_solver->AddRes_RMS(iVar,res_y*res_y);
      // direct_solver->AddRes_Max(iVar,fabs(res_y),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());

      // AddRes_RMS(iVar,res_bary*res_bary);
      // AddRes_Max(iVar,fabs(res_bary),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
    }
  }
}

void COneShotSolver::SetStoreSolutionDelta(){
  unsigned short iVar;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for(iVar = 0; iVar < nVar; iVar++) {
      direct_solver->GetNodes()->SetSolution_DeltaStore(iPoint, iVar, direct_solver->GetNodes()->GetSolution_Delta(iPoint,iVar));
      nodes->SetSolution_DeltaStore(iPoint, iVar, nodes->GetSolution_Delta(iPoint,iVar));
    }
  }
}

void COneShotSolver::SetConstrDerivative(unsigned short iConstr){
  unsigned short iVar;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      DConsVec[iConstr][iPoint][iVar]=nodes->GetSolution(iPoint,iVar);
    }
  }

}

su2double COneShotSolver::MultiplyConstrDerivative(unsigned short iConstr, unsigned short jConstr){
  unsigned short iVar;
  unsigned long iPoint;
  su2double product = 0.0, myProduct=0.0;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      myProduct+= DConsVec[iConstr][iPoint][iVar]*DConsVec[jConstr][iPoint][iVar];
    }
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&myProduct, &product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  product = myProduct;
#endif

  return product;
}
