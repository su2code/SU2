/*!
 * \file COneShotSolver.cpp
 * \brief Main subroutines for solving the one-shot problem.
 * \author T.Dick
 * \version 7.0.6 "Blackbird"
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

 theta = 0.0;
 rho = 0.0;
 nConstr = config->GetnConstr();

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


  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)),
  time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

  unsigned long iPoint;
  unsigned short iVar;

  /*--- For the one-shot solver the solution is not reset in each iteration step to the initial solution ---*/

  if (time_n_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n()[iVar]);
      }
    }
  }
  if (time_n1_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n1()[iVar]);
      }
    }
  }

  /*--- Set the Jacobian to zero since this is not done inside the fluid iteration
   * when running the discrete adjoint solver. ---*/

  direct_solver->Jacobian.SetValZero();

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}

void COneShotSolver::SetGeometrySensitivityLagrangian(CGeometry *geometry){

    unsigned short iDim;
    unsigned long iPoint;

    geometry->InitializeSensitivity();

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        geometry->SetSensitivity(iPoint, iDim, node[iPoint]->GetSensitivity_AugmentedLagrangian(iDim));
      }
    }
}

void COneShotSolver::SetGeometrySensitivityGradient(CGeometry *geometry){

    unsigned short iDim;
    unsigned long iPoint;

    geometry->InitializeSensitivity();

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        geometry->SetSensitivity(iPoint, iDim, node[iPoint]->GetSensitivity_ShiftedLagrangian(iDim));
      }
    }
}

void COneShotSolver::SaveSensitivity(CGeometry *geometry){
    unsigned short iDim;
    unsigned long iPoint;

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        node[iPoint]->SetSensitivity_ShiftedLagrangian(iDim, node[iPoint]->GetSensitivity(iDim));
      }
    }
}

void COneShotSolver::ResetSensitivityLagrangian(CGeometry *geometry){
    unsigned short iDim;
    unsigned long iPoint;

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        node[iPoint]->SetSensitivity_AugmentedLagrangian(iDim, 0.0);
      }
    }
}

void COneShotSolver::UpdateSensitivityLagrangian(CGeometry *geometry, su2double factor){
    unsigned short iDim;
    unsigned long iPoint;

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        node[iPoint]->SetSensitivity_AugmentedLagrangian(iDim, node[iPoint]->GetSensitivity_AugmentedLagrangian(iDim)+factor*node[iPoint]->GetSensitivity(iDim));
      }
    }
}

void COneShotSolver::StoreMeshPoints(CConfig *config, CGeometry *geometry){
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

void COneShotSolver::LoadMeshPoints(CConfig *config, CGeometry *geometry){
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

void COneShotSolver::StoreSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->Set_StoreSolution();
    node[iPoint]->Set_StoreSolution();
  }
}

void COneShotSolver::StoreFormerSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->Set_FormerSolution();
    node[iPoint]->Set_FormerSolution();
  }
}

void COneShotSolver::LoadSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->SetSolution(direct_solver->node[iPoint]->GetSolution_Store());
    node[iPoint]->SetSolution(node[iPoint]->GetSolution_Store());
  }
}

void COneShotSolver::LoadSolutionStep(su2double stepsize){
  unsigned long iPoint, iVar;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      direct_solver->node[iPoint]->SetSolution(iVar, direct_solver->node[iPoint]->GetSolution_Former(iVar)+stepsize*direct_solver->node[iPoint]->GetSolution_Delta_Store(iVar));
      node[iPoint]->SetSolution(iVar, node[iPoint]->GetSolution_Former(iVar)+stepsize*node[iPoint]->GetSolution_Delta_Store(iVar));
    }
  }
}

void COneShotSolver::ShiftFormerSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->SetSolution_Store(direct_solver->node[iPoint]->GetSolution_Former());
    node[iPoint]->SetSolution_Store(node[iPoint]->GetSolution_Former());
  }
}

void COneShotSolver::ShiftStoreSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->SetSolution_Former(direct_solver->node[iPoint]->GetSolution_Store());
    node[iPoint]->SetSolution_Former(node[iPoint]->GetSolution_Store());
  }
}

void COneShotSolver::StoreSaveSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->Set_SaveSolution();
    node[iPoint]->Set_SaveSolution();
  }
}

void COneShotSolver::LoadSaveSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->SetSolution(direct_solver->node[iPoint]->GetSolution_Save());
    node[iPoint]->SetSolution(node[iPoint]->GetSolution_Save());
  }
}

void COneShotSolver::CalculateAlphaBetaGamma(CConfig *config){
  unsigned short iVar;
  unsigned long iPoint;
  su2double normDelta=0.0,    myNormDelta=0.0;
  su2double normDeltaNew=0.0, myNormDeltaNew=0.0;

  /* --- Estimate rho and theta values --- */
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      myNormDelta += direct_solver->node[iPoint]->GetSolution_Delta(iVar)*direct_solver->node[iPoint]->GetSolution_Delta(iVar);
      myNormDeltaNew += (direct_solver->node[iPoint]->GetSolution(iVar)-direct_solver->node[iPoint]->GetSolution_Store(iVar))*(direct_solver->node[iPoint]->GetSolution(iVar)-direct_solver->node[iPoint]->GetSolution_Store(iVar));
    }
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&myNormDelta, &normDelta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&myNormDeltaNew, &normDeltaNew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  normDelta    = myNormDelta;
  normDeltaNew = myNormDeltaNew;
#endif

  rho = min(max(sqrt(normDeltaNew)/sqrt(normDelta), 0.9*rho), 0.9999); // Saturate contractivity
}

void COneShotSolver::SetAlphaBetaGamma(CConfig *config, su2double val_bcheck_norm){

  su2double alpha = 2./((1.-rho)*(1.-rho));
  su2double beta  = 2.;
  su2double gamma = 2./val_bcheck_norm;

  config->SetOneShotAlpha(alpha);
  config->SetOneShotBeta(beta);
  config->SetOneShotGamma(gamma);
}

su2double COneShotSolver::CalculateLagrangianPart(CConfig *config, bool augmented){
  unsigned short iVar;
  unsigned long iPoint;
  su2double Lagrangian=0.0, myLagrangian=0.0;
  su2double helper=0.0;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      direct_solver->node[iPoint]->SetSolution_Delta(iVar, direct_solver->node[iPoint]->GetSolution(iVar)-direct_solver->node[iPoint]->GetSolution_Store(iVar));
    }
    for (iVar = 0; iVar < nVar; iVar++){
      node[iPoint]->SetSolution_Delta(iVar,node[iPoint]->GetSolution(iVar)-node[iPoint]->GetSolution_Store(iVar));
    }
  }

  /* --- Calculate augmented Lagrangian terms (alpha and beta) --- */
  if(augmented){
    for (iPoint = 0; iPoint < nPointDomain; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        helper+=direct_solver->node[iPoint]->GetSolution_Delta(iVar)*direct_solver->node[iPoint]->GetSolution_Delta(iVar);
      }
    }
    myLagrangian+=helper*(config->GetOneShotAlpha()/2);
    helper=0.0;
    for (iPoint = 0; iPoint < nPointDomain; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        helper+=node[iPoint]->GetSolution_Delta(iVar)*node[iPoint]->GetSolution_Delta(iVar);
      }
    }
    myLagrangian+=helper*(config->GetOneShotBeta()/2);
  }

  helper=0.0;
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      helper+=direct_solver->node[iPoint]->GetSolution_Delta(iVar)*node[iPoint]->GetSolution_Store(iVar);
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
    direct_solver->node[iPoint]->SetAdjointSolution(direct_solver->node[iPoint]->GetSolution_Delta());
  }
}

void COneShotSolver::SetAdjoint_OutputZero(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;
  unsigned short iVar;
  su2double * ZeroSolution = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++){
      ZeroSolution[iVar] = 0.0;
  }

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->node[iPoint]->SetAdjointSolution(ZeroSolution);
  }

  delete [] ZeroSolution;
}

void COneShotSolver::ExtractAdjoint_Solution_Clean(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Extract the adjoint solution ---*/

    direct_solver->node[iPoint]->GetAdjointSolution(Solution);

    /*--- Store the adjoint solution ---*/

    node[iPoint]->SetSolution(Solution);
  }

}

void COneShotSolver::UpdateStateVariable(CConfig *config){
    unsigned long iPoint;
    unsigned short iVar;
    su2double fd_step=config->GetFDStep();
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        Solution[iVar] = direct_solver->node[iPoint]->GetSolution_Store(iVar)+fd_step*node[iPoint]->GetSolution_Delta(iVar);
      }
      direct_solver->node[iPoint]->SetSolution(Solution);
    }
}

void COneShotSolver::SetFiniteDifferenceSens(CGeometry *geometry, CConfig* config){
    unsigned short iDim;
    unsigned long iPoint;

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        node[iPoint]->SetSensitivity(iDim, (node[iPoint]->GetSensitivity(iDim)-node[iPoint]->GetSensitivity_ShiftedLagrangian(iDim))*(1./config->GetFDStep()));
      }
    }
}

void COneShotSolver::StoreSolutionDelta(){
  unsigned short iVar;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      direct_solver->node[iPoint]->SetSolution_Delta_Store(iVar, direct_solver->node[iPoint]->GetSolution_Delta(iVar));
    }
    for (iVar = 0; iVar < nVar; iVar++){
      node[iPoint]->SetSolution_Delta_Store(iVar,node[iPoint]->GetSolution_Delta(iVar));
    }
  }
}

void COneShotSolver::SetConstrDerivative(unsigned short iConstr){
  unsigned short iVar;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      DConsVec[iConstr][iPoint][iVar]=node[iPoint]->GetSolution(iVar);
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
