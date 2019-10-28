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

COneShotSolver::COneShotSolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver, unsigned short Kind_Solver, unsigned short iMesh)  : CDiscAdjSolver() {

  unsigned short iVar, iMarker, iDim;
  unsigned long iVertex;
  string text_line, mesh_filename;
  ifstream restart_file;
  string filename, AdjExt;

  adjoint = true;

  nVar = direct_solver->GetnVar();
  nDim = geometry->GetnDim();

  /*--- Initialize arrays to NULL ---*/

  CSensitivity = NULL;

  /*-- Store some information about direct solver ---*/
  this->KindDirect_Solver = Kind_Solver;
  this->direct_solver = direct_solver;


  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 1.0;
  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 1.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 1.0;

  /*--- Define some auxiliary vectors related to the geometry adjoint (nDim) ---*/
  Solution_Geometry = new su2double[nDim];     for (iDim = 0; iDim < nDim; iDim++) Solution_Geometry[iDim] = 1.0;

  /*--- Define some auxiliary vectors related to the residual for problems with a BGS strategy---*/

  if (config->GetMultizone_Residual()){

    Residual_BGS      = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]      = 1.0;
    Residual_Max_BGS  = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 1.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS       = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar] = 0;
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
    }

  }


  /*--- Define some structures for locating max residuals ---*/

  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution   = new su2double[nVar];

  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 1e-16;

  /*--- Sensitivity definition and coefficient in all the markers ---*/

  CSensitivity = new su2double* [nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      CSensitivity[iMarker]        = new su2double [geometry->nVertex[iMarker]];
  }

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          CSensitivity[iMarker][iVertex] = 0.0;
      }
  }

  /*--- Initialize the discrete adjoint solution to zero everywhere. ---*/

  nodes = new CDiscAdjVariable(Solution, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  switch(KindDirect_Solver){
  case RUNTIME_FLOW_SYS:
    SolverName = "ADJ.FLOW";
    break;
  case RUNTIME_HEAT_SYS:
    SolverName = "ADJ.HEAT";
    break;
  case RUNTIME_TURB_SYS:
    SolverName = "ADJ.TURB";
    break;
  default:
    SolverName = "ADJ.SOL";
    break;
  }

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

  unsigned short iMarker;

  if (CSensitivity != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] CSensitivity[iMarker];
    }
    delete [] CSensitivity;
  }

  if (nodes != nullptr) delete nodes;
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

void COneShotSolver::SaveSensitivity(CGeometry *geometry){
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
    direct_solver->GetNodes()->Set_StoreSolution(iPoint);
    nodes->Set_StoreSolution(iPoint);
  }
}

void COneShotSolver::StoreFormerSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->Set_FormerSolution(iPoint);
    nodes->Set_FormerSolution(iPoint);
  }
}

void COneShotSolver::LoadSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution(iPoint,direct_solver->GetNodes()->GetSolution_Store(iPoint));
    nodes->SetSolution(iPoint,nodes->GetSolution_Store(iPoint));
  }
}

void COneShotSolver::LoadSolutionStep(su2double stepsize){
  unsigned long iPoint, iVar;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      direct_solver->GetNodes()->SetSolution(iPoint, iVar, direct_solver->GetNodes()->GetSolution_Former(iPoint,iVar)+stepsize*direct_solver->GetNodes()->GetSolution_Delta_Store(iPoint,iVar));
      nodes->SetSolution(iPoint, iVar, nodes->GetSolution_Former(iPoint,iVar)+stepsize*nodes->GetSolution_Delta_Store(iPoint,iVar));
    }
  }
}

void COneShotSolver::ShiftFormerSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution_Store(iPoint,direct_solver->GetNodes()->GetSolution_Former(iPoint));
    nodes->SetSolution_Store(iPoint,nodes->GetSolution_Former(iPoint));
  }
}

void COneShotSolver::ShiftStoreSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution_Former(iPoint,direct_solver->GetNodes()->GetSolution_Store(iPoint));
    nodes->SetSolution_Former(iPoint,nodes->GetSolution_Store(iPoint));
  }
}

void COneShotSolver::StoreSaveSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->Set_SaveSolution(iPoint);
    nodes->Set_SaveSolution(iPoint);
  }
}

void COneShotSolver::LoadSaveSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution(iPoint,direct_solver->GetNodes()->GetSolution_Save(iPoint));
    nodes->SetSolution(iPoint,nodes->GetSolution_Save(iPoint));
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
      myNormDelta += direct_solver->GetNodes()->GetSolution_Delta(iPoint,iVar)*direct_solver->GetNodes()->GetSolution_Delta(iPoint,iVar);
      myNormDeltaNew += (direct_solver->GetNodes()->GetSolution(iPoint,iVar)-direct_solver->GetNodes()->GetSolution_Store(iPoint,iVar))*(direct_solver->GetNodes()->GetSolution(iPoint,iVar)-direct_solver->GetNodes()->GetSolution_Store(iPoint,iVar));
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
      direct_solver->GetNodes()->SetSolution_Delta(iPoint, iVar, direct_solver->GetNodes()->GetSolution(iPoint,iVar)-direct_solver->GetNodes()->GetSolution_Store(iPoint,iVar));
    }
    for (iVar = 0; iVar < nVar; iVar++){
      nodes->SetSolution_Delta(iPoint,iVar,nodes->GetSolution(iPoint,iVar)-nodes->GetSolution_Store(iPoint,iVar));
    }
  }

  /* --- Calculate augmented Lagrangian terms (alpha and beta) --- */
  if(augmented){
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
  }

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
  su2double * ZeroSolution = new su2double[nVar];
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

void COneShotSolver::StoreSolutionDelta(){
  unsigned short iVar;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      direct_solver->GetNodes()->SetSolution_Delta_Store(iPoint, iVar, direct_solver->GetNodes()->GetSolution_Delta(iPoint,iVar));
    }
    for (iVar = 0; iVar < nVar; iVar++){
      nodes->SetSolution_Delta_Store(iPoint,iVar,nodes->GetSolution_Delta(iPoint,iVar));
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
