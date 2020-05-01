/*!
 * \file CDiscAdjMeshSolver.cpp
 * \brief Main subroutines for solving the discrete adjoint mesh problem.
 * \author Ruben Sanchez
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



#include "../../include/solvers/CDiscAdjMeshSolver.hpp"
#include "../../include/variables/CDiscAdjMeshBoundVariable.hpp"

CDiscAdjMeshSolver::CDiscAdjMeshSolver(void) : CSolver (){

  KindDirect_Solver = 0;

  direct_solver = NULL;

}

CDiscAdjMeshSolver::CDiscAdjMeshSolver(CGeometry *geometry, CConfig *config)  : CSolver(){

  KindDirect_Solver = 0;

}

CDiscAdjMeshSolver::CDiscAdjMeshSolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver)  : CSolver(){

  unsigned short iVar, iDim;

  nVar = geometry->GetnDim();
  nDim = geometry->GetnDim();

  /*-- Store some information about direct solver ---*/
  this->direct_solver = direct_solver;

  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 1.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 1.0;

  /*--- Define some structures for locating max residuals ---*/

  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Define some auxiliary vectors related to the residual for problems with a BGS strategy---*/

  if (config->GetMultizone_Residual()){

    Residual_BGS      = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]      = 1.0;
    Residual_Max_BGS  = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 1.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS       = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar]  = 0;
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
    }

  }

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 1e-16;

  /*--- Initialize the node structure ---*/
  nodes = new CDiscAdjMeshBoundVariable(nPoint,nDim,config);
  SetBaseClassPointerToNodes();

  /*--- Set which points are vertices and allocate boundary data. ---*/

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    nodes->SetSolution(iPoint,Solution);

    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      long iVertex = geometry->node[iPoint]->GetVertex(iMarker);
      if (iVertex >= 0) {
        nodes->Set_isVertex(iPoint,true);
        break;
      }
    }
  }
  static_cast<CDiscAdjMeshBoundVariable*>(nodes)->AllocateBoundaryVariables(config);

}

CDiscAdjMeshSolver::~CDiscAdjMeshSolver(void){
  if (nodes != nullptr) delete nodes;
}


void CDiscAdjMeshSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config_container, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output){


}

void CDiscAdjMeshSolver::SetRecording(CGeometry* geometry, CConfig *config){


  unsigned long iPoint;
  /*--- Reset the solution to the initial (converged) solution ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->GetNodes()->SetBound_Disp(iPoint,nodes->GetBoundDisp_Direct(iPoint));
  }

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}

void CDiscAdjMeshSolver::RegisterSolution(CGeometry *geometry, CConfig *config){

  /*--- Register reference mesh coordinates ---*/
  bool input = true;
  direct_solver->GetNodes()->Register_MeshCoord(input);

}

void CDiscAdjMeshSolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset){

  /*--- Register boundary displacements as input ---*/
  bool input = true;
  direct_solver->GetNodes()->Register_BoundDisp(input);

}

void CDiscAdjMeshSolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;

  /*--- Extract the sensitivities of the mesh coordinates ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /*--- Extract the adjoint solution from the original mesh coordinates ---*/

    direct_solver->GetNodes()->GetAdjoint_MeshCoord(iPoint,Solution);

    /*--- Store the adjoint solution (the container is reused) ---*/

    nodes->SetSolution(iPoint,Solution);

  }

}

void CDiscAdjMeshSolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;

  /*--- Extract the sensitivities of the boundary displacements ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /*--- Extract the adjoint solution of the boundary displacements ---*/

    direct_solver->GetNodes()->GetAdjoint_BoundDisp(iPoint,Solution);

    /*--- Store the sensitivities of the boundary displacements ---*/

    nodes->SetBoundDisp_Sens(iPoint,Solution);

  }

}

void CDiscAdjMeshSolver::SetSensitivity(CGeometry *geometry, CSolver **solver, CConfig *config) {

  unsigned long iPoint;
  unsigned short iDim;
  su2double Sensitivity, eps;
  bool time_stepping = (config->GetTime_Marching() != STEADY);

  /*--- Extract the sensitivities ---*/
  ExtractAdjoint_Solution(geometry, config);

  /*--- Extract the adjoint variables: sensitivities of the boundary displacements ---*/
  ExtractAdjoint_Variables(geometry, config);

  /*--- Store the sensitivities in the flow adjoint container ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    for (iDim = 0; iDim < nDim; iDim++) {

      /*--- The sensitivity was extracted using ExtractAdjoint_Solution ---*/
      Sensitivity = nodes->GetSolution(iPoint,iDim);

      /*--- If sharp edge, set the sensitivity to 0 on that region ---*/
      if (config->GetSens_Remove_Sharp()) {
        eps = config->GetVenkat_LimiterCoeff()*config->GetRefElemLength();
        if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetAdjSharp_LimiterCoeff()*eps )
          Sensitivity = 0.0;
      }

      /*--- Store the sensitivities ---*/
      if (!time_stepping) {
        solver[ADJFLOW_SOL]->GetNodes()->SetSensitivity(iPoint, iDim, Sensitivity);
      } else {
        solver[ADJFLOW_SOL]->GetNodes()->SetSensitivity(iPoint, iDim,
          solver[ADJFLOW_SOL]->GetNodes()->GetSensitivity(iPoint, iDim) + Sensitivity);
      }
    }
  }
  solver[ADJFLOW_SOL]->SetSurface_Sensitivity(geometry, config);

}

void CDiscAdjMeshSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

}
