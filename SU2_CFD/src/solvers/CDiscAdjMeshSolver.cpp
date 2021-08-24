/*!
 * \file CDiscAdjMeshSolver.cpp
 * \brief Main subroutines for solving the discrete adjoint mesh problem.
 * \author Ruben Sanchez
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

#include "../../include/solvers/CDiscAdjMeshSolver.hpp"
#include "../../include/variables/CDiscAdjMeshBoundVariable.hpp"


CDiscAdjMeshSolver::CDiscAdjMeshSolver() : CSolver () {}

CDiscAdjMeshSolver::CDiscAdjMeshSolver(CGeometry *geometry, CConfig *config) : CSolver() {}

CDiscAdjMeshSolver::CDiscAdjMeshSolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver) : CSolver() {

  unsigned short iVar;

  nVar = geometry->GetnDim();
  nDim = geometry->GetnDim();

  /*-- Store some information about direct solver ---*/
  this->direct_solver = direct_solver;

  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual_RMS.resize(nVar,1.0);
  Residual_Max.resize(nVar,1.0);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

  /*--- Define some auxiliary vectors related to the residual for problems with a BGS strategy---*/

  if (config->GetMultizone_Residual()){

    Residual_BGS.resize(nVar,1.0);
    Residual_Max_BGS.resize(nVar,1.0);
    Point_Max_BGS.resize(nVar,0);
    Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
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
      long iVertex = geometry->nodes->GetVertex(iPoint, iMarker);
      if (iVertex >= 0) {
        nodes->Set_isVertex(iPoint,true);
        break;
      }
    }
  }
  static_cast<CDiscAdjMeshBoundVariable*>(nodes)->AllocateBoundaryVariables(config);

}

CDiscAdjMeshSolver::~CDiscAdjMeshSolver(void){
  delete nodes;
}


void CDiscAdjMeshSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config_container,
                                       unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output){
}

void CDiscAdjMeshSolver::SetRecording(CGeometry* geometry, CConfig *config){

  /*--- Reset the solution to the initial (converged) solution ---*/

  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
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

  /*--- Extract the sensitivities of the mesh coordinates ---*/

  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++){

    /*--- Extract the adjoint solution from the original mesh coordinates ---*/

    direct_solver->GetNodes()->GetAdjoint_MeshCoord(iPoint,Solution);

    /*--- Store the adjoint solution (the container is reused) ---*/

    nodes->SetSolution(iPoint,Solution);

  }

}

void CDiscAdjMeshSolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config){

  /*--- Extract the sensitivities of the boundary displacements ---*/

  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++){

    /*--- Extract the adjoint solution of the boundary displacements ---*/

    direct_solver->GetNodes()->GetAdjoint_BoundDisp(iPoint,Solution);

    /*--- Store the sensitivities of the boundary displacements ---*/

    nodes->SetBoundDisp_Sens(iPoint,Solution);

  }

}

void CDiscAdjMeshSolver::SetSensitivity(CGeometry *geometry, CConfig *config, CSolver *solver) {

  const bool time_stepping = (config->GetTime_Marching() != STEADY);
  const auto eps = config->GetAdjSharp_LimiterCoeff()*config->GetRefElemLength();

  /*--- Extract the sensitivities ---*/
  ExtractAdjoint_Solution(geometry, config);

  /*--- Extract the adjoint variables: sensitivities of the boundary displacements ---*/
  ExtractAdjoint_Variables(geometry, config);

  /*--- Store the sensitivities in the flow adjoint container ---*/
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {

    /*--- If sharp edge, set the sensitivity to 0 on that region ---*/
    su2double limiter = 1.0;
    if (config->GetSens_Remove_Sharp() && (geometry->nodes->GetSharpEdge_Distance(iPoint) < eps)) {
      limiter = 0.0;
    }

    for (auto iDim = 0u; iDim < nDim; iDim++) {
      /*--- The sensitivity was extracted using ExtractAdjoint_Solution ---*/
      const su2double Sensitivity = limiter * nodes->GetSolution(iPoint,iDim);

      /*--- Store the sensitivities ---*/
      if (!time_stepping) {
        solver->GetNodes()->SetSensitivity(iPoint, iDim, Sensitivity);
      } else {
        solver->GetNodes()->SetSensitivity(iPoint, iDim,
          solver->GetNodes()->GetSensitivity(iPoint, iDim) + Sensitivity);
      }
    }
  }
  solver->SetSurface_Sensitivity(geometry, config);

}

void CDiscAdjMeshSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

}
