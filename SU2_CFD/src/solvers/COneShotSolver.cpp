/*!
 * \file COneShotSolver.cpp
 * \brief Main subroutines for solving the one-shot problem.
 * \author T.Dick
 * \version 7.1.1 "Blackbird"
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

#include "../../include/solvers/COneShotSolver.hpp"
#include "../../include/variables/CDiscAdjVariable.hpp"

COneShotSolver::COneShotSolver(void) : CDiscAdjSolver () {

}

COneShotSolver::COneShotSolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver, unsigned short Kind_Solver, unsigned short iMesh)  : CDiscAdjSolver(geometry, config, direct_solver, Kind_Solver, iMesh) {

}

COneShotSolver::~COneShotSolver(void) {

}

void COneShotSolver::SetRecording(CGeometry* geometry, CConfig *config){


  bool time_n1_needed = config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND;
  bool time_n_needed = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) || time_n1_needed;

  unsigned long iPoint;
  unsigned short iVar;

  /*--- For the one-shot solver the solution is not reset in each iteration step to the initial solution ---*/

  if (time_n_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->GetNodes()->GetSolution_time_n(iPoint)[iVar]);
      }
    }
  }
  if (time_n1_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->GetNodes()->GetSolution_time_n1(iPoint)[iVar]);
      }
    }
  }

  /*--- Set the Jacobian to zero since this is not done inside the fluid iteration
   * when running the discrete adjoint solver. ---*/

  direct_solver->Jacobian.SetValZero();

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}


void COneShotSolver::StoreMeshPoints(CGeometry *geometry, CConfig *config){
    unsigned long iVertex;
    unsigned short iMarker;

    geometry->nodes->SetCoord_Old();

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          geometry->vertex[iMarker][iVertex]->SetNormal_Old(geometry->vertex[iMarker][iVertex]->GetNormal());
        }
    }
}


void COneShotSolver::LoadMeshPoints(CGeometry *geometry, CConfig *config){
    unsigned long iVertex, iPoint;
    unsigned short iMarker;

    geometry->nodes->GetCoord_Old();

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          geometry->vertex[iMarker][iVertex]->SetNormal(geometry->vertex[iMarker][iVertex]->GetNormal_Old());
        }
    }
}


void  COneShotSolver::UpdateAuxiliaryGeometryVariables(CGeometry **geometry_container, CVolumetricMovement *grid_movement, CConfig *config) {

  su2double MinVolume, MaxVolume;

  /*--- Communicate the updated mesh coordinates. ---*/

  geometry_container[MESH_0]->InitiateComms(geometry_container[MESH_0], config, COORDINATES);
  geometry_container[MESH_0]->CompleteComms(geometry_container[MESH_0], config, COORDINATES);

  /*--- After moving all nodes, update the dual mesh. Recompute the edges and
   dual mesh control volumes in the domain and on the boundaries. ---*/

  grid_movement->UpdateDualGrid(geometry_container[MESH_0], config);

  /*--- Update the multigrid structure after moving the finest grid,
   including computing the grid velocities on the coarser levels
   when the problem is solved in unsteady conditions. ---*/

  grid_movement->UpdateMultiGrid(geometry_container, config);

  /*--- do we need to update the elemnt volumes, etc.??? ---*/

  grid_movement->ComputeDeforming_Element_Volume(geometry_container[MESH_0], MinVolume, MaxVolume, true);
  grid_movement->ComputenNonconvexElements(geometry_container[MESH_0], true);

  if (rank == MASTER_NODE) {
    cout << "Resetting mesh coordinates for linesearch: " << endl;
    cout << "Min. volume: " << MinVolume << ", Max. volume: " << MaxVolume << "." << endl;
  }

}

su2double COneShotSolver::EvaluateGeometryFunction(CGeometry* geometry, CConfig* config, unsigned int iPlane) {

  /*--- Evaluation of the objective function ---*/

  unsigned short nPlane;
  su2double ThicknessValue;
  su2double *Plane_P0, *Plane_Normal;
  su2double Wing_Volume = 0.0, Wing_MinThickness = 0.0, Wing_MaxThickness = 0.0, Wing_MinChord = 0.0, Wing_MaxChord = 0.0, Wing_MinLERadius = 0.0, Wing_MaxLERadius = 0.0, Wing_MinToC = 0.0, Wing_MaxToC = 0.0, Wing_ObjFun_MinToC = 0.0, Wing_MaxTwist = 0.0, Wing_MaxCurvature = 0.0, Wing_MaxDihedral = 0.0;

  std::vector<su2double> Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil, Variable_Airfoil;

  /*--- Set the number of sections, and allocate the memory ---*/

  Plane_P0 = new su2double[3];
  Plane_Normal = new su2double[3];

  /// possible advanced preprocessing???

  geometry->ComputeSurf_Curvature(config);
  geometry->SetPositive_ZArea(config);

  /*--- Create plane structure ---*/

  if (geometry->GetnDim() == 2) {
    Plane_Normal[0] = 0.0;   Plane_P0[0] = 0.0;
    Plane_Normal[1] = 1.0;   Plane_P0[1] = 0.0;
    Plane_Normal[2] = 0.0;   Plane_P0[2] = 0.0;
  }
  else if (geometry->GetnDim() == 3) {
    Plane_Normal[0] = 0.0;    Plane_P0[0] = 0.0;
    Plane_Normal[1] = 1.0;    Plane_P0[1] = config->GetLocationStations(iPlane);
    Plane_Normal[2] = 0.0;    Plane_P0[2] = 0.0;
  }

  /*--- Compute the wing and airfoil description ---*/

  if (geometry->GetnDim() == 3) {
    geometry->Compute_Wing(config, true,
                           Wing_Volume, Wing_MinThickness, Wing_MaxThickness, Wing_MinChord, Wing_MaxChord,
                           Wing_MinLERadius, Wing_MaxLERadius, Wing_MinToC, Wing_MaxToC, Wing_ObjFun_MinToC,
                           Wing_MaxTwist, Wing_MaxCurvature, Wing_MaxDihedral);
  }
  geometry->ComputeAirfoil_Section(Plane_P0, Plane_Normal, -1E6, 1E6, -1E6, 1E6, -1E6, 1E6, nullptr,
                                     Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil,
                                     Variable_Airfoil, true, config);

  if (Xcoord_Airfoil.size() > 1) {
    ThicknessValue  = geometry->Compute_MaxThickness(Plane_P0, Plane_Normal, config, Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil);
  }

  /// MPI broadcast the solution
  SU2_MPI::Bcast(&ThicknessValue, 1, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());

  return ThicknessValue;

}


vector<su2double> COneShotSolver::EvaluateGeometryGradient(CGeometry* geometry, CSurfaceMovement* surface_movement, CConfig *config) {


  ///* TODO: check use of nDV and nDV_Total here? *///

  unsigned short  iDV, iFFDBox, iPlane, nPlane;
  su2double *ThicknessValue, *ThicknessValue_New, *Gradient, delta_eps,
            **Plane_P0, **Plane_Normal;
  su2double Wing_Volume = 0.0, Wing_MinThickness = 0.0, Wing_MaxThickness = 0.0, Wing_MinChord = 0.0, Wing_MaxChord = 0.0, Wing_MinLERadius = 0.0, Wing_MaxLERadius = 0.0, Wing_MinToC = 0.0, Wing_MaxToC = 0.0, Wing_ObjFun_MinToC = 0.0, Wing_MaxTwist = 0.0, Wing_MaxCurvature = 0.0, Wing_MaxDihedral = 0.0,
            Wing_Volume_New = 0.0, Wing_MinThickness_New = 0.0, Wing_MaxThickness_New = 0.0, Wing_MinChord_New = 0.0, Wing_MaxChord_New = 0.0, Wing_MinLERadius_New = 0.0, Wing_MaxLERadius_New = 0.0, Wing_MinToC_New = 0.0, Wing_MaxToC_New = 0.0, Wing_ObjFun_MinToC_New = 0.0, Wing_MaxTwist_New = 0.0, Wing_MaxCurvature_New = 0.0, Wing_MaxDihedral_New = 0.0;
  std::vector<su2double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, *Variable_Airfoil;
  std::vector<su2double> Sum_Gradient(config->GetnDV(), 0.0);
  su2double* SumGradientBuffer;
  bool Local_MoveSurface, MoveSurface = false;

  SumGradientBuffer = new su2double[config->GetnDV()];
  for (auto iDV=0; iDV<config->GetnDV(); iDV++) {
    SumGradientBuffer[iDV]=0.0;
  }

  /// base evaluation

  /*--- Set the number of sections, and allocate the memory ---*/

  if (geometry->GetnDim() == 2) nPlane = 1;
  else nPlane = config->GetnLocationStations();

  Xcoord_Airfoil = new vector<su2double>[nPlane];
  Ycoord_Airfoil = new vector<su2double>[nPlane];
  Zcoord_Airfoil = new vector<su2double>[nPlane];
  Variable_Airfoil = new vector<su2double>[nPlane];

  Plane_P0 = new su2double*[nPlane];
  Plane_Normal = new su2double*[nPlane];
  for(iPlane = 0; iPlane < nPlane; iPlane++ ) {
    Plane_P0[iPlane] = new su2double[3];
    Plane_Normal[iPlane] = new su2double[3];
  }

  ThicknessValue = new su2double[nPlane];
  ThicknessValue_New = new su2double[nPlane];
  Gradient = new su2double[nPlane];

  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    ThicknessValue[iPlane] = 0.0;
    ThicknessValue_New[iPlane] = 0.0;
    Gradient[iPlane] = 0.0;
  }

  geometry->ComputeSurf_Curvature(config);
  geometry->SetPositive_ZArea(config);

  /*--- Create plane structure ---*/

  if (geometry->GetnDim() == 2) {
    Plane_Normal[0][0] = 0.0;   Plane_P0[0][0] = 0.0;
    Plane_Normal[0][1] = 1.0;   Plane_P0[0][1] = 0.0;
    Plane_Normal[0][2] = 0.0;   Plane_P0[0][2] = 0.0;
  }
  else if (geometry->GetnDim() == 3) {
    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      Plane_Normal[iPlane][0] = 0.0;    Plane_P0[iPlane][0] = 0.0;
      Plane_Normal[iPlane][1] = 1.0;    Plane_P0[iPlane][1] = config->GetLocationStations(iPlane);
      Plane_Normal[iPlane][2] = 0.0;    Plane_P0[iPlane][2] = 0.0;
    }
  }

  /*--- Compute the wing and fan description (only 3D). ---*/

  if (geometry->GetnDim() == 3) {
      geometry->Compute_Wing(config, true,
                             Wing_Volume, Wing_MinThickness, Wing_MaxThickness, Wing_MinChord, Wing_MaxChord,
                             Wing_MinLERadius, Wing_MaxLERadius, Wing_MinToC, Wing_MaxToC, Wing_ObjFun_MinToC,
                             Wing_MaxTwist, Wing_MaxCurvature, Wing_MaxDihedral);
  }
  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    geometry->ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane], -1E6, 1E6, -1E6, 1E6, -1E6, 1E6, nullptr,
                                     Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane],
                                     Variable_Airfoil[iPlane], true, config);
  }

  if (rank == MASTER_NODE) {
    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      if (Xcoord_Airfoil[iPlane].size() > 1) {
        ThicknessValue[iPlane]  = geometry->Compute_MaxThickness(Plane_P0[iPlane], Plane_Normal[iPlane], config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
      }
    }
  }

  /// gradient part

  /*--- Copy coordinates to the surface structure ---*/
  surface_movement->CopyBoundary(geometry, config);

  /*--- Definition of the FFD deformation class ---*/
  CFreeFormDefBox** FFDBox = nullptr;
  FFDBox = new CFreeFormDefBox*[MAX_NUMBER_FFD];
  for (iFFDBox = 0; iFFDBox < MAX_NUMBER_FFD; iFFDBox++) FFDBox[iFFDBox] = nullptr;

  for (iDV = 0; iDV < config->GetnDV(); iDV++) {

    /*--- Free Form deformation based ---*/

    if ((config->GetDesign_Variable(iDV) == FFD_CONTROL_POINT_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_CONTROL_POINT) ) {

      /*--- Read the FFD information in the first iteration ---*/

      if (iDV == 0) {
        /*--- Read the FFD information from the grid file ---*/
        surface_movement->ReadFFDInfo(geometry, config, FFDBox, config->GetMesh_FileName());
        /*--- If the FFDBox was not defined in the input file ---*/
        if (!surface_movement->GetFFDBoxDefinition()) {
          SU2_MPI::Error("The input grid doesn't have the entire FFD information!", CURRENT_FUNCTION);
        }
        for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
          surface_movement->CheckFFDDimension(geometry, config, FFDBox[iFFDBox], iFFDBox);
          surface_movement->CheckFFDIntersections(geometry, config, FFDBox[iFFDBox], iFFDBox);
        }
      }

      /*--- Apply the control point change ---*/

      MoveSurface = false;

      for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {

        switch ( config->GetDesign_Variable(iDV) ) {
        case FFD_CONTROL_POINT_2D : Local_MoveSurface = surface_movement->SetFFDCPChange_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
        case FFD_CONTROL_POINT :    Local_MoveSurface = surface_movement->SetFFDCPChange(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
        }

        /*--- Recompute cartesian coordinates using the new control points position ---*/

        if (Local_MoveSurface) {
          MoveSurface = true;
          surface_movement->SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, true);
        }

      }
    }

    /*--- Hicks Henne design variable ---*/

    else if (config->GetDesign_Variable(iDV) == HICKS_HENNE) {
      MoveSurface = true;
      surface_movement->SetHicksHenne(geometry, config, iDV, true);
    }

    if (MoveSurface) {

      /*--- Compute the gradient for the volume. In 2D this is just
         the gradient of the area. ---*/

      if (geometry->GetnDim() == 3) {
        geometry->Compute_Wing(config, false,
                               Wing_Volume_New, Wing_MinThickness_New, Wing_MaxThickness_New, Wing_MinChord_New,
                               Wing_MaxChord_New, Wing_MinLERadius_New, Wing_MaxLERadius_New, Wing_MinToC_New, Wing_MaxToC_New,
                               Wing_ObjFun_MinToC_New, Wing_MaxTwist_New, Wing_MaxCurvature_New, Wing_MaxDihedral_New);

      }
      for (iPlane = 0; iPlane < nPlane; iPlane++) {
        geometry->ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane], -1E6, 1E6, -1E6, 1E6, -1E6, 1E6, nullptr,
                                         Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane],
                                         Variable_Airfoil[iPlane], false, config);
      }
    }

    /*--- Compute gradient ---*/

    if (rank == MASTER_NODE) {

      delta_eps = config->GetDV_Value(iDV);

      if (delta_eps == 0) {
        SU2_MPI::Error("The finite difference steps is zero!!", CURRENT_FUNCTION);
      }

      if (MoveSurface) {
        for (iPlane = 0; iPlane < nPlane; iPlane++) {
          if (Xcoord_Airfoil[iPlane].size() > 1) {
            ThicknessValue_New[iPlane] = geometry->Compute_MaxThickness(Plane_P0[iPlane], Plane_Normal[iPlane], config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
            Gradient[iPlane] = (ThicknessValue_New[iPlane] - ThicknessValue[iPlane]) / delta_eps;
          }
        }
      } else {
        for (iPlane = 0; iPlane < nPlane; iPlane++) {
          Gradient[iPlane] = 0.0;
        }
      }
    }

    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      SumGradientBuffer[iDV] += Gradient[iPlane];
    }

  }

  /// MPI broadcast the solution
  SU2_MPI::Bcast(SumGradientBuffer, config->GetnDV(), MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());

  for (auto iDV=0; iDV<config->GetnDV(); iDV++) {
    Sum_Gradient[iDV] = SumGradientBuffer[iDV];
  }

  /// reset boundaries for safety
  surface_movement->CopyBoundary(geometry, config);

  /// delete variables
  delete [] SumGradientBuffer;
  delete [] Xcoord_Airfoil; delete [] Ycoord_Airfoil; delete [] Zcoord_Airfoil;

  for(iPlane = 0; iPlane < nPlane; iPlane++ ) {
    delete [] Plane_P0[iPlane];
    delete [] Plane_Normal[iPlane];
  }
  delete [] Plane_P0;
  delete [] Plane_Normal;

  return Sum_Gradient;

}
