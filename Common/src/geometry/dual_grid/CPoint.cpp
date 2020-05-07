/*!
 * \file CPoint.cpp
 * \brief Main classes for defining the points of the dual grid
 * \author F. Palacios, T. Economon
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

#include "../../../include/geometry/dual_grid/CPoint.hpp"

CPoint::CPoint(unsigned short val_nDim, unsigned long val_globalindex, CConfig *config) : CDualGrid(val_nDim) {

  unsigned short iDim, jDim;

  /*--- Element, point and edge structures initialization ---*/
  Elem.clear();  nElem  = 0;
  Point.clear(); nPoint = 0;
  Edge.clear();

  Volume            = NULL;           Vertex              = NULL;
  Coord             = NULL;           Coord_Old           = NULL;            Coord_Sum  = NULL;
  Coord_n           = NULL;           Coord_n1            = NULL;            Coord_p1   = NULL;
  GridVel           = NULL;           GridVel_Grad        = NULL;
  AD_InputIndex     = NULL;           AD_OutputIndex      = NULL;

  /*--- Volume (0 -> Vol_nP1, 1-> Vol_n, 2 -> Vol_nM1 ) and coordinates of the control volume ---*/

  if (config->GetTime_Marching() == NO) {
    Volume = new su2double[1];
    Volume[0] = 0.0;
  }
  else {
    Volume = new su2double[3];
    Volume[0] = 0.0;
    Volume[1] = 0.0;
    Volume[2] = 0.0;
  }

  Coord = new su2double[nDim];

  if(config->GetAD_Mode() && config->GetMultizone_Problem()) {
    AD_InputIndex   = new int[nDim];
    AD_OutputIndex  = new int[nDim];
  }

  /*--- Indicator if the control volume has been agglomerated ---*/
  Parent_CV   = 0;
  Agglomerate = false;

  /*--- Flip the normal orientation ---*/
  Flip_Orientation = false;

  /*--- Indicator if the point is going to be moved in a volumetric deformation ---*/
  Move = true;

  /*--- Identify boundaries, physical boundaries (not send-receive
  condition), detect if an element belong to the domain or it must
  be computed with other processor  ---*/
  Domain           = true;
  Boundary         = false;
  SolidBoundary    = false;
  PhysicalBoundary = false;
  PeriodicBoundary = false;

  /*--- Set the global index in the parallel simulation ---*/
  GlobalIndex = val_globalindex;

  /*--- Set the color for mesh partitioning ---*/
  color = 0;

  /*--- For smoothing the numerical grid coordinates ---*/
  if ( config->GetSmoothNumGrid() ) {
    Coord_Old = new su2double[nDim];
    Coord_Sum = new su2double[nDim];
  }

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  bool dynamic_grid = config->GetDynamic_Grid();

  /*--- Grid velocity gradients are only needed for the continuous adjoint ---*/
  bool continuous_adjoint = (config->GetKind_Solver() == ADJ_EULER ||
                             config->GetKind_Solver() == ADJ_NAVIER_STOKES ||
                             config->GetKind_Solver() == ADJ_RANS);

  /*--- Storage of grid velocities for dynamic meshes ---*/

  if ( dynamic_grid ) {
    GridVel  = new su2double[nDim];

    for (iDim = 0; iDim < nDim; iDim++)
      GridVel[iDim] = 0.0;

    if (continuous_adjoint){
    /*--- Gradient of the grid velocity ---*/
      GridVel_Grad = new su2double*[nDim];

      for (iDim = 0; iDim < nDim; iDim++) {
        GridVel_Grad[iDim] = new su2double[nDim];
        for (jDim = 0; jDim < nDim; jDim++)
          GridVel_Grad[iDim][jDim] = 0.0;
      }
    }

    /*--- Structures for storing old node coordinates for computing grid
    velocities via finite differencing with dynamically deforming meshes. ---*/
    /*--- In the case of deformable mesh solver, these coordinates are stored as solutions to the mesh problem ---*/
    if ( config->GetGrid_Movement() && (config->GetTime_Marching() != NO)) {
      Coord_p1 = new su2double[nDim];
      Coord_n  = new su2double[nDim];
      Coord_n1 = new su2double[nDim];
      Coord_Old = new su2double[nDim];
    }
  }

  /*--- Intialize the value of the curvature ---*/
  Curvature = 0.0;

  /*--- Intialize the value of the periodic volume. ---*/
  Periodic_Volume = 0.0;

  /*--- Init walldistance ---*/

  Wall_Distance = 0.0;
}

CPoint::CPoint(su2double val_coord_0, su2double val_coord_1, unsigned long val_globalindex, CConfig *config) : CDualGrid(2) {

  unsigned short iDim, jDim;

  /*--- Element, point and edge structures initialization ---*/
  Elem.clear();  nElem  = 0;
  Point.clear(); nPoint = 0;
  Edge.clear();

  Volume            = NULL;           Vertex              = NULL;
  Coord             = NULL;           Coord_Old           = NULL;            Coord_Sum  = NULL;
  Coord_n           = NULL;           Coord_n1            = NULL;            Coord_p1   = NULL;
  GridVel           = NULL;           GridVel_Grad        = NULL;
  AD_InputIndex     = NULL;           AD_OutputIndex      = NULL;

  /*--- Volume (0 -> Vol_nP1, 1-> Vol_n, 2 -> Vol_nM1 ) and coordinates of the control volume ---*/

  if (config->GetTime_Marching() == NO) {
    Volume = new su2double[1];
    Volume[0] = 0.0;
  }
  else{
    Volume = new su2double[3];
    Volume[0] = 0.0;
    Volume[1] = 0.0;
    Volume[2] = 0.0;
  }

  Coord    = new su2double[nDim];
  Coord[0] = val_coord_0;
  Coord[1] = val_coord_1;

  if(config->GetAD_Mode() && config->GetMultizone_Problem()) {
    AD_InputIndex   = new int[nDim];
    AD_OutputIndex  = new int[nDim];
  }

  /*--- Indicator if the control volume has been agglomerated ---*/
  Parent_CV   = 0;
  Agglomerate = false;

  /*--- Flip the normal orientation ---*/
  Flip_Orientation = false;

  /*--- Indicator if the point is going to be moved in a volumetric deformation ---*/
  Move = true;

  /*--- Identify boundaries, physical boundaries (not send-receive
  condition), detect if an element belong to the domain or it must
  be computed with other processor  ---*/
  Domain           = true;
  Boundary         = false;
  SolidBoundary    = false;
  PhysicalBoundary = false;
  PeriodicBoundary = false;

  /*--- Set the color for mesh partitioning ---*/
  color = 0;

  /*--- Set the global index in the parallel simulation ---*/
  GlobalIndex = val_globalindex;

  /*--- For smoothing the numerical grid coordinates ---*/
  if ( config->GetSmoothNumGrid() ) {
    Coord_Old = new su2double[nDim];
    Coord_Sum = new su2double[nDim];
  }

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  bool dynamic_grid = config->GetDynamic_Grid();

  /*--- Grid velocity gradients are only needed for the continuous adjoint ---*/
  bool continuous_adjoint = (config->GetKind_Solver() == ADJ_EULER ||
                             config->GetKind_Solver() == ADJ_NAVIER_STOKES ||
                             config->GetKind_Solver() == ADJ_RANS);

  /*--- Storage of grid velocities for dynamic meshes ---*/
  if ( dynamic_grid ) {
    GridVel  = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      GridVel[iDim] = 0.0;

    if (continuous_adjoint){
    /*--- Gradient of the grid velocity ---*/
      GridVel_Grad = new su2double*[nDim];
      for (iDim = 0; iDim < nDim; iDim++) {
        GridVel_Grad[iDim] = new su2double[nDim];
        for (jDim = 0; jDim < nDim; jDim++)
          GridVel_Grad[iDim][jDim] = 0.0;
      }
    }

    /*--- Structures for storing old node coordinates for computing grid
    velocities via finite differencing with dynamically deforming meshes. ---*/
    /*--- In the case of deformable mesh solver, these coordinates are stored as solutions to the mesh problem ---*/
    if ( config->GetGrid_Movement() && (config->GetTime_Marching() != NO)) {
      Coord_p1 = new su2double[nDim];
      Coord_n  = new su2double[nDim];
      Coord_n1 = new su2double[nDim];
      Coord_Old = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim ++) {
        Coord_p1[iDim] = Coord[iDim];
        Coord_n[iDim]  = Coord[iDim];
        Coord_n1[iDim] = Coord[iDim];
      }
    }
  }

  /*--- Intialize the value of the curvature ---*/
  Curvature = 0.0;

  /*--- Intialize the value of the periodic volume. ---*/
  Periodic_Volume = 0.0;

}

CPoint::CPoint(su2double val_coord_0, su2double val_coord_1, su2double val_coord_2, unsigned long val_globalindex, CConfig *config) : CDualGrid(3) {

  unsigned short iDim, jDim;

  /*--- Element, point and edge structures initialization ---*/
  Elem.clear();  nElem  = 0;
  Point.clear(); nPoint = 0;
  Edge.clear();

  Volume            = NULL;           Vertex              = NULL;
  Coord             = NULL;           Coord_Old           = NULL;            Coord_Sum  = NULL;
  Coord_n           = NULL;           Coord_n1            = NULL;            Coord_p1   = NULL;
  GridVel           = NULL;           GridVel_Grad        = NULL;
  AD_InputIndex     = NULL;           AD_OutputIndex      = NULL;

  /*--- Volume (0 -> Vol_nP1, 1-> Vol_n, 2 -> Vol_nM1 ) and coordinates of the control volume ---*/
  if ( config->GetTime_Marching() == NO ) {
    Volume = new su2double[1];
    Volume[0] = 0.0;
  }
  else{
    Volume = new su2double[3];
    Volume[0] = 0.0;
    Volume[1] = 0.0;
    Volume[2] = 0.0;
  }

  Coord    = new su2double[nDim];
  Coord[0] = val_coord_0;
  Coord[1] = val_coord_1;
  Coord[2] = val_coord_2;

  if(config->GetAD_Mode() && config->GetMultizone_Problem()) {
    AD_InputIndex   = new int[nDim];
    AD_OutputIndex  = new int[nDim];
  }

  /*--- Indicator if the control volume has been agglomerated ---*/
  Parent_CV = 0;
  Agglomerate = false;

  /*--- Indicator if the point is going to be moved in a volumetric deformation ---*/
  Move = true;

  /*--- Flip the normal orientation ---*/
  Flip_Orientation = false;

  /*--- Identify boundaries, physical boundaries (not send-receive
  condition), detect if an element belong to the domain or it must
  be computed with other processor  ---*/
  Domain           = true;
  Boundary         = false;
  SolidBoundary    = false;
  PhysicalBoundary = false;
  PeriodicBoundary = false;

  /*--- Set the color for mesh partitioning ---*/
  color = 0;

  /*--- Set the global index in the parallel simulation ---*/
  GlobalIndex = val_globalindex;

  /*--- For smoothing the numerical grid coordinates ---*/
  if (config->GetSmoothNumGrid()) {
    Coord_Old = new su2double[nDim];
    Coord_Sum = new su2double[nDim];
  }

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  bool dynamic_grid = config->GetDynamic_Grid();

  /*--- Grid velocity gradients are only needed for the continuous adjoint ---*/
  bool continuous_adjoint = (config->GetKind_Solver() == ADJ_EULER ||
                             config->GetKind_Solver() == ADJ_NAVIER_STOKES ||
                             config->GetKind_Solver() == ADJ_RANS);

  /*--- Storage of grid velocities for dynamic meshes ---*/

  if (dynamic_grid) {
    GridVel = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim ++)
      GridVel[iDim] = 0.0;

    if (continuous_adjoint){
    /*--- Gradient of the grid velocity ---*/
      GridVel_Grad = new su2double*[nDim];
      for (iDim = 0; iDim < nDim; iDim++) {
        GridVel_Grad[iDim] = new su2double[nDim];
        for (jDim = 0; jDim < nDim; jDim++)
          GridVel_Grad[iDim][jDim] = 0.0;
      }
    }

    /*--- Structures for storing old node coordinates for computing grid
    velocities via finite differencing with dynamically deforming meshes. ---*/
    /*--- In the case of deformable mesh solver, these coordinates are stored as solutions to the mesh problem ---*/
    if ( config->GetGrid_Movement() && (config->GetTime_Marching() != NO)) {
      Coord_p1 = new su2double[nDim];
      Coord_n  = new su2double[nDim];
      Coord_n1 = new su2double[nDim];
      Coord_Old = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim ++) {
        Coord_p1[iDim] = Coord[iDim];
        Coord_n[iDim]  = Coord[iDim];
        Coord_n1[iDim] = Coord[iDim];
      }
    }
  }

  /*--- Intialize the value of the curvature ---*/
  Curvature = 0.0;

  /*--- Intialize the value of the periodic volume. ---*/
  Periodic_Volume = 0.0;

}

CPoint::~CPoint() {

  if (Vertex       != NULL && Boundary) delete[] Vertex;
  if (Volume       != NULL) delete[] Volume;
  if (Coord        != NULL) delete[] Coord;
  if (Coord_Old    != NULL) delete[] Coord_Old;
  if (Coord_Sum    != NULL) delete[] Coord_Sum;
  if (Coord_n      != NULL) delete[] Coord_n;
  if (Coord_n1     != NULL) delete[] Coord_n1;
  if (Coord_p1     != NULL) delete[] Coord_p1;
  if (GridVel      != NULL) delete[] GridVel;
  if (GridVel_Grad != NULL) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      delete [] GridVel_Grad[iDim];
    delete [] GridVel_Grad;
  }
  if (AD_InputIndex  != NULL) delete[] AD_InputIndex;
  if (AD_OutputIndex != NULL) delete[] AD_OutputIndex;
 }

void CPoint::SetPoint(unsigned long val_point) {

  unsigned short iPoint;
  bool new_point;

  /*--- Look for the point in the list ---*/
  new_point = true;
  for (iPoint = 0; iPoint < GetnPoint(); iPoint++)
  if (Point[iPoint] == val_point) {
    new_point = false;
    break;
  }

  /*--- Store the point structure and dimensionalizate edge structure ---*/
  if (new_point) {
    Point.push_back(val_point);
    Edge.push_back(-1);
    nPoint = Point.size();
  }

}

void CPoint::SetBoundary(unsigned short val_nmarker) {

  unsigned short imarker;

  /*--- To be sure that we are not goint to initializate twice the same vertex ---*/
  if (!Boundary) {
    Vertex = new long[val_nmarker];

    /*--- The initialization is made with -1 ---*/
    for (imarker = 0; imarker < val_nmarker; imarker++)
      Vertex[imarker] = -1;
  }
  Boundary = true;

}

void CPoint::SetIndex(bool input) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    if(input) {
      AD::SetIndex(AD_InputIndex[iDim], Coord[iDim]);
    }
    else {
      AD::SetIndex(AD_OutputIndex[iDim], Coord[iDim]);
    }
  }
}

void CPoint::SetAdjointSolution(const su2double *adj_sol) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    AD::SetDerivative(AD_OutputIndex[iDim], SU2_TYPE::GetValue(adj_sol[iDim]));
  }
}

su2double CPoint::GetAdjointSolution(unsigned short iDim) const {
  return AD::GetDerivative(AD_InputIndex[iDim]);
}
