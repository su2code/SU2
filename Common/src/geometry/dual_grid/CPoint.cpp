/*!
 * \file CPoint.cpp
 * \brief Main classes for defining the points of the dual grid
 * \author F. Palacios, T. Economon
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

#include "../../../include/geometry/dual_grid/CPoint.hpp"
#include "../../../include/CConfig.hpp"


CPoint::CPoint(unsigned long npoint, unsigned long ndim, const CConfig *config) {

  nDim = ndim;

  /*--- Element, point and edge structures initialization. ---*/

  Elem.resize(npoint);  nElem.resize(npoint) = 0;
  Point.resize(npoint); nPoint.resize(npoint) = 0;
  Edge.resize(npoint);
  Vertex.resize(npoint);

  /*--- Coordinates and volumes. ---*/

  Coord.resize(npoint,nDim) = su2double(0.0);

  Volume.resize(npoint) = su2double(0.0);
  Periodic_Volume.resize(npoint) = su2double(0.0);

  if (config->GetTime_Marching() != NO) {
    Volume_n.resize(npoint) = su2double(0.0);
    Volume_nM1.resize(npoint) = su2double(0.0);
  }

  if(config->GetAD_Mode() && config->GetMultizone_Problem()) {
    AD_InputIndex.resize(npoint,nDim) = 0;
    AD_OutputIndex.resize(npoint,nDim) = 0;
  }

  /*--- Indicator if the control volume has been agglomerated. ---*/
  Parent_CV.resize(npoint) = 0;
  Agglomerate.resize(npoint) = false;
  Agglomerate_Indirect.resize(npoint) = false;

  /*--- Flip the normal orientation ---*/
  Flip_Orientation.resize(npoint) = false;

  /*--- Indicator if the point is going to be moved in a volumetric deformation. ---*/
  Move.resize(npoint) = true;

  /*--- Identify boundaries, physical boundaries (not send-receive
   *    condition), detect if an element belong to the domain or it
   *    must be computed with other processor. ---*/
  Domain.resize(npoint)           = true;
  Boundary.resize(npoint)         = false;
  SolidBoundary.resize(npoint)    = false;
  PhysicalBoundary.resize(npoint) = false;
  PeriodicBoundary.resize(npoint) = false;

  /*--- Set the global index in the parallel simulation. ---*/
  GlobalIndex.resize(npoint) = 0;

  /*--- Set the color for mesh partitioning. ---*/
  Color.resize(npoint) = 0;

  /*--- For smoothing the numerical grid coordinates ---*/
  if (config->GetSmoothNumGrid()) {
    Coord_Old.resize(npoint,nDim) = su2double(0.0);
    Coord_Sum.resize(npoint,nDim) = su2double(0.0);
  }

  /*--- Storage of grid velocities for dynamic meshes. ---*/

  if (config->GetDynamic_Grid()) {
    GridVel.resize(npoint,nDim) = su2double(0.0);

    /*--- Grid velocity gradients are needed for the continuous adjoint. ---*/
    if (config->GetContinuous_Adjoint())
      GridVel_Grad.resize(npoint,nDim,nDim,0.0);

    /*--- Structures for storing old node coordinates for computing grid
     *    velocities via finite differencing with dynamically deforming meshes. ---*/
    /*--- In the case of CMeshSolver, these coordinates are stored as solutions to the mesh problem. ---*/
    if ( config->GetGrid_Movement() && (config->GetTime_Marching() != NO)) {
      Coord_n.resize(npoint,nDim) = su2double(0.0);
      Coord_p1.resize(npoint,nDim) = su2double(0.0);
      Coord_n1.resize(npoint,nDim) = su2double(0.0);
      Coord_Old.resize(npoint,nDim) = su2double(0.0);
    }
  }

  Curvature.resize(npoint) = su2double(0.0);
  Wall_Distance.resize(npoint) = su2double(0.0);

}

void CPoint::SetPoint(unsigned long iPoint, unsigned long point) {

  bool new_point = true;

  for (auto iPt = 0ul; iPt < nPoint(iPoint); iPt++) {
    if (Point[iPoint][iPt] == point) {
      new_point = false;
      break;
    }
  }

  if (new_point) {
    Point[iPoint].push_back(point);
    Edge[iPoint].push_back(-1);
    nPoint(iPoint) = Point[iPoint].size();
  }
}

void CPoint::SetIndex(unsigned long iPoint, bool input) {
  for (unsigned long iDim = 0; iDim < nDim; iDim++) {
    if(input) {
      AD::SetIndex(AD_InputIndex(iPoint,iDim), Coord(iPoint,iDim));
    }
    else {
      AD::SetIndex(AD_OutputIndex(iPoint,iDim), Coord(iPoint,iDim));
    }
  }
}

//CPoint::CPoint(su2double val_coord_0, su2double val_coord_1, unsigned long val_globalindex, CConfig *config) : CDualGrid(2) {
//
//  unsigned short iDim, jDim;
//
//  /*--- Element, point and edge structures initialization ---*/
//  Elem.clear();  nElem  = 0;
//  Point.clear(); nPoint = 0;
//  Edge.clear();
//
//  Volume            = NULL;           Vertex              = NULL;
//  Coord             = NULL;           Coord_Old           = NULL;            Coord_Sum  = NULL;
//  Coord_n           = NULL;           Coord_n1            = NULL;            Coord_p1   = NULL;
//  GridVel           = NULL;           GridVel_Grad        = NULL;
//  AD_InputIndex     = NULL;           AD_OutputIndex      = NULL;
//
//  /*--- Volume (0 -> Vol_nP1, 1-> Vol_n, 2 -> Vol_nM1 ) and coordinates of the control volume ---*/
//
//  if (config->GetTime_Marching() == NO) {
//    Volume = new su2double[1];
//    Volume[0] = 0.0;
//  }
//  else{
//    Volume = new su2double[3];
//    Volume[0] = 0.0;
//    Volume[1] = 0.0;
//    Volume[2] = 0.0;
//  }
//
//  Coord    = new su2double[nDim];
//  Coord[0] = val_coord_0;
//  Coord[1] = val_coord_1;
//
//  if(config->GetAD_Mode() && config->GetMultizone_Problem()) {
//    AD_InputIndex   = new int[nDim];
//    AD_OutputIndex  = new int[nDim];
//  }
//
//  /*--- Indicator if the control volume has been agglomerated ---*/
//  Parent_CV   = 0;
//  Agglomerate = false;
//
//  /*--- Flip the normal orientation ---*/
//  Flip_Orientation = false;
//
//  /*--- Indicator if the point is going to be moved in a volumetric deformation ---*/
//  Move = true;
//
//  /*--- Identify boundaries, physical boundaries (not send-receive
//  condition), detect if an element belong to the domain or it must
//  be computed with other processor  ---*/
//  Domain           = true;
//  Boundary         = false;
//  SolidBoundary    = false;
//  PhysicalBoundary = false;
//  PeriodicBoundary = false;
//
//  /*--- Set the color for mesh partitioning ---*/
//  color = 0;
//
//  /*--- Set the global index in the parallel simulation ---*/
//  GlobalIndex = val_globalindex;
//
//  /*--- For smoothing the numerical grid coordinates ---*/
//  if ( config->GetSmoothNumGrid() ) {
//    Coord_Old = new su2double[nDim];
//    Coord_Sum = new su2double[nDim];
//  }
//
//  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
//  bool dynamic_grid = config->GetDynamic_Grid();
//
//  /*--- Grid velocity gradients are only needed for the continuous adjoint ---*/
//  bool continuous_adjoint = (config->GetKind_Solver() == ADJ_EULER ||
//                             config->GetKind_Solver() == ADJ_NAVIER_STOKES ||
//                             config->GetKind_Solver() == ADJ_RANS);
//
//  /*--- Storage of grid velocities for dynamic meshes ---*/
//  if ( dynamic_grid ) {
//    GridVel  = new su2double[nDim];
//    for (iDim = 0; iDim < nDim; iDim++)
//      GridVel[iDim] = 0.0;
//
//    if (continuous_adjoint){
//    /*--- Gradient of the grid velocity ---*/
//      GridVel_Grad = new su2double*[nDim];
//      for (iDim = 0; iDim < nDim; iDim++) {
//        GridVel_Grad[iDim] = new su2double[nDim];
//        for (jDim = 0; jDim < nDim; jDim++)
//          GridVel_Grad[iDim][jDim] = 0.0;
//      }
//    }
//
//    /*--- Structures for storing old node coordinates for computing grid
//    velocities via finite differencing with dynamically deforming meshes. ---*/
//    /*--- In the case of deformable mesh solver, these coordinates are stored as solutions to the mesh problem ---*/
//    if ( config->GetGrid_Movement() && (config->GetTime_Marching() != NO)) {
//      Coord_p1 = new su2double[nDim];
//      Coord_n  = new su2double[nDim];
//      Coord_n1 = new su2double[nDim];
//      Coord_Old = new su2double[nDim];
//      for (iDim = 0; iDim < nDim; iDim ++) {
//        Coord_p1[iDim] = Coord[iDim];
//        Coord_n[iDim]  = Coord[iDim];
//        Coord_n1[iDim] = Coord[iDim];
//      }
//    }
//  }
//
//  /*--- Intialize the value of the curvature ---*/
//  Curvature = 0.0;
//
//  /*--- Intialize the value of the periodic volume. ---*/
//  Periodic_Volume = 0.0;
//
//}
//
//CPoint::CPoint(su2double val_coord_0, su2double val_coord_1, su2double val_coord_2, unsigned long val_globalindex, CConfig *config) : CDualGrid(3) {
//
//  unsigned short iDim, jDim;
//
//  /*--- Element, point and edge structures initialization ---*/
//  Elem.clear();  nElem  = 0;
//  Point.clear(); nPoint = 0;
//  Edge.clear();
//
//  Volume            = NULL;           Vertex              = NULL;
//  Coord             = NULL;           Coord_Old           = NULL;            Coord_Sum  = NULL;
//  Coord_n           = NULL;           Coord_n1            = NULL;            Coord_p1   = NULL;
//  GridVel           = NULL;           GridVel_Grad        = NULL;
//  AD_InputIndex     = NULL;           AD_OutputIndex      = NULL;
//
//  /*--- Volume (0 -> Vol_nP1, 1-> Vol_n, 2 -> Vol_nM1 ) and coordinates of the control volume ---*/
//  if ( config->GetTime_Marching() == NO ) {
//    Volume = new su2double[1];
//    Volume[0] = 0.0;
//  }
//  else{
//    Volume = new su2double[3];
//    Volume[0] = 0.0;
//    Volume[1] = 0.0;
//    Volume[2] = 0.0;
//  }
//
//  Coord    = new su2double[nDim];
//  Coord[0] = val_coord_0;
//  Coord[1] = val_coord_1;
//  Coord[2] = val_coord_2;
//
//  if(config->GetAD_Mode() && config->GetMultizone_Problem()) {
//    AD_InputIndex   = new int[nDim];
//    AD_OutputIndex  = new int[nDim];
//  }
//
//  /*--- Indicator if the control volume has been agglomerated ---*/
//  Parent_CV = 0;
//  Agglomerate = false;
//
//  /*--- Indicator if the point is going to be moved in a volumetric deformation ---*/
//  Move = true;
//
//  /*--- Flip the normal orientation ---*/
//  Flip_Orientation = false;
//
//  /*--- Identify boundaries, physical boundaries (not send-receive
//  condition), detect if an element belong to the domain or it must
//  be computed with other processor  ---*/
//  Domain           = true;
//  Boundary         = false;
//  SolidBoundary    = false;
//  PhysicalBoundary = false;
//  PeriodicBoundary = false;
//
//  /*--- Set the color for mesh partitioning ---*/
//  color = 0;
//
//  /*--- Set the global index in the parallel simulation ---*/
//  GlobalIndex = val_globalindex;
//
//  /*--- For smoothing the numerical grid coordinates ---*/
//  if (config->GetSmoothNumGrid()) {
//    Coord_Old = new su2double[nDim];
//    Coord_Sum = new su2double[nDim];
//  }
//
//  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
//  bool dynamic_grid = config->GetDynamic_Grid();
//
//  /*--- Grid velocity gradients are only needed for the continuous adjoint ---*/
//  bool continuous_adjoint = (config->GetKind_Solver() == ADJ_EULER ||
//                             config->GetKind_Solver() == ADJ_NAVIER_STOKES ||
//                             config->GetKind_Solver() == ADJ_RANS);
//
//  /*--- Storage of grid velocities for dynamic meshes ---*/
//
//  if (dynamic_grid) {
//    GridVel = new su2double[nDim];
//    for (iDim = 0; iDim < nDim; iDim ++)
//      GridVel[iDim] = 0.0;
//
//    if (continuous_adjoint){
//    /*--- Gradient of the grid velocity ---*/
//      GridVel_Grad = new su2double*[nDim];
//      for (iDim = 0; iDim < nDim; iDim++) {
//        GridVel_Grad[iDim] = new su2double[nDim];
//        for (jDim = 0; jDim < nDim; jDim++)
//          GridVel_Grad[iDim][jDim] = 0.0;
//      }
//    }
//
//    /*--- Structures for storing old node coordinates for computing grid
//    velocities via finite differencing with dynamically deforming meshes. ---*/
//    /*--- In the case of deformable mesh solver, these coordinates are stored as solutions to the mesh problem ---*/
//    if ( config->GetGrid_Movement() && (config->GetTime_Marching() != NO)) {
//      Coord_p1 = new su2double[nDim];
//      Coord_n  = new su2double[nDim];
//      Coord_n1 = new su2double[nDim];
//      Coord_Old = new su2double[nDim];
//      for (iDim = 0; iDim < nDim; iDim ++) {
//        Coord_p1[iDim] = Coord[iDim];
//        Coord_n[iDim]  = Coord[iDim];
//        Coord_n1[iDim] = Coord[iDim];
//      }
//    }
//  }
//
//  /*--- Intialize the value of the curvature ---*/
//  Curvature = 0.0;
//
//  /*--- Intialize the value of the periodic volume. ---*/
//  Periodic_Volume = 0.0;
//
//}
