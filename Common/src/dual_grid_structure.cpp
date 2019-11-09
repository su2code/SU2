/*!
 * \file dual_grid_structure.cpp
 * \brief Main classes for defining the dual grid
 * \author F. Palacios, T. Economon
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

#include "../include/dual_grid_structure.hpp"

unsigned short CDualGrid::nDim = 0;

CDualGrid::CDualGrid(unsigned short val_nDim) { nDim = val_nDim;}

CDualGrid::~CDualGrid() {}

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

su2double CPoint::GetAdjointSolution(unsigned short iDim) {
  return AD::GetDerivative(AD_InputIndex[iDim]);
}

CEdge::CEdge(unsigned long val_iPoint, unsigned long val_jPoint, unsigned short val_nDim) : CDualGrid(val_nDim) {
    
  unsigned short iDim;
    
  /*--- Pointers initialization ---*/
  Coord_CG = NULL;
  Normal   = NULL;
  Nodes    = NULL;

  /*--- Allocate center of gravity coordinates, nodes, and face normal ---*/
  Coord_CG = new su2double [nDim];
  Normal   = new su2double [nDim];
  Nodes    = new unsigned long[2];

  /*--- Initializate the structure ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    Coord_CG[iDim] = 0.0;
    Normal[iDim]   = 0.0;
  }

  Nodes[0] = val_iPoint; 
  Nodes[1] = val_jPoint;

}

CEdge::~CEdge() {
  
  if (Coord_CG != NULL) delete[] Coord_CG;
  if (Normal   != NULL) delete[] Normal;
  if (Nodes    != NULL) delete[] Nodes;
  
}

void CEdge::SetCoord_CG(su2double **val_coord) {

  unsigned short iDim, iNode;

  for (iDim = 0; iDim < nDim; iDim++) {
    Coord_CG[iDim] = 0.0;
    for (iNode = 0; iNode < 2;  iNode++)
      Coord_CG[iDim] += val_coord[iNode][iDim] / 2.0;
  }

}

su2double CEdge::GetVolume(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG, su2double *val_coord_Point) {

  unsigned short iDim;
  su2double vec_a[3] = {0.0,0.0,0.0}, vec_b[3] = {0.0,0.0,0.0}, vec_c[3] = {0.0,0.0,0.0}, vec_d[3] = {0.0,0.0,0.0}, Local_Volume;

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord_Edge_CG, nDim);
  AD::SetPreaccIn(val_coord_Elem_CG, nDim);
  AD::SetPreaccIn(val_coord_FaceElem_CG, nDim);
  AD::SetPreaccIn(val_coord_Point, nDim);

  for (iDim = 0; iDim < nDim; iDim++) {
    vec_a[iDim] = val_coord_Edge_CG[iDim]     - val_coord_Point[iDim];
    vec_b[iDim] = val_coord_FaceElem_CG[iDim] - val_coord_Point[iDim];
    vec_c[iDim] = val_coord_Elem_CG[iDim]     - val_coord_Point[iDim];
  }

  vec_d[0] =   vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1];
  vec_d[1] = -(vec_a[0] * vec_b[2] - vec_a[2] * vec_b[0]);
  vec_d[2] =   vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0];

  Local_Volume = fabs( vec_c[0] * vec_d[0] + vec_c[1] * vec_d[1] + vec_c[2] * vec_d[2] ) / 6.0;

  AD::SetPreaccOut(Local_Volume);
  AD::EndPreacc();

  return Local_Volume;

}

su2double CEdge::GetVolume(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG, su2double *val_coord_Point) {

  unsigned short iDim;
  su2double vec_a[2] = {0.0,0.0}, vec_b[2] = {0.0,0.0}, Local_Volume;

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord_Edge_CG, nDim);
  AD::SetPreaccIn(val_coord_Elem_CG, nDim);
  AD::SetPreaccIn(val_coord_Point, nDim);

  for (iDim = 0; iDim < nDim; iDim++) {
    vec_a[iDim] = val_coord_Elem_CG[iDim] - val_coord_Point[iDim];
    vec_b[iDim] = val_coord_Edge_CG[iDim] - val_coord_Point[iDim];
  }

  Local_Volume = 0.5 * fabs( vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0] );

  AD::SetPreaccOut(Local_Volume);
  AD::EndPreacc();

  return Local_Volume;

}

void CEdge::SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG) {

  unsigned short iDim;
  su2double vec_a[3] = {0.0,0.0,0.0}, vec_b[3] = {0.0,0.0,0.0}, Dim_Normal[3];

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord_Edge_CG, nDim);
  AD::SetPreaccIn(val_coord_Elem_CG, nDim);
  AD::SetPreaccIn(val_coord_FaceElem_CG, nDim);
  AD::SetPreaccIn(Normal, nDim);

  for (iDim = 0; iDim < nDim; iDim++) {
    vec_a[iDim] = val_coord_Elem_CG[iDim]-val_coord_Edge_CG[iDim];
    vec_b[iDim] = val_coord_FaceElem_CG[iDim]-val_coord_Edge_CG[iDim];
  }

  Dim_Normal[0] =  0.5 * ( vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1] );
  Dim_Normal[1] = -0.5 * ( vec_a[0] * vec_b[2] - vec_a[2] * vec_b[0] );
  Dim_Normal[2] =  0.5 * ( vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0] );

  Normal[0] += Dim_Normal[0]; 
  Normal[1] += Dim_Normal[1];       
  Normal[2] += Dim_Normal[2];

  AD::SetPreaccOut(Normal, nDim);
  AD::EndPreacc();
}

void CEdge::SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG) {

  su2double Dim_Normal[2];

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord_Elem_CG, nDim);
  AD::SetPreaccIn(val_coord_Edge_CG, nDim);
  AD::SetPreaccIn(Normal, nDim);

  Dim_Normal[0] =   val_coord_Elem_CG[1] - val_coord_Edge_CG[1];
  Dim_Normal[1] = -(val_coord_Elem_CG[0] - val_coord_Edge_CG[0]);

  Normal[0] += Dim_Normal[0]; 
  Normal[1] += Dim_Normal[1];

  AD::SetPreaccOut(Normal, nDim);
  AD::EndPreacc();

}

CVertex::CVertex(unsigned long val_point, unsigned short val_nDim) : CDualGrid(val_nDim) {

  unsigned short iDim;

  /*--- Set periodic points to zero ---*/
  
  PeriodicPoint[0] = -1; PeriodicPoint[1] = -1; PeriodicPoint[2] = -1;
  PeriodicPoint[3] = -1; PeriodicPoint[4] = -1;
  
  /*--- Identify the points at the perimeter of the actuatrod disk ---*/
  
  ActDisk_Perimeter = false;

  /*--- Pointers initialization ---*/
  
  Nodes  = NULL;
  Normal = NULL;

  /*--- Allocate node, and face normal ---*/
  
  Nodes  = new unsigned long[1]; 
  Normal = new su2double [nDim];

  /*--- Initializate the structure ---*/
  
  Nodes[0] = val_point;
  for (iDim = 0; iDim < nDim; iDim ++) 
    Normal[iDim] = 0.0;

  /*--- Set to zero the variation of the coordinates ---*/
  
  VarCoord[0] = 0.0; 
  VarCoord[1] = 0.0; 
  VarCoord[2] = 0.0;

  /*--- Set to NULL variation of the rotation  ---*/
  
  VarRot = NULL;

  /*--- Set to NULL donor arrays for interpolation ---*/
  
  Donor_Points  = NULL;
  Donor_Proc    = NULL;
  Donor_Coeff   = NULL;
  nDonor_Points = 1;

}

CVertex::~CVertex() {
  
  if (Normal != NULL) delete[] Normal;
  if (Nodes  != NULL) delete[] Nodes;

  /*---  donor arrays for interpolation ---*/
  
  if (VarRot       != NULL) delete[] VarRot;
  if (Donor_Coeff  != NULL) delete[] Donor_Coeff;
  if (Donor_Proc   != NULL) delete[] Donor_Proc;
  if (Donor_Points != NULL) delete[] Donor_Points;

}

void CVertex::SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG) {

  su2double vec_a[3] = {0.0,0.0,0.0}, vec_b[3] = {0.0,0.0,0.0};
  unsigned short iDim;

  assert(nDim == 3);

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord_Edge_CG, nDim);
  AD::SetPreaccIn(val_coord_Elem_CG, nDim);
  AD::SetPreaccIn(val_coord_FaceElem_CG, nDim);
  AD::SetPreaccIn(Normal, nDim);

  for (iDim = 0; iDim < nDim; iDim++) {
    vec_a[iDim] = val_coord_Elem_CG[iDim]-val_coord_Edge_CG[iDim];
    vec_b[iDim] = val_coord_FaceElem_CG[iDim]-val_coord_Edge_CG[iDim];
  }

  Normal[0] += 0.5 * ( vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1]);
  Normal[1] -= 0.5 * ( vec_a[0] * vec_b[2] - vec_a[2] * vec_b[0]);
  Normal[2] += 0.5 * ( vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);

  AD::SetPreaccOut(Normal, nDim);
  AD::EndPreacc();

}

void CVertex::SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG) {

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord_Elem_CG, nDim);
  AD::SetPreaccIn(val_coord_Edge_CG, nDim);
  AD::SetPreaccIn(Normal, nDim);

  Normal[0] += val_coord_Elem_CG[1]-val_coord_Edge_CG[1];
  Normal[1] -= (val_coord_Elem_CG[0]-val_coord_Edge_CG[0]);

  AD::SetPreaccOut(Normal, nDim);
  AD::EndPreacc();
 
}

void CVertex::AddNormal(su2double *val_face_normal) {

  unsigned short i;
  for( i = 0; i < nDim; i++ )
    Normal[i] += val_face_normal[i]; 

}

void CVertex::Allocate_DonorInfo(void){
  
  if( Donor_Points != NULL )  delete [] Donor_Points;
  if( Donor_Proc   != NULL )  delete [] Donor_Proc;
  if( Donor_Coeff  != NULL )  delete [] Donor_Coeff;  
  
  Donor_Points = new unsigned long[nDonor_Points];
  Donor_Proc   = new unsigned long[nDonor_Points];
  Donor_Coeff  = new su2double[nDonor_Points];
}

CTurboVertex::CTurboVertex(unsigned long val_point, unsigned short val_nDim) : CVertex(val_point, val_nDim){
	unsigned short iDim;
 /*--- Pointers initialization ---*/
	TurboNormal = NULL;
	/*--- Allocate node, and face normal ---*/
	TurboNormal = new su2double [nDim];

	/*--- Initializate the structure ---*/
	for (iDim = 0; iDim < nDim; iDim ++) TurboNormal[iDim] = 0.0;

}

CTurboVertex::~CTurboVertex() {

	if (TurboNormal != NULL) delete [] TurboNormal;

}
