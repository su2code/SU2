/*!
 * \file primal_grid_structure.cpp
 * \brief Main classes for defining the primal grid elements (triangle, tetrahedra, etc.).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

#include "../include/primal_grid_structure.hpp"

unsigned short CPrimalGrid::nDim;

CPrimalGrid::CPrimalGrid(void) {
  
  /*--- Set the default values for the pointers ---*/
  Nodes = NULL;
	Neighbor_Elements = NULL;
	Coord_CG = NULL;
	Coord_FaceElems_CG = NULL;
  
}

CPrimalGrid::~CPrimalGrid() {

	if (Nodes != NULL) delete[] Nodes;
	if (Coord_CG != NULL) delete[] Coord_CG;
  if (Neighbor_Elements != NULL) delete[] Neighbor_Elements;
   
}

void CPrimalGrid::SetCG(double **val_coord) {
	unsigned short iDim, iNode, NodeFace, iFace;
	
	for (iDim = 0; iDim < nDim; iDim++) {
		Coord_CG[iDim] = 0.0;
		for (iNode = 0; iNode < GetnNodes();  iNode++)
			Coord_CG[iDim] += val_coord[iNode][iDim]/double(GetnNodes());
	}
	
	for (iFace = 0; iFace < GetnFaces();  iFace++)
		for (iDim = 0; iDim < nDim; iDim++) {
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
			for (iNode = 0; iNode < GetnNodesFace(iFace); iNode++) {
				NodeFace = GetFaces(iFace, iNode);
				Coord_FaceElems_CG[iFace][iDim] += val_coord[NodeFace][iDim]/double(GetnNodesFace(iFace));
			}
		}
}

void CPrimalGrid::GetAllNeighbor_Elements() {
	cout << "( ";
	for (unsigned short iFace = 0; iFace < GetnFaces(); iFace++)
	{
		cout << GetNeighbor_Elements(iFace) << ", ";
	}
	cout << ")"  << endl;
}

unsigned short CVertexMPI::nFaces = 0;

unsigned short CVertexMPI::nNodes = 1;

unsigned short CVertexMPI::nNeighbor_Elements = 0;

unsigned short CVertexMPI::VTK_Type = 1;

unsigned short CVertexMPI::maxNodesFace = 0;

CVertexMPI::CVertexMPI(unsigned long val_point, unsigned short val_ndim) : CPrimalGrid() {
	unsigned short iDim;
	
	/*--- Allocate CG coordinates ---*/
	nDim = val_ndim;
	Coord_CG = new double[nDim];
	for (iDim = 0; iDim < nDim; iDim++) Coord_CG[iDim] = 0.0;
	
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_point;
	
	/*--- By default, no rotation in the solution ---*/
	Rotation_Type = 0;
  
  /*--- By default, assume a single-zone mesh ---*/
  Matching_Zone = ZONE_0;
	
}

CVertexMPI::~CVertexMPI() {
  unsigned short iFaces;
  
    for (iFaces = 0; iFaces < nFaces; iFaces++)
      if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
    if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;

}

void CVertexMPI::Change_Orientation(void) { cout << "Not defined orientation change" << endl; }

unsigned short CLine::Faces[1][2]={{0,1}};

unsigned short CLine::Neighbor_Nodes[2][1]={{1},{0}};

unsigned short CLine::nNodesFace[1]={2};

unsigned short CLine::nNeighbor_Nodes[2]={1,1};

unsigned short CLine::nFaces = 1;

unsigned short CLine::nNodes = 2;

unsigned short CLine::nNeighbor_Elements = 1;

unsigned short CLine::VTK_Type = 3;

unsigned short CLine::maxNodesFace = 2;

CLine::CLine(unsigned long val_iPoint, unsigned long val_jPoint,
             unsigned short val_ndim) : CPrimalGrid() {
	unsigned short iDim, iFace;

	/*--- Allocate CG coordinates ---*/
	nDim = val_ndim;
	Coord_CG = new double[nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Coord_CG[iDim] = 0.0;
	Coord_FaceElems_CG = new double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_iPoint;
	Nodes[1] = val_jPoint;
  
}

CLine::~CLine() {
  unsigned short iFaces;
  
  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;

}

void CLine::Change_Orientation(void) {
	unsigned long iPoint, jPoint;
  
	iPoint = Nodes[0];
	jPoint = Nodes[1];
	Nodes[0] = jPoint;
	Nodes[1] = iPoint;
}

unsigned short CTriangle::Faces[3][2] = {{0,1},{1,2},{2,0}};

unsigned short CTriangle::Neighbor_Nodes[3][2] = {{1,2},{2,0},{0,1}};

unsigned short CTriangle::nNodesFace[3] = {2,2,2};

unsigned short CTriangle::nNeighbor_Nodes[3] = {2,2,2};

unsigned short CTriangle::nFaces = 3;

unsigned short CTriangle::nNodes = 3;

unsigned short CTriangle::nNeighbor_Elements = 3;

unsigned short CTriangle::VTK_Type = 5;

unsigned short CTriangle::maxNodesFace = 2;

CTriangle::CTriangle(unsigned long val_iPoint, unsigned long val_jPoint, 
					 unsigned long val_point_2, unsigned short val_nDim) : CPrimalGrid() {
	unsigned short iDim, iFace, iNeighbor_Elements;

	/*--- Allocate CG coordinates ---*/
	nDim = val_nDim;
	Coord_CG = new double[nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Coord_CG[iDim] = 0.0;
	Coord_FaceElems_CG = new double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_iPoint;
	Nodes[1] = val_jPoint;
	Nodes[2] = val_point_2;
	
	/*--- Allocate and define neighbor elements to a element ---*/
	nNeighbor_Elements = nFaces;
	Neighbor_Elements = new long[nNeighbor_Elements];
	for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
		Neighbor_Elements[iNeighbor_Elements]=-1;
	}
  
}

CTriangle::~CTriangle() {
  unsigned short iFaces;

  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;
  
}

void CTriangle::Change_Orientation(void) {
	unsigned long iPoint, Point_2;
	iPoint = Nodes[0];
	Point_2 = Nodes[2];
	Nodes[0] = Point_2;
	Nodes[2] = iPoint;
  
}

unsigned short CRectangle::Faces[4][2] = {{0,1},{1,2},{2,3},{3,0}};

unsigned short CRectangle::Neighbor_Nodes[4][2] = {{1,3},{2,0},{3,1},{0,2}};

unsigned short CRectangle::nNodesFace[4] = {2,2,2,2};

unsigned short CRectangle::nNeighbor_Nodes[4] = {2,2,2,2};

unsigned short CRectangle::nFaces = 4;

unsigned short CRectangle::nNodes = 4;

unsigned short CRectangle::nNeighbor_Elements = 4;

unsigned short CRectangle::VTK_Type = 9;

unsigned short CRectangle::maxNodesFace = 2;

CRectangle::CRectangle(unsigned long val_iPoint, unsigned long val_jPoint, 
					   unsigned long val_point_2, unsigned long val_point_3, unsigned short val_nDim) 
: CPrimalGrid() {
	unsigned short iDim, iFace, iNeighbor_Elements;

	/*--- Allocate CG coordinates ---*/
	nDim = val_nDim;
	Coord_CG = new double[nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Coord_CG[iDim] = 0.0;
	Coord_FaceElems_CG = new double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_iPoint;
	Nodes[1] = val_jPoint;
	Nodes[2] = val_point_2;
	Nodes[3] = val_point_3;
	
	
	nNeighbor_Elements = nFaces;
	Neighbor_Elements = new long[nNeighbor_Elements];
	for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
		Neighbor_Elements[iNeighbor_Elements]=-1;
	}
  
}

CRectangle::~CRectangle() {
  unsigned short iFaces;
  
  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;
  
}

void CRectangle::Change_Orientation(void) {
	unsigned long jPoint, Point_3;
	jPoint = Nodes[1];
	Point_3 = Nodes[3];
	Nodes[1] = Point_3;
	Nodes[3] = jPoint;
  
}

unsigned short CTetrahedron::Faces[4][3]={{0,2,1},{0,1,3},{0,3,2},{1,2,3}};

unsigned short CTetrahedron::Neighbor_Nodes[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};

unsigned short CTetrahedron::nNodesFace[4]={3,3,3,3};

unsigned short CTetrahedron::nNeighbor_Nodes[4]={3,3,3,3};

unsigned short CTetrahedron::nFaces = 4;

unsigned short CTetrahedron::nNodes = 4;

unsigned short CTetrahedron::nNeighbor_Elements = 4;

unsigned short CTetrahedron::VTK_Type = 10;

unsigned short CTetrahedron::maxNodesFace = 3;

CTetrahedron::CTetrahedron(unsigned long val_iPoint, unsigned long val_jPoint, 
						   unsigned long val_point_2, unsigned long val_point_3) : CPrimalGrid() {
	unsigned short iDim, iFace, iNeighbor_Elements;

	/*--- Allocate CG coordinates ---*/
	nDim = 3;
	Coord_CG = new double[nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Coord_CG[iDim] = 0.0;
	Coord_FaceElems_CG = new double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_iPoint;
	Nodes[1] = val_jPoint;
	Nodes[2] = val_point_2;
	Nodes[3] = val_point_3;
	
	/*--- Allocate and define neighbor elements to a element ---*/
	nNeighbor_Elements = nFaces;
	Neighbor_Elements = new long[nNeighbor_Elements];
	for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
		Neighbor_Elements[iNeighbor_Elements]=-1;
	}
  
}

CTetrahedron::~CTetrahedron() {
  unsigned short iFaces;
  
  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;
  
}

void CTetrahedron::Change_Orientation(void) {
	unsigned long iPoint, jPoint;
	iPoint = Nodes[0];
	jPoint = Nodes[1];
	Nodes[0] = jPoint;
	Nodes[1] = iPoint;
  
}

unsigned short CHexahedron::Faces[6][4] = {{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{0,3,2,1},{4,5,6,7}};

unsigned short CHexahedron::Neighbor_Nodes[8][3] = {{1,3,4},{0,2,5},{1,3,6},{0,2,7},{0,5,7},{4,6,1},{2,5,7},{4,3,6}};

unsigned short CHexahedron::nNodesFace[6] = {4,4,4,4,4,4};

unsigned short CHexahedron::nNeighbor_Nodes[8] = {3,3,3,3,3,3,3,3};

unsigned short CHexahedron::nFaces = 6;

unsigned short CHexahedron::nNodes = 8;

unsigned short CHexahedron::nNeighbor_Elements = 6;

unsigned short CHexahedron::VTK_Type = 12;

unsigned short CHexahedron::maxNodesFace = 4;

CHexahedron::CHexahedron(unsigned long val_iPoint, unsigned long val_jPoint, 
						 unsigned long val_point_2, unsigned long val_point_3, 
						 unsigned long val_point_4, unsigned long val_point_5, 
						 unsigned long val_point_6, unsigned long val_point_7) : CPrimalGrid() {
	unsigned short iDim, iFace, iNeighbor_Elements;

	/*--- Allocate center-of-gravity coordinates ---*/
	nDim = 3;
	Coord_CG = new double[nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Coord_CG[iDim] = 0.0;
	Coord_FaceElems_CG = new double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_iPoint;	Nodes[1] = val_jPoint;
	Nodes[2] = val_point_2;	Nodes[3] = val_point_3;
	Nodes[4] = val_point_4;	Nodes[5] = val_point_5;
	Nodes[6] = val_point_6;	Nodes[7] = val_point_7;
	
	/*--- Allocate and define neighbor elements to a element ---*/
	nNeighbor_Elements = nFaces;
	Neighbor_Elements = new long[nNeighbor_Elements];
	for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
		Neighbor_Elements[iNeighbor_Elements]=-1;
	}
  
}

CHexahedron::~CHexahedron() {
  unsigned short iFaces;
  
  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;
  
}

void CHexahedron::Change_Orientation(void) {
	unsigned long Point_0, Point_1, Point_2, Point_3, Point_4, Point_5, Point_6, Point_7;
	Point_0 = Nodes[0];
	Point_1 = Nodes[1];
	Point_2 = Nodes[2];
	Point_3 = Nodes[3];
	Point_4 = Nodes[4];
	Point_5 = Nodes[5];
	Point_6 = Nodes[6];
	Point_7 = Nodes[7];
	
	Nodes[0] = Point_7;
	Nodes[1] = Point_4;
	Nodes[2] = Point_5;
	Nodes[3] = Point_6;
	Nodes[4] = Point_3;
	Nodes[5] = Point_0;
	Nodes[6] = Point_1;
	Nodes[7] = Point_2;
  
}

unsigned short CWedge::Faces[5][4] = {{3,4,1,0},{5,2,1,4},{2,5,3,0},{0,1,2,2},{5,4,3,3}};

unsigned short CWedge::Neighbor_Nodes[6][3] = {{1,2,3},{0,2,4},{1,0,5},{0,4,5},{3,5,1},{4,3,2}};

unsigned short CWedge::nNodesFace[5] = {4,4,4,3,3};

unsigned short CWedge::nNeighbor_Nodes[6] = {3,3,3,3,3,3};

unsigned short CWedge::nFaces = 5;

unsigned short CWedge::nNodes = 6;

unsigned short CWedge::nNeighbor_Elements = 5;

unsigned short CWedge::VTK_Type = 13;

unsigned short CWedge::maxNodesFace = 4;

CWedge::CWedge(unsigned long val_point_0, unsigned long val_point_1, 
			   unsigned long val_point_2, unsigned long val_point_3, 
			   unsigned long val_point_4, unsigned long val_point_5) : CPrimalGrid() {
	unsigned short iDim, iFace, iNeighbor_Elements;

	/*--- Allocate CG coordinates ---*/
	nDim = 3;
	Coord_CG = new double[nDim];
	for (iDim = 0; iDim < nDim; iDim++) Coord_CG[iDim] = 0.0;
	
	Coord_FaceElems_CG = new double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_point_0;
	Nodes[1] = val_point_1;
	Nodes[2] = val_point_2;
	Nodes[3] = val_point_3;
	Nodes[4] = val_point_4;
	Nodes[5] = val_point_5;
	
	/*--- Allocate and define neighbor elements to a element ---*/
	nNeighbor_Elements = nFaces;
	Neighbor_Elements = new long[nNeighbor_Elements];
	for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
		Neighbor_Elements[iNeighbor_Elements]=-1;
	}
  
}

CWedge::~CWedge() {
  unsigned short iFaces;
  
  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;
  
}

void CWedge::Change_Orientation(void) {
	unsigned long Point_0, Point_1, Point_3, Point_4;
	Point_0 = Nodes[0];
	Point_1 = Nodes[1];
	Point_3 = Nodes[3];
	Point_4 = Nodes[4];
	
	Nodes[0] = Point_1;
	Nodes[1] = Point_0;
	Nodes[3] = Point_4;
	Nodes[4] = Point_3;
}

unsigned short CPyramid::Faces[5][4] = {{0,3,2,1},{4,3,0,0},{4,0,1,1},{2,4,1,1},{3,4,2,2}};

unsigned short CPyramid::Neighbor_Nodes[5][4] = {{1,3,4,4},{0,2,4,4},{1,3,4,4},{2,0,4,4},{0,1,2,3}};

unsigned short CPyramid::nNodesFace[5] = {4,3,3,3,3};

unsigned short CPyramid::nNeighbor_Nodes[5] = {3,3,3,3,4};

unsigned short CPyramid::nFaces = 5;

unsigned short CPyramid::nNodes = 5;

unsigned short CPyramid::nNeighbor_Elements = 5;

unsigned short CPyramid::VTK_Type = 14;

unsigned short CPyramid::maxNodesFace = 4;

CPyramid::CPyramid(unsigned long val_iPoint, unsigned long val_jPoint, 
				   unsigned long val_point_2, unsigned long val_point_3, 
				   unsigned long val_point_4) : CPrimalGrid() {
	unsigned short iDim, iFace, iNeighbor_Elements;

	/*--- Allocate CG coordinates ---*/
	nDim = 3;
	Coord_CG = new double[nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Coord_CG[iDim] = 0.0;
	Coord_FaceElems_CG = new double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_iPoint;
	Nodes[1] = val_jPoint;
	Nodes[2] = val_point_2;
	Nodes[3] = val_point_3;
	Nodes[4] = val_point_4;
	
	/*--- Allocate and define neighbor elements to a element ---*/
	nNeighbor_Elements = nFaces;
	Neighbor_Elements = new long[nNeighbor_Elements];
	for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
		Neighbor_Elements[iNeighbor_Elements]=-1;
	}
  
}

CPyramid::~CPyramid() {
  unsigned short iFaces;
  
  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;
  
}

void CPyramid::Change_Orientation(void) { cout << "Not defined orientation change" << endl; }
