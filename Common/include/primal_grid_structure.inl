/*!
 * \file primal_grid_structure.inl
 * \brief In-Line subroutines of the <i>primal_grid_structure.hpp</i> file.
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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
 
#pragma once

inline unsigned short CPrimalGrid::GetnNodesFace(unsigned short val_face) { return 0; }

inline void CPrimalGrid::SetDomainElement(unsigned long val_domainelement) { DomainElement = val_domainelement; }

inline void CPrimalGrid::SetRotation_Type(unsigned short val_rotation_type) { }

inline unsigned short CPrimalGrid::GetRotation_Type(void) { return 0; }

inline unsigned long CPrimalGrid::GetDomainElement(void) { return DomainElement; }

inline void CPrimalGrid::SetNeighbor_Elements(unsigned long val_elem, unsigned short val_face) { Neighbor_Elements[val_face] = val_elem; }

inline long CPrimalGrid::GetNeighbor_Elements(unsigned short val_face) { return Neighbor_Elements[val_face]; }

inline su2double CPrimalGrid::GetVolume(void) { return Volume; }

inline void CPrimalGrid::SetVolume(su2double val_volume) { Volume = val_volume; }

inline su2double CPrimalGrid::GetCG(unsigned short val_dim) { return Coord_CG[val_dim]; }

inline su2double CPrimalGrid::GetFaceCG(unsigned short val_face, unsigned short val_dim) { return Coord_FaceElems_CG[val_face][val_dim]; }

inline void CPrimalGrid::SetDivide (bool val_divide) {	Divide = val_divide; }

inline bool CPrimalGrid::GetDivide (void) { return Divide; }

inline unsigned long CPrimalGrid::GetGlobalIndex(void) { return GlobalIndex; }

inline void CPrimalGrid::SetGlobalIndex(unsigned long val_globalindex) { GlobalIndex = val_globalindex; }

inline void CPrimalGrid::SetNode(unsigned short val_node, unsigned long val_point) { }

inline unsigned short CVertexMPI::GetnNodes(void) { return nNodes; }

inline unsigned long CVertexMPI::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline void CVertexMPI::SetNode(unsigned short val_node, unsigned long val_point) { Nodes[val_node] = val_point; }

inline unsigned short CVertexMPI::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CVertexMPI::GetRotation_Type(void) { return Rotation_Type; }

inline void CVertexMPI::SetRotation_Type(unsigned short val_rotation_type) { Rotation_Type = val_rotation_type; }

inline unsigned short CVertexMPI::GetnNeighbor_Nodes(unsigned short val_node) { return 0; }

inline unsigned short CVertexMPI::GetnNeighbor_Elements(void) { return 0; }

inline unsigned short CVertexMPI::GetnFaces(void) { return 0; }

inline unsigned short CVertexMPI::GetnNodesFace(unsigned short val_face) { return 0; }

inline unsigned short CVertexMPI::GetMaxNodesFace(void) { return 0; }

inline unsigned short CVertexMPI::GetFaces(unsigned short val_face, unsigned short val_index) { return 0; }

inline unsigned short CVertexMPI::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return 0; }

inline unsigned short CLine::GetFaces(unsigned short val_face, unsigned short val_index) { return Faces[val_face][val_index]; }

inline unsigned short CLine::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return Neighbor_Nodes[val_node][val_index]; }

inline unsigned short CLine::GetnNodesFace(unsigned short val_face) { return nNodesFace[val_face]; }

inline unsigned short CLine::GetnNeighbor_Nodes(unsigned short val_node) { return nNeighbor_Nodes[val_node]; }

inline unsigned long CLine::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline void CLine::SetNode(unsigned short val_node, unsigned long val_point) { Nodes[val_node] = val_point; }

inline unsigned short CLine::GetnNodes(void) { return nNodes; }

inline unsigned short CLine::GetnFaces(void) { return nFaces; }

inline unsigned short CLine::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CLine::GetMaxNodesFace(void) { return maxNodesFace; }

inline unsigned short CLine::GetnNeighbor_Elements(void) { return nNeighbor_Elements; }

inline void CLine::SetDomainElement(unsigned long val_domainelement) {DomainElement = val_domainelement; }

inline unsigned long CLine::GetDomainElement(void) { return DomainElement; }

inline unsigned short CTriangle::GetFaces(unsigned short val_face, unsigned short val_index) { return Faces[val_face][val_index]; }

inline unsigned short CTriangle::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return Neighbor_Nodes[val_node][val_index]; }

inline unsigned short CTriangle::GetnNodesFace(unsigned short val_face) { return nNodesFace[val_face]; }

inline unsigned short CTriangle::GetnNeighbor_Nodes(unsigned short val_node) { return nNeighbor_Nodes[val_node]; }

inline unsigned long CTriangle::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline void CTriangle::SetNode(unsigned short val_node, unsigned long val_point) { Nodes[val_node] = val_point; }

inline unsigned short CTriangle::GetnNodes(void) { return nNodes; }

inline unsigned short CTriangle::GetnFaces(void) { return nFaces; }

inline unsigned short CTriangle::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CTriangle::GetMaxNodesFace(void) { return maxNodesFace; }

inline unsigned short CTriangle::GetnNeighbor_Elements(void) { return nNeighbor_Elements; }

inline void CTriangle::SetDomainElement(unsigned long val_domainelement) { DomainElement = val_domainelement; }

inline unsigned long CTriangle::GetDomainElement(void) { return DomainElement; }
inline unsigned short CQuadrilateral::GetFaces(unsigned short val_face, unsigned short val_index) { return Faces[val_face][val_index]; }

inline unsigned short CQuadrilateral::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return Neighbor_Nodes[val_node][val_index]; }

inline unsigned short CQuadrilateral::GetnNodesFace(unsigned short val_face) { return nNodesFace[val_face]; }

inline unsigned short CQuadrilateral::GetnNeighbor_Nodes(unsigned short val_node) { return nNeighbor_Nodes[val_node]; }

inline unsigned long CQuadrilateral::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline void CQuadrilateral::SetNode(unsigned short val_node, unsigned long val_point) { Nodes[val_node] = val_point; }

inline unsigned short CQuadrilateral::GetnNodes(void) { return nNodes; }

inline unsigned short CQuadrilateral::GetnFaces(void) { return nFaces; }

inline unsigned short CQuadrilateral::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CQuadrilateral::GetMaxNodesFace(void) { return maxNodesFace; }

inline unsigned short CQuadrilateral::GetnNeighbor_Elements(void) { return nNeighbor_Elements; }

inline void CQuadrilateral::SetDomainElement(unsigned long val_domainelement) {	DomainElement = val_domainelement; }

inline unsigned long CQuadrilateral::GetDomainElement(void) { return DomainElement; }

inline unsigned short CTetrahedron::GetFaces(unsigned short val_face, unsigned short val_index) { return Faces[val_face][val_index]; }

inline unsigned short CTetrahedron::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return Neighbor_Nodes[val_node][val_index]; }

inline unsigned short CTetrahedron::GetnNodesFace(unsigned short val_face) { return nNodesFace[val_face]; }

inline unsigned short CTetrahedron::GetnNeighbor_Nodes(unsigned short val_node) { return nNeighbor_Nodes[val_node]; }

inline unsigned long CTetrahedron::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline void CTetrahedron::SetNode(unsigned short val_node, unsigned long val_point) { Nodes[val_node] = val_point; }

inline unsigned short CTetrahedron::GetnNodes(void) { return nNodes; }

inline unsigned short CTetrahedron::GetnFaces(void) { return nFaces; }

inline unsigned short CTetrahedron::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CTetrahedron::GetMaxNodesFace(void) { return maxNodesFace; }

inline unsigned short CTetrahedron::GetnNeighbor_Elements(void) { return nNeighbor_Elements; }

inline unsigned short CHexahedron::GetFaces(unsigned short val_face, unsigned short val_index) { return Faces[val_face][val_index]; }

inline unsigned short CHexahedron::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return Neighbor_Nodes[val_node][val_index]; }

inline unsigned short CHexahedron::GetnNodesFace(unsigned short val_face) { return nNodesFace[val_face]; }

inline unsigned short CHexahedron::GetnNeighbor_Nodes(unsigned short val_node) { return nNeighbor_Nodes[val_node]; }

inline unsigned long CHexahedron::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline void CHexahedron::SetNode(unsigned short val_node, unsigned long val_point) { Nodes[val_node] = val_point; }

inline unsigned short CHexahedron::GetnNodes(void) { return nNodes; }

inline unsigned short CHexahedron::GetnFaces(void) { return nFaces; }

inline unsigned short CHexahedron::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CHexahedron::GetMaxNodesFace(void) { return maxNodesFace; }

inline unsigned short CHexahedron::GetnNeighbor_Elements(void) { return nNeighbor_Elements; }

inline unsigned short CPrism::GetFaces(unsigned short val_face, unsigned short val_index) { return Faces[val_face][val_index]; }

inline unsigned short CPrism::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return Neighbor_Nodes[val_node][val_index]; }

inline unsigned short CPrism::GetnNodesFace(unsigned short val_face) { return nNodesFace[val_face]; }

inline unsigned short CPrism::GetnNeighbor_Nodes(unsigned short val_node) { return nNeighbor_Nodes[val_node]; }

inline unsigned long CPrism::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline void CPrism::SetNode(unsigned short val_node, unsigned long val_point) { Nodes[val_node] = val_point; }

inline unsigned short CPrism::GetnNodes(void) { return nNodes; }

inline unsigned short CPrism::GetnFaces(void) { return nFaces; }

inline unsigned short CPrism::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CPrism::GetMaxNodesFace(void) { return maxNodesFace; }

inline unsigned short CPrism::GetnNeighbor_Elements(void) { return nNeighbor_Elements; }

inline unsigned short CPyramid::GetFaces(unsigned short val_face, unsigned short val_index) { return Faces[val_face][val_index]; }

inline unsigned short CPyramid::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return Neighbor_Nodes[val_node][val_index]; }

inline unsigned short CPyramid::GetnNodesFace(unsigned short val_face) { return nNodesFace[val_face]; }

inline unsigned short CPyramid::GetnNeighbor_Nodes(unsigned short val_node) { return nNeighbor_Nodes[val_node]; }

inline unsigned long CPyramid::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline void CPyramid::SetNode(unsigned short val_node, unsigned long val_point) { Nodes[val_node] = val_point; }

inline unsigned short CPyramid::GetnNodes(void) { return nNodes; }

inline unsigned short CPyramid::GetnFaces(void) { return nFaces; }

inline unsigned short CPyramid::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CPyramid::GetMaxNodesFace(void) { return maxNodesFace; }

inline unsigned short CPyramid::GetnNeighbor_Elements(void) { return nNeighbor_Elements; }
