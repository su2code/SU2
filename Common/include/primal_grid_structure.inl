/*!
 * \file primal_grid_structure.inl
 * \brief In-Line subroutines of the <i>primal_grid_structure.hpp</i> file.
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#pragma once

inline unsigned short CPrimalGrid::GetnNodesFace(unsigned short val_face) { return 0; }

inline void CPrimalGrid::SetDomainElement(unsigned long val_domainelement) { DomainElement = val_domainelement; }

inline void CPrimalGrid::SetRotation_Type(unsigned short val_rotation_type) { }

inline unsigned short CPrimalGrid::GetRotation_Type(void) { return 0; }

inline unsigned short CPrimalGrid::GetMatching_Zone(void) { return 0; }

inline void CPrimalGrid::SetMatching_Zone(unsigned short val_matching_zone) { }

inline unsigned long CPrimalGrid::GetDomainElement(void) { return DomainElement; }

inline void CPrimalGrid::SetNeighbor_Elements(unsigned long val_elem, unsigned short val_face) { Neighbor_Elements[val_face] = val_elem; }

inline long CPrimalGrid::GetNeighbor_Elements(unsigned short val_face) { return Neighbor_Elements[val_face]; }

inline double CPrimalGrid::GetCG(unsigned short val_dim) { return Coord_CG[val_dim]; }

inline double CPrimalGrid::GetFaceCG(unsigned short val_face, unsigned short val_dim) { return Coord_FaceElems_CG[val_face][val_dim]; }

inline void CPrimalGrid::SetDivide (bool val_divide) {	Divide = val_divide; }

inline bool CPrimalGrid::GetDivide (void) { return Divide; }

inline unsigned short CVertexMPI::GetnNodes(void) { return nNodes; }

inline unsigned long CVertexMPI::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline unsigned short CVertexMPI::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CVertexMPI::GetRotation_Type(void) { return Rotation_Type; }

inline void CVertexMPI::SetRotation_Type(unsigned short val_rotation_type) { Rotation_Type = val_rotation_type; }

inline unsigned short CVertexMPI::GetMatching_Zone(void) { return Matching_Zone; }

inline void CVertexMPI::SetMatching_Zone(unsigned short val_matching_zone) { Matching_Zone = val_matching_zone; }

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

inline unsigned short CTriangle::GetnNodes(void) { return nNodes; }

inline unsigned short CTriangle::GetnFaces(void) { return nFaces; }

inline unsigned short CTriangle::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CTriangle::GetMaxNodesFace(void) { return maxNodesFace; }

inline unsigned short CTriangle::GetnNeighbor_Elements(void) { return nNeighbor_Elements; }

inline void CTriangle::SetDomainElement(unsigned long val_domainelement) { DomainElement = val_domainelement; }

inline unsigned long CTriangle::GetDomainElement(void) { return DomainElement; }

inline unsigned short CRectangle::GetFaces(unsigned short val_face, unsigned short val_index) { return Faces[val_face][val_index]; }

inline unsigned short CRectangle::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return Neighbor_Nodes[val_node][val_index]; }

inline unsigned short CRectangle::GetnNodesFace(unsigned short val_face) { return nNodesFace[val_face]; }

inline unsigned short CRectangle::GetnNeighbor_Nodes(unsigned short val_node) { return nNeighbor_Nodes[val_node]; }

inline unsigned long CRectangle::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline unsigned short CRectangle::GetnNodes(void) { return nNodes; }

inline unsigned short CRectangle::GetnFaces(void) { return nFaces; }

inline unsigned short CRectangle::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CRectangle::GetMaxNodesFace(void) { return maxNodesFace; }

inline unsigned short CRectangle::GetnNeighbor_Elements(void) { return nNeighbor_Elements; }

inline void CRectangle::SetDomainElement(unsigned long val_domainelement) {	DomainElement = val_domainelement; }

inline unsigned long CRectangle::GetDomainElement(void) { return DomainElement; }

inline unsigned short CTetrahedron::GetFaces(unsigned short val_face, unsigned short val_index) { return Faces[val_face][val_index]; }

inline unsigned short CTetrahedron::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return Neighbor_Nodes[val_node][val_index]; }

inline unsigned short CTetrahedron::GetnNodesFace(unsigned short val_face) { return nNodesFace[val_face]; }

inline unsigned short CTetrahedron::GetnNeighbor_Nodes(unsigned short val_node) { return nNeighbor_Nodes[val_node]; }

inline unsigned long CTetrahedron::GetNode(unsigned short val_node) { return Nodes[val_node]; }

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

inline unsigned short CHexahedron::GetnNodes(void) { return nNodes; }

inline unsigned short CHexahedron::GetnFaces(void) { return nFaces; }

inline unsigned short CHexahedron::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CHexahedron::GetMaxNodesFace(void) { return maxNodesFace; }

inline unsigned short CHexahedron::GetnNeighbor_Elements(void) { return nNeighbor_Elements; }

inline unsigned short CWedge::GetFaces(unsigned short val_face, unsigned short val_index) { return Faces[val_face][val_index]; }

inline unsigned short CWedge::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return Neighbor_Nodes[val_node][val_index]; }

inline unsigned short CWedge::GetnNodesFace(unsigned short val_face) { return nNodesFace[val_face]; }

inline unsigned short CWedge::GetnNeighbor_Nodes(unsigned short val_node) { return nNeighbor_Nodes[val_node]; }

inline unsigned long CWedge::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline unsigned short CWedge::GetnNodes(void) { return nNodes; }

inline unsigned short CWedge::GetnFaces(void) { return nFaces; }

inline unsigned short CWedge::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CWedge::GetMaxNodesFace(void) { return maxNodesFace; }

inline unsigned short CWedge::GetnNeighbor_Elements(void) { return nNeighbor_Elements; }

inline unsigned short CPyramid::GetFaces(unsigned short val_face, unsigned short val_index) { return Faces[val_face][val_index]; }

inline unsigned short CPyramid::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return Neighbor_Nodes[val_node][val_index]; }

inline unsigned short CPyramid::GetnNodesFace(unsigned short val_face) { return nNodesFace[val_face]; }

inline unsigned short CPyramid::GetnNeighbor_Nodes(unsigned short val_node) { return nNeighbor_Nodes[val_node]; }

inline unsigned long CPyramid::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline unsigned short CPyramid::GetnNodes(void) { return nNodes; }

inline unsigned short CPyramid::GetnFaces(void) { return nFaces; }

inline unsigned short CPyramid::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CPyramid::GetMaxNodesFace(void) { return maxNodesFace; }

inline unsigned short CPyramid::GetnNeighbor_Elements(void) { return nNeighbor_Elements; }
