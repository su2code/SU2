/*!
 * \file primal_grid_structure.inl
 * \brief In-Line subroutines of the <i>primal_grid_structure.hpp</i> file.
 * \author F. Palacios
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

inline void CPrimalGrid::SetColor(unsigned long val_color) { }

inline unsigned long CPrimalGrid::GetColor(void) { return -1; }

inline unsigned short CPrimalGrid::GetRotation_Type(void) { return 0; }

inline unsigned long CPrimalGrid::GetDomainElement(void) { return DomainElement; }

inline void CPrimalGrid::SetNeighbor_Elements(unsigned long val_elem, unsigned short val_face) { Neighbor_Elements[val_face] = val_elem; }

inline su2double CPrimalGrid::GetLengthScale(void) { return LenScale; }

inline void CPrimalGrid::SetLengthScale(su2double val_lenScale) { LenScale = val_lenScale; }

inline unsigned short CPrimalGrid::GetTimeLevel(void) { return TimeLevel; }

inline void CPrimalGrid::SetTimeLevel(unsigned short val_timeLevel) { TimeLevel = val_timeLevel; }

inline bool CPrimalGrid::GetOwnerFace(unsigned short val_face) { return ElementOwnsFace[val_face]; }

inline void CPrimalGrid::SetOwnerFace(bool val_owner, unsigned short val_face) { ElementOwnsFace[val_face] = val_owner; }

inline short CPrimalGrid::GetPeriodicIndex(unsigned short val_face) {return PeriodIndexNeighbors[val_face];}

inline void CPrimalGrid::SetPeriodicIndex(unsigned short val_periodic, unsigned short val_face) {PeriodIndexNeighbors[val_face] = val_periodic; }

inline bool CPrimalGrid::GetJacobianConstantFace(unsigned short val_face) { return JacobianFaceIsConstant[val_face]; }

inline void CPrimalGrid::SetJacobianConstantFace(bool val_JacFaceIsConstant, unsigned short val_face) {JacobianFaceIsConstant[val_face] = val_JacFaceIsConstant; }

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

inline void CPrimalGrid::GetCornerPointsAllFaces(unsigned short &nFaces,
                                                 unsigned short nPointsPerFace[],
                                                 unsigned long  faceConn[6][4]) { }

inline unsigned long CPrimalGrid::GetGlobalElemID(void) { return 0; }

inline unsigned long CPrimalGrid::GetGlobalOffsetDOFsSol(void) { return 0; }

inline unsigned short CPrimalGrid::GetNPolyGrid(void) { return 0; }

inline unsigned short CPrimalGrid::GetNPolySol(void) { return 0; }

inline unsigned short CPrimalGrid::GetNDOFsGrid(void) { return 0; }

inline unsigned short CPrimalGrid::GetNDOFsSol(void) { return 0; }

inline bool CPrimalGrid::GetJacobianConsideredConstant(void) { return false; }

inline void CPrimalGrid::SetJacobianConsideredConstant(bool val_JacobianConsideredConstant) {}

inline void CPrimalGrid::AddOffsetGlobalDOFs(const unsigned long val_offsetRank) {}

inline void CPrimalGrid::AddDonorWallFunctions(const unsigned long donorElement) {}

inline unsigned short CPrimalGrid::GetNDonorsWallFunctions(void) {return 0;}

inline unsigned long *CPrimalGrid::GetDonorsWallFunctions(void) {return NULL;}

inline void CPrimalGrid::SetDonorsWallFunctions(const vector<unsigned long> &donorElements) {}

inline void CPrimalGrid::RemoveMultipleDonorsWallFunctions(void) {}

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

inline unsigned long CPrimalGridFEM::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline unsigned short CPrimalGridFEM::GetnNodesFace(unsigned short val_face) { return -1; }

inline unsigned short CPrimalGridFEM::GetFaces(unsigned short val_face, unsigned short val_index) { return -1; }

inline unsigned short CPrimalGridFEM::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return -1; }

inline unsigned short CPrimalGridFEM::GetnNodes(void) { return nDOFsGrid; }

inline unsigned short CPrimalGridFEM::GetnFaces(void) { return nFaces; }

inline unsigned short CPrimalGridFEM::GetnNeighbor_Nodes(unsigned short val_node) { return -1; }

inline void CPrimalGridFEM::Change_Orientation(void) {}

inline unsigned long CPrimalGridFEM::GetGlobalElemID(void) { return elemIDGlobal; }

inline unsigned long CPrimalGridFEM::GetGlobalOffsetDOFsSol(void) { return offsetDOFsSolGlobal; }

inline unsigned short CPrimalGridFEM::GetnNeighbor_Elements(void) { return nFaces; }

inline unsigned short CPrimalGridFEM::GetMaxNodesFace(void) { return -1; }

inline unsigned short CPrimalGridFEM::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CPrimalGridFEM::GetNPolyGrid(void) { return nPolyGrid; }

inline unsigned short CPrimalGridFEM::GetNPolySol(void) { return nPolySol; }

inline unsigned short CPrimalGridFEM::GetNDOFsGrid(void) { return nDOFsGrid; }

inline unsigned short CPrimalGridFEM::GetNDOFsSol(void) { return nDOFsSol; }

inline bool CPrimalGridFEM::GetJacobianConsideredConstant(void) { return JacobianConsideredConstant; }

inline void CPrimalGridFEM::SetColor(unsigned long val_color) { color = val_color; }

inline unsigned long CPrimalGridFEM::GetColor(void) { return color; }

inline void CPrimalGridFEM::SetJacobianConsideredConstant(bool val_JacobianConsideredConstant) {JacobianConsideredConstant = val_JacobianConsideredConstant;}

inline void CPrimalGridFEM::AddOffsetGlobalDOFs(const unsigned long val_offsetRank) {offsetDOFsSolGlobal += val_offsetRank;}

inline CPrimalGridBoundFEM::~CPrimalGridBoundFEM(){}

inline unsigned long CPrimalGridBoundFEM::GetNode(unsigned short val_node) { return Nodes[val_node]; }

inline unsigned short CPrimalGridBoundFEM::GetnNodesFace(unsigned short val_face) { return -1; }

inline unsigned short CPrimalGridBoundFEM::GetFaces(unsigned short val_face, unsigned short val_index) { return -1; }

inline unsigned short CPrimalGridBoundFEM::GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) { return -1; }

inline unsigned short CPrimalGridBoundFEM::GetnNodes(void) { return nDOFsGrid; }

inline unsigned short CPrimalGridBoundFEM::GetnFaces(void) { return -1; }

inline unsigned short CPrimalGridBoundFEM::GetnNeighbor_Nodes(unsigned short val_node) { return -1; }

inline void CPrimalGridBoundFEM::Change_Orientation(void) {}

inline unsigned long CPrimalGridBoundFEM::GetGlobalElemID(void) { return boundElemIDGlobal; }

inline unsigned short CPrimalGridBoundFEM::GetnNeighbor_Elements(void) { return -1; }

inline unsigned short CPrimalGridBoundFEM::GetMaxNodesFace(void) { return -1; }

inline unsigned short CPrimalGridBoundFEM::GetVTK_Type(void) { return VTK_Type; }

inline unsigned short CPrimalGridBoundFEM::GetNPolyGrid(void) { return nPolyGrid; }

inline unsigned short CPrimalGridBoundFEM::GetNDOFsGrid(void) { return nDOFsGrid; }

inline bool CPrimalGridBoundFEM::GetJacobianConsideredConstant(void) {return JacobianConsideredConstant;}

inline void CPrimalGridBoundFEM::SetJacobianConsideredConstant(bool val_JacobianConsideredConstant) {JacobianConsideredConstant = val_JacobianConsideredConstant;}

inline void CPrimalGridBoundFEM::AddDonorWallFunctions(const unsigned long donorElement) {donorElementsWallFunctions.push_back(donorElement);}

inline unsigned short CPrimalGridBoundFEM::GetNDonorsWallFunctions(void) {return donorElementsWallFunctions.size();}

inline unsigned long *CPrimalGridBoundFEM::GetDonorsWallFunctions(void) {return donorElementsWallFunctions.data();}

inline void CPrimalGridBoundFEM::SetDonorsWallFunctions(const vector<unsigned long> &donorElements) {donorElementsWallFunctions = donorElements;}
