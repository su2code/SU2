/*!
 * \file primal_grid_structure.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>primal_grid_structure.cpp</i> file.
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

#include <iostream>
#include <vector>

#include "dual_grid_structure.hpp"
#include "config_structure.hpp"

using namespace std;

/*!
 * \class CPrimalGrid
 * \brief Class to define the numerical primal grid.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CPrimalGrid {
protected:
	unsigned long *Nodes;         /*!< \brief Vector to store the global nodes of an element. */
	long *Neighbor_Elements;      /*!< \brief Vector to store the elements surronding an element. */
	double *Coord_CG;             /*!< \brief Coordinates of the center-of-gravity of the element. */
	double **Coord_FaceElems_CG;	/*!< \brief Coordinates of the center-of-gravity of the face of the
                                 elements. */
	static unsigned short nDim;		/*!< \brief Dimension of the element (2D or 3D) useful for triangles,
                                 rectangles and edges. */
	unsigned long DomainElement;	/*!< \brief Only for boundaries, in this variable the 3D elements which
                                 correspond with a boundary element is stored. */
	bool Divide;                  /*!< \brief Marker used to know if we are going to divide this element
                                 in the adaptation proccess. */
  
public:
	
	/*!
	 * \brief Constructor of the class.
	 */
	CPrimalGrid(void);
	
	/*!
	 * \overload
	 * \param[in] val_nNodes - Number of nodes of the element.
	 * \param[in] val_nFaces - Number of faces of the element.
	 * \param[in] val_VTK_Type - Type of the element using the vtk nomenclature.
	 */
	CPrimalGrid(unsigned short val_nNodes, unsigned short val_nFaces, unsigned short val_VTK_Type);
	
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CPrimalGrid(void);
	
	/*!
	 * \brief Get the elements that surround an element.
	 * \param[in] val_face - Local index of the face.
	 * \return Global index of the element.
	 */
	long GetNeighbor_Elements(unsigned short val_face);
	
	/*!
	 * \brief Set the elements that surround an element.
	 * \param[in] val_elem - Global index of the element.
	 * \param[in] val_face - Local index of the face.
	 */
	void SetNeighbor_Elements(unsigned long val_elem, unsigned short val_face);
	
	/*!
	 * \brief Set the center of gravity of an element (including edges).
	 * \param[in] val_coord - Coordinates of the element.
	 */
	void SetCG(double **val_coord);
	
	/*!
	 * \brief Get the center of gravity of an element (including edges).
	 * \param[in] val_dim - Coordinate of the center of gravity.
	 * \return Coordinates of the center of gravity.
	 */
	double GetCG(unsigned short val_dim);
	
	/*!
	 * \brief Get the CG of a face of an element.
	 * \param[in] val_face - Local index of the face.
	 * \param[in] val_dim - Coordinate of the center of gravity.
	 * \return Coordinates of the center of gravity.
	 */
	double GetFaceCG(unsigned short val_face, unsigned short val_dim);
	
	/*!
	 * \brief Get all the neighbors of an element.
	 * \return List of all the neighbor of an element.
	 */
	void GetAllNeighbor_Elements(void);
	
	/*!
	 * \brief Set that an element must be divided in the adaptation stage.
	 * \param[in] val_divide - <code>TRUE</code> if the element must be divided; otherwise <code>FALSE</code>.
	 */
	void SetDivide(bool val_divide);
  
	/*!
	 * \brief Get if an element must be divided in the adaptation stage.
	 * \return <code>TRUE</code> if the element must be divided; otherwise <code>FALSE</code>.
	 */
	bool GetDivide(void);
	
	/*!
	 * \brief A virtual member.
	 * \param[in] val_domainelement Index of the domain element which has a face shared by this boundary element.
	 */
	virtual void SetDomainElement(unsigned long val_domainelement);
	
	/*!
	 * \brief A virtual member.
	 * \return Relate the boundary element which a face of a domain element.
	 */
	virtual unsigned long GetDomainElement(void);
  
	/*!
	 * \brief A pure virtual member.
	 */
	virtual void Change_Orientation(void) = 0;
	
	/*!
	 * \brief A pure virtual member.
	 * \return Kind of element using the vtk nomenclature.
	 */
	virtual unsigned short GetVTK_Type(void) = 0;
	
	/*!
	 * \brief A pure virtual member.
	 * \return Type of the element using VTK nomenclature.
	 */
	virtual unsigned short GetRotation_Type(void);
	
	/*!
	 * \brief A pure virtual member.
	 * \param[in] val_rotation_type - Kind of rotation/traslation that must be applied.
	 */
	virtual void SetRotation_Type(unsigned short val_rotation_type);
	
  /*!
	 * \brief A pure virtual member.
	 * \return Matching zone for a sliding interface vertex.
	 */
	virtual unsigned short GetMatching_Zone(void);
	
	/*!
	 * \brief A pure virtual member.
	 * \param[in] val_matching_zone - Matching zone index for a sliding interface vertex.
	 */
	virtual void SetMatching_Zone(unsigned short val_matching_zone);
  
	/*!
	 * \brief A pure virtual member.
	 * \param[in] val_node - Local index of the node.
	 * \return Number of neighbors nodes of a node in the element.
	 */
	virtual unsigned short GetnNeighbor_Nodes(unsigned short val_node) = 0;
	
	/*!
	 * \brief A pure virtual member.
	 * \return Number of neighbors elements of a element.
	 */
	virtual unsigned short GetnNeighbor_Elements(void) = 0;
	
	/*!
	 * \brief A pure virtual member.
	 * \return Number of nodes of an element.
	 */
	virtual unsigned short GetnNodes(void) = 0;
	
	/*!
	 * \brief A pure virtual member.
	 * \return Number of faces of an element.
	 */
	virtual unsigned short GetnFaces(void) = 0;
	
	/*!
	 * \brief A pure virtual member.
	 * \param[in] val_face - Local index of a face.
	 * \return Local index of the nodes that compose a face.
	 */
	virtual unsigned short GetnNodesFace(unsigned short val_face) = 0;
	
	/*!
	 * \brief A pure virtual member.
	 * \return Maximum number of nodes that compose a face.
	 */
	virtual unsigned short GetMaxNodesFace(void) = 0;
	
	/*!
	 * \brief A pure virtual member.
	 * \param[in] val_node - Local index of a node.
	 * \return Global index of the node.
	 */
	virtual unsigned long GetNode(unsigned short val_node) = 0;
	
	/*!
	 * \brief A pure virtual member.
	 * \param[in] val_face - Local index of the face.
	 * \param[in] val_index - Local index of the nodes that compose the face.
	 * \return - Local index of the nodes that compose the face.
	 */
	virtual unsigned short GetFaces(unsigned short val_face, unsigned short val_index) = 0;
	
	/*!
	 * \brief A pure virtual member.
	 * \param[in] val_node - Local index of a node.
	 * \param[in] val_index - Local (to the neighbor nodes of <i>val_node</i>) index of the nodes that
   are neighbor to val_node.
	 * \return Local index of the nodes that are neighbor to <i>val_node</i>.
	 */
	virtual unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) = 0;
};

/*!
 * \class CVertexMPI
 * \brief Class for vertex element definition. This kind
 *        of element is used in the parallelization stuff.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CVertexMPI : public CPrimalGrid {
private:
	static unsigned short nFaces;				/*!< \brief Number of faces of the element. */
	static unsigned short nNodes;				/*!< \brief Number of nodes of the element. */
	static unsigned short VTK_Type;				/*!< \brief Type of element using VTK nomenclature. */
	unsigned short Rotation_Type;			/*!< \brief Definition of the rotation, traslation of the
																		 solution at the vertex. */
  unsigned short Matching_Zone;			/*!< \brief Matching/Donor zone for sliding mesh interfaces. */
	static unsigned short maxNodesFace;			/*!< \brief Maximum number of nodes for a face. */
	static unsigned short nNeighbor_Elements;	/*!< \brief Number of Neighbor_Elements. */
	
public:
  
	/*!
	 * \brief Constructor using the nodes and index.
	 * \param[in] val_point - Index of the 1st triangle point read from the grid file.
	 * \param[in] val_nDim - Number of dimension of the problem (2D or 3D).
	 */
	CVertexMPI(unsigned long val_point, unsigned short val_nDim);
  
  /*!
	 * \brief Destructor of the class.
	 */
	~CVertexMPI(void);
	
	/*!
	 * \brief Get the nodes shared by the line.
	 * \param[in] val_node - Local (to the line) index of the node (a line has 2 nodes).
	 * \return Global index of the line node.
	 */
	unsigned long GetNode(unsigned short val_node);
	
	/*!
	 * \brief Get the number of nodes of an element.
	 * \return Number of nodes that composes an element.
	 */
	unsigned short GetnNodes(void);
	
	/*!
	 * \brief Get the type of the element using VTK nomenclature.
	 * \return Type of the element using VTK nomenclature.
	 */
	unsigned short GetVTK_Type(void);
	
	/*!
	 * \brief Get the type of rotation/traslation that must be applied.
	 * \return Type of the element using VTK nomenclature.
	 */
	unsigned short GetRotation_Type(void);
	
	/*!
	 * \brief Set the type of rotation/traslation that must be applied.
	 * \param[in] val_rotation_type - Kind of rotation/traslation that must be applied.
	 */
	void SetRotation_Type(unsigned short val_rotation_type);
  
  /*!
	 * \brief Get the matching zone index for a sliding interface vertex.
	 * \return Matching zone for a sliding interface vertex.
	 */
  unsigned short GetMatching_Zone(void);
	
	/*!
	 * \brief Set the matching zone index for a sliding interface vertex.
	 * \param[in] val_matching_zone - Matching zone index for a sliding interface vertex.
	 */
  void SetMatching_Zone(unsigned short val_matching_zone);
	
	/*!
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the
	 *        definition of the function in all the derived classes).
	 */
	void Change_Orientation(void);
	
	/*!
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the
	 *        definition of the function in all the derived classes).
	 */
	unsigned short GetnNeighbor_Elements(void);
	
	/*!
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the
	 *        definition of the function in all the derived classes).
	 */
	unsigned short GetnNeighbor_Nodes(unsigned short val_node);
	
	/*!
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the
	 *        definition of the function in all the derived classes).
	 */
	unsigned short GetnFaces(void);
	
	/*!
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the
	 *        definition of the function in all the derived classes).
	 */
	unsigned short GetnNodesFace(unsigned short val_face);
	
	/*!
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the
	 *        definition of the function in all the derived classes).
	 */
	unsigned short GetMaxNodesFace(void);
	
	/*!
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the
	 *        definition of the function in all the derived classes).
	 */
	unsigned short GetFaces(unsigned short val_face, unsigned short val_index);
	
	/*!
	 * \brief This function does nothing (it comes from a pure virtual function, that implies the
	 *        definition of the function in all the derived classes).
	 */
	unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index);
};

/*!
 * \class CLine
 * \brief Class for line element definition.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CLine : public CPrimalGrid {
private:
	static unsigned short Faces[1][2];			/*!< \brief Matrix to store the local nodes of all the faces. */
	static unsigned short Neighbor_Nodes[2][1]; /*!< \brief Neighbor to a nodes in the element. */
	static unsigned short nNodesFace[1];		/*!< \brief Number of nodes of each face of the element. */
	static unsigned short nNeighbor_Nodes[2];	/*!< \brief Number of Neighbor to a nodes in the element. */
	static unsigned short nFaces;				/*!< \brief Number of faces of the element. */
	static unsigned short nNodes;				/*!< \brief Number of nodes of the element. */
	static unsigned short VTK_Type;				/*!< \brief Type of element using VTK nomenclature. */
	static unsigned short maxNodesFace;			/*!< \brief Maximum number of nodes for a face. */
	static unsigned short nNeighbor_Elements;	/*!< \brief Number of Neighbor_Elements. */
  
public:
  
	/*!
	 * \brief Constructor using the nodes and index.
	 * \param[in] val_iPoint - Index of the 1st triangle point read from the grid file.
	 * \param[in] val_jPoint - Index of the 2nd triangle point read from the grid file.
	 * \param[in] val_nDim - Number of dimension of the problem (2D or 3D).
	 */
	CLine(unsigned long val_iPoint, unsigned long val_jPoint, unsigned short val_nDim);
  
  /*!
	 * \brief Destructor of the class.
	 */
	~CLine(void);
  
	/*!
	 * \brief Get the nodes shared by the line.
	 * \param[in] val_node - Local (to the line) index of the node (a line has 2 nodes).
	 * \return Global index of the line node.
	 */
	unsigned long GetNode(unsigned short val_node);
  
	/*!
	 * \brief Get the face index of and element.
	 * \param[in] val_face - Local index of the face.
	 * \param[in] val_index - Local (to the face) index of the nodes that compose the face.
	 * \return Local (to the line) index of the nodes that compose the face.
	 */
	unsigned short GetFaces(unsigned short val_face, unsigned short val_index);
  
	/*!
	 * \brief Get the local index of the neighbors to a node (given the local index).
	 * \param[in] val_node - Local (to the line) index of a node.
	 * \param[in] val_index - Local (to the neighbor nodes of val_node) index of the nodes
   that are neighbor to val_node (each face is composed by 3 nodes).
	 * \return Local (to the line) index of the nodes that are neighbor to val_node.
	 */
	unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index);
	
	/*!
	 * \brief Get the number of neighbors nodes of a node.
	 * \param[in] val_node - Local (to the line) index of a node.
	 * \return Number if neighbors of a node val_node.
	 */
	unsigned short GetnNeighbor_Nodes(unsigned short val_node);
  
	/*!
	 * \brief Get the number of nodes that composes a face of an element.
	 * \param[in] val_face - Local index of the face.
	 * \return Number of nodes that composes a face of an element.
	 */
	unsigned short GetnNodesFace(unsigned short val_face);
	
	/*!
	 * \brief Get the number of nodes of an element.
	 * \return Number of nodes that composes an element.
	 */
	unsigned short GetnNodes(void);
  
	/*!
	 * \brief Get the number of faces of an element.
	 * \return Number of faces of an element.
	 */
	unsigned short GetnFaces(void);
  
	/*!
	 * \brief Get the Maximum number of nodes of a face of an element.
	 * \return Maximum number of nodes of a face of an element.
	 */
	unsigned short GetMaxNodesFace(void);
  
	/*!
	 * \brief Get the type of the element using VTK nomenclature.
	 * \return Type of the element using VTK nomenclature.
	 */
	unsigned short GetVTK_Type(void);
	
	/*!
	 * \brief Get the number of element that are neighbor to this element.
	 * \return Number of neighbor elements.
	 */
	unsigned short GetnNeighbor_Elements(void);
  
	/*!
	 * \brief Set the domain element which shares a face with the boundary element.
	 * \param[in] val_domainelement - Global index of the element.
	 */
	void SetDomainElement(unsigned long val_domainelement);
	
	/*!
	 * \brief Get the domain element which shares a face with the boundary element.
	 * \return Domain element which shares a face with the boundary element.
	 */
	unsigned long GetDomainElement(void);
	
	/*!
	 * \brief Change the orientation of an element.
	 */
	void Change_Orientation(void);
};

/*!
 * \class CTriangle
 * \brief Class for triangle element definition.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CTriangle : public CPrimalGrid {
private:
	static unsigned short Faces[3][2];			/*!< \brief Matrix to store the local nodes of all the faces. */
	static unsigned short Neighbor_Nodes[3][2];	/*!< \brief Neighbor to a nodes in the element. */
	static unsigned short nNodesFace[3];		/*!< \brief Number of nodes of each face of the element. */
	static unsigned short nNeighbor_Nodes[3];	/*!< \brief Number of Neighbor to a nodes in the element. */
	static unsigned short nFaces;				/*!< \brief Number of faces of the element. */
	static unsigned short nNodes;				/*!< \brief Number of nodes of the element. */
	static unsigned short VTK_Type;				/*!< \brief Type of element using VTK nomenclature. */
	static unsigned short maxNodesFace;			/*!< \brief Maximum number of nodes for a face. */
	static unsigned short nNeighbor_Elements;	/*!< \brief Number of Neighbor_Elements. */
  
public:
	
	/*!
	 * \brief Constructor using the nodes and index.
	 * \param[in] val_iPoint - Index of the 1st triangle point read from the grid file.
	 * \param[in] val_jPoint - Index of the 2nd triangle point read from the grid file.
	 * \param[in] val_point_2 - Index of the 3th triangle point read from the grid file.
	 * \param[in] val_nDim - Number of dimension of the problem (2D or 3D), be careful a triangle could be 2D or 3D.
	 */
	CTriangle(unsigned long val_iPoint, unsigned long val_jPoint,
            unsigned long val_point_2, unsigned short val_nDim);
  
  /*!
	 * \brief Destructor of the class.
	 */
	~CTriangle(void);
	
	/*!
	 * \brief Get the nodes shared by the triangle.
	 * \param[in] val_node - Local (to the triangle) index of the node (a triangle has 3 nodes).
	 * \return Global index of the triangle node.
	 */
	unsigned long GetNode(unsigned short val_node);
	
	/*!
	 * \brief Get the face index of and element.
	 * \param[in] val_face - Local index of the face.
	 * \param[in] val_index - Local (to the face) index of the nodes that compose the face.
	 * \return Local (to the triangle) index of the nodes that compose the face.
	 */
	unsigned short GetFaces(unsigned short val_face, unsigned short val_index);
	
	/*!
	 * \brief Get the local index of the neighbors to a node (given the local index).
	 * \param[in] val_node - Local (to the triangle) index of a node.
	 * \param[in] val_index - Local (to the neighbor nodes of val_node) index of the nodes that
   are neighbor to val_node (each face is composed by 3 nodes).
	 * \return Local (to the triangle) index of the nodes that are neighbor to val_node.
	 */
	unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index);
  
	/*!
	 * \brief Get the number of neighbors nodes of a node.
	 * \param[in] val_node - Local (to the triangle) index of a node.
	 * \return Number if neighbors of a node val_node.
	 */
	unsigned short GetnNeighbor_Nodes(unsigned short val_node);
	
	/*!
	 * \brief Get the number of nodes that composes a face of an element.
	 * \param[in] val_face - Local index of the face.
	 * \return Number of nodes that composes a face of an element.
	 */
	unsigned short GetnNodesFace(unsigned short val_face);
	
	/*!
	 * \brief Get the number of nodes of an element.
	 * \return Number of nodes that composes an element.
	 */
	unsigned short GetnNodes(void);
  
	/*!
	 * \brief Get the number of faces of an element.
	 * \return Number of faces of an element.
	 */
	unsigned short GetnFaces(void);
  
	/*!
	 * \brief Get the Maximum number of nodes of a face of an element.
	 * \return Maximum number of nodes of a face of an element.
	 */
	unsigned short GetMaxNodesFace(void);
  
	/*!
	 * \brief Get the type of the element using VTK nomenclature.
	 * \return Type of the element using VTK nomenclature.
	 */
	unsigned short GetVTK_Type(void);
  
	/*!
	 * \brief Get the number of element that are neighbor to this element.
	 * \return Number of neighbor elements.
	 */
	unsigned short GetnNeighbor_Elements(void);
  
	/*!
	 * \brief Change the orientation of an element.
	 */
	void Change_Orientation(void);
	
	/*!
	 * \brief Set the domain element which shares a face with the boundary element.
	 * \param[in] val_domainelement - Global index of the element.
	 */
	void SetDomainElement(unsigned long val_domainelement);
	
	/*!
	 * \brief Get the domain element which shares a face with the boundary element.
	 * \return Domain element which shares a face with the boundary element.
	 */
	unsigned long GetDomainElement(void);
};

/*!
 * \class CRectangle
 * \brief Class for rectangle element definition.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CRectangle : public CPrimalGrid {
private:
	static unsigned short Faces[4][2];			/*!< \brief Matrix to store the local nodes of all the faces. */
	static unsigned short Neighbor_Nodes[4][2];	/*!< \brief Neighbor to a nodes in the element. */
	static unsigned short nNodesFace[4];		/*!< \brief Number of nodes of each face of the element. */
	static unsigned short nNeighbor_Nodes[4];	/*!< \brief Number of Neighbor to a nodes in the element. */
	static unsigned short nFaces;				/*!< \brief Number of faces of the element. */
	static unsigned short nNodes;				/*!< \brief Number of nodes of the element. */
	static unsigned short VTK_Type;				/*!< \brief Type of element using VTK nomenclature. */
	static unsigned short maxNodesFace;			/*!< \brief Maximum number of nodes for a face. */
	static unsigned short nNeighbor_Elements;	/*!< \brief Number of neighbor elements. */
  
public:
  
	/*!
	 * \brief Constructor using the nodes and index.
	 * \param[in] val_iPoint - Index of the 1st point read from the grid file.
	 * \param[in] val_jPoint - Index of the 2nd point read from the grid file.
	 * \param[in] val_point_2 - Index of the 3th point read from the grid file.
	 * \param[in] val_point_3 - Index of the 4th point read from the grid file.
	 * \param[in] val_nDim - Number of dimension of the problem (2D or 3D).
	 */
	CRectangle(unsigned long val_iPoint, unsigned long val_jPoint,
             unsigned long val_point_2, unsigned long val_point_3, unsigned short val_nDim);
  
  /*!
	 * \brief Destructor of the class.
	 */
	~CRectangle(void);
  
	/*!
	 * \brief Get the nodes shared by the triangle.
	 * \param[in] val_node - Local (to the rectangle) index of the node (a rectangle has 4 nodes).
	 * \return Global index of the triangle node.
	 */
	unsigned long GetNode(unsigned short val_node);
  
	/*!
	 * \brief Get the face index of and element.
	 * \param[in] val_face - Local index of the face.
	 * \param[in] val_index - Local (to the face) index of the nodes that compose the face.
	 * \return Local (to the rectangle) index of the nodes that compose the face.
	 */
	unsigned short GetFaces(unsigned short val_face, unsigned short val_index);
  
	/*!
	 * \brief Get the local index of the neighbors to a node (given the local index).
	 * \param[in] val_node - Local (to the element) index of a node.
	 * \param[in] val_index - Local (to the neighbor nodes of val_node) index of the nodes that are neighbor to val_node.
	 * \return Local (to the rectangle) index of the nodes that are neighbor to val_node.
	 */
	unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index);
  
	/*!
	 * \brief Get the number of neighbors nodes of a node.
	 * \param[in] val_node - Local (to the element) index of a node.
	 * \return Number if neighbors of a node val_node.
	 */
	unsigned short GetnNeighbor_Nodes(unsigned short val_node);
	
	/*!
	 * \brief Get the number of nodes that composes a face of an element.
	 * \param[in] val_face - Local index of the face.
	 * \return Number of nodes that composes a face of an element.
	 */
	unsigned short GetnNodesFace(unsigned short val_face);
	
	/*!
	 * \brief Get the number of nodes of an element.
	 * \return Number of nodes that composes an element.
	 */
	unsigned short GetnNodes(void);
	
	/*!
	 * \brief Get the number of faces of an element.
	 * \return Number of faces of an element.
	 */
	unsigned short GetnFaces(void);
	
	/*!
	 * \brief Get the Maximum number of nodes of a face of an element.
	 * \return Maximum number of nodes of a face of an element.
	 */
	unsigned short GetMaxNodesFace(void);
	
	/*!
	 * \brief Get the type of the element using VTK nomenclature.
	 * \return Type of the element using VTK nomenclature.
	 */
	unsigned short GetVTK_Type(void);
	
	/*!
	 * \brief Get the number of element that are neighbor to this element.
	 * \return Number of neighbor elements.
	 */
	unsigned short GetnNeighbor_Elements(void);
	
	/*!
	 * \brief Change the orientation of an element.
	 */
	void Change_Orientation(void);
	
	/*!
	 * \brief Set the domain element which shares a face with the boundary element.
	 * \param[in] val_domainelement - Global index of the element.
	 */
	void SetDomainElement(unsigned long val_domainelement);
	
	/*!
	 * \brief Get the domain element which shares a face with the boundary element.
	 * \return Domain element which shares a face with the boundary element.
	 */
	unsigned long GetDomainElement(void);
};

/*!
 * \class CTetrahedron
 * \brief Class for tetrahedron element definition.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CTetrahedron : public CPrimalGrid {
private:
	static unsigned short Faces[4][3];			/*!< \brief Matrix to store the local nodes of all the faces. */
	static unsigned short Neighbor_Nodes[4][3];	/*!< \brief Neighbor to a nodes in the element. */
	static unsigned short nNodesFace[4];		/*!< \brief Number of nodes of each face of the element. */
	static unsigned short nNeighbor_Nodes[4];	/*!< \brief Number of Neighbor to a nodes in the element. */
	static unsigned short nFaces;				/*!< \brief Number of faces of the element. */
	static unsigned short nNodes;				/*!< \brief Number of nodes of the element. */
	static unsigned short VTK_Type;				/*!< \brief Type of element using VTK nomenclature. */
	static unsigned short maxNodesFace;			/*!< \brief Maximum number of nodes for a face. */
	static unsigned short nNeighbor_Elements;	/*!< \brief Number of neighbor elements. */
	
public:
  
	/*!
	 * \brief Constructor using the nodes and index.
	 * \param[in] val_iPoint - Index of the 1st point read from the grid file.
	 * \param[in] val_jPoint - Index of the 2nd point read from the grid file.
	 * \param[in] val_point_2 - Index of the 3th point read from the grid file.
	 * \param[in] val_point_3 - Index of the 4th point read from the grid file.
	 */
	CTetrahedron(unsigned long val_iPoint, unsigned long val_jPoint,
               unsigned long val_point_2, unsigned long val_point_3);
  
  /*!
	 * \brief Destructor of the class.
	 */
	~CTetrahedron(void);
  
	/*!
	 * \brief Get the nodes shared by the tetrahedron.
	 * \param[in] val_node - Local (to the tetrahedron) index of the node (a tetrahedron has 4 nodes).
	 * \return Global index of the tetrahedron node.
	 */
	unsigned long GetNode(unsigned short val_node);
	
	/*!
	 * \brief Get the face index of and element.
	 * \param[in] val_face - Local index of the face.
	 * \param[in] val_index - Local (to the face) index of the nodes that compose the face.
	 * \return Local (to the element) index of the nodes that compose the face.
	 */
	unsigned short GetFaces(unsigned short val_face, unsigned short val_index);
	
	/*!
	 * \brief Get the local index of the neighbors to a node (given the local index).
	 * \param[in] val_node - Local (to the element) index of a node.
	 * \param[in] val_index - Local (to the neighbor nodes of val_node) index of the nodes that are neighbor to val_node.
	 * \return Local (to the element) index of the nodes that are neighbor to val_node.
	 */
	unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index);
	
	/*!
	 * \brief Get the number of neighbors nodes of a node.
	 * \param[in] val_node - Local (to the element) index of a node.
	 * \return Number if neighbors of a node val_node.
	 */
	unsigned short GetnNeighbor_Nodes(unsigned short val_node);
  
	/*!
	 * \brief Get the number of nodes that composes a face of an element.
	 * \param[in] val_face - Local index of the face.
	 * \return Number of nodes that composes a face of an element.
	 */
	unsigned short GetnNodesFace(unsigned short val_face);
	
	/*!
	 * \brief Get the number of nodes of an element.
	 * \return Number of nodes that composes an element.
	 */
	unsigned short GetnNodes(void);
	
	/*!
	 * \brief Get the number of faces of an element.
	 * \return Number of faces of an element.
	 */
	unsigned short GetnFaces(void);
	
	/*!
	 * \brief Get the Maximum number of nodes of a face of an element.
	 * \return Maximum number of nodes of a face of an element.
	 */
	unsigned short GetMaxNodesFace(void);
	
	/*!
	 * \brief Get the type of the element using VTK nomenclature.
	 * \return Type of the element using VTK nomenclature.
	 */
	unsigned short GetVTK_Type(void);
	
	/*!
	 * \brief Get the number of element that are neighbor to this element.
	 * \return Number of neighbor elements.
	 */
	unsigned short GetnNeighbor_Elements(void);
	
	/*!
	 * \brief Change the orientation of an element.
	 */
	void Change_Orientation(void);
};

/*!
 * \class CHexahedron
 * \brief Class for hexahedron element definition.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CHexahedron : public CPrimalGrid {
private:
	static unsigned short Faces[6][4];			/*!< \brief Matrix to store the local nodes of all the faces. */
	static unsigned short Neighbor_Nodes[8][3];	/*!< \brief Neighbor to a nodes in the element. */
	static unsigned short nNodesFace[6];		/*!< \brief Number of nodes of each face of the element. */
	static unsigned short nNeighbor_Nodes[8];	/*!< \brief Number of Neighbor to a nodes in the element. */
	static unsigned short nFaces;				/*!< \brief Number of faces of the element. */
	static unsigned short nNodes;				/*!< \brief Number of nodes of the element. */
	static unsigned short VTK_Type;				/*!< \brief Type of element using VTK nomenclature. */
	static unsigned short maxNodesFace;			/*!< \brief Maximum number of nodes for a face. */
	static unsigned short nNeighbor_Elements;	/*!< \brief Number of neighbor elements. */
  
public:
	
	/*!
	 * \brief Constructor using the nodes and index.
	 * \param[in] val_iPoint - Index of the 1st point read from the grid file.
	 * \param[in] val_jPoint - Index of the 2nd point read from the grid file.
	 * \param[in] val_point_2 - Index of the 3th point read from the grid file.
	 * \param[in] val_point_3 - Index of the 4th point read from the grid file.
	 * \param[in] val_point_4 - Index of the 5td point read from the grid file.
	 * \param[in] val_point_5 - Index of the 6th point read from the grid file.
	 * \param[in] val_point_6 - Index of the 7th point read from the grid file.
	 * \param[in] val_point_7 - Index of the 8th point read from the grid file.
	 */
	CHexahedron(unsigned long val_iPoint, unsigned long val_jPoint,
              unsigned long val_point_2, unsigned long val_point_3,
              unsigned long val_point_4, unsigned long val_point_5,
              unsigned long val_point_6, unsigned long val_point_7);
  
  /*!
	 * \brief Destructor of the class.
	 */
	~CHexahedron(void);
  
	/*!
	 * \brief Get the nodes shared by the triangle.
	 * \param[in] val_node - Local (to the triangle) index of the node (a triangle has 3 nodes).
	 * \return Global index of the triangle node.
	 */
	unsigned long GetNode(unsigned short val_node);
	
	/*!
	 * \brief Get the face index of and element.
	 * \param[in] val_face - Local index of the face.
	 * \param[in] val_index - Local (to the face) index of the nodes that compose the face.
	 * \return Local (to the element) index of the nodes that compose the face.
	 */
	unsigned short GetFaces(unsigned short val_face, unsigned short val_index);
	
	/*!
	 * \brief Get the local index of the neighbors to a node (given the local index).
	 * \param[in] val_node - Local (to the element) index of a node.
	 * \param[in] val_index - Local (to the neighbor nodes of val_node) index of the nodes
	 *            that are neighbor to val_node.
	 * \return Local (to the element) index of the nodes that are neighbor to val_node.
	 */
	unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index);
  
	/*!
	 * \brief Get the number of neighbors nodes of a node.
	 * \param[in] val_node - Local (to the element) index of a node.
	 * \return Number if neighbors of a node val_node.
	 */
	unsigned short GetnNeighbor_Nodes(unsigned short val_node);
  
	/*!
	 * \brief Get the number of nodes that composes a face of an element.
	 * \param[in] val_face - Local index of the face.
	 * \return Number of nodes that composes a face of an element.
	 */
	unsigned short GetnNodesFace(unsigned short val_face);
	
	/*!
	 * \brief Get the number of nodes of an element.
	 * \return Number of nodes that composes an element.
	 */
	unsigned short GetnNodes(void);
	
	/*!
	 * \brief Get the number of faces of an element.
	 * \return Number of faces of an element.
	 */
	unsigned short GetnFaces(void);
	
	/*!
	 * \brief Get the Maximum number of nodes of a face of an element.
	 * \return Maximum number of nodes of a face of an element.
	 */
	unsigned short GetMaxNodesFace(void);
	
	/*!
	 * \brief Get the type of the element using VTK nomenclature.
	 * \return Type of the element using VTK nomenclature.
	 */
	unsigned short GetVTK_Type(void);
	
	/*!
	 * \brief Get the number of element that are neighbor to this element.
	 * \return Number of neighbor elements.
	 */
	unsigned short GetnNeighbor_Elements(void);
	
	/*!
	 * \brief Change the orientation of an element.
	 */
	void Change_Orientation(void);
};

/*!
 * \class CWedge
 * \brief Class for wedge element definition.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CWedge : public CPrimalGrid {
private:
	static unsigned short Faces[5][4];			/*!< \brief Matrix to store the local nodes of all the faces. */
	static unsigned short Neighbor_Nodes[6][3];	/*!< \brief Neighbor to a nodes in the element. */
	static unsigned short nNodesFace[5];		/*!< \brief Number of nodes of each face of the element. */
	static unsigned short nNeighbor_Nodes[6];	/*!< \brief Number of Neighbor to a nodes in the element. */
	static unsigned short nFaces;				/*!< \brief Number of faces of the element. */
	static unsigned short nNodes;				/*!< \brief Number of nodes of the element. */
	static unsigned short VTK_Type;				/*!< \brief Type of element using VTK nomenclature. */
	static unsigned short maxNodesFace;			/*!< \brief Maximum number of nodes for a face. */
	static unsigned short nNeighbor_Elements;	/*!< \brief Number of neighbor elements. */
	
public:
	
	/*!
	 * \brief Constructor using the nodes and index.
	 * \param[in] val_point_0 - Index of the 1st point read from the grid file.
	 * \param[in] val_point_1 - Index of the 2nd point read from the grid file.
	 * \param[in] val_point_2 - Index of the 3th point read from the grid file.
	 * \param[in] val_point_3 - Index of the 4th point read from the grid file.
	 * \param[in] val_point_4 - Index of the 5th point read from the grid file.
	 * \param[in] val_point_5 - Index of the 6th point read from the grid file.
	 */
	CWedge(unsigned long val_point_0, unsigned long val_point_1,
         unsigned long val_point_2, unsigned long val_point_3,
         unsigned long val_point_4, unsigned long val_point_5);
  
  /*!
	 * \brief Destructor of the class.
	 */
	~CWedge(void);
  
	/*!
	 * \brief Get the nodes shared by the triangle.
	 * \param[in] val_node - Local (to the triangle) index of the node (a wedge has 6 nodes).
	 * \return Global index of the wedge node.
	 */
	unsigned long GetNode(unsigned short val_node);
	
	/*!
	 * \brief Get the face index of and element.
	 * \param[in] val_face - Local index of the face.
	 * \param[in] val_index - Local (to the face) index of the nodes that compose the face.
	 * \return Local (to the element) index of the nodes that compose the face.
	 */
	unsigned short GetFaces(unsigned short val_face, unsigned short val_index);
	
	/*!
	 * \brief Get the local index of the neighbors to a node (given the local index).
	 * \param[in] val_node - Local (to the element) index of a node.
	 * \param[in] val_index - Local (to the neighbor nodes of val_node) index of the nodes that are neighbor to val_node.
	 * \return Local (to the element) index of the nodes that are neighbor to val_node.
	 */
	unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index);
	
	/*!
	 * \brief Get the number of neighbors nodes of a node.
	 * \param[in] val_node - Local (to the element) index of a node.
	 * \return Number if neighbors of a node val_node.
	 */
	unsigned short GetnNeighbor_Nodes(unsigned short val_node);
	
	/*!
	 * \brief Get the number of nodes that composes a face of an element.
	 * \param[in] val_face - Local index of the face.
	 * \return Number of nodes that composes a face of an element.
	 */
	unsigned short GetnNodesFace(unsigned short val_face);
	
	/*!
	 * \brief Get the number of nodes of an element.
	 * \return Number of nodes that composes an element.
	 */
	unsigned short GetnNodes(void);
	
	/*!
	 * \brief Get the number of faces of an element.
	 * \return Number of faces of an element.
	 */
	unsigned short GetnFaces(void);
	
	/*!
	 * \brief Get the Maximum number of nodes of a face of an element.
	 * \return Maximum number of nodes of a face of an element.
	 */
	unsigned short GetMaxNodesFace(void);
	
	/*!
	 * \brief Get the type of the element using VTK nomenclature.
	 * \return Type of the element using VTK nomenclature.
	 */
	unsigned short GetVTK_Type(void);
	
	/*!
	 * \brief Get the number of element that are neighbor to this element.
	 * \return Number of neighbor elements.
	 */
	unsigned short GetnNeighbor_Elements(void);
	
	/*!
	 * \brief Change the orientation of an element.
	 */
	void Change_Orientation(void);
};

/*!
 * \class CPyramid
 * \brief Class for pyramid element definition.
 * \author F. Palacios.
 * \version 2.0.6
 */
class CPyramid : public CPrimalGrid {
private:
	static unsigned short Faces[5][4];			/*!< \brief Matrix to store the local nodes of all the faces. */
	static unsigned short Neighbor_Nodes[5][4];	/*!< \brief Neighbor to a nodes in the element. */
	static unsigned short nNodesFace[5];		/*!< \brief Number of nodes of each face of the element. */
	static unsigned short nNeighbor_Nodes[5];	/*!< \brief Number of Neighbor to a nodes in the element. */
	static unsigned short nFaces;				/*!< \brief Number of faces of the element. */
	static unsigned short nNodes;				/*!< \brief Number of nodes of the element. */
	static unsigned short VTK_Type;				/*!< \brief Type of element using VTK nomenclature. */
	static unsigned short maxNodesFace;			/*!< \brief Maximum number of nodes for a face. */
	static unsigned short nNeighbor_Elements;	/*!< \brief Number of neighbor elements. */
	
public:
	
	/*!
	 * \brief Constructor using the nodes and index.
	 * \param[in] val_iPoint - Index of the 1st point read from the grid file.
	 * \param[in] val_jPoint - Index of the 2nd point read from the grid file.
	 * \param[in] val_point_2 - Index of the 3th point read from the grid file.
	 * \param[in] val_point_3 - Index of the 4th point read from the grid file.
	 * \param[in] val_point_4 - Index of the 5th point read from the grid file.
	 */
	CPyramid(unsigned long val_iPoint, unsigned long val_jPoint,
           unsigned long val_point_2, unsigned long val_point_3,
           unsigned long val_point_4);
  
  /*!
	 * \brief Destructor of the class.
	 */
	~CPyramid(void);
  
	/*!
	 * \brief Get the nodes shared by the pyramid.
	 * \param[in] val_node - Local (to the pyramid) index of the node (a pyramid has 3 nodes).
	 * \return Global index of the pyramid node.
	 */
	unsigned long GetNode(unsigned short val_node);
	
	/*!
	 * \brief Get the face index of and element.
	 * \param[in] val_face - Local index of the face.
	 * \param[in] val_index - Local (to the face) index of the nodes that compose the face.
	 * \return Local (to the element) index of the nodes that compose the face.
	 */
	unsigned short GetFaces(unsigned short val_face, unsigned short val_index);
	
	/*!
	 * \brief Get the local index of the neighbors to a node (given the local index).
	 * \param[in] val_node - Local (to the element) index of a node.
	 * \param[in] val_index - Local (to the neighbor nodes of val_node) index of the nodes that are neighbor to val_node.
	 * \return Local (to the element) index of the nodes that are neighbor to val_node.
	 */
	unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index);
	
	/*!
	 * \brief Get the number of neighbors nodes of a node.
	 * \param[in] val_node - Local (to the element) index of a node.
	 * \return Number if neighbors of a node val_node.
	 */
	unsigned short GetnNeighbor_Nodes(unsigned short val_node);
	
	/*!
	 * \brief Get the number of nodes that composes a face of an element.
	 * \param[in] val_face - Local index of the face.
	 * \return Number of nodes that composes a face of an element.
	 */
	unsigned short GetnNodesFace(unsigned short val_face);
	
	/*!
	 * \brief Get the number of nodes of an element.
	 * \return Number of nodes that composes an element.
	 */
	unsigned short GetnNodes(void);
	
	/*!
	 * \brief Get the number of faces of an element.
	 * \return Number of faces of an element.
	 */
	unsigned short GetnFaces(void);
	
	/*!
	 * \brief Get the Maximum number of nodes of a face of an element.
	 * \return Maximum number of nodes of a face of an element.
	 */
	unsigned short GetMaxNodesFace(void);
	
	/*!
	 * \brief Get the type of the element using VTK nomenclature.
	 * \return Type of the element using VTK nomenclature.
	 */
	unsigned short GetVTK_Type(void);
	
	/*!
	 * \brief Get the number of element that are neighbor to this element.
	 * \return Number of neighbor elements.
	 */
	unsigned short GetnNeighbor_Elements(void);
	
	/*!
	 * \brief Change the orientation of an element.
	 */
	void Change_Orientation(void);
};

#include "primal_grid_structure.inl"
