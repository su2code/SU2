/*!
 * \file primal_grid_structure.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>primal_grid_structure.cpp</i> file.
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

#include "./mpi_structure.hpp"

#include <iostream>
#include <vector>
#include <cstdlib>

#include "dual_grid_structure.hpp"
#include "config_structure.hpp"

using namespace std;

/*!
 * \class CPrimalGrid
 * \brief Class to define the numerical primal grid.
 * \author F. Palacios, T. Economon
 */
class CPrimalGrid {
protected:
  unsigned long *Nodes;         /*!< \brief Vector to store the global nodes of an element. */
  unsigned long GlobalIndex;    /*!< \brief The global index of an element. */
  long *Neighbor_Elements;      /*!< \brief Vector to store the elements surronding an element. */
  short *PeriodIndexNeighbors;  /*!< \brief Vector to store the periodic index of a neighbor.
                                            A -1 indicates no periodic transformation to the neighbor. */
  su2double *Coord_CG;             /*!< \brief Coordinates of the center-of-gravity of the element. */
  su2double **Coord_FaceElems_CG;  /*!< \brief Coordinates of the center-of-gravity of the face of the
                                               elements. */
  static unsigned short nDim;	   /*!< \brief Dimension of the element (2D or 3D) useful for triangles,
                                               quadrilateral and edges. */
  unsigned long DomainElement;	   /*!< \brief Only for boundaries, in this variable the 3D elements which
                                               correspond with a boundary element is stored. */
  bool Divide;                     /*!< \brief Marker used to know if we are going to divide this element
                                               in the adaptation proccess. */
  su2double Volume;                /*!< \brief Volume of the element. */
  bool *JacobianFaceIsConstant;    /*!< \brief Whether or not the Jacobian of the faces can be considered
                                               constant in the transformation to the standard element. */
  bool *ElementOwnsFace;    /*!< \brief Whether or not the element owns the face. */
  su2double LenScale;       /*!< \brief Length scale of the element. */
  unsigned short TimeLevel; /*!< \brief Time level of the element for time accurate local time stepping. */
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
   * \brief Make available the length scale of the element.
   * \return The length scale of the element.
   */
  su2double GetLengthScale(void);

  /*!
   * \brief Set the length scale of the element.
   * \param[in] val_lenScale - Length scale of the element.
   */
  void SetLengthScale(su2double val_lenScale);

  /*!
   * \brief Make available the time level of the element.
   * \return The time level of the element.
   */
  unsigned short GetTimeLevel(void);

  /*!
   * \brief Set the time level of the element.
   * \param[in] val_timeLevel - Time level of the element.
   */
  void SetTimeLevel(unsigned short val_timeLevel);

  /*!
  * \brief Get the boolean to indicate whether or not this element owns the face
           between the current and the adjacent element with index val_face.
  * \param[in] val_face - Local index of the face.
  * \return   Boolean to indicate whether or not the face is owned by this element.
  */
  bool GetOwnerFace(unsigned short val_face);

  /*!
  * \brief Set the boolean to indicate whether or not this element owns the face
           between the current and the adjacent element with index val_face.
  * \param[in] val_owner - Whether or not this element owns the face.
  * \param[in] val_face  - Local index of the face.
  */
  void SetOwnerFace(bool val_owner, unsigned short val_face);

  /*!
  * \brief Get the index of the periodic transformation to the neighboring element.
  * \param[in] val_face - Local index of the face.
  * \return   Index of the periodic transformation to the neighboring element.
  */
  short GetPeriodicIndex(unsigned short val_face);

  /*!
  * \brief Set the index of the periodic transformation to the neighboring element.
  * \param[in] val_periodic - Index of the periodic marker to which the face belongs.
  * \param[in] val_face     - Local index of the face.
  */
  void SetPeriodicIndex(unsigned short val_periodic, unsigned short val_face);

  /*!
  * \brief Get whether or not the Jacobian of the given face is considered constant.
  * \param[in] val_face - Local index of the face.
  * \return  Whether or not the Jacobian of the face is considered constant.
  */
  bool GetJacobianConstantFace(unsigned short val_face);

  /*!
  * \brief Set whether or not the Jacobian of the given face is considered constant.
  * \param[in] val_JacFaceIsConstant - Boolean to indicate whether or not the Jacobian is constant.
  * \param[in] val_face              - Local index of the face.
  */
  void SetJacobianConstantFace(bool val_JacFaceIsConstant, unsigned short val_face);

  /*!
   * \brief Set the center of gravity of an element (including edges).
   * \param[in] val_coord - Coordinates of the element.
   */
  void SetCoord_CG(su2double **val_coord);

  /*!
   * \brief Get the center of gravity of an element (including edges).
   * \param[in] val_dim - Coordinate of the center of gravity.
   * \return Coordinates of the center of gravity.
   */
  su2double GetCG(unsigned short val_dim);

  /*!
   * \brief Set the center of gravity of an element (including edges).
   * \param[in] val_coord - Coordinates of the element.
   */
  void SetVolume(su2double val_volume);
  
  /*!
   * \brief Get the center of gravity of an element (including edges).
   * \param[in] val_dim - Coordinate of the center of gravity.
   * \return Coordinates of the center of gravity.
   */
  su2double GetVolume(void);

  /*!
   * \brief Get the CG of a face of an element.
   * \param[in] val_face - Local index of the face.
   * \param[in] val_dim - Coordinate of the center of gravity.
   * \return Coordinates of the center of gravity.
   */
  su2double GetFaceCG(unsigned short val_face, unsigned short val_dim);

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
  * \brief Initialize the array, which stores whether or not the faces have a constant Jacobian.
  * \param[in] val_nFaces - Number of faces for which Jacobians must be initialized.
  */
  void InitializeJacobianConstantFaces(unsigned short val_nFaces);

  /*!
  * \brief Initialize the information about the neighboring elements.
  * \param[in] val_nFaces - Number of faces for which neighboring information must be initialized.
  */
  void InitializeNeighbors(unsigned short val_nFaces);

  /*!
  * \brief A virtual member.
  * \param[in] val_color - New color of the element.
  */
  virtual void SetColor(unsigned long val_color);

  /*!
  * \brief A virtual member.
  * \return The color of the element in the partitioning.
  */
  virtual unsigned long GetColor(void);
  
  /*!
   * \brief Get the element global index in a parallel computation.
   * \return Global index of the element in a parallel computation.
   */
  unsigned long GetGlobalIndex(void);
  
  /*!
   * \brief Set the global index for an element in a parallel computation.
   * \return Global index of an element in a parallel computation.
   */
  void SetGlobalIndex(unsigned long val_globalindex);

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
   * \param[in] val_node - Local index of a node.
   * \param[in] val_point - Point associated to the node.
   */
  virtual void SetNode(unsigned short val_node, unsigned long val_point);

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

  /*!
  * \brief Virtual function, that must be overwritten by the derived class, if needed.
  * \param[out] nFaces         - Number of faces of this element.
  * \param[out] nPointsPerFace - Number of corner points for each of the faces.
  * \param[out] faceConn       - Global IDs of the corner points of the faces.
  */
  virtual void GetCornerPointsAllFaces(unsigned short &nFaces,
                                       unsigned short nPointsPerFace[],
                                       unsigned long  faceConn[6][4]);

  /*!
  * \brief Virtual function to make available the global ID of this element.
  * \return The global ID of this element.
  */
  virtual unsigned long GetGlobalElemID(void);

  /*!
  * \brief Virtual function to make available the global offset of the solution DOFs.
  * \return The global offset of the solution DOFs.
  */
  virtual unsigned long GetGlobalOffsetDOFsSol(void);

  /*!
  * \brief Virtual function to make available the polynomial degree of the grid.
  * \return The polynomial degree of the grid.
  */
  virtual unsigned short GetNPolyGrid(void);

  /*!
  * \brief Virtual function to make available the polynomial degree of the solution.
  * \return The polynomial degree of the solution.
  */
  virtual unsigned short GetNPolySol(void);

  /*!
  * \brief Virtual function to make available the number of DOFs of the grid in the element.
  * \return The number of DOFs of the Grid in the element.
  */
  virtual unsigned short GetNDOFsGrid(void);

  /*!
  * \brief Virtual function to make available the number of DOFs of the solution in the element.
  * \return The number of DOFs of the solution in the element.
  */
  virtual unsigned short GetNDOFsSol(void);

  /*!
  * \brief Virtual function to get whether or not the Jacobian is considered constant.
  * \return True if the Jacobian is (almost) constant and false otherwise.
  */
  virtual bool GetJacobianConsideredConstant(void);

  /*!
  * \brief Virtual function to set the value of JacobianConsideredConstant.
  * \param[in] val_JacobianConsideredConstant - The value to be set for JacobianConsideredConstant.
  */
  virtual void SetJacobianConsideredConstant(bool val_JacobianConsideredConstant);

  /*!
  * \brief Virtual function to correct the offset of the global DOFs.
  * \param[in] val_offsetRank - The offset that must be added for this rank.
  */
  virtual void AddOffsetGlobalDOFs(const unsigned long val_offsetRank);

  /*!
  * \brief Virtual function to add the given donor ID to the donor elements for the wall function treatment.
  * \param[in] donorElement - Element to be added to donor elements.
  */
  virtual void AddDonorWallFunctions(const unsigned long donorElement);

  /*!
  * \brief Virtual function to make available the number of donor elements for the wall function treatment.
  * \return The number of donor elements.
  */
  virtual unsigned short GetNDonorsWallFunctions(void);

  /*!
  * \brief Virtual function to make available the pointer to the vector for the donor elements
           for the wall function treatment.
  * \return The pointer to the data of donorElementsWallFunctions.
  */
  virtual unsigned long *GetDonorsWallFunctions(void);

  /*!
  * \brief Virtual function to set the global ID's of the donor elements for the wall function treatment.
  * \param[in] donorElements - Vector, which contain the donor elements.
  */
  virtual void SetDonorsWallFunctions(const vector<unsigned long> &donorElements);

  /*!
  * \brief Virtual function to remove the multiple donors for the wall function treatment.
  */
  virtual void RemoveMultipleDonorsWallFunctions(void);
};

/*!
 * \class CVertexMPI
 * \brief Class for vertex element definition. This kind
 *        of element is used in the parallelization stuff.
 * \author F. Palacios
 */
class CVertexMPI : public CPrimalGrid {
private:
  static unsigned short nFaces;				/*!< \brief Number of faces of the element. */
  static unsigned short nNodes;				/*!< \brief Number of nodes of the element. */
  static unsigned short VTK_Type;			/*!< \brief Type of element using VTK nomenclature. */
  unsigned short Rotation_Type;			    /*!< \brief Definition of the rotation, traslation of the
                                                                    solution at the vertex. */
  static unsigned short maxNodesFace;		/*!< \brief Maximum number of nodes for a face. */
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
   * \brief Set the point associated at a node.
   * \param[in] val_node - Local index of a node.
   * \param[in] val_point - Point associated to the node.
   */
  void SetNode(unsigned short val_node, unsigned long val_point);
  
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
 * \author F. Palacios
 */
class CLine : public CPrimalGrid {
private:
  static unsigned short Faces[1][2];		   /*!< \brief Matrix to store the local nodes of all the faces. */
  static unsigned short Neighbor_Nodes[2][1];  /*!< \brief Neighbor to a nodes in the element. */
  static unsigned short nNodesFace[1];		   /*!< \brief Number of nodes of each face of the element. */
  static unsigned short nNeighbor_Nodes[2];	   /*!< \brief Number of Neighbor to a nodes in the element. */
  static unsigned short nFaces;				   /*!< \brief Number of faces of the element. */
  static unsigned short nNodes;				   /*!< \brief Number of nodes of the element. */
  static unsigned short VTK_Type;			   /*!< \brief Type of element using VTK nomenclature. */
  static unsigned short maxNodesFace;		   /*!< \brief Maximum number of nodes for a face. */
  static unsigned short nNeighbor_Elements;	   /*!< \brief Number of Neighbor_Elements. */
  
public:
  
  /*!
   * \brief Constructor using the nodes and index.
   * \param[in] val_point_0 - Index of the 1st triangle point read from the grid file.
   * \param[in] val_point_1 - Index of the 2nd triangle point read from the grid file.
   * \param[in] val_nDim - Number of dimension of the problem (2D or 3D).
   */
  CLine(unsigned long val_point_0, unsigned long val_point_1, unsigned short val_nDim);
  
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
   * \brief Set the point associated at a node.
   * \param[in] val_node - Local index of a node.
   * \param[in] val_point - Point associated to the node.
   */
  void SetNode(unsigned short val_node, unsigned long val_point);
  
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
 * \author F. Palacios
 */
class CTriangle : public CPrimalGrid {
private:
  static unsigned short Faces[3][2];		   /*!< \brief Matrix to store the local nodes of all the faces. */
  static unsigned short Neighbor_Nodes[3][2];  /*!< \brief Neighbor to a nodes in the element. */
  static unsigned short nNodesFace[3];		   /*!< \brief Number of nodes of each face of the element. */
  static unsigned short nNeighbor_Nodes[3];	   /*!< \brief Number of Neighbor to a nodes in the element. */
  static unsigned short nFaces;				   /*!< \brief Number of faces of the element. */
  static unsigned short nNodes;				   /*!< \brief Number of nodes of the element. */
  static unsigned short VTK_Type;			   /*!< \brief Type of element using VTK nomenclature. */
  static unsigned short maxNodesFace;	       /*!< \brief Maximum number of nodes for a face. */
  static unsigned short nNeighbor_Elements;	   /*!< \brief Number of Neighbor_Elements. */
  
public:

  /*!
   * \brief Constructor using the nodes and index.
   * \param[in] val_point_0 - Index of the 1st triangle point read from the grid file.
   * \param[in] val_point_1 - Index of the 2nd triangle point read from the grid file.
   * \param[in] val_point_2 - Index of the 3th triangle point read from the grid file.
   * \param[in] val_nDim - Number of dimension of the problem (2D or 3D), be careful a triangle could be 2D or 3D.
   */
  CTriangle(unsigned long val_point_0, unsigned long val_point_1,
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
   * \brief Set the point associated at a node.
   * \param[in] val_node - Local index of a node.
   * \param[in] val_point - Point associated to the node.
   */
  void SetNode(unsigned short val_node, unsigned long val_point);
  
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
 * \class CQuadrilateral
 * \brief Class for quadrilateral element definition.
 * \author F. Palacios
 */
class CQuadrilateral : public CPrimalGrid {
private:
  static unsigned short Faces[4][2];		   /*!< \brief Matrix to store the local nodes of all the faces. */
  static unsigned short Neighbor_Nodes[4][2];  /*!< \brief Neighbor to a nodes in the element. */
  static unsigned short nNodesFace[4];		   /*!< \brief Number of nodes of each face of the element. */
  static unsigned short nNeighbor_Nodes[4];	   /*!< \brief Number of Neighbor to a nodes in the element. */
  static unsigned short nFaces;				   /*!< \brief Number of faces of the element. */
  static unsigned short nNodes;				   /*!< \brief Number of nodes of the element. */
  static unsigned short VTK_Type;			   /*!< \brief Type of element using VTK nomenclature. */
  static unsigned short maxNodesFace;		   /*!< \brief Maximum number of nodes for a face. */
  static unsigned short nNeighbor_Elements;	   /*!< \brief Number of neighbor elements. */
  
public:
  
  /*!
   * \brief Constructor using the nodes and index.
   * \param[in] val_point_0 - Index of the 1st point read from the grid file.
   * \param[in] val_point_1 - Index of the 2nd point read from the grid file.
   * \param[in] val_point_2 - Index of the 3th point read from the grid file.
   * \param[in] val_point_3 - Index of the 4th point read from the grid file.
   * \param[in] val_nDim - Number of dimension of the problem (2D or 3D).
   */
  CQuadrilateral(unsigned long val_point_0, unsigned long val_point_1,
                 unsigned long val_point_2, unsigned long val_point_3, unsigned short val_nDim);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CQuadrilateral(void);
  
  /*!
   * \brief Get the nodes shared by the triangle.
   * \param[in] val_node - Local (to the quadrilateral) index of the node (a quadrilateral has 4 nodes).
   * \return Global index of the triangle node.
   */
  unsigned long GetNode(unsigned short val_node);
  
  /*!
   * \brief Set the point associated at a node.
   * \param[in] val_node - Local index of a node.
   * \param[in] val_point - Point associated to the node.
   */
  void SetNode(unsigned short val_node, unsigned long val_point);
  
  /*!
   * \brief Get the face index of and element.
   * \param[in] val_face - Local index of the face.
   * \param[in] val_index - Local (to the face) index of the nodes that compose the face.
   * \return Local (to the quadrilateral) index of the nodes that compose the face.
   */
  unsigned short GetFaces(unsigned short val_face, unsigned short val_index);
  
  /*!
   * \brief Get the local index of the neighbors to a node (given the local index).
   * \param[in] val_node - Local (to the element) index of a node.
   * \param[in] val_index - Local (to the neighbor nodes of val_node) index of the nodes that are neighbor to val_node.
   * \return Local (to the quadrilateral) index of the nodes that are neighbor to val_node.
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
 * \author F. Palacios
 */
class CTetrahedron : public CPrimalGrid {
private:
  static unsigned short Faces[4][3];		   /*!< \brief Matrix to store the local nodes of all the faces. */
  static unsigned short Neighbor_Nodes[4][3];  /*!< \brief Neighbor to a nodes in the element. */
  static unsigned short nNodesFace[4];		   /*!< \brief Number of nodes of each face of the element. */
  static unsigned short nNeighbor_Nodes[4];	   /*!< \brief Number of Neighbor to a nodes in the element. */
  static unsigned short nFaces;				   /*!< \brief Number of faces of the element. */
  static unsigned short nNodes;				   /*!< \brief Number of nodes of the element. */
  static unsigned short VTK_Type;			   /*!< \brief Type of element using VTK nomenclature. */
  static unsigned short maxNodesFace;		   /*!< \brief Maximum number of nodes for a face. */
  static unsigned short nNeighbor_Elements;	   /*!< \brief Number of neighbor elements. */

public:
  
  /*!
   * \brief Constructor using the nodes and index.
   * \param[in] val_point_0 - Index of the 1st point read from the grid file.
   * \param[in] val_point_1 - Index of the 2nd point read from the grid file.
   * \param[in] val_point_2 - Index of the 3th point read from the grid file.
   * \param[in] val_point_3 - Index of the 4th point read from the grid file.
   */
  CTetrahedron(unsigned long val_point_0, unsigned long val_point_1,
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
   * \brief Set the point associated at a node.
   * \param[in] val_node - Local index of a node.
   * \param[in] val_point - Point associated to the node.
   */
  void SetNode(unsigned short val_node, unsigned long val_point);

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
 * \author F. Palacios
 */
class CHexahedron : public CPrimalGrid {
private:
  static unsigned short Faces[6][4];			/*!< \brief Matrix to store the local nodes of all the faces. */
  static unsigned short Neighbor_Nodes[8][3];	/*!< \brief Neighbor to a nodes in the element. */
  static unsigned short nNodesFace[6];		    /*!< \brief Number of nodes of each face of the element. */
  static unsigned short nNeighbor_Nodes[8];	    /*!< \brief Number of Neighbor to a nodes in the element. */
  static unsigned short nFaces;				    /*!< \brief Number of faces of the element. */
  static unsigned short nNodes;				    /*!< \brief Number of nodes of the element. */
  static unsigned short VTK_Type;				/*!< \brief Type of element using VTK nomenclature. */
  static unsigned short maxNodesFace;			/*!< \brief Maximum number of nodes for a face. */
  static unsigned short nNeighbor_Elements;	    /*!< \brief Number of neighbor elements. */
  
public:

  /*!
   * \brief Constructor using the nodes and index.
   * \param[in] val_point_0 - Index of the 1st point read from the grid file.
   * \param[in] val_point_1 - Index of the 2nd point read from the grid file.
   * \param[in] val_point_2 - Index of the 3th point read from the grid file.
   * \param[in] val_point_3 - Index of the 4th point read from the grid file.
   * \param[in] val_point_4 - Index of the 5td point read from the grid file.
   * \param[in] val_point_5 - Index of the 6th point read from the grid file.
   * \param[in] val_point_6 - Index of the 7th point read from the grid file.
   * \param[in] val_point_7 - Index of the 8th point read from the grid file.
   */
  CHexahedron(unsigned long val_point_0, unsigned long val_point_1,
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
   * \brief Set the point associated at a node.
   * \param[in] val_node - Local index of a node.
   * \param[in] val_point - Point associated to the node.
   */
  void SetNode(unsigned short val_node, unsigned long val_point);

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
 * \class CPrism
 * \brief Class for prism element definition.
 * \author F. Palacios
 */
class CPrism : public CPrimalGrid {
private:
  static unsigned short Faces[5][4];			/*!< \brief Matrix to store the local nodes of all the faces. */
  static unsigned short Neighbor_Nodes[6][3];	/*!< \brief Neighbor to a nodes in the element. */
  static unsigned short nNodesFace[5];		    /*!< \brief Number of nodes of each face of the element. */
  static unsigned short nNeighbor_Nodes[6];	    /*!< \brief Number of Neighbor to a nodes in the element. */
  static unsigned short nFaces;				    /*!< \brief Number of faces of the element. */
  static unsigned short nNodes;				    /*!< \brief Number of nodes of the element. */
  static unsigned short VTK_Type;				/*!< \brief Type of element using VTK nomenclature. */
  static unsigned short maxNodesFace;			/*!< \brief Maximum number of nodes for a face. */
  static unsigned short nNeighbor_Elements;	    /*!< \brief Number of neighbor elements. */

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
  CPrism(unsigned long val_point_0, unsigned long val_point_1,
         unsigned long val_point_2, unsigned long val_point_3,
         unsigned long val_point_4, unsigned long val_point_5);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CPrism(void);
  
  /*!
   * \brief Get the nodes shared by the triangle.
   * \param[in] val_node - Local (to the triangle) index of the node (a prism has 6 nodes).
   * \return Global index of the prism node.
   */
  unsigned long GetNode(unsigned short val_node);
  
  /*!
   * \brief Set the point associated at a node.
   * \param[in] val_node - Local index of a node.
   * \param[in] val_point - Point associated to the node.
   */
  void SetNode(unsigned short val_node, unsigned long val_point);

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
 * \author F. Palacios
 */
class CPyramid : public CPrimalGrid {
private:
  static unsigned short Faces[5][4];			/*!< \brief Matrix to store the local nodes of all the faces. */
  static unsigned short Neighbor_Nodes[5][4];	/*!< \brief Neighbor to a nodes in the element. */
  static unsigned short nNodesFace[5];		    /*!< \brief Number of nodes of each face of the element. */
  static unsigned short nNeighbor_Nodes[5];	    /*!< \brief Number of Neighbor to a nodes in the element. */
  static unsigned short nFaces;				    /*!< \brief Number of faces of the element. */
  static unsigned short nNodes;				    /*!< \brief Number of nodes of the element. */
  static unsigned short VTK_Type;				/*!< \brief Type of element using VTK nomenclature. */
  static unsigned short maxNodesFace;			/*!< \brief Maximum number of nodes for a face. */
  static unsigned short nNeighbor_Elements;	    /*!< \brief Number of neighbor elements. */

public:

  /*!
   * \brief Constructor using the nodes and index.
   * \param[in] val_point_0 - Index of the 1st point read from the grid file.
   * \param[in] val_point_1 - Index of the 2nd point read from the grid file.
   * \param[in] val_point_2 - Index of the 3th point read from the grid file.
   * \param[in] val_point_3 - Index of the 4th point read from the grid file.
   * \param[in] val_point_4 - Index of the 5th point read from the grid file.
   */
  CPyramid(unsigned long val_point_0, unsigned long val_point_1,
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
   * \brief Set the point associated at a node.
   * \param[in] val_node - Local index of a node.
   * \param[in] val_point - Point associated to the node.
   */
  void SetNode(unsigned short val_node, unsigned long val_point);

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
 * \class CPrimalGridFEM
 * \brief Class to define primal grid element for the FEM solver.
 * \version 7.0.0 "Blackbird"
 */
class CPrimalGridFEM : public CPrimalGrid {
private:
  unsigned short VTK_Type;      /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;     /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nPolySol;      /*!< \brief Polynomial degree for the solution of the element. */
  unsigned short nDOFsGrid;     /*!< \brief Number of DOFs for the geometry of the element. */
  unsigned short nDOFsSol;      /*!< \brief Number of DOFs for the solution of the element. */
  unsigned short nFaces;        /*!< \brief Number of faces of the element. */

  unsigned long elemIDGlobal;        /*!< \brief Global element ID of this element. */
  unsigned long offsetDOFsSolGlobal; /*!< \brief Global offset of the solution DOFs of this element. */
  unsigned long color;               /*!< \brief Color of the element in the partitioning strategy. */

  bool JacobianConsideredConstant;   /*!< \brief Whether or not the Jacobian of the transformation to
                                                 is (almost) constant. */

public:

  /*!
  * \brief Constructor of the class.
  */
  CPrimalGridFEM(void);

  /*!
  * \brief Constructor using data to initialize the element.
  * \param[in] val_elemGlobalID - Global element ID of this element.
  * \param[in] val_VTK_Type     - VTK type to indicate the element type
  * \param[in] val_nPolyGrid    - Polynomial degree to describe the geometry of the element.
  * \param[in] val_nPolySol     - Polynomial degree to describe the solution of the element.
  * \param[in] val_nDOFsGrid    - Number of DOFs used to describe the geometry of the element.
  * \param[in] val_nDOFsSol     - Number of DOFs used to describe the solution of the element.
  * \param[in] val_offDOfsSol   - Global offset of the solution DOFs of the element.
  * \param[in] elem_line        - istringstream, which contains the grid node numbers of the element.
  */
  CPrimalGridFEM(unsigned long  val_elemGlobalID, unsigned short val_VTK_Type,
                 unsigned short val_nPolyGrid,    unsigned short val_nPolySol,
                 unsigned short val_nDOFsGrid,    unsigned short val_nDOFsSol,
                 unsigned long  val_offDOfsSol,   istringstream  &elem_line);

  /*!
  * \brief Constructor using data to initialize the element.
  * \param[in] val_elemGlobalID - Global element ID of this element.
  * \param[in] val_VTK_Type     - VTK type to indicate the element type
  * \param[in] val_nPolyGrid    - Polynomial degree to describe the geometry of the element.
  * \param[in] val_nPolySol     - Polynomial degree to describe the solution of the element.
  * \param[in] val_nDOFsGrid    - Number of DOFs used to describe the geometry of the element.
  * \param[in] val_nDOFsSol     - Number of DOFs used to describe the solution of the element.
  * \param[in] val_offDOfsSol   - Global offset of the solution DOFs of the element.
  * \param[in] connGrid         - Array, which contains the grid node numbers of the element.
  */
  CPrimalGridFEM(unsigned long  val_elemGlobalID, unsigned short val_VTK_Type,
                 unsigned short val_nPolyGrid,    unsigned short val_nPolySol,
                 unsigned short val_nDOFsGrid,    unsigned short val_nDOFsSol,
                 unsigned long  val_offDOfsSol,   const unsigned long *connGrid);

  /*!
  * \brief Destructor of the class.
  */
  ~CPrimalGridFEM(void);

  /*!
  * \brief Get the node shared by the element
  * \param[in] val_node - Local (to the element) index of the node.
  * \return Global index of the node.
  */
  unsigned long GetNode(unsigned short val_node);

  /*!
  * \brief Get the number of nodes that composes a face of an element.
  * \param[in] val_face - Local index of the face.
  * \return Number of nodes that composes a face of an element.
  */
  unsigned short GetnNodesFace(unsigned short val_face);

  /*!
  * \brief Get the face index of an element.
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
  * \brief Get the number of neighbors nodes of a node.
  * \param[in] val_node - Local (to the element) index of a node.
  * \return Number if neighbors of a node val_node.
  */
  unsigned short GetnNeighbor_Nodes(unsigned short val_node);

  /*!
  * \brief Change the orientation of an element.
  */
  void Change_Orientation(void);

  /*!
  * \brief Make available the global ID of this element.
  * \return The global ID of this element.
  */
  unsigned long GetGlobalElemID(void);

  /*!
  * \brief Make available the global offset of the solution DOFs of this element.
  * \return The global offset of the solution DOFs.
  */
  unsigned long GetGlobalOffsetDOFsSol(void);

  /*!
  * \brief Get the number of element that are neighbor to this element.
  * \return Number of neighbor elements.
  */
  unsigned short GetnNeighbor_Elements(void);

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
  * \brief Get the polynomial degree of the grid for this element.
  * \return The polynomial degree of the grid.
  */
  unsigned short GetNPolyGrid(void);

  /*!
  * \brief Get the polynomial degree of the solution for this element.
  * \return The polynomial degree of the solution.
  */
  unsigned short GetNPolySol(void);

  /*!
  * \brief Function to make available the number of DOFs of the grid in the element.
  * \return The number of DOFs of the grid in the element.
  */
  unsigned short GetNDOFsGrid(void);

  /*!
  * \brief Function to make available the number of DOFs of the solution in the element.
  * \return The number of DOFs of the solution in the element.
  */
  unsigned short GetNDOFsSol(void);

  /*!
  * \brief Get all the corner points of all the faces of this element. It must be made sure
           that the numbering of the faces is identical to the numbering used for the
           standard elements.
  * \param[out] nFaces         - Number of faces of this element.
  * \param[out] nPointsPerFace - Number of corner points for each of the faces.
  * \param[out] faceConn       - Global IDs of the corner points of the faces.
  */
  void GetCornerPointsAllFaces(unsigned short &numFaces,
                               unsigned short nPointsPerFace[],
                               unsigned long  faceConn[6][4]);

  /*!
  * \brief Static member function to get the local the corner points of all the faces
           of this element. It must be made sure that the numbering of the faces is
           identical to the numbering used for the standard elements.
  * \param[in]  elementType    - Type of the element using the VTK convention.
  * \param[in]  nPoly          - Polynomial degree of the element.
  * \param[in]  nDOFs          - Number of DOFs of the element.
  * \param[out] nFaces         - Number of faces of this element.
  * \param[out] nPointsPerFace - Number of corner points for each of the faces.
  * \param[out] faceConn       - Global IDs of the corner points of the faces.
  */
  static void GetLocalCornerPointsAllFaces(unsigned short elementType,
                                           unsigned short nPoly,
                                           unsigned short nDOFs,
                                           unsigned short &numFaces,
                                           unsigned short nPointsPerFace[],
                                           unsigned long  faceConn[6][4]);
  /*!
  * \brief Function to get whether or not the Jacobian is considered constant.
  * \return True if the Jacobian is (almost) constant and false otherwise.
  */
  bool GetJacobianConsideredConstant(void);

  /*!
  * \brief Set the color of the element.
  * \param[in] val_color - New color of the element.
  */
  void SetColor(unsigned long val_color);

  /*!
  * \brief Get the color of the element for the partitioning.
  * return - The color of the element in the partitioning.
  */
  unsigned long GetColor(void);

  /*!
  * \brief Function to set the value of JacobianConsideredConstant.
  * \param[in] val_JacobianConsideredConstant - The value to be set for JacobianConsideredConstant.
  */
  void SetJacobianConsideredConstant(bool val_JacobianConsideredConstant);

  /*!
  * \brief Function to correct the offset of the global DOFs.
  * \param[in] val_offsetRank - The offset that must be added for this rank.
  */
  void AddOffsetGlobalDOFs(const unsigned long val_offsetRank);
};

/*!
 * \class CPrimalGridBoundFEM
 * \brief Class to define primal grid boundary element for the FEM solver.
 * \version 7.0.0 "Blackbird"
 */
class CPrimalGridBoundFEM : public CPrimalGrid {
private:
  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;    /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nDOFsGrid;    /*!< \brief Number of DOFs for the geometry of the element. */

  unsigned long boundElemIDGlobal;    /*!< \brief Global boundary element ID of this element. */
  bool JacobianConsideredConstant;    /*!< \brief Whether or not the Jacobian of the transformation to
                                                  is (almost) constant. */

  vector<unsigned long> donorElementsWallFunctions; /*!< \brief The global ID's of the donor elements
                                                                for the wall function treatment. */
public:

  /*!
  * \brief Constructor of the class.
  */
  CPrimalGridBoundFEM(void);

  /*!
  * \brief Constructor using data to initialize the boundary element.
  * \param[in] val_elemGlobalID    - Global boundary element ID of this element.
  * \param[in] val_domainElementID - Global ID of the corresponding domain element.
  * \param[in] val_VTK_Type        - VTK type to indicate the element type
  * \param[in] val_nPolyGrid       - Polynomial degree to describe the geometry of the element.
  * \param[in] val_nDOFsGrid       - Number of DOFs used to describe the geometry of the element.
  * \param[in] val_nodes           - Vector, which contains the global node IDs of the element.
  */
  CPrimalGridBoundFEM(unsigned long         val_elemGlobalID,
                      unsigned long         val_domainElementID,
                      unsigned short        val_VTK_Type,
                      unsigned short        val_nPolyGrid,
                      unsigned short        val_nDOFsGrid,
                      vector<unsigned long> &val_nodes);

  /*!
  * \brief Destructor of the class.
  */
  ~CPrimalGridBoundFEM(void);

  /*!
  * \brief Get the node shared by the element
  * \param[in] val_node - Local (to the element) index of the node.
  * \return Global index of the node.
  */
  unsigned long GetNode(unsigned short val_node);

  /*!
  * \brief Get the number of nodes that composes a face of an element.
  * \param[in] val_face - Local index of the face.
  * \return Number of nodes that composes a face of an element.
  */
  unsigned short GetnNodesFace(unsigned short val_face);

  /*!
  * \brief Get the face index of an element.
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
  * \brief Get the number of neighbors nodes of a node.
  * \param[in] val_node - Local (to the element) index of a node.
  * \return Number if neighbors of a node val_node.
  */
  unsigned short GetnNeighbor_Nodes(unsigned short val_node);

  /*!
  * \brief Change the orientation of an element.
  */
  void Change_Orientation(void);

  /*!
  * \brief Get the number of element that are neighbor to this element.
  * \return Number of neighbor elements.
  */
  unsigned short GetnNeighbor_Elements(void);

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
  * \brief Get the polynomial degree of the grid for this element.
  * \return The polynomial degree of the grid.
  */
  unsigned short GetNPolyGrid(void);

  /*!
  * \brief Function to make available the number of DOFs of the grid in the element.
  * \return The number of DOFs of the grid in the element.
  */
  unsigned short GetNDOFsGrid(void);

  /*!
  * \brief Get the corner points of this boundary element.
  * \param[out] nFaces         - Number of faces of this element, i.e. 1.
  * \param[out] nPointsPerFace - Number of corner points for each of the faces.
  * \param[out] faceConn       - Global IDs of the corner points of the faces.
  */
  void GetCornerPointsAllFaces(unsigned short &nFaces,
                               unsigned short nPointsPerFace[],
                               unsigned long  faceConn[6][4]);

  /*!
  * \brief Static member function to get the local the corner points of all the face
           of this element.
  * \param[in]  elementType    - Type of the element using the VTK convention.
  * \param[in]  nPoly          - Polynomial degree of the element.
  * \param[in]  nDOFs          - Number of DOFs of the element.
  * \param[out] nPointsPerFace - Number of corner points of the face.
  * \param[out] faceConn       - Global IDs of the corner points of the face.
  */
  static void GetLocalCornerPointsFace(unsigned short elementType,
                                       unsigned short nPoly,
                                       unsigned short nDOFs,
                                       unsigned short &nPointsPerFace,
                                       unsigned long  faceConn[]);

  /*!
  * \brief Make available the global ID of this element.
  * \return The global ID of this element.
  */
  unsigned long GetGlobalElemID(void);

  /*!
  * \brief Function to get whether or not the Jacobian is considered constant.
  * \return True if the Jacobian is (almost) constant and false otherwise.
  */
  bool GetJacobianConsideredConstant(void);

  /*!
  * \brief Function to set the value of JacobianConsideredConstant.
  * \param[in] val_JacobianConsideredConstant - The value to be set for JacobianConsideredConstant.
  */
  void SetJacobianConsideredConstant(bool val_JacobianConsideredConstant);

  /*!
  * \brief Add the given donor ID to the donor elements for the wall function treatment.
  * \param[in] donorElement - Element to be added to donor elements.
  */
  void AddDonorWallFunctions(const unsigned long donorElement);

  /*!
  * \brief Make available the number of donor elements for the wall function treatment.
  * \return The number of donor elements.
  */
  unsigned short GetNDonorsWallFunctions(void);

  /*!
  * \brief Make available the pointer to the vector for the donor elements
           for the wall function treatment.
  * \return The pointer to the data of donorElementsWallFunctions.
  */
  unsigned long *GetDonorsWallFunctions(void);

  /*!
  * \brief Set the global ID's of the donor elements for the wall function treatment.
  * \param[in] donorElements - Vector, which contain the donor elements.
  */
  void SetDonorsWallFunctions(const vector<unsigned long> &donorElements);

  /*!
  * \brief Function to remove the multiple donors for the wall function treatment.
  */
  void RemoveMultipleDonorsWallFunctions(void);
};

#include "primal_grid_structure.inl"
