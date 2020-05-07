/*!
 * \file CPrimalGrid.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>primal_grid_structure.cpp</i> file.
 * \author F. Palacios
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

#pragma once

#include "../../mpi_structure.hpp"

#include <iostream>
#include <vector>
#include <cstdlib>

#include "../../CConfig.hpp"


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
  static unsigned short nDim;    /*!< \brief Dimension of the element (2D or 3D) useful for triangles,
                                               quadrilateral and edges. */
  unsigned long DomainElement;     /*!< \brief Only for boundaries, in this variable the 3D elements which
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
  inline long GetNeighbor_Elements(unsigned short val_face) { return Neighbor_Elements[val_face]; }

  /*!
   * \brief Set the elements that surround an element.
   * \param[in] val_elem - Global index of the element.
   * \param[in] val_face - Local index of the face.
   */
  inline void SetNeighbor_Elements(unsigned long val_elem, unsigned short val_face) { Neighbor_Elements[val_face] = val_elem; }

  /*!
   * \brief Make available the length scale of the element.
   * \return The length scale of the element.
   */
  inline su2double GetLengthScale(void) { return LenScale; }

  /*!
   * \brief Set the length scale of the element.
   * \param[in] val_lenScale - Length scale of the element.
   */
  inline void SetLengthScale(su2double val_lenScale) { LenScale = val_lenScale; }

  /*!
   * \brief Make available the time level of the element.
   * \return The time level of the element.
   */
  inline unsigned short GetTimeLevel(void) { return TimeLevel; }

  /*!
   * \brief Set the time level of the element.
   * \param[in] val_timeLevel - Time level of the element.
   */
  inline void SetTimeLevel(unsigned short val_timeLevel) { TimeLevel = val_timeLevel; }

  /*!
  * \brief Get the boolean to indicate whether or not this element owns the face
           between the current and the adjacent element with index val_face.
  * \param[in] val_face - Local index of the face.
  * \return   Boolean to indicate whether or not the face is owned by this element.
  */
  inline bool GetOwnerFace(unsigned short val_face) { return ElementOwnsFace[val_face]; }

  /*!
  * \brief Set the boolean to indicate whether or not this element owns the face
           between the current and the adjacent element with index val_face.
  * \param[in] val_owner - Whether or not this element owns the face.
  * \param[in] val_face  - Local index of the face.
  */
  inline void SetOwnerFace(bool val_owner, unsigned short val_face) { ElementOwnsFace[val_face] = val_owner; }

  /*!
  * \brief Get the index of the periodic transformation to the neighboring element.
  * \param[in] val_face - Local index of the face.
  * \return   Index of the periodic transformation to the neighboring element.
  */
  inline short GetPeriodicIndex(unsigned short val_face) {return PeriodIndexNeighbors[val_face];}

  /*!
  * \brief Set the index of the periodic transformation to the neighboring element.
  * \param[in] val_periodic - Index of the periodic marker to which the face belongs.
  * \param[in] val_face     - Local index of the face.
  */
  inline void SetPeriodicIndex(unsigned short val_periodic, unsigned short val_face) {PeriodIndexNeighbors[val_face] = val_periodic; }

  /*!
  * \brief Get whether or not the Jacobian of the given face is considered constant.
  * \param[in] val_face - Local index of the face.
  * \return  Whether or not the Jacobian of the face is considered constant.
  */
  inline bool GetJacobianConstantFace(unsigned short val_face) { return JacobianFaceIsConstant[val_face]; }

  /*!
  * \brief Set whether or not the Jacobian of the given face is considered constant.
  * \param[in] val_JacFaceIsConstant - Boolean to indicate whether or not the Jacobian is constant.
  * \param[in] val_face              - Local index of the face.
  */
  inline void SetJacobianConstantFace(bool val_JacFaceIsConstant, unsigned short val_face) {JacobianFaceIsConstant[val_face] = val_JacFaceIsConstant; }

  /*!
   * \brief Set the center of gravity of an element (including edges).
   * \param[in] val_coord - Coordinates of the element.
   */
  void SetCoord_CG(const su2double* const* val_coord);

  /*!
   * \brief Get the center of gravity of an element (including edges).
   * \param[in] val_dim - Coordinate of the center of gravity.
   * \return Coordinates of the center of gravity.
   */
  inline su2double GetCG(unsigned short val_dim) { return Coord_CG[val_dim]; }

  /*!
   * \brief Set the center of gravity of an element (including edges).
   * \param[in] val_coord - Coordinates of the element.
   */
  inline void SetVolume(su2double val_volume) { Volume = val_volume; }

  /*!
   * \brief Get the center of gravity of an element (including edges).
   * \param[in] val_dim - Coordinate of the center of gravity.
   * \return Coordinates of the center of gravity.
   */
  inline su2double GetVolume(void) { return Volume; }

  /*!
   * \brief Get the CG of a face of an element.
   * \param[in] val_face - Local index of the face.
   * \param[in] val_dim - Coordinate of the center of gravity.
   * \return Coordinates of the center of gravity.
   */
  inline su2double GetFaceCG(unsigned short val_face, unsigned short val_dim) { return Coord_FaceElems_CG[val_face][val_dim]; }

  /*!
   * \brief Get all the neighbors of an element.
   * \return List of all the neighbor of an element.
   */
  void GetAllNeighbor_Elements(void);

  /*!
   * \brief Set that an element must be divided in the adaptation stage.
   * \param[in] val_divide - <code>TRUE</code> if the element must be divided; otherwise <code>FALSE</code>.
   */
  inline void SetDivide (bool val_divide) { Divide = val_divide; }

  /*!
   * \brief Get if an element must be divided in the adaptation stage.
   * \return <code>TRUE</code> if the element must be divided; otherwise <code>FALSE</code>.
   */
  inline bool GetDivide (void) { return Divide; }

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
  inline virtual void SetColor(unsigned long val_color) { }

  /*!
  * \brief A virtual member.
  * \return The color of the element in the partitioning.
  */
  inline virtual unsigned long GetColor(void) { return -1; }

  /*!
   * \brief Get the element global index in a parallel computation.
   * \return Global index of the element in a parallel computation.
   */
  inline unsigned long GetGlobalIndex(void) { return GlobalIndex; }

  /*!
   * \brief Set the global index for an element in a parallel computation.
   * \return Global index of an element in a parallel computation.
   */
  inline void SetGlobalIndex(unsigned long val_globalindex) { GlobalIndex = val_globalindex; }

  /*!
   * \brief A virtual member.
   * \param[in] val_domainelement Index of the domain element which has a face shared by this boundary element.
   */
  inline virtual void SetDomainElement(unsigned long val_domainelement) { DomainElement = val_domainelement; }

  /*!
   * \brief A virtual member.
   * \return Relate the boundary element which a face of a domain element.
   */
  inline virtual unsigned long GetDomainElement(void) { return DomainElement; }

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
  inline virtual unsigned short GetRotation_Type(void) { return 0; }

  /*!
   * \brief A pure virtual member.
   * \param[in] val_rotation_type - Kind of rotation/traslation that must be applied.
   */
  inline virtual void SetRotation_Type(unsigned short val_rotation_type) { }

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
  inline virtual unsigned short GetnNodesFace(unsigned short val_face) { return 0; }

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
  inline virtual void SetNode(unsigned short val_node, unsigned long val_point) { }

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
   * \param[in] val_index - Local (to the neighbor nodes of <i>val_node</i>) index of the nodes that are neighbor to val_node.
   * \return Local index of the nodes that are neighbor to <i>val_node</i>.
   */
  virtual unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) = 0;

  /*!
   * \brief Virtual function, that must be overwritten by the derived class, if needed.
   * \param[out] nFaces         - Number of faces of this element.
   * \param[out] nPointsPerFace - Number of corner points for each of the faces.
   * \param[out] faceConn       - Global IDs of the corner points of the faces.
   */
  inline virtual void GetCornerPointsAllFaces(unsigned short &nFaces,
                                      unsigned short nPointsPerFace[],
                                      unsigned long  faceConn[6][4]) { }

  /*!
   * \brief Virtual function to make available the global ID of this element.
   * \return The global ID of this element.
   */
  inline virtual unsigned long GetGlobalElemID(void) { return 0; }

  /*!
   * \brief Virtual function to make available the global offset of the solution DOFs.
   * \return The global offset of the solution DOFs.
   */
  inline virtual unsigned long GetGlobalOffsetDOFsSol(void) { return 0; }

  /*!
   * \brief Virtual function to make available the polynomial degree of the grid.
   * \return The polynomial degree of the grid.
   */
  inline virtual unsigned short GetNPolyGrid(void) { return 0; }

  /*!
   * \brief Virtual function to make available the polynomial degree of the solution.
   * \return The polynomial degree of the solution.
   */
  inline virtual unsigned short GetNPolySol(void) { return 0; }

  /*!
   * \brief Virtual function to make available the number of DOFs of the grid in the element.
   * \return The number of DOFs of the Grid in the element.
   */
  inline virtual unsigned short GetNDOFsGrid(void) { return 0; }

  /*!
   * \brief Virtual function to make available the number of DOFs of the solution in the element.
   * \return The number of DOFs of the solution in the element.
   */
  inline virtual unsigned short GetNDOFsSol(void) { return 0; }

  /*!
   * \brief Virtual function to get whether or not the Jacobian is considered constant.
   * \return True if the Jacobian is (almost) constant and false otherwise.
   */
  inline virtual bool GetJacobianConsideredConstant(void) { return false; }

  /*!
   * \brief Virtual function to set the value of JacobianConsideredConstant.
   * \param[in] val_JacobianConsideredConstant - The value to be set for JacobianConsideredConstant.
   */
  inline virtual void SetJacobianConsideredConstant(bool val_JacobianConsideredConstant) {}

  /*!
   * \brief Virtual function to correct the offset of the global DOFs.
   * \param[in] val_offsetRank - The offset that must be added for this rank.
   */
  inline virtual void AddOffsetGlobalDOFs(const unsigned long val_offsetRank) {}

  /*!
   * \brief Virtual function to add the given donor ID to the donor elements for the wall function treatment.
   * \param[in] donorElement - Element to be added to donor elements.
   */
  inline virtual void AddDonorWallFunctions(const unsigned long donorElement) {}

  /*!
   * \brief Virtual function to make available the number of donor elements for the wall function treatment.
   * \return The number of donor elements.
   */
  inline  virtual unsigned short GetNDonorsWallFunctions(void) {return 0;}

  /*!
   * \brief Virtual function to make available the pointer to the vector for the donor elements
            for the wall function treatment.
   * \return The pointer to the data of donorElementsWallFunctions.
   */
  inline virtual unsigned long *GetDonorsWallFunctions(void) {return NULL;}

  /*!
   * \brief Virtual function to set the global ID's of the donor elements for the wall function treatment.
   * \param[in] donorElements - Vector, which contain the donor elements.
   */
  inline virtual void SetDonorsWallFunctions(const vector<unsigned long> &donorElements) {}

  /*!
   * \brief Virtual function to remove the multiple donors for the wall function treatment.
   */
  inline virtual void RemoveMultipleDonorsWallFunctions(void) {}
};
