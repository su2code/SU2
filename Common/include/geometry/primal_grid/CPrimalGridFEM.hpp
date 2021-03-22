/*!
 * \file CPrimalGridFEM.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>CPrimalGridFEM.cpp</i> file.
 * \author F. Palacios
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "CPrimalGrid.hpp"

/*!
 * \class CPrimalGridFEM
 * \brief Class to define primal grid element for the FEM solver.
 * \version 7.1.1 "Blackbird"
 */
class CPrimalGridFEM final: public CPrimalGrid {
private:
  unsigned short VTK_Type;      /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;     /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nPolySol;      /*!< \brief Polynomial degree for the solution of the element. */
  unsigned short nDOFsGrid;     /*!< \brief Number of DOFs for the geometry of the element. */
  unsigned short nDOFsSol;      /*!< \brief Number of DOFs for the solution of the element. */
  unsigned short nFaces;        /*!< \brief Number of faces of the element. */
  unsigned short TimeLevel;     /*!< \brief Time level of the element for time accurate local time stepping. */

  unsigned long elemIDGlobal;   /*!< \brief Global element ID of this element. */
  unsigned long color;          /*!< \brief Color of the element in the partitioning strategy. */

  bool JacobianConsideredConstant; /*!< \brief Whether or not the Jacobian of the transformation
                                               is (almost) constant. */
  bool *JacobianFaceIsConstant;    /*!< \brief Whether or not the Jacobian of the faces can be considered
                                               constant in the transformation to the standard element. */
  bool *ElementOwnsFace;           /*!< \brief Whether or not the element owns the face. */

  short *PeriodIndexNeighbors;     /*!< \brief Vector to store the periodic index of a neighbor.
                                               A -1 indicates no periodic transformation to the neighbor. */

  su2double LenScale;       /*!< \brief Length scale of the element. */
public:

  /*!
   * \brief Constructor using data to initialize the element.
   * \param[in] dataElem - Meta and connectivity data for this element.
   */
  CPrimalGridFEM(const unsigned long *dataElem); 

  /*!
   * \brief Destructor of the class.
   */
  ~CPrimalGridFEM(void) override;

  /*!
   * \brief Change the orientation of an element.
   */
  inline void Change_Orientation(void) override {}

  /*!
   * \brief Get the face index of an element.
   * \param[in] val_face - Local index of the face.
   * \param[in] val_index - Local (to the face) index of the nodes that compose the face.
   * \return Local (to the element) index of the nodes that compose the face.
   */
  inline unsigned short GetFaces(unsigned short val_face, unsigned short val_index) override { return std::numeric_limits<unsigned short>::max(); }  

  /*!
   * \brief Make available the global ID of this element.
   * \return The global ID of this element.
   */
  inline unsigned long GetGlobalElemID(void) override { return elemIDGlobal; }

  /*!
   * \brief Get the Maximum number of nodes of a face of an element.
   * \return Maximum number of nodes of a face of an element.
   */
  inline unsigned short GetMaxNodesFace(void) override { return std::numeric_limits<unsigned short>::max(); }

  /*!
   * \brief Get the number of faces of an element.
   * \return Number of faces of an element.
   */
  inline unsigned short GetnFaces(void) override { return nFaces; }  

  /*!
   * \brief Get the number of element that are neighbor to this element.
   * \return Number of neighbor elements.
   */
  inline unsigned short GetnNeighbor_Elements(void) override { return std::numeric_limits<unsigned short>::max(); }

  /*!
   * \brief Get the number of neighbors nodes of a node.
   * \param[in] val_node - Local (to the element) index of a node.
   * \return Number if neighbors of a node val_node.
   */
  inline unsigned short GetnNeighbor_Nodes(unsigned short val_node) override { return std::numeric_limits<unsigned short>::max(); }  

  /*!
   * \brief Get the local index of the neighbors to a node (given the local index).
   * \param[in] val_node - Local (to the element) index of a node.
   * \param[in] val_index - Local (to the neighbor nodes of val_node) index of the nodes that are neighbor to val_node.
   * \return Local (to the element) index of the nodes that are neighbor to val_node.
   */
  inline unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) override { return std::numeric_limits<unsigned short>::max(); }

  /*!
   * \brief Get the number of nodes of an element.
   * \return Number of nodes that composes an element.
   */
  inline unsigned short GetnNodes(void) override { return nDOFsGrid; }

  /*!
   * \brief Get the number of nodes that composes a face of an element.
   * \param[in] val_face - Local index of the face.
   * \return Number of nodes that composes a face of an element.
   */
  inline unsigned short GetnNodesFace(unsigned short val_face) override { return std::numeric_limits<unsigned short>::max(); }

  /*!
   * \brief Get the node shared by the element
   * \param[in] val_node - Local (to the element) index of the node.
   * \return Global index of the node.
   */
  inline unsigned long GetNode(unsigned short val_node) override { return Nodes[val_node]; }

  /*!
   * \brief Get the type of the element using VTK nomenclature.
   * \return Type of the element using VTK nomenclature.
   */
  inline unsigned short GetVTK_Type(void) override { return VTK_Type; }

  /*!
   * \brief Get the polynomial degree of the grid for this element.
   * \return The polynomial degree of the grid.
   */
  inline unsigned short GetNPolyGrid(void) override { return nPolyGrid; }

  /*!
   * \brief Get the polynomial degree of the solution for this element.
   * \return The polynomial degree of the solution.
   */
  inline unsigned short GetNPolySol(void) override { return nPolySol; }

  /*!
   * \brief Function to make available the number of DOFs of the grid in the element.
   * \return The number of DOFs of the grid in the element.
   */
  inline unsigned short GetNDOFsGrid(void) override { return nDOFsGrid; }

  /*!
   * \brief Function to make available the number of DOFs of the solution in the element.
   * \return The number of DOFs of the solution in the element.
   */
  inline unsigned short GetNDOFsSol(void) override { return nDOFsSol; }

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
                               unsigned long  faceConn[6][4]) override;
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
  inline bool GetJacobianConsideredConstant(void) override { return JacobianConsideredConstant; }

  /*!
   * \brief Set the color of the element.
   * \param[in] val_color - New color of the element.
   */
  inline void SetColor(unsigned long val_color) override { color = val_color; }

  /*!
   * \brief Get the color of the element for the partitioning.
   * return - The color of the element in the partitioning.
   */
  inline unsigned long GetColor(void) override { return color; }

  /*!
   * \brief Function to set the value of JacobianConsideredConstant.
   * \param[in] val_JacobianConsideredConstant - The value to be set for JacobianConsideredConstant.
   */
  inline void SetJacobianConsideredConstant(bool val_JacobianConsideredConstant) override {
    JacobianConsideredConstant = val_JacobianConsideredConstant;
  }

  /*!
   * \brief Make available the length scale of the element.
   * \return The length scale of the element.
   */
  inline su2double GetLengthScale(void) const override { return LenScale; }

  /*!
   * \brief Set the length scale of the element.
   * \param[in] val_lenScale - Length scale of the element.
   */
  inline void SetLengthScale(su2double val_lenScale) override { LenScale = val_lenScale; }

  /*!
   * \brief Make available the time level of the element.
   * \return The time level of the element.
   */
  inline unsigned short GetTimeLevel(void) const override { return TimeLevel; }

  /*!
   * \brief Set the time level of the element.
   * \param[in] val_timeLevel - Time level of the element.
   */
  inline void SetTimeLevel(unsigned short val_timeLevel) override { TimeLevel = val_timeLevel; }

  /*!
   * \brief Get the boolean to indicate whether or not this element owns the face
            between the current and the adjacent element with index val_face.
   * \param[in] val_face - Local index of the face.
   * \return   Boolean to indicate whether or not the face is owned by this element.
   */
  inline bool GetOwnerFace(unsigned short val_face) override { return ElementOwnsFace[val_face]; }

  /*!
   * \brief Set the boolean to indicate whether or not this element owns the face
            between the current and the adjacent element with index val_face.
   * \param[in] val_owner - Whether or not this element owns the face.
   * \param[in] val_face  - Local index of the face.
   */
  inline void SetOwnerFace(bool val_owner, unsigned short val_face) override { ElementOwnsFace[val_face] = val_owner; }

  /*!
   * \brief Get the index of the periodic transformation to the neighboring element.
   * \param[in] val_face - Local index of the face.
   * \return   Index of the periodic transformation to the neighboring element.
   */
  inline short GetPeriodicIndex(unsigned short val_face) override {return PeriodIndexNeighbors[val_face];}

  /*!
   * \brief Set the index of the periodic transformation to the neighboring element.
   * \param[in] val_periodic - Index of the periodic marker to which the face belongs.
   * \param[in] val_face     - Local index of the face.
   */
  inline void SetPeriodicIndex(unsigned short val_periodic, unsigned short val_face) override {
    PeriodIndexNeighbors[val_face] = val_periodic;
  }

  /*!
   * \brief Get whether or not the Jacobian of the given face is considered constant.
   * \param[in] val_face - Local index of the face.
   * \return  Whether or not the Jacobian of the face is considered constant.
   */
  inline bool GetJacobianConstantFace(unsigned short val_face) override { return JacobianFaceIsConstant[val_face]; }

  /*!
   * \brief Set whether or not the Jacobian of the given face is considered constant.
   * \param[in] val_JacFaceIsConstant - Boolean to indicate whether or not the Jacobian is constant.
   * \param[in] val_face              - Local index of the face.
   */
  inline void SetJacobianConstantFace(bool val_JacFaceIsConstant, unsigned short val_face) override {
    JacobianFaceIsConstant[val_face] = val_JacFaceIsConstant;
  }

  /*!
   * \brief Initialize the array, which stores whether or not the faces have a constant Jacobian.
   * \param[in] val_nFaces - Number of faces for which Jacobians must be initialized.
   */
  void InitializeJacobianConstantFaces(unsigned short val_nFaces) override;

  /*!
   * \brief Initialize the information about the neighboring elements.
   * \param[in] val_nFaces - Number of faces for which neighboring information must be initialized.
   */
  void InitializeNeighbors(unsigned short val_nFaces) override;

private:

  /*!
   * \brief Default constructor, disabled.
   */
  CPrimalGridFEM(void);
};
