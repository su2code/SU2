/*!
 * \file CPrimalGridFEM.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>CPrimalGridFEM.cpp</i> file.
 * \author F. Palacios
 * \version 7.0.6 "Blackbird"
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

#include "CPrimalGrid.hpp"

/*!
 * \class CPrimalGridFEM
 * \brief Class to define primal grid element for the FEM solver.
 * \version 7.0.6 "Blackbird"
 */
class CPrimalGridFEM final: public CPrimalGrid {
private:
  unsigned short VTK_Type;      /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;     /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nPolySol;      /*!< \brief Polynomial degree for the solution of the element. */
  unsigned short nDOFsGrid;     /*!< \brief Number of DOFs for the geometry of the element. */
  unsigned short nDOFsSol;      /*!< \brief Number of DOFs for the solution of the element. */

  unsigned long elemIDGlobal;        /*!< \brief Global element ID of this element. */
  unsigned long color;               /*!< \brief Color of the element in the partitioning strategy. */

  bool JacobianConsideredConstant;   /*!< \brief Whether or not the Jacobian of the transformation to
                                                 is (almost) constant. */

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
  inline unsigned short GetFaces(unsigned short val_face, unsigned short val_index) override { return -1; }  

  /*!
   * \brief Get the Maximum number of nodes of a face of an element.
   * \return Maximum number of nodes of a face of an element.
   */
  inline unsigned short GetMaxNodesFace(void) override { return -1; }

  /*!
   * \brief Get the number of faces of an element.
   * \return Number of faces of an element.
   */
  inline unsigned short GetnFaces(void) override { return -1; }  

  /*!
   * \brief Get the number of element that are neighbor to this element.
   * \return Number of neighbor elements.
   */
  inline unsigned short GetnNeighbor_Elements(void) override { return -1; }

  /*!
   * \brief Get the number of neighbors nodes of a node.
   * \param[in] val_node - Local (to the element) index of a node.
   * \return Number if neighbors of a node val_node.
   */
  inline unsigned short GetnNeighbor_Nodes(unsigned short val_node) override { return -1; }  

  /*!
   * \brief Get the local index of the neighbors to a node (given the local index).
   * \param[in] val_node - Local (to the element) index of a node.
   * \param[in] val_index - Local (to the neighbor nodes of val_node) index of the nodes that are neighbor to val_node.
   * \return Local (to the element) index of the nodes that are neighbor to val_node.
   */
  inline unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) override { return -1; }

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
  inline unsigned short GetnNodesFace(unsigned short val_face) override { return -1; }

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
private:

  /*!
   * \brief Default constructor, disabled.
   */
  CPrimalGridFEM(void);
};
