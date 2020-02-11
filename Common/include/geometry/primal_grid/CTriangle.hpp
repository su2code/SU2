/*!
 * \file CTriangle.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>CTriangle.cpp</i> file.
 * \author F. Palacios
 * \version 7.0.1 "Blackbird"
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

#include "CPrimalGrid.hpp"

/*!
 * \class CTriangle
 * \brief Class for triangle element definition.
 * \author F. Palacios
 */
class CTriangle final: public CPrimalGrid {
private:
  static unsigned short Faces[3][2];       /*!< \brief Matrix to store the local nodes of all the faces. */
  static unsigned short Neighbor_Nodes[3][2];  /*!< \brief Neighbor to a nodes in the element. */
  static unsigned short nNodesFace[3];       /*!< \brief Number of nodes of each face of the element. */
  static unsigned short nNeighbor_Nodes[3];    /*!< \brief Number of Neighbor to a nodes in the element. */
  static unsigned short nFaces;          /*!< \brief Number of faces of the element. */
  static unsigned short nNodes;          /*!< \brief Number of nodes of the element. */
  static unsigned short VTK_Type;        /*!< \brief Type of element using VTK nomenclature. */
  static unsigned short maxNodesFace;        /*!< \brief Maximum number of nodes for a face. */
  static unsigned short nNeighbor_Elements;    /*!< \brief Number of Neighbor_Elements. */

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
  ~CTriangle(void) override;

  /*!
   * \brief Get the nodes shared by the triangle.
   * \param[in] val_node - Local (to the triangle) index of the node (a triangle has 3 nodes).
   * \return Global index of the triangle node.
   */
  inline unsigned long GetNode(unsigned short val_node) override { return Nodes[val_node]; }

  /*!
   * \brief Set the point associated at a node.
   * \param[in] val_node - Local index of a node.
   * \param[in] val_point - Point associated to the node.
   */
  inline void SetNode(unsigned short val_node, unsigned long val_point) override { Nodes[val_node] = val_point; }

  /*!
   * \brief Get the face index of and element.
   * \param[in] val_face - Local index of the face.
   * \param[in] val_index - Local (to the face) index of the nodes that compose the face.
   * \return Local (to the triangle) index of the nodes that compose the face.
   */
  inline unsigned short GetFaces(unsigned short val_face, unsigned short val_index) override { return Faces[val_face][val_index]; }

  /*!
   * \brief Get the local index of the neighbors to a node (given the local index).
   * \param[in] val_node - Local (to the triangle) index of a node.
   * \param[in] val_index - Local (to the neighbor nodes of val_node) index of the nodes that
   are neighbor to val_node (each face is composed by 3 nodes).
   * \return Local (to the triangle) index of the nodes that are neighbor to val_node.
   */
  inline unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) override { return Neighbor_Nodes[val_node][val_index]; }

  /*!
   * \brief Get the number of neighbors nodes of a node.
   * \param[in] val_node - Local (to the triangle) index of a node.
   * \return Number if neighbors of a node val_node.
   */
  inline unsigned short GetnNeighbor_Nodes(unsigned short val_node) override { return nNeighbor_Nodes[val_node]; }

  /*!
   * \brief Get the number of nodes that composes a face of an element.
   * \param[in] val_face - Local index of the face.
   * \return Number of nodes that composes a face of an element.
   */
  inline unsigned short GetnNodesFace(unsigned short val_face) override { return nNodesFace[val_face]; }

  /*!
   * \brief Get the number of nodes of an element.
   * \return Number of nodes that composes an element.
   */
  inline unsigned short GetnNodes(void) override { return nNodes; }

  /*!
   * \brief Get the number of faces of an element.
   * \return Number of faces of an element.
   */
  inline unsigned short GetnFaces(void) override { return nFaces; }

  /*!
   * \brief Get the Maximum number of nodes of a face of an element.
   * \return Maximum number of nodes of a face of an element.
   */
  inline unsigned short GetMaxNodesFace(void) override { return maxNodesFace; }

  /*!
   * \brief Get the type of the element using VTK nomenclature.
   * \return Type of the element using VTK nomenclature.
   */
  inline unsigned short GetVTK_Type(void) override { return VTK_Type; }

  /*!
   * \brief Get the number of element that are neighbor to this element.
   * \return Number of neighbor elements.
   */
  inline unsigned short GetnNeighbor_Elements(void) override { return nNeighbor_Elements; }

  /*!
   * \brief Change the orientation of an element.
   */
  inline void Change_Orientation(void) override {
    unsigned long Point_0, Point_2;
    Point_0 = Nodes[0];
    Point_2 = Nodes[2];
    Nodes[0] = Point_2;
    Nodes[2] = Point_0;
  }
};
