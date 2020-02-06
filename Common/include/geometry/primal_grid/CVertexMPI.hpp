/*!
 * \file CVertexMPI.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>primal_grid_structure.cpp</i> file.
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
 * \class CVertexMPI
 * \brief Class for vertex element definition. This kind
 *        of element is used in the parallelization stuff.
 * \author F. Palacios
 */
class CVertexMPI final: public CPrimalGrid {
private:
  static unsigned short nFaces;       /*!< \brief Number of faces of the element. */
  static unsigned short nNodes;       /*!< \brief Number of nodes of the element. */
  static unsigned short VTK_Type;     /*!< \brief Type of element using VTK nomenclature. */
  unsigned short Rotation_Type;         /*!< \brief Definition of the rotation, traslation of the
                                                                    solution at the vertex. */
  static unsigned short maxNodesFace;   /*!< \brief Maximum number of nodes for a face. */
  static unsigned short nNeighbor_Elements; /*!< \brief Number of Neighbor_Elements. */

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
  ~CVertexMPI(void) override;

  /*!
   * \brief Get the nodes shared by the line.
   * \param[in] val_node - Local (to the line) index of the node (a line has 2 nodes).
   * \return Global index of the line node.
   */
  inline unsigned long GetNode(unsigned short val_node) override { return Nodes[val_node]; }

  /*!
   * \brief Set the point associated at a node.
   * \param[in] val_node - Local index of a node.
   * \param[in] val_point - Point associated to the node.
   */
  inline void SetNode(unsigned short val_node, unsigned long val_point) override { Nodes[val_node] = val_point; }

  /*!
   * \brief Get the number of nodes of an element.
   * \return Number of nodes that composes an element.
   */
  inline unsigned short GetnNodes(void) override { return nNodes; }

  /*!
   * \brief Get the type of the element using VTK nomenclature.
   * \return Type of the element using VTK nomenclature.
   */
  inline unsigned short GetVTK_Type(void) override { return VTK_Type; }

  /*!
   * \brief Get the type of rotation/traslation that must be applied.
   * \return Type of the element using VTK nomenclature.
   */
  inline unsigned short GetRotation_Type(void) override { return Rotation_Type; }

  /*!
   * \brief Set the type of rotation/traslation that must be applied.
   * \param[in] val_rotation_type - Kind of rotation/traslation that must be applied.
   */
  inline void SetRotation_Type(unsigned short val_rotation_type) override { Rotation_Type = val_rotation_type; }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  void Change_Orientation(void) override;

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline unsigned short GetnNeighbor_Elements(void) override { return 0; }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline unsigned short GetnNeighbor_Nodes(unsigned short val_node) override { return 0; }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline unsigned short GetnFaces(void) override { return 0; }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline unsigned short GetnNodesFace(unsigned short val_face) override { return 0; }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline unsigned short GetMaxNodesFace(void) override { return 0; }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline unsigned short GetFaces(unsigned short val_face, unsigned short val_index) override { return 0; }

  /*!
   * \brief This function does nothing (it comes from a pure virtual function, that implies the
   *        definition of the function in all the derived classes).
   */
  inline unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index) override { return 0; }
};
