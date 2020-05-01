/*!
 * \file CEdge.hpp
 * \brief Declaration of the edge class <i>CEdge.cpp</i> file.
 * \author F. Palacios, T. Economon
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

#include "../../toolboxes/C2DContainer.hpp"

/*!
 * \class CEdge
 * \brief Class for defining the edges of the dual grid.
 * \author F. Palacios
 */
class CEdge {
  static_assert(su2activematrix::Storage == StorageType::RowMajor, "Needed to return normal as pointer.");

private:
  su2matrix<unsigned long> Nodes; /*!< \brief Vector to store the node indices of the edge. */
  su2activematrix Normal;         /*!< \brief Normal (area) of the edge. */
  su2activematrix Coord_CG;       /*!< \brief Center-of-gravity (mid point) of the edge. */

public:
  enum NodePosition : unsigned long {LEFT = 0, RIGHT = 1};

  /*!
   * \brief Constructor of the class.
   * \param[in] nEdge - Number of edges
   * \param[in] nDim - Number of dimensions of the problem.
   */
  CEdge(unsigned long nEdge, unsigned long nDim);

  /*!
   * \brief No default construction.
   */
  CEdge() = delete;

  /*!
   * \brief Set the center of gravity of the edge.
   * \param[in] iEdge - Edge index.
   * \param[in] nodeCoord - Coordinates of the two nodes.
   */
  template<class T>
  void SetCoord_CG(unsigned long iEdge, const T& nodeCoord) {
    for (auto iDim = 0u; iDim < Coord_CG.cols(); ++iDim)
      Coord_CG(iEdge,iDim) = 0.5 * (nodeCoord[0][iDim] + nodeCoord[1][iDim]);
  }

  /*!
   * \brief Obtain the center of gravity of the edge.
   * \param[in] iEdge - Edge index.
   * \param[in] iDim - Dimension.
   * \return Coordinate of the centre of gravity.
   */
  inline su2double GetCG(unsigned long iEdge, unsigned long iDim) const { return Coord_CG(iEdge,iDim); }

  /*!
   * \brief Get left/right node index defining the edge.
   * \param[in] iEdge - Edge index.
   * \param[in] iNode - Node index 0 or 1, LEFT or RIGHT.
   * \return Index of the node that composes the edge.
   */
  inline unsigned long GetNode(unsigned long iEdge, unsigned long iNode) const { return Nodes(iEdge,iNode); }

  /*!
   * \brief Set the node indices of an edge.
   * \param[in] iEdge - Edge index.
   * \param[in] iPoint - Index of left node.
   * \param[in] jPoint - Index of right node.
   */
  inline void SetNodes(unsigned long iEdge, unsigned long iPoint, unsigned long jPoint) {
    Nodes(iEdge, LEFT) = iPoint;
    Nodes(iEdge, RIGHT) = jPoint;
  }

  /*!
   * \brief Get the number of nodes of an edge (2).
   */
  inline unsigned long GetnNodes() const { return 2; }

  /*!
   * \brief Compute the volume associated with an edge (3D version).
   * \param[in] coord_Edge_CG - Coordinates of the centre of gravity of the edge.
   * \param[in] coord_FaceElem_CG - Coordinates of the centre of gravity of the face of an element.
   * \param[in] coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] coord_Point - Coordinates of the point that form the control volume.
   * \return Local volume associated to the edge.
   */
  static su2double GetVolume(const su2double* coord_Edge_CG,
                             const su2double* coord_FaceElem_CG,
                             const su2double* coord_Elem_CG,
                             const su2double* coord_Point);

  /*!
   * \brief Compute the volume associated with an edge (2D version).
   * \param[in] coord_Edge_CG - Coordinates of the centre of gravity of the edge.
   * \param[in] coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] coord_Point - Coordinates of the point that form the control volume.
   * \return Local volume associated to the edge.
   */
  static su2double GetVolume(const su2double* coord_Edge_CG,
                             const su2double* coord_Elem_CG,
                             const su2double* coord_Point);

  /*!
   * \brief Set the face that corresponds to an edge (3D version).
   * \param[in] iEdge - Edge index.
   * \param[in] coord_Edge_CG - Coordinates of the centre of gravity of the edge.
   * \param[in] coord_FaceElem_CG - Coordinates of the centre of gravity of the face of an element.
   * \param[in] coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
   * \return Compute the normal (dimensional) to the face that makes the control volume boundaries.
   */
  void SetNodes_Coord(unsigned long iEdge,
                      const su2double* coord_Edge_CG,
                      const su2double* coord_FaceElem_CG,
                      const su2double* coord_Elem_CG);

  /*!
   * \brief Set the face that corresponds to an edge (2D version).
   * \param[in] iEdge - Edge index.
   * \param[in] coord_Edge_CG - Coordinates of the centre of gravity of the edge.
   * \param[in] coord_Elem_CG - Coordinates of the centre of gravity of the element.
   * \param[in] config - Definition of the particular problem.
   * \return Compute the normal (dimensional) to the face that makes the contorl volume boundaries.
   */
  void SetNodes_Coord(unsigned long iEdge,
                      const su2double* coord_Edge_CG,
                      const su2double* coord_Elem_CG);

  /*!
   * \brief Copy the the normal vector of a face.
   * \param[in] iEdge - Edge index.
   * \param[out] normal - Object into which the normal (dimensional) will be copied.
   */
  template<class T>
  inline void GetNormal(unsigned long iEdge, T& normal) const {
    for (auto iDim = 0ul; iDim < Normal.cols(); iDim++)
      normal[iDim] = Normal(iEdge,iDim);
  }

  /*!
   * \brief Get the normal to a face of the control volume asociated with an edge.
   * \param[in] iEdge - Edge index.
   * \return Dimensional normal vector, the modulus is the area of the face.
   */
  inline const su2double* GetNormal(unsigned long iEdge) const { return Normal[iEdge]; }

  /*!
   * \brief Initialize normal vector to 0.
   */
  void SetZeroValues(void);

  /*!
   * \brief Set the normal vector of an edge.
   * \param[in] iEdge - Edge index.
   * \param[in] normal - Vector to initialize the normal vector.
   * \return Value of the normal vector.
   */
  template<class T>
  void SetNormal(unsigned long iEdge, const T& normal) {
    for (auto iDim = 0ul; iDim < Normal.cols(); ++iDim)
      Normal(iEdge,iDim) = normal[iDim];
  }

  /*!
   * \brief Add a vector to the normal vector of an edge.
   * \param[in] iEdge - Edge index.
   * \param[in] normal - Vector to add to the normal vector.
   */
  template<class T>
  void AddNormal(unsigned long iEdge, const T& normal) {
    for (auto iDim = 0ul; iDim < Normal.cols(); ++iDim)
      Normal(iEdge,iDim) += normal[iDim];
  }

  /*!
   * \brief Subtract a vector to the normal vector of an edge.
   * \param[in] iEdge - Edge index.
   * \param[in] normal - Vector to add to the normal vector.
   */
  template<class T>
  void SubNormal(unsigned long iEdge, const T& normal) {
    for (auto iDim = 0ul; iDim < Normal.cols(); ++iDim)
      Normal(iEdge,iDim) -= normal[iDim];
  }

};
