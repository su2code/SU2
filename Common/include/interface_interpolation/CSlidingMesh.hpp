/*!
 * \file CSlidingMesh.hpp
 * \brief Sliding mesh interpolation.
 * \author H. Kline
 * \version 7.0.2 "Blackbird"
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

#include "CInterpolator.hpp"

/*!
 * \brief Sliding mesh approach.
 */
class CSlidingMesh final : public CInterpolator {
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - index of the donor zone
   * \param[in] jZone - index of the target zone
   */
  CSlidingMesh(CGeometry ****geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] config - Definition of the particular problem.
   */
  void Set_TransferCoeff(CConfig **config) override;

private:
  /*!
   * \brief For 3-Dimensional grids, build the dual surface element
   * \param[in] map         - array containing the index of the boundary points connected to the node
   * \param[in] startIndex  - for each vertex specifies the corresponding index in the global array containing the indexes of all its neighbouring vertexes
   * \param[in] nNeighbour  - for each vertex specifies the number of its neighbouring vertexes (on the boundary)
   * \param[in] coord       - array containing the coordinates of all the boundary vertexes
   * \param[in] centralNode - label of the vertex around which the dual surface element is built
   * \param[in] element  - double array where element node coordinates will be stored
   */
  int Build_3D_surface_element(unsigned long *map, unsigned long *startIndex, unsigned long* nNeighbor,
                               su2double *coord, unsigned long centralNode, su2double** element);

  /*!
   * \brief For 2-Dimensional grids, compute intersection length of two segments projected along a given direction
   * \param[in] A1 - first  point of segment A
   * \param[in] A2 - second point of segment A
   * \param[in] B1 - first  point of segment B
   * \param[in] B2 - second point of segment B
   * \param[in] Direction - along which segments are projected
   */
  su2double ComputeLineIntersectionLength(su2double* A1, su2double* A2, su2double* B1, su2double* B2, su2double* Direction);

  /*!
   * \brief For 3-Dimensional grids, compute intersection area between two triangle projected on a given plane
   * \param[in] A1 - first  point of triangle A
   * \param[in] A2 - second point of triangle A
   * \param[in] A3 - third  point of triangle A
   * \param[in] B1 - first  point of triangle B
   * \param[in] B2 - second point of triangle B
   * \param[in] B3 - third  point of triangle B
   * \param[in] Direction - vector normal to projection plane
   */
  su2double Compute_Triangle_Intersection(su2double* A1, su2double* A2, su2double* A3, su2double* B1, su2double* B2, su2double* B3, su2double* Direction);

  /*!
   * \brief For 3-Dimensional grids, compute intersection area between two triangle projected on a given plane
   * P1 from triangle P MUST be inside triangle Q, points order doesn't matter
   * \param[in] P1 - first  point of triangle A
   * \param[in] P2 - second point of triangle A
   * \param[in] P3 - third  point of triangle A
   * \param[in] Q1 - first  point of triangle B
   * \param[in] Q2 - second point of triangle B
   * \param[in] Q3 - third  point of triangle B
   */
  su2double ComputeIntersectionArea( su2double* P1, su2double* P2, su2double* P3, su2double* Q1, su2double* Q2, su2double* Q3 );

  /*!
   * \brief For 2-Dimensional grids, check whether, and compute, two lines are intersecting
   * \param[in] A1 - first  defining first line
   * \param[in] A2 - second defining first line
   * \param[in] B1 - first  defining second line
   * \param[in] B2 - second defining second line
   * \param[in] IntersectionPoint - Container for intersection coordinates
   */
  void ComputeLineIntersectionPoint( su2double* A1, su2double* A2, su2double* B1, su2double* B2, su2double* IntersectionPoint );

  /*!
   * \brief For N-Dimensional grids, check whether a point is inside a triangle specified by 3 T points
   * \param[in] Point - query point
   * \param[in] T1 - first  point of triangle T
   * \param[in] T2 - second point of triangle T
   * \param[in] T3 - third  point of triangle T
   */
  bool CheckPointInsideTriangle(su2double* Point, su2double* T1, su2double* T2, su2double* T3);

};
