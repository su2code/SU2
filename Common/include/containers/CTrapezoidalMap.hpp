/*!
 * \file CTrapezoidalMap.hpp
 * \brief Implementation of the trapezoidal map for tabulation and lookup of fluid properties
 * \author D. Mayer, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include <string>
#include <vector>

#include "../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../Common/include/toolboxes/CSquareMatrixCM.hpp"

/*!
 * \class CTrapezoidalMap
 * \ingroup LookUpInterp
 * \brief Construction of trapezoidal map for tabulated lookup
 * \author: D. Mayer, T. Economon
 * \version 8.0.0 "Harrier"
 */
class CTrapezoidalMap {
 protected:
  /* The unique values of x which exist in the data */
  std::vector<su2double> unique_bands_x;

  su2activematrix edge_limits_x;
  su2activematrix edge_limits_y;

  su2vector<std::vector<unsigned long> > edge_to_triangle;

  /* The value that each edge which intersects the band takes within that
   * same band. Used to sort the edges */
  su2vector<std::vector<std::pair<su2double, unsigned long> > > y_edge_at_band_mid;

  double memory_footprint = 0;

 public:
  CTrapezoidalMap() = default;

  CTrapezoidalMap(const su2double* samples_x, const su2double* samples_y, const unsigned long size,
                  const std::vector<std::array<unsigned long, 2> >& edges,
                  const su2vector<std::vector<unsigned long> >& edge_to_triangle, bool display = false);

  /*!
   * \brief return the index to the triangle that contains the coordinates (val_x,val_y)
   * \param[in]  val_x  - x-coordinate or first independent variable
   * \param[in]  val_y  - y-coordinate or second independent variable
   * \param[out] val_index - index to the triangle
   */
  unsigned long GetTriangle(su2double val_x, su2double val_y);

  /*!
   * \brief get the indices of the vertical coordinate band (xmin,xmax) in the 2D search space
   * that contains the coordinate val_x
   * \param[in]  val_x  - x-coordinate or first independent variable
   * \param[out] val_band - a pair(i_low,i_up) , the lower index and upper index between which the value val_x
   * can be found
   */
  std::pair<unsigned long, unsigned long> GetBand(su2double val_x);

  /*!
   * \brief for a given coordinate (val_x,value), known to be in the band (xmin,xmax) with band index (i_low,i_up),
   * find the edges in the band (these edges come from the triangulation) that enclose the coordinate
   * \param[in]  val_band - pair i_low,i_up
   * \param[in]  val_x  - x-coordinate or first independent variable
   * \param[in]  val_y  - y-coordinate or first independent variable
   * \param[out] pair (edge_low,edge_up) - lower edge and upper edge of a triangle that encloses the coordinate
   */
  std::pair<unsigned long, unsigned long> GetEdges(std::pair<unsigned long, unsigned long> val_band, su2double val_x,
                                                   su2double val_y) const;

  /*!
   * \brief determine if the x-coordinate falls within the bounds xmin,xmax of the table
   * \param[in]  val_x  - x-coordinate or first independent variable
   * \param[out] bool - true if val_x is within (xmin,xmax)
   */
  inline bool IsInsideHullX(su2double val_x) {
    return (val_x >= unique_bands_x.front()) && (val_x <= unique_bands_x.back());
  }

  /*!
   * \brief get memory footprint of trapezoidal map.
   * \return - memory footprint in mega bytes.
   */
  double GetMemoryFootprint() const { return memory_footprint; }
};
