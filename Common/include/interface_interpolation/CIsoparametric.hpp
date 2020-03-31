/*!
 * \file CIsoparametric.hpp
 * \brief Isoparametric interpolation using FE shape functions.
 * \author H. Kline
 * \version 7.0.3 "Blackbird"
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
 * \brief Isoparametric interpolation.
 */
class CIsoparametric final : public CInterpolator {
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - index of the donor zone
   * \param[in] jZone - index of the target zone
   */
  CIsoparametric(CGeometry ****geometry_container, const CConfig* const* config, unsigned int iZone, unsigned int jZone);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] config - Definition of the particular problem.
   */
  void Set_TransferCoeff(const CConfig* const* config) override;

private:
  /*!
   * \brief Calculate the isoparametric representation of point iVertex in marker iZone_0 by
   *        nodes of element donor_elem in marker jMarker of zone iZone_1.
   * \param[in] iVertex - vertex index of the point being interpolated.
   * \param[in] nDim - the dimension of the coordinates.
   * \param[in] iZone_1 - zone index of the element to use for interpolation (the DONOR zone)
   * \param[in] donor_elem - element index of the element to use for interpolation (or global index of a point in 2D)
   * \param[in] nDonorPoints - number of donor points in the element.
   * \param[in] xj - point projected onto the plane of the donor element.
   * \param[out] isoparams - isoparametric coefficients. Must be allocated to size nNodes ahead of time. (size> nDonors)
   *
   * \note If the problem is 2D, the 'face' projected onto is actually an edge; the local index
   * of the edge is then stored in iFace, and the global index of the node (from which the edge
   * is referenced)
   */
  void Isoparameters(unsigned short nDim, unsigned short nDonor, const su2double *X,
                     const su2double *xj, su2double* isoparams) const;

};
