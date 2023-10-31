/*!
 * \file CNearestNeighbor.hpp
 * \brief Nearest Neighbor interpolation class.
 * \author H. Kline
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

#include "CInterpolator.hpp"

/*!
 * \brief Nearest Neighbor(s) interpolation.
 * \note The closest k neighbors are used for IDW interpolation, the computational
 * cost of setting up the interpolation is O(N^2 log(k)), this can be improved
 * by using an ADT.
 * \ingroup Interfaces
 */
class CNearestNeighbor final : public CInterpolator {
 private:
  su2double AvgDistance = 0.0, MaxDistance = 0.0;

  /*! \brief Helper struct to (partially) sort neighbours according to distance while
   *         keeping track of the origin of the point (i.e. index and processor). */
  struct DonorInfo {
    su2double dist;
    unsigned pidx;
    int proc;
    DonorInfo(su2double d = 0.0, unsigned i = 0, int p = 0) : dist(d), pidx(i), proc(p) {}
  };

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - index of the donor zone.
   * \param[in] jZone - index of the target zone.
   */
  CNearestNeighbor(CGeometry**** geometry_container, const CConfig* const* config, unsigned int iZone,
                   unsigned int jZone);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes.
   * \param[in] config - Definition of the particular problem.
   */
  void SetTransferCoeff(const CConfig* const* config) override;

  /*!
   * \brief Print interpolation statistics.
   */
  void PrintStatistics(void) const override;
};
