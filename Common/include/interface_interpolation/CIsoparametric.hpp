/*!
 * \file CIsoparametric.hpp
 * \brief Isoparametric interpolation using FE shape functions.
 * \author P. Gomes
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
 * \brief Isoparametric interpolation.
 * \ingroup Interfaces
 */
class CIsoparametric final : public CInterpolator {
 private:
  /*--- Statistics. ---*/
  su2double MaxDistance = 0.0, ErrorRate = 0.0;
  unsigned long ErrorCounter = 0;

  /*! \brief Helper struct to store information about candidate donor elements. */
  struct DonorInfo {
    su2double isoparams[4] = {0.0}; /*!< \brief Interpolation coefficients. */
    su2double distance = 0.0;       /*!< \brief Distance from target to final mapped point on donor plane. */
    unsigned iElem = 0;             /*!< \brief Identification of the element. */
    int error = 0;                  /*!< \brief If the mapped point is "outside" of the donor. */

    /*--- Best donor is one for which the mapped point is closest to target. ---*/
    bool operator<(const DonorInfo& other) const { return distance < other.distance; }
  };

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - index of the donor zone.
   * \param[in] jZone - index of the target zone.
   */
  CIsoparametric(CGeometry**** geometry_container, const CConfig* const* config, unsigned int iZone,
                 unsigned int jZone);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] config - Definition of the particular problem.
   */
  void SetTransferCoeff(const CConfig* const* config) override;

  /*!
   * \brief Print information about the interpolation.
   */
  void PrintStatistics(void) const override;

 private:
  /*!
   * \brief Compute the isoparametric interpolation coefficients for a 2D line element.
   * \param[in] X - Coordinate matrix defining the line.
   * \param[in] xj - Coordinates of the target point.
   * \param[out] isoparams - Isoparametric coefficients.
   * \return 0 on success, 1 if xj is too far outside element bounds.
   */
  static int LineIsoparameters(const su2double X[][3], const su2double* xj, su2double* isoparams);

  /*!
   * \brief Compute the isoparametric interpolation coefficients for a 3D triangle element.
   * \param[in] X - Coordinate matrix defining the triangle.
   * \param[in] xj - Coordinates of the target point.
   * \param[out] isoparams - Isoparametric coefficients.
   * \return 0 on success, 1 if xj is too far outside element bounds.
   */
  static int TriangleIsoparameters(const su2double X[][3], const su2double* xj, su2double* isoparams);

  /*!
   * \brief Compute the isoparametric interpolation coefficients for a 3D quadrilateral element.
   * \param[in] X - Coordinate matrix defining the quadrilateral.
   * \param[in] xj - Coordinates of the target point.
   * \param[out] isoparams - Isoparametric coefficients.
   * \return 0 on success, 1 if xj is too far outside element bounds.
   */
  static int QuadrilateralIsoparameters(const su2double X[][3], const su2double* xj, su2double* isoparams);
};
