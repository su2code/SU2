/*!
 * \file CGridMovement.hpp
 * \brief Headers of the CGridMovement class
 * \author F. Palacios
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

#include "../geometry/CGeometry.hpp"
#include "../CConfig.hpp"

/*!
 * \class CGridMovement
 * \brief Class for moving the surface and volumetric
 *        numerical grid (2D and 3D problems).
 * \author F. Palacios
 */
class CGridMovement {
 protected:
  int rank, /*!< \brief MPI Rank. */
      size; /*!< \brief MPI Size. */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CGridMovement(void);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CGridMovement(void);

  /*!
   * \brief Set the surface/boundary deformation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \return Total deformation applied, which may be less than target if intersection prevention is used.
   */
  inline virtual vector<vector<su2double> > SetSurface_Deformation(CGeometry* geometry, CConfig* config) {
    return vector<vector<su2double> >();
  }
};
