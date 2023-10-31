/*!
 * \file CSobolevSmoothingVariable.hpp
 * \brief Class for defining the variables of the gradient smoothing.
 * \author T.Dick
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

#include "../../../Common/include/containers/CVertexMap.hpp"
#include "CVariable.hpp"

class CSobolevSmoothingVariable final : public CVariable {
 public:
  MatrixType Sensitivity;  /*!< Vector holding the derivative of target functional with respect to the coordinates at this node */

  CVertexMap<unsigned> BoundaryVertexMap; /*!< \brief Stores if a point belongs to the boundary of a boundary. */

  /*!
   * \brief Constructor of the class.
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSobolevSmoothingVariable(unsigned long npoint, unsigned long ndim, const CConfig *config);

  /*!
   * \brief Set the sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline void SetSensitivity(unsigned long iPoint, unsigned long iDim, su2double val) override { Sensitivity(iPoint,iDim) = val;}

  /*!
   * \brief Get the Sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline su2double GetSensitivity(unsigned long iPoint, unsigned long iDim) const override { return Sensitivity(iPoint,iDim); }

  /*!
   * \brief Mark a point as boundary of a boundary
   */
  void MarkAsBoundaryPoint(unsigned long iPoint);

  /*!
   * \brief return wether a point is a boundary of a boundary
   */
  bool GetIsBoundaryPoint(unsigned long iPoint) const;

  /*!
   * \brief Allocate member variables for points marked as vertex (via "MarkAsBoundaryPoint").
   */
  void AllocateBoundaryVariables();
};
