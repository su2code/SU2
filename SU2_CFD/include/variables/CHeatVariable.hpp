/*!
 * \file CHeatVariable.hpp
 * \brief Class for defining the variables of the finite-volume heat equation solver.
 * \author F. Palacios, T. Economon
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

#include "CScalarVariable.hpp"

/*!
 * \class CHeatVariable
 * \brief Class for defining the variables of the finite-volume heat equation solver.
 * \author O. Burghardt
 * \version 8.0.0 "Harrier"
 */
class CHeatVariable final : public CScalarVariable {
public:
  static constexpr size_t MAXNVAR = 1; /*!< \brief Max number of variables, for static arrays. */

  /*!
   * \brief Constructor of the class.
   * \param[in] heat - Values of the Heat solution (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CHeatVariable(su2double heat, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Get the temperature of the point.
   * \return Value of the temperature of the point.
   */
  inline su2double GetTemperature(unsigned long iPoint) const final { return Solution(iPoint, 0); }

};
