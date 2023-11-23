/*!
 * \file CConstantViscosity.hpp
 * \brief Defines a constant laminar viscosity model.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna, T. Economon
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

#include "CViscosityModel.hpp"

/*!
 * \class CConstantViscosity
 * \brief Defines a constant viscosity model.
 * \author S.Vitale, M.Pini
 */
class CConstantViscosity final : public CViscosityModel {
 public:
  /*!
   * \brief Constructor of the class.
   */
  CConstantViscosity(su2double mu_const) {
    mu_ = mu_const;
  }

  /*!
   * \brief Set Viscosity.
   */
  void SetViscosity(su2double t, su2double rho) override {}
};
