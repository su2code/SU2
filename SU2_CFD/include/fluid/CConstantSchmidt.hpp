/*!
 * \file CConstantSchmidt.hpp
 * \brief Defines a mass diffusivity model with constant Schmidt numbers.
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

#include "CDiffusivityModel.hpp"

/*!
 * \class CDiffusivityModel
 * \brief Defines a mass diffusivity model for species equations based on Schmidt number.
 * \author T. Economon
 */
class CConstantSchmidt final : public CDiffusivityModel {
 public:

  /*!
   * \brief Constructor of the class.
   */
  CConstantSchmidt(su2double sc_lam) : sc_lam_(sc_lam) {}

  /*!
   * \brief Set diffusivity.
   */
  void SetDiffusivity(su2double rho, su2double mu_lam, su2double cp, su2double kt) override {
    diff_ = mu_lam / (rho * sc_lam_);
  }

 private:
  su2double sc_lam_{0.0};
};
