/*!
 * \file CConstantLewisDiffusivity.hpp
 * \brief Defines Constant Lewis mass diffusivity.
 * \author M.Heimgartner, C.Morales
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
 * \class CConstantLewisDiffusivity
 * \brief Defines a Constant Lewis mass diffusivity model for species equations.
 * \author M.Heimgartner, C.Morales
 */
class CConstantLewisDiffusivity final : public CDiffusivityModel {
 public:
  /*!
   * \brief Constructor of the class.
   */
  CConstantLewisDiffusivity(su2double Lewis) : Lewis_(Lewis) {}

  /*!
   * \brief Set diffusivity.
   */
  void SetDiffusivity(su2double rho, su2double mu_lam, su2double cp, su2double kt) override {
    diff_ = kt / (Lewis_ * rho * cp);
  }

 private:
  const su2double Lewis_{1.0};
};
