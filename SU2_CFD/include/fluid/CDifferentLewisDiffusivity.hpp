/*!
 * \file CDifferentLewisDiffusivity.hpp
 * \brief Defines different Lewis mass diffusivity.
 * \author Cristopher Morales
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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
 * \class CDifferentLewisDiffusivity
 * \brief Defines different Lewis mass diffusivity model for species equations.
 * \author Cristopher Morales
 */
class CDifferentLewisDiffusivity final : public CDiffusivityModel {
 public:
  /*!
   * \brief Constructor of the class.
   */
  CDifferentLewisDiffusivity(su2double Lewis) : Lewis_(Lewis) {}

  /*!
   * \brief Set diffusivity.
   */
  void SetDiffusivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp,
                      su2double kt) override {
    diff_ = kt / (Lewis_ * rho * cp);
  }

 private:
  const su2double Lewis_{0.0};
};
