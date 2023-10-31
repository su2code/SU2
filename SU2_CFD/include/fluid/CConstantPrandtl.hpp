/*!
 * \file CConstantPrandtl.hpp
 * \brief Defines a non-constant laminar Prandtl number thermal conductivity model.
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

#include "CConductivityModel.hpp"

/*!
 * \class CConstantPrandtl
 * \brief Defines a non-constant laminar thermal conductivity using a constant Prandtl's number
 * \author S.Vitale, M.Pini
 */
class CConstantPrandtl final : public CConductivityModel {
 public:
  /*!
   * \brief Constructor of the class.
   */
  CConstantPrandtl(su2double pr_lam) : pr_lam_(pr_lam) {}

  /*!
   * \brief Set thermal conductivity.
   */
  void SetConductivity(su2double t, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp,
                       su2double dmudrho_t, su2double dmudt_rho) override {
    kt_ = mu_lam * cp / pr_lam_;
    dktdrho_t_ = dmudrho_t * cp / pr_lam_;
    dktdt_rho_ = dmudt_rho * cp / pr_lam_;
  }

 private:
  const su2double pr_lam_{0.0};    /*!< \brief Laminar Prandtl number. */
};
