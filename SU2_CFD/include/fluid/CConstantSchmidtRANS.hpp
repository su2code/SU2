/*!
 * \file CConstantSchmidtRANS.hpp
 * \brief Defines a mass diffusivity model with constant Schmidt numbers for RANS.
 * \author T. Economon, C. Morales
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
 * \class CDiffusivityModelRANS
 * \brief Defines a mass diffusivity model for species equations based on Schmidt numbers for RANS.
 * \author T. Economon
 */
class CConstantSchmidtRANS : public CDiffusivityModel {
 public:

  /*!
   * \brief Constructor of the class.
   */
  CConstantSchmidtRANS(su2double diff_lam_const, su2double sc_turb)
      : diff_lam_const_(diff_lam_const), sc_turb_(sc_turb) {}

  /*!
   * \brief Set mass diffusivity.
   */
  void SetDiffusivity(su2double t, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp,
                      su2double kt) override {
    diff_ = diff_lam_const_ + mu_turb / (rho * sc_turb_);
  }

 private:
  su2double diff_lam_const_{0.0};
  su2double sc_turb_{0.0};
};
