/*!
 * \file CPolynomialConductivityRANS.hpp
 * \brief Defines a non-constant thermal conductivity using a polynomial function of temperature
 *        for RANS problems with the addition of a turbulent component based on a turbulent Prandtl number.
 * \author T. Economon
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
 * \class CPolynomialConductivityRANS
 * \brief Defines a non-constant thermal conductivity using a polynomial function of temperature
 *        for RANS problems with the addition of a turbulent component based on a turbulent Prandtl number.
 *  \author T. Economon
 */
template <int N>
class CPolynomialConductivityRANS final : public CConductivityModel {
 public:
  /*!
   * \brief Constructor of the class.
   */
  CPolynomialConductivityRANS(const su2double* val_coeffs, su2double pr_turb) : pr_turb_(pr_turb) {
    for (int i = 0; i < N; ++i) {
      coeffs_[i] = val_coeffs[i];
    }
  }

  /*!
   * \brief Set thermal conductivity.
   */
  void SetConductivity(su2double t, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp,
                       su2double, su2double) override {
    /* Evaluate the new kt from the coefficients and temperature. */
    kt_ = coeffs_[0];
    su2double t_i = 1.0;
    for (int i = 1; i < N; ++i) {
      t_i *= t;
      kt_ += coeffs_[i] * t_i;
    }

    /* Add the component due to turbulence to compute effective conductivity. */
    kt_ += cp * mu_turb / pr_turb_;
  }

 private:
  array<su2double, N> coeffs_;    /*!< \brief Polynomial coefficients for conductivity as a function of temperature. */
  const su2double pr_turb_{0.0};  /*!< \brief Turbulent Prandtl number. */
};
