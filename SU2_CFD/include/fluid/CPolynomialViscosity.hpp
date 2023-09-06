/*!
 * \file CPolynomialViscosity.hpp
 * \brief Defines a laminar viscosity model as a polynomial function of temperature.
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

#include <array>

#include "CViscosityModel.hpp"

/*!
 * \class CPolynomialViscosity
 * \brief Defines laminar viscosity as a polynomial function of temperature.
 * \author T. Economon
 */
template <int N>
class CPolynomialViscosity final : public CViscosityModel {
 public:
  /*!
   * \brief Constructor of the class.
   */
  CPolynomialViscosity(const su2double* val_coeffs) {
    for (int i = 0; i < N; ++i) {
      coeffs_[i] = val_coeffs[i];
    }
  }

  /*!
   * \brief Set Viscosity.
   */
  void SetViscosity(su2double t, su2double rho) override {
    /* Evaluate the new mu from the coefficients and temperature. */
    mu_ = coeffs_[0];
    su2double t_i = 1.0;
    for (int i = 1; i < N; ++i) {
      t_i *= t;
      mu_ += coeffs_[i] * t_i;
    }
  }

 private:
  array<su2double, N> coeffs_; /*!< \brief Polynomial coefficients for viscosity as a function of temperature. */
};
