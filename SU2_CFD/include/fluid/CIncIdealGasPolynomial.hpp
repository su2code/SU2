/*!
 * \file CIncIdealGasPolynomial.hpp
 * \brief Defines the incompressible Ideal Gas model with polynomial Cp.
 * \author T. Economon
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "CFluidModel.hpp"

/*!
 * \class CIncIdealGasPolynomial
 * \brief Child class for defining a custom incompressible ideal gas model.
 * \author: T. Economon
 */
template <int N>
class CIncIdealGasPolynomial final : public CFluidModel {
 public:
  /*!
   * \brief Constructor of the class.
   */
  CIncIdealGasPolynomial(su2double val_gas_constant, su2double val_operating_pressure) {
    /* In the incompressible ideal gas model, the thermodynamic pressure is decoupled
    from the governing equations and held constant. The density is therefore only a
    function of temperature variations. We also use a molecular weight (g/mol) and the
    universal gas constant to compute the specific gas constant for the fluid. The gas
    is incompressible, so Cp = Cv (gamma = 1). */
    Gas_Constant = val_gas_constant;
    Pressure = val_operating_pressure;
    Gamma = 1.0;
  }

  /*!
   * \brief Set the temperature polynomial coefficients for variable Cp.
   * \param[in] config - configuration container for the problem.
   */
  void SetCpModel(const CConfig* config) override {
    Enthalpy_Ref = 0.0;
    su2double t_i = 1.0;
    for (int i = 0; i < N; ++i) {
      t_i *= STD_REF_TEMP / config->GetInc_Temperature_Ref();
      coeffs_[i] = config->GetCp_PolyCoeffND(i);
      Enthalpy_Ref += coeffs_[i] * t_i / (i + 1);
    }
    Temperature_Min = config->GetTemperatureLimits(0);
  }

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] t - Temperature value at the point.
   */
  void SetTDState_T(su2double t, const su2double *val_scalars = nullptr) override {
    /* The EoS only depends upon temperature. */
    Temperature = t;
    Density = Pressure / (Temperature * Gas_Constant);

    /* Evaluate the new Cp and enthalpy from the coefficients and temperature. */
    Cp = coeffs_[0];
    Enthalpy = coeffs_[0] * t - Enthalpy_Ref;
    su2double t_i = 1.0;
    for (int i = 1; i < N; ++i) {
      t_i *= t;
      Cp += coeffs_[i] * t_i;
      Enthalpy += coeffs_[i] * t_i * t / (i + 1);
    }
    Cv = Cp / Gamma;
  }

  /*!
   * \brief Set the Dimensionless State using enthalpy.
   * \param[in] val_enthalpy - Enthalpy value at the point.
   */
  void SetTDState_h(su2double val_enthalpy, const su2double* val_scalars = nullptr) override {
    Enthalpy = val_enthalpy;
    /*--- convergence criterion for temperature in [K], high accuracy needed for restarts. ---*/
    const su2double toll = 1e-5;
    su2double temp_iter = 300.0;
    su2double Cp_iter = 0.0;
    su2double delta_temp_iter = 1e10;
    su2double delta_enthalpy_iter;
    const int counter_limit = 20;

    int counter = 0;

    /*--- Computing temperature given enthalpy using Newton-Raphson. ---*/
    while ((abs(delta_temp_iter) > toll) && (counter++ < counter_limit)) {
      /* Evaluate the new Cp and enthalpy from the coefficients and temperature. */
      Cp_iter = coeffs_[0];
      su2double Enthalpy_iter = coeffs_[0] * temp_iter - Enthalpy_Ref;
      su2double t_i = 1.0;
      for (int i = 1; i < N; ++i) {
        t_i *= temp_iter;
        Cp_iter += coeffs_[i] * t_i;
        Enthalpy_iter += coeffs_[i] * t_i * temp_iter / (i + 1);
      }

      delta_enthalpy_iter = Enthalpy - Enthalpy_iter;

      delta_temp_iter = delta_enthalpy_iter / Cp_iter;

      temp_iter += delta_temp_iter;
      if (temp_iter < Temperature_Min) {
        cout << "Warning: Negative temperature has been found during Newton-Raphson" << endl;
        temp_iter = Temperature_Min;
        break;
      }
    }
    Temperature = temp_iter;
    Cp = Cp_iter;
    if (counter == counter_limit) {
      cout << "Warning Newton-Raphson exceed number of max iteration in temperature computation" << endl;
    }
    Density = Pressure / (Temperature * Gas_Constant);
    Cv = Cp / Gamma;
  }

 private:
  su2double Gas_Constant{0.0}; /*!< \brief Specific Gas Constant. */
  su2double Gamma{0.0};        /*!< \brief Ratio of specific heats. */
  array<su2double, N> coeffs_; /*!< \brief Polynomial coefficients for heat capacity as a function of temperature. */
  su2double Enthalpy_Ref;      /*!< \brief Enthalpy computed at the reference temperature. */
  su2double Temperature_Min;   /*!< \brief Minimum temperature value allowed in Newton-Raphson iterations. */
};
