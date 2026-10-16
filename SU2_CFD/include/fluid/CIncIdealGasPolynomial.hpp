/*!
 * \file CIncIdealGasPolynomial.hpp
 * \brief Defines the incompressible Ideal Gas model with polynomial Cp.
 * \author T. Economon
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
  CIncIdealGasPolynomial(su2double val_gas_constant, su2double val_operating_pressure, su2double val_Temperature_Ref) {
    /* In the incompressible ideal gas model, the thermodynamic pressure is decoupled
    from the governing equations and held constant. The density is therefore only a
    function of temperature variations. We also use a molecular weight (g/mol) and the
    universal gas constant to compute the specific gas constant for the fluid. The gas
    is incompressible, so Cp = Cv (gamma = 1). */
    Gas_Constant = val_gas_constant;
    Pressure = val_operating_pressure;
    Gamma = 1.0;
    Std_Ref_Temp_ND = val_Temperature_Ref;
  }

  /*!
   * \brief Set the temperature polynomial coefficients for variable Cp.
   * \param[in] config - configuration container for the problem.
   */
  void SetCpModel(const CConfig* config, su2double val_Temperature_Ref) override {
    for (int i = 0; i < N; ++i) {
      /*--- Note that we use cp = dh/dt here. ---*/
      coeffs_[i] = config->GetCp_PolyCoeffND(i);
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
    Enthalpy = coeffs_[0] * (t - Std_Ref_Temp_ND);
    su2double t_i = 1.0;
    su2double tref_i = 1.0;
    for (int i = 1; i < N; ++i) {
      t_i *= t;
      tref_i *= Std_Ref_Temp_ND;
      Cp += coeffs_[i] * t_i;
      Enthalpy += coeffs_[i] * (t_i * t - tref_i * Std_Ref_Temp_ND) / (i + 1);
    }
    Cv = Cp / Gamma;
  }

  /*!
   * \brief Set the Dimensionless State using enthalpy, starting from the reference temperature.
   * \note A failed iteration is not reported, use SolveTDState_h where that matters.
   * \param[in] val_enthalpy - Enthalpy value at the point.
   */
  void SetTDState_h(su2double val_enthalpy, const su2double* val_scalars = nullptr) override {
    SolveTDState_h(val_enthalpy, 0.0, val_scalars);
  }

  /*!
   * \brief Set the Dimensionless State using enthalpy, by Newton-Raphson iteration on temperature.
   * \param[in] val_enthalpy - Enthalpy value at the point.
   * \param[in] val_temperature_guess - Start of the iteration, <= 0 to start from the reference temperature.
   * \param[in] val_scalars - Scalar mass fractions (not used).
   * \return True if the iteration converged above the minimum temperature.
   */
  bool SolveTDState_h(su2double val_enthalpy, su2double val_temperature_guess,
                      const su2double* val_scalars = nullptr) override {
    /*--- Convergence criterion for temperature, high accuracy needed for restarts. ---*/
    constexpr su2double tolerance = 1e-5;
    constexpr int maxIterations = 20;

    Enthalpy = val_enthalpy;
    /*--- The guess is passive, the converged temperature does not depend on it. ---*/
    su2double temperature = (val_temperature_guess > 0.0) ? SU2_TYPE::GetValue(val_temperature_guess)
                                                          : Std_Ref_Temp_ND;
    su2double Cp_iter = coeffs_[0];
    bool converged = false;

    /*--- Compute temperature given enthalpy using Newton-Raphson. ---*/
    for (int iter = 0; iter < maxIterations; ++iter) {
      /*--- Evaluate Cp and the enthalpy from the coefficients at the current temperature. ---*/
      Cp_iter = coeffs_[0];
      su2double Enthalpy_iter = coeffs_[0] * (temperature - Std_Ref_Temp_ND);
      su2double t_i = 1.0;
      su2double tref_i = 1.0;
      for (int i = 1; i < N; ++i) {
        t_i *= temperature;
        tref_i *= Std_Ref_Temp_ND;
        Cp_iter += coeffs_[i] * t_i;
        Enthalpy_iter += coeffs_[i] * (t_i * temperature - tref_i * Std_Ref_Temp_ND) / (i + 1);
      }

      const su2double delta = (Enthalpy - Enthalpy_iter) / Cp_iter;

      /*--- An overshoot below the minimum temperature is clamped and the iteration continues. With a positive Cp
       the enthalpy increases with temperature, so a step below the minimum taken from the minimum itself means
       the solution is below it, stop and report failure. ---*/
      if (temperature + delta < Temperature_Min) {
        if (temperature == Temperature_Min) break;
        temperature = Temperature_Min;
        continue;
      }
      temperature += delta;
      if (fabs(delta) <= tolerance) {
        converged = true;
        break;
      }
    }

    Temperature = temperature;
    Cp = Cp_iter;
    Density = Pressure / (Temperature * Gas_Constant);
    Cv = Cp / Gamma;
    return converged;
  }

 private:
  su2double Gas_Constant{0.0}; /*!< \brief Specific Gas Constant. */
  su2double Gamma{0.0};        /*!< \brief Ratio of specific heats. */
  su2double Std_Ref_Temp_ND{0.0}; /*!< \brief Nondimensional standard reference temperature for enthalpy. */
  array<su2double, N> coeffs_; /*!< \brief Polynomial coefficients for heat capacity as a function of temperature. */
  su2double Temperature_Min;   /*!< \brief Minimum temperature value allowed in Newton-Raphson iterations. */
};
