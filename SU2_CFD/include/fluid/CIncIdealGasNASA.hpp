/*!
 * \file CIncIdealGasNASA.hpp
 * \brief Defines the incompressible Ideal Gas model with NASA polynomials for Cp.
 * \author Pratyksh Gupta
 * \version 8.4.0 "Harrier"
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
#include <iostream>

#include "CFluidModel.hpp"

/*!
 * \class CIncIdealGasNASA
 * \brief Child class for defining an incompressible ideal gas model with NASA polynomials.
 * \author Pratyksh Gupta
 * 
 * Implements NASA 7-coefficient polynomial format (NASA SP-273) for thermodynamic properties:
 * Ref: McBride, B.J., Zehe, M.J., and Gordon, S., "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species", NASA/TP-2002-211556, 2002.
 *   Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
 *   H/(R*T) = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
 *   S/R = a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
 * 
 * Uses a single temperature range provided via CP_POLYCOEFFS (indices 0-6).
 */
template <int N_COEFFS = 7>
class CIncIdealGasNASA final : public CFluidModel {
 public:
  /*!
   * \brief Constructor of the class.
   */
  CIncIdealGasNASA(su2double val_gas_constant, su2double val_operating_pressure, su2double val_Temperature_Ref) {
    /*--- In the incompressible ideal gas model, the thermodynamic pressure is decoupled
    from the governing equations and held constant. The density is therefore only a
    function of temperature variations. The gas is incompressible, so Cp = Cv (gamma = 1). ---*/
    Gas_Constant = val_gas_constant;
    Pressure = val_operating_pressure;
    Gamma = 1.0;
    Std_Ref_Temp_ND = val_Temperature_Ref;
  }

  /*!
   * \brief Set the NASA polynomial coefficients for variable Cp.
   * \param[in] config - configuration container for the problem.
   */
  void SetCpModel(const CConfig* config, su2double val_Temperature_Ref) override {
    
    /*--- Read NASA coefficients from the standard polynomial coefficient array (CP_POLYCOEFFS).
          Indices 0-4: Cp coefficients (a1-a5)
          Index 5: Enthalpy constant (a6)
          Index 6: Entropy constant (a7) ---*/
    for (int i = 0; i < N_COEFFS; ++i) {
      if (i < config->GetnPolyCoeffs()) {
        coeffs_[i] = config->GetCp_PolyCoeff(i);
      } else {
        coeffs_[i] = 0.0; 
      }
    }
    
    Temperature_Min = config->GetTemperatureLimits(0);
    Temperature_Max = config->GetTemperatureLimits(1);
  }

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] t - Temperature value at the point.
   */
  void SetTDState_T(su2double t, const su2double *val_scalars = nullptr) override {
    Temperature = t;
    Density = Pressure / (Temperature * Gas_Constant);

    const su2double a1 = coeffs_[0];
    const su2double a2 = coeffs_[1];
    const su2double a3 = coeffs_[2];
    const su2double a4 = coeffs_[3];
    const su2double a5 = coeffs_[4];
    const su2double a6 = coeffs_[5];

    su2double Cp_over_R = a1 + a2*t + a3*t*t + a4*t*t*t + a5*t*t*t*t;
    
    Cp = Cp_over_R * Gas_Constant;
    
    su2double H_over_RT = a1 + a2*t/2.0 + a3*t*t/3.0 + a4*t*t*t/4.0 + a5*t*t*t*t/5.0 + a6/t;
    
    Enthalpy = H_over_RT * Gas_Constant * t;
    Cv = Cp / Gamma;
  }

  /*!
   * \brief Set the Dimensionless State using enthalpy.
   * \param[in] val_enthalpy - Enthalpy value at the point.
   */
  void SetTDState_h(su2double val_enthalpy, const su2double* val_scalars = nullptr) override {
    Enthalpy = val_enthalpy;
    
    const su2double toll = 1e-5;
    su2double temp_iter = (Temperature_Min + Temperature_Max) / 2.0; /* Start in middle of allowed range */
    if (temp_iter < 1.0) temp_iter = 300.0; /* Fallback if limits are not set or zero */
    
    const su2double a1 = coeffs_[0];
    const su2double a2 = coeffs_[1];
    const su2double a3 = coeffs_[2];
    const su2double a4 = coeffs_[3];
    const su2double a5 = coeffs_[4];
    const su2double a6 = coeffs_[5];

    su2double Cp_iter = 0.0;
    su2double delta_temp_iter = 1e10;
    su2double delta_enthalpy_iter;
    const int counter_limit = 50;
    int counter = 0;
    
    while ((abs(delta_temp_iter) > toll) && (counter++ < counter_limit)) {
      
      su2double Cp_over_R = a1 + a2*temp_iter + a3*temp_iter*temp_iter + 
                            a4*temp_iter*temp_iter*temp_iter + a5*temp_iter*temp_iter*temp_iter*temp_iter;
      
      Cp_iter = Cp_over_R * Gas_Constant;
      
      su2double H_over_RT = a1 + a2*temp_iter/2.0 + a3*temp_iter*temp_iter/3.0 + 
                            a4*temp_iter*temp_iter*temp_iter/4.0 + 
                            a5*temp_iter*temp_iter*temp_iter*temp_iter/5.0 + a6/temp_iter;
      
      su2double Enthalpy_iter = H_over_RT * Gas_Constant * temp_iter;
      
      delta_enthalpy_iter = Enthalpy - Enthalpy_iter;
      delta_temp_iter = delta_enthalpy_iter / Cp_iter;
      temp_iter += delta_temp_iter;
      
      if (temp_iter < Temperature_Min) {
        // Clamp to min but don't break immediately to allow recovery if simple overshoot
         if (abs(delta_temp_iter) < toll) break; 
         temp_iter = Temperature_Min;
      }
      if (temp_iter > Temperature_Max) {
         if (abs(delta_temp_iter) < toll) break;
         temp_iter = Temperature_Max;
      }
    }
    
    Temperature = temp_iter;
    Cp = Cp_iter;
    
    if (counter == counter_limit) {
      if (SU2_MPI::GetRank() == MASTER_NODE)
          std::cout << "Warning: Newton-Raphson exceeds max. iterations in temperature computation (NASA Model)." << std::endl;
    }
    
    Density = Pressure / (Temperature * Gas_Constant);
    Cv = Cp / Gamma;
  }

 private:
  su2double Gas_Constant{0.0};     /*!< \brief Specific Gas Constant. */
  su2double Gamma{0.0};            /*!< \brief Ratio of specific heats. */
  su2double Std_Ref_Temp_ND{0.0};  /*!< \brief Nondimensional standard reference temperature for enthalpy. */
  
  std::array<su2double, N_COEFFS> coeffs_;  /*!< \brief NASA polynomial coefficients. */
  
  su2double Temperature_Min{0.0};  /*!< \brief Minimum temperature allowed. */
  su2double Temperature_Max{1e10}; /*!< \brief Maximum temperature allowed. */
};
