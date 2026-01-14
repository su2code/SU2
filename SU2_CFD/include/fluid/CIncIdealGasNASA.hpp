/*!
 * \file CIncIdealGasNASA.hpp
 * \brief Defines the incompressible Ideal Gas model with NASA polynomials for Cp.
 * \author Pratyksh Gupta (based on CIncIdealGasPolynomial by T. Economon)
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
 * Implements NASA 7-coefficient polynomial format for thermodynamic properties:
 *   Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
 *   H/(R*T) = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
 * 
 * Supports dual temperature ranges (low/high) with separate coefficients.
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
    T_mid = config->GetNASA_TempMid();
    T_low = config->GetNASA_TempLow();
    T_high = config->GetNASA_TempHigh();
    
    for (int i = 0; i < N_COEFFS; ++i) {
      coeffs_low_[i] = config->GetNASA_CoeffLowND(i);
      coeffs_high_[i] = config->GetNASA_CoeffHighND(i);
    }
    
    Temperature_Min = config->GetTemperatureLimits(0);
    Temperature_Max = config->GetTemperatureLimits(1);
    
    if (T_mid <= T_low || T_high <= T_mid) {
      SU2_MPI::Error("Invalid NASA polynomial temperature ranges. Must have T_low < T_mid < T_high.", 
                     CURRENT_FUNCTION);
    }
  }

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] t - Temperature value at the point.
   */
  void SetTDState_T(su2double t, const su2double *val_scalars = nullptr) override {
    Temperature = t;
    Density = Pressure / (Temperature * Gas_Constant);

    const su2double* coeffs = (t < T_mid) ? coeffs_low_.data() : coeffs_high_.data();
    
    su2double Cp_over_R = coeffs[0] + coeffs[1]*t + coeffs[2]*t*t + 
                          coeffs[3]*t*t*t + coeffs[4]*t*t*t*t;
    
    Cp = Cp_over_R * Gas_Constant;
    
    su2double H_over_RT = coeffs[0] + coeffs[1]*t/2.0 + coeffs[2]*t*t/3.0 + 
                          coeffs[3]*t*t*t/4.0 + coeffs[4]*t*t*t*t/5.0 + coeffs[5]/t;
    
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
    su2double temp_iter = (T_low + T_high) / 2.0;
    su2double Cp_iter = 0.0;
    su2double delta_temp_iter = 1e10;
    su2double delta_enthalpy_iter;
    const int counter_limit = 50;
    int counter = 0;
    
    while ((abs(delta_temp_iter) > toll) && (counter++ < counter_limit)) {
      const su2double* coeffs = (temp_iter < T_mid) ? coeffs_low_.data() : coeffs_high_.data();
      
      su2double Cp_over_R = coeffs[0] + coeffs[1]*temp_iter + coeffs[2]*temp_iter*temp_iter + 
                            coeffs[3]*temp_iter*temp_iter*temp_iter + 
                            coeffs[4]*temp_iter*temp_iter*temp_iter*temp_iter;
      
      Cp_iter = Cp_over_R * Gas_Constant;
      
      su2double H_over_RT = coeffs[0] + coeffs[1]*temp_iter/2.0 + coeffs[2]*temp_iter*temp_iter/3.0 + 
                            coeffs[3]*temp_iter*temp_iter*temp_iter/4.0 + 
                            coeffs[4]*temp_iter*temp_iter*temp_iter*temp_iter/5.0 + coeffs[5]/temp_iter;
      
      su2double Enthalpy_iter = H_over_RT * Gas_Constant * temp_iter;
      
      delta_enthalpy_iter = Enthalpy - Enthalpy_iter;
      delta_temp_iter = delta_enthalpy_iter / Cp_iter;
      temp_iter += delta_temp_iter;
      
      if (temp_iter < Temperature_Min) {
        cout << "Warning: Negative temperature found during Newton-Raphson." << endl;
        temp_iter = Temperature_Min;
        break;
      }
      if (temp_iter > Temperature_Max) {
        cout << "Warning: Temperature exceeds maximum during Newton-Raphson." << endl;
        temp_iter = Temperature_Max;
        break;
      }
    }
    
    Temperature = temp_iter;
    Cp = Cp_iter;
    
    if (counter == counter_limit) {
      cout << "Warning: Newton-Raphson exceeds max. iterations in temperature computation." << endl;
    }
    
    Density = Pressure / (Temperature * Gas_Constant);
    Cv = Cp / Gamma;
  }

 private:
  su2double Gas_Constant{0.0};     /*!< \brief Specific Gas Constant. */
  su2double Gamma{0.0};            /*!< \brief Ratio of specific heats. */
  su2double Std_Ref_Temp_ND{0.0};  /*!< \brief Nondimensional standard reference temperature for enthalpy. */
  
  std::array<su2double, N_COEFFS> coeffs_low_;   /*!< \brief NASA polynomial coefficients for low temperature range. */
  std::array<su2double, N_COEFFS> coeffs_high_;  /*!< \brief NASA polynomial coefficients for high temperature range. */
  
  su2double T_low{200.0};          /*!< \brief Lower temperature bound. */
  su2double T_mid{1000.0};         /*!< \brief Mid-point temperature (split between low/high ranges). */
  su2double T_high{6000.0};        /*!< \brief Upper temperature bound. */
  
  su2double Temperature_Min{0.0};  /*!< \brief Minimum temperature allowed in Newton-Raphson iterations. */
  su2double Temperature_Max{1e10}; /*!< \brief Maximum temperature allowed in Newton-Raphson iterations. */
};
