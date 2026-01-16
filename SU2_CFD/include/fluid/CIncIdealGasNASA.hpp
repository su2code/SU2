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
 * Implements NASA 9-coefficient polynomial format for thermodynamic properties with backward compatibility for NASA-7.
 * Ref: McBride, B.J., Zehe, M.J., and Gordon, S., "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species", NASA/TP-2002-211556, 2002.
 * 
 * NASA-9 format (full):
 *   Cp/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
 *   H/(RT) = -a1*T^-2 + a2*ln(T)/T + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + a7*T^4/5 + a8/T
 *   S/R = -a1*T^-2/2 - a2*T^-1 + a3*ln(T) + a4*T + a5*T^2/2 + a6*T^3/3 + a7*T^4/4 + a9
 * 
 * NASA-7 format (subset): Set a1=a2=0 to recover the traditional 7-coefficient format.
 * 
 * Uses a single temperature range provided via CP_POLYCOEFFS (indices 0-8).
 */
template <int N_COEFFS = 9>
class CIncIdealGasNASA final : public CFluidModel {
 public:
  /*!
   * \brief Constructor of the class.
   */
  CIncIdealGasNASA(su2double val_gas_constant, su2double val_operating_pressure, su2double val_Temperature_Ref, su2double val_Ref_Temp_Dim = 1.0) {
    /*--- In the incompressible ideal gas model, the thermodynamic pressure is decoupled
    from the governing equations and held constant. The density is therefore only a
    function of temperature variations. The gas is incompressible, so Cp = Cv (gamma = 1). ---*/
    Gas_Constant = val_gas_constant;
    Pressure = val_operating_pressure;
    Gamma = 1.0;
    Std_Ref_Temp_ND = val_Temperature_Ref;
    Ref_Temp_Dim = val_Ref_Temp_Dim;
  }

  /*!
   * \brief Set the NASA polynomial coefficients for variable Cp.
   * \param[in] config - configuration container for the problem.
   */
  void SetCpModel(const CConfig* config, su2double val_Temperature_Ref) override {
    
    /*--- Read NASA coefficients from the standard polynomial coefficient array (CP_POLYCOEFFS).
          NASA-9 format uses indices 0-8:
          Indices 0-1: Inverse temperature terms (a1*T^-2, a2*T^-1)
          Indices 2-6: Polynomial terms (a3, a4*T, a5*T^2, a6*T^3, a7*T^4)
          Index 7: Enthalpy constant (a8)
          Index 8: Entropy constant (a9)
          
          For NASA-7 compatibility: Set indices 0-1 to zero. ---*/
    for (int i = 0; i < N_COEFFS; ++i) {
      if (i < config->GetnPolyCoeffs()) {
        coeffs_[i] = config->GetCp_PolyCoeff(i);
      } else {
        coeffs_[i] = 0.0; 
      }
    }
    
    Temperature_Min = config->GetTemperatureLimits(0) / Ref_Temp_Dim;
    Temperature_Max = config->GetTemperatureLimits(1) / Ref_Temp_Dim;
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
    const su2double a7 = coeffs_[6];
    const su2double a8 = coeffs_[7];

    // Convert to dimensional temperature for polynomial evaluation (NASA coeffs expect Kelvin)
    const su2double T_dim = t * Ref_Temp_Dim;
    const su2double t_inv = 1.0 / T_dim;
    const su2double t_inv2 = t_inv * t_inv;

    // NASA-9: Cp/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
    su2double Cp_over_R = a1*t_inv2 + a2*t_inv + a3 + a4*T_dim + a5*T_dim*T_dim + a6*T_dim*T_dim*T_dim + a7*T_dim*T_dim*T_dim*T_dim;
    
    Cp = Cp_over_R * Gas_Constant;
    
    // NASA-9: H/(RT) = -a1*T^-2 + a2*ln(T)/T + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + a7*T^4/5 + a8/T
    su2double H_over_RT = -a1*t_inv2 + a2*std::log(T_dim)*t_inv + a3 + a4*T_dim/2.0 + a5*T_dim*T_dim/3.0 + 
                          a6*T_dim*T_dim*T_dim/4.0 + a7*T_dim*T_dim*T_dim*T_dim/5.0 + a8*t_inv;
    
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
    const su2double a7 = coeffs_[6];
    const su2double a8 = coeffs_[7];

    su2double Cp_iter = 0.0;
    su2double delta_temp_iter = 1e10;
    su2double delta_enthalpy_iter;
    const int counter_limit = 50;
    int counter = 0;
    
    while ((abs(delta_temp_iter) > toll) && (counter++ < counter_limit)) {
      
      const su2double T_dim = temp_iter * Ref_Temp_Dim;
      const su2double t_inv = 1.0 / T_dim;
      const su2double t_inv2 = t_inv * t_inv;
      
      // NASA-9: Cp/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
      su2double Cp_over_R = a1*t_inv2 + a2*t_inv + a3 + a4*T_dim + a5*T_dim*T_dim + 
                            a6*T_dim*T_dim*T_dim + a7*T_dim*T_dim*T_dim*T_dim;
      
      Cp_iter = Cp_over_R * Gas_Constant;
      
      // NASA-9: H/(RT) = -a1*T^-2 + a2*ln(T)/T + a3 + a4*T/2 + a5*T^2/3 + a6*T^3/4 + a7*T^4/5 + a8/T
      su2double H_over_RT = -a1*t_inv2 + a2*std::log(T_dim)*t_inv + a3 + a4*T_dim/2.0 + 
                            a5*T_dim*T_dim/3.0 + a6*T_dim*T_dim*T_dim/4.0 + 
                            a7*T_dim*T_dim*T_dim*T_dim/5.0 + a8*t_inv;
      
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
  su2double Ref_Temp_Dim{1.0};     /*!< \brief Dimensional reference temperature for evaluating polynomials. */
  
  std::array<su2double, N_COEFFS> coeffs_;  /*!< \brief NASA polynomial coefficients. */
  
  su2double Temperature_Min{0.0};  /*!< \brief Minimum temperature allowed. */
  su2double Temperature_Max{1e10}; /*!< \brief Maximum temperature allowed. */
};
