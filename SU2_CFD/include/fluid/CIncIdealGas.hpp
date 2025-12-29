/*!
 * \file CIncIdealGas.hpp
 * \brief Defines the incompressible Ideal Gas model.
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

#include "CFluidModel.hpp"

/*!
 * \class CIncIdealGas
 * \brief Child class for defining an incompressible ideal gas model.
 * \author: T. Economon
 */
class CIncIdealGas final : public CFluidModel {
 public:
  /*!
   * \brief Constructor of the class.
   */
  CIncIdealGas(su2double val_Cp, su2double val_gas_constant, su2double val_operating_pressure, su2double val_Temperature_Ref) {
    /*--- In the incompressible ideal gas model, the thermodynamic pressure
  is decoupled from the governing equations and held constant. The
  density is therefore only a function of temperature variations. ---*/
    Gas_Constant = val_gas_constant;
    Pressure = val_operating_pressure;
    Gamma = 1.0;
    Cp = val_Cp;
    Cv = Cp;
    Std_Ref_Temp_ND = STD_REF_TEMP / val_Temperature_Ref;
  }

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] t - Temperature value at the point.
   */
  void SetTDState_T(su2double t, const su2double *val_scalars = nullptr) override {
    /*--- The EoS only depends upon temperature. ---*/
    Temperature = t;
    Density = Pressure / (Temperature * Gas_Constant);
    Enthalpy = Cp * (Temperature - Std_Ref_Temp_ND);  // Sensible enthalpy relative to REF_TEMP
  }

  /*!
   * \brief Set the Dimensionless State using Enthalpy.
   * \param[in] val_enthalpy - Enthalpy value at the point.
   */
  void SetTDState_h(su2double val_enthalpy, const su2double* val_scalars = nullptr) override {
    Enthalpy = val_enthalpy;
    Temperature = Enthalpy / Cp + Std_Ref_Temp_ND;  // Temperature from sensible enthalpy
    Density = Pressure / (Temperature * Gas_Constant);
  }

 private:
  su2double Gas_Constant{0.0}; /*!< \brief Gas Constant. */
  su2double Gamma{0.0};        /*!< \brief Heat Capacity Ratio. */
  su2double Std_Ref_Temp_ND{0.0}; /*!< \brief Nondimensional standard reference temperature for enthalpy. */
};
