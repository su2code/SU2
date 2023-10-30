/*!
 * \file CCoolProp.hpp
 * \brief Defines the state-of-the-art fluid model from CoolProp library.
 * \author P. Yan, G. Gori, A. Guardone
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
#include "CFluidModel.hpp"
#if defined(HAVE_COOLPROP) && !defined(CODI_FORWARD_TYPE) && !defined(CODI_REVERSE_TYPE)
#define USE_COOLPROP
namespace CoolProp {
class AbstractState;
}
#endif
#include <memory>

/*!
 * \class CCoolProp
 * \brief Child class for defining fluid model from CoolProp library.
 * \author: P.Yan
 */
class CCoolProp final : public CFluidModel {
 private:
  su2double Gamma{1.4};                /*!< \brief Ratio of Specific Heats. */
  su2double Gas_Constant{297};         /*!< \brief specific Gas Constant. */
  su2double Pressure_Critical{0.0};    /*!< \brief critical pressure */
  su2double Temperature_Critical{0.0}; /*!< \brief critical temperature */
  su2double acentric_factor{0.0};      /*!< \brief acentric factor */
  const su2double dp{0.01};            /*!< threshold for pressure */
  const su2double dt{0.01};            /*!< threshold for temperature */
#ifdef USE_COOLPROP
  std::unique_ptr<CoolProp::AbstractState> fluid_entity; /*!< \brief fluid entity */
#endif
  /*!
   * \brief Avoid critical pressure
   * \param[in,out] Pressure: Modified so that it is not too close to critical pressure to avoid issues in CoolProp.
   */
  void CheckPressure(su2double& Pressure) const {
    if (Pressure > Pressure_Critical)
      Pressure = fmax(Pressure, (1 + dp) * Pressure_Critical);
    else
      Pressure = fmin(Pressure, (1 - dp) * Pressure_Critical);
  }

  /*!
   * \brief Avoid critical temperature
   * \param[in,out] Temperature: Modified so that it is not too close to critical temperature to avoid issues in
   * CoolProp.
   */
  void CheckTemperature(su2double& Temperature) const {
    if (Temperature > Temperature_Critical)
      Temperature = fmax(Temperature, (1 + dt) * Temperature_Critical);
    else
      Temperature = fmin(Temperature, (1 + dt) * Temperature_Critical);
  }

 public:
  /*!
   * \brief Constructor of the class.
   */
  CCoolProp(const string& fluidname);

#ifdef USE_COOLPROP
  /*!
   * \brief Destructor of the class.
   * \note Needs to be defined in the .cpp to allow using only a forward declaration of CoolProp::AbstractState.
   */
  ~CCoolProp();

  /*!
   * \brief Set the Dimensionless State using Density and Internal Energy
   * \param[in] rho - first thermodynamic variable.
   * \param[in] e - second thermodynamic variable.
   */
  void SetTDState_rhoe(su2double rho, su2double e) override;

  /*!
   * \brief Set the Dimensionless State using Pressure  and Temperature
   * \param[in] P - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  void SetTDState_PT(su2double P, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetTDState_Prho(su2double P, su2double rho) override;

  /*!
   * \brief Set the Dimensionless Internal Energy using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetEnergy_Prho(su2double P, su2double rho) override;

  /*!
   * \brief Set the Dimensionless State using Enthalpy and Entropy
   * \param[in] th1 - first thermodynamic variable (h).
   * \param[in] th2 - second thermodynamic variable (s).
   */
  void SetTDState_hs(su2double h, su2double s) override;

  /*!
   * \brief Set the Dimensionless State using Density and Temperature
   * \param[in] th1 - first thermodynamic variable (rho).
   * \param[in] th2 - second thermodynamic variable (T).
   */
  void SetTDState_rhoT(su2double rho, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Entropy
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (s).
   */
  void SetTDState_Ps(su2double P, su2double s) override;

  /*!
   * \brief compute some derivatives of enthalpy and entropy needed for subsonic inflow BC
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (rho).
   */
  void ComputeDerivativeNRBC_Prho(su2double P, su2double rho) override;
#endif

  /*!
   * \brief Get the value of the critical pressure.
   * \return Critical pressure.
   */
  su2double GetPressure_Critical(void) const { return Pressure_Critical; }

  /*!
   * \brief Get the value of the critical temperature.
   * \return Critical temperature.
   */
  su2double GetTemperature_Critical(void) const { return Temperature_Critical; }

  /*!
   * \brief Get the value of specific gas constant.
   * \return Value of the constant: Gamma
   */
  su2double GetGas_Constant(void) const { return Gas_Constant; }

  /*!
   * \brief Get the value of specific gas constant.
   * \return Value of the constant: Gamma
   */
  su2double GetGamma(void) const { return Gamma; }
};