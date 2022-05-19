/*!
 * \file CFluidScalar.hpp
 * \brief Defines the incompressible Ideal Gas model.
 * \author T. Economon
 * \version 7.3.1 "Blackbird"
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
#include <vector>
#include <memory>

#include "CFluidModel.hpp"

/*!
 * \class CFluidScalar
 * \brief Child class for defining an incompressible ideal gas model.
 * \author: T. Economon
 */
class CFluidScalar final : public CFluidModel {
private:
  unsigned short n_scalars = 0;                       /*!< \brief Number of transported scalars. */
  unsigned short n_species_mixture = 0;               /*!< \brief Number of species in mixture. */
  su2double Gas_Constant = 0.0;                       /*!< \brief Specific gas constant. */
  su2double Gamma = 0.0;                              /*!< \brief Ratio of specific heats of the gas. */
  su2double Pressure_Thermodynamic = 0.0;             /*!< \brief Constant pressure thermodynamic. */

  bool wilke;
  bool davidson;

  std::vector<su2double> massFractions;               /*!< \brief Mass fractions of all species. */
  std::vector<su2double> moleFractions;               /*!< \brief Mole fractions of all species. */
  std::vector<su2double> molarMasses;                 /*!< \brief Molar masses of all species. */
  std::vector<su2double> laminarViscosity;            /*!< \brief Laminar viscosity of all species. */
  std::vector<su2double> specificHeat;                /*!< \brief Specific heat of all species. */
  std::vector<su2double> laminarthermalConductivity;  /*!< \brief Laminar thermal conductivity of all species. */

  static constexpr int ARRAYSIZE = 100;
  std::unique_ptr<CViscosityModel> LaminarViscosityPointers[ARRAYSIZE];
  std::unique_ptr<CConductivityModel> ThermalConductivityPointers[ARRAYSIZE];

 public:
  /*!
   * \brief Constructor of the class.
   */
  CFluidScalar(su2double val_Cp, su2double val_gas_constant, su2double val_operating_pressure, CConfig *config);
  
  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] t - Temperature value at the point.
   */
  void SetTDState_T(su2double t){
    /*--- The EoS only depends upon temperature. ---*/
    Temperature = t;
    Density = Pressure / (Temperature * Gas_Constant);
  }
  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] t - Temperature value at the point.
   */

  void SetTDState_T(su2double val_temperature, const su2double* val_scalars) override;

};
