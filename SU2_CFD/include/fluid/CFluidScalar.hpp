/*!
 * \file CFluidScalar.hpp
 * \brief  Defines the multicomponent incompressible Ideal Gas model for mixtures.
 * \author T. Economon, Mark Heimgartner, Cristopher Morales Ubal
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
#include <memory>
#include <vector>

#include "CFluidModel.hpp"

/*!
 * \class CFluidScalar
 * \brief Child class for defining an incompressible ideal gas model.
 * \author: T. Economon
 */
class CFluidScalar final : public CFluidModel {
 private:
  unsigned short n_species_mixture; /*!< \brief Number of species in mixture. */
  su2double Gas_Constant;           /*!< \brief Specific gas constant. */
  su2double Gamma;                  /*!< \brief Ratio of specific heats of the gas. */
  su2double Pressure_Thermodynamic; /*!< \brief Constant pressure thermodynamic. */

  bool wilke;
  bool davidson;

  std::vector<su2double> massFractions;              /*!< \brief Mass fractions of all species. */
  std::vector<su2double> moleFractions;              /*!< \brief Mole fractions of all species. */
  std::vector<su2double> molarMasses;                /*!< \brief Molar masses of all species. */
  std::vector<su2double> laminarViscosity;           /*!< \brief Laminar viscosity of all species. */
  std::vector<su2double> laminarThermalConductivity; /*!< \brief Laminar thermal conductivity of all species. */

  static const int ARRAYSIZE = 100;
  std::unique_ptr<CViscosityModel> LaminarViscosityPointers[ARRAYSIZE];
  std::unique_ptr<CConductivityModel> ThermalConductivityPointers[ARRAYSIZE];

  /*!
   * \brief Convert mass fractions to mole fractions.
   * \param[in] val_scalars - Scalar mass fraction.
   */
  std::vector<su2double>& massToMoleFractions(const su2double* val_scalars);

  /*!
   * \brief Wilke mixing law for mixture viscosity.
   * \param[in] val_scalars - Scalar mass fraction.
   */
  su2double wilkeViscosity(const su2double* val_scalars);

  /*!
   * \brief Davidson mixing law for mixture viscosity.
   * \param[in] val_scalars - Scalar mass fraction.
   */
  su2double davidsonViscosity(const su2double* val_scalars);

  /*!
   * \brief Wilke mixing law for mixture thermal conductivity.
   * \param[in] val_scalars - Scalar mass fraction.
   */
  su2double wilkeConductivity(const su2double* val_scalars); 

 public:
  /*!
   * \brief Constructor of the class.
   */
  CFluidScalar(su2double val_Cp, su2double val_gas_constant, su2double val_operating_pressure, CConfig* config);
  
  /*!
   * \brief Set viscosity model.
   */
  void SetLaminarViscosityModel(const CConfig* config);

  /*!
   * \brief Set thermal conductivity model.
   */
  void SetThermalConductivityModel(const CConfig* config);

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] t - Temperature value at the point.
   */
  void SetTDState_T(su2double val_temperature, const su2double* val_scalars) override;
};
