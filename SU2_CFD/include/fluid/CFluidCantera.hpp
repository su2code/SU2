/*!
 * \file CFluidCantera.hpp
 * \brief  Defines the multicomponent incompressible Ideal Gas model for reacting flows.
 * \author T. Economon, Cristopher Morales Ubal
 * \version 8.0.1 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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
#include <array>

#include "CFluidModel.hpp"

/*!
 * \class CFluidCantera
 * \brief Child class for defining reacting incompressible ideal gas mixture model.
 * \author: T. Economon
 */
class CFluidCantera final : public CFluidModel {
 private:
  const int n_species_mixture;            /*!< \brief Number of species in mixture. */
  su2double Gas_Constant;           /*!< \brief Specific gas constant. */
  const su2double Gamma;                  /*!< \brief Ratio of specific heats of the gas. */
  const su2double Pressure_Thermodynamic; /*!< \brief Constant pressure thermodynamic. */
  const su2double GasConstant_Ref;        /*!< \brief Gas constant reference needed for Nondimensional problems. */
  const su2double Prandtl_Number;         /*!< \brief Prandlt number.*/

  static constexpr int ARRAYSIZE = 16;

  std::array<su2double, ARRAYSIZE> massFractions;              /*!< \brief Mass fractions of all species. */
  std::array<su2double, ARRAYSIZE> moleFractions;              /*!< \brief Mole fractions of all species. */
  std::array<su2double, ARRAYSIZE> molarMasses;                /*!< \brief Molar masses of all species. */
  std::array<su2double, ARRAYSIZE> massDiffusivity;           /*!< \brief mass diffusivity of all species. */
  std::unique_ptr<CDiffusivityModel> MassDiffusivityPointers[ARRAYSIZE];

  /*!
   * \brief Convert mass fractions to mole fractions.
   * \param[in] val_scalars - Scalar mass fraction.
   */
  void MassToMoleFractions(const su2double* val_scalars);

  /*!
   * \brief Compute gas constant for mixture.
   */
  su2double ComputeGasConstant();

  /*!
   * \brief Compute mass diffusivity for species.
   */
  void ComputeMassDiffusivity();

 public:
  /*!
   * \brief Constructor of the class.
   */
  CFluidCantera(su2double val_Cp, su2double val_gas_constant, su2double val_operating_pressure, const CConfig* config);


  /*!
   * \brief Set mass diffusivity model.
   */
  void SetMassDiffusivityModel(const CConfig* config) override;

  /*!
   * \brief Get fluid laminar viscosity.
   */
  inline su2double GetLaminarViscosity() override { return Mu; }

  /*!
   * \brief Get fluid thermal conductivity.
   */
  inline su2double GetThermalConductivity() override { return Kt + Mu_Turb * Cp / Prandtl_Number; }

  /*!
   * \brief Get fluid mass diffusivity.
   */
  inline su2double GetMassDiffusivity(int ivar) override { return massDiffusivity[ivar]; }

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] t - Temperature value at the point.
   */
  void SetTDState_T(su2double val_temperature, const su2double* val_scalars) override;
};