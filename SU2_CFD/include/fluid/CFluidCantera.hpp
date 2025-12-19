/*!
 * \file CFluidCantera.hpp
 * \brief  Defines the multicomponent incompressible Ideal Gas model for reacting flows.
 * \author T. Economon, Cristopher Morales Ubal
 * \version 8.1.0 "Harrier"
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

#if defined(HAVE_CANTERA)
#define USE_CANTERA
namespace Cantera {
class Solution;
}
#endif


/*!
 * \class CFluidCantera
 * \brief Child class for defining reacting incompressible ideal gas mixture model.
 * \author: T. Economon, Cristopher Morales Ubal
 */
class CFluidCantera final : public CFluidModel {
 private:
 #ifdef USE_CANTERA
  const int n_species_mixture;            /*!< \brief Number of species in mixture. */
  const su2double Pressure_Thermodynamic; /*!< \brief Constant pressure thermodynamic. */
  const su2double GasConstant_Ref;        /*!< \brief Gas constant reference needed for Nondimensional problems. */
  const su2double Prandtl_Turb_Number;    /*!< \brief Prandlt turbulent number.*/
  const su2double Schmidt_Turb_Number;    /*!< \brief Schmidt turbulent number.*/
  const string Transport_Model;           /*!< \brief Transport model used for computing mixture properties*/
  const string Chemical_MechanismFile;    /*!< \brief Chemical reaction mechanism used for in cantera*/
  const string Phase_Name;                /*!< \brief Name of the phase used for in cantera*/
  const bool Combustion;                  /*!< \brief Flag for problems involving combustion.*/

  static constexpr int ARRAYSIZE = 16;

  su2double Heat_Release;                       /*!< \brief heat release due to combustion */
  std::array<string, ARRAYSIZE> gasComposition; /*!< \brief Gas composition. */
  std::shared_ptr<Cantera::Solution> sol;       /*!< \brief Object needed to describe a chemically-reacting solution*/
  std::array<su2double, ARRAYSIZE> chemicalSourceTerm; /*!< \brief chemical source term of all species*/
  std::array<su2double, ARRAYSIZE> molarMasses;        /*!< \brief Molar masses of all species. */
  std::array<su2double, ARRAYSIZE> enthalpyFormation;  /*!< \brief Enthalpy of Formation of all species. */
  std::array<su2double, ARRAYSIZE> massFractions;      /*!< \brief Mass fractions of all species. */
  std::array<su2double, ARRAYSIZE> massDiffusivity;    /*!< \brief mass diffusivity of all species. */
  mutable vector<su2double> enthalpiesSpecies;                 /*!< \brief Molar enthalpies of all species. */
  mutable vector<su2double> specificHeatSpecies;               /*!< \brief Molar heat capacities of all species. */

  /*!
   * \brief Compute mass diffusivity for species.
   */
  void ComputeMassDiffusivity();

  /*!
   * \brief Compute heat release due to combustion.
   */
  void ComputeHeatRelease();

  /*!
   * \brief Set enthalpies of formation.
   */
  void SetEnthalpyFormation(const CConfig* config);
  #endif

 public:
  /*!
   * \brief Constructor of the class.
   */
  CFluidCantera(su2double val_operating_pressure, const CConfig* config);

  #ifdef USE_CANTERA

  /*!
   * \brief Get fluid laminar viscosity.
   */
  inline su2double GetLaminarViscosity() override { return Mu; }

  /*!
   * \brief Get fluid thermal conductivity.
   */
  inline su2double GetThermalConductivity() override { return Kt + Mu_Turb * Cp / Prandtl_Turb_Number; }

  /*!
   * \brief Get fluid mass diffusivity.
   */
  inline su2double GetMassDiffusivity(int ivar) override { return massDiffusivity[ivar]; }

  /*!
   * \brief Compute chemical source term for species.
   * \param[in] val_scalars - Scalar mass fractions.
   */
  void ComputeChemicalSourceTerm(const su2double* val_scalars) override;

  /*!
   * \brief Get Chemical source term species.
   * \param[in] ivar - index of species.
   */
  inline su2double GetChemicalSourceTerm(int ivar) override { return chemicalSourceTerm[ivar]; }

  /*!
   * \brief Get Heat release due to combustion.
   */
  inline su2double GetHeatRelease() override { return Heat_Release; }

  /*!
   * \brief Get enthalpy diffusivity terms.
   */
  void GetEnthalpyDiffusivity(su2double* enthalpy_diffusions) const override;

  /*!
   * \brief Get enthalpy diffusivity terms.
   */
  void GetMassCorrectionDiffusivity(su2double* MassCorrection_diffusions) override;

  /*!
   * \brief Get gradient enthalpy diffusivity terms.
   */
  void GetGradEnthalpyDiffusivity(su2double* grad_enthalpy_diffusions) const override;

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] t - Temperature value at the point.
   */
  void SetTDState_T(su2double val_temperature, const su2double* val_scalars) override;

  /*!
   * \brief Virtual member.
   * \param[in] val_enthalpy - Enthalpy value at the point.
   * \param[in] val_scalars - Scalar mass fractions.
   */
  void SetTDState_h(su2double val_enthalpy, const su2double* val_scalars = nullptr) override;
  #endif
};