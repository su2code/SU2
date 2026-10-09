/*!
 * \file CFluidCantera.hpp
 * \brief  Defines the multicomponent incompressible Ideal Gas model for reacting flows.
 * \author T. Economon, Cristopher Morales Ubal
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

#include <memory>
#include <string>
#include <vector>

#include "CFluidModel.hpp"

#ifdef HAVE_CANTERA
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
#ifdef HAVE_CANTERA
  const int n_species_mixture;            /*!< \brief Number of species in mixture. */
  const su2double Pressure_Thermodynamic; /*!< \brief Constant pressure thermodynamic. */
  const su2double Prandtl_Turb_Number;    /*!< \brief Prandlt turbulent number.*/
  const su2double Schmidt_Turb_Number;    /*!< \brief Schmidt turbulent number.*/
  const string Transport_Model;           /*!< \brief Transport model used for computing mixture properties*/
  const string Chemical_MechanismFile;    /*!< \brief Chemical reaction mechanism used for in cantera*/
  const string Phase_Name;                /*!< \brief Name of the phase used for in cantera*/
  const bool Combustion;                  /*!< \brief Flag for problems involving combustion.*/
  const su2double Chemistry_Min_Temperature; /*!< \brief Temperature below which the chemistry is skipped. */
  const bool Correction_Velocity;            /*!< \brief The diffusive species fluxes are corrected to sum to zero. */

  su2double Heat_Release{0.0};                 /*!< \brief heat release due to combustion */
  bool stateFailed{false};                     /*!< \brief The last state evaluation did not converge. */
  double minTemperature{0.0};                  /*!< \brief Lower temperature bound of the Newton iteration. */
  double maxTemperature{0.0};                  /*!< \brief Upper temperature bound of the Newton iteration. */
  double temperatureGuess{0.0};                /*!< \brief Start of the next temperature iteration, 0 if unset. */
  std::shared_ptr<Cantera::Solution> sol;      /*!< \brief Object needed to describe a chemically-reacting solution*/
  std::vector<size_t> speciesIndices;          /*!< \brief Mechanism index of each transported species. */
  std::vector<su2double> chemicalSourceTerm;   /*!< \brief Chemical source term of the transported species. */
  std::vector<su2double> chemicalSourceJacobian; /*!< \brief Diagonal sink Jacobian of the transported species. */
  std::vector<su2double> enthalpyFormation;    /*!< \brief Enthalpy of formation of the transported species. */
  std::vector<double> molarMasses;             /*!< \brief Molar masses of all mechanism species. */
  std::vector<double> massFractions;           /*!< \brief Mass fractions of all mechanism species. */
  std::vector<double> netProductionRates;      /*!< \brief Net production rates of all mechanism species. */
  std::vector<double> destructionRates;        /*!< \brief Destruction rates of all mechanism species. */
  mutable std::vector<double> massDiffusivity; /*!< \brief Mass diffusivity of all mechanism species. */
  mutable std::vector<double> enthalpiesSpecies;   /*!< \brief Molar enthalpies of all mechanism species. */
  mutable std::vector<double> specificHeatSpecies; /*!< \brief Molar heat capacities of all mechanism species. */
  mutable su2double laminarViscosity{0.0};         /*!< \brief Laminar viscosity at the current state. */
  mutable su2double laminarConductivity{0.0};      /*!< \brief Laminar thermal conductivity at the current state. */
  mutable bool transportValid{false};              /*!< \brief Viscosity and conductivity match the current state. */
  mutable bool diffusivityValid{false};            /*!< \brief Mass diffusivities match the current state. */

  /*!
   * \brief Set enthalpies of formation.
   */
  void SetEnthalpyFormation(const CConfig* config);

  /*!
   * \brief Set species mass fraction array for Cantera.
   * \param[in] val_scalars - Scalar mass fractions.
   */
  void SetMassFractions(const su2double* val_scalars);

  /*!
   * \brief Set the thermodynamic state at the given temperature and the current composition.
   */
  void SetThermoState(double val_temperature);

  /*!
   * \brief Evaluate viscosity and conductivity of the current state when they are out of date.
   */
  void UpdateTransport() const;

  /*!
   * \brief Evaluate the mass diffusivities of the current state when they are out of date.
   */
  void UpdateDiffusivity() const;
#endif

 public:
  /*!
   * \brief Constructor of the class.
   */
  CFluidCantera(su2double val_operating_pressure, const CConfig* config);

#ifdef HAVE_CANTERA
  /*!
   * \brief Get fluid laminar viscosity.
   */
  inline su2double GetLaminarViscosity() override {
    UpdateTransport();
    return laminarViscosity;
  }

  /*!
   * \brief Get fluid thermal conductivity.
   */
  inline su2double GetThermalConductivity() override {
    UpdateTransport();
    return laminarConductivity + Mu_Turb * Cp / Prandtl_Turb_Number;
  }

  /*!
   * \brief Get fluid mass diffusivity.
   * \param[in] ivar - index of species.
   */
  inline su2double GetMassDiffusivity(int ivar) override {
    UpdateDiffusivity();
    return massDiffusivity[speciesIndices[ivar]];
  }

  /*!
   * \brief Compute chemical source term and heat release from the current state, zero below the minimum temperature.
   */
  void ComputeChemicalSourceTerm() override;

  /*!
   * \brief Get Chemical source term species.
   * \param[in] ivar - index of species.
   */
  inline su2double GetChemicalSourceTerm(int ivar) override { return chemicalSourceTerm[ivar]; }

  /*!
   * \brief Get the derivative of the chemical source term with respect to the own density-weighted mass fraction
   *        (always <= 0).
   * \param[in] ivar - index of species.
   */
  inline su2double GetChemicalSourceJacobian(int ivar) override { return chemicalSourceJacobian[ivar]; }

  /*!
   * \brief Get Heat release due to combustion.
   */
  inline su2double GetHeatRelease() override { return Heat_Release; }

  /*!
   * \brief Whether the last state evaluation failed to converge or needed a temperature far outside the data range.
   */
  inline bool GetStateFailed() const override { return stateFailed; }

  /*!
   * \brief Get enthalpy diffusivity terms.
   */
  void GetEnthalpyDiffusivity(su2double* enthalpy_diffusions) const override;

  /*!
   * \brief Get gradient enthalpy diffusivity terms.
   */
  void GetGradEnthalpyDiffusivity(su2double* grad_enthalpy_diffusions) const override;

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] val_temperature - Temperature value at the point.
   * \param[in] val_scalars - Scalar mass fractions.
   */
  void SetTDState_T(su2double val_temperature, const su2double* val_scalars) override;

  /*!
   * \brief Set the state from enthalpy by Newton iteration on temperature.
   * \param[in] val_enthalpy - Enthalpy value at the point.
   * \param[in] val_scalars - Scalar mass fractions.
   */
  void SetTDState_h(su2double val_enthalpy, const su2double* val_scalars) override;

  /*!
   * \brief Set the start of the next temperature iteration of SetTDState_h, used once.
   * \param[in] val_temperature - Temperature guess.
   */
  inline void SetTemperatureGuess(su2double val_temperature) override {
    temperatureGuess = SU2_TYPE::GetValue(val_temperature);
  }
#endif
};
