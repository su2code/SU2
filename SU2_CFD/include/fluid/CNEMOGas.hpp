/*!
 * \file CNEMOGas.hpp
 * \brief Defines the nonequilibrium gas model.
 * \author C. Garbacz, W. Maier, S. R. Copeland
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
#include "../../../Common/include/option_structure.hpp"

/*!
 * \class CNEMOGas
 * \brief Class for defining the 2T (trans-rotational and vibro-electronic) nonequilibrium gas model.
 * \author: C. Garbacz, W. Maier, S. R. Copeland
 */
class CNEMOGas : public CFluidModel {

protected:

  bool frozen,                           /*!< \brief Indicates if mixture is frozen. */
  ionization;                            /*!< \brief Presence of charged species in gas mixture. */

  string gas_model;                      /*!< \brief String gas model. */

  unsigned short nSpecies,               /*!< \brief Number of species in the gas mixture. */
  nHeavy,                                /*!< \brief Number of heavy particles in gas */
  nEl,                                   /*!< \brief Number of electrons in gas */
  nDim,                                  /*!< \brief Number of dimensions. */
  nEnergyEq = 2;                         /*!< \brief Number of energy equations for the 2T model. */
  TRANSCOEFFMODEL Kind_TransCoeffModel;  /*!< \brief Transport coefficients model for NEMO solver. */

  unsigned iSpecies,                     /*!< \brief Common iteration counter for species */
  jSpecies,                              /*!< \brief Common iteration counter for species */
  iDim;                                  /*!< \brief Common iteration counter for dimensions */

  su2double gamma{0.0};                  /*!< \brief Specific heat ratio. */
  su2double omega{0.0};                  /*!< \brief Vibrational energy source term */
  su2double T{0.0};                      /*!< \brief Translational temperature. */
  su2double Tve{0.0};                    /*!< \brief Vibrational temperature. */
  su2double ThermalCond_tr{0.0};         /*!< \brief T-R thermal conductivity of the gas mixture. */
  su2double ThermalCond_ve{0.0};         /*!< \brief V-E thermal conductivity of the gas mixture. */
  su2double RuSI{UNIVERSAL_GAS_CONSTANT};/*!< \brief Universal gas constant [J/(mol*K)] */
  su2double Ru{1000.0*RuSI};             /*!< \brief Universal gas constant [J/(kmol*K)] */
  su2double GasConstant{0.0};            /*!< \brief Universal gas constant [J/(kmol*K)] */
  su2double rhoCvtr{0.0};                /*!< \brief density times T-R specific heat */
  su2double rhoCvve{0.0};                /*!< \brief density times V-E specific heat */

  vector<su2double> MolarMass,           /*!< \brief Species molar mass */
  MassFrac,                              /*!< \brief Species mass fractions */
  rhos,                                  /*!< \brief Species partial densities */
  Cvtrs,                                 /*!< \brief Species T-R specific heats at constant volume */
  Cvves,                                 /*!< \brief Species V-E specific heats at constant volume */
  energies,                              /*!< \brief Mixture energies (total internal energy and vibrational energy) */
  eves,                                  /*!< \brief Vector of species total internal and V-E energies (dimension 2 * nSpecies)*/
  temperatures,                          /*!< \brief translational and vibrational temperatures vector */
  ThermalConductivities,                 /*!< \brief T-R and V-E Thermal conductivities vector*/
  hs,                                    /*!< \brief Species enthalpies */
  MolarFractions,                        /*!< \brief Species molar fractions */
  ws,                                    /*!< \brief Species net production rates */
  taus,                                  /*!< \brief Relaxtion time scales */
  DiffusionCoeff,                        /*!< \brief Species diffusion coefficients*/
  Enthalpy_Formation,                    /*!< \brief Enthalpy of formation */
  Ref_Temperature;                       /*!< \brief Reference temperature for thermodynamic relations */

  su2matrix<int> CatRecombTable;         /*!< \brief Table for catalytic wall recombination pairs. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CNEMOGas(const CConfig* config, unsigned short val_nDim);

  /*!
   * \brief Set mixture thermodynamic state.
   * \param[in] rhos - Species partial densities.
   * \param[in] T    - Translational/Rotational temperature.
   * \param[in] Tve  - Vibrational/Electronic temperature.
   */
  virtual void SetTDStateRhosTTv(vector<su2double>& val_rhos, su2double val_temperature, su2double val_temperature_ve){}

  /*!
   * \brief Set mixture thermodynamic state.
   * \param[in] P    - Pressure.
   * \param[in] Ms   - Mass fractions of the gas.
   * \param[in] T    - Translational/Rotational temperature.
   * \param[in] Tve  - Vibrational/Electronic temperature.
   */
  void SetTDStatePTTv(su2double P, const su2double *val_massfrac, su2double val_temperature, su2double val_temperature_ve);

  /*!
   * \brief Get species T-R specific heats at constant volume.
   */
  virtual const vector<su2double>& GetSpeciesCvTraRot() = 0;

  /*!
   * \brief Compute species V-E specific heats at constant volume.
   */
  virtual vector<su2double>& ComputeSpeciesCvVibEle(su2double val_T) = 0;

  /*!
   * \brief Compute mixture energies (total internal energy and vibrational energy).
   */
  virtual vector<su2double>& ComputeMixtureEnergies() = 0;

  /*!
   * \brief Compute species net production rates.
   */
  virtual vector<su2double>& ComputeNetProductionRates(bool implicit, const su2double *V, const su2double* eve,
                                                       const su2double* cvve, const su2double* dTdU, const su2double* dTvedU,
                                                       su2double **val_jacobian) = 0;

  /*!
   * \brief Populate chemical source term jacobian.
   */
  virtual void ChemistryJacobian(unsigned short iReaction, const su2double *V, const su2double* eve,
                                 const su2double* cvve, const su2double* dTdU, const su2double* dTvedU,
                                 su2double **val_jacobian){};

  /*!
   * \brief Compute vibrational energy source term.
   */
  virtual su2double ComputeEveSourceTerm() { return 0; }

  /*!
   * \brief Compute vibration enery source term jacobian.
   */
  virtual void GetEveSourceTermJacobian(const su2double *V, const su2double *eve, const su2double *cvve,
                                        const su2double *dTdU, const su2double* dTvedU,
                                        su2double **val_jacobian){};

  /*!
   * \brief Compute vector of species V-E energy.
   */
  virtual vector<su2double>& ComputeSpeciesEve(su2double val_T, bool vibe_only = false) = 0;

  /*!
   * \brief Compute species enthalpies.
   */
  virtual vector<su2double>& ComputeSpeciesEnthalpy(su2double val_T, su2double val_Tve, su2double *val_eves) = 0;

  /*!
   * \brief Get species diffusion coefficients.
   */
  virtual vector<su2double>& GetDiffusionCoeff() = 0;

  /*!
   * \brief Get viscosity.
   */
  virtual su2double GetViscosity() { return 0; }

  /*!
   * \brief Get T-R and V-E thermal conductivities vector.
   */
  virtual vector<su2double>& GetThermalConductivities() = 0;

  /*!
   * \brief Compute translational and vibrational temperatures vector.
   */
  virtual vector<su2double>& ComputeTemperatures(vector<su2double>& val_rhos, su2double rhoEmix, su2double rhoEve, su2double rhoEvel, su2double Tve_old) = 0;

  /*!
   * \brief Compute speed of sound.
   */
  su2double ComputeSoundSpeed();

  /*!
   * \brief Compute pressure.
   */
  su2double ComputePressure();

  /*!
   * \brief Compute gas constant.
   */
  su2double ComputeGasConstant();

  /*!
   * \brief Compute ratio of specific heats (Gamma).
   */
  su2double ComputeGamma();

  /*!
   * \brief Compute derivative of pressure w.r.t. conservative variables.
   */
  void ComputedPdU(const su2double *V, const vector<su2double>& val_eves, su2double *val_dPdU);

  /*!
   * \brief Compute derivative of temperature w.r.t. conservative variables.
   */
  void ComputedTdU(const su2double *V, su2double *val_dTdU);

  /*!
   * \brief Compute derivative of vibrational temperature w.r.t. conservative variables.
   */
  void ComputedTvedU(const su2double *V, const vector<su2double>& val_eves, su2double *val_dTvedU);

  /*!
   * \brief Set the translational temperature.
   */
  inline void SetT(su2double val_temperature) { T = val_temperature; }

  /*!
   * \brief Set the vibrational temperature.
   */
  inline void SetTve(su2double val_temperature_ve) { Tve = val_temperature_ve; }

  /*!
   * \brief Set species vibrational energies.
   */
  inline void SetEves(vector<su2double>& val_eves) { eves = val_eves; }

  /*!
   * \brief Compute rhoCvtr.
   */
  inline su2double ComputerhoCvtr() {
    rhoCvtr = 0.0;
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
      rhoCvtr += rhos[iSpecies]*Cvtrs[iSpecies];
    return rhoCvtr;
  }

  /*!
   * \brief Compute rhoCvve.
   */
  su2double ComputerhoCvve();

  /*!
   * \brief Get species molar mass.
   */
  virtual const vector<su2double>& GetSpeciesMolarMass() = 0;

  /*!
   * \brief Get reference temperature.
   */
  virtual const vector<su2double>& GetRefTemperature() = 0;

  /*!
   * \brief Get species formation enthalpy.
   */
  virtual const vector<su2double>& GetSpeciesFormationEnthalpy() = 0;

  /*!
   * \brief Get catalytic wall recombination indices and constants.
   */
  inline const su2matrix<int>& GetCatalyticRecombination() const {return CatRecombTable;}

};
