/*!
 * \file CUserDefinedTCLib.hpp
 * \brief Defines the classes for different user defined ThermoChemistry libraries.
 * \author C. Garbacz, W. Maier, S. R. Copeland
 * \version 7.0.8 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "CNEMOGas.hpp"

/*!
 * \derived class CUserDefinedTCLib
 * \brief Child class for user defined nonequilibrium gas model.
 * \author: C. Garbacz, W. Maier, S. R. Copeland
 */
class CUserDefinedTCLib : public CNEMOGas {

private:

  unsigned short nReactions,        /*!< \brief Number of reactions in chemical model. */
  iEl;                              /*!< \brief Common iteration counter for electrons */

  vector<unsigned short> nElStates; /*!< \brief Number of electron states. */

  C3DIntMatrix Reactions;           /*!</brief reaction map for chemically reacting flows */
  
  vector<su2double>
  ArrheniusCoefficient,             /*!< \brief Arrhenius reaction coefficient */
  ArrheniusEta,                     /*!< \brief Arrhenius reaction temperature exponent */
  ArrheniusTheta,                   /*!< \brief Arrhenius reaction characteristic temperature */
  CharVibTemp,                      /*!< \brief Characteristic vibrational temperature for e_vib */
  RotationModes,	          /*!< \brief Rotational modes of energy storage */
  Tcf_a,                          /*!< \brief Rate controlling temperature exponent (fwd) */
  Tcf_b,                          /*!< \brief Rate controlling temperature exponent (fwd) */
  Tcb_a,                          /*!< \brief Rate controlling temperature exponent (bkw) */
  Tcb_b,                          /*!< \brief Rate controlling temperature exponent (bkw) */
  Diss,                           /*!< \brief Dissociation potential. */
  MassFrac_FreeStream,            /*!< \brief Mixture mass fractions of the fluid. */
  Wall_Catalycity,                /*!< \brief Specified wall species mass-fractions for catalytic boundaries. */
  Particle_Mass,                  /*!< \brief Mass of all particles present in the plasma */
  MolarFracWBE,                   /*!< \brief Molar fractions to be used in Wilke/Blottner/Eucken model */
  phis, mus,                      /*!< \brief Auxiliary vectors to be used in Wilke/Blottner/Eucken model */
  A;                              /*!< \brief Auxiliary vector to be used in net production rate computation */

  su2activematrix CharElTemp,    /*!< \brief Characteristic temperature of electron states. */
  ElDegeneracy,                  /*!< \brief Degeneracy of electron states. */
  RxnConstantTable,              /*!< \brief Table of chemical equiibrium reaction constants */
  Blottner,                      /*!< \brief Blottner viscosity coefficients */
  Dij;                           /*!< \brief Binary diffusion coefficients. */
  
  C3DDoubleMatrix Omega00,       /*!< \brief Collision integrals (Omega(0,0)) */
  Omega11;                       /*!< \brief Collision integrals (Omega(1,1)) */

public:

  /*!
   * \brief Constructor of the class.
   */
  CUserDefinedTCLib(const CConfig* config, unsigned short val_nDim, bool val_viscous);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CUserDefinedTCLib(void);

  /*!
   * \brief Set mixture thermodynamic state.
   * \param[in] rhos - Species partial densities.
   * \param[in] T    - Translational/Rotational temperature.
   * \param[in] Tve  - Vibrational/Electronic temperature.
   */
  void SetTDStateRhosTTv(vector<su2double>& val_rhos, su2double val_temperature, su2double val_temperature_ve) final;

  /*!
   * \brief Get species molar mass.
   */
  vector<su2double>& GetSpeciesMolarMass() final { return MolarMass; }

  /*!
   * \brief Get species T-R specific heats at constant volume.
   */
  vector<su2double>& GetSpeciesCvTraRot() final;

  /*!
   * \brief Get species V-E specific heats at constant volume.
   */
  vector<su2double>& GetSpeciesCvVibEle() final;

  /*!
   * \brief Get mixture energies (total internal energy and vibrational energy).
   */
  vector<su2double>& GetMixtureEnergies() final;

  /*!
   * \brief Get vector of species V-E energy.
   */
  vector<su2double>& GetSpeciesEve(su2double val_T) final;

  /*!
   * \brief Get species net production rates.
   */
  vector<su2double>& GetNetProductionRates() final;

  /*!
   * \brief Get vibrational energy source term.
   */
  su2double GetEveSourceTerm() final;

  /*!
   * \brief Get species enthalpies.
   */
  vector<su2double>& GetSpeciesEnthalpy(su2double val_T, su2double val_Tve, su2double *val_eves) final;

  /*!
   * \brief Get species diffusion coefficients.
   */
  vector<su2double>& GetDiffusionCoeff() final;

  /*!
   * \brief Get viscosity.
   */
  su2double GetViscosity() final;

  /*!
   * \brief Get T-R and V-E thermal conductivities vector.
   */
  vector<su2double>& GetThermalConductivities() final;

  /*!
   * \brief Get translational and vibrational temperatures vector.
   */
  vector<su2double>& GetTemperatures(vector<su2double>& val_rhos, su2double rhoEmix, su2double rhoEve, su2double rhoEvel) final;

  private:

  /*!
  * \brief Set equilibrium reaction constants for finite-rate chemistry
  */
  void GetChemistryEquilConstants(unsigned short iReaction);

  /*!
   * \brief Calculates constants used for Keq correlation.
   * \param[out] A - Reference to coefficient array.
   * \param[in] val_reaction - Reaction number indicator.
   */
  void GetKeqConstants(unsigned short val_Reaction);

  /*!
   * \brief Get species diffusion coefficients with Wilke/Blottner/Eucken transport model.
   */
  void DiffusionCoeffWBE();

  /*!
   * \brief Get viscosity with Wilke/Blottner/Eucken transport model.
   */
  void ViscosityWBE();

  /*!
   * \brief Get T-R and V-E thermal conductivities vector with Wilke/Blottner/Eucken transport model.
   */
  void ThermalConductivitiesWBE();

  /*!
   * \brief Get species diffusion coefficients with Gupta-Yos transport model.
   */
  void DiffusionCoeffGY();

  /*!
   * \brief Get viscosity with Gupta-Yos transport model.
   */
  void ViscosityGY();

  /*!
   * \brief Get T-R and V-E thermal conductivities vector with Gupta-Yos transport model.
   */
  void ThermalConductivitiesGY();

  /*!
   * \brief Get reference temperature.
   */
  vector<su2double>& GetRefTemperature() final { return Ref_Temperature; }

  /*!
   * \brief Get species formation enthalpy.
   */
  vector<su2double>& GetSpeciesFormationEnthalpy() final { return Enthalpy_Formation; }  

  };

