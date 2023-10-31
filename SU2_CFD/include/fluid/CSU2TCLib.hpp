/*!
 * \file CSU2TCLib.hpp
 * \brief Defines the classes for different user defined ThermoChemistry libraries.
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

#include "CNEMOGas.hpp"

/*!
 * \derived class CSU2TCLib
 * \brief Child class for user defined nonequilibrium gas model.
 * \author: C. Garbacz, W. Maier, S. R. Copeland
 */
class CSU2TCLib : public CNEMOGas {

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
  RotationModes,                  /*!< \brief Rotational modes of energy storage */
  Tcf_a,                          /*!< \brief Rate controlling temperature exponent (fwd) */
  Tcf_b,                          /*!< \brief Rate controlling temperature exponent (fwd) */
  Tcb_a,                          /*!< \brief Rate controlling temperature exponent (bkw) */
  Tcb_b,                          /*!< \brief Rate controlling temperature exponent (bkw) */
  taus,                           /*!< \brief Relaxtion time scales */
  Diss,                           /*!< \brief Dissociation potential. */
  MassFrac_FreeStream,            /*!< \brief Mixture mass fractions of the fluid. */
  Wall_Catalycity,                /*!< \brief Specified wall species mass-fractions for catalytic boundaries. */
  Particle_Mass,                  /*!< \brief Mass of all particles present in the plasma */
  MolarFracWBE,                   /*!< \brief Molar fractions to be used in Wilke/Blottner/Eucken model */
  phis, mus,                      /*!< \brief Auxiliary vectors to be used in Wilke/Blottner/Eucken model */
  A;                              /*!< \brief Auxiliary vector to be used in net production rate computation */

  std::array<su2double,1> mu_ref; /*!< \brief Vector containing reference viscosity for Sutherland's law */
  std::array<su2double,1> k_ref;  /*!< \brief Vector containing reference thermal conducivities for Sutherland's law */
  std::array<su2double,1> Sm_ref; /*!< \brief Vector containing Sutherland's constant for viscosity */
  std::array<su2double,1> Sk_ref; /*!< \brief Vector containing Sutherland's constant for thermal conductivities */

  const su2double T_ref_suth = 273.15; /*!<\brief Reference temperature for Sutherland's model [K] */
  
  // Coulomb potential constant values source: Scalabrin, NUMERICAL SIMULATION OF WEAKLY IONIZED HYPERSONIC
  // FLOW OVER REENTRY CAPSULES, 2007.
  /*--- Attractive Coulombic potential constants ---*/
  const su2double D1_a = 0.784;
  const su2double C1_a = -0.476;
  const su2double c1_a = 0.0313;
  const su2double D2_a = 1.262;
  const su2double C2_a = -0.146;
  const su2double c2_a = 0.0377;

  /*--- Repulsive Coulombic potential constants ---*/
  const su2double D1_r = 0.765;
  const su2double C1_r = 0.138;
  const su2double c1_r = 0.0106;
  const su2double D2_r = 1.235;
  const su2double C2_r = 0.157;
  const su2double c2_r = 0.0274;

  su2activematrix CharElTemp,    /*!< \brief Characteristic temperature of electron states. */
  ElDegeneracy,                  /*!< \brief Degeneracy of electron states. */
  RxnConstantTable,              /*!< \brief Table of chemical equiibrium reaction constants */
  Blottner,                      /*!< \brief Blottner viscosity coefficients */
  Dij;                           /*!< \brief Binary diffusion coefficients. */

  C3DDoubleMatrix Omega11,       /*!< \brief Collision integrals (Omega^(1,1)) */
  Omega22;                       /*!< \brief Collision integrals (Omega^(2,2)) */

  /*--- Implicit variables ---*/
  su2double                     /*!< \brief Derivatives w.r.t. conservative variables */
  *dPdU, *dTdU, *dTvedU;

  su2double fwdRxn, bkwRxn,
  kf,kfb,kb,
  coeff, eta, epsilon, T_min,
  Trxnf, Trxnb,
  Thf, Thb, dThf, dThb,
  theta, af, bf, ab, bb;

  vector<su2double>
  dkf, dkb,
  dRfok, dRbok,
  eve, eve_eq, cvve, cvve_eq;

  vector<int>
  alphak, betak;

public:

  /*!
   * \brief Constructor of the class.
   */
  CSU2TCLib(const CConfig* config, unsigned short val_nDim, bool val_viscous);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CSU2TCLib(void);

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
   * \brief Compute species V-E specific heats at constant volume.
   */
  vector<su2double>& ComputeSpeciesCvVibEle(su2double val_T) final;

  /*!
   * \brief Compute mixture energies (total internal energy and vibrational energy).
   */
  vector<su2double>& ComputeMixtureEnergies() final;

  /*!
   * \brief Compute species V-E energy.
   */
  vector<su2double>& ComputeSpeciesEve(su2double val_T, bool vibe_only = false) final;

  /*!
   * \brief Compute species net production rates.
   */
  vector<su2double>& ComputeNetProductionRates(bool implicit, const su2double *V, const su2double* eve,
                                               const su2double* cvve, const su2double* dTdU, const su2double* dTvedU,
                                               su2double **val_jacobian) final;

  /*!
   * \brief Compute chemical source term jacobian.
   */
  void ChemistryJacobian(unsigned short iReaction, const su2double *V, const su2double* eve,
                         const su2double* cvve, const su2double* dTdU, const su2double* dTvedU,
                         su2double **val_jacobian) final;

  /*!
   * \brief Compute vibrational energy source term.
   */
  su2double ComputeEveSourceTerm() final;

  /*!
   * \brief Compute relaxation source term jacobian.
   */
  void GetEveSourceTermJacobian(const su2double *V, const su2double *eve, const su2double *cvve,
                                const su2double *dTdU, const su2double* dTvedU,
                                su2double **val_jacobian) final;

  /*!
   * \brief Compute species enthalpies.
   */
  vector<su2double>& ComputeSpeciesEnthalpy(su2double val_T, su2double val_Tve, su2double *val_eves) final;

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
   * \brief Compute translational and vibrational temperatures vector.
   */
  vector<su2double>& ComputeTemperatures(vector<su2double>& val_rhos, su2double rhoEmix, su2double rhoEve, su2double rhoEvel, su2double Tve_old) final;

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
  void ComputeKeqConstants(unsigned short val_Reaction);

  /*!
   * \brief Calculate species diffusion coefficients with Wilke/Blottner/Eucken transport model.
   */
  void DiffusionCoeffWBE();

  /*!
   * \brief Calculate viscosity with Wilke/Blottner/Eucken transport model.
   */
  void ViscosityWBE();

  /*!
   * \brief Calculate T-R and V-E thermal conductivities vector with Wilke/Blottner/Eucken transport model.
   */
  void ThermalConductivitiesWBE();

  /*!
   * \brief Calculate species diffusion coefficients with Gupta-Yos transport model.
   */
  void DiffusionCoeffGY();

  /*!
   * \brief Calculate viscosity with Gupta-Yos transport model.
   */
  void ViscosityGY();

  /*!
   * \brief Calculate T-R and V-E thermal conductivities vector with Gupta-Yos transport model.
   */
  void ThermalConductivitiesGY();

  /*!
   * \brief Calculate viscosity with Sutherland's transport model.
   */
  void ViscositySuth();

  /*!
   * \brief Calculate T-R and V-E thermal conductivities vector with Sutherland's transport model.
   */
  void ThermalConductivitiesSuth();

  /*!
   *\brief Compute the collision terms (deltas)
   *\param[in] iSpecies - Species of gas
   *\param[in] jSpecies - Collision partner of species i
   *\param[in] Mi - Molar mass of species i
   *\param[in] Mj - Molar mass of species j
   *\param[in] T - Temperature of gas
   *\param[in] d1 - Whether delta_1 or delta_2 is computed
   *\param[out] Delta_1 or Delta_2 - Collision term between species i and j
   */
  su2double ComputeCollisionDelta(unsigned iSpecies, unsigned jSpecies, su2double Mi, su2double Mj, su2double T, bool d1);
  
  /*!
   *\brief Compute the collision cross section
   *\param[in] iSpecies - Species of gas
   *\param[in] jSpecies - Collision partner of species i
   *\param[in] T - Temperature of gas
   *\param[in] d1 - Whether omega^(1,1) or omega^(2,2)
   *\param[in] coulomb - Whether to use Coulomb potential
   *\param[out] Omega_ij - collision cross section
   */
  su2double ComputeCollisionCrossSection(unsigned iSpecies, unsigned jSpecies, su2double T, bool d1, bool coulomb);

  /*!
   * \brief Get reference temperature.
   */
  vector<su2double>& GetRefTemperature() final { return Ref_Temperature; }

  /*!
   * \brief Get species formation enthalpy.
   */
  vector<su2double>& GetSpeciesFormationEnthalpy() final { return Enthalpy_Formation; }

  };

