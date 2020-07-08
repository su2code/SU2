/*!
 * \file CNEMOGas.hpp
 * \brief Defines the nonequilibrium gas model.
 * \author W. Maier, C. Garbacz.
 * \version 7.0.5 "Blackbird"
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

#include "CFluidModel.hpp"
#include "../../../Common/include/option_structure.hpp"

/*!
 * \class CNEMOGas
 * \brief Class for defining the 2T (trans-rotational and vibro-electronic) nonequilibrium gas model.
 * \author: W. Maier, C. Garbacz.
 */
class CNEMOGas : public CFluidModel {

protected:
  
  bool viscous;                          /*!< \brief Presence of viscous effects. */
              
  unsigned short nSpecies,               /*!< \brief Number of species in the gas mixture. */
  nHeavy,                                /*!< \brief Number of heavy particles in gas */
  nEl,                                   /*!< \brief Number of electrons in gas */
  nDim,                                  /*!< \brief Number of dimensions. */ 
  nEnergyEq = 2,                         /*!< \brief Number of energy equations for the 2T model. */
  Kind_TransCoeffModel;                  /*!< \brief Transport coefficients model for NEMO solver. */
              
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
  su2double GasConstant{0.0};             /*!< \brief Universal gas constant [J/(kmol*K)] */

  vector<su2double> MassFrac_Freestream;  /*!< \brief Freestream species mass fractions. */

  vector<su2double> MolarMass,           /*!< \brief Species molar mass */
  MassFrac,                                    /*!< \brief Species mass fractions */
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
  DiffusionCoeff;                        /*!< \brief Species diffusion coefficients*/

public:

  /*!
   * \brief Constructor of the class.
   */
  CNEMOGas(const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CNEMOGas(void);

  /*!
   * \brief Set mixture thermodynamic state.
   * \param[in] rhos - Species partial densities.
   * \param[in] T    - Translational/Rotational temperature.
   * \param[in] Tve  - Vibrational/Electronic temperature.
   */
  virtual void SetTDStateRhosTTv(vector<su2double> val_rhos, su2double val_temperature, su2double val_temperature_ve){}

  /*!
   * \brief Set mixture thermodynamic state.
   * \param[in] P    - Pressure.
   * \param[in] T    - Translational/Rotational temperature.
   * \param[in] Tve  - Vibrational/Electronic temperature.
   */
  void SetTDStatePTTv(su2double P, vector<su2double> val_massfrac, su2double val_temperature, su2double val_temperature_ve);
  
  
  /*!
   * \brief Get species T-R specific heats at constant volume.
   */
  virtual vector<su2double> GetSpeciesCvTraRot(){}
  
  /*!
   * \brief Get species V-E specific heats at constant volume.
   */
  virtual vector<su2double> GetSpeciesCvVibEle(){}

  /*!
   * \brief Get specific heat ratio.
   */
  virtual su2double GetGamma(){}
  
  /*!
   * \brief Get mixture energies (total internal energy and vibrational energy).
   */
  //virtual vector<su2double> GetMixtureEnergies(){}
  
  /*!
   * \brief Get species net production rates.
   */
  virtual vector<su2double> GetNetProductionRates(){}
  
  /*!
   * \brief Get vibrational energy source term.
   */
  virtual su2double GetEveSourceTerm(){}

  /*!
   * \brief Get vector of species V-E energy.
   */
  virtual vector<su2double> GetSpeciesEve(){}
  
  /*!
   * \brief Get species enthalpies.
   */
  virtual vector<su2double> GetSpeciesEnthalpy(){}
  
  /*!
   * \brief Get species diffusion coefficients.
   */
  virtual vector<su2double> GetDiffusionCoeff(){}
  
  /*!
   * \brief Get viscosity.
   */
  virtual su2double GetViscosity(){}
  
  /*!
   * \brief Get T-R and V-E thermal conductivities vector.
   */
  virtual vector<su2double> GetThermalConductivities(){}
  
  /*!
   * \brief Get translational and vibrational temperatures vector.
   */
  virtual vector<su2double> GetTemperatures(su2double *rhos, su2double rhoEmix, su2double rhoEve){}
  
  /*!
   * \brief Get speed of sound.
   */
  su2double GetSoundSpeed();

  /*!
   * \brief Get pressure.
   */
  su2double GetPressure();

  /*!
   * \brief Get derivative of pressure w.r.t. conservative variables.
   */
  virtual void GetdPdU(su2double *V, su2double *val_eves, su2double *val_dPdU){}
  
  /*!
   * \brief Get derivative of temperature w.r.t. conservative variables.
   */
  virtual void GetdTdU(su2double *V, su2double *val_dTdU){}
  
  /*!
   * \brief Get derivative of vibrational temperature w.r.t. conservative variables.
   */
  virtual void GetdTvedU(su2double *V, su2double *val_eves, su2double *val_dTvedU){}

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
  inline void SetEves(vector<su2double> val_eves) { eves = val_eves; }

  /*!
   * \brief Get gas constant.
   */
  su2double GetGasConstant();

};