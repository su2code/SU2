/*!
 * \file CMutationTCLib.hpp
 * \brief Defines the class for the link to Mutation++ ThermoChemistry library.
 * \author W. Maier, C. Garbacz
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

#include "CNEMOGas.hpp"

/*!
 * \derived class CMutationTCLib
 * \brief Child class for Mutation++ nonequilibrium gas model.
 * \author: W. Maier
 */
class CMutationTCLib : public CNEMOGas {

private:

  string GasModel;

public:

  /*!
   * \brief Constructor of the class.
   */
  CMutationTCLib(const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CMutationTCLib(void);

  /*!
   * \brief Set mixture therodynamic state.
   * \param[in] rhos - Species partial densities.
   * \param[in] T    - Translational/Rotational temperature.
   * \param[in] Tve  - Vibrational/Electronic temperature.
   */
  void SetTDStateRhosTTv(vector<su2double> val_rhos, su2double val_temperature, su2double val_temperature_ve);

  /*!
   * \brief Get species T-R specific heats at constant volume.
   */
  vector<su2double> GetSpeciesCvTraRot();

  /*!
   * \brief Get species V-E specific heats at constant volume.
   */
  vector<su2double> GetSpeciesCvVibEle();
  
  /*!
   * \brief Get specific heat ratio.
   */
  su2double GetGamma();
  
  /*!
   * \brief Get mixture energies (total internal energy and vibrational energy).
   */
  vector<su2double> GetMixtureEnergies();

  /*!
   * \brief Get vector of species V-E energy.
   */
  vector<su2double> GetSpeciesEve(su2double val_T);
  
  /*!
   * \brief Get species net production rates.
   */
  vector<su2double> GetNetProductionRates();

  /*!
   * \brief Get vibrational energy source term.
   */
  su2double GetEveSourceTerm();
  
  /*!
   * \brief Get species enthalpies.
   */
  vector<su2double> GetSpeciesEnthalpy(su2double val_T, su2double *val_eves);

  /*!
   * \brief Get species diffusion coefficients.
   */
  vector<su2double> GetDiffusionCoeff();

  /*!
   * \brief Get viscosity.
   */
  su2double GetViscosity();

  
  /*!
   * \brief Get T-R and V-E thermal conductivities vector.
   */
  vector<su2double> GetThermalConductivities();
  
  /*!
   * \brief Get translational and vibrational temperatures vector.
   */
  vector<su2double> GetTemperatures(vector<su2double> rhos, su2double rhoEmix, su2double rhoEve, su2double rhoEvel);
   
  /*!
   * \brief Get derivative of pressure w.r.t. conservative variables.
   */
  void GetdPdU(su2double *V, vector<su2double> val_eves, su2double *val_dPdU);

  /*!
   * \brief Get derivative of temperature w.r.t. conservative variables.
   */
  void GetdTdU(su2double *V, su2double *val_dTdU);

  /*!
   * \brief Get derivative of vibrational temperature w.r.t. conservative variables.
   */
  void GetdTvedU(su2double *V, vector<su2double> val_eves, su2double *val_dTvedU);

};