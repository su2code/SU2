/*!
 * \file CMutationTCLib.cpp
 * \brief Source of the Mutation++ 2T nonequilibrium gas model.
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

#include "../../include/fluid/CMutationTCLib.hpp"

CMutationTCLib::CMutationTCLib(const CConfig* config): CNEMOGas(config){
  
  //cat: if wilke - transportmodel = 'wilke' and so on;

  // nEl = mix.getnumberelectrons; nHeavy = nSpecies-nEl;

	//if frozen send nonchecmttv

}

CMutationTCLib::~CMutationTCLib(){}
  
void CMutationTCLib::SetTDStateRhosTTv(vector<su2double> val_rhos, su2double val_temperature, su2double val_temperature_ve){}

vector<su2double> CMutationTCLib::GetSpeciesCvTraRot(){}

vector<su2double> CMutationTCLib::GetSpeciesCvVibEle(){}

su2double CMutationTCLib::GetGamma(){}

vector<su2double> CMutationTCLib::GetMixtureEnergies(){}

vector<su2double> CMutationTCLib::GetSpeciesEve(su2double val_T){}

vector<su2double> CMutationTCLib::GetNetProductionRates(){}

su2double CMutationTCLib::GetEveSourceTerm(){}

vector<su2double> CMutationTCLib::GetSpeciesEnthalpy(su2double val_T, su2double *val_eves){}

vector<su2double> CMutationTCLib::GetDiffusionCoeff(){}

su2double CMutationTCLib::GetViscosity(){}

vector<su2double> CMutationTCLib::GetThermalConductivities(){}

vector<su2double> CMutationTCLib::GetTemperatures(vector<su2double> rhos, su2double rhoEmix, su2double rhoEve){}

void CMutationTCLib::GetdPdU(su2double *V, vector<su2double> val_eves, su2double *val_dPdU){}

void CMutationTCLib::GetdTdU(su2double *V, su2double *val_dTdU){}
 
void CMutationTCLib::GetdTvedU(su2double *V, vector<su2double> val_eves, su2double *val_dTvedU){}

