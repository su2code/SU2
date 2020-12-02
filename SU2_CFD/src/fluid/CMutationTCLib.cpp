/*!
 * \file CMutationTCLib.cpp
 * \brief Source of the Mutation++ 2T nonequilibrium gas model.
 * \author C. Garbacz
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

#include "../../include/fluid/CMutationTCLib.hpp"

CMutationTCLib::CMutationTCLib(const CConfig* config, unsigned short val_nDim): CNEMOGas(config, val_nDim){
  
  //CGarbacz: if wilke - transportmodel = 'wilke' and so on;

  //CGarbacz: nEl = mix.getnumberelectrons; nHeavy = nSpecies-nEl;

  //CGarbacz: if frozen send nonchecmttv

}

//CGarbacz returning random things to avoid warnings. This will be properly implemented once Mutation++ is in develop

CMutationTCLib::~CMutationTCLib(){}
  
void CMutationTCLib::SetTDStateRhosTTv(vector<su2double>& val_rhos, su2double val_temperature, su2double val_temperature_ve){}

vector<su2double>& CMutationTCLib::GetSpeciesMolarMass(){return MassFrac;}

vector<su2double>& CMutationTCLib::GetSpeciesCvTraRot(){return MassFrac;}

vector<su2double>& CMutationTCLib::GetSpeciesCvVibEle(){return MassFrac;}

vector<su2double>& CMutationTCLib::GetMixtureEnergies(){return MassFrac;}

vector<su2double>& CMutationTCLib::GetSpeciesEve(su2double val_T){return MassFrac;}

vector<su2double>& CMutationTCLib::GetNetProductionRates(){return MassFrac;}

su2double CMutationTCLib::GetEveSourceTerm(){return 0;}

vector<su2double>& CMutationTCLib::GetSpeciesEnthalpy(su2double val_T, su2double val_Tve, su2double *val_eves){return MassFrac;}

vector<su2double>& CMutationTCLib::GetDiffusionCoeff(){return MassFrac;}

su2double CMutationTCLib::GetViscosity(){return 0;}

vector<su2double>& CMutationTCLib::GetThermalConductivities(){return MassFrac;}

vector<su2double>& CMutationTCLib::GetTemperatures(vector<su2double>& rhos, su2double rhoEmix, su2double rhoEve, su2double rhoEvel){return MassFrac;}

vector<su2double>& CMutationTCLib::GetRefTemperature() {return MassFrac;}

vector<su2double>& CMutationTCLib::GetSpeciesFormationEnthalpy() {return MassFrac;}
