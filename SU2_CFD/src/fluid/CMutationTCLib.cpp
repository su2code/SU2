/*!
 * \file CMutationTCLib.cpp
 * \brief Source of the Mutation++ 2T nonequilibrium gas model.
 * \author C. Garbacz
 * \version 7.0.6 "Blackbird"
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
 
  Mutation::MixtureOptions opt(gas_model);
  string transport_model, state_model;	

  /* Allocating memory*/
  Cv_ks.resize(nEnergyEq*nSpecies,0.0);
  h_RT.resize(nSpecies,0.0);
  es.resize(nEnergyEq*nSpecies,0.0);
  omega_vec.resize(1,0.0);

  /*--- Set up inputs to define type of mixture in the Mutation++ library ---*/

  /*--- Define thermochemical nonequilibrium model ---*/
  if(!frozen)
    state_model = "ChemNonEqTTv";
  else 
    state_model = "NonEqTTv";

  /*--- Define transport model ---*/
  if(Kind_TransCoeffModel == WILKE)
    transport_model = "Wilke";
  else if (Kind_TransCoeffModel == GUPTAYOS)
 	transport_model = "Gupta-Yos";
  
  opt.setStateModel(state_model);
  opt.setViscosityAlgorithm(transport_model);
  opt.setThermalConductivityAlgorithm(transport_model); 
  
  /* Initialize mixture object */
  mix.reset(new Mutation::Mixture(opt));

  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) MolarMass[iSpecies] = 1000* mix->speciesMw(iSpecies); // x1000 to have Molar Mass in kg/kmol

  if (mix->hasElectrons()) { nHeavy = nSpecies-1; nEl = 1; }
  else                     { nHeavy = nSpecies;   nEl = 0; }

}

//CGarbacz returning random things to avoid warnings. This will be properly implemented once NEMO is in develop

CMutationTCLib::~CMutationTCLib(){}
  
void CMutationTCLib::SetTDStateRhosTTv(vector<su2double>& val_rhos, su2double val_temperature, su2double val_temperature_ve){

  temperatures[0] = val_temperature;
  temperatures[1] = val_temperature_ve; 

  T   = temperatures[0];
  Tve = temperatures[1];

  rhos = val_rhos;

  mix->setState(rhos.data(), temperatures.data(), 1);

}

vector<su2double>& CMutationTCLib::GetSpeciesCvTraRot(){
 
   mix->getCvsMass(Cv_ks.data());

   for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) Cvtrs[iSpecies] = Cv_ks[iSpecies];

   return Cvtrs;
}


vector<su2double>& CMutationTCLib::GetSpeciesCvVibEle(){

   mix->getCvsMass(Cv_ks.data());

   for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) Cvves[iSpecies] = Cv_ks[nSpecies+iSpecies];

   return Cvves;
}

vector<su2double>& CMutationTCLib::GetMixtureEnergies(){

  SetTDStateRhosTTv(rhos, T, Tve);

  mix->mixtureEnergies(energies.data());

  return energies; 
}

vector<su2double>& CMutationTCLib::GetSpeciesEve(su2double val_T){

  SetTDStateRhosTTv(rhos, T, val_T);

  mix->getEnergiesMass(es.data());

  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) eves[iSpecies] = es[nSpecies+iSpecies];

  return eves; 
}

vector<su2double>& CMutationTCLib::GetNetProductionRates(){

  mix->netProductionRates(ws.data());

  return ws;
}

su2double CMutationTCLib::GetEveSourceTerm(){

 // for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) cout << setprecision(10)<< "rhos[" << iSpecies << "]=" << rhos[iSpecies] << endl;

 //cout << "GetEveSourceTerm()" << endl;
 //cout << "T=" << temperatures[0] << endl;
 //cout << "Tve=" << temperatures[1] << endl;

  mix->energyTransferSource(omega_vec.data());

 // cout << "omega_vec=" << omega_vec[0] << endl;

  omega = omega_vec[0];

  return omega;
}

vector<su2double>& CMutationTCLib::GetSpeciesEnthalpy(su2double val_T, su2double val_Tve, su2double *val_eves){

  su2double RuSI = UNIVERSAL_GAS_CONSTANT;
  su2double Ru   = 1000.0*RuSI;

  mix->speciesHOverRT(val_T, val_Tve, val_T, val_Tve, val_Tve, h_RT.data(), NULL, NULL, NULL, NULL, NULL);

  //std::cout << "calc mut 3"  << std::endl<< std::endl<< std::endl<< std::endl;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) hs[iSpecies] = h_RT[iSpecies]*(RuSI*val_T); 

  return hs;
}

vector<su2double>& CMutationTCLib::GetDiffusionCoeff(){

  mix->averageDiffusionCoeffs(DiffusionCoeff.data());

  return DiffusionCoeff;  
}

su2double CMutationTCLib::GetViscosity(){

  Mu = mix->viscosity();

  return Mu;  
}

vector<su2double>& CMutationTCLib::GetThermalConductivities(){

  mix->frozenThermalConductivityVector(ThermalConductivities.data());

  return ThermalConductivities;  
}

vector<su2double>& CMutationTCLib::GetTemperatures(vector<su2double>& val_rhos, su2double rhoE, su2double rhoEve, su2double rhoEvel){

  rhos = val_rhos;
  
  energies[0] = rhoE - rhoEvel;
  energies[1] = rhoEve;

 // std::cout << "MUTATION rhoEmix="  << energies[0] << std::endl;
 // std::cout << "MUTATION rhoEve="   << energies[1] << std::endl;
 // for (iSpecies=0;iSpecies<nSpecies;iSpecies++) std::cout << "MUTATION rhos[" << iSpecies << "]="   << rhos[iSpecies] << std::endl;

  mix->setState(rhos.data(), energies.data(), 0);

  mix->getTemperatures(temperatures.data());

 // std::cout << "MUTATION temperatures[0]="  << temperatures[0] << std::endl;
 // std::cout << "MUTATION temperatures[1]="  << temperatures[1] << std::endl;

  T   = temperatures[0];
  Tve = temperatures[1];

  return temperatures;  
}

void CMutationTCLib::GetdPdU(su2double *V, vector<su2double>& val_eves, su2double *val_dPdU){}

void CMutationTCLib::GetdTdU(su2double *V, su2double *val_dTdU){}
 
void CMutationTCLib::GetdTvedU(su2double *V, vector<su2double>& val_eves, su2double *val_dTvedU){}

