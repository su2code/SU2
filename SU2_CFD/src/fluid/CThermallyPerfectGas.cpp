/*!
 * \file CThermallyPerfectGas.cpp
 * \brief Source of the ideal gas model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
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

#include "../../include/fluid/CThermallyPerfectGas.hpp"

CThermallyPerfectGas::CThermallyPerfectGas(const CConfig* config, su2double R, bool CompEntropy) : CFluidModel() {
  Mutation::MixtureOptions opt(config->GetGasModel());
  opt.setStateModel("ChemNonEq1T");
  opt.setThermodynamicDatabase("RRHO");
  opt.setMechanism("none");
  mix.reset(new Mutation::Mixture(opt));
  Gas_Constant = R;
  ComputeEntropy = CompEntropy;
  ns = config->GetnSpecies();
  cs = config->GetGas_Composition();
  rhos.resize(ns,0.0);
  energy.resize(1,0.0);
  rhoenergy.resize(1,0.0);
  PT.resize(2,0.0);
  temp.resize(1,0.0);
  MolarMass.resize(ns,0.0);
  Cv_ks.resize(4,0.0);
  Cvtrs.resize(ns,0.0);         
  Cvves.resize(ns,0.0);              
  Gas_Constant = R;
  for(i = 0; i < ns; i++) MolarMass[i] = 1000* mix->speciesMw(i);

  
  Gamma = 1.4;
  Gamma_Minus_One = Gamma - 1.0;
  Gas_Constant = R;
  Cp = Gamma / Gamma_Minus_One * Gas_Constant;
  Cv = Cp - R;
  
}

void CThermallyPerfectGas::SetTDState_rhoe(su2double rho, su2double e) {

  //cout << endl << "SetTDState_rhoe" << endl;

  for (i = 0; i < ns; i++) rhos[i] = rho* cs[i];
  rhoenergy[0] = rho*e; //- 300811.74424616753822
  mix->setState(rhos.data(), rhoenergy.data(), 0);
  Gamma = mix->mixtureFrozenGamma();
  //Gamma = ComputeGamma(rhos);
  Gamma_Minus_One = Gamma - 1;
  Cp = Gamma / Gamma_Minus_One * Gas_Constant;
  Cv = Cp - Gas_Constant;
  Density = rho;
  StaticEnergy = e;
  //cout << endl << "e =" << e << endl;

  mix->getTemperatures(temp.data());
  Temperature = temp[0]; //Gamma_Minus_One * StaticEnergy / Gas_Constant;

  Pressure = mix->pressure(Temperature, rho, cs);//Gamma_Minus_One * Density * StaticEnergy;

  //cout << "SetTDState_rhoe P=" << Pressure << endl;
  ////
  //cout << "SetTDState_rhoe rho=" << rho << endl;
  ////cout << "SetTDState_rhoe Gamma=" << Gamma << endl;
  ////cout << "SetTDState_rhoe e=" << e << endl;
  //cout << "SetTDState_rhoe T=" << Temperature << endl;
  
  SoundSpeed2 = Gamma * Pressure / Density;

  //cout << endl << "SoundSpeed =" << sqrt(SoundSpeed2) << endl;

  dPdrho_e = Gamma_Minus_One * StaticEnergy;
  dPde_rho = Gamma_Minus_One * Density;
  dTdrho_e = 0.0;
  dTde_rho = Gamma_Minus_One / Gas_Constant;

  if (ComputeEntropy) Entropy = (1.0 / Gamma_Minus_One * log(Temperature) + log(1.0 / Density)) * Gas_Constant;
}

void CThermallyPerfectGas::SetTDState_PT(su2double P, su2double T) {
  
  //cout << endl << "SetTDState_rhoe" << endl;


  PT[0] = P;
  PT[1] = T;
  //for (int i = 0; i<ns; i++) cout << "cs[" << i << "]=" << cs[i] <<endl;
  mix->setState(cs, PT.data(), 2);
  
  su2double rho = mix->density(); //P / (T * Gas_Constant);
  for (i = 0; i < ns; i++) rhos[i] = rho* cs[i];

  //cout << endl << "SetTDState_PT P=" << PT[0] << endl;
  //cout << "SetTDState_PT T=" << PT[1] << endl;
  //cout << "SetTDState_PT rho=" << rho << endl;

  Gamma = mix->mixtureFrozenGamma();
  //cout << "SetTDState_PT Gamma=" << Gamma << endl;
  //exit(0);
  //Gamma = ComputeGamma(rhos);
  Gamma_Minus_One = Gamma - 1;
  Cp = Gamma / Gamma_Minus_One * Gas_Constant;
  Cv = Cp - Gas_Constant;

  mix->mixtureEnergies(energy.data());
  su2double e = energy[0];//T * Gas_Constant / Gamma_Minus_One;
  //cout << "SetTDState_PT e=" << e << endl;
  //su2double rho = P / (T * Gas_Constant);
  SetTDState_rhoe(rho, e);
}

// WE CAN IGNORE FOR THERMALLY PERFECT GAS
void CThermallyPerfectGas::SetTDState_Prho(su2double P, su2double rho) {
  su2double e = P / (Gamma_Minus_One * rho);
  SetTDState_rhoe(rho, e);
}

void CThermallyPerfectGas::SetEnergy_Prho(su2double rho, su2double P, su2double T) { 

  //cout << endl << "SetEnergy_Prho" << endl;

  temp[0] = T;
  for (i = 0; i < ns; i++) rhos[i] = rho* cs[i];
  mix->setState(rhos.data(), temp.data(), 1);
  Gamma = mix->mixtureFrozenGamma();
  //Gamma = ComputeGamma(rhos);
  Gamma_Minus_One = Gamma - 1;
  Cp = Gamma / Gamma_Minus_One * Gas_Constant;
  Cv = Cp - Gas_Constant;

  mix->mixtureEnergies(energy.data());
  StaticEnergy = energy[0];//P / (rho * Gamma_Minus_One); 
}

// WE CAN IGNORE FOR THERMALLY PERFECT GAS
void CThermallyPerfectGas::SetTDState_hs(su2double h, su2double s) {
  su2double T = h * Gamma_Minus_One / Gas_Constant / Gamma;
  su2double e = h / Gamma;
  su2double v = exp(-1 / Gamma_Minus_One * log(T) + s / Gas_Constant);

  SetTDState_rhoe(1 / v, e);
}

void CThermallyPerfectGas::SetTDState_Ps(su2double P, su2double s) {
  su2double T = exp(Gamma_Minus_One / Gamma * (s / Gas_Constant + log(P) - log(Gas_Constant)));
  su2double rho = P / (T * Gas_Constant);

  SetTDState_Prho(P, rho);
}

void CThermallyPerfectGas::SetTDState_rhoT(su2double rho, su2double T) {
  su2double e = T * Gas_Constant / Gamma_Minus_One;
  SetTDState_rhoe(rho, e);
}

void CThermallyPerfectGas::ComputeDerivativeNRBC_Prho(su2double P, su2double rho) {
  su2double dPdT_rho, dPdrho_T, dPds_rho;

  SetTDState_Prho(P, rho);

  dPdT_rho = Gas_Constant * rho;
  dPdrho_T = Gas_Constant * Temperature;

  dhdrho_P = -dPdrho_e / dPde_rho - P / rho / rho;
  dhdP_rho = 1.0 / dPde_rho + 1.0 / rho;
  dPds_rho = rho * rho * (SoundSpeed2 - dPdrho_T) / dPdT_rho;
  dsdP_rho = 1.0 / dPds_rho;
  dsdrho_P = -SoundSpeed2 / dPds_rho;
}

//su2double CThermallyPerfectGas::ComputeGamma(vector<su2double> rhos){
//
//  su2double g;
//
//  mix->getCvsMass(Cv_ks.data());
//  for(i = 0; i < ns; i++){
//    Cvtrs[i] = Cv_ks[i];
//    Cvves[i] = Cv_ks[ns+i];
//  }
//
//  /*--- Gamma Computation ---*/
//  su2double rhoR = 0.0, rhoCvtr = 0.0, rhoCvve = 0.0;
//  for(i = 0; i < ns; i++){
//    rhoR += rhos[i]*Ru/MolarMass[i];
//    rhoCvtr += rhos[i]*Cvtrs[i];
//    rhoCvve += rhos[i]*Cvves[i];
//  }
//
//  g = rhoR/(rhoCvtr+rhoCvve)+1;
//
//  return g;
//
//
//}


void CThermallyPerfectGas::SetEnergy_Prho(su2double P, su2double rho) {

  std::cout<<"Wrong energy"<<std::endl; std::exit(1);
 StaticEnergy = P / (rho * Gamma_Minus_One); 

}
