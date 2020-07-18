/*!
 * \file CNEMOGas.cpp
 * \brief Source of the nonequilibrium gas model.
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

#include "../../include/fluid/CNEMOGas.hpp"

#include <iomanip> //cat:delete

CNEMOGas::CNEMOGas(const CConfig* config): CFluidModel(){

  nSpecies            = config->GetnSpecies();
  MassFrac_Freestream = config->GetMassFrac_FreeStream();

  MolarMass.resize(nSpecies,0.0);
  MassFrac.resize(nSpecies,0.0);
  MolarFractions.resize(nSpecies,0.0);
  rhos.resize(nSpecies,0.0);
  Cvtrs.resize(nSpecies,0.0);         
  Cvves.resize(nSpecies,0.0);               
  eves.resize(nSpecies,0.0);            
  hs.resize(nSpecies,0.0);            
  ws.resize(nSpecies,0.0);            
  DiffusionCoeff.resize(nSpecies,0.0);
  temperatures.resize(nEnergyEq,0.0);
  energies.resize(nEnergyEq,0.0);  
  ThermalConductivities.resize(nEnergyEq,0.0);

  Kind_TransCoeffModel = config->GetKind_TransCoeffModel();

  frozen = config->GetFrozen();

}

CNEMOGas::~CNEMOGas(){}


void CNEMOGas::SetTDStatePTTv(su2double val_pressure, vector<su2double> val_massfrac, su2double val_temperature, su2double val_temperature_ve){

  su2double denom;

  MassFrac = val_massfrac;                   
  Pressure = val_pressure;                   
  T        = val_temperature;                
  Tve      = val_temperature_ve; 
  
  denom   = 0.0;   

  /*--- Calculate mixture density from supplied primitive quantities ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    denom += MassFrac[iSpecies] * (Ru/MolarMass[iSpecies]) * T;
  for (iSpecies = 0; iSpecies < nEl; iSpecies++)
    denom += MassFrac[nSpecies-1] * (Ru/MolarMass[nSpecies-1]) * Tve;
  Density = Pressure / denom;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
    rhos[iSpecies]     = MassFrac[iSpecies]*Density;
    MassFrac[iSpecies] = rhos[iSpecies]/Density;
  } 
}


su2double CNEMOGas::GetSoundSpeed(){

  su2double conc, rhoCvtr;

  conc    = 0.0;
  rhoCvtr = 0.0; 
  Density = 0.0;

  vector<su2double> Cvtrs = GetSpeciesCvTraRot();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Density+=rhos[iSpecies];

  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++){
    conc += rhos[iSpecies]/MolarMass[iSpecies];
    rhoCvtr += rhos[iSpecies] * Cvtrs[iSpecies];
  }
  SoundSpeed2 = (1.0 + Ru/rhoCvtr*conc) * Pressure/Density;

 // cout <<setprecision(10)<< "cat: Ru=" << Ru << endl;
 //cout <<setprecision(10)<< "cat: conc=" << conc << endl;

// cout<<endl<<endl;
//
// cout <<setprecision(10)<< "cat: rhoCvtr=" << rhoCvtr << endl;
// cout <<setprecision(10)<< "cat: Pressure=" << Pressure << endl;
// cout <<setprecision(10)<< "cat: Density=" << Density << endl;

  return(sqrt(SoundSpeed2));

}

su2double CNEMOGas::GetPressure(){

  su2double P = 0.0;

 //for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
 //  cout <<setprecision(20)<< "cat: rhos=" << rhos[iSpecies] << endl;
 //  cout <<setprecision(20)<< "cat: MolarMass=" << MolarMass[iSpecies] << endl;
 //}

 //cout <<setprecision(20)<< "cat: Ru=" << Ru << endl;
 //cout <<setprecision(20)<< "cat: T=" << T << endl;
 //cout <<setprecision(20)<< "cat: Tve=" << Tve << endl;

  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    P += rhos[iSpecies] * Ru/MolarMass[iSpecies] * T;
  for (iSpecies = 0; iSpecies < nEl; iSpecies++)
    P += rhos[nSpecies-1] * Ru/MolarMass[nSpecies-1] * Tve;

  Pressure = P;

  return P;

}

su2double CNEMOGas::GetGasConstant(){

  su2double Mass = 0.0;

  // This needs work for Ionization and such
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    Mass += MassFrac[iSpecies] * MolarMass[iSpecies];
  GasConstant = Ru / Mass;
 
  return GasConstant;
}

su2double CNEMOGas::GetrhoCvve() {

    Cvves = GetSpeciesCvVibEle();

    rhoCvve = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      rhoCvve += rhos[iSpecies]*Cvves[iSpecies];

    return rhoCvve;
}


