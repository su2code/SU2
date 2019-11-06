/*!
 * transport_model.cpp
 * \brief Source of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna, T. Economon
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.p
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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
 *make ins
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/reactive_structure.hpp"

CReactive::CReactive() {}

CReactive::~CReactive(void) {}



CReactiveMutation::CReactiveMutation(string Mixfile, string Transport) : CReactive(),opt(Mixfile){

  opt.setMechanism(Mixfile);   
  opt.setStateModel("ChemNonEqTTv");
  opt.setViscosityAlgorithm(Transport);
  opt.setThermalConductivityAlgorithm(Transport);

}

CReactiveMutation::~CReactiveMutation(void) {

   delete [] mutation;
   delete [] Ds;
   delete [] hs;
   delete [] comp;
   delete [] Cp_trs;
   delete [] Cp_ves;
   //delete [] Ws;
 }

void CReactiveMutation::InitializeMixture(CConfig *config) {

  ionization = false;
 
  nSpecies = config->GetnSpecies();

  comp     = new su2double[nSpecies];

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) comp[iSpecies]=config->GetInitial_Gas_Composition(iSpecies);

  //std::cout << "MUTATION 1 reactive"  << std::endl;  
    
  mutation = new CMutation(comp, nSpecies, opt); //Cat1

  //std::cout << "MUTATION 2 reactive"  << std::endl;  

  Ms       = mutation->Mutation_MolarMass(); //Cat1cd -


  hs = new su2double[nSpecies];
  Ds = new su2double[nSpecies];

}

vector<su2double> CReactiveMutation::Get_CvTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve) {

  Cv_ks  = mutation->Mutation_Get_CvModeSpecies(cs, rho, T, Tve);

  Cv_trs.resize(nSpecies);

  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) Cv_trs[iSpecies] = Cv_ks[iSpecies];

  return Cv_trs;

}

 vector<su2double> CReactiveMutation::Get_CvVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve) {

  Cv_ks = mutation->Mutation_Get_CvModeSpecies(cs, rho, T, Tve);

  Cv_ves.resize(nSpecies);

  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) Cv_ves[iSpecies] = Cv_ks[nSpecies+iSpecies];

  return Cv_ves;
}

su2double* CReactiveMutation::Get_CpTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve) {

  Cp_ks  = mutation->Mutation_Get_CpModeSpecies(cs, rho, T, Tve);

  Cp_trs = new su2double[nSpecies];

  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) Cp_trs[iSpecies] = Cp_ks[iSpecies];

  return Cp_trs;

}

su2double* CReactiveMutation::Get_CpVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve) {

  Cp_ks = mutation->Mutation_Get_CpModeSpecies(cs, rho, T, Tve);

  Cp_ves = new su2double[nSpecies];

  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) Cp_ves[iSpecies] = Cp_ks[nSpecies+iSpecies];

  return Cp_ves;
}

su2double CReactiveMutation::Get_Gamma(su2double *cs, su2double rho, su2double T, su2double Tve){

  su2double Cv, sum;

  vector<su2double> Cv_trs, Cv_ves, Cvs;

  Cvs.resize(nSpecies);

  //std::cout << setprecision(6)<< "Mutation: 1" << std::endl << std::endl;

  //Cp_trs = Get_CpTraRotSpecies(cs, rho, T, Tve);
  //Cp_ves = Get_CpVibElSpecies(cs, rho, T, Tve);
  Cv_trs = Get_CvTraRotSpecies(cs, rho, T, Tve);
  Cv_ves = Get_CvVibElSpecies(cs, rho, T, Tve);

  //std::cout << setprecision(6)<< "Mutation: 2" << std::endl << std::endl;

    

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    Cvs[iSpecies] = Cv_trs[iSpecies] + Cv_ves[iSpecies];
    Cv            += cs[iSpecies]*Cvs[iSpecies];
    sum           += cs[iSpecies]/Ms[iSpecies];

  }

  //gamma = 1.4; //Cp/Cv;

  gamma = 1 + (8.314*1000/Cv)*sum  ; //Cp/Cv;

  return gamma;
}

su2double CReactiveMutation::Get_GammaFrozen(su2double *cs, su2double rho, su2double T, su2double Tve){


  gammaFrozen = mutation->Mutation_Get_GammaFrozen(cs, rho, T, Tve);

  return gammaFrozen;
}

su2double CReactiveMutation::Get_GammaEquilibrium(su2double *cs, su2double rho, su2double T, su2double Tve){


  gammaEquilibrium = mutation->Mutation_Get_GammaEquilibrium(cs, rho, T, Tve);

  return gammaEquilibrium;
}





su2double CReactiveMutation::Get_MixtureEnergy(su2double* cs, su2double rho, su2double T, su2double Tve) {


  E = mutation->Mutation_Get_MixtureEnergy(cs, rho, T, Tve);

  return E;
}

vector<su2double> CReactiveMutation::Get_MixtureEnergies(su2double* cs, su2double rho, su2double T, su2double Tve) {

  
  Energies = mutation->Mutation_Get_MixtureEnergies(cs, rho, T, Tve);

  return Energies;
}

vector<su2double> CReactiveMutation::Get_SpeciesEnergies(su2double* cs, su2double rho, su2double T, su2double Tve) {

  
  Energies_Species = mutation->Mutation_Get_SpeciesEnergies(cs, rho, T, Tve);

  return Energies_Species;
}


vector<su2double> CReactiveMutation::Get_NetProductionRates(su2double *cs, su2double rho, su2double T, su2double Tve) {
 

  //Ws = new su2double[nSpecies];

  Ws = mutation->Mutation_Get_NetProductionRates(cs, rho, T, Tve);

  return Ws;

}

su2double CReactiveMutation::Get_VTEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve) {

  
  Omega = mutation->Mutation_Get_EnergysourceTerm(cs, rho, T, Tve);
  OmegaVT = Omega[0];

  //std::cout << "OMEGA Total = "  << OmegaVT <<  std::endl << std::endl << std::endl ;

  return OmegaVT;
}

//su2double CReactiveMutation::Get_CVEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve) {

  
  //Omega = mutation->Mutation_Get_EnergysourceTerm(cs, rho, T, Tve);
  //OmegaCV = Omega[1];

  //std::cout << "OMEGA CV = "  << OmegaCV <<  std::endl << std::endl << std::endl ;

  //return OmegaCV;
//}

su2double CReactiveMutation::Get_ReferenceTemperature(su2double *cs, su2double rho, su2double T, su2double Tve) {


  Tref = mutation->Mutation_Get_ReferenceTemperature(cs, rho, T, Tve);

  return Tref;
}
  
vector<su2double> CReactiveMutation::Get_EnthalpiesFormation(su2double *cs, su2double rho, su2double T, su2double Tve) {


  hf = mutation->Mutation_Get_EnthalpiesFormation(cs, rho, T, Tve);

  return hf;

}

su2double* CReactiveMutation::Get_Enthalpies(su2double *cs, su2double rho, su2double T, su2double Tve) {

  //h.resize(nSpecies);

  

  hs = mutation->Mutation_Get_Enthalpies(cs, rho, T, Tve);

  

  return hs;

}

su2double* CReactiveMutation::Get_DiffusionCoeff(su2double *cs, su2double rho, su2double T, su2double Tve) {

   
   
   Ds = mutation->Mutation_Get_DiffusionCoeff(cs, rho, T, Tve);

   return Ds;


}
su2double  CReactiveMutation::Get_Viscosity(su2double *cs, su2double rho, su2double T, su2double Tve) {

   mu = mutation->Mutation_Get_Viscosity(cs, rho, T, Tve);

   return mu;
}

vector<su2double> CReactiveMutation::Get_ThermalConductivity(su2double *cs, su2double rho, su2double T, su2double Tve) {

   
   lambda = mutation->Mutation_Get_ThermalConductivity(cs, rho, T, Tve);

   return lambda;
}

vector<su2double> CReactiveMutation::Get_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve){

  Temp = mutation->Mutation_Get_Temperatures(cs, rho, rhoE, rhoEve);

  return Temp;
}


su2double CReactiveMutation::Get_SoundSpeedFrozen(su2double *cs, su2double rho, su2double T, su2double Tve){

  a = mutation->Mutation_Get_SoundSpeedFrozen(cs, rho, T, Tve);

  return a;

}

su2double CReactiveMutation::Get_Density(su2double T, su2double *Xs, su2double P){

  Density = mutation->Mutation_Get_Density(T, Xs, P);

  return Density;
}



CReactiveHardCode::CReactiveHardCode() : CReactive() {

  std::cout << "Mutation: new CReactiveHardCode"  << std::endl;  
}

CReactiveHardCode::~CReactiveHardCode(void) {}

void CReactiveHardCode::InitializeMixture(CConfig *config) {

  std::cout << "Mutation: CReactiveHardCode: Initialize"  << std::endl;
  
}

vector<su2double> CReactiveHardCode::Get_CvTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve){}

vector<su2double> CReactiveHardCode::Get_CvVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve){}

su2double* CReactiveHardCode::Get_CpTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve){}

su2double* CReactiveHardCode::Get_CpVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve){}

su2double  CReactiveHardCode::Get_Gamma(su2double *cs, su2double rho, su2double T, su2double Tve){}

su2double  CReactiveHardCode::Get_GammaFrozen(su2double *cs, su2double rho, su2double T, su2double Tve){}

su2double  CReactiveHardCode::Get_GammaEquilibrium(su2double *cs, su2double rho, su2double T, su2double Tve){}

su2double  CReactiveHardCode::Get_MixtureEnergy(su2double* cs, su2double rho, su2double T, su2double Tve) {}

vector<su2double> CReactiveHardCode::Get_MixtureEnergies(su2double* cs, su2double rho, su2double T, su2double Tve) {}

vector<su2double> CReactiveHardCode::Get_SpeciesEnergies(su2double* cs, su2double rho, su2double T, su2double Tve) {}

vector<su2double> CReactiveHardCode::Get_NetProductionRates(su2double *cs, su2double rho, su2double T, su2double Tve) {}

su2double CReactiveHardCode::Get_VTEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve){}

//su2double CReactiveHardCode::Get_CVEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve){}

su2double  CReactiveHardCode::Get_ReferenceTemperature(su2double *cs, su2double rho, su2double T, su2double Tve) {}
  
vector<su2double> CReactiveHardCode::Get_EnthalpiesFormation(su2double *cs, su2double rho, su2double T, su2double Tve) {}

su2double* CReactiveHardCode::Get_Enthalpies(su2double *cs, su2double rho, su2double T, su2double Tve) {}

su2double* CReactiveHardCode::Get_DiffusionCoeff(su2double *cs, su2double rho, su2double T, su2double Tve) {}

su2double  CReactiveHardCode::Get_Viscosity(su2double *cs, su2double rho, su2double T, su2double Tve) {}

vector<su2double> CReactiveHardCode::Get_ThermalConductivity(su2double *cs, su2double rho, su2double T, su2double Tve) {}

vector<su2double> CReactiveHardCode::Get_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve){}

su2double  CReactiveHardCode::Get_SoundSpeedFrozen(su2double *cs, su2double rho, su2double T, su2double Tve){}

su2double CReactiveHardCode::Get_Density(su2double T, su2double *Xs, su2double P){}