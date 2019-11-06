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
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
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

#include "../include/cmutation.hpp"

CMutation::CMutation(su2double *composition, su2double nSpecies, Mutation::MixtureOptions options ):
	 
	mix(options){

  //std::cout << "MUTATION 1 Cmutation"  << std::endl;  

  ns   = nSpecies;

  hs = new su2double[ns];
  Ds = new su2double[ns];
  Xs = new su2double[ns];


  
}

CMutation::~CMutation(void) { 

  
  delete [] Ds;
  delete [] hs;
  delete [] Xs;
  delete [] Cp_ks;
  //delete [] Ws;
}

vector<su2double> CMutation::Mutation_MolarMass(){

  
  
  //MolarMass = new su2double[ns];

  MolarMass.resize(ns);

  for(i = 0; i < ns; i++) MolarMass[i] = 1000* mix.speciesMw(i); // x1000 to have Molar Mass in kg/kmol

  
  return MolarMass;

}

vector<su2double> CMutation::Mutation_SpeciesDensity(su2double *cs, su2double rho){
 
   rhos.resize(ns);

   for(i = 0; i < ns; i++) rhos[i] = cs[i]*rho;

   return rhos; 

}

void CMutation::Mutation_UpdateMixtureState(su2double* cs, su2double rho, su2double T, su2double Tve) {

   rhos  = Mutation_SpeciesDensity(cs, rho);

   temperatures.resize(2);

   temperatures[0] = T;
   temperatures[1] = Tve; 

   mix.setState(rhos.data(), temperatures.data(), 1);

}

vector<su2double> CMutation::Mutation_Get_CvModeSpecies(su2double* cs, su2double rho, su2double T, su2double Tve){

   Cv_ks.resize(2*ns);

   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   mix.getCvsMass(Cv_ks.data());

   return Cv_ks; 

}

su2double* CMutation::Mutation_Get_CpModeSpecies(su2double* cs, su2double rho, su2double T, su2double Tve){

   Cp_ks = new su2double[ns*2];

   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   mix.getCpsMass(Cp_ks);

   //su2double CpMassFrozen, CpMoleFrozen, CpMassEq, CpMoleEq;

   //CpMassFrozen = mix.mixtureFrozenCpMass();

   //CpMoleFrozen = mix.mixtureFrozenCpMole();

   //CpMassEq = mix.mixtureEquilibriumCpMass();

   //CpMoleEq = mix.mixtureEquilibriumCpMole();

   //std::cout << setprecision(6)<< "Mutation: CpMassFrozen=" << CpMassFrozen <<  std::endl;
   //std::cout << setprecision(6)<< "Mutation: CpMoleFrozen=" << CpMoleFrozen <<  std::endl;
   //std::cout << setprecision(6)<< "Mutation: CpMassEq=" << CpMassEq <<  std::endl;
   //std::cout << setprecision(6)<< "Mutation: CpMoleEq=" << CpMoleEq <<  std::endl;

   return Cp_ks; 

}

su2double CMutation::Mutation_Get_GammaFrozen(su2double* cs, su2double rho, su2double T, su2double Tve){
   

   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   gammaFrozen = mix.mixtureFrozenGamma();

   return gammaFrozen; 

}

su2double CMutation::Mutation_Get_GammaEquilibrium(su2double* cs, su2double rho, su2double T, su2double Tve){
   

   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   gammaEquilibrium = mix.mixtureEquilibriumGamma();

   return gammaEquilibrium; 

}

su2double CMutation::Mutation_Get_MixtureEnergy(su2double* cs, su2double rho, su2double T, su2double Tve) {
 
   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   E = mix.mixtureEnergyMass();

   return E; 

}

vector<su2double> CMutation::Mutation_Get_MixtureEnergies(su2double* cs, su2double rho, su2double T, su2double Tve) {
 
   Energies.resize(2);

   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   mix.mixtureEnergies(Energies.data());

   return Energies; 

}

vector<su2double> CMutation::Mutation_Get_SpeciesEnergies(su2double* cs, su2double rho, su2double T, su2double Tve) {
 
   Energies_species.resize(2*ns);

   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   mix.getEnergiesMass(Energies_species.data());

   return Energies_species; 

}


vector<su2double> CMutation::Mutation_Get_NetProductionRates(su2double* cs, su2double rho, su2double T, su2double Tve){
   
   //Ws = new su2double[ns];

  Ws.resize(ns);

   //int nReactions;

   Mutation_UpdateMixtureState(cs, rho, T, Tve);


   //nReactions = mix.nReactions();
   
   
   //su2double *kf = new su2double[nReactions];
   //su2double *kb = new su2double[nReactions];

   //mix.forwardRateCoefficients(kf);
   //mix.backwardRateCoefficients(kb);

   //std::cout << "----- MUTATION -----"<<  std::endl ;

   //for (i = 0; i < nReactions; i++) std::cout << "kf=" << kf[i] << " kb=" << kb[i] << std::endl ;

   mix.netProductionRates(Ws.data());

   //std::cout << "W[N2]=" << Ws[0] <<  std::endl ;
   //std::cout << "W[N]="  << Ws[1] <<  std::endl <<  std::endl;

   


   return Ws; 

}

vector<su2double> CMutation::Mutation_Get_EnergysourceTerm(su2double* cs, su2double rho, su2double T, su2double Tve) {

   Omega.resize(1);
   
   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   mix.energyTransferSource(Omega.data());

   return Omega;
}


su2double CMutation::Mutation_Get_ReferenceTemperature(su2double* cs, su2double rho, su2double T, su2double Tve) {

   
   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   Tref = mix.standardStateT();

   return Tref;
}

vector<su2double> CMutation::Mutation_Get_EnthalpiesFormation(su2double* cs, su2double rho, su2double T, su2double Tve) {

   

   su2double RuSI = UNIVERSAL_GAS_CONSTANT;
   su2double Ru   = 1000.0*RuSI;

   hf.resize(ns);

   vector<su2double> hf_RT; hf_RT.resize(ns);

   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   su2double T_ref = mix.standardStateT();

   mix.speciesHOverRT(T_ref, T_ref, T_ref, T_ref, T_ref, NULL, NULL, NULL, NULL, NULL, hf_RT.data());

   for (i = 0; i < ns; i++) hf[i] = hf_RT[i]*(RuSI*T_ref);

   return hf;
}

su2double* CMutation::Mutation_Get_Enthalpies(su2double* cs, su2double rho, su2double T, su2double Tve) {

   su2double RuSI = UNIVERSAL_GAS_CONSTANT;
   su2double Ru   = 1000.0*RuSI;

   //h.resize(ns);

   //std::cout << "calc mut 1"  << std::endl<< std::endl<< std::endl<< std::endl;


   vector<su2double> h_RT; h_RT.resize(ns);

   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   //std::cout << "calc mut 2"  << std::endl<< std::endl<< std::endl<< std::endl;

   

   mix.speciesHOverRT(T, Tve, T, Tve, Tve, h_RT.data(), NULL, NULL, NULL, NULL, NULL);

   //std::cout << "calc mut 3"  << std::endl<< std::endl<< std::endl<< std::endl;

   for (i = 0; i < ns; i++) hs[i] = h_RT[i]*(RuSI*T); 


   return hs;
}

su2double* CMutation::Mutation_Get_DiffusionCoeff(su2double *cs, su2double rho, su2double T, su2double Tve) {

   //Ds.resize(ns);

   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   mix.averageDiffusionCoeffs(Ds);

   return Ds;


}

su2double  CMutation::Mutation_Get_Viscosity(su2double *cs, su2double rho, su2double T, su2double Tve) {

   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   mu = mix.viscosity();

   return mu;
}

vector<su2double> CMutation::Mutation_Get_ThermalConductivity(su2double *cs, su2double rho, su2double T, su2double Tve) {

   lambda.resize(2);

   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   mix.frozenThermalConductivityVector(lambda.data());

   return lambda;
}

vector<su2double> CMutation::Mutation_Get_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve){

   rhos  = Mutation_SpeciesDensity(cs, rho);

   vector<su2double> en; en.resize(2);

   en[0] = rhoE;
   en[1] = rhoEve; 

   mix.setState(rhos.data(), en.data(), 0);

   Temp.resize(2);

   mix.getTemperatures(Temp.data());

   return Temp;

}

su2double  CMutation::Mutation_Get_SoundSpeedFrozen(su2double *cs, su2double rho, su2double T, su2double Tve){

   Mutation_UpdateMixtureState(cs, rho, T, Tve);

   a = mix.frozenSoundSpeed();

   return a;

}

su2double  CMutation::Mutation_Get_Density(su2double T, su2double *Xs, su2double P){

   
   Density = mix.density(T, P, Xs);

   return Density;

}





