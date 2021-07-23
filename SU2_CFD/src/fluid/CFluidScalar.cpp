#include <cmath>
#include <math.h>
#include <vector>
#include <numeric>
#include <iostream> 

#include "../../include/fluid/CFluidScalar.hpp"
#include "../../include/fluid/CConstantViscosity.hpp"
#include "../../include/fluid/CSutherland.hpp"
#include "../../include/fluid/CPolynomialViscosity.hpp"
#include "../../include/fluid/CConstantConductivity.hpp"
#include "../../include/fluid/CConstantConductivityRANS.hpp"
#include "../../include/fluid/CConstantPrandtl.hpp"
#include "../../include/fluid/CConstantPrandtlRANS.hpp"
#include "../../include/fluid/CPolynomialConductivity.hpp"
#include "../../include/fluid/CPolynomialConductivityRANS.hpp"

#include "../../include/fluid/CIncIdealGas.hpp"

// CFluidScalar::CFluidScalar(CConfig *config, CFluidModel *fluidModel) : CFluidModel(), fluidModel(fluidModel){
CFluidScalar::CFluidScalar(CConfig *config, su2double value_pressure_operating) : CFluidModel() {
  
/*
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
*/

  unsigned short n_scalars = config->GetNScalarsInit(); 
  config->SetNScalarsInit(n_scalars); 

  n_species_mixture = n_scalars + 1; 

  specificHeat.resize(n_species_mixture); 
  molarMasses.resize(n_species_mixture); 
  massFractions.resize(n_species_mixture);
  moleFractions.resize(n_species_mixture);
  laminarViscosity.resize(n_species_mixture);
  laminarthermalConductivity.resize(n_species_mixture); 

  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    molarMasses.at(iVar) = config->GetMolecular_Weight(iVar); 
    specificHeat.at(iVar) = config->GetSpecific_Heat_Cp(iVar); 
  }

  wilke = false;
  davidson = true;

  Pressure = value_pressure_operating;
  Gas_Constant = config->GetGas_Constant();
  Gamma = 1.0;
  // fluidModel = new CIncIdealGas(config->GetSpecific_Heat_Cp(), config->GetGas_Constant(), Pressure);

  SetLaminarViscosityModel(config); 
  SetThermalConductivityModel(config); 
}

std::vector<su2double> CFluidScalar::massToMoleFractions(su2double* val_scalars){
  su2double mixtureMolarMass = 0.0; 
  su2double mixtureFractionDenumerator = 0.0; 
  
  // change if val_scalars becomes array of scalar su2double* 
  massFractions.at(0) = val_scalars[0];
  massFractions.at(1) = 1 - val_scalars[0]; 

  for(int iVar = 0; iVar < n_species_mixture; iVar++){
    mixtureMolarMass += massFractions[iVar] / molarMasses[iVar]; 
  }

  for(int iVar = 0; iVar < n_species_mixture; iVar++){
    moleFractions.at(iVar) = (massFractions[iVar] / molarMasses[iVar]) / mixtureMolarMass;
  }
  
  return moleFractions; 
}

void CFluidScalar::SetLaminarViscosityModel(const CConfig* config) {
  switch (config->GetKind_ViscosityModel()) {
    case VISCOSITYMODEL::CONSTANT:
    // Build a list of LaminarViscosity pointers to be used in e.g. wilkeViscosity to get the species viscosities. 
      for (int iVar = 0; iVar < n_species_mixture; iVar++){
        LaminarViscosityPointers[iVar] = unique_ptr<CConstantViscosity>(new CConstantViscosity(config->GetMu_Constant(iVar)));
      }
      break;
    case VISCOSITYMODEL::SUTHERLAND:
      for (int iVar = 0; iVar < n_species_mixture; iVar++){
        LaminarViscosityPointers[iVar] = unique_ptr<CSutherland>(new CSutherland(config->GetMu_Ref(iVar), config->GetMu_Temperature_Ref(iVar), config->GetMu_S(iVar)));
      }
      break;
    case VISCOSITYMODEL::POLYNOMIAL:
      for (int iVar = 0; iVar < n_species_mixture; iVar++){
        LaminarViscosityPointers[iVar] = unique_ptr<CPolynomialViscosity<N_POLY_COEFFS>>(
          new CPolynomialViscosity<N_POLY_COEFFS>(config->GetMu_PolyCoeffND()));
      }
      break;
    case VISCOSITYMODEL::FLAMELET:
      /* do nothing. Viscosity is obtained from the table and set in setTDState_T */
      break;
    default:
      SU2_MPI::Error("Viscosity model not available.", CURRENT_FUNCTION);
      break;
  }
}

void CFluidScalar::SetThermalConductivityModel(const CConfig* config) {
  switch (config->GetKind_ConductivityModel()) {
    case CONDUCTIVITYMODEL::CONSTANT:
      for(int iVar = 0; iVar < n_species_mixture; iVar++){
        if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
          ThermalConductivityPointers[iVar] = unique_ptr<CConstantConductivityRANS>(
            new CConstantConductivityRANS(config->GetKt_Constant(iVar), config->GetPrandtl_Turb(iVar)));
        } else {
          ThermalConductivityPointers[iVar] = unique_ptr<CConstantConductivity>(new CConstantConductivity(config->GetKt_Constant(iVar)));
        }
      }  
      break;
    case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
      for(int iVar = 0; iVar < n_species_mixture; iVar++){
        if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
          ThermalConductivityPointers[iVar] = unique_ptr<CConstantPrandtlRANS>(
            new CConstantPrandtlRANS(config->GetPrandtl_Lam(iVar), config->GetPrandtl_Turb(iVar)));
        } else {
          ThermalConductivityPointers[iVar] = unique_ptr<CConstantPrandtl>(new CConstantPrandtl(config->GetPrandtl_Lam(iVar)));
        }
      }
      break;
    case CONDUCTIVITYMODEL::POLYNOMIAL:
      for(int iVar = 0; iVar < n_species_mixture; iVar++){
        if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
          ThermalConductivity = unique_ptr<CPolynomialConductivityRANS<N_POLY_COEFFS>>(
            new CPolynomialConductivityRANS<N_POLY_COEFFS>(config->GetKt_PolyCoeffND(), config->GetPrandtl_Turb()));
        } else {
          ThermalConductivity = unique_ptr<CPolynomialConductivity<N_POLY_COEFFS>>(
            new CPolynomialConductivity<N_POLY_COEFFS>(config->GetKt_PolyCoeffND()));
        }
      }
      break;
    case CONDUCTIVITYMODEL::FLAMELET:
      /* do nothing. Conductivity is obtained from the table and set in setTDState_T */
      break;
    default:
      SU2_MPI::Error("Conductivity model not available.", CURRENT_FUNCTION);
      break;
  }
}

su2double CFluidScalar::wilkeViscosity(su2double* val_scalars){
  std::vector<su2double> phi; 
  std::vector<su2double> wilkeNumerator;
  std::vector<su2double> wilkeDenumeratorSum;
  su2double wilkeDenumerator = 0.0; 
  wilkeDenumeratorSum.clear();
  wilkeNumerator.clear(); 
  su2double viscosityMixture = 0.0;
   
  massToMoleFractions(val_scalars); 

  // Density = GetDensity();
  // Temperature = GetTemperature();

  // here laminarViscosity is filled with n_species_mixture viscosity values. Commented because laminarViscosity is hardcoded right now

  for (int iVar = 0; iVar < n_species_mixture; iVar++){
    LaminarViscosityPointers[iVar]->SetViscosity(Temperature, Density);
    laminarViscosity.at(iVar) = LaminarViscosityPointers[iVar]->GetViscosity();
  }

  for(int i = 0; i < n_species_mixture; i++){
    for(int j = 0; j < n_species_mixture; j++){
      phi.push_back(pow((1 + pow(laminarViscosity[i] / laminarViscosity[j],0.5) * pow(molarMasses[j] / molarMasses[i],0.25)),2) / pow(8 * (1 + molarMasses[i] / molarMasses[j]),0.5)); 
      wilkeDenumerator += moleFractions[j] * phi[j];   
    }
    wilkeDenumeratorSum.push_back(wilkeDenumerator); 
    wilkeDenumerator = 0.0; 
    phi.clear();
    wilkeNumerator.push_back(moleFractions[i] * laminarViscosity[i]); 
    viscosityMixture += wilkeNumerator[i] / wilkeDenumeratorSum[i];  
  }
  return viscosityMixture; 
}

su2double CFluidScalar::davidsonViscosity(su2double* val_scalars){
  su2double viscosityMixture = 0.0;
  su2double fluidity = 0.0;
  su2double E = 0.0; 
  su2double mixtureFractionDenumerator = 0.0;
  const su2double A = 0.375;   
  std::vector<su2double> mixtureFractions; 
  mixtureFractions.clear(); 
  
  massToMoleFractions(val_scalars); 

 /* 
  Density = GetDensity();
  Temperature = GetTemperature();
*/

  for (int iVar = 0; iVar < n_species_mixture; iVar++){
    LaminarViscosityPointers[iVar]->SetViscosity(Temperature, Density);
    laminarViscosity.at(iVar) = LaminarViscosityPointers[iVar]->GetViscosity();
  }

  for(int i = 0; i < n_species_mixture; i++){
    mixtureFractionDenumerator += moleFractions[i] * sqrt(molarMasses[i]);
  }

  for(int j = 0; j < n_species_mixture; j++){
    mixtureFractions.push_back((moleFractions[j] * sqrt(molarMasses[j])) / mixtureFractionDenumerator); 
  }

  for(int i = 0; i < n_species_mixture; i++){
    for(int j = 0; j < n_species_mixture; j++){
      E = (2*sqrt(molarMasses[i]) * sqrt(molarMasses[j])) / (molarMasses[i] + molarMasses[j]);
      fluidity += ((mixtureFractions[i] * mixtureFractions[j]) / (sqrt(laminarViscosity[i]) * sqrt(laminarViscosity[j]))) * pow(E, A);  
    }
  }
  return viscosityMixture = 1 / fluidity;
}

su2double CFluidScalar::wilkeConductivity(su2double* val_scalars){
  std::vector<su2double> phi; 
  std::vector<su2double> wilkeNumerator;
  std::vector<su2double> wilkeDenumeratorSum;
  su2double wilkeDenumerator = 0.0; 
  su2double conductivityMixture = 0.0;
  wilkeDenumeratorSum.clear();
  wilkeNumerator.clear(); 
   
  massToMoleFractions(val_scalars); 

  for (int iVar = 0; iVar < n_species_mixture; iVar++){
    ThermalConductivityPointers[iVar]->SetConductivity(Temperature, Density, Mu, Mu_Turb, Cp);
    laminarthermalConductivity.at(iVar) = ThermalConductivityPointers[iVar]->GetConductivity();
  }

  for(int i = 0; i < n_species_mixture; i++){
    for(int j = 0; j < n_species_mixture; j++){
      phi.push_back(pow((1 + pow((laminarViscosity[i]) / (laminarViscosity[j]),0.5) * pow(molarMasses[j] / molarMasses[i],0.25)),2) / pow(8 * (1 + molarMasses[i] / molarMasses[j]),0.5)); 
      wilkeDenumerator += moleFractions[j] * phi[j];   
    }
    wilkeDenumeratorSum.push_back(wilkeDenumerator); 
    wilkeDenumerator = 0.0; 
    phi.clear();
    wilkeNumerator.push_back(moleFractions[i] * laminarthermalConductivity[i]); 
    conductivityMixture += wilkeNumerator[i] / wilkeDenumeratorSum[i];  
  }
  return conductivityMixture; 
}

unsigned long CFluidScalar::SetTDState_T(su2double val_temperature, su2double *val_scalars){
  // fluidModel->SetTDState_T(val_temperature, val_scalars);  // this is probably not needed as well. Relates to the first comment of this document

  //val_scalars moet eigenlijk een pointer array zijn die de massafracties van alle species bevat.
  //dus voor 3 species: 1e array element species 1, 2e array element species 2, 3e array element species 3 = 1 - Yspecies 1 - Yspecies 2
  // bovenste is niet handig. val_scalars moet een pointer array zijn met de massafracties van alle species dus niet N-1.  
  MeanMolecularWeight = 1/(val_scalars[0]/(molarMasses[0]/1000) + (1-val_scalars[0])/(molarMasses[1]/1000)); 
  // FluidModel->SetMeanMolecularWeight(val_scalars); 
   
  // Cp = specificHeat[0] * val_scalars[0] + specificHeat[1] * (1- val_scalars[0]); 
  // Cp = 1000 * val_scalars[0] + 1005 * (1- val_scalars[0]);  
  Cp = 1009.39; 
  Cv = Cp/1.4; 
  Temperature = val_temperature;
  Density = Pressure / ((Temperature * UNIVERSAL_GAS_CONSTANT) / MeanMolecularWeight);

  if(wilke){
    Mu  = wilkeViscosity(val_scalars);
  }
  else if(davidson){
    Mu = davidsonViscosity(val_scalars);  
  }

  Kt = wilkeConductivity(val_scalars); 

  return 0;
}

