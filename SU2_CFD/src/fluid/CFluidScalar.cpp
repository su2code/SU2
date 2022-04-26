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


CFluidScalar::CFluidScalar(su2double val_Cp, su2double val_gas_constant, const su2double value_pressure_operating, CConfig *config) : CFluidModel() {
  // nijso TODO BUG FIXME
  n_scalars = config->GetnSpecies();
  //config->SetNScalarsInit(n_scalars);
  n_species_mixture = n_scalars + 1;

  specificHeat.resize(n_species_mixture);
  molarMasses.resize(n_species_mixture);
  massFractions.resize(n_species_mixture);
  moleFractions.resize(n_species_mixture);
  laminarViscosity.resize(n_species_mixture);
  laminarthermalConductivity.resize(n_species_mixture);

  // nijso FIXME TODO BUG
  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    molarMasses.at(iVar) = config->GetMolecular_Weight(iVar);
    //specificHeat.at(iVar) = config->GetSpecific_Heat_Cp(iVar);
    //molarMasses.at(iVar) = 28.96;    //config->GetMolecular_Weight(iVar);
    specificHeat.at(iVar) = 1000.0 ; //config->GetSpecific_Heat_Cp(iVar);
  }

  wilke = false;
  davidson = true;

  Pressure_Thermodynamic = value_pressure_operating;
  Gas_Constant = val_gas_constant; //config->GetGas_Constant();
  Gamma = 1.0;
  Cp = val_Cp;
  Cv = Cp;
  SetLaminarViscosityModel(config);
  SetThermalConductivityModel(config);
}

void CFluidScalar::SetLaminarViscosityModel(const CConfig* config) {
  switch (config->GetKind_ViscosityModel()) {
    case VISCOSITYMODEL::CONSTANT:
      LaminarViscosity = unique_ptr<CConstantViscosity>(new CConstantViscosity(config->GetMu_ConstantND()));
      break;
    case VISCOSITYMODEL::SUTHERLAND:
      LaminarViscosity = unique_ptr<CSutherland>(
          new CSutherland(config->GetMu_RefND(), config->GetMu_Temperature_RefND(), config->GetMu_SND()));
      break;
    case VISCOSITYMODEL::POLYNOMIAL:
      LaminarViscosity = unique_ptr<CPolynomialViscosity<N_POLY_COEFFS>>(
          new CPolynomialViscosity<N_POLY_COEFFS>(config->GetMu_PolyCoeffND()));
      break;
  }
}

void CFluidScalar::SetThermalConductivityModel(const CConfig* config) {
  switch (config->GetKind_ConductivityModel()) {
    case CONDUCTIVITYMODEL::CONSTANT:
      if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
        ThermalConductivity = unique_ptr<CConstantConductivityRANS>(
            new CConstantConductivityRANS(config->GetThermal_Conductivity_ConstantND(), config->GetPrandtl_Turb()));
      } else {
        ThermalConductivity = unique_ptr<CConstantConductivity>(new CConstantConductivity(config->GetThermal_Conductivity_ConstantND()));
      }
      break;
    case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
      if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
        ThermalConductivity = unique_ptr<CConstantPrandtlRANS>(
            new CConstantPrandtlRANS(config->GetPrandtl_Lam(), config->GetPrandtl_Turb()));
      } else {
        ThermalConductivity = unique_ptr<CConstantPrandtl>(new CConstantPrandtl(config->GetPrandtl_Lam()));
      }
      break;
    case CONDUCTIVITYMODEL::POLYNOMIAL:
      if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
        ThermalConductivity = unique_ptr<CPolynomialConductivityRANS<N_POLY_COEFFS>>(
            new CPolynomialConductivityRANS<N_POLY_COEFFS>(config->GetKt_PolyCoeffND(), config->GetPrandtl_Turb()));
      } else {
        ThermalConductivity = unique_ptr<CPolynomialConductivity<N_POLY_COEFFS>>(
            new CPolynomialConductivity<N_POLY_COEFFS>(config->GetKt_PolyCoeffND()));
      }
      break;
    default:
      SU2_MPI::Error("Conductivity model not available.", CURRENT_FUNCTION);
      break;
  }
}
/*
void CFluidScalar::SetLaminarViscosityModel(const CConfig* config) {
  switch (config->GetKind_ViscosityModel()) {
    case VISCOSITYMODEL::CONSTANT:
      for (int iVar = 0; iVar < n_species_mixture; iVar++){
        LaminarViscosityPointers[iVar] = std::unique_ptr<CConstantViscosity>(new CConstantViscosity(config->GetMu_Constant(iVar)));
      }
      break;
    case VISCOSITYMODEL::SUTHERLAND:
      for (int iVar = 0; iVar < n_species_mixture; iVar++){
        LaminarViscosityPointers[iVar] = std::unique_ptr<CSutherland>(new CSutherland(config->GetMu_Ref(iVar), config->GetMu_Temperature_Ref(iVar), config->GetMu_S(iVar)));
      }
      break;
    case VISCOSITYMODEL::POLYNOMIAL:
      for (int iVar = 0; iVar < n_species_mixture; iVar++){
        LaminarViscosityPointers[iVar] = std::unique_ptr<CPolynomialViscosity<N_POLY_COEFFS>>(
          new CPolynomialViscosity<N_POLY_COEFFS>(config->GetMu_PolyCoeffND()));
      }
      break;
    default:
      SU2_MPI::Error("Viscosity model not available.", CURRENT_FUNCTION);
      break;
  }
}*/

/*
void CFluidScalar::SetThermalConductivityModel(const CConfig* config) {
  switch (config->GetKind_ConductivityModel()) {
    case CONDUCTIVITYMODEL::CONSTANT:
      for(int iVar = 0; iVar < n_species_mixture; iVar++){
        if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
          ThermalConductivityPointers[iVar] = std::unique_ptr<CConstantConductivityRANS>(
            new CConstantConductivityRANS(config->GetThermal_Conductivity_Constant(iVar), config->GetPrandtl_Turb(iVar)));
        } else {
          ThermalConductivityPointers[iVar] = std::unique_ptr<CConstantConductivity>(new CConstantConductivity(config->GetThermal_Conductivity_Constant(iVar)));
        }
      }
      break;
    case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
      for(int iVar = 0; iVar < n_species_mixture; iVar++){
        if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
          ThermalConductivityPointers[iVar] = std::unique_ptr<CConstantPrandtlRANS>(
            new CConstantPrandtlRANS(config->GetPrandtl_Lam(iVar), config->GetPrandtl_Turb(iVar)));
        } else {
          ThermalConductivityPointers[iVar] = std::unique_ptr<CConstantPrandtl>(new CConstantPrandtl(config->GetPrandtl_Lam(iVar)));
        }
      }
      break;
    case CONDUCTIVITYMODEL::POLYNOMIAL:
      for(int iVar = 0; iVar < n_species_mixture; iVar++){
        if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
          ThermalConductivity = std::unique_ptr<CPolynomialConductivityRANS<N_POLY_COEFFS>>(
            new CPolynomialConductivityRANS<N_POLY_COEFFS>(config->GetKt_PolyCoeffND(), config->GetPrandtl_Turb()));
        } else {
          ThermalConductivity = std::unique_ptr<CPolynomialConductivity<N_POLY_COEFFS>>(
            new CPolynomialConductivity<N_POLY_COEFFS>(config->GetKt_PolyCoeffND()));
        }
      }
      break;
    default:
      SU2_MPI::Error("Conductivity model not available.", CURRENT_FUNCTION);
      break;
  }
}*/

/*std::vector<su2double>& CFluidScalar::massToMoleFractions(const su2double * const val_scalars){
  su2double mixtureMolarMass {0.0};
  su2double val_scalars_sum {0.0};

  for(int i_scalar = 0; i_scalar < n_scalars; i_scalar++){
    massFractions[i_scalar] = val_scalars[i_scalar];
    val_scalars_sum += val_scalars[i_scalar];
  }
  massFractions[n_scalars] = 1 - val_scalars_sum;

  for(int iVar = 0; iVar < n_species_mixture; iVar++){
    mixtureMolarMass += massFractions[iVar] / molarMasses[iVar];
  }

  for(int iVar = 0; iVar < n_species_mixture; iVar++){
    moleFractions.at(iVar) = (massFractions[iVar] / molarMasses[iVar]) / mixtureMolarMass;
  }

  return moleFractions;
}*/

/*
su2double CFluidScalar::wilkeViscosity(const su2double * const val_scalars){
  std::vector<su2double> phi;
  std::vector<su2double> wilkeNumerator;
  std::vector<su2double> wilkeDenumeratorSum;
  su2double wilkeDenumerator = 0.0;
  wilkeDenumeratorSum.clear();
  wilkeNumerator.clear();
  su2double viscosityMixture = 0.0;


  for (int iVar = 0; iVar < n_species_mixture; iVar++){
    LaminarViscosityPointers[iVar]->SetViscosity(Temperature, Density);
    laminarViscosity.at(iVar) = LaminarViscosityPointers[iVar]->GetViscosity();
  }

  for(int i = 0; i < n_species_mixture; i++){
    for(int j = 0; j < n_species_mixture; j++){
      phi.push_back(pow(1 + sqrt(laminarViscosity[i] / laminarViscosity[j]) * pow(molarMasses[j] / molarMasses[i], 0.25), 2) / sqrt(8 * (1 + molarMasses[i] / molarMasses[j])));
      wilkeDenumerator += moleFractions[j] * phi[j];
    }
    wilkeDenumeratorSum.push_back(wilkeDenumerator);
    wilkeDenumerator = 0.0;
    phi.clear();
    wilkeNumerator.push_back(moleFractions[i] * laminarViscosity[i]);
    viscosityMixture += wilkeNumerator[i] / wilkeDenumeratorSum[i];
  }
  return viscosityMixture;
}*/

/*
su2double CFluidScalar::davidsonViscosity(const su2double * const val_scalars){
  su2double viscosityMixture = 0.0;
  su2double fluidity = 0.0;
  su2double E = 0.0;
  su2double mixtureFractionDenumerator = 0.0;
  const su2double A = 0.375;
  std::vector<su2double> mixtureFractions;
  mixtureFractions.clear();

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
*/

/*
su2double CFluidScalar::wilkeConductivity(const su2double * const val_scalars){
  std::vector<su2double> phi;
  std::vector<su2double> wilkeNumerator;
  std::vector<su2double> wilkeDenumeratorSum;
  su2double wilkeDenumerator = 0.0;
  su2double conductivityMixture = 0.0;
  wilkeDenumeratorSum.clear();
  wilkeNumerator.clear();

  for (int iVar = 0; iVar < n_species_mixture; iVar++){
    ThermalConductivityPointers[iVar]->SetConductivity(Temperature, Density, Mu, Mu_Turb, Cp);
    laminarthermalConductivity.at(iVar) = ThermalConductivityPointers[iVar]->GetConductivity();
  }

  for(int i = 0; i < n_species_mixture; i++){
    for(int j = 0; j < n_species_mixture; j++){
      phi.push_back(pow(1 + sqrt((laminarViscosity[i]) / (laminarViscosity[j])) * pow(molarMasses[j] / molarMasses[i], 0.25), 2) / sqrt(8 * (1 + molarMasses[i] / molarMasses[j])));
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
*/

void CFluidScalar::SetTDState_T(const su2double val_temperature, su2double * const val_scalars){
  const su2double MeanMolecularWeight = ComputeMeanMolecularWeight(molarMasses, val_scalars);

  // CFluidModel::ComputeMeanSpecificHeatCp(specificHeat, val_scalars);
  //const su2double CpAir300Kelvin = 1009.39;
  //const su2double RatioSpecificHeatsAir = 1.4;
  //Cp = CpAir300Kelvin;
  //Cv = Cp/RatioSpecificHeatsAir;
  Temperature = val_temperature;
  //Density = Pressure_Thermodynamic / (Temperature * Gas_Constant);
  Density = Pressure_Thermodynamic / ((Temperature * UNIVERSAL_GAS_CONSTANT*MeanMolecularWeight)); // /MeanMolecularWeight

  //massToMoleFractions(val_scalars);


  //if(wilke){
    //Mu  = wilkeViscosity(val_scalars);
  //}
  //else if(davidson){
    //Mu = davidsonViscosity(val_scalars);
  //}

  //Kt = wilkeConductivity(val_scalars);
  

  //return 0;
}
