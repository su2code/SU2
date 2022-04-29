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
  
  n_scalars = config->GetnSpecies();
  n_species_mixture = n_scalars + 1;

  specificHeat.resize(n_species_mixture);
  molarMasses.resize(n_species_mixture);
  massFractions.resize(n_species_mixture);
  moleFractions.resize(n_species_mixture);
  laminarViscosity.resize(n_species_mixture);
  laminarthermalConductivity.resize(n_species_mixture);

  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    molarMasses.at(iVar) = config->GetMolecular_Weight(iVar);
    specificHeat.at(iVar) = 1004.703; //config->GetSpecific_Heat_Cp(iVar);
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

void CFluidScalar::SetTDState_T(const su2double val_temperature, su2double * const val_scalars){
  const su2double MeanMolecularWeight = ComputeMeanMolecularWeight(molarMasses, val_scalars);
  Temperature = val_temperature;
  Density = Pressure_Thermodynamic / ((Temperature * UNIVERSAL_GAS_CONSTANT*MeanMolecularWeight)); 
}
