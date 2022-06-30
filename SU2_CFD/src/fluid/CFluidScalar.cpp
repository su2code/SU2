/*!
 * \file CFluidScalar.hpp
 * \brief Defines the multicomponent incompressible Ideal Gas model for mixtures.
 * \author T. Economon, Mark Heimgartner, Cristopher Morales Ubal
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fluid/CFluidScalar.hpp"

#include <math.h>

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>
#include <array>

#include "../../include/fluid/CConstantConductivity.hpp"
#include "../../include/fluid/CConstantConductivityRANS.hpp"
#include "../../include/fluid/CConstantPrandtl.hpp"
#include "../../include/fluid/CConstantPrandtlRANS.hpp"
#include "../../include/fluid/CConstantViscosity.hpp"
#include "../../include/fluid/CPolynomialConductivity.hpp"
#include "../../include/fluid/CPolynomialConductivityRANS.hpp"
#include "../../include/fluid/CPolynomialViscosity.hpp"
#include "../../include/fluid/CSutherland.hpp"

CFluidScalar::CFluidScalar(su2double val_Cp, su2double val_gas_constant, const su2double value_pressure_operating,
                           CConfig* config)
    : CFluidModel() {
  n_species_mixture = config->GetnSpecies() + 1;

  molarMasses.resize(n_species_mixture);
  massFractions.resize(n_species_mixture);
  moleFractions.resize(n_species_mixture);
  laminarViscosity.resize(n_species_mixture);
  laminarThermalConductivity.resize(n_species_mixture);

  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    molarMasses[iVar] = config->GetMolecular_Weight(iVar);
  }

  wilke = (config->GetKind_MixingViscosityModel()==MIXINGVISCOSITYMODEL::WILKE);
  davidson = (config->GetKind_MixingViscosityModel()==MIXINGVISCOSITYMODEL::DAVIDSON);

  Pressure_Thermodynamic = value_pressure_operating;
  Gas_Constant = val_gas_constant;
  Gamma = 1.0;
  Cp = val_Cp;
  Cv = Cp;
  SetLaminarViscosityModel(config);
  SetThermalConductivityModel(config);
}

void CFluidScalar::SetLaminarViscosityModel(const CConfig* config) {
  switch (config->GetKind_ViscosityModel()) {
    case VISCOSITYMODEL::CONSTANT:
      /* Build a list of LaminarViscosity pointers to be used in e.g. wilkeViscosity to get the species viscosities. */
      for (int iVar = 0; iVar < n_species_mixture; iVar++) {
        LaminarViscosityPointers[iVar] =
            std::unique_ptr<CConstantViscosity>(new CConstantViscosity(config->GetMu_Constant(iVar)));
      }
      break;
    case VISCOSITYMODEL::SUTHERLAND:
      for (int iVar = 0; iVar < n_species_mixture; iVar++) {
        LaminarViscosityPointers[iVar] = std::unique_ptr<CSutherland>(
            new CSutherland(config->GetMu_Ref(iVar), config->GetMu_Temperature_Ref(iVar), config->GetMu_S(iVar)));
      }
      break;
    case VISCOSITYMODEL::POLYNOMIAL:
      for (int iVar = 0; iVar < n_species_mixture; iVar++) {
        LaminarViscosityPointers[iVar] = std::unique_ptr<CPolynomialViscosity<N_POLY_COEFFS>>(
            new CPolynomialViscosity<N_POLY_COEFFS>(config->GetMu_PolyCoeffND()));
      }
      break;
    default:
      SU2_MPI::Error("Viscosity model not available.", CURRENT_FUNCTION);
      break;
  }
}

void CFluidScalar::SetThermalConductivityModel(const CConfig* config) {
  switch (config->GetKind_ConductivityModel()) {
    case CONDUCTIVITYMODEL::CONSTANT:
      for (int iVar = 0; iVar < n_species_mixture; iVar++) {
        if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
          ThermalConductivityPointers[iVar] = std::unique_ptr<CConstantConductivityRANS>(new CConstantConductivityRANS(
              config->GetThermal_Conductivity_Constant(iVar), config->GetPrandtl_Turb(iVar)));
        } else {
          ThermalConductivityPointers[iVar] = std::unique_ptr<CConstantConductivity>(
              new CConstantConductivity(config->GetThermal_Conductivity_Constant(iVar)));
        }
      }
      break;
    case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
      for (int iVar = 0; iVar < n_species_mixture; iVar++) {
        if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
          ThermalConductivityPointers[iVar] = std::unique_ptr<CConstantPrandtlRANS>(
              new CConstantPrandtlRANS(config->GetPrandtl_Lam(iVar), config->GetPrandtl_Turb(iVar)));
        } else {
          ThermalConductivityPointers[iVar] =
              std::unique_ptr<CConstantPrandtl>(new CConstantPrandtl(config->GetPrandtl_Lam(iVar)));
        }
      }
      break;
    case CONDUCTIVITYMODEL::POLYNOMIAL:
      for (int iVar = 0; iVar < n_species_mixture; iVar++) {
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
}

std::vector<su2double>& CFluidScalar::MassToMoleFractions(const su2double* val_scalars) {
  su2double mixtureMolarMass{0.0};
  su2double val_scalars_sum{0.0};

  for (int i_scalar = 0; i_scalar < n_species_mixture - 1; i_scalar++) {
    massFractions[i_scalar] = val_scalars[i_scalar];
    val_scalars_sum += val_scalars[i_scalar];
  }
  massFractions[n_species_mixture - 1] = 1 - val_scalars_sum;

  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    mixtureMolarMass += massFractions[iVar] / molarMasses[iVar];
  }

  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    moleFractions[iVar] = (massFractions[iVar] / molarMasses[iVar]) / mixtureMolarMass;
  }

  return moleFractions;
}

su2double CFluidScalar::WilkeViscosity(const su2double* val_scalars) {
  std::vector<su2double> phi;
  std::vector<su2double> wilkeNumerator;
  std::vector<su2double> wilkeDenumeratorSum;
  su2double wilkeDenumerator = 0.0;
  wilkeDenumeratorSum.clear();
  wilkeNumerator.clear();
  su2double viscosityMixture = 0.0;

  /* Fill laminarViscosity with n_species_mixture viscosity values. */
  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    LaminarViscosityPointers[iVar]->SetViscosity(Temperature, Density);
    laminarViscosity[iVar] = LaminarViscosityPointers[iVar]->GetViscosity();
  }

  for (int i = 0; i < n_species_mixture; i++) {
    for (int j = 0; j < n_species_mixture; j++) {
      phi.push_back(
          pow(1 + sqrt(laminarViscosity[i] / laminarViscosity[j]) * pow(molarMasses[j] / molarMasses[i], 0.25), 2) /
          sqrt(8 * (1 + molarMasses[i] / molarMasses[j])));
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

su2double CFluidScalar::DavidsonViscosity(const su2double* val_scalars) {
  su2double viscosityMixture = 0.0;
  su2double fluidity = 0.0;
  su2double E = 0.0;
  su2double mixtureFractionDenumerator = 0.0;
  const su2double A = 0.375;
  std::vector<su2double> mixtureFractions;
  mixtureFractions.clear();

  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    LaminarViscosityPointers[iVar]->SetViscosity(Temperature, Density);
    laminarViscosity[iVar] = LaminarViscosityPointers[iVar]->GetViscosity();
  }

  for (int i = 0; i < n_species_mixture; i++) {
    mixtureFractionDenumerator += moleFractions[i] * sqrt(molarMasses[i]);
  }

  for (int j = 0; j < n_species_mixture; j++) {
    mixtureFractions.push_back((moleFractions[j] * sqrt(molarMasses[j])) / mixtureFractionDenumerator);
  }

  for (int i = 0; i < n_species_mixture; i++) {
    for (int j = 0; j < n_species_mixture; j++) {
      E = (2 * sqrt(molarMasses[i]) * sqrt(molarMasses[j])) / (molarMasses[i] + molarMasses[j]);
      fluidity +=
          ((mixtureFractions[i] * mixtureFractions[j]) / (sqrt(laminarViscosity[i]) * sqrt(laminarViscosity[j]))) *
          pow(E, A);
    }
  }
  return viscosityMixture = 1 / fluidity;
}

su2double CFluidScalar::WilkeConductivity(const su2double* val_scalars) {
  std::vector<su2double> phi;
  std::vector<su2double> wilkeNumerator;
  std::vector<su2double> wilkeDenumeratorSum;
  su2double wilkeDenumerator = 0.0;
  su2double conductivityMixture = 0.0;
  wilkeDenumeratorSum.clear();
  wilkeNumerator.clear();

  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    ThermalConductivityPointers[iVar]->SetConductivity(Temperature, Density, Mu, Mu_Turb, Cp);
    laminarThermalConductivity[iVar] = ThermalConductivityPointers[iVar]->GetConductivity();
  }

  for (int i = 0; i < n_species_mixture; i++) {
    for (int j = 0; j < n_species_mixture; j++) {
      phi.push_back(
          pow(1 + sqrt((laminarViscosity[i]) / (laminarViscosity[j])) * pow(molarMasses[j] / molarMasses[i], 0.25), 2) /
          sqrt(8 * (1 + molarMasses[i] / molarMasses[j])));
      wilkeDenumerator += moleFractions[j] * phi[j];
    }
    wilkeDenumeratorSum.push_back(wilkeDenumerator);
    wilkeDenumerator = 0.0;
    phi.clear();
    wilkeNumerator.push_back(moleFractions[i] * laminarThermalConductivity[i]);
    conductivityMixture += wilkeNumerator[i] / wilkeDenumeratorSum[i];
  }
  return conductivityMixture;
}

void CFluidScalar::SetTDState_T(const su2double val_temperature, const su2double* val_scalars) {
  const su2double MeanMolecularWeight = ComputeMeanMolecularWeight(molarMasses, val_scalars);
  Temperature = val_temperature;
  Density = Pressure_Thermodynamic / ((Temperature * UNIVERSAL_GAS_CONSTANT / MeanMolecularWeight));

  MassToMoleFractions(val_scalars);

  if (wilke) {
    Mu = WilkeViscosity(val_scalars);
  } else if (davidson) {
    Mu = DavidsonViscosity(val_scalars);
  }
}
