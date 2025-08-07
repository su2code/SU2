/*!
 * \file CFluidCantera.cpp
 * \brief Defines the multicomponent incompressible Ideal Gas model for reacting flows.
 * \author T. Economon, Cristopher Morales Ubal
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fluid/CFluidCantera.hpp"

#include <cmath>

#include <numeric>
#ifdef USE_CANTERA
#include "../../Common/include/basic_types/ad_structure.hpp"

#include "../../../subprojects/cantera/include/cantera/core.h"
#include "../../../subprojects/cantera/include/cantera/kinetics/Reaction.h"
#include "../../../subprojects/cantera/include/cantera/zerodim.h"
#include <fstream>
#include <iostream>

using namespace Cantera;
using namespace SU2_TYPE;

CFluidCantera::CFluidCantera(su2double value_pressure_operating, const CConfig* config)
    : CFluidModel(),
      n_species_mixture(config->GetnSpecies() + 1),
      Pressure_Thermodynamic(value_pressure_operating),
      GasConstant_Ref(config->GetGas_Constant_Ref()),
      Prandtl_Number(config->GetPrandtl_Turb()),
      Transport_Model(config->GetTransport_Model()),
      Chemical_MechanismFile(config->GetChemical_MechanismFile()),
      Phase_Name(config->GetPhase_Name()),
      Chemistry_Time_Integration(config->GetChemistryTimeIntegration()){
  if (n_species_mixture > ARRAYSIZE) {
    SU2_MPI::Error("Too many species, increase ARRAYSIZE", CURRENT_FUNCTION);
  }
  sol = std::shared_ptr<Cantera::Solution>(newSolution(Chemical_MechanismFile, Phase_Name, Transport_Model));
  sol->thermo()->getMolecularWeights(&molarMasses[0]);
  combustor = nullptr;
  sim = nullptr;
  if (Chemistry_Time_Integration) {
    combustor = std::static_pointer_cast<Cantera::IdealGasConstPressureReactor>(newReactor("IdealGasConstPressureReactor", sol));
    std::vector<std::shared_ptr<Cantera::ReactorBase>> reactor = { combustor };
    sim = std::shared_ptr<Cantera::ReactorNet>(newReactorNet(reactor));
    su2double Delta_t_max = 0.1;
    combustor->setAdvanceLimit("temperature", GetValue(Delta_t_max));
  }
  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    gasComposition[iVar] = config->GetChemical_GasComposition(iVar);
    // config->SetChemical_GasComposition(iVar, gasComposition[iVar]); //this should be used for later
  }
  SetEnthalpyFormation(config);

  SetMassDiffusivityModel(config);
}

void CFluidCantera::SetMassDiffusivityModel(const CConfig* config) {
  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    MassDiffusivityPointers[iVar] = MakeMassDiffusivityModel(config, iVar);
  }
}

void CFluidCantera::SetEnthalpyFormation(const CConfig* config) {
  /*--- Set mass fractions. ---*/
  const int nsp = sol->thermo()->nSpecies();
  su2double val_scalars_sum{0.0};
  su2double massFractions[nsp]{0.0};
  for (int i_scalar = 0; i_scalar < n_species_mixture - 1; i_scalar++) {
    int speciesIndex = sol->thermo()->speciesIndex(gasComposition[i_scalar]);
    massFractions[speciesIndex] =config->GetSpecies_Init()[i_scalar];
    val_scalars_sum += config->GetSpecies_Init()[i_scalar];
  }
  massFractions[n_species_mixture - 1] = 1.0 - val_scalars_sum;
  sol->thermo()->setMassFractions(massFractions);
  su2double T_ref = 298.15;
  sol->thermo()->setState_TP(T_ref, Pressure_Thermodynamic);
  // The universal gas constant times temperature is retrieved from cantera.
  const su2double uni_gas_constant_temp = sol->thermo()->RT();
  vector<su2double> enthalpiesSpecies(nsp);
  sol->thermo()->getEnthalpy_RT_ref(&enthalpiesSpecies[0]);
  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    int speciesIndex = sol->thermo()->speciesIndex(gasComposition[iVar]);
    enthalpyFormation[iVar] = uni_gas_constant_temp * enthalpiesSpecies[speciesIndex] / molarMasses[speciesIndex];
  }
}

void CFluidCantera::ComputeMassDiffusivity() {
  int nsp = sol->thermo()->nSpecies();
  vector<su2double> diff(nsp);
  sol->transport()->getMixDiffCoeffsMass(&diff[0]);
  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    massDiffusivity[iVar] = diff[sol->thermo()->speciesIndex(gasComposition[iVar])];
  }
}

void CFluidCantera::ComputeChemicalSourceTerm(su2double delta_time, const su2double* val_scalars){
  if (Chemistry_Time_Integration) {
    combustor->syncState();
    //combustor->initialize(); //it might be needed in the future
    su2double Delta_t_max = 0.1;
    combustor->setAdvanceLimit("temperature", GetValue(Delta_t_max));
    sim->setInitialTime(0.0);
    sim->advance(max(1E-15, GetValue(delta_time)));
    sim->setTolerances(1E-12, 1E-12);
    //const su2double density_new = combustor->density();
    // cout << "time simulated " << setprecision(12) << sim->time() << endl;
    // cout << "Temperature before integration: " << setprecision(12) << Temperature << endl;
    // cout << "Temperature after integration: " << setprecision(12) << combustor->temperature() << endl;
    // cout << "Pressure combustor after integration: " << setprecision(12) << combustor->pressure() << endl;
    for (int iVar = 0; iVar < n_species_mixture - 1.0; iVar++) {
      int speciesIndex = sol->thermo()->speciesIndex(gasComposition[iVar]);
      const su2double scalar_new = combustor->massFraction(speciesIndex);
      //const su2double scalat_old = val_scalars[iVar];
      // cout << "Species before integration: " << setprecision(12) << val_scalars[iVar] << endl;
      // cout << "Species after integration: " << setprecision(12) << scalar_new << endl;
      //cout << "Time simulated: " << setprecision(12) << sim->time() << endl;
      const su2double source_term_corr = Density * (scalar_new - val_scalars[iVar]) / abs(sim->time());
      // cout<<"source_term "<<source_term_corr<<endl;
      chemicalSourceTerm[iVar] = source_term_corr;  // molarMasses[speciesIndex]*netProductionRates[speciesIndex];
    }
  } else {
    const int nsp = sol->thermo()->nSpecies();
    vector<su2double> netProductionRates(nsp);
    sol->kinetics()->getNetProductionRates(&netProductionRates[0]);
    for (int iVar = 0; iVar < n_species_mixture - 1.0; iVar++) {
      int speciesIndex = sol->thermo()->speciesIndex(gasComposition[iVar]);
      chemicalSourceTerm[iVar] = molarMasses[speciesIndex] * netProductionRates[speciesIndex];
    }
  }
}

void CFluidCantera::ComputeGradChemicalSourceTerm(const su2double* val_scalars) {
  const int nsp = sol->thermo()->nSpecies();
  //const su2double meanMolecularWeight = sol->thermo()->meanMolecularWeight();
  vector<su2double> netProductionRates_T(nsp);
  // The universal gas constant times temperature is retrieved from cantera.
  const su2double uni_gas_constant_temp = sol->thermo()->RT();
  vector<su2double> enthalpiesSpecies(nsp);
  sol->thermo()->getEnthalpy_RT_ref(&enthalpiesSpecies[0]);
  //Eigen::SparseMatrix<su2double> dW_dC = sol->kinetics()->netProductionRates_ddCi();
  sol->kinetics()->getNetProductionRates_ddT(&netProductionRates_T[0]);
  const int speciesN = sol->thermo()->speciesIndex(gasComposition[n_species_mixture - 1]);
  for (int iVar = 0; iVar < n_species_mixture - 1; iVar++) {
    int speciesIndex_i = sol->thermo()->speciesIndex(gasComposition[iVar]);
    //su2double scalars_sum = 0.0;
    const su2double dT_dYi =
        uni_gas_constant_temp * (enthalpiesSpecies[speciesN] - enthalpiesSpecies[speciesIndex_i]) / Cp;
    gradChemicalSourceTerm[iVar] = molarMasses[speciesIndex_i] * netProductionRates_T[speciesIndex_i] * dT_dYi;
    // for (int jVar = 0; jVar < n_species_mixture - 1; jVar++) {
    //   int speciesIndex_j = sol->thermo()->speciesIndex(gasComposition[jVar]);
    //   su2double dW_dCi = dW_dC.coeff(iVar,jVar) * molarMasses[speciesIndex_i];
    //   if (jVar != iVar) {
    //     const su2double factor = dW_dCi * Density * val_scalars[jVar] / (molarMasses[speciesIndex_j] / 1000);
    //     scalars_sum += val_scalars[jVar];
    //     gradChemicalSourceTerm[iVar] +=
    //         factor * (-dT_dYi / Temperature +
    //                   meanMolecularWeight * (1 / molarMasses[speciesN] - 1 / molarMasses[speciesIndex_i]));
    //   } else {
    //     const su2double factor = dW_dCi * Density / (molarMasses[speciesIndex_j] / 1000);
    //     scalars_sum += val_scalars[jVar];
    //     gradChemicalSourceTerm[iVar] +=
    //         factor * (1.0 + val_scalars[jVar] *
    //                             (-dT_dYi / Temperature +
    //                              meanMolecularWeight * (1 / molarMasses[speciesN] - 1 / molarMasses[speciesIndex_i])));
    //   }
    // }
    // su2double dW_dCN = dW_dC.coeff(iVar, speciesN)* molarMasses[speciesN];
    // const su2double factor = dW_dCN * Density / (molarMasses[speciesN] / 1000);
    // gradChemicalSourceTerm[iVar] +=
    //     factor *
    //     (-1.0 + (1.0 - scalars_sum) * (-dT_dYi / Temperature + meanMolecularWeight * (1 / molarMasses[speciesN] -
    //                                                                                   1 / molarMasses[speciesIndex_i])));
  }
}

void CFluidCantera::ComputeHeatRelease() {
  const int nsp = sol->thermo()->nSpecies();
  vector<su2double> netProductionRates(nsp);
  sol->kinetics()->getNetProductionRates(&netProductionRates[0]);
  Heat_Release = 0.0;
  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    int speciesIndex = sol->thermo()->speciesIndex(gasComposition[iVar]);
    Heat_Release +=
        -1.0 * enthalpyFormation[speciesIndex] * molarMasses[speciesIndex] * netProductionRates[speciesIndex];
  }
}

void CFluidCantera::GetEnthalpyDiffusivity(su2double* enthalpy_diffusions) {
  const int nsp = sol->thermo()->nSpecies();
  // The universal gas constant times temperature is retrieved from cantera.
  const su2double uni_gas_constant_temp = sol->thermo()->RT();
  vector<su2double> enthalpiesSpecies(nsp);
  vector<su2double> diff(nsp);
  sol->thermo()->getEnthalpy_RT_ref(&enthalpiesSpecies[0]);
  sol->transport()->getMixDiffCoeffsMass(&diff[0]);
  const int speciesN = sol->thermo()->speciesIndex(gasComposition[n_species_mixture - 1]);
  for (int iVar = 0; iVar < n_species_mixture - 1; iVar++) {
    int speciesIndex = sol->thermo()->speciesIndex(gasComposition[iVar]);
    enthalpy_diffusions[iVar] = Density * uni_gas_constant_temp *
                                ((enthalpiesSpecies[speciesIndex] * diff[speciesIndex] / molarMasses[speciesIndex]) -
                                 (enthalpiesSpecies[speciesN] * diff[speciesN] / molarMasses[speciesN]));
  }
}

void CFluidCantera::GetMassCorrectionDiffusivity(su2double* massCorrection_diffusions) {
  const int nsp = sol->thermo()->nSpecies();
  vector<su2double> diff(nsp);
  sol->transport()->getMixDiffCoeffsMass(&diff[0]);
  const int speciesN = sol->thermo()->speciesIndex(gasComposition[n_species_mixture - 1]);
  for (int iVar = 0; iVar < n_species_mixture - 1; iVar++) {
    int speciesIndex = sol->thermo()->speciesIndex(gasComposition[iVar]);
    massCorrection_diffusions[iVar] = Density * (diff[speciesIndex] - diff[speciesN]);
  }
}

void CFluidCantera::GetGradEnthalpyDiffusivity(su2double* grad_enthalpy_diffusions) {
  const int nsp = sol->thermo()->nSpecies();
  // The universal gas constant is retrieved from cantera,in order to keep consistency with the values retrieve from it.
  const su2double universal_gas_constant = (sol->thermo()->RT()) / Temperature;
  vector<su2double> specificHeatSpecies(nsp);
  vector<su2double> diff(nsp);
  sol->thermo()->getCp_R_ref(&specificHeatSpecies[0]);
  sol->transport()->getMixDiffCoeffsMass(&diff[0]);
  const int speciesN = sol->thermo()->speciesIndex(gasComposition[n_species_mixture - 1]);
  for (int iVar = 0; iVar < n_species_mixture - 1; iVar++) {
    int speciesIndex = sol->thermo()->speciesIndex(gasComposition[iVar]);
    grad_enthalpy_diffusions[iVar] = Density * universal_gas_constant *
                               ((specificHeatSpecies[speciesIndex] * diff[speciesIndex] / molarMasses[speciesIndex]) -
                                (specificHeatSpecies[speciesN] * diff[speciesN] / molarMasses[speciesN]));
  }
}

void CFluidCantera::SetTDState_T(const su2double val_temperature, const su2double* val_scalars) {
  Temperature = val_temperature;
  /*--- Set mass fractions. ---*/
  const int nsp = sol->thermo()->nSpecies();
  su2double val_scalars_sum{0.0};
  su2double massFractions[nsp]{0.0};
  for (int i_scalar = 0; i_scalar < n_species_mixture - 1; i_scalar++) {
    int speciesIndex = sol->thermo()->speciesIndex(gasComposition[i_scalar]);
    massFractions[speciesIndex] = val_scalars[i_scalar];
    val_scalars_sum += val_scalars[i_scalar];
  }
  massFractions[sol->thermo()->speciesIndex(gasComposition[n_species_mixture - 1])] = 1.0 - val_scalars_sum;
  sol->thermo()->setMassFractions(massFractions);
  sol->thermo()->setState_TP(Temperature, Pressure_Thermodynamic);
  Density = sol->thermo()->density();
  Enthalpy = sol->thermo()->enthalpy_mass();
  Cp = sol->thermo()->cp_mass();
  Cv = sol->thermo()->cv_mass();
  Mu = sol->transport()->viscosity();
  Kt = sol->transport()->thermalConductivity();

  ComputeMassDiffusivity();
  //ComputeChemicalSourceTerm();
  ComputeGradChemicalSourceTerm(val_scalars);
  ComputeHeatRelease();
}

void CFluidCantera::SetTDState_h(const su2double val_enthalpy, const su2double* val_scalars) {
  Enthalpy = val_enthalpy;
  /*--- convergence criterion for temperature in [K], high accuracy needed for restarts. ---*/
  su2double toll = 1e-5;
  su2double temp_iter = 300.0;
  su2double delta_temp_iter = 1e10;
  su2double delta_enthalpy_iter;
  const int counter_limit = 20;

  int counter = 0;

  /*--- Set mass fractions. ---*/
  const int nsp = sol->thermo()->nSpecies();
  su2double val_scalars_sum{0.0};
  su2double massFractions[nsp]{0.0};
  for (int i_scalar = 0; i_scalar < n_species_mixture - 1; i_scalar++) {
    int speciesIndex = sol->thermo()->speciesIndex(gasComposition[i_scalar]);
    massFractions[speciesIndex] = val_scalars[i_scalar];
    val_scalars_sum += val_scalars[speciesIndex];
  }
  massFractions[sol->thermo()->speciesIndex(gasComposition[n_species_mixture - 1])] = 1.0 - val_scalars_sum;
  sol->thermo()->setMassFractions(massFractions);

  /*--- Computing temperature given enthalpy and species mass fractions using Newton-Raphson. ---*/
  while ((abs(delta_temp_iter) > toll) && (counter++ < counter_limit)) {
    /*--- Set thermodynamic state based on the current value of temperature. ---*/
    sol->thermo()->setState_TP(temp_iter, Pressure_Thermodynamic);

    su2double Enthalpy_iter = sol->thermo()->enthalpy_mass();
    su2double Cp_iter = sol->thermo()->cp_mass();

    delta_enthalpy_iter = Enthalpy - Enthalpy_iter;

    delta_temp_iter = delta_enthalpy_iter / Cp_iter;

    temp_iter += delta_temp_iter;
    if (temp_iter < 0.0) {
      cout << "Warning: Negative temperature has been found during Newton-Raphson" << endl;
      break;
    }
  }
  Temperature = temp_iter;
  if (counter == counter_limit) {
    cout << "Warning Newton-Raphson exceed number of max iteration in temperature computation" << endl;
  }
  SetTDState_T(Temperature, val_scalars);
}

#else
CFluidCantera::CFluidCantera(su2double value_pressure_operating, const CConfig* config) {
  SU2_MPI::Error("SU2 was not compiled with Cantera(-Denable-cantera=true)", CURRENT_FUNCTION);
}
#endif