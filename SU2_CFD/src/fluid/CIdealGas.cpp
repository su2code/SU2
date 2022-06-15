/*!
 * \file CIdealGas.cpp
 * \brief Source of the ideal gas model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
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

#include "../../include/fluid/CIdealGas.hpp"
#include "../../src/numerics/MultiLayer_Perceptron/LookUp_MLP.cpp"
#include "../../src/numerics/MultiLayer_Perceptron/NeuralNetwork.cpp"
#include "../../src/numerics/MultiLayer_Perceptron/ReadNeuralNetwork.cpp"
CIdealGas::CIdealGas(su2double gamma, su2double R, bool CompEntropy) : CFluidModel() {
  Gamma = gamma;
  Gamma_Minus_One = Gamma - 1.0;
  Gas_Constant = R;
  Cp = Gamma / Gamma_Minus_One * Gas_Constant;
  Cv = Cp - R;

  ComputeEntropy = CompEntropy;
  ANN = new LookUp_MLP("test_collection.mlp");

}

void CIdealGas::SetTDState_rhoe(su2double rho, su2double e) {

  vector<string> input_names;
  vector<string> output_names;
  vector<su2double> inputs;
  vector<su2double*> outputs;
  su2double S, ds_de, ds_drho, d2s_dedrho, d2s_de2, d2s_drho2;
  input_names.push_back("rho");
  inputs.push_back(rho);
  input_names.push_back("e");
  inputs.push_back(e);
//   s_norm
// ds/de_norm
// ds/drho_norm
// d2s/de.drho_norm
// d2s/de2_norm
// d2s/drho2_norm
  output_names.push_back("s");
  outputs.push_back(&S);
  output_names.push_back("ds/de");
  outputs.push_back(&ds_de);
  output_names.push_back("ds/drho");
  outputs.push_back(&ds_drho);
  output_names.push_back("d2s/de.drho");
  outputs.push_back(&d2s_dedrho);
  output_names.push_back("d2s/de2");
  outputs.push_back(&d2s_de2);
  output_names.push_back("d2s/drho2");
  outputs.push_back(&d2s_drho2);

  ANN->Predict_MLP(input_names, inputs, output_names, outputs);
  su2double specific_volume = 1/rho;
  su2double ds_dv = -ds_drho * rho*rho;

  su2double blue_term = (ds_drho * (2 - rho * (1/ds_de) * d2s_dedrho) + rho * d2s_drho2);
  su2double green_term = (- (1.0/ds_de) * d2s_de2 * ds_drho + d2s_dedrho);
  //SoundSpeed2 = - rho * (1/ds_de) * (blue_term - rho * green_term * (ds_drho / ds_de));
  Temperature = 1.0 / ds_de;
  Pressure = -pow(rho, 2) * Temperature * ds_drho;
  Density = rho;
  StaticEnergy = e;
  Pressure = Gamma_Minus_One * Density * StaticEnergy;
  //Temperature = Gamma_Minus_One * StaticEnergy / Gas_Constant;
  SoundSpeed2 = Gamma * Pressure / Density;
  dPdrho_e = Gamma_Minus_One * StaticEnergy;
  dPde_rho = Gamma_Minus_One * Density;
  dTdrho_e = 0.0;
  dTde_rho = Gamma_Minus_One / Gas_Constant;

  if (ComputeEntropy) Entropy = (1.0 / Gamma_Minus_One * log(Temperature) + log(1.0 / Density)) * Gas_Constant;
}

void CIdealGas::SetTDState_PT(su2double P, su2double T) {
  cout << "Calling TDState_PT" << endl;
  su2double e = T * Gas_Constant / Gamma_Minus_One;
  su2double rho = P / (T * Gas_Constant);
  SetTDState_rhoe(rho, e);
}

void CIdealGas::SetTDState_Prho(su2double P, su2double rho) {
  cout << "Calling TDState_Prho" << endl;
  su2double e = P / (Gamma_Minus_One * rho);
  SetEnergy_Prho(P, rho);
  SetTDState_rhoe(rho, StaticEnergy);
}

//void CIdealGas::SetEnergy_Prho(su2double P, su2double rho) {cout << "Calling SetEnergy_Prho"<<endl; StaticEnergy = P / (rho * Gamma_Minus_One); }

void CIdealGas::SetEnergy_Prho(su2double P, su2double rho){
  cout << "Calling SetEnergy_Prho"<<endl;
  su2double e_current{0.5*(3.32151385e+05 + 4.42868514e+05)};
  su2double tolerance = 10;
  unsigned long maxIter = 200, Iter{0};
  
  vector<string> input_names;
  vector<string> output_names;
  vector<su2double> inputs;
  vector<su2double*> outputs;
  su2double S, ds_de, ds_drho, d2s_dedrho, d2s_de2, d2s_drho2, dP_de;
  input_names.push_back("rho");
  inputs.push_back(rho);
  input_names.push_back("e");
  inputs.push_back(e_current);

  output_names.push_back("s");
  outputs.push_back(&S);
  output_names.push_back("ds/de");
  outputs.push_back(&ds_de);
  output_names.push_back("ds/drho");
  outputs.push_back(&ds_drho);
  output_names.push_back("d2s/de.drho");
  outputs.push_back(&d2s_dedrho);
  output_names.push_back("d2s/de2");
  outputs.push_back(&d2s_de2);
  output_names.push_back("d2s/drho2");
  outputs.push_back(&d2s_drho2);

  ANN->Predict_MLP(input_names, inputs, output_names, outputs); 
  su2double Temperature_current = 1.0 / ds_de;
  su2double Pressure_current = -pow(rho, 2) * Temperature * ds_drho;
  su2double delta_P = P - Pressure_current;
  su2double delta_e;

  while((abs(delta_P) > tolerance) && (Iter < maxIter)){
      dP_de = rho*rho*(pow(ds_de, -2)*d2s_de2*ds_drho - Temperature_current*d2s_dedrho);
      delta_e = (P - Pressure_current)/dP_de;
      e_current = e_current + 0.5*delta_e;
      inputs.at(1) = e_current;

      ANN->Predict_MLP(input_names, inputs, output_names, outputs); 
      Temperature_current = 1.0 / ds_de;
      Pressure_current = -pow(rho, 2) * Temperature * ds_drho;
      delta_P = P - Pressure_current;
      cout << delta_P << " " << dP_de << endl;
      Iter ++;
  }
  if(Iter == maxIter){
    cout << "Maximum iteration count reached!" << endl;
  }else{
    cout << "Successfull thingy!" << endl;
  }
  StaticEnergy = e_current;
}
void CIdealGas::SetTDState_hs(su2double h, su2double s) {
  cout << "Calling TDState_hs" << endl;
  su2double T = h * Gamma_Minus_One / Gas_Constant / Gamma;
  su2double e = h / Gamma;
  su2double v = exp(-1 / Gamma_Minus_One * log(T) + s / Gas_Constant);

  SetTDState_rhoe(1 / v, e);
}

void CIdealGas::SetTDState_Ps(su2double P, su2double s) {
  cout << "Calling TDState_Ps" << endl;
  su2double T = exp(Gamma_Minus_One / Gamma * (s / Gas_Constant + log(P) - log(Gas_Constant)));
  su2double rho = P / (T * Gas_Constant);

  SetTDState_Prho(P, rho);
}

void CIdealGas::SetTDState_rhoT(su2double rho, su2double T) {
  cout << "Calling TDState_rhoT" << endl;
  su2double e = T * Gas_Constant / Gamma_Minus_One;
  SetTDState_rhoe(rho, e);
}

void CIdealGas::ComputeDerivativeNRBC_Prho(su2double P, su2double rho) {
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
