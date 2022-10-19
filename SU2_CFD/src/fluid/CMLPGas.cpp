/*!
 * \file CMLPGas.cpp
 * \brief Source of the data-driven fluid model class using 
 *  multilayer perceptrons for data regression
 * \author E.Bunschoten
 * \version 7.4.0 "Blackbird"
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

#include "../../include/fluid/CMLPGas.hpp"

CMLPGas::CMLPGas(const CConfig* config) : CFluidModel() {

  /* Define a CLookUp_ANN object to allow for data regression */
  ann_input_filename = config->GetMLP_FileName();
  lookup_ann = new MLPToolbox::CLookUp_ANN(ann_input_filename);

  /* Map MLP inputs to outputs for each look-up operation*/
  MapInputs_to_Outputs();
}

CMLPGas::~CMLPGas(){
  delete iomap_rhoe;
  delete iomap_PT;
  delete iomap_Prho;
  delete lookup_ann;
}
void CMLPGas::MapInputs_to_Outputs(){

  input_names_rhoe.push_back("Density");
  input_names_rhoe.push_back("Energy");

  output_names_rhoe.push_back("Temperature");
  outputs_rhoe.push_back(&Temperature);
  output_names_rhoe.push_back("Pressure");
  outputs_rhoe.push_back(&Pressure);
  output_names_rhoe.push_back("SoundSpeed2");
  outputs_rhoe.push_back(&SoundSpeed2);
  output_names_rhoe.push_back("dPdrho_e");
  outputs_rhoe.push_back(&dPdrho_e);
  output_names_rhoe.push_back("dPde_rho");
  outputs_rhoe.push_back(&dPde_rho);
  output_names_rhoe.push_back("dTdrho_e");
  outputs_rhoe.push_back(&dTdrho_e);
  output_names_rhoe.push_back("dTde_rho");
  outputs_rhoe.push_back(&dTde_rho);
  output_names_rhoe.push_back("Entropy");
  outputs_rhoe.push_back(&Entropy);

  iomap_rhoe = new MLPToolbox::CIOMap(lookup_ann, input_names_rhoe, output_names_rhoe);

  input_names_PT.push_back("Pressure");
  input_names_PT.push_back("Temperature");

  output_names_PT.push_back("Density");
  outputs_PT.push_back(&Density);
  output_names_PT.push_back("Energy");
  outputs_PT.push_back(&StaticEnergy);

  iomap_PT = new MLPToolbox::CIOMap(lookup_ann, input_names_PT, output_names_PT);

  input_names_Prho.push_back("Pressure");
  input_names_Prho.push_back("Density");

  output_names_Prho.push_back("Energy");
  outputs_Prho.push_back(&StaticEnergy);
  iomap_Prho = new MLPToolbox::CIOMap(lookup_ann, input_names_Prho, output_names_Prho);
}

void CMLPGas::SetTDState_rhoe(su2double rho, su2double e) {
  Density = rho;
  StaticEnergy = e;
  vector<su2double> ann_inputs;
  ann_inputs.push_back(rho);
  ann_inputs.push_back(e);
  lookup_ann->Predict_ANN(iomap_rhoe, ann_inputs, outputs_rhoe);
  Cp = (Pressure / (R_u * Density * Temperature)) * 1 / (dTde_rho);
  Cv = (1/dTde_rho);
  Gamma = Cp / Cv;
  Gamma_Minus_One = Gamma - 1;
  Gas_Constant = Cp * Gamma_Minus_One;

}

void CMLPGas::SetTDState_PT(su2double P, su2double T) {
  // vector<su2double> ann_inputs;
  // ann_inputs.push_back(P);
  // ann_inputs.push_back(T);
  // lookup_ann->Predict_ANN(iomap_PT, ann_inputs, outputs_PT);

  su2double rho = 0.5*(6.9682399336300005e-01 + 2.4896548789900002e+00);
  su2double e = 0.5*(3.2562734877400001e+05	+4.8566089128400001e+05);
  rho = 1.0;
  e = 3.2562734877400001e+05;
  su2double delta_P, delta_T, delta_rho, delta_e, tolerance_P = 10, tolerance_T = 1.0;
  su2double determinant, relaxation = 0.1;
  unsigned long iter = 0, iter_max = 1000;
  bool converged= false;
  cout << "Target PT: " << P << " " << T << endl;
  while(!converged && (iter < iter_max)){
    SetTDState_rhoe(rho, e);
    delta_P = (P - Pressure);
    delta_T = (T - Temperature);
    cout << delta_P << " " << delta_T << endl;
    if((abs(delta_P)<tolerance_P) && (abs(delta_T) < tolerance_T)){
      converged = true;
      cout << "Converged PT: " << Pressure << " " << Temperature << endl;
      return;
    }else{
      determinant = dPdrho_e * dTde_rho - dPde_rho * dTdrho_e;
      delta_rho = (dTde_rho * delta_P - dPde_rho * delta_T) / determinant;
      delta_T = (-dTdrho_e * delta_P + dPdrho_e * delta_T) / determinant;

      rho -= relaxation * delta_rho;
      e -= relaxation * delta_e;    
    }
    iter++;
  }

  SetTDState_rhoe(rho, e);
  cout << "Final(nonconverged) PT: " << Pressure << " " << Temperature << endl;
}

void CMLPGas::SetTDState_Prho(su2double P, su2double rho) {
  vector<su2double> ann_inputs;
  ann_inputs.push_back(P);
  ann_inputs.push_back(rho);
  lookup_ann->Predict_ANN(iomap_Prho, ann_inputs, outputs_Prho);
  SetTDState_rhoe(rho, StaticEnergy);
}

void CMLPGas::SetEnergy_Prho(su2double P, su2double rho) { 
  vector<su2double> ann_inputs;
  ann_inputs.push_back(P);
  ann_inputs.push_back(rho);
  lookup_ann->Predict_ANN(iomap_Prho, ann_inputs, outputs_Prho);
}
