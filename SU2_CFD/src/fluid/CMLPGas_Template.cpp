/*!
 * \file CMLPGas_Template.cpp
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

#include "../../include/fluid/CMLPGas_Template.hpp"

CMLPGas_Template::CMLPGas_Template(su2double gamma, su2double R, bool CompEntropy) : CFluidModel() {
  Gamma = gamma;
  Gamma_Minus_One = Gamma - 1.0;
  Gas_Constant = R;
  Cp = Gamma / Gamma_Minus_One * Gas_Constant;
  Cv = Cp - R;

  ComputeEntropy = CompEntropy;

  /* Define a CLookUp_ANN object to allow for data regression */
  ann_input_filename = "MLP_collection.mlp";
  lookup_ann = new MLPToolbox::CLookUp_ANN(ann_input_filename);

  /* Map MLP inputs to outputs for each look-up operation*/
  MapInputs_to_Outputs();
}

CMLPGas_Template::~CMLPGas_Template(){
  delete iomap_rhoe;
  delete iomap_PT;
  delete iomap_Prho;
  delete lookup_ann;
}
void CMLPGas_Template::MapInputs_to_Outputs(){

  input_names_rhoe.push_back("rho");
  input_names_rhoe.push_back("e");

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

  input_names_PT.push_back("P");
  input_names_PT.push_back("T");

  output_names_PT.push_back("rho");
  outputs_PT.push_back(&Density);
  output_names_PT.push_back("e");
  outputs_PT.push_back(&StaticEnergy);

  iomap_PT = new MLPToolbox::CIOMap(lookup_ann, input_names_PT, output_names_PT);

  input_names_Prho.push_back("P");
  input_names_Prho.push_back("rho");

  output_names_Prho.push_back("e");
  outputs_Prho.push_back(&StaticEnergy);
  iomap_Prho = new MLPToolbox::CIOMap(lookup_ann, input_names_Prho, output_names_Prho);
}

void CMLPGas_Template::SetTDState_rhoe(su2double rho, su2double e) {
  vector<su2double> ann_inputs;
  ann_inputs.push_back(rho);
  ann_inputs.push_back(e);
  lookup_ann->Predict_ANN(iomap_rhoe, ann_inputs, outputs_rhoe);
}

void CMLPGas_Template::SetTDState_PT(su2double P, su2double T) {
  vector<su2double> ann_inputs;
  ann_inputs.push_back(P);
  ann_inputs.push_back(T);
  lookup_ann->Predict_ANN(iomap_PT, ann_inputs, outputs_PT);

  SetTDState_rhoe(Density, StaticEnergy);
}

void CMLPGas_Template::SetTDState_Prho(su2double P, su2double rho) {
  vector<su2double> ann_inputs;
  ann_inputs.push_back(P);
  ann_inputs.push_back(rho);
  lookup_ann->Predict_ANN(iomap_Prho, ann_inputs, outputs_Prho);
  SetTDState_rhoe(rho, StaticEnergy);
}

void CMLPGas_Template::SetEnergy_Prho(su2double P, su2double rho) { 
  vector<su2double> ann_inputs;
  ann_inputs.push_back(P);
  ann_inputs.push_back(rho);
  lookup_ann->Predict_ANN(iomap_Prho, ann_inputs, outputs_Prho);
}
