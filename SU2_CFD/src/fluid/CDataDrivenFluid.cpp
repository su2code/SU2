/*!
 * \file CDataDrivenFluid.cpp
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

#include "../../include/fluid/CDataDrivenFluid.hpp"
#include "../../../Common/include/toolboxes/multilayer_perceptron/CLookUp_ANN.hpp"

#include <fstream>

CDataDrivenFluid::CDataDrivenFluid(su2double gamma, su2double R, bool CompEntropy) : CFluidModel() {
  Gamma = gamma;
  Gamma_Minus_One = Gamma - 1.0;
  Gas_Constant = R;
  Cp = Gamma / Gamma_Minus_One * Gas_Constant;
  Cv = Cp - R;

  idx_rho = 0;
  idx_e = 1;

  ComputeEntropy = CompEntropy;
  //ANN = new LookUp_MLP("split_dataset.mlp");
  ANN_new = new MLPToolbox::CLookUp_ANN("split_dataset.mlp");

  input_names.push_back("rho");
  input_names_lower_rho.push_back("rho_l");
  inputs.push_back(Density);
  input_names.push_back("e");
  input_names_lower_rho.push_back("e_l");
  inputs.push_back(StaticEnergy);


  output_names.push_back("s");
  outputs.push_back(&Entropy);
  output_names.push_back("ds_de");
  outputs.push_back(&ds_de);
  output_names.push_back("ds_drho");
  outputs.push_back(&ds_drho);
  output_names.push_back("d2s_dedrho");
  outputs.push_back(&d2s_dedrho);
  output_names.push_back("d2s_de2");
  outputs.push_back(&d2s_de2);
  output_names.push_back("d2s_drho2");
  outputs.push_back(&d2s_drho2);

  input_output_map = new MLPToolbox::CIOMap(ANN_new, input_names, output_names);
  input_output_map_lower_rho = new MLPToolbox::CIOMap(ANN_new, input_names_lower_rho, output_names);

  // rho_maxmin = ANN->GetInputNorm(input_names.at(idx_rho));
  // e_maxmin = ANN->GetInputNorm(input_names.at(idx_e));
  rho_maxmin = make_pair(0.80523, 295.251);
  e_maxmin = make_pair(350604.24, 479774.22);
  rho_start = 0.5*(rho_maxmin.first + rho_maxmin.second);
  e_start = 0.5*(e_maxmin.first + e_maxmin.second);

  // ifstream infile("/home/evert/NICFD/testcase/roedata_ws.csv");
  // ofstream outfile("refdata_prediction_newMLP.csv");
  // string line;
  // getline(infile, line);

  // su2double rho_in, e_in;
  // vector<string> input_names;
  // vector<string> output_names;
  // vector<su2double> inputs;
  // vector<su2double*> outputs;
  // su2double S, ds_de, ds_drho, d2s_dedrho, d2s_de2, d2s_drho2;
  // input_names.push_back("rho");
  // inputs.push_back(rho_in);
  // input_names.push_back("e");
  // inputs.push_back(e_in);
  // output_names.push_back("s");
  // outputs.push_back(&S);
  // output_names.push_back("ds/de");
  // outputs.push_back(&ds_de);
  // output_names.push_back("ds/drho");
  // outputs.push_back(&ds_drho);
  // output_names.push_back("d2s/de.drho");
  // outputs.push_back(&d2s_dedrho);
  // output_names.push_back("d2s/de2");
  // outputs.push_back(&d2s_de2);
  // output_names.push_back("d2s/drho2");
  // outputs.push_back(&d2s_drho2);

  // while(getline(infile, line)){
  //   istringstream line_input(line);
  //   line_input >> rho_in >> e_in;
  //   inputs.at(0) = rho_in;
  //   inputs.at(1) = e_in;
  //   Predict_MLP(inputs);
  //   outfile << rho_in <<"\t" << e_in <<"\t"<< Entropy <<"\t" << ds_de << "\t" << ds_drho << "\t" << d2s_dedrho << "\t" << d2s_de2 << "\t" << d2s_drho2 << endl;
  // }
  // infile.close();
  // outfile.close();

  dOutputs_dInputs = new su2double*[6];
  for(size_t i=0; i<6; i++){
    dOutputs_dInputs[i] = new su2double[2];
  }
}


void CDataDrivenFluid::SetTDState_rhoe(su2double rho, su2double e) {

  boundsviolation = 0;

  Predict_MLP(rho, e);

  su2double blue_term = (ds_drho * (2 - rho * pow(ds_de, -1) * d2s_dedrho) + rho * d2s_drho2);
  su2double green_term = (- pow(ds_de, -1) * d2s_de2 * ds_drho + d2s_dedrho);
  
  SoundSpeed2 = - rho * pow(ds_de, -1) * (blue_term - rho * green_term * (ds_drho / ds_de));

  Temperature = 1.0 / ds_de;
  Pressure = -pow(rho, 2) * Temperature * ds_drho;
  Density = rho;
  StaticEnergy = e;
  
  dTde_rho = -pow(ds_de, -2)*d2s_de2;
  dTdrho_e = 0.0;

  dPde_rho = -pow(rho, 2) * dTde_rho * ds_drho;
  dPdrho_e = -2*rho*Temperature*ds_drho - pow(rho, 2)*Temperature*d2s_drho2;
  failedNewtonSolver = false;
}

void CDataDrivenFluid::SetTDState_PT(su2double P, su2double T) {
  //cout << "Start SetTDState_PT" << endl;
  su2double rho = rho_start;
  su2double e = e_start;

  bool converged{false};
  su2double tolerance_P = 10,
            tolerance_T = 1,
            relaxation = 0.1;
  unsigned long iter_max = 1000, Iter{0};

  su2double T_current, P_current, delta_P, delta_T, delta_rho, delta_e;
  su2double dP_de, dP_drho, dT_de, dT_drho, determinant;

  vector<su2double> rho_trend, e_trend, s_trend, ds_de_trend, ds_drho_trend, d2s_dedrho_trend, d2s_de2_trend, d2s_drho2_trend;
  vector<su2double> determinant_trend, dP_de_trend, dP_drho_trend, dT_de_trend, dT_drho_trend, delta_rho_trend, delta_e_trend;

  while(!converged && (Iter < iter_max)){
    // ANN->Predict_MLP(input_names, inputs, output_names, outputs);

    // ds_drho = -exp(ds_drho);
    // d2s_drho2 = exp(d2s_drho2);
    Predict_MLP(rho, e);
    //ds_drho = dOutputs_dInputs[0][idx_rho];
    //ds_de = dOutputs_dInputs[0][idx_e];

    T_current = pow(ds_de, -1);
    P_current = -pow(rho, 2)*T_current*ds_drho;
    delta_P = P_current - P;
    delta_T = T_current - T;
    //if(Iter == 0) cout << "Start: " << rho << " " << e << " " << P_current << " " << P << " " << T_current << " " << T << endl;
    // if(!thingy){
    //   rho_trend.push_back(rho);
    //   e_trend.push_back(e);
    //   s_trend.push_back(Entropy);
    //   ds_de_trend.push_back(ds_de);
    //   ds_drho_trend.push_back(ds_drho);
    //   d2s_dedrho_trend.push_back(d2s_dedrho);
    //   d2s_de2_trend.push_back(d2s_de2);
    //   d2s_drho2_trend.push_back(d2s_drho2);
    // }
    if((abs(delta_P) < tolerance_P) && (abs(delta_T) < tolerance_T)){
      converged = true;
    }else{
      dP_de = pow(rho, 2)*pow(ds_de, -2)*d2s_de2*ds_drho;
      dP_drho = -2*rho*T_current*ds_drho - pow(rho, 2)*T_current*d2s_drho2;

      dT_de = -pow(ds_de, -2)*d2s_de2;
      dT_drho = 0;

      determinant = dP_drho * dT_de - dP_de * dT_drho;
      delta_rho = (dT_de * delta_P - dP_de * delta_T) / determinant;
      delta_e = (-dT_drho * delta_P + dP_drho * delta_T) / determinant;

      rho -= relaxation * delta_rho;
      e -= relaxation * delta_e;
      // rho = max(min(rho, rho_maxmin.second), rho_maxmin.first);
      // e = max(min(e, e_maxmin.second), e_maxmin.first);
    }
    // if(!thingy){
    //   rho_trend.push_back(rho);
    //   e_trend.push_back(e);
    //   s_trend.push_back(Entropy);
    //   ds_de_trend.push_back(ds_de);
    //   ds_drho_trend.push_back(ds_drho);
    //   d2s_dedrho_trend.push_back(d2s_dedrho);
    //   d2s_de2_trend.push_back(d2s_de2);
    //   d2s_drho2_trend.push_back(d2s_drho2);
    //   determinant_trend.push_back(determinant);
    //   dP_de_trend.push_back(dP_de);
    //   dP_drho_trend.push_back(dP_drho);
    //   dT_de_trend.push_back(dT_de);
    //   dT_drho_trend.push_back(dT_drho);
    //   delta_rho_trend.push_back(delta_rho);
    //   delta_e_trend.push_back(delta_e);
    // }
    Iter ++;
  }

  // if(!thingy){
  //   ofstream outfile{"PT_trend.csv"};
  //   for(size_t i=0; i<rho_trend.size(); i++){
  //     outfile << rho_trend.at(i) <<"\t" << e_trend.at(i) <<"\t"<< s_trend.at(i) <<"\t" << ds_de_trend.at(i) << "\t" << ds_drho_trend.at(i) << "\t" << d2s_dedrho_trend.at(i) << "\t" << d2s_de2_trend.at(i) << "\t" << d2s_drho2_trend.at(i) << endl;
  //   }
  //   outfile.close();
  //   ofstream outfile_2{"PT_details_trend.csv"};
  //   for(size_t i=0; i<rho_trend.size(); i++){
  //     outfile_2 << determinant_trend.at(i) <<"\t" << dP_de_trend.at(i) <<"\t"<< dP_drho_trend.at(i) <<"\t" << dT_de_trend.at(i) << "\t" << dT_drho_trend.at(i) << "\t" << delta_rho_trend.at(i) << "\t" << delta_e_trend.at(i) << endl;
  //   }
  //   outfile_2.close();
  // }
    
  

  //cout << "End: " << rho << " " << e << " " << P_current << " " << P << " " << T_current << " " << T << endl;
  Reevaluate_MLP = false;
  SetTDState_rhoe(rho, e);
  Reevaluate_MLP = true;
  failedNewtonSolver = converged;
  nIter_NewtonSolver = Iter;
  //cout << "End SetTDState_PT" << endl;
}

void CDataDrivenFluid::SetTDState_Prho(su2double P, su2double rho) {
  //cout << "SetTDState_Prho" << endl;
  //cout << "Calling SetTDState_Prho" << endl;
  SetEnergy_Prho(P, rho);
  SetTDState_rhoe(rho, StaticEnergy);
  //cout << "Pressure: " << Pressure << " Temperature: " << Temperature << " Soundspeed: " << sqrt(SoundSpeed2) << endl;
}

//void CDataDrivenFluid::SetEnergy_Prho(su2double P, su2double rho) {cout << "Calling SetEnergy_Prho"<<endl; StaticEnergy = P / (rho * Gamma_Minus_One); }

void CDataDrivenFluid::SetEnergy_Prho(su2double P, su2double rho){
  //cout << "Calling SetEnergy_Prho"<<endl;
  su2double e = e_start;

  su2double tolerance = 10;
  su2double relaxation = 0.1;
  bool converged{false};
  unsigned long maxIter = 1000, Iter{0};
  
  su2double T_current, P_current, dP_de, delta_P, delta_e;
  while(!converged && (Iter < maxIter)){
      
      // ANN->Predict_MLP(input_names, inputs, output_names, outputs);
      // ds_drho = -exp(ds_drho);
      // d2s_drho2 = exp(d2s_drho2);
      Predict_MLP(rho, e);
      //ds_drho = dOutputs_dInputs[0][idx_rho];
      //ds_de = dOutputs_dInputs[0][idx_e];

      T_current = 1/ds_de;
      P_current = -pow(rho, 2)*T_current*ds_drho;
      delta_P = P_current - P;
      if(abs(delta_P) < tolerance){
        converged = true;
      }else{
        dP_de = pow(rho, 2)*pow(ds_de, -2)*d2s_de2*ds_drho;
        delta_e = delta_P / dP_de;
        e -= relaxation * delta_e;
        //e = max(min(e, e_maxmin.second), e_maxmin.first);
      }
      Iter ++;
  }
  StaticEnergy = e;
  failedNewtonSolver = converged;
  nIter_NewtonSolver = Iter;
}
void CDataDrivenFluid::SetTDState_hs(su2double h, su2double s) {
  //cout << "Start SetTDState_hs" << endl;
  su2double e = e_start;
  su2double rho = rho_start;

  bool converged{false};
  su2double tolerance_h = 10,
            tolerance_s = 1,
            relaxation = 0.1;
  unsigned long iter_max = 1000, Iter{0};

  su2double T_current, P_current, h_current, s_current;
  su2double delta_h, delta_s, delta_rho, delta_e, dP_de, dP_drho, dT_de, dT_drho, dh_drho, dh_de, determinant;
  while(!converged && (Iter < iter_max)){
    //ANN->Predict_MLP(input_names, inputs, output_names, outputs);
    //ds_drho = -exp(ds_drho);
    //d2s_drho2 = exp(d2s_drho2);
    Predict_MLP(rho, e);
    
    //ds_drho = dOutputs_dInputs[0][idx_rho];
    //ds_de = dOutputs_dInputs[0][idx_e];

    T_current = pow(ds_de, -1);
    P_current = -pow(rho, 2)*T_current*ds_drho;
    h_current = e + P_current / rho;

    s_current = Entropy;
    // if(Iter == 0){
    //   cout << "Start: " << rho << " " << e << " " << h_current << " " << h << " " << s_current << " " << s << endl;
    // }
    delta_h = h_current - h;
    delta_s = s_current - s;
    if((abs(delta_h) < tolerance_h) && (abs(delta_s) < tolerance_s)){
      converged = true;
    }else{
      dP_de = pow(rho, 2)*pow(ds_de, -2)*d2s_de2*ds_drho;
      dP_drho = -2*rho*T_current*ds_drho - pow(rho, 2)*T_current*d2s_drho2;

      dh_de = 1 + dP_de / rho;
      dh_drho = -P_current*pow(rho, -2) + dP_drho / rho;

      determinant = dh_drho * ds_de - dh_de * ds_drho;
      delta_rho = (ds_de * delta_h - dh_de * delta_s) / determinant;
      delta_e = (-ds_drho * delta_h + dh_drho * delta_s) / determinant;

      rho -= relaxation * delta_rho;
      e -= relaxation * delta_e;
      // rho = max(min(rho, rho_maxmin.second), rho_maxmin.first);
      // e = max(min(e, e_maxmin.second), e_maxmin.first);
    }
    Iter ++;
    
  }

  Reevaluate_MLP = false;
  SetTDState_rhoe(rho, e);
  Reevaluate_MLP = true;
  failedNewtonSolver = converged;
  nIter_NewtonSolver = Iter;
  //cout << "End SetTDState_hs" << endl;
}

void CDataDrivenFluid::SetTDState_Ps(su2double P, su2double s) {


  su2double e = e_start;
  su2double rho = rho_start;

  bool converged{false};
  su2double tolerance_P = 10,
            tolerance_s = 1,
            relaxation = 0.1;
  unsigned long iter_max = 1000, Iter{0};

  su2double T_current, P_current, h_current, s_current;
  su2double delta_P, delta_s, delta_rho, delta_e, dP_de, dP_drho, dT_de, dT_drho, dh_drho, dh_de, determinant;

  while(!converged && (Iter < iter_max)){
    //ANN->Predict_MLP(input_names, inputs, output_names, outputs);
    //ds_drho = -exp(ds_drho);
    //d2s_drho2 = exp(d2s_drho2);
    Predict_MLP(rho, e);
    
    //ds_drho = dOutputs_dInputs[0][idx_rho];
    //ds_de = dOutputs_dInputs[0][idx_e];

    T_current = pow(ds_de, -1);
    P_current = -pow(rho, 2)*T_current*ds_drho;
    h_current = e + P_current / rho;

    s_current = Entropy;
    // if(Iter == 0){
    //   cout << "Start: " << rho << " " << e << " " << h_current << " " << h << " " << s_current << " " << s << endl;
    // }
    delta_P = P_current - P;
    delta_s = s_current - s;
    if((abs(delta_P) < tolerance_P) && (abs(delta_s) < tolerance_s)){
      converged = true;
    }else{
      dP_de = pow(rho, 2)*pow(ds_de, -2)*d2s_de2*ds_drho;
      dP_drho = -2*rho*T_current*ds_drho - pow(rho, 2)*T_current*d2s_drho2;

      determinant = dP_drho * ds_de - dP_de * ds_drho;
      delta_rho = (ds_de * delta_P - dP_de * delta_s) / determinant;
      delta_e = (-ds_drho * delta_P + dP_drho * delta_s) / determinant;

      rho -= relaxation * delta_rho;
      e -= relaxation * delta_e;
      // rho = max(min(rho, rho_maxmin.second), rho_maxmin.first);
      // e = max(min(e, e_maxmin.second), e_maxmin.first);
    }
    Iter ++;
    
  }

  Reevaluate_MLP = false;
  SetTDState_rhoe(rho, e);
  Reevaluate_MLP = true;
  failedNewtonSolver = converged;
  nIter_NewtonSolver = Iter;
}

void CDataDrivenFluid::SetTDState_rhoT(su2double rho, su2double T) {
  cout << "Calling TDState_rhoT" << endl;
  su2double e = T * Gas_Constant / Gamma_Minus_One;
  SetTDState_rhoe(rho, e);
}

void CDataDrivenFluid::ComputeDerivativeNRBC_Prho(su2double P, su2double rho) {
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

void CDataDrivenFluid::Predict_MLP(su2double rho, su2double e){
  
  inputs.at(idx_rho) = rho;
  inputs.at(idx_e) = e;
  vector<string> these_outputs;
  if(rho < 10.011693){
    //these_outputs = output_names_lower;
    ANN_new->Predict_ANN(input_output_map_lower_rho, inputs, outputs);
  }else{
    //these_outputs = output_names;
    ANN_new->Predict_ANN(input_output_map, inputs, outputs);
  }
  //ANN->Predict_MLP(input_names, inputs, these_outputs, outputs);
  ds_drho = -exp(ds_drho);
  d2s_drho2 = exp(d2s_drho2);
}