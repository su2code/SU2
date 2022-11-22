/*!
 * \file CDataDrivenFluid.cpp
 * \brief Source of the data-driven fluid model class
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

#include "../../include/fluid/CDataDrivenFluid.hpp"

CDataDrivenFluid::CDataDrivenFluid(const CConfig* config) : CFluidModel() {

  Kind_DataDriven_Method = config->GetKind_DataDriven_Method();

  /*--- For this branch, only MLP's are supported for data-driven fluid model. ---*/
  if(Kind_DataDriven_Method != ENUM_DATADRIVEN_METHOD::MLP)
    SU2_MPI::Error("Only multi-layer perceptrons are currently accepted for data-driven fluid models.", CURRENT_FUNCTION);

  /*--- Retrieve interpolation method file name. ---*/
  input_filename = config->GetDataDriven_Filename();

  /*--- Set up interpolation algorithm according to data-driven method. Currently only MLP's are supported. ---*/
  switch (Kind_DataDriven_Method)
  {
  case ENUM_DATADRIVEN_METHOD::MLP:
    lookup_mlp = new MLPToolbox::CLookUp_ANN(input_filename);
    break;
  case ENUM_DATADRIVEN_METHOD::LUT:
    lookup_table = new CLookUpTable(input_filename, "Density", "Energy");
    break;
  default:
    break;
  }
  
  /*--- Relaxation factor for Newton solvers. ---*/
  Newton_Relaxation = config->GetRelaxation_DataDriven();

  /*--- Preprocessing of inputs and outputs for the interpolation method. ---*/
  MapInputs_to_Outputs();

  /*--- Set initial values for density and energy based on config options ---*/
  rho_start = config->GetDensity_Init_DataDriven();
  e_start   = config->GetEnergy_Init_DataDriven();
  
}

CDataDrivenFluid::~CDataDrivenFluid(){
  switch (Kind_DataDriven_Method)
  {
  case ENUM_DATADRIVEN_METHOD::MLP:
    delete iomap_rhoe;
    delete lookup_mlp;
    break;
  case ENUM_DATADRIVEN_METHOD::LUT:
    delete lookup_table;
    break;
  default:
    break;
  }
  
}
void CDataDrivenFluid::MapInputs_to_Outputs(){

  /*--- Inputs of the data-driven method are density and internal energy. ---*/
  input_names_rhoe.push_back("Density");
  idx_rho = 0;
  input_names_rhoe.push_back("Energy");
  idx_e = 1;

  /*--- Required outputs for the interpolation method are entropy and its partial derivatives with respect to energy and density ---*/
  output_names_rhoe.push_back("s");
  outputs_rhoe.push_back(&Entropy);
  output_names_rhoe.push_back("dsde_rho");
  outputs_rhoe.push_back(&dsde_rho);
  output_names_rhoe.push_back("dsdrho_e");
  outputs_rhoe.push_back(&dsdrho_e);
  output_names_rhoe.push_back("d2sde2");
  outputs_rhoe.push_back(&d2sde2);
  output_names_rhoe.push_back("d2sdedrho");
  outputs_rhoe.push_back(&d2sdedrho);
  output_names_rhoe.push_back("d2sdrho2");
  outputs_rhoe.push_back(&d2sdrho2);

  /*--- Further preprocessing of input and output variables ---*/
  switch (Kind_DataDriven_Method)
  {
  case ENUM_DATADRIVEN_METHOD::MLP:
    /*--- Map MLP inputs to outputs ---*/
    iomap_rhoe = new MLPToolbox::CIOMap(lookup_mlp, input_names_rhoe, output_names_rhoe);
    MLP_inputs.resize(2);
    break;
  
  default:
    break;
  }
}

void CDataDrivenFluid::SetTDState_rhoe(su2double rho, su2double e) {
  Density = rho;
  StaticEnergy = e;
  switch (Kind_DataDriven_Method)
  {
  case ENUM_DATADRIVEN_METHOD::MLP:
    /* --- Set MLP input vector values --- */
    MLP_inputs[idx_rho] = Density;
    MLP_inputs[idx_e]   = StaticEnergy;

    /* --- Evaluate MLP --- */
    lookup_mlp->Predict_ANN(iomap_rhoe, MLP_inputs, outputs_rhoe);
    break;
  case ENUM_DATADRIVEN_METHOD::LUT:
    lookup_table->LookUp_ProgEnth(output_names_rhoe, outputs_rhoe, rho, e, "Density", "Energy");
  default:
    break;
  }

  su2double blue_term = (dsdrho_e * (2 - rho * pow(dsde_rho, -1) * d2sdedrho) + rho * d2sdrho2);
  su2double green_term = (- pow(dsde_rho, -1) * d2sde2 * dsdrho_e + d2sdedrho);
  
  SoundSpeed2 = - rho * pow(dsde_rho, -1) * (blue_term - rho * green_term * (dsdrho_e / dsde_rho));

  Temperature = 1.0 / dsde_rho;
  Pressure = -pow(rho, 2) * Temperature * dsdrho_e;
  Density = rho;
  StaticEnergy = e;
  
  dTde_rho = -pow(dsde_rho, -2)*d2sde2;
  dTdrho_e = 0.0;

  dPde_rho = -pow(rho, 2) * dTde_rho * dsdrho_e;
  dPdrho_e = -2*rho*Temperature*dsdrho_e - pow(rho, 2)*Temperature*d2sdrho2;

  Cp = (1 / dTde_rho) * ( 1 + (1 / Density)*dPde_rho);
  Cv = 1 / dTde_rho;
  Gamma = Cp / Cv;
  Gamma_Minus_One = Gamma - 1;
  Gas_Constant = Cp - Cv;
  
}

void CDataDrivenFluid::SetTDState_PT(su2double P, su2double T) {

  /*--- Setting initial values for density and energy ---*/
  su2double rho = rho_start, 
            e = e_start;

  su2double tolerance_P = 10, // Tolerance for pressure solution
            tolerance_T = 1;  // Tolerance for temperature solution

  bool converged = false;         // Convergence flag
  unsigned long iter_max = 1000,  // Maximum number of iterations
                Iter = 0;       

  su2double delta_P,      // Pressure residual
            delta_T,      // Temperature residual
            delta_rho,    // Density step size
            delta_e,      // Energy step size
            determinant;  // Jacobian determinant
  /*--- Initiating Newton solver ---*/
  while(!converged && (Iter < iter_max)){

    /*--- Determine thermodynamic state based on current density and energy*/
    SetTDState_rhoe(rho ,e);

    /*--- Determine pressure and temperature residuals ---*/
    delta_P = Pressure - P;
    delta_T = Temperature - T;
    /*--- Continue iterative process if residuals are outside tolerances ---*/
    if((abs(delta_P) < tolerance_P) && (abs(delta_T) < tolerance_T)){
      converged = true;
    }else{
      
      /*--- Compute step size for density and energy ---*/
      determinant = dPdrho_e * dTde_rho - dPde_rho * dTdrho_e;

      delta_rho = (dTde_rho * delta_P - dPde_rho * delta_T) / determinant;
      delta_e = (-dTdrho_e * delta_P + dPdrho_e * delta_T) / determinant;

      /*--- Update density and energy values ---*/
      rho -= Newton_Relaxation * delta_rho;
      e -= Newton_Relaxation * delta_e;
      
    }
    Iter ++;
  }

  /*--- Calculate thermodynamic state based on converged values for density and energy ---*/
  SetTDState_rhoe(rho, e);
}

void CDataDrivenFluid::SetTDState_Prho(su2double P, su2double rho) {
  /*--- Computing static energy according to pressure and density ---*/
  SetEnergy_Prho(P, rho);
  

  /*--- Calculate thermodynamic state based on converged value for energy ---*/
  SetTDState_rhoe(rho, StaticEnergy);
}

void CDataDrivenFluid::SetEnergy_Prho(su2double P, su2double rho) { 
  /*--- Setting initial values for energy ---*/
  su2double e = e_start;

  su2double tolerance_P = 10; // Tolerance for pressure solution

  bool converged = false;         // Convergence flag
  unsigned long iter_max = 1000,  // Maximum number of iterations
                Iter = 0;       

  su2double delta_P,      // Pressure residual
            delta_e;      // Energy step size

  while(!converged && (Iter < iter_max)){

    /*--- Determine thermodynamic state based on current density and energy*/
    SetTDState_rhoe(rho ,e);

    /*--- Determine pressure and temperature residuals ---*/
    delta_P = Pressure - P;

    /*--- Continue iterative process if residuals are outside tolerances ---*/
    if((abs(delta_P) < tolerance_P)){
      converged = true;
    }else{
      
      /*--- Compute step size for energy ---*/
      delta_e = delta_P / dPde_rho;

      /*--- energy value ---*/
      e -= Newton_Relaxation * delta_e;
    }
    Iter ++;
  }

  /*--- Setting static energy as the converged value ---*/
  StaticEnergy = e;
}

void CDataDrivenFluid::SetTDState_rhoT(su2double rho, su2double T) {

  /*--- Setting initial values for density and energy ---*/
  su2double e = e_start;

  su2double tolerance_T = 1;  // Tolerance for temperature solution

  bool converged = false;         // Convergence flag
  unsigned long iter_max = 1000,  // Maximum number of iterations
                Iter = 0;       

  su2double delta_T,      // Temperature residual
            delta_e;  // Energy increment

  /*--- Initiating Newton solver ---*/
  while(!converged && (Iter < iter_max)){

    /*--- Determine thermodynamic state based on current density and energy*/
    SetTDState_rhoe(rho ,e);

    /*--- Determine temperature residual ---*/
    delta_T = Temperature - T;

    /*--- Continue iterative process if residuals are outside tolerances ---*/
    if((abs(delta_T) < tolerance_T)){
      converged = true;
    }else{
      

      delta_e = delta_T / dTde_rho;

      /*--- Update energy value ---*/
      e -= Newton_Relaxation * delta_e;
    }
    Iter ++;
  }

  /*--- Calculate thermodynamic state based on converged values for density and energy ---*/
  SetTDState_rhoe(rho, e);
}

void CDataDrivenFluid::SetTDState_hs(su2double h, su2double s) {

  /*--- Setting initial values for density and energy ---*/
  su2double rho = rho_start, 
            e = e_start;

  su2double tolerance_h = 10, // Tolerance for enthalpy solution
            tolerance_s = 1;  // Tolerance for entropy solution

  bool converged = false;         // Convergence flag
  unsigned long iter_max = 1000,  // Maximum number of iterations
                Iter = 0;       

  su2double delta_h,      // Enthalpy
            delta_s,      // Entropy residual
            delta_rho,    // Density step size
            delta_e,      // Energy step size
            determinant;  // Jacobian determinant

  /*--- Initiating Newton solver ---*/
  while(!converged && (Iter < iter_max)){

    /*--- Determine thermodynamic state based on current density and energy*/
    SetTDState_rhoe(rho ,e);

    su2double Enthalpy = e + Pressure / rho;
    /*--- Determine pressure and temperature residuals ---*/
    delta_h = Enthalpy - h;
    delta_s = Entropy - s;

    /*--- Continue iterative process if residuals are outside tolerances ---*/
    if((abs(delta_h) < tolerance_h) && (abs(delta_s) < tolerance_s)){
      converged = true;
    }else{
      
      su2double dh_de = 1 + dPde_rho / rho;
      su2double dh_drho = -Pressure*pow(rho, -2) + dPdrho_e / rho;

      determinant = dh_drho * dsde_rho - dh_de * dsdrho_e;
      delta_rho = (dsde_rho * delta_h - dh_de * delta_s) / determinant;
      delta_e = (-dsdrho_e * delta_h + dh_drho * delta_s) / determinant;

      /*--- Update density and energy values ---*/
      rho -= Newton_Relaxation * delta_rho;
      e -= Newton_Relaxation * delta_e;
    }
    Iter ++;
  }
  /*--- Calculate thermodynamic state based on converged values for density and energy ---*/
  SetTDState_rhoe(rho, e);
}

void CDataDrivenFluid::SetTDState_Ps(su2double P, su2double s) {
  
  /*--- Setting initial values for density and energy ---*/
  su2double rho = rho_start, 
            e = e_start;

  su2double tolerance_P = 10, // Tolerance for pressure solution
            tolerance_s = 1;  // Tolerance for entropy solution

  bool converged = false;         // Convergence flag
  unsigned long iter_max = 1000,  // Maximum number of iterations
                Iter = 0;       

  su2double delta_P,      // Enthalpy
            delta_s,      // Entropy residual
            delta_rho,    // Density step size
            delta_e,      // Energy step size
            determinant;  // Jacobian determinant

  /*--- Initiating Newton solver ---*/
  while(!converged && (Iter < iter_max)){

    /*--- Determine thermodynamic state based on current density and energy*/
    SetTDState_rhoe(rho ,e);

    /*--- Determine pressure and temperature residuals ---*/
    delta_P = Pressure - P;
    delta_s = Entropy - s;

    /*--- Continue iterative process if residuals are outside tolerances ---*/
    if((abs(delta_P) < tolerance_P) && (abs(delta_s) < tolerance_s)){
      converged = true;
    }else{

       /*--- Compute step size for density and energy ---*/
      determinant = dPdrho_e * dsde_rho - dPde_rho * dsdrho_e;

      delta_rho = (dsde_rho * delta_P - dPde_rho * delta_s) / determinant;
      delta_e = (-dsdrho_e * delta_P + dPdrho_e * delta_s) / determinant;

      /*--- Update density and energy values ---*/
      rho -= Newton_Relaxation * delta_rho;
      e -= Newton_Relaxation * delta_e;
    }
    Iter ++;
  }

  /*--- Calculate thermodynamic state based on converged values for density and energy ---*/
  SetTDState_rhoe(rho, e);
}