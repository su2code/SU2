/*!
 * \file variable_direct_2phase.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/variable_structure.hpp"

C2phaseVariable::C2phaseVariable(void) : CVariable() {

	unsigned short iVar;

	Source = 0; Enthalpy_Liquid = 0;
	Radius = 0; Liquid_Fraction = 0;

	Primitive_Liquid = NULL;

}

C2phaseVariable::C2phaseVariable(unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CVariable(val_nDim, val_nvar, config) {
  
  unsigned short iVar;

  /*--- Allocate space for the limiter ---*/
  
  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;
  
  Solution_Max = new su2double [nVar];
  Solution_Min = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }
  
  Primitive_Liquid = new su2double [11];
  for (iVar = 0; iVar < 11; iVar++)
	  Primitive_Liquid[iVar] = 0.0;

  /*--- Allocate residual structures ---*/
  Res_TruncError = new su2double [nVar];

  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }

}

/*
C2phaseVariable::SetSource(su2double S) {Source = S; };

C2phaseVariable::SetLiqEnthalpy(su2double h) {Enthalpy_Liquid = h; };

C2phaseVariable::SetRadius(su2double R)     {Radius = R;};

C2phaseVariable::SetLiquidFrac(su2double Y) {Liquid_Fraction = Y; };
*/


C2phaseVariable::~C2phaseVariable(void) {}



C2phase_HillVariable::C2phase_HillVariable(void) : C2phaseVariable() {

}

C2phase_HillVariable::C2phase_HillVariable(su2double val_R, su2double val_N, su2double rho_m, unsigned short nDim,
		unsigned short nVar, CConfig *config): C2phaseVariable(nDim, nVar, config) {

  unsigned short iVar;
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  /*--- Initialization of variables ---*/

  Solution[0] = rho_m * val_N;                Solution_Old[0] = rho_m * val_N;
  Solution[1] = rho_m * val_N*val_R;          Solution_Old[1] = rho_m * val_N*val_R;
  Solution[2] = rho_m * val_N*pow(val_R, 2);  Solution_Old[2] = rho_m * val_N*pow(val_R, 2);
  Solution[3] = rho_m * val_N*pow(val_R, 3);  Solution_Old[3] = rho_m * val_N*pow(val_R, 3);


  /*--- Allocate and initialize solution for the dual time strategy ---*/
  
  if (dual_time) {
	Solution_time_n[0] = rho_m * val_N;                Solution_time_n1[0] = rho_m * val_N;
	Solution_time_n[1] = rho_m * val_N*val_R;          Solution_time_n1[1] = rho_m * val_N*val_R;
	Solution_time_n[2] = rho_m * val_N*pow(val_R, 2);  Solution_time_n1[2] = rho_m * val_N*pow(val_R, 2);
	Solution_time_n[3] = rho_m * val_N*pow(val_R, 3);  Solution_time_n1[3] = rho_m * val_N*pow(val_R, 3);
  }
    
}

C2phase_HillVariable::~C2phase_HillVariable(void) {

	if (Primitive_Liquid != NULL) delete [] Primitive_Liquid;

}

su2double* C2phase_HillVariable::SetLiquidPrim(su2double *Primitive, su2double *Two_Phase_Var, su2double Rcritical, CFluidModel *FluidModel, CConfig *config) {

	su2double rho_l, rho_m, T_l, h_l, Psat, Tsat, sigma, Rc, Rdroplet, mom3;

	if (Two_Phase_Var[0] != 0 && Two_Phase_Var[3] != 0) Rdroplet = pow(Two_Phase_Var[3]/Two_Phase_Var[0], 1.0/3.0) * config->GetLength_Ref();
	else Rdroplet = 0;

    mom3 = Two_Phase_Var[3];
	P    = Primitive[nDim+1] * config->GetPressure_Ref();
	T    = Primitive[0]      * config->GetTemperature_Ref();
	rho  = Primitive[nDim+2] * config->GetDensity_Ref();
	h    = Primitive[nDim+3] * config->GetEnergy_Ref();
	Rc   = Rcritical         * config->GetLength_Ref();

	FluidModel->SetLiquidProp(P, T, rho, h, Rc, Rdroplet, mom3);

	T_l = FluidModel->GetLiquidTemperature()/config->GetTemperature_Ref();

	rho_l = FluidModel->GetLiquidDensity()/config->GetDensity_Ref();

	h_l = FluidModel->GetLiquidEnthalpy()/config->GetEnergy_Ref();

	Psat = FluidModel->GetPsat()/config->GetPressure_Ref();

	Tsat = FluidModel->GetTsat()/config->GetTemperature_Ref();

	sigma = FluidModel->GetSurfaceTension()/config->GetSurfTension_Ref();

	Rc = FluidModel->GetCriticalRadius()/config->GetLength_Ref();

	Rdroplet = Rdroplet/config->GetLength_Ref();

	rho_m = FluidModel->GetMixtureDensity()/config->GetDensity_Ref();


	Primitive_Liquid[0] = T_l;
	Primitive_Liquid[1] = rho_l;
	Primitive_Liquid[2] = h_l;
	Primitive_Liquid[3] = Psat;
	Primitive_Liquid[4] = Tsat;
	Primitive_Liquid[5] = sigma;
	Primitive_Liquid[6] = Rc;
	Primitive_Liquid[7] = Rdroplet;
	Primitive_Liquid[8] = rho_m;
	Primitive_Liquid[9] = 0.0;
	Primitive_Liquid[10] = 0.0;

	return Primitive_Liquid;

}



