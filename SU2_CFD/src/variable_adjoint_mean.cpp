/*!
 * \file variable_adjoint_mean.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/variable_structure.hpp"

CAdjEulerVariable::CAdjEulerVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	Psi = NULL;
	ForceProj_Vector = NULL;
	ObjFuncSource = NULL;
	IntBoundary_Jump = NULL;
	TS_Source = NULL;
  
}

CAdjEulerVariable::CAdjEulerVariable(double val_psirho, double *val_phi, double val_psie, unsigned short val_ndim,
																		 unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
  bool Incompressible = config->GetIncompressible();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Array initialization ---*/
	Psi = NULL;
	ForceProj_Vector = NULL;
	ObjFuncSource = NULL;
	IntBoundary_Jump = NULL;
	TS_Source = NULL;
  
	/*--- Allocate residual structures ---*/
  Res_TruncError = new double [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
		Res_TruncError[iVar] = 0.0;
	}
  
  /*--- Only for residual smoothing (multigrid) ---*/
	for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++)
		nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
  
  if (nMGSmooth > 0) {
    Residual_Sum = new double [nVar];
    Residual_Old = new double [nVar];
  }
	
	/*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED)
		Undivided_Laplacian = new double [nVar];
	if ((config->GetKind_ConvNumScheme_AdjFlow() == SPACE_UPWIND) &&
			(config->GetKind_SlopeLimit_AdjFlow() != NONE)) {
		Limiter = new double [nVar];
		Solution_Max = new double [nVar];
		Solution_Min = new double [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Limiter[iVar] = 0.0;
			Solution_Max[iVar] = 0.0;
			Solution_Min[iVar] = 0.0;
		}
	}
  
  /*--- Allocate and initialize solution ---*/
	if (Incompressible) {
		Solution[0] = 0.0; 	Solution_Old[0] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Solution[iDim+1] = 0.0;
			Solution_Old[iDim+1] = 0.0;
		}
	}
	else {
		Solution[0] = val_psirho; 	Solution_Old[0] = val_psirho;
		Solution[nVar-1] = val_psie; Solution_Old[nVar-1] = val_psie;
		for (iDim = 0; iDim < nDim; iDim++) {
			Solution[iDim+1] = val_phi[iDim];
			Solution_Old[iDim+1] = val_phi[iDim];
		}
	}
  
  /*--- Allocate and initialize solution for dual time strategy ---*/
	if (dual_time) {
		if (Incompressible) {
			Solution_time_n[0] = 0.0;
			Solution_time_n1[0] = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Solution_time_n[iDim+1] = 0.0;
				Solution_time_n1[iDim+1] = 0.0;
			}
		}
		else {
			Solution_time_n[0] = val_psirho;
			Solution_time_n1[0] = val_psirho;
			for (iDim = 0; iDim < nDim; iDim++) {
				Solution_time_n[iDim+1] = val_phi[iDim];
				Solution_time_n1[iDim+1] = val_phi[iDim];
			}
			Solution_time_n[nVar-1] = val_psie;
			Solution_time_n1[nVar-1] = val_psie;
		}
	}
  
  /*--- Allocate auxiliar vector for sensitivity computation ---*/
	Grad_AuxVar = new double [nDim];
	
	/*--- Allocate and initialize projection vector for wall boundary condition ---*/
	ForceProj_Vector = new double [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		ForceProj_Vector[iDim] = 0.0;
  
  /*--- Allocate and initialize interior boundary jump vector for near field boundary condition ---*/
	IntBoundary_Jump = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		IntBoundary_Jump[iVar] = 0.0;
  
	/*--- Allocate and initialize vector containing objective function sensitivity for discrete adjoint ---*/
	if (config->GetKind_Adjoint() == DISCRETE) {
		ObjFuncSource = new double [nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			ObjFuncSource[iVar] = 0.0;
	}
  
  /*--- Allocate space for the time spectral source terms ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		TS_Source = new double[nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			TS_Source[iVar] = 0.0;
	}
	
}

CAdjEulerVariable::CAdjEulerVariable(double *val_solution, unsigned short val_ndim,
                                     unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Array initialization ---*/
	Psi = NULL;
	ForceProj_Vector = NULL;
	ObjFuncSource = NULL;
	IntBoundary_Jump = NULL;
	TS_Source = NULL;
  
	/*--- Allocate residual structures ---*/
  Res_TruncError = new double [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
		Res_TruncError[iVar] = 0.0;
	}
  
  /*--- Only for residual smoothing (multigrid) ---*/
	for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++)
		nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
  
  if (nMGSmooth > 0) {
    Residual_Sum = new double [nVar];
    Residual_Old = new double [nVar];
  }
  
	/*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED)
		Undivided_Laplacian = new double [nVar];
	if ((config->GetKind_ConvNumScheme_AdjFlow() == SPACE_UPWIND) &&
			(config->GetKind_SlopeLimit_AdjFlow() != NONE)) {
		Limiter = new double [nVar];
		Solution_Max = new double [nVar];
		Solution_Min = new double [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Limiter[iVar] = 0.0;
			Solution_Max[iVar] = 0.0;
			Solution_Min[iVar] = 0.0;
		}
	}
  
	/*--- Solution initialization ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}
  
	/*--- Allocate and initializate solution for dual time strategy ---*/
	if (dual_time) {
		Solution_time_n = new double [nVar];
		Solution_time_n1 = new double [nVar];
    
		for (iVar = 0; iVar < nVar; iVar++) {
			Solution_time_n[iVar] = val_solution[iVar];
			Solution_time_n1[iVar] = val_solution[iVar];
		}
	}
  
  /*--- Allocate auxiliar vector for sensitivity computation ---*/
	Grad_AuxVar = new double [nDim];
	
	/*--- Allocate and initializate projection vector for wall boundary condition ---*/
	ForceProj_Vector = new double [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		ForceProj_Vector[iDim] = 0.0;
  
  /*--- Allocate and initializate interior boundary jump vector for near field boundary condition ---*/
	IntBoundary_Jump = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		IntBoundary_Jump[iVar] = 0.0;
	
	/*--- Allocate and initialize vector containing objective function sensitivity for discrete adjoint ---*/
	if (config->GetKind_Adjoint() == DISCRETE) {
		ObjFuncSource = new double [nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			ObjFuncSource[iVar] = 0.0;
	}
  
	/*--- Allocate space for the time spectral source terms ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		TS_Source = new double[nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			TS_Source[iVar] = 0.0;
	}
  
}

CAdjEulerVariable::~CAdjEulerVariable(void) {
	
	if (Psi               != NULL) delete [] Psi;
	if (ForceProj_Vector  != NULL) delete [] ForceProj_Vector;
	if (ObjFuncSource     != NULL) delete [] ObjFuncSource;
	if (IntBoundary_Jump  != NULL) delete [] IntBoundary_Jump;
	if (TS_Source         != NULL) delete [] TS_Source;
  
}

void CAdjEulerVariable::SetTheta(double val_density, double *val_velocity, double val_enthalpy) {
	unsigned short iDim;
	
	Theta = val_density*Solution[0];
	Theta += val_density*val_enthalpy*Solution[nDim+1];
	
	for (iDim = 0; iDim < nDim; iDim++)
		Theta += val_density*val_velocity[iDim]*Solution[iDim+1];
}

void CAdjEulerVariable::SetPrimVar_Compressible(double val_adjlimit) {
	unsigned short iVar;
  bool check_dens = false;
  
  check_dens = (fabs(Solution[0]) > val_adjlimit);             // Check the adjoint density
  
  /*--- Check that the solution has a physical meaning ---*/
  if (check_dens) {
    
    /*--- Copy the old solution ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    
  }
  
}

void CAdjEulerVariable::SetPrimVar_Incompressible(double val_adjlimit) {
  unsigned short iVar;
  bool check_press = false;
  
  check_press = (fabs(Solution[0]) > val_adjlimit);             // Check the adjoint pressure
  
  /*--- Check that the solution has a physical meaning ---*/
  if (check_press) {
    
    /*--- Copy the old solution ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    
  }
  
}

CAdjNSVariable::CAdjNSVariable(void) : CAdjEulerVariable() { }

CAdjNSVariable::CAdjNSVariable(double *val_solution, unsigned short val_ndim,
                               unsigned short val_nvar, CConfig *config) : CAdjEulerVariable(val_solution, val_ndim, val_nvar, config) {
  
}

CAdjNSVariable::CAdjNSVariable(double val_psirho, double *val_phi, double val_psie,
                               unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CAdjEulerVariable(val_psirho, val_phi, val_psie, val_ndim, val_nvar, config) {
  
  kappapsi_Volume = 0.0;
  
}

CAdjNSVariable::~CAdjNSVariable(void) { }


void SetLaminarViscosity_Jacobian(CConfig *config) {
	// Derivative of Sutherland's Law wrt conservative variables
  
  //	double Temperature, Pressure, Temperature_Dim;
  //	double *U, T_Sens, P_Sens, TDim_Sens;
  //
  //	U = new double[nVar];
  //	T_Sens = new double[nVar];
  //	P_Sens = new double[nVar];
  //	TDim_Sens = new double[nVar];
  //
  //	for (iVar = 0; iVar < nVar; iVar++)
  //		U[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
  //
  //	if (iDim == 2)
  //		Pressure = Gamma_Minus_One*(U[3] - (U[1]*U[1] + U[2]*U[2])/(2*U[0]));
  //	else
  //		Pressure = Gamma_Minus_One*(U[4] - (U[1]*U[1] + U[2]*U[2] + U[3]*U[3])/(2*U[0]));
  //
  //	Temperature = Pressure/(U[0]*Gas_Constant);
  //
  //	if (iDim == 2) {
  //		P_Sens[0] = Gamma_Minus_One*(U[1]*U[1] + U[2]*U[2])/(2*U[0]*U[0]);
  //		P_Sens[1] = -Gamma_Minus_One*U[1]/U[0];
  //		P_Sens[2] = -Gamma_Minus_One*U[2]/U[0];
  //		P_Sens[3] = Gamma_Minus_One;
  //
  //	} else {
  //		P_Sens[0] = Gamma_Minus_One*(U[1]*U[1] + U[2]*U[2] + U[3]*U[3])/(2*U[0]*U[0]);
  //		P_Sens[1] = -Gamma_Minus_One*U[1]/U[0];
  //		P_Sens[2] = -Gamma_Minus_One*U[2]/U[0];
  //		P_Sens[3] = -Gamma_Minus_One*U[3]/U[0];
  //		P_Sens[4] = Gamma_Minus_One;
  //	}
  //
  //	T_Sens[0] = P_Sens[0]/(U[0]*Gas_Constant) - Pressure/(U[0]*U[0]*Gas_Constant);
  //	for (iVar = 1; iVar < nVar; iVar++)
  //		T_Sens[iVar] = P_Sens[iVar]/(U[0]*Gas_Constant);
  //
  //	Temperature_Dim = Temperature*Temperature_Ref;
  //
  //	for (iVar = 0; iVar < nVar; iVar++)
  //		TDim_Sens[iVar] = T_Sens[iVar]*Temperature_Ref;
  //
  //	for (iVar = 0; iVar < nVar; iVar++) {
  //		LaminarViscosity_Jacobian[iVar] = 1.853E-5*(
  //				(3.0/2.0)*(TDim_Sens[iVar]/300.0)*sqrt(Temperature_Dim/300.0) * (300.0+110.3)/(Temperature_Dim+110.3)
  //				- pow(Temperature_Dim/300.0,3.0/2.0) * TDim_Sens[iVar]* (300.0+110.3)/((Temperature_Dim+110.3)*(Temperature_Dim+110.3))
  //				);
  //
  //		LaminarViscosity_Jacobian[iVar] = LaminarViscosity_Jacobian[iVar]/Viscosity_Ref;
  //	}
  
}

void SetEddyViscosity_Jacobian(unsigned short val_Kind_Turb_Model, CVariable *TurbVariable) { }

void CAdjNSVariable::SetTheta(double val_density, double *val_velocity, double val_enthalpy) {
	Theta = val_density*Solution[0];
	Theta += val_density*val_enthalpy*Solution[nDim+1];
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Theta += val_density*val_velocity[iDim]*Solution[iDim+1];
}