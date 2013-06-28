/*!
 * \file variable_direct.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.
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

CPotentialVariable::CPotentialVariable(void) : CVariable() { }

CPotentialVariable::CPotentialVariable(double val_potential, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) 
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

	Residual_Old = new double [nVar];
	Residual_Sum = new double [nVar];

	/*--- Initialization of variables ---*/
	for (iVar = 0; iVar< nVar; iVar++) {
		Solution[iVar] = val_potential; 
		Solution_Old[iVar] = val_potential;
	}
	Charge_Density = new double [2];

	PlasmaRhoUGradient = new double* [3];
	for (iVar = 0; iVar < 3; iVar++)
		PlasmaRhoUGradient[iVar] = new double [nDim];

}

CPotentialVariable::~CPotentialVariable(void) { 

	delete [] Residual_Sum; delete [] Residual_Old;
	delete [] Charge_Density;
	for (unsigned short iVar = 0; iVar < 3; iVar++)
		delete [] PlasmaRhoUGradient[iVar];
	delete [] PlasmaRhoUGradient;

}

CEulerVariable::CEulerVariable(void) : CVariable() { }

CEulerVariable::CEulerVariable(double val_density, double *val_velocity, double val_energy, unsigned short val_ndim, 
		unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim;
	unsigned short iMesh, nMGSmooth = 0;

	bool Incompressible = config->GetIncompressible();
	bool FreeSurface = (config->GetKind_Solver() == FREE_SURFACE_EULER);
	bool AdjEuler = ((config->GetKind_Solver() == ADJ_EULER) || (config->GetKind_Solver() == ADJ_AEROACOUSTIC_EULER));
	bool Viscous = ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver() == RANS) || 
			(config->GetKind_Solver() == ADJ_NAVIER_STOKES) || (config->GetKind_Solver() == ADJ_NAVIER_STOKES) );

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; 
	Res_Visc = new double [nVar]; 
	Res_Sour = new double [nVar];
	Res_TruncError = new double [nVar];

	for (iVar = 0; iVar < nVar; iVar++) {
		Res_Conv[iVar] = 0.0;
		Res_Visc[iVar] = 0.0;
		Res_Sour[iVar] = 0.0;
		Res_TruncError[iVar] = 0.0;
	}

	/*--- Only for residual smoothing (multigrid) ---*/
	for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++)
		nMGSmooth += config->GetMG_CorrecSmooth(iMesh);

	if (nMGSmooth > 0) {
		Residual_Sum = new double [nVar]; 
		Residual_Old = new double [nVar];
	}

	/*--- Only for Runge-Kutta computations ---*/
	if (config->GetKind_TimeIntScheme_Flow() == RUNGE_KUTTA_EXPLICIT) {
		Res_Visc_RK = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			Res_Visc_RK[iVar] = new double [nVar];
	}

	/*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) Limiter = new double [nVar];

	/*--- Solution and old solution initialization ---*/
	if (Incompressible) {
		Solution[0] = config->GetPressure_FreeStreamND();
		Solution_Old[0] = config->GetPressure_FreeStreamND();
		for (iDim = 0; iDim < nDim; iDim++) {
			Solution[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
			Solution_Old[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
		}
	}
	else {
		Solution[0] = val_density;
		Solution_Old[0] = val_density;
		for (iDim = 0; iDim < nDim; iDim++) {
			Solution[iDim+1] = val_density*val_velocity[iDim];
			Solution_Old[iDim+1] = val_density*val_velocity[iDim];
		}
		Solution[nVar-1] = val_density*val_energy; 
		Solution_Old[nVar-1] = val_density*val_energy;
	}

	/*--- Allocate and initialize solution for dual time strategy ---*/
	if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
		if (Incompressible) {
			Solution_time_n[0] = config->GetPressure_FreeStreamND();
			Solution_time_n1[0] = config->GetPressure_FreeStreamND();
			for (iDim = 0; iDim < nDim; iDim++) {
				Solution_time_n[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
				Solution_time_n1[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
			}
		}
		else {
			Solution_time_n[0] = val_density;
			Solution_time_n1[0] = val_density;
			for (iDim = 0; iDim < nDim; iDim++) {
				Solution_time_n[iDim+1] = val_density*val_velocity[iDim];
				Solution_time_n1[iDim+1] = val_density*val_velocity[iDim];
			}
			Solution_time_n[nVar-1] = val_density*val_energy;
			Solution_time_n1[nVar-1] = val_density*val_energy;
		}
	}

	/*--- Allocate space for the time spectral source terms ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		TS_Source = new double[nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			TS_Source[iVar] = 0.0;
	}


	/*--- Allocate and initialize the primitive gradient variable (nVar+1) variables ---*/
	if (AdjEuler || Viscous) {
		Gradient_Primitive = new double* [nVar+1];
		for (iVar = 0; iVar < nVar+1; iVar++)
			Gradient_Primitive[iVar] = new double [nDim];
	}

	/*--- Allocate auxiliary vector for free surface source term ---*/
	if (FreeSurface) Grad_AuxVar = new double [nDim];

}

CEulerVariable::CEulerVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;
	unsigned short iMesh, nMGSmooth = 0;

	bool FreeSurface = (config->GetKind_Solver() == FREE_SURFACE_EULER);
	bool AdjEuler = ((config->GetKind_Solver() == ADJ_EULER) || (config->GetKind_Solver() == ADJ_AEROACOUSTIC_EULER));
	bool Viscous = ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver() == RANS) || 
			(config->GetKind_Solver() == ADJ_NAVIER_STOKES) || (config->GetKind_Solver() == ADJ_NAVIER_STOKES) );

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; 
	Res_Visc = new double [nVar]; 
	Res_Sour = new double [nVar];
	Res_TruncError = new double [nVar];

	for (iVar = 0; iVar < nVar; iVar++) {
		Res_Conv[iVar] = 0.0;
		Res_Visc[iVar] = 0.0;
		Res_Sour[iVar] = 0.0;
		Res_TruncError[iVar] = 0.0;
	}

	/*--- Only for residual smoothing (multigrid) ---*/
	for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++)
		nMGSmooth += config->GetMG_CorrecSmooth(iMesh);

	if (nMGSmooth > 0) {
		Residual_Sum = new double [nVar]; 
		Residual_Old = new double [nVar];
	}

	/*--- Only for Runge-Kutta computations ---*/
	if (config->GetKind_TimeIntScheme_Flow() == RUNGE_KUTTA_EXPLICIT) {
		Res_Visc_RK = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			Res_Visc_RK[iVar] = new double [nVar];
	}	

	/*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) Limiter = new double [nVar];

	/*--- Solution initialization ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}

	/*--- Allocate and initializate solution for dual time strategy ---*/
	if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || 
			(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
		Solution_time_n = new double [nVar];
		Solution_time_n1 = new double [nVar];

		for (iVar = 0; iVar < nVar; iVar++) {
			Solution_time_n[iVar] = val_solution[iVar];
			Solution_time_n1[iVar] = val_solution[iVar];
		}
	}

	/*--- Allocate space for the time spectral source terms ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		TS_Source = new double[nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			TS_Source[iVar] = 0.0;
	}

	/*--- Allocate and initializate the primitive gradient variable (nVar+1) variables ---*/
	if (AdjEuler || Viscous) {
		Gradient_Primitive = new double* [nVar+1];
		for (iVar = 0; iVar < nVar+1; iVar++)
			Gradient_Primitive[iVar] = new double [nDim];
	}

	/*--- Allocate auxiliar vector for free surface source term ---*/
	if (FreeSurface) Grad_AuxVar = new double [nDim];

}

CEulerVariable::~CEulerVariable(void) {
	unsigned short iVar;

	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour;
	delete [] Residual_Sum; delete [] Residual_Old;
	if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
	if (Limiter != NULL) delete [] Limiter;
	if (Solution_time_n != NULL) delete [] Solution_time_n;
	if (Solution_time_n1 != NULL) delete [] Solution_time_n1;
	delete [] Res_TruncError;

	if (Grad_AuxVar != NULL) delete [] Grad_AuxVar;

	if (TS_Source != NULL) delete [] TS_Source;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;

	for (iVar = 0; iVar < nVar+1; iVar++)
		delete [] Gradient_Primitive[iVar];
	delete [] Gradient_Primitive;
}

void CEulerVariable::SetGradient_PrimitiveZero(void) {
	unsigned short iVar, iDim;
	for (iVar = 0; iVar < nVar+1; iVar++)
		for (iDim = 0; iDim < nDim; iDim++)
			Gradient_Primitive[iVar][iDim] = 0.0;
}

double CEulerVariable::GetProjVel(double *val_vector) {
	double ProjVel;
	unsigned short iDim;

	ProjVel = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVel += Solution[iDim+1]*val_vector[iDim]/Solution[0];

	return ProjVel;
}

double CEulerVariable::GetProjVelInc(double *val_vector) {
	double ProjVel;
	unsigned short iDim;

	ProjVel = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVel += (Solution[iDim+1]/Solution[nDim+1])*val_vector[iDim];

	return ProjVel;
}

CNSVariable::CNSVariable(void) : CEulerVariable() { }

CNSVariable::CNSVariable(double val_density, double *val_velocity, double val_energy, 
		unsigned short val_ndim, unsigned short val_nvar,
		CConfig *config) : CEulerVariable(val_density, val_velocity, val_energy, val_ndim, val_nvar, config) {

	Temperature_Ref = config->GetTemperature_Ref();
	Viscosity_Ref   = config->GetViscosity_Ref();
	Viscosity_Inf   = config->GetViscosity_FreeStreamND();
	Prandtl_Lam     = config->GetPrandtl_Lam();
	Prandtl_Turb    = config->GetPrandtl_Turb();

}

CNSVariable::CNSVariable(double *val_solution, unsigned short val_ndim, 
		unsigned short val_nvar, CConfig *config) : CEulerVariable(val_solution, val_ndim, val_nvar, config) {

	Temperature_Ref = config->GetTemperature_Ref();
	Viscosity_Ref   = config->GetViscosity_Ref();
	Viscosity_Inf   = config->GetViscosity_FreeStreamND();
	Prandtl_Lam     = config->GetPrandtl_Lam();
	Prandtl_Turb    = config->GetPrandtl_Turb();
}

CNSVariable::~CNSVariable(void) {
	delete [] Res_Conv;
	delete [] Res_Visc;
	delete [] Res_Sour;
	delete [] Residual_Sum;
	delete [] Residual_Old;
	if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
	if (Limiter != NULL) delete [] Limiter;
	if (Solution_time_n != NULL) delete [] Solution_time_n;
	if (Solution_time_n1 != NULL) delete [] Solution_time_n1;	
	delete [] Res_TruncError;

	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;
}

void CNSVariable::SetLaminarViscosity(CConfig *config) {
	double Temperature_Dim;

	/*--- Calculate viscosity from a non-dim. Sutherland's Law ---*/
	Temperature_Dim = Temperature*Temperature_Ref;
	LaminarViscosity = 1.853E-5*(pow(Temperature_Dim/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_Dim+110.3));
	LaminarViscosity = LaminarViscosity/Viscosity_Ref;

}

void CNSVariable::SetEddyViscosity(unsigned short val_Kind_Turb_Model, CVariable *Turb_Solution) {

	switch (val_Kind_Turb_Model) {
		case NONE :
			EddyViscosity = 0.0;                     break;
		case SA   : case SST  :
			EddyViscosity = Turb_Solution->GetmuT(); break;
	}
}

void CNSVariable::SetVorticity(void) {
	double u_y = Gradient_Primitive[1][1];
	double v_x = Gradient_Primitive[2][0];
	double u_z = 0.0;
	double v_z = 0.0;
	double w_x = 0.0;
	double w_y = 0.0;

	if (nDim == 3) {
		u_z = Gradient_Primitive[1][2];
		v_z = Gradient_Primitive[2][2];
		w_x = Gradient_Primitive[3][0];
		w_y = Gradient_Primitive[3][1];
	}

	Vorticity[0] = w_y-v_z;
	Vorticity[1] = -(w_x-u_z);
	Vorticity[2] = v_x-u_y;

}

void CNSVariable::SetStrainMag(void) {
	double div = Gradient_Primitive[1][0] + Gradient_Primitive[2][1] + Gradient_Primitive[3][2];
	StrainMag = 0.0;
	// add diagonals
	StrainMag += pow(Gradient_Primitive[1][0] - 1.0/3.0*div,2.0);
	StrainMag += pow(Gradient_Primitive[2][1] - 1.0/3.0*div,2.0);
	StrainMag += pow(Gradient_Primitive[3][2] - 1.0/3.0*div,2.0);
	// add off diagonals
	StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][1]+Gradient_Primitive[2][0]),2.0);
	StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][2]+Gradient_Primitive[3][0]),2.0);
	StrainMag += 2.0*pow(0.5*(Gradient_Primitive[2][2]+Gradient_Primitive[3][1]),2.0);

	StrainMag = sqrt(2.0*StrainMag);
}

void CNSVariable::SetPressure(double Gamma, double *Coord, double turb_ke){
	Pressure_Old = Pressure;
	Pressure = (Gamma-1.0)*Solution[0]*(Solution[nVar-1]/Solution[0]-0.5*Velocity2 - turb_ke);
	if (Pressure < 0.0) {
		Pressure = Pressure_Old;
		cout <<"Negative pressure at point: ";
		for (unsigned short iDim = 0; iDim < nDim; iDim++) cout << Coord[iDim] <<" ";
		cout << endl;
		}
}

CTurbVariable::CTurbVariable(void) : CVariable() {}

CTurbVariable::~CTurbVariable(void) {
	if (Limiter != NULL) delete [] Limiter;
}

CTurbVariable::CTurbVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {

	if (config->GetKind_SlopeLimit() != NONE) Limiter = new double [nVar];
}

double CTurbVariable::GetmuT(){ return muT;}

void CTurbVariable::SetmuT(double val_muT){ muT = val_muT;}

CTurbSAVariable::CTurbSAVariable(void) : CTurbVariable() {}

CTurbSAVariable::~CTurbSAVariable(void) {
	if (Limiter != NULL) delete [] Limiter;
  if (TS_Source != NULL) delete [] TS_Source;
}

CTurbSAVariable::CTurbSAVariable(double val_nu_tilde, double val_muT, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_ndim, val_nvar,config) {

	/*--- Initialization of S-A variables ---*/
	Solution[0] = val_nu_tilde;		Solution_Old[0] = val_nu_tilde;

	/*--- Initialization of the eddy viscosity ---*/
	muT = val_muT;
  
  /*--- Allocate and initialize solution for the dual time strategy ---*/
	if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
			Solution_time_n[0]  = val_nu_tilde;
			Solution_time_n1[0] = val_nu_tilde;
  }
  
  /*--- Allocate space for the time spectral source terms ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		TS_Source = new double[nVar];
		for (unsigned short iVar = 0; iVar < nVar; iVar++)
			TS_Source[iVar] = 0.0;
	}
  
}

CTransLMVariable::CTransLMVariable(void) : CTurbVariable() {}

CTransLMVariable::~CTransLMVariable(void) {
	if (Limiter != NULL) delete [] Limiter;
}

CTransLMVariable::CTransLMVariable(double val_nu_tilde, double val_intermittency, double val_REth, double val_muT, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_ndim, val_nvar,config) {

	// Initialization of variables
	Solution[0] = val_nu_tilde;		   Solution_Old[0] = val_nu_tilde;
	Solution[1] = val_intermittency; Solution_Old[1] = val_intermittency;
	Solution[2] = val_REth;          Solution_Old[2] = val_REth;

	// Initialization of eddy viscosity
	muT = val_muT;
}

CTurbSSTVariable::CTurbSSTVariable(void) : CTurbVariable() {}

CTurbSSTVariable::~CTurbSSTVariable(void) {
	if (Limiter != NULL) delete [] Limiter;
}

CTurbSSTVariable::CTurbSSTVariable(double val_kine, double val_omega, double val_muT, unsigned short val_ndim, unsigned short val_nvar,
		double *constants, CConfig *config)
: CTurbVariable(val_ndim, val_nvar,config) {

	// Initialization of variables
	Solution[0] = val_kine;     Solution_Old[0] = val_kine;
	Solution[1] = val_omega;	Solution_Old[1] = val_omega;

	sigma_om2 = constants[3];
	beta_star = constants[6];

	F1   = 1.0;
	F2   = 0.0;
	CDkw = 0.0;

	// Initialization of eddy viscosity
	muT = val_muT;
}

void CTurbSSTVariable::SetBlendingFunc(double val_viscosity, double val_dist, double val_density){
	unsigned short iDim;
	double arg2, arg2A, arg2B, arg1;

	// Cross diffusion
	CDkw = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		CDkw += Gradient[0][iDim]*Gradient[1][iDim];
	CDkw *= 2.0*val_density*sigma_om2/Solution[1];
	CDkw = max(CDkw,pow(10.0,-20.0));

	// F1
	arg2A = sqrt(Solution[0])/(beta_star*Solution[1]*val_dist);
	arg2B = 500.0*val_viscosity/(val_density*val_dist*val_dist*Solution[1]);
	arg2 = max(arg2A,arg2B);
	arg1 = min(arg2,4.0*val_density*sigma_om2*Solution[0]/(CDkw*val_dist*val_dist));
	F1 = tanh(pow(arg1,4.0));

	// F2
	arg2 = max(2.0*arg2A,arg2B);
	F2 = tanh(pow(arg2,2.0));
}

double CTurbSSTVariable::GetF1blending(){
	return F1;
}

double CTurbSSTVariable::GetF2blending(){
	return F2;
}

double CTurbSSTVariable::GetCrossDiff(){
	return CDkw;
}

CPlasmaVariable::CPlasmaVariable(void) : CVariable() { }

CPlasmaVariable::CPlasmaVariable(double *val_density, double **val_velocity, double *val_energy, double *val_energy_vib,
		unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {

	unsigned short iVar, loc = 0, iSpecies, iDim;
	double rho, energy,	energy_vib, energy_el, enthalpy_formation, Gamma;	

	nMonatomics  = config->GetnMonatomics();
	nDiatomics = config->GetnDiatomics();
	nSpecies = config->GetnSpecies();
	Prandtl_Lam     = config->GetPrandtl_Lam();

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar];
	Res_Visc = new double [nVar];
	Res_Sour = new double [nVar];
	/*--- Allocate truncation error for multigrid strategy ---*/
	Res_TruncError = new double [nVar];

	for (iVar = 0; iVar < nVar; iVar++) {
		Res_Conv[iVar] = 0.0;
		Res_Visc[iVar] = 0.0;
		Res_Sour[iVar] = 0.0;
		Res_TruncError[iVar] = 0.0;
	}
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) 
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTERED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_UPWIND) Limiter = new double [nVar];

	/*--- Solution and old solution initialization ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

		Solution[loc + 0] = val_density[iSpecies];
		for (iDim = 0; iDim < nDim; iDim ++ )
			Solution[loc + iDim+1] = val_density[iSpecies]*val_velocity[iSpecies][iDim];
		Solution[loc + nDim + 1] = val_density[iSpecies]*val_energy[iSpecies];
		if ( iSpecies < nDiatomics )
			Solution[loc + nDim+2] = val_density[iSpecies]*val_energy_vib[iSpecies];

		Solution_Old[loc + 0] = val_density[iSpecies];
		for (iDim = 0; iDim < nDim; iDim ++ )
			Solution_Old[loc + iDim+1] = val_density[iSpecies]*val_velocity[iSpecies][iDim];
		Solution_Old[loc + nDim + 1] = val_density[iSpecies]*val_energy[iSpecies];
		if ( iSpecies < nDiatomics )
			Solution_Old[loc + nDim+2] = val_density[iSpecies]*val_energy_vib[iSpecies];

	}

	/*--- Allocate and initializate solution for dual time strategy ---*/
	if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
		/*		Solution_time_n = new double [nVar];
		 Solution_time_n1 = new double [nVar];

		 Solution_time_n[0] = val_density;
		 Solution_time_n1[0] = val_density;
		 for (iDim = 0; iDim < nDim; iDim++) {
		 Solution_time_n[iDim+1] = val_density*val_velocity[iDim];
		 Solution_time_n1[iDim+1] = val_density*val_velocity[iDim];
		 }
		 Solution_time_n[nVar-1] = val_density*val_energy;
		 Solution_time_n1[nVar-1] = val_density*val_energy;*/
	}

	/*--- Allocate and initialize the primitive gradient variable ---*/
	Gradient_Primitive = new double* [nVar+nSpecies];
	for (iVar = 0; iVar < nVar+nSpecies; iVar++)
		Gradient_Primitive[iVar] = new double [nDim];

	/*--- Allocate and initialize the local useful quantities ---*/
	Pressure   	= new double [nSpecies];
	Temperature	= new double [nSpecies];
	Temperature_tr = new double [nSpecies];
	Temperature_vib = new double [nSpecies];
	SoundSpeed  = new double [nSpecies];
	Velocity    = new double*[nSpecies];
	Velocity2   = new double [nSpecies];
	Enthalpy    = new double [nSpecies];
	Max_Lambda_Inv_MultiSpecies = new double [nSpecies];
	Max_Lambda_Visc_MultiSpecies = new double [nSpecies];
	Lambda = new double [nSpecies];
	Species_Delta_Time = new double [nSpecies];
	Sensor_MultiSpecies = new double [nSpecies];
	LaminarViscosity_MultiSpecies = new double [nSpecies];
	EddyViscosity_MultiSpecies = new double [nSpecies];
	ThermalCoeff = new double [nSpecies];

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		rho = Solution[loc + 0];
		energy = Solution[loc + nDim+1];
		energy_vib = 0.0;
		if ( iSpecies < nDiatomics ) energy_vib = Solution[loc + nDim+2];
		energy_el = 0.0;
		Velocity2[iSpecies] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity2[iSpecies] += val_velocity[iSpecies][iDim]*val_velocity[iSpecies][iDim];
		}
		Gamma = config->GetSpecies_Gamma(iSpecies);
		enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Pressure[iSpecies] 		= (Gamma - 1.0)*rho*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - enthalpy_formation);
		SoundSpeed[iSpecies] 	= sqrt(fabs(Gamma*(Gamma - 1.0)*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - enthalpy_formation)));		
		Enthalpy[iSpecies]		= energy + Pressure[iSpecies]/rho;
		Max_Lambda_Inv_MultiSpecies[iSpecies] = 0.0;
		Max_Lambda_Visc_MultiSpecies[iSpecies] = 0.0;
		Lambda[iSpecies] = 0.0;
		LaminarViscosity_MultiSpecies[iSpecies] = 0.0;
		EddyViscosity_MultiSpecies[iSpecies] = 0.0;
		Species_Delta_Time[iSpecies] = 1.0;
	}

	Elec_Field = new double [nDim];
	for(iDim =0; iDim < nDim; iDim ++)
		Elec_Field[iDim] = 0.0;
}

CPlasmaVariable::CPlasmaVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, loc = 0, iSpecies, iDim;

	double rho, energy, energy_vib, energy_el, Gamma, enthalpy_formation;

	nMonatomics  = config->GetnMonatomics();
	nDiatomics = config->GetnDiatomics();
	nSpecies = config->GetnSpecies();
	Prandtl_Lam     = config->GetPrandtl_Lam();

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar];
	Res_Visc = new double [nVar];
	Res_Sour = new double [nVar];
	Res_TruncError = new double [nVar];

	for (iVar = 0; iVar < nVar; iVar++) {
		Res_Conv[iVar] = 0.0;
		Res_Visc[iVar] = 0.0;
		Res_Sour[iVar] = 0.0;
		Res_TruncError[iVar] = 0.0;
	}
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTERED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_UPWIND) Limiter = new double [nVar];

	/*--- Allocate truncation error for multigrid strategy ---*/
	Res_TruncError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_TruncError[iVar] = 0.0;


	/*--- Solution initialization ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}

	/*--- Allocate and initializate the primitive gradient variable ---*/
	Gradient_Primitive = new double* [nVar+nSpecies];
	for (iVar = 0; iVar < nVar+nSpecies; iVar++)
		Gradient_Primitive[iVar] = new double [nDim];

	Pressure   	= new double [nSpecies];
	Temperature	= new double [nSpecies];
	Temperature_tr = new double [nSpecies];
	Temperature_vib = new double [nSpecies];
	SoundSpeed  = new double [nSpecies];
	Velocity    = new double*[nSpecies];
	Velocity2   = new double [nSpecies];
	Enthalpy    = new double [nSpecies];
	Max_Lambda_Inv_MultiSpecies = new double [nSpecies];
	Max_Lambda_Visc_MultiSpecies = new double [nSpecies];
	Lambda = new double [nSpecies];
	Species_Delta_Time = new double [nSpecies];
	Sensor_MultiSpecies = new double [nSpecies];

	LaminarViscosity_MultiSpecies = new double [nSpecies];
	EddyViscosity_MultiSpecies = new double [nSpecies];
	ThermalCoeff = new double [nSpecies];

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		rho = Solution[loc + 0];
		energy = Solution[loc + nDim+1];
		energy_vib = 0.0;
		if ( iSpecies < nDiatomics ) energy_vib = Solution[loc + nDim+2];
		energy_el = 0.0;
		Velocity[iSpecies]  = new double [nDim];
		Velocity2[iSpecies] = 0.0;
		for (iDim = 0; iDim < nDim; iDim ++) {
			Velocity2[iSpecies] 		    += pow((Solution[loc + iDim + 1]/rho),2);
			Velocity[iSpecies][iDim] 	 = Solution[loc + iDim + 1]/rho;
		}
		Gamma = config->GetSpecies_Gamma(iSpecies);
		enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Pressure[iSpecies] 		= (Gamma - 1.0)*rho*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - enthalpy_formation);
		SoundSpeed[iSpecies] 	= sqrt(fabs(Gamma*(Gamma - 1.0)*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - enthalpy_formation)));
		Enthalpy[iSpecies]		= (energy - 0.5*rho*Velocity2[iSpecies] - energy_vib - rho*enthalpy_formation + Pressure[iSpecies])/rho  + enthalpy_formation;
		Max_Lambda_Inv_MultiSpecies[iSpecies] = 0.0;
		Max_Lambda_Visc_MultiSpecies[iSpecies] = 0.0;
		Lambda[iSpecies] = 0.0;
		LaminarViscosity_MultiSpecies[iSpecies] = 0.0;
		EddyViscosity_MultiSpecies[iSpecies] = 0.0;
		Species_Delta_Time[iSpecies] = 1.0;
	}
	Elec_Field = new double [nDim];
	for(iDim =0; iDim < nDim; iDim ++)
		Elec_Field[iDim] = 0.0;
}

CPlasmaVariable::~CPlasmaVariable(void) {
	unsigned short iVar;

	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour; delete [] Res_TruncError;
	delete [] Residual_Sum; delete [] Residual_Old;
	if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
	if (Limiter != NULL) delete [] Limiter;
	if (Solution_time_n != NULL) delete [] Solution_time_n;
	if (Solution_time_n1 != NULL) delete [] Solution_time_n1;

	delete [] Max_Lambda_Inv_MultiSpecies;
	delete [] Max_Lambda_Visc_MultiSpecies;
	delete [] Lambda;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;

	for (iVar = 0; iVar < nVar+nSpecies; iVar++)
		delete [] Gradient_Primitive[iVar];
	delete [] Gradient_Primitive;


	for (iVar = 0; iVar < nSpecies; iVar++) {
		delete [] Velocity[iVar];
	}
	delete [] Velocity;
	delete [] Pressure;
	delete [] Velocity2;
	delete [] Temperature;
	delete [] Temperature_tr;
	delete [] Temperature_vib;
	delete [] SoundSpeed;
	delete [] Enthalpy;

	delete [] Max_Lambda_Inv_MultiSpecies;
	delete [] Max_Lambda_Visc_MultiSpecies; 
	delete [] Lambda;
	delete [] Species_Delta_Time; 
	delete [] Sensor_MultiSpecies;
	delete [] LaminarViscosity_MultiSpecies;
	delete [] EddyViscosity_MultiSpecies; 
	delete [] ThermalCoeff;
}

void CPlasmaVariable::SetGradient_PrimitiveZero(void) {
	unsigned short iVar, iDim;
	for (iVar = 0; iVar < nVar+nSpecies; iVar++)
		for (iDim = 0; iDim < nDim; iDim++)
			Gradient_Primitive[iVar][iDim] = 0.0;
}

double CPlasmaVariable::GetProjVel(double *val_vector, unsigned short val_Species) {
	double ProjVel;
	unsigned short iDim, loc = 0;

	if ( val_Species < nDiatomics ) loc = (nDim+3)*val_Species;
	else loc = (nDim+3)*nDiatomics + (nDim+2)*(val_Species-nDiatomics);

	ProjVel = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVel += Solution[loc + iDim+1]*val_vector[iDim]/Solution[loc + 0];
	return ProjVel;
}

void CPlasmaVariable::SetVelocity2() {

	unsigned short loc = 0, iSpecies, iDim;
	double rho;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {

		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

		rho = Solution[loc + 0];
		Velocity2[iSpecies] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Velocity2[iSpecies] += pow(Solution[loc+iDim+1]/rho, 2);
	}
}

void CPlasmaVariable::SetSoundSpeed(CConfig *config) {

	unsigned short loc = 0, iSpecies;
	double Gamma, Enthalpy_formation;
	double Density, Energy, Energy_vib, Energy_el;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

		Density = Solution[loc+0];
		Energy = Solution[loc+nDim+1] / Density;
		Energy_vib = 0.0;
		Energy_el = 0.0;
		if (iSpecies < nDiatomics )
			Energy_vib = Solution[loc+nDim+2] / Density;
		Gamma = config->GetSpecies_Gamma(iSpecies);
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		SoundSpeed[iSpecies] 	= sqrt(Gamma*(Gamma - 1.0)*(Energy - 0.5*Velocity2[iSpecies] - Energy_vib - Energy_el - Enthalpy_formation));
	}
}

void CPlasmaVariable::SetPressure(CConfig *config, double *Coord) {

	unsigned short loc = 0, iSpecies, iDim;
	double Gamma, Enthalpy_formation;
	double Density, Energy, Energy_vib, Energy_el;
	double Vel2;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

		Density = Solution[loc+0];
		Energy  = Solution[loc+nDim+1] / Density;
		Energy_vib = 0.0;
		if ( iSpecies < nDiatomics )
			Energy_vib = Solution[loc+nDim+2] / Density;
		Energy_el = 0.0;
		Vel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Vel2 += Solution[loc+iDim+1]/Solution[loc+0] * Solution[loc+iDim+1]/Solution[loc+0];
		Gamma = config->GetSpecies_Gamma(iSpecies);
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Pressure[iSpecies] = (Gamma-1.0) * Density * (Energy - 1.0/2.0*Vel2 - Enthalpy_formation - Energy_vib - Energy_el);

		if (Pressure[iSpecies] < 0.0) {
			cout << "Density = " << Density << endl;
			cout << "energy = " << Energy << endl;
			cout << "vel2  = " << Vel2 << endl;

			cout << "Negative pressure of " << Pressure[iSpecies] << " at point: ";
			for (unsigned short iDim = 0; iDim < nDim; iDim++) cout << Coord[iDim] <<" ";
			cout << " in species " << iSpecies;
			cout << endl;
		}

	}
}

void CPlasmaVariable::SetTemperature_TR(CConfig *config, double *Coord) {

	unsigned short loc = 0, iSpecies, iDim;
	double Gamma, Enthalpy_formation, Gas_constant;
	double Density, Energy, Energy_vib, Energy_el;
	double Vel2;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);			
		Density    = Solution[loc+0];
		Energy     = Solution[loc+nDim+1] / Density;
		Energy_vib = 0.0;
		Energy_el  = 0.0;
		Vel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Vel2 += Solution[loc+iDim+1] / Solution[loc+0] * Solution[loc+iDim+1] / Solution[loc+0];
		if (iSpecies < nDiatomics) {
			Energy_vib = Solution[loc+nDim+2] / Density;
		}
		Gamma = config->GetSpecies_Gamma(iSpecies);
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);

		Temperature_tr[iSpecies] = (Gamma-1.0)/Gas_constant 
				* (Energy - 1.0/2.0*Vel2 - Enthalpy_formation - Energy_vib - Energy_el);
		Temperature[iSpecies] = Temperature_tr[iSpecies];
		if (Temperature_tr[iSpecies] < 0.0) {
			cout << "Negative translational-rotational temperature at point: ";
			for (unsigned short iDim = 0; iDim < nDim; iDim++) cout << Coord[iDim] <<" ";
			cout << " in species " << iSpecies;
			cout << endl;
		}
	}
}

void CPlasmaVariable::SetTemperature_vib(CConfig *config, double *Coord) {
	//Based on energy in a harmonic oscillator.  Used by Candler (1988) and Scalabrin (2007).
	unsigned short loc = 0, iSpecies;
	double Density, Energy_vib, CharVibTemp, MolarMass;

	if (nDiatomics == 0) return;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

		if (iSpecies < nDiatomics) {
			Density    = Solution[loc+0];
			Energy_vib = Solution[loc+nDim+2] / Density;		
			CharVibTemp = config->GetCharVibTemp(iSpecies);
			MolarMass = config->GetMolar_Mass(iSpecies);

			Temperature_vib[iSpecies] = CharVibTemp / log(CharVibTemp*UNIVERSAL_GAS_CONSTANT/(Energy_vib*MolarMass) + 1.0);			
			if (Temperature_vib[iSpecies] < 0.0) {
				cout << "Negative vibrational temperature at point: ";
				for (unsigned short iDim = 0; iDim < nDim; iDim++) cout << Coord[iDim] <<" ";
				cout << " in species " << iSpecies;
				cout << endl;
			}
			if (Temperature_vib[iSpecies] != Temperature_vib[iSpecies]) {
				cout << "NaN vibrational temperature at point: ";				
				for (unsigned short iDim = 0; iDim < nDim; iDim++) cout << Coord[iDim] <<" ";
				cout << " in species " << iSpecies;
				cout << endl;
				cout << "log term: " << CharVibTemp*UNIVERSAL_GAS_CONSTANT/(Energy_vib*MolarMass) + 1.0 << endl;
				cout << "energy vib: " << Energy_vib << endl;

			}
		}
		else {
			Temperature_vib[iSpecies] = 0.0;
		}
	}
}


void CPlasmaVariable::SetEnthalpy() {
	unsigned short loc = 0, iSpecies;
	double Density, Energy;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		Density    = Solution[loc+0];
		Energy     = Solution[loc+nDim+1] / Density;
		Enthalpy[iSpecies] = Energy + Pressure[iSpecies] / Density;
	}
}

double CPlasmaVariable::GetVelocity(unsigned short val_dim,unsigned short val_Species) {

	unsigned short loc = 0;
	double rho, rho_vel;

	if ( val_Species < nDiatomics ) loc = (nDim+3)*val_Species;
	else loc = (nDim+3)*nDiatomics + (nDim+2)*(val_Species-nDiatomics);
	rho = Solution[loc + 0];
	rho_vel = Solution[loc + val_dim + 1];
	return (rho_vel/rho);
}

void CPlasmaVariable::SetLaminarViscosity(CConfig *config) {
	double Temperature_Dim;
	unsigned short iSpecies;
	double Temperature_Ref, Viscosity_Ref;

  switch (config->GetKind_GasModel()) {
    case ARGON:
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
        Temperature_Ref	= config->GetTemperature_Ref(iSpecies);
        Viscosity_Ref = config->GetViscosity_Ref(iSpecies);
        Temperature_Dim = Temperature[iSpecies]*Temperature_Ref;
        LaminarViscosity_MultiSpecies[iSpecies] = 2.3E-5*(pow(Temperature_Dim/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_Dim+110.3));
        LaminarViscosity_MultiSpecies[iSpecies] = LaminarViscosity_MultiSpecies[iSpecies]/Viscosity_Ref;
      }
      iSpecies = nSpecies - 1;
      Temperature_Ref	= config->GetTemperature_Ref(iSpecies);
      Viscosity_Ref = config->GetViscosity_Ref(iSpecies);
      Temperature_Dim = Temperature[iSpecies]*Temperature_Ref;
      LaminarViscosity_MultiSpecies[iSpecies] = 2.3E-5*1.3E-7*(pow(Temperature_Dim/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_Dim+110.3));
      LaminarViscosity_MultiSpecies[iSpecies] = LaminarViscosity_MultiSpecies[iSpecies]/Viscosity_Ref;
      //LaminarViscosity_MultiSpecies[iSpecies-1] = 0.0*LaminarViscosity_MultiSpecies[iSpecies]/Viscosity_Ref;
      break;
      
    case N2:
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        
        
        
      }
      break;
  }
}

void CPlasmaVariable ::SetThermalCoeff(CConfig *config) {
	unsigned short iSpecies;
	double Gamma, Gas_constant;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		Gamma = config->GetSpecies_Gamma(iSpecies);
		Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);
		ThermalCoeff[iSpecies] = Gamma * LaminarViscosity_MultiSpecies[iSpecies] / ( (Gamma-1.0)*Gas_constant*Prandtl_Lam );
	}
}

CLevelSetVariable::CLevelSetVariable(void) : CVariable() {}

CLevelSetVariable::CLevelSetVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate limiter (upwind)---*/
	if (config->GetKind_SlopeLimit() != NONE) Limiter = new double [nVar];

}

CLevelSetVariable::CLevelSetVariable(double val_levelset, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar,config) {
	unsigned short iVar;

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate limiter (upwind)---*/
	if (config->GetKind_SlopeLimit() != NONE) Limiter = new double [nVar];

	/*--- Solution and old solution initialization ---*/
	Solution[0] = val_levelset;		Solution_Old[0] = val_levelset;

	if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) 
			|| (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
		Solution_time_n[0] = val_levelset;
		Solution_time_n1[0] = val_levelset;
	}

}

CLevelSetVariable::~CLevelSetVariable(void) {
	unsigned short iVar;

	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour;
	delete [] Residual_Sum; delete [] Residual_Old;
	if (Limiter != NULL) delete [] Limiter;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;
}

CWaveVariable::CWaveVariable(void) : CVariable() { }

CWaveVariable::CWaveVariable(double *val_wave, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) 
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];

	/*--- Allocate direct solution container for adjoint problem ---*/
	Solution_Direct = new double[nVar];

	/*--- Allocate aux gradient vector ---*/
	Grad_AuxVar = new double [nDim];

	/*--- Initialization of variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_wave[iVar]; 
		Solution_Old[iVar] = val_wave[iVar];
		Solution_Direct[iVar] = 0.0;
	}

	/*--- Initialize sources to zero ---*/
	Thickness_Noise  = 0.0;
	Loading_Noise    = 0.0;
	Quadrupole_Noise = 0.0;

}

CWaveVariable::~CWaveVariable(void) { 

	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour;
	delete [] Residual_Sum; delete [] Residual_Old;

}

CFEAVariable::CFEAVariable(void) : CVariable() { }

CFEAVariable::CFEAVariable(double *val_fea, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) 
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];

	/*--- Initialization of variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_fea[iVar]; 
		Solution_Old[iVar] = val_fea[iVar];
	}
}

CFEAVariable::~CFEAVariable(void) { 

	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour;
	delete [] Residual_Sum; delete [] Residual_Old;

}

CHeatVariable::CHeatVariable(void) : CVariable() { }

CHeatVariable::CHeatVariable(double *val_heat, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) 
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;
	
	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	
	/*--- Allocate direct solution container for adjoint problem ---*/
	Solution_Direct = new double[nVar];
	
	/*--- Allocate aux gradient vector ---*/
	Grad_AuxVar = new double [nDim];
	
	/*--- Initialization of variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_heat[iVar]; 
		Solution_Old[iVar] = val_heat[iVar];
		Solution_Direct[iVar] = 0.0;
	}
	
}

CHeatVariable::~CHeatVariable(void) { 
	
	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour;
	delete [] Residual_Sum; delete [] Residual_Old;
	
}

