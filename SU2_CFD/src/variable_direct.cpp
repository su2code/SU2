/*!
 * \file variable_direct.cpp
 * \brief Definition of the solution fields.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
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

	// Initialization of variables
	for (iVar = 0; iVar< nVar; iVar++) {
		Solution[iVar] = val_potential; 
		Solution_Old[iVar] = val_potential;
	}
}

CPotentialVariable::~CPotentialVariable(void) { 

	delete [] Residual_Sum; delete [] Residual_Old;

}

CEulerVariable::CEulerVariable(void) : CVariable() { }

CEulerVariable::CEulerVariable(double val_density, double *val_velocity, double val_energy, unsigned short val_ndim, 
		unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim;

	bool incompressible = config->GetIncompressible();

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate undivided laplacian (centred) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) Limiter = new double [nVar];

	/*--- Allocate truncation error for multigrid strategy ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;

	/*--- Solution and old solution initialization ---*/
	if (incompressible) {
		Solution[0] = config->GetPressure_FreeStreamND()/config->GetDensity_FreeStreamND();
		Solution_Old[0] = config->GetPressure_FreeStreamND()/config->GetDensity_FreeStreamND();
		for (iDim = 0; iDim < nDim; iDim++) {
			Solution[iDim+1] = val_velocity[iDim];
			Solution_Old[iDim+1] = val_velocity[iDim];
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

	/*--- Allocate and initializate solution for dual time strategy ---*/
	if (config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) {
		if (incompressible) {
			Solution_time_n[0] = config->GetPressure_FreeStreamND();
			Solution_time_n1[0] = config->GetPressure_FreeStreamND();
			for (iDim = 0; iDim < nDim; iDim++) {
				Solution_time_n[iDim+1] = val_velocity[iDim];
				Solution_time_n1[iDim+1] = val_velocity[iDim];
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

	/*--- Allocate and initializate the primitive gradient variable (nVar+1) variables ---*/
	Gradient_Primitive = new double* [nVar+1];
	for (iVar = 0; iVar < nVar+1; iVar++)
		Gradient_Primitive[iVar] = new double [nDim];
}

CEulerVariable::CEulerVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate undivided laplacian (centred) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) Limiter = new double [nVar];

	/*--- Allocate truncation error for multigrid strategy ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;

	/*--- Solution initialization ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}

	/*--- Allocate and initializate solution for dual time strategy ---*/
	if (config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) {
		Solution_time_n = new double [nVar];
		Solution_time_n1 = new double [nVar];

		for (iVar = 0; iVar < nVar; iVar++) {
			Solution_time_n[iVar] = val_solution[iVar];
			Solution_time_n1[iVar] = val_solution[iVar];
		}
	}

	/*--- Allocate and initializate the primitive gradient variable ---*/
	Gradient_Primitive = new double* [nVar+1];
	for (iVar = 0; iVar < nVar+1; iVar++)
		Gradient_Primitive[iVar] = new double [nDim];
}

CEulerVariable::~CEulerVariable(void) {
	unsigned short iVar;

	delete [] Res_Conv; delete [] Res_Visc;
	delete [] Residual_Sum; delete [] Residual_Old;
	if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
	if (Limiter != NULL) delete [] Limiter;
	if (Solution_time_n != NULL) delete [] Solution_time_n;
	if (Solution_time_n1 != NULL) delete [] Solution_time_n1;
	delete [] TruncationError;

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
		ProjVel += Solution[iDim+1]*val_vector[iDim];

	return ProjVel;
}

void CEulerVariable::SetSoundSpeedInc(double *val_normal) {
	double ProjVel, Area;
	unsigned short iDim;

	ProjVel = 0.0; Area = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Area = val_normal[iDim]*val_normal[iDim];
		ProjVel += Solution[iDim+1]*val_normal[iDim];
	}
	Area = sqrt(Area);

	SoundSpeedInc = sqrt(ProjVel*ProjVel + BetaInc2*Area*Area) / Area; 
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
	delete [] Residual_Sum;
	delete [] Residual_Old;
	if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
	if (Limiter != NULL) delete [] Limiter;
	if (Solution_time_n != NULL) delete [] Solution_time_n;
	if (Solution_time_n1 != NULL) delete [] Solution_time_n1;	
	delete [] TruncationError;

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

void CNSVariable::SetEddyViscosity(unsigned short val_Kind_Turb_Model, double *Turb_Solution) {
	// Parameters for the Spalart-Allmaras model
	double cv1_3 = 7.1*7.1*7.1;
	double nu, nu_hat, Ji, Ji_3, fv1;
	// Parameters for the Menter SST model
	double kine, omega;

	switch (val_Kind_Turb_Model) {
	case NONE :
		EddyViscosity = 0.0;
		break;
	case SA :
		nu = LaminarViscosity / Solution[0];
		nu_hat = Turb_Solution[0];
		Ji = nu_hat/nu;
		Ji_3 = Ji*Ji*Ji;
		fv1 = Ji_3/(Ji_3+cv1_3);
		EddyViscosity = Solution[0]*fv1*nu_hat;
		break;
	case SA_COMP :
		nu = LaminarViscosity / Solution[0];
		nu_hat = Turb_Solution[0] / Solution[0];
		Ji = nu_hat/nu;
		Ji_3 = Ji*Ji*Ji;
		fv1 = Ji_3/(Ji_3+cv1_3);
		EddyViscosity = Solution[0]*fv1*nu_hat;
		break;
	case SST :
		kine  = Turb_Solution[0] / Solution[0];
		omega = Turb_Solution[0] / Solution[0];
		EddyViscosity = Solution[0]*kine/omega;
		break;
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

CTurbVariable::CTurbVariable(void) : CVariable() {}

CTurbVariable::~CTurbVariable(void) {
	if (Limiter != NULL) delete [] Limiter;
}

CTurbVariable::CTurbVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {

	if (config->GetKind_SlopeLimit_Turb() != NONE) Limiter = new double [nVar];
}

CTurbSAVariable::CTurbSAVariable(void) : CTurbVariable() {}

CTurbSAVariable::~CTurbSAVariable(void) {
	if (Limiter != NULL) delete [] Limiter;
}

CTurbSAVariable::CTurbSAVariable(double val_nu_tilde, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_ndim, val_nvar,config) {

	// Initialization of variables
	Solution[0] = val_nu_tilde;		Solution_Old[0] = val_nu_tilde;

}

CTurbSSTVariable::CTurbSSTVariable(void) : CTurbVariable() {}

CTurbSSTVariable::~CTurbSSTVariable(void) {
	if (Limiter != NULL) delete [] Limiter;
}

CTurbSSTVariable::CTurbSSTVariable(double val_rho_kine, double val_rho_omega, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_ndim, val_nvar,config) {

	// Initialization of variables
	Solution[0] = val_rho_kine;     Solution_Old[0] = val_rho_kine;
	Solution[1] = val_rho_omega;	Solution_Old[1] = val_rho_omega;

}

void CTurbSSTVariable::SetF1blending(double val_viscosity, double val_dist, double val_density){
	unsigned short iDim;
	double sigma_w2 = 0.856;
	double cross_diff = 0.0, CD_kw, arg2, arg1, arg1_4;

	for (iDim = 0; iDim < nDim; iDim++)
		cross_diff += Gradient[0][iDim]*Gradient[1][iDim];
	cross_diff *= 2*val_density*sigma_w2/Solution[1];
	CD_kw = max(cross_diff,pow(10.0,-20.0));

	arg2 = max(sqrt(Solution[0])/(0.09*Solution[1]*val_dist),500*val_viscosity/(val_dist*val_dist*Solution[1]));
	arg1 = min(arg2,4*val_density*sigma_w2*Solution[0]/(CD_kw*val_dist*val_dist));
	arg1_4 = pow(arg1,4.0);
	F1 = tanh(arg1_4);
}

double CTurbSSTVariable::GetF1blending(){
	return F1;
}

CPlasmaMonatomicVariable::CPlasmaMonatomicVariable(void) : CVariable() { }

CPlasmaMonatomicVariable::CPlasmaMonatomicVariable(double *val_density, double **val_velocity, double *val_energy, unsigned short val_ndim, 
		unsigned short val_nvar, unsigned short val_nSpecies, unsigned short  val_nFluids, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, loc = 0, iFluids, iDim;

	nFluids  = val_nFluids;	
	nSpecies = val_nSpecies;

	/*--- Gamma determination ---*/
	double Gamma = config->GetGamma();
	double Gamma_Minus_One = Gamma - 1.0;


	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate undivided laplacian (centred) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) Limiter = new double [nVar];

	/*--- Allocate truncation error for multigrid strategy ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;

	/*--- Solution and old solution initialization ---*/	

	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		loc = (nDim+2)*iFluids;

		Solution[loc + 0] = val_density[iFluids];
		Solution[loc + 1] = val_density[iFluids]*val_velocity[iFluids][0];
		Solution[loc + 2] = val_density[iFluids]*val_velocity[iFluids][1];
		Solution[loc + 3] = val_density[iFluids]*val_energy[iFluids];

		Solution_Old[loc + 0] = val_density[iFluids];
		Solution_Old[loc + 1] = val_density[iFluids]*val_velocity[iFluids][0];
		Solution_Old[loc + 2] = val_density[iFluids]*val_velocity[iFluids][1];
		Solution_Old[loc + 3] = val_density[iFluids]*val_energy[iFluids];

		if (nDim == 3) {
			Solution[loc + 3] = val_density[iFluids]*val_velocity[iFluids][2];
			Solution[loc + 4] = val_density[iFluids]*val_energy[iFluids];

			Solution_Old[loc + 3] = val_density[iFluids]*val_velocity[iFluids][2];
			Solution_Old[loc + 4] = val_density[iFluids]*val_energy[iFluids];
		}
	}


	/*--- Allocate and initializate solution for dual time strategy ---*/
	if (config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) {
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

	/*--- Allocate and initializate the primitive gradient variable ---*/
	Gradient_Primitive = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Gradient_Primitive[iVar] = new double [nDim];

	Pressure	= new double [nFluids];
	SoundSpeed  = new double [nFluids];
	Velocity2   = new double [nFluids];
	Velocity    = new double*[nFluids];
	Enthalpy    = new double [nFluids];
	Source	    = new double [nFluids*(nDim+2)];
	Source_Jacobian = new double* [nFluids*(nDim+2)];
	Max_Lambda_Inv_MultiSpecies = new double [nFluids];

	for (iVar = 0; iVar < nVar; iVar ++) {
		Source[iVar] = 0.0;
		Source_Jacobian[iVar] = new double[nFluids*(nDim+2)];
	}

	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		loc      = (nDim+2)*iFluids;
		double r = Solution[loc + 0];
		double energy      = Solution[loc + nDim + 1];
		Velocity[iFluids]  = new double [nDim];
		Velocity2[iFluids] = 0.0;
		for (iDim = 0; iDim < nDim; iDim ++) {
			Velocity2[iFluids] 		    += pow((Solution[loc + iDim + 1]/r),2);
			Velocity[iFluids][iDim] 	 = Solution[loc + iDim + 1]/r;
		}
		Pressure[iFluids] 		= Gamma_Minus_One*r*(energy/r - 0.5*Velocity2[iFluids]);
		SoundSpeed[iFluids] 	= sqrt(fabs(Gamma*Gamma_Minus_One*(energy/r - 0.5*Velocity2[iFluids])));
		Enthalpy[iFluids]		= (energy + Pressure[iFluids])/r;
		Max_Lambda_Inv_MultiSpecies[iFluids] = 0.0;		
	}
}


CPlasmaMonatomicVariable::CPlasmaMonatomicVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate undivided laplacian (centred) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) Limiter = new double [nVar];

	/*--- Allocate truncation error for multigrid strategy ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;

	/*--- Solution initialization ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}

	/*--- Allocate and initializate solution for dual time strategy ---*/
	if (config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) {
		Solution_time_n = new double [nVar];
		Solution_time_n1 = new double [nVar];

		for (iVar = 0; iVar < nVar; iVar++) {
			Solution_time_n[iVar] = val_solution[iVar];
			Solution_time_n1[iVar] = val_solution[iVar];
		}
	}

	/*--- Allocate and initializate the primitive gradient variable ---*/
	Gradient_Primitive = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Gradient_Primitive[iVar] = new double [nDim];
}

CPlasmaMonatomicVariable::~CPlasmaMonatomicVariable(void) {
	unsigned short iVar;

	delete [] Res_Conv; delete [] Res_Visc;
	delete [] Residual_Sum; delete [] Residual_Old;
	if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
	if (Limiter != NULL) delete [] Limiter;
	if (Solution_time_n != NULL) delete [] Solution_time_n;
	if (Solution_time_n1 != NULL) delete [] Solution_time_n1;
	delete [] TruncationError;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;

	for (iVar = 0; iVar < nVar+1; iVar++)
		delete [] Gradient_Primitive[iVar];
	delete [] Gradient_Primitive;


	for (iVar = 0; iVar < nFluids; iVar++) {
		delete [] Velocity[iVar];
	}
	delete [] Velocity;
	delete [] Pressure;
	delete [] Velocity2;
	delete [] SoundSpeed;
	delete [] Enthalpy;
	delete [] Source;

	for (iVar = 0; iVar < nVar; iVar ++ )
		delete [] Source_Jacobian[iVar];
}

void CPlasmaMonatomicVariable::SetGradient_PrimitiveZero(void) {
	unsigned short iVar, iDim;
	for (iVar = 0; iVar < nVar; iVar++)
		for (iDim = 0; iDim < nDim; iDim++)
			Gradient_Primitive[iVar][iDim] = 0.0;
}


double CPlasmaMonatomicVariable::GetProjVel(double *val_vector, unsigned short iFluids) {
	double ProjVel;
	unsigned short iDim, loc = 0; 
	loc = (nDim+2)*iFluids;
	ProjVel = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVel += Solution[loc + iDim+1]*val_vector[iDim]/Solution[loc + 0];
	return ProjVel;
}

void CPlasmaMonatomicVariable::SetVelocity2() { 

	unsigned short loc = 0, iFluids, iDim;
	double r;
	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		loc = (nDim+2)*iFluids;
		r = Solution[loc + 0];
		Velocity2[iFluids] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Velocity2[iFluids] += pow( (Solution[loc + iDim +1]/r),2);
	}
}

void CPlasmaMonatomicVariable::SetSoundSpeed(double Gamma) {

	unsigned short loc = 0, iFluids, iDim;
	double r, velSqr, energy;
	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		loc = (nDim+2)*iFluids;
		r = Solution[loc + 0];
		velSqr = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			velSqr += pow( (Solution[loc + iDim +1]/r),2);
		energy = Solution[loc + nDim + 1];
		SoundSpeed[iFluids] = sqrt(fabs(Gamma*(Gamma-1.0)*(energy/r-0.5*velSqr)));
	}
}

void CPlasmaMonatomicVariable::SetPressure(double Gamma ) {

	unsigned short loc = 0, iFluids, iDim;
	double r, velSqr, energy;

	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		loc = (nDim+2)*iFluids;
		r = Solution[loc + 0];
		velSqr = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			velSqr += pow( (Solution[loc + iDim +1]/r),2);
		energy = Solution[loc + nDim + 1];
		Pressure[iFluids]  = (Gamma - 1.0)*r*(fabs(energy/r-0.5*velSqr));
	}
}

void CPlasmaMonatomicVariable::SetEnthalpy() { 

	unsigned short loc = 0, iFluids;
	double r, energy;

	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		loc = (nDim+2)*iFluids;
		r = Solution[loc + 0];
		energy = Solution[loc + nDim + 1];
		Enthalpy[iFluids]  = (energy + Pressure[iFluids])/r;
	}
}

double CPlasmaMonatomicVariable::GetVelocity(unsigned short val_dim,unsigned short iFluids) {

	unsigned short loc = 0; 
	loc = (nDim+2)*iFluids;
	double r = Solution[loc + 0];
	double r_vel = Solution[loc + val_dim + 1];
	return (r_vel/r);  
}

void CPlasmaMonatomicVariable::SetSource(double * res) { 

	for (int iVar = 0; iVar < nVar; iVar ++)
		Source[iVar] = res[iVar];

}

void CPlasmaMonatomicVariable::SetSourceJacobian(double ** src_jacobian) { 

	for (int iVar = 0; iVar < nVar; iVar ++) {
		for (int jVar = 0; jVar < nVar; jVar ++)
			Source_Jacobian[iVar][jVar] = src_jacobian[iVar][jVar];
	}
}

CPlasmaDiatomicVariable::CPlasmaDiatomicVariable(void) : CVariable() { }

CPlasmaDiatomicVariable::CPlasmaDiatomicVariable(double *val_density, double **val_velocity, double *val_energy, double *val_energy_vib, 
		double *val_energy_formation, double *val_enthalpy_formation, unsigned short val_ndim,
		unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nMonatomics,
		unsigned short val_nDiatomics, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, loc = 0, iSpecies;

	nMonatomics  = val_nMonatomics;	
	nDiatomics = val_nDiatomics;
	nSpecies = val_nSpecies;

	/*--- Gamma determination ---*/
	double GammaMonatomics = 5.0/3.0;
	double GammaDiatomics = 7.0/5.0;

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate limiter (upwind)---*/
	Limiter = new double [nVar];

	/*--- Allocate truncation error for multigrid strategy ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;

	/*--- Solution and old solution initialization ---*/		
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if (nDim == 2) {

			if ( iSpecies < nDiatomics ) loc = 5*iSpecies;
			else loc = 5*nDiatomics + 4*(iSpecies-nDiatomics);

			Solution[loc + 0] = val_density[iSpecies]; 
			Solution[loc + 1] = val_density[iSpecies]*val_velocity[iSpecies][0]; 
			Solution[loc + 2] = val_density[iSpecies]*val_velocity[iSpecies][1]; 
			Solution[loc + 3] = val_density[iSpecies]*val_energy[iSpecies]; 
			if ( iSpecies < nDiatomics ) Solution[loc + 4] = val_density[iSpecies]*val_energy_vib[iSpecies]; 
			Solution_Old[loc + 0] = val_density[iSpecies]; 
			Solution_Old[loc + 1] = val_density[iSpecies]*val_velocity[iSpecies][0]; 
			Solution_Old[loc + 2] = val_density[iSpecies]*val_velocity[iSpecies][1]; 
			Solution_Old[loc + 3] = val_density[iSpecies]*val_energy[iSpecies];
			if ( iSpecies < nDiatomics ) Solution_Old[loc + 4] = val_density[iSpecies]*val_energy_vib[iSpecies]; 
		}
	}

	/*--- Allocate and initializate the primitive gradient variable ---*/
	Gradient_Primitive = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Gradient_Primitive[iVar] = new double [nDim];

	Pressure   	= new double [nSpecies];
	SoundSpeed  = new double [nSpecies];
	Velocity2   = new double [nSpecies];
	Velocity    = new double*[nSpecies];
	Enthalpy    = new double [nSpecies];

	Energy_Formation   = new double [nSpecies];
	Enthalpy_Formation = new double [nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Energy_Formation[iSpecies]   = val_energy_formation[iSpecies];
		Enthalpy_Formation[iSpecies] = val_enthalpy_formation[iSpecies];
	}

	Max_Lambda_Inv_MultiSpecies = new double [nSpecies];

	if (nDim == 2 ) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

			if ( iSpecies < nDiatomics ) loc = 5*iSpecies;
			else loc = 5*nDiatomics + 4*(iSpecies-nDiatomics);	

			double rho = Solution[loc + 0];
			double u = Solution[loc + 1] / rho;
			double v = Solution[loc + 2] / rho;
			double energy = Solution[loc + 3];
			double energy_vib = 0.0;
			if ( iSpecies < nDiatomics ) energy_vib = Solution[loc + 4];
			double energy_el = 0.0;

			Velocity2[iSpecies] 		= u*u + v*v;
			Velocity[iSpecies] = new double [nDim];
			Velocity[iSpecies][0] 	= u; Velocity[iSpecies][1] 	= v;
			if ( iSpecies < nDiatomics ) {
				Pressure[iSpecies] 		= (GammaDiatomics - 1.0)*rho*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Energy_Formation[iSpecies]);
				SoundSpeed[iSpecies] 	= sqrt(fabs(GammaDiatomics*(GammaDiatomics - 1.0)*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Energy_Formation[iSpecies])));
			}
			else {
				Pressure[iSpecies] 		= (GammaMonatomics - 1.0)*rho*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Energy_Formation[iSpecies]);
				SoundSpeed[iSpecies] 	= sqrt(fabs(GammaMonatomics*(GammaMonatomics - 1.0)*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Energy_Formation[iSpecies])));
			}
			Enthalpy[iSpecies]		= (energy + Pressure[iSpecies] + Enthalpy_Formation[iSpecies] - Energy_Formation[iSpecies])/rho;
			Max_Lambda_Inv_MultiSpecies[iSpecies] = 0.0;		
		}
	}
}

CPlasmaDiatomicVariable::CPlasmaDiatomicVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, 
		CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate undivided laplacian (centred) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) Limiter = new double [nVar];

	/*--- Allocate truncation error for multigrid strategy ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;

	/*--- Solution initialization ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}

	/*--- Allocate and initializate solution for dual time strategy ---*/
	if (config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) {
		Solution_time_n = new double [nVar];
		Solution_time_n1 = new double [nVar];

		for (iVar = 0; iVar < nVar; iVar++) {
			Solution_time_n[iVar] = val_solution[iVar];
			Solution_time_n1[iVar] = val_solution[iVar];
		}
	}

	/*--- Allocate and initializate the primitive gradient variable ---*/
	Gradient_Primitive = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Gradient_Primitive[iVar] = new double [nDim];
}

CPlasmaDiatomicVariable::~CPlasmaDiatomicVariable(void) {
	unsigned short iVar;

	delete [] Res_Conv; delete [] Res_Visc;
	delete [] Residual_Sum; delete [] Residual_Old;
	if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
	if (Limiter != NULL) delete [] Limiter;
	if (Solution_time_n != NULL) delete [] Solution_time_n;
	if (Solution_time_n1 != NULL) delete [] Solution_time_n1;
	delete [] TruncationError;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;

	for (iVar = 0; iVar < nVar+1; iVar++)
		delete [] Gradient_Primitive[iVar];
	delete [] Gradient_Primitive;


	for (iVar = 0; iVar < nFluids; iVar++) {
		delete [] Velocity[iVar];
	}
	delete [] Velocity;
	delete [] Pressure;
	delete [] Velocity2;
	delete [] SoundSpeed;
	delete [] Enthalpy;
	delete [] Source;
	delete [] Energy_Formation;
	delete [] Enthalpy_Formation;

	for (iVar = 0; iVar < nFluids*(nDim+2); iVar ++ )
		delete [] Source_Jacobian[iVar];
}

void CPlasmaDiatomicVariable::SetGradient_PrimitiveZero(void) {
	unsigned short iVar, iDim;
	for (iVar = 0; iVar < nVar; iVar++)
		for (iDim = 0; iDim < nDim; iDim++)
			Gradient_Primitive[iVar][iDim] = 0.0;
}


double CPlasmaDiatomicVariable::GetProjVel(double *val_vector, unsigned short val_Species, unsigned short val_nDiatomics) {
	double ProjVel;
	unsigned short iDim, loc = 0; 

	if (nDim == 2) {
		if (val_Species < val_nDiatomics) loc = 5*val_Species;
		else loc = 5*val_nDiatomics + 4*(val_Species-val_nDiatomics);
	}
	if (nDim == 3) {		
		if (val_Species < val_nDiatomics) loc = 6*val_Species;
		else loc = 6*val_nDiatomics + 5*(val_Species-val_nDiatomics);
	}

	ProjVel = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVel += Solution[loc + iDim+1]*val_vector[iDim]/Solution[loc + 0];
	return ProjVel;
}

void CPlasmaDiatomicVariable::SetVelocity2() { 

	unsigned short loc = 0, iSpecies;
	double rho, u, v, w;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if (nDim == 2) {
			if ( iSpecies < nDiatomics ) loc = 5*iSpecies;
			else loc = 5*nDiatomics + 4*(iSpecies-nDiatomics);	

			rho = Solution[loc + 0];
			u = Solution[loc + 1]/rho;
			v = Solution[loc + 2]/rho;
			Velocity2[iSpecies] = u*u + v*v;

		}
		if (nDim == 3) {
			if ( iSpecies < nDiatomics ) loc = 6*iSpecies;
			else loc = 6*nDiatomics + 5*(iSpecies-nDiatomics);	

			rho = Solution[loc + 0];
			u = Solution[loc + 1]/rho;
			v = Solution[loc + 2]/rho;
			w = Solution[loc + 3]/rho;
			Velocity2[iSpecies] = u*u + v*v + w*w;			
		}		
	}
}

void CPlasmaDiatomicVariable::SetSoundSpeed(double GammaMonatomic, double GammaDiatomic) {

	unsigned short loc = 0, iSpecies;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = 5*iSpecies;
		else loc = 5*nDiatomics + 4*(iSpecies-nDiatomics);	

		double rho = Solution[loc + 0];
		double energy_el = 0.0;

		double energy = Solution[loc + nDim + 1];
		double energy_vib;
		if (iSpecies < nDiatomics )
			energy_vib = Solution[loc + nDim + 2];
		else 
			energy_vib = 0.0;


		if ( iSpecies < nDiatomics )
			SoundSpeed[iSpecies] 	= sqrt(fabs(GammaDiatomic*(GammaDiatomic - 1.0)*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Energy_Formation[iSpecies])));
//			SoundSpeed[iSpecies] = sqrt( GammaDiatomic * Pressure / rho );
		else {
			SoundSpeed[iSpecies] 	= sqrt(fabs(GammaMonatomic*(GammaMonatomic - 1.0)*(energy/rho - 0.5*Velocity2[iSpecies] - energy_el - Energy_Formation[iSpecies])));
//			SoundSpeed[iSpecies] = sqrt( GammaMonatomic * Pressure / rho );
		}
	}
}

void CPlasmaDiatomicVariable::SetPressure(double GammaMonatomic, double GammaDiatomic) {

	unsigned short loc = 0, iSpecies;	

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {

		if ( iSpecies < nDiatomics ) loc = 5*iSpecies;
		else loc = 5*nDiatomics + 4*(iSpecies-nDiatomics);	

		double rho = Solution[loc + 0];
		double energy = Solution[loc + nDim + 1];
		double energy_vib;
		if ( iSpecies < nDiatomics )
			energy_vib = Solution[loc + nDim + 2];
		double energy_el = 0.0;

		if ( iSpecies < nDiatomics )
			Pressure[iSpecies]  = (GammaDiatomic - 1.0)*rho*(fabs(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Energy_Formation[iSpecies]));
		else
			Pressure[iSpecies]  = (GammaMonatomic - 1.0)*rho*(fabs(energy/rho - 0.5*Velocity2[iSpecies] - energy_el - Energy_Formation[iSpecies]));
		
	}
}

void CPlasmaDiatomicVariable::SetEnthalpy() { 

	unsigned short loc = 0, iSpecies;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {

		if ( iSpecies < nDiatomics ) loc = 5*iSpecies;
		else loc = 5*nDiatomics + 4*(iSpecies-nDiatomics);

		double rho = Solution[loc + 0];
		double energy = Solution[loc + nDim + 1];
		Enthalpy[iSpecies]  = (energy + Pressure[iSpecies])/rho;
	}
}

double CPlasmaDiatomicVariable::GetVelocity(unsigned short val_dim,unsigned short val_Species) { 

	unsigned short loc = 0; 

	if ( val_Species < nDiatomics ) loc = (nDim+3)*val_Species;
	else loc = (nDim+3)*nDiatomics + (nDim+2)*(val_Species-nDiatomics);	
	
	double rho = Solution[loc + 0];
	double rho_vel = Solution[loc + val_dim + 1];
	return (rho_vel/rho);  
}

void CPlasmaDiatomicVariable::SetSource(double * res) { 

	for (int iVar = 0; iVar < nFluids*(nDim+2); iVar ++)
		Source[iVar] = res[iVar];

}

void CPlasmaDiatomicVariable::SetSourceJacobian(double ** src_jacobian) { 

	for (int iVar = 0; iVar < nFluids*(nDim+2); iVar ++) {
		for (int jVar = 0; jVar < nFluids*(nDim+2); jVar ++)
			Source_Jacobian[iVar][jVar] = src_jacobian[iVar][jVar];
	}

}

CLevelSetVariable::CLevelSetVariable(void) : CVariable() {}

CLevelSetVariable::CLevelSetVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate limiter (upwind)---*/
	if (config->GetKind_SlopeLimit_LevelSet() != NONE) Limiter = new double [nVar];

}

CLevelSetVariable::CLevelSetVariable(double val_levelset, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar,config) {
	unsigned short iVar;

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate limiter (upwind)---*/
	if (config->GetKind_SlopeLimit_LevelSet() != NONE) Limiter = new double [nVar];

	/*--- Solution and old solution initialization ---*/
	Solution[0] = val_levelset;		Solution_Old[0] = val_levelset;


}

CLevelSetVariable::~CLevelSetVariable(void) {
	unsigned short iVar;

	delete [] Res_Conv; delete [] Res_Visc;
	delete [] Residual_Sum; delete [] Residual_Old;
	if (Limiter != NULL) delete [] Limiter;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;
}

CCombustionVariable::CCombustionVariable(void) : CVariable() {}

CCombustionVariable::CCombustionVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) { 
	unsigned short iVar;

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

}

CCombustionVariable::CCombustionVariable(double val_lambda, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar,config) {
	unsigned short iVar;	

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];


	/*--- Solution and old solution initialization ---*/
	Solution[0] = val_lambda;		Solution_Old[0] = val_lambda;

}

CCombustionVariable::~CCombustionVariable(void) {
}
