/*!
 * \file variable_direct.cpp
 * \brief Definition of the solution fields.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.1.
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

	incompressible = config->GetIncompressible();
	bool freesurface = (config->GetKind_Solver() == FREE_SURF_EULER);

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate auxiliar vector for free surface source term ---*/
	if (freesurface) Grad_AuxVar = new double [nDim];

	/*--- Allocate undivided laplacian (centred) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) Limiter = new double [nVar];

	/*--- Allocate truncation error for multigrid strategy ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;

	/*--- Solution and old solution initialization ---*/
	if (incompressible) {
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

	/*--- Allocate and initializate solution for dual time strategy ---*/
	if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
		if (incompressible) {
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

	/*--- Allocate and initializate the primitive gradient variable (nVar+1) variables ---*/
	Gradient_Primitive = new double* [nVar+1];
	for (iVar = 0; iVar < nVar+1; iVar++)
		Gradient_Primitive[iVar] = new double [nDim];
}

CEulerVariable::CEulerVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

	incompressible = config->GetIncompressible();
	bool freesurface = (config->GetKind_Solver() == FREE_SURF_EULER);

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate auxiliar vector for free surface source term ---*/
	if (freesurface) Grad_AuxVar = new double [nDim];

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
	if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || 
			(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
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

	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour;
	delete [] Residual_Sum; delete [] Residual_Old;
	if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
	if (Limiter != NULL) delete [] Limiter;
	if (Solution_time_n != NULL) delete [] Solution_time_n;
	if (Solution_time_n1 != NULL) delete [] Solution_time_n1;
	delete [] TruncationError;

	if (Grad_AuxVar != NULL) delete [] Grad_AuxVar;

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
		ProjVel += (Solution[iDim+1]/DensityInc)*val_vector[iDim];

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

	incompressible = config->GetIncompressible();

	Temperature_Ref = config->GetTemperature_Ref();
	Viscosity_Ref   = config->GetViscosity_Ref();
	Viscosity_Inf   = config->GetViscosity_FreeStreamND();
	Prandtl_Lam     = config->GetPrandtl_Lam();
	Prandtl_Turb    = config->GetPrandtl_Turb();

}

CNSVariable::CNSVariable(double *val_solution, unsigned short val_ndim, 
		unsigned short val_nvar, CConfig *config) : CEulerVariable(val_solution, val_ndim, val_nvar, config) {

	incompressible = config->GetIncompressible();

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

void CNSVariable::SetEddyViscosity(unsigned short val_Kind_Turb_Model, CVariable *Turb_Solution) {

	switch (val_Kind_Turb_Model) {
	case NONE : EddyViscosity = 0.0; break;
	case SA : case SST: EddyViscosity = Turb_Solution->GetmuT(); break;
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
}

CTurbSAVariable::CTurbSAVariable(double val_nu_tilde, double val_muT, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_ndim, val_nvar,config) {

	// Initialization of variables
	Solution[0] = val_nu_tilde;		Solution_Old[0] = val_nu_tilde;

	// Initialization of eddy viscosity
	muT = val_muT;
}

CTurbSSTVariable::CTurbSSTVariable(void) : CTurbVariable() {}

CTurbSSTVariable::~CTurbSSTVariable(void) {
	if (Limiter != NULL) delete [] Limiter;
}

CTurbSSTVariable::CTurbSSTVariable(double val_kine, double val_omega, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_ndim, val_nvar,config) {

	// Initialization of variables
	Solution[0] = val_kine;     Solution_Old[0] = val_kine;
	Solution[1] = val_omega;	Solution_Old[1] = val_omega;

	F1   = 1.0;
	F2   = 0.0;
	CDkw = 0.0;

}

void CTurbSSTVariable::SetBlendingFunc(double val_viscosity, double val_dist, double val_density){
	unsigned short iDim;
	double sigma_om2 = 0.856, beta_star = 0.09;
	double arg2, arg2A, arg2B, arg1;

	// Cross diffusion
	CDkw = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		CDkw += Gradient[0][iDim]*Gradient[1][iDim];
	CDkw *= 2*val_density*sigma_om2/Solution[1];
	CDkw = max(CDkw,pow(10.0,-20.0));

	// F1
	arg2A = sqrt(Solution[0])/(beta_star*Solution[1]*val_dist);
	arg2B = 500*val_viscosity/(val_density*val_dist*val_dist*Solution[1]);
	arg2 = max(arg2A,arg2B);
	arg1 = min(arg2,4*val_density*sigma_om2*Solution[0]/(CDkw*val_dist*val_dist));
	F1 = tanh(pow(arg1,4.0));

	// F2
	arg2 = max(2*arg2A,arg2B);
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

CPlasmaVariable::CPlasmaVariable(double *val_density, double **val_velocity, double *val_energy, unsigned short val_ndim,
		unsigned short val_nvar, unsigned short val_nSpecies, unsigned short  val_nFluids, unsigned short val_nMonatomics,
		unsigned short val_nDiatomics, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, loc = 0, iSpecies, iDim;

	nFluids  = val_nFluids;
	nSpecies = val_nSpecies;
	nMonatomics  = val_nMonatomics;
		nDiatomics = val_nDiatomics;

	/*--- Gamma determination ---*/
	double Gamma = config->GetGamma();
	double Gamma_Minus_One = Gamma - 1.0;
	Prandtl_Lam     = config->GetPrandtl_Lam();


	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	Source = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate undivided laplacian (centred) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTRED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_UPWIND) Limiter = new double [nVar];

	/*--- Allocate truncation error for multigrid strategy ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;

	/*--- Solution and old solution initialization ---*/

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		loc = (nDim+2)*iSpecies;

		Solution[loc + 0] = val_density[iSpecies];
		Solution[loc + 1] = val_density[iSpecies]*val_velocity[iSpecies][0];
		Solution[loc + 2] = val_density[iSpecies]*val_velocity[iSpecies][1];
		Solution[loc + 3] = val_density[iSpecies]*val_energy[iSpecies];

		Solution_Old[loc + 0] = val_density[iSpecies];
		Solution_Old[loc + 1] = val_density[iSpecies]*val_velocity[iSpecies][0];
		Solution_Old[loc + 2] = val_density[iSpecies]*val_velocity[iSpecies][1];
		Solution_Old[loc + 3] = val_density[iSpecies]*val_energy[iSpecies];

		if (nDim == 3) {
			Solution[loc + 3] = val_density[iSpecies]*val_velocity[iSpecies][2];
			Solution[loc + 4] = val_density[iSpecies]*val_energy[iSpecies];

			Solution_Old[loc + 3] = val_density[iSpecies]*val_velocity[iSpecies][2];
			Solution_Old[loc + 4] = val_density[iSpecies]*val_energy[iSpecies];
		}
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

	/*--- Allocate and initializate the primitive gradient variable ---*/
	Gradient_Primitive = new double* [nVar+nSpecies];
	for (iVar = 0; iVar < nVar+nSpecies; iVar++)
		Gradient_Primitive[iVar] = new double [nDim];

	Pressure	= new double [nSpecies];
	Temperature	= new double [nSpecies];
	SoundSpeed  = new double [nSpecies];
	Velocity2   = new double [nSpecies];
	Velocity    = new double*[nSpecies];
	Enthalpy    = new double [nSpecies];
	Source	    = new double [nSpecies*(nDim+2)];
	Source_Jacobian = new double* [nSpecies*(nDim+2)];
	Max_Lambda_Inv_MultiSpecies = new double [nSpecies];
	Max_Lambda_Visc_MultiSpecies = new double [nSpecies];
	Species_Delta_Time = new double [nSpecies];

	LaminarViscosity_MultiSpecies = new double [nSpecies];
	EddyViscosity_MultiSpecies = new double [nSpecies];
	Kappa = new double [nSpecies];
	for (iVar = 0; iVar < nVar; iVar ++) {
		Source[iVar] = 0.0;
		Source_Jacobian[iVar] = new double[nSpecies*(nDim+2)];
	}

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		loc      = (nDim+2)*iSpecies;
		double r = Solution[loc + 0];
		double energy      = Solution[loc + nDim + 1];

		Velocity[iSpecies]  = new double [nDim];
		Velocity2[iSpecies] = 0.0;
		for (iDim = 0; iDim < nDim; iDim ++) {
			Velocity2[iSpecies] 		    += pow((Solution[loc + iDim + 1]/r),2);
			Velocity[iSpecies][iDim] 	 = Solution[loc + iDim + 1]/r;
		}
		Pressure[iSpecies] 		= Gamma_Minus_One*r*(energy/r - 0.5*Velocity2[iSpecies]);
		SoundSpeed[iSpecies] 	= sqrt(fabs(Gamma*Gamma_Minus_One*(energy/r - 0.5*Velocity2[iSpecies])));
		Enthalpy[iSpecies]		= (energy + Pressure[iSpecies])/r;
		Max_Lambda_Inv_MultiSpecies[iSpecies] = 0.0;
		Max_Lambda_Visc_MultiSpecies[iSpecies] = 0.0;
		LaminarViscosity_MultiSpecies[iSpecies] = 0.0;
		EddyViscosity_MultiSpecies[iSpecies] = 0.0;
		Species_Delta_Time[iSpecies] = 1.0;
	}
}

CPlasmaVariable::CPlasmaVariable(double *val_density, double **val_velocity, double *val_energy, double *val_energy_vib,
		double *val_enthalpy_formation, unsigned short val_ndim,
		unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nMonatomics,
		unsigned short val_nDiatomics, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, loc = 0, iSpecies, iDim;

	nMonatomics  = val_nMonatomics;
	nDiatomics = val_nDiatomics;
	nSpecies = val_nSpecies;
	Prandtl_Lam     = config->GetPrandtl_Lam();

	/*--- Gamma determination ---*/
	double GammaMonatomic = config->GetGammaMonatomic();
	double GammaDiatomic = config->GetGammaDiatomic();

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	Source = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Res_Visc_RK[iVar] = new double [nVar];
		Source[iVar] = 0.0;
	}
	
	/*--- Allocate limiter (upwind)---*/
	Limiter = new double [nVar];

	/*--- Allocate truncation error for multigrid strategy ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;

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

	}

	/*--- Allocate and initializate the primitive gradient variable ---*/
	Gradient_Primitive = new double* [nVar+nSpecies];
	for (iVar = 0; iVar < nVar+nSpecies; iVar++)
		Gradient_Primitive[iVar] = new double [nDim];

	Pressure   	= new double [nSpecies];
	Temperature_tr = new double [nSpecies];
	SoundSpeed  = new double [nSpecies];
	Velocity2   = new double [nSpecies];
	Enthalpy    = new double [nSpecies];
	Enthalpy_Formation = new double [nSpecies];

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Enthalpy_Formation[iSpecies] = val_enthalpy_formation[iSpecies];
	}

	Max_Lambda_Inv_MultiSpecies = new double [nSpecies];
	Max_Lambda_Visc_MultiSpecies = new double [nSpecies];

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

		double rho = Solution[loc + 0];
		double energy = Solution[loc + nDim+1];
		double energy_vib = 0.0;
		if ( iSpecies < nDiatomics ) energy_vib = Solution[loc + nDim+2];
		double energy_el = 0.0;
		
		Velocity2[iSpecies] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity2[iSpecies] += val_velocity[iSpecies][iDim]*val_velocity[iSpecies][iDim];
		}
		
		if ( iSpecies < nDiatomics ) {
		  Pressure[iSpecies] 		= (GammaDiatomic - 1.0)*rho*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Enthalpy_Formation[iSpecies]);
		  SoundSpeed[iSpecies] 	= sqrt(fabs(GammaDiatomic*(GammaDiatomic - 1.0)*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Enthalpy_Formation[iSpecies])));
		}
		else {
		  Pressure[iSpecies] 		= (GammaMonatomic - 1.0)*rho*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Enthalpy_Formation[iSpecies]);
		  SoundSpeed[iSpecies] 	= sqrt(fabs(GammaMonatomic*(GammaMonatomic - 1.0)*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Enthalpy_Formation[iSpecies])));
		}
		Enthalpy[iSpecies]		= (energy - 0.5*rho*Velocity2[iSpecies] - energy_vib - rho*Enthalpy_Formation[iSpecies] + Pressure[iSpecies])/rho  + Enthalpy_Formation[iSpecies];
		Max_Lambda_Inv_MultiSpecies[iSpecies] = 0.0;
		Max_Lambda_Visc_MultiSpecies[iSpecies] = 0.0;
	}

//	Source	    = new double [nSpecies*(nDim+2)];
//	Source_Jacobian = new double* [nSpecies*(nDim+2)];	
	Source	    = new double [nVar];
	Source_Jacobian = new double* [nVar];

	for (iVar = 0; iVar < nVar; iVar ++) {
		Source[iVar] = 0.0;
//		Source_Jacobian[iVar] = new double[nSpecies*(nDim+2)];
		Source_Jacobian[iVar] = new double[nVar];

	}

}

CPlasmaVariable::CPlasmaVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nFluids, unsigned short val_nMonatomics,
		unsigned short val_nDiatomics,	CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, loc = 0, iSpecies, iDim;

	nFluids  = val_nFluids;
	nSpecies = val_nSpecies;
	nMonatomics  = val_nMonatomics;
	nDiatomics = val_nDiatomics;

	/*--- Gamma determination ---*/
	double Gamma = config->GetGamma();
	double Gamma_Minus_One = Gamma - 1.0;
	Prandtl_Lam     = config->GetPrandtl_Lam();


	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];

	/*--- Allocate undivided laplacian (centred) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTRED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_UPWIND) Limiter = new double [nVar];

	/*--- Allocate truncation error for multigrid strategy ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;


	/*--- Solution initialization ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}

	/*--- Allocate and initializate the primitive gradient variable ---*/
	Gradient_Primitive = new double* [nVar+nSpecies];
	for (iVar = 0; iVar < nVar+nSpecies; iVar++)
		Gradient_Primitive[iVar] = new double [nDim];

	Pressure	= new double [nSpecies];
	Temperature	= new double [nSpecies];
	SoundSpeed  = new double [nSpecies];
	Velocity2   = new double [nSpecies];
	Velocity    = new double*[nSpecies];
	Enthalpy    = new double [nSpecies];
	Source	    = new double [nSpecies*(nDim+2)];
	Source_Jacobian = new double* [nSpecies*(nDim+2)];
	Max_Lambda_Inv_MultiSpecies = new double [nSpecies];
	Max_Lambda_Visc_MultiSpecies = new double [nSpecies];
	Species_Delta_Time = new double [nSpecies];

	LaminarViscosity_MultiSpecies = new double [nSpecies];
	EddyViscosity_MultiSpecies = new double [nSpecies];
	Kappa = new double [nSpecies];
	for (iVar = 0; iVar < nVar; iVar ++) {
		Source[iVar] = 0.0;
		Source_Jacobian[iVar] = new double[nSpecies*(nDim+2)];
	}


	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		loc      = (nDim+2)*iSpecies;
		double r = Solution[loc + 0];
		double energy      = Solution[loc + nDim + 1];

		Velocity[iSpecies]  = new double [nDim];
		Velocity2[iSpecies] = 0.0;
		for (iDim = 0; iDim < nDim; iDim ++) {
			Velocity2[iSpecies] 		    += pow((Solution[loc + iDim + 1]/r),2);
			Velocity[iSpecies][iDim] 	 = Solution[loc + iDim + 1]/r;
		}
		Pressure[iSpecies] 		= Gamma_Minus_One*r*(energy/r - 0.5*Velocity2[iSpecies]);
		SoundSpeed[iSpecies] 	= sqrt(fabs(Gamma*Gamma_Minus_One*(energy/r - 0.5*Velocity2[iSpecies])));
		Enthalpy[iSpecies]		= (energy + Pressure[iSpecies])/r;
		Max_Lambda_Inv_MultiSpecies[iSpecies] = 0.0;
		Max_Lambda_Visc_MultiSpecies[iSpecies] = 0.0;
		LaminarViscosity_MultiSpecies[iSpecies] = 0.0;
		EddyViscosity_MultiSpecies[iSpecies] = 0.0;
		Species_Delta_Time[iSpecies] = 1.0;
	}
}

CPlasmaVariable::CPlasmaVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, 
																 unsigned short val_nSpecies, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, loc = 0, iSpecies, iDim;
	
	nSpecies = val_nSpecies;
	nDim = val_ndim;
	nVar = val_nvar;
	
	/*--- Gamma determination ---*/
	double GammaDiatomic = config->GetGammaDiatomic();
	double GammaMonatomic = config->GetGammaMonatomic();
	Prandtl_Lam     = config->GetPrandtl_Lam();
	
	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar]; Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	Res_Visc_RK = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Res_Visc_RK[iVar] = new double [nVar];
	
	/*--- Allocate undivided laplacian (centred) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTRED) Undivided_Laplacian = new double [nVar];
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_UPWIND) Limiter = new double [nVar];
	
	/*--- Allocate truncation error for multigrid strategy ---*/
	TruncationError = new double [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		TruncationError[iVar] = 0.0;
	
	
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
	Temperature_tr = new double [nSpecies];
	SoundSpeed  = new double [nSpecies];
	Velocity2   = new double [nSpecies];
	Enthalpy    = new double [nSpecies];
	Enthalpy_Formation = new double [nSpecies];
	
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Enthalpy_Formation[iSpecies] = 	config->GetEnthalpy_Formation(iSpecies);
	}
	
	Max_Lambda_Inv_MultiSpecies = new double [nSpecies];
	Max_Lambda_Visc_MultiSpecies = new double [nSpecies];
	
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		
		double rho = Solution[loc + 0];
		double energy = Solution[loc + nDim+1];
		double energy_vib = 0.0;
		if ( iSpecies < nDiatomics ) energy_vib = Solution[loc + nDim+2];
		double energy_el = 0.0;
		
		Velocity2[iSpecies] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity2[iSpecies] += Solution[loc+1+iDim]*Solution[loc+1+iDim];
		}
		
		if ( iSpecies < nDiatomics ) {
		  Pressure[iSpecies] 		= (GammaDiatomic - 1.0)*rho*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Enthalpy_Formation[iSpecies]);
		  SoundSpeed[iSpecies] 	= sqrt(fabs(GammaDiatomic*(GammaDiatomic - 1.0)*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Enthalpy_Formation[iSpecies])));
		}
		else {
		  Pressure[iSpecies] 		= (GammaMonatomic - 1.0)*rho*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Enthalpy_Formation[iSpecies]);
		  SoundSpeed[iSpecies] 	= sqrt(fabs(GammaMonatomic*(GammaMonatomic - 1.0)*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - Enthalpy_Formation[iSpecies])));
		}
		Enthalpy[iSpecies]		= (energy - 0.5*rho*Velocity2[iSpecies] - energy_vib - rho*Enthalpy_Formation[iSpecies] + Pressure[iSpecies])/rho  + Enthalpy_Formation[iSpecies];
		Max_Lambda_Inv_MultiSpecies[iSpecies] = 0.0;
		Max_Lambda_Visc_MultiSpecies[iSpecies] = 0.0;
	}
	
	Source	    = new double[nVar];
	Source_Jacobian = new double*[nVar];
	
	for (iVar = 0; iVar < nVar; iVar ++) {
		Source[iVar] = 0.0;
		Source_Jacobian[iVar] = new double[nVar];
	}
}


CPlasmaVariable::~CPlasmaVariable(void) {
	unsigned short iVar;

	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour;
	delete [] Residual_Sum; delete [] Residual_Old;
	if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
	if (Limiter != NULL) delete [] Limiter;
	if (Solution_time_n != NULL) delete [] Solution_time_n;
	if (Solution_time_n1 != NULL) delete [] Solution_time_n1;
	delete [] TruncationError;
	delete [] Max_Lambda_Inv_MultiSpecies;
	delete [] Max_Lambda_Visc_MultiSpecies;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;

	for (iVar = 0; iVar < nVar+nSpecies; iVar++)
		delete [] Gradient_Primitive[iVar];
	delete [] Gradient_Primitive;


	for (iVar = 0; iVar < nFluids; iVar++) {
		delete [] Velocity[iVar];
	}
	delete [] Velocity;
	delete [] Pressure;
	delete [] Velocity2;
	delete [] Temperature_tr;
	delete [] SoundSpeed;
	delete [] Enthalpy;
	delete [] Source;
	delete [] Enthalpy_Formation;

	delete [] Source;
	for (iVar = 0; iVar < nVar; iVar ++ )
		delete [] Source_Jacobian[iVar];
	delete [] Source_Jacobian;
}

void CPlasmaVariable::SetGradient_PrimitiveZero(void) {
	unsigned short iVar, iDim;
	for (iVar = 0; iVar < nVar+nSpecies; iVar++)
		for (iDim = 0; iDim < nDim; iDim++)
			Gradient_Primitive[iVar][iDim] = 0.0;
}

double CPlasmaVariable::GetProjVel(double *val_vector, unsigned short iFluids) {
	double ProjVel;
	unsigned short iDim, loc = 0;
	loc = (nDim+2)*iFluids;
	ProjVel = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVel += Solution[loc + iDim+1]*val_vector[iDim]/Solution[loc + 0];
	return ProjVel;
}

double CPlasmaVariable::GetProjVel(double *val_vector, unsigned short val_Species, unsigned short val_nDiatomics) {
	double ProjVel;
	unsigned short iDim, loc = 0;
		
	if ( val_Species < val_nDiatomics ) loc = (nDim+3)*val_Species;
	else loc = (nDim+3)*val_nDiatomics + (nDim+2)*(val_Species-val_nDiatomics);

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

void CPlasmaVariable::SetSoundSpeed(double Gamma) {

	unsigned short loc = 0, iSpecies, iDim;
	double r, velSqr, energy;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		loc = (nDim+2)*iSpecies;
		r = Solution[loc + 0];
		velSqr = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			velSqr += pow( (Solution[loc + iDim +1]/r),2);
		energy = Solution[loc + nDim + 1];
		SoundSpeed[iSpecies] = sqrt(fabs(Gamma*(Gamma-1.0)*(energy/r-0.5*velSqr)));

	}
}

void CPlasmaVariable::SetSoundSpeed(double GammaMonatomic, double GammaDiatomic) {

	unsigned short loc = 0, iSpecies;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

		double Density = Solution[loc+0];
		double Energy_el = 0.0;

		double Energy = Solution[loc+nDim+1] / Density;
		double Energy_vib;
		if (iSpecies < nDiatomics )
			Energy_vib = Solution[loc+nDim+2] / Density;
		else
			Energy_vib = 0.0;

		if ( iSpecies < nDiatomics )
			SoundSpeed[iSpecies] 	= sqrt(fabs(GammaDiatomic*(GammaDiatomic - 1.0)*(Energy - 0.5*Velocity2[iSpecies] - Energy_vib - Energy_el - Enthalpy_Formation[iSpecies])));
		else {
			SoundSpeed[iSpecies] 	= sqrt(fabs(GammaMonatomic*(GammaMonatomic - 1.0)*(Energy - 0.5*Velocity2[iSpecies] - Energy_el - Enthalpy_Formation[iSpecies])));
		}
	}
}

void CPlasmaVariable::SetPressure(double Gamma ) {

	unsigned short loc = 0, iSpecies, iDim;
	double r, velSqr, energy;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		loc = (nDim+2)*iSpecies;
		r = Solution[loc + 0];
		velSqr = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			velSqr += pow( (Solution[loc + iDim +1]/r),2);
		energy = Solution[loc + nDim + 1];
		Pressure[iSpecies]  = (Gamma - 1.0)*r*(fabs(energy/r-0.5*velSqr));

	}
}

void CPlasmaVariable::SetPressure(double GammaMonatomic, double GammaDiatomic, double *Coord) {

	unsigned short loc = 0, iSpecies, iDim;
	double Density;
	double Energy, Energy_vib, Energy_el;
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
		if ( iSpecies < nDiatomics )
			Pressure[iSpecies] = (GammaDiatomic-1.0) * Density * (Energy - 1.0/2.0*Vel2 - Enthalpy_Formation[iSpecies] - Energy_vib - Energy_el);
		else
			Pressure[iSpecies] = (GammaMonatomic-1.0) * Density * (Energy - 1.0/2.0*Vel2 - Enthalpy_Formation[iSpecies] - Energy_el);
		if (Pressure[iSpecies] < 0.0) {
			cout << "Negative pressure at point: ";
			for (unsigned short iDim = 0; iDim < nDim; iDim++) cout << Coord[iDim] <<" ";
			cout << " in species " << iSpecies;
			cout << endl;
		}
		
	}
}

void CPlasmaVariable::SetTemperature(double* Gas_Constant_MS ) {

	unsigned short loc = 0, iSpecies;
	double r;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		loc = (nDim+2)*iSpecies;
		r = Solution[loc + 0];
		Temperature[iSpecies] = Pressure[iSpecies]/ (Gas_Constant_MS[iSpecies] * r);

	}
}

void CPlasmaVariable::SetTemperature_TR(double* Molar_Mass, double GammaMonatomic, double GammaDiatomic) {
	
	unsigned short loc = 0, iSpecies, iDim;
	double Density;
	double Energy, Energy_vib, Energy_el;
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
			Temperature_tr[iSpecies] =   (GammaDiatomic-1.0)/(UNIVERSAL_GAS_CONSTANT/Molar_Mass[iSpecies]) 
																 * (Energy - 1.0/2.0*Vel2 - Enthalpy_Formation[iSpecies] - Energy_vib - Energy_el);
		} else {
			Temperature_tr[iSpecies] =   (GammaMonatomic-1.0)/(UNIVERSAL_GAS_CONSTANT/Molar_Mass[iSpecies])  
																 * (Energy - 1.0/2.0*Vel2 - Enthalpy_Formation[iSpecies] - Energy_el);
		}
		if (Temperature_tr[iSpecies] < 0.0)
			cout << "Temperature_tr[" << iSpecies << "] is < 0!!!" << endl;
	}
}

void CPlasmaVariable::SetEnthalpy() {
	//COMMENT: Some question as to whether or not to include 1/2 V^2 in the enthalpy...  For now, eliminated from the preprocessing stage

	unsigned short loc = 0, iSpecies, iDim;
	double Density;
	double Energy, Energy_vib;
	double Vel2;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
		Density    = Solution[loc+0];
		Energy     = Solution[loc+nDim+1] / Density;
		Energy_vib = 0.0;
		Vel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Vel2 += Solution[loc+iDim+1]/Solution[loc+0] * Solution[loc+iDim+1]/Solution[loc+0];
//		Enthalpy[iSpecies] = Density*(Energy-1.0/2.0*Vel2) + Pressure[iSpecies];
		Enthalpy[iSpecies] = Energy + Pressure[iSpecies] / Density;
	}
}

double CPlasmaVariable::GetVelocity(unsigned short val_dim,unsigned short val_Species) {

	unsigned short loc = 0;

	if ( val_Species < nDiatomics ) loc = (nDim+3)*val_Species;
	else loc = (nDim+3)*nDiatomics + (nDim+2)*(val_Species-nDiatomics);


	double rho = Solution[loc + 0];
	double rho_vel = Solution[loc + val_dim + 1];
	return (rho_vel/rho);
}

void CPlasmaVariable::SetLaminarViscosity(CConfig *config) {
	double Temperature_Dim;
	unsigned short iSpecies;
	double Temperature_Ref, Viscosity_Ref;

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
}

void CPlasmaVariable ::SetThermalCoeff(double Gamma, double* Gas_Constant_MS ) {
	unsigned short iSpecies;
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
		Kappa[iSpecies] = Gamma * LaminarViscosity_MultiSpecies[iSpecies] / ( (Gamma-1.0)*Gas_Constant_MS[iSpecies]*Prandtl_Lam );
}


void CPlasmaVariable::SetSource(double * res) {

	for (int iVar = 0; iVar < nVar; iVar ++)
		Source[iVar] = res[iVar];
}

void CPlasmaVariable::SetSourceJacobian(double ** src_jacobian) {

	for (int iVar = 0; iVar < nVar; iVar ++) {
		for (int jVar = 0; jVar < nVar; jVar ++)
			Source_Jacobian[iVar][jVar] = src_jacobian[iVar][jVar];
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

CWaveVariable::CWaveVariable(double val_wave, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) 
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;
	
	/*--- Allocate residual structures ---*/
	Res_Visc = new double [nVar]; Res_Sour = new double [nVar];
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];

	/*--- Initialization of variables ---*/
	for (iVar = 0; iVar< nVar; iVar++) {
		Solution[iVar] = val_wave; 
		Solution_Old[iVar] = val_wave;
	}
}

CWaveVariable::~CWaveVariable(void) { 
	
	delete [] Res_Visc; delete [] Res_Sour;
	delete [] Residual_Sum; delete [] Residual_Old;
	
}
