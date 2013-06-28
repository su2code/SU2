/*!
 * \file variable_direct.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.2
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
	unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
	bool Incompressible = config->GetIncompressible();
  bool Freesurface = ((config->GetKind_Solver() == FREE_SURFACE_EULER) ||
                      (config->GetKind_Solver() == FREE_SURFACE_NAVIER_STOKES) ||
                      (config->GetKind_Solver() == FREE_SURFACE_RANS) ||
                      (config->GetKind_Solver() == ADJ_FREE_SURFACE_EULER) ||
                      (config->GetKind_Solver() == ADJ_FREE_SURFACE_NAVIER_STOKES) ||
                      (config->GetKind_Solver() == ADJ_FREE_SURFACE_RANS));
	bool magnet = (config->GetMagnetic_Force() == YES);
  
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
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
		Undivided_Laplacian = new double [nVar];
	if ((config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) &&
			(config->GetKind_SlopeLimit_Flow() != NONE)) {
		Limiter = new double [nVar];
		Solution_Max = new double [nVar];
		Solution_Min = new double [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Limiter[iVar] = 0.0;
      Solution_Max[iVar] = 0.0;
      Solution_Min[iVar] = 0.0;
    }
	}
  
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
		for (iVar = 0; iVar < nVar; iVar++) TS_Source[iVar] = 0.0;
	}
  
  /*--- Allocate auxiliar vector for magnetic field ---*/
  if (magnet) B_Field = new double [3];
  
  /*--- Allocate auxiliar vector for free surface source term ---*/
	if (Freesurface) Grad_AuxVar = new double [nDim];
  
	/*--- Allocate and initialize the primitive variables and gradients ---*/
	if (Incompressible) {
		/*--- Incompressible flow, primitive variables nDim+2, (rho,vx,vy,vz,beta) ---*/
		Primitive = new double [nDim+2];
    for (iVar = 0; iVar < nDim+2; iVar++) Primitive[iVar] = 0.0;
    
    /*--- Compressible flow, gradients primitive variables nDim+2, (rho,vx,vy,vz,beta) ---*/
    Gradient_Primitive = new double* [nDim+2];
    for (iVar = 0; iVar < nDim+2; iVar++) {
      Gradient_Primitive[iVar] = new double [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Gradient_Primitive[iVar][iDim] = 0.0;
    }
	} else {
		/*--- Compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c) ---*/
		Primitive = new double [nDim+5];
    for (iVar = 0; iVar < nDim+5; iVar++) Primitive[iVar] = 0.0;
    
    /*--- Compressible flow, gradients primitive variables nDim+2, (T,vx,vy,vz,P) ---*/
    Gradient_Primitive = new double* [nDim+2];
    for (iVar = 0; iVar < nDim+2; iVar++) {
      Gradient_Primitive[iVar] = new double [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Gradient_Primitive[iVar][iDim] = 0.0;
    }
	}
  
}

CEulerVariable::CEulerVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
	bool Incompressible = config->GetIncompressible();
  bool Freesurface = ((config->GetKind_Solver() == FREE_SURFACE_EULER) ||
                      (config->GetKind_Solver() == FREE_SURFACE_NAVIER_STOKES) ||
                      (config->GetKind_Solver() == FREE_SURFACE_RANS) ||
                      (config->GetKind_Solver() == ADJ_FREE_SURFACE_EULER) ||
                      (config->GetKind_Solver() == ADJ_FREE_SURFACE_NAVIER_STOKES) ||
                      (config->GetKind_Solver() == ADJ_FREE_SURFACE_RANS));
	bool magnet = (config->GetMagnetic_Force() == YES);
  
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
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
		Undivided_Laplacian = new double [nVar];
	if ((config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) &&
			(config->GetKind_SlopeLimit_Flow() != NONE)) {
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
  
  /*--- Allocate auxiliar vector for magnetic field ---*/
  if (magnet) B_Field = new double [3];

  /*--- Allocate auxiliar vector for free surface source term ---*/
	if (Freesurface) Grad_AuxVar = new double [nDim];
  
	/*--- Allocate and initialize the primitive variables and gradients ---*/
	if (Incompressible) {
		/*--- Incompressible flow, primitive variables nDim+2, (rho,vx,vy,vz,beta) ---*/
		Primitive = new double [nDim+2];
    for (iVar = 0; iVar < nDim+2; iVar++) Primitive[iVar] = 0.0;
    
    /*--- Compressible flow, gradients primitive variables nDim+2, (rho,vx,vy,vz,beta) ---*/
    Gradient_Primitive = new double* [nDim+2];
    for (iVar = 0; iVar < nDim+2; iVar++) {
      Gradient_Primitive[iVar] = new double [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Gradient_Primitive[iVar][iDim] = 0.0;
    }
	} else {
		/*--- Compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c) ---*/
		Primitive = new double [nDim+5];
    for (iVar = 0; iVar < nDim+5; iVar++) Primitive[iVar] = 0.0;
    
    /*--- Compressible flow, gradients primitive variables nDim+2, (T,vx,vy,vz,P) ---*/
    Gradient_Primitive = new double* [nDim+2];
    for (iVar = 0; iVar < nDim+2; iVar++) {
      Gradient_Primitive[iVar] = new double [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Gradient_Primitive[iVar][iDim] = 0.0;
    }
	}
  
}

CEulerVariable::~CEulerVariable(void) {
  unsigned short iVar;
  
	delete [] Res_Conv;
  delete [] Res_Visc;
  delete [] Res_Sour;
  delete [] Res_TruncError;

  if (Residual_Sum != NULL) delete [] Residual_Sum;
  if (Residual_Old != NULL) delete [] Residual_Old;
	if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
	if (Limiter != NULL) delete [] Limiter;
  if (Solution_Max != NULL) delete [] Solution_Max;
  if (Solution_Min != NULL) delete [] Solution_Min;
	if (Solution_time_n != NULL) delete [] Solution_time_n;
	if (Solution_time_n1 != NULL) delete [] Solution_time_n1;

	if (Grad_AuxVar != NULL) delete [] Grad_AuxVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
		if (Res_Visc_RK[iVar] != NULL) delete [] Res_Visc_RK[iVar];
	if (Res_Visc_RK != NULL) delete [] Res_Visc_RK;
  
	if (B_Field != NULL) delete [] B_Field;

	if (TS_Source != NULL) delete [] TS_Source;

	delete [] Primitive;
  for (iVar = 0; iVar < nDim+2; iVar++)
    delete [] Gradient_Primitive[iVar];
	delete [] Gradient_Primitive;

}

void CEulerVariable::SetGradient_PrimitiveZero(void) {
	unsigned short iVar, iDim;

	for (iVar = 0; iVar < nDim+2; iVar++)
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
		ProjVel += (Solution[iDim+1]/Primitive[0])*val_vector[iDim];

	return ProjVel;
}

void CEulerVariable::SetPrimVar_Compressible(double Gamma, double Gas_Constant) {
	unsigned short iDim;

	SetVelocity2();									// Compute the modulus of the velocity.
	SetPressure(Gamma);							// Requires Velocity2 computation.
	SetSoundSpeed(Gamma);						// Requires pressure computation.
	SetEnthalpy();									// Requires pressure computation.
	SetTemperature(Gas_Constant);		// Requires pressure computation.

	for (iDim = 0; iDim < nDim; iDim++) 
		Primitive[iDim+1] = Solution[iDim+1] / Solution[0];
	Primitive[nDim+2] = Solution[0];

}

void CEulerVariable::SetPrimVar_Incompressible(double Density_Inf, double ArtComp_Factor, bool freesurface) {
	unsigned short iDim;

	SetBetaInc2(ArtComp_Factor);	// Set the value of beta.
	if (!freesurface) 
		SetDensityInc(Density_Inf);	// Set the value of the density
	SetVelocityInc2();						// Requires density computation.

	for (iDim = 0; iDim < nDim; iDim++) 
		Primitive[iDim+1] = Solution[iDim+1] / Primitive[0];

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

void CNSVariable::SetLaminarViscosity() {
	double Temperature_Dim;

	/*--- Calculate viscosity from a non-dim. Sutherland's Law ---*/
	Temperature_Dim = Primitive[0]*Temperature_Ref;
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

void CNSVariable::SetPrimVar_Compressible(double Gamma, double Gas_Constant, double turb_ke) {
	unsigned short iDim;

	SetVelocity2();								 // Compute the modulus of the velocity.
	SetPressure(Gamma, turb_ke);   // Requires Velocity2 computation.
	SetSoundSpeed(Gamma);					 // Requires pressure computation.
	SetEnthalpy();								 // Requires pressure computation.
	SetTemperature(Gas_Constant);  // Requires pressure computation.
	SetLaminarViscosity();				 // Requires temperature computation.

	for (iDim = 0; iDim < nDim; iDim++)
		Primitive[iDim+1] = Solution[iDim+1] / Solution[0];
	Primitive[nDim+2] = Solution[0];

}

void CNSVariable::SetPrimVar_Incompressible(double Density_Inf, double Viscosity_Inf, double ArtComp_Factor, double turb_ke, bool freesurface) {
	unsigned short iDim;


	SetBetaInc2(ArtComp_Factor);	// Set the value of beta.
	if (!freesurface) {
		SetDensityInc(Density_Inf);							// Set the value of the density
		SetLaminarViscosityInc(Viscosity_Inf);	// Set the value of the viscosity
	}
	SetVelocityInc2();						// Requires density computation.

	for (iDim = 0; iDim < nDim; iDim++) 
		Primitive[iDim+1] = Solution[iDim+1] / GetDensityInc();

}

CTurbVariable::CTurbVariable(void) : CVariable() { }

CTurbVariable::CTurbVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) { }

CTurbVariable::~CTurbVariable(void) { }

double CTurbVariable::GetmuT(){ return muT; }

void CTurbVariable::SetmuT(double val_muT){ muT = val_muT; }

CTurbSAVariable::CTurbSAVariable(void) : CTurbVariable() { }

CTurbSAVariable::CTurbSAVariable(double val_nu_tilde, double val_muT, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

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
		for (iVar = 0; iVar < nVar; iVar++)
			TS_Source[iVar] = 0.0;
	}

	/*--- Allocate space for the limiter ---*/
	if (config->GetKind_SlopeLimit_Turb() != NONE) {
		Limiter = new double [nVar];
		Solution_Max = new double [nVar];
		Solution_Min = new double [nVar];
	}

}

CTurbSAVariable::~CTurbSAVariable(void) {
	if (Limiter != NULL) delete [] Limiter;
	if (Solution_Max != NULL) delete [] Solution_Max;
	if (Solution_Min != NULL) delete [] Solution_Min;
	if (TS_Source != NULL) delete [] TS_Source;
}

CTransLMVariable::CTransLMVariable(void) : CTurbVariable() {}

CTransLMVariable::CTransLMVariable(double val_nu_tilde, double val_intermittency, double val_REth,  unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_ndim, val_nvar,config) {

	// Initialization of variables
	Solution[0] = val_intermittency; Solution_Old[0] = val_intermittency;
	Solution[1] = val_REth;          Solution_Old[1] = val_REth;

}

CTransLMVariable::~CTransLMVariable(void) { }

void CTransLMVariable::SetGammaEff() {

	/* -- Correction for separation-induced transition -- */
	Solution[0] = max(Solution[0],gamma_sep);

}

CTurbSSTVariable::CTurbSSTVariable(void) : CTurbVariable() { }

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

CTurbSSTVariable::~CTurbSSTVariable(void) { }

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

	unsigned short iVar, loc = 0, iSpecies, iDim, iMesh, nMGSmooth = 0;
	double rho, energy,	energy_vib, energy_el, enthalpy_formation, Gamma;	
	bool viscous;
	if ((config->GetKind_Solver() == PLASMA_NAVIER_STOKES) || (config->GetKind_Solver() == ADJ_PLASMA_NAVIER_STOKES)) viscous = true;

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
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTERED) {
		Undivided_Laplacian = new double [nVar];
		Lambda = new double [nSpecies];
	}
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

	Solution_Max = new double[nVar];
	Solution_Min = new double[nVar];

	/*--- Allocate and initialize the local useful quantities ---*/
	Velocity    = new double*[nSpecies];
	Velocity2   = new double [nSpecies];
	Pressure_Old = new double [nSpecies];
	Max_Lambda_Inv_MultiSpecies = new double [nSpecies];
	Species_Delta_Time = new double [nSpecies];
	Sensor_MultiSpecies = new double [nSpecies];

	/*--- Compressible flow, primitive variables nDim+6, (T_tr,vx,vy,vz,T_vi, P, rho, h, c) ---*/
	Primitive = new double*[nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
		Primitive[iSpecies] = new double [nDim+6];

	/*--- Allocate and initialize the primitive gradient variable nDim+2, (T_tr,vx,vy,vz,T_vi) ---*/
	Gradient_Primitive = new double **[nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Gradient_Primitive[iSpecies] = new double *[nDim+3];
		for (iVar = 0; iVar < nDim+3; iVar++)
			Gradient_Primitive[iSpecies][iVar] = new double [nDim];
	}	

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
		Primitive[iSpecies][nDim+2] 	= (Gamma - 1.0)*rho*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - enthalpy_formation);
		Primitive[iSpecies][nDim+5] 	= sqrt(fabs(Gamma*(Gamma - 1.0)*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - enthalpy_formation)));		
		Primitive[iSpecies][nDim+4]		= energy + Primitive[iSpecies][nDim+2]/rho;
		Max_Lambda_Inv_MultiSpecies[iSpecies] = 0.0;
		Species_Delta_Time[iSpecies] = 1.0;
	}


	//	if (viscous) {
	LaminarViscosity_MultiSpecies = new double [nSpecies];
	EddyViscosity_MultiSpecies = new double [nSpecies];
	ThermalCoeff = new double [nSpecies];
	ThermalCoeff_vib = new double[nSpecies];
	Max_Lambda_Visc_MultiSpecies = new double [nSpecies];

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		LaminarViscosity_MultiSpecies[iSpecies] = 0.0;
		EddyViscosity_MultiSpecies[iSpecies] = 0.0;
		ThermalCoeff[iSpecies] = 0.0;
		ThermalCoeff_vib[iSpecies] = 0.0;
		Max_Lambda_Visc_MultiSpecies[iSpecies] = 0.0;
	}
	//	}


	//	if (config->GetElectricSolver()) {
	Elec_Field = new double [nDim];
	for(iDim =0; iDim < nDim; iDim ++)
		Elec_Field[iDim] = 0.0;
	//	}


	B_Field = new double [nDim];
	for(iDim =0; iDim < nDim; iDim ++) {
		Elec_Field[iDim] = 0.0;
		B_Field[iDim] = 0.0;
	}

}
CPlasmaVariable::CPlasmaVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {

	unsigned short iVar, loc = 0, iSpecies, iDim, iMesh, nMGSmooth = 0;
	double rho, energy, energy_vib, energy_el, Gamma, enthalpy_formation;
	bool viscous;
	if ((config->GetKind_Solver() == PLASMA_NAVIER_STOKES) || (config->GetKind_Solver() == ADJ_PLASMA_NAVIER_STOKES)) viscous = true;
	nMonatomics  = config->GetnMonatomics();
	nDiatomics = config->GetnDiatomics();
	nSpecies = config->GetnSpecies();
	Prandtl_Lam     = config->GetPrandtl_Lam();

	/*--- Allocate residual structures ---*/
	Res_Conv = new double [nVar];
	Res_Visc = new double [nVar];
	Res_Sour = new double [nVar];
	Res_TruncError = new double [nVar];

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
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTERED) {
		Undivided_Laplacian = new double [nVar];
		Lambda = new double [nSpecies];
	}
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


	/*--- Allocate and initialize the local useful quantities ---*/
	Velocity    = new double*[nSpecies];
	Velocity2   = new double [nSpecies];
	Pressure_Old = new double [nSpecies];
	Max_Lambda_Inv_MultiSpecies = new double [nSpecies];
	Species_Delta_Time = new double [nSpecies];
	Sensor_MultiSpecies = new double [nSpecies];

	/*--- Compressible flow, primitive variables nDim+6, (T_tr,vx,vy,vz,T_vi, P, rho, h, c) ---*/
	Primitive = new double*[nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
		Primitive[iSpecies] = new double [nDim+6];

	/*--- Allocate and initialize the primitive gradient variable nDim+2, (T_tr,vx,vy,vz,T_vi) ---*/
	Gradient_Primitive = new double **[nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Gradient_Primitive[iSpecies] = new double *[nDim+3];
		for (iVar = 0; iVar < nDim+3; iVar++)
			Gradient_Primitive[iSpecies][iVar] = new double [nDim];
	}

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
			Velocity2[iSpecies] += Solution[loc + iDim+1]/rho*Solution[loc + iDim+1]/rho;
		}
		Gamma = config->GetSpecies_Gamma(iSpecies);
		enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Primitive[iSpecies][nDim+2]     = (Gamma - 1.0)*rho*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - enthalpy_formation);
		Primitive[iSpecies][nDim+5]     = sqrt(fabs(Gamma*(Gamma - 1.0)*(energy/rho - 0.5*Velocity2[iSpecies] - energy_vib/rho - energy_el - enthalpy_formation)));
		Primitive[iSpecies][nDim+4]             = energy + Primitive[iSpecies][nDim+2]/rho;
		Max_Lambda_Inv_MultiSpecies[iSpecies] = 0.0;
		Species_Delta_Time[iSpecies] = 1.0;
	}


	//      if (viscous) {
		LaminarViscosity_MultiSpecies = new double [nSpecies];
		EddyViscosity_MultiSpecies = new double [nSpecies];
		ThermalCoeff = new double [nSpecies];
		ThermalCoeff_vib = new double[nSpecies];
		Max_Lambda_Visc_MultiSpecies = new double [nSpecies];

		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			LaminarViscosity_MultiSpecies[iSpecies] = 0.0;
			EddyViscosity_MultiSpecies[iSpecies] = 0.0;
			ThermalCoeff[iSpecies] = 0.0;
			ThermalCoeff_vib[iSpecies] = 0.0;
			Max_Lambda_Visc_MultiSpecies[iSpecies] = 0.0;
		}
		//      }


		//      if (config->GetElectricSolver()) {
		Elec_Field = new double [nDim];
		for(iDim =0; iDim < nDim; iDim ++)
			Elec_Field[iDim] = 0.0;
		//      }


		B_Field = new double [nDim];
		for(iDim =0; iDim < nDim; iDim ++) {
			Elec_Field[iDim] = 0.0;
			B_Field[iDim] = 0.0;
		}

}







CPlasmaVariable::~CPlasmaVariable(void) {
	unsigned short iVar, iSpecies;

	if (Elec_Field                    !=NULL ) delete [] Elec_Field;
	if (LaminarViscosity_MultiSpecies !=NULL ) delete [] LaminarViscosity_MultiSpecies;
	if (EddyViscosity_MultiSpecies    !=NULL ) delete [] EddyViscosity_MultiSpecies;
	if (ThermalCoeff                  !=NULL ) delete [] ThermalCoeff;
	if (ThermalCoeff_vib              !=NULL ) delete [] ThermalCoeff_vib;
	if (Max_Lambda_Visc_MultiSpecies  !=NULL ) delete [] Max_Lambda_Visc_MultiSpecies;
	if (Lambda                        !=NULL ) delete [] Lambda;

	delete [] Res_Conv; delete [] Res_Visc; delete [] Res_Sour; delete [] Res_TruncError;
	delete [] Residual_Sum; delete [] Residual_Old;
	if (Undivided_Laplacian != NULL) delete [] Undivided_Laplacian;
	if (Limiter != NULL) delete [] Limiter;
	if (Solution_time_n != NULL) delete [] Solution_time_n;
	if (Solution_time_n1 != NULL) delete [] Solution_time_n1;

	delete [] Max_Lambda_Inv_MultiSpecies;
	delete [] Max_Lambda_Visc_MultiSpecies;
	delete [] Lambda;
	delete [] Pressure_Old;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Res_Visc_RK[iVar];
	delete [] Res_Visc_RK;

	for (iSpecies = 0; iSpecies< nSpecies; iSpecies++) {
		for (iVar = 0; iVar < nDim+3; iVar++) {
			delete [] Gradient_Primitive[iSpecies][iVar];
		}
		delete [] Gradient_Primitive[iSpecies];
	}
	delete [] Gradient_Primitive;

	delete [] B_Field;

	for (iVar = 0; iVar < nSpecies; iVar++) {
		delete [] Velocity[iVar];
	}
	delete [] Velocity;
	delete [] Velocity2;

	delete [] Max_Lambda_Inv_MultiSpecies;
	delete [] Max_Lambda_Visc_MultiSpecies;
	delete [] Lambda;
	delete [] Species_Delta_Time;
	delete [] Sensor_MultiSpecies;
	delete [] Solution_Max;
	delete [] Solution_Min;
}

void CPlasmaVariable::SetVel_ResConv_Zero(unsigned short iSpecies) {
	unsigned short loc, iDim;
	if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
	else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

	for (iDim = 0; iDim < nDim; iDim++)
		Res_Conv[loc+iDim+1] = 0.0;
}

void CPlasmaVariable::SetVel_ResVisc_Zero(unsigned short iSpecies) {
	unsigned short loc, iDim;
	if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
	else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

	for (iDim = 0; iDim < nDim; iDim++)
		Res_Visc[loc+iDim+1] = 0.0;
}

void CPlasmaVariable::SetVel_ResSour_Zero(unsigned short iSpecies) {
	unsigned short loc, iDim;
	if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
	else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

	for (iDim = 0; iDim < nDim; iDim++)
		Res_Sour[loc+iDim+1] = 0.0;
}

void CPlasmaVariable::SetVelRes_TruncErrorZero(unsigned short iSpecies) {
	unsigned short loc, iDim;
	if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
	else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

	for (iDim = 0; iDim < nDim; iDim++)
		Res_TruncError[loc+iDim+1] = 0.0;
}

void CPlasmaVariable::SetGradient_PrimitiveZero(void) {
	unsigned short iSpecies, iVar, iDim;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
		for (iVar = 0; iVar < nDim+3; iVar++)
			for (iDim = 0; iDim < nDim; iDim++)
				Gradient_Primitive[iSpecies][iVar][iDim] = 0.0;

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
		Primitive[iSpecies][nDim+5] = sqrt(Gamma*(Gamma - 1.0)*(Energy - 0.5*Velocity2[iSpecies] - Energy_vib - Energy_el - Enthalpy_formation));
	}
}

void CPlasmaVariable::SetPressure(CConfig *config) {

	unsigned short loc = 0, iSpecies, iDim;
	double Gamma, Enthalpy_formation;
	double Density, Energy, Energy_vib, Energy_el;
	double Vel2;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
		Pressure_Old[iSpecies] = Primitive[iSpecies][nDim+2];
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
		Primitive[iSpecies][nDim+2] = (Gamma-1.0) * Density * (Energy - 0.5*Vel2 - Enthalpy_formation - Energy_vib - Energy_el);

		if (Primitive[iSpecies][nDim+2] < 0.0) {
			Primitive[iSpecies][nDim+2] = Pressure_Old[iSpecies];
			cout << "Density = " << Density << endl;
			cout << "energy = " << Energy << endl;
			cout << "vel2  = " << Vel2 << endl;
			cout << "Negative pressure of " << Primitive[iSpecies][nDim+2] << " in species " << iSpecies;
			cout << endl;
		}

	}
}

void CPlasmaVariable::SetTemperature_tr(CConfig *config) {

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

		Primitive[iSpecies][0] = (Gamma-1.0)/Gas_constant
				* (Energy - 1.0/2.0*Vel2 - Enthalpy_formation - Energy_vib - Energy_el);
	}
}

void CPlasmaVariable::SetTemperature_vib(CConfig *config) {
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

			Primitive[iSpecies][nDim+1] = CharVibTemp / log(CharVibTemp*UNIVERSAL_GAS_CONSTANT/(Energy_vib*MolarMass) + 1.0);

		}
		else {
			Primitive[iSpecies][nDim+1] = 0.0;
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
		Primitive[iSpecies][nDim+4] = Energy + Primitive[iSpecies][nDim+2] / Density;
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

void CPlasmaVariable::SetVelocity_Old(double *val_velocity, unsigned short iSpecies) {
	unsigned short loc, iDim;

	if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
	else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

	for (iDim = 0; iDim < nDim; iDim++)
		Solution_Old[loc+iDim+1] = val_velocity[iDim]*Solution[loc+0];
}

void CPlasmaVariable::SetPrimVar(CConfig *config, double *Coord) {
	unsigned short iDim, iSpecies, loc;

	bool viscous = (config->GetKind_ViscNumScheme_Plasma());

	double Gamma, Enthalpy_formation, Gas_constant, Molar_mass, Char_temp_vib;
	double density, energy_tot, energy_vib, energy_el;

	/*--- Primitive variable vector: (Ttr, vx, vy, vz, Tvib, P, rho, h, c) ---*/

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) {
			loc = (nDim+3)*iSpecies;
			energy_vib = Solution[loc+nDim+2]/Solution[loc];
		}	else {
			loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
			energy_vib = 0.0;
		}

		/*--- Store useful quantities ---*/
		density    = Solution[loc+0];
		energy_tot = Solution[loc+nDim+1]/density;
		energy_el  = 0.0;

		Gamma = config->GetSpecies_Gamma(iSpecies);
		Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Molar_mass = config->GetMolar_Mass(iSpecies);

		/*--- Velocity ---*/
		Velocity2[iSpecies] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity2[iSpecies] += Solution[loc+iDim+1]*Solution[loc+iDim+1]/(density*density);
			Primitive[iSpecies][iDim+1] = Solution[loc+iDim+1] / Solution[loc+0];
		}

		/*--- Density ---*/
		Primitive[iSpecies][nDim+3] = Solution[loc+0];

		/*--- Sound speed ---*/
		Primitive[iSpecies][nDim+5] = sqrt(Gamma*(Gamma - 1.0)*(energy_tot - 0.5*Velocity2[iSpecies] - energy_vib - energy_el - Enthalpy_formation));

		/*--- Pressure ---*/
		Primitive[iSpecies][nDim+2] = (Gamma-1.0) * density * (energy_tot - 0.5*Velocity2[iSpecies] - Enthalpy_formation - energy_vib - energy_el);

		/*--- Translational-rotational temperature ---*/
		Primitive[iSpecies][0] = (Gamma-1.0)/Gas_constant * (energy_tot - 0.5*Velocity2[iSpecies] - Enthalpy_formation - energy_vib - energy_el);

		/*--- Vibrational temperature ---*/
		if (iSpecies < nDiatomics) {
			Char_temp_vib = config->GetCharVibTemp(iSpecies);
			Primitive[iSpecies][nDim+1] = Char_temp_vib / log(Char_temp_vib*UNIVERSAL_GAS_CONSTANT/(energy_vib*Molar_mass) + 1.0);
		} else {
			Primitive[iSpecies][nDim+1] = 0.0;
		}

		/*--- Enthalpy ---*/
		Primitive[iSpecies][nDim+4] = energy_tot + Primitive[iSpecies][nDim+2]/density;
    
    
    //ERROR CHECKING
    if (Primitive[iSpecies][0] < 0) {
      cout << "Temp TR[" << iSpecies << "] is is negative! " << Primitive[iSpecies][0] << endl;
      cout << "Coordinate: (" << Coord[0] << ", " << Coord[1] << ", "  << Coord[2] << ")" << endl;
      Primitive[iSpecies][0] = EPS;
      cin.get();
    }
    if (Primitive[iSpecies][nDim+1] < 0) {
      cout << "Temp vib[" << iSpecies << "] is is negative! " << Primitive[iSpecies][nDim+1] << endl;
      Primitive[iSpecies][nDim+1] = EPS;
      cin.get();
    }
    if (Primitive[iSpecies][nDim+2] < 0) {
      cout << "Pressure[" << iSpecies << "] is is negative! " << Primitive[iSpecies][nDim+2] << endl;
      Primitive[iSpecies][nDim+2] = EPS;
      cin.get();
    }
    if (Primitive[iSpecies][nDim+3] < 0) {
      cout << "Density[" << iSpecies << "] is negative! " << Primitive[iSpecies][nDim+3] << endl;
      Primitive[iSpecies][nDim+3] = EPS;
      cin.get();
    }
	}

	if(viscous) {
		SetLaminarViscosity(config);
		SetThermalCoeff(config);
	}
}


void CPlasmaVariable::SetPrimVar(CConfig *config) {
	unsigned short iDim, iSpecies, loc;
  
	bool viscous = (config->GetKind_ViscNumScheme_Plasma());
  
	double Gamma, Enthalpy_formation, Gas_constant, Molar_mass, Char_temp_vib;
	double density, energy_tot, energy_vib, energy_el;
  
	/*--- Primitive variable vector: (Ttr, vx, vy, vz, Tvib, P, rho, h, c) ---*/
  
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) {
			loc = (nDim+3)*iSpecies;
			energy_vib = Solution[loc+nDim+2]/Solution[loc];
		}	else {
			loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
			energy_vib = 0.0;
		}
    
		/*--- Store useful quantities ---*/
		density    = Solution[loc+0];
		energy_tot = Solution[loc+nDim+1]/density;
		energy_el  = 0.0;
    
		Gamma = config->GetSpecies_Gamma(iSpecies);
		Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Molar_mass = config->GetMolar_Mass(iSpecies);
    
		/*--- Velocity ---*/
		Velocity2[iSpecies] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity2[iSpecies] += Solution[loc+iDim+1]*Solution[loc+iDim+1]/(density*density);
			Primitive[iSpecies][iDim+1] = Solution[loc+iDim+1] / Solution[loc+0];
		}
    
		/*--- Density ---*/
		Primitive[iSpecies][nDim+3] = Solution[loc+0];
    
		/*--- Sound speed ---*/
		Primitive[iSpecies][nDim+5] = sqrt(Gamma*(Gamma - 1.0)*(energy_tot - 0.5*Velocity2[iSpecies] - energy_vib - energy_el - Enthalpy_formation));
    
		/*--- Pressure ---*/
		Primitive[iSpecies][nDim+2] = (Gamma-1.0) * density * (energy_tot - 0.5*Velocity2[iSpecies] - Enthalpy_formation - energy_vib - energy_el);
    
		/*--- Translational-rotational temperature ---*/
		Primitive[iSpecies][0] = (Gamma-1.0)/Gas_constant * (energy_tot - 0.5*Velocity2[iSpecies] - Enthalpy_formation - energy_vib - energy_el);
    
		/*--- Vibrational temperature ---*/
		if (iSpecies < nDiatomics) {
			Char_temp_vib = config->GetCharVibTemp(iSpecies);
			Primitive[iSpecies][nDim+1] = Char_temp_vib / log(Char_temp_vib*UNIVERSAL_GAS_CONSTANT/(energy_vib*Molar_mass) + 1.0);
		} else {
			Primitive[iSpecies][nDim+1] = 0.0;
		}
    
		/*--- Enthalpy ---*/
		Primitive[iSpecies][nDim+4] = energy_tot + Primitive[iSpecies][nDim+2]/density;
  
	}
  
	if(viscous) {
		SetLaminarViscosity(config);
		SetThermalCoeff(config);
	}
}



void CPlasmaVariable::SetLaminarViscosity(CConfig *config) {
	double Temperature_Dim;
	unsigned short iSpecies;
	double Temperature_Ref, Viscosity_Ref;
	double As, Bs, Cs;

	switch (config->GetKind_GasModel()) {
	case ARGON:
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Temperature_Ref	= config->GetTemperature_Ref(iSpecies);
			Viscosity_Ref = config->GetViscosity_Ref(iSpecies);
			Temperature_Dim = Primitive[iSpecies][0]*Temperature_Ref;
			LaminarViscosity_MultiSpecies[iSpecies] = 2.3E-5*(pow(Temperature_Dim/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_Dim+110.3));
			LaminarViscosity_MultiSpecies[iSpecies] = LaminarViscosity_MultiSpecies[iSpecies]/Viscosity_Ref;
		}
		iSpecies = nSpecies - 1;
		Temperature_Ref	= config->GetTemperature_Ref(iSpecies);
		Viscosity_Ref = config->GetViscosity_Ref(iSpecies);
		Temperature_Dim = Primitive[iSpecies][0]*Temperature_Ref;
		LaminarViscosity_MultiSpecies[iSpecies] = 2.3E-5*1.3E-7*(pow(Temperature_Dim/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_Dim+110.3));
		LaminarViscosity_MultiSpecies[iSpecies] = LaminarViscosity_MultiSpecies[iSpecies]/Viscosity_Ref;
		break;


	case N2:
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			/*--- Retrieve flow quantities ---*/
			Temperature_Dim = Primitive[iSpecies][0];

			/*--- Calculate species viscosity Blottner et. al. (1971) model ---*/
			As = config->GetBlottnerCoeff(iSpecies,0);
			Bs = config->GetBlottnerCoeff(iSpecies,1);
			Cs = config->GetBlottnerCoeff(iSpecies,2);

			/*--- Calculate the laminar viscosity using the correlation (returns in kg/m*s )---*/
			LaminarViscosity_MultiSpecies[iSpecies] = 0.1*exp( (As*log(Temperature_Dim) + Bs) * log(Temperature_Dim) + Cs );
		}
		break;
	}
}


void CPlasmaVariable ::SetThermalCoeff(CConfig *config) {
	unsigned short iSpecies;
	double Gamma, Gas_constant, cv_t, cv_rot, cv_vib;

	if (config->GetKind_GasModel() == ARGON) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Gamma = config->GetSpecies_Gamma(iSpecies);
			Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);
			ThermalCoeff[iSpecies] = Gamma * LaminarViscosity_MultiSpecies[iSpecies] *Gas_constant/ ( (Gamma-1.0)*Prandtl_Lam );
		}
	} else {
		/*--- Set thermal conductivity on a species-by-species basis ---*/
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
			Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);

			/*--- Set the specific heat at constant volume based available energy storage modes in the species ---*/
			if (iSpecies < nDiatomics) {
				cv_t   = 3.0/2.0 * Gas_constant;
				cv_rot = Gas_constant;
				cv_vib = Gas_constant;
			} else {
				cv_t   = 3.0/2.0 * Gas_constant;
				cv_rot = 0.0;
				cv_vib = 0.0;
			}

			/*--- Set the thermal conductivity according to Eucken's relation ---*/
			ThermalCoeff[iSpecies] = LaminarViscosity_MultiSpecies[iSpecies]*(5.0/2.0*cv_t + cv_rot);
			ThermalCoeff_vib[iSpecies] = LaminarViscosity_MultiSpecies[iSpecies]*cv_vib;
		}
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
	Limiter = new double [nVar];
	Solution_Max = new double [nVar];
	Solution_Min = new double [nVar];

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
	Limiter = new double [nVar];
	Solution_Max = new double [nVar];
	Solution_Min = new double [nVar];

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
	if (Solution_Max != NULL) delete [] Solution_Max;
	if (Solution_Min != NULL) delete [] Solution_Min;

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

