/*!
 * \file variable_direct.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.5
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

CPotentialVariable::CPotentialVariable(void) : CVariable() {

  /*--- Array initialization ---*/
	Charge_Density = NULL;
	PlasmaRhoUGradient = NULL;

}

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
  unsigned short iVar;
  
	if (Charge_Density != NULL) delete [] Charge_Density;
  if (PlasmaRhoUGradient != NULL) {
    for (iVar = 0; iVar < 3; iVar++)
      delete PlasmaRhoUGradient[iVar];
    delete [] PlasmaRhoUGradient;
  }

}

CEulerVariable::CEulerVariable(void) : CVariable() {

  /*--- Array initialization ---*/
	TS_Source = NULL;
  B_Field = NULL;
	Primitive = NULL;
	Gradient_Primitive = NULL;
	Limiter_Primitive = NULL;

}

CEulerVariable::CEulerVariable(double val_density, double *val_velocity, double val_energy, unsigned short val_ndim,
		unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim, iMesh, nMGSmooth = 0;

	bool Incompressible = config->GetIncompressible();
	bool freesurface = config->GetFreeSurface();
	bool magnet = (config->GetMagnetic_Force() == YES);
  bool low_fidelity = config->GetLowFidelitySim();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  /*--- Array initialization ---*/
	TS_Source = NULL;
  B_Field = NULL;
	Primitive = NULL;
	Gradient_Primitive = NULL;
	Limiter_Primitive = NULL;
  
  /*--- Allocate and initialize the primitive variables and gradients ---*/
  if (Incompressible) { nPrimVar = nDim+2; nPrimVarGrad = nDim+2; }
  else { nPrimVar = nDim+5; nPrimVarGrad = nDim+3; }
  
	/*--- Allocate residual structures ---*/
	Res_TruncError = new double [nVar];

	for (iVar = 0; iVar < nVar; iVar++) {
		Res_TruncError[iVar] = 0.0;
	}

	/*--- Only for residual smoothing (multigrid) ---*/
	for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++)
		nMGSmooth += config->GetMG_CorrecSmooth(iMesh);

	if ((nMGSmooth > 0) || low_fidelity || freesurface) {
		Residual_Sum = new double [nVar];
		Residual_Old = new double [nVar];
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
	if (dual_time) {
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
	if (freesurface) Grad_AuxVar = new double [nDim];
  
  /*--- Incompressible flow, primitive variables nDim+2, (rho,vx,vy,vz,beta),
   compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c) ---*/
  Primitive = new double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;
  
  /*--- Incompressible flow, gradients primitive variables nDim+2, (rho,vx,vy,vz,beta), 
   compressible flow, gradients primitive variables nDim+3, (T,vx,vy,vz,P,rho)
   We need P, and rho for running the adjoint problem ---*/
  Gradient_Primitive = new double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }


}

CEulerVariable::CEulerVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim, iMesh, nMGSmooth = 0;

	bool Incompressible = config->GetIncompressible();
	bool freesurface = config->GetFreeSurface();
	bool magnet = (config->GetMagnetic_Force() == YES);
  bool low_fidelity = config->GetLowFidelitySim();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Array initialization ---*/
	TS_Source = NULL;
  B_Field = NULL;
	Primitive = NULL;
	Gradient_Primitive = NULL;
  Limiter_Primitive = NULL;

	/*--- Allocate residual structures ---*/
	Res_TruncError = new double [nVar];

	for (iVar = 0; iVar < nVar; iVar++) {
		Res_TruncError[iVar] = 0.0;
	}

	/*--- Only for residual smoothing (multigrid) ---*/
	for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++)
		nMGSmooth += config->GetMG_CorrecSmooth(iMesh);

	if ((nMGSmooth > 0) || low_fidelity) {
		Residual_Sum = new double [nVar];
		Residual_Old = new double [nVar];
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
	if (dual_time) {
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
	if (freesurface) Grad_AuxVar = new double [nDim];

	/*--- Allocate and initialize the primitive variables and gradients ---*/
  if (Incompressible) { nPrimVar = nDim+2; nPrimVarGrad = nDim+2; }
  else { nPrimVar = nDim+5; nPrimVarGrad = nDim+3; }
  
  /*--- Incompressible flow, primitive variables nDim+2, (rho,vx,vy,vz,beta)
   compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c) ---*/
  Primitive = new double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;
  
  /*--- Incompressible flow, gradients primitive variables nDim+2, (rho,vx,vy,vz,beta)
   Compressible flow, gradients primitive variables nDim+3, (T,vx,vy,vz,P,rho)
   We need P, and rho for running the adjoint problem ---*/
  Gradient_Primitive = new double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }


}

CEulerVariable::~CEulerVariable(void) {
	unsigned short iVar;

	if (B_Field           != NULL) delete [] B_Field;
	if (TS_Source         != NULL) delete [] TS_Source;
  if (Primitive         != NULL) delete [] Primitive;
  if (Limiter_Primitive != NULL) delete [] Limiter_Primitive;

  if (Gradient_Primitive != NULL) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      delete Gradient_Primitive[iVar];
    delete [] Gradient_Primitive;
  }

}

void CEulerVariable::SetGradient_PrimitiveZero(unsigned short val_primvar) {
	unsigned short iVar, iDim;

	for (iVar = 0; iVar < val_primvar; iVar++)
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
		ProjVel += Primitive[iDim+1]*val_vector[iDim];

	return ProjVel;
}

void CEulerVariable::SetPrimVar_Compressible(CConfig *config) {
	unsigned short iDim, iVar;
  bool check_dens = false, check_press = false, check_sos = false, check_temp = false;
  
  double Gas_Constant = config->GetGas_ConstantND();
	double Gamma = config->GetGamma();
  
	SetVelocity2();                               // Compute the modulus of the velocity.
    check_dens = (Solution[0] < 0.0);             // Check the density
	check_press = SetPressure(Gamma);							// Requires Velocity2 computation.
	check_sos = SetSoundSpeed(Gamma);             // Requires pressure computation.
	check_temp = SetTemperature(Gas_Constant);		// Requires pressure computation.
    
    /*--- Check that the solution has a physical meaning ---*/
    if (check_dens || check_press || check_sos || check_temp) {
        
        /*--- Copy the old solution ---*/
        for (iVar = 0; iVar < nVar; iVar++)
            Solution[iVar] = Solution_Old[iVar];
        
        /*--- Recompute the primitive variables ---*/
        SetVelocity2();
        check_press = SetPressure(Gamma);
        check_sos = SetSoundSpeed(Gamma);
        check_temp = SetTemperature(Gas_Constant);
        
    }
    
    SetEnthalpy();                                // Requires pressure computation.
    
	for (iDim = 0; iDim < nDim; iDim++)
		Primitive[iDim+1] = Solution[iDim+1] / Solution[0];
	Primitive[nDim+2] = Solution[0];
    
}

void CEulerVariable::SetPrimVar_Incompressible(double Density_Inf, double levelset, CConfig *config) {
	unsigned short iDim;
  double epsilon, Heaviside, lambda, DensityInc;
  
  double ArtComp_Factor = config->GetArtComp_Factor();
  bool freesurface = config->GetFreeSurface();

  /*--- Set the value of the density ---*/
	if (!freesurface) {
		SetDensityInc(Density_Inf);
  }
  else {
    epsilon = config->GetFreeSurface_Thickness();
    Heaviside = 0.0;
    if (levelset < -epsilon) Heaviside = 1.0;
    if (fabs(levelset) <= epsilon) Heaviside = 1.0 - (0.5*(1.0+(levelset/epsilon)+(1.0/PI_NUMBER)*sin(PI_NUMBER*levelset/epsilon)));
    if (levelset > epsilon) Heaviside = 0.0;

    lambda = config->GetRatioDensity();
    DensityInc = (lambda + (1.0 - lambda)*Heaviside)*config->GetDensity_FreeStreamND();
		SetDensityInc(DensityInc);
  }

  /*--- Set the value of the velocity squared (requires density) ---*/
	SetVelocityInc2();
  
  /*--- Set the value of the artificial compressibility factor ---*/
  SetBetaInc2(ArtComp_Factor);

  /*--- Set the value of the velocity ---*/
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

CNSVariable::~CNSVariable(void) { }

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
	case SA : case SST :
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

void CNSVariable::SetPrimVar_Compressible(CConfig *config, double turb_ke) {
	unsigned short iDim, iVar;
    bool check_dens = false, check_press = false, check_sos = false, check_temp = false;
  
  double Gas_Constant = config->GetGas_ConstantND();
	double Gamma = config->GetGamma();
  
	SetVelocity2();                                 // Compute the modulus of the velocity.
    check_dens = (Solution[0] < 0.0);               // Check the density
	check_press = SetPressure(Gamma, turb_ke);      // Requires Velocity2 computation.
	check_sos = SetSoundSpeed(Gamma);               // Requires pressure computation.
	check_temp = SetTemperature(Gas_Constant);      // Requires pressure computation.
    
    /*--- Check that the solution has a physical meaning ---*/
    if (check_dens || check_press || check_sos || check_temp) {
        
        /*--- Copy the old solution ---*/
        for (iVar = 0; iVar < nVar; iVar++)
            Solution[iVar] = Solution_Old[iVar];
        
        /*--- Recompute the primitive variables ---*/
        SetVelocity2();
        check_press = SetPressure(Gamma, turb_ke);
        check_sos = SetSoundSpeed(Gamma);
        check_temp = SetTemperature(Gas_Constant);
        
    }
    
	SetEnthalpy();                                  // Requires pressure computation.
	SetLaminarViscosity();                          // Requires temperature computation.
    
	for (iDim = 0; iDim < nDim; iDim++)
		Primitive[iDim+1] = Solution[iDim+1] / Solution[0];
	Primitive[nDim+2] = Solution[0];
    
}

void CNSVariable::SetPrimVar_Incompressible(double Density_Inf, double Viscosity_Inf, double turb_ke, double levelset, CConfig *config) {
	unsigned short iDim;
  double epsilon, Heaviside, lambda, DensityInc, ViscosityInc;

	double ArtComp_Factor = config->GetArtComp_Factor();
  bool freesurface = config->GetFreeSurface();

  /*--- Set the value of the density and viscosity ---*/
	if (!freesurface) {
		SetDensityInc(Density_Inf);
		SetLaminarViscosityInc(Viscosity_Inf);
  }
  else {
    epsilon = config->GetFreeSurface_Thickness();
    Heaviside = 0.0;
    if (levelset < -epsilon) Heaviside = 1.0;
    if (fabs(levelset) <= epsilon) Heaviside = 1.0 - (0.5*(1.0+(levelset/epsilon)+(1.0/PI_NUMBER)*sin(PI_NUMBER*levelset/epsilon)));
    if (levelset > epsilon) Heaviside = 0.0;

    lambda = config->GetRatioDensity();
    DensityInc = (lambda + (1.0 - lambda)*Heaviside)*config->GetDensity_FreeStreamND();
		SetDensityInc(DensityInc);

    lambda = config->GetRatioViscosity();
    ViscosityInc = (lambda + (1.0 - lambda)*Heaviside)*config->GetViscosity_FreeStreamND();
    SetLaminarViscosityInc(ViscosityInc);
  }

  /*--- Set the value of the velocity squared (requires density) ---*/
	SetVelocityInc2();
  
  /*--- Set the value of the artificial compressibility factor ---*/
  SetBetaInc2(ArtComp_Factor);
  
  /*--- Set the value of the velocity ---*/
	for (iDim = 0; iDim < nDim; iDim++) 
		Primitive[iDim+1] = Solution[iDim+1] / GetDensityInc();

}

CTurbVariable::CTurbVariable(void) : CVariable() {

  /*--- Array initialization ---*/
	TS_Source = NULL;
  
}

CTurbVariable::CTurbVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {

  /*--- Array initialization ---*/
	TS_Source = NULL;

}

CTurbVariable::~CTurbVariable(void) {

	if (TS_Source != NULL) delete [] TS_Source;

}

double CTurbVariable::GetmuT(){ return muT; }

void CTurbVariable::SetmuT(double val_muT){ muT = val_muT; }

CTurbSAVariable::CTurbSAVariable(void) : CTurbVariable() { }

CTurbSAVariable::CTurbSAVariable(double val_nu_tilde, double val_muT, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
	/*--- Initialization of S-A variables ---*/
	Solution[0] = val_nu_tilde;		Solution_Old[0] = val_nu_tilde;

	/*--- Initialization of the eddy viscosity ---*/
	muT = val_muT;

	/*--- Allocate and initialize solution for the dual time strategy ---*/
	if (dual_time) {
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

CPlasmaVariable::CPlasmaVariable(void) : CVariable() {
	
  /*--- Array initialization ---*/
  Sensor_MultiSpecies = NULL;
	Velocity2 = NULL;
	Max_Lambda_Inv_MultiSpecies = NULL;
	Max_Lambda_Visc_MultiSpecies = NULL;
	Lambda = NULL;
	LaminarViscosity_MultiSpecies = NULL;
	EddyViscosity_MultiSpecies = NULL;
	ThermalCoeff = NULL;
  ThermalCoeff_vib = NULL;
	Species_Delta_Time = NULL;
	Residual_Chemistry = NULL;
	Residual_ElecForce = NULL;
	Residual_MomentumExch = NULL;
	Residual_EnergyExch = NULL;
	Elec_Field = NULL;
	B_Field = NULL;
	Primitive = NULL;
  Gradient_Primitive = NULL;
  LimiterPrimitive = NULL;

}

CPlasmaVariable::CPlasmaVariable(double *val_density, double **val_velocity, double *val_energy, double *val_energy_vib, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {

	unsigned short iVar, loc = 0, iSpecies, iDim, iMesh, nMGSmooth = 0, nPrimVar;
	double rho, energy,	energy_vib, energy_el, enthalpy_formation, Gamma;	
	bool viscous;
	if ((config->GetKind_Solver() == PLASMA_NAVIER_STOKES) || (config->GetKind_Solver() == ADJ_PLASMA_NAVIER_STOKES)) viscous = true;
  
  /*--- Array initialization ---*/
  Sensor_MultiSpecies = NULL;
	Velocity2 = NULL;
	Max_Lambda_Inv_MultiSpecies = NULL;
	Max_Lambda_Visc_MultiSpecies = NULL;
	Lambda = NULL;
	LaminarViscosity_MultiSpecies = NULL;
	EddyViscosity_MultiSpecies = NULL;
	ThermalCoeff = NULL;
  ThermalCoeff_vib = NULL;
	Species_Delta_Time = NULL;
	Residual_Chemistry = NULL;
	Residual_ElecForce = NULL;
	Residual_MomentumExch = NULL;
	Residual_EnergyExch = NULL;
	Elec_Field = NULL;
	B_Field = NULL;
	Primitive = NULL;
  Gradient_Primitive = NULL;
  LimiterPrimitive = NULL;

	nMonatomics  = config->GetnMonatomics();
	nDiatomics = config->GetnDiatomics();
	nSpecies = config->GetnSpecies();
	Prandtl_Lam     = config->GetPrandtl_Lam();
  nPrimVar = nDim+3;

	/*--- Allocate residual structures ---*/
	/*--- Allocate truncation error for multigrid strategy ---*/
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
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTERED) {
		Undivided_Laplacian = new double [nVar];
		Lambda = new double [nSpecies];
	}
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_UPWIND) Limiter = new double [nVar];
  LimiterPrimitive = new double *[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    LimiterPrimitive[iSpecies] = new double[nPrimVar];

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
	Velocity2   = new double [nSpecies];
	Max_Lambda_Inv_MultiSpecies = new double [nSpecies];
	Species_Delta_Time = new double [nSpecies];
	Sensor_MultiSpecies = new double [nSpecies];

	/*--- Compressible flow, primitive variables nDim+6, (T_tr,vx,vy,vz,T_vi, P, rho, h, c) ---*/
	Primitive = new double*[nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Primitive[iSpecies] = new double [nDim+6];
		for (iVar = 0; iVar < nDim+6; iVar++)
			Primitive[iSpecies][iVar] = 1.0;
	}

	/*--- Allocate and initialize the primitive gradient variable nDim+2, (T_tr,vx,vy,vz,T_vi) ---*/
	Gradient_Primitive = new double **[nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Gradient_Primitive[iSpecies] = new double *[nDim+3];
		for (iVar = 0; iVar < nDim+3; iVar++) {
			Gradient_Primitive[iSpecies][iVar] = new double [nDim];
			for (iDim = 0; iDim < nDim; iDim++)
				Gradient_Primitive[iSpecies][iVar][iDim] = 0.0;
		}
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

	Elec_Field  = new double [nDim];
	B_Field 	= new double [nDim];
	for(iDim = 0; iDim < nDim; iDim ++) {
		Elec_Field[iDim] = 0.0;
		B_Field[iDim] = 0.0;
	}

}
CPlasmaVariable::CPlasmaVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {

	unsigned short iVar, loc = 0, iSpecies, iDim, iMesh, nMGSmooth = 0, nPrimVar;
	double rho, energy, energy_vib, energy_el, Gamma, enthalpy_formation;
	bool viscous;
	if ((config->GetKind_Solver() == PLASMA_NAVIER_STOKES) || (config->GetKind_Solver() == ADJ_PLASMA_NAVIER_STOKES)) viscous = true;
  
  /*--- Array initialization ---*/
  Sensor_MultiSpecies = NULL;
	Velocity2 = NULL;
	Max_Lambda_Inv_MultiSpecies = NULL;
	Max_Lambda_Visc_MultiSpecies = NULL;
	Lambda = NULL;
	LaminarViscosity_MultiSpecies = NULL;
	EddyViscosity_MultiSpecies = NULL;
	ThermalCoeff = NULL;
  ThermalCoeff_vib = NULL;
	Species_Delta_Time = NULL;
	Residual_Chemistry = NULL;
	Residual_ElecForce = NULL;
	Residual_MomentumExch = NULL;
	Residual_EnergyExch = NULL;
	Elec_Field = NULL;
	B_Field = NULL;
	Primitive = NULL;
  Gradient_Primitive = NULL;
  LimiterPrimitive = NULL;
  
	nMonatomics  = config->GetnMonatomics();
	nDiatomics = config->GetnDiatomics();
	nSpecies = config->GetnSpecies();
	Prandtl_Lam     = config->GetPrandtl_Lam();
  nPrimVar = nDim+3;

	/*--- Allocate residual structures ---*/
	/*--- Allocate truncation error for multigrid strategy ---*/
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
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTERED) {
		Undivided_Laplacian = new double [nVar];
		Lambda = new double [nSpecies];
	}
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_UPWIND) Limiter = new double [nVar];
  LimiterPrimitive = new double *[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    LimiterPrimitive[iSpecies] = new double[nPrimVar];


	/*--- Solution initialization ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}


	Solution_Max = new double[nVar];
	Solution_Min = new double[nVar];

	/*--- Allocate and initialize the local useful quantities ---*/
	Velocity2   = new double [nSpecies];
	Max_Lambda_Inv_MultiSpecies = new double [nSpecies];
	Species_Delta_Time = new double [nSpecies];
	Sensor_MultiSpecies = new double [nSpecies];

	/*--- Compressible flow, primitive variables nDim+6, (T_tr,vx,vy,vz,T_vi, P, rho, h, c) ---*/
	Primitive = new double*[nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Primitive[iSpecies] = new double [nDim+6];
		for (iVar = 0; iVar < nDim+6; iVar++)
			Primitive[iSpecies][iVar] = 1.0;
	}

	/*--- Allocate and initialize the primitive gradient variable nDim+2, (T_tr,vx,vy,vz,T_vi) ---*/
	Gradient_Primitive = new double **[nSpecies];
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Gradient_Primitive[iSpecies] = new double *[nDim+3];
		for (iVar = 0; iVar < nDim+3; iVar++) {
			Gradient_Primitive[iSpecies][iVar] = new double [nDim];
			for (iDim = 0; iDim < nDim; iDim++)
				Gradient_Primitive[iSpecies][iVar][iDim] = 0.0;
		}
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
		Primitive[iSpecies][nDim+4]     = energy + Primitive[iSpecies][nDim+2]/rho;
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


	Elec_Field = new double [nDim];
	B_Field = new double [nDim];
	for(iDim =0; iDim < nDim; iDim ++) {
		Elec_Field[iDim] = 0.0;
		B_Field[iDim] = 0.0;
	}

}


CPlasmaVariable::~CPlasmaVariable(void) {
	unsigned short iVar, iSpecies;

	if (Elec_Field                    != NULL ) delete [] Elec_Field;
	if (LaminarViscosity_MultiSpecies != NULL ) delete [] LaminarViscosity_MultiSpecies;
	if (EddyViscosity_MultiSpecies    != NULL ) delete [] EddyViscosity_MultiSpecies;
	if (ThermalCoeff                  != NULL ) delete [] ThermalCoeff;
	if (ThermalCoeff_vib              != NULL ) delete [] ThermalCoeff_vib;
	if (Max_Lambda_Visc_MultiSpecies  != NULL ) delete [] Max_Lambda_Visc_MultiSpecies;
	if (Lambda                        != NULL ) delete [] Lambda;
  if (Max_Lambda_Inv_MultiSpecies   != NULL ) delete [] Max_Lambda_Inv_MultiSpecies;
  if (Sensor_MultiSpecies           != NULL ) delete [] Sensor_MultiSpecies;
  if (Velocity2                     != NULL ) delete [] Velocity2;
  if (Species_Delta_Time            != NULL ) delete [] Species_Delta_Time;
  if (B_Field                       != NULL ) delete [] B_Field;
  if (Residual_Chemistry            != NULL ) delete [] Residual_Chemistry;
  if (Residual_ElecForce            != NULL ) delete [] Residual_ElecForce;
  if (Residual_MomentumExch         != NULL ) delete [] Residual_MomentumExch;
  if (Residual_EnergyExch           != NULL ) delete [] Residual_EnergyExch;
  
  if (Primitive != NULL ) {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      delete Primitive[iSpecies];
    delete [] Primitive;
  }
  
  if (Gradient_Primitive != NULL ) {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (iVar = 0; iVar < nDim+3; iVar++)
        delete Gradient_Primitive[iSpecies][iVar];
      delete Gradient_Primitive[iSpecies];
    }
    delete [] Gradient_Primitive;
  }
  
  if (LimiterPrimitive != NULL ) {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      delete LimiterPrimitive[iSpecies];
    delete [] LimiterPrimitive;
  }
  
}

void CPlasmaVariable::SetVel_ResTruncError_Zero(unsigned short iSpecies) {
	unsigned short loc, iDim;
	if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
	else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);

	for (iDim = 0; iDim < nDim; iDim++)
		Res_TruncError[loc+iDim+1] = 0.0;
}

void CPlasmaVariable::SetGradient_PrimitiveZero(unsigned short val_primvar) {
	unsigned short iSpecies, iVar, iDim;

	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
		for (iVar = 0; iVar < val_primvar; iVar++)
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
		Primitive[iSpecies][nDim+5] = sqrt(fabs(Gamma*(Gamma - 1.0)*(Energy - 0.5*Velocity2[iSpecies] - Energy_vib - Energy_el - Enthalpy_formation)));
	}
}

bool CPlasmaVariable::SetPressure(CConfig *config) {

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
		Primitive[iSpecies][nDim+2] = (Gamma-1.0) * Density * (Energy - 0.5*Vel2 - Enthalpy_formation - Energy_vib - Energy_el);

	}
  if (Primitive[iSpecies][nDim+2] >= 0.0) return true;
  else return false;
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
	unsigned short iDim, iSpecies, loc, iVar;

	bool viscous = (config->GetKind_ViscNumScheme_Plasma());

	double Gamma, Enthalpy_formation, Gas_constant, Molar_mass, Char_temp_vib;
	double density, energy_tot, energy_vib, energy_el;

	/*--- Primitive variable vector: (Ttr, vx, vy, vz, Tvib, P, rho, h, c) ---*/
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
		else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
    
		/*--- Store useful quantities ---*/
		density    = Solution[loc+0];
		energy_tot = Solution[loc+nDim+1]/density;
    if ( iSpecies < nDiatomics ) energy_vib = Solution[loc+nDim+2]/Solution[loc+0];
		else energy_vib = 0.0;
		energy_el  = 0.0;

		Gamma = config->GetSpecies_Gamma(iSpecies);
		Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);
		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
		Molar_mass = config->GetMolar_Mass(iSpecies);

    /*--- Basic checkings to determine if the solution has a physical meaning ---*/
    bool check_density = false, check_sos = false, check_press = false,
    check_temp = false, check_vib_temp = false;
    
    if ( density < 0.0 ) check_density = true;

		double Vel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			Vel2 += Solution[loc+iDim+1]*Solution[loc+iDim+1]/(density*density);
    
    double sos = Gamma*(Gamma - 1.0)*(energy_tot - 0.5*Vel2 - energy_vib - energy_el - Enthalpy_formation);
    if ( sos < 0.0 ) check_sos = true;
    
    double press = (Gamma-1.0) * density * (energy_tot - 0.5*Vel2 - Enthalpy_formation - energy_vib - energy_el);
    if ( press < 0.0 ) check_press = true;
    
    double temp = (Gamma-1.0)/Gas_constant * (energy_tot - 0.5*Vel2 - Enthalpy_formation - energy_vib - energy_el);
    if ( temp < 0.0 ) check_temp = true;

    if (iSpecies < nDiatomics) {
			double temp_vib = config->GetCharVibTemp(iSpecies) / log(config->GetCharVibTemp(iSpecies)*UNIVERSAL_GAS_CONSTANT/(energy_vib*Molar_mass) + 1.0);
      if ( temp_vib < 0.0 ) check_vib_temp = true;
		}
    
    /*--- Copy old solution to the current solution because the update was not physical ---*/
    if (check_density || check_sos || check_press || check_temp || check_vib_temp) {
      if ( iSpecies < nDiatomics ) {
        for (iVar = 0; iVar < (nDim+3); iVar++)
          Solution[loc+iVar] = Solution_Old[loc+iVar];
      }
      else {
        for (iVar = 0; iVar < (nDim+2); iVar++)
          Solution[loc+iVar] = Solution_Old[loc+iVar];
      }
    }
    
		/*--- Store useful quantities ---*/
		density    = Solution[loc+0];
		energy_tot = Solution[loc+nDim+1]/density;
    if ( iSpecies < nDiatomics ) energy_vib = Solution[loc+nDim+2]/Solution[loc];
		else energy_vib = 0.0;
		energy_el  = 0.0;
    
    /*--- Velocity ---*/
		Velocity2[iSpecies] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity2[iSpecies] += Solution[loc+iDim+1]*Solution[loc+iDim+1]/(density*density);
			Primitive[iSpecies][iDim+1] = Solution[loc+iDim+1] / density;
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
		}
    else Primitive[iSpecies][nDim+1] = 0.0;

		/*--- Enthalpy ---*/
		Primitive[iSpecies][nDim+4] = energy_tot + Primitive[iSpecies][nDim+2]/density;
    
    
//		/*--- Error checking ---*/
//		if (Primitive[iSpecies][0] < 0) {
//			cout << "Temp TR[" << iSpecies << "] is negative! " << Primitive[iSpecies][0] << endl;
//			cout << "Coordinate: (" << Coord[0] << ", " << Coord[1] << ", "  << Coord[2] << ")" << endl;
//			Primitive[iSpecies][0] = EPS;
//		}
//		if (Primitive[iSpecies][nDim+1] < 0) {
//			cout << "Temp vib[" << iSpecies << "] is is negative! " << Primitive[iSpecies][nDim+1] << endl;
//			Primitive[iSpecies][nDim+1] = EPS;
//		}
//		if (Primitive[iSpecies][nDim+2] < 0) {
//			cout << "Pressure[" << iSpecies << "] is negative! " << Primitive[iSpecies][nDim+2] << endl;
//			cout << "Coordinate: (" << Coord[0] << ", " << Coord[1] << ", "  << Coord[2] << ")" << endl;
//			Primitive[iSpecies][nDim+2] = EPS;
//		}
//		if (Primitive[iSpecies][nDim+3] < 0) {
//			cout << "Density[" << iSpecies << "] is negative! " << Primitive[iSpecies][nDim+3] << endl;
//			cout << "Coordinate: (" << Coord[0] << ", " << Coord[1] << ", "  << Coord[2] << ")" << endl;
//			Primitive[iSpecies][nDim+3] = EPS;
//		}
    
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
		Primitive[iSpecies][nDim+5] = sqrt(fabs(Gamma*(Gamma - 1.0)*(energy_tot - 0.5*Velocity2[iSpecies] - energy_vib - energy_el - Enthalpy_formation)));

		/*--- Pressure ---*/
		Primitive[iSpecies][nDim+2] = (Gamma-1.0) * density * (energy_tot - 0.5*Velocity2[iSpecies] - Enthalpy_formation - energy_vib - energy_el);

		/*--- Translational-rotational temperature ---*/
		Primitive[iSpecies][0] = (Gamma-1.0)/Gas_constant * (energy_tot - 0.5*Velocity2[iSpecies] - Enthalpy_formation - energy_vib - energy_el);

    if (Primitive[iSpecies][0] < 0.0) Primitive[iSpecies][0] = 0;

		/*--- Vibrational temperature ---*/
		if (iSpecies < nDiatomics) {
			Char_temp_vib = config->GetCharVibTemp(iSpecies);
			Primitive[iSpecies][nDim+1] = Char_temp_vib / log(Char_temp_vib*UNIVERSAL_GAS_CONSTANT/(energy_vib*Molar_mass) + 1.0);

      if (Primitive[iSpecies][nDim+1] < 0.0) Primitive[iSpecies][nDim+1] = 0;

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
	unsigned short iSpecies, jSpecies, iLoc, jLoc;
	double Temperature_Ref, Viscosity_Ref;
  double densityMixture, molarConc_i, molarConc_j, Mi, Mj, Tc;
  double delta2_ij, collisionArea, denom;
  double ***Omega11;

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
		LaminarViscosity_MultiSpecies[iSpecies] = 1E-6*2.3E-5*1.3E-7*(pow(Temperature_Dim/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_Dim+110.3));
		LaminarViscosity_MultiSpecies[iSpecies] = LaminarViscosity_MultiSpecies[iSpecies]/Viscosity_Ref;
		break;


	case N2:
            
      /*--- Gupta-Yos approximation ---*/
      densityMixture = 0.0;
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
        else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
        densityMixture += Solution[iLoc];
      }
      
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
        else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
        
        denom = 0.0;
        
        /*--- Calculate molar concentration ---*/
        Mi = config->GetMolar_Mass(iSpecies);
        molarConc_i = Solution[iLoc] / (densityMixture * Mi);
        
        for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
          if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
          else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);
          Mj = config->GetMolar_Mass(jSpecies);
          molarConc_j = Solution[jLoc] / (densityMixture * Mj);
          
          /*--- Calculate "delta" quantities ---*/
          Tc = sqrt(Primitive[iSpecies][0]*Primitive[jSpecies][0]);
          Omega11 = config->GetCollisionIntegral11();
          collisionArea = 1E-20 * Omega11[iSpecies][jSpecies][3] * pow(Tc, Omega11[iSpecies][jSpecies][0]*log(Tc)*log(Tc) + Omega11[iSpecies][jSpecies][1]*log(Tc) + Omega11[iSpecies][jSpecies][2]);
          delta2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (PI_NUMBER * UNIVERSAL_GAS_CONSTANT * Tc * (Mi+Mj))) * collisionArea;

          /*--- Add to denominator of viscosity ---*/
          denom += molarConc_j*delta2_ij;
        }
        
        /*--- Calculate species laminar viscosity ---*/
        LaminarViscosity_MultiSpecies[iSpecies] = (Mi/AVOGAD_CONSTANT * molarConc_i) / denom;
      }
		break;
	}
}


void CPlasmaVariable ::SetThermalCoeff(CConfig *config) {
	unsigned short iSpecies, jSpecies, iLoc, jLoc;
	double Gamma, Gas_constant, cv_vib;
  double densityMixture, molarConc_i, molarConc_j, collisionArea, Mi, Mj, Tc, Tv, Theta_v;
  double denom_t, denom_r, delta1_ij, delta2_ij, a_ij;
  double ***Omega00, ***Omega11;

	if (config->GetKind_GasModel() == ARGON) {
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			Gamma = config->GetSpecies_Gamma(iSpecies);
			Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);
			ThermalCoeff[iSpecies] = Gamma * LaminarViscosity_MultiSpecies[iSpecies] *Gas_constant/ ( (Gamma-1.0)*Prandtl_Lam );
		}
	} else {
    
    /*--- Gupta-Yos approximation ---*/
    densityMixture = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
      else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
      densityMixture += Solution[iLoc];
    }
    
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      if ( iSpecies < nDiatomics ) iLoc = (nDim+3)*iSpecies;
      else iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
      
      /*--- Calculate molar concentration ---*/
      Mi = config->GetMolar_Mass(iSpecies);
      molarConc_i = Solution[iLoc] / (densityMixture * Mi);
      Theta_v = config->GetCharVibTemp(iSpecies);
      Tv = Primitive[iSpecies][nDim+1];
      cv_vib = 0.0;
      if (iSpecies < nDiatomics)
        cv_vib = UNIVERSAL_GAS_CONSTANT/Mi * (Theta_v/Tv)*(Theta_v/Tv) * exp(Theta_v/Tv) / ((exp(Theta_v/Tv) - 1.0)*(exp(Theta_v/Tv) - 1.0));
      
      denom_t = 0.0;
      denom_r = 0.0;
      
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        if ( jSpecies < nDiatomics ) jLoc = (nDim+3)*jSpecies;
        else jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);
        Mj = config->GetMolar_Mass(jSpecies);
        molarConc_j = Solution[jLoc] / (densityMixture * Mj);
        
        /*--- Calculate "delta" quantities ---*/
        Tc = sqrt(Primitive[iSpecies][0]*Primitive[jSpecies][0]);
        Omega00 = config->GetCollisionIntegral00();
        Omega11 = config->GetCollisionIntegral11();
        a_ij = 1.0 + (1.0 - Mi/Mj)*(0.45 - 2.54*Mi/Mj) / ((1.0 + Mi/Mj)*(1.0 + Mi/Mj));
        
        collisionArea = 1E-20 * Omega00[iSpecies][jSpecies][3] * pow(Tc, Omega00[iSpecies][jSpecies][0]*log(Tc)*log(Tc) + Omega00[iSpecies][jSpecies][1]*log(Tc) + Omega00[iSpecies][jSpecies][2]);
        delta1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (PI_NUMBER * UNIVERSAL_GAS_CONSTANT * Tc * (Mi+Mj))) * collisionArea;
        
        collisionArea = 1E-20 * Omega11[iSpecies][jSpecies][3] * pow(Tc, Omega11[iSpecies][jSpecies][0]*log(Tc)*log(Tc) + Omega11[iSpecies][jSpecies][1]*log(Tc) + Omega11[iSpecies][jSpecies][2]);
        delta2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (PI_NUMBER * UNIVERSAL_GAS_CONSTANT * Tc * (Mi+Mj))) * collisionArea;
        
        
        /*--- Add to denominator of translational & rotational thermal conductivity ---*/
        denom_t += a_ij*molarConc_j*delta2_ij;
        denom_r += molarConc_j*delta1_ij;
      }
            
      /*--- Calculate species laminar viscosity ---*/
      ThermalCoeff[iSpecies] = 15.0/4.0 * BOLTZMANN_CONSTANT * molarConc_i / denom_t + BOLTZMANN_CONSTANT * molarConc_i / denom_r;
      ThermalCoeff_vib[iSpecies] = BOLTZMANN_CONSTANT * cv_vib/config->GetSpecies_Gas_Constant(iSpecies) * molarConc_i / denom_r;
    }
	}
}

CLevelSetVariable::CLevelSetVariable(void) : CVariable() {}

CLevelSetVariable::CLevelSetVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {

	/*--- Allocate residual structures ---*/
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];

	/*--- Allocate limiter (upwind) ---*/
	Limiter = new double [nVar];
	Solution_Max = new double [nVar];
	Solution_Min = new double [nVar];
  
  /*--- Allocate primitive solution ---*/
	Primitive = new double [nVar];

}

CLevelSetVariable::CLevelSetVariable(double val_levelset, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar,config) {

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
	/*--- Allocate residual structures ---*/
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];

	/*--- Allocate limiter (upwind)---*/
	Limiter = new double [nVar];
	Solution_Max = new double [nVar];
	Solution_Min = new double [nVar];
  
  /*--- Allocate primitive solution ---*/
	Primitive = new double [nVar];

	/*--- Solution and old solution initialization ---*/
	Solution[0] = val_levelset;		Solution_Old[0] = val_levelset;

	if (dual_time) {
		Solution_time_n[0] = val_levelset;
		Solution_time_n1[0] = val_levelset;
	}
  
  /*--- Set primitive solution ---*/
  Primitive[0] = val_levelset;

}

CLevelSetVariable::~CLevelSetVariable(void) { }

CWaveVariable::CWaveVariable(void) : CVariable() {

  /*--- Array initialization ---*/
	Solution_Direct = NULL;
  
}

CWaveVariable::CWaveVariable(double *val_wave, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

  /*--- Array initialization ---*/
	Solution_Direct = NULL;
  
	/*--- Allocate residual structures ---*/
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

	if (Solution_Direct != NULL) delete [] Solution_Direct;

}

CFEAVariable::CFEAVariable(void) : CVariable() { }

CFEAVariable::CFEAVariable(double *val_fea, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

	/*--- Allocate residual structures ---*/
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];

	/*--- Initialization of variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_fea[iVar];
		Solution_Old[iVar] = val_fea[iVar];
	}
}

CFEAVariable::~CFEAVariable(void) { }

CHeatVariable::CHeatVariable(void) : CVariable() {

  /*--- Array initialization ---*/
	Solution_Direct = NULL;
  
}

CHeatVariable::CHeatVariable(double *val_heat, unsigned short val_ndim, unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar;

  /*--- Array initialization ---*/
	Solution_Direct = NULL;
  
	/*--- Allocate residual structures ---*/
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

  if (Solution_Direct != NULL) delete [] Solution_Direct;

}

CTNE2EulerVariable::CTNE2EulerVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	Primitive = NULL;
	Gradient_Primitive = NULL;
	Limiter_Primitive = NULL;
  
}

CTNE2EulerVariable::CTNE2EulerVariable(double val_density, double *val_massfrac, double *val_velocity, double val_temperature, double val_temperature_ve, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim, iSpecies, iMesh, nMGSmooth = 0;
  unsigned short nDim;
  double energy, energy_v, energy_e, energy_ve, energy_formation, sqvel;
  double *molar_mass, *theta_v, *rotation_modes, *temperature_ref, *enthalpy_formation;
  
  ionization = ((config->GetKind_GasModel() == AIR7) || (config->GetKind_GasModel() == ARGON));
  
  /*--- Array initialization ---*/
	Primitive = NULL;
	Gradient_Primitive = NULL;
	Limiter_Primitive = NULL;
  
  /*--- Acquire parameters from the config class ---*/
  nDim               = val_ndim;
  nSpecies           = config->GetnSpecies();
  molar_mass         = config->GetMolar_Mass();
  theta_v            = config->GetCharVibTemp();
  rotation_modes     = config->GetRotationModes();
  temperature_ref    = config->GetRefTemperature();
  enthalpy_formation = config->GetEnthalpy_Formation();
  
  
  /*--- Allocate and initialize the primitive variables and gradients ---*/
  nPrimVar = nSpecies+nDim+6; nPrimVarGrad = nSpecies+nDim+4;
  
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
  
	/*--- Allocate limiter (upwind)---*/
	if ((config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) &&
			(config->GetKind_SlopeLimit_Flow() != NONE)) {
		Limiter      = new double [nVar];
		Solution_Max = new double [nVar];
		Solution_Min = new double [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Limiter[iVar]      = 0.0;
			Solution_Max[iVar] = 0.0;
			Solution_Min[iVar] = 0.0;
		}
	}
  
  /*--- Calculate energy ---*/
  sqvel     = 0.0;
  energy    = 0.0;
  energy_ve = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    sqvel += val_velocity[iDim]*val_velocity[iDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    energy_v  =  UNIVERSAL_GAS_CONSTANT/molar_mass[iSpecies]
               * theta_v[iSpecies] / ( exp(theta_v[iSpecies]/val_temperature_ve) -1.0 );
    energy_e = 0.0;
    energy_formation = enthalpy_formation[iSpecies] - UNIVERSAL_GAS_CONSTANT/molar_mass[iSpecies]*temperature_ref[iSpecies];
    energy +=  (3.0/2.0 + rotation_modes[iSpecies]/2.0) * (val_temperature - temperature_ref[iSpecies])
             + energy_v + energy_e + energy_formation + 0.5*sqvel;
    energy_ve += energy_v + energy_e;
  }
  if (ionization) {
    iSpecies = nSpecies - 1;
    energy_formation = enthalpy_formation[iSpecies] - UNIVERSAL_GAS_CONSTANT/molar_mass[iSpecies]*temperature_ref[iSpecies];
    energy          -= (3.0/2.0) * (val_temperature-temperature_ref[iSpecies]) + energy_formation + 0.5*sqvel;
    energy_ve       += (3.0/2.0) * (val_temperature-temperature_ref[iSpecies]) + energy_formation + 0.5*sqvel;
  }
  
	/*--- Solution and old solution initialization ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Solution[iSpecies] = val_density*val_massfrac[iSpecies];
    Solution_Old[0]    = val_density*val_massfrac[iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    Solution[nSpecies+iDim+1]     = val_density*val_velocity[iDim];
    Solution_Old[nSpecies+iDim+1] = val_density*val_velocity[iDim];
  }
  Solution[nSpecies+nDim+1]     = val_density*energy;
  Solution_Old[nSpecies+nDim+1] = val_density*energy;
  Solution[nSpecies+nDim+2]     = val_density*energy_ve;
  Solution_Old[nSpecies+nDim+2] = val_density*energy_ve;
  
  /*--- Primitive variables: nSpe+nDim+6 (rho1,...,rhoNs,T,Tve,vx,vy,vz,P,rho,h,c) ---*/
  Primitive = new double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;
  
  /*--- Primitive variable gradients: nSpe+nDim+4, (rho1,rhoNs,T,Tve,vx,vy,vz,P,rho)
   We need P, and rho for running the adjoint problem ---*/
  Gradient_Primitive = new double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }
}

CTNE2EulerVariable::CTNE2EulerVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config) {
	unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Array initialization ---*/
	Primitive = NULL;
	Gradient_Primitive = NULL;
  Limiter_Primitive = NULL;
  
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
  
	/*--- Allocate limiter (upwind)---*/
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
  
	/*--- Allocate and initialize the primitive variables and gradients ---*/
  nSpecies = config->GetnSpecies();
  nPrimVar = nSpecies+nDim+6; nPrimVarGrad = nSpecies+nDim+4;
  
  /*--- Primitive variables nSpe+nDim+6, (rho1,...,rhoNs,T,Tve,vx,vy,vz,P,rho,h,c) ---*/
  Primitive = new double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;
  
  /*--- Primitive variable gradient nSpe+nDim+4, (rho1,...,rhoNs,T,Tve,vx,vy,vz,P,rho)
   We need P, and rho for running the adjoint problem ---*/
  Gradient_Primitive = new double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }
}

CTNE2EulerVariable::~CTNE2EulerVariable(void) {
	unsigned short iVar;
  
  if (Primitive         != NULL) delete [] Primitive;
  if (Limiter_Primitive != NULL) delete [] Limiter_Primitive;
  
  if (Res_TruncError != NULL) delete [] Res_TruncError;
  if (Residual_Old   != NULL) delete [] Residual_Old;
  if (Residual_Sum   != NULL) delete [] Residual_Sum;
  if (Limiter != NULL) delete [] Limiter;
  if (Solution_Max != NULL) delete [] Solution_Max;
  if (Solution_Min != NULL) delete [] Solution_Min;
  
  if (Primitive != NULL) delete [] Primitive;
  if (Gradient_Primitive != NULL) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      delete Gradient_Primitive[iVar];
    delete [] Gradient_Primitive;
  }
  
}

void CTNE2EulerVariable::SetGradient_PrimitiveZero(unsigned short val_primvar) {
	unsigned short iVar, iDim;
  
	for (iVar = 0; iVar < val_primvar; iVar++)
		for (iDim = 0; iDim < nDim; iDim++)
			Gradient_Primitive[iVar][iDim] = 0.0;
}

double CTNE2EulerVariable::GetProjVel(double *val_vector) {
	double ProjVel, density;
	unsigned short iDim, iSpecies;
  
	ProjVel = 0.0;
  density = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    density += Solution[iSpecies];
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVel += Solution[nSpecies+iDim+1]*val_vector[iDim]/density;
  
	return ProjVel;
}

inline bool CTNE2EulerVariable::SetTemperature(CConfig *config) {
  unsigned short iSpecies, iDim;
  double rho, rhoe, rhoev, rhoeform, rhoeref, rhoCvtr, sqvel;
  double *hf, *Tref, *rotmodes, *molarmass;
  
  /*--- Load variables from the config class --*/
  hf        = config->GetEnthalpy_Formation();
  Tref      = config->GetRefTemperature();
  rotmodes  = config->GetRotationModes();
  molarmass = config->GetMolar_Mass();
  
  rhoe     = Solution[nSpecies+nDim+1];
  rhoev    = Solution[nSpecies+nDim+2];
  rho      = 0.0;
  rhoeform = 0.0;
  rhoeref  = 0.0;
  rhoCvtr  = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    rho += Solution[iSpecies];
    rhoeform += Solution[iSpecies] * (hf[iSpecies] - UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies]);
    rhoeref  +=  Solution[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
            * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
    rhoCvtr +=  Solution[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
              * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
  }
  if (ionization) {
    iSpecies = nSpecies-1;
    rhoCvtr    -=  Solution[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
                 * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    rhoeref    -=  Solution[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
                 * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies]*Tref[iSpecies];
  }
  sqvel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    sqvel += (Solution[nSpecies+iDim+1]/rho) * (Solution[nSpecies+iDim+1]/rho);
  
  /*--- Assign translational-rotational temperature ---*/
  Primitive[nSpecies+1] = (rhoe - rhoev - rhoeform + rhoeref - 0.5*rho*sqvel) / rhoCvtr;
  
  if (Primitive[nSpecies+1] > 0.0) return false;
  else return true;
}

bool CTNE2EulerVariable::SetPressure(CConfig *config) {
  unsigned short iSpecies;
  double *molarmass;
  double P;
  
  molarmass = config->GetMolar_Mass();
  P = 0.0;
  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    P += Solution[iSpecies]*UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Primitive[nSpecies+1];
  }
  if (ionization) {
    iSpecies = nSpecies - 1;
    P -= Solution[iSpecies]*UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Primitive[nSpecies+1];
    P += Solution[iSpecies]*UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Primitive[nSpecies+2];
  }
  
  Primitive[nSpecies+nDim+2] = P;
  
  if (Primitive[nSpecies+nDim+2] > 0.0) return false;
  else return true;
}

void CTNE2EulerVariable::SetPrimVar_Compressible(double Gamma, double Gas_Constant) {
	unsigned short iDim, iVar, iSpecies;
  bool check_dens = false, check_press = false, check_sos = false, check_temp = false;
  
	SetVelocity2();                               // Compute the modulus of the velocity.
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    check_dens = ((Solution[iSpecies] < 0.0) || check_dens);  // Check the density
  check_temp  = SetTemperature(config);
	check_press = SetPressure(Temperature, Temperature_ve);	// Requires T & Tve computation.
	check_sos = SetSoundSpeed(Gamma);             // Requires pressure computation.
  
  /*--- Check that the solution has a physical meaning ---*/
  if (check_dens || check_press || check_sos || check_temp) {
    
    /*--- Copy the old solution ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    
    /*--- Recompute the primitive variables ---*/
    SetVelocity2();
    check_press = SetPressure(Gamma);
    check_sos = SetSoundSpeed(Gamma);
    check_temp = SetTemperature(Gas_Constant);
    
  }
  // Primitive variables [rho1,...,rhoNs,T,Tve,u,v,w,P,rho,h,c]
  SetEnthalpy();                                // Requires pressure computation.
  
	for (iDim = 0; iDim < nDim; iDim++)
		Primitive[iDim+1] = Solution[iDim+1] / Solution[0];
	Primitive[nDim+2] = Solution[0];
  
}

CTNE2NSVariable::CTNE2NSVariable(void) : CTNE2EulerVariable() { }

CTNE2NSVariable::CTNE2NSVariable(double val_density, double *val_massfrac, double *val_velocity, double val_temperature, double val_energy_ve, unsigned short val_ndim, unsigned short val_nvar, CConfig *config) : CTNE2EulerVariable(val_density, val_massfrac, val_velocity, val_temperature, val_energy_ve, val_ndim, val_nvar, config) {
  
	Temperature_Ref = config->GetTemperature_Ref();
	Viscosity_Ref   = config->GetViscosity_Ref();
	Viscosity_Inf   = config->GetViscosity_FreeStreamND();
	Prandtl_Lam     = config->GetPrandtl_Lam();
  
}

CTNE2NSVariable::CTNE2NSVariable(double *val_solution, unsigned short val_ndim,
                         unsigned short val_nvar, CConfig *config) : CTNE2EulerVariable(val_solution, val_ndim, val_nvar, config) {
  
	Temperature_Ref = config->GetTemperature_Ref();
	Viscosity_Ref   = config->GetViscosity_Ref();
	Viscosity_Inf   = config->GetViscosity_FreeStreamND();
	Prandtl_Lam     = config->GetPrandtl_Lam();

}

CTNE2NSVariable::~CTNE2NSVariable(void) { }

void CTNE2NSVariable::SetLaminarViscosity() {
	double Temperature_Dim;
  
	/*--- Calculate viscosity from a non-dim. Sutherland's Law ---*/
	Temperature_Dim = Primitive[0]*Temperature_Ref;
	LaminarViscosity = 1.853E-5*(pow(Temperature_Dim/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_Dim+110.3));
	LaminarViscosity = LaminarViscosity/Viscosity_Ref;
  
}

void CTNE2NSVariable::SetVorticity(void) {
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

void CTNE2NSVariable::SetPrimVar_Compressible(double Gamma, double Gas_Constant, double turb_ke) {
	unsigned short iDim, iVar;
  bool check_dens = false, check_press = false, check_sos = false, check_temp = false;
  
	SetVelocity2();                                 // Compute the modulus of the velocity.
  check_dens = (Solution[0] < 0.0);               // Check the density
	check_press = SetPressure(Gamma, turb_ke);      // Requires Velocity2 computation.
	check_sos = SetSoundSpeed(Gamma);               // Requires pressure computation.
	check_temp = SetTemperature(Gas_Constant);      // Requires pressure computation.
  
  /*--- Check that the solution has a physical meaning ---*/
  if (check_dens || check_press || check_sos || check_temp) {
    
    /*--- Copy the old solution ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    
    /*--- Recompute the primitive variables ---*/
    SetVelocity2();
    check_press = SetPressure(Gamma, turb_ke);
    check_sos = SetSoundSpeed(Gamma);
    check_temp = SetTemperature(Gas_Constant);
    
  }
  
	SetEnthalpy();                                  // Requires pressure computation.
	SetLaminarViscosity();                          // Requires temperature computation.
  
	for (iDim = 0; iDim < nDim; iDim++)
		Primitive[iDim+1] = Solution[iDim+1] / Solution[0];
	Primitive[nDim+2] = Solution[0];
  
}
