/*!
 * \file variable_direct_mean.cpp
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
    
    /*--- Output to control the convergence ---*/
    cout << "Negative density, pressure, sos, or temperature. Using old solution." << endl;
    
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

void CNSVariable::SetEddyViscosity(unsigned short val_Kind_Turb_Model, CVariable *TurbVariable) {
  
	switch (val_Kind_Turb_Model) {
    case NONE :
      EddyViscosity = 0.0;                     break;
    case SA : case SST :
      EddyViscosity = TurbVariable->GetmuT(); break;
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
    
    /*--- Output to control the convergence ---*/
    cout << "Negative density, pressure, sos, or temperature. Using old solution." << endl;
    
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