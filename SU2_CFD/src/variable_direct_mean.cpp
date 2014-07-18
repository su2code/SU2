/*!
 * \file variable_direct_mean.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.0 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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

CEulerVariable::CEulerVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	TS_Source = NULL;
	Primitive = NULL;
	Gradient_Primitive = NULL;
	Limiter_Primitive = NULL;
  WindGust = NULL;
  WindGustDer = NULL;
  
}

CEulerVariable::CEulerVariable(double val_density, double *val_velocity, double val_energy, unsigned short val_nDim,
                               unsigned short val_nvar, CConfig *config) : CVariable(val_nDim, val_nvar, config) {
	unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool low_fidelity = config->GetLowFidelitySim();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool viscous = config->GetViscous();
  bool windgust = config->GetWind_Gust();
  
  /*--- Array initialization ---*/
	TS_Source = NULL;
	Primitive = NULL;
	Gradient_Primitive = NULL;
	Limiter_Primitive = NULL;
  WindGust = NULL;
  WindGustDer = NULL;

  /*--- Allocate and initialize the primitive variables and gradients ---*/
  if (incompressible) { nPrimVar = nDim+5; nPrimVarGrad = nDim+3; }
  if (freesurface)    { nPrimVar = nDim+7; nPrimVarGrad = nDim+6; }
  if (compressible)   { nPrimVar = nDim+7; nPrimVarGrad = nDim+4;
    if (viscous) { nSecondaryVar = 8; nSecondaryVarGrad = 2; }
    else { nSecondaryVar = 2; nSecondaryVarGrad = 2; }
  }

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
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
		Undivided_Laplacian = new double [nVar];
  }
  
  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/
  Limiter_Primitive = new double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Limiter_Primitive[iVar] = 0.0;
  
  Limiter_Secondary = new double [nSecondaryVarGrad];
  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
    Limiter_Secondary[iVar] = 0.0;
  
  Limiter = new double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;
  
  Solution_Max = new double [nPrimVarGrad];
  Solution_Min = new double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }
  
	/*--- Solution and old solution initialization ---*/
	if (compressible) {
		Solution[0] = val_density;
		Solution_Old[0] = val_density;
		for (iDim = 0; iDim < nDim; iDim++) {
			Solution[iDim+1] = val_density*val_velocity[iDim];
			Solution_Old[iDim+1] = val_density*val_velocity[iDim];
		}
		Solution[nVar-1] = val_density*val_energy;
		Solution_Old[nVar-1] = val_density*val_energy;
	}
	if (incompressible || freesurface) {
		Solution[0] = config->GetPressure_FreeStreamND();
		Solution_Old[0] = config->GetPressure_FreeStreamND();
		for (iDim = 0; iDim < nDim; iDim++) {
			Solution[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
			Solution_Old[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
		}
	}
  
	/*--- Allocate and initialize solution for dual time strategy ---*/
	if (dual_time) {
    if (compressible) {
			Solution_time_n[0] = val_density;
			Solution_time_n1[0] = val_density;
			for (iDim = 0; iDim < nDim; iDim++) {
				Solution_time_n[iDim+1] = val_density*val_velocity[iDim];
				Solution_time_n1[iDim+1] = val_density*val_velocity[iDim];
			}
			Solution_time_n[nVar-1] = val_density*val_energy;
			Solution_time_n1[nVar-1] = val_density*val_energy;
		}
    if (incompressible || freesurface) {
			Solution_time_n[0] = config->GetPressure_FreeStreamND();
			Solution_time_n1[0] = config->GetPressure_FreeStreamND();
			for (iDim = 0; iDim < nDim; iDim++) {
				Solution_time_n[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
				Solution_time_n1[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
			}
		}
	}
  
	/*--- Allocate space for the time spectral source terms ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		TS_Source = new double[nVar];
		for (iVar = 0; iVar < nVar; iVar++) TS_Source[iVar] = 0.0;
	}
    
  /*--- Allocate vector for wind gust and wind gust derivative field ---*/
	if (windgust) {
    WindGust = new double [nDim];
    WindGustDer = new double [nDim+1];
  }
  
	/*--- Allocate auxiliar vector for free surface source term ---*/
	if (freesurface) Grad_AuxVar = new double [nDim];
  
  /*--- Incompressible flow, primitive variables nDim+3, (P,vx,vy,vz,rho,beta),
        FreeSurface Incompressible flow, primitive variables nDim+4, (P,vx,vy,vz,rho,beta,dist),
        Compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c) ---*/
  Primitive = new double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;
  
  Secondary = new double [nSecondaryVar];
  for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar] = 0.0;

  /*--- Incompressible flow, gradients primitive variables nDim+2, (P,vx,vy,vz,rho),
        FreeSurface Incompressible flow, primitive variables nDim+3, (P,vx,vy,vz,rho,beta,dist),
        Compressible flow, gradients primitive variables nDim+4, (T,vx,vy,vz,P,rho,h)
        We need P, and rho for running the adjoint problem ---*/
  Gradient_Primitive = new double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }
  
  Gradient_Secondary = new double* [nSecondaryVarGrad];
  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
    Gradient_Secondary[iVar] = new double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Secondary[iVar][iDim] = 0.0;
  }
  
}

CEulerVariable::CEulerVariable(double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config) : CVariable(val_nDim, val_nvar, config) {
	unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool low_fidelity = config->GetLowFidelitySim();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool viscous = config->GetViscous();
  bool windgust = config->GetWind_Gust();
  
  /*--- Array initialization ---*/
	TS_Source = NULL;
	Primitive = NULL;
	Gradient_Primitive = NULL;
  Limiter_Primitive = NULL;
  WindGust = NULL;
  WindGustDer = NULL;
  
	/*--- Allocate and initialize the primitive variables and gradients ---*/
  if (incompressible) { nPrimVar = nDim+5; nPrimVarGrad = nDim+3; }
  if (freesurface)    { nPrimVar = nDim+7; nPrimVarGrad = nDim+6; }
  if (compressible)   { nPrimVar = nDim+7; nPrimVarGrad = nDim+4;
    if (viscous) { nSecondaryVar = 8; nSecondaryVarGrad = 2; }
    else { nSecondaryVar = 2; nSecondaryVarGrad = 2; }
  }
  
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
  
  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/
  Limiter_Primitive = new double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Limiter_Primitive[iVar] = 0.0;
  
  Limiter_Secondary = new double [nSecondaryVarGrad];
  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
    Limiter_Secondary[iVar] = 0.0;

  Limiter = new double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;
  
  Solution_Max = new double [nPrimVarGrad];
  Solution_Min = new double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
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
		for (iVar = 0; iVar < nVar; iVar++) TS_Source[iVar] = 0.0;
	}
    
  /*--- Allocate vector for wind gust and wind gust derivative field ---*/
	if (windgust) {
    WindGust = new double [nDim];
    WindGustDer = new double [nDim+1];
  }
  
	/*--- Allocate auxiliar vector for free surface source term ---*/
	if (freesurface) Grad_AuxVar = new double [nDim];

  /*--- Incompressible flow, primitive variables nDim+3, (P,vx,vy,vz,rho,beta),
        FreeSurface Incompressible flow, primitive variables nDim+4, (P,vx,vy,vz,rho,beta,dist),
        Compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c) ---*/
  Primitive = new double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;
  
  Secondary = new double [nSecondaryVar];
  for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar] = 0.0;

  /*--- Incompressible flow, gradients primitive variables nDim+2, (P,vx,vy,vz,rho),
        FreeSurface Incompressible flow, primitive variables nDim+4, (P,vx,vy,vz,rho,beta,dist),
        Compressible flow, gradients primitive variables nDim+4, (T,vx,vy,vz,P,rho,h)
        We need P, and rho for running the adjoint problem ---*/
  Gradient_Primitive = new double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }
  
  Gradient_Secondary = new double* [nSecondaryVarGrad];
  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
    Gradient_Secondary[iVar] = new double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Secondary[iVar][iDim] = 0.0;
  }
  
}

CEulerVariable::~CEulerVariable(void) {
	unsigned short iVar;
  
	if (TS_Source         != NULL) delete [] TS_Source;
  if (Primitive         != NULL) delete [] Primitive;
  if (Limiter_Primitive != NULL) delete [] Limiter_Primitive;
  if (WindGust          != NULL) delete [] WindGust;
  if (WindGustDer       != NULL) delete [] WindGustDer;

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

void CEulerVariable::SetGradient_SecondaryZero(unsigned short val_secondaryvar) {
	unsigned short iVar, iDim;
  
	for (iVar = 0; iVar < val_secondaryvar; iVar++)
		for (iDim = 0; iDim < nDim; iDim++)
			Gradient_Secondary[iVar][iDim] = 0.0;
}

double CEulerVariable::GetProjVel(double *val_vector) {
	double ProjVel;
	unsigned short iDim;
  
	ProjVel = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVel += Primitive[iDim+1]*val_vector[iDim];
  
	return ProjVel;
}

bool CEulerVariable::SetPrimVar_Compressible(CFluidModel *FluidModel) {
	unsigned short iVar;
  bool check_dens = false, check_press = false, check_sos = false, check_temp = false, RightVol = true;
  

  SetVelocity();   // Computes velocity and velocity^2
  double density = GetDensity();
  double staticEnergy = GetEnergy()-0.5*Velocity2;
  /* check will be moved inside fluid model plus error description strings*/
  FluidModel->SetTDState_rhoe(density, staticEnergy);
  
  check_dens = SetDensity();
  check_press = SetPressure(FluidModel->GetPressure());
  check_sos = SetSoundSpeed(FluidModel->GetSoundSpeed2());
  check_temp = SetTemperature(FluidModel->GetTemperature());
  
  /*--- Check that the solution has a physical meaning ---*/
  
  if (check_dens || check_press || check_sos || check_temp) {
    
    /*--- Copy the old solution ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    
    /*--- Recompute the primitive variables ---*/
    
    SetVelocity();   // Computes velocity and velocity^2
    double density = GetDensity();
    double staticEnergy = GetEnergy()-0.5*Velocity2;
    /* check will be moved inside fluid model plus error description strings*/
    FluidModel->SetTDState_rhoe(density, staticEnergy);

    check_dens = SetDensity();
    check_press = SetPressure(FluidModel->GetPressure());
    check_sos = SetSoundSpeed(FluidModel->GetSoundSpeed2());
    check_temp = SetTemperature(FluidModel->GetTemperature());
    
    RightVol = false;
    
  }
  
  /*--- Set enthalpy ---*/
  
  SetEnthalpy();                                // Requires pressure computation.
  
  return RightVol;
  
}

void CEulerVariable::SetSecondaryVar_Compressible(CFluidModel *FluidModel) {
	unsigned short iVar;

   /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/

   SetdPdrho_e(FluidModel->GetdPdrho_e());
   SetdPde_rho(FluidModel->GetdPde_rho());

}

bool CEulerVariable::SetPrimVar_Incompressible(double Density_Inf, CConfig *config) {
  
  double ArtComp_Factor = config->GetArtComp_Factor();
  
  /*--- Set the value of the density ---*/
  
  SetDensityInc(Density_Inf);
  
  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/
  
  SetVelocityInc();
  
  /*--- Set the value of the pressure ---*/
  
	SetPressureInc();
  
  /*--- Set the value of the artificial compressibility factor ---*/
  
  SetBetaInc2(ArtComp_Factor);
  
  return true;
  
}

bool CEulerVariable::SetPrimVar_FreeSurface(CConfig *config) {
  
  double Heaviside, lambda, DensityInc, LevelSet, Distance;

  double ArtComp_Factor = config->GetArtComp_Factor();
  double epsilon = config->GetFreeSurface_Thickness();
  
  /*--- Set the value of the Level Set (already set in SetFreeSurface_Distance(geometry, config)) ---*/
  
  LevelSet = Primitive[nDim+5];

  /*--- Set the value of the Distance (already set in SetFreeSurface_Distance(geometry, config)) ---*/
  
  Distance = Primitive[nDim+6];

  /*--- Set the value of the Heaviside function ---*/

  Heaviside = 0.0;
  if (LevelSet < -epsilon) Heaviside = 1.0;
  if (fabs(LevelSet) <= epsilon) Heaviside = 1.0 - (0.5*(1.0+(LevelSet/epsilon)+(1.0/PI_NUMBER)*sin(PI_NUMBER*LevelSet/epsilon)));
  if (LevelSet > epsilon) Heaviside = 0.0;

  /*--- Set the value of the density ---*/

  lambda = config->GetRatioDensity();
  DensityInc = (lambda + (1.0 - lambda)*Heaviside)*config->GetDensity_FreeStreamND();
  SetDensityInc(DensityInc);
  
  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/
  
  SetVelocityInc();
  
  /*--- Set the value of the pressure ---*/
  
	SetPressureInc();
  
  /*--- Set the value of the artificial compressibility factor ---*/
  
  SetBetaInc2(ArtComp_Factor);
  
  return true;
  
}

CNSVariable::CNSVariable(void) : CEulerVariable() { }

CNSVariable::CNSVariable(double val_density, double *val_velocity, double val_energy,
                         unsigned short val_nDim, unsigned short val_nvar,
                         CConfig *config) : CEulerVariable(val_density, val_velocity, val_energy, val_nDim, val_nvar, config) {
  
	Temperature_Ref = config->GetTemperature_Ref();
	Viscosity_Ref   = config->GetViscosity_Ref();
	Viscosity_Inf   = config->GetViscosity_FreeStreamND();
	Prandtl_Lam     = config->GetPrandtl_Lam();
	Prandtl_Turb    = config->GetPrandtl_Turb();
  
}

CNSVariable::CNSVariable(double *val_solution, unsigned short val_nDim,
                         unsigned short val_nvar, CConfig *config) : CEulerVariable(val_solution, val_nDim, val_nvar, config) {
  
	Temperature_Ref = config->GetTemperature_Ref();
	Viscosity_Ref   = config->GetViscosity_Ref();
	Viscosity_Inf   = config->GetViscosity_FreeStreamND();
	Prandtl_Lam     = config->GetPrandtl_Lam();
	Prandtl_Turb    = config->GetPrandtl_Turb();
}

CNSVariable::~CNSVariable(void) { }

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
  
  double div;
  
  if (nDim == 2) {
    div = Gradient_Primitive[1][0] + Gradient_Primitive[2][1];
    StrainMag = 0.0;
    
    // add diagonals
    StrainMag += pow(Gradient_Primitive[1][0] - 1.0/3.0*div, 2.0);
    StrainMag += pow(Gradient_Primitive[2][1] - 1.0/3.0*div, 2.0);
    
    // add off diagonals
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][1] + Gradient_Primitive[2][0]), 2.0);
    
    StrainMag = sqrt(2.0*StrainMag);
    
  }
  else {
    div = Gradient_Primitive[1][0] + Gradient_Primitive[2][1] + Gradient_Primitive[3][2];
    StrainMag = 0.0;
    
    // add diagonals
    StrainMag += pow(Gradient_Primitive[1][0] - 1.0/3.0*div,2.0);
    StrainMag += pow(Gradient_Primitive[2][1] - 1.0/3.0*div,2.0);
    StrainMag += pow(Gradient_Primitive[3][2] - 1.0/3.0*div,2.0);
    
    // add off diagonals
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][1] + Gradient_Primitive[2][0]), 2.0);
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][2] + Gradient_Primitive[3][0]), 2.0);
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[2][2] + Gradient_Primitive[3][1]), 2.0);
    
    StrainMag = sqrt(2.0*StrainMag);
  }
  
}

bool CNSVariable::SetPrimVar_Compressible(double eddy_visc, double turb_ke, CFluidModel *FluidModel) {
	unsigned short iVar;
  bool check_dens = false, check_press = false, check_sos = false, check_temp = false, RightVol = true;
  
  
  SetVelocity();   // Computes velocity and velocity^2
  double density = GetDensity();
  double staticEnergy = GetEnergy()-0.5*Velocity2 - turb_ke;

  /* check will be moved inside fluid model plus error description strings*/
  FluidModel->SetTDState_rhoe(density, staticEnergy);

  check_dens = SetDensity();
  check_press = SetPressure(FluidModel->GetPressure());
  check_sos = SetSoundSpeed(FluidModel->GetSoundSpeed2());
  check_temp = SetTemperature(FluidModel->GetTemperature());
  
  /*--- Check that the solution has a physical meaning ---*/
  
  if (check_dens || check_press || check_sos || check_temp) {
    
    /*--- Copy the old solution ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    
    /*--- Recompute the primitive variables ---*/
    
    SetVelocity();   // Computes velocity and velocity^2
	double density = GetDensity();
	double staticEnergy = GetEnergy()-0.5*Velocity2 - turb_ke;
	/* check will be moved inside fluid model plus error description strings*/
	FluidModel->SetTDState_rhoe(density, staticEnergy);

	check_dens = SetDensity();
	check_press = SetPressure(FluidModel->GetPressure());
	check_sos = SetSoundSpeed(FluidModel->GetSoundSpeed2());
	check_temp = SetTemperature(FluidModel->GetTemperature());
    
    RightVol = false;
    
  }
  
  /*--- Set enthalpy ---*/
  
	SetEnthalpy();                                  // Requires pressure computation.
  
  /*--- Set laminar viscosity ---*/
  
	SetLaminarViscosity(FluidModel->GetLaminarViscosity(FluidModel->GetTemperature(), GetDensity()));
  
  /*--- Set eddy viscosity ---*/
  
  SetEddyViscosity(eddy_visc);
  
  return RightVol;
  
}

bool CNSVariable::SetPrimVar_Incompressible(double Density_Inf, double Viscosity_Inf, double eddy_visc, double turb_ke, CConfig *config) {
  
	double ArtComp_Factor = config->GetArtComp_Factor();
  
  /*--- Set the value of the density and viscosity ---*/
  
  SetDensityInc(Density_Inf);
  SetLaminarViscosityInc(Viscosity_Inf);
  
  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/
  
  SetVelocityInc();
  
  /*--- Set the value of the pressure ---*/
  
	SetPressureInc();
  
  /*--- Set the value of the artificial compressibility factor ---*/
  
  SetBetaInc2(ArtComp_Factor);
  
  /*--- Set eddy viscosity ---*/
  
  SetEddyViscosityInc(eddy_visc);
  
  return true;
  
}

bool CNSVariable::SetPrimVar_FreeSurface(double eddy_visc, double turb_ke, CConfig *config) {

  double Heaviside, lambda, DensityInc, ViscosityInc, LevelSet, Distance;
  
	double ArtComp_Factor = config->GetArtComp_Factor();
  double epsilon = config->GetFreeSurface_Thickness();

  /*--- Set the value of the Level Set (already set in SetFreeSurface_Distance(geometry, config)) ---*/
  
  LevelSet = Primitive[nDim+5];
  
  /*--- Set the value of the Distance (already set in SetFreeSurface_Distance(geometry, config)) ---*/
  
  Distance = Primitive[nDim+6];
  
  /*--- Set the value of the Heaviside function ---*/

  Heaviside = 0.0;
  if (LevelSet < -epsilon) Heaviside = 1.0;
  if (fabs(LevelSet) <= epsilon) Heaviside = 1.0 - (0.5*(1.0+(LevelSet/epsilon)+(1.0/PI_NUMBER)*sin(PI_NUMBER*LevelSet/epsilon)));
  if (LevelSet > epsilon) Heaviside = 0.0;
  
  /*--- Set the value of the density ---*/

  lambda = config->GetRatioDensity();
  DensityInc = (lambda + (1.0 - lambda)*Heaviside)*config->GetDensity_FreeStreamND();
  SetDensityInc(DensityInc);
  
  /*--- Set the value of the laminar viscosity ---*/

  lambda = config->GetRatioViscosity();
  ViscosityInc = (lambda + (1.0 - lambda)*Heaviside)*config->GetViscosity_FreeStreamND();
  SetLaminarViscosityInc(ViscosityInc);

  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/
  
  SetVelocityInc();
  
  /*--- Set the value of the pressure ---*/
  
	SetPressureInc();
  
  /*--- Set the value of the artificial compressibility factor ---*/
  
  SetBetaInc2(ArtComp_Factor);
  
  /*--- Set eddy viscosity ---*/
  
  SetEddyViscosityInc(eddy_visc);

  return true;
  
}
