/*!
 * \file variable_direct_tne2.cpp
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
  
  /*--- Array initialization ---*/
	Primitive = NULL;
	Gradient_Primitive = NULL;
	Limiter_Primitive = NULL;
  
  /*--- Acquire parameters from the config class ---*/
  nDim               = val_ndim;
  nSpecies           = config->GetnSpecies();
  ionization         = config->GetIonization();
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
    Solution_Old[iSpecies]    = val_density*val_massfrac[iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    Solution[nSpecies+iDim]     = val_density*val_velocity[iDim];
    Solution_Old[nSpecies+iDim] = val_density*val_velocity[iDim];
  }
  Solution[nSpecies+nDim]       = val_density*energy;
  Solution_Old[nSpecies+nDim]   = val_density*energy;
  Solution[nSpecies+nDim+1]     = val_density*energy_ve;
  Solution_Old[nSpecies+nDim+1] = val_density*energy_ve;
  
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

void CTNE2EulerVariable::SetDensity(void) {
  unsigned short iSpecies;
  double Density;
  
  Density = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Density += Solution[iSpecies];
  Primitive[nSpecies+nDim+3] = Density;
}

double CTNE2EulerVariable::GetProjVel(double *val_vector) {
	double ProjVel, density;
	unsigned short iDim, iSpecies;
  
	ProjVel = 0.0;
  density = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    density += Solution[iSpecies];
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVel += Solution[nSpecies+iDim]*val_vector[iDim]/density;
  
	return ProjVel;
}

void CTNE2EulerVariable::SetVelocity(double *val_velocity, bool val_incomp) {
	if (val_incomp) {
		for (unsigned short iDim = 0; iDim < nDim; iDim++)
			Solution[nSpecies+iDim] = val_velocity[iDim]*Primitive[nSpecies+nDim+3];
	}
	else {
		for (unsigned short iDim = 0; iDim < nDim; iDim++)
			Solution[nSpecies+iDim] = val_velocity[iDim]*Primitive[nSpecies+nDim+3];
	}
}

void CTNE2EulerVariable::SetVelocity2(void) {
  unsigned short iDim;
  
  Velocity2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity2 +=  Solution[nSpecies+iDim]*Solution[nSpecies+iDim]
    / (Primitive[nSpecies+nDim+3]*Primitive[nSpecies+nDim+3]);
  }
}

bool CTNE2EulerVariable::SetTemperature(CConfig *config) {
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
  Primitive[nSpecies] = (rhoe - rhoev - rhoeform + rhoeref - 0.5*rho*sqvel) / rhoCvtr;
  
  /*--- Assign vibrational-rotational temperature ---*/
  //NOTE:  NOT CORRECT!!!  NEED TO SOLVE FOR TVE PROPERLY!!!
  Primitive[nSpecies+1] = Primitive[nSpecies];
  
  if ((Primitive[nSpecies] > 0.0) && (Primitive[nSpecies+1])) return false;
  else return true;
}

bool CTNE2EulerVariable::SetPressure(CConfig *config) {
  unsigned short iSpecies;
  double *molarmass;
  double P;
  
  molarmass = config->GetMolar_Mass();
  P = 0.0;
  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    P += Solution[iSpecies]*UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Primitive[nSpecies];
  }
  if (ionization) {
    iSpecies = nSpecies - 1;
    P -= Solution[iSpecies]*UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Primitive[nSpecies];
    P += Solution[iSpecies]*UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Primitive[nSpecies+1];
  }
  
  Primitive[nSpecies+nDim+2] = P;
  
  if (Primitive[nSpecies+nDim+2] > 0.0) return false;
  else return true;
}

bool CTNE2EulerVariable::SetSoundSpeed(CConfig *config) {
  unsigned short iSpecies;
  double dPdrhoE, rhoCvtr, conc, radical;
  double *molarmass, *rotmodes;
  
  molarmass = config->GetMolar_Mass();
  rotmodes  = config->GetRotationModes();
  
  rhoCvtr = 0.0;
  conc    = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    rhoCvtr +=  Solution[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    conc    += Solution[iSpecies] / molarmass[iSpecies];
  }
  if (ionization) {
    iSpecies = nSpecies-1;
    rhoCvtr -=  Solution[iSpecies] * (3.0/2.0 + rotmodes[iSpecies]/2.0)
    * UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies];
    conc    -= Solution[iSpecies] / molarmass[iSpecies];
  }
  dPdrhoE = UNIVERSAL_GAS_CONSTANT/rhoCvtr * conc;
  
  /*--- Calculate a^2 using Gnoffo definition (NASA TP 2867) ---*/
  radical = (1.0 + dPdrhoE) * Primitive[nSpecies+nDim+2]/Primitive[nSpecies+nDim+3];
  
  if (radical < 0.0) return true;
  else {
    Primitive[nSpecies+nDim+5] = sqrt(radical);
    return false;
  }
}

void CTNE2EulerVariable::SetPrimVar_Compressible(CConfig *config) {
	unsigned short iDim, iVar, iSpecies;
  bool check_dens = false, check_press = false, check_sos = false, check_temp = false;
  
  /*--- Primitive variables [rho1,...,rhoNs,T,Tve,u,v,w,P,rho,h,c] ---*/
  
  SetDensity();                             // Compute mixture density
	SetVelocity2();                           // Compute the modulus of the velocity (req. mixture density).
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    check_dens = ((Solution[iSpecies] < 0.0) || check_dens);  // Check the density
  check_temp  = SetTemperature(config);     // Compute temperatures (T & Tve)
 	check_press = SetPressure(config);        // Requires T & Tve computation.
	check_sos   = SetSoundSpeed(config);      // Requires density & pressure computation.
  
  /*--- Check that the solution has a physical meaning ---*/
  if (check_dens || check_press || check_sos || check_temp) {
    
    /*--- Copy the old solution ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    
    /*--- Recompute the primitive variables ---*/
    SetDensity();                           // Compute mixture density
    SetVelocity2();                         // Compute square of the velocity (req. mixture density).
    check_temp  = SetTemperature(config);   // Compute temperatures (T & Tve)
    check_press = SetPressure(config);      // Requires T & Tve computation.
    check_sos   = SetSoundSpeed(config);    // Requires density & pressure computation.
  }
  SetEnthalpy();                            // Requires density & pressure computation.
  
  /*--- Calculate velocities (req. density calculation) ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		Primitive[nSpecies+iDim+2] = Solution[nSpecies+iDim] / Primitive[nSpecies+nDim+3];
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

bool CTNE2NSVariable::SetPressure(CConfig *config) {
  
  /////////////////////////////////////////////////////////
  // SUBTRACT TURBULENT KINETIC ENERGY FROM PRESSURE???? //
  /////////////////////////////////////////////////////////
  
  unsigned short iSpecies;
  double *molarmass;
  double P;
  
  molarmass = config->GetMolar_Mass();
  P = 0.0;
  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    P += Solution[iSpecies]*UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Primitive[nSpecies];
  }
  if (ionization) {
    iSpecies = nSpecies - 1;
    P -= Solution[iSpecies]*UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Primitive[nSpecies];
    P += Solution[iSpecies]*UNIVERSAL_GAS_CONSTANT/molarmass[iSpecies] * Primitive[nSpecies+1];
  }
  
  Primitive[nSpecies+nDim+2] = P;
  
  if (Primitive[nSpecies+nDim+2] > 0.0) return false;
  else return true;
}

void CTNE2NSVariable::SetPrimVar_Compressible(CConfig *config) {
	unsigned short iDim, iVar, iSpecies;
  bool check_dens = false, check_press = false, check_sos = false, check_temp = false;
  
  /*--- Primitive variables [rho1,...,rhoNs,T,Tve,u,v,w,P,rho,h,c] ---*/
  
  SetDensity();                             // Compute mixture density
	SetVelocity2();                           // Compute the modulus of the velocity (req. mixture density).
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    check_dens = ((Solution[iSpecies] < 0.0) || check_dens);  // Check the density
  check_temp  = SetTemperature(config);     // Compute temperatures (T & Tve)
 	check_press = SetPressure(config);        // Requires T & Tve computation.
	check_sos   = SetSoundSpeed(config);      // Requires density & pressure computation.
  
  /*--- Check that the solution has a physical meaning ---*/
  if (check_dens || check_press || check_sos || check_temp) {
    
    /*--- Copy the old solution ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    
    /*--- Recompute the primitive variables ---*/
    SetDensity();                           // Compute mixture density
    SetVelocity2();                         // Compute square of the velocity (req. mixture density).
    check_temp  = SetTemperature(config);   // Compute temperatures (T & Tve)
    check_press = SetPressure(config);      // Requires T & Tve computation.
    check_sos   = SetSoundSpeed(config);    // Requires density & pressure computation.
  }
  SetEnthalpy();                            // Requires density & pressure computation.
  SetLaminarViscosity();                    // Requires temperature computation.
  
  /*--- Calculate velocities (req. density calculation) ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		Primitive[nSpecies+iDim+2] = Solution[nSpecies+iDim] / Primitive[nSpecies+nDim+3];
  
}