/*!
 * \file variable_direct_tne2.cpp
 * \brief Definition of the solution fields.
 * \author S. Copeland
 * \version 4.0.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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
#include <math.h>

CTNE2EulerVariable::CTNE2EulerVariable(void) : CVariable() {  

  /*--- Array initialization ---*/
	Primitive = NULL;
	Gradient_Primitive = NULL;
	Limiter_Primitive = NULL;
  dPdU = NULL;
  dTdU = NULL;
  dTvedU = NULL;
  
  /*--- Define structure of the primtive variable vector ---*/
  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  RHOS_INDEX    = 0;
  T_INDEX       = nSpecies;
  TVE_INDEX     = nSpecies+1;
  VEL_INDEX     = nSpecies+2;
  P_INDEX       = nSpecies+nDim+2;
  RHO_INDEX     = nSpecies+nDim+3;
  H_INDEX       = nSpecies+nDim+4;
  A_INDEX       = nSpecies+nDim+5;
  RHOCVTR_INDEX = nSpecies+nDim+6;
  RHOCVVE_INDEX = nSpecies+nDim+7;

}

CTNE2EulerVariable::CTNE2EulerVariable(unsigned short val_nDim,
                                       unsigned short val_nvar,
                                       unsigned short val_nprimvar,
                                       unsigned short val_nprimvargrad,
                                       CConfig *config) : CVariable(val_nDim,
                                                                    val_nvar,
                                                                    config) {
  
  nDim         = val_nDim;
  nVar         = val_nvar;
  nPrimVar     = val_nprimvar;
  nPrimVarGrad = val_nprimvargrad;
  nSpecies = config->GetnSpecies();
  
  ionization = config->GetIonization();
  
  /*--- Array initialization ---*/
	Primitive = NULL;
	Gradient_Primitive = NULL;
	Limiter_Primitive = NULL;
  dPdU = NULL;
  dTdU = NULL;
  dTvedU = NULL;
  
  /*--- Define structure of the primtive variable vector ---*/
  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  RHOS_INDEX    = 0;
  T_INDEX       = nSpecies;
  TVE_INDEX     = nSpecies+1;
  VEL_INDEX     = nSpecies+2;
  P_INDEX       = nSpecies+nDim+2;
  RHO_INDEX     = nSpecies+nDim+3;
  H_INDEX       = nSpecies+nDim+4;
  A_INDEX       = nSpecies+nDim+5;
  RHOCVTR_INDEX = nSpecies+nDim+6;
  RHOCVVE_INDEX = nSpecies+nDim+7;
  
}


CTNE2EulerVariable::CTNE2EulerVariable(su2double val_pressure,
                                       su2double *val_massfrac,
                                       su2double *val_mach,
                                       su2double val_temperature,
                                       su2double val_temperature_ve,
                                       unsigned short val_nDim,
                                       unsigned short val_nvar,
                                       unsigned short val_nvarprim,
                                       unsigned short val_nvarprimgrad,
                                       CConfig *config) : CVariable(val_nDim,
                                                                    val_nvar,
                                                                    config   ) {
  
  unsigned short iEl, iMesh, iDim, iSpecies, iVar, nDim, nEl, nHeavy, nMGSmooth;
  unsigned short *nElStates;
  su2double *xi, *Ms, *thetav, **thetae, **g, *hf, *Tref;
  su2double rhoE, rhoEve, Ev, Ee, Ef, T, Tve, rho, rhoCvtr, rhos;
  su2double Ru, sqvel, num, denom, conc, soundspeed;
  
  /*--- Get Mutation++ mixture ---*/
  nSpecies     = config->GetnSpecies();
  nDim         = val_nDim;
  nPrimVar     = val_nvarprim;
  nPrimVarGrad = val_nvarprimgrad;
  nMGSmooth    = 0;
  
  /*--- Define structure of the primtive variable vector ---*/
  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  RHOS_INDEX    = 0;
  T_INDEX       = nSpecies;
  TVE_INDEX     = nSpecies+1;
  VEL_INDEX     = nSpecies+2;
  P_INDEX       = nSpecies+nDim+2;
  RHO_INDEX     = nSpecies+nDim+3;
  H_INDEX       = nSpecies+nDim+4;
  A_INDEX       = nSpecies+nDim+5;
  RHOCVTR_INDEX = nSpecies+nDim+6;
  RHOCVVE_INDEX = nSpecies+nDim+7;
  
  /*--- Array initialization ---*/
	Primitive = NULL;
	Gradient_Primitive = NULL;
	Limiter_Primitive = NULL;
  
  /*--- Allocate & initialize residual vectors ---*/
	Res_TruncError = new su2double [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Res_TruncError[iVar] = 0.0;
	}
  
	/*--- If using multigrid, allocate residual-smoothing vectors ---*/
	for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
		nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
	if (nMGSmooth > 0) {
		Residual_Sum = new su2double [nVar];
		Residual_Old = new su2double [nVar];
	}
  
  /*--- If using limiters, allocate the arrays ---*/
  if ((config->GetKind_ConvNumScheme_TNE2() == SPACE_UPWIND) &&
			(config->GetSpatialOrder_TNE2() == SECOND_ORDER_LIMITER)) {
		Limiter      = new su2double [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Limiter[iVar]      = 0.0;
		}
	}
  
  /*--- Allocate & initialize primitive variable & gradient arrays ---*/
  Primitive         = new su2double [nPrimVar];
  Limiter_Primitive = new su2double [nPrimVarGrad];
  Solution_Max      = new su2double [nPrimVarGrad];
  Solution_Min      = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;
  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Limiter_Primitive[iVar] = 0.0;
    Solution_Min[iVar]      = 0.0;
    Solution_Max[iVar]      = 0.0;
    Gradient_Primitive[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }
  
  /*--- Allocate partial derivative vectors ---*/
  dPdU      = new su2double [nVar];
  dTdU      = new su2double [nVar];
  dTvedU    = new su2double [nVar];
  
  /*--- Determine the number of heavy species ---*/
  ionization = config->GetIonization();
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Load variables from the config class --*/
  xi        = config->GetRotationModes();      // Rotational modes of energy storage
  Ms        = config->GetMolar_Mass();         // Species molar mass
  thetav    = config->GetCharVibTemp();        // Species characteristic vib. temperature [K]
  thetae    = config->GetCharElTemp();         // Characteristic electron temperature [K]
  g         = config->GetElDegeneracy();       // Degeneracy of electron states
  nElStates = config->GetnElStates();          // Number of electron states
  Tref      = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf        = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]
  
  /*--- Rename & initialize for convenience ---*/
  Ru      = UNIVERSAL_GAS_CONSTANT;         // Universal gas constant [J/(kmol*K)]
  Tve     = val_temperature_ve;             // Vibrational temperature [K]
  T       = val_temperature;                // Translational-rotational temperature [K]
  sqvel   = 0.0;                            // Velocity^2 [m2/s2]
  rhoE    = 0.0;                            // Mixture total energy per mass [J/kg]
  rhoEve  = 0.0;                            // Mixture vib-el energy per mass [J/kg]
  denom   = 0.0;
  conc    = 0.0;
  rhoCvtr = 0.0;
  
  /*--- Calculate mixture density from supplied primitive quantities ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    denom += val_massfrac[iSpecies] * (Ru/Ms[iSpecies]) * T;
  for (iSpecies = 0; iSpecies < nEl; iSpecies++)
    denom += val_massfrac[nSpecies-1] * (Ru/Ms[nSpecies-1]) * Tve;
  rho = val_pressure / denom;
  
  /*--- Calculate sound speed and extract velocities ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    conc += val_massfrac[iSpecies]*rho/Ms[iSpecies];
    rhoCvtr += rho*val_massfrac[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
  }
  soundspeed = sqrt((1.0 + Ru/rhoCvtr*conc) * val_pressure/rho);
  for (iDim = 0; iDim < nDim; iDim++)
    sqvel += val_mach[iDim]*soundspeed * val_mach[iDim]*soundspeed;
  
  /*--- Calculate energy (RRHO) from supplied primitive quanitites ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    // Species density
    rhos = val_massfrac[iSpecies]*rho;
    
    // Species formation energy
    Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];
    
    // Species vibrational energy
    if (thetav[iSpecies] != 0.0)
      Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Tve)-1.0);
    else
      Ev = 0.0;
    
    // Species electronic energy    
    num = 0.0;
    denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);
    for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
    }
    Ee = Ru/Ms[iSpecies] * (num/denom);
    
    // Mixture total energy
    rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (T-Tref[iSpecies])
                    + Ev + Ee + Ef + 0.5*sqvel);
    
    // Mixture vibrational-electronic energy
    rhoEve += rhos * (Ev + Ee);
  }
  for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
    // Species formation energy
    Ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1] * Tref[nSpecies-1];
    
    // Electron t-r mode contributes to mixture vib-el energy
    rhoEve += (3.0/2.0) * Ru/Ms[nSpecies-1] * (Tve - Tref[nSpecies-1]);
  }
  
  /*--- Initialize Solution & Solution_Old vectors ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Solution[iSpecies]     = rho*val_massfrac[iSpecies];
    Solution_Old[iSpecies] = rho*val_massfrac[iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    Solution[nSpecies+iDim]     = rho*val_mach[iDim]*soundspeed;
    Solution_Old[nSpecies+iDim] = rho*val_mach[iDim]*soundspeed;
  }
  Solution[nSpecies+nDim]       = rhoE;
  Solution_Old[nSpecies+nDim]   = rhoE;
  Solution[nSpecies+nDim+1]     = rhoEve;
  Solution_Old[nSpecies+nDim+1] = rhoEve;
  
  /*--- Assign primitive variables ---*/
  Primitive[T_INDEX]   = val_temperature;
  Primitive[TVE_INDEX] = val_temperature_ve;
  Primitive[P_INDEX]   = val_pressure;
  
}

CTNE2EulerVariable::CTNE2EulerVariable(su2double *val_solution,
                                       unsigned short val_nDim,
                                       unsigned short val_nvar,
                                       unsigned short val_nvarprim,
                                       unsigned short val_nvarprimgrad,
                                       CConfig *config) : CVariable(val_nDim,
                                                                    val_nvar,
                                                                    config) {

	unsigned short iVar, iDim, iMesh, nMGSmooth;
  
  nSpecies     = config->GetnSpecies();
  nDim         = val_nDim;
  nPrimVar     = val_nvarprim;
  nPrimVarGrad = val_nvarprimgrad;
  nMGSmooth    = 0;
  
  /*--- Define structure of the primtive variable vector ---*/
  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  RHOS_INDEX    = 0;
  T_INDEX       = nSpecies;
  TVE_INDEX     = nSpecies+1;
  VEL_INDEX     = nSpecies+2;
  P_INDEX       = nSpecies+nDim+2;
  RHO_INDEX     = nSpecies+nDim+3;
  H_INDEX       = nSpecies+nDim+4;
  A_INDEX       = nSpecies+nDim+5;
  RHOCVTR_INDEX = nSpecies+nDim+6;
  RHOCVVE_INDEX = nSpecies+nDim+7;
  
  /*--- Array initialization ---*/
	Primitive = NULL;
	Gradient_Primitive = NULL;
  Limiter_Primitive = NULL;
  
  /*--- Allocate & initialize residual vectors ---*/
	Res_TruncError = new su2double [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Res_TruncError[iVar] = 0.0;
	}
  
	/*--- If using multigrid, allocate residual-smoothing vectors ---*/
	for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
		nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
	if (nMGSmooth > 0) {
		Residual_Sum = new su2double [nVar];
		Residual_Old = new su2double [nVar];
	}
  
  /*--- If using limiters, allocate the arrays ---*/
  if ((config->GetKind_ConvNumScheme_TNE2() == SPACE_UPWIND) &&
			(config->GetSpatialOrder_TNE2() == SECOND_ORDER_LIMITER)) {
		Limiter           = new su2double [nVar];
    Limiter_Primitive = new su2double [nPrimVarGrad];
		Solution_Max      = new su2double [nPrimVarGrad];
		Solution_Min      = new su2double [nPrimVarGrad];
		for (iVar = 0; iVar < nVar; iVar++)
			Limiter[iVar]      = 0.0;
    for (iVar = 0; iVar < nVar; iVar++) {
      Limiter_Primitive[iVar] = 0.0;
			Solution_Max[iVar]      = 0.0;
			Solution_Min[iVar]      = 0.0;
    }
	}
  
  /*--- Allocate & initialize primitive variable & gradient arrays ---*/
  Primitive = new su2double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;
  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }
  
  /*--- Allocate partial derivative vectors ---*/
  dPdU      = new su2double [nVar];
  dTdU      = new su2double [nVar];
  dTvedU    = new su2double [nVar];
  
  /*--- Determine the number of heavy species ---*/
  ionization = config->GetIonization();
  
	/*--- Initialize Solution & Solution_Old vectors ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar]     = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}
  
  /*--- Initialize Tve to the free stream for Newton-Raphson method ---*/
  Primitive[TVE_INDEX] = config->GetTemperature_FreeStream();
  Primitive[T_INDEX]   = config->GetTemperature_FreeStream();
  Primitive[P_INDEX]   = config->GetPressure_FreeStreamND();
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
  if (dPdU      != NULL) delete [] dPdU;
  if (dTdU      != NULL) delete [] dTdU;
  if (dTvedU    != NULL) delete [] dTvedU;
  
}

void CTNE2EulerVariable::SetGradient_PrimitiveZero(unsigned short val_primvar) {
	unsigned short iVar, iDim;
  
	for (iVar = 0; iVar < val_primvar; iVar++)
		for (iDim = 0; iDim < nDim; iDim++)
			Gradient_Primitive[iVar][iDim] = 0.0;
}

bool CTNE2EulerVariable::SetDensity(void) {
  unsigned short iSpecies;
  su2double Density;
  
  Density = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Primitive[RHOS_INDEX+iSpecies] = Solution[iSpecies];
    Density += Solution[iSpecies];
  }
  Primitive[RHO_INDEX] = Density;
  
  return true;
}

su2double CTNE2EulerVariable::GetProjVel(su2double *val_vector) {
	su2double ProjVel, density;
	unsigned short iDim, iSpecies;
  
	ProjVel = 0.0;
  density = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    density += Solution[iSpecies];
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVel += Solution[nSpecies+iDim]*val_vector[iDim]/density;
  
	return ProjVel;
}

void CTNE2EulerVariable::SetVelocity2(void) {
  unsigned short iDim;
  
  Velocity2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Primitive[VEL_INDEX+iDim] = Solution[nSpecies+iDim] / Primitive[RHO_INDEX];
    Velocity2 +=  Solution[nSpecies+iDim]*Solution[nSpecies+iDim]
    / (Primitive[RHO_INDEX]*Primitive[RHO_INDEX]);
  }
}

bool CTNE2EulerVariable::SetTemperature(CConfig *config) {
  
  // Note: Requires previous call to SetDensity()
  unsigned short iEl, iSpecies, iDim, nHeavy, nEl, iIter, maxIter, *nElStates;
  su2double *xi, *Ms, *thetav, **thetae, **g, *hf, *Tref;
  su2double rho, rhoE, rhoEve, rhoEve_t, rhoE_ref, rhoE_f;
  su2double evs, eels;
  su2double Ru, sqvel, rhoCvtr, rhoCvve = 0.0;
  su2double Cvvs, Cves, Tve, Tve2 = 0.0, Tve_o;
  su2double f, df, tol;
  su2double exptv, thsqr, thoTve;
  su2double num, denom, num2, num3;
  su2double Tmin, Tmax, Tvemin, Tvemax;
  
  /*--- Set tolerance for Newton-Raphson method ---*/
  tol     = 1.0E-4;
  maxIter = 100;
  
  /*--- Set temperature clipping values ---*/
  Tmin   = 100.0;
  Tmax   = 6E4;
  Tvemin = 100.0;
  Tvemax = 4E4;
  
  /*--- Determine the number of heavy species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
  
  /*--- Load variables from the config class --*/
  xi        = config->GetRotationModes();      // Rotational modes of energy storage
  Ms        = config->GetMolar_Mass();         // Species molar mass
  thetav    = config->GetCharVibTemp();        // Species characteristic vib. temperature [K]
  Tref      = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf        = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]
  thetae    = config->GetCharElTemp();
  g         = config->GetElDegeneracy();
  nElStates = config->GetnElStates();
  
  /*--- Rename & initialize for convenience ---*/
  Ru       = UNIVERSAL_GAS_CONSTANT;           // Universal gas constant [J/(kmol*K)]
  rho      = Primitive[RHO_INDEX];             // Mixture density [kg/m3]
  rhoE     = Solution[nSpecies+nDim];          // Density * energy [J/m3]
  rhoEve   = Solution[nSpecies+nDim+1];        // Density * energy_ve [J/m3]
  rhoE_f   = 0.0;                              // Density * formation energy [J/m3]
  rhoE_ref = 0.0;                              // Density * reference energy [J/m3]
  rhoCvtr  = 0.0;                              // Mix spec. heat @ const. volume [J/(kg*K)]
  sqvel    = 0.0;                              // Velocity^2 [m2/s2]
  
  /*--- Error checking ---*/
  if (rhoE < 0.0)
    rhoE = EPS;
  if (rhoEve < 0.0)
    rhoEve = EPS;
  
  /*--- Calculate mixture properties (heavy particles only) ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhoCvtr  += Solution[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
    rhoE_ref += Solution[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * Tref[iSpecies];
    rhoE_f   += Solution[iSpecies] * (hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies]);
  }
  for (iDim = 0; iDim < nDim; iDim++)
    sqvel    += (Solution[nSpecies+iDim]/rho) * (Solution[nSpecies+iDim]/rho);
  
  /*--- Calculate translational-rotational temperature ---*/
  Primitive[T_INDEX] = (rhoE - rhoEve - rhoE_f + rhoE_ref - 0.5*rho*sqvel) / rhoCvtr;
  
  /*--- Calculate vibrational-electronic temperature ---*/
  // NOTE: Cannot write an expression explicitly in terms of Tve
  // NOTE: We use Newton-Raphson to iteratively find the value of Tve
  // NOTE: Use T as an initial guess
  Tve   = Primitive[TVE_INDEX];
  Tve_o = Primitive[TVE_INDEX];
  
  for (iIter = 0; iIter < maxIter; iIter++) {
    rhoEve_t = 0.0;
    rhoCvve  = 0.0;
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      
      /*--- Vibrational energy ---*/
      if (thetav[iSpecies] != 0.0) {
        
        /*--- Rename for convenience ---*/
        thoTve = thetav[iSpecies]/Tve;
        exptv = exp(thetav[iSpecies]/Tve);
        thsqr = thetav[iSpecies]*thetav[iSpecies];
        
        /*--- Calculate vibrational energy ---*/
        evs  = Ru/Ms[iSpecies] * thetav[iSpecies] / (exptv - 1.0);
        
        /*--- Calculate species vibrational specific heats ---*/
        Cvvs  = Ru/Ms[iSpecies] * thoTve*thoTve * exptv / ((exptv-1.0)*(exptv-1.0));
        
        /*--- Add contribution ---*/
        rhoEve_t += Solution[iSpecies] * evs;
        rhoCvve  += Solution[iSpecies] * Cvvs;
      }
      /*--- Electronic energy ---*/
      if (nElStates[iSpecies] != 0) {
        num = 0.0; num2 = 0.0;
        denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);
        num3  = g[iSpecies][0] * (thetae[iSpecies][0]/(Tve*Tve))*exp(-thetae[iSpecies][0]/Tve);
        for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
          thoTve = thetae[iSpecies][iEl]/Tve;
          exptv = exp(-thetae[iSpecies][iEl]/Tve);
          
          num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exptv;
          denom += g[iSpecies][iEl] * exptv;
          num2  += g[iSpecies][iEl] * (thoTve*thoTve) * exptv;
          num3  += g[iSpecies][iEl] * thoTve/Tve * exptv;
        }
        eels = Ru/Ms[iSpecies] * (num/denom);
        Cves = Ru/Ms[iSpecies] * (num2/denom - num*num3/(denom*denom));
        
        rhoEve_t += Solution[iSpecies] * eels;
        rhoCvve  += Solution[iSpecies] * Cves;
      }
    }
    for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
      Cves = 3.0/2.0 * Ru/Ms[nSpecies-1];
      rhoEve_t += Solution[nSpecies-1] * Cves * Tve;
      rhoCvve += Solution[nSpecies-1] * Cves;
    }
    
    /*--- Determine function f(Tve) and df/dTve ---*/
    f  = rhoEve - rhoEve_t;
    df = -rhoCvve;
    Tve2 = Tve - (f/df)*0.5;
    
    /*--- Check for non-physical conditions ---*/
    if ((Tve2 != Tve2) || (Tve2 < 0))
      Tve2 = 1.4*Primitive[T_INDEX];
    
    /*--- Check for convergence ---*/
    if (fabs(Tve2-Tve) < tol) break;
    if (iIter == maxIter-1) {
      cout << "WARNING!!! Tve convergence not reached!" << endl;
      cout << "rhoE: " << rhoE << endl;
      cout << "rhoEve: " << rhoEve << endl;
      cout << "T: " << Primitive[T_INDEX] << endl;
      cout << "Tve2: " << Tve2 << endl;
      cout << "Tve_o: " << Tve_o << endl;
      Tve2 = Tve_o;
      break;
    }
    Tve = Tve2;
  }
  
  Primitive[TVE_INDEX] = Tve2;
  
  /*--- Error checking ---*/
  if (Primitive[T_INDEX] <= Tmin) {
    cout << "WARNING: T = " << Primitive[T_INDEX] << "\t -- Clipping at: " << Tmin << endl;
    Primitive[T_INDEX] = Tmin;
  } else if (Primitive[T_INDEX] >= Tmax) {
    cout << "WARNING: T = " << Primitive[T_INDEX] << "\t -- Clipping at: " << Tmax << endl;
    Primitive[T_INDEX] = Tmax;
  }
  if (Primitive[TVE_INDEX] <= Tvemin) {
    cout << "WARNING: Tve = " << Primitive[TVE_INDEX] << "\t -- Clipping at: " << Tvemin << endl;
    Primitive[TVE_INDEX] = Tvemin;
  } else if (Primitive[TVE_INDEX] >= Tvemax) {
    cout << "WARNING: Tve = " << Primitive[TVE_INDEX] << "\t -- Clipping at: " << Tvemax << endl;
    Primitive[TVE_INDEX] = Tvemax;
  }
    
  
  /*--- Assign Gas Properties ---*/
  Primitive[RHOCVTR_INDEX] = rhoCvtr;
  Primitive[RHOCVVE_INDEX] = rhoCvve;
  
  /*--- Check that the solution is physical ---*/
  if ((Primitive[T_INDEX] > 0.0) && (Primitive[TVE_INDEX] > 0.0)) return false;
  else return true;
}


void CTNE2EulerVariable::SetGasProperties(CConfig *config) {

  // NOTE: Requires computation of vib-el temperature.
  // NOTE: For now, neglecting contribution from electronic excitation
  
  unsigned short iSpecies, nHeavy, nEl;
  su2double rhoCvtr, rhoCvve, Ru, Tve;
  su2double *xi, *Ms, *thetav;
  
  /*--- Determine the number of heavy species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
  
  /*--- Load variables from the config class --*/
  xi     = config->GetRotationModes(); // Rotational modes of energy storage
  Ms     = config->GetMolar_Mass();    // Species molar mass
  thetav = config->GetCharVibTemp();   // Species characteristic vib. temperature [K]
  Ru     = UNIVERSAL_GAS_CONSTANT;     // Universal gas constant [J/(kmol*K)]
  Tve    = Primitive[nSpecies+1];      // Vibrational-electronic temperature [K]
  
  rhoCvtr = 0.0;
  rhoCvve = 0.0;
  
  /*--- Heavy particle contribution ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhoCvtr += Solution[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
    
    if (thetav[iSpecies] != 0.0) {
      rhoCvve += Ru/Ms[iSpecies] * (thetav[iSpecies]/Tve)*(thetav[iSpecies]/Tve) * exp(thetav[iSpecies]/Tve)
                 / ((exp(thetav[iSpecies]/Tve)-1.0)*(exp(thetav[iSpecies]/Tve)-1.0));
    }
    
    //--- [ Electronic energy goes here ] ---//
    
  }
  
  /*--- Free-electron contribution ---*/
  for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
    rhoCvve += Solution[nSpecies-1] * 3.0/2.0 * Ru/Ms[nSpecies-1];
  }
  
  /*--- Store computed values ---*/
  Primitive[RHOCVTR_INDEX] = rhoCvtr;
  Primitive[RHOCVVE_INDEX] = rhoCvve;
}


bool CTNE2EulerVariable::SetPressure(CConfig *config) {
  
  // NOTE: Requires computation of trans-rot & vib-el temperatures.
  
  unsigned short iSpecies, nHeavy, nEl;
  su2double *Ms;
  su2double P, Ru;
  
  /*--- Determine the number of heavy species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
  
  /*--- Read gas mixture properties from config ---*/
  Ms = config->GetMolar_Mass();
  
  /*--- Rename for convenience ---*/
  Ru = UNIVERSAL_GAS_CONSTANT;
  
  /*--- Solve for mixture pressure using ideal gas law & Dalton's law ---*/
  // Note: If free electrons are present, use Tve for their partial pressure
  P = 0.0;
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    P += Solution[iSpecies] * Ru/Ms[iSpecies] * Primitive[T_INDEX];
  
  for (iSpecies = 0; iSpecies < nEl; iSpecies++)
    P += Solution[nSpecies-1] * Ru/Ms[nSpecies-1] * Primitive[TVE_INDEX];
  
  /*--- Store computed values and check for a physical solution ---*/
  Primitive[P_INDEX] = P;
  if (Primitive[P_INDEX] > 0.0) return false;
  else return true;
}

bool CTNE2EulerVariable::SetSoundSpeed(CConfig *config) {
  
  // NOTE: Requires SetDensity(), SetTemperature(), SetPressure(), & SetGasProperties().
  
//  unsigned short iDim, iSpecies, nHeavy, nEl;
//  su2double dPdrhoE, factor, Ru, radical, radical2;
//  su2double *Ms, *xi;
//  
//  /*--- Determine the number of heavy species ---*/
//  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
//  else            { nHeavy = nSpecies;   nEl = 0; }
//  
//  /*--- Read gas mixture properties from config ---*/
//  Ms = config->GetMolar_Mass();
//  xi = config->GetRotationModes();
//  
//  /*--- Rename for convenience ---*/
//  Ru = UNIVERSAL_GAS_CONSTANT;
//  
//  /*--- Calculate partial derivative of pressure w.r.t. rhoE ---*/
//  factor = 0.0;
//  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
//    factor += Solution[iSpecies] / Ms[iSpecies];
//  dPdrhoE = Ru/Primitive[RHOCVTR_INDEX] * factor;
//  
//  /*--- Calculate a^2 using Gnoffo definition (NASA TP 2867) ---*/
//  radical = (1.0+dPdrhoE) * Primitive[P_INDEX]/Primitive[RHO_INDEX];
  
  
  
  /// alternative
  unsigned short iSpecies, iDim;
  su2double radical2;
  
  radical2 = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    radical2 += Primitive[RHOS_INDEX+iSpecies]/Primitive[RHO_INDEX] * dPdU[iSpecies];
  for (iDim = 0; iDim < nDim; iDim++)
    radical2 += Primitive[VEL_INDEX+iDim]*dPdU[nSpecies+iDim];
  radical2 += (Solution[nSpecies+nDim]+Primitive[P_INDEX])/Primitive[RHO_INDEX] * dPdU[nSpecies+nDim];
  radical2 += Solution[nSpecies+nDim+1]/Primitive[RHO_INDEX] * dPdU[nSpecies+nDim+1];
  
  if (radical2 < 0.0) return true;
  else {Primitive[A_INDEX] = sqrt(radical2); return false; }
  
//  if (radical < 0.0) return true;
//  else { Primitive[A_INDEX] = sqrt(radical); return false; }
}

void CTNE2EulerVariable::CalcdPdU(su2double *V, CConfig *config, su2double *val_dPdU) {

  // Note: Requires SetDensity(), SetTemperature(), SetPressure(), & SetGasProperties()
  // Note: Electron energy not included properly.
  
  unsigned short iDim, iSpecies, iEl, nHeavy, nEl, *nElStates;
  su2double *Ms, *Tref, *hf, *xi, *thetav, **thetae, **g;
  su2double Ru, RuBAR, CvtrBAR, rhoCvtr, rhoCvve, Cvtrs, rho_el, sqvel, conc;
  su2double rho, rhos, T, Tve, evibs, eels, ef;
  su2double num, denom;
  
  if (val_dPdU == NULL) {
    cout << "ERROR: CalcdPdU - Array dPdU not allocated!" << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Determine the number of heavy species ---*/
  if (ionization) {
    nHeavy = nSpecies-1;
    nEl    = 1;
    rho_el = V[RHOS_INDEX+nSpecies-1];
  } else {
    nHeavy = nSpecies;
    nEl    = 0;
    rho_el = 0.0;
  }
  
  /*--- Read gas mixture properties from config ---*/
  Ms        = config->GetMolar_Mass();
  Tref      = config->GetRefTemperature();
  hf        = config->GetEnthalpy_Formation();
  xi        = config->GetRotationModes();
  thetav    = config->GetCharVibTemp();
  thetae    = config->GetCharElTemp();
  g         = config->GetElDegeneracy();
  nElStates = config->GetnElStates();
  
  /*--- Rename for convenience ---*/
  Ru      = UNIVERSAL_GAS_CONSTANT;
  T       = V[T_INDEX];
  Tve     = V[TVE_INDEX];
  rho     = V[RHO_INDEX];
  rhoCvtr = V[RHOCVTR_INDEX];
  rhoCvve = V[RHOCVVE_INDEX];
  
  /*--- Pre-compute useful quantities ---*/
  RuBAR   = 0.0;
  CvtrBAR = 0.0;
  sqvel   = 0.0;
  conc    = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    sqvel += V[VEL_INDEX+iDim] * V[VEL_INDEX+iDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    CvtrBAR += V[RHOS_INDEX+iSpecies]*(3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
    conc    += V[RHOS_INDEX+iSpecies]/Ms[iSpecies];
  }
  
  // Species density
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhos  = V[RHOS_INDEX+iSpecies];
    ef    = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];
    Cvtrs = (3.0/2.0+xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
    
    val_dPdU[iSpecies] =  T*Ru/Ms[iSpecies] + Ru*conc/rhoCvtr *
                                          (-Cvtrs*(T-Tref[iSpecies]) -
                                           ef + 0.5*sqvel);
  }
  if (ionization) {
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      evibs = Ru/Ms[iSpecies] * thetav[iSpecies]/(exp(thetav[iSpecies]/Tve)-1.0);
      num = 0.0;
      denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);
      for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
        num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
        denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      }
      eels = Ru/Ms[iSpecies] * (num/denom);
      
      val_dPdU[iSpecies] -= rho_el * Ru/Ms[nSpecies-1] * (evibs + eels)/rhoCvve;
    }
    ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1]*Tref[nSpecies-1];
    val_dPdU[nSpecies-1] = Ru*conc/rhoCvtr * (-ef + 0.5*sqvel)
                         + Ru/Ms[nSpecies-1]*Tve
                         - rho_el*Ru/Ms[nSpecies-1] * (-3.0/2.0*Ru/Ms[nSpecies-1]*Tve)/rhoCvve;
  }
  // Momentum
  for (iDim = 0; iDim < nDim; iDim++)
    val_dPdU[nSpecies+iDim] = -conc*Ru*V[VEL_INDEX+iDim]/rhoCvtr;
  
  // Total energy
  val_dPdU[nSpecies+nDim]   = conc*Ru / rhoCvtr;
  
  // Vib.-el energy
  val_dPdU[nSpecies+nDim+1] = -val_dPdU[nSpecies+nDim]
                            + rho_el*Ru/Ms[nSpecies-1]*1.0/rhoCvve;
}

su2double CTNE2EulerVariable::CalcEve(su2double *V, CConfig *config, unsigned short val_Species) {

  unsigned short iEl, *nElStates;
  su2double *Ms, *thetav, **thetae, **g, *hf, *Tref, Ru;
  su2double Tve, Ev, Eel, Ef;
  su2double num, denom;
  
  /*--- Read gas mixture properties from config ---*/
  Ms        = config->GetMolar_Mass();
  
  /*--- Rename for convenience ---*/
  Ru  = UNIVERSAL_GAS_CONSTANT;
  Tve = V[TVE_INDEX];
  
  /*--- Electron species energy ---*/
  if ((ionization) && (val_Species == nSpecies-1)) {
    
    /*--- Get quantities from CConfig ---*/
    Tref = config->GetRefTemperature();
    hf   = config->GetEnthalpy_Formation();
    
    /*--- Calculate formation energy ---*/
    Ef = hf[val_Species] - Ru/Ms[val_Species] * Tref[val_Species];
    
    /*--- Electron t-r mode contributes to mixture vib-el energy ---*/
    Eel = (3.0/2.0) * Ru/Ms[val_Species] * (Tve - Tref[val_Species]) + Ef;
    Ev  = 0.0;
    
  }

  /*--- Heavy particle energy ---*/
  else {
    
    /*--- Read from CConfig ---*/
    thetav    = config->GetCharVibTemp();
    thetae    = config->GetCharElTemp();
    g         = config->GetElDegeneracy();
    nElStates = config->GetnElStates();
    
    /*--- Calculate vibrational energy (harmonic-oscillator model) ---*/
    if (thetav[val_Species] != 0.0)
      Ev = Ru/Ms[val_Species] * thetav[val_Species] / (exp(thetav[val_Species]/Tve)-1.0);
    else
      Ev = 0.0;
    
    /*--- Calculate electronic energy ---*/
    num = 0.0;
    denom = g[val_Species][0] * exp(thetae[val_Species][0]/Tve);
    for (iEl = 1; iEl < nElStates[val_Species]; iEl++) {
      num   += g[val_Species][iEl] * thetae[val_Species][iEl] * exp(-thetae[val_Species][iEl]/Tve);
      denom += g[val_Species][iEl] * exp(-thetae[val_Species][iEl]/Tve);
    }
    Eel = Ru/Ms[val_Species] * (num/denom);
  }
  
  return Ev + Eel;
}


su2double CTNE2EulerVariable::CalcHs(su2double *V, CConfig *config,
                                  unsigned short val_Species) {

  su2double Ru, *xi, *Ms, *hf, *Tref, T, eve, ef, hs;
  
  /*--- Read from config ---*/
  xi   = config->GetRotationModes();
  Ms   = config->GetMolar_Mass();
  hf   = config->GetEnthalpy_Formation();
  Tref = config->GetRefTemperature();
  
  /*--- Rename for convenience ---*/
  Ru = UNIVERSAL_GAS_CONSTANT;
  T = V[T_INDEX];
  
  /*--- Calculate vibrational-electronic energy per unit mass ---*/
  eve = CalcEve(V, config, val_Species);
  
  /*--- Calculate formation energy ---*/
  ef = hf[val_Species] - Ru/Ms[val_Species]*Tref[val_Species];
  
//  hs = Ru/Ms[val_Species]*T
//     + (3.0/2.0+xi[val_Species]/2.0)*Ru/Ms[val_Species]
//     + hf[val_Species] + eve;
  
  hs = (3.0/2.0+xi[val_Species]/2.0) * Ru/Ms[val_Species] * (T-Tref[val_Species])
   + eve + ef;
  
  return hs;
}

su2double CTNE2EulerVariable::CalcCvve(su2double val_Tve, CConfig *config, unsigned short val_Species) {
  
  unsigned short iEl, *nElStates;
  su2double *Ms, *thetav, **thetae, **g, Ru;
  su2double thoTve, exptv, thsqr, Cvvs, Cves;
  su2double Tve;
  su2double num, num2, num3, denom;
  
  /*--- Read from config ---*/
  thetav    = config->GetCharVibTemp();
  thetae    = config->GetCharElTemp();
  g         = config->GetElDegeneracy();
  Ms        = config->GetMolar_Mass();
  nElStates = config->GetnElStates();
  
  /*--- Rename for convenience ---*/
  Ru  = UNIVERSAL_GAS_CONSTANT;
  Tve = val_Tve;
  
  /*--- Initialize ---*/
  Cves = 0.0;
  Cvvs = 0.0;
  
  /*--- If requesting electron specific heat ---*/
  if (ionization && val_Species == nSpecies-1) {
    Cves = 3.0/2.0 * Ru/Ms[nSpecies-1];
    return Cves;
  }
  
  /*--- Heavy particle specific heat ---*/
  else {
    
    /*--- Vibrational energy ---*/
    if (thetav[val_Species] != 0.0) {
      thoTve = thetav[val_Species]/Tve;
      exptv = exp(thetav[val_Species]/Tve);
      thsqr = thetav[val_Species]*thetav[val_Species];
      Cvvs  = Ru/Ms[val_Species] * thoTve*thoTve * exptv / ((exptv-1.0)*(exptv-1.0));
    }
    
    /*--- Electronic energy ---*/
    if (nElStates[val_Species] != 0) {
      num = 0.0; num2 = 0.0;
      denom = g[val_Species][0] * exp(thetae[val_Species][0]/Tve);
      num3  = g[val_Species][0] * (thetae[val_Species][0]/(Tve*Tve))*exp(-thetae[val_Species][0]/Tve);
      for (iEl = 1; iEl < nElStates[val_Species]; iEl++) {
        thoTve = thetae[val_Species][iEl]/Tve;
        exptv = exp(-thetae[val_Species][iEl]/Tve);
        
        num   += g[val_Species][iEl] * thetae[val_Species][iEl] * exptv;
        denom += g[val_Species][iEl] * exptv;
        num2  += g[val_Species][iEl] * (thoTve*thoTve) * exptv;
        num3  += g[val_Species][iEl] * thoTve/Tve * exptv;
      }
      Cves = Ru/Ms[val_Species] * (num2/denom - num*num3/(denom*denom));
    }
  }
  return Cvvs + Cves;
}



void CTNE2EulerVariable::CalcdTdU(su2double *V, CConfig *config,
                                  su2double *val_dTdU) {
  
  unsigned short iDim, iSpecies, nHeavy, nEl;
  su2double *Ms, *xi, *Tref, *hf;
  su2double v2, ef, T, Cvtrs, rhoCvtr, Ru;
  
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
  
  /*--- Get gas properties from config settings ---*/
  Ms   = config->GetMolar_Mass();
  xi   = config->GetRotationModes();
  hf   = config->GetEnthalpy_Formation();
  Tref = config->GetRefTemperature();
  
  /*--- Rename for convenience ---*/
  rhoCvtr = V[RHOCVTR_INDEX];
  T       = V[T_INDEX];
  Ru      = UNIVERSAL_GAS_CONSTANT;
  
  /*--- Calculate supporting quantities ---*/
  v2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    v2 += V[VEL_INDEX+iDim]*V[VEL_INDEX+iDim];
  
  /*--- Species density derivatives ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    ef    = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];
    Cvtrs = (3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
    val_dTdU[iSpecies]   = (-ef + 0.5*v2 + Cvtrs*(Tref[iSpecies]-T)) / rhoCvtr;
  }
  if (ionization) {
    cout << "CTNE2Variable: NEED TO IMPLEMENT dTdU for IONIZED MIX" << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Momentum derivatives ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    val_dTdU[nSpecies+iDim] = -V[VEL_INDEX+iDim] / V[RHOCVTR_INDEX];
  
  /*--- Energy derivatives ---*/
  val_dTdU[nSpecies+nDim]   = 1.0 / V[RHOCVTR_INDEX];
  val_dTdU[nSpecies+nDim+1] = -1.0 / V[RHOCVTR_INDEX];
  
}

void CTNE2EulerVariable::CalcdTvedU(su2double *V, CConfig *config,
                                    su2double *val_dTvedU) {
  
  unsigned short iDim, iSpecies;
  su2double rhoCvve;
  su2double eve;
  
  /*--- Rename for convenience ---*/
  rhoCvve = V[RHOCVVE_INDEX];
  
  /*--- Species density derivatives ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    eve = CalcEve(V, config, iSpecies);
    val_dTvedU[iSpecies] = -eve/rhoCvve;
  }
  /*--- Momentum derivatives ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    val_dTvedU[nSpecies+iDim] = 0.0;
  
  /*--- Energy derivatives ---*/
  val_dTvedU[nSpecies+nDim]   = 0.0;
  val_dTvedU[nSpecies+nDim+1] = 1.0 / rhoCvve;
  
}

bool CTNE2EulerVariable::SetPrimVar_Compressible(CConfig *config) {

  bool V_err, bkup;
  unsigned short iVar;
  
  /*--- Convert conserved to primitive variables ---*/
  V_err = Cons2PrimVar(config, Solution, Primitive, dPdU, dTdU, dTvedU);
  if (V_err) {
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    bkup = Cons2PrimVar(config, Solution, Primitive, dPdU, dTdU, dTvedU);
  }
  
  SetVelocity2();
  
  return !V_err;
  
//	unsigned short iVar, iSpecies;
//  bool check_dens, check_press, check_sos, check_temp, RightVol;
//  
//  /*--- Initialize booleans that check for physical solutions ---*/
//  check_dens  = false;
//  check_press = false;
//  check_sos   = false;
//  check_temp  = false;
//  RightVol    = true;
//  
//  /*--- Calculate primitive variables ---*/
//  // Solution:  [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
//  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
//  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
//  SetDensity();                             // Compute species & mixture density
//	SetVelocity2();                           // Compute the square of the velocity (req. mixture density).
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    check_dens = ((Solution[iSpecies] < 0.0) || check_dens);  // Check the density
//  check_temp  = SetTemperature(config);     // Compute temperatures (T & Tve) (req. mixture density).
// 	check_press = SetPressure(config);        // Requires T & Tve computation.
//  CalcdPdU(Primitive, config, dPdU);        // Requires density, pressure, rhoCvtr, & rhoCvve.
//  CalcdTdU(Primitive, config, dTdU);
//  CalcdTvedU(Primitive, config, dTvedU);
//  check_sos   = SetSoundSpeed(config);      // Requires density, pressure, rhoCvtr, & rhoCvve.
//  
//  /*--- Check that the solution has a physical meaning ---*/
//  if (check_dens || check_press || check_sos || check_temp) {
//    
//    /*--- Copy the old solution ---*/
//    for (iVar = 0; iVar < nVar; iVar++)
//      Solution[iVar] = Solution_Old[iVar];
//    
//    /*--- Recompute the primitive variables ---*/
//    SetDensity();                           // Compute mixture density
//    SetVelocity2();                         // Compute square of the velocity (req. mixture density).
//    check_temp  = SetTemperature(config);   // Compute temperatures (T & Tve)
//    check_press = SetPressure(config);      // Requires T & Tve computation.
//    CalcdPdU(Primitive, config, dPdU);                       // Requires density, pressure, rhoCvtr, & rhoCvve.
//    CalcdTdU(Primitive, config, dTdU);
//    CalcdTvedU(Primitive, config, dTvedU);
//    check_sos   = SetSoundSpeed(config);    // Requires density & pressure computation.
//  
//    RightVol = false;
//  }
//  
//  SetEnthalpy();                            // Requires density & pressure computation.
//
//  
//  
//  return RightVol;
}


bool CTNE2EulerVariable::Cons2PrimVar(CConfig *config, su2double *U, su2double *V,
                                      su2double *val_dPdU, su2double *val_dTdU,
                                      su2double *val_dTvedU) {
  
  bool ionization, nonphys, nrconvg, converr;
	unsigned short iDim, iEl, iSpecies, nHeavy, nEl, iIter, maxBIter, maxNIter;
  su2double rho, rhoE, rhoEve, rhoE_f, rhoE_ref, rhoEve_t;
  su2double Ru, sqvel, rhoCvtr, rhoCvve;
  su2double Tve, Tve2, Tve_o;
  su2double f, df, tol;
  su2double Tmin, Tmax, Tvemin, Tvemax;
  su2double radical2;
  su2double *xi, *Ms, *hf, *Tref;
  
  /*--- Conserved & primitive vector layout ---*/
  // U:  [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
  // V: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  
  converr = false;
  
  /*--- Read from config ---*/
  xi         = config->GetRotationModes();      // Rotational modes of energy storage
  Ms         = config->GetMolar_Mass();         // Species molar mass
  Tref       = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf         = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]
  ionization = config->GetIonization();
  
  /*--- Rename for convenience ---*/
  Ru     = UNIVERSAL_GAS_CONSTANT;    // Universal gas constant [J/(kmol*K)]
  rhoE   = U[nSpecies+nDim];          // Density * energy [J/m3]
  rhoEve = U[nSpecies+nDim+1];        // Density * energy_ve [J/m3]

  /*--- Determine the number of heavy species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
  
  /*--- Species & mixture density ---*/
  V[RHO_INDEX] = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    V[RHOS_INDEX+iSpecies] = U[iSpecies];
    V[RHO_INDEX]          += U[iSpecies];
    if (V[RHOS_INDEX+iSpecies] < 0.0) {
      V[RHOS_INDEX+iSpecies] = 1E-14;
      U[iSpecies] = 1E-14;
      converr = true;
//      cout << "density < 0" << endl;
    }
  }
  
  /*--- Mixture velocity ---*/
  sqvel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    V[VEL_INDEX+iDim] = U[nSpecies+iDim]/V[RHO_INDEX];
    sqvel += V[VEL_INDEX+iDim]*V[VEL_INDEX+iDim];
  }
  
  
  /*--- T-R Temperature ---*/
  // Set temperature clipping values
  Tmin   = 80.0;
  Tmax   = 6E4;
  
  // Rename for convenience
  rho = V[RHO_INDEX];
  rhoE_f   = 0.0;
  rhoE_ref = 0.0;
  rhoCvtr  = 0.0;
  
  // Calculate mixture properties
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhoCvtr  += U[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
    rhoE_ref += U[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * Tref[iSpecies];
    rhoE_f   += U[iSpecies] * (hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies]);
  }
  
  // Calculate T-R temperature
  V[T_INDEX] = (rhoE - rhoEve - rhoE_f + rhoE_ref - 0.5*rho*sqvel) / rhoCvtr;
  V[RHOCVTR_INDEX] = rhoCvtr;
  
  // Error checking
  if (V[T_INDEX] < Tmin) {
    V[T_INDEX] = Tmin;
    converr = true;
//    cout << "T < Tmin" << endl;
  }
  if (V[T_INDEX] > Tmax) {
    V[T_INDEX] = Tmax;
    converr = true;
//    cout << "T > Tmax" << endl;
  }
  
  Tve = V[T_INDEX];
  V[TVE_INDEX] = V[T_INDEX];
  
  /*--- V-E Temperature ---*/
  // Set temperature clipping values
  Tvemin = 80.0;
  Tvemax = 4E4;
  
  // Set tolerance for Newton-Raphson method
  tol     = 1.0E-4;
  maxNIter = 20;
  maxBIter = 24;
  
  // Check for non-physical solutions
  nonphys = false;
  V[TVE_INDEX] = Tvemin;
  rhoEve_t = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhoEve_t += U[iSpecies]*CalcEve(V, config, iSpecies);
  if (rhoEve < rhoEve_t) {
    nonphys = true;
    converr = true;
    V[TVE_INDEX] = Tvemin;
//    cout << "Tve < Tve min" << endl;
  }
  
  V[TVE_INDEX] = Tvemax;
  rhoEve_t = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhoEve_t += U[iSpecies]*CalcEve(V, config, iSpecies);
  if (rhoEve > rhoEve_t) {
    nonphys = true;
    converr = true;
    V[TVE_INDEX] = Tvemax;
//    cout << "Tve > Tve max" << endl;
  }
  
  // Initialize trial values of Tve for Newton-Raphson method
//  Tve   = Primitive[TVE_INDEX];
//  Tve_o = Primitive[TVE_INDEX];
  Tve   = V[T_INDEX];
  Tve_o = V[T_INDEX];

  // Newton-Raphson
  if (!nonphys) {
    nrconvg = false;
    for (iIter = 0; iIter < maxNIter; iIter++) {
      rhoEve_t     = 0.0;
      rhoCvve      = 0.0;
      V[TVE_INDEX] = Tve;
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        rhoEve_t += U[iSpecies] * CalcEve(V, config, iSpecies);
        rhoCvve  += U[iSpecies] * CalcCvve(Tve, config, iSpecies);
      }
      
      // Find the root
      f  = rhoEve - rhoEve_t;
      df = -rhoCvve;
      Tve2 = Tve - (f/df)*0.5;
      
      // Check for non-physical conditions
      if ((Tve2 != Tve2) || (Tve2 < 0)) {
        nonphys = true;
        break;
      }
      // Check for convergence
      if (fabs(Tve2-Tve) < tol) {
        nrconvg = true;
        break;
      }
      // Update guess and iterate
      Tve = Tve2;
    }
    
    // If converged, load in to the solution
    if (nrconvg == true)
      V[TVE_INDEX] = Tve;
    
    else {
      
      // Bisection method
      Tve_o = Tvemin;
      Tve2  = Tvemax;
      
      for (iIter = 0; iIter < maxBIter; iIter++) {
        Tve = (Tve_o+Tve2)/2.0;
        V[TVE_INDEX] = Tve;
        rhoEve_t = 0.0;
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          rhoEve_t += U[iSpecies] * CalcEve(V, config, iSpecies);
        
        if (fabs(rhoEve_t - rhoEve) < tol) {
          V[TVE_INDEX] = Tve;
          break;
        } else {
          if (rhoEve_t > rhoEve) Tve2 = Tve;
          else                  Tve_o = Tve;
        }
      }
    }
  }

  V[RHOCVVE_INDEX] = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    V[RHOCVVE_INDEX] += U[iSpecies]*CalcCvve(V[TVE_INDEX], config, iSpecies);
  
  /*--- Pressure ---*/
  V[P_INDEX] = 0.0;
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    V[P_INDEX] += U[iSpecies] * Ru/Ms[iSpecies] * V[T_INDEX];
  for (iEl = 0; iEl < nEl; iEl++)
    V[P_INDEX] += U[nSpecies-1] * Ru/Ms[nSpecies-1] * V[TVE_INDEX];

  if (V[P_INDEX] < 0.0) {
//    cout << "P: " << V[P_INDEX] << endl;
    V[P_INDEX] = 1E-14;
    converr = true;
  }
  
  /*--- Partial derivatives of pressure and temperature ---*/
  CalcdPdU(  V, config, val_dPdU  );
  CalcdTdU(  V, config, val_dTdU  );
  CalcdTvedU(V, config, val_dTvedU);
  
  
  /*--- Sound speed ---*/
  radical2 = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    radical2 += V[RHOS_INDEX+iSpecies]/V[RHO_INDEX] * val_dPdU[iSpecies];
  for (iDim = 0; iDim < nDim; iDim++)
    radical2 += V[VEL_INDEX+iDim]*val_dPdU[nSpecies+iDim];
  radical2 += (U[nSpecies+nDim]+V[P_INDEX])/V[RHO_INDEX] * val_dPdU[nSpecies+nDim];
  radical2 += U[nSpecies+nDim+1]/V[RHO_INDEX] * val_dPdU[nSpecies+nDim+1];
  V[A_INDEX] = sqrt(radical2);
  
  if (radical2 < 0.0) {
//    unsigned short iVar;
    converr = true;
    V[A_INDEX] = 1E4;
//    cout << "Soundspeed NaN" << endl;
//    for (iVar = 0; iVar < nVar; iVar++)
//      cout << "dPdU[" << iVar << "]: " << val_dPdU[iVar] << endl;
  }
  
  
  /*--- Enthalpy ---*/
  V[H_INDEX] = (U[nSpecies+nDim] + V[P_INDEX])/V[RHO_INDEX];
  
//  unsigned short iVar;
//  for (iVar = 0; iVar < nPrimVar; iVar++)
//    cout << "V[" << iVar << "]: " << V[iVar] << endl;
//  for (iVar = 0; iVar < nVar; iVar++)
//    cout << "dPdU[" << iVar << "]: " << val_dPdU[iVar] << endl;
//  cin.get();
  
  return converr;
}

void CTNE2EulerVariable::Prim2ConsVar(CConfig *config, su2double *V, su2double *U) {
  unsigned short iDim, iEl, iSpecies, nEl, nHeavy;
  unsigned short *nElStates;
  su2double Ru, Tve, T, sqvel, rhoE, rhoEve, Ef, Ev, Ee, rhos, rhoCvtr, num, denom;
  su2double *thetav, *Ms, *xi, *hf, *Tref;
  su2double **thetae, **g;
  
  /*--- Determine the number of heavy species ---*/
  ionization = config->GetIonization();
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
  
  /*--- Load variables from the config class --*/
  xi        = config->GetRotationModes();      // Rotational modes of energy storage
  Ms        = config->GetMolar_Mass();         // Species molar mass
  thetav    = config->GetCharVibTemp();        // Species characteristic vib. temperature [K]
  thetae    = config->GetCharElTemp();         // Characteristic electron temperature [K]
  g         = config->GetElDegeneracy();       // Degeneracy of electron states
  nElStates = config->GetnElStates();          // Number of electron states
  Tref      = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf        = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]
  
  /*--- Rename & initialize for convenience ---*/
  Ru      = UNIVERSAL_GAS_CONSTANT;         // Universal gas constant [J/(kmol*K)]
  Tve     = V[TVE_INDEX];                   // Vibrational temperature [K]
  T       = V[T_INDEX];                     // Translational-rotational temperature [K]
  sqvel   = 0.0;                            // Velocity^2 [m2/s2]
  rhoE    = 0.0;                            // Mixture total energy per mass [J/kg]
  rhoEve  = 0.0;                            // Mixture vib-el energy per mass [J/kg]
  denom   = 0.0;
  rhoCvtr = 0.0;
  
  for (iDim = 0; iDim < nDim; iDim++)
    sqvel += V[VEL_INDEX+iDim]*V[VEL_INDEX+iDim];
  
  /*--- Set species density ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    U[iSpecies] = V[RHOS_INDEX+iSpecies];
  
  /*--- Set momentum ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    U[nSpecies+iDim] = V[RHO_INDEX]*V[VEL_INDEX+iDim];
  
  /*--- Set the total energy ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhos = U[iSpecies];
    
    // Species formation energy
    Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];
    
    // Species vibrational energy
    if (thetav[iSpecies] != 0.0)
      Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Tve)-1.0);
    else
      Ev = 0.0;
    
    // Species electronic energy
    num = 0.0;
    denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);
    for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
    }
    Ee = Ru/Ms[iSpecies] * (num/denom);
    
    // Mixture total energy
    rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (T-Tref[iSpecies])
                    + Ev + Ee + Ef + 0.5*sqvel);
    
    // Mixture vibrational-electronic energy
    rhoEve += rhos * (Ev + Ee);
  }
  for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
    // Species formation energy
    Ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1] * Tref[nSpecies-1];
    
    // Electron t-r mode contributes to mixture vib-el energy
    rhoEve += (3.0/2.0) * Ru/Ms[nSpecies-1] * (Tve - Tref[nSpecies-1]);
  }
  
  /*--- Set energies ---*/
  U[nSpecies+nDim]   = rhoE;
  U[nSpecies+nDim+1] = rhoEve;
  
//  unsigned short iVar;
//  cout << "In CVariable:" << endl;
//  cout << "U: " << endl;
//  for (iVar = 0; iVar < nVar; iVar++)
//    cout << U[iVar] << endl;
//  cout << "V: " << endl;
//  for (iVar = 0; iVar < nPrimVar; iVar++)
//    cout << V[iVar] << endl;
//  cin.get();
 
  return;
}

CTNE2NSVariable::CTNE2NSVariable(void) : CTNE2EulerVariable() { }


CTNE2NSVariable::CTNE2NSVariable(unsigned short val_nDim,
                                 unsigned short val_nvar,
                                 unsigned short val_nprimvar,
                                 unsigned short val_nprimvargrad,
                                 CConfig *config) : CTNE2EulerVariable(val_nDim,
                                                                       val_nvar,
                                                                       val_nprimvar,
                                                                       val_nprimvargrad,
                                                                       config) {
  
  Temperature_Ref = config->GetTemperature_Ref();
	Viscosity_Ref   = config->GetViscosity_Ref();
	Viscosity_Inf   = config->GetViscosity_FreeStreamND();
	Prandtl_Lam     = config->GetPrandtl_Lam();
  DiffusionCoeff  = new su2double[nSpecies];
}


CTNE2NSVariable::CTNE2NSVariable(su2double val_pressure, su2double *val_massfrac,
                                 su2double *val_mach, su2double val_temperature,
                                 su2double val_temperature_ve,
                                 unsigned short val_nDim,
                                 unsigned short val_nvar,
                                 unsigned short val_nvarprim,
                                 unsigned short val_nvarprimgrad,
                                 CConfig *config) : CTNE2EulerVariable(val_pressure,
                                                                       val_massfrac,
                                                                       val_mach,
                                                                       val_temperature,
                                                                       val_temperature_ve,
                                                                       val_nDim,
                                                                       val_nvar,
                                                                       val_nvarprim,
                                                                       val_nvarprimgrad,
                                                                       config) {
  
	Temperature_Ref = config->GetTemperature_Ref();
	Viscosity_Ref   = config->GetViscosity_Ref();
	Viscosity_Inf   = config->GetViscosity_FreeStreamND();
	Prandtl_Lam     = config->GetPrandtl_Lam();
  DiffusionCoeff  = new su2double[nSpecies];
}

CTNE2NSVariable::CTNE2NSVariable(su2double *val_solution, unsigned short val_nDim,
                                 unsigned short val_nvar,
                                 unsigned short val_nprimvar,
                                 unsigned short val_nprimvargrad,
                                 CConfig *config) : CTNE2EulerVariable(val_solution,
                                                                       val_nDim,
                                                                       val_nvar,
                                                                       val_nprimvar,
                                                                       val_nprimvargrad,
                                                                       config) {  
	Temperature_Ref = config->GetTemperature_Ref();
	Viscosity_Ref   = config->GetViscosity_Ref();
	Viscosity_Inf   = config->GetViscosity_FreeStreamND();
	Prandtl_Lam     = config->GetPrandtl_Lam();
  DiffusionCoeff  = new su2double[nSpecies];
}

CTNE2NSVariable::~CTNE2NSVariable(void) {
  delete [] DiffusionCoeff;
}

void CTNE2NSVariable::SetDiffusionCoeff(CConfig *config) {
  unsigned short iSpecies, jSpecies, nHeavy, nEl;
  su2double rho, T, Tve, P;
  su2double *Ms, Mi, Mj, pi, Ru, kb, gam_i, gam_j, gam_t, Theta_v;
  su2double denom, d1_ij, D_ij;
  su2double ***Omega00, Omega_ij;
  
  /*--- Acquire gas parameters from CConfig ---*/
  Omega00 = config->GetCollisionIntegral00();
  Ms      = config->GetMolar_Mass();
  if (ionization) {nHeavy = nSpecies-1;  nEl = 1;}
  else            {nHeavy = nSpecies;    nEl = 0;}
  
  /*--- Rename for convenience ---*/
  rho  = Primitive[RHO_INDEX];
  T    = Primitive[T_INDEX];
  Tve  = Primitive[TVE_INDEX];
  P    = Primitive[P_INDEX];
  pi   = PI_NUMBER;
  Ru   = UNIVERSAL_GAS_CONSTANT;
  kb   = BOLTZMANN_CONSTANT;
  
  /*--- Calculate mixture gas constant ---*/
  gam_t = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    gam_t += Primitive[RHOS_INDEX+iSpecies] / (rho*Ms[iSpecies]);
  }
  
  /*--- Mixture thermal conductivity via Gupta-Yos approximation ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    
    /*--- Initialize the species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = 0.0;
    
    /*--- Calculate molar concentration ---*/
    Mi      = Ms[iSpecies];
    gam_i   = Primitive[RHOS_INDEX+iSpecies] / (rho*Mi);
    Theta_v = config->GetCharVibTemp(iSpecies);
    
    denom = 0.0;
    for (jSpecies = 0; jSpecies < nHeavy; jSpecies++) {
      if (jSpecies != iSpecies) {
        Mj    = Ms[jSpecies];
        gam_j = Primitive[RHOS_INDEX+iSpecies] / (rho*Mj);
        
        /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
        Omega_ij = 1E-20 * Omega00[iSpecies][jSpecies][3]
        * pow(T, Omega00[iSpecies][jSpecies][0]*log(T)*log(T)
              + Omega00[iSpecies][jSpecies][1]*log(T)
              + Omega00[iSpecies][jSpecies][2]);
        
        /*--- Calculate "delta1_ij" ---*/
        d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;
        
        /*--- Calculate heavy-particle binary diffusion coefficient ---*/
        D_ij = kb*T/(P*d1_ij);
        denom += gam_j/D_ij;
      }
    }
    
    if (ionization) {
      jSpecies = nSpecies-1;
      Mj       = config->GetMolar_Mass(jSpecies);
      gam_j    = Primitive[RHOS_INDEX+iSpecies] / (rho*Mj);
      
      /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
      Omega_ij = 1E-20 * Omega00[iSpecies][jSpecies][3]
      * pow(Tve, Omega00[iSpecies][jSpecies][0]*log(Tve)*log(Tve)
            + Omega00[iSpecies][jSpecies][1]*log(Tve)
            + Omega00[iSpecies][jSpecies][2]);
      
      /*--- Calculate "delta1_ij" ---*/
      d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;
    }
    
    /*--- Assign species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = gam_t*gam_t*Mi*(1-Mi*gam_i) / denom;
  }
  if (ionization) {
    iSpecies = nSpecies-1;
    /*--- Initialize the species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = 0.0;
    
    /*--- Calculate molar concentration ---*/
    Mi      = Ms[iSpecies];
    gam_i   = Primitive[RHOS_INDEX+iSpecies] / (rho*Mi);
    
    denom = 0.0;
    for (jSpecies = 0; jSpecies < nHeavy; jSpecies++) {
      if (iSpecies != jSpecies) {
        Mj    = config->GetMolar_Mass(jSpecies);
        gam_j = Primitive[RHOS_INDEX+iSpecies] / (rho*Mj);
        
        /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
        Omega_ij = 1E-20 * Omega00[iSpecies][jSpecies][3]
        * pow(Tve, Omega00[iSpecies][jSpecies][0]*log(Tve)*log(Tve)
              + Omega00[iSpecies][jSpecies][1]*log(Tve)
              + Omega00[iSpecies][jSpecies][2]);
        
        /*--- Calculate "delta1_ij" ---*/
        d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;
        
        /*--- Calculate heavy-particle binary diffusion coefficient ---*/
        D_ij = kb*Tve/(P*d1_ij);
        denom += gam_j/D_ij;
      }
    }
    DiffusionCoeff[iSpecies] = gam_t*gam_t*Ms[iSpecies]*(1-Ms[iSpecies]*gam_i)
                             / denom;
  }
}


void CTNE2NSVariable::SetLaminarViscosity(CConfig *config) {
  unsigned short iSpecies, jSpecies, nHeavy, nEl;
	su2double rho, T, Tve;
  su2double *Ms, Mi, Mj, pi, Ru, Na, gam_i, gam_j, denom;
  su2double ***Omega11, Omega_ij, d2_ij;
  
  /*--- Acquire gas parameters from CConfig ---*/
  Omega11 = config->GetCollisionIntegral11();
  Ms      = config->GetMolar_Mass();
  if (ionization) {nHeavy = nSpecies-1;  nEl = 1;}
  else            {nHeavy = nSpecies;    nEl = 0;}
  
  /*--- Rename for convenience ---*/
  rho = Primitive[RHO_INDEX];
  T   = Primitive[T_INDEX];
  Tve = Primitive[TVE_INDEX];
  pi  = PI_NUMBER;
  Ru  = UNIVERSAL_GAS_CONSTANT;
  Na  = AVOGAD_CONSTANT;
  
  LaminarViscosity = 0.0;
  
  /*--- Mixture viscosity via Gupta-Yos approximation ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    denom = 0.0;
    
    /*--- Calculate molar concentration ---*/
    Mi    = Ms[iSpecies];
    gam_i = Primitive[RHOS_INDEX+iSpecies] / (rho*Mi);
    for (jSpecies = 0; jSpecies < nHeavy; jSpecies++) {
      Mj    = Ms[jSpecies];
      gam_j = Primitive[RHOS_INDEX+jSpecies] / (rho*Mj);
      
      /*--- Calculate "delta" quantities ---*/
      Omega_ij = 1E-20 * Omega11[iSpecies][jSpecies][3]
                * pow(T, Omega11[iSpecies][jSpecies][0]*log(T)*log(T)
                      + Omega11[iSpecies][jSpecies][1]*log(T)
                      + Omega11[iSpecies][jSpecies][2]);
      d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;
      
      /*--- Add to denominator of viscosity ---*/
      denom += gam_j*d2_ij;
    }
    if (ionization) {
      jSpecies = nSpecies-1;
      Mj    = Ms[jSpecies];
      gam_j = Primitive[RHOS_INDEX+jSpecies] / (rho*Mj);
      
      /*--- Calculate "delta" quantities ---*/
      Omega_ij = 1E-20 * Omega11[iSpecies][jSpecies][3]
      * pow(Tve, Omega11[iSpecies][jSpecies][0]*log(Tve)*log(Tve)
            + Omega11[iSpecies][jSpecies][1]*log(Tve)
            + Omega11[iSpecies][jSpecies][2]);
      d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;

      denom += gam_j*d2_ij;
    }
    
    /*--- Calculate species laminar viscosity ---*/
    LaminarViscosity += (Mi/Na * gam_i) / denom;
  }
  if (ionization) {
    iSpecies = nSpecies-1;
    denom = 0.0;
    
    /*--- Calculate molar concentration ---*/
    Mi    = Ms[iSpecies];
    gam_i = Primitive[RHOS_INDEX+iSpecies] / (rho*Mi);
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      Mj    = Ms[jSpecies];
      gam_j = Primitive[RHOS_INDEX+jSpecies] / (rho*Mj);
      
      /*--- Calculate "delta" quantities ---*/
      Omega_ij = 1E-20 * Omega11[iSpecies][jSpecies][3]
      * pow(Tve, Omega11[iSpecies][jSpecies][0]*log(Tve)*log(Tve)
            + Omega11[iSpecies][jSpecies][1]*log(Tve)
            + Omega11[iSpecies][jSpecies][2]);
      d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;
      
      /*--- Add to denominator of viscosity ---*/
      denom += gam_j*d2_ij;
    }
    LaminarViscosity += (Mi/Na * gam_i) / denom;
  }
}


void CTNE2NSVariable ::SetThermalConductivity(CConfig *config) {
	unsigned short iSpecies, jSpecies, nHeavy, nEl;
  su2double rho, T, Tve, Cvve;
  su2double *xi, *Ms, Mi, Mj, mi, mj, pi, R, Ru, Na, kb, gam_i, gam_j, Theta_v;
  su2double denom_t, denom_r, d1_ij, d2_ij, a_ij;
  su2double ***Omega00, ***Omega11, Omega_ij;
  
  if (ionization) {
    cout << "SetThermalConductivity: NEEDS REVISION w/ IONIZATION" << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Acquire gas parameters from CConfig ---*/
  Omega00 = config->GetCollisionIntegral00();
  Omega11 = config->GetCollisionIntegral11();
  Ms      = config->GetMolar_Mass();
  xi      = config->GetRotationModes();
  if (ionization) {nHeavy = nSpecies-1;  nEl = 1;}
  else            {nHeavy = nSpecies;    nEl = 0;}
  
  /*--- Rename for convenience ---*/
  rho  = Primitive[RHO_INDEX];
  T    = Primitive[T_INDEX];
  Tve  = Primitive[TVE_INDEX];
  Cvve = Primitive[RHOCVVE_INDEX]/rho;
  pi   = PI_NUMBER;
  Ru   = UNIVERSAL_GAS_CONSTANT;
  Na   = AVOGAD_CONSTANT;
  kb   = BOLTZMANN_CONSTANT;
  
  /*--- Calculate mixture gas constant ---*/
  R = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    R += Ru * Primitive[RHOS_INDEX+iSpecies]/rho;
  }
  
  /*--- Mixture thermal conductivity via Gupta-Yos approximation ---*/
  ThermalCond    = 0.0;
  ThermalCond_ve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    
    /*--- Calculate molar concentration ---*/
    Mi      = Ms[iSpecies];
    mi      = Mi/Na;
    gam_i   = Primitive[RHOS_INDEX+iSpecies] / (rho*Mi);
    Theta_v = config->GetCharVibTemp(iSpecies);
    
    denom_t = 0.0;
    denom_r = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      Mj    = config->GetMolar_Mass(jSpecies);
      mj    = Mj/Na;
      gam_j = Primitive[RHOS_INDEX+iSpecies] / (rho*Mj);
      
      a_ij = 1.0 + (1.0 - mi/mj)*(0.45 - 2.54*mi/mj) / ((1.0 + mi/mj)*(1.0 + mi/mj));

      /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
      Omega_ij = 1E-20 * Omega00[iSpecies][jSpecies][3]
                * pow(T, Omega00[iSpecies][jSpecies][0]*log(T)*log(T)
                      + Omega00[iSpecies][jSpecies][1]*log(T)
                      + Omega00[iSpecies][jSpecies][2]);

      /*--- Calculate "delta1_ij" ---*/
      d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;
      
      /*--- Calculate the Omega^(1,1)_ij collision cross section ---*/
      Omega_ij = 1E-20 * Omega11[iSpecies][jSpecies][3]
                * pow(T, Omega11[iSpecies][jSpecies][0]*log(T)*log(T)
                      + Omega11[iSpecies][jSpecies][1]*log(T)
                      + Omega11[iSpecies][jSpecies][2]);
      
      /*--- Calculate "delta2_ij" ---*/
      d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;
      
      denom_t += a_ij*gam_j*d2_ij;
      denom_r += gam_j*d1_ij;
    }
    
    /*--- Translational contribution to thermal conductivity ---*/
    ThermalCond    += (15.0/4.0)*kb*gam_i/denom_t;
    
    /*--- Translational contribution to thermal conductivity ---*/
    if (xi[iSpecies] != 0.0)
      ThermalCond  += kb*gam_i/denom_r;
    
    /*--- Vibrational-electronic contribution to thermal conductivity ---*/
    ThermalCond_ve += kb*Cvve/R*gam_i / denom_r;
  }
}


bool CTNE2NSVariable::SetVorticity(bool val_limiter) {
	su2double u_y = Gradient_Primitive[VEL_INDEX][1];
	su2double v_x = Gradient_Primitive[VEL_INDEX+1][0];
	su2double u_z = 0.0;
	su2double v_z = 0.0;
	su2double w_x = 0.0;
	su2double w_y = 0.0;
  
	if (nDim == 3) {
		u_z = Gradient_Primitive[VEL_INDEX][2];
		v_z = Gradient_Primitive[VEL_INDEX+1][2];
		w_x = Gradient_Primitive[VEL_INDEX+2][0];
		w_y = Gradient_Primitive[VEL_INDEX+2][1];
	}
  
	Vorticity[0] = w_y-v_z;
	Vorticity[1] = -(w_x-u_z);
	Vorticity[2] = v_x-u_y;
  
  return false;
}

bool CTNE2NSVariable::SetPrimVar_Compressible(CConfig *config) {
//	unsigned short iVar, iSpecies;
//  bool check_dens, check_press, check_sos, check_temp, RightVol;
//  
//  /*--- Initialize booleans that check for physical solutions ---*/
//  check_dens  = false;
//  check_press = false;
//  check_sos   = false;
//  check_temp  = false;
//  RightVol    = true;
//  
//  /*--- Calculate primitive variables ---*/
//  // Solution:  [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
//  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
//  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
//  SetDensity();                             // Compute mixture density
//	SetVelocity2();                           // Compute the modulus of the velocity (req. mixture density).
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    check_dens = ((Solution[iSpecies] < 0.0) || check_dens);  // Check the density
//  check_temp  = SetTemperature(config);     // Compute temperatures (T & Tve)
// 	check_press = SetPressure(config);        // Requires T & Tve computation.
//  CalcdPdU(Primitive, config, dPdU);        // Requires density, pressure, rhoCvtr, & rhoCvve.
//  CalcdTdU(Primitive, config, dTdU);
//  CalcdTvedU(Primitive, config, dTvedU);
//  check_sos   = SetSoundSpeed(config);      // Requires density & pressure computation.
//  
//  /*--- Check that the solution has a physical meaning ---*/
//  if (check_dens || check_press || check_sos || check_temp) {
//    
//    /*--- Copy the old solution ---*/
//    for (iVar = 0; iVar < nVar; iVar++)
//      Solution[iVar] = Solution_Old[iVar];
//    
//    /*--- Recompute the primitive variables ---*/
//    SetDensity();                           // Compute mixture density
//    SetVelocity2();                         // Compute square of the velocity (req. mixture density).
//    check_temp  = SetTemperature(config);   // Compute temperatures (T & Tve)
//    check_press = SetPressure(config);      // Requires T & Tve computation.
//    CalcdPdU(Primitive, config, dPdU);      // Requires density, pressure, rhoCvtr, & rhoCvve.
//    CalcdTdU(Primitive, config, dTdU);
//    CalcdTvedU(Primitive, config, dTvedU);
//    check_sos   = SetSoundSpeed(config);    // Requires density & pressure computation.
//    
//    RightVol = false;
//  }

  bool V_err, bkup;
  unsigned short iVar;
  
  V_err = Cons2PrimVar(config, Solution, Primitive, dPdU, dTdU, dTvedU);
  if (V_err) {
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    bkup = Cons2PrimVar(config, Solution, Primitive, dPdU, dTdU, dTvedU);
  }
  
  SetVelocity2();
  
//  SetEnthalpy();                            // Requires density & pressure computation.
  SetDiffusionCoeff(config);
  SetLaminarViscosity(config);              // Requires temperature computation.
  SetThermalConductivity(config);
  
//  return RightVol;
  return !V_err;
}
