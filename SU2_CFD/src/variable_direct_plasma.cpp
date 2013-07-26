/*!
 * \file variable_direct_plasma.cpp
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

void CPlasmaVariable::SetPressure(CConfig *config) {
  
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