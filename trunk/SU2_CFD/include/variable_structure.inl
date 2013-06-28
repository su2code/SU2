/*!
 * \file variable_structure.inl
 * \brief In-Line subroutines of the <i>variable_structure.hpp</i> file.
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

#pragma once

inline void CVariable::SetVelocityInc2(void) { }

inline void CVariable::SetBetaInc2(CConfig *config) { }

inline void CVariable::SetPressureInc(void) { }

inline void CVariable::SetPressure(double val_pressure) { }

inline void CVariable::SetPressureInc(double val_pressure) { }

inline double CVariable::GetBetaInc2(void) { return 0; }

inline double CVariable::GetDiffLevelSet(void) { return 0; }

inline double CVariable::GetDensityInc(void) { return 0; }

inline double CVariable::GetProjVelInc(double *val_vector) { return 0; }

inline void CVariable::SetSolution(unsigned short val_var, double val_solution) { Solution[val_var] = val_solution; }

inline void CVariable::SetUndivided_Laplacian(unsigned short val_var, double val_undivided_laplacian) { Undivided_Laplacian[val_var] = val_undivided_laplacian; }

inline void CVariable::SetResidual(unsigned short val_var, double val_residual) { Residual[val_var] = val_residual; }

inline void CVariable::SetAuxVar(double val_auxvar) { AuxVar = val_auxvar; }

inline void CVariable::SetSolution_Old(unsigned short val_var, double val_solution_old) { Solution_Old[val_var] = val_solution_old; }

inline void CVariable::SetLimiter(unsigned short val_var, double val_limiter) { Limiter[val_var] = val_limiter; }

inline void CVariable::SetAuxVarGradient(unsigned short iDim, double val_gradient) { Grad_AuxVar[iDim] = val_gradient; }

inline double *CVariable::GetSolution(void) { return Solution; }

inline double *CVariable::GetSolution_Old(void) { return Solution_Old; }

inline double *CVariable::GetSolution_time_n(void) { return Solution_time_n; }

inline double *CVariable::GetSolution_time_n1(void) { return Solution_time_n1; }

inline double CVariable::GetAuxVar(void) { return AuxVar; }

inline double *CVariable::GetUnd_Lapl(void) { return Undivided_Laplacian; }

inline double CVariable::GetUnd_Lapl(unsigned short val_var) { return Undivided_Laplacian[val_var]; }

inline double CVariable::GetResidual(unsigned short val_var) { return Residual[val_var]; }

inline void CVariable::AddSolution(unsigned short val_var, double val_solution) {Solution[val_var] = Solution_Old[val_var] + val_solution; }

inline double CVariable::GetSolution(unsigned short val_var) { return Solution[val_var]; }

inline double CVariable::GetSolution_Old(unsigned short val_var) { return Solution_Old[val_var]; }

inline double *CVariable::GetResidual(void) { return Residual; }

inline double *CVariable::GetResConv(void) { return Res_Conv; }

inline double *CVariable::GetResVisc(void) { return Res_Visc; }

inline double *CVariable::GetResSour(void) { return Res_Sour; }

inline double CVariable::GetRes_Visc_RK(unsigned short val_var, unsigned short iRKStep) { return Res_Visc_RK[val_var][iRKStep]; }

inline double *CVariable::GetResidual_Sum(void) { return Residual_Sum; }

inline double *CVariable::GetResidual_Old(void) { return Residual_Old; }

inline void CVariable::SetGradient(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient[val_var][val_dim] = val_value; }

inline void CVariable::AddGradient(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient[val_var][val_dim] += val_value; }

inline void CVariable::SubtractGradient(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient[val_var][val_dim] -= val_value; }

inline void CVariable::AddAuxVarGradient(unsigned short val_dim, double val_value) { Grad_AuxVar[val_dim] += val_value; }

inline void CVariable::SubtractAuxVarGradient(unsigned short val_dim, double val_value) { Grad_AuxVar[val_dim] -= val_value; }

inline double CVariable::GetGradient(unsigned short val_var, unsigned short val_dim) { return Gradient[val_var][val_dim]; }

inline double CVariable::GetLimiter(unsigned short val_var) { return Limiter[val_var]; }

inline double CVariable::GetPreconditioner_Beta() { return 0; }

inline void CVariable::SetPreconditioner_Beta( double val_Beta) { }

inline double **CVariable::GetGradient(void) { return Gradient; }

inline double *CVariable::GetLimiter(void) { return Limiter; }

inline double *CVariable::GetAuxVarGradient(void) { return Grad_AuxVar; }

inline double CVariable::GetAuxVarGradient(unsigned short val_dim) { return Grad_AuxVar[val_dim]; }

inline double *CVariable::GetRes_TruncError(void) { return Res_TruncError; }

inline void CVariable::SetDelta_Time(double val_delta_time) { Delta_Time = val_delta_time; }

inline void CVariable::SetDelta_Time(double val_delta_time, unsigned short iSpecies) {  }

inline double CVariable::GetDelta_Time(void) { return Delta_Time; }

inline double CVariable::GetDelta_Time(unsigned short iSpecies) { return 0;}

inline void CVariable::SetMax_Lambda(double val_max_lambda) { Max_Lambda = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Inv(double val_max_lambda) { Max_Lambda_Inv = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Inv(double val_max_lambda, unsigned short val_species) { }

inline void CVariable::SetMax_Lambda_Visc(double val_max_lambda) { Max_Lambda_Visc = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Visc(double val_max_lambda, unsigned short val_species) { }

inline void CVariable::SetLambda(double val_lambda) { Lambda = val_lambda; }

inline void CVariable::SetLambda(double val_lambda, unsigned short iSpecies) {}

inline void CVariable::AddMax_Lambda(double val_max_lambda) { Max_Lambda += val_max_lambda; }

inline void CVariable::AddMax_Lambda_Inv(double val_max_lambda) { Max_Lambda_Inv += val_max_lambda; }

inline void CVariable::AddMax_Lambda_Visc(double val_max_lambda) { Max_Lambda_Visc += val_max_lambda; }

inline void CVariable::AddLambda(double val_lambda) { Lambda += val_lambda; }

inline void CVariable::AddLambda(double val_lambda, unsigned short iSpecies) {}

inline double CVariable::GetMax_Lambda(void) { return Max_Lambda; }

inline double CVariable::GetMax_Lambda_Inv(void) { return Max_Lambda_Inv; }

inline double CVariable::GetMax_Lambda_Visc(void) { return Max_Lambda_Visc; }

inline double CVariable::GetLambda(void) { return Lambda; }

inline double CVariable::GetLambda(unsigned short iSpecies) { return 0; }

inline double CVariable::GetSensor(void) { return Sensor; }

inline double CVariable::GetSensor(unsigned short iSpecies) { return 0;}

inline double CVariable::GetMax_Lambda_Inv(unsigned short iFluids) { return 0; }

inline double CVariable::GetMax_Lambda_Visc(unsigned short iFluids) { return 0; }

inline void CVariable::AddMax_Lambda_Inv(double val_max_lambda, unsigned short iSpecies) { }

inline void CVariable::AddMax_Lambda_Visc(double val_max_lambda, unsigned short iSpecies) { }

inline void CVariable::SetSensor(double val_sensor) { Sensor = val_sensor; }

inline void CVariable::SetSensor(double val_sensor, unsigned short val_iSpecies) {}

inline void CVariable::AddResidual(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual[iVar] += val_residual[iVar]; }

inline void CVariable::AddResidual(double val_residual) { Residual[0] += val_residual ; }

inline void CVariable::AddRes_Conv(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar] += val_residual[iVar]; }

inline void CVariable::AddRes_Visc(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar] += val_residual[iVar]; }

inline void CVariable::AddRes_Sour(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar] += val_residual[iVar]; }

inline void CVariable::SubtractResidual(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual[iVar] -= val_residual[iVar]; }

inline void CVariable::SubtractRes_Visc(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar] -= val_residual[iVar]; }

inline void CVariable::SubtractRes_Conv(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar] -= val_residual[iVar]; }

inline void CVariable::SubtractRes_Sour(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar] -= val_residual[iVar]; }

inline double CVariable::GetDensity(void) {	return 0; }

inline double CVariable::GetDensity(unsigned short val_iSpecies) {	return 0; }

inline double CVariable::GetEnergy(void) { return 0; }

inline double *CVariable::GetForceProj_Vector(void) { return NULL; }

inline double *CVariable::GetObjFuncSource(void) { return NULL; }

inline double *CVariable::GetIntBoundary_Jump(void) { return NULL; }

inline double CVariable::GetEddyViscosity(void) { return 0; }

inline double CVariable::GetEnthalpy(void) { return 0; }

inline double CVariable::GetEnthalpy(unsigned short iSpecies) { return 0; }

inline double CVariable::GetPressure(void) { return 0; }

inline double CVariable::GetPressure_Old(void) { return 0; }

inline double CVariable::GetDeltaPressure(void) { return 0; }

inline double CVariable::GetProjVel(double *val_vector) { return 0; }

inline double CVariable::GetProjVel(double *val_vector, unsigned short val_species) { return 0; }

inline double CVariable::GetSoundSpeed(void) { return 0; }

inline double CVariable::GetSoundSpeed(unsigned short val_var) { return 0; }

inline double CVariable::GetTemperature(void) { return 0; }

inline double CVariable::GetTemperature(unsigned short val_iSpecies) { return 0; }

inline double CVariable::GetTemperature_TR(unsigned short val_iSpecies) { return 0; }

inline double CVariable::GetTemperature_vib(unsigned short val_iSpecies) { return 0; }

inline double CVariable::GetTheta(void) { return 0; }

inline double CVariable::GetVelocity(unsigned short val_dim, bool val_incomp) { return 0; }

inline double CVariable::GetVelocity2(void) { return 0; }

inline double CVariable::GetPressure(unsigned short val_species) { return 0;}

inline double CVariable::GetVelocity2(unsigned short val_species) { return 0;}

inline double CVariable::GetVelocity(unsigned short val_dim, unsigned short val_species) { return 0;}

inline double CVariable::GetLaminarViscosity(void) { return 0; }

inline double CVariable::GetLaminarViscosityInc(void) { return 0; }

inline double CVariable::GetLaminarViscosity(unsigned short iSpecies) { return 0; }

inline double CVariable::GetEddyViscosity(unsigned short iSpecies) { return 0; }

inline double CVariable::GetVorticity(unsigned short val_dim) { return 0; }

inline double CVariable::GetStrainMag(void) { return 0; }

inline void CVariable::SetForceProj_Vector(double *val_ForceProj_Vector) { }

inline void CVariable::SetObjFuncSource(double *val_ObjFuncSource) { }

inline void CVariable::SetIntBoundary_Jump(double *val_IntBoundary_Jump) { }

inline void CVariable::SetEddyViscosity(unsigned short val_Kind_Turb_Model, CVariable *Turb_Solution) { }

inline void CVariable::SetEnthalpy(void) { }

inline void CVariable::SetBetaInc2(double val_betainc2) { }

inline void CVariable::SetDensityInc(double val_densityinc) { }

inline void CVariable::SetPhi_Old(double *val_phi) { }

inline void CVariable::SetDiffLevelSet(double val_difflevelset) { }

inline void CVariable::SetPressure(double Gamma, double *Coord) { }

inline void CVariable::SetPressure(CConfig *config, double *Coord) { }

inline void CVariable::SetPressure(double Gamma, double *Coord, double turb_ke) { }

inline void CVariable::SetPressure() { }

inline void CVariable::SetDeltaPressure(double *val_velocity, double Gamma) { }

inline void CVariable::SetSoundSpeed(CConfig *config) { }

inline void CVariable::SetSoundSpeed() { }

inline void CVariable::SetSoundSpeed(double Gamma) { }

inline void CVariable::SetTemperature(double Gas_Constant) { }

inline void CVariable::SetTemperature_TR(CConfig *config, double *Coord) { }

inline void CVariable::SetTemperature_vib(CConfig *config, double *Coord) { }

inline void CVariable::SetWallTemperature(double Temperature_Wall) { }

inline void CVariable::SetWallTemperature(double* Temperature_Wall) { }

inline void CVariable::SetThermalCoeff(CConfig *config) { }

inline void CVariable::SetTheta(double val_density, double *val_velocity, double val_enthalpy) { }

inline void CVariable::SetVelocity(double *val_velocity, bool val_incomp) { }

inline void CVariable::SetVelocity2(void) { }

inline void CVariable::SetVelocity_Old(double *val_velocity, bool val_incomp) { }

inline void CVariable::SetVelocity_Old(double *val_velocity, unsigned short iSpecies) { }

inline void CVariable::SetLaminarViscosity(CConfig *config) { }

inline void CVariable::SetLaminarViscosity(double val_laminar_viscosity) { }

inline void CVariable::SetLaminarViscosity(double val_laminar_viscosity, unsigned short iSpecies) { }

inline void CVariable::SetLaminarViscosityInc(double val_laminar_viscosity_inc) { }

inline void CVariable::SetEddyViscosity(double val_eddy_viscosity) { }

inline void CVariable::SetVorticity(void) { }

inline void CVariable::SetStrainMag(void) { }

inline void CVariable::SetGradient_PrimitiveZero(void) { }

inline void CVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { }

inline void CVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { }

inline double CVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return 0; }

inline void CVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { }

inline double **CVariable::GetGradient_Primitive(void) { return NULL; }

inline void CVariable::SetBlendingFunc(double val_viscosity, double val_dist, double val_density) { }

inline double CVariable::GetF1blending(void) { return 0; }

inline double CVariable::GetF2blending(void) { return 0; }

inline double CVariable::GetmuT() { return 0;}

inline void CVariable::SetmuT(double val_muT) { }

inline void CVariable::SetThickness_Noise(double val_thickness_noise) { }

inline void CVariable::SetLoading_Noise(double val_loading_noise) { }

inline void CVariable::SetQuadrupole_Noise(double val_quadrupole_noise) { }

inline double CVariable::GetThickness_Noise() { return 0;}

inline double CVariable::GetLoading_Noise() { return 0;}

inline double CVariable::GetQuadrupole_Noise() { return 0;}

inline double* CVariable::GetSolution_Direct() { return NULL; }

inline void CVariable::SetSolution_Direct(double *val_solution_direct) { }

inline void CVariable::SetChargeDensity(double positive_charge, double negative_charge) { }

inline void CVariable::SetPlasmaRhoUGradient(unsigned short iSpecies, double val_gradient, unsigned short iDim) { }

inline void CVariable::SetPlasmaTimeStep(double dt) { }

inline double* CVariable::GetChargeDensity() { return 0;}

inline double** CVariable::GetPlasmaRhoUGradient() { return 0;}

inline double CVariable::GetPlasmaTimeStep() { return 0;}

inline void CVariable::SetElectricField(double* ElectricField) { }

inline double* CVariable::GetElectricField() { return 0;}

inline void CVariable::SetTimeSpectral_Source(unsigned short val_var, double val_source) { }

inline double CVariable::GetTimeSpectral_Source(unsigned short val_var) { return 0; }

inline double CEulerVariable::GetDensity(void) { return Solution[0]; }

inline double CEulerVariable::GetDensityInc(void) { return Solution[nDim+1]; }

inline double CEulerVariable::GetBetaInc2(void) { return Solution[nDim+2]; }

inline double CEulerVariable::GetEnergy(void) { return Solution[nVar-1]/Solution[0]; };

inline double CEulerVariable::GetEnthalpy(void) { return Enthalpy; }

inline double CEulerVariable::GetPressure(void) { return Pressure; }

inline double CEulerVariable::GetPressure_Old(void) { return Pressure_Old; }

inline double CEulerVariable::GetSoundSpeed(void) { return SoundSpeed; }

inline double CEulerVariable::GetTemperature(void) { return Temperature; }

inline double CEulerVariable::GetVelocity(unsigned short val_dim, bool val_incomp) {
double velocity;
if (val_incomp) velocity = Solution[val_dim+1]/Solution[nDim+1];
else velocity = Solution[val_dim+1]/Solution[0]; 
return velocity;
}

inline double CEulerVariable::GetVelocity2(void) { return Velocity2; }

inline void CEulerVariable::SetEnthalpy(void) { Enthalpy = (Solution[nVar-1] + Pressure) / Solution[0]; }

inline void CEulerVariable::SetDensityInc(double val_density) { Solution[nDim+1] = val_density; }

inline void CEulerVariable::SetBetaInc2(double val_betainc2) { Solution[nDim+2] = val_betainc2; }

inline void CEulerVariable::SetPressure(double Gamma, double *Coord) {
	Pressure_Old = Pressure;
	Pressure = (Gamma-1.0)*Solution[0]*(Solution[nVar-1]/Solution[0]-0.5*Velocity2);
	if (Pressure < 0.0) {
		Pressure = Pressure_Old;
		cout <<"Negative pressure at point: ";
		for (unsigned short iDim = 0; iDim < nDim; iDim++) cout << Coord[iDim] <<" ";
		cout << endl;
		}
}

inline void CEulerVariable::SetSoundSpeed(double Gamma) { SoundSpeed = sqrt(Gamma*Pressure/Solution[0]); }

inline void CEulerVariable::SetTemperature(double Gas_Constant) { Temperature = Pressure / ( Gas_Constant * Solution[0]); }

inline void CEulerVariable::SetVelocity(double *val_velocity, bool val_incomp) { 
	if (val_incomp) {
		for (unsigned short iDim = 0; iDim < nDim; iDim++) 
			Solution[iDim+1] = val_velocity[iDim]*Solution[nDim+1]; 
	}
	else {
		for (unsigned short iDim = 0; iDim < nDim; iDim++) 
			Solution[iDim+1] = val_velocity[iDim]*Solution[0]; 
	}
};

inline void CEulerVariable::SetVelocity2(void) { Velocity2 = 0.0; for (unsigned short iDim = 0; iDim < nDim; iDim++) Velocity2 += Solution[iDim+1]*Solution[iDim+1]/(Solution[0]*Solution[0]); }

inline void CEulerVariable::SetVelocityInc2(void) { Velocity2 = 0.0; for (unsigned short iDim = 0; iDim < nDim; iDim++) Velocity2 += (Solution[iDim+1]/Solution[nDim+1])*(Solution[iDim+1]/Solution[nDim+1]); }

inline void CEulerVariable::SetBetaInc2(CConfig *config) { Solution[nDim+2] = config->GetArtComp_Factor(); }

inline void CEulerVariable::SetPressureInc(void) { Pressure = Solution[0]; }

inline void CEulerVariable::SetPressureInc(double val_pressure) { Solution[0] = val_pressure; Pressure = val_pressure; }

inline void CEulerVariable::SetVelocity_Old(double *val_velocity, bool val_incomp) { 
	if (val_incomp) {
		for (unsigned short iDim = 0; iDim < nDim; iDim++)	
			Solution_Old[iDim+1] = val_velocity[iDim]*Solution[nDim+1]; 
	}
	else {
		for (unsigned short iDim = 0; iDim < nDim; iDim++)	
			Solution_Old[iDim+1] = val_velocity[iDim]*Solution[0]; 
	}
};

inline void CEulerVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] += val_value; }

inline void CEulerVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] -= val_value; }

inline double CEulerVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return Gradient_Primitive[val_var][val_dim]; }

inline void CEulerVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] = val_value; }

inline double **CEulerVariable::GetGradient_Primitive(void) { return Gradient_Primitive; }

inline void CEulerVariable::SetTimeSpectral_Source(unsigned short val_var, double val_source) { TS_Source[val_var] = val_source; }

inline double CEulerVariable::GetTimeSpectral_Source(unsigned short val_var) { return TS_Source[val_var]; }

inline double CEulerVariable::GetPreconditioner_Beta() { return Precond_Beta; }

inline void CEulerVariable::SetPreconditioner_Beta(double val_Beta) { Precond_Beta = val_Beta; }

inline double CNSVariable::GetEddyViscosity(void) { return EddyViscosity; }

inline double CNSVariable::GetLaminarViscosity(void) { return LaminarViscosity; }

inline double CNSVariable::GetLaminarViscosityInc(void) { return LaminarViscosityInc; }

inline double CNSVariable::GetVorticity(unsigned short val_dim) { return Vorticity[val_dim]; }

inline double CNSVariable::GetStrainMag(void) { return StrainMag; }

inline void CNSVariable::SetLaminarViscosity(double val_laminar_viscosity) { LaminarViscosity = val_laminar_viscosity; }

inline void CNSVariable::SetLaminarViscosityInc(double val_laminar_viscosity_inc) { LaminarViscosityInc = val_laminar_viscosity_inc; }

inline void CNSVariable::SetEddyViscosity(double val_eddy_viscosity) { EddyViscosity = val_eddy_viscosity; }

inline void CNSVariable::SetWallTemperature(double Temperature_Wall ) { Temperature = Temperature_Wall; }

inline double *CAdjEulerVariable::GetForceProj_Vector(void) { return ForceProj_Vector; }

inline double *CAdjEulerVariable::GetObjFuncSource(void) { return ObjFuncSource; }

inline double *CAdjEulerVariable::GetIntBoundary_Jump(void) { return IntBoundary_Jump; }

inline double CAdjEulerVariable::GetTheta(void) { return Theta; }

inline void CAdjEulerVariable::SetForceProj_Vector(double *val_ForceProj_Vector) { for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

inline void CAdjEulerVariable::SetObjFuncSource(double *val_ObjFuncSource) { for (unsigned short iVar = 0; iVar < nVar; iVar++) ObjFuncSource[iVar] = val_ObjFuncSource[iVar]; }

inline void CAdjEulerVariable::SetIntBoundary_Jump(double *val_IntBoundary_Jump) { for (unsigned short iVar = 0; iVar < nVar; iVar++) IntBoundary_Jump[iVar] = val_IntBoundary_Jump[iVar]; }

inline void CAdjEulerVariable::SetPhi_Old(double *val_phi) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1]=val_phi[iDim]; };

inline double *CAdjNSVariable::GetForceProj_Vector(void) { return ForceProj_Vector; }

inline double CAdjNSVariable::GetTheta(void) { return Theta; }

inline void CAdjNSVariable::SetForceProj_Vector(double *val_ForceProj_Vector) {	for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

inline void CAdjNSVariable::SetPhi_Old(double *val_phi) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1]=val_phi[iDim]; };

inline void CAdjNSVariable::SetVelSolutionOldDVector(void) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1] = ForceProj_Vector[iDim]; };

inline void CAdjNSVariable::SetVelSolutionDVector(void) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution[iDim+1] = ForceProj_Vector[iDim]; };

inline double *CLinEulerVariable::GetForceProj_Vector(void) { return ForceProj_Vector; }

inline void CLinEulerVariable::SetForceProj_Vector(double *val_ForceProj_Vector) { for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

inline void CLinEulerVariable::SetDeltaVel_Old(double *val_deltavel) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1]=val_deltavel[iDim]; };

inline double CLinEulerVariable::GetDeltaPressure(void) { return DeltaPressure; }

inline double CPlasmaVariable::GetPressure(unsigned short val_iSpecies) { return Pressure[val_iSpecies]; }

inline double CPlasmaVariable::GetVelocity2(unsigned short val_iSpecies) { return Velocity2[val_iSpecies]; }

inline void CPlasmaVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] += val_value; }

inline void CPlasmaVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] -= val_value; }

inline double **CPlasmaVariable::GetGradient_Primitive(void) { return Gradient_Primitive; }

inline double CPlasmaVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return Gradient_Primitive[val_var][val_dim]; }

inline void CPlasmaVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] = val_value; }

inline double CPlasmaVariable::GetSoundSpeed(unsigned short val_species) { return SoundSpeed[val_species]; }

inline double CPlasmaVariable::GetEnthalpy(unsigned short val_species) { return Enthalpy[val_species]; }

inline double CPlasmaVariable::GetMax_Lambda_Inv(unsigned short iSpecies) { return Max_Lambda_Inv_MultiSpecies[iSpecies]; }

inline double CPlasmaVariable::GetMax_Lambda_Visc(unsigned short iSpecies) { return Max_Lambda_Visc_MultiSpecies[iSpecies]; }

inline void CPlasmaVariable::SetMax_Lambda_Inv(double val_max_lambda, unsigned short val_species) { Max_Lambda_Inv_MultiSpecies[val_species] = val_max_lambda; }

inline void CPlasmaVariable::SetMax_Lambda_Visc(double val_max_lambda, unsigned short val_species) { Max_Lambda_Visc_MultiSpecies[val_species] = val_max_lambda;}

inline void CPlasmaVariable::AddMax_Lambda_Inv(double val_max_lambda, unsigned short iSpecies) { Max_Lambda_Inv_MultiSpecies[iSpecies] += val_max_lambda; }

inline void CPlasmaVariable::AddMax_Lambda_Visc(double val_max_lambda, unsigned short iSpecies) { Max_Lambda_Visc_MultiSpecies[iSpecies] += val_max_lambda; }

inline void CPlasmaVariable::SetLambda(double val_lambda, unsigned short iSpecies) { Lambda[iSpecies] = val_lambda; }

inline void CPlasmaVariable::AddLambda(double val_lambda, unsigned short iSpecies) { Lambda[iSpecies] += val_lambda; }

inline double CPlasmaVariable::GetLambda(unsigned short iSpecies) { return Lambda[iSpecies]; }

inline double CPlasmaVariable::GetLaminarViscosity(unsigned short iSpecies) { return LaminarViscosity_MultiSpecies[iSpecies]; }

inline double CPlasmaVariable::GetEddyViscosity(unsigned short iSpecies) { return EddyViscosity_MultiSpecies[iSpecies]; }

inline void CPlasmaVariable::SetLaminarViscosity(double val_laminar_viscosity, unsigned short iSpecies ) { LaminarViscosity_MultiSpecies[iSpecies] = val_laminar_viscosity; }

inline double CPlasmaVariable::GetTemperature(unsigned short iSpecies) { return Temperature[iSpecies]; }

inline double CPlasmaVariable::GetTemperature_TR(unsigned short iSpecies) { return Temperature_tr[iSpecies]; }

inline double CPlasmaVariable::GetTemperature_vib(unsigned short iSpecies) { return Temperature_vib[iSpecies]; }

inline double CPlasmaVariable::GetDensity(unsigned short iSpecies) { return Solution[(nDim+2)*iSpecies + 0]; }

inline void CPlasmaVariable::SetVelocity_Old(double *val_velocity, unsigned short iSpecies) { for (unsigned short iDim = 0; iDim < nDim; iDim++)	Solution_Old[iSpecies*(nDim+2)+ iDim+1] = val_velocity[iDim]*Solution[iSpecies*(nDim+2)+ 0]; }

inline void CPlasmaVariable::SetWallTemperature(double* Temperature_Wall ) { for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) Temperature[iSpecies] = Temperature_Wall[iSpecies]; }

inline void CPlasmaVariable::SetDelta_Time(double val_delta_time, unsigned short iSpecies) { Species_Delta_Time[iSpecies] = val_delta_time;}

inline double CPlasmaVariable::GetDelta_Time(unsigned short iSpecies) { return Species_Delta_Time[iSpecies];}

inline void CPlasmaVariable::SetSensor(double val_sensor, unsigned short iSpecies) {Sensor_MultiSpecies[iSpecies] = val_sensor;}

inline double CPlasmaVariable::GetSensor(unsigned short iSpecies) {return Sensor_MultiSpecies[iSpecies]; }

inline void CPlasmaVariable::SetElectricField(double* ElectricField) {Elec_Field = ElectricField; }

inline double* CPlasmaVariable::GetElectricField() { return Elec_Field;}

inline double *CAdjPlasmaVariable::GetForceProj_Vector(void) { return ForceProj_Vector; }

inline double *CAdjPlasmaVariable::GetIntBoundary_Jump(void) { return IntBoundary_Jump; }

inline double CAdjPlasmaVariable::GetTheta(void) { return Theta; }

inline void CAdjPlasmaVariable::SetForceProj_Vector(double *val_ForceProj_Vector) { for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

inline void CAdjPlasmaVariable::SetIntBoundary_Jump(double *val_IntBoundary_Jump) { for (unsigned short iVar = 0; iVar < nVar; iVar++) IntBoundary_Jump[iVar] = val_IntBoundary_Jump[iVar]; }

inline void CAdjPlasmaVariable::SetPhi_Old(double *val_phi) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1]=val_phi[iDim]; };

inline void CLevelSetVariable::SetDiffLevelSet(double val_difflevelset) { DiffLevelSet = val_difflevelset; }

inline double CLevelSetVariable::GetDiffLevelSet(void) { return DiffLevelSet; }

inline double CFEAVariable::GetPressure(void) { return Pressure; }

inline void CFEAVariable::SetPressure(double val_pressure) { Pressure = val_pressure; }

inline void CWaveVariable::SetThickness_Noise(double val_thickness_noise) { Thickness_Noise = val_thickness_noise; }

inline void CWaveVariable::SetLoading_Noise(double val_loading_noise) { Loading_Noise = val_loading_noise; }

inline void CWaveVariable::SetQuadrupole_Noise(double val_quadrupole_noise) { Quadrupole_Noise = val_quadrupole_noise; }

inline double CWaveVariable::GetThickness_Noise() { return Thickness_Noise;}

inline double CWaveVariable::GetLoading_Noise() { return Loading_Noise;}

inline double CWaveVariable::GetQuadrupole_Noise() { return Quadrupole_Noise;}

inline double* CWaveVariable::GetSolution_Direct() { return Solution_Direct;}

inline void CWaveVariable::SetSolution_Direct(double *val_solution_direct) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Direct[iVar] += val_solution_direct[iVar];}

inline double* CPotentialVariable::GetChargeDensity() { return Charge_Density;}

inline void CPotentialVariable::SetChargeDensity(double positive_charge, double negative_charge) {Charge_Density[0] = positive_charge; Charge_Density[1] = negative_charge;}

inline void CPotentialVariable::SetPlasmaRhoUGradient(unsigned short iSpecies, double val_gradient, unsigned short iDim) { PlasmaRhoUGradient[iSpecies][iDim] = val_gradient;}

inline double** CPotentialVariable::GetPlasmaRhoUGradient() { return PlasmaRhoUGradient;}

inline void CPotentialVariable::SetPlasmaTimeStep(double dt) { PlasmaTimeStep = dt;}

inline double CPotentialVariable::GetPlasmaTimeStep() { return PlasmaTimeStep;}

inline double* CHeatVariable::GetSolution_Direct() { return Solution_Direct;}

inline void CHeatVariable::SetSolution_Direct(double *val_solution_direct) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Direct[iVar] += val_solution_direct[iVar];}

inline void CTurbSAVariable::SetTimeSpectral_Source(unsigned short val_var, double val_source) { TS_Source[val_var] = val_source; }

inline double CTurbSAVariable::GetTimeSpectral_Source(unsigned short val_var) { return TS_Source[val_var]; }
