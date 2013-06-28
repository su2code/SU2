/*!
 * \file variable_structure.inl
 * \brief In-Line subroutines of the <i>variable_structure.hpp</i> file.
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

#pragma once

inline void CVariable::SetVelocityInc2(void) { }

inline void CVariable::SetBetaInc2(CConfig *config) { }

inline void CVariable::SetPressureInc(void) { }

inline void CVariable::SetPressureInc(double val_pressure) { }

inline void CVariable::SetSoundSpeedInc(double *val_normal) { }

inline double CVariable::GetBetaInc2(void) { return 0; }

inline double CVariable::GetDensityInc(void) { return 0; }

inline double CVariable::GetSoundSpeedInc(void) { return 0; }

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

inline double *CVariable::GetResidual(void) { return Residual; }

inline double *CVariable::GetResConv(void) { return Res_Conv; }

inline double *CVariable::GetResVisc(void) { return Res_Visc; }

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

inline double **CVariable::GetGradient(void) { return Gradient; }

inline double *CVariable::GetLimiter(void) { return Limiter; }

inline double *CVariable::GetAuxVarGradient(void) { return Grad_AuxVar; }

inline double CVariable::GetAuxVarGradient(unsigned short val_dim) { return Grad_AuxVar[val_dim]; }

inline double *CVariable::GetTruncationError(void) { return TruncationError; }

inline void CVariable::SetDelta_Time(double val_delta_time) { Delta_Time = val_delta_time; }

inline double CVariable::GetDelta_Time(void) { return Delta_Time; }

inline void CVariable::SetMax_Lambda(double val_max_lambda) { Max_Lambda = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Inv(double val_max_lambda) { Max_Lambda_Inv = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Inv(double val_max_lambda, unsigned short val_Species) { }

inline void CVariable::SetMax_Lambda_Visc(double val_max_lambda) { Max_Lambda_Visc = val_max_lambda; }

inline void CVariable::SetLambda(double val_lambda) { Lambda = val_lambda; }

inline void CVariable::AddMax_Lambda(double val_max_lambda) { Max_Lambda += val_max_lambda; }

inline void CVariable::AddMax_Lambda_Inv(double val_max_lambda) { Max_Lambda_Inv += val_max_lambda; }

inline void CVariable::AddMax_Lambda_Visc(double val_max_lambda) { Max_Lambda_Visc += val_max_lambda; }

inline void CVariable::AddLambda(double val_lambda) { Lambda += val_lambda; }

inline double CVariable::GetMax_Lambda(void) { return Max_Lambda; }

inline double CVariable::GetMax_Lambda_Inv(void) { return Max_Lambda_Inv; }

inline double CVariable::GetMax_Lambda_Visc(void) { return Max_Lambda_Visc; }

inline double CVariable::GetLambda(void) { return Lambda; }

inline double CVariable::GetSensor(void) { return Sensor; }

inline double CVariable::GetMax_Lambda_Inv(unsigned short iFluids) { return 0; }

inline void CVariable::AddMax_Lambda_Inv(double val_max_lambda, unsigned short iFluids) { }

inline void CVariable::SetSensor(double val_sensor) { Sensor = val_sensor; }

inline void CVariable::AddResidual(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual[iVar] += val_residual[iVar]; }

inline void CVariable::AddResidual(double val_residual) { Residual[0] += val_residual ; }

inline void CVariable::AddRes_Conv(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar] += val_residual[iVar]; }

inline void CVariable::AddRes_Visc(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar] += val_residual[iVar]; }

inline void CVariable::SubtractResidual(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual[iVar] -= val_residual[iVar]; }

inline void CVariable::SubtractRes_Visc(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar] -= val_residual[iVar]; }

inline void CVariable::SubtractRes_Conv(double *val_residual) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar] -= val_residual[iVar]; }

inline double CVariable::GetDensity(void) {	return 0; }

inline double CVariable::GetEnergy(void) { return 0; }

inline double *CVariable::GetForceProj_Vector(void) { return NULL; }

inline double *CVariable::GetIntBoundary_Jump(void) { return NULL; }

inline double CVariable::GetEddyViscosity(void) { return 0; }

inline double CVariable::GetEnthalpy(void) { return 0; }

inline double CVariable::GetPressure(void) { return 0; }

inline double CVariable::GetDeltaPressure(void) { return 0; }

inline double CVariable::GetProjVel(double *val_vector) { return 0; }

inline double CVariable::GetProjVel(double *val_vector, unsigned short val_var) { return 0; }

inline double CVariable::GetProjVel(double *val_vector, unsigned short val_Species, unsigned short val_nDiatomics) { return 0; }

inline double CVariable::GetSoundSpeed(void) { return 0; }

inline double CVariable::GetSoundSpeed(unsigned short val_var) { return 0; }

inline double CVariable::GetTemperature(void) { return 0; }

inline double CVariable::GetThermalCoeff(void) { return 0; }

inline double CVariable::GetTheta(void) { return 0; }

inline double CVariable::GetVelocity(unsigned short val_dim) { return 0; }

inline double CVariable::GetVelocity2(void) { return 0; }

inline void CVariable::SetTemperature(unsigned short val_Fluid) { }

inline double CVariable::GetPressure(unsigned short val_Fluid) { return 0;}

inline double CVariable::GetVelocity2(unsigned short val_Fluid) { return 0;}

inline double CVariable::GetVelocity(unsigned short val_dim,unsigned short val_Fluid) { return 0;}

inline double CVariable::GetLaminarViscosity(void) { return 0; }

inline double CVariable::GetLaminarViscosityInc(void) { return 0; }

inline double CVariable::GetVorticity(unsigned short val_dim) { return 0; }

inline void CVariable::SetForceProj_Vector(double *val_ForceProj_Vector) { }

inline void CVariable::SetIntBoundary_Jump(double *val_IntBoundary_Jump) { }

inline void CVariable::SetEddyViscosity(unsigned short val_Kind_Turb_Model, double *Turb_Solution) { }

inline void CVariable::SetEnthalpy(void) { }

inline void CVariable::SetBetaInc2(double val_betainc2) { }

inline void CVariable::SetDensityInc(double val_densityinc) { }

inline void CVariable::SetPhi_Old(double *val_phi) { }

inline void CVariable::SetPressure(double Gamma, double *Coord) { }

inline void CVariable::SetPressure(double Gamma) { }

inline void CVariable::SetPressure(double GammaMonatomic, double GammaDiatomic) { }

inline void CVariable::SetPressure() { }

inline void CVariable::SetDeltaPressure(double *val_velocity, double Gamma) { }

inline void CVariable::SetSoundSpeed(double Gamma) { }

inline void CVariable::SetSoundSpeed(double GammaMonatomic, double GammaDiatomic) { }

inline void CVariable::SetSoundSpeed() { }

inline void CVariable::SetTemperature(double Gas_Constant) { }

inline void CVariable::SetThermalCoeff(double Gamma, double Gas_Constant) { }

inline void CVariable::SetTheta(double val_density, double *val_velocity, double val_enthalpy) { }

inline void CVariable::SetVelocity(double *val_velocity) { }

inline void CVariable::SetVelocity2(void) { }

inline void CVariable::SetVelocity_Old(double *val_velocity) { }

inline void CVariable::SetLaminarViscosity(CConfig *config) { }

inline void CVariable::SetLaminarViscosity(double val_laminar_viscosity) { }

inline void CVariable::SetLaminarViscosityInc(double val_laminar_viscosity_inc) { }

inline void CVariable::SetEddyViscosity(double val_eddy_viscosity) { }

inline void CVariable::SetVorticity(void) { }

inline void CVariable::SetGradient_PrimitiveZero(void) { }

inline void CVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { }

inline void CVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { }

inline double CVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return 0; }

inline void CVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { }

inline double **CVariable::GetGradient_Primitive(void) { return NULL; }

inline double CVariable::GetSource(unsigned short val_var) { return 0; }

inline double CVariable::GetSource_Jacobian(unsigned short val_ivar,unsigned short val_jvar) { return 0;}

inline void CVariable::SetSource(double * res) { }

inline void CVariable::SetSourceJacobian(double ** src_jacobian) { }

inline void CVariable::SetF1blending(double val_viscosity, double val_dist, double val_density) { }

inline double CVariable::GetF1blending(void) { return 0; }

inline double CEulerVariable::GetDensity(void) { return Solution[0]; }

inline double CEulerVariable::GetDensityInc(void) { return DensityInc; }

inline double CEulerVariable::GetBetaInc2(void) { return BetaInc2; }

inline double CEulerVariable::GetEnergy(void) { return Solution[nVar-1]/Solution[0]; };

inline double CEulerVariable::GetEnthalpy(void) { return Enthalpy; }

inline double CEulerVariable::GetPressure(void) { return Pressure; }

inline double CEulerVariable::GetSoundSpeed(void) { return SoundSpeed; }

inline double CEulerVariable::GetTemperature(void) { return Temperature; }

inline double CEulerVariable::GetThermalCoeff(void) { return Kappa; }

inline double CEulerVariable::GetVelocity(unsigned short val_dim) { return Solution[val_dim+1]/Solution[0]; }

inline double CEulerVariable::GetVelocity2(void) { return Velocity2; }

inline void CEulerVariable::SetEnthalpy(void) { Enthalpy = (Solution[nVar-1] + Pressure) / Solution[0]; }

inline void CEulerVariable::SetDensityInc(double val_density) { DensityInc = val_density; }

inline void CEulerVariable::SetBetaInc2(double val_betainc2) { BetaInc2 = val_betainc2; }

inline void CEulerVariable::SetPressure(double Gamma, double *Coord) {
	double OldPressure = Pressure;
	Pressure = (Gamma-1.0)*Solution[0]*(Solution[nVar-1]/Solution[0]-0.5*Velocity2);
	if (Pressure < 0.0) {
		Pressure = OldPressure;
		cout <<"Negative pressure at point: ";
		for (unsigned short iDim = 0; iDim < nDim; iDim++) cout << Coord[iDim] <<" ";
		cout << endl;
		}
}

inline void CEulerVariable::SetSoundSpeed(double Gamma) { SoundSpeed = sqrt(Gamma*Pressure/Solution[0]); }

inline void CEulerVariable::SetTemperature(double Gas_Constant) { Temperature = Pressure / ( Gas_Constant * Solution[0]); }

inline void CEulerVariable::SetVelocity(double *val_velocity) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution[iDim+1] = val_velocity[iDim]*Solution[0]; };

inline void CEulerVariable::SetVelocity2(void) { Velocity2 = 0; for (unsigned short iDim = 0; iDim < nDim; iDim++) Velocity2 += Solution[iDim+1]*Solution[iDim+1]/(Solution[0]*Solution[0]); }

inline void CEulerVariable::SetVelocityInc2(void) { Velocity2 = 0; for (unsigned short iDim = 0; iDim < nDim; iDim++) Velocity2 += Solution[iDim+1]*Solution[iDim+1]; }

inline void CEulerVariable::SetBetaInc2(CConfig *config) { 
//		double rb = config->GetArtComp_Factor(); 
//		double BetaInc2_min = config->GetArtComp_Min(); 
//		BetaInc2 = rb*min(Velocity2, BetaInc2_min);
//		BetaInc2 = fabs(rb*Pressure);
//		BetaInc2 = max(fabs(rb*Pressure), rb);
//		BetaInc2 = max(1.0,rb*config->GetPressure_FreeStreamND());
//		BetaInc2 = config->GetArtComp_Factor()/config->GetDensity_FreeStreamND();
//		BetaInc2 = config->GetArtComp_Factor()/DensityInc;
		BetaInc2 = config->GetArtComp_Factor();
 }

inline void CEulerVariable::SetPressureInc(void) { Pressure = Solution[0]; }

inline void CEulerVariable::SetPressureInc(double val_pressure) { Solution[0] = val_pressure; Pressure = val_pressure; }

inline double CEulerVariable::GetSoundSpeedInc(void) { return SoundSpeedInc; }

inline void CEulerVariable::SetVelocity_Old(double *val_velocity) { for (unsigned short iDim = 0; iDim < nDim; iDim++)	Solution_Old[iDim+1] = val_velocity[iDim]*Solution[0]; };

inline void CEulerVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] += val_value; }

inline void CEulerVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] -= val_value; }

inline double CEulerVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return Gradient_Primitive[val_var][val_dim]; }

inline void CEulerVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] = val_value; }

inline double **CEulerVariable::GetGradient_Primitive(void) { return Gradient_Primitive; }

inline double CNSVariable::GetEddyViscosity(void) { return EddyViscosity; }

inline double CNSVariable::GetLaminarViscosity(void) { return LaminarViscosity; }

inline double CNSVariable::GetLaminarViscosityInc(void) { return LaminarViscosityInc; }

inline double CNSVariable::GetVorticity(unsigned short val_dim) { return Vorticity[val_dim]; }

inline void CNSVariable::SetLaminarViscosity(double val_laminar_viscosity) { LaminarViscosity = val_laminar_viscosity; }

inline void CNSVariable::SetLaminarViscosityInc(double val_laminar_viscosity_inc) { LaminarViscosityInc = val_laminar_viscosity_inc; }

inline void CNSVariable::SetEddyViscosity(double val_eddy_viscosity) { EddyViscosity = val_eddy_viscosity; }

inline void CNSVariable::SetThermalCoeff(double Gamma, double Gas_Constant ) { Kappa = Gamma * LaminarViscosity / ( (Gamma-1.0)*Gas_Constant*Prandtl_Lam ); }

inline double *CAdjEulerVariable::GetForceProj_Vector(void) { return ForceProj_Vector; }

inline double *CAdjEulerVariable::GetIntBoundary_Jump(void) { return IntBoundary_Jump; }

inline double CAdjEulerVariable::GetTheta(void) { return Theta; }

inline void CAdjEulerVariable::SetForceProj_Vector(double *val_ForceProj_Vector) { for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

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

inline void CPlasmaMonatomicVariable::SetTemperature(unsigned short val_Fluid) { }

inline double CPlasmaMonatomicVariable::GetPressure(unsigned short val_Fluid) { return Pressure[val_Fluid]; }

inline double CPlasmaMonatomicVariable::GetVelocity2(unsigned short val_Fluid) { return Velocity2[val_Fluid]; }

inline void CPlasmaMonatomicVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] += val_value; }

inline void CPlasmaMonatomicVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] -= val_value; }

inline double **CPlasmaMonatomicVariable::GetGradient_Primitive(void) { return Gradient_Primitive; }

inline double CPlasmaMonatomicVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return Gradient_Primitive[val_var][val_dim]; }

inline void CPlasmaMonatomicVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] = val_value; }

inline double CPlasmaMonatomicVariable::GetSoundSpeed(unsigned short val_Fluid) { return SoundSpeed[val_Fluid]; }

inline double CPlasmaMonatomicVariable::GetSource(unsigned short val_var) { return Source[val_var]; }

inline double CPlasmaMonatomicVariable::GetSource_Jacobian(unsigned short val_ivar,unsigned short val_jvar) { return Source_Jacobian[val_ivar][val_jvar]; }

inline double CPlasmaMonatomicVariable::GetMax_Lambda_Inv(unsigned short iFluids) { return Max_Lambda_Inv_MultiSpecies[iFluids]; }

inline void CPlasmaMonatomicVariable::AddMax_Lambda_Inv(double val_max_lambda, unsigned short iFluids) { Max_Lambda_Inv_MultiSpecies[iFluids] += val_max_lambda; }

inline void CPlasmaDiatomicVariable::SetTemperature(unsigned short val_Fluid) { }

inline double CPlasmaDiatomicVariable::GetPressure(unsigned short val_Fluid) { return Pressure[val_Fluid]; }

inline double CPlasmaDiatomicVariable::GetVelocity2(unsigned short val_Fluid) { return Velocity2[val_Fluid]; }

inline void CPlasmaDiatomicVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] += val_value; }

inline void CPlasmaDiatomicVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] -= val_value; }

inline double **CPlasmaDiatomicVariable::GetGradient_Primitive(void) { return Gradient_Primitive; }

inline double CPlasmaDiatomicVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return Gradient_Primitive[val_var][val_dim]; }

inline void CPlasmaDiatomicVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] = val_value; }

inline double CPlasmaDiatomicVariable::GetSoundSpeed(unsigned short val_Species) { return SoundSpeed[val_Species]; }

inline double CPlasmaDiatomicVariable::GetSource(unsigned short val_var) { return Source[val_var]; }

inline double CPlasmaDiatomicVariable::GetSource_Jacobian(unsigned short val_ivar,unsigned short val_jvar) { return Source_Jacobian[val_ivar][val_jvar]; }

inline double CPlasmaDiatomicVariable::GetMax_Lambda_Inv(unsigned short iSpecies) { return Max_Lambda_Inv_MultiSpecies[iSpecies]; }

inline void CPlasmaDiatomicVariable::SetMax_Lambda_Inv(double val_max_lambda, unsigned short val_Species) { Max_Lambda_Inv_MultiSpecies[val_Species] = val_max_lambda; }

inline void CPlasmaDiatomicVariable::AddMax_Lambda_Inv(double val_max_lambda, unsigned short iSpecies) { Max_Lambda_Inv_MultiSpecies[iSpecies] += val_max_lambda; }
