/*!
 * \file numerics_structure.inl
 * \brief In-Line subroutines of the <i>numerics_structure.hpp</i> file.
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

#pragma once

inline double CNumerics::Determinant_3x3(double A00, double A01, double A02, double A10, double A11, double A12, double A20, double A21, double A22) {
	return A00*(A11*A22-A12*A21) - A01*(A10*A22-A12*A20) + A02*(A10*A21-A11*A20);
}

inline void CNumerics::ComputeResidual(double *val_residual, CConfig *config) { }

inline void CNumerics::ComputeResidual_MacCormack(double *val_residual, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual_i, double *val_residual_j) { }

inline void CNumerics::ComputeResidual(double *val_residual_i, double *val_residual_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, 
                                   CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                   double **val_JacobianMeanFlow_i, double **val_JacobianMeanFlow_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, 
								   double **val_Jacobian_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual_i, double *val_residual_j, 
								   double **val_Jacobian_ii, double **val_Jacobian_ij, 
								   double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) { }
							
inline void CNumerics::ComputeResidual(double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, 
								   double *val_resvisc_j, double **val_Jacobian_ii, double **val_Jacobian_ij, 
								   double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) { }
							
inline void CNumerics::ComputeResidual(double **val_stiffmatrix_elem, CConfig *config) { }

inline void CNumerics::GetEq_Rxn_Coefficients(double **EqnRxnConstants, CConfig *config) { };
														
inline void CNumerics::ComputeResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::ComputeResidual_TransLM(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config, double &gamma_sep) {}

inline void CNumerics::ComputeResidual_Axisymmetric(double *val_residual, CConfig *config) { }

inline void CNumerics::ComputeResidual_Axisymmetric_ad(double *val_residual, double *val_residuald, CConfig *config) { }

inline void CNumerics::SetJacobian_Axisymmetric(double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::ComputeResidual_Chemistry(double *val_residual, CConfig *config) { }

inline void CNumerics::ComputeResidual_Chemistry_ad(double *val_residual, double *val_residuald, CConfig *config) { }

inline void CNumerics::ComputeResidual_Chemistry(double *val_residual, double **val_Jacobian, CConfig *config) { }

inline void CNumerics::SetJacobian_Chemistry(double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::ComputeResidual_ElecForce(double *val_residual, CConfig *config) { }

inline void CNumerics::ComputeResidual_ElecForce(double *val_residual, double **val_Jacobian, CConfig *config) { }

inline void CNumerics::SetJacobian_ElecForce(double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::ComputeResidual_MomentumExch(double *val_residual, CConfig *config) { }

inline void CNumerics::ComputeResidual_MomentumExch_ad(double *val_residual, double *val_residuald, CConfig *config) { }

inline void CNumerics::ComputeResidual_MomentumExch(double *val_residual, double **val_Jacobian, CConfig *config) { }

inline void CNumerics::SetJacobian_MomentumExch(double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::ComputeResidual_EnergyExch(double *val_residual, double **val_Jacobian, CConfig *config) { }

inline void CNumerics::ComputeResidual_EnergyExch(double *val_residual, double *val_residual_ElecForce, double **val_Jacobian, CConfig *config) { }

inline void CNumerics::ComputeResidual_EnergyExch(double *val_residual, CConfig *config) { }

inline void CNumerics::ComputeResidual_EnergyExch_ad(double *val_residual, double *val_residuald, CConfig *config) { }

inline void CNumerics::SetJacobian_EnergyExch(double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::SetSensor( double val_sensor_i, double val_sensor_j, unsigned short iSpecies) {}

inline void CNumerics::SetPressure(double val_pressure_i, double val_pressure_j, unsigned short iSpecies) { }

inline void CNumerics::SetSoundSpeed(double val_soundspeed_i, double val_soundspeed_j, unsigned short iSpecies) { }

inline void CNumerics::SetEnthalpy(double val_enthalpy_i, double val_enthalpy_j, unsigned short iSpecies) { }

inline void CNumerics::SetLambda(double val_lambda_i, double val_lambda_j, unsigned short iSpecies) { }

inline double CNumerics::GetPrecond_Beta() { return 0; }
	
inline void CNumerics::SetUndivided_Laplacian(double *val_und_lapl_i, double *val_und_lapl_j) {
	Und_Lapl_i = val_und_lapl_i; 
	Und_Lapl_j = val_und_lapl_j; 
}

inline void CNumerics::SetSensor( double val_sensor_i, double val_sensor_j) {
	Sensor_i = val_sensor_i;
	Sensor_j = val_sensor_j;
}

inline void CCentJST_Plasma::SetSensor( double val_sensor_i, double val_sensor_j, unsigned short iSpecies) {
	Sensor_i[iSpecies] = val_sensor_i;
	Sensor_j[iSpecies] = val_sensor_j;
}

inline void CCentJST_PlasmaDiatomic::SetSensor( double val_sensor_i, double val_sensor_j, unsigned short iSpecies) {
	Sensor_i[iSpecies] = val_sensor_i;
	Sensor_j[iSpecies] = val_sensor_j;
}

inline void CCentLax_PlasmaDiatomic::SetSensor( double val_sensor_i, double val_sensor_j, unsigned short iSpecies) {
	Sensor_i[iSpecies] = val_sensor_i;
	Sensor_j[iSpecies] = val_sensor_j;
}

inline void CCentLax_AdjPlasmaDiatomic::SetSensor( double val_sensor_i, double val_sensor_j, unsigned short iSpecies) {
	Sensor_i[iSpecies] = val_sensor_i;
	Sensor_j[iSpecies] = val_sensor_j;
}

inline void CNumerics::SetConservative(double *val_u_i, double *val_u_j) {
	U_i = val_u_i;
	U_j = val_u_j;
}

inline void CNumerics::SetConservative_ZeroOrder(double *val_u_i, double *val_u_j) {
	UZeroOrder_i = val_u_i;
	UZeroOrder_j = val_u_j;
}

inline void CNumerics::SetPrimitive(double *val_v_i, double *val_v_j) {
	V_i = val_v_i;
	V_j = val_v_j;
}

inline void CNumerics::SetPrimitive(double **val_v_i, double**val_v_j) {
  Varray_i = val_v_i;
  Varray_j = val_v_j;
}

inline void CNumerics::SetConservative(double *val_u_0, double *val_u_1, double *val_u_2) {
	U_0 = val_u_0;
	U_1 = val_u_1;
	U_2 = val_u_2;
}

inline void CNumerics::SetConservative(double *val_u_0, double *val_u_1, double *val_u_2, double *val_u_3) {
	U_0 = val_u_0;
	U_1 = val_u_1;
	U_2 = val_u_2;
	U_3 = val_u_3;
}

inline void CNumerics::SetVelocity2_Inf(double velocity2) {
	vel2_inf = velocity2;
}

inline void CNumerics::SetChargeDensity(double *val_u_0, double *val_u_1, double *val_u_2, double *val_u_3) {
	U_0 = val_u_0;
	U_1 = val_u_1;
	U_2 = val_u_2;
	U_3 = val_u_3;
}

inline void CNumerics::SetElecField(double *val_Efield) {}

inline void CNumerics::SetTimeStep(double val_timestep) {TimeStep = val_timestep;}

inline double* CNumerics::GetMagneticField() {return 0;}

inline double CNumerics::GetMagneticForce(unsigned short val_Species, unsigned short val_dim) {return 0;}

inline void CNumerics::SetElec_Cond() {}

inline void CNumerics::SetLaminarViscosity(double val_lam_viscosity_i, double val_lam_viscosity_j) {
	Laminar_Viscosity_i = val_lam_viscosity_i;
	Laminar_Viscosity_j = val_lam_viscosity_j;
}

inline void CNumerics::SetLaminarViscosity(double val_lam_viscosity_i, double val_lam_viscosity_j, unsigned short iSpecies) {
	Laminar_Viscosity_MultipleSpecies_i[iSpecies] = val_lam_viscosity_i;
	Laminar_Viscosity_MultipleSpecies_j[iSpecies] = val_lam_viscosity_j;
}

inline void CNumerics::SetThermalConductivity(double val_therm_conductivity_i, double val_therm_conductivity_j, unsigned short iSpecies) {
	Thermal_Conductivity_MultipleSpecies_i[iSpecies] = val_therm_conductivity_i;
	Thermal_Conductivity_MultipleSpecies_j[iSpecies] = val_therm_conductivity_j;
}

inline void CNumerics::SetThermalConductivity_vib(double val_therm_conductivity_vib_i, double val_therm_conductivity_vib_j, unsigned short iSpecies) {
	Thermal_Conductivity_vib_MultipleSpecies_i[iSpecies] = val_therm_conductivity_vib_i;
	Thermal_Conductivity_vib_MultipleSpecies_j[iSpecies] = val_therm_conductivity_vib_j;
}

inline void CNumerics::SetEddyViscosity(double val_eddy_viscosity_i, double val_eddy_viscosity_j, unsigned short iSpecies) {
	Eddy_Viscosity_MultipleSpecies_i[iSpecies] = val_eddy_viscosity_i;
	Eddy_Viscosity_MultipleSpecies_j[iSpecies] = val_eddy_viscosity_j;
}

inline void CNumerics::SetEddyViscosity(double val_eddy_viscosity_i, double val_eddy_viscosity_j) {
	Eddy_Viscosity_i=val_eddy_viscosity_i;
	Eddy_Viscosity_j=val_eddy_viscosity_j;
}

inline void CNumerics::SetIntermittency(double intermittency_in) { }

inline void CNumerics::SetTurbKineticEnergy(double val_turb_ke_i, double val_turb_ke_j) {
	turb_ke_i = val_turb_ke_i;
	turb_ke_j = val_turb_ke_j;
}

inline void CNumerics::SetDistance(double val_dist_i, double val_dist_j) {
	dist_i = val_dist_i;
	dist_j = val_dist_j;
}

inline void CNumerics::SetAdjointVar(double *val_psi_i, double *val_psi_j) {
	Psi_i = val_psi_i;
	Psi_j = val_psi_j;
}

inline void CNumerics::SetLinearizedVar(double *val_deltau_i, double *val_deltau_j) {
	DeltaU_i = val_deltau_i;
	DeltaU_j = val_deltau_j;
}

inline void CNumerics::SetAdjointVarGradient(double **val_psivar_grad_i, double **val_psivar_grad_j) {
	PsiVar_Grad_i = val_psivar_grad_i;
	PsiVar_Grad_j = val_psivar_grad_j;
}

inline void CNumerics::SetTurbVar(double *val_turbvar_i, double *val_turbvar_j) {
	TurbVar_i = val_turbvar_i;
	TurbVar_j = val_turbvar_j;
}

inline void CNumerics::SetTransVar(double *val_transvar_i, double *val_transvar_j) {
	TransVar_i = val_transvar_i;
	TransVar_j = val_transvar_j;
}

inline void CNumerics::SetTurbVarGradient(double **val_turbvar_grad_i, double **val_turbvar_grad_j) {
	TurbVar_Grad_i = val_turbvar_grad_i;
	TurbVar_Grad_j = val_turbvar_grad_j;
}

inline void CNumerics::SetTransVarGradient(double **val_transvar_grad_i, double **val_transvar_grad_j) {
	TransVar_Grad_i = val_transvar_grad_i;
	TransVar_Grad_j = val_transvar_grad_j;
}

inline void CNumerics::SetLevelSetVar(double *val_levelsetvar_i, double *val_levelsetvar_j) {
	LevelSetVar_i = val_levelsetvar_i;
	LevelSetVar_j = val_levelsetvar_j;
}

inline void CNumerics::SetLevelSetVarGradient(double **val_levelsetvar_grad_i, double **val_levelsetvar_grad_j) {
	LevelSetVar_Grad_i = val_levelsetvar_grad_i;
	LevelSetVar_Grad_j = val_levelsetvar_grad_j;
}

inline void CNumerics::SetPrimVarGradient(double **val_primvar_grad_i, double **val_primvar_grad_j) {
	PrimVar_Grad_i = val_primvar_grad_i;
	PrimVar_Grad_j = val_primvar_grad_j;
}

inline void CNumerics::SetPrimVarGradient(double ***val_primvar_grad_i, double ***val_primvar_grad_j) {
  PrimVar_Grad_i_array = val_primvar_grad_i;
  PrimVar_Grad_j_array = val_primvar_grad_j;
}

inline void CNumerics::SetConsVarGradient(double **val_consvar_grad_i, double **val_consvar_grad_j) {
	ConsVar_Grad_i = val_consvar_grad_i;
	ConsVar_Grad_j = val_consvar_grad_j;
}

inline void CNumerics::SetConsVarGradient(double **val_consvar_grad_0, double **val_consvar_grad_1, double **val_consvar_grad_2) {
	ConsVar_Grad_0 = val_consvar_grad_0;
	ConsVar_Grad_1 = val_consvar_grad_1;
	ConsVar_Grad_2 = val_consvar_grad_2;
}

inline void CNumerics::SetConsVarGradient(double **val_consvar_grad_0, double **val_consvar_grad_1, double **val_consvar_grad_2, double **val_consvar_grad_3) {
	ConsVar_Grad_0 = val_consvar_grad_0;
	ConsVar_Grad_1 = val_consvar_grad_1;
	ConsVar_Grad_2 = val_consvar_grad_2;
	ConsVar_Grad_3 = val_consvar_grad_3;
}

inline void CNumerics::SetConsVarGradient(double **val_consvar_grad) {
	ConsVar_Grad = val_consvar_grad;
}

inline void CNumerics::SetCoord(double *val_coord_i, double *val_coord_j) {
	Coord_i = val_coord_i;
	Coord_j = val_coord_j;
}

inline void CNumerics::SetCoord(double *val_coord_0, double *val_coord_1, 
									 double *val_coord_2) {
	Coord_0 = val_coord_0;
	Coord_1 = val_coord_1;
	Coord_2 = val_coord_2;
}

inline void CNumerics::SetCoord(double *val_coord_0, double *val_coord_1, 
									 double *val_coord_2, double *val_coord_3) {
	Coord_0 = val_coord_0;
	Coord_1 = val_coord_1;
	Coord_2 = val_coord_2;
	Coord_3 = val_coord_3;	
}

inline void CNumerics::SetGridVel(double *val_gridvel_i, double *val_gridvel_j) {
	GridVel_i = val_gridvel_i;
	GridVel_j = val_gridvel_j;
}

inline void CNumerics::SetRotVel(double *val_rotvel_i, double *val_rotvel_j) {
	RotVel_i = val_rotvel_i;
	RotVel_j = val_rotvel_j;
}

inline void CNumerics::SetRotFlux(double val_rot_flux) { Rot_Flux = val_rot_flux; }

inline void CNumerics::SetPressure(double val_pressure_i, double val_pressure_j) {
	Pressure_i = val_pressure_i;
	Pressure_j = val_pressure_j;
}

inline void CCentJST_Plasma::SetPressure(double val_pressure_i, double val_pressure_j, unsigned short iSpecies) {
	Pressure_i[iSpecies] = val_pressure_i;
	Pressure_j[iSpecies] = val_pressure_j;
}

inline void CCentJST_PlasmaDiatomic::SetPressure(double val_pressure_i, double val_pressure_j, unsigned short iSpecies) {
	Pressure_i[iSpecies] = val_pressure_i;
	Pressure_j[iSpecies] = val_pressure_j;
}

inline void CCentLax_PlasmaDiatomic::SetPressure(double val_pressure_i, double val_pressure_j, unsigned short iSpecies) {
	Pressure_i[iSpecies] = val_pressure_i;
	Pressure_j[iSpecies] = val_pressure_j;
}

inline void CCentLax_AdjPlasmaDiatomic::SetPressure(double val_pressure_i, double val_pressure_j, unsigned short iSpecies) {
	Pressure_i[iSpecies] = val_pressure_i;
	Pressure_j[iSpecies] = val_pressure_j;
}

inline void CNumerics::SetDensityInc(double val_densityinc_i, double val_densityinc_j) {
	DensityInc_i = val_densityinc_i;
	DensityInc_j = val_densityinc_j;
}

inline void CNumerics::SetBetaInc2(double val_betainc2_i, double val_betainc2_j) {
	BetaInc2_i = val_betainc2_i;
	BetaInc2_j = val_betainc2_j;
}

inline void CNumerics::SetSoundSpeed(double val_soundspeed_i, double val_soundspeed_j) {
	SoundSpeed_i = val_soundspeed_i;
	SoundSpeed_j = val_soundspeed_j;
}

inline void CCentJST_Plasma::SetSoundSpeed(double val_soundspeed_i, double val_soundspeed_j, unsigned short iSpecies) {
	SoundSpeed_i[iSpecies] = val_soundspeed_i;
	SoundSpeed_j[iSpecies] = val_soundspeed_j;
}

inline void CCentJST_PlasmaDiatomic::SetSoundSpeed(double val_soundspeed_i, double val_soundspeed_j, unsigned short iSpecies) {
	SoundSpeed_i[iSpecies] = val_soundspeed_i;
	SoundSpeed_j[iSpecies] = val_soundspeed_j;
}

inline void CCentLax_PlasmaDiatomic::SetSoundSpeed(double val_soundspeed_i, double val_soundspeed_j, unsigned short iSpecies) {
	SoundSpeed_i[iSpecies] = val_soundspeed_i;
	SoundSpeed_j[iSpecies] = val_soundspeed_j;
}

inline void CCentLax_AdjPlasmaDiatomic::SetSoundSpeed(double val_soundspeed_i, double val_soundspeed_j, unsigned short iSpecies) {
	SoundSpeed_i[iSpecies] = val_soundspeed_i;
	SoundSpeed_j[iSpecies] = val_soundspeed_j;
}

inline void CNumerics::SetEnthalpy(double val_enthalpy_i, double val_enthalpy_j) {
	Enthalpy_i = val_enthalpy_i;
	Enthalpy_j = val_enthalpy_j;
}

inline void CCentJST_Plasma::SetEnthalpy(double val_enthalpy_i, double val_enthalpy_j, unsigned short iSpecies) {
	Enthalpy_i[iSpecies] = val_enthalpy_i;
	Enthalpy_j[iSpecies] = val_enthalpy_j;
}

inline void CCentJST_PlasmaDiatomic::SetEnthalpy(double val_enthalpy_i, double val_enthalpy_j, unsigned short iSpecies) {
	Enthalpy_i[iSpecies] = val_enthalpy_i;
	Enthalpy_j[iSpecies] = val_enthalpy_j;
}

inline void CCentLax_PlasmaDiatomic::SetEnthalpy(double val_enthalpy_i, double val_enthalpy_j, unsigned short iSpecies) {
	Enthalpy_i[iSpecies] = val_enthalpy_i;
	Enthalpy_j[iSpecies] = val_enthalpy_j;
}

inline void CCentLax_AdjPlasmaDiatomic::SetEnthalpy(double val_enthalpy_i, double val_enthalpy_j, unsigned short iSpecies) {
	Enthalpy_i[iSpecies] = val_enthalpy_i;
	Enthalpy_j[iSpecies] = val_enthalpy_j;
}

inline void CNumerics::SetLambda(double val_lambda_i, double val_lambda_j) {
	Lambda_i = val_lambda_i;
	Lambda_j = val_lambda_j;
}

inline void CCentJST_Plasma::SetLambda(double val_lambda_i, double val_lambda_j, unsigned short iSpecies) {
	Lambda_i[iSpecies] = val_lambda_i;
	Lambda_j[iSpecies] = val_lambda_j;
}

inline void CCentJST_PlasmaDiatomic::SetLambda(double val_lambda_i, double val_lambda_j, unsigned short iSpecies) {
	Lambda_i[iSpecies] = val_lambda_i;
	Lambda_j[iSpecies] = val_lambda_j;
}

inline void CCentLax_PlasmaDiatomic::SetLambda(double val_lambda_i, double val_lambda_j, unsigned short iSpecies) {
	Lambda_i[iSpecies] = val_lambda_i;
	Lambda_j[iSpecies] = val_lambda_j;
}

inline void CCentLax_AdjPlasmaDiatomic::SetLambda(double val_lambda_i, double val_lambda_j, unsigned short iSpecies) {
	Lambda_i[iSpecies] = val_lambda_i;
	Lambda_j[iSpecies] = val_lambda_j;
}

inline void CNumerics::SetNeighbor(unsigned short val_neighbor_i, unsigned short val_neighbor_j) {
	Neighbor_i = val_neighbor_i;
	Neighbor_j = val_neighbor_j;
}

inline void CNumerics::SetTurbAdjointVar(double *val_turbpsivar_i, double *val_turbpsivar_j) {
	TurbPsi_i = val_turbpsivar_i;
	TurbPsi_j = val_turbpsivar_j;
}

inline void CNumerics::SetTurbAdjointGradient(double **val_turbpsivar_grad_i, double **val_turbpsivar_grad_j) {
	TurbPsi_Grad_i = val_turbpsivar_grad_i;
	TurbPsi_Grad_j = val_turbpsivar_grad_j;
}

inline void CNumerics::SetTemperature(double val_temp_i, double val_temp_j) {
	Temp_i = val_temp_i;
	Temp_j = val_temp_j;
}

inline void CNumerics::SetTemperature_tr(double* val_temp_i, double* val_temp_j) {
	Temp_tr_i = val_temp_i;
	Temp_tr_j = val_temp_j;
}

inline void CNumerics::SetTemperature_vib(double* val_temp_i, double* val_temp_j) {
	Temp_vib_i = val_temp_i;
	Temp_vib_j = val_temp_j;
}

inline void CNumerics::SetPressure(double* val_pressure_i, double* val_pressure_j) {
	SpeciesPressure_i = val_pressure_i;
	SpeciesPressure_j = val_pressure_j;
}

inline void CNumerics::SetAuxVarGrad(double *val_auxvargrad_i, double *val_auxvargrad_j) {
	AuxVar_Grad_i = val_auxvargrad_i;
	AuxVar_Grad_j = val_auxvargrad_j;
}

inline void CNumerics::SetNormal(double *val_normal) { Normal = val_normal; }

inline void CNumerics::SetVolume(double val_volume) { Volume = val_volume; }

inline void CSourcePieceWise_TurbSST::SetF1blending(double val_F1_i, double val_F1_j){ 
	F1_i = val_F1_i; 
	F1_j = val_F1_j;
}

inline void CSourcePieceWise_TurbSST::SetF2blending(double val_F2_i, double val_F2_j){ 
	F2_i = val_F2_i; 
	F2_j = val_F2_j;
}

inline void CSourcePieceWise_TurbSST::SetStrainMag(double val_StrainMag_i, double val_StrainMag_j){
	StrainMag = val_StrainMag_i;
}

inline void CSourcePieceWise_TurbSST::SetCrossDiff(double val_CDkw_i, double val_CDkw_j){
	CDkw = val_CDkw_i;
}			

inline void CSourcePieceWise_TurbSA::SetIntermittency(double intermittency_in) { intermittency = intermittency_in; }

inline void CSourcePieceWise_Plasma::SetElecField(double *val_Efield) { ElectricField = val_Efield; }

inline double* CSourcePieceWise_Plasma::GetMagneticField() { return JcrossB; }

inline double CSourcePieceWise_Plasma::GetMagneticForce(unsigned short val_Species, unsigned short val_dim) { return Mag_Force[val_Species][val_dim]; }

inline double CUpwRoe_Turkel_Flow::GetPrecond_Beta() { return Beta; }

inline double* CSource_Magnet::GetMagneticField() { return Current_Density; }

inline double CNumerics::GetElec_CondIntegral() {return 0;}

inline double CSource_JouleHeating::GetElec_CondIntegral() {return Elec_Conduct*Coord_i[1]*(Coord_i[1]-Coord_j[1]);}

inline void CNumerics::SetElec_CondIntegralsqr(double val_var) {}

inline void CSource_JouleHeating::SetElec_CondIntegralsqr(double val_var) {Integralsqr = val_var; }

inline double CNumerics::GetKappaPsiVolume() {return 0;}

inline double CSourceViscous_AdjFlow::GetKappaPsiVolume() { return kappapsi_Volume; }

inline void CNumerics::SetKappaPsiVolume(double val_kappapsi_Volume) {}

inline void CSourceViscous_AdjFlow::SetKappaPsiVolume(double val_kappapsi_Volume) { kappapsi_Volume = val_kappapsi_Volume;}

inline void CNumerics::ComputeResidual(double **val_Jacobian_i, double *val_Jacobian_mui, double ***val_Jacobian_gradi, CConfig *config) { }

inline void CNumerics::ComputeResidual(double **val_Jacobian_i, double *val_Jacobian_mui, double ***val_Jacobian_gradi, 
									double **val_Jacobian_j, double *val_Jacobian_muj, double ***val_Jacobian_gradj, CConfig *config) { }
