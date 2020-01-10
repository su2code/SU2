/*!
 * \file numerics_structure.inl
 * \brief In-Line subroutines of the <i>numerics_structure.hpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#pragma once

inline su2double CNumerics::Determinant_3x3(su2double A00,
                                            su2double A01,
                                            su2double A02,
                                            su2double A10,
                                            su2double A11,
                                            su2double A12,
                                            su2double A20,
                                            su2double A21,
                                            su2double A22) {
  
  return A00*(A11*A22-A12*A21) - A01*(A10*A22-A12*A20) + A02*(A10*A21-A11*A20);
  
}

inline void CNumerics::Compute_Mass_Matrix(CElement *element_container, CConfig *config) { }

inline void CNumerics::Compute_Dead_Load(CElement *element_container, CConfig *config) { }

inline void CNumerics::Compute_Tangent_Matrix(CElement *element_container, CConfig *config) { }

inline void CNumerics::Compute_NodalStress_Term(CElement *element_container, CConfig *config) { }

inline void CNumerics::Compute_Averaged_NodalStress(CElement *element_container, CConfig *config) { }

inline void CNumerics::SetMeshElasticProperties(unsigned long iElem, su2double val_E) { }

inline void CNumerics::Set_DV_Val(unsigned short i_DV, su2double val_DV) { }

inline su2double CNumerics::Get_DV_Val(unsigned short i_DV) const { return 0.0; }

inline void CNumerics::Set_ElectricField(unsigned short i_DV, su2double val_EField) { }

inline void CNumerics::SetMaterial_Properties(unsigned short iVal, su2double val_E, su2double val_Nu) { }

inline void CNumerics::SetMaterial_Density(unsigned short iVal, su2double val_Rho, su2double val_Rho_DL) { }

inline void CNumerics::ComputeResidual(su2double *val_residual, CConfig *config) { }

inline void CNumerics::ComputeResidual(su2double *val_residual_i, su2double *val_residual_j) { }

inline void CNumerics::ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, 
                                   CConfig *config) { }

inline void CNumerics::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                                   su2double **val_JacobianMeanFlow_i, su2double **val_JacobianMeanFlow_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(su2double *val_resconv, su2double *val_resvisc, su2double **val_Jacobian_i, 
                   su2double **val_Jacobian_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, 
                   su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, 
                   su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) { }
              
inline void CNumerics::ComputeResidual(su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, 
                   su2double *val_resvisc_j, su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, 
                   su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) { }
              
inline void CNumerics::ComputeResidual(su2double **val_stiffmatrix_elem, CConfig *config) { }

inline void CNumerics::GetEq_Rxn_Coefficients(su2double **EqnRxnConstants, CConfig *config) { };
                            
inline void CNumerics::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::ComputeResidual_TransLM(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config, su2double &gamma_sep) {}

inline void CNumerics::ComputeResidual_Axisymmetric(su2double *val_residual, CConfig *config) { }

inline void CNumerics::ComputeResidual_Axisymmetric_ad(su2double *val_residual, su2double *val_residuald, CConfig *config) { }

inline void CNumerics::SetJacobian_Axisymmetric(su2double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::ComputeVibRelaxation(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::ComputeChemistry(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::GetKeqConstants(su2double *A, unsigned short val_reaction, CConfig *config) { }

inline void CNumerics::SetUndivided_Laplacian(su2double *val_und_lapl_i, su2double *val_und_lapl_j) {
  Und_Lapl_i = val_und_lapl_i; 
  Und_Lapl_j = val_und_lapl_j; 
}

inline void CNumerics::SetSensor( su2double val_sensor_i, su2double val_sensor_j) {
  Sensor_i = val_sensor_i;
  Sensor_j = val_sensor_j;
}

inline void CNumerics::SetConservative(su2double *val_u_i, su2double *val_u_j) {
  U_i = val_u_i;
  U_j = val_u_j;
}

inline void CNumerics::SetConservative_ZeroOrder(su2double *val_u_i, su2double *val_u_j) {
  UZeroOrder_i = val_u_i;
  UZeroOrder_j = val_u_j;
}

inline void CNumerics::SetPrimitive(su2double *val_v_i, su2double *val_v_j) {
  V_i = val_v_i;
  V_j = val_v_j;
}

inline void CNumerics::SetSecondary(su2double *val_s_i, su2double *val_s_j) {
    S_i = val_s_i;
    S_j = val_s_j;
}

inline void CNumerics::SetConservative(su2double *val_u_0, su2double *val_u_1, su2double *val_u_2) {
  U_0 = val_u_0;
  U_1 = val_u_1;
  U_2 = val_u_2;
}

inline void CNumerics::SetConservative(su2double *val_u_0, su2double *val_u_1, su2double *val_u_2, su2double *val_u_3) {
  U_0 = val_u_0;
  U_1 = val_u_1;
  U_2 = val_u_2;
  U_3 = val_u_3;
}

inline void CNumerics::SetVelocity2_Inf(su2double velocity2) {
  vel2_inf = velocity2;
}

inline void CNumerics::SetVorticity(su2double *val_vorticity_i, su2double *val_vorticity_j) {
  Vorticity_i = val_vorticity_i;
  Vorticity_j = val_vorticity_j;
}

inline void CNumerics::SetStrainMag(su2double val_strainmag_i, su2double val_strainmag_j) {
  StrainMag_i = val_strainmag_i;
  StrainMag_j = val_strainmag_j;
}

inline void CNumerics::SetTimeStep(su2double val_timestep) {TimeStep = val_timestep;}

inline void CNumerics::SetLaminarViscosity(su2double val_lam_viscosity_i, su2double val_lam_viscosity_j) {
  Laminar_Viscosity_i = val_lam_viscosity_i;
  Laminar_Viscosity_j = val_lam_viscosity_j;
}

inline void CNumerics::SetThermalConductivity(su2double val_therm_conductivity_i, su2double val_therm_conductivity_j) {
  Thermal_Conductivity_i = val_therm_conductivity_i;
  Thermal_Conductivity_j = val_therm_conductivity_j;
}

inline void CNumerics::SetThermalConductivity_ve(su2double val_therm_conductivity_ve_i, su2double val_therm_conductivity_ve_j) {
  Thermal_Conductivity_ve_i = val_therm_conductivity_ve_i;
  Thermal_Conductivity_ve_j = val_therm_conductivity_ve_j;
}

inline void CNumerics::SetThermalDiffusivity(su2double val_thermal_diffusivity_i,su2double val_thermal_diffusivity_j) {
  Thermal_Diffusivity_i = val_thermal_diffusivity_i;
  Thermal_Diffusivity_j = val_thermal_diffusivity_j;
}

inline void CNumerics::SetDiffusionCoeff(su2double* val_diffusioncoeff_i, su2double* val_diffusioncoeff_j) {
  Diffusion_Coeff_i = val_diffusioncoeff_i;
  Diffusion_Coeff_j = val_diffusioncoeff_j;
}

inline void CNumerics::SetEddyViscosity(su2double val_eddy_viscosity_i, su2double val_eddy_viscosity_j) {
  Eddy_Viscosity_i = val_eddy_viscosity_i;
  Eddy_Viscosity_j = val_eddy_viscosity_j;
}

inline void CNumerics::SetIntermittency(su2double intermittency_in) { }

inline void CNumerics::SetProduction(su2double val_production) { }

inline void CNumerics::SetDestruction(su2double val_destruction) { }

inline void CNumerics::SetCrossProduction(su2double val_crossproduction) { }

inline su2double CNumerics::GetProduction(void) { return 0; }

inline su2double CNumerics::GetDestruction(void) { return 0; }

inline su2double CNumerics::GetCrossProduction(void) { return 0; }

inline su2double CNumerics::GetGammaBC(void) { return 0.0; }

inline void CNumerics::SetTurbKineticEnergy(su2double val_turb_ke_i, su2double val_turb_ke_j) {
  turb_ke_i = val_turb_ke_i;
  turb_ke_j = val_turb_ke_j;
}

inline void CNumerics::SetDistance(su2double val_dist_i, su2double val_dist_j) {
  dist_i = val_dist_i;
  dist_j = val_dist_j;
}

inline void CNumerics::SetAdjointVar(su2double *val_psi_i, su2double *val_psi_j) {
  Psi_i = val_psi_i;
  Psi_j = val_psi_j;
}

inline void CNumerics::SetAdjointVarGradient(su2double **val_psivar_grad_i, su2double **val_psivar_grad_j) {
  PsiVar_Grad_i = val_psivar_grad_i;
  PsiVar_Grad_j = val_psivar_grad_j;
}

inline void CNumerics::SetTurbVar(su2double *val_turbvar_i, su2double *val_turbvar_j) {
  TurbVar_i = val_turbvar_i;
  TurbVar_j = val_turbvar_j;
}

inline void CNumerics::SetTransVar(su2double *val_transvar_i, su2double *val_transvar_j) {
  TransVar_i = val_transvar_i;
  TransVar_j = val_transvar_j;
}

inline void CNumerics::SetTurbVarGradient(su2double **val_turbvar_grad_i, su2double **val_turbvar_grad_j) {
  TurbVar_Grad_i = val_turbvar_grad_i;
  TurbVar_Grad_j = val_turbvar_grad_j;
}

inline void CNumerics::SetTransVarGradient(su2double **val_transvar_grad_i, su2double **val_transvar_grad_j) {
  TransVar_Grad_i = val_transvar_grad_i;
  TransVar_Grad_j = val_transvar_grad_j;
}

inline void CNumerics::SetPrimVarGradient(su2double **val_primvar_grad_i, su2double **val_primvar_grad_j) {
  PrimVar_Grad_i = val_primvar_grad_i;
  PrimVar_Grad_j = val_primvar_grad_j;
}


inline void CNumerics::SetConsVarGradient(su2double **val_consvar_grad_i, su2double **val_consvar_grad_j) {
  ConsVar_Grad_i = val_consvar_grad_i;
  ConsVar_Grad_j = val_consvar_grad_j;
}

inline void CNumerics::SetConsVarGradient(su2double **val_consvar_grad_0, su2double **val_consvar_grad_1, su2double **val_consvar_grad_2) {
  ConsVar_Grad_0 = val_consvar_grad_0;
  ConsVar_Grad_1 = val_consvar_grad_1;
  ConsVar_Grad_2 = val_consvar_grad_2;
}

inline void CNumerics::SetConsVarGradient(su2double **val_consvar_grad_0, su2double **val_consvar_grad_1, su2double **val_consvar_grad_2, su2double **val_consvar_grad_3) {
  ConsVar_Grad_0 = val_consvar_grad_0;
  ConsVar_Grad_1 = val_consvar_grad_1;
  ConsVar_Grad_2 = val_consvar_grad_2;
  ConsVar_Grad_3 = val_consvar_grad_3;
}

inline void CNumerics::SetConsVarGradient(su2double **val_consvar_grad) {
  ConsVar_Grad = val_consvar_grad;
}

inline void CNumerics::SetCoord(su2double *val_coord_i, su2double *val_coord_j) {
  Coord_i = val_coord_i;
  Coord_j = val_coord_j;
}

inline void CNumerics::SetCoord(su2double *val_coord_0, su2double *val_coord_1, 
                   su2double *val_coord_2) {
  Coord_0 = val_coord_0;
  Coord_1 = val_coord_1;
  Coord_2 = val_coord_2;
}

inline void CNumerics::SetCoord(su2double *val_coord_0, su2double *val_coord_1, 
                   su2double *val_coord_2, su2double *val_coord_3) {
  Coord_0 = val_coord_0;
  Coord_1 = val_coord_1;
  Coord_2 = val_coord_2;
  Coord_3 = val_coord_3;  
}

inline void CNumerics::SetGridVel(su2double *val_gridvel_i, su2double *val_gridvel_j) {
  GridVel_i = val_gridvel_i;
  GridVel_j = val_gridvel_j;
}

inline void CNumerics::SetWindGust(su2double *val_windgust_i, su2double *val_windgust_j) {
  WindGust_i = val_windgust_i;
  WindGust_j = val_windgust_j;
}

inline void CNumerics::SetWindGustDer(su2double *val_windgustder_i, su2double *val_windgustder_j) {
  WindGustDer_i = val_windgustder_i;
  WindGustDer_j = val_windgustder_j;
}

inline void CNumerics::SetPressure(su2double val_pressure_i, su2double val_pressure_j) {
  Pressure_i = val_pressure_i;
  Pressure_j = val_pressure_j;
}

inline void CNumerics::SetDensity(su2double val_densityinc_i, su2double val_densityinc_j) {
  DensityInc_i = val_densityinc_i;
  DensityInc_j = val_densityinc_j;
}

inline void CNumerics::SetBetaInc2(su2double val_betainc2_i, su2double val_betainc2_j) {
  BetaInc2_i = val_betainc2_i;
  BetaInc2_j = val_betainc2_j;
}

inline void CNumerics::SetSoundSpeed(su2double val_soundspeed_i, su2double val_soundspeed_j) {
  SoundSpeed_i = val_soundspeed_i;
  SoundSpeed_j = val_soundspeed_j;
}

inline void CNumerics::SetEnthalpy(su2double val_enthalpy_i, su2double val_enthalpy_j) {
  Enthalpy_i = val_enthalpy_i;
  Enthalpy_j = val_enthalpy_j;
}

inline void CNumerics::SetLambda(su2double val_lambda_i, su2double val_lambda_j) {
  Lambda_i = val_lambda_i;
  Lambda_j = val_lambda_j;
}

inline void CNumerics::SetNeighbor(unsigned short val_neighbor_i, unsigned short val_neighbor_j) {
  Neighbor_i = val_neighbor_i;
  Neighbor_j = val_neighbor_j;
}

inline void CNumerics::SetTurbAdjointVar(su2double *val_turbpsivar_i, su2double *val_turbpsivar_j) {
  TurbPsi_i = val_turbpsivar_i;
  TurbPsi_j = val_turbpsivar_j;
}

inline void CNumerics::SetTurbAdjointGradient(su2double **val_turbpsivar_grad_i, su2double **val_turbpsivar_grad_j) {
  TurbPsi_Grad_i = val_turbpsivar_grad_i;
  TurbPsi_Grad_j = val_turbpsivar_grad_j;
}

inline void CNumerics::SetTemperature(su2double val_temp_i, su2double val_temp_j) {
  Temp_i = val_temp_i;
  Temp_j = val_temp_j;
}

inline void CNumerics::SetAuxVarGrad(su2double *val_auxvargrad_i, su2double *val_auxvargrad_j) {
  AuxVar_Grad_i = val_auxvargrad_i;
  AuxVar_Grad_j = val_auxvargrad_j;
}

inline void CNumerics::SetNormal(su2double *val_normal) { Normal = val_normal; }

inline void CNumerics::SetVolume(su2double val_volume) { Volume = val_volume; }

inline void CNumerics::SetDissipation(su2double diss_i, su2double diss_j) {
  Dissipation_i = diss_i;
  Dissipation_j = diss_j;
}

inline su2double CNumerics::GetDissipation(){
  return Dissipation_ij;
}

inline void CSourcePieceWise_TurbSST::SetF1blending(su2double val_F1_i, su2double val_F1_j) { 
  F1_i = val_F1_i; 
  F1_j = val_F1_j;
}

inline void CSourcePieceWise_TurbSST::SetF2blending(su2double val_F2_i, su2double val_F2_j) { 
  F2_i = val_F2_i; 
  F2_j = val_F2_j;
}

inline void CSourcePieceWise_TurbSST::SetCrossDiff(su2double val_CDkw_i, su2double val_CDkw_j) {
  CDkw_i = val_CDkw_i;
  CDkw_j = val_CDkw_j;
}      

inline void CSourcePieceWise_TurbSA::SetIntermittency(su2double intermittency_in) { intermittency = intermittency_in; }

inline void CSourcePieceWise_TurbSA::SetProduction(su2double val_production) { Production = val_production; }

inline void CSourcePieceWise_TurbSA::SetDestruction(su2double val_destruction) { Destruction = val_destruction; }

inline void CSourcePieceWise_TurbSA::SetCrossProduction(su2double val_crossproduction) { CrossProduction = val_crossproduction; }

inline su2double CSourcePieceWise_TurbSA::GetProduction(void) { return Production; }

inline su2double CSourcePieceWise_TurbSA::GetGammaBC(void) { return gamma_BC; }

inline su2double CSourcePieceWise_TurbSA::GetDestruction(void) { return Destruction; }

inline su2double CSourcePieceWise_TurbSA::GetCrossProduction(void) { return CrossProduction; }

inline void CSourcePieceWise_TurbSA_E::SetIntermittency(su2double intermittency_in) { intermittency = intermittency_in; }

inline void CSourcePieceWise_TurbSA_E::SetProduction(su2double val_production) { Production = val_production; }

inline void CSourcePieceWise_TurbSA_E::SetDestruction(su2double val_destruction) { Destruction = val_destruction; }

inline void CSourcePieceWise_TurbSA_E::SetCrossProduction(su2double val_crossproduction) { CrossProduction = val_crossproduction; }

inline su2double CSourcePieceWise_TurbSA_E::GetProduction(void) { return Production; }

inline su2double CSourcePieceWise_TurbSA_E::GetDestruction(void) { return Destruction; }

inline su2double CSourcePieceWise_TurbSA_E::GetCrossProduction(void) { return CrossProduction; }

inline void CSourcePieceWise_TurbSA_E_COMP::SetIntermittency(su2double intermittency_in) { intermittency = intermittency_in; }

inline void CSourcePieceWise_TurbSA_E_COMP::SetProduction(su2double val_production) { Production = val_production; }

inline void CSourcePieceWise_TurbSA_E_COMP::SetDestruction(su2double val_destruction) { Destruction = val_destruction; }

inline void CSourcePieceWise_TurbSA_E_COMP::SetCrossProduction(su2double val_crossproduction) { CrossProduction = val_crossproduction; }

inline su2double CSourcePieceWise_TurbSA_E_COMP::GetProduction(void) { return Production; }

inline su2double CSourcePieceWise_TurbSA_E_COMP::GetDestruction(void) { return Destruction; }

inline su2double CSourcePieceWise_TurbSA_E_COMP::GetCrossProduction(void) { return CrossProduction; }

inline void CSourcePieceWise_TurbSA_COMP::SetIntermittency(su2double intermittency_in) { intermittency = intermittency_in; }

inline void CSourcePieceWise_TurbSA_COMP::SetProduction(su2double val_production) { Production = val_production; }

inline void CSourcePieceWise_TurbSA_COMP::SetDestruction(su2double val_destruction) { Destruction = val_destruction; }

inline void CSourcePieceWise_TurbSA_COMP::SetCrossProduction(su2double val_crossproduction) { CrossProduction = val_crossproduction; }

inline su2double CSourcePieceWise_TurbSA_COMP::GetProduction(void) { return Production; }

inline su2double CSourcePieceWise_TurbSA_COMP::GetDestruction(void) { return Destruction; }

inline su2double CSourcePieceWise_TurbSA_COMP::GetCrossProduction(void) { return CrossProduction; }

inline void CSourcePieceWise_TurbSA_Neg::SetIntermittency(su2double intermittency_in) { intermittency = intermittency_in; }

inline void CSourcePieceWise_TurbSA_Neg::SetProduction(su2double val_production) { Production = val_production; }

inline void CSourcePieceWise_TurbSA_Neg::SetDestruction(su2double val_destruction) { Destruction = val_destruction; }

inline void CSourcePieceWise_TurbSA_Neg::SetCrossProduction(su2double val_crossproduction) { CrossProduction = val_crossproduction; }

inline su2double CSourcePieceWise_TurbSA_Neg::GetProduction(void) { return Production; }

inline su2double CSourcePieceWise_TurbSA_Neg::GetDestruction(void) { return Destruction; }

inline su2double CSourcePieceWise_TurbSA_Neg::GetCrossProduction(void) { return CrossProduction; }

inline void CNumerics::ComputeResidual(su2double **val_Jacobian_i, su2double *val_Jacobian_mui, su2double ***val_Jacobian_gradi, CConfig *config) { }

inline void CNumerics::ComputeResidual(su2double **val_Jacobian_i, su2double *val_Jacobian_mui, su2double ***val_Jacobian_gradi, 
                  su2double **val_Jacobian_j, su2double *val_Jacobian_muj, su2double ***val_Jacobian_gradj, CConfig *config) { }

inline void CNumerics::SetTauWall(su2double val_tauwall_i, su2double val_tauwall_j) { }

inline void CAvgGrad_Base::SetTauWall(su2double val_tauwall_i, su2double val_tauwall_j) {
  TauWall_i = val_tauwall_i;
  TauWall_j = val_tauwall_j;
}

inline su2double CAvgGrad_Base::GetStressTensor(unsigned short iDim, unsigned short jDim) const {
  return tau[iDim][jDim];
}

inline su2double CAvgGrad_Base::GetHeatFluxVector(unsigned short iDim) const {
  return heat_flux_vector[iDim];
}

inline void CNumerics::SetUsing_UQ(bool val_using_uq) { using_uq = val_using_uq; }
