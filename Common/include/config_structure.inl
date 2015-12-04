/*!
 * \file config_structure.inl
 * \brief In-Line subroutines of the <i>config_structure.hpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 4.0.2 "Cardinal"
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

#pragma once

inline su2double CConfig::GetCFL_AdaptParam(unsigned short val_index) { return CFL_AdaptParam[val_index]; }

inline bool CConfig::GetCFL_Adapt(void) { return CFL_Adapt; }

inline void CConfig::SetInflow_Mach(unsigned short val_imarker, su2double val_fanface_mach) { Inflow_Mach[val_imarker] = val_fanface_mach; }

inline void CConfig::SetInflow_Pressure(unsigned short val_imarker, su2double val_fanface_pressure) { Inflow_Pressure[val_imarker] = val_fanface_pressure; }

inline void CConfig::SetBleed_Temperature(unsigned short val_imarker, su2double val_bleed_temp) { Bleed_Temperature[val_imarker] = val_bleed_temp; }

inline void CConfig::SetBleed_MassFlow(unsigned short val_imarker, su2double val_bleed_massflow) { Bleed_MassFlow[val_imarker] = val_bleed_massflow; }

inline void CConfig::SetBleed_Pressure(unsigned short val_imarker, su2double val_bleed_pressure) { Bleed_Pressure[val_imarker] = val_bleed_pressure; }

inline void CConfig::SetExhaust_Pressure(unsigned short val_imarker, su2double val_exhaust_pressure) { Exhaust_Pressure[val_imarker] = val_exhaust_pressure; }

inline void CConfig::SetExhaust_Temperature(unsigned short val_imarker, su2double val_exhaust_temp) { Exhaust_Temperature[val_imarker] = val_exhaust_temp; }

inline unsigned short CConfig::GetnZone(void) { return nZone; }

inline unsigned short CConfig::GetiZone(void) { return iZone; }

inline unsigned short CConfig::GetKind_SU2(void) { return Kind_SU2; }

inline unsigned short CConfig::GetRef_NonDim(void) { return Ref_NonDim; }

inline void CConfig::SetKind_SU2(unsigned short val_kind_su2) { Kind_SU2 = val_kind_su2 ; }

inline bool CConfig::GetAdjoint(void) { return Adjoint; }

inline bool CConfig::GetViscous(void) { return Viscous; }

inline unsigned long CConfig::GetnExtIter(void) { return nExtIter; }

inline unsigned short CConfig::GetnTimeInstances(void) { return nTimeInstances; }

inline su2double CConfig::GetTimeSpectral_Period(void) { return TimeSpectral_Period; }

inline void CConfig::SetExtIter(unsigned long val_iter) { ExtIter = val_iter; }

inline void CConfig::SetIntIter(unsigned long val_iter) { IntIter = val_iter; }

inline unsigned long CConfig::GetExtIter(void) { return ExtIter; }

inline unsigned long CConfig::GetIntIter(void) { return IntIter; }

inline unsigned long CConfig::GetUnst_nIntIter(void) { return Unst_nIntIter; }

inline long CConfig::GetUnst_RestartIter(void) { return Unst_RestartIter; }

inline long CConfig::GetUnst_AdjointIter(void) { return Unst_AdjointIter; }

inline string CConfig::GetPlaneTag(unsigned short index) { return PlaneTag[index]; }

inline su2double CConfig::GetEA_IntLimit(unsigned short index) { return EA_IntLimit[index]; }

inline su2double CConfig::GetEA_ScaleFactor(void) { return EA_ScaleFactor; }

inline su2double CConfig::GetAdjointLimit(void) { return AdjointLimit; }

inline su2double *CConfig::GetHold_GridFixed_Coord(void) { return Hold_GridFixed_Coord; }

inline su2double *CConfig::GetSubsonic_Engine_Box(void) { return Subsonic_Engine_Box; }

inline su2double CConfig::GetRatioDensity(void) { return RatioDensity; }

inline su2double CConfig::GetFreeSurface_Thickness(void) { return FreeSurface_Thickness; }

inline su2double CConfig::GetFreeSurface_Damping_Coeff(void) { return FreeSurface_Damping_Coeff; }

inline su2double CConfig::GetFreeSurface_Damping_Length(void) { return FreeSurface_Damping_Length; }

inline su2double CConfig::GetFreeSurface_Outlet(void) { return FreeSurface_Outlet; }

inline su2double CConfig::GetRatioViscosity(void) { return RatioViscosity; }

inline unsigned short CConfig::GetAnalytical_Surface(void) { return Analytical_Surface; }

inline unsigned short CConfig::GetAxis_Orientation(void) { return Axis_Orientation; }

inline su2double CConfig::GetDualVol_Power(void) { return DualVol_Power; }

inline bool CConfig::GetExtraOutput(void) { return ExtraOutput; }

inline su2double CConfig::GetRefAreaCoeff(void) { return RefAreaCoeff; }

inline su2double CConfig::GetWaveSpeed(void) { return Wave_Speed; }

inline su2double CConfig::GetThermalDiffusivity(void) { return Thermal_Diffusivity; }

inline su2double CConfig::GetElasticyMod(void) { return ElasticyMod; }

inline unsigned short CConfig::GetElas2D_Formulation(void) { return Kind_2DElasForm; }

inline su2double CConfig::GetPoissonRatio(void) { return PoissonRatio; }

inline su2double CConfig::GetMaterialDensity(void) { return MaterialDensity; }

inline su2double CConfig::GetRefLengthMoment(void) { return RefLengthMoment; }

inline su2double CConfig::GetRefElemLength(void) { return RefElemLength; }

inline su2double CConfig::GetRefSharpEdges(void) { return RefSharpEdges; }

inline su2double CConfig::GetDomainVolume(void) { return DomainVolume; }

inline void CConfig::SetRefAreaCoeff(su2double val_area) { RefAreaCoeff = val_area; }

inline void CConfig::SetDomainVolume(su2double val_volume) { DomainVolume = val_volume; }

inline void CConfig::SetnExtIter(unsigned long val_niter) { nExtIter = val_niter; }

inline su2double CConfig::GetMach(void) { return Mach; }

inline su2double CConfig::GetGamma(void) { return Gamma; }

inline su2double CConfig::GetSection_Location(unsigned short val_var) { return Section_Location[val_var]; }

inline su2double CConfig::GetBulk_Modulus(void) { return Bulk_Modulus; }

inline su2double CConfig::GetArtComp_Factor(void) { return ArtComp_Factor; }

inline su2double CConfig::GetGas_Constant(void) { return Gas_Constant; }

inline su2double CConfig::GetGas_ConstantND(void) { return Gas_ConstantND; }

inline su2double CConfig::GetWallTemperature(void) { return Wall_Temperature; }

inline su2double CConfig::GetFreeSurface_Zero(void) { return FreeSurface_Zero; }

inline su2double CConfig::GetFreeSurface_Depth(void) { return FreeSurface_Depth; }

inline su2double CConfig::GetGas_Constant_Ref(void) { return Gas_Constant_Ref; }

inline su2double CConfig::GetTemperature_FreeStream(void) { return Temperature_FreeStream; }

inline su2double CConfig::GetEnergy_FreeStream(void) { return Energy_FreeStream; }

inline su2double CConfig::GetViscosity_FreeStream(void) { return Viscosity_FreeStream; }

inline su2double CConfig::GetDensity_FreeStream(void) { return Density_FreeStream; }

inline su2double CConfig::GetModVel_FreeStream(void) { return ModVel_FreeStream; }

inline su2double CConfig::GetModVel_FreeStreamND(void) { return ModVel_FreeStreamND; }

inline su2double CConfig::GetPressure_FreeStream(void) { return Pressure_FreeStream; }

inline su2double CConfig::GetTemperature_ve_FreeStream(void) { return Temperature_ve_FreeStream; }

inline su2double CConfig::GetPrandtl_Lam(void) { return Prandtl_Lam; }

inline su2double CConfig::GetPrandtl_Turb(void) { return Prandtl_Turb; }

inline su2double CConfig::GetLength_Ref(void) { return Length_Ref; }

inline su2double CConfig::GetPressure_Ref(void) { return Pressure_Ref; }

inline su2double CConfig::GetTemperature_Ref(void) { return Temperature_Ref; }

inline su2double CConfig::GetDensity_Ref(void) { return Density_Ref; }

inline su2double CConfig::GetVelocity_Ref(void) { return Velocity_Ref; }

inline su2double CConfig::GetEnergy_Ref(void) { return Energy_Ref; }

inline su2double CConfig::GetTime_Ref(void) { return Time_Ref; }

inline su2double CConfig::GetViscosity_Ref(void) { return Viscosity_Ref; }

inline su2double CConfig::GetConductivity_Ref(void) { return Conductivity_Ref; }

inline su2double CConfig::GetOmega_Ref(void) { return Omega_Ref; }

inline su2double CConfig::GetForce_Ref(void) { return Force_Ref; }

inline su2double CConfig::GetPressure_FreeStreamND(void) { return Pressure_FreeStreamND; }

inline su2double CConfig::GetTemperature_FreeStreamND(void) { return Temperature_FreeStreamND; }

inline su2double CConfig::GetDensity_FreeStreamND(void) { return Density_FreeStreamND; }

inline su2double* CConfig::GetVelocity_FreeStreamND(void) { return Velocity_FreeStreamND; }

inline su2double* CConfig::GetVelocity_FreeStream(void) { return Velocity_FreeStream; }

inline su2double CConfig::GetEnergy_FreeStreamND(void) { return Energy_FreeStreamND; }

inline su2double CConfig::GetViscosity_FreeStreamND(void) { return Viscosity_FreeStreamND; }

inline su2double CConfig::GetTke_FreeStreamND(void) { return Tke_FreeStreamND; }

inline su2double CConfig::GetOmega_FreeStreamND(void) { return Omega_FreeStreamND; }

inline su2double CConfig::GetTke_FreeStream(void) { return Tke_FreeStream; }

inline su2double CConfig::GetOmega_FreeStream(void) { return Omega_FreeStream; }

inline su2double CConfig::GetNuFactor_FreeStream(void) { return NuFactor_FreeStream; }

inline su2double CConfig::GetNuFactor_Engine(void) { return NuFactor_Engine; }

inline su2double CConfig::GetIntermittency_FreeStream(void) { return Intermittency_FreeStream; }

inline su2double CConfig::GetTurbulenceIntensity_FreeStream(void) { return TurbulenceIntensity_FreeStream; }

inline su2double CConfig::GetTurb2LamViscRatio_FreeStream(void) {return Turb2LamViscRatio_FreeStream;}

inline su2double* CConfig::GetMassFrac_FreeStream(void) { return MassFrac_FreeStream; }

inline su2double CConfig::GetLength_Reynolds(void) { return Length_Reynolds; }

inline unsigned short CConfig::GetnStartUpIter(void) { return nStartUpIter; }

inline su2double *CConfig::GetRefOriginMoment(unsigned short val_marker) {
    RefOriginMoment[0] = RefOriginMoment_X[val_marker];
    RefOriginMoment[1] = RefOriginMoment_Y[val_marker];
    RefOriginMoment[2] = RefOriginMoment_Z[val_marker];
    return RefOriginMoment;
}

inline su2double CConfig::GetRefOriginMoment_X(unsigned short val_marker) { return RefOriginMoment_X[val_marker]; }

inline su2double CConfig::GetRefOriginMoment_Y(unsigned short val_marker) { return RefOriginMoment_Y[val_marker]; }

inline su2double CConfig::GetRefOriginMoment_Z(unsigned short val_marker) { return RefOriginMoment_Z[val_marker]; }

inline void CConfig::SetRefOriginMoment_X(unsigned short val_marker, su2double val_origin) { RefOriginMoment_X[val_marker] = val_origin; }

inline void CConfig::SetRefOriginMoment_Y(unsigned short val_marker, su2double val_origin) { RefOriginMoment_Y[val_marker] = val_origin; }

inline void CConfig::SetRefOriginMoment_Z(unsigned short val_marker, su2double val_origin) { RefOriginMoment_Z[val_marker] = val_origin; }

inline su2double CConfig::GetChargeCoeff(void) { return ChargeCoeff; }

inline su2double CConfig::GetLimiterCoeff(void) { return LimiterCoeff; }

inline unsigned long CConfig::GetLimiterIter(void) { return LimiterIter; }

inline su2double CConfig::GetSharpEdgesCoeff(void) { return SharpEdgesCoeff; }

inline su2double CConfig::GetReynolds(void) { return Reynolds; }

inline su2double CConfig::GetFroude(void) { return Froude; }

inline void CConfig::SetPressure_FreeStreamND(su2double val_pressure_freestreamnd) { Pressure_FreeStreamND = val_pressure_freestreamnd; }

inline void CConfig::SetPressure_FreeStream(su2double val_pressure_freestream) { Pressure_FreeStream = val_pressure_freestream; }

inline void CConfig::SetDensity_FreeStreamND(su2double val_density_freestreamnd) { Density_FreeStreamND = val_density_freestreamnd; }

inline void CConfig::SetDensity_FreeStream(su2double val_density_freestream) { Density_FreeStream = val_density_freestream; }

inline void CConfig::SetViscosity_FreeStream(su2double val_viscosity_freestream) { Viscosity_FreeStream = val_viscosity_freestream; }

inline void CConfig::SetModVel_FreeStream(su2double val_modvel_freestream) { ModVel_FreeStream = val_modvel_freestream; }

inline void CConfig::SetModVel_FreeStreamND(su2double val_modvel_freestreamnd) { ModVel_FreeStreamND = val_modvel_freestreamnd; }

inline void CConfig::SetTemperature_FreeStream(su2double val_temperature_freestream) { Temperature_FreeStream = val_temperature_freestream; }

inline void CConfig::SetTemperature_FreeStreamND(su2double val_temperature_freestreamnd) { Temperature_FreeStreamND = val_temperature_freestreamnd; }

inline void CConfig::SetGas_ConstantND(su2double val_gas_constantnd) { Gas_ConstantND = val_gas_constantnd; }

inline void CConfig::SetVelocity_FreeStreamND(su2double val_velocity_freestreamnd, unsigned short val_dim) { Velocity_FreeStreamND[val_dim] = val_velocity_freestreamnd; }

inline void CConfig::SetViscosity_FreeStreamND(su2double val_viscosity_freestreamnd) { Viscosity_FreeStreamND = val_viscosity_freestreamnd; }

inline void CConfig::SetTke_FreeStreamND(su2double val_tke_freestreamnd) { Tke_FreeStreamND = val_tke_freestreamnd; }

inline void CConfig::SetOmega_FreeStreamND(su2double val_omega_freestreamnd) { Omega_FreeStreamND = val_omega_freestreamnd; }

inline void CConfig::SetTke_FreeStream(su2double val_tke_freestream) { Tke_FreeStream = val_tke_freestream; }

inline void CConfig::SetOmega_FreeStream(su2double val_omega_freestream) { Omega_FreeStream = val_omega_freestream; }

inline void CConfig::SetEnergy_FreeStreamND(su2double val_energy_freestreamnd) { Energy_FreeStreamND = val_energy_freestreamnd; }

inline void CConfig::SetEnergy_FreeStream(su2double val_energy_freestream) { Energy_FreeStream = val_energy_freestream; }

inline void CConfig::SetTotal_UnstTimeND(su2double val_total_unsttimend) { Total_UnstTimeND = val_total_unsttimend; }

inline void CConfig::SetFroude(su2double val_froude) { Froude = val_froude; }

inline void CConfig::SetReynolds(su2double val_reynolds) { Reynolds = val_reynolds; }

inline void CConfig::SetMach(su2double val_mach) { Mach = val_mach; }

inline void CConfig::SetLength_Ref(su2double val_length_ref) { Length_Ref = val_length_ref; }

inline void CConfig::SetVelocity_Ref(su2double val_velocity_ref) { Velocity_Ref = val_velocity_ref; }

inline void CConfig::SetPressure_Ref(su2double val_pressure_ref) { Pressure_Ref = val_pressure_ref; }

inline void CConfig::SetDensity_Ref(su2double val_density_ref) { Density_Ref = val_density_ref; }

inline void CConfig::SetTemperature_Ref(su2double val_temperature_ref) { Temperature_Ref = val_temperature_ref; }

inline void CConfig::SetTime_Ref(su2double val_time_ref) { Time_Ref = val_time_ref; }

inline void CConfig::SetOmega_Ref(su2double val_omega_ref) { Omega_Ref = val_omega_ref; }

inline void CConfig::SetForce_Ref(su2double val_force_ref) { Force_Ref = val_force_ref; }

inline void CConfig::SetGas_Constant_Ref(su2double val_gas_constant_ref) { Gas_Constant_Ref = val_gas_constant_ref; }

inline void CConfig::SetGas_Constant(su2double val_gas_constant) { Gas_Constant = val_gas_constant; }

inline void CConfig::SetViscosity_Ref(su2double val_viscosity_ref) { Viscosity_Ref = val_viscosity_ref; }

inline void CConfig::SetConductivity_Ref(su2double val_conductivity_ref) { Conductivity_Ref = val_conductivity_ref; }

inline void CConfig::SetEnergy_Ref(su2double val_energy_ref) { Energy_Ref = val_energy_ref; }

inline su2double CConfig::GetAoA(void) { return AoA; }

inline void CConfig::SetAoA(su2double val_AoA) { AoA = val_AoA; }

inline void CConfig::SetAoS(su2double val_AoS) { AoS = val_AoS; }

inline su2double CConfig::GetAoS(void) { return AoS; }

inline unsigned short CConfig::GetnMGLevels(void) { return nMGLevels; }

inline void CConfig::SetMGLevels(unsigned short val_nMGLevels) { nMGLevels = val_nMGLevels; }

inline unsigned short CConfig::GetFinestMesh(void) { return FinestMesh; }

inline void CConfig::SetFinestMesh(unsigned short val_finestmesh) { FinestMesh = val_finestmesh; }

inline void CConfig::SubtractFinestMesh(void) { FinestMesh = FinestMesh-1; }

inline unsigned short CConfig::GetDesign_Variable(unsigned short val_dv) { return Design_Variable[val_dv]; }

inline unsigned short CConfig::GetConvCriteria(void) { return ConvCriteria; }

inline unsigned short CConfig::GetMGCycle(void) { return MGCycle; }

inline unsigned short CConfig::GetGeometryMode(void) { return GeometryMode; }

inline su2double CConfig::GetCFL(unsigned short val_mesh) {	return CFL[val_mesh]; }

inline void CConfig::SetCFL(unsigned short val_mesh, su2double val_cfl) { CFL[val_mesh] = val_cfl; }

inline su2double CConfig::GetUnst_CFL(void) {	return Unst_CFL; }

inline su2double CConfig::GetMax_DeltaTime(void) {	return Max_DeltaTime; }

inline su2double CConfig::GetParamDV(unsigned short val_dv, unsigned short val_param) {	return ParamDV[val_dv][val_param]; }

inline su2double CConfig::GetCoordFFDBox(unsigned short val_ffd, unsigned short val_coord) {	return CoordFFDBox[val_ffd][val_coord]; }

inline unsigned short CConfig::GetDegreeFFDBox(unsigned short val_ffd, unsigned short val_degree) {	return DegreeFFDBox[val_ffd][val_degree]; }

inline string CConfig::GetFFDTag(unsigned short val_dv) {	return FFDTag[val_dv]; }

inline string CConfig::GetTagFFDBox(unsigned short val_ffd) {	return TagFFDBox[val_ffd]; }

inline unsigned short CConfig::GetnDV(void) {	return nDV; }

inline unsigned short CConfig::GetnFFDBox(void) {	return nFFDBox; }

inline unsigned short CConfig::GetFFD_Continuity(void) { return FFD_Continuity; }

inline unsigned short CConfig::GetnRKStep(void) { return nRKStep; }

inline su2double CConfig::Get_Alpha_RKStep(unsigned short val_step) { return RK_Alpha_Step[val_step]; }

inline unsigned short CConfig::GetMG_PreSmooth(unsigned short val_mesh) {	
	if (nMG_PreSmooth == 0) return 1;
	else return MG_PreSmooth[val_mesh]; 
}

inline unsigned short CConfig::GetMG_PostSmooth(unsigned short val_mesh) { 
	if (nMG_PostSmooth == 0) return 0;
	else return MG_PostSmooth[val_mesh];
}

inline unsigned short CConfig::GetMG_CorrecSmooth(unsigned short val_mesh) { 
	if (nMG_CorrecSmooth == 0) return 0;
	else return MG_CorrecSmooth[val_mesh]; 
}

inline unsigned long CConfig::GetWrt_Sol_Freq(void) { return Wrt_Sol_Freq; }

inline unsigned long CConfig::GetWrt_Sol_Freq_DualTime(void) { return Wrt_Sol_Freq_DualTime; }

inline unsigned long CConfig::GetWrt_Con_Freq(void) { return Wrt_Con_Freq; }

inline unsigned long CConfig::GetWrt_Con_Freq_DualTime(void) { return Wrt_Con_Freq_DualTime; }

inline bool CConfig::GetWrt_Unsteady(void) { return Wrt_Unsteady; }

inline bool CConfig::GetLowFidelitySim(void) { return LowFidelitySim; }

inline bool CConfig::GetIonization(void) { return ionization; }

inline unsigned short CConfig::GetKind_Solver(void) { return Kind_Solver; }

inline void CConfig::SetKind_Solver(unsigned short val_solver) { Kind_Solver = val_solver; }

inline unsigned short CConfig::GetKind_Regime(void) { return Kind_Regime; }

inline unsigned short CConfig::GetSystemMeasurements(void) { return SystemMeasurements; }

inline unsigned short CConfig::GetKind_GasModel(void) { return Kind_GasModel; }

inline unsigned short CConfig::GetKind_FluidModel(void) { return Kind_FluidModel; }

inline unsigned short CConfig::GetKind_FreeStreamOption(void) {return Kind_FreeStreamOption; } 

inline unsigned short CConfig::GetKind_InitOption(void) {return Kind_InitOption; }

inline su2double CConfig::GetPressure_Critical(void) { return Pressure_Critical; }

inline su2double CConfig::GetTemperature_Critical(void) { return Temperature_Critical; }

inline su2double CConfig::GetAcentric_Factor(void) { return Acentric_Factor; }

inline unsigned short CConfig::GetKind_ViscosityModel(void) { return Kind_ViscosityModel; }

inline unsigned short CConfig::GetKind_ConductivityModel(void) { return Kind_ConductivityModel; }

inline su2double CConfig::GetMu_ConstantND(void) { return Mu_ConstantND; }

inline su2double CConfig::GetKt_ConstantND(void) { return Kt_ConstantND; }

inline su2double CConfig::GetMu_RefND(void) { return Mu_RefND; }

inline su2double CConfig::GetMu_Temperature_RefND(void) { return Mu_Temperature_RefND; }

inline su2double CConfig::GetMu_SND(void) { return Mu_SND; }

inline void CConfig::SetMu_ConstantND(su2double mu_const) { Mu_ConstantND = mu_const; }

inline void CConfig::SetMu_RefND(su2double mu_ref) { Mu_RefND = mu_ref; }

inline void CConfig::SetMu_Temperature_RefND(su2double mu_Tref) {Mu_Temperature_RefND = mu_Tref; }

inline void CConfig::SetMu_SND(su2double mu_s) {Mu_SND = mu_s; }

inline void CConfig::SetKt_ConstantND(su2double kt_const) { Kt_ConstantND = kt_const; }

inline unsigned short CConfig::GetKind_GridMovement(unsigned short val_iZone) { return Kind_GridMovement[val_iZone]; }

inline void CConfig::SetKind_GridMovement(unsigned short val_iZone, unsigned short motion_Type) { Kind_GridMovement[val_iZone] = motion_Type; }

inline su2double CConfig::GetMach_Motion(void) { return Mach_Motion; }

inline su2double CConfig::GetMotion_Origin_X(unsigned short val_iZone) { return Motion_Origin_X[val_iZone]; }

inline su2double CConfig::GetMotion_Origin_Y(unsigned short val_iZone) { return Motion_Origin_Y[val_iZone]; }

inline su2double CConfig::GetMotion_Origin_Z(unsigned short val_iZone) { return Motion_Origin_Z[val_iZone]; }

inline void CConfig::SetMotion_Origin_X(unsigned short val_iZone, su2double val_origin) { Motion_Origin_X[val_iZone] = val_origin; }

inline void CConfig::SetMotion_Origin_Y(unsigned short val_iZone, su2double val_origin) { Motion_Origin_Y[val_iZone] = val_origin; }

inline void CConfig::SetMotion_Origin_Z(unsigned short val_iZone, su2double val_origin) { Motion_Origin_Z[val_iZone] = val_origin; }

inline su2double CConfig::GetTranslation_Rate_X(unsigned short val_iZone) { return  Translation_Rate_X[val_iZone]; }

inline su2double CConfig::GetTranslation_Rate_Y(unsigned short val_iZone) { return  Translation_Rate_Y[val_iZone]; }

inline su2double CConfig::GetTranslation_Rate_Z(unsigned short val_iZone) { return  Translation_Rate_Z[val_iZone]; }

inline su2double CConfig::GetRotation_Rate_X(unsigned short val_iZone) { return  Rotation_Rate_X[val_iZone]; }

inline su2double CConfig::GetRotation_Rate_Y(unsigned short val_iZone) { return  Rotation_Rate_Y[val_iZone]; }

inline su2double CConfig::GetRotation_Rate_Z(unsigned short val_iZone) { return  Rotation_Rate_Z[val_iZone]; }

inline su2double CConfig::GetPitching_Omega_X(unsigned short val_iZone) { return  Pitching_Omega_X[val_iZone]; }

inline su2double CConfig::GetPitching_Omega_Y(unsigned short val_iZone) { return  Pitching_Omega_Y[val_iZone]; }

inline su2double CConfig::GetPitching_Omega_Z(unsigned short val_iZone) { return  Pitching_Omega_Z[val_iZone]; }

inline su2double CConfig::GetPitching_Ampl_X(unsigned short val_iZone) { return  Pitching_Ampl_X[val_iZone]; }

inline su2double CConfig::GetPitching_Ampl_Y(unsigned short val_iZone) { return  Pitching_Ampl_Y[val_iZone]; }

inline su2double CConfig::GetPitching_Ampl_Z(unsigned short val_iZone) { return  Pitching_Ampl_Z[val_iZone]; }

inline su2double CConfig::GetPitching_Phase_X(unsigned short val_iZone) { return  Pitching_Phase_X[val_iZone]; }

inline su2double CConfig::GetPitching_Phase_Y(unsigned short val_iZone) { return  Pitching_Phase_Y[val_iZone]; }

inline su2double CConfig::GetPitching_Phase_Z(unsigned short val_iZone) { return  Pitching_Phase_Z[val_iZone]; }

inline su2double CConfig::GetPlunging_Omega_X(unsigned short val_iZone) { return  Plunging_Omega_X[val_iZone]; }

inline su2double CConfig::GetPlunging_Omega_Y(unsigned short val_iZone) { return  Plunging_Omega_Y[val_iZone]; }

inline su2double CConfig::GetPlunging_Omega_Z(unsigned short val_iZone) { return  Plunging_Omega_Z[val_iZone]; }

inline su2double CConfig::GetPlunging_Ampl_X(unsigned short val_iZone) { return  Plunging_Ampl_X[val_iZone]; }

inline su2double CConfig::GetPlunging_Ampl_Y(unsigned short val_iZone) { return  Plunging_Ampl_Y[val_iZone]; }

inline su2double CConfig::GetPlunging_Ampl_Z(unsigned short val_iZone) { return  Plunging_Ampl_Z[val_iZone]; }

inline unsigned short CConfig::GetMoveMotion_Origin(unsigned short val_marker) {return MoveMotion_Origin[val_marker]; }

inline su2double CConfig::GetminTurkelBeta() { return  Min_Beta_RoeTurkel; }

inline su2double CConfig::GetmaxTurkelBeta() { return  Max_Beta_RoeTurkel; }

inline unsigned short CConfig::GetKind_Gradient_Method(void) { return Kind_Gradient_Method; }

inline unsigned short CConfig::GetKind_Linear_Solver(void) { return Kind_Linear_Solver; }

inline unsigned short CConfig::GetDeform_Linear_Solver(void) { return Deform_Linear_Solver; }

inline unsigned short CConfig::GetKind_Linear_Solver_Prec(void) { return Kind_Linear_Solver_Prec; }

inline void CConfig::SetKind_Linear_Solver_Prec(unsigned short val_kind_prec) { Kind_Linear_Solver_Prec = val_kind_prec; }

inline su2double CConfig::GetLinear_Solver_Error(void) { return Linear_Solver_Error; }

inline unsigned long CConfig::GetLinear_Solver_Iter(void) { return Linear_Solver_Iter; }

inline unsigned long CConfig::GetLinear_Solver_Restart_Frequency(void) { return Linear_Solver_Restart_Frequency; }

inline su2double CConfig::GetRelaxation_Factor_Flow(void) { return Relaxation_Factor_Flow; }

inline su2double CConfig::GetRelaxation_Factor_AdjFlow(void) { return Relaxation_Factor_AdjFlow; }

inline su2double CConfig::GetRelaxation_Factor_Turb(void) { return Relaxation_Factor_Turb; }

inline su2double CConfig::GetRoe_Kappa(void) { return Roe_Kappa; }

inline unsigned short CConfig::GetKind_AdjTurb_Linear_Solver(void) { return Kind_AdjTurb_Linear_Solver; }

inline unsigned short CConfig::GetKind_AdjTurb_Linear_Prec(void) { return Kind_AdjTurb_Linear_Prec; }

inline unsigned short CConfig::GetKind_DiscAdj_Linear_Solver(void) { return Kind_DiscAdj_Linear_Solver; }

inline unsigned short CConfig::GetKind_DiscAdj_Linear_Prec(void) { return Kind_DiscAdj_Linear_Prec; }

inline void CConfig::SetKind_AdjTurb_Linear_Prec(unsigned short val_kind_prec) { Kind_AdjTurb_Linear_Prec = val_kind_prec; }

inline su2double CConfig::GetAdjTurb_Linear_Error(void) { return AdjTurb_Linear_Error; }

inline su2double CConfig::GetEntropyFix_Coeff(void) { return EntropyFix_Coeff; }

inline unsigned short CConfig::GetAdjTurb_Linear_Iter(void) { return AdjTurb_Linear_Iter; }

inline su2double CConfig::GetCFLRedCoeff_AdjTurb(void) { return CFLRedCoeff_AdjTurb; }

inline unsigned long CConfig::GetGridDef_Linear_Iter(void) { return GridDef_Linear_Iter; }

inline unsigned long CConfig::GetGridDef_Nonlinear_Iter(void) { return GridDef_Nonlinear_Iter; }

inline bool CConfig::GetDeform_Output(void) { return Deform_Output; }

inline su2double CConfig::GetDeform_Tol_Factor(void) { return Deform_Tol_Factor; }

inline su2double CConfig::GetDeform_ElasticityMod(void) { return Deform_ElasticityMod; }

inline su2double CConfig::GetDeform_PoissonRatio(void) { return Deform_PoissonRatio; }

inline unsigned short CConfig::GetDeform_Stiffness_Type(void) { return Deform_Stiffness_Type; }

inline bool CConfig::GetVisualize_Deformation(void) { return Visualize_Deformation; }

inline unsigned short CConfig::GetKind_Adaptation(void) { return Kind_Adaptation; }

inline su2double CConfig::GetNew_Elem_Adapt(void) { return New_Elem_Adapt; }

inline unsigned short CConfig::GetKind_TimeIntScheme(void) { return Kind_TimeNumScheme; }

inline unsigned short CConfig::GetKind_ConvNumScheme(void) { return Kind_ConvNumScheme; }

inline unsigned short CConfig::GetKind_Centered(void) { return Kind_Centered; }

inline unsigned short CConfig::GetKind_Upwind(void) { return Kind_Upwind; }

inline unsigned short CConfig::GetSpatialOrder(void) { return SpatialOrder; }

inline unsigned short CConfig::GetSpatialOrder_Flow(void) { return SpatialOrder_Flow; }

inline unsigned short CConfig::GetSpatialOrder_Turb(void) { return SpatialOrder_Turb; }

inline unsigned short CConfig::GetSpatialOrder_AdjLevelSet(void) { return SpatialOrder_AdjLevelSet; }

inline unsigned short CConfig::GetSpatialOrder_AdjFlow(void) { return SpatialOrder_AdjFlow; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Flow(void) { return Kind_TimeIntScheme_Flow; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Wave(void) { return Kind_TimeIntScheme_Wave; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Heat(void) { return Kind_TimeIntScheme_Heat; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Poisson(void) { return Kind_TimeIntScheme_Poisson; }

inline unsigned short CConfig::GetKind_TimeIntScheme_FEA(void) { return Kind_TimeIntScheme_FEA; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Flow(void) { return Kind_ConvNumScheme_Flow; }

inline unsigned short CConfig::GetKind_ConvNumScheme_AdjLevelSet(void) { return Kind_ConvNumScheme_AdjLevelSet; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Template(void) { return Kind_ConvNumScheme_Template; }

inline unsigned short CConfig::GetKind_Centered_Flow(void) { return Kind_Centered_Flow; }

inline unsigned short CConfig::GetKind_Centered_AdjLevelSet(void) { return Kind_Centered_AdjLevelSet; }

inline unsigned short CConfig::GetKind_SlopeLimit(void) { return Kind_SlopeLimit; }

inline unsigned short CConfig::GetKind_SlopeLimit_Flow(void) { return Kind_SlopeLimit_Flow; }

inline unsigned short CConfig::GetKind_SlopeLimit_Turb(void) { return Kind_SlopeLimit_Turb; }

inline unsigned short CConfig::GetKind_SlopeLimit_AdjLevelSet(void) { return Kind_SlopeLimit_AdjLevelSet; }

inline unsigned short CConfig::GetKind_SlopeLimit_AdjTurb(void) { return Kind_SlopeLimit_AdjTurb; }

inline unsigned short CConfig::GetKind_SlopeLimit_AdjFlow(void) { return Kind_SlopeLimit_AdjFlow; }

inline unsigned short CConfig::GetKind_Upwind_Flow(void) { return Kind_Upwind_Flow; }

inline unsigned short CConfig::GetKind_Upwind_AdjLevelSet(void) { return Kind_Upwind_AdjLevelSet; }

inline su2double CConfig::GetKappa_1st_Flow(void) { return Kappa_1st_Flow; }

inline su2double CConfig::GetKappa_2nd_Flow(void) { return Kappa_2nd_Flow; }

inline su2double CConfig::GetKappa_4th_Flow(void) { return Kappa_4th_Flow; }

inline unsigned short CConfig::GetKind_TimeIntScheme_AdjFlow(void) { return Kind_TimeIntScheme_AdjFlow; }

inline unsigned short CConfig::GetKind_ConvNumScheme_AdjFlow(void) { return Kind_ConvNumScheme_AdjFlow; }

inline unsigned short CConfig::GetKind_Centered_AdjFlow(void) { return Kind_Centered_AdjFlow; }

inline unsigned short CConfig::GetKind_Upwind_AdjFlow(void) { return Kind_Upwind_AdjFlow; }

inline su2double CConfig::GetKappa_1st_AdjFlow(void) { return Kappa_1st_AdjFlow; }

inline su2double CConfig::GetKappa_2nd_AdjFlow(void) { return Kappa_2nd_AdjFlow; }

inline su2double CConfig::GetKappa_4th_AdjFlow(void) { return Kappa_4th_AdjFlow; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Turb(void) { return Kind_TimeIntScheme_Turb; }

inline unsigned short CConfig::GetKind_TimeIntScheme_AdjLevelSet(void) { return Kind_TimeIntScheme_AdjLevelSet; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Turb(void) { return Kind_ConvNumScheme_Turb; }

inline unsigned short CConfig::GetKind_Centered_Turb(void) { return Kind_Centered_Turb; }

inline unsigned short CConfig::GetKind_Upwind_Turb(void) {	return Kind_Upwind_Turb; }

inline unsigned short CConfig::GetKind_TimeIntScheme_AdjTurb(void) { return Kind_TimeIntScheme_AdjTurb; }

inline unsigned short CConfig::GetKind_ConvNumScheme_AdjTurb(void) { return Kind_ConvNumScheme_AdjTurb; }

inline unsigned short CConfig::GetKind_Centered_AdjTurb(void) { return Kind_Centered_AdjTurb; }

inline unsigned short CConfig::GetKind_Upwind_AdjTurb(void) { return Kind_Upwind_AdjTurb; }

inline unsigned short CConfig::GetKind_Inlet(void) { return Kind_Inlet; }

inline unsigned short CConfig::GetKind_MixingProcess(void) { return Kind_MixingProcess; }

inline bool CConfig::GetBoolMixingPlane(void) { return (nMarker_MixBound !=0);}

inline bool CConfig::GetBoolTurboPerf(void) { return (nMarker_TurboPerf !=0);}

inline string CConfig::GetMarker_MixingPlane_Bound(unsigned short index) {return Marker_MixBound[index];}

inline string CConfig::GetMarker_MixingPlane_Donor(unsigned short index) {return Marker_MixDonor[index];}

inline unsigned short CConfig::Get_nMarkerMixingPlane(void) { return nMarker_MixBound;}

inline unsigned short CConfig::Get_nMarkerTurboPerf(void) { return nMarker_TurboPerf;}

inline string CConfig::GetMarker_TurboPerf_BoundIn(unsigned short index) {return Marker_TurboBoundIn[index];}

inline string CConfig::GetMarker_TurboPerf_BoundOut(unsigned short index) {return Marker_TurboBoundOut[index];}

inline unsigned short CConfig::GetKind_TurboPerf(unsigned short index) {return Kind_TurboPerformance[index];}

inline unsigned short CConfig::GetnSections(void) { return nSections; }

inline unsigned short CConfig::GetnVolSections(void) { return nVolSections; }

inline void CConfig::SetKind_TimeIntScheme(unsigned short val_kind_timeintscheme) { Kind_TimeNumScheme = val_kind_timeintscheme; }

inline unsigned short CConfig::GetKind_ObjFunc(void) {return Kind_ObjFunc; }

inline su2double CConfig::GetCoeff_ObjChainRule(unsigned short iVar) {return Obj_ChainRuleCoeff[iVar]; }

inline unsigned short CConfig::GetKind_SensSmooth(void) {return Kind_SensSmooth; }

inline unsigned short CConfig::GetUnsteady_Simulation(void) { return Unsteady_Simulation; }

inline bool CConfig::GetRestart(void) {	return Restart; }

inline bool CConfig::GetRestart_Flow(void) { return Restart_Flow; }

inline bool CConfig::GetEquivArea(void) { return EquivArea; }

inline bool CConfig::GetInvDesign_Cp(void) { return InvDesign_Cp; }

inline bool CConfig::GetInvDesign_HeatFlux(void) { return InvDesign_HeatFlux; }

inline void CConfig::SetnMarker_All(unsigned short val_nmarker) { nMarker_All = val_nmarker; }

inline string CConfig::GetMarker_All_TagBound(unsigned short val_marker) { return Marker_All_TagBound[val_marker]; }

inline string CConfig::GetMarker_ActDisk_Inlet(unsigned short val_marker) { return Marker_ActDisk_Inlet[val_marker]; }

inline string CConfig::GetMarker_ActDisk_Outlet(unsigned short val_marker) { return Marker_ActDisk_Outlet[val_marker]; }

inline string CConfig::GetMarker_EngineInflow(unsigned short val_marker) { return Marker_EngineInflow[val_marker]; }

inline string CConfig::GetMarker_EngineBleed(unsigned short val_marker) { return Marker_EngineBleed[val_marker]; }

inline string CConfig::GetMarker_EngineExhaust(unsigned short val_marker) { return Marker_EngineExhaust[val_marker]; }

inline string CConfig::GetMarker_Monitoring(unsigned short val_marker) { return Marker_Monitoring[val_marker]; }

inline string CConfig::GetMarker_Moving(unsigned short val_marker) { return Marker_Moving[val_marker]; }

inline short CConfig::GetMarker_All_TagBound(string val_tag) {
	for (unsigned short iMarker = 0; iMarker < nMarker_All; iMarker++) {
		if (val_tag == Marker_All_TagBound[iMarker])
		return iMarker; 
	}
	return -1;
}

inline unsigned short CConfig::GetMarker_All_KindBC(unsigned short val_marker) { return Marker_All_KindBC[val_marker]; }

inline void CConfig::SetMarker_All_KindBC(unsigned short val_marker, unsigned short val_boundary) { Marker_All_KindBC[val_marker] = val_boundary; }

inline void CConfig::SetMarker_All_TagBound(unsigned short val_marker, string val_index) { Marker_All_TagBound[val_marker] = val_index; }

inline void CConfig::SetMarker_All_Monitoring(unsigned short val_marker, unsigned short val_monitoring) { Marker_All_Monitoring[val_marker] = val_monitoring; }

inline void CConfig::SetMarker_All_GeoEval(unsigned short val_marker, unsigned short val_geoeval) { Marker_All_GeoEval[val_marker] = val_geoeval; }

inline void CConfig::SetMarker_All_Designing(unsigned short val_marker, unsigned short val_designing) { Marker_All_Designing[val_marker] = val_designing; }

inline void CConfig::SetMarker_All_Plotting(unsigned short val_marker, unsigned short val_plotting) { Marker_All_Plotting[val_marker] = val_plotting; }

inline void CConfig::SetMarker_All_FSIinterface(unsigned short val_marker, unsigned short val_fsiinterface) { Marker_All_FSIinterface[val_marker] = val_fsiinterface; }

inline void CConfig::SetMarker_All_TurboPerformance(unsigned short val_marker, unsigned short val_turboperf) { Marker_All_TurboPerformance[val_marker] = val_turboperf; }

inline void CConfig::SetMarker_All_TurboPerformanceFlag(unsigned short val_marker, unsigned short val_turboperflag) { Marker_All_TurboPerformanceFlag[val_marker] = val_turboperflag; }

inline void CConfig::SetMarker_All_DV(unsigned short val_marker, unsigned short val_DV) { Marker_All_DV[val_marker] = val_DV; }

inline void CConfig::SetMarker_All_Moving(unsigned short val_marker, unsigned short val_moving) { Marker_All_Moving[val_marker] = val_moving; }

inline void CConfig::SetMarker_All_PerBound(unsigned short val_marker, short val_perbound) { Marker_All_PerBound[val_marker] = val_perbound; }

inline short CConfig::GetMarker_All_PerBound(unsigned short val_marker) { return Marker_All_PerBound[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Monitoring(unsigned short val_marker) { return Marker_All_Monitoring[val_marker]; }

inline void CConfig::SetMarker_All_Out_1D(unsigned short val_marker, unsigned short val_boundary) { Marker_All_Out_1D[val_marker] = val_boundary; }

inline unsigned short CConfig::GetMarker_All_Out_1D(unsigned short val_marker) { return Marker_All_Out_1D[val_marker]; }

inline unsigned short CConfig::GetMarker_All_GeoEval(unsigned short val_marker) { return Marker_All_GeoEval[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Designing(unsigned short val_marker) { return Marker_All_Designing[val_marker]; }

inline short CConfig::GetMarker_All_SendRecv(unsigned short val_marker) { return Marker_All_SendRecv[val_marker]; }

inline void CConfig::SetMarker_All_SendRecv(unsigned short val_marker, short val_index) { Marker_All_SendRecv[val_marker] = val_index; }

inline unsigned short CConfig::GetMarker_All_Plotting(unsigned short val_marker) { return Marker_All_Plotting[val_marker]; }

inline unsigned short CConfig::GetMarker_All_FSIinterface(unsigned short val_marker) { return Marker_All_FSIinterface[val_marker]; }

inline unsigned short CConfig::GetMarker_All_TurboPerformance(unsigned short val_marker) { return Marker_All_TurboPerformance[val_marker]; }

inline unsigned short CConfig::GetMarker_All_TurboPerformanceFlag(unsigned short val_marker) { return Marker_All_TurboPerformanceFlag[val_marker]; }

inline unsigned short CConfig::GetMarker_n_FSIinterface(void) { return nMarker_FSIinterface; }

inline unsigned short CConfig::GetMarker_All_DV(unsigned short val_marker) { return Marker_All_DV[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Moving(unsigned short val_marker) { return Marker_All_Moving[val_marker]; }

inline unsigned short CConfig::GetnMarker_All(void) { return nMarker_All; }

inline unsigned short CConfig::GetnMarker_Max(void) { return nMarker_Max; }

inline unsigned short CConfig::GetnMarker_EngineInflow(void) {	return nMarker_EngineInflow; }

inline unsigned short CConfig::GetnMarker_EngineBleed(void) {	return nMarker_EngineBleed; }

inline unsigned short CConfig::GetnMarker_EngineExhaust(void) { return nMarker_EngineExhaust; }

inline unsigned short CConfig::GetnMarker_InterfaceBound(void) { return nMarker_InterfaceBound; }

inline unsigned short CConfig::GetnMarker_Monitoring(void) { return nMarker_Monitoring; }

inline unsigned short CConfig::GetnMarker_Out_1D(void) { return nMarker_Out_1D; }

inline unsigned short CConfig::GetnMarker_Moving(void) { return nMarker_Moving; }

inline unsigned short CConfig::GetnMarker_NearFieldBound(void) { return nMarker_NearFieldBound; }

inline unsigned short CConfig::GetnMarker_ActDisk_Inlet(void) { return nMarker_ActDisk_Inlet; }

inline unsigned short CConfig::GetnMarker_ActDisk_Outlet(void) { return nMarker_ActDisk_Outlet; }

inline string CConfig::GetMesh_FileName(void) { return Mesh_FileName; }

inline string CConfig::GetMesh_Out_FileName(void) { return Mesh_Out_FileName; }

inline unsigned short CConfig::GetMesh_FileFormat(void) { return Mesh_FileFormat; }

inline unsigned short CConfig::GetOutput_FileFormat(void) { return Output_FileFormat; }

inline string CConfig::GetConv_FileName(void) { return Conv_FileName; }

inline string CConfig::GetConv_FileName_FSI(void) { return Conv_FileName_FSI; }

inline string CConfig::GetBreakdown_FileName(void) { return Breakdown_FileName; }

inline string CConfig::GetSolution_FlowFileName(void) { return Solution_FlowFileName; }

inline string CConfig::GetSolution_AdjFileName(void) { return Solution_AdjFileName; }

inline string CConfig::GetFlow_FileName(void) { return Flow_FileName; }

inline string CConfig::GetStructure_FileName(void) { return Structure_FileName; }

inline string CConfig::GetSurfStructure_FileName(void) { return SurfStructure_FileName; }

inline string CConfig::GetSurfWave_FileName(void) { return SurfWave_FileName; }

inline string CConfig::GetSurfHeat_FileName(void) { return SurfHeat_FileName; }

inline string CConfig::GetWave_FileName(void) { return Wave_FileName; }

inline string CConfig::GetHeat_FileName(void) { return Heat_FileName; }

inline string CConfig::GetAdjWave_FileName(void) { return AdjWave_FileName; }

inline string CConfig::GetRestart_FlowFileName(void) { return Restart_FlowFileName; }

inline string CConfig::GetRestart_WaveFileName(void) { return Restart_WaveFileName; }

inline string CConfig::GetRestart_HeatFileName(void) { return Restart_HeatFileName; }

inline string CConfig::GetRestart_AdjFileName(void) { return Restart_AdjFileName; }

inline string CConfig::GetAdj_FileName(void) { return Adj_FileName; }

inline string CConfig::GetObjFunc_Grad_FileName(void) { return ObjFunc_Grad_FileName; }

inline string CConfig::GetObjFunc_Value_FileName(void) { return ObjFunc_Value_FileName; }

inline string CConfig::GetSurfFlowCoeff_FileName(void) { return SurfFlowCoeff_FileName; }

inline string CConfig::GetSurfAdjCoeff_FileName(void) { return SurfAdjCoeff_FileName; }

inline unsigned short CConfig::GetResidual_Func_Flow(void) { return Residual_Func_Flow; }

inline unsigned short CConfig::GetCauchy_Func_Flow(void) { return Cauchy_Func_Flow; }

inline unsigned short CConfig::GetCauchy_Func_AdjFlow(void) { return Cauchy_Func_AdjFlow; }

inline unsigned short CConfig::GetCauchy_Elems(void) { return Cauchy_Elems; }

inline unsigned long CConfig::GetStartConv_Iter(void) { return StartConv_Iter; }

inline su2double CConfig::GetCauchy_Eps(void) { return Cauchy_Eps; }

inline su2double CConfig::GetDelta_UnstTimeND(void) { return Delta_UnstTimeND; }

inline su2double CConfig::GetTotal_UnstTimeND(void) { return Total_UnstTimeND; }

inline su2double CConfig::GetDelta_UnstTime(void) { return Delta_UnstTime; }

inline su2double CConfig::GetCurrent_UnstTime(void) { return Current_UnstTime; }

inline void CConfig::SetDelta_UnstTimeND(su2double val_delta_unsttimend) { Delta_UnstTimeND = val_delta_unsttimend; }

inline su2double CConfig::GetTotal_UnstTime(void) { return Total_UnstTime; }

inline bool CConfig::GetEngine_Intake(void) { return Engine_Intake; }

inline su2double CConfig::GetDV_Value(unsigned short val_dv) { return DV_Value[val_dv]; }

inline void CConfig::SetDV_Value(unsigned short val_dv, su2double val) { DV_Value[val_dv] = val; }

inline su2double CConfig::GetOrderMagResidual(void) { return OrderMagResidual; }

inline su2double CConfig::GetMinLogResidual(void) { return MinLogResidual; }

inline su2double CConfig::GetDamp_Engine_Inflow(void) { return Damp_Engine_Inflow; }

inline su2double CConfig::GetDamp_Engine_Bleed(void) { return Damp_Engine_Bleed; }

inline su2double CConfig::GetDamp_Engine_Exhaust(void) { return Damp_Engine_Exhaust; }

inline su2double CConfig::GetDamp_Res_Restric(void) { return Damp_Res_Restric; }

inline su2double CConfig::GetDamp_Correc_Prolong(void) { return Damp_Correc_Prolong; }

inline su2double CConfig::GetPosition_Plane(void) { return Position_Plane; }

inline su2double CConfig::GetWeightCd(void) { return WeightCd; }

inline su2double CConfig::GetFixAzimuthalLine(void) { return FixAzimuthalLine; }

inline su2double CConfig::GetCFLRedCoeff_Turb(void) { return CFLRedCoeff_Turb; }

inline bool CConfig::GetGrid_Movement(void) { return Grid_Movement; }

inline bool CConfig::GetRotating_Frame(void) { return Rotating_Frame; }

inline bool CConfig::GetAxisymmetric(void) { return Axisymmetric; }

inline bool CConfig::GetDebugMode(void) { return DebugMode; }

inline bool CConfig::GetAdaptBoundary(void) { return AdaptBoundary; }

inline bool CConfig::GetPoissonSolver(void) { return PoissonSolver; }

inline bool CConfig::Low_Mach_Preconditioning(void) { return Low_Mach_Precon; }

inline bool CConfig::GetGravityForce(void) { return GravityForce; }

inline bool CConfig::GetSmoothNumGrid(void) { return SmoothNumGrid; }

inline void CConfig::SetSmoothNumGrid(bool val_smoothnumgrid) { SmoothNumGrid = val_smoothnumgrid; }

inline unsigned short CConfig::GetKind_Turb_Model(void) { return Kind_Turb_Model; }

inline unsigned short CConfig::GetKind_Trans_Model(void) { return Kind_Trans_Model; }

inline bool CConfig::GetFrozen_Visc(void) { return Frozen_Visc; }

inline bool CConfig::GetSens_Remove_Sharp(void) { return Sens_Remove_Sharp; }

inline bool CConfig::GetViscous_Limiter_Flow(void) { return Viscous_Limiter_Flow; }

inline bool CConfig::GetViscous_Limiter_Turb(void) { return Viscous_Limiter_Turb; }

inline bool CConfig::GetWrite_Conv_FSI(void) { return Write_Conv_FSI; }

inline bool CConfig::GetHold_GridFixed(void) { return Hold_GridFixed; }

inline unsigned short CConfig::GetnPeriodicIndex(void) { return nPeriodic_Index; }

inline su2double* CConfig::GetPeriodicCenter(unsigned short val_index) { return Periodic_Center[val_index]; }

inline void CConfig::SetPeriodicCenter(unsigned short val_index, su2double* center) { Periodic_Center[val_index] = center; }

inline su2double* CConfig::GetPeriodicRotation(unsigned short val_index) { return Periodic_Rotation[val_index]; }

inline void CConfig::SetPeriodicRotation(unsigned short val_index, su2double* rotation) { Periodic_Rotation[val_index] = rotation; }

inline su2double* CConfig::GetPeriodicTranslate(unsigned short val_index) { return Periodic_Translate[val_index]; }

inline void CConfig::SetPeriodicTranslate(unsigned short val_index, su2double* translate) { Periodic_Translate[val_index] = translate; }

inline su2double CConfig::GetCyclic_Pitch(void) { return Cyclic_Pitch; }

inline su2double CConfig::GetCollective_Pitch(void) { return Collective_Pitch; }

inline string CConfig::GetMotion_FileName(void) { return Motion_Filename; }

inline bool CConfig::GetLow_MemoryOutput(void) { return Low_MemoryOutput; }

inline bool CConfig::GetWrt_Vol_Sol(void) { return Wrt_Vol_Sol; }

inline bool CConfig::GetWrt_Srf_Sol(void) { return Wrt_Srf_Sol; }

inline bool CConfig::GetWrt_Csv_Sol(void) { return Wrt_Csv_Sol; }

inline bool CConfig::GetWrt_Residuals(void) { return Wrt_Residuals; }

inline bool CConfig::GetWrt_Limiters(void) { return Wrt_Limiters; }

inline bool CConfig::GetWrt_SharpEdges(void) { return Wrt_SharpEdges; }

inline bool CConfig::GetWrt_Halo(void) { return Wrt_Halo; }

inline bool CConfig::GetPlot_Section_Forces(void) { return Plot_Section_Forces; }

inline bool CConfig::GetWrt_1D_Output(void) { return Wrt_1D_Output; }

inline vector<vector<su2double> > CConfig::GetAeroelastic_np1(unsigned short iMarker) {return Aeroelastic_np1[iMarker]; }

inline vector<vector<su2double> > CConfig::GetAeroelastic_n(unsigned short iMarker) {return Aeroelastic_n[iMarker]; }

inline vector<vector<su2double> > CConfig::GetAeroelastic_n1(unsigned short iMarker) {return Aeroelastic_n1[iMarker]; }

inline void CConfig::SetAeroelastic_np1(unsigned short iMarker, vector<vector<su2double> > solution) {Aeroelastic_np1[iMarker] = solution;}

inline su2double CConfig::GetAeroelastic_plunge(unsigned short val_marker) {return Aeroelastic_plunge[val_marker]; }

inline su2double CConfig::GetAeroelastic_pitch(unsigned short val_marker) {return Aeroelastic_pitch[val_marker]; }

inline void CConfig::SetAeroelastic_plunge(unsigned short val_marker, su2double val) {Aeroelastic_plunge[val_marker] = val; }

inline void CConfig::SetAeroelastic_pitch(unsigned short val_marker, su2double val) {Aeroelastic_pitch[val_marker] = val; }

inline void CConfig::SetAeroelastic_n1(void) {
        Aeroelastic_n1 = Aeroelastic_n;
}

inline void CConfig::SetAeroelastic_n(void) {
        Aeroelastic_n = Aeroelastic_np1;
}

inline su2double CConfig::GetAeroelastic_Flutter_Speed_Index(void) {return FlutterSpeedIndex; }

inline su2double CConfig::GetAeroelastic_Frequency_Plunge(void) {return PlungeNaturalFrequency; }

inline su2double CConfig::GetAeroelastic_Frequency_Pitch(void) {return PitchNaturalFrequency; }

inline su2double CConfig::GetAeroelastic_Airfoil_Mass_Ratio(void) {return AirfoilMassRatio; }

inline su2double CConfig::GetAeroelastic_CG_Location(void) {return CG_Location; }

inline su2double CConfig::GetAeroelastic_Radius_Gyration_Squared(void) {return RadiusGyrationSquared; }

inline unsigned short CConfig::GetAeroelasticIter(void) {return AeroelasticIter; }

inline bool CConfig::GetWind_Gust(void) { return Wind_Gust; }

inline bool CConfig::GetAeroelastic_Simulation(void) { return Aeroelastic_Simulation; }

inline unsigned short CConfig::GetGust_Type(void) {return Gust_Type; }

inline unsigned short CConfig::GetGust_Dir(void) {return Gust_Dir; }

inline su2double CConfig::GetGust_WaveLength(void) {return Gust_WaveLength; }

inline su2double CConfig::GetGust_Periods(void) {return Gust_Periods; }

inline su2double CConfig::GetGust_Ampl(void) {return Gust_Ampl; }

inline su2double CConfig::GetGust_Begin_Time(void) {return Gust_Begin_Time; }

inline su2double CConfig::GetGust_Begin_Loc(void) {return Gust_Begin_Loc; }

inline unsigned short CConfig::GetnFFD_Iter(void) {return nFFD_Iter; }

inline su2double CConfig::GetFFD_Tol(void) {return FFD_Tol; }

inline long CConfig::GetVisualize_CV(void) {return Visualize_CV; }

inline bool CConfig::GetFixed_CL_Mode(void) { return Fixed_CL_Mode; }

inline su2double CConfig::GetTarget_CL(void) {return Target_CL; }

inline su2double CConfig::GetDamp_Fixed_CL(void) {return Damp_Fixed_CL; }

inline unsigned long CConfig::GetIter_Fixed_CL(void) {return Iter_Fixed_CL; }

inline bool CConfig::GetUpdate_AoA(void) { return Update_AoA; }

inline void CConfig::SetUpdate_AoA(bool val_update) { Update_AoA = val_update; }

inline void CConfig::SetNonphysical_Points(unsigned long val_nonphys_points) { Nonphys_Points = val_nonphys_points; }

inline unsigned long CConfig::GetNonphysical_Points(void) { return Nonphys_Points; }

inline void CConfig::SetNonphysical_Reconstr(unsigned long val_nonphys_reconstr) { Nonphys_Reconstr = val_nonphys_reconstr; }

inline unsigned long CConfig::GetNonphysical_Reconstr(void) { return Nonphys_Reconstr; }

inline unsigned short CConfig::GetConsole_Output_Verb(void) { return Console_Output_Verb; }

inline unsigned short CConfig::GetnIterFSI(void) { return nIterFSI; }

inline su2double CConfig::GetAitkenStatRelax(void) { return AitkenStatRelax; }

inline su2double CConfig::GetAitkenDynMaxInit(void) { return AitkenDynMaxInit; }

inline bool CConfig::GetDeadLoad(void) { return DeadLoad; }

inline unsigned short CConfig::GetDynamic_Analysis(void) { return Dynamic_Analysis; }

inline su2double CConfig::GetDelta_DynTime(void) { return Delta_DynTime; }

inline su2double CConfig::GetTotal_DynTime(void) { return Total_DynTime; }

inline su2double CConfig::GetCurrent_DynTime(void) { return Current_DynTime; }

inline bool CConfig::GetWrt_Dynamic(void) { return Wrt_Dynamic; }

inline su2double CConfig::GetNewmark_alpha(void) { return Newmark_alpha; }

inline su2double CConfig::GetNewmark_delta(void) { return Newmark_delta; }

inline bool CConfig::GetGradual_Load(void) { return Gradual_Load; }

inline bool CConfig::GetRamp_Load(void) { return Ramp_Load; }

inline su2double CConfig::GetRamp_Time(void) { return Ramp_Time; }

inline su2double CConfig::GetStatic_Time(void) { return Static_Time; }

inline unsigned short CConfig::GetPredictorOrder(void) { return Pred_Order; }

inline bool CConfig::GetFSI_Simulation(void) { return FSI_Problem; }

inline unsigned short CConfig::GetRelaxation_Method_FSI(void) { return Kind_BGS_RelaxMethod; }

inline su2double CConfig::GetOrderMagResidualFSI(void) { return OrderMagResidualFSI; }

inline su2double CConfig::GetMinLogResidualFSI(void) { return MinLogResidualFSI; }

inline unsigned short CConfig::GetDirectDiff(){ return DirectDiff;}

inline bool CConfig::GetDiscrete_Adjoint() {return DiscreteAdjoint;}
