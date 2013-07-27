/*!
 * \file config_structure.inl
 * \brief In-Line subroutines of the <i>config_structure.hpp</i> file.
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

inline void CConfig::SetFanFace_Mach(unsigned short val_imarker, double val_fanface_mach) { FanFace_Mach[val_imarker] = val_fanface_mach; }

inline void CConfig::SetFanFace_Pressure(unsigned short val_imarker, double val_fanface_pressure) { FanFace_Pressure[val_imarker] = val_fanface_pressure; }

inline unsigned short CConfig::GetKind_SU2(void) { return Kind_SU2; }

inline bool CConfig::GetAdjoint(void) { return Adjoint; }

inline bool CConfig::GetViscous(void) { return Viscous; }

inline unsigned short CConfig::GetKind_Adjoint(void) { return Kind_Adjoint; }

inline unsigned long CConfig::GetnExtIter(void) { return nExtIter; }

inline unsigned short CConfig::GetnTimeInstances(void) { return nTimeInstances; }

inline double CConfig::GetTimeSpectral_Period(void) { return TimeSpectral_Period; }

inline void CConfig::SetExtIter(unsigned long val_iter) { ExtIter = val_iter; }

inline void CConfig::SetIntIter(unsigned long val_iter) { IntIter = val_iter; }

inline unsigned long CConfig::GetExtIter(void) { return ExtIter; }

inline unsigned long CConfig::GetIntIter(void) { return IntIter; }

inline unsigned long CConfig::GetUnst_nIntIter(void) { return Unst_nIntIter; }

inline unsigned short CConfig::GetMaxChildren(void) { return MaxChildren; }

inline string CConfig::GetPlaneTag(unsigned short index) { return PlaneTag[index]; }

inline double CConfig::GetEA_IntLimit(unsigned short index) { return EA_IntLimit[index]; }

inline double CConfig::GetAdjointLimit(void) { return AdjointLimit; }

inline double *CConfig::GetHold_GridFixed_Coord(void) { return Hold_GridFixed_Coord; }

inline double CConfig::GetMaxDimension(void) { return MaxDimension; }

inline double CConfig::GetRatioDensity(void) { return RatioDensity; }

inline double CConfig::GetFreeSurface_Thickness(void) { return FreeSurface_Thickness; }

inline double CConfig::GetFreeSurface_Damping_Coeff(void) { return FreeSurface_Damping_Coeff; }

inline double CConfig::GetFreeSurface_Damping_Length(void) { return FreeSurface_Damping_Length; }

inline unsigned long CConfig::GetFreeSurface_Reevaluation(void) { return FreeSurface_Reevaluation; }

inline double CConfig::GetFreeSurface_Outlet(void) { return FreeSurface_Outlet; }

inline double CConfig::GetRatioViscosity(void) { return RatioViscosity; }

inline unsigned short CConfig::GetAnalytical_Surface(void) { return Analytical_Surface; }

inline double CConfig::GetDualVol_Power(void) { return DualVol_Power; }

inline bool CConfig::GetVisualize_Partition(void) { return Visualize_Partition; }

inline bool CConfig::GetVisualize_Deformation(void) { return Visualize_Deformation; }

inline double CConfig::GetRefAreaCoeff(void) { return RefAreaCoeff; }

inline double CConfig::GetWaveSpeed(void) { return Wave_Speed; }

inline double CConfig::GetThermalDiffusivity(void) { return Thermal_Diffusivity; }

inline double CConfig::GetElasticyMod(void) { return ElasticyMod; }

inline double CConfig::GetPoissonRatio(void) { return PoissonRatio; }

inline double CConfig::GetMaterialDensity(void) { return MaterialDensity; }

inline double CConfig::GetRefLengthMoment(void) { return RefLengthMoment; }

inline double CConfig::GetRefElemLength(void) { return RefElemLength; }

inline double CConfig::GetDomainVolume(void) { return DomainVolume; }

inline void CConfig::SetRefAreaCoeff(double val_area) { RefAreaCoeff = val_area; }

inline void CConfig::SetDomainVolume(double val_volume) { DomainVolume = val_volume; }

inline void CConfig::SetnExtIter(unsigned long val_niter) { nExtIter = val_niter; }

inline double CConfig::GetMach_FreeStreamND(void) { return Mach; }

inline double CConfig::GetGamma(void) { return Gamma; }

inline double CConfig::GetSpecies_Gamma(unsigned short val_Species) { return Species_Gamma[val_Species]; }

inline int CConfig::GetCharge_Number(unsigned short val_Species) { return Charge_Number[val_Species]; }

inline int ***CConfig::GetReaction_Map(void) { return Reactions; } 

inline double ***CConfig::GetCollisionIntegral00(void) { return Omega00; }

inline double ***CConfig::GetCollisionIntegral11(void) { return Omega11; }

inline double CConfig::GetBulk_Modulus(void) { return Bulk_Modulus; }

inline double CConfig::GetGammaMonatomic(void) { return GammaMonatomic; }

inline double CConfig::GetGammaDiatomic(void) { return GammaDiatomic; }

inline double CConfig::GetArtComp_Factor(void) { return ArtComp_Factor; }

inline double CConfig::GetGas_Constant(void) { return Gas_Constant; }

inline double CConfig::GetGas_ConstantND(void) { return Gas_ConstantND; }

inline double CConfig::GetMixtureMolar_Mass(void) { return Mixture_Molar_mass; }

inline double CConfig::GetSpecies_Gas_Constant(unsigned short val_Species) { return Species_Gas_Constant[val_Species]; }

inline double CConfig::GetBlottnerCoeff(unsigned short val_Species, unsigned short val_Coeff) { return Blottner[val_Species][val_Coeff]; }

inline double CConfig::GetSpecies_Temperature(unsigned short val_Species) { return Species_Temperature_FreeStream[val_Species]; }

inline double CConfig::GetWallTemperature(void) { return Wall_Temperature; }

inline double CConfig::GetFreeSurface_Zero(void) { return FreeSurface_Zero; }

inline double CConfig::GetFreeSurface_Depth(void) { return FreeSurface_Depth; }

inline double CConfig::GetGas_Constant_Ref(void) { return Gas_Constant_Ref; }

inline double CConfig::GetTemperature_FreeStream(void) { return Temperature_FreeStream; }

inline double CConfig::GetPrandtl_Lam(void) { return Prandtl_Lam; }

inline double CConfig::GetPrandtl_Turb(void) { return Prandtl_Turb; }

inline double CConfig::GetLength_Ref(void) { return Length_Ref; }

inline double CConfig::GetPressure_Ref(void) { return Pressure_Ref; }

inline double CConfig::GetTemperature_Ref(void) { return Temperature_Ref; }

inline double CConfig::GetDensity_Ref(void) { return Density_Ref; }

inline double CConfig::GetVelocity_Ref(void) { return Velocity_Ref; }

inline double CConfig::GetTime_Ref(void) { return Time_Ref; }

inline double CConfig::GetViscosity_Ref(void) { return Viscosity_Ref; }

inline double CConfig::GetOmega_Ref(void) { return Omega_Ref; }

inline double CConfig::GetForce_Ref(void) { return Force_Ref; }

inline double CConfig::GetPressure_FreeStreamND(void) { return Pressure_FreeStreamND; }

inline double CConfig::GetPressure_FreeStream(void) { return Pressure_FreeStream; }

inline double CConfig::GetTemperature_FreeStreamND(void) { return Temperature_FreeStreamND; }

inline double CConfig::GetDensity_FreeStreamND(void) { return Density_FreeStreamND; }

inline double* CConfig::GetVelocity_FreeStreamND(void) { return Velocity_FreeStreamND; }

inline double* CConfig::GetVelocity_FreeStream(void) { return Velocity_FreeStream; }

inline double* CConfig::GetOmega_FreeStreamND(void) { return Omega_FreeStreamND; }

inline double CConfig::GetEnergy_FreeStreamND(void) { return Energy_FreeStreamND; }

inline double CConfig::GetViscosity_FreeStreamND(void) { return Viscosity_FreeStreamND; }

inline double CConfig::GetNuFactor_FreeStream(void) { return NuFactor_FreeStream; }

inline double CConfig::GetIntermittency_FreeStream(void) { return Intermittency_FreeStream; }

inline double CConfig::GetTurbulenceIntensity_FreeStream(void) { return TurbulenceIntensity_FreeStream; }

inline double CConfig::GetTurb2LamViscRatio_FreeStream(void) {return Turb2LamViscRatio_FreeStream;}

inline double CConfig::GetLength_Reynolds(void) { return Length_Reynolds; }

inline double CConfig::GetConversion_Factor(void) { return Conversion_Factor; }

inline unsigned short CConfig::GetnStartUpIter(void) { return nStartUpIter; }

inline double *CConfig::GetRefOriginMoment(void) { return RefOriginMoment; }

inline double CConfig::GetChargeCoeff(void) { return ChargeCoeff; }

inline double CConfig::GetLimiterCoeff(void) { return LimiterCoeff; }

inline double CConfig::GetReynolds(void) { return Reynolds; }

inline double CConfig::GetFroude(void) { return Froude; }

inline double CConfig::GetAoA(void) { return AoA; }

inline unsigned short CConfig::GetnDomain(void) { return nDomain; }

inline void CConfig::SetnDomain(unsigned short val_ndomain) { nDomain = val_ndomain; }

inline double CConfig::GetAoS(void) { return AoS; }

inline unsigned short CConfig::GetMGLevels(void) { return nMultiLevel; }

inline void CConfig::SetMGLevels(unsigned short val_nMultiLevel) { nMultiLevel = val_nMultiLevel; }

inline unsigned short CConfig::GetFinestMesh(void) { return FinestMesh; }

inline void CConfig::SetFinestMesh(unsigned short val_finestmesh) { FinestMesh = val_finestmesh; }

inline void CConfig::SubtractFinestMesh(void) { FinestMesh = FinestMesh-1; }

inline unsigned short CConfig::GetDesign_Variable(unsigned short val_dv) { return Design_Variable[val_dv]; }

inline unsigned short CConfig::GetConvCriteria(void) { return ConvCriteria; }

inline unsigned short CConfig::GetMGCycle(void) { return MGCycle; }

inline unsigned short CConfig::GetGeometryMode(void) { return GeometryMode; }

inline double CConfig::GetCFL(unsigned short val_mesh) {	return CFL[val_mesh]; }

inline double CConfig::GetCFL(unsigned short val_mesh, unsigned short val_Species) { return CFL_MS[val_Species][val_mesh]; }

inline double CConfig::GetUnst_CFL(void) {	return Unst_CFL; }

inline double CConfig::GetParamDV(unsigned short val_dv, unsigned short val_param) {	return ParamDV[val_dv][val_param]; }

inline unsigned short CConfig::GetnDV(void) {	return nDV; }

inline unsigned short CConfig::GetnRKStep(void) { return nRKStep; }

inline double CConfig::Get_Alpha_RKStep(unsigned short val_step) { return RK_Alpha_Step[val_step]; }

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

inline unsigned short CConfig::GetKind_Solver(void) { return Kind_Solver; }

inline unsigned short CConfig::GetKind_GasModel(void) { return Kind_GasModel; }

inline unsigned short CConfig::GetKind_GridMovement(unsigned short val_iZone) { return Kind_GridMovement[val_iZone]; }

inline void CConfig::SetKind_GridMovement(unsigned short val_iZone, unsigned short motion_Type) { Kind_GridMovement[val_iZone] = motion_Type; }

inline double CConfig::GetMach_Motion(void) { return Mach_Motion; }

inline double CConfig::GetMotion_Origin_X(unsigned short val_iZone) { return Motion_Origin_X[val_iZone]; }

inline double CConfig::GetMotion_Origin_Y(unsigned short val_iZone) { return Motion_Origin_Y[val_iZone]; }

inline double CConfig::GetMotion_Origin_Z(unsigned short val_iZone) { return Motion_Origin_Z[val_iZone]; }

inline void CConfig::SetMotion_Origin_X(unsigned short val_iZone, double val_origin) { Motion_Origin_X[val_iZone] = val_origin; }

inline void CConfig::SetMotion_Origin_Y(unsigned short val_iZone, double val_origin) { Motion_Origin_Y[val_iZone] = val_origin; }

inline void CConfig::SetMotion_Origin_Z(unsigned short val_iZone, double val_origin) { Motion_Origin_Z[val_iZone] = val_origin; }

inline double CConfig::GetTranslation_Rate_X(unsigned short val_iZone) { return  Translation_Rate_X[val_iZone]; }

inline double CConfig::GetTranslation_Rate_Y(unsigned short val_iZone) { return  Translation_Rate_Y[val_iZone]; }

inline double CConfig::GetTranslation_Rate_Z(unsigned short val_iZone) { return  Translation_Rate_Z[val_iZone]; }

inline double CConfig::GetRotation_Rate_X(unsigned short val_iZone) { return  Rotation_Rate_X[val_iZone]; }

inline double CConfig::GetRotation_Rate_Y(unsigned short val_iZone) { return  Rotation_Rate_Y[val_iZone]; }

inline double CConfig::GetRotation_Rate_Z(unsigned short val_iZone) { return  Rotation_Rate_Z[val_iZone]; }

inline double CConfig::GetPitching_Omega_X(unsigned short val_iZone) { return  Pitching_Omega_X[val_iZone]; }

inline double CConfig::GetPitching_Omega_Y(unsigned short val_iZone) { return  Pitching_Omega_Y[val_iZone]; }

inline double CConfig::GetPitching_Omega_Z(unsigned short val_iZone) { return  Pitching_Omega_Z[val_iZone]; }

inline double CConfig::GetPitching_Ampl_X(unsigned short val_iZone) { return  Pitching_Ampl_X[val_iZone]; }

inline double CConfig::GetPitching_Ampl_Y(unsigned short val_iZone) { return  Pitching_Ampl_Y[val_iZone]; }

inline double CConfig::GetPitching_Ampl_Z(unsigned short val_iZone) { return  Pitching_Ampl_Z[val_iZone]; }

inline double CConfig::GetPitching_Phase_X(unsigned short val_iZone) { return  Pitching_Phase_X[val_iZone]; }

inline double CConfig::GetPitching_Phase_Y(unsigned short val_iZone) { return  Pitching_Phase_Y[val_iZone]; }

inline double CConfig::GetPitching_Phase_Z(unsigned short val_iZone) { return  Pitching_Phase_Z[val_iZone]; }

inline double CConfig::GetPlunging_Omega_X(unsigned short val_iZone) { return  Plunging_Omega_X[val_iZone]; }

inline double CConfig::GetPlunging_Omega_Y(unsigned short val_iZone) { return  Plunging_Omega_Y[val_iZone]; }

inline double CConfig::GetPlunging_Omega_Z(unsigned short val_iZone) { return  Plunging_Omega_Z[val_iZone]; }

inline double CConfig::GetPlunging_Ampl_X(unsigned short val_iZone) { return  Plunging_Ampl_X[val_iZone]; }

inline double CConfig::GetPlunging_Ampl_Y(unsigned short val_iZone) { return  Plunging_Ampl_Y[val_iZone]; }

inline double CConfig::GetPlunging_Ampl_Z(unsigned short val_iZone) { return  Plunging_Ampl_Z[val_iZone]; }

inline double CConfig::GetminTurkelBeta() { return  Min_Beta_RoeTurkel; }

inline double CConfig::GetmaxTurkelBeta() { return  Max_Beta_RoeTurkel; }

inline unsigned short CConfig::GetnSpecies(void) { return nSpecies; }

inline unsigned short CConfig::GetnReactions(void) { return nReactions; }

inline double CConfig::GetArrheniusCoeff(unsigned short iReaction) { return ArrheniusCoefficient[iReaction]; }

inline double CConfig::GetArrheniusEta(unsigned short iReaction) { return ArrheniusEta[iReaction]; }

inline double CConfig::GetArrheniusTheta(unsigned short iReaction) { return ArrheniusTheta[iReaction]; }

inline double CConfig::GetCharVibTemp(unsigned short iSpecies) {return CharVibTemp[iSpecies]; }

inline double CConfig::GetParticle_Mass(unsigned short iSpecies) { return Particle_Mass[iSpecies]; } 

inline double CConfig::GetStagnation_B() { return Stagnation_B; } 

inline double CConfig::GetElec_Conductivity() { return Electric_Cond; } 

inline double CConfig::GetDipoleDist() { return DipoleDist; } 

inline double CConfig::GetMolar_Mass(unsigned short iSpecies) { return Molar_Mass[iSpecies]; } 

inline unsigned short CConfig::GetnMonatomics(void) { return nMonatomics; }
	
inline unsigned short CConfig::GetnDiatomics(void) { return nDiatomics; }	

inline double CConfig::GetMolecular_Diameter(unsigned short iSpecies) { return Molecular_Diameter[iSpecies]; } 

inline int CConfig::GetParticle_ChargeNumber(unsigned short iSpecies) { return Charge_Number[iSpecies]; }

inline double CConfig::GetInitial_Gas_Composition(unsigned short iSpecies) { return Gas_Composition[iSpecies]; }

inline double CConfig::GetEnthalpy_Formation(unsigned short iSpecies) { return Enthalpy_Formation[iSpecies]; }

inline double CConfig::GetTemperature_Ref(unsigned short iSpecies) { return Species_Ref_Temperature[iSpecies]; } 

inline double CConfig::GetViscosity_Ref(unsigned short iSpecies) { return Species_Ref_Viscosity[iSpecies]; } 

inline unsigned short CConfig::GetKind_Gradient_Method(void) { return Kind_Gradient_Method; }

inline unsigned short CConfig::GetKind_GridDef_Method(void) { return Kind_GridDef_Method; }

inline unsigned short CConfig::GetKind_Linear_Solver(void) { return Kind_Linear_Solver; }

inline unsigned short CConfig::GetKind_Linear_Solver_Prec(void) { return Kind_Linear_Solver_Prec; }

inline void CConfig::SetKind_Linear_Solver_Prec(unsigned short val_kind_prec) { Kind_Linear_Solver_Prec = val_kind_prec; }

inline double CConfig::GetLinear_Solver_Error(void) { return Linear_Solver_Error; }

inline unsigned long CConfig::GetLinear_Solver_Iter(void) { return Linear_Solver_Iter; }

inline double CConfig::GetLinear_Solver_Relax(void) { return Linear_Solver_Relax; }

inline unsigned short CConfig::GetKind_AdjTurb_Linear_Solver(void) { return Kind_AdjTurb_Linear_Solver; }

inline unsigned short CConfig::GetKind_AdjTurb_Linear_Prec(void) { return Kind_AdjTurb_Linear_Prec; }

inline void CConfig::SetKind_AdjTurb_Linear_Prec(unsigned short val_kind_prec) { Kind_AdjTurb_Linear_Prec = val_kind_prec; }

inline double CConfig::GetAdjTurb_Linear_Error(void) { return AdjTurb_Linear_Error; }

inline unsigned short CConfig::GetAdjTurb_Linear_Iter(void) { return AdjTurb_Linear_Iter; }

inline double CConfig::GetAdjTurb_CFLRedCoeff(void) { return AdjTurb_CFLRedCoeff; }

inline unsigned long CConfig::GetFEA_Iter(void) { return FEA_Iter; }

inline unsigned short CConfig::GetKind_Adaptation(void) { return Kind_Adaptation; }

inline double CConfig::GetNew_Elem_Adapt(void) { return New_Elem_Adapt; }

inline unsigned short CConfig::GetKind_TimeIntScheme(void) { return Kind_TimeNumScheme; }

inline unsigned short CConfig::GetKind_ConvNumScheme(void) { return Kind_ConvNumScheme; }

inline unsigned short CConfig::GetKind_ViscNumScheme(void) { return Kind_ViscNumScheme; }

inline unsigned short CConfig::GetKind_SourNumScheme(void) { return Kind_SourNumScheme; }

inline unsigned short CConfig::GetKind_Centered(void) { return Kind_Centered; }

inline unsigned short CConfig::GetKind_Upwind(void) { return Kind_Upwind; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Flow(void) { return Kind_TimeIntScheme_Flow; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Plasma(void) { return Kind_TimeIntScheme_Plasma; }

inline unsigned short CConfig::GetKind_TimeIntScheme_AdjPlasma(void) { return Kind_TimeIntScheme_AdjPlasma; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Wave(void) { return Kind_TimeIntScheme_Wave; }

inline unsigned short CConfig::GetKind_TimeIntScheme_FEA(void) { return Kind_TimeIntScheme_FEA; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Flow(void) { return Kind_ConvNumScheme_Flow; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Plasma(void) { return Kind_ConvNumScheme_Plasma; }

inline unsigned short CConfig::GetKind_ConvNumScheme_AdjPlasma(void) { return Kind_ConvNumScheme_AdjPlasma; }

inline unsigned short CConfig::GetKind_ConvNumScheme_LevelSet(void) { return Kind_ConvNumScheme_LevelSet; }

inline unsigned short CConfig::GetKind_ConvNumScheme_AdjLevelSet(void) { return Kind_ConvNumScheme_AdjLevelSet; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Template(void) { return Kind_ConvNumScheme_Template; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Flow(void) { return Kind_ViscNumScheme_Flow; }

inline unsigned short CConfig::GetKind_SourNumScheme_LevelSet(void) { return Kind_SourNumScheme_LevelSet; }

inline unsigned short CConfig::GetKind_SourNumScheme_Wave(void) { return Kind_SourNumScheme_Wave; }

inline unsigned short CConfig::GetKind_SourNumScheme_FEA(void) { return Kind_SourNumScheme_FEA; }

inline unsigned short CConfig::GetKind_ViscNumScheme_AdjLevelSet(void) { return Kind_ViscNumScheme_AdjLevelSet; }

inline unsigned short CConfig::GetKind_SourNumScheme_AdjLevelSet(void) { return Kind_SourNumScheme_AdjLevelSet; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Plasma(void) { return Kind_ViscNumScheme_Plasma; }

inline unsigned short CConfig::GetKind_ViscNumScheme_AdjPlasma(void) { return Kind_ViscNumScheme_AdjPlasma; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Template(void) { return Kind_ViscNumScheme_Template; }

inline unsigned short CConfig::GetKind_SourNumScheme_Flow(void) { return Kind_SourNumScheme_Flow; }

inline unsigned short CConfig::GetKind_SourNumScheme_Plasma(void) { return Kind_SourNumScheme_Plasma; }

inline unsigned short CConfig::GetKind_SourNumScheme_AdjPlasma(void) { return Kind_SourNumScheme_AdjPlasma; }

inline unsigned short CConfig::GetKind_SourJac_Plasma(void) { return Kind_SourJac_Plasma; }

inline unsigned short CConfig::GetKind_SourNumScheme_Template(void) { return Kind_SourNumScheme_Template; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Elec(void) { return Kind_ViscNumScheme_Elec; }

inline unsigned short CConfig::GetKind_SourNumScheme_Elec(void) { return Kind_SourNumScheme_Elec; }

inline unsigned short CConfig::GetKind_Centered_Flow(void) { return Kind_Centered_Flow; }

inline unsigned short CConfig::GetKind_Centered_LevelSet(void) { return Kind_Centered_LevelSet; }

inline unsigned short CConfig::GetKind_Centered_AdjLevelSet(void) { return Kind_Centered_AdjLevelSet; }

inline unsigned short CConfig::GetKind_Centered_Plasma(void) { return Kind_Centered_Plasma; }

inline unsigned short CConfig::GetKind_Centered_AdjPlasma(void) { return Kind_Centered_AdjPlasma; }

inline unsigned short CConfig::GetKind_SlopeLimit(void) { return Kind_SlopeLimit; }

inline unsigned short CConfig::GetKind_SlopeLimit_Flow(void) { return Kind_SlopeLimit_Flow; }

inline unsigned short CConfig::GetKind_SlopeLimit_Plasma(void) { return Kind_SlopeLimit_Plasma; }

inline unsigned short CConfig::GetKind_SlopeLimit_AdjPlasma(void) { return Kind_SlopeLimit_AdjPlasma; }

inline unsigned short CConfig::GetKind_SlopeLimit_Turb(void) { return Kind_SlopeLimit_Turb; }

inline unsigned short CConfig::GetKind_SlopeLimit_LevelSet(void) { return Kind_SlopeLimit_LevelSet; }

inline unsigned short CConfig::GetKind_SlopeLimit_AdjLevelSet(void) { return Kind_SlopeLimit_AdjLevelSet; }

inline unsigned short CConfig::GetKind_SlopeLimit_AdjTurb(void) { return Kind_SlopeLimit_AdjTurb; }

inline unsigned short CConfig::GetKind_SlopeLimit_AdjFlow(void) { return Kind_SlopeLimit_AdjFlow; }

inline unsigned short CConfig::GetKind_Upwind_Flow(void) { return Kind_Upwind_Flow; }

inline unsigned short CConfig::GetKind_Upwind_LevelSet(void) { return Kind_Upwind_LevelSet; }

inline unsigned short CConfig::GetKind_Upwind_AdjLevelSet(void) { return Kind_Upwind_AdjLevelSet; }

inline unsigned short CConfig::GetKind_Upwind_Plasma(void) { return Kind_Upwind_Plasma; }

inline unsigned short CConfig::GetKind_Upwind_AdjPlasma(void) { return Kind_Upwind_AdjPlasma; }

inline double CConfig::GetKappa_1st_Flow(void) { return Kappa_1st_Flow; }

inline double CConfig::GetKappa_2nd_Flow(void) { return Kappa_2nd_Flow; }

inline double CConfig::GetKappa_4th_Flow(void) { return Kappa_4th_Flow; }

inline double CConfig::GetKappa_1st_Plasma(void) { return Kappa_1st_Plasma; }

inline double CConfig::GetKappa_2nd_Plasma(void) { return Kappa_2nd_Plasma; }

inline double CConfig::GetKappa_4th_Plasma(void) { return Kappa_4th_Plasma; }

inline unsigned short CConfig::GetKind_TimeIntScheme_AdjFlow(void) { return Kind_TimeIntScheme_AdjFlow; }

inline unsigned short CConfig::GetKind_ConvNumScheme_AdjFlow(void) { return Kind_ConvNumScheme_AdjFlow; }

inline unsigned short CConfig::GetKind_ViscNumScheme_AdjFlow(void) { return Kind_ViscNumScheme_AdjFlow; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Wave(void) { return Kind_ViscNumScheme_Wave; }

inline unsigned short CConfig::GetKind_ViscNumScheme_FEA(void) { return Kind_ViscNumScheme_FEA; }

inline unsigned short CConfig::GetKind_SourNumScheme_AdjFlow(void) { return Kind_SourNumScheme_AdjFlow; }

inline unsigned short CConfig::GetKind_Centered_AdjFlow(void) { return Kind_Centered_AdjFlow; }

inline unsigned short CConfig::GetKind_Upwind_AdjFlow(void) { return Kind_Upwind_AdjFlow; }

inline double CConfig::GetKappa_1st_AdjFlow(void) { return Kappa_1st_AdjFlow; }

inline double CConfig::GetKappa_2nd_AdjFlow(void) { return Kappa_2nd_AdjFlow; }

inline double CConfig::GetKappa_4th_AdjFlow(void) { return Kappa_4th_AdjFlow; }

inline double CConfig::GetKappa_1st_AdjPlasma(void) { return Kappa_1st_AdjPlasma; }

inline double CConfig::GetKappa_2nd_AdjPlasma(void) { return Kappa_2nd_AdjPlasma; }

inline double CConfig::GetKappa_4th_AdjPlasma(void) { return Kappa_4th_AdjPlasma; }

inline unsigned short CConfig::GetKind_TimeIntScheme_LinFlow(void) { return Kind_TimeIntScheme_LinFlow; }

inline unsigned short CConfig::GetKind_ConvNumScheme_LinFlow(void) { return Kind_ConvNumScheme_LinFlow; }

inline unsigned short CConfig::GetKind_ViscNumScheme_LinFlow(void) { return Kind_ViscNumScheme_LinFlow; }

inline unsigned short CConfig::GetKind_SourNumScheme_LinFlow(void) { return Kind_SourNumScheme_LinFlow; }

inline unsigned short CConfig::GetKind_Centered_LinFlow(void) { return Kind_Centered_LinFlow; }

inline unsigned short CConfig::GetKind_Upwind_LinFlow(void) { return Kind_Upwind_LinFlow; }

inline double CConfig::GetKappa_4th_LinFlow(void) { return Kappa_4th_LinFlow; }

inline double CConfig::GetKappa_1st_LinFlow(void) { return Kappa_1st_LinFlow; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Turb(void) { return Kind_TimeIntScheme_Turb; }

inline unsigned short CConfig::GetKind_TimeIntScheme_LevelSet(void) { return Kind_TimeIntScheme_LevelSet; }

inline unsigned short CConfig::GetKind_TimeIntScheme_AdjLevelSet(void) { return Kind_TimeIntScheme_AdjLevelSet; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Turb(void) { return Kind_ConvNumScheme_Turb; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Turb(void) { return Kind_ViscNumScheme_Turb; }

inline unsigned short CConfig::GetKind_SourNumScheme_Turb(void) { return Kind_SourNumScheme_Turb; }

inline unsigned short CConfig::GetKind_Centered_Turb(void) { return Kind_Centered_Turb; }

inline unsigned short CConfig::GetKind_Upwind_Turb(void) {	return Kind_Upwind_Turb; }

inline unsigned short CConfig::GetKind_TimeIntScheme_AdjTurb(void) { return Kind_TimeIntScheme_AdjTurb; }

inline unsigned short CConfig::GetKind_ConvNumScheme_AdjTurb(void) { return Kind_ConvNumScheme_AdjTurb; }

inline unsigned short CConfig::GetKind_ViscNumScheme_AdjTurb(void) { return Kind_ViscNumScheme_AdjTurb; }

inline unsigned short CConfig::GetKind_SourNumScheme_AdjTurb(void) { return Kind_SourNumScheme_AdjTurb; }

inline unsigned short CConfig::GetKind_Centered_AdjTurb(void) { return Kind_Centered_AdjTurb; }

inline unsigned short CConfig::GetKind_Upwind_AdjTurb(void) { return Kind_Upwind_AdjTurb; }

inline unsigned short CConfig::GetKind_Inlet(void) { return Kind_Inlet; }

inline void CConfig::SetKind_TimeIntScheme(unsigned short val_kind_timeintscheme) { Kind_TimeNumScheme = val_kind_timeintscheme; }

inline void CConfig::SetKind_ViscNumScheme(unsigned short val_kind_viscnumscheme) { Kind_ViscNumScheme = val_kind_viscnumscheme; }

inline void CConfig::SetKind_SourNumScheme(unsigned short val_kind_sournumscheme) { Kind_SourNumScheme = val_kind_sournumscheme; }

inline unsigned short CConfig::GetKind_ObjFunc(void) {return Kind_ObjFunc; }

inline unsigned short CConfig::GetKind_GeoObjFunc(void) {return Kind_GeoObjFunc; }

inline unsigned short CConfig::GetKind_SensSmooth(void) {return Kind_SensSmooth; }

inline unsigned short CConfig::GetContinuous_Eqns(void) {return Continuous_Eqns; }

inline unsigned short CConfig::GetDiscrete_Eqns(void) {return Discrete_Eqns; }

inline unsigned short CConfig::GetUnsteady_Simulation(void) { return Unsteady_Simulation; }

inline bool CConfig::GetRestart(void) {	return Restart; }

inline bool CConfig::GetRestart_Flow(void) { return Restart_Flow; }

inline bool CConfig::GetFullMG(void) { return FullMG; }

inline bool CConfig::GetEquivArea(void) { return EquivArea; }

inline bool CConfig::GetFlowRate(void) { return FlowRate; }

inline double CConfig::GetReduced_Frequency(void) { return Reduced_Frequency; }

inline double CConfig::GetPitching_Amplitude(void) { return Pitching_Amplitude; }

inline void CConfig::SetnMarker_All(unsigned short val_nmarker) { nMarker_All = val_nmarker; }

inline string CConfig::GetMarker_All_Tag(unsigned short val_marker) { return Marker_All_Tag[val_marker]; }

inline string CConfig::GetMarker_NacelleInflow(unsigned short val_marker) { return Marker_NacelleInflow[val_marker]; }

inline string CConfig::GetMarker_NacelleExhaust(unsigned short val_marker) { return Marker_NacelleExhaust[val_marker]; }

inline unsigned short CConfig::GetTag_Marker_All(string val_tag) {
	for (unsigned short iMarker = 0; iMarker < nMarker_All; iMarker++) {
		if (val_tag == Marker_All_Tag[iMarker])
		return iMarker; 
	}
	cout <<"Ups, I don't find the boundary: "<< val_tag << endl; return 0;
}

inline unsigned short CConfig::GetMarker_All_Boundary(unsigned short val_marker) { return Marker_All_Boundary[val_marker]; }

inline void CConfig::SetMarker_All_Boundary(unsigned short val_marker, unsigned short val_boundary) { Marker_All_Boundary[val_marker] = val_boundary; }

inline void CConfig::SetMarker_All_Tag(unsigned short val_marker, string val_index) { Marker_All_Tag[val_marker] = val_index; }

inline void CConfig::SetMarker_All_Monitoring(unsigned short val_marker, unsigned short val_monitoring) { Marker_All_Monitoring[val_marker] = val_monitoring; }

inline void CConfig::SetMarker_All_Designing(unsigned short val_marker, unsigned short val_designing) { Marker_All_Designing[val_marker] = val_designing; }

inline void CConfig::SetMarker_All_Plotting(unsigned short val_marker, unsigned short val_plotting) { Marker_All_Plotting[val_marker] = val_plotting; }

inline void CConfig::SetMarker_All_Moving(unsigned short val_marker, unsigned short val_moving) { Marker_All_Moving[val_marker] = val_moving; }

inline void CConfig::SetMarker_All_PerBound(unsigned short val_marker, short val_perbound) { Marker_All_PerBound[val_marker] = val_perbound; }

inline short CConfig::GetMarker_All_PerBound(unsigned short val_marker) { return Marker_All_PerBound[val_marker]; }

inline void CConfig::SetMarker_All_Sliding(unsigned short val_marker, unsigned short val_slidebound) { Marker_All_Sliding[val_marker] = val_slidebound; }

inline unsigned short CConfig::GetMarker_All_Sliding(unsigned short val_marker) { return Marker_All_Sliding[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Monitoring(unsigned short val_marker) { return Marker_All_Monitoring[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Designing(unsigned short val_marker) { return Marker_All_Designing[val_marker]; }

inline short CConfig::GetMarker_All_SendRecv(unsigned short val_marker) { return Marker_All_SendRecv[val_marker]; }

inline void CConfig::SetMarker_All_SendRecv(unsigned short val_marker, short val_index) { Marker_All_SendRecv[val_marker] = val_index; }

inline unsigned short CConfig::GetMarker_All_Plotting(unsigned short val_marker) { return Marker_All_Plotting[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Moving(unsigned short val_marker) { return Marker_All_Moving[val_marker]; }

inline unsigned short CConfig::GetnMarker_All(void) { return nMarker_All; }

inline unsigned short CConfig::GetnMarker_NacelleInflow(void) {	return nMarker_NacelleInflow; }

inline unsigned short CConfig::GetnMarker_NacelleExhaust(void) { return nMarker_NacelleExhaust; }

inline unsigned short CConfig::GetnMarker_InterfaceBound(void) { return nMarker_InterfaceBound; }

inline unsigned short CConfig::GetnMarker_NearFieldBound(void) { return nMarker_NearFieldBound; }

inline string CConfig::GetMesh_FileName(void) { return Mesh_FileName; }

inline string CConfig::GetMesh_Out_FileName(void) { return Mesh_Out_FileName; }

inline unsigned short CConfig::GetMesh_FileFormat(void) { return Mesh_FileFormat; }

inline unsigned short CConfig::GetOutput_FileFormat(void) { return Output_FileFormat; }

inline string CConfig::GetConv_FileName(void) { return Conv_FileName; }

inline string CConfig::GetSolution_FlowFileName(void) { return Solution_FlowFileName; }

inline string CConfig::GetFarfield_FileName(void) { return Farfield_FileName; }

inline string CConfig::GetSolution_LinFileName(void) { return Solution_LinFileName; }

inline string CConfig::GetSolution_AdjFileName(void) { return Solution_AdjFileName; }

inline string CConfig::GetFlow_FileName(void) { return Flow_FileName; }

inline string CConfig::GetStructure_FileName(void) { return Structure_FileName; }

inline string CConfig::GetWave_FileName(void) { return Wave_FileName; }

inline string CConfig::GetAdjWave_FileName(void) { return AdjWave_FileName; }

inline string CConfig::GetRestart_FlowFileName(void) { return Restart_FlowFileName; }

inline string CConfig::GetRestart_WaveFileName(void) { return Restart_WaveFileName; }

inline string CConfig::GetRestart_HeatFileName(void) { return Restart_HeatFileName; }

inline string CConfig::GetRestart_LinFileName(void) { return Restart_LinFileName; }

inline string CConfig::GetRestart_AdjFileName(void) { return Restart_AdjFileName; }

inline string CConfig::GetAdj_FileName(void) { return Adj_FileName; }

inline string CConfig::GetLin_FileName(void) { return Lin_FileName; }

inline string CConfig::GetObjFunc_Grad_FileName(void) { return ObjFunc_Grad_FileName; }

inline string CConfig::GetObjFunc_Eval_FileName(void) { return ObjFunc_Eval_FileName; }

inline string CConfig::GetSurfFlowCoeff_FileName(void) { return SurfFlowCoeff_FileName; }

inline string CConfig::GetSurfAdjCoeff_FileName(void) { return SurfAdjCoeff_FileName; }

inline string CConfig::GetSurfLinCoeff_FileName(void) { return SurfLinCoeff_FileName; }

inline unsigned short CConfig::GetCauchy_Func_Flow(void) { return Cauchy_Func_Flow; }

inline unsigned short CConfig::GetCauchy_Func_AdjFlow(void) { return Cauchy_Func_AdjFlow; }

inline unsigned short CConfig::GetCauchy_Func_LinFlow(void) { return Cauchy_Func_LinFlow; }

inline unsigned short CConfig::GetCauchy_Elems(void) { return Cauchy_Elems; }

inline unsigned long CConfig::GetStartConv_Iter(void) { return StartConv_Iter; }

inline double CConfig::GetCauchy_Eps(void) { return Cauchy_Eps; }

inline double CConfig::GetCauchy_Eps_OneShot(void) { return Cauchy_Eps_OneShot; }

inline double CConfig::GetCauchy_Eps_FullMG(void) { return Cauchy_Eps_FullMG; }

inline double CConfig::GetDelta_UnstTimeND(void) { return Delta_UnstTimeND; }

inline double CConfig::GetDelta_UnstTime(void) { return Delta_UnstTime; }

inline double CConfig::GetCurrent_UnstTime(void) { return Current_UnstTime; }

inline void CConfig::SetDelta_UnstTimeND(double val_delta_unsttimend) { Delta_UnstTimeND = val_delta_unsttimend; }

inline double CConfig::GetTotal_UnstTime(void) { return Total_UnstTime; }

inline bool CConfig::GetDivide_Element(void) { return Divide_Element; }

inline bool CConfig::GetEngine_Intake(void) { return Engine_Intake; }

inline double CConfig::GetDV_Value_New(unsigned short val_dv) { return DV_Value_New[val_dv]; }

inline double CConfig::GetOrderMagResidual(void) { return OrderMagResidual; }

inline double CConfig::GetMinLogResidual(void) { return MinLogResidual; }

inline double CConfig::GetDamp_Nacelle_Inflow(void) { return Damp_Nacelle_Inflow; }

inline double CConfig::GetDamp_Res_Restric(void) { return Damp_Res_Restric; }

inline double CConfig::GetDamp_Correc_Prolong(void) { return Damp_Correc_Prolong; }

inline double CConfig::GetPosition_Plane(void) { return Position_Plane; }

inline double CConfig::GetWeightCd(void) { return WeightCd; }

inline double CConfig::GetCteViscDrag(void) { return CteViscDrag; }

inline double CConfig::GetLevelSet_CFLRedCoeff(void) { return LevelSet_CFLRedCoeff; }

inline double CConfig::GetTurb_CFLRedCoeff(void) { return Turb_CFLRedCoeff; }

inline double *CConfig::GetRotAxisOrigin(void) { return RotAxisOrigin; }

inline double *CConfig::GetOmega(void) { return Omega; }

inline double CConfig::GetOmegaMag(void) { return Omega_Mag; }

inline double CConfig::GetDV_Value_Old(unsigned short val_dv) { return DV_Value_Old[val_dv]; }

inline bool CConfig::GetGrid_Movement(void) { return Grid_Movement; }

inline bool CConfig::GetRotating_Frame(void) { return Rotating_Frame; }

inline bool CConfig::GetAxisymmetric(void) { return Axisymmetric; }

inline bool CConfig::GetAdaptBoundary(void) { return AdaptBoundary; }

inline bool CConfig::GetIncompressible(void) { return Incompressible; }

inline bool CConfig::GetFreeSurface(void) { return FreeSurface; }

inline bool CConfig::GetAdiabaticWall(void) { return AdiabaticWall; }

inline bool CConfig::GetIsothermalWall(void) { return IsothermalWall; }

inline bool CConfig::GetCatalyticWall(void) { return CatalyticWall; }

inline bool CConfig::GetElectricSolver(void) { return ElectricSolver; }

inline bool CConfig::Low_Mach_Preconditioning(void) { return Low_Mach_Precon; }

inline bool CConfig::GetMacCormackRelaxation(void) { return MacCormackRelaxation; }

inline bool CConfig::GetUnsteady_Farfield(void) { return Unsteady_Farfield; }

inline bool CConfig::MultipleTimeSteps(void) { return PlasmaMultiTimeSteps; }

inline bool CConfig::GetInletConditionsDefined(void) { return Inlet_Outlet_Defined; }

inline double CConfig::GetInlet_Species_Temperature(unsigned short iSpecies) { return Species_Temperature_Inlet[iSpecies]; }

inline double CConfig::GetOutlet_Species_Temperature(unsigned short iSpecies) { return Species_Temperature_Outlet[iSpecies]; }

inline double CConfig::GetInlet_Species_Pressure(unsigned short iSpecies) { return Species_Pressure_Inlet[iSpecies]; }

inline double CConfig::GetOutlet_Species_Pressure(unsigned short iSpecies) { return Species_Pressure_Outlet[iSpecies]; }

inline double CConfig::GetInlet_Species_Velocity(unsigned short iSpecies) { return Species_Velocity_Inlet[iSpecies]; }

inline double CConfig::GetOutlet_Species_Velocity(unsigned short iSpecies) { return Species_Velocity_Outlet[iSpecies]; }

inline bool CConfig::GetGravityForce(void) { return GravityForce; }

inline bool CConfig::GetMagnetic_Force(void) { return MagneticForce; }

inline bool CConfig::GetJouleHeating(void) { return JouleHeating; }

inline bool CConfig::GetSmoothNumGrid(void) { return SmoothNumGrid; }

inline void CConfig::SetSmoothNumGrid(bool val_smoothnumgrid) { SmoothNumGrid = val_smoothnumgrid; }

inline unsigned short CConfig::GetKind_Turb_Model(void) { return Kind_Turb_Model; }

inline unsigned short CConfig::GetKind_Trans_Model(void) { return Kind_Trans_Model; }

inline bool CConfig::GetFrozen_Visc(void) { return Frozen_Visc; }

inline bool CConfig::GetShow_Adj_Sens(void) { return Show_Adj_Sens; }

inline bool CConfig::GetHold_GridFixed(void) { return Hold_GridFixed; }

inline bool CConfig::GetCGNS_To_SU2(void) {return CGNS_To_SU2; }

inline bool CConfig::GetWrite_Converted_Mesh(void) { return Write_Converted_Mesh; }

inline unsigned short CConfig::GetnPeriodicIndex(void) { return nPeriodic_Index; }

inline double* CConfig::GetPeriodicCenter(unsigned short val_index) { return Periodic_Center[val_index]; }

inline void CConfig::SetPeriodicCenter(unsigned short val_index, double* center) { Periodic_Center[val_index] = center; }

inline double* CConfig::GetPeriodicRotation(unsigned short val_index) { return Periodic_Rotation[val_index]; }

inline void CConfig::SetPeriodicRotation(unsigned short val_index, double* rotation) { Periodic_Rotation[val_index] = rotation; }

inline double* CConfig::GetPeriodicTranslate(unsigned short val_index) { return Periodic_Translate[val_index]; }

inline void CConfig::SetPeriodicTranslate(unsigned short val_index, double* translate) { Periodic_Translate[val_index] = translate; }

inline double CConfig::GetCyclic_Pitch(void) { return Cyclic_Pitch; }

inline double CConfig::GetCollective_Pitch(void) { return Collective_Pitch; }

inline string CConfig::GetMotion_FileName(void) { return Motion_Filename; }

inline bool CConfig::GetWrt_Vol_Sol(void) { return Wrt_Vol_Sol; }

inline bool CConfig::GetWrt_Srf_Sol(void) { return Wrt_Srf_Sol; }

inline bool CConfig::GetWrt_Csv_Sol(void) { return Wrt_Csv_Sol; }

inline bool CConfig::GetWrt_Restart(void) { return Wrt_Restart; }

inline bool CConfig::GetWrt_Sol_CGNS(void) { return Wrt_Sol_CGNS; }

inline bool CConfig::GetWrt_Sol_Tec_ASCII(void) { return Wrt_Sol_Tec_ASCII; }

inline bool CConfig::GetWrt_Sol_Tec_Binary(void) { return Wrt_Sol_Tec_Binary; }

inline bool CConfig::GetWrt_Residuals(void) { return Wrt_Residuals; }

inline bool CConfig::GetWrt_Halo(void) { return Wrt_Halo; }

inline bool CConfig::GetRelative_Motion(void) { return Relative_Motion; }

inline double* CConfig::GetAeroelastic_np1(void) {return Aeroelastic_np1; }

inline double* CConfig::GetAeroelastic_n(void) {return Aeroelastic_n; }
    
inline double* CConfig::GetAeroelastic_n1(void) {return Aeroelastic_n1; }

inline void CConfig::SetAeroelastic_np1(unsigned short val_index, double val) {Aeroelastic_np1[val_index] = val;}

inline double CConfig::GetAeroelastic_plunge(void) {return Aeroelastic_plunge;}

inline double CConfig::GetAeroelastic_pitch(void) {return Aeroelastic_pitch;}

inline void CConfig::SetAeroelastic_plunge(double val) {Aeroelastic_plunge = val;}

inline void CConfig::SetAeroelastic_pitch(double val) {Aeroelastic_pitch = val;}

inline void CConfig::SetAeroelastic_n1(void) {
    for (unsigned short i=0; i<4; i++)
        Aeroelastic_n1[i] = Aeroelastic_n[i];
}

inline void CConfig::SetAeroelastic_n(void) {
    for (unsigned short i=0; i<4; i++)
        Aeroelastic_n[i] = Aeroelastic_np1[i];
}
    
inline double CConfig::GetAeroelastic_Frequency_Plunge(void) {return FreqPlungeAeroelastic;}

inline double CConfig::GetAeroelastic_Frequency_Pitch(void) {return FreqPitchAeroelastic;}

inline unsigned short CConfig::GetType_Aeroelastic(void) {return Aeroelastic_Grid_Movement;}

inline unsigned short CConfig::GetAeroelastic_GridVelocity(void) {return Aeroelastic_Grid_Velocity;}

inline	double  CConfig::GetDensity_FreeStreamND_Time(unsigned long val_Ext_Iter) {return Density_FreeStreamND_Time[val_Ext_Iter];}

inline	double  CConfig::GetPressure_FreeStreamND_Time(unsigned long val_Ext_Iter) {return Pressure_FreeStreamND_Time[val_Ext_Iter];}

inline	double  CConfig::GetEnergy_FreeStreamND_Time(unsigned long val_Ext_Iter) {return Energy_FreeStreamND_Time[val_Ext_Iter];}

inline	double  CConfig::GetMach_FreeStreamND_Time(unsigned long val_Ext_Iter) {return Mach_Inf_Time[val_Ext_Iter];}

inline	double*  CConfig::GetVelocity_FreeStreamND_Time(unsigned long val_Ext_Iter) {return Velocity_FreeStreamND_Time[val_Ext_Iter];}

