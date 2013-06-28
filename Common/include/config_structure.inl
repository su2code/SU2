/*!
 * \file config_structure.inl
 * \brief In-Line subroutines of the <i>config_structure.hpp</i> file.
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

inline unsigned short CConfig::GetKind_SU2(void) { return Kind_SU2; }

inline unsigned long CConfig::GetnExtIter(void) { return nExtIter; }

inline unsigned long CConfig::GetnIntIter(void) { return nIntIter; }

inline unsigned short CConfig::GetMaxChildren(void) { return MaxChildren; }

inline string CConfig::GetPlaneTag(unsigned short index) { return PlaneTag[index]; }

inline double CConfig::GetEA_IntLimit(unsigned short index) { return EA_IntLimit[index]; }

inline double CConfig::GetMaxDimension(void) { return MaxDimension; }

inline double CConfig::GetRatioDensity(void) { return RatioDensity; }

inline double CConfig::GetInterfaseThickness(void) { return InterfaseThickness; }

inline double CConfig::GetRatioViscosity(void) { return RatioViscosity; }

inline unsigned short CConfig::GetAnalytical_Surface(void) { return Analytical_Surface; }

inline double CConfig::GetDualVol_Power(void) { return DualVol_Power; }

inline bool CConfig::GetReStartMGCycle(void) { return ReStartMGCycle; }

inline bool CConfig::GetShockTube(void) { return ShockTube; }

inline bool CConfig::GetStretchingFactor_Flow(void) { return StretchingFactor_Flow; }

inline bool CConfig::GetStretchingFactor_Adj(void) { return StretchingFactor_Adj; }

inline bool CConfig::GetStretchingFactor_Lin(void) { return StretchingFactor_Lin; }

inline bool CConfig::GetVisualize_Partition(void) { return Visualize_Partition; }

inline bool CConfig::GetVisualize_Deformation(void) { return Visualize_Deformation; }

inline double CConfig::GetRefAreaCoeff(void) { return RefAreaCoeff; }

inline double CConfig::GetRefLengthMoment(void) { return RefLengthMoment; }

inline double CConfig::GetDomainVolume(void) { return DomainVolume; }

inline void CConfig::SetRefAreaCoeff(double val_area) { RefAreaCoeff = val_area; }

inline void CConfig::SetDomainVolume(double val_volume) { DomainVolume = val_volume; }

inline void CConfig::SetnExtIter(unsigned long val_niter) { nExtIter = val_niter; }

inline double CConfig::GetMach_FreeStreamND(void) { return Mach; }

inline double CConfig::GetGamma(void) { return Gamma; }

inline double CConfig::GetBulk_Modulus(void) { return Bulk_Modulus; }

inline double CConfig::GetGammaMonatomic(void) { return GammaMonatomic; }

inline double CConfig::GetGammaDiatomic(void) { return GammaDiatomic; }

inline double CConfig::GetArtComp_Factor(void) { return ArtComp_Factor; }

inline double CConfig::GetArtComp_Min(void) { return ArtComp_Min; }

inline double CConfig::GetGas_Constant(void) { return Gas_Constant; }

inline double CConfig::GetLevelSet_Zero(void) { return LevelSet_Zero; }

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

inline double CConfig::GetTemperature_FreeStreamND(void) { return Temperature_FreeStreamND; }

inline double CConfig::GetDensity_FreeStreamND(void) { return Density_FreeStreamND; }

inline double* CConfig::GetVelocity_FreeStreamND(void) { return Velocity_FreeStreamND; }

inline double* CConfig::GetOmega_FreeStreamND(void) { return Omega_FreeStreamND; }

inline double CConfig::GetEnergy_FreeStreamND(void) { return Energy_FreeStreamND; }

inline double CConfig::GetViscosity_FreeStreamND(void) { return Viscosity_FreeStreamND; }

inline double CConfig::GetLength_Reynolds(void) { return Length_Reynolds; }

inline double CConfig::GetConversion_Factor(void) { return Conversion_Factor; }

inline double CConfig::GetPrimGrad_Threshold(void) { return PrimGrad_Threshold; }

inline unsigned short CConfig::GetnStartUpIter(void) { return nStartUpIter; }

inline double *CConfig::GetRefOriginMoment(void) { return RefOriginMoment; }

inline double CConfig::GetChargeCoeff(void) { return ChargeCoeff; }

inline double CConfig::GetSmoothCoeff(void) { return SmoothCoeff; }

inline double CConfig::GetLimiterCoeff(void) { return LimiterCoeff; }

inline unsigned short CConfig::GetnSmooth(void) { return nSmooth; }

inline double CConfig::GetReynolds(void) { return Reynolds; }

inline double CConfig::GetAoA(void) { return AoA; }

inline unsigned short CConfig::GetnDomain(void) { return nDomain; }

inline void CConfig::SetnDomain(unsigned short val_ndomain) { nDomain = val_ndomain; }

inline double CConfig::GetAoS(void) { return AoS; }

inline unsigned short CConfig::GetMGLevels(void) { return nMultiLevel; }

inline unsigned short CConfig::GetFinestMesh(void) { return FinestMesh; }

inline void CConfig::SetFinestMesh(unsigned short val_finestmesh) { FinestMesh = val_finestmesh; };

inline void CConfig::SubtractFinestMesh(void) { FinestMesh = FinestMesh-1; }

inline unsigned short CConfig::GetDesign_Variable(unsigned short val_dv) { return Design_Variable[val_dv]; }

inline unsigned short CConfig::GetConvCriteria(void) { return ConvCriteria; }

inline unsigned short CConfig::GetMGCycle(void) { return MGCycle; }

inline double CConfig::GetCFL(unsigned short val_mesh) {	return CFL[val_mesh]; }

inline double CConfig::GetParamDV(unsigned short val_dv, unsigned short val_param) {	return ParamDV[val_dv][val_param]; }

inline unsigned short CConfig::GetnDV(void) {	return nDV; }

inline unsigned short CConfig::GetnRKStep(void) { return nRKStep; }

inline double CConfig::Get_Alpha_RKStep(unsigned short val_step) { return RK_Alpha_Step[val_step]; }

inline double CConfig::Get_Beta_RKStep(unsigned short val_step) { return RK_Beta_Step[val_step]; }

inline unsigned short CConfig::GetMG_PreSmooth(unsigned short val_mesh) {	return MG_PreSmooth[val_mesh]; }

inline unsigned short CConfig::GetMG_PostSmooth(unsigned short val_mesh) { return MG_PostSmooth[val_mesh]; }

inline unsigned short CConfig::GetMG_CorrecSmooth(unsigned short val_mesh) { return MG_CorrecSmooth[val_mesh]; }

inline unsigned long CConfig::GetWrt_Sol_Freq(void) { return Wrt_Sol_Freq; }

inline unsigned long CConfig::GetWrt_Con_Freq(void) { return Wrt_Con_Freq; }

inline unsigned short CConfig::GetKind_Solver(void) { return Kind_Solver; }

inline unsigned short CConfig::GetKind_GasModel(void) { return Kind_GasModel; }

inline unsigned short CConfig::GetnSpecies(void) { return nSpecies; }

inline unsigned short CConfig::GetnFluids(void) {  return nFluids; } 

inline double CConfig::GetParticle_Mass(unsigned short iSpecies) { return Particle_Mass[iSpecies]; } 

inline unsigned short CConfig::GetMPI_Kind_Solver(void) { return Kind_Solver_MPI; }

inline void CConfig::SetKind_Solver(unsigned short val_kind_solver) { Kind_Solver = val_kind_solver; }

inline void CConfig::SetMPI_Kind_Solver(unsigned short val_mpi_kind_solver) { Kind_Solver_MPI = val_mpi_kind_solver; }

inline unsigned short CConfig::GetKind_Gradient_Method(void) { return Kind_Gradient_Method; };

inline unsigned short CConfig::GetKind_GridDef_Method(void) { return Kind_GridDef_Method; };

inline unsigned short CConfig::GetKind_GridDef_Solver(void) { return Kind_GridDef_Solver; };

inline double CConfig::GetGridDef_Error(void) { return GridDef_Error; };

inline unsigned short CConfig::GetKind_Adaptation(void) { return Kind_Adaptation; }

inline double CConfig::GetNew_Elem_Adapt(void) { return New_Elem_Adapt; }

inline unsigned short CConfig::GetKind_TimeIntScheme(void) { return Kind_TimeNumScheme; }

inline unsigned short CConfig::GetKind_ConvNumScheme(void) { return Kind_ConvNumScheme; }

inline unsigned short CConfig::GetKind_ViscNumScheme(void) { return Kind_ViscNumScheme; }

inline unsigned short CConfig::GetKind_SourNumScheme(void) { return Kind_SourNumScheme; }

inline unsigned short CConfig::GetKind_Centred(void) { return Kind_Centred; }

inline unsigned short CConfig::GetKind_Upwind(void) { return Kind_Upwind; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Flow(void) { return Kind_TimeIntScheme_Flow; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Plasma(void) { return Kind_TimeIntScheme_Plasma; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Flow(void) { return Kind_ConvNumScheme_Flow; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Plasma(void) { return Kind_ConvNumScheme_Plasma; }

inline unsigned short CConfig::GetKind_ConvNumScheme_LevelSet(void) { return Kind_ConvNumScheme_LevelSet; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Combustion(void) { return Kind_ConvNumScheme_Combustion; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Template(void) { return Kind_ConvNumScheme_Template; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Flow(void) { return Kind_ViscNumScheme_Flow; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Plasma(void) { return Kind_ViscNumScheme_Plasma; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Template(void) { return Kind_ViscNumScheme_Template; }

inline unsigned short CConfig::GetKind_SourNumScheme_Flow(void) { return Kind_SourNumScheme_Flow; }

inline unsigned short CConfig::GetKind_SourNumScheme_Plasma(void) { return Kind_SourNumScheme_Plasma; }

inline unsigned short CConfig::GetKind_SourNumScheme_Combustion(void) { return Kind_SourNumScheme_Combustion; }

inline unsigned short CConfig::GetKind_SourNumScheme_Template(void) { return Kind_SourNumScheme_Template; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Elec(void) { return Kind_ViscNumScheme_Elec; }

inline unsigned short CConfig::GetKind_SourNumScheme_Elec(void) { return Kind_SourNumScheme_Elec; }

inline unsigned short CConfig::GetKind_Centred_Flow(void) { return Kind_Centred_Flow; }

inline unsigned short CConfig::GetKind_Centred_LevelSet(void) { return Kind_Centred_LevelSet; }

inline unsigned short CConfig::GetKind_Centred_Combustion(void) { return Kind_Centred_Combustion; }

inline unsigned short CConfig::GetKind_Centred_Plasma(void) { return Kind_Centred_Plasma; }

inline unsigned short CConfig::GetKind_SlopeLimit_Flow(void) { return Kind_SlopeLimit_Flow; }

inline unsigned short CConfig::GetKind_SlopeLimit_Turb(void) { return Kind_SlopeLimit_Turb; }

inline unsigned short CConfig::GetKind_SlopeLimit_LevelSet(void) { return Kind_SlopeLimit_LevelSet; }

inline unsigned short CConfig::GetKind_SlopeLimit_AdjTurb(void) { return Kind_SlopeLimit_AdjTurb; }

inline unsigned short CConfig::GetKind_Upwind_Flow(void) { return Kind_Upwind_Flow; }

inline unsigned short CConfig::GetKind_Upwind_LevelSet(void) { return Kind_Upwind_LevelSet; }

inline unsigned short CConfig::GetKind_Upwind_Combustion(void) { return Kind_Upwind_Combustion; }

inline unsigned short CConfig::GetKind_Upwind_Plasma(void) { return Kind_Upwind_Plasma; }

inline double CConfig::GetKappa_1st_Flow(void) { return Kappa_1st_Flow; }

inline double CConfig::GetKappa_2nd_Flow(void) { return Kappa_2nd_Flow; }

inline double CConfig::GetKappa_4th_Flow(void) { return Kappa_4th_Flow; }

inline unsigned short CConfig::GetKind_TimeIntScheme_AdjFlow(void) { return Kind_TimeIntScheme_Adj; }

inline unsigned short CConfig::GetKind_ConvNumScheme_AdjFlow(void) { return Kind_ConvNumScheme_Adj; }

inline unsigned short CConfig::GetKind_ViscNumScheme_AdjFlow(void) { return Kind_ViscNumScheme_Adj; }

inline unsigned short CConfig::GetKind_SourNumScheme_AdjFlow(void) { return Kind_SourNumScheme_Adj; }

inline unsigned short CConfig::GetKind_Centred_AdjFlow(void) { return Kind_Centred_Adj; }

inline unsigned short CConfig::GetKind_Upwind_AdjFlow(void) { return Kind_Upwind_Adj; }

inline double CConfig::GetKappa_4th_Adj(void) { return Kappa_4th_Adj; }

inline double CConfig::GetKappa_1st_Adj(void) { return Kappa_1st_Adj; }

inline unsigned short CConfig::GetKind_TimeIntScheme_LinFlow(void) { return Kind_TimeIntScheme_Lin; }

inline unsigned short CConfig::GetKind_ConvNumScheme_LinFlow(void) { return Kind_ConvNumScheme_Lin; }

inline unsigned short CConfig::GetKind_ViscNumScheme_LinFlow(void) { return Kind_ViscNumScheme_Lin; }

inline unsigned short CConfig::GetKind_SourNumScheme_LinFlow(void) { return Kind_SourNumScheme_Lin; }

inline unsigned short CConfig::GetKind_Centred_LinFlow(void) { return Kind_Centred_Lin; }

inline unsigned short CConfig::GetKind_Upwind_LinFlow(void) { return Kind_Upwind_Lin; }

inline double CConfig::GetKappa_4th_Lin(void) { return Kappa_4th_Lin; }

inline double CConfig::GetKappa_1st_Lin(void) { return Kappa_1st_Lin; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Turb(void) { return Kind_TimeIntScheme_Turb; }

inline unsigned short CConfig::GetKind_TimeIntScheme_LevelSet(void) { return Kind_TimeIntScheme_LevelSet; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Combustion(void) { return Kind_TimeIntScheme_Combustion; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Turb(void) { return Kind_ConvNumScheme_Turb; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Turb(void) { return Kind_ViscNumScheme_Turb; }

inline unsigned short CConfig::GetKind_SourNumScheme_Turb(void) { return Kind_SourNumScheme_Turb; }

inline unsigned short CConfig::GetKind_Centred_Turb(void) { return Kind_Centred_Turb; }

inline unsigned short CConfig::GetKind_Upwind_Turb(void) {	return Kind_Upwind_Turb; }

inline unsigned short CConfig::GetKind_TimeIntScheme_AdjTurb(void) { return Kind_TimeIntScheme_AdjTurb; }

inline unsigned short CConfig::GetKind_ConvNumScheme_AdjTurb(void) { return Kind_ConvNumScheme_AdjTurb; }

inline unsigned short CConfig::GetKind_ViscNumScheme_AdjTurb(void) { return Kind_ViscNumScheme_AdjTurb; }

inline unsigned short CConfig::GetKind_SourNumScheme_AdjTurb(void) { return Kind_SourNumScheme_AdjTurb; }

inline unsigned short CConfig::GetKind_Centred_AdjTurb(void) { return Kind_Centred_AdjTurb; }

inline unsigned short CConfig::GetKind_Upwind_AdjTurb(void) { return Kind_Upwind_AdjTurb; }

inline void CConfig::SetKind_TimeIntScheme(unsigned short val_kind_timeintscheme) { Kind_TimeNumScheme = val_kind_timeintscheme; }

inline void CConfig::SetKind_ViscNumScheme(unsigned short val_kind_viscnumscheme) { Kind_ViscNumScheme = val_kind_viscnumscheme; }

inline void CConfig::SetKind_SourNumScheme(unsigned short val_kind_sournumscheme) { Kind_SourNumScheme = val_kind_sournumscheme; }

inline unsigned short CConfig::GetKind_ObjFunc(void) {return Kind_ObjFunc; }

inline unsigned short CConfig::GetUnsteady_Simulation(void) { return Unsteady_Simulation; }

inline bool CConfig::GetRestart(void) {	return Restart; }

inline bool CConfig::GetRestart_Flow(void) { return Restart_Flow; }

inline bool CConfig::GetFullMG(void) { return FullMG; }

inline bool CConfig::GetEquivArea(void) { return EquivArea; }

inline bool CConfig::GetBiGrid(void) { return BiGrid; }

inline double CConfig::GetReduced_Frequency(void) { return Reduced_Frequency; }

inline double CConfig::GetPitching_Amplitude(void) { return Pitching_Amplitude; }

inline void CConfig::SetnMarker_All(unsigned short val_nmarker) { nMarker_All = val_nmarker; }

inline string CConfig::GetMarker_All_Tag(unsigned short val_marker) { return Marker_All_Tag[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Boundary(unsigned short val_marker) { return Marker_All_Boundary[val_marker]; }

inline void CConfig::SetMarker_All_Boundary(unsigned short val_marker, unsigned short val_boundary) { Marker_All_Boundary[val_marker] = val_boundary; }

inline void CConfig::SetMarker_All_Tag(unsigned short val_marker, string val_index) { Marker_All_Tag[val_marker] = val_index; }

inline void CConfig::SetMarker_All_Monitoring(unsigned short val_marker, unsigned short val_monitoring) { Marker_All_Monitoring[val_marker] = val_monitoring; }

inline void CConfig::SetMarker_All_Plotting(unsigned short val_marker, unsigned short val_plotting) { Marker_All_Plotting[val_marker] = val_plotting; }

inline void CConfig::SetMarker_All_Moving(unsigned short val_marker, unsigned short val_moving) { Marker_All_Moving[val_marker] = val_moving; }

inline void CConfig::SetMarker_All_PerBound(unsigned short val_marker, short val_perbound) { Marker_All_PerBound[val_marker] = val_perbound; }

inline short CConfig::GetMarker_All_PerBound(unsigned short val_marker) { return Marker_All_PerBound[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Monitoring(unsigned short val_marker) { return Marker_All_Monitoring[val_marker]; }

inline short CConfig::GetMarker_All_SendRecv(unsigned short val_marker) { return Marker_All_SendRecv[val_marker]; }

inline void CConfig::SetMarker_All_SendRecv(unsigned short val_marker, short val_index) { Marker_All_SendRecv[val_marker] = val_index; }

inline unsigned short CConfig::GetMarker_All_Plotting(unsigned short val_marker) { return Marker_All_Plotting[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Moving(unsigned short val_marker) { return Marker_All_Moving[val_marker]; }

inline unsigned short CConfig::GetnMarker_All(void) {	return nMarker_All; }

inline string CConfig::GetMesh_FileName(void) { return Mesh_FileName; }

inline string CConfig::GetMesh_Out_FileName(void) { return Mesh_Out_FileName; }

inline unsigned short CConfig::GetMesh_FileFormat(void) { return Mesh_FileFormat; }

inline unsigned short CConfig::GetOutput_FileFormat(void) { return Output_FileFormat; }

inline string CConfig::GetConv_FileName(void) { return Conv_FileName; }

inline string CConfig::GetSolution_FlowFileName(void) { return Solution_FlowFileName; }

inline string CConfig::GetSolution_LinFileName(void) { return Solution_LinFileName; }

inline string CConfig::GetSolution_AdjFileName(void) { return Solution_AdjFileName; }

inline string CConfig::GetFlow_FileName(void) { return Flow_FileName; }

inline string CConfig::GetReStart_FlowFileName(void) { return ReStart_FlowFileName; }

inline string CConfig::GetReStart_LinFileName(void) { return ReStart_LinFileName; }

inline string CConfig::GetReStart_AdjFileName(void) { return ReStart_AdjFileName; }

inline string CConfig::GetAdj_FileName(void) { return Adj_FileName; }

inline string CConfig::GetLin_FileName(void) { return Lin_FileName; }

inline string CConfig::GetObjFunc_Grad_FileName(void) { return ObjFunc_Grad_FileName; }

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

inline double CConfig::GetUnst_TimeStep(void) { return Unst_TimeStep; }

inline double CConfig::GetUnst_Time(void) { return Unst_Time; }

inline bool CConfig::GetDivide_Element(void) { return Divide_Element; }

inline double CConfig::GetDV_Value_New(unsigned short val_dv) { return DV_Value_New[val_dv]; }

inline double CConfig::GetOrderMagResidual(void) { return OrderMagResidual; }

inline double CConfig::GetMinLogResidual(void) { return MinLogResidual; }

inline double CConfig::GetDamp_Res_Restric(void) { return Damp_Res_Restric; }

inline double CConfig::GetDamp_Correc_Prolong(void) { return Damp_Correc_Prolong; }

inline double CConfig::GetPosition_Plane(void) { return Position_Plane; }

inline double CConfig::GetWeightCd(void) { return WeightCd; }

inline double *CConfig::GetRotAxisOrigin(void) { return RotAxisOrigin; }

inline double *CConfig::GetOmega(void) { return Omega; }

inline double CConfig::GetRotVel_Ref(void) { return RotVel_Ref/Velocity_Ref; }

inline double CConfig::GetDV_Value_Old(unsigned short val_dv) { return DV_Value_Old[val_dv]; }

inline double CConfig::GetSens_Limit(void) { return Sens_Limit; }

inline bool CConfig::GetGrid_Movement(void) { return Grid_Movement; }

inline bool CConfig::GetFreeSurface_Movement(void) { return FreeSurface_Movement; }

inline bool CConfig::GetRotating_Frame(void) { return Rotating_Frame; }

inline bool CConfig::GetIncompressible(void) { return Incompressible; }

inline bool CConfig::GetGravityForce(void) { return GravityForce; }

inline bool CConfig::GetSmoothNumGrid(void) { return SmoothNumGrid; }

inline void CConfig::SetSmoothNumGrid(bool val_smoothnumgrid) { SmoothNumGrid = val_smoothnumgrid; }

inline unsigned short CConfig::GetKind_Turb_Model(void) { return Kind_Turb_Model; }

inline bool CConfig::GetFrozen_Visc(void) { return Frozen_Visc; }

inline bool CConfig::GetMove_FFDBox(void) { return Move_FFDBox; }

inline bool CConfig::GetCGNS_To_SU2(void) {return CGNS_To_SU2; }

inline string CConfig::GetNew_SU2_FileName(void) { return New_SU2_FileName; }

inline unsigned short CConfig::GetnPeriodicIndex(void) { return nPeriodic_Index; }

inline double* CConfig::GetPeriodicCenter(unsigned short val_index) { return Periodic_Center[val_index]; }

inline void CConfig::SetPeriodicCenter(unsigned short val_index, double* center) { Periodic_Center[val_index] = center; }

inline double* CConfig::GetPeriodicRotation(unsigned short val_index) { return Periodic_Rotation[val_index]; }

inline void CConfig::SetPeriodicRotation(unsigned short val_index, double* rotation) { Periodic_Rotation[val_index] = rotation; }

inline double* CConfig::GetPeriodicTranslate(unsigned short val_index) { return Periodic_Translate[val_index]; }

inline void CConfig::SetPeriodicTranslate(unsigned short val_index, double* translate) { Periodic_Translate[val_index] = translate; }

