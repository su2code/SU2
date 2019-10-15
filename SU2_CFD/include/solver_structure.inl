/*!
 * \file solver_structure.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

inline void CSolver::SetIterLinSolver(unsigned short val_iterlinsolver) { IterLinSolver = val_iterlinsolver; }

inline void CSolver::SetNondimensionalization(CConfig *config, unsigned short iMesh) { }

inline bool CSolver::GetAdjoint(void) { return adjoint; }

inline unsigned short CSolver::GetIterLinSolver(void) { return IterLinSolver; }

inline su2double CSolver::GetCSensitivity(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep, 
                                     unsigned short iMesh, unsigned short RunTime_EqSystem) { }
                                     
inline void CSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) { }

inline void CSolver::ResetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) { }

inline void CSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) { }

inline void CSolver::LoadRestart_FSI(CGeometry *geometry, CConfig *config, int val_iter) { }

inline void CSolver::PredictStruct_Displacement(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) { }

inline void CSolver::ComputeAitken_Coefficient(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution, unsigned long iOuterIter) { }

inline void CSolver::SetAitken_Relaxation(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) { }

inline void CSolver::Update_StructSolution(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) { }

inline void CSolver::Compute_OFRefGeom(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Compute_OFRefNode(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Compute_OFVolFrac(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Compute_OFCompliance(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::SetForceCoeff(su2double val_forcecoeff_history) { }

inline void CSolver::SetFSI_Residual(su2double val_FSI_residual) { }

inline void CSolver::SetRelaxCoeff(su2double val_relaxecoeff_history) { }

inline su2double CSolver::GetRelaxCoeff(void) { return 0.0; }

inline su2double CSolver::GetForceCoeff(void) { return 0.0; }

inline su2double CSolver::GetFSI_Residual(void) { return 0.0; }

inline void CSolver::Stiffness_Penalty(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container, CConfig *config) { }

inline void CSolver::SetCSensitivity(unsigned short val_marker, unsigned long val_vertex, su2double val_sensitivity) { }

inline void CSolver::Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline void CSolver::Smooth_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline void CSolver::Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline su2double CSolver::GetPhi_Inf(unsigned short val_dim) { return 0; }

inline su2double CSolver::GetPsiRho_Inf(void) { return 0; }

inline su2double* CSolver::GetPsiRhos_Inf(void) { return NULL; }

inline su2double CSolver::GetPsiE_Inf(void) { return 0; }

inline void CSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPrimitive_Limiter_MPI(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPreconditioner(CConfig *config, unsigned long iPoint) { }

inline void CSolver::SetDistance(CGeometry *geometry, CConfig *config) { };

inline su2double CSolver::GetCD_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCL_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CL(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CD(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CSF(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CEff(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFx(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFy(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFz(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMx(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMy(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMz(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CL_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CD_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CSF_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CEff_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFx_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFy_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFz_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMx_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMy_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMz_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CL_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CD_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CSF_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CEff_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFx_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFy_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFz_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMx_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMy_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMz_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_Buffet_Metric(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CL_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CD_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CSF_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CEff_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFx_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFy_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFz_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMx_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMy_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMz_Mnt(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetInflow_MassFlow(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetExhaust_MassFlow(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetInflow_Pressure(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetInflow_Mach(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCSF_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCEff_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_HF_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_MaxHF_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCL_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCSF_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCD_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetAllBound_CL_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CD_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CSF_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CEff_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CMx_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CMy_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CMz_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CoPx_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CoPy_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CoPz_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CFx_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CFy_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CFz_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CL_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CD_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CSF_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CEff_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CMx_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CMy_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CMz_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CoPx_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CoPy_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CoPz_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CFx_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CFy_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CFz_Mnt() { return 0; }

inline su2double CSolver::GetAllBound_CL_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CD_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CSF_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CEff_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CMx_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CMy_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CMz_Visc() { return 0; }

inline su2double CSolver::GetTotal_Buffet_Metric() { return 0; }

inline su2double CSolver::GetAllBound_CoPx_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CoPy_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CoPz_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CFx_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CFy_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CFz_Visc() { return 0; }

inline void CSolver::SetForceProj_Vector(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::SetIntBoundary_Jump(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline su2double CSolver::GetTotal_CL() { return 0; }

inline su2double CSolver::GetTotal_CD() { return 0; }

inline su2double CSolver::GetTotal_NetThrust() { return 0; }

inline su2double CSolver::GetTotal_Power() { return 0; }

inline su2double CSolver::GetTotal_SolidCD() { return 0; }

inline su2double CSolver::GetTotal_ReverseFlow() { return 0; }

inline su2double CSolver::GetTotal_MFR() { return 0; }

inline su2double CSolver::GetTotal_Prop_Eff() { return 0; }

inline su2double CSolver::GetTotal_ByPassProp_Eff() { return 0; }

inline su2double CSolver::GetTotal_Adiab_Eff() { return 0; }

inline su2double CSolver::GetTotal_Poly_Eff() { return 0; }

inline su2double CSolver::GetTotal_IDC_Mach() { return 0; }

inline su2double CSolver::GetTotal_DC60() { return 0; }

inline su2double CSolver::GetTotal_Custom_ObjFunc() { return 0; }

inline su2double CSolver::GetTotal_CMx() { return 0; }

inline su2double CSolver::GetTotal_CMy() { return 0; }

inline su2double CSolver::GetTotal_CMz() { return 0; }

inline su2double CSolver::GetTotal_CoPx() { return 0; }

inline su2double CSolver::GetTotal_CoPy() { return 0; }

inline su2double CSolver::GetTotal_CoPz() { return 0; }

inline su2double CSolver::GetTotal_CFx() { return 0; }

inline su2double CSolver::GetTotal_CFy() { return 0; }

inline su2double CSolver::GetTotal_CFz() { return 0; }

inline su2double CSolver::GetTotal_CSF() { return 0; }

inline su2double CSolver::GetTotal_CEff() { return 0; }

inline su2double CSolver::GetTotal_CT() { return 0; }

inline void CSolver::SetTotal_CT(su2double val_Total_CT) { }

inline su2double CSolver::GetTotal_CQ() { return 0; }

inline su2double CSolver::GetTotal_HeatFlux() { return 0; }

inline su2double CSolver::GetTotal_AvgTemperature() { return 0; }

inline su2double CSolver::GetTotal_MaxHeatFlux() { return 0; }

inline su2double CSolver::Get_PressureDrag() { return 0; }

inline su2double CSolver::Get_ViscDrag() { return 0; }

inline void CSolver::SetTotal_CQ(su2double val_Total_CQ) { }

inline void CSolver::SetTotal_HeatFlux(su2double val_Total_Heat) { }

inline void CSolver::SetTotal_MaxHeatFlux(su2double val_Total_Heat) { }

inline su2double CSolver::GetTotal_CMerit() { return 0; }

inline su2double CSolver::GetTotal_CEquivArea() { return 0; }

inline su2double CSolver::GetTotal_AeroCD() { return 0; }

inline su2double CSolver::GetTotal_IDR() { return 0; }

inline su2double CSolver::GetTotal_IDC() { return 0; }

inline su2double CSolver::GetTotal_CpDiff() { return 0; }

inline su2double CSolver::GetTotal_HeatFluxDiff() { return 0; }

inline su2double CSolver::GetTotal_CFEA() { return 0; }

inline su2double CSolver::GetTotal_CNearFieldOF() { return 0; }

inline su2double CSolver::GetTotal_OFRefGeom() { return 0; }

inline su2double CSolver::GetTotal_OFRefNode() { return 0; }

inline su2double CSolver::GetTotal_OFVolFrac() { return 0; }

inline su2double CSolver::GetTotal_OFCompliance() { return 0; }

inline bool CSolver::IsElementBased(void){ return false; }

inline void CSolver::AddTotal_ComboObj(su2double val_obj) {}

inline void CSolver::SetTotal_CEquivArea(su2double val_cequivarea) { }

inline void CSolver::SetTotal_AeroCD(su2double val_aerocd) { }

inline void CSolver::SetTotal_CpDiff(su2double val_pressure) { }

inline void CSolver::SetTotal_HeatFluxDiff(su2double val_heat) { }

inline void CSolver::SetTotal_CFEA(su2double val_cfea) { }

inline void CSolver::SetTotal_OFRefGeom(su2double val_ofrefgeom) { }

inline void CSolver::SetTotal_OFRefNode(su2double val_ofrefnode) { }

inline su2double CSolver::GetWAitken_Dyn(void) { return 0; }

inline su2double CSolver::GetWAitken_Dyn_tn1(void) { return 0; }

inline void CSolver::SetWAitken_Dyn(su2double waitk) {  }

inline void CSolver::SetWAitken_Dyn_tn1(su2double waitk_tn1) {  }

inline void CSolver::SetLoad_Increment(su2double val_loadIncrement) {  }

inline su2double CSolver::GetLoad_Increment() { return 0; }

inline void CSolver::SetTotal_CNearFieldOF(su2double val_cnearfieldpress) { }

inline su2double CSolver::GetTotal_CWave() { return 0; }

inline su2double CSolver::GetTotal_CHeat() { return 0; }

inline void CSolver::SetTotal_CL(su2double val_Total_CL) { }

inline void CSolver::SetTotal_CD(su2double val_Total_CD) { }

inline void CSolver::SetTotal_NetThrust(su2double val_Total_NetThrust) { }

inline void CSolver::SetTotal_Power(su2double val_Total_Power) { }

inline void CSolver::SetTotal_SolidCD(su2double val_Total_SolidCD) { }

inline void CSolver::SetTotal_ReverseFlow(su2double val_Total_ReverseFlow) { }

inline void CSolver::SetTotal_MFR(su2double val_Total_MFR) { }

inline void CSolver::SetTotal_Prop_Eff(su2double val_Total_Prop_Eff) { }

inline void CSolver::SetTotal_ByPassProp_Eff(su2double val_Total_ByPassProp_Eff) { }

inline void CSolver::SetTotal_Adiab_Eff(su2double val_Total_Adiab_Eff) { }

inline void CSolver::SetTotal_Poly_Eff(su2double val_Total_Poly_Eff) { }

inline void CSolver::SetTotal_IDC(su2double val_Total_IDC) { }

inline void CSolver::SetTotal_IDC_Mach(su2double val_Total_IDC_Mach) { }

inline void CSolver::SetTotal_IDR(su2double val_Total_IDR) { }

inline void CSolver::SetTotal_DC60(su2double val_Total_DC60) { }

inline void CSolver::SetTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { }

inline void CSolver::AddTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { }

inline su2double CSolver::GetCPressure(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetCPressureTarget(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetCPressureTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_pressure) { }

inline void CSolver::SetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_heat) { }

inline su2double *CSolver::GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value) { }

inline su2double *CSolver::GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value) { }

inline void CSolver::SetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value) { }

inline su2double CSolver::GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var) { return 0; }

inline su2double *CSolver::GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var) { return 0; }

inline unsigned long CSolver::GetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex, unsigned long val_index) { }

inline su2double CSolver::GetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex, su2double val_deltap) { }

inline su2double CSolver::GetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex, su2double val_deltat) { }

inline su2double CSolver::GetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim) { return 0; }

inline void CSolver::SetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex, su2double val_ttotal) { }

inline void CSolver::SetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex, su2double val_ptotal) { }

inline void CSolver::SetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_flowdir) { }

inline void CSolver::SetInlet_TurbVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_turb_var) { }

inline void CSolver::SetUniformInlet(CConfig* config, unsigned short iMarker) {};

inline void CSolver::SetInletAtVertex(su2double *val_inlet, unsigned short iMarker, unsigned long iVertex) { };

inline su2double CSolver::GetInletAtVertex(su2double *val_inlet, unsigned long val_inlet_point, unsigned short val_kind_marker, string val_marker, CGeometry *geometry, CConfig *config) { return 0; }

inline void CSolver::UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config) { }

inline su2double CSolver::GetCSkinFriction(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim) { return 0; }

inline su2double CSolver::GetHeatFlux(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetBuffetSensor(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetYPlus(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetStrainMag_Max(void) { return 0; }

inline su2double CSolver::GetOmega_Max(void) { return 0; }

inline void CSolver::SetStrainMag_Max(su2double val_strainmag_max) { }

inline void CSolver::SetOmega_Max(su2double val_omega_max) { }

inline void CSolver::Viscous_Residual(CGeometry *geometry,
                                      CSolver **solver_container,
                                      CNumerics *numerics, CConfig
                                      *config, unsigned short iMesh,
                                      unsigned short iRKstep) { }
                     
inline void CSolver::AddStiffMatrix(su2double ** StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1, unsigned long Point_2, unsigned long Point_3) { }
                     
inline void CSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, 
                          CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) { }
                     
inline void CSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, 
                          CNumerics *numerics, CConfig *config, unsigned short iMesh) { }

inline su2double CSolver::GetTotal_Sens_Geo() { return 0; }

inline su2double CSolver::GetTotal_Sens_Mach() { return 0; }

inline su2double CSolver::GetTotal_Sens_AoA() { return 0; }

inline su2double CSolver::GetTotal_Sens_Press() { return 0; }

inline su2double CSolver::GetTotal_Sens_Temp() { return 0; }

inline su2double CSolver::GetTotal_Sens_BPress() { return 0; }

inline su2double CSolver::GetTotal_Sens_Density() { return 0; }

inline su2double CSolver::GetTotal_Sens_ModVel() { return 0; }

inline su2double CSolver::GetDensity_Inf(void) { return 0; }

inline su2double CSolver::GetDensity_Inf(unsigned short val_var) { return 0; }

inline su2double CSolver::GetModVelocity_Inf(void) { return 0; }

inline su2double CSolver::GetDensity_Energy_Inf(void) { return 0; }

inline su2double CSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return 0; }

inline su2double CSolver::GetDensity_Velocity_Inf(unsigned short val_dim, unsigned short val_var) { return 0; }

inline su2double CSolver::GetVelocity_Inf(unsigned short val_dim) { return 0; }

inline su2double* CSolver::GetVelocity_Inf(void) { return 0; }

inline su2double CSolver::GetPressure_Inf(void) { return 0; }

inline su2double CSolver::GetViscosity_Inf(void) { return 0; }

inline su2double CSolver::GetNuTilde_Inf(void) { return 0; }

inline su2double CSolver::GetTke_Inf(void) { return 0; }

inline su2double CSolver::GetOmega_Inf(void) { return 0; }

inline su2double CSolver::GetTotal_Sens_E(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetTotal_Sens_Nu(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetTotal_Sens_Rho(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetTotal_Sens_Rho_DL(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetTotal_Sens_EField(unsigned short iEField) { return 0.0; }

inline su2double CSolver::GetTotal_Sens_DVFEA(unsigned short iDVFEA) { return 0.0; }

inline su2double CSolver::GetGlobal_Sens_E(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetGlobal_Sens_Nu(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetGlobal_Sens_Rho(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetGlobal_Sens_Rho_DL(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetGlobal_Sens_EField(unsigned short iEField) { return 0.0; }

inline su2double CSolver::GetGlobal_Sens_DVFEA(unsigned short iDVFEA) { return 0.0; }

inline su2double CSolver::GetVal_Young(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetVal_Poisson(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetVal_Rho(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetVal_Rho_DL(unsigned short iVal) { return 0.0; }

inline unsigned short CSolver::GetnEField(void) { return 0; }

inline unsigned short CSolver::GetnDVFEA(void) { return 0; }

inline void CSolver::ReadDV(CConfig *config) { }

inline su2double CSolver::GetVal_EField(unsigned short iVal) { return 0.0; }

inline su2double CSolver::GetVal_DVFEA(unsigned short iVal) { return 0.0; }

inline su2double* CSolver::GetConstants() { return NULL;}

inline void CSolver::SetTotal_ComboObj(su2double ComboObj) {}

inline su2double CSolver::GetTotal_ComboObj(void) { return 0;}

inline void CSolver::Set_Heatflux_Areas(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Evaluate_ObjFunc(CConfig *config) {};

inline void CSolver::Solve_System(CGeometry *geometry, CConfig *config) { }

inline void CSolver::BC_Euler_Wall(CGeometry      *geometry, 
                                   CSolver        **solver_container, 
                                   CNumerics      *conv_numerics, 
                                   CNumerics      *visc_numerics, 
                                   CConfig        *config, 
                                   unsigned short val_marker) { }

inline void CSolver::BC_Clamped(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_DispDir(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Clamped_Post(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }
                   
inline void CSolver::BC_Normal_Displacement(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }
                                                       
inline void CSolver::BC_Normal_Load(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Dir_Load(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }
                   
inline void CSolver::BC_Sine_Load(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }
                   
inline void CSolver::BC_Damper(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }
                  
inline void CSolver::BC_Deforming(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }
                  
inline void CSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
                  unsigned short val_marker) { }

inline void CSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                         CConfig *config) { }

inline void CSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                           CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                           CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config) { }

inline void CSolver::BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                      CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                       CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_ActDisk(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                CConfig *config, unsigned short val_marker, bool val_inlet_surface) { }

inline void CSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                    CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
                  CConfig *config, unsigned short val_marker) { }
                  
inline void CSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                     CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
                     CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_TurboRiemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
										 CConfig *config, unsigned short val_marker) { }

inline void CSolver::PreprocessBC_Giles(CGeometry *geometry, CConfig *config,
																								CNumerics *conv_numerics,unsigned short marker_flag){}

inline void CSolver::BC_Giles(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
                     CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                     CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                     CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
                      CConfig *config, unsigned short val_marker) { }
                      
inline void CSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                     CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                         CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
                      CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
                      CConfig *config, unsigned short val_marker) { }
                                            
inline void CSolver::BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
                      CConfig *config, unsigned short val_marker) { }
                  
inline void CSolver::BC_Dielec(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
                   CConfig *config, unsigned short val_marker) { }
                                         
inline void CSolver::BC_Electrode(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
                  CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                      CConfig *config, unsigned short val_marker) { }

inline void CSolver::GetPower_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) { }

inline void CSolver::GetOutlet_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) { }

inline void CSolver::GetEllipticSpanLoad_Diff(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                     CConfig *config, unsigned short iMesh, bool Output) { }

inline bool CSolver::FixedCL_Convergence(CConfig *config, bool convergence) { return false; }

inline bool CSolver::GetStart_AoA_FD(void) { return false; }

inline bool CSolver::GetEnd_AoA_FD(void) { return false; }

inline unsigned long CSolver::GetIter_Update_AoA(void) { return 0; }

inline su2double CSolver::GetPrevious_AoA(void) { return 0.0; }

inline su2double CSolver::GetAoA_inc(void) { return 0.0; }

inline void CSolver::SetActDisk_BCThrust(CGeometry *geometry, CSolver **solver_container,
                                         CConfig *config, unsigned short iMesh, bool Output) { }

inline void CSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
							        unsigned short iMesh, unsigned long Iteration) { }

inline void CSolver::CheckTimeSynchronization(CConfig         *config,
                                              const su2double TimeSync,
                                              su2double       &timeEvolved,
                                              bool            &syncTimeReached) {}

inline void CSolver::ProcessTaskList_DG(CGeometry *geometry,  CSolver **solver_container,
                                        CNumerics **numerics, CConfig *config,
                                        unsigned short iMesh) {}

inline void CSolver::ADER_SpaceTimeIntegration(CGeometry *geometry,  CSolver **solver_container,
                                               CNumerics **numerics, CConfig *config,
                                               unsigned short iMesh, unsigned short RunTime_EqSystem) {}

inline void CSolver::ComputeSpatialJacobian(CGeometry *geometry,  CSolver **solver_container,
                                            CNumerics **numerics, CConfig *config,
                                            unsigned short iMesh, unsigned short RunTime_EqSystem) {}

inline void CSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                      unsigned short iMesh) { }

inline void CSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
                      unsigned short iMesh) { }  

inline void CSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
                    CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

inline void CSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
                     CConfig *config, unsigned short iMesh) { }

inline void CSolver::Convective_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                       CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

inline void CSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) { }

inline void CSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) { }

inline void CSolver::SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetUpwind_Ducros_Sensor(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_ActDisk(CSolver **solver_container, CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Nearfield(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Pressure_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Momentum_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Friction_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Buffet_Monitoring(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Heat_Fluxes(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Inviscid_DeltaForces(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Viscous_DeltaForces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Wave_Strength(CGeometry *geometry, CConfig *config) { }

inline void CSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, 
                      CConfig *config, unsigned short iRKStep) { }

inline void CSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                            CConfig *config, unsigned short iRKStep) { }

inline void CSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::ImplicitNewmark_Update(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::ImplicitNewmark_Relaxation(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::GeneralizedAlpha_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::GeneralizedAlpha_UpdateDisp(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::GeneralizedAlpha_UpdateSolution(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::GeneralizedAlpha_UpdateLoads(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Compute_Residual(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
                    unsigned short iMesh) { }

inline void CSolver::SetRes_RMS(unsigned short val_var, su2double val_residual) { Residual_RMS[val_var] = val_residual; }

inline void CSolver::AddRes_RMS(unsigned short val_var, su2double val_residual) { Residual_RMS[val_var] += val_residual; }

inline su2double CSolver::GetRes_RMS(unsigned short val_var) { return Residual_RMS[val_var]; }

inline void CSolver::SetRes_Max(unsigned short val_var, su2double val_residual, unsigned long val_point) { Residual_Max[val_var] = val_residual; Point_Max[val_var] = val_point; }

inline void CSolver::AddRes_Max(unsigned short val_var, su2double val_residual, unsigned long val_point, su2double* val_coord) {
  if (val_residual > Residual_Max[val_var]) {
  Residual_Max[val_var] = val_residual;
  Point_Max[val_var] = val_point;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Point_Max_Coord[val_var][iDim] = val_coord[iDim];
  }
}

inline void CSolver::AddRes_Max(unsigned short val_var, su2double val_residual, unsigned long val_point, const su2double* val_coord) {
  if (val_residual > Residual_Max[val_var]) {
    Residual_Max[val_var] = val_residual;
    Point_Max[val_var] = val_point;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Point_Max_Coord[val_var][iDim] = val_coord[iDim];
  }
}

inline su2double CSolver::GetRes_Max(unsigned short val_var) { return Residual_Max[val_var]; }

inline void CSolver::SetRes_BGS(unsigned short val_var, su2double val_residual) { Residual_BGS[val_var] = val_residual; }

inline void CSolver::AddRes_BGS(unsigned short val_var, su2double val_residual) { Residual_BGS[val_var] += val_residual; }

inline su2double CSolver::GetRes_BGS(unsigned short val_var) { return Residual_BGS[val_var]; }

inline void CSolver::SetRes_Max_BGS(unsigned short val_var, su2double val_residual, unsigned long val_point) { Residual_Max_BGS[val_var] = val_residual; Point_Max_BGS[val_var] = val_point; }

inline void CSolver::AddRes_Max_BGS(unsigned short val_var, su2double val_residual, unsigned long val_point, su2double* val_coord) {
  if (val_residual > Residual_Max_BGS[val_var]) {
  Residual_Max_BGS[val_var] = val_residual;
  Point_Max_BGS[val_var] = val_point;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Point_Max_Coord_BGS[val_var][iDim] = val_coord[iDim];
  }
}

inline su2double CSolver::GetRes_Max_BGS(unsigned short val_var) { return Residual_Max_BGS[val_var]; }

inline su2double CSolver::GetRes_FEM(unsigned short val_var) { return 0.0; }

inline unsigned long CSolver::GetPoint_Max(unsigned short val_var) { return Point_Max[val_var]; }

inline su2double* CSolver::GetPoint_Max_Coord(unsigned short val_var) { return Point_Max_Coord[val_var]; }

inline unsigned long CSolver::GetPoint_Max_BGS(unsigned short val_var) { return Point_Max_BGS[val_var]; }

inline su2double* CSolver::GetPoint_Max_Coord_BGS(unsigned short val_var) { return Point_Max_Coord_BGS[val_var]; }

inline void CSolver::Set_OldSolution(CGeometry *geometry) { base_nodes->Set_OldSolution(); }

inline void CSolver::Set_NewSolution(CGeometry *geometry) { }

inline unsigned short CSolver::GetnVar(void) { return nVar; }

inline unsigned short CSolver::GetnOutputVariables(void) { return nOutputVariables; }

inline unsigned short CSolver::GetnPrimVar(void) { return nPrimVar; }

inline unsigned short CSolver::GetnPrimVarGrad(void) { return nPrimVarGrad; }

inline unsigned short CSolver::GetnSecondaryVar(void) { return nSecondaryVar; }

inline unsigned short CSolver::GetnSecondaryVarGrad(void) { return nSecondaryVarGrad; }

inline su2double CSolver::GetMax_Delta_Time(void) { return Max_Delta_Time; }

inline su2double CSolver::GetMin_Delta_Time(void) { return Min_Delta_Time; }

inline su2double CSolver::GetMax_Delta_Time(unsigned short val_Species) { return 0.0; }

inline su2double CSolver::GetMin_Delta_Time(unsigned short val_Species) { return 0.0; }

inline void CSolver::Copy_Zone_Solution(CSolver ***solver1_solution, CGeometry **solver1_geometry, CConfig *solver1_config, 
                      CSolver ***solver2_solution, CGeometry **solver2_geometry, CConfig *solver2_config) {};

inline CFluidModel* CSolver::GetFluidModel(void) { return NULL;}

inline su2double* CSolver::GetVecSolDOFs(void) {return NULL;}

inline unsigned long CSolver::GetnDOFsGlobal(void) {return 0;}

inline void CSolver::Set_ReferenceGeometry(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_Prestretch(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_ElementProperties(CGeometry *geometry, CConfig *config) { }

inline su2double CSolver::Compute_LoadCoefficient(su2double CurrentTime, su2double RampTime, CConfig *config) { return 0.0; }

inline su2double CSolver::Get_ValCoord(CGeometry *geometry, unsigned long indexNode, unsigned short iDim) {return 0.0;}
                      
inline void CSolver::Compute_StiffMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Compute_MassMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Compute_MassRes(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Compute_NodalStressRes(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Compute_NodalStress(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Compute_DeadLoad(CGeometry *geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }  

inline void CSolver::Compute_IntegrationConstants(CConfig *config) { }

inline void CSolver::SetFSI_ConvValue(unsigned short val_index, su2double val_criteria) { };

inline su2double CSolver::GetFSI_ConvValue(unsigned short val_index) { return 0.0; }

inline void CSolver::RegisterSolution(CGeometry *geometry_container, CConfig *config){}

inline void CSolver::RegisterOutput(CGeometry *geometry_container, CConfig *config){}

inline void CSolver::SetAdjoint_Output(CGeometry *geometry, CConfig *config){}

inline void CSolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config){}

inline void CSolver::RegisterObj_Func(CConfig *config){}

inline void CSolver::SetSurface_Sensitivity(CGeometry *geometry, CConfig *config){}

inline void CSolver::SetSensitivity(CGeometry *geometry, CSolver **solver, CConfig *config){}

inline void CSolver::SetAdj_ObjFunc(CGeometry *geometry, CConfig *config){}

inline unsigned long CSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {return 0;}

inline void CSolver::SetRecording(CGeometry *geometry, CConfig *config){}

inline void CSolver::SetPressure_Inf(su2double p_inf){}

inline void CSolver::SetTemperature_Inf(su2double t_inf){}

inline void CSolver::SetDensity_Inf(su2double rho_inf){}

inline void CSolver::SetVelocity_Inf(unsigned short val_dim, su2double val_velocity) { }

inline void CSolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset){}

inline void CSolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config){}

inline void CSolver::SetFreeStream_Solution(CConfig *config){}

inline su2double* CBaselineSolver_FEM::GetVecSolDOFs(void) {return VecSolDOFs.data();}

inline void CSolver::SetTauWall_WF(CGeometry *geometry, CSolver** solver_container, CConfig* config){}

inline void CSolver::SetNuTilde_WF(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                           CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {}

inline void CEulerSolver::Set_NewSolution(CGeometry *geometry) { nodes->SetSolution_New(); }

inline void CSolver::InitTurboContainers(CGeometry *geometry, CConfig *config){}

inline void CSolver::PreprocessAverage(CSolver **solver, CGeometry *geometry, CConfig *config, unsigned short marker_flag){}

inline void CSolver::TurboAverageProcess(CSolver **solver, CGeometry *geometry, CConfig *config, unsigned short marker_flag){}

inline void CSolver::GatherInOutAverageValues(CConfig *config, CGeometry *geometry){ }

inline su2double CSolver::GetAverageDensity(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double CSolver::GetAveragePressure(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double* CSolver::GetAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan){return NULL;}

inline su2double CSolver::GetAverageNu(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double CSolver::GetAverageKine(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double CSolver::GetAverageOmega(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double CSolver::GetExtAverageNu(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double CSolver::GetExtAverageKine(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline su2double CSolver::GetExtAverageOmega(unsigned short valMarker, unsigned short valSpan){return 0.0;}

inline void CSolver::SetExtAverageDensity(unsigned short valMarker, unsigned short valSpan, su2double valDensity){ }

inline void CSolver::SetExtAveragePressure(unsigned short valMarker, unsigned short valSpan, su2double valPressure){ }

inline void CSolver::SetExtAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan, unsigned short valIndex, su2double valTurboVelocity){ }

inline void CSolver::SetExtAverageNu(unsigned short valMarker, unsigned short valSpan, su2double valNu){ }

inline void CSolver::SetExtAverageKine(unsigned short valMarker, unsigned short valSpan, su2double valKine){ }

inline void CSolver::SetExtAverageOmega(unsigned short valMarker, unsigned short valSpan, su2double valOmega){ }

inline su2double CSolver::GetDensityIn(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetPressureIn(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double* CSolver::GetTurboVelocityIn(unsigned short inMarkerTP, unsigned short valSpan){return NULL;}

inline su2double CSolver::GetDensityOut(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetPressureOut(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double* CSolver::GetTurboVelocityOut(unsigned short inMarkerTP, unsigned short valSpan){return NULL;}

inline su2double CSolver::GetKineIn(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetOmegaIn(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetNuIn(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetKineOut(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetOmegaOut(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline su2double CSolver::GetNuOut(unsigned short inMarkerTP, unsigned short valSpan){return 0;}

inline void CSolver::SetDensityIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetPressureIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetTurboVelocityIn(su2double *value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetDensityOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetPressureOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetTurboVelocityOut(su2double *value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetKineIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetOmegaIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetNuIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetKineOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetOmegaOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetNuOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){ }

inline void CSolver::SetFreeStream_TurboSolution(CConfig *config){ }

inline void CSolver::SetBeta_Parameter(CGeometry *geometry, CSolver **solver_container,
                                     CConfig *config, unsigned short iMesh) { }

inline void CSolver::SetRoe_Dissipation(CGeometry *geometry, CConfig *config) {}

inline void CSolver::SetDES_LengthScale(CSolver** solver, CGeometry *geometry, CConfig *config) { }

inline void CSolver::DeformMesh(CGeometry **geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::SetMesh_Stiffness(CGeometry **geometry, CNumerics **numerics, CConfig *config) { }

inline void CSolver::SetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var, su2double relaxation_factor, su2double val_var) { }

inline su2double CSolver::GetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var) { return 0.0; }

inline void CSolver::ComputeVerificationError(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetImplicitPeriodic(bool val_implicit_periodic) { implicit_periodic = val_implicit_periodic; }

inline void CSolver::SetRotatePeriodic(bool val_rotate_periodic) { rotate_periodic = val_rotate_periodic; }

inline string CSolver::GetSolverName(void) {return SolverName;}

inline su2double CEulerSolver::GetDensity_Inf(void) { return Density_Inf; }

inline su2double CEulerSolver::GetModVelocity_Inf(void) { 
  su2double Vel2 = 0; 
  for (unsigned short iDim = 0; iDim < nDim; iDim++) 
    Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim]; 
  return sqrt(Vel2);
}

inline su2double CEulerSolver::GetDensity_Energy_Inf(void) { return Density_Inf*Energy_Inf; }

inline su2double CEulerSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return Density_Inf*Velocity_Inf[val_dim]; }

inline su2double CEulerSolver::GetVelocity_Inf(unsigned short val_dim) { return Velocity_Inf[val_dim]; }

inline su2double *CEulerSolver::GetVelocity_Inf(void) { return Velocity_Inf; }

inline su2double CEulerSolver::GetPressure_Inf(void) { return Pressure_Inf; }

inline su2double CEulerSolver::GetCPressure(unsigned short val_marker, unsigned long val_vertex) { return CPressure[val_marker][val_vertex]; }

inline su2double CEulerSolver::GetCPressureTarget(unsigned short val_marker, unsigned long val_vertex) { return CPressureTarget[val_marker][val_vertex]; }

inline void CEulerSolver::SetCPressureTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_pressure) { CPressureTarget[val_marker][val_vertex] = val_pressure; }

inline su2double *CEulerSolver::GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex) { return CharacPrimVar[val_marker][val_vertex]; }

inline void CEulerSolver::SetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value) { CharacPrimVar[val_marker][val_vertex][val_var] = val_value; }

inline su2double *CEulerSolver::GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex) { return DonorPrimVar[val_marker][val_vertex]; }

inline void CEulerSolver::SetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value) { DonorPrimVar[val_marker][val_vertex][val_var] = val_value; }

inline su2double CEulerSolver::GetDonorPrimVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var) { return DonorPrimVar[val_marker][val_vertex][val_var]; }

inline unsigned long CEulerSolver::GetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex) { return DonorGlobalIndex[val_marker][val_vertex]; }

inline void CEulerSolver::SetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex, unsigned long val_index) { DonorGlobalIndex[val_marker][val_vertex] = val_index; }

inline su2double CEulerSolver::GetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex) { return ActDisk_DeltaP[val_marker][val_vertex]; }

inline void CEulerSolver::SetActDisk_DeltaP(unsigned short val_marker, unsigned long val_vertex, su2double val_deltap) { ActDisk_DeltaP[val_marker][val_vertex] = val_deltap; }

inline su2double CEulerSolver::GetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex) { return ActDisk_DeltaT[val_marker][val_vertex]; }

inline void CEulerSolver::SetActDisk_DeltaT(unsigned short val_marker, unsigned long val_vertex, su2double val_deltat) { ActDisk_DeltaT[val_marker][val_vertex] = val_deltat; }

inline su2double CEulerSolver::GetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex) { return Inlet_Ttotal[val_marker][val_vertex]; }

inline su2double CEulerSolver::GetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex) { return Inlet_Ptotal[val_marker][val_vertex]; }

inline su2double CEulerSolver::GetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim) { return Inlet_FlowDir[val_marker][val_vertex][val_dim]; }

inline void CEulerSolver::SetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex, su2double val_ttotal) {
  /*--- Since this call can be accessed indirectly using python, do some error
   * checking to prevent segmentation faults ---*/
  if (val_marker >= nMarker)
    SU2_MPI::Error("Out-of-bounds marker index used on inlet.", CURRENT_FUNCTION);
  else if (Inlet_Ttotal == NULL || Inlet_Ttotal[val_marker] == NULL)
    SU2_MPI::Error("Tried to set custom inlet BC on an invalid marker.", CURRENT_FUNCTION);
  else if (val_vertex >= nVertex[val_marker])
    SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
  else
    Inlet_Ttotal[val_marker][val_vertex] = val_ttotal;
}

inline void CEulerSolver::SetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex, su2double val_ptotal) {
  /*--- Since this call can be accessed indirectly using python, do some error
   * checking to prevent segmentation faults ---*/
  if (val_marker >= nMarker)
    SU2_MPI::Error("Out-of-bounds marker index used on inlet.", CURRENT_FUNCTION);
  else if (Inlet_Ptotal == NULL || Inlet_Ptotal[val_marker] == NULL)
    SU2_MPI::Error("Tried to set custom inlet BC on an invalid marker.", CURRENT_FUNCTION);
  else if (val_vertex >= nVertex[val_marker])
    SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
  else
    Inlet_Ptotal[val_marker][val_vertex] = val_ptotal;
}

inline void CEulerSolver::SetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_flowdir) {
  /*--- Since this call can be accessed indirectly using python, do some error
   * checking to prevent segmentation faults ---*/
  if (val_marker >= nMarker)
    SU2_MPI::Error("Out-of-bounds marker index used on inlet.", CURRENT_FUNCTION);
  else if (Inlet_FlowDir == NULL || Inlet_FlowDir[val_marker] == NULL)
      SU2_MPI::Error("Tried to set custom inlet BC on an invalid marker.", CURRENT_FUNCTION);
  else if (val_vertex >= nVertex[val_marker])
    SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
  else
    Inlet_FlowDir[val_marker][val_vertex][val_dim] = val_flowdir;
}

inline su2double CEulerSolver::GetCL_Inv(unsigned short val_marker) { return CL_Inv[val_marker]; }

inline su2double CEulerSolver::GetCD_Inv(unsigned short val_marker) { return CD_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CL(unsigned short val_marker) { return Surface_CL[val_marker]; }

inline su2double CEulerSolver::GetSurface_CD(unsigned short val_marker) { return Surface_CD[val_marker]; }

inline su2double CEulerSolver::GetSurface_CSF(unsigned short val_marker) { return Surface_CSF[val_marker]; }

inline su2double CEulerSolver::GetSurface_CEff(unsigned short val_marker) { return Surface_CEff[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFx(unsigned short val_marker) { return Surface_CFx[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFy(unsigned short val_marker) { return Surface_CFy[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFz(unsigned short val_marker) { return Surface_CFz[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMx(unsigned short val_marker) { return Surface_CMx[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMy(unsigned short val_marker) { return Surface_CMy[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMz(unsigned short val_marker) { return Surface_CMz[val_marker]; }

inline su2double CEulerSolver::GetSurface_CL_Inv(unsigned short val_marker) { return Surface_CL_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CD_Inv(unsigned short val_marker) { return Surface_CD_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CSF_Inv(unsigned short val_marker) { return Surface_CSF_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CEff_Inv(unsigned short val_marker) { return Surface_CEff_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFx_Inv(unsigned short val_marker) { return Surface_CFx_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFy_Inv(unsigned short val_marker) { return Surface_CFy_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFz_Inv(unsigned short val_marker) { return Surface_CFz_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMx_Inv(unsigned short val_marker) { return Surface_CMx_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMy_Inv(unsigned short val_marker) { return Surface_CMy_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMz_Inv(unsigned short val_marker) { return Surface_CMz_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CL_Mnt(unsigned short val_marker) { return Surface_CL_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CD_Mnt(unsigned short val_marker) { return Surface_CD_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CSF_Mnt(unsigned short val_marker) { return Surface_CSF_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CEff_Mnt(unsigned short val_marker) { return Surface_CEff_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFx_Mnt(unsigned short val_marker) { return Surface_CFx_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFy_Mnt(unsigned short val_marker) { return Surface_CFy_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFz_Mnt(unsigned short val_marker) { return Surface_CFz_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMx_Mnt(unsigned short val_marker) { return Surface_CMx_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMy_Mnt(unsigned short val_marker) { return Surface_CMy_Mnt[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMz_Mnt(unsigned short val_marker) { return Surface_CMz_Mnt[val_marker]; }

inline su2double CEulerSolver::GetInflow_MassFlow(unsigned short val_marker) { return Inflow_MassFlow[val_marker]; }

inline su2double CEulerSolver::GetExhaust_MassFlow(unsigned short val_marker) { return Exhaust_MassFlow[val_marker]; }

inline su2double CEulerSolver::GetInflow_Pressure(unsigned short val_marker) { return Inflow_Pressure[val_marker]; }

inline su2double CEulerSolver::GetInflow_Mach(unsigned short val_marker) { return Inflow_Mach[val_marker]; }

inline su2double CEulerSolver::GetCSF_Inv(unsigned short val_marker) { return CSF_Inv[val_marker]; }

inline su2double CEulerSolver::GetCEff_Inv(unsigned short val_marker) { return CEff_Inv[val_marker]; }

inline su2double CEulerSolver::GetTotal_CL() { return Total_CL; }

inline void CEulerSolver::SetTotal_ComboObj(su2double ComboObj) {Total_ComboObj = ComboObj; }

inline su2double CEulerSolver::GetTotal_ComboObj() { return Total_ComboObj; }

inline su2double CEulerSolver::GetTotal_CD() { return Total_CD; }

inline su2double CEulerSolver::GetTotal_NetThrust() { return Total_NetThrust; }

inline su2double CEulerSolver::GetTotal_Power() { return Total_Power; }

inline su2double CEulerSolver::GetTotal_SolidCD() { return Total_SolidCD; }

inline su2double CEulerSolver::GetTotal_ReverseFlow() { return Total_ReverseFlow; }

inline su2double CEulerSolver::GetTotal_MFR() { return Total_MFR; }

inline su2double CEulerSolver::GetTotal_Prop_Eff() { return Total_Prop_Eff; }

inline su2double CEulerSolver::GetTotal_ByPassProp_Eff() { return Total_ByPassProp_Eff; }

inline su2double CEulerSolver::GetTotal_Adiab_Eff() { return Total_Adiab_Eff; }

inline su2double CEulerSolver::GetTotal_Poly_Eff() { return Total_Poly_Eff; }

inline su2double CEulerSolver::GetTotal_IDC_Mach() { return Total_IDC_Mach; }

inline su2double CEulerSolver::GetTotal_DC60() { return Total_DC60; }

inline su2double CEulerSolver::GetTotal_Custom_ObjFunc() { return Total_Custom_ObjFunc; }

inline su2double CEulerSolver::GetTotal_CMx() { return Total_CMx; }

inline su2double CEulerSolver::GetTotal_CMy() { return Total_CMy; }

inline su2double CEulerSolver::GetTotal_CMz() { return Total_CMz; }

inline su2double CEulerSolver::GetTotal_CoPx() { return Total_CoPx; }

inline su2double CEulerSolver::GetTotal_CoPy() { return Total_CoPy; }

inline su2double CEulerSolver::GetTotal_CoPz() { return Total_CoPz; }

inline su2double CEulerSolver::GetTotal_CFx() { return Total_CFx; }

inline su2double CEulerSolver::GetTotal_CFy() { return Total_CFy; }

inline su2double CEulerSolver::GetTotal_CFz() { return Total_CFz; }

inline su2double CEulerSolver::GetTotal_CSF() { return Total_CSF; }

inline su2double CEulerSolver::GetTotal_CEff() { return Total_CEff; }

inline su2double CEulerSolver::GetTotal_CT() { return Total_CT; }

inline void CEulerSolver::SetTotal_CT(su2double val_Total_CT) { Total_CT = val_Total_CT; }

inline su2double CEulerSolver::GetTotal_CQ() { return Total_CQ; }

inline su2double CEulerSolver::GetTotal_HeatFlux() { return Total_Heat; }

inline su2double CEulerSolver::GetTotal_MaxHeatFlux() { return Total_MaxHeat; }

inline void CEulerSolver::SetTotal_CQ(su2double val_Total_CQ) { Total_CQ = val_Total_CQ; }

inline void CEulerSolver::SetTotal_HeatFlux(su2double val_Total_Heat) { Total_Heat = val_Total_Heat; }

inline void CEulerSolver::SetTotal_MaxHeatFlux(su2double val_Total_MaxHeat) { Total_MaxHeat = val_Total_MaxHeat; }

inline su2double CEulerSolver::GetTotal_CMerit() { return Total_CMerit; }

inline su2double CEulerSolver::GetTotal_CEquivArea() { return Total_CEquivArea; }

inline su2double CEulerSolver::GetTotal_AeroCD() { return Total_AeroCD; }

inline su2double CEulerSolver::GetTotal_IDR() { return Total_IDR; }

inline su2double CEulerSolver::GetTotal_IDC() { return Total_IDC; }

inline su2double CEulerSolver::GetTotal_CpDiff() { return Total_CpDiff; }

inline su2double CEulerSolver::GetTotal_HeatFluxDiff() { return Total_HeatFluxDiff; }

inline su2double CEulerSolver::GetTotal_CNearFieldOF() { return Total_CNearFieldOF; }

inline void CEulerSolver::AddTotal_ComboObj(su2double val_obj) {Total_ComboObj +=val_obj;}

inline void CEulerSolver::SetTotal_CEquivArea(su2double val_cequivarea) { Total_CEquivArea = val_cequivarea; }

inline void CEulerSolver::SetTotal_AeroCD(su2double val_aerocd) { Total_AeroCD = val_aerocd; }

inline void CEulerSolver::SetTotal_CpDiff(su2double pressure) { Total_CpDiff = pressure; }

inline void CEulerSolver::SetTotal_HeatFluxDiff(su2double heat) { Total_HeatFluxDiff = heat; }

inline void CEulerSolver::SetTotal_CNearFieldOF(su2double cnearfieldpress) { Total_CNearFieldOF = cnearfieldpress; }

inline void CEulerSolver::SetTotal_CL(su2double val_Total_CL) { Total_CL = val_Total_CL; }

inline void CEulerSolver::SetTotal_CD(su2double val_Total_CD) { Total_CD = val_Total_CD; }

inline void CEulerSolver::SetTotal_NetThrust(su2double val_Total_NetThrust) { Total_NetThrust = val_Total_NetThrust; }

inline void CEulerSolver::SetTotal_Power(su2double val_Total_Power) { Total_Power = val_Total_Power; }

inline void CEulerSolver::SetTotal_SolidCD(su2double val_Total_SolidCD) { Total_SolidCD = val_Total_SolidCD; }

inline void CEulerSolver::SetTotal_ReverseFlow(su2double val_Total_ReverseFlow) { Total_ReverseFlow = val_Total_ReverseFlow; }

inline void CEulerSolver::SetTotal_MFR(su2double val_Total_MFR) { Total_MFR = val_Total_MFR; }

inline void CEulerSolver::SetTotal_Prop_Eff(su2double val_Total_Prop_Eff) { Total_Prop_Eff = val_Total_Prop_Eff; }

inline void CEulerSolver::SetTotal_ByPassProp_Eff(su2double val_Total_ByPassProp_Eff) { Total_ByPassProp_Eff = val_Total_ByPassProp_Eff; }

inline void CEulerSolver::SetTotal_Adiab_Eff(su2double val_Total_Adiab_Eff) { Total_Adiab_Eff = val_Total_Adiab_Eff; }

inline void CEulerSolver::SetTotal_Poly_Eff(su2double val_Total_Poly_Eff) { Total_Poly_Eff = val_Total_Poly_Eff; }

inline void CEulerSolver::SetTotal_IDC(su2double val_Total_IDC) { Total_IDC = val_Total_IDC; }

inline void CEulerSolver::SetTotal_IDC_Mach(su2double val_Total_IDC_Mach) { Total_IDC_Mach = val_Total_IDC_Mach; }

inline void CEulerSolver::SetTotal_IDR(su2double val_Total_IDR) { Total_IDR = val_Total_IDR; }

inline void CEulerSolver::SetTotal_DC60(su2double val_Total_DC60) { Total_DC60 = val_Total_DC60; }

inline void CEulerSolver::SetTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { Total_Custom_ObjFunc = val_total_custom_objfunc*val_weight; }

inline void CEulerSolver::AddTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { Total_Custom_ObjFunc += val_total_custom_objfunc*val_weight; }

inline su2double CEulerSolver::GetAllBound_CL_Inv() { return AllBound_CL_Inv; }

inline su2double CEulerSolver::GetAllBound_CD_Inv() { return AllBound_CD_Inv; }

inline su2double CEulerSolver::GetAllBound_CSF_Inv() { return AllBound_CSF_Inv; }

inline su2double CEulerSolver::GetAllBound_CEff_Inv() { return AllBound_CEff_Inv; }

inline su2double CEulerSolver::GetAllBound_CMx_Inv() { return AllBound_CMx_Inv; }

inline su2double CEulerSolver::GetAllBound_CMy_Inv() { return AllBound_CMy_Inv; }

inline su2double CEulerSolver::GetAllBound_CMz_Inv() { return AllBound_CMz_Inv; }

inline su2double CEulerSolver::GetAllBound_CoPx_Inv() { return AllBound_CoPx_Inv; }

inline su2double CEulerSolver::GetAllBound_CoPy_Inv() { return AllBound_CoPy_Inv; }

inline su2double CEulerSolver::GetAllBound_CoPz_Inv() { return AllBound_CoPz_Inv; }

inline su2double CEulerSolver::GetAllBound_CFx_Inv() { return AllBound_CFx_Inv; }

inline su2double CEulerSolver::GetAllBound_CFy_Inv() { return AllBound_CFy_Inv; }

inline su2double CEulerSolver::GetAllBound_CFz_Inv() { return AllBound_CFz_Inv; }

inline su2double CEulerSolver::GetAllBound_CL_Mnt() { return AllBound_CL_Mnt; }

inline su2double CEulerSolver::GetAllBound_CD_Mnt() { return AllBound_CD_Mnt; }

inline su2double CEulerSolver::GetAllBound_CSF_Mnt() { return AllBound_CSF_Mnt; }

inline su2double CEulerSolver::GetAllBound_CEff_Mnt() { return AllBound_CEff_Mnt; }

inline su2double CEulerSolver::GetAllBound_CMx_Mnt() { return AllBound_CMx_Mnt; }

inline su2double CEulerSolver::GetAllBound_CMy_Mnt() { return AllBound_CMy_Mnt; }

inline su2double CEulerSolver::GetAllBound_CMz_Mnt() { return AllBound_CMz_Mnt; }

inline su2double CEulerSolver::GetAllBound_CoPx_Mnt() { return AllBound_CoPx_Mnt; }

inline su2double CEulerSolver::GetAllBound_CoPy_Mnt() { return AllBound_CoPy_Mnt; }

inline su2double CEulerSolver::GetAllBound_CoPz_Mnt() { return AllBound_CoPz_Mnt; }

inline su2double CEulerSolver::GetAllBound_CFx_Mnt() { return AllBound_CFx_Mnt; }

inline su2double CEulerSolver::GetAllBound_CFy_Mnt() { return AllBound_CFy_Mnt; }

inline su2double CEulerSolver::GetAllBound_CFz_Mnt() { return AllBound_CFz_Mnt; }

inline su2double CEulerSolver::GetAverageDensity(unsigned short valMarker, unsigned short valSpan){return AverageDensity[valMarker][valSpan];}

inline su2double CEulerSolver::GetAveragePressure(unsigned short valMarker, unsigned short valSpan){return AveragePressure[valMarker][valSpan];}

inline su2double* CEulerSolver::GetAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan){return AverageTurboVelocity[valMarker][valSpan];}

inline su2double CEulerSolver::GetAverageNu(unsigned short valMarker, unsigned short valSpan){return AverageNu[valMarker][valSpan];}

inline su2double CEulerSolver::GetAverageKine(unsigned short valMarker, unsigned short valSpan){return AverageKine[valMarker][valSpan];}

inline su2double CEulerSolver::GetAverageOmega(unsigned short valMarker, unsigned short valSpan){return AverageOmega[valMarker][valSpan];}

inline su2double CEulerSolver::GetExtAverageNu(unsigned short valMarker, unsigned short valSpan){return ExtAverageNu[valMarker][valSpan];}

inline su2double CEulerSolver::GetExtAverageKine(unsigned short valMarker, unsigned short valSpan){return ExtAverageKine[valMarker][valSpan];}

inline su2double CEulerSolver::GetExtAverageOmega(unsigned short valMarker, unsigned short valSpan){return ExtAverageOmega[valMarker][valSpan];}

inline void CEulerSolver::SetExtAverageDensity(unsigned short valMarker, unsigned short valSpan, su2double valDensity){ExtAverageDensity[valMarker][valSpan] = valDensity;}

inline void CEulerSolver::SetExtAveragePressure(unsigned short valMarker, unsigned short valSpan, su2double valPressure){ExtAveragePressure[valMarker][valSpan] = valPressure;}

inline void CEulerSolver::SetExtAverageTurboVelocity(unsigned short valMarker, unsigned short valSpan, unsigned short valIndex, su2double valTurboVelocity){ExtAverageTurboVelocity[valMarker][valSpan][valIndex] = valTurboVelocity;}

inline void CEulerSolver::SetExtAverageNu(unsigned short valMarker, unsigned short valSpan, su2double valNu){ExtAverageNu[valMarker][valSpan] = valNu;}

inline void CEulerSolver::SetExtAverageKine(unsigned short valMarker, unsigned short valSpan, su2double valKine){ExtAverageKine[valMarker][valSpan] = valKine;}

inline void CEulerSolver::SetExtAverageOmega(unsigned short valMarker, unsigned short valSpan, su2double valOmega){ExtAverageOmega[valMarker][valSpan] = valOmega;}

inline su2double  CEulerSolver::GetDensityIn(unsigned short inMarkerTP, unsigned short valSpan){return DensityIn[inMarkerTP][valSpan];}

inline su2double  CEulerSolver::GetPressureIn(unsigned short inMarkerTP, unsigned short valSpan){return PressureIn[inMarkerTP][valSpan];}

inline su2double* CEulerSolver::GetTurboVelocityIn(unsigned short inMarkerTP, unsigned short valSpan){return TurboVelocityIn[inMarkerTP][valSpan];}

inline su2double  CEulerSolver::GetDensityOut(unsigned short inMarkerTP, unsigned short valSpan){return DensityOut[inMarkerTP][valSpan];}

inline su2double  CEulerSolver::GetPressureOut(unsigned short inMarkerTP, unsigned short valSpan){return PressureOut[inMarkerTP][valSpan];}

inline su2double* CEulerSolver::GetTurboVelocityOut(unsigned short inMarkerTP, unsigned short valSpan){return TurboVelocityOut[inMarkerTP][valSpan];}

inline su2double CEulerSolver::GetKineIn(unsigned short inMarkerTP, unsigned short valSpan){return KineIn[inMarkerTP][valSpan];}

inline su2double CEulerSolver::GetOmegaIn(unsigned short inMarkerTP, unsigned short valSpan){return OmegaIn[inMarkerTP][valSpan];}

inline su2double CEulerSolver::GetNuIn(unsigned short inMarkerTP, unsigned short valSpan){return NuIn[inMarkerTP][valSpan];}

inline su2double CEulerSolver::GetKineOut(unsigned short inMarkerTP, unsigned short valSpan){return KineOut[inMarkerTP][valSpan];}

inline su2double CEulerSolver::GetOmegaOut(unsigned short inMarkerTP, unsigned short valSpan){return OmegaOut[inMarkerTP][valSpan];}

inline su2double CEulerSolver::GetNuOut(unsigned short inMarkerTP, unsigned short valSpan){return NuOut[inMarkerTP][valSpan];}

inline void CEulerSolver::SetDensityIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){DensityIn[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetPressureIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){PressureIn[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetTurboVelocityIn(su2double *value, unsigned short inMarkerTP, unsigned short valSpan){
  unsigned short iDim;

  for(iDim = 0; iDim < nDim; iDim++)
    TurboVelocityIn[inMarkerTP][valSpan][iDim] = value[iDim];
}

inline void CEulerSolver::SetDensityOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){DensityOut[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetPressureOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){PressureOut[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetTurboVelocityOut(su2double *value, unsigned short inMarkerTP, unsigned short valSpan){
  unsigned short iDim;

  for(iDim = 0; iDim < nDim; iDim++)
    TurboVelocityOut[inMarkerTP][valSpan][iDim] = value[iDim];
}

inline void CEulerSolver::SetKineIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){KineIn[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetOmegaIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){OmegaIn[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetNuIn(su2double value, unsigned short inMarkerTP, unsigned short valSpan){NuIn[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetKineOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){KineOut[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetOmegaOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){OmegaOut[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::SetNuOut(su2double value, unsigned short inMarkerTP, unsigned short valSpan){NuOut[inMarkerTP][valSpan] = value;}

inline void CEulerSolver::ComputeTurboVelocity(su2double *cartesianVelocity, su2double *turboNormal, su2double *turboVelocity, unsigned short marker_flag, unsigned short kind_turb) {

  if ((kind_turb == AXIAL && nDim == 3) || (kind_turb == CENTRIPETAL_AXIAL && marker_flag == OUTFLOW) || (kind_turb == AXIAL_CENTRIFUGAL && marker_flag == INFLOW) ){
    turboVelocity[2] =  turboNormal[0]*cartesianVelocity[0] + cartesianVelocity[1]*turboNormal[1];
    turboVelocity[1] =  turboNormal[0]*cartesianVelocity[1] - turboNormal[1]*cartesianVelocity[0];
    turboVelocity[0] = cartesianVelocity[2];
  }
  else{
    turboVelocity[0] =  turboNormal[0]*cartesianVelocity[0] + cartesianVelocity[1]*turboNormal[1];
    turboVelocity[1] =  turboNormal[0]*cartesianVelocity[1] - turboNormal[1]*cartesianVelocity[0];
    if (marker_flag == INFLOW){
      turboVelocity[0] *= -1.0;
      turboVelocity[1] *= -1.0;
    }
    if(nDim == 3)
      turboVelocity[2] = cartesianVelocity[2];
  }
}

inline void CEulerSolver::ComputeBackVelocity(su2double *turboVelocity, su2double *turboNormal, su2double *cartesianVelocity, unsigned short marker_flag, unsigned short kind_turb){

  if ((kind_turb == AXIAL && nDim == 3) || (kind_turb == CENTRIPETAL_AXIAL && marker_flag == OUTFLOW) || (kind_turb == AXIAL_CENTRIFUGAL && marker_flag == INFLOW)){
    cartesianVelocity[0] = turboVelocity[2]*turboNormal[0] - turboVelocity[1]*turboNormal[1];
    cartesianVelocity[1] = turboVelocity[2]*turboNormal[1] + turboVelocity[1]*turboNormal[0];
    cartesianVelocity[2] = turboVelocity[0];
  }
  else{
    cartesianVelocity[0] =  turboVelocity[0]*turboNormal[0] - turboVelocity[1]*turboNormal[1];
    cartesianVelocity[1] =  turboVelocity[0]*turboNormal[1] + turboVelocity[1]*turboNormal[0];

    if (marker_flag == INFLOW){
      cartesianVelocity[0] *= -1.0;
      cartesianVelocity[1] *= -1.0;
    }

    if(nDim == 3)
      cartesianVelocity[2] = turboVelocity[2];
  }
}


inline CFluidModel* CEulerSolver::GetFluidModel(void) { return FluidModel;}

inline void CEulerSolver::SetPressure_Inf(su2double p_inf) {Pressure_Inf = p_inf;}

inline void CEulerSolver::SetTemperature_Inf(su2double t_inf) {Temperature_Inf = t_inf;}

inline bool CEulerSolver::GetStart_AoA_FD(void) { return Start_AoA_FD; }

inline bool CEulerSolver::GetEnd_AoA_FD(void) { return End_AoA_FD; }

inline unsigned long CEulerSolver::GetIter_Update_AoA(void) { return Iter_Update_AoA; }

inline su2double CEulerSolver::GetPrevious_AoA(void) { return AoA_Prev; }

inline su2double CEulerSolver::GetAoA_inc(void) { return AoA_inc; }

inline su2double CNSSolver::GetViscosity_Inf(void) { return Viscosity_Inf; }

inline su2double CNSSolver::GetTke_Inf(void) { return Tke_Inf; }

inline su2double CNSSolver::GetSurface_HF_Visc(unsigned short val_marker) { return Surface_HF_Visc[val_marker]; }

inline su2double CNSSolver::GetSurface_MaxHF_Visc(unsigned short val_marker) { return Surface_MaxHF_Visc[val_marker]; }

inline su2double CNSSolver::GetCL_Visc(unsigned short val_marker) { return CL_Visc[val_marker]; }

inline su2double CNSSolver::GetCSF_Visc(unsigned short val_marker) { return CSF_Visc[val_marker]; }

inline su2double CNSSolver::GetCD_Visc(unsigned short val_marker) { return CD_Visc[val_marker]; }

inline su2double CNSSolver::GetAllBound_CL_Visc() { return AllBound_CL_Visc; }

inline su2double CNSSolver::GetAllBound_CD_Visc() { return AllBound_CD_Visc; }

inline su2double CNSSolver::GetAllBound_CSF_Visc() { return AllBound_CSF_Visc; }

inline su2double CNSSolver::GetAllBound_CEff_Visc() { return AllBound_CEff_Visc; }

inline su2double CNSSolver::GetAllBound_CMx_Visc() { return AllBound_CMx_Visc; }

inline su2double CNSSolver::GetAllBound_CMy_Visc() { return AllBound_CMy_Visc; }

inline su2double CNSSolver::GetAllBound_CMz_Visc() { return AllBound_CMz_Visc; }

inline su2double CNSSolver::GetAllBound_CoPx_Visc() { return AllBound_CoPx_Visc; }

inline su2double CNSSolver::GetAllBound_CoPy_Visc() { return AllBound_CoPy_Visc; }

inline su2double CNSSolver::GetAllBound_CoPz_Visc() { return AllBound_CoPz_Visc; }

inline su2double CNSSolver::GetAllBound_CFx_Visc() { return AllBound_CFx_Visc; }

inline su2double CNSSolver::GetAllBound_CFy_Visc() { return AllBound_CFy_Visc; }

inline su2double CNSSolver::GetAllBound_CFz_Visc() { return AllBound_CFz_Visc; }

inline su2double CNSSolver::GetTotal_Buffet_Metric() { return Total_Buffet_Metric; }

inline su2double CNSSolver::GetSurface_CL_Visc(unsigned short val_marker) { return Surface_CL_Visc[val_marker]; }

inline su2double CNSSolver::GetSurface_CD_Visc(unsigned short val_marker) { return Surface_CD_Visc[val_marker]; }

inline su2double CNSSolver::GetSurface_CSF_Visc(unsigned short val_marker) { return Surface_CSF_Visc[val_marker]; }

inline su2double CNSSolver::GetSurface_CEff_Visc(unsigned short val_marker) { return Surface_CEff_Visc[val_marker]; }

inline su2double CNSSolver::GetSurface_CFx_Visc(unsigned short val_marker) { return Surface_CFx_Visc[val_marker]; }

inline su2double CNSSolver::GetSurface_CFy_Visc(unsigned short val_marker) { return Surface_CFy_Visc[val_marker]; }

inline su2double CNSSolver::GetSurface_CFz_Visc(unsigned short val_marker) { return Surface_CFz_Visc[val_marker]; }

inline su2double CNSSolver::GetSurface_CMx_Visc(unsigned short val_marker) { return Surface_CMx_Visc[val_marker]; }

inline su2double CNSSolver::GetSurface_CMy_Visc(unsigned short val_marker) { return Surface_CMy_Visc[val_marker]; }

inline su2double CNSSolver::GetSurface_CMz_Visc(unsigned short val_marker) { return Surface_CMz_Visc[val_marker]; }

inline su2double CNSSolver::GetSurface_Buffet_Metric(unsigned short val_marker) { return Surface_Buffet_Metric[val_marker]; }

inline su2double CNSSolver::GetCSkinFriction(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim) { return CSkinFriction[val_marker][val_dim][val_vertex]; }

inline su2double CNSSolver::GetHeatFlux(unsigned short val_marker, unsigned long val_vertex) { return HeatFlux[val_marker][val_vertex]; }

inline su2double CNSSolver::GetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex) { return HeatFluxTarget[val_marker][val_vertex]; }

inline void CNSSolver::SetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_heat) { HeatFluxTarget[val_marker][val_vertex] = val_heat; }

inline su2double CNSSolver::GetBuffetSensor(unsigned short val_marker, unsigned long val_vertex) { return Buffet_Sensor[val_marker][val_vertex]; }

inline su2double CNSSolver::GetYPlus(unsigned short val_marker, unsigned long val_vertex) { return YPlus[val_marker][val_vertex]; }

inline su2double CNSSolver::GetStrainMag_Max(void) { return StrainMag_Max; }

inline su2double CNSSolver::GetOmega_Max(void) { return Omega_Max; }

inline void CNSSolver::SetStrainMag_Max(su2double val_strainmag_max) { StrainMag_Max = val_strainmag_max; }

inline void CNSSolver::SetOmega_Max(su2double val_omega_max) { Omega_Max = val_omega_max; }

inline su2double CNSSolver::GetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var) { return HeatConjugateVar[val_marker][val_vertex][pos_var]; }

inline void CNSSolver::SetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var, su2double relaxation_factor, su2double val_var) {
  HeatConjugateVar[val_marker][val_vertex][pos_var] = relaxation_factor*val_var + (1.0-relaxation_factor)*HeatConjugateVar[val_marker][val_vertex][pos_var]; }

inline CFluidModel* CFEM_DG_EulerSolver::GetFluidModel(void) { return FluidModel;}

inline su2double* CFEM_DG_EulerSolver::GetVecSolDOFs(void) {return VecSolDOFs.data();}

inline unsigned long CFEM_DG_EulerSolver::GetnDOFsGlobal(void) {return nDOFsGlobal;}

inline su2double CFEM_DG_EulerSolver::GetDensity_Inf(void) { return Density_Inf; }

inline su2double CFEM_DG_EulerSolver::GetModVelocity_Inf(void) {
  su2double Vel2 = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  return sqrt(Vel2);
}

inline su2double CFEM_DG_EulerSolver::GetDensity_Energy_Inf(void) { return Density_Inf*Energy_Inf; }

inline su2double CFEM_DG_EulerSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return Density_Inf*Velocity_Inf[val_dim]; }

inline su2double CFEM_DG_EulerSolver::GetVelocity_Inf(unsigned short val_dim) { return Velocity_Inf[val_dim]; }

inline su2double *CFEM_DG_EulerSolver::GetVelocity_Inf(void) { return Velocity_Inf; }

inline su2double CFEM_DG_EulerSolver::GetPressure_Inf(void) { return Pressure_Inf; }

inline su2double CFEM_DG_EulerSolver::GetCL_Inv(unsigned short val_marker) { return CL_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetCMz_Inv(unsigned short val_marker) { return CMz_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetCD_Inv(unsigned short val_marker) { return CD_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CL(unsigned short val_marker) { return Surface_CL[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CD(unsigned short val_marker) { return Surface_CD[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CSF(unsigned short val_marker) { return Surface_CSF[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CEff(unsigned short val_marker) { return Surface_CEff[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CFx(unsigned short val_marker) { return Surface_CFx[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CFy(unsigned short val_marker) { return Surface_CFy[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CFz(unsigned short val_marker) { return Surface_CFz[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CMx(unsigned short val_marker) { return Surface_CMx[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CMy(unsigned short val_marker) { return Surface_CMy[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CMz(unsigned short val_marker) { return Surface_CMz[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CL_Inv(unsigned short val_marker) { return Surface_CL_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CD_Inv(unsigned short val_marker) { return Surface_CD_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CSF_Inv(unsigned short val_marker) { return Surface_CSF_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CEff_Inv(unsigned short val_marker) { return Surface_CEff_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CFx_Inv(unsigned short val_marker) { return Surface_CFx_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CFy_Inv(unsigned short val_marker) { return Surface_CFy_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CFz_Inv(unsigned short val_marker) { return Surface_CFz_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CMx_Inv(unsigned short val_marker) { return Surface_CMx_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CMy_Inv(unsigned short val_marker) { return Surface_CMy_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetSurface_CMz_Inv(unsigned short val_marker) { return Surface_CMz_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetCSF_Inv(unsigned short val_marker) { return CSF_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetCEff_Inv(unsigned short val_marker) { return CEff_Inv[val_marker]; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CL() { return Total_CL; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CD() { return Total_CD; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CMx() { return Total_CMx; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CMy() { return Total_CMy; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CMz() { return Total_CMz; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CFx() { return Total_CFx; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CFy() { return Total_CFy; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CFz() { return Total_CFz; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CSF() { return Total_CSF; }

inline su2double CFEM_DG_EulerSolver::GetTotal_CEff() { return Total_CEff; }

inline void CFEM_DG_EulerSolver::SetTotal_CL(su2double val_Total_CL) { Total_CL = val_Total_CL; }

inline void CFEM_DG_EulerSolver::SetTotal_CD(su2double val_Total_CD) { Total_CD = val_Total_CD; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CL_Inv() { return AllBound_CL_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CD_Inv() { return AllBound_CD_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CSF_Inv() { return AllBound_CSF_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CEff_Inv() { return AllBound_CEff_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CMx_Inv() { return AllBound_CMx_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CMy_Inv() { return AllBound_CMy_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CMz_Inv() { return AllBound_CMz_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CFx_Inv() { return AllBound_CFx_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CFy_Inv() { return AllBound_CFy_Inv; }

inline su2double CFEM_DG_EulerSolver::GetAllBound_CFz_Inv() { return AllBound_CFz_Inv; }

inline void CFEM_DG_EulerSolver::SetPressure_Inf(su2double p_inf){Pressure_Inf = p_inf;}

inline void CFEM_DG_EulerSolver::SetTemperature_Inf(su2double t_inf){Temperature_Inf = t_inf;}

inline void CFEM_DG_EulerSolver::BC_HeatFlux_Wall(CConfig                  *config,
                                                  const unsigned long      surfElemBeg,
                                                  const unsigned long      surfElemEnd,
                                                  const CSurfaceElementFEM *surfElem,
                                                  su2double                *resFaces,
                                                  CNumerics                *conv_numerics,
                                                  unsigned short           val_marker,
                                                  su2double                *workArray) {}

inline void CFEM_DG_EulerSolver::BC_Isothermal_Wall(CConfig                  *config,
                                                    const unsigned long      surfElemBeg,
                                                    const unsigned long      surfElemEnd,
                                                    const CSurfaceElementFEM *surfElem,
                                                    su2double                *resFaces,
                                                    CNumerics                *conv_numerics,
                                                    unsigned short           val_marker,
                                                    su2double                *workArray) {}

inline su2double CFEM_DG_NSSolver::GetViscosity_Inf(void) { return Viscosity_Inf; }

inline su2double CFEM_DG_NSSolver::GetTke_Inf(void) { return Tke_Inf; }

inline su2double CFEM_DG_NSSolver::GetCL_Visc(unsigned short val_marker) { return CL_Visc[val_marker]; }

inline su2double CFEM_DG_NSSolver::GetCMz_Visc(unsigned short val_marker) { return CMz_Visc[val_marker]; }

inline su2double CFEM_DG_NSSolver::GetCSF_Visc(unsigned short val_marker) { return CSF_Visc[val_marker]; }

inline su2double CFEM_DG_NSSolver::GetCD_Visc(unsigned short val_marker) { return CD_Visc[val_marker]; }

inline su2double CFEM_DG_NSSolver::GetAllBound_CL_Visc() { return AllBound_CL_Visc; }

inline su2double CFEM_DG_NSSolver::GetAllBound_CSF_Visc() { return AllBound_CSF_Visc; }

inline su2double CFEM_DG_NSSolver::GetAllBound_CD_Visc() { return AllBound_CD_Visc; }

inline su2double CFEM_DG_NSSolver::GetStrainMag_Max(void) { return StrainMag_Max; }

inline su2double CFEM_DG_NSSolver::GetOmega_Max(void) { return Omega_Max; }

inline void CFEM_DG_NSSolver::SetStrainMag_Max(su2double val_strainmag_max) { StrainMag_Max = val_strainmag_max; }

inline void CFEM_DG_NSSolver::SetOmega_Max(su2double val_omega_max) { Omega_Max = val_omega_max; }

inline su2double CAdjEulerSolver::GetCSensitivity(unsigned short val_marker, unsigned long val_vertex) { return CSensitivity[val_marker][val_vertex]; }

inline void CAdjEulerSolver::SetCSensitivity(unsigned short val_marker, unsigned long val_vertex, su2double val_sensitivity) { CSensitivity[val_marker][val_vertex] = val_sensitivity; }

inline unsigned long CAdjEulerSolver::GetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex) { return DonorGlobalIndex[val_marker][val_vertex]; }

inline void CAdjEulerSolver::SetDonorGlobalIndex(unsigned short val_marker, unsigned long val_vertex, unsigned long val_index) { DonorGlobalIndex[val_marker][val_vertex] = val_index; }

inline void CAdjEulerSolver::SetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var, su2double val_value) { DonorAdjVar[val_marker][val_vertex][val_var] = val_value; }

inline su2double CAdjEulerSolver::GetTotal_Sens_Geo() { return Total_Sens_Geo; }

inline su2double CAdjEulerSolver::GetTotal_Sens_Mach() { return Total_Sens_Mach; }

inline su2double CAdjEulerSolver::GetTotal_Sens_AoA() { return Total_Sens_AoA; }

inline su2double CAdjEulerSolver::GetTotal_Sens_Press() { return Total_Sens_Press; }

inline su2double CAdjEulerSolver::GetTotal_Sens_Temp() { return Total_Sens_Temp; }

inline su2double CAdjEulerSolver::GetTotal_Sens_BPress() { return Total_Sens_BPress; }

inline su2double CAdjEulerSolver::GetPsiRho_Inf(void) { return PsiRho_Inf; }

inline su2double CAdjEulerSolver::GetPsiE_Inf(void) { return PsiE_Inf; }

inline su2double *CAdjEulerSolver::GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex) { return DonorAdjVar[val_marker][val_vertex]; }

inline su2double CAdjEulerSolver::GetDonorAdjVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_var) { return DonorAdjVar[val_marker][val_vertex][val_var]; }

inline su2double CAdjEulerSolver::GetPhi_Inf(unsigned short val_dim) { return Phi_Inf[val_dim]; }

inline void CSolver::RefGeom_Sensitivity(CGeometry *geometry, CSolver **solver_container, CConfig *config){ }

inline void CSolver::DE_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){ }

inline void CSolver::Stiffness_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){ }

inline unsigned short CSolver::Get_iElem_iDe(unsigned long iElem){ return 0; }

inline void CSolver::Set_DV_Val(su2double val_EField, unsigned short i_DV){ }

inline su2double CSolver::Get_DV_Val(unsigned short i_DV){ return 0.0; }

inline su2double CSolver::Get_val_I(void){ return 0.0; }

inline su2double CSolver::Get_MassMatrix(unsigned long iPoint, unsigned long jPoint, unsigned short iVar, unsigned short jVar){ return 0.0; }

inline su2double CIncEulerSolver::GetDensity_Inf(void) { return Density_Inf; }

inline su2double CIncEulerSolver::GetModVelocity_Inf(void) {
  su2double Vel2 = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  return sqrt(Vel2);
}

inline CFluidModel* CIncEulerSolver::GetFluidModel(void) { return FluidModel;}

inline su2double CIncEulerSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return Density_Inf*Velocity_Inf[val_dim]; }

inline su2double CIncEulerSolver::GetVelocity_Inf(unsigned short val_dim) { return Velocity_Inf[val_dim]; }

inline void CIncEulerSolver::SetVelocity_Inf(unsigned short val_dim, su2double val_velocity) { Velocity_Inf[val_dim] = val_velocity; }

inline su2double *CIncEulerSolver::GetVelocity_Inf(void) { return Velocity_Inf; }

inline su2double CIncEulerSolver::GetPressure_Inf(void) { return Pressure_Inf; }

inline su2double CIncEulerSolver::GetTemperature_Inf(void) { return Temperature_Inf; }

inline su2double CIncEulerSolver::GetCPressure(unsigned short val_marker, unsigned long val_vertex) { return CPressure[val_marker][val_vertex]; }

inline su2double CIncEulerSolver::GetCPressureTarget(unsigned short val_marker, unsigned long val_vertex) { return CPressureTarget[val_marker][val_vertex]; }

inline void CIncEulerSolver::SetCPressureTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_pressure) { CPressureTarget[val_marker][val_vertex] = val_pressure; }

inline su2double *CIncEulerSolver::GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex) { return CharacPrimVar[val_marker][val_vertex]; }

inline su2double CIncEulerSolver::GetInlet_Ttotal(unsigned short val_marker, unsigned long val_vertex) { return Inlet_Ttotal[val_marker][val_vertex]; }

inline su2double CIncEulerSolver::GetInlet_Ptotal(unsigned short val_marker, unsigned long val_vertex) { return Inlet_Ptotal[val_marker][val_vertex]; }

inline su2double CIncEulerSolver::GetInlet_FlowDir(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim) { return Inlet_FlowDir[val_marker][val_vertex][val_dim]; }

inline su2double CIncEulerSolver::GetCD_Inv(unsigned short val_marker) { return CD_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CL(unsigned short val_marker) { return Surface_CL[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CD(unsigned short val_marker) { return Surface_CD[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CSF(unsigned short val_marker) { return Surface_CSF[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CEff(unsigned short val_marker) { return Surface_CEff[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFx(unsigned short val_marker) { return Surface_CFx[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFy(unsigned short val_marker) { return Surface_CFy[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFz(unsigned short val_marker) { return Surface_CFz[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMx(unsigned short val_marker) { return Surface_CMx[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMy(unsigned short val_marker) { return Surface_CMy[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMz(unsigned short val_marker) { return Surface_CMz[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CL_Inv(unsigned short val_marker) { return Surface_CL_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CD_Inv(unsigned short val_marker) { return Surface_CD_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CSF_Inv(unsigned short val_marker) { return Surface_CSF_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CEff_Inv(unsigned short val_marker) { return Surface_CEff_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFx_Inv(unsigned short val_marker) { return Surface_CFx_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFy_Inv(unsigned short val_marker) { return Surface_CFy_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFz_Inv(unsigned short val_marker) { return Surface_CFz_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMx_Inv(unsigned short val_marker) { return Surface_CMx_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMy_Inv(unsigned short val_marker) { return Surface_CMy_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMz_Inv(unsigned short val_marker) { return Surface_CMz_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetCSF_Inv(unsigned short val_marker) { return CSF_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetCEff_Inv(unsigned short val_marker) { return CEff_Inv[val_marker]; }

inline su2double CIncEulerSolver::GetTotal_CL() { return Total_CL; }

inline su2double CIncEulerSolver::GetTotal_CD() { return Total_CD; }

inline su2double CIncEulerSolver::GetTotal_CMx() { return Total_CMx; }

inline su2double CIncEulerSolver::GetTotal_CMy() { return Total_CMy; }

inline su2double CIncEulerSolver::GetTotal_CMz() { return Total_CMz; }

inline su2double CIncEulerSolver::GetTotal_CoPx() { return Total_CoPx; }

inline su2double CIncEulerSolver::GetTotal_CoPy() { return Total_CoPy; }

inline su2double CIncEulerSolver::GetTotal_CoPz() { return Total_CoPz; }

inline su2double CIncEulerSolver::GetTotal_CFx() { return Total_CFx; }

inline su2double CIncEulerSolver::GetTotal_CFy() { return Total_CFy; }

inline su2double CIncEulerSolver::GetTotal_CFz() { return Total_CFz; }

inline su2double CIncEulerSolver::GetTotal_CSF() { return Total_CSF; }

inline su2double CIncEulerSolver::GetTotal_CEff() { return Total_CEff; }

inline su2double CIncEulerSolver::GetTotal_CT() { return Total_CT; }

inline void CIncEulerSolver::SetTotal_CT(su2double val_Total_CT) { Total_CT = val_Total_CT; }

inline su2double CIncEulerSolver::GetTotal_CQ() { return Total_CQ; }

inline su2double CIncEulerSolver::GetTotal_HeatFlux() { return Total_Heat; }

inline su2double CIncEulerSolver::GetTotal_MaxHeatFlux() { return Total_MaxHeat; }

inline void CIncEulerSolver::SetTotal_CQ(su2double val_Total_CQ) { Total_CQ = val_Total_CQ; }

inline void CIncEulerSolver::SetTotal_HeatFlux(su2double val_Total_Heat) { Total_Heat = val_Total_Heat; }

inline void CIncEulerSolver::SetTotal_MaxHeatFlux(su2double val_Total_MaxHeat) { Total_MaxHeat = val_Total_MaxHeat; }

inline su2double CIncEulerSolver::GetTotal_CMerit() { return Total_CMerit; }

inline su2double CIncEulerSolver::GetTotal_CpDiff() { return Total_CpDiff; }

inline su2double CIncEulerSolver::GetTotal_HeatFluxDiff() { return Total_HeatFluxDiff; }

inline void CIncEulerSolver::SetTotal_CpDiff(su2double pressure) { Total_CpDiff = pressure; }

inline void CIncEulerSolver::SetTotal_HeatFluxDiff(su2double heat) { Total_HeatFluxDiff = heat; }

inline void CIncEulerSolver::SetTotal_CD(su2double val_Total_CD) { Total_CD = val_Total_CD; }

inline su2double CIncEulerSolver::GetTotal_Custom_ObjFunc() { return Total_Custom_ObjFunc; }

inline void CIncEulerSolver::SetTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { Total_Custom_ObjFunc = val_total_custom_objfunc*val_weight; }

inline void CIncEulerSolver::AddTotal_Custom_ObjFunc(su2double val_total_custom_objfunc, su2double val_weight) { Total_Custom_ObjFunc += val_total_custom_objfunc*val_weight; }

inline su2double CIncEulerSolver::GetAllBound_CL_Inv() { return AllBound_CL_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CD_Inv() { return AllBound_CD_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CSF_Inv() { return AllBound_CSF_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CEff_Inv() { return AllBound_CEff_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CMx_Inv() { return AllBound_CMx_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CMy_Inv() { return AllBound_CMy_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CMz_Inv() { return AllBound_CMz_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CoPx_Inv() { return AllBound_CoPx_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CoPy_Inv() { return AllBound_CoPy_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CoPz_Inv() { return AllBound_CoPz_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CFx_Inv() { return AllBound_CFx_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CFy_Inv() { return AllBound_CFy_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CFz_Inv() { return AllBound_CFz_Inv; }

inline su2double CIncEulerSolver::GetAllBound_CL_Mnt() { return AllBound_CL_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CD_Mnt() { return AllBound_CD_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CSF_Mnt() { return AllBound_CSF_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CEff_Mnt() { return AllBound_CEff_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CMx_Mnt() { return AllBound_CMx_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CMy_Mnt() { return AllBound_CMy_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CMz_Mnt() { return AllBound_CMz_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CoPx_Mnt() { return AllBound_CoPx_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CoPy_Mnt() { return AllBound_CoPy_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CoPz_Mnt() { return AllBound_CoPz_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CFx_Mnt() { return AllBound_CFx_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CFy_Mnt() { return AllBound_CFy_Mnt; }

inline su2double CIncEulerSolver::GetAllBound_CFz_Mnt() { return AllBound_CFz_Mnt; }

inline su2double CIncEulerSolver::GetSurface_CL_Mnt(unsigned short val_marker) { return Surface_CL_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CD_Mnt(unsigned short val_marker) { return Surface_CD_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CSF_Mnt(unsigned short val_marker) { return Surface_CSF_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CEff_Mnt(unsigned short val_marker) { return Surface_CEff_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFx_Mnt(unsigned short val_marker) { return Surface_CFx_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFy_Mnt(unsigned short val_marker) { return Surface_CFy_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CFz_Mnt(unsigned short val_marker) { return Surface_CFz_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMx_Mnt(unsigned short val_marker) { return Surface_CMx_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMy_Mnt(unsigned short val_marker) { return Surface_CMy_Mnt[val_marker]; }

inline su2double CIncEulerSolver::GetSurface_CMz_Mnt(unsigned short val_marker) { return Surface_CMz_Mnt[val_marker]; }

inline void CIncEulerSolver::SetPressure_Inf(su2double p_inf){Pressure_Inf = p_inf;}

inline void CIncEulerSolver::SetTemperature_Inf(su2double t_inf){Temperature_Inf = t_inf;}

inline void CIncEulerSolver::SetDensity_Inf(su2double rho_inf){Density_Inf = rho_inf;}

inline void CIncEulerSolver::SetTotal_ComboObj(su2double ComboObj) {Total_ComboObj = ComboObj; }

inline su2double CIncEulerSolver::GetTotal_ComboObj() { return Total_ComboObj; }

inline su2double CIncNSSolver::GetViscosity_Inf(void) { return Viscosity_Inf; }

inline su2double CIncNSSolver::GetTke_Inf(void) { return Tke_Inf; }

inline su2double CIncNSSolver::GetSurface_HF_Visc(unsigned short val_marker) { return Surface_HF_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_MaxHF_Visc(unsigned short val_marker) { return Surface_MaxHF_Visc[val_marker]; }

inline su2double CIncNSSolver::GetCL_Visc(unsigned short val_marker) { return CL_Visc[val_marker]; }

inline su2double CIncNSSolver::GetCSF_Visc(unsigned short val_marker) { return CSF_Visc[val_marker]; }

inline su2double CIncNSSolver::GetCD_Visc(unsigned short val_marker) { return CD_Visc[val_marker]; }

inline su2double CIncNSSolver::GetAllBound_CL_Visc() { return AllBound_CL_Visc; }

inline su2double CIncNSSolver::GetAllBound_CSF_Visc() { return AllBound_CSF_Visc; }

inline su2double CIncNSSolver::GetAllBound_CD_Visc() { return AllBound_CD_Visc; }

inline su2double CIncNSSolver::GetAllBound_CEff_Visc() { return AllBound_CEff_Visc; }

inline su2double CIncNSSolver::GetAllBound_CMx_Visc() { return AllBound_CMx_Visc; }

inline su2double CIncNSSolver::GetAllBound_CMy_Visc() { return AllBound_CMy_Visc; }

inline su2double CIncNSSolver::GetAllBound_CMz_Visc() { return AllBound_CMz_Visc; }

inline su2double CIncNSSolver::GetAllBound_CoPx_Visc() { return AllBound_CoPx_Visc; }

inline su2double CIncNSSolver::GetAllBound_CoPy_Visc() { return AllBound_CoPy_Visc; }

inline su2double CIncNSSolver::GetAllBound_CoPz_Visc() { return AllBound_CoPz_Visc; }

inline su2double CIncNSSolver::GetAllBound_CFx_Visc() { return AllBound_CFx_Visc; }

inline su2double CIncNSSolver::GetAllBound_CFy_Visc() { return AllBound_CFy_Visc; }

inline su2double CIncNSSolver::GetAllBound_CFz_Visc() { return AllBound_CFz_Visc; }

inline su2double CIncNSSolver::GetSurface_CL_Visc(unsigned short val_marker) { return Surface_CL_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CD_Visc(unsigned short val_marker) { return Surface_CD_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CSF_Visc(unsigned short val_marker) { return Surface_CSF_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CEff_Visc(unsigned short val_marker) { return Surface_CEff_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CFx_Visc(unsigned short val_marker) { return Surface_CFx_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CFy_Visc(unsigned short val_marker) { return Surface_CFy_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CFz_Visc(unsigned short val_marker) { return Surface_CFz_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CMx_Visc(unsigned short val_marker) { return Surface_CMx_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CMy_Visc(unsigned short val_marker) { return Surface_CMy_Visc[val_marker]; }

inline su2double CIncNSSolver::GetSurface_CMz_Visc(unsigned short val_marker) { return Surface_CMz_Visc[val_marker]; }

inline su2double CIncNSSolver::GetCSkinFriction(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim) { return CSkinFriction[val_marker][val_dim][val_vertex]; }

inline su2double CIncNSSolver::GetHeatFlux(unsigned short val_marker, unsigned long val_vertex) { return HeatFlux[val_marker][val_vertex]; }

inline su2double CIncNSSolver::GetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex) { return HeatFluxTarget[val_marker][val_vertex]; }

inline void CIncNSSolver::SetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_heat) { HeatFluxTarget[val_marker][val_vertex] = val_heat; }

inline su2double CIncNSSolver::GetYPlus(unsigned short val_marker, unsigned long val_vertex) { return YPlus[val_marker][val_vertex]; }

inline su2double CIncNSSolver::GetStrainMag_Max(void) { return StrainMag_Max; }

inline su2double CIncNSSolver::GetOmega_Max(void) { return Omega_Max; }

inline void CIncNSSolver::SetStrainMag_Max(su2double val_strainmag_max) { StrainMag_Max = val_strainmag_max; }

inline void CIncNSSolver::SetOmega_Max(su2double val_omega_max) { Omega_Max = val_omega_max; }

inline su2double CIncNSSolver::GetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var) { return HeatConjugateVar[val_marker][val_vertex][pos_var]; }

inline void CIncNSSolver::SetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var, su2double relaxation_factor, su2double val_var) {
  HeatConjugateVar[val_marker][val_vertex][pos_var] = relaxation_factor*val_var + (1.0-relaxation_factor)*HeatConjugateVar[val_marker][val_vertex][pos_var]; }

inline su2double CHeatSolverFVM::GetTotal_HeatFlux() { return Total_HeatFlux; }

inline su2double CHeatSolverFVM::GetHeatFlux(unsigned short val_marker, unsigned long val_vertex) { return HeatFlux[val_marker][val_vertex]; }

inline su2double CHeatSolverFVM::GetTotal_AvgTemperature() { return Total_AverageT; }

inline su2double CHeatSolverFVM::GetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var) { return ConjugateVar[val_marker][val_vertex][pos_var]; }

inline void CHeatSolverFVM::SetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var, su2double relaxation_factor, su2double val_var) {
  ConjugateVar[val_marker][val_vertex][pos_var] = relaxation_factor*val_var + (1.0-relaxation_factor)*ConjugateVar[val_marker][val_vertex][pos_var]; }

inline su2double CFEASolver::Get_MassMatrix(unsigned long iPoint, unsigned long jPoint, unsigned short iVar, unsigned short jVar){ 
  return MassMatrix.GetBlock(iPoint, jPoint, iVar, jVar); }

inline unsigned short CFEASolver::Get_iElem_iDe(unsigned long iElem){ return iElem_iDe[iElem]; }

inline su2double CFEASolver::GetRes_FEM(unsigned short val_var) { return Conv_Check[val_var]; }

inline su2double CFEASolver::GetTotal_CFEA() { return Total_CFEA; }

inline void CFEASolver::SetTotal_CFEA(su2double cfea) { Total_CFEA = cfea; }

inline void CFEASolver::SetTotal_OFRefGeom(su2double val_ofrefgeom) { Total_OFRefGeom = val_ofrefgeom; }

inline void CFEASolver::SetTotal_OFRefNode(su2double val_ofrefnode) { Total_OFRefNode = val_ofrefnode; }

inline su2double CFEASolver::GetWAitken_Dyn(void) { return WAitken_Dyn; }

inline su2double CFEASolver::GetWAitken_Dyn_tn1(void) { return WAitken_Dyn_tn1; }

inline void CFEASolver::SetWAitken_Dyn(su2double waitk) { WAitken_Dyn = waitk; }

inline void CFEASolver::SetWAitken_Dyn_tn1(su2double waitk_tn1) { WAitken_Dyn_tn1 = waitk_tn1; }

inline void CFEASolver::SetLoad_Increment(su2double val_loadIncrement) { loadIncrement = val_loadIncrement; }

inline su2double CFEASolver::GetLoad_Increment(void) { return loadIncrement; }

inline void CFEASolver::SetFSI_ConvValue(unsigned short val_index, su2double val_criteria) { FSI_Conv[val_index] = val_criteria; }

inline su2double CFEASolver::GetFSI_ConvValue(unsigned short val_index){ return FSI_Conv[val_index]; }

inline su2double CFEASolver::GetTotal_OFRefGeom(void){ return Total_OFRefGeom; }

inline su2double CFEASolver::GetTotal_OFRefNode(void){ return Total_OFRefNode; }

inline su2double CFEASolver::GetTotal_OFVolFrac(void){ return Total_OFVolFrac; }

inline su2double CFEASolver::GetTotal_OFCompliance(void){ return Total_OFCompliance; }

inline bool CFEASolver::IsElementBased(void){ return element_based; }

inline void CFEASolver::SetForceCoeff(su2double val_forcecoeff_history) { ForceCoeff = val_forcecoeff_history; }

inline void CFEASolver::SetRelaxCoeff(su2double val_relaxecoeff_history) { RelaxCoeff = val_relaxecoeff_history; }

inline void CFEASolver::SetFSI_Residual(su2double val_FSI_residual) { FSI_Residual = val_FSI_residual; }

inline su2double CFEASolver::GetForceCoeff(void) { return ForceCoeff; }

inline su2double CFEASolver::GetRelaxCoeff(void) { return RelaxCoeff; }

inline su2double CFEASolver::GetFSI_Residual(void) { return FSI_Residual; }

inline su2double CFEASolver::Get_ValCoord(CGeometry *geometry, unsigned long indexNode, unsigned short iDim) {return geometry->node[indexNode]->GetCoord(iDim);}

inline void CSolver::SetAdjoint_OutputMesh(CGeometry *geometry, CConfig *config) {}

inline void CSolver::ExtractAdjoint_Geometry(CGeometry *geometry, CConfig *config) {}

inline void CSolver::ExtractAdjoint_CrossTerm(CGeometry *geometry, CConfig *config) {}

inline void CSolver::ExtractAdjoint_CrossTerm_Geometry(CGeometry *geometry, CConfig *config) {}

inline void CSolver::ExtractAdjoint_CrossTerm_Geometry_Flow(CGeometry *geometry, CConfig *config) {}

inline void CSolver::SetMesh_Recording(CGeometry **geometry, CVolumetricMovement *grid_movement, CConfig *config) {}

inline su2double CDiscAdjSolver::GetTotal_Sens_Geo() { return Total_Sens_Geo; }

inline su2double CDiscAdjSolver::GetTotal_Sens_Mach() { return Total_Sens_Mach; }

inline su2double CDiscAdjSolver::GetTotal_Sens_AoA() { return Total_Sens_AoA; }

inline su2double CDiscAdjSolver::GetTotal_Sens_Press() { return Total_Sens_Press; }

inline su2double CDiscAdjSolver::GetTotal_Sens_Temp() { return Total_Sens_Temp; }

inline su2double CDiscAdjSolver::GetTotal_Sens_BPress() { return Total_Sens_BPress; }

inline su2double CDiscAdjSolver::GetTotal_Sens_Density() { return Total_Sens_Density; }

inline su2double CDiscAdjSolver::GetTotal_Sens_ModVel() { return Total_Sens_ModVel; }

inline su2double CDiscAdjSolver::GetCSensitivity(unsigned short val_marker, unsigned long val_vertex) { return CSensitivity[val_marker][val_vertex]; }

inline void CEulerSolver::SetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index, su2double component){ 
  SlidingState[val_marker][val_vertex][val_state][donor_index] = component; 
}

inline void CIncEulerSolver::SetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index, su2double component){
  SlidingState[val_marker][val_vertex][val_state][donor_index] = component;
}

inline void CSolver::SetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index, su2double component){ }

inline su2double CEulerSolver::GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index) { return SlidingState[val_marker][val_vertex][val_state][donor_index]; }

inline su2double CIncEulerSolver::GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index) { return SlidingState[val_marker][val_vertex][val_state][donor_index]; }

inline su2double CSolver::GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index) { return 0; }

inline int CEulerSolver::GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex){ return SlidingStateNodes[val_marker][val_vertex]; }

inline int CIncEulerSolver::GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex){ return SlidingStateNodes[val_marker][val_vertex]; }

inline int CSolver::GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex){ return 0; }

inline void CEulerSolver::SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value){ SlidingStateNodes[val_marker][val_vertex] = value; }

inline void CIncEulerSolver::SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value){ SlidingStateNodes[val_marker][val_vertex] = value; }

inline void CSolver::SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value){}

inline void CSolver::SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex){}

inline void CEulerSolver::SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex){	
  int iVar;

  for( iVar = 0; iVar < nPrimVar+1; iVar++){
    if( SlidingState[val_marker][val_vertex][iVar] != NULL )
      delete [] SlidingState[val_marker][val_vertex][iVar];
  }

  for( iVar = 0; iVar < nPrimVar+1; iVar++)
    SlidingState[val_marker][val_vertex][iVar] = new su2double[ GetnSlidingStates(val_marker, val_vertex) ];
}


inline void CIncEulerSolver::SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex){
  int iVar;

  for( iVar = 0; iVar < nPrimVar+1; iVar++){
    if( SlidingState[val_marker][val_vertex][iVar] != NULL )
      delete [] SlidingState[val_marker][val_vertex][iVar];
  }

  for( iVar = 0; iVar < nPrimVar+1; iVar++)
    SlidingState[val_marker][val_vertex][iVar] = new su2double[ GetnSlidingStates(val_marker, val_vertex) ];
}



inline void CTurbSolver::SetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index, su2double component){ 
  SlidingState[val_marker][val_vertex][val_state][donor_index] = component; 
}

inline int CTurbSolver::GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex){ return SlidingStateNodes[val_marker][val_vertex]; }

inline void CTurbSolver::SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex){	
  int iVar;

  for( iVar = 0; iVar < nVar+1; iVar++){
    if( SlidingState[val_marker][val_vertex][iVar] != NULL )
      delete [] SlidingState[val_marker][val_vertex][iVar];
  }

  for( iVar = 0; iVar < nVar+1; iVar++)
    SlidingState[val_marker][val_vertex][iVar] = new su2double[ GetnSlidingStates(val_marker, val_vertex) ];
}

inline void CTurbSolver::SetnSlidingStates(unsigned short val_marker, unsigned long val_vertex, int value){ SlidingStateNodes[val_marker][val_vertex] = value; }

inline su2double CTurbSolver::GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index) { return SlidingState[val_marker][val_vertex][val_state][donor_index]; }

inline void CTurbSolver::SetInlet_TurbVar(unsigned short val_marker, unsigned long val_vertex, unsigned short val_dim, su2double val_turb_var) {
  /*--- Since this call can be accessed indirectly using python, do some error
   * checking to prevent segmentation faults ---*/
  if (val_marker >= nMarker)
    SU2_MPI::Error("Out-of-bounds marker index used on inlet.", CURRENT_FUNCTION);
  else if (Inlet_TurbVars == NULL || Inlet_TurbVars[val_marker] == NULL)
    SU2_MPI::Error("Tried to set custom inlet BC on an invalid marker.", CURRENT_FUNCTION);
  else if (val_vertex >= nVertex[val_marker])
    SU2_MPI::Error("Out-of-bounds vertex index used on inlet.", CURRENT_FUNCTION);
  else if (val_dim >= nVar)
    SU2_MPI::Error("Out-of-bounds index used for inlet turbulence variable.", CURRENT_FUNCTION);
  else
    Inlet_TurbVars[val_marker][val_vertex][val_dim] = val_turb_var;
}

inline void CTurbSASolver::SetFreeStream_Solution(CConfig *config) {
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) nodes->SetSolution(iPoint, 0, nu_tilde_Inf);
}

inline su2double CTurbSASolver::GetNuTilde_Inf(void) { return nu_tilde_Inf; }

inline void CTurbSSTSolver::SetFreeStream_Solution(CConfig *config){
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
    nodes->SetSolution(iPoint, 0, kine_Inf);
    nodes->SetSolution(iPoint, 1, omega_Inf);
  }
}

inline su2double CTurbSSTSolver::GetTke_Inf(void) { return kine_Inf; }

inline su2double CTurbSSTSolver::GetOmega_Inf(void) { return omega_Inf; }

inline su2double CDiscAdjFEASolver::GetTotal_Sens_E(unsigned short iVal) { return Total_Sens_E[iVal]; }

inline su2double CDiscAdjFEASolver::GetTotal_Sens_Nu(unsigned short iVal) { return Total_Sens_Nu[iVal]; }

inline su2double CDiscAdjFEASolver::GetTotal_Sens_Rho(unsigned short iVal) { return Total_Sens_Rho[iVal]; }

inline su2double CDiscAdjFEASolver::GetTotal_Sens_Rho_DL(unsigned short iVal) { return Total_Sens_Rho_DL[iVal]; }

inline su2double CDiscAdjFEASolver::GetTotal_Sens_EField(unsigned short iEField) { return Total_Sens_EField[iEField]; }

inline su2double CDiscAdjFEASolver::GetTotal_Sens_DVFEA(unsigned short iDVFEA) { return Total_Sens_DV[iDVFEA]; }

inline su2double CDiscAdjFEASolver::GetGlobal_Sens_E(unsigned short iVal) { return Global_Sens_E[iVal]; }

inline su2double CDiscAdjFEASolver::GetGlobal_Sens_Nu(unsigned short iVal) { return Global_Sens_Nu[iVal]; }

inline su2double CDiscAdjFEASolver::GetGlobal_Sens_Rho(unsigned short iVal) { return Global_Sens_Rho[iVal]; }

inline su2double CDiscAdjFEASolver::GetGlobal_Sens_Rho_DL(unsigned short iVal) { return Global_Sens_Rho_DL[iVal]; }

inline su2double CDiscAdjFEASolver::GetGlobal_Sens_EField(unsigned short iEField) { return Global_Sens_EField[iEField]; }

inline su2double CDiscAdjFEASolver::GetGlobal_Sens_DVFEA(unsigned short iDVFEA) { return Global_Sens_DV[iDVFEA]; }

inline su2double CDiscAdjFEASolver::GetVal_Young(unsigned short iVal) { return E_i[iVal]; }

inline su2double CDiscAdjFEASolver::GetVal_Poisson(unsigned short iVal) { return Nu_i[iVal]; }

inline su2double CDiscAdjFEASolver::GetVal_Rho(unsigned short iVal) { return Rho_i[iVal]; }

inline su2double CDiscAdjFEASolver::GetVal_Rho_DL(unsigned short iVal) { return Rho_DL_i[iVal]; }

inline unsigned short CDiscAdjFEASolver::GetnEField(void) { return nEField; }

inline unsigned short CDiscAdjFEASolver::GetnDVFEA(void) { return nDV; }

inline su2double CDiscAdjFEASolver::GetVal_EField(unsigned short iVal) { return EField[iVal]; }

inline su2double CDiscAdjFEASolver::GetVal_DVFEA(unsigned short iVal) { return DV_Val[iVal]; }

inline void CSolver::SetDualTime_Mesh(void){ }

inline vector<string> CSolver::GetSolutionFields(){return fields;}
