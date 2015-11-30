/*!
 * \file solver_structure.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
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

inline void CSolver::SetIterLinSolver(unsigned short val_iterlinsolver) { IterLinSolver = val_iterlinsolver; }

inline unsigned short CSolver::GetnSpecies(void) { return 0; }

inline void CSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Primitive(CGeometry *geometry, CConfig *config) { }

//inline void CSolver::Set_MPI_Secondary(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Primitive_Limiter(CGeometry *geometry, CConfig *config) { }

//inline void CSolver::Set_MPI_Secondary_Limiter(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetNondimensionalization(CGeometry *geometry, CConfig *config, unsigned short iMesh) { }

inline unsigned short CSolver::GetIterLinSolver(void) { return IterLinSolver; }

inline su2double CSolver::GetCSensitivity(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep, 
																		 unsigned short iMesh, unsigned short RunTime_EqSystem) { }
																		 
inline void CSolver::SetFreeSurface_Distance(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry, CGeometry **flow_geometry, CConfig *fea_config, CConfig *flow_config, CNumerics *fea_numerics) { }

inline void CSolver::GetSurface_Pressure(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) { }

inline void CSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter) { }
  
inline void CSolver::SetFlow_Displacement(CGeometry **flow_geometry, CVolumetricMovement *flow_grid_movement, CConfig *flow_config, CConfig *fea_config, CGeometry **fea_geometry, CSolver ***fea_solution) { }

inline void CSolver::SetStruct_Displacement(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) { }

inline void CSolver::PredictStruct_Displacement(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) { }

inline void CSolver::ComputeAitken_Coefficient(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution, unsigned long iFSIIter) { }

inline void CSolver::SetAitken_Relaxation(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) { }

inline void CSolver::Update_StructSolution(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) { }

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

inline void CSolver::Set_MPI_Primitive_Gradient(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPrimitive_Limiter_MPI(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) { }

//inline void CSolver::SetSecondary_Gradient_GG(CGeometry *geometry, CConfig *config) { }

//inline void CSolver::SetSecondary_Gradient_LS(CGeometry *geometry, CConfig *config) { }

//inline void CSolver::Set_MPI_Secondary_Gradient(CGeometry *geometry, CConfig *config) { }

//inline void CSolver::SetSecondary_Limiter_MPI(CGeometry *geometry, CConfig *config) { }

//inline void CSolver::SetSecondary_Limiter(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPreconditioner(CConfig *config, unsigned long iPoint) { }

inline void CSolver::SetDistance(CGeometry *geometry, CConfig *config) { };

inline su2double CSolver::GetCLift_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCMz_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCDrag_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CLift(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CDrag(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CSideForce(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CEff(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFx(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFy(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFz(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMx(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMy(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMz(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CLift_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CDrag_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CSideForce_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CEff_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFx_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFy_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CFz_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMx_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMy_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetSurface_CMz_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetInflow_MassFlow(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetExhaust_MassFlow(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetInflow_Pressure(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetInflow_Mach(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCSideForce_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCEff_Inv(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCLift_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCMz_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCSideForce_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetCDrag_Visc(unsigned short val_marker) { return 0; }

inline su2double CSolver::GetAllBound_CLift_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CDrag_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CSideForce_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CEff_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CMx_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CMy_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CMz_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CFx_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CFy_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CFz_Inv() { return 0; }

inline su2double CSolver::GetAllBound_CLift_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CSideForce_Visc() { return 0; }

inline su2double CSolver::GetAllBound_CDrag_Visc() { return 0; }

inline void CSolver::SetForceProj_Vector(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::SetIntBoundary_Jump(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline su2double CSolver::GetTotal_CLift() { return 0; }

inline su2double CSolver::GetTotal_CDrag() { return 0; }

inline su2double CSolver::GetTotal_CMx() { return 0; }

inline su2double CSolver::GetTotal_CMy() { return 0; }

inline su2double CSolver::GetTotal_CMz() { return 0; }

inline su2double CSolver::GetTotal_CFx() { return 0; }

inline su2double CSolver::GetTotal_CFy() { return 0; }

inline su2double CSolver::GetTotal_CFz() { return 0; }

inline su2double CSolver::GetTotal_CSideForce() { return 0; }

inline su2double CSolver::GetTotal_CEff() { return 0; }

inline su2double CSolver::GetTotal_CT() { return 0; }

inline void CSolver::SetTotal_CT(su2double val_Total_CT) { }

inline su2double CSolver::GetTotal_CQ() { return 0; }

inline su2double CSolver::GetTotal_HeatFlux() { return 0; }

inline su2double CSolver::GetTotal_MaxHeatFlux() { return 0; }

inline su2double CSolver::Get_PressureDrag() { return 0; }

inline su2double CSolver::Get_ViscDrag() { return 0; }

inline void CSolver::SetTotal_CQ(su2double val_Total_CQ) { }

inline void CSolver::SetTotal_HeatFlux(su2double val_Total_Heat) { }

inline void CSolver::SetTotal_MaxHeatFlux(su2double val_Total_Heat) { }

inline su2double CSolver::GetTotal_CMerit() { return 0; }

inline su2double CSolver::GetTotal_CEquivArea() { return 0; }

inline su2double CSolver::GetTotal_CpDiff() { return 0; }

inline su2double CSolver::GetTotal_HeatFluxDiff() { return 0; }

inline su2double CSolver::GetTotal_CFreeSurface() { return 0; }

inline su2double CSolver::GetTotal_CFEA() { return 0; }

inline su2double CSolver::GetTotal_CNearFieldOF() { return 0; }

inline void CSolver::SetTotal_CEquivArea(su2double val_cequivarea) { }

inline void CSolver::SetTotal_CpDiff(su2double val_pressure) { }

inline void CSolver::SetTotal_HeatFluxDiff(su2double val_heat) { }

inline void CSolver::SetTotal_CFEA(su2double val_cfea) { }

inline su2double CSolver::GetWAitken_Dyn(void) { return 0; }

inline su2double CSolver::GetWAitken_Dyn_tn1(void) { return 0; }

inline void CSolver::SetWAitken_Dyn(su2double waitk) {  }

inline void CSolver::SetWAitken_Dyn_tn1(su2double waitk_tn1) {  }

inline void CSolver::SetTotal_CFreeSurface(su2double val_freesurface) { }

inline void CSolver::SetTotal_CNearFieldOF(su2double val_cnearfieldpress) { }

inline su2double CSolver::GetTotal_CWave() { return 0; }

inline su2double CSolver::GetTotal_CHeat() { return 0; }

inline void CSolver::SetTotal_CLift(su2double val_Total_CLift) { }

inline void CSolver::SetTotal_CDrag(su2double val_Total_CDrag) { }

inline su2double CSolver::GetCPressure(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetCPressureTarget(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline void CSolver::SetCPressureTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_pressure) { }

inline void CSolver::SetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_heat) { }

inline su2double *CSolver::GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetCSkinFriction(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetHeatFlux(unsigned short val_marker, unsigned long val_vertex) { return 0; }

inline su2double CSolver::GetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex) { return 0; }

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

inline su2double CSolver::GetTke_Inf(void) { return 0; }

inline su2double* CSolver::GetConstants() {return NULL;}

inline su2double CSolver::GetOneD_TotalPress(void) {return 0;}

inline void CSolver::SetOneD_TotalPress(su2double AveragePressure) { }

inline su2double CSolver::GetOneD_Mach(void) {return 0;}

inline void CSolver::SetOneD_Mach(su2double AverageMach) { }

inline su2double CSolver::GetOneD_Temp(void) {return 0;}

inline void CSolver::SetOneD_Temp(su2double AverageTemperature) { }

inline su2double CSolver::GetOneD_MassFlowRate(void) {return 0;}

inline void CSolver::SetOneD_MassFlowRate(su2double MassFlowRate) { }

inline su2double CSolver::GetOneD_FluxAvgPress(void) {return 0;}

inline void CSolver::SetOneD_FluxAvgPress(su2double PressureRef) { }

inline su2double CSolver::GetOneD_FluxAvgDensity(void) {return 0;}

inline void CSolver::SetOneD_FluxAvgDensity(su2double DensityRef) { }

inline su2double CSolver::GetOneD_FluxAvgVelocity(void) {return 0;}

inline void CSolver::SetOneD_FluxAvgVelocity(su2double VelocityRef) { }

inline su2double CSolver::GetOneD_FluxAvgEntalpy(void) {return 0;}

inline void CSolver::SetOneD_FluxAvgEntalpy(su2double EnthalpyRef) { }

inline su2double CSolver::GetAveragedDensity(unsigned short valMarker){ return 0;}

inline su2double CSolver::GetAveragedPressure(unsigned short valMarker){ return 0;}

inline su2double CSolver::GetAveragedEnthalpy(unsigned short valMarker){ return 0;}

inline su2double CSolver::GetAveragedEntropy(unsigned short valMarker){ return 0;}

inline su2double* CSolver::GetAveragedVelocity(unsigned short valMarker){ return 0;}

inline su2double* CSolver::GetAveragedGridVelocity(unsigned short valMarker){ return 0;}

inline su2double CSolver::GetAveragedNormalVelocity(unsigned short valMarker){ return 0;}

inline su2double CSolver::GetAveragedTangVelocity(unsigned short valMarker){ return 0;}

inline su2double CSolver::GetAveragedTotTemperature(unsigned short valMarker){return 0;}

inline su2double CSolver::GetAveragedTotPressure(unsigned short valMarker){return 0;}

inline su2double CSolver::GetMassFlow(unsigned short valMarker){return 0;}

inline su2double CSolver::GetFlowAngle(unsigned short valMarker){return 0;}

inline su2double CSolver::GetAveragedMach(unsigned short valMarker){return 0;}

inline su2double CSolver::GetAveragedNormalMach(unsigned short valMarker){return 0;}

inline su2double CSolver::GetTotalPressureLoss(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetKineticEnergyLoss(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetTotalTotalEfficiency(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetTotalStaticEfficiency(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetEulerianWork(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetTotalEnthalpyIn(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetFlowAngleIn(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetFlowAngleOut(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetMassFlowIn(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetMassFlowOut(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetMachIn(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetMachOut(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetNormalMachIn(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetNormalMachOut(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetEnthalpyOut(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetVelocityOutIs(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetPressureOut(unsigned short inMarkerTP){return 0;}

inline su2double CSolver::GetPressureRatio(unsigned short inMarkerTP){return 0;}

inline void CSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }


inline void CSolver::BC_Clamped(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }

inline void CSolver::BC_Clamped_Post(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }
									 
inline void CSolver::BC_Normal_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }
									 									 
inline void CSolver::BC_Flow_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }
									 
inline void CSolver::BC_Normal_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }

inline void CSolver::BC_Dir_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }
									 
inline void CSolver::BC_Sine_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }									 
                   
inline void CSolver::BC_Pressure(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
									 unsigned short val_marker) { }
                  
inline void CSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }
									
inline void CSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
								  unsigned short val_marker) { }

inline void CSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									CConfig *config) { }
                  
inline void CSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									CConfig *config) { }

inline void CSolver::BC_ActDisk_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                               CConfig *config) { }
										
inline void CSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
								    CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
									CConfig *config, unsigned short val_marker) { }
									
inline void CSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
										 CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
										 CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_NonReflecting(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker){ }

inline void CSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
										 CConfig *config, unsigned short val_marker) { }
                     
inline void CSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
										  CConfig *config, unsigned short val_marker) { }
                      
inline void CSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
										 CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                         CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
										  CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Engine_Bleed(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                      CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
										  CConfig *config, unsigned short val_marker) { }
																						
inline void CSolver::BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
										  CConfig *config, unsigned short val_marker) { }
								  
inline void CSolver::BC_Dielec(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									 CConfig *config, unsigned short val_marker) { }
									 										  
inline void CSolver::BC_Electrode(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									CConfig *config, unsigned short val_marker) { }

inline void CSolver::Mixing_Process(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker){}

inline void CSolver::MPIMixing_Process(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short marker_flag){}

inline void CSolver::MixedOut_Average (su2double val_init_pressure, su2double *val_Averaged_Flux, su2double *val_normal, su2double *pressure_mix, su2double *density_mix){}

inline void CSolver::MixedOut_Root_Function(su2double *pressure, su2double *val_Averaged_Flux, su2double *val_normal, su2double *valfunc, su2double *density){}

inline void CSolver::Boundary_Fourier(CGeometry *geometry, CSolver **solver_container, CConfig *config,
		                               unsigned short val_Marker, vector<std::complex<su2double> > &c4k,signed long &nboundaryvertex){}

inline void CSolver::Boundary_Fourier(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker,
		                              vector<std::complex<su2double> >& c2k,vector<std::complex<su2double> >& c3k,signed long& nboundaryvertex) {}
inline void CSolver::SetExtAveragedValue(CSolver *solver_container, unsigned short intMarker,  unsigned short extMarker){ }

inline void CSolver::GetEngine_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) { }

inline void CSolver::GetActuatorDisk_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) { }

inline void CSolver::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                     CConfig *config, unsigned short iMesh, bool Output) { }

inline void CSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
							        unsigned short iMesh, unsigned long Iteration) { }	
							        
inline void CSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
							        unsigned short iMesh) { }	

inline void CSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
							        unsigned short iMesh) { }	

inline void CSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
										CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

inline void CSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									   CConfig *config, unsigned short iMesh) { }

inline void CSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) { }

inline void CSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) { }

inline void CSolver::SetDissipation_Switch(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Dissipation_Switch(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::TurboPerformance(CSolver *solver, CConfig *config, unsigned short inMarker,  unsigned short outMarker, unsigned short Kind_TurboPerf, unsigned short inMarkerTP ){ }

inline void CSolver::MPITurboPerformance(CConfig *config){ }

inline void CSolver::StoreTurboPerformance(CSolver *solver, unsigned short inMarkerTP ){ }

inline void CSolver::Viscous_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Inviscid_DeltaForces(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Viscous_DeltaForces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Wave_Strength(CGeometry *geometry, CConfig *config) { }

inline void CSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, 
											CConfig *config, unsigned short iRKStep) { }

inline void CSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

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

inline su2double CSolver::GetRes_Max(unsigned short val_var) { return Residual_Max[val_var]; }

inline unsigned long CSolver::GetPoint_Max(unsigned short val_var) { return Point_Max[val_var]; }

inline su2double* CSolver::GetPoint_Max_Coord(unsigned short val_var) { return Point_Max_Coord[val_var]; }

inline void CSolver::Set_OldSolution(CGeometry *geometry) {
	for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) 
		node[iPoint]->Set_OldSolution(); // The loop should be over nPoints 
                                     //  to guarantee that the boundaries are
                                     //  well updated
}

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

inline CFluidModel* CEulerSolver::GetFluidModel(void) { return FluidModel;}
										  
inline void CSolver::Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline void CSolver::Compute_StiffMassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }										  

inline void CSolver::Compute_StiffMassDampMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline void CSolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }	

inline void CSolver::Compute_IntegrationConstants(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline void CSolver::SetSolution_time_n(CGeometry *geometry, CConfig *config) { }										  
										  
inline su2double CEulerSolver::GetDensity_Inf(void) { return Density_Inf; }

inline su2double CEulerSolver::GetModVelocity_Inf(void) { 
	su2double Vel2 = 0; 
	for (unsigned short iDim = 0; iDim < nDim; iDim++) 
		Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim]; 
	return sqrt(Vel2);
}

inline void CSolver::SetFSI_ConvValue(unsigned short val_index, su2double val_criteria) { };

inline su2double CSolver::GetFSI_ConvValue(unsigned short val_index) { return 0.0; }

inline su2double CEulerSolver::GetDensity_Energy_Inf(void) { return Density_Inf*Energy_Inf; }

inline su2double CEulerSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return Density_Inf*Velocity_Inf[val_dim]; }

inline su2double CEulerSolver::GetVelocity_Inf(unsigned short val_dim) { return Velocity_Inf[val_dim]; }

inline su2double *CEulerSolver::GetVelocity_Inf(void) { return Velocity_Inf; }

inline su2double CEulerSolver::GetPressure_Inf(void) { return Pressure_Inf; }

inline su2double CEulerSolver::GetCPressure(unsigned short val_marker, unsigned long val_vertex) { return CPressure[val_marker][val_vertex]; }

inline su2double CEulerSolver::GetCPressureTarget(unsigned short val_marker, unsigned long val_vertex) { return CPressureTarget[val_marker][val_vertex]; }

inline void CEulerSolver::SetCPressureTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_pressure) { CPressureTarget[val_marker][val_vertex] = val_pressure; }

inline su2double *CEulerSolver::GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex) { return CharacPrimVar[val_marker][val_vertex]; }

inline su2double CEulerSolver::GetCLift_Inv(unsigned short val_marker) { return CLift_Inv[val_marker]; }

inline su2double CEulerSolver::GetCMz_Inv(unsigned short val_marker) { return CMz_Inv[val_marker]; }

inline su2double CEulerSolver::GetCDrag_Inv(unsigned short val_marker) { return CDrag_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CLift(unsigned short val_marker) { return Surface_CLift[val_marker]; }

inline su2double CEulerSolver::GetSurface_CDrag(unsigned short val_marker) { return Surface_CDrag[val_marker]; }

inline su2double CEulerSolver::GetSurface_CSideForce(unsigned short val_marker) { return Surface_CSideForce[val_marker]; }

inline su2double CEulerSolver::GetSurface_CEff(unsigned short val_marker) { return Surface_CEff[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFx(unsigned short val_marker) { return Surface_CFx[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFy(unsigned short val_marker) { return Surface_CFy[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFz(unsigned short val_marker) { return Surface_CFz[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMx(unsigned short val_marker) { return Surface_CMx[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMy(unsigned short val_marker) { return Surface_CMy[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMz(unsigned short val_marker) { return Surface_CMz[val_marker]; }

inline su2double CEulerSolver::GetSurface_CLift_Inv(unsigned short val_marker) { return Surface_CLift_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CDrag_Inv(unsigned short val_marker) { return Surface_CDrag_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CSideForce_Inv(unsigned short val_marker) { return Surface_CSideForce_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CEff_Inv(unsigned short val_marker) { return Surface_CEff_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFx_Inv(unsigned short val_marker) { return Surface_CFx_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFy_Inv(unsigned short val_marker) { return Surface_CFy_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CFz_Inv(unsigned short val_marker) { return Surface_CFz_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMx_Inv(unsigned short val_marker) { return Surface_CMx_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMy_Inv(unsigned short val_marker) { return Surface_CMy_Inv[val_marker]; }

inline su2double CEulerSolver::GetSurface_CMz_Inv(unsigned short val_marker) { return Surface_CMz_Inv[val_marker]; }

inline su2double CEulerSolver::GetInflow_MassFlow(unsigned short val_marker) { return Inflow_MassFlow[val_marker]; }

inline su2double CEulerSolver::GetExhaust_MassFlow(unsigned short val_marker) { return Exhaust_MassFlow[val_marker]; }

inline su2double CEulerSolver::GetInflow_Pressure(unsigned short val_marker) { return Inflow_Pressure[val_marker]; }

inline su2double CEulerSolver::GetInflow_Mach(unsigned short val_marker) { return Inflow_Mach[val_marker]; }

inline su2double CEulerSolver::GetCSideForce_Inv(unsigned short val_marker) { return CSideForce_Inv[val_marker]; }

inline su2double CEulerSolver::GetCEff_Inv(unsigned short val_marker) { return CEff_Inv[val_marker]; }

inline su2double CEulerSolver::GetTotal_CLift() { return Total_CLift; }

inline su2double CEulerSolver::GetTotal_CDrag() { return Total_CDrag; }

inline su2double CEulerSolver::GetTotal_CMx() { return Total_CMx; }

inline su2double CEulerSolver::GetTotal_CMy() { return Total_CMy; }

inline su2double CEulerSolver::GetTotal_CMz() { return Total_CMz; }

inline su2double CEulerSolver::GetTotal_CFx() { return Total_CFx; }

inline su2double CEulerSolver::GetTotal_CFy() { return Total_CFy; }

inline su2double CEulerSolver::GetTotal_CFz() { return Total_CFz; }

inline su2double CEulerSolver::GetTotal_CSideForce() { return Total_CSideForce; }

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

inline su2double CEulerSolver::GetTotal_CpDiff() { return Total_CpDiff; }

inline su2double CEulerSolver::GetTotal_HeatFluxDiff() { return Total_HeatFluxDiff; }

inline su2double CEulerSolver::GetTotal_CNearFieldOF() { return Total_CNearFieldOF; }

inline void CEulerSolver::SetTotal_CEquivArea(su2double cequivarea) { Total_CEquivArea = cequivarea; }

inline void CEulerSolver::SetTotal_CpDiff(su2double pressure) { Total_CpDiff = pressure; }

inline void CEulerSolver::SetTotal_HeatFluxDiff(su2double heat) { Total_HeatFluxDiff = heat; }

inline void CEulerSolver::SetTotal_CNearFieldOF(su2double cnearfieldpress) { Total_CNearFieldOF = cnearfieldpress; }

inline void CEulerSolver::SetTotal_CLift(su2double val_Total_CLift) { Total_CLift = val_Total_CLift; }

inline void CEulerSolver::SetTotal_CDrag(su2double val_Total_CDrag) { Total_CDrag = val_Total_CDrag; }

inline su2double CEulerSolver::GetAllBound_CLift_Inv() { return AllBound_CLift_Inv; }

inline su2double CEulerSolver::GetAllBound_CDrag_Inv() { return AllBound_CDrag_Inv; }

inline su2double CEulerSolver::GetAllBound_CSideForce_Inv() { return AllBound_CSideForce_Inv; }

inline su2double CEulerSolver::GetAllBound_CEff_Inv() { return AllBound_CEff_Inv; }

inline su2double CEulerSolver::GetAllBound_CMx_Inv() { return AllBound_CMx_Inv; }

inline su2double CEulerSolver::GetAllBound_CMy_Inv() { return AllBound_CMy_Inv; }

inline su2double CEulerSolver::GetAllBound_CMz_Inv() { return AllBound_CMz_Inv; }

inline su2double CEulerSolver::GetAllBound_CFx_Inv() { return AllBound_CFx_Inv; }

inline su2double CEulerSolver::GetAllBound_CFy_Inv() { return AllBound_CFy_Inv; }

inline su2double CEulerSolver::GetAllBound_CFz_Inv() { return AllBound_CFz_Inv; }

inline su2double CEulerSolver::GetTotal_CFreeSurface() { return Total_CFreeSurface; }

inline void CEulerSolver::SetTotal_CFreeSurface(su2double cfreesurface) { Total_CFreeSurface = cfreesurface; }

inline su2double CEulerSolver::GetOneD_TotalPress(void) { return OneD_TotalPress; }

inline void CEulerSolver::SetOneD_TotalPress(su2double AveragePressure) { OneD_TotalPress = AveragePressure; }

inline su2double CEulerSolver::GetOneD_Mach(void) {return OneD_Mach;}

inline void CEulerSolver::SetOneD_Mach(su2double AverageMach) { OneD_Mach = AverageMach; }

inline su2double CEulerSolver::GetOneD_Temp(void) {return OneD_Temp;}

inline void CEulerSolver::SetOneD_Temp(su2double AverageTemperature) { OneD_Temp = AverageTemperature; }

inline su2double CEulerSolver::GetOneD_MassFlowRate(void) {return OneD_MassFlowRate;}

inline void CEulerSolver::SetOneD_MassFlowRate(su2double MassFlowRate) { OneD_MassFlowRate = MassFlowRate; }

inline su2double CEulerSolver::GetOneD_FluxAvgPress(void) {return OneD_PressureRef;}

inline void CEulerSolver::SetOneD_FluxAvgPress(su2double PressureRef) {OneD_PressureRef = PressureRef; }

inline su2double CEulerSolver::GetOneD_FluxAvgDensity(void) {return OneD_DensityRef;}

inline void CEulerSolver::SetOneD_FluxAvgDensity(su2double DensityRef) {OneD_DensityRef = DensityRef; }

inline su2double CEulerSolver::GetOneD_FluxAvgVelocity(void) {return OneD_VelocityRef;}

inline void CEulerSolver::SetOneD_FluxAvgVelocity(su2double VelocityRef) {OneD_VelocityRef = VelocityRef; }

inline su2double CEulerSolver::GetOneD_FluxAvgEntalpy(void) {return OneD_EnthalpyRef;}

inline void CEulerSolver::SetOneD_FluxAvgEntalpy(su2double EnthalpyRef) {OneD_EnthalpyRef = EnthalpyRef; }

inline su2double CEulerSolver::GetAveragedDensity(unsigned short valMarker){return AveragedDensity[valMarker];}

inline su2double CEulerSolver::GetAveragedPressure(unsigned short valMarker){return AveragedPressure[valMarker];}

inline su2double CEulerSolver::GetAveragedEnthalpy(unsigned short valMarker){return AveragedEnthalpy[valMarker];}

inline su2double CEulerSolver::GetAveragedEntropy(unsigned short valMarker){return AveragedEntropy[valMarker];}

inline su2double* CEulerSolver::GetAveragedVelocity(unsigned short valMarker){return AveragedVelocity[valMarker];}

inline su2double* CEulerSolver::GetAveragedGridVelocity(unsigned short valMarker){return AveragedGridVel[valMarker];}

inline su2double CEulerSolver::GetAveragedNormalVelocity(unsigned short valMarker){return AveragedNormalVelocity[valMarker];}

inline su2double CEulerSolver::GetAveragedTangVelocity(unsigned short valMarker){return AveragedTangVelocity[valMarker];}

inline su2double CEulerSolver::GetAveragedTotTemperature(unsigned short valMarker){return AveragedTotTemperature[valMarker];}

inline su2double CEulerSolver::GetAveragedTotPressure(unsigned short valMarker){return AveragedTotPressure[valMarker];}

inline su2double CEulerSolver::GetMassFlow(unsigned short valMarker){return MassFlow[valMarker];}

inline su2double CEulerSolver::GetFlowAngle(unsigned short valMarker){return FlowAngle[valMarker];}

inline su2double CEulerSolver::GetAveragedMach(unsigned short valMarker){return AveragedMach[valMarker];}

inline su2double CEulerSolver::GetAveragedNormalMach(unsigned short valMarker){return AveragedNormalMach[valMarker];}

inline su2double CEulerSolver::GetTotalPressureLoss(unsigned short inMarkerTP){return TotalPressureLoss[inMarkerTP];}

inline su2double CEulerSolver::GetKineticEnergyLoss(unsigned short inMarkerTP){return KineticEnergyLoss[inMarkerTP];}

inline su2double CEulerSolver::GetTotalStaticEfficiency(unsigned short inMarkerTP){return TotalStaticEfficiency[inMarkerTP];}

inline su2double CEulerSolver::GetTotalTotalEfficiency(unsigned short inMarkerTP){return TotalTotalEfficiency[inMarkerTP];}

inline su2double CEulerSolver::GetEulerianWork(unsigned short inMarkerTP){return EulerianWork[inMarkerTP];}

inline su2double CEulerSolver::GetTotalEnthalpyIn(unsigned short inMarkerTP){return TotalEnthalpyIn[inMarkerTP];}

inline su2double CEulerSolver::GetFlowAngleIn(unsigned short inMarkerTP){return FlowAngleIn[inMarkerTP];}

inline su2double CEulerSolver::GetFlowAngleOut(unsigned short inMarkerTP){return FlowAngleOut[inMarkerTP];}

inline su2double CEulerSolver::GetMassFlowIn(unsigned short inMarkerTP){return MassFlowIn[inMarkerTP];}

inline su2double CEulerSolver::GetMassFlowOut(unsigned short inMarkerTP){return MassFlowOut[inMarkerTP];}

inline su2double CEulerSolver::GetMachIn(unsigned short inMarkerTP){return MachIn[inMarkerTP];}

inline su2double CEulerSolver::GetMachOut(unsigned short inMarkerTP){return MachOut[inMarkerTP];}

inline su2double CEulerSolver::GetNormalMachIn(unsigned short inMarkerTP){return NormalMachIn[inMarkerTP];}

inline su2double CEulerSolver::GetNormalMachOut(unsigned short inMarkerTP){return NormalMachOut[inMarkerTP];}

inline su2double CEulerSolver::GetEnthalpyOut(unsigned short inMarkerTP){return EnthalpyOut[inMarkerTP];}

inline su2double CEulerSolver::GetVelocityOutIs(unsigned short inMarkerTP){return VelocityOutIs[inMarkerTP];}

inline su2double CEulerSolver::GetPressureOut(unsigned short inMarkerTP){return PressureOut[inMarkerTP];}

inline su2double CEulerSolver::GetPressureRatio(unsigned short inMarkerTP){return PressureRatio[inMarkerTP];}

inline void CEulerSolver::SetExtAveragedValue(CSolver *solver_container, unsigned short intMarker,  unsigned short extMarker){
	ExtAveragedDensity[extMarker]= solver_container->GetAveragedDensity(intMarker);
    ExtAveragedPressure[extMarker]= solver_container->GetAveragedPressure(intMarker);
    ExtAveragedNormalVelocity[extMarker]= solver_container->GetAveragedNormalVelocity(intMarker);
    ExtAveragedTangVelocity[extMarker]= solver_container->GetAveragedTangVelocity(intMarker);
    ExtAveragedTotTemperature[extMarker]= solver_container->GetAveragedTotTemperature(intMarker);
    ExtAveragedTotPressure[extMarker]= solver_container->GetAveragedTotPressure(intMarker);
}

inline void CEulerSolver::StoreTurboPerformance(CSolver *solver, unsigned short inMarkerTP ){
	TotalPressureLoss[inMarkerTP] = solver->GetTotalPressureLoss(inMarkerTP);
	KineticEnergyLoss[inMarkerTP] = solver->GetKineticEnergyLoss(inMarkerTP);
	TotalTotalEfficiency[inMarkerTP] = solver->GetTotalTotalEfficiency(inMarkerTP);
	TotalStaticEfficiency[inMarkerTP]= solver->GetTotalStaticEfficiency(inMarkerTP);
	EulerianWork[inMarkerTP] = solver->GetEulerianWork(inMarkerTP);
	TotalEnthalpyIn[inMarkerTP]= solver->GetTotalEnthalpyIn(inMarkerTP);
	FlowAngleIn[inMarkerTP]= solver->GetFlowAngleIn(inMarkerTP);
	FlowAngleOut[inMarkerTP]= solver->GetFlowAngleOut(inMarkerTP);
	MassFlowIn[inMarkerTP]= solver->GetMassFlowIn(inMarkerTP);
	MassFlowOut[inMarkerTP]= solver->GetMassFlowOut(inMarkerTP);
	MachIn[inMarkerTP]= solver->GetMachIn(inMarkerTP);
	MachOut[inMarkerTP]= solver->GetMachOut(inMarkerTP);
	NormalMachIn[inMarkerTP]= solver->GetNormalMachIn(inMarkerTP);
	NormalMachOut[inMarkerTP]= solver->GetNormalMachOut(inMarkerTP);
	EnthalpyOut[inMarkerTP]= solver->GetEnthalpyOut(inMarkerTP);
	VelocityOutIs[inMarkerTP]= solver->GetVelocityOutIs(inMarkerTP);
	PressureOut[inMarkerTP] = solver->GetPressureOut(inMarkerTP);
	PressureRatio[inMarkerTP] = solver->GetPressureRatio(inMarkerTP);

}

inline su2double CNSSolver::GetViscosity_Inf(void) { return Viscosity_Inf; }

inline su2double CNSSolver::GetTke_Inf(void) { return Tke_Inf; }

inline su2double CNSSolver::GetCLift_Visc(unsigned short val_marker) { return CLift_Visc[val_marker]; }

inline su2double CNSSolver::GetCMz_Visc(unsigned short val_marker) { return CMz_Visc[val_marker]; }

inline su2double CNSSolver::GetCSideForce_Visc(unsigned short val_marker) { return CSideForce_Visc[val_marker]; }

inline su2double CNSSolver::GetCDrag_Visc(unsigned short val_marker) { return CDrag_Visc[val_marker]; }

inline su2double CNSSolver::GetAllBound_CLift_Visc() { return AllBound_CLift_Visc; }

inline su2double CNSSolver::GetAllBound_CSideForce_Visc() { return AllBound_CSideForce_Visc; }

inline su2double CNSSolver::GetAllBound_CDrag_Visc() { return AllBound_CDrag_Visc; }

inline su2double CNSSolver::GetCSkinFriction(unsigned short val_marker, unsigned long val_vertex) { return CSkinFriction[val_marker][val_vertex]; }

inline su2double CNSSolver::GetHeatFlux(unsigned short val_marker, unsigned long val_vertex) { return HeatFlux[val_marker][val_vertex]; }

inline su2double CNSSolver::GetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex) { return HeatFluxTarget[val_marker][val_vertex]; }

inline void CNSSolver::SetHeatFluxTarget(unsigned short val_marker, unsigned long val_vertex, su2double val_heat) { HeatFluxTarget[val_marker][val_vertex] = val_heat; }

inline su2double CNSSolver::GetYPlus(unsigned short val_marker, unsigned long val_vertex) { return YPlus[val_marker][val_vertex]; }

inline su2double CNSSolver::GetStrainMag_Max(void) { return StrainMag_Max; }

inline su2double CNSSolver::GetOmega_Max(void) { return Omega_Max; }

inline void CNSSolver::SetStrainMag_Max(su2double val_strainmag_max) { StrainMag_Max = val_strainmag_max; }

inline void CNSSolver::SetOmega_Max(su2double val_omega_max) { Omega_Max = val_omega_max; }

inline su2double CAdjEulerSolver::GetCSensitivity(unsigned short val_marker, unsigned long val_vertex) { return CSensitivity[val_marker][val_vertex]; }

inline void CAdjEulerSolver::SetCSensitivity(unsigned short val_marker, unsigned long val_vertex, su2double val_sensitivity) {CSensitivity[val_marker][val_vertex]=val_sensitivity; }

inline su2double CAdjEulerSolver::GetTotal_Sens_Geo() { return Total_Sens_Geo; }

inline su2double CAdjEulerSolver::GetTotal_Sens_Mach() { return Total_Sens_Mach; }

inline su2double CAdjEulerSolver::GetTotal_Sens_AoA() { return Total_Sens_AoA; }

inline su2double CAdjEulerSolver::GetTotal_Sens_Press() { return Total_Sens_Press; }

inline su2double CAdjEulerSolver::GetTotal_Sens_Temp() { return Total_Sens_Temp; }

inline su2double CAdjEulerSolver::GetTotal_Sens_BPress() { return Total_Sens_BPress; }

inline su2double CAdjEulerSolver::GetPsiRho_Inf(void) { return PsiRho_Inf; }

inline su2double CAdjEulerSolver::GetPsiE_Inf(void) { return PsiE_Inf; }

inline su2double CAdjEulerSolver::GetPhi_Inf(unsigned short val_dim) { return Phi_Inf[val_dim]; }

inline su2double CFEASolver::GetTotal_CFEA() { return Total_CFEA; }

inline void CFEASolver::SetTotal_CFEA(su2double cfea) { Total_CFEA = cfea; }

inline su2double CFEASolver::GetWAitken_Dyn(void) { return WAitken_Dyn; }

inline su2double CFEASolver::GetWAitken_Dyn_tn1(void) { return WAitken_Dyn_tn1; }

inline void CFEASolver::SetWAitken_Dyn(su2double waitk) { WAitken_Dyn = waitk; }

inline void CFEASolver::SetWAitken_Dyn_tn1(su2double waitk_tn1) { WAitken_Dyn_tn1 = waitk_tn1; }

inline void CFEASolver::SetFSI_ConvValue(unsigned short val_index, su2double val_criteria) { FSI_Conv[val_index] = val_criteria; }

inline su2double CFEASolver::GetFSI_ConvValue(unsigned short val_index){ return FSI_Conv[val_index]; }

inline su2double CWaveSolver::GetTotal_CWave() { return Total_CWave; }

inline su2double CHeatSolver::GetTotal_CHeat() { return Total_CHeat; }

inline void CSolver::RegisterSolution(CGeometry *geometry_container, CConfig *config){}

inline void CSolver::RegisterOutput(CGeometry *geometry_container, CConfig *config){}

inline void CSolver::SetAdjoint_Output(CGeometry *geometry, CConfig *config){}

inline void CSolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config){}

inline void CSolver::RegisterObj_Func(CConfig *config){}

inline void CSolver::SetSurface_Sensitivity(CGeometry *geometry, CConfig *config){}

inline void CSolver::SetSensitivity(CGeometry *geometry, CConfig *config){}

inline void CSolver::SetAdj_ObjFunc(CGeometry *geometry, CConfig *config){}

inline su2double CDiscAdjSolver::GetTotal_Sens_Geo() { return Total_Sens_Geo; }

inline su2double CDiscAdjSolver::GetTotal_Sens_Mach() { return Total_Sens_Mach; }

inline su2double CDiscAdjSolver::GetTotal_Sens_AoA() { return Total_Sens_AoA; }

inline su2double CDiscAdjSolver::GetTotal_Sens_Press() { return Total_Sens_Press; }

inline su2double CDiscAdjSolver::GetTotal_Sens_Temp() { return Total_Sens_Temp; }

inline su2double CDiscAdjSolver::GetTotal_Sens_BPress() { return Total_Sens_BPress; }

inline su2double CDiscAdjSolver::GetCSensitivity(unsigned short val_marker, unsigned long val_vertex) { return CSensitivity[val_marker][val_vertex]; }

inline unsigned long CSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {return 0;}

inline void CSolver::SetRecording(CGeometry *geometry, CConfig *config, unsigned short kind_recording){}

inline void CSolver::SetPressure_Inf(su2double p_inf){}

inline void CSolver::SetTemperature_Inf(su2double t_inf){}

inline void CEulerSolver::SetPressure_Inf(su2double p_inf){Pressure_Inf = p_inf;}

inline void CEulerSolver::SetTemperature_Inf(su2double t_inf){Temperature_Inf = t_inf;}

inline void CSolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset){}

inline void CSolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config){}

