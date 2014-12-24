/*!
 * \file solver_structure.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 3.2.6 "eagle"
 *
 * Copyright (C) 2012-2014 SU2 <https://github.com/su2code>.
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

inline double CSolver::GetCSensitivity(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline void CSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep, 
																		 unsigned short iMesh, unsigned short RunTime_EqSystem) { }
																		 
inline void CSolver::SetFreeSurface_Distance(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry, CGeometry **flow_geometry, CConfig *fea_config, CConfig *flow_config) { }

inline void CSolver::GetSurface_Pressure(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) { }

inline void CSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter) { }
  
inline void CSolver::SetFlow_Displacement(CGeometry **flow_geometry, CVolumetricMovement *flow_grid_movement, CConfig *flow_config, CConfig *fea_config, CGeometry **fea_geometry, CSolver ***fea_solution) { }

inline void CSolver::SetCSensitivity(unsigned short val_marker, unsigned short val_vertex, double val_sensitivity) { }

inline void CSolver::Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline void CSolver::Smooth_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline void CSolver::Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline double CSolver::GetPhi_Inf(unsigned short val_dim) { return 0; }

inline double CSolver::GetPsiRho_Inf(void) { return 0; }

inline double* CSolver::GetPsiRhos_Inf(void) { return NULL; }

inline double CSolver::GetPsiE_Inf(void) { return 0; }

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

inline void CSolver::SetPreconditioner(CConfig *config, unsigned short iPoint) { }

inline void CSolver::SetDistance(CGeometry *geometry, CConfig *config) { };

inline double CSolver::GetCLift_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetCMz_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetCDrag_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CLift(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CDrag(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CSideForce(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CEff(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CFx(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CFy(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CFz(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CMx(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CMy(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CMz(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CLift_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CDrag_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CSideForce_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CEff_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CFx_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CFy_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CFz_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CMx_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CMy_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetSurface_CMz_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetInflow_MassFlow(unsigned short val_marker) { return 0; }

inline double CSolver::GetExhaust_MassFlow(unsigned short val_marker) { return 0; }

inline double CSolver::GetInflow_Pressure(unsigned short val_marker) { return 0; }

inline double CSolver::GetInflow_Mach(unsigned short val_marker) { return 0; }

inline double CSolver::GetCSideForce_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetCEff_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetCLift_Visc(unsigned short val_marker) { return 0; }

inline double CSolver::GetCMz_Visc(unsigned short val_marker) { return 0; }

inline double CSolver::GetCSideForce_Visc(unsigned short val_marker) { return 0; }

inline double CSolver::GetCDrag_Visc(unsigned short val_marker) { return 0; }

inline double CSolver::GetAllBound_CLift_Inv() { return 0; }

inline double CSolver::GetAllBound_CDrag_Inv() { return 0; }

inline double CSolver::GetAllBound_CSideForce_Inv() { return 0; }

inline double CSolver::GetAllBound_CEff_Inv() { return 0; }

inline double CSolver::GetAllBound_CMx_Inv() { return 0; }

inline double CSolver::GetAllBound_CMy_Inv() { return 0; }

inline double CSolver::GetAllBound_CMz_Inv() { return 0; }

inline double CSolver::GetAllBound_CFx_Inv() { return 0; }

inline double CSolver::GetAllBound_CFy_Inv() { return 0; }

inline double CSolver::GetAllBound_CFz_Inv() { return 0; }

inline double CSolver::GetAllBound_CLift_Visc() { return 0; }

inline double CSolver::GetAllBound_CSideForce_Visc() { return 0; }

inline double CSolver::GetAllBound_CDrag_Visc() { return 0; }

inline void CSolver::SetForceProj_Vector(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::SetIntBoundary_Jump(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline double CSolver::GetTotal_CLift() { return 0; }

inline double CSolver::GetTotal_CDrag() { return 0; }

inline double CSolver::GetTotal_CMx() { return 0; }

inline double CSolver::GetTotal_CMy() { return 0; }

inline double CSolver::GetTotal_CMz() { return 0; }

inline double CSolver::GetTotal_CFx() { return 0; }

inline double CSolver::GetTotal_CFy() { return 0; }

inline double CSolver::GetTotal_CFz() { return 0; }

inline double CSolver::GetTotal_CSideForce() { return 0; }

inline double CSolver::GetTotal_CEff() { return 0; }

inline double CSolver::GetTotal_CT() { return 0; }

inline void CSolver::SetTotal_CT(double val_Total_CT) { }

inline double CSolver::GetTotal_CQ() { return 0; }

inline double CSolver::GetTotal_HeatFlux() { return 0; }

inline double CSolver::GetTotal_MaxHeatFlux() { return 0; }

inline double CSolver::Get_PressureDrag() { return 0; }

inline double CSolver::Get_ViscDrag() { return 0; }

inline void CSolver::SetTotal_CQ(double val_Total_CQ) { }

inline void CSolver::SetTotal_HeatFlux(double val_Total_Heat) { }

inline void CSolver::SetTotal_MaxHeatFlux(double val_Total_Heat) { }

inline double CSolver::GetTotal_CMerit() { return 0; }

inline double CSolver::GetTotal_CEquivArea() { return 0; }

inline double CSolver::GetTotal_CpDiff() { return 0; }

inline double CSolver::GetTotal_HeatFluxDiff() { return 0; }

inline double CSolver::GetTotal_CFreeSurface() { return 0; }

inline double CSolver::GetTotal_CFEA() { return 0; }

inline double CSolver::GetTotal_CNearFieldOF() { return 0; }

inline void CSolver::SetTotal_CEquivArea(double val_cequivarea) { }

inline void CSolver::SetTotal_CpDiff(double val_pressure) { }

inline void CSolver::SetTotal_HeatFluxDiff(double val_heat) { }

inline void CSolver::SetTotal_CFEA(double val_cfea) { }

inline void CSolver::SetTotal_CFreeSurface(double val_freesurface) { }

inline void CSolver::SetTotal_CNearFieldOF(double val_cnearfieldpress) { }

inline double CSolver::GetTotal_CDeltaLift() { return 0; }

inline double CSolver::GetTotal_CDeltaDrag() { return 0; }

inline double CSolver::GetTotal_CWave() { return 0; }

inline double CSolver::GetTotal_CHeat() { return 0; }

inline void CSolver::SetTotal_CLift(double val_Total_CLift) { }

inline void CSolver::SetTotal_CDrag(double val_Total_CDrag) { }

inline double CSolver::GetCPressure(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolver::GetCPressureTarget(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline void CSolver::SetCPressureTarget(unsigned short val_marker, unsigned short val_vertex, double val_pressure) { }

inline void CSolver::SetHeatFluxTarget(unsigned short val_marker, unsigned short val_vertex, double val_heat) { }

inline double *CSolver::GetCharacPrimVar(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolver::GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolver::GetHeatFlux(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolver::GetHeatFluxTarget(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolver::GetYPlus(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline void CSolver::Viscous_Residual(CGeometry *geometry,
                                      CSolver **solver_container,
                                      CNumerics *numerics, CConfig
                                      *config, unsigned short iMesh,
                                      unsigned short iRKstep) { }
									   
inline void CSolver::AddStiffMatrix(double ** StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1, unsigned long Point_2, unsigned long Point_3) { }
									   
inline void CSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, 
												  CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) { }
									   
inline void CSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, 
												  CNumerics *numerics, CConfig *config, unsigned short iMesh) { }

inline double CSolver::GetTotal_Sens_Geo() { return 0; }

inline double CSolver::GetTotal_Sens_Mach() { return 0; }

inline double CSolver::GetTotal_Sens_AoA() { return 0; }

inline double CSolver::GetTotal_Sens_Press() { return 0; }

inline double CSolver::GetTotal_Sens_Temp() { return 0; }

inline double CSolver::GetDensity_Inf(void) { return 0; }

inline double CSolver::GetDensity_Inf(unsigned short val_var) { return 0; }

inline double CSolver::GetModVelocity_Inf(void) { return 0; }

inline double CSolver::GetDensity_Energy_Inf(void) { return 0; }

inline double CSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return 0; }

inline double CSolver::GetDensity_Velocity_Inf(unsigned short val_dim, unsigned short val_var) { return 0; }

inline double CSolver::GetVelocity_Inf(unsigned short val_dim) { return 0; }

inline double* CSolver::GetVelocity_Inf(void) { return 0; }

inline double CSolver::GetPressure_Inf(void) { return 0; }

inline double CSolver::GetViscosity_Inf(void) { return 0; }

inline double CSolver::GetTke_Inf(void) { return 0; }

inline double* CSolver::GetConstants() {return NULL;}

inline double CSolver::GetOneD_TotalPress(void){return 0;}

inline void CSolver::SetOneD_TotalPress(double AveragePressure){ }

inline double CSolver::GetOneD_Mach(void){return 0;}

inline void CSolver::SetOneD_Mach(double AverageMach){ }

inline double CSolver::GetOneD_Temp(void){return 0;}

inline void CSolver::SetOneD_Temp(double AverageTemperature){ }

inline double CSolver::GetOneD_MassFlowRate(void){return 0;}

inline void CSolver::SetOneD_MassFlowRate(double MassFlowRate){ }

inline double CSolver::GetOneD_FluxAvgPress(void){return 0;}

inline void CSolver::SetOneD_FluxAvgPress(double PressureRef){ }

inline double CSolver::GetOneD_FluxAvgDensity(void){return 0;}

inline void CSolver::SetOneD_FluxAvgDensity(double DensityRef){ }

inline double CSolver::GetOneD_FluxAvgVelocity(void){return 0;}

inline void CSolver::SetOneD_FluxAvgVelocity(double VelocityRef){ }

inline double CSolver::GetOneD_FluxAvgEntalpy(void){return 0;}

inline void CSolver::SetOneD_FluxAvgEntalpy(double EnthalpyRef){ }

inline void CSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }
									 
inline void CSolver::BC_Normal_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }
									 									 
inline void CSolver::BC_Flow_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }
									 
inline void CSolver::BC_Normal_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }
                   
inline void CSolver::BC_Pressure(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
									 unsigned short val_marker) { }
                  
inline void CSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_IsothermalCatalytic_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_IsothermalNonCatalytic_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }
                  
inline void CSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_HeatFluxCatalytic_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_HeatFluxNonCatalytic_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }
									
inline void CSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
								  unsigned short val_marker) { }

inline void CSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									CConfig *config, unsigned short val_marker) { }
                  
inline void CSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									CConfig *config, unsigned short val_marker) { }

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
										 
inline void CSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
										 CConfig *config, unsigned short val_marker) { }
                     
inline void CSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
										  CConfig *config, unsigned short val_marker) { }
                      
                      
inline void CSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
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
            
inline void CSolver::GetEngine_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) { }

inline void CSolver::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                     CConfig *config, unsigned short iMesh, bool Output) { }

inline void CSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
							        unsigned short iMesh, unsigned long Iteration) { }	
							        
inline void CSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
							        unsigned short iMesh) { }								        

inline void CSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
										CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

inline void CSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									   CConfig *config, unsigned short iMesh) { }

inline void CSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) { }

inline void CSolver::SetDissipation_Switch(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Dissipation_Switch(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Viscous_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Inviscid_DeltaForces(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Viscous_DeltaForces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Wave_Strength(CGeometry *geometry, CConfig *config) { }

inline void CSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, 
											CConfig *config, unsigned short iRKStep) { }

inline void CSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Compute_Residual(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
										unsigned short iMesh) { }

inline void CSolver::SetRes_RMS(unsigned short val_var, double val_residual) { Residual_RMS[val_var] = val_residual; }

inline void CSolver::AddRes_RMS(unsigned short val_var, double val_residual) { Residual_RMS[val_var] += val_residual; }

inline double CSolver::GetRes_RMS(unsigned short val_var) { return Residual_RMS[val_var]; }

inline void CSolver::SetRes_Max(unsigned short val_var, double val_residual, unsigned long val_point) { Residual_Max[val_var] = val_residual; Point_Max[val_var] = val_point; }

inline void CSolver::AddRes_Max(unsigned short val_var, double val_residual, unsigned long val_point, double* val_coord) {
  if (val_residual > Residual_Max[val_var]) {
  Residual_Max[val_var] = val_residual;
  Point_Max[val_var] = val_point;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Point_Max_Coord[val_var][iDim] = val_coord[iDim];
  }
}

inline double CSolver::GetRes_Max(unsigned short val_var) { return Residual_Max[val_var]; }

inline unsigned long CSolver::GetPoint_Max(unsigned short val_var) { return Point_Max[val_var]; }

inline double* CSolver::GetPoint_Max_Coord(unsigned short val_var) { return Point_Max_Coord[val_var]; }

inline void CSolver::Set_OldSolution(CGeometry *geometry) {
	for(unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) 
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

inline double CSolver::GetMax_Delta_Time(void) { return Max_Delta_Time; }

inline double CSolver::GetMin_Delta_Time(void) { return Min_Delta_Time; }

inline double CSolver::GetMax_Delta_Time(unsigned short val_Species) { return 0.0; }

inline double CSolver::GetMin_Delta_Time(unsigned short val_Species) { return 0.0; }

inline void CSolver::Copy_Zone_Solution(CSolver ***solver1_solution, CGeometry **solver1_geometry, CConfig *solver1_config, 
										  CSolver ***solver2_solution, CGeometry **solver2_geometry, CConfig *solver2_config) {};

inline CFluidModel* CSolver::GetFluidModel(void) { return NULL;}

inline CFluidModel* CEulerSolver::GetFluidModel(void) { return FluidModel;}

inline double CEulerSolver::GetDensity_Inf(void) { return Density_Inf; }

inline double CEulerSolver::GetModVelocity_Inf(void) { 
	double Vel2 = 0; 
	for (unsigned short iDim = 0; iDim < nDim; iDim++) 
		Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim]; 
	return sqrt(Vel2);
}

inline double CEulerSolver::GetDensity_Energy_Inf(void) { return Density_Inf*Energy_Inf; }

inline double CEulerSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return Density_Inf*Velocity_Inf[val_dim]; }

inline double CEulerSolver::GetVelocity_Inf(unsigned short val_dim) { return Velocity_Inf[val_dim]; }

inline double *CEulerSolver::GetVelocity_Inf(void) { return Velocity_Inf; }

inline double CEulerSolver::GetPressure_Inf(void) { return Pressure_Inf; }

inline double CEulerSolver::GetCPressure(unsigned short val_marker, unsigned short val_vertex) { return CPressure[val_marker][val_vertex]; }

inline double CEulerSolver::GetCPressureTarget(unsigned short val_marker, unsigned short val_vertex) { return CPressureTarget[val_marker][val_vertex]; }

inline void CEulerSolver::SetCPressureTarget(unsigned short val_marker, unsigned short val_vertex, double val_pressure) { CPressureTarget[val_marker][val_vertex] = val_pressure; }

inline double *CEulerSolver::GetCharacPrimVar(unsigned short val_marker, unsigned short val_vertex) { return CharacPrimVar[val_marker][val_vertex]; }

inline double CEulerSolver::GetCLift_Inv(unsigned short val_marker) { return CLift_Inv[val_marker]; }

inline double CEulerSolver::GetCMz_Inv(unsigned short val_marker) { return CMz_Inv[val_marker]; }

inline double CEulerSolver::GetCDrag_Inv(unsigned short val_marker) { return CDrag_Inv[val_marker]; }

inline double CEulerSolver::GetSurface_CLift(unsigned short val_marker) { return Surface_CLift[val_marker]; }

inline double CEulerSolver::GetSurface_CDrag(unsigned short val_marker) { return Surface_CDrag[val_marker]; }

inline double CEulerSolver::GetSurface_CSideForce(unsigned short val_marker) { return Surface_CSideForce[val_marker]; }

inline double CEulerSolver::GetSurface_CEff(unsigned short val_marker) { return Surface_CEff[val_marker]; }

inline double CEulerSolver::GetSurface_CFx(unsigned short val_marker) { return Surface_CFx[val_marker]; }

inline double CEulerSolver::GetSurface_CFy(unsigned short val_marker) { return Surface_CFy[val_marker]; }

inline double CEulerSolver::GetSurface_CFz(unsigned short val_marker) { return Surface_CFz[val_marker]; }

inline double CEulerSolver::GetSurface_CMx(unsigned short val_marker) { return Surface_CMx[val_marker]; }

inline double CEulerSolver::GetSurface_CMy(unsigned short val_marker) { return Surface_CMy[val_marker]; }

inline double CEulerSolver::GetSurface_CMz(unsigned short val_marker) { return Surface_CMz[val_marker]; }

inline double CEulerSolver::GetSurface_CLift_Inv(unsigned short val_marker) { return Surface_CLift_Inv[val_marker]; }

inline double CEulerSolver::GetSurface_CDrag_Inv(unsigned short val_marker) { return Surface_CDrag_Inv[val_marker]; }

inline double CEulerSolver::GetSurface_CSideForce_Inv(unsigned short val_marker) { return Surface_CSideForce_Inv[val_marker]; }

inline double CEulerSolver::GetSurface_CEff_Inv(unsigned short val_marker) { return Surface_CEff_Inv[val_marker]; }

inline double CEulerSolver::GetSurface_CFx_Inv(unsigned short val_marker) { return Surface_CFx_Inv[val_marker]; }

inline double CEulerSolver::GetSurface_CFy_Inv(unsigned short val_marker) { return Surface_CFy_Inv[val_marker]; }

inline double CEulerSolver::GetSurface_CFz_Inv(unsigned short val_marker) { return Surface_CFz_Inv[val_marker]; }

inline double CEulerSolver::GetSurface_CMx_Inv(unsigned short val_marker) { return Surface_CMx_Inv[val_marker]; }

inline double CEulerSolver::GetSurface_CMy_Inv(unsigned short val_marker) { return Surface_CMy_Inv[val_marker]; }

inline double CEulerSolver::GetSurface_CMz_Inv(unsigned short val_marker) { return Surface_CMz_Inv[val_marker]; }

inline double CEulerSolver::GetInflow_MassFlow(unsigned short val_marker) { return Inflow_MassFlow[val_marker]; }

inline double CEulerSolver::GetExhaust_MassFlow(unsigned short val_marker) { return Exhaust_MassFlow[val_marker]; }

inline double CEulerSolver::GetInflow_Pressure(unsigned short val_marker) { return Inflow_Pressure[val_marker]; }

inline double CEulerSolver::GetInflow_Mach(unsigned short val_marker) { return Inflow_Mach[val_marker]; }

inline double CEulerSolver::GetCSideForce_Inv(unsigned short val_marker) { return CSideForce_Inv[val_marker]; }

inline double CEulerSolver::GetCEff_Inv(unsigned short val_marker) { return CEff_Inv[val_marker]; }

inline double CEulerSolver::GetTotal_CLift() { return Total_CLift; }

inline double CEulerSolver::GetTotal_CDrag() { return Total_CDrag; }

inline double CEulerSolver::GetTotal_CMx() { return Total_CMx; }

inline double CEulerSolver::GetTotal_CMy() { return Total_CMy; }

inline double CEulerSolver::GetTotal_CMz() { return Total_CMz; }

inline double CEulerSolver::GetTotal_CFx() { return Total_CFx; }

inline double CEulerSolver::GetTotal_CFy() { return Total_CFy; }

inline double CEulerSolver::GetTotal_CFz() { return Total_CFz; }

inline double CEulerSolver::GetTotal_CSideForce() { return Total_CSideForce; }

inline double CEulerSolver::GetTotal_CEff() { return Total_CEff; }

inline double CEulerSolver::GetTotal_CT() { return Total_CT; }

inline void CEulerSolver::SetTotal_CT(double val_Total_CT) { Total_CT = val_Total_CT; }

inline double CEulerSolver::GetTotal_CQ() { return Total_CQ; }

inline double CEulerSolver::GetTotal_HeatFlux() { return Total_Heat; }

inline double CEulerSolver::GetTotal_MaxHeatFlux() { return Total_MaxHeat; }

inline void CEulerSolver::SetTotal_CQ(double val_Total_CQ) { Total_CQ = val_Total_CQ; }

inline void CEulerSolver::SetTotal_HeatFlux(double val_Total_Heat) { Total_Heat = val_Total_Heat; }

inline void CEulerSolver::SetTotal_MaxHeatFlux(double val_Total_MaxHeat) { Total_MaxHeat = val_Total_MaxHeat; }

inline double CEulerSolver::GetTotal_CMerit() { return Total_CMerit; }

inline double CEulerSolver::GetTotal_CEquivArea() { return Total_CEquivArea; }

inline double CEulerSolver::GetTotal_CpDiff() { return Total_CpDiff; }

inline double CEulerSolver::GetTotal_HeatFluxDiff() { return Total_HeatFluxDiff; }

inline double CEulerSolver::GetTotal_CNearFieldOF() { return Total_CNearFieldOF; }

inline void CEulerSolver::SetTotal_CEquivArea(double cequivarea) { Total_CEquivArea = cequivarea; }

inline void CEulerSolver::SetTotal_CpDiff(double pressure) { Total_CpDiff = pressure; }

inline void CEulerSolver::SetTotal_HeatFluxDiff(double heat) { Total_HeatFluxDiff = heat; }

inline void CEulerSolver::SetTotal_CNearFieldOF(double cnearfieldpress) { Total_CNearFieldOF = cnearfieldpress; }

inline void CEulerSolver::SetTotal_CLift(double val_Total_CLift) { Total_CLift = val_Total_CLift; }

inline void CEulerSolver::SetTotal_CDrag(double val_Total_CDrag) { Total_CDrag = val_Total_CDrag; }

inline double CEulerSolver::GetAllBound_CLift_Inv() { return AllBound_CLift_Inv; }

inline double CEulerSolver::GetAllBound_CDrag_Inv() { return AllBound_CDrag_Inv; }

inline double CEulerSolver::GetAllBound_CSideForce_Inv() { return AllBound_CSideForce_Inv; }

inline double CEulerSolver::GetAllBound_CEff_Inv() { return AllBound_CEff_Inv; }

inline double CEulerSolver::GetAllBound_CMx_Inv() { return AllBound_CMx_Inv; }

inline double CEulerSolver::GetAllBound_CMy_Inv() { return AllBound_CMy_Inv; }

inline double CEulerSolver::GetAllBound_CMz_Inv() { return AllBound_CMz_Inv; }

inline double CEulerSolver::GetAllBound_CFx_Inv() { return AllBound_CFx_Inv; }

inline double CEulerSolver::GetAllBound_CFy_Inv() { return AllBound_CFy_Inv; }

inline double CEulerSolver::GetAllBound_CFz_Inv() { return AllBound_CFz_Inv; }

inline double CEulerSolver::GetTotal_CFreeSurface() { return Total_CFreeSurface; }

inline void CEulerSolver::SetTotal_CFreeSurface(double cfreesurface) { Total_CFreeSurface = cfreesurface; }

inline double CEulerSolver::GetOneD_TotalPress(void) { return OneD_TotalPress; }

inline void CEulerSolver::SetOneD_TotalPress(double AveragePressure) { OneD_TotalPress = AveragePressure; }

inline double CEulerSolver::GetOneD_Mach(void){return OneD_Mach;}

inline void CEulerSolver::SetOneD_Mach(double AverageMach) { OneD_Mach = AverageMach; }

inline double CEulerSolver::GetOneD_Temp(void){return OneD_Temp;}

inline void CEulerSolver::SetOneD_Temp(double AverageTemperature) { OneD_Temp = AverageTemperature; }

inline double CEulerSolver::GetOneD_MassFlowRate(void){return OneD_MassFlowRate;}

inline void CEulerSolver::SetOneD_MassFlowRate(double MassFlowRate) { OneD_MassFlowRate = MassFlowRate; }

inline double CEulerSolver::GetOneD_FluxAvgPress(void){return OneD_PressureRef;}

inline void CEulerSolver::SetOneD_FluxAvgPress(double PressureRef){OneD_PressureRef = PressureRef; }

inline double CEulerSolver::GetOneD_FluxAvgDensity(void){return OneD_DensityRef;}

inline void CEulerSolver::SetOneD_FluxAvgDensity(double DensityRef){OneD_DensityRef = DensityRef; }

inline double CEulerSolver::GetOneD_FluxAvgVelocity(void){return OneD_VelocityRef;}

inline void CEulerSolver::SetOneD_FluxAvgVelocity(double VelocityRef){OneD_VelocityRef = VelocityRef; }

inline double CEulerSolver::GetOneD_FluxAvgEntalpy(void){return OneD_EnthalpyRef;}

inline void CEulerSolver::SetOneD_FluxAvgEntalpy(double EnthalpyRef){OneD_EnthalpyRef = EnthalpyRef; }

inline double CNSSolver::GetViscosity_Inf(void) { return Viscosity_Inf; }

inline double CNSSolver::GetTke_Inf(void) { return Tke_Inf; }

inline double CNSSolver::GetCLift_Visc(unsigned short val_marker) { return CLift_Visc[val_marker]; }

inline double CNSSolver::GetCMz_Visc(unsigned short val_marker) { return CMz_Visc[val_marker]; }

inline double CNSSolver::GetCSideForce_Visc(unsigned short val_marker) { return CSideForce_Visc[val_marker]; }

inline double CNSSolver::GetCDrag_Visc(unsigned short val_marker) { return CDrag_Visc[val_marker]; }

inline double CNSSolver::GetAllBound_CLift_Visc() { return AllBound_CLift_Visc; }

inline double CNSSolver::GetAllBound_CSideForce_Visc() { return AllBound_CSideForce_Visc; }

inline double CNSSolver::GetAllBound_CDrag_Visc() { return AllBound_CDrag_Visc; }

inline double CNSSolver::GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex) { return CSkinFriction[val_marker][val_vertex]; }

inline double CNSSolver::GetHeatFlux(unsigned short val_marker, unsigned short val_vertex) { return HeatFlux[val_marker][val_vertex]; }

inline double CNSSolver::GetHeatFluxTarget(unsigned short val_marker, unsigned short val_vertex) { return HeatFluxTarget[val_marker][val_vertex]; }

inline void CNSSolver::SetHeatFluxTarget(unsigned short val_marker, unsigned short val_vertex, double val_heat) { HeatFluxTarget[val_marker][val_vertex] = val_heat; }

inline double CNSSolver::GetYPlus(unsigned short val_marker, unsigned short val_vertex) { return YPlus[val_marker][val_vertex]; }

inline double CAdjEulerSolver::GetCSensitivity(unsigned short val_marker, unsigned short val_vertex) { return CSensitivity[val_marker][val_vertex]; }

inline void CAdjEulerSolver::SetCSensitivity(unsigned short val_marker, unsigned short val_vertex, double val_sensitivity) {CSensitivity[val_marker][val_vertex]=val_sensitivity; }

inline double CAdjEulerSolver::GetTotal_Sens_Geo() { return Total_Sens_Geo; }

inline double CAdjEulerSolver::GetTotal_Sens_Mach() { return Total_Sens_Mach; }

inline double CAdjEulerSolver::GetTotal_Sens_AoA() { return Total_Sens_AoA; }

inline double CAdjEulerSolver::GetTotal_Sens_Press() { return Total_Sens_Press; }

inline double CAdjEulerSolver::GetTotal_Sens_Temp() { return Total_Sens_Temp; }

inline double CAdjEulerSolver::GetPsiRho_Inf(void) { return PsiRho_Inf; }

inline double CAdjEulerSolver::GetPsiE_Inf(void) { return PsiE_Inf; }

inline double CAdjEulerSolver::GetPhi_Inf(unsigned short val_dim) { return Phi_Inf[val_dim]; }

inline double CLinEulerSolver::GetTotal_CDeltaLift() { return Total_CDeltaLift; }

inline double CLinEulerSolver::GetTotal_CDeltaDrag() { return Total_CDeltaDrag; }

inline double CFEASolver::GetTotal_CFEA() { return Total_CFEA; }

inline void CFEASolver::SetTotal_CFEA(double cfea) { Total_CFEA = cfea; }

inline double CWaveSolver::GetTotal_CWave() { return Total_CWave; }

inline double CHeatSolver::GetTotal_CHeat() { return Total_CHeat; }
      
inline double CTNE2EulerSolver::GetDensity_Inf(void) { cout << "CTNE2EulerSolver::GetDensity_Inf NOT RETURNING THE CORRECT VALUE!!!" << endl; return 0.0; }

inline double CTNE2EulerSolver::GetVelocity_Inf(unsigned short val_dim) { return Velocity_Inf[val_dim]; }

inline double CTNE2EulerSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { cout << "CTNE2EulerSolver::GetDensity_Velocity_Inf NOT RETURNING THE CORRECT VALUE!!!" << endl; return 0.0; }

inline double CTNE2EulerSolver::GetPressure_Inf(void) { return Pressure_Inf; }

inline double CTNE2EulerSolver::GetDensity_Energy_Inf(void) { cout << "CTNE2EulerSolver::GetDensity_Energy_Inf NOT RETURNING THE CORRECT VALUE!!!" << endl; return 0.0; }

inline double CTNE2EulerSolver::GetModVelocity_Inf(void) { 
	double Vel2 = 0; 
	for (unsigned short iDim = 0; iDim < nDim; iDim++) 
		Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim]; 
	return sqrt(Vel2);
}

inline double CTNE2EulerSolver::GetCPressure(unsigned short val_marker, unsigned short val_vertex) { return CPressure[val_marker][val_vertex]; }

inline double CTNE2EulerSolver::GetCLift_Inv(unsigned short val_marker) { return CLift_Inv[val_marker]; }

inline double CTNE2EulerSolver::GetCDrag_Inv(unsigned short val_marker) { return CDrag_Inv[val_marker]; }

inline double CTNE2EulerSolver::GetCSideForce_Inv(unsigned short val_marker) { return CSideForce_Inv[val_marker]; }

inline double CTNE2EulerSolver::GetCEff_Inv(unsigned short val_marker) { return CEff_Inv[val_marker]; }

inline double CTNE2EulerSolver::GetTotal_CLift() { return Total_CLift; }

inline double CTNE2EulerSolver::GetTotal_CDrag() { return Total_CDrag; }

inline double CTNE2EulerSolver::GetTotal_CMx() { return Total_CMx; }

inline double CTNE2EulerSolver::GetTotal_CMy() { return Total_CMy; }

inline double CTNE2EulerSolver::GetTotal_CMz() { return Total_CMz; }

inline double CTNE2EulerSolver::GetTotal_CFx() { return Total_CFx; }

inline double CTNE2EulerSolver::GetTotal_CFy() { return Total_CFy; }

inline double CTNE2EulerSolver::GetTotal_CFz() { return Total_CFz; }

inline double CTNE2EulerSolver::GetTotal_CSideForce() { return Total_CSideForce; }

inline double CTNE2EulerSolver::GetTotal_CEff() { return Total_CEff; }

inline double CTNE2EulerSolver::GetTotal_HeatFlux() { return Total_Heat; }

inline double CTNE2EulerSolver::GetTotal_MaxHeatFlux() { return Total_MaxHeat; }

inline void CTNE2EulerSolver::SetTotal_HeatFlux(double val_Total_Heat) { Total_Heat = val_Total_Heat; }

inline void CTNE2EulerSolver::SetTotal_MaxHeatFlux(double val_Total_MaxHeat) { Total_MaxHeat = val_Total_MaxHeat; }

inline void CTNE2EulerSolver::SetTotal_CLift(double val_Total_CLift) { Total_CLift = val_Total_CLift; }

inline void CTNE2EulerSolver::SetTotal_CDrag(double val_Total_CDrag) { Total_CDrag = val_Total_CDrag; }

inline double CTNE2EulerSolver::GetAllBound_CLift_Inv() { return AllBound_CLift_Inv; }

inline double CTNE2EulerSolver::GetAllBound_CDrag_Inv() { return AllBound_CDrag_Inv; }

inline double CTNE2EulerSolver::GetAllBound_CSideForce_Inv() { return AllBound_CSideForce_Inv; }

inline double CTNE2EulerSolver::GetAllBound_CEff_Inv() { return AllBound_CEff_Inv; }

inline double CTNE2EulerSolver::GetAllBound_CMx_Inv() { return AllBound_CMx_Inv; }

inline double CTNE2EulerSolver::GetAllBound_CMy_Inv() { return AllBound_CMy_Inv; }

inline double CTNE2EulerSolver::GetAllBound_CMz_Inv() { return AllBound_CMz_Inv; }

inline double CTNE2EulerSolver::GetAllBound_CFx_Inv() { return AllBound_CFx_Inv; }

inline double CTNE2EulerSolver::GetAllBound_CFy_Inv() { return AllBound_CFy_Inv; }

inline double CTNE2EulerSolver::GetAllBound_CFz_Inv() { return AllBound_CFz_Inv; }

inline double CTNE2NSSolver::GetCDrag_Visc(unsigned short val_marker) { return CDrag_Visc[val_marker]; }

inline double CTNE2NSSolver::GetCLift_Visc(unsigned short val_marker) { return CLift_Visc[val_marker]; }

inline double CTNE2NSSolver::GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex) { return CSkinFriction[val_marker][val_vertex]; }

inline double CTNE2NSSolver::GetViscosity_Inf(void) { return Viscosity_Inf; }

inline double CTNE2NSSolver::GetHeatFlux(unsigned short val_marker, unsigned short val_vertex) { return HeatFlux[val_marker][val_vertex]; }

inline double CTNE2NSSolver::GetAllBound_CLift_Visc() { return AllBound_CLift_Visc; }

inline double CTNE2NSSolver::GetAllBound_CDrag_Visc() { return AllBound_CDrag_Visc; }

inline double CAdjTNE2EulerSolver::GetPhi_Inf(unsigned short val_dim) { return Phi_Inf[val_dim]; }

inline double CAdjTNE2EulerSolver::GetPsiE_Inf(void) { return PsiE_Inf; }

inline double* CAdjTNE2EulerSolver::GetPsiRhos_Inf(void) { return PsiRho_Inf; }

inline double CAdjTNE2EulerSolver::GetTotal_Sens_Geo() { return Total_Sens_Geo; }

inline void CAdjTNE2EulerSolver::SetCSensitivity(unsigned short val_marker, unsigned short val_vertex, double val_sensitivity) {CSensitivity[val_marker][val_vertex]=val_sensitivity; }

inline double CAdjTNE2EulerSolver::GetCSensitivity(unsigned short val_marker, unsigned short val_vertex) { return CSensitivity[val_marker][val_vertex]; }

inline double CAdjTNE2EulerSolver::GetTotal_Sens_AoA() { return Total_Sens_AoA; }

inline double CAdjTNE2EulerSolver::GetTotal_Sens_Mach() { return Total_Sens_Mach; }

inline double CAdjTNE2EulerSolver::GetTotal_Sens_Press() { return Total_Sens_Press; }

inline double CAdjTNE2EulerSolver::GetTotal_Sens_Temp() { return Total_Sens_Temp; }

