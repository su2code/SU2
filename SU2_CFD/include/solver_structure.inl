/*!
 * \file solver_structure.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
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

inline void CSolver::SetIterLinSolver(unsigned short val_iterlinsolver) { IterLinSolver = val_iterlinsolver; }

inline unsigned short CSolver::GetnSpecies(void) { return 0; }

inline unsigned short CSolver::GetnMonatomics(void) { return 0; }

inline unsigned short CSolver::GetnDiatomics(void) { return 0; }

inline void CSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) { }

inline unsigned short CSolver::GetIterLinSolver(void) { return IterLinSolver; }

inline double CSolver::GetCSensitivity(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline void CSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep, 
																		 unsigned short iMesh, unsigned short RunTime_EqSystem) { }
																		 
inline void CSolver::SetLevelSet_Distance(CGeometry *geometry, CConfig *config, bool Initialization, bool WriteLevelSet) { }

inline void CSolver::SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry, CGeometry **flow_geometry, CConfig *fea_config, CConfig *flow_config) { }

inline void CSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) { }

inline void CSolver::GetRestart(CGeometry *geometry, CConfig *config, unsigned short val_iZone) { }
  
inline void CSolver::SetNoise_Source(CSolver ***flow_solution, CGeometry **wave_geometry, CConfig *wave_config) { }

inline void CSolver::SetAeroacoustic_Coupling(CSolver ***wave_solution, CSolver ***flow_solution, CNumerics *numerics, CGeometry **flow_geometry, CConfig *flow_config) { }

inline void CSolver::SetFlow_Displacement(CGeometry **flow_geometry, CVolumetricMovement *flow_grid_movement, CConfig *flow_config, CConfig *fea_config, CGeometry **fea_geometry, CSolver ***fea_solution) { }

inline void CSolver::SetCSensitivity(unsigned short val_marker, unsigned short val_vertex, double val_sensitivity) { }

inline void CSolver::Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline void CSolver::Smooth_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline void CSolver::Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

inline double CSolver::GetPhi_Inf(unsigned short val_dim) { return 0; }

inline double CSolver::GetPsiRho_Inf(void) { return 0; }

inline double CSolver::GetPsiE_Inf(void) { return 0; }

inline void CSolver::SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Set_MPI_PrimVar_Gradient(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPrimVar_Limiter_MPI(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPrimVar_Limiter(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPreconditioner(CConfig *config, unsigned short iPoint) { }

inline void CSolver::SetDistance(CGeometry *geometry, CConfig *config) { };

inline double CSolver::GetCLift_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetCDrag_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetFanFace_MassFlow(unsigned short val_marker) { return 0; }

inline double CSolver::GetExhaust_MassFlow(unsigned short val_marker) { return 0; }

inline double CSolver::GetFanFace_Pressure(unsigned short val_marker) { return 0; }

inline double CSolver::GetFanFace_Mach(unsigned short val_marker) { return 0; }

inline double CSolver::GetCSideForce_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetCEff_Inv(unsigned short val_marker) { return 0; }

inline double CSolver::GetCLift_Visc(unsigned short val_marker) { return 0; }

inline double CSolver::GetCDrag_Visc(unsigned short val_marker) { return 0; }

inline double CSolver::GetAllBound_CLift_Inv() { return 0; }

inline double CSolver::GetAllBound_CDrag_Inv() { return 0; }

inline double CSolver::GetAllBound_CSideForce_Inv() { return 0; }

inline double CSolver::GetAllBound_CEff_Inv() { return 0; }

inline double CSolver::GetAllBound_CLift_Visc() { return 0; }

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

inline double CSolver::GetTotal_Q() { return 0; }

inline double CSolver::GetTotal_MaxQ() { return 0; }

inline double CSolver::Get_PressureDrag() { return 0; }

inline double CSolver::Get_ViscDrag() { return 0; }

inline double CSolver::Get_MagnetDrag() { return 0; }

inline void CSolver::SetTotal_CQ(double val_Total_CQ) { }

inline void CSolver::SetTotal_Q(double val_Total_Q) { }

inline void CSolver::SetTotal_MaxQ(double val_Total_Q) { }

inline double CSolver::GetTotal_CMerit() { return 0; }

inline double CSolver::GetTotal_CEquivArea() { return 0; }

inline double CSolver::GetTotal_CFreeSurface() { return 0; }

inline double CSolver::GetTotal_CFEA() { return 0; }

inline double CSolver::GetTotal_CNearFieldOF() { return 0; }

inline void CSolver::SetTotal_CEquivArea(double val_cequivarea) { }

inline void CSolver::SetTotal_CFEA(double val_cfea) { }

inline void CSolver::SetTotal_CFreeSurface(double val_freesurface) { }

inline void CSolver::SetTotal_CNearFieldOF(double val_cnearfieldpress) { }

inline double CSolver::GetTotal_CDeltaLift() { return 0; }

inline double CSolver::GetTotal_CDeltaDrag() { return 0; }

inline double CSolver::GetTotal_CCharge() { return 0; }

inline double CSolver::GetTotal_CWave() { return 0; }

inline double CSolver::GetTotal_CHeat() { return 0; }

inline void CSolver::SetTotal_CLift(double val_Total_CLift) { }

inline void CSolver::SetTotal_CDrag(double val_Total_CDrag) { }

inline void CSolver::SetTotal_CCharge(double val_Total_CCharge) { }

inline double CSolver::GetCPressure(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolver::GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolver::GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolver::GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_iSpecies, unsigned short val_vertex) { return 0; }

inline double CSolver::GetViscForce(unsigned short val_marker, unsigned short val_iSpecies, unsigned short iDim, unsigned short val_vertex) { return 0; }

inline double CSolver::GetPressureForce(unsigned short val_marker, unsigned short val_iSpecies, unsigned short iDim, unsigned short val_vertex) { return 0; }

inline double CSolver::GetYPlus(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline void CSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, 
                                        CNumerics *numerics, CConfig *config, unsigned short iMesh, unsigned short iRKstep) { }

inline void CSolver::Galerkin_Method(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									   CConfig *config, unsigned short iMesh) { }
									   
inline void CSolver::AddStiffMatrix(double ** StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1, unsigned long Point_2, unsigned long Point_3) { }
									   
inline void CSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, 
												  CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) { }
									   
inline void CSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, 
												  CNumerics *numerics, CConfig *config, unsigned short iMesh) { }
												  
inline void CSolver::SourceConserv_Residual(CGeometry *geometry, CSolver **solver_container, 
											   CNumerics *numerics, CConfig *config, unsigned short iMesh) { }

inline void CSolver::Charge_Dist_SourceTerm(CGeometry *geometry, CSolver **solver_container, 
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

inline double CSolver::GetDensity_Energy_Inf(unsigned short val_var) { return 0; }

inline double CSolver::GetDensity_Energy_vib_Inf(unsigned short val_var) { return 0; }

inline double CSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return 0; }

inline double CSolver::GetDensity_Velocity_Inf(unsigned short val_dim, unsigned short val_var) { return 0; }

inline double CSolver::GetVelocity_Inf(unsigned short val_dim) { return 0; }

inline double CSolver::GetPressure_Inf(void) { return 0; }

inline double CSolver::GetViscosity_Inf(void) { return 0; }

inline double CSolver::GetDensity_Inlet(void) { return 0; }

inline double CSolver::GetDensity_Velocity_Inlet(unsigned short val_dim) { return 0; }

inline double CSolver::GetDensity_Energy_Inlet(void) { return 0; }

inline double CSolver::GetDensity_Outlet(void) { return 0; }

inline double CSolver::GetDensity_Outlet(unsigned short val_Fluid) { return 0; }

inline double CSolver::GetDensity_Inlet(unsigned short val_Fluid) { return 0; }

inline double CSolver::GetDensity_Velocity_Outlet(unsigned short val_dim) { return 0; }

inline double CSolver::GetDensity_Velocity_Outlet(unsigned short val_dim,unsigned short val_Fluid) { return 0; }

inline double CSolver::GetDensity_Velocity_Inlet(unsigned short val_dim,unsigned short val_Fluid) { return 0; }

inline double CSolver::GetDensity_Energy_Outlet(void) { return 0; }

inline double CSolver::GetDensity_Energy_Outlet(unsigned short val_Fluid) { return 0; }

inline double CSolver::GetDensity_Energy_Inlet(unsigned short val_Fluid) { return 0; }

inline double* CSolver::GetConstants() {return NULL;}

inline void CSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }
									 
inline void CSolver::BC_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }
									 									 
inline void CSolver::BC_FlowLoad(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }
									 
inline void CSolver::BC_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
									 unsigned short val_marker) { }
                  
inline void CSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }
                  
inline void CSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }
									
inline void CSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
								  unsigned short val_marker) { }

inline void CSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									CConfig *config, unsigned short val_marker) { }
                  
inline void CSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									CConfig *config, unsigned short val_marker) { }
										
inline void CSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
								    CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
									CConfig *config, unsigned short val_marker) { }
									
inline void CSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
										 CConfig *config, unsigned short val_marker) { }
										 
inline void CSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
										 CConfig *config, unsigned short val_marker) { }
                     
inline void CSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
										  CConfig *config, unsigned short val_marker) { }
                      
                      
inline void CSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
										 CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Nacelle_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
										  CConfig *config, unsigned short val_marker) { }
											
inline void CSolver::BC_Nacelle_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
										  CConfig *config, unsigned short val_marker) { }
																						
inline void CSolver::BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
										  CConfig *config, unsigned short val_marker) { }
								  
inline void CSolver::BC_Dielectric(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									 CConfig *config, unsigned short val_marker) { }
									 										  
inline void CSolver::BC_Electrode(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_FWH(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
              unsigned short val_marker) { }
   
inline void CSolver::BC_Observer(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
            unsigned short val_marker) { }
            
inline void CSolver::GetNacelle_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh) { }
                         
inline void CSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
							        unsigned short iMesh, unsigned long Iteration) { }	
							        
inline void CSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
							        unsigned short iMesh) { }								        

inline void CSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
										CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

inline void CSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
									   CConfig *config, unsigned short iMesh) { }

inline void CSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) { }

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

inline void CSolver::Solve_LinearSystem(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
										  unsigned short iMesh) { }

inline void CSolver::Solve_LinearSystem(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::Compute_Residual(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
										unsigned short iMesh) { }

inline void CSolver::SetRes_RMS(unsigned short val_var, double val_residual) { Residual_RMS[val_var] = val_residual; }

inline void CSolver::AddRes_RMS(unsigned short val_var, double val_residual) { Residual_RMS[val_var] += val_residual; }

inline double CSolver::GetRes_RMS(unsigned short val_var) { return Residual_RMS[val_var]; }

inline void CSolver::SetRes_Max(unsigned short val_var, double val_residual, unsigned long val_point) { Residual_Max[val_var] = val_residual; Point_Max[val_var] = val_point; }

inline void CSolver::AddRes_Max(unsigned short val_var, double val_residual, unsigned long val_point) { 
  if (val_residual > Residual_Max[val_var]) {
  Residual_Max[val_var] = val_residual;
  Point_Max[val_var] = val_point; }
}

inline double CSolver::GetRes_Max(unsigned short val_var) { return Residual_Max[val_var]; }

inline unsigned long CSolver::GetPoint_Max(unsigned short val_var) { return Point_Max[val_var]; }

inline void CSolver::Set_OldSolution(CGeometry *geometry) {
	for(unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) 
		node[iPoint]->Set_OldSolution(); // The loop should be over nPoints 
                                     //  to guarantee that the boundaries are
                                     //  well updated
}

inline unsigned short CSolver::GetnVar(void) { return nVar; }

inline unsigned short CSolver::GetnPrimVar(void) { return nPrimVar; }

inline double CSolver::GetMax_Delta_Time(void) { return Max_Delta_Time; }

inline double CSolver::GetMin_Delta_Time(void) { return Min_Delta_Time; }

inline double CSolver::GetMax_Delta_Time(unsigned short val_Species) { return 0.0; }

inline double CSolver::GetMin_Delta_Time(unsigned short val_Species) { return 0.0; }

inline void CSolver::Copy_Zone_Solution(CSolver ***solver1_solution, CGeometry **solver1_geometry, CConfig *solver1_config, 
										  CSolver ***solver2_solution, CGeometry **solver2_geometry, CConfig *solver2_config) {};

inline double CEulerSolver::GetDensity_Inlet(void) { return Density_Inlet; }

inline double CEulerSolver::GetDensity_Energy_Inlet(void) { return Density_Inlet*Energy_Inlet; }

inline double CEulerSolver::GetDensity_Velocity_Inlet(unsigned short val_dim) { return Density_Inlet*Velocity_Inlet[val_dim]; }

inline double CEulerSolver::GetDensity_Outlet(void) { return Density_Outlet; }

inline double CEulerSolver::GetDensity_Energy_Outlet(void) { return Density_Outlet*Energy_Outlet; }

inline double CEulerSolver::GetDensity_Velocity_Outlet(unsigned short val_dim) { return Density_Outlet*Velocity_Outlet[val_dim]; }

inline double CEulerSolver::GetDensity_Inf(void) { return Density_Inf; }

inline double CEulerSolver::GetDensity_Back(void) { return Density_Back; }

inline double CEulerSolver::GetModVelocity_Inf(void) { 
	double Vel2 = 0; 
	for (unsigned short iDim = 0; iDim < nDim; iDim++) 
		Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim]; 
	return sqrt(Vel2);
}

inline double CEulerSolver::GetDensity_Energy_Inf(void) { return Density_Inf*Energy_Inf; }

inline double CEulerSolver::GetDensity_Energy_Back(void) { return Density_Back*Energy_Back; }

inline double CEulerSolver::GetDensity_Velocity_Inf(unsigned short val_dim) { return Density_Inf*Velocity_Inf[val_dim]; }

inline double CEulerSolver::GetDensity_Velocity_Back(unsigned short val_dim) { return Density_Back*Velocity_Back[val_dim]; }

inline double CEulerSolver::GetVelocity_Inf(unsigned short val_dim) { return Velocity_Inf[val_dim]; }

inline double CEulerSolver::GetPressure_Inf(void) { return Pressure_Inf; }

inline double CEulerSolver::GetPressure_Back(void) { return Pressure_Back; }

inline double CEulerSolver::GetCPressure(unsigned short val_marker, unsigned short val_vertex) { return CPressure[val_marker][val_vertex]; }

inline double CEulerSolver::GetCLift_Inv(unsigned short val_marker) { return CLift_Inv[val_marker]; }

inline double CEulerSolver::GetCDrag_Inv(unsigned short val_marker) { return CDrag_Inv[val_marker]; }

inline double CEulerSolver::GetFanFace_MassFlow(unsigned short val_marker) { return FanFace_MassFlow[val_marker]; }

inline double CEulerSolver::GetExhaust_MassFlow(unsigned short val_marker) { return Exhaust_MassFlow[val_marker]; }

inline double CEulerSolver::GetFanFace_Pressure(unsigned short val_marker) { return FanFace_Pressure[val_marker]; }

inline double CEulerSolver::GetFanFace_Mach(unsigned short val_marker) { return FanFace_Mach[val_marker]; }

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

inline double CEulerSolver::GetTotal_Q() { return Total_Q; }

inline double CEulerSolver::GetTotal_MaxQ() { return Total_Maxq; }

inline void CEulerSolver::SetTotal_CQ(double val_Total_CQ) { Total_CQ = val_Total_CQ; }

inline void CEulerSolver::SetTotal_Q(double val_Total_Q) { Total_Q = val_Total_Q; }

inline void CEulerSolver::SetTotal_MaxQ(double val_Total_MaxQ) { Total_Maxq = val_Total_MaxQ; }

inline double CEulerSolver::GetTotal_CMerit() { return Total_CMerit; }

inline double CEulerSolver::GetTotal_CEquivArea() { return Total_CEquivArea; }

inline double CEulerSolver::GetTotal_CNearFieldOF() { return Total_CNearFieldOF; }

inline void CEulerSolver::SetTotal_CEquivArea(double cequivarea) { Total_CEquivArea = cequivarea; }

inline void CEulerSolver::SetTotal_CNearFieldOF(double cnearfieldpress) { Total_CNearFieldOF = cnearfieldpress; }

inline void CEulerSolver::SetTotal_CLift(double val_Total_CLift) { Total_CLift = val_Total_CLift; }

inline void CEulerSolver::SetTotal_CDrag(double val_Total_CDrag) { Total_CDrag = val_Total_CDrag; }

inline double CEulerSolver::GetAllBound_CLift_Inv() { return AllBound_CLift_Inv; }

inline double CEulerSolver::GetAllBound_CDrag_Inv() { return AllBound_CDrag_Inv; }

inline double CEulerSolver::GetAllBound_CSideForce_Inv() { return AllBound_CSideForce_Inv; }

inline double CEulerSolver::GetAllBound_CEff_Inv() { return AllBound_CEff_Inv; }

inline double CNSSolver::GetViscosity_Inf(void) { return Viscosity_Inf; }

inline double CNSSolver::GetCLift_Visc(unsigned short val_marker) { return CLift_Visc[val_marker]; }

inline double CNSSolver::GetCDrag_Visc(unsigned short val_marker) { return CDrag_Visc[val_marker]; }

inline double CNSSolver::GetAllBound_CLift_Visc() { return AllBound_CLift_Visc; }

inline double CNSSolver::GetAllBound_CDrag_Visc() { return AllBound_CDrag_Visc; }

inline double CNSSolver::GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex) { return CSkinFriction[val_marker][val_vertex]; }

inline double CNSSolver::GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_vertex) { return CHeatTransfer[val_marker][val_vertex]; }

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

inline double CElectricSolver::GetTotal_CCharge() { return Total_CCharge; }

inline void CElectricSolver::SetTotal_CCharge(double val_Total_CCharge) {Total_CCharge = val_Total_CCharge; }

inline unsigned short CPlasmaSolver::GetnSpecies(void) { return nSpecies; }

inline unsigned short CPlasmaSolver::GetnMonatomics(void) { return nMonatomics; }

inline unsigned short CPlasmaSolver::GetnDiatomics(void) { return nDiatomics; }

inline double CPlasmaSolver::GetDensity_Energy_vib_Inf(unsigned short  iSpecies) { return Density_Inf[iSpecies]*Energy_vib_Inf[iSpecies]; }

inline double CPlasmaSolver::GetDensity_Energy_Inf(unsigned short  iSpecies) { return Density_Inf[iSpecies]*Energy_Inf[iSpecies]; }

inline double CPlasmaSolver::GetDensity_Velocity_Inf(unsigned short val_dim, unsigned short val_Fluid) { return Density_Inf[val_Fluid]*Velocity_Inf[val_Fluid][val_dim]; }

inline double CPlasmaSolver::GetDensity_Inf(unsigned short iSpecies) { return Density_Inf[iSpecies]; }

inline double CPlasmaSolver::GetDensity_Outlet(unsigned short val_Fluid) { return Density_Outlet[val_Fluid]; }

inline double CPlasmaSolver::GetDensity_Velocity_Outlet(unsigned short val_dim,unsigned short val_Fluid) { return Density_Outlet[val_Fluid]*Velocity_Outlet[val_Fluid][val_dim]; }

inline double CPlasmaSolver::GetDensity_Energy_Outlet(unsigned short val_Fluid) { return Density_Outlet[val_Fluid]*Energy_Outlet[val_Fluid]; }

inline double CPlasmaSolver::GetDensity_Inlet(unsigned short val_Fluid) { return Density_Inlet[val_Fluid]; }

inline double CPlasmaSolver::GetDensity_Velocity_Inlet(unsigned short val_dim,unsigned short val_Fluid) { return Density_Inlet[val_Fluid]*Velocity_Inlet[val_Fluid][val_dim]; }

inline double CPlasmaSolver::GetDensity_Energy_Inlet(unsigned short val_Fluid) { return Density_Inlet[val_Fluid]*Energy_Inlet[val_Fluid]; }

inline double CPlasmaSolver::GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex) { return CSkinFriction[val_marker][val_vertex]; }

inline double CPlasmaSolver::GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_iSpecies, unsigned short val_vertex) { return CHeatTransfer[val_marker][val_iSpecies][val_vertex]; }

inline double CPlasmaSolver::GetViscForce(unsigned short val_marker, unsigned short val_iSpecies, unsigned short iDim, unsigned short val_vertex) { return CViscForce[val_marker][val_iSpecies][iDim][val_vertex]; }

inline double CPlasmaSolver::GetPressureForce(unsigned short val_marker, unsigned short val_iSpecies, unsigned short iDim, unsigned short val_vertex) { return CPressForce[val_marker][val_iSpecies][iDim][val_vertex]; }

inline double CPlasmaSolver::GetTotal_CLift() { return Total_CLift; }

inline double CPlasmaSolver::GetTotal_CDrag() { return Total_CDrag; }

inline double CPlasmaSolver::GetTotal_CMx() { return Total_CMx; }

inline double CPlasmaSolver::GetTotal_CMy() { return Total_CMy; }

inline double CPlasmaSolver::GetTotal_CMz() { return Total_CMz; }

inline double CPlasmaSolver::GetTotal_CFx() { return Total_CFx; }

inline double CPlasmaSolver::GetTotal_CFy() { return Total_CFy; }

inline double CPlasmaSolver::GetTotal_CFz() { return Total_CFz; }

inline double CPlasmaSolver::GetTotal_CSideForce() { return Total_CSideForce; }

inline double CPlasmaSolver::GetTotal_CEff() { return Total_CEff; }

inline double CPlasmaSolver::GetTotal_Q() { return Total_Q; }

inline double CPlasmaSolver::Get_PressureDrag() { return PressureDrag; }

inline double CPlasmaSolver::Get_ViscDrag() { return ViscDrag; }

inline double CPlasmaSolver::Get_MagnetDrag() { return MagnetDrag; }

inline double CPlasmaSolver::GetMax_Delta_Time(unsigned short val_Species) { return Max_Delta_Time[val_Species]; }

inline double CPlasmaSolver::GetMin_Delta_Time(unsigned short val_Species) { return Min_Delta_Time[val_Species]; }

inline double CLevelSetSolver::GetTotal_CFreeSurface() { return Total_CFreeSurface; }

inline void CLevelSetSolver::SetTotal_CFreeSurface(double cfreesurface) { Total_CFreeSurface = cfreesurface; }

inline unsigned short CAdjPlasmaSolver::GetnSpecies(void) { return nSpecies; }

inline unsigned short CAdjPlasmaSolver::GetnMonatomics(void) { return nMonatomics; }

inline unsigned short CAdjPlasmaSolver::GetnDiatomics(void) { return nDiatomics; }

inline double CAdjPlasmaSolver::GetCSensitivity(unsigned short val_marker, unsigned short val_vertex) { return CSensitivity[val_marker][val_vertex]; }

inline void CAdjPlasmaSolver::SetCSensitivity(unsigned short val_marker, unsigned short val_vertex, double val_sensitivity) {CSensitivity[val_marker][val_vertex]=val_sensitivity; }

inline double CAdjPlasmaSolver::GetTotal_Sens_Geo() { return Total_Sens_Geo; }

inline double CAdjPlasmaSolver::GetTotal_Sens_Mach() { return Total_Sens_Mach; }

inline double CAdjPlasmaSolver::GetTotal_Sens_AoA() { return Total_Sens_AoA; }

inline double CAdjPlasmaSolver::GetTotal_Sens_Press() { return Total_Sens_Press; }

inline double CAdjPlasmaSolver::GetTotal_Sens_Temp() { return Total_Sens_Temp; }

inline double CAdjPlasmaSolver::GetPsiRho_Inf(void) { return PsiRho_Inf; }

inline double CAdjPlasmaSolver::GetPsiE_Inf(void) { return PsiE_Inf; }

inline double CAdjPlasmaSolver::GetPhi_Inf(unsigned short val_dim) { return Phi_Inf[val_dim]; }

inline double CFEASolver::GetTotal_CFEA() { return Total_CFEA; }

inline void CFEASolver::SetTotal_CFEA(double cfea) { Total_CFEA = cfea; }

inline double CWaveSolver::GetTotal_CWave() { return Total_CWave; }

inline double CHeatSolver::GetTotal_CHeat() { return Total_CHeat; }

inline void CTurbSolver::CalcEddyViscosity(double *val_FlowVars, double val_laminar_viscosity,
			double *val_TurbVar, double *val_eddy_viscosity) {}


