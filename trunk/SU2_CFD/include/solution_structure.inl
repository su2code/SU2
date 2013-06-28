/*!
 * \file solution_structure.inl
 * \brief In-Line subroutines of the <i>solution_structure.hpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.
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

inline void CSolution::SetIterLinSolver(unsigned short val_iterlinsolver) { IterLinSolver = val_iterlinsolver; }

inline unsigned short CSolution::GetIterLinSolver(void) { return IterLinSolver; }

inline double CSolution::GetCSensitivity(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline void CSolution::SetResidual_DualTime(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep, 
																		 unsigned short iMesh, unsigned short RunTime_EqSystem) { }
																		 
inline void CSolution::SetLevelSet_Distance(CGeometry *geometry, CConfig *config) { }

inline void CSolution::SetFEA_Load(CSolution ***flow_solution, CGeometry **fea_geometry, CGeometry **flow_geometry, CConfig *fea_config, CConfig *flow_config) { }

inline void CSolution::SetInitialCondition(CGeometry **geometry, CSolution ***solution_container, CConfig *config, unsigned long ExtIter) { }

inline void CSolution::GetRestart(CGeometry *geometry, CConfig *config, unsigned short val_iZone) { }

inline void CSolution::SetNoise_Source(CSolution ***flow_solution, CGeometry **wave_geometry, CConfig *wave_config) { }

inline void CSolution::SetAeroacoustic_Coupling(CSolution ***wave_solution, CSolution ***flow_solution, CNumerics *solver, CGeometry **flow_geometry, CConfig *flow_config) { }

inline void CSolution::SetFlow_Displacement(CGeometry **flow_geometry, CVolumetricMovement *flow_grid_movement, CConfig *flow_config, CConfig *fea_config, CGeometry **fea_geometry, CSolution ***fea_solution) { }

inline void CSolution::SetCSensitivity(unsigned short val_marker, unsigned short val_vertex, double val_sensitivity) { }

inline void CSolution::Inviscid_Sensitivity(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config) { }

inline void CSolution::Smooth_Sensitivity(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config) { }

inline void CSolution::Viscous_Sensitivity(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config) { }

inline double CSolution::GetPhi_Inf(unsigned short val_dim) { return 0; }

inline double CSolution::GetPsiRho_Inf(void) { return 0; }

inline double CSolution::GetPsiE_Inf(void) { return 0; }

inline void CSolution::SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config) { }

inline void CSolution::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config) { }

inline void CSolution::SetDistance(CGeometry *geometry, CConfig *config) { };

inline double CSolution::GetCLift_Inv(unsigned short val_marker) { return 0; }

inline double CSolution::GetCDrag_Inv(unsigned short val_marker) { return 0; }

inline double CSolution::GetMassFlow_Rate(unsigned short val_marker) { return 0; }

inline double CSolution::GetFanFace_Pressure(unsigned short val_marker) { return 0; }

inline double CSolution::GetFanFace_Mach(unsigned short val_marker) { return 0; }

inline double CSolution::GetCSideForce_Inv(unsigned short val_marker) { return 0; }

inline double CSolution::GetCEff_Inv(unsigned short val_marker) { return 0; }

inline double CSolution::GetCLift_Visc(unsigned short val_marker) { return 0; }

inline double CSolution::GetCDrag_Visc(unsigned short val_marker) { return 0; }

inline double CSolution::GetAllBound_CLift_Inv() { return 0; }

inline double CSolution::GetAllBound_CDrag_Inv() { return 0; }

inline double CSolution::GetAllBound_CSideForce_Inv() { return 0; }

inline double CSolution::GetAllBound_CEff_Inv() { return 0; }

inline double CSolution::GetAllBound_CLift_Visc() { return 0; }

inline double CSolution::GetAllBound_CDrag_Visc() { return 0; }

inline void CSolution::SetForceProj_Vector(CGeometry *geometry, CSolution **solution_container, CConfig *config) { }

inline void CSolution::SetIntBoundary_Jump(CGeometry *geometry, CSolution **solution_container, CConfig *config) { }

inline double CSolution::GetTotal_CLift() { return 0; }

inline double CSolution::GetTotal_CDrag() { return 0; }

inline double CSolution::GetTotal_CMx() { return 0; }

inline double CSolution::GetTotal_CMy() { return 0; }

inline double CSolution::GetTotal_CMz() { return 0; }

inline double CSolution::GetTotal_CFx() { return 0; }

inline double CSolution::GetTotal_CFy() { return 0; }

inline double CSolution::GetTotal_CFz() { return 0; }

inline double CSolution::GetTotal_CSideForce() { return 0; }

inline double CSolution::GetTotal_CEff() { return 0; }

inline double CSolution::GetTotal_CT() { return 0; }

inline void CSolution::SetTotal_CT(double val_Total_CT) { }

inline double CSolution::GetTotal_CQ() { return 0; }

inline void CSolution::SetTotal_CQ(double val_Total_CQ) { }

inline double CSolution::GetTotal_CMerit() { return 0; }

inline double CSolution::GetTotal_CEquivArea() { return 0; }

inline double CSolution::GetTotal_CFreeSurface() { return 0; }

inline double CSolution::GetTotal_CFEA() { return 0; }

inline double CSolution::GetTotal_CNearFieldOF() { return 0; }

inline void CSolution::SetTotal_CEquivArea(double val_cequivarea) { }

inline void CSolution::SetTotal_CFEA(double val_cfea) { }

inline void CSolution::SetTotal_CFreeSurface(double val_freesurface) { }

inline void CSolution::SetTotal_CNearFieldOF(double val_cnearfieldpress) { }

inline double CSolution::GetTotal_CDeltaLift() { return 0; }

inline double CSolution::GetTotal_CDeltaDrag() { return 0; }

inline double CSolution::GetTotal_CCharge() { return 0; }

inline double CSolution::GetTotal_CWave() { return 0; }

inline double CSolution::GetTotal_CHeat() { return 0; }

inline void CSolution::SetTotal_CLift(double val_Total_CLift) { }

inline void CSolution::SetTotal_CDrag(double val_Total_CDrag) { }

inline void CSolution::SetTotal_CCharge(double val_Total_CCharge) { }

inline double CSolution::GetCPressure(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolution::GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolution::GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline void CSolution::Viscous_Residual(CGeometry *geometry, CSolution **solution_container, 
                                        CNumerics *solver, CConfig *config, unsigned short iMesh, unsigned short iRKstep) { }

inline void CSolution::Galerkin_Method(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
									   CConfig *config, unsigned short iMesh) { }
									   
inline void CSolution::AddStiffMatrix(double ** StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1, unsigned long Point_2, unsigned long Point_3) { }
									   
inline void CSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, 
												  CNumerics *solver, CConfig *config, unsigned short iMesh) { }

inline void CSolution::SourceConserv_Residual(CGeometry *geometry, CSolution **solution_container, 
											   CNumerics *solver, CConfig *config, unsigned short iMesh) { }

inline void CSolution::Charge_Dist_SourceTerm(CGeometry *geometry, CSolution **solution_container, 
											  CNumerics *solver, CConfig *config, unsigned short iMesh) { }

inline double CSolution::GetTotal_Sens_Geo() { return 0; }

inline double CSolution::GetTotal_Sens_Mach() { return 0; }

inline double CSolution::GetTotal_Sens_AoA() { return 0; }

inline double CSolution::GetTotal_Sens_Press() { return 0; }

inline double CSolution::GetTotal_Sens_Temp() { return 0; }

inline double CSolution::GetDensity_Inf(void) { return 0; }

inline double CSolution::GetDensity_Inf(unsigned short val_var) { return 0; }

inline double CSolution::GetModVelocity_Inf(void) { return 0; }

inline double CSolution::GetDensity_Energy_Inf(void) { return 0; }

inline double CSolution::GetDensity_Energy_Inf(unsigned short val_var) { return 0; }

inline double CSolution::GetDensity_Energy_vib_Inf(unsigned short val_var) { return 0; }

inline double CSolution::GetDensity_Velocity_Inf(unsigned short val_dim) { return 0; }

inline double CSolution::GetDensity_Velocity_Inf(unsigned short val_dim, unsigned short val_var) { return 0; }

inline double CSolution::GetVelocity_Inf(unsigned short val_dim) { return 0; }

inline double CSolution::GetPressure_Inf(void) { return 0; }

inline double CSolution::GetViscosity_Inf(void) { return 0; }

inline double CSolution::GetDensity_Inlet(void) { return 0; }

inline double CSolution::GetDensity_Velocity_Inlet(unsigned short val_dim) { return 0; }

inline double CSolution::GetDensity_Energy_Inlet(void) { return 0; }

inline double CSolution::GetDensity_Outlet(void) { return 0; }

inline double CSolution::GetDensity_Outlet(unsigned short val_Fluid) { return 0; }

inline double CSolution::GetDensity_Inlet(unsigned short val_Fluid) { return 0; }

inline double CSolution::GetDensity_Velocity_Outlet(unsigned short val_dim) { return 0; }

inline double CSolution::GetDensity_Velocity_Outlet(unsigned short val_dim,unsigned short val_Fluid) { return 0; }

inline double CSolution::GetDensity_Velocity_Inlet(unsigned short val_dim,unsigned short val_Fluid) { return 0; }

inline double CSolution::GetDensity_Energy_Outlet(void) { return 0; }

inline double CSolution::GetDensity_Energy_Outlet(unsigned short val_Fluid) { return 0; }

inline double CSolution::GetDensity_Energy_Inlet(unsigned short val_Fluid) { return 0; }

inline double* CSolution::GetConstants() {return NULL;}

inline void CSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
									 unsigned short val_marker) { }
									 
inline void CSolution::BC_Displacement(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
									 unsigned short val_marker) { }
									 									 
inline void CSolution::BC_FlowLoad(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
									 unsigned short val_marker) { }
									 
inline void CSolution::BC_Load(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
									 unsigned short val_marker) { }

inline void CSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
								  unsigned short val_marker) { }
									
inline void CSolution::BC_Dirichlet(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
								  unsigned short val_marker) { }

inline void CSolution::MPI_Send_Receive(CGeometry ***geometry, CSolution ****solution_container,
                                        CConfig **config, unsigned short iMGLevel, unsigned short iZone) { }

inline void CSolution::BC_Interface_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
									CConfig *config, unsigned short val_marker) { }
                  
inline void CSolution::BC_NearField_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
									CConfig *config, unsigned short val_marker) { }
										
inline void CSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
								    CConfig *config, unsigned short val_marker) { }

inline void CSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
									CConfig *config, unsigned short val_marker) { cout << "inline" <<endl; }
									
inline void CSolution::BC_Custom(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
										 CConfig *config, unsigned short val_marker) { }
										 
inline void CSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
										 CConfig *config, unsigned short val_marker) { }

inline void CSolution::BC_Supersonic_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
										 CConfig *config, unsigned short val_marker) { }
                     
inline void CSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
										  CConfig *config, unsigned short val_marker) { }

inline void CSolution::BC_NacelleInflow(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
										  CConfig *config, unsigned short val_marker) { }
											
inline void CSolution::BC_NacelleExhaust(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
										  CConfig *config, unsigned short val_marker) { }
																						
inline void CSolution::BC_Neumann(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
										  CConfig *config, unsigned short val_marker) { }
								  
inline void CSolution::BC_Dielectric(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
									 CConfig *config, unsigned short val_marker) { }
									 										  
inline void CSolution::BC_Electrode(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
									CConfig *config, unsigned short val_marker) { }

inline void CSolution::BC_FWH(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
              unsigned short val_marker) { }
   
inline void CSolution::BC_Observer(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
            unsigned short val_marker) { }
                         
inline void CSolution::SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
							        unsigned short iMesh, unsigned long Iteration) { }	
							        
inline void CSolution::Postprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
							        unsigned short iMesh) { }								        

inline void CSolution::Centered_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
										CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

inline void CSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
									   CConfig *config, unsigned short iMesh) { }

inline void CSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, 
									 unsigned short iRKStep) { }

inline void CSolution::SetDissipation_Switch(CGeometry *geometry, CSolution **solution_container, CConfig *config) { }

inline void CSolution::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) { }

inline void CSolution::Inviscid_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolution::Viscous_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolution::Inviscid_DeltaForces(CGeometry *geometry, CSolution **solution_container, CConfig *config) { }

inline void CSolution::Viscous_DeltaForces(CGeometry *geometry, CConfig *config) { }

inline void CSolution::Wave_Strength(CGeometry *geometry, CConfig *config) { }

inline void CSolution::ExplicitRK_Iteration(CGeometry *geometry, CSolution **solution_container, 
											CConfig *config, unsigned short iRKStep) { }

inline void CSolution::ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) { }

inline void CSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) { }

inline void CSolution::Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
										  unsigned short iMesh) { }

inline void CSolution::Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config) { }

inline void CSolution::Compute_Residual(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
										unsigned short iMesh) { }

inline void CSolution::SetRes_Max(unsigned short val_var, double val_residual) { Residual_Max[val_var] = val_residual; }

inline void CSolution::AddRes_Max(unsigned short val_var, double val_residual) { Residual_Max[val_var] += val_residual; }

inline double CSolution::GetRes_Max(unsigned short val_var) { return Residual_Max[val_var]; }

inline void CSolution::Set_OldSolution(CGeometry *geometry) {
	for(unsigned long iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) 
		node[iPoint]->Set_OldSolution(); 
}

inline unsigned short CSolution::GetnVar(void) { return nVar; }

inline double CSolution::GetMax_Delta_Time(void) { return Max_Delta_Time; }

inline double CSolution::GetMin_Delta_Time(void) { return Min_Delta_Time; }

inline unsigned short CSolution::GetnSpecies(void) { return nSpecies; }

inline unsigned short CSolution::GetnMonatomics(void) { return nMonatomics; }

inline unsigned short CSolution::GetnDiatomics(void) { return nDiatomics; }

inline void CSolution::Copy_Zone_Solution(CSolution ***solver1_solution, CGeometry **solver1_geometry, CConfig *solver1_config, 
										  CSolution ***solver2_solution, CGeometry **solver2_geometry, CConfig *solver2_config) {};

inline double CEulerSolution::GetDensity_Inlet(void) { return Density_Inlet; }

inline double CEulerSolution::GetDensity_Energy_Inlet(void) { return Density_Inlet*Energy_Inlet; }

inline double CEulerSolution::GetDensity_Velocity_Inlet(unsigned short val_dim) { return Density_Inlet*Velocity_Inlet[val_dim]; }

inline double CEulerSolution::GetDensity_Outlet(void) { return Density_Outlet; }

inline double CEulerSolution::GetDensity_Energy_Outlet(void) { return Density_Outlet*Energy_Outlet; }

inline double CEulerSolution::GetDensity_Velocity_Outlet(unsigned short val_dim) { return Density_Outlet*Velocity_Outlet[val_dim]; }

inline double CEulerSolution::GetDensity_Inf(void) { return Density_Inf; }

inline double CEulerSolution::GetDensity_Back(void) { return Density_Back; }

inline double CEulerSolution::GetModVelocity_Inf(void) { 
	double Vel2 = 0; 
	for (unsigned short iDim = 0; iDim < nDim; iDim++) 
		Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim]; 
	return sqrt(Vel2);
}

inline double CEulerSolution::GetDensity_Energy_Inf(void) { return Density_Inf*Energy_Inf; }

inline double CEulerSolution::GetDensity_Energy_Back(void) { return Density_Back*Energy_Back; }

inline double CEulerSolution::GetDensity_Velocity_Inf(unsigned short val_dim) { return Density_Inf*Velocity_Inf[val_dim]; }

inline double CEulerSolution::GetDensity_Velocity_Back(unsigned short val_dim) { return Density_Back*Velocity_Back[val_dim]; }

inline double CEulerSolution::GetVelocity_Inf(unsigned short val_dim) { return Velocity_Inf[val_dim]; }

inline double CEulerSolution::GetPressure_Inf(void) { return Pressure_Inf; }

inline double CEulerSolution::GetPressure_Back(void) { return Pressure_Back; }

inline double CEulerSolution::GetCPressure(unsigned short val_marker, unsigned short val_vertex) { return CPressure[val_marker][val_vertex]; }

inline double CEulerSolution::GetCLift_Inv(unsigned short val_marker) { return CLift_Inv[val_marker]; }

inline double CEulerSolution::GetCDrag_Inv(unsigned short val_marker) { return CDrag_Inv[val_marker]; }

inline double CEulerSolution::GetMassFlow_Rate(unsigned short val_marker) { return MassFlow_Rate[val_marker]; }

inline double CEulerSolution::GetFanFace_Pressure(unsigned short val_marker) { return FanFace_Pressure[val_marker]; }

inline double CEulerSolution::GetFanFace_Mach(unsigned short val_marker) { return FanFace_Mach[val_marker]; }

inline double CEulerSolution::GetCSideForce_Inv(unsigned short val_marker) { return CSideForce_Inv[val_marker]; }

inline double CEulerSolution::GetCEff_Inv(unsigned short val_marker) { return CEff_Inv[val_marker]; }

inline double CEulerSolution::GetTotal_CLift() { return Total_CLift; }

inline double CEulerSolution::GetTotal_CDrag() { return Total_CDrag; }

inline double CEulerSolution::GetTotal_CMx() { return Total_CMx; }

inline double CEulerSolution::GetTotal_CMy() { return Total_CMy; }

inline double CEulerSolution::GetTotal_CMz() { return Total_CMz; }

inline double CEulerSolution::GetTotal_CFx() { return Total_CFx; }

inline double CEulerSolution::GetTotal_CFy() { return Total_CFy; }

inline double CEulerSolution::GetTotal_CFz() { return Total_CFz; }

inline double CEulerSolution::GetTotal_CSideForce() { return Total_CSideForce; }

inline double CEulerSolution::GetTotal_CEff() { return Total_CEff; }

inline double CEulerSolution::GetTotal_CT() { return Total_CT; }

inline void CEulerSolution::SetTotal_CT(double val_Total_CT) { Total_CT = val_Total_CT; }

inline double CEulerSolution::GetTotal_CQ() { return Total_CQ; }

inline void CEulerSolution::SetTotal_CQ(double val_Total_CQ) { Total_CQ = val_Total_CQ; }

inline double CEulerSolution::GetTotal_CMerit() { return Total_CMerit; }

inline double CEulerSolution::GetTotal_CEquivArea() { return Total_CEquivArea; }

inline double CEulerSolution::GetTotal_CNearFieldOF() { return Total_CNearFieldOF; }

inline void CEulerSolution::SetTotal_CEquivArea(double cequivarea) { Total_CEquivArea = cequivarea; }

inline void CEulerSolution::SetTotal_CNearFieldOF(double cnearfieldpress) { Total_CNearFieldOF = cnearfieldpress; }

inline void CEulerSolution::SetTotal_CLift(double val_Total_CLift) { Total_CLift = val_Total_CLift; }

inline void CEulerSolution::SetTotal_CDrag(double val_Total_CDrag) { Total_CDrag = val_Total_CDrag; }

inline double CEulerSolution::GetAllBound_CLift_Inv() { return AllBound_CLift_Inv; }

inline double CEulerSolution::GetAllBound_CDrag_Inv() { return AllBound_CDrag_Inv; }

inline double CEulerSolution::GetAllBound_CSideForce_Inv() { return AllBound_CSideForce_Inv; }

inline double CEulerSolution::GetAllBound_CEff_Inv() { return AllBound_CEff_Inv; }

inline double CNSSolution::GetViscosity_Inf(void) { return Viscosity_Inf; }

inline double CNSSolution::GetCLift_Visc(unsigned short val_marker) { return CLift_Visc[val_marker]; }

inline double CNSSolution::GetCDrag_Visc(unsigned short val_marker) { return CDrag_Visc[val_marker]; }

inline double CNSSolution::GetAllBound_CLift_Visc() { return AllBound_CLift_Visc; }

inline double CNSSolution::GetAllBound_CDrag_Visc() { return AllBound_CDrag_Visc; }

inline double CNSSolution::GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex) { return CSkinFriction[val_marker][val_vertex]; }

inline double CNSSolution::GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_vertex) { return CHeatTransfer[val_marker][val_vertex]; }

inline double CAdjEulerSolution::GetCSensitivity(unsigned short val_marker, unsigned short val_vertex) { return CSensitivity[val_marker][val_vertex]; }

inline void CAdjEulerSolution::SetCSensitivity(unsigned short val_marker, unsigned short val_vertex, double val_sensitivity) {CSensitivity[val_marker][val_vertex]=val_sensitivity; }

inline double CAdjEulerSolution::GetTotal_Sens_Geo() { return Total_Sens_Geo; }

inline double CAdjEulerSolution::GetTotal_Sens_Mach() { return Total_Sens_Mach; }

inline double CAdjEulerSolution::GetTotal_Sens_AoA() { return Total_Sens_AoA; }

inline double CAdjEulerSolution::GetTotal_Sens_Press() { return Total_Sens_Press; }

inline double CAdjEulerSolution::GetTotal_Sens_Temp() { return Total_Sens_Temp; }

inline double CAdjEulerSolution::GetPsiRho_Inf(void) { return PsiRho_Inf; }

inline double CAdjEulerSolution::GetPsiE_Inf(void) { return PsiE_Inf; }

inline double CAdjEulerSolution::GetPhi_Inf(unsigned short val_dim) { return Phi_Inf[val_dim]; }

inline double CLinEulerSolution::GetTotal_CDeltaLift() { return Total_CDeltaLift; }

inline double CLinEulerSolution::GetTotal_CDeltaDrag() { return Total_CDeltaDrag; }

inline double CElectricSolution::GetTotal_CCharge() { return Total_CCharge; }

inline void CElectricSolution::SetTotal_CCharge(double val_Total_CCharge) {Total_CCharge = val_Total_CCharge; }

inline double CPlasmaSolution::GetDensity_Energy_vib_Inf(unsigned short  iSpecies) { return Density_Inf[iSpecies]*Energy_vib_Inf[iSpecies]; }

inline double CPlasmaSolution::GetDensity_Energy_Inf(unsigned short  iSpecies) { return Density_Inf[iSpecies]*Energy_Inf[iSpecies]; }

inline double CPlasmaSolution::GetDensity_Velocity_Inf(unsigned short val_dim, unsigned short val_Fluid) { return Density_Inf[val_Fluid]*Velocity_Inf[val_Fluid][val_dim]; }

inline double CPlasmaSolution::GetDensity_Inf(unsigned short iSpecies) { return Density_Inf[iSpecies]; }

inline double CPlasmaSolution::GetDensity_Outlet(unsigned short val_Fluid) { return Density_Outlet[val_Fluid]; }

inline double CPlasmaSolution::GetDensity_Velocity_Outlet(unsigned short val_dim,unsigned short val_Fluid) { return Density_Outlet[val_Fluid]*Velocity_Outlet[val_Fluid][val_dim]; }

inline double CPlasmaSolution::GetDensity_Energy_Outlet(unsigned short val_Fluid) { return Density_Outlet[val_Fluid]*Energy_Outlet[val_Fluid]; }

inline double CPlasmaSolution::GetDensity_Inlet(unsigned short val_Fluid) { return Density_Inlet[val_Fluid]; }

inline double CPlasmaSolution::GetDensity_Velocity_Inlet(unsigned short val_dim,unsigned short val_Fluid) { return Density_Inlet[val_Fluid]*Velocity_Inlet[val_Fluid][val_dim]; }

inline double CPlasmaSolution::GetDensity_Energy_Inlet(unsigned short val_Fluid) { return Density_Inlet[val_Fluid]*Energy_Inlet[val_Fluid]; }

inline double CPlasmaSolution::GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex) { return CSkinFriction[val_marker][val_vertex]; }

inline double CPlasmaSolution::GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_vertex) { return CHeatTransfer[val_marker][val_vertex]; }

inline double CPlasmaSolution::GetTotal_CLift() { return Total_CLift; }

inline double CPlasmaSolution::GetTotal_CDrag() { return Total_CDrag; }

inline double CPlasmaSolution::GetTotal_CMx() { return Total_CMx; }

inline double CPlasmaSolution::GetTotal_CMy() { return Total_CMy; }

inline double CPlasmaSolution::GetTotal_CMz() { return Total_CMz; }

inline double CPlasmaSolution::GetTotal_CFx() { return Total_CFx; }

inline double CPlasmaSolution::GetTotal_CFy() { return Total_CFy; }

inline double CPlasmaSolution::GetTotal_CFz() { return Total_CFz; }

inline double CPlasmaSolution::GetTotal_CSideForce() { return Total_CSideForce; }

inline double CPlasmaSolution::GetTotal_CEff() { return Total_CEff; }

inline double CLevelSetSolution::GetTotal_CFreeSurface() { return Total_CFreeSurface; }

inline void CLevelSetSolution::SetTotal_CFreeSurface(double cfreesurface) { Total_CFreeSurface = cfreesurface; }

inline double CAdjPlasmaSolution::GetCSensitivity(unsigned short val_marker, unsigned short val_vertex) { return CSensitivity[val_marker][val_vertex]; }

inline void CAdjPlasmaSolution::SetCSensitivity(unsigned short val_marker, unsigned short val_vertex, double val_sensitivity) {CSensitivity[val_marker][val_vertex]=val_sensitivity; }

inline double CAdjPlasmaSolution::GetTotal_Sens_Geo() { return Total_Sens_Geo; }

inline double CAdjPlasmaSolution::GetTotal_Sens_Mach() { return Total_Sens_Mach; }

inline double CAdjPlasmaSolution::GetTotal_Sens_AoA() { return Total_Sens_AoA; }

inline double CAdjPlasmaSolution::GetTotal_Sens_Press() { return Total_Sens_Press; }

inline double CAdjPlasmaSolution::GetTotal_Sens_Temp() { return Total_Sens_Temp; }

inline double CAdjPlasmaSolution::GetPsiRho_Inf(void) { return PsiRho_Inf; }

inline double CAdjPlasmaSolution::GetPsiE_Inf(void) { return PsiE_Inf; }

inline double CAdjPlasmaSolution::GetPhi_Inf(unsigned short val_dim) { return Phi_Inf[val_dim]; }

inline double CFEASolution::GetTotal_CFEA() { return Total_CFEA; }

inline void CFEASolution::SetTotal_CFEA(double cfea) { Total_CFEA = cfea; }

inline double CWaveSolution::GetTotal_CWave() { return Total_CWave; }

inline double CHeatSolution::GetTotal_CHeat() { return Total_CHeat; }
