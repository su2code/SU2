/*!
 * \file integration_structure.inl
 * \brief In-Line subroutines of the <i>integration_structure.hpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.5
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

inline double CIntegration::GetCauchy_Value(void) { return Cauchy_Value; }

inline bool CIntegration::GetConvergence(void) { return Convergence; }

inline bool CIntegration::GetConvergence_FullMG(void) { return Convergence_FullMG; }

inline bool CIntegration::GetConvergence_OneShot(void) { return Convergence_OneShot; }

inline void CIntegration::SetConvergence(bool value) { Convergence = value; }

inline void CIntegration::SetMultiGrid_Solver(CGeometry ***geometry, CSolution ****solution_container, CNumerics *****solver_container, 
											  CConfig **config, unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone) { }
	
inline void CIntegration::Multigrid_Iteration(CGeometry ***geometry, CSolution ****solution_container, CNumerics *****solver_container,
							   CConfig **config, unsigned short iMesh, unsigned short mu, unsigned short RunTime_EqSystem,
							   unsigned long Iteration, unsigned short iZone) { }
										
inline void CIntegration::NonDimensional_Parameters(CGeometry **geometry, CSolution ***solution_container, CNumerics ****solver_container, 
																											CConfig *config, unsigned short FinestMesh, unsigned short RunTime_EqSystem, unsigned long Iteration, 
																											double *monitor) { }
	
inline void CIntegration::SetProlongated_Correction(CSolution *sol_fine, CGeometry *geo_fine, CConfig *config) { }

inline void CIntegration::SetProlongated_Solution(unsigned short RunTime_EqSystem, CSolution **sol_fine, CSolution **sol_coarse,
												  CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) { }

inline void CIntegration::SetRestricted_Residual(CSolution *sol_fine, CSolution *sol_coarse, CGeometry *geo_fine, 
											     CGeometry *geo_coarse, CConfig *config)  { }

inline void CIntegration::SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolution **sol_fine, CSolution **sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config)  { }

inline void CIntegration::SetRestricted_Gradient(unsigned short RunTime_EqSystem, CSolution **sol_fine, CSolution **sol_coarse, 
												 CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) { }
	
inline void CIntegration::SetResidual_Term(CGeometry *geometry, CSolution *flow) { }

inline void CIntegration::SetForcing_Term(CSolution *sol_fine, CSolution *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, 
										  CConfig *config) { }

inline void CIntegration::SetSingleGrid_Solver(CGeometry ***geometry, CSolution ****solution_container, CNumerics *****solver_container, 
											  CConfig **config, unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone) { }

inline void CIntegration::SetPotential_Solver(CGeometry ***geometry, CSolution ****solution_container, CNumerics *****solver_container, 
                                              CConfig **config, unsigned short RunTime_EqSystem, unsigned short iMesh, unsigned short iZone) { }
                                              
inline void CIntegration::Smooth_Solution(unsigned short RunTime_EqSystem, CSolution **solution, CGeometry *geometry, unsigned short val_nSmooth, double val_smooth_coeff, CConfig *config) { }

inline void CIntegration::Smooth_PrimVar(CSolution **solution, CGeometry *geometry, CConfig *config) { }

