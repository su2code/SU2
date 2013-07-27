/*!
 * \file solution_template.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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

#include "../include/solver_structure.hpp"

CTemplateSolver::CTemplateSolver(void) : CSolver() { }

CTemplateSolver::CTemplateSolver(CGeometry *geometry, CConfig *config) : CSolver() { }

CTemplateSolver::~CTemplateSolver(void) { }

void CTemplateSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) { }

void CTemplateSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }

void CTemplateSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
																				 CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

void CTemplateSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
																				CConfig *config, unsigned short iMesh) { }

void CTemplateSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
																								 CConfig *config, unsigned short iMesh) { }

void CTemplateSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
																								 CConfig *config, unsigned short iMesh) { }

void CTemplateSolver::Solve_LinearSystem(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
																					 unsigned short iMesh) { }

void CTemplateSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
																			unsigned short val_marker) { }

void CTemplateSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

void CTemplateSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
																		 unsigned short val_marker) { }

void CTemplateSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, 
																 unsigned short val_marker) { }

void CTemplateSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, 
																	unsigned short val_marker) { }

void CTemplateSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, 
																		 unsigned short val_marker) { }

void CTemplateSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

void CTemplateSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, 
																						 CConfig *config, unsigned short iRKStep) { }

void CTemplateSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

void CTemplateSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }
