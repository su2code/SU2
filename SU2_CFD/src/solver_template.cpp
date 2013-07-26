/*!
 * \file solution_template.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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

#include "../include/solver_structure.hpp"

CTemplateSolution::CTemplateSolution(void) : CSolution() { }

CTemplateSolution::CTemplateSolution(CGeometry *geometry, CConfig *config) : CSolution() { }

CTemplateSolution::~CTemplateSolution(void) { }

void CTemplateSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) { }

void CTemplateSolution::SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }

void CTemplateSolution::Centered_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *numerics,
																				 CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

void CTemplateSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *numerics,
																				CConfig *config, unsigned short iMesh) { }

void CTemplateSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *numerics, CNumerics *second_numerics,
																								 CConfig *config, unsigned short iMesh) { }

void CTemplateSolution::Source_Template(CGeometry *geometry, CSolution **solution_container, CNumerics *numerics,
																								 CConfig *config, unsigned short iMesh) { }

void CTemplateSolution::Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
																					 unsigned short iMesh) { }

void CTemplateSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *numerics, CConfig *config, 
																			unsigned short val_marker) { }

void CTemplateSolution::BC_HeatFlux_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

void CTemplateSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
																		 unsigned short val_marker) { }

void CTemplateSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, 
																 unsigned short val_marker) { }

void CTemplateSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, 
																	unsigned short val_marker) { }

void CTemplateSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, 
																		 unsigned short val_marker) { }

void CTemplateSolution::BC_Custom(CGeometry *geometry, CSolution **solution_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

void CTemplateSolution::ExplicitRK_Iteration(CGeometry *geometry, CSolution **solution_container, 
																						 CConfig *config, unsigned short iRKStep) { }

void CTemplateSolution::ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) { }

void CTemplateSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) { }
