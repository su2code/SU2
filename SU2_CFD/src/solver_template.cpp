/*!
 * \file solution_template.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author F. Palacios
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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


#include "../include/solver_structure.hpp"

CTemplateSolver::CTemplateSolver(void) : CSolver() { }

CTemplateSolver::CTemplateSolver(CGeometry *geometry, CConfig *config) : CSolver() { }

CTemplateSolver::~CTemplateSolver(void) { }

void CTemplateSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) { }

void CTemplateSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }

void CTemplateSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

void CTemplateSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                        CConfig *config, unsigned short iMesh) { }

void CTemplateSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                                 CConfig *config, unsigned short iMesh) { }

void CTemplateSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                                 CConfig *config, unsigned short iMesh) { }

void CTemplateSolver::BC_Euler_Wall(CGeometry      *geometry,
                                    CSolver        **solver_container,
                                    CNumerics      *conv_numerics,
                                    CNumerics      *visc_numerics,
                                    CConfig        *config,
                                    unsigned short val_marker) { }

void CTemplateSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

void CTemplateSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                     unsigned short val_marker) { }

void CTemplateSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, 
                                 unsigned short val_marker) { }

void CTemplateSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, 
                                  unsigned short val_marker) { }

void CTemplateSolver::BC_Sym_Plane(CGeometry      *geometry,
                                   CSolver        **solver_container,
                                   CNumerics      *conv_numerics,
                                   CNumerics      *visc_numerics,
                                   CConfig        *config,
                                   unsigned short val_marker) { }

void CTemplateSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

void CTemplateSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, 
                                             CConfig *config, unsigned short iRKStep) { }

void CTemplateSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

void CTemplateSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }
