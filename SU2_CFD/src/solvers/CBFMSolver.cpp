/*!
 * \file CBFMSolver.cpp
 * \brief Subroutines to be implemented for any new solvers
 * \author E.C. Bunschoten
 * \version 7.1.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/solvers/CBFMSolver.hpp"
#include "../../include/variables/CBFMVariable.hpp"
#include "../../include/numerics/ReadBFMInput.hpp"

CBFMSolver::CBFMSolver(void) : CSolver() { }

CBFMSolver::CBFMSolver(CGeometry *geometry, CConfig *config) : CSolver() {
    nodes = new CBFMVariable();
    nPoint = geometry->GetnPoint();
    nDim = geometry->GetnDim();
    
}

CBFMSolver::~CBFMSolver(void) { 
    delete nodes;

}

void CBFMSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) { }

void CBFMSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }

void CBFMSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                         CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

void CBFMSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                        CConfig *config, unsigned short iMesh) { }

void CBFMSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                      CNumerics **numerics_container, CConfig *config, unsigned short iMesh) 
{
    

}

void CBFMSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                                 CConfig *config, unsigned short iMesh) { }

void CBFMSolver::BC_Euler_Wall(CGeometry      *geometry,
                                    CSolver        **solver_container,
                                    CNumerics      *conv_numerics,
                                    CNumerics      *visc_numerics,
                                    CConfig        *config,
                                    unsigned short val_marker) { }

void CBFMSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

void CBFMSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                     unsigned short val_marker) { }

void CBFMSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                 unsigned short val_marker) { }

void CBFMSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                  unsigned short val_marker) { }

void CBFMSolver::BC_Sym_Plane(CGeometry      *geometry,
                                   CSolver        **solver_container,
                                   CNumerics      *conv_numerics,
                                   CNumerics      *visc_numerics,
                                   CConfig        *config,
                                   unsigned short val_marker) { }

void CBFMSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

void CBFMSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                             CConfig *config, unsigned short iRKStep) { }

void CBFMSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

void CBFMSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }
