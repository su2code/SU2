/*!
 * \file iteration_structure.hpp
 * \brief Headers of the main subroutines used by SU2_CFD.
 *        The subroutines and functions are in the <i>definition_structure.cpp</i> file.
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

#ifndef NO_MPI
#include <mpi.h>
#endif
#include <ctime>

#include "solver_structure.hpp"
#include "integration_structure.hpp"
#include "output_structure.hpp"
#include "numerics_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../Common/include/config_structure.hpp"

using namespace std;

/*! 
 * \brief ________________________.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxs of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void MeanFlowIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
											 CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, 
											 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*! 
 * \brief ________________________.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxs of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void PlasmaIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
										 CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, 
										 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*! 
 * \brief ________________________.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxs of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void FreeSurfaceIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
													CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, 
													CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*! 
 * \brief ________________________.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxs of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void FluidStructureIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
														 CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, 
														 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*! 
 * \brief ________________________.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxs of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void AeroacousticIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
													 CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, 
													 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*! 
 * \brief ________________________.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxs of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void WaveIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
									 CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, 
									 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*! 
 * \brief ________________________.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxs of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void FEAIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
									CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, 
									CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*! 
 * \brief ________________________.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxs of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void AdjMeanFlowIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
													CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, 
													CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*! 
 * \brief ________________________.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxs of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void AdjPlasmaIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
												CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, 
												CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*! 
 * \brief ________________________.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxs of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void AdjFreeSurfaceIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
														 CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, 
														 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*! 
 * \brief ________________________.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxs of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void AdjAeroacousticIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
															CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, 
															CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*! 
 * \brief Updates the positions and grid velocities for dynamic meshes between physical time steps.
 * \author T. Economon
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxs of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iZone - Index of the zone.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void SetGrid_Movement(CGeometry **geometry_container, CSurfaceMovement *surface_movement, 
                      CVolumetricMovement *grid_movement, CFreeFormDefBox **FFDBox,
                      CSolver ***solver_container, CConfig *config_container, unsigned short iZone, unsigned long ExtIter);

/*!
 * \brief Computation and storage of the time spectral source terms.
 * \author T. Economon, K. Naik
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] config - Definition of the particular problem.
 * \param[in] nZone - Total number of zones (periodic instances).
 * \param[in] iZone - Current zone number.
 */
void SetTimeSpectral(CGeometry ***geometry_container, CSolver ****solver_container,
		CConfig **config_container, unsigned short nZone, unsigned short iZone);

/*!
 * \brief Computation and storage of the time-spectral mesh velocities.
 * \author T. Economon, K. Naik
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] nZone - Total number of zones (periodic instances).
 */
void SetTimeSpectral_Velocities(CGeometry ***geometry_container,
		CConfig **config_container, unsigned short nZone);

/*!
 * \brief Search and interpolate routine for matching two zones in relative motion (sliding).
 * \author T. Economon
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] config - Definition of the particular problem.
 * \param[in] nZone - Total number of zones (periodic instances).
 */
void SetSliding_Interfaces(CGeometry ***geometry_container, CSolver ****solver_container,
                     CConfig **config_container, unsigned short nZone);
