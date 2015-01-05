/*!
 * \file iteration_structure.hpp
 * \brief Headers of the main subroutines used by SU2_CFD.
 *        The subroutines and functions are in the <i>definition_structure.cpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 3.2.7 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (fpalacios@stanford.edu).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
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

#ifdef HAVE_MPI
  #include "mpi.h"
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
 * \param[in] FFDBox - FFD FFDBoxes of the problem.
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
 * \param[in] FFDBox - FFD FFDBoxes of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void TNE2Iteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
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
 * \param[in] FFDBox - FFD FFDBoxes of the problem.
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
 * \param[in] FFDBox - FFD FFDBoxes of the problem.
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
 * \param[in] FFDBox - FFD FFDBoxes of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void HeatIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
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
 * \param[in] FFDBox - FFD FFDBoxes of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void PoissonIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
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
 * \param[in] FFDBox - FFD FFDBoxes of the problem.
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
 * \param[in] FFDBox - FFD FFDBoxes of the problem.
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
 * \param[in] FFDBox - FFD FFDBoxes of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void AdjTNE2Iteration(COutput *output, CIntegration ***integration_container,
                      CGeometry ***geometry_container,
                      CSolver ****solver_container,
                      CNumerics *****numerics_container,
                      CConfig **config_container,
                      CSurfaceMovement **surface_movement,
                      CVolumetricMovement **grid_movement,
                      CFreeFormDefBox*** FFDBox);

/*!
 * \brief Imposes a gust via the grid velocities.
 * \author S. Padron
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 */
void SetWind_GustField(CConfig *config_container, CGeometry **geometry_container, CSolver ***solver_container);

/*!
 * \brief Reads and initializes the vortex positions, strengths and gradient.
 * \author S. Padron
 * \param[in] nVortex - number of vortices.
 * \param[in] x0 - Vector of x-loc of the vortices.
 * \param[in] y0 - Vector of y-loc of the vortices.
 * \param[in] vort_strength - Vector of vortex strengths.
 * \param[in] r_core - Vector of vortex core size.
 */
void InitializeVortexDistribution(unsigned long &nVortex, vector<double>& x0,vector<double>& y0,vector<double>& vort_strength,vector<double>& r_core);

/*!
 * \brief Updates the positions and grid velocities for dynamic meshes between physical time steps.
 * \author T. Economon
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxes of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iZone - Index of the zone.
 * \param[in] IntIter - Current sudo time iteration number.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void SetGrid_Movement(CGeometry **geometry_container, CSurfaceMovement *surface_movement, 
                      CVolumetricMovement *grid_movement, CFreeFormDefBox **FFDBox,
                      CSolver ***solver_container, CConfig *config_container, unsigned short iZone, unsigned long IntIter, unsigned long ExtIter);

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
 * \brief Computation of the Time-Spectral operator matrix.
 * \author K. Naik
 * \param[in] D - double pointer to the operator matrix.
 * \param[in] nZone - Total number of zones (periodic instances).
 */
void ComputeTimeSpectral_Operator(double **D, double period, unsigned short nZone);

/*!
 * \brief Computation and storage of the time-spectral mesh velocities.
 * \author K. Naik, T. Economon
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] nZone - Total number of zones (periodic instances).
 */
void SetTimeSpectral_Velocities(CGeometry ***geometry_container,
		CConfig **config_container, unsigned short nZone);
