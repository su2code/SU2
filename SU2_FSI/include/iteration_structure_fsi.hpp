/*!
 * \file iteration_structure_fsi.hpp
 * \brief Headers of the main subroutines used by SU2_CFD.
 *        The subroutines and functions are in the <i>definition_structure.cpp</i> file.
 * \author R. Sanchez
 * \version 3.2.9 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
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

#include "../../SU2_CFD/include/solver_structure.hpp"
#include "../../SU2_CFD/include/integration_structure.hpp"
#include "../../SU2_CFD/include/output_structure.hpp"
#include "../../SU2_CFD/include/numerics_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../Common/include/config_structure.hpp"

using namespace std;

/*!
 * \brief Block Gauss-Seidel Iteration function for Fluid-Structure Interaction applications.
 * \author R. Sanchez, F. Palacios.
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
 * \param[in] nFluidIt - Number of fluid iterations within a fixed time step.
 */
void FSI_BGS_Iteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
						  CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
						  CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
						  unsigned long iFluidIt, unsigned long nFluidIt);

/*!
 * \brief CFD Subiteration function for Fluid-Structure Interaction applications.
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
void Flow_Subiteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
						  CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
						  CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*!
 * \brief CFD update function for Fluid-Structure Interaction applications.
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
void Flow_Update(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
				   CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
				   CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
				   unsigned long ExtIter);


/*!
 * \brief FEA Subiteration function for Fluid-Structure Interaction applications (legacy).
 * \author R. Sanchez.
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
void FEA_Subiteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
					     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
						 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*!
 * \brief FEM Subiteration function for Fluid-Structure Interaction applications (structural side).
 * \author R. Sanchez.
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
void FEM_Subiteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
					     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
						 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);


/*!
 * \brief Displacement transfer function for Fluid-Structure Interaction applications.
 * \author R. Sanchez.
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
void FSI_Disp_Transfer(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
					     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
						 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*!
 * \brief Load transfer function for Fluid-Structure Interaction applications.
 * \author R. Sanchez.
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
void FSI_Load_Transfer(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
					     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
						 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
						 unsigned long ExtIter);


/*!
 * \brief FEA update function for Fluid-Structure Interaction applications (legacy).
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void FEA_Update(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
				   CSolver ****solver_container, CConfig **config_container, unsigned long ExtIter);

/*!
 * \brief FEM update function for Fluid-Structure Interaction applications.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void FEM_Update(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
				   CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, unsigned long ExtIter);


/*!
 * \brief Relaxation step for displacement transfer.
 * \author R. Sanchez.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] iFSIIter - Current FSI iteration number.
 */
void FSI_Disp_Relaxation(COutput *output, CGeometry ***geometry_container, CSolver ****solver_container,
					     CConfig **config_container, unsigned long iFSIIter);

/*!
 * \brief Relaxation step for load transfer.
 * \author R. Sanchez.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] iFSIIter - Current FSI iteration number.
 */
void FSI_Load_Relaxation(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
	     	 	 	 	 	CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
	     	 	 	 	 	CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*!
 * \brief Displacement predictor function for Fluid-Structure Interaction applications.
 * \author R. Sanchez.
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
void FSI_Disp_Predictor(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
					     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
						 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);

/*!
 * \brief Load predictor function for Fluid-Structure Interaction applications.
 * \author R. Sanchez.
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
void FSI_Load_Predictor(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
					     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
						 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
						 unsigned long ExtIter);
