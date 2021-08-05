/*!
 * \file COneShotSolver.hpp
 * \brief Main header for defining the OneShot solver.
 * \author L. Kusch, T.Dick
 * \version 7.1.1 "Blackbird"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "CDiscAdjSolver.hpp"

class COneShotSolver : public CDiscAdjSolver {
private:

  // add variables as needed later!

public:

  /*!
   * \brief Constructor of the class.
   */
  COneShotSolver(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Initialize the discrete adjoint solver with the corresponding direct solver.
   * \param[in] Kind_Solver - The kind of direct solver.
   */
  COneShotSolver(CGeometry *geometry, CConfig *config, CSolver* solver, unsigned short Kind_Solver, unsigned short iMesh);

  ~COneShotSolver(void);

  /*!
   * \brief Prepare the solver for a new recording (without setting solution to initial solution).
   * \param[in] geometry - geometry class object
   * \param[in] config - config class object
   */
  void SetRecording(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Store the coordinates for the optimization process.
   * \param[in] geometry - geometry class object
   * \param[in] config - config class object
   */
  void StoreMeshPoints(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Load the stored coordinates for resetting the optimization process.
   * \param[in] geometry - geometry class object
   * \param[in] config - config class object
   */
  void LoadMeshPoints(CGeometry *geometry, CConfig *config);

  /*!
   * \brief After resetting the mesh coordinates do the post processing from mesh deformation.
   * \param[in] geometry_container - we need the geometry on all multigrid mesh levels
   * \param[in] grid_movement - need the class that was used for the original deformation
   * \param[in] config - config class object
   */
  void UpdateAuxiliaryGeometryVariables(CGeometry **geometry_container, CVolumetricMovement *grid_movement, CConfig *config);

  /*!
   * \brief Evaluate a geometric property for optimization.
   */
  su2double EvaluateGeometryFunction(CGeometry* geometry, CConfig *config, unsigned int iPlane);

  /*!
   * \brief Evaluate the gradient of a geometric property for optimization.
   */
  vector<su2double> EvaluateGeometryGradient(CGeometry* geometry, CSurfaceMovement* surface_movement, CConfig* config);

};
