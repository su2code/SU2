/*!
 * \file COneShotSinglezoneDriver.hpp
 * \brief Header to define a driver for OneShot problems.
 *        Logic is based on the CAdjSinglezoneDriver class.
 * \author T.Dick
 * \version 7.0.6 "Blackbird"
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

#include "CDiscAdjSinglezoneDriver.hpp"

/*!
 * \class COneShotSinglezoneDriver
 * \brief Class that extends the AdjSinglezoneDriver class with methods for one-shot optimization
 * \author L. Kusch, T.Dick
 * \version 7.0.6 "Blackbird"
 */
class COneShotSinglezoneDriver : public CDiscAdjSinglezoneDriver {

protected:

  unsigned long nPiggyIter;           /*!< \brief The number of coupled primal and adjoint iterations that are run on the PiggyBack solver.*/
  unsigned short nDV_Total;           /*!< \brief Total number of design variables used in optimization.*/


public:

  /*!
    * \brief Constructor of the class.
    * \param[in] confFile - Configuration file name.
    * \param[in] val_nZone - Total number of zones.
    */
  COneShotSinglezoneDriver(char* confFile,
                   unsigned short val_nZone,
                   SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~COneShotSinglezoneDriver(void);

  /*!
   * \brief Preprocess the one-shot iteration
   * \param[in] TimeIter - index of the current time-step.
   */
  void Preprocess(unsigned long TimeIter) override;

  /*!
   * \brief Runs main routines of the one-shot class.
   */
  void Run(void) override;

  /*!
   * \brief Postprocess the One-Shot driver.
   * \note no secondary recording is needed here like in the DiscAdj driver!
   */
  void Postprocess(void) override;

  /*!
   * \brief Runs an optimization using the one-shot method.
   */
  void RunOneShot();

  /*!
   * \brief Run a Piggyback iteration
   */
  void PiggyBack();

  /*!
   * \brief Executes one primal and dual iteration (in a piggy-back manner).
   * \param[in] current Piggyback iteration number
   */
  void PrimalDualStep(unsigned long iPiggyIter);

  /*!
   * \brief Record one iteration of a flow iteration in within multiple zones.
   * \param[in] kind_recording - Type of recording (either CONS_VARS, MESH_COORDS, COMBINED or NONE)
   */
  void SetRecording(unsigned short kind_recording);

  /*!
   * \brief Projection of the surface sensitivity using algorithmic differentiation (AD) (see also SU2_DOT).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement class of the problem.
   * \param[in] Gradient - Output to store the gradient data.
   */
  void SetProjection_AD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double* Gradient);

  /*!
   * \brief Performs a surface deformation and volumetric deformation (see also SU2_DEF) for a design update.
   * \param[in] geometry - geometry class.
   * \param[in] config - config class.
   * \param[in] surface_movement - surface movement class
   * \param[in] grid_movement - volumetric movement class
   */
  void SurfaceDeformation(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, CVolumetricMovement *grid_movement);

  /*!
   * \brief Compute the search direction using the gradients and the preconditioner.
   */
  void ComputeSearchDirection();

  /*!
   * \brief Check if the search direction is a descent direction.
   */
  bool CheckDescent();

  /*!
   * \brief Store values for the updated design variables.
   */
  void UpdateDesignVariable();

  /*!
   * \brief Project mesh sensitivities to surface and design variables using AD.
   */
  void ProjectMeshSensitivities();

  /*!
   * \brief Compute the inverse preconditioner matrix (BCheck^(-1)) for the multiplier update.
   */
  void ComputePreconditioner();
};
