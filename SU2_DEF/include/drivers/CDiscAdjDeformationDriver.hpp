/*!
* \file CDiscAdjDeformationDriver.cpp
* \brief Main subroutines for driving the projection of sensitivities.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../Common/include/parallelization/mpi_structure.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/fem/fem_geometry_structure.hpp"
#include "../../../Common/include/grid_movement/CSurfaceMovement.hpp"
#include "../../../Common/include/grid_movement/CVolumetricMovement.hpp"
#include "../../../SU2_CFD/include/output/COutput.hpp"
#include "../../../SU2_CFD/include/output/CBaselineOutput.hpp"
#include "../../../SU2_CFD/include/solvers/CBaselineSolver.hpp"

/*!
 * \class CDiscAdjDeformationDriver
 * \brief Class for driving sensitivity DiscAdjDeformations.
 * \author A. Gastaldi, H. Patel
 * \version 7.1.1 "Blackbird"
 */
class CDiscAdjDeformationDriver {
protected:
  char config_file_name[MAX_STRING_SIZE];
  int rank,
      size;
  unsigned short iZone, nZone = SINGLE_ZONE;
  unsigned short iInst;
  unsigned short* nInst;
  su2double StartTime,                          /*!< \brief Start point of the timer for performance benchmarking.*/
            StopTime,                           /*!< \brief Stop point of the timer for performance benchmarking.*/
            UsedTimePreproc,                    /*!< \brief Elapsed time between Start and Stop point of the timer for tracking preprocessing phase.*/
            UsedTimeCompute,                    /*!< \brief Elapsed time between Start and Stop point of the timer for tracking compute phase.*/
            UsedTime;                           /*!< \brief Elapsed time between Start and Stop point of the timer.*/
  su2double** Gradient;
  ofstream Gradient_file;
  CConfig *driver_config;                       /*!< \brief Definition of the driver configuration. */
  CConfig **config_container;                   /*!< \brief Definition of the particular problem. */
  CGeometry ***geometry_container;             /*!< \brief Geometrical definition of the problem. */
  CSurfaceMovement **surface_movement;
  CVolumetricMovement **grid_movement;
  COutput **output_container;                   /*!< \brief Pointer to the COutput class. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CDiscAdjDeformationDriver(char* confFile, SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjDeformationDriver(void);

  /*!
   * \brief [Overload] Launch the computation for single-zone problems.
   */
  void Run();

  /*!
   * \brief Deallocation routine
   */
  void Postprocessing();

protected:
  /*!
   * \brief Init_Containers
   */
  void SetContainers_Null();

  /*!
   * \brief Read in the config and mesh files.
   */
  void Input_Preprocessing();

  /*!
   * \brief Construction of the edge-based data structure.
   */
  void Geometrical_Preprocessing();

  /*!
   * \brief Preprocess the output container.
   */
  void Output_Preprocessing();

  /*!
   * \brief DiscAdjDeformation of the surface sensitivity using finite differences (FD).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement class of the problem.
   * \param[in] Gradient_file - Output file to store the gradient data.
   */

  void SetDiscAdjDeformation_FD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double **Gradient);

  /*!
   * \brief DiscAdjDeformation of the surface sensitivity using algorithmic differentiation (AD).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement class of the problem.
   * \param[in] Gradient_file - Output file to store the gradient data.
   */

  void SetDiscAdjDeformation_AD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double **Gradient);

  /*!
   * \brief Prints the gradient information to a file.
   * \param[in] Gradient - The gradient data.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Gradient_file - Output file to store the gradient data.
   */

  void OutputGradient(su2double** Gradient, CConfig* config, ofstream& Gradient_file);

  /*!
   * \brief Write the sensitivity (including mesh sensitivity) computed with the discrete adjoint method
   *  on the surface and in the volume to a file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - Number of Zones.
   */

  void SetSensitivity_Files(CGeometry ***geometry, CConfig **config, unsigned short val_nZone);

};
