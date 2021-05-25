/*!
 * \file CDeformationDriver.hpp
 * \brief Headers of the main subroutines for driving the mesh deformation.
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

#include "../../../Common/include/grid_movement/CSurfaceMovement.hpp"
#include "../../../Common/include/grid_movement/CVolumetricMovement.hpp"
#include "../../../SU2_CFD/include/output/COutput.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"

/*!
 * \class CDeformationDriver
 * \brief Class for driving mesh deformation solvers.
 * \author A. Gastaldi, H. Patel
 * \version 7.1.1 "Blackbird"
 */
class CDeformationDriver {
protected:
  char config_file_name[MAX_STRING_SIZE];
  int rank,
      size;
  su2double StartTime,                          /*!< \brief Start point of the timer for performance benchmarking.*/
            StopTime,                           /*!< \brief Stop point of the timer for performance benchmarking.*/
            UsedTimePreproc,                    /*!< \brief Elapsed time between Start and Stop point of the timer for tracking preprocessing phase.*/
            UsedTimeCompute,                    /*!< \brief Elapsed time between Start and Stop point of the timer for tracking compute phase.*/
            UsedTime;                           /*!< \brief Elapsed time between Start and Stop point of the timer.*/
  unsigned short iZone, nZone = SINGLE_ZONE;
  CConfig *driver_config;                       /*!< \brief Definition of the driver configuration. */
  CConfig **config_container;                   /*!< \brief Definition of the particular problem. */
  CGeometry **geometry_container;             /*!< \brief Geometrical definition of the problem. */
  CSurfaceMovement **surface_movement;          /*!< \brief Surface movement classes of the problem. */
  CVolumetricMovement **grid_movement;         /*!< \brief Volume grid movement classes of the problem. */
  COutput **output_container;                   /*!< \brief Pointer to the COutput class. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CDeformationDriver(char* confFile, SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CDeformationDriver(void);

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
  void Geometrical_Preprocessing(CConfig *config, CGeometry *&geometry);

  /*!
   * \brief Preprocess the output container.
   */
  void Output_Preprocessing(CConfig *config, CGeometry *geometry, COutput *&output);

};
