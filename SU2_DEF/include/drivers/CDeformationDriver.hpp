/*!
 * \file CDeformationDriver.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
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

/*!
 * \class CDeformationDriver
 * \brief Class for driving single-zone solvers.
 * \author R. Sanchez
 * \version 7.1.1 "Blackbird"
 */
class CDeformationDriver {
protected:
  char config_file_name[MAX_STRING_SIZE];

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CDeformationDriver(char* confFile,
             SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CDeformationDriver(void);

  /*!
   * \brief [Overload] Launch the computation for single-zone problems.
   */
  void RunSolver();

protected:
    /*!
     * \brief Init_Containers
     */
    void SetContainers_Null();

    /*!
     * \brief Read in the config and mesh files.
     */
    void Input_Preprocessing(CConfig **&config, CConfig *&driver_config);

};
