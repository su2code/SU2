/*!
 * \file output_driver.hpp
 * \brief Headers of the main subroutines for screen and history output in multizone problems.
 * \author R. Sanchez, T. Albring
 * \version 6.1.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include "../../../Common/include/mpi_structure.hpp"

#ifdef HAVE_CGNS
  #include "cgnslib.h"
#endif
#ifdef HAVE_TECIO
  #include "TECIO.h"
#endif
#include <fstream>
#include <cmath>
#include <time.h>
#include <fstream>

#include "COutput.hpp"
#include "../../../Common/include/config_structure.hpp"

using namespace std;

/*!
 * \class CDriverOutput
 * \brief Class for writing the multizone output.
 * \author R. Sanchez, T. Albring.
 */
class CDriverOutput : public COutput {

protected:
  unsigned short nZone;
  
  string bgs_res_name;
  bool write_zone;
  
public:

  /*!
   * \brief Constructor of the class.
   */
  CDriverOutput(CConfig *driver_config, CConfig** config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CDriverOutput(void);

  /*!
   * \brief Load the history file fields for multizone computations.
   * \param[in] output - Output information of the particular multizone problem.
   * \param[in] config - Definition of the particular multizone problem.
   */
  void LoadMultizoneHistoryData(COutput **output, CConfig **config);

  /*!
   * \brief Write the history output fields for multizone computations.
   * \param[in] output - Output information of the particular multizone problem.
   * \param[in] config - Definition of the particular multizone problem.
   */
  void SetMultizoneHistoryOutputFields(COutput **output, CConfig **config);

  /*! 
   * \brief Determines if the screen header should be written.
   * \param[in] config - Definition of the particular problem.
   */
  bool WriteScreen_Header(CConfig *config);

  /*!
   * \brief Determines if the screen header should be written.
   * \param[in] config - Definition of the particular problem.
   */
  bool WriteScreen_Output(CConfig *config);
  
  /*!
   * \brief Determines if the history file output.
   * \param[in] config - Definition of the particular problem.
   */
  bool WriteHistoryFile_Output(CConfig *config);

};
