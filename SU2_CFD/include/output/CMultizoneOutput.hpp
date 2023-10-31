/*!
 * \file CMultizoneOutput.hpp
 * \brief Headers of the main subroutines for screen and history output in multizone problems.
 * \author R. Sanchez, T. Albring
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#ifdef HAVE_CGNS
  #include "cgnslib.h"
#endif
#ifdef HAVE_TECIO
  #include "TECIO.h"
#endif
#include <fstream>
#include <cmath>

#include "COutput.hpp"
#include "../../../Common/include/CConfig.hpp"

using namespace std;

/*!
 * \class CMultizoneOutput
 * \brief Class for writing the multizone output.
 * \author R. Sanchez, T. Albring.
 */
class CMultizoneOutput final: public COutput {

protected:
  unsigned short nZone; //!< Number of zones

  string bgs_res_name; //!< Block-Gauss Seidel residual name
  bool write_zone;     //!< Boolean indicating whether the individual zones write to screen

public:

  /*!
   * \brief Constructor of the class.
   */
  CMultizoneOutput(const CConfig *driver_config, const CConfig* const* config, unsigned short nDim);

  /*!
   * \brief Load the multizone history output field values
   * \param[in] output - Container holding the output instances per zone.
   * \param[in] config - Definition of the particular problem.
   */
  void LoadMultizoneHistoryData(const COutput* const* output, const CConfig* const* config) override;

  /*!
   * \brief Set the available multizone history output fields
   * \param[in] output - Container holding the output instances per zone.
   */
  void SetMultizoneHistoryOutputFields(const COutput* const* output, const CConfig* const* config) override;

  /*!
   * \brief Determines if the history file output.
   * \param[in] config - Definition of the particular problem.
   */
  bool WriteHistoryFileOutput(const CConfig *config) override ;

  /*!
   * \brief Determines if the screen header should be written.
   * \param[in] config - Definition of the particular problem.
   */
  bool WriteScreenHeader(const CConfig *config) override ;

  /*!
   * \brief Determines if the screen header should be written.
   * \param[in] config - Definition of the particular problem.
   */
  bool WriteScreenOutput(const CConfig *config) override ;
};
