/*!
 * \file CCGNSFileWriter.hpp
 * \brief Headers fo the tecplot ASCII writer class.
 * \author E. Saetta, L. Russo, R. Tognaccini
 * \version 7.0.8 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
#include "CFileWriter.hpp"
#include "../../solvers/CSolver.hpp"
#include "../COutputLegacy.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
class CGeometry;
class COutputLegacy;


class CCGNSFileWriter final: public CFileWriter{

  unsigned long nElem;  // Global number of elems without halos

  unsigned long timeIter;  //!< Current value of the time iteration
  su2double timeStep;       //!< Current value of the time step

public:

  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] valFileName - The name of the file
   * \param[in] valDataSorter - The parallel sorted data to write
   * \param[in] valTimeIter - The current time iteration
   * \param[in] valTimeStep - The current physical time step value
   */
  CCGNSFileWriter(string valFileName, CParallelDataSorter* valDataSorter);

  /*!
   * \brief Destructor
   */
  ~CCGNSFileWriter() override;

  /*!
   * \brief Write sorted data to file in tecplot ASCII file format
   */
//  void Write_Data_CGNS(CConfig *config);
   void Write_Data_CGNS(CConfig *config);

};

