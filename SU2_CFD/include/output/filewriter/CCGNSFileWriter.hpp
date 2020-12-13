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

class CCGNSFileWriter final: public CFileWriter{

  const CConfig* config;
  unsigned long nElem;  // Global number of elems without halos

public:

  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] valFileName - The name of the file
   * \param[in] valDataSorter - The parallel sorted data to write
   */
  CCGNSFileWriter(string valFileName, CParallelDataSorter* valDataSorter);

  /*!
   * \brief Destructor
   */
  ~CCGNSFileWriter() override;

  /*!
   * \brief Write sorted data to file in CGNS file format
   */
//  void Write_Data();
   void Write_Data();

   void SetConfig(CConfig* configClass) {config = configClass;}

   void callCGNS(int error){if (error) cg_error_print();}

};

