/*!
 * \file CCSVFileWriter.hpp
 * \brief Headers fo the CSV file writer class.
 * \author T. Albring
 * \version 7.0.4 "Blackbird"
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

class CCSVFileWriter final: public CFileWriter{

public:

  /*!
   * \brief Construct a file writer using field name and the data sorter.
   * \param[in] valFileName - The name of the file
   * \param[in] valDataSorter - The parallel sorted data to write
   */
  CCSVFileWriter(string valFileName, CParallelDataSorter* valDataSorter);

  /*!
   * \brief Destructor
   */
  ~CCSVFileWriter() override;

  /*!
   * \brief Write sorted data to file in CSV file format
   * \param[in] - The name of the file
   * \param[in] - The parallel sorted data to write
   */
  void Write_Data() override;

};

