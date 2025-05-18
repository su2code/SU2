/*!
 * \file CParaviewFileWriter.hpp
 * \brief Headers fo the paraview file writer class.
 * \author T. Albring
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

class CParaviewFileWriter final: public CFileWriter{

public:

  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] valDataSorter - The parallel sorted data to write
   */
  CParaviewFileWriter(CParallelDataSorter* valDataSorter);

  /*!
   * \brief Destructor
   */
  ~CParaviewFileWriter() override;

  /*!
   * \brief Write sorted data to file in paraview file format
   * \param[in] val_filename - The name of the file
   */
  void WriteData(string val_filename) override ;

};

