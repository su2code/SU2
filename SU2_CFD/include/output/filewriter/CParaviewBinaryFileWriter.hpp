/*!
 * \file CParaviewBinaryFileWriter.hpp
 * \brief Headers fo paraview binary file writer class.
 * \author T. Albring
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

class CParaviewBinaryFileWriter final: public CFileWriter{

public:

  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names, dimension.
   * \param[in] fields - A list of field names
   * \param[in] nDim - Physical dimension
   * \param[in] fileName - The name of the file
   * \param[in] data_sorter - The parallel sorted data to write
   */
  CParaviewBinaryFileWriter(vector<string> fields, unsigned short nDim,
                            string fileName, CParallelDataSorter* data_sorter);

  /*!
   * \brief Destructor
   */
  ~CParaviewBinaryFileWriter() override;

  /*!
   * \brief Write sorted data to file in paraview binary file format
   */
  void Write_Data() override;

private:

  /*!
   * \brief Change storage of buffer from big endian to little endian
   * \param buffer - Pointer to the beginning of the buffer
   * \param nBytes - The size in bytes of an data entry
   * \param nVar - The number of entries
   */
  void SwapBytes(char *buffer, size_t nBytes, unsigned long nVar);

};

