/*!
 * \file CCSVFileWriter.hpp
 * \brief Headers fo the CSV file writer class.
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

class CSU2MeshFileWriter final: public CFileWriter{

private:
  unsigned short iZone, //!< Index of the current zone
  nZone;                //!< Number of zones

public:

  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names, file extension and dimension.
   * \param[in] fields - A list of field names
   * \param[in] nDim - Physical dimension
   * \param[in] fileName - The name of the file
   * \param[in] data_sorter - The parallel sorted data to write
   * \param[in] iZone - Index of the current zone
   * \param[in] nZone - Number of zones
   */
  CSU2MeshFileWriter(vector<string> fields, unsigned short nDim,
                     string fileName, CParallelDataSorter* data_sorter,
                     unsigned short iZone, unsigned short nZone);

  /*!
   * \brief Destructor
   */
  ~CSU2MeshFileWriter() override;

  /*!
   * \brief Write sorted data to file in SU2 mesh file format
   * \param[in] - The name of the file
   * \param[in] - The parallel sorted data to write
   */
  void Write_Data() override;

};

