/*!
 * \file CTecplotFileWriter.hpp
 * \brief Headers fo the tecplot ASCII writer class.
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

class CTecplotFileWriter final: public CFileWriter{

  unsigned long time_iter;  //!< Current value of the time iteration
  su2double timestep;       //!< Current value of the time step

public:

  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names, file extension and dimension.
   * \param[in] fields - A list of field names
   * \param[in] nDim - Physical dimension
   * \param[in] - The name of the file
   * \param[in] - The parallel sorted data to write
   */
  CTecplotFileWriter(vector<string> fields, unsigned short nDim,
                     string fileName, CParallelDataSorter *dataSorter,
                     unsigned long time_iter, su2double timestep);

  /*!
   * \brief Destructor
   */
  ~CTecplotFileWriter() override;

  /*!
   * \brief Write sorted data to file in tecplot ASCII file format
   */
  void Write_Data() override;

};

