/*!
 * \file CCSVFileWriter.hpp
 * \brief Headers fo the CSV file writer class.
 * \author T. Albring
 * \version 6.2.0 "Falcon"
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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
#include "CFileWriter.hpp"

class CCSVFileWriter final: public CFileWriter{
  
public:
  
  /*!
   * \brief Construct a file writer using field names, file extension and dimension.
   * \param[in] fields - A list of field names
   * \param[in] nDim - Physical dimension
   */  
  CCSVFileWriter(vector<string> fields, unsigned short nDim, string fileName, CParallelDataSorter* data_sorter);
  
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

