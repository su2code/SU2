/*!
 * \file CFileWriter.hpp
 * \brief Headers fo the file writer class.
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
#include "../../../../Common/include/mpi_structure.hpp"
#include "../../../../Common/include/option_structure.hpp"
#include <sys/stat.h>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>

#include "../../output/filewriter/CParallelDataSorter.hpp"

using namespace std;

class CFileWriter{
protected:

  /*!
   * \brief Vector containing the field names
   */
  std::vector<std::string> fieldnames;
  
  /*!
   * \brief The physical dimension of the problem
   */
  unsigned short nDim;
  
  /*!
   * \brief The MPI rank
   */
  int rank;
  
  /*!
   * \brief The MPI size, aka the number of processors.
   */
  int size;
  
  /*!
   * \brief The file extension to be attached to the filename.
   */
  std::string file_ext;
  
  su2double StartTime, /*!< \brief Start time for time measurement  */
            StopTime , /*!< \brief Stop time for time measurement  */
            UsedTime , /*!< \brief Used time for time measurement  */
            Bandwidth, /*!< \brief Used bandwidth  */
            file_size; /*!< \brief Size of the last written file  */
  
  /*!
   * \brief Determine the file size
   * \param[in] filename - Name of the file.
   * \return 
   */
  inline unsigned long Determine_Filesize(std::string filename){
      struct stat stat_buf;
      int rc = stat(filename.c_str(), &stat_buf);
      return rc == 0 ? stat_buf.st_size : -1;
  }
  
  /*!
   * \brief Filename
   */
  string fileName;
  
  /*!
   * \brief The parallel data sorter
   */
  CParallelDataSorter* dataSorter;
  
public:
  /*!
   * \brief Construct a file writer using field names, file extension and dimension.
   * \param[in] fields - A list of field names
   * \param[in] file_ext - The file extension to be attached to the filename
   * \param[in] nDim - Physical dimension
   */  
  CFileWriter(std::vector<std::string> fields, string fileName, CParallelDataSorter* dataSorter, string file_ext, unsigned short nDim);
  
  /*!
   * \brief Destructor
   */
  virtual ~CFileWriter();
  
  /*!
   * \brief Write sorted data to file
   */
  virtual void Write_Data(){}
  
  /*!
   * \brief Get the bandwith used for the last writing
   */
  su2double Get_Bandwidth(){return Bandwidth;}
  
  /*!
   * \brief Get the filesize of the last written file.
   */
  su2double Get_Filesize(){return file_size;}
  
  /*!
   * \brief Get the used time of the last file writing.
   * \return 
   */
  su2double Get_UsedTime(){return UsedTime;}

};

