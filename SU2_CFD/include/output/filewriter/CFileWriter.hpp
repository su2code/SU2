/*!
 * \file CFileWriter.hpp
 * \brief Headers fo the file writer class.
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
#include "../../../../Common/include/mpi_structure.hpp"
#include "../../../../Common/include/option_structure.hpp"
#include <sys/stat.h>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <time.h>

#include "../../output/filewriter/CParallelDataSorter.hpp"

using namespace std;

class CFileWriter{
protected:

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
  std::string fileExt;

  su2double startTime, /*!< \brief Start time for time measurement  */
            stopTime , /*!< \brief Stop time for time measurement  */
            usedTime , /*!< \brief Used time for time measurement  */
            bandwidth, /*!< \brief Used bandwidth  */
            fileSize; /*!< \brief Size of the last written file  */

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
  
#ifdef HAVE_MPI
  /*!
   * \brief The displacement that every process has in the current file view
   */
  MPI_Offset disp;

  /*!
   * \brief The file handle for writing
   */
  MPI_File fhw;
#else

  /*!
   * \brief The file handle for writing
   */
  FILE* fhw;
#endif

public:
  /*!
   * \brief Construct a file writer using field names, the data sorter and the file extension.
   * \param[in] valFileName - The name of the file
   * \param[in] valDataSorter - The parallel sorted data to write
   * \param[in] valFileExt - The file extension.
   */
  CFileWriter(string valFileName, CParallelDataSorter* valDataSorter, string valFileExt);

  /*!
   * \brief Construct a file writer using field names, file extension.
   * \param[in] valFileName - The name of the file
   * \param[in] valFileExt - The file extension to be attached to the filename
   */
  CFileWriter(string valFileName, string valFileExt);

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
  su2double Get_Bandwidth() const {return bandwidth;}

  /*!
   * \brief Get the filesize of the last written file.
   */
  su2double Get_Filesize() const {return fileSize;}

  /*!
   * \brief Get the used time of the last file writing.
   * \return The time used to write to file.
   */
  su2double Get_UsedTime() const {return usedTime;}
  
protected:
  
  /*!
   * \brief Collectively write a binary data array distributed over all processors to file using MPI I/O.
   * \param[in] data - Pointer to the data to write.
   * \param sizeInBytes - The size of the data in bytes on this processor.
   * \param totalSizeInBytes - The total size of the array accumulated over all processors.
   * \param offset - The offset in bytes of the chunk of data the current processor owns within the global array.
   * \return Boolean indicating whether the writing was successful.
   */
  bool WriteMPIBinaryDataAll(const void *data, unsigned long sizeInBytes, unsigned long totalSizeInBytes, unsigned long offset);

  /*!
   * \brief Write a binary data array to a currently opened file using MPI I/O. Note: routine must be called collectively,
   * although only one processor writes its data.
   * \param[in] data - Pointer to the beginning of the data. Can be NULL for all processors that do not write.
   * \param[in] sizeInBytes - The size of the data in bytes.
   * \param[in] processor - Rank of the processor that should write its data.
   * \return Boolean indicating whether the writing was successful.
   */
  bool WriteMPIBinaryData(const void *data, unsigned long sizeInBytes, unsigned short processor);
  
  /*!
   * \brief Write a string to a currently opened file using MPI I/O. Note: routine must be called collectively,
   *  although only one processor writes the string.
   * \param[in] str - The string to write to file.
   * \param[in] processor - Rank of the processor that should the string.
   * \return 
   */
  bool WriteMPIString(const std::string& str, unsigned short processor);
  
  /*!
   * \brief Open a file to write using MPI I/O. Already existing file is deleted.
   * \return Boolean indicating whether the opening was successful.
   */
  bool OpenMPIFile();
  
  /*!
   * \brief Close a file using MPI I/O. 
   * \return Boolean indicating whether the closing was successful.
   */
  bool CloseMPIFile();
  
};

