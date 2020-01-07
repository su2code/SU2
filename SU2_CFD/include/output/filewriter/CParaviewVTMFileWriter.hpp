/*!
 * \file CParaviewVTMFileWriter.hpp
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

class CParaviewVTMFileWriter final: public CFileWriter{

  /*!
   * \brief String stream that stores the content of the vtm file
   */
  stringstream output;
  
  /*!
   * \brief The folder name where all the files associated with the datasets will be stored
   */
  string folderName;
  
  /*!
   * \brief The current zone index
   */
  unsigned short iZone;
  
  /*!
   * \brief The total number of zones
   */
  unsigned short nZone;
  
  /*!
   * \brief Current physical time
   */
  su2double curTime;
  
public:

  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names, dimension.
   * \param[in] fileName - The name of the file
   */
  CParaviewVTMFileWriter(string fileName, string folderName, su2double time, unsigned short iZone, unsigned short nZone);

  /*!
   * \brief Destructor
   */
  ~CParaviewVTMFileWriter() override;

  /*!
   * \brief Write sorted data to file in paraview binary file format
   */
  void Write_Data() override;
  
  /*!
   * \brief Add a new dataset by writing data from a datasorter to file and adding it to the vtm file
   * \param[in] name - The name of the dataset
   * \param[in] file - The name of the vtu dataset file to write
   * \param[in] dataSorter - Datasorter object containing the actual data. Note, data must be sorted.
   */
  void AddDataset(string name, string file, CParallelDataSorter* dataSorter);

  /*!
   * \brief Start a new block
   * \param[in] name - The name of the block
   */
  inline void StartBlock(string name){
    if (rank == MASTER_NODE){
      output << "<Block name=\"" << name << "\">" << endl;
    }   
  }
  
  /*!
   * \brief Close currently opened block
   */
  inline void EndBlock(){
    if (rank == MASTER_NODE){
      output << "</Block>" << endl;
    }
  }
  
  /*!
   * \brief Add a new dataset by writing it to the vtm file
   * \param[in] name - Name of the dataset
   * \param[in] file - vtu file where the data is stored
   */
  inline void AddDataset(string name, string file){
    if (rank == MASTER_NODE){
      output << "<DataSet name=\"" << name <<"\" file=\"" << file << "\"/>" << endl;
    }
  }
  
  /*!
   * \brief Get the name of the folder where the data will be stored
   * \return The folder name
   */
  inline string GetFolderName() const{
    return folderName;
  }
};
