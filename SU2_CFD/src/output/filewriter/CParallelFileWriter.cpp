/*!
 * \file CFileWriter.cpp
 * \brief Filewriter base class.
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

#include "../../../include/output/filewriter/CFileWriter.hpp"


CFileWriter::CFileWriter(string valFileName, CParallelDataSorter *valDataSorter, string valFileExt):
  fileExt(valFileExt),
  fileName(std::move(valFileName)),
  dataSorter(valDataSorter){

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  this->fileName += valFileExt;

  fileSize = 0.0;
  bandwidth = 0.0;
  
}

CFileWriter::CFileWriter(string valFileName, string valFileExt):
  fileExt(valFileExt),
  fileName(std::move(valFileName)){
  
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  this->fileName += valFileExt;
  
  fileSize = 0.0;
  bandwidth = 0.0;
  
}

CFileWriter::~CFileWriter(){

}


