/*!
 * \file CSU2FileWriter.cpp
 * \brief Filewriter class SU2 native ASCII (CSV) format.
 * \author T. Albring
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/output/filewriter/CSU2FileWriter.hpp"

const string CSU2FileWriter::fileExt = ".csv";

CSU2FileWriter::CSU2FileWriter(CParallelDataSorter *valDataSorter) :
  CFileWriter(valDataSorter, fileExt){}

void CSU2FileWriter::WriteData(string val_filename){

  const vector<string> fieldNames = dataSorter->GetRequiredFieldNames();

  OpenMPIFile(val_filename);

  /*--- Write the header. ---*/

  string header = "\"PointID\"";
  for (auto& field : fieldNames) header += ",\"" + field + "\"";
  header += "\n";

  WriteMPIString(header, MASTER_NODE);

  /*--- Each rank formats the data of its own points into a string, and all ranks then write their strings
   to the file at the same time, one after the other in rank order. ---*/

  ostringstream data;
  data.precision(15);

  for (auto iPoint = 0ul; iPoint < dataSorter->GetnPoints(); iPoint++) {

    /*--- Write global index of the current point. ---*/

    data << dataSorter->GetGlobalIndex(iPoint);

    /*--- Loop over the variables and write the values to file. ---*/

    for (size_t iVar = 0; iVar < fieldNames.size(); iVar++)
      data << ", " << scientific << dataSorter->GetData(iVar, iPoint);
    data << "\n";
  }

  WriteMPIStringAll(data.str());

  CloseMPIFile();
}
