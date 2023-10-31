/*!
 * \file CSU2FileWriter.cpp
 * \brief Filewriter class SU2 native ASCII (CSV) format.
 * \author T. Albring
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

  ofstream restart_file;
  const vector<string> fieldNames = dataSorter->GetFieldNames();
  /*--- We append the pre-defined suffix (extension) to the filename (prefix) ---*/
  val_filename.append(fileExt);

  /*--- Set a timer for the file writing. ---*/

  startTime = SU2_MPI::Wtime();

  /*--- Only the FIRST node writes the header (it does not matter if that is the master). ---*/

  if (rank == 0) {
    restart_file.open(val_filename);
    restart_file << "\"PointID\"";
    for (auto& field : fieldNames) restart_file << ",\"" << field << "\"";
    restart_file << "\n";
    restart_file.close();
  }

  /*--- Serialize the writes to the restart file. ---*/

  for (int iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      restart_file.open(val_filename, ios::app);
      restart_file.precision(15);

      for (auto iPoint = 0ul; iPoint < dataSorter->GetnPoints(); iPoint++) {

        /*--- Write global index of the current point. ---*/

        restart_file << dataSorter->GetGlobalIndex(iPoint);

        /*--- Loop over the variables and write the values to file. ---*/

        for (size_t iVar = 0; iVar < fieldNames.size(); iVar++)
          restart_file << ", " << scientific << dataSorter->GetData(iVar, iPoint);
        restart_file << "\n";
      }

      restart_file.close();
    }

    /*--- Wait for iProcessor to finish and close the file. ---*/

    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }

  /*--- Compute and store the write time. ---*/

  stopTime = SU2_MPI::Wtime();

  usedTime = stopTime-startTime;

  /*--- Determine the file size ---*/

  fileSize = DetermineFilesize(val_filename);

  /*--- Compute and store the bandwidth ---*/

  bandwidth = fileSize/(1.0e6)/usedTime;
}
