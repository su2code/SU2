/*!
 * \file CCSVFileWriter.cpp
 * \brief CSV Writer output class
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

#include "../../../include/output/filewriter/CCSVFileWriter.hpp"
#include "../../../include/output/filewriter/CParallelDataSorter.hpp"

CCSVFileWriter::CCSVFileWriter(CParallelDataSorter *valDataSorter) :
  CFileWriter(valDataSorter, ".csv"){}


CCSVFileWriter::~CCSVFileWriter()= default;

void CCSVFileWriter::WriteData(string val_filename){

  /*--- Routine to write the surface CSV files (ASCII). We
   assume here that, as an ASCII file, it is safer to merge the
   surface data onto the master rank for writing for 2 reasons:
   (a) as a surface file, the amount of data should be much less
   than the volume solution, and (b) writing ASCII files in parallel
   requires serializing the IO calls with barriers, which ruins
   the performance at moderate to high rank counts. ---*/

  unsigned short iVar;

  int iProcessor, nProcessor = size;

  unsigned long iPoint, index;
  unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = nullptr;
  unsigned long nLocalVertex_Surface = 0, MaxLocalVertex_Surface = 0;

  const vector<string> fieldNames = dataSorter->GetFieldNames();

  ofstream Surf_file;
  Surf_file.precision(15);

  /*--- Find the max number of surface vertices among all
   partitions so we can set up buffers. The master node will handle
   the writing of the CSV file after gathering all of the data. ---*/

  nLocalVertex_Surface   = dataSorter->GetnPoints();
  Buffer_Send_nVertex[0] = nLocalVertex_Surface;
  if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long[nProcessor];

  /*--- Communicate the number of local vertices on each partition
   to the master node with collective calls. ---*/

  SU2_MPI::Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1,
                     MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());

  SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG,
                  Buffer_Recv_nVertex,  1, MPI_UNSIGNED_LONG,
                  MASTER_NODE, SU2_MPI::GetComm());

  /*--- Allocate buffers for send/recv of the data and global IDs. ---*/

  auto *bufD_Send = new su2double[MaxLocalVertex_Surface*fieldNames.size()]();
  su2double *bufD_Recv = nullptr;

  auto *bufL_Send = new unsigned long [MaxLocalVertex_Surface]();
  unsigned long *bufL_Recv = nullptr;

  /*--- Load send buffers with the local data on this rank. ---*/

  index = 0;
  for (iPoint = 0; iPoint < nLocalVertex_Surface; iPoint++) {

    /*--- Global index values. ---*/

    bufL_Send[iPoint] = dataSorter->GetGlobalIndex(iPoint);

    /*--- Solution data. ---*/

    for (iVar = 0; iVar < fieldNames.size(); iVar++){
      bufD_Send[index] = dataSorter->GetData(iVar, iPoint);
      index++;
    }

  }

  /*--- Only the master rank allocates buffers for the recv. ---*/

  if (rank == MASTER_NODE) {
    bufD_Recv = new su2double[nProcessor*MaxLocalVertex_Surface*fieldNames.size()]();
    bufL_Recv = new unsigned long[nProcessor*MaxLocalVertex_Surface];
  }

  /*--- Collective comms of the solution data and global IDs. ---*/

  SU2_MPI::Gather(bufD_Send, (int)MaxLocalVertex_Surface*fieldNames.size(), MPI_DOUBLE,
                  bufD_Recv, (int)MaxLocalVertex_Surface*fieldNames.size(), MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());

  SU2_MPI::Gather(bufL_Send, (int)MaxLocalVertex_Surface, MPI_UNSIGNED_LONG,
                  bufL_Recv, (int)MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());

  /*--- The master rank alone writes the surface CSV file. ---*/

  if (rank == MASTER_NODE) {

    /*--- We append the pre-defined suffix (extension) to the filename (prefix) ---*/
    val_filename.append(fileExt);

    /*--- Open the CSV file and write the header with variable names. ---*/
    Surf_file.open(val_filename.c_str(), ios::out);
    Surf_file << "\"Point\",";
    for (iVar = 0; iVar < fieldNames.size()-1; iVar++) {
      Surf_file << "\"" << fieldNames[iVar] << "\",";
    }
    Surf_file << "\"" << fieldNames[fieldNames.size()-1] << "\"" << endl;

    /*--- Loop through all of the collected data and write each node's values ---*/

    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nVertex[iProcessor]; iPoint++) {

        /*--- Current index position for global index access. ---*/

        index  = iProcessor*MaxLocalVertex_Surface + iPoint;

        /*--- Write global index values. ---*/

        Surf_file << bufL_Recv[index] << ", ";

        /*--- Reset index for solution data access. ---*/

        index  = (iProcessor*MaxLocalVertex_Surface*fieldNames.size() +
                  iPoint*fieldNames.size());

        /*--- Write the solution data for each field variable. ---*/

        for (iVar = 0; iVar < fieldNames.size(); iVar++){
          Surf_file << scientific << bufD_Recv[index + iVar];
          if (iVar != fieldNames.size() -1) Surf_file << ", ";
        }
        Surf_file << endl;

      }
    }

    /*--- Close the file. ---*/

    Surf_file.close();

  }

  /*--- Free temporary memory. ---*/

  if (rank == MASTER_NODE) {
    delete [] bufL_Recv;
    delete [] bufD_Recv;
    delete [] Buffer_Recv_nVertex;
  }
  delete [] bufL_Send;
  delete [] bufD_Send;
}
