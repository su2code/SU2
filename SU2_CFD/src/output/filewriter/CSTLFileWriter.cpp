#include "../../../include/output/filewriter/CSTLFileWriter.hpp"
#include "../../../include/output/filewriter/CParallelDataSorter.hpp"
#include <iomanip> /*--- for setprecision ---*/

const string CSTLFileWriter::fileExt = ".stl";

CSTLFileWriter::CSTLFileWriter(vector<string> fields, unsigned short nDim, 
                                         string fileName, CParallelDataSorter *dataSorter) : 
  CFileWriter(std::move(fields), std::move(fileName), dataSorter, fileExt, nDim){}


CSTLFileWriter::~CSTLFileWriter(){}


void CSTLFileWriter::Write_Data(){

  /*--- Routine to write the surface STL files (ASCII). We
   assume here that, as an ASCII file, it is safer to merge the
   surface data onto the master rank for writing for 2 reasons:
   (a) as a surface file, the amount of data should be much less
   than the volume solution, and (b) writing ASCII files in parallel
   requires serializing the IO calls with barriers, which ruins
   the performance at moderate to high rank counts. ---*/

  unsigned short iVar,
                 iPoint;

  unsigned long iProcessor,
                nProcessor = size,
                index,
                iElem,
                MaxLocalTriaAll,
                nLocalTria,
                nLocalQuad,
                nLocalTriaAll = 0,
                *Buffer_Recv_nTriaAll = NULL;

  su2double *bufD_Send = NULL,
            *bufD_Recv = NULL;

  vector<unsigned short> Nodelist = {0,1,3, 1,2,3}; // for Quad2Tri, assumes clockwise or counterclockwise rotation

  ofstream Surf_file;

  /*--- Find the max number of surface vertices among all
   partitions so we can set up buffers. The master node will handle
   the writing of the CSV file after gathering all of the data. ---*/

  nLocalTria = dataSorter->GetnElem(TRIANGLE);
  nLocalQuad = dataSorter->GetnElem(QUADRILATERAL);
  nLocalTriaAll = nLocalTria + nLocalQuad*2; // Quad splitted into 2 tris

  if (rank == MASTER_NODE) Buffer_Recv_nTriaAll = new unsigned long[nProcessor];

  /*--- Communicate the maximum of local triangles on any process to each partition and the number of local vertices on each partition
   to the master node with collective calls. ---*/

  SU2_MPI::Allreduce(&nLocalTriaAll, &MaxLocalTriaAll, 1,
                     MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

  cout << "Rank: " << rank << " ltria " << nLocalTria << " lquad " << nLocalQuad << " max " << MaxLocalTriaAll << endl;

  SU2_MPI::Gather(&nLocalTriaAll,       1, MPI_UNSIGNED_LONG,
                  Buffer_Recv_nTriaAll, 1, MPI_UNSIGNED_LONG,
                  MASTER_NODE, MPI_COMM_WORLD);

  if (rank == MASTER_NODE) {
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      cout << "Buffer_Recv_nTriaAll " << iProcessor << " " <<  Buffer_Recv_nTriaAll[iProcessor] << endl;
    }
  }

  /*--- Allocate buffers for send/recv of the data and global IDs. ---*/

  bufD_Send = new su2double[MaxLocalTriaAll*3*3](); // Triangle has 3 Points with 3 coords each, holds all coordinates

  /*--- Only the master rank allocates buffers for the recv. ---*/
  if (rank == MASTER_NODE)
    bufD_Recv = new su2double[nProcessor*MaxLocalTriaAll*3*3];

  /*--- Load send buffers with the local data on this rank. ---*/
  // Tria data
  index = 0;
  for (iElem = 0; iElem < nLocalTria; iElem++) {
    for (iPoint = 0; iPoint < 3; iPoint++) {
      for (iVar = 0; iVar < 3; iVar++){
        bufD_Send[index] = dataSorter->GetData(iVar, dataSorter->GetElem_Connectivity(TRIANGLE, iElem, iPoint) - 1); // (var, GlobalPointindex)
        if ((abs(bufD_Send[index]) < 1e-4 || abs(bufD_Send[index]) > 1e5) && bufD_Send[index] != 0) {
              cout << "Bad bufD_Send value: " << bufD_Send[index] << " at iproc " << rank << " iElem " << iElem << " iPoint " << iPoint << " iVat " << iVar << endl;
        }
        index++;
      }
    }
  }
  // Quad data
  for (iElem = 0; iElem < nLocalQuad; iElem++) {
    for (iPoint = 0; iPoint < Nodelist.size(); iPoint++) {
      for (iVar = 0; iVar < 3; iVar++){
        bufD_Send[index] = dataSorter->GetData(iVar, dataSorter->GetElem_Connectivity(QUADRILATERAL,iElem,Nodelist[iPoint]) - 1);
        index++;
      }
    }
  }

  /*--- Collective comms of the solution data and global IDs. ---*/
  SU2_MPI::Gather(bufD_Send, static_cast<int>(MaxLocalTriaAll*3*3), MPI_DOUBLE,
                  bufD_Recv, static_cast<int>(MaxLocalTriaAll*3*3), MPI_DOUBLE,
                  MASTER_NODE, MPI_COMM_WORLD);

  /*--- The master rank alone writes the surface CSV file. ---*/

  if (rank == MASTER_NODE) {

    /*--- Open the CSV file and write the header with variable names. ---*/
    Surf_file.precision(6);
    Surf_file.open(fileName.c_str(), ios::out);
    Surf_file << "solid SU2_output" << endl;
    /*--- Loop through all of the collected data and write each node's values ---*/

    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      cout << "iProcessor: " << iProcessor << " out of " << nProcessor << endl;
      for (iElem = 0; iElem < Buffer_Recv_nTriaAll[iProcessor]; iElem++) { // loops over nLocalTriaAll

        /*--- Write the solution data for each field variable. ---*/
        Surf_file << "facet normal " << 1 << " " << 2 << " " << 3 << endl;
        Surf_file << "    outer loop" << endl;

        for(iPoint = 0; iPoint < 3; iPoint++) {
          Surf_file << "        vertex";
          index = iProcessor*MaxLocalTriaAll*3*3 + iElem*3*3 + iPoint*3;

          for (iVar = 0; iVar < 3; iVar++) {
            Surf_file << " " <<  bufD_Recv[index + iVar];
            if ((abs(bufD_Recv[index + iVar]) < 1e-4 || abs(bufD_Recv[index + iVar]) > 1e5) && bufD_Recv[index + iVar] != 0) {
              cout << "Bad bufD_Recv value: " << bufD_Recv[index + iVar] << " at iproc " << iProcessor << " iElem " << iElem << " iPoint " << iPoint << " iVat " << iVar << endl;
            } 
          }
          Surf_file << endl;
        }
        Surf_file << "    endloop" << endl;
        Surf_file << "endfacet" << endl;
      }//iElem
    }//iProcessor

    /*--- Close the file. ---*/
    Surf_file << "endsolid SU2_output" << endl;
    Surf_file.close();

  }

  /*--- Free temporary memory. ---*/
  if(bufD_Send != NULL) delete [] bufD_Send;
  if(bufD_Recv != NULL) delete [] bufD_Recv;
  if(Buffer_Recv_nTriaAll != NULL) delete [] Buffer_Recv_nTriaAll;
}
