#include "../../../include/output/filewriter/CSTLFileWriter.hpp"
#include "../../../include/output/filewriter/CParallelDataSorter.hpp"
#include <iomanip> // for setprecision

CSTLFileWriter::CSTLFileWriter(vector<string> fields, unsigned short nDim) : 
  CFileWriter(fields, nDim){

  file_ext = ".stl";

}


CSTLFileWriter::~CSTLFileWriter(){
  
}


void CSTLFileWriter::Write_Data(string filename, CParallelDataSorter *data_sorter){
  cout << "CSTLFileWriter::Write_Data" << endl;
  filename += file_ext;
  
  /*--- Routine to write the surface CSV files (ASCII). We
   assume here that, as an ASCII file, it is safer to merge the
   surface data onto the master rank for writing for 2 reasons:
   (a) as a surface file, the amount of data should be much less
   than the volume solution, and (b) writing ASCII files in parallel
   requires serializing the IO calls with barriers, which ruins
   the performance at moderate to high rank counts. ---*/
  
  unsigned short iVar;
  
  int iProcessor, nProcessor = size;
  
  unsigned long iPoint, index, iElem;
  unsigned long Buffer_Send_nTriaAll[1], *Buffer_Recv_nTriaAll = NULL;
  unsigned long MaxLocalTriaAll = 0, 
                nLocalTria = 0, 
                nLocalQuad = 0, 
                nLocalTriaAll = 0;
    
  ofstream Surf_file;
  Surf_file.precision(6);
  
  /*--- Find the max number of surface vertices among all
   partitions so we can set up buffers. The master node will handle
   the writing of the CSV file after gathering all of the data. ---*/
  
  nLocalTria = data_sorter->GetnElem(TRIANGLE);
  nLocalQuad = data_sorter->GetnElem(QUADRILATERAL);
  nLocalTriaAll = nLocalTria + nLocalQuad*2; // Quad splitted into 2 tris
  cout << "Rank: " << rank << " , nLocalTria: " << nLocalTria << endl;
  cout << "Rank: " << rank << " , nLocalQuad: " << nLocalQuad << endl;
  cout << "Rank: " << rank << " , nLocalTriaAll: " << nLocalTriaAll << endl;

  Buffer_Send_nTriaAll[0] = nLocalTriaAll;
  if (rank == MASTER_NODE) Buffer_Recv_nTriaAll = new unsigned long[nProcessor];
  
  /*--- Communicate the number of local vertices on each partition
   to the master node with collective calls. ---*/
  
  SU2_MPI::Allreduce(&nLocalTriaAll, &MaxLocalTriaAll, 1,
                     MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  
  SU2_MPI::Gather(&Buffer_Send_nTriaAll, 1, MPI_UNSIGNED_LONG,
                  Buffer_Recv_nTriaAll,  1, MPI_UNSIGNED_LONG,
                  MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- Allocate buffers for send/recv of the data and global IDs. ---*/
  
  su2double *bufD_Send = new su2double[MaxLocalTriaAll*3*3]; // Triangle has 3 Points with 3 coords each, holds all coordinates
  su2double *bufD_Recv = NULL;
  
  /*--- Load send buffers with the local data on this rank. ---*/
  // Tria data
  index = 0;
  for (iElem = 0; iElem < nLocalTria; iElem++) {
    for (iPoint = 0; iPoint < 3; iPoint++) {
      /*--- Solution data. ---*/
      for (iVar = 0; iVar < 3; iVar++){
        bufD_Send[index] = data_sorter->GetData(iVar, data_sorter->GetElem_Connectivity(TRIANGLE, iElem, iPoint) - 1); // (var, GlobalPointindex)
        cout << "bufD_Send[index]: " << bufD_Send[index] << endl;
        index++;
      }
    }  
  }
  // Quad data
  for (iElem = 0; iElem < nLocalQuad; iElem++) {

    vector<unsigned short> Nodelist = {0,1,3, 1,2,3}; //assumes clockwise or counterclockwise rotation
    for (unsigned short iPoint = 0; iPoint < Nodelist.size(); iPoint++) {

      /*--- Solution data. ---*/
      for (iVar = 0; iVar < 3; iVar++){
        bufD_Send[index] = data_sorter->GetData(iVar, data_sorter->GetElem_Connectivity(QUADRILATERAL,iElem,Nodelist[iPoint]) - 1); //TK:: Here the data sort is already messy
        cout << "bufD_Send[index]: " << bufD_Send[index] << endl;
        index++;
      }
    }
  }

  /*--- Only the master rank allocates buffers for the recv. ---*/
  
  if (rank == MASTER_NODE) {
    bufD_Recv = new su2double[nProcessor*MaxLocalTriaAll*3*3];
  }
  
  /*--- Collective comms of the solution data and global IDs. ---*/
  
  SU2_MPI::Gather(bufD_Send, (int)MaxLocalTriaAll*3*3, MPI_DOUBLE,
                  bufD_Recv, (int)MaxLocalTriaAll*3*3, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  
  
  /*--- The master rank alone writes the surface CSV file. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*--- Open the CSV file and write the header with variable names. ---*/
    
    Surf_file.open(filename.c_str(), ios::out);
    Surf_file << "solid SU2_output" << endl;
    /*--- Loop through all of the collected data and write each node's values ---*/
    
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nTriaAll[iProcessor]; iElem++) { // loops over nLocalTriaAll

        /*--- Write the solution data for each field variable. ---*/
        Surf_file << "facet normal " << 1 << " " << 2 << " " << 3 << endl;
        Surf_file << "    outer loop" << endl;
        for(iPoint = 0; iPoint < 3; iPoint++) {
          Surf_file << "        vertex";
          for (iVar = 0; iVar < 3; iVar++){
            Surf_file << " " <<  bufD_Recv[iProcessor*MaxLocalTriaAll*3*3 + iElem*3*3 + iPoint*3 + iVar]; // TK:: check, writes the correct data
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
  
  if (rank == MASTER_NODE) {
    delete [] bufD_Recv;
    delete [] Buffer_Recv_nTriaAll;
  }
}
