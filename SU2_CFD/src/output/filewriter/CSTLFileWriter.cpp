#include "../../../include/output/filewriter/CSTLFileWriter.hpp"
#include "../../../include/output/filewriter/CParallelDataSorter.hpp"
#include "../../../Common/include/geometry_structure.hpp"
#include <iomanip>

CSTLFileWriter::CSTLFileWriter(vector<string> fields, unsigned short nDim, CGeometry *geometry, CConfig *config) : 
  CFileWriter(fields, nDim){

  file_ext = ".stl";

  if (geometry != NULL)
    cout << "Geometry unequal NULL in STLWriter." << endl;

  unsigned int iElem, iPoint;
  unsigned short iDim, iNode, NODES_PER_ELEMENT, iMarker;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {

        switch (geometry->bound[iMarker][iElem]->GetVTK_Type()) {
          case TRIANGLE:
            cout << "TRIANGLE found." << endl;
            NODES_PER_ELEMENT = 3;

            // Loop over all Nodes in one element
            for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
              iPoint = geometry->elem[iElem]->GetNode(iNode);

              // Write one point in the vector as three doubles
              for(iDim = 0; iDim < 3; iDim++) {
                Pointlist.push_back(geometry->node[iPoint]->GetCoord(iDim));
              }
            }
            break;

          case QUADRILATERAL:
            cout << "QUADRILATERAL found." << endl;
            // split quad into 2 tris by taking node numbers (0,1,3) and (1,2,3)
            // this assumes the nodes are counted (counter-)clockwise
            
            // Loop over an array of fixed nodes
            vector<unsigned short> Nodelist = {0,1,3, 1,2,3}; 
            for (iNode = 0; iNode < Nodelist.size(); iNode++) {
              iPoint = geometry->bound[iMarker][iElem]->GetNode(Nodelist[iNode]);

              // Write one point in the vector as three doubles
              for(iDim = 0; iDim < 3; iDim++) {
                Pointlist.push_back(geometry->node[iPoint]->GetCoord(iDim));
              }
            }
            break;
        }
    }
}


CSTLFileWriter::~CSTLFileWriter(){
  
}

void CSTLFileWriter::Write_Data(string filename, CParallelDataSorter *data_sorter){
  
  filename += file_ext;
  
  /*--- Routine to write the surface STL files (ASCII). We
   assume here that, as an ASCII file, it is safer to merge the
   surface data onto the master rank for writing for 2 reasons:
   (a) as a surface file, the amount of data should be much less
   than the volume solution, and (b) writing ASCII files in parallel
   requires serializing the IO calls with barriers, which ruins
   the performance at moderate to high rank counts. ---*/
  
  unsigned short iVar;
  
  int iProcessor, nProcessor = size;
  
  unsigned long iPoint, index;
  unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = NULL;
  unsigned long nLocalVertex_Surface = 0, MaxLocalVertex_Surface = 0;
    
  ofstream Surf_file;
  Surf_file.precision(15);
  
  /*--- Find the max number of surface vertices among all
   partitions so we can set up buffers. The master node will handle
   the writing of the STL file after gathering all of the data. ---*/
  
  nLocalVertex_Surface   = data_sorter->GetnPoints();
  Buffer_Send_nVertex[0] = nLocalVertex_Surface;
  if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long[nProcessor];
  
  /*--- Communicate the number of local vertices on each partition
   to the master node with collective calls. ---*/
  
  SU2_MPI::Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1,
                     MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  
  SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG,
                  Buffer_Recv_nVertex,  1, MPI_UNSIGNED_LONG,
                  MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- Allocate buffers for send/recv of the data and global IDs. ---*/
  
  su2double *bufD_Send = new su2double[MaxLocalVertex_Surface*fieldnames.size()];
  su2double *bufD_Recv = NULL;
  
  unsigned long *bufL_Send = new unsigned long [MaxLocalVertex_Surface];
  unsigned long *bufL_Recv = NULL;
  
  /*--- Load send buffers with the local data on this rank. ---*/
  
  index = 0;
  for (iPoint = 0; iPoint < nLocalVertex_Surface; iPoint++) {
    
    /*--- Global index values. ---*/
    
    bufL_Send[iPoint] = data_sorter->GetGlobalIndex(iPoint);
    
    /*--- Solution data. ---*/
    
    for (iVar = 0; iVar < fieldnames.size(); iVar++){
      bufD_Send[index] = data_sorter->GetData(iVar, iPoint);
      index++;
    }
    
  }
  
  /*--- Only the master rank allocates buffers for the recv. ---*/
  
  if (rank == MASTER_NODE) {
    bufD_Recv = new su2double[nProcessor*MaxLocalVertex_Surface*fieldnames.size()];
    bufL_Recv = new unsigned long[nProcessor*MaxLocalVertex_Surface];
  }
  
  /*--- Collective comms of the solution data and global IDs. ---*/
  
  SU2_MPI::Gather(bufD_Send, (int)MaxLocalVertex_Surface*fieldnames.size(), MPI_DOUBLE,
                  bufD_Recv, (int)MaxLocalVertex_Surface*fieldnames.size(), MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  
  SU2_MPI::Gather(bufL_Send, (int)MaxLocalVertex_Surface, MPI_UNSIGNED_LONG,
                  bufL_Recv, (int)MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- The master rank alone writes the surface STL file. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*--- Open the STL file and write the header with variable names. ---*/
    
    Surf_file.open(filename.c_str(), ios::out);
    Surf_file.precision(6);
    Surf_file << "solid SU2_output" << endl;
    
    /*--- Loop through all of the collected data and write each node's values ---*/
      unsigned int iTri;
      for(iTri = 0; iTri < Pointlist.size()/9; iTri++) {
        /*--- Write the solution data for each field variable. ---*/
        Surf_file << "facet normal " << 1 << " " << 2 << " " << 3 << endl;
        Surf_file << "    outer loop" << endl;
        for(iPoint = 0; iPoint < 3; iPoint++) {
          Surf_file << "        vertex";
          for (iVar = 0; iVar < 3; iVar++){
            Surf_file << " " <<  Pointlist[iTri*9 + iPoint*3 + iVar];
          }
          Surf_file << endl;
        }
        Surf_file << "    endloop" << endl;
        Surf_file << "endfacet" << endl;
      }
    /*--- Close the file. ---*/
    Surf_file << "endsolid SU2_output" << endl;
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
