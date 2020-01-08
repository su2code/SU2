/*!
 * \file CSTLFileWriter.cpp
 * \brief STL Writer output class
 * \author T. Kattmann
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

#include "../../../include/output/filewriter/CSTLFileWriter.hpp"
#include "../../../include/output/filewriter/CParallelDataSorter.hpp"
#include <iomanip> /*--- for setprecision ---*/

const string CSTLFileWriter::fileExt = ".stl";

CSTLFileWriter::CSTLFileWriter(vector<string> fields,
                               unsigned short nDim,
                               string fileName,
                               CParallelDataSorter *dataSorter) :
                CFileWriter(std::move(fields),
                            std::move(fileName),
                            dataSorter,
                            fileExt,
                            nDim){}


CSTLFileWriter::~CSTLFileWriter(){}


void CSTLFileWriter::Write_Data(){

  unsigned long nParallel_Tria = dataSorter->GetnElem(TRIANGLE),
                nParallel_Quad = dataSorter->GetnElem(QUADRILATERAL);

  unsigned short iVar;
  //  NodePartitioner node_partitioner(num_nodes, size);
  num_nodes_to_receive.resize(size, 0);
  values_to_receive_displacements.resize(size);
  halo_nodes.clear();
  sorted_halo_nodes.clear();
  
  /* We output a single, partitioned zone where each rank outputs one partition. */
  vector<int32_t> partition_owners;
  partition_owners.reserve(size);
  for (int32_t iRank = 0; iRank < size; ++iRank)
    partition_owners.push_back(iRank);

  /* Gather a list of nodes we refer to but are not outputting. */

  for (unsigned long i = 0; i < nParallel_Tria * N_POINTS_TRIANGLE; ++i)
    if (dataSorter->FindProcessor(dataSorter->GetElem_Connectivity(TRIANGLE, 0, i) -1) != rank)
      halo_nodes.insert(dataSorter->GetElem_Connectivity(TRIANGLE, 0, i)-1);

  for (unsigned long i = 0; i < nParallel_Quad * N_POINTS_QUADRILATERAL; ++i)
    if (dataSorter->FindProcessor(dataSorter->GetElem_Connectivity(QUADRILATERAL, 0, i) -1) != rank)
      halo_nodes.insert(dataSorter->GetElem_Connectivity(QUADRILATERAL, 0, i)-1);
  

  /* Sorted list of halo nodes for this MPI rank. */
  sorted_halo_nodes.assign(halo_nodes.begin(), halo_nodes.end());
      
  /*--- We effectively tack the halo nodes onto the end of the node list for this partition.
    TecIO will later replace them with references to nodes in neighboring partitions. */
  num_halo_nodes = sorted_halo_nodes.size();
  vector<int64_t> halo_node_local_numbers(max((size_t)1, num_halo_nodes)); /* Min size 1 to avoid crashes when we access these vectors below. */
  vector<int32_t> neighbor_partitions(max((size_t)1, num_halo_nodes));
  vector<int64_t> neighbor_nodes(max((size_t)1, num_halo_nodes));
  for(int64_t i = 0; i < static_cast<int64_t>(num_halo_nodes); ++i) {
    halo_node_local_numbers[i] = dataSorter->GetNodeEnd(rank) - dataSorter->GetNodeBegin(rank) + i;
    int owning_rank = dataSorter->FindProcessor(sorted_halo_nodes[i]);
    unsigned long node_number = sorted_halo_nodes[i] - dataSorter->GetNodeBegin(owning_rank);
    neighbor_partitions[i] = owning_rank; /* Partition numbers are 1-based. */
    if (rank == 0) cout << owning_rank << endl;
    neighbor_nodes[i] = static_cast<int64_t>(node_number);
  }

  /* Gather halo node data. First, tell each rank how many nodes' worth of data we need from them. */
  for (size_t i = 0; i < num_halo_nodes; ++i)
    ++num_nodes_to_receive[neighbor_partitions[i]];
  num_nodes_to_send.resize(size);
  SU2_MPI::Alltoall(&num_nodes_to_receive[0], 1, MPI_INT, &num_nodes_to_send[0], 1, MPI_INT, MPI_COMM_WORLD);

  /* Now send the global node numbers whose data we need,
     and receive the same from all other ranks.
     Each rank has globally consecutive node numbers,
     so we can just parcel out sorted_halo_nodes for send. */
  nodes_to_send_displacements.resize(size);
  nodes_to_receive_displacements.resize(size);
  nodes_to_send_displacements[0] = 0;
  nodes_to_receive_displacements[0] = 0;
  for(int iRank = 1; iRank < size; ++iRank) {
    nodes_to_send_displacements[iRank] = nodes_to_send_displacements[iRank - 1] + num_nodes_to_send[iRank - 1];
    nodes_to_receive_displacements[iRank] = nodes_to_receive_displacements[iRank - 1] + num_nodes_to_receive[iRank - 1];
  }
  int total_num_nodes_to_send = nodes_to_send_displacements[size - 1] + num_nodes_to_send[size - 1];

  nodes_to_send.resize(max(1, total_num_nodes_to_send));
  /* The terminology gets a bit confusing here. We're sending the node numbers
     (sorted_halo_nodes) whose data we need to receive, and receiving
     lists of nodes whose data we need to send. */
  if (sorted_halo_nodes.empty()) sorted_halo_nodes.resize(1); /* Avoid crash. */
  SU2_MPI::Alltoallv(&sorted_halo_nodes[0], &num_nodes_to_receive[0], &nodes_to_receive_displacements[0], MPI_UNSIGNED_LONG,
                     &nodes_to_send[0],     &num_nodes_to_send[0],    &nodes_to_send_displacements[0],    MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
  
  /* Now actually send and receive the data */
  data_to_send.resize(max(1, total_num_nodes_to_send * (int)fieldnames.size()));
  halo_var_data.resize(max((size_t)1, fieldnames.size() * num_halo_nodes));
  num_values_to_send.resize(size);
  values_to_send_displacements.resize(size);
  num_values_to_receive.resize(size);
  size_t index = 0;
  for(int iRank = 0; iRank < size; ++iRank) {
    /* We send and receive GlobalField_Counter values per node. */
    num_values_to_send[iRank]              = num_nodes_to_send[iRank] * fieldnames.size();
    values_to_send_displacements[iRank]    = nodes_to_send_displacements[iRank] * fieldnames.size();
    num_values_to_receive[iRank]           = num_nodes_to_receive[iRank] * fieldnames.size();
    values_to_receive_displacements[iRank] = nodes_to_receive_displacements[iRank] * fieldnames.size();
    for(iVar = 0; iVar < fieldnames.size(); ++iVar)
      for(int iNode = 0; iNode < num_nodes_to_send[iRank]; ++iNode) {
        unsigned long node_offset = nodes_to_send[nodes_to_send_displacements[iRank] + iNode] - dataSorter->GetNodeBegin(rank);
        data_to_send[index++] = SU2_TYPE::GetValue(dataSorter->GetData(iVar,node_offset));
      }
  }
  SU2_MPI::Alltoallv(&data_to_send[0],  &num_values_to_send[0],    &values_to_send_displacements[0],    MPI_DOUBLE,
                     &halo_var_data[0], &num_values_to_receive[0], &values_to_receive_displacements[0], MPI_DOUBLE,
                     MPI_COMM_WORLD);

  /*--- Routine to write the surface STL files (ASCII). We
   assume here that, as an ASCII file, it is safer to merge the
   surface data onto the master rank for writing for 2 reasons:
   (a) as a surface file, the amount of data should be much less
   than the volume solution, and (b) writing ASCII files in parallel
   requires serializing the IO calls with barriers, which ruins
   the performance at moderate to high rank counts. ---*/

  unsigned short iPoint;

  unsigned long iProcessor,
                nProcessor = size,
                iElem,
                MaxLocalTriaAll,
                nLocalTria,
                nLocalQuad,
                nLocalTriaAll = 0,
                *Buffer_Recv_nTriaAll = NULL,
                global_node_number;

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

      global_node_number = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, iPoint) - 1; // (var, GlobalPointindex)
      unsigned long local_node_number = global_node_number - dataSorter->GetNodeBegin(rank);

      for (iVar = 0; iVar < 3; iVar++){

        if (dataSorter->FindProcessor(global_node_number) == rank) {
          bufD_Send[index] = dataSorter->GetData(iVar, local_node_number); 
        } else {
          bufD_Send[index] = GetHaloNodeValue(global_node_number, iVar);
        }
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
      global_node_number = dataSorter->GetElem_Connectivity(QUADRILATERAL,iElem,Nodelist[iPoint]) - 1; // (var, GlobalPointindex)
      unsigned long local_node_number = global_node_number - dataSorter->GetNodeBegin(rank); 
      //cout << rank << " " <<  global_node_number << " " << dataSorter->GetNodeBegin(rank) << endl;
      //unsigned long global_node_number_volume = dataSorter->GetGlobalIndex(local_node_number);

      for (iVar = 0; iVar < 3; iVar++){

         if (dataSorter->FindProcessor(global_node_number) == rank) {
        //if (local_node_number < dataSorter->GetnPoints()) {

          bufD_Send[index] = dataSorter->GetData(iVar, local_node_number);
          if ((abs(bufD_Send[index]) < 1e-4 || abs(bufD_Send[index]) > 1e5) && bufD_Send[index] != 0) {
            cout << "Bad bufD_Send value Local: " << bufD_Send[index] << " at iproc " << rank << " iElem " << iElem << " iPoint " << iPoint << " iVat " << iVar <<  endl;
           // cout << "global_node_number global renumered" << global_node_number << " local renumbered " << global_node_number - dataSorter->GetNodeBegin(rank) << " global not renumbered " << dataSorter->GetGlobalIndex(global_node_number - dataSorter->GetNodeBegin(rank)) << endl;
          }
        } else {
           bufD_Send[index] = GetHaloNodeValue(global_node_number, iVar);
          // if ((abs(bufD_Send[index]) < 1e-4 || abs(bufD_Send[index]) > 1e5) ) {
            // cout << "Bad bufD_Send value Halo: " << bufD_Send[index] << " at iproc " << rank << " iElem " << iElem << " iPoint " << iPoint << " iVat " << iVar << endl;
          // }
        }
        
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
            //if ((abs(bufD_Recv[index + iVar]) < 1e-4 || abs(bufD_Recv[index + iVar]) > 1e5) && bufD_Recv[index + iVar] != 0) {
            //  cout << "Bad bufD_Recv value: " << bufD_Recv[index + iVar] << " at iproc " << iProcessor << " iElem " << iElem << " iPoint " << iPoint << " iVat " << iVar << endl;
            //} 
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

double CSTLFileWriter::GetHaloNodeValue(unsigned long global_node_number, unsigned short iVar) {

  vector<unsigned long>::iterator it = lower_bound(sorted_halo_nodes.begin(), sorted_halo_nodes.end(), global_node_number);

  int offset = distance(sorted_halo_nodes.begin(), it);
    cout << offset << " " << sorted_halo_nodes.size() << endl;


  int id = 0;
  for (int iRank = 0; iRank < size; ++iRank) {
    cout << rank << " " << num_nodes_to_receive[iRank] << endl;
    for (int i = 0; i < num_nodes_to_receive[iRank]; i++){  
      int displ = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank]*iVar;
      if (id == offset)
        return halo_var_data[displ+i];
      id++;    
    }
  }

  SU2_MPI::Error("Node not found", CURRENT_FUNCTION);
  

  return 0.0; 
}
