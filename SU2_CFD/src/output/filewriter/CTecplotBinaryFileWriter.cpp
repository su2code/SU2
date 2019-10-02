#include "../../../include/output/filewriter/CTecplotBinaryFileWriter.hpp"
#ifdef HAVE_TECIO
  #include "TECIO.h"
#endif
#include <set>

const string CTecplotBinaryFileWriter::fileExt = ".szplt";

CTecplotBinaryFileWriter::CTecplotBinaryFileWriter(vector<string> fields, unsigned short nDim,
                                                   string fileName, CParallelDataSorter *dataSorter,
                                                   unsigned long time_iter, su2double timestep) : 
  CFileWriter(std::move(fields), std::move(fileName), dataSorter, fileExt, nDim), time_iter(time_iter), timestep(timestep){}

CTecplotBinaryFileWriter::~CTecplotBinaryFileWriter(){}

void CTecplotBinaryFileWriter::Write_Data(){
    
  if (!dataSorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }
  
  /*--- Set a timer for the binary file writing. ---*/
  
#ifndef HAVE_MPI
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StartTime = MPI_Wtime();
#endif  
  
#ifdef HAVE_TECIO
  
  /*--- Reduce the total number of each element. ---*/

  unsigned long nTot_Line, nTot_Tria, nTot_Quad, nTot_Tetr, nTot_Hexa, nTot_Pris, nTot_Pyra;
  unsigned long nParallel_Line = dataSorter->GetnElem(LINE),
                nParallel_Tria = dataSorter->GetnElem(TRIANGLE),
                nParallel_Quad = dataSorter->GetnElem(QUADRILATERAL),
                nParallel_Tetr = dataSorter->GetnElem(TETRAHEDRON),
                nParallel_Hexa = dataSorter->GetnElem(HEXAHEDRON),
                nParallel_Pris = dataSorter->GetnElem(PRISM),
                nParallel_Pyra = dataSorter->GetnElem(PYRAMID);
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nParallel_Line, &nTot_Line, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Tria, &nTot_Tria, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Quad, &nTot_Quad, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Tetr, &nTot_Tetr, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Hexa, &nTot_Hexa, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Pris, &nTot_Pris, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Pyra, &nTot_Pyra, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTot_Line      = nParallel_Line;

  nTot_Tria = nParallel_Tria;
  nTot_Quad = nParallel_Quad;
  nTot_Tetr = nParallel_Tetr;
  nTot_Hexa = nParallel_Hexa;
  nTot_Pris = nParallel_Pris;
  nTot_Pyra = nParallel_Pyra;
#endif
  
  string data_set_title = "Visualization of the solution";

  ostringstream tecplot_variable_names;
  for (size_t iVar = 0; iVar < fieldnames.size()-1; ++iVar) {
    tecplot_variable_names << fieldnames[iVar] << ",";
  }
  tecplot_variable_names << fieldnames[fieldnames.size()-1];

  void* file_handle = NULL;
  int32_t err = tecFileWriterOpen(fileName.c_str(), data_set_title.c_str(), tecplot_variable_names.str().c_str(),
    FILEFORMAT_SZL, FILETYPE_FULL, (int32_t)FieldDataType_Double, NULL, &file_handle);
  if (err) cout << "Error opening Tecplot file '" << fileName << "'" << endl;

#ifdef HAVE_MPI
  err = tecMPIInitialize(file_handle, MPI_COMM_WORLD, MASTER_NODE);
  if (err) cout << "Error initializing Tecplot parallel output." << endl;
#endif
  
  /*--- Define the zone(s). For 2D, and for 3D surfaces, each rank outputs a separate zone. ---*/

  int64_t num_nodes;
  int64_t num_cells;
  int32_t zone_type;
  
  
  num_nodes = static_cast<int64_t>(dataSorter->GetnPointsGlobal());
  num_cells = static_cast<int64_t>(dataSorter->GetnElem());
  if (nDim == 3){
    if ((nTot_Quad > 0 || nTot_Tria > 0) && (nTot_Hexa + nTot_Pris + nTot_Pyra + nTot_Tetr == 0)){
      zone_type = ZONETYPE_FEQUADRILATERAL;
    }
    else {
      zone_type = ZONETYPE_FEBRICK;      
    }
  }
  else {
    if (nTot_Line > 0 && (nTot_Tria + nTot_Quad == 0)){
      zone_type = ZONETYPE_FELINESEG;
      
    }
    else{
      zone_type = ZONETYPE_FEQUADRILATERAL;      
    }
  }

  bool is_unsteady = false;
  passivedouble solution_time = 0.0;

  if (timestep > 0.0){
    is_unsteady = true;
    solution_time = SU2_TYPE::GetValue(timestep)*time_iter;
  }
  
  int32_t zone;
  vector<int32_t> value_locations(fieldnames.size(), 1); /* Nodal variables. */
  err = tecZoneCreateFE(file_handle, "Zone", zone_type, num_nodes, num_cells, NULL, NULL, &value_locations[0], NULL, 0, 0, 0, &zone);
  if (err) cout << rank << ": Error creating Tecplot zone." << endl;
  if (is_unsteady) {
    err = tecZoneSetUnsteadyOptions(file_handle, zone, solution_time, time_iter + 1);
    if (err) cout << rank << ": Error setting Tecplot zone unsteady options." << std::endl;
  }

#ifdef HAVE_MPI

  unsigned short iVar;
  NodePartitioner node_partitioner(num_nodes, size);
  std::set<unsigned long> halo_nodes;
  vector<unsigned long> sorted_halo_nodes;
  vector<passivedouble> halo_var_data;
  vector<int> num_nodes_to_receive(size, 0);
  vector<int> values_to_receive_displacements(size);

  if (zone_type == ZONETYPE_FEBRICK) {

    /* We output a single, partitioned zone where each rank outputs one partition. */
    vector<int32_t> partition_owners;
    partition_owners.reserve(size);
    for (int32_t iRank = 0; iRank < size; ++iRank)
      partition_owners.push_back(iRank);
    err = tecZoneMapPartitionsToMPIRanks(file_handle, zone, size, &partition_owners[0]);
    if (err) cout << rank << ": Error assigning MPI ranks for Tecplot zone partitions." << endl;
  
    /* Gather a list of nodes we refer to but are not outputting. */

    for (unsigned long i = 0; i < nParallel_Tria * N_POINTS_TRIANGLE; ++i)
      if ((unsigned long)dataSorter->GetElem_Connectivity(TRIANGLE, 0, i) <= dataSorter->GetNodeBegin(rank) || 
          dataSorter->GetNodeEnd(rank) < (unsigned long)dataSorter->GetElem_Connectivity(TRIANGLE, 0, i))
        halo_nodes.insert(dataSorter->GetElem_Connectivity(TRIANGLE, 0, i));
  
    for (unsigned long i = 0; i < nParallel_Quad * N_POINTS_QUADRILATERAL; ++i)
      if ((unsigned long)dataSorter->GetElem_Connectivity(QUADRILATERAL, 0, i) <= dataSorter->GetNodeBegin(rank) || 
          dataSorter->GetNodeEnd(rank) < (unsigned long)dataSorter->GetElem_Connectivity(QUADRILATERAL, 0, i))
        halo_nodes.insert(dataSorter->GetElem_Connectivity(QUADRILATERAL, 0, i));

    for (unsigned long i = 0; i < nParallel_Tetr * N_POINTS_TETRAHEDRON; ++i)
      if ((unsigned long)dataSorter->GetElem_Connectivity(TETRAHEDRON, 0, i) <= dataSorter->GetNodeBegin(rank) || 
          dataSorter->GetNodeEnd(rank) < (unsigned long)dataSorter->GetElem_Connectivity(TETRAHEDRON, 0, i))
        halo_nodes.insert(dataSorter->GetElem_Connectivity(TETRAHEDRON, 0, i));

    for (unsigned long i = 0; i < nParallel_Hexa * N_POINTS_HEXAHEDRON; ++i)
      if ((unsigned long)dataSorter->GetElem_Connectivity(HEXAHEDRON, 0, i) <= dataSorter->GetNodeBegin(rank) || 
          dataSorter->GetNodeEnd(rank) < (unsigned long)dataSorter->GetElem_Connectivity(HEXAHEDRON, 0, i))
        halo_nodes.insert(dataSorter->GetElem_Connectivity(HEXAHEDRON, 0, i));
      
    for (unsigned long i = 0; i < nParallel_Pris * N_POINTS_PRISM; ++i)
      if ((unsigned long)dataSorter->GetElem_Connectivity(PRISM, 0, i) <= dataSorter->GetNodeBegin(rank) || 
          dataSorter->GetNodeEnd(rank) < (unsigned long)dataSorter->GetElem_Connectivity(PRISM, 0, i))
        halo_nodes.insert(dataSorter->GetElem_Connectivity(PRISM, 0, i));
    
    for (unsigned long i = 0; i < nParallel_Pyra * N_POINTS_PYRAMID; ++i)
      if ((unsigned long)dataSorter->GetElem_Connectivity(PYRAMID, 0, i) <= dataSorter->GetNodeBegin(rank) || 
          dataSorter->GetNodeEnd(rank) < (unsigned long)dataSorter->GetElem_Connectivity(PYRAMID, 0, i))
        halo_nodes.insert(dataSorter->GetElem_Connectivity(PYRAMID, 0, i));

    /* Sorted list of halo nodes for this MPI rank. */
    sorted_halo_nodes.assign(halo_nodes.begin(), halo_nodes.end());
        
    /* Have to include all nodes our cells refer to or TecIO will barf, so add the halo node count to the number of local nodes. */
    int64_t partition_num_nodes = dataSorter->GetNodeEnd(rank) - dataSorter->GetNodeBegin(rank) + static_cast<int64_t>(halo_nodes.size());
    int64_t partition_num_cells = nParallel_Tetr + nParallel_Hexa + nParallel_Pris + nParallel_Pyra;

    /*--- We effectively tack the halo nodes onto the end of the node list for this partition.
      TecIO will later replace them with references to nodes in neighboring partitions. */
    size_t num_halo_nodes = sorted_halo_nodes.size();
    vector<int64_t> halo_node_local_numbers(max((size_t)1, num_halo_nodes)); /* Min size 1 to avoid crashes when we access these vectors below. */
    vector<int32_t> neighbor_partitions(max((size_t)1, num_halo_nodes));
    vector<int64_t> neighbor_nodes(max((size_t)1, num_halo_nodes));
    for(int64_t i = 0; i < static_cast<int64_t>(num_halo_nodes); ++i) {
      halo_node_local_numbers[i] = dataSorter->GetNodeEnd(rank) - dataSorter->GetNodeBegin(rank) + i + 1;
      int owning_rank;
      unsigned long node_number;
      node_partitioner.GetOwningRankAndNodeNumber(sorted_halo_nodes[i], owning_rank, node_number);
      neighbor_partitions[i] = owning_rank + 1; /* Partition numbers are 1-based. */
      neighbor_nodes[i] = static_cast<int64_t>(node_number);
    }
    err = tecFEPartitionCreate64(file_handle, zone, rank + 1, partition_num_nodes, partition_num_cells,
      static_cast<int64_t>(num_halo_nodes), &halo_node_local_numbers[0], &neighbor_partitions[0], &neighbor_nodes[0], 0, NULL);
    if (err) cout << rank << ": Error creating Tecplot zone partition." << endl;

    /* Gather halo node data. First, tell each rank how many nodes' worth of data we need from them. */
    for (size_t i = 0; i < num_halo_nodes; ++i)
      ++num_nodes_to_receive[neighbor_partitions[i] - 1];
    vector<int> num_nodes_to_send(size);
    SU2_MPI::Alltoall(&num_nodes_to_receive[0], 1, MPI_INT, &num_nodes_to_send[0], 1, MPI_INT, MPI_COMM_WORLD);

    /* Now send the global node numbers whose data we need,
       and receive the same from all other ranks.
       Each rank has globally consecutive node numbers,
       so we can just parcel out sorted_halo_nodes for send. */
    vector<int> nodes_to_send_displacements(size);
    vector<int> nodes_to_receive_displacements(size);
    nodes_to_send_displacements[0] = 0;
    nodes_to_receive_displacements[0] = 0;
    for(int iRank = 1; iRank < size; ++iRank) {
      nodes_to_send_displacements[iRank] = nodes_to_send_displacements[iRank - 1] + num_nodes_to_send[iRank - 1];
      nodes_to_receive_displacements[iRank] = nodes_to_receive_displacements[iRank - 1] + num_nodes_to_receive[iRank - 1];
    }
    int total_num_nodes_to_send = nodes_to_send_displacements[size - 1] + num_nodes_to_send[size - 1];
    vector<unsigned long> nodes_to_send(max(1, total_num_nodes_to_send));

    /* The terminology gets a bit confusing here. We're sending the node numbers
       (sorted_halo_nodes) whose data we need to receive, and receiving
       lists of nodes whose data we need to send. */
    if (sorted_halo_nodes.empty()) sorted_halo_nodes.resize(1); /* Avoid crash. */
    SU2_MPI::Alltoallv(&sorted_halo_nodes[0], &num_nodes_to_receive[0], &nodes_to_receive_displacements[0], MPI_UNSIGNED_LONG,
                       &nodes_to_send[0],     &num_nodes_to_send[0],    &nodes_to_send_displacements[0],    MPI_UNSIGNED_LONG,
                       MPI_COMM_WORLD);
    
    /* Now actually send and receive the data */
    vector<passivedouble> data_to_send(max(1, total_num_nodes_to_send * (int)fieldnames.size()));
    halo_var_data.resize(max((size_t)1, fieldnames.size() * num_halo_nodes));
    vector<int> num_values_to_send(size);
    vector<int> values_to_send_displacements(size);
    vector<int> num_values_to_receive(size);
    size_t index = 0;
    for(int iRank = 0; iRank < size; ++iRank) {
      /* We send and receive GlobalField_Counter values per node. */
      num_values_to_send[iRank]              = num_nodes_to_send[iRank] * fieldnames.size();
      values_to_send_displacements[iRank]    = nodes_to_send_displacements[iRank] * fieldnames.size();
      num_values_to_receive[iRank]           = num_nodes_to_receive[iRank] * fieldnames.size();
      values_to_receive_displacements[iRank] = nodes_to_receive_displacements[iRank] * fieldnames.size();
      for(iVar = 0; iVar < fieldnames.size(); ++iVar)
        for(int iNode = 0; iNode < num_nodes_to_send[iRank]; ++iNode) {
          unsigned long node_offset = nodes_to_send[nodes_to_send_displacements[iRank] + iNode] - dataSorter->GetNodeBegin(rank) - 1;
          data_to_send[index++] =dataSorter->GetData(iVar,node_offset);
        }
    }
    SU2_MPI::Alltoallv(&data_to_send[0],  &num_values_to_send[0],    &values_to_send_displacements[0],    MPI_DOUBLE,
                       &halo_var_data[0], &num_values_to_receive[0], &values_to_receive_displacements[0], MPI_DOUBLE,
                       MPI_COMM_WORLD);
  }
  else {
    /* Zone will be gathered to and output by MASTER_NODE */
    int32_t partition_owner = MASTER_NODE;
    err = tecZoneMapPartitionsToMPIRanks(file_handle, zone, 1, &partition_owner);
  }

  /*--- Write surface and volumetric solution data. ---*/
  
  if (zone_type == ZONETYPE_FEBRICK) {
    std::vector<passivedouble> values_to_write(dataSorter->GetnPoints());
    for (iVar = 0; err == 0 && iVar < fieldnames.size(); iVar++) {
      for(unsigned long i = 0; i < dataSorter->GetnPoints(); ++i)
        values_to_write[i] = dataSorter->GetData(iVar, i);
      err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, rank + 1, dataSorter->GetnPoints(), &values_to_write[0]);
      if (err) cout << rank << ": Error outputting Tecplot variable values." << endl;
      for (int iRank = 0; err == 0 && iRank < size; ++iRank) {
        if (num_nodes_to_receive[iRank] > 0) {
          int var_data_offset = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank] * iVar;
          err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, rank + 1, static_cast<int64_t>(num_nodes_to_receive[iRank]), &halo_var_data[var_data_offset]);
          if (err) cout << rank << ": Error outputting Tecplot halo values." << endl;
        }
      }
    }
  } else {
    if (rank == MASTER_NODE) {
      vector<passivedouble> var_data;
      unsigned long nPoint = dataSorter->GetnPoints();
      vector<unsigned long> num_points(size);
      SU2_MPI::Gather(&nPoint, 1, MPI_UNSIGNED_LONG, &num_points[0], 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      
      for(int iRank = 0; iRank < size; ++iRank) {
        int64_t rank_num_points = num_points[iRank];

        if (rank_num_points > 0) {
          if (iRank == rank) { /* Output local data. */
            std::vector<passivedouble> values_to_write;
            for (iVar = 0; err == 0 && iVar < fieldnames.size(); iVar++) {
              values_to_write.resize(rank_num_points);
              for(unsigned long i = 0; i < (unsigned long)rank_num_points; ++i)
                values_to_write[i] = dataSorter->GetData(iVar,i);
              err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, rank_num_points, &values_to_write[0]); 
              if (err) cout << rank << ": Error outputting Tecplot variable values." << endl;
            }
          }
          else { /* Receive data from other rank. */
            var_data.resize(max((int64_t)1, (int64_t)fieldnames.size() * rank_num_points));
            CBaseMPIWrapper::Recv(&var_data[0], fieldnames.size() * rank_num_points, MPI_DOUBLE, iRank, iRank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (iVar = 0; err == 0 && iVar < fieldnames.size(); iVar++) {
              err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, rank_num_points, &var_data[iVar * rank_num_points]);
              if (err) cout << rank << ": Error outputting Tecplot surface variable values." << endl;
            }
          }
        }
      }
    }
    else { /* Send data to MASTER_NODE */
      unsigned long nPoint = dataSorter->GetnPoints();
      
      SU2_MPI::Gather(&nPoint, 1, MPI_UNSIGNED_LONG, NULL, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    
      vector<passivedouble> var_data;
      size_t var_data_size = fieldnames.size() * dataSorter->GetnPoints();
      var_data.reserve(var_data_size);
      for (iVar = 0; err == 0 && iVar < fieldnames.size() ; iVar++)
          for(unsigned long i = 0; i < dataSorter->GetnPoints(); ++i)
            var_data.push_back(dataSorter->GetData(iVar,i));
      
      if (var_data.size() > 0)
        CBaseMPIWrapper::Send(&var_data[0], static_cast<int>(var_data.size()), MPI_DOUBLE, MASTER_NODE, rank, MPI_COMM_WORLD);
    }
  }

#else

  unsigned short iVar;
  
  vector<passivedouble> var_data;
  size_t var_data_size = fieldnames.size() * dataSorter->GetnPoints();
  var_data.reserve(var_data_size);
  

  for (iVar = 0; err == 0 && iVar <  fieldnames.size(); iVar++) {
    for(unsigned long i = 0; i < dataSorter->GetnPoints(); ++i)
      var_data.push_back(dataSorter->GetData(iVar,i));
    err = tecZoneVarWriteDoubleValues(file_handle, zone, iVar + 1, 0, dataSorter->GetnPoints(), &var_data[iVar * dataSorter->GetnPoints()]);
    if (err) cout << rank << ": Error outputting Tecplot variable value." << endl;
  }
  

#endif /* HAVE_MPI */
  
  /*--- Write connectivity data. ---*/

  unsigned long iElem;
  
#ifdef HAVE_MPI
  if (zone_type == ZONETYPE_FEBRICK) {

    int64_t nodes[8];

    /**
     *  Each rank writes node numbers relative to the partition it is outputting (starting with node number 1).
     *  Ghost (halo) nodes identified above are numbered sequentially just beyond the end of the actual, local nodes.
     *  Note that beg_node and end_node refer to zero-based node numbering, but Conn_* contain one-based node numbers.
     */
#define MAKE_LOCAL(n) dataSorter->GetNodeBegin(rank) < (unsigned long)n && (unsigned long)n <= dataSorter->GetNodeEnd(rank) \
  ? (int64_t)((unsigned long)n - dataSorter->GetNodeBegin(rank)) \
  : GetHaloNodeNumber(n, dataSorter->GetNodeEnd(rank) - dataSorter->GetNodeBegin(rank), sorted_halo_nodes)

    for (iElem = 0; err == 0 && iElem < nParallel_Tetr; iElem++) {
      nodes[0] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 0));
      nodes[1] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 1));
      nodes[2] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 2));
      nodes[3] = nodes[2];
      nodes[4] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3));
      nodes[5] = nodes[4];
      nodes[6] = nodes[4];
      nodes[7] = nodes[4];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank + 1, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }

    for (iElem = 0; err == 0 && iElem < nParallel_Hexa; iElem++) {
      nodes[0] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 0));
      nodes[1] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 1));
      nodes[2] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 2));
      nodes[3] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 3));
      nodes[4] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 4));
      nodes[5] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 5));
      nodes[6] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 6));
      nodes[7] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 7));
      err = tecZoneNodeMapWrite64(file_handle, zone, rank + 1, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
      
    for (iElem = 0; err == 0 && iElem < nParallel_Pris; iElem++) {
      nodes[0] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(PRISM, iElem, 0));
      nodes[1] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(PRISM, iElem, 1));
      nodes[2] = nodes[1];
      nodes[3] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(PRISM, iElem, 2));
      nodes[4] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(PRISM, iElem, 3));
      nodes[5] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(PRISM, iElem, 4));
      nodes[6] = nodes[5];
      nodes[7] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(PRISM, iElem, 5));
      err = tecZoneNodeMapWrite64(file_handle, zone, rank + 1, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
    
    for (iElem = 0; err == 0 && iElem < nParallel_Pyra; iElem++) {
      nodes[0] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(PYRAMID, iElem, 0));
      nodes[1] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(PYRAMID, iElem, 1));
      nodes[2] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(PYRAMID, iElem, 2));
      nodes[3] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(PYRAMID, iElem, 3));
      nodes[4] = MAKE_LOCAL(dataSorter->GetElem_Connectivity(PYRAMID, iElem, 4));
      nodes[5] = nodes[4];
      nodes[6] = nodes[4];
      nodes[7] = nodes[4];
      err = tecZoneNodeMapWrite64(file_handle, zone, rank + 1, 1, 8, nodes);
      if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
    }
  } else {
    if (rank == MASTER_NODE) {

      /* Non-hexahedral output by the master node. Output local data directly, and gather other data from the other ranks. */

      int64_t nodes[4];

      vector<unsigned long> connectivity_sizes(size);
      unsigned long unused = 0;
      SU2_MPI::Gather(&unused, 1, MPI_UNSIGNED_LONG, &connectivity_sizes[0], 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      vector<int64_t> connectivity;
      for(int iRank = 0; iRank < size; ++iRank) {
        if (iRank == rank) {
          for (iElem = 0; err == 0 && iElem < nParallel_Line; iElem++) {
            nodes[0] = dataSorter->GetElem_Connectivity(LINE, iElem, 0);
            nodes[1] = dataSorter->GetElem_Connectivity(LINE, iElem, 1);
            err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 2, nodes);
            if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
          }
          
          for (iElem = 0; err == 0 && iElem < nParallel_Tria; iElem++) {
            nodes[0] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 0);
            nodes[1] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 1);
            nodes[2] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 2);
            nodes[3] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 2);
            err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 4, nodes);
            if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
          }
          
          for (iElem = 0; err == 0 && iElem < nParallel_Quad; iElem++) {
            nodes[0] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 0);
            nodes[1] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 1);
            nodes[2] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 2);
            nodes[3] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 3);
            err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, 4, nodes);
            if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
          }
          
        } else { /* Receive node map and write out. */
          connectivity.resize(max((unsigned long)1, connectivity_sizes[iRank]));
          SU2_MPI::Recv(&connectivity[0], connectivity_sizes[iRank], MPI_UNSIGNED_LONG, iRank, iRank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          err = tecZoneNodeMapWrite64(file_handle, zone, 0, 1, connectivity_sizes[iRank], &connectivity[0]);
          if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
        }
      }
    } else {

      /* Non-hexahedral output by non-master node. Send what we've got to the master node. */

      unsigned long connectivity_size;
      connectivity_size = 2 * nParallel_Line + 4 * (nParallel_Tria + nParallel_Quad);
      SU2_MPI::Gather(&connectivity_size, 1, MPI_UNSIGNED_LONG, NULL, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      vector<int64_t> connectivity;
      connectivity.reserve(connectivity_size);
      for (iElem = 0; err == 0 && iElem < nParallel_Line; iElem++) {
        connectivity.push_back(dataSorter->GetElem_Connectivity(LINE, iElem, 0));
        connectivity.push_back(dataSorter->GetElem_Connectivity(LINE, iElem, 0));
      }
      
      for (iElem = 0; err == 0 && iElem < nParallel_Tria; iElem++) {
        connectivity.push_back(dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 0));
        connectivity.push_back(dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 1));
        connectivity.push_back(dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 2));
        connectivity.push_back(dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 2));
      }
      
      for (iElem = 0; err == 0 && iElem < nParallel_Quad; iElem++) {
        connectivity.push_back(dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 0));
        connectivity.push_back(dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 1));
        connectivity.push_back(dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 2));
        connectivity.push_back(dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 3));
      }
      
      if (connectivity.empty()) connectivity.resize(1); /* Avoid crash */
      SU2_MPI::Send(&connectivity[0], connectivity_size, MPI_UNSIGNED_LONG, MASTER_NODE, rank, MPI_COMM_WORLD);
    }
  }
#else
  
  int64_t nodes[8];
  
  for (iElem = 0; err == 0 && iElem < nParallel_Tria; iElem++) {
    nodes[0] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 0);
    nodes[1] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 1);
    nodes[2] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 2);
    nodes[3] = dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 2);
    err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 4, nodes);
    if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
  }
  
  for (iElem = 0; err == 0 && iElem < nParallel_Quad; iElem++) {
    nodes[0] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 0);
    nodes[1] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 1);
    nodes[2] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 2);
    nodes[3] = dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 3);
    err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 4, nodes);
    if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
  }
  
  for (iElem = 0; err == 0 && iElem < nParallel_Tetr; iElem++) {
    nodes[0] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 0);
    nodes[1] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 1);
    nodes[2] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 2);
    nodes[3] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 2);
    nodes[4] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3);
    nodes[5] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3);
    nodes[6] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3);
    nodes[7] = dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3);
    err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 8, nodes);
    if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
  }
  
  for (iElem = 0; err == 0 && iElem < nParallel_Hexa; iElem++) {
    nodes[0] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 0);
    nodes[1] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 1);
    nodes[2] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 2);
    nodes[3] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 3);
    nodes[4] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 4);
    nodes[5] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 5);
    nodes[6] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 6);
    nodes[7] = dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 7);
    err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 8, nodes);
    if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
  }
  
  for (iElem = 0; err == 0 && iElem < nParallel_Pris; iElem++) {
    nodes[0] = dataSorter->GetElem_Connectivity(PRISM, iElem, 0);
    nodes[1] = dataSorter->GetElem_Connectivity(PRISM, iElem, 1);
    nodes[2] = dataSorter->GetElem_Connectivity(PRISM, iElem, 1);
    nodes[3] = dataSorter->GetElem_Connectivity(PRISM, iElem, 2);
    nodes[4] = dataSorter->GetElem_Connectivity(PRISM, iElem, 3);
    nodes[5] = dataSorter->GetElem_Connectivity(PRISM, iElem, 4);
    nodes[6] = dataSorter->GetElem_Connectivity(PRISM, iElem, 4);
    nodes[7] = dataSorter->GetElem_Connectivity(PRISM, iElem, 5);
    err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 8, nodes);
    if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
  }
  
  for (iElem = 0; err == 0 && iElem < nParallel_Pyra; iElem++) {
    nodes[0] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 0);
    nodes[1] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 1);
    nodes[2] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 2);
    nodes[3] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 3);
    nodes[4] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 4);
    nodes[5] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 4);
    nodes[6] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 4);
    nodes[7] = dataSorter->GetElem_Connectivity(PYRAMID, iElem, 4);
    err = tecZoneNodeMapWrite64(file_handle, zone, rank, 1, 8, nodes);
    if (err) cout << rank << ": Error outputting Tecplot node values." << endl;
  }
  
  

#endif
  
  err = tecFileWriterClose(&file_handle);
  if (err) cout << rank << ": Error finishing Tecplot file output." << endl;
  
#endif /* HAVE_TECIO */
  
  /*--- Compute and store the write time. ---*/
  
#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  UsedTime = StopTime-StartTime;

  file_size = Determine_Filesize(fileName);
  
  /*--- Compute and store the bandwidth ---*/
  
  Bandwidth = file_size/(1.0e6)/UsedTime;
}


int64_t CTecplotBinaryFileWriter::GetHaloNodeNumber(unsigned long global_node_number, unsigned long last_local_node, vector<unsigned long> const &halo_node_list)
{
  vector<unsigned long>::const_iterator it = lower_bound(halo_node_list.begin(), halo_node_list.end(), global_node_number);
  assert(it != halo_node_list.end());
  assert(*it == global_node_number);
  /* When C++11 is universally available, replace the following mouthful with "auto" */
  iterator_traits<vector<unsigned long>::const_iterator>::difference_type offset = distance(halo_node_list.begin(), it);
  assert(offset >= 0);
  return (int64_t)(last_local_node + offset + 1);
}

