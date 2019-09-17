#include "../../../include/output/filewriter/CCatalystWriter.hpp"
#include "../../../include/output/filewriter/CParallelDataSorter.hpp"
#include "../../../Common/include/toolboxes/SU2_LOG.hpp"
#include <set>


CCatalystWriter::CCatalystWriter(vector<string> fields, unsigned short nDim) : CFileWriter(fields, nDim){
  
  Processor = vtkCPProcessor::New();
  Processor->Initialize();
  vtkNew<vtkCPPythonScriptPipeline> pipeline;
  pipeline->Initialize("Test.py");
  Processor->AddPipeline(pipeline.GetPointer());
  VTKGrid = NULL;
}

CCatalystWriter::~CCatalystWriter(){}

void CCatalystWriter::Write_Data(unsigned long TimeStep, double time, CParallelDataSorter *data_sorter){
  LOG_SCOPE_FUNCTION(INFO);
  SU2_INFO << "Time Step: " << TimeStep << ", Time: " << time; 
  vtkNew<vtkCPDataDescription> dataDescription;
  dataDescription->AddInput("input");
  dataDescription->SetTimeData(time, TimeStep);
  if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0){
    vtkCPInputDataDescription* idd = dataDescription->GetInputDescriptionByName("input");
    BuildVTKDataStructures(idd, data_sorter);
    idd->SetGrid(VTKGrid);
    Processor->CoProcess(dataDescription.GetPointer());
  }
}

void CCatalystWriter::BuildVTKGrid(CParallelDataSorter *data_sorter)
{
  LOG_SCOPE_FUNCTION(INFO);
  
  const int NCOORDS = 3;
  
  /*--- Compute our local number of elements, the required storage,
   and reduce the total number of elements and storage globally. ---*/
  
  unsigned long nTot_Line;
  unsigned long nTot_Tria, nTot_Quad;
  unsigned long nTot_Tetr, nTot_Hexa, nTot_Pris, nTot_Pyra;
  unsigned long myElem, myElemStorage, GlobalElem, GlobalElemStorage;
  
  unsigned long nParallel_Line = data_sorter->GetnElem(LINE),
                nParallel_Tria = data_sorter->GetnElem(TRIANGLE),
                nParallel_Quad = data_sorter->GetnElem(QUADRILATERAL),
                nParallel_Tetr = data_sorter->GetnElem(TETRAHEDRON),
                nParallel_Hexa = data_sorter->GetnElem(HEXAHEDRON),
                nParallel_Pris = data_sorter->GetnElem(PRISM),
                nParallel_Pyra = data_sorter->GetnElem(PYRAMID);
  
  SU2_MPI::Allreduce(&nParallel_Line, &nTot_Line, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Tria, &nTot_Tria, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Quad, &nTot_Quad, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Tetr, &nTot_Tetr, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Hexa, &nTot_Hexa, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Pris, &nTot_Pris, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nParallel_Pyra, &nTot_Pyra, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  
  myElem        = (nParallel_Line + nParallel_Tria + nParallel_Quad + nParallel_Tetr +
                   nParallel_Hexa + nParallel_Pris + nParallel_Pyra);
  myElemStorage = (nParallel_Line*3 + nParallel_Tria*4 + nParallel_Quad*5 + nParallel_Tetr*5 +
                   nParallel_Hexa*9 + nParallel_Pris*7 + nParallel_Pyra*6);
  
  GlobalElem        = (nTot_Line + nTot_Tria   + nTot_Quad   + nTot_Tetr   +
                       nTot_Hexa   + nTot_Pris   + nTot_Pyra);
  GlobalElemStorage = (nTot_Line*3 + nTot_Tria*4 + nTot_Quad*5 + nTot_Tetr*5 +
                       nTot_Hexa*9 + nTot_Pris*7 + nTot_Pyra*6);
  

  
  
  unsigned short iVar;
//  NodePartitioner node_partitioner(num_nodes, size);
  vector<passivedouble> halo_var_data;
  vector<int> num_nodes_to_receive(size, 0);
  vector<int> values_to_receive_displacements(size);
  
  halo_nodes.clear();
  sorted_halo_nodes.clear();
  
  /* We output a single, partitioned zone where each rank outputs one partition. */
  vector<int32_t> partition_owners;
  partition_owners.reserve(size);
  for (int32_t iRank = 0; iRank < size; ++iRank)
    partition_owners.push_back(iRank);

  /* Gather a list of nodes we refer to but are not outputting. */

  for (unsigned long i = 0; i < nParallel_Tria * N_POINTS_TRIANGLE; ++i)
    if (data_sorter->FindProcessor(data_sorter->GetElem_Connectivity(TRIANGLE, 0, i) -1) != rank)
      halo_nodes.insert(data_sorter->GetElem_Connectivity(TRIANGLE, 0, i)-1);

  for (unsigned long i = 0; i < nParallel_Quad * N_POINTS_QUADRILATERAL; ++i)
    if (data_sorter->FindProcessor(data_sorter->GetElem_Connectivity(QUADRILATERAL, 0, i) -1) != rank)
      halo_nodes.insert(data_sorter->GetElem_Connectivity(QUADRILATERAL, 0, i)-1);

  for (unsigned long i = 0; i < nParallel_Tetr * N_POINTS_TETRAHEDRON; ++i)
    if ((unsigned long)data_sorter->GetElem_Connectivity(TETRAHEDRON, 0, i)-1 < data_sorter->GetNodeBegin(rank) || 
        data_sorter->GetNodeEnd(rank) <= (unsigned long)data_sorter->GetElem_Connectivity(TETRAHEDRON, 0, i)-1)
      halo_nodes.insert(data_sorter->GetElem_Connectivity(TETRAHEDRON, 0, i)-1);

  for (unsigned long i = 0; i < nParallel_Hexa * N_POINTS_HEXAHEDRON; ++i)
    if ((unsigned long)data_sorter->GetElem_Connectivity(HEXAHEDRON, 0, i)-1 < data_sorter->GetNodeBegin(rank) || 
        data_sorter->GetNodeEnd(rank) <= (unsigned long)data_sorter->GetElem_Connectivity(HEXAHEDRON, 0, i)-1)
      halo_nodes.insert(data_sorter->GetElem_Connectivity(HEXAHEDRON, 0, i)-1);
    
  for (unsigned long i = 0; i < nParallel_Pris * N_POINTS_PRISM; ++i)
    if ((unsigned long)data_sorter->GetElem_Connectivity(PRISM, 0, i)-1 < data_sorter->GetNodeBegin(rank) || 
        data_sorter->GetNodeEnd(rank) <= (unsigned long)data_sorter->GetElem_Connectivity(PRISM, 0, i)-1)
      halo_nodes.insert(data_sorter->GetElem_Connectivity(PRISM, 0, i)-1);
  
  for (unsigned long i = 0; i < nParallel_Pyra * N_POINTS_PYRAMID; ++i)
    if ((unsigned long)data_sorter->GetElem_Connectivity(PYRAMID, 0, i)-1 < data_sorter->GetNodeBegin(rank) || 
        data_sorter->GetNodeEnd(rank) <= (unsigned long)data_sorter->GetElem_Connectivity(PYRAMID, 0, i)-1)
      halo_nodes.insert(data_sorter->GetElem_Connectivity(PYRAMID, 0, i)-1);

  /* Sorted list of halo nodes for this MPI rank. */
  sorted_halo_nodes.assign(halo_nodes.begin(), halo_nodes.end());
      
  /* Have to include all nodes our cells refer to or TecIO will barf, so add the halo node count to the number of local nodes. */
  int64_t partition_num_nodes = data_sorter->GetNodeEnd(rank) - data_sorter->GetNodeBegin(rank) + static_cast<int64_t>(halo_nodes.size());
  int64_t partition_num_cells = nParallel_Tetr + nParallel_Hexa + nParallel_Pris + nParallel_Pyra;

  /*--- We effectively tack the halo nodes onto the end of the node list for this partition.
    TecIO will later replace them with references to nodes in neighboring partitions. */
  size_t num_halo_nodes = sorted_halo_nodes.size();
  vector<int64_t> halo_node_local_numbers(max((size_t)1, num_halo_nodes)); /* Min size 1 to avoid crashes when we access these vectors below. */
  vector<int32_t> neighbor_partitions(max((size_t)1, num_halo_nodes));
  vector<int64_t> neighbor_nodes(max((size_t)1, num_halo_nodes));
  for(int64_t i = 0; i < static_cast<int64_t>(num_halo_nodes); ++i) {
    halo_node_local_numbers[i] = data_sorter->GetNodeEnd(rank) - data_sorter->GetNodeBegin(rank) + i;
    int owning_rank;
    owning_rank = data_sorter->FindProcessor(sorted_halo_nodes[i]);
    unsigned long node_number = sorted_halo_nodes[i] - data_sorter->GetNodeBegin(owning_rank);
    neighbor_partitions[i] = owning_rank + 1; /* Partition numbers are 1-based. */
    neighbor_nodes[i] = static_cast<int64_t>(node_number);
  }

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
        unsigned long node_offset = nodes_to_send[nodes_to_send_displacements[iRank] + iNode] - data_sorter->GetNodeBegin(rank);
        data_to_send[index++] = SU2_TYPE::GetValue(data_sorter->GetData(iVar,node_offset));
      }
  }
  SU2_MPI::Alltoallv(&data_to_send[0],  &num_values_to_send[0],    &values_to_send_displacements[0],    MPI_DOUBLE,
                     &halo_var_data[0], &num_values_to_receive[0], &values_to_receive_displacements[0], MPI_DOUBLE,
                     MPI_COMM_WORLD);
  

  int additional_nodes = 0;
  for (int iRank = 0; iRank < size; ++iRank) {
    additional_nodes += num_nodes_to_receive[iRank];
  }
  for (int i = 0; i < sorted_halo_nodes.size(); i++){
    cout << sorted_halo_nodes[i] << " " <<  data_sorter->GetNodeEnd(rank) << " " <<  data_sorter->GetNodeBegin(rank) << endl;
  }
  // create the points information
  vtkNew<vtkDoubleArray> pointArray;
  pointArray->SetNumberOfComponents(3);
  pointArray->SetNumberOfValues((data_sorter->GetnPoints()+num_halo_nodes)*NCOORDS);
  int id = 0;
  for (int iRank = 0; iRank < size; ++iRank) {
    for (int i = 0; i < num_nodes_to_receive[iRank]; i++){  
      int displ = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank]*0; 
      pointArray->SetValue((data_sorter->GetnPoints() + id)*NCOORDS+0, halo_var_data[displ+i]);
      id++;    
    }
  }
  id = 0;
  for (int iRank = 0; iRank < size; ++iRank) {
    for (int i = 0; i < num_nodes_to_receive[iRank]; i++){  
      int displ = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank]*1; 
      pointArray->SetValue((data_sorter->GetnPoints() + id)*NCOORDS+1, halo_var_data[displ+i]);
      id++;    
    }
  }
  id = 0;
  for (int iRank = 0; iRank < size; ++iRank) {
    for (int i = 0; i < num_nodes_to_receive[iRank]; i++){  
      int displ = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank]*2; 
      if (nDim == 3) pointArray->SetValue((data_sorter->GetnPoints() + id)*NCOORDS+2, halo_var_data[displ+i]);
      else {pointArray->SetValue((data_sorter->GetnPoints() + id)*NCOORDS+2, 0.0);}
      id++;    
    }
  }
  
  for (unsigned long iPoint = 0; iPoint < data_sorter->GetnPoints(); iPoint++){
    pointArray->SetValue(iPoint*NCOORDS+0, data_sorter->GetData(0,iPoint));
    pointArray->SetValue(iPoint*NCOORDS+1, data_sorter->GetData(1,iPoint));
    if (nDim == 3 ) pointArray->SetValue(iPoint*NCOORDS+2, data_sorter->GetData(2,iPoint));
    else pointArray->SetValue(iPoint*NCOORDS+2, 0);
  }
 
    
  vtkNew<vtkPoints> points;
  points->SetData(pointArray.GetPointer());
  VTKGrid->SetPoints(points.GetPointer());

  // create the cells
  VTKGrid->Allocate(static_cast<vtkIdType>(myElemStorage));
  
#define MAKE_LOCAL(n) data_sorter->FindProcessor(n) == rank \
  ? (int64_t)((unsigned long)n - data_sorter->GetNodeBegin(rank)) \
  : GetHaloNodeNumber(n, data_sorter->GetNodeEnd(rank) - data_sorter->GetNodeBegin(rank), sorted_halo_nodes)

#define ISN_LOCAL(n) (data_sorter->FindProcessor(n) != rank)
  
  for (unsigned long iElem = 0; iElem < nParallel_Line; iElem++) {
    vtkIdType tmp[2] = {(vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(LINE, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(LINE, iElem, 1)-1))};
    VTKGrid->InsertNextCell(VTK_LINE, N_POINTS_LINE, tmp);    
  }
  
  
  for (unsigned long iElem = 0; iElem < nParallel_Tria; iElem++) {    
    vtkIdType tmp[3] = {(vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(TRIANGLE, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(TRIANGLE, iElem, 1)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(TRIANGLE, iElem, 2)-1))};
    VTKGrid->InsertNextCell(VTK_TRIANGLE, N_POINTS_TRIANGLE, tmp);   
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Quad; iElem++) {
    vtkIdType tmp[4] = {(vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(QUADRILATERAL, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(QUADRILATERAL, iElem, 1)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(QUADRILATERAL, iElem, 2)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(QUADRILATERAL, iElem, 3)-1))};
    VTKGrid->InsertNextCell(VTK_QUAD, N_POINTS_QUADRILATERAL, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Tetr; iElem++) {
    vtkIdType tmp[4] = {(vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 1)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 2)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3)-1))};
    VTKGrid->InsertNextCell(VTK_TETRA, N_POINTS_TETRAHEDRON, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Hexa; iElem++) {
    vtkIdType tmp[8] = {(vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 1)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 2)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 3)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 4)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 5)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 6)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 7)-1))};
    VTKGrid->InsertNextCell(VTK_HEXAHEDRON, N_POINTS_HEXAHEDRON, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Pris; iElem++) {
    vtkIdType tmp[6] = {(vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(PRISM, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(PRISM, iElem, 1)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(PRISM, iElem, 2)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(PRISM, iElem, 3)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(PRISM, iElem, 4)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(PRISM, iElem, 5)-1))};
    VTKGrid->InsertNextCell(VTK_WEDGE, N_POINTS_PRISM, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Pyra; iElem++) {
    vtkIdType tmp[5] = {(vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(PYRAMID, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(PYRAMID, iElem, 1)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(PYRAMID, iElem, 2)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(PYRAMID, iElem, 3)-1)),
                        (vtkIdType)(MAKE_LOCAL(data_sorter->GetElem_Connectivity(PYRAMID, iElem, 4)-1))};
    VTKGrid->InsertNextCell(VTK_PYRAMID, N_POINTS_PYRAMID, tmp);    
  }
}


void CCatalystWriter::UpdateVTKAttributes(vtkCPInputDataDescription* idd, CParallelDataSorter *data_sorter)
{
  LOG_SCOPE_FUNCTION(INFO);
  
  ERROR_CONTEXT("Number of points", data_sorter->GetnPoints());
  
  bool allocate = false;
  
  if (VTKGrid->GetPointData()->GetNumberOfArrays() == 0){
    
    allocate = true;
    
  }  
  unsigned short iVar;
//  NodePartitioner node_partitioner(num_nodes, size);
  vector<passivedouble> halo_var_data;
  vector<int> num_nodes_to_receive(size, 0);
  vector<int> values_to_receive_displacements(size);
  
  /* We output a single, partitioned zone where each rank outputs one partition. */
  vector<int32_t> partition_owners;
  partition_owners.reserve(size);
  for (int32_t iRank = 0; iRank < size; ++iRank)
    partition_owners.push_back(iRank);

  /*--- We effectively tack the halo nodes onto the end of the node list for this partition.
    TecIO will later replace them with references to nodes in neighboring partitions. */
  size_t num_halo_nodes = sorted_halo_nodes.size();
  vector<int64_t> halo_node_local_numbers(max((size_t)1, num_halo_nodes)); /* Min size 1 to avoid crashes when we access these vectors below. */
  vector<int32_t> neighbor_partitions(max((size_t)1, num_halo_nodes));
  vector<int64_t> neighbor_nodes(max((size_t)1, num_halo_nodes));
  for(int64_t i = 0; i < static_cast<int64_t>(num_halo_nodes); ++i) {
    halo_node_local_numbers[i] = data_sorter->GetNodeEnd(rank) - data_sorter->GetNodeBegin(rank) + i;
    int owning_rank;
    owning_rank = data_sorter->FindProcessor(sorted_halo_nodes[i]);
    unsigned long node_number = sorted_halo_nodes[i] - data_sorter->GetNodeBegin(owning_rank);
    neighbor_partitions[i] = owning_rank + 1; /* Partition numbers are 1-based. */
    neighbor_nodes[i] = static_cast<int64_t>(node_number);
  }

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
        unsigned long node_offset = nodes_to_send[nodes_to_send_displacements[iRank] + iNode] - data_sorter->GetNodeBegin(rank);
        data_to_send[index++] = SU2_TYPE::GetValue(data_sorter->GetData(iVar,node_offset));
      }
  }
  SU2_MPI::Alltoallv(&data_to_send[0],  &num_values_to_send[0],    &values_to_send_displacements[0],    MPI_DOUBLE,
                     &halo_var_data[0], &num_values_to_receive[0], &values_to_receive_displacements[0], MPI_DOUBLE,
                     MPI_COMM_WORLD);
  
  
  for (unsigned short iField = nDim; iField < fieldnames.size(); iField++){
    
    if (idd->IsFieldNeeded(fieldnames[iField].c_str(), vtkDataObject::POINT) == true){
      if (allocate){
        
        vtkNew<vtkDoubleArray> array;
        array->SetName(fieldnames[iField].c_str());
        array->SetNumberOfComponents(1);
        array->SetNumberOfValues(static_cast<vtkIdType>(data_sorter->GetnPoints()+num_halo_nodes));
        
        VTKGrid->GetPointData()->AddArray(array.GetPointer());
        
        SU2_INFO << "Adding catalyst output " << fieldnames[iField];
      }
      
      SU2_INFO << "Setting catalyst output " << fieldnames[iField];
      
      vtkDoubleArray* data =
          vtkDoubleArray::SafeDownCast(VTKGrid->GetPointData()->GetArray(fieldnames[iField].c_str()));
      //      data->SetArray(data_sorter->GetData(iField), static_cast<vtkIdType>(data_sorter->GetnPoints()), 1);
      for (vtkIdType i = 0; i < data_sorter->GetnPoints(); i++)
      {
        su2double val = 0.0;
          val = data_sorter->GetData(iField, i);
        data->SetValue(i, val);
      }
      int id = 0;
      for (int iRank = 0; iRank < size; ++iRank) {
        for (int i = 0; i < num_nodes_to_receive[iRank]; i++){  
          int displ = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank]*iField; 
          data->SetValue(data_sorter->GetnPoints() + id, halo_var_data[displ+i]);
          id++;    
        }
      }
    }
  }
}

void CCatalystWriter::BuildVTKDataStructures(vtkCPInputDataDescription* idd, CParallelDataSorter *data_sorter)
{
  LOG_SCOPE_FUNCTION(INFO);
  
  if (VTKGrid == NULL)
  {
    SU2_INFO << "Creating VTK Grid";
    // The grid structure isn't changing so we only build it
    // the first time it's needed. If we needed the memory
    // we could delete it and rebuild as necessary.
    VTKGrid = vtkUnstructuredGrid::New();
    BuildVTKGrid(data_sorter);
  }
  UpdateVTKAttributes(idd, data_sorter);
  
}




int64_t CCatalystWriter::GetHaloNodeNumber(unsigned long global_node_number, unsigned long last_local_node, vector<unsigned long> const &halo_node_list)
{
  vector<unsigned long>::const_iterator it = lower_bound(halo_node_list.begin(), halo_node_list.end(), global_node_number);
  assert(it != halo_node_list.end());
  assert(*it == global_node_number);
  /* When C++11 is universally available, replace the following mouthful with "auto" */
  iterator_traits<vector<unsigned long>::const_iterator>::difference_type offset = distance(halo_node_list.begin(), it);
  assert(offset >= 0);
  return (int64_t)(last_local_node + offset);
}
