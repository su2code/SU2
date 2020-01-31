#include "../../../include/output/filewriter/CCatalystWriter.hpp"
#include "../../../include/output/filewriter/CParallelDataSorter.hpp"
#include "../../../Common/include/toolboxes/SU2_LOG.hpp"
#include <set>


CCatalystWriter::CCatalystWriter(string fileName, CParallelDataSorter *dataSorter): 
  CFileWriter(std::move(fileName), dataSorter, ""){
  
  Processor = vtkCPProcessor::New();
  Processor->Initialize();
  vtkNew<vtkCPPythonScriptPipeline> pipeline;
  pipeline->Initialize("Test.py");
  Processor->AddPipeline(pipeline.GetPointer());
  VTKGrid = NULL;
}

CCatalystWriter::~CCatalystWriter(){}

void CCatalystWriter::Write_Data(unsigned long TimeStep, double time){
  LOG_SCOPE_FUNCTION(INFO);
  SU2_INFO << "Time Step: " << TimeStep << ", Time: " << time; 
  vtkNew<vtkCPDataDescription> dataDescription;
  dataDescription->AddInput("input");
  dataDescription->SetTimeData(time, TimeStep);
  if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0){
    vtkCPInputDataDescription* idd = dataDescription->GetInputDescriptionByName("input");
    BuildVTKDataStructures(idd, dataSorter);
    idd->SetGrid(VTKGrid);
    Processor->CoProcess(dataDescription.GetPointer());
  }
}

void CCatalystWriter::BuildVTKGrid()
{
  LOG_SCOPE_FUNCTION(INFO);
  
  const vector<string>& fieldNames = dataSorter->GetFieldNames();
  
  const int NCOORDS = 3;
  
  /*--- Compute our local number of elements, the required storage,
   and reduce the total number of elements and storage globally. ---*/

  unsigned long myElem, myElemStorage, GlobalElem, GlobalElemStorage;
  
  unsigned long nParallel_Line = dataSorter->GetnElem(LINE),
                nParallel_Tria = dataSorter->GetnElem(TRIANGLE),
                nParallel_Quad = dataSorter->GetnElem(QUADRILATERAL),
                nParallel_Tetr = dataSorter->GetnElem(TETRAHEDRON),
                nParallel_Hexa = dataSorter->GetnElem(HEXAHEDRON),
                nParallel_Pris = dataSorter->GetnElem(PRISM),
                nParallel_Pyra = dataSorter->GetnElem(PYRAMID);

  myElem            = dataSorter->GetnElem();
  myElemStorage     = dataSorter->GetnConn();
  GlobalElem        = dataSorter->GetnElemGlobal();
  GlobalElemStorage = dataSorter->GetnConnGlobal();
  
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

  for (unsigned long i = 0; i < nParallel_Tetr * N_POINTS_TETRAHEDRON; ++i)
    if ((unsigned long)dataSorter->GetElem_Connectivity(TETRAHEDRON, 0, i)-1 < dataSorter->GetNodeBegin(rank) || 
        dataSorter->GetNodeEnd(rank) <= (unsigned long)dataSorter->GetElem_Connectivity(TETRAHEDRON, 0, i)-1)
      halo_nodes.insert(dataSorter->GetElem_Connectivity(TETRAHEDRON, 0, i)-1);

  for (unsigned long i = 0; i < nParallel_Hexa * N_POINTS_HEXAHEDRON; ++i)
    if ((unsigned long)dataSorter->GetElem_Connectivity(HEXAHEDRON, 0, i)-1 < dataSorter->GetNodeBegin(rank) || 
        dataSorter->GetNodeEnd(rank) <= (unsigned long)dataSorter->GetElem_Connectivity(HEXAHEDRON, 0, i)-1)
      halo_nodes.insert(dataSorter->GetElem_Connectivity(HEXAHEDRON, 0, i)-1);
    
  for (unsigned long i = 0; i < nParallel_Pris * N_POINTS_PRISM; ++i)
    if ((unsigned long)dataSorter->GetElem_Connectivity(PRISM, 0, i)-1 < dataSorter->GetNodeBegin(rank) || 
        dataSorter->GetNodeEnd(rank) <= (unsigned long)dataSorter->GetElem_Connectivity(PRISM, 0, i)-1)
      halo_nodes.insert(dataSorter->GetElem_Connectivity(PRISM, 0, i)-1);
  
  for (unsigned long i = 0; i < nParallel_Pyra * N_POINTS_PYRAMID; ++i)
    if ((unsigned long)dataSorter->GetElem_Connectivity(PYRAMID, 0, i)-1 < dataSorter->GetNodeBegin(rank) || 
        dataSorter->GetNodeEnd(rank) <= (unsigned long)dataSorter->GetElem_Connectivity(PYRAMID, 0, i)-1)
      halo_nodes.insert(dataSorter->GetElem_Connectivity(PYRAMID, 0, i)-1);

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
    int owning_rank;
    owning_rank = dataSorter->FindProcessor(sorted_halo_nodes[i]);
    unsigned long node_number = sorted_halo_nodes[i] - dataSorter->GetNodeBegin(owning_rank);
    neighbor_partitions[i] = owning_rank + 1; /* Partition numbers are 1-based. */
    neighbor_nodes[i] = static_cast<int64_t>(node_number);
  }

  /* Gather halo node data. First, tell each rank how many nodes' worth of data we need from them. */
  for (size_t i = 0; i < num_halo_nodes; ++i)
    ++num_nodes_to_receive[neighbor_partitions[i] - 1];
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
  data_to_send.resize(max(1, total_num_nodes_to_send * (int)fieldNames.size()));
  halo_var_data.resize(max((size_t)1, fieldNames.size() * num_halo_nodes));
  num_values_to_send.resize(size);
  values_to_send_displacements.resize(size);
  num_values_to_receive.resize(size);
  size_t index = 0;
  for(int iRank = 0; iRank < size; ++iRank) {
    /* We send and receive GlobalField_Counter values per node. */
    num_values_to_send[iRank]              = num_nodes_to_send[iRank] * fieldNames.size();
    values_to_send_displacements[iRank]    = nodes_to_send_displacements[iRank] * fieldNames.size();
    num_values_to_receive[iRank]           = num_nodes_to_receive[iRank] * fieldNames.size();
    values_to_receive_displacements[iRank] = nodes_to_receive_displacements[iRank] * fieldNames.size();
    for(iVar = 0; iVar < fieldNames.size(); ++iVar)
      for(int iNode = 0; iNode < num_nodes_to_send[iRank]; ++iNode) {
        unsigned long node_offset = nodes_to_send[nodes_to_send_displacements[iRank] + iNode] - dataSorter->GetNodeBegin(rank);
        data_to_send[index++] = SU2_TYPE::GetValue(dataSorter->GetData(iVar,node_offset));
      }
  }
  SU2_MPI::Alltoallv(&data_to_send[0],  &num_values_to_send[0],    &values_to_send_displacements[0],    MPI_DOUBLE,
                     &halo_var_data[0], &num_values_to_receive[0], &values_to_receive_displacements[0], MPI_DOUBLE,
                     MPI_COMM_WORLD);
  

  // create the points information
  vtkNew<vtkDoubleArray> pointArray;
  pointArray->SetNumberOfComponents(3);
  pointArray->SetNumberOfValues((dataSorter->GetnPoints()+num_halo_nodes)*NCOORDS);
  int id = 0;
  for (int iRank = 0; iRank < size; ++iRank) {
    for (int i = 0; i < num_nodes_to_receive[iRank]; i++){  
      int displ = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank]*0; 
      pointArray->SetValue((dataSorter->GetnPoints() + id)*NCOORDS+0, halo_var_data[displ+i]);
      id++;    
    }
  }
  id = 0;
  for (int iRank = 0; iRank < size; ++iRank) {
    for (int i = 0; i < num_nodes_to_receive[iRank]; i++){  
      int displ = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank]*1; 
      pointArray->SetValue((dataSorter->GetnPoints() + id)*NCOORDS+1, halo_var_data[displ+i]);
      id++;    
    }
  }
  id = 0;
  for (int iRank = 0; iRank < size; ++iRank) {
    for (int i = 0; i < num_nodes_to_receive[iRank]; i++){  
      int displ = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank]*2; 
      if (dataSorter->GetnDim() == 3) pointArray->SetValue((dataSorter->GetnPoints() + id)*NCOORDS+2, halo_var_data[displ+i]);
      else {pointArray->SetValue((dataSorter->GetnPoints() + id)*NCOORDS+2, 0.0);}
      id++;    
    }
  }
  
  for (unsigned long iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++){
    pointArray->SetValue(iPoint*NCOORDS+0, dataSorter->GetData(0,iPoint));
    pointArray->SetValue(iPoint*NCOORDS+1, dataSorter->GetData(1,iPoint));
    if (dataSorter->GetnDim() == 3 ) pointArray->SetValue(iPoint*NCOORDS+2, dataSorter->GetData(2,iPoint));
    else pointArray->SetValue(iPoint*NCOORDS+2, 0);
  }
 
    
  vtkNew<vtkPoints> points;
  points->SetData(pointArray.GetPointer());
  VTKGrid->SetPoints(points.GetPointer());

  // create the cells
  VTKGrid->Allocate(static_cast<vtkIdType>(myElemStorage));
  
#define MAKE_LOCAL(n) dataSorter->FindProcessor(n) == rank \
  ? (int64_t)((unsigned long)n - dataSorter->GetNodeBegin(rank)) \
  : GetHaloNodeNumber(n, dataSorter->GetNodeEnd(rank) - dataSorter->GetNodeBegin(rank), sorted_halo_nodes)

#define ISN_LOCAL(n) (dataSorter->FindProcessor(n) != rank)
  
  for (unsigned long iElem = 0; iElem < nParallel_Line; iElem++) {
    vtkIdType tmp[2] = {(vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(LINE, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(LINE, iElem, 1)-1))};
    VTKGrid->InsertNextCell(VTK_LINE, N_POINTS_LINE, tmp);    
  }
  
  
  for (unsigned long iElem = 0; iElem < nParallel_Tria; iElem++) {    
    vtkIdType tmp[3] = {(vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 1)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 2)-1))};
    VTKGrid->InsertNextCell(VTK_TRIANGLE, N_POINTS_TRIANGLE, tmp);   
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Quad; iElem++) {
    vtkIdType tmp[4] = {(vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 1)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 2)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 3)-1))};
    VTKGrid->InsertNextCell(VTK_QUAD, N_POINTS_QUADRILATERAL, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Tetr; iElem++) {
    vtkIdType tmp[4] = {(vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 1)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 2)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3)-1))};
    VTKGrid->InsertNextCell(VTK_TETRA, N_POINTS_TETRAHEDRON, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Hexa; iElem++) {
    vtkIdType tmp[8] = {(vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 1)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 2)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 3)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 4)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 5)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 6)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 7)-1))};
    VTKGrid->InsertNextCell(VTK_HEXAHEDRON, N_POINTS_HEXAHEDRON, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Pris; iElem++) {
    vtkIdType tmp[6] = {(vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(PRISM, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(PRISM, iElem, 1)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(PRISM, iElem, 2)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(PRISM, iElem, 3)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(PRISM, iElem, 4)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(PRISM, iElem, 5)-1))};
    VTKGrid->InsertNextCell(VTK_WEDGE, N_POINTS_PRISM, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Pyra; iElem++) {
    vtkIdType tmp[5] = {(vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(PYRAMID, iElem, 0)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(PYRAMID, iElem, 1)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(PYRAMID, iElem, 2)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(PYRAMID, iElem, 3)-1)),
                        (vtkIdType)(MAKE_LOCAL(dataSorter->GetElem_Connectivity(PYRAMID, iElem, 4)-1))};
    VTKGrid->InsertNextCell(VTK_PYRAMID, N_POINTS_PYRAMID, tmp);    
  }
}


void CCatalystWriter::UpdateVTKAttributes(vtkCPInputDataDescription* idd, CParallelDataSorter *dataSorter)
{
  LOG_SCOPE_FUNCTION(INFO);
  
  ERROR_CONTEXT("Number of points", dataSorter->GetnPoints());
  
  const vector<string>& fieldNames = dataSorter->GetFieldNames();
  
  bool allocate = false;
  
  if (VTKGrid->GetPointData()->GetNumberOfArrays() == 0){
    
    allocate = true;
    
  }  
  
  unsigned short iVar;
  unsigned long index= 0;
  for(int iRank = 0; iRank < size; ++iRank) {

    for(iVar = 0; iVar < fieldNames.size(); ++iVar)
      for(int iNode = 0; iNode < num_nodes_to_send[iRank]; ++iNode) {
        unsigned long node_offset = nodes_to_send[nodes_to_send_displacements[iRank] + iNode] - dataSorter->GetNodeBegin(rank);
        data_to_send[index++] = SU2_TYPE::GetValue(dataSorter->GetData(iVar,node_offset));
      }
  }
  SU2_MPI::Alltoallv(&data_to_send[0],  &num_values_to_send[0],    &values_to_send_displacements[0],    MPI_DOUBLE,
                     &halo_var_data[0], &num_values_to_receive[0], &values_to_receive_displacements[0], MPI_DOUBLE,
                     MPI_COMM_WORLD);
  
  unsigned long TotalnPoints = dataSorter->GetnPoints()+num_halo_nodes;
  for (unsigned short iField = dataSorter->GetnDim(); iField < fieldNames.size(); iField++){
    string fieldname = fieldNames[iField];
    
    bool output_variable = true, isVector = false;
    size_t found = fieldNames[iField].find("_x");
    if (found!=string::npos) {
      output_variable = true;
      isVector = true;
    }
    found = fieldNames[iField].find("_y");
    if (found!=string::npos) {
      //skip
      output_variable = false;
    }
    found = fieldNames[iField].find("_z");
    if (found!=string::npos) {
      //skip
      output_variable = false;
    }
    
    if (idd->IsFieldNeeded(fieldNames[iField].c_str(), vtkDataObject::POINT) == true && !isVector && output_variable){
      if (allocate){
        
        vtkNew<vtkDoubleArray> array;
        array->SetName(fieldNames[iField].c_str());
        array->SetNumberOfComponents(1);
        array->SetNumberOfValues(static_cast<vtkIdType>(TotalnPoints));
        
        VTKGrid->GetPointData()->AddArray(array.GetPointer());
        
        SU2_INFO << "Adding catalyst output " << fieldNames[iField];
      }
      
      SU2_INFO << "Setting catalyst output " << fieldNames[iField];
      
      vtkDoubleArray* data =
          vtkDoubleArray::SafeDownCast(VTKGrid->GetPointData()->GetArray(fieldNames[iField].c_str()));
      //      data->SetArray(data_sorter->GetData(iField), static_cast<vtkIdType>(data_sorter->GetnPoints()), 1);
      for (vtkIdType i = 0; i < dataSorter->GetnPoints(); i++)
      {
        su2double val = 0.0;
          val = dataSorter->GetData(iField, i);
        data->SetValue(i, val);
      }
      int id = 0;
      for (int iRank = 0; iRank < size; ++iRank) {
        for (int i = 0; i < num_nodes_to_receive[iRank]; i++){  
          int displ = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank]*iField; 
          data->SetValue(dataSorter->GetnPoints() + id, halo_var_data[displ+i]);
          id++;    
        }
      }
    }
    if (isVector && output_variable){
      fieldname.erase(fieldname.end()-2,fieldname.end());
      
      if (idd->IsFieldNeeded(fieldname.c_str(), vtkDataObject::POINT) == true ){
        if (allocate){
          // velocity array
           vtkNew<vtkDoubleArray> array;
           array->SetName(fieldname.c_str());
           array->SetNumberOfComponents(dataSorter->GetnDim());
           array->SetNumberOfTuples(static_cast<vtkIdType>(TotalnPoints));
           VTKGrid->GetPointData()->AddArray(array.GetPointer());
           SU2_INFO << "Adding catalyst output " << fieldname.c_str();
          
        }
        SU2_INFO << "Setting catalyst output " << fieldname;
        
        vtkDoubleArray* data =
            vtkDoubleArray::SafeDownCast(VTKGrid->GetPointData()->GetArray(fieldname.c_str()));
        
        for (vtkIdType i = 0; i < dataSorter->GetnPoints(); i++)
        {            
          data->SetValue(i*dataSorter->GetnDim() + 0, dataSorter->GetData(iField+0, i));
          data->SetValue(i*dataSorter->GetnDim() + 1, dataSorter->GetData(iField+1, i));
          if (dataSorter->GetnDim() == 3){
            data->SetValue(i*dataSorter->GetnDim() + 2, dataSorter->GetData(iField+2, i));            
          } 
        }
        int id = 0;
        for (int iRank = 0; iRank < size; ++iRank) {
          for (int i = 0; i < num_nodes_to_receive[iRank]; i++){  
            int displ = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank]*(iField+0); 
            data->SetValue((dataSorter->GetnPoints() + id)*dataSorter->GetnDim()+0, halo_var_data[displ+i]);
            id++;    
          }
        }
        id = 0;
        for (int iRank = 0; iRank < size; ++iRank) {
          for (int i = 0; i < num_nodes_to_receive[iRank]; i++){  
            int displ = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank]*(iField+1); 
            data->SetValue((dataSorter->GetnPoints() + id)*dataSorter->GetnDim()+1, halo_var_data[displ+i]);
            id++;    
          }
        }
        if (dataSorter->GetnDim() == 3){
          id = 0;
          for (int iRank = 0; iRank < size; ++iRank) {
            for (int i = 0; i < num_nodes_to_receive[iRank]; i++){  
              int displ = values_to_receive_displacements[iRank] + num_nodes_to_receive[iRank]*(iField+2); 
              data->SetValue((dataSorter->GetnPoints() + id)*dataSorter->GetnDim()+2, halo_var_data[displ+i]);
              id++;    
            }
          }
        }
      }
      
    }
  }
}

void CCatalystWriter::BuildVTKDataStructures(vtkCPInputDataDescription* idd, CParallelDataSorter *dataSorter)
{
  LOG_SCOPE_FUNCTION(INFO);
  
  if (VTKGrid == NULL)
  {
    SU2_INFO << "Creating VTK Grid";
    // The grid structure isn't changing so we only build it
    // the first time it's needed. If we needed the memory
    // we could delete it and rebuild as necessary.
    VTKGrid = vtkUnstructuredGrid::New();
    BuildVTKGrid();
  }
  UpdateVTKAttributes(idd, dataSorter);
  
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
