#include "../../../include/output/filewriter/CCatalystWriter.hpp"
#include "../../../include/output/filewriter/CParallelDataSorter.hpp"
#include "../../../Common/include/toolboxes/SU2_LOG.hpp"


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
  
  // create the points information
  vtkNew<vtkDoubleArray> pointArray;
  pointArray->SetNumberOfComponents(3);
  pointArray->SetNumberOfValues(data_sorter->GetnPoints()*NCOORDS);
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

  
  for (unsigned long iElem = 0; iElem < nParallel_Line; iElem++) {
    vtkIdType tmp[2] = {(vtkIdType)data_sorter->GetElem_Connectivity(LINE, iElem, 0)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(LINE, iElem, 1)-1};
    VTKGrid->InsertNextCell(VTK_LINE, N_POINTS_LINE, tmp);    
  }
  
  
  for (unsigned long iElem = 0; iElem < nParallel_Tria; iElem++) {
    vtkIdType tmp[3] = {(vtkIdType)data_sorter->GetElem_Connectivity(TRIANGLE, iElem, 0)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(TRIANGLE, iElem, 1)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(TRIANGLE, iElem, 2)-1};
    VTKGrid->InsertNextCell(VTK_TRIANGLE, N_POINTS_TRIANGLE, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Quad; iElem++) {
    vtkIdType tmp[4] = {(vtkIdType)data_sorter->GetElem_Connectivity(QUADRILATERAL, iElem, 0)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(QUADRILATERAL, iElem, 1)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(QUADRILATERAL, iElem, 2)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(QUADRILATERAL, iElem, 3)-1};
    VTKGrid->InsertNextCell(VTK_QUAD, N_POINTS_QUADRILATERAL, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Tetr; iElem++) {
    vtkIdType tmp[4] = {(vtkIdType)data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 0)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 1)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 2)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3)-1};
    VTKGrid->InsertNextCell(VTK_TETRA, N_POINTS_TETRAHEDRON, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Hexa; iElem++) {
    vtkIdType tmp[8] = {(vtkIdType)data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 0)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 1)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 2)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 3)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 4)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 5)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 6)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(HEXAHEDRON, iElem, 7)-1};
    VTKGrid->InsertNextCell(VTK_HEXAHEDRON, N_POINTS_HEXAHEDRON, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Pris; iElem++) {
    vtkIdType tmp[6] = {(vtkIdType)data_sorter->GetElem_Connectivity(PRISM, iElem, 0)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(PRISM, iElem, 1)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(PRISM, iElem, 2)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(PRISM, iElem, 3)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(PRISM, iElem, 4)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(PRISM, iElem, 5)-1};
    VTKGrid->InsertNextCell(VTK_WEDGE, N_POINTS_PRISM, tmp);    
  }
  
  for (unsigned long iElem = 0; iElem < nParallel_Pyra; iElem++) {
    vtkIdType tmp[6] = {(vtkIdType)data_sorter->GetElem_Connectivity(PYRAMID, iElem, 0)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(PYRAMID, iElem, 1)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(PYRAMID, iElem, 2)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(PYRAMID, iElem, 3)-1,
                        (vtkIdType)data_sorter->GetElem_Connectivity(PYRAMID, iElem, 4)-1};
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
  
  for (unsigned short iField = nDim; iField < fieldnames.size(); iField++){
    
    if (idd->IsFieldNeeded(fieldnames[iField].c_str(), vtkDataObject::POINT) == true){
      if (allocate){
        
        vtkNew<vtkDoubleArray> array;
        array->SetName(fieldnames[iField].c_str());
        array->SetNumberOfComponents(1);
        VTKGrid->GetPointData()->AddArray(array.GetPointer());
        
        SU2_INFO << "Adding catalyst output " << fieldnames[iField];
      }
      
      SU2_INFO << "Setting catalyst output " << fieldnames[iField];
      
      vtkDoubleArray* data =
          vtkDoubleArray::SafeDownCast(VTKGrid->GetPointData()->GetArray(fieldnames[iField].c_str()));
      data->SetArray(data_sorter->GetData(iField), static_cast<vtkIdType>(data_sorter->GetnPoints()), 1);
      
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





