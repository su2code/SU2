#pragma once
#include "CFileWriter.hpp"
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

class CParallelDataSorter;

class CCatalystWriter : public CFileWriter {
  
  vtkCPProcessor* Processor = NULL;
  vtkUnstructuredGrid* VTKGrid;
  
public:
  CCatalystWriter(vector<string> fields, unsigned short nDim);

  ~CCatalystWriter();  
  
  void BuildVTKGrid(CParallelDataSorter *data_sorter);
  
  void UpdateVTKAttributes(vtkCPInputDataDescription* idd, CParallelDataSorter *data_sorter);
  void BuildVTKDataStructures(vtkCPInputDataDescription* idd, CParallelDataSorter *data_sorter);
  
  void Write_Data(unsigned long TimeStep, double time, CParallelDataSorter *data_sorter);
};
