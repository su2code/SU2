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
#include <set>
class CParallelDataSorter;

class CCatalystWriter : public CFileWriter {
  
  vtkCPProcessor* Processor = NULL;
  vtkUnstructuredGrid* VTKGrid;
  std::set<unsigned long> halo_nodes;
  vector<unsigned long> sorted_halo_nodes;
  
  vector<passivedouble> data_to_send;
  vector<passivedouble> halo_var_data;
  vector<int> num_values_to_send;
  vector<int> values_to_send_displacements;
  vector<int> values_to_receive_displacements;
  vector<unsigned long> nodes_to_send;
  vector<int> num_values_to_receive;
  vector<int> nodes_to_send_displacements;
  vector<int> nodes_to_receive_displacements;
  vector<int> num_nodes_to_send;
  size_t num_halo_nodes;
  vector<int> num_nodes_to_receive;
public:
  /*!
   * \brief Construct a file writer using field names, file extension and dimension.
   * \param[in] fields - A list of field names
   * \param[in] nDim - Physical dimension
   * \param[in] fileName - The name of the file
   * \param[in] data_sorter - The parallel sorted data to write
   */  
  CCatalystWriter(string fileName, CParallelDataSorter* data_sorter);

  ~CCatalystWriter();  
  
  void BuildVTKGrid();
  
  void UpdateVTKAttributes(vtkCPInputDataDescription* idd, CParallelDataSorter *dataSorter);
  void BuildVTKDataStructures(vtkCPInputDataDescription* idd, CParallelDataSorter *dataSorter);
  
  void Write_Data(unsigned long TimeStep, double time);
  
  int64_t GetHaloNodeNumber(unsigned long global_node_number, unsigned long last_local_node, vector<unsigned long> const &halo_node_list);
};
