#pragma once

#include "CFileWriter.hpp"

class CTecplotFileWriter : public CFileWriter{
  
  unsigned long time_iter;
  su2double timestep;
  
public:
  
  CTecplotFileWriter(string filename, vector<string> fields, unsigned short nDim, unsigned long time_iter, su2double timestep);
  
  ~CTecplotFileWriter();
  
  void Write_Data(CParallelDataSorter* data_sorter);
  
};

