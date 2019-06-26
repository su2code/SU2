#pragma once

#include "CFileWriter.hpp"

class CParaviewFileWriter : public CFileWriter{
  
public:
  
  CParaviewFileWriter(vector<string> fields, unsigned short nDim);
  
  ~CParaviewFileWriter();
  
  void Write_Data(string filename, CParallelDataSorter* data_sorter);
  
};

