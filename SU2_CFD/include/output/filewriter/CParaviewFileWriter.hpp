#pragma once

#include "CFileWriter.hpp"

class CParaviewFileWriter : public CFileWriter{
  
public:
  
  CParaviewFileWriter(string filename, vector<string> fields, unsigned short nDim);
  
  ~CParaviewFileWriter();
  
  void Write_Data(CParallelDataSorter* data_sorter);
  
};

