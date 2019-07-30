#pragma once
#include "CFileWriter.hpp"

class CSTLFileWriter : public CFileWriter{

public:
  
  CSTLFileWriter(vector<string> fields, unsigned short nDim);
  
  ~CSTLFileWriter();
  
  void Write_Data(string filename, CParallelDataSorter* data_sorter);
  
};

