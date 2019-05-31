#pragma once
#include "CFileWriter.hpp"

class CCSVFileWriter : public CFileWriter{
  

public:
  
  CCSVFileWriter(vector<string> fields, unsigned short nDim);
  
  ~CCSVFileWriter();
  
  void Write_Data(string filename, CParallelDataSorter* data_sorter);
  
};

