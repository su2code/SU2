#pragma once
#include "CFileWriter.hpp"

class CCSVFileWriter : public CFileWriter{
  

public:
  
  CCSVFileWriter(string filename, vector<string> fields, unsigned short nDim);
  
  ~CCSVFileWriter();
  
  void Write_Data(string filename, CParallelDataSorter* data_sorter);
  
};

