#pragma once
#include "CFileWriter.hpp"

class CSU2BinaryFileWriter : public CFileWriter{
  

public:
  
  CSU2BinaryFileWriter(string filename, vector<string> fields, unsigned short nDim);
  
  ~CSU2BinaryFileWriter();
  
  void Write_Data(CParallelDataSorter* data_sorter);
  
};
