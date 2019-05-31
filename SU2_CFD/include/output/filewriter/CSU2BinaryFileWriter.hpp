#pragma once
#include "CFileWriter.hpp"

class CSU2BinaryFileWriter : public CFileWriter{
  

public:
  
  CSU2BinaryFileWriter(vector<string> fields, unsigned short nDim);
  
  ~CSU2BinaryFileWriter();
  
  void Write_Data(string filename, CParallelDataSorter* data_sorter);
  
};
