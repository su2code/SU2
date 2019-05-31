#pragma once

#include "CFileWriter.hpp"

class CParaviewBinaryFileWriter : public CFileWriter{
  
public:
  
  CParaviewBinaryFileWriter(vector<string> fields, unsigned short nDim);
  
  ~CParaviewBinaryFileWriter();
  
  void Write_Data(string filename, CParallelDataSorter* data_sorter);
  
private:
  void SwapBytes(char *buffer, size_t nBytes, unsigned long nVar);
  
};

