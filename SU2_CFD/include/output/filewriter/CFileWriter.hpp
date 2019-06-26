#pragma once
#include "../../../Common/include/geometry_structure.hpp"
#include "CParallelDataSorter.hpp"

class CFileWriter{
protected:

  vector<string> fieldnames;
  unsigned short nDim;
  
  int rank, size;
  
  string file_ext;
  
public:
  
  CFileWriter(vector<string> fields, unsigned short nDim);
  
  virtual ~CFileWriter();
  
  virtual void Write_Data(string filename, CParallelDataSorter* data_sorter){}
  
};

