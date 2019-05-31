#pragma once
#include "../../../Common/include/geometry_structure.hpp"
#include "CParallelDataSorter.hpp"

class CFileWriter{
protected:

  string filename;
  vector<string> fieldnames;
  unsigned short nDim;
  
  int rank, size;
  
  string file_ext;
  
public:
  
  CFileWriter(string filename, vector<string> fields, unsigned short nDim);
  
  virtual ~CFileWriter();
  
  virtual void Write_Data(CParallelDataSorter* data_sorter){}
  
};

