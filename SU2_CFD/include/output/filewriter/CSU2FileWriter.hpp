#pragma once
#include "CFileWriter.hpp"

class CSU2FileWriter : public CFileWriter{
  

public:
  
  CSU2FileWriter(vector<string> fields, unsigned short nDim);
  
  ~CSU2FileWriter();
  
  void Write_Data(string filename, CParallelDataSorter* data_sorter);
  
};

