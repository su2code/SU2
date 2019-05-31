#pragma once
#include "CFileWriter.hpp"

class CSU2FileWriter : public CFileWriter{
  

public:
  
  CSU2FileWriter(string filename, vector<string> fields, unsigned short nDim);
  
  ~CSU2FileWriter();
  
  void Write_Data(CParallelDataSorter* data_sorter);
  
};

