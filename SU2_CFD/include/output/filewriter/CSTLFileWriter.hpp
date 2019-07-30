#pragma once
#include "CFileWriter.hpp"

class CSTLFileWriter : public CFileWriter{
private:
  vector<su2double> Pointlist; // holds doubles. 3doubles for a point and 3 points for triangle. subsequent ordered

public:
  
  CSTLFileWriter(vector<string> fields, unsigned short nDim);
  
  ~CSTLFileWriter();
  
  void Write_Data(string filename, CParallelDataSorter* data_sorter);
  
};

