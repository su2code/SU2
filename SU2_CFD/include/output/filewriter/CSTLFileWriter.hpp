#pragma once
#include "CFileWriter.hpp"

class CSTLFileWriter : public CFileWriter{
private:
  vector<su2double> Pointlist; // holds doubles. 3doubles for a point and 3 points for triangle. subsequent ordered

public:
  
  CSTLFileWriter(vector<string> fields, unsigned short nDim, CGeometry *geometry, CConfig *config);
  
  ~CSTLFileWriter();
  
  void Write_Data(string filename, CParallelDataSorter* data_sorter);

  void Write_Data2(string filename, CParallelDataSorter* data_sorter);
  
};

