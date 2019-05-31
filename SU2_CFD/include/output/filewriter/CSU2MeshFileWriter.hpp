#pragma once
#include "CFileWriter.hpp"

class CSU2MeshFileWriter : public CFileWriter{
  
private:
  unsigned short iZone, nZone;

public:
  
  CSU2MeshFileWriter(string filename, vector<string> fields, unsigned short nDim, unsigned short iZone, unsigned short nZone);
  
  ~CSU2MeshFileWriter();
  
  void Write_Data(CParallelDataSorter* data_sorter);
  
};

