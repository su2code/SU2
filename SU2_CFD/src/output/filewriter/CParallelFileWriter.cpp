#include "../../../include/output/filewriter/CFileWriter.hpp"


CFileWriter::CFileWriter(string filename, vector<string> fields, unsigned short nDim){
  
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  this->nDim = nDim;
  
  this->filename = filename;
  
  this->fieldnames = fields;
  
}


CFileWriter::~CFileWriter(){
  
}
