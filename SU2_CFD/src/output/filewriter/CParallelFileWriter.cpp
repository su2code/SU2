#include "../../../include/output/filewriter/CFileWriter.hpp"


CFileWriter::CFileWriter(vector<string> fields, unsigned short nDim){
  
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  this->nDim = nDim;
    
  this->fieldnames = fields;
  
}


CFileWriter::~CFileWriter(){
  
}
