#include "../../../include/output/filewriter/CFileWriter.hpp"


CFileWriter::CFileWriter(vector<string> fields, string fileName,
                         CParallelDataSorter *dataSorter, string file_ext, unsigned short nDim): 
  fieldnames(fields),
  nDim(nDim),
  file_ext(file_ext),
  fileName(fileName),
  dataSorter(dataSorter){
  
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  this->fileName += file_ext;
  
}


CFileWriter::~CFileWriter(){
  
}
