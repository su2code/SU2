#include "../../../include/output/filewriter/CFileWriter.hpp"


CFileWriter::CFileWriter(vector<string> fields, string fileName,
                         CParallelDataSorter *dataSorter, string file_ext, unsigned short nDim): 
  fieldnames(std::move(fields)),
  nDim(nDim),
  file_ext(file_ext),
  fileName(std::move(fileName)),
  dataSorter(dataSorter){
  
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  this->fileName += file_ext;
  
  file_size = 0.0;
  
}


CFileWriter::~CFileWriter(){
  
}
