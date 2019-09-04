#include "../../../include/output/filewriter/CFileWriter.hpp"


CFileWriter::CFileWriter(vector<string> fields, string file_ext, unsigned short nDim): fieldnames(fields),
                                                                                       nDim(nDim),
                                                                                       file_ext(file_ext){
  
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
}


CFileWriter::~CFileWriter(){
  
}
