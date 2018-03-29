

#include "matrix_interface_structure.hpp"



template<class CalcType, class BaseType>
void TCMatrixInterface<CalcType, BaseType>::Initialize(unsigned short blocksize, bool edgeconnect, CGeometry *geometry, CConfig *config){
  
  BlockSize = blocksize;
  unsigned short MaxSize = BlockSize*BlockSize;
  
  Block_CalcType = new CalcType*[MaxSize];
  Block_BaseType = new BaseType*[MaxSize];

  for (unsigned short i = 0; i < MaxSize; i++){
    Block_CalcType[i] = new CalcType[MaxSize];
    Block_BaseType[i] = new BaseType[MaxSize];
  } 
  
  BlockLin_CalcType = new CalcType[MaxSize];
  BlockLin_BaseType = new BaseType[MaxSize];
  
}

template<class CalcType, class BaseType>
TCMatrixInterface<CalcType, BaseType>::~TCMatrixInterface() {
  
  for (unsigned short i=0; i < BlockSize; i++){
    if (Block_CalcType[i] != NULL) delete [] Block_CalcType[i];
    if (Block_BaseType[i] != NULL) delete [] Block_BaseType[i];
  } 
  
  if (BlockLin_CalcType != NULL) delete [] BlockLin_CalcType;
  if (BlockLin_BaseType != NULL) delete [] BlockLin_BaseType;
  
  if (Block_BaseType != NULL) delete [] Block_BaseType;
  if (Block_CalcType != NULL) delete [] Block_CalcType;
  
}