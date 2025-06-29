/*!
\file GPUComms.cuh
* \brief Header file containing universal functions that provide basic and essential utilities for other GPU processes
* \author A. Raj
* \version 8.2.0 "Harrier"
*
* SU2 Project Website: https://su2code.github.io
*
* The SU2 Project is maintained by the SU2 Foundation
* (http://su2foundation.org)
*
* Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
*
* SU2 is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* SU2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with SU2. If not, see <http://www.gnu.org/licenses/>.
*/

#include<cuda_runtime.h>
#include<iostream>

namespace kernelParameters{

  /*Returns the rounded up value of the decimal quotient to the next integer (in all cases)*/
  inline constexpr int rounded_up_division(const int divisor, int dividend) { return ((dividend + divisor - 1) / divisor); }

  /*Returns the rounded down value of the decimal quotient to the previous integer (in all cases)*/
  inline constexpr int rounded_down_division(const int divisor, int dividend) { return ((dividend - divisor + 1) / divisor); }

  static constexpr short MVP_BLOCK_SIZE = 1024;
  static constexpr short MVP_WARP_SIZE = 32;  

};

struct matrixParameters{

  public:
    unsigned long totalRows;
    unsigned short blockRowSize;
    unsigned short blockColSize;
    unsigned long nPartition; 
    unsigned short blockSize;
    unsigned short rowsPerBlock;
    unsigned short activeThreads;

    matrixParameters(unsigned long nPointDomain, unsigned long nEqn, unsigned long nVar, unsigned long nPartitions){
      totalRows = nPointDomain;
      blockRowSize = nEqn;
      blockColSize = nVar;
      nPartition = nPartitions;
      blockSize = nVar * nEqn;
      rowsPerBlock = kernelParameters::rounded_up_division(kernelParameters::MVP_WARP_SIZE, kernelParameters::MVP_BLOCK_SIZE);
      activeThreads = nVar * (kernelParameters::MVP_WARP_SIZE/nVar);
    }

    __device__ bool validAccess(unsigned long row, unsigned short threadNo, unsigned short threadLimit){
      return (row<totalRows && threadNo<threadLimit);
    }

    __device__ bool validParallelAccess(bool rowInPartition, unsigned short threadNo, unsigned short threadLimit){
      return (rowInPartition && threadNo<activeThreads);
    }

};

/*!
* \brief assert style function that reads return codes after intercepting CUDA API calls.
*        It returns the result code and its location if the call is unsuccessful.
* \param[in] code - result code of CUDA function
* \param[in] file - name of file holding the function
* \param[in] line - line containing the function
*/

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true){
  if (code != cudaSuccess){
     fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
     if (abort) exit(code);
  }
}

#define gpuErrChk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
