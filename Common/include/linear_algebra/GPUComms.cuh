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

/*!< \brief Namespace that contains variables and helper functions that are
    utilized to launch CUDA Kernels. */
namespace kernelParameters{


  /*!
   * \brief Returns the rounded up value of the decimal quotient to the next integer (in all cases).
   */
  inline constexpr int rounded_up_division(const int divisor, int dividend) { return ((dividend + divisor - 1) / divisor); }

  /*!
   * \brief Returns the rounded down value of the decimal quotient to the previous integer (in all cases).
   */
  inline constexpr int rounded_down_division(const int divisor, int dividend) { return ((dividend - divisor + 1) / divisor); }

  static constexpr short BLOCK_SIZE = 640;
  static constexpr short WARP_SIZE = 32;
  static constexpr short ROWS_PER_BLOCK = rounded_up_division(WARP_SIZE, BLOCK_SIZE);


};

/*!< \brief Structure containing information related to the Jacobian Matrix
    which is utilized by any launched Kernel. */
struct matrixParameters{

  public:
    unsigned long totalRows;        /*!< \brief Contains the total number of rows of the Jacbian Matrix. */
    unsigned short blockRowSize;    /*!< \brief Contains the row dimensions of the blocks of the Jacobian Matrix. */
    unsigned short blockColSize;    /*!< \brief Contains the column dimensions of the blocks of the Jacobian Matrix. */
    unsigned int nChainStart;       /*!< \brief Starting partition of the current chain. */
    unsigned int nChainEnd;         /*!< \brief Ending partition of the current chain. */
    unsigned short blockSize;       /*!< \brief Contains the total number of elements in each block of the Jacbian Matrix. */
    unsigned short activeThreads;   /*!< \brief Cotains the number of active threads per iteration during MVP - depending on the
                                        dimensions of the Jacbian Matrix. */

    matrixParameters(unsigned long nPointDomain, unsigned long nEqn, unsigned long nVar, unsigned long nPartitions){
      totalRows = nPointDomain;
      blockRowSize = nEqn;
      blockColSize = nVar;
      nChainStart = 0;
      nChainEnd = 0;
      blockSize = nVar * nEqn;
      activeThreads = nVar * (kernelParameters::WARP_SIZE/nVar);
    }

    /*!
    * \brief Returns the memory index in the shared memory array used by the Symmetric Iteration Kernels.
    */
    __device__ unsigned short shrdMemIndex(unsigned short localRow, unsigned short threadNo){
      return (localRow * blockSize + threadNo);
    }

    /*!
    * \brief Returns a boolean value to check whether the row is under the total number of rows and if the
    *        thread number is within a user-specified thread limit. This is to avoid illegal memory accesses.
    */
    __device__ bool validAccess(unsigned long row, unsigned short threadNo, unsigned short threadLimit){
      return (row<totalRows && threadNo<threadLimit);
    }

    /*!
    * \brief Returns a boolean value to check whether the row is part of the parallel partition being executed and if the
    *        thread number is within a user-specified thread limit. This is to avoid illegal memory accesses.
    * \param[in] rowInPartition - Represents a boolean that indicates the presence/absence of the row in the partition.
    */
    __device__ bool validParallelAccess(bool rowInPartition, unsigned short threadNo, unsigned short threadLimit){
      return (rowInPartition && threadNo<threadLimit);
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
