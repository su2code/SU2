/*!
 * \file GPU_lin_alg.cuh
 * \brief Declaration of the GPU Matrix Vector Product CUDA Kernel.
 *        The implemtation is in <i>GPU_lin_alg.cu</i>.
 * \author A. Raj
 * \version 8.0.1 "Harrier"
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
#include"../../include/linear_algebra/CSysMatrix.hpp"
#include"iostream"

/*!
   * \brief assert style function that reads return codes after intercepting CUDA API calls.
   *        It returns the result code and its location if the call is unsuccessful.
   * \param[in] code - result code of CUDA function
   * \param[in] file - name of file holding the function
   * \param[in] line - line containing the function
   */
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

/*!
   * \brief CUDA Kernel that performs the Matrix Vector Product: All threads execute product += matrix*vector
   * \param[in] matrix - matrix to be multiplied
   * \param[in] vec - vector to multiply the matrix with
   * \param[in] prod - storing the output of the operation
   * \param[in] d_row_ptr - a device array of pointers pointing to the first non-zero element in each row of the block matrix
   * \param[in] d_col_ind - a device array holding the column index of each element of the block matrix
   * \param[in] nPointDomain - number of real points of the mesh
   * \param[in] nVar - number of variables of the problem
   * \param[in] nEqn - number of equations of the problem
   */

template<typename matrixType, typename vectorType>
__global__ void GPUMatrixVectorProductAdd(matrixType* matrix, vectorType* vec, vectorType* prod, unsigned long* d_row_ptr, unsigned long* d_col_ind, unsigned long nPointDomain, unsigned long nVar, unsigned long nEqn);
