/*!
\file GPUComms.cuh
* \brief Header file containing universal functions that provide basic and essential utilities for other GPU processes
* \author A. Raj
* \version 8.3.0 "Harrier"
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

namespace KernelParameters{

inline constexpr int round_up_division(const int multiple, int x) { return ((x + multiple - 1) / multiple); }

static constexpr int MVP_BLOCK_SIZE = 1024;
static constexpr int MVP_WARP_SIZE = 32;
}
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
