/*!
 * \file CSysMatrixGPU.cu
 * \brief Implementations of Kernels and Functions for Matrix Operations on the GPU
 * \author A. Raj, Jesse Li, D. Di Giusto
 * \version 8.5.0 "Harrier"
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

#include "../../include/linear_algebra/CSysMatrix.hpp"
#include "../../include/linear_algebra/CSysMatrix.inl"
#include "../../include/linear_algebra/GPUComms.cuh"

/*!
 * \brief Matrix-vector product kernel.
 */
template<typename matrixType, typename vectorType>
__global__ void GPUMatrixVectorProductKernel(matrixType *invM, vectorType* vec, vectorType* prod, unsigned long nPointDomain, unsigned long nVar)
{
  
  const unsigned long iPoint = blockIdx.x * blockDim.x + threadIdx.x;
  if (iPoint >= nPointDomain) return;

  const auto block = &invM[iPoint * nVar * nVar];
  const auto rhs = &vec[iPoint * nVar];
  auto out = &prod[iPoint * nVar];

  for (auto iVar = 0; iVar < nVar; ++iVar) {
    vectorType sum = vectorType(0);
    for (auto jVar = 0; jVar < nVar; ++jVar) {
      sum += block[iVar * nVar + jVar] * rhs[jVar];
    }
    out[iVar] = sum;
  }
}

/*!
 * \brief Block-LDU SpMV kernel: y[iRow] = (L + D + U) * x per block-row.
 *        One CUDA block per block-row; threadIdx.x indexes output variable (0..nVar-1).
 */
template <class ScalarType>
__global__ void BlockLDU_SpMV_kernel(unsigned long nRows, unsigned long nVar,
                                     const unsigned long* __restrict__ row_ptr_l,
                                     const unsigned long* __restrict__ col_ind_l,
                                     const ScalarType* __restrict__ mat_l,
                                     const ScalarType* __restrict__ mat_d,
                                     const unsigned long* __restrict__ row_ptr_u,
                                     const unsigned long* __restrict__ col_ind_u,
                                     const ScalarType* __restrict__ mat_u,
                                     const ScalarType* __restrict__ x, ScalarType* __restrict__ y) {
  const unsigned long iRow = blockIdx.x;
  const unsigned long iVar = threadIdx.x;
  if (iRow >= nRows || iVar >= nVar) return;

  ScalarType sum = 0;
  /* Lower */
  for (auto k = row_ptr_l[iRow]; k < row_ptr_l[iRow + 1]; ++k) {
    const auto col = col_ind_l[k];
    const ScalarType* blk = mat_l + k * nVar * nVar + iVar * nVar;
    for (unsigned long jVar = 0; jVar < nVar; ++jVar) sum += blk[jVar] * x[col * nVar + jVar];
  }
  /* Diagonal */
  {
    const ScalarType* blk = mat_d + iRow * nVar * nVar + iVar * nVar;
    for (unsigned long jVar = 0; jVar < nVar; ++jVar) sum += blk[jVar] * x[iRow * nVar + jVar];
  }
  /* Upper */
  for (auto k = row_ptr_u[iRow]; k < row_ptr_u[iRow + 1]; ++k) {
    const auto col = col_ind_u[k];
    const ScalarType* blk = mat_u + k * nVar * nVar + iVar * nVar;
    for (unsigned long jVar = 0; jVar < nVar; ++jVar) sum += blk[jVar] * x[col * nVar + jVar];
  }
  y[iRow * nVar + iVar] = sum;
}

template <class ScalarType>
void CSysMatrix<ScalarType>::HtDTransfer(bool trigger) const {
  if (!trigger) return;
  gpuErrChk(cudaMemcpy(gpu.d, mat.d, sizeof(ScalarType) * nPoint * nVar * nEqn, cudaMemcpyHostToDevice));
  gpuErrChk(cudaMemcpy(gpu.l, mat.l, sizeof(ScalarType) * mat.nnz_l * nVar * nEqn, cudaMemcpyHostToDevice));
  gpuErrChk(cudaMemcpy(gpu.u, mat.u, sizeof(ScalarType) * mat.nnz_u * nVar * nEqn, cudaMemcpyHostToDevice));
}

template <class ScalarType>
void CSysMatrix<ScalarType>::GPUMatrixVectorProduct(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                                    CGeometry* geometry, const CConfig* config) const {
  if (nVar != nEqn) {
    SU2_MPI::Error("CUDA CSysMatrix block-LDU SpMV requires square blocks.", CURRENT_FUNCTION);
  }

  ScalarType* d_vec = vec.GetDevicePointer();
  ScalarType* d_prod = prod.GetDevicePointer();
  vec.HtDTransfer();

  dim3 blockDim(static_cast<unsigned>(nVar), 1, 1);
  dim3 gridDim(static_cast<unsigned>(nPointDomain), 1, 1);
  BlockLDU_SpMV_kernel<ScalarType><<<gridDim, blockDim>>>(
      nPointDomain, nVar, gpu.row_ptr_l, gpu.col_ind_l, gpu.l, gpu.d,
      gpu.row_ptr_u, gpu.col_ind_u, gpu.u, d_vec, d_prod);
  gpuErrChk(cudaGetLastError());

  prod.DtHTransfer();
}

template <class ScalarType>
void CSysMatrix<ScalarType>::GPUComputeJacobiPreconditioner(const CSysVector<ScalarType>& vec,
                                                            CSysVector<ScalarType>& prod, CGeometry* geometry,
                                                            const CConfig* config) const {
  SU2_ZONE_SCOPED
  /*--- Apply Jacobi preconditioner, y = D^{-1} * x, the inverse of the diagonal is already known and synced to device ---*/

  ScalarType* d_vec = vec.GetDevicePointer();
  ScalarType* d_prod = prod.GetDevicePointer();

  vec.HtDTransfer(); // this is now the entry point of the cuda section so we always want to copy
  prod.GPUSetVal(0.0);

  dim3 blockDim(KernelParameters::MVP_BLOCK_SIZE,1,1);
  int gridx = KernelParameters::round_up_division(KernelParameters::MVP_BLOCK_SIZE, nPointDomain);
  dim3 gridDim(gridx, 1, 1);

  GPUMatrixVectorProductKernel<<<gridDim, blockDim>>>(d_invM, d_vec, d_prod, nPointDomain, nVar);
  gpuErrChk( cudaPeekAtLastError() );

  prod.DtHTransfer();//forcely copy back prod to host for MPI synch

  /*--- MPI Parallelization ---*/
  CSysMatrixComms::Initiate(prod, geometry, config);
  CSysMatrixComms::Complete(prod, geometry, config);
  prod.HtDTransfer();//forcely copy back prod to device to continue calculations

}

template <class ScalarType>
void CSysMatrix<ScalarType>::GPUBuildJacobiPreconditioner() {
  SU2_ZONE_SCOPED
  /*--- Build Jacobi preconditioner (M = D), compute and store the inverses of the diagonal blocks. ---*/
  SU2_OMP_FOR_DYN(omp_heavy_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
    InverseDiagonalBlock(iPoint, &(invM[iPoint * nVar * nVar]));
  END_SU2_OMP_FOR

  //copy to device or prefetch
  if (invM_is_managed) {
    gpu_um_prefetch(d_invM, nPointDomain * nVar * nVar * sizeof(ScalarType), GPUMemoryAllocation::GetCurrentDevice());
  } else {
    gpuErrChk(cudaMemcpy(d_invM, invM, nPointDomain * nVar * nVar * sizeof(ScalarType), cudaMemcpyHostToDevice));
  }
}

template class CSysMatrix<su2mixedfloat>; //This is a temporary fix for invalid instantiations due to separating the member function from the header file the class is defined in. Will try to rectify it in coming commits.