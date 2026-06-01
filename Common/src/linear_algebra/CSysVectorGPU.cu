/*!
 * \file CSysVectorGPU.cu
 * \brief Implementations of Kernels and Functions for Vector Operations on the GPU
 * \author A. Raj
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

#include "../../include/linear_algebra/CSysVector.hpp"
#include "../../include/linear_algebra/GPUComms.cuh"
template<class ScalarType>
void CSysVector<ScalarType>::HtDTransfer(bool trigger) const
{
   if(trigger) gpuErrChk(cudaMemcpy((void*)(d_vec_val), (void*)&vec_val[0], (sizeof(ScalarType)*nElm), cudaMemcpyHostToDevice));
}

template<class ScalarType>
void CSysVector<ScalarType>::DtHTransfer(bool trigger) const
{
   if(trigger) gpuErrChk(cudaMemcpy((void*)(&vec_val[0]), (void*)d_vec_val, (sizeof(ScalarType)*nElm), cudaMemcpyDeviceToHost));
}

template<class ScalarType>
void CSysVector<ScalarType>::GPUSetVal(ScalarType val, bool trigger) const
{
   if(trigger) gpuErrChk(cudaMemset((void*)(d_vec_val), val, (sizeof(ScalarType)*nElm)));
}

template <class ScalarType>
void CSysVector<ScalarType>::GPUCopy(const CSysVector& src) const {
  (void)src;

  /*--- Implementation area for GPU-to-GPU vector copy used by device-resident Krylov solvers. ---*/
  SU2_MPI::Error("CSysVector::GPUCopy skeleton reached without an implementation.", CURRENT_FUNCTION);
}

template <class ScalarType>
void CSysVector<ScalarType>::GPUScale(ScalarType alpha) const {
  (void)alpha;

  /*--- Implementation area for GPU vector scaling used by device-resident Krylov solvers. ---*/
  SU2_MPI::Error("CSysVector::GPUScale skeleton reached without an implementation.", CURRENT_FUNCTION);
}

template <class ScalarType>
void CSysVector<ScalarType>::GPUAxpy(ScalarType alpha, const CSysVector& x) const {
  (void)alpha;
  (void)x;

  /*--- Implementation area for GPU AXPY used by device-resident Krylov solvers. ---*/
  SU2_MPI::Error("CSysVector::GPUAxpy skeleton reached without an implementation.", CURRENT_FUNCTION);
}

template <class ScalarType>
ScalarType CSysVector<ScalarType>::GPUDot(const CSysVector& other) const {
  (void)other;

  /*--- Implementation area for GPU dot products used by device-resident Krylov solvers. ---*/
  SU2_MPI::Error("CSysVector::GPUDot skeleton reached without an implementation.", CURRENT_FUNCTION);
  return ScalarType(0);
}

template <class ScalarType>
ScalarType CSysVector<ScalarType>::GPUNorm() const {
  /*--- Implementation area for GPU norms used by device-resident Krylov solvers. ---*/
  SU2_MPI::Error("CSysVector::GPUNorm skeleton reached without an implementation.", CURRENT_FUNCTION);
  return ScalarType(0);
}

template void CSysVector<su2mixedfloat>::HtDTransfer(bool trigger) const;
template void CSysVector<su2mixedfloat>::DtHTransfer(bool trigger) const;
template void CSysVector<su2mixedfloat>::GPUSetVal(su2mixedfloat val, bool trigger) const;
template void CSysVector<su2mixedfloat>::GPUCopy(const CSysVector<su2mixedfloat>& src) const;
template void CSysVector<su2mixedfloat>::GPUScale(su2mixedfloat alpha) const;
template void CSysVector<su2mixedfloat>::GPUAxpy(su2mixedfloat alpha, const CSysVector<su2mixedfloat>& x) const;
template su2mixedfloat CSysVector<su2mixedfloat>::GPUDot(const CSysVector<su2mixedfloat>& other) const;
template su2mixedfloat CSysVector<su2mixedfloat>::GPUNorm() const;

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
template void CSysVector<passivedouble>::HtDTransfer(bool trigger) const;
template void CSysVector<passivedouble>::DtHTransfer(bool trigger) const;
template void CSysVector<passivedouble>::GPUSetVal(passivedouble val, bool trigger) const;
template void CSysVector<passivedouble>::GPUCopy(const CSysVector<passivedouble>& src) const;
template void CSysVector<passivedouble>::GPUScale(passivedouble alpha) const;
template void CSysVector<passivedouble>::GPUAxpy(passivedouble alpha, const CSysVector<passivedouble>& x) const;
template passivedouble CSysVector<passivedouble>::GPUDot(const CSysVector<passivedouble>& other) const;
template passivedouble CSysVector<passivedouble>::GPUNorm() const;
#endif

#ifdef CODI_REVERSE_TYPE
template void CSysVector<su2double>::HtDTransfer(bool trigger) const;
template void CSysVector<su2double>::DtHTransfer(bool trigger) const;
template void CSysVector<su2double>::GPUSetVal(su2double val, bool trigger) const;
template void CSysVector<su2double>::GPUCopy(const CSysVector<su2double>& src) const;
template void CSysVector<su2double>::GPUScale(su2double alpha) const;
template void CSysVector<su2double>::GPUAxpy(su2double alpha, const CSysVector<su2double>& x) const;
template su2double CSysVector<su2double>::GPUDot(const CSysVector<su2double>& other) const;
template su2double CSysVector<su2double>::GPUNorm() const;
#endif
