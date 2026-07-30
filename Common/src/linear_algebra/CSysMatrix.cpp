/*!
 * \file CSysMatrix.cpp
 * \brief Implementation of the sparse matrix class.
 * \author F. Palacios, A. Bueno, T. Economon, P. Gomes
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/linear_algebra/CSysMatrix.inl"

#include "../../include/geometry/CGeometry.hpp"
#include "../../include/toolboxes/allocation_toolbox.hpp"

#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

namespace {
/*--- Helper function to regularize small pivots ---*/
template <class ScalarType>
FORCEINLINE void RegularizePivot(ScalarType& pivot, unsigned long row, unsigned long col, const char* context) {
  const float eps = 1e-12;
  if (std::abs(pivot) < eps) {
    pivot = std::copysign(eps, SU2_TYPE::GetValue(pivot));
#ifndef NDEBUG
    std::cout << context << ": Regularized small pivot A(" << row << "," << col << ") to " << pivot << std::endl;
#endif
  }
}
}  // namespace

template <class ScalarType>
CSysMatrix<ScalarType>::CSysMatrix() : rank(SU2_MPI::GetRank()), size(SU2_MPI::GetSize()) {
  SU2_ZONE_SCOPED

  nPoint = nPointDomain = nVar = nEqn = 0;
  mat.nnz_l = mat.nnz_u = 0;
  gpu.nnz_l = gpu.nnz_u = 0;
  ilu.nnz_l = ilu.nnz_u = 0;
  ilu_fill_in = 0;

  omp_partitions = nullptr;

  mat.row_ptr_l = nullptr;
  mat.col_ind_l = nullptr;
  mat.row_ptr_u = nullptr;
  mat.col_ind_u = nullptr;
  l_to_u_transp = nullptr;
  u_to_l_transp = nullptr;
  edge_ptr_l = nullptr;

  mat.d = nullptr;
  mat.l = nullptr;
  mat.u = nullptr;

  gpu.d = nullptr;
  gpu.l = nullptr;
  gpu.u = nullptr;
  gpu.row_ptr_l = nullptr;
  gpu.col_ind_l = nullptr;
  gpu.row_ptr_u = nullptr;
  gpu.col_ind_u = nullptr;

  ilu.l = nullptr;
  ilu.d = nullptr;
  ilu.u = nullptr;

  q_scale_l = nullptr;
  q_blocks_l = nullptr;
  q_scale_u = nullptr;
  q_blocks_u = nullptr;
  q_scale_d = nullptr;
  q_blocks_d = nullptr;

  invM = nullptr;
  d_invM = nullptr;

#ifdef USE_MKL
  MatrixMatrixProductJitter = nullptr;
  MatrixVectorProductJitterBetaOne = nullptr;
  MatrixVectorProductJitterBetaZero = nullptr;
  MatrixVectorProductJitterAlphaMinusOne = nullptr;
#endif
}

template <class ScalarType>
CSysMatrix<ScalarType>::~CSysMatrix() {
  SU2_ZONE_SCOPED

  delete[] omp_partitions;
  MemoryAllocation::aligned_free(ilu.l);
  MemoryAllocation::aligned_free(ilu.d);
  MemoryAllocation::aligned_free(ilu.u);
  MemoryAllocation::aligned_free(mat.d);
  MemoryAllocation::aligned_free(mat.l);
  MemoryAllocation::aligned_free(mat.u);
  MemoryAllocation::aligned_free(invM);
  MemoryAllocation::aligned_free(q_scale_l);
  MemoryAllocation::aligned_free(q_blocks_l);
  MemoryAllocation::aligned_free(q_scale_u);
  MemoryAllocation::aligned_free(q_blocks_u);
  MemoryAllocation::aligned_free(q_scale_d);
  MemoryAllocation::aligned_free(q_blocks_d);

  if (useCuda) {
    GPUMemoryAllocation::gpu_free(gpu.d);
    GPUMemoryAllocation::gpu_free(gpu.l);
    GPUMemoryAllocation::gpu_free(gpu.u);
    GPUMemoryAllocation::gpu_free(gpu.row_ptr_l);
    GPUMemoryAllocation::gpu_free(gpu.col_ind_l);
    GPUMemoryAllocation::gpu_free(gpu.row_ptr_u);
    GPUMemoryAllocation::gpu_free(gpu.col_ind_u);
    GPUMemoryAllocation::gpu_free(d_invM);
  }

#ifdef USE_MKL
  mkl_jit_destroy(MatrixMatrixProductJitter);
  mkl_jit_destroy(MatrixVectorProductJitterBetaZero);
  mkl_jit_destroy(MatrixVectorProductJitterBetaOne);
  mkl_jit_destroy(MatrixVectorProductJitterAlphaMinusOne);
#endif
}

template <class ScalarType>
void CSysMatrix<ScalarType>::Initialize(unsigned long npoint, unsigned long npointdomain, unsigned short nvar,
                                        unsigned short neqn, bool EdgeConnect, CGeometry* geometry,
                                        const CConfig* config, bool needTranspPtr, bool grad_mode, bool allow_quant) {
  SU2_ZONE_SCOPED
  assert(omp_get_thread_num() == 0 && "Only the master thread is allowed to initialize the matrix.");

  if (npoint == 0) return;

  if (mat.d != nullptr) {
    SU2_MPI::Error("CSysMatrix can only be initialized once.", CURRENT_FUNCTION);
  }

  if (nvar > MAXNVAR) {
    SU2_MPI::Error("nVar larger than expected, increase MAXNVAR.", CURRENT_FUNCTION);
  }

  /*--- Application of this matrix, FVM or FEM. ---*/
  const auto type = EdgeConnect ? ConnectivityType::FiniteVolume : ConnectivityType::FiniteElement;

  /*--- Type of preconditioner the matrix will be asked to build. ---*/
  auto prec = config->GetKind_Linear_Solver_Prec();

  if ((!EdgeConnect && !config->GetStructuralProblem()) || (config->GetKind_SU2() == SU2_COMPONENT::SU2_DEF) ||
      (config->GetKind_SU2() == SU2_COMPONENT::SU2_DOT)) {
    /*--- FEM-type connectivity in non-structural context implies mesh deformation. ---*/
    prec = config->GetKind_Deform_Linear_Solver_Prec();
  } else if (config->GetDiscrete_Adjoint() && (prec != ILU)) {
    /*--- Else "upgrade" primal solver settings. ---*/
    prec = config->GetKind_DiscAdj_Linear_Prec();
  }

  /*--- No else if, but separate if case! ---*/
  if (config->GetSmoothGradient() && grad_mode) {
    prec = config->GetKind_Grad_Linear_Solver_Prec();
  }

  useCuda = config->GetCUDA();

  const bool ilu_needed = (prec == ILU);
  const bool diag_needed = (prec == JACOBI) || (prec == LINELET);
#ifndef CODI_REVERSE_TYPE
  const bool q_lus_needed = allow_quant && !useCuda && (prec == Q_LU_SGS);
#else
  /*--- No quantization in adjoint mode for now because TransposeInPlace would get complicated. ---*/
  const bool q_lus_needed = false;
#endif

  /*--- Basic dimensions. ---*/
  nVar = nvar;
  nEqn = neqn;
  nPoint = npoint;
  nPointDomain = npointdomain;

  /*--- Allocate host data. ---*/
  auto allocAndInit = [](ScalarType*& ptr, unsigned long num) {
    ptr = MemoryAllocation::aligned_alloc<ScalarType, true>(64, num * sizeof(ScalarType));
  };

  /*--- L/D/U index structures and value arrays. ---*/
  {
    const auto& pat = geometry->GetSparsePattern(type, 0);
    mat.row_ptr_l = pat.l.outerPtr();
    mat.col_ind_l = pat.l.innerIdx();
    mat.nnz_l = pat.l.getNumNonZeros();
    mat.row_ptr_u = pat.u.outerPtr();
    mat.col_ind_u = pat.u.innerIdx();
    mat.nnz_u = pat.u.getNumNonZeros();
  }
  allocAndInit(mat.d, nPoint * nVar * nEqn);

  if (q_lus_needed) {
    /*--- Q_LU_SGS: no full-precision L/U; off-diagonal blocks live in quantized storage.
     *    L/U are quantized on-the-fly during assembly; diagonal is quantized in Build step. ---*/
#ifndef CODI_REVERSE_TYPE
    quantized_mode = true;
#endif
    auto allocQ = [](QuantType*& ptr, unsigned long n) {
      ptr = MemoryAllocation::aligned_alloc<QuantType, true>(64, n * sizeof(QuantType));
    };
    allocQ(q_scale_l, mat.nnz_l * nVar);
    allocQ(q_blocks_l, mat.nnz_l * nVar * nEqn);
    allocQ(q_scale_u, mat.nnz_u * nVar);
    allocQ(q_blocks_u, mat.nnz_u * nVar * nEqn);
    allocQ(q_scale_d, nPoint * nVar);
    allocQ(q_blocks_d, nPoint * nVar * nEqn);
  } else {
    allocAndInit(mat.l, mat.nnz_l * nVar * nEqn);
    allocAndInit(mat.u, mat.nnz_u * nVar * nEqn);
  }

  if (useCuda) {
    auto GPUAllocAndInit = [](ScalarType*& ptr, unsigned long num) {
      ptr = GPUMemoryAllocation::gpu_alloc<ScalarType, true>(num * sizeof(ScalarType));
    };
    auto GPUAllocAndCopy = [](const su2uint*& ptr, const su2uint* src_ptr, unsigned long num) {
      ptr = GPUMemoryAllocation::gpu_alloc_cpy<su2uint>(src_ptr, num * sizeof(su2uint));
    };
    GPUAllocAndInit(gpu.d, nPoint * nVar * nEqn);
    GPUAllocAndInit(gpu.l, mat.nnz_l * nVar * nEqn);
    GPUAllocAndInit(gpu.u, mat.nnz_u * nVar * nEqn);
    GPUAllocAndCopy(gpu.row_ptr_l, mat.row_ptr_l, nPointDomain + 1);
    GPUAllocAndCopy(gpu.col_ind_l, mat.col_ind_l, mat.nnz_l);
    GPUAllocAndCopy(gpu.row_ptr_u, mat.row_ptr_u, nPointDomain + 1);
    GPUAllocAndCopy(gpu.col_ind_u, mat.col_ind_u, mat.nnz_u);
  }

  if (type == ConnectivityType::FiniteVolume) {
    edge_ptr_l = geometry->GetUToLTransposeSparsePatternMap(type).data();
  }
  if (needTranspPtr) {
    l_to_u_transp = geometry->GetLToUTransposeSparsePatternMap(type).data();
    u_to_l_transp = geometry->GetUToLTransposeSparsePatternMap(type).data();
  }

  /*--- Get ILU sparse pattern, if fill is 0 no new data is allocated. --*/

  if (ilu_needed) {
    ilu_fill_in = config->GetLinear_Solver_ILU_n();

    const auto& pat_ilu = geometry->GetSparsePattern(type, ilu_fill_in);
    ilu.row_ptr_l = pat_ilu.l.outerPtr();
    ilu.col_ind_l = pat_ilu.l.innerIdx();
    ilu.nnz_l = pat_ilu.l.getNumNonZeros();
    ilu.row_ptr_u = pat_ilu.u.outerPtr();
    ilu.col_ind_u = pat_ilu.u.innerIdx();
    ilu.nnz_u = pat_ilu.u.getNumNonZeros();

    if (omp_get_max_threads() > 1 && config->GetLinear_Solver_ILU_levels()) {
      levels_ilu = computeLevels(pat_ilu.l);
    }
  }

  /*--- Preconditioners. ---*/

  if (ilu_needed) {
    allocAndInit(ilu.l, ilu.nnz_l * nVar * nEqn);
    allocAndInit(ilu.d, nPointDomain * nVar * nEqn);
    allocAndInit(ilu.u, ilu.nnz_u * nVar * nEqn);
  }

  if (diag_needed) allocAndInit(invM, nPointDomain * nVar * nEqn);

  if (useCuda && diag_needed) {
    d_invM = GPUMemoryAllocation::gpu_alloc<ScalarType, true>(nPointDomain * nVar * nEqn * sizeof(ScalarType));
  }

  /*--- Thread parallel initialization. ---*/

  int num_threads = omp_get_max_threads();

  /*--- Set suitable chunk sizes for light static for loops, and heavy
   dynamic ones, such that threads are approximately evenly loaded. ---*/
  omp_light_size = computeStaticChunkSize(nPoint * nVar * nEqn, num_threads, OMP_MAX_SIZE_L);
  omp_heavy_size = computeStaticChunkSize(nPointDomain, num_threads, OMP_MAX_SIZE_H);

  omp_num_parts = config->GetLinear_Solver_Prec_Threads();
  if (omp_num_parts == 0) omp_num_parts = num_threads;

  /*--- This is akin to the row_ptr. ---*/
  omp_partitions = new unsigned long[omp_num_parts + 1];
  for (unsigned long i = 0; i <= omp_num_parts; ++i) omp_partitions[i] = nPointDomain;

  /*--- Work estimate based on non-zeros to produce balanced partitions. ---*/

  /*--- Cumulative nnz up to row iPoint for the preconditioner's LDU pattern. ---*/
  auto nnz_up_to = [&](unsigned long iPoint) -> unsigned long {
    if (ilu_needed) return ilu.row_ptr_l[iPoint] + iPoint + ilu.row_ptr_u[iPoint];
    return mat.row_ptr_l[iPoint] + iPoint + mat.row_ptr_u[iPoint];
  };
  const auto nnz_prec = nnz_up_to(nPointDomain);
  const auto nnz_per_part = roundUpDiv(nnz_prec, omp_num_parts);

  for (auto iPoint = 0ul, part = 0ul; iPoint < nPointDomain; ++iPoint) {
    if (nnz_up_to(iPoint) >= part * nnz_per_part) omp_partitions[part++] = iPoint;
  }

  for (unsigned long thread = 0; thread < omp_num_parts; ++thread) {
    const auto begin = omp_partitions[thread];
    const auto end = omp_partitions[thread + 1];
    if (begin == end) {
      cout << "WARNING: Redundant thread has been detected. Performance could be impacted due to low number of nodes "
              "per thread."
           << endl;
      break;
    }
  }

  /*--- Generate MKL Kernels ---*/

#ifdef USE_MKL
  using mkl = mkl_jit_wrapper<ScalarType>;
  mkl::create_gemm(&MatrixMatrixProductJitter, MKL_ROW_MAJOR, MKL_NOTRANS, MKL_NOTRANS, nVar, nVar, nVar, 1.0, nVar,
                   nVar, 0.0, nVar);
  MatrixMatrixProductKernel = mkl::get_gemm(MatrixMatrixProductJitter);

  mkl::create_gemm(&MatrixVectorProductJitterBetaZero, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, 1, nVar, nEqn, 1.0, 1,
                   nEqn, 0.0, 1);
  MatrixVectorProductKernelBetaZero = mkl::get_gemm(MatrixVectorProductJitterBetaZero);

  mkl::create_gemm(&MatrixVectorProductJitterBetaOne, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, 1, nVar, nEqn, 1.0, 1,
                   nEqn, 1.0, 1);
  MatrixVectorProductKernelBetaOne = mkl::get_gemm(MatrixVectorProductJitterBetaOne);

  mkl::create_gemm(&MatrixVectorProductJitterAlphaMinusOne, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, 1, nVar, nEqn,
                   -1.0, 1, nEqn, 1.0, 1);
  MatrixVectorProductKernelAlphaMinusOne = mkl::get_gemm(MatrixVectorProductJitterAlphaMinusOne);
#endif
}

template <class T>
void CSysMatrixComms::Initiate(const CSysVector<T>& x, CGeometry* geometry, const CConfig* config,
                               MPI_QUANTITIES commType) {
  SU2_ZONE_SCOPED
  if (geometry->nP2PSend == 0) return;

  /*--- Local variables ---*/

  const unsigned short COUNT_PER_POINT = x.GetNVar();
  const auto MPI_TYPE = geometry->GetCommType<T>();

  /*--- Create a boolean for reversing the order of comms. ---*/

  const bool reverse = (commType == MPI_QUANTITIES::SOLUTION_MATRIXTRANS);

  /*--- Set the size of the data packet and type depending on quantity. ---*/

  switch (commType) {
    case MPI_QUANTITIES::SOLUTION_MATRIX:
    case MPI_QUANTITIES::SOLUTION_MATRIXTRANS:
      break;
    default:
      SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.", CURRENT_FUNCTION);
      break;
  }

  /*--- Check to make sure we have created a large enough buffer
   for these comms during preprocessing. This is only for the su2double
   buffer. It will be reallocated whenever we find a larger count
   per point. After the first cycle of comms, this should be inactive. ---*/

  geometry->AllocateP2PComms(COUNT_PER_POINT);

  /*--- Load the specified quantity from the solver into the generic
   communication buffer in the geometry class. ---*/

  /*--- Post all non-blocking recvs first before sends. ---*/

  geometry->PostP2PRecvs(geometry, config, MPI_TYPE, COUNT_PER_POINT, reverse);

  for (auto iMessage = 0; iMessage < geometry->nP2PSend; iMessage++) {
    switch (commType) {
      case MPI_QUANTITIES::SOLUTION_MATRIX: {
        auto* bufDSend = geometry->GetP2PSendBuf<T>();

        /*--- Get the offset for the start of this message. ---*/

        const auto msg_offset = geometry->nPoint_P2PSend[iMessage];

        /*--- Total count can include multiple pieces of data per point. ---*/

        const auto nSend = (geometry->nPoint_P2PSend[iMessage + 1] - geometry->nPoint_P2PSend[iMessage]);

        SU2_OMP_FOR_STAT(CSysMatrix<T>::OMP_MIN_SIZE)
        for (auto iSend = 0; iSend < nSend; iSend++) {
          /*--- Get the local index for this communicated data. ---*/

          const auto iPoint = geometry->Local_Point_P2PSend[msg_offset + iSend];

          /*--- Compute the offset in the recv buffer for this point. ---*/

          const auto buf_offset = (msg_offset + iSend) * COUNT_PER_POINT;

          /*--- Load the buffer with the data to be sent. ---*/

          for (auto iVar = 0ul; iVar < x.GetNVar(); iVar++) bufDSend[buf_offset + iVar] = x(iPoint, iVar);
        }
        END_SU2_OMP_FOR
        break;
      }

      case MPI_QUANTITIES::SOLUTION_MATRIXTRANS: {
        /*--- We are going to communicate in reverse, so we use the
         recv buffer for the send instead. Also, all of the offsets
         and counts are derived from the recv data structures. ---*/
        auto* bufDSend = geometry->GetP2PRecvBuf<T>();

        /*--- Get the offset for the start of this message. ---*/

        const auto msg_offset = geometry->nPoint_P2PRecv[iMessage];

        /*--- Total count can include multiple pieces of data per point. ---*/

        const auto nSend = (geometry->nPoint_P2PRecv[iMessage + 1] - geometry->nPoint_P2PRecv[iMessage]);

        SU2_OMP_FOR_STAT(CSysMatrix<T>::OMP_MIN_SIZE)
        for (auto iSend = 0; iSend < nSend; iSend++) {
          /*--- Get the local index for this communicated data. Here we
           again use the recv structure to find the send point, since
           the usual recv points are now the senders in reverse mode. ---*/

          const auto iPoint = geometry->Local_Point_P2PRecv[msg_offset + iSend];

          /*--- Compute the offset in the recv buffer for this point. ---*/

          const auto buf_offset = (msg_offset + iSend) * COUNT_PER_POINT;

          /*--- Load the buffer with the data to be sent. ---*/

          for (auto iVar = 0ul; iVar < x.GetNVar(); iVar++) bufDSend[buf_offset + iVar] = x(iPoint, iVar);
        }
        END_SU2_OMP_FOR
        break;
      }

      default:
        SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.", CURRENT_FUNCTION);
        break;
    }

    /*--- Launch the point-to-point MPI send for this message. ---*/

    geometry->PostP2PSends(geometry, config, MPI_TYPE, COUNT_PER_POINT, iMessage, reverse);
  }
}

template <class T>
void CSysMatrixComms::Complete(CSysVector<T>& x, CGeometry* geometry, const CConfig* config, MPI_QUANTITIES commType) {
  SU2_ZONE_SCOPED
  if (geometry->nP2PRecv == 0) return;

  /*--- Local variables ---*/

  const unsigned short COUNT_PER_POINT = x.GetNVar();

  /*--- Global status so all threads can see the result of Waitany. ---*/
  static typename SelectMPIWrapper<T>::W::Status status;
  int ind;

  /*--- Store the data that was communicated into the appropriate
   location within the local class data structures. ---*/

  for (auto iMessage = 0; iMessage < geometry->nP2PRecv; iMessage++) {
    /*--- For efficiency, recv the messages dynamically based on
     the order they arrive. ---*/

    SU2_OMP_SAFE_GLOBAL_ACCESS(
        SelectMPIWrapper<T>::W::Waitany(geometry->nP2PRecv, geometry->GetP2PRecvReq<T>(), &ind, &status);)

    /*--- Once we have recv'd a message, get the source rank. ---*/

    const auto source = status.MPI_SOURCE;

    switch (commType) {
      case MPI_QUANTITIES::SOLUTION_MATRIX: {
        const auto* bufDRecv = geometry->GetP2PRecvBuf<T>();

        /*--- We know the offsets based on the source rank. ---*/

        const auto jRecv = geometry->P2PRecv2Neighbor[source];

        /*--- Get the offset for the start of this message. ---*/

        const auto msg_offset = geometry->nPoint_P2PRecv[jRecv];

        /*--- Get the number of packets to be received in this message. ---*/

        const auto nRecv = (geometry->nPoint_P2PRecv[jRecv + 1] - geometry->nPoint_P2PRecv[jRecv]);

        SU2_OMP_FOR_STAT(CSysMatrix<T>::OMP_MIN_SIZE)
        for (auto iRecv = 0; iRecv < nRecv; iRecv++) {
          /*--- Get the local index for this communicated data. ---*/

          const auto iPoint = geometry->Local_Point_P2PRecv[msg_offset + iRecv];

          /*--- Compute the offset in the recv buffer for this point. ---*/

          const auto buf_offset = (msg_offset + iRecv) * COUNT_PER_POINT;

          /*--- Store the data correctly depending on the quantity. ---*/

          for (auto iVar = 0ul; iVar < x.GetNVar(); iVar++)
            x(iPoint, iVar) = CSysMatrix<T>::template ActiveAssign<T>(bufDRecv[buf_offset + iVar]);
        }
        END_SU2_OMP_FOR
        break;
      }

      case MPI_QUANTITIES::SOLUTION_MATRIXTRANS: {
        /*--- We are going to communicate in reverse, so we use the
         send buffer for the recv instead. Also, all of the offsets
         and counts are derived from the send data structures. ---*/

        const auto* bufDRecv = geometry->GetP2PSendBuf<T>();

        /*--- We know the offsets based on the source rank. ---*/

        const auto jRecv = geometry->P2PSend2Neighbor[source];

        /*--- Get the offset for the start of this message. ---*/

        const auto msg_offset = geometry->nPoint_P2PSend[jRecv];

        /*--- Get the number of packets to be received in this message. ---*/

        const auto nRecv = (geometry->nPoint_P2PSend[jRecv + 1] - geometry->nPoint_P2PSend[jRecv]);

        SU2_OMP_FOR_STAT(CSysMatrix<T>::OMP_MIN_SIZE)
        for (auto iRecv = 0; iRecv < nRecv; iRecv++) {
          /*--- Get the local index for this communicated data. ---*/

          const auto iPoint = geometry->Local_Point_P2PSend[msg_offset + iRecv];

          /*--- Compute the offset in the recv buffer for this point. ---*/

          const auto buf_offset = (msg_offset + iRecv) * COUNT_PER_POINT;

          /*--- Update receiving point. ---*/

          for (auto iVar = 0ul; iVar < x.GetNVar(); iVar++)
            x(iPoint, iVar) += CSysMatrix<T>::template ActiveAssign<T>(bufDRecv[buf_offset + iVar]);
        }
        END_SU2_OMP_FOR
        break;
      }

      default:
        SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.", CURRENT_FUNCTION);
        break;
    }
  }

  /*--- Verify that all non-blocking point-to-point sends have finished.
   Note that this should be satisfied, as we have received all of the
   data in the loop above at this point. ---*/

#ifdef HAVE_MPI
  SU2_OMP_SAFE_GLOBAL_ACCESS(
      SelectMPIWrapper<T>::W::Waitall(geometry->nP2PSend, geometry->GetP2PSendReq<T>(), MPI_STATUS_IGNORE);)
#endif
}

template <class ScalarType>
void CSysMatrix<ScalarType>::QuantizeBlock(const ScalarType* blk, QuantType* qs, QuantType* qv) const {
  EncodeQuantBlock([&](unsigned long r, unsigned long c) { return blk[r * nVar + c]; }, qs, qv, nVar);
}

template <class ScalarType>
void CSysMatrix<ScalarType>::QuantizeDiagonalBlocks() {
  SU2_ZONE_SCOPED

  if (quantized_mode) {
    /*--- Q_LU_SGS: L/U were quantized during assembly; only the diagonal needs quantization now. ---*/
    SU2_OMP_FOR_DYN(omp_heavy_size)
    for (auto i = 0ul; i < nPointDomain; ++i)
      QuantizeBlock(&mat.d[i * nVar * nVar], &q_scale_d[i * nVar], &q_blocks_d[i * nVar * nVar]);
    END_SU2_OMP_FOR
  }
}

template <class ScalarType>
void CSysMatrix<ScalarType>::SetValZero() {
  SU2_ZONE_SCOPED
  const auto nThreads = static_cast<unsigned long>(omp_get_num_threads());
  const auto iThread = static_cast<unsigned long>(omp_get_thread_num());

  auto zeroChunk = [&](auto* arr, unsigned long n) {
    if (n == 0) return;
    const auto chunk = roundUpDiv(n, nThreads);
    const auto begin = min<size_t>(chunk * iThread, n);
    const auto mySize = min<size_t>(chunk, n - begin) * sizeof(std::remove_pointer_t<decltype(arr)>);
    if (mySize) memset(&arr[begin], 0, mySize);
  };
  zeroChunk(mat.d, nPoint * nVar * nEqn);
  if (!quantized_mode) {
    zeroChunk(mat.l, mat.nnz_l * nVar * nEqn);
    zeroChunk(mat.u, mat.nnz_u * nVar * nEqn);
  } else {
    zeroChunk(q_scale_l, mat.nnz_l * nVar);
    zeroChunk(q_scale_u, mat.nnz_l * nVar);
    zeroChunk(q_blocks_l, mat.nnz_l * nVar * nEqn);
    zeroChunk(q_blocks_u, mat.nnz_u * nVar * nEqn);
  }
  SU2_OMP_BARRIER
}

template <class ScalarType>
void CSysMatrix<ScalarType>::SetValDiagonalZero() {
  SU2_ZONE_SCOPED
  SU2_OMP_FOR_STAT(omp_heavy_size)
  for (auto iVar = 0ul; iVar < nPointDomain * nVar * nEqn; ++iVar) mat.d[iVar] = 0;
  END_SU2_OMP_FOR
}

template <class ScalarType>
void CSysMatrix<ScalarType>::GaussElimination(ScalarType* matrix, ScalarType* vec) const {
#ifdef USE_MKL_LAPACK
  // With MKL_DIRECT_CALL enabled, this is significantly faster than native code on Intel Architectures.
  lapack_int ipiv[MAXNVAR];
  if constexpr (std::is_same_v<ScalarType, double>) {
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, nVar, nVar, matrix, nVar, ipiv);
    LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', nVar, 1, matrix, nVar, ipiv, vec, 1);
  } else {
    static_assert(std::is_same_v<ScalarType, float>, "ScalarType not handled");
    LAPACKE_sgetrf(LAPACK_ROW_MAJOR, nVar, nVar, matrix, nVar, ipiv);
    LAPACKE_sgetrs(LAPACK_ROW_MAJOR, 'N', nVar, 1, matrix, nVar, ipiv, vec, 1);
  }
#else
#define A(I, J) matrix[(I)*nVar + (J)]

  /*--- Transform system in Upper Matrix ---*/

  for (auto iVar = 1ul; iVar < nVar; iVar++) {
    for (auto jVar = 0ul; jVar < iVar; jVar++) {
      /*--- Regularize pivot if too small to prevent divide-by-zero ---*/
      RegularizePivot(A(jVar, jVar), jVar, jVar, "DEBUG GaussElimination");

      ScalarType weight = A(iVar, jVar) / A(jVar, jVar);

      for (auto kVar = jVar; kVar < nVar; kVar++) A(iVar, kVar) -= weight * A(jVar, kVar);
      vec[iVar] -= weight * vec[jVar];
    }
  }

  /*--- Backwards substitution ---*/
  for (auto iVar = nVar; iVar > 0ul;) {
    iVar--;  // unsigned type
    for (auto jVar = iVar + 1; jVar < nVar; jVar++) vec[iVar] -= A(iVar, jVar) * vec[jVar];

    /*--- Regularize diagonal if too small ---*/
    RegularizePivot(A(iVar, iVar), iVar, iVar, "DEBUG GaussElimination backsubst");

    vec[iVar] /= A(iVar, iVar);
  }
#undef A
#endif
}

template <class ScalarType>
void CSysMatrix<ScalarType>::MatrixInverse(ScalarType* matrix, ScalarType* inverse) const {
  /*--- This is a generalization of Gaussian elimination for multiple rhs' (the basis vectors).
   We could call "GaussElimination" multiple times or fully generalize it for multiple rhs,
   the performance of both routines would suffer in both cases without the use of exotic templating.
   And so it feels reasonable to have some duplication here. ---*/

  assert((matrix != inverse) && "Output cannot be the same as the input.");

#define M(I, J) inverse[(I)*nVar + (J)]

  /*--- Initialize the inverse with the identity. ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++)
    for (auto jVar = 0ul; jVar < nVar; jVar++) M(iVar, jVar) = ScalarType(iVar == jVar);

      /*--- Inversion ---*/
#ifdef USE_MKL_LAPACK
  // With MKL_DIRECT_CALL enabled, this is significantly faster than native code on Intel Architectures.
  lapack_int ipiv[MAXNVAR];
  if constexpr (std::is_same_v<ScalarType, double>) {
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, nVar, nVar, matrix, nVar, ipiv);
    LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', nVar, nVar, matrix, nVar, ipiv, inverse, nVar);
  } else {
    static_assert(std::is_same_v<ScalarType, float>, "ScalarType not handled");
    LAPACKE_sgetrf(LAPACK_ROW_MAJOR, nVar, nVar, matrix, nVar, ipiv);
    LAPACKE_sgetrs(LAPACK_ROW_MAJOR, 'N', nVar, nVar, matrix, nVar, ipiv, inverse, nVar);
  }
#else
#define A(I, J) matrix[(I)*nVar + (J)]

  /*--- Transform system in Upper Matrix ---*/
  for (auto iVar = 1ul; iVar < nVar; iVar++) {
    for (auto jVar = 0ul; jVar < iVar; jVar++) {
      /*--- Regularize pivot if too small to prevent divide-by-zero ---*/
      RegularizePivot(A(jVar, jVar), jVar, jVar, "MatrixInverse");

      ScalarType weight = A(iVar, jVar) / A(jVar, jVar);
      for (auto kVar = jVar; kVar < nVar; kVar++) A(iVar, kVar) -= weight * A(jVar, kVar);

      /*--- at this stage M is lower triangular so not all cols need updating ---*/
      for (auto kVar = 0ul; kVar <= jVar; kVar++) M(iVar, kVar) -= weight * M(jVar, kVar);
    }
  }

  /*--- Backwards substitution ---*/
  for (auto iVar = nVar; iVar > 0ul;) {
    iVar--;  // unsigned type
    for (auto jVar = iVar + 1; jVar < nVar; jVar++)
      for (auto kVar = 0ul; kVar < nVar; kVar++) M(iVar, kVar) -= A(iVar, jVar) * M(jVar, kVar);

    /*--- Regularize diagonal if too small ---*/
    RegularizePivot(A(iVar, iVar), iVar, iVar, "DEBUG MatrixInverse backsubst");

    for (auto kVar = 0ul; kVar < nVar; kVar++) {
      M(iVar, kVar) /= A(iVar, iVar);
    }
  }
#undef A
#endif
#undef M
}

template <class ScalarType>
void CSysMatrix<ScalarType>::DeleteValsRowi(unsigned long block_i, unsigned long row) {
  SU2_ZONE_SCOPED
  const auto blkSz = nVar * nEqn;

  auto* d = &mat.d[block_i * blkSz];
  for (auto iVar = 0u; iVar < nEqn; iVar++) d[row * nEqn + iVar] = 0.0;
  d[row * nEqn + row] = 1.0;

  if (quantized_mode) {
    for (auto k = mat.row_ptr_l[block_i]; k < mat.row_ptr_l[block_i + 1]; ++k) {
      for (auto iVar = 0u; iVar < nEqn; iVar++) q_blocks_l[k * blkSz + row * nEqn + iVar] = 0;
    }
    for (auto k = mat.row_ptr_u[block_i]; k < mat.row_ptr_u[block_i + 1]; ++k) {
      for (auto iVar = 0u; iVar < nEqn; iVar++) q_blocks_u[k * blkSz + row * nEqn + iVar] = 0;
    }
  } else {
    for (auto k = mat.row_ptr_l[block_i]; k < mat.row_ptr_l[block_i + 1]; ++k) {
      for (auto iVar = 0u; iVar < nEqn; iVar++) mat.l[k * blkSz + row * nEqn + iVar] = 0;
    }
    for (auto k = mat.row_ptr_u[block_i]; k < mat.row_ptr_u[block_i + 1]; ++k) {
      for (auto iVar = 0u; iVar < nEqn; iVar++) mat.u[k * blkSz + row * nEqn + iVar] = 0;
    }
  }
}

template <class ScalarType>
void CSysMatrix<ScalarType>::MatrixVectorProduct(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                                 CGeometry* geometry, const CConfig* config) const {
  SU2_ZONE_SCOPED
  /*--- Some checks for consistency between CSysMatrix and the CSysVector<ScalarType>s ---*/
#ifndef NDEBUG
  if ((nEqn != vec.GetNVar()) || (nVar != prod.GetNVar())) {
    SU2_MPI::Error("nVar values incompatible.", CURRENT_FUNCTION);
  }
  if (nPoint != prod.GetNBlk()) {
    SU2_MPI::Error("nPoint and nBlk values incompatible.", CURRENT_FUNCTION);
  }
#endif

  /*--- OpenMP parallelization. First need to make view of vectors
   *    consistent, a barrier is implicit at the end of FOR section
   *    (and it is required before master thread communicates). ---*/

  SU2_OMP_BARRIER

  if (quantized_mode) {
    SU2_OMP_FOR_DYN(omp_heavy_size)
    for (auto row_i = 0ul; row_i < nPointDomain; row_i++) {
      QuantizedRowProduct(vec, row_i, &prod[row_i * nVar]);
    }
    END_SU2_OMP_FOR
  } else {
    SU2_OMP_FOR_DYN(omp_heavy_size)
    for (auto row_i = 0ul; row_i < nPointDomain; row_i++) {
      RowProduct(vec, row_i, &prod[row_i * nVar]);
    }
    END_SU2_OMP_FOR
  }

  /*--- MPI Parallelization. ---*/

  CSysMatrixComms::Initiate(prod, geometry, config);
  CSysMatrixComms::Complete(prod, geometry, config);
}

template <class ScalarType>
void CSysMatrix<ScalarType>::BuildJacobiPreconditioner() {
  SU2_ZONE_SCOPED

  if (useCuda) {
#ifdef HAVE_CUDA
    BuildJacobiPreconditionerGPU();
    return;
#else
    SU2_MPI::Error(
        "\nError in building Jacobi preconditioner\nENABLE_CUDA is set to YES\nPlease compile with CUDA options "
        "enabled in Meson to access GPU Functions",
        CURRENT_FUNCTION);
#endif
  }

  /*--- Build Jacobi preconditioner (M = D), compute and store the inverses of the diagonal blocks. ---*/
  SU2_OMP_FOR_DYN(omp_heavy_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
    InverseDiagonalBlock(iPoint, &(invM[iPoint * nVar * nVar]));
  END_SU2_OMP_FOR
}

template <class ScalarType>
void CSysMatrix<ScalarType>::ComputeJacobiPreconditioner(const CSysVector<ScalarType>& vec,
                                                         CSysVector<ScalarType>& prod, CGeometry* geometry,
                                                         const CConfig* config) const {
  SU2_ZONE_SCOPED

  if (config->GetCUDA()) {
#ifdef HAVE_CUDA
    ComputeJacobiPreconditionerGPU(vec, prod, geometry, config);
    return;
#else
    SU2_MPI::Error(
        "\nError in applying Jacobi preconditioner\nENABLE_CUDA is set to YES\nPlease compile with CUDA options "
        "enabled in Meson to access GPU Functions",
        CURRENT_FUNCTION);
#endif
  }

  /*--- Apply Jacobi preconditioner, y = D^{-1} * x, the inverse of the diagonal is already known. ---*/
  SU2_OMP_BARRIER
  SU2_OMP_FOR_DYN(omp_heavy_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
    MatrixVectorProduct(&(invM[iPoint * nVar * nVar]), &vec[iPoint * nVar], &prod[iPoint * nVar]);
  END_SU2_OMP_FOR

  /*--- MPI Parallelization ---*/
  CSysMatrixComms::Initiate(prod, geometry, config);
  CSysMatrixComms::Complete(prod, geometry, config);
}

template <class ScalarType>
void CSysMatrix<ScalarType>::BuildILUPreconditioner() {
  SU2_ZONE_SCOPED
  const auto blockSize = nVar * nVar;
  ScalarType Lij[MAXNVAR * MAXNVAR], Lij_Ujk[MAXNVAR * MAXNVAR];

  /*--- Helper to copy block matrix to compute factorization in-place. ---*/
  auto InitIluRow = [&](const auto iPoint) {
    MatrixCopy(&mat.d[iPoint * blockSize], &ilu.d[iPoint * blockSize]);

    if (ilu_fill_in == 0) {
      /*--- ILU0: Same sparse pattern, copy L and U blocks directly. ---*/
      auto copy = [&](const su2uint* row_ptr, const ScalarType* mat, ScalarType* ilu) {
        const unsigned long begin = row_ptr[iPoint] * blockSize;
        const unsigned long end = row_ptr[iPoint + 1] * blockSize;
        SU2_OMP_SIMD
        for (auto k = begin; k < end; ++k) ilu[k] = mat[k];
      };
      copy(ilu.row_ptr_l, mat.l, ilu.l);
      copy(ilu.row_ptr_u, mat.u, ilu.u);
      return;
    }
    /*--- ILUn: Merge-scan L and U via shared lambda. ---*/
    auto scatterPart = [&](const su2uint* mat_row_ptr, const su2uint* mat_col_ind, const ScalarType* mat_vals,
                           const su2uint* ilu_row_ptr, const su2uint* ilu_col_ind, ScalarType* ilu_vals) {
      auto km = mat_row_ptr[iPoint], km_end = mat_row_ptr[iPoint + 1];
      for (auto k = ilu_row_ptr[iPoint]; k < ilu_row_ptr[iPoint + 1]; ++k) {
        const auto jPoint = ilu_col_ind[k];
        while (km < km_end && mat_col_ind[km] < jPoint) ++km;
        if (km < km_end && mat_col_ind[km] == jPoint) {
          MatrixCopy(&mat_vals[km * blockSize], &ilu_vals[k * blockSize]);
        } else {
          ZeroMatrix(&ilu_vals[k * blockSize]);
        }
      }
    };
    scatterPart(mat.row_ptr_l, mat.col_ind_l, mat.l, ilu.row_ptr_l, ilu.col_ind_l, ilu.l);
    scatterPart(mat.row_ptr_u, mat.col_ind_u, mat.u, ilu.row_ptr_u, ilu.col_ind_u, ilu.u);
  };

  /*--- Update one row of the LU matrix. ---*/
  auto BuildIluRow = [&](const auto iPoint, const auto begin, const auto end) {
    /*--- For this row (unknown), loop over its lower diagonal entries. ---*/

    for (auto kl = ilu.row_ptr_l[iPoint]; kl < ilu.row_ptr_l[iPoint + 1]; ++kl) {
      /*--- jPoint is the column index (jPoint < iPoint). ---*/

      const auto jPoint = ilu.col_ind_l[kl];

      /*--- We only care about the sub matrix within "begin" and "end-1". ---*/

      if (jPoint < begin) continue;

      /*--- Multiply the block by the inverse of the corresponding diagonal block. ---*/

      auto* Block_ij = &ilu.l[kl * blockSize];
      const auto* invUjj = &ilu.d[jPoint * blockSize];
      MatrixMatrixProduct(Block_ij, invUjj, Lij);

      /*--- Lij holds Aij*inv(Ujj). Jump to the upper part of the jPoint row. ---*/

      for (auto ku = ilu.row_ptr_u[jPoint]; ku < ilu.row_ptr_u[jPoint + 1]; ++ku) {
        /*--- Get the column index (kPoint > jPoint). ---*/

        const auto kPoint = ilu.col_ind_u[ku];
        if (kPoint >= end) break;

        /*--- If Aik exists, update it: Aik -= Lij * Ujk ---*/

        auto* Block_ik = GetBlock_ILUMatrix(iPoint, kPoint);
        if (Block_ik == nullptr) continue;

        const auto* Ujk = &ilu.u[ku * blockSize];
        MatrixMatrixProduct(Lij, Ujk, Lij_Ujk);
        MatrixSubtraction(Block_ik, Lij_Ujk, Block_ik);
      }

      /*--- Store Lij in the lower triangular part. ---*/
      MatrixCopy(Lij, Block_ij);
    }

    /*--- Invert the diagonal entry, Uii, for the next rows. ---*/
    InvertDiagonalBlockILUMatrix(iPoint);
  };

  if (levels_ilu.empty()) {
    /*--- Each OMP thread will work on the submatrix defined from
     * row/col "begin" to row/col "end-1" (i.e. the range [begin,end[).
     * Which is exactly what the MPI-only implementation does. ---*/
    SU2_OMP_FOR_STAT(1)
    for (unsigned long thread = 0; thread < omp_num_parts; ++thread) {
      const auto begin = omp_partitions[thread];
      const auto end = omp_partitions[thread + 1];
      if (begin == end) continue;

      InitIluRow(begin);
      InvertDiagonalBlockILUMatrix(begin);

      for (auto iPoint = begin + 1; iPoint < end; iPoint++) {
        InitIluRow(iPoint);
        BuildIluRow(iPoint, begin, end);
      }
    }
    END_SU2_OMP_FOR
  } else {
    /*--- OMP threads work on each level together before moving to the next.
     * Levels are determined such that rows in a level only depend on rows
     * from previous levels. ---*/

    SU2_OMP_FOR_(schedule(static))
    for (auto k = 0ul; k < levels_ilu.getNumNonZeros(0); ++k) {
      const auto iPoint = levels_ilu.getInnerIdx(0, k);
      InitIluRow(iPoint);
      InvertDiagonalBlockILUMatrix(iPoint);
    }
    END_SU2_OMP_FOR

    for (auto level = 1ul; level < levels_ilu.getOuterSize(); ++level) {
      SU2_OMP_FOR_(schedule(static))
      for (auto k = 0ul; k < levels_ilu.getNumNonZeros(level); ++k) {
        const auto iPoint = levels_ilu.getInnerIdx(level, k);
        InitIluRow(iPoint);
        BuildIluRow(iPoint, 0ul, nPointDomain);
      }
      END_SU2_OMP_FOR
    }
  }
}

template <class ScalarType>
void CSysMatrix<ScalarType>::ComputeILUPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                                      CGeometry* geometry, const CConfig* config) const {
  SU2_ZONE_SCOPED
  /*--- Coherent view of vectors. ---*/
  SU2_OMP_BARRIER

  const auto blockSize = nVar * nVar;

  /*--- Forward solve the system using the lower matrix entries
   * that were computed and stored during the ILU preprocessing.
   * Note that we are overwriting the prod vector as we go. ---*/

  auto ForwardSolve = [&](const auto iPoint, const auto begin) {
    for (auto iVar = 0ul; iVar < nVar; ++iVar) prod(iPoint, iVar) = vec(iPoint, iVar);

    for (auto kl = ilu.row_ptr_l[iPoint]; kl < ilu.row_ptr_l[iPoint + 1]; ++kl) {
      const auto jPoint = ilu.col_ind_l[kl];
      if (jPoint < begin) continue;
      MatrixVectorProductSub(&ilu.l[kl * blockSize], &prod[jPoint * nVar], &prod[iPoint * nVar]);
    }
  };

  auto BackwardSolve = [&](const auto iPoint, const auto end) {
    ScalarType aux_vec[MAXNVAR];
    for (auto iVar = 0ul; iVar < nVar; ++iVar) aux_vec[iVar] = prod(iPoint, iVar);

    const auto* invUii = &ilu.d[iPoint * blockSize];

    for (auto ku = ilu.row_ptr_u[iPoint]; ku < ilu.row_ptr_u[iPoint + 1]; ++ku) {
      const auto jPoint = ilu.col_ind_u[ku];
      if (jPoint >= end) break;
      const auto* Block_ij = &ilu.u[ku * blockSize];
      MatrixVectorProductSub(Block_ij, &prod[jPoint * nVar], aux_vec);
    }

    MatrixVectorProduct(invUii, aux_vec, &prod[iPoint * nVar]);
  };

  if (levels_ilu.empty()) {
    SU2_OMP_FOR_STAT(1)
    for (unsigned long thread = 0; thread < omp_num_parts; ++thread) {
      const auto begin = omp_partitions[thread];
      const auto end = omp_partitions[thread + 1];
      if (begin == end) continue;

      for (auto iPoint = begin; iPoint < end; ++iPoint) {
        ForwardSolve(iPoint, begin);
      }
      for (auto iPoint = end; iPoint > begin;) {
        --iPoint;  // unsigned type
        BackwardSolve(iPoint, end);
      }
    }
    END_SU2_OMP_FOR
  } else {
    for (auto level = 0ul; level < levels_ilu.getOuterSize(); ++level) {
      SU2_OMP_FOR_(schedule(static))
      for (auto k = 0ul; k < levels_ilu.getNumNonZeros(level); ++k) {
        const auto iPoint = levels_ilu.getInnerIdx(level, k);
        ForwardSolve(iPoint, 0ul);
      }
      END_SU2_OMP_FOR
    }
    for (auto level = levels_ilu.getOuterSize(); level > 0;) {
      --level;  // unsigned type
      SU2_OMP_FOR_(schedule(static))
      for (auto k = 0ul; k < levels_ilu.getNumNonZeros(level); ++k) {
        const auto iPoint = levels_ilu.getInnerIdx(level, k);
        BackwardSolve(iPoint, nPointDomain);
      }
      END_SU2_OMP_FOR
    }
  }

  /*--- MPI Parallelization ---*/

  CSysMatrixComms::Initiate(prod, geometry, config);
  CSysMatrixComms::Complete(prod, geometry, config);
}

template <class ScalarType>
void CSysMatrix<ScalarType>::ComputeLU_SGSPreconditioner(const CSysVector<ScalarType>& vec,
                                                         CSysVector<ScalarType>& prod, CGeometry* geometry,
                                                         const CConfig* config) const {
  SU2_ZONE_SCOPED
  /*--- First part of the symmetric iteration: (D+L).x* = b ---*/

  /*--- Coherent view of vectors. ---*/
  SU2_OMP_BARRIER

  /*--- OpenMP Parallelization ---*/
  SU2_OMP_FOR_STAT(1)
  for (unsigned long thread = 0; thread < omp_num_parts; ++thread) {
    const auto begin = omp_partitions[thread];
    const auto end = omp_partitions[thread + 1];
    if (begin == end) continue;

    /*--- Each thread will work on the submatrix defined from row/col "begin"
     *    to row/col "end-1", except the last thread that also considers halos.
     *    This is NOT exactly equivalent to the MPI implementation on the same
     *    number of domains, for that we would need to define "thread-halos". ---*/

    ScalarType low_prod[MAXNVAR];

    if (quantized_mode) {
      for (auto iPoint = begin; iPoint < end; ++iPoint) {
        auto idx = iPoint * nVar;
        QuantizedLowerProduct(prod, iPoint, begin, low_prod);
        VectorSubtraction(&vec[idx], low_prod, &prod[idx]);
        QuantizedGaussElimination(iPoint, &prod[idx]);
      }
    } else {
      for (auto iPoint = begin; iPoint < end; ++iPoint) {
        auto idx = iPoint * nVar;
        LowerProduct(prod, iPoint, begin, low_prod);         // Compute L.x*
        VectorSubtraction(&vec[idx], low_prod, &prod[idx]);  // Compute y = b - L.x*
        GaussElimination(iPoint, &prod[idx]);                // Solve D.x* = y
      }
    }
  }
  END_SU2_OMP_FOR

  /*--- MPI Parallelization ---*/

  CSysMatrixComms::Initiate(prod, geometry, config);
  CSysMatrixComms::Complete(prod, geometry, config);

  /*--- Second part of the symmetric iteration: (D+U).x_(1) = D.x* ---*/

  /*--- OpenMP Parallelization ---*/
  SU2_OMP_FOR_STAT(1)
  for (unsigned long thread = 0; thread < omp_num_parts; ++thread) {
    const auto begin = omp_partitions[thread];
    const auto row_end = omp_partitions[thread + 1];
    if (begin == row_end) continue;

    ScalarType up_prod[MAXNVAR], dia_prod[MAXNVAR];

    if (quantized_mode) {
      for (auto iPoint = row_end; iPoint > begin;) {
        iPoint--;
        auto idx = iPoint * nVar;
        QuantizedDiagonalProduct(prod, iPoint, dia_prod);
        QuantizedUpperProduct(prod, iPoint, row_end, up_prod);
        VectorSubtraction(dia_prod, up_prod, &prod[idx]);
        QuantizedGaussElimination(iPoint, &prod[idx]);
      }
    } else {
      for (auto iPoint = row_end; iPoint > begin;) {
        iPoint--;  // because of unsigned type
        auto idx = iPoint * nVar;
        DiagonalProduct(prod, iPoint, dia_prod);           // Compute D.x*
        UpperProduct(prod, iPoint, row_end, up_prod);      // Compute U.x_(n+1)
        VectorSubtraction(dia_prod, up_prod, &prod[idx]);  // Compute y = D.x*-U.x_(n+1)
        GaussElimination(iPoint, &prod[idx]);              // Solve D.x* = y
      }
    }
  }
  END_SU2_OMP_FOR

  /*--- MPI Parallelization ---*/

  CSysMatrixComms::Initiate(prod, geometry, config);
  CSysMatrixComms::Complete(prod, geometry, config);
}

template <class ScalarType>
void CSysMatrix<ScalarType>::BuildLineletPreconditioner(const CGeometry* geometry, const CConfig* config) {
  SU2_ZONE_SCOPED
  BuildJacobiPreconditioner();

  /*--- Allocate working vectors if not done yet. ---*/
  if (!LineletUpper.empty()) return;

  const auto nThreads = omp_get_max_threads();

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    const auto& li = geometry->GetLineletInfo(config);
    if (!li.linelets.empty()) {
      LineletUpper.resize(nThreads);
      LineletVector.resize(nThreads);
      LineletInvDiag.resize(nThreads);
    }
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  SU2_OMP_FOR_STAT(1)
  for (int iThread = 0; iThread < nThreads; ++iThread) {
    const auto size = CGeometry::CLineletInfo::MAX_LINELET_POINTS;
    LineletUpper[iThread].resize(size, nullptr);
    LineletVector[iThread].resize(size * nVar, 0.0);
    LineletInvDiag[iThread].resize(size * nVar * nVar, 0.0);
  }
  END_SU2_OMP_FOR
}

template <class ScalarType>
void CSysMatrix<ScalarType>::ComputeLineletPreconditioner(const CSysVector<ScalarType>& vec,
                                                          CSysVector<ScalarType>& prod, CGeometry* geometry,
                                                          const CConfig* config) const {
  SU2_ZONE_SCOPED
  /*--- Coherent view of vectors. ---*/
  SU2_OMP_BARRIER

  const auto& li = geometry->GetLineletInfo(config);

  /*--- Jacobi preconditioning where there are no linelets. ---*/

  SU2_OMP_FOR_(schedule(dynamic, omp_heavy_size) SU2_NOWAIT)
  for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++)
    if (li.lineletIdx[iPoint] == CGeometry::CLineletInfo::NO_LINELET)
      MatrixVectorProduct(&(invM[iPoint * nVar * nVar]), &vec[iPoint * nVar], &prod[iPoint * nVar]);
  END_SU2_OMP_FOR

  /*--- Solve the tridiagonal systems for the linelets. ---*/

  SU2_OMP_FOR_DYN(1)
  for (auto iLinelet = 0ul; iLinelet < li.linelets.size(); iLinelet++) {
    /*--- Get references to the working vectors allocated for this thread. ---*/

    const int thread = omp_get_thread_num();
    auto& lineletUpper = LineletUpper[thread];
    auto& lineletInvDiag = LineletInvDiag[thread];
    auto& lineletVector = LineletVector[thread];

    /*--- Initialize the solution vector with the rhs ---*/

    const auto nElem = li.linelets[iLinelet].size();

    for (auto iElem = 0ul; iElem < nElem; iElem++) {
      const auto iPoint = li.linelets[iLinelet][iElem];
      for (auto iVar = 0ul; iVar < nVar; iVar++) lineletVector[iElem * nVar + iVar] = vec[iPoint * nVar + iVar];
    }

    /*--- Forward pass, eliminate lower entries, modify diagonal and rhs. ---*/

    /*--- Small temporaries. ---*/
    ScalarType aux_block[MAXNVAR * MAXNVAR], aux_vector[MAXNVAR];

    /*--- Copy diagonal block for first point in this linelet. ---*/
    MatrixCopy(&mat.d[li.linelets[iLinelet][0] * nVar * nEqn], lineletInvDiag.data());

    for (auto iElem = 1ul; iElem < nElem; iElem++) {
      /*--- Setup pointers to required matrices and vectors ---*/
      const auto im1Point = li.linelets[iLinelet][iElem - 1];
      const auto iPoint = li.linelets[iLinelet][iElem];

      const auto* d = GetBlock(iPoint, iPoint);
      const auto* l = GetBlock(iPoint, im1Point);
      const auto* u = GetBlock(im1Point, iPoint);

      auto* inv_dm1 = &lineletInvDiag[(iElem - 1) * nVar * nVar];
      auto* d_prime = &lineletInvDiag[iElem * nVar * nVar];
      auto* b_prime = &lineletVector[iElem * nVar];

      /*--- Invert previous modified diagonal ---*/
      MatrixCopy(inv_dm1, aux_block);
      MatrixInverse(aux_block, inv_dm1);

      /*--- Left-multiply by lower block to obtain the weight ---*/
      MatrixMatrixProduct(l, inv_dm1, aux_block);

      /*--- Multiply weight by upper block to modify current diagonal ---*/
      MatrixMatrixProduct(aux_block, u, d_prime);
      MatrixSubtraction(d, d_prime, d_prime);

      /*--- Update the rhs ---*/
      MatrixVectorProduct(aux_block, &lineletVector[(iElem - 1) * nVar], aux_vector);
      VectorSubtraction(b_prime, aux_vector, b_prime);

      /*--- Cache upper block pointer for the backward substitution phase ---*/
      lineletUpper[iElem - 1] = u;
    }

    /*--- Backwards substitution, LineletVector becomes the solution ---*/

    /*--- x_n = d_n^{-1} * b_n ---*/
    GaussElimination(&lineletInvDiag[(nElem - 1) * nVar * nVar], &lineletVector[(nElem - 1) * nVar]);

    /*--- x_i = d_i^{-1}*(b_i - u_i*x_{i+1}) ---*/
    for (auto iElem = nElem - 1; iElem > 0; --iElem) {
      const auto* inv_dm1 = &lineletInvDiag[(iElem - 1) * nVar * nVar];
      MatrixVectorProduct(lineletUpper[iElem - 1], &lineletVector[iElem * nVar], aux_vector);
      VectorSubtraction(&lineletVector[(iElem - 1) * nVar], aux_vector, aux_vector);
      MatrixVectorProduct(inv_dm1, aux_vector, &lineletVector[(iElem - 1) * nVar]);
    }

    /*--- Copy results to product vector ---*/

    for (auto iElem = 0ul; iElem < nElem; iElem++) {
      const auto iPoint = li.linelets[iLinelet][iElem];
      for (auto iVar = 0ul; iVar < nVar; iVar++) prod[iPoint * nVar + iVar] = lineletVector[iElem * nVar + iVar];
    }
  }
  END_SU2_OMP_FOR

  /*--- MPI Parallelization ---*/

  CSysMatrixComms::Initiate(prod, geometry, config);
  CSysMatrixComms::Complete(prod, geometry, config);
}

template <class ScalarType>
void CSysMatrix<ScalarType>::ComputeResidual(const CSysVector<ScalarType>& sol, const CSysVector<ScalarType>& f,
                                             CSysVector<ScalarType>& res) const {
  SU2_ZONE_SCOPED
  SU2_OMP_BARRIER
  SU2_OMP_FOR_DYN(omp_heavy_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    ScalarType aux_vec[MAXNVAR];
    RowProduct(sol, iPoint, aux_vec);
    VectorSubtraction(aux_vec, &f[iPoint * nVar], &res[iPoint * nVar]);
  }
  END_SU2_OMP_FOR
}

template <class ScalarType>
template <class OtherType>
void CSysMatrix<ScalarType>::EnforceSolutionAtNode(const unsigned long node_i, const OtherType* x_i,
                                                   CSysVector<OtherType>& b) {
  SU2_ZONE_SCOPED
  /*--- Eliminate the row associated with node i (Block_ii = I and all other Block_ij = 0).
   *    To preserve eventual symmetry, also attempt to eliminate the column, if the sparse pattern is not
   *    symmetric the entire column may not be eliminated, the result (matrix and vector) is still correct.
   *    The vector is updated with the product of column i by the known (enforced) solution at node i. ---*/

  /*--- Visit off-diagonal columns (L then U; diagonal is handled by SetVal2Diag outside). ---*/
  auto processOffDiag = [&](unsigned long node_j) {
    auto bij = GetBlock(node_i, node_j);
    auto bji = GetBlock(node_j, node_i);
    if (bji == nullptr) {
      node_j = node_i;
      bji = bij;
    }
    for (auto iVar = 0ul; iVar < nVar; ++iVar) {
      for (auto jVar = 0ul; jVar < nVar; ++jVar) {
        b[node_j * nVar + iVar] -= bji[iVar * nVar + jVar] * x_i[jVar];
        bij[iVar * nVar + jVar] = bji[iVar * nVar + jVar] = 0.0;
      }
    }
  };
  for (auto k = mat.row_ptr_l[node_i]; k < mat.row_ptr_l[node_i + 1]; ++k) processOffDiag(mat.col_ind_l[k]);
  for (auto k = mat.row_ptr_u[node_i]; k < mat.row_ptr_u[node_i + 1]; ++k) processOffDiag(mat.col_ind_u[k]);

  /*--- Set the diagonal block to the identity. ---*/
  SetVal2Diag(node_i, 1.0);

  /*--- Set known solution in rhs vector. ---*/
  b.SetBlock(node_i, x_i);
}

template <class ScalarType>
template <class OtherType>
void CSysMatrix<ScalarType>::EnforceZeroProjection(unsigned long node_i, const OtherType* n, CSysVector<OtherType>& b) {
  SU2_ZONE_SCOPED

  /*--- Visit all columns (L, diagonal, U) of row node_i. ---*/
  auto processCol = [&](unsigned long node_j, bool isDiag) {
    auto bij = GetBlock(node_i, node_j);
    auto bji = GetBlock(node_j, node_i);
    ScalarType nbn{};
    if (bji != nullptr) {
      for (auto iVar = 0ul; iVar < nVar; ++iVar) {
        ScalarType proj{};
        for (auto jVar = 0ul; jVar < nVar; ++jVar) proj += bji[iVar * nVar + jVar] * PassiveAssign(n[jVar]);
        for (auto jVar = 0ul; jVar < nVar; ++jVar) bji[iVar * nVar + jVar] -= proj * PassiveAssign(n[jVar]);
        nbn += proj * PassiveAssign(n[iVar]);
      }
    }
    for (auto jVar = 0ul; jVar < nVar; ++jVar) {
      ScalarType proj{};
      for (auto iVar = 0ul; iVar < nVar; ++iVar) proj += bij[iVar * nVar + jVar] * PassiveAssign(n[iVar]);
      for (auto iVar = 0ul; iVar < nVar; ++iVar) bij[iVar * nVar + jVar] -= proj * PassiveAssign(n[iVar]);
    }
    if (isDiag) {
      for (auto iVar = 0ul; iVar < nVar; ++iVar)
        for (auto jVar = 0ul; jVar < nVar; ++jVar)
          bij[iVar * nVar + jVar] += PassiveAssign(n[iVar]) * nbn * PassiveAssign(n[jVar]);
    }
  };
  for (auto k = mat.row_ptr_l[node_i]; k < mat.row_ptr_l[node_i + 1]; ++k) processCol(mat.col_ind_l[k], false);
  processCol(node_i, true);
  for (auto k = mat.row_ptr_u[node_i]; k < mat.row_ptr_u[node_i + 1]; ++k) processCol(mat.col_ind_u[k], false);

  OtherType proj{};
  for (auto iVar = 0ul; iVar < nVar; ++iVar) proj += b(node_i, iVar) * n[iVar];
  for (auto iVar = 0ul; iVar < nVar; ++iVar) b(node_i, iVar) -= proj * n[iVar];
}

template <class ScalarType>
void CSysMatrix<ScalarType>::SetDiagonalAsColumnSum() {
  SU2_ZONE_SCOPED
  const auto blkSz = nVar * nEqn;

  SU2_OMP_FOR_DYN(omp_heavy_size)
  for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
    auto* d_i = &mat.d[iPoint * blkSz];
    for (auto k = 0ul; k < blkSz; ++k) d_i[k] = 0.0;

    if (!quantized_mode) {
      /*--- For each L entry (iPoint, j): subtract its U-transpose (j, iPoint). ---*/
      for (auto k_l = mat.row_ptr_l[iPoint]; k_l < mat.row_ptr_l[iPoint + 1]; ++k_l)
        MatrixSubtraction(d_i, &mat.u[l_to_u_transp[k_l] * blkSz], d_i);

      /*--- For each U entry (iPoint, j): subtract its L-transpose (j, iPoint). ---*/
      for (auto k_u = mat.row_ptr_u[iPoint]; k_u < mat.row_ptr_u[iPoint + 1]; ++k_u)
        MatrixSubtraction(d_i, &mat.l[u_to_l_transp[k_u] * blkSz], d_i);
    } else {
      auto subtractTransp = [&](su2uint k_transp, const QuantType* qs, const QuantType* qv) {
        const CBlockView<const ScalarType> view{nullptr, &qs[k_transp * nVar], &qv[k_transp * blkSz], nVar};
        for (auto i = 0ul; i < nVar; ++i)
          for (auto j = 0ul; j < nEqn; ++j) d_i[i * nEqn + j] -= view(i, j);
      };
      for (auto k_l = mat.row_ptr_l[iPoint]; k_l < mat.row_ptr_l[iPoint + 1]; ++k_l)
        subtractTransp(l_to_u_transp[k_l], q_scale_u, q_blocks_u);
      for (auto k_u = mat.row_ptr_u[iPoint]; k_u < mat.row_ptr_u[iPoint + 1]; ++k_u)
        subtractTransp(u_to_l_transp[k_u], q_scale_l, q_blocks_l);
    }
  }
  END_SU2_OMP_FOR
}

template <class ScalarType>
void CSysMatrix<ScalarType>::TransposeInPlace() {
  SU2_ZONE_SCOPED
  assert(nVar == nEqn && "Cannot transpose with nVar != nEqn.");

  auto swapAndTransp = [](unsigned long n, ScalarType* a, ScalarType* b) {
    assert(a != b);
    /*--- a=b', b=a' ---*/
    for (auto i = 0ul; i < n; ++i) {
      for (auto j = 0ul; j < i; ++j) {
        const auto lo = i * n + j;
        const auto up = j * n + i;
        std::swap(a[lo], b[up]);
        std::swap(a[up], b[lo]);
      }
      std::swap(a[i * n + i], b[i * n + i]);
    }
  };

  /*--- Swap ij with ji and transpose them. ---*/

  if (edge_ptr_l) {
    /*--- FV path: each edge maps to one U and one L block. ---*/
    SU2_OMP_FOR_DYN(omp_heavy_size * 2)
    for (auto iEdge = 0ul; iEdge < mat.nnz_l; ++iEdge) {
      auto* bij_u = &mat.u[iEdge * nVar * nVar];
      auto* bji_l = &mat.l[edge_ptr_l[iEdge] * nVar * nVar];
      swapAndTransp(nVar, bij_u, bji_l);
    }
    END_SU2_OMP_FOR
  } else if (l_to_u_transp) {
    /*--- FEM/general path: use the L→U transpose map (one L entry per pair). ---*/
    SU2_OMP_FOR_DYN(omp_heavy_size)
    for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
      for (auto k_l = mat.row_ptr_l[iPoint]; k_l < mat.row_ptr_l[iPoint + 1]; ++k_l) {
        const auto k_u = l_to_u_transp[k_l];
        swapAndTransp(nVar, &mat.u[k_u * nVar * nVar], &mat.l[k_l * nVar * nVar]);
      }
    }
    END_SU2_OMP_FOR
  } else {
    /*--- Slow fallback: search for each U entry's L partner via GetBlock. ---*/
    SU2_OMP_FOR_DYN(omp_heavy_size)
    for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
      for (auto k_u = mat.row_ptr_u[iPoint]; k_u < mat.row_ptr_u[iPoint + 1]; ++k_u) {
        const auto jPoint = mat.col_ind_u[k_u];
        auto* bij = &mat.u[k_u * nVar * nVar];
        auto* bji = GetBlock(jPoint, iPoint);
        assert(bji && "Pattern is not symmetric.");
        swapAndTransp(nVar, bij, bji);
      }
    }
    END_SU2_OMP_FOR
  }

  /*--- Transpose the diagonal blocks. ---*/

  SU2_OMP_FOR_STAT(omp_heavy_size)
  for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
    auto bii = &mat.d[iPoint * nVar * nVar];
    for (auto i = 0ul; i < nVar; ++i)
      for (auto j = 0ul; j < i; ++j) std::swap(bii[i * nVar + j], bii[j * nVar + i]);
  }
  END_SU2_OMP_FOR

#ifdef HAVE_PASTIX
  SU2_OMP_MASTER
  pastix_wrapper.SetTransposedSolve();
  END_SU2_OMP_MASTER
#endif
}

template <class ScalarType>
void CSysMatrix<ScalarType>::MatrixMatrixAddition(ScalarType alpha, const CSysMatrix<ScalarType>& B) {
  SU2_ZONE_SCOPED
  /*--- Check that the LDU structure is shared (pointer equality since both come from CGeometry). ---*/
  const bool ok = (mat.row_ptr_l == B.mat.row_ptr_l) && (mat.col_ind_l == B.mat.col_ind_l) &&
                  (mat.row_ptr_u == B.mat.row_ptr_u) && (mat.col_ind_u == B.mat.col_ind_u) && (nVar == B.nVar) &&
                  (nEqn == B.nEqn) && (nPoint == B.nPoint) && (mat.nnz_l == B.mat.nnz_l) && (mat.nnz_u == B.mat.nnz_u);
  if (!ok) SU2_MPI::Error("Matrices do not have compatible sparsity.", CURRENT_FUNCTION);

  SU2_OMP_FOR_STAT(omp_light_size)
  for (auto i = 0ul; i < nPoint * nVar * nEqn; ++i) mat.d[i] += alpha * B.mat.d[i];
  END_SU2_OMP_FOR
  SU2_OMP_FOR_STAT(omp_light_size)
  for (auto i = 0ul; i < mat.nnz_l * nVar * nEqn; ++i) mat.l[i] += alpha * B.mat.l[i];
  END_SU2_OMP_FOR
  SU2_OMP_FOR_STAT(omp_light_size)
  for (auto i = 0ul; i < mat.nnz_u * nVar * nEqn; ++i) mat.u[i] += alpha * B.mat.u[i];
  END_SU2_OMP_FOR
}

template <class ScalarType>
void CSysMatrix<ScalarType>::BuildPastixPreconditioner(CGeometry* geometry, const CConfig* config,
                                                       unsigned short kind_fact) {
  SU2_ZONE_SCOPED
#ifdef HAVE_PASTIX
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    pastix_wrapper.SetLDU(nVar, nPoint, nPointDomain, mat.row_ptr_l, mat.col_ind_l, mat.row_ptr_u, mat.col_ind_u, mat.d,
                          mat.l, mat.u);
    pastix_wrapper.Factorize(geometry, config, kind_fact);
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
#else
  SU2_MPI::Error("SU2 was not compiled with -DHAVE_PASTIX", CURRENT_FUNCTION);
#endif
}

template <class ScalarType>
void CSysMatrix<ScalarType>::ComputePastixPreconditioner(const CSysVector<ScalarType>& vec,
                                                         CSysVector<ScalarType>& prod, CGeometry* geometry,
                                                         const CConfig* config) const {
  SU2_ZONE_SCOPED
#ifdef HAVE_PASTIX
  SU2_OMP_SAFE_GLOBAL_ACCESS(pastix_wrapper.Solve(vec, prod);)

  CSysMatrixComms::Initiate(prod, geometry, config);
  CSysMatrixComms::Complete(prod, geometry, config);
#else
  SU2_MPI::Error("SU2 was not compiled with -DHAVE_PASTIX", CURRENT_FUNCTION);
#endif
}

/*--- Explicit instantiations ---*/

#define INSTANTIATE_COMMS(TYPE)                                                                                       \
  template void CSysMatrixComms::Initiate<TYPE>(const CSysVector<TYPE>&, CGeometry*, const CConfig*, MPI_QUANTITIES); \
  template void CSysMatrixComms::Complete<TYPE>(CSysVector<TYPE>&, CGeometry*, const CConfig*, MPI_QUANTITIES);

#define INSTANTIATE_MATRIX(TYPE)                                                                                  \
  template class CSysMatrix<TYPE>;                                                                                \
  template void CSysMatrix<TYPE>::EnforceSolutionAtNode(unsigned long, const su2double*, CSysVector<su2double>&); \
  template void CSysMatrix<TYPE>::EnforceZeroProjection(unsigned long, const su2double*, CSysVector<su2double>&); \
  INSTANTIATE_COMMS(TYPE)

INSTANTIATE_MATRIX(su2mixedfloat)

#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
INSTANTIATE_MATRIX(passivedouble)
#endif
#ifdef CODI_REVERSE_TYPE
INSTANTIATE_COMMS(su2double)
#endif
