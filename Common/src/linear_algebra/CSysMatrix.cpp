/*!
 * \file CSysMatrix.cpp
 * \brief Implementation of the sparse matrix class.
 * \author F. Palacios, A. Bueno, T. Economon, P. Gomes
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

template <class ScalarType>
CSysMatrix<ScalarType>::CSysMatrix() : rank(SU2_MPI::GetRank()), size(SU2_MPI::GetSize()) {
  nPoint = nPointDomain = nVar = nEqn = 0;
  nnz = nnz_ilu = 0;
  ilu_fill_in = 0;

  omp_partitions = nullptr;

  matrix = nullptr;
  row_ptr = nullptr;
  dia_ptr = nullptr;
  col_ind = nullptr;
  col_ptr = nullptr;

  ILU_matrix = nullptr;
  row_ptr_ilu = nullptr;
  dia_ptr_ilu = nullptr;
  col_ind_ilu = nullptr;

  invM = nullptr;

#ifdef USE_MKL
  MatrixMatrixProductJitter = nullptr;
  MatrixVectorProductJitterBetaOne = nullptr;
  MatrixVectorProductJitterBetaZero = nullptr;
  MatrixVectorProductJitterAlphaMinusOne = nullptr;
#endif
}

template <class ScalarType>
CSysMatrix<ScalarType>::~CSysMatrix() {
  delete[] omp_partitions;
  MemoryAllocation::aligned_free(ILU_matrix);
  MemoryAllocation::aligned_free(matrix);
  MemoryAllocation::aligned_free(invM);

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
                                        const CConfig* config, bool needTranspPtr, bool grad_mode) {
  assert(omp_get_thread_num() == 0 && "Only the master thread is allowed to initialize the matrix.");

  if (npoint == 0) return;

  if (matrix != nullptr) {
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

  /*--- No else if, but separat if case! ---*/
  if (config->GetSmoothGradient() && grad_mode) {
    prec = config->GetKind_Grad_Linear_Solver_Prec();
  }

  const bool ilu_needed = (prec == ILU);
  const bool diag_needed = ilu_needed || (prec == JACOBI) || (prec == LINELET);

  /*--- Basic dimensions. ---*/
  nVar = nvar;
  nEqn = neqn;
  nPoint = npoint;
  nPointDomain = npointdomain;

  /*--- Get sparse structure pointers from geometry,
   *    the data is managed by CGeometry to allow re-use. ---*/

  const auto& csr = geometry->GetSparsePattern(type, 0);

  nnz = csr.getNumNonZeros();
  row_ptr = csr.outerPtr();
  col_ind = csr.innerIdx();
  dia_ptr = csr.diagPtr();

  if (needTranspPtr) col_ptr = geometry->GetTransposeSparsePatternMap(type).data();

  if (type == ConnectivityType::FiniteVolume) {
    edge_ptr.ptr = geometry->GetEdgeToSparsePatternMap().data();
    edge_ptr.nEdge = geometry->GetnEdge();
  }

  /*--- Get ILU sparse pattern, if fill is 0 no new data is allocated. --*/

  if (ilu_needed) {
    ilu_fill_in = config->GetLinear_Solver_ILU_n();

    const auto& csr_ilu = geometry->GetSparsePattern(type, ilu_fill_in);

    row_ptr_ilu = csr_ilu.outerPtr();
    col_ind_ilu = csr_ilu.innerIdx();
    dia_ptr_ilu = csr_ilu.diagPtr();
    nnz_ilu = csr_ilu.getNumNonZeros();
  }

  /*--- Allocate data. ---*/
  auto allocAndInit = [](ScalarType*& ptr, unsigned long num) {
    ptr = MemoryAllocation::aligned_alloc<ScalarType, true>(64, num * sizeof(ScalarType));
  };

  allocAndInit(matrix, nnz * nVar * nEqn);

  /*--- Preconditioners. ---*/

  if (ilu_needed) allocAndInit(ILU_matrix, nnz_ilu * nVar * nEqn);

  if (diag_needed) allocAndInit(invM, nPointDomain * nVar * nEqn);

  /*--- Thread parallel initialization. ---*/

  int num_threads = omp_get_max_threads();

  /*--- Set suitable chunk sizes for light static for loops, and heavy
   dynamic ones, such that threads are approximately evenly loaded. ---*/
  omp_light_size = computeStaticChunkSize(nnz * nVar * nEqn, num_threads, OMP_MAX_SIZE_L);
  omp_heavy_size = computeStaticChunkSize(nPointDomain, num_threads, OMP_MAX_SIZE_H);

  omp_num_parts = config->GetLinear_Solver_Prec_Threads();
  if (omp_num_parts == 0) omp_num_parts = num_threads;

  /*--- This is akin to the row_ptr. ---*/
  omp_partitions = new unsigned long[omp_num_parts + 1];
  for (unsigned long i = 0; i <= omp_num_parts; ++i) omp_partitions[i] = nPointDomain;

  /*--- Work estimate based on non-zeros to produce balanced partitions. ---*/

  const auto row_ptr_prec = ilu_needed ? row_ptr_ilu : row_ptr;
  const auto nnz_prec = row_ptr_prec[nPointDomain];

  const auto nnz_per_part = roundUpDiv(nnz_prec, omp_num_parts);

  for (auto iPoint = 0ul, part = 0ul; iPoint < nPointDomain; ++iPoint) {
    if (row_ptr_prec[iPoint] >= part * nnz_per_part) omp_partitions[part++] = iPoint;
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
                               unsigned short commType) {
  if (geometry->nP2PSend == 0) return;

  /*--- Local variables ---*/

  const unsigned short COUNT_PER_POINT = x.GetNVar();
  const unsigned short MPI_TYPE = COMM_TYPE_DOUBLE;

  /*--- Create a boolean for reversing the order of comms. ---*/

  const bool reverse = (commType == SOLUTION_MATRIXTRANS);

  /*--- Set the size of the data packet and type depending on quantity. ---*/

  switch (commType) {
    case SOLUTION_MATRIX:
    case SOLUTION_MATRIXTRANS:
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
      case SOLUTION_MATRIX: {
        su2double* bufDSend = geometry->bufD_P2PSend;

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

      case SOLUTION_MATRIXTRANS: {
        /*--- We are going to communicate in reverse, so we use the
         recv buffer for the send instead. Also, all of the offsets
         and counts are derived from the recv data structures. ---*/

        su2double* bufDSend = geometry->bufD_P2PRecv;

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
void CSysMatrixComms::Complete(CSysVector<T>& x, CGeometry* geometry, const CConfig* config, unsigned short commType) {
  if (geometry->nP2PRecv == 0) return;

  /*--- Local variables ---*/

  const unsigned short COUNT_PER_POINT = x.GetNVar();

  /*--- Global status so all threads can see the result of Waitany. ---*/
  static SU2_MPI::Status status;
  int ind;

  /*--- Store the data that was communicated into the appropriate
   location within the local class data structures. ---*/

  for (auto iMessage = 0; iMessage < geometry->nP2PRecv; iMessage++) {
    /*--- For efficiency, recv the messages dynamically based on
     the order they arrive. ---*/

    SU2_OMP_SAFE_GLOBAL_ACCESS(SU2_MPI::Waitany(geometry->nP2PRecv, geometry->req_P2PRecv, &ind, &status);)

    /*--- Once we have recv'd a message, get the source rank. ---*/

    const auto source = status.MPI_SOURCE;

    switch (commType) {
      case SOLUTION_MATRIX: {
        const su2double* bufDRecv = geometry->bufD_P2PRecv;

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

      case SOLUTION_MATRIXTRANS: {
        /*--- We are going to communicate in reverse, so we use the
         send buffer for the recv instead. Also, all of the offsets
         and counts are derived from the send data structures. ---*/

        const su2double* bufDRecv = geometry->bufD_P2PSend;

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
  SU2_OMP_SAFE_GLOBAL_ACCESS(SU2_MPI::Waitall(geometry->nP2PSend, geometry->req_P2PSend, MPI_STATUS_IGNORE);)
#endif
}

template <class ScalarType>
void CSysMatrix<ScalarType>::SetValZero() {
  const auto size = nnz * nVar * nEqn;
  const auto chunk = roundUpDiv(size, omp_get_num_threads());
  const auto begin = chunk * omp_get_thread_num();
  const auto mySize = min(chunk, size - begin) * sizeof(ScalarType);
  memset(&matrix[begin], 0, mySize);
  SU2_OMP_BARRIER
}

template <class ScalarType>
void CSysMatrix<ScalarType>::SetValDiagonalZero() {
  SU2_OMP_FOR_STAT(omp_heavy_size)
  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint)
    for (auto index = 0ul; index < nVar * nEqn; ++index) matrix[dia_ptr[iPoint] * nVar * nEqn + index] = 0.0;
  END_SU2_OMP_FOR
}

template <class ScalarType>
void CSysMatrix<ScalarType>::Gauss_Elimination(ScalarType* matrix, ScalarType* vec) const {
#ifdef USE_MKL_LAPACK
  // With MKL_DIRECT_CALL enabled, this is significantly faster than native code on Intel Architectures.
  lapack_int ipiv[MAXNVAR];
  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, nVar, nVar, matrix, nVar, ipiv);
  LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', nVar, 1, matrix, nVar, ipiv, vec, 1);
#else
#define A(I, J) matrix[(I)*nVar + (J)]

  /*--- Transform system in Upper Matrix ---*/
  for (auto iVar = 1ul; iVar < nVar; iVar++) {
    for (auto jVar = 0ul; jVar < iVar; jVar++) {
      ScalarType weight = A(iVar, jVar) / A(jVar, jVar);
      for (auto kVar = jVar; kVar < nVar; kVar++) A(iVar, kVar) -= weight * A(jVar, kVar);
      vec[iVar] -= weight * vec[jVar];
    }
  }

  /*--- Backwards substitution ---*/
  for (auto iVar = nVar; iVar > 0ul;) {
    iVar--;  // unsigned type
    for (auto jVar = iVar + 1; jVar < nVar; jVar++) vec[iVar] -= A(iVar, jVar) * vec[jVar];
    vec[iVar] /= A(iVar, iVar);
  }
#undef A
#endif
}

template <class ScalarType>
void CSysMatrix<ScalarType>::MatrixInverse(ScalarType* matrix, ScalarType* inverse) const {
  /*--- This is a generalization of Gaussian elimination for multiple rhs' (the basis vectors).
   We could call "Gauss_Elimination" multiple times or fully generalize it for multiple rhs,
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
  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, nVar, nVar, matrix, nVar, ipiv);
  LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', nVar, nVar, matrix, nVar, ipiv, inverse, nVar);
#else
#define A(I, J) matrix[(I)*nVar + (J)]

  /*--- Transform system in Upper Matrix ---*/
  for (auto iVar = 1ul; iVar < nVar; iVar++) {
    for (auto jVar = 0ul; jVar < iVar; jVar++) {
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

    for (auto kVar = 0ul; kVar < nVar; kVar++) M(iVar, kVar) /= A(iVar, iVar);
  }
#undef A
#endif
#undef M
}

template <class ScalarType>
void CSysMatrix<ScalarType>::DeleteValsRowi(unsigned long i) {
  const auto block_i = i / nVar;
  const auto row = i % nVar;

  for (auto index = row_ptr[block_i]; index < row_ptr[block_i + 1]; index++) {
    for (auto iVar = 0u; iVar < nVar; iVar++)
      matrix[index * nVar * nVar + row * nVar + iVar] = 0.0;  // Delete row values in the block
    if (col_ind[index] == block_i)
      matrix[index * nVar * nVar + row * nVar + row] = 1.0;  // Set 1 to the diagonal element
  }
}

template <class ScalarType>
void CSysMatrix<ScalarType>::MatrixVectorProduct(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                                 CGeometry* geometry, const CConfig* config) const {
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

  SU2_OMP_FOR_DYN(omp_heavy_size)
  for (auto row_i = 0ul; row_i < nPointDomain; row_i++) {
    RowProduct(vec, row_i, &prod[row_i * nVar]);
  }
  END_SU2_OMP_FOR

  /*--- MPI Parallelization. ---*/

  CSysMatrixComms::Initiate(prod, geometry, config);
  CSysMatrixComms::Complete(prod, geometry, config);
}

template <class ScalarType>
void CSysMatrix<ScalarType>::BuildJacobiPreconditioner() {
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
  /*--- Copy block matrix to compute factorization in-place. ---*/

  if (ilu_fill_in == 0) {
    /*--- ILU0, direct copy. ---*/
    SU2_OMP_FOR_STAT(omp_light_size)
    for (auto iVar = 0ul; iVar < nnz * nVar * nVar; ++iVar) ILU_matrix[iVar] = matrix[iVar];
    END_SU2_OMP_FOR
  } else {
    /*--- ILUn clear the ILU matrix first. ---*/
    SU2_OMP_FOR_STAT(omp_light_size)
    for (auto iVar = 0ul; iVar < nnz_ilu * nVar * nVar; iVar++) ILU_matrix[iVar] = 0.0;
    END_SU2_OMP_FOR

    /*--- ILUn, traverse matrix to access its blocks
     *    sequentially and set them in the ILU matrix. ---*/
    SU2_OMP_FOR_DYN(omp_heavy_size)
    for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {
      for (auto index = row_ptr[iPoint]; index < row_ptr[iPoint + 1]; index++) {
        auto jPoint = col_ind[index];
        SetBlock_ILUMatrix(iPoint, jPoint, &matrix[index * nVar * nVar]);
      }
    }
    END_SU2_OMP_FOR
  }

  /*--- Transform system in Upper Matrix ---*/

  /*--- OpenMP Parallelization, a loop construct is used to ensure
   *    the preconditioner is computed correctly even if called
   *    outside of a parallel section. ---*/

  SU2_OMP_FOR_STAT(1)
  for (unsigned long thread = 0; thread < omp_num_parts; ++thread) {
    const auto begin = omp_partitions[thread];
    const auto end = omp_partitions[thread + 1];
    if (begin == end) continue;

    /*--- Each thread will work on the submatrix defined from row/col "begin"
     *    to row/col "end-1" (i.e. the range [begin,end[). Which is exactly
     *    what the MPI-only implementation does. ---*/

    ScalarType weight[MAXNVAR * MAXNVAR], aux_block[MAXNVAR * MAXNVAR];

    for (auto iPoint = begin + 1; iPoint < end; iPoint++) {
      /*--- Invert and store the previous diagonal block to later compute the weight. ---*/

      InverseDiagonalBlock_ILUMatrix(iPoint - 1, &invM[(iPoint - 1) * nVar * nVar]);

      /*--- For this row (unknown), loop over its lower diagonal entries. ---*/

      for (auto index = row_ptr_ilu[iPoint]; index < dia_ptr_ilu[iPoint]; index++) {
        /*--- jPoint is the column index (jPoint < iPoint). ---*/

        auto jPoint = col_ind_ilu[index];

        /*--- We only care about the sub matrix within "begin" and "end-1". ---*/

        if (jPoint < begin) continue;

        /*--- Multiply the block by the inverse of the corresponding diagonal block. ---*/

        auto Block_ij = &ILU_matrix[index * nVar * nVar];
        MatrixMatrixProduct(Block_ij, &invM[jPoint * nVar * nVar], weight);

        /*--- "weight" holds Aij*inv(Ajj). Jump to the upper part of the jPoint row. ---*/

        for (auto index_ = dia_ptr_ilu[jPoint] + 1; index_ < row_ptr_ilu[jPoint + 1]; index_++) {
          /*--- Get the column index (kPoint > jPoint). ---*/

          auto kPoint = col_ind_ilu[index_];

          if (kPoint >= end) break;

          /*--- If Aik exists, update it: Aik -= Aij*inv(Ajj)*Ajk ---*/

          auto Block_ik = GetBlock_ILUMatrix(iPoint, kPoint);

          if (Block_ik != nullptr) {
            auto Block_jk = &ILU_matrix[index_ * nVar * nVar];
            MatrixMatrixProduct(weight, Block_jk, aux_block);
            MatrixSubtraction(Block_ik, aux_block, Block_ik);
          }
        }

        /*--- Lastly, store "weight" in the lower triangular part, which
         will be reused during the forward solve in the precon/smoother. ---*/

        for (auto iVar = 0ul; iVar < nVar * nVar; ++iVar) Block_ij[iVar] = weight[iVar];
      }
    }
    InverseDiagonalBlock_ILUMatrix(end - 1, &invM[(end - 1) * nVar * nVar]);
  }
  END_SU2_OMP_FOR
}

template <class ScalarType>
void CSysMatrix<ScalarType>::ComputeILUPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                                      CGeometry* geometry, const CConfig* config) const {
  /*--- Coherent view of vectors. ---*/
  SU2_OMP_BARRIER

  /*--- OpenMP Parallelization ---*/
  SU2_OMP_FOR_STAT(1)
  for (unsigned long thread = 0; thread < omp_num_parts; ++thread) {
    const auto begin = omp_partitions[thread];
    const auto end = omp_partitions[thread + 1];
    if (begin == end) continue;

    ScalarType aux_vec[MAXNVAR];

    /*--- Copy vector to then work on prod in place ---*/

    for (auto iVar = begin * nVar; iVar < end * nVar; iVar++) prod[iVar] = vec[iVar];

    /*--- Forward solve the system using the lower matrix entries that
     were computed and stored during the ILU preprocessing. Note
     that we are overwriting the residual vector as we go. ---*/

    for (auto iPoint = begin + 1; iPoint < end; iPoint++) {
      for (auto index = row_ptr_ilu[iPoint]; index < dia_ptr_ilu[iPoint]; index++) {
        auto jPoint = col_ind_ilu[index];
        if (jPoint < begin) continue;
        auto Block_ij = &ILU_matrix[index * nVar * nVar];
        MatrixVectorProductSub(Block_ij, &prod[jPoint * nVar], &prod[iPoint * nVar]);
      }
    }

    /*--- Backwards substitution (starts at the last row) ---*/

    for (auto iPoint = end; iPoint > begin;) {
      iPoint--;  // unsigned type
      for (auto iVar = 0ul; iVar < nVar; iVar++) aux_vec[iVar] = prod[iPoint * nVar + iVar];

      for (auto index = dia_ptr_ilu[iPoint] + 1; index < row_ptr_ilu[iPoint + 1]; index++) {
        auto jPoint = col_ind_ilu[index];
        if (jPoint >= end) break;
        auto Block_ij = &ILU_matrix[index * nVar * nVar];
        MatrixVectorProductSub(Block_ij, &prod[jPoint * nVar], aux_vec);
      }

      MatrixVectorProduct(&invM[iPoint * nVar * nVar], aux_vec, &prod[iPoint * nVar]);
    }
  }
  END_SU2_OMP_FOR

  /*--- MPI Parallelization ---*/

  CSysMatrixComms::Initiate(prod, geometry, config);
  CSysMatrixComms::Complete(prod, geometry, config);
}

template <class ScalarType>
void CSysMatrix<ScalarType>::ComputeLU_SGSPreconditioner(const CSysVector<ScalarType>& vec,
                                                         CSysVector<ScalarType>& prod, CGeometry* geometry,
                                                         const CConfig* config) const {
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

    for (auto iPoint = begin; iPoint < end; ++iPoint) {
      auto idx = iPoint * nVar;
      LowerProduct(prod, iPoint, begin, low_prod);         // Compute L.x*
      VectorSubtraction(&vec[idx], low_prod, &prod[idx]);  // Compute y = b - L.x*
      Gauss_Elimination(iPoint, &prod[idx]);               // Solve D.x* = y
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

    for (auto iPoint = row_end; iPoint > begin;) {
      iPoint--;  // because of unsigned type
      auto idx = iPoint * nVar;
      DiagonalProduct(prod, iPoint, dia_prod);           // Compute D.x*
      UpperProduct(prod, iPoint, row_end, up_prod);      // Compute U.x_(n+1)
      VectorSubtraction(dia_prod, up_prod, &prod[idx]);  // Compute y = D.x*-U.x_(n+1)
      Gauss_Elimination(iPoint, &prod[idx]);             // Solve D.x* = y
    }
  }
  END_SU2_OMP_FOR

  /*--- MPI Parallelization ---*/

  CSysMatrixComms::Initiate(prod, geometry, config);
  CSysMatrixComms::Complete(prod, geometry, config);
}

template <class ScalarType>
void CSysMatrix<ScalarType>::BuildLineletPreconditioner(const CGeometry* geometry, const CConfig* config) {
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
    MatrixCopy(&matrix[dia_ptr[li.linelets[iLinelet][0]] * nVar * nVar], lineletInvDiag.data());

    for (auto iElem = 1ul; iElem < nElem; iElem++) {
      /*--- Setup pointers to required matrices and vectors ---*/
      const auto im1Point = li.linelets[iLinelet][iElem - 1];
      const auto iPoint = li.linelets[iLinelet][iElem];

      const auto* d = &matrix[dia_ptr[iPoint] * nVar * nVar];
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
    Gauss_Elimination(&lineletInvDiag[(nElem - 1) * nVar * nVar], &lineletVector[(nElem - 1) * nVar]);

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
  /*--- Eliminate the row associated with node i (Block_ii = I and all other Block_ij = 0).
   *    To preserve eventual symmetry, also attempt to eliminate the column, if the sparse pattern is not
   *    symmetric the entire column may not be eliminated, the result (matrix and vector) is still correct.
   *    The vector is updated with the product of column i by the known (enforced) solution at node i. ---*/

  for (auto index = row_ptr[node_i]; index < row_ptr[node_i + 1]; ++index) {
    auto node_j = col_ind[index];

    /*--- The diagonal block is handled outside the loop. ---*/
    if (node_j == node_i) continue;

    /*--- Delete block j on row i (bij) and ATTEMPT to delete block i on row j (bji). ---*/
    auto bij = &matrix[index * nVar * nVar];
    auto bji = GetBlock(node_j, node_i);

    /*--- The "attempt" part. ---*/
    if (bji == nullptr) {
      node_j = node_i;
      bji = bij;
    }

    for (auto iVar = 0ul; iVar < nVar; ++iVar) {
      for (auto jVar = 0ul; jVar < nVar; ++jVar) {
        /*--- Column product. ---*/
        b[node_j * nVar + iVar] -= bji[iVar * nVar + jVar] * x_i[jVar];
        /*--- Delete blocks. ---*/
        bij[iVar * nVar + jVar] = bji[iVar * nVar + jVar] = 0.0;
      }
    }
  }

  /*--- Set the diagonal block to the identity. ---*/
  SetVal2Diag(node_i, 1.0);

  /*--- Set known solution in rhs vector. ---*/
  b.SetBlock(node_i, x_i);
}

template <class ScalarType>
template <class OtherType>
void CSysMatrix<ScalarType>::EnforceSolutionAtDOF(unsigned long node_i, unsigned long iVar, OtherType x_i,
                                                  CSysVector<OtherType>& b) {
  for (auto index = row_ptr[node_i]; index < row_ptr[node_i + 1]; ++index) {
    const auto node_j = col_ind[index];

    /*--- Delete row iVar of block j on row i (bij) and ATTEMPT
     *    to delete column iVar block i on row j (bji). ---*/

    auto bij = &matrix[index * nVar * nVar];
    auto bji = GetBlock(node_j, node_i);

    /*--- The "attempt" part. ---*/
    if (bji != nullptr) {
      for (auto jVar = 0ul; jVar < nVar; ++jVar) {
        /*--- Column product. ---*/
        b[node_j * nVar + jVar] -= bji[jVar * nVar + iVar] * x_i;
        /*--- Delete entries. ---*/
        bji[jVar * nVar + iVar] = 0.0;
      }
    }

    /*--- Delete row. ---*/
    for (auto jVar = 0ul; jVar < nVar; ++jVar) bij[iVar * nVar + jVar] = 0.0;

    /*--- Set the diagonal entry of the block to 1. ---*/
    if (node_j == node_i) bij[iVar * (nVar + 1)] = 1.0;
  }

  /*--- Set known solution in rhs vector. ---*/
  b(node_i, iVar) = x_i;
}

template <class ScalarType>
void CSysMatrix<ScalarType>::SetDiagonalAsColumnSum() {
  SU2_OMP_FOR_DYN(omp_heavy_size)
  for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
    auto block_ii = &matrix[dia_ptr[iPoint] * nVar * nEqn];

    for (auto k = 0ul; k < nVar * nEqn; ++k) block_ii[k] = 0.0;

    for (auto k = row_ptr[iPoint]; k < row_ptr[iPoint + 1]; ++k) {
      auto block_ji = &matrix[col_ptr[k] * nVar * nEqn];
      if (block_ji != block_ii) MatrixSubtraction(block_ii, block_ji, block_ii);
    }
  }
  END_SU2_OMP_FOR
}

template <class ScalarType>
void CSysMatrix<ScalarType>::TransposeInPlace() {
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

  if (edge_ptr) {
    /*--- The FV way. ---*/
    SU2_OMP_FOR_DYN(omp_heavy_size * 2)
    for (auto iEdge = 0ul; iEdge < edge_ptr.nEdge; ++iEdge) {
      auto bij = &matrix[edge_ptr(iEdge, 0) * nVar * nVar];
      auto bji = &matrix[edge_ptr(iEdge, 1) * nVar * nVar];

      swapAndTransp(nVar, bij, bji);
    }
    END_SU2_OMP_FOR
  } else if (col_ptr) {
    /*--- If the column pointer was built. ---*/
    SU2_OMP_FOR_DYN(omp_heavy_size)
    for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
      for (auto k = row_ptr[iPoint]; k < dia_ptr[iPoint]; ++k) {
        auto bij = &matrix[k * nVar * nVar];
        auto bji = &matrix[col_ptr[k] * nVar * nVar];

        swapAndTransp(nVar, bij, bji);
      }
    }
    END_SU2_OMP_FOR
  } else {
    /*--- Slow fallback, needs to search for ji. ---*/
    SU2_OMP_FOR_DYN(omp_heavy_size)
    for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
      for (auto k = dia_ptr[iPoint] + 1ul; k < row_ptr[iPoint + 1]; ++k) {
        const auto jPoint = col_ind[k];
        auto bij = &matrix[k * nVar * nVar];
        auto bji = GetBlock(jPoint, iPoint);
        assert(bji && "Pattern is not symmetric.");

        swapAndTransp(nVar, bij, bji);
      }
    }
    END_SU2_OMP_FOR
  }

  /*--- Transpose the diagonal blocks. ---*/

  SU2_OMP_FOR_STAT(omp_heavy_size)
  for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
    auto bii = &matrix[dia_ptr[iPoint] * nVar * nVar];
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
  /*--- Check that the sparse structure is shared between the two matrices,
   *    comparing pointers is ok as they are obtained from CGeometry. ---*/
  bool ok = (row_ptr == B.row_ptr) && (col_ind == B.col_ind) && (nVar == B.nVar) && (nEqn == B.nEqn) && (nnz == B.nnz);

  if (!ok) {
    SU2_MPI::Error("Matrices do not have compatible sparsity.", CURRENT_FUNCTION);
  }

  SU2_OMP_FOR_STAT(omp_light_size)
  for (auto i = 0ul; i < nnz * nVar * nEqn; ++i) matrix[i] += alpha * B.matrix[i];
  END_SU2_OMP_FOR
}

template <class ScalarType>
void CSysMatrix<ScalarType>::BuildPastixPreconditioner(CGeometry* geometry, const CConfig* config,
                                                       unsigned short kind_fact) {
#ifdef HAVE_PASTIX
  /*--- Pastix will launch nested threads. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    pastix_wrapper.SetMatrix(nVar, nPoint, nPointDomain, row_ptr, col_ind, matrix);
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
  template void CSysMatrixComms::Initiate<TYPE>(const CSysVector<TYPE>&, CGeometry*, const CConfig*, unsigned short); \
  template void CSysMatrixComms::Complete<TYPE>(CSysVector<TYPE>&, CGeometry*, const CConfig*, unsigned short);

#define INSTANTIATE_MATRIX(TYPE)                                                                                  \
  template class CSysMatrix<TYPE>;                                                                                \
  template void CSysMatrix<TYPE>::EnforceSolutionAtNode(unsigned long, const su2double*, CSysVector<su2double>&); \
  template void CSysMatrix<TYPE>::EnforceSolutionAtDOF(unsigned long, unsigned long, su2double,                   \
                                                       CSysVector<su2double>&);                                   \
  INSTANTIATE_COMMS(TYPE)

#ifdef CODI_FORWARD_TYPE
/*--- In forward AD only the active type is used. ---*/
INSTANTIATE_MATRIX(su2double)
#else
/*--- Base and reverse AD, matrix is passive. ---*/
INSTANTIATE_MATRIX(su2mixedfloat)
/*--- If using mixed precision (float) instantiate also a version for doubles, and allow cross communications. ---*/
#ifdef USE_MIXED_PRECISION
INSTANTIATE_MATRIX(passivedouble)
#endif
#ifdef CODI_REVERSE_TYPE
INSTANTIATE_COMMS(su2double)
#endif
#endif  // CODI_FORWARD_TYPE
