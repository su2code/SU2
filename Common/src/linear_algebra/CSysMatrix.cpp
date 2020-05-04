/*!
 * \file CSysMatrix.cpp
 * \brief Implementation of the sparse matrix class.
 * \author F. Palacios, A. Bueno, T. Economon
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../include/CConfig.hpp"
#include "../../include/omp_structure.hpp"
#include "../../include/toolboxes/allocation_toolbox.hpp"

#include <cmath>

template<class ScalarType>
CSysMatrix<ScalarType>::CSysMatrix(void) {

  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();

  nPoint = nPointDomain = nVar = nEqn = 0;
  nnz = nnz_ilu = 0;
  ilu_fill_in = 0;
  nLinelet = 0;

  omp_partitions    = nullptr;

  matrix            = nullptr;
  row_ptr           = nullptr;
  dia_ptr           = nullptr;
  col_ind           = nullptr;

  ILU_matrix        = nullptr;
  row_ptr_ilu       = nullptr;
  dia_ptr_ilu       = nullptr;
  col_ind_ilu       = nullptr;

  invM              = nullptr;

#ifdef USE_MKL
  MatrixMatrixProductJitter              = nullptr;
  MatrixVectorProductJitterBetaOne       = nullptr;
  MatrixVectorProductJitterBetaZero      = nullptr;
  MatrixVectorProductJitterAlphaMinusOne = nullptr;
  MatrixVectorProductTranspJitterBetaOne = nullptr;
#endif

}

template<class ScalarType>
CSysMatrix<ScalarType>::~CSysMatrix(void) {

  if (omp_partitions != nullptr) delete [] omp_partitions;
  if (ILU_matrix != nullptr) MemoryAllocation::aligned_free(ILU_matrix);
  if (matrix != nullptr) MemoryAllocation::aligned_free(matrix);
  if (invM != nullptr) MemoryAllocation::aligned_free(invM);

#ifdef USE_MKL
  if ( MatrixMatrixProductJitter != nullptr )              mkl_jit_destroy( MatrixMatrixProductJitter );
  if ( MatrixVectorProductJitterBetaZero != nullptr )      mkl_jit_destroy( MatrixVectorProductJitterBetaZero );
  if ( MatrixVectorProductJitterBetaOne  != nullptr )      mkl_jit_destroy( MatrixVectorProductJitterBetaOne );
  if ( MatrixVectorProductJitterAlphaMinusOne != nullptr ) mkl_jit_destroy( MatrixVectorProductJitterAlphaMinusOne );
  if ( MatrixVectorProductTranspJitterBetaOne != nullptr ) mkl_jit_destroy( MatrixVectorProductTranspJitterBetaOne );
#endif

}

template<class ScalarType>
void CSysMatrix<ScalarType>::Initialize(unsigned long npoint, unsigned long npointdomain,
                                        unsigned short nvar, unsigned short neqn,
                                        bool EdgeConnect, CGeometry *geometry,
                                        CConfig *config, bool needTranspPtr) {

  assert(omp_get_thread_num()==0 && "Only the master thread is allowed to initialize the matrix.");

  if(matrix != nullptr) {
    SU2_OMP_MASTER
    SU2_MPI::Error("CSysMatrix can only be initialized once.", CURRENT_FUNCTION);
  }

  if(nvar > MAXNVAR) {
    SU2_OMP_MASTER
    SU2_MPI::Error("nVar larger than expected, increase MAXNVAR.", CURRENT_FUNCTION);
  }

  /*--- Application of this matrix, FVM or FEM. ---*/
  auto type = EdgeConnect? ConnectivityType::FiniteVolume : ConnectivityType::FiniteElement;

  /*--- Types of preconditioner the matrix will be asked to build. ---*/
  unsigned short sol_prec = config->GetKind_Linear_Solver_Prec();
  unsigned short def_prec = config->GetKind_Deform_Linear_Solver_Prec();
  unsigned short adj_prec = config->GetKind_DiscAdj_Linear_Prec();
  bool adjoint = config->GetDiscrete_Adjoint();

  bool ilu_needed = (sol_prec==ILU) || (def_prec==ILU) || (adjoint && (adj_prec==ILU));

  /*--- Basic dimensions. ---*/
  nVar = nvar;
  nEqn = neqn;
  nPoint = npoint;
  nPointDomain = npointdomain;

  /*--- Get sparse structure pointers from geometry,
   *    the data is managed by CGeometry to allow re-use. ---*/

  const auto& csr = geometry->GetSparsePattern(type,0);

  nnz = csr.getNumNonZeros();
  row_ptr = csr.outerPtr();
  col_ind = csr.innerIdx();
  dia_ptr = csr.diagPtr();

  if (needTranspPtr)
    col_ptr = geometry->GetTransposeSparsePatternMap(type).data();

  if (type == ConnectivityType::FiniteVolume)
    edge_ptr.ptr = geometry->GetEdgeToSparsePatternMap().data();

  /*--- Get ILU sparse pattern, if fill is 0 no new data is allocated. --*/

  if(ilu_needed)
  {
    ilu_fill_in = config->GetLinear_Solver_ILU_n();

    const auto& csr_ilu = geometry->GetSparsePattern(type, ilu_fill_in);

    row_ptr_ilu = csr_ilu.outerPtr();
    col_ind_ilu = csr_ilu.innerIdx();
    dia_ptr_ilu = csr_ilu.diagPtr();
    nnz_ilu = csr_ilu.getNumNonZeros();
  }

  /*--- Allocate data. ---*/
#define ALLOC_AND_INIT(ptr,num) {\
  ptr = MemoryAllocation::aligned_alloc<ScalarType>(64,num*sizeof(ScalarType));\
  for(size_t k=0; k<num; ++k) ptr[k]=0.0; }

  ALLOC_AND_INIT(matrix, nnz*nVar*nEqn)

  /*--- Preconditioners. ---*/

  if (ilu_needed) {
    ALLOC_AND_INIT(ILU_matrix, nnz_ilu*nVar*nEqn)
  }

  if (ilu_needed || (sol_prec==JACOBI) || (sol_prec==LINELET) ||
      (adjoint && (adj_prec==JACOBI)) || (def_prec==JACOBI))
  {
    ALLOC_AND_INIT(invM, nPointDomain*nVar*nEqn);
  }
#undef ALLOC_AND_INIT

  /*--- Thread parallel initialization. ---*/

  int num_threads = omp_get_max_threads();

  /*--- Set suitable chunk sizes for light static for loops, and heavy
   dynamic ones, such that threads are approximately evenly loaded. ---*/
  omp_light_size = computeStaticChunkSize(nnz*nVar*nEqn, num_threads, OMP_MAX_SIZE_L);
  omp_heavy_size = computeStaticChunkSize(nPointDomain, num_threads, OMP_MAX_SIZE_H);

  omp_num_parts = config->GetLinear_Solver_Prec_Threads();
  if (omp_num_parts == 0) omp_num_parts = num_threads;

  /*--- This is akin to the row_ptr. ---*/
  omp_partitions = new unsigned long [omp_num_parts+1];

  /// TODO: Use a work estimate to produce more balanced partitions.
  auto pts_per_part = roundUpDiv(nPointDomain, omp_num_parts);
  for(auto part = 0ul; part < omp_num_parts; ++part)
    omp_partitions[part] = part * pts_per_part;
  omp_partitions[omp_num_parts] = nPointDomain;

  /*--- Generate MKL Kernels ---*/

#ifdef USE_MKL
  mkl_jit_create_dgemm( &MatrixMatrixProductJitter, MKL_ROW_MAJOR, MKL_NOTRANS, MKL_NOTRANS, nVar, nVar, nVar,  1.0, nVar, nVar, 0.0, nVar );
  MatrixMatrixProductKernel = mkl_jit_get_dgemm_ptr( MatrixMatrixProductJitter );

  mkl_jit_create_dgemm( &MatrixVectorProductJitterBetaZero, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, 1, nVar, nVar,  1.0, 1, nVar, 0.0, 1 );
  MatrixVectorProductKernelBetaZero = mkl_jit_get_dgemm_ptr( MatrixVectorProductJitterBetaZero );

  mkl_jit_create_dgemm( &MatrixVectorProductJitterBetaOne, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, 1, nVar, nVar,  1.0, 1, nVar, 1.0, 1 );
  MatrixVectorProductKernelBetaOne = mkl_jit_get_dgemm_ptr( MatrixVectorProductJitterBetaOne );

  mkl_jit_create_dgemm( &MatrixVectorProductJitterAlphaMinusOne, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, 1, nVar, nVar, -1.0, 1, nVar, 1.0, 1 );
  MatrixVectorProductKernelAlphaMinusOne = mkl_jit_get_dgemm_ptr( MatrixVectorProductJitterAlphaMinusOne );

  mkl_jit_create_dgemm( &MatrixVectorProductTranspJitterBetaOne, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, nVar, 1, nVar,  1.0, nVar, nVar, 1.0, nVar );
  MatrixVectorProductTranspKernelBetaOne = mkl_jit_get_dgemm_ptr( MatrixVectorProductTranspJitterBetaOne );
#endif

}

template<class ScalarType>
template<class OtherType>
void CSysMatrix<ScalarType>::InitiateComms(const CSysVector<OtherType> & x,
                                           CGeometry *geometry,
                                           CConfig *config,
                                           unsigned short commType) const {

  /*--- Local variables ---*/

  unsigned short iVar;
  unsigned short COUNT_PER_POINT = 0;
  unsigned short MPI_TYPE        = 0;

  unsigned long iPoint, msg_offset, buf_offset;

  int iMessage, iSend, nSend;

  /*--- Create a boolean for reversing the order of comms. ---*/

  bool reverse = false;

  /*--- Set the size of the data packet and type depending on quantity. ---*/

  switch (commType) {
    case SOLUTION_MATRIX:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      reverse          = false;
      break;
    case SOLUTION_MATRIXTRANS:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      reverse          = true;
      break;
    default:
      SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.",
                     CURRENT_FUNCTION);
      break;
  }

  /*--- Check to make sure we have created a large enough buffer
   for these comms during preprocessing. This is only for the su2double
   buffer. It will be reallocated whenever we find a larger count
   per point. After the first cycle of comms, this should be inactive. ---*/

  if (COUNT_PER_POINT > geometry->countPerPoint) {
    geometry->AllocateP2PComms(COUNT_PER_POINT);
  }

  /*--- Set some local pointers to make access simpler. ---*/

  su2double *bufDSend = geometry->bufD_P2PSend;

  /*--- Load the specified quantity from the solver into the generic
   communication buffer in the geometry class. ---*/

  if (geometry->nP2PSend > 0) {

    /*--- Post all non-blocking recvs first before sends. ---*/

    geometry->PostP2PRecvs(geometry, config, MPI_TYPE, reverse);

    for (iMessage = 0; iMessage < geometry->nP2PSend; iMessage++) {

      switch (commType) {

        case SOLUTION_MATRIX:

          /*--- Get the offset for the start of this message. ---*/

          msg_offset = geometry->nPoint_P2PSend[iMessage];

          /*--- Total count can include multiple pieces of data per point. ---*/

          nSend = (geometry->nPoint_P2PSend[iMessage+1] -
                   geometry->nPoint_P2PSend[iMessage]);

          for (iSend = 0; iSend < nSend; iSend++) {

            /*--- Get the local index for this communicated data. ---*/

            iPoint = geometry->Local_Point_P2PSend[msg_offset + iSend];

            /*--- Compute the offset in the recv buffer for this point. ---*/

            buf_offset = (msg_offset + iSend)*geometry->countPerPoint;

            /*--- Load the buffer with the data to be sent. ---*/

            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = x[iPoint*nVar+iVar];

          }

          break;

        case SOLUTION_MATRIXTRANS:

          /*--- We are going to communicate in reverse, so we use the
           recv buffer for the send instead. Also, all of the offsets
           and counts are derived from the recv data structures. ---*/

          bufDSend = geometry->bufD_P2PRecv;

          /*--- Get the offset for the start of this message. ---*/

          msg_offset = geometry->nPoint_P2PRecv[iMessage];

          /*--- Total count can include multiple pieces of data per point. ---*/

          nSend = (geometry->nPoint_P2PRecv[iMessage+1] -
                   geometry->nPoint_P2PRecv[iMessage]);

          for (iSend = 0; iSend < nSend; iSend++) {

            /*--- Get the local index for this communicated data. Here we
             again use the recv structure to find the send point, since
             the usual recv points are now the senders in reverse mode. ---*/

            iPoint = geometry->Local_Point_P2PRecv[msg_offset + iSend];

            /*--- Compute the offset in the recv buffer for this point. ---*/

            buf_offset = (msg_offset + iSend)*geometry->countPerPoint;

            /*--- Load the buffer with the data to be sent. ---*/

            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = x[iPoint*nVar+iVar];

          }

          break;

        default:
          SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.",
                         CURRENT_FUNCTION);
          break;

      }

      /*--- Launch the point-to-point MPI send for this message. ---*/

      geometry->PostP2PSends(geometry, config, MPI_TYPE, iMessage, reverse);

    }
  }

}

template<class ScalarType>
template<class OtherType>
void CSysMatrix<ScalarType>::CompleteComms(CSysVector<OtherType> & x,
                                           CGeometry *geometry,
                                           CConfig *config,
                                           unsigned short commType) const {

  /*--- Local variables ---*/

  unsigned short iVar;
  unsigned long iPoint, iRecv, nRecv, msg_offset, buf_offset;

  int ind, source, iMessage, jRecv;
  SU2_MPI::Status status;

  /*--- Set some local pointers to make access simpler. ---*/

  const su2double *bufDRecv = geometry->bufD_P2PRecv;

  /*--- Store the data that was communicated into the appropriate
   location within the local class data structures. ---*/

  if (geometry->nP2PRecv > 0) {

    for (iMessage = 0; iMessage < geometry->nP2PRecv; iMessage++) {

      /*--- For efficiency, recv the messages dynamically based on
       the order they arrive. ---*/

      SU2_MPI::Waitany(geometry->nP2PRecv, geometry->req_P2PRecv,
                       &ind, &status);

      /*--- Once we have recv'd a message, get the source rank. ---*/

      source = status.MPI_SOURCE;

      switch (commType) {
        case SOLUTION_MATRIX:

          /*--- We know the offsets based on the source rank. ---*/

          jRecv = geometry->P2PRecv2Neighbor[source];

          /*--- Get the offset for the start of this message. ---*/

          msg_offset = geometry->nPoint_P2PRecv[jRecv];

          /*--- Get the number of packets to be received in this message. ---*/

          nRecv = (geometry->nPoint_P2PRecv[jRecv+1] -
                   geometry->nPoint_P2PRecv[jRecv]);

          for (iRecv = 0; iRecv < nRecv; iRecv++) {

            /*--- Get the local index for this communicated data. ---*/

            iPoint = geometry->Local_Point_P2PRecv[msg_offset + iRecv];

            /*--- Compute the offset in the recv buffer for this point. ---*/

            buf_offset = (msg_offset + iRecv)*geometry->countPerPoint;

            /*--- Store the data correctly depending on the quantity. ---*/

            for (iVar = 0; iVar < nVar; iVar++)
              x[iPoint*nVar+iVar] = ActiveAssign<OtherType,su2double>(bufDRecv[buf_offset+iVar]);

          }
          break;

        case SOLUTION_MATRIXTRANS:

          /*--- We are going to communicate in reverse, so we use the
           send buffer for the recv instead. Also, all of the offsets
           and counts are derived from the send data structures. ---*/

          bufDRecv = geometry->bufD_P2PSend;

          /*--- We know the offsets based on the source rank. ---*/

          jRecv = geometry->P2PSend2Neighbor[source];

          /*--- Get the offset for the start of this message. ---*/

          msg_offset = geometry->nPoint_P2PSend[jRecv];

          /*--- Get the number of packets to be received in this message. ---*/

          nRecv = (geometry->nPoint_P2PSend[jRecv+1] -
                   geometry->nPoint_P2PSend[jRecv]);

          for (iRecv = 0; iRecv < nRecv; iRecv++) {

            /*--- Get the local index for this communicated data. ---*/

            iPoint = geometry->Local_Point_P2PSend[msg_offset + iRecv];

            /*--- Compute the offset in the recv buffer for this point. ---*/

            buf_offset = (msg_offset + iRecv)*geometry->countPerPoint;


            for (iVar = 0; iVar < nVar; iVar++)
              x[iPoint*nVar+iVar] += ActiveAssign<OtherType,su2double>(bufDRecv[buf_offset+iVar]);

          }

          break;
        default:
          SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.",
                         CURRENT_FUNCTION);
          break;
      }
    }

    /*--- Verify that all non-blocking point-to-point sends have finished.
     Note that this should be satisfied, as we have received all of the
     data in the loop above at this point. ---*/

#ifdef HAVE_MPI
    SU2_MPI::Waitall(geometry->nP2PSend, geometry->req_P2PSend, MPI_STATUS_IGNORE);
#endif

  }

}

template<class ScalarType>
void CSysMatrix<ScalarType>::SetValZero() {
  SU2_OMP_FOR_STAT(omp_light_size)
  for (auto index = 0ul; index < nnz*nVar*nEqn; index++)
    matrix[index] = 0.0;
}

template<class ScalarType>
void CSysMatrix<ScalarType>::SetValDiagonalZero() {
  SU2_OMP_FOR_STAT(omp_heavy_size)
  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint)
    for (auto index = 0ul; index < nVar*nEqn; ++index)
      matrix[dia_ptr[iPoint]*nVar*nEqn + index] = 0.0;
}

template<class ScalarType>
void CSysMatrix<ScalarType>::Gauss_Elimination(ScalarType* matrix, ScalarType* vec) const {

#ifdef USE_MKL_LAPACK
  // With MKL_DIRECT_CALL enabled, this is significantly faster than native code on Intel Architectures.
  lapack_int ipiv[MAXNVAR];
  LAPACKE_dgetrf( LAPACK_ROW_MAJOR, nVar, nVar, matrix, nVar, ipiv);
  LAPACKE_dgetrs( LAPACK_ROW_MAJOR, 'N', nVar, 1, matrix, nVar, ipiv, vec, 1 );
#else
#define A(I,J) matrix[(I)*nVar+(J)]

  /*--- Transform system in Upper Matrix ---*/
  for (auto iVar = 1ul; iVar < nVar; iVar++) {
    for (auto jVar = 0ul; jVar < iVar; jVar++) {
      ScalarType weight = A(iVar,jVar) / A(jVar,jVar);
      for (auto kVar = jVar; kVar < nVar; kVar++)
        A(iVar,kVar) -= weight * A(jVar,kVar);
      vec[iVar] -= weight * vec[jVar];
    }
  }

  /*--- Backwards substitution ---*/
  for (auto iVar = nVar; iVar > 0ul;) {
    iVar--; // unsigned type
    for (auto jVar = iVar+1; jVar < nVar; jVar++)
      vec[iVar] -= A(iVar,jVar) * vec[jVar];
    vec[iVar] /= A(iVar,iVar);
  }
#undef A
#endif
}

template<class ScalarType>
void CSysMatrix<ScalarType>::MatrixInverse(ScalarType *matrix, ScalarType *inverse) const {

  /*--- This is a generalization of Gaussian elimination for multiple rhs' (the basis vectors).
   We could call "Gauss_Elimination" multiple times or fully generalize it for multiple rhs,
   the performance of both routines would suffer in both cases without the use of exotic templating.
   And so it feels reasonable to have some duplication here. ---*/

  assert((matrix != inverse) && "Output cannot be the same as the input.");

#define M(I,J) inverse[(I)*nVar+(J)]

  /*--- Initialize the inverse with the identity. ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++)
    for (auto jVar = 0ul; jVar < nVar; jVar++)
      M(iVar,jVar) = ScalarType(iVar==jVar);

  /*--- Inversion ---*/
#ifdef USE_MKL_LAPACK
  // With MKL_DIRECT_CALL enabled, this is significantly faster than native code on Intel Architectures.
  lapack_int ipiv[MAXNVAR];
  LAPACKE_dgetrf( LAPACK_ROW_MAJOR, nVar, nVar, matrix, nVar, ipiv );
  LAPACKE_dgetrs( LAPACK_ROW_MAJOR, 'N', nVar, nVar, matrix, nVar, ipiv, inverse, nVar );
#else
#define A(I,J) matrix[(I)*nVar+(J)]

  /*--- Transform system in Upper Matrix ---*/
  for (auto iVar = 1ul; iVar < nVar; iVar++) {
    for (auto jVar = 0ul; jVar < iVar; jVar++)
    {
      ScalarType weight = A(iVar,jVar) / A(jVar,jVar);

      for (auto kVar = jVar; kVar < nVar; kVar++)
        A(iVar,kVar) -= weight * A(jVar,kVar);

      /*--- at this stage M is lower triangular so not all cols need updating ---*/
      for (auto kVar = 0ul; kVar <= jVar; kVar++)
        M(iVar,kVar) -= weight * M(jVar,kVar);
    }
  }

  /*--- Backwards substitution ---*/
  for (auto iVar = nVar; iVar > 0ul;) {
    iVar--; // unsigned type
    for (auto jVar = iVar+1; jVar < nVar; jVar++)
      for (auto kVar = 0ul; kVar < nVar; kVar++)
        M(iVar,kVar) -= A(iVar,jVar) * M(jVar,kVar);

    for (auto kVar = 0ul; kVar < nVar; kVar++)
      M(iVar,kVar) /= A(iVar,iVar);
  }
#undef A
#endif
#undef M
}

template<class ScalarType>
void CSysMatrix<ScalarType>::DeleteValsRowi(unsigned long i) {

  unsigned long block_i = i/nVar;
  unsigned long row = i - block_i*nVar;
  unsigned long index, iVar;

  for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
    for (iVar = 0; iVar < nVar; iVar++)
      matrix[index*nVar*nVar+row*nVar+iVar] = 0.0; // Delete row values in the block
    if (col_ind[index] == block_i)
      matrix[index*nVar*nVar+row*nVar+row] = 1.0; // Set 1 to the diagonal element
  }

}

template<class ScalarType>
void CSysMatrix<ScalarType>::RowProduct(const CSysVector<ScalarType> & vec,
                                        unsigned long row_i, ScalarType *prod) const {
  unsigned long iVar, index, col_j;

  for (iVar = 0; iVar < nVar; iVar++) prod[iVar] = 0.0;

  for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
    col_j = col_ind[index];
    MatrixVectorProductAdd(&matrix[index*nVar*nVar], &vec[col_j*nVar], prod);
  }

}

template<class ScalarType>
void CSysMatrix<ScalarType>::MatrixVectorProduct(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                                 CGeometry *geometry, CConfig *config) const {

  /*--- Some checks for consistency between CSysMatrix and the CSysVector<ScalarType>s ---*/
#ifndef NDEBUG
  if ( (nVar != vec.GetNVar()) || (nVar != prod.GetNVar()) ) {
    SU2_OMP_MASTER
    SU2_MPI::Error("nVar values incompatible.", CURRENT_FUNCTION);
  }
  if ( (nPoint != vec.GetNBlk()) || (nPoint != prod.GetNBlk()) ) {
    SU2_OMP_MASTER
    SU2_MPI::Error("nPoint and nBlk values incompatible.", CURRENT_FUNCTION);
  }
#endif

  /*--- OpenMP parallelization. First need to make view of vectors
   *    consistent, a barrier is implicit at the end of FOR section
   *    (and it is required before master thread communicates). ---*/

  SU2_OMP_BARRIER

  SU2_OMP_FOR_DYN(omp_heavy_size)
  for (auto row_i = 0ul; row_i < nPointDomain; row_i++) {
    auto prod_begin = row_i*nVar; // offset to beginning of block row_i
    for(auto iVar = 0ul; iVar < nVar; iVar++)
      prod[prod_begin+iVar] = 0.0;
    for (auto index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
      auto vec_begin = col_ind[index]*nVar; // offset to beginning of block col_ind[index]
      auto mat_begin = index*nVar*nVar; // offset to beginning of matrix block[row_i][col_ind[indx]]
      MatrixVectorProductAdd(&matrix[mat_begin], &vec[vec_begin], &prod[prod_begin]);
    }
  }

  /*--- MPI Parallelization by master thread. ---*/

  SU2_OMP_MASTER
  {
    InitiateComms(prod, geometry, config, SOLUTION_MATRIX);
    CompleteComms(prod, geometry, config, SOLUTION_MATRIX);
  }
  SU2_OMP_BARRIER
}

template<class ScalarType>
void CSysMatrix<ScalarType>::MatrixVectorProductTransposed(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                                           CGeometry *geometry, CConfig *config) const {

  unsigned long prod_begin, vec_begin, mat_begin, index, row_i;

  /*--- Some checks for consistency between CSysMatrix and the CSysVector<ScalarType>s ---*/
#ifndef NDEBUG
  if ( (nVar != vec.GetNVar()) || (nVar != prod.GetNVar()) ) {
    SU2_OMP_MASTER
    SU2_MPI::Error("nVar values incompatible.", CURRENT_FUNCTION);
  }
  if ( (nPoint != vec.GetNBlk()) || (nPoint != prod.GetNBlk()) ) {
    SU2_OMP_MASTER
    SU2_MPI::Error("nPoint and nBlk values incompatible.", CURRENT_FUNCTION);
  }
#endif

  /// TODO: The transpose product requires a different thread-parallel strategy.
  prod = ScalarType(0.0); // set all entries of prod to zero
  for (row_i = 0; row_i < nPointDomain; row_i++) {
    vec_begin = row_i*nVar; // offset to beginning of block col_ind[index]
    for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
      prod_begin = col_ind[index]*nVar; // offset to beginning of block row_i
      mat_begin = (index*nVar*nVar); // offset to beginning of matrix block[row_i][col_ind[indx]]
      MatrixVectorProductTransp(&matrix[mat_begin], &vec[vec_begin], &prod[prod_begin]);
    }
  }

  /*--- MPI Parallelization ---*/

  InitiateComms(prod, geometry, config, SOLUTION_MATRIXTRANS);
  CompleteComms(prod, geometry, config, SOLUTION_MATRIXTRANS);

}

template<class ScalarType>
void CSysMatrix<ScalarType>::BuildJacobiPreconditioner(bool transpose) {

  /*--- Build Jacobi preconditioner (M = D), compute and store the inverses of the diagonal blocks. ---*/
  SU2_OMP(for schedule(dynamic,omp_heavy_size) nowait)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
    InverseDiagonalBlock(iPoint, &(invM[iPoint*nVar*nVar]), transpose);

}

template<class ScalarType>
void CSysMatrix<ScalarType>::ComputeJacobiPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                                         CGeometry *geometry, CConfig *config) const {

  /*--- Apply Jacobi preconditioner, y = D^{-1} * x, the inverse of the diagonal is already known. ---*/
  SU2_OMP_BARRIER
  SU2_OMP_FOR_DYN(omp_heavy_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
    MatrixVectorProduct(&(invM[iPoint*nVar*nVar]), &vec[iPoint*nVar], &prod[iPoint*nVar]);

  /*--- MPI Parallelization ---*/
  SU2_OMP_MASTER
  {
    InitiateComms(prod, geometry, config, SOLUTION_MATRIX);
    CompleteComms(prod, geometry, config, SOLUTION_MATRIX);
  }
  SU2_OMP_BARRIER
}

template<class ScalarType>
void CSysMatrix<ScalarType>::BuildILUPreconditioner(bool transposed) {

  /*--- Copy block matrix to compute factorization in-place. ---*/

  if ((ilu_fill_in == 0) && !transposed) {
    /*--- ILU0, direct copy. ---*/
    SU2_OMP_FOR_STAT(omp_light_size)
    for (auto iVar = 0ul; iVar < nnz*nVar*nVar; ++iVar)
      ILU_matrix[iVar] = matrix[iVar];
  }
  else {
    /*--- ILUn clear the ILU matrix first, for ILU0^T
     *    the copy takes care of the clearing. ---*/
    if (ilu_fill_in > 0) {
      SU2_OMP_FOR_STAT(omp_light_size)
      for (auto iVar = 0ul; iVar < nnz_ilu*nVar*nVar; iVar++)
        ILU_matrix[iVar] = 0.0;
    }

    /*--- Transposed or ILUn, traverse matrix to access its blocks
     *    sequentially and set them in the ILU matrix. ---*/
    SU2_OMP_FOR_DYN(omp_heavy_size)
    for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {
      for (auto index = row_ptr[iPoint]; index < row_ptr[iPoint+1]; index++) {
        auto jPoint = col_ind[index];
        if (transposed) {
          SetBlockTransposed_ILUMatrix(jPoint, iPoint, &matrix[index*nVar*nVar]);
        } else {
          SetBlock_ILUMatrix(iPoint, jPoint, &matrix[index*nVar*nVar]);
        }
      }
    }
  }

  /*--- Transform system in Upper Matrix ---*/

  /*--- OpenMP Parallelization, a loop construct is used to ensure
   *    the preconditioner is computed correctly even if called
   *    outside of a parallel section. ---*/

  SU2_OMP_FOR_STAT(1)
  for(unsigned long thread = 0; thread < omp_num_parts; ++thread)
  {
    const auto begin = omp_partitions[thread];
    const auto end = omp_partitions[thread+1];

    /*--- Each thread will work on the submatrix defined from row/col "begin"
     *    to row/col "end-1" (i.e. the range [begin,end[). Which is exactly
     *    what the MPI-only implementation does. ---*/

    ScalarType weight[MAXNVAR*MAXNVAR], aux_block[MAXNVAR*MAXNVAR];

    for (auto iPoint = begin+1; iPoint < end; iPoint++) {

      /*--- Invert and store the previous diagonal block to later compute the weight. ---*/

      InverseDiagonalBlock_ILUMatrix(iPoint-1, &invM[(iPoint-1)*nVar*nVar]);

      /*--- For this row (unknown), loop over its lower diagonal entries. ---*/

      for (auto index = row_ptr_ilu[iPoint]; index < dia_ptr_ilu[iPoint]; index++) {

        /*--- jPoint is the column index (jPoint < iPoint). ---*/

        auto jPoint = col_ind_ilu[index];

        /*--- We only care about the sub matrix within "begin" and "end-1". ---*/

        if (jPoint < begin) continue;

        /*--- Multiply the block by the inverse of the corresponding diagonal block. ---*/

        auto Block_ij = &ILU_matrix[index*nVar*nVar];
        MatrixMatrixProduct(Block_ij, &invM[jPoint*nVar*nVar], weight);

        /*--- "weight" holds Aij*inv(Ajj). Jump to the upper part of the jPoint row. ---*/

        for (auto index_ = dia_ptr_ilu[jPoint]+1; index_ < row_ptr_ilu[jPoint+1]; index_++) {

          /*--- Get the column index (kPoint > jPoint). ---*/

          auto kPoint = col_ind_ilu[index_];

          if (kPoint >= end) break;

          /*--- If Aik exists, update it: Aik -= Aij*inv(Ajj)*Ajk ---*/

          auto Block_ik = GetBlock_ILUMatrix(iPoint, kPoint);

          if (Block_ik != nullptr) {
            auto Block_jk = &ILU_matrix[index_*nVar*nVar];
            MatrixMatrixProduct(weight, Block_jk, aux_block);
            MatrixSubtraction(Block_ik, aux_block, Block_ik);
          }
        }

        /*--- Lastly, store "weight" in the lower triangular part, which
         will be reused during the forward solve in the precon/smoother. ---*/

        for (auto iVar = 0ul; iVar < nVar*nVar; ++iVar)
          Block_ij[iVar] = weight[iVar];
      }
    }
    InverseDiagonalBlock_ILUMatrix(end-1, &invM[(end-1)*nVar*nVar]);

  } // end parallel

}

template<class ScalarType>
void CSysMatrix<ScalarType>::ComputeILUPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                                      CGeometry *geometry, CConfig *config) const {
  /*--- Coherent view of vectors. ---*/
  SU2_OMP_BARRIER

  /*--- OpenMP Parallelization ---*/
  SU2_OMP_FOR_STAT(1)
  for(unsigned long thread = 0; thread < omp_num_parts; ++thread)
  {
    const auto begin = omp_partitions[thread];
    const auto end = omp_partitions[thread+1];

    ScalarType aux_vec[MAXNVAR];

    /*--- Copy vector to then work on prod in place ---*/

    for (auto iVar = begin*nVar; iVar < end*nVar; iVar++)
      prod[iVar] = vec[iVar];

    /*--- Forward solve the system using the lower matrix entries that
     were computed and stored during the ILU preprocessing. Note
     that we are overwriting the residual vector as we go. ---*/

    for (auto iPoint = begin+1; iPoint < end; iPoint++) {
      for (auto index = row_ptr_ilu[iPoint]; index < dia_ptr_ilu[iPoint]; index++) {
        auto jPoint = col_ind_ilu[index];
        if (jPoint < begin) continue;
        auto Block_ij = &ILU_matrix[index*nVar*nVar];
        MatrixVectorProductSub(Block_ij, &prod[jPoint*nVar], &prod[iPoint*nVar]);
      }
    }

    /*--- Backwards substitution (starts at the last row) ---*/

    for (auto iPoint = end; iPoint > begin;) {
      iPoint--; // unsigned type
      for (auto iVar = 0ul; iVar < nVar; iVar++)
        aux_vec[iVar] = prod[iPoint*nVar+iVar];

      for (auto index = dia_ptr_ilu[iPoint]+1; index < row_ptr_ilu[iPoint+1]; index++) {
        auto jPoint = col_ind_ilu[index];
        if (jPoint >= end) break;
        auto Block_ij = &ILU_matrix[index*nVar*nVar];
        MatrixVectorProductSub(Block_ij, &prod[jPoint*nVar], aux_vec);
      }

      MatrixVectorProduct(&invM[iPoint*nVar*nVar], aux_vec, &prod[iPoint*nVar]);
    }
  } // end parallel

  /*--- MPI Parallelization ---*/

  SU2_OMP_MASTER
  {
    InitiateComms(prod, geometry, config, SOLUTION_MATRIX);
    CompleteComms(prod, geometry, config, SOLUTION_MATRIX);
  }
  SU2_OMP_BARRIER
}

template<class ScalarType>
void CSysMatrix<ScalarType>::ComputeLU_SGSPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                                         CGeometry *geometry, CConfig *config) const {

  /*--- First part of the symmetric iteration: (D+L).x* = b ---*/

  /*--- Coherent view of vectors. ---*/
  SU2_OMP_BARRIER

  /*--- OpenMP Parallelization ---*/
  SU2_OMP_FOR_STAT(1)
  for(unsigned long thread = 0; thread < omp_num_parts; ++thread)
  {
    const auto begin = omp_partitions[thread];
    const auto end = omp_partitions[thread+1];

    /*--- Each thread will work on the submatrix defined from row/col "begin"
     *    to row/col "end-1", except the last thread that also considers halos.
     *    This is NOT exactly equivalent to the MPI implementation on the same
     *    number of domains, for that we would need to define "thread-halos". ---*/

    ScalarType low_prod[MAXNVAR];

    for (auto iPoint = begin; iPoint < end; ++iPoint) {
      auto idx = iPoint*nVar;
      LowerProduct(prod, iPoint, begin, low_prod);        // Compute L.x*
      VectorSubtraction(&vec[idx], low_prod, &prod[idx]); // Compute y = b - L.x*
      Gauss_Elimination(iPoint, &prod[idx]);              // Solve D.x* = y
    }
  } // end parallel

  /*--- MPI Parallelization ---*/
  SU2_OMP_MASTER
  {
    InitiateComms(prod, geometry, config, SOLUTION_MATRIX);
    CompleteComms(prod, geometry, config, SOLUTION_MATRIX);
  }
  SU2_OMP_BARRIER

  /*--- Second part of the symmetric iteration: (D+U).x_(1) = D.x* ---*/

  /*--- OpenMP Parallelization ---*/
  SU2_OMP_FOR_STAT(1)
  for(unsigned long thread = 0; thread < omp_num_parts; ++thread)
  {
    const auto begin = omp_partitions[thread];
    const auto row_end = omp_partitions[thread+1];
    /*--- On the last thread partition the upper
     *    product should consider halo columns. ---*/
    const auto col_end = (row_end==nPointDomain)? nPoint : row_end;

    ScalarType up_prod[MAXNVAR], dia_prod[MAXNVAR];

    for (auto iPoint = row_end; iPoint > begin;) {
      iPoint--; // because of unsigned type
      auto idx = iPoint*nVar;
      DiagonalProduct(prod, iPoint, dia_prod);          // Compute D.x*
      UpperProduct(prod, iPoint, col_end, up_prod);     // Compute U.x_(n+1)
      VectorSubtraction(dia_prod, up_prod, &prod[idx]); // Compute y = D.x*-U.x_(n+1)
      Gauss_Elimination(iPoint, &prod[idx]);            // Solve D.x* = y
    }
  } // end parallel

  /*--- MPI Parallelization ---*/
  SU2_OMP_MASTER
  {
    InitiateComms(prod, geometry, config, SOLUTION_MATRIX);
    CompleteComms(prod, geometry, config, SOLUTION_MATRIX);
  }
  SU2_OMP_BARRIER
}

template<class ScalarType>
unsigned long CSysMatrix<ScalarType>::BuildLineletPreconditioner(CGeometry *geometry, CConfig *config) {

  assert(omp_get_thread_num()==0 && "Linelet preconditioner cannot be built by multiple threads.");

  bool add_point;
  unsigned long iEdge, iPoint, jPoint, index_Point, iLinelet, iVertex, next_Point, counter, iElem;
  unsigned short iMarker, iNode;
  su2double alpha = 0.9, weight, max_weight, area, volume_iPoint, volume_jPoint;
  const su2double* normal;
  unsigned long Local_nPoints, Local_nLineLets, Global_nPoints, Global_nLineLets, max_nElem;

  /*--- Memory allocation --*/

  vector<bool> check_Point(nPoint,true);

  LineletBool.clear();
  LineletBool.resize(nPoint,false);

  nLinelet = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == EULER_WALL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == DISPLACEMENT_BOUNDARY)) {
      nLinelet += geometry->nVertex[iMarker];
    }
  }

  /*--- If the domain contains well defined Linelets ---*/

  if (nLinelet != 0) {

    /*--- Basic initial allocation ---*/

    LineletPoint.resize(nLinelet);

    /*--- Define the basic linelets, starting from each vertex ---*/

    iLinelet = 0;

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
          (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
          (config->GetMarker_All_KindBC(iMarker) == EULER_WALL             ) ||
          (config->GetMarker_All_KindBC(iMarker) == DISPLACEMENT_BOUNDARY))
      {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          LineletPoint[iLinelet].push_back(iPoint);
          check_Point[iPoint] = false;
          iLinelet++;
        }
      }
    }

    /*--- Create the linelet structure ---*/

    iLinelet = 0;

    do {

      index_Point = 0;

      do {

        /*--- Compute the value of the max weight ---*/

        iPoint = LineletPoint[iLinelet][index_Point];
        max_weight = 0.0;
        for (iNode = 0; iNode < geometry->node[iPoint]->GetnPoint(); iNode++) {
          jPoint = geometry->node[iPoint]->GetPoint(iNode);
          if ((check_Point[jPoint]) && geometry->node[jPoint]->GetDomain()) {
            iEdge = geometry->FindEdge(iPoint, jPoint);
            normal = geometry->edges->GetNormal(iEdge);
            if (geometry->GetnDim() == 3) area = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
            else area = sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
            volume_iPoint = geometry->node[iPoint]->GetVolume();
            volume_jPoint = geometry->node[jPoint]->GetVolume();
            weight = 0.5*area*((1.0/volume_iPoint)+(1.0/volume_jPoint));
            max_weight = max(max_weight, weight);
          }
        }

        /*--- Verify if any face of the control volume must be added ---*/

        add_point = false;
        counter = 0;
        next_Point = geometry->node[iPoint]->GetPoint(0);
        for (iNode = 0; iNode < geometry->node[iPoint]->GetnPoint(); iNode++) {
          jPoint = geometry->node[iPoint]->GetPoint(iNode);
          iEdge = geometry->FindEdge(iPoint, jPoint);
          normal = geometry->edges->GetNormal(iEdge);
          if (geometry->GetnDim() == 3) area = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
          else area = sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
          volume_iPoint = geometry->node[iPoint]->GetVolume();
          volume_jPoint = geometry->node[jPoint]->GetVolume();
          weight = 0.5*area*((1.0/volume_iPoint)+(1.0/volume_jPoint));
          if (((check_Point[jPoint]) && (weight/max_weight > alpha) && (geometry->node[jPoint]->GetDomain())) &&
              ((index_Point == 0) || ((index_Point > 0) && (jPoint != LineletPoint[iLinelet][index_Point-1])))) {
            add_point = true;
            next_Point = jPoint;
            counter++;
          }
        }

        /*--- We have arrived to an isotropic zone ---*/

        if (counter > 1) add_point = false;

        /*--- Add a typical point to the linelet, no leading edge ---*/

        if (add_point) {
          LineletPoint[iLinelet].push_back(next_Point);
          check_Point[next_Point] = false;
          index_Point++;
        }

      } while (add_point);
      iLinelet++;
    } while (iLinelet < nLinelet);

    /*--- Identify the points that belong to a Linelet ---*/

    for (iLinelet = 0; iLinelet < nLinelet; iLinelet++) {
      for (iElem = 0; iElem < LineletPoint[iLinelet].size(); iElem++) {
        iPoint = LineletPoint[iLinelet][iElem];
        LineletBool[iPoint] = true;
      }
    }

    /*--- Identify the maximum number of elements in a Linelet ---*/

    max_nElem = LineletPoint[0].size();
    for (iLinelet = 1; iLinelet < nLinelet; iLinelet++)
      if (LineletPoint[iLinelet].size() > max_nElem)
        max_nElem = LineletPoint[iLinelet].size();

  }

  /*--- The domain doesn't have well defined linelets ---*/

  else {

    max_nElem = 0;

  }

  /*--- Screen output ---*/

  Local_nPoints = 0;
  for (iLinelet = 0; iLinelet < nLinelet; iLinelet++) {
    Local_nPoints += LineletPoint[iLinelet].size();
  }
  Local_nLineLets = nLinelet;

  SU2_MPI::Allreduce(&Local_nPoints, &Global_nPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nLineLets, &Global_nLineLets, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  /*--- Memory allocation --*/

  LineletUpper.resize(omp_get_max_threads(), vector<const ScalarType*>(max_nElem,nullptr));
  LineletVector.resize(omp_get_max_threads(), vector<ScalarType>(max_nElem*nVar,0.0));
  LineletInvDiag.resize(omp_get_max_threads(), vector<ScalarType>(max_nElem*nVar*nVar,0.0));

  return (unsigned long)(passivedouble(Global_nPoints) / Global_nLineLets);

}

template<class ScalarType>
void CSysMatrix<ScalarType>::ComputeLineletPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                                          CGeometry *geometry, CConfig *config) const {
  /*--- Coherent view of vectors. ---*/
  SU2_OMP_BARRIER

  /*--- Jacobi preconditioning where there is no linelet ---*/

  SU2_OMP(for schedule(dynamic,omp_heavy_size) nowait)
  for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++)
    if (!LineletBool[iPoint])
      MatrixVectorProduct(&(invM[iPoint*nVar*nVar]), &vec[iPoint*nVar], &prod[iPoint*nVar]);

  /*--- Solve each linelet using the Thomas algorithm ---*/

  SU2_OMP_FOR_DYN(1)
  for (auto iLinelet = 0ul; iLinelet < nLinelet; iLinelet++) {

    /*--- Get references to the working vectors allocated for this thread. ---*/

    int thread = omp_get_thread_num();
    vector<const ScalarType*>& lineletUpper = LineletUpper[thread];
    vector<ScalarType>& lineletInvDiag = LineletInvDiag[thread];
    vector<ScalarType>& lineletVector = LineletVector[thread];

    /*--- Initialize the solution vector with the rhs ---*/

    auto nElem = LineletPoint[iLinelet].size();

    for (auto iElem = 0ul; iElem < nElem; iElem++) {
      auto iPoint = LineletPoint[iLinelet][iElem];
      for (auto iVar = 0ul; iVar < nVar; iVar++)
        lineletVector[iElem*nVar+iVar] = vec[iPoint*nVar+iVar];
    }

    /*--- Forward pass, eliminate lower entries, modify diagonal and rhs. ---*/

    /*--- Small temporaries. ---*/
    ScalarType aux_block[MAXNVAR*MAXNVAR], aux_vector[MAXNVAR];

    /*--- Copy diagonal block for first point in this linelet. ---*/
    MatrixCopy(&matrix[dia_ptr[LineletPoint[iLinelet][0]]*nVar*nVar],
               lineletInvDiag.data());

    for (auto iElem = 1ul; iElem < nElem; iElem++) {

      /*--- Setup pointers to required matrices and vectors ---*/
      auto im1Point = LineletPoint[iLinelet][iElem-1];
      auto iPoint = LineletPoint[iLinelet][iElem];

      auto d = &matrix[dia_ptr[iPoint]*nVar*nVar];
      auto l = GetBlock(iPoint, im1Point);
      auto u = GetBlock(im1Point, iPoint);

      auto inv_dm1 = &lineletInvDiag[(iElem-1)*nVar*nVar];
      auto d_prime = &lineletInvDiag[iElem*nVar*nVar];
      auto b_prime = &lineletVector[iElem*nVar];

      /*--- Invert previous modified diagonal ---*/
      MatrixCopy(inv_dm1, aux_block);
      MatrixInverse(aux_block, inv_dm1);

      /*--- Left-multiply by lower block to obtain the weight ---*/
      MatrixMatrixProduct(l, inv_dm1, aux_block);

      /*--- Multiply weight by upper block to modify current diagonal ---*/
      MatrixMatrixProduct(aux_block, u, d_prime);
      MatrixSubtraction(d, d_prime, d_prime);

      /*--- Update the rhs ---*/
      MatrixVectorProduct(aux_block, &lineletVector[(iElem-1)*nVar], aux_vector);
      VectorSubtraction(b_prime, aux_vector, b_prime);

      /*--- Cache upper block pointer for the backward substitution phase ---*/
      lineletUpper[iElem-1] = u;
    }

    /*--- Backwards substitution, LineletVector becomes the solution ---*/

    /*--- x_n = d_n^{-1} * b_n ---*/
    Gauss_Elimination(&lineletInvDiag[(nElem-1)*nVar*nVar], &lineletVector[(nElem-1)*nVar]);

    /*--- x_i = d_i^{-1}*(b_i - u_i*x_{i+1}) ---*/
    for (auto iElem = nElem-1; iElem > 0; --iElem) {
      auto inv_dm1 = &lineletInvDiag[(iElem-1)*nVar*nVar];
      MatrixVectorProduct(lineletUpper[iElem-1], &lineletVector[iElem*nVar], aux_vector);
      VectorSubtraction(&lineletVector[(iElem-1)*nVar], aux_vector, aux_vector);
      MatrixVectorProduct(inv_dm1, aux_vector, &lineletVector[(iElem-1)*nVar]);
    }

    /*--- Copy results to product vector ---*/

    for (auto iElem = 0ul; iElem < nElem; iElem++) {
      auto iPoint = LineletPoint[iLinelet][iElem];
      for (auto iVar = 0ul; iVar < nVar; iVar++)
        prod[iPoint*nVar+iVar] = lineletVector[iElem*nVar+iVar];
    }

  }

  /*--- MPI Parallelization ---*/

  SU2_OMP_MASTER
  {
    InitiateComms(prod, geometry, config, SOLUTION_MATRIX);
    CompleteComms(prod, geometry, config, SOLUTION_MATRIX);
  }
  SU2_OMP_BARRIER

}

template<class ScalarType>
void CSysMatrix<ScalarType>::ComputeResidual(const CSysVector<ScalarType> & sol, const CSysVector<ScalarType> & f,
                                             CSysVector<ScalarType> & res) const {
  SU2_OMP_FOR_DYN(omp_heavy_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    ScalarType aux_vec[MAXNVAR];
    RowProduct(sol, iPoint, aux_vec);
    VectorSubtraction(aux_vec, &f[iPoint*nVar], &res[iPoint*nVar]);
  }
}

template<class ScalarType>
template<class OtherType>
void CSysMatrix<ScalarType>::EnforceSolutionAtNode(const unsigned long node_i, const OtherType *x_i, CSysVector<OtherType> & b) {

  /*--- Eliminate the row associated with node i (Block_ii = I and all other Block_ij = 0).
   *    To preserve eventual symmetry, also attempt to eliminate the column, if the sparse pattern is not
   *    symmetric the entire column may not be eliminated, the result (matrix and vector) is still correct.
   *    The vector is updated with the product of column i by the known (enforced) solution at node i. ---*/

  for (auto index = row_ptr[node_i]; index < row_ptr[node_i+1]; ++index) {

    auto node_j = col_ind[index];

    /*--- The diagonal block is handled outside the loop. ---*/
    if (node_j == node_i) continue;

    /*--- Delete block j on row i (bij) and ATTEMPT to delete block i on row j (bji). ---*/
    auto bij = &matrix[index*nVar*nVar];
    auto bji = GetBlock(node_j, node_i);

    /*--- The "attempt" part. ---*/
    if (bji == nullptr) {
      node_j = node_i;
      bji = bij;
    }

    for(auto iVar = 0ul; iVar < nVar; ++iVar) {
      for(auto jVar = 0ul; jVar < nVar; ++jVar) {
        /*--- Column product. ---*/
        b[node_j*nVar+iVar] -= bji[iVar*nVar+jVar] * x_i[jVar];
        /*--- Delete blocks. ---*/
        bij[iVar*nVar+jVar] = bji[iVar*nVar+jVar] = 0.0;
      }
    }

  }

  /*--- Set the diagonal block to the identity. ---*/
  SetVal2Diag(node_i, 1.0);

  /*--- Set know solution in rhs vector. ---*/
  b.SetBlock(node_i, x_i);

}

template<class ScalarType>
template<class OtherType>
void CSysMatrix<ScalarType>::EnforceSolutionAtDOF(unsigned long node_i, unsigned long iVar,
                                                  OtherType x_i, CSysVector<OtherType> & b) {

  for (auto index = row_ptr[node_i]; index < row_ptr[node_i+1]; ++index) {

    const auto node_j = col_ind[index];

    /*--- Delete row iVar of block j on row i (bij) and ATTEMPT
     *    to delete column iVar block i on row j (bji). ---*/

    auto bij = &matrix[index*nVar*nVar];
    auto bji = GetBlock(node_j, node_i);

    /*--- The "attempt" part. ---*/
    if (bji != nullptr) {
      for(auto jVar = 0ul; jVar < nVar; ++jVar) {
        /*--- Column product. ---*/
        b[node_j*nVar+jVar] -= bji[jVar*nVar+iVar] * x_i;
        /*--- Delete entries. ---*/
        bji[jVar*nVar+iVar] = 0.0;
      }
    }

    /*--- Delete row. ---*/
    for(auto jVar = 0ul; jVar < nVar; ++jVar)
      bij[iVar*nVar+jVar] = 0.0;

    /*--- Set the diagonal entry of the block to 1. ---*/
    if (node_j == node_i)
      bij[iVar*(nVar+1)] = 1.0;
  }

  /*--- Set know solution in rhs vector. ---*/
  b(node_i, iVar) = x_i;

}

template<class ScalarType>
void CSysMatrix<ScalarType>::SetDiagonalAsColumnSum() {

  SU2_OMP_FOR_DYN(omp_heavy_size)
  for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {

    auto block_ii = &matrix[dia_ptr[iPoint]*nVar*nEqn];

    for (auto k = 0ul; k < nVar*nEqn; ++k) block_ii[k] = 0.0;

    for (auto k = row_ptr[iPoint]; k < row_ptr[iPoint+1]; ++k) {
      auto block_ji = &matrix[col_ptr[k]*nVar*nEqn];
      if (block_ji != block_ii) MatrixSubtraction(block_ii, block_ji, block_ii);
    }
  }
}

template<class ScalarType>
void CSysMatrix<ScalarType>::MatrixMatrixAddition(ScalarType alpha, const CSysMatrix<ScalarType>& B) {

  /*--- Check that the sparse structure is shared between the two matrices,
   *    comparing pointers is ok as they are obtained from CGeometry. ---*/
  bool ok = (row_ptr == B.row_ptr) && (col_ind == B.col_ind) &&
            (nVar == B.nVar) && (nEqn == B.nEqn) && (nnz == B.nnz);

  if (!ok) {
    SU2_OMP_MASTER
    SU2_MPI::Error("Matrices do not have compatible sparsity.", CURRENT_FUNCTION);
  }

  SU2_OMP_FOR_STAT(omp_light_size)
  for (auto i = 0ul; i < nnz*nVar*nEqn; ++i)
    matrix[i] += alpha*B.matrix[i];

}

template<class ScalarType>
void CSysMatrix<ScalarType>::BuildPastixPreconditioner(CGeometry *geometry, CConfig *config,
                                                       unsigned short kind_fact, bool transposed) {
#ifdef HAVE_PASTIX
  /*--- Pastix will launch nested threads. ---*/
  SU2_OMP_MASTER
  {
    pastix_wrapper.SetMatrix(nVar,nPoint,nPointDomain,row_ptr,col_ind,matrix);
    pastix_wrapper.Factorize(geometry, config, kind_fact, transposed);
  }
  SU2_OMP_BARRIER
#else
  SU2_OMP_MASTER
  SU2_MPI::Error("SU2 was not compiled with -DHAVE_PASTIX", CURRENT_FUNCTION);
#endif
}

template<class ScalarType>
void CSysMatrix<ScalarType>::ComputePastixPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                                         CGeometry *geometry, CConfig *config) const {
#ifdef HAVE_PASTIX
  SU2_OMP_BARRIER
  SU2_OMP_MASTER
  {
    pastix_wrapper.Solve(vec,prod);
    InitiateComms(prod, geometry, config, SOLUTION_MATRIX);
    CompleteComms(prod, geometry, config, SOLUTION_MATRIX);
  }
  SU2_OMP_BARRIER
#else
  SU2_OMP_MASTER
  SU2_MPI::Error("SU2 was not compiled with -DHAVE_PASTIX", CURRENT_FUNCTION);
#endif
}

#if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)
template<>
void CSysMatrix<su2double>::BuildPastixPreconditioner(CGeometry *geometry, CConfig *config,
                                                      unsigned short kind_fact, bool transposed) {
  SU2_OMP_MASTER
  SU2_MPI::Error("The PaStiX preconditioner is only available in CSysMatrix<passivedouble>", CURRENT_FUNCTION);
}
template<>
void CSysMatrix<su2double>::ComputePastixPreconditioner(const CSysVector<su2double> & vec, CSysVector<su2double> & prod,
                                                        CGeometry *geometry, CConfig *config) const {
  SU2_OMP_MASTER
  SU2_MPI::Error("The PaStiX preconditioner is only available in CSysMatrix<passivedouble>", CURRENT_FUNCTION);
}
#endif

/*--- Explicit instantiations ---*/
template class CSysMatrix<su2double>;
template void  CSysMatrix<su2double>::InitiateComms(const CSysVector<su2double>&, CGeometry*, CConfig*, unsigned short) const;
template void  CSysMatrix<su2double>::CompleteComms(CSysVector<su2double>&, CGeometry*, CConfig*, unsigned short) const;
template void  CSysMatrix<su2double>::EnforceSolutionAtNode(unsigned long, const su2double*, CSysVector<su2double>&);
template void  CSysMatrix<su2double>::EnforceSolutionAtDOF(unsigned long, unsigned long, su2double, CSysVector<su2double>&);

#ifdef CODI_REVERSE_TYPE
template class CSysMatrix<passivedouble>;
template void  CSysMatrix<passivedouble>::InitiateComms(const CSysVector<passivedouble>&, CGeometry*, CConfig*, unsigned short) const;
template void  CSysMatrix<passivedouble>::InitiateComms(const CSysVector<su2double>&, CGeometry*, CConfig*, unsigned short) const;
template void  CSysMatrix<passivedouble>::CompleteComms(CSysVector<passivedouble>&, CGeometry*, CConfig*, unsigned short) const;
template void  CSysMatrix<passivedouble>::CompleteComms(CSysVector<su2double>&, CGeometry*, CConfig*, unsigned short) const;
template void  CSysMatrix<passivedouble>::EnforceSolutionAtNode(unsigned long, const passivedouble*, CSysVector<passivedouble>&);
template void  CSysMatrix<passivedouble>::EnforceSolutionAtNode(unsigned long, const su2double*, CSysVector<su2double>&);
template void  CSysMatrix<passivedouble>::EnforceSolutionAtDOF(unsigned long, unsigned long, passivedouble, CSysVector<passivedouble>&);
template void  CSysMatrix<passivedouble>::EnforceSolutionAtDOF(unsigned long, unsigned long, su2double, CSysVector<su2double>&);
#endif
