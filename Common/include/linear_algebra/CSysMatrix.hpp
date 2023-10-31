/*!
 * \file CSysMatrix.hpp
 * \brief Declaration of the block-sparse matrix class.
 *        The implemtation is in <i>CSysMatrix.cpp</i>.
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

#pragma once

#include "../../include/CConfig.hpp"
#include "CSysVector.hpp"
#include "CPastixWrapper.hpp"

#include <cstdlib>
#include <vector>
#include <cassert>

/*--- In forward mode the matrix is not of a built-in type. ---*/
#if defined(HAVE_MKL) && !defined(CODI_FORWARD_TYPE)
#include "mkl.h"
#ifndef __INTEL_MKL__
#error Could not determine the MKL version
#endif
/*--- JIT is only available since 2019. ---*/
#if __INTEL_MKL__ >= 2019
#define USE_MKL
/*---
 Lapack direct calls only seem to be created for Intel compilers, and it is not worthwhile
 making "getrf" and "getrs" compatible with AD since they are not used as often as "gemm".
---*/
#if defined(__INTEL_COMPILER) && defined(MKL_DIRECT_CALL_SEQ) && !defined(CODI_REVERSE_TYPE)
#define USE_MKL_LAPACK
#endif
template <class T>
struct mkl_jit_wrapper {
  using gemm_t = dgemm_jit_kernel_t;
  template <class... Ts>
  static void create_gemm(Ts&&... args) {
    mkl_jit_create_dgemm(args...);
  }
  static gemm_t get_gemm(void* jitter) { return mkl_jit_get_dgemm_ptr(jitter); }
};
template <>
struct mkl_jit_wrapper<float> {
  using gemm_t = sgemm_jit_kernel_t;
  template <class... Ts>
  static void create_gemm(Ts&&... args) {
    mkl_jit_create_sgemm(args...);
  }
  static gemm_t get_gemm(void* jitter) { return mkl_jit_get_sgemm_ptr(jitter); }
};
#else
#warning The current version of MKL does not support JIT gemm kernels
#endif
#endif

class CGeometry;

/*!
 * \brief Helper to communicate distributed vectors.
 * \ingroup SpLinSys
 */
struct CSysMatrixComms {
  /*!
   * \brief Routine to load a vector quantity into the data structures for MPI point-to-point
   *        communication and to launch non-blocking sends and recvs.
   * \param[in] x        - CSysVector holding the array of data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   */
  template <class T>
  static void Initiate(const CSysVector<T>& x, CGeometry* geometry, const CConfig* config,
                       unsigned short commType = SOLUTION_MATRIX);

  /*!
   * \brief Routine to complete the set of non-blocking communications launched by
   *        Initiate() and unpacking of the data in the vector.
   * \param[in] x        - CSysVector holding the array of data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be unpacked.
   */
  template <class T>
  static void Complete(CSysVector<T>& x, CGeometry* geometry, const CConfig* config,
                       unsigned short commType = SOLUTION_MATRIX);
};

/*!
 * \class CSysMatrix
 * \ingroup SpLinSys
 * \brief Main class for defining block-compressed-row-storage sparse matrices.
 */
template <class ScalarType>
class CSysMatrix {
 private:
  friend struct CSysMatrixComms;

  const int rank; /*!< \brief MPI Rank. */
  const int size; /*!< \brief MPI Size. */

  enum : size_t {
    MAXNVAR = 20
  }; /*!< \brief Maximum number of variables the matrix can handle. The static
                 size is needed for fast, per-thread, static memory allocation. */

  enum { OMP_MAX_SIZE_L = 8192 }; /*!< \brief Max. chunk size used in light parallel for loops. */
  enum { OMP_MAX_SIZE_H = 512 };  /*!< \brief Max. chunk size used in heavy parallel for loops. */
  enum { OMP_MIN_SIZE = 32 };     /*!< \brief Chunk size for finer grain operations. */
  unsigned long omp_light_size;   /*!< \brief Actual chunk size used in light loops (e.g. over non zeros). */
  unsigned long omp_heavy_size;   /*!< \brief Actual chunk size used in heavy loops (e.g. over rows). */
  unsigned long omp_num_parts;    /*!< \brief Number of threads used in thread-parallel LU_SGS and ILU. */
  unsigned long* omp_partitions;  /*!< \brief Point indexes of LU_SGS and ILU thread-parallel sub partitioning. */

  unsigned long nPoint;       /*!< \brief Number of points in the grid. */
  unsigned long nPointDomain; /*!< \brief Number of points in the grid (excluding halos). */
  unsigned long nVar;         /*!< \brief Number of variables (and rows of the blocks). */
  unsigned long nEqn;         /*!< \brief Number of equations (and columns of the blocks). */

  ScalarType* matrix;           /*!< \brief Entries of the sparse matrix. */
  unsigned long nnz;            /*!< \brief Number of possible nonzero entries in the matrix. */
  const unsigned long* row_ptr; /*!< \brief Pointers to the first element in each row. */
  const unsigned long* dia_ptr; /*!< \brief Pointers to the diagonal element in each row. */
  const unsigned long* col_ind; /*!< \brief Column index for each of the elements in val(). */
  const unsigned long* col_ptr; /*!< \brief The transpose of col_ind, pointer to blocks with the same column index. */

  ScalarType* ILU_matrix;           /*!< \brief Entries of the ILU sparse matrix. */
  unsigned long nnz_ilu;            /*!< \brief Number of possible nonzero entries in the matrix (ILU). */
  const unsigned long* row_ptr_ilu; /*!< \brief Pointers to the first element in each row (ILU). */
  const unsigned long* dia_ptr_ilu; /*!< \brief Pointers to the diagonal element in each row (ILU). */
  const unsigned long* col_ind_ilu; /*!< \brief Column index for each of the elements in val() (ILU). */
  unsigned short ilu_fill_in;       /*!< \brief Fill in level for the ILU preconditioner. */

  ScalarType* invM; /*!< \brief Inverse of (Jacobi) preconditioner, or diagonal of ILU. */

  /*--- Temporary (hence mutable) working memory used in the Linelet preconditioner, outer vector is for threads ---*/
  mutable vector<vector<const ScalarType*> >
      LineletUpper; /*!< \brief Pointers to the upper blocks of the tri-diag system (working memory). */
  mutable vector<vector<ScalarType> >
      LineletInvDiag; /*!< \brief Inverse of the diagonal blocks of the tri-diag system (working memory). */
  mutable vector<vector<ScalarType> >
      LineletVector; /*!< \brief Solution and RHS of the tri-diag system (working memory). */

#ifdef USE_MKL
  using gemm_t = typename mkl_jit_wrapper<ScalarType>::gemm_t;
  void* MatrixMatrixProductJitter;               /*!< \brief Jitter handle for MKL JIT based GEMM. */
  gemm_t MatrixMatrixProductKernel;              /*!< \brief MKL JIT based GEMM kernel. */
  void* MatrixVectorProductJitterBetaZero;       /*!< \brief Jitter handle for MKL JIT based GEMV. */
  gemm_t MatrixVectorProductKernelBetaZero;      /*!< \brief MKL JIT based GEMV kernel. */
  void* MatrixVectorProductJitterBetaOne;        /*!< \brief Jitter handle for MKL JIT based GEMV with BETA=1.0. */
  gemm_t MatrixVectorProductKernelBetaOne;       /*!< \brief MKL JIT based GEMV kernel with BETA=1.0. */
  void* MatrixVectorProductJitterAlphaMinusOne;  /*!< \brief Jitter handle for MKL JIT based GEMV with ALPHA=-1.0 and
                                                    BETA=1.0. */
  gemm_t MatrixVectorProductKernelAlphaMinusOne; /*!< \brief MKL JIT based GEMV kernel with ALPHA=-1.0 and BETA=1.0. */
#endif

#ifdef HAVE_PASTIX
  mutable CPastixWrapper<ScalarType> pastix_wrapper;
#endif

  /*!
   * \brief Auxilary object to wrap the edge map pointer used in fast block updates, i.e. without linear searches.
   */
  struct {
    const unsigned long* ptr = nullptr;
    unsigned long nEdge = 0;

    operator bool() { return nEdge != 0; }

    inline unsigned long operator()(unsigned long edge, unsigned long node) const { return ptr[2 * edge + node]; }
    inline unsigned long ij(unsigned long edge) const { return ptr[2 * edge]; }
    inline unsigned long ji(unsigned long edge) const { return ptr[2 * edge + 1]; }

  } edge_ptr;

  /*!
   * \brief Handle type conversion for when we Set, Add, etc. blocks, preserving derivative information (if supported by
   * types).
   */
  template <class DstType, class SrcType, su2enable_if<std::is_arithmetic<DstType>::value> = 0>
  FORCEINLINE static DstType ActiveAssign(const SrcType& val) {
    return SU2_TYPE::GetValue(val);
  }

  template <class DstType, class SrcType, su2enable_if<!std::is_arithmetic<DstType>::value> = 0>
  FORCEINLINE static DstType ActiveAssign(const SrcType& val) {
    return val;
  }

  /*!
   * \brief Handle type conversion for when we Set, Add, etc. blocks, discarding derivative information.
   */
  template <class SrcType>
  FORCEINLINE static ScalarType PassiveAssign(const SrcType& val) {
    return SU2_TYPE::GetValue(val);
  }

  /*!
   * \brief Calculates the matrix-vector product: product = matrix*vector
   * \param[in] matrix
   * \param[in] vector
   * \param[out] product
   */
  void MatrixVectorProduct(const ScalarType* matrix, const ScalarType* vector, ScalarType* product) const;

  /*!
   * \brief Calculates the matrix-vector product: product += matrix*vector
   * \param[in] matrix
   * \param[in] vector
   * \param[in,out] product
   */
  void MatrixVectorProductAdd(const ScalarType* matrix, const ScalarType* vector, ScalarType* product) const;

  /*!
   * \brief Calculates the matrix-vector product: product -= matrix*vector
   * \param[in] matrix
   * \param[in] vector
   * \param[in,out] product
   */
  void MatrixVectorProductSub(const ScalarType* matrix, const ScalarType* vector, ScalarType* product) const;

  /*!
   * \brief Calculates the matrix-matrix product
   */
  void MatrixMatrixProduct(const ScalarType* matrix_a, const ScalarType* matrix_b, ScalarType* product) const;

  /*!
   * \brief Subtract b from a and store the result in c.
   */
  FORCEINLINE void VectorSubtraction(const ScalarType* a, const ScalarType* b, ScalarType* c) const {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) c[iVar] = a[iVar] - b[iVar];
  }

  /*!
   * \brief Subtract b from a and store the result in c.
   */
  FORCEINLINE void MatrixSubtraction(const ScalarType* a, const ScalarType* b, ScalarType* c) const {
    SU2_OMP_SIMD
    for (unsigned long iVar = 0; iVar < nVar * nEqn; iVar++) c[iVar] = a[iVar] - b[iVar];
  }

  /*!
   * \brief Copy matrix src into dst, transpose if required.
   */
  FORCEINLINE void MatrixCopy(const ScalarType* src, ScalarType* dst) const {
    SU2_OMP_SIMD
    for (auto iVar = 0ul; iVar < nVar * nEqn; ++iVar) dst[iVar] = src[iVar];
  }

  /*!
   * \brief Solve a small (nVar x nVar) linear system using Gaussian elimination.
   * \param[in,out] matrix - On entry the system matrix, on exit the factorized matrix.
   * \param[in,out] vec - On entry the rhs, on exit the solution.
   */
  void Gauss_Elimination(ScalarType* matrix, ScalarType* vec) const;

  /*!
   * \brief Invert a small dense matrix.
   * \param[in,out] matrix - On entry the system matrix, on exit the factorized matrix.
   * \param[out] inverse - the matrix inverse.
   */
  void MatrixInverse(ScalarType* matrix, ScalarType* inverse) const;

  /*!
   * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i,i) subblock and rhs.
   * \param[in] block_i - Index of the (i,i) diagonal block.
   * \param[in] rhs - Right-hand-side of the linear system.
   * \return Solution of the linear system (overwritten on rhs).
   */
  inline void Gauss_Elimination(unsigned long block_i, ScalarType* rhs) const;

  /*!
   * \brief Inverse diagonal block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[out] invBlock - Inverse block.
   */
  inline void InverseDiagonalBlock(unsigned long block_i, ScalarType* invBlock) const;

  /*!
   * \brief Inverse diagonal block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[out] invBlock - Inverse block.
   */
  inline void InverseDiagonalBlock_ILUMatrix(unsigned long block_i, ScalarType* invBlock) const;

  /*!
   * \brief Copies the block (i, j) of the matrix-by-blocks structure in the internal variable *block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   */
  inline ScalarType* GetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j);

  /*!
   * \brief Set the value of a block in the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to set to A(i, j).
   */
  inline void SetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, ScalarType* val_block);

  /*!
   * \brief Performs the product of i-th row of the upper part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the upper part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \param[in] col_ub - Exclusive upper bound for column indices considered in multiplication.
   * \param[out] prod - Result of the product U(A)*vec.
   */
  inline void UpperProduct(const CSysVector<ScalarType>& vec, unsigned long row_i, unsigned long col_ub,
                           ScalarType* prod) const;

  /*!
   * \brief Performs the product of i-th row of the lower part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the lower part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \param[in] col_lb - Inclusive lower bound for column indices considered in multiplication.
   * \param[out] prod - Result of the product L(A)*vec.
   */
  inline void LowerProduct(const CSysVector<ScalarType>& vec, unsigned long row_i, unsigned long col_lb,
                           ScalarType* prod) const;

  /*!
   * \brief Performs the product of i-th row of the diagonal part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the diagonal part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return prod Result of the product D(A)*vec (stored at *prod_row_vector).
   */
  inline void DiagonalProduct(const CSysVector<ScalarType>& vec, unsigned long row_i, ScalarType* prod) const;

  /*!
   * \brief Performs the product of i-th row of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the row of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return Result of the product (stored at *prod_row_vector).
   */
  void RowProduct(const CSysVector<ScalarType>& vec, unsigned long row_i, ScalarType* prod) const;

 public:
  /*!
   * \brief Constructor of the class.
   */
  CSysMatrix(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CSysMatrix(void);

  /*!
   * \brief Initializes the sparse matrix.
   * \note The preconditioners require nVar == nEqn (square blocks).
   * \param[in] npoint - Number of points including halos.
   * \param[in] npointdomain - Number of points excluding halos.
   * \param[in] nvar - Number of variables (and rows of the blocks).
   * \param[in] neqn - Number of equations (and columns of the blocks).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] needTranspPtr - If "col_ptr" should be created, used for "SetDiagonalAsColumnSum".
   */
  void Initialize(unsigned long npoint, unsigned long npointdomain, unsigned short nvar, unsigned short neqn,
                  bool EdgeConnect, CGeometry* geometry, const CConfig* config, bool needTranspPtr = false,
                  bool grad_mode = false);

  /*!
   * \brief Sets to zero all the entries of the sparse matrix.
   */
  void SetValZero(void);

  /*!
   * \brief Sets to zero all the block diagonal entries of the sparse matrix.
   */
  void SetValDiagonalZero(void);

  /*!
   * \brief Get a pointer to the start of block "ij"
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \return Pointer to location in memory where the block starts.
   */
  FORCEINLINE const ScalarType* GetBlock(unsigned long block_i, unsigned long block_j) const {
    /*--- The position of the diagonal block is known which allows halving the search space. ---*/
    const auto end = (block_j < block_i) ? dia_ptr[block_i] : row_ptr[block_i + 1];
    for (auto index = (block_j < block_i) ? row_ptr[block_i] : dia_ptr[block_i]; index < end; ++index)
      if (col_ind[index] == block_j) return &matrix[index * nVar * nEqn];
    return nullptr;
  }

  /*!
   * \brief Get a pointer to the start of block "ij", non-const version
   */
  FORCEINLINE ScalarType* GetBlock(unsigned long block_i, unsigned long block_j) {
    const CSysMatrix& const_this = *this;
    return const_cast<ScalarType*>(const_this.GetBlock(block_i, block_j));
  }

  /*!
   * \brief Gets the value of a particular entry in block "ij".
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] iVar - Row of the block.
   * \param[in] jVar - Column of the block.
   * \return Value of the block entry.
   */
  FORCEINLINE ScalarType GetBlock(unsigned long block_i, unsigned long block_j, unsigned short iVar,
                                  unsigned short jVar) const {
    auto mat_ij = GetBlock(block_i, block_j);
    if (!mat_ij) return 0.0;
    return mat_ij[iVar * nEqn + jVar];
  }

  /*!
   * \brief Set the value of a block (in flat format) in the sparse matrix with scaling.
   * \note If the template param Overwrite is false we add to the block (bij += alpha*b).
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to set to A(i, j).
   * \param[in] alpha - Scale factor.
   */
  template <class OtherType, bool Overwrite = true, su2enable_if<!is_pointer<OtherType>::value> = 0>
  inline void SetBlock(unsigned long block_i, unsigned long block_j, const OtherType* val_block,
                       OtherType alpha = 1.0) {
    auto mat_ij = GetBlock(block_i, block_j);
    if (!mat_ij) return;
    SU2_OMP_SIMD
    for (auto iVar = 0ul; iVar < nVar * nEqn; ++iVar) {
      mat_ij[iVar] = (Overwrite ? ScalarType(0) : mat_ij[iVar]) + PassiveAssign(alpha * val_block[iVar]);
    }
  }

  /*!
   * \brief Add a scaled block (in flat format) to the sparse matrix (see SetBlock).
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to set to A(i, j).
   * \param[in] alpha - Scale factor.
   */
  template <class OtherType, su2enable_if<!is_pointer<OtherType>::value> = 0>
  inline void AddBlock(unsigned long block_i, unsigned long block_j, const OtherType* val_block,
                       OtherType alpha = 1.0) {
    SetBlock<OtherType, false>(block_i, block_j, val_block, alpha);
  }

  /*!
   * \brief Set the value of a scaled block in the sparse matrix.
   * \note If the template param Overwrite is false we add to the block (bij += alpha*b).
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to set to A(i, j).
   * \param[in] alpha - Scale factor.
   */
  template <class OtherType, bool Overwrite = true>
  inline void SetBlock(unsigned long block_i, unsigned long block_j, const OtherType* const* val_block,
                       OtherType alpha = 1.0) {
    auto mat_ij = GetBlock(block_i, block_j);
    if (!mat_ij) return;
    for (auto iVar = 0ul; iVar < nVar; ++iVar) {
      for (auto jVar = 0ul; jVar < nEqn; ++jVar) {
        *mat_ij = (Overwrite ? ScalarType(0) : *mat_ij) + PassiveAssign(alpha * val_block[iVar][jVar]);
        ++mat_ij;
      }
    }
  }

  /*!
   * \brief Adds a scaled block to the sparse matrix (see SetBlock).
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to add to A(i, j).
   * \param[in] alpha - Scale factor.
   */
  template <class OtherType>
  inline void AddBlock(unsigned long block_i, unsigned long block_j, const OtherType* const* val_block,
                       OtherType alpha = 1.0) {
    SetBlock<OtherType, false>(block_i, block_j, val_block, alpha);
  }

  /*!
   * \brief Subtracts the specified block to the sparse matrix (see AddBlock).
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to subtract to A(i, j).
   */
  template <class OtherType>
  inline void SubtractBlock(unsigned long block_i, unsigned long block_j, const OtherType* const* val_block) {
    AddBlock(block_i, block_j, val_block, OtherType(-1));
  }

  /*!
   * \brief Set the value of a scaled block in the sparse matrix.
   * \note This is an templated overload for C2Dcontainer specialization su2matrix.
   *       It assumes that MatrixType supports a member type Scalar and access operator[][].
   *       If the template param Overwrite is false we add to the block (bij += alpha*b).
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to set to A(i, j).
   * \param[in] alpha - Scale factor.
   */
  template <class MatrixType, bool Overwrite = true>
  inline void SetBlock(unsigned long block_i, unsigned long block_j, MatrixType& val_block,
                       typename MatrixType::Scalar alpha = 1.0) {
    auto mat_ij = GetBlock(block_i, block_j);
    if (!mat_ij) return;
    for (auto iVar = 0ul; iVar < nVar; ++iVar) {
      for (auto jVar = 0ul; jVar < nEqn; ++jVar) {
        *mat_ij = (Overwrite ? ScalarType(0) : *mat_ij) + PassiveAssign(alpha * val_block(iVar, jVar));
        ++mat_ij;
      }
    }
  }

  /*!
   * \brief Adds a scaled block to the sparse matrix (see SetBlock).
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to add to A(i, j).
   * \param[in] alpha - Scale factor.
   */
  template <class MatrixType>
  inline void AddBlock(unsigned long block_i, unsigned long block_j, MatrixType& val_block,
                       typename MatrixType::Scalar alpha = 1.0) {
    SetBlock<MatrixType, false>(block_i, block_j, val_block, alpha);
  }

  /*!
   * \brief Update 4 blocks ii, ij, ji, jj (add to i* sub from j*).
   * \note This method assumes an FVM-type sparse pattern.
   * \param[in] edge - Index of edge that connects iPoint and jPoint.
   * \param[in] iPoint - Row to which we add the blocks.
   * \param[in] jPoint - Row from which we subtract the blocks.
   * \param[in] block_i - Adds to ii, subs from ji.
   * \param[in] block_j - Adds to ij, subs from jj.
   * \param[in] scale - Scale blocks during update (axpy type op).
   */
  template <class MatrixType, class OtherType = ScalarType>
  inline void UpdateBlocks(unsigned long iEdge, unsigned long iPoint, unsigned long jPoint, const MatrixType& block_i,
                           const MatrixType& block_j, OtherType scale = 1) {
    ScalarType* bii = &matrix[dia_ptr[iPoint] * nVar * nEqn];
    ScalarType* bjj = &matrix[dia_ptr[jPoint] * nVar * nEqn];
    ScalarType* bij = &matrix[edge_ptr(iEdge, 0) * nVar * nEqn];
    ScalarType* bji = &matrix[edge_ptr(iEdge, 1) * nVar * nEqn];

    unsigned long iVar, jVar, offset = 0;

    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nEqn; jVar++) {
        bii[offset] += PassiveAssign(block_i[iVar][jVar] * scale);
        bij[offset] += PassiveAssign(block_j[iVar][jVar] * scale);
        bji[offset] -= PassiveAssign(block_i[iVar][jVar] * scale);
        bjj[offset] -= PassiveAssign(block_j[iVar][jVar] * scale);
        ++offset;
      }
    }
  }

  /*!
   * \brief Short-hand for the "subtractive" version (sub from i* add to j*) of UpdateBlocks.
   */
  template <class MatrixType>
  inline void UpdateBlocksSub(unsigned long iEdge, unsigned long iPoint, unsigned long jPoint,
                              const MatrixType& block_i, const MatrixType& block_j) {
    UpdateBlocks<MatrixType, ScalarType>(iEdge, iPoint, jPoint, block_i, block_j, -1);
  }

  /*!
   * \brief SIMD version, does the update for multiple edges and points.
   * \note Nothing is updated if the mask is 0.
   */
  template <class MatTypeSIMD, size_t N, class I, class F = ScalarType>
  FORCEINLINE void UpdateBlocks(simd::Array<I, N> iEdge, simd::Array<I, N> iPoint, simd::Array<I, N> jPoint,
                                const MatTypeSIMD& block_i, const MatTypeSIMD& block_j, simd::Array<F, N> mask = 1) {
    static_assert(MatTypeSIMD::StaticSize, "This method requires static size blocks.");
    static_assert(MatTypeSIMD::IsRowMajor, "Block storage is not compatible with matrix.");
    constexpr size_t blkSz = MatTypeSIMD::StaticSize;
    assert(blkSz == nVar * nEqn);

    /*--- "Transpose" the blocks, scale, and possibly convert types,
     * giving the compiler the chance to vectorize all of these. ---*/
    ScalarType blk_i[N][blkSz], blk_j[N][blkSz];

    for (size_t i = 0; i < blkSz; ++i) {
      SU2_OMP_SIMD_IF_NOT_AD
      for (size_t k = 0; k < N; ++k) {
        blk_i[k][i] = PassiveAssign(-mask[k] * block_i.data()[i][k]);
        blk_j[k][i] = PassiveAssign(mask[k] * block_j.data()[i][k]);
      }
    }

    /*--- Update one by one skipping if mask is 0. ---*/
    for (size_t k = 0; k < N; ++k) {
      if (mask[k] == 0) continue;

      /*--- Fetch the blocks. ---*/
      auto bii = &matrix[dia_ptr[iPoint[k]] * blkSz];
      auto bjj = &matrix[dia_ptr[jPoint[k]] * blkSz];
      auto bij = &matrix[edge_ptr(iEdge[k], 0) * blkSz];
      auto bji = &matrix[edge_ptr(iEdge[k], 1) * blkSz];

      /*--- Update, block i was negated during transpose in the
       * hope the assignments below become non-temporal stores. ---*/
      SU2_OMP_SIMD
      for (size_t i = 0; i < blkSz; ++i) {
        bii[i] -= blk_i[k][i];
        bjj[i] -= blk_j[k][i];
        bij[i] = blk_j[k][i];
        bji[i] = blk_i[k][i];
      }
    }
  }

  /*!
   * \brief Sets 2 blocks ij and ji (add to i* sub from j*) associated with
   *        one edge of an FVM-type sparse pattern.
   * \note The parameter Overwrite allows completely writing over the
   *       current values held by the matrix (true), or updating them (false).
   * \param[in] edge - Index of edge that connects iPoint and jPoint.
   * \param[in] block_i - Subs from ji.
   * \param[in] block_j - Adds to ij.
   * \param[in] scale - Scale blocks during update (axpy type op).
   */
  template <class MatrixType, class OtherType = ScalarType, bool Overwrite = true>
  inline void SetBlocks(unsigned long iEdge, const MatrixType& block_i, const MatrixType& block_j,
                        OtherType scale = 1) {
    ScalarType* bij = &matrix[edge_ptr(iEdge, 0) * nVar * nEqn];
    ScalarType* bji = &matrix[edge_ptr(iEdge, 1) * nVar * nEqn];

    unsigned long iVar, jVar, offset = 0;

    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nEqn; jVar++) {
        bij[offset] = (Overwrite ? ScalarType(0) : bij[offset]) + PassiveAssign(block_j[iVar][jVar] * scale);
        bji[offset] = (Overwrite ? ScalarType(0) : bji[offset]) - PassiveAssign(block_i[iVar][jVar] * scale);
        ++offset;
      }
    }
  }

  /*!
   * \brief Short-hand for the "additive overwrite" version of SetBlocks.
   */
  template <class MatrixType, class OtherType = ScalarType>
  inline void UpdateBlocks(unsigned long iEdge, const MatrixType& block_i, const MatrixType& block_j,
                           OtherType scale = 1) {
    SetBlocks<MatrixType, OtherType, false>(iEdge, block_i, block_j, scale);
  }

  /*!
   * \brief Short-hand for the "subtractive" version (sub from i* add to j*) of SetBlocks.
   */
  template <class MatrixType>
  inline void UpdateBlocksSub(unsigned long iEdge, const MatrixType& block_i, const MatrixType& block_j) {
    SetBlocks<MatrixType, ScalarType, false>(iEdge, block_i, block_j, -1);
  }

  /*!
   * \brief SIMD version, does the update for multiple edges.
   * \note Nothing is updated if the mask is 0.
   */
  template <class MatTypeSIMD, size_t N, class I, class F = ScalarType>
  FORCEINLINE void SetBlocks(simd::Array<I, N> iEdge, const MatTypeSIMD& block_i, const MatTypeSIMD& block_j,
                             simd::Array<F, N> mask = 1) {
    static_assert(MatTypeSIMD::StaticSize, "This method requires static size blocks.");
    static_assert(MatTypeSIMD::IsRowMajor, "Block storage is not compatible with matrix.");
    constexpr size_t blkSz = MatTypeSIMD::StaticSize;
    assert(blkSz == nVar * nEqn);

    /*--- "Transpose" the blocks, scale, and possibly convert types,
     * giving the compiler the chance to vectorize all of these. ---*/
    ScalarType blk_i[N][blkSz], blk_j[N][blkSz];

    for (size_t i = 0; i < blkSz; ++i) {
      SU2_OMP_SIMD_IF_NOT_AD
      for (size_t k = 0; k < N; ++k) {
        blk_i[k][i] = PassiveAssign(-mask[k] * block_i.data()[i][k]);
        blk_j[k][i] = PassiveAssign(mask[k] * block_j.data()[i][k]);
      }
    }

    /*--- Update one by one skipping if mask is 0. ---*/
    for (size_t k = 0; k < N; ++k) {
      if (mask[k] == 0) continue;

      /*--- Fetch the blocks. ---*/
      auto bij = &matrix[edge_ptr(iEdge[k], 0) * blkSz];
      auto bji = &matrix[edge_ptr(iEdge[k], 1) * blkSz];

      /*--- Update, block i was negated during transpose in the
       * hope the assignments below become non-temporal stores. ---*/
      SU2_OMP_SIMD
      for (size_t i = 0; i < blkSz; ++i) {
        bij[i] = blk_j[k][i];
        bji[i] = blk_i[k][i];
      }
    }
  }

  /*!
   * \brief Sets the specified block to the (i, i) subblock of the sparse matrix.
   *        Scales the input block by factor alpha. If the Overwrite parameter is
   *        false we update instead (bii += alpha*b).
   * \param[in] block_i - Diagonal index.
   * \param[in] val_block - Block to add to the diagonal of the matrix.
   * \param[in] alpha - Scale factor.
   */
  template <class OtherType, bool Overwrite = true, class T = ScalarType>
  inline void SetBlock2Diag(unsigned long block_i, const OtherType& val_block, T alpha = 1.0) {
    auto mat_ii = &matrix[dia_ptr[block_i] * nVar * nEqn];

    for (auto iVar = 0ul; iVar < nVar; iVar++)
      for (auto jVar = 0ul; jVar < nEqn; jVar++) {
        *mat_ii = (Overwrite ? ScalarType(0) : *mat_ii) + PassiveAssign(alpha * val_block[iVar][jVar]);
        ++mat_ii;
      }
  }

  /*!
   * \brief Non overwrite version of SetBlock2Diag, also with scaling.
   */
  template <class OtherType, class T = ScalarType>
  inline void AddBlock2Diag(unsigned long block_i, const OtherType& val_block, T alpha = 1.0) {
    SetBlock2Diag<OtherType, false>(block_i, val_block, alpha);
  }

  /*!
   * \brief Short-hand to AddBlock2Diag with alpha = -1, i.e. subtracts from the current diagonal.
   */
  template <class OtherType>
  inline void SubtractBlock2Diag(unsigned long block_i, const OtherType& val_block) {
    AddBlock2Diag(block_i, val_block, -1.0);
  }

  /*!
   * \brief Adds the specified value to the diagonal of the (i, i) subblock
   *        of the matrix-by-blocks structure.
   * \param[in] block_i - Diagonal index.
   * \param[in] val_matrix - Value to add to the diagonal elements of A(i, i).
   */
  template <class OtherType>
  inline void AddVal2Diag(unsigned long block_i, OtherType val_matrix) {
    for (auto iVar = 0ul; iVar < nVar; iVar++)
      matrix[dia_ptr[block_i] * nVar * nVar + iVar * (nVar + 1)] += PassiveAssign(val_matrix);
  }

  /*!
   * \brief Adds the specified value to the diagonal of the (i, i) subblock
   *        of the matrix-by-blocks structure.
   * \param[in] block_i - Diagonal index.
   * \param[in] iVar - Variable index.
   * \param[in] val - Value to add to the diagonal elements of A(i, i).
   */
  template <class OtherType>
  inline void AddVal2Diag(unsigned long block_i, unsigned long iVar, OtherType val) {
    matrix[dia_ptr[block_i] * nVar * nVar + iVar * (nVar + 1)] += PassiveAssign(val);
  }

  /*!
   * \brief Sets the specified value to the diagonal of the (i, i) subblock
   *        of the matrix-by-blocks structure.
   * \param[in] block_i - Diagonal index.
   * \param[in] val_matrix - Value to add to the diagonal elements of A(i, i).
   */
  template <class OtherType>
  inline void SetVal2Diag(unsigned long block_i, OtherType val_matrix) {
    unsigned long iVar, index = dia_ptr[block_i] * nVar * nVar;

    /*--- Clear entire block before setting its diagonal. ---*/
    SU2_OMP_SIMD
    for (iVar = 0; iVar < nVar * nVar; iVar++) matrix[index + iVar] = 0.0;

    for (iVar = 0; iVar < nVar; iVar++) matrix[index + iVar * (nVar + 1)] = PassiveAssign(val_matrix);
  }

  /*!
   * \brief Deletes the values of the row i of the sparse matrix.
   * \param[in] i - Index of the row.
   */
  void DeleteValsRowi(unsigned long i);

  /*!
   * \brief Modifies this matrix (A) and a rhs vector (b) such that (A^-1 * b)_i = x_i.
   * \param[in] node_i - Index of the node for which to enforce the solution of all DOF's.
   * \param[in] x_i - Values to enforce (nVar sized).
   * \param[in,out] b - The rhs vector (b := b - A_{*,i} * x_i;  b_i = x_i).
   */
  template <class OtherType>
  void EnforceSolutionAtNode(unsigned long node_i, const OtherType* x_i, CSysVector<OtherType>& b);

  /*!
   * \brief Version of EnforceSolutionAtNode for a single degree of freedom.
   */
  template <class OtherType>
  void EnforceSolutionAtDOF(unsigned long node_i, unsigned long iVar, OtherType x_i, CSysVector<OtherType>& b);

  /*!
   * \brief Sets the diagonal entries of the matrix as the sum of the blocks in the corresponding column.
   */
  void SetDiagonalAsColumnSum();

  /*!
   * \brief Transposes the matrix, any preconditioner that was computed may be invalid.
   */
  void TransposeInPlace();

  /*!
   * \brief Add a scaled sparse matrix to "this" (axpy-type operation, A = A+alpha*B).
   * \note Matrices must have the same sparse pattern.
   * \param[in] alpha - The scaling constant.
   * \param[in] B - Matrix being.
   */
  void MatrixMatrixAddition(ScalarType alpha, const CSysMatrix& B);

  /*!
   * \brief Performs the product of a sparse matrix by a CSysVector.
   * \param[in] vec - CSysVector to be multiplied by the sparse matrix A.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] prod - Result of the product.
   */
  void MatrixVectorProduct(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod, CGeometry* geometry,
                           const CConfig* config) const;

  /*!
   * \brief Build the Jacobi preconditioner.
   */
  void BuildJacobiPreconditioner();

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeJacobiPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod, CGeometry* geometry,
                                   const CConfig* config) const;

  /*!
   * \brief Build the ILU preconditioner.
   */
  void BuildILUPreconditioner();

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeILUPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod, CGeometry* geometry,
                                const CConfig* config) const;

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   */
  void ComputeLU_SGSPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod, CGeometry* geometry,
                                   const CConfig* config) const;

  /*!
   * \brief Build the Linelet preconditioner.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void BuildLineletPreconditioner(const CGeometry* geometry, const CConfig* config);

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   */
  void ComputeLineletPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                    CGeometry* geometry, const CConfig* config) const;

  /*!
   * \brief Compute the linear residual.
   * \param[in] sol - Solution (x).
   * \param[in] f - Right hand side (b).
   * \param[out] res - Residual (Ax-b).
   */
  void ComputeResidual(const CSysVector<ScalarType>& sol, const CSysVector<ScalarType>& f,
                       CSysVector<ScalarType>& res) const;

  /*!
   * \brief Factorize matrix using PaStiX.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] kind_fact - Type of factorization.
   */
  void BuildPastixPreconditioner(CGeometry* geometry, const CConfig* config, unsigned short kind_fact);

  /*!
   * \brief Apply the PaStiX factorization to CSysVec.
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product M*vec.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputePastixPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod, CGeometry* geometry,
                                   const CConfig* config) const;
};
