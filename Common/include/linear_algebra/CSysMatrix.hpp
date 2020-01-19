/*!
 * \file CSysMatrix.hpp
 * \brief Declaration of the block-sparse matrix class.
 *        The implemtation is in <i>CSysMatrix.cpp</i>.
 * \author F. Palacios, A. Bueno, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/mpi_structure.hpp"
#include "CSysVector.hpp"
#include "CPastixWrapper.hpp"

#include <cstdlib>
#include <vector>

using namespace std;

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
#else
  #warning The current version of MKL does not support JIT gemm kernels
#endif
#endif

class CConfig;
class CGeometry;

/*!
 * \class CSysMatrix
 * \brief Main class for defining block-compressed-row-storage sparse matrices.
 * \author A. Bueno, F. Palacios
 */
template<class ScalarType>
class CSysMatrix {
private:
  /*--- We are friends with all other possible CSysMatrices. ---*/
  template<class T> friend class CSysMatrix;

  int rank;     /*!< \brief MPI Rank. */
  int size;     /*!< \brief MPI Size. */

  enum : size_t { MAXNVAR = 8 };    /*!< \brief Maximum number of variables the matrix can handle. The static
                                                size is needed for fast, per-thread, static memory allocation. */

  enum { OMP_MAX_SIZE_L = 8192 };   /*!< \brief Max. chunk size used in light parallel for loops. */
  enum { OMP_MAX_SIZE_H = 512 };    /*!< \brief Max. chunk size used in heavy parallel for loops. */
  unsigned long omp_light_size;     /*!< \brief Actual chunk size used in light loops (e.g. over non zeros). */
  unsigned long omp_heavy_size;     /*!< \brief Actual chunk size used in heavy loops (e.g. over rows). */
  unsigned long omp_num_parts;      /*!< \brief Number of threads used in thread-parallel LU_SGS and ILU. */
  unsigned long *omp_partitions;    /*!< \brief Point indexes of LU_SGS and ILU thread-parallel sub partitioning. */

  unsigned long nPoint;             /*!< \brief Number of points in the grid. */
  unsigned long nPointDomain;       /*!< \brief Number of points in the grid (excluding halos). */
  unsigned long nVar;               /*!< \brief Number of variables. */
  unsigned long nEqn;               /*!< \brief Number of equations. */

  ScalarType *matrix;               /*!< \brief Entries of the sparse matrix. */
  unsigned long nnz;                /*!< \brief Number of possible nonzero entries in the matrix. */
  const unsigned long *row_ptr;     /*!< \brief Pointers to the first element in each row. */
  const unsigned long *dia_ptr;     /*!< \brief Pointers to the diagonal element in each row. */
  const unsigned long *col_ind;     /*!< \brief Column index for each of the elements in val(). */

  ScalarType *ILU_matrix;           /*!< \brief Entries of the ILU sparse matrix. */
  unsigned long nnz_ilu;            /*!< \brief Number of possible nonzero entries in the matrix (ILU). */
  const unsigned long *row_ptr_ilu; /*!< \brief Pointers to the first element in each row (ILU). */
  const unsigned long *dia_ptr_ilu; /*!< \brief Pointers to the diagonal element in each row (ILU). */
  const unsigned long *col_ind_ilu; /*!< \brief Column index for each of the elements in val() (ILU). */
  unsigned short ilu_fill_in;       /*!< \brief Fill in level for the ILU preconditioner. */

  ScalarType *invM;                 /*!< \brief Inverse of (Jacobi) preconditioner, or diagonal of ILU. */

  unsigned long nLinelet;                      /*!< \brief Number of Linelets in the system. */
  vector<bool> LineletBool;                    /*!< \brief Identify if a point belong to a Linelet. */
  vector<vector<unsigned long> > LineletPoint; /*!< \brief Linelet structure. */

  /*--- Temporary (hence mutable) working memory used in the Linelet preconditioner, outer vector is for threads ---*/
  mutable vector<vector<const ScalarType*> > LineletUpper; /*!< \brief Pointers to the upper blocks of the tri-diag system (working memory). */
  mutable vector<vector<ScalarType> > LineletInvDiag;      /*!< \brief Inverse of the diagonal blocks of the tri-diag system (working memory). */
  mutable vector<vector<ScalarType> > LineletVector;       /*!< \brief Solution and RHS of the tri-diag system (working memory). */

#ifdef USE_MKL
  void * MatrixMatrixProductJitter;                            /*!< \brief Jitter handle for MKL JIT based GEMM. */
  dgemm_jit_kernel_t MatrixMatrixProductKernel;                /*!< \brief MKL JIT based GEMM kernel. */
  void * MatrixVectorProductJitterBetaZero;                    /*!< \brief Jitter handle for MKL JIT based GEMV. */
  dgemm_jit_kernel_t MatrixVectorProductKernelBetaZero;        /*!< \brief MKL JIT based GEMV kernel. */
  void * MatrixVectorProductJitterBetaOne;                     /*!< \brief Jitter handle for MKL JIT based GEMV with BETA=1.0. */
  dgemm_jit_kernel_t MatrixVectorProductKernelBetaOne;         /*!< \brief MKL JIT based GEMV kernel with BETA=1.0. */
  void * MatrixVectorProductJitterAlphaMinusOne;               /*!< \brief Jitter handle for MKL JIT based GEMV with ALPHA=-1.0 and BETA=1.0. */
  dgemm_jit_kernel_t MatrixVectorProductKernelAlphaMinusOne;   /*!< \brief MKL JIT based GEMV kernel with ALPHA=-1.0 and BETA=1.0. */
  void * MatrixVectorProductTranspJitterBetaOne;               /*!< \brief Jitter handle for MKL JIT based GEMV (transposed) with BETA=1.0. */
  dgemm_jit_kernel_t MatrixVectorProductTranspKernelBetaOne;   /*!< \brief MKL JIT based GEMV (transposed) kernel with BETA=1.0. */
#endif

#ifdef HAVE_PASTIX
  mutable CPastixWrapper pastix_wrapper;
#endif

  /*!
   * \brief Auxilary object to wrap the edge map pointer used in fast block updates, i.e. without linear searches.
   */
  struct {
    const unsigned long *ptr = nullptr;

    inline unsigned long operator() (unsigned long edge, unsigned long node) const {
      return ptr[2*edge+node];
    }
    inline unsigned long ij(unsigned long edge) const { return ptr[2*edge]; }
    inline unsigned long ji(unsigned long edge) const { return ptr[2*edge+1]; }

  } edge_ptr;

  /*!
   * \brief Handle type conversion for when we Set, Add, etc. blocks, preserving derivative information (if supported by types).
   * \note See specializations for discrete adjoint right outside this class's declaration.
   */
  template<class DstType, class SrcType>
  inline DstType ActiveAssign(const SrcType & val) const { return val; }

  /*!
   * \brief Handle type conversion for when we Set, Add, etc. blocks, discarding derivative information.
   */
  template<class DstType, class SrcType>
  inline DstType PassiveAssign(const SrcType & val) const {
#if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)
    return SU2_TYPE::GetValue(val);
#else
    return val;
#endif
  }

  /*!
   * \brief Calculates the matrix-vector product: product = matrix*vector
   * \param[in] matrix
   * \param[in] vector
   * \param[out] product
   */
  inline void MatrixVectorProduct(const ScalarType *matrix, const ScalarType *vector, ScalarType *product) const;

  /*!
   * \brief Calculates the matrix-vector product: product += matrix*vector
   * \param[in] matrix
   * \param[in] vector
   * \param[in,out] product
   */
  inline void MatrixVectorProductAdd(const ScalarType *matrix, const ScalarType *vector, ScalarType *product) const;

  /*!
   * \brief Calculates the matrix-vector product: product -= matrix*vector
   * \param[in] matrix
   * \param[in] vector
   * \param[in,out] product
   */
  inline void MatrixVectorProductSub(const ScalarType *matrix, const ScalarType *vector, ScalarType *product) const;

  /*!
   * \brief Calculates the matrix-vector product: product += matrix^T * vector
   * \param[in] matrix
   * \param[in] vector
   * \param[in,out] product
   */
  inline void MatrixVectorProductTransp(const ScalarType *matrix, const ScalarType *vector, ScalarType *product) const;

  /*!
   * \brief Calculates the matrix-matrix product
   */
  inline void MatrixMatrixProduct(const ScalarType *matrix_a, const ScalarType *matrix_b, ScalarType *product) const;

  /*!
   * \brief Subtract b from a and store the result in c.
   */
  inline void VectorSubtraction(const ScalarType *a, const ScalarType *b, ScalarType *c) const {
    for(unsigned long iVar = 0; iVar < nVar; iVar++)
      c[iVar] = a[iVar] - b[iVar];
  }

  /*!
   * \brief Subtract b from a and store the result in c.
   */
  inline void MatrixSubtraction(const ScalarType *a, const ScalarType *b, ScalarType *c) const {
    for(unsigned long iVar = 0; iVar < nVar*nEqn; iVar++)
      c[iVar] = a[iVar] - b[iVar];
  }

  /*!
   * \brief Copy matrix src into dst, transpose is required.
   */
  inline void MatrixCopy(const ScalarType *src, ScalarType *dst, bool transposed = false) const {
    if (!transposed) {
      for(auto iVar = 0ul; iVar < nVar*nVar; ++iVar)
        dst[iVar] = src[iVar];
    }
    else {
      for (auto iVar = 0ul; iVar < nVar; ++iVar)
        for (auto jVar = 0ul; jVar < nVar; ++jVar)
          dst[iVar*nVar+jVar] = src[jVar*nVar+iVar];
    }
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
  void MatrixInverse(ScalarType *matrix, ScalarType *inverse) const;

  /*!
   * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i,i) subblock and rhs.
   * \param[in] block_i - Index of the (i,i) diagonal block.
   * \param[in] rhs - Right-hand-side of the linear system.
   * \param[in] transposed - If true the transposed of the block is used (default = false).
   * \return Solution of the linear system (overwritten on rhs).
   */
  inline void Gauss_Elimination(unsigned long block_i, ScalarType* rhs, bool transposed = false) const;

  /*!
   * \brief Inverse diagonal block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[out] invBlock - Inverse block.
   */
  inline void InverseDiagonalBlock(unsigned long block_i, ScalarType *invBlock, bool transposed = false) const;

  /*!
   * \brief Inverse diagonal block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[out] invBlock - Inverse block.
   */
  inline void InverseDiagonalBlock_ILUMatrix(unsigned long block_i, ScalarType *invBlock) const;

  /*!
   * \brief Copies the block (i, j) of the matrix-by-blocks structure in the internal variable *block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   */
  inline ScalarType *GetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j);

  /*!
   * \brief Set the value of a block in the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to set to A(i, j).
   */
  inline void SetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, ScalarType *val_block);

  /*!
   * \brief Set the transposed value of a block in the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to set to A(i, j).
   */
  inline void SetBlockTransposed_ILUMatrix(unsigned long block_i, unsigned long block_j, ScalarType *val_block);

  /*!
   * \brief Performs the product of i-th row of the upper part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the upper part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \param[in] col_ub - Exclusive upper bound for column indices considered in multiplication.
   * \param[out] prod - Result of the product U(A)*vec.
   */
  inline void UpperProduct(const CSysVector<ScalarType> & vec, unsigned long row_i,
                           unsigned long col_ub, ScalarType *prod) const;

  /*!
   * \brief Performs the product of i-th row of the lower part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the lower part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \param[in] col_lb - Inclusive lower bound for column indices considered in multiplication.
   * \param[out] prod - Result of the product L(A)*vec.
   */
  inline void LowerProduct(const CSysVector<ScalarType> & vec, unsigned long row_i,
                           unsigned long col_lb, ScalarType *prod) const;

  /*!
   * \brief Performs the product of i-th row of the diagonal part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the diagonal part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return prod Result of the product D(A)*vec (stored at *prod_row_vector).
   */
  inline void DiagonalProduct(const CSysVector<ScalarType> & vec, unsigned long row_i, ScalarType *prod) const;

  /*!
   * \brief Performs the product of i-th row of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the row of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return Result of the product (stored at *prod_row_vector).
   */
  void RowProduct(const CSysVector<ScalarType> & vec, unsigned long row_i, ScalarType *prod) const;

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
   * \brief Initializes sparse matrix system.
   * \param[in] npoint - Number of points including halos.
   * \param[in] npointdomain - Number of points excluding halos.
   * \param[in] nvar - Number of variables.
   * \param[in] neqn - Number of equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Initialize(unsigned long npoint, unsigned long npointdomain,
                  unsigned short nvar, unsigned short neqn,
                  bool EdgeConnect, CGeometry *geometry, CConfig *config);

  /*!
   * \brief Sets to zero all the entries of the sparse matrix.
   */
  void SetValZero(void);

  /*!
   * \brief Sets to zero all the block diagonal entries of the sparse matrix.
   */
  void SetValDiagonalZero(void);

  /*!
   * \brief Routine to load a vector quantity into the data structures for MPI point-to-point
   *        communication and to launch non-blocking sends and recvs.
   * \param[in] x        - CSysVector holding the array of data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   */
  template<class OtherType>
  void InitiateComms(const CSysVector<OtherType> & x,
                     CGeometry *geometry,
                     CConfig *config,
                     unsigned short commType) const;

  /*!
   * \brief Routine to complete the set of non-blocking communications launched by
   *        InitiateComms() and unpacking of the data in the vector.
   * \param[in] x        - CSysVector holding the array of data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be unpacked.
   */
  template<class OtherType>
  void CompleteComms(CSysVector<OtherType> & x,
                     CGeometry *geometry,
                     CConfig *config,
                     unsigned short commType) const;

  /*!
   * \brief Get a pointer to the start of block "ij"
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \return Pointer to location in memory where the block starts.
   */
  inline ScalarType *GetBlock(unsigned long block_i, unsigned long block_j) {

    for (unsigned long index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++)
      if (col_ind[index] == block_j)
        return &(matrix[index*nVar*nEqn]);
    return nullptr;
  }

  /*!
   * \brief Get a pointer to the start of block "ij", const version
   */
  inline const ScalarType *GetBlock(unsigned long block_i, unsigned long block_j) const {

    for (unsigned long index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++)
      if (col_ind[index] == block_j)
        return &(matrix[index*nVar*nEqn]);
    return nullptr;
  }

  /*!
   * \brief Gets the value of a particular entry in block "ij".
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] iVar - Row of the block.
   * \param[in] jVar - Column of the block.
   * \return Value of the block entry.
   */
  inline ScalarType GetBlock(unsigned long block_i, unsigned long block_j,
                             unsigned short iVar, unsigned short jVar) const {

    for (unsigned long index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++)
      if (col_ind[index] == block_j)
        return matrix[index*nVar*nEqn+iVar*nEqn+jVar];
    return 0.0;
  }

  /*!
   * \brief Set the value of a block in the sparse matrix.
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to set to A(i, j).
   */
  template<class OtherType>
  inline void SetBlock(unsigned long block_i, unsigned long block_j, OtherType **val_block) {

    unsigned long iVar, jVar, index;

    for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
      if (col_ind[index] == block_j) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nEqn; jVar++)
            matrix[index*nVar*nEqn+iVar*nEqn+jVar] = PassiveAssign<ScalarType,OtherType>(val_block[iVar][jVar]);
        break;
      }
    }
  }

  /*!
   * \brief Set the value of a block in the sparse matrix.
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to set to A(i, j).
   */
  template<class OtherType>
  inline void SetBlock(unsigned long block_i, unsigned long block_j, OtherType *val_block) {

    unsigned long iVar, index;

    for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
      if (col_ind[index] == block_j) {
        for (iVar = 0; iVar < nVar*nEqn; iVar++)
          matrix[index*nVar*nEqn+iVar] = PassiveAssign<ScalarType,OtherType>(val_block[iVar]);
        break;
      }
    }
  }

  /*!
   * \brief Adds the specified block to the sparse matrix.
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to add to A(i, j).
   */
  template<class OtherType>
  inline void AddBlock(unsigned long block_i, unsigned long block_j, OtherType **val_block) {

    unsigned long iVar, jVar, index;

    for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
      if (col_ind[index] == block_j) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nEqn; jVar++)
            matrix[index*nVar*nEqn+iVar*nEqn+jVar] += PassiveAssign<ScalarType,OtherType>(val_block[iVar][jVar]);
        break;
      }
    }
  }

  /*!
   * \brief Subtracts the specified block to the sparse matrix.
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to subtract to A(i, j).
   */
  template<class OtherType>
  inline void SubtractBlock(unsigned long block_i, unsigned long block_j, OtherType **val_block) {

    unsigned long iVar, jVar, index;

    for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
      if (col_ind[index] == block_j) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nEqn; jVar++)
            matrix[index*nVar*nEqn+iVar*nEqn+jVar] -= PassiveAssign<ScalarType,OtherType>(val_block[iVar][jVar]);
        break;
      }
    }
  }

  /*!
   * \brief Update 4 blocks ii, ij, ji, jj (add to i* sub from j*).
   * \note The template parameter Sign, can be used create a "subtractive"
   *       update i.e. subtract from row i and add to row j instead.
   * \param[in] edge - Index of edge that connects iPoint and jPoint.
   * \param[in] iPoint - Row to which we add the blocks.
   * \param[in] jPoint - Row from which we subtract the blocks.
   * \param[in] block_i - Adds to ii, subs from ji.
   * \param[in] block_j - Adds to ij, subs from jj.
   */
  template<class OtherType, int Sign = 1>
  inline void UpdateBlocks(unsigned long iEdge, unsigned long iPoint, unsigned long jPoint,
                           OtherType **block_i, OtherType **block_j) {

    ScalarType *bii = &matrix[dia_ptr[iPoint]*nVar*nEqn];
    ScalarType *bjj = &matrix[dia_ptr[jPoint]*nVar*nEqn];
    ScalarType *bij = &matrix[edge_ptr(iEdge,0)*nVar*nEqn];
    ScalarType *bji = &matrix[edge_ptr(iEdge,1)*nVar*nEqn];

    unsigned long iVar, jVar, offset = 0;

    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nEqn; jVar++) {
        bii[offset] += PassiveAssign<ScalarType,OtherType>(block_i[iVar][jVar]) * Sign;
        bij[offset] += PassiveAssign<ScalarType,OtherType>(block_j[iVar][jVar]) * Sign;
        bji[offset] -= PassiveAssign<ScalarType,OtherType>(block_i[iVar][jVar]) * Sign;
        bjj[offset] -= PassiveAssign<ScalarType,OtherType>(block_j[iVar][jVar]) * Sign;
        ++offset;
      }
    }
  }

  /*!
   * \brief Adds the specified value to the diagonal of the (i, i) subblock
   *        of the matrix-by-blocks structure.
   * \param[in] block_i - Diagonal index.
   * \param[in] val_matrix - Value to add to the diagonal elements of A(i, i).
   */
  template<class OtherType>
  inline void AddVal2Diag(unsigned long block_i, OtherType val_matrix) {

    unsigned long iVar, index;

    for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
      if (col_ind[index] == block_i) { // Only elements on the diagonal
        for (iVar = 0; iVar < nVar; iVar++)
          matrix[index*nVar*nVar+iVar*nVar+iVar] += PassiveAssign<ScalarType,OtherType>(val_matrix);
        break;
      }
    }
  }

  /*!
   * \brief Sets the specified value to the diagonal of the (i, i) subblock
   *        of the matrix-by-blocks structure.
   * \param[in] block_i - Diagonal index.
   * \param[in] val_matrix - Value to add to the diagonal elements of A(i, i).
   */
  template<class OtherType>
  inline void SetVal2Diag(unsigned long block_i, OtherType val_matrix) {

    unsigned long iVar, jVar, index;

    for (index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++) {
      if (col_ind[index] == block_i) { // Only elements on the diagonal

        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            matrix[index*nVar*nVar+iVar*nVar+jVar] = 0.0;

        for (iVar = 0; iVar < nVar; iVar++)
          matrix[index*nVar*nVar+iVar*nVar+iVar] = PassiveAssign<ScalarType,OtherType>(val_matrix);

        break;
      }
    }
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
  template<class OtherType>
  void EnforceSolutionAtNode(const unsigned long node_i, const OtherType *x_i, CSysVector<OtherType> & b);

  /*!
   * \brief Add a scaled sparse matrix to "this" (axpy-type operation, A = A+alpha*B).
   * \note Matrices must have the same sparse pattern.
   * \param[in] alpha - The scaling constant.
   * \param[in] B - Matrix being.
   */
  template<class OtherType>
  void MatrixMatrixAddition(OtherType alpha, const CSysMatrix<OtherType>& B);

  /*!
   * \brief Performs the product of a sparse matrix by a CSysVector.
   * \param[in] vec - CSysVector to be multiplied by the sparse matrix A.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] prod - Result of the product.
   */
  void MatrixVectorProduct(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                           CGeometry *geometry, CConfig *config) const;

  /*!
   * \brief Performs the product of a sparse matrix by a CSysVector.
   * \param[in] vec - CSysVector to be multiplied by the sparse matrix A.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] prod - Result of the product.
   */
  void MatrixVectorProductTransposed(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                     CGeometry *geometry, CConfig *config) const;

  /*!
   * \brief Build the Jacobi preconditioner.
   */
  void BuildJacobiPreconditioner(bool transpose = false);

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeJacobiPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                   CGeometry *geometry, CConfig *config) const;

  /*!
   * \brief Build the ILU preconditioner.
   * \param[in] transposed - Flag to use the transposed matrix to construct the preconditioner.
   */
  void BuildILUPreconditioner(bool transposed = false);

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeILUPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                CGeometry *geometry, CConfig *config) const;

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   */
  void ComputeLU_SGSPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                   CGeometry *geometry, CConfig *config) const;

  /*!
   * \brief Build the Linelet preconditioner.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \return Average number of points per linelet.
   */
  unsigned long BuildLineletPreconditioner(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   */
  void ComputeLineletPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                    CGeometry *geometry, CConfig *config) const;

  /*!
   * \brief Compute the linear residual.
   * \param[in] sol - Solution (x).
   * \param[in] f - Right hand side (b).
   * \param[out] res - Residual (Ax-b).
   */
  void ComputeResidual(const CSysVector<ScalarType> & sol, const CSysVector<ScalarType> & f,
                       CSysVector<ScalarType> & res) const;

  /*!
   * \brief Factorize matrix using PaStiX.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] kind_fact - Type of factorization.
   * \param[in] transposed - Flag to use the transposed matrix during application of the preconditioner.
   */
  void BuildPastixPreconditioner(CGeometry *geometry, CConfig *config, unsigned short kind_fact, bool transposed = false);

  /*!
   * \brief Apply the PaStiX factorization to CSysVec.
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product M*vec.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputePastixPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod,
                                   CGeometry *geometry, CConfig *config) const;

};

#ifdef CODI_REVERSE_TYPE
template<> template<>
inline passivedouble CSysMatrix<passivedouble>::ActiveAssign(const su2double & val) const { return SU2_TYPE::GetValue(val); }

template<> template<>
inline passivedouble CSysMatrix<su2double>::ActiveAssign(const su2double & val) const { return SU2_TYPE::GetValue(val); }
#endif
