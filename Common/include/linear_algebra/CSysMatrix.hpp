/*!
 * \file matrix_structure.hpp
 * \brief Headers of the main subroutines for creating the sparse matrices-by-blocks.
 *        The subroutines and functions are in the <i>matrix_structure.cpp</i> file.
 * \author F. Palacios, A. Bueno, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../mpi_structure.hpp"
#include <limits>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "../config_structure.hpp"
#include "../geometry_structure.hpp"
#include "CSysVector.hpp"
#include "CPastixWrapper.hpp"

#if defined(HAVE_MKL) && !defined(CODI_FORWARD_TYPE)
#include "mkl.h"
#ifndef __INTEL_MKL__
  #error Could not determine the MKL version
#endif
/*--- JIT is only available since 2019 ---*/
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

using namespace std;

const su2double eps = numeric_limits<passivedouble>::epsilon(); /*!< \brief machine epsilon */


/*!
 * \class CSysMatrix
 * \brief Main class for defining sparse matrices-by-blocks
 with compressed row format.
 * \author A. Bueno, F. Palacios
 */
template<class ScalarType>
class CSysMatrix {
private:
  int rank;     /*!< \brief MPI Rank. */
  int size;     /*!< \brief MPI Size. */
  unsigned long nPoint,   /*!< \brief Number of points in the grid. */
  nPointDomain,           /*!< \brief Number of points in the grid. */
  nVar,                   /*!< \brief Number of variables. */
  nEqn;                   /*!< \brief Number of equations. */
  ScalarType *matrix;            /*!< \brief Entries of the sparse matrix. */
  ScalarType *ILU_matrix;         /*!< \brief Entries of the ILU sparse matrix. */
  unsigned long nnz;                 /*!< \brief Number of possible nonzero entries in the matrix. */
  unsigned long *row_ptr;            /*!< \brief Pointers to the first element in each row. */
  unsigned long *col_ind;            /*!< \brief Column index for each of the elements in val(). */
  unsigned long nnz_ilu;             /*!< \brief Number of possible nonzero entries in the matrix (ILU). */
  unsigned long *row_ptr_ilu;        /*!< \brief Pointers to the first element in each row (ILU). */
  unsigned long *col_ind_ilu;        /*!< \brief Column index for each of the elements in val() (ILU). */
  unsigned short ilu_fill_in;        /*!< \brief Fill in level for the ILU preconditioner. */

  ScalarType *block;             /*!< \brief Internal array to store a subblock of the matrix. */
  ScalarType *block_inverse;     /*!< \brief Internal array to store a subblock of the matrix. */
  ScalarType *block_weight;      /*!< \brief Internal array to store a subblock of the matrix. */
  ScalarType *prod_row_vector;   /*!< \brief Internal array to store the product of a matrix-by-blocks "row" with a vector. */
  ScalarType *aux_vector;        /*!< \brief Auxiliary array to store intermediate results. */
  ScalarType *sum_vector;        /*!< \brief Auxiliary array to store intermediate results. */
  ScalarType *invM;              /*!< \brief Inverse of (Jacobi) preconditioner, or diagonal of ILU. */

  unsigned long nLinelet;                       /*!< \brief Number of Linelets in the system. */
  vector<bool> LineletBool;                     /*!< \brief Identify if a point belong to a Linelet. */
  vector<vector<unsigned long> > LineletPoint;  /*!< \brief Linelet structure. */
  vector<const ScalarType*> LineletUpper;       /*!< \brief Pointers to the upper blocks of the tri-diag system. */
  vector<ScalarType> LineletInvDiag;            /*!< \brief Inverse of the diagonal blocks of the tri-diag system. */
  vector<ScalarType> LineletVector;             /*!< \brief Solution and RHS of the tri-diag system. */

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
  lapack_int * mkl_ipiv;
#endif

#ifdef HAVE_PASTIX
  CPastixWrapper pastix_wrapper;
#endif

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
   * \brief Assigns values to the sparse-matrix structure (used in Initialize).
   * \param[in] val_nPoint - Number of points in the nPoint x nPoint block structure
   * \param[in] val_nVar - Number of nVar x nVar variables in each subblock of the matrix-by-block structure.
   * \param[in] val_nEq - Number of nEqn x nVar variables in each subblock of the matrix-by-block structure.
   * \param[in] val_row_ptr - Pointers to the first element in each row.
   * \param[in] val_col_ind - Column index for each of the elements in val().
   * \param[in] val_nnz - Number of possible nonzero entries in the matrix.
   * \param[in] config - Definition of the particular problem.
   */
  void SetIndexes(unsigned long val_nPoint, unsigned long val_nPointDomain, unsigned short val_nVar, unsigned short val_nEq, unsigned long* val_row_ptr, unsigned long* val_col_ind, unsigned long val_nnz, CConfig *config);

  /*!
   * \brief Assigns values to the sparse-matrix structure (used in Initialize).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] iPoint - Base point to compute neighbours.
   * \param[in] deep_level - Deep level for the recursive algorithm.
   * \param[in] fill_level - ILU fill in level.
   * \param[in] EdgeConnect - There is (or not) an edge structure).
   * \param[in] vneighs - Storage the neighbours points to iPoint.
   */
  void SetNeighbours(CGeometry *geometry, unsigned long iPoint, unsigned short deep_level, unsigned short fill_level, bool EdgeConnect, vector<unsigned long> & vneighs);

  /*!
   * \brief Calculates the matrix-vector product: product = matrix*vector
   * \param[in] matrix
   * \param[in] vector
   * \param[out] product
   */
  inline void MatrixVectorProduct(const ScalarType *matrix, const ScalarType *vector, ScalarType *product);

  /*!
   * \brief Calculates the matrix-vector product: product += matrix*vector
   * \param[in] matrix
   * \param[in] vector
   * \param[in,out] product
   */
  inline void MatrixVectorProductAdd(const ScalarType *matrix, const ScalarType *vector, ScalarType *product);

  /*!
   * \brief Calculates the matrix-vector product: product -= matrix*vector
   * \param[in] matrix
   * \param[in] vector
   * \param[in,out] product
   */
  inline void MatrixVectorProductSub(const ScalarType *matrix, const ScalarType *vector, ScalarType *product);

  /*!
   * \brief Calculates the matrix-vector product: product += matrix^T * vector
   * \param[in] matrix
   * \param[in] vector
   * \param[in,out] product
   */
  inline void MatrixVectorProductTransp(const ScalarType *matrix, const ScalarType *vector, ScalarType *product);

  /*!
   * \brief Calculates the matrix-matrix product
   * \param[in] matrix_a
   * \param[in] matrix_b
   * \param[out] product
   */
  inline void MatrixMatrixProduct(const ScalarType *matrix_a, const ScalarType *matrix_b, ScalarType *product);

  /*!
   * \brief Subtract b from a and store the result in c.
   */
  inline void VectorSubtraction(const ScalarType *a, const ScalarType *b, ScalarType *c) {
    for(unsigned long iVar = 0; iVar < nVar; iVar++)
      c[iVar] = a[iVar] - b[iVar];
  }

  /*!
   * \brief Subtract b from a and store the result in c.
   */
  inline void MatrixSubtraction(const ScalarType *a, const ScalarType *b, ScalarType *c) {
    for(unsigned long iVar = 0; iVar < nVar*nEqn; iVar++)
      c[iVar] = a[iVar] - b[iVar];
  }

  /*!
   * \brief Solve a small (nVar x nVar) linear system using Gaussian elimination.
   * \param[in,out] matrix - On entry the system matrix, on exit the factorized matrix.
   * \param[in,out] vec - On entry the rhs, on exit the solution.
   */
  inline void Gauss_Elimination(ScalarType* matrix, ScalarType* vec);

  /*!
   * \brief Invert a small dense matrix.
   * \param[in] matrix - the matrix.
   * \param[out] inverse - the matrix inverse.
   */
  inline void MatrixInverse(const ScalarType *matrix, ScalarType *inverse);

  /*!
   * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i, i) subblock and rhs.
   * \param[in] block_i - Index of the (i, i) subblock in the matrix-by-blocks structure.
   * \param[in] rhs - Right-hand-side of the linear system.
   * \param[in] transposed - If true the transposed of the block is used (default = false).
   * \return Solution of the linear system (overwritten on rhs).
   */
  inline void Gauss_Elimination(unsigned long block_i, ScalarType* rhs, bool transposed = false);

  /*!
   * \brief Inverse diagonal block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[out] invBlock - Inverse block.
   */
  inline void InverseDiagonalBlock(unsigned long block_i, ScalarType *invBlock, bool transpose = false);

  /*!
   * \brief Inverse diagonal block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[out] invBlock - Inverse block.
   */
  inline void InverseDiagonalBlock_ILUMatrix(unsigned long block_i, ScalarType *invBlock);

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
   * \brief Subtracts the specified block to the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to subtract to A(i, j).
   */
  inline void SubtractBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, ScalarType *val_block);

  /*!
   * \brief Performs the product of i-th row of the upper part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the upper part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return prod Result of the product U(A)*vec (stored at *prod_row_vector).
   */
  void UpperProduct(const CSysVector<ScalarType> & vec, unsigned long row_i);

  /*!
   * \brief Performs the product of i-th row of the lower part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the lower part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return prod Result of the product L(A)*vec (stored at *prod_row_vector).
   */
  void LowerProduct(const CSysVector<ScalarType> & vec, unsigned long row_i);

  /*!
   * \brief Performs the product of i-th row of the diagonal part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the diagonal part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return prod Result of the product D(A)*vec (stored at *prod_row_vector).
   */
  void DiagonalProduct(const CSysVector<ScalarType> & vec, unsigned long row_i);

  /*!
   * \brief Performs the product of i-th row of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the row of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return Result of the product (stored at *prod_row_vector).
   */
  void RowProduct(const CSysVector<ScalarType> & vec, unsigned long row_i);

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
   * \param[in] nVar - Number of variables.
   * \param[in] nEqn - Number of equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Initialize(unsigned long nPoint, unsigned long nPointDomain, unsigned short nVar, unsigned short nEqn,
                  bool EdgeConnect, CGeometry *geometry, CConfig *config);

  /*!
   * \brief Sets to zero all the entries of the sparse matrix.
   */
  inline void SetValZero(void) {
    if(matrix != NULL)
      for (unsigned long index = 0; index < nnz*nVar*nEqn; index++)
        matrix[index] = 0.0;
  }

  /*!
   * \brief Routine to load a vector quantity into the data structures for MPI point-to-point communication and to launch non-blocking sends and recvs.
   * \param[in] x        - CSysVector holding the array of data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   */
  template<class OtherType>
  void InitiateComms(CSysVector<OtherType> & x,
                     CGeometry *geometry,
                     CConfig *config,
                     unsigned short commType);

  /*!
   * \brief Routine to complete the set of non-blocking communications launched by InitiateComms() and unpacking of the data in the vector.
   * \param[in] x        - CSysVector holding the array of data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be unpacked.
   */
  template<class OtherType>
  void CompleteComms(CSysVector<OtherType> & x,
                     CGeometry *geometry,
                     CConfig *config,
                     unsigned short commType);

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

    return NULL;
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
                             unsigned short iVar, unsigned short jVar) {

    for (unsigned long index = row_ptr[block_i]; index < row_ptr[block_i+1]; index++)
      if (col_ind[index] == block_j)
        return matrix[index*nVar*nEqn+iVar*nEqn+jVar];

    return 0.0;
  }

  /*!
   * \brief Set the value of a block in the sparse matrix.
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] **val_block - Block to set to A(i, j).
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
   * \param[in] *val_block - Block to set to A(i, j).
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
   * \param[in] **val_block - Block to add to A(i, j).
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
   * \param[in] **val_block - Block to subtract to A(i, j).
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
   * \brief Performs the product of a sparse matrix by a CSysVector.
   * \param[in] vec - CSysVector to be multiplied by the sparse matrix A.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] prod - Result of the product.
   */
  void MatrixVectorProduct(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod, CGeometry *geometry, CConfig *config);

  /*!
   * \brief Performs the product of a sparse matrix by a CSysVector.
   * \param[in] vec - CSysVector to be multiplied by the sparse matrix A.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] prod - Result of the product.
   */
  void MatrixVectorProductTransposed(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod, CGeometry *geometry, CConfig *config);

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
  void ComputeJacobiPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod, CGeometry *geometry, CConfig *config);

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
  void ComputeILUPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod, CGeometry *geometry, CConfig *config);

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   */
  void ComputeLU_SGSPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod, CGeometry *geometry, CConfig *config);

  /*!
   * \brief Build the Linelet preconditioner.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  unsigned short BuildLineletPreconditioner(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   */
  void ComputeLineletPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod, CGeometry *geometry, CConfig *config);

  /*!
   * \brief Compute the residual Ax-b
   * \param[in] sol - CSysVector to be multiplied by the preconditioner.
   * \param[in] f - Result of the product A*vec.
   * \param[out] res - Result of the product A*vec.
   */
  void ComputeResidual(const CSysVector<ScalarType> & sol, const CSysVector<ScalarType> & f, CSysVector<ScalarType> & res);

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
  void ComputePastixPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod, CGeometry *geometry, CConfig *config);

};

#ifdef CODI_REVERSE_TYPE
template<> template<>
inline passivedouble CSysMatrix<passivedouble>::ActiveAssign(const su2double & val) const { return SU2_TYPE::GetValue(val); }

template<> template<>
inline passivedouble CSysMatrix<su2double>::ActiveAssign(const su2double & val) const { return SU2_TYPE::GetValue(val); }
#endif
