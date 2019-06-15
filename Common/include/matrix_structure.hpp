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

#include "./mpi_structure.hpp"
#include <limits>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "vector_structure.hpp"

#ifdef HAVE_MKL
#include "mkl.h"
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
  int rank, 	/*!< \brief MPI Rank. */
  size;       	/*!< \brief MPI Size. */
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
  ScalarType *block_inverse;             /*!< \brief Internal array to store a subblock of the matrix. */
  ScalarType *block_weight;             /*!< \brief Internal array to store a subblock of the matrix. */
  ScalarType *prod_block_vector; /*!< \brief Internal array to store the product of a subblock with a vector. */
  ScalarType *prod_row_vector;   /*!< \brief Internal array to store the product of a matrix-by-blocks "row" with a vector. */
  ScalarType *aux_vector;         /*!< \brief Auxiliary array to store intermediate results. */
  ScalarType *sum_vector;         /*!< \brief Auxiliary array to store intermediate results. */
  ScalarType *invM;              /*!< \brief Inverse of (Jacobi) preconditioner. */
  
  bool *LineletBool;                          /*!< \brief Identify if a point belong to a linelet. */
  vector<unsigned long> *LineletPoint;        /*!< \brief Linelet structure. */
  unsigned long nLinelet;                     /*!< \brief Number of Linelets in the system. */
  ScalarType **UBlock, **invUBlock, **LBlock,
  **yVector, **zVector, **rVector, *LFBlock,
  *LyVector, *FzVector;           /*!< \brief Arrays of the Linelet preconditioner methodology. */
  unsigned long max_nElem;

#if defined(HAVE_MKL) && !(defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))
  void * MatrixMatrixProductJitter;                   		/*!< \brief Jitter handle for MKL JIT based GEMM. */
  dgemm_jit_kernel_t MatrixMatrixProductKernel;               	/*!< \brief MKL JIT based GEMM kernel. */
  void * MatrixVectorProductJitterBetaZero;           		/*!< \brief Jitter handle for MKL JIT based GEMV. */
  dgemm_jit_kernel_t MatrixVectorProductKernelBetaZero;       	/*!< \brief MKL JIT based GEMV kernel. */
  void * MatrixVectorProductJitterBetaOne;            		/*!< \brief Jitter handle for MKL JIT based GEMV with BETA=1.0. */
  dgemm_jit_kernel_t MatrixVectorProductKernelBetaOne;        	/*!< \brief MKL JIT based GEMV kernel with BETA=1.0. */
  bool useMKL;
#endif

  /*!
   * \brief Handle type conversion for when we Set, Add, etc. blocks, preserving derivative information (if supported by types).
   */
  template<class DstType, class SrcType>
  DstType ActiveAssign(const SrcType & val) const;

  /*!
   * \brief Handle type conversion for when we Set, Add, etc. blocks, discarding derivative information.
   */
  template<class DstType, class SrcType>
  DstType PassiveAssign(const SrcType & val) const;

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
   * \brief Initializes space matrix system.
   * \param[in] nVar - Number of variables.
   * \param[in] nEqn - Number of equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Initialize(unsigned long nPoint, unsigned long nPointDomain, unsigned short nVar, unsigned short nEqn,
                  bool EdgeConnect, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Assigns values to the sparse-matrix structure.
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
   * \brief Assigns values to the sparse-matrix structure.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] iPoint - Base point to compute neighbours.
   * \param[in] deep_level - Deep level for the recursive algorithm.
   * \param[in] fill_level - ILU fill in level.
   * \param[in] EdgeConnect - There is (or not) an edge structure).
   * \param[in] vneighs - Storage the neighbours points to iPoint.
   */
  void SetNeighbours(CGeometry *geometry, unsigned long iPoint, unsigned short deep_level, unsigned short fill_level, bool EdgeConnect, vector<unsigned long> & vneighs);
  
  /*!
   * \brief Sets to zero all the entries of the sparse matrix.
   */
  void SetValZero(void);
  
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
   * \brief Copies the block (i, j) of the matrix-by-blocks structure in the internal variable *block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   */
  ScalarType *GetBlock(unsigned long block_i, unsigned long block_j);
  
  /*!
   * \brief Copies the block (i, j) of the matrix-by-blocks structure in the internal variable *block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   */
  ScalarType GetBlock(unsigned long block_i, unsigned long block_j, unsigned short iVar, unsigned short jVar);
  
  /*!
   * \brief Set the value of a block in the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to set to A(i, j).
   */
  template<class OtherType>
  void SetBlock(unsigned long block_i, unsigned long block_j, OtherType **val_block);
  
  /*!
   * \brief Set the value of a block in the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to set to A(i, j).
   */
  template<class OtherType>
  void SetBlock(unsigned long block_i, unsigned long block_j, OtherType *val_block);
  
  /*!
   * \brief Adds the specified block to the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to add to A(i, j).
   */
  template<class OtherType>
  void AddBlock(unsigned long block_i, unsigned long block_j, OtherType **val_block);
  
  /*!
   * \brief Subtracts the specified block to the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to subtract to A(i, j).
   */
  template<class OtherType>
  void SubtractBlock(unsigned long block_i, unsigned long block_j, OtherType **val_block);
  
  /*!
   * \brief Copies the block (i, j) of the matrix-by-blocks structure in the internal variable *block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   */
  ScalarType *GetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j);
  
  /*!
   * \brief Set the value of a block in the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to set to A(i, j).
   */
  void SetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, ScalarType *val_block);
  
  
  /*!
   * \brief Set the transposed value of a block in the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to set to A(i, j).
   */
  void SetBlockTransposed_ILUMatrix(unsigned long block_i, unsigned long block_j, ScalarType *val_block);
  
  /*!
   * \brief Subtracts the specified block to the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to subtract to A(i, j).
   */
  void SubtractBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, ScalarType *val_block);
  
  /*!
   * \brief Adds the specified value to the diagonal of the (i, i) subblock
   *        of the matrix-by-blocks structure.
   * \param[in] block_i - Index of the block in the matrix-by-blocks structure.
   * \param[in] val_matrix - Value to add to the diagonal elements of A(i, i).
   */
  template<class OtherType>
  void AddVal2Diag(unsigned long block_i, OtherType val_matrix);
  
  /*!
   * \brief Sets the specified value to the diagonal of the (i, i) subblock
   *        of the matrix-by-blocks structure.
   * \param[in] block_i - Index of the block in the matrix-by-blocks structure.
   * \param[in] val_matrix - Value to add to the diagonal elements of A(i, i).
   */
  template<class OtherType>
  void SetVal2Diag(unsigned long block_i, OtherType val_matrix);
  
  /*!
   * \brief Calculates the matrix-vector product
   * \param[in] matrix
   * \param[in] vector
   * \param[out] product
   */
  void MatrixVectorProduct(ScalarType *matrix, ScalarType *vector, ScalarType *product);
  
  /*!
   * \brief Calculates the matrix-matrix product
   * \param[in] matrix_a
   * \param[in] matrix_b
   * \param[out] product
   */
  void MatrixMatrixProduct(ScalarType *matrix_a, ScalarType *matrix_b, ScalarType *product);
  
  /*!
   * \brief Deletes the values of the row i of the sparse matrix.
   * \param[in] i - Index of the row.
   */
  void DeleteValsRowi(unsigned long i);

  /*!
   * \brief Recursive definition of determinate using expansion by minors. Written by Paul Bourke
   * \param[in] a - Matrix to compute the determinant.
   * \param[in] n - Size of the quare matrix.
   * \return Value of the determinant.
   */
  ScalarType MatrixDeterminant(ScalarType **a, unsigned long n);
  
  /*!
   * \brief Find the cofactor matrix of a square matrix. Written by Paul Bourke
   * \param[in] a - Matrix to compute the determinant.
   * \param[in] n - Size of the quare matrix.
   * \param[out] b - cofactor matrix
   */
  void MatrixCoFactor(ScalarType **a, unsigned long n, ScalarType **b) ;
  
  /*!
   * \brief Transpose of a square matrix, do it in place. Written by Paul Bourke
   * \param[in] a - Matrix to compute the determinant.
   * \param[in] n - Size of the quare matrix.
   */
  void MatrixTranspose(ScalarType **a, unsigned long n) ;
  
  /*!
   * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i, i) subblock and rhs.
   * \param[in] block_i - Index of the (i, i) subblock in the matrix-by-blocks structure.
   * \param[in] rhs - Right-hand-side of the linear system.
   * \param[in] transposed - If true the transposed of the block is used (default = false).
   * \return Solution of the linear system (overwritten on rhs).
   */
  void Gauss_Elimination(unsigned long block_i, ScalarType* rhs, bool transposed = false);
  
  /*!
   * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i, i) subblock and rhs.
   * \param[in] Block - matrix-by-blocks structure.
   * \param[in] rhs - Right-hand-side of the linear system.
   * \return Solution of the linear system (overwritten on rhs).
   */
  void Gauss_Elimination(ScalarType* Block, ScalarType* rhs);
  
  /*!
   * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i, i) subblock and rhs.
   * \param[in] block_i - Index of the (i, i) subblock in the matrix-by-blocks structure.
   * \param[in] rhs - Right-hand-side of the linear system.
   * \return Solution of the linear system (overwritten on rhs).
   */
  void Gauss_Elimination_ILUMatrix(unsigned long block_i, ScalarType* rhs);
  
  /*!
   * \fn void CSysMatrix::ProdBlockVector(unsigned long block_i, unsigned long block_j, su2double* vec);
   * \brief Performs the product of the block (i, j) by vector vec.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] vec - Vector to be multiplied by the block (i, j) of the sparse matrix A.
   * \return Product of A(i, j) by vector *vec (stored at *prod_block_vector).
   */
  void ProdBlockVector(unsigned long block_i, unsigned long block_j, const CSysVector<ScalarType> & vec);
  
  /*!
   * \brief Performs the product of i-th row of the upper part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the upper part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return prod Result of the product U(A)*vec (stored at *prod_row_vector).
   */
  void UpperProduct(CSysVector<ScalarType> & vec, unsigned long row_i);
  
  /*!
   * \brief Performs the product of i-th row of the lower part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the lower part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return prod Result of the product L(A)*vec (stored at *prod_row_vector).
   */
  void LowerProduct(CSysVector<ScalarType> & vec, unsigned long row_i);
  
  /*!
   * \brief Performs the product of i-th row of the diagonal part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the diagonal part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return prod Result of the product D(A)*vec (stored at *prod_row_vector).
   */
  void DiagonalProduct(CSysVector<ScalarType> & vec, unsigned long row_i);
  
  /*!
   * \brief Performs the product of i-th row of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the row of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return Result of the product (stored at *prod_row_vector).
   */
  void RowProduct(const CSysVector<ScalarType> & vec, unsigned long row_i);
  
  /*!
   * \brief Performs the product of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the sparse matrix A.
   * \param[out] prod - Result of the product.
   * \return Result of the product A*vec.
   */
  void MatrixVectorProduct(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod);
  
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
   * \brief Performs the product of two block matrices.
   */
  void GetMultBlockBlock(ScalarType *c, ScalarType *a, ScalarType *b);
  
  /*!
   * \brief Performs the product of a block matrices by a vector.
   */
  void GetMultBlockVector(ScalarType *c, ScalarType *a, ScalarType *b);
  
  /*!
   * \brief Performs the subtraction of two matrices.
   */
  void GetSubsBlock(ScalarType *c, ScalarType *a, ScalarType *b);
  
  /*!
   * \brief Performs the subtraction of two vectors.
   */
  void GetSubsVector(ScalarType *c, ScalarType *a, ScalarType *b);
  
  /*!
   * \brief Inverse diagonal block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[out] invBlock - Inverse block.
   */
  void InverseDiagonalBlock(unsigned long block_i, ScalarType *invBlock, bool transpose = false);
  
 	/*!
   * \brief Inverse diagonal block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[out] invBlock - Inverse block.
   */
  void InverseDiagonalBlock_ILUMatrix(unsigned long block_i, ScalarType *invBlock);
  
  /*!
   * \brief Inverse a block.
   * \param[in] Block - block matrix.
   * \param[out] invBlock - Inverse block.
   */
  void InverseBlock(ScalarType *Block, ScalarType *invBlock);
  
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
   * \brief Apply Jacobi as a classical iterative smoother
   * \param[in] b - CSysVector containing the residual (b)
   * \param[in] x - CSysVector containing the solution (x^k)
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] residual
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] x - CSysVector containing the result of the smoothing (x^k+1 = x^k + M^-1*(b - A*x^k).
   */
  unsigned long Jacobi_Smoother(const CSysVector<ScalarType> & b, CSysVector<ScalarType> & x, CMatrixVectorProduct<ScalarType> & mat_vec, ScalarType tol, unsigned long m, ScalarType *residual, bool monitoring, CGeometry *geometry, CConfig *config);
  
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
   * \brief Apply ILU as a classical iterative smoother
   * \param[in] b - CSysVector containing the residual (b)
   * \param[in] x - CSysVector containing the solution (x^k)
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] residual
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] x - CSysVector containing the result of the smoothing (x^k+1 = x^k + M^-1*(b - A*x^k).
   */
  unsigned long ILU_Smoother(const CSysVector<ScalarType> & b, CSysVector<ScalarType> & x, CMatrixVectorProduct<ScalarType> & mat_vec, ScalarType tol, unsigned long m, ScalarType *residual, bool monitoring, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   */
  void ComputeLU_SGSPreconditioner(const CSysVector<ScalarType> & vec, CSysVector<ScalarType> & prod, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Apply LU_SGS as a classical iterative smoother
   * \param[in] b - CSysVector containing the residual (b)
   * \param[in] x - CSysVector containing the solution (x^k)
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] residual
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] x - CSysVector containing the result of the smoothing (x^k+1 = x^k + M^-1*(b - A*x^k).
   */
  unsigned long LU_SGS_Smoother(const CSysVector<ScalarType> & b, CSysVector<ScalarType> & x, CMatrixVectorProduct<ScalarType> & mat_vec, ScalarType tol, unsigned long m, ScalarType *residual, bool monitoring, CGeometry *geometry, CConfig *config);
  
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
  
};

/*!
 * \class CSysMatrixVectorProduct
 * \brief specialization of matrix-vector product that uses CSysMatrix class
 */
template<class ScalarType>
class CSysMatrixVectorProduct : public CMatrixVectorProduct<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the product. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the products
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  CSysMatrixVectorProduct(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~CSysMatrixVectorProduct() {}
  
  /*!
   * \brief operator that defines the CSysMatrix-CSysVector product
   * \param[in] u - CSysVector that is being multiplied by the sparse matrix
   * \param[out] v - CSysVector that is the result of the product
   */
  void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const;
};

/*!
 * \class CSysMatrixVectorProduct
 * \brief specialization of matrix-vector product that uses CSysMatrix class
 */
template<class ScalarType>
class CSysMatrixVectorProductTransposed : public CMatrixVectorProduct<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the product. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the products
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  CSysMatrixVectorProductTransposed(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~CSysMatrixVectorProductTransposed() {}
  
  /*!
   * \brief operator that defines the CSysMatrix-CSysVector product
   * \param[in] u - CSysVector that is being multiplied by the sparse matrix
   * \param[out] v - CSysVector that is the result of the product
   */
  void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const;
};

/*!
 * \class CJacobiPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class ScalarType>
class CJacobiPreconditioner : public CPreconditioner<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  CJacobiPreconditioner(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~CJacobiPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - CSysVector that is being preconditioned
   * \param[out] v - CSysVector that is the result of the preconditioning
   */
  void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const;
};

/*!
 * \class CJacobiTransposedPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class ScalarType>
class CJacobiTransposedPreconditioner : public CPreconditioner<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  CJacobiTransposedPreconditioner(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~CJacobiTransposedPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - CSysVector that is being preconditioned
   * \param[out] v - CSysVector that is the result of the preconditioning
   */
  void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const;
};

/*!
 * \class CILUPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class ScalarType>
class CILUPreconditioner : public CPreconditioner<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  CILUPreconditioner(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~CILUPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - CSysVector that is being preconditioned
   * \param[out] v - CSysVector that is the result of the preconditioning
   */
  void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const;
};

/*!
 * \class CLU_SGSPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class ScalarType>
class CLU_SGSPreconditioner : public CPreconditioner<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  CLU_SGSPreconditioner(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~CLU_SGSPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - CSysVector that is being preconditioned
   * \param[out] v - CSysVector that is the result of the preconditioning
   */
  void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const;
};

/*!
 * \class CLineletPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class ScalarType>
class CLineletPreconditioner : public CPreconditioner<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  CLineletPreconditioner(CSysMatrix<ScalarType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~CLineletPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - CSysVector that is being preconditioned
   * \param[out] v - CSysVector that is the result of the preconditioning
   */
  void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const;
};

#include "matrix_structure.inl"
