/*!
 * \file matrix_structure.hpp
 * \brief Headers of the main subroutines for creating the sparse matrices-by-blocks.
 *        The subroutines and functions are in the <i>matrix_structure.cpp</i> file.
 * \author F. Palacios, A. Bueno, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

using namespace std;

const su2double eps = numeric_limits<passivedouble>::epsilon(); /*!< \brief machine epsilon */


/*!
 * \class CSysMatrix
 * \brief Main class for defining sparse matrices-by-blocks
 with compressed row format.
 * \author A. Bueno, F. Palacios
 * \version 5.0.0 "Raven"
 */
template<class CalcType>
class TCSysMatrix {
private:
  int rank, 	/*!< \brief MPI Rank. */
  size;       	/*!< \brief MPI Size. */
  unsigned long nPoint,   /*!< \brief Number of points in the grid. */
  nPointDomain,           /*!< \brief Number of points in the grid. */
  nVar,                   /*!< \brief Number of variables. */
  nEqn;                   /*!< \brief Number of equations. */
  CalcType *matrix;            /*!< \brief Entries of the sparse matrix. */
  CalcType *ILU_matrix;         /*!< \brief Entries of the ILU sparse matrix. */
  unsigned long nnz;                 /*!< \brief Number of possible nonzero entries in the matrix. */
  unsigned long *row_ptr;            /*!< \brief Pointers to the first element in each row. */
  unsigned long *col_ind;            /*!< \brief Column index for each of the elements in val(). */
  unsigned long nnz_ilu;             /*!< \brief Number of possible nonzero entries in the matrix (ILU). */
  unsigned long *row_ptr_ilu;        /*!< \brief Pointers to the first element in each row (ILU). */
  unsigned long *col_ind_ilu;        /*!< \brief Column index for each of the elements in val() (ILU). */
  unsigned short ilu_fill_in;        /*!< \brief Fill in level for the ILU preconditioner. */
  
  CalcType *block;             /*!< \brief Internal array to store a subblock of the matrix. */
  CalcType *block_inverse;             /*!< \brief Internal array to store a subblock of the matrix. */
  CalcType *block_weight;             /*!< \brief Internal array to store a subblock of the matrix. */
  CalcType *prod_block_vector; /*!< \brief Internal array to store the product of a subblock with a vector. */
  CalcType *prod_row_vector;   /*!< \brief Internal array to store the product of a matrix-by-blocks "row" with a vector. */
  CalcType *aux_vector;         /*!< \brief Auxiliary array to store intermediate results. */
  CalcType *sum_vector;         /*!< \brief Auxiliary array to store intermediate results. */
  CalcType *invM;              /*!< \brief Inverse of (Jacobi) preconditioner. */
  
  bool *LineletBool;                          /*!< \brief Identify if a point belong to a linelet. */
  vector<unsigned long> *LineletPoint;        /*!< \brief Linelet structure. */
  unsigned long nLinelet;                     /*!< \brief Number of Linelets in the system. */
  CalcType **UBlock, **invUBlock, **LBlock,
  **yVector, **zVector, **rVector, *LFBlock,
  *LyVector, *FzVector;           /*!< \brief Arrays of the Linelet preconditioner methodology. */
  unsigned long max_nElem;
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  TCSysMatrix(void);
  
  /*!
   * \brief Destructor of the class.
   */
  ~TCSysMatrix(void);
  
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
   * \brief Copies the block (i, j) of the matrix-by-blocks structure in the internal variable *block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   */
  CalcType *GetBlock(unsigned long block_i, unsigned long block_j);
  
  /*!
   * \brief Copies the block (i, j) of the matrix-by-blocks structure in the internal variable *block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   */
  CalcType GetBlock(unsigned long block_i, unsigned long block_j, unsigned short iVar, unsigned short jVar);
  
  /*!
   * \brief Set the value of a block in the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to set to A(i, j).
   */
  void SetBlock(unsigned long block_i, unsigned long block_j, CalcType **val_block);
  
  /*!
   * \brief Set the value of a block in the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to set to A(i, j).
   */
  void SetBlock(unsigned long block_i, unsigned long block_j, CalcType *val_block);
  
  /*!
   * \brief Adds the specified block to the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to add to A(i, j).
   */
  void AddBlock(unsigned long block_i, unsigned long block_j, CalcType **val_block);
  
  /*!
   * \brief Subtracts the specified block to the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to subtract to A(i, j).
   */
  void SubtractBlock(unsigned long block_i, unsigned long block_j, CalcType **val_block);
  
  /*!
   * \brief Copies the block (i, j) of the matrix-by-blocks structure in the internal variable *block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   */
  CalcType *GetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j);
  
  /*!
   * \brief Set the value of a block in the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to set to A(i, j).
   */
  void SetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, CalcType *val_block);
  
  
  /*!
   * \brief Set the transposed value of a block in the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to set to A(i, j).
   */
  void SetBlockTransposed_ILUMatrix(unsigned long block_i, unsigned long block_j, CalcType *val_block);
  
  /*!
   * \brief Subtracts the specified block to the sparse matrix.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] **val_block - Block to subtract to A(i, j).
   */
  void SubtractBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, CalcType *val_block);
  
  /*!
   * \brief Adds the specified value to the diagonal of the (i, i) subblock
   *        of the matrix-by-blocks structure.
   * \param[in] block_i - Index of the block in the matrix-by-blocks structure.
   * \param[in] val_matrix - Value to add to the diagonal elements of A(i, i).
   */
  void AddVal2Diag(unsigned long block_i, CalcType val_matrix);
  
  /*!
   * \brief Sets the specified value to the diagonal of the (i, i) subblock
   *        of the matrix-by-blocks structure.
   * \param[in] block_i - Index of the block in the matrix-by-blocks structure.
   * \param[in] val_matrix - Value to add to the diagonal elements of A(i, i).
   */
  void SetVal2Diag(unsigned long block_i, CalcType val_matrix);
  
  /*!
   * \brief Calculates the matrix-vector product
   * \param[in] matrix
   * \param[in] vector
   * \param[out] product
   */
  void MatrixVectorProduct(CalcType *matrix, CalcType *vector, CalcType *product);
  
  /*!
   * \brief Calculates the matrix-matrix product
   * \param[in] matrix_a
   * \param[in] matrix_b
   * \param[out] product
   */
  void MatrixMatrixProduct(CalcType *matrix_a, CalcType *matrix_b, CalcType *product);
  
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
  CalcType MatrixDeterminant(CalcType **a, unsigned long n);
  
  /*!
   * \brief Find the cofactor matrix of a square matrix. Written by Paul Bourke
   * \param[in] a - Matrix to compute the determinant.
   * \param[in] n - Size of the quare matrix.
   * \param[out] b - cofactor matrix
   */
  void MatrixCoFactor(CalcType **a, unsigned long n, CalcType **b) ;
  
  /*!
   * \brief Transpose of a square matrix, do it in place. Written by Paul Bourke
   * \param[in] a - Matrix to compute the determinant.
   * \param[in] n - Size of the quare matrix.
   */
  void MatrixTranspose(CalcType **a, unsigned long n) ;
  
  /*!
   * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i, i) subblock and rhs.
   * \param[in] block_i - Index of the (i, i) subblock in the matrix-by-blocks structure.
   * \param[in] rhs - Right-hand-side of the linear system.
   * \param[in] transposed - If true the transposed of the block is used (default = false).
   * \return Solution of the linear system (overwritten on rhs).
   */
  void Gauss_Elimination(unsigned long block_i, CalcType* rhs, bool transposed = false);
  
  /*!
   * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i, i) subblock and rhs.
   * \param[in] Block - matrix-by-blocks structure.
   * \param[in] rhs - Right-hand-side of the linear system.
   * \return Solution of the linear system (overwritten on rhs).
   */
  void Gauss_Elimination(CalcType* Block, CalcType* rhs);
  
  /*!
   * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i, i) subblock and rhs.
   * \param[in] block_i - Index of the (i, i) subblock in the matrix-by-blocks structure.
   * \param[in] rhs - Right-hand-side of the linear system.
   * \return Solution of the linear system (overwritten on rhs).
   */
  void Gauss_Elimination_ILUMatrix(unsigned long block_i, CalcType* rhs);
  
  /*!
   * \fn void CSysMatrix::ProdBlockVector(unsigned long block_i, unsigned long block_j, CalcType* vec);
   * \brief Performs the product of the block (i, j) by vector vec.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
   * \param[in] vec - Vector to be multiplied by the block (i, j) of the sparse matrix A.
   * \return Product of A(i, j) by vector *vec (stored at *prod_block_vector).
   */
  void ProdBlockVector(unsigned long block_i, unsigned long block_j, const TCSysVector<CalcType> & vec);
  
  /*!
   * \brief Performs the product of i-th row of the upper part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the upper part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return prod Result of the product U(A)*vec (stored at *prod_row_vector).
   */
  void UpperProduct(TCSysVector<CalcType> & vec, unsigned long row_i);
  
  /*!
   * \brief Performs the product of i-th row of the lower part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the lower part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return prod Result of the product L(A)*vec (stored at *prod_row_vector).
   */
  void LowerProduct(TCSysVector<CalcType> & vec, unsigned long row_i);
  
  /*!
   * \brief Performs the product of i-th row of the diagonal part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the diagonal part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return prod Result of the product D(A)*vec (stored at *prod_row_vector).
   */
  void DiagonalProduct(TCSysVector<CalcType> & vec, unsigned long row_i);
  
  /*!
   * \brief Send receive the solution using MPI.
   * \param[in] x - Solution..
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SendReceive_Solution(TCSysVector<CalcType> & x, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Send receive the solution using MPI and the transposed structure of the matrix.
   * \param[in] x - Solution..
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SendReceive_SolutionTransposed(TCSysVector<CalcType> & x, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Performs the product of i-th row of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the row of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return Result of the product (stored at *prod_row_vector).
   */
  void RowProduct(const TCSysVector<CalcType> & vec, unsigned long row_i);
  
  /*!
   * \brief Performs the product of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the sparse matrix A.
   * \param[out] prod - Result of the product.
   * \return Result of the product A*vec.
   */
  void MatrixVectorProduct(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod);
  
  /*!
   * \brief Performs the product of a sparse matrix by a CSysVector.
   * \param[in] vec - TCSysVector<CalcType> to be multiplied by the sparse matrix A.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] prod - Result of the product.
   */
  void MatrixVectorProduct(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Performs the product of a sparse matrix by a CSysVector.
   * \param[in] vec - TCSysVector<CalcType> to be multiplied by the sparse matrix A.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] prod - Result of the product.
   */
  void MatrixVectorProductTransposed(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Performs the product of two block matrices.
   */
  void GetMultBlockBlock(CalcType *c, CalcType *a, CalcType *b);
  
  /*!
   * \brief Performs the product of a block matrices by a vector.
   */
  void GetMultBlockVector(CalcType *c, CalcType *a, CalcType *b);
  
  /*!
   * \brief Performs the subtraction of two matrices.
   */
  void GetSubsBlock(CalcType *c, CalcType *a, CalcType *b);
  
  /*!
   * \brief Performs the subtraction of two vectors.
   */
  void GetSubsVector(CalcType *c, CalcType *a, CalcType *b);
  
  /*!
   * \brief Inverse diagonal block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[out] invBlock - Inverse block.
   */
  void InverseDiagonalBlock(unsigned long block_i, CalcType *invBlock, bool transpose = false);
  
 	/*!
   * \brief Inverse diagonal block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[out] invBlock - Inverse block.
   */
  void InverseDiagonalBlock_ILUMatrix(unsigned long block_i, CalcType *invBlock);
  
  /*!
   * \brief Inverse a block.
   * \param[in] Block - block matrix.
   * \param[out] invBlock - Inverse block.
   */
  void InverseBlock(CalcType *Block, CalcType *invBlock);
  
  /*!
   * \brief Build the Jacobi preconditioner.
   */
  void BuildJacobiPreconditioner(bool transpose = false);
  
  /*!
   * \brief Multiply TCSysVector<CalcType> by the preconditioner
   * \param[in] vec - TCSysVector<CalcType> to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeJacobiPreconditioner(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Apply Jacobi as a classical iterative smoother
   * \param[in] b - TCSysVector<CalcType> containing the residual (b)
   * \param[in] x - TCSysVector<CalcType> containing the solution (x^k)
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] residual
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] x - TCSysVector<CalcType> containing the result of the smoothing (x^k+1 = x^k + M^-1*(b - A*x^k).
   */
  unsigned long Jacobi_Smoother(const TCSysVector<CalcType> & b, TCSysVector<CalcType> & x, TCMatrixVectorProduct<CalcType> & mat_vec, CalcType tol, unsigned long m, CalcType *residual, bool monitoring, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Build the ILU preconditioner.
   * \param[in] transposed - Flag to use the transposed matrix to construct the preconditioner.
   */
  void BuildILUPreconditioner(bool transposed = false);
  
  /*!
   * \brief Multiply TCSysVector<CalcType> by the preconditioner
   * \param[in] vec - TCSysVector<CalcType> to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeILUPreconditioner(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Apply ILU as a classical iterative smoother
   * \param[in] b - TCSysVector<CalcType> containing the residual (b)
   * \param[in] x - TCSysVector<CalcType> containing the solution (x^k)
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] residual
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] x - TCSysVector<CalcType> containing the result of the smoothing (x^k+1 = x^k + M^-1*(b - A*x^k).
   */
  unsigned long ILU_Smoother(const TCSysVector<CalcType> & b, TCSysVector<CalcType> & x, TCMatrixVectorProduct<CalcType> & mat_vec, CalcType tol, unsigned long m, CalcType *residual, bool monitoring, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Multiply TCSysVector<CalcType> by the preconditioner
   * \param[in] vec - TCSysVector<CalcType> to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   */
  void ComputeLU_SGSPreconditioner(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Apply LU_SGS as a classical iterative smoother
   * \param[in] b - TCSysVector<CalcType> containing the residual (b)
   * \param[in] x - TCSysVector<CalcType> containing the solution (x^k)
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] residual
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] x - TCSysVector<CalcType> containing the result of the smoothing (x^k+1 = x^k + M^-1*(b - A*x^k).
   */
  unsigned long LU_SGS_Smoother(const TCSysVector<CalcType> & b, TCSysVector<CalcType> & x, TCMatrixVectorProduct<CalcType> & mat_vec, CalcType tol, unsigned long m, CalcType *residual, bool monitoring, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Build the Linelet preconditioner.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  unsigned short BuildLineletPreconditioner(CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Multiply TCSysVector<CalcType> by the preconditioner
   * \param[in] vec - TCSysVector<CalcType> to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   */
  void ComputeLineletPreconditioner(const TCSysVector<CalcType> & vec, TCSysVector<CalcType> & prod, CGeometry *geometry, CConfig *config);
  
  /*!
   * \brief Compute the residual Ax-b
   * \param[in] sol - TCSysVector<CalcType> to be multiplied by the preconditioner.
   * \param[in] f - Result of the product A*vec.
   * \param[out] res - Result of the product A*vec.
   */
  void ComputeResidual(const TCSysVector<CalcType> & sol, const TCSysVector<CalcType> & f, TCSysVector<CalcType> & res);
  
};

typedef TCSysMatrix<su2double> CSysMatrix;

/*!
 * \class TCSysMatrixVectorProduct
 * \brief specialization of matrix-vector product that uses CSysMatrix class
 */
template<class CalcType>
class TCSysMatrixVectorProduct : public TCMatrixVectorProduct<CalcType> {
private:
  TCSysMatrix<CalcType>* sparse_matrix; /*!< \brief pointer to matrix that defines the product. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the products
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  TCSysMatrixVectorProduct(TCSysMatrix<CalcType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~TCSysMatrixVectorProduct() {}
  
  /*!
   * \brief operator that defines the CSysMatrix-TCSysVector<CalcType> product
   * \param[in] u - TCSysVector<CalcType> that is being multiplied by the sparse matrix
   * \param[out] v - TCSysVector<CalcType> that is the result of the product
   */
  void operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const;
};

typedef TCSysMatrixVectorProduct<su2double> CSysMatrixVectorProduct;

/*!
 * \class TCSysMatrixVectorProduct
 * \brief specialization of matrix-vector product that uses CSysMatrix class
 */
template<class CalcType>
class TCSysMatrixVectorProductTransposed : public TCMatrixVectorProduct<CalcType> {
private:
  TCSysMatrix<CalcType>* sparse_matrix; /*!< \brief pointer to matrix that defines the product. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the products
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  TCSysMatrixVectorProductTransposed(TCSysMatrix<CalcType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~TCSysMatrixVectorProductTransposed() {}
  
  /*!
   * \brief operator that defines the CSysMatrix-TCSysVector<CalcType> product
   * \param[in] u - TCSysVector<CalcType> that is being multiplied by the sparse matrix
   * \param[out] v - TCSysVector<CalcType> that is the result of the product
   */
  void operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const;
};

typedef TCSysMatrixVectorProductTransposed<su2double> CSysMatrixVectorProductTransposed;

/*!
 * \class TCJacobiPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class CalcType>
class TCJacobiPreconditioner : public TCPreconditioner<CalcType> {
private:
  TCSysMatrix<CalcType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  TCJacobiPreconditioner(TCSysMatrix<CalcType>& matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~TCJacobiPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - TCSysVector<CalcType> that is being preconditioned
   * \param[out] v - TCSysVector<CalcType> that is the result of the preconditioning
   */
  void operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const;
};

typedef TCJacobiPreconditioner<su2double> CJacobiPreconditioner;

/*!
 * \class TCJacobiTransposedPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class CalcType>
class TCJacobiTransposedPreconditioner : public TCPreconditioner<CalcType> {
private:
  TCSysMatrix<CalcType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  TCJacobiTransposedPreconditioner(TCSysMatrix<CalcType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~TCJacobiTransposedPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - TCSysVector<CalcType> that is being preconditioned
   * \param[out] v - TCSysVector<CalcType> that is the result of the preconditioning
   */
  void operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const;
};

typedef TCJacobiTransposedPreconditioner<su2double> CJacobiTransposedPreconditioner;

/*!
 * \class CILUPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class CalcType>
class TCILUPreconditioner : public TCPreconditioner<CalcType> {
private:
  TCSysMatrix<CalcType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  TCILUPreconditioner(TCSysMatrix<CalcType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~TCILUPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - TCSysVector<CalcType> that is being preconditioned
   * \param[out] v - TCSysVector<CalcType> that is the result of the preconditioning
   */
  void operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const;
};

typedef TCILUPreconditioner<su2double> CILUPreconditioner;

/*!
 * \class TCLU_SGSPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class CalcType>
class TCLU_SGSPreconditioner : public TCPreconditioner<CalcType> {
private:
  TCSysMatrix<CalcType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  TCLU_SGSPreconditioner(TCSysMatrix<CalcType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~TCLU_SGSPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - TCSysVector<CalcType> that is being preconditioned
   * \param[out] v - TCSysVector<CalcType> that is the result of the preconditioning
   */
  void operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const;
};

typedef TCLU_SGSPreconditioner<su2double> CLU_SGSPreconditioner;

/*!
 * \class TCLineletPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class CalcType>
class TCLineletPreconditioner : public TCPreconditioner<CalcType> {
private:
  TCSysMatrix<CalcType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref -
   * \param[in] config_ref -
   */
  TCLineletPreconditioner(TCSysMatrix<CalcType> & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
  /*!
   * \brief destructor of the class
   */
  ~TCLineletPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - TCSysVector<CalcType> that is being preconditioned
   * \param[out] v - TCSysVector<CalcType> that is the result of the preconditioning
   */
  void operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v) const;
};

typedef TCLineletPreconditioner<su2double> CLineletPreconditioner;

template class TCSysMatrix<su2double>;
#ifdef CODI_REVERSE_TYPE
template class TCSysMatrix<passivedouble>;
#endif
#include "matrix_structure.inl"

