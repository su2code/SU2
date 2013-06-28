/*!
 * \file sparse_structure.hpp
 * \brief Headers of the main subroutines for creating the sparse matrices-by-blocks.
 *        The subroutines and functions are in the <i>sparse_structure.cpp</i> file.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <iostream>
#include <cmath>

#include "config_structure.hpp"
#include "linear_solvers_structure.hpp"

using namespace std;

/*! 
 * \class CSparseMatrix
 * \brief Main class for defining sparse matrices-by-blocks 
          with compressed row format.
 * \author A. Bueno.
 * \version 1.0.
 */
class CSparseMatrix {
private:
	unsigned long nPoint;      /*!< \brief Number of points in the grid. */
	unsigned long nVar;        /*!< \brief Number of variables. */
	double *val;               /*!< \brief Entries of the sparse matrix. */
	unsigned long *row_ptr;    /*!< \brief Pointers to the first element in each row. */
	unsigned long *col_ind;    /*!< \brief Column index for each of the elements in val(). */
	unsigned long nnz;         /*!< \brief Number of possible nonzero entries in the matrix. */
	double *block;             /*!< \brief Internal array to store a subblock of the matrix. */
	double *prod_block_vector; /*!< \brief Internal array to store the product of a subblock with a vector. */
	double *prod_row_vector;   /*!< \brief Internal array to store the product of a matrix-by-blocks "row" with a vector. */
	double *aux_vector;		   /*!< \brief Auxilar array to store intermediate results. */	
	double *invM;              /*!< \brief Inverse of (Jacobi) preconditioner. */
	
public:
	
	/*! 
	 * \brief Constructor of the class. 
	 */
	CSparseMatrix(void);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSparseMatrix(void);
	
	/*! 
	 * \brief Assings values to the sparse-matrix structure. 
	 * \param[in] val_nPoint - Number of points in the nPoint x nPoint block structure
	 * \param[in] val_nVar - Number of nVar x nVar variables in each subblock of the matrix-by-block structure.
	 * \param[in] val_row_ptr - Pointers to the first element in each row.
	 * \param[in] val_col_ind - Column index for each of the elements in val().
	 * \param[in] val_nnz - Number of possible nonzero entries in the matrix.
	 * \param[in] preconditioner - If <code>TRUE</code> then it use a preconditioner.
	 */
	void SetIndexes(unsigned long val_nPoint, unsigned short val_nVar, unsigned long* val_row_ptr, 
					unsigned long* val_col_ind, unsigned long val_nnz, bool preconditioner);
	
	/*! 
	 * \brief Sets to zero all the entries of the sparse matrix. 
	 */
	void SetValZero(void);
	
	/*! 
	 * \brief Scales the entries of the sparse matrix. 
	 * \param[in] val_scale - Factor of scaling.
	 */
	void ScaleVals(double val_scale);
	
	/*! 
	 * \brief Copies the block (i,j) of the matrix-by-blocks structure in the internal variable *block.
	 * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
	 */
	void GetBlock(unsigned long block_i, unsigned long block_j);
	
	/*! 
	 * \brief Displays the content of the internal variable <i>*block</i> (for debug purposes).
	 */
	void DisplayBlock(void);
	
	/*! 
	 * \brief Adds the specified block to the sparse matrix. 
	 * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] **val_block - Block to add to A(i,j).
	 */
	void AddBlock(unsigned long block_i, unsigned long block_j, double **val_block);
	
	/*! 
	 * \brief Subtracts the specified block to the sparse matrix. 
	 * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] **val_block - Block to subtract to A(i,j).
	 */
	void SubtractBlock(unsigned long block_i, unsigned long block_j, double **val_block);
	
	/*! 
	 * \brief Adds the specified value to the diagonal of the (i,i) subblock 
	 *        of the matrix-by-blocks structure. 
	 * \param[in] block_i - Index of the block in the matrix-by-blocks structure.
	 * \param[in] val_val - Value to add to the diagonal elements of A(i,i).
	 */
	void AddVal2Diag(unsigned long block_i, double val_val);
	
	/*! 
	 * \brief Deletes the values of the row i of the sparse matrix.
	 * \param[in] i - Index of the row.
	 */
	void DeleteValsRowi(unsigned long i);
	
	/*! 
	 * \brief Returns the sum of the row i.
	 * \param[in] i - Index of the row.
	 * \return The sum of the row i.
	 */
	double SumAbsRowi(unsigned long i);
	
	/*! 
	 * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i,i) subblock and rhs. 
	 * \param[in] block_i - Index of the (i,i) subblock in the matrix-by-blocks structure.
	 * \param[in] rhs - Right-hand-side of the linear system.
	 * \return Solution of the linear system (overwritten on rhs).
	 */
	void Gauss_Elimination(unsigned long block_i, double* rhs);
	
	/*! 
	 * \fn void CSparseMatrix::ProdBlockVector(unsigned long block_i, unsigned long block_j, double* vec);
	 * \brief Performs the product of the block (i,j) by vector vec. 
	 * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] vec - Vector to be multiplied by the block (i,j) of the sparse matrix A.
	 * \return Product of A(i,j) by vector *vec (stored at *prod_block_vector).
	 */
	void ProdBlockVector(unsigned long block_i, unsigned long block_j, double* vec);
	
	/*! 
	 * \brief Performs the product of i-th row of the upper part of a sparse matrix by a vector.
	 * \param[in] vec - Vector to be multiplied by the upper part of the sparse matrix A.
	 * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
	 * \return prod Result of the product U(A)*vec (stored at *prod_row_vector).
	 */
	void UpperProduct(double* vec, unsigned long row_i);
	
	/*! 
	 * \brief Performs the product of i-th row of the lower part of a sparse matrix by a vector.
	 * \param[in] vec - Vector to be multiplied by the lower part of the sparse matrix A.
	 * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
	 * \return prod Result of the product L(A)*vec (stored at *prod_row_vector).
	 */
	void LowerProduct(double* vec, unsigned long row_i);
	
	/*! 
	 * \brief Performs the product of i-th row of the diagonal part of a sparse matrix by a vector.
	 * \param[in] vec - Vector to be multiplied by the diagonal part of the sparse matrix A.
	 * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
	 * \return prod Result of the product D(A)*vec (stored at *prod_row_vector).
	 */
	void DiagonalProduct(double* vec, unsigned long row_i);
	
	/*! 
	 * \brief Performs a single LU-Symmetric Gauss Seidel (SGS) iteration over vector x (overwriten on x_n).
	 * \param[in] b - RHS of the equation.
	 * \param[in] x_n - Approximate solution at iteration (1). 
	 */
	void LU_SGSIteration(double* b, double* x_n);

	/*! 
	 * \brief Performs a Symmetric Gauss Seidel (SGS) iteration over vector x (overwriten on x_n).
	 * \param[in] b - RHS of the equation.
	 * \param[in] x_n - Approximate solution at iteration (n). 
	 * \return Norm of the residual at iteration (n+1). 
	 */
	double SGSIteration(double* b, double* x_n);
	
	/*! 
	 * \brief Solves the linear system Ax = b using the Symmetric Gauss Seidel (SGS) algorithm. 
	 * \param[in] b - RHS of the equation.
	 * \param[in] x_i - Initial candidate for the solution (x_i is overwritten).
	 * \param[in] tol - Tolerance in order to stop the iterative proccess.
	 * \param[in] max_it - Maximum number of iterations.
	 * \param[in] monitoring - Boolean variable in order to monitore the convergence proccess.
	 * \return Approximate solution to the linear system Ax = b with tolerance tol (overwritten on x_i).
	 */
	void SGSSolution(double* b, double* x_i, double tol, int max_it, bool monitoring);
	
	/*! 
	 * \brief Performs the product of i-th row of a sparse matrix by a vector.
	 * \param[in] vec - Vector to be multiplied by the row of the sparse matrix A.
	 * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
	 * \return Result of the product (stored at *prod_row_vector).
	 */
	void RowProduct(double* vec, unsigned long row_i);
	
	/*!
	 * \brief Performs the product of a sparse matrix by a vector.
	 * \param[in] vec - Vector to be multiplied by the sparse matrix A.
	 * \param[out] prod - Result of the product.
	 * \return Result of the product A*vec.
	 */
	void MatrixVectorProduct(double* vec, double* prod);

	/*!
	 * \brief Performs the product of a sparse matrix by a CSysVector.
	 * \param[in] vec - CSysVector to be multiplied by the sparse matrix A.
	 * \param[out] prod - Result of the product.
	 */
  void MatrixVectorProduct(const CSysVector & vec, CSysVector & prod);
	
	/*! 
	 * \brief Solves the linear system Ax = b using the Conjugate Gradient (CG) algorithm. 
	 * \param[in] b - RHS of the equation.
	 * \param[in] x_i - Initial candidate for the solution (x_i is overwritten).
	 * \param[in] tol - Tolerance in order to stop the iterative proccess.
	 * \param[in] max_it - Maximum number of iterations.
	 * \param[in] monitoring - Boolean variable in order to monitore the convergence proccess.
	 * \return Approximate solution to the linear system Ax = b with tolerance tol (overwritten on x_i).
	 */
	void CGSolution(double* b, double* x_i, double tol, int max_it, bool monitoring);
	
	/*! 
	 * \brief Solves the linear system Ax = b using a preconditioned Conjugate Gradient (CG) algorithm. 
	 * \param[in] b - RHS of the equation.
	 * \param[in] x_i - Initial candidate for the solution (x_i is overwritten).
	 * \param[in] tol - Tolerance in order to stop the iterative proccess.
	 * \param[in] max_it - Maximum number of iterations.
	 * \param[in] monitoring - Boolean variable in order to monitore the convergence proccess.
	 */
	void PreconditionedCGSolution(double* b, double* x_i, double tol, int max_it, bool monitoring); 
	
	/*! 
	 * \brief Inverse diagonal block. 
	 * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
	 * \param[out] invBlock - Inverse block.
	 */
	void InverseDiagonalBlock(unsigned long block_i, double **invBlock);
	
	/*! 
	 * \brief Bouild the Jacobi preconditioner. 
	 */
	void BuildJacobiPreconditioner(void);
	
	/*! 
	 * \brief Multiply the preconditioner by a vector. 
	 * \param[in] vec - Vector to be multiplied by the preconditioner.
	 * \param[out] prod - Result of the product A*vec.
	 */
	void PrecondVectorProduct(double* vec, double* prod);

	/*! 
	 * \brief Multiply CSysVector by the preconditioner
	 * \param[in] vec - CSysVector to be multiplied by the preconditioner.
	 * \param[out] prod - Result of the product A*vec.
	 */
	void PrecondVectorProduct(const CSysVector & vec, CSysVector & prod);
};

/*!
 * \class CSparseMatrixVectorProduct
 * \brief specialization of matrix-vector product that uses CSparseMatrix class
 */
class CSparseMatrixVectorProduct : public CMatrixVectorProduct {
private:
  CSparseMatrix* sparse_matrix; /*!< \brief pointer to matrix that defines the product. */

 public:

  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the products
   */
  CSparseMatrixVectorProduct(CSparseMatrix & matrix_ref);

  /*!
   * \brief destructor of the class
   */
  ~CSparseMatrixVectorProduct(){}

  /*!
   * \brief operator that defines the CSparseMatrix-CSysVector product
   * \param[in] u - CSysVector that is being multiplied by the sparse matrix
   * \param[out] v - CSysVector that is the result of the product
   */
  void operator()(const CSysVector & u, CSysVector & v) const;
};

/*!
 * \class CSparseMatrixPreconditioner
 * \brief specialization of preconditioner that uses CSparseMatrix class
 */
class CSparseMatrixPreconditioner : public CPreconditioner {
private:
  CSparseMatrix* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */

 public:

  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   */
  CSparseMatrixPreconditioner(CSparseMatrix & matrix_ref);

  /*!
   * \brief destructor of the class
   */
  ~CSparseMatrixPreconditioner() {}

  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - CSysVector that is being preconditioned
   * \param[out] v - CSysVector that is the result of the preconditioning
   */
  void operator()(const CSysVector & u, CSysVector & v) const;
};

/*!
 * \class CIdentityPreconditioner
 * \brief specialization of preconditioner that does nothing (leaves vector unchanged)
 */
class CIdentityPreconditioner : public CPreconditioner {
 public:
  
  /*!
   * \brief default constructor of the class
   */
  CIdentityPreconditioner() {}

  /*!
   * \brief destructor of the class
   */
  ~CIdentityPreconditioner() {}

  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - CSysVector that is being preconditioned
   * \param[out] v - CSysVector that is the result of the preconditioning
   */
  void operator()(const CSysVector & u, CSysVector & v) const;
};

#include "sparse_structure.inl"
