/*!
 * \file sparse_structure.hpp
 * \brief Headers of the main subroutines for creating the sparse matrices-by-blocks.
 *        The subroutines and functions are in the <i>sparse_structure.cpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.1
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
#include "geometry_structure.hpp"
#include "linear_solvers_structure.hpp"

using namespace std;

/*! 
 * \class CSparseMatrix
 * \brief Main class for defining sparse matrices-by-blocks 
          with compressed row format.
 * \author A. Bueno.
 * \version 2.0.1
 */
class CSparseMatrix {
private:
	unsigned long nPoint;      /*!< \brief Number of points in the grid. */
	unsigned long nPointDomain;/*!< \brief Number of points in the grid. */
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
	unsigned short nSub_blocks; /*!< \brief  Number of sub-blocks in the nVar*nVar block structure */
	unsigned short *Sub_block_sizes;		/*!< \brief  Size of each sub-block in the nVar*nVar block structure */
	bool blockDiagonalJacobian; /*!< \brief flag if the Jacobian has a block diagonal structure like in multi species flow */
	bool *LineletBool;						 /*!< \brief Identify if a point belong to a linelet. */
	vector<unsigned long> *LineletPoint;	 /*!< \brief Linelet structure. */
	unsigned long nLinelet;							 /*!< \brief Number of Linelets in the system. */
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
	 * \param[in] blockDiagonal - If <code>TRUE</code> then modified Gauss Elimination is used for matrix inversion.
	 * \param[in] nSub_blocks - Number of sub-blocks in the nVar*nVar block structure.
	 * \param[in] Sub_block_sizes -Size of each sub-block in the nVar*nVar block structure.
	 */
	void SetIndexes(unsigned long val_nPoint, unsigned long val_nPointDomain, unsigned short val_nVar, unsigned long* val_row_ptr, 
			unsigned long* val_col_ind, unsigned long val_nnz, bool preconditioner, bool blockDiagonal, unsigned short nSub_blocks, 
									unsigned short *Sub_block_sizes);

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
	 * \brief Returns the content of the internal variable <i>*block</i> (for debug purposes).
	 */
	void ReturnBlock(double **val_block);

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
	 * \brief Adds the specified value to the diagonal of the (i,i) subblock
	 *        of the matrix-by-blocks structure.
	 * \param[in] block_i - Index of the block in the matrix-by-blocks structure.
	 * \param[in] *val_val - Values to add to the diagonal elements of A(i,i).
	 * \param[in] num_dim - number of dimensions
	 *
	 */
	void AddVal2Diag(unsigned long block_i, double* val_val, unsigned short num_dim);

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
	 * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i,i) subblock and rhs. 
	 * \param[in] A - matrix-by-blocks structure.
	 * \param[in] rhs - Right-hand-side of the linear system.
	 * \return Solution of the linear system (overwritten on rhs).
	 */
	void Gauss_Elimination(double* Block, double* rhs);

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
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void LU_SGSIteration(double* b, double* x_n, CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Send receive the solution using MPI.
	 * \param[in] x_n - Solution..
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SendReceive_Solution(double* x, CGeometry *geometry, CConfig *config);
	
	/*!
	 * \brief Send receive the solution using MPI.
	 * \param[in] vec - CSysVector to be multiplied by the sparse matrix A.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SendReceive_Solution(CSysVector & vec, CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Solves the linear system Ax = b using the Symmetric Gauss Seidel (SGS) algorithm. 
	 * \param[in] b - RHS of the equation.
	 * \param[in] x_i - Initial candidate for the solution (x_i is overwritten).
	 * \param[in] tol - Tolerance in order to stop the iterative proccess.
	 * \param[in] max_it - Maximum number of iterations.
	 * \param[in] monitoring - Boolean variable in order to monitore the convergence proccess.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \return Approximate solution to the linear system Ax = b with tolerance tol (overwritten on x_i).
	 */
	void SGSSolution(double* b, double* x_i, double tol, int max_it, bool monitoring, CGeometry *geometry, CConfig *config);
	
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
	void MatrixVectorProduct(double* vec, double* prod, double nPoint);

	/*!
	 * \brief Performs the product of a sparse matrix by a CSysVector.
	 * \param[in] vec - CSysVector to be multiplied by the sparse matrix A.
	 * \param[out] prod - Result of the product.
	 */
	void MatrixVectorProduct(const CSysVector & vec, CSysVector & prod);
	
	/*!
	 * \brief Performs the product of two block matrices.
	 */
	void GetMultBlockBlock(double *c, double *a, double *b);
	
	/*!
	 * \brief Performs the product of a block matrices by a vector.
	 */
	void GetMultBlockVector(double *c, double *a, double *b);
	
	/*!
	 * \brief Performs the substraction of two matrices.
	 */
	void GetSubsBlock(double *c, double *a, double *b);
	
	/*!
	 * \brief Performs the substraction of two vectors.
	 */
	void GetSubsVector(double *c, double *a, double *b);
	
	/*! 
	 * \brief Inverse diagonal block. 
	 * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
	 * \param[out] invBlock - Inverse block.
	 */
	void InverseDiagonalBlock(unsigned long block_i, double **invBlock);
	
	/*! 
	 * \brief Inverse a block. 
	 * \param[in] Block - block matrix.
	 * \param[out] invBlock - Inverse block.
	 */
	void InverseBlock(double *Block, double *invBlock);

	/*! 
	 * \brief Build the Jacobi preconditioner. 
	 */
	void BuildJacobiPreconditioner(void);

	/*! 
	 * \brief Build the Linelet preconditioner. 
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void BuildLineletPreconditioner(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Multiply CSysVector by the preconditioner
	 * \param[in] vec - CSysVector to be multiplied by the preconditioner.
	 * \param[out] prod - Result of the product A*vec.
	 */
	void ComputeJacobiPreconditioner(const CSysVector & vec, CSysVector & prod);
	
	/*! 
	 * \brief Multiply CSysVector by the preconditioner
	 * \param[in] vec - CSysVector to be multiplied by the preconditioner.
	 * \param[out] prod - Result of the product A*vec.
	 */
	void ComputeLineletPreconditioner(const CSysVector & vec, CSysVector & prod);

	/*!
	 * \brief Write the header of the history file.
	 * \param[in] LinConvHist_file - Pointer to the convergence history file (which is defined in the main subroutine).
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetLin_Sol_History_Header(ofstream *LinConvHist_file, CConfig *config);

	/*!
	 * \brief Write the history file and the convergence on the screen for serial computations.
	 * \param[in] LinConvHist_file - Pointer to the convergence history file (which is defined in the main subroutine).
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] iExtIter - Current external iteration.
	 * \param[in] linear_resid - Current linear residual.
	 */
	void SetLin_Sol_History_Iter(ofstream *LinConvHist_file, CConfig *config, unsigned long iExtIter, double linear_residual);
	
	/*! 
	 * \brief Solves the linear system Ax = b using a preconditioned Conjugate Gradient (CG) algorithm. 
	 * \param[in] b - RHS of the equation.
	 * \param[in] x_i - Initial candidate for the solution (x_i is overwritten).
	 * \param[in] tol - Tolerance in order to stop the iterative proccess.
	 * \param[in] max_it - Maximum number of iterations.
	 * \param[in] monitoring - Boolean variable in order to monitore the convergence proccess.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void CGSolution(double* b, double* x_i, double tol, int max_it, bool monitoring, 
									CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Multiply the preconditioner by a vector. 
	 * \param[in] vec - Vector to be multiplied by the preconditioner.
	 * \param[out] prod - Result of the product A*vec.
	 */
	void PrecondVectorProduct(double* vec, double* prod, double nPoint);
	
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
 * \class CSparseMatrixSolMPI
 * \brief specialization of matrix-vector product that uses CSparseMatrix class
 */
class CSparseMatrixSolMPI : public CSolutionSendReceive {
private:
	CSparseMatrix* sparse_matrix; /*!< \brief pointer to matrix that defines the product. */
	CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
	CConfig* config; /*!< \brief pointer to matrix that defines the config. */
	
public:
	
	/*!
	 * \brief constructor of the class
	 * \param[in] matrix_ref - matrix reference that will be used to define the products
	 */
	CSparseMatrixSolMPI(CSparseMatrix & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
	
	/*!
	 * \brief destructor of the class
	 */
	CSparseMatrixSolMPI(){}
	
	/*!
	 * \brief operator that defines the CSparseMatrix-CSysVector product
	 * \param[in] u - CSysVector that is being multiplied by the sparse matrix
	 */
	void operator()(CSysVector & u) const;
};

/*!
 * \class CJacobiPreconditioner
 * \brief specialization of preconditioner that uses CSparseMatrix class
 */
class CJacobiPreconditioner : public CPreconditioner {
private:
	CSparseMatrix* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */

public:

	/*!
	 * \brief constructor of the class
	 * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
	 */
	CJacobiPreconditioner(CSparseMatrix & matrix_ref);

	/*!
	 * \brief destructor of the class
	 */
	~CJacobiPreconditioner() {}

	/*!
	 * \brief operator that defines the preconditioner operation
	 * \param[in] u - CSysVector that is being preconditioned
	 * \param[out] v - CSysVector that is the result of the preconditioning
	 */
	void operator()(const CSysVector & u, CSysVector & v) const;
};

/*!
 * \class CLineletPreconditioner
 * \brief specialization of preconditioner that uses CSparseMatrix class
 */
class CLineletPreconditioner : public CPreconditioner {
private:
	CSparseMatrix* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
	
public:
	
	/*!
	 * \brief constructor of the class
	 * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
	 */
	CLineletPreconditioner(CSparseMatrix & matrix_ref);
	
	/*!
	 * \brief destructor of the class
	 */
	~CLineletPreconditioner() {}
	
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
