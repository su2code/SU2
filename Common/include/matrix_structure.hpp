/*!
 * \file matrix_structure.hpp
 * \brief Headers of the main subroutines for creating the sparse matrices-by-blocks.
 *        The subroutines and functions are in the <i>matrix_structure.cpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
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
 * \class CSysMatrix
 * \brief Main class for defining sparse matrices-by-blocks
 with compressed row format.
 * \author A. Bueno, F. Palacios.
 * \version 2.0.6
 */
class CSysMatrix {
private:
	unsigned long nPoint,   /*!< \brief Number of points in the grid. */
	nPointDomain,           /*!< \brief Number of points in the grid. */
	nVar,                   /*!< \brief Number of variables. */
	nEqn;                   /*!< \brief Number of equations. */
	double *matrix;               /*!< \brief Entries of the sparse matrix. */
	unsigned long *row_ptr;    /*!< \brief Pointers to the first element in each row. */
	unsigned long *col_ind;    /*!< \brief Column index for each of the elements in val(). */
	unsigned long nnz;         /*!< \brief Number of possible nonzero entries in the matrix. */
	double *block;             /*!< \brief Internal array to store a subblock of the matrix. */
	double *prod_block_vector; /*!< \brief Internal array to store the product of a subblock with a vector. */
	double *prod_row_vector;   /*!< \brief Internal array to store the product of a matrix-by-blocks "row" with a vector. */
	double *aux_vector;		   /*!< \brief Auxilar array to store intermediate results. */
	double *invM;              /*!< \brief Inverse of (Jacobi) preconditioner. */
	bool *LineletBool;						 /*!< \brief Identify if a point belong to a linelet. */
	vector<unsigned long> *LineletPoint;	 /*!< \brief Linelet structure. */
	unsigned long nLinelet;							 /*!< \brief Number of Linelets in the system. */
  
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
  void Initialize(unsigned long nPoint, unsigned long nPointDomain, unsigned short nVar, unsigned short nEqn, CGeometry *geometry);
  
  /*!
	 * \brief Assings values to the sparse-matrix structure.
	 * \param[in] val_nPoint - Number of points in the nPoint x nPoint block structure
	 * \param[in] val_nVar - Number of nVar x nVar variables in each subblock of the matrix-by-block structure.
   * \param[in] val_nEq - Number of nEqn x nVar variables in each subblock of the matrix-by-block structure.
	 * \param[in] val_row_ptr - Pointers to the first element in each row.
	 * \param[in] val_col_ind - Column index for each of the elements in val().
	 * \param[in] val_nnz - Number of possible nonzero entries in the matrix.
	 * \param[in] preconditioner - If <code>TRUE</code> then it use a preconditioner.
	 */
	void SetIndexes(unsigned long val_nPoint, unsigned long val_nPointDomain, unsigned short val_nVar, unsigned short val_nEq, unsigned long* val_row_ptr, unsigned long* val_col_ind, unsigned long val_nnz);
  
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
	 * \brief Adds the specified block to the sparse matrix.
	 * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] **val_block - Block to add to A(i,j).
	 */
	void AddBlock(unsigned long point_i, double **val_block);
  
	/*!
	 * \brief Subtracts the specified block to the sparse matrix.
	 * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] **val_block - Block to subtract to A(i,j).
	 */
	void SubtractBlock(unsigned long point_i, double **val_block);
  
	/*!
	 * \brief Adds the specified value to the diagonal of the (i,i) subblock
	 *        of the matrix-by-blocks structure.
	 * \param[in] block_i - Index of the block in the matrix-by-blocks structure.
	 * \param[in] val_matrix - Value to add to the diagonal elements of A(i,i).
	 */
	void AddVal2Diag(unsigned long block_i, double val_matrix);
  
	/*!
	 * \brief Adds the specified value to the diagonal of the (i,i) subblock
	 *        of the matrix-by-blocks structure.
	 * \param[in] block_i - Index of the block in the matrix-by-blocks structure.
	 * \param[in] *val_matrix - Values to add to the diagonal elements of A(i,i).
	 * \param[in] num_dim - number of dimensions
	 *
	 */
	void AddVal2Diag(unsigned long block_i, double* val_matrix, unsigned short num_dim);
  
  /*!
	 * \brief Adds the specified value to the diagonal of the (i,i) subblock
	 *        of the matrix-by-blocks structure.
	 * \param[in] block_i - Index of the block in the matrix-by-blocks structure.
	 * \param[in] *val_matrix - Values to add to the diagonal elements of A(i,i).
	 * \param[in] val_nDim - number of dimensions
   * \param[in] val_nDiatomics - number of diatomic species
	 *
	 */
	void AddVal2Diag(unsigned long block_i, double* val_matrix, unsigned short val_nDim, unsigned short val_nDiatomics);
  
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
	 * \fn void CSysMatrix::ProdBlockVector(unsigned long block_i, unsigned long block_j, double* vec);
	 * \brief Performs the product of the block (i,j) by vector vec.
	 * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
	 * \param[in] vec - Vector to be multiplied by the block (i,j) of the sparse matrix A.
	 * \return Product of A(i,j) by vector *vec (stored at *prod_block_vector).
	 */
	void ProdBlockVector(unsigned long block_i, unsigned long block_j, const CSysVector & vec);
  
  /*!
	 * \brief Performs the product of i-th row of the upper part of a sparse matrix by a vector.
	 * \param[in] vec - Vector to be multiplied by the upper part of the sparse matrix A.
	 * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
	 * \return prod Result of the product U(A)*vec (stored at *prod_row_vector).
	 */
	void UpperProduct(CSysVector & vec, unsigned long row_i);
  
  /*!
	 * \brief Performs the product of i-th row of the lower part of a sparse matrix by a vector.
	 * \param[in] vec - Vector to be multiplied by the lower part of the sparse matrix A.
	 * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
	 * \return prod Result of the product L(A)*vec (stored at *prod_row_vector).
	 */
	void LowerProduct(CSysVector & vec, unsigned long row_i);
  
  /*!
	 * \brief Performs the product of i-th row of the diagonal part of a sparse matrix by a vector.
	 * \param[in] vec - Vector to be multiplied by the diagonal part of the sparse matrix A.
	 * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
	 * \return prod Result of the product D(A)*vec (stored at *prod_row_vector).
	 */
	void DiagonalProduct(CSysVector & vec, unsigned long row_i);
	
  /*!
	 * \brief Send receive the solution using MPI.
	 * \param[in] x - Solution..
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SendReceive_Solution(CSysVector & x, CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Performs the product of i-th row of a sparse matrix by a vector.
	 * \param[in] vec - Vector to be multiplied by the row of the sparse matrix A.
	 * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
	 * \return Result of the product (stored at *prod_row_vector).
	 */
	void RowProduct(const CSysVector & vec, unsigned long row_i);
  
  /*!
	 * \brief Performs the product of a sparse matrix by a vector.
	 * \param[in] vec - Vector to be multiplied by the sparse matrix A.
	 * \param[out] prod - Result of the product.
	 * \return Result of the product A*vec.
	 */
	void MatrixVectorProduct(const CSysVector & vec, CSysVector & prod);
  
	/*!
	 * \brief Performs the product of a sparse matrix by a CSysVector.
	 * \param[in] vec - CSysVector to be multiplied by the sparse matrix A.
	 * \param[out] prod - Result of the product.
	 */
	void MatrixVectorProduct(const CSysVector & vec, CSysVector & prod, CGeometry *geometry, CConfig *config);
	
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
	void ComputeJacobiPreconditioner(const CSysVector & vec, CSysVector & prod, CGeometry *geometry, CConfig *config);
	
  /*!
	 * \brief Multiply CSysVector by the preconditioner
	 * \param[in] vec - CSysVector to be multiplied by the preconditioner.
	 * \param[out] prod - Result of the product A*vec.
	 */
	void ComputeLU_SGSPreconditioner(const CSysVector & vec, CSysVector & prod, CGeometry *geometry, CConfig *config);
  
	/*!
	 * \brief Multiply CSysVector by the preconditioner
	 * \param[in] vec - CSysVector to be multiplied by the preconditioner.
	 * \param[out] prod - Result of the product A*vec.
	 */
	void ComputeLineletPreconditioner(const CSysVector & vec, CSysVector & prod, CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Multiply CSysVector by the preconditioner
	 * \param[in] vec - CSysVector to be multiplied by the preconditioner.
	 * \param[out] prod - Result of the product A*vec.
	 */
	void ComputeIdentityPreconditioner(const CSysVector & vec, CSysVector & prod, CGeometry *geometry, CConfig *config);
	
};

/*!
 * \class CSysMatrixVectorProduct
 * \brief specialization of matrix-vector product that uses CSysMatrix class
 */
class CSysMatrixVectorProduct : public CMatrixVectorProduct {
private:
	CSysMatrix* sparse_matrix; /*!< \brief pointer to matrix that defines the product. */
	CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
	CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
	/*!
	 * \brief constructor of the class
	 * \param[in] matrix_ref - matrix reference that will be used to define the products
	 */
	CSysMatrixVectorProduct(CSysMatrix & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
	/*!
	 * \brief destructor of the class
	 */
	~CSysMatrixVectorProduct(){}
  
	/*!
	 * \brief operator that defines the CSysMatrix-CSysVector product
	 * \param[in] u - CSysVector that is being multiplied by the sparse matrix
	 * \param[out] v - CSysVector that is the result of the product
	 */
	void operator()(const CSysVector & u, CSysVector & v) const;
};

/*!
 * \class CJacobiPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
class CJacobiPreconditioner : public CPreconditioner {
private:
	CSysMatrix* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
	CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
	CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
  
	/*!
	 * \brief constructor of the class
	 * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
	 */
	CJacobiPreconditioner(CSysMatrix & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
  
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
 * \class CLU_SGSPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
class CLU_SGSPreconditioner : public CPreconditioner {
private:
	CSysMatrix* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
	CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
	
	/*!
	 * \brief constructor of the class
	 * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
	 */
	CLU_SGSPreconditioner(CSysMatrix & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
	
	/*!
	 * \brief destructor of the class
	 */
	~CLU_SGSPreconditioner() {}
	
	/*!
	 * \brief operator that defines the preconditioner operation
	 * \param[in] u - CSysVector that is being preconditioned
	 * \param[out] v - CSysVector that is the result of the preconditioning
	 */
	void operator()(const CSysVector & u, CSysVector & v) const;
};

/*!
 * \class CLineletPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
class CLineletPreconditioner : public CPreconditioner {
private:
	CSysMatrix* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
	CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
public:
	
	/*!
	 * \brief constructor of the class
	 * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
	 */
	CLineletPreconditioner(CSysMatrix & matrix_ref, CGeometry *geometry_ref, CConfig *config_ref);
	
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

#include "matrix_structure.inl"
