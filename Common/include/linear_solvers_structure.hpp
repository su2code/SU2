/*!
 * \file linear_solvers_structure.hpp
 * \brief Headers for the classes related to linear solvers (LU_SGS, CG, GMRES, etc)
 *        The subroutines and functions are in the <i>linear_solvers_structure.cpp</i> file.
 * \author Current Development: Stanford University.
 * \version 1.1.
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

#ifndef NO_MPI
#include <mpi.h>
#endif
#include <climits>
#include <limits>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

const double eps = numeric_limits<double>::epsilon(); /*!< \brief machine epsilon */

/*! 
 * \class CSysVector
 * \brief Class for holding and manipulating vectors needed by linear solvers
 * \author J. Hicken.
 * \version 1.1.
 *
 * We could use the STL vector as a base class here, but this gives us
 * more flexibility with the underlying data (e.g. we may decide to
 * use a block storage scheme rather than a continuous storage
 * scheme).
 */
class CSysVector {
  
private:

	unsigned int nElm; /*!< \brief total number of elements (or number elements on this processor) */
  unsigned long nElmGlobal; /*!< \brief total number of elements over all processors */
	unsigned short nVar; /*!< \brief number of elements in a block */
	unsigned int nBlk; /*!< \brief number of blocks (or number of blocks on this processor) */
  unsigned int myrank; /*!< \brief processor rank (only used for parallel runs) */
  double* vec_val; /*!< \brief storage for the element values */

public:
  
  /*! 
   * \brief constructor of the class. 
   * \param[in] size - number of elements locally 
   * \param[in] val - default value for elements
   */
  CSysVector(const unsigned int & size, const double & val = 0.0);
  
  /*! 
   * \brief constructor of the class. 
   * \param[in] numBlk - number of blocks locally 
   * \param[in] numVar - number of variables in each block
   * \param[in] val - default value for elements
   */		      
  CSysVector(const unsigned int & numBlk, const unsigned short & numVar, 
	  const double & val = 0.0);
  
  /*! 
   * \brief class destructor 
   */
  virtual ~CSysVector();
  
  /*! 
   * \brief copy constructor of the class. 
   * \param[in] u - CSysVector that is being copied
   */
  CSysVector(const CSysVector & u);

  /*!
   * \brief constructor from array
   * \param[in] size - number of elements locally 
   * \param[in] u_array - vector stored as array being copied 
   */
  explicit CSysVector(const unsigned int & size, const double* u_array);

  /*!
   * \brief constructor from array
   * \param[in] numBlk - number of blocks locally 
   * \param[in] numVar - number of variables in each block
   * \param[in] u_array - vector stored as array being copied 
   */
  explicit CSysVector(const unsigned int & numBlk, const unsigned short & numVar, 
                      const double* u_array);

  /*!
   * \brief return the number of local elements in the CSysVector
   */
  unsigned int GetLocSize() const;

  /*!
   * \brief return the size of the CSysVector (over all processors)
   */
  unsigned long GetSize() const;

  /*!
   * \brief return the number of variables at each block (typically number per node)
   */
  unsigned short GetNVar() const;

  /*!
   * \brief return the number of blocks (typically number of nodes locally)
   */
  unsigned int GetNBlk() const;  
  
  /*! 
   * \brief set calling CSysVector to scaling of another CSysVector
   * \param[in] a - scalar factor for x
   * \param[in] x - CSysVector that is being scaled
   */
  void Equals_AX(const double & a, CSysVector & x);

  /*! 
   * \brief adds a scaled CSysVector to calling CSysVector
   * \param[in] a - scalar factor for x
   * \param[in] x - CSysVector that is being scaled
   */
  void Plus_AX(const double & a, CSysVector & x);

  /*! 
   * \brief general linear combination of two CSysVectors
   * \param[in] a - scalar factor for x
   * \param[in] x - first CSysVector in linear combination
   * \param[in] b - scalar factor for y
   * \param[in] y - second CSysVector in linear combination
   */
  void Equals_AX_Plus_BY(const double & a, CSysVector & x, const double & b, CSysVector & y);

  /*! 
   * \brief assignment operator with deep copy 
   * \param[in] u - CSysVector whose values are being assigned
   */
  CSysVector & operator=(const CSysVector & u);
  
  /*!
   * \brief CSysVector=double assignment operator
   * \param[in] val - value assigned to each element of CSysVector
   */
  CSysVector & operator=(const double & val);

  /*! 
   * \brief addition operator
   * \param[in] u - CSysVector being added to *this
   */
  CSysVector operator+(const CSysVector & u) const;
  
  /*! 
   * \brief compound addition-assignment operator
   * \param[in] u - CSysVector being added to calling object
   */
  CSysVector & operator+=(const CSysVector & u);

  /*! 
   * \brief subtraction operator
   * \param[in] u - CSysVector being subtracted from *this
   */
  CSysVector operator-(const CSysVector & u) const;

  /*! 
   * \brief compound subtraction-assignment operator
   * \param[in] u - CSysVector being subtracted from calling object
   */
  CSysVector & operator-=(const CSysVector & u);

  /*! 
   * \brief vector * scalar multiplication operator
   * \param[in] val - value to multiply *this by
   */
  CSysVector operator*(const double & val) const;

  /*! 
   * \brief scalar * vector multiplication operator
   * \param[in] val - scalar value to multiply by
   * \param[in] u - CSysVector having its elements scaled
   */
  friend CSysVector operator*(const double & val, const CSysVector & u);

  /*! 
   * \brief compound scalar multiplication-assignment operator
   * \param[in] val - value to multiply calling object by
   */
  CSysVector & operator*=(const double & val);

  /*! 
   * \brief vector-scalar division operator (no scalar/vector operator)
   * \param[in] val - value to divide elements of *this by
   */
  CSysVector operator/(const double & val) const;

  /*! 
   * \brief compound scalar division-assignment operator
   * \param[in] val - value to divide elements of calling object by
   */
  CSysVector & operator/=(const double & val);

  /*! 
   * \brief indexing operator with assignment permitted
   * \param[in] i = local index to access
   */
  double & operator[](const unsigned int & i);

  /*! 
   * \brief indexing operator with assignment not permitted
   * \param[in] i = local index to access
   */
  const double & operator[](const unsigned int & i) const;

  // consider operator()(block,variable) as well?

  /*! 
   * \brief the L2 norm of the CSysVector
   * \result the L2 norm
   */
  double norm() const;

  /*! 
   * \brief copies the contents of the calling CSysVector into an array
   * \param[out] u_array - array into which information is being copied
   * \pre u_array must be allocated and have the same size as CSysVector
   */
  void CopyToArray(double* u_array);

  /*! 
   * \brief dot-product between two CSysVectors 
   * \param[in] u - first CSysVector in dot product
   * \param[in] v - second CSysVector in dot product
   */
  friend double dotProd(const CSysVector & u, const CSysVector & v);

};

/*!
 * \class CMatrixVectorProduct
 * \brief abstract base class for defining matrix-vector products
 * \author J. Hicken.
 * \version 1.1.
 *
 * The Krylov-subspace solvers require only matrix-vector products and
 * not the actual matrix/Jacobian.  We need some way to indicate which
 * function will perform the product.  However, sometimes the
 * functions that define the product will require different numbers
 * and types of inputs.  For example, the forward-difference
 * approximation to a Jacobian-vector product requires the vector that
 * defines the Jacobian and a perturbation parameter.  The
 * CMatrixVectorProduct class is used to derive child classes that can
 * handle the different types of matrix-vector products and still be
 * passed to a single implementation of the Krylov solvers.
 */
class CMatrixVectorProduct {
 public:
  virtual ~CMatrixVectorProduct() = 0; ///< class destructor
  virtual void operator()(const CSysVector & u, CSysVector & v)
      const = 0; ///< matrix-vector product operation
};
inline CMatrixVectorProduct::~CMatrixVectorProduct() {}

/*!
 * \class CPreconditioner
 * \brief abstract base class for defining preconditioning operation
 * \author J. Hicken.
 * \version 1.1.
 *
 * See the remarks regarding the CMatrixVectorProduct class.  The same
 * idea applies here to the preconditioning operation.
 */
class CPreconditioner {
 public:
  virtual ~CPreconditioner() = 0; ///< class destructor
  virtual void operator()(const CSysVector & u, CSysVector & v)
      const = 0; ///< preconditioning operation
};
inline CPreconditioner::~CPreconditioner() {}

/*! 
 * \class CSysSolve
 * \brief Class for solving linear systems using classical and Krylov-subspace iterative methods
 * \author J. Hicken.
 * \version 1.1.
 *
 * The individual solvers could be stand-alone subroutines, but by
 * creating CSysSolve objects we can more easily assign different
 * matrix-vector products and preconditioners to different problems
 * that may arise in a hierarchical solver (i.e. multigrid).
 */
class CSysSolve {

private:
  
  /*!
   * \brief sign transfer function
   * \param[in] x - value having sign prescribed
   * \param[in] y - value that defined the sign
   *
   * this may already be defined as a global function somewhere, if
   * so, feel free to delete this and replace it as needed with the
   * appropriate global function
   */
  double sign(const double & x, const double & y) const;

  /*! 
   * \brief applys a Givens rotation to a 2-vector
   * \param[in] s - sine of the Givens rotation angle
   * \param[in] c - cosine of the Givens rotation angle
   * \param[in,out] h1 - first element of 2x1 vector being transformed
   * \param[in,out] h2 - second element of 2x1 vector being transformed
   */
  void applyGivens(const double & s, const double & c, double & h1, double & h2);

  /*! 
   * \brief generates the Givens rotation matrix for a given 2-vector
   * \param[in,out] dx - element of 2x1 vector being transformed
   * \param[in,out] dy - element of 2x1 vector being set to zero
   * \param[in,out] s - sine of the Givens rotation angle
   * \param[in,out] c - cosine of the Givens rotation angle
   *
   * Based on givens() of SPARSKIT, which is based on p.202 of
   * "Matrix Computations" by Golub and van Loan.
   */
  void generateGivens(double & dx, double & dy, double & s, double & c);

  /*! 
   * \brief finds the solution of the upper triangular system Hsbg*x = rhs
   *
   * \param[in] n - size of the reduced system
   * \param[in] Hsbg - upper triangular matrix
   * \param[in] rhs - right-hand side of the reduced system
   * \param[out] x - solution of the reduced system
   *
   * \pre the upper Hessenberg matrix has been transformed into a 
   * triangular matrix.
   */
  void solveReduced(const int & n, const vector<vector<double> > & Hsbg, 
		    const vector<double> & rhs, vector<double> & x);

  /*! 
   * \brief Modified Gram-Schmidt orthogonalization 
   * \author Based on Kesheng John Wu's mgsro subroutine in Saad's SPARSKIT
   * 
   * \tparam Vec - a generic vector class
   * \param[in] i - index indicating which vector in w is being orthogonalized 
   * \param[in,out] Hsbg - the upper Hessenberg begin updated 
   * \param[in,out] w - the (i+1)th vector of w is orthogonalized against the 
   *                    previous vectors in w
   *
   * \pre the vectors w[0:i] are orthonormal
   * \post the vectors w[0:i+1] are orthonormal
   *
   * Reothogonalization is performed if the cosine of the angle between 
   * w[i+1] and w[k], k < i+1, is greater than 0.98.  The norm of the "new"
   * vector is kept in nrm0 and updated after operating with each vector
   * 
   */
  void modGramSchmidt(int i, vector<vector<double> > & Hsbg, vector<CSysVector> & w);

  /*! 
   * \brief writes header information for a CSysSolve residual history
   * \param[in,out] os - ostream class object for output
   * \param[in] solver - string describing the solver
   * \param[in] restol - the target tolerance to solve to 
   * \param[in] resinit - the initial residual norm (absolute)
   *
   * \pre the ostream object os should be open
   */
  void writeHeader(ostream & os, const std::string & solver, 
		   const double & restol, const double & resinit);
  
  /*! 
   * \brief writes residual convergence data for one iteration to a stream
   * \param[in,out] os - ostream class object for output
   * \param[in] iter - current iteration
   * \param[in] res - the (absolute) residual norm value
   * \param[in] resinit - the initial residual norm 
   *
   * \pre the ostream object os should be open
   */
  void writeHistory(ostream & os, const int & iter, 
		    const double & res, const double & resinit);
  
public:

  /*! \brief Conjugate Gradient method
   * \param[in] b - the right hand size vector
   * \param[in,out] x - on entry the intial guess, on exit the solution
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] precond - object that defines preconditioner
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   */
  void ConjugateGradient(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                         CPreconditioner & precond, double tol, int m, bool monitoring); 

  /*! 
   * \brief Flexible Generalized Minimal RESidual method
   * \param[in] b - the right hand size vector
   * \param[in,out] x - on entry the intial guess, on exit the solution
   * \param[in] mat_vec - object that defines matrix-vector product
   * \param[in] precond - object that defines preconditioner
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] m - maximum size of the search subspace
   * \param[in] monitoring - turn on priting residuals from solver to screen.
   */
  void FlexibleGMRES(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                     CPreconditioner & precond, double tol, int m, bool monitoring); 

};

#include "linear_solvers_structure.inl"
