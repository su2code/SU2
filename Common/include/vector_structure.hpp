/*!
 * \file vector_structure.hpp
 * \brief Headers for the classes related to linear solvers (CG, FGMRES, etc)
 *        The subroutines and functions are in the <i>linear_solvers_structure.cpp</i> file.
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

#ifndef NO_MPI
#include <mpi.h>
#endif
#include <climits>
#include <limits>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

const double eps = numeric_limits<double>::epsilon(); /*!< \brief machine epsilon */

/*!
 * \class CSysVector
 * \brief Class for holding and manipulating vectors needed by linear solvers
 * \author J. Hicken.
 * \version 2.0.6
 *
 * We could use the STL vector as a base class here, but this gives us
 * more flexibility with the underlying data (e.g. we may decide to
 * use a block storage scheme rather than a continuous storage
 * scheme).
 */
class CSysVector {
  
private:
  
	unsigned int nElm; /*!< \brief total number of elements (or number elements on this processor) */
	unsigned int nElmDomain; /*!< \brief total number of elements (or number elements on this processor without Ghost cells) */
  unsigned long nElmGlobal; /*!< \brief total number of elements over all processors */
	unsigned short nVar; /*!< \brief number of elements in a block */
	unsigned int nBlk; /*!< \brief number of blocks (or number of blocks on this processor) */
	unsigned int nBlkDomain; /*!< \brief number of blocks (or number of blocks on this processor without Ghost cells) */
  unsigned int myrank; /*!< \brief processor rank (only used for parallel runs) */
  double* vec_val; /*!< \brief storage for the element values */
  
public:
  
  /*!
   * \brief default constructor of the class.
   */
  CSysVector(void);
  
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
  CSysVector(const unsigned int & numBlk, const unsigned int & numBlkDomain, const unsigned short & numVar, const double & val = 0.0);
  
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
   * \param[in] numBlkDomain - number of blocks locally (without g cells)
   * \param[in] numVar - number of variables in each block
   * \param[in] u_array - vector stored as array being copied
   */
  explicit CSysVector(const unsigned int & numBlk, const unsigned int & numBlkDomain, const unsigned short & numVar,
                      const double* u_array);
  
  /*!
   * \brief class destructor
   */
  virtual ~CSysVector();
  
  /*!
   * \brief Initialize the class.
   * \param[in] numBlk - number of blocks locally
   * \param[in] numVar - number of variables in each block
   * \param[in] val - default value for elements
   */
  void Initialize(const unsigned int & numBlk, const unsigned int & numBlkDomain, const unsigned short & numVar, const double & val = 0.0);
  
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
   * \brief return the number of blocks (typically number of nodes locally)
   */
  unsigned int GetNBlkDomain() const;
  
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
	 * \brief Subtract val_residual to the residual.
	 * \param[in] val_ipoint - index of the point where subtract the residual.
   * \param[in] val_residual - Value to subtract to the residual.
	 */
  void SubtractBlock(unsigned long val_ipoint, double *val_residual);
  
  /*!
	 * \brief Add val_residual to the residual.
	 * \param[in] val_ipoint - index of the point where add the residual.
   * \param[in] val_residual - Value to add to the residual.
	 */
  void AddBlock(unsigned long val_ipoint, double *val_residual);
  
  /*!
	 * \brief Set val_residual to the residual.
	 * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_var - inde of the residual to be set.
   * \param[in] val_residual - Value to set to the residual.
	 */
  void SetBlock(unsigned long val_ipoint, unsigned short val_var, double val_residual);
  
  /*!
	 * \brief Set val_residual to the residual.
	 * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_residual - Value to set to the residual.
	 */
  void SetBlock(unsigned long val_ipoint, double *val_residual);
  
  /*!
	 * \brief Set the residual to zero.
	 * \param[in] val_ipoint - index of the point where set the residual.
	 */
  void SetBlock_Zero(unsigned long val_ipoint);
  
  /*!
	 * \brief Set the velocity residual to zero.
	 * \param[in] val_ipoint - index of the point where set the residual.
	 */
  void SetBlock_Zero(unsigned long val_ipoint, unsigned short val_var);
	
  /*!
	 * \brief Get the value of the residual.
	 * \param[in] val_ipoint - index of the point where set the residual.
   * \return Pointer to the residual.
	 */
  double *GetBlock(unsigned long val_ipoint);
	
  /*!
	 * \brief Get the value of the residual.
	 * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_var - inde of the residual to be set.
   * \return Value of the residual.
	 */
  double GetBlock(unsigned long val_ipoint, unsigned short val_var);
  
  
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
 * \version 2.0.6
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
 * \version 2.0.6
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

#include "vector_structure.inl"
