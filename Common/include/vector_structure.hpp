/*!
 * \file vector_structure.hpp
 * \brief Headers for the classes related to linear solvers (CG, FGMRES, etc)
 *        The subroutines and functions are in the <i>linear_solvers_structure.cpp</i> file.
 * \author F. Palacios, J. Hicken, T. Economon
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

#include <climits>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>

using namespace std;

/*--- Forward declaration of template friend functions. ---*/
template<class T> class CSysVector;
template<class T> T dotProd(const CSysVector<T> & u, const CSysVector<T> & v);

/*!
 * \class CSysVector
 * \brief Class for holding and manipulating vectors needed by linear solvers
 * \author J. Hicken.
 *
 * We could use the STL vector as a base class here, but this gives us
 * more flexibility with the underlying data (e.g. we may decide to
 * use a block storage scheme rather than a continuous storage
 * scheme).
 */
template<class ScalarType>
class CSysVector {
  
private:
	unsigned long nElm; /*!< \brief total number of elements (or number elements on this processor) */
	unsigned long nElmDomain; /*!< \brief total number of elements (or number elements on this processor without Ghost cells) */
#ifdef HAVE_MPI
  unsigned long nElmGlobal; /*!< \brief total number of elements over all processors */
#endif
	unsigned short nVar; /*!< \brief number of elements in a block */
	unsigned long nBlk; /*!< \brief number of blocks (or number of blocks on this processor) */
	unsigned long nBlkDomain; /*!< \brief number of blocks (or number of blocks on this processor without Ghost cells) */
  ScalarType* vec_val; /*!< \brief storage for the element values */
  
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
  CSysVector(const unsigned long & size, const ScalarType & val = 0.0);
  
  /*!
   * \brief constructor of the class.
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain
   * \param[in] numVar - number of variables in each block
   * \param[in] val - default value for elements
   */
  CSysVector(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const ScalarType & val = 0.0);
  
  /*!
   * \brief copy constructor of the class.
   * \param[in] u - CSysVector that is being copied
   */
  CSysVector(const CSysVector & u);
  
  /*!
	 * \brief Sets to zero all the entries of the vector.
	 */
	void SetValZero(void);
  
  /*!
   * \brief constructor from array
   * \param[in] size - number of elements locally
   * \param[in] u_array - vector stored as array being copied
   */
  explicit CSysVector(const unsigned long & size, const ScalarType* u_array);
  
  /*!
   * \brief constructor from array
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain - number of blocks locally (without g cells)
   * \param[in] numVar - number of variables in each block
   * \param[in] u_array - vector stored as array being copied
   */
  explicit CSysVector(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar,
                      const ScalarType* u_array);
  
  /*!
   * \brief class destructor
   */
  virtual ~CSysVector();
  
  /*!
   * \brief Initialize the class.
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain
   * \param[in] numVar - number of variables in each block
   * \param[in] val - default value for elements
   */
  void Initialize(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const ScalarType & val = 0.0);
  
  /*!
   * \brief return the number of local elements in the CSysVector
   */
  unsigned long GetLocSize() const;
  
  /*!
   * \brief return the number of local elements in the CSysVector without ghost cells
   */
  unsigned long GetNElmDomain() const;
  
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
  unsigned long GetNBlk() const;
	
	/*!
   * \brief return the number of blocks (typically number of nodes locally)
   */
  unsigned long GetNBlkDomain() const;
  
  /*!
   * \brief set calling CSysVector to scaling of another CSysVector
   * \param[in] a - scalar factor for x
   * \param[in] x - CSysVector that is being scaled
   */
  void Equals_AX(const ScalarType & a, CSysVector & x);
  
  /*!
   * \brief adds a scaled CSysVector to calling CSysVector
   * \param[in] a - scalar factor for x
   * \param[in] x - CSysVector that is being scaled
   */
  void Plus_AX(const ScalarType & a, CSysVector & x);
  
  /*!
   * \brief general linear combination of two CSysVectors
   * \param[in] a - scalar factor for x
   * \param[in] x - first CSysVector in linear combination
   * \param[in] b - scalar factor for y
   * \param[in] y - second CSysVector in linear combination
   */
  void Equals_AX_Plus_BY(const ScalarType & a, CSysVector & x, const ScalarType & b, CSysVector & y);
  
  /*!
   * \brief assignment operator with deep copy
   * \param[in] u - CSysVector whose values are being assigned
   */
  CSysVector & operator=(const CSysVector & u);
  
  /*!
   * \brief CSysVector=su2double assignment operator
   * \param[in] val - value assigned to each element of CSysVector
   */
  CSysVector & operator=(const ScalarType & val);
  
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
  CSysVector operator*(const ScalarType & val) const;
  
  /*!
   * \brief compound scalar multiplication-assignment operator
   * \param[in] val - value to multiply calling object by
   */
  CSysVector & operator*=(const ScalarType & val);
  
  /*!
   * \brief vector-scalar division operator (no scalar/vector operator)
   * \param[in] val - value to divide elements of *this by
   */
  CSysVector operator/(const ScalarType & val) const;
  
  /*!
   * \brief compound scalar division-assignment operator
   * \param[in] val - value to divide elements of calling object by
   */
  CSysVector & operator/=(const ScalarType & val);
  
  /*!
   * \brief indexing operator with assignment permitted
   * \param[in] i = local index to access
   */
  ScalarType & operator[](const unsigned long & i);
  
  /*!
   * \brief indexing operator with assignment not permitted
   * \param[in] i = local index to access
   */
  const ScalarType & operator[](const unsigned long & i) const;
    
  /*!
   * \brief the L2 norm of the CSysVector
   * \result the L2 norm
   */
  ScalarType norm() const;
  
  /*!
   * \brief copies the contents of the calling CSysVector into an array
   * \param[out] u_array - array into which information is being copied
   * \pre u_array must be allocated and have the same size as CSysVector
   */
  void CopyToArray(ScalarType* u_array);
  
  /*!
	 * \brief Subtract val_residual to the residual.
	 * \param[in] val_ipoint - index of the point where subtract the residual.
   * \param[in] val_residual - Value to subtract to the residual.
	 */
  void SubtractBlock(unsigned long val_ipoint, ScalarType *val_residual);
  
  /*!
	 * \brief Add val_residual to the residual.
	 * \param[in] val_ipoint - index of the point where add the residual.
   * \param[in] val_residual - Value to add to the residual.
	 */
  void AddBlock(unsigned long val_ipoint, ScalarType *val_residual);
  
  /*!
	 * \brief Set val_residual to the residual.
	 * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_var - inde of the residual to be set.
   * \param[in] val_residual - Value to set to the residual.
	 */
  void SetBlock(unsigned long val_ipoint, unsigned short val_var, ScalarType val_residual);
  
  /*!
	 * \brief Set val_residual to the residual.
	 * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_residual - Value to set to the residual.
	 */
  void SetBlock(unsigned long val_ipoint, ScalarType *val_residual);
  
  /*!
	 * \brief Set the residual to zero.
	 * \param[in] val_ipoint - index of the point where set the residual.
	 */
  void SetBlock_Zero(unsigned long val_ipoint);
  
  /*!
	 * \brief Set the velocity residual to zero.
	 * \param[in] val_ipoint - index of the point where set the residual.
	 * \param[in] val_var - inde of the residual to be set.
	 */
  void SetBlock_Zero(unsigned long val_ipoint, unsigned short val_var);
	
  /*!
	 * \brief Get the value of the residual.
	 * \param[in] val_ipoint - index of the point where set the residual.
   * \return Pointer to the residual.
	 */
  ScalarType *GetBlock(unsigned long val_ipoint);
	
  /*!
	 * \brief Get the value of the residual.
	 * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_var - inde of the residual to be set.
   * \return Value of the residual.
	 */
  ScalarType GetBlock(unsigned long val_ipoint, unsigned short val_var);
  
  /*!
   * \brief dot-product between two CSysVectors
   * \param[in] u - first CSysVector in dot product
   * \param[in] v - second CSysVector in dot product
   */
  friend ScalarType dotProd<ScalarType>(const CSysVector & u, const CSysVector & v);
  
  /*!
   * \brief Set our values (resizing if required) by copying from other, the derivative information is lost.
   * \param[in] other - source CSysVector
   */
  template<class T>
  void PassiveCopy(const CSysVector<T>& other);
};

/*!
 * \brief scalar * vector multiplication operator
 * \param[in] val - scalar value to multiply by
 * \param[in] u - CSysVector having its elements scaled
 */
template<class ScalarType>
CSysVector<ScalarType> operator*(const ScalarType & val, const CSysVector<ScalarType> & u);

/*!
 * \class CMatrixVectorProduct
 * \brief abstract base class for defining matrix-vector products
 * \author J. Hicken.
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
template<class ScalarType>
class CMatrixVectorProduct {
public:
  virtual ~CMatrixVectorProduct() = 0; ///< class destructor
  virtual void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v)
  const = 0; ///< matrix-vector product operation
};
template<class ScalarType>
inline CMatrixVectorProduct<ScalarType>::~CMatrixVectorProduct() {}

/*!
 * \class CPreconditioner
 * \brief abstract base class for defining preconditioning operation
 * \author J. Hicken.
 *
 * See the remarks regarding the CMatrixVectorProduct class.  The same
 * idea applies here to the preconditioning operation.
 */
template<class ScalarType>
class CPreconditioner {
public:
  virtual ~CPreconditioner() = 0; ///< class destructor
  virtual void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v)
  const = 0; ///< preconditioning operation
};
template<class ScalarType>
inline CPreconditioner<ScalarType>::~CPreconditioner() {}

#include "vector_structure.inl"
