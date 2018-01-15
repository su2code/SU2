/*!
 * \file vector_structure.hpp
 * \brief Headers for the classes related to linear solvers (CG, FGMRES, etc)
 *        The subroutines and functions are in the <i>linear_solvers_structure.cpp</i> file.
 * \author F. Palacios, J. Hicken, T. Economon
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

#include <climits>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>

using namespace std;

/*!
 * \class TCSysVector
 * \brief Class for holding and manipulating vectors needed by linear solvers
 * \author J. Hicken.
 * \version 5.0.0 "Raven"
 *
 * We could use the STL vector as a base class here, but this gives us
 * more flexibility with the underlying data (e.g. we may decide to
 * use a block storage scheme rather than a continuous storage
 * scheme).
 */
template<class CalcType>
class TCSysVector {
  
private:
	unsigned long nElm; /*!< \brief total number of elements (or number elements on this processor) */
	unsigned long nElmDomain; /*!< \brief total number of elements (or number elements on this processor without Ghost cells) */
#ifdef HAVE_MPI
  unsigned long nElmGlobal; /*!< \brief total number of elements over all processors */
#endif
	unsigned short nVar; /*!< \brief number of elements in a block */
	unsigned long nBlk; /*!< \brief number of blocks (or number of blocks on this processor) */
	unsigned long nBlkDomain; /*!< \brief number of blocks (or number of blocks on this processor without Ghost cells) */
  
public:
  CalcType* vec_val; /*!< \brief storage for the element values */
  
  /*!
   * \brief default constructor of the class.
   */
  TCSysVector(void);
  
  /*!
   * \brief constructor of the class.
   * \param[in] size - number of elements locally
   * \param[in] val - default value for elements
   */
  TCSysVector(const unsigned long & size, const CalcType & val = 0.0);
  
  /*!
   * \brief constructor of the class.
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain
   * \param[in] numVar - number of variables in each block
   * \param[in] val - default value for elements
   */
  TCSysVector(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const CalcType & val = 0.0);
  
  /*!
   * \brief copy constructor of the class.
   * \param[in] u - CSysVector that is being copied
   */
  TCSysVector(const TCSysVector<CalcType> & u);
  
  /*!
	 * \brief Sets to zero all the entries of the vector.
	 */
	void SetValZero(void);
  
  /*!
   * \brief constructor from array
   * \param[in] size - number of elements locally
   * \param[in] u_array - vector stored as array being copied
   */
  explicit TCSysVector(const unsigned long & size, const CalcType* u_array);
  
  /*!
   * \brief constructor from array
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain - number of blocks locally (without g cells)
   * \param[in] numVar - number of variables in each block
   * \param[in] u_array - vector stored as array being copied
   */
  explicit TCSysVector(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar,
                      const CalcType* u_array);
  
  /*!
   * \brief class destructor
   */
  virtual ~TCSysVector();
  
  /*!
   * \brief Initialize the class.
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain
   * \param[in] numVar - number of variables in each block
   * \param[in] val - default value for elements
   */
  void Initialize(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const CalcType & val = 0.0);
  
  /*!
   * \brief return the number of local elements in the CSysVector
   */
  unsigned long GetLocSize() const;
  
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
  void Equals_AX(const CalcType & a, TCSysVector<CalcType> & x);
  
  /*!
   * \brief adds a scaled CSysVector to calling CSysVector
   * \param[in] a - scalar factor for x
   * \param[in] x - CSysVector that is being scaled
   */
  void Plus_AX(const CalcType & a, TCSysVector<CalcType> & x);
  
  /*!
   * \brief general linear combination of two CSysVectors
   * \param[in] a - scalar factor for x
   * \param[in] x - first CSysVector in linear combination
   * \param[in] b - scalar factor for y
   * \param[in] y - second CSysVector in linear combination
   */
  void Equals_AX_Plus_BY(const CalcType & a, TCSysVector<CalcType> & x, const CalcType & b, TCSysVector<CalcType> & y);
  
  /*!
   * \brief assignment operator with deep copy
   * \param[in] u - CSysVector whose values are being assigned
   */
  TCSysVector<CalcType>& operator=(const TCSysVector<CalcType> & u);
  
  /*!
   * \brief CSysVector=CalcType assignment operator
   * \param[in] val - value assigned to each element of CSysVector
   */
  TCSysVector<CalcType> & operator=(const CalcType & val);
  
  /*!
   * \brief addition operator
   * \param[in] u - CSysVector being added to *this
   */
  TCSysVector<CalcType> operator+(const TCSysVector<CalcType> & u) const;
  
  /*!
   * \brief compound addition-assignment operator
   * \param[in] u - CSysVector being added to calling object
   */
  TCSysVector<CalcType> & operator+=(const TCSysVector<CalcType> & u);
  
  /*!
   * \brief subtraction operator
   * \param[in] u - CSysVector being subtracted from *this
   */
  TCSysVector<CalcType> operator-(const TCSysVector<CalcType> & u) const;
  
  /*!
   * \brief compound subtraction-assignment operator
   * \param[in] u - CSysVector being subtracted from calling object
   */
  TCSysVector<CalcType> & operator-=(const TCSysVector<CalcType> & u);
  
  /*!
   * \brief vector * scalar multiplication operator
   * \param[in] val - value to multiply *this by
   */
  TCSysVector<CalcType> operator*(const CalcType & val) const;
  
  /*!
   * \brief scalar * vector multiplication operator
   * \param[in] val - scalar value to multiply by
   * \param[in] u - CSysVector having its elements scaled
   */
  template<class T>
  friend TCSysVector<T> operator*(const CalcType & val, const TCSysVector<T> & u);
  
  /*!
   * \brief compound scalar multiplication-assignment operator
   * \param[in] val - value to multiply calling object by
   */
  TCSysVector<CalcType> & operator*=(const CalcType & val);
  
  /*!
   * \brief vector-scalar division operator (no scalar/vector operator)
   * \param[in] val - value to divide elements of *this by
   */
  TCSysVector<CalcType> operator/(const CalcType & val) const;
  
  /*!
   * \brief compound scalar division-assignment operator
   * \param[in] val - value to divide elements of calling object by
   */
  TCSysVector<CalcType>& operator/=(const CalcType & val);
  
  /*!
   * \brief indexing operator with assignment permitted
   * \param[in] i = local index to access
   */
  CalcType & operator[](const unsigned long & i);
  
  /*!
   * \brief indexing operator with assignment not permitted
   * \param[in] i = local index to access
   */
  const CalcType & operator[](const unsigned long & i) const;
    
  /*!
   * \brief the L2 norm of the CSysVector
   * \result the L2 norm
   */
  CalcType norm() const;
  
  /*!
   * \brief copies the contents of the calling CSysVector into an array
   * \param[out] u_array - array into which information is being copied
   * \pre u_array must be allocated and have the same size as CSysVector
   */
  void CopyToArray(CalcType* u_array);
  
  /*!
	 * \brief Subtract val_residual to the residual.
	 * \param[in] val_ipoint - index of the point where subtract the residual.
   * \param[in] val_residual - Value to subtract to the residual.
	 */
  void SubtractBlock(unsigned long val_ipoint, CalcType *val_residual);
  
  /*!
	 * \brief Add val_residual to the residual.
	 * \param[in] val_ipoint - index of the point where add the residual.
   * \param[in] val_residual - Value to add to the residual.
	 */
  void AddBlock(unsigned long val_ipoint, CalcType *val_residual);
  
  /*!
	 * \brief Set val_residual to the residual.
	 * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_var - inde of the residual to be set.
   * \param[in] val_residual - Value to set to the residual.
	 */
  void SetBlock(unsigned long val_ipoint, unsigned short val_var, CalcType val_residual);
  
  /*!
	 * \brief Set val_residual to the residual.
	 * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_residual - Value to set to the residual.
	 */
  void SetBlock(unsigned long val_ipoint, CalcType *val_residual);
  
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
  CalcType *GetBlock(unsigned long val_ipoint);
	
  /*!
	 * \brief Get the value of the residual.
	 * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_var - inde of the residual to be set.
   * \return Value of the residual.
	 */
  CalcType GetBlock(unsigned long val_ipoint, unsigned short val_var);
  
  
  /*!
   * \brief dot-product between two CSysVectors
   * \param[in] u - first CSysVector in dot product
   * \param[in] v - second CSysVector in dot product
   */
  template<class T>
  friend T dotProd(const TCSysVector<T> & u, const TCSysVector<T> & v);
  
};

typedef TCSysVector<su2double> CSysVector;

/*!
 * \class TCMatrixVectorProduct
 * \brief abstract base class for defining matrix-vector products
 * \author J. Hicken.
 * \version 5.0.0 "Raven"
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
template<class CalcType>
class TCMatrixVectorProduct {
public:
  virtual ~TCMatrixVectorProduct() = 0; ///< class destructor
  virtual void operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v)
  const = 0; ///< matrix-vector product operation
};
template<class CalcType>
inline TCMatrixVectorProduct<CalcType>::~TCMatrixVectorProduct() {}

typedef TCMatrixVectorProduct<su2double> CMatrixVectorProduct;

/*!
 * \class TCPreconditioner
 * \brief abstract base class for defining preconditioning operation
 * \author J. Hicken.
 * \version 5.0.0 "Raven"
 *
 * See the remarks regarding the CMatrixVectorProduct class.  The same
 * idea applies here to the preconditioning operation.
 */
template<class CalcType>
class TCPreconditioner {
public:
  virtual ~TCPreconditioner() = 0; ///< class destructor
  virtual void operator()(const TCSysVector<CalcType> & u, TCSysVector<CalcType> & v)
  const = 0; ///< preconditioning operation
};
template<class CalcType>
inline TCPreconditioner<CalcType>::~TCPreconditioner() {}

typedef TCPreconditioner<su2double> CPreconditioner;

template class TCSysVector<su2double>;
#ifdef CODI_REVERSE_TYPE
template class TCSysVector<passivedouble>;
#endif
#include "vector_structure.inl"
