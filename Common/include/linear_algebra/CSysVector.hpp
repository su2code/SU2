/*!
 * \file CSysVector.hpp
 * \brief Declararion of the vector class used in the solution of
 *        large, distributed, sparse linear systems.
 * \author F. Palacios, J. Hicken, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include <cmath>
#include <cstdlib>


/*!
 * \class CSysVector
 * \brief Class for holding and manipulating vectors needed by linear solvers
 * \author J. Hicken.
 *
 * We could use the STL vector as a base class here, but this gives us
 * more flexibility with the underlying data (e.g. we may decide to
 * use a block storage scheme rather than a continuous storage scheme).
 */
template<class ScalarType>
class CSysVector {

private:
  enum { OMP_MAX_SIZE = 4096 };   /*!< \brief Maximum chunk size used in parallel for loops. */

  unsigned long omp_chunk_size;   /*!< \brief Static chunk size used in loop, determined at initialization. */
  ScalarType* vec_val;            /*!< \brief storage for the element values, 64 byte aligned (do not use normal new/delete) */
  unsigned long nElm;             /*!< \brief total number of elements (or number elements on this processor) */
  unsigned long nElmDomain;       /*!< \brief total number of elements (or number elements on this processor without Ghost cells) */
  unsigned long nVar;             /*!< \brief number of elements in a block */
  mutable ScalarType dotRes;      /*!< \brief result of dot product. to perform a reduction with OpenMP the
                                              variable needs to be declared outside the parallel region */

  /*!
   * \brief Generic initialization from a scalar or array.
   * \note If val==nullptr vec_val is not initialized, only allocated.
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain - number of blocks locally (without g cells)
   * \param[in] numVar - number of variables in each block
   * \param[in] val - default value for elements
   * \param[in] valIsArray - if true val is treated as array
   */
  void Initialize(unsigned long numBlk, unsigned long numBlkDomain, unsigned long numVar,
                  const ScalarType* val, bool valIsArray);

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
  CSysVector(unsigned long size, ScalarType val = 0.0) {
    nElm = 0; vec_val = nullptr;
    Initialize(size, size, 1, &val, false);
  }

  /*!
   * \brief constructor of the class.
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain - number of blocks locally (without g cells)
   * \param[in] numVar - number of variables in each block
   * \param[in] val - default value for elements
   */
  CSysVector(unsigned long numBlk, unsigned long numBlkDomain, unsigned long numVar, ScalarType val = 0.0) {
    nElm = 0; vec_val = nullptr;
    Initialize(numBlk, numBlkDomain, numVar, &val, false);
  }

  /*!
   * \brief constructor from array
   * \param[in] size - number of elements locally
   * \param[in] u_array - vector stored as array being copied
   */
  explicit CSysVector(unsigned long size, const ScalarType* u_array) {
    nElm = 0; vec_val = nullptr;
    Initialize(size, size, 1, u_array, true);
  }

  /*!
   * \brief constructor from array
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain - number of blocks locally (without g cells)
   * \param[in] numVar - number of variables in each block
   * \param[in] u_array - vector stored as array being copied
   */
  explicit CSysVector(unsigned long numBlk, unsigned long numBlkDomain, unsigned long numVar, const ScalarType* u_array) {
    nElm = 0; vec_val = nullptr;
    Initialize(numBlk, numBlkDomain, numVar, u_array, true);
  }

  /*!
   * \brief copy constructor of the class.
   * \param[in] u - CSysVector that is being copied
   */
  CSysVector(const CSysVector & u) {
    nElm = 0; vec_val = nullptr;
    Initialize(u.GetNBlk(), u.GetNBlkDomain(), u.nVar, u.vec_val, true);
  }

  /*!
   * \brief Set our values (resizing if required) by copying from other, the derivative information is lost.
   * \param[in] other - source CSysVector
   */
  template<class T>
  void PassiveCopy(const CSysVector<T>& other);

  /*!
   * \brief class destructor
   */
  ~CSysVector();

  /*!
   * \brief Initialize the class with a scalar.
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain - number of blocks locally (without g cells)
   * \param[in] numVar - number of variables in each block
   * \param[in] val - default value for elements
   */
  void Initialize(unsigned long numBlk, unsigned long numBlkDomain, unsigned long numVar, ScalarType val = 0.0) {
    Initialize(numBlk, numBlkDomain, numVar, &val, false);
  }

  /*!
   * \brief Initialize the class with an array.
   * \note If ptr==nullptr no copy occurs.
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain - number of blocks locally (without g cells)
   * \param[in] numVar - number of variables in each block
   * \param[in] ptr - pointer to data with which to initialize the vector
   */
  void Initialize(unsigned long numBlk, unsigned long numBlkDomain, unsigned long numVar, const ScalarType* ptr) {
    Initialize(numBlk, numBlkDomain, numVar, ptr, true);
  }

  /*!
   * \brief return the number of local elements in the CSysVector
   */
  inline unsigned long GetLocSize() const { return nElm; }

  /*!
   * \brief return the number of local elements in the CSysVector without ghost cells
   */
  inline unsigned long GetNElmDomain() const { return nElmDomain; }

  /*!
   * \brief return the number of variables at each block (typically number per node)
   */
  inline unsigned long GetNVar() const { return nVar; }

  /*!
   * \brief return the number of blocks (typically number of nodes locally)
   */
  inline unsigned long GetNBlk() const { return nElm/nVar; }

  /*!
   * \brief return the number of blocks (typically number of nodes locally)
   */
  inline unsigned long GetNBlkDomain() const { return nElmDomain/nVar; }

  /*!
   * \brief set calling CSysVector to scaling of another CSysVector
   * \param[in] a - scalar factor for x
   * \param[in] x - CSysVector that is being scaled
   */
  void Equals_AX(ScalarType a, const CSysVector & x);

  /*!
   * \brief adds a scaled CSysVector to calling CSysVector
   * \param[in] a - scalar factor for x
   * \param[in] x - CSysVector that is being scaled
   */
  void Plus_AX(ScalarType a, const CSysVector & x);

  /*!
   * \brief general linear combination of two CSysVectors
   * \param[in] a - scalar factor for x
   * \param[in] x - first CSysVector in linear combination
   * \param[in] b - scalar factor for y
   * \param[in] y - second CSysVector in linear combination
   */
  void Equals_AX_Plus_BY(ScalarType a, const CSysVector & x, ScalarType b, const CSysVector & y);

  /*!
   * \brief assignment operator with deep copy
   * \param[in] u - CSysVector whose values are being assigned
   */
  CSysVector & operator=(const CSysVector & u);

  /*!
   * \brief CSysVector=su2double assignment operator
   * \param[in] val - value assigned to each element of CSysVector
   */
  CSysVector & operator=(ScalarType val);

  /*!
   * \brief Sets to zero all the entries of the vector.
   */
  inline void SetValZero(void) { *this = ScalarType(0.0); }

  /*!
   * \brief compound addition-assignment operator
   * \param[in] u - CSysVector being added to calling object
   */
  CSysVector & operator+=(const CSysVector & u);

  /*!
   * \brief compound subtraction-assignment operator
   * \param[in] u - CSysVector being subtracted from calling object
   */
  CSysVector & operator-=(const CSysVector & u);

  /*!
   * \brief compound scalar multiplication-assignment operator
   * \param[in] val - value to multiply calling object by
   */
  CSysVector & operator*=(ScalarType val);

  /*!
   * \brief compound scalar division-assignment operator
   * \param[in] val - value to divide elements of calling object by
   */
  CSysVector & operator/=(ScalarType val);

  /*!
   * \brief Dot product between "this" and another vector
   * \param[in] u - Another vector.
   * \return result of dot product
   */
  ScalarType dot(const CSysVector & u) const;

  /*!
   * \brief squared L2 norm of the vector (via dot with self)
   * \return squared L2 norm
   */
  inline ScalarType squaredNorm() const { return dot(*this); }

  /*!
   * \brief L2 norm of the vector
   * \return L2 norm
   */
  inline ScalarType norm() const { return sqrt(squaredNorm()); }

  /*!
   * \brief indexing operator with assignment permitted
   * \param[in] i = local index to access
   */
  inline ScalarType & operator[](const unsigned long & i) { return vec_val[i]; }

  /*!
   * \brief indexing operator with assignment not permitted
   * \param[in] i = local index to access
   */
  inline const ScalarType & operator[](const unsigned long & i) const { return vec_val[i]; }

  /*!
   * \brief copies the contents of the calling CSysVector into an array
   * \param[out] u_array - array into which information is being copied
   * \pre u_array must be allocated and have the same size as CSysVector
   */
  void CopyToArray(ScalarType* u_array) const;

  /*!
   * \brief Subtract val_residual to the residual.
   * \param[in] val_ipoint - index of the point where subtract the residual.
   * \param[in] val_residual - Value to subtract to the residual.
   */
  inline void SubtractBlock(unsigned long val_ipoint, const ScalarType *val_residual) {
    for (auto iVar = 0ul; iVar < nVar; iVar++)
      vec_val[val_ipoint*nVar+iVar] -= val_residual[iVar];
  }

  /*!
   * \brief Add val_residual to the residual.
   * \param[in] val_ipoint - index of the point where add the residual.
   * \param[in] val_residual - Value to add to the residual.
   */
  inline void AddBlock(unsigned long val_ipoint, const ScalarType *val_residual) {
    for (auto iVar = 0ul; iVar < nVar; iVar++)
      vec_val[val_ipoint*nVar+iVar] += val_residual[iVar];
  }

  /*!
   * \brief Set val_residual to the residual.
   * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_var - inde of the residual to be set.
   * \param[in] val_residual - Value to set to the residual.
   */
  inline void SetBlock(unsigned long val_ipoint, unsigned long val_var, ScalarType val_residual) {
    vec_val[val_ipoint*nVar+val_var] = val_residual;
  }

  /*!
   * \brief Set val_residual to the residual.
   * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_residual - Value to set to the residual.
   */
  inline void SetBlock(unsigned long val_ipoint, const ScalarType *val_residual) {
    for (auto iVar = 0ul; iVar < nVar; iVar++)
      vec_val[val_ipoint*nVar+iVar] = val_residual[iVar];
  }

  /*!
   * \brief Set the residual to zero.
   * \param[in] val_ipoint - index of the point where set the residual.
   */
  inline void SetBlock_Zero(unsigned long val_ipoint) {
    for (auto iVar = 0ul; iVar < nVar; iVar++)
      vec_val[val_ipoint*nVar+iVar] = 0.0;
  }

  /*!
   * \brief Set the velocity residual to zero.
   * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_var - inde of the residual to be set.
   */
  inline void SetBlock_Zero(unsigned long val_ipoint, unsigned long val_var) {
    vec_val[val_ipoint*nVar+val_var] = 0.0;
  }

  /*!
   * \brief Get the value of the residual.
   * \param[in] val_ipoint - index of the point where set the residual.
   * \return Pointer to the residual.
   */
  inline ScalarType *GetBlock(unsigned long val_ipoint) { return &vec_val[val_ipoint*nVar]; }

  /*!
   * \brief Get the value of the residual.
   * \param[in] val_ipoint - index of the point where set the residual.
   * \param[in] val_var - inde of the residual to be set.
   * \return Value of the residual.
   */
  inline ScalarType GetBlock(unsigned long val_ipoint, unsigned long val_var) const {
    return vec_val[val_ipoint*nVar+val_var];
  }
};
