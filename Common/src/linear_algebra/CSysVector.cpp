/*!
 * \file CSysVector.cpp
 * \brief Main classes required for solving linear systems of equations
 * \author F. Palacios, J. Hicken
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

#include "../../include/linear_algebra/CSysVector.hpp"
#include "../../include/mpi_structure.hpp"
#include "../../include/toolboxes/allocation_toolbox.hpp"


template<class ScalarType>
CSysVector<ScalarType>::CSysVector(void) {

  vec_val = nullptr;
  nElm = 0;
  nElmDomain = 0;
  nElmGlobal = 0;
  nVar = 0;
  nBlk = 0;
  nBlkDomain = 0;
}

template<class ScalarType>
void CSysVector<ScalarType>::Initialize(unsigned long numBlk, unsigned long numBlkDomain,
                                        unsigned long numVar, const ScalarType* val, bool valIsArray) {

  if ((nElm != numBlk*numVar) && (vec_val != nullptr)) {
    MemoryAllocation::aligned_free(vec_val);
    vec_val = nullptr;
  }

  nBlk = numBlk;
  nElm = numBlk*numVar;
  nBlkDomain = numBlkDomain;
  nElmDomain = numBlkDomain*numVar;
  nVar = numVar;

  SU2_MPI::Allreduce(&nElm, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  if (vec_val == nullptr)
    vec_val = MemoryAllocation::aligned_alloc<ScalarType>(64, nElm*sizeof(ScalarType));

  if(!valIsArray)
    for(auto i=0ul; i<nElm; i++) vec_val[i] = *val;
  else
    for(auto i=0ul; i<nElm; i++) vec_val[i] = val[i];
}

template<class ScalarType>
template<class T>
void CSysVector<ScalarType>::PassiveCopy(const CSysVector<T>& other) {

  /*--- This is a method and not the overload of an operator to make sure who
   calls it knows the consequence to the derivative information (lost) ---*/

  /*--- check if self-assignment, otherwise perform deep copy ---*/
  if ((const void*)this == (const void*)&other) return;

  /*--- determine if (re-)allocation is needed ---*/
  if (nElm != other.GetLocSize() && vec_val != nullptr) {
    MemoryAllocation::aligned_free(vec_val);
    vec_val = nullptr;
  }

  /*--- copy ---*/
  nElm = other.GetLocSize();
  nElmDomain = other.GetNElmDomain();
  nBlk = other.GetNBlk();
  nBlkDomain = other.GetNBlkDomain();
  nVar = other.GetNVar();
  nElmGlobal = other.GetSize();

  if (vec_val == nullptr)
    vec_val = MemoryAllocation::aligned_alloc<ScalarType>(64, nElm*sizeof(ScalarType));

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = SU2_TYPE::GetValue(other[i]);
}

template<class ScalarType>
CSysVector<ScalarType>::~CSysVector() {

  if (vec_val != nullptr)
    MemoryAllocation::aligned_free(vec_val);
}

template<class ScalarType>
void CSysVector<ScalarType>::Equals_AX(ScalarType a, const CSysVector<ScalarType> & x) {

  assert(nElm == x.nElm && "Sizes do not match");

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = a * x.vec_val[i];
}

template<class ScalarType>
void CSysVector<ScalarType>::Plus_AX(ScalarType a, const CSysVector<ScalarType> & x) {

  assert(nElm == x.nElm && "Sizes do not match");

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] += a * x.vec_val[i];
}

template<class ScalarType>
void CSysVector<ScalarType>::Equals_AX_Plus_BY(ScalarType a, const CSysVector<ScalarType> & x,
                                               ScalarType b, const CSysVector<ScalarType> & y) {
  assert(nElm == x.nElm && nElm == y.Elm && "Sizes do not match");

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = a * x.vec_val[i] + b * y.vec_val[i];
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator=(const CSysVector<ScalarType> & u) {
  /*--- check if self-assignment, otherwise perform deep copy ---*/
  if (this != &u)
    Initialize(u.nBlk, u.nBlkDomain, u.nVar, u.vec_val, true);
  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator=(ScalarType val) {
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = val;
  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator+=(const CSysVector<ScalarType> & u) {

  assert(nElm == u.nElm && "Sizes do not match");

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] += u.vec_val[i];
  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator-=(const CSysVector<ScalarType> & u) {

  assert(nElm == u.nElm && "Sizes do not match");

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] -= u.vec_val[i];
  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator*=(ScalarType val) {

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] *= val;
  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator/=(ScalarType val) {

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] /= val;
  return *this;
}

template<class ScalarType>
void CSysVector<ScalarType>::CopyToArray(ScalarType* u_array) const {

  for (unsigned long i = 0; i < nElm; i++)
    u_array[i] = vec_val[i];
}

template<class ScalarType>
ScalarType dotProd(const CSysVector<ScalarType> & u, const CSysVector<ScalarType> & v) {

  assert(u.nElm == v.nElm && "Sizes do not match");

  /*--- find local inner product and, if a parallel run, sum over all
   processors (we use nElemDomain instead of nElem) ---*/
  ScalarType loc_prod = 0.0;
  for (unsigned long i = 0; i < u.nElmDomain; i++)
    loc_prod += u.vec_val[i]*v.vec_val[i];
  ScalarType prod = 0.0;

  SelectMPIWrapper<ScalarType>::W::Allreduce(&loc_prod, &prod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return prod;
}

/*--- Explicit instantiations ---*/
template class CSysVector<su2double>;
template void CSysVector<su2double>::PassiveCopy(const CSysVector<su2double>&);
template su2double dotProd<su2double>(const CSysVector<su2double> & u, const CSysVector<su2double> & v);

#ifdef CODI_REVERSE_TYPE
template class CSysVector<passivedouble>;
template void CSysVector<su2double>::PassiveCopy(const CSysVector<passivedouble>&);
template void CSysVector<passivedouble>::PassiveCopy(const CSysVector<su2double>&);
template passivedouble dotProd<passivedouble>(const CSysVector<passivedouble> & u, const CSysVector<passivedouble> & v);
#endif
