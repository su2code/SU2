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
#include "../../include/omp_structure.hpp"
#include "../../include/toolboxes/allocation_toolbox.hpp"

/*!
 * \brief OpenMP worksharing construct used in CSysVector for loops.
 * \note The loop will only run in parallel if methods are called from a
 * parallel region (if not the results will still be correct).
 * Static schedule to reduce overhead, chunk size determined at initialization.
 * "nowait" clause is safe when calling CSysVector methods after each other
 * as the loop size is the same. Methods of other classes that operate on a
 * CSysVector and do not have the same work scheduling must use a
 * SU2_OMP_BARRIER before using the vector.
 */
#define PARALLEL_FOR SU2_OMP(for schedule(static,omp_chunk_size) nowait)

template<class ScalarType>
CSysVector<ScalarType>::CSysVector(void) {

  vec_val = nullptr;
  nElm = 0;
  nElmDomain = 0;
  nVar = 0;
  omp_chunk_size = OMP_MAX_SIZE;
  dotRes = 0.0;
}

template<class ScalarType>
void CSysVector<ScalarType>::Initialize(unsigned long numBlk, unsigned long numBlkDomain,
                                        unsigned long numVar, const ScalarType* val, bool valIsArray) {

  /*--- Assert that this method is only called by one thread. ---*/
  assert(omp_get_thread_num()==0 && "Only the master thread is allowed to initialize the vector.");

  if ((nElm != numBlk*numVar) && (vec_val != nullptr)) {
    MemoryAllocation::aligned_free(vec_val);
    vec_val = nullptr;
  }

  nElm = numBlk*numVar;
  nElmDomain = numBlkDomain*numVar;
  nVar = numVar;

  omp_chunk_size = computeStaticChunkSize(nElm, omp_get_max_threads(), OMP_MAX_SIZE);

  if (vec_val == nullptr)
    vec_val = MemoryAllocation::aligned_alloc<ScalarType>(64, nElm*sizeof(ScalarType));

  if(val != nullptr) {
    if(!valIsArray) {
      for(auto i=0ul; i<nElm; i++) vec_val[i] = *val;
    }
    else {
      for(auto i=0ul; i<nElm; i++) vec_val[i] = val[i];
    }
  }
}

template<class ScalarType>
template<class T>
void CSysVector<ScalarType>::PassiveCopy(const CSysVector<T>& other) {

  /*--- This is a method and not the overload of an operator to make sure who
   calls it knows the consequence to the derivative information (lost) ---*/

  /*--- check if self-assignment, otherwise perform deep copy ---*/
  if ((const void*)this == (const void*)&other) return;

  SU2_OMP_MASTER
  Initialize(other.GetNBlk(), other.GetNBlkDomain(), other.GetNVar(), nullptr, true);
  SU2_OMP_BARRIER

  PARALLEL_FOR
  for(auto i=0ul; i<nElm; i++)
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

  PARALLEL_FOR
  for(auto i=0ul; i<nElm; i++) vec_val[i] = a * x.vec_val[i];
}

template<class ScalarType>
void CSysVector<ScalarType>::Plus_AX(ScalarType a, const CSysVector<ScalarType> & x) {

  assert(nElm == x.nElm && "Sizes do not match");

  PARALLEL_FOR
  for(auto i=0ul; i<nElm; i++) vec_val[i] += a * x.vec_val[i];
}

template<class ScalarType>
void CSysVector<ScalarType>::Equals_AX_Plus_BY(ScalarType a, const CSysVector<ScalarType> & x,
                                               ScalarType b, const CSysVector<ScalarType> & y) {
  assert(nElm == x.nElm && nElm == y.nElm && "Sizes do not match");

  PARALLEL_FOR
  for(auto i=0ul; i<nElm; i++)
    vec_val[i] = a * x.vec_val[i] + b * y.vec_val[i];
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator=(const CSysVector<ScalarType> & u) {

  assert(nElm == u.nElm && "Sizes do not match");

  PARALLEL_FOR
  for(auto i=0ul; i<nElm; i++) vec_val[i] = u.vec_val[i];

  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator=(ScalarType val) {

  PARALLEL_FOR
  for(auto i=0ul; i<nElm; i++) vec_val[i] = val;

  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator+=(const CSysVector<ScalarType> & u) {

  assert(nElm == u.nElm && "Sizes do not match");

  PARALLEL_FOR
  for(auto i=0ul; i<nElm; i++) vec_val[i] += u.vec_val[i];

  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator-=(const CSysVector<ScalarType> & u) {

  assert(nElm == u.nElm && "Sizes do not match");

  PARALLEL_FOR
  for(auto i=0ul; i<nElm; i++) vec_val[i] -= u.vec_val[i];

  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator*=(ScalarType val) {

  PARALLEL_FOR
  for(auto i=0ul; i<nElm; i++) vec_val[i] *= val;

  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator/=(ScalarType val) {

  PARALLEL_FOR
  for(auto i=0ul; i<nElm; i++) vec_val[i] /= val;

  return *this;
}

template<class ScalarType>
void CSysVector<ScalarType>::CopyToArray(ScalarType* u_array) const {

  PARALLEL_FOR
  for(auto i=0ul; i<nElm; i++) u_array[i] = vec_val[i];
}

template<class ScalarType>
ScalarType CSysVector<ScalarType>::dot(const CSysVector<ScalarType> & u) const {
#if !defined(CODI_FORWARD_TYPE) && !defined(CODI_REVERSE_TYPE)

  /*--- All threads get the same "view" of the vectors and shared variable. ---*/
  SU2_OMP_BARRIER
  dotRes = 0.0;
  SU2_OMP_BARRIER

  /*--- Reduction over all threads in this mpi rank using the shared variable. ---*/
  ScalarType sum = 0.0;

  PARALLEL_FOR
  for(auto i=0ul; i<nElmDomain; ++i)
    sum += vec_val[i]*u.vec_val[i];

  SU2_OMP(atomic)
  dotRes += sum;

  /*--- Wait for all atomic updates. ---*/
  SU2_OMP_BARRIER

#ifdef HAVE_MPI
  /*--- Reduce across all mpi ranks, only master thread communicates. ---*/
  SU2_OMP_MASTER
  {
    sum = dotRes;
    SelectMPIWrapper<ScalarType>::W::Allreduce(&sum, &dotRes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  /*--- Make view of result consistent across threads. ---*/
  SU2_OMP_BARRIER
#endif // MPI
#else // CODI_TYPE
  /*--- Compatible version, no OMP reductions, no atomics, master does everything. ---*/
  SU2_OMP_BARRIER
  SU2_OMP_MASTER
  {
    ScalarType sum = 0.0;
    for(auto i=0ul; i<nElmDomain; ++i)
      sum += vec_val[i]*u.vec_val[i];
#ifdef HAVE_MPI
    /*--- Reduce across all mpi ranks. ---*/
    SelectMPIWrapper<ScalarType>::W::Allreduce(&sum, &dotRes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    dotRes = sum;
#endif // MPI
  }
  SU2_OMP_BARRIER
#endif // CODI
  return dotRes;
}

/*--- Explicit instantiations ---*/
template class CSysVector<su2double>;
template void CSysVector<su2double>::PassiveCopy(const CSysVector<su2double>&);

#ifdef CODI_REVERSE_TYPE
template class CSysVector<passivedouble>;
template void CSysVector<su2double>::PassiveCopy(const CSysVector<passivedouble>&);
template void CSysVector<passivedouble>::PassiveCopy(const CSysVector<su2double>&);
#endif
