/*!
 * \file CSysVector.cpp
 * \brief Implementation and explicit instantiations of CSysVector.
 * \author P. Gomes, F. Palacios, J. Hicken, T. Economon
 * \version 7.0.6 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../include/toolboxes/allocation_toolbox.hpp"

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
   * calls it knows the consequence to the derivative information (lost) ---*/

  /*--- check if self-assignment, otherwise perform deep copy ---*/
  if ((const void*)this == (const void*)&other) return;

  SU2_OMP_MASTER
  Initialize(other.GetNBlk(), other.GetNBlkDomain(), other.GetNVar(), nullptr, true);
  SU2_OMP_BARRIER

  CSYSVEC_PARFOR
  for(auto i=0ul; i<nElm; i++)
    vec_val[i] = SU2_TYPE::GetValue(other[i]);
}

template<class ScalarType>
CSysVector<ScalarType>::~CSysVector() {

  MemoryAllocation::aligned_free(vec_val);
}

/*--- Explicit instantiations ---*/
/*--- We allways need su2double (regardless if it is passive or active). ---*/
template class CSysVector<su2double>;
#if defined(CODI_REVERSE_TYPE) || defined(USE_MIXED_PRECISION)
/*--- In reverse AD (or with mixed precision) we will also have passive (or float) vectors,
 *    and copy operations between them and active (or double) vectors, respectively. ---*/
template class CSysVector<su2mixedfloat>;
template void CSysVector<su2mixedfloat>::PassiveCopy(const CSysVector<su2double>&);
template void CSysVector<su2double>::PassiveCopy(const CSysVector<su2mixedfloat>&);
#endif
