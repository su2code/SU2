/*!
 * \file CSysVector.cpp
 * \brief Implementation and explicit instantiations of CSysVector.
 * \author P. Gomes, F. Palacios, J. Hicken, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

template <class ScalarType>
void CSysVector<ScalarType>::Initialize(unsigned long numBlk, unsigned long numBlkDomain, unsigned long numVar,
                                        const ScalarType* val, bool valIsArray, bool errorIfParallel) {
  if (errorIfParallel && omp_in_parallel()) {
    assert(false);
    SU2_MPI::Error("If this class were constructed in parallel its operations would be incorrect.", CURRENT_FUNCTION);
  }

  if (omp_get_thread_num())
    SU2_MPI::Error("Only the master thread is allowed to initialize the vector.", CURRENT_FUNCTION);

  if (nElm != numBlk * numVar) {
    MemoryAllocation::aligned_free(vec_val);
    vec_val = nullptr;
  }

  nElm = numBlk * numVar;
  nElmDomain = numBlkDomain * numVar;
  nVar = numVar;

  omp_chunk_size = computeStaticChunkSize(nElm, omp_get_max_threads(), OMP_MAX_SIZE);

  if (vec_val == nullptr) vec_val = MemoryAllocation::aligned_alloc<ScalarType, true>(64, nElm * sizeof(ScalarType));

  if (val != nullptr) {
    if (!valIsArray) {
      for (auto i = 0ul; i < nElm; i++) vec_val[i] = *val;
    } else {
      for (auto i = 0ul; i < nElm; i++) vec_val[i] = val[i];
    }
  }
}

template <class ScalarType>
CSysVector<ScalarType>::~CSysVector() {
  if (!std::is_trivial<ScalarType>::value)
    for (auto i = 0ul; i < nElm; i++) vec_val[i].~ScalarType();
  MemoryAllocation::aligned_free(vec_val);
}

/*--- Explicit instantiations ---*/
/*--- We allways need su2double (regardless if it is passive or active). ---*/
template class CSysVector<su2double>;
#ifdef USE_MIXED_PRECISION
/*--- In reverse AD (or with mixed precision) we will also have passive (or float) vectors. ---*/
template class CSysVector<su2mixedfloat>;
#endif
#ifdef CODI_REVERSE_TYPE
template class CSysVector<passivedouble>;
#endif
