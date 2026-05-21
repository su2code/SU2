/*!
 * \file CSysVector.cpp
 * \brief Implementation and explicit instantiations of CSysVector.
 * \author P. Gomes, F. Palacios, J. Hicken, T. Economon
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

  d_vec_val = GPUMemoryAllocation::gpu_alloc<ScalarType, true>(nElm * sizeof(ScalarType));

#ifdef HAVE_OMP
  dot_scratch.reset(new ScalarType[omp_get_max_threads()]);
#endif

  if (val != nullptr) {
    if (!valIsArray) {
      for (auto i = 0ul; i < nElm; i++) vec_val[i] = *val;
    } else {
      for (auto i = 0ul; i < nElm; i++) vec_val[i] = val[i];
    }
  }
}

template <class ScalarType>
const su2matrix<ScalarType>& CSysVector<ScalarType>::multiDot(const std::vector<CSysVector<ScalarType>>& V,
                                                              const size_t i0, const size_t n,
                                                              const std::vector<CSysVector<ScalarType>>& W,
                                                              const size_t m) {
  SU2_ZONE_SCOPED
  static constexpr size_t BLOCK_SIZE = 1024;
  static su2matrix<ScalarType> shared;

  if (n == 0 || m == 0) return shared;

  SU2_OMP_BARRIER
  const size_t size = V[0].nElmDomain;

  su2matrix<ScalarType> local(n, m);
  local.setConstant(0);

  SU2_OMP_FOR_(schedule(static) SU2_NOWAIT)
  for (size_t offset = 0; offset < size; offset += BLOCK_SIZE) {
    const auto limit = std::min(offset + BLOCK_SIZE, size);
    for (size_t i = 0; i < n; ++i) {
      const auto& vi = V[i0 + i];
      for (size_t j = 0; j < m; ++j) {
        const auto& wj = W[j];
        ScalarType sum = 0.0;
        SU2_OMP_SIMD
        for (auto k = offset; k < limit; ++k) {
          sum += vi[k] * wj[k];
        }
        local(i, j) += sum;
      }
    }
  }
  END_SU2_OMP_FOR

  /*--- Reduce over all threads in an ordered way to ensure a deterministic result. ---*/
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      W[j].dot_scratch[omp_get_thread_num()] = local(i, j);
    }
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    for (size_t j = 0; j < m; ++j) {
      for (int t = 1; t < omp_get_num_threads(); ++t) {
        local(i, j) += W[j].dot_scratch[t];
      }
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Single AllReduce of the result, only the master thread communicates. ---*/
  SU2_OMP_MASTER {
    shared.resize(n, m);

    const auto mpi_type = (sizeof(ScalarType) < sizeof(double)) ? MPI_FLOAT : MPI_DOUBLE;
    SelectMPIWrapper<ScalarType>::W::Allreduce(local.data(), shared.data(), n * m, mpi_type, MPI_SUM,
                                               SU2_MPI::GetComm());
  }
  END_SU2_OMP_MASTER

  /*--- All threads have the same view of the result. ---*/
  SU2_OMP_BARRIER

  return shared;
}

template <class ScalarType>
CSysVector<ScalarType>::~CSysVector() {
  if constexpr (!std::is_trivial_v<ScalarType>) {
    for (auto i = 0ul; i < nElm; i++) vec_val[i].~ScalarType();
  }
  MemoryAllocation::aligned_free(vec_val);

  GPUMemoryAllocation::gpu_free(d_vec_val);
}

/*--- Explicit instantiations ---*/
template class CSysVector<su2mixedfloat>;
#ifdef USE_MIXED_PRECISION
template class CSysVector<passivedouble>;
#endif
#ifdef CODI_REVERSE_TYPE
template class CSysVector<su2double>;
#endif
