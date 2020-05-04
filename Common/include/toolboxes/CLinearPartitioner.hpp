/*!
 * \file CLinearPartitioner.hpp
 * \brief Header file for the class CLinearPartitioner.
 *        The implementations are in the <i>CLinearPartitioner.cpp</i> file.
 * \author T. Economon
 * \version 7.0.4 "Blackbird"
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

#pragma once

#include "../mpi_structure.hpp"

#include <vector>

using namespace std;

/*!
 * \class CLinearPartitioner
 * \brief Helper class that provides the counts for each rank in a linear
 *        partitioning given the global count as input.
 * \author: T. Economon
 */
class CLinearPartitioner {
  
protected:
  
  int size;  /*!< \brief MPI Size. */

  vector<unsigned long> firstIndex;               /*!< \brief Vector containing the first index on each rank due to a linear partitioning by global count. */
  vector<unsigned long> lastIndex;                /*!< \brief Vector containing the last index on each rank due to a linear partitioning by global count. */
  vector<unsigned long> sizeOnRank;               /*!< \brief Vector containing the total size of the current rank's linear partition. */
  vector<unsigned long> cumulativeSizeBeforeRank; /*!< \brief Vector containing the cumulative size of all linear partitions before the current rank. */
  
public:

  /*!
   * \brief Constructor of the CLinearPartitioner class.
   * \param[in] val_global_count - global count to be linearly partitioned.
   * \param[in] val_offset - offset from 0 for the first index on rank 0 (typically 0).
   * \param[in] isDisjoint - boolean controlling whether the linear partitions should be disjoint (default is false).
   */
  CLinearPartitioner(unsigned long val_global_count,
                     unsigned long val_offset,
                     bool          isDisjoint = false);

  /*!
   * \brief Destructor of the CLinearPartitioner class.
   */
  ~CLinearPartitioner(void);
  
  /*!
   * \brief Get the rank that owns the index based on the linear partitioning.
   * \param[in] val_index - Current index.
   * \returns Owning rank for the current index based on linear partitioning.
   */
  unsigned long GetRankContainingIndex(unsigned long val_index);
  
  /*!
   * \brief Get the first index of the current rank's linear partition.
   * \param[in] val_rank - MPI rank identifier.
   * \returns First index of the current rank's linear partition.
   */
  inline unsigned long GetFirstIndexOnRank(int val_rank) {
    return firstIndex[val_rank];
  }
  
  /*!
   * \brief Get the last index of the current rank's linear partition.
   * \param[in] val_rank - MPI rank identifier.
   * \returns Last index of the current rank's linear partition.
   */
  inline unsigned long GetLastIndexOnRank(int val_rank) {
    return lastIndex[val_rank];
  }
  
  /*!
   * \brief Get the total size of the current rank's linear partition.
   * \param[in] val_rank - MPI rank identifier.
   * \returns Size of the current rank's linear partition.
   */
  inline unsigned long GetSizeOnRank(int val_rank) {
    return sizeOnRank[val_rank];
  }
  
  /*!
   * \brief Get the cumulative size of all linear partitions before the current rank.
   * \param[in] val_rank - MPI rank identifier.
   * \returns Cumulative size of all linear partitions before the current rank.
   */
  inline unsigned long GetCumulativeSizeBeforeRank(int val_rank) {
    return cumulativeSizeBeforeRank[val_rank];
  }
  
};
