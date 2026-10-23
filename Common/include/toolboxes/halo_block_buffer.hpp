/*!
 * \file halo_block_buffer.hpp
 * \brief Message buffers for matrix blocks, with integer indices kept apart from scalar values.
 * \author SU2 Contributors
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

#pragma once

#include <vector>

/*!
 * \class CHaloBlockBuffer
 * \brief Per point: a block count, then a global column index and a block per entry, in separate arrays.
 */
template <class ScalarType>
struct CHaloBlockBuffer {
  std::vector<unsigned long> index;
  std::vector<ScalarType> value;

  /*!
   * \brief Starts a new point.
   * \return Position of the point's block count.
   */
  size_t BeginPoint() {
    index.push_back(0);
    return index.size() - 1;
  }

  /*!
   * \brief Appends a block of n scalars to the point started at countPos.
   */
  void AddBlock(size_t countPos, unsigned long column, const ScalarType* block, size_t n) {
    index.push_back(column);
    value.insert(value.end(), block, block + n);
    index[countPos] += 1;
  }
};

/*!
 * \class CHaloBlockReader
 * \brief Walks the arrays written by CHaloBlockBuffer.
 */
template <class ScalarType>
class CHaloBlockReader {
  const unsigned long* idx;
  const ScalarType* val;
  size_t blkSize;

 public:
  CHaloBlockReader(const unsigned long* index, const ScalarType* value, size_t blockSize)
      : idx(index), val(value), blkSize(blockSize) {}

  /*!
   * \return Number of blocks of the next point.
   */
  unsigned long NextPoint() { return *idx++; }

  /*!
   * \brief Reads the next block of the current point.
   * \param[out] column - Its global column index.
   * \return Pointer to its scalars.
   */
  const ScalarType* NextBlock(unsigned long& column) {
    column = *idx++;
    const auto* blk = val;
    val += blkSize;
    return blk;
  }
};
