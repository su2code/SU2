/*!
 * \file CFastFindAndEraseQueue.hpp
 * \brief A queue-type container (push back, pop front), but with
 *        fast deletion of arbitrary items (possibly in the middle).
 * \author P. Gomes
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
#pragma once

#include <cassert>
#include <stdlib.h>
#include <limits>
#include <vector>
#include <unordered_map>

/*!
 * \class CFastFindAndEraseQueue
 * \ingroup Containers
 * \brief A queue-type container (push back, pop front), but with
 *        fast deletion of arbitrary items (possibly in the middle).
 * \param[in] ItemType_ - Type of the stored items.
 * \param[in] ErasedValue_ - Value used to mark items for erasure.
 * \param[in] CleanupThreshold - Number of marked items that triggers full cleanup.
 * \note It would not be a good idea to use non-trivial item types.
 */
template <typename ItemType_ = unsigned long, ItemType_ ErasedValue_ = std::numeric_limits<ItemType_>::max(),
          size_t CleanupThreshold = 8192>
class CFastFindAndEraseQueue {
 public:
  enum { ErasedValue = ErasedValue_ };
  using ItemType = ItemType_;
  using Iterator = typename std::vector<ItemType>::const_iterator;

 private:
  size_t erasedCounter = 0;                     /*!< \brief How many items have been marked since last cleanup. */
  std::vector<ItemType> items;                  /*!< \brief The stored items. */
  std::unordered_map<ItemType, size_t> indexes; /*!< \brief Map items to their location. */

  /*!
   * \brief Cleanup, shifts non-erased items forward, re-mapping them.
   * \note This is the performance hotspot (due to the re-mapping).
   */
  void cleanup() {
    size_t idx = 0;
    for (auto value : items) {
      if (value != ErasedValue) {
        items[idx] = value;
        indexes[value] = idx;
        ++idx;
      }
    }
    items.resize(idx);
    erasedCounter = 0;
  }

 public:
  /*!
   * \brief Default construct.
   */
  CFastFindAndEraseQueue() = default;

  /*!
   * \brief Construct and initialize with [0,N) range.
   */
  CFastFindAndEraseQueue(size_t N) {
    items.resize(N);
    for (size_t i = 0; i < N; ++i) {
      items[i] = i;
      indexes[i] = i;
    }
  }

  /*--- Wrapper methods for stl-type operations. ---*/
  size_t size() const { return indexes.size(); }
  bool empty() const { return indexes.empty(); }
  Iterator begin() const { return items.cbegin(); }
  Iterator end() const { return items.cend(); }

  /*!
   * \brief Push back an item, that should be unique.
   * \param[in] value - The new item.
   */
  void push_back(ItemType value) {
    assert(value != ErasedValue && "#bad_idea");
    assert(!indexes.count(value) && "Items must be unique.");
    indexes[value] = items.size();
    items.push_back(value);
  }

  /*!
   * \brief Get the item in front of the queue.
   */
  ItemType front() const {
    for (auto value : items)
      if (value != ErasedValue) return value;
    assert(false && "Queue is empty.");
    return ErasedValue;
  }

  /*!
   * \brief Main event, look for an item and erase it if present.
   * \param[in] value - Item to look for.
   * \return True, if something was erased, false otherwise.
   */
  bool findAndErase(const ItemType value) {
    if (empty()) return false;

    /*--- Get the item position. ---*/
    auto it = indexes.find(value);
    if (it == indexes.end()) return false;  // value is not here

    const auto idx = it->second;
    assert(items[idx] != ErasedValue && "The item was already erased?!");

    /*--- Erase and cleanup items. ---*/
    indexes.erase(it);
    items[idx] = ErasedValue;
    if (++erasedCounter == CleanupThreshold) cleanup();

    return true;
  }
};
