/*!
 * \file partition_load.hpp
 * \brief Load balance measure for graph partitions.
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

#include <algorithm>
#include <numeric>
#include <vector>

/*!
 * \brief Ratio of the largest partition load to the mean load.
 * \param[in] load - Total vertex weight of each partition.
 */
inline double PartitionImbalance(const std::vector<unsigned long>& load) {
  const auto sum = std::accumulate(load.begin(), load.end(), 0ul);
  if (sum == 0) return 1.0;
  return static_cast<double>(*std::max_element(load.begin(), load.end())) * load.size() / sum;
}
