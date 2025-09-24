/*!
 * \file random_toolbox.hpp
 * \brief Collection of utility functions for random number generation.
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include <random>
#include <functional>

namespace RandomToolbox {
/// \addtogroup RandomToolbox
/// @{

/*!
 * \brief Combine two 64-bit integers into a single hash value.
 * \param[in] v1 Current hash value.
 * \param[in] v2 Value to mix in.
 * \return Combined hash value.
 */
inline uint64_t HashCombine(uint64_t v1, uint64_t v2) {
  const uint64_t prime = 1099511628211ULL;
  v1 ^= v2;
  v1 *= prime;
  return v1;
}

/*!
 * \brief Convert a double to a 64-bit integer suitable for hashing.
 * \param[in] x Double to integer.
 * \return Hash value of the double (not portable).
 */
inline uint64_t ToUInt64(double x) { return std::hash<double>{}(x); }

/*!
 * \brief Build a deterministic seed from physical time.
 * \param[in] x First integer value.
 * \param[in] y Second integer value.
 * \return 64-bit seed value.
 */
inline uint64_t GetSeed(uint64_t x, uint64_t y) { return HashCombine(x, y); }

/*!
 * \brief Generate a standard normally-distributed random number.
 * \param[in] seed Seed for the random number generator.
 * \param[in] mean Mean of the normal distribution (default 0).
 * \param[in] stddev Standard deviation of the normal distribution (default 1).
 * \return Normally-distributed random number.
 */
inline double GetRandomNormal(uint64_t seed, double mean = 0.0, double stddev = 1.0) {
  std::mt19937_64 gen(seed);
  std::normal_distribution<double> rnd(mean, stddev);
  return rnd(gen);
}

/// @}
}  // namespace RandomToolbox
