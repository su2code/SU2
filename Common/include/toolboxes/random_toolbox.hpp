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
#include <cstdint>

namespace RandomToolbox {
/// \addtogroup RandomToolbox
/// @{

/*!
 * \brief SplitMix64 hash function for 64-bit integers.
 * \param[in] x Input value to hash.
 * \return Hashed 64-bit output.
 */
static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;                  // golden ratio offset
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL; // first mixing step
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; // second mixing step
    return x ^ (x >> 31);                        // final avalanche
}

/*!
 * \brief Generate a deterministic 64-bit hash from three integers.
 * \param[in] nodeIndex Global node index.
 * \param[in] iDim Dimension index.
 * \param[in] timeIter Current time iteration of the simulation.
 * \return 64-bit hash value.
 */
inline uint64_t GetHash(unsigned long nodeIndex, unsigned short iDim, unsigned long timeIter) {
    uint64_t x = nodeIndex;
    x ^= splitmix64(iDim);
    x ^= splitmix64(timeIter);
    return splitmix64(x);
}

/*!
 * \brief Convert a 64-bit hash into a uniform double in (0,1].
 * Uses the top 53 bits of the hash to fill the mantissa of a double.
 * Ensures the result is never zero, suitable for Box-Muller transform.
 * \param[in] x 64-bit hash.
 * \return Uniform double in the interval (0,1].
 */
inline double HashToUniform(uint64_t x) {
    constexpr double inv53 = 1.0 / 9007199254740992.0; // 1/2^53
    uint64_t uInt = x >> 11;                           // top 53 bits
    return (uInt + 1) * inv53;                         // map to (0,1]
}

/*!
 * \brief Generate a standard normal random number from a 64-bit hash.
 * Uses two deterministic uniforms derived from the hash and its bitwise NOT
 * as inputs to the Box-Muller transform.
 * \param[in] x 64-bit hash.
 * \return Standard normal random number (mean=0, stddev=1).
 */
inline double HashToNormal(uint64_t x) {
    constexpr double pi = 3.14159265358979323846;
    double u = HashToUniform(x);    // first uniform
    double v = HashToUniform(~x);   // second uniform (bitwise NOT)
    double r = sqrt(-2.0 * log(u));  
    double theta = 2.0 * pi * v;
    return r * cos(theta);          // one normal sample
}

/*!
 * \brief Generate a deterministic standard normal number for a cell, dimension, and timestep.
 * 
 * Combines hashing and Box-Muller in one function.
 * 
 * \param[in] nodeIndex Global node index.
 * \param[in] dim Dimension index.
 * \param[in] timeIter Simulation timestep (1-based).
 * \return Standard normal random number.
 */
inline double GetNormal(unsigned long nodeIndex, unsigned long dim, unsigned long timeIter) {
    uint64_t hash = GetHash(nodeIndex, dim, timeIter);
    return HashToNormal(hash);
}

/*!
 * \brief Compute modified bessel function of first kind (order 0).
 * \param[in] x Argument of Bessel funtion.
 * \return Value of Bessel function.
 */
inline double GetBesselZero(double x) {
  double abx = fabs(x);
  if (abx < 3.75) {
    double t = abx / 3.75;
    double p = 1.0 + t * t * (3.5156229 + t * t * (3.0899424 + t * t * (1.2067492 + t * t * (0.2659732 + t * t * (0.0360768 + t * t * 0.0045813)))));
    return log(p);
  } else {
    double t = 3.75 / abx;
    double poly = 0.39894228 + t * (0.01328592 + t * (0.00225319 + t * (-0.00157565 + t * (0.00916281 + t * (-0.02057706 + t * (0.02635537 + t * (-0.01647633 + t * 0.00392377)))))));
    return abx - 0.5 * log(abx) + log(poly);
  }
}

/*!
 * \brief Compute integral involving the product of three modified Bessel functions. 
 *  Useful for scaling the smoothed stochastic source terms in Langevin equations.
 * \param[in] beta_x Argument in x-direction.
 * \param[in] beta_y Argument in y-direction.
 * \param[in] beta_z Argument in z-direction.
 * \return Value of the integral.
 */
inline double GetBesselIntegral(double beta_x, double beta_y, double beta_z) {
  const double A = 1.0 + 2.0 * (beta_x + beta_y + beta_z);
  const double Bx = 2.0 * beta_x;
  const double By = 2.0 * beta_y;
  const double Bz = 2.0 * beta_z;
  const int N = 4000;
  const double t_max = 20.0;
  const double dt = t_max / N;
  double sum = 0.0;
  for (int i = 1; i <= N; i++) {
    double t = i * dt;
    double lx = GetBesselZero(Bx * t);
    double ly = GetBesselZero(By * t);
    double lz = GetBesselZero(Bz * t);
    double lin = log(t) - A * t + lx + ly + lz;
    double integrand = exp(lin);
    if (i==N) integrand *= 0.5;
    sum += integrand;
  }
  return sum * dt;
}

/// @}
}  // namespace RandomToolbox