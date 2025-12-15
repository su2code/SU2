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
inline unsigned long HashCombine(unsigned long v1, unsigned long v2) {
  const unsigned long prime = 1099511628211ULL;
  v1 ^= v2;
  v1 *= prime;
  return v1;
}

/*!
 * \brief Convert a double to a 64-bit integer suitable for hashing.
 * \param[in] x Double to integer.
 * \return Hash value of the double (not portable).
 */
inline unsigned long ToUInt64(double x) { return std::hash<double>{}(x); }

/*!
 * \brief Build a deterministic seed from physical time.
 * \param[in] x First integer value.
 * \param[in] y Second integer value.
 * \return 64-bit seed value.
 */
inline unsigned long GetSeed(unsigned long x, unsigned long y) { return HashCombine(x, y); }

/*!
 * \brief Generate a standard normally-distributed random number.
 * \param[in] gen Pseudo-random number generator.
 * \param[in] mean Mean of the normal distribution (default 0).
 * \param[in] stddev Standard deviation of the normal distribution (default 1).
 * \return Normally-distributed random number.
 */
inline double GetRandomNormal(std::mt19937 gen, double mean = 0.0, double stddev = 1.0) {
  std::normal_distribution<double> rnd(mean, stddev);
  return rnd(gen);
}

/*!
 * \brief Generate a uniformly-distributed random number.
 * \param[in] gen Pseudo-random number generator.
 * \param[in] xmin Lower boundary of the interval (default 0).
 * \param[in] xmax Upper bounary of the interval (default 1).
 * \return Uniformly-distributed random number.
 */
inline double GetRandomUniform(std::mt19937 gen, double xmin = 0.0, double xmax = 1.0) {
  std::uniform_real_distribution<double> rnd(xmin, xmax);
  return rnd(gen);
}

/*!
 * \brief Compute modified bessel function of first kind (order 0).
 * \param[in] x Argument of Bessel funtion.
 * \return Value of Bessel function.
 */
inline double GetBesselZero(double x) {
  double abx = fabs(x);

  if (abx < 3.75) {
    double t = x / 3.75;
    double p =
        1.0 +
        t * t *
            (3.5156229 +
             t * t * (3.0899424 + t * t * (1.2067492 + t * t * (0.2659732 + t * t * (0.0360768 + t * t * 0.0045813)))));
    return log(p);
  } else {
    double t = 3.75 / abx;
    double poly =
        0.39894228 +
        t * (0.01328592 +
             t * (0.00225319 +
                  t * (-0.00157565 +
                       t * (0.00916281 + t * (-0.02057706 + t * (0.02635537 + t * (-0.01647633 + t * 0.00392377)))))));
    return abx - 0.5 * log(abx) + log(poly);
  }
}

/*!
 * \brief Compute integral involving product of three modified Bessel functions.
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
  const double t_max = 40.0;
  const double dt = t_max / N;

  double sum = 0.0;

  for (int i = 1; i < N; i++) {
    double t = i * dt;

    double e = exp(-A * t);

    double lx = GetBesselZero(Bx * t);
    double ly = GetBesselZero(By * t);
    double lz = GetBesselZero(Bz * t);

    double lin = log(t) - A * t + lx + ly + lz;

    double integrand = exp(lin);
    sum += integrand;
  }

  return sum * dt;
}

/// @}
}  // namespace RandomToolbox
