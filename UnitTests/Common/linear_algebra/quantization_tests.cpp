/*!
 * \file Quantization_tests.cpp
 * \brief Unit tests for the int8 row-scaled block quantization used by Q_LU_SGS.
 * \author P. Gomes
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

#include "catch.hpp"
#include "../../../Common/include/linear_algebra/CSysMatrix.hpp"

/*--- Row-major access into a flat quantized block, matching CBlockView's decode. ---*/
static su2double Decode(const int8_t* qs, const int8_t* qv, unsigned long nVar, unsigned long r, unsigned long c) {
  return static_cast<su2double>(qv[r * nVar + c]) * DecodeQuantScale(qs[r]);
}

TEST_CASE("Quantization round-trip is stable for well-behaved values", "[LinearAlgebra]") {
  /*--- EncodeQuantBlock always encodes a full nVar x nVar block (one scale per row),
   * so every row is set to the same values here and only row 0 is inspected. Values
   * are moderate: none close enough to the row's own max to saturate, and none land
   * exactly on the int8 min (-128), so the derived per-row scale does not shift
   * between successive encode/decode passes. ---*/
  constexpr unsigned long nVar = 4;
  const su2double row[nVar] = {3.0, -2.0, 1.5, 0.75};
  auto f = [&](unsigned long, unsigned long c) { return row[c]; };

  int8_t qs1[nVar], qv1[nVar * nVar];
  EncodeQuantBlock(f, qs1, qv1, nVar);

  auto decoded = [&](unsigned long, unsigned long c) { return Decode(qs1, qv1, nVar, 0, c); };

  int8_t qs2[nVar], qv2[nVar * nVar];
  EncodeQuantBlock(decoded, qs2, qv2, nVar);

  CHECK(qs2[0] == qs1[0]);
  for (unsigned long c = 0; c < nVar; ++c) CHECK(qv2[c] == qv1[c]);

  /*--- A second decode/encode cycle from the now-stable representation must
   * reproduce the exact same codes again. ---*/
  auto decoded2 = [&](unsigned long, unsigned long c) { return Decode(qs2, qv2, nVar, 0, c); };
  int8_t qs3[nVar], qv3[nVar * nVar];
  EncodeQuantBlock(decoded2, qs3, qv3, nVar);

  CHECK(qs3[0] == qs2[0]);
  for (unsigned long c = 0; c < nVar; ++c) CHECK(qv3[c] == qv2[c]);
}

TEST_CASE("Quantization saturates and truncates with a known, bounded error", "[LinearAlgebra]") {
  /*--- One row engineered so that, at the scale it forces (2^0 = 1 here, since the
   * largest magnitude in the row is in [64, 128)):
   *   col 0:  127.9 -> rounds to 128, clamped to  127 (int8 max), error ~0.9 (close to
   *           the theoretical worst case: clamping can only ever push a value that
   *           rounds to +128 down to +127, an error that approaches but never reaches
   *           1 full scale unit).
   *   col 1: -127.9 -> rounds to -128 (int8 min), no clamping needed since -128 is a
   *           representable int8 value; error ~0.1 (near-exact). Clamping can only
   *           ever trigger on the positive side because int8 is asymmetric
   *           ([-128, 127]) while the encoding is symmetric around the row's max abs.
   *   col 2:   0.2  -> rounds to 0: small values are truncated away entirely.
   *   col 3:  60.0  -> reconstructs exactly, since the row's scale (2^0 = 1) is an
   *           exact power of two and 60 fits in int8. ---*/
  constexpr unsigned long nVar = 4;
  const su2double row[nVar] = {127.9, -127.9, 0.2, 60.0};
  auto f = [&](unsigned long, unsigned long c) { return row[c]; };

  int8_t qs[nVar], qv[nVar * nVar];
  EncodeQuantBlock(f, qs, qv, nVar);

  CHECK(qs[0] == 0);
  CHECK(static_cast<int>(qv[0]) == 127);   // Saturated at the int8 maximum.
  CHECK(static_cast<int>(qv[1]) == -128);  // Hits the int8 minimum exactly.
  CHECK(static_cast<int>(qv[2]) == 0);     // Truncated to zero.
  CHECK(static_cast<int>(qv[3]) == 60);    // Exact, mid-range value.

  const su2double scale = DecodeQuantScale(qs[0]);
  CHECK(scale == Approx(1.0));

  const su2double error0 = std::abs(row[0] - Decode(qs, qv, nVar, 0, 0));
  const su2double error1 = std::abs(row[1] - Decode(qs, qv, nVar, 0, 1));
  const su2double error2 = std::abs(row[2] - Decode(qs, qv, nVar, 0, 2));
  const su2double error3 = std::abs(row[3] - Decode(qs, qv, nVar, 0, 3));

  /*--- The quantization error for any entry is strictly bounded by one scale unit:
   * 0.5 from rounding to the nearest representable level, plus (only for values that
   * saturate on the positive side) less than another 0.5 from clamping 128 down to 127. ---*/
  CHECK(error0 < scale);
  CHECK(error1 < scale);
  CHECK(error2 < scale);
  CHECK(error3 == Approx(0.0).margin(1e-12));

  CHECK(error0 == Approx(0.9).margin(1e-5));
  CHECK(error1 == Approx(0.1).margin(1e-5));
  CHECK(error2 == Approx(0.2).margin(1e-5));
}
