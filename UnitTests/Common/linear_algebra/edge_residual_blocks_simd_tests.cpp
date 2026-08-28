/*!
 * \file edge_residual_blocks_simd_tests.cpp
 * \brief Unit tests for the SIMD overloads of CSysMatrix::SetBlocks (four independent blocks),
 *        SetOffDiagBlocks, and CSysVector::AddBlock, added for the vectorized scalar edge loop.
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

#include <set>

#include "catch.hpp"
#include "../../UnitQuadTestCase.hpp"
#include "../../../SU2_CFD/include/numerics/util.hpp"

/*--- UnitQuadTestCase's fixed config: 3D compressible Euler, nVar = nDim + 2. ---*/
constexpr size_t nVar = 5;
using Double = simd::Array<su2double>;
constexpr size_t N = Double::Size;

static void FillBlock(Matrix<Double, nVar, nVar>& block, size_t lane, su2double base) {
  for (auto i = 0u; i < nVar; ++i)
    for (auto j = 0u; j < nVar; ++j) block(i, j)[lane] = base + 0.1 * i + 0.01 * j + lane;
}

static void CheckBlock(const CSysMatrix<su2mixedfloat>& matrix, unsigned long i, unsigned long j,
                       const Matrix<Double, nVar, nVar>& block, size_t lane, double tol) {
  auto view = matrix.GetBlockView(i, j);
  for (auto iVar = 0u; iVar < nVar; ++iVar)
    for (auto jVar = 0u; jVar < nVar; ++jVar)
      CHECK(SU2_TYPE::GetValue(view(iVar, jVar)) == Approx(SU2_TYPE::GetValue(block(iVar, jVar)[lane])).margin(tol));
}

TEST_CASE("SIMD SetBlocks, SetOffDiagBlocks and AddBlock match their scalar counterparts", "[LinearAlgebra]") {
  cout.rdbuf(nullptr);

  UnitQuadTestCase testCase;
  testCase.InitConfig();
  testCase.InitGeometry();
  testCase.InitSolver();

  cout.rdbuf(testCase.orig_buf);

  auto* solver = testCase.solver[FLOW_SOL];
  auto& matrix = solver->Jacobian;
  REQUIRE(solver->GetnVar() == nVar);

  /*--- N edges with pairwise disjoint endpoints, so a lane's diagonal block never accumulates
   * a second edge's contribution and CheckBlock can compare it to a single lane's own fill.
   * (The production code relies on edges within a SIMD group being contiguous for the SIMD
   * GetNode gather, but this test builds iPoint/jPoint with the scalar GetNode overload, one
   * edge at a time, so that constraint does not apply here.) ---*/
  simd::Array<unsigned long, N> iEdge, iPoint, jPoint;
  std::set<unsigned long> used;
  size_t found = 0;
  for (unsigned long e = 0; e < testCase.geometry->GetnEdge() && found < N; ++e) {
    const auto p0 = testCase.geometry->edges->GetNode(e, 0);
    const auto p1 = testCase.geometry->edges->GetNode(e, 1);
    if (used.count(p0) || used.count(p1)) continue;
    iEdge[found] = e;
    iPoint[found] = p0;
    jPoint[found] = p1;
    used.insert(p0);
    used.insert(p1);
    ++found;
  }
  REQUIRE(found == N);

  Matrix<Double, nVar, nVar> jac_ii, jac_ij, jac_ji, jac_jj;
  for (size_t k = 0; k < N; ++k) {
    FillBlock(jac_ii, k, 1.0);
    FillBlock(jac_ij, k, 2.0);
    FillBlock(jac_ji, k, 3.0);
    FillBlock(jac_jj, k, 4.0);
  }

  SECTION("SetBlocks: diagonal accumulates, off-diagonal is set, one lane at a time via the scalar overload") {
    matrix.SetValZero();
    matrix.SetBlocks(iEdge, iPoint, jPoint, jac_ii, jac_ij, jac_ji, jac_jj);

    for (size_t k = 0; k < N; ++k) {
      CheckBlock(matrix, iPoint[k], iPoint[k], jac_ii, k, 1e-6);
      CheckBlock(matrix, jPoint[k], jPoint[k], jac_jj, k, 1e-6);
      CheckBlock(matrix, iPoint[k], jPoint[k], jac_ij, k, 1e-6);
      CheckBlock(matrix, jPoint[k], iPoint[k], jac_ji, k, 1e-6);
    }
  }

  SECTION("SetBlocks: a lane whose mask is 0 is left untouched") {
    matrix.SetValZero();
    matrix.SetBlocks(iEdge, iPoint, jPoint, jac_ii, jac_ij, jac_ji, jac_jj);

    Matrix<Double, nVar, nVar> jac_ii2 = jac_ii, jac_ij2 = jac_ij, jac_ji2 = jac_ji, jac_jj2 = jac_jj;
    for (size_t k = 0; k < N; ++k) FillBlock(jac_ii2, k, 100.0);

    simd::Array<su2double, N> mask = 1;
    mask[0] = 0;
    matrix.SetBlocks(iEdge, iPoint, jPoint, jac_ii2, jac_ij2, jac_ji2, jac_jj2, mask);

    /*--- Lane 0: SetBlocks was skipped, diagonal is the single accumulation from the first
     * call above, off-diagonal is still what that first call set. ---*/
    CheckBlock(matrix, iPoint[0], iPoint[0], jac_ii, 0, 1e-6);
    CheckBlock(matrix, iPoint[0], jPoint[0], jac_ij, 0, 1e-6);

    /*--- Lane 1 (points disjoint from every other lane's, by construction): diagonal
     * accumulated twice, off-diagonal overwritten by the second call. ---*/
    if (N > 1) {
      Matrix<Double, nVar, nVar> jac_ii_2x = jac_ii;
      for (auto i = 0u; i < nVar; ++i)
        for (auto j = 0u; j < nVar; ++j) jac_ii_2x(i, j)[1] = jac_ii(i, j)[1] + jac_ii2(i, j)[1];
      CheckBlock(matrix, iPoint[1], iPoint[1], jac_ii_2x, 1, 1e-6);
      CheckBlock(matrix, iPoint[1], jPoint[1], jac_ij2, 1, 1e-6);
    }
  }

  SECTION("SetOffDiagBlocks leaves the diagonal untouched") {
    matrix.SetValZero();
    matrix.SetBlocks(iEdge, iPoint, jPoint, jac_ii, jac_ij, jac_ji, jac_jj);

    Matrix<Double, nVar, nVar> jac_ij_new, jac_ji_new;
    for (size_t k = 0; k < N; ++k) {
      FillBlock(jac_ij_new, k, 5.0);
      FillBlock(jac_ji_new, k, 6.0);
    }
    matrix.SetOffDiagBlocks(iEdge, jac_ij_new, jac_ji_new);

    for (size_t k = 0; k < N; ++k) {
      CheckBlock(matrix, iPoint[k], iPoint[k], jac_ii, k, 1e-6);
      CheckBlock(matrix, jPoint[k], jPoint[k], jac_jj, k, 1e-6);
      CheckBlock(matrix, iPoint[k], jPoint[k], jac_ij_new, k, 1e-6);
      CheckBlock(matrix, jPoint[k], iPoint[k], jac_ji_new, k, 1e-6);
    }
  }
}

TEST_CASE("SIMD CSysVector::AddBlock writes independent values to independent points", "[LinearAlgebra]") {
  cout.rdbuf(nullptr);

  UnitQuadTestCase testCase;
  testCase.InitConfig();
  testCase.InitGeometry();
  testCase.InitSolver();

  cout.rdbuf(testCase.orig_buf);

  auto* solver = testCase.solver[FLOW_SOL];
  auto& vector = solver->LinSysRes;
  REQUIRE(solver->GetnVar() == nVar);
  REQUIRE(testCase.geometry->GetnPoint() >= static_cast<long>(N));

  simd::Array<unsigned long, N> iPoint;
  Vector<Double, nVar> block;
  for (size_t k = 0; k < N; ++k) {
    iPoint[k] = k;
    for (auto i = 0u; i < nVar; ++i) block(i)[k] = 1.0 + 0.1 * i + k;
  }

  vector.SetValZero();
  vector.AddBlock(iPoint, block);
  vector.AddBlock(iPoint, block);

  for (size_t k = 0; k < N; ++k) {
    const auto* row = vector.GetBlock(iPoint[k]);
    for (auto i = 0u; i < nVar; ++i)
      CHECK(SU2_TYPE::GetValue(row[i]) == Approx(2.0 * SU2_TYPE::GetValue(block(i)[k])).margin(1e-6));
  }
}
