/*!
 * \file edge_residual_blocks_tests.cpp
 * \brief Unit tests for CSysMatrix::SetBlocks (four independent blocks) and SetOffDiagBlocks.
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
#include "../../UnitQuadTestCase.hpp"

/*--- A block whose entries are all different, so a mixed-up index reads back wrong. ---*/
static void FillBlock(su2double block[][8], unsigned long nVar, su2double base) {
  for (auto i = 0u; i < nVar; ++i)
    for (auto j = 0u; j < nVar; ++j) block[i][j] = base + 0.1 * i + 0.01 * j;
}

static void CheckBlock(const CSysMatrix<su2mixedfloat>& matrix, unsigned long i, unsigned long j, unsigned long nVar,
                       const su2double block[][8], double tol) {
  auto view = matrix.GetBlockView(i, j);
  for (auto iVar = 0u; iVar < nVar; ++iVar)
    for (auto jVar = 0u; jVar < nVar; ++jVar)
      CHECK(SU2_TYPE::GetValue(view(iVar, jVar)) == Approx(block[iVar][jVar]).margin(tol));
}

TEST_CASE("SetBlocks and SetOffDiagBlocks assemble four independent blocks", "[LinearAlgebra]") {
  cout.rdbuf(nullptr);

  UnitQuadTestCase testCase;
  testCase.InitConfig();
  testCase.InitGeometry();
  testCase.InitSolver();

  cout.rdbuf(testCase.orig_buf);

  auto* solver = testCase.solver[FLOW_SOL];
  auto& matrix = solver->Jacobian;
  const auto nVar = solver->GetnVar();
  REQUIRE(testCase.geometry->GetnEdge() > 0);

  const auto iEdge = 0ul;
  const auto iPoint = testCase.geometry->edges->GetNode(iEdge, 0);
  const auto jPoint = testCase.geometry->edges->GetNode(iEdge, 1);

  su2double jac_ii[8][8], jac_ij[8][8], jac_ji[8][8], jac_jj[8][8];
  FillBlock(jac_ii, nVar, 1.0);
  FillBlock(jac_ij, nVar, 2.0);
  FillBlock(jac_ji, nVar, 3.0);
  FillBlock(jac_jj, nVar, 4.0);

  SECTION("SetBlocks: diagonal accumulates, off-diagonal is set") {
    matrix.SetValZero();

    matrix.SetBlocks(iEdge, iPoint, jPoint, jac_ii, jac_ij, jac_ji, jac_jj);
    CheckBlock(matrix, iPoint, iPoint, nVar, jac_ii, 1e-6);
    CheckBlock(matrix, jPoint, jPoint, nVar, jac_jj, 1e-6);
    CheckBlock(matrix, iPoint, jPoint, nVar, jac_ij, 1e-6);
    CheckBlock(matrix, jPoint, iPoint, nVar, jac_ji, 1e-6);

    /*--- A second call must double the diagonal (accumulated) and leave the
     * off-diagonal exactly as set (overwritten, not doubled). ---*/
    matrix.SetBlocks(iEdge, iPoint, jPoint, jac_ii, jac_ij, jac_ji, jac_jj);

    su2double jac_ii_2x[8][8], jac_jj_2x[8][8];
    for (auto i = 0u; i < nVar; ++i)
      for (auto j = 0u; j < nVar; ++j) {
        jac_ii_2x[i][j] = 2 * jac_ii[i][j];
        jac_jj_2x[i][j] = 2 * jac_jj[i][j];
      }
    CheckBlock(matrix, iPoint, iPoint, nVar, jac_ii_2x, 1e-6);
    CheckBlock(matrix, jPoint, jPoint, nVar, jac_jj_2x, 1e-6);
    CheckBlock(matrix, iPoint, jPoint, nVar, jac_ij, 1e-6);
    CheckBlock(matrix, jPoint, iPoint, nVar, jac_ji, 1e-6);
  }

  SECTION("SetOffDiagBlocks leaves the diagonal untouched") {
    matrix.SetValZero();
    matrix.SetBlocks(iEdge, iPoint, jPoint, jac_ii, jac_ij, jac_ji, jac_jj);

    su2double jac_ij_new[8][8], jac_ji_new[8][8];
    FillBlock(jac_ij_new, nVar, 5.0);
    FillBlock(jac_ji_new, nVar, 6.0);
    matrix.SetOffDiagBlocks(iEdge, jac_ij_new, jac_ji_new);

    CheckBlock(matrix, iPoint, iPoint, nVar, jac_ii, 1e-6);
    CheckBlock(matrix, jPoint, jPoint, nVar, jac_jj, 1e-6);
    CheckBlock(matrix, iPoint, jPoint, nVar, jac_ij_new, 1e-6);
    CheckBlock(matrix, jPoint, iPoint, nVar, jac_ji_new, 1e-6);
  }
}

TEST_CASE("SetBlocks and SetOffDiagBlocks with quantized off-diagonal storage", "[LinearAlgebra]") {
  cout.rdbuf(nullptr);

  UnitQuadTestCase testCase;
  /*--- Q_JACOBI keeps the off-diagonal blocks in int8 storage, which the block writers encode on
   * the fly instead of storing full precision values. ---*/
  testCase.AddOption("LINEAR_SOLVER_PREC= Q_JACOBI");
  testCase.InitConfig();
  testCase.InitGeometry();
  testCase.InitSolver();

  cout.rdbuf(testCase.orig_buf);

  auto* solver = testCase.solver[FLOW_SOL];
  auto& matrix = solver->Jacobian;
  const auto nVar = solver->GetnVar();
  REQUIRE(testCase.geometry->GetnEdge() > 0);

  const auto iEdge = 0ul;
  const auto iPoint = testCase.geometry->edges->GetNode(iEdge, 0);
  const auto jPoint = testCase.geometry->edges->GetNode(iEdge, 1);

  su2double jac_ii[8][8], jac_ij[8][8], jac_ji[8][8], jac_jj[8][8];
  FillBlock(jac_ii, nVar, 1.0);
  FillBlock(jac_ij, nVar, 2.0);
  FillBlock(jac_ji, nVar, 3.0);
  FillBlock(jac_jj, nVar, 4.0);

  /*--- Int8 with a per-row exponent scale, so the off-diagonal blocks come back to within about
   * one part in a hundred of the largest entry of their row; the diagonal is full precision. ---*/
  const double quantTol = 0.05;

  SECTION("SetBlocks") {
    matrix.SetValZero();
    matrix.SetBlocks(iEdge, iPoint, jPoint, jac_ii, jac_ij, jac_ji, jac_jj);

    CheckBlock(matrix, iPoint, iPoint, nVar, jac_ii, 1e-6);
    CheckBlock(matrix, jPoint, jPoint, nVar, jac_jj, 1e-6);
    CheckBlock(matrix, iPoint, jPoint, nVar, jac_ij, quantTol);
    CheckBlock(matrix, jPoint, iPoint, nVar, jac_ji, quantTol);
  }

  SECTION("SetOffDiagBlocks") {
    matrix.SetValZero();
    matrix.SetBlocks(iEdge, iPoint, jPoint, jac_ii, jac_ij, jac_ji, jac_jj);

    su2double jac_ij_new[8][8], jac_ji_new[8][8];
    FillBlock(jac_ij_new, nVar, 5.0);
    FillBlock(jac_ji_new, nVar, 6.0);
    matrix.SetOffDiagBlocks(iEdge, jac_ij_new, jac_ji_new);

    CheckBlock(matrix, iPoint, iPoint, nVar, jac_ii, 1e-6);
    CheckBlock(matrix, jPoint, jPoint, nVar, jac_jj, 1e-6);
    CheckBlock(matrix, iPoint, jPoint, nVar, jac_ij_new, quantTol);
    CheckBlock(matrix, jPoint, iPoint, nVar, jac_ji_new, quantTol);
  }
}
