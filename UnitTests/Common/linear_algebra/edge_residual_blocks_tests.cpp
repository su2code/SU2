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
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
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

TEST_CASE("Complete periodic implicit operator and transpose", "[Periodic][LinearAlgebra]") {
  const auto kind = GENERATE(0u, 1u, 2u, 3u);
  const bool rotation = kind == 1;
  const auto nPairs = kind > 1 ? kind : 1u;
  UnitQuadTestCase field;
  const auto start = field.config_options.find("MARKER_HEATFLUX=");
  const auto end = field.config_options.find("VISCOSITY_MODEL=");
  field.config_options.replace(start, end - start,
                               nPairs == 1 ? "MARKER_CUSTOM= (y_minus,y_plus,z_plus,z_minus)\n"
                                           : (nPairs == 2 ? "MARKER_CUSTOM= (z_plus,z_minus)\n" : ""));
  std::string periodic = rotation ? "MARKER_PERIODIC= (x_minus,x_plus, 0,0.5,0.5, 90,0,0, 1,0,0"
                                  : "MARKER_PERIODIC= (x_minus,x_plus, 0,0,0, 0,0,0, 1,0,0";
  if (nPairs > 1) periodic += ", y_minus,y_plus, 0,0,0, 0,0,0, 0,1,0";
  if (nPairs > 2) periodic += ", z_minus,z_plus, 0,0,0, 0,0,0, 0,0,1";
  field.AddOption(periodic + ")");
  /*--- Avoid asking a float Krylov solver to converge below roundoff. ---*/
  field.AddOption("LINEAR_SOLVER_PREC= JACOBI\nLINEAR_SOLVER_ITER= 150");
  field.AddOption(sizeof(su2mixedfloat) == sizeof(float) ? "LINEAR_SOLVER_ERROR= 1e-7" : "LINEAR_SOLVER_ERROR= 1e-10");
  field.InitConfig();
  field.InitGeometry(true);
  for (auto pair = 1u; pair <= nPairs; ++pair) field.geometry->MatchPeriodic(field.config.get(), pair);
  field.geometry->PreprocessPeriodicComms(field.geometry.get(), field.config.get());
  field.InitSolver();
  auto& matrix = field.solver[FLOW_SOL]->Jacobian;
  const auto nVar = field.solver[FLOW_SOL]->GetnVar();
  const auto nPoint = field.geometry->GetnPoint();
  const auto nDomain = field.geometry->GetnPointDomain();
  unsigned long globalPoints = 0;
  SU2_MPI::Allreduce(&nDomain, &globalPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  REQUIRE(nDomain > 0);
  const auto size = globalPoints * nVar;
  std::vector<su2double> partialA(size * size, 0), partialP(size * size, 0), A(size * size), P(size * size);
  matrix.SetValZero();
  for (auto i = 0ul; i < nPoint; ++i) {
    const auto globalI = field.geometry->nodes->GetGlobalIndex(i);
    std::vector<unsigned long> columns = {i};
    for (auto j : field.geometry->nodes->GetPoints(i)) columns.push_back(j);
    for (auto j : columns) {
      const auto globalJ = field.geometry->nodes->GetGlobalIndex(j);
      auto* block = matrix.GetBlock(i, j);
      REQUIRE(block != nullptr);
      for (auto a = 0u; a < nVar; ++a)
        for (auto b = 0u; b < nVar; ++b) {
          /*--- Nonsymmetric blocks expose mistakes in the reverse operator. ---*/
          const auto value =
              (i == j && a == b ? 3.0 : 0.0) + 0.001 * (1 + a + 2 * b) + (i == j ? 0 : 0.002 * (1 + globalI));
          block[a * nVar + b] = value;
          if (i < nDomain)
            partialA[(globalI * nVar + a) * size + globalJ * nVar + b] = SU2_TYPE::GetValue(block[a * nVar + b]);
        }
    }
    if (i < nDomain)
      for (auto a = 0u; a < nVar; ++a) partialP[(globalI * nVar + a) * size + globalI * nVar + a] = 1;
  }
  SU2_MPI::Allreduce(partialA.data(), A.data(), A.size(), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(partialP.data(), P.data(), P.size(), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  for (auto pair = 1u; pair <= nPairs; ++pair) {
    std::fill(partialP.begin(), partialP.end(), 0);
    for (auto i = 0ul; i < nDomain; ++i) {
      const auto globalI = field.geometry->nodes->GetGlobalIndex(i);
      for (auto a = 0u; a < nVar; ++a)
        std::copy_n(P.data() + (globalI * nVar + a) * size, size, partialP.data() + (globalI * nVar + a) * size);
    }
    for (auto marker = 0u; marker < field.geometry->GetnMarker(); ++marker) {
      if (field.config->GetMarker_All_KindBC(marker) != PERIODIC_BOUNDARY) continue;
      const auto index = static_cast<unsigned short>(field.config->GetMarker_All_PerBound(marker));
      if (index != pair && index != pair + nPairs) continue;
      const auto* angles = field.config->GetPeriodicRotAngles(field.config->GetMarker_All_TagBound(marker));
      su2double q[3][3];
      GeometryToolbox::RotationMatrix(angles[0], angles[1], angles[2], q);
      for (auto vertex = 0ul; vertex < field.geometry->GetnVertex(marker); ++vertex) {
        const auto* point = field.geometry->vertex[marker][vertex];
        const auto i = point->GetNode();
        if (i >= nDomain) continue;
        const auto globalI = field.geometry->nodes->GetGlobalIndex(i);
        const auto donor = point->GetDonorGlobalIndex();
        for (auto a = 0u; a < nVar; ++a)
          for (auto j = 0ul; j < size; ++j) {
            auto value = P[(globalI * nVar + a) * size + j];
            if (a >= 1 && a <= 3) {
              for (auto b = 1u; b <= 3; ++b) value += q[b - 1][a - 1] * P[(donor * nVar + b) * size + j];
            } else {
              value += P[(donor * nVar + a) * size + j];
            }
            partialP[(globalI * nVar + a) * size + j] = 0.5 * value;
          }
      }
    }
    SU2_MPI::Allreduce(partialP.data(), P.data(), P.size(), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  }
  auto product = [&](const std::vector<su2double>& mat, const std::vector<su2double>& x, bool transpose = false) {
    std::vector<su2double> result(size, 0);
    for (auto i = 0ul; i < size; ++i)
      for (auto j = 0ul; j < size; ++j) result[i] += (transpose ? mat[j * size + i] : mat[i * size + j]) * x[j];
    return result;
  };
  std::vector<su2double> x(size), y(size);
  for (auto i = 0ul; i < size; ++i) {
    x[i] = sin(0.1 * i);
    y[i] = cos(0.07 * i);
  }
  const auto px = product(P, x), py = product(P, y);
  auto reference = product(P, product(A, px));
  auto transposeReference = product(P, product(A, py, true));
  for (auto i = 0ul; i < size; ++i) {
    reference[i] += x[i] - px[i];
    transposeReference[i] += y[i] - py[i];
  }
  CSysVector<su2mixedfloat> input(nPoint, nDomain, nVar), output(nPoint, nDomain, nVar),
      transposeOutput(nPoint, nDomain, nVar);
  for (auto i = 0ul; i < nPoint; ++i)
    for (auto a = 0u; a < nVar; ++a) input(i, a) = x[field.geometry->nodes->GetGlobalIndex(i) * nVar + a];
  SU2_OMP_PARALLEL { matrix.MatrixVectorProduct(input, output, field.geometry.get(), field.config.get()); }
  END_SU2_OMP_PARALLEL
  su2double error = 0;
  for (auto i = 0ul; i < nDomain; ++i)
    for (auto a = 0u; a < nVar; ++a)
      error = std::max(error, fabs(SU2_TYPE::GetValue(output(i, a)) -
                                   reference[field.geometry->nodes->GetGlobalIndex(i) * nVar + a]));
  CHECK(error < 1e-5);
  /*--- Krylov solution accuracy depends on the arithmetic used by its basis. ---*/
  const auto solveTolerance =
      sizeof(su2mixedfloat) == sizeof(float) ? sqrt(std::numeric_limits<float>::epsilon()) : 1e-5;
  CSysSolve<su2mixedfloat> system;
  CSysVector<su2double> rhs(nPoint, nDomain, nVar), solution(nPoint, nDomain, nVar);
  solution = su2double(0);
  for (auto i = 0ul; i < nPoint; ++i)
    for (auto a = 0u; a < nVar; ++a) rhs(i, a) = reference[field.geometry->nodes->GetGlobalIndex(i) * nVar + a];
  SU2_OMP_PARALLEL { system.Solve(matrix, rhs, solution, field.geometry.get(), field.config.get()); }
  END_SU2_OMP_PARALLEL
  error = 0;
  for (auto i = 0ul; i < nDomain; ++i)
    for (auto a = 0u; a < nVar; ++a)
      error = std::max(
          error, fabs(SU2_TYPE::GetValue(solution(i, a)) - x[field.geometry->nodes->GetGlobalIndex(i) * nVar + a]));
  CHECK(error < solveTolerance);
  SU2_OMP_PARALLEL { matrix.TransposeInPlace(); }
  END_SU2_OMP_PARALLEL
  for (auto i = 0ul; i < nPoint; ++i)
    for (auto a = 0u; a < nVar; ++a) input(i, a) = y[field.geometry->nodes->GetGlobalIndex(i) * nVar + a];
  SU2_OMP_PARALLEL { matrix.MatrixVectorProduct(input, transposeOutput, field.geometry.get(), field.config.get()); }
  END_SU2_OMP_PARALLEL
  error = 0;
  su2double dotForward = 0, dotTranspose = 0;
  for (auto i = 0ul; i < nDomain; ++i)
    for (auto a = 0u; a < nVar; ++a) {
      const auto global = field.geometry->nodes->GetGlobalIndex(i) * nVar + a;
      error = std::max(error, fabs(SU2_TYPE::GetValue(transposeOutput(i, a)) - transposeReference[global]));
      dotForward += y[global] * SU2_TYPE::GetValue(output(i, a));
      dotTranspose += x[global] * SU2_TYPE::GetValue(transposeOutput(i, a));
    }
  CHECK(error < 1e-5);
  su2double dots[] = {dotForward, dotTranspose}, globalDots[2] = {};
  SU2_MPI::Allreduce(dots, globalDots, 2, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  CHECK(globalDots[0] == Approx(globalDots[1]).margin(1e-4));
  SU2_OMP_PARALLEL { matrix.TransposeInPlace(); }
  END_SU2_OMP_PARALLEL
  solution = su2double(0);
  for (auto i = 0ul; i < nPoint; ++i)
    for (auto a = 0u; a < nVar; ++a)
      rhs(i, a) = transposeReference[field.geometry->nodes->GetGlobalIndex(i) * nVar + a];
  SU2_OMP_PARALLEL { system.Solve_b(matrix, rhs, solution, field.geometry.get(), field.config.get()); }
  END_SU2_OMP_PARALLEL
  error = 0;
  for (auto i = 0ul; i < nDomain; ++i)
    for (auto a = 0u; a < nVar; ++a)
      error = std::max(
          error, fabs(SU2_TYPE::GetValue(solution(i, a)) - y[field.geometry->nodes->GetGlobalIndex(i) * nVar + a]));
  CHECK(error < solveTolerance);
}
