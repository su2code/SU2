/*!
 * \file periodic_ad_tests.cpp
 * \brief Finite-difference check of the periodic linear-solve external adjoint.
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

TEST_CASE("Periodic linear solve external adjoint", "[Periodic][AD tests]") {
  const bool rotation = GENERATE(false, true);
  AD::Reset();
  UnitQuadTestCase field;
  const auto start = field.config_options.find("MARKER_HEATFLUX=");
  const auto end = field.config_options.find("VISCOSITY_MODEL=");
  field.config_options.replace(start, end - start, "MARKER_CUSTOM= (y_minus,y_plus,z_plus,z_minus)\n");
  field.SetOption("SOLVER= EULER");
  field.SetOption("KIND_VERIFICATION_SOLUTION= NO_VERIFICATION_SOLUTION");
  field.AddOption("MATH_PROBLEM= DISCRETE_ADJOINT");
  field.AddOption(rotation ? "MARKER_PERIODIC= (x_minus,x_plus, 0,0.5,0.5, 90,0,0, 1,0,0)"
                           : "MARKER_PERIODIC= (x_minus,x_plus, 0,0,0, 0,0,0, 1,0,0)");
  field.AddOption("LINEAR_SOLVER= FGMRES\nLINEAR_SOLVER_PREC= JACOBI\nLINEAR_SOLVER_ERROR= 1e-12");
  field.AddOption("LINEAR_SOLVER_ITER= 200\nDISCADJ_LIN_SOLVER= FGMRES\nDISCADJ_LIN_PREC= JACOBI");
  field.InitConfig();
  field.InitGeometry(true);
  field.geometry->MatchPeriodic(field.config.get(), 1);
  field.geometry->PreprocessPeriodicComms(field.geometry.get(), field.config.get());
  const auto nPoint = field.geometry->GetnPoint(), nDomain = field.geometry->GetnPointDomain();
  REQUIRE(nDomain > 0);
  constexpr unsigned short nVar = 5;
  CSysMatrix<su2mixedfloat> matrix;
  matrix.Initialize(nPoint, nDomain, nVar, nVar, true, field.geometry.get(), field.config.get());
  matrix.SetPeriodicProjection(1);
  matrix.SetValZero();
  for (auto i = 0ul; i < nPoint; ++i) {
    std::vector<unsigned long> columns = {i};
    for (auto j : field.geometry->nodes->GetPoints(i)) columns.push_back(j);
    for (auto j : columns) {
      auto* block = matrix.GetBlock(i, j);
      for (auto a = 0u; a < nVar; ++a)
        for (auto b = 0u; b < nVar; ++b)
          block[a * nVar + b] = (i == j && a == b ? 3.0 : 0.0) + 0.001 * (1 + a + 2 * b) +
                                (i == j ? 0 : 0.002 * (1 + field.geometry->nodes->GetGlobalIndex(i)));
    }
  }
  CSysSolve<su2mixedfloat> system;
  CSysVector<su2double> rhs(nPoint, nDomain, nVar), solution(nPoint, nDomain, nVar);
  auto solve = [&](const su2double& parameter) {
    rhs = su2double(0);
    solution = su2double(0);
    for (auto i = 0ul; i < nDomain; ++i)
      for (auto a = 0u; a < nVar; ++a) {
        const auto global = field.geometry->nodes->GetGlobalIndex(i) * nVar + a;
        rhs(i, a) = cos(0.03 * global) + parameter * sin(0.1 * global);
      }
    system.Solve(matrix, rhs, solution, field.geometry.get(), field.config.get());
    su2double objective = 0;
    for (auto i = 0ul; i < nDomain; ++i)
      for (auto a = 0u; a < nVar; ++a)
        objective += cos(0.07 * (field.geometry->nodes->GetGlobalIndex(i) * nVar + a)) * solution(i, a);
    return objective;
  };
  su2double parameter = 0.2;
  AD::StartRecording();
  AD::RegisterInput(parameter);
  auto objective = solve(parameter);
  AD::RegisterOutput(objective);
  AD::StopRecording();
  SU2_TYPE::SetDerivative(objective, 1.0);
  AD::ComputeAdjoint();
  su2double localDerivative = SU2_TYPE::GetDerivative(parameter), derivative = 0;
  SU2_MPI::Allreduce(&localDerivative, &derivative, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  /*--- The external solve transposes its matrix during the reverse sweep. ---*/
  matrix.TransposeInPlace();
  AD::Reset();
  constexpr double step = 1e-5;
  const auto plus = SU2_TYPE::GetValue(solve(su2double(0.2 + step)));
  const auto minus = SU2_TYPE::GetValue(solve(su2double(0.2 - step)));
  su2double localDifference = (plus - minus) / (2 * step), difference = 0;
  SU2_MPI::Allreduce(&localDifference, &difference, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  CHECK(SU2_TYPE::GetValue(derivative) == Approx(SU2_TYPE::GetValue(difference)).epsilon(1e-6).margin(1e-6));
  AD::Reset();
}
