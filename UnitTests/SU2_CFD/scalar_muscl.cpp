/*!
 * \file scalar_muscl.cpp
 * \brief Tests for scalar reconstruction, boundary fluxes and gradient storage.
 * \author rois1995
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
#include "../UnitQuadTestCase.hpp"
#include "../../SU2_CFD/include/solvers/CTurbSolver.hpp"
#include "../../SU2_CFD/include/numerics/scalar/scalar_edge_flux.hpp"
#include "../../SU2_CFD/include/variables/CSpeciesVariable.hpp"
#include "../../SU2_CFD/include/variables/CHeatVariable.hpp"
#include "../../SU2_CFD/include/variables/CPoissonVariable.hpp"

TEST_CASE("Edge reconstruction of scalars and of the flow", "[ScalarMuscl]") {
  auto edgeFace = [](auto value, auto edgeLimiter, su2double delta, su2double projection) {
    using Double = decltype(value);
    using Int = typename CLaneTraits<Double>::Int;
    C3DDoubleMatrix gradient(2, 1, 2);
    gradient(0, 0, 0) = projection;
    gradient(1, 0, 0) = -projection;
    Vector<Double, 2> edge;
    edge(0) = Double(1.0);
    edge(1) = Double(0.0);
    CPair<CScalarValues<Double, 1>> face;
    face.i.all(0) = Double(0.0);
    face.j.all(0) = Double(delta);
    su2activematrix unusedLimiter(2, 1);
    reconstruct<1, decltype(edgeLimiter)::value>(Int(0UL), Int(1UL), edge, gradient, unusedLimiter,
                                                 LIMITER::VAN_ALBADA_EDGE, 0, face, Double(0.0), Double(1.0));
    return face;
  };
  /*--- Scalars use the van Albada function of the point limiters, the flow keeps its own form. ---*/
  const std::integral_constant<MusclLimiter, MusclLimiter::EDGE_SCALAR> scalarEdge{};
  const std::integral_constant<MusclLimiter, MusclLimiter::EDGE> flowEdge{};
  auto scalarLimiter = [](su2double proj, su2double delta) {
    return LimiterHelpers<>::vanAlbadaFunction(proj, delta, EPS);
  };
  auto flowLimiter = [](su2double proj, su2double delta) {
    const su2double delta2 = pow(delta, 2) + 1e-6;
    return (delta2 + proj * delta) / (pow(proj, 2) + delta2);
  };
  for (const su2double delta : {0.0, 1e-8, -1e-8, 1e-3, -1e-3}) {
    for (const su2double projection : {0.0, 2e-8, -2e-8, 2e-3}) {
      const auto checkFace = [&](auto edgeLimiter, auto limiter) {
        const auto scalar = edgeFace(su2double(0.0), edgeLimiter, delta, projection);
        const auto vector = edgeFace(simd::Array<su2double, 4>(0.0), edgeLimiter, delta, projection);
        const auto left = 0.5 * projection * limiter(projection, delta);
        const auto right = delta + 0.5 * projection * limiter(-projection, delta);
        CHECK(scalar.i.all(0) == Approx(left).margin(1e-20));
        CHECK(scalar.j.all(0) == Approx(right).margin(1e-20));
        for (size_t lane = 0; lane < 4; ++lane) {
          CHECK(vector.i.all(0)[lane] == Approx(left).margin(1e-20));
          CHECK(vector.j.all(0)[lane] == Approx(right).margin(1e-20));
        }
      };
      checkFace(scalarEdge, scalarLimiter);
      checkFace(flowEdge, flowLimiter);
    }
  }
}

struct ScalarBoundaryTestSolver : CScalarSolver<CTurbVariable> {
  using CScalarSolver::BoundedScalarBCFlux;
  void LoadRestart(CGeometry**, CSolver***, CConfig*, int, bool) override {}
  ScalarBoundaryTestSolver(CGeometry* geometry, CConfig* config, const CSolver* flow, bool conservative)
      : CScalarSolver(geometry, config, flow, conservative, true) {
    nPoint = geometry->GetnPoint();
    nPointDomain = geometry->GetnPointDomain();
    nVar = 1;
    nodes = new CTurbVariable(nPoint, nDim, nVar, config);
    SetBaseClassPointerToNodes();
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  }
};

TEST_CASE("Bounded scalar boundary flux upwinds density and weights the correction", "[ScalarMuscl]") {
  UnitQuadTestCase test;
  test.InitConfig();
  test.InitGeometry();
  test.InitSolver();
  for (const bool conservative : {false, true}) {
    ScalarBoundaryTestSolver solver(test.geometry.get(), test.config.get(), test.solver[FLOW_SOL], conservative);
    solver.GetNodes()->GetSolution(0)[0] = 0.5;
    const su2double rho = 2.0, ghostRho = 5.0, normal[3] = {1.0, 0.0, 0.0};
    for (const bool implicit : {false, true}) {
      for (const su2double speed : {-3.0, 0.0, 3.0}) {
        for (const su2double weight : {0.25, 1.0}) {
          for (const bool ghost : {false, true}) {
            solver.LinSysRes.SetValZero();
            solver.Jacobian.SetValZero();
            const su2double velocity[3] = {speed, 0.0, 0.0};
            const auto mass =
                solver.BoundedScalarBCFlux(0, implicit, rho, velocity, normal, ghost ? &ghostRho : nullptr, weight);
            const auto expectedMass = speed * (ghost && speed < 0.0 ? ghostRho : rho);
            const auto q = weight * (conservative ? expectedMass : speed);
            CHECK(mass == Approx(expectedMass));
            CHECK(solver.LinSysRes(0, 0) == Approx(-0.5 * q));
            CHECK(solver.Jacobian.GetBlockView(0, 0)(0, 0) == Approx(implicit ? -q / (conservative ? rho : 1.0) : 0.0));
          }
        }
      }
    }
  }
}

TEST_CASE("Scalar reconstruction gradients follow the constructor's solver-specific flag", "[ScalarMuscl]") {
  struct Variable : CScalarVariable {
    using CScalarVariable::CScalarVariable;
    bool HasAuxiliaryGradient() const { return Gradient_Aux.size() != 0; }
  };
  for (const bool separate : {false, true}) {
    for (const bool scalarMuscl : {false, true}) {
      std::stringstream options;
      options << "SOLVER= EULER\nMUSCL_FLOW= YES\nNUM_METHOD_GRAD= GREEN_GAUSS\n"
              << "MUSCL_TURB= " << (scalarMuscl ? "YES" : "NO") << "\n"
              << "MUSCL_SPECIES= " << (scalarMuscl ? "YES" : "NO") << "\n"
              << "MUSCL_HEAT= " << (scalarMuscl ? "YES" : "NO") << "\n"
              << "NUM_METHOD_GRAD_RECON= " << (separate ? "LEAST_SQUARES" : "GREEN_GAUSS") << "\n";
      CConfig config(options, SU2_COMPONENT::SU2_CFD, false);
      REQUIRE(config.GetReconstructionGradientRequired() == separate);
      for (const bool muscl : {false, true}) {
        Variable variable(2, 2, 1, &config, muscl);
        CHECK(variable.HasAuxiliaryGradient() == (separate && muscl));
        CHECK((&variable.GetGradient_Reconstruction() != &variable.GetGradient()) == (separate && muscl));
      }
      CTurbVariable turb(2, 2, 1, &config);
      const su2double speciesInf[1] = {0.5};
      CSpeciesVariable species(speciesInf, 2, 2, 1, &config);
      CHeatVariable heat(1.0, 2, 2, 1, &config);
      CPoissonVariable poisson(0.0, 2, 2, 1, &config);
      CHECK((&turb.GetGradient_Reconstruction() != &turb.GetGradient()) == (separate && scalarMuscl));
      CHECK((&species.GetGradient_Reconstruction() != &species.GetGradient()) == (separate && scalarMuscl));
      CHECK((&heat.GetGradient_Reconstruction() != &heat.GetGradient()) == (separate && scalarMuscl));
      CHECK(&poisson.GetGradient_Reconstruction() == &poisson.GetGradient());
    }
  }
}
