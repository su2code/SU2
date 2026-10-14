/*!
 * \file scalar_limiter_eps.cpp
 * \brief Tests for turbulence limiter units, zero references and derivatives.
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
#include <cmath>

namespace {
struct ScalarLimiterTestCase : UnitQuadTestCase {
  ScalarLimiterTestCase(const std::string& model, const std::string& limiter, unsigned short dimension) {
    config_options =
        "SOLVER= RANS\nKIND_TURB_MODEL= " + model + "\nSLOPE_LIMITER_TURB= " + limiter +
        "\nMUSCL_TURB= YES\nVENKAT_LIMITER_COEFF= 0.2\nMACH_NUMBER= 0.5\nREYNOLDS_NUMBER= 1e6\nMESH_FORMAT= " +
        (dimension == 2 ? "RECTANGLE" : "BOX") +
        "\nMARKER_HEATFLUX= (y_minus,0,y_plus,0)\nMESH_BOX_SIZE= 5,5,5\nMARKER_FAR= (x_minus,x_plus" +
        (dimension == 2 ? ")\n" : ",z_minus,z_plus)\n");
    InitConfig();
    InitGeometry();
    InitSolver();
    config->SetGlobalParam(config->GetKind_Solver(), RUNTIME_TURB_SYS);
  }
  ~ScalarLimiterTestCase() {
    delete solver[TURB_SOL];
    solver[TURB_SOL] = nullptr;
  }
  CTurbSolver& Turb() { return *static_cast<CTurbSolver*>(solver[TURB_SOL]); }
  su2double Evaluate(const su2double& amplitude, const su2double& reference) {
    auto* nodes = Turb().GetNodes();
    auto& gradient = nodes->GetGradient_Reconstruction();
    auto* ref = const_cast<su2double*>(Turb().GetLimiterReference());
    for (auto var = 0u; var < Turb().GetnVar(); ++var) ref[var] = reference * (var + 1.0);
    for (auto point = 0ul; point < geometry->GetnPoint(); ++point) {
      const auto x = geometry->nodes->GetCoord(point, 0);
      for (auto var = 0u; var < Turb().GetnVar(); ++var) {
        nodes->GetSolution(point)[var] = amplitude * (1.0 + var + x * x);
        for (auto dim = 0u; dim < geometry->GetnDim(); ++dim) {
          su2double slope = 0.0;
          if (dim == 0) slope = 4.0 + 2.0 * x;
          gradient(point, var, dim) = amplitude * slope;
        }
      }
    }
    SU2_OMP_PARALLEL { Turb().SetSolution_Limiter(geometry.get(), config.get()); }
    su2double result = 0.0;
    for (auto point = 0ul; point < geometry->GetnPointDomain(); ++point)
      for (auto var = 0u; var < Turb().GetnVar(); ++var) result += nodes->GetLimiter()(point, var);
    return result;
  }
  std::vector<double> Snapshot() {
    std::vector<double> values;
    for (auto point = 0ul; point < geometry->GetnPointDomain(); ++point)
      for (auto var = 0u; var < Turb().GetnVar(); ++var)
        values.push_back(SU2_TYPE::GetValue(Turb().GetNodes()->GetLimiter()(point, var)));
    return values;
  }
};
const std::vector<std::string> scaledLimiters = {"VENKATAKRISHNAN", "NISHIKAWA_R3", "NISHIKAWA_R4", "NISHIKAWA_R5"};
}  // namespace

TEST_CASE("Turbulence point limiters preserve their values across field units", "[ScalarLimiterEps]") {
  for (const std::string model : {"SA", "SST"}) {
    for (const unsigned short dim : {2, 3}) {
      for (const auto& kind : scaledLimiters) {
        ScalarLimiterTestCase test(model, kind, dim);
        CHECK(test.solver[FLOW_SOL]->GetLimiterReference() == nullptr);
        test.Evaluate(1.0, 10.0);
        const auto reference = test.Snapshot();
        REQUIRE_FALSE(reference.empty());
        for (const double units : {1e-12, 1e-6, 1.0, 1e6, 1e12}) {
          for (const double sign : {-1.0, 1.0}) {
            CAPTURE(model, dim, kind, units, sign);
            test.Evaluate(units, sign * 10.0 * units);
            const auto values = test.Snapshot();
            REQUIRE(values.size() == reference.size());
            for (size_t i = 0; i < values.size(); ++i) CHECK(values[i] == Approx(reference[i]).margin(1e-12));
          }
        }
        test.Evaluate(1.0, 0.0);
        const auto zeroReference = test.Snapshot();
        test.Evaluate(1.0, 1e-300);
        const auto tinyReference = test.Snapshot();
        for (size_t i = 0; i < zeroReference.size(); ++i) {
          CHECK(std::isfinite(zeroReference[i]));
          CHECK(tinyReference[i] == Approx(zeroReference[i]).margin(1e-12));
        }
        test.Evaluate(0.0, 0.0);
        for (const auto value : test.Snapshot()) CHECK(value == Approx(1.0));
      }
    }
  }
}

TEST_CASE("Wang's range-relative turbulence limiter ignores the reference magnitude", "[ScalarLimiterEps]") {
  for (const unsigned short dim : {2, 3}) {
    ScalarLimiterTestCase test("SST", "VENKATAKRISHNAN_WANG", dim);
    test.Evaluate(1.0, 10.0);
    const auto reference = test.Snapshot();
    for (const double value : {0.0, -10.0, 1e-300, 1e12}) {
      test.Evaluate(1.0, value);
      const auto values = test.Snapshot();
      for (size_t i = 0; i < values.size(); ++i) CHECK(values[i] == Approx(reference[i]).margin(1e-12));
    }
  }
}

#if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)
TEST_CASE("Scalar limiter reference derivatives match finite differences", "[ScalarLimiterEps][AD]") {
  for (const auto& kind : scaledLimiters) {
    AD::Reset();
    {
      ScalarLimiterTestCase test("SST", kind, 2);
      su2double amplitude = 1.0, reference = 10.0;
      double dAmplitude, dReference;
#ifdef CODI_REVERSE_TYPE
      AD::StartRecording();
      AD::RegisterInput(amplitude);
      AD::RegisterInput(reference);
      su2double result = test.Evaluate(amplitude, reference);
      AD::RegisterOutput(result);
      AD::StopRecording();
      SU2_TYPE::SetDerivative(result, 1.0);
      AD::ComputeAdjoint();
      dAmplitude = SU2_TYPE::GetDerivative(amplitude);
      dReference = SU2_TYPE::GetDerivative(reference);
#else
      SU2_TYPE::SetDerivative(amplitude, 1.0);
      SU2_TYPE::SetDerivative(reference, 0.0);
      dAmplitude = SU2_TYPE::GetDerivative(test.Evaluate(amplitude, reference));
      SU2_TYPE::SetDerivative(amplitude, 0.0);
      SU2_TYPE::SetDerivative(reference, 1.0);
      dReference = SU2_TYPE::GetDerivative(test.Evaluate(amplitude, reference));
#endif
      const double h = 1e-4;
      const double plusAmplitude = SU2_TYPE::GetValue(test.Evaluate(su2double(1.0 + h), su2double(10.0)));
      const double minusAmplitude = SU2_TYPE::GetValue(test.Evaluate(su2double(1.0 - h), su2double(10.0)));
      const double plusReference = SU2_TYPE::GetValue(test.Evaluate(su2double(1.0), su2double(10.0 + h)));
      const double minusReference = SU2_TYPE::GetValue(test.Evaluate(su2double(1.0), su2double(10.0 - h)));
      CHECK(dAmplitude == Approx((plusAmplitude - minusAmplitude) / (2 * h)).epsilon(1e-5).margin(1e-7));
      CHECK(dReference == Approx((plusReference - minusReference) / (2 * h)).epsilon(1e-5).margin(1e-7));
      CHECK(std::abs(dReference) > 1e-8);
      CHECK(dAmplitude + 10.0 * dReference == Approx(0.0).margin(1e-8));
    }
    AD::Reset();
  }
}
#endif
