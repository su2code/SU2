/*!
 * \file moving_wall_tests.cpp
 * \brief Regression for compressible wall projection in fixed and uniformly translating frames.
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
#include <array>
#include <memory>
#include "../../UnitQuadTestCase.hpp"
#include "../../../SU2_CFD/include/numerics/flow/convection/hllc.hpp"
#include "../../../SU2_CFD/include/fluid/CPengRobinson.hpp"

TEST_CASE("Compressible wall projection commutes with a Galilean boost", "[HLLC-wall]") {
  for (const auto nDim : {2u, 3u}) {
    for (const bool realGas : {false, true}) {
      UnitQuadTestCase fixed, moving;
      const std::string markers =
          nDim == 2 ? "(x_minus,x_plus,y_minus,y_plus)" : "(x_minus,x_plus,y_minus,y_plus,z_minus,z_plus)";
      const std::string options = "SOLVER= EULER\nMESH_FORMAT= " + std::string(nDim == 2 ? "RECTANGLE" : "BOX") +
                                  "\nMESH_BOX_SIZE= 5,5,5\nREF_DIMENSIONALIZATION= DIMENSIONAL\n"
                                  "MACH_NUMBER= 0.1\nFREESTREAM_OPTION= DENSITY_FS\nFREESTREAM_DENSITY= 230\n"
                                  "FREESTREAM_PRESSURE= 1500112\nGAS_CONSTANT= 35.152\n"
                                  "CONV_NUM_METHOD_FLOW= HLLC\nTIME_DISCRE_FLOW= EULER_IMPLICIT\n"
                                  "GRID_MOVEMENT= STEADY_TRANSLATION\nTRANSLATION_RATE= 0,0,0\nMARKER_EULER= " +
                                  markers + "\nMARKER_MONITORING= (x_minus)\n" +
                                  (realGas ? "FLUID_MODEL= PR_GAS\nGAMMA_VALUE= 1.0125\nCRITICAL_PRESSURE= "
                                             "1415200\nCRITICAL_TEMPERATURE= 564.1\nACENTRIC_FACTOR= 0.529\n"
                                           : "FLUID_MODEL= IDEAL_GAS\nGAMMA_VALUE= 1.4\n");
      fixed.config_options = moving.config_options = options;
      fixed.InitConfig();
      fixed.InitGeometry();
      fixed.InitSolver();
      moving.InitConfig();
      moving.InitGeometry();
      moving.InitSolver();
      const auto nVar = nDim + 2;
      const su2double boost[3] = {20, -7, 3};
      const su2double velocity[3] = {11, 4, -2};
      su2double B[5][5] = {}, T[5][5] = {}, baseJac[5][5] = {}, boostedJac[5][5] = {};
      su2double c2 = 0;
      for (auto i = 0u; i < nVar; ++i) {
        B[i][i] = T[i][i] = 1;
        for (auto j = 0u; j < nVar; ++j) baseJac[i][j] = 1 + i + 2 * j;
      }
      for (auto i = 0u; i < nDim; ++i) {
        c2 += boost[i] * boost[i];
        B[i + 1][0] = boost[i];
        T[i + 1][0] = -boost[i];
        B[nVar - 1][i + 1] = boost[i];
        T[nVar - 1][i + 1] = -boost[i];
      }
      B[nVar - 1][0] = T[nVar - 1][0] = 0.5 * c2;
      for (auto i = 0u; i < nVar; ++i)
        for (auto j = 0u; j < nVar; ++j)
          for (auto k = 0u; k < nVar; ++k)
            for (auto l = 0u; l < nVar; ++l) boostedJac[i][j] += B[i][k] * baseJac[k][l] * T[l][j];
      std::unique_ptr<CFluidModel> fluid =
          realGas ? std::unique_ptr<CFluidModel>(new CPengRobinson(1.0125, 35.152, 1415200, 564.1, 0.529))
                  : std::unique_ptr<CFluidModel>(new CIdealGas(1.4, 35.152));
      fluid->SetTDState_Prho(1500112, 230);
      std::array<su2double, 5> baseState{};
      baseState[0] = 230;
      su2double speed2 = 0;
      for (auto i = 0u; i < nDim; ++i) {
        baseState[i + 1] = 230 * velocity[i];
        speed2 += velocity[i] * velocity[i];
      }
      baseState[nVar - 1] = 230 * (fluid->GetStaticEnergy() + 0.5 * speed2);
      const su2double baseResidual[5] = {0.3, -0.5, 0.2, 0.7, -0.1};
      for (auto frame = 0u; frame < 2; ++frame) {
        auto& test = frame ? moving : fixed;
        auto* solver = test.solver[FLOW_SOL];
        auto* nodes = static_cast<CEulerVariable*>(solver->GetNodes());
        for (auto p = 0ul; p < test.geometry->GetnPoint(); ++p) {
          su2double state[5] = {}, residual[5] = {};
          for (auto i = 0u; i < nVar; ++i)
            for (auto j = 0u; j < nVar; ++j) {
              state[i] += (frame ? B[i][j] : (i == j)) * baseState[j];
              residual[i] += (frame ? B[i][j] : (i == j)) * baseResidual[j];
            }
          nodes->SetSolution(p, state);
          nodes->SetSolution_Old(p, state);
          REQUIRE(nodes->SetPrimVar(p, fluid.get()));
          nodes->SetSecondaryVar(p, fluid.get());
          for (auto i = 0u; i < nDim; ++i) test.geometry->nodes->SetGridVel(p, i, frame ? boost[i] : 0);
          for (auto i = 0u; i < nVar; ++i) {
            solver->LinSysRes(p, i) = residual[i];
            nodes->GetResTruncError(p)[i] = residual[i];
          }
          su2double block[25] = {};
          for (auto i = 0u; i < nVar; ++i)
            for (auto j = 0u; j < nVar; ++j) block[i * nVar + j] = frame ? boostedJac[i][j] : baseJac[i][j];
          solver->Jacobian.SetBlock(p, p, block);
          for (auto k = 0u; k < test.geometry->nodes->GetnPoint(p); ++k)
            solver->Jacobian.SetBlock(p, test.geometry->nodes->GetPoint(p, k), block);
        }
        std::unique_ptr<CNumerics> numerics =
            realGas ? std::unique_ptr<CNumerics>(new CUpwGeneralHLLC_Flow(nDim, nVar, test.config.get()))
                    : std::unique_ptr<CNumerics>(new CUpwHLLC_Flow(nDim, nVar, test.config.get()));
        for (auto marker = 0u; marker < test.config->GetnMarker_All(); ++marker)
          if (test.config->GetMarker_All_TagBound(marker) == "x_minus")
            solver->BC_Sym_Plane(test.geometry.get(), test.solver, numerics.get(), nullptr, test.config.get(), marker);
      }
      for (auto p = 0ul; p < fixed.geometry->GetnPoint(); ++p) {
        auto* f = fixed.solver[FLOW_SOL];
        auto* m = moving.solver[FLOW_SOL];
        for (auto i = 0u; i < nVar; ++i) {
          su2double state = 0, residual = 0, truncation = 0;
          for (auto j = 0u; j < nVar; ++j) {
            state += B[i][j] * f->GetNodes()->GetSolution_Old(p)[j];
            residual += B[i][j] * f->LinSysRes(p, j);
            truncation += B[i][j] * f->GetNodes()->GetResTruncError(p)[j];
          }
          CHECK(m->GetNodes()->GetSolution_Old(p)[i] == Approx(state).margin(1e-6));
          CHECK(m->LinSysRes(p, i) == Approx(residual).margin(1e-8));
          CHECK(m->GetNodes()->GetResTruncError(p)[i] == Approx(truncation).margin(1e-8));
          for (auto block = 0u; block <= fixed.geometry->nodes->GetnPoint(p); ++block) {
            const auto neighbor = block ? fixed.geometry->nodes->GetPoint(p, block - 1) : p;
            for (auto j = 0u; j < nVar; ++j) {
              su2double expected = 0;
              for (auto k = 0u; k < nVar; ++k)
                for (auto l = 0u; l < nVar; ++l)
                  expected += B[i][k] * f->Jacobian.GetBlock(p, neighbor)[k * nVar + l] * T[l][j];
              CHECK(m->Jacobian.GetBlock(p, neighbor)[i * nVar + j] == Approx(expected).margin(1e-5));
            }
          }
        }
      }

      for (const bool correction : {false, true}) {
        const su2double delta[5] = {0.3, -0.5, 0.2, 0.7, -0.1};
        for (auto frame = 0u; frame < 2; ++frame) {
          auto& test = frame ? moving : fixed;
          auto* solver = test.solver[FLOW_SOL];
          for (auto p = 0ul; p < test.geometry->GetnPoint(); ++p) {
            su2double state[5] = {};
            for (auto i = 0u; i < nVar; ++i)
              for (auto j = 0u; j < nVar; ++j)
                state[i] += (frame ? B[i][j] : (i == j)) * (correction ? delta[j] : baseState[j]);
            if (correction)
              solver->GetNodes()->SetSolution_Old(p, state);
            else
              solver->GetNodes()->SetSolution(p, state);
          }
          solver->MultigridProjectEulerWall(test.geometry.get(), test.config.get(), correction);
        }
        for (auto p = 0ul; p < fixed.geometry->GetnPoint(); ++p) {
          const auto* f = correction ? fixed.solver[FLOW_SOL]->GetNodes()->GetSolution_Old(p)
                                     : fixed.solver[FLOW_SOL]->GetNodes()->GetSolution(p);
          const auto* m = correction ? moving.solver[FLOW_SOL]->GetNodes()->GetSolution_Old(p)
                                     : moving.solver[FLOW_SOL]->GetNodes()->GetSolution(p);
          for (auto i = 0u; i < nVar; ++i) {
            su2double expected = 0;
            for (auto j = 0u; j < nVar; ++j) expected += B[i][j] * f[j];
            CHECK(m[i] == Approx(expected).margin(1e-6));
          }
        }
      }
    }
  }
}
