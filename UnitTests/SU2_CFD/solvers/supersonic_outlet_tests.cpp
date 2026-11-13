/*!
 * \file supersonic_outlet_tests.cpp
 * \brief Regression for extrapolated supersonic outlet thermodynamic data.
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
#include "../../../SU2_CFD/include/variables/CNSVariable.hpp"
#include "../../../SU2_CFD/include/fluid/CPengRobinson.hpp"
#include "../../../SU2_CFD/include/numerics/flow/convection/hllc.hpp"
#include "../../../SU2_CFD/include/numerics/flow/flow_diffusion.hpp"

TEST_CASE("Supersonic outlet supplies general-gas thermodynamic data", "[Supersonic-outlet]") {
  for (const auto nDim : {2u, 3u}) {
    for (const bool realGas : {true, false}) {
      for (const bool moving : {false, true}) {
        for (const bool viscous : {false, true}) {
          UnitQuadTestCase test;
          test.config_options =
              std::string(viscous ? "SOLVER= NAVIER_STOKES\n" : "SOLVER= EULER\n") +
              "MESH_FORMAT= " + std::string(nDim == 2 ? "RECTANGLE" : "BOX") +
              "\nMESH_BOX_SIZE= 5,5,5\nREF_DIMENSIONALIZATION= DIMENSIONAL\nMACH_NUMBER= 0.1\n"
              "FREESTREAM_OPTION= DENSITY_FS\nFREESTREAM_DENSITY= 230\nFREESTREAM_PRESSURE= 1500112\n"
              "GAS_CONSTANT= 35.152\nCONV_NUM_METHOD_FLOW= HLLC\nTIME_DISCRE_FLOW= EULER_IMPLICIT\n"
              "REYNOLDS_NUMBER= 1E6\nREYNOLDS_LENGTH= 1\n"
              "VISCOSITY_MODEL= CONSTANT_VISCOSITY\nCONDUCTIVITY_MODEL= CONSTANT_CONDUCTIVITY\n"
              "GRID_MOVEMENT= STEADY_TRANSLATION\nTRANSLATION_RATE= 0,0,0\n"
              "MARKER_SUPERSONIC_OUTLET= (x_plus)\nMARKER_MONITORING= (x_plus)\nMARKER_EULER= " +
              std::string(nDim == 2 ? "(x_minus,y_minus,y_plus)\n" : "(x_minus,y_minus,y_plus,z_minus,z_plus)\n") +
              (realGas ? "FLUID_MODEL= PR_GAS\nGAMMA_VALUE= 1.0125\nCRITICAL_PRESSURE= 1415200\n"
                         "CRITICAL_TEMPERATURE= 564.1\nACENTRIC_FACTOR= 0.529\n"
                       : "FLUID_MODEL= IDEAL_GAS\nGAMMA_VALUE= 1.4\n");
          test.InitConfig();
          test.InitGeometry();
          test.InitSolver();
          std::unique_ptr<CFluidModel> gas =
              realGas ? std::unique_ptr<CFluidModel>(new CPengRobinson(1.0125, 35.152, 1415200, 564.1, .529))
                      : std::unique_ptr<CFluidModel>(new CIdealGas(1.4, 35.152));
          gas->SetLaminarViscosityModel(test.config.get());
          gas->SetThermalConductivityModel(test.config.get());
          gas->SetTDState_Prho(1500112, 230);
          const auto nVar = nDim + 2;
          const su2double grid[3] = {moving ? 20.0 : 0.0, 0, 0};
          const su2double velocity[3] = {2 * gas->GetSoundSpeed() + grid[0], 2, 0};
          const su2double density = gas->GetDensity(), pressure = gas->GetPressure();
          const auto energy = gas->GetStaticEnergy() + .5 * GeometryToolbox::SquaredNorm(nDim, velocity);
          const su2double enthalpy = energy + pressure / density;
          su2double state[5] = {density};
          for (auto i = 0u; i < nDim; ++i) state[i + 1] = density * velocity[i];
          state[nVar - 1] = density * energy;
          auto* solver = test.solver[FLOW_SOL];
          auto* nodes = static_cast<CEulerVariable*>(solver->GetNodes());
          for (auto p = 0ul; p < test.geometry->GetnPoint(); ++p) {
            nodes->SetSolution(p, state);
            nodes->SetSolution_Old(p, state);
            if (viscous) {
              REQUIRE(static_cast<CNSVariable*>(nodes)->SetPrimVar(p, 0, 0, gas.get()));
              static_cast<CNSVariable*>(nodes)->SetSecondaryVar(p, gas.get());
            } else {
              REQUIRE(nodes->SetPrimVar(p, gas.get()));
              nodes->SetSecondaryVar(p, gas.get());
            }
            test.geometry->nodes->SetGridVel(p, grid);
            for (auto i = 0u; i < nVar; ++i) solver->LinSysRes(p, i) = 0;
          }
          std::unique_ptr<CNumerics> numerics =
              realGas ? std::unique_ptr<CNumerics>(new CUpwGeneralHLLC_Flow(nDim, nVar, test.config.get()))
                      : std::unique_ptr<CNumerics>(new CUpwHLLC_Flow(nDim, nVar, test.config.get()));
          std::unique_ptr<CNumerics> diffusion;
          if (viscous) {
            if (realGas)
              diffusion.reset(new CGeneralAvgGrad_Flow(nDim, nVar, true, test.config.get()));
            else
              diffusion.reset(new CAvgGrad_Flow(nDim, nVar, true, test.config.get()));
          }
          for (auto marker = 0u; marker < test.config->GetnMarker_All(); ++marker) {
            if (test.config->GetMarker_All_TagBound(marker) != "x_plus") continue;
            solver->BC_Supersonic_Outlet(test.geometry.get(), test.solver, numerics.get(), diffusion.get(),
                                         test.config.get(), marker);
            for (auto v = 0ul; v < test.geometry->nVertex[marker]; ++v) {
              const auto p = test.geometry->vertex[marker][v]->GetNode();
              su2double normal[3], flux[5];
              test.geometry->vertex[marker][v]->GetNormal(normal);
              for (auto i = 0u; i < nDim; ++i) normal[i] = -normal[i];
              numerics->GetInviscidProjFlux(&density, velocity, &pressure, &enthalpy, normal, flux);
              numerics->SetNormal(normal);
              numerics->SetPrimitive(nodes->GetPrimitive(p), nodes->GetPrimitive(p));
              if (realGas) numerics->SetSecondary(nodes->GetSecondary(p), nodes->GetSecondary(p));
              const auto reference = numerics->ComputeResidual(test.config.get());
              for (auto i = 0u; i < nVar; ++i) {
                // The boundary must use the same extrapolated state as an interior equal-state flux.
                CHECK(solver->LinSysRes(p, i) == Approx(SU2_TYPE::GetValue(reference.residual[i])));
                if (!moving) CHECK(solver->LinSysRes(p, i) == Approx(SU2_TYPE::GetValue(flux[i])));
                for (auto j = 0u; j < nVar; ++j) CHECK(std::isfinite(solver->Jacobian.GetBlock(p, p)[i * nVar + j]));
              }
            }
          }
        }
      }
    }
  }
}
