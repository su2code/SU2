/*!
 * \file CNumerics_tests.cpp
 * \brief Unit tests for the numerics classes.
 * \author C. Pederson
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
#include <sstream>
#include "../../../SU2_CFD/include/numerics/CNumerics.hpp"
#include "../../../SU2_CFD/include/numerics/turbulent/turb_sources.hpp"
#include "../../../SU2_CFD/include/variables/CEulerVariable.hpp"

TEST_CASE("NTS blending has a minimum of 0.05", "[Upwind/central blending]") {
  std::stringstream config_options;

  config_options << "SOLVER= NAVIER_STOKES" << std::endl;
  config_options << "ROE_LOW_DISSIPATION= "
                 << "NTS" << std::endl;
  config_options << "REYNOLDS_NUMBER= 5" << std::endl;

  /*--- Setup ---*/

  CConfig* config = new CConfig(config_options, SU2_COMPONENT::SU2_CFD, false);

  const su2double dissipation_i = 0;
  const su2double dissipation_j = 0;
  const su2double sensor_i = 0;
  const su2double sensor_j = 0;

  /*--- Test ---*/

  CNumerics numerics;
  const su2double dissipation_ij =
      numerics.GetRoe_Dissipation(dissipation_i, dissipation_j, sensor_i, sensor_j, config);

  REQUIRE(dissipation_ij >= 0.05);

  /*--- Teardown ---*/

  delete config;
}

TEST_CASE("SST pressure dilatation uses kinematic eddy viscosity", "[SST][Compressibility]") {
  constexpr unsigned short nDim = 2;
  using Indices = CEulerVariable::CIndices<unsigned short>;
  const Indices idx(nDim, 0);
  const std::array<su2double, 11> constants = {0.85, 1.0, 0.5, 0.856, 0.075, 0.0828, 0.09, 0.31, 5.0 / 9.0, 0.44, 10.0};
  const std::array<su2double, 2> turbulence = {2.0, 5.0};
  su2activematrix primitive_gradient(nDim + 1, nDim), turbulence_gradient(2, nDim);
  turbulence_gradient = su2double(0.0);

  for (const bool pressure_dilatation : {false, true}) {
    std::stringstream options;
    options << "SOLVER= RANS\nKIND_TURB_MODEL= SST\nREYNOLDS_NUMBER= 5e6\n"
            << "SST_OPTIONS= V2003m, COMPRESSIBILITY-SARKAR";
    if (pressure_dilatation) options << ", COMPRESSIBILITY-PRESS-DIL";
    options << '\n';
    CConfig config(options, SU2_COMPONENT::SU2_CFD, false);
    CSourcePieceWise_TurbSST<Indices> numerics(nDim, 2, constants.data(), 0.0, 0.0, &config);

    for (const su2double density : {0.5, 1.0, 2.0}) {
      std::array<su2double, nDim + 9> primitive{};
      primitive[idx.Density()] = density;
      primitive[idx.SoundSpeed()] = 10.0;
      primitive[idx.LaminarViscosity()] = 1e-5;
      primitive[idx.EddyViscosity()] = 0.4 * density;
      numerics.SetPrimitive(primitive.data(), nullptr);
      numerics.SetScalarVar(turbulence.data(), nullptr);
      numerics.SetScalarVarGradient(CMatrixView<const su2double>(turbulence_gradient), nullptr);
      numerics.SetVolume(1.0);
      numerics.SetDistance(1.0, 1.0);
      numerics.SetF1blending(0.0, 0.0);
      numerics.SetF2blending(0.0);
      numerics.SetCrossDiff(0.0);

      for (const su2double strain : {0.0, 3.0}) {
        CAPTURE(pressure_dilatation, density, strain);
        primitive_gradient = su2double(0.0);
        primitive_gradient(idx.Velocity(), 1) = strain;
        const std::array<su2double, 3> vorticity = {0.0, 0.0, -strain};
        numerics.SetPrimVarGradient(CMatrixView<const su2double>(primitive_gradient), nullptr);
        numerics.SetVorticity(vorticity.data(), nullptr);
        numerics.SetStrainMag(strain, 0.0);
        const auto residual = numerics.ComputeResidual(&config);

        // Unit-density sources use -Pi/nu_t from Tang and Gao (2015), Eq. (11), for both signs of Pi.
        const su2double expected_k =
            strain == 0.0 ? (pressure_dilatation ? -0.9288 : -0.936) : (pressure_dilatation ? 2.6136 : 2.664);
        const su2double expected_omega =
            strain == 0.0 ? (pressure_dilatation ? -1.998 : -1.98) : (pressure_dilatation ? 2.106 : 1.98);
        CHECK(residual.residual[0] == Approx(density * expected_k).epsilon(1e-12));
        CHECK(residual.residual[1] == Approx(density * expected_omega).epsilon(1e-12));
      }
    }
  }
}
