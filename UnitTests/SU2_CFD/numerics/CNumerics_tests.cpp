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
#include <sstream>
#include "../../../SU2_CFD/include/numerics/CNumerics.hpp"

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

TEST_CASE("QCR2000 corrects only the turbulent stress", "[QCR]") {
  constexpr size_t nDim = 3;
  su2double gradvel[3][3] = {{0.0}};
  gradvel[0][1] = 100.0;
  const su2double mu_lam = 2.0;
  const su2double mu_turb = 3.0;

  /*--- At a no-slip wall the eddy viscosity is zero, so the correction
   * must leave the (purely laminar) stress tensor unchanged. ---*/

  su2double tau_wall[3][3] = {{0.0}}, tau_wall_ref[3][3] = {{0.0}};
  CNumerics::ComputeStressTensor(nDim, tau_wall, gradvel, mu_lam);
  CNumerics::ComputeStressTensor(nDim, tau_wall_ref, gradvel, mu_lam);
  CNumerics::AddQCR(nDim, gradvel, tau_wall, 0.0);
  for (size_t iDim = 0; iDim < nDim; iDim++)
    for (size_t jDim = 0; jDim < nDim; jDim++) REQUIRE(tau_wall[iDim][jDim] == tau_wall_ref[iDim][jDim]);

  /*--- Scaling the correction of the total stress by mu_t / (mu_l + mu_t)
   * is equivalent to correcting the turbulent stress tensor alone. ---*/

  su2double tau_total[3][3] = {{0.0}}, tau_turb[3][3] = {{0.0}}, tau_lam[3][3] = {{0.0}};
  CNumerics::ComputeStressTensor(nDim, tau_total, gradvel, mu_lam + mu_turb);
  CNumerics::ComputeStressTensor(nDim, tau_turb, gradvel, mu_turb);
  CNumerics::ComputeStressTensor(nDim, tau_lam, gradvel, mu_lam);
  CNumerics::AddQCR(nDim, gradvel, tau_total, mu_turb / (mu_lam + mu_turb));
  CNumerics::AddQCR(nDim, gradvel, tau_turb, 1.0);
  for (size_t iDim = 0; iDim < nDim; iDim++)
    for (size_t jDim = 0; jDim < nDim; jDim++)
      REQUIRE(tau_total[iDim][jDim] == Approx(tau_lam[iDim][jDim] + tau_turb[iDim][jDim]).margin(1e-12));
}
