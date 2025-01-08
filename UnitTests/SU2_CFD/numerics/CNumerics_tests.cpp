/*!
 * \file CNumerics_tests.cpp
 * \brief Unit tests for the numerics classes.
 * \author C. Pederson
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
