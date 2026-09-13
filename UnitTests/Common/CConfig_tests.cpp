/*!
 * \file CConfig_tests.cpp
 * \brief Unit tests for option defaults that depend on other options.
 * \author SU2 Contributors
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
#include <string>
#include "../../Common/include/CConfig.hpp"

namespace {

/*--- Minimal incompressible case on a box mesh, the options under test are appended to it. ---*/
const std::string base_options =
    "SOLVER= INC_NAVIER_STOKES\n"
    "MESH_FORMAT= BOX\n"
    "MESH_BOX_SIZE= 5,5,5\n"
    "MESH_BOX_LENGTH= 1,1,1\n"
    "MESH_BOX_OFFSET= 0,0,0\n"
    "MARKER_HEATFLUX= (y_minus, 0.0, y_plus, 0.0)\n"
    "MARKER_CUSTOM= (x_minus, x_plus, z_plus, z_minus)\n"
    "INC_DENSITY_INIT= 1.1766\n"
    "INC_TEMPERATURE_INIT= 300.0\n";

INIT_OPTION_INC GetInitOptionInc(const std::string& options) {
  std::stringstream ss(base_options + options);
  auto* orig_buf = std::cout.rdbuf(nullptr);
  CConfig config(ss, SU2_COMPONENT::SU2_CFD, false);
  std::cout.rdbuf(orig_buf);
  return config.GetKind_InitOption_Inc();
}

const std::string ideal_gas_options =
    "INC_DENSITY_MODEL= VARIABLE\n"
    "INC_ENERGY_EQUATION= YES\n"
    "FLUID_MODEL= INC_IDEAL_GAS\n";

const std::string mixture_options =
    "INC_DENSITY_MODEL= VARIABLE\n"
    "INC_ENERGY_EQUATION= YES\n"
    "FLUID_MODEL= FLUID_MIXTURE\n"
    "KIND_SCALAR_MODEL= SPECIES_TRANSPORT\n"
    "SPECIES_INIT= 1.0\n"
    "MOLECULAR_WEIGHT= 28.96, 16.043\n"
    "SPECIFIC_HEAT_CP= 1009.39, 2225.0\n"
    "MARKER_SPECIES_STRONG_BC= (x_minus)\n";

}  // namespace

TEST_CASE("INIT_OPTION_INC defaults", "[Config]") {
  /*--- Without the option each fluid model keeps the density it used before INIT_OPTION_INC existed. ---*/

  CHECK(GetInitOptionInc(ideal_gas_options) == INIT_OPTION_INC::DENSITY_INIT);
  CHECK(GetInitOptionInc(mixture_options) == INIT_OPTION_INC::OPERATING_PRESSURE);
  CHECK(GetInitOptionInc("INC_DENSITY_MODEL= CONSTANT\n") == INIT_OPTION_INC::DENSITY_INIT);

  /*--- An explicit setting overrides the default. ---*/

  CHECK(GetInitOptionInc(ideal_gas_options + "INIT_OPTION_INC= OPERATING_PRESSURE\n") ==
        INIT_OPTION_INC::OPERATING_PRESSURE);
  CHECK(GetInitOptionInc(ideal_gas_options + "INIT_OPTION_INC= DENSITY_INIT\n") == INIT_OPTION_INC::DENSITY_INIT);
}
