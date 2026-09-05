/*!
 * \file CNEMOGas_composition_tests.cpp
 * \brief Unit tests for the freestream mass fraction validation of the
 *        SU2 two-temperature thermochemistry library.
 * \author J. Bellon
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
#include <cmath>
#include <sstream>
#include <string>
#include "../../../Common/include/CConfig.hpp"
#include "../../../SU2_CFD/include/fluid/CSU2TCLib.hpp"

namespace {

/*--- Minimal dimensional NEMO configuration for a given gas model and
 * freestream composition. ---*/
std::stringstream NEMOConfigOptions(const std::string& gasModel, const std::string& composition) {
  std::stringstream cfg;
  cfg << "SOLVER= NEMO_EULER\n"
      << "MATH_PROBLEM= DIRECT\n"
      << "SYSTEM_MEASUREMENTS= SI\n"
      << "GAS_MODEL= " << gasModel << "\n"
      << "FLUID_MODEL= SU2_NONEQ\n"
      << "GAS_COMPOSITION= " << composition << "\n"
      << "FROZEN_MIXTURE= NO\n"
      << "TRANSPORT_COEFF_MODEL= WILKE\n"
      << "INIT_OPTION= TD_CONDITIONS\n"
      << "FREESTREAM_OPTION= TEMPERATURE_FS\n"
      << "REF_DIMENSIONALIZATION= DIMENSIONAL\n"
      << "MACH_NUMBER= 5\n"
      << "FREESTREAM_PRESSURE= 476.0\n"
      << "FREESTREAM_TEMPERATURE= 300.0\n"
      << "FREESTREAM_TEMPERATURE_VE= 300.0\n"
      << "CONV_NUM_METHOD_FLOW= AUSM\n"
      << "TIME_DISCRE_FLOW= EULER_IMPLICIT\n";
  return cfg;
}

/*--- Accumulate the parsed composition the same way the constructor does. ---*/
double AccumulatedMassFraction(const CConfig& config) {
  double sum = 0.0;
  const su2double* massFrac = config.GetGas_Composition();
  for (unsigned short iSpecies = 0; iSpecies < config.GetnSpecies(); ++iSpecies)
    sum += SU2_TYPE::GetValue(massFrac[iSpecies]);
  return sum;
}

}  // namespace

TEST_CASE("NEMO freestream mass fractions that sum to unity within roundoff are accepted", "[NEMO][fluid]") {
  SECTION("AIR-5 with four trace species, decimal literals sum to 1 minus one ulp") {
    auto options = NEMOConfigOptions("AIR-5", "( 0.999, 0.00025, 0.00025, 0.00025, 0.00025 )");
    CConfig config(options, SU2_COMPONENT::SU2_CFD, false);
    REQUIRE(config.GetnSpecies() == 5);

    /*--- The accumulated double precision sum is not bit-exact unity; an
     * exact equality check would reject this physically exact input. ---*/
    const double sum = AccumulatedMassFraction(config);
    CHECK(sum != 1.0);
    CHECK(std::abs(sum - 1.0) <= 1.0e-10);

    /*--- The library constructor validates the composition and exits the
     * process on rejection, so reaching the assertions is the test. ---*/
    CSU2TCLib model(&config, 3, false);
    REQUIRE(model.GetSpeciesMolarMass().size() == 5);
  }

  SECTION("N2 with a trace atomic species") {
    auto options = NEMOConfigOptions("N2", "( 0.999, 0.001 )");
    CConfig config(options, SU2_COMPONENT::SU2_CFD, false);
    REQUIRE(config.GetnSpecies() == 2);
    CHECK(std::abs(AccumulatedMassFraction(config) - 1.0) <= 1.0e-10);

    CSU2TCLib model(&config, 3, false);
    REQUIRE(model.GetSpeciesMolarMass().size() == 2);
  }
}
