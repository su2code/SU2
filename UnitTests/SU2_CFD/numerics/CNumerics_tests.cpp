/*!
 * \file CNumerics_tests.cpp
 * \brief Unit tests for the numerics classes.
 * \author C. Pederson
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../SU2_CFD/include/numerics_structure.hpp"

static void WriteCfgFile(unsigned short nDim, const char* filename,
                  std::string blending) {
  std::ofstream cfg_file;

  cfg_file.open(filename, ios::out);
  cfg_file << "SOLVER= NAVIER_STOKES" << std::endl;
  cfg_file << "ROE_LOW_DISSIPATION= " << blending << std::endl;
  cfg_file << "REYNOLDS_NUMBER= 5" << std::endl;
  cfg_file.close();
}

TEST_CASE("NTS blending has a minimum of 0.05", "[Upwind/central blending]") {

  /*--- Setup ---*/

  const unsigned short nDim = 3;

  char cfg_filename[100] = "convective_blending_test.cfg";
  WriteCfgFile(nDim, cfg_filename, "NTS");
  CConfig* config = new CConfig(cfg_filename, SU2_CFD, false);
  std::remove(cfg_filename);

  const su2double dissipation_i = 0;
  const su2double dissipation_j = 0;
  const su2double sensor_i = 0;
  const su2double sensor_j = 0;
  su2double dissipation_ij;

  /*--- Test ---*/

  CNumerics numerics;
  numerics.SetRoe_Dissipation(dissipation_i, dissipation_j,
                              sensor_i, sensor_j,
                              dissipation_ij, config);

  REQUIRE(dissipation_ij >= 0.05);

  /*--- Teardown ---*/

  delete config;
}

