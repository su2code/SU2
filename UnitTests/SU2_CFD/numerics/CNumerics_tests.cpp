#include "catch.hpp"
#include <sstream>
#include "../../../SU2_CFD/include/numerics_structure.hpp"

TEST_CASE("NTS blending has a minimum of 0.05", "[Upwind/central blending]") {

  std::stringstream config_options;
  
  config_options << "SOLVER= NAVIER_STOKES" << std::endl;
  config_options << "ROE_LOW_DISSIPATION= " << "NTS" << std::endl;
  config_options << "REYNOLDS_NUMBER= 5" << std::endl;
  
  /*--- Setup ---*/

  const unsigned short nDim = 3;

  CConfig* config = new CConfig(config_options, SU2_CFD, false);

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

