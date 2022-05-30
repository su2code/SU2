/*!
 * \file CLookupTable_tests.cpp
 * \brief Unit tests for the lookup table.
 * \author N. Beishuizen
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../Common/include/CConfig.hpp"
#include "../../../SU2_CFD/include/numerics/CTrapezoidalMap.hpp"
#include "../../../SU2_CFD/include/numerics/CLookUpTable.hpp"
#include "../../../SU2_CFD/include/numerics/CFileReaderLUT.hpp"

TEST_CASE("LUTreader", "[tabulated chemistry]") {

  CLookUpTable *look_up_table;

  /* string names of the controlling variables */
  string name_prog = "PROGVAR";
  string name_enth = "ENTHALPY";
  
  /* values at which we perform the lookup */
  su2double prog = 0.90;
  su2double enth = -220000.0;

  /* values to lookup based on column keys */
  vector<string> look_up_tags={"ViscosityDyn", "Density"};
  vector<su2double> look_up_data={0.0, 0.0}; 
  
  look_up_table = new CLookUpTable("methane_air_mixing.drg","PROGVAR","ENTHALPY");
  look_up_table->LookUp_ProgEnth(look_up_tags, look_up_data, prog,enth, name_prog, name_enth); 
  CHECK(look_up_data[0] == Approx(1.19152e-5));
  CHECK(look_up_data[1] == Approx(0.682905));

  /* value lookup based on string*/
  string look_up_tag = "Density";
  su2double look_up_dat;
  look_up_table->LookUp_ProgEnth(look_up_tag, &look_up_dat, prog,enth, name_prog, name_enth); 
  CHECK(look_up_dat == Approx(0.682905));

  auto limitsEnth = look_up_table->GetTableLimitsEnth();
  CHECK(limitsEnth.first == Approx(-1.0));
  CHECK(limitsEnth.second == Approx(1.0));

  auto limitsProgvar = look_up_table->GetTableLimitsProg();
  CHECK(limitsProgvar.first == Approx(0.0));
  CHECK(limitsProgvar.second == Approx(1.0));


}

