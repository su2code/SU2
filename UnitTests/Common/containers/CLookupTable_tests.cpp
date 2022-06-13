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
#include "../../../Common/include/containers/CTrapezoidalMap.hpp"
#include "../../../Common/include/containers/CLookUpTable.hpp"
#include "../../../Common/include/containers/CFileReaderLUT.hpp"

using namespace std;

TEST_CASE("LUTreader", "[tabulated chemistry]") {

  CLookUpTable *look_up_table1, *look_up_table2;

  /* string names of the controlling variables */
  string name_prog = "PROGVAR";
  string name_enth = "ENTHALPY";
  
  /* values at which we perform the lookup */
  su2double prog = 0.90;
  su2double enth = -220000.0;

  /* values to lookup based on column keys */
  vector<string> look_up_tags={"ViscosityDyn", "Density"};
  vector<su2double> look_up_data={0.0, 0.0}; 
  
  look_up_table1 = new CLookUpTable("../UnitTests/Common/containers/methane_air_mixing.drg","PROGVAR","ENTHALPY");

  look_up_table1->LookUp_ProgEnth(look_up_tags, look_up_data, prog,enth, name_prog, name_enth); 
  cout << "check 1" << endl;
  CHECK(look_up_data[0] == Approx(1.19152e-5));
  cout << "check 2" << endl;
  CHECK(look_up_data[1] == Approx(0.682905));

  /* value lookup based on string*/
  string look_up_tag = "Density";
  su2double look_up_dat;
  look_up_table1->LookUp_ProgEnth(look_up_tag, &look_up_dat, prog,enth, name_prog, name_enth); 
  cout << "check 3" << endl;
  CHECK(look_up_dat == Approx(0.682905));

  // find the table limits
  auto limitsEnth = look_up_table1->GetTableLimitsEnth();
  cout << "check 4" << endl;
  CHECK(limitsEnth.first == Approx(-1.0));
  cout << "check 5" << endl;
  CHECK(limitsEnth.second == Approx(1.0));

  auto limitsProgvar = look_up_table1->GetTableLimitsProg();
  cout << "check 6" << endl;
  CHECK(limitsProgvar.first == Approx(0.0));
  cout << "check 7" << endl;
  CHECK(limitsProgvar.second == Approx(1.0));


  cout << "checking value outside lookup table" << endl;
  /* lookup value outside of lookup table */
  prog = 1.10;
  enth = -220000.0;
  look_up_tag = "Density";
  look_up_table1->LookUp_ProgEnth(look_up_tag, &look_up_dat, prog, enth, name_prog, name_enth); 
  cout << "check 8: " << look_up_dat << endl;
  CHECK(look_up_dat == Approx(0.6516888435));
 
  delete look_up_table1; 


  look_up_table2 = new CLookUpTable("../UnitTests/Common/containers/lookuptable.drg","PROGVAR","ENTHALPY");

  prog = 0.55;
  enth = 0.25; 
  look_up_tag = "Density";
  look_up_table2->LookUp_ProgEnth(look_up_tag, &look_up_dat, prog,enth, name_prog, name_enth); 
  cout << "check 9" << endl;
  CHECK(look_up_dat == Approx(1.05));


  prog = 0.65;
  enth = 0.95; 
  look_up_tag = "Density";
  look_up_table2->LookUp_ProgEnth(look_up_tag, &look_up_dat, prog,enth, name_prog, name_enth); 
  cout << "check 10" << endl;
  CHECK(look_up_dat == Approx(1.19));

  delete look_up_table2; 

}

