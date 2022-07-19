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
#include <filesystem>

#include <dirent.h>
#include <stdio.h>

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/containers/CTrapezoidalMap.hpp"
#include "../../../Common/include/containers/CLookUpTable.hpp"
#include "../../../Common/include/containers/CFileReaderLUT.hpp"

using namespace std;


TEST_CASE("LUTreader", "[tabulated chemistry]") {

  /*--- smaller and trivial lookup table ---*/

  CLookUpTable *look_up_table;


  /*--- 2D lookup table with progress variable and enthalpy as variables ---*/

  look_up_table = new CLookUpTable("src/SU2/UnitTests/Common/containers/lookuptable.drg","PROGVAR","ENTHALPY");


  /*--- string names of the controlling variables ---*/

  string name_CV1 = "PROGVAR";
  string name_CV2 = "ENTHALPY";


 /*--- look up a single value for density ---*/

  su2double prog = 0.55;
  su2double enth = -0.5; 
  string look_up_tag = "Density";
  su2double look_up_dat;
  look_up_table->LookUp_ProgEnth(look_up_tag, &look_up_dat, prog, enth, name_CV1, name_CV2); 
  CHECK(look_up_dat == Approx(1.02));

  /*--- look up a single value for viscosity ---*/

  prog = 0.65;
  enth = 0.95; 
  look_up_tag = "Viscosity";
  look_up_table->LookUp_ProgEnth(look_up_tag, &look_up_dat, prog, enth, name_CV1, name_CV2); 
  CHECK(look_up_dat == Approx(0.0000715714));


  /* find the table limits */
  
  auto limitsEnth = look_up_table->GetTableLimitsEnth();
  CHECK(limitsEnth.first == Approx(-1.0));
  CHECK(limitsEnth.second == Approx(1.0));

  auto limitsProgvar = look_up_table->GetTableLimitsProg();
  CHECK(limitsProgvar.first == Approx(0.0));
  CHECK(limitsProgvar.second == Approx(1.0));


  /* lookup value outside of lookup table */

  prog = 1.10;
  enth = 1.1;
  look_up_tag = "Density";
  look_up_table->LookUp_ProgEnth(look_up_tag, &look_up_dat, prog, enth, name_CV1, name_CV2); 
  CHECK(look_up_dat == Approx(1.2));
 
  delete look_up_table; 

}

