/*!
 * \file CLookupTable_tests.cpp
 * \brief Unit tests for the lookup table.
 * \author N. Beishuizen
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
#include <stdio.h>

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/containers/CTrapezoidalMap.hpp"
#include "../../../Common/include/containers/CLookUpTable.hpp"
#include "../../../Common/include/containers/CFileReaderLUT.hpp"

TEST_CASE("LUTreader", "[tabulated chemistry]") {
  /*--- smaller and trivial lookup table ---*/

  CLookUpTable look_up_table("src/SU2/UnitTests/Common/containers/lookuptable.drg", "ProgressVariable", "EnthalpyTot");

  /*--- string names of the controlling variables ---*/

  string name_CV1 = "ProgressVariable";
  string name_CV2 = "EnthalpyTot";

  /*--- look up a single value for density ---*/

  su2double prog = 0.55;
  su2double enth = -0.5;
  string look_up_tag = "Density";
  su2double look_up_dat;
  look_up_table.LookUp_XY(look_up_tag, &look_up_dat, prog, enth);
  CHECK(look_up_dat == Approx(1.02));

  /*--- look up a single value for viscosity ---*/

  prog = 0.6;
  enth = 0.9;
  look_up_tag = "Viscosity";
  look_up_table.LookUp_XY(look_up_tag, &look_up_dat, prog, enth);
  CHECK(look_up_dat == Approx(0.0000674286));

  /* find the table limits */

  auto limitsEnth = look_up_table.GetTableLimitsY();
  CHECK(SU2_TYPE::GetValue(*limitsEnth.first) == Approx(-1.0));
  CHECK(SU2_TYPE::GetValue(*limitsEnth.second) == Approx(1.0));

  auto limitsProgvar = look_up_table.GetTableLimitsX();
  CHECK(SU2_TYPE::GetValue(*limitsProgvar.first) == Approx(0.0));
  CHECK(SU2_TYPE::GetValue(*limitsProgvar.second) == Approx(1.0));

  /* lookup value outside of lookup table */

  prog = 1.10;
  enth = 1.1;
  look_up_tag = "Density";
  look_up_table.LookUp_XY(look_up_tag, &look_up_dat, prog, enth);
  CHECK(look_up_dat == Approx(1.1738796125));
}

TEST_CASE("LUTreader_3D", "[tabulated chemistry]") {
  /*--- smaller and trivial lookup table ---*/

  CLookUpTable look_up_table("src/SU2/UnitTests/Common/containers/lookuptable_3D.drg", "ProgressVariable",
                             "EnthalpyTot");

  /*--- string names of the controlling variables ---*/

  string name_CV1 = "ProgressVariable";
  string name_CV2 = "EnthalpyTot";

  /*--- look up a single value for density ---*/

  su2double prog = 0.55;
  su2double enth = -0.5;
  su2double mfrac = 0.5;
  string look_up_tag = "Density";
  su2double look_up_dat;
  look_up_table.LookUp_XYZ(look_up_tag, &look_up_dat, prog, enth, mfrac);
  CHECK(look_up_dat == Approx(1.02));

  /*--- look up a single value for viscosity ---*/

  prog = 0.6;
  enth = 0.9;
  mfrac = 0.8;
  look_up_tag = "Viscosity";
  look_up_table.LookUp_XYZ(look_up_tag, &look_up_dat, prog, enth, mfrac);
  CHECK(look_up_dat == Approx(0.0000674286));

  /* find the table limits */

  auto limitsEnth = look_up_table.GetTableLimitsY();
  CHECK(SU2_TYPE::GetValue(*limitsEnth.first) == Approx(-1.0));
  CHECK(SU2_TYPE::GetValue(*limitsEnth.second) == Approx(1.0));

  auto limitsProgvar = look_up_table.GetTableLimitsX();
  CHECK(SU2_TYPE::GetValue(*limitsProgvar.first) == Approx(0.0));
  CHECK(SU2_TYPE::GetValue(*limitsProgvar.second) == Approx(1.0));

  /* lookup value outside of lookup table */

  prog = 1.10;
  enth = 1.1;
  mfrac = 2.0;
  look_up_tag = "Density";
  look_up_table.LookUp_XYZ(look_up_tag, &look_up_dat, prog, enth, mfrac);
  CHECK(look_up_dat == Approx(1.1738796125));
}
