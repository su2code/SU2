/*!
 * \file CSurfaceFVMDataSorter_tests.cpp
 * \brief Tests of the surface connectivity sorting of the FVM output.
 * \version 8.6.0 "Harrier"
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

#include <string>
#include <vector>

#include "../../UnitQuadTestCase.hpp"
#include "../../../SU2_CFD/include/output/filewriter/CFVMDataSorter.hpp"
#include "../../../SU2_CFD/include/output/filewriter/CSurfaceFVMDataSorter.hpp"

namespace {

/*--- Number of quadrilaterals in the sorted surface connectivity of the given markers of a box mesh. ---*/
unsigned long CountSurfaceQuads(const std::string& boxSize, const std::vector<std::string>& markers) {
  UnitQuadTestCase test;
  test.SetOption("MESH_BOX_SIZE=" + boxSize);
  test.InitConfig();
  test.InitGeometry();

  CFVMDataSorter volumeSorter(test.config.get(), test.geometry.get(), {"x", "y", "z"});
  CSurfaceFVMDataSorter surfaceSorter(test.config.get(), test.geometry.get(), &volumeSorter);
  surfaceSorter.SortConnectivity(test.config.get(), test.geometry.get(), markers);

  return surfaceSorter.GetnElemGlobal(QUADRILATERAL);
}

const std::vector<std::string> allMarkers = {"x_minus", "x_plus", "y_minus", "y_plus", "z_minus", "z_plus"};

}  // namespace

TEST_CASE("Surface sorting keeps the elements of every marker", "[Output]") {
  /*--- One element per marker: the element index is 0 on every marker. ---*/
  CHECK(CountSurfaceQuads("2,2,2", allMarkers) == 6);
  CHECK(CountSurfaceQuads("2,2,2", {"x_minus", "z_plus"}) == 2);
  CHECK(CountSurfaceQuads("2,2,2", {"y_plus"}) == 1);

  /*--- 4x4 elements per marker. ---*/
  CHECK(CountSurfaceQuads("5,5,5", allMarkers) == 6 * 16);
}
