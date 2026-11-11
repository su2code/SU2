/*!
 * \file reference_area_tests.cpp
 * \brief Analytical checks for automatically computed projected reference area.
 * \version 8.6.0 "Harrier"
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
 * Distributed under the GNU Lesser General Public License, version 2.1 or later.
 */
#include "catch.hpp"
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/geometry/CPhysicalGeometry.hpp"

namespace {

/* A diamond of width 2, or the same diamond extruded through a unit span. */
void WriteDiamond(const std::string& filename, unsigned short nDim, bool split) {
  const std::array<std::array<double, 2>, 5> points = {{{0, 0}, {1, 1}, {2, 0}, {1, -1}, {1, 0}}};
  std::ofstream mesh(filename);
  mesh << "NDIME= " << nDim << "\nNELEM= 4\n";
  for (int i = 0; i < 4; ++i) {
    const int j = (i + 1) % 4;
    if (nDim == 2) mesh << "5 4 " << i << " " << j << " " << i << "\n";
    else mesh << "13 4 " << i << " " << j << " 9 " << i + 5 << " " << j + 5 << " " << i << "\n";
  }
  mesh << "NPOIN= " << (nDim == 2 ? 5 : 10) << "\n";
  for (int layer = 0; layer < (nDim == 2 ? 1 : 2); ++layer)
    for (int i = 0; i < 5; ++i) {
      if (nDim == 2) mesh << points[i][0] << " " << points[i][1] << " " << i << "\n";
      else mesh << points[i][0] << " " << layer << " " << points[i][1] << " " << i + 5 * layer << "\n";
    }
  mesh << "NMARK= " << (split ? 2 : 1) + (nDim == 3 ? 1 : 0) << "\n";
  const auto wall = [&](const std::string& name, int first, int last) {
    mesh << "MARKER_TAG= " << name << "\nMARKER_ELEMS= " << last - first << "\n";
    for (int i = first; i < last; ++i) {
      const int j = (i + 1) % 4;
      if (nDim == 2) mesh << "3 " << i << " " << j << "\n";
      else mesh << "9 " << i << " " << j << " " << j + 5 << " " << i + 5 << "\n";
    }
  };
  if (split) { wall("upper", 0, 2); wall("lower", 2, 4); }
  else wall("wall", 0, 4);
  if (nDim == 3) {
    mesh << "MARKER_TAG= ends\nMARKER_ELEMS= 8\n";
    for (int i = 0; i < 4; ++i) {
      const int j = (i + 1) % 4;
      mesh << "5 4 " << j << " " << i << "\n";
      mesh << "5 9 " << i + 5 << " " << j + 5 << "\n";
    }
  }
}

void WriteRefinedDiamond(const std::string& filename, unsigned int n) {
  const std::array<std::array<double, 2>, 4> corners = {{{0, 0}, {1, 1}, {2, 0}, {1, -1}}};
  const unsigned int count = 4 * n;
  std::ofstream mesh(filename);
  mesh << "NDIME= 2\nNELEM= " << count << "\n";
  for (unsigned int i = 0; i < count; ++i)
    mesh << "5 " << count << " " << i << " " << (i + 1) % count << " " << i << "\n";
  mesh << "NPOIN= " << count + 1 << "\n";
  for (unsigned int edge = 0; edge < 4; ++edge) {
    for (unsigned int j = 0; j < n; ++j) {
      const double t = double(j) / n;
      mesh << (1-t)*corners[edge][0] + t*corners[(edge+1)%4][0] << " "
           << (1-t)*corners[edge][1] + t*corners[(edge+1)%4][1] << " " << edge*n+j << "\n";
    }
  }
  mesh << "1 0 " << count << "\nNMARK= 1\nMARKER_TAG= wall\nMARKER_ELEMS= " << count << "\n";
  for (unsigned int i = 0; i < count; ++i) mesh << "3 " << i << " " << (i+1)%count << "\n";
}

double ReferenceArea(unsigned short nDim, bool split, unsigned int refinements = 0) {
  const std::string file = "reference_diamond_" + std::to_string(nDim) + (split ? "_split.su2" : ".su2");
  if (SU2_MPI::GetRank() == 0) {
    if (refinements) WriteRefinedDiamond(file, refinements);
    else WriteDiamond(file, nDim, split);
  }
  SU2_MPI::Barrier(SU2_MPI::GetComm());
  const auto markers = split ? "upper, lower" : "wall";
  std::stringstream options;
  options << "SOLVER= EULER\nMESH_FORMAT= SU2\nMESH_FILENAME= " << file << "\nMACH_NUMBER= 0.3\n"
          << "REF_AREA= 0\nMARKER_EULER= ( " << markers << " )\nMARKER_MONITORING= ( " << markers << " )\n";
  if (nDim == 3) options << "MARKER_FAR= ( ends )\n";
  auto config = std::make_unique<CConfig>(options, SU2_COMPONENT::SU2_CFD, false);
  auto aux = std::make_unique<CPhysicalGeometry>(config.get(), 0, 1);
  aux->SetColorGrid_Parallel(config.get());
  auto geometry = std::make_unique<CPhysicalGeometry>(aux.get(), config.get());
  geometry->SetSendReceive(config.get());
  geometry->SetBoundaries(config.get());
  geometry->SetPoint_Connectivity();
  geometry->SetElement_Connectivity();
  geometry->SetBoundVolume();
  geometry->Check_IntElem_Orientation(config.get());
  geometry->Check_BoundElem_Orientation(config.get());
  geometry->SetEdges();
  geometry->SetVertex(config.get());
  geometry->SetControlVolume(config.get(), ALLOCATE);
  geometry->SetBoundControlVolume(config.get(), ALLOCATE);
  geometry->SetPositive_ZArea(config.get());
  return SU2_TYPE::GetValue(config->GetRefArea());
}

}  // namespace

TEST_CASE("Automatic reference area matches a diamond projection", "[ReferenceArea]") {
  for (unsigned short nDim : {2, 3}) {
    CAPTURE(nDim);
    CHECK(ReferenceArea(nDim, false) == Approx(2.0).margin(1e-12));
    CHECK(ReferenceArea(nDim, true) == Approx(2.0).margin(1e-12));
  }
}

TEST_CASE("Reference area under diamond edge refinement", "[ReferenceAreaRefinement]") {
  for (unsigned int n : {1, 2, 4, 10}) {
    const auto area = ReferenceArea(2, false, n);
    if (SU2_MPI::GetRank() == 0) std::cout << "REFINEMENT_RESULT n=" << n << " area=" << area << std::endl;
    CAPTURE(n);
    CHECK(area == Approx(2.0).margin(1e-12));
  }
}
