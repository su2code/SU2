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
    if (nDim == 2)
      mesh << "5 4 " << i << " " << j << " " << i << "\n";
    else
      mesh << "13 4 " << i << " " << j << " 9 " << i + 5 << " " << j + 5 << " " << i << "\n";
  }
  mesh << "NPOIN= " << (nDim == 2 ? 5 : 10) << "\n";
  for (int layer = 0; layer < (nDim == 2 ? 1 : 2); ++layer)
    for (int i = 0; i < 5; ++i) {
      if (nDim == 2)
        mesh << points[i][0] << " " << points[i][1] << " " << i << "\n";
      else
        mesh << points[i][0] << " " << layer << " " << points[i][1] << " " << i + 5 * layer << "\n";
    }
  mesh << "NMARK= " << (split ? 2 : 1) + (nDim == 3 ? 1 : 0) << "\n";
  const auto wall = [&](const std::string& name, int first, int last) {
    mesh << "MARKER_TAG= " << name << "\nMARKER_ELEMS= " << last - first << "\n";
    for (int i = first; i < last; ++i) {
      const int j = (i + 1) % 4;
      if (nDim == 2)
        mesh << "3 " << i << " " << j << "\n";
      else
        mesh << "9 " << i << " " << j << " " << j + 5 << " " << i + 5 << "\n";
    }
  };
  if (split) {
    wall("upper", 0, 2);
    wall("lower", 2, 4);
  } else
    wall("wall", 0, 4);
  if (nDim == 3) {
    mesh << "MARKER_TAG= ends\nMARKER_ELEMS= 8\n";
    for (int i = 0; i < 4; ++i) {
      const int j = (i + 1) % 4;
      mesh << "5 4 " << j << " " << i << "\n";
      mesh << "5 9 " << i + 5 << " " << j + 5 << "\n";
    }
  }
}

void WriteRefinedDiamond(const std::string& filename, unsigned short nDim, unsigned int n) {
  const std::array<std::array<double, 2>, 4> corners = {{{0, 0}, {1, 1}, {2, 0}, {1, -1}}};
  const unsigned int count = 4 * n, layerSize = count + 1;
  std::ofstream mesh(filename);
  mesh << "NDIME= " << nDim << "\nNELEM= " << count << "\n";
  for (unsigned int i = 0; i < count; ++i) {
    const auto j = (i + 1) % count;
    if (nDim == 2)
      mesh << "5 " << count << " " << i << " " << j << " " << i << "\n";
    else
      mesh << "13 " << count << " " << i << " " << j << " " << count + layerSize << " " << i + layerSize << " "
           << j + layerSize << " " << i << "\n";
  }
  mesh << "NPOIN= " << layerSize * (nDim == 2 ? 1 : 2) << "\n";
  for (unsigned int layer = 0; layer < (nDim == 2 ? 1u : 2u); ++layer) {
    for (unsigned int edge = 0; edge < 4; ++edge) {
      for (unsigned int j = 0; j < n; ++j) {
        const double t = double(j) / n;
        const auto x = (1 - t) * corners[edge][0] + t * corners[(edge + 1) % 4][0];
        const auto z = (1 - t) * corners[edge][1] + t * corners[(edge + 1) % 4][1];
        mesh << x << " ";
        if (nDim == 3) mesh << layer << " ";
        mesh << z << " " << layer * layerSize + edge * n + j << "\n";
      }
    }
    mesh << "1 ";
    if (nDim == 3) mesh << layer << " ";
    mesh << "0 " << layer * layerSize + count << "\n";
  }
  mesh << "NMARK= " << (nDim == 3 ? 2 : 1) << "\nMARKER_TAG= wall\nMARKER_ELEMS= " << count << "\n";
  for (unsigned int i = 0; i < count; ++i) {
    const auto j = (i + 1) % count;
    if (nDim == 2)
      mesh << "3 " << i << " " << j << "\n";
    else
      mesh << "9 " << i << " " << j << " " << j + layerSize << " " << i + layerSize << "\n";
  }
  if (nDim == 3) {
    mesh << "MARKER_TAG= ends\nMARKER_ELEMS= " << 2 * count << "\n";
    for (unsigned int i = 0; i < count; ++i) {
      const auto j = (i + 1) % count;
      mesh << "5 " << count << " " << j << " " << i << "\n";
      mesh << "5 " << count + layerSize << " " << i + layerSize << " " << j + layerSize << "\n";
    }
  }
}

double ReferenceArea(unsigned short nDim, bool split, unsigned int refinements = 0, double configuredArea = 0.0) {
  const std::string file = "reference_diamond_" + std::to_string(nDim) + (split ? "_split.su2" : ".su2");
  if (SU2_MPI::GetRank() == 0) {
    if (refinements)
      WriteRefinedDiamond(file, nDim, refinements);
    else
      WriteDiamond(file, nDim, split);
  }
  SU2_MPI::Barrier(SU2_MPI::GetComm());
  const auto markers = split ? "upper, lower" : "wall";
  std::stringstream options;
  options << "SOLVER= EULER\nMESH_FORMAT= SU2\nMESH_FILENAME= " << file << "\nMACH_NUMBER= 0.3\n"
          << "REF_AREA= " << configuredArea << "\nMARKER_EULER= ( " << markers << " )\nMARKER_MONITORING= ( " << markers
          << " )\n";
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
  /* Splitting the wall at the sharp corners must not change the projected area. */
  for (unsigned short nDim : {2, 3}) {
    CAPTURE(nDim);
    CHECK(ReferenceArea(nDim, false) == Approx(2.0).margin(1e-12));
    CHECK(ReferenceArea(nDim, true) == Approx(2.0).margin(1e-12));
    CHECK(ReferenceArea(nDim, false, 0, 7.0) == Approx(7.0).margin(1e-12));
  }
}

TEST_CASE("Reference area under diamond edge refinement", "[ReferenceAreaRefinement]") {
  for (unsigned short nDim : {2, 3}) {
    for (unsigned int n : {1, 2, 4, 10}) {
      const auto area = ReferenceArea(nDim, false, n);
      if (SU2_MPI::GetRank() == 0)
        std::cout << "REFINEMENT_RESULT dimension=" << nDim << " n=" << n << " area=" << area << std::endl;
      CAPTURE(nDim, n);
      CHECK(area == Approx(2.0).margin(1e-12));
    }
  }
}

#if defined(CODI_REVERSE_TYPE)
TEST_CASE("Reference area derivative matches the diamond chord", "[ReferenceArea][AD]") {
  if (SU2_MPI::GetSize() != 1) return;
  const std::string file = "reference_diamond_ad.su2";
  WriteDiamond(file, 2, false);
  std::stringstream options;
  options << "SOLVER= EULER\nMESH_FORMAT= SU2\nMESH_FILENAME= " << file
          << "\nMACH_NUMBER= 0.3\nREF_AREA= 0\nMARKER_EULER= ( wall )\nMARKER_MONITORING= ( wall )\n";
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

  su2double* trailingX = nullptr;
  for (unsigned long i = 0; i < geometry->GetnPoint(); ++i) {
    if (geometry->nodes->GetGlobalIndex(i) == 2) trailingX = &geometry->nodes->GetCoord(i)[0];
  }
  REQUIRE(trailingX != nullptr);
  /* Moving the trailing edge changes the projected chord one-for-one. */
  AD::Reset();
  AD::StartRecording();
  AD::RegisterInput(*trailingX);
  geometry->SetBoundControlVolume(config.get(), UPDATE);
  geometry->SetPositive_ZArea(config.get());
  su2double area = config->GetRefArea();
  AD::RegisterOutput(area);
  AD::StopRecording();
  SU2_TYPE::SetDerivative(area, 1.0);
  AD::ComputeAdjoint();
  CHECK(SU2_TYPE::GetValue(area) == Approx(2.0).margin(1e-12));
  CHECK(SU2_TYPE::GetDerivative(*trailingX) == Approx(1.0).margin(1e-12));
  AD::Reset();
}
#endif
