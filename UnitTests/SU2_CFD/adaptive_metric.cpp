/*!
 * \file adaptive_metric.cpp
 * \brief Unit tests for the adaptive metric normalization (mesh adaptation).
 * \author B. Munguía
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
#include <memory>
#include <sstream>

#include "../../Common/include/CConfig.hpp"
#include "../../Common/include/containers/container_decorators.hpp"
#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../SU2_CFD/include/adaptive_metric/computeAdaptiveMetric.hpp"

namespace {

/*--- Normalization parameters used in the checks. They must match the ADAPTIVE_METRIC_* options in
 *    AdaptiveMetricTestCase::configOptions. The mesh is a unit cube, so the domain volume is 1. ---*/
constexpr unsigned long COMPLEXITY = 1000;
constexpr passivedouble HMIN = 1e-8, HMAX = 100.0, ARMAX = 1e6;

/*!
 * \brief Unit cube box mesh with three sensors, each given a constant Hessian:
 *        0: identity, 1: rank-one v v^T with v = (1,1,0), 2: zero.
 */
struct AdaptiveMetricTestCase {
  /*--- The ADAPTIVE_METRIC_* values must match the constants above. ---*/
  const std::string configOptions =
      "SOLVER= EULER\n"
      "MESH_FORMAT= BOX\n"
      "INIT_OPTION= TD_CONDITIONS\n"
      "MARKER_EULER= (y_minus, y_plus)\n"
      "MARKER_FAR= (x_minus, x_plus, z_plus, z_minus)\n"
      "MESH_BOX_SIZE= 5,5,5\n"
      "MESH_BOX_LENGTH= 1,1,1\n"
      "MESH_BOX_OFFSET= 0,0,0\n"
      "COMPUTE_ADAPTIVE_METRIC= YES\n"
      "ADAPTIVE_METRIC_SENSOR= (IDENTITY, RANK_ONE, ZERO)\n"
      "ADAPTIVE_METRIC_NORM= 2\n"
      "ADAPTIVE_METRIC_COMPLEXITY= 1000\n"
      "ADAPTIVE_METRIC_HMIN= 1E-8\n"
      "ADAPTIVE_METRIC_HMAX= 100.0\n"
      "ADAPTIVE_METRIC_ARMAX= 1E6\n";

  const std::vector<std::string> sensorNames = {"IDENTITY", "RANK_ONE", "ZERO"};

  std::unique_ptr<CConfig> config;
  std::unique_ptr<CGeometry> geometry;

  AdaptiveMetricTestCase() {
    auto origBuf = cout.rdbuf();
    cout.rdbuf(nullptr);
    stringstream ss(configOptions);
    config = std::unique_ptr<CConfig>(new CConfig(ss, SU2_COMPONENT::SU2_CFD, false));
    {
      auto aux_geometry = std::unique_ptr<CGeometry>(new CPhysicalGeometry(config.get(), 0, 1));
      geometry = std::unique_ptr<CGeometry>(new CPhysicalGeometry(aux_geometry.get(), config.get()));
    }
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
    cout.rdbuf(origBuf);
  }
};

/*!
 * \brief Load the lower-triangular (xx, xy, yy, xz, yz, zz) metric of a sensor at a point.
 */
void loadMetric(const C3DPassiveMatrix& metric, unsigned long iPoint, unsigned long iSensor,
                passivedouble (&mat)[3][3]) {
  for (auto i = 0ul; i < 3; ++i)
    for (auto j = 0ul; j <= i; ++j) mat[i][j] = mat[j][i] = metric(iPoint, iSensor, i * (i + 1) / 2 + j);
}

}  // namespace

TEST_CASE("Adaptive metric normalization", "[AdaptiveMetric]") {
  AdaptiveMetricTestCase test;
  const auto& geometry = *test.geometry;
  const auto& config = *test.config;
  const auto nPoint = geometry.GetnPoint();
  const auto nPointDomain = geometry.GetnPointDomain();
  constexpr unsigned long nSensor = 3, nSymMat = 6;

  /*--- Constant Hessians: identity, rank-one v v^T with v = (1,1,0), and zero. ---*/
  const passivedouble identity[nSymMat] = {1, 0, 1, 0, 0, 1};
  const passivedouble rankOne[nSymMat] = {1, 1, 1, 0, 0, 0};
  CVectorOfMatrix hessian(nPoint, nSensor, nSymMat, 0.0);
  for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
    for (auto iMat = 0ul; iMat < nSymMat; ++iMat) {
      hessian(iPoint, 0, iMat) = identity[iMat];
      hessian(iPoint, 1, iMat) = rankOne[iMat];
    }
  }

  C3DPassiveMatrix metric(nPoint, nSensor, nSymMat, 0.0);
  accumulateAdaptiveMetric(geometry, hessian, 1.0, metric);
  const auto integrals = integrateAdaptiveMetric(geometry, config, metric);
  normalizeAdaptiveMetric(geometry, config, integrals, test.sensorNames, metric);

  /*--- Identity: |H| = I and the integral is the domain volume (1), so M = N^(2/3) I. ---*/
  CHECK(integrals[0] == Approx(1.0).epsilon(1e-12));
  const passivedouble uniform = std::pow(static_cast<passivedouble>(COMPLEXITY), 2.0 / 3.0);

  const passivedouble eigMin = 1.0 / (HMAX * HMAX), eigMax = 1.0 / (HMIN * HMIN);

  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
    passivedouble A[3][3], eigVec[3][3], eigVal[3], work[3];

    /*--- Identity and zero Hessians give the same uniform isotropic metric: for a zero Hessian the
     *    eigenvalue floor cancels out of the normalization. ---*/
    for (const auto iSensor : {0ul, 2ul}) {
      loadMetric(metric, iPoint, iSensor, A);
      for (auto i = 0ul; i < 3; ++i) {
        for (auto j = 0ul; j < 3; ++j) {
          if (i == j) {
            CHECK(A[i][j] == Approx(uniform).epsilon(1e-8));
          } else {
            CHECK(std::abs(A[i][j]) < 1e-8 * uniform);
          }
        }
      }
    }

    /*--- Rank-one: finite, SPD, within the size and aspect ratio constraints, and stretched along v. ---*/
    loadMetric(metric, iPoint, 1, A);
    CBlasStructure::EigenDecomposition(A, eigVec, eigVal, 3, work);

    passivedouble lamMin = eigVal[0], lamMax = eigVal[0];
    unsigned long iMax = 0;
    for (auto iDim = 0ul; iDim < 3; ++iDim) {
      CHECK(std::isfinite(eigVal[iDim]));
      CHECK(eigVal[iDim] >= eigMin * (1 - 1e-10));
      CHECK(eigVal[iDim] <= eigMax * (1 + 1e-10));
      lamMin = std::min(lamMin, eigVal[iDim]);
      if (eigVal[iDim] > lamMax) {
        lamMax = eigVal[iDim];
        iMax = iDim;
      }
    }
    CHECK(lamMin > 0.0);
    CHECK(lamMax / lamMin <= ARMAX * ARMAX * (1 + 1e-8));

    /*--- The largest eigenvalue's eigenvector is along v = (1,1,0)/sqrt(2). ---*/
    const passivedouble alignment = std::abs(eigVec[0][iMax] + eigVec[1][iMax]) / std::sqrt(2.0);
    CHECK(alignment == Approx(1.0).epsilon(1e-10));
  }
}
