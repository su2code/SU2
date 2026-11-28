/*!
 * \file adaptive_metric.cpp
 * \brief Unit tests for the Green-Gauss Hessians and the adaptive metric normalization (mesh adaptation).
 * \author B. Munguía
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

#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

#include "../../Common/include/CConfig.hpp"
#include "../../Common/include/containers/container_decorators.hpp"
#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../SU2_CFD/include/solvers/CSolver.hpp"
#include "../../SU2_CFD/include/gradients/computeGradientsGreenGauss.hpp"
#include "../../SU2_CFD/include/adaptive_metric/computeAdaptiveMetric.hpp"

namespace {

/*--- Normalization parameters used in the checks. They must match the ADAPTIVE_METRIC_* options in
 *    configOptions. The mesh is the unit square (2D) or unit cube (3D), so the domain volume is 1.
 *    ARMAX bounds the eigenvalue ratio of the stored metric by ARMAX^2, and re-decomposing the stored
 *    metric recovers its smallest eigenvalue only to ~ARMAX^2 * machine epsilon; ARMAX = 1e3 keeps that
 *    roundoff (~1e-10) well inside the tolerance of the aspect ratio check. ---*/
constexpr unsigned long COMPLEXITY = 1000;
constexpr passivedouble HMIN = 1e-8, HMAX = 100.0, ARMAX = 1e3;

/*!
 * \brief Config of a uniform 9 x 9 (x 9) node mesh of the unit square (cube), with three sensors.
 */
std::string configOptions(unsigned short nDim) {
  const std::string mesh = (nDim == 2) ? "MESH_FORMAT= RECTANGLE\n"
                                         "MARKER_FAR= (x_minus, x_plus, y_minus, y_plus)\n"
                                       : "MESH_FORMAT= BOX\n"
                                         "MARKER_FAR= (x_minus, x_plus, y_minus, y_plus, z_minus, z_plus)\n";
  /*--- The ADAPTIVE_METRIC_* values must match the constants above. ---*/
  return "SOLVER= EULER\n"
         "INIT_OPTION= TD_CONDITIONS\n" +
         mesh +
         "MESH_BOX_SIZE= 9,9,9\n"
         "MESH_BOX_LENGTH= 1,1,1\n"
         "MESH_BOX_OFFSET= 0,0,0\n"
         "COMPUTE_ADAPTIVE_METRIC= YES\n"
         "ADAPTIVE_METRIC_SENSOR= (IDENTITY, RANK_ONE, ZERO)\n"
         "ADAPTIVE_METRIC_NORM= 2\n"
         "ADAPTIVE_METRIC_COMPLEXITY= 1000\n"
         "ADAPTIVE_METRIC_HMIN= 1E-8\n"
         "ADAPTIVE_METRIC_HMAX= 100.0\n"
         "ADAPTIVE_METRIC_ARMAX= 1E3\n";
}

/*!
 * \brief Config and geometry of the test mesh.
 */
struct AdaptiveMetricTestCase {
  const std::vector<std::string> sensorNames = {"IDENTITY", "RANK_ONE", "ZERO"};

  std::unique_ptr<CConfig> config;
  std::unique_ptr<CGeometry> geometry;

  explicit AdaptiveMetricTestCase(unsigned short nDim) {
    auto origBuf = cout.rdbuf();
    cout.rdbuf(nullptr);
    stringstream ss(configOptions(nDim));
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
 * \brief Index of component (i,j), i >= j, in the lower-triangular storage xx, xy, yy[, xz, yz, zz].
 */
constexpr unsigned long symIndex(unsigned long i, unsigned long j) { return i * (i + 1) / 2 + j; }

/*!
 * \brief Load the symmetric tensor of a sensor at a point into a full matrix.
 */
template <size_t nDim, class Container>
void loadSym(const Container& field, unsigned long iPoint, unsigned long iSensor, passivedouble (&mat)[nDim][nDim]) {
  for (auto i = 0ul; i < nDim; ++i)
    for (auto j = 0ul; j <= i; ++j) mat[i][j] = mat[j][i] = field(iPoint, iSensor, symIndex(i, j));
}

/*!
 * \brief Quadratic field f = 1/2 x^T A x + b.x + 1, whose Hessian is A.
 */
struct QuadraticField {
  static constexpr passivedouble A[3][3] = {{2.0, 0.5, -0.3}, {0.5, -1.0, 0.7}, {-0.3, 0.7, 1.5}};
  static constexpr passivedouble b[3] = {0.3, -0.2, 0.1};

  const CGeometry& geometry;

  su2double operator()(unsigned long iPoint, unsigned long) const {
    const auto nDim = geometry.GetnDim();
    const auto x = geometry.nodes->GetCoord(iPoint);
    su2double f = 1.0;
    for (auto i = 0ul; i < nDim; ++i) {
      f += b[i] * x[i];
      for (auto j = 0ul; j < nDim; ++j) f += 0.5 * x[i] * A[i][j] * x[j];
    }
    return f;
  }
};

/*!
 * \brief The Green-Gauss gradient of a quadratic field is exact at the interior points of a uniform Cartesian
 *        mesh, so its Green-Gauss Hessian is exact at the points at least two layers away from the boundary.
 */
template <size_t nDim>
void testGreenGaussHessian() {
  AdaptiveMetricTestCase test(nDim);
  auto& geometry = *test.geometry;
  const auto& config = *test.config;
  const auto nPoint = geometry.GetnPoint();
  constexpr unsigned long nSymMat = 3 * (nDim - 1);

  const QuadraticField field{geometry};
  C3DDoubleMatrix gradient(nPoint, 1, nDim), hessian(nPoint, 1, nSymMat);
  computeGradientsGreenGauss(nullptr, MPI_QUANTITIES::SOLUTION, PERIODIC_NONE, geometry, config, field, 0, 1, -1,
                             gradient);
  computeHessiansGreenGauss(nullptr, MPI_QUANTITIES::HESSIAN, PERIODIC_NONE, geometry, config, gradient, 0, 1, -1,
                            hessian);

  /*--- Mesh lines in each direction, to find the points at least two layers inside the domain. ---*/
  std::vector<passivedouble> lines[nDim];
  for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint)
    for (auto iDim = 0ul; iDim < nDim; ++iDim) lines[iDim].push_back(geometry.nodes->GetCoord(iPoint, iDim));
  for (auto& line : lines) {
    std::sort(line.begin(), line.end());
    line.erase(
        std::unique(line.begin(), line.end(), [](passivedouble a, passivedouble b) { return std::abs(a - b) < 1e-12; }),
        line.end());
  }

  unsigned long nChecked = 0;
  for (auto iPoint = 0ul; iPoint < geometry.GetnPointDomain(); ++iPoint) {
    bool inner = true;
    for (auto iDim = 0ul; iDim < nDim; ++iDim) {
      const auto x = geometry.nodes->GetCoord(iPoint, iDim);
      const auto& line = lines[iDim];
      inner = inner && (x > line[1] + 1e-12) && (x < line[line.size() - 2] - 1e-12);
    }
    if (!inner) continue;

    passivedouble H[nDim][nDim];
    loadSym<nDim>(hessian, iPoint, 0, H);
    for (auto i = 0ul; i < nDim; ++i)
      for (auto j = 0ul; j < nDim; ++j) CHECK(H[i][j] == Approx(QuadraticField::A[i][j]).margin(1e-10));
    ++nChecked;
  }

  /*--- 9 nodes per direction, minus two layers on each side. ---*/
  CHECK(nChecked == static_cast<unsigned long>(std::pow(5, nDim)));
}

/*!
 * \brief Normalize three sensors with constant Hessians: 0: identity, 1: rank-one v v^T with v = (1,1[,0]), 2: zero.
 */
template <size_t nDim>
void testNormalization() {
  AdaptiveMetricTestCase test(nDim);
  const auto& geometry = *test.geometry;
  const auto& config = *test.config;
  const auto nPoint = geometry.GetnPoint();
  constexpr unsigned long nSensor = 3, nSymMat = 3 * (nDim - 1);

  const passivedouble v[3] = {1.0, 1.0, 0.0};
  CVectorOfMatrix hessian(nPoint, nSensor, nSymMat, 0.0);
  for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
    for (auto i = 0ul; i < nDim; ++i) {
      for (auto j = 0ul; j <= i; ++j) {
        hessian(iPoint, 0, symIndex(i, j)) = (i == j) ? 1.0 : 0.0;
        hessian(iPoint, 1, symIndex(i, j)) = v[i] * v[j];
      }
    }
  }

  C3DPassiveMatrix metric(nPoint, nSensor, nSymMat, 0.0);
  accumulateAdaptiveMetric(geometry, hessian, 1.0, metric);
  const auto integrals = integrateAdaptiveMetric(geometry, config, metric);
  normalizeAdaptiveMetric(geometry, config, integrals, test.sensorNames, metric);

  /*--- Identity: |H| = I and the integral is the domain volume (1), so M = N^(2/nDim) I. ---*/
  CHECK(integrals[0] == Approx(1.0).epsilon(1e-12));
  const passivedouble uniform = std::pow(static_cast<passivedouble>(COMPLEXITY), 2.0 / nDim);

  const passivedouble eigMin = 1.0 / (HMAX * HMAX), eigMax = 1.0 / (HMIN * HMIN);

  for (auto iPoint = 0ul; iPoint < geometry.GetnPointDomain(); ++iPoint) {
    passivedouble A[nDim][nDim], eigVec[nDim][nDim], eigVal[nDim], work[nDim];

    /*--- Identity and zero Hessians give the same uniform isotropic metric: for a zero Hessian the
     *    eigenvalue floor cancels out of the normalization. ---*/
    for (const auto iSensor : {0ul, 2ul}) {
      loadSym<nDim>(metric, iPoint, iSensor, A);
      for (auto i = 0ul; i < nDim; ++i) {
        for (auto j = 0ul; j < nDim; ++j) {
          if (i == j) {
            CHECK(A[i][j] == Approx(uniform).epsilon(1e-8));
          } else {
            CHECK(std::abs(A[i][j]) < 1e-8 * uniform);
          }
        }
      }
    }

    /*--- Rank-one: finite, SPD, within the size and aspect ratio constraints, and stretched along v. ---*/
    loadSym<nDim>(metric, iPoint, 1, A);
    CBlasStructure::EigenDecomposition(A, eigVec, eigVal, nDim, work);

    passivedouble lamMin = eigVal[0], lamMax = eigVal[0];
    unsigned long iMax = 0;
    for (auto iDim = 0ul; iDim < nDim; ++iDim) {
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

    /*--- The largest eigenvalue's eigenvector is along v / |v| = (1,1[,0]) / sqrt(2). ---*/
    const passivedouble alignment = std::abs(eigVec[0][iMax] + eigVec[1][iMax]) / std::sqrt(2.0);
    CHECK(alignment == Approx(1.0).epsilon(1e-10));
  }
}

}  // namespace

TEST_CASE("Green-Gauss Hessian of a quadratic field 2D", "[AdaptiveMetric]") { testGreenGaussHessian<2>(); }

TEST_CASE("Green-Gauss Hessian of a quadratic field 3D", "[AdaptiveMetric]") { testGreenGaussHessian<3>(); }

TEST_CASE("Adaptive metric normalization 2D", "[AdaptiveMetric]") { testNormalization<2>(); }

TEST_CASE("Adaptive metric normalization 3D", "[AdaptiveMetric]") { testNormalization<3>(); }
