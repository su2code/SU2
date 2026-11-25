/*!
 * \file computeAdaptiveMetric.cpp
 * \brief Computation of the adaptive metric tensor for anisotropic mesh adaptation.
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

#include "../../include/adaptive_metric/computeAdaptiveMetric.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../../Common/include/parallelization/mpi_structure.hpp"

namespace {

/*--- Lower bound on |lambda|, applied only to the local eigenvalues within a pass (never stored).
 *    It keeps det(|M|) > 0 so the negative power in the normalization stays finite; the size
 *    constraints then clip any floored eigenvalue to 1/hmax^2. ---*/
constexpr passivedouble EIGENVALUE_FLOOR = 1e-20;

/*!
 * \brief Index of component (i,j), i >= j, in the lower-triangular storage of a symmetric tensor.
 */
constexpr unsigned long symIndex(unsigned long i, unsigned long j) { return i * (i + 1) / 2 + j; }

/*!
 * \brief Load the symmetric tensor of a sensor at a point into a full matrix.
 */
template <size_t nDim, class Container>
void loadSymTensor(const Container& field, unsigned long iPoint, unsigned long iSensor,
                   passivedouble (&mat)[nDim][nDim]) {
  for (auto i = 0ul; i < nDim; ++i) {
    for (auto j = 0ul; j <= i; ++j) {
      passivedouble value;
      if constexpr (std::is_same<typename Container::Scalar, passivedouble>::value) {
        value = field(iPoint, iSensor, symIndex(i, j));
      } else {
        value = SU2_TYPE::GetValue(field(iPoint, iSensor, symIndex(i, j)));
      }
      mat[i][j] = value;
      mat[j][i] = value;
    }
  }
}

/*!
 * \brief Store the lower triangle of a full symmetric matrix as the tensor of a sensor at a point.
 */
template <size_t nDim>
void storeSymTensor(const passivedouble (&mat)[nDim][nDim], unsigned long iPoint, unsigned long iSensor,
                    C3DPassiveMatrix& field) {
  for (auto i = 0ul; i < nDim; ++i)
    for (auto j = 0ul; j <= i; ++j) field(iPoint, iSensor, symIndex(i, j)) = mat[i][j];
}

/*!
 * \brief Eigen-decompose a symmetric matrix and replace its eigenvalues by max(|lambda|, floor).
 * \note Non-finite eigenvalues are replaced by the floor.
 */
template <size_t nDim>
void absEigenDecomposition(const passivedouble (&mat)[nDim][nDim], passivedouble (&eigVec)[nDim][nDim],
                           passivedouble (&eigVal)[nDim], passivedouble floor) {
  passivedouble work[nDim];
  CBlasStructure::EigenDecomposition(mat, eigVec, eigVal, nDim, work);
  for (auto iDim = 0ul; iDim < nDim; ++iDim)
    eigVal[iDim] = std::isfinite(eigVal[iDim]) ? std::max(std::fabs(eigVal[iDim]), floor) : floor;
}

template <size_t nDim>
passivedouble productOf(const passivedouble (&eigVal)[nDim]) {
  passivedouble prod = 1.0;
  for (auto iDim = 0ul; iDim < nDim; ++iDim) prod *= eigVal[iDim];
  return prod;
}

template <size_t nDim>
void accumulate(const CGeometry& geometry, const CVectorOfMatrix& hessian, passivedouble coeff,
                C3DPassiveMatrix& metric) {
  const auto nPointDomain = geometry.GetnPointDomain();
  const auto nSensor = metric.rows();

  passivedouble A[nDim][nDim], eigVec[nDim][nDim], eigVal[nDim];

  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
    for (auto iSensor = 0ul; iSensor < nSensor; ++iSensor) {
      loadSymTensor<nDim>(hessian, iPoint, iSensor, A);
      absEigenDecomposition<nDim>(A, eigVec, eigVal, 0.0);
      CBlasStructure::EigenRecomposition(A, eigVec, eigVal, nDim);

      for (auto i = 0ul; i < nDim; ++i)
        for (auto j = 0ul; j <= i; ++j) metric(iPoint, iSensor, symIndex(i, j)) += coeff * A[i][j];
    }
  }
}

template <size_t nDim>
std::vector<passivedouble> integrate(const CGeometry& geometry, const CConfig& config,
                                     const C3DPassiveMatrix& metric) {
  const auto nPointDomain = geometry.GetnPointDomain();
  const auto nSensor = metric.rows();

  const passivedouble p = config.GetAdaptiveMetricParsedOptions().norm;
  const passivedouble normExp = p / (2.0 * p + nDim);

  passivedouble A[nDim][nDim], eigVec[nDim][nDim], eigVal[nDim];
  std::vector<passivedouble> localIntegral(nSensor, 0.0);

  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
    const passivedouble volume = SU2_TYPE::GetValue(geometry.nodes->GetVolume(iPoint));
    for (auto iSensor = 0ul; iSensor < nSensor; ++iSensor) {
      loadSymTensor<nDim>(metric, iPoint, iSensor, A);
      absEigenDecomposition<nDim>(A, eigVec, eigVal, EIGENVALUE_FLOOR);
      localIntegral[iSensor] += std::pow(productOf<nDim>(eigVal), normExp) * volume;
    }
  }

  /*--- SU2_MPI with MPI_DOUBLE expects su2double buffers, which are AD types in AD and direct differentiation
   *    builds, so reduce through su2double copies of the (passive) integrals. ---*/
  std::vector<su2double> localBuffer(localIntegral.begin(), localIntegral.end()), globalBuffer(nSensor, 0.0);
  SU2_MPI::Allreduce(localBuffer.data(), globalBuffer.data(), static_cast<int>(nSensor), MPI_DOUBLE, MPI_SUM,
                     SU2_MPI::GetComm());

  std::vector<passivedouble> globalIntegral(nSensor);
  for (auto iSensor = 0ul; iSensor < nSensor; ++iSensor) globalIntegral[iSensor] = SU2_TYPE::GetValue(globalBuffer[iSensor]);

  return globalIntegral;
}

template <size_t nDim>
void normalize(const CGeometry& geometry, const CConfig& config, const std::vector<passivedouble>& integrals,
               const std::vector<std::string>& sensorNames, C3DPassiveMatrix& metric) {
  const auto nPointDomain = geometry.GetnPointDomain();
  const auto nSensor = metric.rows();
  const auto& adapt = config.GetAdaptiveMetricParsedOptions();

  /*--- Constants defining the normalization ---*/
  const passivedouble p = adapt.norm;
  const passivedouble N = static_cast<passivedouble>(adapt.complexity) * adapt.n_time_subintervals;
  const passivedouble normExp = -1.0 / (2.0 * p + nDim);

  std::vector<passivedouble> globalFactor(nSensor);
  for (auto iSensor = 0ul; iSensor < nSensor; ++iSensor) {
    if (!(integrals[iSensor] > std::numeric_limits<passivedouble>::min()) || !std::isfinite(integrals[iSensor])) {
      SU2_MPI::Error("Cannot normalize the adaptive metric of sensor " + sensorNames[iSensor] +
                     ": its global integral is zero or not finite (is the sensor field constant?).",
                     CURRENT_FUNCTION);
    }
    globalFactor[iSensor] = std::pow(N / integrals[iSensor], 2.0 / nDim);
  }

  /*--- Size and aspect ratio constraints of each sensor ---*/
  std::vector<passivedouble> eigMin(nSensor), eigMax(nSensor), armax2(nSensor);
  for (auto iSensor = 0ul; iSensor < nSensor; ++iSensor) {
    const passivedouble hmin = SU2_TYPE::GetValue(adapt.hmin[iSensor]);
    const passivedouble hmax = SU2_TYPE::GetValue(adapt.hmax[iSensor]);
    const passivedouble armax = SU2_TYPE::GetValue(adapt.armax[iSensor]);
    eigMax[iSensor] = 1.0 / (hmin * hmin);
    eigMin[iSensor] = 1.0 / (hmax * hmax);
    armax2[iSensor] = armax * armax;
  }

  passivedouble A[nDim][nDim], eigVec[nDim][nDim], eigVal[nDim];

  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
    for (auto iSensor = 0ul; iSensor < nSensor; ++iSensor) {
      loadSymTensor<nDim>(metric, iPoint, iSensor, A);
      absEigenDecomposition<nDim>(A, eigVec, eigVal, EIGENVALUE_FLOOR);

      /*--- Local scaling of the Lp-norm normalization ---*/
      const passivedouble factor = globalFactor[iSensor] * std::pow(productOf<nDim>(eigVal), normExp);

      /*--- Scale and clip by the size constraints ---*/
      for (auto iDim = 0ul; iDim < nDim; ++iDim)
        eigVal[iDim] = std::min(std::max(factor * eigVal[iDim], eigMin[iSensor]), eigMax[iSensor]);

      /*--- Clip by the aspect ratio constraint ---*/
      const passivedouble eigLargest = *std::max_element(eigVal, eigVal + nDim);
      for (auto iDim = 0ul; iDim < nDim; ++iDim) eigVal[iDim] = std::max(eigVal[iDim], eigLargest / armax2[iSensor]);

      CBlasStructure::EigenRecomposition(A, eigVec, eigVal, nDim);
      storeSymTensor<nDim>(A, iPoint, iSensor, metric);
    }
  }
}

}  // namespace

void resetAdaptiveMetric(C3DPassiveMatrix& metric) { std::fill(metric.data(), metric.data() + metric.size(), 0.0); }

void accumulateAdaptiveMetric(const CGeometry& geometry, const CVectorOfMatrix& hessian, passivedouble coeff,
                              C3DPassiveMatrix& metric) {
  switch (geometry.GetnDim()) {
    case 2:
      accumulate<2>(geometry, hessian, coeff, metric);
      break;
    case 3:
      accumulate<3>(geometry, hessian, coeff, metric);
      break;
    default:
      SU2_MPI::Error("The adaptive metric is only available in 2D and 3D.", CURRENT_FUNCTION);
      break;
  }
}

std::vector<passivedouble> integrateAdaptiveMetric(const CGeometry& geometry, const CConfig& config,
                                                   const C3DPassiveMatrix& metric) {
  switch (geometry.GetnDim()) {
    case 2:
      return integrate<2>(geometry, config, metric);
    case 3:
      return integrate<3>(geometry, config, metric);
    default:
      SU2_MPI::Error("The adaptive metric is only available in 2D and 3D.", CURRENT_FUNCTION);
      return {};
  }
}

void normalizeAdaptiveMetric(const CGeometry& geometry, const CConfig& config,
                             const std::vector<passivedouble>& integrals,
                             const std::vector<std::string>& sensorNames, C3DPassiveMatrix& metric) {
  switch (geometry.GetnDim()) {
    case 2:
      normalize<2>(geometry, config, integrals, sensorNames, metric);
      break;
    case 3:
      normalize<3>(geometry, config, integrals, sensorNames, metric);
      break;
    default:
      SU2_MPI::Error("The adaptive metric is only available in 2D and 3D.", CURRENT_FUNCTION);
      break;
  }
}
