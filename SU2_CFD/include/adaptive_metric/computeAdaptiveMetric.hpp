/*!
 * \file computeAdaptiveMetric.hpp
 * \brief Computation of the adaptive metric tensor for anisotropic mesh adaptation.
 * \note The metric of each sensor is built from the absolute value of its Hessian (accumulated in
 *       time for unsteady flows), then optionally normalized in the Lp norm. The absolute value and
 *       the eigenvalue floor are only applied in eigenspace, inside each pass that needs them, and
 *       the final sizes are set by the size and aspect ratio constraints.
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

#pragma once

#include <string>
#include <vector>
#include "../../../Common/include/containers/container_decorators.hpp"

class CGeometry;
class CConfig;

/*!
 * \brief Set all components of the adaptive metric field to zero.
 * \param[in,out] metric - Adaptive metric field (nPoint, nSensor, nSymMat).
 */
void resetAdaptiveMetric(C3DPassiveMatrix& metric);

/*!
 * \brief Add coeff * |H| of each sensor to its adaptive metric, where |H| is the matrix absolute value of the Hessian.
 * \note No eigenvalue floor is applied here; normalization floors the eigenvalues of the accumulated metric.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] hessian - Hessian field of the sensors (nPoint, nSensor, nSymMat).
 * \param[in] coeff - Weight of this contribution (e.g. the time step for unsteady flows, negative to remove one).
 * \param[in,out] metric - Adaptive metric field (nPoint, nSensor, nSymMat).
 */
void accumulateAdaptiveMetric(const CGeometry& geometry, const CVectorOfMatrix& hessian, passivedouble coeff,
                              C3DPassiveMatrix& metric);

/*!
 * \brief Integrate det(|M|)^(p/(2p+n)) of each sensor's metric over the domain, for the Lp-norm normalization.
 * \note Nothing is stored; the eigenvalues are floored in eigenspace only for this computation.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] metric - Adaptive metric field (nPoint, nSensor, nSymMat).
 * \return Global integral for each sensor.
 */
std::vector<passivedouble> integrateAdaptiveMetric(const CGeometry& geometry, const CConfig& config,
                                                   const C3DPassiveMatrix& metric);

/*!
 * \brief Perform the Lp-norm normalization of each sensor's metric and apply the size and aspect ratio constraints.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] integrals - Global integrals from integrateAdaptiveMetric.
 * \param[in] sensorNames - Name of each sensor, for error messages.
 * \param[in,out] metric - Adaptive metric field (nPoint, nSensor, nSymMat).
 */
void normalizeAdaptiveMetric(const CGeometry& geometry, const CConfig& config,
                             const std::vector<passivedouble>& integrals,
                             const std::vector<std::string>& sensorNames, C3DPassiveMatrix& metric);
