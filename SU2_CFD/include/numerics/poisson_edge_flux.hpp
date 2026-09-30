/*!
 * \file poisson_edge_flux.hpp
 * \brief Pressure correction (Poisson) equation as a third-layer scalar flux,
 *        see numerics/scalar/scalar_edge_flux.hpp.
 * \author T. Aalbers, P. Gomes
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

#include "scalar/scalar_edge_flux.hpp"

/*!
 * \class CScalarFlux_Poisson
 * \ingroup ViscDiscr
 * \brief Diffusion of the pressure correction, with a diagonal (single-equation) coefficient.
 * \note The equation has no convective term at all, so the solver runs this kernel with
 *       ScalarFluxOptions::convective cleared and the inherited CUpwScalarFlux::finalizeFlux
 *       is never reached; see CPoissonSolver::Viscous_Residual.
 * \note The diffusion coefficient is the momentum coefficient vol/A_p carried by the solver's
 *       own variables, not a flow primitive, so nothing here reads the flow's primitive row;
 *       the flow variables are only asked which points carry a strong velocity BC.
 */
template <class Double, class FlowIndices, int nDim, size_t nVar = 1>
class CScalarFlux_Poisson final
    : public CUpwScalarBase<Double, CScalarFlux_Poisson<Double, FlowIndices, nDim, nVar>, FlowIndices, nDim, nVar> {
 public:
  static constexpr bool Conservative = false;
  static constexpr bool DiagonalDiffusion = true;

  using Base = CUpwScalarBase<Double, CScalarFlux_Poisson, FlowIndices, nDim, nVar>;
  using Int = typename Base::Int;

  using Base::Base;

  /*!
   * \brief Momentum coefficient, an i/j average, identical for both edge sides (TSL Jacobian).
   * \note A point under a strong velocity BC has no momentum equation, and so no momentum
   *       coefficient of its own; the edge uses that of its other node instead.
   */
  template <class VariableType>
  FORCEINLINE CPair<Vector<Double, nVar>> coefficients(const FlowIndices&, Int iPoint,
                                                       const EdgeSide<VariableType>& side_i, Int jPoint,
                                                       const EdgeSide<VariableType>& side_j,
                                                       const CPair<Double>&) const {
    const bool strong_i = side_i.flowNodes->GetStrongBC(iPoint);
    const bool strong_j = side_j.flowNodes->GetStrongBC(jPoint);

    const Double coeff_i = strong_i ? gatherVariables(jPoint, side_j.scalarNodes.GetMomCoeff())
                                    : gatherVariables(iPoint, side_i.scalarNodes.GetMomCoeff());
    const Double coeff_j = strong_j ? gatherVariables(iPoint, side_i.scalarNodes.GetMomCoeff())
                                    : gatherVariables(jPoint, side_j.scalarNodes.GetMomCoeff());

    Vector<Double, nVar> D;
    D(0) = 0.5 * (coeff_i + coeff_j);
    return {D, D};
  }
};
