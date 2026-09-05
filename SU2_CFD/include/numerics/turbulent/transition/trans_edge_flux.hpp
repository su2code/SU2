/*!
 * \file trans_edge_flux.hpp
 * \brief Langtry-Menter transition model as a third-layer scalar flux, see numerics/scalar/scalar_edge_flux.hpp.
 * \author P. Gomes
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

#include "../../scalar/scalar_edge_flux.hpp"

/*!
 * \class CScalarFlux_TransLM
 * \ingroup ViscDiscr
 * \brief Convection and diffusion of the Langtry-Menter transition model, conservative with a
 *        diagonal, i/j-symmetric diffusion matrix. The coefficients depend only on the flow's
 *        mu/mu_t, not on the transported gamma/Re_theta, so no coefficientJacobians override
 *        is needed.
 * \note LM writes no finalizeFlux of its own: the inherited CUpwScalarFlux one is exactly
 *       flux(iVar) = a0*rho_i*phi_i(iVar) + a1*rho_j*phi_j(iVar), Conservative weighting by
 *       density, which is the model's whole convective term.
 */
template <class Double, class FlowIndices, int nDim, size_t nVar = 2>
class CScalarFlux_TransLM
    : public CUpwScalarBase<Double, CScalarFlux_TransLM<Double, FlowIndices, nDim, nVar>, FlowIndices, nDim, nVar> {
 public:
  static constexpr bool Conservative = true;
  static constexpr bool DiagonalDiffusion = true;

  using Base = CUpwScalarBase<Double, CScalarFlux_TransLM, FlowIndices, nDim, nVar>;
  using Int = typename Base::Int;
  using Base::Base;

  /*!
   * \brief Diffusion coefficients, an i/j average of (mu+mu_t) for intermittency and of
   *        2*(mu+mu_t) for the momentum-thickness Reynolds number; identical for both edge sides.
   */
  template <class VariableType>
  FORCEINLINE CPair<Vector<Double, nVar>> coefficients(const FlowIndices& idx, Int iPoint,
                                                       const EdgeSide<VariableType>& side_i, Int jPoint,
                                                       const EdgeSide<VariableType>& side_j,
                                                       const CPair<Double>&) const {
    const Double mu_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.LaminarViscosity());
    const Double mu_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.LaminarViscosity());
    const Double muT_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.EddyViscosity());
    const Double muT_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.EddyViscosity());

    const Double diff_i_gamma = mu_i + muT_i;
    const Double diff_j_gamma = mu_j + muT_j;
    const Double diff_i_ReThetaT = 2.0 * (mu_i + muT_i);
    const Double diff_j_ReThetaT = 2.0 * (mu_j + muT_j);

    Vector<Double, nVar> D;
    D(0) = 0.5 * (diff_i_gamma + diff_j_gamma);
    D(1) = 0.5 * (diff_i_ReThetaT + diff_j_ReThetaT);
    return {D, D};
  }
};
