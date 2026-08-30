/*!
 * \file species_edge_flux.hpp
 * \brief Species transport model as a third-layer scalar flux, see numerics/scalar/scalar_edge_flux.hpp.
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

#include "../scalar/scalar_edge_flux.hpp"

/*!
 * \class CScalarFlux_Species
 * \ingroup ViscDiscr
 * \brief Convection and diffusion of the species transport model, conservative with a diagonal,
 *        i/j-symmetric diffusion matrix. Unlike SA/SST/LM, the equation count is only known at
 *        runtime (one per transported species), so this is the framework's first Dynamic-nVar
 *        model: nEqn is passed to the base explicitly, and coefficients() loops to it rather than
 *        to a compile-time nVar.
 * \note Species writes no finalizeFlux of its own: the inherited CUpwScalarFlux one is exactly
 *       flux(iVar) = a0*rho_i*Y_i(iVar) + a1*rho_j*Y_j(iVar), Conservative weighting by density,
 *       which is the model's whole convective term.
 */
template <class Double, class FlowIndices, int nDim, size_t nVar = Dynamic>
class CScalarFlux_Species
    : public CUpwScalarBase<Double, CScalarFlux_Species<Double, FlowIndices, nDim, nVar>, FlowIndices, nDim, nVar> {
 public:
  static constexpr bool Conservative = true;
  static constexpr bool DiagonalDiffusion = true;

  using Base = CUpwScalarBase<Double, CScalarFlux_Species, FlowIndices, nDim, nVar>;
  using Int = typename Base::Int;

  explicit CScalarFlux_Species(const CConfig& config)
      : Base(config, config.GetnSpecies()),
        turbulence(config.GetKind_Turb_Model() != TURB_MODEL::NONE),
        Sc_t(config.GetSchmidt_Number_Turbulent()) {}

  /*!
   * \brief Diffusion coefficients, an i/j average of (rho * mass diffusivity) per species, plus a
   *        turbulent (mu_t/Sc_t) contribution shared by every species, when a turbulence model is
   *        active; identical for both edge sides.
   * \note The laminar and turbulent averages are kept as two separate 0.5*(...) terms summed at
   *       the end, matching CAvgGrad_Species::FinishResidualCalc's exact operation order, rather
   *       than folding the turbulent term into the same average as the laminar one.
   */
  template <class VariableType>
  FORCEINLINE CPair<Vector<Double, Base::Size>> coefficients(const FlowIndices& idx, Int iPoint,
                                                             const EdgeSide<VariableType>& side_i, Int jPoint,
                                                             const EdgeSide<VariableType>& side_j) const {
    const Double rho_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.Density());
    const Double rho_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.Density());

    Double diffTurb = 0.0;
    if (turbulence) {
      const Double muT_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.EddyViscosity());
      const Double muT_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.EddyViscosity());
      diffTurb = 0.5 * (muT_i / Sc_t + muT_j / Sc_t);
    }

    Vector<Double, Base::Size> D;
    for (size_t iVar = 0; iVar < this->nEqn; ++iVar) {
      const Double D_lam_i = gatherVariables(iPoint, side_i.scalarNodes.GetDiffusivity(), iVar);
      const Double D_lam_j = gatherVariables(jPoint, side_j.scalarNodes.GetDiffusivity(), iVar);
      const Double diffLam = 0.5 * (rho_i * D_lam_i + rho_j * D_lam_j);
      D(iVar) = diffLam + diffTurb;
    }
    return {D, D};
  }

 private:
  const bool turbulence;
  const su2double Sc_t;
};
