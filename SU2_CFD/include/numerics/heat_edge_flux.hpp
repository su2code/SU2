/*!
 * \file heat_edge_flux.hpp
 * \brief Heat transport as a third-layer scalar flux, see numerics/scalar/scalar_edge_flux.hpp.
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

#include "scalar/scalar_edge_flux.hpp"

/*!
 * \class CScalarFlux_Heat
 * \ingroup ViscDiscr
 * \brief Convection and diffusion of temperature, non-conservative, with a diagonal
 *        (single-equation) diffusion coefficient.
 * \note The temperature has no notion of density weighting, so Conservative is false and the
 *       inherited CUpwScalarFlux::finalizeFlux, flux(0) = a0*phi_i(0) + a1*phi_j(0), is exactly
 *       the model's whole convective term; no override is needed.
 * \note The solver runs in two modes, weakly-coupled energy equation on a fluid zone or standalone
 *       conduction on a solid one (CHeatSolver::flow); the diffusion coefficient is the flow's
 *       thermal conductivity over specific heat, plus a turbulent contribution, in the former, and
 *       the configured constant thermal diffusivity in the latter. EdgeSide::flowNodes is null in
 *       the solid case, so this is the only place that may read it, and only when flow is set.
 */
template <class Double, class FlowIndices, int nDim, size_t nVar = 1>
class CScalarFlux_Heat final
    : public CUpwScalarBase<Double, CScalarFlux_Heat<Double, FlowIndices, nDim, nVar>, FlowIndices, nDim, nVar> {
 public:
  static constexpr bool Conservative = false;
  static constexpr bool DiagonalDiffusion = true;

  using Base = CUpwScalarBase<Double, CScalarFlux_Heat, FlowIndices, nDim, nVar>;
  using Int = typename Base::Int;

  explicit CScalarFlux_Heat(const CConfig& config)
      : Base(config),
        flow(config.GetFluidProblem()),
        prandtlTurb(config.GetPrandtl_Turb()),
        constDiffusivity(config.GetThermalDiffusivity()) {}

  /*!
   * \brief Thermal diffusivity, an i/j average, identical for both edge sides (TSL Jacobian).
   */
  template <class VariableType>
  FORCEINLINE CPair<Vector<Double, nVar>> coefficients(const FlowIndices&, Int iPoint,
                                                       const EdgeSide<VariableType>& side_i, Int jPoint,
                                                       const EdgeSide<VariableType>& side_j,
                                                       const CPair<Double>&) const {
    Vector<Double, nVar> D;
    if (flow) {
      const Double diff_i = side_i.flowNodes->GetThermalConductivity(iPoint) / side_i.flowNodes->GetSpecificHeatCp(iPoint) +
                            side_i.flowNodes->GetEddyViscosity(iPoint) / prandtlTurb;
      const Double diff_j = side_j.flowNodes->GetThermalConductivity(jPoint) / side_j.flowNodes->GetSpecificHeatCp(jPoint) +
                            side_j.flowNodes->GetEddyViscosity(jPoint) / prandtlTurb;
      D(0) = 0.5 * (diff_i + diff_j);
    } else {
      D(0) = constDiffusivity;
    }
    return {D, D};
  }

 private:
  const bool flow;
  const su2double prandtlTurb;
  const su2double constDiffusivity;
};
