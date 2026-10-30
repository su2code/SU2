/*!
 * \file turb_sa_edge_flux.hpp
 * \brief Spalart-Allmaras model as a third-layer scalar flux, see numerics/scalar/scalar_edge_flux.hpp.
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
 * \class CScalarFlux_SA
 * \ingroup ConvDiscr
 * \ingroup ViscDiscr
 * \brief Convection and diffusion of the Spalart-Allmaras model, non-conservative and with a
 *        diagonal (but asymmetric) diffusion coefficient.
 * \note SA writes its own convective term rather than using the inherited CUpwScalarFlux one,
 *       because with stochastic backscatter active (nVar 4) the three Langevin equations are
 *       advected with a centered flux, unlike the plain upwind SA equation itself.
 */
template <class Double, class FlowIndices, int nDim, size_t nVar>
class CScalarFlux_SA
    : public CUpwScalarBase<Double, CScalarFlux_SA<Double, FlowIndices, nDim, nVar>, FlowIndices, nDim, nVar> {
 public:
  static constexpr bool Conservative = false;
  static constexpr bool DiagonalDiffusion = true;
  static constexpr bool DiffusionReadsDensity = true; /*!< \brief The kinematic viscosities below. */

  using Base = CUpwScalarBase<Double, CScalarFlux_SA, FlowIndices, nDim, nVar>;
  using Int = typename Base::Int;

  explicit CScalarFlux_SA(const CConfig& config)
      : Base(config),
        negativeSA(config.GetSAParsedOptions().version == SA_OPTIONS::NEG),
        accurateJacobians(config.GetUse_Accurate_Turb_Jacobians()) {}

 private:
  static constexpr passivedouble sigma = 2.0 / 3.0; /*!< \brief Constant of the diffusion term. */
  static constexpr passivedouble cb2 = 0.622;       /*!< \brief Constant of the diffusion term. */
  static constexpr passivedouble cn1 = 16.0;        /*!< \brief Constant of the SA-neg diffusion correction. */

  /*!< \brief Whether nu_tilde may go negative (SA_OPTIONS= NEGATIVE), which needs the fn-corrected
   *          diffusion coefficient below to keep the diffusion term from turning anti-diffusive. */
  const bool negativeSA;
  const bool accurateJacobians;

 public:
  /*!
   * \brief SA convection, plus the centered advection of the backscatter equations when nVar > 1.
   */
  template <class VariableType, size_t Size>
  FORCEINLINE void finalizeFlux(const FlowIndices&, const ScalarFluxOptions& opt, Int, const EdgeSide<VariableType>&,
                                Int, const EdgeSide<VariableType>&, const Double& a0, const Double& a1,
                                const CPair<Double>&, const CPair<CScalarValues<Double, Size>>& phi,
                                EdgeResidual<Double, nVar>& res) const {
    const Double flux = a0 * phi.i.all(0) + a1 * phi.j.all(0);

    res.flux_i(0) += flux;
    if (!opt.oneSided) res.flux_j(0) -= flux;

    if (opt.implicit) {
      res.jac_ii(0, 0) += a0;
      if (!opt.oneSided) {
        res.jac_ij(0, 0) += a1;
        res.jac_ji(0, 0) -= a0;
        res.jac_jj(0, 0) -= a1;
      }
    }

    /*--- Stochastic backscatter: three Langevin equations, advected with the mean of the two
     * upwinding weights and with no diffusion. ---*/
    const Double avg = 0.5 * (a0 + a1);
    for (size_t iVar = 1; iVar < res.nVar; ++iVar) {
      const Double flux_bs = avg * (phi.i.all(iVar) + phi.j.all(iVar));

      res.flux_i(iVar) += flux_bs;
      if (!opt.oneSided) res.flux_j(iVar) -= flux_bs;

      if (opt.implicit) {
        res.jac_ii(iVar, iVar) += avg;
        if (!opt.oneSided) {
          res.jac_ij(iVar, iVar) += avg;
          res.jac_ji(iVar, iVar) -= avg;
          res.jac_jj(iVar, iVar) -= avg;
        }
      }
    }
  }

  /*!
   * \brief fn, the positivity-preserving correction to the SA-neg diffusion coefficient
   *        (Allmaras, Johnson & Spalart), 1 when nu_tilde is not negative enough to need it.
   */
  FORCEINLINE Double fn(const Double& zeta) const {
    if (!negativeSA || zeta >= 0.0) return 1.0;
    const Double zeta3 = zeta * zeta * zeta;
    return (cn1 + zeta3) / (cn1 - zeta3);
  }

  /*!
   * \brief Diffusion coefficients of both orientations of the edge.
   * \note The coefficient is not symmetric: it uses the transported variable of the row it is
   *       going to be used for (the quadratic, non-conservative part of the diffusion term).
   *       Coefficients past index 0 are left at zero, the backscatter equations have no diffusion.
   */
  template <class VariableType>
  FORCEINLINE CPair<Vector<Double, nVar>> coefficients(const FlowIndices& idx, Int iPoint,
                                                       const EdgeSide<VariableType>& side_i, Int jPoint,
                                                       const EdgeSide<VariableType>& side_j,
                                                       const CPair<Double>& rho) const {
    const Double nu_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.LaminarViscosity()) / rho.i;
    const Double nu_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.LaminarViscosity()) / rho.j;

    const Double nuTilde_i = gatherVariables(iPoint, side_i.scalarNodes.GetSolution(), 0);
    const Double nuTilde_j = gatherVariables(jPoint, side_j.scalarNodes.GetSolution(), 0);

    const Double nu_ij = 0.5 * (nu_i + nu_j);
    const Double nuTilde_ij = 0.5 * (nuTilde_i + nuTilde_j);

    /*--- fn is only ever != 1 under SA_OPTIONS= NEGATIVE, and then only where the row's own
     * nu_tilde pulls the coefficient negative; without it (nu + nu_tilde going anti-diffusive)
     * the equation diverges, see Allmaras, Johnson & Spalart's negative SA modification. ---*/
    const Double fn_i = fn(((1.0 + cb2) * nuTilde_ij - cb2 * nuTilde_i) / nu_ij);
    const Double fn_j = fn(((1.0 + cb2) * nuTilde_ij - cb2 * nuTilde_j) / nu_ij);

    Vector<Double, nVar> D_i, D_j;
    D_i(0) = (nu_ij + (1.0 + cb2) * nuTilde_ij * fn_i - cb2 * nuTilde_i * fn_i) / sigma;
    D_j(0) = (nu_ij + (1.0 + cb2) * nuTilde_ij * fn_j - cb2 * nuTilde_j * fn_j) / sigma;
    for (size_t iVar = 1; iVar < nVar; ++iVar) {
      D_i(iVar) = 0.0;
      D_j(iVar) = 0.0;
    }
    return {D_i, D_j};
  }

  /*!
   * \brief Extra Jacobian terms from the dependence of the diffusion coefficient on nu_tilde.
   * \note Skipped for SA-neg, whose fn-corrected coefficient was found (upstream, pre-migration)
   *       to diverge with exact Jacobians; frozen (TSL-only) Jacobians are used there instead.
   *       Skipped for standard SA too unless USE_ACCURATE_TURB_JACOBIANS is set, matching the
   *       pre-migration default of using frozen Jacobians there as well.
   */
  template <class Coefficients, size_t Size>
  FORCEINLINE void coefficientJacobians(const ScalarFluxOptions& opt, const Coefficients&,
                                        const Vector<Double, Size>& projGrad, EdgeResidual<Double, nVar>& res) const {
    if (negativeSA || !accurateJacobians) return;

    /*--- d(diffusion coefficient of i)/d(nu_tilde_i), and its counterpart w.r.t. nu_tilde_j;
     * the coefficient of j is the same expression with i and j swapped, so the same two
     * derivatives apply to both orientations. Both are per-edge constants, so the coefficients
     * themselves are not read here. ---*/
    const Double dDC_dNuTilde_i = ((1.0 + cb2) * 0.5 - cb2) / sigma;
    const Double dDC_dNuTilde_j = (1.0 + cb2) * 0.5 / sigma;

    res.jac_ii(0, 0) -= dDC_dNuTilde_i * projGrad(0);
    if (opt.oneSided) return;

    res.jac_ij(0, 0) -= dDC_dNuTilde_j * projGrad(0);
    res.jac_ji(0, 0) += dDC_dNuTilde_j * projGrad(0);
    res.jac_jj(0, 0) += dDC_dNuTilde_i * projGrad(0);
  }
};
