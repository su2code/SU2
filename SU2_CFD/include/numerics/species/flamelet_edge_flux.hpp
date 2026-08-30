/*!
 * \file flamelet_edge_flux.hpp
 * \brief Flamelet transport as a third-layer scalar flux, see numerics/scalar/scalar_edge_flux.hpp.
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

#include "species_edge_flux.hpp"

/*!
 * \class CScalarFlux_Flamelet
 * \ingroup ViscDiscr
 * \brief Convection and diffusion of the flamelet controlling variables and passive species,
 *        which is species transport plus two preferential diffusion terms.
 * \note The preferential diffusion terms read the beta scalars and their gradients from the
 *       auxiliary variables of the solver's own containers, which the per-marker ghost containers
 *       of a boundary do not carry. They are an interior edge term: boundaries instantiate
 *       CScalarFlux_Species, as they did before this model existed.
 */
template <class Double, class FlowIndices, int nDim, size_t nVar = Dynamic>
class CScalarFlux_Flamelet final
    : public CScalarFluxSpeciesBase<Double, CScalarFlux_Flamelet<Double, FlowIndices, nDim, nVar>, FlowIndices, nDim,
                                    nVar> {
 public:
  using Base = CScalarFluxSpeciesBase<Double, CScalarFlux_Flamelet, FlowIndices, nDim, nVar>;
  using Int = typename Base::Int;

  explicit CScalarFlux_Flamelet(const CConfig& config)
      : Base(config),
        preferentialDiffusion(config.GetFlameletParsedOptions().preferential_diffusion),
        nControlVars(config.GetFlameletParsedOptions().n_control_vars) {}

  /*!
   * \brief Preferential diffusion, two terms with the shape of the ordinary diffusion but of
   *        states the model synthesises: div(D grad(beta - phi)) for each controlling variable,
   *        and a thermal term div(beta_T D grad(T)) on the enthalpy equation.
   * \note The thermal term has no implicit part, matching the treatment of the heat flux it
   *       models; the first term has the same thin shear layer Jacobian as the ordinary diffusion,
   *       because it is the same operator applied to a shifted state.
   */
  template <class VariableType>
  FORCEINLINE void extraDiffusionTerms(const FlowIndices& idx, const ScalarFluxOptions& opt, Int iPoint,
                                       const EdgeSide<VariableType>& side_i, Int jPoint,
                                       const EdgeSide<VariableType>& side_j, const CPair<Double>& rho,
                                       const Vector<Double, nDim>& normal, const Vector<Double, nDim>& vector_ij,
                                       EdgeResidual<Double, nVar>& res) const {
    if (!preferentialDiffusion) return;

    const Double dist2_ij = fmax(squaredNorm(vector_ij), EPS);
    const Double proj_vector_ij = dot(vector_ij, normal) / dist2_ij;
    const Double proj_on_rho_i = proj_vector_ij / rho.i;
    const Double proj_on_rho_j = proj_vector_ij / rho.j;

    const Double diffTurb = Base::turbulentDiffusivity(idx, iPoint, side_i, jPoint, side_j);

    /*--- The gradient of a controlling variable is subtracted from that of its beta scalar, so
     * that what is added here is the difference from the ordinary diffusion already applied. ---*/
    for (auto iScalar = 0u; iScalar < nControlVars; ++iScalar) {
      const auto iBeta = betaIndex(iScalar);

      const Double phi_i = gatherVariables(iPoint, side_i.scalarNodes.GetAuxVar(), iBeta) -
                           gatherVariables(iPoint, side_i.scalarNodes.GetSolution(), iScalar);
      const Double phi_j = gatherVariables(jPoint, side_j.scalarNodes.GetAuxVar(), iBeta) -
                           gatherVariables(jPoint, side_j.scalarNodes.GetSolution(), iScalar);

      auto grad_i = gatherVariables<nDim>(iPoint, side_i.scalarNodes.GetAuxVarGradient(), iBeta);
      auto grad_j = gatherVariables<nDim>(jPoint, side_j.scalarNodes.GetAuxVarGradient(), iBeta);
      const auto gradPhi_i = gatherVariables<nDim>(iPoint, side_i.scalarNodes.GetGradient(), iScalar);
      const auto gradPhi_j = gatherVariables<nDim>(jPoint, side_j.scalarNodes.GetGradient(), iScalar);
      for (int iDim = 0; iDim < nDim; ++iDim) {
        grad_i(iDim) -= gradPhi_i(iDim);
        grad_j(iDim) -= gradPhi_j(iDim);
      }

      const Double D_i = gatherVariables(iPoint, side_i.scalarNodes.GetDiffusivity(), iScalar);
      const Double D_j = gatherVariables(jPoint, side_j.scalarNodes.GetDiffusivity(), iScalar);
      const Double D = 0.5 * (rho.i * D_i + rho.j * D_j) + diffTurb;

      const Double projGrad = projectedGradient(opt, grad_i, grad_j, phi_i, phi_j, normal, vector_ij, dist2_ij);

      res.flux_i(iScalar) -= D * projGrad;
      res.flux_j(iScalar) += D * projGrad;

      if (opt.implicit) {
        res.jac_ii(iScalar, iScalar) += D * proj_on_rho_i;
        res.jac_ij(iScalar, iScalar) -= D * proj_on_rho_j;
        res.jac_ji(iScalar, iScalar) -= D * proj_on_rho_i;
        res.jac_jj(iScalar, iScalar) += D * proj_on_rho_j;
      }
    }

    /*--- Thermal term, on the enthalpy equation alone, driven by the temperature gradient. ---*/
    if (nControlVars <= I_ENTH) return;

    const Double T_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.Temperature());
    const Double T_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.Temperature());

    const auto gradT_i = gatherVariables<nDim>(iPoint, side_i.flowNodes->GetGradient_Primitive(), idx.Temperature());
    const auto gradT_j = gatherVariables<nDim>(jPoint, side_j.flowNodes->GetGradient_Primitive(), idx.Temperature());

    const Double Dth_i = gatherVariables(iPoint, side_i.scalarNodes.GetAuxVar(), I_BETA_ENTH_THERMAL) *
                         gatherVariables(iPoint, side_i.scalarNodes.GetDiffusivity(), I_ENTH);
    const Double Dth_j = gatherVariables(jPoint, side_j.scalarNodes.GetAuxVar(), I_BETA_ENTH_THERMAL) *
                         gatherVariables(jPoint, side_j.scalarNodes.GetDiffusivity(), I_ENTH);
    const Double Dth = 0.5 * (rho.i * Dth_i + rho.j * Dth_j) + diffTurb;

    const Double projGradT = projectedGradient(opt, gradT_i, gradT_j, T_i, T_j, normal, vector_ij, dist2_ij);

    res.flux_i(I_ENTH) -= Dth * projGradT;
    res.flux_j(I_ENTH) += Dth * projGradT;
  }

 private:
  const bool preferentialDiffusion;
  const unsigned short nControlVars;

  /*!
   * \brief Auxiliary variable holding the beta scalar of a controlling variable.
   */
  static FORCEINLINE unsigned short betaIndex(unsigned short iScalar) {
    switch (iScalar) {
      case I_PROGVAR:
        return I_BETA_PROGVAR;
      case I_ENTH:
        return I_BETA_ENTH;
      default:
        return I_BETA_MIXFRAC;
    }
  }

  /*!
   * \brief Average gradient of one synthesised state projected on the normal, corrected for
   *        skewness when asked, which is what the ordinary diffusion does for a transported one.
   */
  FORCEINLINE Double projectedGradient(const ScalarFluxOptions& opt, const Vector<Double, nDim>& grad_i,
                                       const Vector<Double, nDim>& grad_j, const Double& phi_i, const Double& phi_j,
                                       const Vector<Double, nDim>& normal, const Vector<Double, nDim>& vector_ij,
                                       const Double& dist2_ij) const {
    Vector<Double, nDim> avgGrad;
    for (int iDim = 0; iDim < nDim; ++iDim) avgGrad(iDim) = 0.5 * (grad_i(iDim) + grad_j(iDim));

    if (opt.correctGradient) {
      const Double corr = (dot(avgGrad, vector_ij) - phi_j + phi_i) / dist2_ij;
      for (int iDim = 0; iDim < nDim; ++iDim) avgGrad(iDim) -= corr * vector_ij(iDim);
    }
    return dot(avgGrad, normal);
  }
};
