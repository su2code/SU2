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
#include "../../variables/CSpeciesFlameletVariable.hpp"

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
        nControlVars(config.GetFlameletParsedOptions().n_control_vars),
        pdMethod(config.GetFlameletParsedOptions().pd_method),
        nMajorSpecies(config.GetFlameletParsedOptions().n_pd_major_species),
        pdTermsPerCV(FlameletPDTermsPerCV(config.GetFlameletParsedOptions().n_pd_major_species)) {}

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

    if (pdMethod == FLAMELET_PD_METHOD::SOURCE_TERM) {
      sourceTermFluxes(idx, opt, iPoint, side_i, jPoint, side_j, rho, normal, vector_ij, res);
      return;
    }

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
      if (!opt.oneSided) res.flux_j(iScalar) += D * projGrad;

      if (opt.implicit) {
        res.jac_ii(iScalar, iScalar) += D * proj_on_rho_i;
        if (!opt.oneSided) {
          res.jac_ij(iScalar, iScalar) -= D * proj_on_rho_j;
          res.jac_ji(iScalar, iScalar) -= D * proj_on_rho_i;
          res.jac_jj(iScalar, iScalar) += D * proj_on_rho_j;
        }
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
    if (!opt.oneSided) res.flux_j(I_ENTH) += Dth * projGradT;
  }

 private:
  const bool preferentialDiffusion;
  const unsigned short nControlVars;
  const FLAMELET_PD_METHOD pdMethod;
  const unsigned short nMajorSpecies;
  const unsigned short pdTermsPerCV;

  /*!
   * \brief Resolved preferential diffusion fluxes of the SOURCE_TERM method, Eq. (14) of
   *        Schepers & van Oijen (C&F 280, 2025): for each controlling variable phi_k,
   *          J_{phi_k} = - sum_sp D_{phi_k,sp} grad(Y_sp) - (D^T_{phi_k}/T) grad(T),
   *        entering the transport equation as -div(J). The major species mass fractions are the
   *        solver's auxiliary variables, so their CFD-resolved gradients drive the flux, which is
   *        the point of tabulating a coefficient per species instead of pre-contracting them
   *        against the one-dimensional flamelet gradients.
   * \note The tabulated coefficients are already density weighted (kg/(m s), or W/m for the
   *       enthalpy row), unlike the kinematic diffusivity the ordinary term uses, so they are
   *       averaged directly and NOT multiplied by rho. The 1/T of the thermal flux is not
   *       tabulated and is applied here with the local temperature.
   * \note Explicit: the coefficients and the major species are manifold lookups, so there is no
   *       exact Jacobian with respect to the transported scalars. The implicit part is a
   *       stabilising self-diffusion added to the Jacobian alone, below.
   */
  template <class VariableType>
  FORCEINLINE void sourceTermFluxes(const FlowIndices& idx, const ScalarFluxOptions& opt, Int iPoint,
                                    const EdgeSide<VariableType>& side_i, Int jPoint,
                                    const EdgeSide<VariableType>& side_j, const CPair<Double>& rho,
                                    const Vector<Double, nDim>& normal, const Vector<Double, nDim>& vector_ij,
                                    EdgeResidual<Double, nVar>& res) const {
    const Double dist2_ij = fmax(squaredNorm(vector_ij), EPS);
    const Double proj_vector_ij = dot(vector_ij, normal) / dist2_ij;

    const Double T_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.Temperature());
    const Double T_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.Temperature());
    const auto gradT_i = gatherVariables<nDim>(iPoint, side_i.flowNodes->GetGradient_Primitive(), idx.Temperature());
    const auto gradT_j = gatherVariables<nDim>(jPoint, side_j.flowNodes->GetGradient_Primitive(), idx.Temperature());

    /*--- The scalar solvers are templated on CSpeciesVariable, which does not carry the Eq. (14)
     coefficients. Only the flamelet solver instantiates this flux, and only on its own interior
     edges, so its container is what these sides hold. ---*/
    const auto& pdCoeff_i = static_cast<const CSpeciesFlameletVariable&>(side_i.scalarNodes).GetPDFluxCoeffs();
    const auto& pdCoeff_j = static_cast<const CSpeciesFlameletVariable&>(side_j.scalarNodes).GetPDFluxCoeffs();

    for (auto iScalar = 0u; iScalar < nControlVars; ++iScalar) {
      /*--- Magnitude of the Eq. (14) flux on each side, accumulated as the terms are applied, for
       the stabilising diffusivity below. ---*/
      Double fluxMag_i = 0.0, fluxMag_j = 0.0;

      /*--- Molecular terms, one per major species. ---*/
      for (auto iSp = 0u; iSp < nMajorSpecies; ++iSp) {
        const auto col = iScalar * pdTermsPerCV + iSp;

        const Double Y_i = gatherVariables(iPoint, side_i.scalarNodes.GetAuxVar(), iSp);
        const Double Y_j = gatherVariables(jPoint, side_j.scalarNodes.GetAuxVar(), iSp);
        const auto gradY_i = gatherVariables<nDim>(iPoint, side_i.scalarNodes.GetAuxVarGradient(), iSp);
        const auto gradY_j = gatherVariables<nDim>(jPoint, side_j.scalarNodes.GetAuxVarGradient(), iSp);

        const Double c_i = gatherVariables(iPoint, pdCoeff_i, col);
        const Double c_j = gatherVariables(jPoint, pdCoeff_j, col);
        const Double D = 0.5 * (c_i + c_j);

        const Double projGrad = projectedGradient(opt, gradY_i, gradY_j, Y_i, Y_j, normal, vector_ij, dist2_ij);

        res.flux_i(iScalar) -= D * projGrad;
        if (!opt.oneSided) res.flux_j(iScalar) += D * projGrad;

        if (opt.implicit) {
          fluxMag_i += fabs(c_i) * sqrt(fmax(squaredNorm(gradY_i), EPS));
          fluxMag_j += fabs(c_j) * sqrt(fmax(squaredNorm(gradY_j), EPS));
        }
      }

      /*--- Thermal (Soret) term. ---*/
      const auto colT = iScalar * pdTermsPerCV + FlameletPDThermalTerm(nMajorSpecies);
      const Double Dth = 0.5 * (gatherVariables(iPoint, pdCoeff_i, colT) / T_i +
                                gatherVariables(jPoint, pdCoeff_j, colT) / T_j);

      const Double projGradT = projectedGradient(opt, gradT_i, gradT_j, T_i, T_j, normal, vector_ij, dist2_ij);

      res.flux_i(iScalar) -= Dth * projGradT;
      if (!opt.oneSided) res.flux_j(iScalar) += Dth * projGradT;

      /*--- Stabilisation. The fluxes above are explicit, and explicit diffusion carries a limit on
       the pseudo time step that the flame front cells reach first. A self-diffusion of the
       controlling variable is added to the Jacobian ALONE: the residual still holds exactly the
       Eq. (14) fluxes, so the converged solution is untouched, while the positive definite
       operator lifts that limit the way an implicit treatment of a real self-diffusion would.
       D_stab is the baseline diffusivity scaled by how much larger the Eq. (14) coefficients are,
       capped so a vanishing baseline cannot make it unbounded. ---*/
      if (opt.implicit) {
        const Double Dbase_i = rho.i * gatherVariables(iPoint, side_i.scalarNodes.GetDiffusivity(), iScalar);
        const Double Dbase_j = rho.j * gatherVariables(jPoint, side_j.scalarNodes.GetDiffusivity(), iScalar);

        /*--- Recast the flux bound as a self diffusivity, |J| <= D_stab |grad phi|, so the implicit
         operator dominates the explicit flux without damping more than the flux itself warrants.
         Scaling by the control variable's own gradient is what keeps it that tight: dropping it
         leaves D_stab at the raw coefficient magnitude, which over damps and cancels the very CFL
         gain the term exists to provide. Where grad(phi) vanishes while the cross gradients do not
         the bound degenerates, so it is capped against the baseline diffusion. ---*/
        fluxMag_i += fabs(gatherVariables(iPoint, pdCoeff_i, colT) / T_i) * sqrt(fmax(squaredNorm(gradT_i), EPS));
        fluxMag_j += fabs(gatherVariables(jPoint, pdCoeff_j, colT) / T_j) * sqrt(fmax(squaredNorm(gradT_j), EPS));

        const auto gradPhi_i = gatherVariables<nDim>(iPoint, side_i.scalarNodes.GetGradient(), iScalar);
        const auto gradPhi_j = gatherVariables<nDim>(jPoint, side_j.scalarNodes.GetGradient(), iScalar);

        const Double Dstab_i = fmin(fluxMag_i / sqrt(fmax(squaredNorm(gradPhi_i), EPS)), C_STAB_MAX * Dbase_i);
        const Double Dstab_j = fmin(fluxMag_j / sqrt(fmax(squaredNorm(gradPhi_j), EPS)), C_STAB_MAX * Dbase_j);
        const Double Dstab = 0.5 * (Dstab_i + Dstab_j);

        res.jac_ii(iScalar, iScalar) += Dstab * proj_vector_ij / rho.i;
        if (!opt.oneSided) {
          res.jac_ij(iScalar, iScalar) -= Dstab * proj_vector_ij / rho.j;
          res.jac_ji(iScalar, iScalar) -= Dstab * proj_vector_ij / rho.i;
          res.jac_jj(iScalar, iScalar) += Dstab * proj_vector_ij / rho.j;
        }
      }
    }
  }

  /*!
   * \brief Cap on the stabilising diffusivity, as a multiple of the baseline diffusion.
   */
  static constexpr passivedouble C_STAB_MAX = 100.0;

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
