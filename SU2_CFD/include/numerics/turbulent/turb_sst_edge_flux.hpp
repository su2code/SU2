/*!
 * \file turb_sst_edge_flux.hpp
 * \brief Menter SST model as a third-layer scalar flux, see numerics/scalar/scalar_edge_flux.hpp.
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
 * \class CScalarFlux_SST
 * \ingroup ViscDiscr
 * \brief Convection and diffusion of the Menter SST model, conservative with a coupled (but
 *        neither symmetric nor diagonal) 2x2 diffusion matrix.
 * \note SST writes no finalizeFlux of its own: the inherited CUpwScalarFlux one is exactly
 *       flux(iVar) = a0*rho_i*phi_i(iVar) + a1*rho_j*phi_j(iVar), Conservative weighting by
 *       density, which is the model's whole convective term.
 */
template <class Double, class FlowIndices, int nDim, size_t nVar = 2>
class CScalarFlux_SST
    : public CUpwScalarBase<Double, CScalarFlux_SST<Double, FlowIndices, nDim, nVar>, FlowIndices, nDim, nVar> {
 public:
  static constexpr bool Conservative = true;
  static constexpr bool DiagonalDiffusion = false;

  using Base = CUpwScalarBase<Double, CScalarFlux_SST, FlowIndices, nDim, nVar>;
  using Int = typename Base::Int;
  using Base::Base;

 private:
  /*--- Fixed regardless of SST_OPTIONS::version: only the production-limiter and source-term
   * constants (alfa/gamma) differ by version, not these. ---*/
  static constexpr passivedouble sigma_k1 = 0.85;
  static constexpr passivedouble sigma_k2 = 1.0;
  static constexpr passivedouble sigma_om1 = 0.5;
  static constexpr passivedouble sigma_om2 = 0.856;

 public:
  /*!
   * \brief Diffusion coefficients of both orientations of the edge, and the terms of the cross
   *        diffusion that the Jacobian correction below needs, so that neither the gathers nor
   *        the blending are repeated for it.
   * \note The cross term reads the transported omega of whichever point its row is being written
   *       for, so it is not symmetric: i, read by i's row, uses omega at i; j, read by j's row,
   *       uses omega at j. Every other entry is an i/j average, so it is the same in both.
   */
  struct CCoefficients {
    Matrix<Double, nVar, nVar> i, j;
    Double lambda_ij, omega_i, omega_j;
  };

  template <class VariableType>
  FORCEINLINE CCoefficients coefficients(const FlowIndices& idx, Int iPoint, const EdgeSide<VariableType>& side_i,
                                         Int jPoint, const EdgeSide<VariableType>& side_j,
                                         const CPair<Double>& rho) const {
    const Double mu_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.LaminarViscosity());
    const Double mu_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.LaminarViscosity());
    const Double muT_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.EddyViscosity());
    const Double muT_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.EddyViscosity());

    const Double F1_i = gatherVariables(iPoint, side_i.scalarNodes.GetF1blending());
    const Double F1_j = gatherVariables(jPoint, side_j.scalarNodes.GetF1blending());

    CCoefficients D;
    D.omega_i = gatherVariables(iPoint, side_i.scalarNodes.GetSolution(), 1);
    D.omega_j = gatherVariables(jPoint, side_j.scalarNodes.GetSolution(), 1);

    const Double sigma_kine_i = F1_i * sigma_k1 + (1.0 - F1_i) * sigma_k2;
    const Double sigma_kine_j = F1_j * sigma_k1 + (1.0 - F1_j) * sigma_k2;
    const Double sigma_omega_i = F1_i * sigma_om1 + (1.0 - F1_i) * sigma_om2;
    const Double sigma_omega_j = F1_j * sigma_om1 + (1.0 - F1_j) * sigma_om2;

    const Double diff_kine = 0.5 * ((mu_i + sigma_kine_i * muT_i) + (mu_j + sigma_kine_j * muT_j));
    const Double diff_omega = 0.5 * ((mu_i + sigma_omega_i * muT_i) + (mu_j + sigma_omega_j * muT_j));

    const Double lambda_i = 2.0 * (1.0 - F1_i) * rho.i * sigma_omega_i;
    const Double lambda_j = 2.0 * (1.0 - F1_j) * rho.j * sigma_omega_j;
    D.lambda_ij = 0.5 * (lambda_i + lambda_j);
    const Double w_ij = 0.5 * (D.omega_i + D.omega_j);

    /*--- Cross-diffusion coefficient: a divergence-theorem term (diff_omega_T2) plus a cell
     * centre correction (diff_omega_T3) that reads the transported omega of the row's own point. ---*/
    const Double diff_omega_T2 = D.lambda_ij;
    const Double diff_omega_T3_i = -D.omega_i * D.lambda_ij / w_ij;
    const Double diff_omega_T3_j = -D.omega_j * D.lambda_ij / w_ij;

    /*--- D.i(0,1) and D.j(0,1) are left zero: there is no diffusive coupling from omega into
     * the k row. ---*/
    D.i = Double(0.0);
    D.j = Double(0.0);
    D.i(0, 0) = diff_kine;
    D.i(1, 1) = diff_omega;
    D.i(1, 0) = diff_omega_T2 + diff_omega_T3_i;

    D.j(0, 0) = diff_kine;
    D.j(1, 1) = diff_omega;
    D.j(1, 0) = diff_omega_T2 + diff_omega_T3_j;

    return D;
  }

  /*!
   * \brief Extra Jacobian terms from the dependence of the cross-diffusion coefficient on omega.
   * \note diff_omega_T3_i and diff_omega_T3_j both depend on omega_i and omega_j through w_ij, so
   *       each of the four blocks needs a correction beyond the one diffusionTerms already applies
   *       through projGrad. The correction only depends on which point's omega is being
   *       differentiated against, not on which row it lands in: differentiating against omega_i
   *       gives +E_j in both jac_ii and jac_ji, differentiating against omega_j gives -E_i in both
   *       jac_ij and jac_jj.
   */
  template <size_t Size>
  FORCEINLINE void coefficientJacobians(const ScalarFluxOptions& opt, const CCoefficients& D,
                                        const Vector<Double, Size>& projGrad, EdgeResidual<Double, nVar>& res) const {
    const Double denom = pow(D.omega_i + D.omega_j, 2.0);
    const Double E_i = 2.0 * D.lambda_ij * D.omega_i / denom * projGrad(0);
    const Double E_j = 2.0 * D.lambda_ij * D.omega_j / denom * projGrad(0);

    res.jac_ii(1, 1) += E_j;
    if (opt.oneSided) return;

    res.jac_ij(1, 1) -= E_i;
    res.jac_ji(1, 1) += E_j;
    res.jac_jj(1, 1) -= E_i;
  }
};
