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
 * \brief Convection and diffusion of the Menter SST model, conservative with a coupled (but
 *        neither symmetric nor diagonal) 2x2 diffusion matrix.
 * \note SST writes no finalizeFlux of its own: the inherited CUpwScalarFlux one is exactly
 *       flux(iVar) = a0*rho_i*phi_i(iVar) + a1*rho_j*phi_j(iVar), Conservative weighting by
 *       density, which is the whole of the model's old convective term.
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
   * \brief Diffusion coefficients of both orientations of the edge, see CAvgGrad_TurbSST.
   * \note The old discretization evaluates the diffusion numerics twice per edge, once with i
   *       first and once with j first (CScalarSolver::Viscous_Residual_NonCons), because the
   *       cross term below reads the transported omega of whichever point was passed first. That
   *       is exactly what returning two different matrices here reproduces: D.i, read by i's row,
   *       uses omega at i; D.j, read by j's row, uses omega at j. Every other entry is symmetric
   *       (an i/j average), so it is the same in both matrices.
   */
  template <class VariableType>
  FORCEINLINE CPair<Matrix<Double, nVar, nVar>> coefficients(const FlowIndices& idx, Int iPoint,
                                                             const EdgeSide<VariableType>& side_i, Int jPoint,
                                                             const EdgeSide<VariableType>& side_j) const {
    const Double rho_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.Density());
    const Double rho_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.Density());
    const Double mu_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.LaminarViscosity());
    const Double mu_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.LaminarViscosity());
    const Double muT_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.EddyViscosity());
    const Double muT_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.EddyViscosity());

    const Double F1_i = side_i.scalarNodes.GetF1blending(iPoint);
    const Double F1_j = side_j.scalarNodes.GetF1blending(jPoint);
    const Double omega_i = gatherVariables(iPoint, side_i.scalarNodes.GetSolution(), 1);
    const Double omega_j = gatherVariables(jPoint, side_j.scalarNodes.GetSolution(), 1);

    const Double sigma_kine_i = F1_i * sigma_k1 + (1.0 - F1_i) * sigma_k2;
    const Double sigma_kine_j = F1_j * sigma_k1 + (1.0 - F1_j) * sigma_k2;
    const Double sigma_omega_i = F1_i * sigma_om1 + (1.0 - F1_i) * sigma_om2;
    const Double sigma_omega_j = F1_j * sigma_om1 + (1.0 - F1_j) * sigma_om2;

    const Double diff_kine = 0.5 * ((mu_i + sigma_kine_i * muT_i) + (mu_j + sigma_kine_j * muT_j));
    const Double diff_omega = 0.5 * ((mu_i + sigma_omega_i * muT_i) + (mu_j + sigma_omega_j * muT_j));

    const Double lambda_i = 2.0 * (1.0 - F1_i) * rho_i * sigma_omega_i;
    const Double lambda_j = 2.0 * (1.0 - F1_j) * rho_j * sigma_omega_j;
    const Double lambda_ij = 0.5 * (lambda_i + lambda_j);
    const Double w_ij = 0.5 * (omega_i + omega_j);

    /*--- Same two terms as CAvgGrad_TurbSST's old diff_omega_T2/diff_omega_T3, kept as the exact
     * same two additions (not algebraically simplified) so this rounds identically. ---*/
    const Double diff_omega_T2 = lambda_ij;
    const Double diff_omega_T3_i = -omega_i * lambda_ij / w_ij;
    const Double diff_omega_T3_j = -omega_j * lambda_ij / w_ij;

    Matrix<Double, nVar, nVar> D_i = Double(0.0), D_j = Double(0.0);
    D_i(0, 0) = diff_kine;
    D_i(1, 1) = diff_omega;
    D_i(1, 0) = diff_omega_T2 + diff_omega_T3_i;

    D_j(0, 0) = diff_kine;
    D_j(1, 1) = diff_omega;
    D_j(1, 0) = diff_omega_T2 + diff_omega_T3_j;

    return {D_i, D_j};
  }

  /*!
   * \brief Extra Jacobian terms from the dependence of the cross-diffusion coefficient on omega.
   * \note Unlike SA's, this correction is not a per-edge constant: it comes out of the same
   *       twice-per-edge evaluation coefficients() reproduces, so it lands on all four Jacobian
   *       blocks rather than mirroring the D.i/D.j split (jac_ii and jac_ji share one term, jac_ij
   *       and jac_jj the other) -- worked out by hand against the old i->j / j->i pair of
   *       Viscous_Residual_NonCons calls, not read off the single-evaluation shape most other
   *       models have.
   */
  template <class VariableType, size_t Size>
  FORCEINLINE void coefficientJacobians(const FlowIndices& idx, Int iPoint, const EdgeSide<VariableType>& side_i,
                                        Int jPoint, const EdgeSide<VariableType>& side_j,
                                        const Vector<Double, Size>& projGrad, EdgeResidual<Double, nVar>& res) const {
    const Double rho_i = gatherVariables(iPoint, side_i.flowNodes->GetPrimitive(), idx.Density());
    const Double rho_j = gatherVariables(jPoint, side_j.flowNodes->GetPrimitive(), idx.Density());
    const Double F1_i = side_i.scalarNodes.GetF1blending(iPoint);
    const Double F1_j = side_j.scalarNodes.GetF1blending(jPoint);
    const Double omega_i = gatherVariables(iPoint, side_i.scalarNodes.GetSolution(), 1);
    const Double omega_j = gatherVariables(jPoint, side_j.scalarNodes.GetSolution(), 1);

    const Double sigma_omega_i = F1_i * sigma_om1 + (1.0 - F1_i) * sigma_om2;
    const Double sigma_omega_j = F1_j * sigma_om1 + (1.0 - F1_j) * sigma_om2;
    const Double lambda_i = 2.0 * (1.0 - F1_i) * rho_i * sigma_omega_i;
    const Double lambda_j = 2.0 * (1.0 - F1_j) * rho_j * sigma_omega_j;
    const Double lambda_ij = 0.5 * (lambda_i + lambda_j);

    const Double denom = pow(omega_i + omega_j, 2.0);
    const Double E_i = 2.0 * lambda_ij * omega_i / denom * projGrad(0);
    const Double E_j = 2.0 * lambda_ij * omega_j / denom * projGrad(0);

    res.jac_ii(1, 1) += E_j;
    res.jac_ij(1, 1) -= E_i;
    res.jac_ji(1, 1) += E_j;
    res.jac_jj(1, 1) -= E_i;
  }
};
