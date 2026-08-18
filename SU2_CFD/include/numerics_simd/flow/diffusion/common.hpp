/*!
 * \file common.hpp
 * \brief Helper functions for viscous methods.
 * \author P. Gomes, C. Pederson, A. Bueno, F. Palacios, T. Economon
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

#include "../../CNumericsSIMD.hpp"
#include "../../util.hpp"
#include "../variables.hpp"
#include "../../../numerics/CNumerics.hpp"

/*!
 * \brief Average gradients at i/j points.
 */
template<size_t nVar, size_t nDim, class GradientType>
FORCEINLINE MatrixDbl<nVar,nDim> averageGradient(const Int& iPoint, const Int& jPoint,
                                                 const GradientType& gradient) {
  auto avgGrad = gatherVariables<nVar,nDim>(iPoint, gradient);
  auto grad_j = gatherVariables<nVar,nDim>(jPoint, gradient);
  for (size_t iVar = 0; iVar < nVar; ++iVar) {
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      avgGrad(iVar,iDim) *= 0.5;
      avgGrad(iVar,iDim) += 0.5 * grad_j(iVar,iDim);
    }
  }
  return avgGrad;
}

template<size_t nSecVar, class SecondaryType>
FORCEINLINE CCompressibleSecondary<nSecVar> averageSecondary(Int iPoint, Int jPoint,
                                                             const SecondaryType& secondary) {

    CCompressibleSecondary<nSecVar> out;
    auto second_i = gatherVariables<nSecVar>(iPoint, secondary);
    auto second_j = gatherVariables<nSecVar>(jPoint, secondary);
    for (size_t iVar = 0; iVar < nSecVar; ++iVar) {
        out.all(iVar) = 0.5 * (second_i(iVar) + second_j(iVar));
    }
    return out;
}

/*!
 * \brief Correct average gradient with the directional derivative to avoid decoupling.
 */
template<size_t nVar, size_t nDim, class PrimitiveType>
FORCEINLINE void correctGradient(const PrimitiveType& V,
                                 const VectorDbl<nDim>& vector_ij,
                                 const VectorDbl<nDim>& diss,
                                 MatrixDbl<nVar,nDim>& avgGrad) {
  for (size_t iVar = 0; iVar < nVar; ++iVar) {
//    Double corr = ( V.j.all(iVar) - V.i.all(iVar));
    Double corr = ( V.j.all(iVar) - V.i.all(iVar) - dot(avgGrad[iVar],vector_ij));
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      avgGrad(iVar,iDim) += corr * diss(iDim);
    }
  }
}

/*!
 * \brief Compute the stress tensor.
 * \note Second viscosity term ignored.
 */
template<size_t nVar, size_t nDim>
FORCEINLINE MatrixDbl<nDim> stressTensor(const Double& viscosity,
                                         const MatrixDbl<nVar,nDim>& grad) {
  /*--- Hydrostatic term. ---*/
  Double velDiv = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    velDiv += grad(iDim+1,iDim);
  }
  Double pTerm = 2.0/3.0 * viscosity * velDiv;

  MatrixDbl<nDim> tau;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    /*--- Deviatoric term. ---*/
    for (size_t jDim = 0; jDim < nDim; ++jDim) {
      tau(iDim,jDim) = viscosity * (grad(jDim+1,iDim) + grad(iDim+1,jDim));
    }
    tau(iDim,iDim) -= pTerm;
  }
  return tau;
}

/*!
 * \brief Add perturbed stress tensor.
 * \note Not inlined because it is not easy to vectorize properly, due to tred2 and tql2.
 */
template<class PrimitiveType, class MatrixType, size_t nDim, class... Ts>
NEVERINLINE void addPerturbedRSM(const PrimitiveType& V,
                                 const MatrixType& grad,
                                 const Double& turb_ke,
                                 MatrixDbl<nDim,nDim>& tau,
                                 Ts... uq_args) {
  /*--- Handle SIMD dimensions 1 by 1. ---*/
  for (size_t k = 0; k < Double::Size; ++k) {
    su2double velgrad[nDim][nDim];
    for (size_t iVar = 0; iVar < nDim; ++iVar)
      for (size_t iDim = 0; iDim < nDim; ++iDim)
        velgrad[iVar][iDim] = grad(iVar+1,iDim)[k];

    su2double rsm[3][3];
    CNumerics::ComputePerturbedRSM(nDim, uq_args..., velgrad, V.density()[k],
                                   V.eddyVisc()[k], turb_ke[k], rsm);

    for (size_t iDim = 0; iDim < nDim; ++iDim)
      for (size_t jDim = 0; jDim < nDim; ++jDim)
        tau(iDim,jDim)[k] -= V.density()[k] * rsm[iDim][jDim];
  }
}

/*!
 * \brief SA-QCR2000 modification of the stress tensor.
 */
template<class MatrixType, size_t nDim>
FORCEINLINE void addQCR(const MatrixType& grad, MatrixDbl<nDim>& tau) {
  constexpr passivedouble c_cr1 = 0.3;

  /*--- Denominator, antisymmetric normalized rotation tensor. ---*/
  Double denom = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim)
    for (size_t jDim = 0; jDim < nDim; ++jDim)
      denom += grad(iDim+1,jDim) * grad(iDim+1,jDim);

  const Double factor = 1 / sqrt(fmax(denom,1e-10));

  /*--- Compute the QCR term, and update the stress tensor. ---*/
  MatrixDbl<nDim> qcr;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    for (size_t jDim = 0; jDim < nDim; ++jDim) {
      qcr(iDim,jDim) = 0.0;
      for (size_t kDim = 0; kDim < nDim; ++kDim) {
        Double O_ik = (grad(iDim+1,kDim) - grad(kDim+1,iDim)) * factor;
        Double O_jk = (grad(jDim+1,kDim) - grad(kDim+1,jDim)) * factor;
        qcr(iDim,jDim) += O_ik*tau(jDim,kDim) + O_jk*tau(iDim,kDim);
      }
    }
  }
  for (size_t iDim = 0; iDim < nDim; ++iDim)
    for (size_t jDim = 0; jDim < nDim; ++jDim)
      tau(iDim,jDim) -= c_cr1 * qcr(iDim,jDim);
}

/*!
 * \brief Scale the stress tensor according to the target (from a
 *        wall function) magnitude in the tangential direction.
 */
template<class Container, size_t nDim>
FORCEINLINE void addTauWall(const Int& iPoint, const Int& jPoint,
                            const Container& tauWall,
                            const VectorDbl<nDim>& unitNormal,
                            MatrixDbl<nDim>& tau) {

  Double tauWall_i = fmax(gatherVariables(iPoint, tauWall), 0.0);
  Double tauWall_j = fmax(gatherVariables(jPoint, tauWall), 0.0);

  Double isWall_i = tauWall_i > 0.0;
  Double isWall_j = tauWall_j > 0.0;
  /*--- Arithmetic xor. ---*/
  Double isNormalEdge = isWall_i+isWall_j - 2*isWall_i*isWall_j;

  /*--- Tau wall is 0 for edges that are not normal to walls. ---*/
  Double tauWall_ij = (tauWall_i+tauWall_j) * isNormalEdge;

  /*--- Scale is 1 for those edges, i.e. tau is not changed. ---*/
  Double scale =
      tauWall_ij / fmax(norm(tangentProjection(tau,unitNormal)), EPS) + (1.0-isNormalEdge);

  for (size_t iDim = 0; iDim < nDim; ++iDim)
    for (size_t jDim = 0; jDim < nDim; ++jDim)
      tau(iDim,jDim) *= scale;
}

/*!
 * \brief Jacobian of the stress tensor (compressible flow).
 */
template<size_t nVar, size_t nDim, class PrimitiveType>
FORCEINLINE MatrixDbl<nDim,nVar> stressTensorJacobian(const PrimitiveType& V,
                                                      const VectorDbl<nDim>& normal,
                                                      const Double& dist_ij) {
  Double viscosity = V.laminarVisc() + V.eddyVisc();
  Double xi = viscosity / (V.density() * dist_ij);
  MatrixDbl<nDim,nVar> jac;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    /*--- Momentum. ---*/
    for (size_t jDim = 0; jDim < nDim; ++jDim) {
      jac(iDim,jDim+1) = (-1/3.0) * xi * normal(iDim) * normal(jDim);
    }
    jac(iDim,iDim+1) -= xi;
    /*--- Density. ---*/
    jac(iDim,0) = -dot<nDim>(&jac(iDim,1), V.velocity());
    /*--- Energy. ---*/
    jac(iDim,nDim+1) = 0.0;
  }
  return jac;
}

/*!
 * \brief Viscous flux for compressible flows.
 */
template<size_t nVar, size_t nDim, class PrimitiveType>
FORCEINLINE VectorDbl<nVar> viscousFlux(const PrimitiveType& V,
                                        const MatrixDbl<nDim>& tau,
                                        const VectorDbl<nDim>& heatFlux,
                                        const VectorDbl<nDim>& normal) {
  VectorDbl<nVar> flux;
  flux(0) = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    /*--- Using the symmetry of the tensor. ---*/
    flux(iDim+1) = dot(tau[iDim], normal);
  }
  flux(nDim+1) = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    auto viscWork = dot<nDim>(tau[iDim], V.velocity());
    flux(nDim+1) += normal(iDim) * (heatFlux(iDim) + viscWork);
  }
  return flux;
}


/*!
 * \brief Exact Jacobian of the projected viscous flux (ideal gas), with the average LSQ
 * gradient contribution treated as frozen: the state dependence enters through the
 * corrected-gradient terms diss*(V_j - V_i), the averaged velocity in the viscous work,
 * and (via the averaged secondary variables) mu(T) and kt(T). Exact when the LSQ part of
 * the corrected gradient is zero.
 * \note Outputs d(total flux)/dU = -d(viscous flux)/dU, the caller adds area*jac to the
 * convective Jacobians. This must remain algebraically identical to the scalar twin
 * (CAvgGrad_Flow::SetExactViscousProjJacs in numerics/flow/flow_diffusion.cpp).
 */
template<size_t nVar, size_t nDim, class PrimitiveType>
FORCEINLINE void viscousFluxJacobian(const Double& gamma,
                                     const PrimitiveType& V,
                                     const PrimitiveType& Vi,
                                     const PrimitiveType& Vj,
                                     const VectorDbl<nDim>& unitNormal,
                                     const VectorDbl<nDim>& diss,
                                     const MatrixDbl<nVar,nDim>& grad,
                                     const Double& Pr_t,
                                     const CCompressibleSecondary<8>& secondary,
                                     MatrixDbl<nDim+2>& jac_i,
                                     MatrixDbl<nDim+2>& jac_j) {

  constexpr size_t nEqn = nDim+2;
  const auto& n = unitNormal;
  const Double two_third = 2.0/3.0;

  const Double mu = V.laminarVisc() + V.eddyVisc();
  const Double cond = V.thermalCond() + V.cp()*V.eddyVisc()/Pr_t;

  /*--- Corrected mean strain tensor S_kl = du_k/dx_l + du_l/dx_k - 2/3 delta_kl div(u),
   * its normal projection sn, viscous work factor snv (with the averaged velocity),
   * normal T-gradient, and normal component of the correction direction. ---*/

  Double trace = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) trace += grad(1+iDim,iDim);

  VectorDbl<nDim> sn;
  Double snv = 0.0, gradT_n = 0.0, dn = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    sn(iDim) = 0.0;
    for (size_t jDim = 0; jDim < nDim; ++jDim) {
      Double S_ij = grad(1+iDim,jDim) + grad(1+jDim,iDim);
      if (iDim == jDim) S_ij -= two_third*trace;
      sn(iDim) += S_ij*n(jDim);
    }
    snv += sn(iDim)*V.velocity(iDim);
    gradT_n += grad(0,iDim)*n(iDim);
    dn += diss(iDim)*n(iDim);
  }

  /*--- Averaged temperature-derivatives of the transport properties (zero for
   * constant-property models). ---*/

  const Double dmudT = secondary.dmudT_rho();
  const Double dktdT = secondary.dktdT_rho();

  for (int side = 0; side < 2; ++side) {

    const auto& Vn = (side == 0) ? Vi : Vj;
    auto& jac = (side == 0) ? jac_i : jac_j;
    const Double sgn = (side == 0) ? -1.0 : 1.0;

    const Double T = Vn.temperature();
    const Double p = Vn.pressure();
    const Double rho = Vn.density();
    Double q2 = 0.0;
    for (size_t iDim = 0; iDim < nDim; ++iDim) q2 += Vn.velocity(iDim)*Vn.velocity(iDim);

    /*--- Nodal T = p/(R*rho) exactly for the ideal gas. ---*/
    const Double dT_drho = -T/rho;
    const Double dT_dp = T/p;

    /*--- The 0.5 accounts for d(mean value)/d(nodal value). ---*/
    const Double dmu_dT_side = 0.5*dmudT;
    const Double dcond_dT_side = 0.5*dktdT;

    /*--- Jacobian of the projected flux F = [0, mu*sn, mu*snv + cond*gradT_n] with
     * respect to this side's primitives w = [rho, u_k, p]. ---*/

    MatrixDbl<nEqn> dFdw;
    for (size_t r = 0; r < nEqn; ++r)
      for (size_t c = 0; c < nEqn; ++c)
        dFdw(r,c) = 0.0;

    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      dFdw(iDim+1,0) = dmu_dT_side*dT_drho*sn(iDim);
      for (size_t kDim = 0; kDim < nDim; ++kDim) {
        const Double delta = (iDim == kDim) ? 1.0 : 0.0;
        const Double dsn_du = sgn*(dn*delta + diss(iDim)*n(kDim) - two_third*n(iDim)*diss(kDim));
        dFdw(iDim+1,kDim+1) = mu*dsn_du;
      }
      dFdw(iDim+1,nDim+1) = dmu_dT_side*dT_dp*sn(iDim);
    }

    const Double dqn_dT = sgn*cond*dn + dcond_dT_side*gradT_n;
    dFdw(nDim+1,0) = (dmu_dT_side*snv + dqn_dT)*dT_drho;
    for (size_t kDim = 0; kDim < nDim; ++kDim) {
      dFdw(nDim+1,kDim+1) = 0.5*mu*sn(kDim);
      for (size_t iDim = 0; iDim < nDim; ++iDim)
        dFdw(nDim+1,kDim+1) += dFdw(iDim+1,kDim+1)*V.velocity(iDim);
    }
    dFdw(nDim+1,nDim+1) = (dmu_dT_side*snv + dqn_dT)*dT_dp;

    /*--- Transform to conservative variables using this side's state and negate (the
     * viscous flux is subtracted from the total flux). dW/dU rows: rho: (1, 0, 0);
     * u_k: (-u_k/rho, delta/rho, 0); p: ((g-1)*q2/2, -(g-1)*u, (g-1)).
     * Row 0 of the flux is zero. ---*/

    for (size_t r = 0; r < nEqn; ++r) {
      Double col0 = dFdw(r,0) + dFdw(r,nDim+1)*0.5*(gamma-1)*q2;
      for (size_t kDim = 0; kDim < nDim; ++kDim) col0 -= dFdw(r,kDim+1)*Vn.velocity(kDim)/rho;
      jac(r,0) = -col0;
      for (size_t c = 0; c < nDim; ++c)
        jac(r,c+1) = -(dFdw(r,c+1)/rho - dFdw(r,nDim+1)*(gamma-1)*Vn.velocity(c));
      jac(r,nDim+1) = -dFdw(r,nDim+1)*(gamma-1);
    }
  }
}
