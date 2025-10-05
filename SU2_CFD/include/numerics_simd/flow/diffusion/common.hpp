/*!
 * \file common.hpp
 * \brief Helper functions for viscous methods.
 * \author P. Gomes, C. Pederson, A. Bueno, F. Palacios, T. Economon
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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
FORCEINLINE MatrixDbl<nVar,nDim> averageGradient(Int iPoint, Int jPoint,
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
FORCEINLINE MatrixDbl<nDim> stressTensor(Double viscosity,
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
FORCEINLINE void addTauWall(Int iPoint, Int jPoint,
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
                                                      Double dist_ij) {
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
 * \brief Viscous flux for compressible flows.
 */
template<size_t nVar, size_t nDim, class PrimitiveType>
FORCEINLINE void viscousFluxJacobian(const PrimitiveType& V,
                                     const PrimitiveType& Vi,
                                     const PrimitiveType& Vj,
                             const VectorDbl<nDim>& unitNormal,
                             const VectorDbl<nDim>& diss,
                             const MatrixDbl<nVar,nDim>& grad,
                             const Double& cp, const Double& Pr_l, const Double& Pr_t,
                             const CCompressibleSecondary<8>& secondary,
                             MatrixDbl<nDim+2>& jac_i,
                             MatrixDbl<nDim+2>& jac_j) {

  //! Arithmetic averages of velocities and temperature.
  //
  //        u = half*(uL + uR);
  //        v = half*(vL + vR);
  //        w = half*(wL + wR);
  //        T = half*(TL + TR);
  int iDim, jDim, iVar, jVar;

  const Double gamma = 1.4;

  VectorDbl<nDim> vel;
  for(iDim = 0; iDim < nDim; ++iDim) vel(iDim) = V.velocity(iDim);
  const auto& rho = V.density();
  const auto& T = V.temperature();
  const auto& p = V.pressure();
  const auto& H = V.enthalpy();
  const auto& a = V.speedSound();
  const auto& q2 = squaredNorm<nDim>(V.velocity());
  const auto& qn  = dot(V.velocity(),unitNormal);

  VectorDbl<nDim> velL;
  for(iDim = 0; iDim < nDim; ++iDim) velL(iDim) = Vi.velocity(iDim);
  const auto& rhoL = Vi.density();
  const auto& TL = Vi.temperature();
  const auto& pL = Vi.pressure();
  const auto& HL = Vi.enthalpy();
  const auto& aL = Vi.speedSound();
  const auto& q2L = squaredNorm<nDim>(Vi.velocity());
  const auto& qnL  = dot(Vi.velocity(),unitNormal);

  VectorDbl<nDim> velR;
  for(iDim = 0; iDim < nDim; ++iDim) velL(iDim) = Vj.velocity(iDim);
  const auto& rhoR = Vj.density();
  const auto& TR = Vj.temperature();
  const auto& pR = Vj.pressure();
  const auto& HR = Vj.enthalpy();
  const auto& aR = Vj.speedSound();
  const auto& q2R = squaredNorm<nDim>(Vj.velocity());
  const auto& qnR  = dot(Vj.velocity(),unitNormal);

  const auto& njk = unitNormal;

  MatrixDbl<nDim,nDim> grad_Vel;
  for (iDim = 0; iDim < nDim; ++iDim)
    for (jDim = 0; jDim < nDim; ++jDim)
//      grad_Vel(iDim,jDim) = diss(jDim)*(velR(iDim) - velL(iDim));
      grad_Vel(iDim,jDim) = grad(iDim+1,jDim);

  VectorDbl<nDim> grad_T, grad_p, grad_rho;
  for (iDim = 0; iDim < nDim; ++iDim) {
//    grad_T(iDim) = diss(iDim)*(TR - TL);
    grad_T(iDim) = grad(0,iDim);
//    grad_p(iDim) = diss(iDim)*(pR - pL);
    grad_p(iDim) = grad(nDim+1,iDim);
//    grad_rho(iDim) = diss(iDim)*(rhoR - rhoL);
    grad_rho(iDim) = grad(nDim+2,iDim);
  }

//  const Double Rt = p/rho;
//  for (iDim = 0; iDim < nDim; ++iDim)
//      grad_rho(iDim) = grad_p(iDim) / Rt - rho*grad_T(iDim)/T;

  const Double four_third = 4.0 / 3.0;
  const Double two_third = 2.0 / 3.0;
  const Double half = 0.5;
  const Double zero = 0.0;
  const Double one = 1.0;

  const Double& mu = V.laminarVisc() + V.eddyVisc();
  const auto& mu_l = V.laminarVisc();
  const auto& mu_t = V.eddyVisc();

  //! Viscous stresses (Stokes' hypothesis is assumed)
  //
  //!    tauxx =  mu*sxx
  //!    tauyy =  mu*sxy
  //!    tauzz =  mu*szz
//  const Double sxx = four_third*grad_u(ix) - two_third*grad_v(iy) - two_third*grad_w(iz);
//  const Double syy = four_third*grad_v(iy) - two_third*grad_u(ix) - two_third*grad_w(iz);
//  const Double szz = four_third*grad_w(iz) - two_third*grad_u(ix) - two_third*grad_v(iy);

  /** diagonal part **/
  VectorDbl<nDim> sDiag;
  for (iDim = 0; iDim < nDim; ++iDim) {
    sDiag(iDim) = 2.0 * grad_Vel(iDim,iDim);
    for (jDim = 0; jDim < nDim; ++jDim)
      sDiag(iDim) -= two_third * grad_Vel(jDim,jDim);
  }
  //
  //!    tauxy =  mu*sxy
  //!    tauxz =  mu*sxz
  //!    tauyz =  mu*syz
  //       sxy = grad_u(iy) + grad_v(ix);
  //       sxz = grad_u(iz) + grad_w(ix);
  //       syz = grad_v(iz) + grad_w(iy);

  /** symmetric part **/
  VectorDbl<nDim> sSymm;
  sSymm(0) = grad_Vel(0,1) + grad_Vel(1,0);
#if nDim == 3
    sSymm(1) = grad_Vel(0,2) + grad_Vel(2,0);
    sSymm(2) = grad_Vel(1,2) + grad_Vel(2,1);
#endif

  MatrixDbl<nDim,nDim> sTens;
  for (iDim = 0; iDim < nDim; ++iDim) sTens(iDim,iDim) = sDiag(iDim);
  for (iDim = 1; iDim < nDim; ++iDim) sTens(0,iDim) = sSymm(iDim-1);
  for (iDim = 2; iDim < nDim; ++iDim) sTens(1,iDim) = sSymm(iDim);

  for (iDim = 0; iDim < nDim; ++iDim)
    for (jDim = 0; jDim < iDim; ++jDim)
      sTens(iDim,jDim) = sTens(jDim,iDim);

  //
  //! Normal components
  //
  //!    taunx = mu*snx
  //!    tauny = mu*sny
  //!    taunz = mu*snz
  //       snx = sxx*njk(ix) + sxy*njk(iy) + sxz*njk(iz);
  //       sny = sxy*njk(ix) + syy*njk(iy) + syz*njk(iz);
  //       snz = sxz*njk(ix) + syz*njk(iy) + szz*njk(iz);

  VectorDbl<nDim> sn;
  for (iDim = 0; iDim < nDim; ++iDim)
    for (jDim = 0; jDim < nDim; ++jDim)
      sn(iDim) = sTens(iDim,jDim) * njk(jDim);
  //
  //! Viscous work
  //!    taunv = mu*snv = mu* ( snx*u + sny*v + snz*w )
  //       snv = snx*u + sny*v + snz*w

  Double snv = 0.0;
  for (iDim = 0; iDim < nDim; ++iDim)
    snv += sn(iDim) * vel(iDim);



  //!--------------------------------------------------------------------------------
  //! Part 1. Compute dFn_alpha/ducL
  //!
  //!  dFn_alpha/ducL = [ zero, -d(taunx)/ducL, -d(tauny)/ducL, -d(taunz)/ducL, d(-taunv+qn)/ducL]
  //!                 = (dFn_alpha/dwL)*(dwL/ducL)
  //!
  //!  where wL = [rhoL, uL, vL, wL, pL] - primitive variables.
  //!
  //!  So, we proceed as follows:
  //!
  //!  1.1 Compute the derivatives of the viscosity: dmu/drhoL, dmu/dpL
  //!  1.2 Compute d(taunx)/dwL, -d(tauny)/dwL, -d(taunz)/dwL, -d(taunv)/dwL
  //!  1.3 Compute d(qn)/dwL
  //!  1.4 Compute the flux Jacbian with respect to the primitive variables:
  //!      -> dFndwL
  //!  1.5 Compute the flux Jacbian with respect to the conservative variables:
  //!      -> dFnducL = dFndwL*(dWL/ducL)
  //
  //!----------------------------------------------------------------------
  //! 1.1 Compute the derivatives of the viscosity: dmu/drhoL, dmu/dpL
  //
  //! dmu/dTL = d( (one+C/T_inf)/(T+C/T_inf)*T**(three*half) * M_inf/Re_inf )/dTL

//    dmu_dTL = half*mu/(T+C/T_inf)*(one + three*(C/T_inf)/T) * ( half );
  const Double& dmul_dT = secondary.dmudT_rho();

  const Double dmu_dTL = 0.5 * dmul_dT;

  //! Note: dmu/dW = [dmu/drho, dmu/du, dmu/dv, dmu/dw, dmu/dp ] = [dmu/drho, 0, 0, 0, dmu/dp ]
  //!       So, we only need dmu/drho and dmu/dp.
  //
  //! dmu/drhoL = dmu/dTL*(dTL/drhoL); dT/drho = -gamma*p/rho**2 because T = p/(R*rho).

//    dmu_drhoL = dmu_dTL * (-gamma*pL/rhoL**2);
  const Double dmu_drhoL = dmu_dTL * (-TL/rhoL);

  //!   dmu/dpL = dmu/dTL*(dTL/dpL); dT/dp = gamma/rho because T = gamma*p/rho.

//      dmu_dpL = dmu_dTL * ( gamma/rhoL      );
      const Double dmu_dpL = dmu_dTL * (-TL/pL);

  //!----------------------------------------------------------------------
  //! 1.2 Compute d(taunx)/dwL, -d(tauny)/dwL, -d(taunz)/dwL, -d(taunv)/dwL
  //!
  //! Derivatives of the interface velocity gradients:
  //!
  //! Interface gradients are given by
  //!
  //!   grad_u = (Average LSQ gradient of u) + alpha/Lr*(uR - uL)*njk
  //!   grad_v = (Average LSQ gradient of v) + alpha/Lr*(vR - vL)*njk
  //!   grad_w = (Average LSQ gradient of w) + alpha/Lr*(wR - wL)*njk
  //!
  //! but for simplicity and compantness, we ignore the average LSQ gradients.
  //! Note: If the average LSQ gradients are computed on a compact stencil,
  //!       we may add their contributions.
  //
  //  dux_duL = alpha/Lr*(zero-one)*njk(ix)
  //  duy_duL = alpha/Lr*(zero-one)*njk(iy)
  //  duz_duL = alpha/Lr*(zero-one)*njk(iz)
  //
  //  dvx_dvL = alpha/Lr*(zero-one)*njk(ix)
  //  dvy_dvL = alpha/Lr*(zero-one)*njk(iy)
  //  dvz_dvL = alpha/Lr*(zero-one)*njk(iz)
  //
  //  dwx_dwL = alpha/Lr*(zero-one)*njk(ix)
  //  dwy_dwL = alpha/Lr*(zero-one)*njk(iy)
  //  dwz_dwL = alpha/Lr*(zero-one)*njk(iz)

      MatrixDbl<nDim,nDim> dGradVel_dVelL;
      for (iDim = 0; iDim < nDim; ++iDim)
        for (jDim = 0; jDim < nDim; ++jDim)
          dGradVel_dVelL(iDim,jDim) = - diss(jDim);

  //
  //
  //! Derivatives of the stress tensor.
  //! Note: many terms are zero, e.g., dux/dvL=0, dwz/dvL=0, etc.

      MatrixDbl<nDim> dsDiag_dVelL;
      for (iDim = 0; iDim < nDim; ++iDim) {
        for (jDim = 0; jDim < nDim; ++jDim)
          dsDiag_dVelL(iDim,jDim) = - two_third * dGradVel_dVelL(jDim,jDim);
        dsDiag_dVelL(iDim,iDim) = four_third * dGradVel_dVelL(iDim,iDim);
      }

  //
  //!  dsxy_duL = duy_duL + dvx_duL
  //!  dsxz_duL = duz_duL + dwx_duL
  //!  dsyz_duL = dvz_duL + dwy_duL
  //!  Ignoring terms that are zero, we obtain the following.
  VectorDbl<nDim> dsxy_dvelL, dsxz_dvelL, dsyz_dvelL;
  dsxy_dvelL(0) = dGradVel_dVelL(0,1);
  //
  //
  //!  dsxy_dvL = duy_dvL + dvx_dvL
  //!  dsxz_dvL = duz_dvL + dwx_dvL
  //!  dsyz_dvL = dvz_dvL + dwy_dvL
  //!  Ignoring terms that are zero, we obtain the following.
     dsxy_dvelL(1) = dGradVel_dVelL(1,0);
  //
  //
  //!  dsxy_dwL = duy_dwL + dvx_dwL
  //!  dsxz_dwL = duz_dwL + dwx_dwL
  //!  dsyz_dwL = dvz_dwL + dwy_dwL
  //!  Ignoring terms that are zero, we obtain the following.
#if nDim == 3
    dsxy_dvelL(2) = zero;

    dsxz_dvelL(0) = dGradVel_dVelL(0,2);
    dsxz_dvelL(1) = zero;
    dsxz_dvelL(2) = dGradVel_dVelL(2, 0);

    dsyz_dvelL(0) = zero;
    dsyz_dvelL(1) = dGradVel_dVelL(1,2);
    dsyz_dvelL(2) = dGradVel_dVelL(2, 1);
#endif


    MatrixDbl<nDim,nDim> dsn_dvelL;

    dsn_dvelL(0,0) = dsDiag_dVelL(0,0)*njk(0) + dsxy_dvelL(0)*njk(1);
    dsn_dvelL(1,0) = dsxy_dvelL(0)*njk(0) + dsDiag_dVelL(1,0)*njk(1);
#if nDim == 3
      dsn_dvelL(2, 0) = dsxz_dvelL(0) * njk(0) + dsyz_dvelL(0) * njk(1) + dsDiag_dVelL(2, 0) * njk(2);
      dsn_dvelL(0,0) += dsxz_dvelL(0)*njk(2);
      dsn_dvelL(1,0) += dsyz_dvelL(0)*njk(2);
#endif

    dsn_dvelL(0,1) = dsDiag_dVelL(0,1)*njk(0) + dsxy_dvelL(1)*njk(1);
    dsn_dvelL(1,1) = dsxy_dvelL(1)*njk(0) + dsDiag_dVelL(1,1)*njk(1);
#if nDim == 3
      dsn_dvelL(2, 1) = dsxz_dvelL(1) * njk(0) + dsyz_dvelL(1) * njk(1) + dsDiag_dVelL(2, 1) * njk(2);
      dsn_dvelL(0,1) += dsxz_dvelL(1)*njk(2);
      dsn_dvelL(1,1) += dsyz_dvelL(1)*njk(2);
#endif

#if nDim == 3
      dsn_dvelL(0,2) = dsDiag_dVelL(0,2)*njk(0) + dsxy_dvelL(2)*njk(1) + dsxz_dvelL(2)*njk(2);
      dsn_dvelL(1,2) = dsxy_dvelL(2)*njk(0) + dsDiag_dVelL(1,2)*njk(1) + dsyz_dvelL(2)*njk(2);
      dsn_dvelL(2,2) = dsxz_dvelL(2) * njk(0) + dsyz_dvelL(2) * njk(1) + dsDiag_dVelL(2, 2) * njk(2);
#endif
  //
  //! Derivatives of the viscous stresses.
  //! Note: mu does not depend on the velocity.
  //
    MatrixDbl<nDim> dtauDiag_dvelL;
    for (iDim = 0; iDim < nDim; ++iDim)
      for (jDim = 0; jDim < nDim; ++jDim)
        dtauDiag_dvelL(iDim,jDim) = mu*dsDiag_dVelL(iDim,jDim);


    VectorDbl<nDim> dtauxy_dvelL, dtauxz_dvelL, dtauyz_dvelL;
    for (iDim = 0; iDim < nDim; ++iDim)
      dtauxy_dvelL(iDim) = mu * dsxy_dvelL(iDim);

    if (nDim == 3) {
      for (iDim = 0; iDim < nDim; ++iDim)
        dtauxz_dvelL(iDim) = mu * dsxz_dvelL(iDim);

      for (iDim = 0; iDim < nDim; ++iDim)
        dtauyz_dvelL(iDim) = mu * dsyz_dvelL(iDim);
    }
  //
  //! Derivatives of the viscous stresses in the normal direction.
  //
  //! Derivatives of taunx
  const Double dtaunx_drhoL = dmu_drhoL * sn(0);

  MatrixDbl<nDim,nDim> dtaun_dvelL;
  dtaun_dvelL(0,0) = dtauDiag_dvelL(0,0)*njk(0) + dtauxy_dvelL(0)*njk(1);
  dtaun_dvelL(0,1) = dtauDiag_dvelL(0,1)*njk(0) + dtauxy_dvelL(1)*njk(1);
#if nDim == 3
    dtaun_dvelL(0,2) = dtauDiag_dvelL(0,2)*njk(0) + dtauxy_dvelL(2)*njk(1) + dtauxz_dvelL(2)*njk(2);
    dtaun_dvelL(0,0) += dtauxz_dvelL(0)*njk(2);
    dtaun_dvelL(0,1) += dtauxz_dvelL(1)*njk(2);
#endif

  const Double dtaunx_dpL   = dmu_dpL   * sn(0);
  //
  //! Derivatives of tauny
  const Double dtauny_drhoL = dmu_drhoL * sn(1);

    dtaun_dvelL(1,0) = dtauDiag_dvelL(1,0)*njk(1) + dtauxy_dvelL(0)*njk(0);
    dtaun_dvelL(1,1) = dtauDiag_dvelL(1,1)*njk(1) + dtauxy_dvelL(1)*njk(0);
#if nDim == 3
        dtaun_dvelL(1,2) = dtauDiag_dvelL(1,2)*njk(1) + dtauxy_dvelL(2)*njk(0) + dtauyz_dvelL(2)*njk(2);
        dtaun_dvelL(1,0) += dtauyz_dvelL(0)*njk(2);
        dtaun_dvelL(1,1) += dtauyz_dvelL(1)*njk(2);
#endif

  const Double dtauny_dpL   = dmu_dpL   * sn(1);
  //
  //! Derivatives of taunz
  Double dtaunz_drhoL, dtaunz_dpL;
#if nDim == 3
  dtaunz_drhoL = dmu_drhoL * sn(2);

  dtaun_dvelL(2,0) = dtauDiag_dvelL(2,0)*njk(2) + dtauxz_dvelL(0)*njk(0) + dtauyz_dvelL(0)*njk(1);
  dtaun_dvelL(2,1) = dtauDiag_dvelL(2,1)*njk(2) + dtauxz_dvelL(1)*njk(0) + dtauyz_dvelL(1)*njk(1);;
  dtaun_dvelL(2,2) = dtauDiag_dvelL(2,2)*njk(2) + dtauxz_dvelL(2)*njk(0) + dtauyz_dvelL(2)*njk(1);
  dtaunz_dpL   = dmu_dpL   * sn(2);
#endif
  //
  //! Derivatives of taunv = mu* ( snx*u + sny*v + snz*w )
  VectorDbl<nDim> dtaunv_dvelL;
  const Double dtaunv_drhoL = dmu_drhoL * snv;
  for (iDim = 0; iDim < nDim; ++iDim) {
    dtaunv_dvelL(iDim) = mu * sn(iDim) * half;
    for (jDim = 0; jDim < nDim; ++jDim)
      dtaunv_dvelL(iDim) += mu * dsn_dvelL(jDim,iDim)*vel(jDim);
  }
  const Double dtaunv_dpL   = dmu_dpL   * snv;
  //
  //!----------------------------------------------------------------------

  //!----------------------------------------------------------------------
  //! 1.3 Compute d(qn)/dwL
  //
  //! First compute the derivatives of the interface gradients:
  //!
  //!       rho = half*(rhoR + rhoL)
  //!        a2 = gamma*half*(pR + pL)/rho
  //!  grad_rho = (Average LSQ gradient of rho) + alpha/Lr*(rhoR - rhoL)*njk
  //!  grad_p   = (Average LSQ gradient of p)   + alpha/Lr*(  pR - pL  )*njk
  //!  grad_T   = ( gamma*grad_p - a2*grad_rho) /rho
  //!
  //!  Again, we ignore the average LSQ gradients.
  //
  //  drhox_drhoL = alpha/Lr*(zero-one)*njk(ix) ! Average LSQ gradient ignored
  //  drhoy_drhoL = alpha/Lr*(zero-one)*njk(iy) ! Average LSQ gradient ignored
  //  drhoz_drhoL = alpha/Lr*(zero-one)*njk(iz) ! Average LSQ gradient ignored
  VectorDbl<nDim> drho_drhoL;
  for (iDim = 0; iDim < nDim; ++iDim)
    drho_drhoL(iDim) = - diss(iDim);

  VectorDbl<nDim> dp_dpL;
  for (iDim = 0; iDim < nDim; ++iDim)
    dp_dpL(iDim) = - diss(iDim);

  VectorDbl<nDim> dT_drhoL;
  const Double Rrho = p / T;
  for (iDim = 0; iDim < nDim; ++iDim)
    dT_drhoL(iDim) = half * ((-grad_p(iDim) + 2.0*p*grad_rho(iDim)/rho)/(Rrho*rho) - p*drho_drhoL(iDim)/(Rrho*rho));
//    dT_drhoL(iDim) = (- gamma*grad_p(iDim)/(rho*rho)
//                 + 2.0*gamma*p/(rho*rho*rho)*grad_rho(iDim) ) * half
//                 - a*a/rho*drho_drhoL(iDim);

  VectorDbl<nDim> dT_dpL;
  for (iDim = 0; iDim < nDim; ++iDim)
    dT_dpL(iDim) = half * (dp_dpL(iDim)/(Rrho) - (grad_rho(iDim) / (Rrho*rho)));
//    dT_dpL(iDim) = gamma*dp_dpL(iDim)/rho - ( gamma/(rho*rho)*grad_rho(iDim) ) * half;

  //
  //! Derivatives of the heat fluxes:
  //!       q = -cp*(mu_lam/Pr_lam + mu_tur/Pr_tur)*grad_T
  //!       qx = - mu*grad_T(ix)/(Prandtl*(gamma-one))
  //!       qy = - mu*grad_T(iy)/(Prandtl*(gamma-one))
  //!       qz = - mu*grad_T(iz)/(Prandtl*(gamma-one))

  //

  VectorDbl<nDim> dq_drhoL, dq_dpL;
  for (iDim = 0; iDim < nDim; ++iDim)
    dq_drhoL(iDim) = - cp*(mu_l/Pr_l + mu_t/Pr_t)*dT_drhoL(iDim) -cp*dmu_drhoL/Pr_l * grad_T(iDim);
  for (iDim = 0; iDim < nDim; ++iDim)
    dq_dpL(iDim) = - cp*(mu_l/Pr_l + mu_t/Pr_t)*dT_dpL(iDim) -cp*dmu_dpL/Pr_l * grad_T(iDim);

  //
  //! Derivatives of the normal heat flux: qn = qx*njk(ix) + qy*njk(iy) + qz*njk(iz)
  //
  Double dqn_drhoL = 0.0, dqn_dpL = 0.0;
  for (iDim = 0; iDim < nDim; ++iDim) {
    dqn_drhoL += dq_drhoL(iDim)*njk(iDim);
    dqn_dpL += dq_dpL(iDim)*njk(iDim);
  }

  //
  //!----------------------------------------------------------------------
  //! 1.4 Compute the flux Jacbian with respect to the primitive variables:
  //!     -> dFndwL
  //
  MatrixDbl<nDim+2,nDim+2> dFndwL;
  //! 1st row: d(zero)/dW=zero (no viscous term in the continuity equation)
  dFndwL(0,0) = zero;
  for (iDim = 0; iDim < nDim; ++iDim)
    dFndwL(0,iDim+1) = zero;
  dFndwL(0,nDim+1) = zero;
  //
  //! 2nd row: d(-taunx)/dW; taunx = tauxx*nx + tauxy*ny + tauxz*nz
    dFndwL(1,0) = - dtaunx_drhoL;
    for (iDim = 0; iDim < nDim; ++iDim)
      dFndwL(1,iDim+1) = - dtaun_dvelL(0,iDim);
    dFndwL(1,nDim+1) = - dtaunx_dpL  ;
  //
  //! 3rd row: d(-tauny)/dW; tauny = tauyx*nx + tauyy*ny + tauyz*nz
    dFndwL(2,0) = - dtauny_drhoL;
    for (iDim = 0; iDim < nDim; ++iDim)
      dFndwL(2,iDim+1) = - dtaun_dvelL(1,iDim);
    dFndwL(2,nDim+1) = - dtauny_dpL;
  //
  //! 4th row: d(-taunz)/dW; taunz = tauzx*nx + tauzy*ny + tauzz*nz
#if nDim == 3
    dFndwL(3, 0) = -dtaunz_drhoL;
    for (iDim = 0; iDim < nDim; ++iDim)
      dFndwL(3, iDim + 1) = -dtaun_dvelL(2,iDim);
    dFndwL(3, nDim + 1) = -dtaunz_dpL;
#endif
  //
  //! 5th row: d(-taunv + qn)/dW; taunv = taunx*u + tauny*v + taunz*nz; qn=qx*nx+qy*ny+qz*nz
    dFndwL(nDim+1,0) = - dtaunv_drhoL + dqn_drhoL;
    for (iDim = 0; iDim < nDim; ++iDim)
      dFndwL(nDim+1,iDim+1) = - dtaunv_dvelL(iDim);
    dFndwL(nDim+1,nDim+1) = - dtaunv_dpL   + dqn_dpL;
  //
  //!----------------------------------------------------------------------
  //! 1.5 Compute the flux Jacbian with respect to the conservative variables:
  //!     -> dFnducL = dFndwL*(dWL/ducL)
  //!
  //!   dW/duc(1,:) = (     1,     0,     0,     0, 0 )
  //!   dW/duc(2,:) = (-u/rho, 1/rho,     0,     0, 0 )
  //!   dW/duc(3,:) = (-v/rho,     0, 1/rho,     0, 0 )
  //!   dW/duc(4,:) = (-w/rho,     0,     0, 1/rho, 0 )
  //!   dW/duc(5,:) = ( q^2/2,    -u,    -v,    -w, 1 )*(gamma-1)
  //
  //! 1st row
    jac_i(0,0) = zero;
    for (iDim = 0; iDim < nDim; ++iDim)
      jac_i(0,iDim+1) = zero;
    jac_i(0,nDim+1) = zero;
  //
  //! 2nd row

    for (jDim = 0; jDim < nDim; ++jDim) {
      jac_i(jDim+1, 0) = dFndwL(jDim+1, 0);
      for (iDim = 0; iDim < nDim; ++iDim)
        jac_i(jDim+1, 0) += dFndwL(jDim+1, iDim + 1) * (-velL(iDim) / rhoL);
      jac_i(jDim+1, 0) += dFndwL(jDim+1, nDim + 1) * half * (gamma - one) * q2L;
      for (iDim = 0; iDim < nDim; ++iDim)
        jac_i(jDim+1, iDim + 1) = dFndwL(jDim+1, iDim + 1) * (one / rhoL) + dFndwL(jDim+1, nDim + 1) * (-(gamma - one) * velL(iDim));
      jac_i(jDim+1, nDim + 1) = dFndwL(jDim+1, nDim + 1) * (gamma - one);
    }
  //
  //! 5th row
    jac_i(nDim+1,0) = dFndwL(nDim+1,0);
    for (iDim = 0; iDim < nDim; ++iDim)
      jac_i(nDim+1,0) += dFndwL(nDim+1,iDim+1) * (-velL(iDim)/rhoL);
    jac_i(nDim+1,0) += dFndwL(nDim+1,nDim+1) * half*(gamma-one)*q2L;
    for (iDim = 0; iDim < nDim; ++iDim)
      jac_i(nDim+1,iDim+1) = dFndwL(nDim+1,iDim+1) * (one/rhoL) + dFndwL(nDim+1,nDim+1) * (-(gamma-one)*velL(iDim));
    jac_i(nDim+1,nDim+1) = dFndwL(nDim+1,nDim+1) *   (gamma-one);


//!--------------------------------------------------------------------------------
//! Part 2. Compute dFn_alpha/ducR
//!
//!  dFn_alpha/ducR = [ zero, -d(taunx)/ducR, -d(tauny)/ducR, -d(taunz)/ducR, d(-taunv+qn)/ducR]
//!                 = (dFn_alpha/dwR)*(dwR/ducR)
//!
//!  where wR = [rhoR, uR, vR, wR, pR] - primitive variables.
//!
//!  So, we proceed as follows:
//!
//!  1.1 Compute the derivatives of the viscosity: dmu/drhoR, dmu/dpR
//!  1.2 Compute d(taunx)/dwR, -d(tauny)/dwR, -d(taunz)/dwR, -d(taunv)/dwR
//!  1.3 Compute d(qn)/dwR
//!  1.4 Compute the flux Jacbian with respect to the primitive variables:
//!      -> dFndwR
//!  1.5 Compute the flux Jacbian with respect to the conservative variables:
//!      -> dFnducR = dFndwR*(dWR/ducR)
//
//!----------------------------------------------------------------------
//! 1.1 Compute the derivatives of the viscosity: dmu/drhoR, dmu/dpR
//
//! dmu/dTR = d( (one+C/T_inf)/(T+C/T_inf)*T**(three*half) * M_inf/Re_inf )/dTR

//    dmu_dTR = half*mu/(T+C/T_inf)*(one + three*(C/T_inf)/T) * ( half );
const Double dmu_dTR = half * dmul_dT;

//! Note: dmu/dW = [dmu/drho, dmu/du, dmu/dv, dmu/dw, dmu/dp ] = [dmu/drho, 0, 0, 0, dmu/dp ]
//!       So, we only need dmu/drho and dmu/dp.
//
//! dmu/drhoR = dmu/dTR*(dTR/drhoR); dT/drho = -gamma*p/rho**2 because T = gamma*p/rho.

//    dmu_drhoR = dmu_dTR * (-gamma*pR/rhoR**2);
    const Double dmu_drhoR = dmu_dTR * (-TR/rhoR);

    //!   dmu/dpL = dmu/dTL*(dTL/dpL); dT/dp = gamma/rho because T = gamma*p/rho.

//      dmu_dpL = dmu_dTL * ( gamma/rhoL      );
    const Double dmu_dpR = dmu_dTR * (-TR/pR);

//!----------------------------------------------------------------------
//! 2.2 Compute d(taunx)/dwR, -d(tauny)/dwR, -d(taunz)/dwR, -d(taunv)/dwR
//!
//! Derivatives of the interface velocity gradients:
//!
//! Interface gradients are given by
//!
//!   grad_u = (Average RSQ gradient of u) + alpha/Rr*(uR - uR)*njk
//!   grad_v = (Average RSQ gradient of v) + alpha/Rr*(vR - vR)*njk
//!   grad_w = (Average RSQ gradient of w) + alpha/Rr*(wR - wR)*njk
//!
//! but for simplicity and compantness, we ignore the average RSQ gradients.
//! Note: If the average RSQ gradients are computed on a compact stencil,
//!       we may add their contributions.


MatrixDbl<nDim,nDim> dGradVel_dVelR;
for (iDim = 0; iDim < nDim; ++iDim)
  for (jDim = 0; jDim < nDim; ++jDim)
    dGradVel_dVelR(iDim,jDim) = diss(jDim);

//
//
//! Derivatives of the stress tensor.
//! Note: many terms are zero, e.g., dux/dvR=0, dwz/dvR=0, etc.

MatrixDbl<nDim> dsDiag_dVelR;
for (iDim = 0; iDim < nDim; ++iDim) {
  for (jDim = 0; jDim < nDim; ++jDim)
    dsDiag_dVelR(iDim,jDim) = - two_third * dGradVel_dVelR(jDim,jDim);
  dsDiag_dVelR(iDim,iDim) = four_third * dGradVel_dVelR(iDim,iDim);
}

//
//!  dsxy_duR = duy_duR + dvx_duR
//!  dsxz_duR = duz_duR + dwx_duR
//!  dsyz_duR = dvz_duR + dwy_duR
//!  Ignoring terms that are zero, we obtain the following.
VectorDbl<nDim> dsxy_dvelR, dsxz_dvelR, dsyz_dvelR;
dsxy_dvelR(0) = dGradVel_dVelR(0,1);
//
//
//!  dsxy_dvR = duy_dvR + dvx_dvR
//!  dsxz_dvR = duz_dvR + dwx_dvR
//!  dsyz_dvR = dvz_dvR + dwy_dvR
//!  Ignoring terms that are zero, we obtain the following.
dsxy_dvelR(1) = dGradVel_dVelR(1,0);
//
//
//!  dsxy_dwR = duy_dwR + dvx_dwR
//!  dsxz_dwR = duz_dwR + dwx_dwR
//!  dsyz_dwR = dvz_dwR + dwy_dwR
//!  Ignoring terms that are zero, we obtain the following.
#if nDim == 3
  dsxy_dvelR(2) = zero;

  dsxz_dvelR(0) = dGradVel_dVelR(0,2);
  dsxz_dvelR(1) = zero;
  dsxz_dvelR(2) = dGradVel_dVelR(2, 0);

  dsyz_dvelR(0) = zero;
  dsyz_dvelR(1) = dGradVel_dVelR(1,2);
  dsyz_dvelR(2) = dGradVel_dVelR(2, 1);

#endif


MatrixDbl<nDim,nDim> dsn_dvelR;

dsn_dvelR(0,0) = dsDiag_dVelR(0,0)*njk(0) + dsxy_dvelR(0)*njk(1);
dsn_dvelR(1,0) = dsxy_dvelR(0)*njk(0) + dsDiag_dVelR(1,0)*njk(1);
#if nDim == 3
  dsn_dvelR(2, 0) = dsxz_dvelR(0) * njk(0) + dsyz_dvelR(0) * njk(1) + dsDiag_dVelR(2, 0) * njk(2);
  dsn_dvelR(0,0) += dsxz_dvelR(0)*njk(2);
  dsn_dvelR(1,0) += dsyz_dvelR(0)*njk(2);
#endif

dsn_dvelR(0,1) = dsDiag_dVelR(0,1)*njk(0) + dsxy_dvelR(1)*njk(1);
dsn_dvelR(1,1) = dsxy_dvelR(1)*njk(0) + dsDiag_dVelR(1,1)*njk(1);
#if nDim == 3
  dsn_dvelR(2, 1) = dsxz_dvelR(1) * njk(0) + dsyz_dvelR(1) * njk(1) + dsDiag_dVelR(2, 1) * njk(2);
  dsn_dvelR(0,1) += dsxz_dvelR(1)*njk(2);
  dsn_dvelR(1,1) += dsyz_dvelR(1)*njk(2);
#endif

#if nDim == 3
  dsn_dvelR(0,2) = dsDiag_dVelR(0,2)*njk(0) + dsxy_dvelR(2)*njk(1);
  dsn_dvelR(1,2) = dsxy_dvelR(2)*njk(0) + dsDiag_dVelR(1,2)*njk(1);
  dsn_dvelR(2, 2) = dsxz_dvelR(2) * njk(0) + dsyz_dvelR(2) * njk(1) + dsDiag_dVelR(2, 2) * njk(2);
  dsn_dvelR(0,2) += dsxz_dvelR(2)*njk(2);
  dsn_dvelR(1,2) += dsyz_dvelR(2)*njk(2);
#endif
//
//! Derivatives of the viscous stresses.
//! Note: mu does not depend on the velocity.
//
MatrixDbl<nDim> dtauDiag_dvelR;
for (iDim = 0; iDim < nDim; ++iDim)
  for (jDim = 0; jDim < nDim; ++jDim)
    dtauDiag_dvelR(iDim,jDim) = mu*dsDiag_dVelR(iDim,jDim);


VectorDbl<nDim> dtauxy_dvelR, dtauxz_dvelR, dtauyz_dvelR;
for (iDim = 0; iDim < nDim; ++iDim)
  dtauxy_dvelR(iDim) = mu * dsxy_dvelR(iDim);

#if nDim == 3
  for (iDim = 0; iDim < nDim; ++iDim)
    dtauxz_dvelR(iDim) = mu * dsxz_dvelR(iDim);

  for (iDim = 0; iDim < nDim; ++iDim)
    dtauyz_dvelR(iDim) = mu * dsyz_dvelR(iDim);
#endif
//
//! Derivatives of the viscous stresses in the normal direction.
//
//! Derivatives of taunx
const Double dtaunx_drhoR = dmu_drhoR * sn(0);

MatrixDbl<nDim,nDim> dtaun_dvelR;
dtaun_dvelR(0,0) = dtauDiag_dvelR(0,0)*njk(0) + dtauxy_dvelR(0)*njk(1);
dtaun_dvelR(0,1) = dtauDiag_dvelR(0,1)*njk(0) + dtauxy_dvelR(1)*njk(1);
#if nDim == 3
  dtaun_dvelR(0,2) = dtauDiag_dvelR(0,2)*njk(0) + dtauxy_dvelR(2)*njk(1) + dtauxz_dvelR(2)*njk(2);
  dtaun_dvelR(0,0) += dtauxz_dvelR(0)*njk(2);
  dtaun_dvelR(0,1) += dtauxz_dvelR(1)*njk(2);
#endif

const Double dtaunx_dpR   = dmu_dpR   * sn(0);
//
//! Derivatives of tauny
const Double dtauny_drhoR = dmu_drhoR * sn(1);

dtaun_dvelR(1,0) = dtauDiag_dvelR(1,0)*njk(1) + dtauxy_dvelR(0)*njk(0);
dtaun_dvelR(1,1) = dtauDiag_dvelR(1,1)*njk(1) + dtauxy_dvelR(1)*njk(0);
#if nDim == 3
  dtaun_dvelR(1,2) = dtauDiag_dvelR(1,2)*njk(1) + dtauxy_dvelR(2)*njk(0) + dtauyz_dvelR(2)*njk(2);
  dtaun_dvelR(1,0) += dtauyz_dvelR(0)*njk(2);
  dtaun_dvelR(1,1) += dtauyz_dvelR(1)*njk(2);
#endif

const Double dtauny_dpR   = dmu_dpR   * sn(1);
//
//! Derivatives of taunz
Double dtaunz_drhoR, dtaunz_dpR;
#if nDim == 3
  dtaunz_drhoR = dmu_drhoR * sn(2);

  dtaun_dvelR(2,0) = dtauDiag_dvelR(2,0)*njk(2) + dtauxz_dvelR(0)*njk(0);
  dtaun_dvelR(2,1) = dtauDiag_dvelR(2,1)*njk(2) + dtauxz_dvelR(1)*njk(0);
  dtaun_dvelR(2,2) = dtauDiag_dvelR(2,2)*njk(2) + dtauxz_dvelR(2)*njk(0) + dtauyz_dvelR(2)*njk(1);
  dtaun_dvelR(2,0) += dtauyz_dvelR(0)*njk(1);
  dtaun_dvelR(2,1) += dtauyz_dvelR(1)*njk(1);
  dtaunz_dpR   = dmu_dpR   * sn(2);
#endif
//
//! Derivatives of taunv = mu* ( snx*u + sny*v + snz*w )
VectorDbl<nDim> dtaunv_dvelR;
const Double dtaunv_drhoR = dmu_drhoR * snv;
for (iDim = 0; iDim < nDim; ++iDim) {
  dtaunv_dvelR(iDim) = mu * sn(iDim) * half;
  for (jDim = 0; jDim < nDim; ++jDim)
    dtaunv_dvelR(iDim) += mu * dsn_dvelR(jDim,iDim)*vel(jDim);
}
const Double dtaunv_dpR   = dmu_dpR   * snv;
//
//!----------------------------------------------------------------------

//!----------------------------------------------------------------------
//! 1.3 Compute d(qn)/dwR
//
//! First compute the derivatives of the interface gradients:
//!
//!       rho = half*(rhoR + rhoR)
//!        a2 = gamma*half*(pR + pR)/rho
//!  grad_rho = (Average RSQ gradient of rho) + alpha/Rr*(rhoR - rhoR)*njk
//!  grad_p   = (Average RSQ gradient of p)   + alpha/Rr*(  pR - pR  )*njk
//!  grad_T   = ( gamma*grad_p - a2*grad_rho) /rho
//!
//!  Again, we ignore the average RSQ gradients.
//
//  drhox_drhoR = alpha/Rr*(zero-one)*njk(ix) ! Average RSQ gradient ignored
//  drhoy_drhoR = alpha/Rr*(zero-one)*njk(iy) ! Average RSQ gradient ignored
//  drhoz_drhoR = alpha/Rr*(zero-one)*njk(iz) ! Average RSQ gradient ignored
VectorDbl<nDim> drho_drhoR;
for (iDim = 0; iDim < nDim; ++iDim)
  drho_drhoR(iDim) = diss(iDim);

VectorDbl<nDim> dp_dpR;
for (iDim = 0; iDim < nDim; ++iDim)
  dp_dpR(iDim) = diss(iDim);

VectorDbl<nDim> dT_drhoR;
for (iDim = 0; iDim < nDim; ++iDim)
  dT_drhoR(iDim) = half * ((-grad_p(iDim) + 2.0*p*grad_rho(iDim)/rho)/(Rrho*rho) - p*drho_drhoR(iDim)/(Rrho*rho));
//    dT_drhoL(iDim) = (- gamma*grad_p(iDim)/(rho*rho)
//                 + 2.0*gamma*p/(rho*rho*rho)*grad_rho(iDim) ) * half
//                 - a*a/rho*drho_drhoL(iDim);

VectorDbl<nDim> dT_dpR;
for (iDim = 0; iDim < nDim; ++iDim)
  dT_dpR(iDim) = half * (dp_dpR(iDim)/(Rrho) - (grad_rho(iDim) / (Rrho*rho)));

//
//! Derivatives of the heat fluxes:
//!       q = -cp*(mu_lam/Pr_lam + mu_tur/Pr_tur)*grad_T
//!       qx = - mu*grad_T(ix)/(Prandtl*(gamma-one))
//!       qy = - mu*grad_T(iy)/(Prandtl*(gamma-one))
//!       qz = - mu*grad_T(iz)/(Prandtl*(gamma-one))

//

VectorDbl<nDim> dq_drhoR, dq_dpR;
for (iDim = 0; iDim < nDim; ++iDim)
  dq_drhoR(iDim) = - cp*(mu_l/Pr_l + mu_t/Pr_t)*dT_drhoR(iDim) -cp*dmu_drhoR/Pr_l * grad_T(iDim);
for (iDim = 0; iDim < nDim; ++iDim)
  dq_dpR(iDim) = - cp*(mu_l/Pr_l + mu_t/Pr_t)*dT_dpR(iDim) -cp*dmu_dpR/Pr_l * grad_T(iDim);

//
//! Derivatives of the normal heat flux: qn = qx*njk(ix) + qy*njk(iy) + qz*njk(iz)
//
Double dqn_drhoR = 0.0, dqn_dpR = 0.0;
for (iDim = 0; iDim < nDim; ++iDim) {
  dqn_drhoR += dq_drhoR(iDim)*njk(iDim);
  dqn_dpR += dq_dpR(iDim)*njk(iDim);
}

//
//!----------------------------------------------------------------------
//! 1.4 Compute the flux Jacbian with respect to the primitive variables:
//!     -> dFndwR
//
MatrixDbl<nDim+2,nDim+2> dFndwR;
//! 1st row: d(zero)/dW=zero (no viscous term in the continuity equation)
dFndwR(0,0) = zero;
dFndwR(0,1) = zero;
dFndwR(0,2) = zero;
dFndwR(0,3) = zero;
dFndwR(0,nDim+1) = zero;
//
//! 2nd row: d(-taunx)/dW; taunx = tauxx*nx + tauxy*ny + tauxz*nz
dFndwR(1,0) = - dtaunx_drhoR;
for (iDim = 0; iDim < nDim; ++iDim)
  dFndwR(1,iDim+1) = - dtaun_dvelR(0,iDim);
dFndwR(1,nDim+1) = - dtaunx_dpR  ;
//
//! 3rd row: d(-tauny)/dW; tauny = tauyx*nx + tauyy*ny + tauyz*nz
dFndwR(2,0) = - dtauny_drhoR;
for (iDim = 0; iDim < nDim; ++iDim)
  dFndwR(2,iDim+1) = - dtaun_dvelR(1,iDim);
dFndwR(2,nDim+1) = - dtauny_dpR;
//
//! 4th row: d(-taunz)/dW; taunz = tauzx*nx + tauzy*ny + tauzz*nz
if (nDim == 3) {
  dFndwR(3, 0) = -dtaunz_drhoR;
  for (iDim = 0; iDim < nDim; ++iDim)
    dFndwR(3, iDim + 1) = -dtaun_dvelR(2,iDim);
  dFndwR(3, nDim + 1) = -dtaunz_dpR;
}
//
//! 5th row: d(-taunv + qn)/dW; taunv = taunx*u + tauny*v + taunz*nz; qn=qx*nx+qy*ny+qz*nz
dFndwR(nDim+1,0) = - dtaunv_drhoR + dqn_drhoR;
for (iDim = 0; iDim < nDim; ++iDim)
  dFndwR(nDim+1,iDim+1) = - dtaunv_dvelR(iDim);
dFndwR(nDim+1,nDim+1) = - dtaunv_dpR   + dqn_dpR;
//
//!----------------------------------------------------------------------
//! 1.5 Compute the flux Jacbian with respect to the conservative variables:
//!     -> dFnducR = dFndwR*(dWR/ducR)
//!
//!   dW/duc(1,:) = (     1,     0,     0,     0, 0 )
//!   dW/duc(2,:) = (-u/rho, 1/rho,     0,     0, 0 )
//!   dW/duc(3,:) = (-v/rho,     0, 1/rho,     0, 0 )
//!   dW/duc(4,:) = (-w/rho,     0,     0, 1/rho, 0 )
//!   dW/duc(5,:) = ( q^2/2,    -u,    -v,    -w, 1 )*(gamma-1)
//
//! 1st row
jac_j(0,0) = zero;
for (iDim = 0; iDim < nDim; ++iDim)
  jac_j(0,iDim+1) = zero;
jac_j(0,nDim+1) = zero;
//
//! 2nd row

for (jDim = 0; jDim < nDim; ++jDim) {
  jac_j(jDim+1, 0) = dFndwR(jDim+1, 0);
  for (iDim = 0; iDim < nDim; ++iDim)
    jac_j(jDim+1, 0) += dFndwR(jDim+1, iDim + 1) * (-velR(iDim) / rhoR);
  jac_j(jDim+1, 0) += dFndwR(jDim+1, nDim + 1) * half * (gamma - one) * q2R;
  for (iDim = 0; iDim < nDim; ++iDim)
    jac_j(jDim+1, iDim + 1) = dFndwR(jDim+1, iDim + 1) * (one / rhoR) + dFndwR(jDim+1, nDim + 1) * (-(gamma - one) * velR(iDim));
  jac_j(jDim+1, nDim + 1) = dFndwR(jDim+1, nDim + 1) * (gamma - one);
}
//
//! 5th row
jac_j(nDim+1,0) = dFndwR(nDim+1,0);
for (iDim = 0; iDim < nDim; ++iDim)
  jac_j(nDim+1,0) += dFndwR(nDim+1,iDim+1) * (-velR(iDim)/rhoR);
jac_j(nDim+1,0) += dFndwR(nDim+1,nDim+1) * half*(gamma-one)*q2R;
for (iDim = 0; iDim < nDim; ++iDim)
  jac_j(nDim+1,iDim+1) = dFndwR(nDim+1,iDim+1) * (one/rhoR) + dFndwR(nDim+1,nDim+1) * (-(gamma-one)*velR(iDim));
jac_j(nDim+1,nDim+1) = dFndwR(nDim+1,nDim+1) *   (gamma-one);

}
