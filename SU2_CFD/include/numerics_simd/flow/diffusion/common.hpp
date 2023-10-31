/*!
 * \file common.hpp
 * \brief Helper functions for viscous methods.
 * \author P. Gomes, C. Pederson, A. Bueno, F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

/*!
 * \brief Correct average gradient with the directional derivative to avoid decoupling.
 */
template<size_t nVar, size_t nDim, class PrimitiveType>
FORCEINLINE void correctGradient(const PrimitiveType& V,
                                 const VectorDbl<nDim>& vector_ij,
                                 Double dist2_ij,
                                 MatrixDbl<nVar,nDim>& avgGrad) {
  for (size_t iVar = 0; iVar < nVar; ++iVar) {
    Double corr = (dot(avgGrad[iVar],vector_ij) - V.j.all(iVar) + V.i.all(iVar)) / dist2_ij;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      avgGrad(iVar,iDim) -= corr * vector_ij(iDim);
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
