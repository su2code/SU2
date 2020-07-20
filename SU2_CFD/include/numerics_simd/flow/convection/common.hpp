/*!
 * \file common.hpp
 * \brief Common convection-related methods.
 * \author P. Gomes, F. Palacios, T. Economon
 * \version 7.0.6 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../../variables/CNSVariable.hpp"

/*!
 * \brief Unlimited reconstruction.
 */
template<size_t nVar, size_t nDim, class Field_t, class Gradient_t>
FORCEINLINE VectorDbl<nVar> musclUnlimited(Int iPoint,
                                           const VectorDbl<nDim>& vector_ij,
                                           Double direction,
                                           const Field_t& field,
                                           const Gradient_t& gradient) {
  auto vars = gatherVariables<nVar>(iPoint, field);
  auto grad = gatherVariables<nVar,nDim>(iPoint, gradient);
  for (size_t iVar = 0; iVar < nVar; ++iVar) {
    vars(iVar) += direction * dot(grad[iVar], vector_ij);
  }
  return vars;
}

/*!
 * \brief Limited reconstruction with point-based limiter.
 */
template<size_t nVar, size_t nDim, class Field_t, class Gradient_t>
FORCEINLINE VectorDbl<nVar> musclPointLimited(Int iPoint,
                                              const VectorDbl<nDim>& vector_ij,
                                              Double direction,
                                              const Field_t& field,
                                              const Field_t& limiter,
                                              const Gradient_t& gradient) {
  auto vars = gatherVariables<nVar>(iPoint, field);
  auto lim = gatherVariables<nVar>(iPoint, limiter);
  auto grad = gatherVariables<nVar,nDim>(iPoint, gradient);
  for (size_t iVar = 0; iVar < nVar; ++iVar) {
    vars(iVar) += lim(iVar) * direction * dot(grad[iVar], vector_ij);
  }
  return vars;
}

/*!
 * \brief Limited reconstruction with edge-based limiter.
 */
template<size_t nDim, size_t nVar, class Field_t, class Gradient_t>
FORCEINLINE void musclEdgeLimited(Int iPoint,
                                  Int jPoint,
                                  const VectorDbl<nDim>& vector_ij,
                                  const Field_t& field,
                                  const Gradient_t& gradient,
                                  VectorDbl<nVar>& vars_i,
                                  VectorDbl<nVar>& vars_j) {
  vars_i = gatherVariables<nVar>(iPoint, field);
  vars_j = gatherVariables<nVar>(jPoint, field);

  auto grad_i = gatherVariables<nVar,nDim>(iPoint, gradient);
  auto grad_j = gatherVariables<nVar,nDim>(jPoint, gradient);

  for (size_t iVar = 0; iVar < nVar; ++iVar) {
    auto proj_i = dot(grad_i[iVar], vector_ij);
    auto proj_j = dot(grad_j[iVar], vector_ij);
    auto delta_ij = vars_j(iVar) - vars_i(iVar);
    auto delta_ij_2 = pow(delta_ij,2);
    /// TODO: Customize the limiter function (template functor).
    auto lim_i = (delta_ij_2 + 2*proj_i*delta_ij) / (4*pow(proj_i,2) + delta_ij_2 + EPS);
    auto lim_j = (delta_ij_2 + 2*proj_j*delta_ij) / (4*pow(proj_j,2) + delta_ij_2 + EPS);
    vars_i(iVar) += lim_i * proj_i;
    vars_j(iVar) -= lim_j * proj_j;
  }
}

/*!
 * \brief Retrieve primitive variables for points i/j, reconstructing them if needed.
 * \param[in] iPoint, jPoint - Nodes of the edge.
 * \param[in] muscl - If true, reconstruct, else simply fetch.
 * \param[in] vector_ij - Distance vector from i to j.
 * \param[in] solution - Entire solution container (a derived CVariable).
 * \return Pair of primitive variables.
 */
template<class PrimitiveType, size_t nDim, class VariableType>
FORCEINLINE CPair<PrimitiveType> reconstructPrimitives(Int iPoint, Int jPoint, bool muscl,
                                                       ENUM_LIMITER limiterType,
                                                       const VectorDbl<nDim>& vector_ij,
                                                       const VariableType& solution) {
  CPair<PrimitiveType> V;
  constexpr size_t nVar = PrimitiveType::nVar;

  const auto& primitives = solution.GetPrimitive();
  const auto& gradients = solution.GetGradient_Reconstruction();
  const auto& limiters = solution.GetLimiter_Primitive();

  if (!muscl) {
    V.i.all = gatherVariables<nVar>(iPoint, primitives);
    V.j.all = gatherVariables<nVar>(jPoint, primitives);
  }
  else {
    switch (limiterType) {
    case NO_LIMITER:
      V.i.all = musclUnlimited<nVar>(iPoint, vector_ij, 1, primitives, gradients);
      V.j.all = musclUnlimited<nVar>(jPoint, vector_ij,-1, primitives, gradients);
      break;
    case VAN_ALBADA_EDGE:
      musclEdgeLimited(iPoint, jPoint, vector_ij, primitives, gradients, V.i.all, V.j.all);
      break;
    default:
      V.i.all = musclPointLimited<nVar>(iPoint, vector_ij, 1, primitives, limiters, gradients);
      V.j.all = musclPointLimited<nVar>(jPoint, vector_ij,-1, primitives, limiters, gradients);
      break;
    }
    /// TODO: Extra reconstruction checks needed.
  }
  return V;
}

/*!
 * \brief Compute and return the P tensor (compressible flow, ideal gas).
 */
template<size_t nDim>
FORCEINLINE MatrixDbl<nDim+2> pMatrix(Double gamma, Double density, const VectorDbl<nDim>& velocity,
                                      Double projVel, Double speedSound, const VectorDbl<nDim>& normal) {
  MatrixDbl<nDim+2> pMat;
  const Double vel2 = 0.5*squaredNorm(velocity);

  if (nDim == 2) {
    pMat(0,0) = 1.0;
    pMat(0,1) = 0.0;

    pMat(1,0) = velocity(0);
    pMat(1,1) = density*normal(1);

    pMat(2,0) = velocity(1);
    pMat(2,1) = -density*normal(0);

    pMat(3,0) = vel2;
    pMat(3,1) = density*(velocity(0)*normal(1) - velocity(1)*normal(0));
  }
  else {
    pMat(0,0) = normal(0);
    pMat(0,1) = normal(1);
    pMat(0,2) = normal(2);

    pMat(1,0) = velocity(0)*normal(0);
    pMat(1,1) = velocity(0)*normal(1) - density*normal(2);
    pMat(1,2) = velocity(0)*normal(2) + density*normal(1);

    pMat(2,0) = velocity(1)*normal(0) + density*normal(2);
    pMat(2,1) = velocity(1)*normal(1);
    pMat(2,2) = velocity(1)*normal(2) - density*normal(0);

    pMat(3,0) = velocity(2)*normal(0) - density*normal(1);
    pMat(3,1) = velocity(2)*normal(1) + density*normal(0);
    pMat(3,2) = velocity(2)*normal(2);

    pMat(4,0) = vel2*normal(0) + density*(velocity(1)*normal(2) - velocity(2)*normal(1));
    pMat(4,1) = vel2*normal(1) - density*(velocity(0)*normal(2) - velocity(2)*normal(0));
    pMat(4,2) = vel2*normal(2) + density*(velocity(0)*normal(1) - velocity(1)*normal(0));
  }

  /*--- Last two columns. ---*/

  const Double rhoOn2 = 0.5*density;
  const Double rhoOnTwoC = rhoOn2 / speedSound;
  pMat(0,nDim) = rhoOnTwoC;
  pMat(0,nDim+1) = rhoOnTwoC;

  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    pMat(iDim+1,nDim) = rhoOnTwoC * velocity(iDim) + rhoOn2 * normal(iDim);
    pMat(iDim+1,nDim+1) = rhoOnTwoC * velocity(iDim) - rhoOn2 * normal(iDim);
  }

  pMat(nDim+1,nDim) = rhoOnTwoC * vel2 + rhoOn2 * (projVel + speedSound/(gamma-1));
  pMat(nDim+1,nDim+1) = rhoOnTwoC * vel2 - rhoOn2 * (projVel - speedSound/(gamma-1));

  return pMat;
}

/*!
 * \brief Compute and return the inverse P tensor (compressible flow, ideal gas).
 */
template<size_t nDim>
FORCEINLINE MatrixDbl<nDim+2> pMatrixInv(Double gamma, Double density, const VectorDbl<nDim>& velocity,
                                         Double projVel, Double speedSound, const VectorDbl<nDim>& normal) {
  MatrixDbl<nDim+2> pMatInv;

  const Double c2 = pow(speedSound,2);
  const Double vel2 = 0.5*squaredNorm(velocity);
  const Double oneOnRho = 1 / density;

  if (nDim == 2) {
    Double tmp = (gamma-1)/c2;
    pMatInv(0,0) = 1.0 - tmp*vel2;
    pMatInv(0,1) = tmp*velocity(0);
    pMatInv(0,2) = tmp*velocity(1);
    pMatInv(0,3) = -tmp;

    pMatInv(1,0) = (normal(0)*velocity(1)-normal(1)*velocity(0))*oneOnRho;
    pMatInv(1,1) = normal(1)*oneOnRho;
    pMatInv(1,2) = -normal(0)*oneOnRho;
    pMatInv(1,3) = 0.0;
  }
  else {
    Double tmp = (gamma-1)/c2 * normal(0);
    pMatInv(0,0) = normal(0) - tmp*vel2 - (normal(2)*velocity(1)-normal(1)*velocity(2))*oneOnRho;
    pMatInv(0,1) = tmp*velocity(0);
    pMatInv(0,2) = tmp*velocity(1) + normal(2)*oneOnRho;
    pMatInv(0,3) = tmp*velocity(2) - normal(1)*oneOnRho;
    pMatInv(0,4) = -tmp;

    tmp = (gamma-1)/c2 * normal(1);
    pMatInv(1,0) = normal(1) - tmp*vel2 + (normal(2)*velocity(0)-normal(0)*velocity(2))*oneOnRho;
    pMatInv(1,1) = tmp*velocity(0) - normal(2)*oneOnRho;
    pMatInv(1,2) = tmp*velocity(1);
    pMatInv(1,3) = tmp*velocity(2) + normal(0)*oneOnRho;
    pMatInv(1,4) = -tmp;

    tmp = (gamma-1)/c2 * normal(2);
    pMatInv(2,0) = normal(2) - tmp*vel2 - (normal(1)*velocity(0)-normal(0)*velocity(1))*oneOnRho;
    pMatInv(2,1) = tmp*velocity(0) + normal(1)*oneOnRho;
    pMatInv(2,2) = tmp*velocity(1) - normal(0)*oneOnRho;
    pMatInv(2,3) = tmp*velocity(2);
    pMatInv(2,4) = -tmp;
  }

  /*--- Last two rows. ---*/

  const Double gamma_minus_1_on_rho_times_c = (gamma-1) / (density*speedSound);

  for (size_t iVar = nDim; iVar < nDim+2; ++iVar) {
    Double sign = (iVar==nDim)? 1 : -1;
    pMatInv(iVar,0) = -sign*projVel*oneOnRho + gamma_minus_1_on_rho_times_c * vel2;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      pMatInv(iVar,iDim+1) = sign*normal(iDim)*oneOnRho - gamma_minus_1_on_rho_times_c * velocity(iDim);
    }
    pMatInv(iVar,nDim+1) = gamma_minus_1_on_rho_times_c;
  }

  return pMatInv;
}

/*!
 * \brief Convective projected (onto normal) flux (compressible flow).
 */
template<size_t nDim, size_t N>
FORCEINLINE VectorDbl<nDim+2> inviscidProjFlux(const CCompressiblePrimitives<nDim,N>& V,
                                               const CCompressibleConservatives<nDim>& U,
                                               const VectorDbl<nDim>& normal) {
  Double mdot = dot(U.momentum(), normal);
  VectorDbl<nDim+2> flux;
  flux(0) = mdot;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    flux(iDim+1) = mdot*V.velocity(iDim) + normal(iDim)*V.pressure();
  }
  flux(nDim+1) = mdot*V.enthalpy();
  return flux;
}

/*!
 * \brief Jacobian of the convective flux (compressible flow, ideal gas).
 */
template<size_t nDim, class RandomAccessIterator>
FORCEINLINE MatrixDbl<nDim+2> inviscidProjJac(Double gamma, RandomAccessIterator velocity,
                                              Double energy, const VectorDbl<nDim>& normal,
                                              Double scale) {
  MatrixDbl<nDim+2> jac;

  Double projVel = dot(velocity, normal);
  Double gamma_m_1 = gamma-1;
  Double phi = 0.5*gamma_m_1*squaredNorm<nDim>(velocity);
  Double a1 = gamma*energy - phi;

  jac(0,0) = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    jac(0,iDim+1) = scale * normal(iDim);
  }
  jac(0,nDim+1) = 0.0;

  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    jac(iDim+1,0) = scale * (normal(iDim)*phi - velocity[iDim]*projVel);
    for (size_t jDim = 0; jDim < nDim; ++jDim) {
      jac(iDim+1,jDim+1) = scale * (normal(jDim)*velocity[iDim] - gamma_m_1*normal(iDim)*velocity[jDim]);
    }
    jac(iDim+1,iDim+1) += scale * projVel;
    jac(iDim+1,nDim+1) = scale * gamma_m_1 * normal(iDim);
  }

  jac(nDim+1,0) = scale * projVel * (phi-a1);
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    jac(nDim+1,iDim+1) = scale * (normal(iDim)*a1 - gamma_m_1*velocity[iDim]*projVel);
  }
  jac(nDim+1,nDim+1) = scale * gamma * projVel;

  return jac;
}

/*!
 * \brief (Low) Dissipation coefficient for Roe schemes.
 */
template<class VariableType>
FORCEINLINE Double roeDissipation(Int iPoint,
                                  Int jPoint,
                                  ENUM_ROELOWDISS type,
                                  const VariableType& solution) {
  if (type == NO_ROELOWDISS) {
    return 1.0;
  }

  const auto& sol = static_cast<const CNSVariable&>(solution);
  const auto& sensor = sol.GetSensor();
  const auto& dissip = sol.GetRoe_Dissipation();

  const Double si = gatherVariables(iPoint, sensor);
  const Double sj = gatherVariables(jPoint, sensor);
  const Double avgSensor = 0.5 * (si + sj);

  const Double di = gatherVariables(iPoint, dissip);
  const Double dj = gatherVariables(jPoint, dissip);
  const Double avgDissip = 0.5 * (di + dj);

  /*--- A minimum level of upwinding is used to enhance stability. ---*/
  constexpr passivedouble minDissip = 0.05;

  switch (type) {
    case FD:
    case FD_DUCROS: {
      Double d = max(minDissip, 1.0 - avgDissip);

      if (type == FD_DUCROS) {
        /*--- See Jonhsen et al. JCP 229 (2010) pag. 1234 ---*/
        d = max(d, 0.05 + 0.95*(avgSensor > 0.65));
      }
      return d;
    }
    case NTS:
      return max(minDissip, avgDissip);

    case NTS_DUCROS:
      /*--- See Xiao et al. INT J HEAT FLUID FL 51 (2015) pag. 141
       * https://doi.org/10.1016/j.ijheatfluidflow.2014.10.007 ---*/
      return max(minDissip, avgSensor+avgDissip - avgSensor*avgDissip);

    default:
      assert(false);
      return 1.0;
  }
}

/*!
 * \brief Update of a flux Jacobian due to a scalar dissipation term.
 */
template<class VariableType, size_t nVar>
FORCEINLINE void scalarDissipationJacobian(const VariableType& V,
                                           Double gamma,
                                           Double dissipConst,
                                           MatrixDbl<nVar>& jac) {
  /*--- Diagonal entries. ---*/
  for (size_t iVar = 0; iVar < nVar-1; ++iVar) {
    jac(iVar,iVar) += dissipConst;
  }
  jac(nVar-1,nVar-1) += dissipConst * gamma;

  /*--- N-1 columns of last row. ---*/
  dissipConst *= gamma-1.0;
  for (size_t iDim = 0; iDim < VariableType::nDim; ++iDim) {
    jac(nVar-1,iDim+1) -= dissipConst * V.velocity(iDim);
    jac(nVar-1,0) += dissipConst * pow(V.velocity(iDim), 2);
  }
}
