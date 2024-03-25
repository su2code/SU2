/*!
 * \file roe.hpp
 * \brief Roe-family of convective schemes.
 * \author P. Gomes, A. Bueno, F. Palacios
 * \version 7.5.1 "Blackbird"
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
#include "common.hpp"
#include "../../../variables/CEulerVariable.hpp"
#include "../../../../../Common/include/geometry/CGeometry.hpp"

/*!
 * \class CRoeBase
 * \ingroup ConvDiscr
 * \brief Base class for Roe schemes, derived classes implement
 * the dissipation term in a const "finalizeFlux" method.
 * A base class implementing "viscousTerms" is accepted as template parameter.
 * Similarly to derived, that method should update the flux and Jacobians, but
 * whereas "finalizeFlux" is passed data prepared by CRoeBase, "viscousTerms"
 * takes the same input arguments as "ComputeFlux", i.e. it can fetch more
 * data from CVariable. Derived is meant to implement small details,
 * Base is meant to do heavy lifting.
 */
template<class Derived, class Base>
class CRoeBase : public Base {
protected:
  using Base::nDim;
  static constexpr size_t nVar = CCompressibleConservatives<nDim>::nVar;
  static constexpr size_t nPrimVarGrad = nDim+4;
  static constexpr size_t nPrimVar = Max(Base::nPrimVar, nPrimVarGrad);

  const su2double kappa;
  const su2double gamma;
  const su2double entropyFix;
  const bool finestGrid;
  const bool dynamicGrid;
  const bool muscl;
  const bool correction;
  const LIMITER typeLimiter;

  /*!
   * \brief Constructor, store some constants and forward args to base.
   */
  template<class... Ts>
  CRoeBase(const CConfig& config, unsigned iMesh, Ts&... args) : Base(config, iMesh, args...),
    kappa(config.GetRoe_Kappa()),
    gamma(config.GetGamma()),
    entropyFix(config.GetEntropyFix_Coeff()),
    finestGrid(iMesh == MESH_0),
    dynamicGrid(config.GetDynamic_Grid()),
    muscl(finestGrid && config.GetMUSCL_Flow()),
    correction(muscl && config.GetFluxCorrection()),
    typeLimiter(config.GetKind_SlopeLimit_Flow()) {
  }

public:
  /*!
   * \brief Implementation of the base Roe flux.
   */
  void ComputeFlux(Int iEdge,
                   const CConfig& config,
                   const CGeometry& geometry,
                   const CVariable& solution_,
                   UpdateType updateType,
                   Double updateMask,
                   CSysVector<su2double>& vector,
                   SparseMatrixType& matrix) const final {

    /*--- Start preaccumulation, inputs are registered
     *    automatically in "gatherVariables". ---*/
    AD::StartPreacc();

    const bool implicit = (config.GetKind_TimeIntScheme() == EULER_IMPLICIT);
    const auto& solution = static_cast<const CEulerVariable&>(solution_);

    const auto iPoint = geometry.edges->GetNode(iEdge,0);
    const auto jPoint = geometry.edges->GetNode(iEdge,1);

    /*--- Geometric properties. ---*/

    const auto vector_ij = distanceVector<nDim>(iPoint, jPoint, geometry.nodes->GetCoord());

    const auto normal = gatherVariables<nDim>(iEdge, geometry.edges->GetNormal());
    const auto area = norm(normal);
    VectorDbl<nDim> unitNormal;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      unitNormal(iDim) = normal(iDim) / area;
    }

    /*--- Reconstructed primitives. ---*/

    CPair<CCompressiblePrimitives<nDim,nPrimVar> > V1st;
    V1st.i.all = gatherVariables<nPrimVar>(iPoint, solution.GetPrimitive());
    V1st.j.all = gatherVariables<nPrimVar>(jPoint, solution.GetPrimitive());

    auto V = reconstructPrimitives<CCompressiblePrimitives<nDim,nPrimVarGrad> >(
                 iEdge, iPoint, jPoint, muscl, typeLimiter, V1st, vector_ij, solution);

    /*--- Compute conservative variables. ---*/

    CPair<CCompressibleConservatives<nDim> > U;
    U.i = compressibleConservatives(V.i);
    U.j = compressibleConservatives(V.j);

    /*--- Roe averaged variables. ---*/

    auto roeAvg = roeAveragedVariables(gamma, V, unitNormal);

    /*--- P tensor. ---*/

    auto pMat = pMatrix(gamma, roeAvg.density, roeAvg.velocity,
                        roeAvg.projVel, roeAvg.speedSound, unitNormal);

    /*--- Grid motion. ---*/

    Double projGridVel = 0.0, projVel = roeAvg.projVel;
    if (dynamicGrid) {
      const auto& gridVel = geometry.nodes->GetGridVel();
      projGridVel = 0.5*(dot(gatherVariables<nDim>(iPoint,gridVel), unitNormal)+
                         dot(gatherVariables<nDim>(jPoint,gridVel), unitNormal));
      projVel -= projGridVel;
    }

    /*--- Convective eigenvalues. ---*/

    VectorDbl<nVar> lambda;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      lambda(iDim) = projVel;
    }
    lambda(nDim) = projVel + roeAvg.speedSound;
    lambda(nDim+1) = projVel - roeAvg.speedSound;

    /*--- Apply Mavriplis' entropy correction to eigenvalues. ---*/

    Double maxLambda = abs(projVel) + roeAvg.speedSound;

    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      lambda(iVar) = fmax(abs(lambda(iVar)), entropyFix*maxLambda);
    }

    /*--- Inviscid fluxes and Jacobians. ---*/

    auto flux_i = inviscidProjFlux(V.i, U.i, normal);
    auto flux_j = inviscidProjFlux(V.j, U.j, normal);

    VectorDbl<nVar> flux;
    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      flux(iVar) = 0.5 * (flux_i(iVar) + flux_j(iVar));
    }

    MatrixDbl<nVar> jac_i, jac_j;
    if (implicit) {
      jac_i = inviscidProjJac(gamma, V.i.velocity(), U.i.energy(), normal, kappa);
      jac_j = inviscidProjJac(gamma, V.j.velocity(), U.j.energy(), normal, kappa);
    }

    /*--- Correct for grid motion. ---*/

    if (dynamicGrid) {
      for (size_t iVar = 0; iVar < nVar; ++iVar) {
        Double dFdU = projGridVel * area * 0.5;
        flux(iVar) -= dFdU * (U.i.all(iVar) + U.j.all(iVar));

        if (implicit) {
          jac_i(iVar,iVar) -= dFdU;
          jac_j(iVar,iVar) -= dFdU;
        }
      }
    }

    /*--- Finalize in derived class (static polymorphism). ---*/

    const auto derived = static_cast<const Derived*>(this);

    derived->finalizeFlux(flux, jac_i, jac_j, implicit, area, unitNormal, V,
                          U, roeAvg, lambda, pMat, iPoint, jPoint, solution);

    /*--- Correct for second order on arbitrary grids ---*/
    if (correction) {
      const auto corrX = gatherVariables<nDim>(iEdge, geometry.edges->GetCorrection_X());
      const auto corrY = gatherVariables<nDim>(iEdge, geometry.edges->GetCorrection_Y());
      const auto corrZ = gatherVariables<nDim>(iEdge, geometry.edges->GetCorrection_Z());
      CorrectFlux(iPoint, jPoint, solution, V, flux, corrX, corrY, corrZ, gamma);
//      reconstructPrimitives<CCompressiblePrimitives<nDim,nPrimVarGrad> >(iEdge, iPoint, jPoint, muscl, typeLimiter, V1st, vector_ij, solution);
    }

    /*--- Add the contributions from the base class (static decorator). ---*/

    Base::viscousTerms(iEdge, iPoint, jPoint, V1st, solution_, vector_ij, geometry,
                       config, area, unitNormal, implicit, flux, jac_i, jac_j);

    /*--- Stop preaccumulation. ---*/

    stopPreacc(flux);

    /*--- Update the vector and system matrix. ---*/

    updateLinearSystem(iEdge, iPoint, jPoint, implicit, updateType,
                       updateMask, flux, jac_i, jac_j, vector, matrix);
  }
};

/*!
 * \class CRoeScheme
 * \ingroup ConvDiscr
 * \brief Classical Roe scheme.
 */
template<class Decorator>
class CRoeScheme : public CRoeBase<CRoeScheme<Decorator>,Decorator> {
private:
  using Base = CRoeBase<CRoeScheme<Decorator>,Decorator>;
  using Base::nDim;
  using Base::nVar;
  using Base::gamma;
  using Base::kappa;
  const ENUM_ROELOWDISS typeDissip;

public:
  /*!
   * \brief Constructor, store some constants and forward to base.
   */
  template<class... Ts>
  CRoeScheme(const CConfig& config, Ts&... args) : Base(config, args...),
    typeDissip(static_cast<ENUM_ROELOWDISS>(config.GetKind_RoeLowDiss())) {
  }

  /*!
   * \brief Updates flux and Jacobians with standard Roe dissipation.
   * \note "Ts" is here just in case other schemes in the family need extra args.
   */
  template<class PrimVarType, class ConsVarType, class... Ts>
  FORCEINLINE void finalizeFlux(VectorDbl<nVar>& flux,
                                MatrixDbl<nVar>& jac_i,
                                MatrixDbl<nVar>& jac_j,
                                bool implicit,
                                Double area,
                                const VectorDbl<nDim>& unitNormal,
                                const CPair<PrimVarType>& V,
                                const CPair<ConsVarType>& U,
                                const CRoeVariables<nDim>& roeAvg,
                                const VectorDbl<nVar>& lambda,
                                const MatrixDbl<nVar>& pMat,
                                Int iPoint,
                                Int jPoint,
                                const CEulerVariable& solution,
                                Ts&...) const {
    /*--- Inverse P tensor. ---*/

    auto pMatInv = pMatrixInv(gamma, roeAvg.density, roeAvg.velocity,
                              roeAvg.projVel, roeAvg.speedSound, unitNormal);

    /*--- Diference between conservative variables at jPoint and iPoint. ---*/

    VectorDbl<nVar> deltaU;
    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      deltaU(iVar) = U.j.all(iVar) - U.i.all(iVar);
    }

    /*--- Dissipation terms. ---*/

    Double dissipation = roeDissipation(iPoint, jPoint, typeDissip, solution);

    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      for (size_t jVar = 0; jVar < nVar; ++jVar) {
        /*--- Compute |projModJacTensor| = P x |Lambda| x P^-1. ---*/

        Double projModJacTensor = 0.0;
        for (size_t kVar = 0; kVar < nVar; ++kVar) {
          projModJacTensor += pMat(iVar,kVar) * lambda(kVar) * pMatInv(kVar,jVar);
        }

        Double dDdU = projModJacTensor * (1-kappa) * area * dissipation;

        /*--- Update flux and Jacobians. ---*/

        flux(iVar) -= dDdU * deltaU(jVar);

        if(implicit) {
          jac_i(iVar,jVar) += dDdU;
          jac_j(iVar,jVar) -= dDdU;
        }
      }
    }
  }
};

template<class Derived, class Base>
class CRoeNewBase : public Base {
 protected:
  using Base::nDim;
  static constexpr size_t nVar = CCompressibleConservatives<nDim>::nVar;
  static constexpr size_t nPrimVarGrad = nDim+4;
  static constexpr size_t nPrimVar = Max(Base::nPrimVar, nPrimVarGrad);

  const su2double kappa;
  const su2double gamma;
  const su2double entropyFix;
  const bool finestGrid;
  const bool dynamicGrid;
  const bool muscl;
  const bool correction;
  const LIMITER typeLimiter;

  /*!
   * \brief Constructor, store some constants and forward args to base.
   */
  template<class... Ts>
  CRoeNewBase(const CConfig& config, unsigned iMesh, Ts&... args) : Base(config, iMesh, args...),
                                                                    kappa(config.GetRoe_Kappa()),
                                                                    gamma(config.GetGamma()),
                                                                    entropyFix(config.GetEntropyFix_Coeff()),
                                                                    finestGrid(iMesh == MESH_0),
                                                                    dynamicGrid(config.GetDynamic_Grid()),
                                                                    muscl(finestGrid && config.GetMUSCL_Flow()),
                                                                    correction(muscl && config.GetFluxCorrection()),
                                                                    typeLimiter(config.GetKind_SlopeLimit_Flow()) {
  }

 public:

  /*!
   * \brief Implementation of the base Roe flux.
   */
  void ComputeFlux(Int iEdge,
                   const CConfig& config,
                   const CGeometry& geometry,
                   const CVariable& solution_,
                   UpdateType updateType,
                   Double updateMask,
                   CSysVector<su2double>& vector,
                   SparseMatrixType& matrix) const final {

    /*--- Start preaccumulation, inputs are registered
     *    automatically in "gatherVariables". ---*/
    AD::StartPreacc();

    const bool implicit = (config.GetKind_TimeIntScheme() == EULER_IMPLICIT);
    const auto& solution = static_cast<const CEulerVariable&>(solution_);

    const auto iPoint = geometry.edges->GetNode(iEdge,0);
    const auto jPoint = geometry.edges->GetNode(iEdge,1);

    /*--- Geometric properties. ---*/

    const auto vector_ij = distanceVector<nDim>(iPoint, jPoint, geometry.nodes->GetCoord());

    const auto normal = gatherVariables<nDim>(iEdge, geometry.edges->GetNormal());
    const auto area = norm(normal);
    VectorDbl<nDim> unitNormal;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      unitNormal(iDim) = normal(iDim) / area;
    }

    /*--- Reconstructed primitives. ---*/

    CPair<CCompressiblePrimitives<nDim,nPrimVar> > V1st;
    V1st.i.all = gatherVariables<nPrimVar>(iPoint, solution.GetPrimitive());
    V1st.j.all = gatherVariables<nPrimVar>(jPoint, solution.GetPrimitive());

    auto V = reconstructPrimitives<CCompressiblePrimitives<nDim,nPrimVarGrad> >(
        iEdge, iPoint, jPoint, muscl, typeLimiter, V1st, vector_ij, solution);

    /*--- Compute conservative variables. ---*/

    CPair<CCompressibleConservatives<nDim> > U;
    U.i = compressibleConservatives(V.i);
    U.j = compressibleConservatives(V.j);

    /*--- Roe averaged variables. ---*/

    auto roeAvg = roeAveragedVariables(gamma, V, unitNormal);

    /*--- Grid motion. ---*/

    Double projGridVel = 0.0, projVel = roeAvg.projVel;
    if (dynamicGrid) {
      const auto& gridVel = geometry.nodes->GetGridVel();
      projGridVel = 0.5*(dot(gatherVariables<nDim>(iPoint,gridVel), unitNormal)+
                           dot(gatherVariables<nDim>(jPoint,gridVel), unitNormal));
      projVel -= projGridVel;
    }

    const auto& vel = roeAvg.velocity;
    const auto& rho = roeAvg.density;
    const auto& H = roeAvg.enthalpy;
    const auto& qn = projVel;
    const auto& a = roeAvg.speedSound;
    const auto& a2 = a*a;
    const auto& q2 = squaredNorm<nDim>(roeAvg.velocity);

    VectorDbl<nDim> velL;
    for (int iDim = 0; iDim < nDim; ++iDim) velL(iDim) = V.i.velocity(iDim);
    const auto& rhoL = V.i.density();
    const auto& pL = V.i.pressure();
    const auto& HL = V.i.enthalpy();
//    const auto& aL = V.i.speedSound();
    const auto& aL = sqrt(gamma*pL/rhoL);
    const auto& q2L = squaredNorm<nDim>(V.i.velocity());
    const auto& qnL  = dot(V.i.velocity(),unitNormal);


    VectorDbl<nDim> velR;
    for (int iDim = 0; iDim < nDim; ++iDim) velR(iDim) = V.j.velocity(iDim);
    const auto& rhoR = V.j.density();
    const auto& pR = V.j.pressure();
    const auto& HR = V.j.enthalpy();
//    const auto& aR = V.j.speedSound();
    const auto& aR = sqrt(gamma*pR/rhoR);
    const auto& q2R = squaredNorm<nDim>(V.j.velocity());
    const auto& qnR  = dot(V.j.velocity(),unitNormal);

    const auto& n = unitNormal;

    //!Wave Strengths
    const Double drho = rhoR - rhoL; //!Density difference
    const Double dp   = pR - pL;   //!Pressure difference
    const Double dqn  = qnR - qnL;  //!Normal velocity difference

    VectorDbl<4> LdU;
    LdU(0) = (dp - rho*a*dqn )/(2*a2); //!Left-moving acoustic wave strength
    LdU(1) =  drho - dp/(a2);            //!Entropy wave strength
    LdU(2) = (dp + rho*a*dqn )/(2*a2); //!Right-moving acoustic wave strength
    LdU(3) = rho;                         //!Shear wave strength (not really, just a factor)

    //!Absolute values of the wave Speeds

    VectorDbl<4> ws;
    ws(0) = abs(qn-a); //!Left-moving acoustic wave
    ws(1) = abs(qn);   //!Entropy wave
    ws(2) = abs(qn+a); //!Right-moving acoustic wave
    ws(3) = abs(qn);   //!Shear waves

    auto ws_orig = ws;

    //!Harten's entropy fix JCP(1983), 49, pp357-393 is applied to avoid vanishing
    //!wave speeds by making a parabolic fit near ws = 0. It is an entropy fix for
    //!nonlinear waves (avoids expnsion shock), but applied here as an eigenvalue
    //!limiting, for the pusrpose of avoiding 0 eigenvalues (wave speeds).

    //!Nonlinear fields
    VectorDbl<4> dws;
    const Double elimc_nonlinear = 0.25;
    const Double elimc_linear = 0.05;

    dws(0) = elimc_nonlinear *a;
    dws(2) = elimc_nonlinear *a;
    dws(1) = elimc_linear *a;
    dws(3) = elimc_linear *a;

    for(int k = 0; k < 4; ++k) {
      const Double fix = ws(k) < dws(k);
      ws(k) = fix * (0.5 * (ws(k) * ws(k) / dws(k) + dws(k))) + (1-fix) * ws(k);
    }

    //!Right Eigenvectors
    //!Note: Two shear wave components are combined into 1, so that tangent vectors
    //!      are not required. And that's why there are only 4 vectors here.
    //!      See "I do like CFD, VOL.1" about how tangent vectors are eliminated.

    MatrixDbl<nDim+2> R;
    //! Left-moving acoustic wave
    R(0,0) = 1;
    for (int iDim = 0; iDim < nDim; ++iDim)
        R(iDim+1,0) = vel(iDim) - a*n(iDim);
    R(nDim+1,0) = H - a*qn;

    //! Entropy wave
    R(0,1) = 1;
    for (int iDim = 0; iDim < nDim; ++iDim)
        R(iDim+1,1) = vel(iDim);
    R(nDim+1,1) = 0.5*q2;

    //! Right-moving acoustic wave
    R(0,2) = 1;
    for (int iDim = 0; iDim < nDim; ++iDim)
        R(iDim+1,2) = vel(iDim) + a*n(iDim);
    R(nDim+1,2) = H + a*qn;

    //! Two shear wave components combined into 1 (wave strength incorporated).
    VectorDbl<nDim> dVel;
    for (int iDim = 0; iDim < nDim; ++iDim)
        dVel(iDim) = velR(iDim) - velL(iDim);

    const Double velDotdVel = dot(vel,dVel);

    R(0,3) = 0;
      for (int iDim = 0; iDim < nDim; ++iDim)
    R(iDim+1,3) = dVel(iDim) - dqn*n(iDim);

    R(nDim+1,3) = velDotdVel - qn*dqn;

    //! We are now ready to compute the Jacobians for the Roe flux:
    //!
    //! Roe flux function -> Fn_Roe = 0.5*[Fn(ucL)+Fn(ucR)] - 0.5*sum_{k=1,4}|lambda_k|*(LdU)_k*r_k

    //!--------------------------------------------------------------------------------
    //! Part 1. Compute dFn_Roe/ducL
    //!
    //!  dFn_Roe/ducL =   d(0.5*Fn(ucL))/duL
    //!                 - 0.5*sum_{k=1,4} [d(|lambda_k|)/ducL]*(LdU)_k*r_k
    //!                 - 0.5*sum_{k=1,4} |lambda_k|*[d(LdU)_k/ducL]*r_k
    //!                 - 0.5*sum_{k=1,4} |lambda_k|*(LdU)_k*[dr_k/ducL]
    //!
    //!  So, we proceed as follows:
    //!
    //!  1.1 Compute                d(0.5*Fn(ucL))/duL
    //!  1.2 Compute various deriavives that will be used in the following steps.
    //!  1.3 Add the second term, - 0.5*sum_{k=1,4} [d(|lambda_k|)/ducL]*(LdU)_k*r_k
    //!  1.4 Add the  third term, - 0.5*sum_{k=1,4} |lambda_k|*[d(LdU)_k/ducL]*r_k
    //!  1.5 Add the fourth term, - 0.5*sum_{k=1,4} |lambda_k|*(LdU)_k*[dr_k/ducL]
    //!

    //!--------------------------------------
    //! 1.1 Compute "d(0.5*Fn(ucL))/ducL"
    //!
    //!     (See "I Do Like CFD, VOL.1", page 55, for the analytical Jacobian, dFn(u)/du)

    MatrixDbl<nDim+2> dFnducL;

    //!  1st column
    dFnducL(0,0) = 0;
    for (int iDim = 0; iDim < nDim; ++iDim)
        dFnducL(iDim+1,0) = 0.5*(gamma-1)*q2L*n(iDim)  - velL(iDim)*qnL;
    dFnducL(nDim+1,0) = 0.5*(gamma-1)*q2L*qnL - HL*qnL;

    //!  2nd column
    for (int jDim = 0; jDim < nDim; ++jDim) {
        dFnducL(0,jDim+1) =     n(jDim);
        dFnducL(nDim+1,jDim+1) =  HL*n(jDim) - (gamma-1)*velL(jDim)*qnL;
        for (int iDim = 0; iDim < nDim; ++iDim)
            dFnducL(iDim+1,jDim+1) =  velL(iDim)*n(jDim) - (gamma-1)*velL(jDim)*n(iDim);
        dFnducL(jDim+1,jDim+1) += qnL;
    }

    //!  5th column
    dFnducL(0,nDim+1) =  0;
    for (int iDim = 0; iDim < nDim; ++iDim)
        dFnducL(iDim+1,nDim+1) = (gamma-1)*n(iDim);
    dFnducL(nDim+1,nDim+1) =  gamma*qnL;

    //!  Factor 1/2
    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) *= 0.5;

    //!--------------------------------------
    //! 1.2 Compute various deriavives that will be used in the following steps.
    //!     (See "I Do Like CFD, VOL.1" for details)

    //! dqn/ducL

    VectorDbl<nVar> dqn_ducL;
    dqn_ducL(0) = -0.5*(qnL+qn) / (rhoL+rho);
    for (int iDim = 0; iDim < nDim; ++iDim)
        dqn_ducL(iDim+1) =            n(iDim)  / (rhoL+rho);
    dqn_ducL(nDim+1) =          0;

    //! d(|qn|)/ducL

    VectorDbl<nVar> dabs_qn_ducL;
    Double mask = qn < 0;
    Double sigQN = mask * -1 + (1-mask) * 1;
    dabs_qn_ducL(0) = -0.5*sigQN*(qnL+qn) / (rhoL+rho);
    for (int iDim = 0; iDim < nDim; ++iDim)
        dabs_qn_ducL(iDim+1) =       sigQN*     n(iDim)  / (rhoL+rho);
    dabs_qn_ducL(nDim+1) =  0;

    //! da/ducL
    VectorDbl<nVar> da_ducL;
    const Double velDotVelL = dot(vel,velL);
    da_ducL(0) =  0.5*(gamma-1)/a*( 0.5*( velDotVelL + q2 )
                                          +  0.5*(HL-H) - aL*aL/(gamma-1) + 0.5*(gamma-2)*q2L )/(rhoL+rho);
    for (int iDim = 0; iDim < nDim; ++iDim)
        da_ducL(iDim+1) = -0.5*(gamma-1)*(vel(iDim)+(gamma-1)*velL(iDim))/a  / (rhoL+rho);
    da_ducL(nDim+1) =  0.5*gamma*(gamma-1)/a               / (rhoL+rho);

    //! drho/ducL
    VectorDbl<nVar> drho_ducL;
    drho_ducL(0) = 0.5*rho/rhoL;
    for (int iDim = 0; iDim < nDim; ++iDim)
        drho_ducL(iDim+1) = 0;
    drho_ducL(nDim+1) = 0;

    // dVel/ducL
    MatrixDbl<nDim,nVar> dvel_ducL;

    for (int iDim = 0; iDim < nDim; ++iDim) {
        dvel_ducL(iDim,0) = -0.5*(velL(iDim)+vel(iDim)) / (rhoL+rho);
        dvel_ducL(iDim,nDim+1) = 0;
        for (int jDim = 0; jDim < nDim; ++jDim) dvel_ducL(iDim,jDim+1) = 0;
        dvel_ducL(iDim,iDim+1) = 1.0/(rhoL+rho);
    }

    //! dH/ducL
    VectorDbl<nVar> dH_ducL;
    dH_ducL(0) = ( 0.5*(HL-H) - aL*aL/(gamma-1) + 0.5*(gamma-2)*q2L ) / (rhoL+rho);
    for (int iDim = 0; iDim < nDim; ++iDim)
        dH_ducL(iDim+1) = ( 1 - gamma )*velL(iDim) / (rhoL+rho);
    dH_ducL(nDim+1) =              gamma / (rhoL+rho);

    //! d(rhoR-rhoL)/ducL = - drhoL/ducL = - (drhoL/dWL)*(dWL/ducL) = -(1,0,0,0,0)*dW/dU
    VectorDbl<nVar> ddrho_ducL;
    ddrho_ducL(0) = - (  1 );
    for (int iDim = 0; iDim < nDim; ++iDim)
        ddrho_ducL(iDim+1) = 0;
    ddrho_ducL(nDim+1) = 0;

    //! d(pR-pL)/ducL = - dpL/ducL = - (dpL/dWL)*(dWL/ducL) = -(0,0,0,0,1)*dW/dU
    VectorDbl<nVar> ddp_ducL;
    ddp_ducL(0)   = - ( 0.5*(gamma-1)*q2L );
    for (int iDim = 0; iDim < nDim; ++iDim)
        ddp_ducL(iDim+1)   = - (    - (gamma-1)*velL(iDim)  );
    ddp_ducL(nDim+1)   = - (       gamma-1      );

    //! d(qnR-qnL)/ducL = - dqnL/ducL = - (dqnL/dWL)*(dWL/ducL) = -(0,nx,ny,nz,0)*dW/dU
    VectorDbl<nVar> ddqn_ducL;
    ddqn_ducL(0)  = - (-qnL/rhoL);
    for (int iDim = 0; iDim < nDim; ++iDim)
        ddqn_ducL(iDim+1)  = - (  n(iDim)/rhoL);
    ddqn_ducL(nDim+1)  = - ( 0    );

    //! d(uR-uL)/ducL = - duL/ducL = - (duL/dWL)*(dWL/ducL) = -(0,1,0,0,0)*dW/dU
    MatrixDbl<nDim,nVar> ddvel_ducL;
    for (int iDim = 0; iDim < nDim; ++iDim) {
        ddvel_ducL(iDim,0) = velL(iDim)/rhoL;
        ddvel_ducL(iDim,nDim+1) = 0;
        for (int jDim = 0; jDim < nDim; ++jDim) ddvel_ducL(iDim,jDim+1) = 0;
        ddvel_ducL(iDim,iDim+1) = -1.0/rhoL;
    }


    VectorDbl<nVar> ddws1_ducL;
    for(int k = 0; k < nVar; ++k)
      ddws1_ducL(k) = elimc_nonlinear*da_ducL(k);

    VectorDbl<nVar> ddws3_ducL;
    for(int k = 0; k < nVar; ++k)
      ddws3_ducL(k) = elimc_nonlinear*da_ducL(k);

    VectorDbl<nVar> ddws2_ducL;
    for(int k = 0; k < nVar; ++k)
      ddws2_ducL(k) = elimc_linear*da_ducL(k);

    VectorDbl<nVar> ddws4_ducL;
    for(int k = 0; k < nVar; ++k)
      ddws4_ducL(k) = elimc_linear*da_ducL(k);

    //!--------------------------------------
    //! 1.3 Differentiate the absolute values of the wave speeds, and
    //!     add the second term, - 0.5*sum_{k=1,4} [d(|lambda_k|)/ducL]*(LdU)_k*r_k

    //!  dws1_ducL = d(|qn-a|)/dUcL
    //!
    //!  Note on entropy fix:
    //!　 Let   dws1' = 0.5 * ( ws(1)*ws(1)/dws(1)+dws(1) ) <- Entropy fix
    ////! 　Then, dws1'/dUcL = (dws1'/dws1) * dws1_ducL
    //!                    = ws(1)/dws(1) * dws1_ducL

    //!  Absolute value

    VectorDbl<nVar> dws1_ducL;
    mask = qn-a > 0;
    Double s = (mask) * 1 + (1-mask) * -1;
    for (int iVar = 0; iVar < nVar; ++iVar) {
      dws1_ducL(iVar) = s * (dqn_ducL(iVar) - da_ducL(iVar));
    }


    //!  Entropy fix/Eigenvalue limiting
    Double fix = ws_orig(0) < dws(0);
    for (int iVar = 0; iVar < nVar; ++iVar)
      dws1_ducL(iVar) = fix * (ws_orig(0)/dws(0) * dws1_ducL(iVar) + 0.5 * (-ws_orig(0)*ws_orig(0)/(dws(0)*dws(0)) + 1)*ddws1_ducL(iVar)) +
                        (1-fix) * dws1_ducL(iVar);


    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) -= 0.5 * dws1_ducL(jVar)*LdU(0)*R(iVar,0);


    //!  dws2_ducL = d(|qn|)/dUcL

    auto dws2_ducL = dabs_qn_ducL;

    //!  Entropy fix/Eigenvalue limiting
    fix = ws_orig(0) < dws(0);
    for (int iVar = 0; iVar < nVar; ++iVar)
      dws2_ducL(iVar) = fix * (ws_orig(1)/dws(1) * dws2_ducL(iVar) + 0.5 * (-ws_orig(1)*ws_orig(1)/(dws(1)*dws(1)) + 1)*ddws2_ducL(iVar)) +
                        (1-fix) * dws2_ducL(iVar);


    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) -= 0.5 * dws2_ducL(jVar)*LdU(1)*R(iVar,1);


    //!  dws3_ducL = d(|qn+a|)/dUcL
    //!
    //!  Note on entropy fix:
    //!　 Let   dws3' = 0.5 * ( ws(3)*ws(3)/dws(3)+dws(3) ) <- Entropy fix
    //! 　Then, dws3'/dUcL = (dws3'/dws3) * dws3_ducL
    //!                    = ws(3)/dws(3) * dws3_ducL

    //!  Absolute value
    VectorDbl<nVar> dws3_ducL;
    mask = qn+a > 0;
    s = mask * 1 + (1-mask) * -1;
    for (int iVar = 0; iVar < nVar; ++iVar) {
      dws3_ducL(iVar) = s * (dqn_ducL(iVar) + da_ducL(iVar));
    }


    //!  Entropy fix/Eigenvalue limiting
    fix = ws_orig(2) < dws(2);
    for (int iVar = 0; iVar < nVar; ++iVar)
      dws3_ducL(iVar) = fix * (ws_orig(2)/dws(2) * dws3_ducL(iVar) + 0.5 * (-ws_orig(2)*ws_orig(2)/(dws(2)*dws(2)) + 1)*ddws3_ducL(iVar))
                        + (1-fix) * dws3_ducL(iVar);

    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) -= 0.5 * dws3_ducL(jVar)*LdU(2)*R(iVar,2);

    //!  dws4_ducL = d(|qn|)/dUcL = dws1_ducL

    auto dws4_ducL = dabs_qn_ducL;

    //!  Entropy fix/Eigenvalue limiting
    fix = ws_orig(3) < dws(3);
    for (int iVar = 0; iVar < nVar; ++iVar)
      dws4_ducL(iVar) = fix * (ws_orig(3)/dws(3) * dws4_ducL(iVar) + 0.5 * (-ws_orig(3)*ws_orig(3)/(dws(3)*dws(3)) + 1)*ddws4_ducL(iVar)) +
                        (1-fix) * dws4_ducL(iVar);


    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) -= 0.5 * dws4_ducL(jVar)*LdU(3)*R(iVar,3);

    //!--------------------------------------
    //! 1.4 Differentiate the wave strength, and
    //!     add the third term, - 0.5*sum_{k=1,4} |lambda_k|*[d(LdU)_k/ducL]*r_k.

    //!  dLdU1_ducL = d( dp/(2a^2) - rho/(2a)*dqn )/dUcL
    //!
    //!             = dp*[d(1/(2a^2))/dUcL] - dqn*[d(rho/(2a))/dUcL]
    //!             + [d(dp)/dUcL]/(2a^2)   - [d(dqn)/dUcL]*rho/(2a)
    //!
    //!             = -a^(-3)*dp*[da/dUcL]  - dqn*( [drho/dUcL]/(2a) - rho/(2*a^2)*[da/dUcL] )
    //!             + [d(dp)/dUcL]/(2a^2)   - [d(dqn)/dUcL]*rho/(2a)
    //!
    //!             = ( -2*dp + rho*a*dqn )/(2a^3)*[da/dUcL] - dqn*[drho/dUcL]/(2a)
    //!             + [d(dp)/dUcL]/(2a^2)   - [d(dqn)/dUcL]*rho/(2a)

    VectorDbl<nVar> dLdU1_ducL;
    for(int iVar = 0; iVar < nVar; ++iVar) {
      dLdU1_ducL(iVar) = 0.5*(-2*dp+rho*a*dqn )/(a2*a) * (da_ducL(iVar))
                         - 0.5*dqn/a * (drho_ducL(iVar))
                         + 0.5*(ddp_ducL(iVar))/a2
                         - 0.5*rho*(ddqn_ducL(iVar))/a;
    }

    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) -= 0.5 * ws(0)*dLdU1_ducL(jVar)*R(iVar,0);

    //!  dLdU2_ducL = d( drho - dp/a^2 )/dUcL
    //!              = [d(drho)/dUcL] - [d(dp)/dUcL]/a^2 + 2*dp/a^3*[da/dUcL]

    VectorDbl<nVar> dLdU2_ducL;
    for(int iVar = 0; iVar < nVar; ++iVar)
      dLdU2_ducL(iVar) = (ddrho_ducL(iVar)) - (ddp_ducL(iVar))/a2 + 2*dp/(a2*a) * (da_ducL(iVar));

    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) -= 0.5 * ws(1)*dLdU2_ducL(jVar)*R(iVar,1);

    //!  dLdU3_ducL = d( dp/(2a^2) + rho/(2a)*dqn )/dUcL
    //!
    //!             = dp*[d(1/(2a^2))/dUcL] + dqn*[d(rho/(2a))/dUcL]
    //!             + [d(dp)/dUcL]/(2a^2)   + [d(dqn)/dUcL]*rho/(2a)
    //!
    //!             = -a^(-3)*dp*[da/dUcL]  + dqn*( [drho/dUcL]/(2a) - rho/(2*a^2)*[da/dUcL] )
    //!             + [d(dp)/dUcL]/(2a^2)   + [d(dqn)/dUcL]*rho/(2a)
    //!
    //!             = ( -2*dp - rho*a*dqn )/(2a^3)*[da/dUcL]  + dqn*[drho/dUcL]/(2a)
    //!             + [d(dp)/dUcL]/(2a^2)   + [d(dqn)/dUcL]*rho/(2a)

    VectorDbl<nVar> dLdU3_ducL;
    for(int iVar = 0; iVar < nVar; ++iVar) {
      dLdU3_ducL(iVar) = 0.5 * (-2 * dp - rho * a * dqn) / (a2*a) * (da_ducL(iVar)) + 0.5 * dqn / a * (drho_ducL(iVar)) +
                         0.5 * (ddp_ducL(iVar)) / a2 + 0.5 * rho * (ddqn_ducL(iVar)) / a;
    }

    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) -= 0.5 * ws(2)*dLdU3_ducL(jVar)*R(iVar,2);

    //!  dLdU4_ducL = d(rho)/dUcL

    auto dLdU4_ducL = drho_ducL;

    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) -= 0.5 * ws(3)*dLdU4_ducL(jVar)*R(iVar,3);

    //!--------------------------------------
    //! 1.5 Differentiate the right-eigenvectors, and
    //!     add the fourth term, - 0.5*sum_{k=1,4} |lambda_k|*(LdU)_k*[dr_k/ducL]

    //! dR1_ducL = dR(:,1)/dUcL
    //!
    //! Left-moving acoustic wave
    //!
    //!        Eigenvector -> Differentiated
    //!
    //!  R(1,1) = 1      -> 0
    //!  R(2,1) = u - a*nx -> du/dUcL - da/dUcL*nx
    //!  R(3,1) = v - a*ny -> dv/dUcL - da/dUcL*ny
    //!  R(4,1) = w - a*nz -> dw/dUcL - da/dUcL*nz
    //!  R(5,1) = H - a*qn -> dH/dUcL - da/dUcL*qn - dqn/dUcL*a
    MatrixDbl<nDim+2> dR1_ducL;
    for (int iVar = 0; iVar < nVar; ++iVar) dR1_ducL(0,iVar) = 0;
    for (int iDim = 0; iDim < nDim; ++iDim)
        for (int iVar = 0; iVar < nVar; ++iVar)
            dR1_ducL(iDim+1,iVar) = (dvel_ducL(iDim,iVar)) - (da_ducL(iVar)) * n(iDim);
    for (int iVar = 0; iVar < nVar; ++iVar) dR1_ducL(nDim+1,iVar) = (dH_ducL(iVar)) - (da_ducL(iVar)) * qn - (dqn_ducL(iVar)) * a;

    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) -= 0.5 * ws(0)*LdU(0)*dR1_ducL(iVar,jVar);

    //! dR2_ducL = dR(:,2)/dUcL
    //!
    //! Entropy wave
    //!
    //!                  Eigenvector -> Differentiated
    //!
    //!  R(1,2) = 1                -> 0
    //!  R(2,2) = u                  -> du/dUcL
    //!  R(3,2) = v                  -> dv/dUcL
    //!  R(4,2) = w                  -> dw/dUcL
    //!  R(5,2) = 0.5*(u*u+v*v+w*w) -> u*du/dUcL + v*dv/dUcL + w*dw/dUcL
    MatrixDbl<nDim+2> dR2_ducL;
    for (int iVar = 0; iVar < nVar; ++iVar) dR2_ducL(0,iVar) = 0;
    for (int iDim = 0; iDim < nDim; ++iDim)
        for (int iVar = 0; iVar < nVar; ++iVar)
            dR2_ducL(iDim+1,iVar) = (dvel_ducL(iDim,iVar));
    for (int iVar = 0; iVar < nVar; ++iVar) {
        dR2_ducL(nDim+1,iVar) = 0;
        for (int iDim = 0; iDim < nDim; ++iDim)
            dR2_ducL(nDim+1,iVar) += vel(iDim)*(dvel_ducL(iDim,iVar));
    }

    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) -= 0.5 * ws(1)*LdU(1)*dR2_ducL(iVar,jVar);

    //! dR3_ducL = dR(:,3)/dUcL
    //!
    //! Right-moving acoustic wave
    //!
    //!        Eigenvector -> Differentiated
    //!
    //!  R(1,3) = 1      -> 0
    //!  R(2,3) = u + a*nx -> du/dUcL + da/dUcL*nx
    //!  R(3,3) = v + a*ny -> dv/dUcL + da/dUcL*ny
    //!  R(4,3) = w + a*nz -> dw/dUcL + da/dUcL*nz
    //!  R(5,3) = H + a*qn -> dH/dUcL + da/dUcL*qn + dqn/dUcL*a

    MatrixDbl<nDim+2> dR3_ducL;
    for (int iVar = 0; iVar < nVar; ++iVar) dR3_ducL(0,iVar) = 0;
    for (int iDim = 0; iDim < nDim; ++iDim)
        for (int iVar = 0; iVar < nVar; ++iVar)
            dR3_ducL(iDim+1,iVar) = (dvel_ducL(iDim,iVar)) + (da_ducL(iVar)) * n(iDim);
    for (int iVar = 0; iVar < nVar; ++iVar) dR3_ducL(nDim+1,iVar) = (dH_ducL(iVar)) + (da_ducL(iVar)) * qn + (dqn_ducL(iVar)) * a;

    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) -= 0.5 * ws(2)*LdU(2)*dR3_ducL(iVar,jVar);

    //! dR4_ducL = dR(:,4)/dUcL
    //! Two shear wave components combined into 1 (wave strength incorporated).
    //! So, it is not really an eigenvector.
    //!
    //!                  Combined vector -> Differentiated
    //!
    //!  R(1,4) = 0                   -> 0
    //!  R(2,4) = du - dqn*nx            -> d(du)/dUcL - d(dqn)/dUcL*nx
    //!  R(3,4) = dv - dqn*ny            -> d(dv)/dUcL - d(dqn)/dUcL*ny
    //!  R(4,4) = dw - dqn*nz            -> d(dw)/dUcL - d(dqn)/dUcL*nz
    //!  R(5,4) = u*du+v*dv+w*dw-qn*dqn  -> du/dUcL*du     + d(du)/dUcL*u
    //!                                   + dv/dUcL*dv     + d(dv)/dUcL*v
    //!                                   + dw/dUcL*dw     + d(dw)/dUcL*w
    //!                                   - d(qn)/dUcL*dqn - d(dqn)/dUcL*qn

    MatrixDbl<nDim+2> dR4_ducL;
    for (int iVar = 0; iVar < nVar; ++iVar) dR4_ducL(0,iVar) = 0;
    for (int iDim = 0; iDim < nDim; ++iDim)
        for (int iVar = 0; iVar < nVar; ++iVar)
            dR4_ducL(iDim+1,iVar) = (ddvel_ducL(iDim,iVar)) - (ddqn_ducL(iVar))*n(iDim);

    for (int iVar = 0; iVar < nVar; ++iVar) {
        dR4_ducL(nDim+1,iVar) = - dqn*( dqn_ducL(iVar)) - qn*(ddqn_ducL(iVar));
        for (int iDim = 0; iDim < nDim; ++iDim)
            dR4_ducL(nDim+1,iVar) += dVel(iDim)*( dvel_ducL(iDim,iVar)) + vel(iDim)*(ddvel_ducL(iDim,iVar));
    }

    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducL(iVar,jVar) -= 0.5 * ws(3)*LdU(3)*dR4_ducL(iVar,jVar);


    //!--------------------------------------------------------------------------------
    //! Part 2. Compute dFn_Roe/ducR
    //!
    //!  dFn_Roe/ducR =   d(0.5*Fn(ucR))/duR
    //!                 - 0.5*sum_{k=1,4} [d(|lambda_k|)/ducR]*(LdU)_k*r_k
    //!                 - 0.5*sum_{k=1,4} |lambda_k|*[d(LdU)_k/ducR]*r_k
    //!                 - 0.5*sum_{k=1,4} |lambda_k|*(LdU)_k*[dr_k/ducR]
    //!
    //!  So, we proceed as follows:
    //!
    //!  1.1 Compute                d(0.5*Fn(ucR))/duR
    //!  1.2 Compute various deriavives that will be used in the following steps.
    //!  1.3 Add the second term, - 0.5*sum_{k=1,4} [d(|lambda_k|)/ducR]*(LdU)_k*r_k
    //!  1.4 Add the  third term, - 0.5*sum_{k=1,4} |lambda_k|*[d(LdU)_k/ducR]*r_k
    //!  1.5 Add the fourth term, - 0.5*sum_{k=1,4} |lambda_k|*(LdU)_k*[dr_k/ducR]
    //!

    //!--------------------------------------
    //! 2.1 Compute "d(0.5*Fn(ucR))/ducR"
    //!
    //!     (See "I Do Like CFD, VOL.1", page 55, for the analytical Jacobian, dFn(u)/du)

    //!  1st column
    MatrixDbl<nVar> dFnducR;
    dFnducR(0,0) = 0;
    for (int iDim = 0; iDim < nDim; ++iDim)
        dFnducR(iDim+1,0) = 0.5*(gamma-1)*q2R*n(iDim)  - velR(iDim)*qnR;
    dFnducR(nDim+1,0) = 0.5*(gamma-1)*q2R*qnR - HR*qnR;

    //!  2nd column
    for (int jDim = 0; jDim < nDim; ++jDim) {
        dFnducR(0,jDim+1) =     n(jDim);
        dFnducR(nDim+1,jDim+1) =  HR*n(jDim) - (gamma-1)*velR(jDim)*qnR;
        for (int iDim = 0; iDim < nDim; ++iDim)
            dFnducR(iDim+1,jDim+1) =  velR(iDim)*n(jDim) - (gamma-1)*velR(jDim)*n(iDim);
        dFnducR(jDim+1,jDim+1) += qnR;
    }

    //!  5th column
    dFnducR(0,nDim+1) =  0;
    for (int iDim = 0; iDim < nDim; ++iDim)
        dFnducR(iDim+1,nDim+1) = (gamma-1)*n(iDim);
    dFnducR(nDim+1,nDim+1) =  gamma*qnR;

    //!  Factor 1/2
    for(int iVar = 0; iVar < nVar; ++iVar)
      for(int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) *= 0.5;

    //!--------------------------------------
    //! 2.2 Compute various deriavives that will be used in the following steps.
    //!     (See "I Do Like CFD, VOL.1" for details)

    //! dqn/ducR

    VectorDbl<nVar> dqn_ducR;
    dqn_ducR(0) = -0.5*(qnR+qn) / (rhoR+rho);
    for (int iDim = 0; iDim < nDim; ++iDim)
        dqn_ducR(iDim+1) =            n(iDim)  / (rhoR+rho);
    dqn_ducR(nDim+1) =          0;

    //! d(|qn|)/ducR
    VectorDbl<nVar> dabs_qn_ducR;
    mask = qn < 0;
    sigQN = mask * -1 + (1-mask) * 1;
    dabs_qn_ducR(0) = -0.5*sigQN*(qnR+qn) / (rhoR+rho);
    for (int iDim = 0; iDim < nDim; ++iDim)
        dabs_qn_ducR(iDim+1) =       sigQN*     n(iDim)  / (rhoR+rho);
    dabs_qn_ducR(nDim+1) =  0;

    //! da/ducR
    VectorDbl<nVar> da_ducR;
    const Double velDotVelR = dot(vel,velR);
    da_ducR(0) =  0.5*(gamma-1)/a*( 0.5*( velDotVelR + q2 )
                                          +  0.5*(HR-H) - aR*aR/(gamma-1) + 0.5*(gamma-2)*q2R )/(rhoR+rho);
    for (int iDim = 0; iDim < nDim; ++iDim)
        da_ducR(iDim+1) = -0.5*(gamma-1)*(vel(iDim)+(gamma-1)*velR(iDim))/a  / (rhoR+rho);
    da_ducR(nDim+1) =  0.5*gamma*(gamma-1)/a  / (rhoR+rho);

    //! drho/ducR
    VectorDbl<nVar> drho_ducR;
    drho_ducR(0) = 0.5*rho/rhoR;
    for (int iDim = 0; iDim < nDim; ++iDim) drho_ducR(iDim+1) = 0;
    drho_ducR(nDim+1) = 0;

    //! du/ducR
    MatrixDbl<nDim,nVar> dvel_ducR;

    for (int iDim = 0; iDim < nDim; ++iDim) {
      dvel_ducR(iDim,0) = -0.5*(velR(iDim)+vel(iDim)) / (rhoR+rho);
      dvel_ducR(iDim,nDim+1) = 0;
      for (int jDim = 0; jDim < nDim; ++jDim) dvel_ducR(iDim,jDim+1) = 0;
      dvel_ducR(iDim,iDim+1) = 1.0 / (rhoR+rho);
    }

    //! dH/ducR
    VectorDbl<nVar> dH_ducR;
    dH_ducR(0) = ( 0.5*(HR-H) - aR*aR/(gamma-1) + 0.5*(gamma-2)*q2R ) / (rhoR+rho);
    for (int iDim = 0; iDim < nDim; ++iDim)
      dH_ducR(iDim+1) = ( 1 - gamma )*velR(iDim) / (rhoR+rho);
    dH_ducR(nDim+1) =              gamma / (rhoR+rho);

    //! d(rhoR-rhoR)/ducR = drhoR/ducR = (drhoR/dWR)*(dWR/ducR) = (1,0,0,0,0)*dW/dU
    VectorDbl<nVar> ddrho_ducR;
    ddrho_ducR(0) =  (  1 );
    for (int iDim = 0; iDim < nDim; ++iDim)
      ddrho_ducR(iDim+1) =  ( 0 );
    ddrho_ducR(nDim+1) =  ( 0 );

    //! d(pR-pR)/ducR = dpR/ducR = (dpR/dWR)*(dWR/ducR) = (0,0,0,0,1)*dW/dU
    VectorDbl<nVar> ddp_ducR;
    ddp_ducR(0)   =  ( 0.5*(gamma-1)*q2R );
    for (int iDim = 0; iDim < nDim; ++iDim)
      ddp_ducR(iDim+1)   =  (    - (gamma-1)*velR(iDim)  );
    ddp_ducR(nDim+1)   =  (       gamma-1      );

    //! d(qnR-qnR)/ducR = dqnR/ducR = (dqnR/dWR)*(dWR/ducR) = (0,nx,ny,nz,0)*dW/dU
    VectorDbl<nVar> ddqn_ducR;
    ddqn_ducR(0)  =  (-qnR/rhoR);
    for (int iDim = 0; iDim < nDim; ++iDim)
      ddqn_ducR(iDim+1)  =  (  n(iDim)/rhoR);
    ddqn_ducR(nDim+1)  =  ( 0    );

    //! d(uR-uR)/ducR = duR/ducR = (duR/dWR)*(dWR/ducR) = -(0,1,0,0,0)*dW/dU
    MatrixDbl<nDim,nVar> ddvel_ducR;
    for (int iDim = 0; iDim < nDim; ++iDim) {
      ddvel_ducR(iDim,0) = -velR(iDim) / rhoR;
      ddvel_ducR(iDim,nDim+1) = 0;
      for (int jDim = 0; jDim < nDim; ++jDim) ddvel_ducR(iDim,jDim+1) = 0;
      ddvel_ducR(iDim,iDim+1) = 1.0 / rhoR;
    }


    VectorDbl<nVar> ddws1_ducR;
    for (int iVar = 0; iVar < nVar; ++iVar) ddws1_ducR(iVar) = elimc_nonlinear*da_ducR(iVar);

    VectorDbl<nVar> ddws3_ducR;
    for (int iVar = 0; iVar < nVar; ++iVar) ddws3_ducR(iVar) = elimc_nonlinear*da_ducR(iVar);

    VectorDbl<nVar> ddws2_ducR;
    for (int iVar = 0; iVar < nVar; ++iVar) ddws2_ducR(iVar) = elimc_linear*da_ducR(iVar);

    VectorDbl<nVar> ddws4_ducR;
    for (int iVar = 0; iVar < nVar; ++iVar) ddws4_ducR(iVar) = elimc_linear*da_ducR(iVar);

    //!--------------------------------------
    //! 2.3 Differentiate the absolute values of the wave speeds, and
    //!     add the second term, - 0.5*sum_{k=1,4} [d(|lambda_k|)/ducR]*(LdU)_k*r_k
    //
    //!  dws1_ducR = d(|qn-a|)/dUcR
    //!
    //!  Note on entropy fix:
    //!　 Let   dws1' = 0.5 * ( ws(1)*ws(1)/dws(1)+dws(1) ) <- Entropy fix
    //! 　Then, dws1'/dUcR = (dws1'/dws1) * dws1_ducR    + ( -ws(1)*ws(1)/dws(1)**2*ddws(1) + ddws(1) )
    //!                    = ws(1)/dws(1) * dws1_ducR

    //!  Absolute value

    VectorDbl<nVar> dws1_ducR;
    mask = qn-a > 0;
    s = (mask) * 1 + (1-mask) * -1;
    for (int iVar = 0; iVar < nVar; ++iVar)
      dws1_ducR(iVar) =  s * (dqn_ducR(iVar) - da_ducR(iVar));


    //!  Entropy fix/Eigenvalue limiting

    fix = ws_orig(0) < dws(0);
    for (int iVar = 0; iVar < nVar; ++iVar)
      dws1_ducR(iVar) = fix * (ws_orig(0)/dws(0) * dws1_ducR(iVar) + 0.5 * (-ws_orig(0)*ws_orig(0)/(dws(0)*dws(0)) + 1)*ddws1_ducR(iVar)) +
                        (1-fix) * dws1_ducR(iVar);

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) -= 0.5 * dws1_ducR(jVar)*LdU(0)*R(iVar,0);

    //!  dws2_ducR = d(|qn|)/dUcR

    auto dws2_ducR = dabs_qn_ducR;

    //!  Entropy fix/Eigenvalue limiting
    fix = ws_orig(1) < dws(1);
    for (int iVar = 0; iVar < nVar; ++iVar)
      dws2_ducR(iVar) = fix * (ws_orig(1)/dws(1) * dws2_ducR(iVar) + 0.5 * (-ws_orig(1)*ws_orig(1)/(dws(1)*dws(1)) + 1)*ddws2_ducR(iVar)) +
                        (1-fix) * dws2_ducR(iVar);

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) -= 0.5 * dws2_ducR(jVar)*LdU(1)*R(iVar,1);

    //!  dws3_ducR = d(|qn+a|)/dUcR
    //!
    //!  Note on entropy fix:
    //!　 Let   dws3' = 0.5 * ( ws(3)*ws(3)/dws(3)+dws(3) ) <- Entropy fix
    //! 　Then, dws3'/dUcR = (dws3'/dws3) * dws3_ducR
    //!                    = ws(3)/dws(3) * dws3_ducR

    //!  Absolute value
    VectorDbl<nVar> dws3_ducR;
    mask = qn+a > 0;
    s = mask * 1 + (1-mask) * -1;
    for (int iVar = 0; iVar < nVar; ++iVar)
      dws3_ducR(iVar) =  s * (dqn_ducR(iVar) + da_ducR(iVar));


    //!  Entropy fix/Eigenvalue limiting

    fix = ws_orig(2) < dws(2);
    for (int iVar = 0; iVar < nVar; ++iVar)
      dws3_ducR(iVar) = fix * (ws_orig(2)/dws(2) * dws3_ducR(iVar) + 0.5 * (-ws_orig(2)*ws_orig(2)/(dws(2)*dws(2)) + 1)*ddws3_ducR(iVar)) +
                        (1-fix) * dws3_ducR(iVar);

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) -= 0.5 * dws3_ducR(jVar)*LdU(2)*R(iVar,2);

    //!  dws4_ducR = d(|qn|)/dUcR = dws1_ducR

    auto dws4_ducR = dabs_qn_ducR;

    //!  Entropy fix/Eigenvalue limiting
    fix = ws_orig(3) < dws(3);
    for (int iVar = 0; iVar < nVar; ++iVar)
      dws4_ducR(iVar) = fix * (ws_orig(3)/dws(3) * dws4_ducR(iVar) + 0.5 * (-ws_orig(3)*ws_orig(3)/(dws(3)*dws(3)) + 1)*ddws4_ducR(iVar)) +
                        (1-fix) * dws4_ducR(iVar);

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) -= 0.5 * dws4_ducR(jVar)*LdU(3)*R(iVar,3);

    //!--------------------------------------
    //! 2.4 Differentiate the wave strength, and
    //!     add the third term, - 0.5*sum_{k=1,4} |lambda_k|*[d(LdU)_k/ducR]*r_k.
    //
    //!  dLdU1_ducR = d( dp/(2a^2) - rho/(2a)*dqn )/dUcR
    //!
    //!             = dp*[d(1/(2a^2))/dUcR] - dqn*[d(rho/(2a))/dUcR]
    //!             + [d(dp)/dUcR]/(2a^2)   - [d(dqn)/dUcR]*rho/(2a)
    //!
    //!             = -a^(-3)*dp*[da/dUcR]  - dqn*( [drho/dUcR]/(2a) - rho/(2*a^2)*[da/dUcR] )
    //!             + [d(dp)/dUcR]/(2a^2)   - [d(dqn)/dUcR]*rho/(2a)
    //!
    //!             = ( -2*dp + rho*a*dqn )/(2a^3)*[da/dUcR] - dqn*[drho/dUcR]/(2a)
    //!             + [d(dp)/dUcR]/(2a^2)   - [d(dqn)/dUcR]*rho/(2a)

    VectorDbl<nVar> dLdU1_ducR;
    for (int iVar = 0; iVar < nVar; ++iVar) {
      dLdU1_ducR(iVar) = 0.5 * (-2 * dp + rho * a * dqn) / (a2 * a) * (da_ducR(iVar)) -
                         0.5 * dqn / a * (drho_ducR(iVar)) + 0.5 * (ddp_ducR(iVar)) / a2 -
                         0.5 * rho * (ddqn_ducR(iVar)) / a;
    }

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) -= 0.5 * ws(0)*dLdU1_ducR(jVar)*R(iVar,0);

    //!  dLdU2_ducR = d( drho - dp/a^2 )/dUcR
    //!              = [d(drho)/dUcR] - [d(dp)/dUcR]/a^2 + 2*dp/a^3*[da/dUcR]

    VectorDbl<nVar> dLdU2_ducR;
    for (int iVar = 0; iVar < nVar; ++iVar)
      dLdU2_ducR(iVar) = (ddrho_ducR(iVar)) - (ddp_ducR(iVar))/a2 + 2*dp/(a2*a) * (da_ducR(iVar));

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) = dFnducR(iVar,jVar) - 0.5 * ws(1)*dLdU2_ducR(jVar)*R(iVar,1);

    //!  dLdU3_ducR = d( dp/(2a^2) + rho/(2a)*dqn )/dUcR
    //!
    //!             = dp*[d(1/(2a^2))/dUcR] + dqn*[d(rho/(2a))/dUcR]
    //!             + [d(dp)/dUcR]/(2a^2)   + [d(dqn)/dUcR]*rho/(2a)
    //!
    //!             = -a^(-3)*dp*[da/dUcR]  + dqn*( [drho/dUcR]/(2a) - rho/(2*a^2)*[da/dUcR] )
    //!             + [d(dp)/dUcR]/(2a^2)   + [d(dqn)/dUcR]*rho/(2a)
    //!
    //!             = ( -2*dp - rho*a*dqn )/(2a^3)*[da/dUcR]  + dqn*[drho/dUcR]/(2a)
    //!             + [d(dp)/dUcR]/(2a^2)   + [d(dqn)/dUcR]*rho/(2a)

    VectorDbl<nVar> dLdU3_ducR;
    for (int iVar = 0; iVar < nVar; ++iVar) {
      dLdU3_ducR(iVar) = 0.5 * (-2 * dp - rho * a * dqn) / (a2 * a) * (da_ducR(iVar)) +
                         0.5 * dqn / a * (drho_ducR(iVar)) + 0.5 * (ddp_ducR(iVar)) / a2 +
                         0.5 * rho * (ddqn_ducR(iVar)) / a;
    }

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) = dFnducR(iVar,jVar) - 0.5 * ws(2)*dLdU3_ducR(jVar)*R(iVar,2);

    //!  dLdU4_ducR = d(rho)/dUcR

    auto dLdU4_ducR = drho_ducR;

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) = dFnducR(iVar,jVar) - 0.5 * ws(3)*dLdU4_ducR(jVar)*R(iVar,3);

    //!--------------------------------------
    //! 2.5 Differentiate the right-eigenvectors, and
    //!     add the fourth term, - 0.5*sum_{k=1,4} |lambda_k|*(LdU)_k*[dr_k/ducL]
    //
    //! dR1_ducR = dR(:,1)/dUcR
    //!
    //! Reft-moving acoustic wave
    //!
    //!        Eigenvector -> Differentiated
    //!
    //!  R(1,1) = 1      -> 0
    //!  R(2,1) = u - a*nx -> du/dUcR - da/dUcR*nx
    //!  R(3,1) = v - a*ny -> dv/dUcR - da/dUcR*ny
    //!  R(4,1) = w - a*nz -> dw/dUcR - da/dUcR*nz
    //!  R(5,1) = H - a*qn -> dH/dUcR - da/dUcR*qn - dqn/dUcR*a
    MatrixDbl<nVar> dR1_ducR;
    for (int iVar = 0; iVar < nVar; ++iVar) dR1_ducR(0,iVar) = 0;
    for (int iDim = 0; iDim < nDim; ++iDim)
      for (int iVar = 0; iVar < nVar; ++iVar)
        dR1_ducR(iDim+1,iVar) = (dvel_ducR(iDim,iVar)) - (da_ducR(iVar)) * n(iDim);
    for (int iVar = 0; iVar < nVar; ++iVar) dR1_ducR(nDim+1,iVar) = (dH_ducR(iVar)) - (da_ducR(iVar)) * qn - (dqn_ducR(iVar)) * a;

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) -= 0.5 * ws(0)*LdU(0)*dR1_ducR(iVar,jVar);

    //! dR2_ducR = dR(:,2)/dUcR
    //!
    //! Entropy wave
    //!
    //!                  Eigenvector -> Differentiated
    //!
    //!  R(1,2) = 1                -> 0
    //!  R(2,2) = u                  -> du/dUcR
    //!  R(3,2) = v                  -> dv/dUcR
    //!  R(4,2) = w                  -> dw/dUcR
    //!  R(5,2) = 0.5*(u*u+v*v+w*w) -> u*du/dUcR + v*dv/dUcR + w*dw/dUcR

    MatrixDbl<nVar> dR2_ducR;
    for (int iVar = 0; iVar < nVar; ++iVar) dR2_ducR(0,iVar) = 0;
    for (int iDim = 0; iDim < nDim; ++iDim)
      for (int iVar = 0; iVar < nVar; ++iVar)
        dR2_ducR(iDim+1,iVar) = (dvel_ducR(iDim,iVar));
    for (int iVar = 0; iVar < nVar; ++iVar) {
      dR2_ducR(nDim+1, iVar)= 0;
      for (int iDim = 0; iDim < nDim; ++iDim)
        dR2_ducR(nDim+1, iVar) += vel(iDim) * (dvel_ducR(iDim,iVar));
    }

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) -= 0.5 * ws(1)*LdU(1)*dR2_ducR(iVar,jVar);

    //! dR3_ducR = dR(:,3)/dUcR
    //!
    //! Right-moving acoustic wave
    //!
    //!        Eigenvector -> Differentiated
    //!
    //!  R(1,3) = 1      -> 0
    //!  R(2,3) = u + a*nx -> du/dUcR + da/dUcR*nx
    //!  R(3,3) = v + a*ny -> dv/dUcR + da/dUcR*ny
    //!  R(4,3) = w + a*nz -> dw/dUcR + da/dUcR*nz
    //!  R(5,3) = H + a*qn -> dH/dUcR + da/dUcR*qn + dqn/dUcR*a

    MatrixDbl<nVar> dR3_ducR;
    for (int iVar = 0; iVar < nVar; ++iVar) dR3_ducR(0,iVar) = 0;
    for (int iDim = 0; iDim < nDim; ++iDim)
      for (int iVar = 0; iVar < nVar; ++iVar)
        dR3_ducR(iDim+1,iVar) = (dvel_ducR(iDim,iVar)) + (da_ducR(iVar)) * n(iDim);
    for (int iVar = 0; iVar < nVar; ++iVar) dR3_ducR(nDim+1,iVar) = (dH_ducR(iVar)) + (da_ducR(iVar)) * qn + (dqn_ducR(iVar)) * a;

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) -= 0.5 * ws(2)*LdU(2)*dR3_ducR(iVar,jVar);

    //! dR4_ducR = dR(:,4)/dUcR
    //! Two shear wave components combined into 1 (wave strength incorporated).
    //! So, it is not really an eigenvector.
    //!
    //!                  Combined vector -> Differentiated
    //!
    //!  R(1,4) = 0                   -> 0
    //!  R(2,4) = du - dqn*nx            -> d(du)/dUcR - d(dqn)/dUcR*nx
    //!  R(3,4) = dv - dqn*ny            -> d(dv)/dUcR - d(dqn)/dUcR*ny
    //!  R(4,4) = dw - dqn*nz            -> d(dw)/dUcR - d(dqn)/dUcR*nz
    //!  R(5,4) = u*du+v*dv+w*dw-qn*dqn  -> du/dUcR*du     + d(du)/dUcR*u
    //!                                   + dv/dUcR*dv     + d(dv)/dUcR*v
    //!                                   + dw/dUcR*dw     + d(dw)/dUcR*w
    //!                                   - d(qn)/dUcR*dqn - d(dqn)/dUcR*qn

    MatrixDbl<nVar> dR4_ducR;
    for (int iVar = 0; iVar < nVar; ++iVar) dR4_ducR(0,iVar) = 0;
    for (int iDim = 0; iDim < nDim; ++iDim)
      for (int iVar = 0; iVar < nVar; ++iVar)
        dR4_ducR(iDim+1,iVar) = (ddvel_ducR(iDim,iVar)) - (ddqn_ducR(iVar))*n(iDim);
    for (int iVar = 0; iVar < nVar; ++iVar){
      dR4_ducR(nDim+1,iVar) = - dqn*( dqn_ducR(iVar)) - qn*(ddqn_ducR(iVar));
      for (int iDim = 0; iDim < nDim; ++iDim)
        dR4_ducR(nDim+1,iVar) += dVel(iDim)*dvel_ducR(iDim,iVar) + vel(iDim)*ddvel_ducR(iDim,iVar);
    }

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar)
        dFnducR(iVar,jVar) -= 0.5 * ws(3)*LdU(3)*dR4_ducR(iVar,jVar);

    for (int iVar = 0; iVar < nVar; ++iVar)
      for (int jVar = 0; jVar < nVar; ++jVar) {
        dFnducL(iVar, jVar) *= area;
        dFnducR(iVar, jVar) *= area;
      }


    //!Dissipation Term: |An|(UR-UL) = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ]
    VectorDbl<nVar> diss;
    for (int iVar = 0; iVar < nVar; ++iVar) {
      diss(iVar) = 0;
      for (int iWave = 0; iWave < 4; ++iWave)
        diss(iVar) += ws(iWave)*LdU(iWave)*R(iVar,iWave);
    }

    //!Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)

    VectorDbl<nVar> flux, flux_i, flux_j;

    flux_i(0) = rhoL*qnL;
    for (int iDim = 0; iDim < nDim; ++iDim)
      flux_i(iDim+1) = rhoL*qnL * velL(iDim) + pL*n(iDim);
    flux_i(nDim+1) = rhoL*qnL * HL;

    flux_j(0) = rhoR*qnR;
    for (int iDim = 0; iDim < nDim; ++iDim)
      flux_j(iDim+1) = rhoR*qnR * velR(iDim) + pR*n(iDim);
    flux_j(nDim+1) = rhoR*qnR * HR;

    //! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]
    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      flux(iVar) = 0.5 * (flux_i(iVar) + flux_j(iVar) - diss(iVar)) * area;
    }

    /*--- Inviscid fluxes and Jacobians. ---*/

    MatrixDbl<nVar> jac_i, jac_j;
    if (implicit) {
      jac_i = dFnducL;
      jac_j = dFnducR;
    }

    /*--- Correct for grid motion. ---*/

    if (dynamicGrid) {
      for (size_t iVar = 0; iVar < nVar; ++iVar) {
        Double dFdU = projGridVel * area * 0.5;
        flux(iVar) -= dFdU * (U.i.all(iVar) + U.j.all(iVar));

        if (implicit) {
          jac_i(iVar,iVar) -= dFdU;
          jac_j(iVar,iVar) -= dFdU;
        }
      }
    }

    /*--- Correct for second order on arbitrary grids ---*/
    if (correction) {
      const auto corrX = gatherVariables<nDim>(iEdge, geometry.edges->GetCorrection_X());
      const auto corrY = gatherVariables<nDim>(iEdge, geometry.edges->GetCorrection_Y());
      const auto corrZ = gatherVariables<nDim>(iEdge, geometry.edges->GetCorrection_Z());
      CorrectFlux(iPoint, jPoint, solution, V, flux, corrX, corrY, corrZ, gamma);
      //      reconstructPrimitives<CCompressiblePrimitives<nDim,nPrimVarGrad> >(iEdge, iPoint, jPoint, muscl, typeLimiter, V1st, vector_ij, solution);
    }

    /*--- Add the contributions from the base class (static decorator). ---*/

    Base::viscousTerms(iEdge, iPoint, jPoint, V1st, solution_, vector_ij, geometry,
                       config, area, unitNormal, implicit, flux, jac_i, jac_j);

    /*--- Stop preaccumulation. ---*/

    stopPreacc(flux);

    /*--- Update the vector and system matrix. ---*/

    updateLinearSystem(iEdge, iPoint, jPoint, implicit, updateType,
                       updateMask, flux, jac_i, jac_j, vector, matrix);
  }

};

template<class Decorator>
class CRoeSchemeNew : public CRoeNewBase<CRoeSchemeNew<Decorator>,Decorator> {
 private:
  using Base = CRoeNewBase<CRoeSchemeNew<Decorator>, Decorator>;
  using Base::nDim;
  using Base::nVar;
  using Base::gamma;
  using Base::kappa;
  const ENUM_ROELOWDISS typeDissip;

 public:
  /*!
   * \brief Constructor, store some constants and forward to base.
   */
  template <class... Ts>
  CRoeSchemeNew(const CConfig& config, Ts&... args)
      : Base(config, args...), typeDissip(static_cast<ENUM_ROELOWDISS>(config.GetKind_RoeLowDiss())) {}

  /*!
   * \brief Updates flux and Jacobians with standard Roe dissipation.
   * \note "Ts" is here just in case other schemes in the family need extra args.
   */
  template <class PrimVarType, class ConsVarType, class... Ts>
  FORCEINLINE void finalizeFlux(VectorDbl<nVar>& flux, MatrixDbl<nVar>& jac_i, MatrixDbl<nVar>& jac_j, bool implicit,
                                Double area, const VectorDbl<nDim>& unitNormal, const CPair<PrimVarType>& V,
                                const CPair<ConsVarType>& U, const CRoeVariables<nDim>& roeAvg,
                                const VectorDbl<nVar>& lambda, const MatrixDbl<nVar>& pMat, Int iPoint, Int jPoint,
                                const CEulerVariable& solution, Ts&...) const {};
};



/********************************************************************************/
//!* -- 3D Roe Flux Jacobian --
//!*
//!* This subroutine computes the left and right Jacobians for the Roe flux
//!* for the Euler equations in the direction, njk=[nx,ny,nz].
//!*
//!* Conservative form of the Euler equations:
//!*
//!*     dU/dt + dF/dx + dG/dy + dH/dz = 0
//!*
//!* The normal flux is defined by
//!*
//!*     Fn = F*nx + G*ny + H*nz = | rho*qn          |
//!*                               | rho*qn*u + p*nx |
//!*                               | rho*qn*v + p*ny |
//!*                               | rho*qn*w + p*nz |
//!*                               | rho*qn*H        |    (qn = u*nx + v*ny + w*nz)
//!*
//!* The Roe flux is given by
//!*
//!*   Fn_Roe = 1/2 [ Fn(UR) + Fn(UL) - |An|dU ],
//!*
//!*  where
//!*
//!*    An = dFn/dU,  |An| = R|Lambda|L, dU = UR - UL.
//!*
//!* The dissipation term, |An|dU, is actually computed as
//!*
//!*     sum_{k=1,4} |lambda_k| * (LdU)_k * r_k,
//!*
//!* where lambda_k is the k-th eigenvalue, (LdU)_k is the k-th wave strength,
//!* and r_k is the k-th right-eigenvector evaluated at the Roe-average state.
//!*
//!* Note: The 4th component is a combined contribution from 2 shear waves.
//!*       They are combined to eliminate the tangent vectors.
//!*       So, (LdU)_4 is not really a wave strength, and
//!*       r_4 is not really an eigenvector.
//!*       See "I do like CFD, VOL.1" about how tangent vectors are eliminated.
//!*
//!* Note: In the code, the vector of conserative variables are denoted by uc.
//!*
//!* ------------------------------------------------------------------------------
//!*  Input: ucL(1:5) =  Left state (rhoL, rhoL*uL, rhoL*vL, rhoL*wL, rhoL*EL)
//!*         ucR(1:5) = Right state (rhoR, rhoR*uR, rhoR*vR, rhoR*wR, rhoR*ER)
//!*         njk(1:3) = unit face normal vector (nx, ny, nz), pointing from Left to Right.
//!*
//!*           njk
//!*  Face normal ^   o Right data point
//!*              |  .
//!*              | .
//!*              |.
//!*       -------x-------- Face
//!*             .                 Left and right states are
//!*            .                   1. Values at data points for 1st-order accuracy
//!*           .                    2. Extrapolated values at the face midpoint 'x'
//!*          o Left data point        for 2nd/higher-order accuracy.
//!*
//!*
//!* Output: dFnducL: Derivative of the Roe flux w.r.t. the left state ucL
//!*         dFnducR: Derivative of the Roe flux w.r.t. the left state ucR
//!* ------------------------------------------------------------------------------
//!*
//!* Katate Masatsuka, November 2012. http://www.cfdbooks.com
//!********************************************************************************


