/*!
 * \file upwind.hpp
 * \brief Roe-family of convective schemes.
 * \author P. Gomes, A. Bueno, F. Palacios
 * \version 8.4.0 "Harrier"
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
#include "common.hpp"
#include "../../../variables/CEulerVariable.hpp"
#include "../../../../../Common/include/geometry/CGeometry.hpp"

/*!
 * \class CUpwindBase
 * \ingroup ConvDiscr
 * \brief Base class for upwind schemes, derived classes implement
 * details of the convective flux in a const "finalizeFlux" method.
 * A base class implementing "viscousTerms" is accepted as template parameter.
 * Similarly to derived, that method should update the flux and Jacobians, but
 * whereas "finalizeFlux" is passed data prepared by CUpwindBase, "viscousTerms"
 * takes the same input arguments as "ComputeFlux", i.e. it can fetch more
 * data from CVariable. Derived is meant to implement small details,
 * Base is meant to do heavy lifting.
 */
template<class Derived, class Base>
class CUpwindBase : public Base {
protected:
  using Base::nDim;
  static constexpr size_t nVar = CCompressibleConservatives<nDim>::nVar;
  static constexpr size_t nPrimVarGrad = nDim+4;
  static constexpr size_t nPrimVar = Max(Base::nPrimVar, nPrimVarGrad);

  const su2double gamma;
  const su2double gasConst;
  const bool finestGrid;
  const bool dynamicGrid;
  const bool muscl;
  const su2double umusclKappa;
  const LIMITER typeLimiter;

  /*!
   * \brief Constructor, store some constants and forward args to base.
   */
  template<class... Ts>
  CUpwindBase(const CConfig& config, unsigned iMesh, Ts&... args) : Base(config, iMesh, args...),
    gamma(config.GetGamma()),
    gasConst(config.GetGas_ConstantND()),
    finestGrid(iMesh == MESH_0),
    dynamicGrid(config.GetDynamic_Grid()),
    muscl(finestGrid && config.GetMUSCL_Flow()),
    umusclKappa(config.GetMUSCL_Kappa_Flow()),
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
    const su2double umusclRamp = config.GetMUSCLRampValue() * config.GetNewtonKrylovRelaxation();
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

    /*--- Recompute density and enthalpy instead of reconstructing. ---*/
    auto V = reconstructPrimitives<CCompressiblePrimitives<nDim,nPrimVarGrad> >(
        iEdge, iPoint, jPoint, gamma, gasConst, muscl, umusclKappa, umusclRamp, typeLimiter, V1st, vector_ij, solution);

    /*--- Compute conservative variables. ---*/

    CPair<CCompressibleConservatives<nDim> > U;
    U.i = compressibleConservatives(V.i);
    U.j = compressibleConservatives(V.j);

    /*--- Finalize in derived class (static polymorphism). ---*/

    const auto derived = static_cast<const Derived*>(this);
    VectorDbl<nVar> flux;
    MatrixDbl<nVar> jac_i, jac_j;
    derived->finalizeFlux(flux, jac_i, jac_j, implicit, area, unitNormal,
                          normal, V, U, iPoint, jPoint, solution, geometry);

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
class CRoeScheme : public CUpwindBase<CRoeScheme<Decorator>, Decorator> {
private:
  using Base = CUpwindBase<CRoeScheme<Decorator>, Decorator>;
  using Base::nDim;
  using Base::nVar;
  using Base::nPrimVarGrad;
  using Base::nPrimVar;

  const su2double kappa;
  const su2double entropyFix;
  const ENUM_ROELOWDISS typeDissip;
  using Base::gamma;
  using Base::gasConst;
  using Base::dynamicGrid;

public:
  /*!
   * \brief Constructor, store some constants and forward to base.
   */
  template<class... Ts>
  CRoeScheme(const CConfig& config, Ts&... args) : Base(config, args...),
    kappa(config.GetRoe_Kappa()),
    entropyFix(config.GetEntropyFix_Coeff()),
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
                                const VectorDbl<nDim>& normal,
                                const CPair<PrimVarType>& V,
                                const CPair<ConsVarType>& U,
                                Int iPoint,
                                Int jPoint,
                                const CEulerVariable& solution,
                                const CGeometry& geometry,
                                Ts&...) const {
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

    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      flux(iVar) = 0.5 * (flux_i(iVar) + flux_j(iVar));
    }
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

    /*--- P tensor. ---*/

    auto pMat = pMatrix(gamma, roeAvg.density, roeAvg.velocity,
                        roeAvg.projVel, roeAvg.speedSound, unitNormal);

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

        if (implicit) {
          jac_i(iVar,jVar) += dDdU;
          jac_j(iVar,jVar) -= dDdU;
        }
      }
    }
  }
};

/*!
 * \class CMSWScheme
 * \ingroup ConvDiscr
 * \brief MSW scheme with switch to SW based on a pressure sensor.
 */
template<class Decorator>
class CMSWScheme : public CUpwindBase<CMSWScheme<Decorator>, Decorator> {
protected:
  using Base = CUpwindBase<CMSWScheme<Decorator>, Decorator>;
  using Base::nDim;
  using Base::nVar;
  using Base::nPrimVarGrad;
  using Base::nPrimVar;

  const su2double alpha;
  using Base::gamma;
  using Base::gasConst;
  using Base::dynamicGrid;

public:
  /*!
   * \brief Constructor, store some constants and forward args to base.
   */
  template<class... Ts>
  CMSWScheme(const CConfig& config, unsigned iMesh, Ts&... args) : Base(config, iMesh, args...),
    alpha(config.GetMSW_Alpha()) {
  }

  /*!
   * \brief Implementation of the flux.
   */
  template<class PrimVarType, class ConsVarType, class... Ts>
  FORCEINLINE void finalizeFlux(VectorDbl<nVar>& flux,
                                MatrixDbl<nVar>& jac_i,
                                MatrixDbl<nVar>& jac_j,
                                bool implicit,
                                Double area,
                                const VectorDbl<nDim>& unitNormal,
                                const VectorDbl<nDim>& normal,
                                const CPair<PrimVarType>& V,
                                const CPair<ConsVarType>& U,
                                Int iPoint,
                                Int jPoint,
                                const CEulerVariable& solution,
                                const CGeometry& geometry,
                                Ts&...) const {

    /*--- Weighted states for flux-vector spliting (see the non-SIMD version in fvs.cpp for notes). ---*/

    const auto si = gatherVariables(iPoint, solution.GetSensor());
    const auto sj = gatherVariables(jPoint, solution.GetSensor());

    const Double dp = fmax(si, sj) - alpha * 0.06;
    const Double w = 0.25 * (1 - sign(dp)) * (1 - exp(-100 * abs(dp)));
    const Double onemw = 1 - w;

    CPair<CCompressiblePrimitives<nDim, nPrimVarGrad>> Vweighted;
    for (size_t iVar = 0; iVar < nPrimVarGrad; ++iVar) {
      Vweighted.i.all(iVar) = onemw * V.i.all(iVar) + w * V.j.all(iVar);
      Vweighted.j.all(iVar) = onemw * V.j.all(iVar) + w * V.i.all(iVar);
    }
    const Double soundSpeed_i = sqrt(gamma * gasConst * Vweighted.i.temperature());
    const Double soundSpeed_j = sqrt(gamma * gasConst * Vweighted.j.temperature());

    Double projGridVel_i = 0.0, projGridVel_j = 0.0;
    if (dynamicGrid) {
      const auto& gridVel = geometry.nodes->GetGridVel();
      const Double vn_i = dot(gatherVariables<nDim>(iPoint, gridVel), unitNormal);
      const Double vn_j = dot(gatherVariables<nDim>(jPoint, gridVel), unitNormal);
      projGridVel_i = onemw * vn_i + w * vn_j;
      projGridVel_j = onemw * vn_j + w * vn_i;
    }
    Double projVel_i = dot(Vweighted.i.velocity(), unitNormal) - projGridVel_i;
    Double projVel_j = dot(Vweighted.j.velocity(), unitNormal) - projGridVel_j;

    /*--- Lambda+ ---*/

    VectorDbl<nVar> lambda;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      lambda(iDim) = fmax(projVel_i, 0);
    }
    lambda(nDim) = fmax(projVel_i + soundSpeed_i, 0);
    lambda(nDim+1) = fmax(projVel_i - soundSpeed_i, 0);

    auto pMat = pMatrix(gamma, Vweighted.i.density(), Vweighted.i.velocity(),
                        projVel_i, soundSpeed_i, unitNormal);
    auto pMatInv = pMatrixInv(gamma, Vweighted.i.density(), Vweighted.i.velocity(),
                              projVel_i, soundSpeed_i, unitNormal);

    auto updateFlux = [&](const auto& u, auto& jac) {
      for (size_t iVar = 0; iVar < nVar; ++iVar) {
        for (size_t jVar = 0; jVar < nVar; ++jVar) {
          /*--- Compute |projModJacTensor| = P x |Lambda| x P^-1. ---*/

          Double projModJacTensor = 0.0;
          for (size_t kVar = 0; kVar < nVar; ++kVar) {
            projModJacTensor += pMat(iVar,kVar) * lambda(kVar) * pMatInv(kVar,jVar);
          }
          Double dFdU = projModJacTensor * area;

          flux(iVar) += dFdU * u(jVar);
          jac(iVar, jVar) = dFdU;
        }
      }
    };
    updateFlux(U.i.all, jac_i);

    /*--- Lambda- ---*/

    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      lambda(iDim) = fmin(projVel_j, 0);
    }
    lambda(nDim) = fmin(projVel_j + soundSpeed_j, 0);
    lambda(nDim+1) = fmin(projVel_j - soundSpeed_j, 0);

    pMat = pMatrix(gamma, Vweighted.j.density(), Vweighted.j.velocity(),
                   projVel_j, soundSpeed_j, unitNormal);
    pMatInv = pMatrixInv(gamma, Vweighted.j.density(), Vweighted.j.velocity(),
                         projVel_j, soundSpeed_j, unitNormal);
    updateFlux(U.j.all, jac_j);
  }
};
