/*!
 * \file roe.hpp
 * \brief Roe-family of convective schemes.
 * \author P. Gomes, A. Bueno, F. Palacios
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
