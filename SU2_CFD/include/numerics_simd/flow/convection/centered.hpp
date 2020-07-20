/*!
 * \file centered.hpp
 * \brief Centered convective schemes.
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
#include "common.hpp"
#include "../../../variables/CEulerVariable.hpp"
#include "../../../../../Common/include/geometry/CGeometry.hpp"

/*!
 * \class CCenteredBase
 * \brief Base class for Centered schemes, derived classes implement
 * the dissipation term in a const "finalizeFlux" method.
 * \note See CRoeBase for the role of Base.
 */
template<class Derived, class Base>
class CCenteredBase : public Base {
protected:
  using Base::nDim;
  enum: size_t {nVar = CCompressibleConservatives<nDim>::nVar};
  enum: size_t {nPrimVar = Base::IsViscous? Base::nPrimVar : nDim+5};

  const su2double gamma;
  const su2double fixFactor;
  const bool dynamicGrid;
  static constexpr passivedouble stretchParam = 0.3;

  /*!
   * \brief Constructor, store some constants and forward args to base.
   */
  template<class... Ts>
  CCenteredBase(const CConfig& config, Ts&... args) : Base(config, args...),
    gamma(config.GetGamma()),
    fixFactor(config.GetCent_Jac_Fix_Factor()),
    dynamicGrid(config.GetDynamic_Grid()) {
  }

  /*!
   * \brief Special treatment needed to fetch integer data.
   */
  template<class T, size_t N>
  FORCEINLINE static Double numNeighbor(simd::Array<T,N> idx, const CGeometry& geometry) {
    Double n;
    for (size_t k=0; k<N; ++k) n[k] = geometry.nodes->GetnNeighbor(idx[k]);
    return n;
  }
  FORCEINLINE static Double numNeighbor(unsigned long idx, const CGeometry& geometry) {
    return geometry.nodes->GetnNeighbor(idx);
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

    const auto normal = gatherVariables<nDim>(iEdge, geometry.edges->GetNormal());
    const auto area = norm(normal);
    VectorDbl<nDim> unitNormal;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      unitNormal(iDim) = normal(iDim) / area;
    }

    /*--- Primitive variables. ---*/

    CPair<CCompressiblePrimitives<nDim,nPrimVar> > V;
    V.i.all = gatherVariables<nPrimVar>(iPoint, solution.GetPrimitive());
    V.j.all = gatherVariables<nPrimVar>(jPoint, solution.GetPrimitive());

    CCompressiblePrimitives<nDim,nPrimVar> avgV;
    for (size_t iVar = 0; iVar < nPrimVar; ++iVar) {
      avgV.all(iVar) = 0.5 * (V.i.all(iVar) + V.j.all(iVar));
    }

    /*--- Compute conservative variables. ---*/

    CPair<CCompressibleConservatives<nDim> > U;
    U.i = compressibleConservatives(V.i);
    U.j = compressibleConservatives(V.j);

    auto avgU = compressibleConservatives(avgV);

    VectorDbl<nVar> diffU;
    for (size_t iVar = 0; iVar < nVar-1; ++iVar) {
      diffU(iVar) = U.i.all(iVar) - U.j.all(iVar);
    }
    diffU(nVar-1) = V.i.density()*V.i.enthalpy() - V.j.density()*V.j.enthalpy();

    /*--- Inviscid fluxes and Jacobians. ---*/

    auto flux = inviscidProjFlux(avgV, avgU, normal);

    MatrixDbl<nVar> jac_i, jac_j;
    if (implicit) {
      jac_i = inviscidProjJac(gamma, avgV.velocity(), avgU.energy(), normal, 0.5);
      jac_j = jac_i;
    }

    /*--- Grid motion. ---*/

    Double projGridVel = 0.0;
    if (dynamicGrid) {
      const auto& gridVel = geometry.nodes->GetGridVel();
      projGridVel = 0.5*(dot(gatherVariables<nDim>(iPoint,gridVel), normal)+
                         dot(gatherVariables<nDim>(jPoint,gridVel), normal));

      for (size_t iVar = 0; iVar < nVar; ++iVar) {
        flux(iVar) -= projGridVel * avgU.all(iVar);
        if (implicit) {
          jac_i(iVar,iVar) -= 0.5 * projGridVel;
          jac_j(iVar,iVar) -= 0.5 * projGridVel;
        }
      }
    }

    /*--- Compute the local spectral radius and the stretching factor. ---*/

    const Double projVel_i = dot(V.i.velocity(), normal) - projGridVel;
    const Double projVel_j = dot(V.j.velocity(), normal) - projGridVel;

    const Double locLambda_i = abs(projVel_i) + V.i.speedSound()*area;
    const Double locLambda_j = abs(projVel_j) + V.j.speedSound()*area;
    const Double avgLambda = 0.5 * (locLambda_i + locLambda_j);

    const auto lambda_i = gatherVariables(iPoint, solution.GetLambda());
    const auto lambda_j = gatherVariables(jPoint, solution.GetLambda());
    const Double phi_i = pow(0.25*lambda_i/avgLambda, stretchParam);
    const Double phi_j = pow(0.25*lambda_j/avgLambda, stretchParam);
    const Double stretchFactor = 4*phi_i*phi_j / (phi_i + phi_j) * avgLambda;

    /*--- Finalize in derived class (static polymorphism). ---*/

    const auto derived = static_cast<const Derived*>(this);

    derived->finalizeFlux(flux, jac_i, jac_j, implicit, stretchFactor,
                          V, diffU, iPoint, jPoint, geometry, solution);

    /*--- Add the contributions from the base class (static decorator). ---*/

    Base::updateFlux(iEdge, iPoint, jPoint, V, avgV, solution_, geometry,
                     config, area, unitNormal, implicit, flux, jac_i, jac_j);

    /*--- Stop preaccumulation. ---*/

    AD::SetPreaccOut(flux, nVar, Double::Size);
    AD::EndPreacc();

    /*--- Update the vector and system matrix. ---*/

    updateLinearSystem(iEdge, iPoint, jPoint, implicit, updateType,
                       updateMask, flux, jac_i, jac_j, vector, matrix);
  }
};

/*!
 * \class CJSTScheme
 * \brief Classical JST scheme with scalar dissipation.
 */
template<class Decorator>
class CJSTScheme : public CCenteredBase<CJSTScheme<Decorator>,Decorator> {
private:
  using Base = CCenteredBase<CJSTScheme<Decorator>,Decorator>;
  using Base::nDim;
  using Base::nVar;
  using Base::gamma;
  using Base::fixFactor;
  const su2double kappa2;
  const su2double kappa4;

public:
  /*!
   * \brief Constructor, forward everything to base.
   */
  template<class... Ts>
  CJSTScheme(const CConfig& config, Ts&... args) : Base(config, args...),
    kappa2(config.GetKappa_2nd_Flow()),
    kappa4(config.GetKappa_4th_Flow()) {
  }

  /*!
   * \brief Updates flux and Jacobians with JST dissipation.
   * \note "Ts" is here just in case other schemes in the family need extra args.
   */
  template<class PrimVarType, class... Ts>
  FORCEINLINE void finalizeFlux(VectorDbl<nVar>& flux,
                                MatrixDbl<nVar>& jac_i,
                                MatrixDbl<nVar>& jac_j,
                                bool implicit,
                                Double stretchFactor,
                                const CPair<PrimVarType>& V,
                                const VectorDbl<nVar>& diffU,
                                Int iPoint,
                                Int jPoint,
                                const CGeometry& geometry,
                                const CEulerVariable& solution,
                                Ts&...) const {

    /*--- Compute dissipation coefficients. ---*/

    const auto ni = Base::numNeighbor(iPoint, geometry);
    const auto nj = Base::numNeighbor(jPoint, geometry);
    const Double sc2 = 3 * (ni+nj) / (ni*nj);
    const Double sc4 = 0.25*pow(sc2, 2);

    const auto si = gatherVariables(iPoint, solution.GetSensor());
    const auto sj = gatherVariables(jPoint, solution.GetSensor());
    const Double eps2 = kappa2 * 0.5*(si+sj) * sc2;
    const Double eps4 = max(0.0, kappa4-eps2) * sc4;

    /*--- Update flux and Jacobians with dissipation terms. ---*/

    const auto lapl_i = gatherVariables<nVar>(iPoint, solution.GetUndivided_Laplacian());
    const auto lapl_j = gatherVariables<nVar>(jPoint, solution.GetUndivided_Laplacian());

    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      flux(iVar) += (eps2*diffU(iVar) - eps4*(lapl_i(iVar)-lapl_j(iVar))) * stretchFactor;
    }

    if (implicit) {
      const Double dissip_i = fixFactor*(eps2 + eps4*(ni+1)) * stretchFactor;
      const Double dissip_j = -fixFactor*(eps2 + eps4*(nj+1)) * stretchFactor;
      scalarDissipationJacobian(V.i, gamma, dissip_i, jac_i);
      scalarDissipationJacobian(V.j, gamma, dissip_j, jac_j);
    }
  }
};

/*!
 * \class CJSTkeScheme
 * \brief JST scheme without 4th order dissipation.
 */
template<class Decorator>
class CJSTkeScheme : public CCenteredBase<CJSTkeScheme<Decorator>,Decorator> {
private:
  using Base = CCenteredBase<CJSTkeScheme<Decorator>,Decorator>;
  using Base::nDim;
  using Base::nVar;
  using Base::gamma;
  using Base::fixFactor;
  const su2double kappa2;

public:
  /*!
   * \brief Constructor, forward everything to base.
   */
  template<class... Ts>
  CJSTkeScheme(const CConfig& config, Ts&... args) : Base(config, args...),
    kappa2(config.GetKappa_2nd_Flow()) {
  }

  /*!
   * \brief Updates flux and Jacobians with 2nd order dissipation.
   * \note "Ts" is here just in case other schemes in the family need extra args.
   */
  template<class PrimVarType, class... Ts>
  FORCEINLINE void finalizeFlux(VectorDbl<nVar>& flux,
                                MatrixDbl<nVar>& jac_i,
                                MatrixDbl<nVar>& jac_j,
                                bool implicit,
                                Double stretchFactor,
                                const CPair<PrimVarType>& V,
                                const VectorDbl<nVar>& diffU,
                                Int iPoint,
                                Int jPoint,
                                const CGeometry& geometry,
                                const CEulerVariable& solution,
                                Ts&...) const {

    /*--- Compute dissipation coefficient. ---*/

    const auto ni = Base::numNeighbor(iPoint, geometry);
    const auto nj = Base::numNeighbor(jPoint, geometry);
    const Double sc2 = 3 * (ni+nj) / (ni*nj);

    const auto si = gatherVariables(iPoint, solution.GetSensor());
    const auto sj = gatherVariables(jPoint, solution.GetSensor());
    const Double dissip = kappa2 * 0.5*(si+sj) * sc2 * stretchFactor;

    /*--- Update flux and Jacobians with dissipation term. ---*/

    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      flux(iVar) += dissip * diffU(iVar);
    }

    if (implicit) {
      scalarDissipationJacobian(V.i, gamma, fixFactor*dissip, jac_i);
      scalarDissipationJacobian(V.j, gamma, -fixFactor*dissip, jac_j);
    }
  }
};

/*!
 * \class CLaxScheme
 * \brief Lax–Friedrichs 1st order scheme.
 */
template<class Decorator>
class CLaxScheme : public CCenteredBase<CLaxScheme<Decorator>,Decorator> {
private:
  using Base = CCenteredBase<CLaxScheme<Decorator>,Decorator>;
  using Base::nDim;
  using Base::nVar;
  using Base::gamma;
  using Base::fixFactor;
  const su2double kappa0;

public:
  /*!
   * \brief Constructor, forward everything to base.
   */
  template<class... Ts>
  CLaxScheme(const CConfig& config, Ts&... args) : Base(config, args...),
    kappa0(config.GetKappa_1st_Flow()) {
  }

  /*!
   * \brief Updates flux and Jacobians with 1st order scalar dissipation.
   * \note "Ts" is here just in case other schemes in the family need extra args.
   */
  template<class PrimVarType, class... Ts>
  FORCEINLINE void finalizeFlux(VectorDbl<nVar>& flux,
                                MatrixDbl<nVar>& jac_i,
                                MatrixDbl<nVar>& jac_j,
                                bool implicit,
                                Double stretchFactor,
                                const CPair<PrimVarType>& V,
                                const VectorDbl<nVar>& diffU,
                                Int iPoint,
                                Int jPoint,
                                const CGeometry& geometry,
                                Ts&...) const {

    /*--- Compute dissipation coefficient. ---*/

    const auto ni = Base::numNeighbor(iPoint, geometry);
    const auto nj = Base::numNeighbor(jPoint, geometry);
    const Double dissip = kappa0 * nDim * (ni+nj) / (ni*nj) * stretchFactor;

    /*--- Update flux and Jacobians with dissipation term. ---*/

    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      flux(iVar) += dissip * diffU(iVar);
    }

    if (implicit) {
      scalarDissipationJacobian(V.i, gamma, fixFactor*dissip, jac_i);
      scalarDissipationJacobian(V.j, gamma, -fixFactor*dissip, jac_j);
    }
  }
};
