/*!
 * \file viscous_fluxes.hpp
 * \brief Decorator classes for computation of viscous fluxes.
 * \author P. Gomes, C. Pederson, A. Bueno, F. Palacios, T. Economon
 * \version 7.0.8 "Blackbird"
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

/*!
 * \class CNoViscousFlux
 * \brief Numerics classes that accept a compile-time decorator should use this
 * class template as a "do-nothing" decorator and as a link to the interface when
 * they are not being decorated.
 * Compile-time decoration works by specifying the base class as a template parameter.
 * Then the class being decorated should call a method of its base class to add some
 * contribution to the flux/source and Jacobians just before writting the results to
 * CSysVector and CSysMatrix. The mechanism can be used to chain any number of classes
 * at compile-time, but its main purpose is to combine convective and viscous fluxes
 * in the most (nearly) efficient way.
 */
template<size_t NDIM>
class CNoViscousFlux : public CNumericsSIMD {
protected:
  static constexpr size_t nDim = NDIM;
  static constexpr size_t nPrimVar = 0;

  template<class... Ts>
  CNoViscousFlux(Ts&...) {}

  /*!
   * \brief Empty method, real decorators should take as arguments whatever
   * the decorated class can pass them to avoid expensive data accesses.
   */
  template<class... Ts>
  void viscousTerms(Ts&...) const {}
};

/*!
 * \class CCompressibleViscousFlux
 * \brief Decorator class to add viscous fluxes (compressible flow, ideal gas).
 */
template<size_t NDIM>
class CCompressibleViscousFlux : public CNumericsSIMD {
protected:
  static constexpr size_t nDim = NDIM;
  static constexpr size_t nPrimVar = nDim+7;
  static constexpr size_t nPrimVarGrad = nDim+1;

  const su2double gamma;
  const su2double gasConst;
  const su2double prandtlLam;
  const su2double prandtlTurb;
  const su2double cp;
  const bool correct;
  const bool useSA_QCR;

  /*!
   * \brief Constructor, initialize constants and booleans.
   */
  template<class... Ts>
  CCompressibleViscousFlux(const CConfig& config, int iMesh, Ts&...) :
    gamma(config.GetGamma()),
    gasConst(config.GetGas_ConstantND()),
    prandtlLam(config.GetPrandtl_Lam()),
    prandtlTurb(config.GetPrandtl_Turb()),
    cp(gamma * gasConst / (gamma - 1)),
    correct(iMesh == MESH_0),
    useSA_QCR(config.GetQCR()) {
  }

  /*!
   * \brief Add viscous contributions to flux and jacobians.
   */
  template<class PrimVarType, size_t nVar>
  FORCEINLINE void viscousTerms(Int iEdge,
                                Int iPoint,
                                Int jPoint,
                                const PrimVarType& avgV,
                                const CPair<PrimVarType>& V,
                                const CVariable& solution_,
                                const VectorDbl<nDim>& vector_ij,
                                const CGeometry& geometry,
                                const CConfig& config,
                                Double area,
                                const VectorDbl<nDim>& unitNormal,
                                bool implicit,
                                VectorDbl<nVar>& flux,
                                MatrixDbl<nVar>& jac_i,
                                MatrixDbl<nVar>& jac_j) const {

    static_assert(PrimVarType::nVar <= nPrimVar,"");

    const auto& solution = static_cast<const CNSVariable&>(solution_);
    const auto& gradient = solution.GetGradient_Primitive();

    /*--- Compute distance and handle zero without "ifs" by making it large. ---*/

    auto dist2_ij = squaredNorm(vector_ij);
    Double mask = dist2_ij < EPS*EPS;
    dist2_ij += mask / (EPS*EPS);

    /*--- Compute the corrected mean gradient. ---*/

    auto avgGrad = averageGradient<nPrimVarGrad,nDim>(iPoint, jPoint, gradient);
    if(correct) correctGradient(V, vector_ij, dist2_ij, avgGrad);

    /// TODO: Uncertainty quantification (needs a way to access tke, maybe in ctor).

    /*--- Stress and heat flux tensors. ---*/

    auto tau = stressTensor(avgV, avgGrad);
    if(useSA_QCR) addQCR(avgGrad, tau);

    Double cond = cp * (avgV.laminarVisc()/prandtlLam + avgV.eddyVisc()/prandtlTurb);
    VectorDbl<nDim> heatFlux;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      heatFlux(iDim) = cond * avgGrad(0,iDim);
    }

    /*--- Projected flux. ---*/

    auto viscFlux = viscousFlux<nVar>(avgV, tau, heatFlux, unitNormal);
    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      viscFlux(iVar) *= area;
      flux(iVar) -= viscFlux(iVar);
    }

    if (!implicit) return;

    /*--- Flux Jacobians. ---*/

    Double dist_ij = sqrt(dist2_ij);
    auto dtau = stressTensorJacobian<nVar>(avgV, unitNormal, dist_ij);
    Double contraction = 0.0;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      contraction += dtau(iDim,0) * avgV.velocity(iDim);
    }

    /*--- Energy flux Jacobian. ---*/
    VectorDbl<nVar> dEdU;
    Double vel2 = 0.5 * squaredNorm<nDim>(avgV.velocity());
    Double phi = (gamma-1) / avgV.density();
    Double RdTdrho = phi*vel2 - avgV.pressure() / pow(avgV.density(),2);
    Double condOnRd = cond / (gasConst * dist_ij);

    dEdU(0) = area * (condOnRd * RdTdrho - contraction);
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      dEdU(iDim+1) = area * (condOnRd*phi*avgV.velocity(iDim) + dtau(iDim,0));
    }
    dEdU(nDim+1) = area * condOnRd * phi;

    /*--- Update momentum and energy terms ("symmetric" part). ---*/
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      for (size_t iVar = 0; iVar < nVar; ++iVar) {
        jac_i(iDim+1,iVar) -= area * dtau(iDim,iVar);
        jac_j(iDim+1,iVar) += area * dtau(iDim,iVar);
      }
    }
    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      jac_i(nDim+1,iVar) += dEdU(iVar);
      jac_j(nDim+1,iVar) -= dEdU(iVar);
    }
    /*--- "Non-symmetric" energy terms. ---*/
    Double proj = dot<nDim>(&viscFlux(1), avgV.velocity());
    Double halfOnRho = 0.5/avgV.density();
    jac_i(nDim+1,0) += halfOnRho * proj;
    jac_j(nDim+1,0) += halfOnRho * proj;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      jac_i(nDim+1,iDim+1) -= halfOnRho * viscFlux(iDim+1);
      jac_j(nDim+1,iDim+1) -= halfOnRho * viscFlux(iDim+1);
    }
  }

  /*!
   * \overload Average primitives if not provided yet.
   */
  template<class PrimVarType, class... Ts>
  FORCEINLINE void viscousTerms(Int iEdge,
                                Int iPoint,
                                Int jPoint,
                                const CPair<PrimVarType>& V,
                                Ts&... args) const {
    PrimVarType avgV;
    for (size_t iVar = 0; iVar < PrimVarType::nVar; ++iVar) {
      avgV.all(iVar) = 0.5 * (V.i.all(iVar) + V.j.all(iVar));
    }

    /*--- Continue calculation. ---*/
    viscousTerms(iEdge, iPoint, jPoint, avgV, V, args...);
  }

  /*!
   * \overload Compute the i-j vector if not provided yet.
   */
  template<class PrimVarType, class... Ts>
  FORCEINLINE void viscousTerms(Int iEdge,
                                Int iPoint,
                                Int jPoint,
                                const PrimVarType& avgV,
                                const CPair<PrimVarType>& V,
                                const CVariable& solution_,
                                const CGeometry& geometry,
                                Ts&... args) const {

    const auto vector_ij = distanceVector<nDim>(iPoint, jPoint, geometry.nodes->GetCoord());

    /*--- Continue calculation. ---*/
    viscousTerms(iEdge, iPoint, jPoint, avgV, V, solution_, vector_ij, geometry, args...);
  }
};
