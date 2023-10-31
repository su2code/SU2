/*!
 * \file viscous_fluxes.hpp
 * \brief Decorator classes for computation of viscous fluxes.
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
#include "common.hpp"

/*!
 * \class CNoViscousFlux
 * \ingroup ViscDiscr
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
 * \class CCompressibleViscousFluxBase
 * \ingroup ViscDiscr
 * \brief Decorator class to add viscous fluxes (compressible flow).
 */
template<size_t NDIM, class Derived>
class CCompressibleViscousFluxBase : public CNumericsSIMD {
protected:
  static constexpr size_t nDim = NDIM;
  static constexpr size_t nPrimVarGrad = nDim+1;

  const su2double gamma;
  const su2double gasConst;
  const su2double prandtlLam;
  const su2double prandtlTurb;
  const su2double cp;
  const bool correct;
  const bool useSA_QCR;
  const bool wallFun;
  const bool uq;
  const bool uq_permute;
  const size_t uq_eigval_comp;
  const su2double uq_delta_b;
  const su2double uq_urlx;

  const CVariable* turbVars;

  /*!
   * \brief Constructor, initialize constants and booleans.
   */
  template<class... Ts>
  CCompressibleViscousFluxBase(const CConfig& config, int iMesh,
                               const CVariable* turbVars_, Ts&...) :
    gamma(config.GetGamma()),
    gasConst(config.GetGas_ConstantND()),
    prandtlLam(config.GetPrandtl_Lam()),
    prandtlTurb(config.GetPrandtl_Turb()),
    cp(gamma * gasConst / (gamma - 1)),
    correct(iMesh == MESH_0),
    useSA_QCR(config.GetSAParsedOptions().qcr2000),
    wallFun(config.GetWall_Functions()),
    uq(config.GetSSTParsedOptions().uq),
    uq_permute(config.GetUQ_Permute()),
    uq_eigval_comp(config.GetEig_Val_Comp()),
    uq_delta_b(config.GetUQ_Delta_B()),
    uq_urlx(config.GetUQ_URLX()),
    turbVars(turbVars_) {
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

    static_assert(PrimVarType::nVar <= Derived::nPrimVar,"");

    /*--- Pointer on which to call the "compile-time virtual" methods. ---*/

    const auto derived = static_cast<const Derived*>(this);

    const auto& solution = static_cast<const CNSVariable&>(solution_);
    const auto& gradient = solution.GetGradient_Primitive();

    /*--- Compute distance and handle zero without "ifs" by making it large. ---*/

    auto dist2_ij = squaredNorm(vector_ij);
    Double mask = dist2_ij < EPS*EPS;
    dist2_ij += mask / (EPS*EPS);

    /*--- Compute the corrected mean gradient. ---*/

    auto avgGrad = averageGradient<nPrimVarGrad,nDim>(iPoint, jPoint, gradient);
    if(correct) correctGradient(V, vector_ij, dist2_ij, avgGrad);

    /*--- Stress and heat flux tensors. ---*/

    auto tau = stressTensor(avgV.laminarVisc() + (uq? Double(0.0) : avgV.eddyVisc()), avgGrad);
    if(useSA_QCR) addQCR(avgGrad, tau);
    if(uq) {
      Double turb_ke = 0.5*(gatherVariables(iPoint, turbVars->GetSolution()) +
                            gatherVariables(jPoint, turbVars->GetSolution()));
      addPerturbedRSM(avgV, avgGrad, turb_ke, tau,
                      uq_eigval_comp, uq_permute, uq_delta_b, uq_urlx);
    }

    if(wallFun) addTauWall(iPoint, jPoint, solution.GetTau_Wall(), unitNormal, tau);

    Double cond = derived->thermalConductivity(avgV);
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

    /*--- Energy flux Jacobian. ---*/
    auto dEdU = derived->energyJacobian(avgV, dtau, cond, area, dist_ij, iPoint, jPoint, solution);

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

/*!
 * \class CCompressibleViscousFlux
 * \ingroup ViscDiscr
 * \brief Decorator class to add viscous fluxes (compressible flow, ideal gas).
 */
template<size_t NDIM>
class CCompressibleViscousFlux : public CCompressibleViscousFluxBase<NDIM, CCompressibleViscousFlux<NDIM> > {
public:
  static constexpr size_t nPrimVar = NDIM+7;
  using Base = CCompressibleViscousFluxBase<NDIM, CCompressibleViscousFlux<NDIM> >;
  using Base::gamma;
  using Base::gasConst;
  using Base::prandtlLam;
  using Base::prandtlTurb;
  using Base::cp;

  /*!
   * \brief Constructor, initialize constants and booleans.
   */
  template<class... Ts>
  CCompressibleViscousFlux(Ts&... args) : Base(args...) {}

  /*!
   * \brief Compute the thermal conductivity.
   */
  template<class PrimitiveType>
  FORCEINLINE Double thermalConductivity(const PrimitiveType& V) const {
    return cp * (V.laminarVisc()/prandtlLam + V.eddyVisc()/prandtlTurb);
  }

  /*!
   * \brief Compute Jacobian of the energy flux, except the part due to the work of viscous forces.
   */
  template<size_t nVar, size_t nDim, class PrimitiveType, class... Ts>
  FORCEINLINE VectorDbl<nVar> energyJacobian(const PrimitiveType& V,
                                             const MatrixDbl<nDim,nVar>& dtau,
                                             Double thermalCond,
                                             Double area,
                                             Double dist_ij,
                                             Ts&... args) const {
    Double vel2 = 0.5 * squaredNorm<nDim>(V.velocity());
    Double phi = (gamma-1) / V.density();
    Double RdTdrho = phi*vel2 - V.pressure() / pow(V.density(),2);
    Double condOnRd = thermalCond / (gasConst * dist_ij);
    Double contraction = 0.0;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      contraction += dtau(iDim,0) * V.velocity(iDim);
    }
    VectorDbl<nVar> dEdU;
    dEdU(0) = area * (condOnRd * RdTdrho - contraction);
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      dEdU(iDim+1) = area * (dtau(iDim,0) - condOnRd*phi*V.velocity(iDim));
    }
    dEdU(nDim+1) = area * condOnRd * phi;

    return dEdU;
  }
};

/*!
 * \class CGeneralCompressibleViscousFlux
 * \ingroup ViscDiscr
 * \brief Decorator class to add viscous fluxes (compressible flow, real gas).
 */
template<size_t NDIM>
class CGeneralCompressibleViscousFlux : public CCompressibleViscousFluxBase<NDIM, CGeneralCompressibleViscousFlux<NDIM> > {
public:
  static constexpr size_t nPrimVar = NDIM+9;
  static constexpr size_t nSecVar = 4;
  using Base = CCompressibleViscousFluxBase<NDIM, CGeneralCompressibleViscousFlux<NDIM> >;
  using Base::prandtlTurb;

  /*!
   * \brief Constructor, initialize constants and booleans.
   */
  template<class... Ts>
  CGeneralCompressibleViscousFlux(Ts&... args) : Base(args...) {}

  /*!
   * \brief Compute the thermal conductivity.
   */
  template<class PrimitiveType>
  FORCEINLINE Double thermalConductivity(const PrimitiveType& V) const {
    return V.thermalCond() + V.cp()*V.eddyVisc()/prandtlTurb;
  }

  /*!
   * \brief Compute Jacobian of the energy flux, except the part due to the work of viscous forces.
   */
  template<size_t nVar, size_t nDim, class PrimitiveType, class VariableType>
  FORCEINLINE VectorDbl<nVar> energyJacobian(const PrimitiveType& V,
                                             const MatrixDbl<nDim,nVar>& dtau,
                                             Double thermalCond,
                                             Double area,
                                             Double dist_ij,
                                             Int iPoint,
                                             Int jPoint,
                                             const VariableType& solution) const {
    Double vel2 = squaredNorm<nDim>(V.velocity());
    Double contraction = 0.0;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      contraction += dtau(iDim,0) * V.velocity(iDim);
    }

    auto secVar_i = gatherVariables<nSecVar>(iPoint, solution.GetSecondary());
    auto secVar_j = gatherVariables<nSecVar>(jPoint, solution.GetSecondary());

    Double dTdrho_e = 0.5 * (secVar_i(2) + secVar_j(2));
    Double dTde_rho = 0.5 * (secVar_i(3) + secVar_j(3)) / V.density();

    Double condOnDist = thermalCond / dist_ij;
    Double dTdrho = dTdrho_e + dTde_rho*(vel2-V.enthalpy()+V.pressure()/V.density());

    VectorDbl<nVar> dEdU;
    dEdU(0) = area * (condOnDist * dTdrho - contraction);
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      dEdU(iDim+1) = area * (dtau(iDim,0) - condOnDist*dTde_rho*V.velocity(iDim));
    }
    dEdU(nDim+1) = area * condOnDist * dTde_rho;

    return dEdU;
  }
};
