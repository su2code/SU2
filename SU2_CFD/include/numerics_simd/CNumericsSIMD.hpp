/*!
 * \file CNumericsSIMD.hpp
 * \brief Vectorized (SIMD) numerics classes.
 * \author P. Gomes
 * \version 7.0.5 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is masize_tained by the SU2 Foundation
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

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/parallelization/vectorization.hpp"
#include "../../../Common/include/toolboxes/C2DContainer.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

/*!
 * \enum UpdateType
 * \brief Ways to update vectors and system matrices.
 * COLORING is the typical i/j update, whereas for REDUCTION
 * the fluxes are stored and the matrix diagonal is not modified.
 */
enum class UpdateType {COLORING, REDUCTION};

/*--- Scalar types (SIMD). ---*/
using Double = simd::Array<su2double>
using Int = simd::Array<unsigned long, Double::Size>

/*--- Static vector types. ---*/
template<class Type, size_t Size>
using Vector<Type,Size> = C2DContainer<unsigned long, Type, StorageType::ColumnMajor, Type::Align, Size, 1>;

template<size_t Size> using VectorInt<Size> = Vector<Int, Size>;
template<size_t Size> using VectorDbl<Size> = Vector<Double, Size>;

/*--- Static matrix types. ---*/
template<class Type, size_t Rows, size_t Cols>
using Matrix<Type,Rows,Cols> = C2DContainer<unsigned long, Type, StorageType::RowMajor, Type::Align, Rows, Cols>;

template<size_t Rows, size_t Cols = Rows> using MatrixInt<Rows,Cols> = Matrix<Int, Rows, Cols>;
template<size_t Rows, size_t Cols = Rows> using MatrixDbl<Rows,Cols> = Matrix<Double, Rows, Cols>;

/*!
 * \class CNumericsSIMD
 * \brief Base class to define the size_terface.
 */
class CNumericsSIMD {
public:
  /*!
   * \brief Interface for edge flux computation.
   * \param[in] iEdge - The starting edge for flux computation.
   * \param[in] config - Problem definitions.
   * \param[in] geometry - Problem geometry.
   * \param[in] solution - Solution variables.
   * \param[in] updateType - Type of update done on vector and matrix.
   * \param[in,out] vector - Target for the fluxes.
   * \param[in,out] matrix - Target for the flux Jacobians.
   */
  virtual void ComputeFlux(unsigned long iEdge,
                           const CConfig& config,
                           const CGeometry& geometry,
                           const CVariable& solution,
                           UpdateType updateType,
                           CSysVector& vector,
                           CSysMatrix& matrix) const = 0;

  /*! \brief Destructor of the class. */
  virtual ~CNumerics(void) = default;
};

template<size_t nDim>
FORCEINLINE Double dot(const VectorDbl<nDim>& a,
                       const VectorDbl<nDim>& b,) {
  Double sum = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    sum += a(iDim) * b(iDim);
  }
  return sum;
}

template<size_t nDim, class ForwardIterator>
FORCEINLINE Double dot(ForwardIterator iterator,
                       const VectorDbl<nDim>& vector) {
  Double sum = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    sum += *(iterator++) * vector(iDim);
  }
  return sum;
}

template<size_t nDim>
FORCEINLINE Double squaredNorm(const VectorDbl<nDim>& vector) {
  Double sum = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    sum += pow(vector(iDim),2);
  }
  return sum;
}

template<size_t nDim, class ForwardIterator>
FORCEINLINE Double squaredNorm(ForwardIterator iterator) {
  Double sum = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    sum += pow(*(iterator++),2);
  }
  return sum;
}

template<size_t nDim>
FORCEINLINE Double norm(const VectorDbl<nDim>& vector) { return sqrt(squaredNorm(vector)); }

template<size_t nDim, class Coord_t>
FORCEINLINE VectorDbl<nDim> distanceVector(Int iPosize_t, Int jPosize_t, const Coord_t& coords) {
  VectorDbl<nDim> vector_ij;
  auto coord_i = coords.innerIter<nDim>(iPosize_t);
  auto coord_j = coords.innerIter<nDim>(jPosize_t);
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    vector_ij(iDim) = 0.5 * (*(coord_j++) - *(coord_i++));
  }
  return vector_ij;
}

template<size_t nVar, template Int_t, class Field_t>
FORCEINLINE VectorDbl<nVar> gatherVariables(Int_t idx, const Field_t& vars) {
  VectorDbl<nVar> v;
  auto it = vars.innerIter<nVar>(idx);
  while (it.go()) v(iVar) = *(it++);
  return v;
}

template<size_t nVar, size_t nDim, class Field_t, class Gradient_t>
FORCEINLINE VectorDbl<nVar> musclUnlimited(Int iPosize_t,
                                           const VectorDbl<nDim>& vector_ij,
                                           const Field_t& field,
                                           const Gradient_t& gradient) {
  auto vars = gatherVariables<nVar>(iPosize_t, field);
  for (size_t iVar = 0; iVar < nVar; ++iVar) {
    vars(iVar) += dot(gradient.innerIter<nVar,nDim>(iPosize_t,iVar), vector_ij);
  }
  return vars;
}

template<size_t nVar, size_t nDim, class Field_t, class Gradient_t>
FORCEINLINE VectorDbl<nVar> musclPosize_tLimited(Int iPosize_t,
                                              const VectorDbl<nDim>& vector_ij,
                                              const Field_t& field,
                                              const Field_t& limiter,
                                              const Gradient_t& gradient) {
  auto vars = gatherVariables<nVar>(iPosize_t, field);
  auto itLim = limiter.innerIter<nVar>(iPosize_t);
  for (size_t iVar = 0; iVar < nVar; ++iVar) {
    vars(iVar) += *(itLim++) * dot(gradient.innerIter<nVar,nDim>(iPosize_t,iVar), vector_ij);
  }
  return vars;
}

template<size_t nDim, size_t nVar, class Field_t, class Gradient_t>
FORCEINLINE void musclEdgeLimited(Int iPosize_t,
                                  Int jPosize_t,
                                  const VectorDbl<nDim>& vector_ij,
                                  const Field_t& field,
                                  const Gradient_t& gradient,
                                  VectorDbl<nVar>& vars_i,
                                  VectorDbl<nVar>& vars_j) {
  vars_i = gatherVariables<nVar>(iPosize_t, field);
  vars_j = gatherVariables<nVar>(jPosize_t, field);

  for (size_t iVar = 0; iVar < nVar; ++iVar) {
    auto proj_i = dot(gradient.innerIter<nVar,nDim>(iPosize_t,iVar), vector_ij);
    auto proj_j = dot(gradient.innerIter<nVar,nDim>(jPosize_t,iVar), vector_ij);
    auto delta_ij = vars_j(iVar) - vars_i(iVar);
    auto delta_ij_2 = pow(delta_ij,2);
    /// TODO: Customize the limiter function
    auto lim_i = (delta_ij_2 + 2*proj_i*delta_ij) / (4*pow(proj_i,2) + delta_ij_2 + eps);
    auto lim_j = (delta_ij_2 + 2*proj_j*delta_ij) / (4*pow(proj_j,2) + delta_ij_2 + eps);
    vars_i(iVar) += lim_i * proj_i;
    vars_j(iVar) -= lim_j * proj_j;
  }
}

template<class T>
struct CPair {
  T i, j;
};

template<size_t nDim>
struct CCompressiblePrimitives {
  enum : size_t {nVar = nDim+4};
  VectorDbl<nVar> all;
  Double& temperature = all(0);
  Double& pressure = all(nDim+1);
  Double& density = all(nDim+2);
  Double& enthalpy = all(nDim+3);
  FORCEINLINE Double& velocity(size_t iDim) { return all(iDim+1); }
  FORCEINLINE Double* velocity() { return &velocity(0); }
};

template<size_t nDim, class Variable_t>
FORCEINLINE CPair<CCompressibleConservatives<nDim> > getCompressiblePrimitives(Int iPosize_t, Int jPosize_t, bool muscl,
                                                                               ENUM_LIMITER limiterType,
                                                                               const VectorDbl<nDim>& vector_ij,
                                                                               const Variable_t& solution) {
  CPair<CCompressiblePrimitives<nDim> > V;
  constexpr size_t nVar = V::nVar;

  const auto& primitives = solution.Primitive;
  const auto& gradients = solution.Gradient_Reconstruction;
  const auto& limiters = solution.Limiter_Primitive;

  if (!muscl) {
    V.i.all = gatherVariables<nVar>(iPosize_t, primitives);
    V.j.all = gatherVariables<nVar>(jPosize_t, primitives);
  }
  else {
    switch (limiterType) {
    case NO_LIMITER:
      V.i.all = musclUnlimited<nVar>(iPosize_t, vector_ij, primitives, gradients);
      V.j.all = musclUnlimited<nVar>(jPosize_t,-vector_ij, primitives, gradients);
    case VAN_ALBADA_EDGE:
      musclEdgeLimited(iPosize_t, jPosize_t, vector_ij, primitives, gradients, V.i.all, V.j.all);
    default:
      V.i.all = musclPosize_tLimited<nVar>(iPosize_t, vector_ij, primitives, limiters, gradients);
      V.j.all = musclPosize_tLimited<nVar>(jPosize_t,-vector_ij, primitives, limiters, gradients);
    }
    /// TODO: Extra reconstruction checks needed.
  }
  return V;
}

template<size_t nDim>
struct CCompressibleConservatives {
  enum : size_t {nVar = nDim+2};
  VectorDbl<nVar> all;
  Double& density = all(0);
  Double& rhoEnergy = all(nDim+1);
  FORCEINLINE Double& momentum(size_t iDim) { return all(iDim+1); }
  FORCEINLINE Double* momentum() { return &momentum(0); }
};

template<size_t nDim>
FORCEINLINE CCompressibleConservatives<nDim> getCompressibleConservatives(const CCompressiblePrimitives<nDim>& V) {
  CCompressibleConservatives<nDim> U;
  U.density = V.density;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    U.momentum(iDim) = V.density * V.velocity(iDim);
  }
  U.rhoEnergy = V.density * V.enthalpy - V.pressure;
}

template<size_t nDim>
struct CRoeVariables {
  Double density;
  VectorDbl<nDim> velocity;
  Double enthalpy;
  Double speedSound;
  Double projVel;
};

template<size_t nDim>
FORCEINLINE CRoeVariables<nDim> getRoeAveragedVariables(Double gamma,
                                                        const CPair<CCompressibleConservatives<nDim> >& V,
                                                        const VectorDbl<nDim>& normal) {
  CRoeVariables<nDim> roeAvg;
  auto R = sqrt(abs(V.j.density / V.i.density));
  roeAvg.density = R * V.i.density;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    roeAvg.velocity(iDim) = (R*V.j.velocity(iDim) + V.i.velocity(iDim)) / (R+1);
  }
  roeAvg.enthalpy = (R*V.j.enthalpy + V.i.enthalpy) / (R+1);
  roeAvg.soundSpeed = sqrt((gamma-1) * (roeAvg.enthalpy - 0.5*squaredNorm(roeAvg.velocity)));
  roeAvg.projVel = dot(roeAvg.velocity, normal);
  return roeAvg;
}

template<size_t nDim>
FORCEINLINE MatrixDbl<nDim+2> pMatrix(Double gamma, Double density, const VectorDbl<nDim>& velocity,
                                      Double projVel, Double soundSpeed, const VectorDbl<nDim>& normal) {
  MatrixDbl<nDim+2> pMat;

  auto oneOnC = 1/soundSpeed;
  auto vel2 = 0.5*squaredNorm(velocity);

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
    pMat(4,1) = vel2*normal(1) - density*(velocity(0)*normal(2) + velocity(2)*normal(0));
    pMat(4,2) = vel2*normal(2) + density*(velocity(0)*normal(1) - velocity(1)*normal(0));
  }

  /*--- Last two columns. ---*/

  pMat(0,nDim) = 0.5*density*oneOnC;
  pMat(0,nDim+1) = 0.5*density*oneOnC;

  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    pMat(iDim+1,nDim) = 0.5*density*(velocity(iDim)*oneOnC + normal(iDim));
    pMat(iDim+1,nDim+1) = 0.5*density*(velocity(iDim)*oneOnC - normal(iDim));
  }

  pMat(nDim+1,nDim) = 0.5*density*(vel2*oneOnC + projVel + soundSpeed/(gamma-1));
  pMat(nDim+1,nDim+1) = 0.5*density*(vel2*oneOnC - projVel + soundSpeed/(gamma-1));

  return pMat;
}

template<size_t nDim>
FORCEINLINE MatrixDbl<nDim+2> pMatrixInv(Double gamma, Double density, const VectorDbl<nDim>& velocity,
                                         Double projVel, Double soundSpeed, const VectorDbl<nDim>& normal) {
  MatrixDbl<nDim+2> pMatInv;

  auto c2 = pow(soundSpeed,2);
  auto vel2 = 0.5*squaredNorm(velocity);

  if (nDim == 2) {
    auto tmp = (gamma-1)/c2;
    pMatInv(0,0) =  1.0 - tmp*vel2;
    pMatInv(0,1) = tmp*velocity(0);
    pMatInv(0,2) = tmp*velocity(1);
    pMatInv(0,3) = -tmp;

    pMatInv(1,0) = (normal(0)*velocity(1)-normal(1)*velocity(0))/density;
    pMatInv(1,1) = normal(1)/density;
    pMatInv(1,2) = -normal(0)/density;
    pMatInv(1,3) = 0.0;
  }
  else {
    auto tmp = (gamma-1)/c2 * normal(0);
    pMatInv(0,0) = normal(0) - tmp*vel2 - (normal(2)*velocity(1)-normal(1)*velocity(2))/density;
    pMatInv(0,1) = tmp*velocity(0);
    pMatInv(0,2) = tmp*velocity(1) + normal(2)/density;
    pMatInv(0,3) = tmp*velocity(2) - normal(1)/density;
    pMatInv(0,4) = -tmp;

    tmp = (gamma-1)/c2 * normal(1);
    pMatInv(1,0) = normal(1) - tmp*vel2 + (normal(2)*velocity(0)-normal(0)*velocity(2))/density;
    pMatInv(1,1) = tmp*velocity(0) - normal(2)/density;
    pMatInv(1,2) = tmp*velocity(1);
    pMatInv(1,3) = tmp*velocity(2) + normal(0)/density;
    pMatInv(1,4) = -tmp;

    tmp = (gamma-1)/c2 * normal(2);
    pMatInv(2,0) = normal(2) - tmp*vel2 - (normal(1)*velocity(0)-normal(0)*velocity(1))/density;
    pMatInv(2,1) = tmp*velocity(0) + normal(1)/density;
    pMatInv(2,2) = tmp*velocity(1) - normal(0)/density+;
    pMatInv(2,3) = tmp*velocity(2);
    pMatInv(2,4) = -tmp;
  }

  /*--- Last two rows. ---*/

  auto gamma_minus_1_on_rho_times_c = (gamma-1) / (density*speedSound);

  for (size_t iVar = nDim; iVar < nDim+2; ++iVar) {
    Double sign = (iVar==nDim)? 1 : -1;
    pMatInv(iVar,0) = -sign*projVel/density + gamma_minus_1_on_rho_times_c * vel2;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      pMatInv(iVar,iDim+1) = sign*normal(iDim)/density - gamma_minus_1_on_rho_times_c * velocity(iDim);
    }
    pMatInv(iVar,nDim+1) = gamma_minus_1_on_rho_times_c;
  }

  return pMatInv;
}

template<size_t nDim>
FORCEINLINE VectorDbl<nDim+2> inviscidProjFlux(const CCompressiblePrimitives<nDim>& V
                                               const CCompressibleConservatives<nDim>& U
                                               const VectorDbl<nDim>& normal) {
  VectorDbl<nDim+2> flux = dot(U.momentum(), normal);
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    flux(iDim+1) *= V.velocity(iDim);
    flux(iDim+1) += V.pressure * normal(iDim);
  }
  flux(nDim+1) *= V.enthalpy;
  return flux;
}

template<size_t nDim, class RandomIterator>
FORCEINLINE MatrixDbl<nDim+2> inviscidProjJac(Double gamma, RandomIterator velocity,
                                              Double energy, const VectorDbl<nDim>& normal,
                                              Double scale) {
  MatrixDbl<nDim+2> jac;

  auto projVel = dot(velocity, normal);
  auto gamma_m_1 = gamma-1;
  auto phi = 0.5*gamma_m_1*squaredNorm<nDim>(velocity);
  auto a1 = gamma*energy - phi;

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
 * \class CRoeBase
 * \brief Base class for Roe schemes, derived classes implement
 * the dissipation term in a static method "finalizeResidual".
 */
template<class Derived, int NDIM>
class CRoeBase : public CNumericsSIMD {
protected:
  enum: size_t {nDim = NDIM};
  enum: size_t {nVar = CCompressibleConservatives<nDim>::nVar};

  const su2double kappa;
  const su2double gamma;
  const su2double entropyFix;
  const bool dynamicGrid;

  /*!
   * \brief Constructor.
   */
  CRoeScheme(const CConfig& config) :
    kappa(config.GetRoe_Kappa()),
    gamma(config.GetGamma()),
    entropyFix(config.GetEntropyFix_Coeff())
    dynamicGrid(config.GetDynamic_Grid()) {

  }

public:
  /*!
   * \brief Implementation of the base Roe flux.
   */
  void ComputeFlux(unsigned long iEdge,
                   const CConfig& config,
                   const CGeometry& geometry,
                   const CVariable& solution_,
                   UpdateType updateType,
                   CSysVector& vector,
                   CSysMatrix& matrix) override final {

    const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
    const auto solution = static_cast<CEulerVariable&>(solution_);

    const auto iPosize_t = geometry.edges->GetNodes(iEdge,0);
    const auto jPosize_t = geometry.edges->GetNodes(iEdge,1);

    /*--- Geometrix properties. ---*/

    const auto vector_ij = distanceVector(iPosize_t, jPosize_t, geometry.nodes->Coord);

    const auto normal = loadVariables<nDim>(iEdge, geometry.edges->Normal);
    const auto area = norm(normal);
    VectorDbl<nDim> unitNormal;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      unitNormal(iDim) = normal(iDim) / area;
    }

    /*--- Reconstructed primitives. ---*/

    auto V = getCompressiblePrimitives(iPosize_t, jPosize_t,
                                       config.GetKind_SlopeLimit_Flow(),
                                       vector_ij, solution);

    /*--- Compute conservative variables. ---*/

    CPair<CCompressibleConservatives<nDim> > U;
    U.i = getCompressibleConservatives(V.i);
    U.j = getCompressibleConservatives(V.j);

    /*--- Roe averaged variables. ---*/

    auto roeAvg = getRoeAveragedVariables(gamma, V, unitNormal);

    /*--- P tensor. ---*/

    auto pMat = pMatrix(gamma, roeAvg.density, roeAvg.velocity,
                        roeAvg.projVel, roeAvg.soundSpeed, unitNormal);

    /*--- Grid motion. ---*/

    Double projGridVel = 0.0, projVel = roeAvg.projVel;
    if (dynamicGrid) {
      const auto& gridVel = geometry.nodes->GridVel;
      projGridVel = 0.5*(dot(gridVel.innerIter<nDim>(iPosize_t), unitNormal)+
                         dot(gridVel.innerIter<nDim>(jPosize_t), unitNormal));
      projVel -= projGridVel;
    }

    /*--- Convective eigenvalues. ---*/

    VectorDbl<nVar> lambda = projVel;
    lambda(nDim) += roeAvg.soundSpeed;
    lambda(nDim+1) -= roeAvg.soundSpeed;

    /*--- Apply Mavriplis' entropy correction to eigenvalues. ---*/

    auto maxLambda = abs(projVel) + roeAvg.soundSpeed;

    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      lambda(iVar) = max(abs(lambda(iVar)), entropyFix*maxLambda);
    }

    /*--- Inviscid fluxes and Jacobians. ---*/

    auto flux_i = inviscidProjFlux(V.i, U.i, unitNormal);
    auto flux_j = inviscidProjFlux(V.j, U.j, unitNormal);

    VectorDbl<nVar> flux;
    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      flux(iVar) = kappa * (flux_i(iVar) + flux_j(iVar));
    }

    MatrixDbl<nVar> jac_i, jac_j;
    if (implicit) {
      jac_i = inviscidProjJac(gamma, V.i.velocity(), U.i.rhoEnergy/U.i.density, normal, kappa);
      jac_j = inviscidProjJac(gamma, V.j.velocity(), U.j.rhoEnergy/U.j.density, normal, kappa);
    }

    /*--- Finalize in derived class (static polymorphism). ---*/

    Derived::finalizeResidual(flux, jac_i, jac_j, implicit, gamma, kappa, area,
                              unitNormal, V, U, roeAvg, lambda, pMat);

    /*--- Correct for grid motion ---*/

    if (dynamicGrid) {
      for (size_t iVar = 0; iVar < nVar; iVar++) {
        Double dFdU = projGridVel * area * 0.5;
        flux(iVar) -= dFdU * (U.i.all(iVar) + U.j.all(iVar));

        if (implicit) {
          jac_i(iVar,iVar) -= dFdU;
          jac_j(iVar,iVar) -= dFdU;
        }
      }
    }

    /*--- Update the vector and system matrix. ---*/


    if (implicit) {


    }
  }
};

/*!
 * \class CRoeScheme
 * \brief Classical Roe scheme.
 */
template<int nDim>
class CRoeScheme : public CRoeBase<CRoeScheme,nDim> {
public:
  CRoeScheme(const CConfig& config) : CRoeBase(config) {}

  template<class... Ts>
  FORCEINLINE static void finalizeResidual(VectorDbl<nVar>& flux,
                                           MatrixDbl<nVar>& jac_i,
                                           MatrixDbl<nVar>& jac_j,
                                           bool implicit;
                                           Double gamma,
                                           Double kappa,
                                           Double area,
                                           const VectorDbl<nDim>& unitNormal,
                                           const CPair<CCompressiblePrimitives<nDim> >& V,
                                           const CPair<CCompressibleConservatives<nDim> >& U,
                                           const CRoeVariables& roeAvg,
                                           const VectorDbl<nVar>& lambda,
                                           const MatrixDbl<nVar>& pMat,
                                           Ts&...) {
    /*--- Inverse P tensor. ---*/

    auto pMatInv = pMatrixInv(gamma, roeAvg.density, roeAvg.velocity,
                              roeAvg.projVel, roeAvg.soundSpeed, unitNormal);

    /*--- Diference between conservative variables at jPosize_t and iPosize_t. ---*/

    VectorDbl<nVar> deltaU;
    for (iVar = 0; iVar < nVar; iVar++) {
      deltaU(iVar) = U.j.all(iVar) - U.i.all(iVar);
    }

    /*--- Dissipation terms. ---*/

    /// TODO: Low dissipation computation would go here.
    Double dissipation = 1.0;

    for (size_t iVar = 0; iVar < nVar; iVar++) {
      for (size_t jVar = 0; jVar < nVar; jVar++) {
        /*--- Compute |projModJacTensor| = P x |Lambda| x P^-1. ---*/

        Double projModJacTensor = 0.0;
        for (size_t kVar = 0; kVar < nVar; kVar++) {
          projModJacTensor += pMat(iVar,kVar) * lambda(kVar) * pMatInv(kVar,jVar);
        }

        auto dDdU = projModJacTensor * (1-kappa) * area * dissipation;

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
