/*!
 * \file CNumericsSIMD.hpp
 * \brief Vectorized (SIMD) numerics classes.
 * \author P. Gomes
 * \version 7.0.5 "Blackbird"
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

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/parallelization/vectorization.hpp"
#include "../../../Common/include/containers/C2DContainer.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

/*!
 * \enum UpdateType
 * \brief Ways to update vectors and system matrices.
 * COLORING is the typical i/j update, whereas for REDUCTION
 * the fluxes are stored and the matrix diagonal is not modified.
 */
enum class UpdateType {COLORING, REDUCTION};

#if !defined(CODI_REVERSE_TYPE) && !defined(CODI_FORWARD_TYPE)
/*--- Scalar types (SIMD). ---*/
using Double = simd::Array<su2double>;
using Int = simd::Array<unsigned long, Double::Size>;
constexpr size_t Int_Size = Int::Size;
#else
/*--- No SIMD with AD. ---*/
using Double = su2double;
using Int = unsigned long;
constexpr size_t Int_Size = 1;
#endif

/*--- Static vector types. ---*/
template<class Type, size_t Size>
using Vector = C2DContainer<unsigned long, Type, StorageType::ColumnMajor, simd::SIMD_SIZE, Size, 1>;

template<size_t Size> using VectorInt = Vector<Int, Size>;
template<size_t Size> using VectorDbl = Vector<Double, Size>;

/*--- Static matrix types. ---*/
template<class Type, size_t Rows, size_t Cols>
using Matrix = C2DContainer<unsigned long, Type, StorageType::RowMajor, simd::SIMD_SIZE, Rows, Cols>;

template<size_t Rows, size_t Cols = Rows> using MatrixInt = Matrix<Int, Rows, Cols>;
template<size_t Rows, size_t Cols = Rows> using MatrixDbl = Matrix<Double, Rows, Cols>;

#ifdef CODI_FORWARD_TYPE
using SparseMatrixType = CSysMatrix<su2double>;
#else
using SparseMatrixType = CSysMatrix<su2mixedfloat>;
#endif

/*!
 * \class CNumericsSIMD
 * \brief Base class to define the interface.
 * \note See CNumericsEmptyDecorator.
 */
class CNumericsSIMD {
public:
  /*!
   * \brief Interface for edge flux computation.
   * \param[in] iEdge - The edges for flux computation.
   * \param[in] config - Problem definitions.
   * \param[in] geometry - Problem geometry.
   * \param[in] solution - Solution variables.
   * \param[in] updateType - Type of update done on vector and matrix.
   * \param[in] updateMask - SIMD array of 1's and 0's, the latter prevent the update.
   * \param[in,out] vector - Target for the fluxes.
   * \param[in,out] matrix - Target for the flux Jacobians.
   * \note The update mask is used to handle "remainder" edges (nEdge mod simdSize).
   */
  virtual void ComputeFlux(Int iEdge,
                           const CConfig& config,
                           const CGeometry& geometry,
                           const CVariable& solution,
                           UpdateType updateType,
                           Double updateMask,
                           CSysVector<su2double>& vector,
                           SparseMatrixType& matrix) const = 0;

  /*! \brief Destructor of the class. */
  virtual ~CNumericsSIMD(void) = default;
};

/*!
 * \class CNumericsEmptyDecorator
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
class CNumericsEmptyDecorator : public CNumericsSIMD {
protected:
  enum: size_t {nDim = NDIM};

  template<class... Ts>
  CNumericsEmptyDecorator(Ts&...) {}

  /*!
   * \brief Empty method, real decorators should take as arguments whatever
   * the decorated class can pass them to avoid expensive data accesses.
   */
  template<class... Ts>
  void updateFlux(Ts&...) const {}
};

template<size_t nDim, class ForwardIterator, class T>
FORCEINLINE Double dot(ForwardIterator iterator, const T* ptr) {
  Double sum = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    sum += *(iterator++) * ptr[iDim];
  }
  return sum;
}

template<size_t nDim, class ForwardIterator>
FORCEINLINE Double dot(ForwardIterator iterator,
                       const VectorDbl<nDim>& vector) {
  return dot<nDim>(iterator, vector.data());
}

template<size_t nDim>
FORCEINLINE Double dot(const VectorDbl<nDim>& a, const VectorDbl<nDim>& b) {
  return dot<nDim>(a.data(), b.data());
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
FORCEINLINE Double squaredNorm(const VectorDbl<nDim>& vector) {
  return squaredNorm<nDim>(vector.data());
}

template<size_t nDim>
FORCEINLINE Double norm(const VectorDbl<nDim>& vector) { return sqrt(squaredNorm(vector)); }

template<size_t nDim, class Coord_t>
FORCEINLINE VectorDbl<nDim> distanceVector(Int iPoint, Int jPoint, const Coord_t& coords) {
  using T = VectorDbl<nDim>;
  auto coord_i = coords.template get<T>(iPoint);
  auto coord_j = coords.template get<T>(jPoint);
  T vector_ij;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    vector_ij(iDim) = 0.5 * (coord_j(iDim) - coord_i(iDim));
  }
  return vector_ij;
}

template<size_t nVar, class Field_t>
FORCEINLINE VectorDbl<nVar> gatherVariables(Int idx, const Field_t& vars) {
  return vars.template get<VectorDbl<nVar> >(idx);
}

template<size_t nRows, size_t nCols, class Field_t>
FORCEINLINE MatrixDbl<nRows,nCols> gatherVariables(Int idx, const Field_t& vars) {
  return vars.template get<MatrixDbl<nRows,nCols> >(idx);
}

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

template<class T>
struct CPair {
  T i, j;
};

template<size_t nDim, size_t nVar_= nDim+4>
struct CCompressiblePrimitives {
  enum : size_t {nVar = nVar_};
  VectorDbl<nVar> all;
  FORCEINLINE Double& temperature() { return all(0); }
  FORCEINLINE Double& pressure() { return all(nDim+1); }
  FORCEINLINE Double& density() { return all(nDim+2); }
  FORCEINLINE Double& enthalpy() { return all(nDim+3); }
  FORCEINLINE Double& velocity(size_t iDim) { return all(iDim+1); }
  FORCEINLINE const Double& temperature() const { return all(0); }
  FORCEINLINE const Double& pressure() const { return all(nDim+1); }
  FORCEINLINE const Double& density() const { return all(nDim+2); }
  FORCEINLINE const Double& enthalpy() const { return all(nDim+3); }
  FORCEINLINE const Double& velocity(size_t iDim) const { return all(iDim+1); }
  FORCEINLINE const Double* velocity() const { return &velocity(0); }

  FORCEINLINE Double& speedSound() { return all(nDim+4); }
  FORCEINLINE Double& laminarVisc() { return all(nDim+5); }
  FORCEINLINE Double& eddyVisc() { return all(nDim+6); }
  FORCEINLINE const Double& speedSound() const { return all(nDim+4); }
  FORCEINLINE const Double& laminarVisc() const { return all(nDim+5); }
  FORCEINLINE const Double& eddyVisc() const { return all(nDim+6); }
};

template<class PrimitiveType, size_t nDim, class Variable_t>
FORCEINLINE CPair<PrimitiveType> reconstructPrimitives(Int iPoint, Int jPoint, bool muscl,
                                                       ENUM_LIMITER limiterType,
                                                       const VectorDbl<nDim>& vector_ij,
                                                       const Variable_t& solution) {
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

template<size_t nDim>
struct CCompressibleConservatives {
  enum : size_t {nVar = nDim+2};
  VectorDbl<nVar> all;

  FORCEINLINE Double& density() { return all(0); }
  FORCEINLINE Double& rhoEnergy() { return all(nDim+1); }
  FORCEINLINE Double& momentum(size_t iDim) { return all(iDim+1); }
  FORCEINLINE const Double& density() const { return all(0); }
  FORCEINLINE const Double& rhoEnergy() const { return all(nDim+1); }
  FORCEINLINE const Double& momentum(size_t iDim) const { return all(iDim+1); }

  FORCEINLINE Double energy() const { return rhoEnergy() / density(); }
  FORCEINLINE const Double* momentum() const { return &momentum(0); }
};

template<size_t nDim>
FORCEINLINE CCompressibleConservatives<nDim> compressibleConservatives(const CCompressiblePrimitives<nDim>& V) {
  CCompressibleConservatives<nDim> U;
  U.density() = V.density();
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    U.momentum(iDim) = V.density() * V.velocity(iDim);
  }
  U.rhoEnergy() = V.density() * V.enthalpy() - V.pressure();
  return U;
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
FORCEINLINE CRoeVariables<nDim> roeAveragedVariables(Double gamma,
                                                     const CPair<CCompressiblePrimitives<nDim> >& V,
                                                     const VectorDbl<nDim>& normal) {
  CRoeVariables<nDim> roeAvg;
  Double R = sqrt(V.j.density() / V.i.density());
  Double D = 1 / (R+1);
  roeAvg.density = R * V.i.density();
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    roeAvg.velocity(iDim) = (R*V.j.velocity(iDim) + V.i.velocity(iDim)) * D;
  }
  roeAvg.enthalpy = (R*V.j.enthalpy() + V.i.enthalpy()) * D;
  roeAvg.speedSound = sqrt((gamma-1) * (roeAvg.enthalpy - 0.5*squaredNorm(roeAvg.velocity)));
  roeAvg.projVel = dot(roeAvg.velocity, normal);
  return roeAvg;
}

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
    pMat(2,1) = -1*density*normal(0);

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
    pMatInv(0,3) = -1*tmp;

    pMatInv(1,0) = (normal(0)*velocity(1)-normal(1)*velocity(0))*oneOnRho;
    pMatInv(1,1) = normal(1)*oneOnRho;
    pMatInv(1,2) = -1*normal(0)*oneOnRho;
    pMatInv(1,3) = 0.0;
  }
  else {
    Double tmp = (gamma-1)/c2 * normal(0);
    pMatInv(0,0) = normal(0) - tmp*vel2 - (normal(2)*velocity(1)-normal(1)*velocity(2))*oneOnRho;
    pMatInv(0,1) = tmp*velocity(0);
    pMatInv(0,2) = tmp*velocity(1) + normal(2)*oneOnRho;
    pMatInv(0,3) = tmp*velocity(2) - normal(1)*oneOnRho;
    pMatInv(0,4) = -1*tmp;

    tmp = (gamma-1)/c2 * normal(1);
    pMatInv(1,0) = normal(1) - tmp*vel2 + (normal(2)*velocity(0)-normal(0)*velocity(2))*oneOnRho;
    pMatInv(1,1) = tmp*velocity(0) - normal(2)*oneOnRho;
    pMatInv(1,2) = tmp*velocity(1);
    pMatInv(1,3) = tmp*velocity(2) + normal(0)*oneOnRho;
    pMatInv(1,4) = -1*tmp;

    tmp = (gamma-1)/c2 * normal(2);
    pMatInv(2,0) = normal(2) - tmp*vel2 - (normal(1)*velocity(0)-normal(0)*velocity(1))*oneOnRho;
    pMatInv(2,1) = tmp*velocity(0) + normal(1)*oneOnRho;
    pMatInv(2,2) = tmp*velocity(1) - normal(0)*oneOnRho;
    pMatInv(2,3) = tmp*velocity(2);
    pMatInv(2,4) = -1*tmp;
  }

  /*--- Last two rows. ---*/

  const Double gamma_minus_1_on_rho_times_c = (gamma-1) / (density*speedSound);

  for (size_t iVar = nDim; iVar < nDim+2; ++iVar) {
    Double sign = (iVar==nDim)? 1 : -1;
    pMatInv(iVar,0) = -1*sign*projVel*oneOnRho + gamma_minus_1_on_rho_times_c * vel2;
    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      pMatInv(iVar,iDim+1) = sign*normal(iDim)*oneOnRho - gamma_minus_1_on_rho_times_c * velocity(iDim);
    }
    pMatInv(iVar,nDim+1) = gamma_minus_1_on_rho_times_c;
  }

  return pMatInv;
}

template<size_t nDim>
FORCEINLINE VectorDbl<nDim+2> inviscidProjFlux(const CCompressiblePrimitives<nDim>& V,
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

template<size_t nVar>
FORCEINLINE void updateLinearSystem(Int iEdge,
                                    Int iPoint,
                                    Int jPoint,
                                    bool implicit,
                                    UpdateType updateType,
                                    Double updateMask,
                                    const VectorDbl<nVar>& flux,
                                    const MatrixDbl<nVar>& jac_i,
                                    const MatrixDbl<nVar>& jac_j,
                                    CSysVector<su2double>& vector,
                                    SparseMatrixType& matrix) {
  if (updateType == UpdateType::COLORING) {
    vector.UpdateBlocks(iPoint, jPoint, flux, updateMask);
    if(implicit) matrix.UpdateBlocks(iEdge, iPoint, jPoint, jac_i, jac_j, updateMask);
  }
  else {
    vector.SetBlock(iEdge, flux, updateMask);
    if(implicit) matrix.SetBlocks(iEdge, jac_i, jac_j, updateMask);
  }
}

/*!
 * \class CRoeBase
 * \brief Base class for Roe schemes, derived classes implement
 * the dissipation term in a static "finalizeFlux" method.
 * A base class implementing "updateFlux" is accepted as template parameter.
 * Similarly to derived, that method should update the flux and Jacobians, but
 * whereas "finalizeFlux" is passed data prepared by CRoeBase, "updateFlux"
 * takes the same input arguments as "ComputeFlux", i.e. it must fetch the data
 * it needs from CVariable. Derived is meant to implement small details,
 * Base is meant to do heavy lifting like computing viscous fluxes.
 */
template<class Derived, class Base>
class CRoeBase : public Base {
protected:
  enum: size_t {nDim = Base::nDim};
  enum: size_t {nVar = CCompressibleConservatives<nDim>::nVar};

  const su2double kappa;
  const su2double gamma;
  const su2double entropyFix;
  const bool finestGrid;
  const bool dynamicGrid;

  /*!
   * \brief Constructor, store some constants and forward args to base.
   */
  template<class... Ts>
  CRoeBase(const CConfig& config, unsigned iMesh, Ts&... args) : Base(config, iMesh, args...),
    kappa(config.GetRoe_Kappa()),
    gamma(config.GetGamma()),
    entropyFix(config.GetEntropyFix_Coeff()),
    finestGrid(iMesh == MESH_0),
    dynamicGrid(config.GetDynamic_Grid()) {
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

    const bool implicit = (config.GetKind_TimeIntScheme() == EULER_IMPLICIT);
    const bool muscl = finestGrid && config.GetMUSCL_Flow();
    const auto limiter = static_cast<ENUM_LIMITER>(config.GetKind_SlopeLimit_Flow());
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

    auto V = reconstructPrimitives<CCompressiblePrimitives<nDim> >(
               iPoint, jPoint, muscl, limiter, vector_ij, solution);

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
      lambda(iVar) = max(abs(lambda(iVar)), entropyFix*maxLambda);
    }

    /*--- Inviscid fluxes and Jacobians. ---*/

    auto flux_i = inviscidProjFlux(V.i, U.i, normal);
    auto flux_j = inviscidProjFlux(V.j, U.j, normal);

    VectorDbl<nVar> flux;
    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      flux(iVar) = kappa * (flux_i(iVar) + flux_j(iVar));
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

    Derived::finalizeFlux(flux, jac_i, jac_j, implicit, gamma, kappa,
                          area, unitNormal, V, U, roeAvg, lambda, pMat);

    /*--- Add the contributions from the base class (static decorator). ---*/

    Base::updateFlux(iEdge, iPoint, jPoint, config, geometry, solution_,
                     area, unitNormal, vector_ij, implicit, flux, jac_i, jac_j);

    /*--- Update the vector and system matrix. ---*/

    updateLinearSystem(iEdge, iPoint, jPoint, implicit, updateType,
                       updateMask, flux, jac_i, jac_j, vector, matrix);
  }
};

/*!
 * \class CRoeScheme
 * \brief Classical Roe scheme.
 */
template<class Decorator>
class CRoeScheme : public CRoeBase<CRoeScheme<Decorator>,Decorator> {
private:
  using Base = CRoeBase<CRoeScheme<Decorator>,Decorator>;
  enum: size_t {nDim = Base::nDim};
  enum: size_t {nVar = Base::nVar};
public:
  /*!
   * \brief Constructor, forward everything to base.
   */
  template<class... Ts>
  CRoeScheme(Ts&... args) : Base(args...) {}

  /*!
   * \brief Updates flux and Jacobians with standard Roe dissipation.
   * \note "Ts" is here just in case other schemes in the family need extra args.
   */
  template<class... Ts>
  FORCEINLINE static void finalizeFlux(VectorDbl<nVar>& flux,
                                       MatrixDbl<nVar>& jac_i,
                                       MatrixDbl<nVar>& jac_j,
                                       bool implicit,
                                       Double gamma,
                                       Double kappa,
                                       Double area,
                                       const VectorDbl<nDim>& unitNormal,
                                       const CPair<CCompressiblePrimitives<nDim> >& V,
                                       const CPair<CCompressibleConservatives<nDim> >& U,
                                       const CRoeVariables<nDim>& roeAvg,
                                       const VectorDbl<nVar>& lambda,
                                       const MatrixDbl<nVar>& pMat,
                                       Ts&...) {
    /*--- Inverse P tensor. ---*/

    auto pMatInv = pMatrixInv(gamma, roeAvg.density, roeAvg.velocity,
                              roeAvg.projVel, roeAvg.speedSound, unitNormal);

    /*--- Diference between conservative variables at jPoint and iPoint. ---*/

    VectorDbl<nVar> deltaU;
    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      deltaU(iVar) = U.j.all(iVar) - U.i.all(iVar);
    }

    /*--- Dissipation terms. ---*/

    /// TODO: Low dissipation computation would go here.
    Double dissipation = 1.0;

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

template<size_t nVar, size_t nDim, class PrimitiveType>
FORCEINLINE MatrixDbl<nDim> stressTensor(const PrimitiveType& V,
                                         const MatrixDbl<nVar,nDim> grad) {
  Double viscosity = V.laminarVisc() + V.eddyVisc();

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

template<size_t nVar, size_t nDim, class PrimitiveType>
FORCEINLINE VectorDbl<nVar> viscousFlux(const PrimitiveType& V,
                                        const MatrixDbl<nDim>& tau,
                                        const VectorDbl<nDim>& heatFlux,
                                        const VectorDbl<nDim>& normal) {
  VectorDbl<nVar> flux;
  flux(0) = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    /*--- Assuming symmetric stress tensor. ---*/
    flux(iDim+1) = dot(tau[iDim], normal);
  }
  flux(nDim+1) = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    auto viscWork = dot<nDim>(tau[iDim], V.velocity());
    flux(nDim+1) += normal(iDim) * (heatFlux(iDim) + viscWork);
  }
  return flux;
}

template<size_t nVar, size_t nDim, class PrimitiveType>
FORCEINLINE MatrixDbl<nDim,nVar> stressTensorJacobian(const PrimitiveType& V,
                                                      const VectorDbl<nDim> normal,
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
    jac(iDim,0) = -1 * dot<nDim>(&jac(iDim,1), V.velocity());
    /*--- Energy. ---*/
    jac(iDim,nDim+1) = 0.0;
  }
  return jac;
}

/*!
 * \class CCompressibleViscousFlux
 * \brief Decorator class to add viscous fluxes (compressible ideal gas).
 */
template<size_t NDIM>
class CCompressibleViscousFlux : public CNumericsSIMD {
protected:
  enum: size_t {nDim = NDIM};

  const su2double gamma;
  const su2double gasConst;
  const su2double prandtlLam;
  const su2double prandtlTurb;
  const su2double cp;
  const bool correct;

  template<class... Ts>
  CCompressibleViscousFlux(const CConfig& config, int iMesh, Ts&...) :
    gamma(config.GetGamma()),
    gasConst(config.GetGas_ConstantND()),
    prandtlLam(config.GetPrandtl_Lam()),
    prandtlTurb(config.GetPrandtl_Turb()),
    cp(gamma * gasConst / (gamma - 1)),
    correct(iMesh == MESH_0) {
  }

  /*!
   * \brief Add viscous contributions to flux and jacobians.
   */
  template<size_t nVar>
  FORCEINLINE void updateFlux(Int iEdge,
                              Int iPoint,
                              Int jPoint,
                              const CConfig& config,
                              const CGeometry& geometry,
                              const CVariable& solution_,
                              Double area,
                              const VectorDbl<nDim>& unitNormal,
                              VectorDbl<nDim> vector_ij,
                              bool implicit,
                              VectorDbl<nVar>& flux,
                              MatrixDbl<nVar>& jac_i,
                              MatrixDbl<nVar>& jac_j) const {

    constexpr size_t nPrimVar = nDim+7;
    constexpr size_t nPrimVarGrad = nDim+1;

    const auto& solution = static_cast<const CNSVariable&>(solution_);
    const auto& primitives = solution.GetPrimitive();
    const auto& gradient = solution.GetGradient_Primitive();

    /*--- (Re)compute some geometric properties. ---*/

    for (size_t iDim = 0; iDim < nDim; ++iDim) {
      /*--- Double distance from i to j as it was computed to face. ---*/
      vector_ij(iDim) *= 2;
    }
    auto dist2_ij = squaredNorm(vector_ij);
    /*--- Handle zero distance without "ifs" by making it large. ---*/
    Double mask = dist2_ij < EPS*EPS;
    dist2_ij += mask / (EPS*EPS);

    /*--- Get primitive variables and average them. ---*/

    CPair<CCompressiblePrimitives<nDim,nPrimVar> > V;
    V.i.all = gatherVariables<nPrimVar>(iPoint, primitives);
    V.j.all = gatherVariables<nPrimVar>(jPoint, primitives);

    CCompressiblePrimitives<nDim,nPrimVar> avgV;
    for (size_t iVar = 0; iVar < nPrimVar; ++iVar) {
      avgV.all(iVar) = 0.5 * (V.i.all(iVar) + V.j.all(iVar));
    }

    /*--- Compute the corrected mean gradient. ---*/

    auto avgGrad = averageGradient<nPrimVarGrad,nDim>(iPoint, jPoint, gradient);
    if(correct) correctGradient(V, vector_ij, dist2_ij, avgGrad);

    /// TODO: Wall shear stress (from wall functions).

    /// TODO: Uncertainty quantification (needs a way to access tke, maybe in ctor).

    /*--- Stress and heat flux tensors. ---*/

    auto tau = stressTensor(avgV, avgGrad);

    /// TODO: SA QCR.

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
};
