/*!
 * \file variables.hpp
 * \brief Collection of types to store physical variables.
 * \author P. Gomes
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

#include "../CNumericsSIMD.hpp"
#include "../util.hpp"

/*!
 * \brief Type to store compressible primitive variables and access them by name.
 */
template<size_t nDim_, size_t nVar_>
struct CCompressiblePrimitives {
  static constexpr size_t nDim = nDim_;
  static constexpr size_t nVar = nVar_;
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

  /*--- Un-reconstructed variables. ---*/
  FORCEINLINE Double& speedSound() { return all(nDim+4); }
  FORCEINLINE Double& laminarVisc() { return all(nDim+5); }
  FORCEINLINE Double& eddyVisc() { return all(nDim+6); }
  FORCEINLINE Double& thermalCond() { return all(nDim+7); }
  FORCEINLINE Double& cp() { return all(nDim+8); }
  FORCEINLINE const Double& speedSound() const { return all(nDim+4); }
  FORCEINLINE const Double& laminarVisc() const { return all(nDim+5); }
  FORCEINLINE const Double& eddyVisc() const { return all(nDim+6); }
  FORCEINLINE const Double& thermalCond() const { return all(nDim+7); }
  FORCEINLINE const Double& cp() const { return all(nDim+8); }
};

/*!
 * \brief Type to store compressible conservative (i.e. solution) variables.
 */
template<size_t nDim_>
struct CCompressibleConservatives {
  static constexpr size_t nDim = nDim_;
  static constexpr size_t nVar = nDim+2;
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

/*!
 * \brief Primitive to conservative conversion.
 */
template<size_t nDim, size_t N>
FORCEINLINE CCompressibleConservatives<nDim> compressibleConservatives(const CCompressiblePrimitives<nDim,N>& V) {
  CCompressibleConservatives<nDim> U;
  U.density() = V.density();
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    U.momentum(iDim) = V.density() * V.velocity(iDim);
  }
  U.rhoEnergy() = V.density() * V.enthalpy() - V.pressure();
  return U;
}

/*!
 * \brief Roe-averaged variables.
 */
template<size_t nDim>
struct CRoeVariables {
  Double density;
  VectorDbl<nDim> velocity;
  Double enthalpy;
  Double speedSound;
  Double projVel;
};

/*!
 * \brief Compute Roe-averaged variables from pair of primitive variables.
 */
template<size_t nDim, class PrimVarType>
FORCEINLINE CRoeVariables<nDim> roeAveragedVariables(Double gamma,
                                                     const CPair<PrimVarType>& V,
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
