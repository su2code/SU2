/*!
 * \file CNumericsSIMD.hpp
 * \brief Vectorized (SIMD) numerics classes.
 * \note This should be the only cpp for this family of classes
 * (which are all templates). All compilation takes place here.
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

#include "CNumericsSIMD.hpp"
#include "flow/convection/roe.hpp"
#include "flow/convection/centered.hpp"
#include "flow/diffusion/viscous_fluxes.hpp"

namespace {

/*!
 * \brief Upwind factory implementation for ideal gas.
 */
template<class ViscousDecorator>
CNumericsSIMD* createUpwindIdealNumerics(const CConfig& config, int iMesh, const CVariable* turbVars) {
  CNumericsSIMD* obj = nullptr;
  switch (config.GetKind_Upwind_Flow()) {
    case UPWIND::ROE:
      obj = new CRoeScheme<ViscousDecorator>(config, iMesh, turbVars);
      break;
    default:
      break;
  }
  return obj;
}

/*!
 * \brief Upwind factory implementation for real gas.
 */
template<class ViscousDecorator>
CNumericsSIMD* createUpwindGeneralNumerics(const CConfig& config, int iMesh, const CVariable* turbVars) {
  return nullptr;
}

/*!
 * \brief Centered factory implementation.
 */
template<class ViscousDecorator>
CNumericsSIMD* createCenteredNumerics(const CConfig& config, int iMesh, const CVariable* turbVars) {
  CNumericsSIMD* obj = nullptr;
  switch ((iMesh==MESH_0)? config.GetKind_Centered_Flow() : CENTERED::LAX) {
    case CENTERED::NONE:
      break;
    case CENTERED::LAX:
      obj = new CLaxScheme<ViscousDecorator>(config, iMesh, turbVars);
      break;
    case CENTERED::JST:
      obj = new CJSTScheme<ViscousDecorator>(config, iMesh, turbVars);
      break;
    case CENTERED::JST_KE:
      obj = new CJSTkeScheme<ViscousDecorator>(config, iMesh, turbVars);
      break;
    case CENTERED::JST_MAT:
      obj = new CJSTmatScheme<ViscousDecorator>(config, iMesh, turbVars);
      break;
  }
  return obj;
}

/*!
 * \brief Generic factory implementation.
 */
template<int nDim>
CNumericsSIMD* createNumerics(const CConfig& config, int iMesh, const CVariable* turbVars) {
  CNumericsSIMD* obj = nullptr;
  const bool ideal_gas = (config.GetKind_FluidModel() == STANDARD_AIR) ||
                         (config.GetKind_FluidModel() == IDEAL_GAS);

  switch (config.GetKind_ConvNumScheme_Flow()) {
    case SPACE_UPWIND:
      if (config.GetViscous()) {
        if (ideal_gas)
          obj = createUpwindIdealNumerics<CCompressibleViscousFlux<nDim> >(config, iMesh, turbVars);
        else
          obj = createUpwindGeneralNumerics<CGeneralCompressibleViscousFlux<nDim> >(config, iMesh, turbVars);
      }
      else {
        if (ideal_gas)
          obj = createUpwindIdealNumerics<CNoViscousFlux<nDim> >(config, iMesh, turbVars);
        else
          obj = createUpwindGeneralNumerics<CNoViscousFlux<nDim> >(config, iMesh, turbVars);
      }
      break;

    case SPACE_CENTERED:
      if (config.GetViscous()) {
        if (ideal_gas)
          obj = createCenteredNumerics<CCompressibleViscousFlux<nDim> >(config, iMesh, turbVars);
        else
          obj = createCenteredNumerics<CGeneralCompressibleViscousFlux<nDim> >(config, iMesh, turbVars);
      }
      else {
        obj = createCenteredNumerics<CNoViscousFlux<nDim> >(config, iMesh, turbVars);
      }
      break;
    default:
    break;
  }

  return obj;
}

} // namespace

/*!
 * \brief This function instantiates both 2D and 3D versions of the implementation in
 * createNumerics, which in turn instantiates the class templates of the different
 * numerical methods.
 */
CNumericsSIMD* CNumericsSIMD::CreateNumerics(const CConfig& config, int nDim, int iMesh, const CVariable* turbVars) {
#ifndef CODI_REVERSE_TYPE
  if ((Double::Size < 4) && (SU2_MPI::GetRank() == MASTER_NODE)) {
    cout << "WARNING: SU2 was not compiled for an AVX-capable architecture. Performance could be better,\n"
            "         see https://su2code.github.io/docs_v7/Build-SU2-Linux-MacOS/#compiler-optimizations" << endl;
  }
#endif
  if (nDim == 2) return createNumerics<2>(config, iMesh, turbVars);
  if (nDim == 3) return createNumerics<3>(config, iMesh, turbVars);

  return nullptr;
}
