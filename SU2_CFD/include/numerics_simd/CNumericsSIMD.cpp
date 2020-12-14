/*!
 * \file CNumericsSIMD.hpp
 * \brief Vectorized (SIMD) numerics classes.
 * \note This should be the only cpp for this family of classes
 * (which are all templates). All compilation takes place here.
 * \author P. Gomes
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

#include "CNumericsSIMD.hpp"
#include "flow/convection/roe.hpp"
#include "flow/convection/centered.hpp"
#include "flow/diffusion/viscous_fluxes.hpp"

/*!
 * \brief Generic factory implementation.
 */
template<class ViscousDecorator>
CNumericsSIMD* createNumerics(const CConfig& config, int iMesh) {
  CNumericsSIMD* obj = nullptr;
  switch (config.GetKind_ConvNumScheme_Flow()) {
    case SPACE_UPWIND:
      switch (config.GetKind_Upwind_Flow()) {
        case ROE:
          obj = new CRoeScheme<ViscousDecorator>(config, iMesh);
          break;
      }
      break;

    case SPACE_CENTERED:
      switch ((iMesh==MESH_0)? config.GetKind_Centered_Flow() : LAX) {
        case NO_CENTERED:
          break;
        case LAX:
          obj = new CLaxScheme<ViscousDecorator>(config, iMesh);
          break;
        case JST:
          obj = new CJSTScheme<ViscousDecorator>(config, iMesh);
          break;
        case JST_KE:
          obj = new CJSTkeScheme<ViscousDecorator>(config, iMesh);
          break;
        case JST_MAT:
          obj = new CJSTmatScheme<ViscousDecorator>(config, iMesh);
          break;
      }
      break;
  }
  return obj;
}

/*!
 * \brief This function instantiates both 2D and 3D versions of the implementation in
 * createNumerics, which in turn instantiates the class templates of the different
 * numerical methods.
 */
CNumericsSIMD* CNumericsSIMD::CreateNumerics(const CConfig& config, int nDim, int iMesh) {
  if ((Double::Size < 4) && (SU2_MPI::GetRank() == MASTER_NODE)) {
    cout << "WARNING: SU2 was not compiled for an AVX-capable architecture." << endl;
  }
  CNumericsSIMD* obj = nullptr;
  if (config.GetViscous()) {
    if (nDim == 2) obj = createNumerics<CCompressibleViscousFlux<2> >(config, iMesh);
    if (nDim == 3) obj = createNumerics<CCompressibleViscousFlux<3> >(config, iMesh);
  } else {
    if (nDim == 2) obj = createNumerics<CNoViscousFlux<2> >(config, iMesh);
    if (nDim == 3) obj = createNumerics<CNoViscousFlux<3> >(config, iMesh);
  }
  return obj;
}
