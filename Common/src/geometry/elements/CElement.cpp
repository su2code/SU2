/*!
 * \file CElement.cpp
 * \brief Definition of the Finite Element structure (elements)
 * \author R. Sanchez
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

#include "../../../include/geometry/elements/CElement.hpp"

CElement::CElement(unsigned short ngauss, unsigned short nnodes, unsigned short ndim) {
  nGaussPoints = ngauss;
  nNodes = nnodes;
  nDim = ndim;

  /*--- Allocate structures. ---*/

  /*--- Always allocate MAXNDIM to have enough space in all cases.
   * e.g. elements embedded on a curved surfaces. ---*/
  CurrentCoord.resize(nNodes, MAXNDIM) = su2double(0.0);
  RefCoord.resize(nNodes, MAXNDIM) = su2double(0.0);

  GaussPoint.reserve(nGaussPoints);
  for (unsigned short iGauss = 0; iGauss < nGaussPoints; ++iGauss) GaussPoint.emplace_back(iGauss, nDim, nNodes);

  GaussWeight.resize(nGaussPoints) = su2double(0.0);
  NodalExtrap.resize(nNodes, nGaussPoints) = su2double(0.0);

  NodalStress.resize(nNodes, 6) = su2double(0.0);

  Mab.resize(nNodes, nNodes);
  Ks_ab.resize(nNodes, nNodes);
  Kab.resize(nNodes);
  for (auto& kab_i : Kab) kab_i.resize(nNodes, nDim * nDim);

  Kt_a.resize(nNodes, nDim);
  FDL_a.resize(nNodes, nDim);

  HiHj.resize(nNodes, nNodes);
  DHiDHj.resize(nNodes);
  for (auto& DHiDHj_a : DHiDHj) {
    DHiDHj_a.resize(nNodes);
    for (auto& DHiDHj_ab : DHiDHj_a) DHiDHj_ab.resize(nDim, nDim);
  }

  ClearElement();
}

void CElement::ClearElement() {
  Mab.setConstant(0.0);
  Kt_a.setConstant(0.0);
  FDL_a.setConstant(0.0);
  Ks_ab.setConstant(0.0);

  HiHj.setConstant(0.0);
  for (auto& DHiDHj_a : DHiDHj) {
    for (auto& DHiDHj_ab : DHiDHj_a) DHiDHj_ab.setConstant(0.0);
  }
  for (auto& kab_i : Kab) kab_i.setConstant(0.0);
}
