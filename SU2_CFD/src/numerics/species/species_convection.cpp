/*!
 * \file species_convection.cpp
 * \brief Implementation of numerics classes to compute convective
 *        fluxes in turbulence problems.
 * \author T. Kattmann
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/species/species_convection.hpp"

CUpwSca_Species::CUpwSca_Species(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
    : CUpwScalar(val_nDim, val_nVar, config) {}

void CUpwSca_Species::ExtraADPreaccIn() {
  AD::SetPreaccIn(V_i, nDim + 3);  /// Note TK::Check the +3 addition!
  AD::SetPreaccIn(V_j, nDim + 3);
}

void CUpwSca_Species::FinishResidualCalc(const CConfig* config) {
  for (auto iVar = 0u; iVar < nVar; iVar++) {
    Flux[iVar] = a0 * Density_i * ScalarVar_i[iVar] + a1 * Density_j * ScalarVar_j[iVar];

    if (implicit) {
      for (auto jVar = 0u; jVar < nVar; jVar++) {
        if (iVar == jVar) {
          // note that Jacobians are taken wrt rho*Y not Y alone.
          Jacobian_i[iVar][jVar] = a0;
          Jacobian_j[iVar][jVar] = a1;
        } else {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }  // jVar
    }
  }  // iVar
}
