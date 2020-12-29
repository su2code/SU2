/*!
 * \file CInternalFaceFEM_DG.cpp
 * \brief Implementations of the member functions of CInternalFaceFEM_DG.
 * \author E. van der Weide
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

#include "../../../include/geometry/fem_grid/CInternalFaceFEM_DG.hpp"
#include "../../../include/geometry/fem_grid/CVolumeElementFEM_DG.hpp"

/*---------------------------------------------------------------------*/
/*---        Public member functions of CInternalFaceFEM_DG.        ---*/
/*---------------------------------------------------------------------*/

void CInternalFaceFEM_DG::InitGridVelocities(const unsigned short nDim) {

  /*--- Determine the padded number of integration points. ---*/
  const unsigned short nIntPad = standardElemGrid->GetNIntegrationPad();

  /*--- Allocate the memory for the grid velocities and initialize
        them to zero. ---*/
  gridVelocities.resize(nIntPad, nDim);
  gridVelocities.setConstant(0.0);
}

void CInternalFaceFEM_DG::MetricTermsIntegrationPoints(const bool                         viscousTerms,
                                                       const unsigned short               nDim,
                                                       const vector<CVolumeElementFEM_DG> &volElem) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}
