/*!
 * \file fem_work_estimate_metis.cpp
 * \brief This file contains the implementation of the member functions WorkEstimateMetis
          for the FEM standard elements.
 * \author E. van der Weide
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

#include "../../include/fem/fem_standard_element.hpp"

su2double CFEMStandardElement::WorkEstimateMetis(CConfig* config) {
  /* TEMPORARY IMPLEMENTATION. */
  return nIntegration + 0.1 * nDOFs;
}

su2double CFEMStandardInternalFace::WorkEstimateMetis(CConfig* config) {
  /* TEMPORARY IMPLEMENTATION. */
  return 2.0 * nIntegration + 0.05 * (nDOFsFaceSide0 + nDOFsFaceSide1);
}

su2double CFEMStandardBoundaryFace::WorkEstimateMetis(CConfig* config) {
  /* TEMPORARY IMPLEMENTATION. */
  return nIntegration + 0.05 * nDOFsFace;
}

su2double CFEMStandardBoundaryFace::WorkEstimateMetisWallFunctions(CConfig* config, const unsigned short nPointsWF) {
  /* TEMPORARY IMPLEMENTATION. */
  return 0.25 * nIntegration * nPointsWF;
}
