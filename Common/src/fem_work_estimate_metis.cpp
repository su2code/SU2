/*!
 * \file fem_work_estimate_metis.cpp
 * \brief This file contains the implementation of the member functions WorkEstimateMetis
          for the FEM standard elements.
 * \author E. van der Weide
 * \version 6.2.0 "Falcon"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/fem_standard_element.hpp"

su2double CFEMStandardElement::WorkEstimateMetis(CConfig *config) {

  /* TEMPORARY IMPLEMENTATION. */
  return nIntegration + 0.1*nDOFs;
}

su2double CFEMStandardInternalFace::WorkEstimateMetis(CConfig *config) {

  /* TEMPORARY IMPLEMENTATION. */
  return 2.0*nIntegration + 0.05*(nDOFsFaceSide0 + nDOFsFaceSide1);
}

su2double CFEMStandardBoundaryFace::WorkEstimateMetis(CConfig *config) {

  /* TEMPORARY IMPLEMENTATION. */
  return nIntegration + 0.05*nDOFsFace;
}

su2double CFEMStandardBoundaryFace::WorkEstimateMetisWallFunctions(
                                           CConfig              *config,
                                           const unsigned short nPointsWF) {

  /* TEMPORARY IMPLEMENTATION. */
  return 0.25*nIntegration*nPointsWF;
}
