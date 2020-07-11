/*!
 * \file geometry_structure_fem_part.cpp
 * \brief Main subroutines for distributin the grid for the Fluid FEM solver.
 * \author F. Palacios, T. Economon
 * \version 7.0.6 "Blackbird"
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

#include "../../include/geometry/CPhysicalGeometry.hpp"
#include "../../include/fem/fem_standard_element.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridFEM.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridBoundFEM.hpp"

#include "../../include/adt_structure.hpp"
#include "../../include/blas_structure.hpp"
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>

void CPhysicalGeometry::SetColorFEMGrid_Parallel(CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CPhysicalGeometry::DeterminePeriodicFacesFEMGrid(CConfig                *config,
                                                      vector<CFaceOfElement> &localFaces) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CPhysicalGeometry::DetermineFEMConstantJacobiansAndLenScale(CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CPhysicalGeometry::DetermineDonorElementsWallFunctions(CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CPhysicalGeometry::DetermineTimeLevelElements(
                          CConfig                              *config,
                          const vector<CFaceOfElement>         &localFaces,
                          map<unsigned long, CUnsignedShort2T> &mapExternalElemIDToTimeLevel) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CPhysicalGeometry::ComputeFEMGraphWeights(
              CConfig                                    *config,
              const vector<CFaceOfElement>               &localFaces,
              const vector<vector<unsigned long> >       &adjacency,
              const map<unsigned long, CUnsignedShort2T> &mapExternalElemIDToTimeLevel,
                    vector<su2double>                    &vwgt,
                    vector<vector<su2double> >           &adjwgt) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}
