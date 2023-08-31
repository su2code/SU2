/*!
 * \file CUserDefinedSolution.cpp
 * \brief Implementations of the member functions of CUserDefinedSolution.
 * \author T. Economon, E. van der Weide
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

#include "../../../include/toolboxes/MMS/CUserDefinedSolution.hpp"

CUserDefinedSolution::CUserDefinedSolution() : CVerificationSolution() {}

CUserDefinedSolution::CUserDefinedSolution(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_iMesh,
                                           CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  /*--- Write a message that the solution is initialized for a
   user-defined verification case. ---*/

  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for a user-defined verification case!!!" << endl;
    cout << endl << flush;
  }

  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
}

CUserDefinedSolution::~CUserDefinedSolution() = default;

void CUserDefinedSolution::GetBCState(const su2double* val_coords, const su2double val_t,
                                      su2double* val_solution) const {
  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
}

void CUserDefinedSolution::GetSolution(const su2double* val_coords, const su2double val_t,
                                       su2double* val_solution) const {
  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
}

void CUserDefinedSolution::GetMMSSourceTerm(const su2double* val_coords, const su2double val_t,
                                            su2double* val_source) const {
  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
}

bool CUserDefinedSolution::IsManufacturedSolution() const {
  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
  return false; /* True if manufactured. */
}
