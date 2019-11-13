/*!
 * \file CUserDefinedSolution.cpp
 * \brief Implementations of the member functions of CUserDefinedSolution.
 * \author T. Economon, E. van der Weide
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

CUserDefinedSolution::CUserDefinedSolution(void) : CVerificationSolution() { }

CUserDefinedSolution::CUserDefinedSolution(unsigned short val_nDim,
                                           unsigned short val_nVar,
                                           unsigned short val_iMesh,
                                           CConfig*       config)
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

CUserDefinedSolution::~CUserDefinedSolution(void) { }

void CUserDefinedSolution::GetBCState(const su2double *val_coords,
                                      const su2double val_t,
                                      su2double       *val_solution) {

  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
}

void CUserDefinedSolution::GetSolution(const su2double *val_coords,
                                       const su2double val_t,
                                       su2double       *val_solution) {

  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
}

void CUserDefinedSolution::GetMMSSourceTerm(const su2double *val_coords,
                                            const su2double val_t,
                                            su2double       *val_source) {

  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
}

bool CUserDefinedSolution::IsManufacturedSolution(void) {
  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
  return false;  /* True if manufactured. */
}
