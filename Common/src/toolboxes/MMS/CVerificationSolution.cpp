/*!
 * \file CVerificationSolution.cpp
 * \brief Implementations of the member functions of CVerificationSolution.
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

#include "../../../include/toolboxes/MMS/CVerificationSolution.hpp"

CVerificationSolution::CVerificationSolution(void) { }

CVerificationSolution::CVerificationSolution(unsigned short val_nDim,
                                             unsigned short val_nVar,
                                             unsigned short val_iMesh,
                                             CConfig*       config) {
  
  /*--- Store the rank and size for the calculation. ---*/
  
  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();
  
  /*--- Store the dimension and number of variables. ---*/
  
  nDim = val_nDim;
  nVar = val_nVar;
  
}

CVerificationSolution::~CVerificationSolution(void) { }

void CVerificationSolution::GetSolution(const su2double *val_coords,
                                        const su2double val_t,
                                        su2double       *val_solution) {

  SU2_MPI::Error("Function must be overwritten by the derived class", CURRENT_FUNCTION);
}

void CVerificationSolution::GetInitialCondition(const su2double *val_coords,
                                                su2double       *val_solution) {
  
  /*--- Initial conditions call the GetSolution() method at t = 0. ---*/
  GetSolution(val_coords, 0.0, val_solution);
}

void CVerificationSolution::GetBCState(const su2double *val_coords,
                                       const su2double val_t,
                                       su2double       *val_solution) {

  SU2_MPI::Error("Function must be overwritten by the derived class", CURRENT_FUNCTION);
}

void CVerificationSolution::GetMMSSourceTerm(const su2double *val_coords,
                                             const su2double val_t,
                                             su2double       *val_source) {

  /* Default implementation of the source terms for the method of manufactured
     solutions. Simply set them to zero. */
  for(unsigned short iVar=0; iVar<nVar; ++iVar)
    val_source[iVar] = 0.0;
}

bool CVerificationSolution::IsManufacturedSolution(void) {return false;}

bool CVerificationSolution::ExactSolutionKnown(void) {return true;}

void CVerificationSolution::GetLocalError(const su2double *val_coords,
                                          const su2double val_t,
                                          const su2double *val_solution,
                                          su2double       *val_error) {
  
  /*--- Get the value of the verification solution first.
        Use val_error to store this solution. ---*/
  
  GetSolution(val_coords, val_t, val_error);
  
  /*--- Compute the local error as the difference between the current
   numerical solution and the verification solution. ---*/
  
  for (unsigned short iVar=0; iVar<nVar; ++iVar)
    val_error[iVar] = val_solution[iVar] - val_error[iVar];
  
}
