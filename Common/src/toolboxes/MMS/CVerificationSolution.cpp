/*!
 * \file CVerificationSolution.cpp
 * \brief Implementations of the member functions of CVerificationSolution.
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

#include "../../../include/toolboxes/MMS/CVerificationSolution.hpp"

CVerificationSolution::CVerificationSolution() {
  /*--- Initialize the pointers to NULL. ---*/
  Error_RMS = nullptr;
  Error_Max = nullptr;
  Error_Point_Max = nullptr;
  Error_Point_Max_Coord = nullptr;
}

CVerificationSolution::CVerificationSolution(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_iMesh,
                                             CConfig* config) {
  /*--- Store the kind of solver ---*/

  Kind_Solver = config->GetKind_Solver();

  /*--- Store the rank and size for the calculation. ---*/

  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();

  /*--- Store the dimension and number of variables. ---*/

  nDim = val_nDim;
  nVar = val_nVar;

  /*--- Allocate space for global error metrics for verification. ---*/
  Error_RMS = new su2double[nVar];
  Error_Max = new su2double[nVar];

  Error_Point_Max = new unsigned long[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) Error_Point_Max[iVar] = 0;

  Error_Point_Max_Coord = new su2double*[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Error_Point_Max_Coord[iVar] = new su2double[nDim];
    for (unsigned short iDim = 0; iDim < nDim; iDim++) Error_Point_Max_Coord[iVar][iDim] = 0.0;
  }
}

CVerificationSolution::~CVerificationSolution() {
  /*--- Release the memory of the pointers, if allocated. ---*/
  delete[] Error_RMS;
  delete[] Error_Max;

  delete[] Error_Point_Max;

  if (Error_Point_Max_Coord != nullptr) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete[] Error_Point_Max_Coord[iVar];
    }
    delete[] Error_Point_Max_Coord;
  }
}

void CVerificationSolution::GetSolution(const su2double* val_coords, const su2double val_t,
                                        su2double* val_solution) const {
  SU2_MPI::Error("Function must be overwritten by the derived class", CURRENT_FUNCTION);
}

void CVerificationSolution::GetInitialCondition(const su2double* val_coords, su2double* val_solution) const {
  /*--- Initial conditions call the GetSolution() method at t = 0. ---*/
  GetSolution(val_coords, 0.0, val_solution);
}

void CVerificationSolution::GetBCState(const su2double* val_coords, const su2double val_t,
                                       su2double* val_solution) const {
  SU2_MPI::Error("Function must be overwritten by the derived class", CURRENT_FUNCTION);
}

void CVerificationSolution::GetMMSSourceTerm(const su2double* val_coords, const su2double val_t,
                                             su2double* val_source) const {
  /* Default implementation of the source terms for the method of manufactured
     solutions. Simply set them to zero. */
  for (unsigned short iVar = 0; iVar < nVar; ++iVar) val_source[iVar] = 0.0;
}

bool CVerificationSolution::IsManufacturedSolution() const { return false; }

bool CVerificationSolution::ExactSolutionKnown() const { return true; }

void CVerificationSolution::GetLocalError(const su2double* val_coords, const su2double val_t,
                                          const su2double* val_solution, su2double* val_error) const {
  /*--- Get the value of the verification solution first.
        Use val_error to store this solution. ---*/

  GetSolution(val_coords, val_t, val_error);

  /*--- Compute the local error as the difference between the current
   numerical solution and the verification solution. ---*/

  for (unsigned short iVar = 0; iVar < nVar; ++iVar) val_error[iVar] = val_solution[iVar] - val_error[iVar];
}

void CVerificationSolution::SetVerificationError(unsigned long nDOFsGlobal, CConfig* config) {
  /* Disable the reduce for the error to avoid overhead if requested. */
  if (config->GetComm_Level() == COMM_FULL) {
    /*--- Get the number of ranks and the MPI communicator. ---*/
    int size = SU2_MPI::GetSize();
    SU2_MPI::Comm comm = SU2_MPI::GetComm();

    /*--- The local L2 norms must be added to obtain the global value. ---*/
    vector<su2double> rbufError(nVar, 0.0);
    SU2_MPI::Allreduce(Error_RMS, rbufError.data(), nVar, MPI_DOUBLE, MPI_SUM, comm);

    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      SetError_RMS(iVar, max(EPS * EPS, sqrt(rbufError[iVar] / nDOFsGlobal)));
    }

    /*--- The global maximum norms must be obtained. ---*/
    rbufError.resize(nVar * size, 0.0);
    SU2_MPI::Allgather(Error_Max, nVar, MPI_DOUBLE, rbufError.data(), nVar, MPI_DOUBLE, comm);

    vector<unsigned long> rbufPoint(nVar * size);
    SU2_MPI::Allgather(Error_Point_Max, nVar, MPI_UNSIGNED_LONG, rbufPoint.data(), nVar, MPI_UNSIGNED_LONG, comm);

    vector<su2double> sbufCoor(nDim * nVar, 0.0);
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      for (unsigned short iDim = 0; iDim < nDim; ++iDim)
        sbufCoor[iVar * nDim + iDim] = Error_Point_Max_Coord[iVar][iDim];
    }

    vector<su2double> rbufCoor(nDim * nVar * size, 0.0);
    SU2_MPI::Allgather(sbufCoor.data(), nVar * nDim, MPI_DOUBLE, rbufCoor.data(), nVar * nDim, MPI_DOUBLE, comm);

    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      for (int proc = 0; proc < size; ++proc)
        AddError_Max(iVar, rbufError[proc * nVar + iVar], rbufPoint[proc * nVar + iVar],
                     &rbufCoor[proc * nVar * nDim + iVar * nDim]);
    }
  }
}
