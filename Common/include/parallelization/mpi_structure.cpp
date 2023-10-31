/*!
 * \file mpi_structure.cpp
 * \brief Main subroutines for the mpi structures.
 * \author T. Albring
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

#include "mpi_structure.hpp"
#include <cstring>  // memcpy

/* Initialise the MPI Communicator Rank and Size */
int CBaseMPIWrapper::Rank = 0;
int CBaseMPIWrapper::Size = 1;

/* Set the default MPI Communicator */
#ifdef HAVE_MPI
CBaseMPIWrapper::Comm CBaseMPIWrapper::currentComm = MPI_COMM_WORLD;
#else
CBaseMPIWrapper::Comm CBaseMPIWrapper::currentComm = 0;  // dummy value
#endif

#ifdef HAVE_MPI
int CBaseMPIWrapper::MinRankError;
bool CBaseMPIWrapper::winMinRankErrorInUse = false;
CBaseMPIWrapper::Win CBaseMPIWrapper::winMinRankError;

void CBaseMPIWrapper::Error(const std::string& ErrorMsg, const std::string& FunctionName) {
  /* Set MinRankError to Rank, as the error message is called on this rank. */
  MinRankError = Rank;
  int flag = 0;

#if MPI_VERSION >= 3
  /* Find out whether the error call is collective via MPI_Ibarrier. */
  Request barrierRequest;
  MPI_Ibarrier(currentComm, &barrierRequest);

  /* Try to complete the non-blocking barrier call for a second. */
  double startTime = SU2_MPI::Wtime();
  while (true) {
    MPI_Test(&barrierRequest, &flag, MPI_STATUS_IGNORE);
    if (flag) break;

    double currentTime = SU2_MPI::Wtime();
    if (currentTime > startTime + 1.0) break;
  }
#else
  /* MPI_Ibarrier function is not supported. Simply wait for one
     second to give other ranks the opportunity to reach this point. */
#ifdef _MSC_VER
  _sleep(1);
#else
  sleep(1);
#endif
#endif

  if (flag) {
    /* The barrier is completed and hence the error call is collective.
       Set MinRankError to 0. */
    MinRankError = 0;
  } else {
    /* The error call is not collective and the minimum rank must be
       determined by one sided communication. Loop over the lower numbered
       ranks to check if they participate in the error message. */
    for (int i = 0; i < Rank; ++i) {
      int MinRankErrorOther;
      MPI_Win_lock(MPI_LOCK_SHARED, i, 0, winMinRankError);
      MPI_Get(&MinRankErrorOther, 1, MPI_INT, i, 0, 1, MPI_INT, winMinRankError);
      MPI_Win_unlock(i, winMinRankError);

      if (MinRankErrorOther < MinRankError) {
        MinRankError = MinRankErrorOther;
        break;
      }
    }
  }

  /* Check if this rank must write the error message and do so. */
  if (Rank == MinRankError) {
    std::cout << std::endl << std::endl;
    std::cout << "Error in \"" << FunctionName << "\": " << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << ErrorMsg << std::endl;
    std::cout << "------------------------------ Error Exit -------------------------------" << std::endl;
    std::cout << std::endl << std::endl;
  }
  Abort(currentComm, EXIT_FAILURE);
}

void CBaseMPIWrapper::CopyData(const void* sendbuf, void* recvbuf, int size, Datatype datatype, int recvshift,
                               int sendshift) {
  if (datatype == MPI_DOUBLE) {
    for (int i = 0; i < size; i++) {
      static_cast<su2double*>(recvbuf)[i + recvshift] = static_cast<const su2double*>(sendbuf)[i + sendshift];
    }
  } else {
    int scalarsize;
    MPI_Type_size(datatype, &scalarsize);
    const char* src = static_cast<const char*>(sendbuf) + sendshift * static_cast<size_t>(scalarsize);
    char* dest = static_cast<char*>(recvbuf) + recvshift * static_cast<size_t>(scalarsize);
    std::memcpy(static_cast<void*>(dest), static_cast<const void*>(src), size * static_cast<size_t>(scalarsize));
  }
}
#else  // HAVE_MPI

void CBaseMPIWrapper::Error(std::string ErrorMsg, std::string FunctionName) {
  if (Rank == 0) {
    std::cout << std::endl << std::endl;
    std::cout << "Error in \"" << FunctionName << "\": " << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << ErrorMsg << std::endl;
    std::cout << "------------------------------ Error Exit -------------------------------" << std::endl;
    std::cout << std::endl << std::endl;
  }
  Abort(currentComm, 0);
}

void CBaseMPIWrapper::CopyData(const void* sendbuf, void* recvbuf, int size, Datatype datatype, int recvshift,
                               int sendshift) {
  switch (datatype) {
    case MPI_DOUBLE:
      for (int i = 0; i < size; i++) {
        static_cast<su2double*>(recvbuf)[i + recvshift] = static_cast<const su2double*>(sendbuf)[i + sendshift];
      }
      break;
    case MPI_UNSIGNED_LONG:
      for (int i = 0; i < size; i++) {
        static_cast<unsigned long*>(recvbuf)[i + recvshift] = static_cast<const unsigned long*>(sendbuf)[i + sendshift];
      }
      break;
    case MPI_LONG:
      for (int i = 0; i < size; i++) {
        static_cast<long*>(recvbuf)[i + recvshift] = static_cast<const long*>(sendbuf)[i + sendshift];
      }
      break;
    case MPI_UNSIGNED_SHORT:
      for (int i = 0; i < size; i++) {
        static_cast<unsigned short*>(recvbuf)[i + recvshift] =
            static_cast<const unsigned short*>(sendbuf)[i + sendshift];
      }
      break;
    case MPI_CHAR:
      for (int i = 0; i < size; i++) {
        static_cast<char*>(recvbuf)[i + recvshift] = static_cast<const char*>(sendbuf)[i + sendshift];
      }
      break;
    case MPI_SHORT:
      for (int i = 0; i < size; i++) {
        static_cast<short*>(recvbuf)[i + recvshift] = static_cast<const short*>(sendbuf)[i + sendshift];
      }
      break;
    case MPI_INT:
      for (int i = 0; i < size; i++) {
        static_cast<int*>(recvbuf)[i + recvshift] = static_cast<const int*>(sendbuf)[i + sendshift];
      }
      break;
    default:
      Error("Unknown type", CURRENT_FUNCTION);
      break;
  };
}
#endif

#ifdef HAVE_MPI
#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE
MediTypes* mediTypes;
#include <medi/medi.cpp>
#endif  // defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE

#endif  // HAVE_MPI
