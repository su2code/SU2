/*!
 * \file mpi_structure.hpp
 * \brief Headers of the mpi interface for generalized datatypes.
 *        The subroutines and functions are in the <i>mpi_structure.cpp</i> file.
 * \author T. Albring
 * \version 3.2.9 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (francisco.palacios@boeing.com).
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

#pragma once

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "./datatype_structure.hpp"

#ifdef HAVE_MPI

/*!
 * \class CMPIWrapper
 * \brief Class for defining the MPI wrapper routines; this class features as a base class for
 * MPI interfaces for non-primitive dataypes e.g. used by AD, complex etc.
 * \author T. Albring
 * \version 3.2.9 "eagle"
 */

class CMPIWrapper {
public:
  static void Init(int *argc, char***argv);

  static void Isend(void *buf, int count, MPI_Datatype datatype, int dest,
                    int tag, MPI_Comm comm, MPI_Request* request);

  static void Irecv(void *buf, int count, MPI_Datatype datatype, int source,
                    int tag, MPI_Comm comm, MPI_Request* request);

  static void Wait(MPI_Request *request, MPI_Status *status);

  static void Waitall(int nrequests, MPI_Request *request, MPI_Status *status);

  static void Send(void *buf, int count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm);

  static void Recv(void *buf, int count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm, MPI_Status *status);

  static void Bcast(void *buf, int count, MPI_Datatype datatype, int root,
                    MPI_Comm comm);

  static void Bsend(void *buf, int count, MPI_Datatype datatype, int dest,
                    int tag, MPI_Comm comm);

  static void Reduce(void *sendbuf, void *recvbuf, int count,
                     MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

  static void Allreduce(void *sendbuf, void *recvbuf, int count,
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

  static void Gather(void *sendbuf, int sendcnt,MPI_Datatype sendtype,
                     void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);

  static void Allgather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                        void *recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm);

  static void Sendrecv(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                       int dest, int sendtag, void *recvbuf, int recvcnt,
                       MPI_Datatype recvtype,int source, int recvtag,
                       MPI_Comm comm, MPI_Status *status);

protected:
  static char* buff;

};
#endif
#include "mpi_structure.inl"
