/*!
 * \file mpi_structure.hpp
 * \brief Headers of the mpi interface for generalized datatypes.
 *        The subroutines and functions are in the <i>mpi_structure.cpp</i> file.
 * \author T. Albring
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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
#include <map>
#endif

#include "./datatype_structure.hpp"
#include <stdlib.h>

#ifdef HAVE_MPI

/*--- Depending on the datatype used, the correct MPI wrapper class is defined.
 * For the default (double type) case this results in using the normal MPI routines. ---*/
#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE

#include <medi/medi.hpp>
using namespace medi;

class CMediMPIWrapper;
typedef CMediMPIWrapper SU2_MPI;

#if defined CODI_REVERSE_TYPE
#include <medi/codiMediPackTypes.hpp>
#if CODI_PRIMAL_INDEX_TAPE
typedef CoDiPackToolPrimalRestore<su2double> MediTool;
#else
typedef CoDiPackTool<su2double> MediTool;
#endif // defined CODI_REVERSE_TYPE
#elif defined CODI_FORWARD_TYPE
#include <medi/codiForwardMediPackTypes.hpp>
typedef CoDiPackForwardTool<su2double> MediTool;
#endif // defined CODI_FORWARD_TYPE

typedef medi::MpiTypeDefault<MediTool> AMPI_ADOUBLE_TYPE;
extern AMPI_ADOUBLE_TYPE* AMPI_ADOUBLE;

#else
class CMPIWrapper;
typedef CMPIWrapper SU2_MPI;
#endif // defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE

/*!
 * \class CMPIWrapper
 * \brief Class for defining the MPI wrapper routines; this class features as a base class for
 * MPI interfaces for non-primitive dataypes e.g. used by AD, complex etc.
 * \author T. Albring
 * \version 5.0.0 "Raven"
 */

class CMPIWrapper {

public:

  typedef MPI_Comm Comm;
  
  typedef MPI_Status Status;

#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE
  typedef AMPI_Request Request;
#else
  typedef MPI_Request Request;
#endif

  static void Init(int *argc, char***argv);
  
  static void Buffer_attach(void *buffer, int size);
  
  static void Buffer_detach(void *buffer, int *size);
  
  static void Finalize();
  
  static void Comm_rank(Comm comm, int* rank);
  
  static void Comm_size(Comm comm, int* size);
  
  static void Barrier(Comm comm);
  
  static void Abort(Comm comm, int error);
  
  static inline MPI_Request* convertRequest(Request* request);

  static void Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);

  static void Isend(void *buf, int count, MPI_Datatype datatype, int dest,
                    int tag, MPI_Comm comm, Request* request);

  static void Irecv(void *buf, int count, MPI_Datatype datatype, int source,
                    int tag, MPI_Comm comm, Request* request);

  static void Wait(Request *request, MPI_Status *status);

  static void Waitall(int nrequests, Request *request, MPI_Status *status);

  static void Waitany(int nrequests, Request *request,
                      int *index, MPI_Status *status);

  static void Testall(int count, Request* array_of_requests, int *flag, MPI_Status* array_of_statuses);

  static void Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);

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

  static void Scatter(void *sendbuf, int sendcnt,MPI_Datatype sendtype,
                      void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);

  static void Allgather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                        void *recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm);

  static void Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                         void *recvbuf, int *recvcounts, int *displs,
                         MPI_Datatype recvtype, MPI_Comm comm);

  static void Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                       void *recvbuf, int recvcount, MPI_Datatype recvtype,
                       MPI_Comm comm);

  static void Sendrecv(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                       int dest, int sendtag, void *recvbuf, int recvcnt,
                       MPI_Datatype recvtype,int source, int recvtag,
                       MPI_Comm comm, MPI_Status *status);
};

typedef MPI_Comm SU2_Comm;

#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE

/*!
 * \class CMediMPIWrapper
 * \brief MPI wrapper functions for MediPack tool.
 * \author T. Albring
 * \version 5.0.0 "Raven"
 */

class CMediMPIWrapper: public CMPIWrapper {
public:
  static void Init(int *argc, char***argv);

  static AMPI_Comm convertComm(MPI_Comm comm);

  static AMPI_Status* convertStatus(MPI_Status* status);

  static AMPI_Op convertOp(MPI_Op op);

  static void Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);

  static void Isend(void *buf, int count, MPI_Datatype datatype, int dest,
                    int tag, MPI_Comm comm, Request* request);

  static void Irecv(void *buf, int count, MPI_Datatype datatype, int source,
                    int tag, MPI_Comm comm, Request* request);

  static void Wait(SU2_MPI::Request *request, MPI_Status *status);

  static void Waitall(int nrequests, Request *request, MPI_Status *status);

  static void Waitany(int nrequests, SU2_MPI::Request *request,
                      int *index, MPI_Status *status);

  static void Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);

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

  static void Scatter(void *sendbuf, int sendcnt,MPI_Datatype sendtype,
                      void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);

  static void Allgather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                        void *recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm);

  static void Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                         void *recvbuf, int *recvcounts, int *displs,
                         MPI_Datatype recvtype, MPI_Comm comm);

  static void Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                       void *recvbuf, int recvcount, MPI_Datatype recvtype,
                       MPI_Comm comm);

  static void Sendrecv(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                       int dest, int sendtag, void *recvbuf, int recvcnt,
                       MPI_Datatype recvtype,int source, int recvtag,
                       MPI_Comm comm, MPI_Status *status);

};
#endif

#else //HAVE_MPI

#define MPI_COMM_WORLD 0
#define MPI_UNSIGNED_LONG 1
#define MPI_LONG 2
#define MPI_UNSIGNED_SHORT 3
#define MPI_DOUBLE 4
#define MPI_ANY_SOURCE 5
#define MPI_SUM 6
#define MPI_CHAR 7
#define MPI_SHORT 8
#define MPI_MIN 9
#define MPI_MAX 10
#define MPI_INT 11
class CMPIWrapper {
  
public:
  typedef double Comm;
  typedef int Datatype;
  typedef int Request;
  typedef int Op;
  
  struct Status {
    int MPI_TAG;
    int MPI_SOURCE;
    Status(): MPI_TAG(0), MPI_SOURCE(0){}
  };
    
  static void Init(int *argc, char***argv);
  
  static void Buffer_attach(void *buffer, int size);
  
  static void Buffer_detach(void *buffer, int *size);
  
  static void Finalize();
  
  static void Comm_rank(Comm comm, int* rank);
  
  static void Comm_size(Comm comm, int* size);
  
  static void Barrier(Comm comm);
  
  static void Abort(Comm comm, int error);
  
  static void Get_count(Status *status, Datatype datatype, int *count);

  static void Isend(void *buf, int count, Datatype datatype, int dest,
                    int tag, Comm comm, Request* request);

  static void Irecv(void *buf, int count, Datatype datatype, int source,
                    int tag, Comm comm, Request* request);

  static void Wait(Request *request, Status *status);

  static void Waitall(int nrequests, Request *request, Status *status);

  static void Waitany(int nrequests, Request *request,
                      int *index, Status *status);

  static void Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);

  static void Send(void *buf, int count, Datatype datatype, int dest,
                   int tag, Comm comm);

  static void Recv(void *buf, int count, Datatype datatype, int dest,
                   int tag, Comm comm, Status *status);

  static void Bcast(void *buf, int count, Datatype datatype, int root,
                    Comm comm);

  static void Bsend(void *buf, int count, Datatype datatype, int dest,
                    int tag, Comm comm);

  static void Reduce(void *sendbuf, void *recvbuf, int count,
                     Datatype datatype, Op op, int root, Comm comm);

  static void Allreduce(void *sendbuf, void *recvbuf, int count,
                        Datatype datatype, Op op, Comm comm);

  static void Gather(void *sendbuf, int sendcnt, Datatype sendtype,
                     void *recvbuf, int recvcnt, Datatype recvtype, int root, Comm comm);

  static void Scatter(void *sendbuf, int sendcnt, Datatype sendtype,
                      void *recvbuf, int recvcnt, Datatype recvtype, int root, Comm comm);

  static void Allgather(void *sendbuf, int sendcnt, Datatype sendtype,
                        void *recvbuf, int recvcnt, Datatype recvtype, Comm comm);

  static void Allgatherv(void *sendbuf, int sendcnt, Datatype sendtype,
                         void *recvbuf, int recvcnt, int *displs, Datatype recvtype, Comm comm);

  static void Sendrecv(void *sendbuf, int sendcnt, Datatype sendtype,
                       int dest, int sendtag, void *recvbuf, int recvcnt,
                       Datatype recvtype,int source, int recvtag,
                       Comm comm, Status *status);
  
  static void Alltoall(void *sendbuf, int sendcount, Datatype sendtype,
                           void *recvbuf, int recvcount, Datatype recvtype,
                           Comm comm);
  
  static int Probe(int source, int tag, Comm comm, Status *status);
  
  static void CopyData(void *sendbuf, void *recvbuf, int size, Datatype datatype);
  
};
typedef double SU2_Comm;
typedef CMPIWrapper SU2_MPI;
extern CMPIWrapper::Status* MPI_STATUS_IGNORE;

#endif
#include "mpi_structure.inl"
