/*!
 * \file mpi_structure.hpp
 * \brief Headers of the mpi interface for generalized datatypes.
 *        The subroutines and functions are in the <i>mpi_structure.cpp</i> file.
 * \author T. Albring
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

#pragma once

#ifdef HAVE_MPI
#include "mpi.h"
#include <map>
#endif

#include "./datatype_structure.hpp"
#include <stdlib.h>
#include <unistd.h>

#ifdef HAVE_MPI

/*--- Depending on the datatype used, the correct MPI wrapper class is defined.
 * For the default (double type) case this results in using the normal MPI routines. ---*/
#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE

#include <medi/medi.hpp>
using namespace medi;

class CMediMPIWrapper;
typedef CMediMPIWrapper SU2_MPI;

#if defined CODI_REVERSE_TYPE
#include <codi/externals/codiMediPackTypes.hpp>
#if CODI_PRIMAL_INDEX_TAPE
typedef CoDiPackToolPrimalRestore<su2double> MediTool;
#else
typedef CoDiPackTool<su2double> MediTool;
#endif // defined CODI_REVERSE_TYPE
#elif defined CODI_FORWARD_TYPE
#include <codi/externals/codiForwardMediPackTypes.hpp>
typedef CoDiPackForwardTool<su2double> MediTool;
#endif // defined CODI_FORWARD_TYPE
#define AMPI_ADOUBLE ((medi::MpiTypeInterface*)MediTool::MPI_TYPE)

#else
class CBaseMPIWrapper;
typedef CBaseMPIWrapper SU2_MPI;
#endif // defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE

/*--- Select the appropriate MPI wrapper based on datatype, to use in templated classes. ---*/
template<class T> struct SelectMPIWrapper { typedef SU2_MPI W; };

/*--- In AD we specialize for the passive wrapper. ---*/
#if defined CODI_REVERSE_TYPE
class CBaseMPIWrapper;
template<> struct SelectMPIWrapper<passivedouble> { typedef CBaseMPIWrapper W; };
#endif

/*!
 * \class CMPIWrapper
 * \brief Class for defining the MPI wrapper routines; this class features as a base class for
 * MPI interfaces for non-primitive dataypes e.g. used by AD, complex etc.
 * \author T. Albring
 */

class CBaseMPIWrapper {

public:
  
  typedef MPI_Request  Request;
  typedef MPI_Status   Status;
  typedef MPI_Datatype Datatype;
  typedef MPI_Op       Op;
  typedef MPI_Comm     Comm;
  typedef MPI_Win      Win;
  
protected:
  
  static int Rank, Size, MinRankError;
  static Comm currentComm;
  static bool winMinRankErrorInUse;
  static Win  winMinRankError;
  
public:
  
  static int GetRank();
  
  static int GetSize();
  
  static Comm GetComm();
  
  static void SetComm(Comm NewComm);

  static void Error(std::string ErrorMsg, std::string FunctionName);

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

  static void Testall(int count, Request* array_of_requests, int *flag, Status* array_of_statuses);

  static void Probe(int source, int tag, Comm comm, Status *status);

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

  static void Gather(void *sendbuf, int sendcnt,Datatype sendtype,
                     void *recvbuf, int recvcnt, Datatype recvtype, int root, Comm comm);

  static void Scatter(void *sendbuf, int sendcnt,Datatype sendtype,
                      void *recvbuf, int recvcnt, Datatype recvtype, int root, Comm comm);

  static void Allgather(void *sendbuf, int sendcnt, Datatype sendtype,
                        void *recvbuf, int recvcnt, Datatype recvtype, Comm comm);

  static void Allgatherv(void *sendbuf, int sendcount, Datatype sendtype,
                         void *recvbuf, int *recvcounts, int *displs,
                         Datatype recvtype, Comm comm);

  static void Alltoall(void *sendbuf, int sendcount, Datatype sendtype,
                       void *recvbuf, int recvcount, Datatype recvtype,
                       Comm comm);

  static void Alltoallv(void *sendbuf, int *sendcounts, int *sdispls, Datatype sendtype,
                        void *recvbuf, int *recvcounts, int *recvdispls, Datatype recvtype,
                        Comm comm);

  static void Sendrecv(void *sendbuf, int sendcnt, Datatype sendtype,
                       int dest, int sendtag, void *recvbuf, int recvcnt,
                       Datatype recvtype,int source, int recvtag,
                       Comm comm, Status *status);

  static void Reduce_scatter(void *sendbuf, void *recvbuf, int *recvcounts,
                             Datatype datatype, Op op, Comm comm);
};

typedef MPI_Comm SU2_Comm;

#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE

/*!
 * \class CMediMPIWrapper
 * \brief MPI wrapper functions for MediPack tool.
 * \author T. Albring
 */

class CMediMPIWrapper: public CBaseMPIWrapper {
public:

  typedef AMPI_Request Request;
  typedef AMPI_Status Status;

  static AMPI_Comm     convertComm(MPI_Comm comm);

  static AMPI_Op       convertOp(MPI_Op op);

  static AMPI_Datatype convertDatatype(MPI_Datatype datatype);

  static void SetComm(Comm NewComm);
  
  static void Init(int *argc, char***argv);

  static void Init_AMPI(void);

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

  static void Testall(int count, Request* array_of_requests, int *flag, Status* array_of_statuses);

  static void Probe(int source, int tag, Comm comm, Status *status);

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

  static void Gather(void *sendbuf, int sendcnt,Datatype sendtype,
                     void *recvbuf, int recvcnt, Datatype recvtype, int root, Comm comm);

  static void Scatter(void *sendbuf, int sendcnt,Datatype sendtype,
                      void *recvbuf, int recvcnt, Datatype recvtype, int root, Comm comm);

  static void Allgather(void *sendbuf, int sendcnt, Datatype sendtype,
                        void *recvbuf, int recvcnt, Datatype recvtype, Comm comm);

  static void Allgatherv(void *sendbuf, int sendcount, Datatype sendtype,
                         void *recvbuf, int *recvcounts, int *displs,
                         Datatype recvtype, Comm comm);

  static void Alltoall(void *sendbuf, int sendcount, Datatype sendtype,
                       void *recvbuf, int recvcount, Datatype recvtype,
                       Comm comm);

  static void Alltoallv(void *sendbuf, int *sendcounts, int *sdispls, Datatype sendtype,
                        void *recvbuf, int *recvcounts, int *rdispls, Datatype recvtype,
                        Comm comm);

  static void Sendrecv(void *sendbuf, int sendcnt, Datatype sendtype,
                       int dest, int sendtag, void *recvbuf, int recvcnt,
                       Datatype recvtype,int source, int recvtag,
                       Comm comm, Status *status);

  static void Reduce_scatter(void *sendbuf, void *recvbuf, int *recvcounts,
                             Datatype datatype, Op op, Comm comm);

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
#define MPI_PROD 12
class CBaseMPIWrapper {
  
public:
  typedef int Comm;
  typedef int Datatype;
  typedef int Request;
  typedef int Op;
  
  struct Status {
    int MPI_TAG;
    int MPI_SOURCE;
    Status(): MPI_TAG(0), MPI_SOURCE(0){}
  };

private:
  static int Rank, Size;
  static Comm currentComm;

public:
  static int GetRank();
  
  static int GetSize();  
  
  static Comm GetComm();
  
  static void SetComm(Comm NewComm);
  
  static void Error(std::string ErrorMsg, std::string FunctionName);
    
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

  static void Probe(int source, int tag, Comm comm, Status *status);

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

  static void Alltoallv(void *sendbuf, int *sendcounts, int *sdispls, Datatype sendtype,
                        void *recvbuf, int *recvcounts, int *rdispls, Datatype recvtype,
                        Comm comm);

  static void Reduce_scatter(void *sendbuf, void *recvbuf, int *recvcounts,
                             Datatype datatype, Op op, Comm comm);
    
  static void CopyData(void *sendbuf, void *recvbuf, int size, Datatype datatype);
  
};
typedef int SU2_Comm;
typedef CBaseMPIWrapper SU2_MPI;
extern CBaseMPIWrapper::Status* MPI_STATUS_IGNORE;

#endif

/* Depending on the compiler, define the correct macro to get the current function name */

#if defined(__GNUC__) || (defined(__ICC) && (__ICC >= 600))
# define CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__FUNCSIG__)
# define CURRENT_FUNCTION __FUNCSIG__
#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600))
# define CURRENT_FUNCTION __FUNCTION__
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
# define CURRENT_FUNCTION __func__
#elif defined(__cplusplus) && (__cplusplus >= 201103)
# define CURRENT_FUNCTION __func__
#else
# define CURRENT_FUNCTION "(unknown)"
#endif

#include "mpi_structure.inl"
