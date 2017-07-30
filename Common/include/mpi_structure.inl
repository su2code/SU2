/*!
 * \file mpi_structure.hpp
 * \brief In-Line subroutines of the <i>mpi_structure.hpp</i> file.
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
inline void CMPIWrapper::Init(int *argc, char ***argv) {
  MPI_Init(argc,argv);
}

inline MPI_Request* CMPIWrapper::convertRequest(SU2_MPI::Request* request) {
#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE
  return &request->request;
#else
  return request;
#endif
}

inline void CMPIWrapper::Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count) {
  MPI_Get_count(status, datatype, count);
}

inline void CMPIWrapper::Isend(void *buf, int count, MPI_Datatype datatype,
                               int dest, int tag, MPI_Comm comm, SU2_MPI::Request *request) {
  MPI_Isend(buf,count,datatype,dest,tag,comm,convertRequest(request));
}

inline void CMPIWrapper::Irecv(void *buf, int count, MPI_Datatype datatype,
                               int dest, int tag, MPI_Comm comm, SU2_MPI::Request *request) {
  MPI_Irecv(buf,count,datatype,dest,tag,comm,convertRequest(request));
}

inline void CMPIWrapper::Wait(SU2_MPI::Request *request, MPI_Status *status) {
#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE
  AMPI_Wait(request,status);
#else
  MPI_Wait(request,status);
#endif
}

inline void CMPIWrapper::Waitall(int nrequests, SU2_MPI::Request *request, MPI_Status *status) {
#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE
  AMPI_Waitall(nrequests,request,status);
#else
  MPI_Waitall(nrequests,request,status);
#endif
}

inline void CMPIWrapper::Send(void *buf, int count, MPI_Datatype datatype,
                              int dest, int tag, MPI_Comm comm) {
  MPI_Send(buf,count,datatype,dest,tag,comm);
}

inline void CMPIWrapper::Recv(void *buf, int count, MPI_Datatype datatype,
                              int dest,int tag, MPI_Comm comm, MPI_Status *status) {
  MPI_Recv(buf,count,datatype,dest,tag,comm,status);
}

inline void CMPIWrapper::Bcast(void *buf, int count, MPI_Datatype datatype,
                               int root, MPI_Comm comm) {
  MPI_Bcast(buf,count,datatype,root,comm);
}

inline void CMPIWrapper::Bsend(void *buf, int count, MPI_Datatype datatype,
                               int dest, int tag, MPI_Comm comm) {
  MPI_Bsend(buf,count,datatype,dest,tag,comm);
}

inline void CMPIWrapper::Reduce(void *sendbuf, void *recvbuf, int count,
                                MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
  MPI_Reduce(sendbuf, recvbuf,count,datatype,op,root,comm);
}

inline void CMPIWrapper::Allreduce(void *sendbuf, void *recvbuf, int count,
                                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  MPI_Allreduce(sendbuf,recvbuf,count,datatype,op,comm);
}

inline void CMPIWrapper::Gather(void *sendbuf, int sendcnt,MPI_Datatype sendtype,
                                void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  MPI_Gather(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,root,comm);
}

inline void CMPIWrapper::Scatter(void *sendbuf, int sendcnt,MPI_Datatype sendtype,
                                 void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
}

inline void CMPIWrapper::Allgather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                                   void *recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm) {
  MPI_Allgather(sendbuf,sendcnt,sendtype, recvbuf, recvcnt, recvtype, comm);
}


inline void CMPIWrapper::Sendrecv(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                                  int dest, int sendtag, void *recvbuf, int recvcnt,
                                  MPI_Datatype recvtype,int source, int recvtag,
                                  MPI_Comm comm, MPI_Status *status) {
  MPI_Sendrecv(sendbuf,sendcnt,sendtype,dest,sendtag,recvbuf,recvcnt,recvtype,source,recvtag,comm,status);
}

inline void CMPIWrapper::Waitany(int nrequests, SU2_MPI::Request *request,
                                 int *index, MPI_Status *status) {
#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE
  AMPI_Waitany(nrequests, request, index, status);
#else
  MPI_Waitany(nrequests, request, index, status);
#endif
}
  

#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE


inline void CMediMPIWrapper::Init(int *argc, char ***argv) {
  AMPI_Init(argc,argv);
  MediTool::init();
  AMPI_ADOUBLE = new AMPI_ADOUBLE_TYPE();
}

inline AMPI_Comm CMediMPIWrapper::convertComm(MPI_Comm comm) {
  return AMPI_COMM_WORLD;
}

inline AMPI_Status* CMediMPIWrapper::convertStatus(MPI_Status* status) {
  return status;
}

inline AMPI_Op CMediMPIWrapper::convertOp(MPI_Op op) {
  if(MPI_SUM == op) {
    return MediTool::OP_SUM;
  } else if(MPI_PROD == op) {
    return MediTool::OP_PROD;
  } else if(MPI_MIN == op) {
    return MediTool::OP_MIN;
  } else if(MPI_MAX == op) {
    return MediTool::OP_MAX;
  } else {
    std::cerr << "Implement conversion." << std::endl;
    exit(-1);
    return MediTool::OP_SUM;
  }
}

inline void CMediMPIWrapper::Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count) {
  if (datatype != MPI_DOUBLE) {
    MPI_Get_count(status, datatype, count);
  } else {
    AMPI_Get_count(status, AMPI_ADOUBLE, count);
  }
}


inline void CMediMPIWrapper::Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, SU2_MPI::Request *request) {
  if (datatype != MPI_DOUBLE) {
    MPI_Isend(buf,count,datatype,dest,tag,comm,convertRequest(request));
  } else {
    AMPI_Isend((su2double*)buf,count,AMPI_ADOUBLE,dest,tag,convertComm(comm),request);
  }
}

inline void CMediMPIWrapper::Irecv(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, SU2_MPI::Request *request) {
  if (datatype != MPI_DOUBLE) {
    MPI_Irecv(buf,count,datatype,dest,tag,comm,convertRequest(request));
  } else {
    AMPI_Irecv((su2double*)buf,count,AMPI_ADOUBLE,dest,tag,convertComm(comm),request);
  }
}

inline void CMediMPIWrapper::Wait(SU2_MPI::Request *request, MPI_Status *status) {
  AMPI_Wait(request,convertStatus(status));
}

inline void CMediMPIWrapper::Waitall(int nrequests, SU2_MPI::Request *request, MPI_Status *status) {
  AMPI_Waitall(nrequests,request,convertStatus(status));
}

inline void CMediMPIWrapper::Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Send(buf,count,datatype,dest,tag,comm);
  } else {
    AMPI_Send((su2double*)buf,count,AMPI_ADOUBLE,dest,tag,convertComm(comm));
  }
}

inline void CMediMPIWrapper::Recv(void *buf, int count, MPI_Datatype datatype,int dest,int tag, MPI_Comm comm, MPI_Status *status) {
  if (datatype != MPI_DOUBLE) {
    MPI_Recv(buf,count,datatype,dest,tag,comm,status);
  } else {
    AMPI_Recv((su2double*)buf,count,AMPI_ADOUBLE,dest,tag,convertComm(comm),convertStatus(status));
  }
}

inline void CMediMPIWrapper::Bcast(void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Bcast(buf,count,datatype,root,comm);
  } else {
    AMPI_Bcast((su2double*)buf,count,AMPI_ADOUBLE,root,convertComm(comm));
  }
}

inline void CMediMPIWrapper::Bsend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Bsend(buf,count,datatype,dest,tag,comm);
  } else {
    AMPI_Bsend((su2double*)buf,count,AMPI_ADOUBLE,dest,tag,convertComm(comm));
  }
}

inline void CMediMPIWrapper::Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Reduce(sendbuf, recvbuf,count,datatype,op,root,comm);
  } else {
    AMPI_Reduce((su2double*)sendbuf,(su2double*)recvbuf,count,AMPI_ADOUBLE,convertOp(op),root,convertComm(comm));
  }
}

inline void CMediMPIWrapper::Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Allreduce(sendbuf,recvbuf,count,datatype,op,comm);
  } else {
    AMPI_Allreduce((su2double*)sendbuf,(su2double*)recvbuf,count,AMPI_ADOUBLE,convertOp(op),convertComm(comm));
  }
}

inline void CMediMPIWrapper::Gather(void *sendbuf, int sendcnt,MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  if (sendtype != MPI_DOUBLE) {
    MPI_Gather(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,root,comm);
  } else {
    AMPI_Gather((su2double*)sendbuf,sendcnt,AMPI_ADOUBLE,(su2double*)recvbuf,recvcnt,AMPI_ADOUBLE,root,convertComm(comm));
  }
}

inline void CMediMPIWrapper::Scatter(void *sendbuf, int sendcnt,MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  if (sendtype != MPI_DOUBLE) {
    MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
  }else {
    AMPI_Scatter((su2double*)sendbuf, sendcnt, AMPI_ADOUBLE, (su2double*)recvbuf, recvcnt, AMPI_ADOUBLE, root, convertComm(comm));
  }
}

inline void CMediMPIWrapper::Allgather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm) {
  if (sendtype != MPI_DOUBLE) {
    MPI_Allgather(sendbuf,sendcnt,sendtype, recvbuf, recvcnt, recvtype, comm);
  } else {
    AMPI_Allgather((su2double*)sendbuf,sendcnt,AMPI_ADOUBLE, (su2double*)recvbuf, recvcnt, AMPI_ADOUBLE, convertComm(comm));
  }
}

inline void CMediMPIWrapper::Sendrecv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcnt, MPI_Datatype recvtype,int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
  if(sendtype != MPI_DOUBLE) {
    MPI_Sendrecv(sendbuf,sendcnt,sendtype,dest,sendtag,recvbuf,recvcnt,recvtype,source,recvtag,comm,status);
  } else {
    AMPI_Sendrecv((su2double*)sendbuf,sendcnt,AMPI_ADOUBLE,dest,sendtag,(su2double*)recvbuf,recvcnt,AMPI_ADOUBLE,source,recvtag,convertComm(comm),convertStatus(status));
  }
}

inline void CMediMPIWrapper::Waitany(int nrequests, SU2_MPI::Request *request,
                                        int *index, MPI_Status *status) {

  AMPI_Waitany(nrequests, request, index, convertStatus(status));
}
#endif
#endif
