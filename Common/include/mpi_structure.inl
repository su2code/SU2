/*!
 * \file mpi_structure.hpp
 * \brief In-Line subroutines of the <i>mpi_structure.hpp</i> file.
 * \author T. Albring
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
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

inline void CMPIWrapper::Isend(void *buf, int count, MPI_Datatype datatype,
                               int dest, int tag, MPI_Comm comm, MPI_Request *request) {
  MPI_Isend(buf,count,datatype,dest,tag,comm,request);
}

inline void CMPIWrapper::Irecv(void *buf, int count, MPI_Datatype datatype,
                               int dest, int tag, MPI_Comm comm, MPI_Request *request) {
  MPI_Irecv(buf,count,datatype,dest,tag,comm,request);
}

inline void CMPIWrapper::Wait(MPI_Request *request, MPI_Status *status) {
  MPI_Wait(request,status);
}

inline void CMPIWrapper::Waitall(int nrequests, MPI_Request *request, MPI_Status *status) {
  MPI_Waitall(nrequests,request,status);
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

inline void CMPIWrapper::Waitany(int nrequests, MPI_Request *request,
                                 int *index, MPI_Status *status) {
  MPI_Waitany(nrequests, request, index, status);
}
  
#if defined COMPLEX_TYPE || defined ADOLC_FORWARD_TYPE || defined CODI_FORWARD_TYPE
inline void CAuxMPIWrapper::Allgather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm) {
  if (sendtype != MPI_DOUBLE) {
    MPI_Allgather(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,comm);
  } else {
    int size;
    MPI_Comm_size(comm,&size);
    Gather(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,0,comm);
    Bcast(recvbuf,recvcnt*size,recvtype,0,comm);
  }
}

inline void CAuxMPIWrapper::Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Allreduce(sendbuf,recvbuf,count,datatype,op,comm);
  } else {
    Reduce(sendbuf,recvbuf, count,datatype,op,0,comm);
    Bcast(recvbuf,count,datatype,0,comm);
  }
}

inline void CAuxMPIWrapper::Sendrecv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcnt, MPI_Datatype recvtype,int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
  if (sendtype != MPI_DOUBLE) {
    MPI_Sendrecv(sendbuf,sendcnt,sendtype,dest,sendtag,recvbuf,recvcnt,recvtype,source,recvtag,comm,status);
  } else {
    Bsend(sendbuf,sendcnt,sendtype,dest,sendtag,comm);
    Recv(recvbuf,recvcnt,recvtype,source,recvtag,comm,status);
  }
}

#endif
#ifdef CODI_REVERSE_TYPE
inline void CAdjointMPIWrapper::Init(int *argc, char ***argv) {
  AMPI_Init(argc,argv);
//  buff = new char[1000000];
//  MPI_Buffer_attach(buff, 1000000);
}

inline void CAdjointMPIWrapper::Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
  if (datatype != MPI_DOUBLE) {
    MPI_Isend(buf,count,datatype,dest,tag,comm,request);
  } else {
    AMPI_Isend(buf,count,datatype,dest,tag,comm,request);
  }
}

inline void CAdjointMPIWrapper::Irecv(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
  if (datatype != MPI_DOUBLE) {
    MPI_Irecv(buf,count,datatype,dest,tag,comm,request);
  } else {
    AMPI_Irecv(buf,count,datatype,dest,tag,comm,request);
  }
}

inline void CAdjointMPIWrapper::Wait(MPI_Request *request, MPI_Status *status) {
  AMPI_Wait(request,status);
}

inline void CAdjointMPIWrapper::Waitall(int nrequests, MPI_Request *request, MPI_Status *status) {
  AMPI_Waitall(nrequests,request,status);
}

inline void CAdjointMPIWrapper::Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Send(buf,count,datatype,dest,tag,comm);
  } else {
    AMPI_Send(buf,count,datatype,dest,tag,comm);
  }
}

inline void CAdjointMPIWrapper::Recv(void *buf, int count, MPI_Datatype datatype,int dest,int tag, MPI_Comm comm, MPI_Status *status) {
  if (datatype != MPI_DOUBLE) {
    MPI_Recv(buf,count,datatype,dest,tag,comm,status);
  } else {
    AMPI_Recv(buf,count,datatype,dest,tag,comm,status);
  }
}

inline void CAdjointMPIWrapper::Bcast(void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Bcast(buf,count,datatype,root,comm);
  } else {
    AMPI_Bcast(buf,count,datatype,root,comm);
  }
}

inline void CAdjointMPIWrapper::Bsend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Bsend(buf,count,datatype,dest,tag,comm);
  } else {
    AMPI_Bsend(buf,count,datatype,dest,tag,comm);
  }
}

inline void CAdjointMPIWrapper::Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Reduce(sendbuf, recvbuf,count,datatype,op,root,comm);
  } else {
    AMPI_Reduce(sendbuf,recvbuf,count,datatype,op,root,comm);
  }
}

inline void CAdjointMPIWrapper::Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Allreduce(sendbuf,recvbuf,count,datatype,op,comm);
  } else {
    AMPI_Allreduce(sendbuf,recvbuf,count,datatype,op,comm);
  }
}

inline void CAdjointMPIWrapper::Gather(void *sendbuf, int sendcnt,MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  if (sendtype != MPI_DOUBLE) {
    MPI_Gather(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,root,comm);
  } else {
    AMPI_Gather(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,root,comm);
  }
}

inline void CAdjointMPIWrapper::Scatter(void *sendbuf, int sendcnt,MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) {
  if (sendtype != MPI_DOUBLE) {
    MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
  }else {
    AMPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
  }
}

inline void CAdjointMPIWrapper::Allgather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm) {
  if (sendtype != MPI_DOUBLE) {
    MPI_Allgather(sendbuf,sendcnt,sendtype, recvbuf, recvcnt, recvtype, comm);
  } else {
    int size;
    MPI_Comm_size(comm,&size);
    AMPI_Gather(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,0,comm);
    AMPI_Bcast(recvbuf,recvcnt*size,recvtype,0,comm);
  }
}

inline void CAdjointMPIWrapper::Sendrecv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcnt, MPI_Datatype recvtype,int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
  if(sendtype != MPI_DOUBLE) {
    MPI_Sendrecv(sendbuf,sendcnt,sendtype,dest,sendtag,recvbuf,recvcnt,recvtype,source,recvtag,comm,status);
  } else {
    AMPI_Sendrecv(sendbuf,sendcnt,sendtype,dest,sendtag,recvbuf,recvcnt,recvtype,source,recvtag,comm,status);
  }
}

inline void CAdjointMPIWrapper::Waitany(int nrequests, MPI_Request *request,
                                        int *index, MPI_Status *status) {
  AMPI_Waitany(nrequests, request, index, status);
}
#endif
#endif
