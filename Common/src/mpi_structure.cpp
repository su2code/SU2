/*!
 * \file mpi_structure.cpp
 * \brief Main subroutines for the mpi structures.
 * \author T. Albring
 * \version 4.2.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/mpi_structure.hpp"

#ifdef HAVE_MPI

#if defined COMPLEX_TYPE || defined ADOLC_FORWARD_TYPE || defined CODI_FORWARD_TYPE
std::map<MPI_Request*, CAuxMPIWrapper::CommInfo>
CAuxMPIWrapper::CommInfoMap;
std::map<MPI_Request*, CAuxMPIWrapper::CommInfo>::iterator
CAuxMPIWrapper::CommInfoIterator;

void CAuxMPIWrapper::Isend(void *buf, int count, MPI_Datatype datatype,
                               int dest, int tag, MPI_Comm comm, MPI_Request *request) {
  if (datatype != MPI_DOUBLE) {
    MPI_Isend(buf,count,datatype,dest,tag,comm,request);
  } else {

    int iVal;

    /*--- Create request object for Complex communication ---*/

    MPI_Request* RequestAux = new MPI_Request;

    /*--- Create buffer objects (Note: they will be deleted in the wait routine!) ---*/

    double *ValueBuffer = new double[count];
    double *AuxBuffer   = new double[count];
    su2double *SendBuffer = static_cast<su2double*>(buf);


    /*--- Extract real value and imag value ---*/

    for (iVal = 0; iVal < count; iVal++) {
      ValueBuffer[iVal] = SU2_TYPE::GetValue(SendBuffer[iVal]);
      AuxBuffer[iVal]   = SU2_TYPE::GetSecondary(SendBuffer[iVal]);
    }

    /*---  Send real value and imag value ---*/

    MPI_Isend(ValueBuffer,count,datatype,dest,tag,comm,request);
    MPI_Isend(AuxBuffer,count,datatype,dest,tag+100,comm,RequestAux);

    /*--- Create info object for wait routine to find the request for the aux var ---*/

    CommInfo info;
    info.ValueBuffer      = ValueBuffer;
    info.AuxBuffer        = AuxBuffer;
    info.RequestAux       = RequestAux;
    info.su2doubleBuffer  = SendBuffer;
    info.count            = count;
    info.Type             = ISEND;

    /*--- Insert info object into global map -- */

    CommInfoMap.insert(std::pair<MPI_Request*, CommInfo>(request,info));

  }
}


void CAuxMPIWrapper::Irecv(void *buf, int count, MPI_Datatype datatype,
                               int source, int tag, MPI_Comm comm, MPI_Request *request) {
  if (datatype != MPI_DOUBLE) {
    MPI_Irecv(buf,count,datatype,source,tag,comm,request);
  } else {

    /*--- Create request object for Complex communication ---*/

    MPI_Request* RequestAux = new MPI_Request;

    /*--- Create buffer objects (Note: they will be deleted in the wait routine!) ---*/

    double *ValueBuffer = new double[count];
    double *AuxBuffer   = new double[count];
    su2double *RecvBuffer = static_cast<su2double*>(buf);

    /*---  Recv real value and imag value ---*/

    MPI_Irecv(ValueBuffer,count,datatype,source,tag,comm,request);
    MPI_Irecv(AuxBuffer,count,datatype,source,tag+100,comm,RequestAux);

    /*--- Create info object for wait routine to find the request for the aux var ---*/

    CommInfo info;
    info.ValueBuffer      = ValueBuffer;
    info.AuxBuffer        = AuxBuffer;
    info.RequestAux       = RequestAux;
    info.su2doubleBuffer  = RecvBuffer;
    info.count            = count;
    info.Type             = IRECV;

    /*--- Insert info object into global map -- */

    CommInfoMap.insert(std::pair<MPI_Request*, CommInfo>(request,info));

  }
}

void CAuxMPIWrapper::Wait(MPI_Request *request, MPI_Status *status) {

  /*--- First wait for send/recv of normal value operation to finish ---*/

  MPI_Wait(request, status);

  /* Search for request in case there is also a aux. request and finalize this comm. ---*/

  if((CommInfoIterator = CommInfoMap.find(request)) != CommInfoMap.end()) {
    FinalizeCommunication(CommInfoIterator);
  }
}


void CAuxMPIWrapper::FinalizeCommunication(
  std::map<MPI_Request *, CommInfo>::iterator &CommInfoIt) {

  /*--- Get info about communication ---*/

  CommInfo info = CommInfoIt->second;
  double *ValueBuffer   = info.ValueBuffer;
  double *AuxBuffer     = info.AuxBuffer;
  su2double *Buffer     = info.su2doubleBuffer;
  MPI_Status status;
  MPI_Request *RequestAux = info.RequestAux;
  int iVal, count = info.count;

  /*--- Wait for aux. request ---*/

  MPI_Wait(RequestAux,&status);

  switch(info.Type) {
  case ISEND:
    break;
  case IRECV:
    for (iVal = 0; iVal < count; iVal++) {
      SU2_TYPE::SetValue(Buffer[iVal], ValueBuffer[iVal]);
      SU2_TYPE::SetSecondary(Buffer[iVal], AuxBuffer[iVal]);
    }
    break;
  }

  /*--- Delete buffer ---*/

  delete [] ValueBuffer;
  delete [] AuxBuffer;

  /*--- Delete Request ---*/

  delete RequestAux;
  CommInfoMap.erase(CommInfoIt);

}

void CAuxMPIWrapper::Waitall(int nrequests, MPI_Request *request,
                                 MPI_Status *status) {

  /*--- Wait for normal requests to finish ---*/

  MPI_Waitall(nrequests, request, status);

  /*--- Wait for aux. requests and finish communication ---*/

  for (int iVal = 0; iVal < nrequests; iVal++) {
    if((CommInfoIterator = CommInfoMap.find(&request[iVal])) != CommInfoMap.end()) {
      FinalizeCommunication(CommInfoIterator);
    }
  }
}

void CAuxMPIWrapper::Waitany(int nrequests, MPI_Request *request,
                             int *index, MPI_Status *status){

  /*--- Wait for any normal request to finish ---*/

  MPI_Waitany(nrequests, request, index, status);

  /*--- Wait for particular aux. request and finish communication ---*/

  if((CommInfoIterator = CommInfoMap.find(&request[*index])) != CommInfoMap.end()) {
    FinalizeCommunication(CommInfoIterator);
  }

}
void CAuxMPIWrapper::Send(void *buf, int count, MPI_Datatype datatype,
                              int dest, int tag, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Send(buf,count,datatype,dest,tag,comm);
  } else {

    double *AuxBuffer      = new double[count];
    double *ValueBuffer    = new double[count];
    su2double *SendBuffer  = static_cast<su2double*>(buf);

    int iVal;

    for (iVal = 0; iVal < count; iVal++) {
      ValueBuffer[iVal] = SU2_TYPE::GetValue(SendBuffer[iVal]);
      AuxBuffer[iVal]   = SU2_TYPE::GetSecondary(SendBuffer[iVal]);
    }

    MPI_Send(ValueBuffer,count,datatype,dest,tag,comm);
    MPI_Send(AuxBuffer,count,datatype,dest,tag+100,comm);

    delete [] ValueBuffer;
    delete [] AuxBuffer;
  }
}

void CAuxMPIWrapper::Recv(void *buf, int count, MPI_Datatype datatype,
                              int dest, int tag, MPI_Comm comm, MPI_Status* status) {
  if (datatype != MPI_DOUBLE) {
    MPI_Recv(buf,count,datatype,dest,tag,comm,status);
  } else {
    double *AuxBuffer      = new double[count];
    double *ValueBuffer    = new double[count];
    su2double *RecvBuffer  = static_cast<su2double*>(buf);

    int iVal;

    MPI_Recv(ValueBuffer,count,datatype,dest,tag,comm,status);
    MPI_Recv(AuxBuffer,count,datatype,dest,tag+100,comm,status);

    for (iVal = 0; iVal < count; iVal++) {
      SU2_TYPE::SetValue(RecvBuffer[iVal], ValueBuffer[iVal]);
      SU2_TYPE::SetSecondary(RecvBuffer[iVal], AuxBuffer[iVal]);
    }

    delete [] AuxBuffer;
    delete [] ValueBuffer;
  }
}

void CAuxMPIWrapper::Bsend(void *buf, int count, MPI_Datatype datatype,
                               int dest, int tag, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Bsend(buf,count,datatype,dest,tag,comm);
  } else {

    double *AuxBuffer      = new double[count];
    double *ValueBuffer    = new double[count];
    su2double *SendBuffer  = static_cast<su2double*>(buf);

    int iVal;

    for (iVal = 0; iVal < count; iVal++) {
      ValueBuffer[iVal] = SU2_TYPE::GetValue(SendBuffer[iVal]);
      AuxBuffer[iVal]   = SU2_TYPE::GetSecondary(SendBuffer[iVal]);
    }

    MPI_Bsend(ValueBuffer,count,datatype,dest,tag,comm);
    MPI_Bsend(AuxBuffer,count,datatype,dest,tag+100,comm);

    delete [] ValueBuffer;
    delete [] AuxBuffer;
  }
}

void CAuxMPIWrapper::Reduce(void *sendbuf, void *recvbuf, int count,
                                MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Reduce(sendbuf, recvbuf,count,datatype,op,root,comm);
  } else {
    int rank;
    MPI_Comm_rank(comm, &rank);
    su2double* SendBuffer = static_cast< su2double*>(sendbuf);
    su2double* RecvBuffer = static_cast< su2double*>(recvbuf);

    double  *SendAuxBuffer = NULL, *SendValueBuffer = NULL,
            *RecvAuxBuffer = NULL, *RecvValueBuffer = NULL;

    int iVal = 0;

    if (op == MPI_SUM) {
      SendAuxBuffer = new double[count];
      SendValueBuffer = new double[count];

      if (rank == root) {
        RecvAuxBuffer = new double[count];
        RecvValueBuffer = new double[count];
      }

      for (iVal = 0; iVal < count; iVal++) {
        SendValueBuffer[iVal] = SU2_TYPE::GetValue(SendBuffer[iVal]);
        SendAuxBuffer[iVal]   = SU2_TYPE::GetSecondary(SendBuffer[iVal]);
      }

      MPI_Reduce(SendValueBuffer, RecvValueBuffer, count, datatype, op, root, comm);
      MPI_Reduce(SendAuxBuffer, RecvAuxBuffer, count, datatype, op, root,
                 comm);
      if (rank == root) {
        for (iVal = 0; iVal < count; iVal++) {
            SU2_TYPE::SetValue(RecvBuffer[iVal], RecvValueBuffer[iVal]);
            SU2_TYPE::SetSecondary(RecvBuffer[iVal], RecvAuxBuffer[iVal]);
        }
      }

      delete [] SendValueBuffer;
      delete [] SendAuxBuffer;
      if (rank == root) {
        delete [] RecvValueBuffer;
        delete [] RecvAuxBuffer;
      }
    } else if(op == MPI_MAX || op == MPI_MIN) {
      ValRank* SendValLoc = new ValRank[count];
      ValRank* RecvValLoc = new ValRank[count];
      double temp = 0;

      for (iVal = 0; iVal < count; iVal++) {
        SendValLoc[iVal].value = SU2_TYPE::GetValue(SendBuffer[iVal]);
        SendValLoc[iVal].rank = rank;
      }

      if (op == MPI_MAX)
        MPI_Allreduce(SendValLoc, RecvValLoc, count, MPI_DOUBLE_INT,MPI_MAXLOC,comm);
      else if (op == MPI_MIN)
        MPI_Allreduce(SendValLoc, RecvValLoc, count, MPI_DOUBLE_INT,MPI_MINLOC,comm);

      MPI_Status status;
      for (iVal = 0; iVal < count; iVal++) {
        if (rank == RecvValLoc[iVal].rank) {
          temp = SU2_TYPE::GetSecondary(SendBuffer[iVal]);
          MPI_Bsend(&temp ,1 , MPI_DOUBLE, root, rank, comm);
        }
        if (rank == root) {
          MPI_Recv(&temp,1, MPI_DOUBLE, RecvValLoc[iVal].rank, RecvValLoc[iVal].rank,
                   comm,
                   &status);
          SU2_TYPE::SetValue(RecvBuffer[iVal], RecvValLoc[iVal].value);
          SU2_TYPE::SetSecondary(RecvBuffer[iVal], temp);
        }
      }

      delete [] SendValLoc;
      delete [] RecvValLoc;
    } else {
      if (rank == root)
        std::cout << "Reduce operation not implemented for this kind of operation" <<
                  std::endl;
      MPI_Abort(comm,1);
    }
  }
}

void CAuxMPIWrapper::Gather(void *sendbuf, int sendcnt,
                                MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype,
                                int root, MPI_Comm comm) {
  if (sendtype != MPI_DOUBLE) {
    MPI_Gather(sendbuf,sendcnt,sendtype, recvbuf, recvcnt, recvtype,root,comm);
  } else {

    double* SendValueBuffer = new double[sendcnt];
    double* SendAuxBuffer   = new double[sendcnt];

    int iVal;

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    double *RecvValueBuffer = NULL, *RecvAuxBuffer = NULL;

    if (rank == root) {
      RecvValueBuffer = new double[recvcnt*size];
      RecvAuxBuffer   = new double[recvcnt*size];
    }

    su2double *SendBuffer = static_cast< su2double* >(sendbuf);
    su2double *RecvBuffer = static_cast< su2double* >(recvbuf);

    for (iVal = 0; iVal < sendcnt; iVal++) {
      SendValueBuffer[iVal] = SU2_TYPE::GetValue(SendBuffer[iVal]);
      SendAuxBuffer[iVal]   = SU2_TYPE::GetSecondary(SendBuffer[iVal]);
    }

    MPI_Gather(SendValueBuffer, sendcnt, sendtype, RecvValueBuffer, recvcnt,
               recvtype, root, comm);
    MPI_Gather(SendAuxBuffer, sendcnt, sendtype, RecvAuxBuffer, recvcnt,
               recvtype, root, comm);

    if (rank == root) {
      for (iVal = 0; iVal < recvcnt*size; iVal++) {
        SU2_TYPE::SetValue(RecvBuffer[iVal],  RecvValueBuffer[iVal]);
        SU2_TYPE::SetSecondary(RecvBuffer[iVal], RecvAuxBuffer[iVal]);
      }
      delete [] RecvValueBuffer;
      delete [] RecvAuxBuffer;
    }

    delete [] SendValueBuffer;
    delete [] SendAuxBuffer;
  }
}

void CAuxMPIWrapper::Scatter(void *sendbuf, int sendcnt,
                            MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype,
                            int root, MPI_Comm comm) {
  if (sendtype != MPI_DOUBLE) {
    MPI_Scatter(sendbuf,sendcnt,sendtype, recvbuf, recvcnt, recvtype,root,comm);
  } else {

    double* SendValueBuffer = NULL;
    double* SendAuxBuffer   = NULL;

    int iVal;

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    double *RecvValueBuffer = NULL, *RecvAuxBuffer = NULL;

    RecvValueBuffer = new double[recvcnt*size];
    RecvAuxBuffer   = new double[recvcnt*size];

    su2double *SendBuffer =  static_cast< su2double* >(sendbuf);
    su2double *RecvBuffer =  static_cast< su2double* >(recvbuf);

    if (rank == root){
      SendValueBuffer = new double[sendcnt];
      SendAuxBuffer    = new double[sendcnt];

      for (iVal = 0; iVal < sendcnt; iVal++) {
        SendValueBuffer[iVal] = SU2_TYPE::GetValue(SendBuffer[iVal]);
        SendAuxBuffer[iVal]   = SU2_TYPE::GetSecondary(SendBuffer[iVal]);
      }
    }

    MPI_Scatter(SendValueBuffer, sendcnt, sendtype, RecvValueBuffer, recvcnt,
               recvtype, root, comm);
    MPI_Scatter(SendAuxBuffer, sendcnt, sendtype, RecvAuxBuffer, recvcnt,
               recvtype, root, comm);

    for (iVal = 0; iVal < recvcnt*size; iVal++) {
      SU2_TYPE::SetValue(RecvBuffer[iVal],  RecvValueBuffer[iVal]);
      SU2_TYPE::SetSecondary(RecvBuffer[iVal], RecvAuxBuffer[iVal]);
    }
    delete [] RecvValueBuffer;
    delete [] RecvAuxBuffer;

    if (rank == root){
      delete [] SendValueBuffer;
      delete [] SendAuxBuffer;
    }
  }
}

void CAuxMPIWrapper::Bcast(void *buf, int count, MPI_Datatype datatype,
                               int root, MPI_Comm comm) {
  if (datatype != MPI_DOUBLE) {
    MPI_Bcast(buf,count,datatype,root,comm);
  } else {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    double *AuxBuffer = new double[count];
    double *ValueBuffer = new double[count];

    int iVal = 0;

    su2double *Buffer = static_cast< su2double* >(buf);

    if (rank == root) {
      for (iVal = 0; iVal < count; iVal++) {
        ValueBuffer[iVal] = SU2_TYPE::GetValue(Buffer[iVal]);
        AuxBuffer[iVal]   = SU2_TYPE::GetSecondary(Buffer[iVal]);
      }
    }

    MPI_Bcast(ValueBuffer, count, datatype, root, comm);
    MPI_Bcast(AuxBuffer, count, datatype, root, comm);

    if (rank != root) {
      for (iVal = 0; iVal < count; iVal++) {
        SU2_TYPE::SetValue(Buffer[iVal], ValueBuffer[iVal]);
        SU2_TYPE::SetSecondary(Buffer[iVal], AuxBuffer[iVal]);
      }
    }

    delete [] ValueBuffer;
    delete [] AuxBuffer;
  }
}
#endif
#ifdef CODI_REVERSE_TYPE
#define AD_TYPE su2double
#include "externals/ampi_interface_realreverse_old.cpp"
#undef AD_TYPE
#endif
#endif
