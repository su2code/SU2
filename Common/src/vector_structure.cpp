/*!
 * \file vector_structure.cpp
 * \brief Main classes required for solving linear systems of equations
 * \author F. Palacios, J. Hicken
 * \version 6.0.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include "../include/vector_structure.hpp"

template<class CalcType>
TCSysVector<CalcType>::TCSysVector(void) {
  
  vec_val = NULL;
  
}

template<class CalcType>
TCSysVector<CalcType>::TCSysVector(const unsigned long & size, const CalcType & val) {
  
  nElm = size; nElmDomain = size;
  nBlk = nElm; nBlkDomain = nElmDomain;
  nVar = 1;
  
  /*--- Check for invalid size, then allocate memory and initialize values ---*/
  if ( (nElm <= 0) || (nElm >= UINT_MAX) ) {
    char buf[100];
    SPRINTF(buf, "Invalid input: size = %lu", size );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  vec_val = new CalcType[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = val;
  
#ifdef HAVE_MPI
  unsigned long nElmLocal = (unsigned long)nElm;
  SU2_MPI::Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
}

template<class CalcType>
TCSysVector<CalcType>::TCSysVector(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar,
                       const CalcType & val) {

  nElm = numBlk*numVar; nElmDomain = numBlkDomain*numVar;
  nBlk = numBlk; nBlkDomain = numBlkDomain;
  nVar = numVar;
  
  /*--- Check for invalid input, then allocate memory and initialize values ---*/
  if ( (nElm <= 0) || (nElm >= ULONG_MAX) ) {
    char buf[100];
    SPRINTF(buf, "invalid input: numBlk, numVar = %lu, %u", numBlk, numVar );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }
	
  vec_val = new CalcType[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = val;
  
#ifdef HAVE_MPI
  unsigned long nElmLocal = (unsigned long)nElm;
  SU2_MPI::Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
}

template<class CalcType>
TCSysVector<CalcType>::TCSysVector(const TCSysVector<CalcType> & u) {
  
  /*--- Copy size information, allocate memory, and initialize values ---*/
  nElm = u.nElm; nElmDomain = u.nElmDomain;
  nBlk = u.nBlk; nBlkDomain = u.nBlkDomain;
  nVar = u.nVar;
  
  vec_val = new CalcType[nElm];
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = u.vec_val[i];
  
#ifdef HAVE_MPI
  nElmGlobal = u.nElmGlobal;
#endif
  
}

template<class CalcType>
TCSysVector<CalcType>::TCSysVector(const unsigned long & size, const CalcType* u_array) {
  
  nElm = size; nElmDomain = size;
  nBlk = nElm; nBlkDomain = nElmDomain;
  nVar = 1;
  
  /*--- Check for invalid size, then allocate memory and initialize values ---*/
  if ( (nElm <= 0) || (nElm >= ULONG_MAX) ) {
    char buf[100];    
    SPRINTF(buf, "Invalid input: size = %lu", size );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  vec_val = new CalcType[nElm];
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = u_array[i];

#ifdef HAVE_MPI
  unsigned long nElmLocal = (unsigned long)nElm;
  SU2_MPI::Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
}

template<class CalcType>
TCSysVector<CalcType>::TCSysVector(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar,
                       const CalcType* u_array) {

  nElm = numBlk*numVar; nElmDomain = numBlkDomain*numVar;
  nBlk = numBlk; nBlkDomain = numBlkDomain;
  nVar = numVar;
  
  /*--- check for invalid input, then allocate memory and initialize values ---*/
  if ( (nElm <= 0) || (nElm >= ULONG_MAX) ) {
    char buf[100];
    SPRINTF(buf, "invalid input: numBlk, numVar = %lu, %u", numBlk, numVar );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  vec_val = new CalcType[nElm];
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = u_array[i];

#ifdef HAVE_MPI
  unsigned long nElmLocal = (unsigned long)nElm;
  SU2_MPI::Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
}

template<class CalcType>
TCSysVector<CalcType>::~TCSysVector() {
  delete [] vec_val;
  
  nElm = 0; nElmDomain = 0;
  nBlk = 0; nBlkDomain = 0;
  nVar = 0;
  
}

template<class CalcType>
void TCSysVector<CalcType>::Initialize(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const CalcType & val) {
  
  nElm = numBlk*numVar; nElmDomain = numBlkDomain*numVar;
  nBlk = numBlk; nBlkDomain = numBlkDomain;
  nVar = numVar;
  
  /*--- Check for invalid input, then allocate memory and initialize values ---*/
  if ( (nElm <= 0) || (nElm >= ULONG_MAX) ) {
    char buf[100];
    SPRINTF(buf, "invalid input: numBlk, numVar = %lu, %u", numBlk, numVar );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }
	
  vec_val = new CalcType[nElm];
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = val;
  
#ifdef HAVE_MPI
  unsigned long nElmLocal = (unsigned long)nElm;
  SU2_MPI::Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
}

template<class CalcType>
void TCSysVector<CalcType>::Equals_AX(const CalcType & a, TCSysVector<CalcType> & x) {
  /*--- check that *this and x are compatible ---*/
  if (nElm != x.nElm) {
    cerr << "TCSysVector<CalcType>::Equals_AX(): " << "sizes do not match";
    throw(-1);
  }
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = a * x.vec_val[i];
}

template<class CalcType>
void TCSysVector<CalcType>::Plus_AX(const CalcType & a, TCSysVector<CalcType> & x) {
  /*--- check that *this and x are compatible ---*/
  if (nElm != x.nElm) {
    SU2_MPI::Error("Sizes do not match", CURRENT_FUNCTION);
  }
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] += a * x.vec_val[i];
}

template<class CalcType>
void TCSysVector<CalcType>::Equals_AX_Plus_BY(const CalcType & a, TCSysVector<CalcType> & x, const CalcType & b, TCSysVector<CalcType> & y) {
  /*--- check that *this, x and y are compatible ---*/
  if ((nElm != x.nElm) || (nElm != y.nElm)) {
    SU2_MPI::Error("Sizes do not match", CURRENT_FUNCTION);
  }
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = a * x.vec_val[i] + b * y.vec_val[i];
}

template<class CalcType>
TCSysVector<CalcType> & TCSysVector<CalcType>::operator=(const TCSysVector<CalcType> & u) {
  
  /*--- check if self-assignment, otherwise perform deep copy ---*/
  if (this == &u) return *this;
  
  if (nElm != u.nElm){
    SU2_MPI::Error("Sizes do not match", CURRENT_FUNCTION);
  }
  
  nBlk = u.nBlk;
	nBlkDomain = u.nBlkDomain;
  
  nVar = u.nVar;
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = u.vec_val[i];
  
#ifdef HAVE_MPI
  nElmGlobal = u.nElmGlobal;
#endif
  
  return *this;
}

template<class CalcType>
TCSysVector<CalcType> & TCSysVector<CalcType>::operator=(const CalcType & val) {
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = val;
  return *this;
}

template<class CalcType>
TCSysVector<CalcType> TCSysVector<CalcType>::operator+(const TCSysVector<CalcType> & u) const {
  
  /*--- Use copy constructor and compound addition-assignment ---*/
  TCSysVector<CalcType> sum(*this);
  sum += u;
  return sum;
}

template<class CalcType>
TCSysVector<CalcType> & TCSysVector<CalcType>::operator+=(const TCSysVector<CalcType> & u) {
  
  /*--- Check for consistent sizes, then add elements ---*/
  if (nElm != u.nElm) {
    SU2_MPI::Error("Sizes do not match", CURRENT_FUNCTION);
  }
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] += u.vec_val[i];
  return *this;
}

template<class CalcType>
TCSysVector<CalcType> TCSysVector<CalcType>::operator-(const TCSysVector<CalcType> & u) const {
  
  /*--- Use copy constructor and compound subtraction-assignment ---*/
  TCSysVector<CalcType> diff(*this);
  diff -= u;
  return diff;
}

template<class CalcType>
TCSysVector<CalcType> & TCSysVector<CalcType>::operator-=(const TCSysVector<CalcType> & u) {
  
  /*--- Check for consistent sizes, then subtract elements ---*/
  if (nElm != u.nElm) {
    SU2_MPI::Error("Sizes do not match", CURRENT_FUNCTION);
  }
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] -= u.vec_val[i];
  return *this;
}

template<class CalcType>
TCSysVector<CalcType> TCSysVector<CalcType>::operator*(const CalcType & val) const {
  
  /*--- use copy constructor and compound scalar
   multiplication-assignment ---*/
  TCSysVector<CalcType> prod(*this);
  prod *= val;
  return prod;
}

template<class CalcType>
TCSysVector<CalcType> operator*(const CalcType & val, const CSysVector & u) {
  
  /*--- use copy constructor and compound scalar
   multiplication-assignment ---*/
  TCSysVector<CalcType> prod(u);
  prod *= val;
  return prod;
}

template<class CalcType>
TCSysVector<CalcType> & TCSysVector<CalcType>::operator*=(const CalcType & val) {
  
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] *= val;
  return *this;
}

template<class CalcType>
TCSysVector<CalcType> TCSysVector<CalcType>::operator/(const CalcType & val) const {
  
  /*--- use copy constructor and compound scalar
   division-assignment ---*/
  TCSysVector<CalcType> quotient(*this);
  quotient /= val;
  return quotient;
}

template<class CalcType>
TCSysVector<CalcType> & TCSysVector<CalcType>::operator/=(const CalcType & val) {
  
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] /= val;
  return *this;
}

template<class CalcType>
CalcType TCSysVector<CalcType>::norm() const {
  
  /*--- just call dotProd on this*, then sqrt ---*/
  CalcType val = dotProd(*this, *this);
  if (val < 0.0) {
    SU2_MPI::Error("Inner product of CSysVector is negative", CURRENT_FUNCTION);
  }
  return sqrt(val);
}

template<class CalcType>
void TCSysVector<CalcType>::CopyToArray(CalcType* u_array) {
  
  for (unsigned long i = 0; i < nElm; i++)
    u_array[i] = vec_val[i];
}

template<class CalcType>
void TCSysVector<CalcType>::AddBlock(unsigned long val_ipoint, CalcType *val_residual) {
  unsigned short iVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    vec_val[val_ipoint*nVar+iVar] += val_residual[iVar];
}

template<class CalcType>
void TCSysVector<CalcType>::SubtractBlock(unsigned long val_ipoint, CalcType *val_residual) {
  unsigned short iVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    vec_val[val_ipoint*nVar+iVar] -= val_residual[iVar];
}

template<class CalcType>
void TCSysVector<CalcType>::SetBlock(unsigned long val_ipoint, CalcType *val_residual) {
  unsigned short iVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    vec_val[val_ipoint*nVar+iVar] = val_residual[iVar];
}

template<class CalcType>
void TCSysVector<CalcType>::SetBlock(unsigned long val_ipoint, unsigned short val_var, CalcType val_residual) {

  vec_val[val_ipoint*nVar+val_var] = val_residual;
}

template<class CalcType> 
void TCSysVector<CalcType>::SetBlock_Zero(unsigned long val_ipoint) {
  unsigned short iVar;

  for (iVar = 0; iVar < nVar; iVar++)
    vec_val[val_ipoint*nVar+iVar] = 0.0;
}

template<class CalcType>
void TCSysVector<CalcType>::SetBlock_Zero(unsigned long val_ipoint, unsigned short val_var) {
    vec_val[val_ipoint*nVar+val_var] = 0.0;
}

template<class CalcType>
CalcType TCSysVector<CalcType>::GetBlock(unsigned long val_ipoint, unsigned short val_var) {
  return vec_val[val_ipoint*nVar + val_var];
}

template<class CalcType>
CalcType *TCSysVector<CalcType>::GetBlock(unsigned long val_ipoint) {
  return &vec_val[val_ipoint*nVar];
}

template<class CalcType>
CalcType dotProd(const TCSysVector<CalcType> & u, const TCSysVector<CalcType> & v) {
  
  /*--- check for consistent sizes ---*/
  if (u.nElm != v.nElm) {
    SU2_MPI::Error("Sizes do not match", CURRENT_FUNCTION);
  }
  
  /*--- find local inner product and, if a parallel run, sum over all
   processors (we use nElemDomain instead of nElem) ---*/
  CalcType loc_prod = 0.0;
  for (unsigned long i = 0; i < u.nElmDomain; i++)
    loc_prod += u.vec_val[i]*v.vec_val[i];
  CalcType prod = 0.0;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&loc_prod, &prod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  prod = loc_prod;
#endif
  
  return prod;
}

template<>
void TCSysVector<su2double>::SendReceive(CGeometry *geometry, CConfig *config) {
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive = NULL, *Buffer_Send = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
      
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      
      Buffer_Receive = new su2double [nBufferR_Vector];
      Buffer_Send = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sended ---*/
      
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send[iVertex*nVar+iVar] = vec_val[iPoint*nVar+iVar];
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      
      SU2_MPI::Sendrecv(Buffer_Send, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive[iVar*nVertexR+iVertex] = Buffer_Send[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      
      delete [] Buffer_Send;
      
      /*--- Do the coordinate transformation ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        
        for (iVar = 0; iVar < nVar; iVar++)
          vec_val[iPoint*nVar+iVar] = Buffer_Receive[iVertex*nVar+iVar];
        
      }
      
      /*--- Deallocate receive buffer ---*/
      
      delete [] Buffer_Receive;
      
    }
    
  }
  
}


template<>
void TCSysVector<su2double>::SendReceive_Reverse(CGeometry *geometry, CConfig *config) {

  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive = NULL, *Buffer_Send = NULL;

#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker + 1;  MarkerR = iMarker;

#ifdef HAVE_MPI

      receive_from = config->GetMarker_All_SendRecv(MarkerR)-1;
      send_to = abs(config->GetMarker_All_SendRecv(MarkerS))-1;

#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/

      Buffer_Receive = new su2double [nBufferR_Vector];
      Buffer_Send = new su2double[nBufferS_Vector];

      /*--- Copy the solution that should be sended ---*/

      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send[iVertex*nVar+iVar] = vec_val[iPoint*nVar+iVar];
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/

      SU2_MPI::Sendrecv(Buffer_Send, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/

      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive[iVar*nVertexR+iVertex] = Buffer_Send[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/

      delete [] Buffer_Send;

      /*--- Do the coordinate transformation ---*/

      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/

        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

        /*--- Copy transformed conserved variables back into buffer. ---*/

        for (iVar = 0; iVar < nVar; iVar++)
          vec_val[iPoint*nVar+iVar] += Buffer_Receive[iVertex*nVar+iVar];

      }

      /*--- Deallocate receive buffer ---*/

      delete [] Buffer_Receive;

    }

  }

}

#ifdef CODI_REVERSE_TYPE
template<>
void TCSysVector<passivedouble>::SendReceive(CGeometry *geometry, CConfig *config) {
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  passivedouble *Buffer_Receive = NULL, *Buffer_Send = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
      
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      
      Buffer_Receive = new passivedouble [nBufferR_Vector];
      Buffer_Send = new passivedouble[nBufferS_Vector];
      
      /*--- Copy the solution that should be sended ---*/
      
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send[iVertex*nVar+iVar] = vec_val[iPoint*nVar+iVar];
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      
      CBaseMPIWrapper::Sendrecv(Buffer_Send, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive[iVar*nVertexR+iVertex] = Buffer_Send[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      
      delete [] Buffer_Send;
      
      /*--- Do the coordinate transformation ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        
        for (iVar = 0; iVar < nVar; iVar++)
          vec_val[iPoint*nVar+iVar] = Buffer_Receive[iVertex*nVar+iVar];
        
      }
      
      /*--- Deallocate receive buffer ---*/
      
      delete [] Buffer_Receive;
      
    }
    
  }
  
}



template<>
void TCSysVector<passivedouble>::SendReceive_Reverse(CGeometry *geometry, CConfig *config) {

  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  passivedouble *Buffer_Receive = NULL, *Buffer_Send = NULL;

#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker + 1;  MarkerR = iMarker;

#ifdef HAVE_MPI

      receive_from = config->GetMarker_All_SendRecv(MarkerR)-1;
      send_to = abs(config->GetMarker_All_SendRecv(MarkerS))-1;

#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/

      Buffer_Receive = new passivedouble [nBufferR_Vector];
      Buffer_Send = new passivedouble[nBufferS_Vector];

      /*--- Copy the solution that should be sended ---*/

      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send[iVertex*nVar+iVar] = vec_val[iPoint*nVar+iVar];
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/

      CBaseMPIWrapper::Sendrecv(Buffer_Send, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/

      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive[iVar*nVertexR+iVertex] = Buffer_Send[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/

      delete [] Buffer_Send;

      /*--- Do the coordinate transformation ---*/

      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/

        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

        /*--- Copy transformed conserved variables back into buffer. ---*/

        for (iVar = 0; iVar < nVar; iVar++)
          vec_val[iPoint*nVar+iVar] += Buffer_Receive[iVertex*nVar+iVar];

      }

      /*--- Deallocate receive buffer ---*/

      delete [] Buffer_Receive;

    }

  }

}
#endif
