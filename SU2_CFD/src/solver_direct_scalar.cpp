/*!
 * \file solution_direct_scalar.cpp
 * \brief Main subroutines for solving scalar transport equations.
 * \author T. Economon
 * \version 6.1.0 "Falcon"
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

#include "../include/solver_structure.hpp"

CScalarSolver::CScalarSolver(void) : CSolver() {
  
  FlowPrimVar_i    = NULL;
  FlowPrimVar_j    = NULL;
  lowerlimit       = NULL;
  upperlimit       = NULL;
  nVertex          = NULL;
  nMarker          = 0;
  Inlet_ScalarVars = NULL;
  Scalar_Inf       = NULL;
  
}

CScalarSolver::CScalarSolver(CGeometry* geometry, CConfig *config) : CSolver() {

  FlowPrimVar_i    = NULL;
  FlowPrimVar_j    = NULL;
  lowerlimit       = NULL;
  upperlimit       = NULL;
  nMarker          = config->GetnMarker_All();
  Inlet_ScalarVars = NULL;
  Scalar_Inf       = NULL;
  
  /*--- Store the number of vertices on each marker for deallocation later ---*/
  
  nVertex = new unsigned long[nMarker];
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];
  
}

CScalarSolver::~CScalarSolver(void) {
  
  if (Inlet_ScalarVars != NULL) {
    for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
      if (Inlet_ScalarVars[iMarker] != NULL) {
        for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          delete [] Inlet_ScalarVars[iMarker][iVertex];
        }
        delete [] Inlet_ScalarVars[iMarker];
      }
    }
    delete [] Inlet_ScalarVars;
  }
  
  if (FlowPrimVar_i != NULL) delete [] FlowPrimVar_i;
  if (FlowPrimVar_j != NULL) delete [] FlowPrimVar_j;
  if (lowerlimit != NULL)    delete [] lowerlimit;
  if (upperlimit != NULL)    delete [] upperlimit;
  if (nVertex != NULL)       delete [] nVertex;
  if (Scalar_Inf != NULL)    delete [] Scalar_Inf;

}

void CScalarSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector, nBufferS_Scalar, nBufferR_Scalar;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL, *Buffer_Receive_muT = NULL, *Buffer_Send_muT = NULL;
  
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
      nBufferS_Scalar = nVertexS;             nBufferR_Scalar = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      Buffer_Receive_muT = new su2double [nBufferR_Scalar];
      Buffer_Send_muT = new su2double[nBufferS_Scalar];
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_muT[iVertex] = 0.0; //node[iPoint]->GetmuT();   // TDE
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      SU2_MPI::Sendrecv(Buffer_Send_muT, nBufferS_Scalar, MPI_DOUBLE, send_to, 1,
                        Buffer_Receive_muT, nBufferR_Scalar, MPI_DOUBLE, receive_from, 1, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        Buffer_Receive_muT[iVertex] = 0.0; //node[iPoint]->GetmuT();  // TDE
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      delete [] Buffer_Send_muT;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy conservative variables. ---*/
        //node[iPoint]->SetmuT(Buffer_Receive_muT[iVertex]);  // TDE
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_muT;
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CScalarSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
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
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Old(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
}

void CScalarSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
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
      nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
      Buffer_Send_Gradient = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
        
        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
        }
        
        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient(iVar, iDim, Gradient[iVar][iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;
      
    }
    
  }
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}

void CScalarSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
  su2double *Limiter = new su2double [nVar];
  
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
      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
      Buffer_Send_Limit = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetLimiter(iVar, Buffer_Receive_Limit[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
  }
  
  delete [] Limiter;
  
}

void CScalarSolver::Upwind_Residual(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CNumerics *numerics,
                                    CConfig *config,
                                    unsigned short iMesh) {
  
  su2double *Scalar_i, *Scalar_j, *Limiter_i = NULL, *Limiter_j = NULL;
  su2double *V_i, *V_j, **Gradient_i, **Gradient_j;
  su2double Project_Grad_i, Project_Grad_j;
  
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim, iVar;
  
  bool muscl         = (config->GetMUSCL_Scalar() && (iMesh == MESH_0));
  bool limiter       = (config->GetKind_SlopeLimit_Scalar() != NO_LIMITER);
  bool grid_movement = config->GetGrid_Movement();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive variables w/o reconstruction ---*/
    
    V_i = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
    V_j = solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive();
    numerics->SetPrimitive(V_i, V_j);
    
    /*--- Scalar variables w/o reconstruction ---*/
    
    Scalar_i = node[iPoint]->GetSolution();
    Scalar_j = node[jPoint]->GetSolution();
    numerics->SetScalarVar(Scalar_i, Scalar_j);
    
    /*--- Grid Movement ---*/
    
    if (grid_movement)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
    if (muscl) {
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      /*--- Mean flow primitive variables using gradient reconstruction and limiters ---*/
      
      Gradient_i = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
      Gradient_j = solver_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive();
      if (limiter) {
        Limiter_i = solver_container[FLOW_SOL]->node[iPoint]->GetLimiter_Primitive();
        Limiter_j = solver_container[FLOW_SOL]->node[jPoint]->GetLimiter_Primitive();
      }
      
      for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnPrimVarGrad(); iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        if (limiter) {
          FlowPrimVar_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          FlowPrimVar_i[iVar] = V_i[iVar] + Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }
      
      numerics->SetPrimitive(FlowPrimVar_i, FlowPrimVar_j);
      
      /*--- Scalar variables using gradient reconstruction and limiters ---*/
      
      Gradient_i = node[iPoint]->GetGradient();
      Gradient_j = node[jPoint]->GetGradient();
      if (limiter) {
        Limiter_i = node[iPoint]->GetLimiter();
        Limiter_j = node[jPoint]->GetLimiter();
      }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        if (limiter) {
          Solution_i[iVar] = Scalar_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          Solution_j[iVar] = Scalar_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          Solution_i[iVar] = Scalar_i[iVar] + Project_Grad_i;
          Solution_j[iVar] = Scalar_j[iVar] + Project_Grad_j;
        }
      }
      
      numerics->SetScalarVar(Solution_i, Solution_j);
      
    }
    
    /*--- Add and subtract residual ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    LinSysRes.AddBlock(iPoint, Residual);
    LinSysRes.SubtractBlock(jPoint, Residual);
    
    /*--- Implicit part ---*/
    
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    
  }
  
}

void CScalarSolver::Viscous_Residual(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics *numerics,
                                   CConfig *config,
                                     unsigned short iMesh,
                                     unsigned short iRKStep) {
  
  unsigned long iEdge, iPoint, jPoint;
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Points coordinates, and normal vector ---*/
    
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(),
                           solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive());
    
    /*--- Scalar variables w/o reconstruction. ---*/
    
    numerics->SetScalarVar(node[iPoint]->GetSolution(),
                           node[jPoint]->GetSolution());
    
    /*--- Scalar variable gradients. ---*/

    numerics->SetScalarVarGradient(node[iPoint]->GetGradient(),
                                   node[jPoint]->GetGradient());
    
    /*--- Mass diffusivity coefficients. ---*/
    
    numerics->SetDiffusionCoeff(node[iPoint]->GetDiffusivity(),
                                node[jPoint]->GetDiffusivity());
    
    /*--- Compute residuals and Jacobians ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    /*--- Add/subtract residual and update Jacobians ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    
  }
  
}

void CScalarSolver::ImplicitEuler_Iteration(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, *local_Res_TruncError, Vol;
  
  bool scalar_clipping = config->GetScalar_Clipping();
  su2double scalar_clipping_min = config->GetScalar_Clipping_Min();
  su2double scalar_clipping_max = config->GetScalar_Clipping_Max();
  
  bool adjoint = ( config->GetContinuous_Adjoint() ||
                  (config->GetDiscrete_Adjoint() && config->GetFrozen_Visc_Disc()));
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  if (incompressible) SetPreconditioner(geometry, solver_container, config);

  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Read the residual ---*/

    local_Res_TruncError = node[iPoint]->GetResTruncError();

    /*--- Read the volume ---*/
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
    Delta = Vol / (config->GetCFLRedCoeff_Scalar()*solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
    Jacobian.AddVal2Diag(iPoint, Delta);
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] +  + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]),
                 geometry->node[iPoint]->GetGlobalIndex(),
                 geometry->node[iPoint]->GetCoord());
    }
  }
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }
  
  /*--- Solve or smooth the linear system ---*/
  
  CSysSolve system;
  system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  if (!adjoint) {
    
    /*--- Update the scalar solution. ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (scalar_clipping) {
        node[iPoint]->AddClippedSolution(0, config->GetRelaxation_Factor_Scalar()*LinSysSol[iPoint],
                                         scalar_clipping_min, scalar_clipping_max);
      }
      else {
        node[iPoint]->AddSolution(0, config->GetRelaxation_Factor_Scalar()*LinSysSol[iPoint]);
      }
    }
    
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CScalarSolver::BC_Sym_Plane(CGeometry *geometry,
                                 CSolver **solver_container,
                                 CNumerics *conv_numerics,
                                 CNumerics *visc_numerics,
                                 CConfig *config,
                                 unsigned short val_marker) {
  
  /*--- Convective fluxes across symmetry plane are equal to zero. ---*/
  
}

void CScalarSolver::BC_Euler_Wall(CGeometry *geometry,
                                  CSolver **solver_container,
                                CNumerics *numerics,
                                  CConfig *config,
                                  unsigned short val_marker) {
  
  /*--- Convective fluxes across euler wall are equal to zero. ---*/
  
}

void CScalarSolver::BC_HeatFlux_Wall(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CNumerics *conv_numerics,
                                            CNumerics *visc_numerics,
                                            CConfig *config,
                                            unsigned short val_marker) {
  
  /*--- Convective fluxes across viscous walls are equal to zero. ---*/
  
}

void CScalarSolver::BC_Isothermal_Wall(CGeometry *geometry,
                                              CSolver **solver_container,
                                              CNumerics *conv_numerics,
                                              CNumerics *visc_numerics,
                                              CConfig *config,
                                              unsigned short val_marker) {
  
  /*--- Convective fluxes across viscous walls are equal to zero. ---*/
  
}

void CScalarSolver::BC_Far_Field(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CNumerics *conv_numerics,
                                        CNumerics *visc_numerics,
                                        CConfig *config,
                                        unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  su2double *Normal, *V_infty, *V_domain;
  
  bool grid_movement  = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the infinity ---*/
      
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Set scalar variable at the wall, and at infinity ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = Scalar_Inf[iVar];
      }
      conv_numerics->SetScalarVar(Solution_i, Solution_j);
      
      /*--- Set Normal (it is necessary to change the sign) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute residuals and Jacobians ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add residuals and Jacobians ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  delete [] Normal;
  
}

void CScalarSolver::BC_Inlet(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CNumerics *conv_numerics,
                                    CNumerics *visc_numerics,
                                    CConfig *config,
                                    unsigned short val_marker) {
  
  unsigned short iDim, iVar;
  unsigned long iVertex, iPoint;
  su2double *V_inlet, *V_domain, *Normal;
  
  Normal = new su2double[nDim];
  
  bool grid_movement  = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Allocate the value at the inlet ---*/
      
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Set the scalar variable states (prescribed for an inflow) ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = Inlet_ScalarVars[val_marker][iVertex][iVar];
      }
      conv_numerics->SetScalarVar(Solution_i, Solution_j);
      
      /*--- Set various other quantities in the conv_numerics class ---*/
      
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      //      /*--- Viscous contribution, commented out because serious convergence problems ---*/
      //
      //      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      //      visc_numerics->SetNormal(Normal);
      //
      //      /*--- Conservative variables w/o reconstruction ---*/
      //
      //      visc_numerics->SetPrimitive(V_domain, V_inlet);
      //
      //      /*--- Scalar variables w/o reconstruction, and its gradients ---*/
      //
      //      visc_numerics->SetScalarVar(Solution_i, Solution_j);
      //      visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      //
      //      /*--- Compute residual, and Jacobians ---*/
      //
      //      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      //
      //      /*--- Subtract residual, and update Jacobians ---*/
      //
      //      LinSysRes.SubtractBlock(iPoint, Residual);
      //      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}

void CScalarSolver::BC_Outlet(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics *conv_numerics,
                                     CNumerics *visc_numerics,
                                     CConfig *config,
                                     unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  su2double *V_outlet, *V_domain, *Normal;
  
  bool grid_movement  = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the outlet ---*/
      
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Set the scalar variables. Here we use a Neumann BC such
       that the scalar variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> ScalarVar_internal,
       Solution_j --> ScalarVar_outlet ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
      }
      conv_numerics->SetScalarVar(Solution_i, Solution_j);
      
      /*--- Set Normal (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      //      /*--- Viscous contribution, commented out because serious convergence problems ---*/
      //
      //      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      //      visc_numerics->SetNormal(Normal);
      //
      //      /*--- Conservative variables w/o reconstruction ---*/
      //
      //      visc_numerics->SetPrimitive(V_domain, V_outlet);
      //
      //      /*--- Scalar variables w/o reconstruction, and its gradients ---*/
      //
      //      visc_numerics->SetScalarVar(Solution_i, Solution_j);
      //      visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      //
      //      /*--- Compute residual, and Jacobians ---*/
      //
      //      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      //
      //      /*--- Subtract residual, and update Jacobians ---*/
      //
      //      LinSysRes.SubtractBlock(iPoint, Residual);
      //      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  
}

void CScalarSolver::SetInletAtVertex(su2double *val_inlet,
                                            unsigned short iMarker,
                                            unsigned long iVertex) {
  unsigned short iVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Inlet_ScalarVars[iMarker][iVertex][iVar] = val_inlet[Inlet_Position+iVar];
  }
  
}

su2double CScalarSolver::GetInletAtVertex(su2double *val_inlet,
                                          unsigned long val_inlet_point,
                                          unsigned short val_kind_marker,
                                          string val_marker,
                                          CGeometry *geometry,
                                          CConfig *config) {
  
  /*--- Local variables ---*/
  
  unsigned short iMarker, iDim, iVar;
  unsigned long iPoint, iVertex;
  su2double Area = 0.0;
  su2double Normal[3] = {0.0,0.0,0.0};
  
  /*--- Alias positions within inlet file for readability ---*/
  
  if (val_kind_marker == INLET_FLOW) {
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) &&
          (config->GetMarker_All_TagBound(iMarker) == val_marker)) {
        
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++){
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (iPoint == val_inlet_point) {
            
            /*-- Compute boundary face area for this vertex. ---*/
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            Area = sqrt(Area);
            
            /*--- Access and store the inlet variables for this vertex. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++)
              val_inlet[Inlet_Position+iVar] = Inlet_ScalarVars[iMarker][iVertex][iVar];
            
            /*--- Exit once we find the point. ---*/
            
            return Area;
            
          }
        }
      }
    }
    
  }
  
  /*--- If we don't find a match, then the child point is not on the
   current inlet boundary marker. Return zero area so this point does
   not contribute to the restriction operator and continue. ---*/
  
  return Area;
  
}

void CScalarSolver::SetUniformInlet(CConfig* config, unsigned short iMarker) {
  
  for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Inlet_ScalarVars[iMarker][iVertex][iVar] = Scalar_Inf[iVar];
  }
  
}



void CScalarSolver::SetResidual_DualTime(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CConfig *config,
                                       unsigned short iRKStep,
                                         unsigned short iMesh,
                                         unsigned short RunTime_EqSystem) {
  
  /*--- Local variables ---*/
  
  unsigned short iVar, jVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;
  
  su2double *U_time_nM1, *U_time_n, *U_time_nP1;
  su2double Volume_nM1, Volume_nP1, TimeStep;
  su2double Density_nM1, Density_n, Density_nP1;
  su2double *Normal = NULL, *GridVel_i = NULL, *GridVel_j = NULL, Residual_GCL;
  
  bool implicit      = (config->GetKind_TimeIntScheme_Scalar() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  /*--- Store the physical time step ---*/
  
  TimeStep = config->GetDelta_UnstTimeND();
  
  /*--- Compute the dual time-stepping source term for static meshes ---*/
  
  if (!grid_movement) {
    
    /*--- Loop over all nodes (excluding halos) ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/
      
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/
      
       /*--- Get the density to compute the conservative variables. ---*/
      
        if (incompressible){
          /*--- This is temporary and only valid for constant-density problems:
           density could also be temperature dependent, but as it is not a part
           of the solution vector it's neither stored for previous time steps
           nor updated with the solution at the end of each iteration. */
          Density_nM1 = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
          Density_n   = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
          Density_nP1 = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
        }
        else{
          Density_nM1 = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[0];
          Density_n   = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n()[0];
          Density_nP1 = solver_container[FLOW_SOL]->node[iPoint]->GetSolution()[0];
        }
        
        for (iVar = 0; iVar < nVar; iVar++) {
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Residual[iVar] = ( Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*Volume_nP1 / TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Residual[iVar] = ( 3.0*Density_nP1*U_time_nP1[iVar] - 4.0*Density_n*U_time_n[iVar]
                              +1.0*Density_nM1*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
        }
      
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
    
  } else {
    
    /*--- For unsteady flows on dynamic meshes (rigidly transforming or
     dynamically deforming), the Geometric Conservation Law (GCL) should be
     satisfied in conjunction with the ALE formulation of the governing
     equations. The GCL prevents accuracy issues caused by grid motion, i.e.
     a uniform free-stream should be preserved through a moving grid. First,
     we will loop over the edges and boundaries to compute the GCL component
     of the dual time source term that depends on grid velocities. ---*/
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      /*--- Get indices for nodes i & j plus the face normal ---*/
      
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      Normal = geometry->edge[iEdge]->GetNormal();
      
      /*--- Grid velocities stored at nodes i & j ---*/
      
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      
      /*--- Compute the GCL term by averaging the grid velocities at the
       edge mid-point and dotting with the face normal. ---*/
      
      Residual_GCL = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual_GCL += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      
      /*--- Compute the GCL component of the source term for node i ---*/
      
      U_time_n = node[iPoint]->GetSolution_time_n();
      
      /*--- Multiply by density at node i  ---*/
      
        if (incompressible) Density_n = solver_container[FLOW_SOL]->node[iPoint]->GetDensity(); // Temporary fix
        else Density_n = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n()[0];
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = Density_n*U_time_n[iVar]*Residual_GCL;

      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Compute the GCL component of the source term for node j ---*/
      
      U_time_n = node[jPoint]->GetSolution_time_n();
      
      /*--- Multiply by density at node j ---*/
      
        if (incompressible) Density_n = solver_container[FLOW_SOL]->node[jPoint]->GetDensity(); // Temporary fix
        else Density_n = solver_container[FLOW_SOL]->node[jPoint]->GetSolution_time_n()[0];
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = Density_n*U_time_n[iVar]*Residual_GCL;

      LinSysRes.SubtractBlock(jPoint, Residual);
      
    }
    
    /*---  Loop over the boundary edges ---*/
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          
          /*--- Get the index for node i plus the boundary face normal ---*/
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          
          /*--- Grid velocities stored at boundary node i ---*/
          
          GridVel_i = geometry->node[iPoint]->GetGridVel();
          
          /*--- Compute the GCL term by dotting the grid velocity with the face
           normal. The normal is negated to match the boundary convention. ---*/
          
          Residual_GCL = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Residual_GCL -= 0.5*(GridVel_i[iDim]+GridVel_i[iDim])*Normal[iDim];
          
          /*--- Compute the GCL component of the source term for node i ---*/
          
          U_time_n = node[iPoint]->GetSolution_time_n();
          
          /*--- Multiply by density at node i  ---*/
          
            if (incompressible) Density_n = solver_container[FLOW_SOL]->node[iPoint]->GetDensity(); // Temporary fix
            else Density_n = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n()[0];
            for (iVar = 0; iVar < nVar; iVar++)
              Residual[iVar] = Density_n*U_time_n[iVar]*Residual_GCL;

          LinSysRes.AddBlock(iPoint, Residual);
        }
    }
    
    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/
      
      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();

        if (incompressible){
          /*--- This is temporary and only valid for constant-density problems:
           density could also be temperature dependent, but as it is not a part
           of the solution vector it's neither stored for previous time steps
           nor updated with the solution at the end of each iteration. */
          Density_nM1 = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
          Density_n   = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
          Density_nP1 = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
        }
        else{
          Density_nM1 = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[0];
          Density_n   = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_time_n()[0];
          Density_nP1 = solver_container[FLOW_SOL]->node[iPoint]->GetSolution()[0];
        }
        
        for (iVar = 0; iVar < nVar; iVar++) {
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Residual[iVar] = (Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*(Volume_nP1/TimeStep);
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Residual[iVar] = (Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
            + (Density_nM1*U_time_nM1[iVar] - Density_n*U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
        }
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {  // TDE density in Jacobian
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1/TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (3.0*Volume_nP1)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }
  
}

void CScalarSolver::LoadRestart(CGeometry **geometry,
                                CSolver ***solver,
                                CConfig *config,
                                int val_iter,
                                bool val_update_geo) {
  
  /*--- Restart the solution from file information ---*/
  
  unsigned short iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine;
  su2double Area_Children, Area_Parent, *Solution_Fine;
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();
  
  string UnstExt, text_line;
  ifstream restart_file;
  string restart_filename = config->GetSolution_FlowFileName();
  
  bool turbulent = (config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS);
  bool turb_SST  = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool turb_SA   = ((turbulent) && (config->GetKind_Turb_Model() == SA));
  
  /*--- Modify file name for multizone problems ---*/
  if (nZone >1)
    restart_filename = config->GetMultizone_FileName(restart_filename, iZone);
  
  /*--- Modify file name for an unsteady restart ---*/
  
  if (dual_time|| time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);
  
  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/
  
  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }
  
  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
  
  /*--- Skip flow variables ---*/
  
  unsigned short skipVars = 0;
  
  if (nDim == 2) skipVars += 6;
  if (nDim == 3) skipVars += 8;
  
  if (turbulent) {
    if (turb_SA) skipVars += 1;
    else if (turb_SST) skipVars += 2;
  }
  
  /*--- Load data from the restart into correct containers. ---*/
  
  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {
    
    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/
    
    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);
    
    if (iPoint_Local > -1) {
      
      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/
      
      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      node[iPoint_Local]->SetSolution(Solution);
      iPoint_Global_Local++;
      
      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }
    
  }
  
  /*--- Detect a wrong solution file ---*/
  
  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }
  
#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1,
                     MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
    SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }
  
  /*--- MPI solution ---*/
  
  //TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart.
  solver[MESH_0][SCALAR_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][SCALAR_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  
  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  solver[MESH_0][SCALAR_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
  
  /*--- Interpolate the solution down to the coarse multigrid levels ---*/
  
  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][SCALAR_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][SCALAR_SOL]->node[iPoint]->SetSolution(Solution);
    }
    solver[iMesh][SCALAR_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    solver[iMesh][SCALAR_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
  }
  
  /*--- Delete the class memory that is used to load the restart. ---*/
  
  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;
  
}

CPassiveScalarSolver::CPassiveScalarSolver(void) : CScalarSolver() {
  
  Inlet_ScalarVars = NULL;
  
}

CPassiveScalarSolver::CPassiveScalarSolver(CGeometry *geometry,
                                           CConfig *config,
                                           unsigned short iMesh)
: CScalarSolver(geometry, config) {
  
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  su2double Density_Inf, Viscosity_Inf;
  
  bool turbulent = ((config->GetKind_Solver() == RANS) ||
  (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool turb_SST  = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool turb_SA   = ((turbulent) && (config->GetKind_Turb_Model() == SA));
  
  /*--- Dimension of the problem --> passive scalar will only ever
   have a single equation. Other child classes of CScalarSolver
   will have variable numbers of equations. ---*/
  
  nVar     = 1;
  nPrimVar = 1;
  
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  /*--- Define geometry constants in the solver structure ---*/
  
  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];
  
  /*--- Fluid model pointer initialization ---*/
  
  FluidModel = NULL;
  
  /*--- Define some auxiliar vector related with the residual ---*/
  
  Residual = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]  = 0.0;
  Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_i = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]  = 0.0;
  Residual_j = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]  = 0.0;
  Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliar vector related with the solution ---*/
  
  Solution = new su2double[nVar];
  Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
  
  /*--- Define some auxiliar vector related with the geometry ---*/
  
  Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
  
  /*--- Define some auxiliar vector related with the flow solution ---*/
  
  FlowPrimVar_i = new su2double [nDim+9]; FlowPrimVar_j = new su2double [nDim+9];
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  /*--- Initialization of the structure of the whole Jacobian ---*/
  
  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Passive Scalar). MG level: " << iMesh <<"." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  
  if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
      (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
    nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
    if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
  }
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Computation of gradients by least squares ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    
    /*--- c vector := transpose(WA)*(Wb) ---*/
    
    Cvector = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Initialize lower and upper limits---*/
  
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  
  lowerlimit[0] = -1.0e15;
  upperlimit[0] =  1.0e15;
  
  /*--- Read farfield conditions from config ---*/
  
  Density_Inf   = config->GetDensity_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();

  /*--- Set up fluid model for the diffusivity ---*/
  
  su2double Diffusivity_Ref = 1.0;
  
  su2double DiffusivityND = config->GetDiffusivity_Constant()/Diffusivity_Ref;
  
  config->SetDiffusivity_ConstantND(DiffusivityND);
  //config->SetDiffusivity_Ref(Diffusivity_Ref);
  
  FluidModel = new CFluidModel();
  FluidModel->SetMassDiffusivityModel(config);
  
  /*--- Scalar variable state at the far-field. ---*/
  
  Scalar_Inf = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Scalar_Inf[iVar] = config->GetScalar_Init();
  
  /*--- Initialize the solution to the far-field state everywhere. ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CScalarVariable(Scalar_Inf, nDim, nVar, config);
  
  /*--- MPI solution ---*/
  
  //TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart
  Set_MPI_Solution(geometry, config);
  Set_MPI_Solution(geometry, config);
  
  /*--- Initializate quantities for SlidingMesh Interface ---*/
  
  unsigned long iMarker;
  
  SlidingState       = new su2double*** [nMarker];
  SlidingStateNodes  = new int*         [nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++){
    
    SlidingState[iMarker]      = NULL;
    SlidingStateNodes[iMarker] = NULL;
    
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){
      
      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)];
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)];
      
      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++){
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1];
        
        SlidingStateNodes[iMarker][iPoint] = 0;
        for (iVar = 0; iVar < nPrimVar+1; iVar++)
          SlidingState[iMarker][iPoint][iVar] = NULL;
      }
      
    }
  }
  
  /*-- Allocation of inlets has to happen in derived classes
   (not CScalarSolver), due to arbitrary number of scalar variables.
   First, we also set the column index for any inlet profiles. ---*/
  
  Inlet_Position = nDim*2+2;
  if (turbulent) {
    if (turb_SA) Inlet_Position += 1;
    else if (turb_SST) Inlet_Position += 2;
  }
  
  Inlet_ScalarVars = new su2double**[nMarker];
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_ScalarVars[iMarker] = new su2double*[nVertex[iMarker]];
    for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
      Inlet_ScalarVars[iMarker][iVertex] = new su2double[nVar];
      for (unsigned short iVar = 0; iVar < nVar; iVar++)
        Inlet_ScalarVars[iMarker][iVertex][0] = Scalar_Inf[iVar];
    }
  }
  
}

CPassiveScalarSolver::~CPassiveScalarSolver(void) {
  
  unsigned long iMarker, iVertex;
  unsigned short iVar;
  
  if ( SlidingState != NULL ) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if ( SlidingState[iMarker] != NULL ) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          if ( SlidingState[iMarker][iVertex] != NULL ){
            for (iVar = 0; iVar < nPrimVar+1; iVar++)
              delete [] SlidingState[iMarker][iVertex][iVar];
            delete [] SlidingState[iMarker][iVertex];
          }
        delete [] SlidingState[iMarker];
      }
    }
    delete [] SlidingState;
  }
  
  if ( SlidingStateNodes != NULL ){
    for (iMarker = 0; iMarker < nMarker; iMarker++){
      if (SlidingStateNodes[iMarker] != NULL)
        delete [] SlidingStateNodes[iMarker];
    }
    delete [] SlidingStateNodes;
  }
  
  if (FluidModel != NULL) delete FluidModel;

}

void CPassiveScalarSolver::Preprocessing(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CConfig *config,
                                         unsigned short iMesh,
                                         unsigned short iRKStep,
                                         unsigned short RunTime_EqSystem,
                                         bool Output) {
  
  unsigned long ErrorCounter = 0;
  unsigned long ExtIter = config->GetExtIter();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool limiter_flow     = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) &&
                           (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_scalar   = ((config->GetKind_SlopeLimit_Scalar() != NO_LIMITER) &&
                           (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  /*--- Initialize the Jacobian matrices ---*/
  
  if (!disc_adjoint) Jacobian.SetValZero();
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  
  /*--- Upwind second-order reconstruction ---*/
  
  if ((limiter_scalar) && (iMesh == MESH_0)) SetSolution_Limiter(geometry, config);
  
  if ((limiter_flow) && (iMesh == MESH_0)) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);
  
}

void CPassiveScalarSolver::Postprocessing(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CConfig *config,
                                          unsigned short iMesh) { }

unsigned long CPassiveScalarSolver::SetPrimitive_Variables(CSolver **solver_container,
                                                           CConfig *config,
                                                           bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double Density, Temperature, lam_visc = 0.0, eddy_visc = 0.0, Cp;
  unsigned short iVar, turb_model = config->GetKind_Turb_Model();
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /*--- Retrieve the density, temperature, Cp, and laminar viscosity. ---*/

    Density     = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
    Cp          = solver_container[FLOW_SOL]->node[iPoint]->GetSpecificHeatCp();
    Temperature = solver_container[FLOW_SOL]->node[iPoint]->GetTemperature();
    lam_visc    = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    
    /*--- Retrieve the value of the kinetic energy (if needed) ---*/
    
    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->node[iPoint]->GetmuT();
    }
    
    /*--- Compute and store the mass diffusivity. ---*/

    FluidModel->SetDiffusivityState(Temperature, Density, lam_visc, eddy_visc, Cp);
    
    for (iVar = 0; iVar < nVar; iVar++)
      node[iPoint]->SetDiffusivity(FluidModel->GetMassDiffusivity(), iVar);
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  return ErrorCounter;
  
}

void CPassiveScalarSolver::SetInitialCondition(CGeometry **geometry,
                                               CSolver ***solver_container,
                                               CConfig *config,
                                               unsigned long ExtIter) {
  
  unsigned short iMesh, iVar;
  unsigned long iPoint;
  su2double *Coords;
  
  if (ExtIter == 0) {
    
    if (rank == MASTER_NODE)
      cout << "Setting the convection-diffusion initial condition." << endl;
    
    su2double *Init_Sol = new su2double[nVar];
    
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
        
        /*--- Set the initial condition for the conv-diff problem. ---*/
        Coords = geometry[iMesh]->node[iPoint]->GetCoord();
        if (Coords[1] <= 1.0)
          Init_Sol[0] = 1.0;
        else
          Init_Sol[0] = 0.0;
        
        solver_container[iMesh][SCALAR_SOL]->node[iPoint]->SetSolution(Init_Sol);
      }
      solver_container[iMesh][SCALAR_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    }
    delete [] Init_Sol;
    
  }
  
}


void CPassiveScalarSolver::SetPreconditioner(CGeometry *geometry, CSolver **solver_container,  CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  
  su2double  BetaInc2, Density, dRhodT, dRhodC, Temperature, Cp, Delta;
  
  bool variable_density = (config->GetKind_DensityModel() == VARIABLE);
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Access the primitive variables at this node. ---*/
    
    Density     = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
    BetaInc2    = solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2();
    Cp          = solver_container[FLOW_SOL]->node[iPoint]->GetSpecificHeatCp();
    Temperature = solver_container[FLOW_SOL]->node[iPoint]->GetTemperature();
    
    unsigned short nVar_Flow = solver_container[FLOW_SOL]->GetnVar();
    
    su2double SolP = solver_container[FLOW_SOL]->LinSysSol[iPoint*nVar_Flow+0];
    su2double SolT = solver_container[FLOW_SOL]->LinSysSol[iPoint*nVar_Flow+nDim+1];
    
    /*--- We need the derivative of the equation of state to build the
     preconditioning matrix. For now, the only option is the ideal gas
     law, but in the future, dRhodT should be in the fluid model. ---*/
    
    if (variable_density) {
      dRhodT = -Density/Temperature;
    } else {
      dRhodT = 0.0;
    }
    
    /*--- Passive scalars have no impact on the density. ---*/
    
    dRhodC = 0.0;
    
    /*--- Modify matrix diagonal with term including volume and time step. ---*/
    
    su2double Vol = geometry->node[iPoint]->GetVolume();
    Delta = Vol / (config->GetCFLRedCoeff_Scalar()*
                   solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
    
    /*--- Calculating the inverse of the preconditioning matrix
     that multiplies the time derivative during time integration. ---*/
    
    if (implicit) {
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        total_index = iPoint*nVar+iVar;
        
        su2double c = node[iPoint]->GetSolution(0);
        
        /*--- Compute the lag terms for the decoupled linear system from
         the mean flow equations and add to the residual for the scalar.
         In short, we are effectively making these terms explicit. ---*/
        
        su2double artcompc1 = SolP * c/(Density*BetaInc2);
        su2double artcompc2 = SolT * dRhodT * c/(Density);
        
        LinSysRes[total_index] += artcompc1 + artcompc2;
        
        /*--- Add the extra Jacobian term to the scalar system. ---*/
        
        su2double Jaccomp = c * dRhodC + Density; //This is Gamma
        su2double JacTerm = Jaccomp*Delta;
        
        Jacobian.AddVal2Diag(iPoint, JacTerm);
        
      }
      
    }
    
  }
  
}

void CPassiveScalarSolver::Source_Residual(CGeometry *geometry,
                                           CSolver **solver_container,
                                           CNumerics *numerics,
                                           CNumerics *second_numerics,
                                    CConfig *config,
                                           unsigned short iMesh) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool axisymmetric   = config->GetAxisymmetric();
  bool viscous        = config->GetViscous();
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(), NULL);
    
    /*--- Gradient of the primitive and conservative variables ---*/
    
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
    
    /*--- Scalar variables w/o reconstruction, and its gradient ---*/
    
    numerics->SetScalarVar(node[iPoint]->GetSolution(), NULL);
    numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), NULL);
    
    /*--- Set volume ---*/
    
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    
    /*--- Compute the source term ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);
    
    /*--- Subtract residual and the Jacobian ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    
  }
  
  /*--- Axisymmetry source term for the scalar equation. ---*/
  
  if (axisymmetric) {
    
    /*--- Zero out Jacobian structure ---*/
    
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar ++)
        for (unsigned short jVar = 0; jVar < nVar; jVar ++)
          Jacobian_i[iVar][jVar] = 0.0;
    }
    
    /*--- loop over points ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Primitive variables w/o reconstruction ---*/
      
      second_numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(), NULL);
      
      /*--- Scalar variables w/o reconstruction ---*/
      
      second_numerics->SetScalarVar(node[iPoint]->GetSolution(), NULL);
      
      /*--- Mass diffusivity coefficients. ---*/
      
      second_numerics->SetDiffusionCoeff(node[iPoint]->GetDiffusivity(),
                                  NULL);
      
      /*--- Set control volume ---*/
      
      second_numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Set y coordinate ---*/
      
      second_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                         NULL);
      
      /*--- If viscous, we need gradients for extra terms. ---*/
      
      if (viscous) {
        
        /*--- Gradient of the scalar variables ---*/
        
        second_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), NULL);
        
      }
      
      /*--- Compute Source term Residual ---*/
      
      second_numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add Residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
      
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
}

CCombustionScalarSolver::CCombustionScalarSolver(void) : CScalarSolver() {
  
  Inlet_ScalarVars = NULL;
  
}

CCombustionScalarSolver::CCombustionScalarSolver(CGeometry *geometry,
                                           CConfig *config,
                                           unsigned short iMesh)
: CScalarSolver(geometry, config) {
  
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  su2double Density_Inf, Viscosity_Inf;
  
  bool turbulent = ((config->GetKind_Solver() == RANS) ||
                    (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool turb_SST  = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool turb_SA   = ((turbulent) && (config->GetKind_Turb_Model() == SA));
  
  /*--- Dimension of the problem --> passive scalar will only ever
   have a single equation. Other child classes of CScalarSolver
   will have variable numbers of equations. ---*/
  
  nVar     = 1;
  nPrimVar = 1;
  
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  /*--- Define geometry constants in the solver structure ---*/
  
  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];
  
  /*--- Fluid model pointer initialization ---*/
  
  FluidModel = NULL;
  
  /*--- Define some auxiliar vector related with the residual ---*/
  
  Residual = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]  = 0.0;
  Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_i = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]  = 0.0;
  Residual_j = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]  = 0.0;
  Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliar vector related with the solution ---*/
  
  Solution = new su2double[nVar];
  Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
  
  /*--- Define some auxiliar vector related with the geometry ---*/
  
  Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
  
  /*--- Define some auxiliar vector related with the flow solution ---*/
  
  FlowPrimVar_i = new su2double [nDim+9]; FlowPrimVar_j = new su2double [nDim+9];
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  /*--- Initialization of the structure of the whole Jacobian ---*/
  
  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Combustion Scalar). MG level: " << iMesh <<"." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  
  if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
      (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
    nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
    if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
  }
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Computation of gradients by least squares ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    
    /*--- c vector := transpose(WA)*(Wb) ---*/
    
    Cvector = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Initialize lower and upper limits---*/
  
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  
  lowerlimit[0] = -1.0e15;
  upperlimit[0] =  1.0e15;
  
  /*--- Read farfield conditions from config ---*/
  
  Density_Inf   = config->GetDensity_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  
  /*--- Set up fluid model for the diffusivity ---*/
  
  su2double Diffusivity_Ref = 1.0;
  
  su2double DiffusivityND = config->GetDiffusivity_Constant()/Diffusivity_Ref;
  
  config->SetDiffusivity_ConstantND(DiffusivityND);
  //config->SetDiffusivity_Ref(Diffusivity_Ref);
  
  FluidModel = new CFluidModel();
  FluidModel->SetMassDiffusivityModel(config);
  
  /*--- Scalar variable state at the far-field. ---*/
  
  Scalar_Inf = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Scalar_Inf[iVar] = config->GetScalar_Init();
  
  /*--- Initialize the solution to the far-field state everywhere. ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CScalarVariable(Scalar_Inf, nDim, nVar, config);
  
  /*--- MPI solution ---*/
  
  //TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart
  Set_MPI_Solution(geometry, config);
  Set_MPI_Solution(geometry, config);
  
  /*--- Initializate quantities for SlidingMesh Interface ---*/
  
  unsigned long iMarker;
  
  SlidingState       = new su2double*** [nMarker];
  SlidingStateNodes  = new int*         [nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++){
    
    SlidingState[iMarker]      = NULL;
    SlidingStateNodes[iMarker] = NULL;
    
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){
      
      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)];
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)];
      
      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++){
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1];
        
        SlidingStateNodes[iMarker][iPoint] = 0;
        for (iVar = 0; iVar < nPrimVar+1; iVar++)
          SlidingState[iMarker][iPoint][iVar] = NULL;
      }
      
    }
  }
  
  /*-- Allocation of inlets has to happen in derived classes
   (not CScalarSolver), due to arbitrary number of scalar variables.
   First, we also set the column index for any inlet profiles. ---*/
  
  Inlet_Position = nDim*2+2;
  if (turbulent) {
    if (turb_SA) Inlet_Position += 1;
    else if (turb_SST) Inlet_Position += 2;
  }
  
  Inlet_ScalarVars = new su2double**[nMarker];
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_ScalarVars[iMarker] = new su2double*[nVertex[iMarker]];
    for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
      Inlet_ScalarVars[iMarker][iVertex] = new su2double[nVar];
      for (unsigned short iVar = 0; iVar < nVar; iVar++)
        Inlet_ScalarVars[iMarker][iVertex][0] = Scalar_Inf[iVar];
    }
  }
  
}

CCombustionScalarSolver::~CCombustionScalarSolver(void) {
  
  unsigned long iMarker, iVertex;
  unsigned short iVar;
  
  if ( SlidingState != NULL ) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if ( SlidingState[iMarker] != NULL ) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          if ( SlidingState[iMarker][iVertex] != NULL ){
            for (iVar = 0; iVar < nPrimVar+1; iVar++)
              delete [] SlidingState[iMarker][iVertex][iVar];
            delete [] SlidingState[iMarker][iVertex];
          }
        delete [] SlidingState[iMarker];
      }
    }
    delete [] SlidingState;
  }
  
  if ( SlidingStateNodes != NULL ){
    for (iMarker = 0; iMarker < nMarker; iMarker++){
      if (SlidingStateNodes[iMarker] != NULL)
        delete [] SlidingStateNodes[iMarker];
    }
    delete [] SlidingStateNodes;
  }
  
}

void CCombustionScalarSolver::Preprocessing(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CConfig *config,
                                         unsigned short iMesh,
                                         unsigned short iRKStep,
                                         unsigned short RunTime_EqSystem,
                                         bool Output) {
  
  unsigned long ErrorCounter = 0;
  unsigned long ExtIter = config->GetExtIter();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool limiter_flow     = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) &&
                           (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_scalar   = ((config->GetKind_SlopeLimit_Scalar() != NO_LIMITER) &&
                           (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  /*--- Initialize the Jacobian matrices ---*/
  
  if (!disc_adjoint) Jacobian.SetValZero();
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  
  /*--- Upwind second order reconstruction ---*/
  
  if (limiter_scalar) SetSolution_Limiter(geometry, config);
  
  if (limiter_flow) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);
  
}

void CCombustionScalarSolver::Postprocessing(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CConfig *config,
                                          unsigned short iMesh) { }

unsigned long CCombustionScalarSolver::SetPrimitive_Variables(CSolver **solver_container,
                                                           CConfig *config,
                                                           bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double Density, Temperature, lam_visc = 0.0, eddy_visc = 0.0, Cp;
  unsigned short iVar, turb_model = config->GetKind_Turb_Model();
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /*--- Retrieve the density, temperature, Cp, and laminar viscosity. ---*/
    
    Density     = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
    Cp          = solver_container[FLOW_SOL]->node[iPoint]->GetSpecificHeatCp();
    Temperature = solver_container[FLOW_SOL]->node[iPoint]->GetTemperature();
    lam_visc    = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    
    /*--- Retrieve the value of the kinetic energy (if needed) ---*/
    
    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->node[iPoint]->GetmuT();
    }
    
    /*--- Compute and store the mass diffusivity. ---*/
    
    FluidModel->SetDiffusivityState(Temperature, Density, lam_visc, eddy_visc, Cp);
    
    for (iVar = 0; iVar < nVar; iVar++)
      node[iPoint]->SetDiffusivity(FluidModel->GetMassDiffusivity(), iVar);
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  return ErrorCounter;
  
}

void CCombustionScalarSolver::SetInitialCondition(CGeometry **geometry,
                                               CSolver ***solver_container,
                                               CConfig *config,
                                               unsigned long ExtIter) {
  
  if (ExtIter == 0) {
    
    if (rank == MASTER_NODE)
      cout << "Setting the combustion initial condition." << endl;
    
  }
  
}


void CCombustionScalarSolver::SetPreconditioner(CGeometry *geometry, CSolver **solver_container,  CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  
  su2double  BetaInc2, Density, dRhodT, dRhodC, Temperature, Cp, Delta;
  
  bool variable_density = (config->GetKind_DensityModel() == VARIABLE);
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Access the primitive variables at this node. ---*/
    
    Density     = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
    BetaInc2    = solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2();
    Cp          = solver_container[FLOW_SOL]->node[iPoint]->GetSpecificHeatCp();
    Temperature = solver_container[FLOW_SOL]->node[iPoint]->GetTemperature();
    
    unsigned short nVar_Flow = solver_container[FLOW_SOL]->GetnVar();
    
    su2double SolP = solver_container[FLOW_SOL]->LinSysSol[iPoint*nVar_Flow+0];
    su2double SolT = solver_container[FLOW_SOL]->LinSysSol[iPoint*nVar_Flow+nDim+1];
    
    /*--- We need the derivative of the equation of state to build the
     preconditioning matrix. For now, the only option is the ideal gas
     law, but in the future, dRhodT should be in the fluid model. ---*/
    
    if (variable_density) {
      dRhodT = -Density/Temperature;
    } else {
      dRhodT = 0.0;
    }
    
    /*--- Passive scalars have no impact on the density. ---*/
    
    dRhodC = 0.0;
    
    /*--- Modify matrix diagonal with term including volume and time step. ---*/
    
    su2double Vol = geometry->node[iPoint]->GetVolume();
    Delta = Vol / (config->GetCFLRedCoeff_Scalar()*
                   solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
    
    /*--- Calculating the inverse of the preconditioning matrix
     that multiplies the time derivative during time integration. ---*/
    
    if (implicit) {
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        total_index = iPoint*nVar+iVar;
        
        su2double c = node[iPoint]->GetSolution(0);
        
        /*--- Compute the lag terms for the decoupled linear system from
         the mean flow equations and add to the residual for the scalar.
         In short, we are effectively making these terms explicit. ---*/
        
        su2double artcompc1 = SolP * c/(Density*BetaInc2);
        su2double artcompc2 = SolT * dRhodT * c/(Density);
        
        LinSysRes[total_index] += artcompc1 + artcompc2;
        
        /*--- Add the extra Jacobian term to the scalar system. ---*/
        
        su2double Jaccomp = c * dRhodC + Density; //This is Gamma
        su2double JacTerm = Jaccomp*Delta;
        
        Jacobian.AddVal2Diag(iPoint, JacTerm);
        
      }
      
    }
    
  }
  
}

void CCombustionScalarSolver::Source_Residual(CGeometry *geometry,
                                           CSolver **solver_container,
                                           CNumerics *numerics,
                                           CNumerics *second_numerics,
                                           CConfig *config,
                                           unsigned short iMesh) {
  
  unsigned short iVar, iDim;
  unsigned long iPoint;
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool axisymmetric   = config->GetAxisymmetric();
  bool viscous        = config->GetViscous();
  
  /*--- Combustion source term implementation. (active by default)  ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    su2double rho_u = 1.0;    //config->GetTransportedScalar_UnburntDensity();
    //su2double alpha_u = 0.01; //config->GetTransportedScalar_UnburntThermalDiffusion();
    su2double Sl = 0.5;       //config->GetLaminarFlamespeed();
    
    su2double **ScalarVar_Grad, GradNorm2 = 0.0;
    
    /*--- Get volume of the dual cell. ---*/
    
    su2double Volume = geometry->node[iPoint]->GetVolume();
    
    /*--- Get the gradient and compute its magnitude squared. ---*/
    
    ScalarVar_Grad = node[iPoint]->GetGradient();
    GradNorm2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      GradNorm2 += ScalarVar_Grad[0][iDim]*ScalarVar_Grad[0][iDim];
    
    /*--- Compute the production term. ---*/
    
    Residual[0] = rho_u*Sl*sqrt(GradNorm2)*Volume;
    
    /*--- Implicit part for production term (to do). ---*/

    Jacobian_i[0][0] = 0.0;
    
    /*--- Add Residual ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    
    /*--- Implicit part ---*/
    
    if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    
  }
  
  /*--- Axisymmetry source term for the scalar equation. ---*/
  
  if (axisymmetric) {
    
    /*--- Zero out Jacobian structure ---*/
    
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar ++)
        for (unsigned short jVar = 0; jVar < nVar; jVar ++)
          Jacobian_i[iVar][jVar] = 0.0;
    }
    
    /*--- loop over points ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Primitive variables w/o reconstruction ---*/
      
      second_numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(), NULL);
      
      /*--- Scalar variables w/o reconstruction ---*/
      
      second_numerics->SetScalarVar(node[iPoint]->GetSolution(), NULL);
      
      /*--- Mass diffusivity coefficients. ---*/
      
      second_numerics->SetDiffusionCoeff(node[iPoint]->GetDiffusivity(),
                                         NULL);
      
      /*--- Set control volume ---*/
      
      second_numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Set y coordinate ---*/
      
      second_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                NULL);
      
      /*--- If viscous, we need gradients for extra terms. ---*/
      
      if (viscous) {
        
        /*--- Gradient of the scalar variables ---*/
        
        second_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), NULL);
        
      }
      
      /*--- Compute Source term Residual ---*/
      
      second_numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add Residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
      
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
}
