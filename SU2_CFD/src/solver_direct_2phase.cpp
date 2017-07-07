/*!
 * \file solution_direct_2phase.cpp
 * \brief Main subrotuines for solving direct problems
 * \author F. Palacios, A. Bueno
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

#include "../include/solver_structure.hpp"

C2phaseSolver::C2phaseSolver(void) : CSolver() {
  
  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;
  lowerlimit    = NULL;
  upperlimit    = NULL;

}

C2phaseSolver::C2phaseSolver(CConfig *config) : CSolver() {

 /*
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  */

  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;
  lowerlimit    = NULL;
  upperlimit    = NULL;

}

C2phaseSolver::~C2phaseSolver(void) {
  
  if (FlowPrimVar_i != NULL) delete [] FlowPrimVar_i;
  if (FlowPrimVar_j != NULL) delete [] FlowPrimVar_j;
  if (lowerlimit != NULL) delete [] lowerlimit;
  if (upperlimit != NULL) delete [] upperlimit;
  
//  if (Primitive_Liquid != NULL) delete [] Primitive_Liquid;

}

void C2phaseSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector, nBufferS_Scalar, nBufferR_Scalar;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL, *Buffer_Receive_muT = NULL, *Buffer_Send_muT = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
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
      
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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
        
        /*--- Copy conservative variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
}

void C2phaseSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
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

void C2phaseSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
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

void C2phaseSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
  su2double *Limiter = new su2double [nVar];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
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



void C2phaseSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {
  
	  su2double *Two_phase_i, *Two_phase_j, *Limiter_i = NULL, *Limiter_j = NULL, *V_i, *V_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
	  su2double *Liquid_vec;
	  unsigned long iEdge, iPoint, jPoint;
	  unsigned short iDim, iVar;

	  bool second_order  = ((config->GetSpatialOrder() == SECOND_ORDER) || (config->GetSpatialOrder() == SECOND_ORDER_LIMITER));
	  bool limiter       = (config->GetSpatialOrder() == SECOND_ORDER_LIMITER);
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

	    /*--- 2phase variables w/o reconstruction ---*/

	    Two_phase_i = node[iPoint]->GetSolution();
	    Two_phase_j = node[jPoint]->GetSolution();
	    numerics->Set2phaseVar(Two_phase_i, Two_phase_j);

	    Liquid_vec = node[iPoint]->SetLiquidPrim(V_i, Two_phase_i, numerics->Primitive_Liquid[6], FluidModel, config);
	    numerics->SetPrimitive_Liquid(Liquid_vec);

	    /*--- Grid Movement ---*/

	    if (grid_movement)
	      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());

	    if (second_order) {

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

	      /*--- Turbulent variables using gradient reconstruction and limiters ---*/

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
	          Solution_i[iVar] = Two_phase_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
	          Solution_j[iVar] = Two_phase_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
	        }
	        else {
	          Solution_i[iVar] = Two_phase_i[iVar] + Project_Grad_i;
	          Solution_j[iVar] = Two_phase_j[iVar] + Project_Grad_j;
	        }
	      }

	      numerics->Set2phaseVar(Solution_i, Solution_j);

	//      numerics->SetPrimitive_Liquid(node[iPoint]->SetLiquidPrim(FlowPrimVar_i, Solution_i, config));

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




void C2phaseSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  

	  /*--- Convective fluxes across symmetry plane are equal to zero. ---*/

}

void C2phaseSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                CNumerics *numerics, CConfig *config, unsigned short val_marker) {

	 /*--- Convective fluxes across wall are equal to zero. ---*/

  
}

void C2phaseSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res, *Sol;
  unsigned short iVar;
  unsigned long iPoint;

  bool autoreset = config->Get_AutoReset_NegativeSol();

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() / Vol;

    local_Res_TruncError = node[iPoint]->GetResTruncError();
    local_Residual = LinSysRes.GetBlock(iPoint);

    for (iVar = 0; iVar < nVar; iVar++) {
      Res = local_Residual[iVar] + local_Res_TruncError[iVar];
      node[iPoint]->AddSolution(iVar, -Res*Delta);
      AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());

    }

    if (autoreset) {
		Sol = node[iPoint]->GetSolution();

		for (iVar = 0; iVar < nVar; iVar++) {
		   Sol[iVar] = max(0.0, Sol[iVar]);
		}

		node[iPoint]->SetSolution(Sol);
    }
  }

  /*--- MPI solution ---*/

  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);


}

void C2phaseSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, Vol, density_old = 0.0, density = 0.0;
  su2double *Sol;
  

  bool adjoint = config->GetContinuous_Adjoint();
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool autoreset = config->Get_AutoReset_NegativeSol();
  
  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Read the volume ---*/
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    Delta = Vol / (config->GetCFLRedCoeff_2phase()*solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
    Jacobian.AddVal2Diag(iPoint, Delta);
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - LinSysRes[total_index];
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
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
    
       for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

        	for (iVar = 0; iVar < nVar; iVar++) {
        		node[iPoint]->AddSolution(iVar, config->GetRelaxation_Factor_2phase()*LinSysSol[iPoint*nVar+iVar]);
        	}

        	if (autoreset) {
				Sol = node[iPoint]->GetSolution();

				for (iVar = 0; iVar < nVar; iVar++) {
					Sol[iVar] = max(Sol[iVar], 0.0);
				}

				node[iPoint]-> SetSolution(Sol);
        	}

       }
  }
  



  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void C2phaseSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                       unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
  
	  /*--- Local variables ---*/

	  unsigned short iVar, jVar, iMarker, iDim;
	  unsigned long iPoint, jPoint, iEdge, iVertex;

	  su2double *U_time_nM1, *U_time_n, *U_time_nP1;
	  su2double Volume_nM1, Volume_nP1, TimeStep;
	  su2double Density_nM1, Density_n, Density_nP1;
	  su2double *Normal = NULL, *GridVel_i = NULL, *GridVel_j = NULL, Residual_GCL;

	  bool implicit      = (config->GetKind_TimeIntScheme_2phase() == EULER_IMPLICIT);
	  bool grid_movement = config->GetGrid_Movement();

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

	        for (iVar = 0; iVar < nVar; iVar++) {
	          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
	            Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
	          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
	            Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
	                              +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
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

          for (iVar = 0; iVar < nVar; iVar++)
	          Residual[iVar] = U_time_n[iVar]*Residual_GCL;

	      LinSysRes.AddBlock(iPoint, Residual);

	      /*--- Compute the GCL component of the source term for node j ---*/

	      U_time_n = node[jPoint]->GetSolution_time_n();

          for (iVar = 0; iVar < nVar; iVar++)
	          Residual[iVar] = U_time_n[iVar]*Residual_GCL;

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

	          for (iVar = 0; iVar < nVar; iVar++)
	            Residual[iVar] = U_time_n[iVar]*Residual_GCL;

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

	      /*--- Compute the dual time-stepping source residual. Due to the
	       introduction of the GCL term above, the remainder of the source residual
	       due to the time discretization has a new form.---*/

	        for (iVar = 0; iVar < nVar; iVar++) {
	          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
	            Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
	          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
	            Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
	            + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
	        }

	      /*--- Store the residual and compute the Jacobian contribution due
	       to the dual time source term. ---*/

	      LinSysRes.AddBlock(iPoint, Residual);
	      if (implicit) {
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


void C2phaseSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {


	  /*--- Restart the solution from file information ---*/

	  unsigned short iVar, iMesh;
	  unsigned long iPoint, index, iChildren, Point_Fine;
	  su2double Area_Children, Area_Parent, *Solution_Fine, *Sol;
	  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
	  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
	  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
	                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
	  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
	  unsigned short iZone = config->GetiZone();
	  unsigned short nZone = config->GetnZone();

	  string UnstExt, text_line;
	  ifstream restart_file;
	  string restart_filename = config->GetSolution_FlowFileName();

	  int rank = MASTER_NODE;
	#ifdef HAVE_MPI
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

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

	  if (compressible) {
	    if (nDim == 2) skipVars += 6;
	    if (nDim == 3) skipVars += 8;
	  }
	  if (incompressible) {
	    if (nDim == 2) skipVars += 5;
	    if (nDim == 3) skipVars += 7;
	  }

	  if (config->GetKind_Turb_Model() != NONE) {
		  if (config->GetKind_Turb_Model() == SST)
			  skipVars += 2;
		  else
			  skipVars += 1;
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

	      index = counter*Restart_Vars[0] + skipVars;
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
	  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
	#endif
	  if (rbuf_NotMatching != 0) {
	    if (rank == MASTER_NODE) {
	      cout << endl << "The solution file " << restart_filename.data() << " doesn't match with the mesh file!" << endl;
	      cout << "It could be empty lines at the end of the file." << endl << endl;
	    }
	#ifndef HAVE_MPI
	    exit(EXIT_FAILURE);
	#else
	    MPI_Barrier(MPI_COMM_WORLD);
	    MPI_Abort(MPI_COMM_WORLD,1);
	    MPI_Finalize();
	#endif
	  }


	  solver[MESH_0][TWO_PHASE_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
	  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
	  solver[MESH_0][TWO_PHASE_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);

	  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

	  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
	    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
	      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
	      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
	      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
	        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
	        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
	        Solution_Fine = solver[iMesh-1][TWO_PHASE_SOL]->node[Point_Fine]->GetSolution();
	        for (iVar = 0; iVar < nVar; iVar++) {
	          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
	        }
	      }
	      solver[iMesh][TWO_PHASE_SOL]->node[iPoint]->SetSolution(Solution);
	    }
	    solver[iMesh][TWO_PHASE_SOL]->Set_MPI_Solution(geometry[iMesh], config);
	    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
	    solver[iMesh][TWO_PHASE_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);

	    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
	    	Sol = node[iPoint]->SetLiquidPrim(solver[iMesh][FLOW_SOL]->node[iPoint]->GetPrimitive(),
	    			           node[iPoint]->GetSolution(), 1e-10,
	    			           solver[iMesh][FLOW_SOL]->GetFluidModel(), config);
	    }
	  }


	  /*--- Delete the class memory that is used to load the restart. ---*/

	  if (Restart_Vars != NULL) delete [] Restart_Vars;
	  if (Restart_Data != NULL) delete [] Restart_Data;
	  Restart_Vars = NULL; Restart_Data = NULL;



	}


C2phase_HillSolver::C2phase_HillSolver(void) : C2phaseSolver() {
  
}

C2phase_HillSolver::C2phase_HillSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh, CFluidModel* FluidModel) : C2phaseSolver() {
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  ifstream restart_file;
  string text_line;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Dimension of the problem --> dependent on the 2phase model. ---*/
  nVar = 4;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  /*--- Define geometry constants in the solver structure ---*/
  
  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];
  
  /*--- Single grid simulation ---*/
  
  if (iMesh == MESH_0) {
    
    /*--- Define some auxiliary vector related with the residual ---*/
    
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
    
    /*--- Define some auxiliary vector related with the solution ---*/
    Solution = new su2double[nVar]; Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
    for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]  = 0.0;

    /*--- Define some auxiliary vectors related to the liquid phase ---*/
    /*Primitive_Liquid = new su2double[9];*/
    //for (iVar = 0; iVar < 10; iVar++) Primitive_Liquid[iVar]  = 0.0;


    /*--- Define some auxiliary vector related with the geometry ---*/
    Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
    
    /*--- Define some auxiliary vector related with the flow solution ---*/
    FlowPrimVar_i = new su2double [nDim+9]; FlowPrimVar_j = new su2double [nDim+9];
    
    /*--- Jacobians and vector structures for implicit computations ---*/
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }
    
    /*--- Initialization of the structure of the whole Jacobian ---*/
    
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (MOM model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  }
  
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
  
  lowerlimit[0] = 1.0e-50;
  upperlimit[0] = 1.0e40;

  lowerlimit[1] = 1.0e-50;
  upperlimit[1] = 1.0e20;

  lowerlimit[2] = 1.0e-60;
  upperlimit[2] = 1.0e10;
  
  lowerlimit[3] = 1.0e-70;
  upperlimit[3] = 1.0e10;

  
  /*--- Far-field flow state quantities and initialization. ---*/
  su2double RInf, NInf, DInf;

  DInf    = 1;
  RInf    = 0;
  NInf    = 0;

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint] = new C2phase_HillVariable(RInf, NInf, DInf, nDim, nVar, config);
  }

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
}

C2phase_HillSolver::~C2phase_HillSolver(void) {
  
}



void C2phase_HillSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint;

  unsigned long ExtIter      = config->GetExtIter();
  bool limiter_flow          = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    
    /*--- Initialize the residual vector ---*/
    LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  /*--- Initialize the Jacobian matrices ---*/
  
  Jacobian.SetValZero();

  /*--- Upwind second order reconstruction ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

  if (config->GetSpatialOrder_2phase() == SECOND_ORDER_LIMITER) SetSolution_Limiter(geometry, config);
  
  if (limiter_flow) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);

}

void C2phase_HillSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {

  su2double *Two_phase_i, *Two_phase_j, *Limiter_i = NULL, *Limiter_j = NULL, *V_i, *V_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
  su2double V_Left, V_Right, Sol_Left, Sol_Right, Delta, Delta_Left, Delta_Right;
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim, iVar, jVar;

  bool second_order  = ((config->GetSpatialOrder_2phase() == SECOND_ORDER) || (config->GetSpatialOrder_2phase() == SECOND_ORDER_LIMITER));
  bool limiter       = (config->GetSpatialOrder_2phase() == SECOND_ORDER_LIMITER);
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

    /*--- 2phase variables w/o reconstruction ---*/

    Two_phase_i = node[iPoint]->GetSolution();
    Two_phase_j = node[jPoint]->GetSolution();

    numerics->Set2phaseVar(Two_phase_i, Two_phase_j);

    /*--- Grid Movement ---*/

    if (grid_movement)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());

    if (second_order) {

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

        } else {
          FlowPrimVar_i[iVar] = V_i[iVar] + Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }

      numerics->SetPrimitive(FlowPrimVar_i, FlowPrimVar_j);

      /*--- Turbulent variables using gradient reconstruction and limiters ---*/

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

        	Sol_Left  = max(Two_phase_i[iVar] - Project_Grad_i, 0.0);
        	Sol_Right = max(Two_phase_j[iVar] - Project_Grad_j, 0.0);

        	if (config->GetKind_SlopeLimit_2phase() == MINMOD) {

        		Delta_Left = Two_phase_i[iVar]-Sol_Left; Delta_Right = (Two_phase_j[iVar]-Two_phase_i[iVar])/2;

        		if (Delta_Left > 0 && Delta_Right > 0) {
        		    Delta    = 1.0;
        		} else if (Delta_Left < 0 && Delta_Right < 0) {
        			Delta = -1.0;
        		} else {
        			Delta = 0;
        		}
        		Delta = Delta *min(fabs(Delta_Left),fabs(Delta_Right));
        		Solution_i[iVar] = Two_phase_i[iVar] + config->GetLimiterCoeff_2phase() * Delta;

        		Delta_Left = Sol_Right-Two_phase_j[iVar]; Delta_Right = (Two_phase_j[iVar]-Two_phase_i[iVar])/2;

        		if (Delta_Left > 0 && Delta_Right > 0) {
        		    Delta    = 1.0;
        		} else if (Delta_Left < 0 && Delta_Right < 0) {
        			Delta = -1.0;
        		} else {
        			Delta = 0;
        		}
        		Delta = Delta *min(fabs(Delta_Left),fabs(Delta_Right));
        		Solution_j[iVar] = Two_phase_j[iVar] - config->GetLimiterCoeff_2phase() * Delta;


        	} else if (config->GetKind_SlopeLimit_2phase() == VAN_ALBADA) {
                // VAN_ALBADA LIMITER

        	   Delta_Left = Two_phase_j[iVar]-Two_phase_i[iVar];

               Delta = Project_Grad_i*Delta_Left* (2*Project_Grad_i + Delta_Left);
               Delta = Delta/(4*Project_Grad_i*Project_Grad_i+ Delta_Left*Delta_Left+1e-30);

               Solution_i[iVar] = Two_phase_i[iVar] + config->GetLimiterCoeff_2phase() * Delta;

               Delta = - Project_Grad_j*Delta_Left*(-2*Project_Grad_j + Delta_Left);
               Delta = Delta/(4*Project_Grad_j*Project_Grad_j+Delta_Left*Delta_Left+1e-30);

               Solution_j[iVar] = Two_phase_j[iVar] - config->GetLimiterCoeff_2phase() * Delta;


        	} else if (config->GetKind_SlopeLimit_2phase() == PUT) {

           		Delta_Left = Two_phase_i[iVar]-Sol_Left; Delta_Right = (Two_phase_j[iVar]-Two_phase_i[iVar])/2;

           		Delta =  (2* Delta_Left* Delta_Right +1e-30)/(pow(Delta_Left,2) + pow(Delta_Right, 2)+1e-30);
           		Delta = 0.25 * Delta * ( (1-0.33*Delta) * Delta_Right + (1+0.33*Delta) * Delta_Left) ;
        		Solution_i[iVar] = Two_phase_i[iVar] + config->GetLimiterCoeff_2phase() * Delta;

        		Delta_Left = Sol_Right-Two_phase_j[iVar]; Delta_Right = (Two_phase_j[iVar]-Two_phase_i[iVar])/2;

        		Delta =  (2* Delta_Left* Delta_Right +1e-30)/(pow(Delta_Left,2) + pow(Delta_Right, 2)+1e-30);
        		Delta = 0.25 * Delta * ( (1-0.33*Delta) * Delta_Left + (1+0.33*Delta) * Delta_Right) ;
        		Solution_j[iVar] = Two_phase_j[iVar] - config->GetLimiterCoeff_2phase() * Delta;


        	} else if (config->GetKind_SlopeLimit_2phase() == SUPERBEE) {

        		Delta_Left = Two_phase_i[iVar]-Sol_Left; Delta_Right = (Two_phase_j[iVar]-Two_phase_i[iVar])/2;

        		if (Delta_Left > 0 && Delta_Right > 0) {
        		    Delta    = 1.0;
        		} else if (Delta_Left < 0 && Delta_Right < 0) {
        			Delta = -1.0;
        		} else {
        			Delta = 0;
        		}
        		Delta = Delta *max( min(fabs(Delta_Left), 2*fabs(Delta_Right)) , min(fabs(Delta_Right), 2*fabs(Delta_Left)));
        		Solution_i[iVar] = Two_phase_i[iVar] + config->GetLimiterCoeff_2phase() * Delta;

        		Delta_Left = Sol_Right-Two_phase_j[iVar]; Delta_Right = (Two_phase_j[iVar]-Two_phase_i[iVar])/2;

        		if (Delta_Left > 0 && Delta_Right > 0) {
        		    Delta    = 1.0;
        		} else if (Delta_Left < 0 && Delta_Right < 0) {
        			Delta = -1.0;
        		} else {
        			Delta = 0;
        		}
        		Delta = Delta *max( min(fabs(Delta_Left), 2*fabs(Delta_Right)) , min(fabs(Delta_Right), 2*fabs(Delta_Left)));
        		Solution_j[iVar] = Two_phase_j[iVar] - config->GetLimiterCoeff_2phase() * Delta;

        	} else {
        		Solution_i[iVar] = Two_phase_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
        		Solution_j[iVar] = Two_phase_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        	}

        }        else {
          Solution_i[iVar] = Two_phase_i[iVar] + Project_Grad_i;
          Solution_j[iVar] = Two_phase_j[iVar] + Project_Grad_j;
        }

        Solution_i[iVar] = max(Solution_i[iVar], 0.0);
        Solution_j[iVar] = max(Solution_j[iVar], 0.0);

      }



      numerics->Set2phaseVar(Solution_i, Solution_j);
    }

    /*--- Add and subtract residual if in the metastable region---*/

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

void C2phase_HillSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

  unsigned long iPoint, iNode;
  su2double r, s, y, *Liquid_vec, rho_v, rho_m, mom0, mom1, mom3;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
     rho_v   = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
     mom0    = solver_container[TWO_PHASE_SOL]->node[iPoint]->GetSolution(0);
     mom1    = solver_container[TWO_PHASE_SOL]->node[iPoint]->GetSolution(1);
     mom3    = solver_container[TWO_PHASE_SOL]->node[iPoint]->GetSolution(3);

     Liquid_vec = node[iPoint]->GetLiquidPrim();

     if (mom0 != 0.0 && mom1!=0.0 && Liquid_vec[1]!=0) {
    	 r = mom1 / mom0;

    	 y = mom3*(Liquid_vec[1] - rho_v);
     	 y = y + 0.75 * rho_v / 3.1415;
     	 y = mom3*Liquid_vec[1] / y;

     	 rho_m = y/ Liquid_vec[1] + (1.0 - y)/ rho_v;
     	 rho_m = 1.0/ rho_m;

     	 s = 3.0 * rho_m * y * Liquid_vec[9];
     	 s = s/r;

     } else {
    	 r = 0; y = 0; rho_m = rho_v; s = 0;
     }

	 node[iPoint]->SetSource(s);
	 node[iPoint]->SetLiqEnthalpy(Liquid_vec[2]);
	 node[iPoint]->SetRadius(r);
     node[iPoint]->SetLiquidFrac(y);
    
  }

}


void C2phase_HillSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {

}

void C2phase_HillSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) {
  
  su2double *Two_phase_sol, *Prim_vec, *Liquid_vec;
  unsigned long iPoint;
  unsigned short iVar, jVar;


  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

	Prim_vec      = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
	Two_phase_sol = solver_container[TWO_PHASE_SOL]->node[iPoint]->GetSolution();

	Liquid_vec = node[iPoint]->GetLiquidPrim();
	Liquid_vec = node[iPoint]->SetLiquidPrim(Prim_vec, Two_phase_sol, Liquid_vec[6], solver_container[FLOW_SOL]->GetFluidModel(), config);

	numerics->Set2phaseVar(Two_phase_sol, NULL);

	numerics->SetPrimitive(Prim_vec, NULL);

	numerics->SetVolume(geometry->node[iPoint]->GetVolume());

    numerics->ComputeResidual(Residual, Jacobian_i, Liquid_vec, NULL, config);
    node[iPoint]->SetLiquidPrim(Liquid_vec);

    /*--- Subtract residual and the Jacobian ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

  }

}


void C2phase_HillSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex, total_index;
  unsigned short iDim, iVar;
/*
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    //--- Check if the node belongs to the domain (i.e, not a halo node) ---
    if (geometry->node[iPoint]->GetDomain()) {
      
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        }

      //--- Set the solution values and zero the residual ---
       node[iPoint]->SetSolution_Old(Solution_i);
       node[iPoint]->SetSolution(Solution_i);

//      LinSysRes.SetBlock_Zero(iPoint);

      //--- Change rows of the Jacobian (includes 1 in the diagonal) ---
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }

    }
  }
*/
}

void C2phase_HillSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                        unsigned short val_marker) {
  
  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar;
/*
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    //--- Check if the node belongs to the domain (i.e, not a halo node) ---
    if (geometry->node[iPoint]->GetDomain()) {

      Solution_i = node[iPoint]-> GetSolution();
      
      //--- Set the solution values and zero the residual ---
      node[iPoint]->SetSolution_Old(Solution_i);
      node[iPoint]->SetSolution(Solution_i);
      LinSysRes.SetBlock_Zero(iPoint);
      
      //--- Change rows of the Jacobian (includes 1 in the diagonal) ---
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
      
    }
  }
*/
}

void C2phase_HillSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  su2double *Normal, *V_infty, *V_domain;
  unsigned short iVar, iDim;
  
  bool grid_movement = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the infinity ---*/
      
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Set turbulent variable at the wall, and at infinity ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
      Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      Solution_j[iVar] = node[iPoint]->GetSolution(iVar);

      
      conv_numerics->Set2phaseVar(Solution_i, Solution_j);
      
      /*--- Set Normal (it is necessary to change the sign) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
      Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute residuals and Jacobians ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add residuals and Jacobians ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  delete [] Normal;
  
}

void C2phase_HillSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Riemann(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT: case STATIC_SUPERSONIC_INFLOW_PT: case STATIC_SUPERSONIC_INFLOW_PD: case DENSITY_VELOCITY:
    BC_Inlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  case STATIC_PRESSURE:
    BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}


void C2phase_HillSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                              unsigned short val_marker) {
  
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_inlet, *V_domain, *Normal;

  Normal = new su2double[nDim];

  bool grid_movement  = config->GetGrid_Movement();

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = 0.0;
      }

      /*--- Allocate the value at the inlet ---*/
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      conv_numerics->Set2phaseVar(Solution_i, Solution_j);
      
      /*--- Set various other quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}

void C2phase_HillSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex, Point_Normal;
  unsigned short iVar, iDim;
  su2double *V_outlet, *V_domain, *Normal;
  
  bool grid_movement  = config->GetGrid_Movement();

  Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      V_outlet = V_domain;

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*Solution_i --> 2phaseVar_internal,
      Solution_j --> 2phaseVar_outlet ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      }
      conv_numerics->Set2phaseVar(Solution_i, Solution_i);
      
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

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}




C2phase_QMOMSolver::C2phase_QMOMSolver(void) : C2phaseSolver() {

}

C2phase_QMOMSolver::C2phase_QMOMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : C2phaseSolver() {
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  ifstream restart_file;
  string text_line;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Array initialization ---*/


  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Dimension of the problem --> dependent on the 2phase model. ---*/

  nVar = 4;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];

  /*--- Single grid simulation ---*/

  if (iMesh == MESH_0) {

    /*--- Define some auxiliary vector related with the residual ---*/

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

    /*--- Define some auxiliary vector related with the solution ---*/

    Solution = new su2double[nVar];
    Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];


    /*--- Define some auxiliary vectors related to the liquid phase ---*/

    Primitive_Liquid = new su2double[10];

    /*--- Define some auxiliary vector related with the geometry ---*/

    Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];

    /*--- Define some auxiliary vector related with the flow solution ---*/

    FlowPrimVar_i = new su2double [nDim+9]; FlowPrimVar_j = new su2double [nDim+9];

    /*--- Jacobians and vector structures for implicit computations ---*/

    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (QMOM model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  }

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

  lowerlimit[0] = 1.0e-10;
  upperlimit[0] = 1.0e10;

  lowerlimit[1] = 1.0e-4;
  upperlimit[1] = 1.0e15;

  /*--- Far-field flow state quantities and initialization. ---*/
  su2double RInf, NInf, DInf;

  DInf    = config->GetDensity_FreeStreamND();
  RInf    = 0;
  NInf    = 0;

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new C2phase_HillVariable(RInf, NInf, DInf, nDim, nVar, config);

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);

}

C2phase_QMOMSolver::~C2phase_QMOMSolver(void) {

}

void C2phase_QMOMSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint;

  unsigned long ExtIter      = config->GetExtIter();
  bool limiter_flow          = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Initialize the residual vector ---*/

    LinSysRes.SetBlock_Zero(iPoint);

  }

  /*--- Initialize the Jacobian matrices ---*/

  Jacobian.SetValZero();

  /*--- Upwind second order reconstruction ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

  if (config->GetSpatialOrder() == SECOND_ORDER_LIMITER) SetSolution_Limiter(geometry, config);

  if (limiter_flow) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);

}

void C2phase_QMOMSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {

	  su2double *Two_phase_i, *Two_phase_j, *Limiter_i = NULL, *Limiter_j = NULL, *V_i, *V_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
	  su2double *Liquid_vec, *Solution_L, *Solution_Right;
	  unsigned long iEdge, iPoint, jPoint;
	  unsigned short iDim, iVar;

	  bool second_order  = ((config->GetSpatialOrder() == SECOND_ORDER) || (config->GetSpatialOrder() == SECOND_ORDER_LIMITER));
	  bool limiter       = (config->GetSpatialOrder() == SECOND_ORDER_LIMITER);
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

	    /*--- 2phase variables w/o reconstruction ---*/

	    Two_phase_i = node[iPoint]->GetSolution();
	    Two_phase_j = node[jPoint]->GetSolution();
	    numerics->Set2phaseVar(Two_phase_i, Two_phase_j);

	    Liquid_vec = node[iPoint]->SetLiquidPrim(V_i, Two_phase_i, numerics->Primitive_Liquid[6], FluidModel, config);
	    numerics->SetPrimitive_Liquid(Liquid_vec);

	    /*--- Grid Movement ---*/

	    if (grid_movement)
	      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());

	    if (second_order) {

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


	        } else {
	          FlowPrimVar_i[iVar] = V_i[iVar] + Project_Grad_i;
	          FlowPrimVar_j[iVar] = V_j[iVar] + Project_Grad_j;
	        }
	      }

	      numerics->SetPrimitive(FlowPrimVar_i, FlowPrimVar_j);

	      /*--- Turbulent variables using gradient reconstruction and limiters ---*/

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
					Solution_i[iVar] = Two_phase_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
					Solution_j[iVar] = Two_phase_j[iVar] + Limiter_j[iVar]*Project_Grad_j;

	        } else {
	          Solution_i[iVar] = Two_phase_i[iVar] + Project_Grad_i;
	          Solution_j[iVar] = Two_phase_j[iVar] + Project_Grad_j;
	        }
	      }

	      numerics->Set2phaseVar(Solution_i, Solution_j);

	//      numerics->SetPrimitive_Liquid(node[iPoint]->SetLiquidPrim(FlowPrimVar_i, Solution_i, config));

	    }

	    /*--- Add and subtract residual ---*/

	    cout << "Solution before upwind " << Two_phase_i[0] << endl;
	    getchar();

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


void C2phase_QMOMSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

	 unsigned long iPoint, iNode, rho_m;
	  su2double R, S, y, *Liquid_vec, rho_v, rhoN, rhoNR, rhoNR3;

	  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
	  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

	  /*--- Compute mean flow and turbulence gradients ---*/

	  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
	//    solver_container[FLOW_SOL]->SetPrimitive_Gradient_GG(geometry, config);
	    SetSolution_Gradient_GG(geometry, config);
	  }
	  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
	//    solver_container[FLOW_SOL]->SetPrimitive_Gradient_LS(geometry, config);
	    SetSolution_Gradient_LS(geometry, config);
	  }

	  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

	     rho_v  = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
	     rhoN      = node[iPoint]->GetSolution(0);
	     rhoNR     = node[iPoint]->GetSolution(1);
	     rhoNR3    = node[iPoint]->GetSolution(3);

	     Liquid_vec= node[iPoint]->GetLiquidPrim();

	     if (rhoN != 0) {

	    	 R = rhoNR / rhoN;
	     	 y = rhoNR3*(Liquid_vec[1] - rho_v);
	     	 y = y + 0.75 * rho_v / 3.14;
	     	 y = y*Liquid_vec[1] / y;

	     	 rho_m = y/ Liquid_vec[1] + (1.0 - y)/ rho_v;
	     	 rho_m = 1.0/ rho_m;

	     	 S = rho_m * 3 * y / R * Liquid_vec[9];
	     } else {
	    	 R = 0; y = 0; rho_m = rho_v; S = 0;
	     }

		 node[iPoint]->SetSource(S);
		 node[iPoint]->SetLiqEnthalpy(Liquid_vec[2]);
		 node[iPoint]->SetRadius(R);
	     node[iPoint]->SetLiquidFrac(y);

	  }

}

void C2phase_QMOMSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) {

  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(), NULL);


    /*--- Set volume ---*/

    numerics->SetVolume(geometry->node[iPoint]->GetVolume());


    /*--- Compute the source term ---*/

    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);


    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

  }

}


void C2phase_QMOMSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex, total_index;
  unsigned short iDim, iVar;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {


      /*--- Set the solution values and zero the residual ---*/
      node[iPoint]->SetSolution_Old(Solution_j);
      node[iPoint]->SetSolution(Solution_j);
      LinSysRes.SetBlock_Zero(iPoint);

      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }

    }
  }

}

void C2phase_QMOMSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                        unsigned short val_marker) {

  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Set the solution values and zero the residual ---*/
      node[iPoint]->SetSolution_Old(Solution_j);
      node[iPoint]->SetSolution(Solution_j);
      LinSysRes.SetBlock_Zero(iPoint);

      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }

    }
  }

}

void C2phase_QMOMSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;
  su2double *Normal, *V_infty, *V_domain;
  unsigned short iVar, iDim;

  bool grid_movement = config->GetGrid_Movement();

  Normal = new su2double[nDim];

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Set turbulent variable at the wall, and at infinity ---*/

      for (iVar = 0; iVar < nVar; iVar++)
      Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      Solution_j[iVar] = node[iPoint]->GetSolution(iVar);


      conv_numerics->Set2phaseVar(Solution_i, Solution_j);

      /*--- Set Normal (it is necessary to change the sign) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
      Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Grid Movement ---*/

      if (grid_movement)
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute residuals and Jacobians ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Add residuals and Jacobians ---*/

      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  delete [] Normal;

}

void C2phase_QMOMSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Riemann(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT: case STATIC_SUPERSONIC_INFLOW_PT: case STATIC_SUPERSONIC_INFLOW_PD: case DENSITY_VELOCITY:
    BC_Inlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  case STATIC_PRESSURE:
    BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}

void C2phase_QMOMSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                              unsigned short val_marker) {

  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_inlet, *V_domain, *Normal;

  Normal = new su2double[nDim];

  bool grid_movement  = config->GetGrid_Movement();

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      for (iVar = 0; iVar < nVar; iVar++)
      Solution_i[iVar] = node[iPoint]->GetSolution(iVar);

      Solution_j[0]= 0;
      Solution_j[1]= 0;
      Solution_j[2]= 0;
      Solution_j[3]= 0;

      conv_numerics->Set2phaseVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);

      if (grid_movement)
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}

void C2phase_QMOMSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex, Point_Normal;
  unsigned short iVar, iDim;
  su2double *V_outlet, *V_domain, *Normal;

  bool grid_movement  = config->GetGrid_Movement();

  Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Allocate the value at the outlet ---*/
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);

      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> TurbVar_internal,
       Solution_j --> TurbVar_outlet ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
      }
      conv_numerics->Set2phaseVar(Solution_i, Solution_j);

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

      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_outlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->Set2phaseVar(Solution_i, Solution_j);
      visc_numerics->Set2phaseVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());

      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}



