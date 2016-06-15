/*!
 * \file solution_adjoint_levelset.cpp
 * \brief Main subrotuines for solving the level set problem.
 * \author F. Palacios
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

#include "../include/solver_structure.hpp"

CAdjLevelSetSolver::CAdjLevelSetSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
	unsigned short iVar, iDim, nLineLets;
	unsigned long iPoint, index;
	su2double dull_val;
	ifstream restart_file;
	string text_line, mesh_filename, filename, AdjExt;
	
  bool restart = config->GetRestart();
	
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	node = new CVariable*[geometry->GetnPoint()];
	
	/*--- Dimension of the problem ---*/
	nVar = 1;
  
  /*--- Single grid simulation ---*/
	if (iMesh == MESH_0) {
	
    /*--- Define some auxiliar vector related with the residual ---*/
    Residual      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
    Residual_RMS  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
    Residual_Max  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
    Point_Max  = new unsigned long[nVar]; for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
    Point_Max_Coord = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
    }
    Residual_i    = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
    Residual_j    = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
	
    /*--- Define some auxiliar vector related with the solution ---*/
    Solution    = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]    = 0.0;
    Solution_i  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar]  = 0.0;
    Solution_j  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar]  = 0.0;
    
    /*--- Define some auxiliar vector related with the geometry ---*/
    Vector    = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]    = 0.0;
    Vector_i  = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim]  = 0.0;
    Vector_j  = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim]  = 0.0;
	
    /*--- Define some auxiliar vector related with the flow solution ---*/
    FlowSolution_i = new su2double[nDim+1]; for (iVar = 0; iVar < nDim+1; iVar++) FlowSolution_i[iVar]  = 0.0;
    FlowSolution_j = new su2double[nDim+1]; for (iVar = 0; iVar < nDim+1; iVar++) FlowSolution_j[iVar]  = 0.0;
	
    /*--- Solution and residual vectors ---*/
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    
    /*--- Jacobians and vector structures for implicit computations ---*/
    if (config->GetKind_TimeIntScheme_AdjLevelSet() == EULER_IMPLICIT) {
      
      /*--- Point to point Jacobians ---*/
      Jacobian_i = new su2double* [nVar];
      Jacobian_j = new su2double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[iVar] = new su2double [nVar];
        Jacobian_j[iVar] = new su2double [nVar];
      }
      
      Jacobian_ii = new su2double* [nVar];
      Jacobian_ij = new su2double* [nVar];
      Jacobian_ji = new su2double* [nVar];
      Jacobian_jj = new su2double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_ii[iVar] = new su2double [nVar];
        Jacobian_ij[iVar] = new su2double [nVar];
        Jacobian_ji[iVar] = new su2double [nVar];
        Jacobian_jj[iVar] = new su2double [nVar];
      }
      
      /*--- Initialization of the structure of the whole Jacobian ---*/
      if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Adj. Level Set). MG level: 0." << endl;
      Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
      
      if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
          (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
        nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
        if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
      }
      
    }
	
    /*--- Computation of gradients by least squares ---*/
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
      Smatrix = new su2double* [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Smatrix[iDim] = new su2double [nDim];
      /*--- c vector := transpose(WA)*(Wb) ---*/
      cvector = new su2double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++)
        cvector[iVar] = new su2double [nDim];
    }
    
  }
	
	/*--- Restart the solution from file information ---*/
	if (!restart || (iMesh != MESH_0)) {
    
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			node[iPoint] = new CAdjLevelSetVariable(0.0, nDim, nVar, config);
		}
    
	}
	else {

		/*--- Restart the solution from file information ---*/
		mesh_filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(mesh_filename);
    
		restart_file.open(filename.data(), ios::in);		
		
		if (restart_file.fail()) {
			cout << "There is no level set restart file!!" << endl;
			exit(EXIT_FAILURE);
		}
    
    /*--- In case this is a parallel simulation, we need to perform the 
     Global2Local index transformation first. ---*/
    long *Global2Local;
    Global2Local = new long[geometry->GetGlobal_nPointDomain()];
    /*--- First, set all indices to a negative value by default ---*/
    for (iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
      Global2Local[iPoint] = -1;
    }
    /*--- Now fill array with the transform values only for local points ---*/
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    }

    /*--- Read all lines in the restart file ---*/
    long iPoint_Local; unsigned long iPoint_Global = 0;
    
    /*--- The first line is the header ---*/
    getline (restart_file, text_line);
    
    while (getline (restart_file, text_line)) {
			istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives 
       on a different processor, the value of iPoint_Local will be -1. 
       Otherwise, the local index for this node on the current processor 
       will be returned and used to instantiate the vars. ---*/
      iPoint_Local = Global2Local[iPoint_Global];
      if (iPoint_Local >= 0) {
        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        node[iPoint_Local] = new CAdjLevelSetVariable(Solution[0], nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
      node[iPoint] = new CAdjLevelSetVariable(Solution[0], nDim, nVar, config);
    }
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
	}
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the gradient (Mean flow source term) ---*/
  if (iMesh == MESH_0) {
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  }
  
}

CAdjLevelSetSolver::~CAdjLevelSetSolver(void) { }

void CAdjLevelSetSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
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
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
	}
  
}

void CAdjLevelSetSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
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
          Solution[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetLimiter(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
	}
}

void CAdjLevelSetSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;

#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
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

void CAdjLevelSetSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
	unsigned long iPoint;
	
	bool implicit      = (config->GetKind_TimeIntScheme_AdjLevelSet() == EULER_IMPLICIT);
	bool second_order  = ((config->GetSpatialOrder_AdjLevelSet() == SECOND_ORDER) || (config->GetSpatialOrder_AdjLevelSet() == SECOND_ORDER_LIMITER));

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		LinSysRes.SetBlock_Zero(iPoint);
	
	/*--- Implicit part ---*/
	if (implicit)
		Jacobian.SetValZero();
	
  if (second_order) {
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  }
  
}

void CAdjLevelSetSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {
	su2double *LevelSet_var_i, *LevelSet_var_j, *U_i, *U_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, DensityInc_i, DensityInc_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	
	bool implicit      = (config->GetKind_TimeIntScheme_AdjLevelSet() == EULER_IMPLICIT);
	bool second_order  = ((config->GetSpatialOrder_AdjLevelSet() == SECOND_ORDER) || (config->GetSpatialOrder_AdjLevelSet() == SECOND_ORDER_LIMITER));
		
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
		
		/*--- Conservative variables w/o reconstruction ---*/
		U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
		U_j = solver_container[FLOW_SOL]->node[jPoint]->GetSolution();
		numerics->SetConservative(U_i, U_j);
		
		/*--- Level Set variables w/o reconstruction ---*/
		LevelSet_var_i = node[iPoint]->GetSolution();
		LevelSet_var_j = node[jPoint]->GetSolution();
		numerics->SetLevelSetVar(LevelSet_var_i, LevelSet_var_j);
		
		/*--- Set the value of the density ---*/
		DensityInc_i = solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
		DensityInc_j = solver_container[FLOW_SOL]->node[jPoint]->GetDensityInc();
		numerics->SetDensityInc(DensityInc_i, DensityInc_j);
				
		if (second_order) {

			for (iDim = 0; iDim < nDim; iDim++) {
				Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
				Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}
			
			/*--- Level Set variables using gradient reconstruction ---*/
			Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0.0; Project_Grad_j = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				Solution_i[iVar] = LevelSet_var_i[iVar] + Project_Grad_i;
				Solution_j[iVar] = LevelSet_var_j[iVar] + Project_Grad_j;
			}
			numerics->SetLevelSetVar(Solution_i, Solution_j);
      
		}
		
		numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);	
		
		/*--- Add and subtract Residual ---*/
		LinSysRes.AddBlock(iPoint, Residual_i);
		LinSysRes.AddBlock(jPoint, Residual_j);
    		
		/*--- Implicit part ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
			Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
			Jacobian.AddBlock(jPoint, iPoint, Jacobian_ji);
			Jacobian.AddBlock(jPoint, jPoint, Jacobian_jj);
		}		
		
	}
}

void CAdjLevelSetSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
																								 CConfig *config, unsigned short iMesh) {
//	unsigned short iVar, iDim;
//	unsigned long iPoint;
//  su2double epsilon, DeltaDirac, lambda, dRho_dPhi, dMud_Phi, Vol, DiffLevelSet, LevelSet, *AdjMeanFlow, **AdjLevelSetGradient, **AdjMeanFlowGradient,
//  Density, Velocity[3] = {0.0,0.0,0.0}, ProjAdj, dFc_dRho[3][4], ProjFlux;
//  
//	su2double Froude2 = config->GetFroude()*config->GetFroude();
//  
//	for (iVar = 0; iVar < nVar; iVar++)
//		Residual[iVar] = 0;
//
//	/*--- loop over points ---*/
//  
//	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 
//		
//		Vol = geometry->node[iPoint]->GetVolume();
//
//    /*--- Direct problem quantities ---*/
//    
//		DiffLevelSet = solver_container[LEVELSET_SOL]->node[iPoint]->GetDiffLevelSet();
//		LevelSet = solver_container[LEVELSET_SOL]->node[iPoint]->GetSolution(0);
//		MeanFlow = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
//    Density = solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
//    
//    /*--- Adjoint problem quantities ---*/
//    
//		AdjMeanFlow = solver_container[ADJFLOW_SOL]->node[iPoint]->GetSolution();
//    AdjLevelSetGradient = solver_container[ADJLEVELSET_SOL]->node[iPoint]->GetGradient();
//    AdjMeanFlowGradient = solver_container[ADJFLOW_SOL]->node[iPoint]->GetGradient();
//    
//    /*--- Projected adjoint velocity ---*/
//    
//    ProjAdj = 0.0;
//    for (iDim = 0; iDim < nDim; iDim++) {
//      Velocity[iDim] = solver_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim);
//      ProjAdj += Velocity[iDim]*AdjLevelSetGradient[0][iDim];
//    }
//    
//    /*--- Compute the flow solution using the level set value. ---*/
//    
//    epsilon = config->GetFreeSurface_Thickness();
//    DeltaDirac = 0.0;
//    if (fabs(LevelSet) <= epsilon) DeltaDirac = 1.0 - (0.5*(1.0+(LevelSet/epsilon)+(1.0/PI_NUMBER)*sin(PI_NUMBER*LevelSet/epsilon)));
//    
//    /*--- Set the value of the incompressible density for free surface flows (density ratio g/l) ---*/
//    
//    lambda = config->GetRatioDensity();
//    dRho_dPhi = (1.0 - lambda)*DeltaDirac*config->GetDensity_FreeStreamND();
//    
//    /*--- Set the value of the incompressible viscosity for free surface flows (viscosity ratio g/l) ---*/
//    
//    lambda = config->GetRatioViscosity();
//    dMud_Phi = (1.0 - lambda)*DeltaDirac*config->GetViscosity_FreeStreamND();
//
//    /*--- Flux derivative ---*/
//    
//    for (iDim = 0; iDim < nDim; iDim++) {
//      dFc_dRho[iDim][0] = 0.0;
//      dFc_dRho[iDim][1] = Velocity[iDim]*Velocity[0];
//      dFc_dRho[iDim][2] = Velocity[iDim]*Velocity[1];
//      dFc_dRho[iDim][3] = Velocity[iDim]*Velocity[2];
//    }
//    
//    /*--- Projected flux derivative ---*/
//    
//    ProjFlux = 0.0;
//    for (iDim = 0; iDim < nDim; iDim++) {
//      for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
//        ProjFlux += AdjMeanFlowGradient[iVar][iDim]*dFc_dRho[iDim][iVar];
//      }
//    }
//    
// 		Residual[0] = ( ProjFlux -(LevelSet*ProjAdj/Density) - (AdjMeanFlow[nDim]/Froude2))*dRho_dPhi*Vol;
//		
//    /*--- The Free surface objective function requires the levelt set difference ---*/
//    
//    if (config->GetKind_ObjFunc() == FREE_SURFACE) Residual[0] +=  DiffLevelSet* Vol;
//
//		/*--- Add Residual ---*/
//		LinSysRes.AddBlock(iPoint, Residual);
//	}
//	
}

void CAdjLevelSetSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) { }

void CAdjLevelSetSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	su2double *U_domain, *U_wall, *LevelSet_domain, *LevelSet_wall;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new su2double[solver_container[FLOW_SOL]->GetnVar()]; 
	U_wall = new su2double[solver_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new su2double[1];
	LevelSet_wall = new su2double[1];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
						
			for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
				U_domain[iVar] = solver_container[FLOW_SOL]->node[Point_Normal]->GetSolution(iVar);
				U_wall[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			}
			
			LevelSet_domain[0] = node[Point_Normal]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			numerics->SetNormal(Vector);
			
			numerics->SetConservative(U_wall, U_wall);
			numerics->SetLevelSetVar(LevelSet_wall, LevelSet_wall);
			numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);	
			
			LinSysRes.AddBlock(iPoint, Residual_i);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
			
		}
	}
}

void CAdjLevelSetSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                        CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	su2double *U_domain, *U_wall, *LevelSet_domain, *LevelSet_wall;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new su2double[solver_container[FLOW_SOL]->GetnVar()];
	U_wall = new su2double[solver_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new su2double[1];
	LevelSet_wall = new su2double[1];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
				U_domain[iVar] = solver_container[FLOW_SOL]->node[Point_Normal]->GetSolution(iVar);
				U_wall[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			}
			
			LevelSet_domain[0] = node[Point_Normal]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			conv_numerics->SetNormal(Vector);
			
			conv_numerics->SetConservative(U_wall, U_wall);
			conv_numerics->SetLevelSetVar(LevelSet_wall, LevelSet_wall);
			conv_numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			LinSysRes.AddBlock(iPoint, Residual_i);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
			
		}
	}
}

void CAdjLevelSetSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	su2double *U_domain, *U_wall, *LevelSet_domain, *LevelSet_wall;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new su2double[solver_container[FLOW_SOL]->GetnVar()];
	U_wall = new su2double[solver_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new su2double[1];
	LevelSet_wall = new su2double[1];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
				U_domain[iVar] = solver_container[FLOW_SOL]->node[Point_Normal]->GetSolution(iVar);
				U_wall[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			}
			
			LevelSet_domain[0] = node[Point_Normal]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			conv_numerics->SetNormal(Vector);
			
			conv_numerics->SetConservative(U_wall, U_wall);
			conv_numerics->SetLevelSetVar(LevelSet_wall, LevelSet_wall);
			conv_numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			LinSysRes.AddBlock(iPoint, Residual_i);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
			
		}
	}
}

void CAdjLevelSetSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	su2double *U_domain, *U_wall, *LevelSet_domain, *LevelSet_wall;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new su2double[solver_container[FLOW_SOL]->GetnVar()];
	U_wall = new su2double[solver_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new su2double[1];
	LevelSet_wall = new su2double[1];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
				U_domain[iVar] = solver_container[FLOW_SOL]->node[Point_Normal]->GetSolution(iVar);
				U_wall[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			}
			
			LevelSet_domain[0] = node[Point_Normal]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			conv_numerics->SetNormal(Vector);
			
			conv_numerics->SetConservative(U_wall, U_wall);
			conv_numerics->SetLevelSetVar(LevelSet_wall, LevelSet_wall);
			conv_numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			LinSysRes.AddBlock(iPoint, Residual_i);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
			
		}
	}
}

void CAdjLevelSetSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	su2double *U_domain, *U_wall, *LevelSet_domain, *LevelSet_wall;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new su2double[solver_container[FLOW_SOL]->GetnVar()];
	U_wall = new su2double[solver_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new su2double[1];
	LevelSet_wall = new su2double[1];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
				U_domain[iVar] = solver_container[FLOW_SOL]->node[Point_Normal]->GetSolution(iVar);
				U_wall[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			}
			
			LevelSet_domain[0] = node[Point_Normal]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			conv_numerics->SetNormal(Vector);
			
			conv_numerics->SetConservative(U_wall, U_wall);
			conv_numerics->SetLevelSetVar(LevelSet_wall, LevelSet_wall);
			conv_numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			LinSysRes.AddBlock(iPoint, Residual_i);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
			
		}
	}
}

void CAdjLevelSetSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	
	bool implicit = (config->GetKind_TimeIntScheme_AdjLevelSet() == EULER_IMPLICIT);
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set the solution to the original value ---*/
		Solution[0] = node[iPoint]->GetSolution(0);
		
		node[iPoint]->SetSolution_Old(Solution);
		LinSysRes.SetBlock_Zero(iPoint);
		
		/*--- Includes 1 in the diagonal ---*/
		if (implicit)
			Jacobian.DeleteValsRowi(iPoint);
	}
}

void CAdjLevelSetSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
	unsigned long iPoint;
	su2double Delta = 0.0, Vol;
	
	/*--- Set maximum residual to zero ---*/
  
    SetRes_RMS(0, 0.0);
    SetRes_Max(0, 0.0, 0);
    
	/*--- Build implicit system ---*/
  
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		/*--- Read the volume ---*/
    
		Vol = geometry->node[iPoint]->GetVolume();
		
		/*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
		Delta = Vol / solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time();
        
		Jacobian.AddVal2Diag(iPoint, Delta);
        
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
		LinSysRes[iPoint] = -LinSysRes[iPoint];
		LinSysSol[iPoint] = 0.0;
		AddRes_RMS(0, LinSysRes[iPoint]*LinSysRes[iPoint]);
    AddRes_Max(0, fabs(LinSysRes[iPoint]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
	}
    
    /*--- Initialize residual and solution at the ghost points ---*/
  
    for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
        LinSysRes[iPoint] = 0.0;
        LinSysSol[iPoint] = 0.0;
    }
	
  /*--- Solve or smooth the linear system ---*/
  
  CSysSolve system;
  system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
	
	/*--- Update solution (system written in terms of increments), be careful with the update of the
	 scalar equations which includes the density ---*/
  
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->AddSolution(0, LinSysSol[iPoint]);
    
    /*--- MPI solution ---*/
  
    Set_MPI_Solution(geometry, config);
    
    /*--- Compute the root mean square residual ---*/
  
    SetResidual_RMS(geometry, config);
  
}

void CAdjLevelSetSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
	su2double *U_time_nM1, *U_time_n, *U_time_nP1, Volume_nM1, Volume_n, Volume_nP1, TimeStep;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool Grid_Movement = config->GetGrid_Movement();
  
	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1  = node[iPoint]->GetSolution_time_n1();
		U_time_n    = node[iPoint]->GetSolution_time_n();
		U_time_nP1  = node[iPoint]->GetSolution();
    
		/*--- Volume at time n-1 and n ---*/
		if (Grid_Movement) {
			Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
			Volume_n = geometry->node[iPoint]->GetVolume_n();
			Volume_nP1 = geometry->node[iPoint]->GetVolume();
		}
		else {
			Volume_nM1 = geometry->node[iPoint]->GetVolume();
			Volume_n = geometry->node[iPoint]->GetVolume();
			Volume_nP1 = geometry->node[iPoint]->GetVolume();
		}
		    
		/*--- Time Step ---*/
		TimeStep = config->GetDelta_UnstTimeND();
		
		/*--- Compute Residual ---*/
    if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
      Residual[0] = ( U_time_nP1[0]*Volume_nP1 - U_time_n[0]*Volume_n ) / TimeStep;
    if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
      Residual[0] = ( 3.0*U_time_nP1[0]*Volume_nP1 - 4.0*U_time_n[0]*Volume_n
                     +  1.0*U_time_nM1[0]*Volume_nM1 ) / (2.0*TimeStep);
    
		/*--- Add Residual ---*/
		LinSysRes.AddBlock(iPoint, Residual);
		
		if (implicit) {
      Jacobian_i[0][0] = 0.0;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Jacobian_i[0][0] = Volume_nP1 / TimeStep;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
        Jacobian_i[0][0] = (Volume_nP1*3.0)/(2.0*TimeStep);
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
    
	}
}
