/*!
 * \file solution_adjoint_mean.cpp
 * \brief Main subrotuines for solving adjoint problems (Euler, Navier-Stokes, etc.).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/solver_structure.hpp"

CAdjTNE2EulerSolver::CAdjTNE2EulerSolver(void) : CSolver() {
  
  /*--- Array initialization ---*/
	Phi_Inf = NULL;
	Sens_Mach = NULL;
	Sens_AoA = NULL;
	Sens_Geo = NULL;
	Sens_Press = NULL;
	Sens_Temp = NULL;
	p1_Und_Lapl = NULL;
	p2_Und_Lapl = NULL;
  Jacobian_Axisymmetric = NULL;
	CSensitivity = NULL;
  
}

CAdjTNE2EulerSolver::CAdjTNE2EulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
	unsigned long iPoint, index, iVertex;
	string text_line, mesh_filename;
	unsigned short iDim, iVar, iMarker;
	ifstream restart_file;
	string filename, AdjExt;
  
	bool restart = config->GetRestart();
	bool axisymmetric = config->GetAxisymmetric();
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Array initialization ---*/
	Phi_Inf = NULL;
	Sens_Mach = NULL;
	Sens_AoA = NULL;
	Sens_Geo = NULL;
	Sens_Press = NULL;
	Sens_Temp = NULL;
	p1_Und_Lapl = NULL;
	p2_Und_Lapl = NULL;
  Jacobian_Axisymmetric = NULL;
	CSensitivity = NULL;
  
	/*--- Set the gamma value ---*/
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
  nMarker = config->GetnMarker_All();
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
	nVar = nDim + 2;
	node = new CVariable*[nPoint];
  
	/*--- Define some auxiliary vectors related to the residual ---*/
	Residual = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
	Residual_i = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
	Res_Conv_i = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_i[iVar]    = 0.0;
  Res_Visc_i = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_i[iVar]    = 0.0;
	Res_Conv_j = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_j[iVar]    = 0.0;
  Res_Visc_j = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_j[iVar]    = 0.0;
  
	/*--- Define some auxiliary vectors related to the solution ---*/
	Solution   = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
	Solution_i = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar]   = 0.0;
  Solution_j = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar]   = 0.0;
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
	Vector   = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
	Vector_i = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
	Vector_j = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
	/*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
	if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED) {
		p1_Und_Lapl = new double [nPoint];
		p2_Und_Lapl = new double [nPoint];
	}
  
	/*--- Define some auxiliary vectors related to the geometry ---*/
	Vector_i = new double[nDim]; Vector_j = new double[nDim];
  
  /*--- Point to point Jacobians. These are always defined because
   they are also used for sensitivity calculations. ---*/
  Jacobian_i = new double* [nVar];
  Jacobian_j = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new double [nVar];
    Jacobian_j[iVar] = new double [nVar];
  }
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
	/*--- Jacobians and vector structures for implicit computations ---*/
  if (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT) {
		Jacobian_ii = new double* [nVar];
		Jacobian_ij = new double* [nVar];
		Jacobian_ji = new double* [nVar];
		Jacobian_jj = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_ii[iVar] = new double [nVar];
			Jacobian_ij[iVar] = new double [nVar];
			Jacobian_ji[iVar] = new double [nVar];
			Jacobian_jj[iVar] = new double [nVar];
		}
    
    if (rank == MASTER_NODE)
      cout << "Initialize jacobian structure (Adjoint Euler). MG level: " << iMesh <<"." << endl;
		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
    
    if (axisymmetric) {
      Jacobian_Axisymmetric = new double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_Axisymmetric[iVar] = new double [nVar];
    }
  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No jacobian structure (Adjoint Euler). MG level: " << iMesh <<"." << endl;
  }
  
	/*--- Computation of gradients by least squares ---*/
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			cvector[iVar] = new double [nDim];
	}
  
	/*--- Sensitivity definition and coefficient in all the markers ---*/
	CSensitivity = new double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		CSensitivity[iMarker] = new double [geometry->nVertex[iMarker]];
	}
	Sens_Geo  = new double[nMarker];
	Sens_Mach = new double[nMarker];
	Sens_AoA  = new double[nMarker];
	Sens_Press = new double[nMarker];
	Sens_Temp  = new double[nMarker];
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		Sens_Geo[iMarker]  = 0.0;
		Sens_Mach[iMarker] = 0.0;
		Sens_AoA[iMarker]  = 0.0;
		Sens_Press[iMarker] = 0.0;
		Sens_Temp[iMarker]  = 0.0;
		for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
			CSensitivity[iMarker][iVertex] = 0.0;
	}
  
	/*--- Adjoint flow at the inifinity, initialization stuff ---*/
	PsiRho_Inf = 0.0; PsiE_Inf   = 0.0;
	Phi_Inf    = new double [nDim];
	Phi_Inf[0] = 0.0; Phi_Inf[1] = 0.0;
	if (nDim == 3) Phi_Inf[2] = 0.0;
  
	if (!restart || geometry->GetFinestMGLevel() == false) {
		/*--- Restart the solution from infinity ---*/
		for (iPoint = 0; iPoint < nPoint; iPoint++)
			node[iPoint] = new CAdjTNE2EulerVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, nDim, nVar, config);
	}
	else {
    
		/*--- Restart the solution from file information ---*/
		mesh_filename = config->GetSolution_AdjFileName();
    
		/*--- Change the name, depending of the objective function ---*/
		filename.assign(mesh_filename);
		filename.erase (filename.end()-4, filename.end());
		switch (config->GetKind_ObjFunc()) {
      case DRAG_COEFFICIENT: AdjExt = "_cd.dat"; break;
      case LIFT_COEFFICIENT: AdjExt = "_cl.dat"; break;
      case SIDEFORCE_COEFFICIENT: AdjExt = "_csf.dat"; break;
      case PRESSURE_COEFFICIENT: AdjExt = "_cp.dat"; break;
      case MOMENT_X_COEFFICIENT: AdjExt = "_cmx.dat"; break;
      case MOMENT_Y_COEFFICIENT: AdjExt = "_cmy.dat"; break;
      case MOMENT_Z_COEFFICIENT: AdjExt = "_cmz.dat"; break;
      case EFFICIENCY: AdjExt = "_eff.dat"; break;
      case EQUIVALENT_AREA: AdjExt = "_ea.dat"; break;
      case NEARFIELD_PRESSURE: AdjExt = "_nfp.dat"; break;
      case FORCE_X_COEFFICIENT: AdjExt = "_cfx.dat"; break;
      case FORCE_Y_COEFFICIENT: AdjExt = "_cfy.dat"; break;
      case FORCE_Z_COEFFICIENT: AdjExt = "_cfz.dat"; break;
      case THRUST_COEFFICIENT: AdjExt = "_ct.dat"; break;
      case TORQUE_COEFFICIENT: AdjExt = "_cq.dat"; break;
      case FIGURE_OF_MERIT: AdjExt = "_merit.dat"; break;
      case FREE_SURFACE: AdjExt = "_fs.dat"; break;
      case NOISE: AdjExt = "_fwh.dat"; break;
      case HEAT_LOAD: AdjExt = "_Q.dat"; break;
		}
		filename.append(AdjExt);
		restart_file.open(filename.data(), ios::in);
    
		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
			cout << "There is no adjoint restart file!! " << filename.data() << "."<< endl;
			cout << "Press any key to exit..." << endl;
			cin.get(); exit(1);
		}
    
		/*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
		long *Global2Local;
		Global2Local = new long[geometry->GetGlobal_nPointDomain()];
		/*--- First, set all indices to a negative value by default ---*/
		for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
			Global2Local[iPoint] = -1;
		}
		/*--- Now fill array with the transform values only for local points ---*/
		for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
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
        if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
        if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
				node[iPoint_Local] = new CAdjTNE2EulerVariable(Solution, nDim, nVar, config);
			}
			iPoint_Global++;
		}
    
		/*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CAdjTNE2EulerVariable(Solution, nDim, nVar, config);
		}
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
	}
  
	/*--- Define solver parameters needed for execution of destructor ---*/
	if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED) space_centered = true;
	else space_centered = false;
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
}

CAdjTNE2EulerSolver::~CAdjTNE2EulerSolver(void) {
  unsigned short iVar, iMarker;
  
  if (Phi_Inf != NULL) delete [] Phi_Inf;
	if (Sens_Mach != NULL) delete [] Sens_Mach;
	if (Sens_AoA != NULL) delete [] Sens_AoA;
	if (Sens_Geo != NULL) delete [] Sens_Geo;
	if (Sens_Press != NULL) delete [] Sens_Press;
	if (Sens_Temp != NULL) delete [] Sens_Temp;
	if (p1_Und_Lapl != NULL) delete [] p1_Und_Lapl;
	if (p2_Und_Lapl != NULL) delete [] p2_Und_Lapl;
  
  
	if (Jacobian_Axisymmetric != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Jacobian_Axisymmetric[iVar];
    delete [] Jacobian_Axisymmetric;
  }
  
	if (CSensitivity != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete CSensitivity[iMarker];
    delete [] CSensitivity;
  }
  
}

void CAdjTNE2EulerSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new double [nBufferR_Vector];
      Buffer_Send_U = new double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
#ifndef NO_MPI
      
      //      /*--- Send/Receive using non-blocking communications ---*/
      //      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
      //      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
      //      send_request.Wait(status);
      //      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
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

void CAdjTNE2EulerSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new double [nBufferR_Vector];
      Buffer_Send_U = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
      }
      
#ifndef NO_MPI
      
      //      /*--- Send/Receive using non-blocking communications ---*/
      //      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
      //      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
      //      send_request.Wait(status);
      //      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
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
          node[iPoint]->SetSolution_Old(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
	}
}

void CAdjTNE2EulerSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new double [nBufferR_Vector];
      Buffer_Send_Limit = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter(iVar);
      }
      
#ifndef NO_MPI
      
      //      /*--- Send/Receive using non-blocking communications ---*/
      //      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Limit, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
      //      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Limit, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
      //      send_request.Wait(status);
      //      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Limit, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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

void CAdjTNE2EulerSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
  double **Gradient = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new double[nDim];
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new double [nBufferR_Vector];
      Buffer_Send_Gradient = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
      }
      
#ifndef NO_MPI
      
      //      /*--- Send/Receive using non-blocking communications ---*/
      //      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Gradient, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
      //      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Gradient, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
      //      send_request.Wait(status);
      //      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Gradient, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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


void CAdjTNE2EulerSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Send_Undivided_Laplacian = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Undivided_Laplacian = new double [nBufferR_Vector];
      Buffer_Send_Undivided_Laplacian = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Undivided_Laplacian[iVar*nVertexS+iVertex] = node[iPoint]->GetUndivided_Laplacian(iVar);
      }
      
#ifndef NO_MPI
      
      //      /*--- Send/Receive using non-blocking communications ---*/
      //      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
      //      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
      //      send_request.Wait(status);
      //      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex] = Buffer_Send_Undivided_Laplacian[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Undivided_Laplacian;
      
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
          Solution[iVar] = Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetUndivided_Laplacian(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Undivided_Laplacian;
      
    }
    
	}
  
}

void CAdjTNE2EulerSolver::SetForceProj_Vector(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	double *ForceProj_Vector, x = 0.0, y = 0.0, z = 0.0, *Normal, C_d, C_l, C_t, C_q;
	double x_origin, y_origin, z_origin, WDrag, Area;
	unsigned short iMarker, iDim;
	unsigned long iVertex, iPoint;
	double Alpha      = (config->GetAoA()*PI_NUMBER)/180.0;
	double Beta       = (config->GetAoS()*PI_NUMBER)/180.0;
	double RefAreaCoeff    = config->GetRefAreaCoeff();
	double RefLengthMoment  = config->GetRefLengthMoment();
	double *RefOriginMoment = config->GetRefOriginMoment();
	double RefVel2, RefDensity;
  
  bool grid_movement = config->GetGrid_Movement();
  
	ForceProj_Vector = new double[nDim];
  
	/*--- If we have a rotating frame problem or an unsteady problem with
   mesh motion, use special reference values for the force coefficients.
   Otherwise, use the freestream values, which is the standard convention. ---*/
	if (grid_movement) {
    double Gas_Constant = config->GetGas_ConstantND();
    double Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    double Mach_Motion = config->GetMach_Motion();
		RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  } else {
		double *Velocity_Inf = config->GetVelocity_FreeStreamND();
		RefVel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
	}
  
  RefDensity  = config->GetDensity_FreeStreamND();
  
	/*--- In parallel computations the Cd, and Cl must be recomputed using all the processors ---*/
#ifdef NO_MPI
	C_d = solver_container[FLOW_SOL]->GetTotal_CDrag();
	C_l = solver_container[FLOW_SOL]->GetTotal_CLift();
	C_t = solver_container[FLOW_SOL]->GetTotal_CT();
	C_q = solver_container[FLOW_SOL]->GetTotal_CQ();
#else
	double *sbuf_force = new double[4];
	double *rbuf_force = new double[4];
	sbuf_force[0] = solver_container[FLOW_SOL]->GetTotal_CDrag();
	sbuf_force[1] = solver_container[FLOW_SOL]->GetTotal_CLift();
	sbuf_force[2] = solver_container[FLOW_SOL]->GetTotal_CT();
	sbuf_force[3] = solver_container[FLOW_SOL]->GetTotal_CQ();
	MPI::COMM_WORLD.Reduce(sbuf_force, rbuf_force, 4, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
	MPI::COMM_WORLD.Bcast(rbuf_force, 4, MPI::DOUBLE, MASTER_NODE);
	C_d = rbuf_force[0];
	C_l = rbuf_force[1];
	C_t = rbuf_force[2];
	C_q = rbuf_force[3];
	delete [] sbuf_force;
	delete [] rbuf_force;
#endif
  
	/*--- Compute coefficients needed for objective function evaluation. ---*/
	C_d += config->GetCteViscDrag();
	double C_p    = 1.0/(0.5*RefDensity*RefAreaCoeff*RefVel2);
	double invCD  = 1.0 / C_d;
	double CLCD2  = C_l / (C_d*C_d);
	double invCQ  = 1.0/C_q;
	double CTRCQ2 = C_t/(RefLengthMoment*C_q*C_q);
  
	x_origin = RefOriginMoment[0]; y_origin = RefOriginMoment[1]; z_origin = RefOriginMoment[2];
  
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		if ((config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE) &&
        (config->GetMarker_All_Monitoring(iMarker) == YES))
			for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
				x = geometry->node[iPoint]->GetCoord(0);
				y = geometry->node[iPoint]->GetCoord(1);
				if (nDim == 3) z = geometry->node[iPoint]->GetCoord(2);
        
				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				switch (config->GetKind_ObjFunc()) {
          case DRAG_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = C_p*cos(Alpha); ForceProj_Vector[1] = C_p*sin(Alpha); }
            if (nDim == 3) { ForceProj_Vector[0] = C_p*cos(Alpha)*cos(Beta); ForceProj_Vector[1] = C_p*sin(Beta); ForceProj_Vector[2] = C_p*sin(Alpha)*cos(Beta); }
            break;
          case LIFT_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = -C_p*sin(Alpha); ForceProj_Vector[1] = C_p*cos(Alpha); }
            if (nDim == 3) { ForceProj_Vector[0] = -C_p*sin(Alpha); ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = C_p*cos(Alpha); }
            break;
          case SIDEFORCE_COEFFICIENT :
            if (nDim == 2) { cout << "This functional is not possible in 2D!!" << endl;
              cout << "Press any key to exit..." << endl; cin.get(); exit(1);
            }
            if (nDim == 3) { ForceProj_Vector[0] = -C_p*sin(Beta) * cos(Alpha); ForceProj_Vector[1] = C_p*cos(Beta); ForceProj_Vector[2] = -C_p*sin(Beta) * sin(Alpha); }
            break;
          case PRESSURE_COEFFICIENT :
            if (nDim == 2) {
              Area = sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1]);
              ForceProj_Vector[0] = -C_p*Normal[0]/Area; ForceProj_Vector[1] = -C_p*Normal[1]/Area;
            }
            if (nDim == 3) {
              Area = sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);
              ForceProj_Vector[0] = -C_p*Normal[0]/Area; ForceProj_Vector[1] = -C_p*Normal[1]/Area; ForceProj_Vector[2] = -C_p*Normal[2]/Area;
            }
            break;
          case MOMENT_X_COEFFICIENT :
            if (nDim == 2) { cout << "This functional is not possible in 2D!!" << endl; cout << "Press any key to exit..." << endl; cin.get(); exit(1); }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = -C_p*(z - z_origin)/RefLengthMoment; ForceProj_Vector[2] = C_p*(y - y_origin)/RefLengthMoment; }
            break;
          case MOMENT_Y_COEFFICIENT :
            if (nDim == 2) { cout << "This functional is not possible in 2D!!" << endl; cout << "Press any key to exit..." << endl; cin.get(); exit(1); }
            if (nDim == 3) { ForceProj_Vector[0] = C_p*(z - z_origin)/RefLengthMoment; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = -C_p*(x - x_origin)/RefLengthMoment; }
            break;
          case MOMENT_Z_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = -C_p*(y - y_origin)/RefLengthMoment; ForceProj_Vector[1] = C_p*(x - x_origin)/RefLengthMoment; }
            if (nDim == 3) { ForceProj_Vector[0] = -C_p*(y - y_origin)/RefLengthMoment; ForceProj_Vector[1] = C_p*(x - x_origin)/RefLengthMoment; ForceProj_Vector[2] = 0; }
            break;
          case EFFICIENCY :
            if (nDim == 2) { ForceProj_Vector[0] = -C_p*(invCD*sin(Alpha)+CLCD2*cos(Alpha)); ForceProj_Vector[1] = C_p*(invCD*cos(Alpha)-CLCD2*sin(Alpha)); }
            if (nDim == 3) { ForceProj_Vector[0] = -C_p*(invCD*sin(Alpha)+CLCD2*cos(Alpha)*cos(Beta)); ForceProj_Vector[1] = -C_p*CLCD2*sin(Beta); ForceProj_Vector[2] = C_p*(invCD*cos(Alpha)-CLCD2*sin(Alpha)*cos(Beta)); }
            break;
          case EQUIVALENT_AREA :
            WDrag = config->GetWeightCd();
            if (nDim == 2) { ForceProj_Vector[0] = C_p*cos(Alpha)*WDrag; ForceProj_Vector[1] = C_p*sin(Alpha)*WDrag; }
            if (nDim == 3) { ForceProj_Vector[0] = C_p*cos(Alpha)*cos(Beta)*WDrag; ForceProj_Vector[1] = C_p*sin(Beta)*WDrag; ForceProj_Vector[2] = C_p*sin(Alpha)*cos(Beta)*WDrag; }
            break;
          case NEARFIELD_PRESSURE :
            WDrag = config->GetWeightCd();
            if (nDim == 2) { ForceProj_Vector[0] = C_p*cos(Alpha)*WDrag; ForceProj_Vector[1] = C_p*sin(Alpha)*WDrag; }
            if (nDim == 3) { ForceProj_Vector[0] = C_p*cos(Alpha)*cos(Beta)*WDrag; ForceProj_Vector[1] = C_p*sin(Beta)*WDrag; ForceProj_Vector[2] = C_p*sin(Alpha)*cos(Beta)*WDrag; }
            break;
          case FORCE_X_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = C_p; ForceProj_Vector[1] = 0.0; }
            if (nDim == 3) { ForceProj_Vector[0] = C_p; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 0.0; }
            break;
          case FORCE_Y_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = C_p; }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = C_p; ForceProj_Vector[2] = 0.0; }
            break;
          case FORCE_Z_COEFFICIENT :
            if (nDim == 2) {cout << "This functional is not possible in 2D!!" << endl;
              cout << "Press any key to exit..." << endl;
              cin.get(); exit(1);
            }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = C_p; }
            break;
          case THRUST_COEFFICIENT :
            if (nDim == 2) {cout << "This functional is not possible in 2D!!" << endl;
              cout << "Press any key to exit..." << endl;
              cin.get(); exit(1);
            }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = C_p; }
            break;
          case TORQUE_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = C_p*(y - y_origin)/RefLengthMoment; ForceProj_Vector[1] = -C_p*(x - x_origin)/RefLengthMoment; }
            if (nDim == 3) { ForceProj_Vector[0] = C_p*(y - y_origin)/RefLengthMoment; ForceProj_Vector[1] = -C_p*(x - x_origin)/RefLengthMoment; ForceProj_Vector[2] = 0; }
            break;
          case FIGURE_OF_MERIT :
            if (nDim == 2) {cout << "This functional is not possible in 2D!!" << endl;
              cout << "Press any key to exit..." << endl;
              cin.get(); exit(1);
            }
            if (nDim == 3) {
              ForceProj_Vector[0] = -C_p*invCQ;
              ForceProj_Vector[1] = -C_p*CTRCQ2*(z - z_origin);
              ForceProj_Vector[2] =  C_p*CTRCQ2*(y - y_origin);
            }
            break;
          case FREE_SURFACE :
            if (nDim == 2) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 0.0; }
            break;
          case NOISE:
            if (nDim == 2) { ForceProj_Vector[0] = 0.0;
              ForceProj_Vector[1] = 0.0; }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0;
              ForceProj_Vector[1] = 0.0;
              ForceProj_Vector[2] = 0.0; }
            break;
          case HEAT_LOAD:
            if (nDim == 2) { ForceProj_Vector[0] = 0.0;
              ForceProj_Vector[1] = 0.0; }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0;
              ForceProj_Vector[1] = 0.0;
              ForceProj_Vector[2] = 0.0; }
            break;
				}
        
				/*--- Store the force projection vector at this node ---*/
				node[iPoint]->SetForceProj_Vector(ForceProj_Vector);
        
			}
  
	delete [] ForceProj_Vector;
}


void CAdjTNE2EulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
	unsigned long iPoint, Point_Fine;
	unsigned short iMesh, iChildren, iVar;
	double LevelSet, Area_Children, Area_Parent, LevelSet_Fine, *Solution, *Solution_Fine;
  
	bool restart = config->GetRestart();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- If restart solution, then interpolate the flow solution to
   all the multigrid levels, this is important with the dual time strategy ---*/
  if (restart) {
    Solution = new double[nVar];
    for (iMesh = 1; iMesh <= config->GetMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
        for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
        for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
          Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
          Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
          Solution_Fine = solver_container[iMesh-1][ADJFLOW_SOL]->node[Point_Fine]->GetSolution();
          for (iVar = 0; iVar < nVar; iVar++) {
            Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
          }
        }
        solver_container[iMesh][ADJFLOW_SOL]->node[iPoint]->SetSolution(Solution);
        
      }
      solver_container[iMesh][ADJFLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    }
    delete [] Solution;
  }
  
  /*--- The value of the solution for the first iteration of the dual time ---*/
  for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      if ((ExtIter == 0) && (dual_time)) {
        solver_container[iMesh][ADJFLOW_SOL]->node[iPoint]->Set_Solution_time_n();
        solver_container[iMesh][ADJFLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
      }
    }
  }
  
}

void CAdjTNE2EulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
  
  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the flow problem may use different methods). ---*/
  bool implicit       = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool upwind_2nd     = (config->GetKind_Upwind_AdjFlow() == ROE_2ND);
  bool center         = (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst     = (config->GetKind_Centered_AdjFlow() == JST);
  bool limiter        = (config->GetKind_SlopeLimit() != NONE);
  double adj_limit    = config->GetAdjointLimit();
  
  /*--- Compute nacelle inflow and exhaust properties ---*/
  GetNacelle_Properties(geometry, config, iMesh);
  
	/*--- Residual initialization ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Set the primitive variables incompressible and compressible
     adjoint variables ---*/
    //		node[iPoint]->SetPrimVar_Compressible(adj_limit);
    
		/*--- Initialize the convective residual vector ---*/
		LinSysRes.SetBlock_Zero(iPoint);
    
	}
  
  /*--- Upwind second order reconstruction ---*/
  if ((upwind_2nd) && (iMesh == MESH_0)) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    
    /*--- Limiter computation ---*/
		if (limiter) SetSolution_Limiter(geometry, config);
	}
  
  /*--- Artificial dissipation ---*/
  if (center) {
    if ((center_jst) && (iMesh == MESH_0)) {
      SetDissipation_Switch(geometry, config);
      SetUndivided_Laplacian(geometry, config);
      if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
      if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    }
  }
  
	/*--- Implicit solution ---*/
	if (implicit) Jacobian.SetValZero();
  
}

void CAdjTNE2EulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                            CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
	unsigned long iEdge, iPoint, jPoint;
  
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool high_order_diss = ((config->GetKind_Centered_AdjFlow() == JST) && (iMesh == MESH_0));
	bool grid_movement  = config->GetGrid_Movement();
  
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Points in edge, normal, and neighbors---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
		numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());
    
		/*--- Adjoint variables w/o reconstruction ---*/
		numerics->SetAdjointVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    
		/*--- Conservative variables w/o reconstruction ---*/
		numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                              solver_container[FLOW_SOL]->node[jPoint]->GetSolution());
    
    numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(),
                            solver_container[FLOW_SOL]->node[jPoint]->GetSoundSpeed());
    numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(),
                          solver_container[FLOW_SOL]->node[jPoint]->GetEnthalpy());
    
		numerics->SetLambda(solver_container[FLOW_SOL]->node[iPoint]->GetLambda(),
                        solver_container[FLOW_SOL]->node[jPoint]->GetLambda());
    
		if (high_order_diss) {
			numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(), node[jPoint]->GetUndivided_Laplacian());
			numerics->SetSensor(solver_container[FLOW_SOL]->node[iPoint]->GetSensor(),
                          solver_container[FLOW_SOL]->node[jPoint]->GetSensor());
		}
    
		/*--- Mesh motion ---*/
		if (grid_movement) {
			numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
		}
    
		/*--- Compute residuals ---*/
		numerics->ComputeResidual(Res_Conv_i, Res_Visc_i, Res_Conv_j, Res_Visc_j,
                              Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
    
		/*--- Update convective and artificial dissipation residuals ---*/
		LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
		LinSysRes.SubtractBlock(jPoint, Res_Conv_j);
    LinSysRes.SubtractBlock(iPoint, Res_Visc_i);
    LinSysRes.SubtractBlock(jPoint, Res_Visc_j);
    
		/*--- Implicit contribution to the residual ---*/
		if (implicit) {
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
			Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
			Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ji);
			Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_jj);
		}
	}
}


void CAdjTNE2EulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {
	double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, *Limiter_i = NULL,
  *Limiter_j = NULL, *Psi_i = NULL, *Psi_j = NULL, *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
  
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool high_order_diss = (((config->GetKind_Upwind_AdjFlow() == ROE_2ND) ||
                           (config->GetKind_Upwind_AdjFlow() == SW_2ND)) && (iMesh == MESH_0));
	bool grid_movement  = config->GetGrid_Movement();
	bool limiter = (config->GetKind_SlopeLimit() != NONE);
  
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Adjoint variables w/o reconstruction ---*/
    Psi_i = node[iPoint]->GetSolution(); Psi_j = node[jPoint]->GetSolution();
    numerics->SetAdjointVar(Psi_i, Psi_j);
    
		/*--- Conservative variables w/o reconstruction ---*/
		U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
		U_j = solver_container[FLOW_SOL]->node[jPoint]->GetSolution();
		numerics->SetConservative(U_i, U_j);
    
    numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(),
                            solver_container[FLOW_SOL]->node[jPoint]->GetSoundSpeed());
    numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(),
                          solver_container[FLOW_SOL]->node[jPoint]->GetEnthalpy());
    
		/*--- Mesh motion ---*/
		if (grid_movement) {
			numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
		}
    
		/*--- High order reconstruction using MUSCL strategy ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
      Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
    }
    
    Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
    if (limiter) { Limiter_i = node[iPoint]->GetLimiter(); Limiter_j = node[jPoint]->GetLimiter(); }
    
    for (iVar = 0; iVar < nVar; iVar++) {
      Project_Grad_i = 0; Project_Grad_j = 0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
        Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
      }
      if (limiter) {
        Solution_i[iVar] = Psi_i[iVar] + Project_Grad_i*Limiter_i[iDim];
        Solution_j[iVar] = Psi_j[iVar] + Project_Grad_j*Limiter_j[iDim];
      }
      else {
        Solution_i[iVar] = Psi_i[iVar] + Project_Grad_i;
        Solution_j[iVar] = Psi_j[iVar] + Project_Grad_j;
        
      }
    }
    /*--- Set conservative variables with reconstruction ---*/
    numerics->SetAdjointVar(Solution_i, Solution_j);
    
    
		/*--- Compute the residual---*/
    numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii,
                              Jacobian_ij, Jacobian_ji, Jacobian_jj,
                              config);
    
		/*--- Add and Subtract Residual ---*/
    LinSysRes.SubtractBlock(iPoint, Residual_i);
    LinSysRes.SubtractBlock(jPoint, Residual_j);
    
    /*--- Implicit contribution to the residual ---*/
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ji);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_jj);
    }
	}
}

void CAdjTNE2EulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                          CConfig *config, unsigned short iMesh) {
  
	unsigned short iVar, iDim;
	unsigned long iPoint;
  
	bool axisymmetric = config->GetAxisymmetric();
  
	for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
  
	if (axisymmetric) {
    
		bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
    /*--- Zero out Jacobian structure ---*/
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar ++)
        for (unsigned short jVar = 0; jVar < nVar; jVar ++)
          Jacobian_i[iVar][jVar] = 0.0;
    }
    
		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
			
			/*--- Set solution ---*/
			numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(), solver_container[FLOW_SOL]->node[iPoint]->GetSolution());
      
			/*--- Set adjoint variables ---*/
			numerics->SetAdjointVar(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
      
			/*--- Set control volume ---*/
			numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
			/*--- Set coordinate ---*/
			numerics->SetCoord(geometry->node[iPoint]->GetCoord(),geometry->node[iPoint]->GetCoord());
      
			/*--- Compute Source term Residual ---*/
			numerics->ComputeResidual(Residual, Jacobian_i, config);
      
			/*--- Add Residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
		}
	}
  
}

void CAdjTNE2EulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint, iEdge;
	unsigned short iVar;
	double *Diff;
  
  Diff = new double[nVar];
  
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		node[iPoint]->SetUnd_LaplZero();
  
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
    
		for (iVar = 0; iVar < nVar; iVar++)
			Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);
    
#ifdef STRUCTURED_GRID
    
		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    
#else
    
    bool boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    bool boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both in the boundary ---*/
		if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
		}
		
		/*--- iPoint inside the domain, jPoint on the boundary ---*/
		if (!boundary_i && boundary_j)
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
		
		/*--- jPoint inside the domain, iPoint on the boundary ---*/
		if (boundary_i && !boundary_j)
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    
#endif
    
	}
  
#ifdef STRUCTURED_GRID
  
  unsigned long Point_Normal = 0, iVertex;
	unsigned short iMarker;
  double *Psi_mirror;
  
  Psi_mirror = new double[nVar];
  
	/*--- Loop over all boundaries and include an extra contribution
	 from a halo node. Find the nearest normal, interior point
	 for a boundary node and make a linear approximation. ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
				config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
				config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY &&
				config->GetMarker_All_Boundary(iMarker) != PERIODIC_BOUNDARY) {
      
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
				if (geometry->node[iPoint]->GetDomain()) {
          
					Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
          
					/*--- Interpolate & compute difference in the conserved variables ---*/
					for (iVar = 0; iVar < nVar; iVar++) {
						Psi_mirror[iVar] = 2.0*node[iPoint]->GetSolution(iVar) - node[Point_Normal]->GetSolution(iVar);
						Diff[iVar]   = node[iPoint]->GetSolution(iVar) - Psi_mirror[iVar];
					}
          
					/*--- Subtract contribution at the boundary node only ---*/
					node[iPoint]->SubtractUnd_Lapl(Diff);
				}
			}
    }
  }
  
	delete [] Psi_mirror;
  
#endif
  
  delete [] Diff;
  
  /*--- MPI parallelization ---*/
  Set_MPI_Undivided_Laplacian(geometry, config);
  
}

void CAdjTNE2EulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta, *local_Res_TruncError, Vol;
  
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
		Delta = Vol / solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time();
    
		Jacobian.AddVal2Diag(iPoint, Delta);
    
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			LinSysRes[total_index] = -(LinSysRes[total_index] + local_Res_TruncError[iVar]);
			LinSysSol[total_index] = 0.0;
			AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex());
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
  
	/*--- Solve the linear system (Krylov subspace methods) ---*/
  CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
  
  CPreconditioner* precond = NULL;
  if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
    Jacobian.BuildJacobiPreconditioner();
    precond = new CJacobiPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LU_SGS) {
    precond = new CLU_SGSPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    Jacobian.BuildJacobiPreconditioner();
    Jacobian.BuildLineletPreconditioner(geometry, config);
    precond = new CLineletPreconditioner(Jacobian, geometry, config);
  }
  
  CSysSolve system;
  if (config->GetKind_Linear_Solver() == BCGSTAB)
    system.BCGSTAB(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                   config->GetLinear_Solver_Iter(), false);
  else if (config->GetKind_Linear_Solver() == FGMRES)
    system.FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                  config->GetLinear_Solver_Iter(), false);
  
  delete mat_vec;
  delete precond;
  
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		for (iVar = 0; iVar < nVar; iVar++)
			node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar]);
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CAdjTNE2EulerSolver::Solve_LinearSystem(CGeometry *geometry, CSolver **solver_container, CConfig *config){
	unsigned long iPoint;
	unsigned long total_index;
	unsigned short iVar;
	double *ObjFuncSource;
  
	/*--- Build linear system ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		ObjFuncSource = node[iPoint]->GetObjFuncSource();
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			LinSysRes[total_index] = ObjFuncSource[iVar];
			LinSysSol[total_index] = 0.0;
		}
	}
  
	/*--- Solve the linear system (Krylov subspace methods) ---*/
  CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
  
  CPreconditioner* precond = NULL;
  if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
    Jacobian.BuildJacobiPreconditioner();
    precond = new CJacobiPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LU_SGS) {
    precond = new CLU_SGSPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    Jacobian.BuildJacobiPreconditioner();
    Jacobian.BuildLineletPreconditioner(geometry, config);
    precond = new CLineletPreconditioner(Jacobian, geometry, config);
  }
  
  CSysSolve system;
  if (config->GetKind_Linear_Solver() == BCGSTAB)
    system.BCGSTAB(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                   config->GetLinear_Solver_Iter(), true);
  else if (config->GetKind_Linear_Solver() == FGMRES)
    system.FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                  config->GetLinear_Solver_Iter(), true);
  
  delete mat_vec;
  delete precond;
  
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		for (iVar = 0; iVar < nVar; iVar++)
			node[iPoint]->SetSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar]);
  
}

void CAdjTNE2EulerSolver::Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
	unsigned long iVertex, iPoint, Neigh;
	unsigned short iPos, jPos;
	unsigned short iDim, iMarker, iNeigh;
	double *d = NULL, *Normal = NULL, *Psi = NULL, *U = NULL, Enthalpy, conspsi, Mach_Inf,
  Area, **PrimVar_Grad = NULL, **ConsVar_Grad = NULL, *ConsPsi_Grad = NULL,
  ConsPsi, d_press, grad_v, Beta2, v_gradconspsi, UnitaryNormal[3], *RotVel = NULL, *GridVel = NULL,
  LevelSet, Target_LevelSet;
  
	bool grid_movement  = config->GetGrid_Movement();
  
	/*--- Initialize sensitivities to zero ---*/
	Total_Sens_Geo = 0.0; Total_Sens_Mach = 0.0; Total_Sens_AoA = 0.0;
	Total_Sens_Press = 0.0; Total_Sens_Temp = 0.0;
	//	Total_Sens_Far = 0.0;
  
	/*--- Compute surface sensitivity ---*/
  
  /*--- Loop over boundary markers to select those for Euler walls ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    if (config->GetMarker_All_Boundary(iMarker) == EULER_WALL)
      
    /*--- Loop over points on the surface to store the auxiliary variable ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Psi = node[iPoint]->GetSolution();
          U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
          Enthalpy = solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy();
          conspsi = U[0]*Psi[0] + U[0]*Enthalpy*Psi[nDim+1];
          for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];
          
          node[iPoint]->SetAuxVar(conspsi);
          
          /*--- Also load the auxiliary variable for first neighbors ---*/
          for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
            Neigh = geometry->node[iPoint]->GetPoint(iNeigh);
            Psi = node[Neigh]->GetSolution();
            U = solver_container[FLOW_SOL]->node[Neigh]->GetSolution();
            Enthalpy = solver_container[FLOW_SOL]->node[Neigh]->GetEnthalpy();
            conspsi = U[0]*Psi[0] + U[0]*Enthalpy*Psi[nDim+1];
            for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];
            node[Neigh]->SetAuxVar(conspsi);
          }
        }
      }
  
  /*--- Compute surface gradients of the auxiliary variable ---*/
  SetAuxVar_Surface_Gradient(geometry, config);
  
  /*--- Evaluate the shape sensitivity ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Sens_Geo[iMarker] = 0.0;
    
    if (config->GetMarker_All_Boundary(iMarker) == EULER_WALL) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          
          d = node[iPoint]->GetForceProj_Vector();
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
          
          PrimVar_Grad = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
          ConsVar_Grad = solver_container[FLOW_SOL]->node[iPoint]->GetGradient();
          ConsPsi_Grad = node[iPoint]->GetAuxVarGradient();
          ConsPsi = node[iPoint]->GetAuxVar();
          
          /*--- Adjustment for grid movement - double check this ---*/
          if (grid_movement) GridVel = geometry->node[iPoint]->GetGridVel();
          
          d_press = 0.0; grad_v = 0.0; v_gradconspsi = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            
            /*-- Retrieve the value of the pressure gradient ---*/
            d_press += d[iDim]*PrimVar_Grad[nDim+1][iDim];
            
            /*-- Retrieve the value of the velocity gradient ---*/
            grad_v += PrimVar_Grad[iDim+1][iDim]*ConsPsi;
            
            /*-- Retrieve the value of the theta gradient ---*/
            v_gradconspsi += solver_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim, false) * ConsPsi_Grad[iDim];
            if (grid_movement) v_gradconspsi -= GridVel[iDim] * ConsPsi_Grad[iDim];
          }
          
          /*--- Compute sensitivity for each surface point ---*/
          CSensitivity[iMarker][iVertex] = (d_press + grad_v + v_gradconspsi) * Area;
          Sens_Geo[iMarker] -= CSensitivity[iMarker][iVertex] * Area;
        }
      }
      Total_Sens_Geo += Sens_Geo[iMarker];
    }
  }
  
}

void CAdjTNE2EulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	double *d = NULL, *Normal, *U, *Psi_Aux, ProjVel = 0.0, bcn, vn = 0.0, Area, *UnitaryNormal, *Coord;
	double *Velocity, *Psi, *ObjFuncSource, Enthalpy = 0.0, sq_vel, phin, phis1, phis2, DensityInc = 0.0, BetaInc2 = 0.0;
	unsigned short iDim, iVar, jDim;
	double *dPressure;
  
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool grid_movement = config->GetGrid_Movement();
  
	UnitaryNormal = new double[nDim];
	Velocity = new double[nDim];
	Psi      = new double[nVar];
	ObjFuncSource = new double[nVar];
	dPressure = new double[nVar];
  
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		if (geometry->node[iPoint]->GetDomain()) {
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Coord = geometry->node[iPoint]->GetCoord();
      
      /*--- Create a copy of the adjoint solution ---*/
      Psi_Aux = node[iPoint]->GetSolution();
      for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];
      
			/*--- Flow solution ---*/
			U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
      
      /*--- Read the value of the objective function ---*/
      d = node[iPoint]->GetForceProj_Vector();
      
      /*--- Normal vector computation ---*/
			Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
			Area = sqrt(Area);
			for (iDim = 0; iDim < nDim; iDim++) UnitaryNormal[iDim] = -Normal[iDim]/Area;
      
      for (iDim = 0; iDim < nDim; iDim++)
        Velocity[iDim] = U[iDim+1] / U[0];
      
      Enthalpy = solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy();
      sq_vel   = 0.5*solver_container[FLOW_SOL]->node[iPoint]->GetVelocity2();
      
      /*--- Compute projections ---*/
      ProjVel = 0.0; bcn = 0.0; vn = 0.0, phin = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel -= Velocity[iDim]*Normal[iDim];
        bcn     += d[iDim]*UnitaryNormal[iDim];
        vn      += Velocity[iDim]*UnitaryNormal[iDim];
        phin    += Psi[iDim+1]*UnitaryNormal[iDim];
      }
      
      /*--- Extra boundary term for grid movement ---*/
      if (grid_movement) {
        double ProjGridVel = 0.0;
        double *GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitaryNormal[iDim];
        phin -= Psi[nVar-1]*ProjGridVel;
      }
      
      /*--- Introduce the boundary condition ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Psi[iDim+1] -= ( phin - bcn ) * UnitaryNormal[iDim];
      
      /*--- Inner products after introducing BC (Psi has changed) ---*/
      phis1 = 0.0; phis2 = Psi[0] + Enthalpy * Psi[nVar-1];
      for (iDim = 0; iDim < nDim; iDim++) {
        phis1 -= Normal[iDim]*Psi[iDim+1];
        phis2 += Velocity[iDim]*Psi[iDim+1];
      }
      
      /*--- Flux of the Euler wall ---*/
      Residual[0] = ProjVel * Psi[0] - phis2 * ProjVel + phis1 * Gamma_Minus_One * sq_vel;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual[iDim+1] = ProjVel * Psi[iDim+1] - phis2 * Normal[iDim] - phis1 * Gamma_Minus_One * Velocity[iDim];
      Residual[nVar-1] = ProjVel * Psi[nVar-1] + phis1 * Gamma_Minus_One;
      
      /*--- Flux adjustment for grid movement (TDE) ---*/
      if (grid_movement) {
        double ProjGridVel = 0.0;
        double *GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel -= GridVel[iDim]*Normal[iDim];
        Residual[0] -= ProjGridVel*Psi[0];
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[iDim+1] -= ProjGridVel*Psi[iDim+1];
        Residual[nVar-1] -= ProjGridVel*Psi[nVar-1];
      }
      
      if (implicit) {
        
        /*--- Adjoint density ---*/
        Jacobian_ii[0][0] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_ii[0][iDim+1] = -ProjVel * (Velocity[iDim] - UnitaryNormal[iDim] * vn);
        Jacobian_ii[0][nVar-1] = -ProjVel * Enthalpy;
        
        /*--- Adjoint velocities ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          Jacobian_ii[iDim+1][0] = -Normal[iDim];
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_ii[iDim+1][jDim+1] = -ProjVel*(UnitaryNormal[jDim]*UnitaryNormal[iDim] - Normal[iDim] * (Velocity[jDim] - UnitaryNormal[jDim] * vn));
          Jacobian_ii[iDim+1][iDim+1] += ProjVel;
          Jacobian_ii[iDim+1][nVar-1] = -Normal[iDim] * Enthalpy;
        }
        
        /*--- Adjoint energy ---*/
        Jacobian_ii[nVar-1][0] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_ii[nVar-1][iDim+1] = 0.0;
        Jacobian_ii[nVar-1][nVar-1] = ProjVel;
        
        /*--- Jacobian contribution due to grid movement (TDE) ---*/
        if (grid_movement) {
          double ProjGridVel = 0.0;
          double *GridVel = geometry->node[iPoint]->GetGridVel();
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel -= GridVel[iDim]*Normal[iDim];
          Jacobian_ii[0][0] -= ProjGridVel;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_ii[iDim+1][iDim+1] -= ProjGridVel;
          Jacobian_ii[nVar-1][nVar-1] -= ProjGridVel;
        }
        
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
        
      }
      
      /*--- Update residual ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      
		}
	}
  
	delete [] Velocity;
	delete [] UnitaryNormal;
	delete [] Psi;
  
}

void CAdjTNE2EulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                       CConfig *config, unsigned short val_marker) {
  
	unsigned long iVertex, iPoint;
	double *Normal, *U, *Psi_Aux, ProjVel = 0.0, vn = 0.0, Area, *UnitaryNormal, *Coord;
	double *Velocity, *Psi, Enthalpy = 0.0, sq_vel, phin, phis1, phis2, DensityInc = 0.0, BetaInc2 = 0.0;
	unsigned short iDim, iVar, jDim;
  
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool grid_movement = config->GetGrid_Movement();
  
	UnitaryNormal = new double[nDim];
	Velocity = new double[nDim];
	Psi      = new double[nVar];
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		if (geometry->node[iPoint]->GetDomain()) {
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Coord = geometry->node[iPoint]->GetCoord();
      
			/*--- Create a copy of the adjoint solution ---*/
			Psi_Aux = node[iPoint]->GetSolution();
			for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];
      
			/*--- Flow solution ---*/
			U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
      
			Area = 0;
			for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
			Area = sqrt(Area);
      
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim]   = -Normal[iDim]/Area;
      
      for (iDim = 0; iDim < nDim; iDim++)
        Velocity[iDim] = U[iDim+1] / U[0];
      
      Enthalpy = solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy();
      sq_vel   = 0.5*solver_container[FLOW_SOL]->node[iPoint]->GetVelocity2();
      
      /*--- Compute projections ---*/
      ProjVel = 0.0; vn = 0.0, phin = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel -= Velocity[iDim]*Normal[iDim];
        vn      += Velocity[iDim]*UnitaryNormal[iDim];
        phin    += Psi[iDim+1]*UnitaryNormal[iDim];
      }
      
      /*--- Grid Movement ---*/
      if (grid_movement) {
        double ProjGridVel = 0.0;
        double *GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) {
          ProjGridVel += GridVel[iDim]*UnitaryNormal[iDim];
        }
        phin -= Psi[nVar-1]*ProjGridVel;
      }
      
      /*--- Introduce the boundary condition ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Psi[iDim+1] -= phin * UnitaryNormal[iDim];
      
      /*--- Inner products after introducing BC (Psi has changed) ---*/
      phis1 = 0.0; phis2 = Psi[0] + Enthalpy * Psi[nVar-1];
      for (iDim = 0; iDim < nDim; iDim++) {
        phis1 -= Normal[iDim]*Psi[iDim+1];
        phis2 += Velocity[iDim]*Psi[iDim+1];
      }
      
      /*--- Flux of the Euler wall ---*/
      Residual[0] = ProjVel * Psi[0] - phis2 * ProjVel + phis1 * Gamma_Minus_One * sq_vel;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual[iDim+1] = ProjVel * Psi[iDim+1] - phis2 * Normal[iDim] - phis1 * Gamma_Minus_One * Velocity[iDim];
      Residual[nVar-1] = ProjVel * Psi[nVar-1] + phis1 * Gamma_Minus_One;
      
      /*--- Grid Movement ---*/
      if (grid_movement) {
        double ProjGridVel = 0.0;
        double *GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel -= GridVel[iDim]*Normal[iDim];
        Residual[0] -= ProjGridVel*Psi[0];
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[iDim+1] -= ProjGridVel*Psi[iDim+1];
        Residual[nVar-1] -= ProjGridVel*Psi[nVar-1];
      }
      
			/*--- Update residual ---*/
			LinSysRes.SubtractBlock(iPoint, Residual);
      
      
			/*--- Implicit stuff ---*/
			if (implicit) {
        
        /*--- Adjoint density ---*/
        Jacobian_ii[0][0] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_ii[0][iDim+1] = -ProjVel * (Velocity[iDim] - UnitaryNormal[iDim] * vn);
        Jacobian_ii[0][nVar-1] = -ProjVel * Enthalpy;
        
        /*--- Adjoint velocities ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          Jacobian_ii[iDim+1][0] = -Normal[iDim];
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_ii[iDim+1][jDim+1] = -ProjVel*(UnitaryNormal[jDim]*UnitaryNormal[iDim] - Normal[iDim] * (Velocity[jDim] - UnitaryNormal[jDim] * vn));
          Jacobian_ii[iDim+1][iDim+1] += ProjVel;
          Jacobian_ii[iDim+1][nVar-1] = -Normal[iDim] * Enthalpy;
        }
        
        /*--- Adjoint energy ---*/
        Jacobian_ii[nVar-1][0] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_ii[nVar-1][iDim+1] = 0.0;
        Jacobian_ii[nVar-1][nVar-1] = ProjVel;
        
        /*--- Contribution from grid movement ---*/
        if (grid_movement) {
          double ProjGridVel = 0.0;
          double *GridVel = geometry->node[iPoint]->GetGridVel();
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel -= GridVel[iDim]*Normal[iDim];
          Jacobian_ii[0][0] -= ProjGridVel;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_ii[iDim+1][iDim+1] -= ProjGridVel;
          Jacobian_ii[nVar-1][nVar-1] -= ProjGridVel;
        }
        
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
			}
		}
	}
  
	delete [] Velocity;
	delete [] UnitaryNormal;
	delete [] Psi;
}


void CAdjTNE2EulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                       CConfig *config, unsigned short val_marker) {
  
	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	double *Normal, *U_domain, *U_infty, *Psi_domain, *Psi_infty;
  
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool grid_movement = config->GetGrid_Movement();
  
	Normal = new double[nDim];
	U_domain = new double[nVar]; U_infty = new double[nVar];
	Psi_domain = new double[nVar]; Psi_infty = new double[nVar];
  
	/*--- Loop over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- If the node belongs to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
			conv_numerics->SetNormal(Normal);
      
			/*--- Flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
      
			/*--- Solution at infinity ---*/
      /*--- Flow Solution at infinity ---*/
      U_infty[0] = solver_container[FLOW_SOL]->GetDensity_Inf();
      U_infty[1] = solver_container[FLOW_SOL]->GetDensity_Velocity_Inf(0);
      U_infty[2] = solver_container[FLOW_SOL]->GetDensity_Velocity_Inf(1);
      U_infty[3] = solver_container[FLOW_SOL]->GetDensity_Energy_Inf();
      if (nDim == 3) {
        U_infty[3] = solver_container[FLOW_SOL]->GetDensity_Velocity_Inf(2);
        U_infty[4] = solver_container[FLOW_SOL]->GetDensity_Energy_Inf();
      }
			conv_numerics->SetConservative(U_domain, U_infty);
      
      /*--- Adjoint flow solution at the wall ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
        Psi_infty[iVar] = 0.0;
      }
      conv_numerics->SetAdjointVar(Psi_domain, Psi_infty);
      
      conv_numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(),
                                   solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed());
      conv_numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(),
                                 solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy());
      
			/*--- Grid Movement ---*/
			if (grid_movement)
				conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
			/*--- Compute the upwind flux ---*/
      conv_numerics->ComputeResidual(Residual_i, Residual_j,
                                     Jacobian_ii, Jacobian_ij,
                                     Jacobian_ji, Jacobian_jj, config);
      
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      
      /*--- Implicit contribution to the residual ---*/
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      
    }
  }
	delete [] Normal;
	delete [] U_domain; delete [] U_infty;
	delete [] Psi_domain; delete [] Psi_infty;
}

CAdjTNE2NSSolver::CAdjTNE2NSSolver(void) : CAdjTNE2EulerSolver() { }

CAdjTNE2NSSolver::CAdjTNE2NSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CAdjTNE2EulerSolver() {
	unsigned long iPoint, index, iVertex;
	string text_line, mesh_filename;
	unsigned short iDim, iVar, iMarker;
	ifstream restart_file;
	string filename, AdjExt;
  
	bool restart = config->GetRestart();
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
	/*--- Set the gamma value ---*/
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	/*--- Define geometry constants in the solver structure ---*/
	nDim    = geometry->GetnDim();
  nMarker = config->GetnMarker_All();
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
	nVar = nDim + 2;
	node = new CVariable*[nPoint];
  
	/*--- Define some auxiliary arrays related to the residual ---*/
	Residual      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
	Residual_RMS  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
	Residual_Max  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
	Point_Max  = new unsigned long[nVar]; for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
	Residual_i    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
	Residual_j    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
	Res_Conv_i = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Res_Conv_i[iVar]    = 0.0;
  Res_Visc_i   = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Res_Visc_i[iVar]    = 0.0;
	Res_Conv_j = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Res_Conv_j[iVar]    = 0.0;
  Res_Visc_j   = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Res_Visc_j[iVar]    = 0.0;
  
	/*--- Define some auxiliary arrays related to the solution ---*/
	Solution   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
	Solution_i = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
	Solution_j = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;
  
	/*--- Define some auxiliary vectors related to the geometry ---*/
	Vector   = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
	Vector_i = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
	Vector_j = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Point to point Jacobians. These are always defined because
   they are also used for sensitivity calculations. ---*/
  Jacobian_i = new double* [nVar];
  Jacobian_j = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new double [nVar];
    Jacobian_j[iVar] = new double [nVar];
  }
  
  /*--- Solution and residual vectors ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT) {
    
		Jacobian_ii = new double* [nVar];
		Jacobian_ij = new double* [nVar];
		Jacobian_ji = new double* [nVar];
		Jacobian_jj = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_ii[iVar] = new double [nVar];
			Jacobian_ij[iVar] = new double [nVar];
			Jacobian_ji[iVar] = new double [nVar];
			Jacobian_jj[iVar] = new double [nVar];
		}
    if (rank == MASTER_NODE)
      cout << "Initialize jacobian structure (Adjoint N-S). MG level: " << iMesh <<"." << endl;
		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
    
  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No jacobian structure (Adjoint N-S). MG level: " << iMesh <<"." << endl;
  }
  
	/*--- Array structures for computation of gradients by least squares ---*/
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			cvector[iVar] = new double [nDim];
	}
  
	/*--- Sensitivity definition and coefficient on all markers ---*/
	CSensitivity = new double* [nMarker];
	for (iMarker=0; iMarker<nMarker; iMarker++) {
		CSensitivity[iMarker] = new double [geometry->nVertex[iMarker]];
	}
	Sens_Geo   = new double[nMarker];
	Sens_Mach  = new double[nMarker];
	Sens_AoA   = new double[nMarker];
	Sens_Press = new double[nMarker];
	Sens_Temp  = new double[nMarker];
  
  /*--- Initialize sensitivities to zero ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		Sens_Geo[iMarker]   = 0.0;
		Sens_Mach[iMarker]  = 0.0;
		Sens_AoA[iMarker]   = 0.0;
		Sens_Press[iMarker] = 0.0;
		Sens_Temp[iMarker]  = 0.0;
		for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
			CSensitivity[iMarker][iVertex] = 0.0;
	}
  
	/*--- Initialize the adjoint variables to zero (infinity state) ---*/
	PsiRho_Inf = 0.0;
  if (config->GetKind_ObjFunc() == HEAT_LOAD)
    PsiE_Inf = -1.0;
  else
    PsiE_Inf = 0.0;
	Phi_Inf = new double [nDim];
	Phi_Inf[0] = 0.0; Phi_Inf[1] = 0.0;
	if (nDim == 3) Phi_Inf[2] = 0.0;
  
	if (!restart || geometry->GetFinestMGLevel() == false) {
		/*--- Restart the solution from infinity ---*/
		for (iPoint = 0; iPoint < nPoint; iPoint++)
			node[iPoint] = new CAdjNSVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, nDim, nVar, config);
	}
	else {
    
		/*--- Restart the solution from file information ---*/
		mesh_filename = config->GetSolution_AdjFileName();
    
		/*--- Change the name, depending of the objective function ---*/
		filename.assign(mesh_filename);
		filename.erase (filename.end()-4, filename.end());
		switch (config->GetKind_ObjFunc()) {
      case DRAG_COEFFICIENT: AdjExt = "_cd.dat"; break;
      case LIFT_COEFFICIENT: AdjExt = "_cl.dat"; break;
      case SIDEFORCE_COEFFICIENT: AdjExt = "_csf.dat"; break;
      case PRESSURE_COEFFICIENT: AdjExt = "_cp.dat"; break;
      case MOMENT_X_COEFFICIENT: AdjExt = "_cmx.dat"; break;
      case MOMENT_Y_COEFFICIENT: AdjExt = "_cmy.dat"; break;
      case MOMENT_Z_COEFFICIENT: AdjExt = "_cmz.dat"; break;
      case EFFICIENCY: AdjExt = "_eff.dat"; break;
      case EQUIVALENT_AREA: AdjExt = "_ea.dat"; break;
      case NEARFIELD_PRESSURE: AdjExt = "_nfp.dat"; break;
      case FORCE_X_COEFFICIENT: AdjExt = "_cfx.dat"; break;
      case FORCE_Y_COEFFICIENT: AdjExt = "_cfy.dat"; break;
      case FORCE_Z_COEFFICIENT: AdjExt = "_cfz.dat"; break;
      case THRUST_COEFFICIENT: AdjExt = "_ct.dat"; break;
      case TORQUE_COEFFICIENT: AdjExt = "_cq.dat"; break;
      case FIGURE_OF_MERIT: AdjExt = "_merit.dat"; break;
      case FREE_SURFACE: AdjExt = "_fs.dat"; break;
      case NOISE: AdjExt = "_fwh.dat"; break;
      case HEAT_LOAD: AdjExt = "_Q.dat"; break;
		}
		filename.append(AdjExt);
		restart_file.open(filename.data(), ios::in);
    
		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
			cout << "There is no adjoint restart file!! " << filename.data() << "."<< endl;
			cout << "Press any key to exit..." << endl;
			cin.get(); exit(1);
		}
    
		/*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
		long *Global2Local;
		Global2Local = new long[geometry->GetGlobal_nPointDomain()];
		/*--- First, set all indices to a negative value by default ---*/
		for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
			Global2Local[iPoint] = -1;
		}
		/*--- Now fill array with the transform values only for local points ---*/
		for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
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
        if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
        if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
				node[iPoint_Local] = new CAdjNSVariable(Solution, nDim, nVar, config);
			}
			iPoint_Global++;
		}
    
		/*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CAdjNSVariable(Solution, nDim, nVar, config);
		}
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
	}
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
}

CAdjTNE2NSSolver::~CAdjTNE2NSSolver(void) {
  
}


void CAdjTNE2NSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
  
  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the flow problem may use different methods). ---*/
  bool implicit       = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool upwind_2nd     = (config->GetKind_Upwind_AdjFlow() == ROE_2ND) ||
  (config->GetKind_Upwind_AdjFlow() == SW_2ND);
  bool center         = (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst     = (config->GetKind_Centered_AdjFlow() == JST);
  bool limiter        = (config->GetKind_SlopeLimit() != NONE);
  double adj_limit    = config->GetAdjointLimit();
  
  /*--- Compute nacelle inflow and exhaust properties ---*/
  GetNacelle_Properties(geometry, config, iMesh);
  
	/*--- Residual initialization ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Set the primitive variables incompressible and compressible
     adjoint variables ---*/
    //		node[iPoint]->SetPrimVar_Compressible(adj_limit);
    
		/*--- Initialize the convective residual vector ---*/
		LinSysRes.SetBlock_Zero(iPoint);
    
	}
  
  /*--- Compute gradients for upwind second-order reconstruction ---*/
  if ((upwind_2nd) && (iMesh == MESH_0)) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    
    /*--- Limiter computation ---*/
		if (limiter) SetSolution_Limiter(geometry, config);
	}
  
  /*--- Artificial dissipation for centered schemes ---*/
  if (center) {
    if ((center_jst) && (iMesh == MESH_0)) {
      SetDissipation_Switch(geometry, config);
      SetUndivided_Laplacian(geometry, config);
      if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
      if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    }
  }
  
	/*--- Compute gradients for solution reconstruction and viscous term
   (be careful, if an upwind strategy is used, then we compute the gradient twice) ---*/
	switch (config->GetKind_Gradient_Method()) {
    case GREEN_GAUSS :
      SetSolution_Gradient_GG(geometry, config);
      if ((config->GetKind_Solver() == ADJ_RANS) && (!config->GetFrozen_Visc()))
        solver_container[ADJTURB_SOL]->SetSolution_Gradient_GG(geometry, config);
      break;
    case WEIGHTED_LEAST_SQUARES :
      SetSolution_Gradient_LS(geometry, config);
      if ((config->GetKind_Solver() == ADJ_RANS) && (!config->GetFrozen_Visc()))
        solver_container[ADJTURB_SOL]->SetSolution_Gradient_LS(geometry, config);
      break;
	}
  
	/*--- Initialize the Jacobian for implicit integration ---*/
	if (implicit) Jacobian.SetValZero();
  
}

void CAdjTNE2NSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                        CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iPoint, jPoint, iEdge;
  
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge, coordinates and normal vector---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction and adjoint variables w/o reconstruction---*/
    numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                              solver_container[FLOW_SOL]->node[jPoint]->GetSolution());
    numerics->SetAdjointVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    
    /*--- Gradient of Adjoint Variables ---*/
    numerics->SetAdjointVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
    
    /*--- Viscosity and eddy viscosity---*/
    numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                  solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
    
    numerics->SetEddyViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
                               solver_container[FLOW_SOL]->node[jPoint]->GetEddyViscosity());
    
    /*--- Compute residual in a non-conservative way, and update ---*/
    numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
    
    /*--- Update adjoint viscous residual ---*/
    LinSysRes.SubtractBlock(iPoint, Residual_i);
    LinSysRes.AddBlock(jPoint, Residual_j);
    
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
      Jacobian.AddBlock(jPoint, iPoint, Jacobian_ji);
      Jacobian.AddBlock(jPoint, jPoint, Jacobian_jj);
    }
    
  }
  
}

void CAdjTNE2NSSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                       CConfig *config, unsigned short iMesh) {
	unsigned long iPoint;
  unsigned short iVar, iDim;
  
  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
  
	/*--- Loop over all the points, note that we are supposing that primitive and
	 adjoint gradients have been computed previously ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
		/*--- Conservative variables w/o reconstruction ---*/
		numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);
    
		/*--- Gradient of primitive and adjoint variables ---*/
		numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
		numerics->SetAdjointVarGradient(node[iPoint]->GetGradient(), NULL);
    
		/*--- Laminar viscosity, and eddy viscosity (adjoint with frozen viscosity) ---*/
		numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
		numerics->SetEddyViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(), 0.0);
    
		/*--- Set temperature of the fluid ---*/
		numerics->SetTemperature(solver_container[FLOW_SOL]->node[iPoint]->GetTemperature(), 0.0);
    
		/*--- Set volume ---*/
		numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    
		/*--- If turbulence computation we must add some coupling terms to the NS adjoint eq. ---*/
		if ((config->GetKind_Solver() == ADJ_RANS) && (!config->GetFrozen_Visc())) {
      
			/*--- Turbulent variables w/o reconstruction ---*/
			numerics->SetTurbVar(solver_container[TURB_SOL]->node[iPoint]->GetSolution(), NULL);
      
			/*--- Gradient of Turbulent Variables w/o reconstruction ---*/
			numerics->SetTurbVarGradient(solver_container[TURB_SOL]->node[iPoint]->GetGradient(), NULL);
      
			/*--- Turbulent adjoint variables w/o reconstruction ---*/
			numerics->SetTurbAdjointVar(solver_container[ADJTURB_SOL]->node[iPoint]->GetSolution(), NULL);
      
			/*--- Gradient of Adjoint turbulent variables w/o reconstruction ---*/
			numerics->SetTurbAdjointGradient(solver_container[ADJTURB_SOL]->node[iPoint]->GetGradient(), NULL);
      
			/*--- Set distance to the surface ---*/
      //			numerics->SetDistance(geometry->node[iPoint]->GetWallDistance(), 0.0);
      
		}
    
		/*--- Compute residual ---*/
		numerics->ComputeResidual(Residual, config);
    
    /*--- Add and substract to the residual ---*/
		LinSysRes.AddBlock(iPoint, Residual);
	}
  
  // WARNING: The rotating frame source term has been placed in the second
  // source term container since the section below is commented. This needs a
  // permanent fix asap!
  
  //  if ((config->GetKind_Solver() == ADJ_RANS) && (!config->GetFrozen_Visc()) && (config->GetKind_Adjoint() == CONTINUOUS)) {
  //    unsigned long jPoint, iEdge;
  //
  //    /*--- Gradient of primitive variables already computed in the previous step ---*/
  //    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
  //
  //      /*--- Points in edge, and normal vector ---*/
  //      iPoint = geometry->edge[iEdge]->GetNode(0);
  //      jPoint = geometry->edge[iEdge]->GetNode(1);
  //      second_numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
  //
  //      /*--- Conservative variables w/o reconstruction ---*/
  //      second_numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(),
  //                                     solver_container[FLOW_SOL]->node[jPoint]->GetSolution());
  //
  //      /*--- Gradient of primitive variables w/o reconstruction ---*/
  //      second_numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(),
  //                                        solver_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive());
  //
  //      /*--- Viscosity ---*/
  //      second_numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
  //                                         solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
  //
  //      /*--- Turbulent variables w/o reconstruction ---*/
  //      second_numerics->SetTurbVar(solver_container[TURB_SOL]->node[iPoint]->GetSolution(),
  //                                solver_container[TURB_SOL]->node[jPoint]->GetSolution());
  //
  //      /*--- Turbulent adjoint variables w/o reconstruction ---*/
  //      second_numerics->SetTurbAdjointVar(solver_container[ADJTURB_SOL]->node[iPoint]->GetSolution(),
  //                                       solver_container[ADJTURB_SOL]->node[jPoint]->GetSolution());
  //
  //      /*--- Set distance to the surface ---*/
  //      second_numerics->SetDistance(geometry->node[iPoint]->GetWallDistance(), geometry->node[jPoint]->GetWallDistance());
  //
  //      /*--- Add and Subtract Residual ---*/
  //      for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
  //      second_numerics->ComputeResidual(Residual, config);
  //      LinSysRes.AddBlock(iPoint, Residual);
  //      LinSysRes.SubtractBlock(jPoint, Residual);
  //    }
  //
  //  }
  
}

void CAdjTNE2NSSolver::Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
	unsigned long iVertex, iPoint;
	unsigned short iDim, jDim, iMarker;
	double **PsiVar_Grad, **PrimVar_Grad, div_phi, *Normal, Area,
	normal_grad_psi5, normal_grad_T, sigma_partial,
  cp, Laminar_Viscosity, heat_flux_factor, LevelSet, Target_LevelSet, temp_sens;
  
  double *Psi, *U, Enthalpy, **GridVel_Grad, gradPsi5_v, psi5_tau_partial, psi5_tau_grad_vel, source_v;
  double Density, Pressure, div_vel, val_turb_ke, vartheta, vartheta_partial, psi5_p_div_vel, *Omega, rho_v[3], CrossProduct[3];
  double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  
	double Gas_Constant = config->GetGas_ConstantND();
  
	cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  
	double *UnitaryNormal = new double[nDim];
	double *normal_grad_vel = new double[nDim];
	double *tang_deriv_psi5 = new double[nDim];
	double *tang_deriv_T = new double[nDim];
	double **Sigma = new double* [nDim];
  
	for (iDim = 0; iDim < nDim; iDim++)
		Sigma[iDim] = new double [nDim];
  
  double *normal_grad_gridvel = new double[nDim];
  double *normal_grad_v_ux =new double[nDim];
  double **Sigma_Psi5v = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Sigma_Psi5v[iDim] = new double [nDim];
  double **tau = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		tau[iDim] = new double [nDim];
  double *Velocity = new double[nDim];
  
  /*--- Compute gradient of adjoint variables on the surface ---*/
  SetSurface_Gradient(geometry, config);
  
  Total_Sens_Geo = 0.0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    Sens_Geo[iMarker] = 0.0;
    
    if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL)) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          
          PsiVar_Grad = node[iPoint]->GetGradient();
          PrimVar_Grad = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
          
          Laminar_Viscosity  = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
          
          heat_flux_factor = cp * Laminar_Viscosity / PRANDTL;
          
          /*--- Compute face area and the nondimensional normal to the surface ---*/
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) { Area += Normal[iDim]*Normal[iDim]; } Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++) { UnitaryNormal[iDim] = Normal[iDim] / Area; }
          
          /*--- Compute the sensitivity related to the temperature ---*/
          normal_grad_psi5 = 0.0; normal_grad_T = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            normal_grad_psi5 += PsiVar_Grad[nVar-1][iDim]*UnitaryNormal[iDim];
            normal_grad_T += PrimVar_Grad[0][iDim]*UnitaryNormal[iDim];
          }
          
          temp_sens = 0.0;
          if (config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) {
            
            /*--- Heat Flux Term: temp_sens = (\partial_tg \psi_5)\cdot (k \partial_tg T) ---*/
            for (iDim = 0; iDim < nDim; iDim++) {
              tang_deriv_psi5[iDim] = PsiVar_Grad[nVar-1][iDim] - normal_grad_psi5*UnitaryNormal[iDim];
              tang_deriv_T[iDim] = PrimVar_Grad[0][iDim] - normal_grad_T*UnitaryNormal[iDim];
            }
            for (iDim = 0; iDim < nDim; iDim++)
              temp_sens += heat_flux_factor * tang_deriv_psi5[iDim] * tang_deriv_T[iDim];
            
          } else if (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) {
            
            /*--- Isothermal Term: temp_sens = - k * \partial_n(\psi_5) * \partial_n(T) ---*/
            temp_sens = - heat_flux_factor * normal_grad_psi5 * normal_grad_T;
            
          }
          
          /*--- Term: sigma_partial = \Sigma_{ji} n_i \partial_n v_j ---*/
          div_phi = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            div_phi += PsiVar_Grad[iDim+1][iDim];
            for (jDim = 0; jDim < nDim; jDim++)
              Sigma[iDim][jDim] = Laminar_Viscosity * (PsiVar_Grad[iDim+1][jDim]+PsiVar_Grad[jDim+1][iDim]);
          }
          for (iDim = 0; iDim < nDim; iDim++)
            Sigma[iDim][iDim] -= TWO3*Laminar_Viscosity * div_phi;
          
          for (iDim = 0; iDim < nDim; iDim++) {
            normal_grad_vel[iDim] = 0.0;
            for (jDim = 0; jDim < nDim; jDim++)
              normal_grad_vel[iDim] += PrimVar_Grad[iDim+1][jDim]*UnitaryNormal[jDim];
          }
          
          sigma_partial = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            for (jDim = 0; jDim < nDim; jDim++)
              sigma_partial += UnitaryNormal[iDim]*Sigma[iDim][jDim]*normal_grad_vel[jDim];
          
          /*--- Compute sensitivity for each surface point ---*/
          CSensitivity[iMarker][iVertex] = (sigma_partial - temp_sens)*Area;
          Sens_Geo[iMarker] -= CSensitivity[iMarker][iVertex]*Area;
        }
      }
      Total_Sens_Geo += Sens_Geo[iMarker];
    }
  }
  
	delete [] UnitaryNormal;
	delete [] normal_grad_vel;
	delete [] tang_deriv_psi5;
	delete [] tang_deriv_T;
	for (iDim = 0; iDim < nDim; iDim++)
		delete Sigma[iDim];
	delete [] Sigma;
  
  delete [] normal_grad_gridvel;
  delete [] normal_grad_v_ux;
  for (iDim = 0; iDim < nDim; iDim++)
		delete Sigma_Psi5v[iDim];
  for (iDim = 0; iDim < nDim; iDim++)
		delete tau[iDim];
	delete [] tau;
  delete [] Velocity;
}

void CAdjTNE2NSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, total_index, Point_Normal;
	unsigned short iDim, iVar, jVar, jDim;
	double *d, *U, l1psi, mu_dyn, Temp, dVisc_T, rho, pressure, div_phi,
  force_stress, Sigma_5, **PsiVar_Grad, phi[3];
  double phis1, phis2, sq_vel, ProjVel, Enthalpy, *GridVel, phi_u, d_n;
  double Energy, ViscDens, XiDens, Density, SoundSpeed, Pressure, dPhiE_dn, Laminar_Viscosity, Eddy_Viscosity,
  Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz,
  Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5,
  Sigma_yz5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;
  
  
  double *Psi = new double[nVar];
	double **Tau = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Tau[iDim] = new double [nDim];
  double *Velocity = new double[nDim];
  double *Normal = new double[nDim];
  
  double **GradPhi = new double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    GradPhi[iDim] = new double [nDim];
  double *GradPsiE = new double [nDim];
  
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  
	double Gas_Constant = config->GetGas_ConstantND();
	double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Retrieve adjoint solution at the wall boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Psi[iVar] = node[iPoint]->GetSolution(iVar);
      
			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Get the force projection vector (based on the objective function) ---*/
			d = node[iPoint]->GetForceProj_Vector();
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv_i[iVar] = 0.0;
        Res_Visc_i[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar ++)
            Jacobian_ii[iVar][jVar] = 0.0;
        }
      }
      
      /*--- Adjustments to strong boundary condition for dynamic meshes ---*/
      if ( grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) {
          phi[iDim] = d[iDim] - Psi[nVar-1]*GridVel[iDim];
        }
      } else {
        for (iDim = 0; iDim < nDim; iDim++) {
          phi[iDim] = d[iDim];
        }
      }
      
			/*--- Strong BC imposition for the adjoint velocity equations ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      
      node[iPoint]->SetVel_ResTruncError_Zero();
			for (iDim = 0; iDim < nDim; iDim++)
				node[iPoint]->SetSolution_Old(iDim+1, phi[iDim]);
      
			/*--- Modify the velocity rows of the Jacobian ---*/
			if (implicit) {
				for (iVar = 1; iVar <= nDim; iVar++) {
					total_index = iPoint*nVar+iVar;
					Jacobian.DeleteValsRowi(total_index);
				}
			}
      
      /*--- Additional contributions to adjoint density and energy (weak imposition) ---*/
      
      
      /*--- Energy resiudal due to the convective term ---*/
      l1psi = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        l1psi += Normal[iDim]*d[iDim];
      Res_Conv_i[nVar-1] = l1psi*Gamma_Minus_One;
      
      /*--- Components of the effective and adjoint stress tensors ---*/
      PsiVar_Grad = node[iPoint]->GetGradient();
      div_phi = 0;
      for (iDim = 0; iDim < nDim; iDim++) {
        div_phi += PsiVar_Grad[iDim+1][iDim];
        for (jDim = 0; jDim < nDim; jDim++)
          Tau[iDim][jDim] = (PsiVar_Grad[iDim+1][jDim]+PsiVar_Grad[jDim+1][iDim]);
      }
      for (iDim = 0; iDim < nDim; iDim++)
        Tau[iDim][iDim] -= TWO3*div_phi;
      
      /*--- force_stress = n_i \Tau_{ij} d_j ---*/
      force_stress = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          force_stress += Normal[iDim]*Tau[iDim][jDim]*d[jDim];
      
      /*--- \partial \mu_dyn \partial T ---*/
      mu_dyn = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
      Temp = solver_container[FLOW_SOL]->node[iPoint]->GetTemperature();
      dVisc_T = 0.0;  // dVisc_T = mu_dyn*(Temp+3.0*mu2)/(2.0*Temp*(Temp+mu2));
      
      /*--- \Sigma_5 ---*/
      Sigma_5 = (Gamma/Cp)*dVisc_T*force_stress;
      
      /*--- Imposition of residuals ---*/
      rho = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
      pressure = solver_container[FLOW_SOL]->node[iPoint]->GetPressure(false);
      Res_Conv_i[0] = pressure*Sigma_5/(Gamma_Minus_One*rho*rho);
      Res_Conv_i[nVar-1] -= Sigma_5/rho;
      
      /*--- Flux contribution and Jacobian contributions for moving
       walls. Note that these are only for the adjoint density and
       adjoint energy equations (the adjoint vel. uses a strong BC). ---*/
      if (grid_movement) {
        
        /*--- Get the appropriate grid velocity at this node ---*/
        if (grid_movement)
          GridVel = geometry->node[iPoint]->GetGridVel();
        
        /*--- Get the enthalpy from the direct solution ---*/
        Enthalpy = solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy();
        
        /*--- Compute projections, velocity squared divided by two, and
         other inner products. Note that we are imposing v = u_wall from
         the direct problem and that phi = d - \psi_5 * v ---*/
        ProjVel = 0.0; sq_vel = 0.0; phi_u = 0.0; d_n = 0.0;
        phis1 = 0.0; phis2 = Psi[0] + Enthalpy * Psi[nVar-1];
        for (iDim = 0; iDim < nDim; iDim++) {
          ProjVel += GridVel[iDim]*Normal[iDim];
          sq_vel  += 0.5*GridVel[iDim]*GridVel[iDim];
          phis1   += Normal[iDim]*phi[iDim];
          phis2   += GridVel[iDim]*phi[iDim];
          phi_u   += GridVel[iDim]*phi[iDim];
          d_n     += d[iDim]*Normal[iDim];
        }
        phis1 += ProjVel * Psi[nVar-1];
        
        /*--- Convective flux at the wall node (adjoint density & energy only) ---*/
        
        /*--- Version 1 (full) ---*/
        //Res_Conv_i[0] = ProjVel * Psi[0] - phis2 * ProjVel + phis1 * Gamma_Minus_One * sq_vel - ProjVel*Psi[0];
        //Res_Conv_i[nVar-1] = ProjVel * Psi[nVar-1] + phis1 * Gamma_Minus_One - ProjVel*Psi[nVar-1];
        
        /*--- Simplified version ---*/
        Res_Conv_i[0] = -(Psi[0] + phi_u + Psi[nVar-1]*Enthalpy)*ProjVel + d_n*Gamma_Minus_One*sq_vel;
        Res_Conv_i[nVar-1] = d_n * Gamma_Minus_One;
        
        /*--- TO DO: Implicit contributions for convective part ---*/
        
        
        /*--- Viscous flux contributions at the wall node ---*/
        U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
        Laminar_Viscosity = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
        Eddy_Viscosity = solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(); // Should be zero at the wall
        Density = U[0];
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = GridVel[iDim];
        }
        Energy = U[nDim+1] / Density;
        SoundSpeed = sqrt(Gamma*Gamma_Minus_One*(Energy-sq_vel));
        Pressure = (SoundSpeed * SoundSpeed * Density) / Gamma;
        ViscDens = (Laminar_Viscosity + Eddy_Viscosity) / Density;
        XiDens = Gamma * (Laminar_Viscosity/PRANDTL + Eddy_Viscosity/PRANDTL_TURB) / Density;
        
        /*--- Average of the derivatives of the adjoint variables ---*/
        PsiVar_Grad = node[iPoint]->GetGradient();
        
        for (iDim = 0; iDim < nDim; iDim++) {
          GradPsiE[iDim] =  PsiVar_Grad[nVar-1][iDim];
          for (jDim = 0; jDim < nDim; jDim++)
            GradPhi[iDim][jDim] =  PsiVar_Grad[iDim+1][jDim];
        }
        
        /*--- Impose dPhiE_dn = 0 (adiabatic walls with frozen viscosity). Note
         that this is where a different adjoint boundary condition for temperature
         could be imposed. ---*/
        dPhiE_dn = 0.0;
        
        if (nDim ==2) {
          
          /*--- Compute the adjoint stress tensor ---*/
          Sigma_xx  = ViscDens * (FOUR3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1]);
          Sigma_yy  = ViscDens * (-TWO3 * GradPhi[0][0] + FOUR3 * GradPhi[1][1]);
          Sigma_xy  = ViscDens * (GradPhi[1][0] + GradPhi[0][1]);
          Sigma_xx5 = ViscDens * ( FOUR3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1]);
          Sigma_yy5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] + FOUR3 * Velocity[1] * GradPsiE[1]);
          Sigma_xy5 = ViscDens * (Velocity[0] * GradPsiE[1] + Velocity[1] * GradPsiE[0]);
          Sigma_5   = XiDens * dPhiE_dn;
          eta_xx    = Sigma_xx + Sigma_xx5;
          eta_yy    = Sigma_yy + Sigma_yy5;
          eta_xy    = Sigma_xy + Sigma_xy5;
          
          /*--- Viscous flux at the wall node (adjoint density & energy only) ---*/
          Res_Visc_i[0] = - (Velocity[0] * Normal[0] * eta_xx  + Velocity[1] * Normal[1] * eta_yy
                             + (Velocity[0] * Normal[1] + Velocity[1] * Normal[0]) * eta_xy
                             - (sq_vel - Pressure/(Density*Gamma_Minus_One)) * Sigma_5);
          Res_Visc_i[1] = 0.0;
          Res_Visc_i[2] = 0.0;
          Res_Visc_i[3] = Sigma_5;
          
        } else if (nDim == 3) {
          
          /*--- Compute the adjoint stress tensor ---*/
          Sigma_xx  = ViscDens * (FOUR3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1] - TWO3  * GradPhi[2][2]);
          Sigma_yy  = ViscDens * (-TWO3 * GradPhi[0][0] + FOUR3 * GradPhi[1][1] - TWO3  * GradPhi[2][2]);
          Sigma_zz  = ViscDens * (-TWO3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1] + FOUR3 * GradPhi[2][2]);
          Sigma_xy  = ViscDens * (GradPhi[1][0] + GradPhi[0][1]);
          Sigma_xz  = ViscDens * (GradPhi[2][0] + GradPhi[0][2]);
          Sigma_yz  = ViscDens * (GradPhi[2][1] + GradPhi[1][2]);
          Sigma_xx5 = ViscDens * ( FOUR3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1] -  TWO3 * Velocity[2] * GradPsiE[2]);
          Sigma_yy5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] + FOUR3 * Velocity[1] * GradPsiE[1] -  TWO3 * Velocity[2] * GradPsiE[2]);
          Sigma_zz5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1] + FOUR3 * Velocity[2] * GradPsiE[2]);
          Sigma_xy5 = ViscDens * (Velocity[0] * GradPsiE[1] + Velocity[1] * GradPsiE[0]);
          Sigma_xz5 = ViscDens * (Velocity[0] * GradPsiE[2] + Velocity[2] * GradPsiE[0]);
          Sigma_yz5 = ViscDens * (Velocity[1] * GradPsiE[2] + Velocity[2] * GradPsiE[1]);
          Sigma_5   = XiDens * dPhiE_dn;
          eta_xx    = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
          eta_xy    = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;
          
          /*--- Viscous flux at the wall node (adjoint density & energy only) ---*/
          Res_Visc_i[0] = - (Velocity[0] * Normal[0] * eta_xx  + Velocity[1] * Normal[1] * eta_yy + Velocity[2] * Normal[2] * eta_zz
                             + (Velocity[0] * Normal[1] + Velocity[1] * Normal[0]) * eta_xy
                             + (Velocity[0] * Normal[2] + Velocity[2] * Normal[0]) * eta_xz
                             + (Velocity[2] * Normal[1] + Velocity[1] * Normal[2]) * eta_yz
                             - (sq_vel - Pressure/(Density*Gamma_Minus_One)) * Sigma_5);
          Res_Visc_i[1] = 0.0;
          Res_Visc_i[2] = 0.0;
          Res_Visc_i[3] = 0.0;
          Res_Visc_i[4] = Sigma_5;
        }
      }
      
      /*--- Update convective and viscous residuals ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
      LinSysRes.SubtractBlock(iPoint, Res_Visc_i);
      if (implicit) {
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      }
      
		}
    
	}
  
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Tau[iDim];
	delete [] Tau;
  delete [] Psi;
  delete [] Velocity;
  delete [] Normal;
  delete [] GradPsiE;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] GradPhi[iDim];
  delete [] GradPhi;
  
}


void CAdjTNE2NSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, total_index, Point_Normal;
	unsigned short iDim, iVar, jVar, jDim;
	double *d, q, *U, l1psi, mu_dyn, Temp, dVisc_T, rho, pressure, div_phi,
  force_stress, Sigma_5, **PsiVar_Grad, phi[3];
  double phis1, phis2, sq_vel, ProjVel, Enthalpy, *GridVel, phi_u, d_n;
  double Energy, ViscDens, XiDens, Density, SoundSpeed, Pressure, dPhiE_dn, Laminar_Viscosity, Eddy_Viscosity,
  Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz,
  Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5,
  Sigma_yz5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;
  
  double *Psi = new double[nVar];
	double **Tau = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Tau[iDim] = new double [nDim];
  double *Velocity = new double[nDim];
  double *Normal = new double[nDim];
  
  double **GradPhi = new double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    GradPhi[iDim] = new double [nDim];
  double *GradPsiE = new double [nDim];
  
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  bool heat_flux_obj  = (config->GetKind_ObjFunc() == HEAT_FLUX);
  
	double Gas_Constant = config->GetGas_ConstantND();
	double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Retrieve adjoint solution at the wall boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Psi[iVar] = node[iPoint]->GetSolution(iVar);
      
			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Get the force projection vector (based on the objective function) ---*/
			d = node[iPoint]->GetForceProj_Vector();
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv_i[iVar] = 0.0;
        Res_Visc_i[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar ++)
            Jacobian_ii[iVar][jVar] = 0.0;
        }
      }
      
      /*--- Adjustments to strong boundary condition for dynamic meshes ---*/
      if ( grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) {
          phi[iDim] = d[iDim] - Psi[nVar-1]*GridVel[iDim];
        }
      } else {
        for (iDim = 0; iDim < nDim; iDim++) {
          phi[iDim] = d[iDim];
        }
      }
      
			/*--- Strong BC imposition for the adjoint velocity equations ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();
			for (iDim = 0; iDim < nDim; iDim++)
				node[iPoint]->SetSolution_Old(iDim+1, phi[iDim]);
      
			/*--- Modify the velocity rows of the Jacobian ---*/
			if (implicit) {
				for (iVar = 1; iVar <= nDim; iVar++) {
					total_index = iPoint*nVar+iVar;
					Jacobian.DeleteValsRowi(total_index);
				}
			}
      
      /*--- Strong enforcement of the adjoint energy equation for temperature-based objectives ---*/
      if (heat_flux_obj) {
        // This is the derivative of our objective function j = kdndT for heat flux with a negative sign from the formulation of the adjoint boundary conditions.
        q = -1.0;
      } else {
        q = 0.0;
      }
      
      /*--- Strong enforcement of the energy equations ---*/
      LinSysRes.SetBlock_Zero(iPoint, nVar-1);
      node[iPoint]->SetEnergy_ResTruncError_Zero();
      node[iPoint]->SetSolution_Old(nDim+1, q);
      
      /*--- Modify the energy rows of the Jacobian ---*/
      if (implicit) {
        iVar = nDim+1;
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
      
      /*--- Additional contributions to adjoint density and energy (weak imposition) ---*/
      /*--- Components of the effective and adjoint stress tensors ---*/
      PsiVar_Grad = node[iPoint]->GetGradient();
      div_phi = 0;
      for (iDim = 0; iDim < nDim; iDim++) {
        div_phi += PsiVar_Grad[iDim+1][iDim];
        for (jDim = 0; jDim < nDim; jDim++)
          Tau[iDim][jDim] = (PsiVar_Grad[iDim+1][jDim]+PsiVar_Grad[jDim+1][iDim]);
      }
      for (iDim = 0; iDim < nDim; iDim++)
        Tau[iDim][iDim] -= TWO3*div_phi;
      
      /*--- force_stress = n_i \Tau_{ij} d_j ---*/
      force_stress = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          force_stress += Normal[iDim]*Tau[iDim][jDim]*d[jDim];
      
      /*--- \partial \mu_dyn \partial T ---*/
      mu_dyn = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
      Temp = solver_container[FLOW_SOL]->node[iPoint]->GetTemperature();
      dVisc_T = 0.0;  // dVisc_T = mu_dyn*(Temp+3.0*mu2)/(2.0*Temp*(Temp+mu2));
      
      /*--- \Sigma_5 ---*/
      Sigma_5 = (Gamma/Cp)*dVisc_T*force_stress;
      
      /*--- Imposition of residuals ---*/
      rho = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
      pressure = solver_container[FLOW_SOL]->node[iPoint]->GetPressure(false);
      Res_Conv_i[0] = pressure*Sigma_5/(Gamma_Minus_One*rho*rho);
      
      /*--- Flux contribution and Jacobian contributions for moving
       walls. Note that these are only for the adjoint density and
       adjoint energy equations (the adjoint vel. uses a strong BC). ---*/
      if (grid_movement) {
        
        /*--- Get the appropriate grid velocity at this node ---*/
        if (grid_movement)
          GridVel = geometry->node[iPoint]->GetGridVel();
        
        /*--- Get the enthalpy from the direct solution ---*/
        Enthalpy = solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy();
        
        /*--- Compute projections, velocity squared divided by two, and
         other inner products. Note that we are imposing v = u_wall from
         the direct problem and that phi = d - \psi_5 * v ---*/
        ProjVel = 0.0; sq_vel = 0.0; phi_u = 0.0; d_n = 0.0;
        phis1 = 0.0; phis2 = Psi[0] + Enthalpy * Psi[nVar-1];
        for (iDim = 0; iDim < nDim; iDim++) {
          ProjVel += GridVel[iDim]*Normal[iDim];
          sq_vel  += 0.5*GridVel[iDim]*GridVel[iDim];
          phis1   += Normal[iDim]*phi[iDim];
          phis2   += GridVel[iDim]*phi[iDim];
          phi_u   += GridVel[iDim]*phi[iDim];
          d_n     += d[iDim]*Normal[iDim];
        }
        phis1 += ProjVel * Psi[nVar-1];
        
        /*--- Convective flux at the wall node (adjoint density & energy only) ---*/
        
        /*--- Version 1 (full) ---*/
        //Res_Conv_i[0] = ProjVel * Psi[0] - phis2 * ProjVel + phis1 * Gamma_Minus_One * sq_vel - ProjVel*Psi[0];
        //Res_Conv_i[nVar-1] = ProjVel * Psi[nVar-1] + phis1 * Gamma_Minus_One - ProjVel*Psi[nVar-1];
        
        /*--- Simplified version ---*/
        Res_Conv_i[0] = -(Psi[0] + phi_u + Psi[nVar-1]*Enthalpy)*ProjVel + d_n*Gamma_Minus_One*sq_vel;
        
        /*--- TO DO: Implicit contributions for convective part ---*/
        
        
        /*--- Viscous flux contributions at the wall node ---*/
        U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
        Laminar_Viscosity = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
        Eddy_Viscosity = solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(); // Should be zero at the wall
        Density = U[0];
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = GridVel[iDim];
        }
        Energy = U[nDim+1] / Density;
        SoundSpeed = sqrt(Gamma*Gamma_Minus_One*(Energy-sq_vel));
        Pressure = (SoundSpeed * SoundSpeed * Density) / Gamma;
        ViscDens = (Laminar_Viscosity + Eddy_Viscosity) / Density;
        XiDens = Gamma * (Laminar_Viscosity/PRANDTL + Eddy_Viscosity/PRANDTL_TURB) / Density;
        
        /*--- Average of the derivatives of the adjoint variables ---*/
        PsiVar_Grad = node[iPoint]->GetGradient();
        
        for (iDim = 0; iDim < nDim; iDim++) {
          GradPsiE[iDim] =  PsiVar_Grad[nVar-1][iDim];
          for (jDim = 0; jDim < nDim; jDim++)
            GradPhi[iDim][jDim] =  PsiVar_Grad[iDim+1][jDim];
        }
        
        /*--- Impose dPhiE_dn = 0 (adiabatic walls with frozen viscosity). Note
         that this is where a different adjoint boundary condition for temperature
         could be imposed. ---*/
        dPhiE_dn = 0.0;
        
        if (nDim ==2) {
          
          /*--- Compute the adjoint stress tensor ---*/
          Sigma_xx  = ViscDens * (FOUR3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1]);
          Sigma_yy  = ViscDens * (-TWO3 * GradPhi[0][0] + FOUR3 * GradPhi[1][1]);
          Sigma_xy  = ViscDens * (GradPhi[1][0] + GradPhi[0][1]);
          Sigma_xx5 = ViscDens * ( FOUR3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1]);
          Sigma_yy5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] + FOUR3 * Velocity[1] * GradPsiE[1]);
          Sigma_xy5 = ViscDens * (Velocity[0] * GradPsiE[1] + Velocity[1] * GradPsiE[0]);
          Sigma_5   = XiDens * dPhiE_dn;
          eta_xx    = Sigma_xx + Sigma_xx5;
          eta_yy    = Sigma_yy + Sigma_yy5;
          eta_xy    = Sigma_xy + Sigma_xy5;
          
          /*--- Viscous flux at the wall node (adjoint density & energy only) ---*/
          Res_Visc_i[0] = - (Velocity[0] * Normal[0] * eta_xx  + Velocity[1] * Normal[1] * eta_yy
                             + (Velocity[0] * Normal[1] + Velocity[1] * Normal[0]) * eta_xy
                             - (sq_vel - Pressure/(Density*Gamma_Minus_One)) * Sigma_5);
          Res_Visc_i[1] = 0.0;
          Res_Visc_i[2] = 0.0;
          
        } else if (nDim == 3) {
          
          /*--- Compute the adjoint stress tensor ---*/
          Sigma_xx  = ViscDens * (FOUR3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1] - TWO3  * GradPhi[2][2]);
          Sigma_yy  = ViscDens * (-TWO3 * GradPhi[0][0] + FOUR3 * GradPhi[1][1] - TWO3  * GradPhi[2][2]);
          Sigma_zz  = ViscDens * (-TWO3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1] + FOUR3 * GradPhi[2][2]);
          Sigma_xy  = ViscDens * (GradPhi[1][0] + GradPhi[0][1]);
          Sigma_xz  = ViscDens * (GradPhi[2][0] + GradPhi[0][2]);
          Sigma_yz  = ViscDens * (GradPhi[2][1] + GradPhi[1][2]);
          Sigma_xx5 = ViscDens * ( FOUR3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1] -  TWO3 * Velocity[2] * GradPsiE[2]);
          Sigma_yy5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] + FOUR3 * Velocity[1] * GradPsiE[1] -  TWO3 * Velocity[2] * GradPsiE[2]);
          Sigma_zz5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1] + FOUR3 * Velocity[2] * GradPsiE[2]);
          Sigma_xy5 = ViscDens * (Velocity[0] * GradPsiE[1] + Velocity[1] * GradPsiE[0]);
          Sigma_xz5 = ViscDens * (Velocity[0] * GradPsiE[2] + Velocity[2] * GradPsiE[0]);
          Sigma_yz5 = ViscDens * (Velocity[1] * GradPsiE[2] + Velocity[2] * GradPsiE[1]);
          Sigma_5   = XiDens * dPhiE_dn;
          eta_xx    = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
          eta_xy    = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;
          
          /*--- Viscous flux at the wall node (adjoint density & energy only) ---*/
          Res_Visc_i[0] = - (Velocity[0] * Normal[0] * eta_xx  + Velocity[1] * Normal[1] * eta_yy + Velocity[2] * Normal[2] * eta_zz
                             + (Velocity[0] * Normal[1] + Velocity[1] * Normal[0]) * eta_xy
                             + (Velocity[0] * Normal[2] + Velocity[2] * Normal[0]) * eta_xz
                             + (Velocity[2] * Normal[1] + Velocity[1] * Normal[2]) * eta_yz
                             - (sq_vel - Pressure/(Density*Gamma_Minus_One)) * Sigma_5);
          Res_Visc_i[1] = 0.0;
          Res_Visc_i[2] = 0.0;
          Res_Visc_i[3] = 0.0;
        }
      }
      
      /*--- Update convective and viscous residuals ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
      LinSysRes.SubtractBlock(iPoint, Res_Visc_i);
      if (implicit) {
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      }
      
		}
    
	}
  
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Tau[iDim];
	delete [] Tau;
  delete [] Psi;
  delete [] Velocity;
  delete [] Normal;
  delete [] GradPsiE;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] GradPhi[iDim];
  delete [] GradPhi;
  
}

