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

CAdjEulerSolver::CAdjEulerSolver(void) : CSolver() {

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

CAdjEulerSolver::CAdjEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
	unsigned long iPoint, index, iVertex;
	string text_line, mesh_filename;
	unsigned short iDim, iVar, iMarker;
	ifstream restart_file;
	string filename, AdjExt;

	bool restart = config->GetRestart();
	bool incompressible = config->GetIncompressible();
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

	if (incompressible) nVar = nDim + 1;
	else nVar = nDim + 2;
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

	/*--- Jacobians and vector structures for discrete computations ---*/
	if (config->GetKind_Adjoint() == DISCRETE) {

		/*--- Point to point Jacobians ---*/
		Jacobian_i = new double* [nVar];
		Jacobian_j = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new double [nVar];
			Jacobian_j[iVar] = new double [nVar];
		}

		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
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
			node[iPoint] = new CAdjEulerVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, nDim, nVar, config);
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
				if (incompressible) {
					if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2];
					if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
				}
				else {
					if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
					if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
				}
				node[iPoint_Local] = new CAdjEulerVariable(Solution, nDim, nVar, config);
			}
			iPoint_Global++;
		}

		/*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CAdjEulerVariable(Solution, nDim, nVar, config);
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

CAdjEulerSolver::~CAdjEulerSolver(void) {
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

void CAdjEulerSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
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

void CAdjEulerSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
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

void CAdjEulerSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
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

void CAdjEulerSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
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


void CAdjEulerSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) {
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

void CAdjEulerSolver::Set_MPI_Dissipation_Switch(CGeometry *geometry, CConfig *config) {
	unsigned short iMarker, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
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
			nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new double [nBufferR_Vector];
      Buffer_Send_Lambda = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetSensor();
      }
      
#ifndef NO_MPI
      
      //      /*--- Send/Receive using non-blocking communications ---*/
      //      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Lambda, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
      //      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Lambda, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
      //      send_request.Wait(status);
      //      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Lambda, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        Buffer_Receive_Lambda[iVertex] = Buffer_Send_Lambda[iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Lambda;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        node[iPoint]->SetSensor(Buffer_Receive_Lambda[iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Lambda;
      
    }
    
	}
}

void CAdjEulerSolver::SetForceProj_Vector(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
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
  
	bool rotating_frame = config->GetRotating_Frame();
  bool grid_movement = config->GetGrid_Movement();
  
	ForceProj_Vector = new double[nDim];
  
	/*--- If we have a rotating frame problem or an unsteady problem with
   mesh motion, use special reference values for the force coefficients.
   Otherwise, use the freestream values, which is the standard convention. ---*/
	if (rotating_frame || grid_movement) {
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

void CAdjEulerSolver::SetIntBoundary_Jump(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	unsigned short iMarker, iVar, jVar, kVar, jc, jrjc, jrjcm1, jrjcp1, jr, jm, jrm1, jrjr, jrp1, jmjm, iDim, jDim;
	unsigned long iVertex, iPoint, iPointNearField, nPointNearField = 0;
	double factor = 1.0, AngleDouble, data, aux, *IntBound_Vector, *coord, u, v, sq_vel, *FlowSolution, A[5][5], M[5][5], AM[5][5], b[5], WeightSB, sum, MinDist = 1E6,
			Dist, DerivativeOF = 0.0, *Normal;
	short iPhiAngle = 0, IndexNF_inv[180], iColumn;
	ifstream index_file;
	string text_line;
	vector<vector<double> > NearFieldWeight;
	vector<double> CoordNF;
	vector<short> IndexNF;

	bool incompressible = config->GetIncompressible();

	IntBound_Vector = new double [nVar];

	/*--- If equivalent area objective function, read the value of
	 the derivative from a file, this is a preprocess of the direct solution ---*/

	if (config->GetKind_ObjFunc() == EQUIVALENT_AREA) {

		/*--- Read derivative of the objective function at the NearField from file ---*/
		index_file.open("WeightNF.dat", ios::in);
		if (index_file.fail()) {
			cout << "There is no Weight Nearfield Pressure file (WeightNF.dat)." << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}

		nPointNearField = 0;

		while (index_file) {
			string line;
			getline(index_file, line);
			istringstream is(line);

			/*--- The first row provides the azimuthal angle ---*/
			if (nPointNearField == 0) {
				is >> data; // The first column is related with the coordinate
				while (is.good()) { is >> data; IndexNF.push_back(int(data)); }
			}
			else {
				is >> data; CoordNF.push_back(data); // The first column is the point coordinate
				vector<double> row;  
				while (is.good()) { is >> data; row.push_back(data); }
				NearFieldWeight.push_back(row);				
			}
			nPointNearField++;
		}

		/*--- Note tha the first row is the azimuthal angle ---*/
		nPointNearField = nPointNearField - 1;

		for (iPhiAngle = 0; iPhiAngle < 180; iPhiAngle++)
			IndexNF_inv[iPhiAngle] = -1;

		for (unsigned short iIndex = 0; iIndex < IndexNF.size(); iIndex++)
			IndexNF_inv[IndexNF[iIndex]] = iIndex;

	}

	/*--- Compute the jump on the adjoint variables for the upper and the lower side ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

				double Area = 0.0; double UnitaryNormal[3];
				for (iDim = 0; iDim < nDim; iDim++)
					Area += Normal[iDim]*Normal[iDim];
				Area = sqrt (Area);

				for (iDim = 0; iDim < nDim; iDim++)
					UnitaryNormal[iDim] = Normal[iDim]/Area;

				if (geometry->node[iPoint]->GetDomain()) {

					coord = geometry->node[iPoint]->GetCoord();
					DerivativeOF = 0.0;

					/*--- Just in case the functional depend also on the surface pressure ---*/
					WeightSB = 1.0-config->GetWeightCd(); 

					double AoA, XcoordRot = 0.0, YcoordRot = 0.0, ZcoordRot = 0.0;

					if (nDim == 2) XcoordRot = coord[0];			
					if (nDim == 3) {
						/*--- Rotate the nearfield cylinder  ---*/
						AoA = -(config->GetAoA()*PI_NUMBER/180.0);
						XcoordRot = coord[0]*cos(AoA) - coord[2]*sin(AoA);
						YcoordRot = coord[1];
						ZcoordRot = coord[0]*sin(AoA) + coord[2]*cos(AoA);
					}

					switch (config->GetKind_ObjFunc()) {	
					case EQUIVALENT_AREA :

						if (nDim == 2) iPhiAngle = 0;
						if (nDim == 3) {
							/*--- Compute the azimuthal angle of the iPoint ---*/
							AngleDouble = atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER;
							iPhiAngle = (short) floor(AngleDouble + 0.5);
							if (iPhiAngle < 0) iPhiAngle = 180 + iPhiAngle;
						}

						if (iPhiAngle <= 60) {
							iColumn = IndexNF_inv[iPhiAngle];
							/*--- An azimuthal angle is not defined... this happens with MG levels ---*/
							if (iColumn < 0.0) {
								if (IndexNF_inv[iPhiAngle+1] > 0) { iColumn = IndexNF_inv[iPhiAngle+1]; goto end; }
								if (IndexNF_inv[iPhiAngle-1] > 0) { iColumn = IndexNF_inv[iPhiAngle-1]; goto end; }
								if (IndexNF_inv[iPhiAngle+2] > 0) { iColumn = IndexNF_inv[iPhiAngle+2]; goto end; }
								if (IndexNF_inv[iPhiAngle-2] > 0) { iColumn = IndexNF_inv[iPhiAngle-2]; goto end; }
								if (IndexNF_inv[iPhiAngle+3] > 0) { iColumn = IndexNF_inv[iPhiAngle+3]; goto end; }
								if (IndexNF_inv[iPhiAngle-3] > 0) { iColumn = IndexNF_inv[iPhiAngle-3]; goto end; }
								if (IndexNF_inv[iPhiAngle+4] > 0) { iColumn = IndexNF_inv[iPhiAngle+4]; goto end; }
								if (IndexNF_inv[iPhiAngle-4] > 0) { iColumn = IndexNF_inv[iPhiAngle-4]; goto end; }
							}

							end:

							if (iColumn < 0.0) { cout <<" An azimuthal angle is not defined..." << endl; }

							/*--- Find the value of the weight in the table, using the azimuthal angle  ---*/
							MinDist = 1E6;
							for (iPointNearField = 0; iPointNearField < nPointNearField; iPointNearField++) {
								Dist = fabs(CoordNF[iPointNearField] - XcoordRot);
								if (Dist <= MinDist) {
									MinDist = Dist;
									DerivativeOF = factor*WeightSB*NearFieldWeight[iPointNearField][iColumn];
								}
							}
						}
						else DerivativeOF = 0.0;

						if ((MinDist > 1E-6) || (coord[nDim-1] > 0.0)) DerivativeOF = 0.0;

						break;

					case NEARFIELD_PRESSURE :
						DerivativeOF = factor*WeightSB*(solver_container[FLOW_SOL]->node[iPoint]->GetPressure(incompressible)
								- solver_container[FLOW_SOL]->GetPressure_Inf());
						break;
					}

					/*--- Compute the jump of the adjoint variables (2D, and 3D problems) --*/
					if (nDim == 2) {

						FlowSolution = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
						u = FlowSolution[1]/FlowSolution[0]; v = FlowSolution[2]/FlowSolution[0];
						double Energy = FlowSolution[3]/FlowSolution[0];
						double Rho = FlowSolution[0];

						sq_vel = u*u+v*v;
						A[0][0] = 0.0;																					A[0][1] = 0.0;									A[0][2] = 1.0;																				A[0][3] = 0.0;
						A[1][0] = -u*v;																					A[1][1] = v;										A[1][2] = u;																					A[1][3] = 0.0;
						A[2][0] = 0.5*(Gamma-3.0)*v*v+0.5*Gamma_Minus_One*u*u;	A[2][1] = -Gamma_Minus_One*u;		A[2][2] = (3.0-Gamma)*v;															A[2][3] = Gamma_Minus_One;
						A[3][0] = -Gamma*v*Energy+Gamma_Minus_One*v*sq_vel;					A[3][1] = -Gamma_Minus_One*u*v; A[3][2] = Gamma*Energy-0.5*Gamma_Minus_One*(u*u+3.0*v*v);	A[3][3] = Gamma*v;

						M[0][0] = 1.0;				M[0][1] = 0.0;		M[0][2] = 0.0;		M[0][3] = 0.0;
						M[1][0] = u;					M[1][1] = Rho;		M[1][2] = 0.0;		M[1][3] = 0.0;
						M[2][0] = v;					M[2][1] = 0.0;		M[2][2] = Rho;		M[2][3] = 0.0;
						M[3][0] = 0.5*sq_vel;	M[3][1] = Rho*u;	M[3][2] = Rho*v;	M[3][3] = 1.0/Gamma_Minus_One;

						for (iVar = 0; iVar < 4; iVar++)
							for (jVar = 0; jVar < 4; jVar++) {
								aux = 0.0;
								for (kVar = 0; kVar < 4; kVar++)
									aux += A[iVar][kVar]*M[kVar][jVar];
								AM[iVar][jVar] = aux;
							}

						for (iVar = 0; iVar < nVar; iVar++)
							for (jVar = 0; jVar < nVar; jVar++)
								A[iVar][jVar] = AM[jVar][iVar];

						b[0] = 0.0; b[1] = 0.0; b[2] = 0.0; b[3] = DerivativeOF; 
					}

					if (nDim == 3) {

						FlowSolution = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
						double velocity[3];
						double Energy = FlowSolution[4]/FlowSolution[0];
						double Rho = FlowSolution[0];

						double sqvel = 0.0; double proj_vel = 0.0;
						for (iDim = 0; iDim < nDim; iDim++) {
							velocity[iDim] = FlowSolution[iDim+1]/FlowSolution[0];
							sqvel    += velocity[iDim]*velocity[iDim];
							proj_vel += velocity[iDim]*UnitaryNormal[iDim];
						}

						double phi = 0.5*Gamma_Minus_One*sqvel;
						double a1 = Gamma*Energy-phi;
						double a2 = Gamma-1.0;

						/*--- Compute the projected Jacobian ---*/
						A[0][0] = 0.0;
						for (iDim = 0; iDim < nDim; iDim++) A[0][iDim+1] = UnitaryNormal[iDim];
						A[0][nDim+1] = 0.0;

						for (iDim = 0; iDim < nDim; iDim++) {
							A[iDim+1][0] = (UnitaryNormal[iDim]*phi - velocity[iDim]*proj_vel);
							for (jDim = 0; jDim < nDim; jDim++)
								A[iDim+1][jDim+1] = (UnitaryNormal[jDim]*velocity[iDim]-a2*UnitaryNormal[iDim]*velocity[jDim]);
							A[iDim+1][iDim+1] += proj_vel;
							A[iDim+1][nDim+1] = a2*UnitaryNormal[iDim];
						}

						A[nDim+1][0] = proj_vel*(phi-a1);
						for (iDim = 0; iDim < nDim; iDim++)
							A[nDim+1][iDim+1] = (UnitaryNormal[iDim]*a1-a2*velocity[iDim]*proj_vel);
						A[nDim+1][nDim+1] = Gamma*proj_vel;

						/*--- Compute the transformation matrix ---*/
						M[0][0] = 1.0;					M[0][1] = 0.0;							M[0][2] = 0.0;							M[0][3] = 0.0;							M[0][4] = 0.0;
						M[1][0] = velocity[0];	M[1][1] = Rho;							M[1][2] = 0.0;							M[1][3] = 0.0;							M[1][4] = 0.0;
						M[2][0] = velocity[1];	M[2][1] = 0.0;							M[2][2] = Rho;							M[2][3] = 0.0;							M[2][4] = 0.0;
						M[3][0] = velocity[2];	M[3][1] = 0.0;							M[3][2] = 0.0;							M[3][3] = Rho;							M[3][4] = 0.0;
						M[4][0] = 0.5*sqvel;		M[4][1] = Rho*velocity[0];	M[4][2] = Rho*velocity[1];	M[4][3] = Rho*velocity[2];	M[4][4] = 1.0/Gamma_Minus_One;

						/*--- Compute A times M ---*/
						for (iVar = 0; iVar < 5; iVar++)
							for (jVar = 0; jVar < 5; jVar++) {
								aux = 0.0;
								for (kVar = 0; kVar < 5; kVar++)
									aux += A[iVar][kVar]*M[kVar][jVar];
								AM[iVar][jVar] = aux;
							}

						/*--- Compute the transpose matrix ---*/
						for (iVar = 0; iVar < nVar; iVar++)
							for (jVar = 0; jVar < nVar; jVar++)
								A[iVar][jVar] = AM[jVar][iVar];

						/*--- Create the soruce term (AM)^T X = b ---*/
						b[0] = 0.0; b[1] = 0.0; b[2] = 0.0; b[3] = 0.0; b[4] = DerivativeOF;
					}

					/*--- Solve the linear system using a LU descomposition --*/
					for (jc = 1; jc < nVar; jc++)
						A[0][jc] /= A[0][0];

					jrjc = 0;						
					for (;;) {
						jrjc++; jrjcm1 = jrjc-1; jrjcp1 = jrjc+1;
						for (jr = jrjc; jr < nVar; jr++) {
							sum = A[jr][jrjc];
							for (jm = 0; jm <= jrjcm1; jm++)
								sum -= A[jr][jm]*A[jm][jrjc];
							A[jr][jrjc] = sum;
						}
						if ( jrjc == (nVar-1) ) goto stop;
						for (jc = jrjcp1; jc<nVar; jc++) {
							sum = A[jrjc][jc];
							for (jm = 0; jm <= jrjcm1; jm++)
								sum -= A[jrjc][jm]*A[jm][jc];
							A[jrjc][jc] = sum/A[jrjc][jrjc];
						}
					}

					stop:

					b[0] = b[0]/A[0][0];
					for (jr = 1; jr<nVar; jr++) {
						jrm1 = jr-1;
						sum = b[jr];
						for (jm = 0; jm<=jrm1; jm++)
							sum -= A[jr][jm]*b[jm];
						b[jr] = sum/A[jr][jr];
					}

					for (jrjr = 1; jrjr<nVar; jrjr++) {
						jr = (nVar-1)-jrjr;
						jrp1 = jr+1;
						sum = b[jr];
						for (jmjm = jrp1; jmjm<nVar; jmjm++) {
							jm = (nVar-1)-jmjm+jrp1;
							sum -= A[jr][jm]*b[jm];
						}
						b[jr] = sum;
					}

					/*--- Update the internal boundary jump --*/
					for (iVar = 0; iVar < nVar; iVar++)
						IntBound_Vector[iVar] = b[iVar];

					node[iPoint]->SetIntBoundary_Jump(IntBound_Vector);
				}
			}

	delete [] IntBound_Vector;

}

void CAdjEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
	unsigned long iPoint, Point_Fine;
	unsigned short iMesh, iChildren, iVar;
	double LevelSet, Area_Children, Area_Parent, LevelSet_Fine, *Solution, *Solution_Fine;
  
	bool restart = config->GetRestart();
	bool freesurface = config->GetFreeSurface();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
	if (freesurface) {
    
    for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
      
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        
        /*--- Set initial boundary condition at iter 0 ---*/
        if ((ExtIter == 0) && (!restart)) {
          
          /*--- Compute the adjoint level set value in all the MG levels ---*/
          if (iMesh == MESH_0) {
            solver_container[iMesh][ADJLEVELSET_SOL]->node[iPoint]->SetSolution(0, 0.0);
          }
          else {
            Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
            LevelSet = 0.0;
            for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
              Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
              Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
              LevelSet_Fine = solver_container[iMesh-1][ADJLEVELSET_SOL]->node[Point_Fine]->GetSolution(0);
              LevelSet += LevelSet_Fine*Area_Children/Area_Parent;
            }
            solver_container[iMesh][ADJLEVELSET_SOL]->node[iPoint]->SetSolution(0, LevelSet);
          }
          
          /*--- Compute the flow solution using the level set value. ---*/
					for (iVar = 0; iVar < nVar; iVar++)
            solver_container[iMesh][ADJFLOW_SOL]->node[iPoint]->SetSolution(iVar, 0.0);
          
        }
      }
      
      /*--- Set the MPI communication ---*/
      solver_container[iMesh][ADJFLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
      solver_container[iMesh][ADJLEVELSET_SOL]->Set_MPI_Solution(geometry[iMesh], config);
      
      /*--- The value of the solution for the first iteration of the dual time ---*/
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        if ((ExtIter == 0) && (dual_time)) {
          solver_container[iMesh][ADJLEVELSET_SOL]->node[iPoint]->Set_Solution_time_n();
          solver_container[iMesh][ADJLEVELSET_SOL]->node[iPoint]->Set_Solution_time_n1();
        }
      }
      
    }
    
	}
  
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

void CAdjEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
  
  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the flow problem may use different methods). ---*/
  bool implicit       = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool upwind_2nd     = (config->GetKind_Upwind_AdjFlow() == ROE_2ND);
  bool center         = (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst     = (config->GetKind_Centered_AdjFlow() == JST);
  bool limiter        = (config->GetKind_SlopeLimit() != NONE);
  bool incompressible = config->GetIncompressible();
  double adj_limit    = config->GetAdjointLimit();

  /*--- Compute nacelle inflow and exhaust properties ---*/
  GetNacelle_Properties(geometry, config, iMesh);
  
	/*--- Residual initialization ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Set the primitive variables incompressible and compressible 
     adjoint variables ---*/
		if (incompressible) node[iPoint]->SetPrimVar_Incompressible(adj_limit);
		else node[iPoint]->SetPrimVar_Compressible(adj_limit);
    
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
	if ((implicit) || (config->GetKind_Adjoint() == DISCRETE) ) Jacobian.SetValZero();
  
}

void CAdjEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
		CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

	unsigned long iEdge, iPoint, jPoint;

	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool high_order_diss = ((config->GetKind_Centered_AdjFlow() == JST) && (iMesh == MESH_0));
	
	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();
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

		if (incompressible) {
			numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solver_container[FLOW_SOL]->node[jPoint]->GetDensityInc());
			numerics->SetBetaInc2(solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2(), solver_container[FLOW_SOL]->node[jPoint]->GetBetaInc2());
			numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
		}
		else {
			numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
					solver_container[FLOW_SOL]->node[jPoint]->GetSoundSpeed());
			numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
					solver_container[FLOW_SOL]->node[jPoint]->GetEnthalpy());
		}

		numerics->SetLambda(solver_container[FLOW_SOL]->node[iPoint]->GetLambda(), 
				solver_container[FLOW_SOL]->node[jPoint]->GetLambda());

		if (high_order_diss) {
			numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(), node[jPoint]->GetUndivided_Laplacian());
			numerics->SetSensor(solver_container[FLOW_SOL]->node[iPoint]->GetSensor(),
					solver_container[FLOW_SOL]->node[jPoint]->GetSensor());
		}

		/*--- Rotating frame - use lower-order dissipation sensor ---*/
		if (rotating_frame) {
			numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());
			numerics->SetRotFlux(geometry->edge[iEdge]->GetRotFlux());
			numerics->SetSensor(node[iPoint]->GetSensor(),node[jPoint]->GetSensor());
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


void CAdjEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {  
	double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, *Limiter_i = NULL,
			*Limiter_j = NULL, *Psi_i = NULL, *Psi_j = NULL, *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;

	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool high_order_diss = (((config->GetKind_Upwind_AdjFlow() == ROE_2ND) ||
                           (config->GetKind_Upwind_AdjFlow() == SW_2ND)) && (iMesh == MESH_0));
	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();
	bool limiter = (config->GetKind_SlopeLimit() != NONE);

	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

		if(config->GetKind_Adjoint() != DISCRETE) {
			/*--- Adjoint variables w/o reconstruction ---*/
			Psi_i = node[iPoint]->GetSolution(); Psi_j = node[jPoint]->GetSolution();
			numerics->SetAdjointVar(Psi_i, Psi_j);
		}

		/*--- Conservative variables w/o reconstruction ---*/
		U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
		U_j = solver_container[FLOW_SOL]->node[jPoint]->GetSolution();
		numerics->SetConservative(U_i, U_j);

		if (incompressible) {
			numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
					solver_container[FLOW_SOL]->node[jPoint]->GetDensityInc());
			numerics->SetBetaInc2(solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2(),
					solver_container[FLOW_SOL]->node[jPoint]->GetBetaInc2());
			numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
		}
		else {		
			numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
					solver_container[FLOW_SOL]->node[jPoint]->GetSoundSpeed());
			numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
					solver_container[FLOW_SOL]->node[jPoint]->GetEnthalpy());
		}

		/*--- Rotating frame ---*/
		if (rotating_frame) {
			numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());
			numerics->SetRotFlux(geometry->edge[iEdge]->GetRotFlux());
		}

		/*--- Mesh motion ---*/
		if (grid_movement) {
			numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
		}

		/*--- High order reconstruction using MUSCL strategy ---*/
		if ((high_order_diss) && (config->GetKind_Adjoint() != DISCRETE)) {
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
		}

		/*--- Compute the residual---*/
		if (config->GetKind_Adjoint() == DISCRETE)
			numerics->ComputeResidual(Jacobian_i, Jacobian_j, config);
		else
			numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

		/*--- Add and Subtract Residual ---*/
		if (config->GetKind_Adjoint() == DISCRETE) {
			if (!high_order_diss) {
				// Transpose of block positions
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
				Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_i);
				Jacobian.AddBlock(jPoint, iPoint, Jacobian_j);
				Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
			} else { // include effect of reconstruction

				// NOT MPI READY (nodes in separate domain may not get picked up)

				// get list of Normals and solution values
				//				double **Normals, **U_js;
				//
				//				nNeigh = node[iPoint]->GetnPoint();
				//
				//				Normals = new double*[nNeigh];
				//				for (iNeigh = 0; iNeigh < nNeigh; iNeigh++)
				//					Normals[iNeigh] = new double*[nDim];
				//
				//				U_js = new double*[nNeigh];
				//				for (iNeigh = 0; iNeigh < nNeigh; iNeigh++)
				//					U_js[iNeigh] = new double*[nVar];
				//
				//				for (iNeigh = 0; iNeigh < nNeigh; iNeigh++) {
				//
				//					kPoint = node[iPoint]->GetPoint(iNeigh);
				//
				//					kEdge = geometry->FindEdge(iPoint, kPoint);
				//
				//					kNormal = geometry->edge[kEdge]->GetNormal();
				//
				//					for (iDim = 0; iDim < nDim; iDim++)
				//						Normals[iNeigh][iDim] = kNormal[iDim];
				//
				//					U_k = solver_container[FLOW_SOL]->node[kPoint]->GetSolution();
				//
				//					for (iVar = 0; iVar < nVar; iVar++)
				//						U_js[iNeigh][iVar] = U_k[iVar];
				//
				//				}
				//
				//				for (iNeigh = 0; iNeigh < nNeigh; iNeigh++)
				//					delete [] Normals[iNeigh];
				//
				//				delete [] Normals;
				//
				//				for (iNeigh = 0; iNeigh < nNeigh; iNeigh++)
				//					delete [] U_js[iNeigh];
				//
				//				delete [] U_js;
				//
				//				nNeigh = node[jPoint]->GetnPoint();
				//
				//				for (iNeigh = 0; iNeigh < nNeigh; iNeigh++) {
				//
				//					kPoint = node[jPoint]->GetPoint(iNeigh);
				//
				//					kEdge = geometry->FindEdge(jPoint, kPoint);
				//
				//					kNormal = geometry->edge[kEdge]->GetNormal();
				//
				//					for (iDim = 0; iDim < nDim; iDim++)
				//						Normals[iNeigh][iDim] = kNormal[iDim];
				//
				//					U_k = solver_container[FLOW_SOL]->node[kPoint]->GetSolution();
				//
				//					for (iVar = 0; iVar < nVar; iVar++)
				//						U_js[iNeigh][iVar] = U_k[iVar];
				//
				//				}

			}
		}
		else {

			LinSysRes.SubtractBlock(iPoint, Residual_i);
			LinSysRes.SubtractBlock(jPoint, Residual_j);

			/*--- Implicit contribution to the residual ---*/
			if ((implicit) && (config->GetKind_Adjoint() != DISCRETE)) {
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
				Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
				Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ji);
				Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_jj);
			}
		}
	}
}

void CAdjEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
		CConfig *config, unsigned short iMesh) {

	unsigned short iVar, iDim;
	unsigned long iPoint;

	bool rotating_frame = config->GetRotating_Frame();
	bool axisymmetric = config->GetAxisymmetric();
	bool gravity = (config->GetGravityForce() == YES);
	bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
	bool freesurface = config->GetFreeSurface();

	for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;

	if (rotating_frame) {

		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) { 

			/*--- Set solution  ---*/
			numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

			/*--- Set control volume ---*/
			numerics->SetVolume(geometry->node[iPoint]->GetVolume());

			/*--- Set rotational velocity ---*/
			numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());

			/*--- Compute Residual ---*/
			numerics->ComputeResidual(Residual, Jacobian_i, config);

			/*--- Add Residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);
		}
	}

	if (time_spectral) {

		double Volume, Source;

		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

			/*--- Get control volume ---*/
			Volume = geometry->node[iPoint]->GetVolume();

			/*--- Get stored time spectral source term ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Source = node[iPoint]->GetTimeSpectral_Source(iVar);
				Residual[iVar] = Source*Volume;
			}

			/*--- Add Residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

		}
	}

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

	if (gravity) {

	}

	if (freesurface) {
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      double Volume = geometry->node[iPoint]->GetVolume();
      double **Gradient = solver_container[ADJLEVELSET_SOL]->node[iPoint]->GetGradient();
      double coeff = solver_container[LEVELSET_SOL]->node[iPoint]->GetSolution(0) / solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
      
      Residual[0] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Residual[iDim+1] = coeff*Gradient[0][iDim]*Volume;
      }
      
      LinSysRes.AddBlock(iPoint, Residual);
      
		}		
	}

}

void CAdjEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
		CConfig *config, unsigned short iMesh) {
}

void CAdjEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
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

void CAdjEulerSolver::SetDissipation_Switch(CGeometry *geometry, CConfig *config) {

	double dx = 0.1;
	double LimK = 0.03;
	double eps2 =  pow((LimK*dx),3);

	unsigned long iPoint, jPoint;
	unsigned short iNeigh, nNeigh, iDim;
	double **Gradient_i, *Coord_i, *Coord_j, diff_coord, dist_ij, r_u, r_u_ij, 
	du_max, du_min, u_ij, *Solution_i, *Solution_j, dp, dm;


	for (iPoint = 0; iPoint < nPoint; iPoint++) 

		if (geometry->node[iPoint]->GetDomain()) {

			Solution_i = node[iPoint]->GetSolution();
			Gradient_i = node[iPoint]->GetGradient();
			Coord_i = geometry->node[iPoint]->GetCoord();
			nNeigh = geometry->node[iPoint]->GetnPoint();

			/*--- Find max and min value of the variable in the control volume around the mesh point ---*/
			du_max = 1.0E-8; du_min = -1.0E-8;
			for (iNeigh = 0; iNeigh < nNeigh; iNeigh++) {
				jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
				Solution_j = node[jPoint]->GetSolution();
				du_max = max(du_max, Solution_j[0] - Solution_i[0]);
				du_min = min(du_min, Solution_j[0] - Solution_i[0]);
			}

			r_u = 1.0;
			for (iNeigh = 0; iNeigh < nNeigh; iNeigh++) {

				/*--- Unconstrained reconstructed solution ---*/
				jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
				Solution_j = node[jPoint]->GetSolution();
				Coord_j = geometry->node[jPoint]->GetCoord();
				u_ij = Solution_i[0]; dist_ij = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					diff_coord = Coord_j[iDim]-Coord_i[iDim];
					u_ij += 0.5*diff_coord*Gradient_i[0][iDim];
				}

				/*--- Venkatakrishnan limiter ---*/
				if ((u_ij - Solution_i[0]) >= 0.0) dp = du_max;
				else	dp = du_min;
				dm = u_ij - Solution_i[0];
				r_u_ij = (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2);

				/*--- Take the smallest value of the limiter ---*/
				r_u = min(r_u, r_u_ij);

			}
			node[iPoint]->SetSensor(1.0-r_u);
		}
  
  /*--- MPI parallelization ---*/
  Set_MPI_Dissipation_Switch(geometry, config);
  
}

void CAdjEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, 
		CConfig *config, unsigned short iRKStep) {
	double *Residual, *Res_TruncError, Vol, Delta, Res;
	unsigned short iVar;
	unsigned long iPoint;
  
	double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);

	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Delta = solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() / Vol;
    
		Res_TruncError = node[iPoint]->GetResTruncError();
		Residual = LinSysRes.GetBlock(iPoint);
    
		for (iVar = 0; iVar < nVar; iVar++) {
      Res = Residual[iVar] + Res_TruncError[iVar];
			node[iPoint]->AddSolution(iVar, -Res*Delta*RK_AlphaCoeff);
			AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex());
		}
    
	}

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CAdjEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
	unsigned short iVar;
	unsigned long iPoint;

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
      AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex());
		}
    
	}

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CAdjEulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
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

void CAdjEulerSolver::Solve_LinearSystem(CGeometry *geometry, CSolver **solver_container, CConfig *config){
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

void CAdjEulerSolver::Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
	unsigned long iVertex, iPoint, Neigh;
	unsigned short iPos, jPos;
	unsigned short iDim, iMarker, iNeigh;
	double *d = NULL, *Normal = NULL, *Psi = NULL, *U = NULL, Enthalpy, conspsi, Mach_Inf,
  Area, **PrimVar_Grad = NULL, **ConsVar_Grad = NULL, *ConsPsi_Grad = NULL,
  ConsPsi, d_press, grad_v, Beta2, v_gradconspsi, UnitaryNormal[3], *RotVel = NULL, *GridVel = NULL,
  LevelSet, Target_LevelSet;
	//double RefDensity, *RefVelocity = NULL, RefPressure;

	double r, ru, rv, rw, rE, p, T; // used in farfield sens
	double dp_dr, dp_dru, dp_drv, dp_drw, dp_drE; // used in farfield sens
	double dH_dr, dH_dru, dH_drv, dH_drw, dH_drE, H; // used in farfield sens
	//	double alpha, beta;
	double *USens, *U_infty;

	double Gas_Constant = config->GetGas_ConstantND();

	double **D, *Dd;
	D = new double*[nDim];
	for (iPos=0; iPos<nDim; iPos++)
		D[iPos] = new double[nDim];

	Dd = new double[nDim];

	USens = new double[nVar];
	U_infty = new double[nVar];

	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();
	bool grid_movement  = config->GetGrid_Movement();
    bool freesurface = config->GetFreeSurface();
  
	/*--- Initialize sensitivities to zero ---*/
	Total_Sens_Geo = 0.0; Total_Sens_Mach = 0.0; Total_Sens_AoA = 0.0;
	Total_Sens_Press = 0.0; Total_Sens_Temp = 0.0;
	//	Total_Sens_Far = 0.0;

	/*--- Compute surface sensitivity ---*/
	if (config->GetKind_Adjoint() != DISCRETE) {
    
		/*--- Loop over boundary markers to select those for Euler walls ---*/
		for (iMarker = 0; iMarker < nMarker; iMarker++)
			if (config->GetMarker_All_Boundary(iMarker) == EULER_WALL)

				/*--- Loop over points on the surface to store the auxiliary variable ---*/
				for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					if (geometry->node[iPoint]->GetDomain()) {
						Psi = node[iPoint]->GetSolution();
						U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
						if (incompressible) {
							Beta2 = solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2();
							conspsi = Beta2*Psi[0];
						} else {
							Enthalpy = solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy();
							conspsi = U[0]*Psi[0] + U[0]*Enthalpy*Psi[nDim+1];
						}
						for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];

						node[iPoint]->SetAuxVar(conspsi);

						/*--- Also load the auxiliary variable for first neighbors ---*/
						for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
							Neigh = geometry->node[iPoint]->GetPoint(iNeigh);
							Psi = node[Neigh]->GetSolution();
							U = solver_container[FLOW_SOL]->node[Neigh]->GetSolution();
							if (incompressible) {
								Beta2 = solver_container[FLOW_SOL]->node[Neigh]->GetBetaInc2();
								conspsi = Beta2*Psi[0];
							} else {
								Enthalpy = solver_container[FLOW_SOL]->node[Neigh]->GetEnthalpy();
								conspsi = U[0]*Psi[0] + U[0]*Enthalpy*Psi[nDim+1];
							}
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

						/*--- Adjustment for a rotating frame ---*/
						if (rotating_frame) RotVel = geometry->node[iPoint]->GetRotVel();

						/*--- Adjustment for grid movement - double check this ---*/
						if (grid_movement) GridVel = geometry->node[iPoint]->GetGridVel();

						d_press = 0.0; grad_v = 0.0; v_gradconspsi = 0.0;
						for (iDim = 0; iDim < nDim; iDim++) {
              
              /*-- Retrieve the value of the pressure gradient ---*/
              if (incompressible) d_press += d[iDim]*ConsVar_Grad[0][iDim];
							else d_press += d[iDim]*PrimVar_Grad[nDim+1][iDim];

              /*-- Retrieve the value of the velocity gradient ---*/
							grad_v += PrimVar_Grad[iDim+1][iDim]*ConsPsi;
              
              /*-- Retrieve the value of the theta gradient ---*/
							v_gradconspsi += solver_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim, incompressible) * ConsPsi_Grad[iDim];
							if (rotating_frame) v_gradconspsi -= RotVel[iDim] * ConsPsi_Grad[iDim];
							if (grid_movement) v_gradconspsi -= GridVel[iDim] * ConsPsi_Grad[iDim];
						}

            /*--- Compute additional term in the surface sensitivity for
             free surface problem. ---*/
            if (freesurface) {
              LevelSet = solver_container[LEVELSET_SOL]->node[iPoint]->GetSolution(0);
              Target_LevelSet = geometry->node[iPoint]->GetCoord(nDim-1);
              d_press += 0.5*(Target_LevelSet - LevelSet)*(Target_LevelSet - LevelSet);
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

	/*--- Farfield Sensitivity, only for compressible flows ---*/
  if (!incompressible) {
    
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (config->GetMarker_All_Boundary(iMarker) == FAR_FIELD) {
        Sens_Mach[iMarker] = 0.0;
        Sens_AoA[iMarker] = 0.0;
        Sens_Press[iMarker] = 0.0;
        Sens_Temp[iMarker] = 0.0;
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (geometry->node[iPoint]->GetDomain()) {
            Psi = node[iPoint]->GetSolution();
            U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            
            Mach_Inf   = config->GetMach_FreeStreamND();
            if (grid_movement || rotating_frame)
              Mach_Inf = config->GetMach_Motion();
            
            r = U[0]; ru = U[1]; rv = U[2];
            if (nDim == 2) { rw = 0.0; rE = U[3]; }
            else { rw = U[3]; rE = U[4]; }
            
            p = Gamma_Minus_One*(rE-(ru*ru + rv*rv + rw*rw)/(2*r));
            
            Area = 0.0; for (iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            Area = sqrt(Area);
            
            for (iDim = 0; iDim < nDim; iDim++)
              UnitaryNormal[iDim] = -Normal[iDim]/Area;
            
            if (config->GetKind_Adjoint() == CONTINUOUS) {
              
              H = (rE + p)/r;
              
              dp_dr = Gamma_Minus_One*(ru*ru + rv*rv + rw*rw)/(2*r*r);
              dp_dru = -Gamma_Minus_One*ru/r;
              dp_drv = -Gamma_Minus_One*rv/r;
              if (nDim == 2) {
                dp_drw = 0.0;
                dp_drE = Gamma_Minus_One;
              } else {
                dp_drw = -Gamma_Minus_One*rw/r;
                dp_drE = Gamma_Minus_One;
              }
              
              
              dH_dr = (-H + dp_dr)/r;
              dH_dru = dp_dru/r;
              dH_drv = dp_drv/r;
              if (nDim == 2) {
                dH_drw = 0.0;
                dH_drE = (1 + dp_drE)/r;
              } else {
                dH_drw = dp_drw/r;
                dH_drE = (1 + dp_drE)/r;
              }
              
              if (nDim == 2) {
                Jacobian_j[0][0] = 0.0;
                Jacobian_j[1][0] = Area*UnitaryNormal[0];
                Jacobian_j[2][0] = Area*UnitaryNormal[1];
                Jacobian_j[3][0] = 0.0;
                
                Jacobian_j[0][1] = (-(ru*ru)/(r*r) + dp_dr)*Area*UnitaryNormal[0] +
                (-(ru*rv)/(r*r))*Area*UnitaryNormal[1];
                Jacobian_j[1][1] = (2*ru/r + dp_dru)*Area*UnitaryNormal[0] +
                (rv/r)*Area*UnitaryNormal[1];
                Jacobian_j[2][1] = (dp_drv)*Area*UnitaryNormal[0] +
                (ru/r)*Area*UnitaryNormal[1];
                Jacobian_j[3][1] = (dp_drE)*Area*UnitaryNormal[0];
                
                Jacobian_j[0][2] = (-(ru*rv)/(r*r))*Area*UnitaryNormal[0] +
                (-(rv*rv)/(r*r) + dp_dr)*Area*UnitaryNormal[1];
                Jacobian_j[1][2] = (rv/r)*Area*UnitaryNormal[0] +
                (dp_dru)*Area*UnitaryNormal[1];
                Jacobian_j[2][2] = (ru/r)*Area*UnitaryNormal[0] +
                (2*rv/r + dp_drv)*Area*UnitaryNormal[1];
                Jacobian_j[3][2] = (dp_drE)*Area*UnitaryNormal[1];
                
                Jacobian_j[0][3] = (ru*dH_dr)*Area*UnitaryNormal[0] +
                (rv*dH_dr)*Area*UnitaryNormal[1];
                Jacobian_j[1][3] = (H + ru*dH_dru)*Area*UnitaryNormal[0] +
                (rv*dH_dru)*Area*UnitaryNormal[1];
                Jacobian_j[2][3] = (ru*dH_drv)*Area*UnitaryNormal[0] +
                (H + rv*dH_drv)*Area*UnitaryNormal[1];
                Jacobian_j[3][3] = (ru*dH_drE)*Area*UnitaryNormal[0] +
                (rv*dH_drE)*Area*UnitaryNormal[1];
              } else {
                Jacobian_j[0][0] = 0.0;
                Jacobian_j[1][0] = Area*UnitaryNormal[0];
                Jacobian_j[2][0] = Area*UnitaryNormal[1];
                Jacobian_j[3][0] = Area*UnitaryNormal[2];
                Jacobian_j[4][0] = 0.0;
                
                Jacobian_j[0][1] = (-(ru*ru)/(r*r) + dp_dr)*Area*UnitaryNormal[0] +
                (-(ru*rv)/(r*r))*Area*UnitaryNormal[1] +
                (-(ru*rw)/(r*r))*Area*UnitaryNormal[2];
                Jacobian_j[1][1] = (2*ru/r + dp_dru)*Area*UnitaryNormal[0] +
                (rv/r)*Area*UnitaryNormal[1] +
                (rw/r)*Area*UnitaryNormal[2];
                Jacobian_j[2][1] = (dp_drv)*Area*UnitaryNormal[0] +
                (ru/r)*Area*UnitaryNormal[1];
                Jacobian_j[3][1] = (dp_drw)*Area*UnitaryNormal[0] +
                (ru/r)*Area*UnitaryNormal[2];
                Jacobian_j[4][1] = (dp_drE)*Area*UnitaryNormal[0];
                
                Jacobian_j[0][2] = (-(ru*rv)/(r*r))*Area*UnitaryNormal[0] +
                (-(rv*rv)/(r*r) + dp_dr)*Area*UnitaryNormal[1] +
                (-(rv*rw)/(r*r))*Area*UnitaryNormal[2];
                Jacobian_j[1][2] = (rv/r)*Area*UnitaryNormal[0] +
                (dp_dru)*Area*UnitaryNormal[1];
                Jacobian_j[2][2] = (ru/r)*Area*UnitaryNormal[0] +
                (2*rv/r + dp_drv)*Area*UnitaryNormal[1] +
                (rw/r)*Area*UnitaryNormal[2];
                Jacobian_j[3][2] = (dp_drw)*Area*UnitaryNormal[1] +
                (rv/r)*Area*UnitaryNormal[2];
                Jacobian_j[4][2] = (dp_drE)*Area*UnitaryNormal[1];
                
                Jacobian_j[0][3] = (-(ru*rw)/(r*r))*Area*UnitaryNormal[0] +
                (-(rv*rw)/(r*r))*Area*UnitaryNormal[1] +
                (-(rw*rw)/(r*r) + dp_dr)*Area*UnitaryNormal[2];
                Jacobian_j[1][3] = (rw/r)*Area*UnitaryNormal[0] +
                (dp_dru)*Area*UnitaryNormal[2];
                Jacobian_j[2][3] = (rw/r)*Area*UnitaryNormal[1] +
                (dp_drv)*Area*UnitaryNormal[2];
                Jacobian_j[3][3] = (ru/r)*Area*UnitaryNormal[0] +
                (rv/r)*Area*UnitaryNormal[1] +
                (2*rw/r + dp_drw)*Area*UnitaryNormal[2];
                Jacobian_j[4][3] = (dp_drE)*Area*UnitaryNormal[2];
                
                Jacobian_j[0][4] = (ru*dH_dr)*Area*UnitaryNormal[0] +
                (rv*dH_dr)*Area*UnitaryNormal[1] +
                (rw*dH_dr)*Area*UnitaryNormal[2];
                Jacobian_j[1][4] = (H + ru*dH_dru)*Area*UnitaryNormal[0] +
                (rv*dH_dru)*Area*UnitaryNormal[1] +
                (rw*dH_dru)*Area*UnitaryNormal[2];
                Jacobian_j[2][4] = (ru*dH_drv)*Area*UnitaryNormal[0] +
                (H + rv*dH_drv)*Area*UnitaryNormal[1] +
                (rw*dH_drv)*Area*UnitaryNormal[2];
                Jacobian_j[3][4] = (ru*dH_drw)*Area*UnitaryNormal[0] +
                (rv*dH_drw)*Area*UnitaryNormal[1] +
                (H + rw*dH_drw)*Area*UnitaryNormal[2];
                Jacobian_j[4][4] = (ru*dH_drE)*Area*UnitaryNormal[0] +
                (rv*dH_drE)*Area*UnitaryNormal[1] +
                (rw*dH_drE)*Area*UnitaryNormal[2];
              }
              
            }
            else if (config->GetKind_Adjoint() == DISCRETE) {
              
							/*--- Flow Solution at infinity ---*/
							U_infty[0] = solver_container[FLOW_SOL]->GetDensity_Inf();
							U_infty[1] = solver_container[FLOW_SOL]->GetDensity_Velocity_Inf(0);
							U_infty[2] = solver_container[FLOW_SOL]->GetDensity_Velocity_Inf(1);
							U_infty[3] = solver_container[FLOW_SOL]->GetDensity_Energy_Inf();
							if (nDim == 3) {
								U_infty[3] = solver_container[FLOW_SOL]->GetDensity_Velocity_Inf(2);
								U_infty[4] = solver_container[FLOW_SOL]->GetDensity_Energy_Inf();
							}
              numerics->SetConservative(U, U_infty);
              for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
              numerics->SetNormal(Normal);
              
              if (incompressible) {
                numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                                      solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
                numerics->SetBetaInc2(solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2(),
                                    solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2());
                numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
              }
              else {
                numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(),
                                      solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed());
                numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(),
                                    solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy());
              }
              
              /*--- Compute the upwind flux ---*/
              numerics->ComputeResidual(Jacobian_i, Jacobian_j, config);
              
            }
            
            // Mach
            USens[0] = 0.0; USens[1] = ru/Mach_Inf; USens[2] = rv/Mach_Inf;
            if (nDim == 2) { USens[3] = Gamma*Mach_Inf*p; }
            else { USens[3] = rw/Mach_Inf; USens[4] = Gamma*Mach_Inf*p; }
            
            for (iPos = 0; iPos < nVar; iPos++)
              for (jPos = 0; jPos < nVar; jPos++) {
                Sens_Mach[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
              }
            
            // Alpha
            USens[0] = 0.0;
            if (nDim == 2) { USens[1] = -rv; USens[2] = ru; USens[3] = 0.0; }
            else { USens[1] = -rw; USens[2] = 0.0; USens[3] = ru; USens[4] = 0.0; }
            
            for (iPos = 0; iPos < nVar; iPos++)
              for (jPos = 0; jPos < nVar; jPos++) {
                Sens_AoA[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
              }
            
            // Pressure
            USens[0] = r/p;
            USens[1] = ru/p;
            USens[2] = rv/p;
            if (nDim == 2) {
              USens[3] = rE/p;
            } else {
              USens[3] = rw/p;
              USens[4] = rE/p;
            }
            
            for (iPos = 0; iPos < nVar; iPos++)
              for (jPos = 0; jPos < nVar; jPos++) {
                Sens_Press[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
              }
            
            // Temperature
            
            T = p/(r*Gas_Constant);
            USens[0] = -r/T;
            USens[1] = 0.5*ru/T;
            USens[2] = 0.5*rv/T;
            if (nDim == 2) {
              USens[3] = (ru*ru + rv*rv + rw*rw)/(r*T);
            } else {
              USens[3] = 0.5*rw/T;
              USens[4] = (ru*ru + rv*rv + rw*rw)/(r*T);
            }
            
            for (iPos = 0; iPos < nVar; iPos++)
              for (jPos = 0; jPos < nVar; jPos++) {
                Sens_Temp[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
              }
            
          }
        }
        Total_Sens_Mach -= Sens_Mach[iMarker];
        Total_Sens_AoA -= Sens_AoA[iMarker];
        Total_Sens_Press -= Sens_Press[iMarker];
        Total_Sens_Temp -= Sens_Temp[iMarker];
      }
    }
    
    // Explicit contribution from farfield quantity (Cl or Cd)
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (config->GetMarker_All_Boundary(iMarker) == EULER_WALL) {
        
        //Sens_Far = 0.0;
        Sens_Mach[iMarker] = 0.0;
        Sens_AoA[iMarker] = 0.0;
        Sens_Press[iMarker] = 0.0;
        Sens_Temp[iMarker] = 0.0;
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (geometry->node[iPoint]->GetDomain()) {
            U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            p = solver_container[FLOW_SOL]->node[iPoint]->GetPressure(incompressible);
            
            Mach_Inf   = config->GetMach_FreeStreamND();
            
            d = node[iPoint]->GetForceProj_Vector();
            
            Area = 0.0; for (iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            Area = sqrt(Area);
            
            for (iDim = 0; iDim < nDim; iDim++)
              UnitaryNormal[iDim] = -Normal[iDim]/Area;
            
            
            // Mach
            for (iPos=0; iPos<nDim; iPos++)
              Dd[iPos] = -(2/Mach_Inf)*d[iPos];
            
            for (iPos=0; iPos<nDim; iPos++)
              Sens_Mach[iMarker] += p*Dd[iPos]*Area*UnitaryNormal[iPos];
            
            // Alpha
            if (nDim == 2) {
              D[0][0] = 0.0;
              D[0][1] = -1.0;
              
              D[1][0] = 1.0;
              D[1][1] = 0.0;
            } else {
              D[0][0] = 0.0;
              D[0][1] = 0.0;
              D[0][2] = -1.0;
              
              D[1][0] = 0.0;
              D[1][1] = 0.0;
              D[1][2] = 0.0;
              
              D[2][0] = 1.0;
              D[2][1] = 0.0;
              D[2][2] = 0.0;
            }
            
            for (iPos=0; iPos<nDim; iPos++)
              Dd[iPos] = 0.0;
            for (iPos=0; iPos<nDim; iPos++)
              for (jPos=0; jPos<nDim; jPos++)
                Dd[iPos] += D[iPos][jPos]*d[jPos];
            
            for (iPos=0; iPos<nDim; iPos++)
              Sens_AoA[iMarker] += p*Dd[iPos]*Area*UnitaryNormal[iPos];
            
            // Pressure
            for (iPos=0; iPos<nDim; iPos++)
              Dd[iPos] = -(1/p)*d[iPos];
            
            for (iPos=0; iPos<nDim; iPos++)
              Sens_Press[iMarker] += p*Dd[iPos]*Area*UnitaryNormal[iPos];
            
            // Temperature
            for (iPos=0; iPos<nDim; iPos++)
              Dd[iPos] = 0.0;
            
            for (iPos=0; iPos<nDim; iPos++)
              Sens_Temp[iMarker] += p*Dd[iPos]*Area*UnitaryNormal[iPos];
            
          }
        }
        
        
        Total_Sens_Mach += Sens_Mach[iMarker];
        Total_Sens_AoA += Sens_AoA[iMarker];
        Total_Sens_Press += Sens_Press[iMarker];
        Total_Sens_Temp += Sens_Temp[iMarker];
      }
    }
  }

	for (iPos=0; iPos<nDim; iPos++)
		delete [] D[iPos];
	delete [] D;
	delete [] Dd;

	delete [] USens;
	delete [] U_infty;

}

void CAdjEulerSolver::Smooth_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
	unsigned short iMarker;
	unsigned long iVertex, jVertex, nVertex, iPoint;
	double **A, *b, Sens, *ArchLength, *Coord_begin, *Coord_end, dist;

	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == EULER_WALL) {
			nVertex = geometry->nVertex[iMarker];

			/*--- Allocate the linear system ---*/
			A = new double* [nVertex]; 
			b = new double [nVertex]; 
			ArchLength = new double [nVertex];
			for (iVertex = 0; iVertex < nVertex; iVertex++) {
				A[iVertex] = new double [nVertex];
			}

			/*--- Initialization ---*/
			for (iVertex = 0; iVertex < nVertex; iVertex++) {
				b[iVertex] = 0.0; ArchLength[iVertex] = 0.0;
				for (jVertex = 0; jVertex < nVertex; jVertex++)
					A[iVertex][jVertex] = 0.0;
			}

			/*--- Set the arch length ---*/
			ArchLength[0] = 0.0;
			for (iVertex = 1; iVertex < nVertex; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex-1]->GetNode();
				Coord_begin = geometry->node[iPoint]->GetCoord();
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Coord_end = geometry->node[iPoint]->GetCoord();
				dist = sqrt (pow( Coord_end[0]-Coord_begin[0], 2.0) + pow( Coord_end[1]-Coord_begin[1], 2.0));
				ArchLength[iVertex] = ArchLength[iVertex-1] + dist;
			}

			/*--- Remove the trailing edge effect ---*/
			double MinPosSens = 0.0; double MinNegSens = 0.0;
			for (iVertex = 0; iVertex < nVertex; iVertex++) {
				Sens = CSensitivity[iMarker][iVertex];
				if (ArchLength[iVertex] > ArchLength[nVertex-1]*0.01) { MinNegSens = Sens; break; }
			}

			for (iVertex = 0; iVertex < nVertex; iVertex++) {
				Sens = CSensitivity[iMarker][iVertex];
				if (ArchLength[iVertex] > ArchLength[nVertex-1]*0.99) { MinPosSens = Sens; break; }
			}

			for (iVertex = 0; iVertex < nVertex; iVertex++) {
				if (ArchLength[iVertex] < ArchLength[nVertex-1]*0.01)
					CSensitivity[iMarker][iVertex] = MinNegSens;
				if (ArchLength[iVertex] > ArchLength[nVertex-1]*0.99)
					CSensitivity[iMarker][iVertex] = MinPosSens;
			}

			/*--- Set the right hand side of the system ---*/
			for (iVertex = 0; iVertex < nVertex; iVertex++) {
				b[iVertex] = CSensitivity[iMarker][iVertex];
			}

			/*--- Set the mass matrix ---*/
			double Coeff = 0.0, BackDiff = 0.0, ForwDiff = 0.0, CentDiff = 0.0;
			double epsilon = 5E-5;
			for (iVertex = 0; iVertex < nVertex; iVertex++) {

				if ((iVertex != nVertex-1) && (iVertex != 0)) {
					BackDiff = (ArchLength[iVertex]-ArchLength[iVertex-1]);
					ForwDiff = (ArchLength[iVertex+1]-ArchLength[iVertex]);
					CentDiff = (ArchLength[iVertex+1]-ArchLength[iVertex-1]);
				}
				if (iVertex == nVertex-1) {
					BackDiff = (ArchLength[nVertex-1]-ArchLength[nVertex-2]);
					ForwDiff = (ArchLength[0]-ArchLength[nVertex-1]);
					CentDiff = (ArchLength[0]-ArchLength[nVertex-2]);					
				}
				if (iVertex == 0) {
					BackDiff = (ArchLength[0]-ArchLength[nVertex-1]);
					ForwDiff = (ArchLength[1]-ArchLength[0]);
					CentDiff = (ArchLength[1]-ArchLength[nVertex-1]);					
				}

				Coeff = epsilon*2.0/(BackDiff*ForwDiff*CentDiff);		

				A[iVertex][iVertex] = Coeff*CentDiff;

				if (iVertex != 0) A[iVertex][iVertex-1] = -Coeff*ForwDiff;
				else A[iVertex][nVertex-1] = -Coeff*ForwDiff;

				if (iVertex != nVertex-1) A[iVertex][iVertex+1] = -Coeff*BackDiff;
				else A[iVertex][0] = -Coeff*BackDiff;

			}

			/*--- Add the gradient value in the main diagonal ---*/
			for (iVertex = 0; iVertex < nVertex; iVertex++)
				A[iVertex][iVertex] += 1.0;			

			/*--- Dirichlet boundary condition ---*/
			unsigned long iVertex = int(nVertex/2);
			A[iVertex][iVertex] = 1.0; 
			A[iVertex][iVertex+1] = 0.0; 
			A[iVertex][iVertex-1] = 0.0; 

			Gauss_Elimination(A, b, nVertex);

			/*--- Set the new value of the sensitiviy ---*/
			for (iVertex = 0; iVertex < nVertex; iVertex++)
				CSensitivity[iMarker][iVertex] = b[iVertex];

			/*--- Deallocate the linear system ---*/
			for (iVertex = 0; iVertex < nVertex; iVertex++)
				delete [] A[iVertex];
			delete [] A;	
			delete [] b;
			delete [] ArchLength;

		}
	}


}

void CAdjEulerSolver::GetNacelle_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh) {
	unsigned short iDim, iMarker, iVar;
	unsigned long iVertex, iPoint;
	double Area, Flow_Dir[3], alpha;
  
  unsigned short nMarker_NacelleInflow = config->GetnMarker_NacelleInflow();
  unsigned short nMarker_NacelleExhaust = config->GetnMarker_NacelleExhaust();
  
  if ((nMarker_NacelleInflow != 0) || (nMarker_NacelleExhaust != 0)) {
    
    /*--- Check the flow orientation in the nacelle inflow ---*/
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_Boundary(iMarker) == NACELLE_EXHAUST) {
        
        /*--- Loop over all the vertices on this boundary marker ---*/
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          /*--- Normal vector for this vertex (negate for outward convention) ---*/
          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
          
          for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
          
          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Vector[iDim]*Vector[iDim];
          Area = sqrt (Area);
          
          /*--- Compute unitary vector ---*/
          for (iDim = 0; iDim < nDim; iDim++)
            Vector[iDim] /= Area;
          
          /*--- The flow direction is defined by the local velocity on the surface ---*/
          for (iDim = 0; iDim < nDim; iDim++)
            Flow_Dir[iDim] = node[iPoint]->GetSolution(iDim+1) / node[iPoint]->GetSolution(0);
          
          /*--- Dot product of normal and flow direction. ---*/
          alpha = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            alpha += Vector[iDim]*Flow_Dir[iDim];
          
          /*--- Flow in the wrong direction. ---*/
          if (alpha < 0.0) {
            
            /*--- Copy the old solution ---*/
            for (iVar = 0; iVar < nVar; iVar++)
              node[iPoint]->SetSolution(iVar, node[iPoint]->GetSolution_Old(iVar));
            
          }
          
        }
      }
    }
    
  }
  
}

void CAdjEulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	double *d = NULL, *Normal, *U, *Psi_Aux, ProjVel = 0.0, bcn, vn = 0.0, Area, *UnitaryNormal, *Coord;
	double *Velocity, *Psi, *ObjFuncSource, Enthalpy = 0.0, sq_vel, phin, phis1, phis2, DensityInc = 0.0, BetaInc2 = 0.0;
	unsigned short iDim, iVar, jDim;
	double *dPressure;
  
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();
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
			if(config->GetKind_Adjoint() != DISCRETE) {
				Psi_Aux = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];
			}
      
			/*--- Flow solution ---*/
			U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
      
      /*--- Read the value of the objective function ---*/
      d = node[iPoint]->GetForceProj_Vector();
      
      /*--- Normal vector computation ---*/
			Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
			Area = sqrt(Area);
			for (iDim = 0; iDim < nDim; iDim++) UnitaryNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Incompressible solver ---*/
			if (incompressible) {
        
        DensityInc = solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
        BetaInc2 = solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2();
        
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = U[iDim+1] / solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
        
        /*--- Compute projections ---*/
        bcn = 0.0; phin = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          bcn += d[iDim]*UnitaryNormal[iDim];
          phin += Psi[iDim+1]*UnitaryNormal[iDim];
        }
        
        /*--- Introduce the boundary condition ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Psi[iDim+1] -= ( phin - bcn ) * UnitaryNormal[iDim];
        
        /*--- Inner products after introducing BC (Psi has changed) ---*/
        phis1 = 0.0; phis2 = Psi[0] * (BetaInc2 / DensityInc);
        for (iDim = 0; iDim < nDim; iDim++) {
          phis1 -= Normal[iDim]*Psi[iDim+1];
          phis2 += Velocity[iDim]*Psi[iDim+1];
        }
        
        /*--- Flux of the Euler wall ---*/
        Residual[0] = phis1;
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[iDim+1] = - phis2 * Normal[iDim];
        
        /*--- Update residual ---*/
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        if (implicit) {
          
          /*--- Adjoint density ---*/
          Jacobian_ii[0][0] = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_ii[0][iDim+1] = - Normal[iDim];
          
          /*--- Adjoint velocities ---*/
          for (iDim = 0; iDim < nDim; iDim++) {
            Jacobian_ii[iDim+1][0] = -Normal[iDim] * (BetaInc2 / DensityInc) ;
            for (jDim = 0; jDim < nDim; jDim++)
              Jacobian_ii[iDim+1][jDim+1] = - Normal[iDim] * Velocity[jDim];
          }
          
          /*--- Update Jacobian ---*/
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
        }
        
			}
      
      /*--- Compressible solver ---*/
      else {
        
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
        
        /*--- Extra boundary term for a rotating frame ---*/
        if (rotating_frame) {
          double ProjRotVel = 0.0;
          double *RotVel = geometry->node[iPoint]->GetRotVel();
          for (iDim = 0; iDim < nDim; iDim++) {
            ProjRotVel += RotVel[iDim]*UnitaryNormal[iDim];
          }
          ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux()/Area;
          phin -= Psi[nVar-1]*ProjRotVel;
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
        
        /*--- Flux adjustment for a rotating Frame ---*/
        if (rotating_frame) {
          double ProjRotVel = 0.0;
          double *RotVel = geometry->node[iPoint]->GetRotVel();
          for (iDim = 0; iDim < nDim; iDim++)
            ProjRotVel -= RotVel[iDim]*Normal[iDim];
          ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux();
          Residual[0] -= ProjRotVel*Psi[0];
          for (iDim = 0; iDim < nDim; iDim++)
            Residual[iDim+1] -= ProjRotVel*Psi[iDim+1];
          Residual[nVar-1] -= ProjRotVel*Psi[nVar-1];
        }
        
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
          
          /*--- Jacobian contribution due to a rotating frame ---*/
          if (rotating_frame) {
            double ProjRotVel = 0.0;
            double *RotVel = geometry->node[iPoint]->GetRotVel();
            for (iDim = 0; iDim < nDim; iDim++)
              ProjRotVel -= RotVel[iDim]*Normal[iDim];
            ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux();
            Jacobian_ii[0][0] -= ProjRotVel;
            for (iDim = 0; iDim < nDim; iDim++)
              Jacobian_ii[iDim+1][iDim+1] -= ProjRotVel;
            Jacobian_ii[nVar-1][nVar-1] -= ProjRotVel;
          }
          
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
	}
  
	delete [] Velocity;
	delete [] UnitaryNormal;
	delete [] Psi;
  
}

void CAdjEulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
		CConfig *config, unsigned short val_marker) {

	unsigned long iVertex, iPoint;
	double *Normal, *U, *Psi_Aux, ProjVel = 0.0, vn = 0.0, Area, *UnitaryNormal, *Coord;
	double *Velocity, *Psi, Enthalpy = 0.0, sq_vel, phin, phis1, phis2, DensityInc = 0.0, BetaInc2 = 0.0;
	unsigned short iDim, iVar, jDim;

	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();
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

			if (incompressible) {

				DensityInc = solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
				BetaInc2 = solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2();

				for (iDim = 0; iDim < nDim; iDim++)
					Velocity[iDim] = U[iDim+1] / solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc();

				/*--- Compute projections ---*/
				phin = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					phin += Psi[iDim+1]*UnitaryNormal[iDim];

				/*--- Introduce the boundary condition ---*/
				for (iDim = 0; iDim < nDim; iDim++) 
					Psi[iDim+1] -= phin * UnitaryNormal[iDim];

				/*--- Inner products after introducing BC (Psi has changed) ---*/
				phis1 = 0.0; phis2 = Psi[0] * (BetaInc2 / DensityInc);
				for (iDim = 0; iDim < nDim; iDim++) {
					phis1 -= Normal[iDim]*Psi[iDim+1];
					phis2 += Velocity[iDim]*Psi[iDim+1];
				}

				/*--- Flux of the Euler wall ---*/
				Residual[0] = phis1;
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[iDim+1] = - phis2 * Normal[iDim];

			}

			else {

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

				/*--- Rotating Frame ---*/
				if (rotating_frame) {
					double ProjRotVel = 0.0;
					double *RotVel = geometry->node[iPoint]->GetRotVel();
					for (iDim = 0; iDim < nDim; iDim++) {
						ProjRotVel += RotVel[iDim]*UnitaryNormal[iDim];
					}
					ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux()/Area;
					phin -= Psi[nVar-1]*ProjRotVel;				
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

				/*--- Rotating Frame ---*/
				if (rotating_frame) {
					double ProjRotVel = 0.0;
					double *RotVel = geometry->node[iPoint]->GetRotVel();
					for (iDim = 0; iDim < nDim; iDim++)
						ProjRotVel -= RotVel[iDim]*Normal[iDim];
					ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux();
					Residual[0] -= ProjRotVel*Psi[0];
					for (iDim = 0; iDim < nDim; iDim++)
						Residual[iDim+1] -= ProjRotVel*Psi[iDim+1];
					Residual[nVar-1] -= ProjRotVel*Psi[nVar-1];
				}

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
			}

			/*--- Update residual ---*/
			LinSysRes.SubtractBlock(iPoint, Residual);


			/*--- Implicit stuff ---*/
			if (implicit) {

				if (incompressible) {

					/*--- Adjoint density ---*/
					Jacobian_ii[0][0] = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_ii[0][iDim+1] = - Normal[iDim];

					/*--- Adjoint velocities ---*/
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_ii[iDim+1][0] = -Normal[iDim] * (BetaInc2 / DensityInc) ;
						for (jDim = 0; jDim < nDim; jDim++)
							Jacobian_ii[iDim+1][jDim+1] = - Normal[iDim] * Velocity[jDim];
					}

				}

				else {

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

					/*--- Contribution from a rotating frame ---*/
					if (rotating_frame) {
						double ProjRotVel = 0.0;
						double *RotVel = geometry->node[iPoint]->GetRotVel();
						for (iDim = 0; iDim < nDim; iDim++)
							ProjRotVel -= RotVel[iDim]*Normal[iDim];
						ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux();
						Jacobian_ii[0][0] -= ProjRotVel;
						for (iDim = 0; iDim < nDim; iDim++)
							Jacobian_ii[iDim+1][iDim+1] -= ProjRotVel;
						Jacobian_ii[nVar-1][nVar-1] -= ProjRotVel;
					}

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
				}

				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
			}
		}
	}

	delete [] Velocity;
	delete [] UnitaryNormal;
	delete [] Psi;
}

void CAdjEulerSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
		CConfig *config, unsigned short val_marker) {

#ifdef NO_MPI

	unsigned long iVertex, iPoint, jPoint;
	unsigned short iDim;
	double *Psi_i, *Psi_j, *U_i, *U_j, *Coord;

	double  *Normal = new double[nDim];

	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		jPoint = geometry->vertex[val_marker][iVertex]->GetDonorPoint();
		Coord = geometry->node[iPoint]->GetCoord();

		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Adjoint variables w/o reconstruction ---*/
			Psi_i = node[iPoint]->GetSolution();
			Psi_j = node[jPoint]->GetSolution();

			/*--- Conservative variables w/o reconstruction ---*/
			U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
			U_j = solver_container[FLOW_SOL]->node[jPoint]->GetSolution();
			numerics->SetConservative(U_i, U_j);

			/*--- SoundSpeed enthalpy and lambda variables w/o reconstruction ---*/
			numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
					solver_container[FLOW_SOL]->node[jPoint]->GetSoundSpeed());
			numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
					solver_container[FLOW_SOL]->node[jPoint]->GetEnthalpy());

			/*--- Set face vector, and area ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = - Normal[iDim];
			numerics->SetNormal(Normal);

			/*--- Just do a periodic BC ---*/
			numerics->SetAdjointVar(Psi_i, Psi_j);

			/*--- Compute residual ---*/			
			numerics->ComputeResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

			LinSysRes.SubtractBlock(iPoint, Res_Conv_i);

		}
	}

	delete[] Normal;

#else

	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	unsigned long iVertex, iPoint, jPoint;
	unsigned short iVar, iDim;
	double *Adjoint_Var, Psi_i[5], Psi_j[5], *U_i, *U_j;

	double *Normal = new double [nDim]; 
	double *Buffer_Send_Psi = new double[nVar];
	double *Buffer_Receive_Psi = new double[nVar];

	/*--- Do the send process, by the moment we are sending each 
	 node individually, this must be changed ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];

			/*--- We only send the information that belong to other boundary ---*/
			if (jProcessor != rank) {
				Adjoint_Var = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Send_Psi[iVar] = Adjoint_Var[iVar];
				MPI::COMM_WORLD.Bsend(Buffer_Send_Psi, nVar, MPI::DOUBLE, jProcessor, iPoint);
			}
		}
	}


	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];

			/*--- We only receive the information that belong to other boundary ---*/
			if (jProcessor != rank)
				MPI::COMM_WORLD.Recv(Buffer_Receive_Psi, nVar, MPI::DOUBLE, jProcessor, jPoint);
			else {
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Receive_Psi[iVar] = node[jPoint]->GetSolution(iVar); 
			}

			/*--- Store the solution for both points ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Psi_i[iVar] = node[iPoint]->GetSolution(iVar); 
				Psi_j[iVar] = Buffer_Receive_Psi[iVar]; 
			}

			/*--- Conservative variables w/o reconstruction (the same at both points) ---*/
			U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
			U_j = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
			numerics->SetConservative(U_i, U_j);

			/*--- SoundSpeed enthalpy and lambda variables w/o reconstruction (the same at both points) ---*/
			numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
					solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed());
			numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
					solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy());

			/*--- Set face vector, and area ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = - Normal[iDim];
			numerics->SetNormal(Normal);

			/*--- Just do a periodic BC ---*/
			numerics->SetAdjointVar(Psi_i, Psi_j);

			/*--- Compute residual ---*/			
			numerics->ComputeResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
		}
	}

	delete[] Buffer_Send_Psi;
	delete[] Buffer_Receive_Psi;
	delete[] Normal;
#endif

}

void CAdjEulerSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
		CConfig *config, unsigned short val_marker) {

#ifdef NO_MPI

	unsigned long iVertex, iPoint, jPoint, Pin, Pout;
	unsigned short iVar, iDim;
	double  Psi_out[5], Psi_in[5], Psi_out_ghost[5], Psi_in_ghost[5], 
	MeanPsi[5], *Psi_i, *Psi_j, *U_i, *U_j, *IntBoundary_Jump, *Coord;

	double  *Normal = new double[nDim];

	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		jPoint = geometry->vertex[val_marker][iVertex]->GetDonorPoint();
		Coord = geometry->node[iPoint]->GetCoord();

		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Adjoint variables w/o reconstruction ---*/
			Psi_i = node[iPoint]->GetSolution();
			Psi_j = node[jPoint]->GetSolution();

			/*--- Conservative variables w/o reconstruction ---*/
			U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
			U_j = solver_container[FLOW_SOL]->node[jPoint]->GetSolution();
			numerics->SetConservative(U_i, U_j);

			/*--- SoundSpeed enthalpy and lambda variables w/o reconstruction ---*/
			numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
					solver_container[FLOW_SOL]->node[jPoint]->GetSoundSpeed());
			numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
					solver_container[FLOW_SOL]->node[jPoint]->GetEnthalpy());

			/*--- Set face vector, and area ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = - Normal[iDim];
			numerics->SetNormal(Normal);

			/*--- If equivalent area or nearfield pressure condition ---*/
			if ((config->GetKind_ObjFunc() == EQUIVALENT_AREA) || 
					(config->GetKind_ObjFunc() == NEARFIELD_PRESSURE)) {

				if (Normal[nDim-1] < 0.0) { Pin = iPoint; Pout = jPoint; }
				else { Pout = iPoint; Pin = jPoint; }

				for (iVar = 0; iVar < nVar; iVar++) {
					Psi_out[iVar] = node[Pout]->GetSolution(iVar);
					Psi_in[iVar] = node[Pin]->GetSolution(iVar);	
					MeanPsi[iVar] = 0.5*(Psi_out[iVar] + Psi_in[iVar]);
				}

				IntBoundary_Jump = node[iPoint]->GetIntBoundary_Jump();

				/*--- Inner point ---*/
				if (iPoint == Pin) {
					for (iVar = 0; iVar < nVar; iVar++)
						Psi_in_ghost[iVar] = 2.0*MeanPsi[iVar] - Psi_in[iVar] - IntBoundary_Jump[iVar];
					numerics->SetAdjointVar(Psi_in, Psi_in_ghost);
				}

				/*--- Outer point ---*/
				if (iPoint == Pout) {
					for (iVar = 0; iVar < nVar; iVar++)
						Psi_out_ghost[iVar] = 2.0*MeanPsi[iVar] - Psi_out[iVar] + IntBoundary_Jump[iVar];
					numerics->SetAdjointVar(Psi_out, Psi_out_ghost);
				}
			}
			else {
				/*--- Just do a periodic BC ---*/
				numerics->SetAdjointVar(Psi_i, Psi_j);
			}

			/*--- Compute residual ---*/			
			numerics->ComputeResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

			LinSysRes.SubtractBlock(iPoint, Res_Conv_i);

		}
	}

	delete[] Normal;

#else

	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	unsigned long iVertex, iPoint, jPoint, Pin, Pout;
	unsigned short iVar, iDim;
	double *Adjoint_Var, 
	Psi_out[5], Psi_in[5], Psi_i[5], Psi_j[5], Psi_in_ghost[5], Psi_out_ghost[5], MeanPsi[5], *U_i, *U_j, 
	*IntBoundary_Jump;

	double *Normal = new double [nDim]; 
	double *Buffer_Send_Psi = new double[nVar];
	double *Buffer_Receive_Psi = new double[nVar];

	/*--- Do the send process, by the moment we are sending each 
	 node individually, this must be changed ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];

			/*--- We only send the information that belong to other boundary ---*/
			if (jProcessor != rank) {
				Adjoint_Var = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Send_Psi[iVar] = Adjoint_Var[iVar];
				MPI::COMM_WORLD.Bsend(Buffer_Send_Psi, nVar, MPI::DOUBLE, jProcessor, iPoint);
			}
		}
	}


	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];

			/*--- We only receive the information that belong to other boundary ---*/
			if (jProcessor != rank)
				MPI::COMM_WORLD.Recv(Buffer_Receive_Psi, nVar, MPI::DOUBLE, jProcessor, jPoint);
			else {
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Receive_Psi[iVar] = node[jPoint]->GetSolution(iVar); 
			}

			/*--- Store the solution for both points ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Psi_i[iVar] = node[iPoint]->GetSolution(iVar); 
				Psi_j[iVar] = Buffer_Receive_Psi[iVar]; 
			}

			/*--- Conservative variables w/o reconstruction (the same at both points) ---*/
			U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
			U_j = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
			numerics->SetConservative(U_i, U_j);

			/*--- SoundSpeed enthalpy and lambda variables w/o reconstruction (the same at both points) ---*/
			numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
					solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed());
			numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
					solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy());

			/*--- Set face vector, and area ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = - Normal[iDim];
			numerics->SetNormal(Normal);

			/*--- If equivalent area or nearfield pressure condition ---*/
			if ((config->GetKind_ObjFunc() == EQUIVALENT_AREA) || 
					(config->GetKind_ObjFunc() == NEARFIELD_PRESSURE)) {

				/*--- Inner nearfield boundary ---*/
				if (Normal[nDim-1] < 0.0)  { 
					Pin = iPoint; Pout = jPoint;
					for (iVar = 0; iVar < nVar; iVar++) {
						Psi_in[iVar] = Psi_i[iVar];
						Psi_out[iVar] = Psi_j[iVar];
						MeanPsi[iVar] = 0.5*(Psi_out[iVar] + Psi_in[iVar]);
					}
				}
				/*--- Outer nearfield boundary ---*/
				else { 
					Pout = iPoint; Pin = jPoint; 
					for (iVar = 0; iVar < nVar; iVar++) {
						Psi_in[iVar] = Psi_j[iVar];
						Psi_out[iVar] = Psi_i[iVar];
						MeanPsi[iVar] = 0.5*(Psi_out[iVar] + Psi_in[iVar]);
					}
				}

				IntBoundary_Jump = node[iPoint]->GetIntBoundary_Jump();

				/*--- Inner point ---*/
				if (iPoint == Pin) {
					for (iVar = 0; iVar < nVar; iVar++)
						Psi_in_ghost[iVar] = 2.0*MeanPsi[iVar] - Psi_in[iVar] - IntBoundary_Jump[iVar];
					numerics->SetAdjointVar(Psi_in, Psi_in_ghost);
				}

				/*--- Outer point ---*/
				if (iPoint == Pout) {
					for (iVar = 0; iVar < nVar; iVar++)
						Psi_out_ghost[iVar] = 2.0*MeanPsi[iVar] - Psi_out[iVar] + IntBoundary_Jump[iVar];
					numerics->SetAdjointVar(Psi_out, Psi_out_ghost);	
				}
			}
			else {
				/*--- Just do a periodic BC ---*/
				numerics->SetAdjointVar(Psi_i, Psi_j);
			}

			/*--- Compute residual ---*/			
			numerics->ComputeResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
		}
	}

	delete[] Buffer_Send_Psi;
	delete[] Buffer_Receive_Psi;
	delete[] Normal;
#endif	
}

void CAdjEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
		CConfig *config, unsigned short val_marker) {

	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	double *Normal, *U_domain, *U_infty, *Psi_domain, *Psi_infty;

	bool rotating_frame = config->GetRotating_Frame();
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
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
			if (incompressible) {
				U_infty[0] = solver_container[FLOW_SOL]->GetPressure_Inf();
				U_infty[1] = solver_container[FLOW_SOL]->GetVelocity_Inf(0)*config->GetDensity_FreeStreamND();
				U_infty[2] = solver_container[FLOW_SOL]->GetVelocity_Inf(1)*config->GetDensity_FreeStreamND();
				if (nDim == 3) U_infty[3] = solver_container[FLOW_SOL]->GetVelocity_Inf(2)*config->GetDensity_FreeStreamND();
			}
			else {
				/*--- Flow Solution at infinity ---*/
				U_infty[0] = solver_container[FLOW_SOL]->GetDensity_Inf();
				U_infty[1] = solver_container[FLOW_SOL]->GetDensity_Velocity_Inf(0);
				U_infty[2] = solver_container[FLOW_SOL]->GetDensity_Velocity_Inf(1);
				U_infty[3] = solver_container[FLOW_SOL]->GetDensity_Energy_Inf();
				if (nDim == 3) {
					U_infty[3] = solver_container[FLOW_SOL]->GetDensity_Velocity_Inf(2);
					U_infty[4] = solver_container[FLOW_SOL]->GetDensity_Energy_Inf();
				}
			}
			conv_numerics->SetConservative(U_domain, U_infty);

			if(config->GetKind_Adjoint() != DISCRETE) {
				/*--- Adjoint flow solution at the wall ---*/
				for (iVar = 0; iVar < nVar; iVar++) {
					Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
					Psi_infty[iVar] = 0.0;
				}
				conv_numerics->SetAdjointVar(Psi_domain, Psi_infty);
			}

			if (incompressible) {
				conv_numerics->SetDensityInc(config->GetDensity_FreeStreamND(), config->GetDensity_FreeStreamND());
				conv_numerics->SetBetaInc2(solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2(), 
						solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2());
				conv_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
			}
			else {		
				conv_numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
						solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed());
				conv_numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
						solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy());
			}

			/*--- Rotating Frame ---*/
			if (rotating_frame) {
				conv_numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
				conv_numerics->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}

			/*--- Grid Movement ---*/
			if (grid_movement)
				conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

			/*--- Compute the upwind flux ---*/
			if (config->GetKind_Adjoint() == DISCRETE)
				conv_numerics->ComputeResidual(Jacobian_i, Jacobian_j, config);
			else
				conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

			/*--- Add and Subtract Residual ---*/
			if(config->GetKind_Adjoint() == DISCRETE) {
				//Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

			}
			else {
				LinSysRes.SubtractBlock(iPoint, Residual_i);

				/*--- Implicit contribution to the residual ---*/
				if ((implicit) && (config->GetKind_Adjoint() != DISCRETE))
					Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
			}
		}
	}

	delete [] Normal;
	delete [] U_domain; delete [] U_infty;
	delete [] Psi_domain; delete [] Psi_infty;
}

void CAdjEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iVar, iDim, Kind_Inlet = config->GetKind_Inlet();
	unsigned long iVertex, iPoint, Point_Normal;
	double P_Total, T_Total, Velocity[3], Density_Inlet, Velocity2, H_Total,
  Temperature, Riemann, Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2,
  SoundSpeed_Total2, Vel_Mag, alpha, aa, bb, cc, dd, bcn, phin, Area, UnitaryNormal[3],
  ProjGridVel, *GridVel, ProjRotVel, *RotVel;
  
	double Two_Gamma_M1 = 2.0/Gamma_Minus_One;
	double Gas_Constant = config->GetGas_ConstantND();
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
	bool grid_movement = config->GetGrid_Movement();
	bool rotating_frame = config->GetRotating_Frame();
    bool freesurface = config->GetFreeSurface();
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);

	double *Normal = new double[nDim];
	double *U_domain   = new double[nVar]; double *U_inlet   = new double[nVar];
	double *Psi_domain = new double[nVar]; double *Psi_inlet = new double[nVar];

	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check that the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
			conv_numerics->SetNormal(Normal);

			Area = 0.0; for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;

			/*--- Set the normal point ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

			/*--- Flow solution at the boundary ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);

			/*--- Adjoint flow solution at the boundary ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Construct the flow & adjoint states at the inlet ---*/
			if (incompressible) {

//        /*--- Flow solution at the boundary ---*/
//        Density_Inlet = solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
//        for (iVar = 0; iVar < nVar; iVar++)
//          U_inlet[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
        
				/*--- Pressure and density using the internal value ---*/
				U_inlet[0] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(0);
				Density_Inlet = solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
        
				/*--- The velocity is computed from the infinity values ---*/
				for (iDim = 0; iDim < nDim; iDim++)
					U_inlet[iDim+1] = solver_container[FLOW_SOL]->GetVelocity_Inf(iDim)*Density_Inlet;
        
				/*--- The y/z velocity is interpolated due to the
         free surface effect on the pressure ---*/
				if (freesurface) U_inlet[nDim] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(nDim);
        
				/*--- Adjoint solution at the inlet ---*/
				Psi_inlet[0] = node[iPoint]->GetSolution(0);
        for (iDim = 0; iDim < nDim; iDim++)
          Psi_inlet[iDim+1] = 0.0;

			}	else {

				/*--- Subsonic, compressible inflow: first build the flow state
         using the same method as the direct problem. Then, based on
         those conservative values, compute the characteristic-based 
         adjoint boundary condition. The boundary update to be applied
         depends on whether total conditions or mass flow are specified. ---*/

				switch (Kind_Inlet) {

				/*--- Total properties have been specified at the inlet. ---*/
				case TOTAL_CONDITIONS:

					/*--- Retrieve the specified total conditions for this inlet. ---*/
					P_Total  = config->GetInlet_Ptotal(Marker_Tag);
					T_Total  = config->GetInlet_Ttotal(Marker_Tag);
					Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

					/*--- Non-dim. the inputs if necessary. ---*/
					P_Total /= config->GetPressure_Ref();
					T_Total /= config->GetTemperature_Ref();

					/*--- Store primitives and set some variables for clarity. ---*/
					Density = U_domain[0];
					Velocity2 = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						Velocity[iDim] = U_domain[iDim+1]/Density;
						Velocity2 += Velocity[iDim]*Velocity[iDim];
					}
					Energy      = U_domain[nVar-1]/Density;
					Pressure    = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
					H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
					SoundSpeed2 = Gamma*Pressure/Density;

					/*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/
					Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
					for (iDim = 0; iDim < nDim; iDim++)
						Riemann += Velocity[iDim]*UnitaryNormal[iDim];

					/*--- Total speed of sound ---*/
					SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy
							+ Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

					/*--- Dot product of normal and flow direction. This should
             be negative due to outward facing boundary normal convention. ---*/
					alpha = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						alpha += UnitaryNormal[iDim]*Flow_Dir[iDim];

					/*--- Coefficients in the quadratic equation for the velocity ---*/
					aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
					bb = -1.0*Gamma_Minus_One*alpha*Riemann;
					cc =  0.5*Gamma_Minus_One*Riemann*Riemann
							-2.0*SoundSpeed_Total2/Gamma_Minus_One;

					/*--- Solve quadratic equation for velocity magnitude. Value must
             be positive, so the choice of root is clear. ---*/
					dd = bb*bb - 4.0*aa*cc;
					dd = sqrt(max(0.0,dd));
					Vel_Mag   = (-bb + dd)/(2.0*aa);
					Vel_Mag   = max(0.0,Vel_Mag);
					Velocity2 = Vel_Mag*Vel_Mag;

					/*--- Compute speed of sound from total speed of sound eqn. ---*/
					SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

					/*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
					Mach2 = Velocity2/SoundSpeed2;
					Mach2 = min(1.0,Mach2);
					Velocity2   = Mach2*SoundSpeed2;
					Vel_Mag     = sqrt(Velocity2);
					SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

					/*--- Compute new velocity vector at the inlet ---*/
					for (iDim = 0; iDim < nDim; iDim++)
						Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

					/*--- Static temperature from the speed of sound relation ---*/
					Temperature = SoundSpeed2/(Gamma*Gas_Constant);

					/*--- Static pressure using isentropic relation at a point ---*/
					Pressure = P_Total*pow((Temperature/T_Total),Gamma/Gamma_Minus_One);

					/*--- Density at the inlet from the gas law ---*/
					Density = Pressure/(Gas_Constant*Temperature);

					/*--- Using pressure, density, & velocity, compute the energy ---*/
					Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;

					/*--- Conservative variables, using the derived quantities ---*/
					U_inlet[0] = Density;
					U_inlet[1] = Velocity[0]*Density;
					U_inlet[2] = Velocity[1]*Density;
					U_inlet[3] = Energy*Density;
					if (nDim == 3) {
						U_inlet[3] = Velocity[2]*Density;
						U_inlet[4] = Energy*Density;
					}

					/*--- Adjoint solution at the inlet. Set to zero for now
             but should be replaced with derived expression for this type of
             inlet. ---*/
					for (iVar = 0; iVar < nVar; iVar++)
						Psi_inlet[iVar] = 0.0;

					break;

					/*--- Mass flow has been specified at the inlet. ---*/
				case MASS_FLOW:

					/*--- Retrieve the specified mass flow for the inlet. ---*/
					Density  = config->GetInlet_Ttotal(Marker_Tag);
					Vel_Mag  = config->GetInlet_Ptotal(Marker_Tag);
					Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

					/*--- Non-dim. the inputs if necessary. ---*/
					Density /= config->GetDensity_Ref();
					Vel_Mag /= config->GetVelocity_Ref();

					/*--- Get primitives from current inlet state. ---*/
					for (iDim = 0; iDim < nDim; iDim++)
						Velocity[iDim] = solver_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim, incompressible);
					Pressure    = solver_container[FLOW_SOL]->node[iPoint]->GetPressure(incompressible);
					SoundSpeed2 = Gamma*Pressure/U_domain[0];

					/*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/
					Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
					for (iDim = 0; iDim < nDim; iDim++)
						Riemann += Velocity[iDim]*UnitaryNormal[iDim];

					/*--- Speed of sound squared for fictitious inlet state ---*/
					SoundSpeed2 = Riemann;
					for (iDim = 0; iDim < nDim; iDim++)
						SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitaryNormal[iDim];

					SoundSpeed2 = max(0.0,0.5*Gamma_Minus_One*SoundSpeed2);
					SoundSpeed2 = SoundSpeed2*SoundSpeed2;

					/*--- Pressure for the fictitious inlet state ---*/
					Pressure = SoundSpeed2*Density/Gamma;

					/*--- Energy for the fictitious inlet state ---*/
					Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Vel_Mag*Vel_Mag;

					/*--- Conservative variables, using the derived quantities ---*/
					U_inlet[0] = Density;
					U_inlet[1] = Vel_Mag*Flow_Dir[0]*Density;
					U_inlet[2] = Vel_Mag*Flow_Dir[1]*Density;
					U_inlet[3] = Energy*Density;
					if (nDim == 3) {
						U_inlet[3] = Vel_Mag*Flow_Dir[2]*Density;
						U_inlet[4] = Energy*Density;
					}

					/*--- Retrieve current adjoint solution values at the boundary ---*/
					for (iVar = 0; iVar < nVar; iVar++)
						Psi_inlet[iVar] = node[iPoint]->GetSolution(iVar);

					/*--- Some terms needed for the adjoint BC ---*/
					bcn = 0.0; phin = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						bcn  -= (Gamma/Gamma_Minus_One)*Velocity[iDim]*UnitaryNormal[iDim];
						phin += Psi_domain[iDim+1]*UnitaryNormal[iDim];
					}

					/*--- Extra boundary term for a rotating frame ---*/
					if (rotating_frame) {
						ProjRotVel = 0.0;
						RotVel = geometry->node[iPoint]->GetRotVel();
						for (iDim = 0; iDim < nDim; iDim++) {
							ProjRotVel += RotVel[iDim]*UnitaryNormal[iDim];
						}
						ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux()/Area;
						bcn -= (1.0/Gamma_Minus_One)*ProjRotVel;
					}

					/*--- Extra boundary term for grid movement ---*/
					if (grid_movement) {
						ProjGridVel = 0.0;
						GridVel = geometry->node[iPoint]->GetGridVel();
						for (iDim = 0; iDim < nDim; iDim++)
							ProjGridVel += GridVel[iDim]*UnitaryNormal[iDim];
						bcn -= (1.0/Gamma_Minus_One)*ProjGridVel;
					}

					/*--- Impose value for PsiE based on hand-derived expression. ---*/
					Psi_inlet[nVar-1] = -phin*(1.0/bcn);

					break;
				}
			}

			/*--- Set the flow and adjoint states in the solver ---*/
			conv_numerics->SetConservative(U_domain, U_inlet);
			conv_numerics->SetAdjointVar(Psi_domain, Psi_inlet);

			if (incompressible) {
				conv_numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), Density_Inlet);
				conv_numerics->SetBetaInc2(solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2(),
                                 solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2());
				conv_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
			}
			else {		
				conv_numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
						solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed());
				conv_numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
						solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy());
			}

			/*--- Rotational frame ---*/
			if (rotating_frame) {
				conv_numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(),
						geometry->node[iPoint]->GetRotVel());
				conv_numerics->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}

			/*--- Grid Movement ---*/
			if (grid_movement)
				conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
						geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
			conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
                               Jacobian_ji, Jacobian_jj, config);

			/*--- Add and Subtract Residual ---*/
			LinSysRes.SubtractBlock(iPoint, Residual_i);

			/*--- Implicit contribution to the residual ---*/
			if (implicit) 
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}

	/*--- Free locally allocated memory ---*/
	delete [] Normal;
	delete [] U_domain;   delete [] U_inlet;
	delete [] Psi_domain; delete [] Psi_inlet;

}

void CAdjEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

	/*--- Local variables and initialization. ---*/
	unsigned short iVar, iDim;

	unsigned long iVertex, iPoint, Point_Normal;

	double Pressure, P_Exit, Velocity[3], Velocity2, Entropy;
	double Density, Energy, Riemann, Height;
	double Vn, SoundSpeed, Mach_Exit, Vn_Exit, Ubn, a1, LevelSet, Density_Outlet;
	double *U_domain = new double[nVar]; double *U_outlet = new double[nVar];
	double *Psi_domain = new double [nVar]; double *Psi_outlet = new double [nVar];
	double *Normal = new double[nDim];

	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();
    double FreeSurface_Zero = config->GetFreeSurface_Zero();
	double PressFreeSurface = solver_container[FLOW_SOL]->GetPressure_Inf();
    double epsilon          = config->GetFreeSurface_Thickness();
    double RatioDensity     = config->GetRatioDensity();
    double Froude           = config->GetFroude();
    bool freesurface = config->GetFreeSurface();

  
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);

	/*--- Loop over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

			double Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;

			/*--- Set the normal point ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

			/*--- Flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);

			/*--- Adjoint flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Construct the flow & adjoint states at the outlet ---*/
			if (incompressible) {

        if (freesurface) {
          
//          /*--- Flow solution at the boundary ---*/
//          Density_Outlet = solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
//          for (iVar = 0; iVar < nVar; iVar++)
//            U_outlet[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
          
					/*--- Density computation at the exit using the level set function ---*/
					Height = geometry->node[iPoint]->GetCoord(nDim-1);
					LevelSet = Height - FreeSurface_Zero;
          
					/*--- Pressure computation the density at the exit (imposed) ---*/
					if (LevelSet < -epsilon) Density_Outlet = config->GetDensity_FreeStreamND();
					if (LevelSet > epsilon) Density_Outlet = RatioDensity*config->GetDensity_FreeStreamND();
					U_outlet[0] = PressFreeSurface + Density_Outlet*((FreeSurface_Zero-Height)/(Froude*Froude));
          
					/*--- Neumman condition in the interface for the pressure and density ---*/
					if (fabs(LevelSet) <= epsilon) {
						U_outlet[0] = solver_container[FLOW_SOL]->node[Point_Normal]->GetSolution(0);
						Density_Outlet = solver_container[FLOW_SOL]->node[Point_Normal]->GetDensityInc();
					}
          
        }
        
        else {
          
					/*--- Imposed pressure and density ---*/
					Density_Outlet = solver_container[FLOW_SOL]->GetDensity_Inf();
					U_outlet[0] = solver_container[FLOW_SOL]->GetPressure_Inf();
          
				}
          
        /*--- Neumman condition for the velocity ---*/
				for (iDim = 0; iDim < nDim; iDim++)
					U_outlet[iDim+1] = node[Point_Normal]->GetSolution(iDim+1);

				/*--- Adjoint flow solution at the outlet (hard-coded for size[3] again?) ---*/
				Psi_outlet[2] = 0.0;
				double coeff = (2.0*U_domain[1])/ solver_container[FLOW_SOL]->node[Point_Normal]->GetBetaInc2();
				Psi_outlet[1] = node[Point_Normal]->GetSolution(1);
				Psi_outlet[0] = -coeff*Psi_outlet[1];
        
//				/*--- Adjoint solution at the inlet ---*/
//        for (iVar = 0; iVar < nVar; iVar++)
//          Psi_outlet[iVar] = 0.0;
        
        
			} else {

				/*--- Retrieve the specified back pressure for this outlet. ---*/
				P_Exit = config->GetOutlet_Pressure(Marker_Tag);

				/*--- Non-dim. the inputs if necessary. ---*/
				P_Exit = P_Exit/config->GetPressure_Ref();

				/*--- Check whether the flow is supersonic at the exit. The type
         of boundary update depends on this. ---*/
				Density = U_domain[0];
				Velocity2 = 0.0; Vn = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Velocity[iDim] = U_domain[iDim+1]/Density;
					Velocity2 += Velocity[iDim]*Velocity[iDim];
					Vn += Velocity[iDim]*UnitaryNormal[iDim];
				}
				Energy     = U_domain[nVar-1]/Density;
				Pressure   = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
				SoundSpeed = sqrt(Gamma*Pressure/Density);
				Mach_Exit  = sqrt(Velocity2)/SoundSpeed;

				if (Mach_Exit >= 1.0) {

					/*--- Supersonic exit flow: there are no incoming characteristics,
           so no boundary condition is necessary. Set outlet state to current
           state so that upwinding handles the direction of propagation. This
           means that all variables can be imposed for the adjoint problem,
           so set them all to zero to remove contributions from the boundary
           integral in the adjoint formulation (impose orthogonality). ---*/

					for (iVar = 0; iVar < nVar; iVar++) {
						U_outlet[iVar] = U_domain[iVar];
						Psi_outlet[iVar] = 0.0;
					}

				} else {

					/*--- Subsonic exit flow: there is one incoming characteristic (u-c),
           therefore one variable can be specified (back pressure) and is used
           to update the conservative variables. Compute the entropy and the
           acoustic Riemann variable. These invariants, as well as the
           tangential velocity components, are extrapolated. Adapted from an
           original implementation in the Stanford University multi-block
           (SUmb) solver in the routine bcSubsonicOutflow.f90 by Edwin van
           der Weide, last modified 09-10-2007. ---*/

					Entropy = Pressure*pow(1.0/Density,Gamma);
					Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

					/*--- Compute the new fictitious state at the outlet ---*/
					Density    = pow(P_Exit/Entropy,1.0/Gamma);
					Pressure   = P_Exit;
					SoundSpeed = sqrt(Gamma*P_Exit/Density);
					Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
					Velocity2  = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitaryNormal[iDim];
						Velocity2 += Velocity[iDim]*Velocity[iDim];
					}
					Energy  = P_Exit/(Density*Gamma_Minus_One) + 0.5*Velocity2;

					/*--- Conservative variables, using the derived quantities ---*/
					U_outlet[0] = Density;
					U_outlet[1] = Velocity[0]*Density;
					U_outlet[2] = Velocity[1]*Density;
					U_outlet[3] = Energy*Density;
					if (nDim == 3) {
						U_outlet[3] = Velocity[2]*Density;
						U_outlet[4] = Energy*Density;
					}

					/*--- One condition is imposed at the exit (back pressure) for the
           flow problem, so nVar-1 conditions are imposed at the outlet on
           the adjoint variables. Choose PsiE as the free variable and compute 
           PsiRho & Phi using hand-derived expressions. ---*/

					/*--- Compute (Vn - Ubn).n term for use in the BC. ---*/
					Vn = 0.0; Ubn = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Vn += Velocity[iDim]*UnitaryNormal[iDim];

					/*--- Extra boundary term for a rotating frame ---*/
					if (rotating_frame) {
						double ProjRotVel = 0.0;
						double *RotVel = geometry->node[iPoint]->GetRotVel();
						for (iDim = 0; iDim < nDim; iDim++) {
							ProjRotVel += RotVel[iDim]*UnitaryNormal[iDim];
						}
						ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux()/Area;
						Ubn = ProjRotVel;
					}

					/*--- Extra boundary term for grid movement ---*/
					if (grid_movement) {
						double ProjGridVel = 0.0;
						double *GridVel = geometry->node[iPoint]->GetGridVel();
						for (iDim = 0; iDim < nDim; iDim++)
							ProjGridVel += GridVel[iDim]*UnitaryNormal[iDim];
						Ubn = ProjGridVel;
					}

					/*--- Shorthand for repeated term in the boundary conditions ---*/
					a1 = Gamma*(P_Exit/(Density*Gamma_Minus_One))/(Vn-Ubn);

					/*--- Impose values for PsiRho & Phi using PsiE from domain. ---*/
					Psi_outlet[nVar-1] = Psi_domain[nVar-1];
					Psi_outlet[0] = 0.5*Psi_outlet[nVar-1]*Velocity2;
					for (iDim = 0; iDim < nDim; iDim++) {
						Psi_outlet[0]   += Psi_outlet[nVar-1]*a1*Velocity[iDim]*UnitaryNormal[iDim];
						Psi_outlet[iDim+1] = -Psi_outlet[nVar-1]*(a1*UnitaryNormal[iDim] + Velocity[iDim]);
					}

					//          /*--- Giles & Pierce ---*/
					//          Psi_outlet[nVar-1] = Psi_domain[nVar-1];
					//          Psi_outlet[0] = 0.5*Psi_outlet[nVar-1]*Velocity2;
					//          for (iDim = 0; iDim < nDim; iDim++) {
					//            Psi_outlet[0]   += Psi_outlet[nVar-1]*(SoundSpeed/Gamma_Minus_One)*Velocity[iDim]*UnitaryNormal[iDim];
					//            Psi_outlet[iDim+1] = -Psi_outlet[nVar-1]*((SoundSpeed/Gamma_Minus_One)*UnitaryNormal[iDim] + Velocity[iDim]);
					//          }


				}
			}

			/*--- Set the flow and adjoint states in the solver ---*/
			conv_numerics->SetNormal(Normal);
			conv_numerics->SetConservative(U_domain, U_outlet);
			conv_numerics->SetAdjointVar(Psi_domain, Psi_outlet);

			if (incompressible) {
				conv_numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), Density_Outlet);
				conv_numerics->SetBetaInc2(solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2(),
						solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2());
				conv_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
						geometry->node[iPoint]->GetCoord());
			}
			else {
				conv_numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(),
						solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed());
				conv_numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(),
						solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy());
			}

			/*--- Rotational frame ---*/
			if (rotating_frame) {
				conv_numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(),
						geometry->node[iPoint]->GetRotVel());
				conv_numerics->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}

			/*--- Grid Movement ---*/
			if (grid_movement)
				conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
						geometry->node[iPoint]->GetGridVel());

			conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
					Jacobian_ji, Jacobian_jj, config);

			/*--- Add and Subtract Residual ---*/
			LinSysRes.SubtractBlock(iPoint, Residual_i);

			/*--- Implicit contribution to the residual ---*/
			if (implicit)
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}

	/*--- Free locally allocated memory ---*/
	delete [] Normal;
	delete [] U_domain; delete [] U_outlet;
	delete [] Psi_domain; delete [] Psi_outlet;

}

void CAdjEulerSolver::BC_Nacelle_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	double *Normal, *U_domain, *U_inflow, *Psi_domain, *Psi_inflow, Pressure, P_Fan, Velocity[3], Velocity2, Entropy, Density, Energy, Riemann, Enthalpy, Vn, SoundSpeed, Mach_Exit, Vn_Exit, UnitaryNormal[3], Area, a1;
	unsigned short iVar, iDim;
	unsigned long iVertex, iPoint;

	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);

	Normal = new double[nDim];
	U_domain = new double[nVar]; U_inflow = new double[nVar];
	Psi_domain = new double[nVar]; Psi_inflow = new double[nVar];

	/*--- Loop over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

			Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;

			/*--- Current solution at this boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);

			/*--- Retrieve the specified back pressure for this outlet. ---*/
			P_Fan = config->GetFanFace_Pressure(Marker_Tag);

			/*--- Check whether the flow is supersonic at the exit. The type
			 of boundary update depends on this. ---*/
			Density = U_domain[0];
			Velocity2 = 0.0; Vn = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Velocity[iDim] = U_domain[iDim+1]/Density;
				Velocity2 += Velocity[iDim]*Velocity[iDim];
				Vn += Velocity[iDim]*UnitaryNormal[iDim];
			}
			Energy     = U_domain[nVar-1]/Density;
			Pressure   = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
			SoundSpeed = sqrt(Gamma*Pressure/Density);
			Mach_Exit  = sqrt(Velocity2)/SoundSpeed;

			/*--- Subsonic exit flow: there is one incoming characteristic,
			 therefore one variable can be specified (back pressure) and is used
			 to update the conservative variables. Compute the entropy and the
			 acoustic variable. These riemann invariants, as well as the tangential
			 velocity components, are extrapolated. ---*/
			Entropy = Pressure*pow(1.0/Density,Gamma);
			Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

			/*--- Compute the new fictious state at the outlet ---*/
			Density    = pow(P_Fan/Entropy,1.0/Gamma);
			Pressure   = P_Fan;
			SoundSpeed = sqrt(Gamma*P_Fan/Density);
			Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
			Velocity2  = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitaryNormal[iDim];
				Velocity2 += Velocity[iDim]*Velocity[iDim];
			}
			Energy  = P_Fan/(Density*Gamma_Minus_One) + 0.5*Velocity2;
			Enthalpy = (Energy*Density + Pressure) / Density;

			/*--- Conservative variables, using the derived quantities ---*/
			U_inflow[0] = Density;
			U_inflow[1] = Velocity[0]*Density;
			U_inflow[2] = Velocity[1]*Density;
			U_inflow[3] = Energy*Density;
			if (nDim == 3) {
				U_inflow[3] = Velocity[2]*Density;
				U_inflow[4] = Energy*Density;
			}

			conv_numerics->SetConservative(U_domain, U_inflow);
			conv_numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), SoundSpeed);		
			conv_numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), Enthalpy);

			/*--- Adjoint flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
      
      /*--- Shorthand for repeated term in the boundary conditions ---*/
      a1 = Gamma*(P_Fan/(Density*Gamma_Minus_One))/(Vn);
      
      /*--- Impose values for PsiRho & Phi using PsiE from domain. ---*/
      Psi_inflow[nVar-1] = Psi_domain[nVar-1];
      Psi_inflow[0] = 0.5*Psi_inflow[nVar-1]*Velocity2;
      for (iDim = 0; iDim < nDim; iDim++) {
        Psi_inflow[0]   += Psi_inflow[nVar-1]*a1*Velocity[iDim]*UnitaryNormal[iDim];
        Psi_inflow[iDim+1] = -Psi_inflow[nVar-1]*(a1*UnitaryNormal[iDim] + Velocity[iDim]);
      }
      
			conv_numerics->SetAdjointVar(Psi_domain, Psi_inflow);

			/*--- Set the normal vector ---*/
			conv_numerics->SetNormal(Normal);

			/*--- Compute the residual ---*/
			conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

			/*--- Add and Subtract Residual ---*/
			LinSysRes.SubtractBlock(iPoint, Residual_i);

			/*--- Implicit contribution to the residual ---*/
			if (implicit) 
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}

	delete [] Normal;
	delete [] U_domain; delete [] U_inflow;
	delete [] Psi_domain; delete [] Psi_inflow;

}

void CAdjEulerSolver::BC_Nacelle_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, Point_Normal;
	double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann, Enthalpy, Area, UnitaryNormal[3], Pressure, Density, Energy, Mach2, SoundSpeed2, SoundSpeed_Total2, SoundSpeed, Vel_Mag, alpha, aa, bb, cc, dd, Flow_Dir[3], *Normal, *U_domain, *U_exhaust, *Psi_domain, *Psi_exhaust;
	unsigned short iVar, iDim;

	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
	double Gas_Constant = config->GetGas_ConstantND();

	Normal = new double[nDim];
	U_domain = new double[nVar]; U_exhaust = new double[nVar];
	Psi_domain = new double[nVar]; Psi_exhaust = new double[nVar];

	/*--- Loop over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

			Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;

			/*--- Current solution at this boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
      
			/*--- Subsonic inflow: there is one outgoing characteristic (u-c),
			 therefore we can specify all but one state variable at the inlet.
			 The outgoing Riemann invariant provides the final piece of info. ---*/

			/*--- Retrieve the specified total conditions for this inlet. ---*/
			P_Total  = config->GetNozzle_Ptotal(Marker_Tag);
			T_Total  = config->GetNozzle_Ttotal(Marker_Tag);

			/*--- Non-dim. the inputs if necessary. ---*/
			P_Total /= config->GetPressure_Ref();
			T_Total /= config->GetTemperature_Ref();

			/*--- Store primitives and set some variables for clarity. ---*/
			Density = U_domain[0];
			Velocity2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Velocity[iDim] = U_domain[iDim+1]/Density;
				Velocity2 += Velocity[iDim]*Velocity[iDim];
			}
			Energy      = U_domain[nVar-1]/Density;
			Pressure    = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
			H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
			SoundSpeed2 = Gamma*Pressure/Density;

			/*--- Compute the acoustic Riemann invariant that is extrapolated
			 from the domain interior. ---*/
			Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
			for (iDim = 0; iDim < nDim; iDim++)
				Riemann += Velocity[iDim]*UnitaryNormal[iDim];

			/*--- Total speed of sound ---*/
			SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

			/*--- The flow direction is defined by the surface normal ---*/
			for (iDim = 0; iDim < nDim; iDim++)
				Flow_Dir[iDim] = -UnitaryNormal[iDim];

			/*--- Dot product of normal and flow direction. This should
			 be negative due to outward facing boundary normal convention. ---*/
			alpha = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				alpha += UnitaryNormal[iDim]*Flow_Dir[iDim];

			/*--- Coefficients in the quadratic equation for the velocity ---*/
			aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
			bb = -1.0*Gamma_Minus_One*alpha*Riemann;
			cc =  0.5*Gamma_Minus_One*Riemann*Riemann -2.0*SoundSpeed_Total2/Gamma_Minus_One;

			/*--- Solve quadratic equation for velocity magnitude. Value must
			 be positive, so the choice of root is clear. ---*/
			dd = bb*bb - 4.0*aa*cc;
			dd = sqrt(max(0.0,dd));
			Vel_Mag   = (-bb + dd)/(2.0*aa);
			Vel_Mag   = max(0.0,Vel_Mag);
			Velocity2 = Vel_Mag*Vel_Mag;

			/*--- Compute speed of sound from total speed of sound eqn. ---*/
			SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

			/*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
			Mach2 = Velocity2/SoundSpeed2;
			Mach2 = min(1.0,Mach2);
			Velocity2   = Mach2*SoundSpeed2;
			Vel_Mag     = sqrt(Velocity2);
			SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
			SoundSpeed = sqrt(SoundSpeed2);

			/*--- Compute new velocity vector at the inlet ---*/
			for (iDim = 0; iDim < nDim; iDim++)
				Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

			/*--- Static temperature from the speed of sound relation ---*/
			Temperature = SoundSpeed2/(Gamma*Gas_Constant);

			/*--- Static pressure using isentropic relation at a point ---*/
			Pressure = P_Total*pow((Temperature/T_Total),Gamma/Gamma_Minus_One);

			/*--- Density at the inlet from the gas law ---*/
			Density = Pressure/(Gas_Constant*Temperature);

			/*--- Using pressure, density, & velocity, compute the energy ---*/
			Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;

			/*--- Conservative variables, using the derived quantities ---*/
			U_exhaust[0] = Density;
      for (iDim = 0; iDim < nDim; iDim++)
        U_exhaust[iDim+1] = Velocity[iDim]*Density;
			U_exhaust[nDim+1] = Energy*Density;

			conv_numerics->SetConservative(U_domain, U_exhaust);
			conv_numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), SoundSpeed);
			conv_numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), Enthalpy);

			/*--- Adjoint flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Adjoint flow solution at the exhaust (this should be improved using characteristics bc) ---*/
      Psi_exhaust[0] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
				Psi_exhaust[iDim+1] = node[Point_Normal]->GetSolution(iDim+1);
      Psi_exhaust[nDim+1] = 0.0;
      
			conv_numerics->SetAdjointVar(Psi_domain, Psi_exhaust);

			/*--- Set the normal vector ---*/
			conv_numerics->SetNormal(Normal);

			/*--- Compute the residual ---*/
			conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

			/*--- Add and Subtract Residual ---*/
			LinSysRes.SubtractBlock(iPoint, Residual_i);

			/*--- Implicit contribution to the residual ---*/
			if (implicit)
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}

	delete [] Normal;
	delete [] U_domain; delete [] U_exhaust;
	delete [] Psi_domain; delete [] Psi_exhaust;

}

void CAdjEulerSolver::BC_FWH(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

	/*--- Dirichlet BC to set the solution from the adjoint coupling terms ---*/

	unsigned long iPoint, iVertex, total_index;
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Store solution computed in coupling routine ---*/
		node[iPoint]->SetSolution(node[iPoint]->GetIntBoundary_Jump());
		node[iPoint]->SetSolution_Old(node[iPoint]->GetIntBoundary_Jump());
		LinSysRes.SetBlock_Zero(iPoint);
		node[iPoint]->SetRes_TruncErrorZero();

		/*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
		if (implicit)
			for (unsigned short iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				Jacobian.DeleteValsRowi(total_index);
			}
	}

}

void CAdjEulerSolver::SetAeroacoustic_Coupling(CSolver ***wave_solution, CSolver ***flow_solution, CNumerics *numerics, CGeometry **flow_geometry, CConfig *flow_config) {

	/*--- Local variables and initialization ---*/

	unsigned short iMarker, iVar, jVar, kVar, iDim;
	unsigned short jc, jrjc, jrjcm1, jrjcp1, jr, jm, jrm1, jrjr, jrp1, jmjm;
	unsigned long iVertex, iPoint;
	double aux, *coord, u, v, w = 0.0, sq_vel, E = 0.0;
	double *U_i, M[5][5], AM[5][5], b[5], sum, rho;
	double *Phi = NULL, *U_i_old = NULL, *Normal = NULL;
	double delta_T = flow_config->GetDelta_UnstTimeND();
	double *Psi_New = new double [nVar];
	double *Velocity;
	Velocity= new double[nDim];
	double **A;
	A = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		A[iVar] = new double[nVar];
	ifstream index_file;
	string text_line;


	/*--- Compute the value of the adjoint variables for the coupled adjoint.
   This requires solving a small system at every node on the FWH surface. ---*/

	for (iMarker = 0; iMarker < flow_config->GetnMarker_All(); iMarker++)
		if (flow_config->GetMarker_All_Boundary(iMarker) == FWH_SURFACE)
			for(iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
				iPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();

				if (flow_geometry[MESH_0]->node[iPoint]->GetDomain()) {

					/*--- Some geometry information for this boundary node ---*/
					coord = flow_geometry[MESH_0]->node[iPoint]->GetCoord();
					Normal = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNormal();
					double Area = 0.0; double UnitaryNormal[3];
					for (iDim = 0; iDim < nDim; iDim++)
						Area += Normal[iDim]*Normal[iDim];
					Area = sqrt (Area);

					/*--- Flip sign for the adjoint - should be opposite of the direct solution? ---*/
					for (iDim = 0; iDim < nDim; iDim++)
						UnitaryNormal[iDim] = Normal[iDim]/Area;

					/*--- Direct solution at this point for building the inviscid Jacobian --*/
					U_i = flow_solution[MESH_0][FLOW_SOL]->node[iPoint]->GetSolution();
					u = U_i[1]/U_i[0]; v = U_i[2]/U_i[0];
					Velocity[0] = U_i[1]/U_i[0]; Velocity[1] = U_i[2]/U_i[0];
					rho = U_i[0];
					sq_vel = u*u+v*v;
					if (nDim == 2)	E = U_i[3]/U_i[0];
					if (nDim == 3) { w = U_i[3]/U_i[0]; E = U_i[4]/U_i[0]; }

					if (nDim == 2) {

						/*--- Build matrix for inviscid Jacobian projected in the local normal direction ---*/
						numerics->GetInviscidProjJac(Velocity, &E, Normal, 1.0, A);

						/*--- M = dV/dU for converting from conservative to primitive variables ---*/
						M[0][0] = 1.0;				M[0][1] = 0.0;		M[0][2] = 0.0;		M[0][3] = 0.0;
						M[1][0] = u;					M[1][1] = rho;		M[1][2] = 0.0;		M[1][3] = 0.0;
						M[2][0] = v;					M[2][1] = 0.0;		M[2][2] = rho;		M[2][3] = 0.0;
						M[3][0] = 0.5*sq_vel;	M[3][1] = rho*u;	M[3][2] = rho*v;	M[3][3] = 1.0/Gamma_Minus_One;

						/*--- Multiply [A] by [M] since we are taking the Jacobian of
             the r.h.s. w.r.t the primitive variables ---*/
						for (iVar = 0; iVar < nVar; iVar++)
							for (jVar = 0; jVar < nVar; jVar++) {
								aux = 0.0;
								for (kVar = 0; kVar < nVar; kVar++)
									aux += A[iVar][kVar]*M[kVar][jVar];
								AM[iVar][jVar] = aux;
							}

						/*--- Transpose the product of the A & M matrices ---*/
						for (iVar = 0; iVar < nVar; iVar++)
							for (jVar = 0; jVar < nVar; jVar++)
								A[iVar][jVar] = AM[jVar][iVar];

						/*--- Build right hand side: phi*dQ/dU ---*/
						/*--- Don't have the mesh velocity term yet - assume fixed ---*/

						Phi = wave_solution[MESH_0][WAVE_SOL]->node[iPoint]->GetSolution();
						U_i_old = flow_solution[MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n();

						/*--- Try to flip the normal here w.r.t the adjoint? not likely ---*/
						b[0] = 0.0;
						for (iDim = 0; iDim < nDim; iDim++) {
							b[0] += Phi[0]*(U_i[iDim+1]/U_i[0] - U_i_old[iDim+1]/U_i_old[0])*(UnitaryNormal[iDim]*Area)/delta_T;
							b[iDim+1] = Phi[0]*(U_i[0] - U_i_old[0])*(UnitaryNormal[iDim]*Area)/delta_T;
						}
						b[3] = 0.0;

					}

					if (nDim == 3) {
						// Do nothing in 3-D at the moment
					}

					/*--- Solve the system using a LU decomposition --*/
					for (jc = 1; jc < nVar; jc++)
						A[0][jc] /= A[0][0];

					jrjc = 0;
					for (;;) {
						jrjc++; jrjcm1 = jrjc-1; jrjcp1 = jrjc+1;
						for (jr = jrjc; jr < nVar; jr++) {
							sum = A[jr][jrjc];
							for (jm = 0; jm <= jrjcm1; jm++)
								sum -= A[jr][jm]*A[jm][jrjc];
							A[jr][jrjc] = sum;
						}
						if ( jrjc == (nVar-1) ) goto stop;
						for (jc = jrjcp1; jc<nVar; jc++) {
							sum = A[jrjc][jc];
							for (jm = 0; jm <= jrjcm1; jm++)
								sum -= A[jrjc][jm]*A[jm][jc];
							A[jrjc][jc] = sum/A[jrjc][jrjc];
						}
					}

					stop:

					b[0] = b[0]/A[0][0];
					for (jr = 1; jr<nVar; jr++) {
						jrm1 = jr-1;
						sum = b[jr];
						for (jm = 0; jm<=jrm1; jm++)
							sum -= A[jr][jm]*b[jm];
						b[jr] = sum/A[jr][jr];
					}

					for (jrjr = 1; jrjr<nVar; jrjr++) {
						jr = (nVar-1)-jrjr;
						jrp1 = jr+1;
						sum = b[jr];
						for (jmjm = jrp1; jmjm<nVar; jmjm++) {
							jm = (nVar-1)-jmjm+jrp1;
							sum -= A[jr][jm]*b[jm];
						}
						b[jr] = sum;
					}

					for (iVar = 0; iVar < nVar; iVar++)
						Psi_New[iVar] = b[iVar];

					/*--- Store new adjoint solution for use in BC_FWH() ---*/
					node[iPoint]->SetIntBoundary_Jump(Psi_New);

					/*--- Store old direct solution for computing derivative on the next step ---*/
					/*--- Get old solution for computing time derivative ---*/
					flow_solution[MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();

				}

			}

	/*--- Free memory and exit ---*/
	delete [] Psi_New;
	delete [] Velocity;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] A[iVar];
	delete [] A;


}

void CAdjEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep,
		unsigned short iMesh, unsigned short RunTime_EqSystem) {
	unsigned short iVar, jVar;
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1, Volume_nM1, Volume_n, Volume_nP1, TimeStep;

	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool FlowEq = (RunTime_EqSystem == RUNTIME_FLOW_SYS);
	bool AdjEq = (RunTime_EqSystem == RUNTIME_ADJFLOW_SYS);
	bool incompressible = config->GetIncompressible();
	bool Grid_Movement = config->GetGrid_Movement();

	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n   = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();

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
		for(iVar = 0; iVar < nVar; iVar++) {
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				Residual[iVar] = ( U_time_nP1[iVar]*Volume_nP1 - U_time_n[iVar]*Volume_n ) / TimeStep;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				Residual[iVar] = ( 3.0*U_time_nP1[iVar]*Volume_nP1 - 4.0*U_time_n[iVar]*Volume_n
						+  1.0*U_time_nM1[iVar]*Volume_nM1 ) / (2.0*TimeStep);
		}

		if ((incompressible && FlowEq) || (incompressible && AdjEq)) Residual[0] = 0.0;

		/*--- Add Residual ---*/
		LinSysRes.AddBlock(iPoint, Residual);

		if (implicit) {
			for (iVar = 0; iVar < nVar; iVar++) {
				for (jVar = 0; jVar < nVar; jVar++)
					Jacobian_i[iVar][jVar] = 0.0;

				if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
					Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
				if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
					Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
			}
			if ((incompressible && FlowEq) ||
					(incompressible && AdjEq)) Jacobian_i[0][0] = 0.0;
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}
  
}

CAdjNSSolver::CAdjNSSolver(void) : CAdjEulerSolver() { }

CAdjNSSolver::CAdjNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CAdjEulerSolver() {
	unsigned long iPoint, index, iVertex;
	string text_line, mesh_filename;
	unsigned short iDim, iVar, iMarker;
	ifstream restart_file;
	string filename, AdjExt;

	bool restart = config->GetRestart();
	bool incompressible = config->GetIncompressible();

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
  
	if (incompressible) nVar = nDim + 1;
	else nVar = nDim + 2;
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
				if (incompressible) {
					if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2];
					if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
				}
				else {
					if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
					if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
				}
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

CAdjNSSolver::~CAdjNSSolver(void) {

}


void CAdjNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
  
  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the flow problem may use different methods). ---*/
  bool implicit       = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool upwind_2nd     = (config->GetKind_Upwind_AdjFlow() == ROE_2ND) ||
  (config->GetKind_Upwind_AdjFlow() == SW_2ND);
  bool center         = (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst     = (config->GetKind_Centered_AdjFlow() == JST);
  bool limiter        = (config->GetKind_SlopeLimit() != NONE);
  bool incompressible = config->GetIncompressible();
  double adj_limit    = config->GetAdjointLimit();

  /*--- Compute nacelle inflow and exhaust properties ---*/
  GetNacelle_Properties(geometry, config, iMesh);
  
	/*--- Residual initialization ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Set the primitive variables incompressible and compressible
     adjoint variables ---*/
		if (incompressible) node[iPoint]->SetPrimVar_Incompressible(adj_limit);
		else node[iPoint]->SetPrimVar_Compressible(adj_limit);
    
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
	if ((implicit) || (config->GetKind_Adjoint() == DISCRETE) ) Jacobian.SetValZero();
  
}

void CAdjNSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
		CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iPoint, jPoint, iEdge;

	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
  
	if (true) {

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
			if (incompressible) {
        numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                              solver_container[FLOW_SOL]->node[jPoint]->GetDensityInc());
				numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
						solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosityInc());
      }
			else {
				numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
						solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
      }

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

}

void CAdjNSSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
		CConfig *config, unsigned short iMesh) {
	unsigned long iPoint;
  unsigned short iVar, iDim;
  
  bool rotating_frame = config->GetRotating_Frame();
    bool freesurface = config->GetFreeSurface();
  
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
			numerics->SetDistance(geometry->node[iPoint]->GetWallDistance(), 0.0);

		}

		/*--- Compute residual ---*/
		numerics->ComputeResidual(Residual, config);
    
    /*--- Add and substract to the residual ---*/
		LinSysRes.AddBlock(iPoint, Residual);
	}

  // WARNING: The rotating frame source term has been placed in the second
  // source term container since the section below is commented. This needs a
  // permanent fix asap!
  
  if (rotating_frame) {
    
		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
			/*--- Set solution  ---*/
			second_numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
      
			/*--- Set control volume ---*/
			second_numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
			/*--- Set rotational velocity ---*/
			second_numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
      
			/*--- Compute Residual ---*/
			second_numerics->ComputeResidual(Residual, Jacobian_i, config);
      
			/*--- Add Residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);
		}
	}
  
	if (freesurface) {
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      double Volume = geometry->node[iPoint]->GetVolume();
      double **Gradient = solver_container[ADJLEVELSET_SOL]->node[iPoint]->GetGradient();
      double coeff = solver_container[LEVELSET_SOL]->node[iPoint]->GetSolution(0) / solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
      
      Residual[0] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Residual[iDim+1] = coeff*Gradient[0][iDim]*Volume;
      }
      
      LinSysRes.AddBlock(iPoint, Residual);
      
		}
	}
  
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

void CAdjNSSolver::Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
	unsigned long iVertex, iPoint;
	unsigned short iDim, jDim, iMarker;
	double **PsiVar_Grad, **PrimVar_Grad, div_phi, *Normal, Area,
	normal_grad_psi5, normal_grad_T, sigma_partial,
  cp, Laminar_Viscosity, heat_flux_factor, LevelSet, Target_LevelSet, temp_sens;
  
  double *Psi, *U, Enthalpy, **GridVel_Grad, gradPsi5_v, psi5_tau_partial, psi5_tau_grad_vel, source_v;
  double Density, Pressure, div_vel, val_turb_ke, vartheta, vartheta_partial, psi5_p_div_vel, *Omega, rho_v[3], CrossProduct[3];
  double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  
	double Gas_Constant = config->GetGas_ConstantND();
	bool incompressible = config->GetIncompressible();
    bool rotating_frame = config->GetRotating_Frame();
    bool freesurface = config->GetFreeSurface();

  
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
  
	if (config->GetKind_Adjoint() != DISCRETE) {
    
		/*--- Compute gradient of adjoint variables on the surface ---*/
		SetSurface_Gradient(geometry, config);
    
    /*--- Compute gradient of the grid velocity if applicable - only need to do this once at beginning ---*/
    if (rotating_frame)
      SetRotVel_Gradient(geometry, config);
    
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
            
            if (incompressible) Laminar_Viscosity  = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc();
            else Laminar_Viscosity  = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
            
            heat_flux_factor = cp * Laminar_Viscosity / PRANDTL;
            
            /*--- Compute face area and the nondimensional normal to the surface ---*/
            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) { Area += Normal[iDim]*Normal[iDim]; } Area = sqrt(Area);
            for (iDim = 0; iDim < nDim; iDim++) { UnitaryNormal[iDim] = Normal[iDim] / Area; }
            
            /*--- Compute the sensitivity related to the temperature ---*/
            if (!incompressible) {
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
            } else {
              
              /*--- Incompressible case ---*/
              temp_sens = 0.0;
            }            
            
            /*--- Term: sigma_partial = \Sigma_{ji} n_i \partial_n v_j ---*/
            if (!incompressible) {
              div_phi = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) {
                div_phi += PsiVar_Grad[iDim+1][iDim];
                for (jDim = 0; jDim < nDim; jDim++)
                  Sigma[iDim][jDim] = Laminar_Viscosity * (PsiVar_Grad[iDim+1][jDim]+PsiVar_Grad[jDim+1][iDim]);
              }
              for (iDim = 0; iDim < nDim; iDim++)
                Sigma[iDim][iDim] -= TWO3*Laminar_Viscosity * div_phi;
            }
            else {
              for (iDim = 0; iDim < nDim; iDim++) {
                for (jDim = 0; jDim < nDim; jDim++)
                  Sigma[iDim][jDim] = Laminar_Viscosity * PsiVar_Grad[jDim+1][iDim];
              }
            }
            
            for (iDim = 0; iDim < nDim; iDim++) {
              normal_grad_vel[iDim] = 0.0;
              for (jDim = 0; jDim < nDim; jDim++)
                normal_grad_vel[iDim] += PrimVar_Grad[iDim+1][jDim]*UnitaryNormal[jDim];
            }
            
            sigma_partial = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              for (jDim = 0; jDim < nDim; jDim++)
                sigma_partial += UnitaryNormal[iDim]*Sigma[iDim][jDim]*normal_grad_vel[jDim];
            
            /*--- Compute additional terms in the surface sensitivity for
             moving walls in a rotating frame or dynamic mesh problem. ---*/
            if (rotating_frame) {
              
              Psi = node[iPoint]->GetSolution();
              U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
              Density = U[0];
              Pressure = solver_container[FLOW_SOL]->node[iPoint]->GetPressure(incompressible);
							Enthalpy = solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy();
              
              /*--- Turbulent kinetic energy ---*/
              if (config->GetKind_Turb_Model() == SST)
                val_turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
              else
                val_turb_ke = 0.0;
              
              div_vel = 0.0;
              for (iDim = 0 ; iDim < nDim; iDim++) {
                Velocity[iDim] = U[iDim+1]/Density;
                div_vel += PrimVar_Grad[iDim+1][iDim];
              }
              
              for (iDim = 0 ; iDim < nDim; iDim++)
                for (jDim = 0 ; jDim < nDim; jDim++)
                  tau[iDim][jDim] = Laminar_Viscosity*(PrimVar_Grad[jDim+1][iDim] + PrimVar_Grad[iDim+1][jDim])
                  - TWO3*Laminar_Viscosity*div_vel*delta[iDim][jDim]
                  - TWO3*Density*val_turb_ke*delta[iDim][jDim];
                            
              /*--- Form normal_grad_gridvel = \partial_n (u_x) ---*/
              GridVel_Grad = geometry->node[iPoint]->GetRotVel_Grad();
              for (iDim = 0; iDim < nDim; iDim++) {
                normal_grad_gridvel[iDim] = 0.0;
                for (jDim = 0; jDim < nDim; jDim++)
                  normal_grad_gridvel[iDim] += GridVel_Grad[iDim][jDim]*UnitaryNormal[jDim];
              }
              
              /*--- Form normal_grad_v_ux = \partial_n (v - u_x) ---*/
              for (iDim = 0; iDim < nDim; iDim++)
                normal_grad_v_ux[iDim] = normal_grad_vel[iDim] - normal_grad_gridvel[iDim];
              
              //if (iVertex == 10) cout << normal_grad_gridvel[0] << "  "<<normal_grad_vel[0]<<endl;
              /*--- Form Sigma_Psi5v ---*/
              gradPsi5_v = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) {
                gradPsi5_v += PsiVar_Grad[nDim+1][iDim]*Velocity[iDim];
                for (jDim = 0; jDim < nDim; jDim++)
                  Sigma_Psi5v[iDim][jDim] = Laminar_Viscosity * (PsiVar_Grad[nDim+1][iDim]*Velocity[jDim]+PsiVar_Grad[nDim+1][jDim]*Velocity[iDim]);
              }
              for (iDim = 0; iDim < nDim; iDim++)
                Sigma_Psi5v[iDim][iDim] -= TWO3*Laminar_Viscosity * gradPsi5_v;
             
              
              /*--- Now compute various terms of surface sensitivity ---*/
              
              /*--- Form vartheta_partial = \vartheta * \partial_n (v - u_x) . n ---*/
							vartheta = Density*Psi[0] + Density*Enthalpy*Psi[nDim+1];
              for (iDim = 0; iDim < nDim; iDim++) {
                vartheta += U[iDim+1]*Psi[iDim+1];
              }
              vartheta_partial = 0.0;
              for (iDim = 0; iDim < nDim; iDim++)
                vartheta_partial += vartheta * normal_grad_v_ux[iDim] * UnitaryNormal[iDim];
              
//              /*--- Form sigma_partial = n_i ( \Sigma_Phi_{ij} + \Sigma_Psi5v_{ij} ) \partial_n (v - u_x)_j ---*/
//              sigma_partial = 0.0;
//              for (iDim = 0; iDim < nDim; iDim++)
//                for (jDim = 0; jDim < nDim; jDim++)
//                  sigma_partial += UnitaryNormal[iDim]*(Sigma[iDim][jDim]+Sigma_Psi5v[iDim][jDim])*normal_grad_v_ux[jDim];
              
              /*--- Form psi5_tau_partial = \Psi_5 * \partial_n (v - u_x)_i * tau_{ij} * n_j ---*/
              psi5_tau_partial = 0.0;
              for (iDim = 0; iDim < nDim; iDim++)
                for (jDim = 0; jDim < nDim; jDim++)
                  psi5_tau_partial -= Psi[nDim+1]*normal_grad_v_ux[iDim]*tau[iDim][jDim]*UnitaryNormal[jDim];
              
              /*--- Form psi5_p_div_vel = ---*/
              psi5_p_div_vel = -Psi[nDim+1]*Pressure*div_vel;
              
              /*--- Form psi5_tau_grad_vel = \Psi_5 * tau_{ij} : \nabla( v ) ---*/
              psi5_tau_grad_vel = 0.0;
              for (iDim = 0; iDim < nDim; iDim++)
                for (jDim = 0; jDim < nDim; jDim++)
                  psi5_tau_grad_vel += Psi[nDim+1]*tau[iDim][jDim]*PrimVar_Grad[iDim+1][jDim];
                            
              /*--- Retrieve the angular velocity vector ---*/
              Omega = config->GetOmega_FreeStreamND();
              
              /*--- Calculate momentum source terms as: rho * ( Omega X V ) ---*/
              for(iDim = 0; iDim < nDim; iDim++)
                rho_v[iDim] = U[iDim+1];
              
              CrossProduct[0] = Omega[1]*rho_v[2] - Omega[2]*rho_v[1];
              CrossProduct[1] = Omega[2]*rho_v[0] - Omega[0]*rho_v[2];
              CrossProduct[2] = Omega[0]*rho_v[1] - Omega[1]*rho_v[0];
              
              source_v = 0.0;
              for(iDim = 0; iDim < nDim; iDim++)
                source_v -= Psi[nDim+1]*CrossProduct[iDim]*Velocity[iDim];
              
              /*--- For simplicity, store all additional terms within sigma_partial ---*/
              //if (iVertex==10)cout << sigma_partial << "\t" << vartheta_partial << "\t" << psi5_tau_partial << "\t" << psi5_p_div_vel << "\t" << psi5_tau_grad_vel << "\t" <<tang_psi_5<<endl;
              //sigma_partial = sigma_partial + vartheta_partial + psi5_tau_partial + psi5_p_div_vel + psi5_tau_grad_vel + source_v;
              //sigma_partial = sigma_partial + vartheta_partial + psi5_p_div_vel + psi5_tau_grad_vel;
              
            }
            
            /*--- Compute additional term in the surface sensitivity for
             free surface problem. ---*/
            if (freesurface) {
              LevelSet = solver_container[LEVELSET_SOL]->node[iPoint]->GetSolution(0);
              Target_LevelSet = geometry->node[iPoint]->GetCoord(nDim-1);
              sigma_partial += 0.5*(Target_LevelSet - LevelSet)*(Target_LevelSet - LevelSet);
            }
            
            /*--- Compute sensitivity for each surface point ---*/
            CSensitivity[iMarker][iVertex] = (sigma_partial - temp_sens)*Area;
            Sens_Geo[iMarker] -= CSensitivity[iMarker][iVertex]*Area;
          }
				}
				Total_Sens_Geo += Sens_Geo[iMarker];
			}
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

void CAdjNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
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
	bool incompressible = config->GetIncompressible();
    bool rotating_frame = config->GetRotating_Frame();
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
            if (rotating_frame) {
                GridVel = geometry->node[iPoint]->GetRotVel();
                for (iDim = 0; iDim < nDim; iDim++) {
                    phi[iDim] = d[iDim] - Psi[nVar-1]*GridVel[iDim];
                }
            } else if ( grid_movement) {
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
            if (incompressible) {
                
                /*--- Pressure residual due to the convective term ---*/
                l1psi = 0.0;
                for (iDim = 0; iDim < nDim; iDim++)
                    l1psi += Normal[iDim]*d[iDim];
                Res_Conv_i[0] = l1psi;
                
            } else {
                
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
                pressure = solver_container[FLOW_SOL]->node[iPoint]->GetPressure(incompressible);
                Res_Conv_i[0] = pressure*Sigma_5/(Gamma_Minus_One*rho*rho);
                Res_Conv_i[nVar-1] -= Sigma_5/rho;
                
                /*--- Flux contribution and Jacobian contributions for moving
                 walls. Note that these are only for the adjoint density and
                 adjoint energy equations (the adjoint vel. uses a strong BC). ---*/
                if (rotating_frame || grid_movement) {
                    
                    /*--- Get the appropriate grid velocity at this node ---*/
                    if (rotating_frame)
                        GridVel = geometry->node[iPoint]->GetRotVel();
                    else if (grid_movement)
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


void CAdjNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
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
	bool incompressible = config->GetIncompressible();
  bool rotating_frame = config->GetRotating_Frame();
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
      if (rotating_frame) {
        GridVel = geometry->node[iPoint]->GetRotVel();
        for (iDim = 0; iDim < nDim; iDim++) {
          phi[iDim] = d[iDim] - Psi[nVar-1]*GridVel[iDim];
        }
      } else if ( grid_movement) {
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
      if (incompressible) {
        
        /*--- Pressure residual due to the convective term ---*/
        l1psi = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          l1psi += Normal[iDim]*d[iDim];
        Res_Conv_i[0] = l1psi;
        
      } else {
        
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
        pressure = solver_container[FLOW_SOL]->node[iPoint]->GetPressure(incompressible);
        Res_Conv_i[0] = pressure*Sigma_5/(Gamma_Minus_One*rho*rho);
        
        /*--- Flux contribution and Jacobian contributions for moving
         walls. Note that these are only for the adjoint density and
         adjoint energy equations (the adjoint vel. uses a strong BC). ---*/
        if (rotating_frame || grid_movement) {
          
          /*--- Get the appropriate grid velocity at this node ---*/
          if (rotating_frame)
            GridVel = geometry->node[iPoint]->GetRotVel();
          else if (grid_movement)
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

