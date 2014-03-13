/*!
 * \file solution_adjoint_mean.cpp
 * \brief Main subrotuines for solving adjoint problems (Euler, Navier-Stokes, etc.).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.0.0 "eagle"
 *
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
  PsiRho_Inf   = NULL;
	Phi_Inf      = NULL;
	Sens_Mach    = NULL;
	Sens_AoA     = NULL;
	Sens_Geo     = NULL;
	Sens_Press   = NULL;
	Sens_Temp    = NULL;
	iPoint_UndLapl  = NULL;
	jPoint_UndLapl  = NULL;
	CSensitivity = NULL;
  Jacobian_Axisymmetric = NULL;
  
}

CAdjTNE2EulerSolver::CAdjTNE2EulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {

  bool restart;
	unsigned long iPoint, index, iVertex;
	unsigned short iDim, iSpecies, iVar, iMarker, nLineLets;
  double dull_val;
  string text_line, mesh_filename;
  string filename, AdjExt;
	ifstream restart_file;
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  /*--- Array initialization ---*/
  PsiRho_Inf   = NULL;
	Phi_Inf      = NULL;
	Sens_Mach    = NULL;
	Sens_AoA     = NULL;
	Sens_Geo     = NULL;
	Sens_Press   = NULL;
	Sens_Temp    = NULL;
	iPoint_UndLapl  = NULL;
	jPoint_UndLapl  = NULL;
	CSensitivity = NULL;
  
  /*--- Set booleans for solver settings ---*/
  restart      = config->GetRestart();
  
	/*--- Define constants in the solver structure ---*/
  nSpecies     = config->GetnSpecies();
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nDim         = geometry->GetnDim();
  
  /*--- Set the size of the solution array ---*/
	nVar         = nSpecies + nDim + 2;
  nPrimVar     = nSpecies + nDim + 8;
  nPrimVarGrad = nSpecies + nDim + 8;
  
  /*--- Allocate a CVariable array for each node of the mesh ---*/
	node = new CVariable*[nPoint];
  
	/*--- Define some auxiliary vectors related to the residual ---*/
	Residual     = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max    = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
	Residual_i   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
	Res_Conv_i   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_i[iVar]    = 0.0;
  Res_Visc_i   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_i[iVar]    = 0.0;
	Res_Conv_j   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_j[iVar]    = 0.0;
  Res_Visc_j   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_j[iVar]    = 0.0;
  
	/*--- Define some auxiliary vectors related to the solution ---*/
	Solution   = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
	Solution_i = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar]   = 0.0;
  Solution_j = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar]   = 0.0;
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
	Vector   = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
	Vector_i = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
	Vector_j = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
	/*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
	if (config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_CENTERED) {
		iPoint_UndLapl = new double [nPoint];
		jPoint_UndLapl = new double [nPoint];
	}
  
	/*--- Define some auxiliary vectors related to the geometry ---*/
	Vector_i = new double[nDim];
  Vector_j = new double[nDim];
  
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
  
	/*--- Allocate Jacobians for implicit time-stepping ---*/
  if (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT) {
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
      cout << "Initialize Jacobian structure (Adjoint Euler). MG level: " << iMesh <<"." << endl;
		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry);
    
    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
  
  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No jacobian structure (Adjoint Euler). MG level: " << iMesh <<"." << endl;
  }
  
	/*--- Allocate arrays for gradient computation by least squares ---*/
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
  
	/*--- Sensitivity arrays on boundaries ---*/
	CSensitivity = new double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		CSensitivity[iMarker] = new double [geometry->nVertex[iMarker]];
	}
	Sens_Geo   = new double[nMarker];
	Sens_Mach  = new double[nMarker];
	Sens_AoA   = new double[nMarker];
	Sens_Press = new double[nMarker];
	Sens_Temp  = new double[nMarker];
  
  /*--- Initialize sensitivity arrays ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		Sens_Geo[iMarker]   = 0.0;
		Sens_Mach[iMarker]  = 0.0;
		Sens_AoA[iMarker]   = 0.0;
		Sens_Press[iMarker] = 0.0;
		Sens_Temp[iMarker]  = 0.0;
		for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
			CSensitivity[iMarker][iVertex] = 0.0;
	}
  
	/*--- Adjoint flow at the inifinity, initialization stuff ---*/
  PsiRho_Inf = new double [nSpecies];
  Phi_Inf    = new double [nDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    PsiRho_Inf[iSpecies] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Phi_Inf[iDim] = 0.0;
  PsiE_Inf   = 0.0;
  PsiEve_Inf = 0.0;
  
	if (!restart || geometry->GetFinestMGLevel() == false) {
		/*--- Restart the solution from infinity ---*/
		for (iPoint = 0; iPoint < nPoint; iPoint++)
			node[iPoint] = new CAdjTNE2EulerVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, PsiEve_Inf, nDim, nVar, config);
	}
	else {
    
		/*--- Restart the solution from file information ---*/
		mesh_filename = config->GetSolution_AdjFileName();
    
		/*--- Change the name, depending of the objective function ---*/
		filename.assign(mesh_filename);
		filename.erase (filename.end()-4, filename.end());
		switch (config->GetKind_ObjFunc()) {
      case DRAG_COEFFICIENT:      AdjExt = "_cd.dat"; break;
      case LIFT_COEFFICIENT:      AdjExt = "_cl.dat"; break;
      case SIDEFORCE_COEFFICIENT: AdjExt = "_csf.dat"; break;
      case PRESSURE_COEFFICIENT:  AdjExt = "_cp.dat"; break;
      case MOMENT_X_COEFFICIENT:  AdjExt = "_cmx.dat"; break;
      case MOMENT_Y_COEFFICIENT:  AdjExt = "_cmy.dat"; break;
      case MOMENT_Z_COEFFICIENT:  AdjExt = "_cmz.dat"; break;
      case EFFICIENCY:            AdjExt = "_eff.dat"; break;
      case EQUIVALENT_AREA:       AdjExt = "_ea.dat"; break;
      case NEARFIELD_PRESSURE:    AdjExt = "_nfp.dat"; break;
      case FORCE_X_COEFFICIENT:   AdjExt = "_cfx.dat"; break;
      case FORCE_Y_COEFFICIENT:   AdjExt = "_cfy.dat"; break;
      case FORCE_Z_COEFFICIENT:   AdjExt = "_cfz.dat"; break;
      case THRUST_COEFFICIENT:    AdjExt = "_ct.dat"; break;
      case TORQUE_COEFFICIENT:    AdjExt = "_cq.dat"; break;
      case FIGURE_OF_MERIT:       AdjExt = "_merit.dat"; break;
      case FREE_SURFACE:          AdjExt = "_fs.dat"; break;
      case HEAT_LOAD:             AdjExt = "_Q.dat"; break;
		}
		filename.append(AdjExt);
		restart_file.open(filename.data(), ios::in);
    
		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
			cout << "There is no adjoint restart file!! " << filename.data() << "."<< endl;
			exit(1);
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
        
        point_line >> index;
        for (iDim = 0; iDim < nDim; iDim++)
          point_line >> dull_val;
        for (iVar = 0; iVar < nVar; iVar++)
          point_line >> Solution[iVar];
        
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
	if (config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_CENTERED) space_centered = true;
	else space_centered = false;
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
}

CAdjTNE2EulerSolver::~CAdjTNE2EulerSolver(void) {
  unsigned short iMarker;
  
  if (PsiRho_Inf  != NULL) delete [] PsiRho_Inf;
  if (Phi_Inf     != NULL) delete [] Phi_Inf;
	if (Sens_Mach   != NULL) delete [] Sens_Mach;
	if (Sens_AoA    != NULL) delete [] Sens_AoA;
	if (Sens_Geo    != NULL) delete [] Sens_Geo;
	if (Sens_Press  != NULL) delete [] Sens_Press;
	if (Sens_Temp   != NULL) delete [] Sens_Temp;
	if (iPoint_UndLapl != NULL) delete [] iPoint_UndLapl;
	if (jPoint_UndLapl != NULL) delete [] jPoint_UndLapl;
  
	if (CSensitivity != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete CSensitivity[iMarker];
    delete [] CSensitivity;
  }
  
}

void CAdjTNE2EulerSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR;
  unsigned long nBufferS_Vector, nBufferR_Vector;
  int send_to, receive_from;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta;
  double phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
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
#ifdef WINDOWS
	  MPI_Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                               Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, NULL);
#else
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
#endif
      
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
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        
        rotMatrix[0][1] = cosPhi*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        
        rotMatrix[0][2] = -sinPhi;
        rotMatrix[1][2] = sinTheta*cosPhi;
        rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex];
        }
        else {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+2] = rotMatrix[2][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
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
	unsigned long iVertex, iPoint, nVertexS, nVertexR;
  unsigned long nBufferS_Vector, nBufferR_Vector;
  int send_to, receive_from;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta;
  double phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
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
#ifdef WINDOWS
	  MPI_Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                               Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, NULL);
#else
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
#endif
      
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
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        
        rotMatrix[0][1] = cosPhi*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        
        rotMatrix[0][2] = -sinPhi;
        rotMatrix[1][2] = sinTheta*cosPhi;
        rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex];
        }
        else {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+2] = rotMatrix[2][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
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
	unsigned long iVertex, iPoint, nVertexS, nVertexR;
  unsigned long nBufferS_Vector, nBufferR_Vector;
  int send_to, receive_from;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta;
  double phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  double *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
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
      
      /*--- Copy the limiter that should be sent ---*/
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
#ifdef WINDOWS
	  MPI_Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                               Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, NULL);
#else
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Limit, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
#endif
      
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
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        
        rotMatrix[0][1] = cosPhi*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        
        rotMatrix[0][2] = -sinPhi;
        rotMatrix[1][2] = sinTheta*cosPhi;
        rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex];
        }
        else {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_Limit[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_Limit[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+2] = rotMatrix[2][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_Limit[(nSpecies+2)*nVertexR+iVertex];
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
	unsigned long iVertex, iPoint, nVertexS, nVertexR;
  unsigned long nBufferS_Vector, nBufferR_Vector;
  int send_to, receive_from;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta;
  double phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  double *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
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
			nBufferS_Vector = nVertexS*nVar*nDim;   nBufferR_Vector = nVertexR*nVar*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new double [nBufferR_Vector];
      Buffer_Send_Gradient = new double[nBufferS_Vector];
      
      /*--- Copy the gradient that should be sent ---*/
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
#ifdef WINDOWS
	  MPI_Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                               Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, NULL);
#else
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Gradient, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
#endif
      
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
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        
        rotMatrix[0][1] = cosPhi*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        
        rotMatrix[0][2] = -sinPhi;
        rotMatrix[1][2] = sinTheta*cosPhi;
        rotMatrix[2][2] = cosTheta*cosPhi;
        
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
	unsigned long iVertex, iPoint, nVertexS, nVertexR;
  unsigned long nBufferS_Vector, nBufferR_Vector;
  int send_to, receive_from;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta;
  double phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  double *Buffer_Receive_Undivided_Laplacian = NULL;
  double *Buffer_Send_Undivided_Laplacian = NULL;
  
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
#ifdef WINDOWS
	  MPI_Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                               Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, NULL);
#else
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
#endif
      
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
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        
        rotMatrix[0][1] = cosPhi*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        
        rotMatrix[0][2] = -sinPhi;
        rotMatrix[1][2] = sinTheta*cosPhi;
        rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex];
        }
        else {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_Undivided_Laplacian[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_Undivided_Laplacian[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+2] = rotMatrix[2][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_Undivided_Laplacian[(nSpecies+2)*nVertexR+iVertex];
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

void CAdjTNE2EulerSolver::SetForceProj_Vector(CGeometry *geometry,
                                              CSolver **solver_container,
                                              CConfig *config) {
  bool ionization;
	unsigned short iMarker, iDim, nHeavy, nEl;
	unsigned long iVertex, iPoint;
	double Alpha      = (config->GetAoA()*PI_NUMBER)/180.0;
	double Beta       = (config->GetAoS()*PI_NUMBER)/180.0;
	double RefAreaCoeff    = config->GetRefAreaCoeff();
	double RefLengthMoment  = config->GetRefLengthMoment();
	double *RefOriginMoment = config->GetRefOriginMoment(0);
  double *ForceProj_Vector, x = 0.0, y = 0.0, z = 0.0, *Normal, C_d, C_l, C_t, C_q;
	double x_origin, y_origin, z_origin, WDrag, Area;
	double RefVel2, RefDensity;
  int rank = MASTER_NODE;
  
#ifndef NO_MPI
#ifdef WINDOWS
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
  rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
	ForceProj_Vector = new double[nDim];
  
  /*--- Determine the number of heavy species ---*/
  ionization = config->GetIonization();
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
  
	/*--- Acquire free stream velocity & density ---*/
  RefVel2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    RefVel2  += solver_container[TNE2_SOL]->GetVelocity_Inf(iDim)
               *solver_container[TNE2_SOL]->GetVelocity_Inf(iDim);
  RefDensity = solver_container[TNE2_SOL]->node_infty->GetDensity();
  
	/*--- In parallel computations the Cd, and Cl must be recomputed using all the processors ---*/
#ifdef NO_MPI
	C_d = solver_container[TNE2_SOL]->GetTotal_CDrag();
	C_l = solver_container[TNE2_SOL]->GetTotal_CLift();
	C_t = solver_container[TNE2_SOL]->GetTotal_CT();
	C_q = solver_container[TNE2_SOL]->GetTotal_CQ();
#else
	double *sbuf_force = new double[4];
	double *rbuf_force = new double[4];
	sbuf_force[0] = solver_container[TNE2_SOL]->GetTotal_CDrag();
	sbuf_force[1] = solver_container[TNE2_SOL]->GetTotal_CLift();
	sbuf_force[2] = solver_container[TNE2_SOL]->GetTotal_CT();
	sbuf_force[3] = solver_container[TNE2_SOL]->GetTotal_CQ();
#ifdef WINDOWS
	MPI_Reduce(sbuf_force, rbuf_force, 4, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
	MPI_Bcast(rbuf_force, 4, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
	MPI::COMM_WORLD.Reduce(sbuf_force, rbuf_force, 4, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
	MPI::COMM_WORLD.Bcast(rbuf_force, 4, MPI::DOUBLE, MASTER_NODE);
#endif
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
  
	x_origin = RefOriginMoment[0];
  y_origin = RefOriginMoment[1];
  z_origin = RefOriginMoment[2];
  
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
            if ((nDim == 2) && (rank == MASTER_NODE)) { cout << "This functional is not possible in 2D!!" << endl;
              exit(1);
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
            if ((nDim == 2) && (rank == MASTER_NODE)) { cout << "This functional is not possible in 2D!!" << endl; exit(1); }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = -C_p*(z - z_origin)/RefLengthMoment; ForceProj_Vector[2] = C_p*(y - y_origin)/RefLengthMoment; }
            break;
          case MOMENT_Y_COEFFICIENT :
            if ((nDim == 2) && (rank == MASTER_NODE)) { cout << "This functional is not possible in 2D!!" << endl; exit(1); }
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
            if ((nDim == 2) && (rank == MASTER_NODE)) {cout << "This functional is not possible in 2D!!" << endl;
              exit(1);
            }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = C_p; }
            break;
          case THRUST_COEFFICIENT :
            if ((nDim == 2) && (rank == MASTER_NODE)) {cout << "This functional is not possible in 2D!!" << endl;
              exit(1);
            }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = C_p; }
            break;
          case TORQUE_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = C_p*(y - y_origin)/RefLengthMoment; ForceProj_Vector[1] = -C_p*(x - x_origin)/RefLengthMoment; }
            if (nDim == 3) { ForceProj_Vector[0] = C_p*(y - y_origin)/RefLengthMoment; ForceProj_Vector[1] = -C_p*(x - x_origin)/RefLengthMoment; ForceProj_Vector[2] = 0; }
            break;
          case FIGURE_OF_MERIT :
            if ((nDim == 2) && (rank == MASTER_NODE)) {cout << "This functional is not possible in 2D!!" << endl;
              exit(1);
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


void CAdjTNE2EulerSolver::SetInitialCondition(CGeometry **geometry,
                                              CSolver ***solver_container,
                                              CConfig *config,
                                              unsigned long ExtIter) {

	unsigned long iPoint, Point_Fine;
	unsigned short iMesh, iChildren, iVar;
	double Area_Children, Area_Parent;
  double *Solution, *Solution_Fine;
  
	bool restart = config->GetRestart();
  
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
          Solution_Fine = solver_container[iMesh-1][ADJTNE2_SOL]->node[Point_Fine]->GetSolution();
          for (iVar = 0; iVar < nVar; iVar++) {
            Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
          }
        }
        solver_container[iMesh][ADJTNE2_SOL]->node[iPoint]->SetSolution(Solution);
        
      }
      solver_container[iMesh][ADJTNE2_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    }
    delete [] Solution;
  }
}

void CAdjTNE2EulerSolver::Preprocessing(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CConfig *config,
                                        unsigned short iMesh,
                                        unsigned short iRKStep,
                                        unsigned short RunTime_EqSystem, bool Output) {

  bool implicit, upwind_2nd, center, center_jst, limiter, RightSol;
	unsigned long iPoint, ErrorCounter = 0;
  double SharpEdge_Distance;
  int rank;
  
#ifdef NO_MPI
	rank = MASTER_NODE;
#else
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the direct problem may use different methods). ---*/
  upwind_2nd = ((config->GetKind_Upwind_AdjTNE2() == ROE_2ND)  ||
                (config->GetKind_Upwind_AdjTNE2() == SW_2ND)   ||
                (config->GetKind_Upwind_AdjTNE2() == MSW_2ND)  ||
                (config->GetKind_Upwind_AdjTNE2() == AUSM_2ND) ||
                (config->GetKind_Upwind_AdjTNE2() == AUSMPWPLUS_2ND));
  center     = (config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_CENTERED );
  center_jst = (config->GetKind_Centered_AdjTNE2()      == JST            );
  implicit   = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT );
  limiter    = (config->GetKind_SlopeLimit()            != NONE           );
  
	/*--- Residual initialization ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Get the distance form a sharp edge ---*/
    SharpEdge_Distance = geometry->node[iPoint]->GetSharpEdge_Distance();
    
    /*--- Set the primitive variables incompressible and compressible
     adjoint variables ---*/
    RightSol = node[iPoint]->SetPrimVar_Compressible(SharpEdge_Distance,
                                                     false, config);
    if (!RightSol) ErrorCounter++;
    
		/*--- Initialize the convective residual vector ---*/
		LinSysRes.SetBlock_Zero(iPoint);
	}
  
  /*--- Upwind second order reconstruction ---*/
  if ((upwind_2nd) && (iMesh == MESH_0)) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
      SetSolution_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
      SetSolution_Gradient_LS(geometry, config);
    
    /*--- Limiter computation ---*/
		if (limiter) SetSolution_Limiter(geometry, config);
	}
  
  /*--- Artificial dissipation ---*/
  if (center) {
    if ((center_jst) && (iMesh == MESH_0)) {
      SetDissipation_Switch(geometry, config);
      SetUndivided_Laplacian(geometry, config);
      if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
        SetSolution_Gradient_GG(geometry, config);
      if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
        SetSolution_Gradient_LS(geometry, config);
    }
  }
  
	/*--- Implicit solution ---*/
	if (implicit) Jacobian.SetValZero();
  
  /*--- Error message ---*/
#ifndef NO_MPI
  unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
#ifdef WINDOWS
  MPI_Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI::UNSIGNED_LONG, MPI::SUM);
#endif
#endif
  if ((ErrorCounter != 0) && (rank == MASTER_NODE) && (iMesh == MESH_0))
    cout <<"The solution contains "<< ErrorCounter << " non-physical points." << endl;
  
}

void CAdjTNE2EulerSolver::Centered_Residual(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CNumerics *numerics,
                                            CConfig *config,
                                            unsigned short iMesh,
                                            unsigned short iRKStep) {
  bool implicit, high_order_diss;
  unsigned short iVar, jVar;
	unsigned long iEdge, iPoint, jPoint;
  
  /*--- Set booleans from configuration settings ---*/
	implicit        = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
	high_order_diss = ((config->GetKind_Centered_AdjTNE2() == JST)
                     && (iMesh == MESH_0));
  
  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( solver_container[TNE2_SOL]->node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( solver_container[TNE2_SOL]->node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( solver_container[TNE2_SOL]->node[0]->GetPIndex()       );
  numerics->SetTIndex      ( solver_container[TNE2_SOL]->node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( solver_container[TNE2_SOL]->node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( solver_container[TNE2_SOL]->node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( solver_container[TNE2_SOL]->node[0]->GetHIndex()       );
  numerics->SetAIndex      ( solver_container[TNE2_SOL]->node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvveIndex() );
  
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Identify points on the edge, normal, and neighbors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
		numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(),
                          geometry->node[jPoint]->GetnNeighbor());
    
		/*--- Pass the adjoint variables w/o reconstruction to CNumerics ---*/
		numerics->SetAdjointVar(node[iPoint]->GetSolution(),
                            node[jPoint]->GetSolution());
    
		/*--- Pass conservative & primitive variables w/o reconstruction ---*/
		numerics->SetConservative(solver_container[TNE2_SOL]->node[iPoint]->GetSolution(),
                              solver_container[TNE2_SOL]->node[jPoint]->GetSolution());
    numerics->SetPrimitive(solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar(),
                           solver_container[TNE2_SOL]->node[jPoint]->GetPrimVar());    

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU(  solver_container[TNE2_SOL]->node[iPoint]->GetdPdU(),
                        solver_container[TNE2_SOL]->node[jPoint]->GetdPdU());
    numerics->SetdTdU(  solver_container[TNE2_SOL]->node[iPoint]->GetdTdU(),
                        solver_container[TNE2_SOL]->node[jPoint]->GetdTdU());
    numerics->SetdTvedU(solver_container[TNE2_SOL]->node[iPoint]->GetdTvedU(),
                        solver_container[TNE2_SOL]->node[jPoint]->GetdTvedU());
    
    /*--- Set the value of the largest eigenvalue ---*/
    numerics->SetLambda(solver_container[TNE2_SOL]->node[iPoint]->GetLambda(),
                        solver_container[TNE2_SOL]->node[jPoint]->GetLambda());
    
		if (high_order_diss) {
			numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(),
                                       node[jPoint]->GetUndivided_Laplacian());
			numerics->SetSensor(solver_container[TNE2_SOL]->node[iPoint]->GetSensor(),
                          solver_container[TNE2_SOL]->node[jPoint]->GetSensor());
		}  
    
		/*--- Compute residuals ---*/
		numerics->ComputeResidual(Res_Conv_i, Res_Visc_i, Res_Conv_j, Res_Visc_j,
                              Jacobian_ii, Jacobian_ij, Jacobian_ji,
                              Jacobian_jj, config);
    
    /*--- Error checking ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      if ((Res_Conv_i[iVar] != Res_Conv_i[iVar]) ||
          (Res_Visc_i[iVar] != Res_Visc_i[iVar])) {
        cout << "NaN in Centered Residual" << endl;
      }
      for (jVar = 0; jVar < nVar; jVar++) {
        if (Jacobian_ii[iVar][jVar] != Jacobian_ii[iVar][jVar])
          cout << "NaN in Centered Jacobian i" << endl;
        if (Jacobian_jj[iVar][jVar] != Jacobian_jj[iVar][jVar])
          cout << "NaN in Centered Jacobian j" << endl;
      }
    }
    
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


void CAdjTNE2EulerSolver::Upwind_Residual(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CNumerics *numerics, CConfig *config,
                                          unsigned short iMesh) {
	
  bool implicit, high_order_diss, limiter;
	unsigned short iDim, iVar;
    unsigned long iEdge, iPoint, jPoint;
  double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
  double *Limiter_i, *Limiter_j;
  double *Psi_i, *Psi_j;
  double *U_i, *U_j, *V_i, *V_j;
  
  /*--- Initialize ---*/
  Limiter_i = NULL;
  Limiter_j = NULL;
  Psi_i     = NULL;
  Psi_j     = NULL;
  
  /*--- Set booleans based on config settings ---*/
	implicit        = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  limiter         = (config->GetKind_SlopeLimit_AdjTNE2() != NONE);
	high_order_diss = (((config->GetKind_Upwind_AdjTNE2() == ROE_2ND)  ||
                      (config->GetKind_Upwind_AdjTNE2() == SW_2ND)   ||
                      (config->GetKind_Upwind_AdjTNE2() == AUSM_2ND) ||
                      (config->GetKind_Upwind_AdjTNE2() == MSW_2ND)    )
                     && (iMesh == MESH_0));
  
  if (high_order_diss)
    cout << "WARNING!!! Upwind_Residual: 2nd order accuracy not in place!" << endl;
  if (limiter)
    cout << "WARNING!!! Upwind_Residual: Limiter not in place!" << endl;
  
  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( solver_container[TNE2_SOL]->node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( solver_container[TNE2_SOL]->node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( solver_container[TNE2_SOL]->node[0]->GetPIndex()       );
  numerics->SetTIndex      ( solver_container[TNE2_SOL]->node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( solver_container[TNE2_SOL]->node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( solver_container[TNE2_SOL]->node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( solver_container[TNE2_SOL]->node[0]->GetHIndex()       );
  numerics->SetAIndex      ( solver_container[TNE2_SOL]->node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvveIndex() );
  
  /*--- Loop over edges and calculate convective fluxes ---*/
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Retrieve node numbers and pass edge normal to CNumerics ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Pass conserved and primitive variables from CVariable to CNumerics class ---*/
    U_i = solver_container[TNE2_SOL]->node[iPoint]->GetSolution();
    U_j = solver_container[TNE2_SOL]->node[jPoint]->GetSolution();
    V_i = solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar();
    V_j = solver_container[TNE2_SOL]->node[jPoint]->GetPrimVar();
    numerics->SetPrimitive(V_i, V_j);
    numerics->SetConservative(U_i, U_j);
    
    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU(solver_container[TNE2_SOL]->node[iPoint]->GetdPdU(),
                      solver_container[TNE2_SOL]->node[jPoint]->GetdPdU());
    numerics->SetdTdU(solver_container[TNE2_SOL]->node[iPoint]->GetdTdU(),
                      solver_container[TNE2_SOL]->node[jPoint]->GetdTdU());
    numerics->SetdTvedU(solver_container[TNE2_SOL]->node[iPoint]->GetdTvedU(),
                        solver_container[TNE2_SOL]->node[jPoint]->GetdTvedU());
    
    /*--- Adjoint variables w/o reconstruction ---*/
    Psi_i = solver_container[ADJTNE2_SOL]->node[iPoint]->GetSolution();
    Psi_j = solver_container[ADJTNE2_SOL]->node[jPoint]->GetSolution();
    numerics->SetAdjointVar(Psi_i, Psi_j);    
    
		/*--- High order reconstruction using MUSCL strategy ---*/
    if (high_order_diss){
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(  geometry->node[jPoint]->GetCoord(iDim)
                              - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(  geometry->node[iPoint]->GetCoord(iDim)
                              - geometry->node[jPoint]->GetCoord(iDim));
      }
      Gradient_i = node[iPoint]->GetGradient();
      Gradient_j = node[jPoint]->GetGradient();
      if (limiter) {
        Limiter_i = node[iPoint]->GetLimiter();
        Limiter_j = node[jPoint]->GetLimiter();
      }
      
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
    numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
                              Jacobian_ji, Jacobian_jj, config);
    
    /*--- Error checking ---*/
//    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
//      if (Residual_i[iVar] != Residual_i[iVar]) {
//        cout << "NaN in Convective Residual" << endl;
//      }
//      for (unsigned short jVar = 0; jVar < nVar; jVar++) {
//        if (Jacobian_ii[iVar][jVar] != Jacobian_ii[iVar][jVar])
//          cout << "NaN in Convective Jacobian i" << endl;
//        if (Jacobian_jj[iVar][jVar] != Jacobian_jj[iVar][jVar])
//          cout << "NaN in Convective Jacobian j" << endl;
//      }
//    }
    
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

void CAdjTNE2EulerSolver::Source_Residual(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CNumerics *numerics,
                                          CNumerics *second_numerics,
                                          CConfig *config,
                                          unsigned short iMesh) {
  bool implicit;
	unsigned short iVar, jVar;
	unsigned long iPoint;
  
  /*--- Assign booleans from config settings ---*/
  implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( solver_container[TNE2_SOL]->node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( solver_container[TNE2_SOL]->node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( solver_container[TNE2_SOL]->node[0]->GetPIndex()       );
  numerics->SetTIndex      ( solver_container[TNE2_SOL]->node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( solver_container[TNE2_SOL]->node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( solver_container[TNE2_SOL]->node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( solver_container[TNE2_SOL]->node[0]->GetHIndex()       );
  numerics->SetAIndex      ( solver_container[TNE2_SOL]->node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvveIndex() );  
  
  /*--- Loop over points in the domain ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Set conserved & primitive variables at point i ---*/
    numerics->SetConservative(solver_container[TNE2_SOL]->node[iPoint]->GetSolution(),
                              solver_container[TNE2_SOL]->node[iPoint]->GetSolution());
    numerics->SetPrimitive(solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar(),
                           solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar());
    
    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU(solver_container[TNE2_SOL]->node[iPoint]->GetdPdU(),
                      solver_container[TNE2_SOL]->node[iPoint]->GetdPdU());
    numerics->SetdTdU(solver_container[TNE2_SOL]->node[iPoint]->GetdTdU(),
                      solver_container[TNE2_SOL]->node[iPoint]->GetdTdU());
    numerics->SetdTvedU(solver_container[TNE2_SOL]->node[iPoint]->GetdTvedU(),
                        solver_container[TNE2_SOL]->node[iPoint]->GetdTvedU());
    
    /*--- Set adjoint variables at point i ---*/
    numerics->SetAdjointVar(node[iPoint]->GetSolution(),
                            node[iPoint]->GetSolution());
    
    /*--- Set volume of the dual grid cell ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    
    /*--- Compute chemistry source terms ---*/
    numerics->ComputeChemistry(Residual_i, Jacobian_i, config);
    
    /*--- Error checking ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      if (Residual_i[iVar] != Residual_i[iVar]) {
        cout << "NaN in Chemistry Residual" << endl;
      }
      for (jVar = 0; jVar < nVar; jVar++) {
        if (Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar])
          cout << "NaN in Chemistry Jacobian i" << endl;
      }
    }
    
    /*--- Take the transpose of the source Jacobian matrix ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        Jacobian_ii[iVar][jVar] = Jacobian_i[jVar][iVar];
    
    /*--- Compute the adjoint source term residual (dQ/dU^T * Psi) ---*/
    for (iVar = 0; iVar < nVar; iVar ++) {
      Residual[iVar] = 0.0;
      for (jVar = 0; jVar < nVar; jVar++) {
        Residual[iVar] += Jacobian_ii[iVar][jVar] * node[iPoint]->GetSolution(jVar);
      }
    }
    
    /*--- Subtract Residual (and Jacobian) ---*/
    LinSysRes.SubtractBlock(iPoint, Residual);
    if (implicit)
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
    
    /*--- Compute vibrational relaxation source terms ---*/
    numerics->ComputeVibRelaxation(Residual_i, Jacobian_i, config);
    
    /*--- Error checking ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      if (Residual_i[iVar] != Residual_i[iVar]) {
        cout << "NaN in Energy Exchange Residual" << endl;
      }
      for (jVar = 0; jVar < nVar; jVar++) {
        if (Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar])
          cout << "NaN in Energy Exchange Jacobian i" << endl;
      }
    }
    
    /*--- Take the transpose of the source Jacobian matrix ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        Jacobian_ii[iVar][jVar] = Jacobian_i[jVar][iVar];
    
    /*--- Compute the adjoint source term residual (dQ/dU^T * Psi) ---*/
    for (iVar = 0; iVar < nVar; iVar ++) {
      Residual[iVar] = 0.0;
      for (jVar = 0; jVar < nVar; jVar++) {
        Residual[iVar] = Jacobian_ii[iVar][jVar] * node[iPoint]->GetSolution(jVar);
      }
    }
    
    /*--- Subtract Residual (and Jacobian) ---*/
    LinSysRes.SubtractBlock(iPoint, Residual);
    if (implicit)
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);

  }
}

void CAdjTNE2EulerSolver::SetUndivided_Laplacian(CGeometry *geometry,
                                                 CConfig *config) {
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

void CAdjTNE2EulerSolver::ExplicitEuler_Iteration(CGeometry *geometry,
                                                  CSolver **solver_container,
                                                  CConfig *config) {
	double *local_Residual, Vol, Delta, Res;
	unsigned short iVar;
	unsigned long iPoint;
  
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Delta = solver_container[TNE2_SOL]->node[iPoint]->GetDelta_Time() / Vol;
    
//		local_Res_TruncError = node[iPoint]->GetResTruncError();
		local_Residual = LinSysRes.GetBlock(iPoint);
    
		for (iVar = 0; iVar < nVar; iVar++) {
      Res = local_Residual[iVar];// + local_Res_TruncError[iVar];
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


void CAdjTNE2EulerSolver::ImplicitEuler_Iteration(CGeometry *geometry,
                                                  CSolver **solver_container,
                                                  CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index, IterLinSol=0;
	double Delta, Vol;
  
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
		Delta = Vol / solver_container[TNE2_SOL]->node[iPoint]->GetDelta_Time();
		Jacobian.AddVal2Diag(iPoint, Delta);
    
    if (Delta <= 0 || Delta != Delta) {
      cout << "NaN in Timestep" << endl;
    }
    
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			LinSysRes[total_index] = -(LinSysRes[total_index]);
			LinSysSol[total_index] = 0.0;
			AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex());
      
      if (LinSysRes[total_index] != LinSysRes[total_index])
        cout << "Linsysres NaN!" << endl;
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
    precond = new CLineletPreconditioner(Jacobian, geometry, config);
  }
  
  CSysSolve system;
  if (config->GetKind_Linear_Solver() == BCGSTAB)
    IterLinSol = system.BCGSTAB(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                                config->GetLinear_Solver_Iter(), false);
  else if (config->GetKind_Linear_Solver() == FGMRES)
    IterLinSol = system.FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                               config->GetLinear_Solver_Iter(), false);
  
  /*--- The the number of iterations of the linear solver ---*/
  SetIterLinSolver(IterLinSol);
  
  /*--- Deallocate memory ---*/
  delete mat_vec;
  delete precond;
  
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar]);
    }
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CAdjTNE2EulerSolver::Inviscid_Sensitivity(CGeometry *geometry,
                                               CSolver **solver_container,
                                               CNumerics *numerics,
                                               CConfig *config) {
	unsigned long iSpecies, iVar, jVar, iVertex, iPoint, Neigh;
	unsigned short iDim, iMarker, iNeigh;
  unsigned short P_INDEX, VEL_INDEX, RHO_INDEX, H_INDEX;
	double *d, *Normal, *UnitNormal;
  double *Psi, *U, *V, *dPdU, *USens;
  double rho, rhou, rhov, rhow, rhoE, rhoEve, H, p;
  double conspsi, Area;
  double Mach_Inf;
  double **PrimVar_Grad, **ConsVar_Grad, *ConsPsi_Grad;
  double ConsPsi, d_press, grad_v, v_gradconspsi;
  
  /*--- Initialization ---*/
  d            = NULL;
  Normal       = NULL;
  Psi          = NULL;
  U            = NULL;
  USens        = NULL;
  PrimVar_Grad = NULL;
  ConsVar_Grad = NULL;
  ConsPsi_Grad = NULL;
  
  /*--- Allocate arrays ---*/
  UnitNormal = new double[nDim];
  USens      = new double[nVar];
  
	/*--- Initialize sensitivities to zero ---*/
	Total_Sens_Geo   = 0.0;
  Total_Sens_Mach  = 0.0;
  Total_Sens_AoA   = 0.0;
	Total_Sens_Press = 0.0;
  Total_Sens_Temp  = 0.0;
	//	Total_Sens_Far = 0.0;
  
  RHO_INDEX = solver_container[TNE2_SOL]->node[0]->GetRhoIndex();
  H_INDEX   = solver_container[TNE2_SOL]->node[0]->GetHIndex();
  P_INDEX   = solver_container[TNE2_SOL]->node[0]->GetPIndex();
  VEL_INDEX = solver_container[TNE2_SOL]->node[0]->GetVelIndex();
  
	/*--- Compute surface sensitivity ---*/
  
  /*--- Loop over boundary markers to select those for Euler walls ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_Boundary(iMarker) == EULER_WALL) {
      
      /*--- Loop over points on the surface to store the auxiliary variable ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Psi     = node[iPoint]->GetSolution();
          U       = solver_container[TNE2_SOL]->node[iPoint]->GetSolution();
          H       = solver_container[TNE2_SOL]->node[iPoint]->GetEnthalpy();
          rho     = solver_container[TNE2_SOL]->node[iPoint]->GetDensity();
          conspsi = 0.0;
          for (iVar = 0; iVar < nSpecies+nDim; iVar++)
            conspsi += U[iVar]*Psi[iVar];
          conspsi   += rho*H*Psi[nSpecies+nDim];
          conspsi   += U[nSpecies+nDim+1]*Psi[nSpecies+nDim+1];
          
          node[iPoint]->SetAuxVar(conspsi);
          
          /*--- Also load the auxiliary variable for first neighbors ---*/
          for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
            Neigh = geometry->node[iPoint]->GetPoint(iNeigh);
            Psi   = node[Neigh]->GetSolution();
            U     = solver_container[TNE2_SOL]->node[Neigh]->GetSolution();
            H     = solver_container[TNE2_SOL]->node[Neigh]->GetEnthalpy();
            
            rho     = solver_container[TNE2_SOL]->node[Neigh]->GetDensity();
            conspsi = 0.0;
            for (iVar = 0; iVar < nSpecies+nDim; iVar++)
              conspsi += U[iVar]*Psi[iVar];
            conspsi   += rho*H*Psi[nSpecies+nDim];
            conspsi   += U[nSpecies+nDim+1]*Psi[nSpecies+nDim+1];
            node[Neigh]->SetAuxVar(conspsi);
          }
        }
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
          
          PrimVar_Grad = solver_container[TNE2_SOL]->node[iPoint]->GetGradient_Primitive();
          ConsVar_Grad = solver_container[TNE2_SOL]->node[iPoint]->GetGradient();
          ConsPsi_Grad = node[iPoint]->GetAuxVarGradient();
          ConsPsi = node[iPoint]->GetAuxVar();          
          
          d_press = 0.0; grad_v = 0.0; v_gradconspsi = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            
            /*-- Retrieve the value of the pressure gradient ---*/
            d_press += d[iDim]*PrimVar_Grad[P_INDEX][iDim];
            
            /*-- Retrieve the value of the velocity gradient ---*/
            grad_v += PrimVar_Grad[VEL_INDEX+iDim][iDim]*ConsPsi;
            
            /*-- Retrieve the value of the theta gradient ---*/
            v_gradconspsi += solver_container[TNE2_SOL]->node[iPoint]->GetVelocity(iDim) * ConsPsi_Grad[iDim];
          }
          
          /*--- Compute sensitivity for each surface point ---*/
          CSensitivity[iMarker][iVertex] = (d_press + grad_v + v_gradconspsi) * Area;
          Sens_Geo[iMarker] -= CSensitivity[iMarker][iVertex] * Area;
        }
      }
      Total_Sens_Geo += Sens_Geo[iMarker];
    }
  }
  
  /*--- Farfield Sensitivity (Mach, AoA, Press, Temp) ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if (config->GetMarker_All_Boundary(iMarker) == FAR_FIELD) {
      
      Sens_Mach[iMarker]  = 0.0;
      Sens_AoA[iMarker]   = 0.0;
      Sens_Press[iMarker] = 0.0;
      Sens_Temp[iMarker]  = 0.0;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          Psi      = node[iPoint]->GetSolution();
          U        = solver_container[TNE2_SOL]->node[iPoint]->GetSolution();
          V        = solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar();
          dPdU     = solver_container[TNE2_SOL]->node[iPoint]->GetdPdU();
          Normal   = geometry->vertex[iMarker][iVertex]->GetNormal();
          Mach_Inf = config->GetMach_FreeStreamND();
          
          rho = V[RHO_INDEX];
          rhou = U[nSpecies];
          rhov = U[nSpecies+1];
          if (nDim == 2) {
            rhow   = 0.0;
            rhoE   = U[nSpecies+nDim];
            rhoEve = U[nSpecies+nDim+1];
          }
          else {
            rhow   = U[nSpecies+2];
            rhoE   = U[nSpecies+nDim];
            rhoEve = U[nSpecies+nDim+1];
          }
          p = V[P_INDEX];
          H = V[H_INDEX];
          
          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++)
            UnitNormal[iDim] = -Normal[iDim]/Area;
          
          /*--- Get the inviscid projected Jacobian ---*/
          numerics->GetInviscidProjJac(U, V, dPdU, UnitNormal, 1.0, Jacobian_j);
          
          /*--- Take the transpose & integrate over dual-face area ---*/
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_jj[iVar][jVar] = Jacobian_j[jVar][iVar]*Area;
          
          
          /*--- Mach number sensitivity ---*/
          for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
            USens[iSpecies] = 0.0;
          USens[nSpecies]   = rhou/Mach_Inf;
          USens[nSpecies+1] = rhov/Mach_Inf;
          if (nDim == 2)
            USens[nSpecies+2] = Gamma*Mach_Inf*p;
          else {
            USens[nSpecies+2] = rhow/Mach_Inf;
            USens[nSpecies+3] = Gamma*Mach_Inf*p; }
          for (iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++) {
              Sens_Mach[iMarker] += Psi[iVar]*Jacobian_j[jVar][iVar]*USens[jVar];
            }
          }
          
          /*--- AoA sensitivity ---*/
          for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
            USens[iSpecies] = 0.0;
          if (nDim == 2) {
            USens[nSpecies]   = -rhov;
            USens[nSpecies+1] = rhou;
            USens[nSpecies+2] = 0.0;
          } else {
            USens[nSpecies]   = -rhow;
            USens[nSpecies+1] = 0.0;
            USens[nSpecies+2] = rhou;
            USens[nSpecies+3] = 0.0;
          }
          for (iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++) {
              Sens_AoA[iMarker] += Psi[iVar]*Jacobian_j[jVar][iVar]*USens[jVar];
            }
          }
          
//          /*--- Pressure sensitivity ---*/
//          USens[0] = r/p; USens[1] = ru/p; USens[2] = rv/p;
//          if (nDim == 2) { USens[3] = rE/p; }
//          else { USens[3] = rw/p; USens[4] = rE/p; }
//          for (iPos = 0; iPos < nVar; iPos++) {
//            for (jPos = 0; jPos < nVar; jPos++) {
//              Sens_Press[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
//            }
//          }
//          
//          /*--- Temperature sensitivity ---*/
//          T = p/(r*Gas_Constant);
//          USens[0] = -r/T; USens[1] = 0.5*ru/T; USens[2] = 0.5*rv/T;
//          if (nDim == 2) { USens[3] = (ru*ru + rv*rv + rw*rw)/(r*T); }
//          else { USens[3] = 0.5*rw/T; USens[4] = (ru*ru + rv*rv + rw*rw)/(r*T); }
//          for (iPos = 0; iPos < nVar; iPos++) {
//            for (jPos = 0; jPos < nVar; jPos++) {
//              Sens_Temp[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
//            }
//          }
        }
      }
      Total_Sens_Mach -= Sens_Mach[iMarker];
      Total_Sens_AoA -= Sens_AoA[iMarker];
      Total_Sens_Press -= Sens_Press[iMarker];
      Total_Sens_Temp -= Sens_Temp[iMarker];
    }
  }
  
  /*--- Explicit contribution from objective function quantity ---*/
//  for (iMarker = 0; iMarker < nMarker; iMarker++) {
//    
//    if (config->GetMarker_All_Boundary(iMarker) == EULER_WALL) {
//      
//      Sens_Mach[iMarker]  = 0.0;
//      Sens_AoA[iMarker]   = 0.0;
//      Sens_Press[iMarker] = 0.0;
//      Sens_Temp[iMarker]  = 0.0;
//      
//      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
//        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
//        
//        if (geometry->node[iPoint]->GetDomain()) {
//          
//          U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
//          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
//          p = solver_container[FLOW_SOL]->node[iPoint]->GetPressure();
//          
//          Mach_Inf   = config->GetMach_FreeStreamND();
//          if (grid_movement) Mach_Inf = config->GetMach_Motion();
//          
//          d = node[iPoint]->GetForceProj_Vector();
//          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
//          Area = sqrt(Area);
//          for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;
//          
//          /*--- Mach number sensitivity ---*/
//          for (iPos = 0; iPos < nDim; iPos++) Dd[iPos] = -(2/Mach_Inf)*d[iPos];
//          for (iPos = 0; iPos < nDim; iPos++) Sens_Mach[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];
//          
//          /*--- AoA sensitivity ---*/
//          if (nDim == 2) {
//            D[0][0] = 0.0; D[0][1] = -1.0;
//            D[1][0] = 1.0; D[1][1] = 0.0;
//          }
//          else {
//            D[0][0] = 0.0; D[0][1] = 0.0; D[0][2] = -1.0;
//            D[1][0] = 0.0; D[1][1] = 0.0; D[1][2] = 0.0;
//            D[2][0] = 1.0; D[2][1] = 0.0; D[2][2] = 0.0;
//          }
//          
//          for (iPos = 0; iPos < nDim; iPos++) Dd[iPos] = 0.0;
//          for (iPos = 0; iPos < nDim; iPos++) {
//            for (jPos = 0; jPos < nDim; jPos++)
//              Dd[iPos] += D[iPos][jPos]*d[jPos];
//          }
//          
//          for (iPos = 0; iPos < nDim; iPos++)
//            Sens_AoA[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];
//          
//          /*--- Pressure sensitivity ---*/
//          for (iPos = 0; iPos<nDim; iPos++) Dd[iPos] = -(1/p)*d[iPos];
//          for (iPos = 0; iPos<nDim; iPos++)
//            Sens_Press[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];
//          
//          /*--- Temperature sensitivity ---*/
//          for (iPos = 0; iPos<nDim; iPos++) Dd[iPos] = 0.0;
//          for (iPos = 0; iPos<nDim; iPos++)
//            Sens_Temp[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];
//          
//        }
//      }
//      
//      Total_Sens_Mach   += Sens_Mach[iMarker];
//      Total_Sens_AoA    += Sens_AoA[iMarker];
//      Total_Sens_Press  += Sens_Press[iMarker];
//      Total_Sens_Temp   += Sens_Temp[iMarker];
//      
//    }
//  }
  
  delete [] USens;
  delete [] UnitNormal;
}

void CAdjTNE2EulerSolver::BC_Euler_Wall(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) {
  
  bool implicit;
  unsigned short iDim, iVar, jVar;
  unsigned long iVertex, iPoint;
	double *d, *Normal, Area, *UnitNormal, *Coord;
  double *Psi, *Psi_Aux, phin, bcn;
  double *U, *V, *dPdU;
  
  /*--- Set booleans from config ---*/
	implicit   = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   (solver_container[TNE2_SOL]->node[0]->GetRhosIndex()   );
  numerics->SetRhoIndex    (solver_container[TNE2_SOL]->node[0]->GetRhoIndex()    );
  numerics->SetPIndex      (solver_container[TNE2_SOL]->node[0]->GetPIndex()      );
  numerics->SetTIndex      (solver_container[TNE2_SOL]->node[0]->GetTIndex()      );
  numerics->SetTveIndex    (solver_container[TNE2_SOL]->node[0]->GetTveIndex()    );
  numerics->SetVelIndex    (solver_container[TNE2_SOL]->node[0]->GetVelIndex()    );
  numerics->SetHIndex      (solver_container[TNE2_SOL]->node[0]->GetHIndex()      );
  numerics->SetAIndex      (solver_container[TNE2_SOL]->node[0]->GetAIndex()      );
  numerics->SetRhoCvtrIndex(solver_container[TNE2_SOL]->node[0]->GetRhoCvtrIndex());
  numerics->SetRhoCvveIndex(solver_container[TNE2_SOL]->node[0]->GetRhoCvveIndex());

  /*--- Initialize ---*/
  d = NULL;
  
  /*--- Allocate arrays ---*/
	UnitNormal    = new double[nDim];
	Psi           = new double[nVar];
  
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Get node information ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Coord = geometry->node[iPoint]->GetCoord();
      
      /*--- Compute geometry parameters ---*/
			Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
			Area = sqrt(Area);
			for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Create a copy of the adjoint solution ---*/
      Psi_Aux = node[iPoint]->GetSolution();
      for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];
      
			/*--- Set the direct solution ---*/
			U    = solver_container[TNE2_SOL]->node[iPoint]->GetSolution();
      V    = solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar();
      dPdU = solver_container[TNE2_SOL]->node[iPoint]->GetdPdU();
      
      /*--- Get the force projection vector, d ---*/
      d = node[iPoint]->GetForceProj_Vector();
      
      /*--- Compute projections ---*/
      bcn = 0.0;
      phin = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        bcn     += d[iDim]*UnitNormal[iDim];
        phin    += Psi[nSpecies+iDim]*UnitNormal[iDim];
      }

      /*--- Introduce the boundary condition ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Psi[nSpecies+iDim] -= ( phin - bcn ) * UnitNormal[iDim];

      numerics->GetInviscidProjJac(U, V, dPdU, UnitNormal, 1.0, Jacobian_i);
      
      /*--- Flux of the Euler wall: (Adotn)^T * Psi ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
        for (jVar = 0; jVar < nVar; jVar++) {
          Residual[iVar] += Jacobian_i[jVar][iVar]*Psi[jVar]*Area;
        }
      }

      if (implicit)
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_ii[iVar][jVar] = Jacobian_i[jVar][iVar]*Area;
      
      /*--- Update residual ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}
	delete [] UnitNormal;
	delete [] Psi;
}

void CAdjTNE2EulerSolver::BC_Sym_Plane(CGeometry *geometry,
                                       CSolver **solver_container,
                                       CNumerics *conv_numerics,
                                       CNumerics *visc_numerics,
                                       CConfig *config,
                                       unsigned short val_marker) {
	
  bool implicit;
  unsigned short iDim, iVar, jVar;
  unsigned long iVertex, iPoint;
	double *Normal, Area, *UnitNormal, *Coord;
  double *Psi, *Psi_Aux, phin;
  double *U, *V, *dPdU;
  
  /*--- Set booleans from config ---*/
	implicit   = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  conv_numerics->SetRhosIndex   (solver_container[TNE2_SOL]->node[0]->GetRhosIndex()   );
  conv_numerics->SetRhoIndex    (solver_container[TNE2_SOL]->node[0]->GetRhoIndex()    );
  conv_numerics->SetPIndex      (solver_container[TNE2_SOL]->node[0]->GetPIndex()      );
  conv_numerics->SetTIndex      (solver_container[TNE2_SOL]->node[0]->GetTIndex()      );
  conv_numerics->SetTveIndex    (solver_container[TNE2_SOL]->node[0]->GetTveIndex()    );
  conv_numerics->SetVelIndex    (solver_container[TNE2_SOL]->node[0]->GetVelIndex()    );
  conv_numerics->SetHIndex      (solver_container[TNE2_SOL]->node[0]->GetHIndex()      );
  conv_numerics->SetAIndex      (solver_container[TNE2_SOL]->node[0]->GetAIndex()      );
  conv_numerics->SetRhoCvtrIndex(solver_container[TNE2_SOL]->node[0]->GetRhoCvtrIndex());
  conv_numerics->SetRhoCvveIndex(solver_container[TNE2_SOL]->node[0]->GetRhoCvveIndex());
  
  /*--- Allocate arrays ---*/
	UnitNormal    = new double[nDim];
	Psi           = new double[nVar];
  
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Get node information ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Coord = geometry->node[iPoint]->GetCoord();
      
      /*--- Compute geometry parameters ---*/
			Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
			Area = sqrt(Area);
			for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Create a copy of the adjoint solution ---*/
      Psi_Aux = node[iPoint]->GetSolution();
      for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];
      
			/*--- Set the direct solution ---*/
			U = solver_container[TNE2_SOL]->node[iPoint]->GetSolution();
      V = solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar();
      dPdU = solver_container[TNE2_SOL]->node[iPoint]->GetdPdU();
      
      /*--- Compute projections ---*/
      phin = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        phin    += Psi[nSpecies+iDim]*UnitNormal[iDim];

      /*--- Introduce the boundary condition ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Psi[nSpecies+iDim] -= phin * UnitNormal[iDim];
      
      conv_numerics->GetInviscidProjJac(U, V, dPdU, UnitNormal, 1.0, Jacobian_i);
      
      /*--- Flux of the Euler wall: (Adotn)^T * Psi ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
        for (jVar = 0; jVar < nVar; jVar++) {
          Residual[iVar] += Jacobian_i[jVar][iVar]*Psi[jVar]*Area;
        }
      }
      
      if (implicit)
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_ii[iVar][jVar] = Jacobian_i[jVar][iVar]*Area;
      
      /*--- Update residual ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}
	delete [] UnitNormal;
	delete [] Psi;
}


void CAdjTNE2EulerSolver::BC_Far_Field(CGeometry *geometry,
                                       CSolver **solver_container,
                                       CNumerics *conv_numerics,
                                       CNumerics *visc_numerics,
                                       CConfig *config,
                                       unsigned short val_marker) {
  bool implicit;
	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	double *Normal, *U_domain, *U_infty, *V_domain, *V_infty;
  double *Psi_domain, *Psi_infty;
  
  /*--- Set booleans from config settings ---*/
	implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  /*--- Allocate arrays ---*/
	Normal     = new double[nDim];
	Psi_domain = new double[nVar];
  Psi_infty  = new double[nVar];
  
  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  conv_numerics->SetRhosIndex   ( solver_container[TNE2_SOL]->node[0]->GetRhosIndex()    );
  conv_numerics->SetRhoIndex    ( solver_container[TNE2_SOL]->node[0]->GetRhoIndex()     );
  conv_numerics->SetPIndex      ( solver_container[TNE2_SOL]->node[0]->GetPIndex()       );
  conv_numerics->SetTIndex      ( solver_container[TNE2_SOL]->node[0]->GetTIndex()       );
  conv_numerics->SetTveIndex    ( solver_container[TNE2_SOL]->node[0]->GetTveIndex()     );
  conv_numerics->SetVelIndex    ( solver_container[TNE2_SOL]->node[0]->GetVelIndex()     );
  conv_numerics->SetHIndex      ( solver_container[TNE2_SOL]->node[0]->GetHIndex()       );
  conv_numerics->SetAIndex      ( solver_container[TNE2_SOL]->node[0]->GetAIndex()       );
  conv_numerics->SetRhoCvtrIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvtrIndex() );
  conv_numerics->SetRhoCvveIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvveIndex() );
  
	/*--- Loop over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- If the node belongs to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
			conv_numerics->SetNormal(Normal);
      
      /*--- Retrieve adjoint solution from boundary and free stream ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
        Psi_infty[iVar] = 0.0;
      }
      
			/*--- Retrieve solution from boundary & free stream ---*/
      U_domain = solver_container[TNE2_SOL]->node[iPoint]->GetSolution();
      U_infty  = solver_container[TNE2_SOL]->node_infty->GetSolution();
      V_domain = solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar();
      V_infty  = solver_container[TNE2_SOL]->node_infty->GetPrimVar();
    
      /*--- Pass conserved & primitive variables to CNumerics ---*/
			conv_numerics->SetConservative(U_domain, U_infty);
      conv_numerics->SetPrimitive(V_domain,V_infty);
      
      /*--- Pass supplementary information to CNumerics ---*/
      conv_numerics->SetdPdU(solver_container[TNE2_SOL]->node[iPoint]->GetdPdU(),
                             solver_container[TNE2_SOL]->node_infty->GetdPdU());
      conv_numerics->SetdTdU(solver_container[TNE2_SOL]->node[iPoint]->GetdTdU(),
                             solver_container[TNE2_SOL]->node_infty->GetdTdU());
      conv_numerics->SetdTvedU(solver_container[TNE2_SOL]->node[iPoint]->GetdTvedU(),
                               solver_container[TNE2_SOL]->node_infty->GetdTvedU());
      
      /*--- Pass adjoint solution to CNumerics ---*/
      conv_numerics->SetAdjointVar(Psi_domain, Psi_infty);
      
			/*--- Compute the convective residual (and Jacobian) ---*/
      conv_numerics->ComputeResidual(Residual_i, Residual_j,
                                     Jacobian_ii, Jacobian_ij,
                                     Jacobian_ji, Jacobian_jj, config);
      
      /*--- Apply contribution to the linear system ---*/
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
    }
  }
	delete [] Normal;
	delete [] Psi_domain;
  delete [] Psi_infty;
}

CAdjTNE2NSSolver::CAdjTNE2NSSolver(void) : CAdjTNE2EulerSolver() { }

CAdjTNE2NSSolver::CAdjTNE2NSSolver(CGeometry *geometry,
                                   CConfig *config,
                                   unsigned short iMesh) : CAdjTNE2EulerSolver() {

  bool restart;
	unsigned short iDim, iMarker, iSpecies, iVar, nLineLets;
  unsigned long iPoint, index, iVertex;
  double dull_val;
  string text_line, mesh_filename;
  string filename, AdjExt;
	ifstream restart_file;
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  /*--- Set booleans for solver settings ---*/
  restart = config->GetRestart();
  
	/*--- Define constants in the solver structure ---*/
  nSpecies     = config->GetnSpecies();
	nDim         = geometry->GetnDim();
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Set the size of the solution array ---*/
	nVar         = nSpecies+nDim+2;
  nPrimVar     = nSpecies+nDim+8;
  nPrimVarGrad = nSpecies+nDim+8;
  
  /*--- Allocate teh CVariable array for each node in the mesh ---*/
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
	if (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT) {
    
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
		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry);
    
    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
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
	PsiRho_Inf = new double [nSpecies];
  Phi_Inf    = new double [nDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    PsiRho_Inf[iSpecies] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Phi_Inf[iDim] = 0.0;
  if (config->GetKind_ObjFunc() == HEAT_LOAD) {
    PsiE_Inf = -1.0;
    PsiEve_Inf = -1.0;
  } else {
    PsiE_Inf = 0.0;
    PsiEve_Inf = 0.0;
  }
  
	if (!restart || geometry->GetFinestMGLevel() == false) {
		/*--- Restart the solution from infinity ---*/
		for (iPoint = 0; iPoint < nPoint; iPoint++)
			node[iPoint] = new CAdjTNE2NSVariable(PsiRho_Inf, Phi_Inf,
                                            PsiE_Inf, PsiEve_Inf,
                                            nDim, nVar, config);
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
      case HEAT_LOAD: AdjExt = "_Q.dat"; break;
		}
		filename.append(AdjExt);
		restart_file.open(filename.data(), ios::in);
    
		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
			cout << "There is no adjoint restart file!! " << filename.data() << "."<< endl;
			exit(1);
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
        point_line >> index;
        for (iDim = 0; iDim < nDim; iDim++)
          point_line >> dull_val;
        for (iVar = 0; iVar < nVar; iVar++)
          point_line >> Solution[iVar];
        
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
  unsigned short iMarker;
  
  if (PsiRho_Inf  != NULL) delete [] PsiRho_Inf;
  if (Phi_Inf     != NULL) delete [] Phi_Inf;
	if (Sens_Mach   != NULL) delete [] Sens_Mach;
	if (Sens_AoA    != NULL) delete [] Sens_AoA;
	if (Sens_Geo    != NULL) delete [] Sens_Geo;
	if (Sens_Press  != NULL) delete [] Sens_Press;
	if (Sens_Temp   != NULL) delete [] Sens_Temp;
	if (iPoint_UndLapl != NULL) delete [] iPoint_UndLapl;
	if (jPoint_UndLapl != NULL) delete [] jPoint_UndLapl;
  
	if (CSensitivity != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete CSensitivity[iMarker];
    delete [] CSensitivity;
  }
  
}


void CAdjTNE2NSSolver::Preprocessing(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CConfig *config,
                                     unsigned short iMesh,
                                     unsigned short iRKStep,
                                     unsigned short RunTime_EqSystem, bool Output) {
  bool implicit, upwind_2nd, center, center_jst, limiter, RightSol;
	unsigned long iPoint, ErrorCounter;
  double SharpEdge_Distance;
  int rank;
  
#ifdef NO_MPI
	rank = MASTER_NODE;
#else
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the flow problem may use different methods). ---*/
  implicit   = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  upwind_2nd = ((config->GetKind_Upwind_AdjTNE2() == ROE_2ND) ||
                     (config->GetKind_Upwind_AdjTNE2() == MSW_2ND) ||
                     (config->GetKind_Upwind_AdjTNE2() == AUSM_2ND) ||
                     (config->GetKind_Upwind_AdjTNE2() == AUSMPWPLUS_2ND));
  center     = (config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_CENTERED);
  center_jst = (config->GetKind_Centered_AdjTNE2() == JST);
  limiter    = (config->GetKind_SlopeLimit() != NONE);
  
  /*--- Initialize marker for tracking non-physical solutions ---*/
  ErrorCounter = 0;
  
	/*--- Residual initialization ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Get the distance form a sharp edge ---*/
    SharpEdge_Distance = geometry->node[iPoint]->GetSharpEdge_Distance();
    
    /*--- Set the primitive variables incompressible and compressible
     adjoint variables ---*/
    RightSol = node[iPoint]->SetPrimVar_Compressible(SharpEdge_Distance,
                                                      false, config);
    if (!RightSol) ErrorCounter++;
    
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
  
  /*--- Error message ---*/
#ifndef NO_MPI
  unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
#ifdef WINDOWS
  MPI_Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI::UNSIGNED_LONG, MPI::SUM);
#endif
#endif
  if ((ErrorCounter != 0) && (rank == MASTER_NODE) && (iMesh == MESH_0))
    cout <<"The solution contains "<< ErrorCounter << " non-physical points." << endl;
  
}

void CAdjTNE2NSSolver::Viscous_Residual(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CNumerics *numerics,
                                        CConfig *config,
                                        unsigned short iMesh,
                                        unsigned short iRKStep) {
  bool implicit;
	unsigned long iPoint, jPoint, iEdge;
  
  /*--- Determine if using implicit time-stepping ---*/
	implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( solver_container[TNE2_SOL]->node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( solver_container[TNE2_SOL]->node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( solver_container[TNE2_SOL]->node[0]->GetPIndex()       );
  numerics->SetTIndex      ( solver_container[TNE2_SOL]->node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( solver_container[TNE2_SOL]->node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( solver_container[TNE2_SOL]->node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( solver_container[TNE2_SOL]->node[0]->GetHIndex()       );
  numerics->SetAIndex      ( solver_container[TNE2_SOL]->node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvveIndex() );
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Identify points on the edge, normal, and neighbors ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[jPoint]->GetCoord());
    
    /*--- Pass the adjoint variables w/o reconstruction to CNumerics ---*/
    numerics->SetAdjointVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    
    /*--- Pass conservative & primitive variables w/o reconstruction ---*/
    numerics->SetConservative(solver_container[TNE2_SOL]->node[iPoint]->GetSolution(),
                              solver_container[TNE2_SOL]->node[jPoint]->GetSolution());
    numerics->SetPrimitive(solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar(),
                           solver_container[TNE2_SOL]->node[jPoint]->GetPrimVar());
    
    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU(  solver_container[TNE2_SOL]->node[iPoint]->GetdPdU(),
                        solver_container[TNE2_SOL]->node[jPoint]->GetdPdU());
    numerics->SetdTdU(  solver_container[TNE2_SOL]->node[iPoint]->GetdTdU(),
                        solver_container[TNE2_SOL]->node[jPoint]->GetdTdU());
    numerics->SetdTvedU(solver_container[TNE2_SOL]->node[iPoint]->GetdTvedU(),
                        solver_container[TNE2_SOL]->node[jPoint]->GetdTvedU());
    
    /*--- Pass transport coefficients to CNumerics ---*/
    numerics->SetDiffusionCoeff(solver_container[TNE2_SOL]->node[iPoint]->GetDiffusionCoeff(),
                                solver_container[TNE2_SOL]->node[jPoint]->GetDiffusionCoeff() );
    numerics->SetLaminarViscosity(solver_container[TNE2_SOL]->node[iPoint]->GetLaminarViscosity(),
                                  solver_container[TNE2_SOL]->node[jPoint]->GetLaminarViscosity() );
    numerics->SetThermalConductivity(solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity(),
                                     solver_container[TNE2_SOL]->node[jPoint]->GetThermalConductivity());
    numerics->SetThermalConductivity_ve(solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve(),
                                        solver_container[TNE2_SOL]->node[jPoint]->GetThermalConductivity_ve() );
    
    /*--- Gradient of Adjoint Variables ---*/
    numerics->SetAdjointVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
    
//    /*--- Compute residual in a non-conservative way, and update ---*/
//    numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
//    
//    /*--- Update adjoint viscous residual ---*/
//    LinSysRes.SubtractBlock(iPoint, Residual_i);
//    LinSysRes.AddBlock(jPoint, Residual_j);
//    
//    if (implicit) {
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
//      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
//      Jacobian.AddBlock(jPoint, iPoint, Jacobian_ji);
//      Jacobian.AddBlock(jPoint, jPoint, Jacobian_jj);
//    }
  }
}

void CAdjTNE2NSSolver::Source_Residual(CGeometry *geometry,
                                       CSolver **solver_container,
                                       CNumerics *numerics,
                                       CNumerics *second_numerics,
                                       CConfig *config, unsigned short iMesh) {
	bool implicit;
  unsigned short iVar, jVar;
  unsigned long iPoint;

  /*--- Numerics notes ---*/
  // numerics -- Container for TNE2 direct source terms (chemistry, etc.)
  // second_numerics -- Container for AdjTNE2 source terms (cons. & visc.)
  //
  // Note: We need to pass the appropriate primitive variable structure to both
  //       the numerics & second_numerics containers.  Chemistry & vib.-el.
  //       sources are handled exactly the same as the direct solver, and we
  //       calculate the source residual by transposing the source Jacobian
  //       matrix & multiply by the adjoint variables.
  
  /*--- Assign booleans from config settings ---*/
  implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( solver_container[TNE2_SOL]->node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( solver_container[TNE2_SOL]->node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( solver_container[TNE2_SOL]->node[0]->GetPIndex()       );
  numerics->SetTIndex      ( solver_container[TNE2_SOL]->node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( solver_container[TNE2_SOL]->node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( solver_container[TNE2_SOL]->node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( solver_container[TNE2_SOL]->node[0]->GetHIndex()       );
  numerics->SetAIndex      ( solver_container[TNE2_SOL]->node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvveIndex() );
  
  second_numerics->SetRhosIndex   ( solver_container[TNE2_SOL]->node[0]->GetRhosIndex()    );
  second_numerics->SetRhoIndex    ( solver_container[TNE2_SOL]->node[0]->GetRhoIndex()     );
  second_numerics->SetPIndex      ( solver_container[TNE2_SOL]->node[0]->GetPIndex()       );
  second_numerics->SetTIndex      ( solver_container[TNE2_SOL]->node[0]->GetTIndex()       );
  second_numerics->SetTveIndex    ( solver_container[TNE2_SOL]->node[0]->GetTveIndex()     );
  second_numerics->SetVelIndex    ( solver_container[TNE2_SOL]->node[0]->GetVelIndex()     );
  second_numerics->SetHIndex      ( solver_container[TNE2_SOL]->node[0]->GetHIndex()       );
  second_numerics->SetAIndex      ( solver_container[TNE2_SOL]->node[0]->GetAIndex()       );
  second_numerics->SetRhoCvtrIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvtrIndex() );
  second_numerics->SetRhoCvveIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvveIndex() );
  
	/*--- Loop over all the points, note that we are supposing that primitive and
	 adjoint gradients have been computed previously ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*---+++ Direct problem +++---*/
		/*--- Set conserved & primitive variables at point i ---*/
		numerics->SetConservative(solver_container[TNE2_SOL]->node[iPoint]->GetSolution(),
                              solver_container[TNE2_SOL]->node[iPoint]->GetSolution());
    numerics->SetPrimitive(solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar(),
                           solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar());
    second_numerics->SetConservative(solver_container[TNE2_SOL]->node[iPoint]->GetSolution(),
                                     solver_container[TNE2_SOL]->node[iPoint]->GetSolution());
    second_numerics->SetPrimitive(solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar(),
                                  solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar());
    
    /*--- Pass the adjoint variables to CNumerics ---*/
    second_numerics->SetAdjointVar(node[iPoint]->GetSolution(),
                                   node[iPoint]->GetSolution());
    
    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU(solver_container[TNE2_SOL]->node[iPoint]->GetdPdU(),
                      solver_container[TNE2_SOL]->node[iPoint]->GetdPdU());
    numerics->SetdTdU(solver_container[TNE2_SOL]->node[iPoint]->GetdTdU(),
                      solver_container[TNE2_SOL]->node[iPoint]->GetdTdU());
    numerics->SetdTvedU(solver_container[TNE2_SOL]->node[iPoint]->GetdTvedU(),
                        solver_container[TNE2_SOL]->node[iPoint]->GetdTvedU());
    second_numerics->SetdPdU(solver_container[TNE2_SOL]->node[iPoint]->GetdPdU(),
                      solver_container[TNE2_SOL]->node[iPoint]->GetdPdU());
    second_numerics->SetdTdU(solver_container[TNE2_SOL]->node[iPoint]->GetdTdU(),
                      solver_container[TNE2_SOL]->node[iPoint]->GetdTdU());
    second_numerics->SetdTvedU(solver_container[TNE2_SOL]->node[iPoint]->GetdTvedU(),
                        solver_container[TNE2_SOL]->node[iPoint]->GetdTvedU());
    
		/*--- Gradient of primitive and adjoint variables ---*/
		numerics->SetPrimVarGradient(solver_container[TNE2_SOL]->node[iPoint]->GetGradient_Primitive(),
                                 solver_container[TNE2_SOL]->node[iPoint]->GetGradient_Primitive());
		numerics->SetAdjointVarGradient(node[iPoint]->GetGradient(),
                                    node[iPoint]->GetGradient());
    second_numerics->SetPrimVarGradient(solver_container[TNE2_SOL]->node[iPoint]->GetGradient_Primitive(),
                                        solver_container[TNE2_SOL]->node[iPoint]->GetGradient_Primitive());
		second_numerics->SetAdjointVarGradient(node[iPoint]->GetGradient(),
                                           node[iPoint]->GetGradient());
    
    /*--- Pass transport coefficients to CNumerics ---*/
    numerics->SetDiffusionCoeff(solver_container[TNE2_SOL]->node[iPoint]->GetDiffusionCoeff(),
                                solver_container[TNE2_SOL]->node[iPoint]->GetDiffusionCoeff());
    numerics->SetLaminarViscosity(solver_container[TNE2_SOL]->node[iPoint]->GetLaminarViscosity(),
                                  solver_container[TNE2_SOL]->node[iPoint]->GetLaminarViscosity());
    numerics->SetThermalConductivity(solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity(),
                                     solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity());
    numerics->SetThermalConductivity_ve(solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve(),
                                        solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve());
    second_numerics->SetDiffusionCoeff(solver_container[TNE2_SOL]->node[iPoint]->GetDiffusionCoeff(),
                                solver_container[TNE2_SOL]->node[iPoint]->GetDiffusionCoeff());
    second_numerics->SetLaminarViscosity(solver_container[TNE2_SOL]->node[iPoint]->GetLaminarViscosity(),
                                  solver_container[TNE2_SOL]->node[iPoint]->GetLaminarViscosity());
    second_numerics->SetThermalConductivity(solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity(),
                                     solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity());
    second_numerics->SetThermalConductivity_ve(solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve(),
                                        solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve());
    
		/*--- Set volume ---*/
		numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    second_numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    
    
    /*--- Compute chemistry source terms ---*/
    numerics->ComputeChemistry(Residual_i, Jacobian_i, config);
    
    /*--- Take the transpose of the source Jacobian matrix ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        Jacobian_ii[iVar][jVar] = Jacobian_i[jVar][iVar];
    
    /*--- Compute the adjoint source term residual (dQ/dU^T * Psi) ---*/
    for (iVar = 0; iVar < nVar; iVar ++) {
      Residual[iVar] = 0.0;
      for (jVar = 0; jVar < nVar; jVar++) {
        Residual[iVar] += Jacobian_ii[iVar][jVar] * node[iPoint]->GetSolution(jVar);
      }
    }
    
    /*--- Subtract Residual (and Jacobian) ---*/
    LinSysRes.SubtractBlock(iPoint, Residual);
    if (implicit)
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
    
    /*--- Compute vibrational relaxation source terms ---*/
    numerics->ComputeVibRelaxation(Residual_i, Jacobian_i, config);
    
    /*--- Take the transpose of the source Jacobian matrix ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        Jacobian_ii[iVar][jVar] = Jacobian_i[jVar][iVar];
    
    /*--- Compute the adjoint source term residual (dQ/dU^T * Psi) ---*/
    for (iVar = 0; iVar < nVar; iVar ++) {
      Residual[iVar] = 0.0;
      for (jVar = 0; jVar < nVar; jVar++) {
        Residual[iVar] += Jacobian_ii[iVar][jVar] * node[iPoint]->GetSolution(jVar);
      }
    }
    
    /*--- Subtract Residual (and Jacobian) ---*/
    LinSysRes.SubtractBlock(iPoint, Residual);
    if (implicit)
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);

    
		/*--- Compute viscous source term residual ---*/
		second_numerics->ComputeSourceViscous(Residual_i, config);
    
    /*--- Add and substract to the residual ---*/
		LinSysRes.AddBlock(iPoint, Residual_i);
	}
}

void CAdjTNE2NSSolver::Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
//	unsigned long iVertex, iPoint;
//	unsigned short iDim, jDim, iMarker;
//	double **PsiVar_Grad, **PrimVar_Grad, div_phi, *Normal, Area,
//	normal_grad_psi5, normal_grad_T, sigma_partial,
//  cp, Laminar_Viscosity, heat_flux_factor, temp_sens;
//  
//	double Gas_Constant = config->GetGas_ConstantND();
//  
//	cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
//  
//	double *UnitNormal = new double[nDim];
//	double *normal_grad_vel = new double[nDim];
//	double *tang_deriv_psi5 = new double[nDim];
//	double *tang_deriv_T = new double[nDim];
//	double **Sigma = new double* [nDim];
//  
//	for (iDim = 0; iDim < nDim; iDim++)
//		Sigma[iDim] = new double [nDim];
//  
//  double *normal_grad_gridvel = new double[nDim];
//  double *normal_grad_v_ux =new double[nDim];
//  double **Sigma_Psi5v = new double* [nDim];
//	for (iDim = 0; iDim < nDim; iDim++)
//		Sigma_Psi5v[iDim] = new double [nDim];
//  double **tau = new double* [nDim];
//	for (iDim = 0; iDim < nDim; iDim++)
//		tau[iDim] = new double [nDim];
//  double *Velocity = new double[nDim];
//  
//  /*--- Compute gradient of adjoint variables on the surface ---*/
//  SetSurface_Gradient(geometry, config);
//  
//  Total_Sens_Geo = 0.0;
//  for (iMarker = 0; iMarker < nMarker; iMarker++) {
//    
//    Sens_Geo[iMarker] = 0.0;
//    
//    if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
//        (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL)) {
//      
//      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
//        
//        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
//        if (geometry->node[iPoint]->GetDomain()) {
//          
//          PsiVar_Grad = node[iPoint]->GetGradient();
//          PrimVar_Grad = solver_container[TNE2_SOL]->node[iPoint]->GetGradient_Primitive();
//          
//          Laminar_Viscosity  = solver_container[TNE2_SOL]->node[iPoint]->GetLaminarViscosity();
//          
//          heat_flux_factor = cp * Laminar_Viscosity / PRANDTL;
//          
//          /*--- Compute face area and the nondimensional normal to the surface ---*/
//          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
//          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) { Area += Normal[iDim]*Normal[iDim]; } Area = sqrt(Area);
//          for (iDim = 0; iDim < nDim; iDim++) { UnitNormal[iDim] = Normal[iDim] / Area; }
//          
//          /*--- Compute the sensitivity related to the temperature ---*/
//          normal_grad_psi5 = 0.0; normal_grad_T = 0.0;
//          for (iDim = 0; iDim < nDim; iDim++) {
//            normal_grad_psi5 += PsiVar_Grad[nVar-1][iDim]*UnitNormal[iDim];
//            normal_grad_T += PrimVar_Grad[0][iDim]*UnitNormal[iDim];
//          }
//          
//          temp_sens = 0.0;
//          if (config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) {
//            
//            /*--- Heat Flux Term: temp_sens = (\partial_tg \psi_5)\cdot (k \partial_tg T) ---*/
//            for (iDim = 0; iDim < nDim; iDim++) {
//              tang_deriv_psi5[iDim] = PsiVar_Grad[nVar-1][iDim] - normal_grad_psi5*UnitNormal[iDim];
//              tang_deriv_T[iDim] = PrimVar_Grad[0][iDim] - normal_grad_T*UnitNormal[iDim];
//            }
//            for (iDim = 0; iDim < nDim; iDim++)
//              temp_sens += heat_flux_factor * tang_deriv_psi5[iDim] * tang_deriv_T[iDim];
//            
//          } else if (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) {
//            
//            /*--- Isothermal Term: temp_sens = - k * \partial_n(\psi_5) * \partial_n(T) ---*/
//            temp_sens = - heat_flux_factor * normal_grad_psi5 * normal_grad_T;
//            
//          }
//          
//          /*--- Term: sigma_partial = \Sigma_{ji} n_i \partial_n v_j ---*/
//          div_phi = 0.0;
//          for (iDim = 0; iDim < nDim; iDim++) {
//            div_phi += PsiVar_Grad[iDim+1][iDim];
//            for (jDim = 0; jDim < nDim; jDim++)
//              Sigma[iDim][jDim] = Laminar_Viscosity * (PsiVar_Grad[iDim+1][jDim]+PsiVar_Grad[jDim+1][iDim]);
//          }
//          for (iDim = 0; iDim < nDim; iDim++)
//            Sigma[iDim][iDim] -= TWO3*Laminar_Viscosity * div_phi;
//          
//          for (iDim = 0; iDim < nDim; iDim++) {
//            normal_grad_vel[iDim] = 0.0;
//            for (jDim = 0; jDim < nDim; jDim++)
//              normal_grad_vel[iDim] += PrimVar_Grad[iDim+1][jDim]*UnitNormal[jDim];
//          }
//          
//          sigma_partial = 0.0;
//          for (iDim = 0; iDim < nDim; iDim++)
//            for (jDim = 0; jDim < nDim; jDim++)
//              sigma_partial += UnitNormal[iDim]*Sigma[iDim][jDim]*normal_grad_vel[jDim];
//          
//          /*--- Compute sensitivity for each surface point ---*/
//          CSensitivity[iMarker][iVertex] = (sigma_partial - temp_sens)*Area;
//          Sens_Geo[iMarker] -= CSensitivity[iMarker][iVertex]*Area;
//        }
//      }
//      Total_Sens_Geo += Sens_Geo[iMarker];
//    }
//  }
//  
//	delete [] UnitNormal;
//	delete [] normal_grad_vel;
//	delete [] tang_deriv_psi5;
//	delete [] tang_deriv_T;
//	for (iDim = 0; iDim < nDim; iDim++)
//		delete Sigma[iDim];
//	delete [] Sigma;
//  
//  delete [] normal_grad_gridvel;
//  delete [] normal_grad_v_ux;
//  for (iDim = 0; iDim < nDim; iDim++)
//		delete Sigma_Psi5v[iDim];
//  for (iDim = 0; iDim < nDim; iDim++)
//		delete tau[iDim];
//	delete [] tau;
//  delete [] Velocity;
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
  
	bool implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
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
      mu_dyn = solver_container[TNE2_SOL]->node[iPoint]->GetLaminarViscosity();
      Temp = solver_container[TNE2_SOL]->node[iPoint]->GetTemperature();
      dVisc_T = 0.0;  // dVisc_T = mu_dyn*(Temp+3.0*mu2)/(2.0*Temp*(Temp+mu2));
      
      /*--- \Sigma_5 ---*/
      Sigma_5 = (Gamma/Cp)*dVisc_T*force_stress;
      
      /*--- Imposition of residuals ---*/
      rho = solver_container[TNE2_SOL]->node[iPoint]->GetDensity();
      pressure = solver_container[TNE2_SOL]->node[iPoint]->GetPressure();
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
        Enthalpy = solver_container[TNE2_SOL]->node[iPoint]->GetEnthalpy();
        
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
        U = solver_container[TNE2_SOL]->node[iPoint]->GetSolution();
        Laminar_Viscosity = solver_container[TNE2_SOL]->node[iPoint]->GetLaminarViscosity();
        Eddy_Viscosity = solver_container[TNE2_SOL]->node[iPoint]->GetEddyViscosity(); // Should be zero at the wall
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


void CAdjTNE2NSSolver::BC_Isothermal_Wall(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CNumerics *conv_numerics,
                                          CNumerics *visc_numerics,
                                          CConfig *config,
                                          unsigned short val_marker) {

  ////////////////////
  bool implicit, heat_flux_obj;
	unsigned long iVertex, iPoint, total_index, Point_Normal;
	unsigned short iDim, iVar, jVar;
  unsigned short RHOS_INDEX, RHO_INDEX;
	double *V, *dPdU, *d, q, dn;
  double phi[3];
  
  //  double *U, mu_dyn, Temp, dVisc_T, rho, pressure, div_phi,
  //  force_stress, Sigma_5, **PsiVar_Grad, phi[3];
  //  double phis1, phis2, sq_vel, ProjVel, Enthalpy, *GridVel, phi_u, d_n;
  //  double Energy, ViscDens, XiDens, Density, SoundSpeed, Pressure, dPhiE_dn, Laminar_Viscosity, Eddy_Viscosity,
  //  Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz,
  //  Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5,
  //  Sigma_yz5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;
  
  double *Psi, **Tau, *Velocity, *Normal, **GradPhi, *GradPsiE;
  
  /*--- Set booleans from CConfig specifications ---*/
  implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  heat_flux_obj  = (config->GetKind_ObjFunc() == HEAT_FLUX);
  
  /*--- Allocate arrays ---*/
  Psi = new double[nVar];
	Tau = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Tau[iDim] = new double [nDim];
  Velocity = new double[nDim];
  Normal = new double[nDim];
  GradPhi = new double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    GradPhi[iDim] = new double [nDim];
  GradPsiE = new double [nDim];
  
  /*--- Get primitive vector locators ---*/
  RHOS_INDEX = solver_container[TNE2_SOL]->node[0]->GetRhosIndex();
  RHO_INDEX  = solver_container[TNE2_SOL]->node[0]->GetRhoIndex();
  
  
  /*--- Loop over all boundary points ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Get node and neighbor information ---*/
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv_i[iVar] = 0.0;
        Res_Visc_i[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar ++)
            Jacobian_ii[iVar][jVar] = 0.0;
        }
      }
      
      /*--- Retrieve adjoint solution at the boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Psi[iVar] = node[iPoint]->GetSolution(iVar);
      
      /*--- Retrieve primitive variables at the boundary node ---*/
      V = solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar();
      dPdU = solver_container[TNE2_SOL]->node[iPoint]->GetdPdU();
      
			/*--- Normal vector for this vertex ---*/
      // Note: Convention is outward facing normal
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Get the force projection vector (based on the objective function) ---*/
			d = node[iPoint]->GetForceProj_Vector();
      
      /*--- Apply the momentum boundary condition ---*/
      dn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        phi[iDim] = d[iDim];
        dn += d[iDim]*Normal[iDim];
      }
      
      /*--- Apply the B.C. to the linear system ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
      node[iPoint]->SetVel_ResTruncError_Zero();
			for (iDim = 0; iDim < nDim; iDim++)
				node[iPoint]->SetSolution_Old(nSpecies+iDim, phi[iDim]);
			if (implicit) {
				for (iVar = nSpecies; iVar < nSpecies+nDim; iVar++) {
					total_index = iPoint*nVar+iVar;
					Jacobian.DeleteValsRowi(total_index);
				}
			}
      
      /*--- Apply the energy & vib. el energy boundary condition ---*/
      if (heat_flux_obj) {
        // Note: This is the derivative of our objective function j = kdndT
        //       for heat flux with a negative sign from the formulation of the
        //       adjoint boundary conditions.
        q = -1.0;
      } else {
        q = 0.0;
      }
      
      /*--- Apply the boundary condition to the linear system ---*/
      LinSysRes.SetBlock_Zero(iPoint, nSpecies+nDim);
      LinSysRes.SetBlock_Zero(iPoint, nSpecies+nDim+1);
      node[iPoint]->SetSolution_Old(nSpecies+nDim,   q);
      node[iPoint]->SetSolution_Old(nSpecies+nDim+1, q);
      if (implicit) {
        iVar = nSpecies+nDim;
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
        Jacobian.DeleteValsRowi(total_index+1);
      }
      
      /*--- Determine contribution to adjoint density ---*/
//      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//        Res_Conv_i[iSpecies] = -dn*dPdU[iSpecies];
//      }
//      
//      /*--- Update convective and viscous residuals ---*/
//      LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
//      LinSysRes.SubtractBlock(iPoint, Res_Visc_i);
//      if (implicit) {
//        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
//      }
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

  
  
  //////////////////////////
//  bool implicit;
//	unsigned long iVertex, iPoint, total_index;
//	unsigned short iVar, iDim;
//	double *Normal, *U_domain, *U_infty, *V_domain, *V_infty;
//  double *Psi_domain, *Psi_infty;
//  double phi[3], *d, dn, q;
//  
//  /*--- Set booleans from config settings ---*/
//	implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
//  
//  /*--- Allocate arrays ---*/
//	Normal     = new double[nDim];
//	Psi_domain = new double[nVar];
//  Psi_infty  = new double[nVar];
//  
//  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
//  conv_numerics->SetRhosIndex   ( solver_container[TNE2_SOL]->node[0]->GetRhosIndex()    );
//  conv_numerics->SetRhoIndex    ( solver_container[TNE2_SOL]->node[0]->GetRhoIndex()     );
//  conv_numerics->SetPIndex      ( solver_container[TNE2_SOL]->node[0]->GetPIndex()       );
//  conv_numerics->SetTIndex      ( solver_container[TNE2_SOL]->node[0]->GetTIndex()       );
//  conv_numerics->SetTveIndex    ( solver_container[TNE2_SOL]->node[0]->GetTveIndex()     );
//  conv_numerics->SetVelIndex    ( solver_container[TNE2_SOL]->node[0]->GetVelIndex()     );
//  conv_numerics->SetHIndex      ( solver_container[TNE2_SOL]->node[0]->GetHIndex()       );
//  conv_numerics->SetAIndex      ( solver_container[TNE2_SOL]->node[0]->GetAIndex()       );
//  conv_numerics->SetRhoCvtrIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvtrIndex() );
//  conv_numerics->SetRhoCvveIndex( solver_container[TNE2_SOL]->node[0]->GetRhoCvveIndex() );
//  
//	/*--- Loop over all the vertices ---*/
//	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//    
//		/*--- If the node belongs to the domain ---*/
//		if (geometry->node[iPoint]->GetDomain()) {
//      
//			/*--- Set the normal vector ---*/
//			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
//			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
//			conv_numerics->SetNormal(Normal);
//      
//      /*--- Retrieve adjoint solution from boundary and free stream ---*/
//      for (iVar = 0; iVar < nSpecies; iVar++) {
//        Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
//        Psi_infty[iVar] = Psi_domain[iVar];
//      }
//      for (iVar = nSpecies; iVar < nVar; iVar++) {
//        Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
//        Psi_infty[iVar] = 0.0;
//      }
//      
//			/*--- Retrieve solution from boundary & free stream ---*/
//      U_domain = solver_container[TNE2_SOL]->node[iPoint]->GetSolution();
//      U_infty  = solver_container[TNE2_SOL]->node_infty->GetSolution();
//      V_domain = solver_container[TNE2_SOL]->node[iPoint]->GetPrimVar();
//      V_infty  = solver_container[TNE2_SOL]->node_infty->GetPrimVar();
//      
//      /*--- Pass conserved & primitive variables to CNumerics ---*/
//			conv_numerics->SetConservative(U_domain, U_infty);
//      conv_numerics->SetPrimitive(V_domain,V_infty);
//      
//      /*--- Pass supplementary information to CNumerics ---*/
//      conv_numerics->SetdPdU(solver_container[TNE2_SOL]->node[iPoint]->GetdPdU(),
//                             solver_container[TNE2_SOL]->node_infty->GetdPdU());
//      conv_numerics->SetdTdU(solver_container[TNE2_SOL]->node[iPoint]->GetdTdU(),
//                             solver_container[TNE2_SOL]->node_infty->GetdTdU());
//      conv_numerics->SetdTvedU(solver_container[TNE2_SOL]->node[iPoint]->GetdTvedU(),
//                               solver_container[TNE2_SOL]->node_infty->GetdTvedU());
//      
//      /*--- Pass adjoint solution to CNumerics ---*/
//      conv_numerics->SetAdjointVar(Psi_domain, Psi_infty);
//      
//			/*--- Compute the convective residual (and Jacobian) ---*/
//      conv_numerics->ComputeResidual(Residual_i, Residual_j,
//                                     Jacobian_ii, Jacobian_ij,
//                                     Jacobian_ji, Jacobian_jj, config);
//      
//      /*--- Apply contribution to the linear system ---*/
//      LinSysRes.SubtractBlock(iPoint, Residual_i);
//      if (implicit)
//        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
//      
//      /*--- Get the force projection vector (based on the objective function) ---*/
//			d = node[iPoint]->GetForceProj_Vector();
//
//      /*--- Apply the momentum boundary condition ---*/
//      for (iDim = 0; iDim < nDim; iDim++) {
//        phi[iDim] = d[iDim];
//      }
//
//      /*--- Apply the B.C. to the linear system ---*/
//      for (iDim = 0; iDim < nDim; iDim++)
//        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
//			for (iDim = 0; iDim < nDim; iDim++)
//				node[iPoint]->SetSolution_Old(nSpecies+iDim, phi[iDim]);
//			if (implicit) {
//				for (iVar = nSpecies; iVar < nSpecies+nDim; iVar++) {
//					total_index = iPoint*nVar+iVar;
//					Jacobian.DeleteValsRowi(total_index);
//				}
//			}
//
//      /*--- Set the energy B.C. ---*/
//      q = 0.0;
//
//      /*--- Apply the boundary condition to the linear system ---*/
//      LinSysRes.SetBlock_Zero(iPoint, nSpecies+nDim);
//      LinSysRes.SetBlock_Zero(iPoint, nSpecies+nDim+1);
//      node[iPoint]->SetSolution_Old(nSpecies+nDim,   q);
//      node[iPoint]->SetSolution_Old(nSpecies+nDim+1, q);
//      if (implicit) {
//        iVar = nSpecies+nDim;
//        total_index = iPoint*nVar+iVar;
//        Jacobian.DeleteValsRowi(total_index);
//        Jacobian.DeleteValsRowi(total_index+1);
//      }
//    }
//  }
//	delete [] Normal;
//	delete [] Psi_domain;
//  delete [] Psi_infty;
  
  
//      /*--- Additional contributions to adjoint density and energy (weak imposition) ---*/
//      /*--- Components of the effective and adjoint stress tensors ---*/
////      PsiVar_Grad = node[iPoint]->GetGradient();
////      div_phi = 0;
////      for (iDim = 0; iDim < nDim; iDim++) {
////        div_phi += PsiVar_Grad[iDim+1][iDim];
////        for (jDim = 0; jDim < nDim; jDim++)
////          Tau[iDim][jDim] = (PsiVar_Grad[iDim+1][jDim]+PsiVar_Grad[jDim+1][iDim]);
////      }
////      for (iDim = 0; iDim < nDim; iDim++)
////        Tau[iDim][iDim] -= TWO3*div_phi;
////      
////      /*--- force_stress = n_i \Tau_{ij} d_j ---*/
////      force_stress = 0.0;
////      for (iDim = 0; iDim < nDim; iDim++)
////        for (jDim = 0; jDim < nDim; jDim++)
////          force_stress += Normal[iDim]*Tau[iDim][jDim]*d[jDim];
////      
////      /*--- \partial \mu_dyn \partial T ---*/
////      mu_dyn = solver_container[TNE2_SOL]->node[iPoint]->GetLaminarViscosity();
////      Temp = solver_container[TNE2_SOL]->node[iPoint]->GetTemperature();
////      dVisc_T = 0.0;  // dVisc_T = mu_dyn*(Temp+3.0*mu2)/(2.0*Temp*(Temp+mu2));
////      
////      /*--- \Sigma_5 ---*/
////      Sigma_5 = (Gamma/Cp)*dVisc_T*force_stress;
////      
////      /*--- Imposition of residuals ---*/
////      rho = solver_container[TNE2_SOL]->node[iPoint]->GetDensity();
////      pressure = solver_container[TNE2_SOL]->node[iPoint]->GetPressure(false);
////      Res_Conv_i[0] = pressure*Sigma_5/(Gamma_Minus_One*rho*rho);
////      
////      /*--- Flux contribution and Jacobian contributions for moving
////       walls. Note that these are only for the adjoint density and
////       adjoint energy equations (the adjoint vel. uses a strong BC). ---*/
////      if (grid_movement) {
////        
////        /*--- Get the appropriate grid velocity at this node ---*/
////        if (grid_movement)
////          GridVel = geometry->node[iPoint]->GetGridVel();
////        
////        /*--- Get the enthalpy from the direct solution ---*/
////        Enthalpy = solver_container[TNE2_SOL]->node[iPoint]->GetEnthalpy();
////        
////        /*--- Compute projections, velocity squared divided by two, and
////         other inner products. Note that we are imposing v = u_wall from
////         the direct problem and that phi = d - \psi_5 * v ---*/
////        ProjVel = 0.0; sq_vel = 0.0; phi_u = 0.0; d_n = 0.0;
////        phis1 = 0.0; phis2 = Psi[0] + Enthalpy * Psi[nVar-1];
////        for (iDim = 0; iDim < nDim; iDim++) {
////          ProjVel += GridVel[iDim]*Normal[iDim];
////          sq_vel  += 0.5*GridVel[iDim]*GridVel[iDim];
////          phis1   += Normal[iDim]*phi[iDim];
////          phis2   += GridVel[iDim]*phi[iDim];
////          phi_u   += GridVel[iDim]*phi[iDim];
////          d_n     += d[iDim]*Normal[iDim];
////        }
////        phis1 += ProjVel * Psi[nVar-1];
////        
////        /*--- Convective flux at the wall node (adjoint density & energy only) ---*/
////        
////        /*--- Version 1 (full) ---*/
////        //Res_Conv_i[0] = ProjVel * Psi[0] - phis2 * ProjVel + phis1 * Gamma_Minus_One * sq_vel - ProjVel*Psi[0];
////        //Res_Conv_i[nVar-1] = ProjVel * Psi[nVar-1] + phis1 * Gamma_Minus_One - ProjVel*Psi[nVar-1];
////        
////        /*--- Simplified version ---*/
////        Res_Conv_i[0] = -(Psi[0] + phi_u + Psi[nVar-1]*Enthalpy)*ProjVel + d_n*Gamma_Minus_One*sq_vel;
////        
////        /*--- TO DO: Implicit contributions for convective part ---*/
////        
////        
////        /*--- Viscous flux contributions at the wall node ---*/
////        U = solver_container[TNE2_SOL]->node[iPoint]->GetSolution();
////        Laminar_Viscosity = solver_container[TNE2_SOL]->node[iPoint]->GetLaminarViscosity();
////        Eddy_Viscosity = solver_container[TNE2_SOL]->node[iPoint]->GetEddyViscosity(); // Should be zero at the wall
////        Density = U[0];
////        for (iDim = 0; iDim < nDim; iDim++) {
////          Velocity[iDim] = GridVel[iDim];
////        }
////        Energy = U[nDim+1] / Density;
////        SoundSpeed = sqrt(Gamma*Gamma_Minus_One*(Energy-sq_vel));
////        Pressure = (SoundSpeed * SoundSpeed * Density) / Gamma;
////        ViscDens = (Laminar_Viscosity + Eddy_Viscosity) / Density;
////        XiDens = Gamma * (Laminar_Viscosity/PRANDTL + Eddy_Viscosity/PRANDTL_TURB) / Density;
////        
////        /*--- Average of the derivatives of the adjoint variables ---*/
////        PsiVar_Grad = node[iPoint]->GetGradient();
////        
////        for (iDim = 0; iDim < nDim; iDim++) {
////          GradPsiE[iDim] =  PsiVar_Grad[nVar-1][iDim];
////          for (jDim = 0; jDim < nDim; jDim++)
////            GradPhi[iDim][jDim] =  PsiVar_Grad[iDim+1][jDim];
////        }
////        
////        /*--- Impose dPhiE_dn = 0 (adiabatic walls with frozen viscosity). Note
////         that this is where a different adjoint boundary condition for temperature
////         could be imposed. ---*/
////        dPhiE_dn = 0.0;
////
////        if (nDim ==2) {
////          
////          /*--- Compute the adjoint stress tensor ---*/
////          Sigma_xx  = ViscDens * (FOUR3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1]);
////          Sigma_yy  = ViscDens * (-TWO3 * GradPhi[0][0] + FOUR3 * GradPhi[1][1]);
////          Sigma_xy  = ViscDens * (GradPhi[1][0] + GradPhi[0][1]);
////          Sigma_xx5 = ViscDens * ( FOUR3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1]);
////          Sigma_yy5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] + FOUR3 * Velocity[1] * GradPsiE[1]);
////          Sigma_xy5 = ViscDens * (Velocity[0] * GradPsiE[1] + Velocity[1] * GradPsiE[0]);
////          Sigma_5   = XiDens * dPhiE_dn;
////          eta_xx    = Sigma_xx + Sigma_xx5;
////          eta_yy    = Sigma_yy + Sigma_yy5;
////          eta_xy    = Sigma_xy + Sigma_xy5;
////          
////          /*--- Viscous flux at the wall node (adjoint density & energy only) ---*/
////          Res_Visc_i[0] = - (Velocity[0] * Normal[0] * eta_xx  + Velocity[1] * Normal[1] * eta_yy
////                             + (Velocity[0] * Normal[1] + Velocity[1] * Normal[0]) * eta_xy
////                             - (sq_vel - Pressure/(Density*Gamma_Minus_One)) * Sigma_5);
////          Res_Visc_i[1] = 0.0;
////          Res_Visc_i[2] = 0.0;
////          
////        } else if (nDim == 3) {
////          
////          /*--- Compute the adjoint stress tensor ---*/
////          Sigma_xx  = ViscDens * (FOUR3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1] - TWO3  * GradPhi[2][2]);
////          Sigma_yy  = ViscDens * (-TWO3 * GradPhi[0][0] + FOUR3 * GradPhi[1][1] - TWO3  * GradPhi[2][2]);
////          Sigma_zz  = ViscDens * (-TWO3 * GradPhi[0][0] -  TWO3 * GradPhi[1][1] + FOUR3 * GradPhi[2][2]);
////          Sigma_xy  = ViscDens * (GradPhi[1][0] + GradPhi[0][1]);
////          Sigma_xz  = ViscDens * (GradPhi[2][0] + GradPhi[0][2]);
////          Sigma_yz  = ViscDens * (GradPhi[2][1] + GradPhi[1][2]);
////          Sigma_xx5 = ViscDens * ( FOUR3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1] -  TWO3 * Velocity[2] * GradPsiE[2]);
////          Sigma_yy5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] + FOUR3 * Velocity[1] * GradPsiE[1] -  TWO3 * Velocity[2] * GradPsiE[2]);
////          Sigma_zz5 = ViscDens * (- TWO3 * Velocity[0] * GradPsiE[0] -  TWO3 * Velocity[1] * GradPsiE[1] + FOUR3 * Velocity[2] * GradPsiE[2]);
////          Sigma_xy5 = ViscDens * (Velocity[0] * GradPsiE[1] + Velocity[1] * GradPsiE[0]);
////          Sigma_xz5 = ViscDens * (Velocity[0] * GradPsiE[2] + Velocity[2] * GradPsiE[0]);
////          Sigma_yz5 = ViscDens * (Velocity[1] * GradPsiE[2] + Velocity[2] * GradPsiE[1]);
////          Sigma_5   = XiDens * dPhiE_dn;
////          eta_xx    = Sigma_xx + Sigma_xx5; eta_yy = Sigma_yy + Sigma_yy5; eta_zz = Sigma_zz + Sigma_zz5;
////          eta_xy    = Sigma_xy + Sigma_xy5; eta_xz = Sigma_xz + Sigma_xz5; eta_yz = Sigma_yz + Sigma_yz5;
////          
////          /*--- Viscous flux at the wall node (adjoint density & energy only) ---*/
////          Res_Visc_i[0] = - (Velocity[0] * Normal[0] * eta_xx  + Velocity[1] * Normal[1] * eta_yy + Velocity[2] * Normal[2] * eta_zz
////                             + (Velocity[0] * Normal[1] + Velocity[1] * Normal[0]) * eta_xy
////                             + (Velocity[0] * Normal[2] + Velocity[2] * Normal[0]) * eta_xz
////                             + (Velocity[2] * Normal[1] + Velocity[1] * Normal[2]) * eta_yz
////                             - (sq_vel - Pressure/(Density*Gamma_Minus_One)) * Sigma_5);
////          Res_Visc_i[1] = 0.0;
////          Res_Visc_i[2] = 0.0;
////          Res_Visc_i[3] = 0.0;
////        }
////      }
//		}
//	}

