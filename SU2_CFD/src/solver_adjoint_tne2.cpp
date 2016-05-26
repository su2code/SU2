/*!
 * \file solution_adjoint_tne2.cpp
 * \brief Main subrotuines for solving adjoint problems (Euler, Navier-Stokes, etc.).
 * \author S. Copeland
 * \version 4.0.0 "Cardinal"
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
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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
  su2double dull_val;
  string text_line, mesh_filename;
  string filename, AdjExt;
	ifstream restart_file;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
	Residual     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max    = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
	Residual_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
	Res_Conv_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_i[iVar]    = 0.0;
  Res_Visc_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_i[iVar]    = 0.0;
	Res_Conv_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_j[iVar]    = 0.0;
  Res_Visc_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_j[iVar]    = 0.0;
  
	/*--- Define some auxiliary vectors related to the solution ---*/
	Solution   = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
	Solution_i = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar]   = 0.0;
  Solution_j = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar]   = 0.0;
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
	Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
	Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
	Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
	/*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
	if (config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_CENTERED) {
		iPoint_UndLapl = new su2double [nPoint];
		jPoint_UndLapl = new su2double [nPoint];
	}
  
	/*--- Define some auxiliary vectors related to the geometry ---*/
	Vector_i = new su2double[nDim];
  Vector_j = new su2double[nDim];
  
  /*--- Point to point Jacobians. These are always defined because
   they are also used for sensitivity calculations. ---*/
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
	/*--- Allocate Jacobians for implicit time-stepping ---*/
  if (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT) {
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
    
    if (rank == MASTER_NODE)
      cout << "Initialize Jacobian structure (Adjoint Euler). MG level: " << iMesh <<"." << endl;
		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
  
  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (Adjoint Euler). MG level: " << iMesh <<"." << endl;
  }
  
	/*--- Allocate arrays for gradient computation by least squares ---*/
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
  
	/*--- Sensitivity arrays on boundaries ---*/
	CSensitivity = new su2double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		CSensitivity[iMarker] = new su2double [geometry->nVertex[iMarker]];
	}
	Sens_Geo   = new su2double[nMarker];
	Sens_Mach  = new su2double[nMarker];
	Sens_AoA   = new su2double[nMarker];
	Sens_Press = new su2double[nMarker];
	Sens_Temp  = new su2double[nMarker];
  
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
  PsiRho_Inf = new su2double [nSpecies];
  Phi_Inf    = new su2double [nDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    PsiRho_Inf[iSpecies] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Phi_Inf[iDim] = 0.0;
  PsiE_Inf   = 0.0;
  PsiEve_Inf = 0.0;
  
	if (!restart || (iMesh != MESH_0)) {
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
      case DRAG_COEFFICIENT:        AdjExt = "_cd.dat";       break;
      case LIFT_COEFFICIENT:        AdjExt = "_cl.dat";       break;
      case SIDEFORCE_COEFFICIENT:   AdjExt = "_csf.dat";      break;
      case INVERSE_DESIGN_PRESSURE: AdjExt = "_invpress.dat"; break;
      case INVERSE_DESIGN_HEATFLUX: AdjExt = "_invheat.dat";  break;
      case MOMENT_X_COEFFICIENT:    AdjExt = "_cmx.dat";      break;
      case MOMENT_Y_COEFFICIENT:    AdjExt = "_cmy.dat";      break;
      case MOMENT_Z_COEFFICIENT:    AdjExt = "_cmz.dat";      break;
      case EFFICIENCY:              AdjExt = "_eff.dat";      break;
      case EQUIVALENT_AREA:         AdjExt = "_ea.dat";       break;
      case NEARFIELD_PRESSURE:      AdjExt = "_nfp.dat";      break;
      case FORCE_X_COEFFICIENT:     AdjExt = "_cfx.dat";      break;
      case FORCE_Y_COEFFICIENT:     AdjExt = "_cfy.dat";      break;
      case FORCE_Z_COEFFICIENT:     AdjExt = "_cfz.dat";      break;
      case THRUST_COEFFICIENT:      AdjExt = "_ct.dat";       break;
      case TORQUE_COEFFICIENT:      AdjExt = "_cq.dat";       break;
      case FIGURE_OF_MERIT:         AdjExt = "_merit.dat";    break;
      case FREE_SURFACE:            AdjExt = "_fs.dat";       break;
      case TOTAL_HEATFLUX:          AdjExt = "_totheat.dat";  break;
      case MAXIMUM_HEATFLUX:        AdjExt = "_maxheat.dat";  break;
		}
		filename.append(AdjExt);
		restart_file.open(filename.data(), ios::in);
    
		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
		  if (rank == MASTER_NODE)
		    cout << "There is no adjoint restart file!! " << filename.data() << "."<< endl;
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
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
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
		for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
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
	su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta;
  su2double phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  MPI_Status status;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
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
	su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta;
  su2double phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  MPI_Status status;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
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
	su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta;
  su2double phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  su2double *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
#ifdef HAVE_MPI
  MPI_Status status;
#endif

	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
      Buffer_Send_Limit = new su2double[nBufferS_Vector];
      
      /*--- Copy the limiter that should be sent ---*/
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
	su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta;
  su2double phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  su2double *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
#ifdef HAVE_MPI
  MPI_Status status;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar*nDim;   nBufferR_Vector = nVertexR*nVar*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
      Buffer_Send_Gradient = new su2double[nBufferS_Vector];
      
      /*--- Copy the gradient that should be sent ---*/
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
	su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta;
  su2double phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  su2double *Buffer_Receive_Undivided_Laplacian = NULL;
  su2double *Buffer_Send_Undivided_Laplacian = NULL;
  
#ifdef HAVE_MPI
  MPI_Status status;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Undivided_Laplacian = new su2double [nBufferR_Vector];
      Buffer_Send_Undivided_Laplacian = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Undivided_Laplacian[iVar*nVertexS+iVertex] = node[iPoint]->GetUndivided_Laplacian(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
	  SU2_MPI::Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                               Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
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
	su2double Alpha      = (config->GetAoA()*PI_NUMBER)/180.0;
	su2double Beta       = (config->GetAoS()*PI_NUMBER)/180.0;
	su2double RefAreaCoeff    = config->GetRefAreaCoeff();
	su2double RefLengthMoment  = config->GetRefLengthMoment();
	su2double *RefOriginMoment = config->GetRefOriginMoment(0);
  su2double *ForceProj_Vector, x = 0.0, y = 0.0, z = 0.0, *Normal, C_d, C_l, C_t, C_q;
	su2double x_origin, y_origin, z_origin, Area;
	su2double RefVel2, RefDensity;
  int rank = MASTER_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
	ForceProj_Vector = new su2double[nDim];
  
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
#ifndef HAVE_MPI
	C_d = solver_container[TNE2_SOL]->GetTotal_CDrag();
	C_l = solver_container[TNE2_SOL]->GetTotal_CLift();
	C_t = solver_container[TNE2_SOL]->GetTotal_CT();
	C_q = solver_container[TNE2_SOL]->GetTotal_CQ();
#else
	su2double *sbuf_force = new su2double[4];
	su2double *rbuf_force = new su2double[4];
	sbuf_force[0] = solver_container[TNE2_SOL]->GetTotal_CDrag();
	sbuf_force[1] = solver_container[TNE2_SOL]->GetTotal_CLift();
	sbuf_force[2] = solver_container[TNE2_SOL]->GetTotal_CT();
	sbuf_force[3] = solver_container[TNE2_SOL]->GetTotal_CQ();

	SU2_MPI::Reduce(sbuf_force, rbuf_force, 4, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
	SU2_MPI::Bcast(rbuf_force, 4, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

	C_d = rbuf_force[0];
	C_l = rbuf_force[1];
	C_t = rbuf_force[2];
	C_q = rbuf_force[3];
	delete [] sbuf_force;
	delete [] rbuf_force;
#endif
  
	/*--- Compute coefficients needed for objective function evaluation. ---*/
	su2double C_p    = 1.0/(0.5*RefDensity*RefAreaCoeff*RefVel2);
	su2double invCD  = 1.0 / C_d;
	su2double CLCD2  = C_l / (C_d*C_d);
  
	x_origin = RefOriginMoment[0];
  y_origin = RefOriginMoment[1];
  z_origin = RefOriginMoment[2];
  
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		if ((config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) &&
        (config->GetMarker_All_Monitoring(iMarker) == YES))
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
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
              exit(EXIT_FAILURE);
            }
            if (nDim == 3) { ForceProj_Vector[0] = -C_p*sin(Beta) * cos(Alpha); ForceProj_Vector[1] = C_p*cos(Beta); ForceProj_Vector[2] = -C_p*sin(Beta) * sin(Alpha); }
            break;
          case INVERSE_DESIGN_PRESSURE :
            if (nDim == 2) {
              Area = sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1]);
              ForceProj_Vector[0] = -C_p*Normal[0]/Area; ForceProj_Vector[1] = -C_p*Normal[1]/Area;
            }
            if (nDim == 3) {
              Area = sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);
              ForceProj_Vector[0] = -C_p*Normal[0]/Area; ForceProj_Vector[1] = -C_p*Normal[1]/Area; ForceProj_Vector[2] = -C_p*Normal[2]/Area;
            }
            break;
          case INVERSE_DESIGN_HEATFLUX:
            if (nDim == 2) { ForceProj_Vector[0] = 0.0;
              ForceProj_Vector[1] = 0.0; }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0;
              ForceProj_Vector[1] = 0.0;
              ForceProj_Vector[2] = 0.0; }
            break;
          case MOMENT_X_COEFFICIENT :
            if ((nDim == 2) && (rank == MASTER_NODE)) { cout << "This functional is not possible in 2D!!" << endl; exit(EXIT_FAILURE); }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = -C_p*(z - z_origin)/RefLengthMoment; ForceProj_Vector[2] = C_p*(y - y_origin)/RefLengthMoment; }
            break;
          case MOMENT_Y_COEFFICIENT :
            if ((nDim == 2) && (rank == MASTER_NODE)) { cout << "This functional is not possible in 2D!!" << endl; exit(EXIT_FAILURE); }
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
              exit(EXIT_FAILURE);
            }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = C_p; }
            break;
          case TOTAL_HEATFLUX:
            if (nDim == 2) { ForceProj_Vector[0] = 0.0;
              ForceProj_Vector[1] = 0.0; }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0;
              ForceProj_Vector[1] = 0.0;
              ForceProj_Vector[2] = 0.0; }
            break;
          case MAXIMUM_HEATFLUX:
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
	su2double Area_Children, Area_Parent;
  su2double *Solution, *Solution_Fine;
  
	bool restart = config->GetRestart();
  
  /*--- If restart solution, then interpolate the flow solution to
   all the multigrid levels, this is important with the dual time strategy ---*/
  if (restart) {
    Solution = new su2double[nVar];
    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
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

  bool implicit, second_order, center, center_jst, limiter, RightSol;
	unsigned long iPoint, ErrorCounter = 0;
  su2double SharpEdge_Distance;
  int rank;
  
#ifndef HAVE_MPI
	rank = MASTER_NODE;
#else
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the direct problem may use different methods). ---*/
  second_order = ((config->GetSpatialOrder() == SECOND_ORDER) || (config->GetSpatialOrder() == SECOND_ORDER_LIMITER));
  limiter      = (config->GetSpatialOrder() == SECOND_ORDER_LIMITER);
  center       = (config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_CENTERED);
  center_jst   = (config->GetKind_Centered_AdjTNE2() == JST);
  implicit     = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
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
  if ((second_order) && (iMesh == MESH_0)) {
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
#ifdef HAVE_MPI
  unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
  SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
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
  bool implicit, second_order;
  unsigned short iVar, jVar;
	unsigned long iEdge, iPoint, jPoint;
  
  /*--- Set booleans from configuration settings ---*/
	implicit        = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
	second_order = ((config->GetKind_Centered_AdjTNE2() == JST) && (iMesh == MESH_0));
  
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
    numerics->SetPrimitive(solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive(),
                           solver_container[TNE2_SOL]->node[jPoint]->GetPrimitive());    

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
    
		if (second_order) {
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
	
  bool implicit, second_order, limiter;
	unsigned short iDim, iVar;
    unsigned long iEdge, iPoint, jPoint;
  su2double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
  su2double *Limiter_i, *Limiter_j;
  su2double *Psi_i, *Psi_j;
  su2double *U_i, *U_j, *V_i, *V_j;
  
  /*--- Initialize ---*/
  Limiter_i = NULL;
  Limiter_j = NULL;
  Psi_i     = NULL;
  Psi_j     = NULL;
  
  /*--- Set booleans based on config settings ---*/
	implicit        = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
	second_order    = (((config->GetSpatialOrder_AdjTNE2() == SECOND_ORDER) || (config->GetSpatialOrder_AdjTNE2() == SECOND_ORDER_LIMITER)) && (iMesh == MESH_0));
  limiter         = (config->GetSpatialOrder_AdjTNE2() == SECOND_ORDER_LIMITER);

  if (second_order)
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
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Retrieve node numbers and pass edge normal to CNumerics ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Pass conserved and primitive variables from CVariable to CNumerics class ---*/
    U_i = solver_container[TNE2_SOL]->node[iPoint]->GetSolution();
    U_j = solver_container[TNE2_SOL]->node[jPoint]->GetSolution();
    V_i = solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive();
    V_j = solver_container[TNE2_SOL]->node[jPoint]->GetPrimitive();
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
    if (second_order) {
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
    numerics->SetPrimitive(solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive(),
                           solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive());
    
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
	su2double *Diff;
  
  Diff = new su2double[nVar];
  
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
  su2double *Psi_mirror;
  
  Psi_mirror = new su2double[nVar];
  
	/*--- Loop over all boundaries and include an extra contribution
	 from a halo node. Find the nearest normal, interior point
	 for a boundary node and make a linear approximation. ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
				config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
				config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
				config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY) {
      
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
	su2double *local_Residual, Vol, Delta, Res;
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
      AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
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
	unsigned long iPoint, total_index, IterLinSol = 0;
	su2double Delta, Vol;
  
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
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      
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
  
  /*--- Solve or smooth the linear system ---*/
  
  CSysSolve system;
  IterLinSol = system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
	/*--- Update solution (system written in terms of increments) ---*/
  
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
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
	su2double *d, *Normal, *UnitNormal;
  su2double *Psi, *U, *V, *dPdU, *USens;
  su2double rho, rhou, rhov, rhow, rhoE, rhoEve, H, p;
  su2double conspsi, Area;
  su2double Mach_Inf;
  su2double **PrimVar_Grad, **ConsVar_Grad, *ConsPsi_Grad;
  su2double ConsPsi, d_press, grad_v, v_gradconspsi;
  
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
  UnitNormal = new su2double[nDim];
  USens      = new su2double[nVar];
  
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
    if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {
      
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
    
    if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {
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
    
    if (config->GetMarker_All_KindBC(iMarker) == FAR_FIELD) {
      
      Sens_Mach[iMarker]  = 0.0;
      Sens_AoA[iMarker]   = 0.0;
      Sens_Press[iMarker] = 0.0;
      Sens_Temp[iMarker]  = 0.0;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          Psi      = node[iPoint]->GetSolution();
          U        = solver_container[TNE2_SOL]->node[iPoint]->GetSolution();
          V        = solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive();
          dPdU     = solver_container[TNE2_SOL]->node[iPoint]->GetdPdU();
          Normal   = geometry->vertex[iMarker][iVertex]->GetNormal();
          Mach_Inf = config->GetMach();
          
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
//    if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {
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
//          Mach_Inf   = config->GetMach();
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
	su2double *d, *Normal, Area, *UnitNormal, *Coord;
  su2double *Psi, *Psi_Aux, phin, bcn;
  su2double *U, *V, *dPdU;
  
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
	UnitNormal    = new su2double[nDim];
	Psi           = new su2double[nVar];
  
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
      V    = solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive();
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
	su2double *Normal, Area, *UnitNormal, *Coord;
  su2double *Psi, *Psi_Aux, phin;
  su2double *U, *V, *dPdU;
  
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
	UnitNormal    = new su2double[nDim];
	Psi           = new su2double[nVar];
  
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
      V = solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive();
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
	su2double *Normal, *U_domain, *U_infty, *V_domain, *V_infty;
  su2double *Psi_domain, *Psi_infty;
  
  /*--- Set booleans from config settings ---*/
	implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  /*--- Allocate arrays ---*/
	Normal     = new su2double[nDim];
	Psi_domain = new su2double[nVar];
  Psi_infty  = new su2double[nVar];
  
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
      V_domain = solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive();
      V_infty  = solver_container[TNE2_SOL]->node_infty->GetPrimitive();
    
      /*--- Pass conserved & primitive variables to CNumerics ---*/
			conv_numerics->SetConservative(U_domain, U_infty);
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
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
  su2double dull_val;
  string text_line, mesh_filename;
  string filename, AdjExt;
	ifstream restart_file;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
	Res_Conv_i = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Res_Conv_i[iVar]    = 0.0;
  Res_Visc_i   = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Res_Visc_i[iVar]    = 0.0;
	Res_Conv_j = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Res_Conv_j[iVar]    = 0.0;
  Res_Visc_j   = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Res_Visc_j[iVar]    = 0.0;
  
	/*--- Define some auxiliary arrays related to the solution ---*/
	Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
	Solution_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
	Solution_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;
  
	/*--- Define some auxiliary vectors related to the geometry ---*/
	Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
	Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
	Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Point to point Jacobians. These are always defined because
   they are also used for sensitivity calculations. ---*/
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  /*--- Solution and residual vectors ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT) {
    
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
    if (rank == MASTER_NODE)
      cout << "Initialize Jacobian structure (Adjoint N-S). MG level: " << iMesh <<"." << endl;
		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (Adjoint N-S). MG level: " << iMesh <<"." << endl;
  }
  
	/*--- Array structures for computation of gradients by least squares ---*/
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
  
	/*--- Sensitivity definition and coefficient on all markers ---*/
	CSensitivity = new su2double* [nMarker];
	for (iMarker=0; iMarker<nMarker; iMarker++) {
		CSensitivity[iMarker] = new su2double [geometry->nVertex[iMarker]];
	}
	Sens_Geo   = new su2double[nMarker];
	Sens_Mach  = new su2double[nMarker];
	Sens_AoA   = new su2double[nMarker];
	Sens_Press = new su2double[nMarker];
	Sens_Temp  = new su2double[nMarker];
  
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
	PsiRho_Inf = new su2double [nSpecies];
  Phi_Inf    = new su2double [nDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    PsiRho_Inf[iSpecies] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Phi_Inf[iDim] = 0.0;
  if ((config->GetKind_ObjFunc() == TOTAL_HEATFLUX) ||
      (config->GetKind_ObjFunc() == MAXIMUM_HEATFLUX) ||
      (config->GetKind_ObjFunc() == INVERSE_DESIGN_HEATFLUX)) {
    PsiE_Inf = -1.0;
    PsiEve_Inf = -1.0;
  } else {
    PsiE_Inf = 0.0;
    PsiEve_Inf = 0.0;
  }
  
	if (!restart || (iMesh != MESH_0)) {
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
      case DRAG_COEFFICIENT:        AdjExt = "_cd.dat";       break;
      case LIFT_COEFFICIENT:        AdjExt = "_cl.dat";       break;
      case SIDEFORCE_COEFFICIENT:   AdjExt = "_csf.dat";      break;
      case INVERSE_DESIGN_PRESSURE: AdjExt = "_invpress.dat"; break;
      case INVERSE_DESIGN_HEATFLUX: AdjExt = "_invheat.dat";  break;
      case MOMENT_X_COEFFICIENT:    AdjExt = "_cmx.dat";      break;
      case MOMENT_Y_COEFFICIENT:    AdjExt = "_cmy.dat";      break;
      case MOMENT_Z_COEFFICIENT:    AdjExt = "_cmz.dat";      break;
      case EFFICIENCY:              AdjExt = "_eff.dat";      break;
      case EQUIVALENT_AREA:         AdjExt = "_ea.dat";       break;
      case NEARFIELD_PRESSURE:      AdjExt = "_nfp.dat";      break;
      case FORCE_X_COEFFICIENT:     AdjExt = "_cfx.dat";      break;
      case FORCE_Y_COEFFICIENT:     AdjExt = "_cfy.dat";      break;
      case FORCE_Z_COEFFICIENT:     AdjExt = "_cfz.dat";      break;
      case THRUST_COEFFICIENT:      AdjExt = "_ct.dat";       break;
      case TORQUE_COEFFICIENT:      AdjExt = "_cq.dat";       break;
      case FIGURE_OF_MERIT:         AdjExt = "_merit.dat";    break;
      case FREE_SURFACE:            AdjExt = "_fs.dat";       break;
      case TOTAL_HEATFLUX:          AdjExt = "_totheat.dat";  break;
      case MAXIMUM_HEATFLUX:        AdjExt = "_maxheat.dat";  break;
		}
		filename.append(AdjExt);
		restart_file.open(filename.data(), ios::in);
    
		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
		  if (rank == MASTER_NODE)
		    cout << "There is no adjoint restart file!! " << filename.data() << "."<< endl;
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
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
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
		for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
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
  bool implicit, second_order, center, center_jst, limiter, RightSol;
	unsigned long iPoint, ErrorCounter;
  su2double SharpEdge_Distance;
  int rank;
  
#ifndef HAVE_MPI
	rank = MASTER_NODE;
#else
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the flow problem may use different methods). ---*/
  implicit     = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  second_order = ((config->GetSpatialOrder() == SECOND_ORDER) || (config->GetSpatialOrder() == SECOND_ORDER_LIMITER));
  limiter      = (config->GetSpatialOrder() == SECOND_ORDER_LIMITER);
  center       = (config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_CENTERED);
  center_jst   = (config->GetKind_Centered_AdjTNE2() == JST);
  
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
  if ((second_order) && (iMesh == MESH_0)) {
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
#ifdef HAVE_MPI
  unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
  SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
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
    numerics->SetPrimitive(solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive(),
                           solver_container[TNE2_SOL]->node[jPoint]->GetPrimitive());
    
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
    numerics->SetPrimitive(solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive(),
                           solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive());
    second_numerics->SetConservative(solver_container[TNE2_SOL]->node[iPoint]->GetSolution(),
                                     solver_container[TNE2_SOL]->node[iPoint]->GetSolution());
    second_numerics->SetPrimitive(solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive(),
                                  solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive());
    
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
    
    /*--- Add and subtract to the residual ---*/
		LinSysRes.AddBlock(iPoint, Residual_i);
	}
}

void CAdjTNE2NSSolver::Viscous_Sensitivity(CGeometry *geometry,
                                           CSolver **solver_container,
                                           CNumerics *numerics,
                                           CConfig *config) {
  // Boolean declarations
  
  // Counter & iterator declarations
	unsigned long iVertex, iPoint;
	unsigned short iDim, jDim, iMarker, iVar;
  
  //Geometry declarations
  su2double Area, eps;
  su2double *Normal, *UnitNormal;
  
  Normal     = NULL;
  UnitNormal = NULL;
  
  UnitNormal = new su2double[nDim];
  
  // Direct problem declarations
  unsigned short T_INDEX, TVE_INDEX, VEL_INDEX;
  su2double rho, H;
  su2double div_vel, dnT, dnTve;
  su2double mu, ktr, kve;
  su2double *U, *dnvel;
  su2double **GradV, **sigma;
  
  U     = NULL;
  GradV = NULL;
  
  dnvel = new su2double [nDim];
  sigma = new su2double *[nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    sigma[iDim] = new su2double[nDim];
    
  
  // Adjoint problem declarations
  su2double vartheta, dnPsiE, dnPsiEve, div_phi, GPsiEdotVel;
  su2double B1, B21, B22, B23, B24, B31, B32, B33, B34;
  su2double *Psi, *d;
  su2double **GradPsi, **SigmaPhi, **SigmaPsiE;
  
  Psi     = NULL;
  d       = NULL;
  GradPsi = NULL;
  
  SigmaPhi = new su2double *[nDim];
  SigmaPsiE = new su2double *[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    SigmaPhi[iDim] = new su2double[nDim];
    SigmaPsiE[iDim] = new su2double[nDim];
  }
  
  /*--- Initialize total sensitivites ---*/
  Total_Sens_Geo   = 0.0;
  Total_Sens_Mach  = 0.0;
  Total_Sens_AoA   = 0.0;
  Total_Sens_Press = 0.0;
  Total_Sens_Temp  = 0.0;
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    /*--- Initialize marker geometric sensitivity ---*/
    Sens_Geo[iMarker] = 0.0;
    
    
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC    ) ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC   ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC)) {
      
      /*--- Loop over all boundary nodes ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Calculate geometric quantities ---*/
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++)
            UnitNormal[iDim] = Normal[iDim] / Area;
          
          /*--- Load flow quantities ---*/
          U     = solver_container[TNE2_SOL]->node[iPoint]->GetSolution();
          rho   = solver_container[TNE2_SOL]->node[iPoint]->GetDensity();
          H     = solver_container[TNE2_SOL]->node[iPoint]->GetEnthalpy();
          mu    = solver_container[TNE2_SOL]->node[iPoint]->GetLaminarViscosity();
          ktr   = solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity();
          kve   = solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve();
          GradV = solver_container[TNE2_SOL]->node[iPoint]->GetGradient_Primitive();
          
          VEL_INDEX = solver_container[TNE2_SOL]->node[iPoint]->GetVelIndex();
          T_INDEX   = solver_container[TNE2_SOL]->node[iPoint]->GetTIndex();
          TVE_INDEX = solver_container[TNE2_SOL]->node[iPoint]->GetTveIndex();
          
          /*--- Load adjoint quantities ---*/
          Psi     = node[iPoint]->GetSolution();
          GradPsi = node[iPoint]->GetGradient();
          
          /*--- Calculate support quantities ---*/
          
          // vartheta
          vartheta = 0.0;
          for (iVar = 0; iVar < nSpecies+nDim; iVar++)
            vartheta += U[iVar]*Psi[iVar];
          vartheta += rho*H*Psi[nSpecies+nDim];
          vartheta += U[nSpecies+nDim+1]*Psi[nSpecies+nDim+1];
          
          // Velocity normal derivatives
          for (iDim = 0; iDim < nDim; iDim++) {
            dnvel[iDim] = 0.0;
            for (jDim = 0; jDim < nDim; jDim++)
              dnvel[iDim] += GradV[VEL_INDEX+iDim][jDim]*UnitNormal[jDim];
          }
          
          // Temperature normal derivatives
          dnT   = 0.0;
          dnTve = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            dnT   += GradV[T_INDEX][iDim]   * UnitNormal[iDim];
            dnTve += GradV[TVE_INDEX][iDim] * UnitNormal[iDim];
          }
          
          // Adjoint energy normal derivatives
          dnPsiE   = 0.0;
          dnPsiEve = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            dnPsiE   += GradPsi[nSpecies+nDim][iDim]   * UnitNormal[iDim];
            dnPsiEve += GradPsi[nSpecies+nDim+1][iDim] * UnitNormal[iDim];
          }
          
          // Viscous stress tensor
          div_vel = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            div_vel += GradV[VEL_INDEX+iDim][iDim];
            for (jDim =0 ; jDim < nDim; jDim++)
              sigma[iDim][jDim] = 0.0;
          }
          for (iDim = 0; iDim < nDim; iDim++) {
            for (jDim = 0; jDim < nDim; jDim++) {
              sigma[iDim][jDim] += mu*(GradV[VEL_INDEX+iDim][jDim] +
                                      GradV[VEL_INDEX+jDim][iDim]  );
            }
            sigma[iDim][iDim] -= 2.0/3.0 * div_vel;
          }
          
          // SigmaPhi
          div_phi = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            div_phi += GradPsi[nSpecies+iDim][iDim];
            for (jDim =0 ; jDim < nDim; jDim++)
              SigmaPhi[iDim][jDim] = 0.0;
          }
          for (iDim = 0; iDim < nDim; iDim++) {
            for (jDim = 0; jDim < nDim; jDim++) {
              SigmaPhi[iDim][jDim] += GradPsi[nSpecies+iDim][jDim] +
                                      GradPsi[nSpecies+jDim][iDim];
            }
            SigmaPhi[iDim][iDim] -= 2.0/3.0 * div_phi;
          }
          
          // SigmaPsi5
          GPsiEdotVel = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            GPsiEdotVel += GradPsi[nSpecies+nDim][iDim] *
                           U[nSpecies+iDim]/rho;
            for (jDim = 0; jDim < nDim; jDim++)
              SigmaPsiE[iDim][jDim] = 0.0;
          }
          for (iDim = 0; iDim < nDim; iDim++) {
            for (jDim = 0; jDim < nDim; jDim++) {
              SigmaPsiE[iDim][jDim] += GradPsi[nSpecies+nDim][iDim] *
                                       U[nSpecies+jDim]/rho +
                                       GradPsi[nSpecies+nDim][jDim] *
                                       U[nSpecies+iDim]/rho;
            }
            SigmaPsiE[iDim][iDim] -= 2.0/3.0*GPsiEdotVel;
          }
          
          
          /*--- B1: Convective sensitivity ---*/
          // Note: The deltaP term is always canceled using adjoint B.C.'s
          // and is not included in the surface sensitivities.
          B1 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            B1 = -vartheta * dnvel[iDim] * UnitNormal[iDim];
          
          /*--- B21: 1st order viscous sensitivity, diffusion ---*/
          B21 = 0.0;
          
          /*--- B22: 1st-viscous sensitivity, viscosity ---*/
          // Note: The deltaSigma term is always canceled using adjoint B.C.'s
          // and is not included in the surface sensitivities.
          B22 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            for (jDim = 0; jDim < nDim; jDim++)
              B22 -= Psi[nSpecies+nDim]*sigma[iDim][jDim]*dnvel[jDim]*UnitNormal[iDim];
          
          /*--- B23: 1st-viscous sensitivity, tr heat flux ---*/
          // Note: This term is used in adjoint B.C. for isothermal cases
          B23 = 0.0;
          
          /*--- B24: 1st-viscous sensitivity, ve heat flux ---*/
          // Note: This term is used in adjoint B.C. for isothermal cases
          B24 = 0.0;
          
          /*--- B31: 2nd-viscous sensitivity, diffusion ---*/
          // Note: This term is used in adjoint B.C. for non-catalytic cases
          B31 = 0.0;
          
          /*--- B32: 2nd-viscous sensitivity, viscosity ---*/
          B32 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            for (jDim = 0; jDim < nDim; jDim++) {
              B32 -= mu*(SigmaPhi[iDim][jDim] +
                         SigmaPsiE[iDim][jDim] ) * dnvel[jDim]*UnitNormal[iDim];
            }
          }
          
          /*--- B33: 2nd-viscous sensitivity, tr heat flux ---*/
          B33 = 0.0;
          if (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)
            B33 = -ktr * dnPsiE * dnT;
          
          /*--- B34: 2nd-viscous sensitivity, ve heat flux ---*/
          B34 = 0.0;
          if (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)
            B34 = -kve * (dnPsiE+dnPsiEve) * dnTve;
          
          /*--- Gather all sensitivities for dJ/dS ---*/
          CSensitivity[iMarker][iVertex] = (B1 - (B21+B22+B23+B24) +
                                            (B31+B32+B33+B34)       )*Area;
          
          /*--- If sharp edge, set the sensitivity to 0 on that region ---*/
          
          if (config->GetSens_Remove_Sharp()) {
            eps = config->GetLimiterCoeff()*config->GetRefElemLength();
            if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetSharpEdgesCoeff()*eps )
              CSensitivity[iMarker][iVertex] = 0.0;
          }
          
          /*--- Accumulate for total geometric sensitivity on the marker ---*/
          Sens_Geo[iMarker] -= CSensitivity[iMarker][iVertex];
          
        }
      }
      
      /*--- Accumulate all marker geometric sensitivities ---*/
      Total_Sens_Geo += Sens_Geo[iMarker];
    }
  }
  
  /*--- Farfield Sensitivity (Mach, AoA, Press, Temp), only for compressible flows ---*/
  
//  if (compressible) {
//    
//    for (iMarker = 0; iMarker < nMarker; iMarker++) {
//      
//      if (config->GetMarker_All_KindBC(iMarker) == FAR_FIELD) {
//        
//        Sens_Mach[iMarker]  = 0.0;
//        Sens_AoA[iMarker]   = 0.0;
//        Sens_Press[iMarker] = 0.0;
//        Sens_Temp[iMarker]  = 0.0;
//        
//        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
//          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
//          
//          if (geometry->node[iPoint]->GetDomain()) {
//            Psi = node[iPoint]->GetSolution();
//            U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
//            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
//            
//            Mach_Inf   = config->GetMach();
//            if (grid_movement) Mach_Inf = config->GetMach_Motion();
//            
//            r = U[0]; ru = U[1]; rv = U[2];
//            if (nDim == 2) { rw = 0.0; rE = U[3]; }
//            else { rw = U[3]; rE = U[4]; }
//            p = Gamma_Minus_One*(rE-(ru*ru + rv*rv + rw*rw)/(2*r));
//            
//            Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
//            Area = sqrt(Area);
//            for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;
//            
//            H = (rE + p)/r;
//            
//            dp_dr = Gamma_Minus_One*(ru*ru + rv*rv + rw*rw)/(2*r*r);
//            dp_dru = -Gamma_Minus_One*ru/r;
//            dp_drv = -Gamma_Minus_One*rv/r;
//            if (nDim == 2) { dp_drw = 0.0; dp_drE = Gamma_Minus_One; }
//            else { dp_drw = -Gamma_Minus_One*rw/r; dp_drE = Gamma_Minus_One; }
//            
//            dH_dr = (-H + dp_dr)/r; dH_dru = dp_dru/r; dH_drv = dp_drv/r;
//            if (nDim == 2) { dH_drw = 0.0; dH_drE = (1 + dp_drE)/r; }
//            else { dH_drw = dp_drw/r; dH_drE = (1 + dp_drE)/r; }
//            
//            if (nDim == 2) {
//              Jacobian_j[0][0] = 0.0;
//              Jacobian_j[1][0] = Area*UnitNormal[0];
//              Jacobian_j[2][0] = Area*UnitNormal[1];
//              Jacobian_j[3][0] = 0.0;
//              
//              Jacobian_j[0][1] = (-(ru*ru)/(r*r) + dp_dr)*Area*UnitNormal[0] + (-(ru*rv)/(r*r))*Area*UnitNormal[1];
//              Jacobian_j[1][1] = (2*ru/r + dp_dru)*Area*UnitNormal[0] + (rv/r)*Area*UnitNormal[1];
//              Jacobian_j[2][1] = (dp_drv)*Area*UnitNormal[0] + (ru/r)*Area*UnitNormal[1];
//              Jacobian_j[3][1] = (dp_drE)*Area*UnitNormal[0];
//              
//              Jacobian_j[0][2] = (-(ru*rv)/(r*r))*Area*UnitNormal[0] + (-(rv*rv)/(r*r) + dp_dr)*Area*UnitNormal[1];
//              Jacobian_j[1][2] = (rv/r)*Area*UnitNormal[0] + (dp_dru)*Area*UnitNormal[1];
//              Jacobian_j[2][2] = (ru/r)*Area*UnitNormal[0] + (2*rv/r + dp_drv)*Area*UnitNormal[1];
//              Jacobian_j[3][2] = (dp_drE)*Area*UnitNormal[1];
//              
//              Jacobian_j[0][3] = (ru*dH_dr)*Area*UnitNormal[0] + (rv*dH_dr)*Area*UnitNormal[1];
//              Jacobian_j[1][3] = (H + ru*dH_dru)*Area*UnitNormal[0] + (rv*dH_dru)*Area*UnitNormal[1];
//              Jacobian_j[2][3] = (ru*dH_drv)*Area*UnitNormal[0] + (H + rv*dH_drv)*Area*UnitNormal[1];
//              Jacobian_j[3][3] = (ru*dH_drE)*Area*UnitNormal[0] + (rv*dH_drE)*Area*UnitNormal[1];
//            }
//            else {
//              Jacobian_j[0][0] = 0.0;
//              Jacobian_j[1][0] = Area*UnitNormal[0];
//              Jacobian_j[2][0] = Area*UnitNormal[1];
//              Jacobian_j[3][0] = Area*UnitNormal[2];
//              Jacobian_j[4][0] = 0.0;
//              
//              Jacobian_j[0][1] = (-(ru*ru)/(r*r) + dp_dr)*Area*UnitNormal[0] + (-(ru*rv)/(r*r))*Area*UnitNormal[1] + (-(ru*rw)/(r*r))*Area*UnitNormal[2];
//              Jacobian_j[1][1] = (2*ru/r + dp_dru)*Area*UnitNormal[0] + (rv/r)*Area*UnitNormal[1] + (rw/r)*Area*UnitNormal[2];
//              Jacobian_j[2][1] = (dp_drv)*Area*UnitNormal[0] + (ru/r)*Area*UnitNormal[1];
//              Jacobian_j[3][1] = (dp_drw)*Area*UnitNormal[0] + (ru/r)*Area*UnitNormal[2];
//              Jacobian_j[4][1] = (dp_drE)*Area*UnitNormal[0];
//              
//              Jacobian_j[0][2] = (-(ru*rv)/(r*r))*Area*UnitNormal[0] + (-(rv*rv)/(r*r) + dp_dr)*Area*UnitNormal[1] + (-(rv*rw)/(r*r))*Area*UnitNormal[2];
//              Jacobian_j[1][2] = (rv/r)*Area*UnitNormal[0] + (dp_dru)*Area*UnitNormal[1];
//              Jacobian_j[2][2] = (ru/r)*Area*UnitNormal[0] + (2*rv/r + dp_drv)*Area*UnitNormal[1] + (rw/r)*Area*UnitNormal[2];
//              Jacobian_j[3][2] = (dp_drw)*Area*UnitNormal[1] + (rv/r)*Area*UnitNormal[2];
//              Jacobian_j[4][2] = (dp_drE)*Area*UnitNormal[1];
//              
//              Jacobian_j[0][3] = (-(ru*rw)/(r*r))*Area*UnitNormal[0] + (-(rv*rw)/(r*r))*Area*UnitNormal[1] + (-(rw*rw)/(r*r) + dp_dr)*Area*UnitNormal[2];
//              Jacobian_j[1][3] = (rw/r)*Area*UnitNormal[0] + (dp_dru)*Area*UnitNormal[2];
//              Jacobian_j[2][3] = (rw/r)*Area*UnitNormal[1] + (dp_drv)*Area*UnitNormal[2];
//              Jacobian_j[3][3] = (ru/r)*Area*UnitNormal[0] + (rv/r)*Area*UnitNormal[1] + (2*rw/r + dp_drw)*Area*UnitNormal[2];
//              Jacobian_j[4][3] = (dp_drE)*Area*UnitNormal[2];
//              
//              Jacobian_j[0][4] = (ru*dH_dr)*Area*UnitNormal[0] + (rv*dH_dr)*Area*UnitNormal[1] + (rw*dH_dr)*Area*UnitNormal[2];
//              Jacobian_j[1][4] = (H + ru*dH_dru)*Area*UnitNormal[0] + (rv*dH_dru)*Area*UnitNormal[1] + (rw*dH_dru)*Area*UnitNormal[2];
//              Jacobian_j[2][4] = (ru*dH_drv)*Area*UnitNormal[0] + (H + rv*dH_drv)*Area*UnitNormal[1] + (rw*dH_drv)*Area*UnitNormal[2];
//              Jacobian_j[3][4] = (ru*dH_drw)*Area*UnitNormal[0] + (rv*dH_drw)*Area*UnitNormal[1] + (H + rw*dH_drw)*Area*UnitNormal[2];
//              Jacobian_j[4][4] = (ru*dH_drE)*Area*UnitNormal[0] + (rv*dH_drE)*Area*UnitNormal[1] + (rw*dH_drE)*Area*UnitNormal[2];
//            }
//            
//            /*--- Mach number sensitivity ---*/
//            
//            USens[0] = 0.0; USens[1] = ru/Mach_Inf; USens[2] = rv/Mach_Inf;
//            if (nDim == 2) { USens[3] = Gamma*Mach_Inf*p; }
//            else { USens[3] = rw/Mach_Inf; USens[4] = Gamma*Mach_Inf*p; }
//            for (iPos = 0; iPos < nVar; iPos++) {
//              for (jPos = 0; jPos < nVar; jPos++) {
//                Sens_Mach[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
//              }
//            }
//            
//            /*--- AoA sensitivity ---*/
//            
//            USens[0] = 0.0;
//            if (nDim == 2) { USens[1] = -rv; USens[2] = ru; USens[3] = 0.0; }
//            else { USens[1] = -rw; USens[2] = 0.0; USens[3] = ru; USens[4] = 0.0; }
//            for (iPos = 0; iPos < nVar; iPos++) {
//              for (jPos = 0; jPos < nVar; jPos++) {
//                Sens_AoA[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
//              }
//            }
//            
//            /*--- Pressure sensitivity ---*/
//            
//            USens[0] = r/p; USens[1] = ru/p; USens[2] = rv/p;
//            if (nDim == 2) { USens[3] = rE/p; }
//            else { USens[3] = rw/p; USens[4] = rE/p; }
//            for (iPos = 0; iPos < nVar; iPos++) {
//              for (jPos = 0; jPos < nVar; jPos++) {
//                Sens_Press[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
//              }
//            }
//            
//            /*--- Temperature sensitivity ---*/
//            
//            T = p/(r*Gas_Constant);
//            USens[0] = -r/T; USens[1] = 0.5*ru/T; USens[2] = 0.5*rv/T;
//            if (nDim == 2) { USens[3] = (ru*ru + rv*rv + rw*rw)/(r*T); }
//            else { USens[3] = 0.5*rw/T; USens[4] = (ru*ru + rv*rv + rw*rw)/(r*T); }
//            for (iPos = 0; iPos < nVar; iPos++) {
//              for (jPos = 0; jPos < nVar; jPos++) {
//                Sens_Temp[iMarker] += Psi[iPos]*Jacobian_j[jPos][iPos]*USens[jPos];
//              }
//            }
//          }
//        }
//        Total_Sens_Mach -= Sens_Mach[iMarker];
//        Total_Sens_AoA -= Sens_AoA[iMarker];
//        Total_Sens_Press -= Sens_Press[iMarker];
//        Total_Sens_Temp -= Sens_Temp[iMarker];
//      }
//    }
//    
//    /*--- Explicit contribution from objective function quantity ---*/
//    
//    for (iMarker = 0; iMarker < nMarker; iMarker++) {
//      
//      if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {
//        
//        Sens_Mach[iMarker]  = 0.0;
//        Sens_AoA[iMarker]   = 0.0;
//        Sens_Press[iMarker] = 0.0;
//        Sens_Temp[iMarker]  = 0.0;
//        
//        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
//          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
//          
//          if (geometry->node[iPoint]->GetDomain()) {
//            
//            U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
//            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
//            p = solver_container[FLOW_SOL]->node[iPoint]->GetPressure();
//            
//            Mach_Inf   = config->GetMach();
//            if (grid_movement) Mach_Inf = config->GetMach_Motion();
//            
//            d = node[iPoint]->GetForceProj_Vector();
//            Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
//            Area = sqrt(Area);
//            for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;
//            
//            /*--- Mach number sensitivity ---*/
//            
//            for (iPos = 0; iPos < nDim; iPos++) Dd[iPos] = -(2/Mach_Inf)*d[iPos];
//            for (iPos = 0; iPos < nDim; iPos++) Sens_Mach[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];
//            
//            /*--- AoA sensitivity ---*/
//            /* Coefficients with an explicit AoA dependence - NOTE: Still need to implement right dependency for EFFICIENCY */
//            if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT || config->GetKind_ObjFunc() == LIFT_COEFFICIENT || config->GetKind_ObjFunc() == SIDEFORCE_COEFFICIENT || config->GetKind_ObjFunc() == EQUIVALENT_AREA || config->GetKind_ObjFunc() == NEARFIELD_PRESSURE) {
//            	if (nDim == 2) {
//            		D[0][0] = 0.0; D[0][1] = -1.0;
//            		D[1][0] = 1.0; D[1][1] = 0.0;
//            	}
//            	else {
//            		D[0][0] = 0.0; D[0][1] = 0.0; D[0][2] = -1.0;
//            		D[1][0] = 0.0; D[1][1] = 0.0; D[1][2] = 0.0;
//            		D[2][0] = 1.0; D[2][1] = 0.0; D[2][2] = 0.0;
//            	}
//            	for (iPos = 0; iPos < nDim; iPos++) Dd[iPos] = 0.0;
//            	for (iPos = 0; iPos < nDim; iPos++) {
//            		for (jPos = 0; jPos < nDim; jPos++)
//                  Dd[iPos] += D[iPos][jPos]*d[jPos];
//            	}
//            }
//            /* Coefficients with no explicit AoA dependece */
//            else {
//            	for (iPos = 0; iPos<nDim; iPos++) Dd[iPos] = 0.0;
//            }
//            
//            for (iPos = 0; iPos < nDim; iPos++)
//              Sens_AoA[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];
//            
//            /*--- Pressure sensitivity ---*/
//            
//            for (iPos = 0; iPos<nDim; iPos++) Dd[iPos] = -(1/p)*d[iPos];
//            for (iPos = 0; iPos<nDim; iPos++)
//              Sens_Press[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];
//            
//            /*--- Temperature sensitivity ---*/
//            
//            for (iPos = 0; iPos<nDim; iPos++) Dd[iPos] = 0.0;
//            for (iPos = 0; iPos<nDim; iPos++)
//              Sens_Temp[iMarker] += p*Dd[iPos]*Area*UnitNormal[iPos];
//            
//          }
//        }
//        
//        Total_Sens_Mach   += Sens_Mach[iMarker];
//        Total_Sens_AoA    += Sens_AoA[iMarker];
//        Total_Sens_Press  += Sens_Press[iMarker];
//        Total_Sens_Temp   += Sens_Temp[iMarker];
//        
//      }
//    }
//  }
  
#ifdef HAVE_MPI
  
  su2double MyTotal_Sens_Geo   = Total_Sens_Geo;     Total_Sens_Geo = 0.0;
//  su2double MyTotal_Sens_Mach  = Total_Sens_Mach;    Total_Sens_Mach = 0.0;
//  su2double MyTotal_Sens_AoA   = Total_Sens_AoA;     Total_Sens_AoA = 0.0;
//  su2double MyTotal_Sens_Press = Total_Sens_Press;   Total_Sens_Press = 0.0;
//  su2double MyTotal_Sens_Temp  = Total_Sens_Temp;    Total_Sens_Temp = 0.0;
  
  SU2_MPI::Allreduce(&MyTotal_Sens_Geo, &Total_Sens_Geo, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyTotal_Sens_Mach, &Total_Sens_Mach, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyTotal_Sens_AoA, &Total_Sens_AoA, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyTotal_Sens_Press, &Total_Sens_Press, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyTotal_Sens_Temp, &Total_Sens_Temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#endif
  
  /*--- Deallocate arrays ---*/
  
  // Geometric arrays
  delete [] UnitNormal;
  
  // Flow solution arrays
  delete [] dnvel;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] sigma[iDim];
  delete [] sigma;

  // Adjoint solution arrays
  for (iDim = 0; iDim < nDim; iDim++) {
    delete [] SigmaPhi[iDim];
    delete [] SigmaPsiE[iDim];
  }
  delete [] SigmaPhi;
  delete [] SigmaPsiE;
  
}

void CAdjTNE2NSSolver::BC_HeatFlux_Wall(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CNumerics *conv_numerics,
                                        CNumerics *visc_numerics,
                                        CConfig *config,
                                        unsigned short val_marker) {
  
  bool implicit;
  unsigned short iDim, iVar, jVar;
  unsigned long iPoint, iVertex, total_index;
  su2double ktr, kve;
  su2double dnPsiE, dnPsiEve;
  su2double *Normal, *d;
  su2double *dPdU, *dTdU, *dTvedU;
  su2double *Psi, *phi;
  su2double **GradPsi;
  
//	unsigned long iVertex, iPoint, total_index, Point_Normal;
//	unsigned short iDim, iVar, jVar, jDim;
//	su2double *d, *U, l1psi, mu_dyn, Temp, dVisc_T, rho, pressure, div_phi,
//  force_stress, Sigma_5, **PsiVar_Grad, phi[3];
//  su2double phis1, phis2, sq_vel, ProjVel, Enthalpy, *GridVel, phi_u, d_n;
//  su2double Energy, ViscDens, XiDens, Density, SoundSpeed, Pressure, dPhiE_dn, Laminar_Viscosity, Eddy_Viscosity,
//  Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_xz, Sigma_yz,
//  Sigma_xx5, Sigma_yy5, Sigma_zz5, Sigma_xy5, Sigma_xz5,
//  Sigma_yz5, eta_xx, eta_yy, eta_zz, eta_xy, eta_xz, eta_yz;
//  
//  
//  su2double *Psi = new su2double[nVar];
//	su2double **Tau = new su2double* [nDim];
//	for (iDim = 0; iDim < nDim; iDim++)
//		Tau[iDim] = new su2double [nDim];
//  su2double *Velocity = new su2double[nDim];
//  su2double *Normal = new su2double[nDim];
//  
//  su2double **GradPhi = new su2double* [nDim];
//  for (iDim = 0; iDim < nDim; iDim++)
//    GradPhi[iDim] = new su2double [nDim];
//  su2double *GradPsiE = new su2double [nDim];
  
  
  /*--- Set booleans ---*/
  implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  /*--- Allocate arrays ---*/
  Normal = new su2double[nDim];
  phi    = new su2double[nDim];
  Psi    = new su2double[nVar];

  /*--- Loop over boundary points ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
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
      
      /*--- Retrieve adjoint solution at the wall boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Psi[iVar] = node[iPoint]->GetSolution(iVar);
      GradPsi = node[iPoint]->GetGradient();
      
			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Get the force projection vector ---*/
      // Note: For force-based objective functions, this will be non-zero and
      //       for energy-based objectives, it will be zero.
			d = node[iPoint]->GetForceProj_Vector();
      for (iDim = 0; iDim < nDim; iDim++)
        phi[iDim] = d[iDim];
      
      /*--- Acquire flow quantities ---*/
      ktr = solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity();
      kve = solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve();
      dPdU   = solver_container[TNE2_SOL]->node[iPoint]->GetdPdU();
      dTdU   = solver_container[TNE2_SOL]->node[iPoint]->GetdTdU();
      dTvedU = solver_container[TNE2_SOL]->node[iPoint]->GetdTvedU();

      /*--- Weak imposition of the energy equation ---*/
      // Note: dn PsiE & dn PsiEve = 0.  We apply 'proportional control' to
      //       drive the boundary condition to satisfaction.
      dnPsiE   = 0.0;
      dnPsiEve = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        dnPsiE   += GradPsi[nSpecies+nDim][iDim]  *Normal[iDim];
        dnPsiEve += GradPsi[nSpecies+nDim+1][iDim]*Normal[iDim];
      }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Visc_i[iVar] = ktr*dTdU[iVar]  *dnPsiE +
                           kve*dTvedU[iVar]*(dnPsiE + dnPsiEve);
      }
      
      /*--- Convective terms ---*/
      // Energy
      for (iDim = 0; iDim < nDim; iDim++) {
        Res_Conv_i[nSpecies+nDim]   += phi[iDim]*Normal[iDim]*dPdU[nSpecies+nDim];
        Res_Conv_i[nSpecies+nDim+1] += phi[iDim]*Normal[iDim]*dPdU[nSpecies+nDim+1];
      }
      
      /*--- Apply the viscous residual ---*/
      LinSysRes.AddBlock(iPoint, Res_Conv_i);
      LinSysRes.SubtractBlock(iPoint, Res_Visc_i);
      
      
      
      /*--- Strong BC imposition for the adjoint velocity equations ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
      node[iPoint]->SetVel_ResTruncError_Zero();
			for (iDim = 0; iDim < nDim; iDim++)
				node[iPoint]->SetSolution_Old(nSpecies+iDim, phi[iDim]);
			if (implicit) {
				for (iVar = 0; iVar < nDim; iVar++) {
					total_index = iPoint*nVar+(nSpecies+iVar);
					Jacobian.DeleteValsRowi(total_index);
				}
			}
      
      
//      
//      
//      
//      
//      /*--- Energy resiudal due to the convective term ---*/
//      l1psi = 0.0;
//      for (iDim = 0; iDim < nDim; iDim++)
//        l1psi += Normal[iDim]*d[iDim];
//      Res_Conv_i[nVar-1] = l1psi*Gamma_Minus_One;
//      
//      /*--- Components of the effective and adjoint stress tensors ---*/
//      PsiVar_Grad = node[iPoint]->GetGradient();
//      div_phi = 0;
//      for (iDim = 0; iDim < nDim; iDim++) {
//        div_phi += PsiVar_Grad[iDim+1][iDim];
//        for (jDim = 0; jDim < nDim; jDim++)
//          Tau[iDim][jDim] = (PsiVar_Grad[iDim+1][jDim]+PsiVar_Grad[jDim+1][iDim]);
//      }
//      for (iDim = 0; iDim < nDim; iDim++)
//        Tau[iDim][iDim] -= TWO3*div_phi;
//      
//      /*--- force_stress = n_i \Tau_{ij} d_j ---*/
//      force_stress = 0.0;
//      for (iDim = 0; iDim < nDim; iDim++)
//        for (jDim = 0; jDim < nDim; jDim++)
//          force_stress += Normal[iDim]*Tau[iDim][jDim]*d[jDim];
//      
//      /*--- \partial \mu_dyn \partial T ---*/
//      mu_dyn = solver_container[TNE2_SOL]->node[iPoint]->GetLaminarViscosity();
//      Temp = solver_container[TNE2_SOL]->node[iPoint]->GetTemperature();
//      dVisc_T = 0.0;  // dVisc_T = mu_dyn*(Temp+3.0*mu2)/(2.0*Temp*(Temp+mu2));
//      
//      /*--- \Sigma_5 ---*/
//      Sigma_5 = (Gamma/Cp)*dVisc_T*force_stress;
//      
//      /*--- Imposition of residuals ---*/
//      rho = solver_container[TNE2_SOL]->node[iPoint]->GetDensity();
//      pressure = solver_container[TNE2_SOL]->node[iPoint]->GetPressure();
//      Res_Conv_i[0] = pressure*Sigma_5/(Gamma_Minus_One*rho*rho);
//      Res_Conv_i[nVar-1] -= Sigma_5/rho;
//      
//      /*--- Update convective and viscous residuals ---*/
//      LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
//      LinSysRes.SubtractBlock(iPoint, Res_Visc_i);
//      if (implicit) {
//        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
//      }
		}
	}
  
  
  
//	for (iDim = 0; iDim < nDim; iDim++)
//		delete [] Tau[iDim];
//	delete [] Tau;
//  delete [] Psi;
//  delete [] Velocity;
//  delete [] Normal;
//  delete [] GradPsiE;
//  for (iDim = 0; iDim < nDim; iDim++)
//    delete [] GradPhi[iDim];
//  delete [] GradPhi;
  
}


void CAdjTNE2NSSolver::BC_Isothermal_Wall(CGeometry *geometry,
                                          CSolver **solver_container,
                                          CNumerics *conv_numerics,
                                          CNumerics *visc_numerics,
                                          CConfig *config,
                                          unsigned short val_marker) {


  bool implicit, heat_flux_obj;
	unsigned long iVertex, iPoint, total_index, Point_Normal;
	unsigned short iDim, iVar, jVar;
  unsigned short RHOS_INDEX, RHO_INDEX, T_INDEX, TVE_INDEX;
	su2double *V, *dPdU, *d, q, dn;
  su2double *GradT, *GradTve;
  su2double ktr, kve, qtr, qve;
  su2double Area;
  su2double phi[3];
  su2double pnorm;
  su2double *Psi, *Normal, UnitNormal[3];
  
  /*--- Set booleans from CConfig specifications ---*/
  implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  heat_flux_obj  = ((config->GetKind_ObjFunc() == INVERSE_DESIGN_HEATFLUX) ||
                    (config->GetKind_ObjFunc() == TOTAL_HEATFLUX) ||
                    (config->GetKind_ObjFunc() == MAXIMUM_HEATFLUX));
  
  /*--- Allocate arrays ---*/
  Psi = new su2double[nVar];
  Normal = new su2double[nDim];
  
  /*--- Get primitive vector locators ---*/
  RHOS_INDEX = solver_container[TNE2_SOL]->node[0]->GetRhosIndex();
  RHO_INDEX  = solver_container[TNE2_SOL]->node[0]->GetRhoIndex();
  T_INDEX = solver_container[TNE2_SOL]->node[0]->GetTIndex();
  TVE_INDEX = solver_container[TNE2_SOL]->node[0]->GetTveIndex();
  
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
      V = solver_container[TNE2_SOL]->node[iPoint]->GetPrimitive();
      dPdU = solver_container[TNE2_SOL]->node[iPoint]->GetdPdU();
      
			/*--- Normal vector for this vertex ---*/
      // Note: Convention is outward facing normal
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
        Normal[iDim] = -Normal[iDim];
        Area += Normal[iDim]*Normal[iDim];
      }
      Area = sqrt(Area);
      for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = Normal[iDim]/Area;
      
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
      
      /*--- If heat-flux objective, determine appropriate energy B.C. ---*/
      if (heat_flux_obj) {
        
        /*--- Read from config file ---*/
        pnorm = config->GetPnormHeat();
        
        /*--- Determine local heat flux ---*/
        ktr = solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity();
        kve = solver_container[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve();
        GradT   = solver_container[TNE2_SOL]->node[iPoint]->GetGradient_Primitive()[T_INDEX];
        GradTve = solver_container[TNE2_SOL]->node[iPoint]->GetGradient_Primitive()[TVE_INDEX];
        qtr = 0.0;
        qve = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          qtr += ktr*GradT[iDim]*UnitNormal[iDim];
          qve += kve*GradTve[iDim]*UnitNormal[iDim];
        }
        q = -pnorm * pow(qtr+qve, pnorm-1.0) * Area;
        
      } else {
        q = 0.0;
      }
      
      /*--- Apply the boundary condition to the linear system ---*/
      LinSysRes.SetBlock_Zero(iPoint, nSpecies+nDim);
      LinSysRes.SetBlock_Zero(iPoint, nSpecies+nDim+1);
      node[iPoint]->SetSolution_Old(nSpecies+nDim,   q);
      node[iPoint]->SetSolution_Old(nSpecies+nDim+1, 0.0);
      if (implicit) {
        iVar = nSpecies+nDim;
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
        Jacobian.DeleteValsRowi(total_index+1);
      }
    }
  }
  
  
  delete [] Psi;
  delete [] Normal;
}



