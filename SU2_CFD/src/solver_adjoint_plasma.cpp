/*!
 * \file solution_adjoint_plasma.cpp
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

CAdjPlasmaSolver::CAdjPlasmaSolver(void) : CSolver() { }

CAdjPlasmaSolver::CAdjPlasmaSolver(CGeometry *geometry, CConfig *config) : CSolver() {
	unsigned long iPoint, iVertex;
	string text_line, mesh_filename;
	unsigned short iDim = 0, iVar, iMarker;
	ifstream restart_file;
	string filename, AdjExt;

	bool restart = config->GetRestart();
	bool axisymmetric = config->GetAxisymmetric();

	/*--- Define geometry constans in the solver structure ---*/
	nDim        = geometry->GetnDim();
	nMonatomics = config->GetnMonatomics();
	nDiatomics  = config->GetnDiatomics();
	nSpecies    = config->GetnSpecies();
	nVar        = nMonatomics*(nDim+2) + nDiatomics*(nDim+3);
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
	node = new CVariable*[geometry->GetnPoint()];

	/*--- Define some auxiliary vectors related to the residual ---*/
	Residual = new double[nVar];	  Residual_RMS = new double[nVar];
	Residual_i = new double[nVar];	Residual_j = new double[nVar];
	Res_Conv_i = new double[nVar];	Res_Visc_i = new double[nVar];
	Res_Conv_j = new double[nVar];	Res_Visc_j = new double[nVar];
  Residual_Max = new double[nVar]; Point_Max = new unsigned long[nVar];

	Residual_Chemistry = new double[nVar]; Residual_MomentumExch = new double[nVar];
	Residual_ElecForce = new double[nVar]; Residual_EnergyExch = new double[nVar];
	Res_Conv = new double[nVar];	Res_Visc = new double[nVar];	Res_Sour = new double[nVar];
	if (axisymmetric) {
		Residual_Axisymmetric = new double[nVar];
	}
	

	/*--- Define some auxiliary vectors related to the solution ---*/
	Solution   = new double[nVar];
	Solution_i = new double[nVar]; Solution_j = new double[nVar];
  
	/*--- Define some auxiliary vectors related to the geometry ---*/
//  Vector = new double[nDim];
	Vector_i = new double[nDim]; Vector_j = new double[nDim];

  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT) {
		Jacobian_ii = new double* [nVar];
		Jacobian_ij = new double* [nVar];
		Jacobian_ji = new double* [nVar];
		Jacobian_jj = new double* [nVar];
		Jacobian_Chemistry    = new double* [nVar];
		Jacobian_ElecForce    = new double* [nVar];
		Jacobian_MomentumExch = new double* [nVar];
		Jacobian_EnergyExch   = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_ii[iVar] = new double [nVar];
			Jacobian_ij[iVar] = new double [nVar];
			Jacobian_ji[iVar] = new double [nVar];
			Jacobian_jj[iVar] = new double [nVar];
			Jacobian_Chemistry[iVar]    = new double [nVar];
			Jacobian_ElecForce[iVar]    = new double [nVar];
			Jacobian_MomentumExch[iVar] = new double [nVar];
			Jacobian_EnergyExch[iVar]   = new double [nVar];
			
		}
		if (axisymmetric) {
			Jacobian_Axisymmetric = new double*[nVar];
			for (iVar = 0; iVar < nVar; iVar++)
				Jacobian_Axisymmetric[iVar] = new double[nVar];
		}
		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);

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
	CSensitivity = new double* [config->GetnMarker_All()];
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		CSensitivity[iMarker] = new double [geometry->nVertex[iMarker]];
	}
	Sens_Geo  = new double[config->GetnMarker_All()];
	Sens_Mach = new double[config->GetnMarker_All()];
	Sens_AoA  = new double[config->GetnMarker_All()];
	Sens_Press = new double[config->GetnMarker_All()];
	Sens_Temp  = new double[config->GetnMarker_All()];
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Sens_Geo[iMarker]  = 0.0;
		Sens_Mach[iMarker] = 0.0;
		Sens_AoA[iMarker]  = 0.0;
		Sens_Press[iMarker] = 0.0;
		Sens_Temp[iMarker]  = 0.0;
		for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
			CSensitivity[iMarker][iVertex] = 0.0;
	}

	/*--- Adjoint flow at the inifinity, initialization stuff ---*/
	PsiRho_Inf = 0.0; PsiE_Inf   = 0.0; PsiEvib_Inf = 0.0;
	Phi_Inf    = new double [nDim];
	Phi_Inf[0] = 0.0; Phi_Inf[1] = 0.0;
	if (nDim == 3) Phi_Inf[2] = 0.0;		
  
	if (!restart || geometry->GetFinestMGLevel() == false) {
		/*--- Restart the solution from infinity ---*/
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CAdjPlasmaVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, PsiEvib_Inf, nDim, nVar, config);
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
    for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    }
    
		/*--- Read all lines in the restart file ---*/
    long iPoint_Local; unsigned long iPoint_Global = 0; unsigned long index;
    
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

        /*--- First value is the point index, then the conservative vars. ---*/
        point_line >> index;
        for (iVar = 0; iVar < nVar; iVar++)
          point_line >> Solution[iVar];
        
        node[iPoint_Local] = new CAdjPlasmaVariable(Solution, nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
      node[iPoint] = new CAdjPlasmaVariable(Solution, nDim, nVar, config);
    }
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
	}
  
  /*--- Define solver parameters needed for execution of destructor ---*/
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTERED) space_centered = true;
	else space_centered = false;
  
  /*--- Send solution information using MPI ---*/
  Set_MPI_Solution(geometry, config);
  
}

CAdjPlasmaSolver::~CAdjPlasmaSolver(void) {
	unsigned short iVar, iDim;

	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Jacobian_ii[iVar]; delete [] Jacobian_ij[iVar];
		delete [] Jacobian_ji[iVar]; delete [] Jacobian_jj[iVar];
	}
	delete [] Jacobian_ii; delete [] Jacobian_ij;
	delete [] Jacobian_ji; delete [] Jacobian_jj;
	
	delete [] Residual; delete [] Residual_Max;
	delete [] Residual_i; delete [] Residual_j;
	delete [] Res_Conv_i; delete [] Res_Visc_i;
	delete [] Res_Conv_j; delete [] Res_Visc_j;
	delete [] Solution; 
	delete [] Solution_i; delete [] Solution_j;
//  delete [] Vector;
	delete [] Vector_i; delete [] Vector_j;
	delete [] Sens_Geo; delete [] Sens_Mach;
	delete [] Sens_AoA; delete [] Sens_Press;
	delete [] Sens_Temp; delete [] Phi_Inf;
	
  if (space_centered) {
		delete [] p1_Und_Lapl;
		delete [] p2_Und_Lapl;
	}
  
	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;
	
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
	
/*	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		delete [] CSensitivity[iMarker];
	 delete [] CSensitivity; */
	
/*	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		delete [] node[iPoint];
	delete [] node; */
}


void CAdjPlasmaSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  
	unsigned short iVar, iMarker, iPeriodic_Index, iSpecies, loc;
	unsigned long iVertex, iPoint, nVertex, nBuffer_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi,
	sinPsi, *newSolution = NULL, *Buffer_Receive_U = NULL;
	short SendRecv;
	int send_to, receive_from;
  
#ifndef NO_MPI
	double *Buffer_Send_U = NULL, *Buffer_Send_LaminarViscosity = NULL,
  *Buffer_Send_EddyViscosity = NULL, *Buffer_Send_VGrad = NULL, *Buffer_Send_UGrad = NULL, *Buffer_Send_Limit = NULL, *Buffer_Send_Undivided_Laplacian = NULL,
  *Buffer_Send_Sensor = NULL, *Buffer_Send_Lambda = NULL;
	unsigned short *Buffer_Send_Neighbor = NULL;
#endif
  
	newSolution = new double[nVar];
  
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
      
			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			nVertex = geometry->nVertex[iMarker];
      
			nBuffer_Vector			= nVertex*nVar;
      
			send_to = SendRecv-1;
			receive_from = abs(SendRecv)-1;
      
#ifndef NO_MPI
      
			/*--- Send information using MPI  ---*/
			if (SendRecv > 0) {
				/*--- Allocate upwind variables ---*/
				Buffer_Send_U = new double[nBuffer_Vector];
        
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
          
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
					/*--- Copy data ---*/
					for (iVar = 0; iVar < nVar; iVar++) {
						Buffer_Send_U[iVar*nVertex+iVertex] = node[iPoint]->GetSolution(iVar);
					}
				}
        
				/*--- Send the buffer, and deallocate information using MPI ---*/
				MPI::COMM_WORLD.Bsend(Buffer_Send_U, nBuffer_Vector, MPI::DOUBLE, send_to, 0); delete [] Buffer_Send_U;
			}
      
#endif
      
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
        
				/*--- Allocate upwind variables ---*/
				Buffer_Receive_U = new double [nBuffer_Vector];
        
        
#ifdef NO_MPI
        
				/*--- Get the information from the donor point directly. This is a
				 serial computation with access to all nodes. Note that there is an
				 implicit ordering in the list. ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
          
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
					/*--- Copy data ---*/
					for (iVar = 0; iVar < nVar; iVar++) {
						Buffer_Receive_U[iVar*nVertex+iVertex] = node[iPoint]->GetSolution(iVar);
					}
				}
        
#else
				/*--- Receive the information using MPI---*/
				MPI::COMM_WORLD.Recv(Buffer_Receive_U, nBuffer_Vector, MPI::DOUBLE, receive_from, 0);
        
#endif
        
				/*--- Do the coordinate transformation ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
          
					/*--- Find point and its type of transformation ---*/
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					iPeriodic_Index = geometry->vertex[iMarker][iVertex]->GetRotation_Type();
          
					/*--- Retrieve the supplied periodic information. ---*/
					angles = config->GetPeriodicRotation(iPeriodic_Index);
          
					/*--- Store angles separately for clarity. ---*/
					theta    = angles[0];   phi    = angles[1]; psi    = angles[2];
					cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
					sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
          
					/*--- Compute the rotation matrix. Note that the implicit
					 ordering is rotation about the x-axis, y-axis,
					 then z-axis. Note that this is the transpose of the matrix
					 used during the preprocessing stage. ---*/
					rotMatrix[0][0] = cosPhi*cosPsi; rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi; rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
					rotMatrix[0][1] = cosPhi*sinPsi; rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi; rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
					rotMatrix[0][2] = -sinPhi; rotMatrix[1][2] = sinTheta*cosPhi; rotMatrix[2][2] = cosTheta*cosPhi;
          
					/*--- Copy conserved variables before performing transformation. ---*/
					for (iVar = 0; iVar < nVar; iVar++)
						newSolution[iVar] = Buffer_Receive_U[iVar*nVertex+iVertex];
          
          
					/*--- Rotate the momentum components. ---*/
					if (nDim == 2) {
						for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
							if (iSpecies < nDiatomics) loc = (nDim+3)*iSpecies;
							else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
							newSolution[loc + 1] = rotMatrix[0][0]*Buffer_Receive_U[(loc + 1)*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_U[(loc + 2)*nVertex+iVertex];
							newSolution[loc + 2] = rotMatrix[1][0]*Buffer_Receive_U[(loc + 1)*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_U[(loc + 2)*nVertex+iVertex];
						}
					}
					else {
						for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
							if (iSpecies < nDiatomics) loc = (nDim+3)*iSpecies;
							else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
							newSolution[loc + 1] = rotMatrix[0][0]*Buffer_Receive_U[(loc + 1)*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_U[(loc + 2)*nVertex+iVertex] + rotMatrix[0][2]*Buffer_Receive_U[(loc + 3)*nVertex+iVertex];
							newSolution[loc + 2] = rotMatrix[1][0]*Buffer_Receive_U[(loc + 1)*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_U[(loc + 2)*nVertex+iVertex] + rotMatrix[1][2]*Buffer_Receive_U[(loc + 3)*nVertex+iVertex];
							newSolution[loc + 3] = rotMatrix[2][0]*Buffer_Receive_U[(loc + 1)*nVertex+iVertex] + rotMatrix[2][1]*Buffer_Receive_U[(loc + 2)*nVertex+iVertex] + rotMatrix[2][2]*Buffer_Receive_U[(loc + 3)*nVertex+iVertex];
						}
					}
					/*--- Copy transformed conserved variables back into buffer. ---*/
					for (iVar = 0; iVar < nVar; iVar++)
						Buffer_Receive_U[iVar*nVertex+iVertex] = newSolution[iVar];
          
          
					/*--- Upwind method. Store the received information ---*/
					for (iVar = 0; iVar < nVar; iVar++) {
						node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertex+iVertex]);
					}
				}
				delete [] Buffer_Receive_U;
			}
		}
	}
	delete [] newSolution;
}

void CAdjPlasmaSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, iMarker, iPeriodic_Index;
	unsigned long iVertex, iPoint, nVertex, nBuffer_VectorGrad;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi,
	sinPsi, **newGradient = NULL, *Buffer_Receive_UGrad = NULL;
	short SendRecv;
	int send_to, receive_from;
    
#ifndef NO_MPI
    
    MPI::COMM_WORLD.Barrier();
	double *Buffer_Send_UGrad = NULL;
    
#endif
    
	newGradient = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		newGradient[iVar] = new double[3];
    
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			nVertex = geometry->nVertex[iMarker];
			nBuffer_VectorGrad = nVertex*nVar*nDim;
			send_to = SendRecv-1;
			receive_from = abs(SendRecv)-1;
            
#ifndef NO_MPI
            
			/*--- Send information using MPI  ---*/
			if (SendRecv > 0) {
                Buffer_Send_UGrad = new double[nBuffer_VectorGrad];
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    for (iVar = 0; iVar < nVar; iVar++)
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Send_UGrad[iDim*nVar*nVertex+iVar*nVertex+iVertex] = node[iPoint]->GetGradient(iVar,iDim);
				}
                MPI::COMM_WORLD.Bsend(Buffer_Send_UGrad, nBuffer_VectorGrad, MPI::DOUBLE, send_to, 0); delete [] Buffer_Send_UGrad;
			}
            
#endif
            
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
                Buffer_Receive_UGrad = new double [nBuffer_VectorGrad];
                
#ifdef NO_MPI
                
				/*--- Receive information without MPI ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    for (iVar = 0; iVar < nVar; iVar++)
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Receive_UGrad[iDim*nVar*nVertex+iVar*nVertex+iVertex] = node[iPoint]->GetGradient(iVar,iDim);
				}
                
#else
                
                MPI::COMM_WORLD.Recv(Buffer_Receive_UGrad, nBuffer_VectorGrad, MPI::DOUBLE, receive_from, 0);
                
#endif
                
				/*--- Do the coordinate transformation ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
                    
					/*--- Find point and its type of transformation ---*/
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					iPeriodic_Index = geometry->vertex[iMarker][iVertex]->GetRotation_Type();
                    
					/*--- Retrieve the supplied periodic information. ---*/
					angles = config->GetPeriodicRotation(iPeriodic_Index);
                    
					/*--- Store angles separately for clarity. ---*/
					theta    = angles[0];   phi    = angles[1]; psi    = angles[2];
					cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
					sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
                    
					/*--- Compute the rotation matrix. Note that the implicit
					 ordering is rotation about the x-axis, y-axis,
					 then z-axis. Note that this is the transpose of the matrix
					 used during the preprocessing stage. ---*/
					rotMatrix[0][0] = cosPhi*cosPsi; rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi; rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
					rotMatrix[0][1] = cosPhi*sinPsi; rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi; rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
					rotMatrix[0][2] = -sinPhi; rotMatrix[1][2] = sinTheta*cosPhi; rotMatrix[2][2] = cosTheta*cosPhi;
                    
                    for (iVar = 0; iVar < nVar; iVar++)
                        for (iDim = 0; iDim < nDim; iDim++)
                            newGradient[iVar][iDim] = Buffer_Receive_UGrad[iDim*nVar*nVertex+iVar*nVertex+iVertex];
                    
                    /*--- Need to rotate the gradients for all conserved variables. ---*/
                    for (iVar = 0; iVar < nVar; iVar++) {
                        if (nDim == 2) {
                            newGradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_UGrad[0*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_UGrad[1*nVar*nVertex+iVar*nVertex+iVertex];
                            newGradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_UGrad[0*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_UGrad[1*nVar*nVertex+iVar*nVertex+iVertex];
                        }
                        else {
                            newGradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_UGrad[0*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_UGrad[1*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[0][2]*Buffer_Receive_UGrad[2*nVar*nVertex+iVar*nVertex+iVertex];
                            newGradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_UGrad[0*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_UGrad[1*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[1][2]*Buffer_Receive_UGrad[2*nVar*nVertex+iVar*nVertex+iVertex];
                            newGradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_UGrad[0*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[2][1]*Buffer_Receive_UGrad[1*nVar*nVertex+iVar*nVertex+iVertex] + rotMatrix[2][2]*Buffer_Receive_UGrad[2*nVar*nVertex+iVar*nVertex+iVertex];
                        }
                    }
                    
                    /*--- Copy transformed gradients back into buffer. ---*/
                    for (iVar = 0; iVar < nVar; iVar++)
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Receive_UGrad[iDim*nVar*nVertex+iVar*nVertex+iVertex] = newGradient[iVar][iDim];
                    
                    
					/*--- Store the received information ---*/
                    for (iVar = 0; iVar < nVar; iVar++)
                        for (iDim = 0; iDim < nDim; iDim++)
                            node[iPoint]->SetGradient(iVar, iDim, Buffer_Receive_UGrad[iDim*nVar*nVertex+iVar*nVertex+iVertex]);
				}
                delete [] Buffer_Receive_UGrad;
			}
		}
	}
    
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] newGradient[iVar];
	delete [] newGradient;
    
#ifndef NO_MPI
    
    MPI::COMM_WORLD.Barrier();
    
#endif
    
}

void CAdjPlasmaSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
	
	bool implicit   = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  bool upwind_2nd = ( (config->GetKind_Upwind_AdjPlasma() == ROE_2ND) ||
                      (config->GetKind_Upwind_AdjPlasma() == SW_2ND)  ||
                      (config->GetKind_Upwind_AdjPlasma() == MSW_2ND)    );
  bool limiter    = (config->GetKind_SlopeLimit_AdjPlasma() != NONE);
  
	/*--- Residual initialization ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		
		/*--- Initialize the convective residual vector ---*/
		LinSysRes.SetBlock_Zero(iPoint);
		
	}
  
  /*--- Upwind 2nd order flux reconstruction ---*/
  if ( (upwind_2nd) && ((iMesh == MESH_0))) {
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
    else if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    
    /*--- Limiter computation ---*/
    if (limiter) SetSolution_Limiter(geometry, config);
  }
	
	/*--- Implicit solution ---*/
	if (implicit) Jacobian.SetValZero();
	
}

void CAdjPlasmaSolver::SetForceProj_Vector(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	double *ForceProj_Vector, x = 0.0, y = 0.0, z = 0.0, *Normal;
	unsigned short iMarker;
	unsigned long iVertex, iPoint;
	double Alpha      = (config->GetAoA()*PI_NUMBER)/180.0;
	double Beta       = (config->GetAoS()*PI_NUMBER)/180.0;

	ForceProj_Vector = new double [nDim];

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
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
						if (nDim == 2) { ForceProj_Vector[0] = cos(Alpha); ForceProj_Vector[1] = sin(Alpha); }
						if (nDim == 3) { ForceProj_Vector[0] = cos(Alpha)*cos(Beta); ForceProj_Vector[1] = sin(Beta); ForceProj_Vector[2] = sin(Alpha)*cos(Beta); }
						break;
					case LIFT_COEFFICIENT :
						if (nDim == 2) { ForceProj_Vector[0] = -sin(Alpha); ForceProj_Vector[1] = cos(Alpha); }
						if (nDim == 3) { ForceProj_Vector[0] = -sin(Alpha); ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = cos(Alpha); }
						break;
					case SIDEFORCE_COEFFICIENT :
						if (nDim == 2) { cout << "This functional is not possible in 2D!!" << endl;
							cout << "Press any key to exit..." << endl; cin.get(); exit(1);
						}
						if (nDim == 3) { ForceProj_Vector[0] = -sin(Beta) * cos(Alpha); ForceProj_Vector[1] = cos(Beta); ForceProj_Vector[2] = -sin(Beta) * sin(Alpha); }
						break;
          case FORCE_X_COEFFICIENT :
						if (nDim == 2) { ForceProj_Vector[0] = 1.0; ForceProj_Vector[1] = 0.0; }
						if (nDim == 3) { ForceProj_Vector[0] = 1.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 0.0; }
						break;
          case FORCE_Y_COEFFICIENT :
						if (nDim == 2) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 1.0; }
						if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 1.0; ForceProj_Vector[2] = 0.0; }
						break;
          case FORCE_Z_COEFFICIENT :
						if (nDim == 2) {cout << "This functional is not possible in 2D!!" << endl;
							cout << "Press any key to exit..." << endl;
							cin.get(); exit(1);
            }
						if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 1.0; }
						break;
				}
        
        /*--- Store the force projection vector at this node ---*/
				node[iPoint]->SetForceProj_Vector(ForceProj_Vector);
			}
	
	delete [] ForceProj_Vector;
}

void CAdjPlasmaSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
																			 CConfig *config, unsigned short iMesh, unsigned short iRKStep) { 	
	unsigned long iEdge, iPoint, jPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	unsigned short iSpecies;
	
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		
		/*--- Points in edge, set normal vectors, and number of neighbors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
		
		numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());
		
		/*--- Adjoint variables w/o reconstruction ---*/
		numerics->SetAdjointVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
		
		/*--- Set conservative variables w/o reconstruction ---*/
		numerics->SetConservative(solver_container[PLASMA_SOL]->node[iPoint]->GetSolution(), 
														solver_container[PLASMA_SOL]->node[jPoint]->GetSolution());
		
		
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			numerics->SetPressure(solver_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies), 
													solver_container[PLASMA_SOL]->node[jPoint]->GetPressure(iSpecies), iSpecies);
			numerics->SetSoundSpeed(solver_container[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies), 
														solver_container[PLASMA_SOL]->node[jPoint]->GetSoundSpeed(iSpecies),iSpecies);
			numerics->SetEnthalpy(solver_container[PLASMA_SOL]->node[iPoint]->GetEnthalpy(iSpecies), 
													solver_container[PLASMA_SOL]->node[jPoint]->GetEnthalpy(iSpecies),iSpecies);

      switch (config->GetKind_ConvNumScheme_Plasma()) {
        case SPACE_CENTERED:
          numerics->SetLambda(solver_container[PLASMA_SOL]->node[iPoint]->GetLambda(iSpecies), solver_container[PLASMA_SOL]->node[jPoint]->GetLambda(iSpecies),iSpecies);
          break;
        case SPACE_UPWIND:
          numerics->SetLambda(solver_container[PLASMA_SOL]->node[iPoint]->GetMax_Lambda_Inv(iSpecies), solver_container[PLASMA_SOL]->node[jPoint]->GetMax_Lambda_Inv(iSpecies),iSpecies);
          break;
      }
    }
				
		/*--- Compute residuals ---*/
		numerics->ComputeResidual(Res_Conv_i, Res_Visc_i, Res_Conv_j, Res_Visc_j, 
												Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
				
//		cout << Res_Visc_i[0] <<" "<< Res_Visc_i[1] <<" "<< Res_Visc_i[2] <<" "<< Res_Visc_i[3] <<" "<< Res_Visc_i[4] <<" "<< Res_Visc_i[5] <<" "<< Res_Visc_i[6] <<" "<< Res_Visc_i[7] <<endl;
//		cout << Res_Visc_j[0] <<" "<< Res_Visc_j[1] <<" "<< Res_Visc_j[2] <<" "<< Res_Visc_j[3] <<" "<< Res_Visc_j[4] <<" "<< Res_Visc_j[5] <<" "<< Res_Visc_j[6] <<" "<< Res_Visc_j[7] <<endl;

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

void CAdjPlasmaSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
																				 CConfig *config, unsigned short iMesh) {
	double *Psi_i, *Psi_j, *U_i, *U_j;
  double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
  double *Limiter_i = NULL, *Limiter_j = NULL;
	unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim, iVar;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  bool high_order_diss = (( (config->GetKind_Upwind_AdjPlasma() == ROE_2ND)  ||
                            (config->GetKind_Upwind_AdjPlasma() == SW_2ND)   ||
                            (config->GetKind_Upwind_AdjPlasma() == MSW_2ND)) &&
                            (iMesh == MESH_0));
  
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

		/*--- Adjoint variables w/o reconstruction ---*/
		Psi_i = node[iPoint]->GetSolution(); Psi_j = node[jPoint]->GetSolution();
		numerics->SetAdjointVar(Psi_i, Psi_j);
		
		/*--- Conservative variables w/o reconstruction ---*/
		U_i = solver_container[PLASMA_SOL]->node[iPoint]->GetSolution();
		U_j = solver_container[PLASMA_SOL]->node[jPoint]->GetSolution();
		numerics->SetConservative(U_i, U_j);
    
    
    if ((high_order_diss) && (config->GetKind_Adjoint() != DISCRETE) ) {
      
      cout << "CAdjPlasmaSolver::Upwind_Residual - FLUX RECONSTRUCTION NOT VERIFIED!!!" << endl;
			for (iDim = 0; iDim < nDim; iDim++) {
				Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
				Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}
      
			Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
			if (config->GetKind_SlopeLimit() != NONE) {
				Limiter_j = node[jPoint]->GetLimiter(); Limiter_i = node[iPoint]->GetLimiter();
			}
      
			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				if (config->GetKind_SlopeLimit() == NONE) {
					Solution_i[iVar] = U_i[iVar] + Project_Grad_i;
					Solution_j[iVar] = U_j[iVar] + Project_Grad_j;
				}
				else {
					Solution_i[iVar] = U_i[iVar] + Project_Grad_i*Limiter_i[iVar];
					Solution_j[iVar] = U_j[iVar] + Project_Grad_j*Limiter_j[iVar];
				}
			}
			numerics->SetAdjointVar(Solution_i, Solution_j);
		}
		
		numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
		
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

void CAdjPlasmaSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
																				 CConfig *config, unsigned short iMesh) {
	
	unsigned short iVar, jVar;
	unsigned long iPoint;
	double *Psi_i;
	bool implicit = (config->GetKind_TimeIntScheme_AdjPlasma() == EULER_IMPLICIT);
	bool axisymmetric = config->GetAxisymmetric();	
	
	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		for (iVar = 0; iVar < nVar; iVar++) {
			Residual[iVar] = 0.0;
			Residual_Chemistry[iVar] = 0.0;
			Residual_MomentumExch[iVar] = 0.0;
			Residual_ElecForce[iVar] = 0.0;
			Residual_EnergyExch[iVar] = 0.0;
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_Chemistry[iVar][jVar] = 0.0;
				Jacobian_MomentumExch[iVar][jVar] = 0.0;
				Jacobian_ElecForce[iVar][jVar] = 0.0;
				Jacobian_EnergyExch[iVar][jVar] = 0.0;
			}
		}
		if (axisymmetric) {
			for (iVar = 0; iVar < nVar; iVar++) {
				Residual_Axisymmetric[iVar] = 0.0;
				for (jVar = 0; jVar < nVar; jVar++) 
					Jacobian_Axisymmetric[iVar][jVar] = 0.0;
			}
		}
    
    /*--- Set y coordinate ---*/
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),geometry->node[iPoint]->GetCoord());
    
    /*--- Set solution  ---*/
    numerics->SetConservative(solver_container[PLASMA_SOL]->node[iPoint]->GetSolution(), solver_container[PLASMA_SOL]->node[iPoint]->GetSolution());
    
    /*--- Set control volume ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    
    /*--- Set all primitive variables and gradients ---*/
		numerics->SetPrimitive(solver_container[PLASMA_SOL]->node[iPoint]->GetPrimVar_Plasma(), solver_container[PLASMA_SOL]->node[iPoint]->GetPrimVar_Plasma());
		numerics->SetPrimVarGradient(solver_container[PLASMA_SOL]->node[iPoint]->GetGradient_Primitive_Plasma(),solver_container[PLASMA_SOL]->node[iPoint]->GetGradient_Primitive_Plasma());
    
    /*--- Load auxiliary vector with local adjoint variables ---*/
    Psi_i = node[iPoint]->GetSolution();		    
    
    /*--- Axisymmetric source terms ---*/
    if (axisymmetric) {
      numerics->SetJacobian_Axisymmetric(Jacobian_Axisymmetric, config);			
      for (iVar = 0; iVar < nVar; iVar ++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          Residual_Axisymmetric[iVar] += Jacobian_Axisymmetric[jVar][iVar]*Psi_i[jVar];
          Jacobian_ii[iVar][jVar] = Jacobian_Axisymmetric[jVar][iVar];
        }
      }
      LinSysRes.AddBlock(iPoint, Residual_Axisymmetric);
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);	
    }
    
    /*--- Chemistry source terms ---*/
    numerics->SetJacobian_Chemistry(Jacobian_Chemistry, config);
    for (iVar = 0; iVar < nVar; iVar ++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Residual_Chemistry[iVar] += Jacobian_Chemistry[jVar][iVar]*Psi_i[jVar];
        Jacobian_ii[iVar][jVar] = Jacobian_Chemistry[jVar][iVar];
      }
    }
    LinSysRes.SubtractBlock(iPoint, Residual_Chemistry);
    if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);	
    
    /*--- Momentum exchange source terms ---*/
    numerics->SetJacobian_MomentumExch(Jacobian_MomentumExch, config);			
    for (iVar = 0; iVar < nVar; iVar ++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Residual_MomentumExch[iVar] += Jacobian_MomentumExch[jVar][iVar]*Psi_i[jVar];
        Jacobian_ii[iVar][jVar] = Jacobian_MomentumExch[jVar][iVar];
      }
    }
    LinSysRes.SubtractBlock(iPoint, Residual_MomentumExch);
    if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
    
    /*--- Energy exchange source terms ---*/
    numerics->SetJacobian_EnergyExch(Jacobian_EnergyExch, config);
    for (iVar = 0; iVar < nVar; iVar ++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Residual_EnergyExch[iVar] += Jacobian_EnergyExch[jVar][iVar]*Psi_i[jVar];
        Jacobian_ii[iVar][jVar] = Jacobian_EnergyExch[jVar][iVar];
      }
    }
    LinSysRes.SubtractBlock(iPoint, Residual_EnergyExch);
    if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
  }
}

void CAdjPlasmaSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
																				 CConfig *config, unsigned short iMesh) {
}


void CAdjPlasmaSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	unsigned short iVar, iSpecies;
	unsigned long iPoint, total_index;
	double Delta, Delta_flow, *local_Res_TruncError, Vol;
    bool MultipleTimeSteps = config->MultipleTimeSteps();
    double *Species_Delta;
    
    Species_Delta = new double [nSpecies];
    
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
    }
	
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
        local_Res_TruncError = node[iPoint]->GetResTruncError();
        Vol = geometry->node[iPoint]->GetVolume();
        
        /*--- Modify matrix diagonal to assure diagonal dominance ---*/
        if (!MultipleTimeSteps) {
            Delta_flow = Vol/(solver_container[PLASMA_SOL]->node[iPoint]->GetDelta_Time());
            Delta = Delta_flow;
            Jacobian.AddVal2Diag(iPoint, Delta);
        } else {
            for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
                Species_Delta[iSpecies] = Vol/(solver_container[PLASMA_SOL]->node[iPoint]->GetDelta_Time(iSpecies));
            Jacobian.AddVal2Diag(iPoint, Species_Delta, nDim, nDiatomics);
        }
        
        
        for (iVar = 0; iVar < nVar; iVar++) {
            total_index = iPoint*nVar+iVar;
            /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
            LinSysRes[total_index] = -(LinSysRes[total_index]+local_Res_TruncError[iVar]);
            LinSysSol[total_index] = 0.0;
            AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
            AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex());
        }
    }
    
    /*--- Initialize residual and solution at the ghost points ---*/
    for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
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
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		for (iVar = 0; iVar < nVar; iVar++)
            node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar]);
	
    /*--- MPI solution ---*/
    Set_MPI_Solution(geometry, config);
    
    /*--- Compute the root mean square residual ---*/
    SetResidual_RMS(geometry, config);
	
    delete [] Species_Delta;
}

void CAdjPlasmaSolver::Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
  unsigned long iVertex, iPoint, Neigh;
	//unsigned short iPos, jPos;
	unsigned short iDim, iMarker, iNeigh;
  double *d = NULL, *Normal = NULL, *Psi = NULL, *U = NULL, Enthalpy, conspsi,
	Area, **PrimVar_Grad = NULL, *ConsPsi_Grad = NULL, ConsPsi, d_press, grad_v,
	v_gradconspsi;
  //double UnitaryNormal[3], *RotVel = NULL, *GridVel = NULL;
  //double Mach_Inf, Beta2;
	//double RefDensity, *RefVelocity = NULL, RefPressure;
  
	//double r, ru, rv, rw, rE, p; // used in farfield sens
	//double dp_dr, dp_dru, dp_drv, dp_drw, dp_drE; // used in farfield sens
	//double dH_dr, dH_dru, dH_drv, dH_drw, dH_drE, H; // used in farfield sens
	//double alpha, beta;
	//double *USens, *U_infty;
  
  /*--- Loop over boundary markers to select those for Euler walls ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
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
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
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
          ConsPsi_Grad = node[iPoint]->GetAuxVarGradient();
          ConsPsi = node[iPoint]->GetAuxVar();
          
          d_press = 0; grad_v = 0; v_gradconspsi = 0;
          for (iDim = 0; iDim < nDim; iDim++) {
            d_press += d[iDim]*PrimVar_Grad[nDim+1][iDim];
            grad_v += PrimVar_Grad[iDim+1][iDim]*ConsPsi;
//            v_gradconspsi += solver_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim, incompressible) * ConsPsi_Grad[iDim];
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

void CAdjPlasmaSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {
	
	unsigned long iVertex, iPoint;
	double *d, *Normal, *U, *Psi_Aux, ProjVel, bcn, Area, *UnitaryNormal, *Coord, Gamma_Minus_One;
  double *Velocity, *Psi, Enthalpy, hf, Energy_vib, sq_vel, phin;
	unsigned short iDim, iVar, jVar, jDim, loc, iSpecies;
  double phidotu, phidotn, Energy_el, dPdrho, *dPdrhou, dPdrhoE, dPdrhoEv;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	
	UnitaryNormal = new double[nDim];
	Velocity = new double[nDim];
	Psi      = new double[nVar];
  dPdrhou = new double[nDim];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		if (geometry->node[iPoint]->GetDomain()) {
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Coord = geometry->node[iPoint]->GetCoord();
			
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {
				
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				
				/*--- Create a copy of the adjoint solution ---*/
				Psi_Aux = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];			
				
				/*--- Flow solution and projected force vector ---*/
				U = solver_container[PLASMA_SOL]->node[iPoint]->GetSolution();
				d = node[iPoint]->GetForceProj_Vector();
				
        /*--- Geometry parameters ---*/
				Area = 0; 
				for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
				Area = sqrt(Area);
				for (iDim = 0; iDim < nDim; iDim++)
					UnitaryNormal[iDim]   = -Normal[iDim]/Area;
								

        Gamma_Minus_One = config->GetSpecies_Gamma(iSpecies);
				Enthalpy = solver_container[PLASMA_SOL]->node[iPoint]->GetEnthalpy(iSpecies);
        hf = config->GetEnthalpy_Formation(iSpecies);
        sq_vel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
					Velocity[iDim] = U[loc+iDim+1] / U[loc+0];
          sq_vel += Velocity[iDim]*Velocity[iDim];
        }
        Energy_el = 0;
        Energy_vib = 0.0;
        if (iSpecies < nDiatomics)
          Energy_vib = U[loc+nDim+2]/U[loc+0];
				
				/*--- Compute projections ---*/
				ProjVel = 0.0; bcn = 0.0; phin = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
//					ProjVel -= Velocity[iDim]*Normal[iDim];
					ProjVel += Velocity[iDim]*UnitaryNormal[iDim];
					bcn     += d[iDim]*UnitaryNormal[iDim];
					phin    += Psi[loc+iDim+1]*UnitaryNormal[iDim];
				}

				/*--- Introduce the boundary condition ---*/
				for (iDim = 0; iDim < nDim; iDim++) 
					Psi[loc+iDim+1] -= ( phin - bcn ) * UnitaryNormal[iDim];
        
        /*--- Pre-compute useful quantities ---*/
        phidotu = 0.0;
        phidotn = 0.0;
        for (iDim = 0; iDim < nDim; iDim ++) {
          phidotu += Velocity[iDim] * Psi[loc+iDim+1];
          phidotn += Psi[loc+iDim+1] * UnitaryNormal[iDim];
        }
        dPdrho = Gamma_Minus_One * (0.5*sq_vel - hf - Energy_el);
        for (iDim = 0; iDim < nDim; iDim++)
          dPdrhou[iDim] = -Gamma_Minus_One*Velocity[iDim];
        dPdrhoE = Gamma_Minus_One;
        dPdrhoEv = -Gamma_Minus_One;
        
        /*--- Flux of the Euler wall: Psi^T * (dF/dU dot n) ---*/
        Residual[loc+0] = dPdrho*phidotn - ProjVel*phidotu + ProjVel*(dPdrho-Enthalpy)*Psi[loc+nDim+1];
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[loc+iDim+1] = UnitaryNormal[iDim]*Psi[loc] + ProjVel*Psi[loc+iDim+1] + phidotu*UnitaryNormal[iDim] + dPdrhou[iDim]*phidotn 
          + (dPdrhou[iDim]*ProjVel+Enthalpy*UnitaryNormal[iDim])*Psi[loc+nDim+1];
        Residual[loc+nDim+1] = dPdrhoE*phidotn + ProjVel*(1.0+dPdrhoE)*Psi[loc+nDim+1];
        
        if (iSpecies < nDiatomics) {
          Residual[loc+0] -= ProjVel*Energy_vib*Psi[loc+nDim+2];
          for (iDim = 0; iDim < nDim; iDim++)
            Residual[loc+iDim+1] += Energy_vib*UnitaryNormal[iDim]*Psi[loc+nDim+2];
          Residual[loc+nDim+1] += 0.0;
          Residual[loc+nDim+2] = dPdrhoEv*phidotn + ProjVel*dPdrhoEv*Psi[loc+nDim+1] + ProjVel*Psi[loc+nDim+2];
        }       
        
        /*--- Calculate Jacobians for implicit time marching ---*/
        if (implicit) {
          
          /*--- Adjoint density ---*/
          Jacobian_ii[loc+0][loc+0] = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_ii[loc+0][loc+iDim+1] = dPdrho*UnitaryNormal[iDim] - ProjVel*Velocity[iDim];
          Jacobian_ii[loc+0][loc+nDim+1] = ProjVel*(dPdrho - Enthalpy);
          
          /*--- Adjoint velocity ---*/
          for (iDim = 0; iDim < nDim; iDim++) {
            Jacobian_ii[loc+iDim+1][0] = UnitaryNormal[iDim];
            for (jDim = 0; jDim < nDim; jDim++)
              Jacobian_ii[loc+iDim+1][loc+jDim+1] = Velocity[jDim]*UnitaryNormal[iDim] + dPdrhou[iDim]*UnitaryNormal[jDim];
            Jacobian_ii[loc+iDim+1][loc+iDim+1] += ProjVel;
            Jacobian_ii[loc+iDim+1][loc+nDim+1] = dPdrhou[iDim]*ProjVel + Enthalpy*UnitaryNormal[iDim];
          }
          
          /*--- Adjoint energy ---*/
          Jacobian_ii[loc+nDim+1][loc+0] = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) 
            Jacobian_ii[loc+nDim+1][loc+iDim+1] = dPdrhoE*UnitaryNormal[iDim];
          Jacobian_ii[loc+nDim+1][loc+nDim+1] = ProjVel*(1.0+dPdrhoE);
          
          /*--- Adjoint vibrational energy ---*/
          if (iSpecies < nDiatomics) {
            Jacobian_ii[loc+0][loc+nDim+2] = -ProjVel*Energy_vib;
            for (iDim = 0; iDim < nDim; iDim++)
              Jacobian_ii[loc+iDim+1][loc+nDim+2] = Energy_vib*UnitaryNormal[iDim];
            Jacobian_ii[loc+nDim+1][loc+nDim+2] = 0.0;
            
            Jacobian_ii[loc+nDim+2][loc+0] = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Jacobian_ii[loc+nDim+2][loc+iDim+1] = dPdrhoEv*UnitaryNormal[iDim];
            Jacobian_ii[loc+nDim+2][loc+nDim+1] = ProjVel*dPdrhoEv;
            Jacobian_ii[loc+nDim+2][loc+nDim+2] = ProjVel;
          }  
        }
              
        /*--- Integrate over the area --*/
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] = Residual[iVar]*Area;
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_ii[iVar][jVar] = Jacobian_ii[iVar][jVar]*Area;          
        }        
			}			
			/*--- Update residual ---*/
			LinSysRes.SubtractBlock(iPoint, Residual);
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}	
	delete [] Velocity;
	delete [] UnitaryNormal;
	delete [] Psi;
  delete [] dPdrhou;
}

void CAdjPlasmaSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
																		 CConfig *config, unsigned short val_marker) {
	
	unsigned long iVertex, iPoint;
	double *Normal, *U, *Psi_Aux, ProjVel, bcn, Area, *UnitaryNormal, *Coord, Gamma_Minus_One;
  double *Velocity, *Psi, Enthalpy, hf, Energy_vib, sq_vel, phin;
	unsigned short iDim, iVar, jVar, jDim, loc, iSpecies;
  double phidotu, phidotn, Energy_el, dPdrho, *dPdrhou, dPdrhoE, dPdrhoEv;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	
	UnitaryNormal = new double[nDim];
	Velocity = new double[nDim];
	Psi      = new double[nVar];
  dPdrhou = new double[nDim];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		if (geometry->node[iPoint]->GetDomain()) {
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Coord = geometry->node[iPoint]->GetCoord();
			
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {
				
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				
				/*--- Create a copy of the adjoint solution ---*/
				Psi_Aux = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];
				
				/*--- Flow solution and projected force vector ---*/
				U = solver_container[PLASMA_SOL]->node[iPoint]->GetSolution();
				
        /*--- Geometry parameters ---*/
				Area = 0;
				for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
				Area = sqrt(Area);
				for (iDim = 0; iDim < nDim; iDim++)
					UnitaryNormal[iDim]   = -Normal[iDim]/Area;
        
        Gamma_Minus_One = config->GetSpecies_Gamma(iSpecies);
				Enthalpy = solver_container[PLASMA_SOL]->node[iPoint]->GetEnthalpy(iSpecies);
        hf = config->GetEnthalpy_Formation(iSpecies);
        sq_vel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
					Velocity[iDim] = U[loc+iDim+1] / U[loc+0];
          sq_vel += Velocity[iDim]*Velocity[iDim];
        }
        Energy_el = 0;
        Energy_vib = 0.0;
        if (iSpecies < nDiatomics)
          Energy_vib = U[loc+nDim+2]/U[loc+0];
				
				/*--- Compute projections ---*/
				ProjVel = 0.0; bcn = 0.0; phin = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
          //					ProjVel -= Velocity[iDim]*Normal[iDim];
					ProjVel += Velocity[iDim]*UnitaryNormal[iDim];
					phin    += Psi[loc+iDim+1]*UnitaryNormal[iDim];
				}
        
				/*--- Introduce the boundary condition ---*/
				for (iDim = 0; iDim < nDim; iDim++)
					Psi[loc+iDim+1] -= phin * UnitaryNormal[iDim];
              
        /*--- Pre-compute useful quantities ---*/
        phidotu = 0.0;
        phidotn = 0.0;
        for (iDim = 0; iDim < nDim; iDim ++) {
          phidotu += Velocity[iDim] * Psi[loc+iDim+1];
          phidotn += Psi[loc+iDim+1] * UnitaryNormal[iDim];
        }
        dPdrho = Gamma_Minus_One * (0.5*sq_vel - hf - Energy_el);
        for (iDim = 0; iDim < nDim; iDim++)
          dPdrhou[iDim] = -Gamma_Minus_One*Velocity[iDim];
        dPdrhoE = Gamma_Minus_One;
        dPdrhoEv = -Gamma_Minus_One;
        
        /*--- Flux of the Euler wall: Psi^T * (dF/dU dot n) ---*/
        Residual[loc+0] = dPdrho*phidotn - ProjVel*phidotu + ProjVel*(dPdrho-Enthalpy)*Psi[loc+nDim+1];
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[loc+iDim+1] = UnitaryNormal[iDim]*Psi[loc] + ProjVel*Psi[loc+iDim+1] + phidotu*UnitaryNormal[iDim] + dPdrhou[iDim]*phidotn
          + (dPdrhou[iDim]*ProjVel+Enthalpy*UnitaryNormal[iDim])*Psi[loc+nDim+1];
        Residual[loc+nDim+1] = dPdrhoE*phidotn + ProjVel*(1.0+dPdrhoE)*Psi[loc+nDim+1];
        
        if (iSpecies < nDiatomics) {
          Residual[loc+0] -= ProjVel*Energy_vib*Psi[loc+nDim+2];
          for (iDim = 0; iDim < nDim; iDim++)
            Residual[loc+iDim+1] += Energy_vib*UnitaryNormal[iDim]*Psi[loc+nDim+2];
          Residual[loc+nDim+1] += 0.0;
          Residual[loc+nDim+2] = dPdrhoEv*phidotn + ProjVel*dPdrhoEv*Psi[loc+nDim+1] + ProjVel*Psi[loc+nDim+2];
        }
        
        /*--- Calculate Jacobians for implicit time marching ---*/
        if (implicit) {
          
          /*--- Adjoint density ---*/
          Jacobian_ii[loc+0][loc+0] = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_ii[loc+0][loc+iDim+1] = dPdrho*UnitaryNormal[iDim] - ProjVel*Velocity[iDim];
          Jacobian_ii[loc+0][loc+nDim+1] = ProjVel*(dPdrho - Enthalpy);
          
          /*--- Adjoint velocity ---*/
          for (iDim = 0; iDim < nDim; iDim++) {
            Jacobian_ii[loc+iDim+1][0] = UnitaryNormal[iDim];
            for (jDim = 0; jDim < nDim; jDim++)
              Jacobian_ii[loc+iDim+1][loc+jDim+1] = Velocity[jDim]*UnitaryNormal[iDim] + dPdrhou[iDim]*UnitaryNormal[jDim];
            Jacobian_ii[loc+iDim+1][loc+iDim+1] += ProjVel;
            Jacobian_ii[loc+iDim+1][loc+nDim+1] = dPdrhou[iDim]*ProjVel + Enthalpy*UnitaryNormal[iDim];
          }
          
          /*--- Adjoint energy ---*/
          Jacobian_ii[loc+nDim+1][loc+0] = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_ii[loc+nDim+1][loc+iDim+1] = dPdrhoE*UnitaryNormal[iDim];
          Jacobian_ii[loc+nDim+1][loc+nDim+1] = ProjVel*(1.0+dPdrhoE);
          
          /*--- Adjoint vibrational energy ---*/
          if (iSpecies < nDiatomics) {
            Jacobian_ii[loc+0][loc+nDim+2] = -ProjVel*Energy_vib;
            for (iDim = 0; iDim < nDim; iDim++)
              Jacobian_ii[loc+iDim+1][loc+nDim+2] = Energy_vib*UnitaryNormal[iDim];
            Jacobian_ii[loc+nDim+1][loc+nDim+2] = 0.0;
            
            Jacobian_ii[loc+nDim+2][loc+0] = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Jacobian_ii[loc+nDim+2][loc+iDim+1] = dPdrhoEv*UnitaryNormal[iDim];
            Jacobian_ii[loc+nDim+2][loc+nDim+1] = ProjVel*dPdrhoEv;
            Jacobian_ii[loc+nDim+2][loc+nDim+2] = ProjVel;
          }
        }
        
        /*--- Integrate over the area --*/
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] = Residual[iVar]*Area;
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_ii[iVar][jVar] = Jacobian_ii[iVar][jVar]*Area;
        }
			}
			/*--- Update residual ---*/
			LinSysRes.SubtractBlock(iPoint, Residual);
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}
	delete [] Velocity;
	delete [] UnitaryNormal;
	delete [] Psi;
  delete [] dPdrhou;
}

void CAdjPlasmaSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	double *Normal, *U_domain, *U_infty, *Psi_domain, *Psi_infty;
	unsigned short loc, iSpecies;

	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	Normal = new double[nDim];
	U_domain = new double[nVar]; U_infty = new double[nVar];
	Psi_domain = new double[nVar]; Psi_infty = new double[nVar];
		
	/*--- Loop over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
			conv_numerics->SetNormal(Normal);
			
			/*--- Flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = solver_container[PLASMA_SOL]->node[iPoint]->GetSolution(iVar);
			
			/*--- Solution at the infinity ---*/
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				U_infty[loc + 0] = solver_container[PLASMA_SOL]->GetDensity_Inf(iSpecies);
        for (iDim = 0; iDim < nDim; iDim++)
          U_infty[loc+iDim+1] = solver_container[PLASMA_SOL]->GetDensity_Velocity_Inf(iDim, iSpecies);
				U_infty[loc+nDim+1] = solver_container[PLASMA_SOL]->GetDensity_Energy_Inf(iSpecies);
        if (iSpecies < nDiatomics)
          U_infty[loc+nDim+2] = solver_container[PLASMA_SOL]->GetDensity_Energy_vib_Inf(iSpecies);          
			}      
			
			conv_numerics->SetConservative(U_domain, U_infty);

			/*--- Adjoint flow solution at the farfield ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
				Psi_infty[iVar] = 0.0;
			}
			conv_numerics->SetAdjointVar(Psi_domain, Psi_infty);
			
			/*--- Compute the upwind flux ---*/
			conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			/*--- Add and Subtract Residual ---*/
			LinSysRes.SubtractBlock(iPoint, Residual_i);
			
			/*--- Implicit contribution to the residual ---*/
			if (implicit) 
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}
	
	delete [] Normal;
	delete [] U_domain; delete [] U_infty;
	delete [] Psi_domain; delete [] Psi_infty;
}
