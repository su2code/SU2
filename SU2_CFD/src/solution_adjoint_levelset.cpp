/*!
 * \file solution_adjoint_levelset.cpp
 * \brief Main subrotuines for solving the level set problem.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.5
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

#include "../include/solution_structure.hpp"

CAdjLevelSetSolution::CAdjLevelSetSolution(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolution() {
	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double dull_val;
	ifstream restart_file;
	string text_line, mesh_filename, filename, AdjExt;
	
  bool restart = config->GetRestart();
	
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	node = new CVariable*[geometry->GetnPoint()];
	
	/*--- Dimension of the problem ---*/
	nVar = 1;
  
  /*--- Single grid simulation ---*/
	if (iMesh == MESH_0) {
	
    /*--- Define some auxiliar vector related with the residual ---*/
    Residual      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
    Residual_RMS  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
    Residual_Max  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
    Point_Max  = new unsigned long[nVar]; for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
    Residual_i    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
    Residual_j    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
	
    /*--- Define some auxiliar vector related with the solution ---*/
    Solution    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]    = 0.0;
    Solution_i  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar]  = 0.0;
    Solution_j  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar]  = 0.0;
    
    /*--- Define some auxiliar vector related with the geometry ---*/
    Vector    = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]    = 0.0;
    Vector_i  = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim]  = 0.0;
    Vector_j  = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim]  = 0.0;
	
    /*--- Define some auxiliar vector related with the flow solution ---*/
    FlowSolution_i = new double[nDim+1]; for (iVar = 0; iVar < nDim+1; iVar++) FlowSolution_i[iVar]  = 0.0;
    FlowSolution_j = new double[nDim+1]; for (iVar = 0; iVar < nDim+1; iVar++) FlowSolution_j[iVar]  = 0.0;
	
    /*--- Solution and residual vectors ---*/
    xsol = new double [geometry->GetnPoint()*nVar];
    xres = new double [geometry->GetnPoint()*nVar];
    
    /*--- Jacobians and vector structures for implicit computations ---*/
    if (config->GetKind_TimeIntScheme_AdjLevelSet() == EULER_IMPLICIT) {
      
      /*--- Point to point Jacobians ---*/
      Jacobian_i = new double* [nVar];
      Jacobian_j = new double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[iVar] = new double [nVar];
        Jacobian_j[iVar] = new double [nVar];
      }
      
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
      
      /*--- Initialization of the structure of the whole Jacobian ---*/
      if (rank == MASTER_NODE) cout << "Initialize jacobian structure (Adj. Level Set). MG level: 0." << endl;
      Initialize_SparseMatrix_Structure(&Jacobian, nVar, nVar, geometry, config);
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
    
  }
	
	/*--- Restart the solution from file information ---*/
	if (!restart || geometry->GetFinestMGLevel() == false) {
    
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			node[iPoint] = new CAdjLevelSetVariable(0.0, nDim, nVar, config);
		}
    
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
		
		if (restart_file.fail()) {
			cout << "There is no level set restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
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
    for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
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
        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        node[iPoint_Local] = new CAdjLevelSetVariable(Solution[0], nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
      node[iPoint] = new CAdjLevelSetVariable(Solution[0], nDim, nVar, config);
    }
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
	}
  
  /*--- MPI solution ---*/
  SetSolution_MPI(geometry, config);
  
  /*--- Compute the gradient (Mean flow source term) ---*/
  if (iMesh == MESH_0) {
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  }
  
}

CAdjLevelSetSolution::~CAdjLevelSetSolution(void) { }

void CAdjLevelSetSolution::SetSolution_MPI(CGeometry *geometry, CConfig *config) {
	unsigned short iMarker;
	unsigned long iVertex, iPoint, nVertex, nBuffer_Vector;
	double *Buffer_Receive_U = NULL;
	short SendRecv;
	int send_to, receive_from;
  
#ifndef NO_MPI
	double *Buffer_Send_U = NULL;
#endif
  
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
      
			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			nVertex = geometry->nVertex[iMarker];
			nBuffer_Vector = nVertex*nVar;
      
			send_to = SendRecv-1;
			receive_from = abs(SendRecv)-1;
      
#ifndef NO_MPI
      
			/*--- Send information using MPI  ---*/
			if (SendRecv > 0) {
        Buffer_Send_U = new double[nBuffer_Vector];
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Buffer_Send_U[iVertex] = node[iPoint]->GetSolution(0);
				}
        MPI::COMM_WORLD.Bsend(Buffer_Send_U, nBuffer_Vector, MPI::DOUBLE, send_to, 0); delete [] Buffer_Send_U;
			}
      
#endif
      
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
        Buffer_Receive_U = new double [nBuffer_Vector];
        
#ifdef NO_MPI
				/*--- Get the information from the donor point directly. This is a
				 serial computation with access to all nodes. Note that there is an
				 implicit ordering in the list. ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Buffer_Receive_U[iVertex] = node[iPoint]->GetSolution(0);
        }
#else
        MPI::COMM_WORLD.Recv(Buffer_Receive_U, nBuffer_Vector, MPI::DOUBLE, receive_from, 0);
#endif
        
				/*--- Do the coordinate transformation ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
          
					/*--- Find point and its type of transformation ---*/
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          node[iPoint]->SetSolution(0, Buffer_Receive_U[iVertex]);
          
				}
        delete [] Buffer_Receive_U;
			}
		}
	}
}

void CAdjLevelSetSolution::SetSolution_Limiter_MPI(CGeometry *geometry, CConfig *config) {
	unsigned short iMarker;
	unsigned long iVertex, iPoint, nVertex, nBuffer_Vector;
	double *Buffer_Receive_Limit = NULL;
	short SendRecv;
	int send_to, receive_from;
  
#ifndef NO_MPI
	double *Buffer_Send_Limit;
#endif
  
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
				Buffer_Send_Limit = new double[nBuffer_Vector];
        
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Buffer_Send_Limit[iVertex] = node[iPoint]->GetLimiter(0);
          
				}
				/*--- Send the buffer, and deallocate information using MPI ---*/
        MPI::COMM_WORLD.Bsend(Buffer_Send_Limit, nBuffer_Vector, MPI::DOUBLE, send_to, 2); delete [] Buffer_Send_Limit;
			}
      
#endif
      
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
        
				/*--- Allocate upwind variables ---*/
				Buffer_Receive_Limit = new double [nBuffer_Vector];
        
#ifdef NO_MPI
        
				/*--- Get the information from the donor point directly. This is a
				 serial computation with access to all nodes. Note that there is an
				 implicit ordering in the list. ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Buffer_Receive_Limit[iVertex] = node[iPoint]->GetLimiter(0);
				}
        
#else
				/*--- Receive the information using MPI---*/
        MPI::COMM_WORLD.Recv(Buffer_Receive_Limit, nBuffer_Vector, MPI::DOUBLE, receive_from, 2);
#endif
        
				/*--- Do the coordinate transformation ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          node[iPoint]->SetLimiter(0, Buffer_Receive_Limit[iVertex]);
				}
				delete [] Buffer_Receive_Limit;
			}
		}
	}
  
}

void CAdjLevelSetSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_AdjLevelSet() == EULER_IMPLICIT);
	bool high_order_diss = (config->GetKind_Upwind_AdjLevelSet() == SCALAR_UPWIND_2ND);

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		Set_Residual_Zero(iPoint);
	
	/*--- Implicit part ---*/
	if (implicit)
		Jacobian.SetValZero();
	
  if (high_order_diss) {
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  }
  
}

void CAdjLevelSetSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {
	double *LevelSet_var_i, *LevelSet_var_j, *U_i, *U_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, DensityInc_i, DensityInc_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	
	bool implicit = (config->GetKind_TimeIntScheme_AdjLevelSet() == EULER_IMPLICIT);
	bool high_order_diss = (config->GetKind_Upwind_AdjLevelSet() == SCALAR_UPWIND_2ND);
		
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());
		
		/*--- Conservative variables w/o reconstruction ---*/
		U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
		U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
		solver->SetConservative(U_i, U_j);
		
		/*--- Level Set variables w/o reconstruction ---*/
		LevelSet_var_i = node[iPoint]->GetSolution();
		LevelSet_var_j = node[jPoint]->GetSolution();
		solver->SetLevelSetVar(LevelSet_var_i, LevelSet_var_j);
		
		/*--- Set the value of the density ---*/
		DensityInc_i = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
		DensityInc_j = solution_container[FLOW_SOL]->node[jPoint]->GetDensityInc();
		solver->SetDensityInc(DensityInc_i, DensityInc_j);
				
		if (high_order_diss) {

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
			solver->SetLevelSetVar(Solution_i, Solution_j);
      
		}
		
		solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);	
		
		/*--- Add and subtract Residual ---*/
		AddResidual(iPoint, Residual_i);
		AddResidual(jPoint, Residual_j);
    		
		/*--- Implicit part ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
			Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
			Jacobian.AddBlock(jPoint, iPoint, Jacobian_ji);
			Jacobian.AddBlock(jPoint, jPoint, Jacobian_jj);
		}		
		
	}
}

void CAdjLevelSetSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CNumerics *second_solver,
																								 CConfig *config, unsigned short iMesh) {
	unsigned short iVar, iDim;
	unsigned long iPoint;
  double epsilon, DeltaDirac, lambda, dRho_dPhi, dMud_Phi, Vol, DiffLevelSet, LevelSet, *MeanFlow, *AdjMeanFlow, **AdjLevelSetGradient, **AdjMeanFlowGradient, Density, Velocity[3], ProjAdj, dFc_dRho[3][4], ProjFlux;
  
	double Froude2 = config->GetFroude()*config->GetFroude();
  bool incompressible = config->GetIncompressible();

	for (iVar = 0; iVar < nVar; iVar++)
		Residual[iVar] = 0;

	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 
		
		Vol = geometry->node[iPoint]->GetVolume();

    /*--- Direct problem quantities ---*/
		DiffLevelSet = solution_container[LEVELSET_SOL]->node[iPoint]->GetDiffLevelSet();
		LevelSet = solution_container[LEVELSET_SOL]->node[iPoint]->GetSolution(0);
		MeanFlow = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
    Density = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
    
    /*--- Adjoint problem quantities ---*/
		AdjMeanFlow = solution_container[ADJFLOW_SOL]->node[iPoint]->GetSolution();
    AdjLevelSetGradient = solution_container[ADJLEVELSET_SOL]->node[iPoint]->GetGradient();
    AdjMeanFlowGradient = solution_container[ADJFLOW_SOL]->node[iPoint]->GetGradient();
    
    /*--- Projected adjoint velocity ---*/
    ProjAdj = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity[iDim] = solution_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim, incompressible);
      ProjAdj += Velocity[iDim]*AdjLevelSetGradient[0][iDim];
    }
    
    /*--- Compute the flow solution using the level set value. ---*/
    epsilon = config->GetFreeSurface_Thickness();
    DeltaDirac = 0.0;
    if (fabs(LevelSet) <= epsilon) DeltaDirac = 1.0 - (0.5*(1.0+(LevelSet/epsilon)+(1.0/PI_NUMBER)*sin(PI_NUMBER*LevelSet/epsilon)));
    
    /*--- Set the value of the incompressible density for free surface flows (density ratio g/l) ---*/
    lambda = config->GetRatioDensity();
    dRho_dPhi = (1.0 - lambda)*DeltaDirac*config->GetDensity_FreeStreamND();
    
    /*--- Set the value of the incompressible viscosity for free surface flows (viscosity ratio g/l) ---*/
    lambda = config->GetRatioViscosity();
    dMud_Phi = (1.0 - lambda)*DeltaDirac*config->GetViscosity_FreeStreamND();

    /*--- Flux derivative ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      dFc_dRho[iDim][0] = 0.0;
      dFc_dRho[iDim][1] = Velocity[iDim]*Velocity[0];
      dFc_dRho[iDim][2] = Velocity[iDim]*Velocity[1];
      dFc_dRho[iDim][3] = Velocity[iDim]*Velocity[2];
    }
    
    /*--- Projected flux derivative ---*/
    ProjFlux = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
        ProjFlux += AdjMeanFlowGradient[iVar][iDim]*dFc_dRho[iDim][iVar];
      }
    }
    
 		Residual[0] = ( ProjFlux -(LevelSet*ProjAdj/Density) - (AdjMeanFlow[nDim]/Froude2))*dRho_dPhi*Vol;
		
    /*--- The Free surface objective function requires the levelt set difference ---*/
    if (config->GetKind_ObjFunc() == FREE_SURFACE) Residual[0] +=  DiffLevelSet* Vol;

		/*--- Add Residual ---*/
		AddResidual(iPoint, Residual);
	}
	
}

void CAdjLevelSetSolution::Source_Template(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) { }

void CAdjLevelSetSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	double *U_domain, *U_wall, *LevelSet_domain, *LevelSet_wall;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new double[solution_container[FLOW_SOL]->GetnVar()]; 
	U_wall = new double[solution_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new double[1];
	LevelSet_wall = new double[1];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
						
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
				U_domain[iVar] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(iVar);
				U_wall[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			}
			
			LevelSet_domain[0] = node[Point_Normal]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			solver->SetNormal(Vector);
			
			solver->SetConservative(U_wall, U_wall);
			solver->SetLevelSetVar(LevelSet_wall, LevelSet_wall);
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);	
			
			AddResidual(iPoint, Residual_i);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
			
		}
	}
}

void CAdjLevelSetSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver,
                                        CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	double *U_domain, *U_wall, *LevelSet_domain, *LevelSet_wall;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new double[solution_container[FLOW_SOL]->GetnVar()];
	U_wall = new double[solution_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new double[1];
	LevelSet_wall = new double[1];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
				U_domain[iVar] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(iVar);
				U_wall[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			}
			
			LevelSet_domain[0] = node[Point_Normal]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			conv_solver->SetNormal(Vector);
			
			conv_solver->SetConservative(U_wall, U_wall);
			conv_solver->SetLevelSetVar(LevelSet_wall, LevelSet_wall);
			conv_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      conv_solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			AddResidual(iPoint, Residual_i);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
			
		}
	}
}

void CAdjLevelSetSolution::BC_HeatFlux_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	double *U_domain, *U_wall, *LevelSet_domain, *LevelSet_wall;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new double[solution_container[FLOW_SOL]->GetnVar()];
	U_wall = new double[solution_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new double[1];
	LevelSet_wall = new double[1];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
				U_domain[iVar] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(iVar);
				U_wall[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			}
			
			LevelSet_domain[0] = node[Point_Normal]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			conv_solver->SetNormal(Vector);
			
			conv_solver->SetConservative(U_wall, U_wall);
			conv_solver->SetLevelSetVar(LevelSet_wall, LevelSet_wall);
			conv_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      conv_solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			AddResidual(iPoint, Residual_i);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
			
		}
	}
}

void CAdjLevelSetSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	double *U_domain, *U_wall, *LevelSet_domain, *LevelSet_wall;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new double[solution_container[FLOW_SOL]->GetnVar()];
	U_wall = new double[solution_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new double[1];
	LevelSet_wall = new double[1];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
				U_domain[iVar] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(iVar);
				U_wall[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			}
			
			LevelSet_domain[0] = node[Point_Normal]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			conv_solver->SetNormal(Vector);
			
			conv_solver->SetConservative(U_wall, U_wall);
			conv_solver->SetLevelSetVar(LevelSet_wall, LevelSet_wall);
			conv_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      conv_solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			AddResidual(iPoint, Residual_i);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
			
		}
	}
}

void CAdjLevelSetSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	double *U_domain, *U_wall, *LevelSet_domain, *LevelSet_wall;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new double[solution_container[FLOW_SOL]->GetnVar()];
	U_wall = new double[solution_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new double[1];
	LevelSet_wall = new double[1];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
				U_domain[iVar] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(iVar);
				U_wall[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			}
			
			LevelSet_domain[0] = node[Point_Normal]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			conv_solver->SetNormal(Vector);
			
			conv_solver->SetConservative(U_wall, U_wall);
			conv_solver->SetLevelSetVar(LevelSet_wall, LevelSet_wall);
			conv_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      conv_solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			AddResidual(iPoint, Residual_i);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
			
		}
	}
}

void CAdjLevelSetSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	
	bool implicit = (config->GetKind_TimeIntScheme_AdjLevelSet() == EULER_IMPLICIT);
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set the solution to the original value ---*/
		Solution[0] = node[iPoint]->GetSolution(0);
		
		node[iPoint]->SetSolution_Old(Solution);
		Set_Residual_Zero(iPoint);
		
		/*--- Includes 1 in the diagonal ---*/
		if (implicit)
			Jacobian.DeleteValsRowi(iPoint);
	}
}

void CAdjLevelSetSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned long iPoint;
	double Delta = 0.0, Vol;
	
	/*--- Set maximum residual to zero ---*/
  SetRes_RMS(0, 0.0);
  SetRes_Max(0, 0.0, 0);
  
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		/*--- Read the volume ---*/
		Vol = geometry->node[iPoint]->GetVolume();
		
		/*--- Modify matrix diagonal to assure diagonal dominance ---*/
		Delta = Vol / solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time();
    
		Jacobian.AddVal2Diag(iPoint,Delta);
    
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/    
		xres[iPoint] = -xres[iPoint];
		xsol[iPoint] = 0.0;
		AddRes_RMS(0, xres[iPoint]*xres[iPoint]);
    AddRes_Max(0, fabs(xres[iPoint]), geometry->node[iPoint]->GetGlobalIndex());
	}
  
  /*--- Initialize residual and solution at the ghost points ---*/
  for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
    xres[iPoint] = 0.0;
    xsol[iPoint] = 0.0;
  }
	
  /*--- Solve the linear system (Stationary iterative methods) ---*/
	if (config->GetKind_Linear_Solver() == SYM_GAUSS_SEIDEL)
		Jacobian.SGSSolution(xres, xsol, config->GetLinear_Solver_Error(),
												 config->GetLinear_Solver_Iter(), false, geometry, config);
	
	if (config->GetKind_Linear_Solver() == LU_SGS)
    Jacobian.LU_SGSIteration(xres, xsol, geometry, config);
	
	/*--- Solve the linear system (Krylov subspace methods) ---*/
	if ((config->GetKind_Linear_Solver() == BCGSTAB) ||
			(config->GetKind_Linear_Solver() == GMRES)) {
		
		CSysVector rhs_vec((const unsigned int)geometry->GetnPoint(),
                       (const unsigned int)geometry->GetnPointDomain(), nVar, xres);
		CSysVector sol_vec((const unsigned int)geometry->GetnPoint(),
                       (const unsigned int)geometry->GetnPointDomain(), nVar, xsol);
		
		CMatrixVectorProduct* mat_vec = new CSparseMatrixVectorProduct(Jacobian);
		CSolutionSendReceive* sol_mpi = new CSparseMatrixSolMPI(Jacobian, geometry, config);
    
		CPreconditioner* precond = NULL;
		if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
			Jacobian.BuildJacobiPreconditioner();
			precond = new CJacobiPreconditioner(Jacobian);
		}
		else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
			Jacobian.BuildJacobiPreconditioner();
			precond = new CLineletPreconditioner(Jacobian);
		}
		else if (config->GetKind_Linear_Solver_Prec() == NO_PREC)
			precond = new CIdentityPreconditioner();
		
		CSysSolve system;
		if (config->GetKind_Linear_Solver() == BCGSTAB)
			system.BCGSTAB(rhs_vec, sol_vec, *mat_vec, *precond, *sol_mpi, config->GetLinear_Solver_Error(),
										 config->GetLinear_Solver_Iter(), false);
		else if (config->GetKind_Linear_Solver() == GMRES)
			system.GMRES(rhs_vec, sol_vec, *mat_vec, *precond, *sol_mpi, config->GetLinear_Solver_Error(),
                   config->GetLinear_Solver_Iter(), false);
		
		sol_vec.CopyToArray(xsol);
		delete mat_vec;
		delete precond;
		delete sol_mpi;
	}
	
	/*--- Update solution (system written in terms of increments), be careful with the update of the
	 scalar equations which includes the density ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->AddSolution(0, config->GetLinear_Solver_Relax()*xsol[iPoint]);
  
  /*--- MPI solution ---*/
  SetSolution_MPI(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CAdjLevelSetSolution::SetResidual_DualTime(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1, Volume_nM1, Volume_n, Volume_nP1, TimeStep;
	
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
		AddResidual(iPoint, Residual);
		
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
