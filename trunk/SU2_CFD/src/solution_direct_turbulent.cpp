/*!
 * \file solution_direct_turbulent.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.1
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

CTurbSolution::CTurbSolution(void) : CSolution() { }

CTurbSolution::CTurbSolution(CConfig *config) : CSolution() {
	
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
}
CTurbSolution::~CTurbSolution(void) {}

void CTurbSolution::MPI_Send_Receive(CGeometry ***geometry, CSolution ****solution_container,
                                     CConfig **config, unsigned short iMGLevel, unsigned short iZone) {

#ifndef NO_MPI
	unsigned short iVar, iMarker;
	double *Turb_Var, **Turb_Grad;
	unsigned long iVertex, iPoint;
	
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config[iZone]->GetnMarker_All(); iMarker++)
		if (config[iZone]->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			
			short SendRecv = config[iZone]->GetMarker_All_SendRecv(iMarker);
			unsigned long nVertex = geometry[iZone][iMGLevel]->nVertex[iMarker];
			
			/*--- Send information  ---*/
			if (SendRecv > 0) {
				/*--- Dimensionalization ---*/
				unsigned long nBuffer_Vector = geometry[iZone][iMGLevel]->nVertex[iMarker]*nVar;
				int send_to = SendRecv-1;
				
				double *Buffer_Send_Turb = new double[nBuffer_Vector];
				double *Buffer_Send_Turbx = new double[nBuffer_Vector];
				double *Buffer_Send_Turby = new double[nBuffer_Vector];
				double *Buffer_Send_Turbz = new double[nBuffer_Vector];
				
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();

					Turb_Var = node[iPoint]->GetSolution();
					Turb_Grad = node[iPoint]->GetGradient();
					
					for (iVar = 0; iVar < nVar; iVar++) {
						Buffer_Send_Turb[iVar*nVertex+iVertex] = Turb_Var[iVar];
						Buffer_Send_Turbx[iVar*nVertex+iVertex] = Turb_Grad[iVar][0];
						Buffer_Send_Turby[iVar*nVertex+iVertex] = Turb_Grad[iVar][1];
						if (nDim == 3) Buffer_Send_Turbz[iVar*nVertex+iVertex] = Turb_Grad[iVar][2];
					}
				}
				
				MPI::COMM_WORLD.Bsend(Buffer_Send_Turb, nBuffer_Vector, MPI::DOUBLE, send_to, 0);
				MPI::COMM_WORLD.Bsend(Buffer_Send_Turbx, nBuffer_Vector, MPI::DOUBLE, send_to, 1);
				MPI::COMM_WORLD.Bsend(Buffer_Send_Turby, nBuffer_Vector, MPI::DOUBLE, send_to, 2);
				if (nDim == 3) MPI::COMM_WORLD.Bsend(Buffer_Send_Turbz, nBuffer_Vector, MPI::DOUBLE, send_to, 3);
				
				delete [] Buffer_Send_Turb;
				delete [] Buffer_Send_Turbx;
				delete [] Buffer_Send_Turby;
				delete [] Buffer_Send_Turbz;
				
			}
			
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
				/*--- Dimensionalization ---*/
				unsigned long nBuffer_Vector = geometry[iZone][iMGLevel]->nVertex[iMarker]*nVar;
				int receive_from = abs(SendRecv)-1;
				
				double *Buffer_Receive_Turb = new double [nBuffer_Vector];
				double *Buffer_Receive_Turbx = new double [nBuffer_Vector];
				double *Buffer_Receive_Turby = new double [nBuffer_Vector];
				double *Buffer_Receive_Turbz = new double [nBuffer_Vector];
				
				MPI::COMM_WORLD.Recv(Buffer_Receive_Turb, nBuffer_Vector, MPI::DOUBLE, receive_from, 0);
				MPI::COMM_WORLD.Recv(Buffer_Receive_Turbx, nBuffer_Vector, MPI::DOUBLE, receive_from, 1);
				MPI::COMM_WORLD.Recv(Buffer_Receive_Turby, nBuffer_Vector, MPI::DOUBLE, receive_from, 2);
				if (nDim == 3) MPI::COMM_WORLD.Recv(Buffer_Receive_Turbz, nBuffer_Vector, MPI::DOUBLE, receive_from, 3);
				
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();
					for (iVar = 0; iVar < nVar; iVar++) {
						node[iPoint]->SetSolution(iVar, Buffer_Receive_Turb[iVar*nVertex+iVertex]);
						node[iPoint]->SetGradient(iVar, 0, Buffer_Receive_Turbx[iVar*nVertex+iVertex]);
						node[iPoint]->SetGradient(iVar, 1, Buffer_Receive_Turby[iVar*nVertex+iVertex]);
						if (nDim == 3) node[iPoint]->SetGradient(iVar, 2, Buffer_Receive_Turbz[iVar*nVertex+iVertex]);
					}
				}
				
				delete [] Buffer_Receive_Turb;
				delete [] Buffer_Receive_Turbx;
				delete [] Buffer_Receive_Turby;
				delete [] Buffer_Receive_Turbz;
				
			}
		}
#endif
}

void CTurbSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
								 CConfig *config, unsigned short val_marker) {
	/*--- Convective fluxes across symmetry plane are equal to zero. ---*/
}

void CTurbSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta, Delta_flow, *local_Residual, Vol;
	double density_old, density;

	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar,0.0);

	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			local_Residual = node[iPoint]->GetResidual();
			Vol = geometry->node[iPoint]->GetVolume();

			/*--- Modify matrix diagonal to assure diagonal dominance ---*/
			Delta_flow = Vol / (solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() + EPS);
			Delta = Delta_flow;
			Jacobian.AddVal2Diag(iPoint,Delta);

			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;

				/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
				rhs[total_index] = -local_Residual[iVar];
				xsol[total_index] = 0.0;
				AddRes_Max( iVar, local_Residual[iVar]*local_Residual[iVar]*Vol );
			}
		}

	/*--- Solve the linear system (Stationary iterative methods) ---*/
	if (config->GetKind_Linear_Solver() == SYM_GAUSS_SEIDEL) 
		Jacobian.SGSSolution(rhs, xsol, config->GetLinear_Solver_Error(), 
												 config->GetLinear_Solver_Iter(), false, geometry, config);
	
	if (config->GetKind_Linear_Solver() == LU_SGS) 
    Jacobian.LU_SGSIteration(rhs, xsol, geometry, config);
	
	/*--- Solve the linear system (Krylov subspace methods) ---*/
	if ((config->GetKind_Linear_Solver() == BCGSTAB) || 
			(config->GetKind_Linear_Solver() == GMRES)) {
		
		CSysVector rhs_vec(geometry->GetnPoint(), geometry->GetnPointDomain(), nVar, rhs);
		CSysVector sol_vec(geometry->GetnPoint(), geometry->GetnPointDomain(), nVar, xsol);
		
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
			system.FlexibleGMRES(rhs_vec, sol_vec, *mat_vec, *precond, *sol_mpi, config->GetLinear_Solver_Error(), 
													 config->GetLinear_Solver_Iter(), false);		
		
		sol_vec.CopyToArray(xsol);
		delete mat_vec; 
		delete precond;  
    delete sol_mpi;
	}
	
	/*--- Update solution (system written in terms of increments) ---*/
	switch (config->GetKind_Turb_Model()){
	case SA:
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
			for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);
		break;
	case SST:
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++){
			density_old = solution_container[FLOW_SOL]->node[iPoint]->GetSolution_Old(0);
			density     = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(0);

			for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->AddConservativeSolution(iVar,xsol[iPoint*nVar+iVar],
						density,density_old,lowerlimit[iVar],upperlimit[iVar]);
		}
		break;
	}

#ifdef NO_MPI
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, sqrt(GetRes_Max(iVar)));
#endif
	
}

CTurbSASolution::CTurbSASolution(void) : CTurbSolution() {}

CTurbSASolution::CTurbSASolution(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTurbSolution() {
	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double Density_Inf, Viscosity_Inf, Factor_nu_Inf;
	
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	bool incompressible = config->GetIncompressible();

	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
	
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	node = new CVariable*[geometry->GetnPoint()];

	/*--- Dimension of the problem --> dependent of the turbulent model ---*/
	nVar = 1;

	/*--- Single grid simulation ---*/
	if (iMesh == MESH_0) {

		/*--- Define some auxiliar vector related with the residual ---*/
		Residual = new double[nVar]; Residual_Max = new double[nVar];
		Residual_i = new double[nVar]; Residual_j = new double[nVar];
		
		/*--- Define some auxiliar vector related with the solution ---*/
		Solution = new double[nVar];
		Solution_i = new double[nVar]; Solution_j = new double[nVar];
		
		/*--- Define some auxiliar vector related with the geometry ---*/
		Vector_i = new double[nDim]; Vector_j = new double[nDim];
		
		/*--- Define some auxiliar vector related with the flow solution ---*/
		FlowSolution_i = new double [nDim+2]; FlowSolution_j = new double [nDim+2];
				
		/*--- Jacobians and vector structures for implicit computations ---*/
		if (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT) {
			/*--- Point to point Jacobians ---*/
			Jacobian_i = new double* [nVar];
			Jacobian_j = new double* [nVar];
			for (iVar = 0; iVar < nVar; iVar++) {
				Jacobian_i[iVar] = new double [nVar];
				Jacobian_j[iVar] = new double [nVar];
			}
			/*--- Initialization of the structure of the whole Jacobian ---*/
			if (rank == MASTER_NODE) cout << "Initialize jacobian structure (SA model)." << endl;
			Initialize_Jacobian_Structure(geometry, config);
			xsol = new double [geometry->GetnPoint()*nVar];
			rhs = new double [geometry->GetnPoint()*nVar];
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

	/*--- Read farfield conditions from config ---*/
	Density_Inf   = config->GetDensity_FreeStreamND();
	Viscosity_Inf = config->GetViscosity_FreeStreamND();
	
	/*--- Factor_nu_Inf in [3.0, 5.0] ---*/
	Factor_nu_Inf = 3.0;
	nu_tilde_Inf  = Factor_nu_Inf*Viscosity_Inf/Density_Inf;

	/*--- Eddy viscosity ---*/
	double Ji, Ji_3, fv1, cv1_3 = 7.1*7.1*7.1;
	double muT_Inf;
	Ji = nu_tilde_Inf/Viscosity_Inf*Density_Inf;
	Ji_3 = Ji*Ji*Ji;
	fv1 = Ji_3/(Ji_3+cv1_3);
	muT_Inf = Density_Inf*fv1*nu_tilde_Inf;

	/*--- Restart the solution from file information ---*/
	if (!restart || geometry->GetFinestMGLevel() == false) {
		if (config->GetKind_Turb_Model() == SA)
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
				node[iPoint] = new CTurbSAVariable(nu_tilde_Inf,muT_Inf, nDim, nVar, config);
	}
	else {
    
		/*--- Restart the solution from file information ---*/
    ifstream restart_file;
		string filename = config->GetSolution_FlowFileName();
		
    /*--- Append time step for unsteady restart ---*/
    if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
      char buffer[50];
      unsigned long flowIter = config->GetnExtIter() - 1;
      filename.erase (filename.end()-4, filename.end());
      if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
      if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
      if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
      if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
      if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
      string UnstExt = string(buffer);
      filename.append(UnstExt);
    }
    restart_file.open(filename.data(), ios::in);
    
		if (restart_file.fail()) {
			cout << "There is no turbulent restart file!!" << endl;
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
    long iPoint_Local; unsigned long iPoint_Global = 0; string text_line;
    double U[5], Velocity2, Pressure, Temperature, Temperature_Dim, LaminarViscosity, muT = 0.0;
    double Gas_Constant = config->GetGas_Constant()/config->GetGas_Constant_Ref();
    while (getline (restart_file,text_line)) {
			istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives 
       on a different processor, the value of iPoint_Local will be -1. 
       Otherwise, the local index for this node on the current processor 
       will be returned and used to instantiate the vars. ---*/
      iPoint_Local = Global2Local[iPoint_Global];
      if (iPoint_Local >= 0) {
        
				if (incompressible) {
					if (nDim == 2) point_line >> index >> U[0] >> U[1] >> U[2] >> Solution[0];
					if (nDim == 3) point_line >> index >> U[0] >> U[1] >> U[2] >> U[3] >> Solution[0];	
				}
				else {
					if (nDim == 2) point_line >> index >> U[0] >> U[1] >> U[2] >> U[3] >> Solution[0];
					if (nDim == 3) point_line >> index >> U[0] >> U[1] >> U[2] >> U[3] >> U[4] >> Solution[0];
				}
				
        /*--- Compute the eddy viscosity from the conservative variables and nu ---*/
        Velocity2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity2 += (U[iDim+1]/U[0])*(U[iDim+1]/U[0]);
        Pressure    = Gamma_Minus_One*U[0]*(U[nDim+1]/U[0]-0.5*Velocity2);
        Temperature = Pressure / ( Gas_Constant * U[0]);
        /*--- Calculate lam. viscosity from a non-dim. Sutherland's Law ---*/
        Temperature_Dim  = Temperature*config->GetTemperature_Ref();
        LaminarViscosity = 1.853E-5*(pow(Temperature_Dim/300.0,3.0/2.0) 
                                    *(300.0+110.3)/(Temperature_Dim+110.3));
        LaminarViscosity = LaminarViscosity/config->GetViscosity_Ref();
        Ji   = Solution[0]/(LaminarViscosity/U[0]);
        Ji_3 = Ji*Ji*Ji;
        fv1  = Ji_3/(Ji_3+cv1_3);
        muT  = U[0]*fv1*Solution[0];
        
        /*--- Instantiate the solution at this node ---*/
        node[iPoint_Local] = new CTurbSAVariable(Solution[0], muT, nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
      node[iPoint] = new CTurbSAVariable(Solution[0], muT, nDim, nVar, config);
    }
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
	}

}

CTurbSASolution::~CTurbSASolution(void) {
	unsigned short iVar, iDim;

	delete [] Residual; delete [] Residual_Max;
	delete [] Residual_i; delete [] Residual_j;
	delete [] Solution;
	delete [] Solution_i; delete [] Solution_j;
	delete [] Vector_i; delete [] Vector_j;

	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Jacobian_i[iVar];
		delete [] Jacobian_j[iVar];
	}
	delete [] Jacobian_i; delete [] Jacobian_j;

	delete [] xsol; delete [] rhs;

	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
	
}


void CTurbSASolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		node[iPoint]->SetResidualZero();
	Jacobian.SetValZero();

	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
}

void CTurbSASolution::Postprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {
	double rho, mu, nu, *nu_hat, muT, Ji, Ji_3, fv1;
	double cv1_3 = 7.1*7.1*7.1;
	unsigned long iPoint;

	bool incompressible = config->GetIncompressible();

	/*--- Compute eddy viscosity ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		
		if (incompressible) {
			rho = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
			mu  = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc();
		}
		else {
			rho = solution_container[FLOW_SOL]->node[iPoint]->GetDensity();
			mu  = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
		}

		nu  = mu/rho;
		nu_hat = node[iPoint]->GetSolution();
		
		Ji   = nu_hat[0]/nu;
		Ji_3 = Ji*Ji*Ji;
		fv1  = Ji_3/(Ji_3+cv1_3);
		
		muT = rho*fv1*nu_hat[0];
		node[iPoint]->SetmuT(muT);
	}
}

void CTurbSASolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {

	double *turb_var_i, *turb_var_j, *U_i, *U_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	
	bool high_order_diss = (config->GetKind_Upwind_Turb() == SCALAR_UPWIND_2ND);
	bool rotating_frame = config->GetRotating_Frame();
  bool grid_movement = config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();

  if (high_order_diss) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solution_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solution_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
		if (config->GetKind_SlopeLimit() != NONE) SetSolution_Limiter(geometry, config);
	}

	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());

		/*--- Conservative variables w/o reconstruction ---*/
		U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
		U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
		solver->SetConservative(U_i, U_j);

		/*--- Turbulent variables w/o reconstruction ---*/
		turb_var_i = node[iPoint]->GetSolution();
		turb_var_j = node[jPoint]->GetSolution();
		solver->SetTurbVar(turb_var_i,turb_var_j);

		/*--- Incompresible density w/o reconstruction ---*/
		if (incompressible)
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 
														solution_container[FLOW_SOL]->node[jPoint]->GetDensityInc());
		
		/*--- Rotating Frame ---*/
		if (rotating_frame)
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());
    
    /*--- Grid Movement ---*/
		if (grid_movement)
			solver->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
		if (high_order_diss) {
			/*--- Conservative solution using gradient reconstruction ---*/
			for (iDim = 0; iDim < nDim; iDim++) {
				Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
				Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}
			Gradient_i = solution_container[FLOW_SOL]->node[iPoint]->GetGradient();
			Gradient_j = solution_container[FLOW_SOL]->node[jPoint]->GetGradient();
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				FlowSolution_i[iVar] = U_i[iVar] + Project_Grad_i;
				FlowSolution_j[iVar] = U_j[iVar] + Project_Grad_j;
			}
			solver->SetConservative(FlowSolution_i, FlowSolution_j);

			/*--- Turbulent variables using gradient reconstruction ---*/
			Gradient_i = node[iPoint]->GetGradient();
			Gradient_j = node[jPoint]->GetGradient();

			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				Solution_i[iVar] = turb_var_i[iVar] + Project_Grad_i;
				Solution_j[iVar] = turb_var_j[iVar] + Project_Grad_j;
			}
			solver->SetTurbVar(Solution_i, Solution_j);
		}

		/*--- Add and subtract Residual ---*/
		solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
		node[iPoint]->AddResidual(Residual);
		node[jPoint]->SubtractResidual(Residual);

		/*--- Implicit part ---*/
		Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
		Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
		Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
		Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);
	}

}

void CTurbSASolution::Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
										   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

	unsigned long iEdge, iPoint, jPoint;
	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();

	if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit) {
		
		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

			/*--- Points in edge ---*/
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);

			/*--- Points coordinates, and normal vector ---*/
			solver->SetCoord(geometry->node[iPoint]->GetCoord(),
														geometry->node[jPoint]->GetCoord());
			solver->SetNormal(geometry->edge[iEdge]->GetNormal());

			/*--- Conservative variables w/o reconstruction ---*/
			solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
															solution_container[FLOW_SOL]->node[jPoint]->GetSolution());

			/*--- Laminar Viscosity, and density (incompresible solver)  ---*/
			if (incompressible) {
				solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
																		solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosityInc());
				solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 
															solution_container[FLOW_SOL]->node[jPoint]->GetDensityInc());
			}
			else {
				solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
																		solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
			}
			
			/*--- Eddy Viscosity ---*/
			solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
                               solution_container[FLOW_SOL]->node[jPoint]->GetEddyViscosity());

			/*--- Turbulent variables w/o reconstruction, and its gradients ---*/
			solver->SetTurbVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
			solver->SetTurbVarGradient(node[iPoint]->GetGradient(),node[jPoint]->GetGradient());

			/*--- Compute residual, and Jacobians ---*/
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);

			/*--- Add and subtract residual, and update Jacobians ---*/
			node[iPoint]->SubtractResidual(Residual);
			node[jPoint]->AddResidual(Residual);

			Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
			Jacobian.SubtractBlock(iPoint,jPoint,Jacobian_j);
			Jacobian.AddBlock(jPoint,iPoint,Jacobian_i);
			Jacobian.AddBlock(jPoint,jPoint,Jacobian_j);
		}
	}

}

void CTurbSASolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
											   CConfig *config, unsigned short iMesh) {
	unsigned long iPoint;

	bool incompressible = config->GetIncompressible();
  bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
  bool transition = (config->GetKind_Trans_Model() == LM);
//  transition = false;
  
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

		/*--- Conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);

		/*--- Gradient of the primitive and conservative variables ---*/
		solver->SetPrimVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);

		/*--- Laminar viscosity and density (incompressible solver) ---*/
		if (incompressible) {
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(), 0.0);
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 0.0);
		}
		else {
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
		}

    if (transition) {
      solver->SetIntermittency(solution_container[TRANS_SOL]->node[iPoint]->GetIntermittency() );
    }

		/*--- Turbulent variables w/o reconstruction, and its gradient ---*/
		solver->SetTurbVar(node[iPoint]->GetSolution(), NULL);
		solver->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);

		/*--- Set volume ---*/
		solver->SetVolume(geometry->node[iPoint]->GetVolume());

		/*--- Set distance to the surface ---*/
		solver->SetDistance(geometry->node[iPoint]->GetWallDistance(), 0.0);

		/*--- Compute the source term ---*/
		solver->SetResidual(Residual, Jacobian_i, NULL, config);

		/*--- Subtract residual and the jacobian ---*/
		node[iPoint]->SubtractResidual(Residual);

		Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
	}
  
  if (time_spectral) {
    
		double Volume, Source;
    unsigned short nVar_Turb = solution_container[TURB_SOL]->GetnVar();
    
		/*--- Loop over points ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
      
			/*--- Get control volume ---*/
			Volume = geometry->node[iPoint]->GetVolume();
      
			/*--- Access stored time spectral source term ---*/
			for (unsigned short iVar = 0; iVar < nVar_Turb; iVar++) {
				Source = node[iPoint]->GetTimeSpectral_Source(iVar);
				Residual[iVar] = Source*Volume;
			}
      
			/*--- Add Residual ---*/
			node[iPoint]->AddResidual(Residual);
      
		}
	}
}


void CTurbSASolution::Source_Template(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
											   CConfig *config, unsigned short iMesh) {

}

void CTurbSASolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {

	unsigned long iPoint, iVertex;
	unsigned short iVar;
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Get the velocity vector ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Solution[iVar] = 0.0;
			
			node[iPoint]->SetSolution_Old(Solution);
			node[iPoint]->SetResidualZero();
			
			/*--- includes 1 in the diagonal ---*/
			Jacobian.DeleteValsRowi(iPoint);
		}
	}

}

void CTurbSASolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {

	unsigned long iPoint, iVertex;
	unsigned short iVar, iDim;
	double *Normal;

	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement	= config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();
  bool gravity        = config->GetGravityForce();

	Normal = new double[nDim];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Retrieve solution at the farfield boundary node ---*/
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			
			/*--- Construct solution state at infinity (farfield) ---*/
			if (incompressible) {
				double yFreeSurface = config->GetFreeSurface_Zero();
				double PressFreeSurface = solution_container[FLOW_SOL]->GetPressure_Inf();
				double Density = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
				double Froude = config->GetFroude();
				double yCoord = geometry->node[iPoint]->GetCoord(1);
				double Pressure;
				if (gravity) Pressure = PressFreeSurface + Density*(yFreeSurface-yCoord)/(Froude*Froude);
				else Pressure = PressFreeSurface;
				
				FlowSolution_j[0] = Pressure;
				FlowSolution_j[1] = solution_container[FLOW_SOL]->GetVelocity_Inf(0)*Density;
				FlowSolution_j[2] = solution_container[FLOW_SOL]->GetVelocity_Inf(1)*Density;
				if (nDim == 3) FlowSolution_j[3] = solution_container[FLOW_SOL]->GetVelocity_Inf(2)*Density;			
			}
			else {
				FlowSolution_j[0] = solution_container[FLOW_SOL]->GetDensity_Inf(); 
				FlowSolution_j[nDim+1] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
				for (iDim = 0; iDim < nDim; iDim++)
					FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(iDim);
			}
			
			/*--- Rotating Frame ---*/
			if (rotating_frame)
				solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
			
			/*--- Grid Movement ---*/
			if (grid_movement)
				solver->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
			
			/*--- Incompressible --*/
			if (incompressible)
				solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 
															solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			
			solver->SetConservative(FlowSolution_i, FlowSolution_j); 
			
			/*--- Set turbulent variable at the wall, and at infinity ---*/
			for (iVar = 0; iVar < nVar; iVar++) 
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
			Solution_j[0] = nu_tilde_Inf;
			solver->SetTurbVar(Solution_i, Solution_j);
			
			/*--- Set Normal (it is necessary to change the sign) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = -Normal[iDim]; 
			solver->SetNormal(Normal);
			
			/*--- Compute residuals and jacobians ---*/
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			
			/*--- Add residuals and jacobians ---*/
			node[iPoint]->AddResidual(Residual);
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			
		}
	}
	
	delete [] Normal; 

}

void CTurbSASolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
  
  /*--- Local variables and initialization. ---*/
	unsigned short iVar, iDim, Kind_Inlet = config->GetKind_Inlet();
  
	unsigned long iVertex, iPoint, Point_Normal;
  
	double P_Total, T_Total, Velocity[3];
	double Velocity2, H_Total, Temperature, Riemann;
	double Pressure, Density, Energy, *Flow_Dir, Mach2;
	double SoundSpeed2, SoundSpeed_Total2, Vel_Mag;
  double alpha, aa, bb, cc, dd;
  double Two_Gamma_M1 = 2.0/Gamma_Minus_One;
  double Gas_Constant = config->GetGas_Constant()/config->GetGas_Constant_Ref();
  double *Normal = new double[nDim];
  
	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	bool rotating_frame = config->GetRotating_Frame();
  bool grid_movement  = config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();
  
  string Marker_Tag = config->GetMarker_All_Tag(val_marker);

  /*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      double Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
      
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;
      
			/*--- Current conservative variables at this boundary node (U_domain) ---*/
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			
			/*--- Build the fictitious intlet state based on characteristics ---*/
			if (incompressible) {
				
				/*--- Index of the closest interior node ---*/
				Point_Normal = iPoint; //geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor();
				
				/*--- Pressure computation using the internal value ---*/
				FlowSolution_j[0] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(0);
				
				/*--- The velocity is computed from the interior, and normal 
         derivative for the density ---*/
				for (iDim = 0; iDim < nDim; iDim++)
					FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetVelocity_Inf(iDim)*solution_container[FLOW_SOL]->node[Point_Normal]->GetDensityInc();
				
				/*--- Note that the y velocity is recomputed due to the 
         free surface effect on the pressure ---*/
				FlowSolution_j[nDim] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(nDim);
				
			}
			else {
        
        /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
         therefore we can specify all but one state variable at the inlet.
         The outgoing Riemann invariant provides the final piece of info.
         Adapted from an original implementation in the Stanford University
         multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
         written by Edwin van der Weide, last modified 04-20-2009. ---*/
        
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
            Density = FlowSolution_i[0];
            Velocity2 = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Velocity[iDim] = FlowSolution_i[iDim+1]/Density;
              Velocity2 += Velocity[iDim]*Velocity[iDim];
            }
            Energy      = FlowSolution_i[nVar-1]/Density;
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
            FlowSolution_j[0] = Density;
            FlowSolution_j[1] = Velocity[0]*Density;
            FlowSolution_j[2] = Velocity[1]*Density;
            FlowSolution_j[3] = Energy*Density;
            if (nDim == 3) {
              FlowSolution_j[3] = Velocity[2]*Density;
              FlowSolution_j[4] = Energy*Density;
            }
            
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
              Velocity[iDim] = node[iPoint]->GetVelocity(iDim, incompressible);
            Pressure    = node[iPoint]->GetPressure(incompressible);
            SoundSpeed2 = Gamma*Pressure/FlowSolution_i[0];
            
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
            FlowSolution_j[0] = Density;
            FlowSolution_j[1] = Vel_Mag*Flow_Dir[0]*Density;
            FlowSolution_j[2] = Vel_Mag*Flow_Dir[1]*Density;
            FlowSolution_j[3] = Energy*Density;
            if (nDim == 3) {
              FlowSolution_j[3] = Vel_Mag*Flow_Dir[2]*Density;
              FlowSolution_j[4] = Energy*Density;
            }
            
            break;
        }
			}
      
      /*--- Set the conservative variable states. Note that we really only
       need the density and momentum components so that we can compute
       the velocity for the convective part of the upwind residual. ---*/
			solver->SetConservative(FlowSolution_i,FlowSolution_j);
      
      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      Solution_j[iVar]= nu_tilde_Inf;
			solver->SetTurbVar(Solution_i,Solution_j);

			/*--- Set various other quantities in the solver class ---*/
			solver->SetNormal(Normal);
      if (incompressible)
				solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
															solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      if (rotating_frame) {
        solver->SetRotVel(geometry->node[iPoint]->GetRotVel(),
                          geometry->node[iPoint]->GetRotVel());
        solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
      }
      if (grid_movement)
        solver->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                           geometry->node[iPoint]->GetGridVel());
      
			/*--- Compute the residual using an upwind scheme ---*/
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddResidual(Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
      
		}
	}
  
	/*--- Free locally allocated memory ---*/
  delete[] Normal;

}

void CTurbSASolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																CConfig *config, unsigned short val_marker) {
	
  /*--- Local variables and initialization. ---*/
  unsigned long iPoint, iVertex;
	unsigned short iVar, iDim;
  
  bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();
  
	double *Normal = new double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Set the conservative variables, density & Velocity same as in 
			 the interior, don't need to modify pressure for the convection. ---*/
			solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), 
															solution_container[FLOW_SOL]->node[iPoint]->GetSolution()); 
			
      /*--- Set the turbulent variables. Here we use a Neumann BC such
			 that the turbulent variable is copied from the interior of the
			 domain to the outlet before computing the residual.
			 Solution_i --> TurbVar_internal,
			 Solution_j --> TurbVar_outlet ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
				Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
			}
			solver->SetTurbVar(Solution_i, Solution_j);
      
			/*--- Set Normal (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = -Normal[iDim];
			solver->SetNormal(Normal);
      
			/*--- Set various quantities in the solver class ---*/
			if (incompressible)
				solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 
															solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      if (rotating_frame) {
				solver->SetRotVel(geometry->node[iPoint]->GetRotVel(),
                          geometry->node[iPoint]->GetRotVel());
				solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
      if (grid_movement)
        solver->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                           geometry->node[iPoint]->GetGridVel());
			
			/*--- Compute the residual using an upwind scheme ---*/
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddResidual(Residual);
			
      /*--- Jacobian contribution for implicit integration ---*/
			if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
		}
	}
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}

void CTurbSASolution::SetResidual_DualTime(CGeometry *geometry, CSolution **solution_container, CConfig *config,
                                          unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {

  unsigned short iVar, jVar;
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1;
  double Volume_nM1, Volume_n, Volume_nP1, TimeStep;
	
	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	bool grid_movement = config->GetGrid_Movement();
  
	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n   = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();

		/*--- Volume at time n-1 and n ---*/
		if (grid_movement) {
			Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
			Volume_n   = geometry->node[iPoint]->GetVolume_n();
			Volume_nP1 = geometry->node[iPoint]->GetVolume();
		}
		else {
			Volume_nM1 = geometry->node[iPoint]->GetVolume();
			Volume_n   = geometry->node[iPoint]->GetVolume();
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

		/*--- Add Residual ---*/
		node[iPoint]->AddResidual(Residual);

		if (implicit) {
			for (iVar = 0; iVar < nVar; iVar++) {
				for (jVar = 0; jVar < nVar; jVar++)
					Jacobian_i[iVar][jVar] = 0.0;
				
				if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
					Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
				if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
					Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
			}
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}
}

CTurbSSTSolution::CTurbSSTSolution(void) : CTurbSolution() {}

CTurbSSTSolution::CTurbSSTSolution(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTurbSolution() {

	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double dull_val;
	ifstream restart_file;
	string text_line;
	
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	bool incompressible = config->GetIncompressible();

	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	node = new CVariable*[geometry->GetnPoint()];

	/*--- Dimension of the problem --> dependent of the turbulent model ---*/
	nVar = 2;

	/*--- Single grid simulation ---*/
	if (iMesh == MESH_0) {
		
		/*--- Define some auxiliary vector related with the residual ---*/
		Residual = new double[nVar]; Residual_Max = new double[nVar];
		Residual_i = new double[nVar]; Residual_j = new double[nVar];
		
		/*--- Define some auxiliary vector related with the solution ---*/
		Solution = new double[nVar];
		Solution_i = new double[nVar]; Solution_j = new double[nVar];
		
		/*--- Define some auxiliary vector related with the geometry ---*/
		Vector_i = new double[nDim]; Vector_j = new double[nDim];
		
		/*--- Define some auxiliary vector related with the flow solution ---*/
		FlowSolution_i = new double [nDim+2]; FlowSolution_j = new double [nDim+2];
		
		/*--- Jacobians and vector structures for implicit computations ---*/
		if (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT) {
			/*--- Point to point Jacobians ---*/
			Jacobian_i = new double* [nVar];
			Jacobian_j = new double* [nVar];
			for (iVar = 0; iVar < nVar; iVar++) {
				Jacobian_i[iVar] = new double [nVar];
				Jacobian_j[iVar] = new double [nVar];
			}
			/*--- Initialization of the structure of the whole Jacobian ---*/
			if (rank == MASTER_NODE) cout << "Initialize jacobian structure (SST model)." << endl;
			Initialize_Jacobian_Structure(geometry, config);
			xsol = new double [geometry->GetnPoint()*nVar];
			rhs = new double [geometry->GetnPoint()*nVar];
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

	/* --- Initialize value for model constants --- */
	constants = new double[10];
	constants[0] = 0.85;   //sigma_k1
	constants[1] = 1.0;    //sigma_k2
	constants[2] = 0.5;    //sigma_om1
	constants[3] = 0.856;  //sigma_om2
	constants[4] = 0.075;  //beta_1
	constants[5] = 0.0828; //beta_2
	constants[6] = 0.09;   //betaStar
	constants[7] = 0.31;   //a1
	constants[8] = constants[4]/constants[6] - constants[2]*0.41*0.41/sqrt(constants[6]);  //alfa_1
	constants[9] = constants[5]/constants[6] - constants[3]*0.41*0.41/sqrt(constants[6]);  //alfa_2

	/* --- Initialize lower and upper limits--- */
	lowerlimit = new double[nVar];
	upperlimit = new double[nVar];

	lowerlimit[0] = 1.0e-10;
	upperlimit[0] = 1.0e10;

	lowerlimit[1] = 1.0e-4;
	upperlimit[1] = 1.0e15;

	/*--- Flow infinity initialization stuff ---*/
	double rhoInf, *VelInf, muLamInf, Intensity, viscRatio, muT_Inf;

	rhoInf    = config->GetDensity_FreeStreamND();
	VelInf    = config->GetVelocity_FreeStreamND();
	muLamInf  = config->GetViscosity_FreeStreamND();
	Intensity = config->GetTurbulenceIntensity_FreeStream();
	viscRatio = config->GetTurb2LamViscRatio_FreeStream();

	double VelMag = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		VelMag += VelInf[iDim]*VelInf[iDim];
	VelMag = sqrt(VelMag);

	kine_Inf  = 3.0/2.0*(VelMag*VelMag*Intensity*Intensity);
	omega_Inf = rhoInf*kine_Inf/(muLamInf*viscRatio);

	/*--- Eddy viscosity, initialized without stress limiter ---*/
	muT_Inf = rhoInf*kine_Inf/omega_Inf;

	/*--- Restart the solution from file information ---*/
	if (!restart || geometry->GetFinestMGLevel() == false) {
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CTurbSSTVariable(kine_Inf,omega_Inf, muT_Inf, nDim, nVar, constants, config);
	}
	else {
		/*--- Restart the solution from file information ---*/
		string mesh_filename = config->GetSolution_FlowFileName();
		restart_file.open(mesh_filename.data(), ios::in);
    
		if (restart_file.fail()) {
			cout << "There is no turbulent restart file!!" << endl;
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
		while (getline (restart_file,text_line)) {
			istringstream point_line(text_line);
      
			/*--- Retrieve local index. If this node from the restart file lives
       	   on a different processor, the value of iPoint_Local will be -1.
       	   Otherwise, the local index for this node on the current processor
       	   will be returned and used to instantiate the vars. ---*/
			iPoint_Local = Global2Local[iPoint_Global];
			if (iPoint_Local >= 0) {
        
				if (incompressible) {
					if (nDim == 2) point_line >> index >> rhoInf >> dull_val >> dull_val >> Solution[0] >> Solution[1];
					if (nDim == 3) point_line >> index >> rhoInf >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
				}
				else {
					if (nDim == 2) point_line >> index >> rhoInf >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
					if (nDim == 3) point_line >> index >> rhoInf >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
				}
				
        /*--- Eddy viscosity, initialized without stress limiter ---*/
        muT_Inf = rhoInf*Solution[0]/Solution[1];
        
				/*--- Instantiate the solution at this node ---*/
				node[iPoint_Local] = new CTurbSSTVariable(Solution[0], Solution[1], muT_Inf, nDim, nVar, constants, config);
			}
			iPoint_Global++;
		}
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
      node[iPoint] = new CTurbSSTVariable(Solution[0], Solution[1], muT_Inf, nDim, nVar, constants, config);
    }
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
	}

}

CTurbSSTSolution::~CTurbSSTSolution(void) {
	unsigned short iVar, iDim;

	delete [] Residual; delete [] Residual_Max;
	delete [] Residual_i; delete [] Residual_j;
	delete [] Solution;
	delete [] Solution_i; delete [] Solution_j;
	delete [] Vector_i; delete [] Vector_j;
	delete [] FlowSolution_i; delete [] FlowSolution_j;

	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Jacobian_i[iVar];
		delete [] Jacobian_j[iVar];
	}
	delete [] Jacobian_i; delete [] Jacobian_j;

	delete [] xsol; delete [] rhs;

	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;

  delete [] constants;
  delete [] upperlimit;
  delete [] lowerlimit;
  
}

void CTurbSSTSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++){
		node[iPoint]->SetResidualZero();
	}
	Jacobian.SetValZero();

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  
  // Compute mean flow gradients
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
		solution_container[FLOW_SOL]->SetPrimVar_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
		solution_container[FLOW_SOL]->SetPrimVar_Gradient_LS(geometry, config);
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
    solution_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
    solution_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
  
}

void CTurbSSTSolution::Postprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {
	double rho, mu, dist, omega, kine, vorticity[3], vortMag, strMag, F2, muT, zeta;
	double a1 = constants[7];
	unsigned long iPoint;
	
	bool incompressible = config->GetIncompressible();
	
	// Compute mean flow gradients
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
		solution_container[FLOW_SOL]->SetPrimVar_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
		solution_container[FLOW_SOL]->SetPrimVar_Gradient_LS(geometry, config);

	// Compute turbulence variable gradients
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
		SetSolution_Gradient_GG(geometry);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
		SetSolution_Gradient_LS(geometry, config);

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		// Compute vorticity and rate of strain magnitude
		solution_container[FLOW_SOL]->node[iPoint]->SetVorticity();
    vorticity[0] = solution_container[FLOW_SOL]->node[iPoint]->GetVorticity(0);
    vorticity[1] = solution_container[FLOW_SOL]->node[iPoint]->GetVorticity(1);
    vorticity[2] = solution_container[FLOW_SOL]->node[iPoint]->GetVorticity(2);
    vortMag = sqrt(vorticity[0]*vorticity[0] + vorticity[1]*vorticity[1] + vorticity[2]*vorticity[2]);
		//vort   = 0.0; //not yet implemented
    solution_container[FLOW_SOL]->node[iPoint]->SetStrainMag();
		strMag = solution_container[FLOW_SOL]->node[iPoint]->GetStrainMag();

		// Compute blending functions and cross diffusion
		if (incompressible) {
			rho  = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
			mu   = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc();
		}
		else {
			rho  = solution_container[FLOW_SOL]->node[iPoint]->GetDensity();
			mu   = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
		}
		dist = geometry->node[iPoint]->GetWallDistance();

		node[iPoint]->SetBlendingFunc(mu,dist,rho);
		F2 = node[iPoint]->GetF2blending();

		// Compute the eddy viscosity
		kine  = node[iPoint]->GetSolution(0);
		omega = node[iPoint]->GetSolution(1);

    //zeta = min(1.0/omega,a1/(vortMag*F2));
		zeta = min(1.0/omega,a1/(strMag*F2));
		//zeta = 1.0/omega;
		muT = min(max(rho*kine*zeta,0.0),1.0);
		node[iPoint]->SetmuT(muT);
	}
}

void CTurbSSTSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {
	double *turb_var_i, *turb_var_j, *U_i, *U_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	
	bool high_order_diss = (config->GetKind_Upwind_Turb() == SCALAR_UPWIND_2ND);
	bool rotating_frame = config->GetRotating_Frame();
  bool grid_movement = config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();

	if (high_order_diss) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solution_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solution_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
		if (config->GetKind_SlopeLimit() != NONE) SetSolution_Limiter(geometry, config);
	}

	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		
		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());

		/*--- Conservative variables w/o reconstruction ---*/
		U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
		U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
		solver->SetConservative(U_i, U_j);

		/*--- Turbulent variables w/o reconstruction ---*/
		turb_var_i = node[iPoint]->GetSolution();
		turb_var_j = node[jPoint]->GetSolution();
		solver->SetTurbVar(turb_var_i,turb_var_j);

		/*--- Incompresible density w/o reconstruction ---*/
		if (incompressible)
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 
														solution_container[FLOW_SOL]->node[jPoint]->GetDensityInc());
		
		/*--- Rotating Frame ---*/
		if (rotating_frame)
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());
    
    /*--- Grid Movement ---*/
		if (grid_movement)
			solver->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
		if (high_order_diss) {
			/*--- Conservative solution using gradient reconstruction ---*/
			for (iDim = 0; iDim < nDim; iDim++) {
				Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
				Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}
			Gradient_i = solution_container[FLOW_SOL]->node[iPoint]->GetGradient();
			Gradient_j = solution_container[FLOW_SOL]->node[jPoint]->GetGradient();
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				FlowSolution_i[iVar] = U_i[iVar] + Project_Grad_i;
				FlowSolution_j[iVar] = U_j[iVar] + Project_Grad_j;
			}
			solver->SetConservative(FlowSolution_i, FlowSolution_j);

			/*--- Turbulent variables using gradient reconstruction ---*/
			Gradient_i = node[iPoint]->GetGradient();
			Gradient_j = node[jPoint]->GetGradient();
			
			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				Solution_i[iVar] = turb_var_i[iVar] + Project_Grad_i;
				Solution_j[iVar] = turb_var_j[iVar] + Project_Grad_j;

			}
			solver->SetTurbVar(Solution_i, Solution_j);
		}

		/*--- Add and subtract Residual ---*/
		solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
		node[iPoint]->AddResidual(Residual);
		node[jPoint]->SubtractResidual(Residual);

		/*--- Implicit part ---*/
		Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
		Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
		Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
		Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);
	}
}

void CTurbSSTSolution::Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
										   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint;
	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();

	if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit) {

		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

			/*--- Points in edge ---*/
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);

			/*--- Points coordinates, and normal vector ---*/
			solver->SetCoord(geometry->node[iPoint]->GetCoord(),
														geometry->node[jPoint]->GetCoord());
			solver->SetNormal(geometry->edge[iEdge]->GetNormal());

			/*--- Conservative variables w/o reconstruction ---*/
			solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
															solution_container[FLOW_SOL]->node[jPoint]->GetSolution());

			/*--- Laminar Viscosity, and density (incompresible solver)  ---*/
			if (incompressible) {
				solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
																		solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosityInc());
				solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 
															solution_container[FLOW_SOL]->node[jPoint]->GetDensityInc());
			}
			else {
				solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
																		solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
			}

			/*--- Eddy Viscosity ---*/
			solver->SetEddyViscosity(node[iPoint]->GetmuT(), node[jPoint]->GetmuT());

			/*--- Turbulent variables w/o reconstruction, and its gradients ---*/
			solver->SetTurbVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
			solver->SetTurbVarGradient(node[iPoint]->GetGradient(),
																 node[jPoint]->GetGradient());

			/*--- Menter's first blending function ---*/
			solver->SetF1blending(node[iPoint]->GetF1blending(),node[jPoint]->GetF1blending());

			/*--- Compute residual, and Jacobians ---*/
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);

			/*--- Add and subtract residual, and update Jacobians ---*/
			node[iPoint]->SubtractResidual(Residual);
			node[jPoint]->AddResidual(Residual);
			Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
			Jacobian.SubtractBlock(iPoint,jPoint,Jacobian_j);
			Jacobian.AddBlock(jPoint,iPoint,Jacobian_i);
			Jacobian.AddBlock(jPoint,jPoint,Jacobian_j);
		}
	}
}

void CTurbSSTSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
											   CConfig *config, unsigned short iMesh) {
	unsigned long iPoint;
	bool incompressible = config->GetIncompressible();

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

		/*--- Conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);

		/*--- Gradient of the primitive and conservative variables ---*/
		solver->SetPrimVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);

		/*--- Laminar viscosity, and density (incompresible solver) ---*/
		if (incompressible) {
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(), 0.0);
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 0.0);
		}
		else {
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
		}
		
		/*--- Eddy Viscosity ---*/
		solver->SetEddyViscosity(node[iPoint]->GetmuT(), 0.0);

		/*--- Turbulent variables w/o reconstruction, and its gradient ---*/
		solver->SetTurbVar(node[iPoint]->GetSolution(), NULL);
		solver->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);

		/*--- Set volume ---*/
		solver->SetVolume(geometry->node[iPoint]->GetVolume());

		/*--- Set distance to the surface ---*/
		solver->SetDistance(geometry->node[iPoint]->GetWallDistance(), 0.0);

		/*--- Menter's first blending function ---*/
		solver->SetF1blending(node[iPoint]->GetF1blending(),0.0);

		/*--- Menter's second blending function ---*/
		solver->SetF2blending(node[iPoint]->GetF2blending(),0.0);

		/*--- Rate of strain magnitude ---*/
		solver->SetStrainMag(solution_container[FLOW_SOL]->node[iPoint]->GetStrainMag(),0.0);

		/*--- Cross diffusion ---*/
		solver->SetCrossDiff(node[iPoint]->GetCrossDiff(),0.0);

		/*--- Compute the source term ---*/
		solver->SetResidual(Residual, Jacobian_i, NULL, config);

		/*--- Subtract residual and the jacobian ---*/
		node[iPoint]->SubtractResidual(Residual);
		Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
	}
}

void CTurbSSTSolution::Source_Template(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
											   CConfig *config, unsigned short iMesh) {

}

void CTurbSSTSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, jPoint, iVertex, total_index;
	unsigned short iDim, iVar;
	double distance, density, laminar_viscosity, beta_1;

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- distance to closest neighbor ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor();
			distance = 0.0;
			for(iDim = 0; iDim < nDim; iDim++){
				distance += (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim))*
						        (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}
			distance = sqrt(distance);

			/*--- Set wall values ---*/
			density = solution_container[FLOW_SOL]->node[jPoint]->GetSolution(0);
			beta_1 = constants[4];
			laminar_viscosity = solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();

			Solution[0] = 0.0;
			Solution[1] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);
			//Solution[1] = 60.0*laminar_viscosity/(beta_1*distance*distance);

      /*--- Set the solution values and zero the residual ---*/
			node[iPoint]->SetSolution_Old(Solution);
			node[iPoint]->SetSolution(Solution);
      node[iPoint]->SetResidualZero();
      
      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
      
		}
	}
}

void CTurbSSTSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	double *Normal;
	unsigned short iVar, iDim;

	Normal = new double[nDim];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Set conservative variables at the wall, and at infinity ---*/
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			
			FlowSolution_j[0] = solution_container[FLOW_SOL]->GetDensity_Inf();
			FlowSolution_j[nDim+1] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
			for (iDim = 0; iDim < nDim; iDim++)
				FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(iDim);
			
			solver->SetConservative(FlowSolution_i, FlowSolution_j);
			
			/*--- Set turbulent variable at the wall, and at infinity ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
			Solution_j[0] = kine_Inf;
			Solution_j[1] = omega_Inf;
			solver->SetTurbVar(Solution_i, Solution_j);
			
			/*--- Set Normal (it is necessary to change the sign) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = -Normal[iDim];
			solver->SetNormal(Normal);
			
			/*--- Compute residuals and jacobians ---*/
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			
			/*--- Add residuals and jacobians ---*/
			node[iPoint]->AddResidual(Residual);
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}
	
	delete [] Normal;
}

void CTurbSSTSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
									  unsigned short val_marker) {

	/*--- Local variables and initialization. ---*/
	unsigned short iVar, iDim, Kind_Inlet = config->GetKind_Inlet();
  
	unsigned long iVertex, iPoint, Point_Normal;
  
	double P_Total, T_Total, Velocity[3];
	double Velocity2, H_Total, Temperature, Riemann;
	double Pressure, Density, Energy, *Flow_Dir, Mach2;
	double SoundSpeed2, SoundSpeed_Total2, Vel_Mag;
    double alpha, aa, bb, cc, dd;
    double Two_Gamma_M1 = 2.0/Gamma_Minus_One;
    double Gas_Constant = config->GetGas_Constant()/config->GetGas_Constant_Ref();
    double *Normal = new double[nDim];
  
	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();
  
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
  
	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
			double Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
      
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;
      
			/*--- Current conservative variables at this boundary node (U_domain) ---*/
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			
			/*--- Build the fictitious intlet state based on characteristics ---*/
			if (incompressible) {
				
				/*--- Index of the closest interior node ---*/
				Point_Normal = iPoint; //geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor();
				
				/*--- Pressure computation using the internal value ---*/
				FlowSolution_j[0] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(0);
				
				/*--- The velocity is computed from the interior, and normal derivative for the density ---*/
				for (iDim = 0; iDim < nDim; iDim++)
					FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetVelocity_Inf(iDim)*
					solution_container[FLOW_SOL]->node[Point_Normal]->GetDensityInc();
				
				/*--- Note that the y velocity is recomputed due to the
                free surface effect on the pressure ---*/
				FlowSolution_j[nDim] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(nDim);
				
			}
			else {
        
				/*--- Subsonic inflow: there is one outgoing characteristic (u-c),
                therefore we can specify all but one state variable at the inlet.
                The outgoing Riemann invariant provides the final piece of info.
                Adapted from an original implementation in the Stanford University
                multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
                written by Edwin van der Weide, last modified 04-20-2009. ---*/
        
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
					Density = FlowSolution_i[0];
					Velocity2 = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						Velocity[iDim] = FlowSolution_i[iDim+1]/Density;
						Velocity2 += Velocity[iDim]*Velocity[iDim];
					}
					Energy      = FlowSolution_i[nVar-1]/Density;
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
					FlowSolution_j[0] = Density;
					FlowSolution_j[1] = Velocity[0]*Density;
					FlowSolution_j[2] = Velocity[1]*Density;
					FlowSolution_j[3] = Energy*Density;
					if (nDim == 3) {
						FlowSolution_j[3] = Velocity[2]*Density;
						FlowSolution_j[4] = Energy*Density;
					}
            
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
						Velocity[iDim] = node[iPoint]->GetVelocity(iDim, incompressible);
					Pressure    = node[iPoint]->GetPressure(incompressible);
					SoundSpeed2 = Gamma*Pressure/FlowSolution_i[0];
            
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
					FlowSolution_j[0] = Density;
					FlowSolution_j[1] = Vel_Mag*Flow_Dir[0]*Density;
					FlowSolution_j[2] = Vel_Mag*Flow_Dir[1]*Density;
					FlowSolution_j[3] = Energy*Density;
					if (nDim == 3) {
						FlowSolution_j[3] = Vel_Mag*Flow_Dir[2]*Density;
						FlowSolution_j[4] = Energy*Density;
					}
            
					break;
				}
			}
      
			/*--- Set the conservative variable states. Note that we really only
            need the density and momentum components so that we can compute
            the velocity for the convective part of the upwind residual. ---*/
			solver->SetConservative(FlowSolution_i,FlowSolution_j);
      
			/*--- Set the turbulent variable states. Use free-stream SST
            values for the turbulent state at the inflow. ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
			Solution_j[0]= kine_Inf;
			Solution_j[1]= omega_Inf;

			solver->SetTurbVar(Solution_i,Solution_j);
      
			/*--- Set various other quantities in the solver class ---*/
			solver->SetNormal(Normal);
			if (incompressible)
				solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
															solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			if (rotating_frame) {
				solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
				solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}

			if (grid_movement)
				solver->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
			/*--- Compute the residual using an upwind scheme ---*/
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddResidual(Residual);
      
			/*--- Jacobian contribution for implicit integration ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
		}
	}
  
	/*--- Free locally allocated memory ---*/
	delete[] Normal;
}

void CTurbSSTSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
									   CConfig *config, unsigned short val_marker) {

	/*--- Local variables and initialization. ---*/
	unsigned long iPoint, iVertex;
	unsigned short iVar, iDim;
  
	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();
  
	double *Normal = new double[nDim];
  
	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Set the conservative variables, density & Velocity same as in
			 the interior, don't need to modify pressure for the convection. ---*/
			solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
															solution_container[FLOW_SOL]->node[iPoint]->GetSolution());
			
			/*--- Set the turbulent variables. Here we use a Neumann BC such
			 that the turbulent variable is copied from the interior of the
			 domain to the outlet before computing the residual.
			 Solution_i --> TurbVar_internal,
			 Solution_j --> TurbVar_outlet ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
				Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
			}
			solver->SetTurbVar(Solution_i, Solution_j);
      
			/*--- Set Normal (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = -Normal[iDim];
			solver->SetNormal(Normal);
      
			/*--- Set various quantities in the solver class ---*/
			if (incompressible)
				solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
															solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			if (rotating_frame) {
				solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
				solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
			
			if (grid_movement)
				solver->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

			/*--- Compute the residual using an upwind scheme ---*/
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddResidual(Residual);
			
			/*--- Jacobian contribution for implicit integration ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}
  
	/*--- Free locally allocated memory ---*/
	delete[] Normal;
}

void CTurbSSTSolution::SetResidual_DualTime(CGeometry *geometry, CSolution **solution_container, CConfig *config,
                                            unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
  
  unsigned short iVar, jVar;
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1;
  double Volume_nM1, Volume_n, Volume_nP1, TimeStep;
	
	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	bool grid_movement = config->GetGrid_Movement();
  
	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n   = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();
    
		/*--- Volume at time n-1 and n ---*/
		if (grid_movement) {
			Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
			Volume_n   = geometry->node[iPoint]->GetVolume_n();
			Volume_nP1 = geometry->node[iPoint]->GetVolume();
		}
		else {
			Volume_nM1 = geometry->node[iPoint]->GetVolume();
			Volume_n   = geometry->node[iPoint]->GetVolume();
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
    
		/*--- Add Residual ---*/
		node[iPoint]->AddResidual(Residual);
    
		if (implicit) {
			for (iVar = 0; iVar < nVar; iVar++) {
				for (jVar = 0; jVar < nVar; jVar++)
					Jacobian_i[iVar][jVar] = 0.0;
				
				if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
					Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
				if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
					Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
			}
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}
}

double* CTurbSSTSolution::GetConstants() {
	return constants;
}
