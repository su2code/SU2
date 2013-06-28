/*!
 * \file solution_direct_turbulent.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.1.
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

void CTurbSolution::MPI_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config,
																		unsigned short val_mesh) {
	
#ifndef NO_MPI
	unsigned short iVar, iMarker;
	double *Turb_Var, **Turb_Grad;
	unsigned long iVertex, iPoint;
	
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			
			short SendRecv = config->GetMarker_All_SendRecv(iMarker);
			unsigned long nVertex = geometry->nVertex[iMarker];
			
			/*--- Send information  ---*/
			if (SendRecv > 0) {
				/*--- Dimensionalization ---*/
				unsigned long nBuffer_Vector = geometry->nVertex[iMarker]*nVar;
				int send_to = SendRecv-1;
				
				double *Buffer_Send_Turb = new double[nBuffer_Vector];
				double *Buffer_Send_Turbx = new double[nBuffer_Vector];
				double *Buffer_Send_Turby = new double[nBuffer_Vector];
				double *Buffer_Send_Turbz = new double[nBuffer_Vector];
				
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					
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
				unsigned long nBuffer_Vector = geometry->nVertex[iMarker]*nVar;
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
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
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
								 CConfig *config, unsigned short val_marker) { }

void CTurbSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta, Delta_flow, *local_Residual, Vol;

	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar,0.0);

	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			local_Residual = node[iPoint]->GetResidual();
			Vol = geometry->node[iPoint]->GetVolume();

			/*--- Modify matrix diagonal to assure diagonal dominance ---*/
			Delta_flow = Vol/(solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
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

	/*--- Solve the system ---*/
	if (config->GetKind_Linear_Solver() == SYM_GAUSS_SEIDEL) Jacobian.SGSSolution(rhs, xsol, config->GetLinear_Solver_Error(), 
																																								config->GetLinear_Solver_Iter(), false, geometry, config);
	if (config->GetKind_Linear_Solver() == BCGSTAB) Jacobian.BCGSTABSolution(rhs, xsol, config->GetLinear_Solver_Error(), 
																																					 config->GetLinear_Solver_Iter(), false, geometry, config);
	if (config->GetKind_Linear_Solver() == LU_SGS) Jacobian.LU_SGSIteration(rhs, xsol, geometry, config);
	
	
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		for (iVar = 0; iVar < nVar; iVar++)
			node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);

#ifdef NO_MPI
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, sqrt(GetRes_Max(iVar)));
#endif
	
}

CTurbSASolution::CTurbSASolution(void) : CTurbSolution() {}

CTurbSASolution::CTurbSASolution(CGeometry *geometry, CConfig *config) : CTurbSolution() {
	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double Density_Inf, Viscosity_Inf, Factor_nu_Inf, dull_val;
	ifstream restart_file;
	char *cstr;
	string text_line;
	bool restart = (config->GetRestart() || config->GetRestart_Flow());

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	node = new CVariable*[geometry->GetnPoint()];

	/*--- Dimension of the problem --> dependent of the turbulent model ---*/
	nVar = 1;

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
		Initialize_Jacobian_Structure(geometry, config);
		xsol = new double [geometry->GetnPoint()*nVar];
		rhs = new double [geometry->GetnPoint()*nVar];
	}

	/*--- Computation of gradients by least squares ---*/
	if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
		(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) {
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			cvector[iVar] = new double [nDim];
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
	if (!restart) {
		if (config->GetKind_Turb_Model() == SA)
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
				node[iPoint] = new CTurbSAVariable(nu_tilde_Inf,muT_Inf, nDim, nVar, config);
		if (config->GetKind_Turb_Model() == SA_COMP)
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
				node[iPoint] = new CTurbSAVariable(Density_Inf*nu_tilde_Inf,muT_Inf, nDim, nVar, config);
	}
	else {
		string mesh_filename = config->GetSolution_FlowFileName();
		cstr = new char [mesh_filename.size()+1];
		strcpy (cstr, mesh_filename.c_str());
		restart_file.open(cstr, ios::in);
		if (restart_file.fail()) {
			cout << "There is no turbulent restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}

		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
			if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
			node[iPoint] = new CTurbSAVariable(Solution[0], 0, nDim, nVar, config);
		}
		restart_file.close();
	}

}

CTurbSASolution::~CTurbSASolution(void) {
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
	
}

void CTurbSASolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		node[iPoint]->SetResidualZero();
	Jacobian.SetValZero();

	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry);
	if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
		(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) SetSolution_Gradient_LS(geometry, config);
}

void CTurbSASolution::Postprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {
	// Compute eddy viscosity
	double rho, mu, nu, *nu_hat, muT;
	double Ji, Ji_3, fv1;
	double cv1_3 = 7.1*7.1*7.1;
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++){
		rho = solution_container[FLOW_SOL]->node[iPoint]->GetDensity();
		mu  = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
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
	double *turb_var_i, *turb_var_j, *Limiter_i = NULL, *Limiter_j = NULL, *U_i, *U_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	bool high_order_diss = (config->GetKind_Upwind_Turb() == SCALAR_UPWIND_2ND);
	bool rotating_frame = config->GetRotating_Frame();

	if (high_order_diss) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solution_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry);
		if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
			(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) solution_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
		if (config->GetKind_SlopeLimit() == VENKATAKRISHNAN) SetSolution_Limiter(geometry, config);
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

		/*--- Rotational frame ---*/
		if (rotating_frame)
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());

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
			if (config->GetKind_SlopeLimit() != NONE) {
				Limiter_i = node[iPoint]->GetLimiter();
				Limiter_j = node[jPoint]->GetLimiter();
			}
			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				if (config->GetKind_SlopeLimit() == NONE) {
					Solution_i[iVar] = turb_var_i[iVar] + Project_Grad_i;
					Solution_j[iVar] = turb_var_j[iVar] + Project_Grad_j;
				}
				else {
					Solution_i[iVar] = turb_var_i[iVar] + Project_Grad_i*Limiter_i[iVar];
					Solution_j[iVar] = turb_var_j[iVar] + Project_Grad_j*Limiter_j[iVar];
				}
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

	if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit) {
		/*--- If SA_COMP --> Need gradient of flow conservative variables ---*/
		if (config->GetKind_Turb_Model() == SA_COMP) {
			if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solution_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry);
			if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
					(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) solution_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
		}

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

			/*--- Laminar Viscosity ---*/
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
																	solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
			/*--- Eddy Viscosity ---*/
			solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
																				solution_container[FLOW_SOL]->node[jPoint]->GetEddyViscosity());

			/*--- Turbulent variables w/o reconstruction, and its gradients ---*/
			solver->SetTurbVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
			solver->SetTurbVarGradient(node[iPoint]->GetGradient(),
																 node[jPoint]->GetGradient());

			if (config->GetKind_Turb_Model() == SA_COMP) {
				solver->SetConsVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient(),
																	 solution_container[FLOW_SOL]->node[jPoint]->GetGradient());
			}

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

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

		/*--- Conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);

		/*--- Gradient of the primitive and conservative variables ---*/
		solver->SetPrimVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
		if (config->GetKind_Turb_Model() == SA_COMP || config->GetKind_Turb_Model() == SST)
			solver->SetConsVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient(), NULL);

		/*--- Laminar viscosity ---*/
		solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);

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
		Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
	}
}

void CTurbSASolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar;

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Get the velocity vector ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Solution[iVar] = 0.0;

		node[iPoint]->SetSolution_Old(Solution);
		node[iPoint]->SetResidualZero();
		
		/*--- includes 1 in the diagonal ---*/
		Jacobian.DeleteValsRowi(iPoint);
	}
}

void CTurbSASolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar, iDim;
	double *Normal;
	bool rotating_frame = config->GetRotating_Frame();
	
	Normal = new double[nDim];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set conservative variables at the wall, and at the infinity ---*/
		for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
			FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);

		FlowSolution_j[0] = solution_container[FLOW_SOL]->GetDensity_Inf(); 
		FlowSolution_j[nDim+1] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
		for (iDim = 0; iDim < nDim; iDim++)
			FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(iDim);
		
		/*--- Rotational frame ---*/
		if (rotating_frame)
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());

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
		solver->SetResidual(Residual, Jacobian_i, NULL, config);
		
		/*--- Add residuals and jacobians ---*/
		node[iPoint]->AddResidual(Residual);
		Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
	}
	
	delete [] Normal; 
}

void CTurbSASolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
									  unsigned short val_marker) {		
	unsigned long iPoint, iVertex;
	unsigned short iVar, iDim;
	double Pressure_Total, Temperature_Total, Velocity[3], Velocity2, Enthalpy_Total, 
	Enthalpy_Static, Temperature_Static, Pressure_Static, Density, Energy, *Flow_Direction;
  double Gas_Constant = config->GetGas_Constant();
  double *Normal = new double[3];
  
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- FlowSolution_i -> U_internal ---*/
		for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
			FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
    
    /*--- Retrieve the specified flow quantities for this inlet boundary. ---*/
    Pressure_Total    = config->GetInlet_Ptotal(config->GetMarker_All_Tag(val_marker));
    Temperature_Total = config->GetInlet_Ttotal(config->GetMarker_All_Tag(val_marker));
    Flow_Direction    = config->GetInlet_FlowDir(config->GetMarker_All_Tag(val_marker));
    
    /*--- Interpolate the velocity from the interior of the domain.  ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity[iDim] = solution_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim);
    Velocity2 = solution_container[FLOW_SOL]->node[iPoint]->GetVelocity2();
    
    /*--- Correct the flow velocity to be in the specified direction at the inlet. ---*/
    Velocity[0] = Flow_Direction[0]*sqrt(Velocity2); 
    Velocity[1] = Flow_Direction[1]*sqrt(Velocity2); 
    if (nDim == 3) Velocity[2] = Flow_Direction[2]*sqrt(Velocity2); 
    
    /*--- With the total temperature we compute the total enthalpy ---*/
    Enthalpy_Total = (Gamma * Gas_Constant / Gamma_Minus_One) * Temperature_Total;
    
    /*--- If we substract the $0.5*velocity^2$ we obtain the static enthalpy (at the inlet) ---*/
    Enthalpy_Static = Enthalpy_Total - 0.5*Velocity2;
    
    /*--- With the static enthalpy (at the inlet) it is possible to compute the static temperature (at the inlet) ---*/
    Temperature_Static = Enthalpy_Static * Gamma_Minus_One / (Gamma * Gas_Constant);
    
    /*--- With the static temperature (at the inlet), and the total temperature (in the stagnation point), using the 
     isentropic relations it is possible to compute the static pressure at the inlet ---*/
    Pressure_Static = Pressure_Total * pow((Temperature_Static/Temperature_Total), Gamma/Gamma_Minus_One);
    
    /*--- Using the static pressure (at the inlet) and static temperature (at the inlet) we will compute 
     the density (at the inlet) ---*/
    Density = Pressure_Static / (Gas_Constant * Temperature_Static);
    
    /*--- Using pressure, density, and velocity we can compute the energy at the inlet ---*/
    Energy = Pressure_Static/(Density*Gamma_Minus_One)+0.5*Velocity2;
    
    /*--- Conservative variables, using the derived quatities ---*/
    FlowSolution_j[0] = Density;
    FlowSolution_j[1] = Velocity[0]*Density;
    FlowSolution_j[2] = Velocity[1]*Density;
    FlowSolution_j[3] = Energy*Density;
    if (nDim == 3) {
      FlowSolution_j[3] = Velocity[2]*Density;
      FlowSolution_j[4] = Energy*Density;
    }
		
		/*--- Set the conservative variables ---*/
		solver->SetConservative(FlowSolution_i,FlowSolution_j); 
		
    /*--- Set Normal (it is necessary to change the sign) ---*/
		geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
		for (unsigned short iDim = 0; iDim < nDim; iDim++)
			Normal[iDim] = -Normal[iDim]; 
		solver->SetNormal(Normal);
    
		/*--- Set the turbulent variable ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
		Solution_j[iVar]= nu_tilde_Inf;
		solver->SetTurbVar(Solution_i,Solution_j);
		
		/*--- Add Residual and Jacobians ---*/
		solver->SetResidual(Residual, Jacobian_i, NULL, config);
		node[iPoint]->AddResidual(Residual);
		Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
	}
  delete[] Normal;
  

}

void CTurbSASolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
									   CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar;
  double *Normal = new double[3];
  
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set the conservative variables, density & Velocity same as in 
		 the interior, don't need to modify pressure for the convection. ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), 
								solution_container[FLOW_SOL]->node[iPoint]->GetSolution()); 
		
		/*--- Set the turbulent variables. Here we use a Neumann BC such
     that the turbulent variable is copied from the interior of the 
     domain to the outlet before computing the redidual.
     Solution_i --> TurbVar_internal, 
		 Solution_j --> TurbVar_outlet ---*/
		for (iVar = 0; iVar < nVar; iVar++) 
			Solution_i[iVar] = node[iPoint]->GetSolution(iVar); 
		Solution_j[iVar]= node[iPoint]->GetSolution(iVar);
		solver->SetTurbVar(Solution_i,Solution_j);
    
    /*--- Set Normal (it is necessary to change the sign) ---*/
		geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
		for (unsigned short iDim = 0; iDim < nDim; iDim++)
			Normal[iDim] = -Normal[iDim]; 
		solver->SetNormal(Normal);
		
		/*--- Add Residual and Jacobians ---*/
		solver->SetResidual(Residual, Jacobian_i, NULL, config);
		node[iPoint]->AddResidual(Residual);
		Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
	}
  delete[] Normal;
}

CTurbSSTSolution::CTurbSSTSolution(void) : CTurbSolution() {}

CTurbSSTSolution::CTurbSSTSolution(CGeometry *geometry, CConfig *config) : CTurbSolution() {
	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double dull_val;
	ifstream restart_file;
	char *cstr;
	string text_line;
	bool restart = (config->GetRestart() || config->GetRestart_Flow());

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	node = new CVariable*[geometry->GetnPoint()];

	/*--- Dimension of the problem --> dependent of the turbulent model ---*/
	nVar = 2;

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
		Initialize_Jacobian_Structure(geometry, config);
		xsol = new double [geometry->GetnPoint()*nVar];
		rhs = new double [geometry->GetnPoint()*nVar];
	}

	/*--- Computation of gradients by least squares ---*/
	if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
		(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) {
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			cvector[iVar] = new double [nDim];
	}

	/* --- Initialize value for model constants --- */
	/*sigma_k1  = 0.85;
	sigma_k2  = 1.0;
	sigma_om1 = 0.5;
	sigma_om2 = 0.856;
	beta_1    = 0.075;
	beta_2    = 0.0828;
	beta_star = 0.09;
	a1        = 0.31;
	gamma_1   = beta_1/beta_star - sigma_om1*0.41*0.41/sqrt(beta_star);
	gamma_2   = beta_2/beta_star - sigma_om2*0.41*0.41/sqrt(beta_star);*/

	/*--- Flow infinity initialization stuff ---*/
	double rhoInf   = config->GetDensity_FreeStreamND();
	double *VelInf  = config->GetVelocity_FreeStreamND();
	double muLamInf = config->GetViscosity_FreeStreamND();

	double VelMag = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		VelMag += VelInf[iDim]*VelInf[iDim];
	VelMag = sqrt(VelMag);

  /*--- Warning: These have been initialized to zero because
   they were being used uninitialized in the next couple of lines. ---*/
	double Intensity = 0.0;
	double viscRatio = 0.0;

	kine_Inf  = 3.0/2.0*(VelMag*VelMag*Intensity*Intensity);
	omega_Inf = rhoInf*kine_Inf/(muLamInf*viscRatio);

	/*--- Restart the solution from file information ---*/
	if (!restart) {
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CTurbSSTVariable(kine_Inf,omega_Inf, nDim, nVar, config);
	}
	else {
		string mesh_filename = config->GetSolution_FlowFileName();
		cstr = new char [mesh_filename.size()+1];
		strcpy (cstr, mesh_filename.c_str());
		restart_file.open(cstr, ios::in);
		if (restart_file.fail()) {
			cout << "There is no turbulent restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}

		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
			if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
			node[iPoint] = new CTurbSSTVariable(Solution[0],0, nDim, nVar, config);
		}
		restart_file.close();
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

}

void CTurbSSTSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++){
		node[iPoint]->SetResidualZero();
	}
	Jacobian.SetValZero();

}

void CTurbSSTSolution::Postprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {

	double rho, mu, dist, omega, kine, vort, strMag, F2, muT;
	double a1 = 0.31;
	double zeta;
	unsigned long iPoint;
	
	// Compute mean flow gradients
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
		solution_container[FLOW_SOL]->SetPrimVar_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
		solution_container[FLOW_SOL]->SetPrimVar_Gradient_LS(geometry, config);

	// Compute turbulence variable gradients
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
		SetSolution_Gradient_GG(geometry);
	if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) || (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES))
		SetSolution_Gradient_LS(geometry, config);

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		// Compute vorticity and rate of strain magnitude
		solution_container[FLOW_SOL]->node[iPoint]->SetVorticity();
		solution_container[FLOW_SOL]->node[iPoint]->SetStrainMag();
		vort   = 0.0; //not yet implemented
		strMag = solution_container[FLOW_SOL]->node[iPoint]->GetStrainMag();

		// Compute blending functions and cross diffusion
		rho  = solution_container[FLOW_SOL]->node[iPoint]->GetDensity();
		mu   = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
		dist = geometry->node[iPoint]->GetWallDistance();

		node[iPoint]->SetBlendingFunc(mu,dist,rho);
		F2 = node[iPoint]->GetF2blending();

		// Compute the eddy viscosity
		kine  = node[iPoint]->GetSolution(0);
		omega = node[iPoint]->GetSolution(1);

		zeta = min(1.0/omega,a1/strMag*F2);
		muT = min(max(rho*kine*zeta,0.0),1.0);
		node[iPoint]->SetmuT(muT);
	}
}

void CTurbSSTSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {
	double *turb_var_i, *turb_var_j, *Limiter_i = NULL, *Limiter_j = NULL, *U_i, *U_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	bool high_order_diss = (config->GetKind_Upwind_Turb() == SCALAR_UPWIND_2ND);
	bool rotating_frame = config->GetRotating_Frame();

	if (high_order_diss) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solution_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry);
		if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
			(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) solution_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
		if (config->GetKind_SlopeLimit() == VENKATAKRISHNAN) SetSolution_Limiter(geometry, config);
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

		/*--- Rotational frame ---*/
		if (rotating_frame)
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());

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
			if (config->GetKind_SlopeLimit() != NONE) {
				Limiter_i = node[iPoint]->GetLimiter();
				Limiter_j = node[jPoint]->GetLimiter();
			}
			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				if (config->GetKind_SlopeLimit() == NONE) {
					Solution_i[iVar] = turb_var_i[iVar] + Project_Grad_i;
					Solution_j[iVar] = turb_var_j[iVar] + Project_Grad_j;
				}
				else {
					Solution_i[iVar] = turb_var_i[iVar] + Project_Grad_i*Limiter_i[iVar];
					Solution_j[iVar] = turb_var_j[iVar] + Project_Grad_j*Limiter_j[iVar];
				}
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

			/*--- Laminar Viscosity ---*/
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
																	solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
			/*--- Eddy Viscosity ---*/
			solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
																				solution_container[FLOW_SOL]->node[jPoint]->GetEddyViscosity());

			/*--- Turbulent variables w/o reconstruction, and its gradients ---*/
			solver->SetTurbVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
			solver->SetTurbVarGradient(node[iPoint]->GetGradient(),
																 node[jPoint]->GetGradient());

			/*--- Menter's first blending function ---*/
			solver->SetF1blending(node[iPoint]->GetF1blending(),node[jPoint]->GetF1blending());

			/*--- Rate of strain magnitude ---*/
			solver->SetStrainMag(solution_container[FLOW_SOL]->node[iPoint]->GetStrainMag(),0.0);

			/*--- Cross diffusion ---*/
			solver->SetCrossDiff(node[iPoint]->GetCrossDiff(),0.0);

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

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

		/*--- Conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);

		/*--- Gradient of the primitive and conservative variables ---*/
		solver->SetPrimVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);

		/*--- Laminar viscosity ---*/
		solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);

		/*--- Eddy Viscosity ---*/
		solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(), 0.0);

		/*--- Turbulent variables w/o reconstruction, and its gradient ---*/
		solver->SetTurbVar(node[iPoint]->GetSolution(), NULL);
		solver->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);

		/*--- Set volume ---*/
		solver->SetVolume(geometry->node[iPoint]->GetVolume());

		/*--- Set distance to the surface ---*/
		solver->SetDistance(geometry->node[iPoint]->GetWallDistance(), 0.0);

		/*--- Menter's first blending function ---*/
		solver->SetF1blending(node[iPoint]->GetF1blending(),0.0);

		/*--- Compute the source term ---*/
		solver->SetResidual(Residual, Jacobian_i, NULL, config);

		/*--- Subtract residual and the jacobian ---*/
		node[iPoint]->SubtractResidual(Residual);
		Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
	}
}

void CTurbSSTSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar;

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Get the velocity vector ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Solution[iVar] = 0.0;

		node[iPoint]->SetSolution_Old(Solution);
		node[iPoint]->SetResidualZero();

		/*--- includes 1 in the diagonal ---*/
		Jacobian.DeleteValsRowi(iPoint);
	}
}

void CTurbSSTSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	double *Normal;
	unsigned short iVar, iDim;

	Normal = new double[nDim];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Set conservative variables at the wall, and at the infinity ---*/
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
//		Solution_j[0] = nu_tilde_Inf;
		solver->SetTurbVar(Solution_i, Solution_j);

		/*--- Set Normal (it is necessary to change the sign) ---*/
		geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
		for (iDim = 0; iDim < nDim; iDim++)
			Normal[iDim] = -Normal[iDim];
		solver->SetNormal(Normal);

		/*--- Compute residuals and jacobians ---*/
		solver->SetResidual(Residual, Jacobian_i, NULL, config);

		/*--- Add residuals and jacobians ---*/
		node[iPoint]->AddResidual(Residual);
		Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
	}
	
	delete [] Normal;
}

void CTurbSSTSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
									  unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar, iDim;
	double SoundSpeed_internal, Vel2_internal, Mach_internal, Mrel, P_internal, SoundSpeed_inlet, DensVel2;
	double AoA = (config->GetAoA()*PI_NUMBER) / 180.0;
	double AoS = (config->GetAoS()*PI_NUMBER) / 180.0;

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- FlowSolution_i -> U_internal ---*/
		for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
			FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);

		/*--- Compute the Mach number at the internal nodes of the inlet ---*/
		SoundSpeed_internal = solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
		Vel2_internal = 0;
		for ( iDim = 0; iDim < nDim; iDim++)
			Vel2_internal += (FlowSolution_i[iDim+1]*FlowSolution_i[iDim+1])/(FlowSolution_i[0]*FlowSolution_i[0]);
		Mach_internal = sqrt(Vel2_internal)/SoundSpeed_internal;

		/*--- Compute the relative Mach number (isentropic relation) ---*/
		Mrel = (1.0+(Gamma_Minus_One*config->GetMach_FreeStreamND()*config->GetMach_FreeStreamND())/2.0)/
		(1.+(Gamma_Minus_One*Mach_internal*Mach_internal)/2.0);

		/*--- Compute the pressure at the internal nodes of the inlet ---*/
		P_internal = pow(Mrel,Gamma/Gamma_Minus_One)/(Gamma*config->GetMach_FreeStreamND()*config->GetMach_FreeStreamND());

		/*--- Compute the SoundSpeed at the inlet ---*/
		SoundSpeed_inlet = (1.0/config->GetMach_FreeStreamND())*sqrt(Mrel);

		/*--- Compute the solution at the inlet, FlowSolution_j --> U_inlet ---*/
		FlowSolution_j[0] = Gamma*P_internal/(SoundSpeed_inlet*SoundSpeed_inlet);
		if (nDim == 2) {
			FlowSolution_j[1] = SoundSpeed_inlet*Mach_internal*FlowSolution_j[0]*cos(AoA);
			FlowSolution_j[2] = SoundSpeed_inlet*Mach_internal*FlowSolution_j[0]*sin(AoA);
		}
		if (nDim == 3) {
			FlowSolution_j[1] = SoundSpeed_inlet*Mach_internal*FlowSolution_j[0]*cos(AoA)*cos(AoS);
			FlowSolution_j[2] = SoundSpeed_inlet*Mach_internal*FlowSolution_j[0]*sin(AoS);
			FlowSolution_j[3] = SoundSpeed_inlet*Mach_internal*FlowSolution_j[0]*sin(AoA)*cos(AoS);
		}

		DensVel2 = 0;
		for (iDim = 0; iDim < nDim; iDim++)
			DensVel2 += 0.5*FlowSolution_j[iDim+1]*FlowSolution_j[iDim+1]/FlowSolution_j[0];
		FlowSolution_j[nDim+1] = P_internal/Gamma_Minus_One + DensVel2;

		/*--- Set the conservative variables ---*/
		solver->SetConservative(FlowSolution_i,FlowSolution_j);

		/*--- Set the turbulent variable ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
//		Solution_j[iVar]= nu_tilde_Inf;
		solver->SetTurbVar(Solution_i,Solution_j);

		/*--- Set normal vector ---*/
		solver->SetNormal(geometry->vertex[val_marker][iVertex]->GetNormal());

		/*--- Add Residual and Jacobians ---*/
		solver->SetResidual(Residual, Jacobian_i, NULL, config);
		node[iPoint]->AddResidual(Residual);
		Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
	}
}

void CTurbSSTSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
									   CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar;

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Set the conservative variables, density & Velocity same as in
		 the interior, don't need to modify pressure for the convection. ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
								solution_container[FLOW_SOL]->node[iPoint]->GetSolution());

		/*--- Set the turbulent variables, Solution_i --> TurbVar_internal,
		 Solution_j --> TurbVar_outlet ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
//		Solution_j[iVar]= nu_tilde_Inf;
		solver->SetTurbVar(Solution_i,Solution_j);

		/*--- Set normal vector ---*/
		solver->SetNormal(geometry->vertex[val_marker][iVertex]->GetNormal());

		/*--- Add Residual and Jacobians ---*/
		solver->SetResidual(Residual, Jacobian_i, NULL, config);
		node[iPoint]->AddResidual(Residual);
		Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
	}
}
