/*!
 * \file solution_direct_levelset.cpp
 * \brief Main subrotuines for solving the level set problem.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
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

CLevelSetSolution::CLevelSetSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double dull_val;
	ifstream restart_file;
	char *cstr;
	string text_line;
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	
	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	node = new CVariable*[geometry->GetnPoint()];
	
	/*--- Dimension of the problem ---*/
	nVar = 1;
	
	/*--- Define some auxiliar vector related with the residual ---*/
	Residual = new double[nVar]; Residual_Max = new double[nVar];
	Residual_i = new double[nVar]; Residual_j = new double[nVar];
	
	/*--- Define some auxiliar vector related with the solution ---*/
	Solution = new double[nVar];
	Solution_i = new double[nVar]; Solution_j = new double[nVar];
	
	/*--- Define some auxiliar vector related with the geometry ---*/
	Vector = new double[nDim];
	Vector_i = new double[nDim]; Vector_j = new double[nDim];
	
	/*--- Define some auxiliar vector related with the flow solution ---*/
	FlowSolution_i = new double [nDim+2]; FlowSolution_j = new double [nDim+2];
	
	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT) {
		/*--- Point to point Jacobians ---*/
		Jacobian_i = new double* [nVar];
		Jacobian_j = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new double [nVar];
			Jacobian_j[iVar] = new double [nVar];
		}
		/*--- Initialization of the structure of the whole Jacobian ---*/
		InitializeJacobianStructure(geometry, config);
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
	
	/*--- levelset_Inf at the infinity ---*/
	levelset_Inf = 0.0;
	
	/*--- Restart the solution from file information ---*/
	if (!restart) {
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			double levelset = geometry->node[iPoint]->GetCoord(1) - config->GetLevelSet_Zero();
			node[iPoint] = new CLevelSetVariable(levelset, nDim, nVar, config);
		}
	}
	else {
		string mesh_filename = config->GetSolution_FlowFileName();
		cstr = new char [mesh_filename.size()+1];
		strcpy (cstr, mesh_filename.c_str());
		restart_file.open(cstr, ios::in);
		if (restart_file.fail()) {
			cout << "There is no level set restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}
		
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0];
			if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
			node[iPoint] = new CLevelSetVariable(Solution[0], nDim, nVar, config);
		}
		restart_file.close();
	}
}

CLevelSetSolution::~CLevelSetSolution(void) {
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

void CLevelSetSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		node[iPoint]->Set_ResConv_Zero();
		node[iPoint]->SetResidualZero();
	}
	
	/*--- Implicit part ---*/
	if (implicit)
		Jacobian.SetValZero();
	
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
}

void CLevelSetSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {
	double *LevelSet_var_i, *LevelSet_var_j, *Limiter_i = NULL, *Limiter_j = NULL, *U_i, *U_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	bool high_order_diss = (config->GetKind_Upwind_LevelSet() == SCALAR_UPWIND_2ND);
	
	if (high_order_diss) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
			solution_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry);
			SetSolution_Gradient_GG(geometry);
		}
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
			solution_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
			SetSolution_Gradient_LS(geometry, config);
		}
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
		
		/*--- Level Set variables w/o reconstruction ---*/
		LevelSet_var_i = node[iPoint]->GetSolution();
		LevelSet_var_j = node[jPoint]->GetSolution();
		solver->SetLevelSetVar(LevelSet_var_i,LevelSet_var_j);
				
		if (high_order_diss) {

			for (iDim = 0; iDim < nDim; iDim++) {
				Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
				Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}
			
			/*--- Conservative solution using gradient reconstruction ---*/
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
			
			/*--- Level Set variables using gradient reconstruction ---*/
			Gradient_i = node[iPoint]->GetGradient();
			Gradient_j = node[jPoint]->GetGradient();
			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				if (config->GetKind_SlopeLimit_LevelSet() == NONE) {
					Solution_i[iVar] = LevelSet_var_i[iVar] + Project_Grad_i;
					Solution_j[iVar] = LevelSet_var_j[iVar] + Project_Grad_j;
				}
				else {
					Solution_i[iVar] = LevelSet_var_i[iVar] + Project_Grad_i*Limiter_i[iVar];
					Solution_j[iVar] = LevelSet_var_j[iVar] + Project_Grad_j*Limiter_j[iVar];
				}
			}
			solver->SetLevelSetVar(Solution_i, Solution_j);
		}
		
		/*--- Add and subtract Residual ---*/
		solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
		node[iPoint]->AddRes_Conv(Residual);
		node[jPoint]->SubtractRes_Conv(Residual);
		
		/*--- Implicit part ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
			Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
			Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);
		}
	}
}

void CLevelSetSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar;
	
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set the solution to the original value ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Solution[iVar] = node[iPoint]->GetSolution(0);
		
		node[iPoint]->SetSolution_Old(Solution);
		node[iPoint]->SetResidualZero();
		node[iPoint]->Set_ResConv_Zero();

		/*--- Includes 1 in the diagonal ---*/
		if (implicit)
			Jacobian.DeleteValsRowi(iPoint);
	}
}

void CLevelSetSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar;
	
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set the solution to the original value ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Solution[iVar] = node[iPoint]->GetSolution(0);
		
		node[iPoint]->SetSolution_Old(Solution);
		node[iPoint]->SetResidualZero();
		node[iPoint]->Set_ResConv_Zero();

		/*--- Includes 1 in the diagonal ---*/
		if (implicit)
			Jacobian.DeleteValsRowi(iPoint);
	}
}

void CLevelSetSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar, iDim;
	double *Normal;
	
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	
	Normal = new double[nDim];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set flow variables at the wall ---*/
		for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
			FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
		
		/*--- Set flow variables at the infinity ---*/
		FlowSolution_j[0] = solution_container[FLOW_SOL]->GetPressure_Inf();
		
		/*--- For the water, the velocity is defined at the infinity ---*/
		for (iDim = 0; iDim < nDim; iDim++)
			FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetVelocity_Inf(iDim);				
		solver->SetConservative(FlowSolution_i, FlowSolution_j); 
		
		/*--- Set level set variable at the wall, and at infinity, basically we maintain the distance ---*/
		for (iVar = 0; iVar < nVar; iVar++) 
			Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
		Solution_j[0] = geometry->node[iPoint]->GetCoord(1) - config->GetLevelSet_Zero();
		solver->SetLevelSetVar(Solution_i, Solution_j);
		
		/*--- Set the normal vector ---*/
		geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
		for (iDim = 0; iDim < nDim; iDim++)
			Normal[iDim] = -Normal[iDim]; 
		solver->SetNormal(Normal);
		
		/*--- Compute residuals and jacobians ---*/
		solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
		
		/*--- Add residuals and jacobians ---*/
		node[iPoint]->AddRes_Conv(Residual);
		if (implicit)
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
	}
	
	delete [] Normal; 
}

void CLevelSetSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar;

	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set the solution to the original value ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Solution[iVar] = node[iPoint]->GetSolution(0);
		
		node[iPoint]->SetSolution_Old(Solution);
		node[iPoint]->SetResidualZero();
		node[iPoint]->Set_ResConv_Zero();
		
		/*--- Includes 1 in the diagonal ---*/
		if (implicit)
			Jacobian.DeleteValsRowi(iPoint);
	}
}

void CLevelSetSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar;
	
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set the solution to the original value ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Solution[iVar] = node[iPoint]->GetSolution(0);
		
		node[iPoint]->SetSolution_Old(Solution);
		node[iPoint]->SetResidualZero();
		node[iPoint]->Set_ResConv_Zero();
		
		/*--- Includes 1 in the diagonal ---*/
		if (implicit)
			Jacobian.DeleteValsRowi(iPoint);
	}
}

void CLevelSetSolution::BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config,
																				unsigned short val_marker, unsigned short val_mesh) {
	
#ifndef NO_MPI
	unsigned short iVar;
	double *LevelSet_Var, **LevelSet_Grad;
	unsigned long iVertex, iPoint;
	short SendRecv = config->GetMarker_All_SendRecv(val_marker);
	
	/*--- Send information  ---*/
	if (SendRecv > 0) {
		
		double Buffer_Send_LevelSet[geometry->nVertex[val_marker]][nVar];
//		double Buffer_Send_LevelSetx[geometry->nVertex[val_marker]][nVar];
//		double Buffer_Send_LevelSety[geometry->nVertex[val_marker]][nVar];
//		double Buffer_Send_LevelSetz[geometry->nVertex[val_marker]][nVar];
		
		/*--- Dimensionalization ---*/
		unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
		int send_to = SendRecv;
		
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
			
			LevelSet_Var = node[iPoint]->GetSolution();
			LevelSet_Grad = node[iPoint]->GetGradient();
			
			for (iVar = 0; iVar < nVar; iVar++) {
				Buffer_Send_LevelSet[iVertex][iVar] = LevelSet_Var[iVar];
//				Buffer_Send_LevelSetx[iVertex][iVar] = LevelSet_Grad[iVar][0];
//				Buffer_Send_LevelSety[iVertex][iVar] = LevelSet_Grad[iVar][1];
//				if (nDim == 3) Buffer_Send_LevelSetz[iVertex][iVar] = LevelSet_Grad[iVar][2];
			}
		}
		
		MPI::COMM_WORLD.Bsend(&Buffer_Send_LevelSet,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
//		MPI::COMM_WORLD.Bsend(&Buffer_Send_LevelSetx,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
//		MPI::COMM_WORLD.Bsend(&Buffer_Send_LevelSety,nBuffer_Vector,MPI::DOUBLE,send_to, 2);
//		if (nDim == 3) MPI::COMM_WORLD.Bsend(&Buffer_Send_LevelSetz,nBuffer_Vector,MPI::DOUBLE,send_to, 3);
	}
	
	/*--- Receive information  ---*/
	if (SendRecv < 0) {
		
		double Buffer_Receive_LevelSet[geometry->nVertex[val_marker]][nVar];
//		double Buffer_Receive_LevelSetx[geometry->nVertex[val_marker]][nVar];
//		double Buffer_Receive_LevelSety[geometry->nVertex[val_marker]][nVar];
//		double Buffer_Receive_LevelSetz[geometry->nVertex[val_marker]][nVar];
		
		/*--- Dimensionalization ---*/
		unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
		int receive_from = abs(SendRecv);
		
		MPI::COMM_WORLD.Recv(&Buffer_Receive_LevelSet,nBuffer_Vector,MPI::DOUBLE,receive_from, 0);
//		MPI::COMM_WORLD.Recv(&Buffer_Receive_LevelSetx,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
//		MPI::COMM_WORLD.Recv(&Buffer_Receive_LevelSety,nBuffer_Vector,MPI::DOUBLE,receive_from, 2);
//		if (nDim == 3) MPI::COMM_WORLD.Recv(&Buffer_Receive_LevelSetz,nBuffer_Vector,MPI::DOUBLE,receive_from, 3);
		
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
			for (iVar = 0; iVar < nVar; iVar++) {
				node[iPoint]->SetSolution(iVar, Buffer_Receive_LevelSet[iVertex][iVar]);
//				node[iPoint]->SetGradient(iVar, 0, Buffer_Receive_LevelSetx[iVertex][iVar]);
//				node[iPoint]->SetGradient(iVar, 1, Buffer_Receive_LevelSety[iVertex][iVar]);
//				if (nDim == 3) node[iPoint]->SetGradient(iVar, 2, Buffer_Receive_LevelSetz[iVertex][iVar]);
			}
		}
	}
#endif
}

void CLevelSetSolution::BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config,
																					unsigned short val_marker, unsigned short val_mesh) {
	
#ifndef NO_MPI
	unsigned short iVar;
	unsigned long iVertex, iPoint;

	short SendRecv = config->GetMarker_All_SendRecv(val_marker);
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);

	
	/*--- Receive information  ---*/
	if (SendRecv < 0) {
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
			node[iPoint]->SetResidualZero();
			if (implicit) {
				for (iVar = 0; iVar < nVar; iVar++)
					Jacobian.DeleteValsRowi(iPoint);
			}
		}
	}
#endif
}

void CLevelSetSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																		 CConfig *config, unsigned short val_marker) { }

void CLevelSetSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta = 0.0, *local_Residual, Vol;
	
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar,0.0);
	
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		local_Residual = node[iPoint]->GetResidual();
		Vol = geometry->node[iPoint]->GetVolume();
		
		/*--- Modify matrix diagonal to assure diagonal dominance ---*/
		if (config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) {
			/*--- If the flow solver is doing a dual time stepping, the time step 
			 is the unsteady time step (be careful, this could be very large, and the 
			 iteration can be unstable ---*/
			Delta = Vol/config->GetUnst_TimeStep();
		}
		if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
			/*--- For a normal unsteady simulatio we use the global time step ---*/
			Delta = Vol/solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time();
		}
		
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
	if ((config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) || 
			(config->GetUnsteady_Simulation() == TIME_STEPPING)) {
		/*--- If it is a unsteady problem, the linear system must be converged ---*/
		Jacobian.SGSSolution(rhs, xsol, 1e-9, 99999, true);
	}
	else {
		Jacobian.LU_SGSIteration(rhs, xsol);
	}
	
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		for (iVar = 0; iVar < nVar; iVar++)
			node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);
	
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_Max(iVar, sqrt(GetRes_Max(iVar)));
	}
}

void CLevelSetSolution::ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	double *Residual, Vol, Delta, *Coord, LevelSetZero;
	unsigned short iVar;
	unsigned long iPoint;
	
	LevelSetZero = config->GetLevelSet_Zero();
	
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);
	
	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Coord = geometry->node[iPoint]->GetCoord();

		if (config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) {
			/*--- If the flow solver is doing a dual time stepping, the time step 
			 is the unsteady time step (be careful, this could be very large, and the 
			 iteration can be unstable... CFL criteria ---*/
			Delta = config->GetUnst_TimeStep() / Vol;
		}
		else {
			/*--- For a normal unsteady simulatio we use the global time step ---*/
			Delta = solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() / Vol;
		}
		
		Residual = node[iPoint]->GetResidual();
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, -Residual[iVar]*Delta);
			AddRes_Max( iVar, Residual[iVar]*Residual[iVar]*Vol );	
		}
	}
	
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, sqrt(GetRes_Max(iVar)) );
}

void CLevelSetSolution::RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, 
																					CConfig *config, unsigned short iRKStep) {
	double *Residual, Vol, Delta;
	unsigned short iVar;
	unsigned long iPoint;
	double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
	
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);
	
	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();

		if (config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) {
			/*--- If the flow solver is doing a dual time stepping, the time step 
			 is the unsteady time step (be careful, this could be very large, and the 
			 iteration can be unstable ---*/
			Delta = config->GetUnst_TimeStep() / Vol;
		}
		else {
			/*--- For a normal unsteady simulatio we use the global time step ---*/
			Delta = solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() / Vol;
		}
		
		Residual = node[iPoint]->GetResidual();
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, -Residual[iVar]*Delta*RK_AlphaCoeff);
			AddRes_Max( iVar, Residual[iVar]*Residual[iVar]*Vol );
		}
	}
	
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, sqrt(GetRes_Max(iVar)) );
}
