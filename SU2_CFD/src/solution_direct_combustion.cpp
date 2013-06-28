/*!
 * \file solution_direct_combustion.cpp
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

CCombustionSolution::CCombustionSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned short iVar, iDim;
	unsigned long iPoint;
	double Lambda_Ini;
	ifstream restart_file;
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
	Vector_i = new double[nDim]; Vector_j = new double[nDim];
	
	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
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
	
	/*--- Lambda_Ini at the initial ---*/
	Lambda_Ini = 0.0;
	
	if (!restart) {
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			node[iPoint] = new CCombustionVariable(Lambda_Ini, nDim, nVar, config);
		}
	}
}

CCombustionSolution::~CCombustionSolution(void) {
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

void CCombustionSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		node[iPoint]->Set_ResConv_Zero();
		node[iPoint]->Set_ResVisc_Zero();
		node[iPoint]->SetResidualZero();
	}
		
}

void CCombustionSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {
	double *LambdaComb_var_i, *LambdaComb_var_j, *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		
		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());
		
		/*--- Conservative variables w/o reconstruction ---*/
		U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
		U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
		solver->SetConservative(U_i, U_j);
		
		/*--- Lambda combustion variables w/o reconstruction ---*/
		LambdaComb_var_i = node[iPoint]->GetSolution();
		LambdaComb_var_j = node[jPoint]->GetSolution();
		solver->SetLambdaComb(LambdaComb_var_i[0],LambdaComb_var_j[0]);
		
		/*--- Add and subtract Residual ---*/
		solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
		node[iPoint]->AddResidual(Residual);
		node[jPoint]->SubtractResidual(Residual);
		
	}
}

void CCombustionSolution::SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
											  CConfig *config, unsigned short iMesh) {
	unsigned short iVar;
	unsigned long iPoint;
	double *U, *LambdaComb;	
	
	for (iVar = 0; iVar < nVar; iVar++)
		Residual[iVar] = 0;
	
	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 
		
		/*--- Set solution  ---*/
		U = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
		solver->SetConservative(U, U);
		
		/*--- Lambda combustion variables w/o reconstruction ---*/
		LambdaComb = node[iPoint]->GetSolution();
		solver->SetLambdaComb(LambdaComb[0], LambdaComb[0]);

		/*--- Set control volume ---*/
		solver->SetVolume(geometry->node[iPoint]->GetVolume());
		
		/*--- Compute Residual ---*/
		solver->SetResidual(Residual, config);
		
		/*--- Add Residual ---*/
		node[iPoint]->AddResidual(Residual);
	}
	
}


void CCombustionSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set the solution to the original value ---*/
		Solution[0] = node[iPoint]->GetSolution(0);
		
		node[iPoint]->SetSolution_Old(Solution);
		node[iPoint]->SetResidualZero();
		
	}
}

void CCombustionSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set the solution to the original value ---*/
		Solution[0] = node[iPoint]->GetSolution(0);
		
		node[iPoint]->SetSolution_Old(Solution);
		node[iPoint]->SetResidualZero();
		
	}
}

void CCombustionSolution::BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config,
																				unsigned short val_marker, unsigned short val_mesh) {
	
#ifndef NO_MPI
	unsigned short iVar;
	double *Combustion_Var, **Combustion_Grad;
	unsigned long iVertex, iPoint;
	short SendRecv = config->GetMarker_All_SendRecv(val_marker);
	
	/*--- Send information  ---*/
	if (SendRecv > 0) {
		
		double Buffer_Send_Combustion[geometry->nVertex[val_marker]][nVar];
		double Buffer_Send_Combustionx[geometry->nVertex[val_marker]][nVar];
		double Buffer_Send_Combustiony[geometry->nVertex[val_marker]][nVar];
		double Buffer_Send_Combustionz[geometry->nVertex[val_marker]][nVar];
		
		/*--- Dimensionalization ---*/
		unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
		int send_to = SendRecv;
		
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
			
			Combustion_Var = node[iPoint]->GetSolution();
			Combustion_Grad = node[iPoint]->GetGradient();
			
			for (iVar = 0; iVar < nVar; iVar++) {
				Buffer_Send_Combustion[iVertex][iVar] = Combustion_Var[iVar];
				Buffer_Send_Combustionx[iVertex][iVar] = Combustion_Grad[iVar][0];
				Buffer_Send_Combustiony[iVertex][iVar] = Combustion_Grad[iVar][1];
				if (nDim == 3) Buffer_Send_Combustionz[iVertex][iVar] = Combustion_Grad[iVar][2];
			}
		}
		
		MPI::COMM_WORLD.Bsend(&Buffer_Send_Combustion,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
		MPI::COMM_WORLD.Bsend(&Buffer_Send_Combustionx,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
		MPI::COMM_WORLD.Bsend(&Buffer_Send_Combustiony,nBuffer_Vector,MPI::DOUBLE,send_to, 2);
		if (nDim == 3) MPI::COMM_WORLD.Bsend(&Buffer_Send_Combustionz,nBuffer_Vector,MPI::DOUBLE,send_to, 3);
	}
	
	/*--- Receive information  ---*/
	if (SendRecv < 0) {
		
		double Buffer_Receive_Combustion[geometry->nVertex[val_marker]][nVar];
		double Buffer_Receive_Combustionx[geometry->nVertex[val_marker]][nVar];
		double Buffer_Receive_Combustiony[geometry->nVertex[val_marker]][nVar];
		double Buffer_Receive_Combustionz[geometry->nVertex[val_marker]][nVar];
		
		/*--- Dimensionalization ---*/
		unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
		int receive_from = abs(SendRecv);
		
		MPI::COMM_WORLD.Recv(&Buffer_Receive_Combustion,nBuffer_Vector,MPI::DOUBLE,receive_from, 0);
		MPI::COMM_WORLD.Recv(&Buffer_Receive_Combustionx,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
		MPI::COMM_WORLD.Recv(&Buffer_Receive_Combustiony,nBuffer_Vector,MPI::DOUBLE,receive_from, 2);
		if (nDim == 3) MPI::COMM_WORLD.Recv(&Buffer_Receive_Combustionz,nBuffer_Vector,MPI::DOUBLE,receive_from, 3);
		
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
			for (iVar = 0; iVar < nVar; iVar++) {
				node[iPoint]->SetSolution(iVar, Buffer_Receive_Combustion[iVertex][iVar]);
				node[iPoint]->SetGradient(iVar, 0, Buffer_Receive_Combustionx[iVertex][iVar]);
				node[iPoint]->SetGradient(iVar, 1, Buffer_Receive_Combustiony[iVertex][iVar]);
				if (nDim == 3) node[iPoint]->SetGradient(iVar, 2, Buffer_Receive_Combustionz[iVertex][iVar]);
			}
		}
	}
#endif
}

void CCombustionSolution::BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config,
																					unsigned short val_marker, unsigned short val_mesh) {
	
#ifndef NO_MPI
	unsigned short iVar;
	unsigned long iVertex, iPoint;
	short SendRecv = config->GetMarker_All_SendRecv(val_marker);
	
	/*--- Receive information  ---*/
	if (SendRecv < 0) {
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
			node[iPoint]->SetResidualZero();
			for (iVar = 0; iVar < nVar; iVar++)
				Jacobian.DeleteValsRowi(iPoint);
		}
	}
#endif
}

void CCombustionSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																		 CConfig *config, unsigned short val_marker) { }

void CCombustionSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
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
	//	Jacobian.SGSSolution(rhs, xsol, 1e-9, 10, false);
	Jacobian.LU_SGSIteration(rhs, xsol);
	
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		for (iVar = 0; iVar < nVar; iVar++)
			node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);
	
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_Max(iVar, sqrt(GetRes_Max(iVar)));
	}
}

void CCombustionSolution::ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	
	double *Residual, Vol, Delta;
	unsigned short iVar;
	unsigned long iPoint;
	
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);
	
	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Delta = solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() / Vol;
		
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

void CCombustionSolution::RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, 
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
		Delta = solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() / Vol;
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
