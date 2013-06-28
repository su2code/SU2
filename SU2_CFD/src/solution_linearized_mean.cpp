/*!
 * \file solution_linearized_mean.cpp
 * \brief Main subrotuines for solving linearized problems (Euler, Navier-Stokes, etc.).
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

CLinEulerSolution::CLinEulerSolution(void) : CSolution() { }

CLinEulerSolution::CLinEulerSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned long iPoint, index;
	string text_line, mesh_filename;
	unsigned short iDim, iVar;
	ifstream restart_file;
	char *cstr;
	bool restart = config->GetRestart();
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	
	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	nVar = geometry->GetnDim()+2;
	node = new CVariable*[geometry->GetnPoint()];
	
	/*--- Define some auxiliar vector related with the residual ---*/
	Residual = new double[nVar];	Residual_Max = new double[nVar];
	Residual_i = new double[nVar];	Residual_j = new double[nVar];
	Res_Conv = new double[nVar];	Res_Visc = new double[nVar];
	
	/*--- Define some auxiliar vector related with the solution ---*/
	Solution   = new double[nVar];
	Solution_i = new double[nVar];	Solution_j = new double[nVar];
	
	/*--- Define some auxiliar vector related with the geometry ---*/
	Vector_i = new double[nDim];	Vector_j = new double[nDim];
	
	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT) {
		/*--- Point to point Jacobians ---*/
		Jacobian_i = new double* [nVar]; Jacobian_j = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new double [nVar]; Jacobian_j[iVar] = new double [nVar]; }
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
	
	/*--- Delta forces definition and coefficient in all the markers ---*/
	DeltaForceInviscid = new double[nDim];
	CDeltaDrag_Inv = new double[config->GetnMarker_All()];
	CDeltaLift_Inv = new double[config->GetnMarker_All()];
	
	/*--- Linearized flow at the inifinity, inizialization stuff ---*/
	DeltaRho_Inf = 0.0;	
	DeltaE_Inf = 0.0;
	DeltaVel_Inf = new double [nDim];
	if (nDim == 2) {
		DeltaVel_Inf[0] = 0.0;	
		DeltaVel_Inf[1] = 0.0;
	}
	if (nDim == 3) {
		DeltaVel_Inf[0] = 0.0;
		DeltaVel_Inf[1] = 0.0;
		DeltaVel_Inf[2] = 0.0;		
	}
	
	/*--- Restart the solution from file information ---*/
	if (!restart) {
		for (iPoint=0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CLinEulerVariable(DeltaRho_Inf, DeltaVel_Inf, DeltaE_Inf, nDim, nVar, config);
	}
	else {
		mesh_filename = config->GetSolution_LinFileName();
		cstr = new char [mesh_filename.size()+1];
		strcpy (cstr, mesh_filename.c_str());
		restart_file.open(cstr, ios::in);
		if (restart_file.fail()) {
			cout << "There is no linearized restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}
		
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
			if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
			node[iPoint] = new CLinEulerVariable(Solution, nDim, nVar, config);
		}
		restart_file.close();
	}
}

CLinEulerSolution::~CLinEulerSolution(void) {
	unsigned short iVar, iDim;
	
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Jacobian_i[iVar];
		delete [] Jacobian_j[iVar];
	}
	delete [] Jacobian_i; delete [] Jacobian_j;
	
	delete [] Residual; delete [] Residual_Max;
	delete [] Residual_i; delete [] Residual_j;
	delete [] Res_Conv_i; delete [] Res_Visc_i;
	delete [] Res_Conv_j; delete [] Res_Visc_j;
	delete [] Solution; 
	delete [] Solution_i; delete [] Solution_j;
	delete [] Vector_i; delete [] Vector_j;
	delete [] xsol; delete [] rhs;
	delete [] DeltaForceInviscid; delete [] DeltaVel_Inf;
	delete [] CDeltaDrag_Inv; delete [] CDeltaLift_Inv;
	
	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;
	
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
	
	/*	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
	 delete [] node[iPoint];
	 delete [] node; */
}

void CLinEulerSolution::Centred_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
											 CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint;
	bool implicit = (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT);
	bool dissipation = ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit);
	bool high_order_diss = ((config->GetKind_Centred() == JST) && (iMesh == MESH_0));
	
	if (dissipation && high_order_diss) 
		SetUndivided_Laplacian(geometry, config);
	
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		/*--- Points in edge, normal, and neighbors---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());
		solver->SetNeighbor(geometry->node[iPoint]->GetnPoint(), geometry->node[jPoint]->GetnPoint());

		/*--- Linearized variables w/o reconstruction ---*/
		solver->SetLinearizedVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());

		/*--- Conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), 
								solution_container[FLOW_SOL]->node[jPoint]->GetSolution());

		/*--- SoundSpeed enthalpy and lambda variables w/o reconstruction ---*/
		solver->SetSoundSpeed(solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
							  solution_container[FLOW_SOL]->node[jPoint]->GetSoundSpeed());
		solver->SetEnthalpy(solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
							solution_container[FLOW_SOL]->node[jPoint]->GetEnthalpy());
		solver->SetLambda(solution_container[FLOW_SOL]->node[iPoint]->GetLambda(), 
						  solution_container[FLOW_SOL]->node[jPoint]->GetLambda());

		/*--- Undivided laplacian ---*/
		if (dissipation && high_order_diss) 
			solver->SetUndivided_Laplacian(node[iPoint]->GetUnd_Lapl(),node[jPoint]->GetUnd_Lapl());
		
		/*--- Compute residual ---*/
		solver->SetResidual(Res_Conv, Res_Visc, Jacobian_i, Jacobian_j, dissipation, config);
		
		/*--- Update convective and artificial dissipation residuals ---*/
		node[iPoint]->AddRes_Conv(Res_Conv);
		node[jPoint]->SubtractRes_Conv(Res_Conv);
		if (dissipation) {
			node[iPoint]->AddRes_Visc(Res_Visc);
			node[jPoint]->SubtractRes_Visc(Res_Visc);
		}

		/*--- Set implicit stuff ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
			Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
			Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
		}
	}
}

void CLinEulerSolution::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iVar;
	double *Diff = new double[nVar];
	bool boundary_i, boundary_j;
		
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
			node[iPoint]->SetUnd_LaplZero();
	
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		
		for (iVar = 0; iVar < nVar; iVar++)
			Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);
				
		boundary_i = geometry->node[iPoint]->GetBoundary_Physical();
		boundary_j = geometry->node[jPoint]->GetBoundary_Physical();
		
		/*--- Both points inside Omega ---*/
		if (!boundary_i && !boundary_j) {
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
		}
		
		/*--- iPoint inside Omega, jPoint on the boundary ---*/
		if (!boundary_i && boundary_j)
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
		
		/*--- jPoint inside Omega, iPoint on the boundary ---*/
		if (boundary_i && !boundary_j)
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
		
		/*--- Both points on the boundary ---*/
		if (boundary_i && boundary_j) {
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
		}
	}
	delete [] Diff;
}

void CLinEulerSolution::RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, 
											CConfig *config, unsigned short iRKStep) {
	double *Residual, Vol, Delta, *TruncationError;
	unsigned short iVar;
	unsigned long iPoint;
	double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);

	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, 0.0 );
	
	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			Vol = geometry->node[iPoint]->GetVolume();
			Delta = solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() / Vol;
			TruncationError = node[iPoint]->GetTruncationError();
			Residual = node[iPoint]->GetResidual();
			for (iVar = 0; iVar < nVar; iVar++) {
				node[iPoint]->AddSolution(iVar, -(Residual[iVar]+TruncationError[iVar])*Delta*RK_AlphaCoeff);
				AddRes_Max( iVar, Residual[iVar]*Residual[iVar]*Vol );
			}
		}
	
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, sqrt(GetRes_Max(iVar)) );
}

void CLinEulerSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
									  unsigned short iRKStep) {
	unsigned long iPoint;
	bool implicit = (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT);
	
	/*--- Residual inicialization ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		node[iPoint]->Set_ResConv_Zero();
		if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit)
			node[iPoint]->Set_ResVisc_Zero();
	}
	
	/*--- Inicialize the jacobian matrices ---*/
	if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT)
		Jacobian.SetValZero();
}

void CLinEulerSolution::Inviscid_DeltaForces(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned long iVertex, Point;
	unsigned short iDim, iMarker, Boundary, Monitoring;
	double  *Face_Normal, dS, DeltaPressure, *Velocity;
	double Alpha = config->GetAoA()*PI_NUMBER / 180.0;
	double Beta  = config->GetAoS()*PI_NUMBER / 180.0;
	double RefAreaCoeff = config->GetRefAreaCoeff();
	double Density_Inf = solution_container[FLOW_SOL]->GetDensity_Inf();
	double ModVelocity_Inf = solution_container[FLOW_SOL]->GetModVelocity_Inf();
	double C_p = 1.0/(0.5*Density_Inf*RefAreaCoeff*ModVelocity_Inf*ModVelocity_Inf);
	
	/*-- Inicialization ---*/
	Total_CDeltaDrag = 0.0; Total_CDeltaLift = 0.0;
	AllBound_CDeltaDrag_Inv = 0.0; AllBound_CDeltaLift_Inv = 0.0;
	Velocity = new double[nDim];

	/*--- Loop over the Euler markers ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);
		if ((Boundary == EULER_WALL) || (Boundary == NO_SLIP_WALL)) {
			for (iDim = 0; iDim < nDim; iDim++) DeltaForceInviscid[iDim] = 0.0;
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				Point = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[Point]->GetDomain()) {
					
					/*--- Compute pressure on the boundary ---*/
					for (iDim = 0; iDim < nDim; iDim++) 
						Velocity[iDim] = solution_container[FLOW_SOL]->node[Point]->GetVelocity(iDim);
					
					double rho = solution_container[FLOW_SOL]->node[Point]->GetSolution(0) + node[Point]->GetSolution(0);
					double rhoE = solution_container[FLOW_SOL]->node[Point]->GetSolution(nVar-1) + node[Point]->GetSolution(nVar-1);
					double Pressure = solution_container[FLOW_SOL]->node[Point]->GetPressure();
					double rhoVel[3];
					double sqr_vel = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						rhoVel[iDim] = solution_container[FLOW_SOL]->node[Point]->GetSolution(iDim+1) + node[Point]->GetSolution(iDim+1);
						sqr_vel += rhoVel[iDim]*rhoVel[iDim]/(rho*rho);
					}
					DeltaPressure = Gamma_Minus_One*rho*(rhoE/rho-0.5*sqr_vel)-Pressure;
					
					if (Monitoring == YES) {
						Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
						dS = 0.0; for (iDim = 0; iDim < nDim; iDim++) dS += Face_Normal[iDim]*Face_Normal[iDim]; dS = sqrt(dS);
						for (iDim = 0; iDim < nDim; iDim++)
							DeltaForceInviscid[iDim] -= C_p*DeltaPressure*Face_Normal[iDim];
					}
				}
			}
			
			/*--- Transform ForceInviscid into CLift and CDrag ---*/
			if  (Monitoring == YES) {
				if (nDim == 2) {
					CDeltaDrag_Inv[iMarker] =  DeltaForceInviscid[0]*cos(Alpha) + DeltaForceInviscid[1]*sin(Alpha);
					CDeltaLift_Inv[iMarker] = -DeltaForceInviscid[0]*sin(Alpha) + DeltaForceInviscid[1]*cos(Alpha);
				}
				if (nDim == 3) {
					CDeltaDrag_Inv[iMarker] =  DeltaForceInviscid[0]*cos(Alpha)*cos(Beta) + DeltaForceInviscid[1]*sin(Beta) + DeltaForceInviscid[2]*sin(Alpha)*cos(Beta);
					CDeltaLift_Inv[iMarker] = -DeltaForceInviscid[0]*sin(Alpha) + DeltaForceInviscid[2]*cos(Alpha);
				}
				
				AllBound_CDeltaDrag_Inv += CDeltaDrag_Inv[iMarker];
				AllBound_CDeltaLift_Inv += CDeltaLift_Inv[iMarker];
			}
		}
	}
	
	Total_CDeltaDrag += AllBound_CDeltaDrag_Inv;
	Total_CDeltaLift += AllBound_CDeltaLift_Inv;
	
	delete [] Velocity;
}


void CLinEulerSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {
	unsigned short iDim, iVar;
	unsigned long iVertex, iPoint;
	double phi, a1, a2, sqvel, d, *Face_Normal, *U, dS, VelxDeltaRhoVel, Energy, Rho, 
		   *Normal, *Delta_RhoVel, *Vel, *DeltaU, Delta_RhoE, Delta_Rho;

	Normal = new double[nDim];
	Delta_RhoVel = new double[nDim];
	Vel = new double[nDim];
	
	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		/*--- Point, normal vector, and area ---*/
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			
			Face_Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			dS = 0; for (iDim = 0; iDim < nDim; iDim++) dS += Face_Normal[iDim]*Face_Normal[iDim]; dS = sqrt(dS);
			
			/*--- Linearized and flow solution ---*/
			DeltaU = node[iPoint]->GetSolution();
			U = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
			
			/*--- Value of the linearized velocity projection ---*/
			d = 0.0;
			
			Rho = U[0];
			Energy = U[nDim+1] / Rho;
			Delta_Rho = DeltaU[0];
			Delta_RhoE = DeltaU[nDim+1];
			sqvel = 0.0; VelxDeltaRhoVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Normal[iDim] = -Face_Normal[iDim]/dS;
				Vel[iDim] = U[iDim+1] / Rho;
				Delta_RhoVel[iDim] = DeltaU[iDim+1];
				sqvel += Vel[iDim]*Vel[iDim];
				VelxDeltaRhoVel += Vel[iDim]*Delta_RhoVel[iDim];
			}  
			phi = 0.5*Gamma_Minus_One*sqvel;
			a1 = Gamma*Energy-phi; a2 = Gamma-1.0;
			
			Residual[0] = Rho * d;
			for (iDim = 0; iDim < nDim; iDim++)
				Residual[iDim+1] = Normal[iDim]*(phi*Delta_Rho+a2*(Delta_RhoE-VelxDeltaRhoVel)) + Vel[iDim] * Rho * d;
			Residual[nVar-1] = a1 * Rho * d;
			
			for (iVar = 0; iVar < nVar; iVar++)
				Residual[iVar] = dS * Residual[iVar];
			
			/*--- Add value to the residual ---*/
			node[iPoint]->AddRes_Conv(Residual);
		}
	}
	
	delete [] Vel;
	delete [] Normal;
	delete [] Delta_RhoVel;
}

void CLinEulerSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iVar, jVar, iDim;
	
	double *Face_Normal;
	double *kappa = new double[nDim];
	double *velocity = new double[nDim];
	double *U_wall = new double[nVar];
	double *U_infty = new double[nVar];
	double *DeltaU_wall = new double[nVar];
	double *DeltaU_infty = new double[nVar];
	double *DeltaU_update = new double[nVar];
	double *W_wall = new double[nVar];
	double *W_infty = new double[nVar];
	double *W_update = new double[nVar];
	double Mach = config->GetMach_FreeStreamND();
	double DeltaMach = 1.0;

	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_wall[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			
			/*--- Flow Solution at the infinity ---*/
			U_infty[0] = solution_container[FLOW_SOL]->GetDensity_Inf();
			U_infty[1] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(0);
			U_infty[2] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(1);
			U_infty[3] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
			if (nDim == 3) {
				U_infty[3] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(2);
				U_infty[4] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
			}
			
			/*--- Delta flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				DeltaU_wall[iVar] = node[iPoint]->GetSolution(iVar);
			
			/*--- Linearized flow solution at the far-field ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				DeltaU_infty[iVar] = 0.0;			
			DeltaU_infty[nVar-1] = -2.0*DeltaMach/(Gamma*Gamma_Minus_One*Mach*Mach*Mach);			

			/*--- Normal vector ---*/
			Face_Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			double dS = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				dS += Face_Normal[iDim]*Face_Normal[iDim];
			dS = sqrt (dS);
			
			for (iDim = 0; iDim < nDim; iDim++)
				kappa[iDim] = -Face_Normal[iDim]/dS;
			
			/*--- Computation of P and inverse P matrix using values at the infinity ---*/
			double sq_vel = 0.0, vn = 0.0;
			double rho = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(0);
			double rhoE = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(nVar-1);
			for (iDim = 0; iDim < nDim; iDim++) {
				velocity[iDim] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iDim+1)/rho;
				sq_vel +=velocity[iDim]*velocity[iDim];
				vn += velocity[iDim]*kappa[iDim]*dS;
			}
			
/*			double sq_vel = 0.0, vn = 0.0;
			double rho = U_infty[0];
			double rhoE = U_infty[nVar-1];
			for (iDim = 0; iDim < nDim; iDim++) {
				velocity[iDim] = U_infty[iDim+1]/rho;
				sq_vel +=velocity[iDim]*velocity[iDim];
				vn += velocity[iDim]*kappa[iDim]*dS;
			}*/
			
			double c = sqrt(Gamma*Gamma_Minus_One*(rhoE/rho-0.5*sq_vel));
			double energy = rhoE / rho;
			
			double **P_Matrix, **invP_Matrix;
			P_Matrix = new double* [nVar];
			invP_Matrix = new double* [nVar];
			for (iVar = 0; iVar < nVar; iVar++) {
				P_Matrix[iVar] = new double [nVar];
				invP_Matrix[iVar] = new double [nVar];
			}
			
			solver->GetPMatrix_inv(&rho, velocity, &c, kappa, invP_Matrix);
			solver->GetPMatrix(&rho, velocity, &c, kappa, P_Matrix);
			
			/*--- computation of characteristics variables at the wall ---*/			
			for (iVar=0; iVar < nVar; iVar++) {
				W_wall[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_wall[iVar] +=invP_Matrix[iVar][jVar]*DeltaU_wall[jVar];
			}
			
			/*--- computation of characteristics variables at the far-field ---*/			
			for (iVar=0; iVar < nVar; iVar++) {
				W_infty[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_infty[iVar] +=invP_Matrix[iVar][jVar]*DeltaU_infty[jVar];
			}
			
			/*--- fix characteristics value ---*/			
			if (nDim == 2) {
				if(vn > 0.0) { 
					W_update[0] = W_wall[0];
					W_update[1] = W_wall[1];
				}
				else {
					W_update[0] = W_infty[0];
					W_update[1] = W_infty[1];
				}
				
				if(vn+c*dS > 0.0) W_update[2] = W_wall[2];
				else W_update[2] = W_infty[2];
				
				if(vn-c*dS > 0.0) W_update[3] = W_wall[3];
				else W_update[3] = W_infty[3];
			}
			
			if (nDim == 3) {
				if(vn > 0.0) { 
					W_update[0] = W_wall[0];
					W_update[1] = W_wall[1];
					W_update[2] = W_wall[2];
				}
				else {
					W_update[0] = W_infty[0];
					W_update[1] = W_infty[1];
					W_update[2] = W_infty[2];
				}
				
				if(vn+c*dS > 0.0) W_update[3] = W_wall[3];
				else W_update[3] = W_infty[3];
				
				if(vn-c*dS > 0.0) W_update[4] = W_wall[4];
				else W_update[4] = W_infty[4];
			}
			
			/*--- conservative variables using characteristics ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				DeltaU_update[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					DeltaU_update[iVar] +=P_Matrix[iVar][jVar]*W_update[jVar];
			}
			
			/*--- Residual computation ---*/
			double **Jac_Matrix;
			Jac_Matrix = new double* [nVar];
			for (iVar = 0; iVar < nVar; iVar++) {
				Jac_Matrix[iVar] = new double [nVar];
			}			
			solver->GetInviscidProjJac(velocity, energy, kappa, 1.0, Jac_Matrix);
			for (iVar = 0; iVar < nVar; iVar++) {
				Residual[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					Residual[iVar] +=Jac_Matrix[iVar][jVar]*DeltaU_update[jVar]*dS;
			}
			
			node[iPoint]->AddRes_Conv(Residual);
			
		}
	}
}

void CLinEulerSolution::BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
										unsigned short val_marker, unsigned short val_mesh) {
#ifndef NO_MPI
	unsigned short iVar;
	unsigned long iVertex, Point;
	double *Linearized_Var, *Linearized_Undivided_Laplacian = NULL, **Linearized_Grad = NULL;
	short SendRecv = config->GetMarker_All_SendRecv(val_marker);
	
	/*--- Send information  ---*/
	if (SendRecv > 0) {
		
		/*--- Upwind scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {
			double Buffer_Send_DeltaU[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_DeltaUx[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_DeltaUy[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_DeltaUz[geometry->nVertex[val_marker]][nVar];
			unsigned long nBuffer = geometry->nVertex[val_marker]*nVar;
			int send_to = SendRecv;
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				Point = geometry->vertex[val_marker][iVertex]->GetNode();
				Linearized_Var = node[Point]->GetSolution();
				if (val_mesh == MESH_0) Linearized_Grad = node[Point]->GetGradient();
				
				for (iVar = 0; iVar < nVar; iVar++) {
					Buffer_Send_DeltaU[iVertex][iVar] = Linearized_Var[iVar];
					if (val_mesh == MESH_0) {
						Buffer_Send_DeltaUx[iVertex][iVar] = Linearized_Grad[iVar][0];					
						Buffer_Send_DeltaUy[iVertex][iVar] = Linearized_Grad[iVar][1];
						if (nDim == 3) Buffer_Send_DeltaUz[iVertex][iVar] = Linearized_Grad[iVar][2];
					}
				}
			}
			
			MPI::COMM_WORLD.Bsend(&Buffer_Send_DeltaU,nBuffer,MPI::DOUBLE,send_to, 0);
			if (val_mesh == MESH_0) {
				MPI::COMM_WORLD.Bsend(&Buffer_Send_DeltaUx,nBuffer,MPI::DOUBLE,send_to, 1);
				MPI::COMM_WORLD.Bsend(&Buffer_Send_DeltaUy,nBuffer,MPI::DOUBLE,send_to, 2);
				if (nDim == 3) MPI::COMM_WORLD.Bsend(&Buffer_Send_DeltaUz,nBuffer,MPI::DOUBLE,send_to, 3);
			}
		}
		
		/*--- Centered scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {
			
			double Buffer_Send_DeltaU[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Undivided_Laplacian[geometry->nVertex[val_marker]][nVar];
			unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
			int send_to = SendRecv;
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				Point = geometry->vertex[val_marker][iVertex]->GetNode();
				Linearized_Var = node[Point]->GetSolution();
				if (val_mesh == MESH_0) Linearized_Undivided_Laplacian = node[Point]->GetUnd_Lapl();
				for (iVar = 0; iVar < nVar; iVar++) {
					Buffer_Send_DeltaU[iVertex][iVar] = Linearized_Var[iVar];
					if (val_mesh == MESH_0) Buffer_Send_Undivided_Laplacian[iVertex][iVar] = Linearized_Undivided_Laplacian[iVar];					
				}
			}
			
			MPI::COMM_WORLD.Bsend(&Buffer_Send_DeltaU,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Bsend(&Buffer_Send_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
			
		}
	}
	
	/*--- Receive information  ---*/
	if (SendRecv < 0) {
		/*--- Upwind scheme (Not validated)---*/
		if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {
			double Buffer_Receive_DeltaU[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_DeltaUx[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_DeltaUy[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_DeltaUz[geometry->nVertex[val_marker]][nVar];
			unsigned long nBuffer = geometry->nVertex[val_marker]*nVar;
			int receive_from = abs(SendRecv);
			
			MPI::COMM_WORLD.Recv(&Buffer_Receive_DeltaU,nBuffer,MPI::DOUBLE,receive_from, 0);
			if (val_mesh == MESH_0) {
				MPI::COMM_WORLD.Recv(&Buffer_Receive_DeltaUx,nBuffer,MPI::DOUBLE,receive_from, 1);
				MPI::COMM_WORLD.Recv(&Buffer_Receive_DeltaUy,nBuffer,MPI::DOUBLE,receive_from, 2);
				if (nDim == 3) MPI::COMM_WORLD.Recv(&Buffer_Receive_DeltaUz,nBuffer,MPI::DOUBLE,receive_from, 3);
			}
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				Point = geometry->vertex[val_marker][iVertex]->GetNode();	
				for (iVar = 0; iVar < nVar; iVar++) {
					node[Point]->SetSolution(iVar, Buffer_Receive_DeltaU[iVertex][iVar]);
					if (val_mesh == MESH_0) {
						node[Point]->SetGradient(iVar, 0, Buffer_Receive_DeltaUx[iVertex][iVar]);
						node[Point]->SetGradient(iVar, 1, Buffer_Receive_DeltaUy[iVertex][iVar]);
						if (nDim == 3) node[Point]->SetGradient(iVar, 2, Buffer_Receive_DeltaUz[iVertex][iVar]);
					}
					if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT)
						Jacobian.DeleteValsRowi(Point*nVar+iVar);
				}
				if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT) {
					node[Point]->Set_ResConv_Zero();
					node[Point]->Set_ResVisc_Zero();
					node[Point]->SetResidualZero();
					node[Point]->SetTruncationErrorZero();
				}
			}
		}
		
		/*--- Centered scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {
			double Buffer_Receive_DeltaU[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Undivided_Laplacian[geometry->nVertex[val_marker]][nVar];
			unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
			int receive_from = abs(SendRecv);
			
			MPI::COMM_WORLD.Recv(&Buffer_Receive_DeltaU,nBuffer_Vector,MPI::DOUBLE,receive_from, 0);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Recv(&Buffer_Receive_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				Point = geometry->vertex[val_marker][iVertex]->GetNode();
				for (iVar = 0; iVar < nVar; iVar++) {
					node[Point]->SetSolution(iVar, Buffer_Receive_DeltaU[iVertex][iVar]);
					if (val_mesh == MESH_0) node[Point]->SetUndivided_Laplacian(iVar, Buffer_Receive_Undivided_Laplacian[iVertex][iVar]);
					if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT)
						Jacobian.DeleteValsRowi(Point*nVar+iVar);
				}
				if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT) {
					node[Point]->Set_ResConv_Zero();
					node[Point]->Set_ResVisc_Zero();
					node[Point]->SetResidualZero();
					node[Point]->SetTruncationErrorZero();
				}
			}
		}			
		
	}
#endif
}

void CLinEulerSolution::BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
																				unsigned short val_marker, unsigned short val_mesh) {
#ifndef NO_MPI
	unsigned short iVar;
	unsigned long iVertex, Point;
	double *Linearized_Var, *Linearized_Undivided_Laplacian = NULL, **Linearized_Grad = NULL;
	short SendRecv = config->GetMarker_All_SendRecv(val_marker);
	
	/*--- Send information  ---*/
	if (SendRecv > 0) {
		
		/*--- Upwind scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {
			double Buffer_Send_DeltaU[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_DeltaUx[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_DeltaUy[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_DeltaUz[geometry->nVertex[val_marker]][nVar];
			unsigned long nBuffer = geometry->nVertex[val_marker]*nVar;
			int send_to = SendRecv;
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				Point = geometry->vertex[val_marker][iVertex]->GetNode();
				Linearized_Var = node[Point]->GetSolution();
				if (val_mesh == MESH_0) Linearized_Grad = node[Point]->GetGradient();
				
				for (iVar = 0; iVar < nVar; iVar++) {
					Buffer_Send_DeltaU[iVertex][iVar] = Linearized_Var[iVar];
					if (val_mesh == MESH_0) {
						Buffer_Send_DeltaUx[iVertex][iVar] = Linearized_Grad[iVar][0];					
						Buffer_Send_DeltaUy[iVertex][iVar] = Linearized_Grad[iVar][1];
						if (nDim == 3) Buffer_Send_DeltaUz[iVertex][iVar] = Linearized_Grad[iVar][2];
					}
				}
			}
			
			MPI::COMM_WORLD.Bsend(&Buffer_Send_DeltaU,nBuffer,MPI::DOUBLE,send_to, 0);
			if (val_mesh == MESH_0) {
				MPI::COMM_WORLD.Bsend(&Buffer_Send_DeltaUx,nBuffer,MPI::DOUBLE,send_to, 1);
				MPI::COMM_WORLD.Bsend(&Buffer_Send_DeltaUy,nBuffer,MPI::DOUBLE,send_to, 2);
				if (nDim == 3) MPI::COMM_WORLD.Bsend(&Buffer_Send_DeltaUz,nBuffer,MPI::DOUBLE,send_to, 3);
			}
		}
		
		/*--- Centered scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {
			
			double Buffer_Send_DeltaU[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Undivided_Laplacian[geometry->nVertex[val_marker]][nVar];
			unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
			int send_to = SendRecv;
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				Point = geometry->vertex[val_marker][iVertex]->GetNode();
				Linearized_Var = node[Point]->GetSolution();
				if (val_mesh == MESH_0) Linearized_Undivided_Laplacian = node[Point]->GetUnd_Lapl();
				for (iVar = 0; iVar < nVar; iVar++) {
					Buffer_Send_DeltaU[iVertex][iVar] = Linearized_Var[iVar];
					if (val_mesh == MESH_0) Buffer_Send_Undivided_Laplacian[iVertex][iVar] = Linearized_Undivided_Laplacian[iVar];					
				}
			}
			
			MPI::COMM_WORLD.Bsend(&Buffer_Send_DeltaU,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Bsend(&Buffer_Send_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
			
		}
	}
	
	/*--- Receive information  ---*/
	if (SendRecv < 0) {
		/*--- Upwind scheme (Not validated)---*/
		if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {
			double Buffer_Receive_DeltaU[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_DeltaUx[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_DeltaUy[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_DeltaUz[geometry->nVertex[val_marker]][nVar];
			unsigned long nBuffer = geometry->nVertex[val_marker]*nVar;
			int receive_from = abs(SendRecv);
			
			MPI::COMM_WORLD.Recv(&Buffer_Receive_DeltaU,nBuffer,MPI::DOUBLE,receive_from, 0);
			if (val_mesh == MESH_0) {
				MPI::COMM_WORLD.Recv(&Buffer_Receive_DeltaUx,nBuffer,MPI::DOUBLE,receive_from, 1);
				MPI::COMM_WORLD.Recv(&Buffer_Receive_DeltaUy,nBuffer,MPI::DOUBLE,receive_from, 2);
				if (nDim == 3) MPI::COMM_WORLD.Recv(&Buffer_Receive_DeltaUz,nBuffer,MPI::DOUBLE,receive_from, 3);
			}
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				Point = geometry->vertex[val_marker][iVertex]->GetNode();	
				for (iVar = 0; iVar < nVar; iVar++) {
					node[Point]->SetSolution(iVar, Buffer_Receive_DeltaU[iVertex][iVar]);
					if (val_mesh == MESH_0) {
						node[Point]->SetGradient(iVar, 0, Buffer_Receive_DeltaUx[iVertex][iVar]);
						node[Point]->SetGradient(iVar, 1, Buffer_Receive_DeltaUy[iVertex][iVar]);
						if (nDim == 3) node[Point]->SetGradient(iVar, 2, Buffer_Receive_DeltaUz[iVertex][iVar]);
					}
					if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT)
						Jacobian.DeleteValsRowi(Point*nVar+iVar);
				}
				if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT) {
					node[Point]->Set_ResConv_Zero();
					node[Point]->Set_ResVisc_Zero();
					node[Point]->SetResidualZero();
					node[Point]->SetTruncationErrorZero();
				}
			}
		}
		
		/*--- Centered scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {
			double Buffer_Receive_DeltaU[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Undivided_Laplacian[geometry->nVertex[val_marker]][nVar];
			unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
			int receive_from = abs(SendRecv);
			
			MPI::COMM_WORLD.Recv(&Buffer_Receive_DeltaU,nBuffer_Vector,MPI::DOUBLE,receive_from, 0);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Recv(&Buffer_Receive_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				Point = geometry->vertex[val_marker][iVertex]->GetNode();
				for (iVar = 0; iVar < nVar; iVar++) {
					node[Point]->SetSolution(iVar, Buffer_Receive_DeltaU[iVertex][iVar]);
					if (val_mesh == MESH_0) node[Point]->SetUndivided_Laplacian(iVar, Buffer_Receive_Undivided_Laplacian[iVertex][iVar]);
					if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT)
						Jacobian.DeleteValsRowi(Point*nVar+iVar);
				}
				if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT) {
					node[Point]->Set_ResConv_Zero();
					node[Point]->Set_ResVisc_Zero();
					node[Point]->SetResidualZero();
					node[Point]->SetTruncationErrorZero();
				}
			}
		}			
		
	}
#endif
}

CLinElectricSolution::CLinElectricSolution(void) : CSolution() { }

CLinElectricSolution::~CLinElectricSolution(void) {
	
	unsigned short iVar, iDim;
	
	delete [] Residual;
	delete [] Residual_Max;
	delete [] Solution;
	
	if (nDim == 2) {
		for (iVar = 0; iVar < 3; iVar++)
			delete [] StiffMatrix_Elem[iVar];
	}
	
	if (nDim == 3) {
		for (iVar = 0; iVar < 4; iVar++)
			delete [] StiffMatrix_Elem[iVar];
	}
	
	delete [] StiffMatrix_Elem;
	delete [] StiffMatrix_Node;
	delete [] xsol;
	delete [] rhs;
	
	// Computation of gradients by least-squares
	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;
	
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
	
}

CLinElectricSolution::CLinElectricSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	
	unsigned long nPoint;
	unsigned short nMarker;
	
	nPoint = geometry->GetnPoint();
	nDim = geometry->GetnDim();
	nMarker = config->GetnMarker_All(); 
	node = new CVariable*[nPoint];
	nVar = 1;		
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	
	Residual = new double[nVar];
	Residual_Max = new double[nVar];
	Solution = new double[nVar];
	
	
	// - - - - STRUCTURES FOR SOLVING THE LINEAR SYSTEM - - - - - - 
	// Point to point stiffness matrix
	if (nDim == 2) {
		StiffMatrix_Elem = new double* [3];
		for (unsigned short iVar = 0; iVar < 3; iVar++) {
			StiffMatrix_Elem[iVar] = new double [3];
		}
	}
	
	if (nDim == 3) {
		StiffMatrix_Elem = new double* [4];
		for (unsigned short iVar = 0; iVar < 4; iVar++) {
			StiffMatrix_Elem[iVar] = new double [4];
		}
	}	
	
	StiffMatrix_Node = new double* [1];
	for (unsigned short iVar = 0; iVar < 1; iVar++) {
		StiffMatrix_Node[iVar] = new double [1];
	}
	
	// Initialization of the structure of the whole Jacobian
	InitializeStiffMatrixStructure(geometry, config);
	xsol = new double [nPoint*nVar];
	xres = new double [nPoint*nVar];
	rhs = new double [nPoint*nVar];
	
	// - - - - STRUCTURES LEAST SQUARES GRADIENTS - - - - - - - 
	// Computation of gradients by least squares
	Smatrix = new double* [nDim]; // S matrix := inv(R)*traspose(inv(R))
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Smatrix[iDim] = new double [nDim];
	
	cvector = new double* [nVar]; // c vector := transpose(WA)*(Wb)
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		cvector[iVar] = new double [nDim];
	
	
	// - - - - INITIALIZATION - - - - - - - 
	for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
		node[iPoint] = new CPotentialVariable(0.0, nDim, nVar, config);
	
	bool restart = config->GetRestart();
	
	if (!restart) {
		for (unsigned long iPoint=0; iPoint < nPoint; iPoint++)
			node[iPoint] = new CPotentialVariable(0.0, nDim, nVar, config);
	}
	else {
		string mesh_filename = config->GetSolution_LinFileName();
		ifstream restart_file;	// definition of file
		
		char *cstr; cstr = new char [mesh_filename.size()+1];
		strcpy (cstr, mesh_filename.c_str());
		restart_file.open(cstr, ios::in);	// Open the file
		if (restart_file.fail()) {
			cout << "There is no linear restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}
		unsigned long index;
		string text_line;	// definition of the text line
		
		for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			point_line >> index >> Solution[0];
			node[iPoint] = new CPotentialVariable(Solution[0], nDim, nVar, config);
		}
		restart_file.close();
	}		
	
}

void CLinElectricSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		node[iPoint]->SetResidualZero(); // Inicialize the residual vector
	StiffMatrix.SetValZero();
}

void CLinElectricSolution::Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
											  unsigned short iMesh) {
	
	
	unsigned long iPoint, nPoint = geometry->GetnPoint();
	unsigned short iVar = 0, iter = 0;
	double norm;
	
	// Build lineal system
	for (iPoint = 0; iPoint < nPoint; iPoint++) {		
		rhs[iPoint] = node[iPoint]->GetResidual(iVar);
		xsol[iPoint] = node[iPoint]->GetSolution(iVar);
	}
	
	// Solve the system
	norm = 1; iter = 0;
	
	while (norm > config->GetCauchy_Eps()) {
		iter++;
		norm = StiffMatrix.SGSIteration(rhs,xsol);
		if (iter == config->GetnExtIter()) break;		
	}
	
	SetRes_Max(0, norm);
	
	// Update solution
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		node[iPoint]->SetSolution(0,xsol[iPoint]);
	
}

void CLinElectricSolution::Compute_Residual(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
											unsigned short iMesh) {
	
	unsigned long iPoint;
	unsigned short iVar = 0;
	
	// Build lineal system
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {		
		rhs[iPoint] = node[iPoint]->GetResidual(iVar);
		xsol[iPoint] = node[iPoint]->GetSolution(iVar);
		xres[iPoint] = 0.0;
	}
	
	StiffMatrix.MatrixVectorProduct(xsol,xres);
	
	// Update residual
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		node[iPoint]->SetResidual(0,xres[iPoint]-rhs[iPoint]);
	
	SetResidual_Smoothing(geometry, 10, 10.0);
	
}

void CLinElectricSolution::SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
												  CConfig *config, unsigned short iMesh) {
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, iPoint;
	double *Coord_i, *Coord_j, *Coord_2, a[3], b[3], Area_Local;
	unsigned short iDim, iNode;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		
		Coord_i = geometry->node[Point_0]->GetCoord();
		Coord_j = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		
		for (iDim=0; iDim < nDim; iDim++) {
			a[iDim] = Coord_i[iDim]-Coord_2[iDim];
			b[iDim] = Coord_j[iDim]-Coord_2[iDim];
		}
		
		Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0])/3.0;
		
		for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
			iPoint = geometry->elem[iElem]->GetNode(iNode);
			solver->SetCoord(geometry->node[iPoint]->GetCoord(), NULL, NULL, NULL);
			solver->SetVolume(Area_Local);
			solver->SetResidual(Residual, config);
			node[iPoint]->AddResidual(Residual);		
		}
		
	}
}

void CLinElectricSolution::Galerkin_Method(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
										   CConfig *config, unsigned short iMesh) {
	
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
		/*--- Points in edge ---*/
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		if (nDim == 3) Point_3 = geometry->elem[iElem]->GetNode(3);
		
		/*--- Points coordinates ---*/
		Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		if (nDim == 3) Coord_3 = geometry->node[Point_3]->GetCoord();
		
		if (nDim == 2) solver->SetCoord(Coord_0, Coord_1, Coord_2);
		if (nDim == 3) solver->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
		
		solver->SetResidual(StiffMatrix_Elem, config);
		
		if (nDim == 2) {
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]; StiffMatrix.AddBlock(Point_0,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]; StiffMatrix.AddBlock(Point_0,Point_1,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]; StiffMatrix.AddBlock(Point_0,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]; StiffMatrix.AddBlock(Point_1,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]; StiffMatrix.AddBlock(Point_1,Point_1,StiffMatrix_Node);		
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]; StiffMatrix.AddBlock(Point_1,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]; StiffMatrix.AddBlock(Point_2,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]; StiffMatrix.AddBlock(Point_2,Point_1,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]; StiffMatrix.AddBlock(Point_2,Point_2,StiffMatrix_Node);			
		}
		
		if (nDim == 3) {
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]; StiffMatrix.AddBlock(Point_0,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]; StiffMatrix.AddBlock(Point_0,Point_1,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]; StiffMatrix.AddBlock(Point_0,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][3]; StiffMatrix.AddBlock(Point_0,Point_3,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]; StiffMatrix.AddBlock(Point_1,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]; StiffMatrix.AddBlock(Point_1,Point_1,StiffMatrix_Node);		
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]; StiffMatrix.AddBlock(Point_1,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][3]; StiffMatrix.AddBlock(Point_1,Point_3,StiffMatrix_Node);
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]; StiffMatrix.AddBlock(Point_2,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]; StiffMatrix.AddBlock(Point_2,Point_1,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]; StiffMatrix.AddBlock(Point_2,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][3]; StiffMatrix.AddBlock(Point_2,Point_3,StiffMatrix_Node);
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][0]; StiffMatrix.AddBlock(Point_3,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][1]; StiffMatrix.AddBlock(Point_3,Point_1,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][2]; StiffMatrix.AddBlock(Point_3,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][3]; StiffMatrix.AddBlock(Point_3,Point_3,StiffMatrix_Node);
		}
	}
}

void CLinElectricSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {
	unsigned long Point, iVertex;
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		Point = geometry->vertex[val_marker][iVertex]->GetNode();
		Solution[0]= 0.0;
		node[Point]->SetSolution(Solution);
		node[Point]->SetResidual(Solution);
		StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
	}
	
}

void CLinElectricSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long Point, iVertex;
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		Point = geometry->vertex[val_marker][iVertex]->GetNode();
		Solution[0]= 0.0;
		node[Point]->SetSolution(Solution);
		node[Point]->SetResidual(Solution);
		StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
	}
}
