/*!
 * \file solution_linearized_mean.cpp
 * \brief Main subrotuines for solving linearized problems (Euler, Navier-Stokes, etc.).
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

CLinEulerSolver::CLinEulerSolver(void) : CSolver() { }

CLinEulerSolver::CLinEulerSolver(CGeometry *geometry, CConfig *config) : CSolver() {
	unsigned long iPoint, index;
	string text_line, mesh_filename;
	unsigned short iDim, iVar;
	ifstream restart_file;
	bool restart = config->GetRestart();
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	
	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	nVar = geometry->GetnDim()+2;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
	node = new CVariable*[geometry->GetnPoint()];
	
	/*--- Define some auxiliar vector related with the residual ---*/
	Residual = new double[nVar];	Residual_RMS = new double[nVar];  
	Residual_i = new double[nVar]; Residual_j = new double[nVar];
	Res_Conv = new double[nVar];	Res_Visc = new double[nVar]; Res_Sour = new double[nVar];
  Residual_Max = new double[nVar]; Point_Max = new unsigned long[nVar];

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
	if (!restart || geometry->GetFinestMGLevel() == false) {
		for (iPoint=0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CLinEulerVariable(DeltaRho_Inf, DeltaVel_Inf, DeltaE_Inf, nDim, nVar, config);
	}
	else {
    
    /*--- Restart the solution from file information ---*/
		mesh_filename = config->GetSolution_FlowFileName();
    restart_file.open(mesh_filename.data(), ios::in);
		
    /*--- In case there is no file ---*/
		if (restart_file.fail()) {
			cout << "There is no linearized restart file!!" << endl;
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
    long iPoint_Local; unsigned long iPoint_Global = 0;
    
    /*--- The first line is the header ---*/
    getline (restart_file, text_line);
    
    while (getline (restart_file,text_line)) {
			istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives 
       on a different processor, the value of iPoint_Local will be -1. 
       Otherwise, the local index for this node on the current processor 
       will be returned and used to instantiate the vars. ---*/
      iPoint_Local = Global2Local[iPoint_Global];
      if (iPoint_Local >= 0) {
          if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
          if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
        node[iPoint_Local] = new CLinEulerVariable(Solution, nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
      node[iPoint] = new CLinEulerVariable(Solution, nDim, nVar, config);
    }
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
	}
}

CLinEulerSolver::~CLinEulerSolver(void) {
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

void CLinEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
											 CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint;
	bool implicit = (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT);
	bool high_order_diss = ((config->GetKind_Centered() == JST) && (iMesh == MESH_0));
	
	if (high_order_diss) 
		SetUndivided_Laplacian(geometry, config);
	
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		/*--- Points in edge, normal, and neighbors---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
		numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());

		/*--- Linearized variables w/o reconstruction ---*/
		numerics->SetLinearizedVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());

		/*--- Conservative variables w/o reconstruction ---*/
		numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(), 
								solver_container[FLOW_SOL]->node[jPoint]->GetSolution());

		/*--- SoundSpeed enthalpy and lambda variables w/o reconstruction ---*/
		numerics->SetSoundSpeed(solver_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
							  solver_container[FLOW_SOL]->node[jPoint]->GetSoundSpeed());
		numerics->SetEnthalpy(solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
							solver_container[FLOW_SOL]->node[jPoint]->GetEnthalpy());
		numerics->SetLambda(solver_container[FLOW_SOL]->node[iPoint]->GetLambda(), 
						  solver_container[FLOW_SOL]->node[jPoint]->GetLambda());

		/*--- Undivided laplacian ---*/
		if (high_order_diss) 
			numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(),node[jPoint]->GetUndivided_Laplacian());
		
		/*--- Compute residual ---*/
		numerics->ComputeResidual(Res_Conv, Res_Visc, Jacobian_i, Jacobian_j, config);
		
		/*--- Update convective and artificial dissipation residuals ---*/
		LinSysRes.AddBlock(iPoint, Res_Conv);
		LinSysRes.SubtractBlock(jPoint, Res_Conv);
    LinSysRes.AddBlock(iPoint, Res_Visc);
    LinSysRes.SubtractBlock(jPoint, Res_Visc);

		/*--- Set implicit stuff ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
			Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
			Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
		}
	}
}

void CLinEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) { }

void CLinEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, 
											CConfig *config, unsigned short iRKStep) {
	double *Residual, Vol, Delta, *Res_TruncError;
	unsigned short iVar;
	unsigned long iPoint;
	double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);

	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
	
	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			Vol = geometry->node[iPoint]->GetVolume();
			Delta = solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() / Vol;
			Res_TruncError = node[iPoint]->GetResTruncError();
			Residual = LinSysRes.GetBlock(iPoint);
			for (iVar = 0; iVar < nVar; iVar++) {
				node[iPoint]->AddSolution(iVar, -(Residual[iVar]+Res_TruncError[iVar])*Delta*RK_AlphaCoeff);
				AddRes_RMS(iVar, Residual[iVar]*Residual[iVar]);
        AddRes_Max(iVar, fabs(Residual[iVar]), geometry->node[iPoint]->GetGlobalIndex());
			}
		}
	
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);

}

void CLinEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
	
	/*--- Residual inicialization ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		LinSysRes.SetBlock_Zero(iPoint);
	}
	
	/*--- Inicialize the jacobian matrices ---*/
	if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT)
		Jacobian.SetValZero();
}

void CLinEulerSolver::Inviscid_DeltaForces(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	unsigned long iVertex, Point;
	unsigned short iDim, iMarker, Boundary, Monitoring;
	double  *Face_Normal, dS, DeltaPressure, *Velocity;
	double Alpha = config->GetAoA()*PI_NUMBER / 180.0;
	double Beta  = config->GetAoS()*PI_NUMBER / 180.0;
	double RefAreaCoeff = config->GetRefAreaCoeff();
	double Density_Inf = solver_container[FLOW_SOL]->GetDensity_Inf();
	double ModVelocity_Inf = solver_container[FLOW_SOL]->GetModVelocity_Inf();
	double C_p = 1.0/(0.5*Density_Inf*RefAreaCoeff*ModVelocity_Inf*ModVelocity_Inf);
	bool incompressible = config->GetIncompressible();

	/*-- Inicialization ---*/
	Total_CDeltaDrag = 0.0; Total_CDeltaLift = 0.0;
	AllBound_CDeltaDrag_Inv = 0.0; AllBound_CDeltaLift_Inv = 0.0;
	Velocity = new double[nDim];

	/*--- Loop over the Euler markers ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);
		if ((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL)) {
			for (iDim = 0; iDim < nDim; iDim++) DeltaForceInviscid[iDim] = 0.0;
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				Point = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[Point]->GetDomain()) {
					
					/*--- Compute pressure on the boundary ---*/
					for (iDim = 0; iDim < nDim; iDim++) 
						Velocity[iDim] = solver_container[FLOW_SOL]->node[Point]->GetVelocity(iDim, config->GetIncompressible());
					
					double rho = solver_container[FLOW_SOL]->node[Point]->GetSolution(0) + node[Point]->GetSolution(0);
					double rhoE = solver_container[FLOW_SOL]->node[Point]->GetSolution(nVar-1) + node[Point]->GetSolution(nVar-1);
					double Pressure = solver_container[FLOW_SOL]->node[Point]->GetPressure(incompressible);
					double rhoVel[3];
					double sqr_vel = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						rhoVel[iDim] = solver_container[FLOW_SOL]->node[Point]->GetSolution(iDim+1) + node[Point]->GetSolution(iDim+1);
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


void CLinEulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iDim, iVar;
	unsigned long iVertex, iPoint;
	double phi, a1, a2, sqvel, d, *Face_Normal, *U, dS, VelxDeltaRhoVel, Energy, Rho, 
		   *Normal, *Delta_RhoVel, *Vel, *DeltaU, Delta_RhoE, Delta_Rho;

	Normal = new double[nDim];
	Delta_RhoVel = new double[nDim];
	Vel = new double[nDim];
	
	/*--- Loop over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		/*--- Point, normal vector, and area ---*/
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			
			Face_Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			dS = 0; for (iDim = 0; iDim < nDim; iDim++) dS += Face_Normal[iDim]*Face_Normal[iDim]; dS = sqrt(dS);
			
			/*--- Linearized and flow solution ---*/
			DeltaU = node[iPoint]->GetSolution();
			U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
			
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
			LinSysRes.AddBlock(iPoint, Residual);
		}
	}
	
	delete [] Vel;
	delete [] Normal;
	delete [] Delta_RhoVel;
}

void CLinEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
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

	/*--- Loop over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_wall[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			
			/*--- Flow Solution at the infinity ---*/
			U_infty[0] = solver_container[FLOW_SOL]->GetDensity_Inf();
			U_infty[1] = solver_container[FLOW_SOL]->GetDensity_Velocity_Inf(0);
			U_infty[2] = solver_container[FLOW_SOL]->GetDensity_Velocity_Inf(1);
			U_infty[3] = solver_container[FLOW_SOL]->GetDensity_Energy_Inf();
			if (nDim == 3) {
				U_infty[3] = solver_container[FLOW_SOL]->GetDensity_Velocity_Inf(2);
				U_infty[4] = solver_container[FLOW_SOL]->GetDensity_Energy_Inf();
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
			double rho = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(0);
			double rhoE = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(nVar-1);
			for (iDim = 0; iDim < nDim; iDim++) {
				velocity[iDim] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iDim+1)/rho;
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
			
			conv_numerics->GetPMatrix_inv(&rho, velocity, &c, kappa, invP_Matrix);
			conv_numerics->GetPMatrix(&rho, velocity, &c, kappa, P_Matrix);
			
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
			conv_numerics->GetInviscidProjJac(velocity, &energy, kappa, 1.0, Jac_Matrix);
			for (iVar = 0; iVar < nVar; iVar++) {
				Residual[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					Residual[iVar] +=Jac_Matrix[iVar][jVar]*DeltaU_update[jVar]*dS;
			}
			
			LinSysRes.AddBlock(iPoint, Residual);
			
		}
	}
}
