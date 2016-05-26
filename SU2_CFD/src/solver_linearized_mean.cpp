/*!
 * \file solution_linearized_mean.cpp
 * \brief Main subrotuines for solving linearized problems (Euler, Navier-Stokes, etc.).
 * \author F. Palacios
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

CLinEulerSolver::CLinEulerSolver(void) : CSolver() { }

CLinEulerSolver::CLinEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
	unsigned long iPoint, index;
	string text_line, mesh_filename;
	unsigned short iDim, iVar, nLineLets;
	ifstream restart_file;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
	bool restart = config->GetRestart();
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	su2double dull_val;
  
	/*--- Define geometry constants in the solver structure ---*/
	nDim = geometry->GetnDim();
	nVar = geometry->GetnDim()+2;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
	node = new CVariable*[geometry->GetnPoint()];
	
	/*--- Define some auxiliar vector related with the residual ---*/
	Residual = new su2double[nVar];	Residual_RMS = new su2double[nVar];  
	Residual_i = new su2double[nVar]; Residual_j = new su2double[nVar];
	Res_Conv = new su2double[nVar];	Res_Visc = new su2double[nVar]; Res_Sour = new su2double[nVar];
  Residual_Max = new su2double[nVar];
  
  /*--- Define some structures for locating max residuals ---*/
  Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

	/*--- Define some auxiliar vector related with the solution ---*/
	Solution   = new su2double[nVar];
	Solution_i = new su2double[nVar];	Solution_j = new su2double[nVar];
	
	/*--- Define some auxiliar vector related with the geometry ---*/
	Vector_i = new su2double[nDim];	Vector_j = new su2double[nDim];
	
	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT) {
		/*--- Point to point Jacobians ---*/
		Jacobian_i = new su2double* [nVar]; Jacobian_j = new su2double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new su2double [nVar]; Jacobian_j[iVar] = new su2double [nVar]; }
		/*--- Initialization of the structure of the whole Jacobian ---*/
		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
	}
	
	/*--- Computation of gradients by least squares ---*/
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
	
	/*--- Delta forces definition and coefficient in all the markers ---*/
	DeltaForceInviscid = new su2double[nDim];
	CDeltaDrag_Inv = new su2double[config->GetnMarker_All()];
	CDeltaLift_Inv = new su2double[config->GetnMarker_All()];
	
	/*--- Linearized flow at the inifinity, inizialization stuff ---*/
	DeltaRho_Inf = 0.0;	
	DeltaE_Inf = 0.0;
	DeltaVel_Inf = new su2double [nDim];
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
	if (!restart || (iMesh != MESH_0)) {
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
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
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
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
        node[iPoint_Local] = new CLinEulerVariable(Solution, nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
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
	bool second_order = ((config->GetKind_Centered() == JST) && (iMesh == MESH_0));
	
	if (second_order) 
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
		if (second_order) 
			numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(), node[jPoint]->GetUndivided_Laplacian());
		
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
	su2double *Residual, Vol, Delta, *Res_TruncError;
	unsigned short iVar;
	unsigned long iPoint;
	su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);

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
        AddRes_Max(iVar, fabs(Residual[iVar]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
			}
		}
	
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);

}

void CLinEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
	unsigned long iPoint;
	
	/*--- Residual inicialization ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		LinSysRes.SetBlock_Zero(iPoint);
	}
	
	/*--- Inicialize the Jacobian matrices ---*/
	if (config->GetKind_TimeIntScheme_LinFlow() == EULER_IMPLICIT)
		Jacobian.SetValZero();
}

void CLinEulerSolver::Inviscid_DeltaForces(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	unsigned long iVertex, Point;
	unsigned short iDim, iMarker, Boundary, Monitoring;
	su2double  *Face_Normal, DeltaPressure, *Velocity;
	su2double Alpha = config->GetAoA()*PI_NUMBER / 180.0;
	su2double Beta  = config->GetAoS()*PI_NUMBER / 180.0;
	su2double RefAreaCoeff = config->GetRefAreaCoeff();
	su2double Density_Inf = solver_container[FLOW_SOL]->GetDensity_Inf();
	su2double ModVelocity_Inf = solver_container[FLOW_SOL]->GetModVelocity_Inf();
	su2double C_p = 1.0/(0.5*Density_Inf*RefAreaCoeff*ModVelocity_Inf*ModVelocity_Inf);

	/*-- Inicialization ---*/
	Total_CDeltaDrag = 0.0; Total_CDeltaLift = 0.0;
	AllBound_CDeltaDrag_Inv = 0.0; AllBound_CDeltaLift_Inv = 0.0;
	Velocity = new su2double[nDim];

	/*--- Loop over the Euler markers ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_KindBC(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);
		if ((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL)) {
			for (iDim = 0; iDim < nDim; iDim++) DeltaForceInviscid[iDim] = 0.0;
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				Point = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[Point]->GetDomain()) {
					
					/*--- Compute pressure on the boundary ---*/
					for (iDim = 0; iDim < nDim; iDim++) 
						Velocity[iDim] = solver_container[FLOW_SOL]->node[Point]->GetVelocity(iDim);
					
					su2double rho = solver_container[FLOW_SOL]->node[Point]->GetSolution(0) + node[Point]->GetSolution(0);
					su2double rhoE = solver_container[FLOW_SOL]->node[Point]->GetSolution(nVar-1) + node[Point]->GetSolution(nVar-1);
					su2double Pressure = solver_container[FLOW_SOL]->node[Point]->GetPressure();
					su2double rhoVel[3];
					su2double sqr_vel = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						rhoVel[iDim] = solver_container[FLOW_SOL]->node[Point]->GetSolution(iDim+1) + node[Point]->GetSolution(iDim+1);
						sqr_vel += rhoVel[iDim]*rhoVel[iDim]/(rho*rho);
					}
					DeltaPressure = Gamma_Minus_One*rho*(rhoE/rho-0.5*sqr_vel)-Pressure;
					
					if (Monitoring == YES) {
						Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
						for (iDim = 0; iDim < nDim; iDim++)
							DeltaForceInviscid[iDim] -= C_p*DeltaPressure*Face_Normal[iDim];
					}
				}
			}
			
			/*--- Transform ForceInviscid into CLift and CDrag ---*/
			if (Monitoring == YES) {
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
	su2double phi, a1, a2, sqvel, d, *Face_Normal, *U, dS, VelxDeltaRhoVel, Energy, Rho, 
		   *Normal, *Delta_RhoVel, *Vel, *DeltaU, Delta_RhoE, Delta_Rho;

	Normal = new su2double[nDim];
	Delta_RhoVel = new su2double[nDim];
	Vel = new su2double[nDim];
	
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
	
	su2double *Face_Normal;
	su2double *kappa = new su2double[nDim];
	su2double *velocity = new su2double[nDim];
	su2double *U_wall = new su2double[nVar];
	su2double *U_infty = new su2double[nVar];
	su2double *DeltaU_wall = new su2double[nVar];
	su2double *DeltaU_infty = new su2double[nVar];
	su2double *DeltaU_update = new su2double[nVar];
	su2double *W_wall = new su2double[nVar];
	su2double *W_infty = new su2double[nVar];
	su2double *W_update = new su2double[nVar];
	su2double Mach = config->GetMach();
	su2double DeltaMach = 1.0;

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
			su2double dS = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				dS += Face_Normal[iDim]*Face_Normal[iDim];
			dS = sqrt (dS);
			
			for (iDim = 0; iDim < nDim; iDim++)
				kappa[iDim] = -Face_Normal[iDim]/dS;
			
			/*--- Computation of P and inverse P matrix using values at the infinity ---*/
			su2double sq_vel = 0.0, vn = 0.0;
			su2double rho = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(0);
			su2double rhoE = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(nVar-1);
			for (iDim = 0; iDim < nDim; iDim++) {
				velocity[iDim] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iDim+1)/rho;
				sq_vel +=velocity[iDim]*velocity[iDim];
				vn += velocity[iDim]*kappa[iDim]*dS;
			}
			
/*			su2double sq_vel = 0.0, vn = 0.0;
			su2double rho = U_infty[0];
			su2double rhoE = U_infty[nVar-1];
			for (iDim = 0; iDim < nDim; iDim++) {
				velocity[iDim] = U_infty[iDim+1]/rho;
				sq_vel +=velocity[iDim]*velocity[iDim];
				vn += velocity[iDim]*kappa[iDim]*dS;
			}*/
			
			su2double c = sqrt(Gamma*Gamma_Minus_One*(rhoE/rho-0.5*sq_vel));
			su2double energy = rhoE / rho;
			
			su2double **P_Matrix, **invP_Matrix;
			P_Matrix = new su2double* [nVar];
			invP_Matrix = new su2double* [nVar];
			for (iVar = 0; iVar < nVar; iVar++) {
				P_Matrix[iVar] = new su2double [nVar];
				invP_Matrix[iVar] = new su2double [nVar];
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
				if (vn > 0.0) { 
					W_update[0] = W_wall[0];
					W_update[1] = W_wall[1];
				}
				else {
					W_update[0] = W_infty[0];
					W_update[1] = W_infty[1];
				}
				
				if (vn+c*dS > 0.0) W_update[2] = W_wall[2];
				else W_update[2] = W_infty[2];
				
				if (vn-c*dS > 0.0) W_update[3] = W_wall[3];
				else W_update[3] = W_infty[3];
			}
			
			if (nDim == 3) {
				if (vn > 0.0) { 
					W_update[0] = W_wall[0];
					W_update[1] = W_wall[1];
					W_update[2] = W_wall[2];
				}
				else {
					W_update[0] = W_infty[0];
					W_update[1] = W_infty[1];
					W_update[2] = W_infty[2];
				}
				
				if (vn+c*dS > 0.0) W_update[3] = W_wall[3];
				else W_update[3] = W_infty[3];
				
				if (vn-c*dS > 0.0) W_update[4] = W_wall[4];
				else W_update[4] = W_infty[4];
			}
			
			/*--- conservative variables using characteristics ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				DeltaU_update[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					DeltaU_update[iVar] +=P_Matrix[iVar][jVar]*W_update[jVar];
			}
			
			/*--- Residual computation ---*/
			su2double **Jac_Matrix;
			Jac_Matrix = new su2double* [nVar];
			for (iVar = 0; iVar < nVar; iVar++) {
				Jac_Matrix[iVar] = new su2double [nVar];
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
