/*!
 * \file solution_direct_levelset.cpp
 * \brief Main subrotuines for solving the level set problem.
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

CLevelSetSolution::CLevelSetSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double dull_val;
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
	Vector = new double[nDim];
	Vector_i = new double[nDim]; Vector_j = new double[nDim];
	
	/*--- Define some auxiliar vector related with the flow solution ---*/
	FlowSolution_i = new double[nDim+3]; FlowSolution_j = new double[nDim+3];
	
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
	
	/*--- levelset_Inf at the infinity ---*/
	levelset_Inf = 0.0;
	
	/*--- Restart the solution from file information ---*/
	if (!restart || geometry->GetFinestMGLevel() == false) {
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			double levelset = geometry->node[iPoint]->GetCoord(nDim-1) - config->GetFreeSurface_Zero();
			node[iPoint] = new CLevelSetVariable(levelset, nDim, nVar, config);
		}	
	}
	else {
    
    /*--- Restart the solution from file information ---*/
		string mesh_filename = config->GetSolution_FlowFileName();
		restart_file.open(mesh_filename.data(), ios::in);
    
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
    while (getline (restart_file,text_line)) {
			istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives 
       on a different processor, the value of iPoint_Local will be -1. 
       Otherwise, the local index for this node on the current processor 
       will be returned and used to instantiate the vars. ---*/
      iPoint_Local = Global2Local[iPoint_Global];
      if (iPoint_Local >= 0) {
        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        node[iPoint_Local] = new CLevelSetVariable(Solution[0], nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
      node[iPoint] = new CLevelSetVariable(Solution[0], nDim, nVar, config);
    }
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
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

void CLevelSetSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		node[iPoint]->Set_ResConv_Zero();
		node[iPoint]->Set_ResSour_Zero();
		node[iPoint]->SetResidualZero();
	}
	
	/*--- Implicit part ---*/
	if (implicit)
		Jacobian.SetValZero();
	
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
}

void CLevelSetSolution::SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
																	unsigned short iMesh, unsigned long Iteration) {
	double *Normal, Vol, Mean_ProjVel, Lambda, Local_Delta_Time,
	Global_Delta_Time = 1E6, Global_Delta_UnstTimeND;
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iMarker;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetMax_Lambda_Inv(0.0);
	
	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		
		/*--- Point identification, Normal vector and area ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); 
		jPoint = geometry->edge[iEdge]->GetNode(1);
		Normal = geometry->edge[iEdge]->GetNormal();
		
		/*--- Mean Values ---*/
		Mean_ProjVel = 0.5 * (solution_container[FLOW_SOL]->node[iPoint]->GetProjVelInc(Normal) + 
													solution_container[FLOW_SOL]->node[jPoint]->GetProjVelInc(Normal));
		
		/*--- Inviscid contribution ---*/
		Lambda = fabs(Mean_ProjVel);
		node[iPoint]->AddMax_Lambda_Inv(Lambda);
		node[jPoint]->AddMax_Lambda_Inv(Lambda);
	}
	
	/*--- Loop boundary edges ---*/
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) { 
		for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			
			/*--- Point identification, Normal vector and area ---*/
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			
			/*--- Mean Values ---*/
			Mean_ProjVel = solution_container[FLOW_SOL]->node[iPoint]->GetProjVelInc(Normal);
			
			/*--- Inviscid contribution ---*/
			Lambda = fabs(Mean_ProjVel);
			node[iPoint]->AddMax_Lambda_Inv(Lambda);
		}
	}
	
	/*--- Each element uses their own speed, steaday state simulation ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Local_Delta_Time = config->GetLevelSet_CFLRedCoeff()*config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
		Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
		Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
		Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
		node[iPoint]->SetDelta_Time(Local_Delta_Time);
	}
		
	/*--- For exact time solution use the minimum delta time of the whole mesh ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
#ifndef NO_MPI
		double rbuf_time, sbuf_time;
		sbuf_time = Global_Delta_Time;
		MPI::COMM_WORLD.Reduce(&sbuf_time, &rbuf_time, 1, MPI::DOUBLE, MPI::MIN, MASTER_NODE);
		MPI::COMM_WORLD.Bcast(&rbuf_time, 1, MPI::DOUBLE, MASTER_NODE);
		MPI::COMM_WORLD.Barrier();
		Global_Delta_Time = rbuf_time;
#endif
		for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
			node[iPoint]->SetDelta_Time(Global_Delta_Time);
	}
	
	/*--- Recompute the unsteady time step for the dual time stratey 
	 if the unsteady CFL is diferent from 0 ---*/
	if (((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) && 
			(Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
		Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);
		
#ifndef NO_MPI
		double rbuf_time, sbuf_time;
		sbuf_time = Global_Delta_UnstTimeND;
		MPI::COMM_WORLD.Reduce(&sbuf_time, &rbuf_time, 1, MPI::DOUBLE, MPI::MIN, MASTER_NODE);
		MPI::COMM_WORLD.Bcast(&rbuf_time, 1, MPI::DOUBLE, MASTER_NODE);
		MPI::COMM_WORLD.Barrier();
		Global_Delta_UnstTimeND = rbuf_time;
#endif
		config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
	}
	
	/*--- The pseudo local time cannot be greater than the physical time ---*/
	if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
			|| (config->GetUnsteady_Simulation() == DT_STEPPING_2ND))
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
			if (!implicit) { 
				Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
				node[iPoint]->SetDelta_Time(Local_Delta_Time);
			}
	
}

void CLevelSetSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {
	double *LevelSet_var_i, *LevelSet_var_j, *U_i, *U_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, DensityInc_i, DensityInc_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	double factor, x, x_o, x_od, x_i, x_id, DampingFactor;

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
			
			/*--- Damping factor in the flow variables ---*/			
			factor = config->GetFreeSurface_Damping_Length();
			x = 0.5*(geometry->node[iPoint]->GetCoord()[0] + geometry->node[jPoint]->GetCoord()[0]);
			
			x_o = config->GetFreeSurface_Outlet();
			x_od = x_o - factor*2.0*PI_NUMBER*config->GetFroude()*config->GetFroude();
			
			x_i = config->GetFreeSurface_Inlet();
			x_id = x_i + factor*2.0*PI_NUMBER*config->GetFroude()*config->GetFroude();
			
			DampingFactor = 1.0;
			if (x >= x_od) DampingFactor = (1.0-pow((x-x_od)/(x_o-x_od), 2.0));	
			if (x <= x_id) DampingFactor = (1.0-pow((x-x_id)/(x_i-x_id), 2.0));	
			
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
				Solution_i[iVar] = LevelSet_var_i[iVar] + Project_Grad_i;
				Solution_j[iVar] = LevelSet_var_j[iVar] + Project_Grad_j;
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

void CLevelSetSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																										CConfig *config, unsigned short iMesh) {
	unsigned long iPoint;
	double Vol, x_o, x_od, x_i, x_id, x, z, levelset, DampingFactor;
	double factor = config->GetFreeSurface_Damping_Length();
	
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);

	x_o = config->GetFreeSurface_Outlet();
	x_od = x_o - factor*2.0*PI_NUMBER*config->GetFroude()*config->GetFroude();

	x_i = config->GetFreeSurface_Inlet();
	x_id = x_i + factor*2.0*PI_NUMBER*config->GetFroude()*config->GetFroude();

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 
		
		Vol = geometry->node[iPoint]->GetVolume();
		x = geometry->node[iPoint]->GetCoord()[0];
		z = geometry->node[iPoint]->GetCoord()[nDim-1]-config->GetFreeSurface_Zero();
		levelset = node[iPoint]->GetSolution(0);

		DampingFactor = 0.0;
		if (x >= x_od) 
			DampingFactor = config->GetFreeSurface_Damping_Coeff()*pow((x-x_od)/(x_o-x_od), 2.0);	
		if (x <= x_id) 
			DampingFactor = config->GetFreeSurface_Damping_Coeff()*pow((x-x_id)/(x_i-x_id), 2.0);		
			
		Residual[0] = Vol*(levelset-z)*DampingFactor;
		Jacobian_i[0][0] = Vol*DampingFactor;
		
		node[iPoint]->AddRes_Sour(Residual);
		if (implicit)
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		
	}
	
}

void CLevelSetSolution::Source_Template(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																										CConfig *config, unsigned short iMesh) {
}


void CLevelSetSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar, iDim;
	double *U_domain, *U_wall, *LevelSet_domain, *LevelSet_wall, DensityInc_domain, DensityInc_wall;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new double[solution_container[FLOW_SOL]->GetnVar()]; 
	U_wall = new double[solution_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new double[1];
	LevelSet_wall = new double[1];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
					
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
				U_domain[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
				U_wall[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			}
			DensityInc_domain = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
			DensityInc_wall = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
			
			LevelSet_domain[0] = node[iPoint]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			solver->SetNormal(Vector);
			
			solver->SetConservative(U_domain, U_wall);
			solver->SetDensityInc(DensityInc_domain, DensityInc_wall);
			solver->SetLevelSetVar(LevelSet_domain, LevelSet_wall);
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			
			node[iPoint]->AddRes_Conv(Residual);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			
		}
	}
}

void CLevelSetSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
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
		node[iPoint]->Set_ResSour_Zero();

		/*--- Includes 1 in the diagonal ---*/
		if (implicit)
			Jacobian.DeleteValsRowi(iPoint);
	}
}

void CLevelSetSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set the solution to the original value ---*/
		Solution[0] = node[iPoint]->GetSolution(0);
		
		node[iPoint]->SetSolution_Old(Solution);
		node[iPoint]->SetResidualZero();
		node[iPoint]->Set_ResConv_Zero();
		node[iPoint]->Set_ResSour_Zero();
		
		/*--- Includes 1 in the diagonal ---*/
		if (implicit)
			Jacobian.DeleteValsRowi(iPoint);
	}
}

void CLevelSetSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	double *U_domain, *U_inlet, *LevelSet_domain, *LevelSet_inlet;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new double[solution_container[FLOW_SOL]->GetnVar()]; 
	U_inlet = new double[solution_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new double[nVar];
	LevelSet_inlet = new double[nVar];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Compute closest normal neighbor ---*/
			Point_Normal = iPoint; //geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor();
			
			/*--- Interpolated solution interior to the inlet ---*/
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				U_domain[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			
			/*--- Pressure computation using the internal value ---*/
			U_inlet[0] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(0);
			
			/*--- The velocity is computed from the interior, 
			 and normal derivative for the density ---*/
			for (iDim = 0; iDim < nDim; iDim++)
				U_inlet[iDim+1] = solution_container[FLOW_SOL]->GetVelocity_Inf(iDim)*solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
			
			/*--- Note that the y velocity is recomputed due to the 
			 free surface effect on the pressure ---*/
			U_inlet[nDim] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(nDim); 
			
			LevelSet_domain[0] = node[iPoint]->GetSolution(0);
//			LevelSet_inlet[0] = 2.0*node[iPoint]->GetSolution(0) - node[geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor()]->GetSolution(0);
			LevelSet_inlet[0] = node[Point_Normal]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			solver->SetNormal(Vector);
			
			solver->SetConservative(U_domain, U_inlet);
			solver->SetLevelSetVar(LevelSet_domain, LevelSet_inlet);
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			
			node[iPoint]->AddRes_Conv(Residual);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			
		}
	}
}

void CLevelSetSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	double *U_domain, *U_outlet, *LevelSet_domain, *LevelSet_outlet, PressFreeSurface, 
	Froude, yFreeSurface, yCoord;
	double Heaviside, LevelSet, epsilon, Density_Exit;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool gravity = config->GetGravityForce();

	U_domain = new double[solution_container[FLOW_SOL]->GetnVar()]; 
	U_outlet = new double[solution_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new double[nVar];
	LevelSet_outlet = new double[nVar];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Compute closest normal neighbor ---*/
			Point_Normal = iPoint; //geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor();
			
			/*--- Interpolated solution interior to the outlet ---*/
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				U_domain[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);

			/*--- Density computation at the exit ---*/
			Heaviside = 0.0;
			yFreeSurface = config->GetFreeSurface_Zero();
			yCoord = geometry->node[iPoint]->GetCoord(nDim-1);
			LevelSet = yCoord - yFreeSurface;
			epsilon = config->GetFreeSurface_Thickness();
			if (LevelSet < -epsilon) Heaviside = 1.0;
			if (fabs(LevelSet) <= epsilon) Heaviside = 1.0 - (0.5*(1.0+(LevelSet/epsilon)+(1.0/PI_NUMBER)*sin(PI_NUMBER*LevelSet/epsilon)));
			if (LevelSet > epsilon) Heaviside = 0.0;
			Density_Exit = (config->GetRatioDensity() + (1.0 - config->GetRatioDensity())*Heaviside)*config->GetDensity_FreeStreamND();
			
			/*--- Pressure computation using density at the exit ---*/
			PressFreeSurface = solution_container[FLOW_SOL]->GetPressure_Inf();
			Froude = config->GetFroude();
			if (gravity) U_outlet[0] = PressFreeSurface + Density_Exit*((yFreeSurface-yCoord)/(Froude*Froude));
			else U_outlet[0] = PressFreeSurface;	

			/*--- Neumman condition in the interface ---*/
			if (fabs(LevelSet) <= epsilon) {
				U_outlet[0] = node[Point_Normal]->GetSolution(0);
				Density_Exit = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
			}
			
			/*--- Velocity interpolation ---*/
			for (iDim = 0; iDim < nDim; iDim++)
				U_outlet[iDim+1] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(iDim+1); 
			
			LevelSet_domain[0] = node[iPoint]->GetSolution(0);
//			LevelSet_outlet[0] = 2.0*node[iPoint]->GetSolution(0) - node[geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor()]->GetSolution(0);
			LevelSet_outlet[0] = node[Point_Normal]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			solver->SetNormal(Vector);
			
			solver->SetConservative(U_domain, U_outlet);
			solver->SetLevelSetVar(LevelSet_domain, LevelSet_outlet);
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), Density_Exit);
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
		
			node[iPoint]->AddRes_Conv(Residual);
			
			if (implicit) Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			
		}
	}
}

void CLevelSetSolution::MPI_Send_Receive(CGeometry ***geometry, CSolution ****solution_container,
                                         CConfig **config, unsigned short iMGLevel, unsigned short iZone) {
	
#ifndef NO_MPI
	unsigned short iVar, iMarker;
	double *LevelSet_Var, **LevelSet_Grad, *Conserv_Var;
	unsigned long iVertex, iPoint;
	
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config[iZone]->GetnMarker_All(); iMarker++)
		if (config[iZone]->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			
			short SendRecv = config[iZone]->GetMarker_All_SendRecv(iMarker);
			unsigned long nVertex = geometry[iZone][iMGLevel]->nVertex[iMarker];
			
			/*--- Send information  ---*/
			if (SendRecv > 0) {
				/*--- Dimensionalization ---*/
				unsigned long nBuffer_Scalar = nVertex*nVar;
				unsigned long nBuffer_Vector = nVertex*solution_container[iZone][iMGLevel][FLOW_SOL]->GetnVar();
				
				int send_to = SendRecv-1;
				
				double *Buffer_Send_LevelSet = new double [nBuffer_Scalar];
				double *Buffer_Send_LevelSetx = new double [nBuffer_Scalar];
				double *Buffer_Send_LevelSety = new double [nBuffer_Scalar];
				double *Buffer_Send_LevelSetz = new double [nBuffer_Scalar];
				double *Buffer_Send_U = new double[nBuffer_Vector];
				
				for (iVertex = 0; iVertex < geometry[iZone][iMGLevel]->nVertex[iMarker]; iVertex++) {
					iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();
					
					Conserv_Var = solution_container[iZone][iMGLevel][FLOW_SOL]->node[iPoint]->GetSolution();
					LevelSet_Var = node[iPoint]->GetSolution();
					LevelSet_Grad = node[iPoint]->GetGradient();
					
					for (iVar = 0; iVar < nVar; iVar++) {
						Buffer_Send_LevelSet[iVar*nVertex+iVertex] = LevelSet_Var[iVar];
						Buffer_Send_LevelSetx[iVar*nVertex+iVertex] = LevelSet_Grad[iVar][0];
						Buffer_Send_LevelSety[iVar*nVertex+iVertex] = LevelSet_Grad[iVar][1];
						if (nDim == 3) Buffer_Send_LevelSetz[iVar*nVertex+iVertex] = LevelSet_Grad[iVar][2];
					}
					
					for (iVar = 0; iVar < solution_container[iZone][iMGLevel][FLOW_SOL]->GetnVar(); iVar++)
						Buffer_Send_U[iVar*nVertex+iVertex] = Conserv_Var[iVar];
					
				}
				
				MPI::COMM_WORLD.Bsend(Buffer_Send_LevelSet,nBuffer_Scalar,MPI::DOUBLE,send_to, 0);
				MPI::COMM_WORLD.Bsend(Buffer_Send_LevelSetx,nBuffer_Scalar,MPI::DOUBLE,send_to, 1);
				MPI::COMM_WORLD.Bsend(Buffer_Send_LevelSety,nBuffer_Scalar,MPI::DOUBLE,send_to, 2);
				if (nDim == 3) MPI::COMM_WORLD.Bsend(Buffer_Send_LevelSetz,nBuffer_Scalar,MPI::DOUBLE,send_to, 3);
				MPI::COMM_WORLD.Bsend(Buffer_Send_U, nBuffer_Vector, MPI::DOUBLE, send_to, 4);
				
				delete[] Buffer_Send_LevelSet;
				delete[] Buffer_Send_LevelSetx;
				delete[] Buffer_Send_LevelSety;
				delete[] Buffer_Send_LevelSetz;
				delete[] Buffer_Send_U;
				
			}
			
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
				
				/*--- Dimensionalization ---*/
				unsigned long nBuffer_Scalar = nVertex*nVar;
				unsigned long nBuffer_Vector = nVertex*solution_container[iZone][iMGLevel][FLOW_SOL]->GetnVar();
				
				int receive_from = abs(SendRecv)-1;
				
				double *Buffer_Receive_LevelSet = new double [nBuffer_Scalar];
				double *Buffer_Receive_LevelSetx = new double [nBuffer_Scalar];
				double *Buffer_Receive_LevelSety = new double [nBuffer_Scalar];
				double *Buffer_Receive_LevelSetz = new double [nBuffer_Scalar];
				double *Buffer_Receive_U = new double [nBuffer_Vector];
				
				MPI::COMM_WORLD.Recv(Buffer_Receive_LevelSet,nBuffer_Scalar,MPI::DOUBLE,receive_from, 0);
				MPI::COMM_WORLD.Recv(Buffer_Receive_LevelSetx,nBuffer_Scalar,MPI::DOUBLE,receive_from, 1);
				MPI::COMM_WORLD.Recv(Buffer_Receive_LevelSety,nBuffer_Scalar,MPI::DOUBLE,receive_from, 2);
				if (nDim == 3) MPI::COMM_WORLD.Recv(Buffer_Receive_LevelSetz,nBuffer_Scalar,MPI::DOUBLE,receive_from, 3);
				MPI::COMM_WORLD.Recv(Buffer_Receive_U, nBuffer_Vector, MPI::DOUBLE, receive_from, 4);
				
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();
					for (iVar = 0; iVar < nVar; iVar++) {
						node[iPoint]->SetSolution(iVar, Buffer_Receive_LevelSet[iVar*nVertex+iVertex]);
						node[iPoint]->SetGradient(iVar, 0, Buffer_Receive_LevelSetx[iVar*nVertex+iVertex]);
						node[iPoint]->SetGradient(iVar, 1, Buffer_Receive_LevelSety[iVar*nVertex+iVertex]);
						if (nDim == 3) node[iPoint]->SetGradient(iVar, 2, Buffer_Receive_LevelSetz[iVar*nVertex+iVertex]);
					}
					
					for (iVar = 0; iVar < solution_container[iZone][iMGLevel][FLOW_SOL]->GetnVar(); iVar++)
						solution_container[iZone][iMGLevel][FLOW_SOL]->node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertex+iVertex]);
					
				}
				
				delete[] Buffer_Receive_LevelSet;
				delete[] Buffer_Receive_LevelSetx;
				delete[] Buffer_Receive_LevelSety;
				delete[] Buffer_Receive_LevelSetz;
				delete[] Buffer_Receive_U;
				
			}
		}
#endif
}

void CLevelSetSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																		 CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar, iDim;
	double *U_domain, *U_wall, *LevelSet_domain, *LevelSet_wall, DensityInc_domain, DensityInc_wall;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	U_domain = new double[solution_container[FLOW_SOL]->GetnVar()]; 
	U_wall = new double[solution_container[FLOW_SOL]->GetnVar()];
	LevelSet_domain = new double[1];
	LevelSet_wall = new double[1];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
				U_domain[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
				U_wall[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			}
			DensityInc_domain = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
			DensityInc_wall = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
			
			LevelSet_domain[0] = node[iPoint]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			solver->SetNormal(Vector);
			
			solver->SetConservative(U_domain, U_wall);
			solver->SetDensityInc(DensityInc_domain, DensityInc_wall);
			solver->SetLevelSetVar(LevelSet_domain, LevelSet_wall);
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			
			node[iPoint]->AddRes_Conv(Residual);
			
			if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			
		}
	}
}

void CLevelSetSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned long iPoint;
	double Delta = 0.0, Res, *local_ResConv, *local_ResSour, Vol;
	
	/*--- Set maximum residual to zero ---*/
	SetRes_Max(0,0.0);
	
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		/*--- Read the residual ---*/
		local_ResConv = node[iPoint]->GetResConv();
		local_ResSour = node[iPoint]->GetResSour();
		
		/*--- Read the volume ---*/
		Vol = geometry->node[iPoint]->GetVolume();
		
		/*--- Modify matrix diagonal to assure diagonal dominance ---*/
		Delta = Vol / node[iPoint]->GetDelta_Time();
		
		Jacobian.AddVal2Diag(iPoint,Delta);
			
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		Res = local_ResConv[0]+local_ResSour[0];
		
		rhs[iPoint] = -Res;
		xsol[iPoint] = 0.0;
		AddRes_Max( 0, Res*Res*Vol );
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
	
	/*--- Update solution (system written in terms of increments), be careful with the update of the 
	 scalar equations which includes the density ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->AddSolution(0, xsol[iPoint]);

#ifdef NO_MPI
	/*--- Compute the norm-2 of the residual ---*/
	SetRes_Max(0, sqrt(GetRes_Max(0)));
#endif

}

void CLevelSetSolution::ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	double *Residual, Vol, Delta, *Coord;
	unsigned long iPoint;
	
	/*--- Set maximum residual to zero ---*/
	SetRes_Max(0, 0.0);
	
	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Coord = geometry->node[iPoint]->GetCoord();

		Delta = node[iPoint]->GetDelta_Time() / Vol;
		
		Residual = node[iPoint]->GetResidual();
		node[iPoint]->AddSolution(0, - Residual[0]*Delta);

		AddRes_Max( 0, Residual[0]*Residual[0]*Vol );	
	}
	
#ifdef NO_MPI
	/*--- Compute the norm-2 of the residual ---*/
	SetRes_Max( 0, sqrt(GetRes_Max(0)) );
#endif
	
}

void CLevelSetSolution::ExplicitRK_Iteration(CGeometry *geometry, CSolution **solution_container, 
																					CConfig *config, unsigned short iRKStep) {
	double *Residual, Vol, Delta;
	unsigned long iPoint;
	double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
	
	/*--- Set maximum residual to zero ---*/
	SetRes_Max(0, 0.0);
	
	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();

		Delta = node[iPoint]->GetDelta_Time() / Vol;
		
		Residual = node[iPoint]->GetResidual();
		node[iPoint]->AddSolution(0, - Residual[0]*Delta*RK_AlphaCoeff);

		AddRes_Max( 0, Residual[0]*Residual[0]*Vol );
	}
	
#ifdef NO_MPI
	/*--- Compute the norm-2 of the residual ---*/
	SetRes_Max( 0, sqrt(GetRes_Max(0)) );
#endif
	
}

void CLevelSetSolution::SetLevelSet_Distance(CGeometry *geometry, CConfig *config) {
	
	double *coord, dist2, dist, *iCoord, *jCoord, *U_i, *U_j;
	unsigned short iDim;
	unsigned long iPoint, jPoint, iVertex, nVertex_LevelSet, iEdge;
	
#ifdef NO_MPI
	
	/*--- Identification of the 0 level set points and coordinates ---*/
	nVertex_LevelSet = 0;
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0); U_i = node[iPoint]->GetSolution();
		jPoint = geometry->edge[iEdge]->GetNode(1); U_j = node[jPoint]->GetSolution();
		if (U_i[0]*U_j[0] < 0.0) nVertex_LevelSet ++;
	}
	
	/*--- Allocate vector of boundary coordinates ---*/
	double **Coord_LevelSet;
	Coord_LevelSet = new double* [nVertex_LevelSet];
	for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
		Coord_LevelSet[iVertex] = new double [nDim];
	
	/*--- Get coordinates of the points of the surface ---*/
	nVertex_LevelSet = 0;
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0); U_i = node[iPoint]->GetSolution(); iCoord = geometry->node[iPoint]->GetCoord();
		jPoint = geometry->edge[iEdge]->GetNode(1); U_j = node[jPoint]->GetSolution(); jCoord = geometry->node[jPoint]->GetCoord();
		if (U_i[0]*U_j[0] < 0.0) {
			for (iDim = 0; iDim < nDim; iDim++)
				Coord_LevelSet[nVertex_LevelSet][iDim] = iCoord[iDim]-U_i[0]*(jCoord[iDim]-iCoord[iDim])/(U_j[0]-U_i[0]);
			nVertex_LevelSet++;
		}
	}
	
	/*--- Get coordinates of the points and compute distances to the surface ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		coord = geometry->node[iPoint]->GetCoord();
		/*--- Compute the squared distance to the rest of points, and get the minimum ---*/
		dist = 1E20;
		for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++) {
			dist2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				dist2 += (coord[iDim]-Coord_LevelSet[iVertex][iDim])*(coord[iDim]-Coord_LevelSet[iVertex][iDim]);
			if (dist2 < dist) dist = dist2;
		}
		double NumberSign = 1.0;
		if (node[iPoint]->GetSolution(0) != 0.0) NumberSign = node[iPoint]->GetSolution(0)/fabs(node[iPoint]->GetSolution(0));
		node[iPoint]->SetSolution(0, sqrt(dist)*NumberSign);
	}
	
	/*--- Deallocate vector of boundary coordinates ---*/
	for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
		delete [] Coord_LevelSet[iVertex];
	delete [] Coord_LevelSet;
	
#else 
		
	/*--- Update the Level Set solution ---*/
	for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {

			double *LevelSet_Var;
			unsigned long iVertex, iPoint;
			
			short SendRecv = config->GetMarker_All_SendRecv(iMarker);
			unsigned long nVertex = geometry->nVertex[iMarker];
			
			/*--- Send information  ---*/
			if (SendRecv > 0) {
				
				/*--- Dimensionalization ---*/
				unsigned long nBuffer_Scalar = nVertex*nVar;
				int send_to = SendRecv-1;
				double *Buffer_Send_LevelSet = new double [nBuffer_Scalar];
				
				for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					LevelSet_Var = node[iPoint]->GetSolution();
					Buffer_Send_LevelSet[iVertex] = LevelSet_Var[0];
				}
				
				MPI::COMM_WORLD.Bsend(Buffer_Send_LevelSet,nBuffer_Scalar,MPI::DOUBLE,send_to, 0);
				
				delete[] Buffer_Send_LevelSet;
			}
			
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
				
				/*--- Dimensionalization ---*/
				unsigned long nBuffer_Scalar = nVertex*nVar;
				int receive_from = abs(SendRecv)-1;
				double *Buffer_Receive_LevelSet = new double [nBuffer_Scalar];
				
				MPI::COMM_WORLD.Recv(Buffer_Receive_LevelSet,nBuffer_Scalar,MPI::DOUBLE,receive_from, 0);
				
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					node[iPoint]->SetSolution(0, Buffer_Receive_LevelSet[iVertex]);
				}
				
				delete[] Buffer_Receive_LevelSet;
			}
		}
		
	int iProcessor;
	
	/*--- Count the number of wall nodes in the whole mesh ---*/
	unsigned long nLocalVertex_LevelSet = 0, nGlobalVertex_LevelSet = 0, MaxLocalVertex_LevelSet = 0;
	
	int nProcessor = MPI::COMM_WORLD.Get_size();
	
	unsigned long *Buffer_Send_nVertex = new unsigned long [1];
	unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
	
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0); U_i = node[iPoint]->GetSolution();
		jPoint = geometry->edge[iEdge]->GetNode(1); U_j = node[jPoint]->GetSolution();
		if (U_i[0]*U_j[0] < 0.0) nLocalVertex_LevelSet ++;
	}
	
	Buffer_Send_nVertex[0] = nLocalVertex_LevelSet;	
	
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_LevelSet, &nGlobalVertex_LevelSet, 1, MPI::UNSIGNED_LONG, MPI::SUM); 	
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_LevelSet, &MaxLocalVertex_LevelSet, 1, MPI::UNSIGNED_LONG, MPI::MAX); 	
	MPI::COMM_WORLD.Allgather(Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG);
	
	double *Buffer_Send_Coord = new double [MaxLocalVertex_LevelSet*nDim];
	double *Buffer_Receive_Coord = new double [nProcessor*MaxLocalVertex_LevelSet*nDim];
	unsigned long nBuffer = MaxLocalVertex_LevelSet*nDim;
	
	for (iVertex = 0; iVertex < MaxLocalVertex_LevelSet; iVertex++)
		for (iDim = 0; iDim < nDim; iDim++)
			Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
	
	nVertex_LevelSet = 0;
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0); U_i = node[iPoint]->GetSolution(); iCoord = geometry->node[iPoint]->GetCoord();
		jPoint = geometry->edge[iEdge]->GetNode(1); U_j = node[jPoint]->GetSolution(); jCoord = geometry->node[jPoint]->GetCoord();

		if (U_i[0]*U_j[0] < 0.0) {
			for (iDim = 0; iDim < nDim; iDim++)
				Buffer_Send_Coord[nVertex_LevelSet*nDim+iDim] = iCoord[iDim]-U_i[0]*(jCoord[iDim]-iCoord[iDim])/(U_j[0]-U_i[0]);
			nVertex_LevelSet++;
		}
	}
	
	MPI::COMM_WORLD.Allgather(Buffer_Send_Coord, nBuffer, MPI::DOUBLE, Buffer_Receive_Coord, nBuffer, MPI::DOUBLE);
	
	/*--- Get coordinates of the points and compute distances to the surface ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		coord = geometry->node[iPoint]->GetCoord();
		
		/*--- Compute the squared distance and get the minimum ---*/
		dist = 1E20;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
			for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
				dist2 = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					dist2 += (coord[iDim]-Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_LevelSet+iVertex)*nDim+iDim])*
					(coord[iDim]-Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_LevelSet+iVertex)*nDim+iDim]);
				if (dist2 < dist) dist = dist2;
			}
		double NumberSign = 1.0;
		if (node[iPoint]->GetSolution(0) != 0.0) NumberSign = node[iPoint]->GetSolution(0)/fabs(node[iPoint]->GetSolution(0));
		node[iPoint]->SetSolution(0, sqrt(dist)*NumberSign);
	}
	
	delete [] Buffer_Send_Coord;
	delete [] Buffer_Receive_Coord;
	delete [] Buffer_Send_nVertex;
	delete [] Buffer_Receive_nVertex;
	
#endif
	
}

void CLevelSetSolution::SetResidual_DualTime(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep, 
																		 unsigned short iMesh, unsigned short RunTime_EqSystem) {

	unsigned short iVar, jVar;
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1, Volume_nM1, Volume_n, Volume_nP1, TimeStep;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool Grid_Movement = config->GetGrid_Movement();

	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 
		
		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();
		
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
		for(iVar = 0; iVar < nVar; iVar++) {
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				Residual[iVar] = ( U_time_nP1[iVar]*Volume_nP1 - U_time_n[iVar]*Volume_n ) / TimeStep;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				Residual[iVar] = ( 3.0*U_time_nP1[iVar]*Volume_nP1 - 4.0*U_time_n[iVar]*Volume_n
													+  1.0*U_time_nM1[iVar]*Volume_nM1 ) / (2.0*TimeStep);
		}
				
		/*--- Add Residual ---*/
		node[iPoint]->AddRes_Conv(Residual);
		
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

