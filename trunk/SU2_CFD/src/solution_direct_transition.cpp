/*!
 * \file solution_direct_turbulent.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.
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

CTransLMSolution::CTransLMSolution(void) : CTurbSolution() {}

CTransLMSolution::CTransLMSolution(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTurbSolution() {
	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double Density_Inf, Viscosity_Inf, Intermittency_Inf, tu_Inf, REth_Inf, Factor_nu_Inf, dull_val;
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
	
	if (iMesh == MESH_0) {
		
		/*--- Define some auxiliar vector related with the residual ---*/
		Residual = new double[nVar]; Residual_Max = new double[nVar];
		Residual_i = new double[nVar]; Residual_j = new double[nVar];
		
		/*--- Define some auxiliar vector related with the solution ---*/
		Solution   = new double[nVar];
		Solution_i = new double[nVar]; Solution_j = new double[nVar];
		
		/*--- Define some auxiliar vector related with the geometry ---*/
		Vector_i = new double[nDim]; Vector_j = new double[nDim];
		
		/*--- Define some auxiliar vector related with the flow solution ---*/
		FlowSolution_i = new double [nDim+2]; FlowSolution_j = new double [nDim+2];
		
		/*--- Jacobians and vector structures for implicit computations ---*/
		if (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT) {
			cout << "No implicit SAGT yet!!\n";   // TODO, Aniket
			int j;
			cin >> j;
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
		
	}
	
	/*--- Read farfield conditions from config ---*/
	Density_Inf       = config->GetDensity_FreeStreamND();
  Viscosity_Inf     = config->GetViscosity_FreeStreamND();
  Intermittency_Inf = config->GetIntermittency_FreeStream();
  tu_Inf            = config->GetTurbulenceIntensity_FreeStream();
	
  /*-- Initialize REth from correlation --*/
  if (tu_Inf <= 1.3) {
    REth_Inf = (1173.51-589.428*tu_Inf+0.2196/(tu_Inf*tu_Inf));
  } else {
    REth_Inf = 331.5*pow(tu_Inf-0.5658,-0.671);
  }
	
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
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      node[iPoint] = new CTransLMVariable(nu_tilde_Inf, Intermittency_Inf, REth_Inf, muT_Inf, nDim, nVar, config);
	}
	else {
    cout << "No SAGT restart yet!!" << endl; // TODO, Aniket
    int j;
    cin >> j;
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

CTransLMSolution::~CTransLMSolution(void){
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


void CTransLMSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iRKStep) {
	// TODO, Aniket
}

void CTransLMSolution::Postprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {
	// TODO, Aniket
}

void CTransLMSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {
	double *trans_var_i, *trans_var_j, *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	bool rotating_frame = config->GetRotating_Frame();
  bool grid_movement = config->GetGrid_Movement();
	
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
		trans_var_i = node[iPoint]->GetSolution();
		trans_var_j = node[jPoint]->GetSolution();
		solver->SetTransVar(trans_var_i,trans_var_j);
		
		/*--- Rotating Frame ---*/
		if (rotating_frame)
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());
		
    /*--- Grid Movement ---*/
		if (grid_movement)
			solver->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    		
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


void CTransLMSolution::Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																				CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint;
	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	
	if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit) {
		
    /*---  Need gradient of flow conservative variables ---*/
	  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solution_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry);
	  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solution_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
		
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
			solver->SetTransVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
			solver->SetTransVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
			
		  solver->SetConsVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient(),
																 solution_container[FLOW_SOL]->node[jPoint]->GetGradient());
			
			
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

void CTransLMSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																			 CConfig *config, unsigned short iMesh) {
  unsigned long iPoint;
	
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		/*--- Conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);
		
		/*--- Gradient of the primitive and conservative variables ---*/
		solver->SetPrimVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
		
		/*--- Laminar viscosity ---*/
		solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
		
		/*--- Turbulent variables w/o reconstruction, and its gradient ---*/
		solver->SetTransVar(node[iPoint]->GetSolution(), NULL);
		solver->SetTransVarGradient(node[iPoint]->GetGradient(), NULL);
		
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

void CTransLMSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
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

void CTransLMSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar, iDim;
	double *Normal;
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement = config->GetGrid_Movement();
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
		
		/*--- Rotating Frame ---*/
		if (rotating_frame)
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
		
    /*--- Grid Movement ---*/
		if (grid_movement)
			solver->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
    
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

void CTransLMSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
																unsigned short val_marker) { }

void CTransLMSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																 CConfig *config, unsigned short val_marker) { }
