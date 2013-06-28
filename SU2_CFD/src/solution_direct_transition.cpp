/*!
 * \file solution_direct_turbulent.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.3
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
	double Density_Inf, Viscosity_Inf, tu_Inf, nu_tilde_Inf, Factor_nu_Inf, dull_val,rey;
	ifstream restart_file;
	char *cstr;
	string text_line;
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	
    cout << "Entered constructor for CTransLMSolution -AA\n";
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	
	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	node = new CVariable*[geometry->GetnPoint()];
	
	/*--- Dimension of the problem --> dependent of the turbulent model ---*/
	nVar = 2;
	
	if (iMesh == MESH_0) {
		
		/*--- Define some auxiliar vector related with the residual ---*/
		Residual = new double[nVar]; Residual_RMS = new double[nVar];
		Residual_i = new double[nVar]; Residual_j = new double[nVar];
    Residual_Max = new double[nVar]; Point_Max = new unsigned long[nVar];
    

		/*--- Define some auxiliar vector related with the solution ---*/
		Solution   = new double[nVar];
		Solution_i = new double[nVar]; Solution_j = new double[nVar];
		
		/*--- Define some auxiliar vector related with the geometry ---*/
		Vector_i = new double[nDim]; Vector_j = new double[nDim];
		
		/*--- Define some auxiliar vector related with the flow solution ---*/
		FlowSolution_i = new double [nDim+2]; FlowSolution_j = new double [nDim+2];
		
		/*--- Jacobians and vector structures for implicit computations ---*/
		if (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT) {
//      cout << "Jacobian stuff -AA\n";
			/*--- Point to point Jacobians ---*/
			Jacobian_i = new double* [nVar];
			Jacobian_j = new double* [nVar];
			for (iVar = 0; iVar < nVar; iVar++) {
				Jacobian_i[iVar] = new double [nVar];
				Jacobian_j[iVar] = new double [nVar];
			}
			/*--- Initialization of the structure of the whole Jacobian ---*/
			Initialize_SparseMatrix_Structure(&Jacobian, nVar, nVar, geometry, config);

			xsol = new double [geometry->GetnPoint()*nVar];
			rhs = new double [geometry->GetnPoint()*nVar];
			}
		
		/*--- Computation of gradients by least squares ---*/
//    cout << "Gradients -AA\n";
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
  rey      = config->GetReynolds();
	
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
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      if (iPoint == 9745||iPoint == 9746||iPoint == 9608||iPoint == 9609) {
        node[iPoint] = new CTransLMVariable(nu_tilde_Inf, 0.0, 1100.0, nDim, nVar, config);
//      } else if (iPoint == 12328||iPoint == 12191||iPoint == 12190||iPoint == 12327) {
//        node[iPoint] = new CTransLMVariable(nu_tilde_Inf, 0.0, 100.0, nDim, nVar, config);
      } else {
        node[iPoint] = new CTransLMVariable(nu_tilde_Inf, Intermittency_Inf, REth_Inf, nDim, nVar, config);
      }
    }

//    cout << "Finished creating CTransLMVariable's -AA\n";
	}
	else {
    cout << "No LM restart yet!!" << endl; // TODO, Aniket
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

//  cout << "Completed constructor for CTransLMSolution -AA\n";
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

void CTransLMSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iPoint;

//  cout << "Preprocessing LM -AA\n";
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		node[iPoint]->SetResidualZero();
  Jacobian.SetValZero();

	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
}

void CTransLMSolution::Postprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {

  // Replace intermittency with gamma_eff
	//for (unsigned int iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
    //node[iPoint]->SetGammaEff();
}

void CTransLMSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta, Delta_flow, *local_Residual, Vol;


  /*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

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
				AddRes_RMS(iVar, local_Residual[iVar]*local_Residual[iVar]*Vol);
        AddRes_Max(iVar, fabs(local_Residual[iVar]), geometry->node[iPoint]->GetGlobalIndex());
			}
		}
  
  /*--- Initialize residual and solution at the ghost points ---*/
  for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      rhs[total_index] = 0.0;
      xsol[total_index] = 0.0;
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
		
		CSysVector rhs_vec((const unsigned int)geometry->GetnPoint(),
                       (const unsigned int)geometry->GetnPointDomain(), nVar, rhs);
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
			system.FlexibleGMRES(rhs_vec, sol_vec, *mat_vec, *precond, *sol_mpi, config->GetLinear_Solver_Error(), 
													 config->GetLinear_Solver_Iter(), false);		
		
		sol_vec.CopyToArray(xsol);
		delete mat_vec; 
		delete precond;  
    delete sol_mpi;

	}
	
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);

//			node[iPoint]->AddSolution(0,xsol[iPoint]);
//			node[iPoint]->AddSolution(1,xsol[iPoint]);
  }

  /*--- MPI solution ---*/
  SetSolution_MPI(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CTransLMSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {
//  double *trans_var_i, *trans_var_j, *U_i, *U_j;
//  unsigned long iEdge, iPoint, jPoint;
//  bool rotating_frame = config->GetRotating_Frame();
//  bool grid_movement = config->GetGrid_Movement();

//  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
//
//    /*--- Points in edge and normal vectors ---*/
//    iPoint = geometry->edge[iEdge]->GetNode(0);
//    jPoint = geometry->edge[iEdge]->GetNode(1);
//    solver->SetNormal(geometry->edge[iEdge]->GetNormal());
//
//    /*--- Conservative variables w/o reconstruction ---*/
//    U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
//    U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
//    solver->SetConservative(U_i, U_j);
//
//    /*--- Transition variables w/o reconstruction ---*/
//    trans_var_i = node[iPoint]->GetSolution();
//    trans_var_j = node[jPoint]->GetSolution();
//    solver->SetTransVar(trans_var_i,trans_var_j);
//
//    /*--- Add and subtract Residual ---*/
//    solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
//    node[iPoint]->AddResidual(Residual);
//    node[jPoint]->SubtractResidual(Residual);
//
//    /*--- Implicit part ---*/
//    Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
//    Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
//    Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
//    Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);
//
//  }

}


void CTransLMSolution::Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																				CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
//	unsigned long iEdge, iPoint, jPoint;
//	bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);

//	if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit) {
//	
//  /*---  Need gradient of flow conservative variables ---*/
//	  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solution_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry, config);
//	  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solution_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
//		
//		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
//			
//			/*--- Points in edge ---*/
//			iPoint = geometry->edge[iEdge]->GetNode(0);
//			jPoint = geometry->edge[iEdge]->GetNode(1);
//			
//			/*--- Points coordinates, and normal vector ---*/
//			solver->SetCoord(geometry->node[iPoint]->GetCoord(),
//											 geometry->node[jPoint]->GetCoord());
//			
//			solver->SetNormal(geometry->edge[iEdge]->GetNormal());
//			
//			/*--- Conservative variables w/o reconstruction ---*/
//			solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
//															solution_container[FLOW_SOL]->node[jPoint]->GetSolution());
//			
//			/*--- Laminar Viscosity ---*/
//			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
//																	solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
//			/*--- Eddy Viscosity ---*/
//			solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
//															 solution_container[FLOW_SOL]->node[jPoint]->GetEddyViscosity());
//			
//			/*--- Transition variables w/o reconstruction, and its gradients ---*/
//			solver->SetTransVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
//			solver->SetTransVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
//			
//		  solver->SetConsVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient(),
//																 solution_container[FLOW_SOL]->node[jPoint]->GetGradient());
//			
//			
//			/*--- Compute residual, and Jacobians ---*/
//			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
//			
//			/*--- Add and subtract residual, and update Jacobians ---*/
//			node[iPoint]->SubtractResidual(Residual);
//      node[jPoint]->AddResidual(Residual);
//
//      Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
//      Jacobian.SubtractBlock(iPoint,jPoint,Jacobian_j);
//      Jacobian.AddBlock(jPoint,iPoint,Jacobian_i);
//      Jacobian.AddBlock(jPoint,jPoint,Jacobian_j);
//
//		}
//	}
}

void CTransLMSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CNumerics *second_solver,
																			 CConfig *config, unsigned short iMesh) {
  unsigned long iPoint;
  double gamma_sep;

  //cout << "Setting Trans residual -AA " << endl;
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
   // cout << "iPoint: " << iPoint << endl;
		
		/*--- Conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);
		
		/*--- Gradient of the primitive and conservative variables ---*/
		solver->SetPrimVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
		
		/*--- Laminar viscosity ---*/
		solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
    solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),0.0);
		
		/*--- Turbulent variables w/o reconstruction, and its gradient ---*/
		solver->SetTransVar(node[iPoint]->GetSolution(), NULL);
		solver->SetTransVarGradient(node[iPoint]->GetGradient(), NULL);
		
		/*--- Set volume ---*/
		solver->SetVolume(geometry->node[iPoint]->GetVolume());
		
		/*--- Set distance to the surface ---*/
		solver->SetDistance(geometry->node[iPoint]->GetWallDistance(), 0.0);
		
		/*--- Compute the source term ---*/
		solver->SetResidual_TransLM(Residual, Jacobian_i, NULL, config, gamma_sep);
		
    /*-- Store gamma_sep in variable class --*/
    node[iPoint]->SetGammaSep(gamma_sep);

		/*--- Subtract residual and the jacobian ---*/
		node[iPoint]->SubtractResidual(Residual);
    Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);

	}
}

void CTransLMSolution::Source_Template(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																			 CConfig *config, unsigned short iMesh) {
}

void CTransLMSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
  double *U_i;
  double *U_domain = new double[nVar];
  double *U_wall   = new double[nVar];
  double *Normal   = new double[nDim];
  double *Residual = new double[nVar];

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT); 

//  cout << "Setting wall BC -AA\n";
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      /*--- Set both interior and exterior point to current value ---*/
      for (iVar=0; iVar < nVar; iVar++) {
        U_domain[iVar] = node[iPoint]->GetSolution(iVar);
        U_wall[iVar]   = node[iPoint]->GetSolution(iVar); 	
      }

      /*--- Set various quantities in the solver class ---*/
      solver->SetNormal(Normal);
      solver->SetTransVar(U_domain,U_wall);
  		U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
  		solver->SetConservative(U_i, U_i);

      /*--- Compute the residual using an upwind scheme ---*/
//      cout << "BC calling SetResidual: -AA" << endl;
      solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
//      cout << "Returned from BC call of SetResidual: -AA" << endl;
//      cout << "Residual[0] = " << Residual[0] << endl;
//      cout << "Residual[1] = " << Residual[1] << endl;
      node[iPoint]->AddResidual(Residual);

//      cout << "Implicit part -AA" << endl;
			/*--- Jacobian contribution for implicit integration ---*/
			if (implicit) {
        Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
      }
  }
}
}

void CTransLMSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
  double *U_FF   = new double[nVar];
  double *U_domain = new double[nVar];
  double *Residual = new double[nVar];
  double *Normal = new double[nDim];
  double FlowSolution_i[nDim+2], FlowSolution_j[nDim+2];

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
  //  cout << "Setting FF BC -AA\n";
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //  cout << "iVertex="<<iVertex<<", " << "iPoint=" << iPoint <<endl;
    if (geometry->node[iPoint]->GetDomain()) {


      //cout << "Get flow solution" <<endl;
			/*--- Set conservative variables at the wall, and at infinity ---*/
			for (unsigned short iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);

			FlowSolution_j[0] = solution_container[FLOW_SOL]->GetDensity_Inf();
			FlowSolution_j[nDim+1] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
			for (unsigned short iDim = 0; iDim < nDim; iDim++)
				FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(iDim);

      //cout << "Set flow solution" <<endl;
			conv_solver->SetConservative(FlowSolution_i, FlowSolution_j);

      /*--- Get this node's solution and freestream values, tell them to solver ---*/
      U_domain[0] = node[iPoint]->GetSolution(0);
      U_domain[1] = node[iPoint]->GetSolution(1);
      U_FF[0] = Intermittency_Inf;
      U_FF[1] = REth_Inf;
   //   cout << "itmc: " << U_domain[0] << ", " << U_FF[0] << endl;
   //   cout << "Reth: " << U_domain[1] << ", " << U_FF[1] << endl;
      conv_solver->SetTransVar(U_domain,U_FF);

      //cout << "Get normal" <<endl;
      /*--- Get normal and reverse it due to outward facing convention---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (unsigned short idim=0;idim<nDim;idim++)
        Normal[idim] = -Normal[idim];
      conv_solver->SetNormal(Normal);

      conv_solver->SetResidual(Residual, Jacobian_i, NULL, config);
      node[iPoint]->AddResidual(Residual);

      if (implicit) {
    	Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
      }
    }

  }
}

void CTransLMSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config,
																unsigned short val_marker) {
//cout << "BC_Inlet reached -AA" << endl;
//BC_Far_Field(geometry, solution_container,solver,config,val_marker);
}

void CTransLMSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver,
																 CConfig *config, unsigned short val_marker) {
//cout << "BC_Outlet reached -AA" << endl;
//BC_Far_Field(geometry, solution_container,solver,config,val_marker);
}

void CTransLMSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																 CConfig *config, unsigned short val_marker) {
//cout << "BC_Outlet reached -AA" << endl;
BC_NS_Wall(geometry, solution_container,solver,config,val_marker);
}
