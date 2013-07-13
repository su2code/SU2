/*!
 * \file solution_direct_heat.cpp
 * \brief Main subrotuines for solving the heat equation.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.5
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

CHeatSolution::CHeatSolution(void) : CSolution() { }

CHeatSolution::CHeatSolution(CGeometry *geometry, 
                             CConfig *config) : CSolution() {
  
	unsigned long nPoint;
	unsigned short nMarker, iVar;
  
	nPoint  = geometry->GetnPoint();
	nDim    = geometry->GetnDim();
	nMarker = config->GetnMarker_All(); 
	node    = new CVariable*[nPoint];
	nVar    = 2; // solve as a 2 eq. system		
  
	Residual     = new double[nVar]; Residual_RMS = new double[nVar];
	Solution     = new double[nVar];
  Res_Sour     = new double[nVar];
  Residual_Max = new double[nVar]; Point_Max = new unsigned long[nVar];
  

	/*--- Point to point stiffness matrix (only for triangles)---*/	
	StiffMatrix_Elem = new double*[nDim+1];
	for (iVar = 0; iVar < nDim+1; iVar++) {
		StiffMatrix_Elem[iVar] = new double [nDim+1];
	}
  
	StiffMatrix_Node = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		StiffMatrix_Node[iVar] = new double [nVar];
	}
  
	/*--- Initialization of matrix structures ---*/
	Initialize_SparseMatrix_Structure(&StiffMatrixSpace, nVar, nVar, geometry, config);
	Initialize_SparseMatrix_Structure(&StiffMatrixTime, nVar, nVar, geometry, config);
	Initialize_SparseMatrix_Structure(&Jacobian, nVar, nVar, geometry, config);
  
  /*--- Initialization of linear solver structures ---*/
	xsol = new double [nPoint*nVar];
	xres = new double [nPoint*nVar];
  
  /* Heat strength coefficient for all of the markers */
  
  CHeat = new double[config->GetnMarker_All()];
  Total_CHeat = 0.0;
  
  /* Check for a restart (not really used), initialize from zero otherwise */
  
	bool restart = (config->GetRestart());
	if (!restart) {
    
    double *Coord, Radius;
		for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
      
      /*--- Set up the initial condition for the drum problem ---*/
      Coord = geometry->node[iPoint]->GetCoord();
      Radius = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        Radius += Coord[iDim]*Coord[iDim];
      Radius = sqrt(Radius);
      
      /*--- Symmetrically plucked drum ---*/
      //      Solution[0] = 1.0 - Radius; 
      //      Solution[1] = 0.0;
      
      /*--- Off-center strike ---*/
			//      if ((Radius > 0.4) && (Radius < 0.6)) {
			//        Solution[0] = 10.0; 
			//      } else
			//          Solution[0] = 0.0;
			//      Solution[1] = 0.0;
      
      /*--- Struck drum ---*/
			//      Solution[0] = 0.0; 
			//      Solution[1] = -1.0;
      
      /*--- Zero initial condition for testing source terms & forcing BCs ---*/
      Solution[0] = 0.0; 
      Solution[1] = 0.0;
      
			node[iPoint] = new CHeatVariable(Solution, nDim, nVar, config);
      
      /* Copy solution to old containers if using dual time */
      
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST) {
        node[iPoint]->Set_Solution_time_n(); 
      } else if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) {
        node[iPoint]->Set_Solution_time_n();
        node[iPoint]->Set_Solution_time_n1();
      }
      
    }
	} else {
    
    cout << "Heat restart file not currently configured!!" << endl;
    cout << "Press any key to exit..." << endl;
    cin.get();
    exit(1);
    
		string mesh_filename = config->GetSolution_FlowFileName();
		ifstream restart_file;
    
		char *cstr; cstr = new char [mesh_filename.size()+1];
		strcpy (cstr, mesh_filename.c_str());
		restart_file.open(cstr, ios::in);
        
		if (restart_file.fail()) {
			cout << "There is no Heat restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}
		unsigned long index;
		string text_line;
    
		for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			point_line >> index >> Solution[0] >> Solution[1];
			node[iPoint] = new CHeatVariable(Solution, nDim, nVar, config);
		}
		restart_file.close();
	}
  
  /* Load all of the matrices and Jacobians that are fixed */
  /* during the simulation in order to save effort.        */
  //SetTime_Matrix(geometry, config);
  
}

CHeatSolution::~CHeatSolution(void) {
  
	unsigned short iVar;
  
	delete [] Residual;
	delete [] Residual_Max;
	delete [] Solution;
  
	for (iVar = 0; iVar < nDim+1; iVar++)
		delete [] StiffMatrix_Elem[iVar];
  
  
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] StiffMatrix_Node[iVar];
  
	delete [] StiffMatrix_Elem;
	delete [] StiffMatrix_Node;
	
	delete [] xsol;
	delete [] xres;

}

void CHeatSolution::Preprocessing(CGeometry *geometry, 
                                  CSolution **solution_container,
                                  CConfig   *config, 
                                  unsigned short iMesh,
                                  unsigned short iRKStep,
                                  unsigned short RunTime_EqSystem) {
  
  /* Set residuals and matrix entries to zero */
  
	for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    Set_Residual_Zero(iPoint);
	}
	
  /* Zero out the entries in the various matrices */
  
	StiffMatrixSpace.SetValZero();
	StiffMatrixTime.SetValZero();
	Jacobian.SetValZero();
  
}

void CHeatSolution::Source_Residual(CGeometry *geometry, 
                                    CSolution **solution_container, 
                                    CNumerics *solver, CNumerics *second_solver,
                                    CConfig   *config, 
                                    unsigned short iMesh) { }


void CHeatSolution::Source_Template(CGeometry *geometry,
                                    CSolution **solution_container,
                                    CNumerics *solver,
                                    CConfig   *config,
                                    unsigned short iMesh) { }

void CHeatSolution::Galerkin_Method(CGeometry *geometry, 
                                    CSolution **solution_container, 
                                    CNumerics *solver,
                                    CConfig   *config, 
                                    unsigned short iMesh) {
  
  /* Local variables and initialization */
  
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, iPoint, total_index;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, Thermal_Diffusivity ;
	unsigned short iVar;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
    
		Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
    
		solver->SetCoord(Coord_0, Coord_1, Coord_2);
    
		solver->SetResidual(StiffMatrix_Elem, config);
    
    /*--- Compute the square of the Heat speed ---*/
		Thermal_Diffusivity  = config->GetThermalDiffusivity();
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[0][0]*Thermal_Diffusivity ; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_0, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[0][1]*Thermal_Diffusivity ; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_0, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[0][2]*Thermal_Diffusivity ; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_0, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][0]*Thermal_Diffusivity ; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_1, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][1]*Thermal_Diffusivity ; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_1, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][2]*Thermal_Diffusivity ; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_1, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[2][0]*Thermal_Diffusivity ; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_2, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[2][1]*Thermal_Diffusivity ; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_2, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[2][2]*Thermal_Diffusivity ; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_2, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
    
	}
	
  /* Prepare solution vector for multiplication */
  
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			xsol[total_index] = node[iPoint]->GetSolution(iVar);
			xres[total_index] = 0.0;
		}
	
	StiffMatrixSpace.MatrixVectorProduct(xsol, xres);
	
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			Residual[iVar] = xres[total_index];
		}
		SubtractResidual(iPoint, Residual);
	}
}


void CHeatSolution::BC_Euler_Wall(CGeometry *geometry, 
                                  CSolution **solution_container, 
                                  CNumerics *solver, 
                                  CConfig   *config, 
																	unsigned short val_marker) {
	
  /* Local variables */
  
  unsigned long iPoint, iVertex, total_index, iter;
	double deltaT, omega, time, ampl, Heat_sol[2];
  
  /* Set the values needed for periodic forcing */
  
  deltaT = config->GetDelta_UnstTimeND();
  iter   = config->GetExtIter();  
  time   = static_cast<double>(iter)*deltaT;
  omega  = 1.0;
  ampl   = 1.0;
  
    /* Compute sin Heat forcing at the boundary */
  
  Heat_sol[0] = ampl*sin(omega*time);
  Heat_sol[1] = 0.0; //-ampl*cos(omega*time_new)*omega;
  
  /* Set the solution at the boundary nodes and zero the residual */
	
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
			Solution[iVar] = Heat_sol[iVar];
			Residual[iVar] = 0.0;
		}
		node[iPoint]->SetSolution(Solution);
		node[iPoint]->SetSolution_Old(Solution);
    SetResidual(iPoint, Residual);
		SetResidual(iPoint, Residual);
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      Jacobian.DeleteValsRowi(total_index);
    }
	}
  
}

void CHeatSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver,
                                 CNumerics *visc_solver, CConfig *config, unsigned short val_marker) { }

void CHeatSolution::SetResidual_DualTime(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep, 
                                        unsigned short iMesh, unsigned short RunTime_EqSystem) {
	
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double a[3], b[3], Area_Local, Time_Num;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
	unsigned short iDim, iVar, jVar;
	double TimeJac = 0.0;
	
	/*--- Numerical time step (this system is uncoditional stable... a very big number can be used) ---*/
	Time_Num = config->GetDelta_UnstTimeND();
	
	/*--- Loop through elements to compute contributions from the matrix
   blocks involving time. These contributions are also added to the 
   Jacobian w/ the time step. Spatial source terms are also computed. ---*/
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
    /* Get node numbers and their coordinate vectors */
    
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		
    Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		
    /* Compute triangle area (2-D) */
		
    for (iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
		Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

		
		/*----------------------------------------------------------------*/
		/*--- Block contributions to the Jacobian (includes time step) ---*/
		/*----------------------------------------------------------------*/
		
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				StiffMatrix_Node[iVar][jVar] = 0.0;	
		
		if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST) TimeJac = 1.0/Time_Num;
		if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) TimeJac = 3.0/(2.0*Time_Num);
		
		/*--- Diagonal value identity matrix ---*/
			StiffMatrix_Node[0][0] = 1.0*TimeJac; 

		
		/*--- Diagonal value ---*/
			StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local*TimeJac);

    Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_1, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_2, StiffMatrix_Node);
		
		/*--- Off Diagonal value ---*/
			StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local*TimeJac);

		Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_1, StiffMatrix_Node);
		Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_2, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_2, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_1, StiffMatrix_Node);

    
	}
	
	unsigned long iPoint, total_index;
	double *U_time_nM1, *U_time_n, *U_time_nP1;
	
	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 
		
		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n   = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();
		
		/*--- Compute Residual ---*/
		for(iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			xres[total_index] = 0.0;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				xsol[total_index] = (U_time_nP1[iVar] - U_time_n[iVar]);
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				xsol[total_index] = (U_time_nP1[iVar] - (4.0/3.0)*U_time_n[iVar] +  (1.0/3.0)*U_time_nM1[iVar]);
		}
	}
	
	StiffMatrixTime.MatrixVectorProduct(xsol, xres);		
  
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			Residual[iVar] = xres[total_index];
		}
		SubtractResidual(iPoint, Residual);
	}
	
}


void CHeatSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
    unsigned short iVar;
	unsigned long iPoint, total_index;
    
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
    }
	
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
        
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			xsol[total_index] = 0.0;
			AddRes_RMS(iVar, xres[total_index]*xres[total_index]);
            AddRes_Max(iVar, fabs(xres[total_index]), geometry->node[iPoint]->GetGlobalIndex());
		}
	}
    
    /*--- Initialize residual and solution at the ghost points ---*/
    for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
            total_index = iPoint*nVar + iVar;
            xres[total_index] = 0.0;
            xsol[total_index] = 0.0;
        }
    }
    
	/*--- Solve the linear system (Stationary iterative methods) ---*/
	if (config->GetKind_Linear_Solver() == SYM_GAUSS_SEIDEL)
		Jacobian.SGSSolution(xres, xsol, config->GetLinear_Solver_Error(),
                             config->GetLinear_Solver_Iter(), false, geometry, config);
	
	if (config->GetKind_Linear_Solver() == LU_SGS)
        Jacobian.LU_SGSIteration(xres, xsol, geometry, config);
	
	/*--- Solve the linear system (Krylov subspace methods) ---*/
	if ((config->GetKind_Linear_Solver() == BCGSTAB) ||
        (config->GetKind_Linear_Solver() == GMRES)) {
		
		CSysVector rhs_vec((const unsigned int)geometry->GetnPoint(),
                           (const unsigned int)geometry->GetnPointDomain(), nVar, xres);
		CSysVector sol_vec((const unsigned int)geometry->GetnPoint(),
                           (const unsigned int)geometry->GetnPointDomain(), nVar, xsol);
		
		CMatrixVectorProduct* mat_vec = new CSparseMatrixVectorProduct(Jacobian, geometry, config);

		CPreconditioner* precond = NULL;
		if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
			Jacobian.BuildJacobiPreconditioner();
			precond = new CJacobiPreconditioner(Jacobian, geometry, config);
		}
		else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
			Jacobian.BuildJacobiPreconditioner();
			precond = new CLineletPreconditioner(Jacobian, geometry, config);
		}
		else if (config->GetKind_Linear_Solver_Prec() == NO_PREC) {
			precond = new CIdentityPreconditioner(Jacobian, geometry, config);
        }
		
		CSysSolve system;
		if (config->GetKind_Linear_Solver() == BCGSTAB)
			system.BCGSTAB(rhs_vec, sol_vec, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                           config->GetLinear_Solver_Iter(), false);
		else if (config->GetKind_Linear_Solver() == GMRES)
			system.GMRES(rhs_vec, sol_vec, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                         config->GetLinear_Solver_Iter(), false);
		
		sol_vec.CopyToArray(xsol);
		delete mat_vec;
		delete precond;
	}
	
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*xsol[iPoint*nVar+iVar]);
		}
	}
	
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CHeatSolution::SetTime_Matrix(CGeometry *geometry, 
                                   CConfig   *config) {
  
  
  /* Local variables and initialization */
  
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double a[3], b[3], Area_Local, Time_Num, Time_Phys;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
	unsigned short iDim;
	
  /*--- Numerical time step. This system is unconditionally stable,
   so a very big step can be used. ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING) 
    Time_Num = config->GetDelta_UnstTimeND();
	else Time_Num = 1E+30;
  
  /*--- Physical timestep for source terms ---*/
  Time_Phys = config->GetDelta_UnstTimeND();
  
	/* Loop through elements to compute contributions from the matrix     */
  /* blocks involving time. These contributions are also added to the   */
  /* Jacobian w/ the time step. Spatial source terms are also computed. */
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
    /* Get node numbers and their coordinate vectors */
    
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		
    Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		
    /* Compute triangle area (2-D) */
		
    for (iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
		Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
		
    /*--- Block contributions to the Jacobian (includes time step) ---*/
    
		StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] =  1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] =  1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
    
  }
  
  //unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	//double a[3], b[3], Area_Local, Time_Num;
	//double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
	unsigned short iVar, jVar;
	double TimeJac = 0.0;
	
	/*--- Numerical time step (this system is uncoditional stable... a very big number can be used) ---*/
  if (config->GetUnsteady_Simulation() == TIME_STEPPING) 
    Time_Num = 1E+30;
	else Time_Num = config->GetDelta_UnstTimeND();
	
	/*--- Loop through elements to compute contributions from the matrix
   blocks involving time. These contributions are also added to the 
   Jacobian w/ the time step. Spatial source terms are also computed. ---*/
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
    /* Get node numbers and their coordinate vectors */
    
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		
    Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		
    /* Compute triangle area (2-D) */
		
    for (iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
		Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
    
		
		/*----------------------------------------------------------------*/
		/*--- Block contributions to the Jacobian (includes time step) ---*/
		/*----------------------------------------------------------------*/
		
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				StiffMatrix_Node[iVar][jVar] = 0.0;	
		
		if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST) TimeJac = 1.0/Time_Num;
		if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) TimeJac = 3.0/(2.0*Time_Num);
		
		/*--- Diagonal value identity matrix ---*/
    StiffMatrix_Node[0][0] = 1.0*TimeJac; 
    
		
		/*--- Diagonal value ---*/
    StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local*TimeJac);
    
    Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_1, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_2, StiffMatrix_Node);
		
		/*--- Off Diagonal value ---*/
    StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local*TimeJac);
    
		Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_1, StiffMatrix_Node);
		Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_2, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_2, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_1, StiffMatrix_Node);
    
	}
}

void CHeatSolution::GetRestart(CGeometry *geometry, CConfig *config) {
  
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif
  
  /*--- Restart the solution from file information ---*/
  string restart_filename = config->GetRestart_HeatFileName();
  unsigned long iPoint, index, nFlowIter, adjIter, flowIter;
  char buffer[50];
  string UnstExt, text_line;
  ifstream restart_file;
  
  /*--- For the unsteady adjoint, we integrate backwards through
   physical time, so load in the direct solution files in reverse. ---*/  
  if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    nFlowIter = config->GetnExtIter();
    adjIter   = config->GetExtIter();
    flowIter  = nFlowIter - adjIter - 1;
    restart_filename.erase (restart_filename.end()-4, restart_filename.end());
    if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
    if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
    if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
    if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
    if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
    UnstExt = string(buffer);
    restart_filename.append(UnstExt);
  } else {
    flowIter  =config->GetExtIter();
    restart_filename.erase (restart_filename.end()-4, restart_filename.end());
    if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
    if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
    if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
    if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
    if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
    UnstExt = string(buffer);
    restart_filename.append(UnstExt);
  }
  
  /*--- Open the flow solution from the restart file ---*/
  if (rank == MASTER_NODE)
    cout << "Reading in the direct heat solution from iteration " << flowIter << "." << endl;;
  restart_file.open(restart_filename.data(), ios::in);
  
  /*--- In case there is no file ---*/
  if (restart_file.fail()) {
    cout << "There is no Heat restart file!!" << endl;
    cout << "Press any key to exit..." << endl;
    cin.get(); exit(1);
  }
  
  /*--- Read the restart file ---*/
  for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    getline(restart_file,text_line);
    istringstream point_line(text_line);
    point_line >> index >> Solution[0] >> Solution[1];
    node[iPoint]->SetSolution_Direct(Solution);
  }
  
  /*--- Close the restart file ---*/
  restart_file.close();
  
}
