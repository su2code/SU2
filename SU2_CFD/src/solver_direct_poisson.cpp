/*!
 * \file solution_direct_poisson.cpp
 * \brief Main subrotuines for solving direct problems
 * \author F. Palacios
 * \version 4.2.0 "Cardinal"
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
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

CPoissonSolver::CPoissonSolver(void) : CSolver() { }

CPoissonSolver::CPoissonSolver(CGeometry *geometry, CConfig *config) : CSolver() {
  
	unsigned long nPoint, iPoint;
	unsigned short iVar, iDim;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
	nDim =          geometry->GetnDim();
  nPoint =        geometry->GetnPoint();
  nPointDomain =  geometry->GetnPointDomain();
	nVar =          1;
	node =          new CVariable*[nPoint];
  
	Residual = new su2double[nVar]; Residual_RMS = new su2double[nVar];
	Solution = new su2double[nVar];
  Residual_Max = new su2double[nVar];
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
	/*--- Point to point stiffness matrix ---*/
  
	if (nDim == 2) {
		StiffMatrix_Elem = new su2double* [3];
		Source_Vector	 = new su2double [3];
		for (unsigned short iVar = 0; iVar < 3; iVar++) {
			StiffMatrix_Elem[iVar] = new su2double [3];
		}
	}
  
	if (nDim == 3) {
		StiffMatrix_Elem = new su2double* [4];
		Source_Vector	 = new su2double [4];
		for (unsigned short iVar = 0; iVar < 4; iVar++) {
			StiffMatrix_Elem[iVar] = new su2double [4];
		}
	}
  
	StiffMatrix_Node = new su2double* [1];
	for (unsigned short iVar = 0; iVar < 1; iVar++) {
		StiffMatrix_Node[iVar] = new su2double [1];
	}
  
	/*--- Initialization of the structure of the whole Jacobian ---*/
  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Poisson equation)." << endl;
	StiffMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  
  /*--- Solution and residual vectors ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);

	/*--- Computation of gradients by least squares ---*/
  
	Smatrix = new su2double* [nDim]; // S matrix := inv(R)*traspose(inv(R))
	for (iDim = 0; iDim < nDim; iDim++)
		Smatrix[iDim] = new su2double [nDim];
  
	cvector = new su2double* [nVar]; // c vector := transpose(WA)*(Wb)
	for (iVar = 0; iVar < nVar; iVar++)
		cvector[iVar] = new su2double [nDim];
  
	/*--- Always instantiate and initialize the variable to a zero value. ---*/
  
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		node[iPoint] = new CPotentialVariable(0.0, nDim, nVar, config);
  
}

CPoissonSolver::~CPoissonSolver(void) {
  
	unsigned short iVar;
  
	delete [] Source_Vector;
  
	if (nDim == 2) {
		for (iVar = 0; iVar < 3; iVar++)
			delete [] StiffMatrix_Elem[iVar];
	}
  
	if (nDim == 3) {
		for (iVar = 0; iVar < 4; iVar++)
			delete [] StiffMatrix_Elem[iVar];
	}
  
	for (iVar = 0; iVar < 1; iVar++)
		delete [] StiffMatrix_Node[iVar];
  
	delete [] StiffMatrix_Elem;
	delete [] StiffMatrix_Node;
  
}

void CPoissonSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
	unsigned long iPoint;
  
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		LinSysRes.SetBlock_Zero(iPoint);
  
	StiffMatrix.SetValZero();
  
}

void CPoissonSolver::Compute_Residual(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                      unsigned short iMesh) {
  
	unsigned long iPoint;
	unsigned short iVar = 0;
  
	/*--- Build linear system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		LinSysRes[iPoint] = LinSysRes.GetBlock(iPoint, iVar);
		LinSysSol[iPoint] = node[iPoint]->GetSolution(iVar);
	}
  
	StiffMatrix.MatrixVectorProduct(LinSysSol, LinSysRes);
  
	/*--- Update residual ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		LinSysRes.SetBlock(iPoint, 0, LinSysRes[iPoint]);
	}
}

/*!
 * \method Source_Residual
 * \brief Source terms of the poisson solver
 * \author A. Lonkar
 */
void CPoissonSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                     CConfig *config, unsigned short iMesh) {
//  
//	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
//	su2double a[3], b[3], c[3], d[3], Area_Local, Volume_Local;
//	//	su2double Local_Delta_Time;
//	su2double **Gradient_0, **Gradient_1, **Gradient_2, **Gradient_3;
//	su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;;
//	unsigned short iDim;
//	su2double  dt;
//	//	su2double  dx, u, c;
//	bool MacCormack_relaxation = (config->GetMacCormackRelaxation());
//  
//	if (nDim == 2) {
//		if (config->GetPoissonSolver()) {
//      
//			for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
//        
//				Point_0 = geometry->elem[iElem]->GetNode(0);
//				Point_1 = geometry->elem[iElem]->GetNode(1);
//				Point_2 = geometry->elem[iElem]->GetNode(2);
//        
//				Coord_0 = geometry->node[Point_0]->GetCoord();
//				Coord_1 = geometry->node[Point_1]->GetCoord();
//				Coord_2 = geometry->node[Point_2]->GetCoord();
//
//				for (iDim = 0; iDim < nDim; iDim++) {
//					a[iDim] = Coord_0[iDim]-Coord_2[iDim];
//					b[iDim] = Coord_1[iDim]-Coord_2[iDim];
//				}
//        
//				Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
//        
//				Gradient_0 = node[Point_0]->GetPlasmaRhoUGradient();
//				Gradient_1 = node[Point_1]->GetPlasmaRhoUGradient();
//				Gradient_2 = node[Point_2]->GetPlasmaRhoUGradient();
//        
//				numerics->SetVolume(Area_Local);
//        
//				dt = node[Point_0]->GetPlasmaTimeStep();
//        
//				/*		u = 4800;
//				 c = 87110;
//				 c = 800;
//				 dx = 0.004/81;
//				 Local_Delta_Time = config->GetCFL(iMesh) * dx/(u+c);
//				 numerics->SetTimeStep(Local_Delta_Time);
//				 */
//        
//				numerics->SetCoord(Coord_0, Coord_1, Coord_2);
//				numerics->SetTimeStep(dt);
//				numerics->SetChargeDensity(node[Point_0]->GetChargeDensity(), node[Point_1]->GetChargeDensity(), node[Point_2]->GetChargeDensity(), node[Point_3]->GetChargeDensity());
//				numerics->SetConsVarGradient(Gradient_0, Gradient_1, Gradient_2 );
//				numerics->ComputeResidual_MacCormack(Source_Vector, config);
//        
//				LinSysRes.AddBlock(Point_0, &Source_Vector[0]);
//				LinSysRes.AddBlock(Point_1, &Source_Vector[1]);
//				LinSysRes.AddBlock(Point_2, &Source_Vector[2]);
//
//				if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
//          
//					Point_0 = geometry->elem[iElem]->GetNode(3);
//					Point_1 = geometry->elem[iElem]->GetNode(0);
//					Point_2 = geometry->elem[iElem]->GetNode(2);
//          
//					Coord_0 = geometry->node[Point_0]->GetCoord();
//					Coord_1 = geometry->node[Point_1]->GetCoord();
//					Coord_2 = geometry->node[Point_2]->GetCoord();
//
//					for (iDim = 0; iDim < nDim; iDim++) {
//						a[iDim] = Coord_0[iDim]-Coord_2[iDim];
//						b[iDim] = Coord_1[iDim]-Coord_2[iDim];
//					}
//          
//					Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
//          
//					Gradient_0 = node[Point_0]->GetPlasmaRhoUGradient();
//					Gradient_1 = node[Point_1]->GetPlasmaRhoUGradient();
//					Gradient_2 = node[Point_2]->GetPlasmaRhoUGradient();
//          
//					numerics->SetVolume(Area_Local);
//          
//					/*		u = 4800;
//					 c = 87110;
//					 c = 732.0;
//					 dx = 0.004/81;
//					 Local_Delta_Time = config->GetCFL(iMesh) * dx/(u+c);
//					 numerics->SetTimeStep(Local_Delta_Time);
//					 */
//          
//					dt = node[Point_0]->GetPlasmaTimeStep();
//					numerics->SetCoord(Coord_0, Coord_1, Coord_2);
//					numerics->SetTimeStep(dt);
//					numerics->SetChargeDensity(node[Point_0]->GetChargeDensity(), node[Point_1]->GetChargeDensity(), node[Point_2]->GetChargeDensity(), node[Point_3]->GetChargeDensity());
//					numerics->SetConsVarGradient(Gradient_0, Gradient_1, Gradient_2 );
//					numerics->ComputeResidual_MacCormack(Source_Vector, config);
//					LinSysRes.AddBlock(Point_0, &Source_Vector[0]);
//					LinSysRes.AddBlock(Point_1, &Source_Vector[1]);
//					LinSysRes.AddBlock(Point_2, &Source_Vector[2]);
//				}
//			}
//		}
//	}
//	if (nDim == 3) {
//		if (config->GetPoissonSolver()) {
//			for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
//				Point_0 = geometry->elem[iElem]->GetNode(0);	Coord_0 = geometry->node[Point_0]->GetCoord();
//				Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1 = geometry->node[Point_1]->GetCoord();
//				Point_2 = geometry->elem[iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord();
//				Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3 = geometry->node[Point_3]->GetCoord();
//        
//				for (iDim = 0; iDim < nDim; iDim++) {
//					a[iDim] = Coord_0[iDim]-Coord_2[iDim];
//					b[iDim] = Coord_1[iDim]-Coord_2[iDim];
//					c[iDim] = Coord_3[iDim]-Coord_2[iDim];
//				}
//        
//				d[0] = a[1]*b[2]-a[2]*b[1];
//				d[1] = -(a[0]*b[2]-a[2]*b[0]);
//				d[2] = a[0]*b[1]-a[1]*b[0];
//        
//				/*--- Compute element volume ---*/
//				Volume_Local = fabs(c[0]*d[0] + c[1]*d[1] + c[2]*d[2])/6.0;
//				numerics->SetVolume(Volume_Local);
//				numerics->SetChargeDensity(node[Point_0]->GetChargeDensity(), node[Point_1]->GetChargeDensity(), node[Point_2]->GetChargeDensity(), node[Point_3]->GetChargeDensity());
//        
//				if (MacCormack_relaxation) {
//          
//					Gradient_0 = node[Point_0]->GetPlasmaRhoUGradient();
//					Gradient_1 = node[Point_1]->GetPlasmaRhoUGradient();
//					Gradient_2 = node[Point_2]->GetPlasmaRhoUGradient();
//					Gradient_3 = node[Point_3]->GetPlasmaRhoUGradient();
//					numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
//					numerics->SetConsVarGradient(Gradient_0, Gradient_1, Gradient_2, Gradient_3 );
//					numerics->ComputeResidual_MacCormack(Source_Vector, config);
//				}
//				else numerics->ComputeResidual(Source_Vector, config);
//        
//				LinSysRes.AddBlock(Point_0, &Source_Vector[0]);
//				LinSysRes.AddBlock(Point_1, &Source_Vector[1]);
//				LinSysRes.AddBlock(Point_2, &Source_Vector[2]);
//				LinSysRes.AddBlock(Point_3, &Source_Vector[3]);
//			}
//		}
//	}
}

void CPoissonSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {
}


/*!
 * \method Copy_Zone_Solution
 * \brief Copy solution from solver 1 into solver 2
 * \author A. Lonkar
 */
void CPoissonSolver::Copy_Zone_Solution(CSolver ***solver1_solution,
                                        CGeometry **solver1_geometry,
                                        CConfig *solver1_config,
                                        CSolver ***solver2_solution,
                                        CGeometry **solver2_geometry,
                                        CConfig *solver2_config) {
	unsigned long iPoint;
	unsigned short iDim;
	su2double neg_EFvalue;
	su2double *E_field = new su2double [nDim];
  
	for (iPoint = 0; iPoint < solver1_geometry[MESH_0]->GetnPointDomain(); iPoint++) {
		for (iDim =0; iDim < nDim; iDim ++) {
			neg_EFvalue = solver1_solution[MESH_0][POISSON_SOL]->node[iPoint]->GetGradient(0, iDim);
			E_field[iDim] = -1.0*neg_EFvalue;
		}
	}
};

/*!
 * \method Galerkin_Method
 * \brief calculate the element stiffness matrix
 * \author A. Lonkar
 */
void CPoissonSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
	su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3 = NULL;
  
	if (nDim == 2 ) {
		for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      
			Point_0 = geometry->elem[iElem]->GetNode(0);
			Point_1 = geometry->elem[iElem]->GetNode(1);
			Point_2 = geometry->elem[iElem]->GetNode(2);
      
			Coord_0 = geometry->node[Point_0]->GetCoord();
			Coord_1 = geometry->node[Point_1]->GetCoord();
			Coord_2 = geometry->node[Point_2]->GetCoord();
      
			numerics->SetCoord(Coord_0, Coord_1, Coord_2);
			numerics->ComputeResidual(StiffMatrix_Elem, config);
			AddStiffMatrix(StiffMatrix_Elem, Point_0, Point_1, Point_2, Point_3);
		}
	}
  
	if (nDim == 3 ) {
    
		for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      
      Point_0 = geometry->elem[iElem]->GetNode(0); 	Coord_0 = geometry->node[Point_0]->GetCoord();
      Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1 = geometry->node[Point_1]->GetCoord();
      Point_2 = geometry->elem[iElem]->GetNode(2); 	Coord_2 = geometry->node[Point_2]->GetCoord();
      Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3 = geometry->node[Point_3]->GetCoord();
      
      numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
      numerics->ComputeResidual(StiffMatrix_Elem, config);
      AddStiffMatrix(StiffMatrix_Elem, Point_0, Point_1, Point_2, Point_3);
      
		}
	}
}

/*!
 * \method AddStiffMatrix
 * \brief Assemble the Global stiffness matrix
 * \author A. Lonkar
 */
void CPoissonSolver::AddStiffMatrix(su2double **StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1, unsigned long Point_2, unsigned long Point_3) {
  
	if (nDim == 2 ) {
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]; StiffMatrix.AddBlock(Point_0, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]; StiffMatrix.AddBlock(Point_0, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]; StiffMatrix.AddBlock(Point_0, Point_2, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]; StiffMatrix.AddBlock(Point_1, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]; StiffMatrix.AddBlock(Point_1, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]; StiffMatrix.AddBlock(Point_1, Point_2, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]; StiffMatrix.AddBlock(Point_2, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]; StiffMatrix.AddBlock(Point_2, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]; StiffMatrix.AddBlock(Point_2, Point_2, StiffMatrix_Node);
	}
	if (nDim == 3) {
    
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]; StiffMatrix.AddBlock(Point_0, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]; StiffMatrix.AddBlock(Point_0, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]; StiffMatrix.AddBlock(Point_0, Point_2, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][3]; StiffMatrix.AddBlock(Point_0, Point_3, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]; StiffMatrix.AddBlock(Point_1, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]; StiffMatrix.AddBlock(Point_1, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]; StiffMatrix.AddBlock(Point_1, Point_2, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][3]; StiffMatrix.AddBlock(Point_1, Point_3, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]; StiffMatrix.AddBlock(Point_2, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]; StiffMatrix.AddBlock(Point_2, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]; StiffMatrix.AddBlock(Point_2, Point_2, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][3]; StiffMatrix.AddBlock(Point_2, Point_3, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][0]; StiffMatrix.AddBlock(Point_3, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][1]; StiffMatrix.AddBlock(Point_3, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][2]; StiffMatrix.AddBlock(Point_3, Point_2, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][3]; StiffMatrix.AddBlock(Point_3, Point_3, StiffMatrix_Node);
    
	}
}

void CPoissonSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container,
                                  CConfig *config, unsigned short val_marker) {
  unsigned long Point, iVertex;
  
	/*--- Identify if a boundary is Dirichlet or Neumman ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    Point = geometry->vertex[val_marker][iVertex]->GetNode();
    Solution[0]= 10.0;
    node[Point]->SetSolution(Solution);

    LinSysRes.SetBlock(Point, Solution);
    LinSysSol.SetBlock(Point, Solution);

    StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
  }
  
}

void CPoissonSolver::BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                unsigned short val_marker) { }

void CPoissonSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
	unsigned long iPoint, total_index;
  unsigned short iVar;
	
	/*--- Build implicit system ---*/
  
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			LinSysSol[total_index] = 0.0;
		}
    
	}
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }
	
  /*--- Solve or smooth the linear system ---*/
  
  CSysSolve system;
  system.Solve(StiffMatrix, LinSysRes, LinSysSol, geometry, config);
  
	/*--- Update solution (system written in terms of increments) ---*/
  
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->SetSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
		}
	}
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*---  Compute the residual Ax-f ---*/
  
	StiffMatrix.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);
  
  /*--- Set maximum residual to zero ---*/
  
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Compute the residual ---*/
  
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
      AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
		}
	}
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}
