/*!
 * \file solution_direct_elasticity.cpp
 * \brief Main subrotuines for solving the linear elasticity equation.
 * \author F. Palacios, R. Sanchez
 * \version 4.0.1 "Cardinal"
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

CFEASolver::CFEASolver(void) : CSolver() { }

CFEASolver::CFEASolver(CGeometry *geometry, CConfig *config) : CSolver() {
  
  unsigned long iPoint;
  unsigned short iVar, jVar, NodesElement = 0, nLineLets;
  unsigned short iDim;
  su2double dull_val;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  nPoint        = geometry->GetnPoint();
  nPointDomain  = geometry->GetnPointDomain();
  nDim          = geometry->GetnDim();
  node          = new CVariable*[nPoint];
  
  
  WAitken_Dyn       = 0.0;
  WAitken_Dyn_tn1   = 0.0;
  
  SetFSI_ConvValue(0,0.0);
  SetFSI_ConvValue(1,0.0);
  
  nVar = nDim;
  
  if (nDim == 2) NodesElement = 4;
  if (nDim == 3) NodesElement = 8;
  
  /*--- Define some auxiliary vectors related to the residual ---*/
  
  Residual = new su2double[nVar];          for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS = new su2double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max = new su2double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  
  Solution   = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  
  /*--- Element aux stiffness matrix definition ---*/
  
  StiffMatrix_Elem = new su2double*[NodesElement*nDim];
  for (iVar = 0; iVar < NodesElement*nDim; iVar++) {
    StiffMatrix_Elem[iVar] = new su2double [NodesElement*nDim];
    for (jVar = 0; jVar < NodesElement*nDim; jVar++) {
      StiffMatrix_Elem[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Node aux stiffness matrix definition ---*/
  
  StiffMatrix_Node = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    StiffMatrix_Node[iVar] = new su2double [nVar];
    for (jVar = 0; jVar < nVar; jVar++) {
      StiffMatrix_Node[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Element aux mass matrix definition ---*/
  
  MassMatrix_Elem = new su2double*[NodesElement*nDim];
  for (iVar = 0; iVar < NodesElement*nDim; iVar++) {
    MassMatrix_Elem[iVar] = new su2double [NodesElement*nDim];
    for (jVar = 0; jVar < NodesElement*nDim; jVar++) {
      MassMatrix_Elem[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Node aux mass matrix definition ---*/
  
  MassMatrix_Node = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    MassMatrix_Node[iVar] = new su2double [nVar];
    for (jVar = 0; jVar < nVar; jVar++) {
      MassMatrix_Node[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Node aux mass matrix definition ---*/
  
  MassMatrix_Node_Int = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    MassMatrix_Node_Int[iVar] = new su2double [nVar];
    for (jVar = 0; jVar < nVar; jVar++) {
      MassMatrix_Node_Int[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Element aux damping matrix definition ---*/
  
  DampMatrix_Elem = new su2double*[NodesElement*nDim];
  for (iVar = 0; iVar < NodesElement*nDim; iVar++) {
    DampMatrix_Elem[iVar] = new su2double [NodesElement*nDim];
    for (jVar = 0; jVar < NodesElement*nDim; jVar++) {
      DampMatrix_Elem[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Node aux damping matrix definition ---*/
  
  DampMatrix_Node = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    DampMatrix_Node[iVar] = new su2double [nVar];
    for (jVar = 0; jVar < nVar; jVar++) {
      DampMatrix_Node[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Initialization of integration constants ---*/
  
  for (iVar = 0; iVar < 8; iVar++){
    a_dt[iVar]=0.0;
  }
  
  
  
  /*--- DESTRUCT THIS! ---*/
  
  /*--- Element aux dead load vector definition ---*/
  DeadLoadVector_Elem = new su2double [NodesElement*nDim];
  
  /*--- Node aux dead load vector definition ---*/
  DeadLoadVector_Node = new su2double [nVar];
  
  
  /*--- Initialization of matrix structures ---*/
  
  if (rank == MASTER_NODE) cout << "Initialize Stiffness structure (Linear Elasticity)." << endl;
  
  if (nDim==2){
    unsigned short form2d;
    form2d=config->GetElas2D_Formulation();
    if (form2d==0) cout << "Plane stress model for 2D structural analysis (Linear Elasticity)." << endl;
    if (form2d==1) cout << "Plane strain model for 2D structural analysis (Linear Elasticity)." << endl;
  }
  
  StiffMatrixSpace.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
  StiffMatrixTime.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
  MassMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
  DampMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
  
  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Linear Elasticity)." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
  
  if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
      (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
    nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
    if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
  }
  
  /*--- Initialization of linear solver structures ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  TimeRes_Aux.Initialize(nPoint, nPointDomain, nVar, 0.0);
  TimeRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  
  /*--- Computation of gradients by least squares ---*/
  
  Smatrix = new su2double* [nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Smatrix[iDim] = new su2double [nDim];
  
  cvector = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    cvector[iVar] = new su2double [nDim];
  
  /*--- Check for a restart, initialize from zero otherwise ---*/
  
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  
  if (!restart) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      node[iPoint] = new CFEAVariable(Solution, nDim, nVar, config);
    }
  }
  else {
    unsigned long index;
    string text_line, mesh_filename;
    ifstream restart_file;
    
    /*--- Restart the solution from file information ---*/
    
    mesh_filename = config->GetSolution_FlowFileName();
    restart_file.open(mesh_filename.data(), ios::in);
    
    /*--- In case there is no file ---*/
    
    if (restart_file.fail()) {
      cout << "There is no fea restart file!!" << endl;
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
        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2];
        node[iPoint_Local] = new CFEAVariable(Solution, nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    
    for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
      node[iPoint] = new CFEAVariable(Solution, nDim, nVar, config);
    }
    
    /*--- Close the restart file ---*/
    
    restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    
    delete [] Global2Local;
  }
  
}

CFEASolver::~CFEASolver(void) {
  
	unsigned short iVar, iDim, NodesElement = 0;
  
  if (nDim == 2) NodesElement = 4;
  if (nDim == 3) NodesElement = 8;
  
  delete [] Residual;
  delete [] Residual_Max;
  delete [] Solution;
  
  for (iVar = 0; iVar < NodesElement*nDim; iVar++)
    delete [] StiffMatrix_Elem[iVar];
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] StiffMatrix_Node[iVar];
  
  for (iVar = 0; iVar < NodesElement*nDim; iVar++)
    delete [] MassMatrix_Elem[iVar];
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] MassMatrix_Elem[iVar];
  
  for (iVar = 0; iVar < NodesElement*nDim; iVar++)
    delete [] DampMatrix_Elem[iVar];
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] DampMatrix_Elem[iVar];
  
  delete [] StiffMatrix_Elem;
  delete [] StiffMatrix_Node;
  delete [] MassMatrix_Elem;
  delete [] MassMatrix_Node;
  delete [] DampMatrix_Elem;
  delete [] DampMatrix_Node;
  
  /*--- Computation of gradients by least-squares ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] Smatrix[iDim];
  delete [] Smatrix;
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] cvector[iVar];
  delete [] cvector;
  
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] MassMatrix_Elem[iVar];

	for (iVar = 0; iVar < NodesElement*nDim; iVar++)
		delete [] DampMatrix_Elem[iVar];

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] DampMatrix_Elem[iVar];

	delete [] StiffMatrix_Elem;
	delete [] StiffMatrix_Node;
	delete [] MassMatrix_Elem;
	delete [] MassMatrix_Node;
	delete [] DampMatrix_Elem;
	delete [] DampMatrix_Node;
  
	/*--- Computation of gradients by least-squares ---*/
  
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;
  
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;

}



void CFEASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {
  GetSurface_Pressure(geometry, config);
  
  unsigned long ExtIter = config->GetExtIter();
  
  /*--- Set residuals and auxiliar variables to zero ---*/
  
  Initialize_SystemMatrix(geometry, solver_container, config);
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  if (ExtIter == 0){
    
    if (!dynamic){
      Compute_StiffMatrix(geometry, solver_container, numerics[VISC_TERM], config);
    }
    else if (dynamic){
      /*--- Compute the integration constants ---*/
      Compute_IntegrationConstants(geometry, solver_container, numerics[VISC_TERM], config);
      
      /*--- Compute the stiffness and mass matrices ---*/
      Compute_StiffMassMatrix(geometry, solver_container, numerics[VISC_TERM], config);
      
      //	  Compute_StiffMassDampMatrix(geometry, solver_container, numerics[VISC_TERM], config);
    }
    
  }
  
}

void CFEASolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

	unsigned long iPoint;
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
	unsigned long ExtIter = config->GetExtIter();
	bool fsi = config->GetFSI_Simulation();

	  if (!fsi) {

		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
			LinSysRes.SetBlock_Zero(iPoint);
			LinSysAux.SetBlock_Zero(iPoint);
		}

		/*--- Set matrix entries to zero ---*/

		/*--- Static calculation ---*/

		if (dynamic && ExtIter == 0){

			StiffMatrixSpace.SetValZero();
			StiffMatrixTime.SetValZero();

			MassMatrix.SetValZero();
			DampMatrix.SetValZero();

		}

		/*--- Dynamic calculation ---*/
		else if (!dynamic){

			StiffMatrixSpace.SetValZero();
			Jacobian.SetValZero();

		}

	  }

		else {

			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
				LinSysAux.SetBlock_Zero(iPoint);
			}

			/*--- Set matrix entries to zero ---*/

			/*--- Static calculation ---*/

			if (dynamic && ExtIter == 0){

				StiffMatrixSpace.SetValZero();
				StiffMatrixTime.SetValZero();

				MassMatrix.SetValZero();
				DampMatrix.SetValZero();

			}

			/*--- Dynamic calculation ---*/
			else if (!dynamic){

				StiffMatrixSpace.SetValZero();
				Jacobian.SetValZero();

			}

		}

}

void CFEASolver::Compute_IntegrationConstants(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
  
  su2double Delta_t= config->GetDelta_DynTime();
  su2double delta = config->GetNewmark_delta(), alpha = config->GetNewmark_alpha();
  
  /*--- Integration constants for Newmark scheme ---*/
  
  a_dt[0]= 1 / (alpha*pow(Delta_t,2.0));
  a_dt[1]= delta / (alpha*Delta_t);
  a_dt[2]= 1 / (alpha*Delta_t);
  a_dt[3]= 1 /(2*alpha) - 1;
  a_dt[4]= delta/alpha - 1;
  a_dt[5]= (Delta_t/2) * (delta/alpha - 2);
  a_dt[6]= Delta_t * (1-delta);
  a_dt[7]= delta * Delta_t;
  
}

void CFEASolver::Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
  
  unsigned short iVar, jVar, nNodes = 0, iNodes, iDim, jDim, form2d;
  unsigned long iElem, PointCorners[8];
  su2double CoordCorners[8][3];
  
  form2d=config->GetElas2D_Formulation();
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     nNodes = 3;
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL)    nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      nNodes = 5;
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        nNodes = 6;
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   nNodes = 8;
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNodes = 0; iNodes < nNodes; iNodes++) {
      PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
      for (iDim = 0; iDim < nDim; iDim++) {
        CoordCorners[iNodes][iDim] = geometry->node[PointCorners[iNodes]]->GetCoord(iDim);
      }
    }
    
    /*--- We set the element stiffness matrix ---*/
    
    if (nDim == 2) numerics->SetFEA_StiffMatrix2D(StiffMatrix_Elem, CoordCorners, nNodes, form2d);
    if (nDim == 3) numerics->SetFEA_StiffMatrix3D(StiffMatrix_Elem, CoordCorners, nNodes);
    
    /*--- Initialization of the auxiliar matrix ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        StiffMatrix_Node[iVar][jVar] = 0.0;
    
    /*--- Transform the stiffness matrix into the
     contributions for the individual nodes relative to each other. ---*/
    
    for (iVar = 0; iVar < nNodes; iVar++) {
      for (jVar = 0; jVar < nNodes; jVar++) {
        for (iDim = 0; iDim < nVar; iDim++) {
          for (jDim = 0; jDim < nVar; jDim++) {
            StiffMatrix_Node[iDim][jDim] = StiffMatrix_Elem[(iVar*nDim)+iDim][(jVar*nDim)+jDim];
          }
        }
        StiffMatrixSpace.AddBlock(PointCorners[iVar], PointCorners[jVar], StiffMatrix_Node);
      }
    }
    
  }
  
}

void CFEASolver::Compute_StiffMassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
  
  unsigned short iVar, jVar, nNodes = 0, iNodes, iDim, jDim, form2d;
  unsigned long iElem, PointCorners[8];
  su2double CoordCorners[8][3];
  
  form2d=config->GetElas2D_Formulation();
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     nNodes = 3;
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL)    nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      nNodes = 5;
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        nNodes = 6;
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   nNodes = 8;
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNodes = 0; iNodes < nNodes; iNodes++) {
      PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
      for (iDim = 0; iDim < nDim; iDim++) {
        CoordCorners[iNodes][iDim] = geometry->node[PointCorners[iNodes]]->GetCoord(iDim);
      }
    }
    
    /*--- We set the element stiffness matrix ---*/
    
    /*--- This solves the problem but... why? ---*/
    for (iVar = 0; iVar < nNodes*nDim; iVar++) {
      StiffMatrix_Elem[iVar] = new su2double [nNodes*nDim];
      for (jVar = 0; jVar < nNodes*nDim; jVar++) {
        StiffMatrix_Elem[iVar][jVar] = 0.0;
      }
    }
    
    if (nDim == 2) numerics->SetFEA_StiffMassMatrix2D(StiffMatrix_Elem, MassMatrix_Elem, CoordCorners, nNodes, form2d);
    if (nDim == 3) numerics->SetFEA_StiffMassMatrix3D(StiffMatrix_Elem, MassMatrix_Elem, CoordCorners, nNodes);
    
    /*--- Initialization of the auxiliar matrix ---*/
    
    for (iVar = 0; iVar < nVar; iVar++){
      for (jVar = 0; jVar < nVar; jVar++){
        StiffMatrix_Node[iVar][jVar] = 0.0;
        MassMatrix_Node[iVar][jVar] = 0.0;
      }
    }
    
    /*--- Transform the stiffness and mass matrices into the
     contributions for the individual nodes relative to each other. ---*/
    
    for (iVar = 0; iVar < nNodes; iVar++) {
      for (jVar = 0; jVar < nNodes; jVar++) {
        for (iDim = 0; iDim < nVar; iDim++) {
          for (jDim = 0; jDim < nVar; jDim++) {
            StiffMatrix_Node[iDim][jDim] = StiffMatrix_Elem[(iVar*nDim)+iDim][(jVar*nDim)+jDim];
            MassMatrix_Node[iDim][jDim] = MassMatrix_Elem[(iVar*nDim)+iDim][(jVar*nDim)+jDim];
            MassMatrix_Node_Int[iDim][jDim] = a_dt[0] * MassMatrix_Elem[(iVar*nDim)+iDim][(jVar*nDim)+jDim];
          }
        }
        MassMatrix.AddBlock(PointCorners[iVar], PointCorners[jVar], MassMatrix_Node);
        StiffMatrixTime.AddBlock(PointCorners[iVar], PointCorners[jVar], StiffMatrix_Node);
        StiffMatrixTime.AddBlock(PointCorners[iVar], PointCorners[jVar], MassMatrix_Node_Int);
      }
    }
    
  }
  
}

void CFEASolver::Compute_StiffMassDampMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
  
  cout << "Here we will compute the damping matrix." << endl;
}


void CFEASolver::SetSolution_time_n(CGeometry *geometry, CConfig *config) {
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  if (dynamic){
    for(unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      // The loop is over nPoints so the boundaries are also updated
      node[iPoint]->SetSolution_time_n();
      node[iPoint]->SetSolution_Vel_time_n();
      node[iPoint]->SetSolution_Accel_time_n();
    }
  }
  
}





void CFEASolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                 CConfig *config, unsigned short iMesh) {
  
  /*--- Compute body forces load vector ---*/
  
  
  
  /*--- Compute initial stresses effect ---*/
  
  
}

void CFEASolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                  CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
}


/*--------------------------------------------------------------------------------------------------
 * Definition of new boundary conditions
 ---------------------------------------------------------------------------------------------------*/

void CFEASolver::BC_Clamped(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) {

	unsigned long iPoint, iVertex;
	unsigned short iVar, jVar;

	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

	su2double **mIdentity, **mZeros;  // Variables to delete blocks in the jacobian

	mIdentity = new su2double *[nDim]; // Number of rows, allocate memory for each
	for(int iMat=0; iMat<nDim; iMat++) // i < Number of rows
		mIdentity[iMat] = new su2double[nDim]; // Number of columns, allocate memory for each

	mZeros = new su2double *[nDim]; // Number of rows, allocate memory for each
	for(int iMat=0; iMat<nDim; iMat++) // i < Number of rows
		mZeros[iMat] = new su2double[nDim]; // Number of columns, allocate memory for each

	// Initialise matrices

	for(int iMat=0; iMat<nDim; iMat++){
		for(int jMat=0; jMat<nDim; jMat++){
			mZeros[iMat][jMat]=0.0;
//			if (iMat==jMat) mIdentity[iMat][jMat]=config->GetElasticyMod();
			if (iMat==jMat) mIdentity[iMat][jMat]=1.0;
			else mIdentity[iMat][jMat]=0;
		}
	}


	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		/*--- Get node index ---*/

		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		if (nDim == 2) {
			Solution[0] = 0.0;  Solution[1] = 0.0;
			Residual[0] = 0.0;  Residual[1] = 0.0;
		}
		else {
			Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
			Residual[0] = 0.0;  Residual[1] = 0.0;  Residual[2] = 0.0;
		}

		node[iPoint]->SetSolution(Solution);

		if (dynamic){
			node[iPoint]->SetSolution_Vel(Solution);
			node[iPoint]->SetSolution_Accel(Solution);
		}

		LinSysRes.SetBlock(iPoint, Residual);

		/*--- Set the boundary clamped condition ---*/

		/*--- If the problem is dynamic ---*/

		if(dynamic){

			/*--- Enforce that in the previous time step all nodes had 0 U, U', U'' ---*/

			node[iPoint]->SetSolution_time_n(Solution);
			node[iPoint]->SetSolution_Vel_time_n(Solution);
			node[iPoint]->SetSolution_Accel_time_n(Solution);

			/*--- Delete the rows for a particular node ---*/
			for (jVar = 0; jVar < nPoint; jVar++){
				if (iPoint==jVar) {
					StiffMatrixTime.SetBlock(iPoint,jVar,mIdentity);
					MassMatrix.SetBlock(iPoint,jVar,mIdentity);
				}
				else {
					StiffMatrixTime.SetBlock(iPoint,jVar,mZeros);
					MassMatrix.SetBlock(iPoint,jVar,mZeros);
				}
			}

			/*--- Delete the columns for a particular node ---*/
			for (iVar = 0; iVar < nPoint; iVar++){
				if (iVar==iPoint) {
					StiffMatrixTime.SetBlock(iVar,iPoint,mIdentity);
					MassMatrix.SetBlock(iVar,iPoint,mIdentity);
				}
				else {
					StiffMatrixTime.SetBlock(iVar,iPoint,mZeros);
					MassMatrix.SetBlock(iVar,iPoint,mZeros);
				}
			}

		}

		/*--- If the problem is static ---*/

		else{

			/*--- Delete the rows for a particular node ---*/
			for (jVar = 0; jVar < nPoint; jVar++){
				if (iPoint==jVar) {
					Jacobian.SetBlock(iPoint,jVar,mIdentity);
					StiffMatrixSpace.SetBlock(iPoint,jVar,mIdentity);
				}
				else {
					Jacobian.SetBlock(iPoint,jVar,mZeros);
					StiffMatrixSpace.SetBlock(iPoint,jVar,mZeros);
				}
			}

			/*--- Delete the columns for a particular node ---*/
			for (iVar = 0; iVar < nPoint; iVar++){
				if (iVar==iPoint) {
					Jacobian.SetBlock(iVar,iPoint,mIdentity);
					StiffMatrixSpace.SetBlock(iPoint,jVar,mIdentity);
				}
				else {
					Jacobian.SetBlock(iVar,iPoint,mZeros);
					StiffMatrixSpace.SetBlock(iPoint,jVar,mZeros);
				}
			}

		}

	}


}


void CFEASolver::BC_Clamped_Post(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) {

	unsigned long iPoint, iVertex;

	su2double **mIdentity, **mZeros;  // Variables to delete blocks in the jacobian


	mIdentity = new su2double *[nDim]; // Number of rows, allocate memory for each
	for(int iMat=0; iMat<nDim; iMat++) // i < Number of rows
		mIdentity[iMat] = new su2double[nDim]; // Number of columns, allocate memory for each

	mZeros = new su2double *[nDim]; // Number of rows, allocate memory for each
	for(int iMat=0; iMat<nDim; iMat++) // i < Number of rows
		mZeros[iMat] = new su2double[nDim]; // Number of columns, allocate memory for each


	// Initialise matrices

	for(int iMat=0; iMat<nDim; iMat++){
		for(int jMat=0; jMat<nDim; jMat++){
			mZeros[iMat][jMat]=0;
			if (iMat==jMat) mIdentity[iMat][jMat]=1;
			else mIdentity[iMat][jMat]=0;
		}
	}


	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
	iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

	if (nDim == 2) {
		Solution[0] = 0.0;  Solution[1] = 0.0;
		Residual[0] = 0.0;  Residual[1] = 0.0;
	}
	else {
		Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
		Residual[0] = 0.0;  Residual[1] = 0.0;  Residual[2] = 0.0;
	}

	node[iPoint]->SetSolution(Solution);
	node[iPoint]->SetSolution_Old(Solution);


    /*--- Re-set the displacement condition ---*/
    LinSysRes.SetBlock(iPoint, Residual);
    
  }
  
}


void CFEASolver::BC_Normal_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) {
  unsigned long iPoint, iVertex, total_index;
  unsigned short iVar, iDim;
  su2double *Normal, Area, UnitaryNormal[3] = {0.0,0.0,0.0};
  
  su2double TotalDispl = config->GetDispl_Value(config->GetMarker_All_TagBound(val_marker));
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
    
    /*--- Compute area, and unitary normal ---*/
    
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    for (iDim = 0; iDim < nDim; iDim++) UnitaryNormal[iDim] = Normal[iDim]/Area;
    
    if (config->GetUnsteady_Simulation() == STEADY) {
      if (nDim == 2) {
        Solution[0] = TotalDispl*UnitaryNormal[0];  Solution[1] = TotalDispl*UnitaryNormal[1];
        Residual[0] = TotalDispl*UnitaryNormal[0];  Residual[1] = TotalDispl*UnitaryNormal[1];
      }
      else {
        Solution[0] = TotalDispl*UnitaryNormal[0];  Solution[1] = TotalDispl*UnitaryNormal[1];  Solution[2] = TotalDispl*UnitaryNormal[2];
        Residual[0] = TotalDispl*UnitaryNormal[0];  Residual[1] = TotalDispl*UnitaryNormal[1];  Residual[2] = TotalDispl*UnitaryNormal[2];
      }
    }
    else {
      if (nDim == 2) {
        Solution[0] = TotalDispl*UnitaryNormal[0];  Solution[1] = TotalDispl*UnitaryNormal[1];
        Solution[2] = 0.0;                          Solution[3] = 0.0;
        Residual[0] = 0.0;                          Residual[1] = 0.0;
        Residual[2] = 0.0;                          Residual[3] = 0.0;
      }
      else {
        Solution[0] = TotalDispl*UnitaryNormal[0];  Solution[1] = TotalDispl*UnitaryNormal[1];  Solution[2] = TotalDispl*UnitaryNormal[2];
        Solution[3] = 0.0;                          Solution[4] = 0.0;                          Solution[5] = 0.0;
        Residual[0] = 0.0;                          Residual[1] = 0.0;                          Residual[2] = 0.0;
        Residual[3] = 0.0;                          Residual[4] = 0.0;                          Residual[5] = 0.0;
      }
    }
    
    node[iPoint]->SetSolution(Solution);
    node[iPoint]->SetSolution_Old(Solution);
    
    LinSysRes.SetBlock(iPoint, Residual);
    
    /*--- Set the dirichlet condition ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      Jacobian.DeleteValsRowi(total_index);
    }
    
  }
}

void CFEASolver::BC_Normal_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                unsigned short val_marker) {
  
  su2double a[3], b[3];
  unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
  su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
  su2double Normal_Elem[3] = {0.0, 0.0, 0.0};
  //su2double Length_Elem = 0.0, Area_Elem = 0.0;
  unsigned short iDim;
  
  su2double TotalLoad = config->GetLoad_Value(config->GetMarker_All_TagBound(val_marker));
  
  for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
    
    Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);                   Coord_0 = geometry->node[Point_0]->GetCoord();
    Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);                   Coord_1 = geometry->node[Point_1]->GetCoord();
    if (nDim == 3) { Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord(); }
    
    /*--- Compute area (3D), and length of the surfaces (2D) ---*/
    
    if (nDim == 2) {
      
      for (iDim = 0; iDim < nDim; iDim++) a[iDim] = Coord_0[iDim]-Coord_1[iDim];
      
      //Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);
      Normal_Elem[0] =   a[1];
      Normal_Elem[1] = -(a[0]);
      
    }
    
    if (nDim == 3) {
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = Coord_0[iDim]-Coord_2[iDim];
        b[iDim] = Coord_1[iDim]-Coord_2[iDim];
      }
      
      //Area_Elem = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
      
      Normal_Elem[0] = -(0.5*(a[1]*b[2]-a[2]*b[1]));
      Normal_Elem[1] = -(-0.5*(a[0]*b[2]-a[2]*b[0]));
      Normal_Elem[2] = -(0.5*(a[0]*b[1]-a[1]*b[0]));
      
    }
    
    
    if (nDim == 2) {
      Residual[0] = (1.0/2.0)*TotalLoad*Normal_Elem[0]; Residual[1] = (1.0/2.0)*TotalLoad*Normal_Elem[1];
      LinSysRes.AddBlock(Point_0, Residual);
      Residual[0] = (1.0/2.0)*TotalLoad*Normal_Elem[0]; Residual[1] = (1.0/2.0)*TotalLoad*Normal_Elem[1];
      LinSysRes.AddBlock(Point_1, Residual);
    }
    
    else {
      Residual[0] = (1.0/3.0)*TotalLoad*Normal_Elem[0]; Residual[1] = (1.0/3.0)*TotalLoad*Normal_Elem[1]; Residual[2] = (1.0/3.0)*TotalLoad*Normal_Elem[2];
      LinSysRes.AddBlock(Point_0, Residual);
      
      Residual[0] = (1.0/3.0)*TotalLoad*Normal_Elem[0]; Residual[1] = (1.0/3.0)*TotalLoad*Normal_Elem[1]; Residual[2] = (1.0/3.0)*TotalLoad*Normal_Elem[2];
      LinSysRes.AddBlock(Point_1, Residual);
      
      Residual[0] = (1.0/3.0)*TotalLoad*Normal_Elem[0]; Residual[1] = (1.0/3.0)*TotalLoad*Normal_Elem[1]; Residual[2] = (1.0/3.0)*TotalLoad*Normal_Elem[2];
      LinSysRes.AddBlock(Point_2, Residual);
    }
    
  }
  
}


void CFEASolver::BC_Dir_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                             unsigned short val_marker) {
  
  su2double a[3] = {0.0, 0.0, 0.0}, b[3] = {0.0, 0.0, 0.0};
  su2double AC[3] = {0.0, 0.0, 0.0}, BD[3] = {0.0, 0.0, 0.0};
  unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3=0;
  su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
  su2double Length_Elem = 0.0, Area_Elem = 0.0;
//  su2double Normal_Elem[3] = {0.0, 0.0, 0.0};
  unsigned short iDim;
  
  su2double LoadDirVal = config->GetLoad_Dir_Value(config->GetMarker_All_TagBound(val_marker));
  su2double LoadDirMult = config->GetLoad_Dir_Multiplier(config->GetMarker_All_TagBound(val_marker));
  su2double *Load_Dir_Local= config->GetLoad_Dir(config->GetMarker_All_TagBound(val_marker));
  
  su2double TotalLoad;
  
  bool Gradual_Load = config->GetGradual_Load();
  su2double CurrentTime=config->GetCurrent_DynTime();
  su2double ModAmpl, NonModAmpl;
  
  bool Ramp_Load = config->GetRamp_Load();
  su2double Ramp_Time = config->GetRamp_Time();
  
  if (Ramp_Load){
    ModAmpl=LoadDirVal*LoadDirMult*CurrentTime/Ramp_Time;
    NonModAmpl=LoadDirVal*LoadDirMult;
    TotalLoad=min(ModAmpl,NonModAmpl);
  }
  else if (Gradual_Load){
    ModAmpl=2*((1/(1+exp(-1*CurrentTime)))-0.5);
    TotalLoad=ModAmpl*LoadDirVal*LoadDirMult;
  }
  else{
    TotalLoad=LoadDirVal*LoadDirMult;
  }
  
  /*--- Compute the norm of the vector that was passed in the config file ---*/
  su2double Norm = 0.0;
  if (nDim==2) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]);
  if (nDim==3) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]+Load_Dir_Local[2]*Load_Dir_Local[2]);
  
  for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
    
    Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);     Coord_0 = geometry->node[Point_0]->GetCoord();
    Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);     Coord_1 = geometry->node[Point_1]->GetCoord();
    
    /*--- Compute area (3D), and length of the surfaces (2D) ---*/
    
    if (nDim == 2) {
      
      for (iDim = 0; iDim < nDim; iDim++) a[iDim] = Coord_0[iDim]-Coord_1[iDim];
      
      Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);
//      Normal_Elem[0] =   a[1];
//      Normal_Elem[1] = -(a[0]);
      
    } else {
      
      Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);
      Coord_2 = geometry->node[Point_2]->GetCoord();
      
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE){
        
        for (iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = Coord_1[iDim]-Coord_0[iDim];
          b[iDim] = Coord_2[iDim]-Coord_0[iDim];
        }
        
        su2double Ni=0 , Nj=0, Nk=0;
        
        Ni=a[1]*b[2]-a[2]*b[1];
        Nj=-a[0]*b[2]+a[2]*b[0];
        Nk=a[0]*b[1]-a[1]*b[0];
        
        Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
        
        //Area_Elem = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
        
      } else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL){
        
        Point_3 = geometry->bound[val_marker][iElem]->GetNode(3);
        Coord_3 = geometry->node[Point_3]->GetCoord();
        
        for (iDim = 0; iDim < nDim; iDim++) {
          AC[iDim] = Coord_2[iDim]-Coord_0[iDim];
          BD[iDim] = Coord_3[iDim]-Coord_1[iDim];
        }
        
        su2double Ni=0 , Nj=0, Nk=0;
        
        Ni=AC[1]*BD[2]-AC[2]*BD[1];
        Nj=-AC[0]*BD[2]+AC[2]*BD[0];
        Nk=AC[0]*BD[1]-AC[1]*BD[0];
        
        Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
        
      }
    }
    
    if (nDim == 2) {
      
      Residual[0] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
      Residual[1] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
      
      LinSysRes.AddBlock(Point_0, Residual);
      LinSysRes.AddBlock(Point_1, Residual);
      
    }
    
    else {
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE){
        
        Residual[0] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
        Residual[2] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;
        
        LinSysRes.AddBlock(Point_0, Residual);
        LinSysRes.AddBlock(Point_1, Residual);
        LinSysRes.AddBlock(Point_2, Residual);
      }
      else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL){
        
        Residual[0] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
        Residual[2] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;
        
        LinSysRes.AddBlock(Point_0, Residual);
        LinSysRes.AddBlock(Point_1, Residual);
        LinSysRes.AddBlock(Point_2, Residual);
        LinSysRes.AddBlock(Point_3, Residual);
        
      }
      
    }
    
  }
  
}

void CFEASolver::BC_Sine_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                              unsigned short val_marker) {
  
  su2double a[3] = {0.0, 0.0, 0.0}, b[3] = {0.0, 0.0, 0.0};
  su2double AC[3] = {0.0, 0.0, 0.0}, BD[3] = {0.0, 0.0, 0.0};
  unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3=0;
  su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
  su2double Length_Elem = 0.0, Area_Elem = 0.0;
//  su2double Normal_Elem[3] = {0.0, 0.0, 0.0};
  unsigned short iDim;
  
  su2double LoadAmplitude = config->GetLoad_Sine_Amplitude(config->GetMarker_All_TagBound(val_marker));
  su2double LoadFrequency = config->GetLoad_Sine_Frequency(config->GetMarker_All_TagBound(val_marker));
  su2double *Load_Dir_Local= config->GetLoad_Sine_Dir(config->GetMarker_All_TagBound(val_marker));
  
  su2double CurrentTime=config->GetCurrent_DynTime();
  
  su2double TotalLoad = LoadAmplitude*sin(2*PI_NUMBER*LoadFrequency*CurrentTime);
  
  /*--- Compute the norm of the vector that was passed in the config file ---*/
  su2double Norm = 0.0;
  if (nDim==2) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]);
  if (nDim==3) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]+Load_Dir_Local[2]*Load_Dir_Local[2]);
  
  for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
    
    Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);     Coord_0 = geometry->node[Point_0]->GetCoord();
    Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);     Coord_1 = geometry->node[Point_1]->GetCoord();
    
    /*--- Compute area (3D), and length of the surfaces (2D) ---*/
    
    if (nDim == 2) {
      
      for (iDim = 0; iDim < nDim; iDim++) a[iDim] = Coord_0[iDim]-Coord_1[iDim];
      
      Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);
//      Normal_Elem[0] =   a[1];
//      Normal_Elem[1] = -(a[0]);
      
    } else { // 3D
      
      Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);
      Coord_2 = geometry->node[Point_2]->GetCoord();
      
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE){
        
        for (iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = Coord_1[iDim]-Coord_0[iDim];
          b[iDim] = Coord_2[iDim]-Coord_0[iDim];
        }
        
        su2double Ni=0 , Nj=0, Nk=0;
        
        Ni=a[1]*b[2]-a[2]*b[1];
        Nj=-a[0]*b[2]+a[2]*b[0];
        Nk=a[0]*b[1]-a[1]*b[0];
        
        Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
        
        
        //Area_Elem = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
        
      } else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL){
        
        Point_3 = geometry->bound[val_marker][iElem]->GetNode(3);
        Coord_3 = geometry->node[Point_3]->GetCoord();
        
        for (iDim = 0; iDim < nDim; iDim++) {
          AC[iDim] = Coord_2[iDim]-Coord_0[iDim];
          BD[iDim] = Coord_3[iDim]-Coord_1[iDim];
        }
        
        su2double Ni=0 , Nj=0, Nk=0;
        
        Ni=AC[1]*BD[2]-AC[2]*BD[1];
        Nj=-AC[0]*BD[2]+AC[2]*BD[0];
        Nk=AC[0]*BD[1]-AC[1]*BD[0];
        
        Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
        
      }
    }
    
    if (nDim == 2) {
      
      Residual[0] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
      Residual[1] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
      
      LinSysRes.AddBlock(Point_0, Residual);
      LinSysRes.AddBlock(Point_1, Residual);
      
    }
    
    else {
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE){
        
        Residual[0] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
        Residual[2] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;
        
        LinSysRes.AddBlock(Point_0, Residual);
        LinSysRes.AddBlock(Point_1, Residual);
        LinSysRes.AddBlock(Point_2, Residual);
        
        
      }
      else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL){
        
        Residual[0] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
        Residual[2] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;
        
        LinSysRes.AddBlock(Point_0, Residual);
        LinSysRes.AddBlock(Point_1, Residual);
        LinSysRes.AddBlock(Point_2, Residual);
        LinSysRes.AddBlock(Point_3, Residual);
        
      }
      
    }
    
  }
  
}

void CFEASolver::BC_Pressure(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                             unsigned short val_marker) { }

void CFEASolver::BC_Flow_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                              unsigned short val_marker) { }


void CFEASolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics_container, unsigned short iMesh) {

  unsigned long iPoint, iElem;
  su2double **Stress, VonMises_Stress=0.0, MaxVonMises_Stress = 0.0;
  su2double Sxx=0.0,Syy=0.0,Szz=0.0,Sxy=0.0,Sxz=0.0,Syz=0.0,S1,S2;
  
  unsigned long PointCorners[8];
  unsigned short nNodes=0, iNodes, iDim, jDim, form2d;
  su2double CoordCorners[8][3];
//  su2double CoordGauss[8][3];
  
  /*--- Container of the shape functions ---*/
  CNumerics *numerics;
  numerics=numerics_container[VISC_TERM];
  
  /*--- Enforcement of displacement boundary conditions ---*/
  unsigned short MainSolver = config->GetContainerPosition(RUNTIME_FEA_SYS);
  unsigned int iMarker;
  
  form2d=config->GetElas2D_Formulation();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case CLAMPED_BOUNDARY:
        solver_container[MainSolver]->BC_Clamped(geometry, solver_container, numerics_container[CONV_BOUND_TERM], config, iMarker);
        break;
    }
  }
  
  /* --- Initialize the stress and the number of elements connected to each node ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->Initialize_Connectivity();
    for (iDim = 0; iDim < nDim; iDim++){
      for (jDim = 0; jDim < nDim; jDim++){
        node[iPoint]->SetStress(iDim, jDim, 0);
      }
    }
  }
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE){     nNodes = 3;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL){    nNodes = 4;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON){  nNodes = 4;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID){      nNodes = 5;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM){        nNodes = 6;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON){   nNodes = 8;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNodes = 0; iNodes < nNodes; iNodes++) {
      
      /*--- Get the index of the nodes and saves it in PointCorners ---*/
      PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
      
      for (iDim = 0; iDim < nDim; iDim++) {
        
        /*--- Get the coordinates of the nodes and saves it in CoordCorners ---*/
        CoordCorners[iNodes][iDim] = geometry->node[PointCorners[iNodes]]->GetCoord(iDim);
        
        /*--- Initialization of the gauss coordinate matrix ---*/
//        CoordGauss[iNodes][iDim] = 0.0;
      }
    }
    
    /*----------------------------------------------------------------------------------*/
    /*--- We obtain the stresses in the element, from the finite element formulation ---*/
    /*----------------------------------------------------------------------------------*/
    
    /*--- For a 2D element ---*/
    
    if (nDim == 2) {
      
      su2double StressNodal[8][3], DispElement[8];
      
      /*--- Set the element displacements vector, from the global solution ---*/
      
      for (iNodes = 0; iNodes < nNodes; iNodes++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          DispElement[nDim*iNodes+iDim]=node[PointCorners[iNodes]]->GetSolution(iDim);
        }
      }
      
      /*--- Obtain the stresses in the nodes ---*/
      numerics->GetFEA_StressNodal2D(StressNodal, DispElement, CoordCorners, nNodes, form2d);
      
      /*--- Add the value of the element stress extrapolated to the node to the value stored on the nodes  ---*/
      /*--- At the same point, add a counter to take into account on how many elements connect to the node ---*/
      
      for (iPoint = 0; iPoint < nNodes; iPoint++) {
        
        node[PointCorners[iPoint]]->AddStress(0, 0, StressNodal[iPoint][0]);
        node[PointCorners[iPoint]]->AddStress(1, 1, StressNodal[iPoint][1]);
        node[PointCorners[iPoint]]->AddStress(0, 1, StressNodal[iPoint][2]);
        node[PointCorners[iPoint]]->AddStress(1, 0, StressNodal[iPoint][2]);
        
        node[PointCorners[iPoint]]->Upgrade_Connectivity();
        
      }
      
    }
    
    /*--- For a 3D element ---*/
    
    if (nDim == 3) {
      
      su2double StressNodal[8][6], DispElement[24];
      
      /*--- Set the element displacements vector, from the global solution ---*/
      
      for (iNodes = 0; iNodes < nNodes; iNodes++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          DispElement[nDim*iNodes+iDim]=node[PointCorners[iNodes]]->GetSolution(iDim);
        }
      }
      
      /*--- Obtain the stresses in the gaussian points ---*/
      numerics->GetFEA_StressNodal3D(StressNodal, DispElement, CoordCorners, nNodes);
      
      /*--- Add the value of the element stress extrapolated to the node to the value stored on the nodes  ---*/
      /*--- At the same point, add a counter to take into account on how many elements connect to the node ---*/
      
      for (iPoint = 0; iPoint < nNodes; iPoint++) {
        
        node[PointCorners[iPoint]]->AddStress(0, 0, StressNodal[iPoint][0]);
        node[PointCorners[iPoint]]->AddStress(1, 1, StressNodal[iPoint][1]);
        node[PointCorners[iPoint]]->AddStress(2, 2, StressNodal[iPoint][2]);
        
        node[PointCorners[iPoint]]->AddStress(0, 1, StressNodal[iPoint][3]);
        node[PointCorners[iPoint]]->AddStress(1, 0, StressNodal[iPoint][3]);
        
        node[PointCorners[iPoint]]->AddStress(0, 2, StressNodal[iPoint][4]);
        node[PointCorners[iPoint]]->AddStress(2, 0, StressNodal[iPoint][4]);
        
        node[PointCorners[iPoint]]->AddStress(1, 2, StressNodal[iPoint][5]);
        node[PointCorners[iPoint]]->AddStress(2, 1, StressNodal[iPoint][5]);
        
        node[PointCorners[iPoint]]->Upgrade_Connectivity();
        
      }
      
    }
    
  }
  
  
  /* --- Variable to store the number of elements connected to each node ---*/
  
  su2double nElPerNode=0;
  
  /* --- For the number of nodes in the mesh ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /* --- Get the stresses, added up from all the elements that connect to the node ---*/
    
    Stress     = node[iPoint]->GetStress();
    nElPerNode = node[iPoint]->Get_Connectivity();
    
    /* --- Compute the stress averaged from all the elements connecting to the node and the Von Mises stress ---*/
    
    if (geometry->GetnDim() == 2) {
      
      Sxx=Stress[0][0]/nElPerNode;
      Syy=Stress[1][1]/nElPerNode;
      Sxy=Stress[0][1]/nElPerNode;
      
      S1=(Sxx+Syy)/2+sqrt(((Sxx-Syy)/2)*((Sxx-Syy)/2)+Sxy*Sxy);
      S2=(Sxx+Syy)/2-sqrt(((Sxx-Syy)/2)*((Sxx-Syy)/2)+Sxy*Sxy);
      
      VonMises_Stress = sqrt(S1*S1+S2*S2-2*S1*S2);
      
    }
    else if (geometry->GetnDim() == 3) {
      
      Sxx = Stress[0][0]/nElPerNode;
      Syy = Stress[1][1]/nElPerNode;
      Szz = Stress[2][2]/nElPerNode;
      
      Sxy = Stress[0][1]/nElPerNode;
      Sxz = Stress[0][2]/nElPerNode;
      Syz = Stress[1][2]/nElPerNode;
      
      VonMises_Stress = sqrt(0.5*(    pow(Sxx - Syy, 2.0)
                                  + pow(Syy - Szz, 2.0)
                                  + pow(Szz - Sxx, 2.0)
                                  + 6.0*(Sxy*Sxy+Sxz*Sxz+Syz*Syz)
                                  ));
      
    }
    
    node[iPoint]->SetVonMises_Stress(VonMises_Stress);
    
    /*--- Compute the maximum value of the Von Mises Stress ---*/
    
    MaxVonMises_Stress = max(MaxVonMises_Stress, VonMises_Stress);
    
    /*--- Set the new value of the stress, averaged from the number of elements ---*/
    
    node[iPoint]->SetStress(0, 0, Sxx);
    node[iPoint]->SetStress(1, 1, Syy);
    node[iPoint]->SetStress(0, 1, Sxy);
    node[iPoint]->SetStress(1, 0, Sxy);
    
    if (geometry->GetnDim() == 3) {
      node[iPoint]->SetStress(2, 2, Szz);
      node[iPoint]->SetStress(0, 2, Sxz);
      node[iPoint]->SetStress(2, 0, Sxz);
      node[iPoint]->SetStress(1, 2, Syz);
      node[iPoint]->SetStress(2, 1, Syz);
    }
    
  }
  
#ifdef HAVE_MPI
  
  /*--- Compute MaxVonMises_Stress using all the nodes ---*/
  
  su2double MyMaxVonMises_Stress = MaxVonMises_Stress; MaxVonMises_Stress = 0.0;
  SU2_MPI::Allreduce(&MyMaxVonMises_Stress, &MaxVonMises_Stress, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
#endif
  
  /*--- Set the value of the MaxVonMises_Stress as the CFEA coeffient ---*/
  
  Total_CFEA = MaxVonMises_Stress;
  
}

void CFEASolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep,
                                      unsigned short iMesh, unsigned short RunTime_EqSystem) {  }

void CFEASolver::ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  
  unsigned short iVar;
  unsigned long iPoint, total_index, IterLinSol;
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  unsigned long ExtIter = config->GetExtIter();
  
  su2double *PointTimeRes = NULL;
  
  unsigned short check=0;
  
  if ((dynamic) && (ExtIter == 0) && (check==0)){
    
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
    
    
    /*--- Solve the linear dynamic system ---*/
    
    CSysSolve femSystem;
    IterLinSol = femSystem.Solve(MassMatrix, LinSysRes, LinSysSol, geometry, config);
    SetIterLinSolver(IterLinSol);
    
    /*--- Update solution and (if dynamic) advance velocity and acceleration vectors in time ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        /*--- Acceleration for t=0 ---*/
        node[iPoint]->SetSolution(iVar, 0.0);
        node[iPoint]->SetSolution_Vel(iVar, 0.0);
        node[iPoint]->SetSolution_Accel(iVar, LinSysSol[iPoint*nVar+iVar]);
        
      }
      
    }
    
  }
  else{
    
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
    
    /*--- If dynamic analysis ---*/
    
    if (dynamic){
      
      /*--- Get mass term ---*/
      
      /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        for (iVar = 0; iVar < nVar; iVar++){
          
          Residual[iVar] = a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)+		//a0*U(t)
          a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)+	//a2*U'(t)
          a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);	//a3*U''(t)
          
        }
        
        TimeRes_Aux.SetBlock(iPoint, Residual);
        
      }
      
      /*--- Once computed, compute M*TimeRes_Aux ---*/
      
      MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);
      
      /*--- Add the components of M*TimeRes_Aux to the residual R(t+dt) ---*/
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        PointTimeRes = TimeRes.GetBlock(iPoint);
        
        LinSysRes.AddBlock(iPoint, PointTimeRes);
        
      }
      
      //			su2double *check;
      //			for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
      //			check = LinSysRes.GetBlock(iPoint);	// This avoids the problem in the corner, but...
      //			cout << check[0] << "\t" << check[1] << endl;
      //			}
      
      /*--- Solve the linear dynamic system ---*/
      
      CSysSolve femSystem;
      IterLinSol = femSystem.Solve(StiffMatrixTime, LinSysRes, LinSysSol, geometry, config);
      SetIterLinSolver(IterLinSol);
    }
    else {
      
      /*--- Solve the linear static system ---*/
      
      CSysSolve femSystem;
      IterLinSol = femSystem.Solve(StiffMatrixSpace, LinSysRes, LinSysSol, geometry, config);
      SetIterLinSolver(IterLinSol);
    }
    
    
    
    /*--- Update solution and (if dynamic) advance velocity and acceleration vectors in time ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        /*--- Displacements component of the solution ---*/
        
        node[iPoint]->SetSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
        
      }
      
      if (dynamic){
        
        for (iVar = 0; iVar < nVar; iVar++) {
          
          /*--- Acceleration component of the solution ---*/
          /*--- U''(t+dt) = a0*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/
          
          Solution[iVar]=a_dt[0]*(node[iPoint]->GetSolution(iVar) -
                                  node[iPoint]->GetSolution_time_n(iVar)) -
          a_dt[2]* node[iPoint]->GetSolution_Vel_time_n(iVar) -
          a_dt[3]* node[iPoint]->GetSolution_Accel_time_n(iVar);
          
        }
        
        /*--- Set the acceleration in the node structure ---*/
        
        node[iPoint]->SetSolution_Accel(Solution);
        
        for (iVar = 0; iVar < nVar; iVar++) {
          
          /*--- Velocity component of the solution ---*/
          /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/
          
          Solution[iVar]=node[iPoint]->GetSolution_Vel_time_n(iVar)+
          a_dt[6]* node[iPoint]->GetSolution_Accel_time_n(iVar) +
          a_dt[7]* node[iPoint]->GetSolution_Accel(iVar);
          
        }
        
        /*--- Set the velocity in the node structure ---*/
        
        node[iPoint]->SetSolution_Vel(Solution);
        
      }
      
    }
    
    /*--- MPI solution ---*/
    
    Set_MPI_Solution(geometry, config);
    
    /*---  Compute the residual Ax-f ---*/
    
    if (dynamic){
      
      StiffMatrixTime.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);
      
    }
    else {
      
      StiffMatrixSpace.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);
      
    }
    
    /*--- Compute the reactions ---*/
    
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
  
}

void CFEASolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index, IterLinSol;
  
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
  IterLinSol = system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  SetIterLinSolver(IterLinSol);
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*---  Compute the residual Ax-f ---*/
  
  Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);
  
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

void CFEASolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  unsigned long iPoint;
  
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Problem dimension and physical time step ---*/
  unsigned short nDim = geometry[MESH_0]->GetnDim();
  
  for (iPoint = 0; iPoint < geometry[MESH_0]->GetnPoint(); iPoint++) {
    
    /*--- Set initial boundary condition at the first iteration ---*/
    if ((ExtIter == 0) && (!restart)) {
      
      if (config->GetUnsteady_Simulation() == STEADY) {
        if (nDim == 2) { Solution[0] = 0.0;	Solution[1] = 0.0; }
        else { Solution[0] = 0.0;	Solution[1] = 0.0;	Solution[2] = 0.0; }
      }
      else {
        if (nDim == 2) { Solution[0] = 0.0;	Solution[1] = 0.0; Solution[2] = 0.0;	Solution[3] = 0.0; }
        else { Solution[0] = 0.0;	Solution[1] = 0.0;	Solution[2] = 0.0; Solution[3] = 0.0;	Solution[4] = 0.0;	Solution[5] = 0.0; }
      }
      
      node[iPoint]->SetSolution(Solution);
      
      if (dual_time) {
        node[iPoint]->Set_Solution_time_n();
        node[iPoint]->Set_Solution_time_n1();
      }
      
    }
  }
}

void CFEASolver::GetSurface_Pressure(CGeometry *geometry, CConfig *config) {
  
  unsigned short iMarker, icommas, iDim;
  unsigned long iVertex, iPoint, iExtIter;
  su2double Pressure = 0.0, Dist, Coord[3];
  string text_line;
  string::size_type position;
  ifstream Surface_file;
  char buffer[50], cstr[200];
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  /*--- Reset the value of the Flow_Pressure ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    node[iPoint]->SetFlow_Pressure(0.0);
  
  for (iExtIter = 0; iExtIter < config->GetnExtIter(); iExtIter++) {
    
    /*--- Prepare to read surface sensitivity files (CSV) ---*/
    
    string surfadj_filename = config->GetSurfFlowCoeff_FileName();
    
    /*--- Remove the domain number from the surface csv filename ---*/
    
    if (size > SINGLE_NODE) {
      if ((rank+1 >= 0) && (rank+1 < 10)) surfadj_filename.erase (surfadj_filename.end()-2, surfadj_filename.end());
      if ((rank+1 >= 10) && (rank+1 < 100)) surfadj_filename.erase (surfadj_filename.end()-3, surfadj_filename.end());
      if ((rank+1 >= 100) && (rank+1 < 1000)) surfadj_filename.erase (surfadj_filename.end()-4, surfadj_filename.end());
      if ((rank+1 >= 1000) && (rank+1 < 10000)) surfadj_filename.erase (surfadj_filename.end()-5, surfadj_filename.end());
    }
    strcpy (cstr, surfadj_filename.c_str());
    
    /*--- Write file name with extension if unsteady or steady ---*/
    
    if ((config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) ||
        (config->GetUnsteady_Simulation() == TIME_SPECTRAL)) {
      if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.csv",  SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.csv",   SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.csv",    SU2_TYPE::Int(iExtIter));
      if  (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iExtIter));
    }
    else SPRINTF (buffer, ".csv");
    
    strcat (cstr, buffer);
    
    /*--- Open the surface file ---*/
    
    Surface_file.open(cstr, ios::in);
    
    getline(Surface_file, text_line);
    
    /*--- Res the surface file ---*/
    
    while (getline(Surface_file, text_line)) {
      
      /*--- Remove commas from the surface file ---*/
      
      for (icommas = 0; icommas < 100; icommas++) {
        position = text_line.find( ",", 0 );
        if (position!=string::npos) text_line.erase (position, 1);
      }
      
      /*--- Read the file ---*/
      istringstream point_line(text_line);
      if (nDim == 2) { point_line >> iPoint >> Coord[0] >> Coord[1] >> Pressure; }
      if (nDim == 3) { point_line >> iPoint >> Coord[0] >> Coord[1] >> Coord[2] >> Pressure; }
      
      /*--- Compute the distance from the surface to the points in the .csv files ---*/
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == PRESSURE_BOUNDARY) {
          for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            
            /*--- Compute the distance between the point and the grid points ---*/
            Dist = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Dist += (Coord[iDim]-geometry->node[iPoint]->GetCoord(iDim))*(Coord[iDim]-geometry->node[iPoint]->GetCoord(iDim));
            Dist = sqrt(Dist);
            
            /*--- Check the distance and set the pressure ---*/
            if (Dist < 1E-10) { node[iPoint]->SetFlow_Pressure(Pressure); }
            
          }
        }
      }
      
    }
    
    Surface_file.close();
    
  }
  
}

void CFEASolver::SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry, CGeometry **flow_geometry,
                             CConfig *fea_config, CConfig *flow_config, CNumerics *fea_numerics) {

	unsigned short nVertexFEA, nVertexFlow, iVertex, nMarkerFSIint, iDim, jDim;
	unsigned short markFEA, markFlow, iPoint, iMarkerFSIint;
	unsigned short nMarkerFEA, nMarkerFlow, iMarkerFEA, iMarkerFlow;
	unsigned long *nodeVertex, *donorVertex;
	su2double *nodePress, *nodeShearStress, **normalsVertex, **normalsVertex_Unit, **tn_f, *tn_e;
	su2double **tractionPrev;
	su2double Pdiff, maxPdiff=0.0;
	su2double FSDensity, *FSVelocity, FSVelocity2, factorForces;
	su2double Viscosity_Ref, Velocity_Ref, Density_Ref, Pressure_Ref;

	su2double *Velocity_ND, Density_ND, *Velocity_Real, Density_Real, Velocity2_Real, Velocity2_ND;

	bool compressible       = (flow_config->GetKind_Regime() == COMPRESSIBLE);
	bool incompressible     = (flow_config->GetKind_Regime() == INCOMPRESSIBLE);
//	bool freesurface        = (flow_config->GetKind_Regime() == FREESURFACE);

	bool viscous_flow        = ((flow_config->GetKind_Solver() == NAVIER_STOKES) ||
			(flow_config->GetKind_Solver() == RANS) );

	su2double Pinf;


	su2double ModAmpl;
	su2double CurrentTime=fea_config->GetCurrent_DynTime();
	su2double Static_Time=fea_config->GetStatic_Time();

  bool Ramp_Load = fea_config->GetRamp_Load();
	su2double Ramp_Time = fea_config->GetRamp_Time();

	if (CurrentTime <= Static_Time){
		ModAmpl=0.0;
	}
	else if((CurrentTime > Static_Time) &&
			(CurrentTime <= (Static_Time + Ramp_Time)) &&
			(Ramp_Load)){
		ModAmpl=(CurrentTime-Static_Time)/Ramp_Time;
		ModAmpl=max(ModAmpl,0.0);
		ModAmpl=min(ModAmpl,1.0);
	}
	else{
		ModAmpl=1.0;
	}

	/*--- Number of markers in the FSI interface ---*/
	nMarkerFSIint = (fea_config->GetMarker_n_FSIinterface())/2;

	/*--- Initialization of vectors of residuals ---*/
	/*--- WATCH OUT! This Shouldn't be here I think... For the dead load */

	for (iPoint = 0; iPoint < fea_geometry[MESH_0]->GetnPoint(); iPoint ++) {
		LinSysRes.SetBlock_Zero(iPoint);
	}

  	/*--- Redimensionalize the pressure ---*/

    Velocity_Real = flow_config->GetVelocity_FreeStream();
    Density_Real = flow_config->GetDensity_FreeStream();

    Velocity_ND = flow_config->GetVelocity_FreeStreamND();
    Density_ND = flow_config->GetDensity_FreeStreamND();


	Velocity2_Real = 0.0;
	Velocity2_ND = 0.0;
    for (iDim = 0; iDim < nDim; iDim++){
    	Velocity2_Real += Velocity_Real[iDim]*Velocity_Real[iDim];
    	Velocity2_ND += Velocity_ND[iDim]*Velocity_ND[iDim];
    }

  	Velocity_Ref  = flow_config->GetVelocity_Ref();
  	Viscosity_Ref = flow_config->GetViscosity_Ref();
  	Density_Ref   = flow_solution[MESH_0][FLOW_SOL]->GetDensity_Inf();
  	Pressure_Ref  = flow_config->GetPressure_Ref();

    factorForces = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);

	/*--- Loop over all the markers on the interface ---*/

	for (iMarkerFSIint=0; iMarkerFSIint < nMarkerFSIint; iMarkerFSIint++){

		nMarkerFEA=fea_geometry[MESH_0]->GetnMarker();
		nMarkerFlow=flow_geometry[MESH_0]->GetnMarker();

		/*--- Identification of the markers ---*/

		for (iMarkerFEA=0; iMarkerFEA < nMarkerFEA; iMarkerFEA++){
			if ( fea_config->GetMarker_All_FSIinterface(iMarkerFEA) == (iMarkerFSIint+1)){
				markFEA=iMarkerFEA;
			}
		}

		for (iMarkerFlow=0; iMarkerFlow < nMarkerFlow; iMarkerFlow++){
			if (flow_config->GetMarker_All_FSIinterface(iMarkerFlow) == (iMarkerFSIint+1)){
				markFlow=iMarkerFlow;
			}
		}


		nVertexFEA = fea_geometry[MESH_0]->GetnVertex(markFEA);
		nVertexFlow = flow_geometry[MESH_0]->GetnVertex(markFlow);

		nodePress = new su2double [nVertexFlow];
		nodeShearStress = new su2double [nVertexFlow];
		nodeVertex = new unsigned long [nVertexFlow];
		donorVertex = new unsigned long [nVertexFlow];

		tn_e = new su2double [nVar*nDim];

		tn_f = new su2double* [nVertexFlow];
		for (iVertex = 0; iVertex < nVertexFlow; iVertex++) {
			tn_f[iVertex] = new su2double[nDim];
		}

		normalsVertex = new su2double* [nVertexFlow];
		for (iVertex = 0; iVertex < nVertexFlow; iVertex++) {
			normalsVertex[iVertex] = new su2double[nDim];
		}

		normalsVertex_Unit = new su2double* [nVertexFlow];
		for (iVertex = 0; iVertex < nVertexFlow; iVertex++) {
			normalsVertex_Unit[iVertex] = new su2double[nDim];
		}

		su2double a[3], b[3];
		unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
		su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, **Grad_PrimVar;
		su2double Length_Elem = 0.0, Area_Elem = 0.0, Normal_Elem[3] = {0.0, 0.0, 0.0};
		su2double Viscosity = 0.0, Density = 0.0;
		su2double TauElem_0[3], Tau[3][3], TauElem_1[3], Tau_1[3][3];

		su2double Force_0, Force_1, Force_0_comp, Force_1_comp, div_vel, Delta;
		su2double Area;

		su2double CoordCorners[4][3], Fnodal[12], FnodalRelax[12];

		/*--- Loop over the nodes in the fluid mesh, calculate the tf vector (unitary) ---*/

		su2double Pn, Pnm1, check1, check2;

		/*--- Here, we are looping over the fluid, and we find the pointer to the structure (donorVertex) ---*/
		for (iVertex=0; iVertex < nVertexFlow; iVertex++){

			// Node from the flow mesh
			nodeVertex[iVertex]=flow_geometry[MESH_0]->vertex[markFlow][iVertex]->GetNode();

			// Normals at the vertex: these normals go inside the fluid domain.
			normalsVertex[iVertex]=flow_geometry[MESH_0]->vertex[markFlow][iVertex]->GetNormal();

			// Unit normals
	        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
	        	Area += normalsVertex[iVertex][iDim]*normalsVertex[iVertex][iDim];
	        }
	      	Area = sqrt(Area);

	        for (iDim = 0; iDim < nDim; iDim++) {
	          normalsVertex_Unit[iVertex][iDim] = normalsVertex[iVertex][iDim]/Area;
	        }

			// Corresponding node on the structural mesh
			donorVertex[iVertex]=flow_geometry[MESH_0]->vertex[markFlow][iVertex]->GetDonorPoint();

			// Retrieve the values of pressure, viscosity and density
			if (incompressible){

				Pn=flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetPressureInc();
				Pinf=flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();

				if (viscous_flow){

					Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetGradient_Primitive();
					Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetLaminarViscosityInc();
					Density = flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetDensityInc();

				}
			}
			else if (compressible){

				Pn=flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetPressure();
				Pinf=flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();

				if (viscous_flow){

					Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetGradient_Primitive();
					Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetLaminarViscosity();
					Density = flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetDensity();

				}
			}

			// Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
			for (iDim = 0; iDim < nDim; iDim++) {
				tn_f[iVertex][iDim] = -(Pn-Pinf)*normalsVertex[iVertex][iDim];
			}

			// Calculate tn in the fluid nodes for the viscous term

			if (viscous_flow){

				// Divergence of the velocity
				div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_PrimVar[iDim+1][iDim];
				if (incompressible) div_vel = 0.0;

				for (iDim = 0; iDim < nDim; iDim++) {

					for (jDim = 0 ; jDim < nDim; jDim++) {
						// Dirac delta
						Delta = 0.0; if (iDim == jDim) Delta = 1.0;

						// Viscous stress
						Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
								TWO3*Viscosity*div_vel*Delta;

						// Viscous component in the tn vector --> Units of force (non-dimensional).
						tn_f[iVertex][iDim] += Tau[iDim][jDim]*normalsVertex[iVertex][jDim];
					}
				}
			}

			// Rescale tn to SI units

			for (iDim = 0; iDim < nDim; iDim++) {
				tn_f[iVertex][iDim] = tn_f[iVertex][iDim]*factorForces;
			}

			// Apply time-dependent coefficient (static structure, ramp load, full load)

			for (iDim = 0; iDim < nDim; iDim++) {
				tn_f[iVertex][iDim] = tn_f[iVertex][iDim]*ModAmpl;
			}

			// This works only for matching meshes

			for (iDim=0; iDim < nDim; iDim++){
				Residual[iDim]=tn_f[iVertex][iDim];
			}

			LinSysRes.AddBlock(donorVertex[iVertex], Residual);

		}

	}

}


void CFEASolver::SetStruct_Displacement(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) {
  
  
  unsigned long iPoint, iDim;
  unsigned long nPoint, nDim;
  su2double *Coord, *VarCoord, *Displacement;
  
  
  nPoint = fea_geometry[MESH_0]->GetnPoint();
  nDim = fea_geometry[MESH_0]->GetnDim();
  
  VarCoord = new su2double [nDim];
  
  for (iPoint=0; iPoint < nPoint; iPoint++){
    
    Coord = fea_geometry[MESH_0]->node[iPoint]->GetCoord();
    
    Displacement = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
    
    for (iDim = 0; iDim < nDim; iDim++)
      VarCoord[iDim] = (Coord[iDim]+Displacement[iDim]);
    
    fea_geometry[MESH_0]->node[iPoint]->SetCoord(VarCoord);
    
  }
  
}


void CFEASolver::PredictStruct_Displacement(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) {
  
  unsigned short predOrder=fea_config->GetPredictorOrder();
  su2double Delta_t= fea_config->GetDelta_DynTime();
  unsigned long iPoint, iDim;
  unsigned long nPoint, nDim;
  su2double *solDisp, *solVel, *solVel_tn, *valPred, *checkPred;
  
  nPoint = fea_geometry[MESH_0]->GetnPoint();
  nDim = fea_geometry[MESH_0]->GetnDim();
  
  solDisp=new su2double [nDim];
  solVel=new su2double [nDim];
  solVel_tn=new su2double [nDim];
  valPred=new su2double [nDim];
  checkPred=new su2double [nDim];
  
  for (iPoint=0; iPoint<nPoint; iPoint++){
    if (predOrder==0) fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred();
    else if (predOrder==1) {
      
      solDisp = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
      solVel = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel();
      valPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
      
      for (iDim=0; iDim<nDim; iDim++){
        valPred[iDim] = solDisp[iDim] + Delta_t*solVel[iDim];
      }
      
      //			fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred(valPred);
      
      
    }
    else if (predOrder==2) {
      
      solDisp = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
      solVel = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel();
      solVel_tn = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n();
      valPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
      
      for (iDim=0; iDim<nDim; iDim++){
        valPred[iDim] = solDisp[iDim] + 0.5*Delta_t*(3*solVel[iDim]-solVel_tn[iDim]);
      }
      
      //			fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred(valPred);
      
    }
    else {
      cout<< "Higher order predictor not implemented. Solving with order 0." << endl;
      fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred();
    }
  }
  
  delete [] solDisp;
  delete [] solVel;
  delete [] solVel_tn;
  delete [] valPred;
  delete [] checkPred;
  
}

void CFEASolver::ComputeAitken_Coefficient(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution, unsigned long iFSIIter) {
  unsigned long iPoint, iDim;
  unsigned long nPoint, nDim;
  su2double *dispPred, *dispCalc, *dispPred_Old, *dispCalc_Old;
  su2double deltaU[3] = {0.0, 0.0, 0.0}, deltaU_p1[3] = {0.0, 0.0, 0.0};
  su2double delta_deltaU[3] = {0.0, 0.0, 0.0};
  su2double numAitk, denAitk, cocAitk, WAitken;
  su2double CurrentTime=fea_config->GetCurrent_DynTime();
  su2double Static_Time=fea_config->GetStatic_Time();
  su2double WAitkDyn_tn1, WAitkDyn_Max, WAitkDyn;

  nPoint = fea_geometry[MESH_0]->GetnPoint();
  nDim = fea_geometry[MESH_0]->GetnDim();

  WAitken=fea_config->GetAitkenStatRelax();

  dispPred	=new su2double [iDim];
  dispPred_Old=new su2double [iDim];
  dispCalc	=new su2double [iDim];
  dispCalc_Old=new su2double [iDim];

	numAitk = 0.0;
	denAitk = 0.0;

	ofstream historyFile_FSI;
	bool writeHistFSI = fea_config->GetWrite_Conv_FSI();
	if (writeHistFSI){
		char cstrFSI[200];
		string filenameHistFSI = fea_config->GetConv_FileName_FSI();
		strcpy (cstrFSI, filenameHistFSI.data());
		historyFile_FSI.open (cstrFSI, std::ios_base::app);
	}


	/*--- Only when there is movement, and a dynamic coefficient is requested, it makes sense to compute the Aitken's coefficient ---*/

	if (CurrentTime > Static_Time) {

		if (iFSIIter == 0){

			WAitkDyn_tn1 = GetWAitken_Dyn_tn1();
			WAitkDyn_Max = fea_config->GetAitkenDynMaxInit();

			WAitkDyn = min(WAitkDyn_tn1, WAitkDyn_Max);

			/*--- Temporal fix, only for now ---*/
			WAitkDyn = max(WAitkDyn, 0.1);

			SetWAitken_Dyn(WAitkDyn);
			if (writeHistFSI){
				historyFile_FSI << " " << endl ;
				historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
				historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
				historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << WAitkDyn ;
			}

		}
		else{

			for (iPoint=0; iPoint<nPoint; iPoint++){

				dispPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
				dispPred_Old = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred_Old();
				dispCalc = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
				dispCalc_Old = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Old();

				for (iDim=0; iDim < nDim; iDim++){

					/*--- Compute the deltaU and deltaU_n+1 ---*/
					deltaU[iDim] = dispCalc_Old[iDim] - dispPred_Old[iDim];
					deltaU_p1[iDim] = dispCalc[iDim] - dispPred[iDim];

					/*--- Compute the difference ---*/
					delta_deltaU[iDim] = deltaU_p1[iDim] - deltaU[iDim];

					/*--- Add numerator and denominator ---*/
					numAitk += deltaU[iDim] * delta_deltaU[iDim];
					denAitk += delta_deltaU[iDim] * delta_deltaU[iDim];

				}

			}

				WAitkDyn = GetWAitken_Dyn();

			if (denAitk > 1E-8){
				WAitkDyn = - 1.0 * WAitkDyn * numAitk / denAitk ;
			}

				WAitkDyn = max(WAitkDyn, 0.1);
				WAitkDyn = min(WAitkDyn, 1.0);

				SetWAitken_Dyn(WAitkDyn);

				if (writeHistFSI){
					historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
					historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
					historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << WAitkDyn << "," ;
				}

		}

	}

	if (writeHistFSI){historyFile_FSI.close();}


}

void CFEASolver::SetAitken_Relaxation(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) {
  
  unsigned long iPoint, iDim;
  unsigned long nPoint, nDim;
  unsigned short RelaxMethod_FSI;
  su2double *dispPred, *dispCalc;
  su2double WAitken;
  su2double CurrentTime=fea_config->GetCurrent_DynTime();
  su2double Static_Time=fea_config->GetStatic_Time();
  
  nPoint = fea_geometry[MESH_0]->GetnPoint();
  nDim = fea_geometry[MESH_0]->GetnDim();
  
  dispPred=new su2double [nDim];
  dispCalc=new su2double [nDim];
  
  RelaxMethod_FSI = fea_config->GetRelaxation_Method_FSI();
  
  /*--- Only when there is movement it makes sense to update the solutions... ---*/
  
  if (CurrentTime > Static_Time) {
    
    if (RelaxMethod_FSI == NO_RELAXATION){
      WAitken = 1.0;
    }
    else if (RelaxMethod_FSI == FIXED_PARAMETER){
      WAitken = fea_config->GetAitkenStatRelax();
    }
    else if (RelaxMethod_FSI == AITKEN_DYNAMIC){
      WAitken = GetWAitken_Dyn();
    }
    else {
      WAitken = 1.0;
      cout << "No relaxation parameter used. " << endl;
    }
    
    
    for (iPoint=0; iPoint<nPoint; iPoint++){
      
      /*--- Retrieve pointers to the predicted and calculated solutions ---*/
      dispPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
      dispCalc = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
      
      /*--- Set predicted solution as the old predicted solution ---*/
      fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred_Old();
      
      /*--- Set calculated solution as the old solution (needed for dynamic Aitken relaxation) ---*/
      fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Old(dispCalc);
      
      /*--- Apply the Aitken relaxation ---*/
      for (iDim=0; iDim < nDim; iDim++){
        dispPred[iDim] = (1.0 - WAitken)*dispPred[iDim] + WAitken*dispCalc[iDim];
      }
      
      /*--- Set obtained solution as the new predicted solution ---*/
      /*--- As dispPred is the pointer to the solution_Pred, we don't need to do this... ---*/
      //fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred(dispPred);
      
    }
    
  }
  
  delete [] dispCalc;
  delete [] dispPred;
  
}

void CFEASolver::Update_StructSolution(CGeometry **fea_geometry, CConfig *fea_config, CSolver ***fea_solution) {
  
  su2double *valSolutionPred;
  
  for (unsigned long iPoint=0; iPoint<fea_geometry[MESH_0]->GetnPoint(); iPoint++){
    
    valSolutionPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
    
    fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution(valSolutionPred);
    
  }
  
}




