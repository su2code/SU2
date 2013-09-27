/*!
 * \file solution_direct_fea.cpp
 * \brief Main subrotuines for solving the FEA equation.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.7
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
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
	unsigned short nMarker, iVar, NodesElement;
  double dull_val;
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
	nDim    = geometry->GetnDim();
	nMarker = config->GetnMarker_All();
	node    = new CVariable*[nPoint];
  if (config->GetUnsteady_Simulation() == STEADY) nVar = nDim;
  else nVar = 2*nDim;
  
	if (nDim == 2) NodesElement = 3;	// Triangles in 2D
	if (nDim == 3) NodesElement = 4;	// Tets in 3D
	
  /*--- Define some auxiliary vectors related to the residual ---*/
  Residual = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
  
  /*--- Define some auxiliary vectors related to the solution ---*/
	Solution   = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  
	/*--- Element aux stiffness matrix definition ---*/
	StiffMatrix_Elem = new double*[NodesElement*nDim];
	for (iVar = 0; iVar < NodesElement*nDim; iVar++) {
		StiffMatrix_Elem[iVar] = new double [NodesElement*nDim];
	}
	
	/*--- Node aux stiffness matrix definition ---*/
	StiffMatrix_Node = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		StiffMatrix_Node[iVar] = new double [nVar];
	}
  
	/*--- Initialization of matrix structures ---*/
  StiffMatrixSpace.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
	StiffMatrixTime.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
  if (rank == MASTER_NODE) cout << "Initialize jacobian structure (Linear Elasticity)." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
  
  /*--- Initialization of linear solver structures ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
	/*--- Computation of gradients by least squares ---*/
	Smatrix = new double* [nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Smatrix[iDim] = new double [nDim];
  
	cvector = new double* [nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		cvector[iVar] = new double [nDim];
  
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
    for(iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
      node[iPoint] = new CFEAVariable(Solution, nDim, nVar, config);
    }
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
	}
  
}

CFEASolver::~CFEASolver(void) {
  
	unsigned short iVar, iDim, NodesElement;
  
	if (nDim == 2) NodesElement = 3;	// Triangles in 2D
	if (nDim == 3) NodesElement = 4;	// Tets in 3D
  
	delete [] Residual;
	delete [] Residual_Max;
	delete [] Solution;
  
	for (iVar = 0; iVar < NodesElement*nDim; iVar++)
		delete [] StiffMatrix_Elem[iVar];
  
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] StiffMatrix_Node[iVar];
  
	delete [] StiffMatrix_Elem;
	delete [] StiffMatrix_Node;
  
	/*--- Computation of gradients by least-squares ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;
  
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
  
}

void CFEASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
	
  /*--- Set residuals and auxiliar variable to zero ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		LinSysRes.SetBlock_Zero(iPoint);
    LinSysAux.SetBlock_Zero(iPoint);
	}
	
	/*--- Set matrix entries to zero ---*/
	StiffMatrixSpace.SetValZero();
	StiffMatrixTime.SetValZero();
	Jacobian.SetValZero();
  
}

void CFEASolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                 CConfig *config, unsigned short iMesh) {
  
  if (config->GetUnsteady_Simulation() != STEADY) {
    
    unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
    double a[3], b[3], c[3], d[3], Area_Local = 0.0, Volume_Local = 0.0, Time_Num;
    double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
    unsigned short iDim, iVar, jVar;
//    double MassMatrix_Elem_2D [6][6] =
//    {{ 2, 0, 1, 0, 1, 0 },
//      { 0, 2, 0, 1, 0, 1 },
//      { 1, 0, 2, 0, 1, 0 },
//      { 0, 1, 0, 2, 0, 1 },
//      { 1, 0, 1, 0, 2, 0 },
//      { 0, 1, 0, 1, 0, 2 }};
//    
//    double MassMatrix_Elem_3D [12][12] =
//    {{ 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0 },
//      { 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0 },
//      { 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1 },
//      { 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0 },
//      { 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0 },
//      { 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1 },
//      { 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0 },
//      { 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0 },
//      { 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1 },
//      { 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0 },
//      { 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0 },
//      { 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2 }};
    
    double Density = config->GetMaterialDensity();
    
    /*--- Numerical time step (this system is uncoditional stable... a very big number can be used) ---*/
    if (config->GetUnsteady_Simulation() == TIME_STEPPING) Time_Num = config->GetDelta_UnstTimeND();
    else Time_Num = 1.0;
    if (config->GetUnsteady_Simulation() == STEADY) Time_Num = 1.0;
    
    /*--- Loop through elements to compute contributions from the matrix
     blocks involving time. These contributions are also added to the
     Jacobian w/ the time step. Spatial source terms are also computed. ---*/
    
    for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      
      /*--- Get node numbers and their coordinate vectors ---*/
      Point_0 = geometry->elem[iElem]->GetNode(0);	Coord_0 = geometry->node[Point_0]->GetCoord();
      Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1 = geometry->node[Point_1]->GetCoord();
      Point_2 = geometry->elem[iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord();
      if (nDim == 3) { Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3 = geometry->node[Point_3]->GetCoord(); }
      
      if (nDim == 2) {
        for (iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = Coord_0[iDim]-Coord_2[iDim];
          b[iDim] = Coord_1[iDim]-Coord_2[iDim];
        }
        /*--- Compute element area ---*/
        Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
      }
      if (nDim == 3) {
        for (iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = Coord_0[iDim]-Coord_2[iDim];
          b[iDim] = Coord_1[iDim]-Coord_2[iDim];
          c[iDim] = Coord_3[iDim]-Coord_2[iDim];
        }
        d[0] = a[1]*b[2]-a[2]*b[1];
        d[1] = -(a[0]*b[2]-a[2]*b[0]);
        d[2] = a[0]*b[1]-a[1]*b[0];
        /*--- Compute element volume ---*/
        Volume_Local = fabs(c[0]*d[0] + c[1]*d[1] + c[2]*d[2])/6.0;
      }
      
      /*----------------------------------------------------------------*/
      /*--- Block contributions to the Jacobian (includes time step) ---*/
      /*----------------------------------------------------------------*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          StiffMatrix_Node[iVar][jVar] = 0.0;
      
      /*--- Diagonal value identity matrix ---*/
      if (nDim == 2) {
        StiffMatrix_Node[0][0] = 1.0/Time_Num;
        StiffMatrix_Node[1][1] = 1.0/Time_Num;
      }
      else {
        StiffMatrix_Node[0][0] = 1.0/Time_Num;
        StiffMatrix_Node[1][1] = 1.0/Time_Num;
        StiffMatrix_Node[2][2] = 1.0/Time_Num;
      }
      
      /*--- Diagonal value ---*/
      if (nDim == 2) {
        StiffMatrix_Node[2][2] = Density*(2.0/12.0)*(Area_Local/Time_Num);
        StiffMatrix_Node[3][3] = Density*(2.0/12.0)*(Area_Local/Time_Num);
      }
      else {
        StiffMatrix_Node[3][3] = Density*(2.0/20.0)*(Volume_Local/Time_Num);
        StiffMatrix_Node[4][4] = Density*(2.0/20.0)*(Volume_Local/Time_Num);
        StiffMatrix_Node[5][5] = Density*(2.0/20.0)*(Volume_Local/Time_Num);
      }
      Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
      Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
      Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
      if (nDim == 3) Jacobian.AddBlock(Point_3, Point_3, StiffMatrix_Node);
      
      /*--- Off Diagonal value ---*/
      if (nDim == 2) {
        StiffMatrix_Node[2][2] = Density*(1.0/12.0)*(Area_Local/Time_Num);
        StiffMatrix_Node[3][3] = Density*(1.0/12.0)*(Area_Local/Time_Num);
      }
      else {
        StiffMatrix_Node[3][3] = Density*(1.0/20.0)*(Volume_Local/Time_Num);
        StiffMatrix_Node[4][4] = Density*(1.0/20.0)*(Volume_Local/Time_Num);
        StiffMatrix_Node[5][5] = Density*(1.0/20.0)*(Volume_Local/Time_Num);
      }
      Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
      Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
      Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
      Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
      Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
      Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
      
      if (nDim == 3) {
        Jacobian.AddBlock(Point_0, Point_3, StiffMatrix_Node);
        Jacobian.AddBlock(Point_1, Point_3, StiffMatrix_Node);
        Jacobian.AddBlock(Point_2, Point_3, StiffMatrix_Node);
        Jacobian.AddBlock(Point_3, Point_0, StiffMatrix_Node);
        Jacobian.AddBlock(Point_3, Point_1, StiffMatrix_Node);
        Jacobian.AddBlock(Point_3, Point_2, StiffMatrix_Node);
      }
      
//      /*--- Add volumetric forces as source term (gravity and coriollis) ---*/
//      double RhoG= Density*STANDART_GRAVITY;
//      
//      unsigned short NodesElement = 3; // Triangles in 2D
//      if (nDim == 3) NodesElement = 4;	// Tets in 3D
//      
//      double *ElemForce = new double [nDim*NodesElement];
//      double *ElemResidual = new double [nDim*NodesElement];
//      
//      if (nDim == 2) {
//        ElemForce[0] = 0.0;	ElemForce[1] = -RhoG;
//        ElemForce[2] = 0.0;	ElemForce[3] = -RhoG;
//        ElemForce[4] = 0.0;	ElemForce[5] = -RhoG;
//        
//        for (iVar = 0; iVar < nDim*NodesElement; iVar++) {
//          ElemResidual[iVar] = 0.0;
//          for (jVar = 0; jVar < nDim*NodesElement; jVar++) {
//            ElemResidual[iVar] += MassMatrix_Elem_2D[iVar][jVar]*ElemForce[jVar]*Area_Local/12.0;
//          }
//        }
//        
//        Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = ElemResidual[0]; Residual[3] = ElemResidual[1];
//        LinSysRes.AddBlock(Point_0, Residual);
//        
//        Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = ElemResidual[2]; Residual[3] = ElemResidual[3];
//        LinSysRes.AddBlock(Point_1, Residual);
//        
//        Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = ElemResidual[4]; Residual[3] = ElemResidual[5];
//        LinSysRes.AddBlock(Point_2, Residual);
//        
//      }
//      
//      else {
//        
//        ElemForce[0] = 0.0; ElemForce[1] = 0.0;		ElemForce[2] = -RhoG;
//        ElemForce[3] = 0.0; ElemForce[4] = 0.0;		ElemForce[5] = -RhoG;
//        ElemForce[6] = 0.0; ElemForce[7] = 0.0;		ElemForce[8] = -RhoG;
//        ElemForce[9] = 0.0; ElemForce[10] = 0.0;	ElemForce[11] = -RhoG;
//        
//        for (iVar = 0; iVar < nDim*NodesElement; iVar++) {
//          ElemResidual[iVar] = 0.0;
//          for (jVar = 0; jVar < nDim*NodesElement; jVar++) {
//            ElemResidual[iVar] += MassMatrix_Elem_3D[iVar][jVar]*ElemForce[jVar]*Volume_Local/20.0;
//          }
//        }
//        
//        Residual[0] = 0.0;							Residual[1] = 0.0;							Residual[2] = 0.0;
//        Residual[3] = ElemResidual[0];	Residual[4] = ElemResidual[1];	Residual[5] = ElemResidual[2];
//        LinSysRes.AddBlock(Point_0, Residual);
//        
//        Residual[0] = 0.0;							Residual[1] = 0.0;							Residual[2] = 0.0;
//        Residual[3] = ElemResidual[3];	Residual[4] = ElemResidual[4];	Residual[5] = ElemResidual[5];
//        LinSysRes.AddBlock(Point_1, Residual);
//        
//        Residual[0] = 0.0;							Residual[1] = 0.0;							Residual[2] = 0.0;
//        Residual[3] = ElemResidual[6];	Residual[4] = ElemResidual[7];	Residual[5] = ElemResidual[8];
//        LinSysRes.AddBlock(Point_2, Residual);
//        
//        Residual[0] = 0.0;							Residual[1] = 0.0;							Residual[2] = 0.0;
//        Residual[3] = ElemResidual[9];	Residual[4] = ElemResidual[10];	Residual[5] = ElemResidual[11];
//        LinSysRes.AddBlock(Point_3, Residual);
//        
//      }
//      
//      delete [] ElemForce;
//      delete [] ElemResidual;
      
    }
    
  }
  
}

void CFEASolver::Galerkin_Method(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                 CConfig *config, unsigned short iMesh) {
  
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0, iPoint, total_index;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
	unsigned short iVar, jVar;
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
		Point_0 = geometry->elem[iElem]->GetNode(0);	Coord_0 = geometry->node[Point_0]->GetCoord();
		Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1 = geometry->node[Point_1]->GetCoord();
		Point_2 = geometry->elem[iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord();
		if (nDim == 3) { Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3 = geometry->node[Point_3]->GetCoord(); }
    
		if (nDim == 2) numerics->SetCoord(Coord_0, Coord_1, Coord_2);
		if (nDim == 3) numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
		
		numerics->ComputeResidual(StiffMatrix_Elem, config);
    
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				StiffMatrix_Node[iVar][jVar] = 0.0;
    
    if (config->GetUnsteady_Simulation() == STEADY) {
      
      if (nDim == 2) {
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[0][1];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][0];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[1][1];
        Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[0][3];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][2];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[1][3];
        Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][4];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[0][5];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][4];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[1][5];
        Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[2][1];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[3][0];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[3][1];
        Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[2][3];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[3][2];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[3][3];
        Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][4];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[2][5];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[3][4];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[3][5];
        Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[4][0];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[4][1];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[5][0];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[5][1];
        Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[4][2];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[4][3];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[5][2];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[5][3];
        Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[4][4];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[4][5];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[5][4];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[5][5];
        Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
        
      }
      
      else {
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[0][1];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[0][2];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][0];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[1][1];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[1][2];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[2][0];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[2][1];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[2][2];
        Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][3];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[0][4];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[0][5];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][3];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[1][4];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[1][5];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[2][3];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[2][4];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[2][5];
        Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][6];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[0][7];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[0][8];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][6];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[1][7];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[1][8];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[2][6];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[2][7];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[2][8];
        Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][9];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[0][10];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[0][11];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][9];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[1][10];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[1][11];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[2][9];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[2][10];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[2][11];
        Jacobian.AddBlock(Point_0, Point_3, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][0];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[3][1];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[3][2];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[4][0];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[4][1];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[4][2];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[5][0];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[5][1];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[5][2];
        Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][3];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[3][4];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[3][5];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[4][3];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[4][4];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[4][5];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[5][3];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[5][4];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[5][5];
        Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][6];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[3][7];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[3][8];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[4][6];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[4][7];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[4][8];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[5][6];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[5][7];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[5][8];
        Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][9];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[3][10];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[3][11];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[4][9];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[4][10];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[4][11];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[5][9];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[5][10];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[5][11];
        Jacobian.AddBlock(Point_1, Point_3, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[6][0];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[6][1];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[6][2];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[7][0];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[7][1];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[7][2];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[8][0];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[8][1];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[8][2];
        Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[6][3];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[6][4];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[6][5];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[7][3];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[7][4];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[7][5];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[8][3];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[8][4];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[8][5];
        Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[6][6];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[6][7];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[6][8];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[7][6];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[7][7];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[7][8];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[8][6];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[8][7];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[8][8];
        Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[6][9];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[6][10];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[6][11];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[7][9];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[7][10];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[7][11];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[8][9];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[8][10];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[8][11];
        Jacobian.AddBlock(Point_2, Point_3, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[9][0];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[9][1];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[9][2];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[10][0];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[10][1];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[10][2];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[11][0];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[11][1];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[11][2];
        Jacobian.AddBlock(Point_3, Point_0, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[9][3];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[9][4];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[9][5];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[10][3];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[10][4];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[10][5];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[11][3];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[11][4];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[11][5];
        Jacobian.AddBlock(Point_3, Point_1, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[9][6];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[9][7];	StiffMatrix_Node[0][2] = StiffMatrix_Elem[9][8];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[10][6];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[10][7];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[10][8];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[11][6];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[11][7];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[11][8];
        Jacobian.AddBlock(Point_3, Point_2, StiffMatrix_Node);
        
        StiffMatrix_Node[0][0] = StiffMatrix_Elem[9][9];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[9][10];		StiffMatrix_Node[0][2] = StiffMatrix_Elem[9][11];
        StiffMatrix_Node[1][0] = StiffMatrix_Elem[10][9];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[10][10];	StiffMatrix_Node[1][2] = StiffMatrix_Elem[10][11];
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[11][9];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[11][10];	StiffMatrix_Node[2][2] = StiffMatrix_Elem[11][11];
        Jacobian.AddBlock(Point_3, Point_3, StiffMatrix_Node);
        
      }
      
    }
    
    else {
            
      if (nDim == 2) {
        StiffMatrix_Node[0][2] = -1.0;
        StiffMatrix_Node[1][3] = -1.0;

        StiffMatrix_Node[2][0] = StiffMatrix_Elem[0][0];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[0][1];
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[1][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[1][1];
        StiffMatrixSpace.AddBlock(Point_0, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
        
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[0][2];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[0][3];
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[1][2];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[1][3];
        StiffMatrixSpace.AddBlock(Point_0, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
        
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[0][4];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[0][5];
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[1][4];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[1][5];
        StiffMatrixSpace.AddBlock(Point_0, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
        
        StiffMatrix_Node[2][0] = StiffMatrix_Elem[2][0];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[2][1];
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][1];
        StiffMatrixSpace.AddBlock(Point_1, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);

        StiffMatrix_Node[2][0] = StiffMatrix_Elem[2][2];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[2][3];
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][2];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][3];
        StiffMatrixSpace.AddBlock(Point_1, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);

        StiffMatrix_Node[2][0] = StiffMatrix_Elem[2][4];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[2][5];
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][4];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][5];
        StiffMatrixSpace.AddBlock(Point_1, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);

        StiffMatrix_Node[2][0] = StiffMatrix_Elem[4][0];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[4][1];
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[5][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[5][1];
        StiffMatrixSpace.AddBlock(Point_2, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);

        StiffMatrix_Node[2][0] = StiffMatrix_Elem[4][2];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[4][3];
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[5][2];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[5][3];
        StiffMatrixSpace.AddBlock(Point_2, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);

        StiffMatrix_Node[2][0] = StiffMatrix_Elem[4][4];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[4][5];
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[5][4];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[5][5];
        StiffMatrixSpace.AddBlock(Point_2, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);

      }
      else {
        StiffMatrix_Node[0][3] = -1.0;
        StiffMatrix_Node[1][4] = -1.0;
        StiffMatrix_Node[2][5] = -1.0;

        StiffMatrix_Node[3][0] = StiffMatrix_Elem[0][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[0][1];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[0][2];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[1][0];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[1][1];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[1][2];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[2][0];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[2][1];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[2][2];
        StiffMatrixSpace.AddBlock(Point_0, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
        
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[0][3];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[0][4];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[0][5];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[1][3];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[1][4];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[1][5];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[2][3];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[2][4];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[2][5];
        StiffMatrixSpace.AddBlock(Point_0, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
        
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[0][6];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[0][7];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[0][8];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[1][6];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[1][7];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[1][8];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[2][6];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[2][7];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[2][8];
        StiffMatrixSpace.AddBlock(Point_0, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
        
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[0][9];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[0][10];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[0][11];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[1][9];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[1][10];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[1][11];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[2][9];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[2][10];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[2][11];
        StiffMatrixSpace.AddBlock(Point_0, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_3, StiffMatrix_Node);
        
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][1];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[3][2];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[4][0];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[4][1];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[4][2];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[5][0];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[5][1];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[5][2];
        StiffMatrixSpace.AddBlock(Point_1, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
        
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][3];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][4];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[3][5];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[4][3];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[4][4];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[4][5];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[5][3];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[5][4];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[5][5];
        StiffMatrixSpace.AddBlock(Point_1, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);

        StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][6];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][7];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[3][8];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[4][6];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[4][7];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[4][8];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[5][6];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[5][7];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[5][8];
        StiffMatrixSpace.AddBlock(Point_1, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
        
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][9];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][10];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[3][11];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[4][9];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[4][10];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[4][11];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[5][9];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[5][10];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[5][11];
        StiffMatrixSpace.AddBlock(Point_1, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_3, StiffMatrix_Node);

        StiffMatrix_Node[3][0] = StiffMatrix_Elem[6][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[6][1];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[6][2];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[7][0];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[7][1];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[7][2];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[8][0];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[8][1];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[8][2];
        StiffMatrixSpace.AddBlock(Point_2, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);

        StiffMatrix_Node[3][0] = StiffMatrix_Elem[6][3];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[6][4];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[6][5];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[7][3];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[7][4];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[7][5];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[8][3];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[8][4];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[8][5];
        StiffMatrixSpace.AddBlock(Point_2, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);

        StiffMatrix_Node[3][0] = StiffMatrix_Elem[6][6];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[6][7];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[6][8];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[7][6];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[7][7];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[7][8];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[8][6];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[8][7];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[8][8];
        StiffMatrixSpace.AddBlock(Point_2, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);

        StiffMatrix_Node[3][0] = StiffMatrix_Elem[6][9];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[6][10];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[6][11];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[7][9];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[7][10];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[7][11];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[8][9];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[8][10];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[8][11];
        StiffMatrixSpace.AddBlock(Point_2, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_3, StiffMatrix_Node);

        StiffMatrix_Node[3][0] = StiffMatrix_Elem[9][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[9][1];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[9][2];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[10][0];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[10][1];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[10][2];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[11][0];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[11][1];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[11][2];
        StiffMatrixSpace.AddBlock(Point_3, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_3, Point_0, StiffMatrix_Node);
        
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[9][3];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[9][4];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[9][5];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[10][3];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[10][4];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[10][5];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[11][3];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[11][4];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[11][5];
        StiffMatrixSpace.AddBlock(Point_3, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_3, Point_1, StiffMatrix_Node);
        
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[9][6];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[9][7];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[9][8];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[10][6];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[10][7];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[10][8];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[11][6];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[11][7];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[11][8];
        StiffMatrixSpace.AddBlock(Point_3, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_3, Point_2, StiffMatrix_Node);
        
        StiffMatrix_Node[3][0] = StiffMatrix_Elem[9][9];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[9][10];		StiffMatrix_Node[3][2] = StiffMatrix_Elem[9][11];
        StiffMatrix_Node[4][0] = StiffMatrix_Elem[10][9];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[10][10];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[10][11];
        StiffMatrix_Node[5][0] = StiffMatrix_Elem[11][9];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[11][10];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[11][11];
        StiffMatrixSpace.AddBlock(Point_3, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_3, Point_3, StiffMatrix_Node);

      }
      
    }
    
	}
  
  if (config->GetUnsteady_Simulation() != STEADY) {
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        LinSysSol[total_index] = node[iPoint]->GetSolution(iVar);
        LinSysAux[total_index] = 0.0;
      }
    
    StiffMatrixSpace.MatrixVectorProduct(LinSysSol, LinSysAux);
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Residual[iVar] = LinSysAux[total_index];
      }
      LinSysRes.SubtractBlock(iPoint, Residual);
    }
  }
  
}

void CFEASolver::BC_Normal_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) {
	unsigned long iPoint, iVertex, total_index;
	unsigned short iVar, iDim;
  double *Normal, Area, UnitaryNormal[3];
	
	double TotalDispl = config->GetDispl_Value(config->GetMarker_All_Tag(val_marker));
	
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
	
	double a[3], b[3];
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
	double Length_Elem = 0.0, Area_Elem = 0.0, Normal_Elem[3] = {0.0, 0.0, 0.0};
	unsigned short iDim;
	
	double TotalLoad = 100*config->GetLoad_Value(config->GetMarker_All_Tag(val_marker));
	
	for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
		Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);                   Coord_0 = geometry->node[Point_0]->GetCoord();
		Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);                   Coord_1 = geometry->node[Point_1]->GetCoord();
		if (nDim == 3) { Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord(); }
    
		/*--- Compute area (3D), and length of the surfaces (2D) ---*/
		if (nDim == 2) {
			for (iDim = 0; iDim < nDim; iDim++)
				a[iDim] = Coord_0[iDim]-Coord_1[iDim];
			Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);
      
      Normal_Elem[0] = -a[1];
			Normal_Elem[1] = a[0];
      
		}
		if (nDim == 3) {
			for (iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}
			Area_Elem = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
      
      Normal_Elem[0] = 0.5*(a[1]*b[2]-a[2]*b[1]);
			Normal_Elem[1] = -0.5*(a[0]*b[2]-a[2]*b[0]);
			Normal_Elem[2] = 0.5*(a[0]*b[1]-a[1]*b[0]);
		}
		
    if (config->GetUnsteady_Simulation() == STEADY) {
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
    else {
      if (nDim == 2) {
        Residual[0] = 0.0; Residual[1] = 0.0;
        Residual[2] = (1.0/2.0)*TotalLoad*Normal_Elem[0]; Residual[3] = (1.0/2.0)*TotalLoad*Normal_Elem[1];
        LinSysRes.AddBlock(Point_0, Residual);
        Residual[0] = 0.0; Residual[1] = 0.0;
        Residual[2] = (1.0/2.0)*TotalLoad*Normal_Elem[0]; Residual[3] = (1.0/2.0)*TotalLoad*Normal_Elem[1];
        LinSysRes.AddBlock(Point_1, Residual);
      }
      else {
        Residual[0] = 0.0;                                          Residual[1] = 0.0;                                          Residual[2] = 0.0;
        Residual[3] = (1.0/3.0)*TotalLoad*Normal_Elem[0]; Residual[4] = (1.0/3.0)*TotalLoad*Normal_Elem[1]; Residual[5] = (1.0/3.0)*TotalLoad*Normal_Elem[2];
        LinSysRes.AddBlock(Point_0, Residual);
        
        Residual[0] = 0.0;                                          Residual[1] = 0.0;                                          Residual[2] = 0.0;
        Residual[3] = (1.0/3.0)*TotalLoad*Normal_Elem[0]; Residual[4] = (1.0/3.0)*TotalLoad*Normal_Elem[1]; Residual[5] = (1.0/3.0)*TotalLoad*Normal_Elem[2];
        LinSysRes.AddBlock(Point_1, Residual);
        
        Residual[0] = 0.0;                                          Residual[1] = 0.0;                                          Residual[2] = 0.0;
        Residual[3] = (1.0/3.0)*TotalLoad*Normal_Elem[0]; Residual[4] = (1.0/3.0)*TotalLoad*Normal_Elem[1]; Residual[5] = (1.0/3.0)*TotalLoad*Normal_Elem[2];
        LinSysRes.AddBlock(Point_2, Residual);
      }
    }
		
	}
	
}

void CFEASolver::BC_Flow_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                              unsigned short val_marker) {
	
  //	double a[3], b[3], Press_0 = 0.0, Press_1 = 0.0, Press_2 = 0.0, Press_Elem;
  //	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
  //	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
  //	double Length_Elem, Area_Elem, Normal_Elem[3] = {0.0,0.0,0.0};
  //	unsigned short iDim;
  //
  //	for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
  //		Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);	Coord_0 = geometry->node[Point_0]->GetCoord(); Press_0 = 1.0; //node[Point_0]->GetPressure();
  //		Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);	Coord_1 = geometry->node[Point_1]->GetCoord(); Press_1 = 1.0; //node[Point_1]->GetPressure();
  //		if (nDim == 3) { Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord(); Press_2 = 1.0; }//node[Point_2]->GetPressure(); }
  //
  //		/*--- Compute area (3D), and length of the surfaces (2D) ---*/
  //		if (nDim == 2) {
  //			for (iDim = 0; iDim < nDim; iDim++)
  //				a[iDim] = Coord_0[iDim]-Coord_1[iDim];
  //			Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);
  //
  //			Normal_Elem[0] = -a[1];
  //			Normal_Elem[1] = a[0];
  //
  //			Press_Elem = 0.5*(Press_0 + Press_1);
  //		}
  //		else {
  //
  //			for (iDim = 0; iDim < nDim; iDim++) {
  //				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
  //				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
  //			}
  //			Area_Elem = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  //
  //			Normal_Elem[0] = 0.5*(a[1]*b[2]-a[2]*b[1]);
  //			Normal_Elem[1] = -0.5*(a[0]*b[2]-a[2]*b[0]);
  //			Normal_Elem[2] = 0.5*(a[0]*b[1]-a[1]*b[0]);
  //
  //			Press_Elem = (Press_0 + Press_1 + Press_2)/3.0;
  //		}
  //
  //		if (nDim == 2) {
  //			Residual[0] = 0.0; Residual[1] = 0.0;
  //			Residual[2] = (1.0/2.0)*Press_Elem*Normal_Elem[0]; Residual[3] = (1.0/2.0)*Press_Elem*Normal_Elem[1];
  //			LinSysRes.AddBlock(Point_0, Residual);
  //
  //			Residual[0] = 0.0; Residual[1] = 0.0;
  //			Residual[2] = (1.0/2.0)*Press_Elem*Normal_Elem[0]; Residual[3] = (1.0/2.0)*Press_Elem*Normal_Elem[1];
  //			LinSysRes.AddBlock(Point_1, Residual);
  //		}
  //		else {
  //			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0;
  //			Residual[3] = (1.0/3.0)*Press_Elem*Normal_Elem[0]; Residual[4] = (1.0/3.0)*Press_Elem*Normal_Elem[1]; Residual[5] = (1.0/3.0)*Press_Elem*Normal_Elem[2];
  //			LinSysRes.AddBlock(Point_0, Residual);
  //
  //			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0;
  //			Residual[3] = (1.0/3.0)*Press_Elem*Normal_Elem[0]; Residual[4] = (1.0/3.0)*Press_Elem*Normal_Elem[1]; Residual[5] = (1.0/3.0)*Press_Elem*Normal_Elem[2];
  //			LinSysRes.AddBlock(Point_1, Residual);
  //
  //			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0;
  //			Residual[3] = (1.0/3.0)*Press_Elem*Normal_Elem[0]; Residual[4] = (1.0/3.0)*Press_Elem*Normal_Elem[1]; Residual[5] = (1.0/3.0)*Press_Elem*Normal_Elem[2];
  //			LinSysRes.AddBlock(Point_2, Residual);
  //		}
  //
  //    /*
  //
  //     double StiffMatrix_BoundElem[9][9], Vector_BoundElem[9], a[3], b[3], LoadNode[9], Press_0, Press_1, Press_2, Press_Elem;
  //     unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Vertex_0, Vertex_1, Vertex_2;
  //     double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, Length_Elem, Area_Elem, *Normal, Area_0, Area_1, Area_2, Normal_Elem[3];
  //     unsigned short iVar, jVar, iDim;
  //
  //     double LocalLoad = config->GetLoad_Value(config->GetMarker_All_Tag(val_marker));
  //
  //     if (nDim == 2) {
  //
  //     StiffMatrix_BoundElem[0][0] = (2.0/6.0)*Length_Elem;	StiffMatrix_BoundElem[0][1] = 0.0;										StiffMatrix_BoundElem[0][2] = (1.0/6.0)*Length_Elem;			StiffMatrix_BoundElem[0][3] = 0.0;
  //     StiffMatrix_BoundElem[1][0] = 0.0;										StiffMatrix_BoundElem[1][1] = (2.0/6.0)*Length_Elem;	StiffMatrix_BoundElem[1][2] = 0.0;												StiffMatrix_BoundElem[1][3] = (1.0/6.0)*Length_Elem;
  //     StiffMatrix_BoundElem[2][0] = (1.0/6.0)*Length_Elem;	StiffMatrix_BoundElem[2][1] = 0.0;										StiffMatrix_BoundElem[2][2] = (2.0/6.0)*Length_Elem;			StiffMatrix_BoundElem[2][3] = 0.0;
  //     StiffMatrix_BoundElem[3][0] = 0.0;										StiffMatrix_BoundElem[3][1] = (1.0/6.0)*Length_Elem;	StiffMatrix_BoundElem[3][2] = 0.0;												StiffMatrix_BoundElem[3][3] = (2.0/6.0)*Length_Elem;
  //
  //     LoadNode[0] = 0.0; LoadNode[1] = LocalLoad; LoadNode[2] = 0.0;	LoadNode[3] = LocalLoad;
  //
  //     for (iVar = 0; iVar < 4; iVar++) {
  //     Vector_BoundElem[iVar] = 0.0;
  //     for (jVar = 0; jVar < 4; jVar++) {
  //     Vector_BoundElem[iVar] += StiffMatrix_BoundElem[iVar][jVar]*LoadNode[jVar];
  //     }
  //     }
  //
  //     Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = Vector_BoundElem[0]; Residual[3] = Vector_BoundElem[1];
  //     node[Point_0]->LinSysRes.AddBlock(Residual);
  //
  //     Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = Vector_BoundElem[2]; Residual[3] = Vector_BoundElem[3];
  //     node[Point_1]->LinSysRes.AddBlock(Residual);
  //
  //     }
  //
  //     if (nDim == 3) {
  //
  //     StiffMatrix_BoundElem[0][0] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[0][1] = 0.0;										StiffMatrix_BoundElem[0][2] = 0.0;										StiffMatrix_BoundElem[0][3] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[0][4] = 0.0;										StiffMatrix_BoundElem[0][5] = 0.0;										StiffMatrix_BoundElem[0][6] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[0][7] = 0.0;										StiffMatrix_BoundElem[0][8] = 0.0;
  //     StiffMatrix_BoundElem[1][0] = 0.0;										StiffMatrix_BoundElem[1][1] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[1][2] = 0.0;										StiffMatrix_BoundElem[1][3] = 0.0;										StiffMatrix_BoundElem[1][4] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[1][5] = 0.0;										StiffMatrix_BoundElem[1][6] = 0.0;										StiffMatrix_BoundElem[1][7] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[1][8] = 0.0;
  //     StiffMatrix_BoundElem[2][0] = 0.0;										StiffMatrix_BoundElem[2][1] = 0.0;										StiffMatrix_BoundElem[2][2] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[2][3] = 0.0;										StiffMatrix_BoundElem[2][4] = 0.0;										StiffMatrix_BoundElem[2][5] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[2][6] = 0.0;										StiffMatrix_BoundElem[2][7] = 0.0;										StiffMatrix_BoundElem[2][8] = (1.0/12.0)*Area_Elem;
  //     StiffMatrix_BoundElem[3][0] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[3][1] = 0.0;										StiffMatrix_BoundElem[3][2] = 0.0;										StiffMatrix_BoundElem[3][3] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[3][4] = 0.0;										StiffMatrix_BoundElem[3][5] = 0.0;										StiffMatrix_BoundElem[3][6] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[3][7] = 0.0;										StiffMatrix_BoundElem[3][8] = 0.0;
  //     StiffMatrix_BoundElem[4][0] = 0.0;										StiffMatrix_BoundElem[4][1] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[4][2] = 0.0;										StiffMatrix_BoundElem[4][3] = 0.0;										StiffMatrix_BoundElem[4][4] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[4][5] = 0.0;										StiffMatrix_BoundElem[4][6] = 0.0;										StiffMatrix_BoundElem[4][7] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[4][8] = 0.0;
  //     StiffMatrix_BoundElem[5][0] = 0.0;										StiffMatrix_BoundElem[5][1] = 0.0;										StiffMatrix_BoundElem[5][2] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[5][3] = 0.0;										StiffMatrix_BoundElem[5][4] = 0.0;										StiffMatrix_BoundElem[5][5] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[5][6] = 0.0;										StiffMatrix_BoundElem[5][7] = 0.0;										StiffMatrix_BoundElem[5][8] = (1.0/12.0)*Area_Elem;
  //     StiffMatrix_BoundElem[6][0] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[6][1] = 0.0;										StiffMatrix_BoundElem[6][2] = 0.0;										StiffMatrix_BoundElem[6][3] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[6][4] = 0.0;										StiffMatrix_BoundElem[6][5] = 0.0;										StiffMatrix_BoundElem[6][6] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[6][7] = 0.0;										StiffMatrix_BoundElem[6][8] = 0.0;
  //     StiffMatrix_BoundElem[7][0] = 0.0;										StiffMatrix_BoundElem[7][1] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[7][2] = 0.0;										StiffMatrix_BoundElem[7][3] = 0.0;										StiffMatrix_BoundElem[7][4] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[7][5] = 0.0;										StiffMatrix_BoundElem[7][6] = 0.0;										StiffMatrix_BoundElem[7][7] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[7][8] = 0.0;
  //     StiffMatrix_BoundElem[8][0] = 0.0;										StiffMatrix_BoundElem[8][1] = 0.0;										StiffMatrix_BoundElem[8][2] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[8][3] = 0.0;										StiffMatrix_BoundElem[8][4] = 0.0;										StiffMatrix_BoundElem[8][5] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[8][6] = 0.0;										StiffMatrix_BoundElem[8][7] = 0.0;										StiffMatrix_BoundElem[8][8] = (2.0/12.0)*Area_Elem;
  //
  //     LoadNode[0] = 0.0;	LoadNode[1] = 0.0;	LoadNode[2] = LocalLoad;
  //     LoadNode[3] = 0.0;	LoadNode[4] = 0.0;	LoadNode[5] = LocalLoad;
  //     LoadNode[6] = 0.0;	LoadNode[7] = 0.0;	LoadNode[8] = LocalLoad;
  //
  //     for (iVar = 0; iVar < 9; iVar++) {
  //     Vector_BoundElem[iVar] = 0.0;
  //     for (jVar = 0; jVar < 9; jVar++) {
  //     Vector_BoundElem[iVar] += StiffMatrix_BoundElem[iVar][jVar]*LoadNode[jVar];
  //     }
  //     }
  //
  //     Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0; Residual[3] = Vector_BoundElem[0]; Residual[4] = Vector_BoundElem[1]; Residual[5] = Vector_BoundElem[2];
  //     node[Point_0]->LinSysRes.AddBlock(Residual);
  //
  //     Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0; Residual[3] = Vector_BoundElem[3]; Residual[4] = Vector_BoundElem[4]; Residual[5] = Vector_BoundElem[5];
  //     node[Point_1]->LinSysRes.AddBlock(Residual);
  //
  //     Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0; Residual[3] = Vector_BoundElem[6]; Residual[4] = Vector_BoundElem[7]; Residual[5] = Vector_BoundElem[8];
  //     node[Point_2]->LinSysRes.AddBlock(Residual);
  //
  //     }
  //     */
  //
  //	}
	
}


void CFEASolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  unsigned long iPoint;
  double Strain_xx, Strain_yy, Strain_xy, Strain_zz, Strain_xz, Strain_yz, **Stress, VonMises_Stress, MaxVonMises_Stress = 0.0, Strain_Trace;
  
  double E = config->GetElasticyMod();
  double Nu = config->GetPoissonRatio();
  double Mu = E / (2.0*(1.0 + Nu));
  double Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));		// For plane strain and 3-D
  
	/*--- Compute the gradient of the displacement ---*/
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  
  /*--- Compute the strains and streeses ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    Strain_xx = node[iPoint]->GetGradient(0,0);
    Strain_yy = node[iPoint]->GetGradient(1,1);
    Strain_xy = 0.5*(node[iPoint]->GetGradient(0,1) + node[iPoint]->GetGradient(1,0));
    Strain_Trace = Strain_xx + Strain_yy;
    
    if (geometry->GetnDim() == 3) {
      Strain_zz = node[iPoint]->GetGradient(2,2);
      Strain_xz = 0.5*(node[iPoint]->GetGradient(0,2) + node[iPoint]->GetGradient(2,0));
      Strain_yz = 0.5*(node[iPoint]->GetGradient(1,2) + node[iPoint]->GetGradient(2,1));
      Strain_Trace += Strain_zz;
    }
    
    node[iPoint]->SetStress(0, 0, 2.0*Mu*Strain_xx + Lambda*Strain_Trace);
    node[iPoint]->SetStress(1, 1, 2.0*Mu*Strain_yy + Lambda*Strain_Trace);
    node[iPoint]->SetStress(0, 1, 2.0*Mu*Strain_xy);
    node[iPoint]->SetStress(1, 0, 2.0*Mu*Strain_xy);
    
    if (geometry->GetnDim() == 3) {
      node[iPoint]->SetStress(2, 2, 2.0*Mu*Strain_zz + Lambda*Strain_Trace);
      node[iPoint]->SetStress(0, 2, 2.0*Mu*Strain_xz);
      node[iPoint]->SetStress(2, 0, 2.0*Mu*Strain_xz);
      node[iPoint]->SetStress(1, 2, 2.0*Mu*Strain_yz);
      node[iPoint]->SetStress(2, 1, 2.0*Mu*Strain_yz);
    }
    
    /*---Compute Von Mises criteria ---*/
    Stress = node[iPoint]->GetStress();
    
    if (geometry->GetnDim() == 2) {
      VonMises_Stress = sqrt(  Stress[0][0]*Stress[0][0]
                             - Stress[0][0]*Stress[1][1]
                             + Stress[1][1]*Stress[1][1]
                             + 3.0*Stress[0][1]*Stress[0][1]
                             );
    }
    else {
      VonMises_Stress = sqrt(0.5*(  pow(Stress[0][0] - Stress[1][1], 2.0)
                                  + pow(Stress[1][1] - Stress[2][2], 2.0)
                                  + pow(Stress[2][2] - Stress[0][0], 2.0)
                                  + 6.0*(Stress[0][1]*Stress[0][1]+Stress[1][2]*Stress[1][2]+Stress[2][0]*Stress[2][0])
                                  ));
    }
    
    /*---Store the Von Mises criteria ---*/
    node[iPoint]->SetVonMises_Stress(VonMises_Stress);
    
    /*--- Compute the maximum value of the Von Mises Stress ---*/
    MaxVonMises_Stress = max(MaxVonMises_Stress, VonMises_Stress);
    
  }
  
#ifndef NO_MPI
  
  /*--- Compute MaxVonMises_Stress using all the nodes ---*/
  double MyMaxVonMises_Stress = MaxVonMises_Stress; MaxVonMises_Stress = 0.0;
  MPI::COMM_WORLD.Allreduce(&MyMaxVonMises_Stress, &MaxVonMises_Stress, 1, MPI::DOUBLE, MPI::MAX);
  
#endif
  
  /*--- Set the value of the MaxVonMises_Stress as the CFEA coeffient ---*/
  Total_CFEA = MaxVonMises_Stress;
  
}

void CFEASolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep,
                                      unsigned short iMesh, unsigned short RunTime_EqSystem) {
	
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
	double a[3], b[3], c[3], d[3], Area_Local = 0.0, Volume_Local = 0.0, Time_Num;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
	unsigned short iDim, iVar, jVar;
	double Density = config->GetMaterialDensity(), TimeJac = 0.0;
	
	/*--- Numerical time step (this system is uncoditional stable... a very big number can be used) ---*/
	Time_Num = config->GetDelta_UnstTimeND();
	
	/*--- Loop through elements to compute contributions from the matrix
   blocks involving time. These contributions are also added to the
   Jacobian w/ the time step. Spatial source terms are also computed. ---*/
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
    /*--- Get node numbers and their coordinate vectors ---*/
		Point_0 = geometry->elem[iElem]->GetNode(0);	Coord_0 = geometry->node[Point_0]->GetCoord();
		Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1 = geometry->node[Point_1]->GetCoord();
		Point_2 = geometry->elem[iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord();
		if (nDim == 3) { Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3 = geometry->node[Point_3]->GetCoord(); }
		
		if (nDim == 2) {
			
			for (iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}
			
			/*--- Compute element area ---*/
			Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
		}
		else {
			
			for (iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
				c[iDim] = Coord_3[iDim]-Coord_2[iDim];
			}
			d[0] = a[1]*b[2]-a[2]*b[1]; d[1] = -(a[0]*b[2]-a[2]*b[0]); d[2] = a[0]*b[1]-a[1]*b[0];
			
			/*--- Compute element volume ---*/
			Volume_Local = fabs(c[0]*d[0] + c[1]*d[1] + c[2]*d[2])/6.0;
		}
		
		/*----------------------------------------------------------------*/
		/*--- Block contributions to the Jacobian (includes time step) ---*/
		/*----------------------------------------------------------------*/
		
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				StiffMatrix_Node[iVar][jVar] = 0.0;
		
		if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST) TimeJac = 1.0/Time_Num;
		if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) TimeJac = 3.0/(2.0*Time_Num);
		
		/*--- Diagonal value identity matrix ---*/
		if (nDim == 2) {
			StiffMatrix_Node[0][0] = 1.0*TimeJac;
			StiffMatrix_Node[1][1] = 1.0*TimeJac;
		}
		else {
			StiffMatrix_Node[0][0] = 1.0*TimeJac;
			StiffMatrix_Node[1][1] = 1.0*TimeJac;
			StiffMatrix_Node[2][2] = 1.0*TimeJac;
		}
		
		/*--- Diagonal value ---*/
		if (nDim == 2) {
			StiffMatrix_Node[2][2] = Density*(2.0/12.0)*(Area_Local*TimeJac);
			StiffMatrix_Node[3][3] = Density*(2.0/12.0)*(Area_Local*TimeJac);
		}
		else {
			StiffMatrix_Node[3][3] = Density*(2.0/20.0)*(Volume_Local*TimeJac);
			StiffMatrix_Node[4][4] = Density*(2.0/20.0)*(Volume_Local*TimeJac);
			StiffMatrix_Node[5][5] = Density*(2.0/20.0)*(Volume_Local*TimeJac);
		}
    Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_1, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_2, StiffMatrix_Node);
		if (nDim == 3) { Jacobian.AddBlock(Point_3, Point_3, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_2, StiffMatrix_Node); }
		
		/*--- Off Diagonal value ---*/
		if (nDim == 2) {
			StiffMatrix_Node[2][2] = Density*(1.0/12.0)*(Area_Local*TimeJac);
			StiffMatrix_Node[3][3] = Density*(1.0/12.0)*(Area_Local*TimeJac);
		}
		else {
			StiffMatrix_Node[3][3] = Density*(1.0/20.0)*(Volume_Local*TimeJac);
			StiffMatrix_Node[4][4] = Density*(1.0/20.0)*(Volume_Local*TimeJac);
			StiffMatrix_Node[5][5] = Density*(1.0/20.0)*(Volume_Local*TimeJac);
		}
		Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_1, StiffMatrix_Node);
		Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_2, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_2, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_1, StiffMatrix_Node);
		if (nDim == 3) {
			Jacobian.AddBlock(Point_0, Point_3, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_3, StiffMatrix_Node);
			Jacobian.AddBlock(Point_1, Point_3, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_3, StiffMatrix_Node);
			Jacobian.AddBlock(Point_2, Point_3, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_3, StiffMatrix_Node);
			Jacobian.AddBlock(Point_3, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_3, Point_0, StiffMatrix_Node);
			Jacobian.AddBlock(Point_3, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_3, Point_1, StiffMatrix_Node);
			Jacobian.AddBlock(Point_3, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_3, Point_2, StiffMatrix_Node);
		}
	}
	
	unsigned long iPoint, total_index;
	double *U_time_nM1, *U_time_n, *U_time_nP1;
	
	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		
		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n   = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();
		
		/*--- Compute Residual ---*/
		for(iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				LinSysSol[total_index] = ( U_time_nP1[iVar] - U_time_n[iVar] );
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				LinSysSol[total_index] = ( U_time_nP1[iVar] - (4.0/3.0)*U_time_n[iVar] + (1.0/3.0)*U_time_nM1[iVar] );
		}
	}
	
  /*--- Contribution to the residual ---*/
	StiffMatrixTime.MatrixVectorProduct(LinSysSol, LinSysAux);
  
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			Residual[iVar] = LinSysAux[total_index];
		}
		LinSysRes.SubtractBlock(iPoint, Residual);
	}
	
}

void CFEASolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	
  unsigned short iVar;
	unsigned long iPoint, total_index;
  bool output;
	
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
	
	/*--- Solve the linear system (Krylov subspace methods) ---*/
  CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
  
  CPreconditioner* precond = NULL;
  if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
    Jacobian.BuildJacobiPreconditioner();
    precond = new CJacobiPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LU_SGS) {
    precond = new CLU_SGSPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    Jacobian.BuildJacobiPreconditioner();
    Jacobian.BuildLineletPreconditioner(geometry, config);
    precond = new CLineletPreconditioner(Jacobian, geometry, config);
  }
  
  CSysSolve system;
  
  if (config->GetUnsteady_Simulation() == STEADY) output = true;
  else output = false;
  
  if (config->GetKind_Linear_Solver() == BCGSTAB)
    system.BCGSTAB(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(), config->GetLinear_Solver_Iter(), output);
  else if (config->GetKind_Linear_Solver() == FGMRES)
    system.FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(), config->GetLinear_Solver_Iter(), output);
  
  delete mat_vec;
  delete precond;
  
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			if (config->GetUnsteady_Simulation() == STEADY) node[iPoint]->SetSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      else node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar]);
		}
	}
	
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
      AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex());
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

void CFEASolver::SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry, CGeometry **flow_geometry,
                             CConfig *fea_config, CConfig *flow_config) {
  //	unsigned short iMarker;
  //	unsigned long iVertex, iPoint;
  //	double Pressure;
  //
  //#ifdef NO_MPI
  //
  //	unsigned long iPoint_Donor;
  //
  //	for (iMarker = 0; iMarker < fea_config->GetnMarker_All(); iMarker++) {
  //		if (fea_config->GetMarker_All_Boundary(iMarker) == FLOWLOAD_BOUNDARY) {
  //			for(iVertex = 0; iVertex < fea_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
  //				iPoint = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
  //				iPoint_Donor = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetDonorPoint();
  //				Pressure = 1.0; //(flow_solution[MESH_0][FLOW_SOL]->node[iPoint_Donor]->GetPressure()-101325.0);
  //				node[iPoint]->SetPressureValue(Pressure);
  //			}
  //		}
  //	}
  //
  //#else
  //
  //	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
  //	double *Buffer_Send_U = new double [1];
  //	double *Buffer_Receive_U = new double [1];
  //	unsigned long jPoint;
  //
  //	/*--- Do the send process, by the moment we are sending each
  //	 node individually, this must be changed ---*/
  //	for (iMarker = 0; iMarker < flow_config->GetnMarker_All(); iMarker++) {
  //		/*--- There must be a better way to identify the marker that correspond with a fluid structure interation!!! ---*/
  //		if ((flow_config->GetMarker_All_Boundary(iMarker) == EULER_WALL) &&
  //        (flow_config->GetMarker_All_Moving(iMarker) == YES)) {
  //			for(iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
  //				iPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
  //
  //				if (flow_geometry[MESH_0]->node[iPoint]->GetDomain()) {
  //
  //					/*--- Find the associate pair to the original node (index and processor) ---*/
  //					jPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
  //					jProcessor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
  //
  //					/*--- We only send the pressure that belong to other boundary ---*/
  //					if (jProcessor != rank) {
  //						Buffer_Send_U[0] = 1.0; //(flow_solution[MESH_0][FLOW_SOL]->node[iPoint]->GetPressure()-101325.0);;
  //						MPI::COMM_WORLD.Bsend(Buffer_Send_U, 1, MPI::DOUBLE, jProcessor, iPoint);
  //					}
  //
  //				}
  //			}
  //		}
  //	}
  //
  //	/*--- Now the loop is over the fea points ---*/
  //	for (iMarker = 0; iMarker < fea_config->GetnMarker_All(); iMarker++) {
  //		if (fea_config->GetMarker_All_Boundary(iMarker) == LOAD_BOUNDARY) {
  //			for(iVertex = 0; iVertex < fea_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
  //				iPoint = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
  //				if (fea_geometry[MESH_0]->node[iPoint]->GetDomain()) {
  //
  //					/*--- Find the associate pair to the original node ---*/
  //					jPoint = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
  //					jProcessor = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
  //
  //					/*--- We only receive the information that belong to other boundary ---*/
  //					if (jProcessor != rank)
  //						MPI::COMM_WORLD.Recv(Buffer_Receive_U, 1, MPI::DOUBLE, jProcessor, jPoint);
  //					else
  //						Buffer_Receive_U[0] = 1.0; //(flow_solution[MESH_0][FLOW_SOL]->node[jPoint]->GetPressure()-101325.0);
  //					
  //					/*--- Store the solution for both points ---*/
  //					Pressure = Buffer_Receive_U[1];
  //					node[iPoint]->SetPressureValue(Pressure);
  //					
  //				}
  //			}
  //		}
  //	}
  //  
  //	delete[] Buffer_Send_U;
  //	delete[] Buffer_Receive_U;
  //  
  //#endif
}
