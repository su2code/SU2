/*!
 * \file solution_direct_elasticity.cpp
 * \brief Main subrotuines for solving the linear elasticity equation.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.9
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
  
  nPoint =        geometry->GetnPoint();
  nPointDomain =  geometry->GetnPointDomain();
	nDim =          geometry->GetnDim();
	nMarker =       config->GetnMarker_All();
	node =          new CVariable*[nPoint];
  if (config->GetUnsteady_Simulation() == STEADY) nVar = nDim;
  else nVar = 2*nDim;
  
	if (nDim == 2) NodesElement = 3;
	if (nDim == 3) NodesElement = 4;
	
  /*--- Define some auxiliary vectors related to the residual ---*/
  Residual = new double[nVar];          for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS = new double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max = new double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  
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
  StiffMatrixSpace.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry);
	StiffMatrixTime.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry);
  if (rank == MASTER_NODE) cout << "Initialize jacobian structure (Linear Elasticity)." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry);
  
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
	
  
  GetSurface_Pressure(geometry, config);

  
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
  
  unsigned short iVar, jVar, nNodes, iNodes, iDim, jDim;
	unsigned long iElem, PointCorners[8], iPoint, total_index;
	double CoordCorners[8][3];
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     nNodes = 3;
    if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      nNodes = 5;
    if (geometry->elem[iElem]->GetVTK_Type() == WEDGE)        nNodes = 6;
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   nNodes = 8;
    
    for (iNodes = 0; iNodes < nNodes; iNodes++) {
      PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
      for (iDim = 0; iDim < nDim; iDim++) {
        CoordCorners[iNodes][iDim] = geometry->node[PointCorners[iNodes]]->GetCoord(iDim);
      }
    }
    
    if (nDim == 2) numerics->SetFEA_StiffMatrix2D(StiffMatrix_Elem, CoordCorners, nNodes);
    else numerics->SetFEA_StiffMatrix3D(StiffMatrix_Elem, CoordCorners, nNodes);
    
    /*--- Initialization of the auxiliar matrix ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        StiffMatrix_Node[iVar][jVar] = 0.0;
    
    if (config->GetUnsteady_Simulation() == STEADY) {

      /*--- Transform the stiffness matrix into the
       contributions for the individual nodes relative to each other. ---*/
      
      for (iVar = 0; iVar < nNodes; iVar++) {
        for (jVar = 0; jVar < nNodes; jVar++) {
          for (iDim = 0; iDim < nVar; iDim++) {
            for (jDim = 0; jDim < nVar; jDim++) {
              StiffMatrix_Node[iDim][jDim] = StiffMatrix_Elem[(iVar*nVar)+iDim][(jVar*nVar)+jDim];
            }
          }
          Jacobian.AddBlock(PointCorners[iVar], PointCorners[jVar], StiffMatrix_Node);
        }
      }
    }
    
    else {
      
      if (nDim == 2) {
        StiffMatrix_Node[0][2] = -1.0;
        StiffMatrix_Node[1][3] = -1.0;
      }
      else {
        StiffMatrix_Node[0][3] = -1.0;
        StiffMatrix_Node[1][4] = -1.0;
        StiffMatrix_Node[2][5] = -1.0;
      }
      
      for (iVar = 0; iVar < nNodes; iVar++) {
        for (jVar = 0; jVar < nNodes; jVar++) {
          for (iDim = 0; iDim < nVar; iDim++) {
            for (jDim = 0; jDim < nVar; jDim++) {
              StiffMatrix_Node[nVar+iDim][jDim] = StiffMatrix_Elem[(iVar*nVar)+iDim][(jVar*nVar)+jDim];
            }
          }
          Jacobian.AddBlock(PointCorners[iVar], PointCorners[jVar], StiffMatrix_Node);
        }
      }
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          LinSysSol[total_index] = node[iPoint]->GetSolution(iVar);
          LinSysAux[total_index] = 0.0;
        }
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
  
  /*--- Deallocate memory and exit ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete StiffMatrix_Node[iVar];
  delete [] StiffMatrix_Node;
  
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
      
      Normal_Elem[0] = -(-a[1]);
			Normal_Elem[1] = -(a[0]);
      
		}
		if (nDim == 3) {
			for (iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}
			Area_Elem = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
      
      Normal_Elem[0] = -(0.5*(a[1]*b[2]-a[2]*b[1]));
			Normal_Elem[1] = -(-0.5*(a[0]*b[2]-a[2]*b[0]));
			Normal_Elem[2] = -(0.5*(a[0]*b[1]-a[1]*b[0]));
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

void CFEASolver::BC_Pressure(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                unsigned short val_marker) {
	
#ifndef DEBUG
  
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double *Coord_0 = NULL, *Coord_1 = NULL, *Coord_2 = NULL, Length_Elem = 0.0, Area_Elem = 0.0,
  Normal_Elem[3] = {0.0, 0.0, 0.0}, Pressure[3] = {0.0, 0.0, 0.0}, a[3], b[3];
	unsigned short iDim;
		
	for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
		Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);
    Coord_0 = geometry->node[Point_0]->GetCoord();
    Pressure[0] = node[Point_0]->GetFlow_Pressure();
		Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);
    Coord_1 = geometry->node[Point_1]->GetCoord();
    Pressure[1] = node[Point_1]->GetFlow_Pressure();
		if (nDim == 3) {
      Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);
      Coord_2 = geometry->node[Point_2]->GetCoord();
      Pressure[2] = node[Point_2]->GetFlow_Pressure();
    }
    
		/*--- Compute area (3D), and length of the surfaces (2D) ---*/
    
		if (nDim == 2) {
			for (iDim = 0; iDim < nDim; iDim++)
				a[iDim] = Coord_0[iDim]-Coord_1[iDim];
			Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);
      
      Normal_Elem[0] = -(-a[1]);
			Normal_Elem[1] = -(a[0]);
      
		}
		if (nDim == 3) {
			for (iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}
			Area_Elem = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
      
      Normal_Elem[0] = -(0.5*(a[1]*b[2]-a[2]*b[1]));
			Normal_Elem[1] = -(-0.5*(a[0]*b[2]-a[2]*b[0]));
			Normal_Elem[2] = -(0.5*(a[0]*b[1]-a[1]*b[0]));
		}
		
    /*--- Add the residual corresponding to the force on the surface ---*/

    if (config->GetUnsteady_Simulation() == STEADY) {
      if (nDim == 2) {
        Residual[0] = (1.0/2.0)*Pressure[0]*Normal_Elem[0];
        Residual[1] = (1.0/2.0)*Pressure[0]*Normal_Elem[1];
        LinSysRes.AddBlock(Point_0, Residual);
        Residual[0] = (1.0/2.0)*Pressure[1]*Normal_Elem[0];
        Residual[1] = (1.0/2.0)*Pressure[1]*Normal_Elem[1];
        LinSysRes.AddBlock(Point_1, Residual);
      }
      else {
        Residual[0] = (1.0/3.0)*Pressure[0]*Normal_Elem[0];
        Residual[1] = (1.0/3.0)*Pressure[0]*Normal_Elem[1];
        Residual[2] = (1.0/3.0)*Pressure[0]*Normal_Elem[2];
        LinSysRes.AddBlock(Point_0, Residual);
        
        Residual[0] = (1.0/3.0)*Pressure[1]*Normal_Elem[0];
        Residual[1] = (1.0/3.0)*Pressure[1]*Normal_Elem[1];
        Residual[2] = (1.0/3.0)*Pressure[1]*Normal_Elem[2];
        LinSysRes.AddBlock(Point_1, Residual);
        
        Residual[0] = (1.0/3.0)*Pressure[2]*Normal_Elem[0];
        Residual[1] = (1.0/3.0)*Pressure[2]*Normal_Elem[1];
        Residual[2] = (1.0/3.0)*Pressure[2]*Normal_Elem[2];
        LinSysRes.AddBlock(Point_2, Residual);
      }
    }
    else {
      if (nDim == 2) {
        Residual[0] = 0.0; Residual[1] = 0.0;
        Residual[2] = (1.0/2.0)*Pressure[0]*Normal_Elem[0];
        Residual[3] = (1.0/2.0)*Pressure[0]*Normal_Elem[1];
        LinSysRes.AddBlock(Point_0, Residual);
        Residual[0] = 0.0; Residual[1] = 0.0;
        Residual[2] = (1.0/2.0)*Pressure[1]*Normal_Elem[0];
        Residual[3] = (1.0/2.0)*Pressure[1]*Normal_Elem[1];
        LinSysRes.AddBlock(Point_1, Residual);
      }
      else {
        Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0;
        Residual[3] = (1.0/3.0)*Pressure[0]*Normal_Elem[0];
        Residual[4] = (1.0/3.0)*Pressure[0]*Normal_Elem[1];
        Residual[5] = (1.0/3.0)*Pressure[0]*Normal_Elem[2];
        LinSysRes.AddBlock(Point_0, Residual);
        
        Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0;
        Residual[3] = (1.0/3.0)*Pressure[1]*Normal_Elem[0];
        Residual[4] = (1.0/3.0)*Pressure[1]*Normal_Elem[1];
        Residual[5] = (1.0/3.0)*Pressure[1]*Normal_Elem[2];
        LinSysRes.AddBlock(Point_1, Residual);
        
        Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0;
        Residual[3] = (1.0/3.0)*Pressure[2]*Normal_Elem[0];
        Residual[4] = (1.0/3.0)*Pressure[2]*Normal_Elem[1];
        Residual[5] = (1.0/3.0)*Pressure[2]*Normal_Elem[2];
        LinSysRes.AddBlock(Point_2, Residual);
      }
    }
		
	}
  
#else

  double StiffMatrix_BoundElem[9][9], Vector_BoundElem[9], LoadNode[9];
  unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
  double Length_Elem, Area_Elem, Pressure[3];
  unsigned short iVar, jVar;
  
  if (nDim == 2) {
    
    Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);
    Pressure[0] = node[Point_0]->GetFlow_Pressure();
		Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);
    Pressure[1] = node[Point_1]->GetFlow_Pressure();
    
    StiffMatrix_BoundElem[0][0] = (2.0/6.0)*Length_Elem;
    StiffMatrix_BoundElem[0][1] = 0.0;
    StiffMatrix_BoundElem[0][2] = (1.0/6.0)*Length_Elem;
    StiffMatrix_BoundElem[0][3] = 0.0;
    StiffMatrix_BoundElem[1][0] = 0.0;
    StiffMatrix_BoundElem[1][1] = (2.0/6.0)*Length_Elem;
    StiffMatrix_BoundElem[1][2] = 0.0;
    StiffMatrix_BoundElem[1][3] = (1.0/6.0)*Length_Elem;
    StiffMatrix_BoundElem[2][0] = (1.0/6.0)*Length_Elem;
    StiffMatrix_BoundElem[2][1] = 0.0;
    StiffMatrix_BoundElem[2][2] = (2.0/6.0)*Length_Elem;
    StiffMatrix_BoundElem[2][3] = 0.0;
    StiffMatrix_BoundElem[3][0] = 0.0;
    StiffMatrix_BoundElem[3][1] = (1.0/6.0)*Length_Elem;
    StiffMatrix_BoundElem[3][2] = 0.0;
    StiffMatrix_BoundElem[3][3] = (2.0/6.0)*Length_Elem;
    
    LoadNode[0] = 0.0; LoadNode[1] = Pressure[0]; LoadNode[2] = 0.0;	LoadNode[3] = Pressure[1];
    
    for (iVar = 0; iVar < 4; iVar++) {
      Vector_BoundElem[iVar] = 0.0;
      for (jVar = 0; jVar < 4; jVar++) {
        Vector_BoundElem[iVar] += StiffMatrix_BoundElem[iVar][jVar]*LoadNode[jVar];
      }
    }
    
    Residual[0] = 0.0;
    Residual[1] = 0.0;
    Residual[2] = Vector_BoundElem[0];
    Residual[3] = Vector_BoundElem[1];
    LinSysRes.AddBlock(Point_0, Residual);
    
    Residual[0] = 0.0;
    Residual[1] = 0.0;
    Residual[2] = Vector_BoundElem[2];
    Residual[3] = Vector_BoundElem[3];
    LinSysRes.AddBlock(Point_1, Residual);
    
  }
  
  if (nDim == 3) {
    
    Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);
    Pressure[0] = node[Point_0]->GetFlow_Pressure();
		Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);
    Pressure[1] = node[Point_1]->GetFlow_Pressure();
    Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);
    Pressure[2] = node[Point_2]->GetFlow_Pressure();
    
    StiffMatrix_BoundElem[0][0] = (2.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[0][1] = 0.0;
    StiffMatrix_BoundElem[0][2] = 0.0;
    StiffMatrix_BoundElem[0][3] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[0][4] = 0.0;
    StiffMatrix_BoundElem[0][5] = 0.0;
    StiffMatrix_BoundElem[0][6] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[0][7] = 0.0;
    StiffMatrix_BoundElem[0][8] = 0.0;
    StiffMatrix_BoundElem[1][0] = 0.0;
    StiffMatrix_BoundElem[1][1] = (2.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[1][2] = 0.0;
    StiffMatrix_BoundElem[1][3] = 0.0;
    StiffMatrix_BoundElem[1][4] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[1][5] = 0.0;
    StiffMatrix_BoundElem[1][6] = 0.0;
    StiffMatrix_BoundElem[1][7] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[1][8] = 0.0;
    StiffMatrix_BoundElem[2][0] = 0.0;
    StiffMatrix_BoundElem[2][1] = 0.0;
    StiffMatrix_BoundElem[2][2] = (2.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[2][3] = 0.0;
    StiffMatrix_BoundElem[2][4] = 0.0;
    StiffMatrix_BoundElem[2][5] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[2][6] = 0.0;
    StiffMatrix_BoundElem[2][7] = 0.0;
    StiffMatrix_BoundElem[2][8] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[3][0] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[3][1] = 0.0;
    StiffMatrix_BoundElem[3][2] = 0.0;
    StiffMatrix_BoundElem[3][3] = (2.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[3][4] = 0.0;
    StiffMatrix_BoundElem[3][5] = 0.0;
    StiffMatrix_BoundElem[3][6] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[3][7] = 0.0;
    StiffMatrix_BoundElem[3][8] = 0.0;
    StiffMatrix_BoundElem[4][0] = 0.0;
    StiffMatrix_BoundElem[4][1] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[4][2] = 0.0;
    StiffMatrix_BoundElem[4][3] = 0.0;
    StiffMatrix_BoundElem[4][4] = (2.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[4][5] = 0.0;
    StiffMatrix_BoundElem[4][6] = 0.0;
    StiffMatrix_BoundElem[4][7] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[4][8] = 0.0;
    StiffMatrix_BoundElem[5][0] = 0.0;
    StiffMatrix_BoundElem[5][1] = 0.0;
    StiffMatrix_BoundElem[5][2] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[5][3] = 0.0;
    StiffMatrix_BoundElem[5][4] = 0.0;
    StiffMatrix_BoundElem[5][5] = (2.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[5][6] = 0.0;
    StiffMatrix_BoundElem[5][7] = 0.0;
    StiffMatrix_BoundElem[5][8] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[6][0] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[6][1] = 0.0;
    StiffMatrix_BoundElem[6][2] = 0.0;
    StiffMatrix_BoundElem[6][3] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[6][4] = 0.0;
    StiffMatrix_BoundElem[6][5] = 0.0;
    StiffMatrix_BoundElem[6][6] = (2.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[6][7] = 0.0;
    StiffMatrix_BoundElem[6][8] = 0.0;
    StiffMatrix_BoundElem[7][0] = 0.0;
    StiffMatrix_BoundElem[7][1] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[7][2] = 0.0;
    StiffMatrix_BoundElem[7][3] = 0.0;
    StiffMatrix_BoundElem[7][4] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[7][5] = 0.0;
    StiffMatrix_BoundElem[7][6] = 0.0;
    StiffMatrix_BoundElem[7][7] = (2.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[7][8] = 0.0;
    StiffMatrix_BoundElem[8][0] = 0.0;
    StiffMatrix_BoundElem[8][1] = 0.0;
    StiffMatrix_BoundElem[8][2] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[8][3] = 0.0;
    StiffMatrix_BoundElem[8][4] = 0.0;
    StiffMatrix_BoundElem[8][5] = (1.0/12.0)*Area_Elem;
    StiffMatrix_BoundElem[8][6] = 0.0;
    StiffMatrix_BoundElem[8][7] = 0.0;
    StiffMatrix_BoundElem[8][8] = (2.0/12.0)*Area_Elem;
    
    LoadNode[0] = 0.0;	LoadNode[1] = 0.0;	LoadNode[2] = Pressure[0];
    LoadNode[3] = 0.0;	LoadNode[4] = 0.0;	LoadNode[5] = Pressure[1];
    LoadNode[6] = 0.0;	LoadNode[7] = 0.0;	LoadNode[8] = Pressure[2];
    
    for (iVar = 0; iVar < 9; iVar++) {
      Vector_BoundElem[iVar] = 0.0;
      for (jVar = 0; jVar < 9; jVar++) {
        Vector_BoundElem[iVar] += StiffMatrix_BoundElem[iVar][jVar]*LoadNode[jVar];
      }
    }
    
    Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0;
    Residual[3] = Vector_BoundElem[0]; Residual[4] = Vector_BoundElem[1]; Residual[5] = Vector_BoundElem[2];
    LinSysRes.AddBlock(Point_0, Residual);
    
    Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0;
    Residual[3] = Vector_BoundElem[3]; Residual[4] = Vector_BoundElem[4]; Residual[5] = Vector_BoundElem[5];
    LinSysRes.AddBlock(Point_1, Residual);
    
    Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0;
    Residual[3] = Vector_BoundElem[6]; Residual[4] = Vector_BoundElem[7]; Residual[5] = Vector_BoundElem[8];
    LinSysRes.AddBlock(Point_2, Residual);
    
  }

#endif

}

void CFEASolver::BC_Flow_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                              unsigned short val_marker) { }


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

void CFEASolver::GetSurface_Pressure(CGeometry *geometry, CConfig *config) {
  
  unsigned short iMarker, icommas, iDim;
  unsigned long iVertex, iPoint, iExtIter;
  double Pressure, Dist, Coord[3];
  string text_line;
  string::size_type position;
  ifstream Surface_file;
  char buffer[50], cstr[200];
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifndef NO_MPI
#ifdef WINDOWS
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
  rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();
#endif
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
      if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
      if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
      if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
      if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
      if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
    }
    else sprintf (buffer, ".csv");
    
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
        if (config->GetMarker_All_Boundary(iMarker) == PRESSURE_BOUNDARY) {
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
                             CConfig *fea_config, CConfig *flow_config) { }
