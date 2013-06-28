/*!
 * \file solution_direct_fea.cpp
 * \brief Main subrotuines for solving the FEA equation.
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

CFEASolution::CFEASolution(void) : CSolution() { }

CFEASolution::CFEASolution(CGeometry *geometry, CConfig *config) : CSolution() {

	unsigned long nPoint, iPoint;
	unsigned short nMarker, iVar, NodesElement;
  
	nPoint  = geometry->GetnPoint();
	nDim    = geometry->GetnDim();
	nMarker = config->GetnMarker_All(); 
	node    = new CVariable*[nPoint];
	nVar    = 2*nDim;
	if (nDim == 2) NodesElement = 3;	// Trianles in 2D
	if (nDim == 3) NodesElement = 4;	// Tets in 3D
	
	Residual     = new double[nVar]; Residual_RMS = new double[nVar];
	Solution     = new double[nVar];
  Residual_Max = new double[nVar]; Point_Max = new unsigned long[nVar];
  

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
    Initialize_SparseMatrix_Structure(&StiffMatrixSpace, nVar, nVar, geometry, config);
	Initialize_SparseMatrix_Structure(&StiffMatrixTime, nVar, nVar, geometry, config);
    Initialize_SparseMatrix_Structure(&Jacobian, nVar, nVar, geometry, config);

  /*--- Initialization of linear solver structures ---*/
	xsol = new double [nPoint*nVar];
	xres = new double [nPoint*nVar];
	rhs  = new double [nPoint*nVar];

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
        if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1];
        if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2];
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

CFEASolution::~CFEASolution(void) {

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
	
	delete [] xsol;
	delete [] xres;
	delete [] rhs;

	/*--- Computation of gradients by least-squares ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
}

void CFEASolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iPoint;
	
  /*--- Set residuals to zero ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		node[iPoint]->Set_ResConv_Zero();
		node[iPoint]->Set_ResVisc_Zero();
		node[iPoint]->Set_ResSour_Zero();
	}
	
	/*--- Set matrix entries to zero ---*/
	StiffMatrixSpace.SetValZero();
	StiffMatrixTime.SetValZero();
	Jacobian.SetValZero();

}

void CFEASolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CNumerics *second_solver,
																	 CConfig *config, unsigned short iMesh) {
		
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
	double a[3], b[3], c[3], d[3], Area_Local = 0.0, Volume_Local = 0.0, Time_Num;
	double *Coord_0_ = NULL, *Coord_1_= NULL, *Coord_2_= NULL, *Coord_3_= NULL;
	double Coord_0[3], Coord_1[3], Coord_2[3], Coord_3[3];
	unsigned short iDim, iVar, jVar;
	double MassMatrix_Elem_2D [6][6] = 
	 {{ 2, 0, 1, 0, 1, 0 },
		{ 0, 2, 0, 1, 0, 1 },
		{ 1, 0, 2, 0, 1, 0 },
		{ 0, 1, 0, 2, 0, 1 },
		{ 1, 0, 1, 0, 2, 0 },
		{ 0, 1, 0, 1, 0, 2 }}; 
	
	double MassMatrix_Elem_3D [12][12] = 
	 {{ 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0 },
		{ 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0 },
		{ 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1 },
		{ 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0 },
		{ 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0 },
		{ 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1 },
		{ 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0 },
		{ 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0 },
		{ 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1 },
		{ 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0 },
		{ 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0 },
		{ 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2 }}; 
	
	double Density = config->GetMaterialDensity();
	
	/*--- Numerical time step (this system is uncoditional stable... a very big number can be used) ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING) Time_Num = config->GetDelta_UnstTimeND();
	else Time_Num = 1E30;
	if (config->GetUnsteady_Simulation() == NO) Time_Num = 0.01;
			
	/*--- Loop through elements to compute contributions from the matrix
   blocks involving time. These contributions are also added to the 
   Jacobian w/ the time step. Spatial source terms are also computed. ---*/
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
    /*--- Get node numbers and their coordinate vectors ---*/
		Point_0 = geometry->elem[iElem]->GetNode(0);	Coord_0_ = geometry->node[Point_0]->GetCoord();
		Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1_ = geometry->node[Point_1]->GetCoord();	
		Point_2 = geometry->elem[iElem]->GetNode(2);	Coord_2_ = geometry->node[Point_2]->GetCoord();
		if (nDim == 3) { Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3_ = geometry->node[Point_3]->GetCoord(); }
					
		/*--- Modification of the goemtry due to the current displacement (Value of the solution) ---*/
		for (iDim = 0; iDim < nDim; iDim++) {
			Coord_0[iDim] = Coord_0_[iDim] + node[Point_0]->GetSolution()[iDim];
			Coord_1[iDim] = Coord_1_[iDim] + node[Point_1]->GetSolution()[iDim];
			Coord_2[iDim] = Coord_2_[iDim] + node[Point_2]->GetSolution()[iDim];
			if (nDim == 3) 
				Coord_3[iDim] = Coord_3_[iDim] + node[Point_3]->GetSolution()[iDim];	
		}
		/*--- Modification of the goemtry due to the current displacement (Value of the solution) ---*/
		
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
		
		/*--- Add volumetric forces as source term (gravity and coriollis) ---*/
		double RhoG= Density*STANDART_GRAVITY;

		unsigned short NodesElement = 3; // Triangles in 2D
		if (nDim == 3) NodesElement = 4;	// Tets in 3D

		double *ElemForce = new double [nDim*NodesElement];
		double *ElemResidual = new double [nDim*NodesElement];
		
		if (nDim == 2) {
			ElemForce[0] = 0.0;	ElemForce[1] = -RhoG;	
			ElemForce[2] = 0.0;	ElemForce[3] = -RhoG; 
			ElemForce[4] = 0.0;	ElemForce[5] = -RhoG;
			
			for (iVar = 0; iVar < nDim*NodesElement; iVar++) {
				ElemResidual[iVar] = 0.0;
				for (jVar = 0; jVar < nDim*NodesElement; jVar++) {
					ElemResidual[iVar] += MassMatrix_Elem_2D[iVar][jVar]*ElemForce[jVar]*Area_Local/12.0;
				}
			}
			
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = ElemResidual[0]; Residual[3] = ElemResidual[1]; 
//			node[Point_0]->AddRes_Sour(Residual);

			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = ElemResidual[2]; Residual[3] = ElemResidual[3]; 
//			node[Point_1]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = ElemResidual[4]; Residual[3] = ElemResidual[5]; 
//			node[Point_2]->AddRes_Sour(Residual);
			
		}
		
		else {
			
			ElemForce[0] = 0.0; ElemForce[1] = 0.0;		ElemForce[2] = -RhoG;	
			ElemForce[3] = 0.0; ElemForce[4] = 0.0;		ElemForce[5] = -RhoG;	
			ElemForce[6] = 0.0; ElemForce[7] = 0.0;		ElemForce[8] = -RhoG;	
			ElemForce[9] = 0.0; ElemForce[10] = 0.0;	ElemForce[11] = -RhoG;
			
			for (iVar = 0; iVar < nDim*NodesElement; iVar++) {
				ElemResidual[iVar] = 0.0;
				for (jVar = 0; jVar < nDim*NodesElement; jVar++) {
					ElemResidual[iVar] += MassMatrix_Elem_3D[iVar][jVar]*ElemForce[jVar]*Volume_Local/20.0;
				}
			}

			Residual[0] = 0.0;							Residual[1] = 0.0;							Residual[2] = 0.0;
			Residual[3] = ElemResidual[0];	Residual[4] = ElemResidual[1];	Residual[5] = ElemResidual[2];
//			node[Point_0]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0;							Residual[1] = 0.0;							Residual[2] = 0.0;
			Residual[3] = ElemResidual[3];	Residual[4] = ElemResidual[4];	Residual[5] = ElemResidual[5]; 
//			node[Point_1]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0;							Residual[1] = 0.0;							Residual[2] = 0.0;
			Residual[3] = ElemResidual[6];	Residual[4] = ElemResidual[7];	Residual[5] = ElemResidual[8]; 
//			node[Point_2]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0;							Residual[1] = 0.0;							Residual[2] = 0.0;
			Residual[3] = ElemResidual[9];	Residual[4] = ElemResidual[10];	Residual[5] = ElemResidual[11]; 
//			node[Point_3]->AddRes_Sour(Residual);
			
		}
		
		delete [] ElemForce;
		delete [] ElemResidual;
		
		
	}

}

void CFEASolution::Source_Template(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																	 CConfig *config, unsigned short iMesh) {

}

void CFEASolution::Galerkin_Method(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh) {

	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0, iPoint, total_index;
	double *Coord_0_ = NULL, *Coord_1_= NULL, *Coord_2_= NULL, *Coord_3_= NULL;
	double Coord_0[3], Coord_1[3], Coord_2[3], Coord_3[3];
	unsigned short iVar, jVar, iDim;

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		Point_0 = geometry->elem[iElem]->GetNode(0);	Coord_0_ = geometry->node[Point_0]->GetCoord();
		Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1_ = geometry->node[Point_1]->GetCoord();	
		Point_2 = geometry->elem[iElem]->GetNode(2);	Coord_2_ = geometry->node[Point_2]->GetCoord();
		if (nDim == 3) { Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3_ = geometry->node[Point_3]->GetCoord(); }
		
		/*--- Modification of the goemtry due to the current displacement (Value of the solution) ---*/
		for (iDim = 0; iDim < nDim; iDim++) {
			Coord_0[iDim] = Coord_0_[iDim] + node[Point_0]->GetSolution()[iDim];
			Coord_1[iDim] = Coord_1_[iDim] + node[Point_1]->GetSolution()[iDim];
			Coord_2[iDim] = Coord_2_[iDim] + node[Point_2]->GetSolution()[iDim];
			if (nDim == 3) 
				Coord_3[iDim] = Coord_3_[iDim] + node[Point_3]->GetSolution()[iDim];	
		}
		/*--- Modification of the goemtry due to the current displacement (Value of the solution) ---*/
		
		if (nDim == 2) solver->SetCoord(Coord_0, Coord_1, Coord_2);
		if (nDim == 3) solver->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
		
		solver->SetResidual(StiffMatrix_Elem, config);
    
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				StiffMatrix_Node[iVar][jVar] = 0.0;	
		
		if (nDim == 2) {
			StiffMatrix_Node[0][2] = -1.0;	
			StiffMatrix_Node[1][3] = -1.0;
		}
		else {
			StiffMatrix_Node[0][3] = -1.0;	
			StiffMatrix_Node[1][4] = -1.0;	
			StiffMatrix_Node[2][5] = -1.0;
		}
		
		if (nDim == 2) {
			StiffMatrix_Node[2][0] = StiffMatrix_Elem[0][0];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[0][1];
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[1][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[1][1];
		}
		else {
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[0][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[0][1];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[0][2];
			StiffMatrix_Node[4][0] = StiffMatrix_Elem[1][0];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[1][1];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[1][2];
			StiffMatrix_Node[5][0] = StiffMatrix_Elem[2][0];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[2][1];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[2][2];			
		}
		StiffMatrixSpace.AddBlock(Point_0, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
    
		if (nDim == 2) {
			StiffMatrix_Node[2][0] = StiffMatrix_Elem[0][2];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[0][3];
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[1][2];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[1][3];
		}
		else {
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[0][3];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[0][4];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[0][5];
			StiffMatrix_Node[4][0] = StiffMatrix_Elem[1][3];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[1][4];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[1][5];
			StiffMatrix_Node[5][0] = StiffMatrix_Elem[2][3];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[2][4];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[2][5];			
		}
    StiffMatrixSpace.AddBlock(Point_0, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
    
		if (nDim == 2) {
			StiffMatrix_Node[2][0] = StiffMatrix_Elem[0][4];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[0][5];
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[1][4];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[1][5];
		}
		else {
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[0][6];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[0][7];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[0][8];
			StiffMatrix_Node[4][0] = StiffMatrix_Elem[1][6];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[1][7];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[1][8];
			StiffMatrix_Node[5][0] = StiffMatrix_Elem[2][6];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[2][7];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[2][8];			
		}
    StiffMatrixSpace.AddBlock(Point_0, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
    
		if (nDim == 3) {
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[0][9];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[0][10];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[0][11];
			StiffMatrix_Node[4][0] = StiffMatrix_Elem[1][9];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[1][10];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[1][11];
			StiffMatrix_Node[5][0] = StiffMatrix_Elem[2][9];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[2][10];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[2][11];
			StiffMatrixSpace.AddBlock(Point_0, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_3, StiffMatrix_Node);
		}
		
		if (nDim == 2) {
			StiffMatrix_Node[2][0] = StiffMatrix_Elem[2][0];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[2][1];
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][1];
		}
		else {
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][1];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[3][2];
			StiffMatrix_Node[4][0] = StiffMatrix_Elem[4][0];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[4][1];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[4][2];
			StiffMatrix_Node[5][0] = StiffMatrix_Elem[5][0];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[5][1];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[5][2];			
		}
    StiffMatrixSpace.AddBlock(Point_1, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
    
		if (nDim == 2) {
			StiffMatrix_Node[2][0] = StiffMatrix_Elem[2][2];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[2][3];
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][2];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][3];
		}
		else {
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][3];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][4];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[3][5];
			StiffMatrix_Node[4][0] = StiffMatrix_Elem[4][3];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[4][4];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[4][5];
			StiffMatrix_Node[5][0] = StiffMatrix_Elem[5][3];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[5][4];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[5][5];			
		}
    StiffMatrixSpace.AddBlock(Point_1, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
    
		if (nDim == 2) {
			StiffMatrix_Node[2][0] = StiffMatrix_Elem[2][4];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[2][5];
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][4];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][5];
		}
		else {
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][6];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][7];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[3][8];
			StiffMatrix_Node[4][0] = StiffMatrix_Elem[4][6];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[4][7];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[4][8];
			StiffMatrix_Node[5][0] = StiffMatrix_Elem[5][6];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[5][7];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[5][8];			
		}
    StiffMatrixSpace.AddBlock(Point_1, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
    
		if (nDim == 3) {
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[3][9];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[3][10];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[3][11];
			StiffMatrix_Node[4][0] = StiffMatrix_Elem[4][9];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[4][10];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[4][11];
			StiffMatrix_Node[5][0] = StiffMatrix_Elem[5][9];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[5][10];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[5][11];
			StiffMatrixSpace.AddBlock(Point_1, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_3, StiffMatrix_Node);
		}
		
		if (nDim == 2) {
			StiffMatrix_Node[2][0] = StiffMatrix_Elem[4][0];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[4][1];
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[5][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[5][1];
		}
		else {
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[6][0];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[6][1];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[6][2];
			StiffMatrix_Node[4][0] = StiffMatrix_Elem[7][0];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[7][1];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[7][2];
			StiffMatrix_Node[5][0] = StiffMatrix_Elem[8][0];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[8][1];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[8][2];			
		}
    StiffMatrixSpace.AddBlock(Point_2, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
    
		if (nDim == 2) {
			StiffMatrix_Node[2][0] = StiffMatrix_Elem[4][2];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[4][3];
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[5][2];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[5][3];
		}
		else {
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[6][3];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[6][4];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[6][5];
			StiffMatrix_Node[4][0] = StiffMatrix_Elem[7][3];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[7][4];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[7][5];
			StiffMatrix_Node[5][0] = StiffMatrix_Elem[8][3];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[8][4];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[8][5];			
		}
    StiffMatrixSpace.AddBlock(Point_2, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
    
		if (nDim == 2) {
			StiffMatrix_Node[2][0] = StiffMatrix_Elem[4][4];	StiffMatrix_Node[2][1] = StiffMatrix_Elem[4][5];	StiffMatrix_Node[2][2] = 0.0;		StiffMatrix_Node[2][3] = 0.0;
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[5][4];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[5][5];	StiffMatrix_Node[3][2] = 0.0;		StiffMatrix_Node[3][3] = 0.0;
		}
		else {
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[6][6];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[6][7];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[6][8];
			StiffMatrix_Node[4][0] = StiffMatrix_Elem[7][6];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[7][7];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[7][8];
			StiffMatrix_Node[5][0] = StiffMatrix_Elem[8][6];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[8][7];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[8][8];			
		}
		StiffMatrixSpace.AddBlock(Point_2, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
		
		if (nDim == 3) {
			StiffMatrix_Node[3][0] = StiffMatrix_Elem[6][9];	StiffMatrix_Node[3][1] = StiffMatrix_Elem[6][10];	StiffMatrix_Node[3][2] = StiffMatrix_Elem[6][11];
			StiffMatrix_Node[4][0] = StiffMatrix_Elem[7][9];	StiffMatrix_Node[4][1] = StiffMatrix_Elem[7][10];	StiffMatrix_Node[4][2] = StiffMatrix_Elem[7][11];
			StiffMatrix_Node[5][0] = StiffMatrix_Elem[8][9];	StiffMatrix_Node[5][1] = StiffMatrix_Elem[8][10];	StiffMatrix_Node[5][2] = StiffMatrix_Elem[8][11];
			StiffMatrixSpace.AddBlock(Point_2, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_3, StiffMatrix_Node);
		}
		
		if (nDim == 3) {
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
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			xsol[total_index] = node[iPoint]->GetSolution(iVar);
			xres[total_index] = 0.0;
		}
	
	StiffMatrixSpace.MatrixVectorProduct(xsol, xres);
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			Residual[iVar] = xres[total_index];
		}
		node[iPoint]->SubtractRes_Visc(Residual);
	} 
}

/*void CFEASolution::BC_Displacement(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
																	 unsigned short val_marker) {
	unsigned long iPoint, iVertex, total_index;
	unsigned short iVar;
	
	double TotalDispl = config->GetDispl_Value(config->GetMarker_All_Tag(val_marker));
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		if (nDim == 2) {
			Solution[0] = TotalDispl;		Solution[1] = 0.0;	Solution[2] = 0.0;	Solution[3] = 0.0;
			Residual[0] = 0.0;					Residual[1] = 0.0;	Residual[2] = 0.0;	Residual[3] = 0.0;
		}
		else {
			Solution[0] = 0.0; Solution[1] = TotalDispl;	Solution[2] = 0.0;	Solution[3] = 0.0;	Solution[4] = 0.0;	Solution[5] = 0.0;
			Residual[0] = 0.0; Residual[1] = 0.0;					Residual[2] = 0.0;	Residual[3] = 0.0;	Residual[4] = 0.0;	Residual[5] = 0.0;
		}
		
		node[iPoint]->SetSolution(Solution);
		node[iPoint]->SetSolution_Old(Solution);
    node[iPoint]->SetRes_Visc(Residual); 
		node[iPoint]->SetRes_Sour(Residual);
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      Jacobian.DeleteValsRowi(total_index);
    }
	}
}*/

void CFEASolution::BC_Displacement(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
																	 unsigned short val_marker) {
	unsigned long iVertex, total_index;
	unsigned short iVar, iDim, nDim; 
	unsigned long iPoint;
	double r[3], rotCoord[3], rotMatrix[3][3], vel[3];
  double *Coord_, Coord[3], *Center, *Omega, Lref, dt;
	double dtheta, dphi, dpsi, cosTheta, sinTheta;
  double cosPhi, sinPhi, cosPsi, sinPsi;
	
	/*--- Problem dimension and physical time step ---*/
	nDim = geometry->GetnDim();
	dt   = config->GetDelta_UnstTimeND();
  
	/*--- Center of rotation & angular velocity vector from config ---*/
	Center = config->GetRotAxisOrigin();
	Omega  = config->GetOmega_FreeStreamND();
	Lref   = config->GetLength_Ref();
  
	/*--- Compute delta change in the angle about the x, y, & z axes. ---*/
	dtheta = Omega[0]*dt; dphi   = Omega[1]*dt; dpsi   = Omega[2]*dt;
  
	/*--- Store angles separately for clarity. Compute sines/cosines. ---*/
	cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
	sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
  
	/*--- Compute the rotation matrix. Note that the implicit
   ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
	rotMatrix[0][0] = cosPhi*cosPsi;
	rotMatrix[1][0] = cosPhi*sinPsi;
	rotMatrix[2][0] = -sinPhi;
  
	rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
	rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
	rotMatrix[2][1] = sinTheta*cosPhi;
  
	rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
	rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
	rotMatrix[2][2] = cosTheta*cosPhi;
	
	/*--- Loop over and rotate each node in the marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Coordinates of the current point ---*/
    Coord_ = geometry->node[iPoint]->GetCoord();
		
		/*--- If we use dual time stepping the displacement if fixed during the inner iteration ---*/
		if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || 
				(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {		
			for (iDim = 0; iDim < nDim; iDim++)
				Coord[iDim] = Coord_[iDim] + node[iPoint]->GetSolution_time_n()[iDim];
		}
		else {
			for (iDim = 0; iDim < nDim; iDim++)
				Coord[iDim] = Coord_[iDim] + node[iPoint]->GetSolution()[iDim];
		}
		
    /*--- Calculate non-dim. position from rotation center ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
    if (nDim == 2) r[nDim] = 0.0;
    
    /*--- Compute transformed point coordinates ---*/
    rotCoord[0] = rotMatrix[0][0]*r[0] + rotMatrix[0][1]*r[1] + rotMatrix[0][2]*r[2] + Center[0];
    rotCoord[1] = rotMatrix[1][0]*r[0] + rotMatrix[1][1]*r[1] + rotMatrix[1][2]*r[2] + Center[1];
    rotCoord[2] = rotMatrix[2][0]*r[0] + rotMatrix[2][1]*r[1] + rotMatrix[2][2]*r[2] + Center[2];
		
		/*--- Compute the linear velocity v = Omega*r ---*/
		vel[0] = Omega[1]*r[2]-Omega[2]*r[1]; 
		vel[1] = -(Omega[0]*r[2]-Omega[2]*r[0]); 
		vel[2] = Omega[0]*r[1]-Omega[1]*r[0];
		
		
		/*--- Note that the displacement is computed with respect the original location ---*/
		if (nDim == 2) {
			Solution[0] = rotCoord[0]-Coord_[0];	Solution[1] = rotCoord[1]-Coord_[1];	
			Solution[2] = vel[0];	Solution[3] = vel[1];
			Residual[0] = 0.0;	Residual[1] = 0.0;
			Residual[2] = 0.0;	Residual[3] = 0.0;
		}
		else {
			Solution[0] = rotCoord[0]-Coord_[0];	Solution[1] = rotCoord[1]-Coord_[1];	Solution[2] = rotCoord[2]-Coord_[2];	
			Solution[3] = vel[0];	Solution[4] = vel[1];	Solution[5] = vel[2];
			Residual[0] = 0.0;	Residual[1] = 0.0;	Residual[2] = 0.0;										
			Residual[3] = 0.0;	Residual[4] = 0.0;	Residual[5] = 0.0;
		}
				
		
#ifdef CHECK
		
		double Cyclic_Pitch, Cyclic_Omega;
    double Cyclic_Mag, VarCoord[3];
    double Cyclic_Origin_New[3], Cyclic_Axis_New[3];
    double Time_New, Time_Old, Phase_Lag, Alpha_New, Alpha_Old, dalpha;
    double DEG2RAD = PI_NUMBER/180.0;
    unsigned long iter, iVertex;
    unsigned short iMarker;
    string Marker_Tag;
		
		/*--- Cyclic pitch information from config. Convert degrees to
		 radians and note that our pitching frequency should match the 
		 rotation frequency of the rotor, i.e. time of one revolution 
		 equals one cyclic pitch period. ---*/
		Cyclic_Pitch  = config->GetCyclic_Pitch()*DEG2RAD;
		Cyclic_Omega  = config->GetOmegaMag();
		
		/*--- Compute physical time based on iteration number ---*/
		iter     = config->GetExtIter();
		Time_New = static_cast<double>(iter)*dt;
		Time_Old = Time_New;
		if (iter != 0) Time_Old = (static_cast<double>(iter)-1.0)*dt;
		
		/*--- Retrieve user specified cyclic pitch origin & axis ---*/
		double Cyclic_Offset = 0.0*DEG2RAD;
		double Cyclic_Origin[3] = {0.0, 0.0675, 0.0081};
		double Cyclic_Axis[3]   = {1.0, 0.0, 0.0};
		
		/*--- Compute cyclic pitch angle delta at this time step. ---*/
		Phase_Lag = PI_NUMBER/2.0 + Cyclic_Offset;
		Alpha_New = Cyclic_Pitch*cos(Cyclic_Omega*Time_New - Phase_Lag);
		Alpha_Old = Cyclic_Pitch*cos(Cyclic_Omega*Time_Old - Phase_Lag);
		dalpha    = (1E-10 + (Alpha_New - Alpha_Old));
		
		/*--- Normalize the cyclic pitching axis, just in case ---*/		
		Cyclic_Mag = 0.0;
		for (iDim = 0; iDim < 3; iDim++)
			Cyclic_Mag += Cyclic_Axis[iDim]*Cyclic_Axis[iDim];
		Cyclic_Mag = sqrt(Cyclic_Mag);
		for (iDim = 0; iDim < 3; iDim++)
			Cyclic_Axis[iDim] = Cyclic_Axis[iDim]/Cyclic_Mag;
		
		/*--- Get total rotation since start of simulation (t = 0) ---*/
		dtheta = Omega[0]*Time_New; dphi   = Omega[1]*Time_New; dpsi   = Omega[2]*Time_New;
		
		/*--- Store angles separately for clarity. Compute sines/cosines. ---*/
		cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
		sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
		
		/*--- Compute the rotation matrix. Note that the implicit
		 ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
		rotMatrix[0][0] = cosPhi*cosPsi;
		rotMatrix[1][0] = cosPhi*sinPsi;
		rotMatrix[2][0] = -sinPhi;
		
		rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
		rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
		rotMatrix[2][1] = sinTheta*cosPhi;
		
		rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
		rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
		rotMatrix[2][2] = cosTheta*cosPhi;
		
		/*--- Transform origin and axis with the rotation matrix. ---*/
		Cyclic_Origin_New[0] = rotMatrix[0][0]*Cyclic_Origin[0] + rotMatrix[0][1]*Cyclic_Origin[1] + rotMatrix[0][2]*Cyclic_Origin[2];
		Cyclic_Origin_New[1] = rotMatrix[1][0]*Cyclic_Origin[0] + rotMatrix[1][1]*Cyclic_Origin[1] + rotMatrix[1][2]*Cyclic_Origin[2];
		Cyclic_Origin_New[2] = rotMatrix[2][0]*Cyclic_Origin[0] + rotMatrix[2][1]*Cyclic_Origin[1] + rotMatrix[2][2]*Cyclic_Origin[2];
		Cyclic_Axis_New[0]   = rotMatrix[0][0]*Cyclic_Axis[0] + rotMatrix[0][1]*Cyclic_Axis[1] + rotMatrix[0][2]*Cyclic_Axis[2];
		Cyclic_Axis_New[1]   = rotMatrix[1][0]*Cyclic_Axis[0] + rotMatrix[1][1]*Cyclic_Axis[1] + rotMatrix[1][2]*Cyclic_Axis[2];
		Cyclic_Axis_New[2]   = rotMatrix[2][0]*Cyclic_Axis[0] + rotMatrix[2][1]*Cyclic_Axis[1] + rotMatrix[2][2]*Cyclic_Axis[2];
		
		/*--- Project pitching delta onto cyclic axis ---*/		
		dtheta = Cyclic_Axis_New[0]*dalpha;  
		dphi   = Cyclic_Axis_New[1]*dalpha; 
		dpsi   = Cyclic_Axis_New[2]*dalpha;
		
		/*--- Store angles separately for clarity. Compute sines/cosines. ---*/			
		cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
		sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
		
		/*--- Compute the rotation matrix. Note that the implicit
		 ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
		rotMatrix[0][0] = cosPhi*cosPsi;
		rotMatrix[1][0] = cosPhi*sinPsi;
		rotMatrix[2][0] = -sinPhi;
		
		rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
		rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
		rotMatrix[2][1] = sinTheta*cosPhi;
		
		rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
		rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
		rotMatrix[2][2] = cosTheta*cosPhi;
		
		/*--- Coordinates of the current point ---*/
    Coord_ = geometry->node[iPoint]->GetCoord();
		
		/*--- If we use dual time stepping the displacement if fixed during the inner iteration ---*/
		if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || 
				(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {		
			for (iDim = 0; iDim < nDim; iDim++)
				Coord[iDim] = Coord_[iDim] + node[iPoint]->GetSolution_time_n()[iDim];
		}
		else {
			for (iDim = 0; iDim < nDim; iDim++)
				Coord[iDim] = Coord_[iDim] + node[iPoint]->GetSolution()[iDim];
		}
		
//		/*--- Coordinates of the current point, note that we are adding the previous movement ---*/
//		for (iDim = 0; iDim < nDim; iDim++)
//			Coord[iDim] = Coord_[iDim] + Solution[iDim];
		
		/*--- Calculate non-dim. position from rotation center ---*/
		for (iDim = 0; iDim < nDim; iDim++)
			r[iDim] = (Coord[iDim]-Cyclic_Origin_New[iDim])/Lref;
		if (nDim == 2) r[nDim] = 0.0;
		
		/*--- Compute transformed point coordinates. ---*/
		rotCoord[0] = rotMatrix[0][0]*r[0] + rotMatrix[0][1]*r[1] + rotMatrix[0][2]*r[2] + Cyclic_Origin_New[0];
		rotCoord[1] = rotMatrix[1][0]*r[0] + rotMatrix[1][1]*r[1] + rotMatrix[1][2]*r[2] + Cyclic_Origin_New[1];
		rotCoord[2] = rotMatrix[2][0]*r[0]  + rotMatrix[2][1]*r[1] + rotMatrix[2][2]*r[2] + Cyclic_Origin_New[2];
		
		/*--- Note that the displacement is computed with respect the original location ---*/
		if (nDim == 2) {
			Solution[0] = rotCoord[0]-Coord_[0];	Solution[1] = rotCoord[1]-Coord_[1];	
			Solution[2] = vel[0];	Solution[3] = vel[1];
			Residual[0] = 0.0;	Residual[1] = 0.0;
			Residual[2] = 0.0;	Residual[3] = 0.0;
		}
		else {
			Solution[0] = rotCoord[0]-Coord_[0];	Solution[1] = rotCoord[1]-Coord_[1];	Solution[2] = rotCoord[2]-Coord_[2];	
			Solution[3] = vel[0];	Solution[4] = vel[1];	Solution[5] = vel[2];
			Residual[0] = 0.0;	Residual[1] = 0.0;	Residual[2] = 0.0;										
			Residual[3] = 0.0;	Residual[4] = 0.0;	Residual[5] = 0.0;
		}
#endif
		
		node[iPoint]->SetSolution(Solution);
		node[iPoint]->SetSolution_Old(Solution);
    node[iPoint]->SetRes_Visc(Residual); 
		node[iPoint]->SetRes_Sour(Residual);
		
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      Jacobian.DeleteValsRowi(total_index);
    }
		
	}
	

}

void CFEASolution::BC_FlowLoad(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
													 unsigned short val_marker) {
	
	double a[3], b[3], Press_0 = 0.0, Press_1 = 0.0, Press_2 = 0.0, Press_Elem;
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double *Coord_0_ = NULL, *Coord_1_= NULL, *Coord_2_= NULL;
	double Coord_0[3], Coord_1[3], Coord_2[3];
	double Length_Elem, Area_Elem, Normal_Elem[3] = {0.0,0.0,0.0};
	unsigned short iDim;
		
	for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {		
		Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);	Coord_0_ = geometry->node[Point_0]->GetCoord(); Press_0 = 1.0; //node[Point_0]->GetPressure();
		Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);	Coord_1_ = geometry->node[Point_1]->GetCoord(); Press_1 = 1.0; //node[Point_1]->GetPressure();
		if (nDim == 3) { Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);	Coord_2_ = geometry->node[Point_2]->GetCoord(); Press_2 = 1.0; }//node[Point_2]->GetPressure(); }
		
		/*--- Modification of the goemtry due to the current displacement (Value of the solution) ---*/
		for (iDim = 0; iDim < nDim; iDim++) {
			Coord_0[iDim] = Coord_0_[iDim] + node[Point_0]->GetSolution()[iDim];
			Coord_1[iDim] = Coord_1_[iDim] + node[Point_1]->GetSolution()[iDim];
			if (nDim == 3) 
				Coord_2[iDim] = Coord_2_[iDim] + node[Point_2]->GetSolution()[iDim];			
		}
		/*--- Modification of the goemtry due to the current displacement (Value of the solution) ---*/		
		
		/*--- Compute area (3D), and length of the surfaces (2D) ---*/
		if (nDim == 2) {
			for (iDim = 0; iDim < nDim; iDim++)
				a[iDim] = Coord_0[iDim]-Coord_1[iDim];
			Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);
			
			Normal_Elem[0] = -a[1];
			Normal_Elem[1] = a[0];

			Press_Elem = 0.5*(Press_0 + Press_1);
		}
		else {
			
			for (iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}
			Area_Elem = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
			
			Normal_Elem[0] = 0.5*(a[1]*b[2]-a[2]*b[1]);
			Normal_Elem[1] = -0.5*(a[0]*b[2]-a[2]*b[0]);
			Normal_Elem[2] = 0.5*(a[0]*b[1]-a[1]*b[0]);
			
			Press_Elem = (Press_0 + Press_1 + Press_2)/3.0;
		}	
		
		if (nDim == 2) {
			Residual[0] = 0.0; Residual[1] = 0.0; 
			Residual[2] = (1.0/2.0)*Press_Elem*Normal_Elem[0]; Residual[3] = (1.0/2.0)*Press_Elem*Normal_Elem[1];
			node[Point_0]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0; Residual[1] = 0.0; 
			Residual[2] = (1.0/2.0)*Press_Elem*Normal_Elem[0]; Residual[3] = (1.0/2.0)*Press_Elem*Normal_Elem[1];
			node[Point_1]->AddRes_Sour(Residual);
		}
		else {
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0; 
			Residual[3] = (1.0/3.0)*Press_Elem*Normal_Elem[0]; Residual[4] = (1.0/3.0)*Press_Elem*Normal_Elem[1]; Residual[5] = (1.0/3.0)*Press_Elem*Normal_Elem[2];
			node[Point_0]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0; 
			Residual[3] = (1.0/3.0)*Press_Elem*Normal_Elem[0]; Residual[4] = (1.0/3.0)*Press_Elem*Normal_Elem[1]; Residual[5] = (1.0/3.0)*Press_Elem*Normal_Elem[2];
			node[Point_1]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0; 
			Residual[3] = (1.0/3.0)*Press_Elem*Normal_Elem[0]; Residual[4] = (1.0/3.0)*Press_Elem*Normal_Elem[1]; Residual[5] = (1.0/3.0)*Press_Elem*Normal_Elem[2];
			node[Point_2]->AddRes_Sour(Residual);
		}
		
	/*
	 
	 double StiffMatrix_BoundElem[9][9], Vector_BoundElem[9], a[3], b[3], LoadNode[9], Press_0, Press_1, Press_2, Press_Elem;
	 unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Vertex_0, Vertex_1, Vertex_2;
	 double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, Length_Elem, Area_Elem, *Normal, Area_0, Area_1, Area_2, Normal_Elem[3];
	 unsigned short iVar, jVar, iDim;
	 
	 double LocalLoad = config->GetLoad_Value(config->GetMarker_All_Tag(val_marker));
	 
		if (nDim == 2) {
			
			StiffMatrix_BoundElem[0][0] = (2.0/6.0)*Length_Elem;	StiffMatrix_BoundElem[0][1] = 0.0;										StiffMatrix_BoundElem[0][2] = (1.0/6.0)*Length_Elem;			StiffMatrix_BoundElem[0][3] = 0.0;
			StiffMatrix_BoundElem[1][0] = 0.0;										StiffMatrix_BoundElem[1][1] = (2.0/6.0)*Length_Elem;	StiffMatrix_BoundElem[1][2] = 0.0;												StiffMatrix_BoundElem[1][3] = (1.0/6.0)*Length_Elem;
			StiffMatrix_BoundElem[2][0] = (1.0/6.0)*Length_Elem;	StiffMatrix_BoundElem[2][1] = 0.0;										StiffMatrix_BoundElem[2][2] = (2.0/6.0)*Length_Elem;			StiffMatrix_BoundElem[2][3] = 0.0;
			StiffMatrix_BoundElem[3][0] = 0.0;										StiffMatrix_BoundElem[3][1] = (1.0/6.0)*Length_Elem;	StiffMatrix_BoundElem[3][2] = 0.0;												StiffMatrix_BoundElem[3][3] = (2.0/6.0)*Length_Elem;
			
			LoadNode[0] = 0.0; LoadNode[1] = LocalLoad; LoadNode[2] = 0.0;	LoadNode[3] = LocalLoad;
			
			for (iVar = 0; iVar < 4; iVar++) {
				Vector_BoundElem[iVar] = 0.0;
				for (jVar = 0; jVar < 4; jVar++) {
					Vector_BoundElem[iVar] += StiffMatrix_BoundElem[iVar][jVar]*LoadNode[jVar];
				}
			}
			
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = Vector_BoundElem[0]; Residual[3] = Vector_BoundElem[1];
			node[Point_0]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = Vector_BoundElem[2]; Residual[3] = Vector_BoundElem[3];
			node[Point_1]->AddRes_Sour(Residual);
			
		}

		if (nDim == 3) {
			
			StiffMatrix_BoundElem[0][0] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[0][1] = 0.0;										StiffMatrix_BoundElem[0][2] = 0.0;										StiffMatrix_BoundElem[0][3] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[0][4] = 0.0;										StiffMatrix_BoundElem[0][5] = 0.0;										StiffMatrix_BoundElem[0][6] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[0][7] = 0.0;										StiffMatrix_BoundElem[0][8] = 0.0;
			StiffMatrix_BoundElem[1][0] = 0.0;										StiffMatrix_BoundElem[1][1] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[1][2] = 0.0;										StiffMatrix_BoundElem[1][3] = 0.0;										StiffMatrix_BoundElem[1][4] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[1][5] = 0.0;										StiffMatrix_BoundElem[1][6] = 0.0;										StiffMatrix_BoundElem[1][7] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[1][8] = 0.0;
			StiffMatrix_BoundElem[2][0] = 0.0;										StiffMatrix_BoundElem[2][1] = 0.0;										StiffMatrix_BoundElem[2][2] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[2][3] = 0.0;										StiffMatrix_BoundElem[2][4] = 0.0;										StiffMatrix_BoundElem[2][5] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[2][6] = 0.0;										StiffMatrix_BoundElem[2][7] = 0.0;										StiffMatrix_BoundElem[2][8] = (1.0/12.0)*Area_Elem;
			StiffMatrix_BoundElem[3][0] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[3][1] = 0.0;										StiffMatrix_BoundElem[3][2] = 0.0;										StiffMatrix_BoundElem[3][3] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[3][4] = 0.0;										StiffMatrix_BoundElem[3][5] = 0.0;										StiffMatrix_BoundElem[3][6] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[3][7] = 0.0;										StiffMatrix_BoundElem[3][8] = 0.0;
			StiffMatrix_BoundElem[4][0] = 0.0;										StiffMatrix_BoundElem[4][1] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[4][2] = 0.0;										StiffMatrix_BoundElem[4][3] = 0.0;										StiffMatrix_BoundElem[4][4] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[4][5] = 0.0;										StiffMatrix_BoundElem[4][6] = 0.0;										StiffMatrix_BoundElem[4][7] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[4][8] = 0.0;
			StiffMatrix_BoundElem[5][0] = 0.0;										StiffMatrix_BoundElem[5][1] = 0.0;										StiffMatrix_BoundElem[5][2] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[5][3] = 0.0;										StiffMatrix_BoundElem[5][4] = 0.0;										StiffMatrix_BoundElem[5][5] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[5][6] = 0.0;										StiffMatrix_BoundElem[5][7] = 0.0;										StiffMatrix_BoundElem[5][8] = (1.0/12.0)*Area_Elem;
			StiffMatrix_BoundElem[6][0] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[6][1] = 0.0;										StiffMatrix_BoundElem[6][2] = 0.0;										StiffMatrix_BoundElem[6][3] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[6][4] = 0.0;										StiffMatrix_BoundElem[6][5] = 0.0;										StiffMatrix_BoundElem[6][6] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[6][7] = 0.0;										StiffMatrix_BoundElem[6][8] = 0.0;
			StiffMatrix_BoundElem[7][0] = 0.0;										StiffMatrix_BoundElem[7][1] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[7][2] = 0.0;										StiffMatrix_BoundElem[7][3] = 0.0;										StiffMatrix_BoundElem[7][4] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[7][5] = 0.0;										StiffMatrix_BoundElem[7][6] = 0.0;										StiffMatrix_BoundElem[7][7] = (2.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[7][8] = 0.0;
			StiffMatrix_BoundElem[8][0] = 0.0;										StiffMatrix_BoundElem[8][1] = 0.0;										StiffMatrix_BoundElem[8][2] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[8][3] = 0.0;										StiffMatrix_BoundElem[8][4] = 0.0;										StiffMatrix_BoundElem[8][5] = (1.0/12.0)*Area_Elem;		StiffMatrix_BoundElem[8][6] = 0.0;										StiffMatrix_BoundElem[8][7] = 0.0;										StiffMatrix_BoundElem[8][8] = (2.0/12.0)*Area_Elem;
			
			LoadNode[0] = 0.0;	LoadNode[1] = 0.0;	LoadNode[2] = LocalLoad;
			LoadNode[3] = 0.0;	LoadNode[4] = 0.0;	LoadNode[5] = LocalLoad;
			LoadNode[6] = 0.0;	LoadNode[7] = 0.0;	LoadNode[8] = LocalLoad;
			
			for (iVar = 0; iVar < 9; iVar++) {
				Vector_BoundElem[iVar] = 0.0;
				for (jVar = 0; jVar < 9; jVar++) {
					Vector_BoundElem[iVar] += StiffMatrix_BoundElem[iVar][jVar]*LoadNode[jVar];
				}
			}
			
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0; Residual[3] = Vector_BoundElem[0]; Residual[4] = Vector_BoundElem[1]; Residual[5] = Vector_BoundElem[2];
			node[Point_0]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0; Residual[3] = Vector_BoundElem[3]; Residual[4] = Vector_BoundElem[4]; Residual[5] = Vector_BoundElem[5];
			node[Point_1]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0; Residual[3] = Vector_BoundElem[6]; Residual[4] = Vector_BoundElem[7]; Residual[5] = Vector_BoundElem[8];
			node[Point_2]->AddRes_Sour(Residual);
						
		}
	 */
		
	}
	
}

void CFEASolution::BC_Load(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
															 unsigned short val_marker) {
	
	double a[3], b[3];
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double *Coord_0_ = NULL, *Coord_1_= NULL, *Coord_2_= NULL;
	double Coord_0[3], Coord_1[3], Coord_2[3];
	double Length_Elem = 0.0, Area_Elem = 0.0;
	unsigned short iDim;
	
	double TotalLoad = config->GetLoad_Value(config->GetMarker_All_Tag(val_marker));
	
	for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {		
		Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);	Coord_0_ = geometry->node[Point_0]->GetCoord();
		Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);	Coord_1_ = geometry->node[Point_1]->GetCoord();
		if (nDim == 3) { Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);	Coord_2_ = geometry->node[Point_2]->GetCoord();}
		
		/*--- Modification of the goemtry due to the current displacement (Value of the solution) ---*/
		for (iDim = 0; iDim < nDim; iDim++) {
			Coord_0[iDim] = Coord_0_[iDim] + node[Point_0]->GetSolution()[iDim];
			Coord_1[iDim] = Coord_1_[iDim] + node[Point_1]->GetSolution()[iDim];
			if (nDim == 3) 
				Coord_2[iDim] = Coord_2_[iDim] + node[Point_2]->GetSolution()[iDim];			
		}
		/*--- Modification of the goemtry due to the current displacement (Value of the solution) ---*/		
		
		/*--- Compute area (3D), and length of the surfaces (2D) ---*/
		if (nDim == 2) {
			for (iDim = 0; iDim < nDim; iDim++)
				a[iDim] = Coord_0[iDim]-Coord_1[iDim];
			Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);			
		}
		if (nDim == 3) {
			for (iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}
			Area_Elem = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
		}	
		
		if (nDim == 2) {
			Residual[0] = 0.0; Residual[1] = 0.0; 
			Residual[2] = (1.0/2.0)*TotalLoad*Length_Elem; Residual[3] = 0.0;
			node[Point_0]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0; Residual[1] = 0.0; 
			Residual[2] = (1.0/2.0)*TotalLoad*Length_Elem; Residual[3] = 0.0;
			node[Point_1]->AddRes_Sour(Residual);
		}
		else {
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0; 
			Residual[3] = 0.0; Residual[4] = 0.0; Residual[5] = (1.0/3.0)*TotalLoad*Area_Elem;
			node[Point_0]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0; 
			Residual[3] = 0.0; Residual[4] = 0.0; Residual[5] = (1.0/3.0)*TotalLoad*Area_Elem;
			node[Point_1]->AddRes_Sour(Residual);
			
			Residual[0] = 0.0; Residual[1] = 0.0; Residual[2] = 0.0; 
			Residual[3] = 0.0; Residual[4] = 0.0; Residual[5] = (1.0/3.0)*TotalLoad*Area_Elem;
			node[Point_2]->AddRes_Sour(Residual);
		}
		

	}
	
}

void CFEASolution::MPI_Send_Receive(CGeometry ***geometry, CSolution ****solution_container,
                                    CConfig **config, unsigned short iMGLevel, unsigned short iZone) {
	
#ifndef NO_MPI
	unsigned short iVar, iMarker;
	double *Displacement_Var;
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
				
				int send_to = SendRecv-1;
				
				double *Buffer_Send_Displacement = new double [nBuffer_Scalar];
				
				for (iVertex = 0; iVertex < geometry[iZone][iMGLevel]->nVertex[iMarker]; iVertex++) {
					iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();
					
					Displacement_Var = node[iPoint]->GetSolution();
					
					for (iVar = 0; iVar < nVar; iVar++)
						Buffer_Send_Displacement[iVar*nVertex+iVertex] = Displacement_Var[iVar];
					
				}
				
				MPI::COMM_WORLD.Bsend(Buffer_Send_Displacement,nBuffer_Scalar,MPI::DOUBLE,send_to, 0);
				
				delete[] Buffer_Send_Displacement;
				
			}
			
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
				
				/*--- Dimensionalization ---*/
				unsigned long nBuffer_Scalar = nVertex*nVar;
				
				int receive_from = abs(SendRecv)-1;
				
				double *Buffer_Receive_Displacement = new double [nBuffer_Scalar];
				
				MPI::COMM_WORLD.Recv(Buffer_Receive_Displacement,nBuffer_Scalar,MPI::DOUBLE,receive_from, 0);
				
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();
					for (iVar = 0; iVar < nVar; iVar++)
						node[iPoint]->SetSolution(iVar, Buffer_Receive_Displacement[iVar*nVertex+iVertex]);
										
				}
				
				delete[] Buffer_Receive_Displacement;
				
			}
		}
#endif
}

void CFEASolution::Postprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {

	/*--- Compute the gradient of the displacement ---*/
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

}

void CFEASolution::SetResidual_DualTime(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep, 
																		 unsigned short iMesh, unsigned short RunTime_EqSystem) {
	
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
	double a[3], b[3], c[3], d[3], Area_Local = 0.0, Volume_Local = 0.0, Time_Num;
	double *Coord_0_ = NULL, *Coord_1_= NULL, *Coord_2_= NULL, *Coord_3_= NULL;
	double Coord_0[3], Coord_1[3], Coord_2[3], Coord_3[3];
	unsigned short iDim, iVar, jVar;
	double Density = config->GetMaterialDensity(), TimeJac = 0.0;
	
	/*--- Numerical time step (this system is uncoditional stable... a very big number can be used) ---*/
	Time_Num = config->GetDelta_UnstTimeND();
	
	/*--- Loop through elements to compute contributions from the matrix
   blocks involving time. These contributions are also added to the 
   Jacobian w/ the time step. Spatial source terms are also computed. ---*/
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
    /*--- Get node numbers and their coordinate vectors ---*/
		Point_0 = geometry->elem[iElem]->GetNode(0);	Coord_0_ = geometry->node[Point_0]->GetCoord();
		Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1_ = geometry->node[Point_1]->GetCoord();	
		Point_2 = geometry->elem[iElem]->GetNode(2);	Coord_2_ = geometry->node[Point_2]->GetCoord();
		if (nDim == 3) { Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3_ = geometry->node[Point_3]->GetCoord(); }
		
		/*--- Modification of the goemtry due to the current displacement (Value of the solution) ---*/
		for (iDim = 0; iDim < nDim; iDim++) {
			Coord_0[iDim] = Coord_0_[iDim] + node[Point_0]->GetSolution()[iDim];
			Coord_1[iDim] = Coord_1_[iDim] + node[Point_1]->GetSolution()[iDim];
			Coord_2[iDim] = Coord_2_[iDim] + node[Point_2]->GetSolution()[iDim];
			if (nDim == 3) 
				Coord_3[iDim] = Coord_3_[iDim] + node[Point_3]->GetSolution()[iDim];	
		}
		/*--- Modification of the goemtry due to the current displacement (Value of the solution) ---*/
		
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
			xres[total_index] = 0.0;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				xsol[total_index] = ( U_time_nP1[iVar] - U_time_n[iVar] );
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				xsol[total_index] = ( U_time_nP1[iVar] - (4.0/3.0)*U_time_n[iVar] + (1.0/3.0)*U_time_nM1[iVar] );
		}
	}
	
	StiffMatrixTime.MatrixVectorProduct(xsol, xres);		
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			Residual[iVar] = xres[total_index];
		}
		node[iPoint]->SubtractRes_Visc(Residual);
	}
	
}

void CFEASolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	
  unsigned short iVar;
	unsigned long iPoint, total_index;
	double *local_ResVisc, *local_ResSour, Norm;
		
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
	
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		
		/*--- Read the residual ---*/
		local_ResVisc = node[iPoint]->GetResVisc();
		local_ResSour = node[iPoint]->GetResSour();
				
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			rhs[total_index] = local_ResVisc[iVar] + local_ResSour[iVar];
			xsol[total_index] = 0.0;
			AddRes_RMS(iVar, rhs[total_index]*rhs[total_index]);
      AddRes_Max(iVar, fabs(rhs[total_index]), geometry->node[iPoint]->GetGlobalIndex());
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
	Norm = 0.0;
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, xsol[iPoint*nVar+iVar]);
			Norm += xsol[iPoint*nVar+iVar]*xsol[iPoint*nVar+iVar];
		}
	}
	Norm = sqrt(Norm);
	SetTotal_CFEA(Norm);
	
  /*--- MPI solution ---*/
  SetSolution_MPI(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CFEASolution::SetInitialCondition(CGeometry **geometry, CSolution ***solution_container, CConfig *config, unsigned long ExtIter) {
	unsigned long iPoint;
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	unsigned short iDim; 
	double r[3], vel[3], *Coord, *Center, *Omega, Lref;
	
	/*--- Problem dimension and physical time step ---*/
	unsigned short nDim = geometry[MESH_0]->GetnDim();

	/*--- Center of rotation & angular velocity vector from config ---*/
	Center = config->GetRotAxisOrigin();
	Omega  = config->GetOmega_FreeStreamND();
	Lref   = config->GetLength_Ref();
	
	for (iPoint = 0; iPoint < geometry[MESH_0]->GetnPoint(); iPoint++) {
		
		/*--- Set initial boundary condition at the first iteration ---*/
		if ((ExtIter == 0) && (!restart)) {
			
			/*--- Coordinates of the current point ---*/
			Coord = geometry[MESH_0]->node[iPoint]->GetCoord();
			
			/*--- Calculate non-dim. position from rotation center ---*/
			for (iDim = 0; iDim < nDim; iDim++)
				r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
			if (nDim == 2) r[nDim] = 0.0;
			
			
			/*--- Compute the linear velocity v = Omega*r ---*/
			vel[0] = Omega[1]*r[2]-Omega[2]*r[1]; 
			vel[1] = -(Omega[0]*r[2]-Omega[2]*r[0]); 
			vel[2] = Omega[0]*r[1]-Omega[1]*r[0];
			
			
			if (nDim == 2) {
				Solution[0] = 0.0;	Solution[1] = 0.0;	
				Solution[2] = vel[0];	Solution[3] = vel[1];
			}
			else {
				Solution[0] = 0.0;	Solution[1] = 0.0;	Solution[2] = 0.0;	
				Solution[3] = vel[0];	Solution[4] = vel[1];	Solution[5] = vel[2];
			}
			
			node[iPoint]->SetSolution(Solution);

			if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
				node[iPoint]->Set_Solution_time_n();
				node[iPoint]->Set_Solution_time_n1();
			}
			
		}
	}
}

void CFEASolution::SetFEA_Load(CSolution ***flow_solution, CGeometry **fea_geometry, CGeometry **flow_geometry, 
															 CConfig *fea_config, CConfig *flow_config) {
	unsigned short iMarker;
	unsigned long iVertex, iPoint;
	double Pressure;	

#ifdef NO_MPI
	
	unsigned long iPoint_Donor;

	for (iMarker = 0; iMarker < fea_config->GetnMarker_All(); iMarker++) {
		if (fea_config->GetMarker_All_Boundary(iMarker) == FLOWLOAD_BOUNDARY) {
			for(iVertex = 0; iVertex < fea_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
				iPoint = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
				iPoint_Donor = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetDonorPoint();
				Pressure = 1.0; //(flow_solution[MESH_0][FLOW_SOL]->node[iPoint_Donor]->GetPressure()-101325.0);
				node[iPoint]->SetPressureValue(Pressure);
			}
		}
	}
	
#else
	
	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	double *Buffer_Send_U = new double [1];
	double *Buffer_Receive_U = new double [1];
	unsigned long jPoint;
	
	/*--- Do the send process, by the moment we are sending each 
	 node individually, this must be changed ---*/
	for (iMarker = 0; iMarker < flow_config->GetnMarker_All(); iMarker++) {
		/*--- There must be a better way to identify the marker that correspond with a fluid structure interation!!! ---*/
		if ((flow_config->GetMarker_All_Boundary(iMarker) == EULER_WALL) && 
			(flow_config->GetMarker_All_Moving(iMarker) == YES)) {
			for(iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
				iPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
				
				if (flow_geometry[MESH_0]->node[iPoint]->GetDomain()) {
					
					/*--- Find the associate pair to the original node (index and processor) ---*/
					jPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
					jProcessor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
					
					/*--- We only send the pressure that belong to other boundary ---*/
					if (jProcessor != rank) {
						Buffer_Send_U[0] = 1.0; //(flow_solution[MESH_0][FLOW_SOL]->node[iPoint]->GetPressure()-101325.0);;
						MPI::COMM_WORLD.Bsend(Buffer_Send_U, 1, MPI::DOUBLE, jProcessor, iPoint);
					}
					
				}		
			}
		}
	}
	
	/*--- Now the loop is over the fea points ---*/
	for (iMarker = 0; iMarker < fea_config->GetnMarker_All(); iMarker++) {
		if (fea_config->GetMarker_All_Boundary(iMarker) == LOAD_BOUNDARY) {
			for(iVertex = 0; iVertex < fea_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
				iPoint = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
				if (fea_geometry[MESH_0]->node[iPoint]->GetDomain()) {
					
					/*--- Find the associate pair to the original node ---*/
					jPoint = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
					jProcessor = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
					
					/*--- We only receive the information that belong to other boundary ---*/
					if (jProcessor != rank)
						MPI::COMM_WORLD.Recv(Buffer_Receive_U, 1, MPI::DOUBLE, jProcessor, jPoint);
					else
						Buffer_Receive_U[0] = 1.0; //(flow_solution[MESH_0][FLOW_SOL]->node[jPoint]->GetPressure()-101325.0);
					
					/*--- Store the solution for both points ---*/
					Pressure = Buffer_Receive_U[1];
					node[iPoint]->SetPressureValue(Pressure);
					
				}
			}
		}
	}
		
	delete[] Buffer_Send_U;
	delete[] Buffer_Receive_U;
		
#endif
}

