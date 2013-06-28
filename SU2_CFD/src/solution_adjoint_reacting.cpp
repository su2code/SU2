/*!
 * \file solution_adjoint_reacting.cpp
 * \brief Main subrotuines for solving adjoint problems (Euler, Navier-Stokes, etc.).
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
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

CAdjElectricSolution::CAdjElectricSolution(void) : CSolution() { }

CAdjElectricSolution::~CAdjElectricSolution(void) {
	
	unsigned short iVar, iDim;
	
	delete [] Residual;
	delete [] Residual_Max;
	delete [] Solution;
	
	if (nDim == 2) {
		for (iVar = 0; iVar < 3; iVar++)
			delete [] StiffMatrix_Elem[iVar];
	}
	
	if (nDim == 3) {
		for (iVar = 0; iVar < 4; iVar++)
			delete [] StiffMatrix_Elem[iVar];
	}
	
	delete [] StiffMatrix_Elem;
	delete [] StiffMatrix_Node;
	delete [] xsol;
	delete [] rhs;
	
	// Computation of gradients by least-squares
	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;
	
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
	
}

CAdjElectricSolution::CAdjElectricSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	
	unsigned long nPoint;
	unsigned short nMarker;
	
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	
	nPoint = geometry->GetnPoint();
	nDim = geometry->GetnDim();
	nMarker = config->GetnMarker_All(); 
	node = new CVariable*[nPoint];
	nVar = 1;		
	
	Residual = new double[nVar];
	Residual_Max = new double[nVar];
	Solution = new double[nVar];
	
	
	// - - - - STRUCTURES FOR SOLVING THE LINEAR SYSTEM - - - - - - 
	// Point to point stiffness matrix
	if (nDim == 2) {
		StiffMatrix_Elem = new double* [3];
		for (unsigned short iVar = 0; iVar < 3; iVar++) {
			StiffMatrix_Elem[iVar] = new double [3];
		}
	}
	
	if (nDim == 3) {
		StiffMatrix_Elem = new double* [4];
		for (unsigned short iVar = 0; iVar < 4; iVar++) {
			StiffMatrix_Elem[iVar] = new double [4];
		}
	}	
	
	StiffMatrix_Node = new double* [1];
	for (unsigned short iVar = 0; iVar < 1; iVar++) {
		StiffMatrix_Node[iVar] = new double [1];
	}
	
	// Initialization of the structure of the whole Jacobian
	InitializeStiffMatrixStructure(geometry, config);
	xsol = new double [nPoint*nVar];
	xres = new double [nPoint*nVar];
	rhs = new double [nPoint*nVar];
	
	// - - - - STRUCTURES LEAST SQUARES GRADIENTS - - - - - - - 
	// Computation of gradients by least squares
	Smatrix = new double* [nDim]; // S matrix := inv(R)*traspose(inv(R))
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Smatrix[iDim] = new double [nDim];
	
	cvector = new double* [nVar]; // c vector := transpose(WA)*(Wb)
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		cvector[iVar] = new double [nDim];
	
	
	// - - - - INITIALIZATION - - - - - - - 
	for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
		node[iPoint] = new CPotentialVariable(0.0, nDim, nVar, config);
	
	bool restart = config->GetRestart();
	
	if (!restart) {
		for (unsigned long iPoint=0; iPoint < nPoint; iPoint++)
			node[iPoint] = new CPotentialVariable(0.0, nDim, nVar, config);
	}
	else {
		string copy, mesh_filename;
		char buffer[50], cstr[200];
		mesh_filename = config->GetSolution_AdjFileName();
		copy.assign(mesh_filename);
		copy.erase (copy.end()-4, copy.end());
		strcpy (cstr, copy.c_str()); 
		if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT) sprintf (buffer, "_cd.dat"); 
		if (config->GetKind_ObjFunc() == LIFT_COEFFICIENT) sprintf (buffer, "_cl.dat");
		if (config->GetKind_ObjFunc() == SIDEFORCE_COEFFICIENT) sprintf (buffer, "_csf.dat"); 
		if (config->GetKind_ObjFunc() == PRESSURE_COEFFICIENT) sprintf (buffer, "_cp.dat"); 
		if (config->GetKind_ObjFunc() == MOMENT_X_COEFFICIENT) sprintf (buffer, "_cmx.dat"); 
		if (config->GetKind_ObjFunc() == MOMENT_Y_COEFFICIENT) sprintf (buffer, "_cmy.dat"); 
		if (config->GetKind_ObjFunc() == MOMENT_Z_COEFFICIENT) sprintf (buffer, "_cmz.dat"); 
		if (config->GetKind_ObjFunc() == EFFICIENCY) sprintf (buffer, "_eff.dat"); 
		if (config->GetKind_ObjFunc() == ELECTRIC_CHARGE) sprintf (buffer, "_ec.dat");
		if (config->GetKind_ObjFunc() == EQUIVALENT_AREA) sprintf (buffer, "_ea.dat"); 
		if (config->GetKind_ObjFunc() == NEARFIELD_PRESSURE) sprintf (buffer, "_nfp.dat"); 
		strcat(cstr, buffer);
		ifstream restart_file;
		restart_file.open(cstr, ios::in);
		
		if (restart_file.fail()) {
			cout << "There is no adjoint restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}
		unsigned long index;
		string text_line;	// definition of the text line
		
		for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			point_line >> index >> Solution[0];
			node[iPoint] = new CPotentialVariable(Solution[0], nDim, nVar, config);
		}
		restart_file.close();
	}		
	
}

void CAdjElectricSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		node[iPoint]->SetResidualZero(); // Inicialize the residual vector
	StiffMatrix.SetValZero();
}

void CAdjElectricSolution::Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
											  unsigned short iMesh) {
	
	
	unsigned long iPoint, nPoint = geometry->GetnPoint();
	unsigned short iVar = 0, iter = 0;
	double norm;
	
	// Build lineal system
	for (iPoint = 0; iPoint < nPoint; iPoint++) {		
		rhs[iPoint] = node[iPoint]->GetResidual(iVar);
		xsol[iPoint] = node[iPoint]->GetSolution(iVar);
	}
	
	// Solve the system
	norm = 1; iter = 0;
	while (norm > config->GetCauchy_Eps()) {
		iter++;
		norm = StiffMatrix.SGSIteration(rhs,xsol);
		if (iter == config->GetnExtIter()) break;		
	}
	
	SetRes_Max(0, norm);
	
	// Update solution
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		node[iPoint]->SetSolution(0,xsol[iPoint]);
	
}

void CAdjElectricSolution::Compute_Residual(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
											unsigned short iMesh) {
	
	unsigned long iPoint;
	unsigned short iVar = 0;
	
	// Build lineal system
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {		
		rhs[iPoint] = node[iPoint]->GetResidual(iVar);
		xsol[iPoint] = node[iPoint]->GetSolution(iVar);
		xres[iPoint] = 0.0;
	}
	
	StiffMatrix.MatrixVectorProduct(xsol,xres);
	
	// Update residual
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetResidual(0,xres[iPoint]-rhs[iPoint]);
	
	SetResidual_Smoothing(geometry, 10, 10.0);
	
	
}

void CAdjElectricSolution::SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
													CConfig *config, unsigned short iMesh) {
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, iPoint;
	double *Coord_0, *Coord_1, *Coord_2, a[3], b[3], Area_Local;
	unsigned short iDim, iNode;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		
		Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		
		for (iDim=0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
		
		Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0])/3.0;
		
		for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
			iPoint = geometry->elem[iElem]->GetNode(iNode);
			solver->SetCoord(geometry->node[iPoint]->GetCoord(), NULL, NULL, NULL);
			solver->SetVolume(Area_Local);
			solver->SetResidual(Residual, config);
			node[iPoint]->AddResidual(Residual);		
		}
		
	}
}

void CAdjElectricSolution::Galerkin_Method(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
										   CConfig *config, unsigned short iMesh) {
	
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
	double *Coord_0 = NULL, *Coord_1 = NULL, *Coord_2 = NULL, *Coord_3 = NULL;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) { // loop over edges
		// Points in edge 
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		if (nDim == 3) Point_3 = geometry->elem[iElem]->GetNode(3);
		
		// Points coordinates 
		Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		if (nDim == 3) Coord_3 = geometry->node[Point_3]->GetCoord();
		
		if (nDim == 2) solver->SetCoord(Coord_0, Coord_1, Coord_2);
		if (nDim == 3) solver->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
		
		solver->SetResidual(StiffMatrix_Elem, config);
		
		if (nDim == 2) {
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]; StiffMatrix.AddBlock(Point_0,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]; StiffMatrix.AddBlock(Point_0,Point_1,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]; StiffMatrix.AddBlock(Point_0,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]; StiffMatrix.AddBlock(Point_1,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]; StiffMatrix.AddBlock(Point_1,Point_1,StiffMatrix_Node);		
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]; StiffMatrix.AddBlock(Point_1,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]; StiffMatrix.AddBlock(Point_2,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]; StiffMatrix.AddBlock(Point_2,Point_1,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]; StiffMatrix.AddBlock(Point_2,Point_2,StiffMatrix_Node);			
		}
		
		if (nDim == 3) {
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]; StiffMatrix.AddBlock(Point_0,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]; StiffMatrix.AddBlock(Point_0,Point_1,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]; StiffMatrix.AddBlock(Point_0,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][3]; StiffMatrix.AddBlock(Point_0,Point_3,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]; StiffMatrix.AddBlock(Point_1,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]; StiffMatrix.AddBlock(Point_1,Point_1,StiffMatrix_Node);		
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]; StiffMatrix.AddBlock(Point_1,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][3]; StiffMatrix.AddBlock(Point_1,Point_3,StiffMatrix_Node);
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]; StiffMatrix.AddBlock(Point_2,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]; StiffMatrix.AddBlock(Point_2,Point_1,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]; StiffMatrix.AddBlock(Point_2,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][3]; StiffMatrix.AddBlock(Point_2,Point_3,StiffMatrix_Node);
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][0]; StiffMatrix.AddBlock(Point_3,Point_0,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][1]; StiffMatrix.AddBlock(Point_3,Point_1,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][2]; StiffMatrix.AddBlock(Point_3,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][3]; StiffMatrix.AddBlock(Point_3,Point_3,StiffMatrix_Node);
		}
	}
}

void CAdjElectricSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {
	unsigned long Point, iVertex;
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		Point = geometry->vertex[val_marker][iVertex]->GetNode();
		Solution[0]= 0.0;
		node[Point]->SetSolution(Solution);
		node[Point]->SetResidual(Solution);
		StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
	}
	
}

void CAdjElectricSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long Point, iVertex;
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		Point = geometry->vertex[val_marker][iVertex]->GetNode();
		Solution[0]= 0.0;
		node[Point]->SetSolution(Solution);
		node[Point]->SetResidual(Solution);
		StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
	}
}
