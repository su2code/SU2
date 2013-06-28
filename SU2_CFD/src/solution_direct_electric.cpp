/*!
 * \file solution_direct_reacting_diatomic.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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

CElectricSolution::CElectricSolution(void) : CSolution() { }

CElectricSolution::CElectricSolution(CGeometry *geometry, CConfig *config) : CSolution() {

	unsigned long nPoint;
	unsigned short nMarker;

	nPoint = geometry->GetnPoint();
	nDim = geometry->GetnDim();
	nMarker = config->GetnMarker_All(); 
	node = new CVariable*[nPoint];
	nVar = 1;		
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	Residual = new double[nVar];
	Residual_Max = new double[nVar];
	Solution = new double[nVar];

	/*--- Point to point stiffness matrix ---*/
	if (nDim == 2) {
		StiffMatrix_Elem = new double* [3];
		Source_Vector	 = new double [3];
		for (unsigned short iVar = 0; iVar < 3; iVar++) {
			StiffMatrix_Elem[iVar] = new double [3];
		}
	}

	if (nDim == 3) {
		StiffMatrix_Elem = new double* [4];
		Source_Vector	 = new double [4];
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
	bool restart = (config->GetRestart() || config->GetRestart_Flow());

	if (!restart) {
		for (unsigned long iPoint=0; iPoint < nPoint; iPoint++)
			node[iPoint] = new CPotentialVariable(0.0, nDim, nVar, config);
	}
	else {
		string mesh_filename = config->GetSolution_FlowFileName();
		ifstream restart_file;	// definition of file

		char *cstr; cstr = new char [mesh_filename.size()+1];
		strcpy (cstr, mesh_filename.c_str());
		restart_file.open(cstr, ios::in);	// Open the file
		if (restart_file.fail()) {
			cout << "There is no flow restart file!!" << endl;
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

CElectricSolution::~CElectricSolution(void) {

	unsigned short iVar, iDim;

	delete [] Residual;
	delete [] Residual_Max;
	delete [] Solution;
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
	delete [] xsol;
	delete [] xres;
	delete [] rhs;

	// Computation of gradients by least-squares
	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
}

void CElectricSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		node[iPoint]->SetResidualZero();

	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry); 
	if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
			(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) SetSolution_Gradient_LS(geometry, config);

	StiffMatrix.SetValZero();
}

void CElectricSolution::Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
		unsigned short iMesh) {

	unsigned long iPoint;
	unsigned short iVar = 0;
	double norm = 1E6;
	unsigned long iter_max = 100000;

	/*--- Build lineal system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {		
		rhs[iPoint] = node[iPoint]->GetResidual(iVar);
		xsol[iPoint] = node[iPoint]->GetSolution(iVar);
	}

	/*--- Solve the system ---*/

	/*	CSysVector rhs_vec(geometry->GetnPoint(), 1, rhs);
	 CSysVector sol_vec(geometry->GetnPoint(), 1, xsol);
	 //StiffMatrix.BuildJacobiPreconditioner();
	 CMatrixVectorProduct* mat_vec = new CSparseMatrixVectorProduct(StiffMatrix);
	 //CPreconditioner* precond = new CSparseMatrixPreconditioner(StiffMatrix);
	 CPreconditioner* precond = new CIdentityPreconditioner();
	 CSysSolve system;
	 system.ConjugateGradient(rhs_vec, sol_vec, *mat_vec, *precond, 1E-10, iter_max, true);
	 sol_vec.CopyToArray(xsol);
	 delete mat_vec;
	 delete precond;*/

	StiffMatrix.CGSolution(rhs, xsol, 1E-14, iter_max, true);
	SetRes_Max(0, norm);
	/*--- Update solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetSolution(0,xsol[iPoint]);

}

void CElectricSolution::Compute_Residual(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
		unsigned short iMesh) {

	unsigned long iPoint;
	unsigned short iVar = 0;


	// Build linear system
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

void CElectricSolution::SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh) {
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double a[3], b[3], Area_Local;
	double Local_Delta_Time_0, Local_Delta_Time_1, Local_Delta_Time_2, Local_Delta_Time;
	double **Gradient_0, **Gradient_1, **Gradient_2;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
	unsigned short iDim;
	double Lambda_iFluid;
	double  dx, u, c;

	if (nDim !=3) {
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

			Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

			Gradient_0 = solution_container[PLASMA_SOL]->node[Point_0]->GetGradient();
			Gradient_1 = solution_container[PLASMA_SOL]->node[Point_1]->GetGradient();
			Gradient_2 = solution_container[PLASMA_SOL]->node[Point_2]->GetGradient();

			solver->SetVolume(Area_Local);

			Lambda_iFluid = max(solution_container[PLASMA_SOL]->node[Point_0]->GetMax_Lambda_Inv(0), solution_container[PLASMA_SOL]->node[Point_0]->GetMax_Lambda_Inv(1));
			Lambda_iFluid = max(solution_container[PLASMA_SOL]->node[Point_0]->GetMax_Lambda_Inv(2), Lambda_iFluid);
			Local_Delta_Time_0 = config->GetCFL(iMesh)*Area_Local / Lambda_iFluid;

			Lambda_iFluid = max(solution_container[PLASMA_SOL]->node[Point_1]->GetMax_Lambda_Inv(0), solution_container[PLASMA_SOL]->node[Point_1]->GetMax_Lambda_Inv(1));
			Lambda_iFluid = max(solution_container[PLASMA_SOL]->node[Point_1]->GetMax_Lambda_Inv(2), Lambda_iFluid);
			Local_Delta_Time_1 = config->GetCFL(iMesh)*Area_Local / Lambda_iFluid;

			Lambda_iFluid = max(solution_container[PLASMA_SOL]->node[Point_2]->GetMax_Lambda_Inv(0), solution_container[PLASMA_SOL]->node[Point_2]->GetMax_Lambda_Inv(1));
			Lambda_iFluid = max(solution_container[PLASMA_SOL]->node[Point_2]->GetMax_Lambda_Inv(2), Lambda_iFluid);
			Local_Delta_Time_2 = config->GetCFL(iMesh)*Area_Local / Lambda_iFluid;

			Local_Delta_Time   = (Local_Delta_Time_0 + Local_Delta_Time_1 + Local_Delta_Time_2 ) / 3.0;

			u = 4800;
			c = 87110;
			c = 800;
			dx = 0.004/81;
			Local_Delta_Time = config->GetCFL(iMesh) * dx/(u+c);
			if (nDim == 2) solver->SetCoord(Coord_0, Coord_1, Coord_2);
			solver->SetTimeStep(Local_Delta_Time);
			solver->SetConservative(solution_container[PLASMA_SOL]->node[Point_0]->GetSolution(), solution_container[PLASMA_SOL]->node[Point_1]->GetSolution(), solution_container[PLASMA_SOL]->node[Point_2]->GetSolution());
			solver->SetConsVarGradient(Gradient_0, Gradient_1, Gradient_2 );
			solver->SetResidual(Source_Vector, config);

			node[Point_0]->AddResidual(Source_Vector[0]);
			node[Point_1]->AddResidual(Source_Vector[1]);
			node[Point_2]->AddResidual(Source_Vector[2]);

		}

		for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

			Point_0 = geometry->elem[iElem]->GetNode(3);
			Point_1 = geometry->elem[iElem]->GetNode(1);
			Point_2 = geometry->elem[iElem]->GetNode(2);

			Coord_0 = geometry->node[Point_0]->GetCoord();
			Coord_1 = geometry->node[Point_1]->GetCoord();
			Coord_2 = geometry->node[Point_2]->GetCoord();

			for (iDim=0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}

			Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

			Gradient_0 = solution_container[PLASMA_SOL]->node[Point_0]->GetGradient();
			Gradient_1 = solution_container[PLASMA_SOL]->node[Point_1]->GetGradient();
			Gradient_2 = solution_container[PLASMA_SOL]->node[Point_2]->GetGradient();

			solver->SetVolume(Area_Local);

			Lambda_iFluid = max(solution_container[PLASMA_SOL]->node[Point_0]->GetMax_Lambda_Inv(0), solution_container[PLASMA_SOL]->node[Point_0]->GetMax_Lambda_Inv(1));
			Lambda_iFluid = max(solution_container[PLASMA_SOL]->node[Point_0]->GetMax_Lambda_Inv(2), Lambda_iFluid);
			Local_Delta_Time_0 = config->GetCFL(iMesh)*Area_Local / Lambda_iFluid;

			Lambda_iFluid = max(solution_container[PLASMA_SOL]->node[Point_1]->GetMax_Lambda_Inv(0), solution_container[PLASMA_SOL]->node[Point_1]->GetMax_Lambda_Inv(1));
			Lambda_iFluid = max(solution_container[PLASMA_SOL]->node[Point_1]->GetMax_Lambda_Inv(2), Lambda_iFluid);
			Local_Delta_Time_1 = config->GetCFL(iMesh)*Area_Local / Lambda_iFluid;

			Lambda_iFluid = max(solution_container[PLASMA_SOL]->node[Point_2]->GetMax_Lambda_Inv(0), solution_container[PLASMA_SOL]->node[Point_2]->GetMax_Lambda_Inv(1));
			Lambda_iFluid = max(solution_container[PLASMA_SOL]->node[Point_2]->GetMax_Lambda_Inv(2), Lambda_iFluid);
			Local_Delta_Time_2 = config->GetCFL(iMesh)*Area_Local / Lambda_iFluid;

			Local_Delta_Time   = (Local_Delta_Time_0 + Local_Delta_Time_1 + Local_Delta_Time_2 ) / 3.0;

			u = 4800;
			c = 87110;
			c = 732.0;
			dx = 0.004/81;
			Local_Delta_Time = config->GetCFL(iMesh) * dx/(u+c);
			if (nDim == 2) solver->SetCoord(Coord_0, Coord_1, Coord_2);
			solver->SetTimeStep(Local_Delta_Time);
			solver->SetConservative(solution_container[PLASMA_SOL]->node[Point_0]->GetSolution(), solution_container[PLASMA_SOL]->node[Point_1]->GetSolution(), solution_container[PLASMA_SOL]->node[Point_2]->GetSolution());
			solver->SetConsVarGradient(Gradient_0, Gradient_1, Gradient_2 );
			solver->SetResidual(Source_Vector, config);
			node[Point_0]->AddResidual(Source_Vector[0]);
			node[Point_1]->AddResidual(Source_Vector[1]);
			node[Point_2]->AddResidual(Source_Vector[2]);

		}
	}

}

void CElectricSolution::Galerkin_Method(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh) {

	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);

		Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();

		solver->SetCoord(Coord_0, Coord_1, Coord_2);

		solver->SetResidual(StiffMatrix_Elem, config);

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

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		Point_0 = geometry->elem[iElem]->GetNode(3);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);

		Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();

		solver->SetCoord(Coord_0, Coord_1, Coord_2);

		solver->SetResidual(StiffMatrix_Elem, config);

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
}

void CElectricSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {
	/*	unsigned long Point, iVertex;

	 for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
	 Point = geometry->vertex[val_marker][iVertex]->GetNode();
	 Solution[0]= 0.0;
	 node[Point]->SetSolution(Solution);
	 node[Point]->SetResidual(Solution);
	 StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
	 }*/

}

void CElectricSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
		unsigned short val_marker) {
	/*	unsigned long iPoint, iVertex;

	 for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
	 iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
	 Solution[0]= 0;
	 node[iPoint]->SetSolution(Solution);
	 node[iPoint]->SetResidual(Solution);
	 StiffMatrix.DeleteValsRowi(iPoint); // & includes 1 in the diagonal
	 }*/

}

void CElectricSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
		unsigned short val_marker) {
	unsigned long Point, iVertex;

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		Point = geometry->vertex[val_marker][iVertex]->GetNode();
		Solution[0]= 0;
		node[Point]->SetSolution(Solution);
		node[Point]->SetResidual(Solution);
		StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
	}
}

void CElectricSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
		unsigned short val_marker) {
	/*	unsigned long iPoint, iVertex;

	 for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
	 iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
	 Solution[0]= 0.0;
	 node[iPoint]->SetSolution(Solution);
	 node[iPoint]->SetResidual(Solution);
	 StiffMatrix.DeleteValsRowi(iPoint); // & includes 1 in the diagonal
	 }
	 */
}


