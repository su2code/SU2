/*!
 * \file solution_direct_wave.cpp
 * \brief Main subrotuines for solving the wave equation.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.1.
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

CWaveSolution::CWaveSolution(void) : CSolution() { }

CWaveSolution::CWaveSolution(CGeometry *geometry, CConfig *config) : CSolution() {

	unsigned long nPoint;
	unsigned short nMarker, iVar;
	
	nPoint = geometry->GetnPoint();
	nDim = geometry->GetnDim();
	nMarker = config->GetnMarker_All(); 
	node = new CVariable*[nPoint];
	nVar = 1;		

	Residual = new double[nVar];
	Residual_Max = new double[nVar];
	Solution = new double[nVar];

	/*--- Point to point stiffness matrix (only for tets and hexa)---*/
	Source_Vector	 = new double [nDim+1];
	
	StiffMatrix_Elem = new double* [nDim+1];
	for (iVar = 0; iVar < nDim+1; iVar++) {
		StiffMatrix_Elem[iVar] = new double [nDim+1];
	}

	StiffMatrix_Node = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		StiffMatrix_Node[iVar] = new double [nVar];
	}

	/*--- Initialization of the structure of the whole Jacobian ---*/
	Initialize_StiffMatrixSpace_Structure(geometry, config);
	Initialize_StiffMatrixTime_Structure(geometry, config);
	Initialize_Jacobian_Structure(geometry, config);

	xsol = new double [nPoint*nVar];
	xres = new double [nPoint*nVar];
	rhs = new double [nPoint*nVar];

	/*--- Computation of gradients by least squares ---*/
	Smatrix = new double* [nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Smatrix[iDim] = new double [nDim];

	cvector = new double* [nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		cvector[iVar] = new double [nDim];

	bool restart = (config->GetRestart() || config->GetRestart_Flow());

	if (!restart) {
		for (unsigned long iPoint=0; iPoint < nPoint; iPoint++)
			node[iPoint] = new CWaveVariable(0.0, nDim, nVar, config);
	}
	else {
		string mesh_filename = config->GetSolution_FlowFileName();
		ifstream restart_file;

		char *cstr; cstr = new char [mesh_filename.size()+1];
		strcpy (cstr, mesh_filename.c_str());
		restart_file.open(cstr, ios::in);
		if (restart_file.fail()) {
			cout << "There is no flow restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}
		unsigned long index;
		string text_line;

		for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			point_line >> index >> Solution[0];
			node[iPoint] = new CWaveVariable(Solution[0], nDim, nVar, config);
		}
		restart_file.close();
	}

}

CWaveSolution::~CWaveSolution(void) {

	unsigned short iVar, iDim;

	delete [] Residual;
	delete [] Residual_Max;
	delete [] Solution;

	for (iVar = 0; iVar < nDim+1; iVar++)
		delete [] StiffMatrix_Elem[iVar];


	for (iVar = 0; iVar < 1; iVar++)
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

void CWaveSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		node[iPoint]->Set_ResVisc_Zero();
		node[iPoint]->Set_ResSour_Zero();
	}
	
	StiffMatrixSpace.SetValZero();
	StiffMatrixTime.SetValZero();
	Jacobian.SetValZero();

}

void CWaveSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh) {
	
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double a[3], b[3], Area_Local, Time_Local, Delta;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
	unsigned short iDim;
	
	/*--- Modify the jacobian with the time contribution ---*/
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		
		Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		
		for (iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
		
		Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
		Time_Local = 1.0;
		Delta = (Area_Local/Time_Local);
		
		StiffMatrix_Node[0][0] = (2.0/12.0)*(Area_Local/Time_Local); Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (1.0/12.0)*(Area_Local/Time_Local); Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (1.0/12.0)*(Area_Local/Time_Local); Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (1.0/12.0)*(Area_Local/Time_Local); Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (2.0/12.0)*(Area_Local/Time_Local); Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (1.0/12.0)*(Area_Local/Time_Local); Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (1.0/12.0)*(Area_Local/Time_Local); Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (1.0/12.0)*(Area_Local/Time_Local); Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (2.0/12.0)*(Area_Local/Time_Local); Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
		
		StiffMatrix_Node[0][0] = (2.0/12.0)*Area_Local; StiffMatrixTime.AddBlock(Point_0, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (1.0/12.0)*Area_Local; StiffMatrixTime.AddBlock(Point_0, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (1.0/12.0)*Area_Local; StiffMatrixTime.AddBlock(Point_0, Point_2, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (1.0/12.0)*Area_Local; StiffMatrixTime.AddBlock(Point_1, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (2.0/12.0)*Area_Local; StiffMatrixTime.AddBlock(Point_1, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (1.0/12.0)*Area_Local; StiffMatrixTime.AddBlock(Point_1, Point_2, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (1.0/12.0)*Area_Local; StiffMatrixTime.AddBlock(Point_2, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (1.0/12.0)*Area_Local; StiffMatrixTime.AddBlock(Point_2, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = (2.0/12.0)*Area_Local; StiffMatrixTime.AddBlock(Point_2, Point_2, StiffMatrix_Node);
		
//		solver->SetVolume(Area_Local);
//		solver->SetCoord(Coord_0, Coord_1, Coord_2);
//		solver->SetResidual(Residual, config);
		
//		node[Point_0]->AddRes_Sour(Residual);
//		node[Point_1]->AddRes_Sour(Residual);
//		node[Point_2]->AddRes_Sour(Residual);

	}

}

void CWaveSolution::Galerkin_Method(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh) {

	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, iPoint, total_index;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, SoundSpeed;
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
		SoundSpeed = 1.0; //343.2;
		
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]*SoundSpeed; 
		StiffMatrixSpace.AddBlock(Point_0, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]*SoundSpeed;  
		StiffMatrixSpace.AddBlock(Point_0, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]*SoundSpeed;  
		StiffMatrixSpace.AddBlock(Point_0, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]*SoundSpeed;  
		StiffMatrixSpace.AddBlock(Point_1, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]*SoundSpeed;  
		StiffMatrixSpace.AddBlock(Point_1, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]*SoundSpeed;  
		StiffMatrixSpace.AddBlock(Point_1, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]*SoundSpeed;  
		StiffMatrixSpace.AddBlock(Point_2, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]*SoundSpeed;  
		StiffMatrixSpace.AddBlock(Point_2, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
		StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]*SoundSpeed;  
		StiffMatrixSpace.AddBlock(Point_2, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);

	}
	
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			xsol[total_index] = node[iPoint]->GetSolution(iVar);
			xres[total_index] = 0.0;
		}
	
	StiffMatrixSpace.MatrixVectorProduct(xsol, xres, geometry->GetnPointDomain());
	
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			Residual[iVar] = xres[total_index];
		}
		node[iPoint]->AddRes_Visc(Residual);
	}
}

void CWaveSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
																	unsigned short val_marker) {
	unsigned long Point, iVertex;
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		Point = geometry->vertex[val_marker][iVertex]->GetNode();
		Solution[0] = 1.0; 
		node[Point]->SetSolution(Solution);
		node[Point]->SetSolution_Old(Solution);
		Residual[0] = 0.0; node[Point]->SetRes_Visc(Residual); 
		Residual[0] = 0.0; node[Point]->SetRes_Sour(Residual);
		Jacobian.DeleteValsRowi(Point);
	}

}

void CWaveSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
		unsigned short val_marker) {
	unsigned long Point, iVertex;

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		Point = geometry->vertex[val_marker][iVertex]->GetNode();
		Solution[0]= 0.0;
		node[Point]->SetSolution(Solution);
		node[Point]->SetSolution_Old(Solution);
		Residual[0] = 0.0; node[Point]->SetRes_Visc(Residual); 
		Residual[0]= 0.0; node[Point]->SetRes_Sour(Residual);
		Jacobian.DeleteValsRowi(Point);
	}

}

void CWaveSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Res, *local_ResVisc, *local_ResSour;
		
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);
	
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		/*--- Read the residual ---*/
		local_ResVisc = node[iPoint]->GetResVisc();
		local_ResSour = node[iPoint]->GetResSour();
				
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			Res = local_ResVisc[iVar] + local_ResSour[iVar];
			rhs[total_index] = -Res;
			xsol[total_index] = 0.0;
			AddRes_Max( iVar, Res*Res );
		}
	}
	
	/*--- Solve the system ---*/
	if (config->GetKind_Linear_Solver() == SYM_GAUSS_SEIDEL) Jacobian.SGSSolution(rhs, xsol, config->GetLinear_Solver_Error(), 
																																								config->GetLinear_Solver_Iter(), true, geometry, config);
	if (config->GetKind_Linear_Solver() == BCGSTAB) Jacobian.BCGSTABSolution(rhs, xsol, config->GetLinear_Solver_Error(), 
																																					 config->GetLinear_Solver_Iter(), true, geometry, config);
	if (config->GetKind_Linear_Solver() == LU_SGS) Jacobian.LU_SGSIteration(rhs, xsol, geometry, config);
	
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, xsol[iPoint*nVar+iVar]);
		}
	}
	
#ifdef NO_MPI
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_Max(iVar, sqrt(GetRes_Max(iVar)));
	}
#endif

}
