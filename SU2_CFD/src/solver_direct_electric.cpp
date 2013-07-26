/*!
 * \file solution_direct_electric.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
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

#include "../include/solver_structure.hpp"

CElectricSolver::CElectricSolver(void) : CSolver() { }

CElectricSolver::CElectricSolver(CGeometry *geometry, CConfig *config) : CSolver() {

	unsigned long nPoint;
	unsigned short nMarker;

	nDim = geometry->GetnDim();
	nMarker = config->GetnMarker_All(); 
	node = new CVariable*[nPoint];
	nVar = 1;		
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
	Residual = new double[nVar]; Residual_RMS = new double[nVar];
	Solution = new double[nVar];
  Residual_Max = new double[nVar]; Point_Max = new unsigned long[nVar];
 
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

	/*--- Initialization of the structure of the whole Jacobian ---*/
	StiffMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);

  /*--- Solution and residual vectors ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

	/*--- Computation of gradients by least squares ---*/
	Smatrix = new double* [nDim]; // S matrix := inv(R)*traspose(inv(R))
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Smatrix[iDim] = new double [nDim];

	cvector = new double* [nVar]; // c vector := transpose(WA)*(Wb)
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		cvector[iVar] = new double [nDim];


	/*--- Always instantiate and initialize the variable to a zero value. ---*/
	for (unsigned long iPoint=0; iPoint < nPoint; iPoint++)
		node[iPoint] = new CPotentialVariable(0.0, nDim, nVar, config);

}

CElectricSolver::~CElectricSolver(void) {

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

	/*--- Computation of gradients by least-squares ---*/
	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
}

void CElectricSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container,
																			CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		LinSysRes.SetBlock_Zero(iPoint);

	StiffMatrix.SetValZero();
}

void CElectricSolver::Solve_LinearSystem(CGeometry *geometry, CSolver **solver_container, 
		CConfig *config, unsigned short iMesh) {
	unsigned long iPoint;
	unsigned short iVar = 0;
	double norm = 1E6;
	int iter_max = 10000;

	/*--- Build lineal system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		LinSysRes[iPoint] = LinSysRes.GetBlock(iPoint, iVar);
		LinSysSol[iPoint] = node[iPoint]->GetSolution(iVar);
	}

	/*--- Solve the system ---*/
	CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(StiffMatrix, geometry, config);

	StiffMatrix.BuildJacobiPreconditioner();

	CPreconditioner* precond = NULL;
	Jacobian.BuildJacobiPreconditioner();
	precond = new CJacobiPreconditioner(Jacobian, geometry, config);			

	CSysSolve system;
	system.ConjugateGradient(LinSysRes, LinSysSol, *mat_vec, *precond, 1E-12, iter_max, true);	

	delete mat_vec; 
	delete precond;

	SetRes_RMS(0, norm);
	/*--- Update solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		node[iPoint]->SetSolution(0,LinSysSol[iPoint]);

	}

}

void CElectricSolver::Compute_Residual(CGeometry *geometry, CSolver **solver_container, CConfig *config, 
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
 * \brief Source terms of the electric solver
 * \author A. Lonkar
 */
void CElectricSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
		CConfig *config, unsigned short iMesh) {

	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
	double a[3], b[3],c[3], d[3], Area_Local,Volume_Local;
	//	double Local_Delta_Time;
	double **Gradient_0, **Gradient_1, **Gradient_2, **Gradient_3;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;;
	unsigned short iDim;
	double  dt;
	//	double  dx, u, c;
	bool MacCormack_relaxation = (config->GetMacCormackRelaxation());

	if (nDim == 2) {
		if (config->GetElectricSolver()) {

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

				Gradient_0 = node[Point_0]->GetPlasmaRhoUGradient();
				Gradient_1 = node[Point_1]->GetPlasmaRhoUGradient();
				Gradient_2 = node[Point_2]->GetPlasmaRhoUGradient();

				numerics->SetVolume(Area_Local);

				dt = node[Point_0]->GetPlasmaTimeStep();

				/*		u = 4800;
				 c = 87110;
				 c = 800;
				 dx = 0.004/81;
				 Local_Delta_Time = config->GetCFL(iMesh) * dx/(u+c);
				 numerics->SetTimeStep(Local_Delta_Time);
				 */

				numerics->SetCoord(Coord_0, Coord_1, Coord_2);
				numerics->SetTimeStep(dt);
				numerics->SetChargeDensity(node[Point_0]->GetChargeDensity(), node[Point_1]->GetChargeDensity(), node[Point_2]->GetChargeDensity(), node[Point_3]->GetChargeDensity());
				numerics->SetConsVarGradient(Gradient_0, Gradient_1, Gradient_2 );
				numerics->ComputeResidual_MacCormack(Source_Vector, config);

				LinSysRes.AddBlock(Point_0, &Source_Vector[0]);
				LinSysRes.AddBlock(Point_1, &Source_Vector[1]);
				LinSysRes.AddBlock(Point_2, &Source_Vector[2]);

				if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE) {

					Point_0 = geometry->elem[iElem]->GetNode(3);
					Point_1 = geometry->elem[iElem]->GetNode(0);
					Point_2 = geometry->elem[iElem]->GetNode(2);

					Coord_0 = geometry->node[Point_0]->GetCoord();
					Coord_1 = geometry->node[Point_1]->GetCoord();
					Coord_2 = geometry->node[Point_2]->GetCoord();

					for (iDim=0; iDim < nDim; iDim++) {
						a[iDim] = Coord_0[iDim]-Coord_2[iDim];
						b[iDim] = Coord_1[iDim]-Coord_2[iDim];
					}

					Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

					Gradient_0 = node[Point_0]->GetPlasmaRhoUGradient();
					Gradient_1 = node[Point_1]->GetPlasmaRhoUGradient();
					Gradient_2 = node[Point_2]->GetPlasmaRhoUGradient();

					numerics->SetVolume(Area_Local);

					/*		u = 4800;
					 c = 87110;
					 c = 732.0;
					 dx = 0.004/81;
					 Local_Delta_Time = config->GetCFL(iMesh) * dx/(u+c);
					 numerics->SetTimeStep(Local_Delta_Time);
					 */

					dt = node[Point_0]->GetPlasmaTimeStep();
					numerics->SetCoord(Coord_0, Coord_1, Coord_2);
					numerics->SetTimeStep(dt);
					numerics->SetChargeDensity(node[Point_0]->GetChargeDensity(), node[Point_1]->GetChargeDensity(), node[Point_2]->GetChargeDensity(), node[Point_3]->GetChargeDensity());
					numerics->SetConsVarGradient(Gradient_0, Gradient_1, Gradient_2 );
					numerics->ComputeResidual_MacCormack(Source_Vector, config);
					LinSysRes.AddBlock(Point_0, &Source_Vector[0]);
					LinSysRes.AddBlock(Point_1, &Source_Vector[1]);
					LinSysRes.AddBlock(Point_2, &Source_Vector[2]);
				}
			}
		}
	}
	if(nDim == 3) {
		if (config->GetElectricSolver()) {
			for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
				Point_0 = geometry->elem[iElem]->GetNode(0);	Coord_0 = geometry->node[Point_0]->GetCoord();
				Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1 = geometry->node[Point_1]->GetCoord();
				Point_2 = geometry->elem[iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord();
				Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3 = geometry->node[Point_3]->GetCoord();

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
				numerics->SetVolume(Volume_Local);
				numerics->SetChargeDensity(node[Point_0]->GetChargeDensity(), node[Point_1]->GetChargeDensity(), node[Point_2]->GetChargeDensity(), node[Point_3]->GetChargeDensity());

				if (MacCormack_relaxation) {

					Gradient_0 = node[Point_0]->GetPlasmaRhoUGradient();
					Gradient_1 = node[Point_1]->GetPlasmaRhoUGradient();
					Gradient_2 = node[Point_2]->GetPlasmaRhoUGradient();
					Gradient_3 = node[Point_3]->GetPlasmaRhoUGradient();
					numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
					numerics->SetConsVarGradient(Gradient_0, Gradient_1, Gradient_2, Gradient_3 );
					numerics->ComputeResidual_MacCormack(Source_Vector, config);
				}
				else numerics->ComputeResidual(Source_Vector, config);

				LinSysRes.AddBlock(Point_0, &Source_Vector[0]);
				LinSysRes.AddBlock(Point_1, &Source_Vector[1]);
				LinSysRes.AddBlock(Point_2, &Source_Vector[2]);
				LinSysRes.AddBlock(Point_3, &Source_Vector[3]);
			}
		}
	}
}

void CElectricSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
		CConfig *config, unsigned short iMesh) {
}


/*!
 * \method Copy_Zone_Solution
 * \brief Copy solution from solver 1 into solver 2
 * \author A. Lonkar
 */
void CElectricSolver::Copy_Zone_Solution(CSolver ***solver1_solution, CGeometry **solver1_geometry, CConfig *solver1_config,
		CSolver ***solver2_solution, CGeometry **solver2_geometry, CConfig *solver2_config) {
	unsigned long iPoint;
	unsigned short iDim;
	double neg_EFvalue;
	double *E_field = new double [nDim];

	for (iPoint = 0; iPoint < solver1_geometry[MESH_0]->GetnPointDomain(); iPoint++) {
		for (iDim =0; iDim < nDim; iDim ++) {
			neg_EFvalue = solver1_solution[MESH_0][ELEC_SOL]->node[iPoint]->GetGradient(0,iDim);
			E_field[iDim] = -1.0*neg_EFvalue;
		}
		solver2_solution[MESH_0][PLASMA_SOL]->node[iPoint]->SetElectricField(E_field);
	}
};

/*!
 * \method Galerkin_Method
 * \brief calculate the element stiffness matrix
 * \author A. Lonkar
 */
void CElectricSolver::Galerkin_Method(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
		CConfig *config, unsigned short iMesh) {

	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3 = NULL;
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

		for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
			if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE) {

				Point_0 = geometry->elem[iElem]->GetNode(3);
				Point_1 = geometry->elem[iElem]->GetNode(0);
				Point_2 = geometry->elem[iElem]->GetNode(2);

				Coord_0 = geometry->node[Point_0]->GetCoord();
				Coord_1 = geometry->node[Point_1]->GetCoord();
				Coord_2 = geometry->node[Point_2]->GetCoord();

				numerics->SetCoord(Coord_0, Coord_1, Coord_2);
				numerics->ComputeResidual(StiffMatrix_Elem, config);
				AddStiffMatrix(StiffMatrix_Elem,Point_0, Point_1, Point_2, Point_3);
			}
		}
	}

	if (nDim == 3 ) {

		for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {

				Point_0 = geometry->elem[iElem]->GetNode(0); 	Coord_0 = geometry->node[Point_0]->GetCoord();
				Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1 = geometry->node[Point_1]->GetCoord();
				Point_2 = geometry->elem[iElem]->GetNode(2); 	Coord_2 = geometry->node[Point_2]->GetCoord();
				Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3 = geometry->node[Point_3]->GetCoord();

				numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
				numerics->ComputeResidual(StiffMatrix_Elem, config);
				AddStiffMatrix(StiffMatrix_Elem, Point_0, Point_1, Point_2, Point_3);
			}

			if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {

				/* Tetrahedron: 1, nodes: [0,2,3,6] */
				Point_0 = geometry->elem[iElem]->GetNode(0); 	Coord_0 = geometry->node[Point_0]->GetCoord();
				Point_1 = geometry->elem[iElem]->GetNode(2);	Coord_1 = geometry->node[Point_1]->GetCoord();
				Point_2 = geometry->elem[iElem]->GetNode(3); 	Coord_2 = geometry->node[Point_2]->GetCoord();
				Point_3 = geometry->elem[iElem]->GetNode(6);	Coord_3 = geometry->node[Point_3]->GetCoord();

				numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
				numerics->ComputeResidual(StiffMatrix_Elem, config);
				AddStiffMatrix(StiffMatrix_Elem, Point_0, Point_1, Point_2, Point_3);
				/******************/

				/* Tetrahedron: 2, nodes: [0,3,7,6]  */
				Point_0 = geometry->elem[iElem]->GetNode(0); 	Coord_0 = geometry->node[Point_0]->GetCoord();
				Point_1 = geometry->elem[iElem]->GetNode(3);	Coord_1 = geometry->node[Point_1]->GetCoord();
				Point_2 = geometry->elem[iElem]->GetNode(7); 	Coord_2 = geometry->node[Point_2]->GetCoord();
				Point_3 = geometry->elem[iElem]->GetNode(6);	Coord_3 = geometry->node[Point_3]->GetCoord();

				numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
				numerics->ComputeResidual(StiffMatrix_Elem, config);
				AddStiffMatrix(StiffMatrix_Elem, Point_0, Point_1, Point_2, Point_3);
				/******************/

				/* Tetrahedron: 3, nodes: [0,7,4, 6]  */
				Point_0 = geometry->elem[iElem]->GetNode(0); 	Coord_0 = geometry->node[Point_0]->GetCoord();
				Point_1 = geometry->elem[iElem]->GetNode(7);	Coord_1 = geometry->node[Point_1]->GetCoord();
				Point_2 = geometry->elem[iElem]->GetNode(4); 	Coord_2 = geometry->node[Point_2]->GetCoord();
				Point_3 = geometry->elem[iElem]->GetNode(6);	Coord_3 = geometry->node[Point_3]->GetCoord();

				numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
				numerics->ComputeResidual(StiffMatrix_Elem, config);
				AddStiffMatrix(StiffMatrix_Elem, Point_0, Point_1, Point_2, Point_3);
				/******************/

				/* Tetrahedron: 4, nodes: [0,5,6,4]  */
				Point_0 = geometry->elem[iElem]->GetNode(0); 	Coord_0 = geometry->node[Point_0]->GetCoord();
				Point_1 = geometry->elem[iElem]->GetNode(5);	Coord_1 = geometry->node[Point_1]->GetCoord();
				Point_2 = geometry->elem[iElem]->GetNode(6); 	Coord_2 = geometry->node[Point_2]->GetCoord();
				Point_3 = geometry->elem[iElem]->GetNode(4);	Coord_3 = geometry->node[Point_3]->GetCoord();

				numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
				numerics->ComputeResidual(StiffMatrix_Elem, config);
				AddStiffMatrix(StiffMatrix_Elem, Point_0, Point_1, Point_2, Point_3);
				/******************/

				/* Tetrahedron: 5, nodes: [1,5,6,0]  */
				Point_0 = geometry->elem[iElem]->GetNode(1); 	Coord_0 = geometry->node[Point_0]->GetCoord();
				Point_1 = geometry->elem[iElem]->GetNode(5);	Coord_1 = geometry->node[Point_1]->GetCoord();
				Point_2 = geometry->elem[iElem]->GetNode(6); 	Coord_2 = geometry->node[Point_2]->GetCoord();
				Point_3 = geometry->elem[iElem]->GetNode(0);	Coord_3 = geometry->node[Point_3]->GetCoord();

				numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
				numerics->ComputeResidual(StiffMatrix_Elem, config);
				AddStiffMatrix(StiffMatrix_Elem, Point_0, Point_1, Point_2, Point_3);
				/******************/

				/* Tetrahedron: 6, nodes: [1,6,2,0]  */
				Point_0 = geometry->elem[iElem]->GetNode(1); 	Coord_0 = geometry->node[Point_0]->GetCoord();
				Point_1 = geometry->elem[iElem]->GetNode(6);	Coord_1 = geometry->node[Point_1]->GetCoord();
				Point_2 = geometry->elem[iElem]->GetNode(2); 	Coord_2 = geometry->node[Point_2]->GetCoord();
				Point_3 = geometry->elem[iElem]->GetNode(0);	Coord_3 = geometry->node[Point_3]->GetCoord();

				numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
				numerics->ComputeResidual(StiffMatrix_Elem, config);
				AddStiffMatrix(StiffMatrix_Elem, Point_0, Point_1, Point_2, Point_3);
				/******************/
			}

		}
	}
}

/*!
 * \method AddStiffMatrix
 * \brief Assemble the Global stiffness matrix
 * \author A. Lonkar
 */
void CElectricSolver::AddStiffMatrix(double **StiffMatrix_Elem, unsigned long Point_0, unsigned long Point_1,unsigned long Point_2, unsigned long Point_3) {

	if (nDim == 2 ) {
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

/*!
 * \method BC_Euler_Wall
 * \brief Dirichlet/Neumann boundary condition
 * \author A. Lonkar
 */
void CElectricSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
		unsigned short val_marker) {
	unsigned long Point, iVertex;

	/*--- Identify if a boundary is Dirichlet or Neumman ---*/
	bool Dirichlet = config->GetDirichlet_Boundary(config->GetMarker_All_Tag(val_marker));
	if (Dirichlet) {
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			Point = geometry->vertex[val_marker][iVertex]->GetNode();
			Solution[0]= config->GetDirichlet_Value(config->GetMarker_All_Tag(val_marker));
			node[Point]->SetSolution(Solution);
			LinSysRes.SetBlock(Point, Solution);
			StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
		}
	}
}

/*!
 * \method BC_Sym_Plane
 * \brief Dirichlet/Neumann boundary condition
 * \author A. Lonkar
 */
void CElectricSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
		unsigned short val_marker) {
	unsigned long Point, iVertex;

	/*--- Identify if a boundary is Dirichlet or Neumman ---*/
	bool Dirichlet = config->GetDirichlet_Boundary(config->GetMarker_All_Tag(val_marker));
	if (Dirichlet) {
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			Point = geometry->vertex[val_marker][iVertex]->GetNode();
			Solution[0]= config->GetDirichlet_Value(config->GetMarker_All_Tag(val_marker));
			node[Point]->SetSolution(Solution);
			LinSysRes.SetBlock(Point, Solution);
			StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
		}
	}
}

/*!
 * \method BC_HeatFlux_Wall
 * \brief Dirichlet/Neumann boundary condition
 * \author A. Lonkar
 */
void CElectricSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
		unsigned short val_marker) {
	unsigned long Point, iVertex;

	/*--- Identify if a boundary is Dirichlet or Neumman ---*/
	bool Dirichlet = config->GetDirichlet_Boundary(config->GetMarker_All_Tag(val_marker));
	if (Dirichlet) {
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			Point = geometry->vertex[val_marker][iVertex]->GetNode();
			Solution[0]= config->GetDirichlet_Value(config->GetMarker_All_Tag(val_marker));
			node[Point]->SetSolution(Solution);
			LinSysRes.SetBlock(Point, Solution);
			StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
		}
	}
}

/*!
 * \method BC_Outlet
 * \brief Dirichlet/Neumann boundary condition
 * \author A. Lonkar
 */
void CElectricSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
		unsigned short val_marker) {
	unsigned long Point, iVertex;

	/*--- Identify if a boundary is Dirichlet or Neumman ---*/
	bool Dirichlet = config->GetDirichlet_Boundary(config->GetMarker_All_Tag(val_marker));
	if (Dirichlet) {
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			Point = geometry->vertex[val_marker][iVertex]->GetNode();
			Solution[0]= config->GetDirichlet_Value(config->GetMarker_All_Tag(val_marker));
			node[Point]->SetSolution(Solution);
			LinSysRes.SetBlock(Point, Solution);
			StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
		}
	}
}

/*!
 * \method BC_Inlet
 * \brief Dirichlet/Neumann boundary condition
 * \author A. Lonkar
 */
void CElectricSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
		unsigned short val_marker) {
	unsigned long Point, iVertex;

	/*--- Identify if a boundary is Dirichlet or Neumman ---*/
	bool Dirichlet = config->GetDirichlet_Boundary(config->GetMarker_All_Tag(val_marker));
	if (Dirichlet) {
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			Point = geometry->vertex[val_marker][iVertex]->GetNode();
			Solution[0]= config->GetDirichlet_Value(config->GetMarker_All_Tag(val_marker));
			node[Point]->SetSolution(Solution);
			LinSysRes.SetBlock(Point, Solution);
			StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
		}
	}
}

/*!
 * \method BC_Far_Field
 * \brief Dirichlet/Neumann boundary condition
 * \author A. Lonkar
 */
void CElectricSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
		unsigned short val_marker) {
	unsigned long Point, iVertex;

	/*--- Identify if a boundary is Dirichlet or Neumman ---*/
	bool Dirichlet = config->GetDirichlet_Boundary(config->GetMarker_All_Tag(val_marker));
	if (Dirichlet) {
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			Point = geometry->vertex[val_marker][iVertex]->GetNode();
			Solution[0]= config->GetDirichlet_Value(config->GetMarker_All_Tag(val_marker));
			node[Point]->SetSolution(Solution);
			LinSysRes.SetBlock(Point, Solution);
			StiffMatrix.DeleteValsRowi(Point); // & includes 1 in the diagonal
		}
	}
}
