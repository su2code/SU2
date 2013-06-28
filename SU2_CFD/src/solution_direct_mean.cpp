/*!
 * \file solution_direct_mean.cpp
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

CPotentialSolution::CPotentialSolution(void) : CSolution() { }

CPotentialSolution::CPotentialSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	
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
	
	Vector = new double[nDim];
	
	/*--- Point to point stiffness matrix ---*/
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
	
	/*--- Initialization of the structure of the whole Jacobian ---*/
	InitializeStiffMatrixStructure(geometry, config);
	xsol = new double [nPoint*nVar];
	rhs = new double [nPoint*nVar];
	
	/*--- Computation of gradients by least squares ---*/
	Smatrix = new double* [nDim];
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Smatrix[iDim] = new double [nDim];
	
	cvector = new double* [nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		cvector[iVar] = new double [nDim];
	
	
	for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
		node[iPoint] = new CPotentialVariable(0.0, nDim, nVar, config);
	
	Density_Inf = 1.0;	
	Energy_Inf = ((1.0/(Gamma*(Gamma-1)*config->GetMach_FreeStreamND()*config->GetMach_FreeStreamND() )) + 0.5);
	Velocity_Inf = new double [nDim];
	
	if (nDim == 2) {
		Velocity_Inf[0] = cos( (config->GetAoA()*PI_NUMBER) /180.0 );	
		Velocity_Inf[1] = sin( (config->GetAoA()*PI_NUMBER) /180.0 );
	}
	if (nDim == 3) {
		Velocity_Inf[0] = cos( (config->GetAoA()*PI_NUMBER) /180.0 ) * cos( (config->GetAoS()*PI_NUMBER) /180.0 );
		Velocity_Inf[1] = sin( (config->GetAoS()*PI_NUMBER) /180.0 );
		Velocity_Inf[2] = sin( (config->GetAoA()*PI_NUMBER) /180.0 ) * cos( (config->GetAoS()*PI_NUMBER) /180.0 );		
	}
	
	double Vel2 = 0;
	for (unsigned short iDim = 0; iDim < nDim; iDim++)
		Vel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
	Pressure_Inf = Gamma_Minus_One*Density_Inf*(Energy_Inf-0.5*Vel2);
	
}

CPotentialSolution::~CPotentialSolution(void) {
	
	delete [] Residual;
	delete [] Residual_Max;
	delete [] Residual_i;
	delete [] Residual_j;
	
	delete [] Solution;
	delete [] Solution_i;
	delete [] Solution_j;
	
	delete [] Vector_i;
	delete [] Vector_j;
	
	delete [] FlowSolution_i;
	delete [] FlowSolution_j;
	
	for (unsigned short iVar = 0; iVar < nVar; iVar++) {
		delete [] Jacobian_i[iVar];
		delete [] Jacobian_j[iVar];
	}
	delete [] Jacobian_i;
	delete [] Jacobian_j;
	delete [] xsol;
	delete [] rhs;
	
	for (unsigned short iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;
	
	for (unsigned short iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
}

void CPotentialSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		node[iPoint]->SetResidualZero();
	
	StiffMatrix.SetValZero();
}

void CPotentialSolution::Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
											unsigned short iMesh) {
	unsigned long iPoint;
	
	/*--- Build lineal system ---*/
	unsigned short var = 0;
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {		
		rhs[iPoint] = node[iPoint]->GetResidual(var);
		xsol[iPoint] = node[iPoint]->GetSolution(var);
	}
	
	/*--- Solve the system ---*/
	StiffMatrix.SGSSolution(rhs, xsol, 1e-6, 999999, true);
	
	SetRes_Max(0, 1e-6);
	
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetSolution(0,xsol[iPoint]);
	
}

void CPotentialSolution::Compute_Residual(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
										  unsigned short iMesh) { }

void CPotentialSolution::Galerkin_Method(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
										 CConfig *config, unsigned short iMesh) {
	
	unsigned long iElem, iPoint = 0, jPoint = 0, Point_2 = 0, Point_3 = 0;
	double *Coord_i = NULL, *Coord_j = NULL, *Coord_2 = NULL, *Coord_3 = NULL;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
		/*--- Points in edge ---*/
		iPoint = geometry->elem[iElem]->GetNode(0);
		jPoint = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		if (nDim == 3) Point_3 = geometry->elem[iElem]->GetNode(3);
		
		/*--- Points coordinates  ---*/
		Coord_i = geometry->node[iPoint]->GetCoord();
		Coord_j = geometry->node[jPoint]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		if (nDim == 3) Coord_3 = geometry->node[Point_3]->GetCoord();
		
		if (nDim == 2) solver->SetCoord(Coord_i, Coord_j, Coord_2);
		if (nDim == 3) solver->SetCoord(Coord_i, Coord_j, Coord_2, Coord_3);
		
		solver->SetResidual(StiffMatrix_Elem, config);
		
		if (nDim == 2) {
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]; StiffMatrix.AddBlock(iPoint,iPoint,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]; StiffMatrix.AddBlock(iPoint,jPoint,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]; StiffMatrix.AddBlock(iPoint,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]; StiffMatrix.AddBlock(jPoint,iPoint,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]; StiffMatrix.AddBlock(jPoint,jPoint,StiffMatrix_Node);		
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]; StiffMatrix.AddBlock(jPoint,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]; StiffMatrix.AddBlock(Point_2,iPoint,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]; StiffMatrix.AddBlock(Point_2,jPoint,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]; StiffMatrix.AddBlock(Point_2,Point_2,StiffMatrix_Node);			
		}
		
		if (nDim == 3) {
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]; StiffMatrix.AddBlock(iPoint,iPoint,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]; StiffMatrix.AddBlock(iPoint,jPoint,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]; StiffMatrix.AddBlock(iPoint,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][3]; StiffMatrix.AddBlock(iPoint,Point_3,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]; StiffMatrix.AddBlock(jPoint,iPoint,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]; StiffMatrix.AddBlock(jPoint,jPoint,StiffMatrix_Node);		
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]; StiffMatrix.AddBlock(jPoint,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][3]; StiffMatrix.AddBlock(jPoint,Point_3,StiffMatrix_Node);
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]; StiffMatrix.AddBlock(Point_2,iPoint,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]; StiffMatrix.AddBlock(Point_2,jPoint,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]; StiffMatrix.AddBlock(Point_2,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][3]; StiffMatrix.AddBlock(Point_2,Point_3,StiffMatrix_Node);
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][0]; StiffMatrix.AddBlock(Point_3,iPoint,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][1]; StiffMatrix.AddBlock(Point_3,jPoint,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][2]; StiffMatrix.AddBlock(Point_3,Point_2,StiffMatrix_Node);			
			StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][3]; StiffMatrix.AddBlock(Point_3,Point_3,StiffMatrix_Node);
		}
	}
}

void CPotentialSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
									   unsigned short val_marker) { }

void CPotentialSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
									  unsigned short val_marker) {
	unsigned long Point, iVertex;
	unsigned short iDim;
	double vel_proj;
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		Point = geometry->vertex[val_marker][iVertex]->GetNode();
		geometry->vertex[val_marker][iVertex]->GetNormal(Vector); // Vector --> Face_Normal
		vel_proj = 0; 
		for (iDim = 0; iDim < nDim; iDim++)
			vel_proj += GetDensity_Velocity_Inf(iDim)*Vector[iDim];
		
		Residual[0] = -0.5*vel_proj;
		node[Point]->AddResidual(Residual);
	}
}

void CPotentialSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
								  unsigned short val_marker) {
	unsigned long Point, iVertex;
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		Point = geometry->vertex[val_marker][iVertex]->GetNode();
		Solution[0]= 0.0;
		node[Point]->SetSolution_Old(Solution);
		StiffMatrix.DeleteValsRowi(Point);
	}
}

void CPotentialSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
								   unsigned short val_marker) {
	unsigned long Point, iVertex;
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		Point = geometry->vertex[val_marker][iVertex]->GetNode();
		Solution[0]= 0.0;
		node[Point]->SetSolution_Old(Solution);
		StiffMatrix.DeleteValsRowi(Point);
	}
}

void CPotentialSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
									  unsigned short val_marker) { }

CEulerSolution::CEulerSolution(void) : CSolution() { }

CEulerSolution::CEulerSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned long iPoint, index;
	unsigned short iVar, iDim, iMarker;
	string mesh_filename, text_line;
	ifstream restart_file;

	double AoA = (config->GetAoA()*PI_NUMBER) / 180.0;
	double AoS = (config->GetAoS()*PI_NUMBER) / 180.0;
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();
	bool shocktube = config->GetShockTube();
	
	/*--- Set the gamma value ---*/
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Define geometry constant in the solver structure ---*/
	nDim = geometry->GetnDim();	
	if (incompressible) nVar = nDim + 1;
	else nVar = nDim + 2;
	nMarker = config->GetnMarker_All();
	nPoint = geometry->GetnPoint();
	node = new CVariable*[geometry->GetnPoint()];
	Sum_Delta_Time = 0.0;
	
	/*--- Define some auxiliary vectors related to the residual ---*/
	Residual   = new double[nVar]; Residual_Max = new double[nVar];
	Residual_i = new double[nVar]; Residual_j   = new double[nVar];
	Res_Conv   = new double[nVar]; Res_Visc     = new double[nVar];
	
	/*--- Define some auxiliary vectors related to the solution ---*/
	Solution   = new double[nVar];
	Solution_i = new double[nVar]; Solution_j = new double[nVar];
	
	/*--- Define some auxiliary vectors related to the geometry ---*/
	Vector   = new double[nDim];
	Vector_i = new double[nDim]; Vector_j = new double[nDim];
	
	/*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED) {
		p1_Und_Lapl = new double [geometry->GetnPoint()]; 
		p2_Und_Lapl = new double [geometry->GetnPoint()]; 
  }
	
	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {

		/*--- Point to point Jacobians ---*/
		Jacobian_i = new double* [nVar]; 
		Jacobian_j = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new double [nVar]; 
			Jacobian_j[iVar] = new double [nVar]; 
		}
		
		/*--- Initialization of the structure for the whole Jacobian ---*/
		InitializeJacobianStructure(geometry, config);
		xsol = new double [geometry->GetnPoint()*nVar];
		rhs  = new double [geometry->GetnPoint()*nVar];
	}
	
	/*--- Computation of gradients by least squares ---*/
	if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) || 
			(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) {
		least_squares = true;
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim]; 
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
		
    /*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nVar+1]; 
		for (iVar = 0; iVar < nVar+1; iVar++)
			cvector[iVar] = new double [nDim];
	}
	
	/*--- Forces definition and coefficient in all the markers ---*/
	CPressure = new double* [config->GetnMarker_All()];
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		CPressure[iMarker] = new double [geometry->nVertex[iMarker]];
	
	ForceInviscid  = new double[nDim];
	MomentInviscid = new double[3];
	CDrag_Inv      = new double[config->GetnMarker_All()];
	CLift_Inv      = new double[config->GetnMarker_All()];
	CSideForce_Inv = new double[config->GetnMarker_All()];
	CPress_Inv     = new double[config->GetnMarker_All()];
	CMx_Inv        = new double[config->GetnMarker_All()];
	CMy_Inv        = new double[config->GetnMarker_All()];
	CMz_Inv        = new double[config->GetnMarker_All()];
	CEff_Inv       = new double[config->GetnMarker_All()];
	CEquivArea_Inv = new double[config->GetnMarker_All()];
	CNearFieldPress_Inv = new double[config->GetnMarker_All()];
	Total_CDrag = 0.0; Total_CLift = 0.0; Total_CSideForce = 0.0;
	Total_CPress = 0.0; Total_CMx = 0.0; Total_CMy = 0.0; Total_CMz = 0.0;
	Total_CEff = 0.0; Total_CEquivArea = 0.0; Total_CNearFieldPress = 0.0;
	
	/*--- Read farfield conditions from config ---*/
	Density_Inf  = config->GetDensity_FreeStreamND();
	Pressure_Inf = config->GetPressure_FreeStreamND();
	Velocity_Inf = config->GetVelocity_FreeStreamND();
	Energy_Inf   = config->GetEnergy_FreeStreamND();
	Mach_Inf     = config->GetMach_FreeStreamND();
	
	/*--- Inlet/Outlet boundary conditions, using infinity values ---*/
	Density_Inlet = Density_Inf;		Density_Outlet = Density_Inf;
	Pressure_Inlet = Pressure_Inf;	Pressure_Outlet = Pressure_Inf;
	Energy_Inlet = Energy_Inf;			Energy_Outlet = Energy_Inf;
	Mach_Inlet = Mach_Inf;					Mach_Outlet = Mach_Inf;
	Velocity_Inlet  = new double [nDim]; Velocity_Outlet = new double [nDim];
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_Inlet[iDim] = Velocity_Inf[iDim];
		Velocity_Outlet[iDim] = Velocity_Inf[iDim];
	}		
	
	/*--- Inlet/Outlet boundary conditions, using shock tube problem ---*/
	if (shocktube) {		
		double Temperature = 300.0;
		double Ru = 8314.462175; 
		double M1 = AVOGAD_CONSTANT*config->GetParticle_Mass(0);     
		double R1 = Ru/M1;
		Density_Inlet = 2.1331E-1;													Density_Outlet = 0.841949;
		Pressure_Inlet = Density_Inlet * R1 * Temperature;	Pressure_Outlet = 3.72271E6;
		Mach_Inlet = 14.95;																	Mach_Outlet =  0.4504;
		
		if (nDim == 2) {
			Velocity_Inlet[0]  = cos(AoA) * Mach_Inlet * sqrt(Gamma*Pressure_Inlet/Density_Inlet);	
			Velocity_Inlet[1]  = sin(AoA) * Mach_Inlet * sqrt(Gamma*Pressure_Inlet/Density_Inlet);
			Velocity_Outlet[0] = cos(AoA) * Mach_Outlet * sqrt(Gamma*Pressure_Outlet/Density_Outlet);	
			Velocity_Outlet[1] = sin(AoA) * Mach_Outlet * sqrt(Gamma*Pressure_Outlet/Density_Outlet);
		}
		if (nDim == 3) {
			Velocity_Inlet[0]  = cos(AoA) * cos(AoS) * Mach_Inlet * sqrt(Gamma*Pressure_Inlet/Density_Inlet);
			Velocity_Inlet[1]  = sin(AoS) * Mach_Inlet * sqrt(Gamma*Pressure_Inlet/Density_Inlet);
			Velocity_Inlet[2]  = sin(AoA) * cos(AoS) * Mach_Inlet * sqrt(Gamma*Pressure_Inlet/Density_Inlet);
			Velocity_Outlet[0] = cos(AoA) * cos(AoS) * Mach_Outlet * sqrt(Gamma*Pressure_Outlet/Density_Outlet);
			Velocity_Outlet[1] = sin(AoS) * Mach_Outlet * sqrt(Gamma*Pressure_Outlet/Density_Outlet);
			Velocity_Outlet[2] = sin(AoA) * cos(AoS) * Mach_Outlet * sqrt(Gamma*Pressure_Outlet/Density_Outlet);
		}
	
		Density_Inf  = Density_Inlet;
		Pressure_Inf = Pressure_Inlet;
		Velocity_Inf = Velocity_Inlet;
		Energy_Inf   = Energy_Inlet;
		Mach_Inf     = Mach_Inlet;
		
		cout << " Pressure_Inf = " << Pressure_Inf << endl;
		cout << " Density_Inf = " << Density_Inf << endl;
		cout << " Energy_Inf = " << Energy_Inf << endl;
		cout << " Mach_Inf = " << Mach_Inf << endl;
		cout << " Velocity_Inlet = " << Velocity_Inlet[0] << endl;	
		cout << " Pressure_Inlet = " << Pressure_Inlet << endl;
		cout << " Density_Inlet = " << Density_Inlet << endl;
		cout << " Energy_Inlet = " << Energy_Inlet << endl;
		cout << " Mach_Inlet = " << Mach_Inlet << endl;
		cout << " Pressure_Outlet = " << Pressure_Outlet << endl;
		cout << " Density_Outlet = " << Density_Outlet << endl;
		cout << " Energy_Outlet = " << Energy_Outlet << endl;
		cout << " Mach_Outlet = " << Mach_Outlet << endl;
		cout << " Velocity_Outlet = " << Velocity_Outlet[0] << endl;
	}
	
	/*--- For a rotating frame, set the velocity due to rotation
	 at each point just once, and store it in the geometry class. ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		if (rotating_frame) {
			double RotVel[3], distance[3];
			double *coord = geometry->node[iPoint]->GetCoord();
			double *axis  = config->GetRotAxisOrigin();
			double *omega = config->GetOmega_FreeStreamND();
			double Lref   = config->GetLength_Ref();
			
			/*--- Calculate non-dim distance fron rotation center ---*/
			distance[0] = (coord[0]-axis[0])/Lref;
			distance[1] = (coord[1]-axis[1])/Lref;
			distance[2] = (coord[2]-axis[2])/Lref;
			
			/*--- Calculate the angular velocity as omega X r ---*/
			RotVel[0] = omega[1]*(distance[2]) - omega[2]*(distance[1]);
			RotVel[1] = omega[2]*(distance[0]) - omega[0]*(distance[2]);
			RotVel[2] = omega[0]*(distance[1]) - omega[1]*(distance[0]);
			
			geometry->node[iPoint]->SetRotVel(RotVel);
		}
	}
	
	if (!restart) {
		/*--- Restart the solution from infinity ---*/
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			if (!shocktube) {
				/*--- Initialization using value at the infinity (normal case) ---*/
				node[iPoint] = new CEulerVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
			}
			else {
				/*--- Initialization of the shock tube problem, jump in the middle of the tube ---*/
				if(geometry->node[iPoint]->GetCoord(0) <= 0.0) 
					node[iPoint] = new CEulerVariable(Density_Inlet, Velocity_Inlet, Energy_Inlet, nDim, nVar, config);
				else
					node[iPoint] = new CEulerVariable(Density_Outlet, Velocity_Outlet, Energy_Outlet, nDim, nVar, config);
			}
		}
	}
	else {
		
		/*--- Restart the solution from file information ---*/
		mesh_filename = config->GetSolution_FlowFileName();
		restart_file.open(mesh_filename.data(), ios::in);
		
		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
			cout << "There is no flow restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get(); exit(1);
		}
		
		/*--- Read the restart file ---*/
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			if (incompressible) {
				if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2];
				if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
			}
			else {
				if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
				if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
			}
			node[iPoint] = new CEulerVariable(Solution, nDim, nVar, config);
		}
		
		/*--- Close the restart file ---*/
		restart_file.close();

	}
	
	/*--- Define solver parameters needed for execution of destructor ---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED) space_centered = true;
	else space_centered = false;
	if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
	else euler_implicit = false;
	if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||	
			(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) least_squares = true;
	else least_squares = false;
}

CEulerSolution::~CEulerSolution(void) {
	unsigned short iVar, iDim, iMarker, iPoint;
	
	delete [] Residual; delete [] Residual_Max;
	delete [] Residual_i; delete [] Residual_j;
	delete [] Res_Conv; delete [] Res_Visc;
	delete [] Solution; delete [] Solution_i; delete [] Solution_j;
	delete [] Vector_i; delete [] Vector_j;
	
	if (space_centered) {
		delete [] p1_Und_Lapl;
		delete [] p2_Und_Lapl;
	}
	
	delete [] ForceInviscid; delete [] MomentInviscid;
	delete [] CDrag_Inv; delete [] CLift_Inv; delete [] CSideForce_Inv;
	delete [] CMx_Inv; delete [] CMy_Inv; delete [] CMz_Inv;
	delete [] CPress_Inv; delete [] CEff_Inv; delete [] CEquivArea_Inv; delete [] CNearFieldPress_Inv;
	
	if (euler_implicit){
		for (iVar = 0; iVar < nVar; iVar++) {
			delete [] Jacobian_i[iVar];
			delete [] Jacobian_j[iVar];
		}
		delete [] Jacobian_i;
		delete [] Jacobian_j;
		delete [] xsol;
		delete [] rhs;
	}
	
	if (least_squares){
		for (iDim = 0; iDim < nDim; iDim++)
			delete [] Smatrix[iDim];
		delete [] Smatrix;
		
		for (iVar = 0; iVar < nVar+1; iVar++)
			delete [] cvector[iVar];
		delete [] cvector;
	}
	
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		delete [] CPressure[iMarker];
	delete [] CPressure;
	
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		delete node[iPoint];
	delete [] node;
}

void CEulerSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
	bool twophase = (config->GetKind_Solver() == TWO_PHASE_FLOW);
	double Gas_Constant = config->GetGas_Constant();

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {

		/*--- Compute squared velocity, sound velocity, pressure, and enthalpy ---*/
		if (incompressible) {
			node[iPoint]->SetPressureInc();
			node[iPoint]->SetVelocityInc2();
			node[iPoint]->SetBetaInc2(config);
			if (!twophase) node[iPoint]->SetDensityInc(Density_Inf);
		}
		else {
			node[iPoint]->SetVelocity2();
			node[iPoint]->SetPressure(Gamma, geometry->node[iPoint]->GetCoord());
			node[iPoint]->SetTemperature(Gas_Constant);
			node[iPoint]->SetSoundSpeed(Gamma);
			node[iPoint]->SetEnthalpy();
		}
		
		/*--- Initialize the convective residual vector ---*/
		node[iPoint]->Set_ResConv_Zero();

		/*--- Initialize the viscous residual vector ---*/
		if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit) 
			node[iPoint]->Set_ResVisc_Zero();
	}
	
	/*--- Initialize the jacobian matrices ---*/
	if (implicit) Jacobian.SetValZero();
}

void CEulerSolution::SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {
	double *Normal, Area, dV, Mean_SoundSpeed, Mean_ProjVel, Mean_BetaInc2, Lambda, Local_Delta_Time, 
	Global_Delta_Time = 1E6;
  double ProjVel, ProjVel_i, ProjVel_j;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;
	
  bool rotating_frame = config->GetRotating_Frame();
  bool incompressible = config->GetIncompressible();
  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
	
	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		node[iPoint]->SetMax_Lambda_Inv(0.0);
		if (incompressible) {
			node[iPoint]->SetVelocityInc2();
			node[iPoint]->SetBetaInc2(config);
		}
		else {
			node[iPoint]->SetVelocity2();
			node[iPoint]->SetPressure(Gamma, geometry->node[iPoint]->GetCoord());
			node[iPoint]->SetSoundSpeed(Gamma);
		}
	}
	
	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		
		/*--- Point identification, Normal vector and area ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); 
		jPoint = geometry->edge[iEdge]->GetNode(1);
		Normal = geometry->edge[iEdge]->GetNormal();
		Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
		
		/*--- Mean Values ---*/
		if (incompressible) {
			Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVelInc(Normal) + node[jPoint]->GetProjVelInc(Normal));
			Mean_BetaInc2 = 0.5 * (node[iPoint]->GetBetaInc2() + node[jPoint]->GetBetaInc2());
			Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + Mean_BetaInc2*Area*Area);
		}
		else {
			Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
			Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
		}
		
		/*--- Contribution from rotating framework ---*/
    if (rotating_frame) {
			double *RotVel_i = geometry->node[iPoint]->GetRotVel();
			double *RotVel_j = geometry->node[jPoint]->GetRotVel();
			ProjVel_i = 0.0; ProjVel_j =0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				ProjVel_i += RotVel_i[iDim]*Normal[iDim];
				ProjVel_j += RotVel_j[iDim]*Normal[iDim];
			}
			Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
		}
    
		/*--- Inviscid contribution ---*/
		Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed ;
		node[iPoint]->AddMax_Lambda_Inv(Lambda);
		node[jPoint]->AddMax_Lambda_Inv(Lambda);
	}
	
	/*--- Loop boundary edges ---*/
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) { 
		for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			
			/*--- Point identification, Normal vector and area ---*/
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

			/*--- Mean Values ---*/
			if (incompressible) {
				Mean_ProjVel = node[iPoint]->GetProjVelInc(Normal);
				Mean_BetaInc2 = node[iPoint]->GetBetaInc2();
				Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + Mean_BetaInc2*Area*Area); 
			}
			else {
				Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
				Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;		
			}
			
			/*--- Contribution from rotating framework ---*/
      if (rotating_frame) {
				double *RotVel = geometry->node[iPoint]->GetRotVel();
				ProjVel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					ProjVel += RotVel[iDim]*Normal[iDim];
				Mean_ProjVel -= ProjVel;
			}
      
			/*--- Inviscid contribution ---*/
			Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
			node[iPoint]->AddMax_Lambda_Inv(Lambda);
		}
	}
	
	/*--- Each element uses their own speed ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		dV = geometry->node[iPoint]->GetVolume();
		Local_Delta_Time = config->GetCFL(iMesh)*dV / node[iPoint]->GetMax_Lambda_Inv();
		Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
		Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
		Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
/*		if (incompressible) {
			double DensictyInc = node[iPoint]->GetDensityInc();
			DensictyInc = max(DensictyInc, config->GetRatioDensity());
			node[iPoint]->SetDelta_Time(Local_Delta_Time*DensictyInc);
		}
		else {*/
			node[iPoint]->SetDelta_Time(Local_Delta_Time);
//		}
	}
	
	Sum_Delta_Time += Min_Delta_Time;
	
	/*--- For exact time solution use the minimum delta time of the whole mesh ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING)
		for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
				node[iPoint]->SetDelta_Time(Global_Delta_Time);
}

void CEulerSolution::Centred_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
										  CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool dissipation = ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit);
	bool high_order_diss = ((config->GetKind_Centred() == JST) && (iMesh == MESH_0));
	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();

	/*--- Artificial dissipation preprocessing ---*/
	if (dissipation) {
		SetSpectral_Radius(geometry, config);
		if (high_order_diss) {
			if (!incompressible) SetPress_Switch(geometry); 
			SetUndivided_Laplacian(geometry, config); 
		}
	}
	
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		/*--- Points in edge, Set normal vectors, and number of neighbors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());
		solver->SetNeighbor(geometry->node[iPoint]->GetnPoint(), geometry->node[jPoint]->GetnPoint());

		/*--- Set conservative variables w/o reconstruction ---*/
		solver->SetConservative(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
		
		/*--- Set pressure, soundSpeed, and enthalpy variable w/o reconstruction ---*/
		if (incompressible) {
			solver->SetDensityInc(node[iPoint]->GetDensityInc(), node[jPoint]->GetDensityInc());
			solver->SetBetaInc2(node[iPoint]->GetBetaInc2(), node[jPoint]->GetBetaInc2());
			double *iCoord = NULL, *jCoord= NULL;
			iCoord = geometry->node[iPoint]->GetCoord();
			jCoord = geometry->node[jPoint]->GetCoord();
			solver->SetCoord(iCoord, jCoord);
		}
		else {
			solver->SetPressure(node[iPoint]->GetPressure(), node[jPoint]->GetPressure());
			solver->SetSoundSpeed(node[iPoint]->GetSoundSpeed(), node[jPoint]->GetSoundSpeed());
			solver->SetEnthalpy(node[iPoint]->GetEnthalpy(), node[jPoint]->GetEnthalpy());
		}
		
		/*--- Set lambda, undivided laplacian and sensor variables ---*/
		solver->SetLambda(node[iPoint]->GetLambda(), node[jPoint]->GetLambda());
		
		/*--- Undivided laplacian and sensor ---*/
		if (dissipation && high_order_diss) {
			solver->SetUndivided_Laplacian(node[iPoint]->GetUnd_Lapl(), node[jPoint]->GetUnd_Lapl());
			if (!incompressible)
				solver->SetSensor(node[iPoint]->GetSensor(), node[jPoint]->GetSensor()); 
		}
		
		/*--- Rotational frame ---*/
		if (rotating_frame)
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());
		
		/*--- Compute residuals ---*/
		solver->SetResidual(Res_Conv, Res_Visc, Jacobian_i, Jacobian_j, dissipation, config);
		
		/*--- Update convective and artificial dissipation residuals ---*/
		node[iPoint]->AddRes_Conv(Res_Conv);
		node[jPoint]->SubtractRes_Conv(Res_Conv);
		if (dissipation) {
			node[iPoint]->AddRes_Visc(Res_Visc);
			node[jPoint]->SubtractRes_Visc(Res_Visc); }

		/*--- Set implicit stuff ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
			Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
			Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j); }
	}
}

void CEulerSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
									 CConfig *config, unsigned short iMesh) {
	double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, *Limiter_i = NULL, *Limiter_j = NULL, *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool high_order_diss = (((config->GetKind_Upwind() == ROE_2ND) || (config->GetKind_Upwind() == AUSM_2ND) 
													 || (config->GetKind_Upwind() == HLLC_2ND)) && (iMesh == MESH_0));
	bool incompressible = config->GetIncompressible();

	if (high_order_diss) { 
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry); 
		if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) || 
			(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) SetSolution_Gradient_LS(geometry, config);
		if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN) SetSolution_Limiter(geometry, config);
	}
	
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());
		
		/*--- Set conservative variables w/o reconstruction ---*/
		solver->SetConservative(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
		
		/*--- Set the density and beta w/o reconstruction (incompressible flows) ---*/
		if (incompressible) {
			solver->SetDensityInc(node[iPoint]->GetDensityInc(), node[jPoint]->GetDensityInc());
			solver->SetBetaInc2(node[iPoint]->GetBetaInc2(), node[jPoint]->GetBetaInc2());
			double *iCoord = NULL, *jCoord= NULL;
			iCoord = geometry->node[iPoint]->GetCoord();
			jCoord = geometry->node[jPoint]->GetCoord();
			solver->SetCoord(iCoord, jCoord);
		}
		
		if (high_order_diss) {
			U_i = node[iPoint]->GetSolution();
			U_j = node[jPoint]->GetSolution();

			for (iDim = 0; iDim < nDim; iDim++) {
				Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
				Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}
			Gradient_i = node[iPoint]->GetGradient(); 
			Gradient_j = node[jPoint]->GetGradient(); 
			if (config->GetKind_SlopeLimit_Flow() != NONE) {
				Limiter_j = node[jPoint]->GetLimiter();
				Limiter_i = node[iPoint]->GetLimiter();				
			}
			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				if (config->GetKind_SlopeLimit_Flow() == NONE) {
					Solution_i[iVar] = U_i[iVar] + Project_Grad_i;
					Solution_j[iVar] = U_j[iVar] + Project_Grad_j;
				}
				else {
					Solution_i[iVar] = U_i[iVar] + Project_Grad_i*Limiter_i[iVar];
					Solution_j[iVar] = U_j[iVar] + Project_Grad_j*Limiter_j[iVar];
				}
			}
			/*--- Set conservative variables with reconstruction ---*/
			solver->SetConservative(Solution_i, Solution_j);
		}

		/*--- Compute the residual ---*/
		solver->SetResidual(Res_Conv, Jacobian_i, Jacobian_j, config); 
		
		U_i = node[iPoint]->GetSolution();
		U_j = node[jPoint]->GetSolution();
		
		/*--- Update residual value ---*/
		node[iPoint]->AddRes_Conv(Res_Conv);
		node[jPoint]->SubtractRes_Conv(Res_Conv);
		
		/*--- Set implicit jacobians ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
			Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
			Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);
		}
	}
	
}

void CEulerSolution::SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
												CConfig *config, unsigned short iMesh) {
	unsigned short iVar;
	unsigned long iPoint;
	
	bool gravity = config->GetGravityForce();
	bool rotating_frame = config->GetRotating_Frame();

	if (rotating_frame) {
		for (iVar = 0; iVar < nVar; iVar++)
			Residual[iVar] = 0;
		
		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 
			
			/*--- Set solution  ---*/
			solver->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
			
			/*--- Set control volume ---*/
			solver->SetVolume(geometry->node[iPoint]->GetVolume());
			
			/*--- Set rotational velocity ---*/
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
			
			/*--- Compute Residual ---*/
			solver->SetResidual(Residual, config);
			
			/*--- Add Residual ---*/
			node[iPoint]->AddRes_Conv(Residual);
		}
	}
	
	if (gravity) {
		for (iVar = 0; iVar < nVar; iVar++)
			Residual[iVar] = 0;
		
		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 
			
			/*--- Set solution  ---*/
			solver->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
			
			/*--- Set control volume ---*/
			solver->SetVolume(geometry->node[iPoint]->GetVolume());
			
			/*--- Compute Residual ---*/
			solver->SetResidual(Residual, config);
			
			/*--- Add Residual ---*/
			node[iPoint]->AddRes_Conv(Residual);
		}
	}
	
}

void CEulerSolution::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint, iEdge;
	double Pressure_i = 0, Pressure_j = 0;
	unsigned short iVar;
	bool boundary_i, boundary_j;
	double *Diff = new double[nVar];

	bool incompressible = config->GetIncompressible();
	
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetUnd_LaplZero();
	
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		
		/*--- Solution differences ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);
		
		/*--- Correction for compressible flows ---*/
		if (!incompressible) {
			Pressure_i = node[iPoint]->GetPressure(); Pressure_j = node[jPoint]->GetPressure();
			Diff[nVar-1] = (node[iPoint]->GetSolution(nVar-1) + Pressure_i) - (node[jPoint]->GetSolution(nVar-1) + Pressure_j);
		}
		
		boundary_i = geometry->node[iPoint]->GetBoundary_Physical();
		boundary_j = geometry->node[jPoint]->GetBoundary_Physical();
			
		/*--- Both points inside Omega ---*/
		if (!boundary_i && !boundary_j) {
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
		}
		
		/*--- iPoint inside Omega, jPoint on the boundary ---*/
		if (!boundary_i && boundary_j)
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
		
		/*--- jPoint inside Omega, iPoint on the boundary ---*/
		if (boundary_i && !boundary_j)
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
		
		/*--- Both points on the boundary ---*/
		if (boundary_i && boundary_j) {
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
		}
	}
	
	
	/*--- Prof. Jameson's idea for the conservation equation ---*/
	if (incompressible) {
		double *PressureLaplacian;
		PressureLaplacian = new double [geometry->GetnPoint()];
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			PressureLaplacian[iPoint] = 0.0;

		SetPressureLaplacian(geometry, PressureLaplacian);
		
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			node[iPoint]->SetUnd_Lapl(0, PressureLaplacian[iPoint]);
		}
		
		delete[] PressureLaplacian;
	}
	
	
	delete [] Diff;
}

void CEulerSolution::SetSpectral_Radius(CGeometry *geometry, CConfig *config) {
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker;
	double *Normal, Area, ProjVel, Lambda, SoundSpeed, BetaInc2;
	
	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
			node[iPoint]->SetLambda(0.0);
	
	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		
		Normal = geometry->edge[iEdge]->GetNormal();
		Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
		
		/*--- Inviscid contribution to the Point i ---*/
		if (incompressible) {
			BetaInc2 = node[iPoint]->GetBetaInc2();
			ProjVel = node[iPoint]->GetProjVelInc(Normal);
			SoundSpeed = sqrt(ProjVel*ProjVel + BetaInc2*Area*Area); 
		}
		else {
			SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
			ProjVel = node[iPoint]->GetProjVel(Normal);
		}
			
		/*--- Contribution from rotational frame ---*/
		if (rotating_frame) {
			double *RotVel = geometry->node[iPoint]->GetRotVel();
			double ProjRotVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjRotVel += RotVel[iDim]*Normal[iDim];
			ProjVel -= ProjRotVel;
		}
		
		Lambda = 0.5*(fabs(ProjVel) + SoundSpeed);
		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);
		
		/*--- Inviscid contribution to the Point j ---*/
		if (incompressible) {
			BetaInc2 = node[jPoint]->GetBetaInc2();
			ProjVel = node[jPoint]->GetProjVelInc(Normal);
			SoundSpeed = sqrt(ProjVel*ProjVel + BetaInc2*Area*Area); 
		}
		else {
			SoundSpeed = node[jPoint]->GetSoundSpeed() * Area;
			ProjVel = node[jPoint]->GetProjVel(Normal);
		}
    
		/*--- Contribution from rotational frame ---*/
		if (rotating_frame) {
			double *RotVel = geometry->node[jPoint]->GetRotVel();
			double ProjRotVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjRotVel += RotVel[iDim]*Normal[iDim];
			ProjVel -= ProjRotVel;
		}
		
		Lambda = 0.5*(fabs(ProjVel) + SoundSpeed);
		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);
	}
	
	/*--- Loop boundary edges ---*/
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			if (geometry->node[iPoint]->GetDomain()) {
				
				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

				/*--- Inviscid contribution to the iPoint ---*/
				if (incompressible) {
					BetaInc2 = node[iPoint]->GetBetaInc2();
					ProjVel = node[iPoint]->GetProjVelInc(Normal);
					SoundSpeed = sqrt(ProjVel*ProjVel + BetaInc2*Area*Area); 
				}
				else {
					SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
					ProjVel = node[iPoint]->GetProjVel(Normal);
				}
				
				/*--- Contribution from rotational frame ---*/
				if (rotating_frame) {
					double *RotVel = geometry->node[iPoint]->GetRotVel();
					double ProjRotVel = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						ProjRotVel += RotVel[iDim]*Normal[iDim];
					ProjVel -= ProjRotVel;
				}
				
				Lambda = fabs(ProjVel) + SoundSpeed;
				node[iPoint]->AddLambda(Lambda);
			}
		}
}

void CEulerSolution::SetPress_Switch(CGeometry *geometry) {
	unsigned long iEdge, iPoint, jPoint;
	double Pressure_i, Pressure_j;
	
	/*--- Reset variables to store the undivided pressure ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		p1_Und_Lapl[iPoint] = 0.0;
		p2_Und_Lapl[iPoint] = 0.0;
	}
	
	/*--- Evaluate the pressure sensor ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		
		Pressure_i = node[iPoint]->GetPressure();
		Pressure_j = node[jPoint]->GetPressure();
		
		if (geometry->node[iPoint]->GetDomain()) p1_Und_Lapl[iPoint] += (Pressure_j - Pressure_i);
		if (geometry->node[jPoint]->GetDomain()) p1_Und_Lapl[jPoint] += (Pressure_i - Pressure_j);
		
		if (geometry->node[iPoint]->GetDomain()) p2_Und_Lapl[iPoint] += (Pressure_i + Pressure_j);
		if (geometry->node[jPoint]->GetDomain()) p2_Und_Lapl[jPoint] += (Pressure_i + Pressure_j);
		
	}
	
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetSensor(fabs(p1_Und_Lapl[iPoint]) / p2_Und_Lapl[iPoint]);
}

void CEulerSolution::Inviscid_Forces(CGeometry *geometry, CConfig *config) {
	unsigned long iVertex, iPoint;
	unsigned short iDim, iMarker, Boundary, Monitoring;
	double Pressure, *Normal = NULL, dist[3], *Coord, Face_Area, PressInviscid;
	double factor, NFPressOF, RefVel2, RefDensity, RefPressure;
	bool rotating_frame = config->GetRotating_Frame();
  
	double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
	double Beta            = config->GetAoS()*PI_NUMBER/180.0;
	double RefAreaCoeff    = config->GetRefAreaCoeff();
	double RefLengthMoment = config->GetRefLengthMoment();
	double *Origin         = config->GetRefOriginMoment();
	double Pressure_Inf    = config->GetPressure_FreeStreamND();
    
  /*--- If this is a rotating frame problem, use the specified speed
        for computing the force coefficients. Otherwise, use the 
        freestream values which is the standard convention. ---*/
  if (rotating_frame) {
    RefVel2     = config->GetRotVel_Ref()*config->GetRotVel_Ref();
    RefDensity  = config->GetDensity_FreeStreamND();
    RefPressure = config->GetPressure_FreeStreamND();
  } else {
    double *Velocity_Inf = config->GetVelocity_FreeStreamND();
    RefVel2 = 0.0; 
    for (iDim = 0; iDim < nDim; iDim++) 
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    RefDensity  = config->GetDensity_FreeStreamND();
    RefPressure = config->GetPressure_FreeStreamND();
  }
  
	/*-- Initialization ---*/
	Total_CDrag = 0.0; Total_CLift = 0.0; Total_CSideForce = 0.0; 
	Total_CMx = 0.0; Total_CMy = 0.0; Total_CMz = 0.0;
	Total_CPress = 0.0; Total_CEff = 0.0;
	AllBound_CDrag_Inv = 0.0; AllBound_CLift_Inv = 0.0; AllBound_CSideForce_Inv = 0.0; 
	AllBound_CMx_Inv = 0.0; AllBound_CMy_Inv = 0.0; AllBound_CMz_Inv = 0.0;
	AllBound_CPress_Inv = 0.0; AllBound_CEff_Inv = 0.0; 
	AllBound_CNearFieldPress_Inv = 0.0;
	
	/*--- Loop over the Euler and Navier-Stokes markers ---*/
	factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary   = config->GetMarker_All_Boundary(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);
		
		if ((Boundary == EULER_WALL) || (Boundary == NO_SLIP_WALL) || (Boundary == NEARFIELD_BOUNDARY)) {
			
			for (iDim = 0; iDim < nDim; iDim++) ForceInviscid[iDim] = 0.0;
			MomentInviscid[0] = 0.0; MomentInviscid[1] = 0.0; MomentInviscid[2] = 0.0;
			NFPressOF = 0.0; PressInviscid = 0.0;

			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Pressure = node[iPoint]->GetPressure();
				CPressure[iMarker][iVertex] = (Pressure - RefPressure)*factor*RefAreaCoeff;
				
				/*--- Note that the pressure coefficient is computed at the 
				 halo cells (for visualization purposes), but not the forces ---*/
				if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					Coord = geometry->node[iPoint]->GetCoord();

					/*--- Quadratic objective function for the near field.
           This uses the infinity pressure regardless of Mach number. ---*/
					NFPressOF += 0.5*pow((Pressure - Pressure_Inf), 2.0)*Normal[nDim-1];

					Face_Area = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						/*--- Total force, distance computation, and face area ---*/
						ForceInviscid[iDim] -= Pressure*Normal[iDim]*factor;
						dist[iDim] = Coord[iDim] - Origin[iDim];
						Face_Area += Normal[iDim]*Normal[iDim];
					}
					Face_Area = sqrt(Face_Area);
					PressInviscid += CPressure[iMarker][iVertex]*Face_Area;
					
					/*--- Moment with respect to the reference axis ---*/
					if (iDim == 3) {
						MomentInviscid[0] -= Pressure*(Normal[2]*dist[1]-Normal[1]*dist[2])*factor/RefLengthMoment;
						MomentInviscid[1] -= Pressure*(Normal[0]*dist[2]-Normal[2]*dist[0])*factor/RefLengthMoment;
					}
					MomentInviscid[2] -= Pressure*(Normal[1]*dist[0]-Normal[0]*dist[1])*factor/RefLengthMoment;
				}
			}
			
			/*--- Transform ForceInviscid and MomentInviscid into non-dimensional coefficient ---*/
			if  (Monitoring == YES) {
				if (nDim == 2) {
					if (Boundary != NEARFIELD_BOUNDARY) {
						CDrag_Inv[iMarker] =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
						CLift_Inv[iMarker] = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
						CSideForce_Inv[iMarker] = 0.0;
						CPress_Inv[iMarker] = PressInviscid;
						CMx_Inv[iMarker] = 0.0;
						CMy_Inv[iMarker] = 0.0;
						CMz_Inv[iMarker] = MomentInviscid[2];
						CEff_Inv[iMarker] = CLift_Inv[iMarker]/CDrag_Inv[iMarker];
						CNearFieldPress_Inv[iMarker] = 0.0;
					}
					else {
						CDrag_Inv[iMarker] = 0.0; CLift_Inv[iMarker] = 0.0; CSideForce_Inv[iMarker] = 0.0;
						CPress_Inv[iMarker] = 0.0;
						CMx_Inv[iMarker] = 0.0; CMy_Inv[iMarker] = 0.0; CMz_Inv[iMarker] = 0.0;
						CEff_Inv[iMarker] = 0.0;
						CNearFieldPress_Inv[iMarker] = NFPressOF;
					}
				}
				if (nDim == 3) {
					if (Boundary != NEARFIELD_BOUNDARY) {
						CDrag_Inv[iMarker] =  ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta);
						CLift_Inv[iMarker] = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
						CSideForce_Inv[iMarker] = -ForceInviscid[0]*sin(Beta)*cos(Alpha) + ForceInviscid[1]*cos(Beta) - ForceInviscid[2]*sin(Beta)*sin(Alpha);
						CPress_Inv[iMarker] = PressInviscid;
						CMx_Inv[iMarker] = MomentInviscid[0];
						CMy_Inv[iMarker] = MomentInviscid[1];
						CMz_Inv[iMarker] = MomentInviscid[2];
						CEff_Inv[iMarker] = CLift_Inv[iMarker]/CDrag_Inv[iMarker];
						CNearFieldPress_Inv[iMarker] = 0.0;
					}
					else {
						CDrag_Inv[iMarker] = 0.0; CLift_Inv[iMarker] = 0.0; CSideForce_Inv[iMarker] = 0.0;
						CPress_Inv[iMarker] = 0.0;
						CMx_Inv[iMarker] = 0.0; CMy_Inv[iMarker] = 0.0; CMz_Inv[iMarker] = 0.0;
						CEff_Inv[iMarker] = 0.0;
						CNearFieldPress_Inv[iMarker] = NFPressOF;
					}
				}
				
				AllBound_CDrag_Inv += CDrag_Inv[iMarker];
				AllBound_CLift_Inv += CLift_Inv[iMarker];
				AllBound_CSideForce_Inv += CSideForce_Inv[iMarker];
				AllBound_CPress_Inv += CPress_Inv[iMarker];
				AllBound_CMx_Inv += CMx_Inv[iMarker];
				AllBound_CMy_Inv += CMy_Inv[iMarker];
				AllBound_CMz_Inv += CMz_Inv[iMarker];
				AllBound_CEff_Inv += CEff_Inv[iMarker];
				AllBound_CNearFieldPress_Inv += CNearFieldPress_Inv[iMarker];
			}
		}
	}
	
	Total_CDrag += AllBound_CDrag_Inv;
	Total_CLift += AllBound_CLift_Inv;
	Total_CSideForce += AllBound_CSideForce_Inv;
	Total_CPress += AllBound_CPress_Inv;
	Total_CMx += AllBound_CMx_Inv;
	Total_CMy += AllBound_CMy_Inv;
	Total_CMz += AllBound_CMz_Inv;
	Total_CEff += AllBound_CEff_Inv;
	Total_CNearFieldPress += AllBound_CNearFieldPress_Inv;
}

void CEulerSolution::RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, 
																					CConfig *config, unsigned short iRKStep) {
	double *Residual, *TruncationError, Vol, Delta;
	unsigned short iVar;
	unsigned long iPoint;
	double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
	
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);
	
	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Delta = node[iPoint]->GetDelta_Time() / Vol;
		TruncationError = node[iPoint]->GetTruncationError();
		Residual = node[iPoint]->GetResidual();
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, -(Residual[iVar]+TruncationError[iVar])*Delta*RK_AlphaCoeff);
			AddRes_Max( iVar, Residual[iVar]*Residual[iVar]*Vol );
		}
	}
	
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, sqrt(GetRes_Max(iVar)) );
}

void CEulerSolution::ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	double *Residual, *TruncationError, Vol, Delta;
	unsigned short iVar;
	unsigned long iPoint;
	
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);
	
	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Delta = node[iPoint]->GetDelta_Time() / Vol;
		TruncationError = node[iPoint]->GetTruncationError();
		Residual = node[iPoint]->GetResidual();
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, -(Residual[iVar]+TruncationError[iVar])*Delta);
			AddRes_Max( iVar, Residual[iVar]*Residual[iVar]*Vol );	
		}
	}
	
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, sqrt(GetRes_Max(iVar)) );
}

void CEulerSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta, Res, *local_ResConv, *local_ResVisc, *local_TruncationError, Vol;
	
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);
	
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			local_TruncationError = node[iPoint]->GetTruncationError();
			local_ResConv = node[iPoint]->GetResConv();
			local_ResVisc = node[iPoint]->GetResVisc();
			Vol = geometry->node[iPoint]->GetVolume();
			
			/*--- Modify matrix diagonal to assure diagonal dominance ---*/
			Delta = geometry->node[iPoint]->GetVolume()/node[iPoint]->GetDelta_Time();
		
			Jacobian.AddVal2Diag(iPoint,Delta);
			
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
				Res = local_ResConv[iVar]+local_ResVisc[iVar];
				rhs[total_index] = -(Res+local_TruncationError[iVar]);
				xsol[total_index] = 0.0;
				AddRes_Max( iVar, Res*Res*Vol );
			}
		}
	
	/*--- Solve the system ---*/
//	Jacobian.SGSSolution(rhs, xsol, 1e-9, 10, false);
	Jacobian.LU_SGSIteration(rhs, xsol);

		
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);	
		}
	}
	
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, sqrt(GetRes_Max(iVar)));
	
}

void CEulerSolution::SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint, iEdge, iVertex;
	unsigned short iDim, iVar, iMarker;
	double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average, 
	Partial_Gradient, Partial_Res, *Normal;
	
	bool incompressible = config->GetIncompressible();

	PrimVar_Vertex = new double [nVar+1];
	PrimVar_i = new double [nVar+1];
	PrimVar_j = new double [nVar+1]; 
	
	/*--- Set Gradient_Primitive to Zero (nVar+1) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetGradient_PrimitiveZero();

	/*--- Loop interior edges ---*/ 
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		
		if (incompressible) {
			PrimVar_i[0] = node[iPoint]->GetSolution(0);
			for (iDim = 0; iDim < nDim; iDim++)
				PrimVar_i[iDim+1] = node[iPoint]->GetSolution(iDim+1);
			PrimVar_i[nDim+1] = 0.0;
			PrimVar_j[0] = node[jPoint]->GetSolution(0);
			for (iDim = 0; iDim < nDim; iDim++)
				PrimVar_j[iDim+1] = node[jPoint]->GetSolution(iDim+1);
			PrimVar_j[nDim+1] = 0.0;
		}
		else {
			PrimVar_i[0] = node[iPoint]->GetTemperature();
			for (iDim = 0; iDim < nDim; iDim++)
				PrimVar_i[iDim+1] = node[iPoint]->GetVelocity(iDim);
			PrimVar_i[nDim+1] = node[iPoint]->GetPressure();
			PrimVar_i[nDim+2] = node[iPoint]->GetDensity();
			
			PrimVar_j[0] = node[jPoint]->GetTemperature();
			for (iDim = 0; iDim < nDim; iDim++)
				PrimVar_j[iDim+1] = node[jPoint]->GetVelocity(iDim);
			PrimVar_j[nDim+1] = node[jPoint]->GetPressure();
			PrimVar_j[nDim+2] = node[jPoint]->GetDensity();
		}

		Normal = geometry->edge[iEdge]->GetNormal();		
		for (iVar = 0; iVar < nVar+1; iVar++) {
			PrimVar_Average =  0.5 * ( PrimVar_i[iVar] + PrimVar_j[iVar] );
			for (iDim = 0; iDim < nDim; iDim++) {
				Partial_Res = PrimVar_Average*Normal[iDim];
				if (geometry->node[iPoint]->GetDomain())
					node[iPoint]->AddGradient_Primitive(iVar, iDim, Partial_Res);
				if (geometry->node[jPoint]->GetDomain())
					node[jPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
			}				
		}
	}
	
	/*--- Loop boundary edges ---*/
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
		for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			if (geometry->node[iPoint]->GetDomain()) {
				if (incompressible) {
					PrimVar_Vertex[0] = node[iPoint]->GetSolution(0);
					for (iDim = 0; iDim < nDim; iDim++)
						PrimVar_Vertex[iDim+1] = node[iPoint]->GetSolution(iDim+1);
					PrimVar_Vertex[nDim+1] = 0.0;
				}
				else {
					PrimVar_Vertex[0] = node[iPoint]->GetTemperature();
					for (iDim = 0; iDim < nDim; iDim++)
						PrimVar_Vertex[iDim+1] = node[iPoint]->GetVelocity(iDim);
					PrimVar_Vertex[nDim+1] = node[iPoint]->GetPressure();
					PrimVar_Vertex[nDim+2] = node[iPoint]->GetDensity();
				}

				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				for (iVar = 0; iVar < nVar+1; iVar++)
					for (iDim = 0; iDim < nDim; iDim++) {
						Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
						node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
					}
			}
		}
	}
	
	/*--- Update gradient value ---*/ 
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			for (iVar = 0; iVar < nVar+1; iVar++) {
				for (iDim = 0; iDim < nDim; iDim++) {
					Partial_Gradient = node[iPoint]->GetGradient_Primitive(iVar,iDim) / (geometry->node[iPoint]->GetVolume()+EPS);
//					if (Partial_Gradient > config->GetPrimGrad_Threshold()) Partial_Gradient = config->GetPrimGrad_Threshold();
//					if (Partial_Gradient < -config->GetPrimGrad_Threshold()) Partial_Gradient = -config->GetPrimGrad_Threshold();
					node[iPoint]->SetGradient_Primitive(iVar, iDim, Partial_Gradient);
				}
			}		
		}

	delete [] PrimVar_Vertex;
	delete [] PrimVar_i;
	delete [] PrimVar_j;

}

void CEulerSolution::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, jDim, iNeigh;
	unsigned long iPoint, jPoint;
	double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a, 
	r23_b, r33, weight, product;
	
	PrimVar_i = new double [nVar+1];
	PrimVar_j = new double [nVar+1]; 
	
	/*--- Loop over points of the grid ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			Coord_i = geometry->node[iPoint]->GetCoord();
			PrimVar_i[0] = node[iPoint]->GetTemperature();
			for (iDim = 0; iDim < nDim; iDim++)
				PrimVar_i[iDim+1] = node[iPoint]->GetVelocity(iDim);
			PrimVar_i[nDim+1] = node[iPoint]->GetPressure();
			PrimVar_i[nDim+2] = node[iPoint]->GetDensity();
			
			/*--- Inizialization of variables ---*/
			for (iVar = 0; iVar < nVar+1; iVar++)
				for (iDim = 0; iDim < nDim; iDim++)
					cvector[iVar][iDim] = 0.0;
			r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0; r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;
			
			for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
				jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
				Coord_j = geometry->node[jPoint]->GetCoord();
				
				PrimVar_j[0] = node[jPoint]->GetTemperature();
				for (iDim = 0; iDim < nDim; iDim++)
					PrimVar_j[iDim+1] = node[jPoint]->GetVelocity(iDim);
				PrimVar_j[nDim+1] = node[jPoint]->GetPressure();
				PrimVar_j[nDim+2] = node[jPoint]->GetDensity();
				
				weight = 0.0;
//				if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
					for (iDim = 0; iDim < nDim; iDim++)
						weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
//				}
//				else weight = 1.0;
				
				/*--- Sumations for entries of upper triangular matrix R ---*/
				r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/(weight);
				r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/(weight);
				r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/(weight);
				if (nDim == 3) {
					r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
					r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/(weight);
					r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
					r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/(weight);
				}
				
				/*--- Entries of c:= transpose(A)*b ---*/
				for (iVar = 0; iVar < nVar+1; iVar++)
					for (iDim = 0; iDim < nDim; iDim++)
						cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/(weight);
			}
			
			/*--- Entries of upper triangular matrix R ---*/
			r11 = sqrt(r11);
			r12 = r12/(r11);
			r22 = sqrt(r22-r12*r12);
			if (nDim == 3) {
				r13 = r13/(r11);
				r23 = r23_a/(r22) - r23_b*r12/(r11*r22);
				r33 = sqrt(r33-r23*r23-r13*r13);
			}
			/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
			if (nDim == 2) {
				double detR2 = (r11*r22)*(r11*r22);
				Smatrix[0][0] = (r12*r12+r22*r22)/(detR2);
				Smatrix[0][1] = -r11*r12/(detR2);
				Smatrix[1][0] = Smatrix[0][1];
				Smatrix[1][1] = r11*r11/(detR2);
			}
			else {
				double detR2 = (r11*r22*r33)*(r11*r22*r33);
				double z11, z12, z13, z22, z23, z33;
				z11 = r22*r33;
				z12 = -r12*r33;
				z13 = r12*r23-r13*r22;
				z22 = r11*r33;
				z23 = -r11*r23;
				z33 = r11*r22;
				Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/(detR2);
				Smatrix[0][1] = (z12*z22+z13*z23)/(detR2);
				Smatrix[0][2] = (z13*z33)/(detR2);
				Smatrix[1][0] = Smatrix[0][1];
				Smatrix[1][1] = (z22*z22+z23*z23)/(detR2);
				Smatrix[1][2] = (z23*z33)/(detR2);
				Smatrix[2][0] = Smatrix[0][2];
				Smatrix[2][1] = Smatrix[1][2];
				Smatrix[2][2] = (z33*z33)/(detR2);
			}
			/*--- Computation of the gradient: S*c ---*/
			for (iVar = 0; iVar < nVar+1; iVar++) {
				for (iDim = 0; iDim < nDim; iDim++) {
					product = 0.0;
					for (jDim = 0; jDim < nDim; jDim++)
						product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
//					if (product > config->GetPrimGrad_Threshold()) product = config->GetPrimGrad_Threshold();
//					if (product < -config->GetPrimGrad_Threshold()) product = -config->GetPrimGrad_Threshold();
					node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
				}
			}
		}
	
	delete [] PrimVar_i;
	delete [] PrimVar_j;
}

void CEulerSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iDim, iVar, jVar;
	double Pressure, *Normal, *Solution, ProjGridVel = 0.0, *GridVel, ProjRotVel = 0.0, *RotVel, Froude;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool GridMovement = config->GetGrid_Movement();
	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();
	bool gravity = config->GetGravityForce();


	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Compute the projected residual ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      
      double Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
			
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;
			
			if (incompressible) {
				Pressure = node[iPoint]->GetSolution(0); ///config->GetDensity_FreeStreamND();
				if (gravity) {
					Froude = config->GetVelocity_Ref() / sqrt( config->GetLength_Ref() * STANDART_GRAVITY);
					Pressure += geometry->node[iPoint]->GetCoord(1)/(Froude*Froude);
				}
			}
			else Pressure = node[iPoint]->GetPressure();
			
			/*--- Set the residual using the pressure ---*/
			Residual[0] = 0;
			for (iDim = 0; iDim < nDim; iDim++)
				Residual[iDim+1] = Pressure*UnitaryNormal[iDim]*Area;
			if (!incompressible) Residual[nVar-1] = 0;

			/*--- If there is any grid movement ---*/
			if (GridMovement) {
				Solution = node[iPoint]->GetSolution();
				GridVel = geometry->node[iPoint]->GetGridVel();
				ProjGridVel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					ProjGridVel -= GridVel[iDim]*Normal[iDim];
				for (iVar = 0; iVar < nVar; iVar++)
					Residual[iVar] -= ProjGridVel*Solution[iVar];
			}
			
			/*--- If there is a rotational frame ---*/
			if (rotating_frame) {
				RotVel = geometry->node[iPoint]->GetRotVel();
				ProjRotVel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					ProjRotVel += RotVel[iDim]*UnitaryNormal[iDim]*Area;
				Residual[nVar-1] = Pressure*ProjRotVel;
			}
			
			/*--- Add value to the residual ---*/
			node[iPoint]->AddRes_Conv(Residual);
						
			/*--- In case we are doing a implicit computation ---*/
			if (implicit) {
				if (incompressible)  {
					for (iVar = 0; iVar < nVar; iVar++) {
						for (jVar = 0; jVar < nVar; jVar++)
							Jacobian_i[iVar][jVar] = 0.0;
					}
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_i[iDim+1][0] = -Normal[iDim]; ///config->GetDensity_FreeStreamND();
					Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
				}
				else {
					double a2 = Gamma-1.0;
					double phi = 0.5*a2*node[iPoint]->GetVelocity2();

					for (iVar = 0; iVar < nVar; iVar++) {
						Jacobian_i[0][iVar] = 0.0;
						Jacobian_i[nDim+1][iVar] = 0.0;
					}
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_i[iDim+1][0] = -phi*Normal[iDim];
						for (unsigned short jDim = 0; jDim < nDim; jDim++)
							Jacobian_i[iDim+1][jDim+1] = a2*node[iPoint]->GetVelocity(jDim)*Normal[iDim];
						Jacobian_i[iDim+1][nDim+1] = -a2*Normal[iDim];
					}
					if (config->GetGrid_Movement()) {
						for(iVar = 0; iVar < nVar; iVar++)
							Jacobian_i[iVar][iVar] -= ProjGridVel;
					}
					
					Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
				}
			}
		}
	}
}

void CEulerSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iVar, jVar, iDim;
	double froude, gravityforce, Area, sq_vel, vn, rho, rhoE, c, betainc2, densityinc, **P_Matrix, **invP_Matrix, pressure, enthalpy, *Normal, 
	*UnitaryNormal, *velocity, *U_domain, *U_infty, *U_update, *W_domain, *W_infty, *W_update;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();
	
	UnitaryNormal = new double[nDim]; velocity = new double[nDim];
	U_domain = new double[nVar]; U_infty = new double[nVar]; U_update = new double[nVar];
	W_domain = new double[nVar]; W_infty = new double[nVar]; W_update = new double[nVar];
	P_Matrix = new double* [nVar]; invP_Matrix = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Matrix[iVar] = new double [nVar];
		invP_Matrix[iVar] = new double [nVar];
	}
		
	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Solution at the infinity ---*/
			if (incompressible) {
				/*--- One phase fluid ---*/
				U_infty[0] = node[iPoint]->GetDensityInc()*GetPressure_Inf();
				U_infty[1] = GetVelocity_Inf(0);
				U_infty[2] = GetVelocity_Inf(1);
				if (nDim == 3) U_infty[3] = GetVelocity_Inf(2);
			}
			else {
				U_infty[0] = GetDensity_Inf();
				U_infty[1] = GetDensity_Velocity_Inf(0);
				U_infty[2] = GetDensity_Velocity_Inf(1);
				U_infty[3] = GetDensity_Energy_Inf();
				if (nDim == 3) {
					U_infty[3] = GetDensity_Velocity_Inf(2);
					U_infty[4] = GetDensity_Energy_Inf();
				}
			}
					
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
			
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;
			
			/*--- Computation of P and inverse P matrix using values at the infinity ---*/
			sq_vel = 0.0; vn = 0.0;
			rho = U_infty[0]; rhoE = U_infty[nVar-1];
			betainc2 = node[iPoint]->GetBetaInc2();
			densityinc = node[iPoint]->GetDensityInc();

			for (iDim = 0; iDim < nDim; iDim++) {
				if (incompressible) velocity[iDim] = U_infty[iDim+1];
				else velocity[iDim] = U_infty[iDim+1]/rho;
				sq_vel += velocity[iDim]*velocity[iDim];
				vn += velocity[iDim]*UnitaryNormal[iDim]*Area;
			}
			
			/*--- For a rotational frame ---*/
			if (rotating_frame) {
				double ProjRotVel = 0.0;
				double *RotVel = geometry->node[iPoint]->GetRotVel();
				for (iDim = 0; iDim < nDim; iDim++)
					ProjRotVel += RotVel[iDim]*UnitaryNormal[iDim]*Area;
				vn -= ProjRotVel;
			}
			
			if (incompressible) c = sqrt(vn*vn + betainc2 * Area * Area)/Area;
			else c = sqrt(Gamma*Gamma_Minus_One*(rhoE/rho-0.5*sq_vel));
			
			if (incompressible) {
				solver->GetPArtCompMatrix_inv(&densityinc, velocity, &betainc2, UnitaryNormal, invP_Matrix);
				solver->GetPArtCompMatrix(&densityinc, velocity, &betainc2, UnitaryNormal, P_Matrix);
				
				
				/*--- Check factorization ---*/
/*				double **A_Matrix, **A_Matrix_Projected, *Lambda;
				A_Matrix = new double* [nVar]; A_Matrix_Projected = new double* [nVar];
				Lambda = new double [nVar];
				for (iVar = 0; iVar < nVar; iVar++) {
					A_Matrix[iVar] = new double [nVar];
					A_Matrix_Projected[iVar] = new double [nVar];
				}
				solver->GetInviscidArtCompProjJac(densityinc, velocity, betainc2, Normal, 1.0, A_Matrix);

				Lambda[0] = vn;
				Lambda[1] = vn + c*Area;
				Lambda[2] = vn - c*Area;
				
				double Proj_ModJac_Tensor_ij;
				for (int iVar = 0; iVar < nVar; iVar++) {
					for (int jVar = 0; jVar < nVar; jVar++) { 
						Proj_ModJac_Tensor_ij = 0.0;
						for (int kVar = 0; kVar < nVar; kVar++) 
							Proj_ModJac_Tensor_ij += P_Matrix[iVar][kVar]*Lambda[kVar]*invP_Matrix[kVar][jVar];
						A_Matrix_Projected[iVar][jVar] += Proj_ModJac_Tensor_ij;
					}
				}
							
				cout << A_Matrix[0][0] <<" "<< A_Matrix[0][1] <<" "<< A_Matrix[0][2] << endl;
				cout << A_Matrix[1][0] <<" "<< A_Matrix[1][1] <<" "<< A_Matrix[1][2] << endl;
				cout << A_Matrix[2][0] <<" "<< A_Matrix[2][1] <<" "<< A_Matrix[2][2] << endl<< endl;
				
				cout << A_Matrix_Projected[0][0] <<" "<< A_Matrix_Projected[0][1] <<" "<< A_Matrix_Projected[0][2] << endl;
				cout << A_Matrix_Projected[1][0] <<" "<< A_Matrix_Projected[1][1] <<" "<< A_Matrix_Projected[1][2] << endl;
				cout << A_Matrix_Projected[2][0] <<" "<< A_Matrix_Projected[2][1] <<" "<< A_Matrix_Projected[2][2] << endl;

				cin.get();*/

				
			}
			else {
				solver->GetPMatrix_inv(&rho, velocity, &c, UnitaryNormal, invP_Matrix);
				solver->GetPMatrix(&rho, velocity, &c, UnitaryNormal, P_Matrix);
			}
			
			/*--- computation of characteristics variables at the wall ---*/			
			for (iVar=0; iVar < nVar; iVar++) {
				W_domain[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_domain[iVar] +=invP_Matrix[iVar][jVar]*U_domain[jVar];
			}
			
			/*--- computation of characteristics variables at the infinity ---*/			
			for (iVar=0; iVar < nVar; iVar++) {
				W_infty[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_infty[iVar] +=invP_Matrix[iVar][jVar]*U_infty[jVar];
			}
						
			/*--- fix characteristics value ---*/
			if (nDim == 2) {
				if (incompressible) {
					if(vn > 0.0) W_update[0] = W_domain[0];
					else W_update[0] = W_infty[0];
					if(vn+c*Area > 0.0) W_update[1] = W_domain[1];
					else W_update[1] = W_infty[1];
					if(vn-c*Area > 0.0) W_update[2] = W_domain[2];
					else W_update[2] = W_infty[2];	
				}
				else {
					if(vn > 0.0) { 
						W_update[0] = W_domain[0];
						W_update[1] = W_domain[1];
					}
					else {
						W_update[0] = W_infty[0];
						W_update[1] = W_infty[1];
					}
					
					if(vn+c*Area > 0.0) W_update[2] = W_domain[2];
					else W_update[2] = W_infty[2];
					
					if(vn-c*Area > 0.0) W_update[3] = W_domain[3];
					else W_update[3] = W_infty[3];
				}
			}
			
			if (nDim == 3) {
				if(vn > 0.0) { 
					W_update[0] = W_domain[0];
					W_update[1] = W_domain[1];
					W_update[2] = W_domain[2];
				}
				else {
					W_update[0] = W_infty[0];
					W_update[1] = W_infty[1];
					W_update[2] = W_infty[2];
				}
				
				if(vn+c*Area > 0.0) W_update[3] = W_domain[3];
				else W_update[3] = W_infty[3];
				
				if(vn-c*Area > 0.0) W_update[4] = W_domain[4];
				else W_update[4] = W_infty[4];
			}
			
			/*--- conservative variables using characteristics ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				U_update[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					U_update[iVar] +=P_Matrix[iVar][jVar]*W_update[jVar];
			}
			
			if (incompressible) {
				/*--- Residual computation incompressible flow ---*/
				rho = node[iPoint]->GetDensityInc();
				betainc2 = node[iPoint]->GetBetaInc2();
				froude = config->GetVelocity_Ref() / sqrt( config->GetLength_Ref() * STANDART_GRAVITY);
				gravityforce = geometry->node[iPoint]->GetCoord(nDim-1) / (froude*froude);
				pressure = U_update[0];
				for (iDim = 0; iDim < nDim; iDim++)
					velocity[iDim] = U_update[iDim+1];
				
				solver->GetInviscidArtCompProjFlux(&rho, velocity, &pressure, &gravityforce, &betainc2, UnitaryNormal, Residual);
			}
			else {
				/*--- Residual computation compressible flow ---*/
				rho = U_update[0]; rhoE = U_update[nVar-1];
				sq_vel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					velocity[iDim] = U_update[iDim+1]/rho;
					sq_vel +=velocity[iDim]*velocity[iDim];
				}
				c = sqrt(Gamma*Gamma_Minus_One*(rhoE/rho - 0.5*sq_vel));
				pressure = (c * c * rho) / Gamma;
				enthalpy = (rhoE + pressure) / rho;
				
				solver->GetInviscidProjFlux(&rho, velocity, &pressure, &enthalpy, UnitaryNormal, Residual);
			}
			
			if (rotating_frame) {
				double ProjRotVel = 0.0;
				double *RotVel = geometry->node[iPoint]->GetRotVel();
				for (iDim = 0; iDim < nDim; iDim++)
					ProjRotVel += RotVel[iDim]*UnitaryNormal[iDim];
				for (iVar = 0; iVar < nVar; iVar++)
					Residual[iVar] -= ProjRotVel * U_update[iVar];
			}
			
			for	(iVar = 0; iVar < nVar; iVar++) 
				Residual[iVar] = Residual[iVar]*Area;
			node[iPoint]->AddRes_Conv(Residual);
			
			/*--- In case we are doing a implicit computation ---*/
			if (implicit) {
				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
				solver->SetNormal(Vector);
				solver->SetConservative(U_domain, U_infty);
				if (incompressible) {
					solver->SetDensityInc(node[iPoint]->GetDensityInc(), node[iPoint]->GetDensityInc());
					solver->SetBetaInc2(node[iPoint]->GetBetaInc2(), node[iPoint]->GetBetaInc2());
					double *iCoord = NULL;
					iCoord = geometry->node[iPoint]->GetCoord();
					solver->SetCoord(iCoord, iCoord);
				}
				solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			}
			
		}
	}
	
	delete [] UnitaryNormal; delete [] velocity;
	delete [] U_domain; delete [] U_infty; delete [] U_update;
	delete [] W_domain; delete [] W_infty; delete [] W_update;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Matrix[iVar];
		delete [] invP_Matrix[iVar];
	}
	delete [] P_Matrix;
	delete [] invP_Matrix;
}

void CEulerSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	double *U_domain, *U_inlet, Pressure_Total, Temperature_Total, Velocity[3], Velocity2, Enthalpy_Total, 
	Enthalpy_Static, Temperature_Static, Pressure_Static, Density, Energy, *Flow_Direction, Mach_Inlet, *Velocity_Inlet, Velocity_Inlet2;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
	double Gas_Constant = config->GetGas_Constant();
	
	U_domain = new double[nVar]; U_inlet = new double[nVar];
	
	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);
			
			if (incompressible) {
				
				U_inlet[0] = node[iPoint]->GetSolution(0);
				for (iDim = 0; iDim < nDim; iDim++)
					U_inlet[iDim+1] = GetVelocity_Inf(iDim);
				
			}
			else {
        
        /*--- Retrieve the specified flow quantities for this inlet boundary. ---*/
				Pressure_Total    = config->GetInlet_Ptotal(config->GetMarker_All_Tag(val_marker));
				Temperature_Total = config->GetInlet_Ttotal(config->GetMarker_All_Tag(val_marker));
        Flow_Direction    = config->GetInlet_FlowDir(config->GetMarker_All_Tag(val_marker));
        Velocity_Inlet    = config->GetInlet_FlowDir(config->GetMarker_All_Tag(val_marker));

        
				/*--- Interpolate the velocity from the interior of the domain.  ---*/
				for (iDim = 0; iDim < nDim; iDim++)
					Velocity[iDim] = node[iPoint]->GetVelocity(iDim);
				Velocity2 = node[iPoint]->GetVelocity2();
				
        /*--- Check whether the flow is supersonic at the inlet. ---*/
        Mach_Inlet = sqrt(Velocity2)/node[iPoint]->GetSoundSpeed();
        
//        if (Mach_Inlet >= 1.0) {
//          cout << "Supersonic Inflow: Mach number = " << Mach_Inlet << endl;
//          /*--- Supseronic inflow: there are no characteristics moving upstream,
//           so all flow quantities can be specified at the inlet. We use the
//           specified total quantites to find the total density, then the 
//           density from the isentropic relations. ---*/
//          
//          /*--- Use the specified velocity components.  ---*/
//          Velocity2 = 0.0;
//          for (iDim = 0; iDim < nDim; iDim++) {
//            Velocity[iDim] = Velocity_Inlet[iDim];
//            Velocity2 += Velocity_Inlet[iDim]*Velocity_Inlet[iDim];
//          }
//          
//          /*--- With the total temperature we compute the total enthalpy ---*/
//          Enthalpy_Total = (Gamma * Gas_Constant / Gamma_Minus_One) * Temperature_Total;
//          
//          /*--- If we subtract the $0.5*velocity^2$ we obtain the static enthalpy (at the inlet) ---*/
//          Enthalpy_Static = Enthalpy_Total - 0.5*Velocity2;
//          
//          /*--- With the static enthalpy (at the inlet) it is possible to compute the static temperature (at the inlet) ---*/
//          Temperature_Static = Enthalpy_Static * Gamma_Minus_One / (Gamma * Gas_Constant);
//          
//          /*--- With the static temperature (at the inlet), and the total temperature (in the stagnation point), using the 
//           isentropic relations it is possible to compute the static pressure at the inlet ---*/
//          Pressure_Static = Pressure_Total * pow((Temperature_Static/Temperature_Total), Gamma/Gamma_Minus_One);
//          
//          /*--- Using the static pressure (at the inlet) and static temperature (at the inlet) we will compute 
//           the density (at the inlet) ---*/
//          Density = Pressure_Static / (Gas_Constant * Temperature_Static);
//          
//        } else {
          
          /*--- Subsonic inflow: one characteristic is moving upstream,
           therefore we cannot specify the velocity, and have to use the
           other flow variables from the interior for the update. ---*/
          
          /*--- Find the unit vector of the specified velocity.  ---*/
          Velocity_Inlet2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim]   = Velocity_Inlet[iDim];
            Velocity_Inlet2 += Velocity_Inlet[iDim]*Velocity_Inlet[iDim];
          }
          for (iDim = 0; iDim < nDim; iDim++) {
            Flow_Direction[iDim] = Velocity_Inlet[iDim]/sqrt(Velocity_Inlet2);
          }
          
          /*--- Correct the interior flow velocity to be in the specified direction at the inlet. ---*/
          Velocity[0] = Flow_Direction[0]*sqrt(Velocity2); 
          Velocity[1] = Flow_Direction[1]*sqrt(Velocity2); 
          if (nDim == 3) Velocity[2] = Flow_Direction[2]*sqrt(Velocity2); 
          
          /*--- With the total temperature we compute the total enthalpy ---*/
          Enthalpy_Total = (Gamma * Gas_Constant / Gamma_Minus_One) * Temperature_Total;
          
          /*--- If we subtract the $0.5*velocity^2$ we obtain the static enthalpy (at the inlet) ---*/
          Enthalpy_Static = Enthalpy_Total - 0.5*Velocity2;
          
          /*--- With the static enthalpy (at the inlet) it is possible to compute the static temperature (at the inlet) ---*/
          Temperature_Static = Enthalpy_Static * Gamma_Minus_One / (Gamma * Gas_Constant);
          
          /*--- With the static temperature (at the inlet), and the total temperature (in the stagnation point), using the 
           isentropic relations it is possible to compute the static pressure at the inlet ---*/
          Pressure_Static = Pressure_Total * pow((Temperature_Static/Temperature_Total), Gamma/Gamma_Minus_One);
          
          /*--- Using the static pressure (at the inlet) and static temperature (at the inlet) we will compute 
           the density (at the inlet) ---*/
          Density = Pressure_Static / (Gas_Constant * Temperature_Static);
          
//        }
        
				/*--- Using pressure, density, and velocity we can compute the energy at the inlet ---*/
				Energy = Pressure_Static/(Density*Gamma_Minus_One)+0.5*Velocity2;
				
				/*--- Conservative variables, using the derived quatities ---*/
				U_inlet[0] = Density;
				U_inlet[1] = Velocity[0]*Density;
				U_inlet[2] = Velocity[1]*Density;
				U_inlet[3] = Energy*Density;
				if (nDim == 3) {
					U_inlet[3] = Velocity[2]*Density;
					U_inlet[4] = Energy*Density;
				}
			}
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			solver->SetNormal(Vector);
			
			/*--- Compute the residual using an upwind scheme ---*/
			solver->SetConservative(U_domain, U_inlet);
			if (incompressible) {
				solver->SetDensityInc(node[iPoint]->GetDensityInc(), node[iPoint]->GetDensityInc());
				solver->SetBetaInc2(node[iPoint]->GetBetaInc2(), node[iPoint]->GetBetaInc2());
				double *iCoord = NULL;
				iCoord = geometry->node[iPoint]->GetCoord();
				solver->SetCoord(iCoord, iCoord);
			}
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddRes_Conv(Residual);
			
			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			
		}
	}
}

void CEulerSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	double *U_domain, *U_outlet, Pressure, Velocity[3], Velocity2, Mach_Exit;
  double Density, Energy, Temperature, yFreeSurface, PressFreeSurface;
	double DensityFreeSurface, Froude, yCoord;
	
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	double Gas_Constant = config->GetGas_Constant();
	bool incompressible = config->GetIncompressible();
	bool gravity = config->GetGravityForce();
	
	U_domain = new double[nVar]; U_outlet = new double[nVar];
	
	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Interpolated solution interior to the outlet ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);
			
			if (incompressible) {
				
				/*--- The pressure is defined with the free-stream values, and the 
				 total velocity (corrected) is computed from the interior---*/
				for (iDim = 0; iDim < nDim; iDim++)
					Velocity[iDim] = node[iPoint]->GetSolution(iDim+1);
				Velocity2 = node[iPoint]->GetVelocity2();
				Velocity[0] = sqrt(Velocity2); Velocity[1] = 0.0;
				if (nDim == 3) Velocity[2] = 0.0;
				
				/*--- Bernoulli principle will be used to compute the pressure at each point, we are supposing that  ---*/				
				yFreeSurface = config->GetLevelSet_Zero();
				PressFreeSurface = GetPressure_Inf();
				DensityFreeSurface = GetDensity_Inf();
				Density = node[iPoint]->GetDensityInc();
				Froude = config->GetVelocity_Ref() / sqrt( config->GetLength_Ref() * STANDART_GRAVITY);
				yCoord = geometry->node[iPoint]->GetCoord(1);
				if (gravity) Pressure = Density*( PressFreeSurface/DensityFreeSurface + (yFreeSurface-yCoord)/(Froude*Froude));
				else Pressure = Density*( PressFreeSurface/DensityFreeSurface );
				
//				U_outlet[0] = node[iPoint]->GetSolution(0);
				U_outlet[0] = Pressure/DensityFreeSurface;
//				U_outlet[0] = PressFreeSurface/DensityFreeSurface;
//				U_outlet[0] = PressFreeSurface / Density;
				
				U_outlet[1] = Velocity[0];
				U_outlet[2] = Velocity[1];
				if (nDim == 3) U_outlet[3] = Velocity[2];
			} 
			
			else {
        
				/*--- Retrieve the specified back pressure for this outlet boundary. ---*/
				Pressure = config->GetOutlet_Pressure(config->GetMarker_All_Tag(val_marker));
				
				/*--- The velocity, and the temperature are interpolated (zero order),
				 from the domain, be careful everything is precomputed at the preprocesing stage  ---*/
				Temperature = node[iPoint]->GetTemperature();
				for (iDim = 0; iDim < nDim; iDim++)
					Velocity[iDim] = node[iPoint]->GetVelocity(iDim);
				Velocity2 = node[iPoint]->GetVelocity2();
				
        /*--- Check whether the flow is supersonic at the exit. ---*/
        Mach_Exit = sqrt(Velocity2)/node[iPoint]->GetSoundSpeed();
        
        if (Mach_Exit >= 1.0) {
          
          /*--- Supseronic exit flow: there are no characteristics entering the 
           domain, and the solution can simply be copied from the interior. ---*/
          
          /*--- Conservative variables, using the derived quantities ---*/
          U_outlet[0] = U_domain[0];
          U_outlet[1] = U_domain[1];
          U_outlet[2] = U_domain[2];
          U_outlet[3] = U_domain[3];
          if (nDim == 3) {
            U_outlet[3] = U_domain[3];
            U_outlet[4] = U_domain[4];
          }
        
        } else {
          
          /*--- Subsonic exit flow: one characteristic is entering the domain,
           therefore one variable can be specified (back pressure) and is used
           to update the conservative variables. ---*/
          
          /*--- Compute the density, and the energy at the outlet ---*/
          Density = Pressure/(Gas_Constant * Temperature);
          Energy  = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
          
          /*--- Conservative variables, using the derived quantities ---*/
          U_outlet[0] = Density;
          U_outlet[1] = Velocity[0]*Density;
          U_outlet[2] = Velocity[1]*Density;
          U_outlet[3] = Energy*Density;
          if (nDim == 3) {
            U_outlet[3] = Velocity[2]*Density;
            U_outlet[4] = Energy*Density;
          }
        }
			}
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			solver->SetNormal(Vector);
			
			/*--- Compute the residual using an upwind scheme ---*/
			solver->SetConservative(U_domain, U_outlet);
			if (incompressible) {
				solver->SetDensityInc(node[iPoint]->GetDensityInc(), node[iPoint]->GetDensityInc());
				solver->SetBetaInc2(node[iPoint]->GetBetaInc2(), node[iPoint]->GetBetaInc2());
				double *iCoord = NULL;
				iCoord = geometry->node[iPoint]->GetCoord();
				solver->SetCoord(iCoord, iCoord);
			}
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddRes_Conv(Residual);
			
			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			
		}
	}
	
}

void CEulerSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iDim, iVar;
	double Pressure, *Normal;
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			/*--- Compute the projected residual ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Pressure = node[iPoint]->GetPressure();
			Residual[0] = 0; Residual[nDim+1] = 0;
			for (iDim = 0; iDim < nDim; iDim++)
				Residual[iDim+1] = -Pressure*Normal[iDim];
					
			/*--- Add value to the residual ---*/
			node[iPoint]->AddRes_Conv(Residual);
//			cout << "Resid "<< Pressure << " " << Residual[0] << " " << Residual[1] << " " <<Residual[2] << " " <<Residual[3] << endl<< endl;

		
			/*--- In case we are doing a implicit computation ---*/
			if (implicit) {
				double a2 = Gamma-1.0;
				double phi = 0.5*a2*node[iPoint]->GetVelocity2();
				for (iVar = 0; iVar < nVar; iVar++) {
					Jacobian_i[0][iVar] = 0.0;
					Jacobian_i[nDim+1][iVar] = 0.0;
				}
				for (iDim = 0; iDim < nDim; iDim++) {
					Jacobian_i[iDim+1][0] = -phi*Normal[iDim];
					for (unsigned short jDim = 0; jDim < nDim; jDim++)
						Jacobian_i[iDim+1][jDim+1] = a2*node[iPoint]->GetVelocity(jDim)*Normal[iDim];
					Jacobian_i[iDim+1][nDim+1] = -a2*Normal[iDim];
				}
				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}
		}
	}
}

void CEulerSolution::BC_Interface_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
																					 CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, jPoint;
	unsigned short iVar, iDim;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
#ifdef NO_MPI
	
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPoint();
			
			/*--- Store the solution for both points ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar); 
				Solution_j[iVar] = node[jPoint]->GetSolution(iVar); 
			}
			
			/*--- Set Conservative Variables ---*/
			solver->SetConservative(Solution_i, Solution_j);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++)
				Vector[iDim] = -Vector[iDim]; 
			solver->SetNormal(Vector);
			
			/*--- Add Residuals and Jacobians ---*/
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddRes_Conv(Residual);
			if (implicit) Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);	
		}
	}
	
#else
	
	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	double Buffer_Send_U[nVar], *Conserv_Var, Buffer_Receive_U[nVar];
	
	/*--- Do the send process, by the moment we are sending each 
	 node individually, this must be changed ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
			
			/*--- We only send the information that belong to other boundary ---*/
			if (jProcessor != rank) {
				Conserv_Var = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Send_U[iVar] = Conserv_Var[iVar];
				MPI::COMM_WORLD.Bsend(&Buffer_Send_U, nVar, MPI::DOUBLE, jProcessor, iPoint);
			}
		}		
	}
	
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
			
			/*--- We only receive the information that belong to other boundary ---*/
			if (jProcessor != rank)
				MPI::COMM_WORLD.Recv(&Buffer_Receive_U, nVar, MPI::DOUBLE, jProcessor, jPoint);
			else {
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Receive_U[iVar] = node[jPoint]->GetSolution(iVar); 
			}

			/*--- Store the solution for both points ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar); 
				Solution_j[iVar] = Buffer_Receive_U[iVar];
				Solution[iVar] = 0.5*(Solution_i[iVar]+Solution_j[iVar]);
			}
			
			solver->SetConservative(Solution_i, Solution_j);
//			solver->SetConservative(Solution_i, Solution);
			
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++)
				Vector[iDim] = -Vector[iDim]; 
			solver->SetNormal(Vector);
			
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddRes_Conv(Residual);
			if (implicit) Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
		}

	}
	
#endif
	
}

void CEulerSolution::BC_NearField_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
																		 CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, jPoint;
	unsigned short iVar, iDim;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
#ifdef NO_MPI
	
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPoint();
			
			/*--- Store the solution for both points ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar); 
				Solution_j[iVar] = node[jPoint]->GetSolution(iVar); 
			}
			
			/*--- Set Conservative Variables ---*/
			solver->SetConservative(Solution_i, Solution_j);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++)
				Vector[iDim] = -Vector[iDim]; 
			solver->SetNormal(Vector);
			
			/*--- Add Residuals and Jacobians ---*/
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddRes_Conv(Residual);
			if (implicit) Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);	
		}
	}
	
#else
	
	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	double Buffer_Send_U[nVar], *Conserv_Var, Buffer_Receive_U[nVar];

	/*--- Do the send process, by the moment we are sending each 
	 node individually, this must be changed ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
			/*--- We only send the information that belong to other boundary, -1 processor 
			 means that the boundary condition is not applied ---*/
			if ((jProcessor != rank) && (jProcessor != -1)) {
				Conserv_Var = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Send_U[iVar] = Conserv_Var[iVar];
				MPI::COMM_WORLD.Bsend(&Buffer_Send_U, nVar, MPI::DOUBLE, jProcessor, iPoint);
			}
		}
	}
	
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
			
			/*--- -1 processor means that the boundary condition is not applied ---*/
			if (jProcessor != -1) {
				/*--- We only receive the information that belong to other boundary ---*/
				if (jProcessor != rank)
					MPI::COMM_WORLD.Recv(&Buffer_Receive_U, nVar, MPI::DOUBLE, jProcessor, jPoint);
				else {
					for (iVar = 0; iVar < nVar; iVar++)
						Buffer_Receive_U[iVar] = node[jPoint]->GetSolution(iVar); 
				}
				
				/*--- Store the solution for both points ---*/
				for (iVar = 0; iVar < nVar; iVar++) {
					Solution_i[iVar] = node[iPoint]->GetSolution(iVar); 
					Solution_j[iVar] = Buffer_Receive_U[iVar]; 
				}
				
				/*--- Set Conservative Variables ---*/
				solver->SetConservative(Solution_i, Solution_j);
				
				/*--- Set Normal ---*/
				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = -Vector[iDim]; 
				solver->SetNormal(Vector);
				
				/*--- Add Residuals and Jacobians ---*/
				solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
				node[iPoint]->AddRes_Conv(Residual);
				if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);	
			}
		}
	}
	
#endif
	
}

void CEulerSolution::BC_Dirichlet(CGeometry *geometry, CSolution **solution_container, 
																					 CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, jPoint, total_index;
	unsigned short iVar;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
#ifdef NO_MPI
	
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPoint();
			
			/*--- Store the solution for both points ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar); 
				Solution_j[iVar] = node[jPoint]->GetSolution(iVar); 
			}
			
			/*--- Set the residual, truncation error and solution value ---*/
			node[iPoint]->SetSolution_Old(Solution_j);
			node[iPoint]->SetSolution(Solution_j);
			node[iPoint]->Set_ResConv_Zero();
			node[iPoint]->Set_ResVisc_Zero();
			node[iPoint]->SetTruncationErrorZero();
			
			/*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
			if (implicit)
				for (iVar = 0; iVar < nVar; iVar++) {
					total_index = iPoint*nVar+iVar;
					Jacobian.DeleteValsRowi(total_index);
				}
			
		}
	}
	
#else
	
	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	double Buffer_Send_U[nVar], *Conserv_Var, Buffer_Receive_U[nVar];
	
	/*--- Do the send process, by the moment we are sending each 
	 node individually, this must be changed ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
			
			/*--- We only send the information that belong to other boundary ---*/
			if (jProcessor != rank) {
				Conserv_Var = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Send_U[iVar] = Conserv_Var[iVar];
				MPI::COMM_WORLD.Bsend(&Buffer_Send_U, nVar, MPI::DOUBLE, jProcessor, iPoint);
			}
		}		
	}
	
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
			
			/*--- We only receive the information that belong to other boundary ---*/
			if (jProcessor != rank)
				MPI::COMM_WORLD.Recv(&Buffer_Receive_U, nVar, MPI::DOUBLE, jProcessor, jPoint);
			else {
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Receive_U[iVar] = node[jPoint]->GetSolution(iVar); 
			}
			
			/*--- Store the solution for both points ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar); 
				Solution_j[iVar] = Buffer_Receive_U[iVar]; 
			}
			
			/*--- Set the residual, truncation error and solution value ---*/
			node[iPoint]->SetSolution_Old(Solution_j);
			node[iPoint]->SetSolution(Solution_j);
			node[iPoint]->Set_ResConv_Zero();
			node[iPoint]->Set_ResVisc_Zero();
			node[iPoint]->SetTruncationErrorZero();
			
			/*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
			if (implicit)
				for (iVar = 0; iVar < nVar; iVar++) {
					total_index = iPoint*nVar+iVar;
					Jacobian.DeleteValsRowi(total_index);
				}
			
		}
	}
	
#endif
	
}


void CEulerSolution::BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
									 unsigned short val_marker, unsigned short val_mesh) {
	unsigned short iVar;
	unsigned long iVertex, iPoint;
	double *Conserv_Var, *Conserv_Undivided_Laplacian = NULL;
	
	short SendRecv = config->GetMarker_All_SendRecv(val_marker);
	bool incompressible = config->GetIncompressible();
	
	/*--- Send information  ---*/
	if (SendRecv > 0) {
		
		/*--- Sending only occurs with MPI ---*/
#ifndef NO_MPI
		
		double **Conserv_Grad = NULL, *Grad_Limit = NULL;
		
		/*--- Upwind scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {
			double Buffer_Send_U[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Ux[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Uy[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Uz[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Limit[geometry->nVertex[val_marker]][nVar];
			
			unsigned long nBuffer = geometry->nVertex[val_marker]*nVar;
			int send_to = SendRecv;
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
				Conserv_Var = node[iPoint]->GetSolution();
				if (val_mesh == MESH_0) Conserv_Grad = node[iPoint]->GetGradient();
				if ((val_mesh == MESH_0) && (config->GetKind_SlopeLimit_Flow() != NONE)) 
					Grad_Limit = node[iPoint]->GetLimiter();
				
				for (iVar = 0; iVar < nVar; iVar++) {
					Buffer_Send_U[iVertex][iVar] = Conserv_Var[iVar];
					if (val_mesh == MESH_0) {
						Buffer_Send_Ux[iVertex][iVar] = Conserv_Grad[iVar][0];					
						Buffer_Send_Uy[iVertex][iVar] = Conserv_Grad[iVar][1];
						if (nDim == 3) Buffer_Send_Uz[iVertex][iVar] = Conserv_Grad[iVar][2];
						if (config->GetKind_SlopeLimit_Flow() != NONE) 
							Buffer_Send_Limit[iVertex][iVar] = Grad_Limit[iVar];
					}
				}
				
			}
			MPI::COMM_WORLD.Bsend(&Buffer_Send_U,nBuffer,MPI::DOUBLE,send_to, 0);
			if (val_mesh == MESH_0) {
				MPI::COMM_WORLD.Bsend(&Buffer_Send_Ux,nBuffer,MPI::DOUBLE,send_to, 1);
				MPI::COMM_WORLD.Bsend(&Buffer_Send_Uy,nBuffer,MPI::DOUBLE,send_to, 2);
				if (nDim == 3) MPI::COMM_WORLD.Bsend(&Buffer_Send_Uz,nBuffer,MPI::DOUBLE,send_to, 3);
				if (config->GetKind_SlopeLimit_Flow() != NONE) MPI::COMM_WORLD.Bsend(&Buffer_Send_Limit,nBuffer,MPI::DOUBLE,send_to, 4);
			}
		}
		
		/*--- Centered scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {
			
			double Buffer_Send_U[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Undivided_Laplacian[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Sensor[geometry->nVertex[val_marker]];
			double Buffer_Send_Lambda[geometry->nVertex[val_marker]];
			double Buffer_Send_BetaInc2[geometry->nVertex[val_marker]];
			double Buffer_Send_DensityInc[geometry->nVertex[val_marker]];
			unsigned short Buffer_Send_Neighbor[geometry->nVertex[val_marker]];

			unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
			unsigned long nBuffer_Scalar = geometry->nVertex[val_marker];
			int send_to = SendRecv;
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
				Conserv_Var = node[iPoint]->GetSolution();
				if (val_mesh == MESH_0) Conserv_Undivided_Laplacian = node[iPoint]->GetUnd_Lapl();
				for (iVar = 0; iVar < nVar; iVar++) {
					Buffer_Send_U[iVertex][iVar] = Conserv_Var[iVar];
					if (val_mesh == MESH_0) Buffer_Send_Undivided_Laplacian[iVertex][iVar] = Conserv_Undivided_Laplacian[iVar];					
				}
				if (val_mesh == MESH_0) Buffer_Send_Sensor[iVertex] = node[iPoint]->GetSensor();
				Buffer_Send_Lambda[iVertex] = node[iPoint]->GetLambda();
				Buffer_Send_Neighbor[iVertex] = geometry->node[iPoint]->GetnPoint();
				
				if (incompressible) {
					Buffer_Send_BetaInc2[iVertex] = node[iPoint]->GetBetaInc2();
					Buffer_Send_DensityInc[iVertex] = node[iPoint]->GetDensityInc();
				}
				
			}
			
			MPI::COMM_WORLD.Bsend(&Buffer_Send_U,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Bsend(&Buffer_Send_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Bsend(&Buffer_Send_Sensor,nBuffer_Scalar,MPI::DOUBLE,send_to, 2);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_Lambda,nBuffer_Scalar,MPI::DOUBLE,send_to, 3);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_Neighbor,nBuffer_Scalar,MPI::UNSIGNED_SHORT,send_to, 4);

			if (incompressible) {
				MPI::COMM_WORLD.Bsend(&Buffer_Send_BetaInc2,nBuffer_Scalar,MPI::DOUBLE,send_to, 5);
				MPI::COMM_WORLD.Bsend(&Buffer_Send_DensityInc,nBuffer_Scalar,MPI::DOUBLE,send_to, 6);
			}
		}
#endif
	}
	
	/*--- Receive information  ---*/
	if (SendRecv < 0) {
		/*--- Upwind scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {
			double Buffer_Receive_U[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Ux[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Uy[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Uz[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Limit[geometry->nVertex[val_marker]][nVar];
			
#ifdef NO_MPI
			cout << "Upwinding for periodic boundaries in serial not yet supported." << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(0);
#else
			unsigned long nBuffer = geometry->nVertex[val_marker]*nVar;
			int receive_from = abs(SendRecv);
			
			MPI::COMM_WORLD.Recv(&Buffer_Receive_U,nBuffer,MPI::DOUBLE,receive_from, 0);
			if (val_mesh == MESH_0) {
				MPI::COMM_WORLD.Recv(&Buffer_Receive_Ux,nBuffer,MPI::DOUBLE,receive_from, 1);
				MPI::COMM_WORLD.Recv(&Buffer_Receive_Uy,nBuffer,MPI::DOUBLE,receive_from, 2);
				if (nDim == 3) MPI::COMM_WORLD.Recv(&Buffer_Receive_Uz,nBuffer,MPI::DOUBLE,receive_from, 3);
				if (config->GetKind_SlopeLimit_Flow() != NONE) MPI::COMM_WORLD.Recv(&Buffer_Receive_Limit,nBuffer,MPI::DOUBLE,receive_from, 4);
			}
#endif
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				iPoint = geometry->vertex[val_marker][iVertex]->GetNode();	
				for (iVar = 0; iVar < nVar; iVar++) {
					node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVertex][iVar]);
					if (val_mesh == MESH_0) {
						node[iPoint]->SetGradient(iVar, 0, Buffer_Receive_Ux[iVertex][iVar]);
						node[iPoint]->SetGradient(iVar, 1, Buffer_Receive_Uy[iVertex][iVar]);
						if (nDim == 3) node[iPoint]->SetGradient(iVar, 2, Buffer_Receive_Uz[iVertex][iVar]);
						if (config->GetKind_SlopeLimit_Flow() != NONE) node[iPoint]->SetLimiter(iVar, Buffer_Receive_Limit[iVertex][iVar]);
					}
				}
			}
		}
		
		/*--- Centered scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {
			double Buffer_Receive_U[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Undivided_Laplacian[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Sensor[geometry->nVertex[val_marker]];
			double Buffer_Receive_Lambda[geometry->nVertex[val_marker]];
			unsigned short Buffer_Receive_Neighbor[geometry->nVertex[val_marker]];
			double Buffer_Receive_BetaInc2[geometry->nVertex[val_marker]];
			double Buffer_Receive_DensityInc[geometry->nVertex[val_marker]];
			
			int receive_from = abs(SendRecv);
			double rotMatrix[3][3], *angles, newSolution[nVar];
			double theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
			unsigned short iPeriodic_Index;
			
			/*--- Allow for periodic boundaries to use SEND_RECEIVE in serial.
			 Serial computations will only use the BC in receive mode, as
			 the proc will always be sending information to itself. ---*/
#ifdef NO_MPI
			/*--- Retrieve the donor information from the matching marker ---*/

			unsigned short donor_marker = 1;
			unsigned long donorPoint;
			for (unsigned short iMark = 0; iMark < config->GetnMarker_All(); iMark++){
				if (config->GetMarker_All_SendRecv(iMark) == receive_from) donor_marker = iMark;

			}
			
			/*--- Get the information from the donor directly. This is a serial
			 computation with access to all nodes. Note that there is an
			 implicit ordering in the list. ---*/
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				donorPoint = geometry->vertex[donor_marker][iVertex]->GetNode();
				iPeriodic_Index = geometry->vertex[val_marker][iVertex]->GetRotation_Type();

				Conserv_Var = node[donorPoint]->GetSolution();
				if (val_mesh == MESH_0) Conserv_Undivided_Laplacian = node[donorPoint]->GetUnd_Lapl();
				for (iVar = 0; iVar < nVar; iVar++) {
					Buffer_Receive_U[iVertex][iVar] = Conserv_Var[iVar];
					if (val_mesh == MESH_0) Buffer_Receive_Undivided_Laplacian[iVertex][iVar] = Conserv_Undivided_Laplacian[iVar];					
				}
				if (val_mesh == MESH_0) Buffer_Receive_Sensor[iVertex] = node[donorPoint]->GetSensor();
				Buffer_Receive_Lambda[iVertex] = node[donorPoint]->GetLambda();

				if (incompressible) {
					Buffer_Receive_BetaInc2[iVertex] = node[donorPoint]->GetBetaInc2();
					Buffer_Receive_DensityInc[iVertex] = node[donorPoint]->GetDensityInc();
				}
				
				Buffer_Receive_Neighbor[iVertex] = geometry->node[donorPoint]->GetnPoint();
			}
#else
			
			unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
			unsigned long nBuffer_Scalar = geometry->nVertex[val_marker];
			
			MPI::COMM_WORLD.Recv(&Buffer_Receive_U,nBuffer_Vector,MPI::DOUBLE,receive_from, 0);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Recv(&Buffer_Receive_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Recv(&Buffer_Receive_Sensor,nBuffer_Scalar,MPI::DOUBLE,receive_from, 2);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_Lambda,nBuffer_Scalar,MPI::DOUBLE,receive_from, 3);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_Neighbor,nBuffer_Scalar,MPI::UNSIGNED_SHORT,receive_from, 4);
			
			if (incompressible) {
				MPI::COMM_WORLD.Recv(&Buffer_Receive_BetaInc2,nBuffer_Scalar,MPI::DOUBLE,receive_from, 5);
				MPI::COMM_WORLD.Recv(&Buffer_Receive_DensityInc,nBuffer_Scalar,MPI::DOUBLE,receive_from, 6);
			}
#endif
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				
				iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
				iPeriodic_Index = geometry->vertex[val_marker][iVertex]->GetRotation_Type();

				/*--- Retrieve the supplied periodic information. ---*/
				angles = config->GetPeriodicRotation(iPeriodic_Index);

				/*--- Store angles separately for clarity. ---*/
				theta    = angles[0];   phi    = angles[1]; psi    = angles[2];
				cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
				sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);

				/*--- Compute the rotation matrix. Note that the implicit
				 ordering is rotation about the x-axis, y-axis,
				 then z-axis. Note that this is the transpose of the matrix
				 used during the preprocessing stage. ---*/
				rotMatrix[0][0] = cosPhi*cosPsi;
				rotMatrix[0][1] = cosPhi*sinPsi;
				rotMatrix[0][2] = -sinPhi;
				
				rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
				rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
				rotMatrix[1][2] = sinTheta*cosPhi;
				
				rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
				rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
				rotMatrix[2][2] = cosTheta*cosPhi;
				
				/*--- Copy solution before performing transformation. ---*/
				for (iVar = 0; iVar < nVar; iVar++) {
					newSolution[iVar] = Buffer_Receive_U[iVertex][iVar];
				}

				/*--- Need to rotate the momentum components. ---*/
				if (nDim == 2) {
					newSolution[1] = rotMatrix[0][0]*Buffer_Receive_U[iVertex][1] + rotMatrix[0][1]*Buffer_Receive_U[iVertex][2];
					newSolution[2] = rotMatrix[1][0]*Buffer_Receive_U[iVertex][1] + rotMatrix[1][1]*Buffer_Receive_U[iVertex][2];
				} else {
					newSolution[1] = rotMatrix[0][0]*Buffer_Receive_U[iVertex][1] + rotMatrix[0][1]*Buffer_Receive_U[iVertex][2] + rotMatrix[0][2]*Buffer_Receive_U[iVertex][3];
					newSolution[2] = rotMatrix[1][0]*Buffer_Receive_U[iVertex][1] + rotMatrix[1][1]*Buffer_Receive_U[iVertex][2] + rotMatrix[1][2]*Buffer_Receive_U[iVertex][3];
					newSolution[3] = rotMatrix[2][0]*Buffer_Receive_U[iVertex][1] + rotMatrix[2][1]*Buffer_Receive_U[iVertex][2] + rotMatrix[2][2]*Buffer_Receive_U[iVertex][3];
				}

				for (iVar = 0; iVar < nVar; iVar++) {
					node[iPoint]->SetSolution(iVar, newSolution[iVar]);
					if (val_mesh == MESH_0) node[iPoint]->SetUndivided_Laplacian(iVar, Buffer_Receive_Undivided_Laplacian[iVertex][iVar]);
				}
				if (val_mesh == MESH_0) node[iPoint]->SetSensor(Buffer_Receive_Sensor[iVertex]);
				node[iPoint]->SetLambda(Buffer_Receive_Lambda[iVertex]);
				
				if (incompressible) {
					node[iPoint]->SetBetaInc2(Buffer_Receive_BetaInc2[iVertex]);
					node[iPoint]->SetDensityInc(Buffer_Receive_DensityInc[iVertex]);
				}
				
				geometry->node[iPoint]->SetnPoint(Buffer_Receive_Neighbor[iVertex]);
				
			}
		}			
		
	}
}

void CEulerSolution::BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
																		 unsigned short val_marker, unsigned short val_mesh) {
	
	unsigned short iVar;
	unsigned long iVertex, iPoint;
	short SendRecv = config->GetMarker_All_SendRecv(val_marker);
	
	/*--- Receive information  ---*/
	if (SendRecv < 0) {
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
			iPoint = geometry->vertex[val_marker][iVertex]->GetNode();	
			node[iPoint]->Set_ResConv_Zero();
			node[iPoint]->Set_ResVisc_Zero();
			node[iPoint]->SetResidualZero();
			node[iPoint]->SetTruncationErrorZero();
			if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT)
				for (iVar = 0; iVar < nVar; iVar++)
					Jacobian.DeleteValsRowi(iPoint*nVar+iVar);
		}
	}

}

void CEulerSolution::BC_Custom(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	double *U_domain, press, tempe, velo[3], b, rho, L, V, psi, Vpsi, theta, xx, yy; 
	U_domain = new double[nVar];
	
	/*-------- Set Problem Parameters for Ringleb flow -------------*/
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);	
	unsigned short imax = 500;
	double tol = 1.e-10;
	double Vmin = 0.;
	double Vmax = 2.;
	double Vp = 0.5*( Vmin + Vmax );
	
	/*--- Loop over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Copy the solution on the boundary ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Solution[iVar] = node[iPoint]->GetSolution(iVar);
			
			xx = geometry->node[iPoint]->GetCoord(0); 
			yy = geometry->node[iPoint]->GetCoord(1);
			
			/*--- Iterate to find solution (Adapted from Matsuzaka) ---*/
			for(unsigned short i = 1; i <= imax; i++) {
				b  = sqrt(1. - 0.2*Vp*Vp);
				rho = pow(b,5);
				L  = 1/b + 1/(3*pow(b,3)) + 1/(5*pow(b,5)) - 0.5*log((1+b)/(1-b));
				V   = sqrt(1/sqrt(fabs((xx-0.5*L)*(xx-0.5*L)+yy*yy) )/2/rho);
				if (fabs(V-Vp) < tol) break;
				if (i == imax)
					cout<<"did not converge... "<<fabs(V-Vp)<<" "<<xx<<" "<<yy<<endl;
				Vp = V;
			}
			
			b = sqrt(1. - 0.2*V*V);
			rho = pow(b,5);
			press = rho*b*b/Gamma;
			L = 1/b + 1/(3*pow(b,3)) + 1/(5*pow(b,5))-0.5*log((1+b)/(1-b));
			psi = sqrt(0.5/V/V-(xx-0.5*L)*rho);    Vpsi= sqrt(0.5-V*V*(xx-0.5*L)*rho);
			if ( fabs( Vpsi - 1. ) < tol  ) theta = 0.5*PI_NUMBER;
			else theta = asin( psi*V );
			velo[0] = -V*cos(theta); velo[1] = -V*sin(theta); velo[2] = 0.;
			tempe = press/rho;
			
			/*--- Load up the solution vector ---*/
			Solution[0] = rho; 
			for (iVar = 1; iVar< nVar-1; iVar++) 
				Solution[iVar] = rho*velo[iVar-1]; 
			Solution[nVar-1] = press/(Gamma-1.)+0.5*rho*(velo[0]*velo[0]+velo[1]*velo[1]); 
			
			/*--- Weak boundary imposition ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar]= node[iPoint]->GetSolution(iVar);
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			solver->SetNormal(Vector);
			
			/*--- Compute the residual using an upwind scheme ---*/
			solver->SetConservative(U_domain, Solution);
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddRes_Conv(Residual);
			
			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			
		}
	}
}

CNSSolution::CNSSolution(void) : CEulerSolution() { }

CNSSolution::CNSSolution(CGeometry *geometry, CConfig *config) : CEulerSolution() {
	unsigned long iPoint, index;
	unsigned short iVar, iDim, iMarker;
	ifstream restart_file;
	string mesh_filename, text_line;
	
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();

	/*--- Set the gamma value ---*/
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Define geometry constants in the solver structure ---*/
	nDim = geometry->GetnDim();
	if (incompressible) nVar = nDim + 1;
	else nVar = nDim + 2;
	nMarker = config->GetnMarker_All();
	nPoint = geometry->GetnPoint();
	node = new CVariable*[geometry->GetnPoint()];

	/*--- Define some auxiliar vector related with the residual ---*/
	Residual = new double[nVar]; Residual_Max = new double[nVar];
	Residual_i = new double[nVar]; Residual_j = new double[nVar];
	Res_Conv = new double[nVar]; Res_Visc = new double[nVar];
	
	/*--- Define some auxiliar vector related with the solution ---*/
	Solution = new double[nVar];
	Solution_i = new double[nVar]; Solution_j = new double[nVar];
	
	/*--- Define some auxiliar vector related with the primitive variables ---*/
	PrimVar = new double[nVar];
	PrimVar_i = new double[nVar]; PrimVar_j = new double[nVar];
	
	/*--- Define some auxiliar vector related with the geometry ---*/
	Vector = new double[nDim];
	Vector_i = new double[nDim]; Vector_j = new double[nDim];

	/*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/
	if ((config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED) ||
		(config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTRED) ) {
		p1_Und_Lapl = new double [geometry->GetnPoint()]; 
		p2_Und_Lapl = new double [geometry->GetnPoint()]; 
	}
	
	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
		
		/*--- Point to point Jacobians ---*/
		Jacobian_i = new double* [nVar]; 
		Jacobian_j = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new double [nVar]; 
			Jacobian_j[iVar] = new double [nVar]; 
		}
		
		/*--- Initialization of the structure of the whole Jacobian ---*/
		InitializeJacobianStructure(geometry, config);
		xsol = new double [geometry->GetnPoint()*nVar];
		rhs = new double [geometry->GetnPoint()*nVar];
	}

	/*--- Computation of gradients by least squares ---*/
	if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) || 
		(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) {
		
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim]; 
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
		
		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nVar+1]; 
		for (iVar = 0; iVar < nVar+1; iVar++)
			cvector[iVar] = new double [nDim];
	}

	/*--- Inviscid forces definition and coefficient in all the markers ---*/
	CPressure = new double* [config->GetnMarker_All()];
	for (iMarker=0; iMarker<config->GetnMarker_All(); iMarker++)
		CPressure[iMarker] = new double [geometry->nVertex[iMarker]];
	
	ForceInviscid = new double[nDim];
	MomentInviscid = new double[3];
	CDrag_Inv = new double[config->GetnMarker_All()];
	CLift_Inv = new double[config->GetnMarker_All()];
	CSideForce_Inv = new double[config->GetnMarker_All()];
	CPress_Inv = new double[config->GetnMarker_All()];
	CMx_Inv = new double[config->GetnMarker_All()];
	CMy_Inv = new double[config->GetnMarker_All()];
	CMz_Inv = new double[config->GetnMarker_All()];
	CEff_Inv = new double[config->GetnMarker_All()];
	CEquivArea_Inv = new double[config->GetnMarker_All()];
	CNearFieldPress_Inv = new double[config->GetnMarker_All()];
	Total_CDrag = 0.0; Total_CLift = 0.0; Total_CSideForce = 0.0;
	Total_CPress = 0.0; Total_CMx = 0.0; Total_CMy = 0.0; Total_CMz = 0.0;
	Total_CEff = 0.0; Total_CEquivArea = 0.0; Total_CNearFieldPress = 0.0;
	
	/*--- Viscous forces definition and coefficient in all the markers ---*/
	CSkinFriction = new double* [config->GetnMarker_All()];
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		CSkinFriction[iMarker] = new double [geometry->nVertex[iMarker]];	
	ForceViscous = new double[nDim];
	MomentViscous = new double[nDim];
	CDrag_Visc = new double[config->GetnMarker_All()];
	CLift_Visc = new double[config->GetnMarker_All()];
	CMx_Visc = new double[config->GetnMarker_All()];
	CMy_Visc = new double[config->GetnMarker_All()];
	CMz_Visc = new double[config->GetnMarker_All()];
	CEff_Visc = new double[config->GetnMarker_All()];

  /*--- Read farfield conditions from config ---*/
  Density_Inf   = config->GetDensity_FreeStreamND();
  Pressure_Inf  = config->GetPressure_FreeStreamND();
  Velocity_Inf  = config->GetVelocity_FreeStreamND();
	Energy_Inf    = config->GetEnergy_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  Prandtl_Lam   = config->GetPrandtl_Lam();
  Prandtl_Turb  = config->GetPrandtl_Turb();
	

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		/*--- For a rotating frame, set the velocity due to rotation
		 at each point just once, and store it in the geometry class. ---*/
		if (rotating_frame) {
			double RotVel[3], distance[3];
			double *coord = geometry->node[iPoint]->GetCoord();
			double *axis  = config->GetRotAxisOrigin();
			double *omega = config->GetOmega_FreeStreamND();
			double Lref   = config->GetLength_Ref();
			
			/*--- Calculate non-dim distance fron rotation center ---*/
			distance[0] = (coord[0]-axis[0])/Lref;
			distance[1] = (coord[1]-axis[1])/Lref;
			distance[2] = (coord[2]-axis[2])/Lref;
			
			/*--- Calculate the angular velocity as omega X r ---*/
			RotVel[0] = omega[1]*(distance[2]) - omega[2]*(distance[1]);
			RotVel[1] = omega[2]*(distance[0]) - omega[0]*(distance[2]);
			RotVel[2] = omega[0]*(distance[1]) - omega[1]*(distance[0]);
			
			geometry->node[iPoint]->SetRotVel(RotVel);
		}
		
	}
	
	/*--- Restart the solution from file information ---*/
	if (!restart) {
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
	}
	else {
	
		/*--- Restart the solution from file information ---*/
		mesh_filename = config->GetSolution_FlowFileName();
		restart_file.open(mesh_filename.data(), ios::in);

		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
			cout << "There is no flow restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get(); exit(1);
		}

		/*--- Read the restart file ---*/
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			if (incompressible) {
				if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2];
				if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
			}
			else {
				if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
				if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
			}
			node[iPoint] = new CNSVariable(Solution, nDim, nVar, config);
		}
		
		/*--- Close the restart file ---*/
		restart_file.close();
	}
	
	/*--- For incompressible solver set the initial values for the density and viscosity,
	 unless a twophase problem, this must be constant during the computation ---*/
	if (incompressible)
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		node[iPoint]->SetDensityInc(Density_Inf);
		node[iPoint]->SetLaminarViscosityInc(Viscosity_Inf);
	}
	
	/*--- Define solver parameters needed for execution of destructor ---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED || (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTRED))
		space_centered = true;
	else
		space_centered = false;
	if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT)
		euler_implicit = true;
	else
		euler_implicit = false;
	if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||	(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES))
		least_squares = true;
	else
		least_squares = false;
}

CNSSolution::~CNSSolution(void) {
	unsigned short iVar, iDim;
	
	delete [] Residual; delete [] Residual_Max;
	delete [] Residual_i; delete [] Residual_j;
	delete [] Res_Conv; delete [] Res_Visc;
	delete [] Solution; delete [] Solution_i; delete [] Solution_j;
	delete [] PrimVar; delete [] PrimVar_i; delete [] PrimVar_j;
	delete [] Vector_i; delete [] Vector_j;
	delete [] p1_Und_Lapl; delete [] p2_Und_Lapl;
	delete [] xsol; delete [] rhs;
	delete [] ForceInviscid; delete [] MomentInviscid;
	delete [] CDrag_Inv; delete [] CLift_Inv;
	delete [] CMx_Inv; delete [] CMy_Inv; delete [] CMz_Inv;
	delete [] ForceViscous; delete [] MomentViscous;
	delete [] CDrag_Visc; delete [] CLift_Visc;
	delete [] CMx_Visc; delete [] CMy_Visc; delete [] CMz_Visc;
	delete [] Velocity_Inf;
	
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Jacobian_i[iVar];
		delete [] Jacobian_j[iVar];
	}
	delete [] Jacobian_i; delete [] Jacobian_j;
	
	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;
	
	for (iVar = 0; iVar < nVar+1; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
	
	/*	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
	 delete [] CPressure[iMarker];
	 delete [] CPressure; */
	
	/*	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
	 delete [] CSkinFriction[iMarker];
	 delete [] CSkinFriction; */
	
	/*	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
	 delete [] node[iPoint];
	 delete [] node; */
}

void CNSSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
	double Gas_Constant = config->GetGas_Constant();

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		/*--- Compute squared velocity, sound velocity, pressure, enthalpy, vorticity, and temperature ---*/
		if (incompressible) {			
			node[iPoint]->SetPressureInc();
			node[iPoint]->SetVelocityInc2();
			node[iPoint]->SetBetaInc2(config);
		}
		else {
			node[iPoint]->SetVelocity2();
			node[iPoint]->SetPressure(Gamma, geometry->node[iPoint]->GetCoord());
      node[iPoint]->SetSoundSpeed(Gamma);
			node[iPoint]->SetEnthalpy();
			node[iPoint]->SetVorticity();
			node[iPoint]->SetTemperature(Gas_Constant);
			node[iPoint]->SetLaminarViscosity(config);	
		}
		
		/*--- Compute laminar, turbulent, and termal coefficient ---*/

		if (config->GetKind_Turb_Model() != NONE) 
			node[iPoint]->SetEddyViscosity(config->GetKind_Turb_Model(), 
																		 solution_container[TURB_SOL]->node[iPoint]->GetSolution());
		else node[iPoint]->SetEddyViscosity(NONE,NULL);
		
		node[iPoint]->SetThermalCoeff(Gamma, Gas_Constant);
		
		/*--- Initialize the convective and viscous residual vector ---*/
		node[iPoint]->Set_ResConv_Zero();
		if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit)
			node[iPoint]->Set_ResVisc_Zero();
	}
	
	/*--- Initialize the jacobian matrices ---*/
	if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT)
		Jacobian.SetValZero();
}

void CNSSolution::SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {
	double Mean_BetaInc2, *Normal, Area, dV, Mean_SoundSpeed, Mean_ProjVel, Lambda, Local_Delta_Time, Local_Delta_Time_Visc, 
	Global_Delta_Time = 1E6, Mean_LaminarVisc, Mean_EddyVisc, Mean_Density, Lambda_1, Lambda_2, K_v = 0.25;
	unsigned long iEdge, iVertex, iPoint = 0, jPoint = 0;
	unsigned short iDim, iMarker;
	double ProjVel, ProjVel_i, ProjVel_j;

	bool rotating_frame = config->GetRotating_Frame();
  bool incompressible = config->GetIncompressible();
	double Gas_Constant = config->GetGas_Constant();

	Min_Delta_Time = 1.E6;
	Max_Delta_Time = 0.0;
	
	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed and viscosity ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		node[iPoint]->SetMax_Lambda_Inv(0.0);
		node[iPoint]->SetMax_Lambda_Visc(0.0); 
		if (incompressible) {
			node[iPoint]->SetVelocityInc2();
			node[iPoint]->SetBetaInc2(config);
		}
		else {
			node[iPoint]->SetVelocity2();
      node[iPoint]->SetPressure(Gamma, geometry->node[iPoint]->GetCoord());
			node[iPoint]->SetSoundSpeed(Gamma);
		}

		node[iPoint]->SetTemperature(Gas_Constant);
		node[iPoint]->SetLaminarViscosity(config);
		if (config->GetKind_Turb_Model() != NONE)
			node[iPoint]->SetEddyViscosity(config->GetKind_Turb_Model(),solution_container[TURB_SOL]->node[iPoint]->GetSolution());
		else
			node[iPoint]->SetEddyViscosity(NONE, NULL);
	}
	
	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) { 
		
		/*--- Point identification, Normal vector and area ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); 
		jPoint = geometry->edge[iEdge]->GetNode(1);
		Normal = geometry->edge[iEdge]->GetNormal();
		Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
		
		/*--- Mean Values ---*/
		if (incompressible) {
			Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVelInc(Normal) + node[jPoint]->GetProjVelInc(Normal));
			Mean_BetaInc2 = 0.5 * (node[iPoint]->GetBetaInc2() + node[jPoint]->GetBetaInc2());
			Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + Mean_BetaInc2*Area*Area); 
		}
		else {
			Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
			Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
		}
		
		/*--- Contribution from rotational framework ---*/
		if (rotating_frame) {
			double *RotVel_i = geometry->node[iPoint]->GetRotVel();
			double *RotVel_j = geometry->node[jPoint]->GetRotVel();
			ProjVel_i = 0.0; ProjVel_j =0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				ProjVel_i += RotVel_i[iDim]*Normal[iDim];
				ProjVel_j += RotVel_j[iDim]*Normal[iDim];
			}
			Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j) ;
		}
		
		/*--- Inviscid contribution ---*/
		Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed ;
		node[iPoint]->AddMax_Lambda_Inv(Lambda);
		node[jPoint]->AddMax_Lambda_Inv(Lambda);
		
		/*--- Viscous contribution ---*/
		if (incompressible) {
			Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosityInc() + node[jPoint]->GetLaminarViscosityInc());
			Mean_EddyVisc    = 0.5*(node[iPoint]->GetEddyViscosity() + node[jPoint]->GetEddyViscosity());
			Mean_Density     = 0.5*(node[iPoint]->GetDensityInc() + node[jPoint]->GetDensityInc());
		}
		else {
			Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
			Mean_EddyVisc    = 0.5*(node[iPoint]->GetEddyViscosity() + node[jPoint]->GetEddyViscosity());
			Mean_Density     = 0.5*(node[iPoint]->GetSolution(0) + node[jPoint]->GetSolution(0));
		}
		
		Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
		Lambda_2 = (1.0 + (Prandtl_Lam/Prandtl_Turb)*(Mean_EddyVisc/Mean_LaminarVisc))*(Gamma*Mean_LaminarVisc/Prandtl_Lam);
		Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
		
		node[iPoint]->AddMax_Lambda_Visc(Lambda);
		node[jPoint]->AddMax_Lambda_Visc(Lambda);
	}
	
	/*--- Loop boundary edges ---*/
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) { 
		for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			
			/*--- Point identification, Normal vector and area ---*/
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
			
			/*--- Mean Values ---*/
			if (incompressible) {
				Mean_ProjVel = node[iPoint]->GetProjVelInc(Normal);
				Mean_BetaInc2 = node[iPoint]->GetBetaInc2();
				Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + Mean_BetaInc2*Area*Area); 
			}
			else {
				Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
				Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
			}
			
			if (rotating_frame) {
				double *RotVel = geometry->node[iPoint]->GetRotVel();
				ProjVel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					ProjVel += RotVel[iDim]*Normal[iDim];
				Mean_ProjVel -= ProjVel;
			}
			
			/*--- Inviscid contribution ---*/
			Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
			node[iPoint]->AddMax_Lambda_Inv(Lambda);            
			
			/*--- Viscous contribution ---*/
			if (incompressible) {
				Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosityInc() + node[jPoint]->GetLaminarViscosityInc());
				Mean_EddyVisc    = 0.5*(node[iPoint]->GetEddyViscosity() + node[jPoint]->GetEddyViscosity());
				Mean_Density     = 0.5*(node[iPoint]->GetDensityInc() + node[jPoint]->GetDensityInc());
			}
			else {
				Mean_LaminarVisc = node[iPoint]->GetLaminarViscosity();
				Mean_EddyVisc    = node[iPoint]->GetEddyViscosity();
				Mean_Density     = node[iPoint]->GetSolution(0);
			}
			
			Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
			Lambda_2 = (1.0 + (Prandtl_Lam/Prandtl_Turb)*(Mean_EddyVisc/Mean_LaminarVisc))*(Gamma*Mean_LaminarVisc/Prandtl_Lam);
			Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
			node[iPoint]->AddMax_Lambda_Visc(Lambda);
		}
	}
	
	/*--- Each element uses their own speed ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		dV = geometry->node[iPoint]->GetVolume();
		Local_Delta_Time = config->GetCFL(iMesh)*dV/node[iPoint]->GetMax_Lambda_Inv();
		Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*dV*dV/node[iPoint]->GetMax_Lambda_Visc();
		Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
		Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
		Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
		Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
		node[iPoint]->SetDelta_Time(Local_Delta_Time);
	}
	
	Sum_Delta_Time += Min_Delta_Time;

	/*--- For exact time solution use the minimum delta time of the whole mesh ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING)
		for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
			node[iPoint]->SetDelta_Time(Global_Delta_Time);
}

void CNSSolution::Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
								   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iPoint, jPoint, iEdge;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();

	if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit) {
		
		/*--- Compute gradients (not in the preprocessing, because we want to be sure that we have
		 the gradient of the primitive varibles, not the conservative) ---*/
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetPrimVar_Gradient_GG(geometry, config); 
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetPrimVar_Gradient_LS(geometry, config);
		
		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
			
			/*--- Points, coordinates and normal vector in edge ---*/
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);
			solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
			solver->SetNormal(geometry->edge[iEdge]->GetNormal());
			
			/*--- Conservative variables, and primitive Variables w/o reconstruction ---*/
			solver->SetConservative(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
			solver->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(),node[jPoint]->GetGradient_Primitive());
			
			/*--- Laminar and eddy viscosity ---*/
			if (incompressible) {
				solver->SetDensityInc(node[iPoint]->GetDensityInc(), node[jPoint]->GetDensityInc());
				solver->SetLaminarViscosity(node[iPoint]->GetLaminarViscosityInc(), node[jPoint]->GetLaminarViscosityInc());
			}
			else
				solver->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[jPoint]->GetLaminarViscosity());
				
			solver->SetEddyViscosity(node[iPoint]->GetEddyViscosity(), node[jPoint]->GetEddyViscosity());
			
			/*--- Compute and update residual ---*/
			solver->SetResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
			node[iPoint]->SubtractRes_Visc(Res_Visc);
			node[jPoint]->AddRes_Visc(Res_Visc);
			
			/*--- Implicit part ---*/
			if (implicit) {
				Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
				Jacobian.SubtractBlock(iPoint,jPoint,Jacobian_j);
				Jacobian.AddBlock(jPoint,iPoint,Jacobian_i);
				Jacobian.AddBlock(jPoint,jPoint,Jacobian_j);
			}
		}
	}
}

void CNSSolution::Viscous_Forces(CGeometry *geometry, CConfig *config) {
	unsigned long iVertex, iPoint;
	unsigned short Boundary, Monitoring, iMarker, iDim, jDim;
	double **Tau, Delta, Viscosity, **Grad_PrimVar, div_vel, *Normal, *TauElem;
  double dist[3], *Coord, *Kappa, *TauTangent, Area, WallShearStress, TauNormal;
  double factor, RefVel2, RefDensity;
	double Alpha        = config->GetAoA()*PI_NUMBER/180.0;
	double Beta         = config->GetAoS()*PI_NUMBER/180.0;
	double RefAreaCoeff = config->GetRefAreaCoeff();
	double *Origin      = config->GetRefOriginMoment();
  bool rotating_frame = config->GetRotating_Frame();

  /*--- If this is a rotating frame problem, use the specified speed
        for computing the force coefficients. Otherwise, use the 
        freestream values which is the standard convention. ---*/
  if (rotating_frame) {
    RefVel2     = config->GetRotVel_Ref()*config->GetRotVel_Ref();
    RefDensity  = config->GetDensity_FreeStreamND();
  } else {
    double *Velocity_Inf = config->GetVelocity_FreeStreamND();
    RefVel2 = 0.0; 
    for (iDim = 0; iDim < nDim; iDim++) 
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    RefDensity  = config->GetDensity_FreeStreamND();
  }

	/*-- Initialization --*/
	AllBound_CDrag_Visc = 0.0; AllBound_CLift_Visc = 0.0;
	AllBound_CMx_Visc = 0.0; AllBound_CMy_Visc = 0.0; AllBound_CMz_Visc = 0.0;
	AllBound_CEff_Visc = 0.0;
	/*--- Vector and variables initialization ---*/
	Kappa      = new double [nDim];
	TauElem    = new double [nDim];
	TauTangent = new double [nDim];
	Tau        = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Tau[iDim]   = new double [nDim]; 

	/*--- Loop over the Navier-Stokes markers ---*/
	factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);
		
		if (Boundary == NO_SLIP_WALL) {
			
			for (iDim = 0; iDim < nDim; iDim++) ForceViscous[iDim] = 0.0;
			MomentViscous[0] = 0.0; MomentViscous[1] = 0.0; MomentViscous[2] = 0.0;

			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					Coord = geometry->node[iPoint]->GetCoord();
					for (iDim = 0; iDim < nDim; iDim++) dist[iDim] = Coord[iDim] - Origin[iDim];

					Viscosity = node[iPoint]->GetLaminarViscosity();
					Grad_PrimVar = node[iPoint]->GetGradient_Primitive();
					
					div_vel = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						div_vel += Grad_PrimVar[iDim+1][iDim];
					
					for (iDim = 0; iDim < nDim; iDim++)
						for (jDim = 0 ; jDim < nDim; jDim++) {
							Delta = 0.0; if (iDim == jDim) Delta = 1.0;
							Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + 
														 Grad_PrimVar[iDim+1][jDim]) - TWO3*Viscosity*div_vel*Delta;
						}
					
					/*--- Compute viscous foces ---*/
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
					for (iDim = 0; iDim < nDim; iDim++) Kappa[iDim] = Normal[iDim]/Area;	
					
					for (iDim = 0; iDim < nDim; iDim++) {
						TauElem[iDim] = 0.0;
						for (jDim = 0; jDim < nDim; jDim++)
							TauElem[iDim] += Tau[iDim][jDim]*Kappa[jDim] ;
					}
					
					if ((geometry->node[iPoint]->GetDomain()) && (Monitoring == YES)) {
						for (iDim = 0; iDim < nDim; iDim++)
							ForceViscous[iDim] += TauElem[iDim]*Area*factor;
						
						/*--- Moment with respect to the reference axis ---*/
						if (iDim == 3) MomentViscous[0] += ForceViscous[2]*dist[1] - ForceViscous[1]*dist[2];
						if (iDim == 3) MomentViscous[1] += ForceViscous[0]*dist[2] - ForceViscous[2]*dist[0];
						MomentViscous[2] += ForceViscous[1]*dist[0] - ForceViscous[0]*dist[1];
					}

					/*--- Compute wall shear stress, and skin friction coefficient ---*/
					TauNormal = 0.0; 
					for (iDim = 0; iDim < nDim; iDim++) 
						TauNormal += TauElem[iDim] * Kappa[iDim];
					
					for (iDim = 0; iDim < nDim; iDim++) 
						TauTangent[iDim] = TauElem[iDim] - TauNormal * Kappa[iDim];
					
					WallShearStress = 0.0; 
					for (iDim = 0; iDim < nDim; iDim++) 
						WallShearStress += TauTangent[iDim]*TauTangent[iDim]; 
					
					/*--- Note that the wall Shear Stress is just mu(delta u/delta y)---*/

					CSkinFriction[iMarker][iVertex] = sqrt(WallShearStress) / (0.5*RefDensity*RefVel2);
					

				}
			}
			
			/*--- Transform ForceInviscid into CLift and CDrag ---*/
			if  (Monitoring == YES) {
				if (nDim == 2) {
					CDrag_Visc[iMarker] =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
					CLift_Visc[iMarker] = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
					CMx_Visc[iMarker] = 0.0;
					CMy_Visc[iMarker] = 0.0;
					CMz_Visc[iMarker] = MomentViscous[2];
					CEff_Visc[iMarker] = CLift_Visc[iMarker]/CDrag_Visc[iMarker];
				}
				if (nDim == 3) {
					CDrag_Visc[iMarker] =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
					CLift_Visc[iMarker] = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
					CMx_Visc[iMarker] = MomentViscous[0];
					CMy_Visc[iMarker] = MomentViscous[1];
					CMz_Visc[iMarker] = MomentViscous[2];
					CEff_Visc[iMarker] = CLift_Visc[iMarker]/CDrag_Visc[iMarker];
				}
				
				AllBound_CDrag_Visc += CDrag_Visc[iMarker];
				AllBound_CLift_Visc += CLift_Visc[iMarker];
				AllBound_CMx_Visc += CMx_Visc[iMarker];
				AllBound_CMy_Visc += CMy_Visc[iMarker];
				AllBound_CMz_Visc += CMz_Visc[iMarker];
				AllBound_CEff_Visc += CEff_Visc[iMarker];
			}
		}
	}

	Total_CDrag += AllBound_CDrag_Visc;
	Total_CLift += AllBound_CLift_Visc;
	Total_CMx += AllBound_CMx_Visc;
	Total_CMy += AllBound_CMy_Visc;
	Total_CMz += AllBound_CMz_Visc;	
	Total_CEff = Total_CLift/Total_CDrag;

//cout << "Total_CEff: " << Total_CEff << endl;
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Tau[iDim];
	delete [] Tau;
	delete [] Kappa;
	delete [] TauTangent;
	delete [] TauElem;
}

void CNSSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, total_index;
	unsigned short iVar, iDim;
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Vector --> Velocity_corrected ---*/
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
			
			/*--- Set the residual, truncation error and velocity value ---*/
			node[iPoint]->SetVelocity_Old(Vector);
			node[iPoint]->SetVel_ResConv_Zero();
			node[iPoint]->SetVel_ResVisc_Zero();
			node[iPoint]->SetVelTruncationErrorZero();
			
			/*--- Only change velocity-rows of the Jacobian (includes 1 in the diagonal) ---*/
			if (implicit)
				for (iVar = 1; iVar <= nDim; iVar++) {
					total_index = iPoint*nVar+iVar;
					Jacobian.DeleteValsRowi(total_index);
				}
		}
	}
}

void CNSSolution::BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
									 unsigned short val_marker, unsigned short val_mesh) {
	
#ifndef NO_MPI
	unsigned short iVar;
	unsigned long iVertex, iPoint;
	double *Conserv_Var, *Conserv_Undivided_Laplacian = NULL, **Conserv_Grad = NULL, **PrimVar_Grad, *Grad_Limit = NULL;
	
	short SendRecv = config->GetMarker_All_SendRecv(val_marker);
	bool incompressible = config->GetIncompressible();
		
	/*--- Send information  ---*/
	if (SendRecv > 0) {
		/*--- Upwind scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {
			
			/*--- Inviscid part ---*/
			double Buffer_Send_U[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Ux[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Uy[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Uz[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Limit[geometry->nVertex[val_marker]][nVar];
			
			/*--- Viscous part ---*/
			double Buffer_Send_LaminarViscosity[geometry->nVertex[val_marker]];
			double Buffer_Send_EddyViscosity[geometry->nVertex[val_marker]];
			double Buffer_Send_Vx[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Vy[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Vz[geometry->nVertex[val_marker]][nVar];
			
			/*--- Dimensionalization ---*/
			unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
			unsigned long nBuffer_Scalar = geometry->nVertex[val_marker];
			int send_to = SendRecv;
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
				Conserv_Var = node[iPoint]->GetSolution();
				PrimVar_Grad = node[iPoint]->GetGradient_Primitive();
				if (val_mesh == MESH_0) Conserv_Grad = node[iPoint]->GetGradient();
				if ((val_mesh == MESH_0) && (config->GetKind_SlopeLimit_Flow() != NONE)) 
					Grad_Limit = node[iPoint]->GetLimiter();
				
				for (iVar = 0; iVar < nVar; iVar++) {
					Buffer_Send_U[iVertex][iVar] = Conserv_Var[iVar];
					if (val_mesh == MESH_0) {
						Buffer_Send_Ux[iVertex][iVar] = Conserv_Grad[iVar][0];					
						Buffer_Send_Uy[iVertex][iVar] = Conserv_Grad[iVar][1];
						if (nDim == 3) Buffer_Send_Uz[iVertex][iVar] = Conserv_Grad[iVar][2];
						if (config->GetKind_SlopeLimit_Flow() != NONE) 
							Buffer_Send_Limit[iVertex][iVar] = Grad_Limit[iVar];
					}
					Buffer_Send_Vx[iVertex][iVar] = PrimVar_Grad[iVar][0];
					Buffer_Send_Vy[iVertex][iVar] = PrimVar_Grad[iVar][1];					
					if (nDim == 3) Buffer_Send_Vz[iVertex][iVar] = PrimVar_Grad[iVar][2];
				}
				Buffer_Send_LaminarViscosity[iVertex] = node[iPoint]->GetLaminarViscosity();
				Buffer_Send_EddyViscosity[iVertex] = node[iPoint]->GetEddyViscosity();
			}
			MPI::COMM_WORLD.Bsend(&Buffer_Send_U,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
			if (val_mesh == MESH_0) {
				MPI::COMM_WORLD.Bsend(&Buffer_Send_Ux,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
				MPI::COMM_WORLD.Bsend(&Buffer_Send_Uy,nBuffer_Vector,MPI::DOUBLE,send_to, 2);
				if (nDim == 3) MPI::COMM_WORLD.Bsend(&Buffer_Send_Uz,nBuffer_Vector,MPI::DOUBLE,send_to, 3);
				if (config->GetKind_SlopeLimit_Flow() != NONE) MPI::COMM_WORLD.Bsend(&Buffer_Send_Limit,nBuffer_Vector,MPI::DOUBLE,send_to, 4);
			}
			MPI::COMM_WORLD.Bsend(&Buffer_Send_LaminarViscosity,nBuffer_Scalar,MPI::DOUBLE,send_to, 5);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_EddyViscosity,nBuffer_Scalar,MPI::DOUBLE,send_to, 6);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_Vx,nBuffer_Vector,MPI::DOUBLE,send_to, 7);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_Vy,nBuffer_Vector,MPI::DOUBLE,send_to, 8);
			if (nDim == 3) MPI::COMM_WORLD.Bsend(&Buffer_Send_Vz,nBuffer_Vector,MPI::DOUBLE,send_to, 9);
		}
		
		/*--- Centered scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {
			
			/*--- Inviscid part ---*/
			double Buffer_Send_U[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Undivided_Laplacian[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Sensor[geometry->nVertex[val_marker]];
			double Buffer_Send_Lambda[geometry->nVertex[val_marker]];
			double Buffer_Send_BetaInc2[geometry->nVertex[val_marker]];
			double Buffer_Send_DensityInc[geometry->nVertex[val_marker]];
			double Buffer_Send_LaminarViscosityInc[geometry->nVertex[val_marker]];
			unsigned short Buffer_Send_Neighbor[geometry->nVertex[val_marker]];

			/*--- Viscous part ---*/
			double Buffer_Send_LaminarViscosity[geometry->nVertex[val_marker]];
			double Buffer_Send_EddyViscosity[geometry->nVertex[val_marker]];
			double Buffer_Send_Vx[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Vy[geometry->nVertex[val_marker]][nVar];
			double Buffer_Send_Vz[geometry->nVertex[val_marker]][nVar];

			/*--- Dimensionalization ---*/
			unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
			unsigned long nBuffer_Scalar = geometry->nVertex[val_marker];
			int send_to = SendRecv;
			
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
				Conserv_Var = node[iPoint]->GetSolution();
				PrimVar_Grad = node[iPoint]->GetGradient_Primitive();
				if (val_mesh == MESH_0) Conserv_Undivided_Laplacian = node[iPoint]->GetUnd_Lapl();
				for (iVar = 0; iVar < nVar; iVar++) {
					Buffer_Send_U[iVertex][iVar] = Conserv_Var[iVar];
					if (val_mesh == MESH_0) Buffer_Send_Undivided_Laplacian[iVertex][iVar] = Conserv_Undivided_Laplacian[iVar];
					Buffer_Send_Vx[iVertex][iVar] = PrimVar_Grad[iVar][0];
					Buffer_Send_Vy[iVertex][iVar] = PrimVar_Grad[iVar][1];					
					if (nDim == 3) Buffer_Send_Vz[iVertex][iVar] = PrimVar_Grad[iVar][2];
				}
				if (val_mesh == MESH_0) Buffer_Send_Sensor[iVertex] = node[iPoint]->GetSensor();
				Buffer_Send_Lambda[iVertex] = node[iPoint]->GetLambda();
				Buffer_Send_Neighbor[iVertex] = geometry->node[iPoint]->GetnPoint();
				
				if (incompressible) {
					Buffer_Send_BetaInc2[iVertex] = node[iPoint]->GetBetaInc2();
					Buffer_Send_DensityInc[iVertex] = node[iPoint]->GetDensityInc();
					Buffer_Send_LaminarViscosityInc[iVertex] = node[iPoint]->GetLaminarViscosityInc();
				}
				
				Buffer_Send_LaminarViscosity[iVertex] = node[iPoint]->GetLaminarViscosity();
				Buffer_Send_EddyViscosity[iVertex] = node[iPoint]->GetEddyViscosity();
			}
			
			MPI::COMM_WORLD.Bsend(&Buffer_Send_U,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Bsend(&Buffer_Send_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Bsend(&Buffer_Send_Sensor,nBuffer_Scalar,MPI::DOUBLE,send_to, 2);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_Lambda,nBuffer_Scalar,MPI::DOUBLE,send_to, 3);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_Neighbor,nBuffer_Scalar,MPI::UNSIGNED_SHORT,send_to, 4);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_LaminarViscosity,nBuffer_Scalar,MPI::DOUBLE,send_to, 5);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_EddyViscosity,nBuffer_Scalar,MPI::DOUBLE,send_to, 6);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_Vx,nBuffer_Vector,MPI::DOUBLE,send_to, 7);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_Vy,nBuffer_Vector,MPI::DOUBLE,send_to, 8);
			if (nDim == 3) MPI::COMM_WORLD.Bsend(&Buffer_Send_Vz,nBuffer_Vector,MPI::DOUBLE,send_to, 9);
			
			if (incompressible) {
				MPI::COMM_WORLD.Bsend(&Buffer_Send_BetaInc2,nBuffer_Scalar,MPI::DOUBLE,send_to, 10);
				MPI::COMM_WORLD.Bsend(&Buffer_Send_DensityInc,nBuffer_Scalar,MPI::DOUBLE,send_to, 11);
				MPI::COMM_WORLD.Bsend(&Buffer_Send_LaminarViscosityInc,nBuffer_Scalar,MPI::DOUBLE,send_to, 12);
			}
		}
	}
	
	/*--- Receive information  ---*/
	if (SendRecv < 0) {
		/*--- Upwind scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {
			
			/*--- Inviscid part ---*/
			double Buffer_Receive_U[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Ux[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Uy[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Uz[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Limit[geometry->nVertex[val_marker]][nVar];
			
			/*--- Viscous part ---*/
			double Buffer_Receive_LaminarViscosity[geometry->nVertex[val_marker]];
			double Buffer_Receive_EddyViscosity[geometry->nVertex[val_marker]];
			double Buffer_Receive_Vx[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Vy[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Vz[geometry->nVertex[val_marker]][nVar];
			
			/*--- Dimensionalization ---*/
			unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
			unsigned long nBuffer_Scalar = geometry->nVertex[val_marker];
			int receive_from = abs(SendRecv);
			
			MPI::COMM_WORLD.Recv(&Buffer_Receive_U,nBuffer_Vector,MPI::DOUBLE,receive_from, 0);
			if (val_mesh == MESH_0) {
				MPI::COMM_WORLD.Recv(&Buffer_Receive_Ux,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
				MPI::COMM_WORLD.Recv(&Buffer_Receive_Uy,nBuffer_Vector,MPI::DOUBLE,receive_from, 2);
				if (nDim == 3) MPI::COMM_WORLD.Recv(&Buffer_Receive_Uz,nBuffer_Vector,MPI::DOUBLE,receive_from, 3);
				if (config->GetKind_SlopeLimit_Flow() != NONE) MPI::COMM_WORLD.Recv(&Buffer_Receive_Limit,nBuffer_Vector,MPI::DOUBLE,receive_from, 4);
			}
			MPI::COMM_WORLD.Recv(&Buffer_Receive_LaminarViscosity,nBuffer_Scalar,MPI::DOUBLE,receive_from, 5);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_EddyViscosity,nBuffer_Scalar,MPI::DOUBLE,receive_from, 6);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_Vx,nBuffer_Vector,MPI::DOUBLE,receive_from, 7);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_Vy,nBuffer_Vector,MPI::DOUBLE,receive_from, 8);
			if (nDim == 3) MPI::COMM_WORLD.Recv(&Buffer_Receive_Vz,nBuffer_Vector,MPI::DOUBLE,receive_from, 9);
						
			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				iPoint = geometry->vertex[val_marker][iVertex]->GetNode();	
				for (iVar = 0; iVar < nVar; iVar++) {
					node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVertex][iVar]);
					node[iPoint]->SetGradient_Primitive(iVar, 0, Buffer_Receive_Vx[iVertex][iVar]);
					node[iPoint]->SetGradient_Primitive(iVar, 1, Buffer_Receive_Vy[iVertex][iVar]);
					if (nDim == 3) node[iPoint]->SetGradient_Primitive(iVar, 2, Buffer_Receive_Vz[iVertex][iVar]);
					if (val_mesh == MESH_0) {
						node[iPoint]->SetGradient(iVar, 0, Buffer_Receive_Ux[iVertex][iVar]);
						node[iPoint]->SetGradient(iVar, 1, Buffer_Receive_Uy[iVertex][iVar]);
						if (nDim == 3) node[iPoint]->SetGradient(iVar, 2, Buffer_Receive_Uz[iVertex][iVar]);
						if (config->GetKind_SlopeLimit_Flow() != NONE) node[iPoint]->SetLimiter(iVar, Buffer_Receive_Limit[iVertex][iVar]);
					}
				}
				node[iPoint]->SetLaminarViscosity(Buffer_Receive_LaminarViscosity[iVertex]);
				node[iPoint]->SetEddyViscosity(Buffer_Receive_EddyViscosity[iVertex]);
			}
		}
		
		/*--- Centered scheme ---*/
		if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {

			/*--- Inviscid part ---*/
			double Buffer_Receive_U[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Undivided_Laplacian[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Sensor[geometry->nVertex[val_marker]];
			double Buffer_Receive_Lambda[geometry->nVertex[val_marker]];
			unsigned short Buffer_Receive_Neighbor[geometry->nVertex[val_marker]];
			double Buffer_Receive_BetaInc2[geometry->nVertex[val_marker]];
			double Buffer_Receive_DensityInc[geometry->nVertex[val_marker]];
			double Buffer_Receive_LaminarViscosityInc[geometry->nVertex[val_marker]];

			/*--- Viscous part ---*/
			double Buffer_Receive_LaminarViscosity[geometry->nVertex[val_marker]];
			double Buffer_Receive_EddyViscosity[geometry->nVertex[val_marker]];
			double Buffer_Receive_Vx[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Vy[geometry->nVertex[val_marker]][nVar];
			double Buffer_Receive_Vz[geometry->nVertex[val_marker]][nVar];

			/*--- Dimensionalization ---*/
			unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
			unsigned long nBuffer_Scalar = geometry->nVertex[val_marker];
			int receive_from = abs(SendRecv);
			
			MPI::COMM_WORLD.Recv(&Buffer_Receive_U,nBuffer_Vector,MPI::DOUBLE,receive_from, 0);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Recv(&Buffer_Receive_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Recv(&Buffer_Receive_Sensor,nBuffer_Scalar,MPI::DOUBLE,receive_from, 2);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_Lambda,nBuffer_Scalar,MPI::DOUBLE,receive_from, 3);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_Neighbor,nBuffer_Scalar,MPI::UNSIGNED_SHORT,receive_from, 4);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_LaminarViscosity,nBuffer_Scalar,MPI::DOUBLE,receive_from, 5);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_EddyViscosity,nBuffer_Scalar,MPI::DOUBLE,receive_from, 6);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_Vx,nBuffer_Vector,MPI::DOUBLE,receive_from, 7);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_Vy,nBuffer_Vector,MPI::DOUBLE,receive_from, 8);
			if (nDim == 3) MPI::COMM_WORLD.Recv(&Buffer_Receive_Vz,nBuffer_Vector,MPI::DOUBLE,receive_from, 9);
			
			if (incompressible) {
				MPI::COMM_WORLD.Recv(&Buffer_Receive_BetaInc2,nBuffer_Scalar,MPI::DOUBLE,receive_from, 10);
				MPI::COMM_WORLD.Recv(&Buffer_Receive_DensityInc,nBuffer_Scalar,MPI::DOUBLE,receive_from, 11);
				MPI::COMM_WORLD.Recv(&Buffer_Receive_LaminarViscosityInc,nBuffer_Scalar,MPI::DOUBLE,receive_from, 12);
			}

			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
				for (iVar = 0; iVar < nVar; iVar++) {
					node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVertex][iVar]);
					node[iPoint]->SetGradient_Primitive(iVar, 0, Buffer_Receive_Vx[iVertex][iVar]);
					node[iPoint]->SetGradient_Primitive(iVar, 1, Buffer_Receive_Vy[iVertex][iVar]);
					if (nDim == 3) node[iPoint]->SetGradient_Primitive(iVar, 2, Buffer_Receive_Vz[iVertex][iVar]);
					if (val_mesh == MESH_0) node[iPoint]->SetUndivided_Laplacian(iVar, Buffer_Receive_Undivided_Laplacian[iVertex][iVar]);
				}
				if (val_mesh == MESH_0) node[iPoint]->SetSensor(Buffer_Receive_Sensor[iVertex]);
				node[iPoint]->SetLambda(Buffer_Receive_Lambda[iVertex]);
				if (incompressible) {
					node[iPoint]->SetBetaInc2(Buffer_Receive_BetaInc2[iVertex]);
					node[iPoint]->SetDensityInc(Buffer_Receive_DensityInc[iVertex]);
					node[iPoint]->SetLaminarViscosityInc(Buffer_Receive_LaminarViscosityInc[iVertex]);
				}
				geometry->node[iPoint]->SetnPoint(Buffer_Receive_Neighbor[iVertex]);
				node[iPoint]->SetLaminarViscosity(Buffer_Receive_LaminarViscosity[iVertex]);
				node[iPoint]->SetEddyViscosity(Buffer_Receive_EddyViscosity[iVertex]);
			}
		}			
	}
#endif
}
