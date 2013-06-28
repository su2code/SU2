/*!
 * \file solution_direct_reacting_monatomic.cpp
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

CPlasmaMonatomicSolution::CPlasmaMonatomicSolution(void) : CSolution() { }

CPlasmaMonatomicSolution::CPlasmaMonatomicSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned long iPoint;
	unsigned short iVar, iDim, iSpecies,iFluids;
	double Vel2, Ru;
	double M1, M2, M3;
	double R1, R2, R3;

	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	nFluids = config->GetnFluids() ; 
	nSpecies = config->GetnSpecies(); 
	nVar =  nSpecies + nFluids*nDim + nFluids;	// Continuity + nDim*Momentum + Energy;
	node = new CVariable*[geometry->GetnPoint()];

	/*--- Define  some auxiliar vector related with the residual ---*/
	Residual = new double[nVar];	Residual_Max = new double[nVar];
	Residual_i = new double[nVar];	Residual_j = new double[nVar];
	Res_Conv = new double[nVar];	Res_Visc = new double[nVar];

	/*--- Define some auxiliar vector related with the solution ---*/
	Solution = new double[nVar];
	Solution_i = new double[nVar]; Solution_j = new double[nVar];

	/*--- Define some auxiliar vector related with the geometry ---*/
	Vector = new double[nDim];
	Vector_i = new double[nDim]; Vector_j = new double[nDim];

	/*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTRED) {
		p1_Und_Lapl = new double [geometry->GetnPoint()]; p2_Und_Lapl = new double [geometry->GetnPoint()]; }

	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT) {
		/*--- Point to point Jacobians ---*/
		Jacobian_i = new double* [nVar]; Jacobian_j = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new double [nVar]; Jacobian_j[iVar] = new double [nVar]; }
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
		cvector = new double* [nVar]; 
		for (iVar = 0; iVar < nVar; iVar++)
			cvector[iVar] = new double [nDim];
	}

	M1 = AVOGAD_CONSTANT*config->GetParticle_Mass(0);     
	M3 = AVOGAD_CONSTANT*config->GetParticle_Mass(2);   
	M2 = M1-M3;

	Ru = 8314.462175;
	R1 = Ru/M1;
	R2 = Ru/M2;
	R3 = Ru/M3;

	/*--- Flow infinity initialization stuff ---*/
	Density_Inf = new double [nSpecies];
	Density_Inf[0] = 2.1331E-1;
	Density_Inf[1] = 1E-3*Density_Inf[0]*M2/M1;
	Density_Inf[2] = 1E-3*Density_Inf[0]*M3/M1;

	Pressure_Inf = new double [nFluids];
	Pressure_Inf[0] = .1335E5;
	Pressure_Inf[1] = .1335E2;
	Pressure_Inf[2] = .1335E2;

	Velocity_Inf = new double*[nFluids];
	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		Velocity_Inf[iFluids] = new double [nDim];
		Velocity_Inf[iFluids][0] = 4800.0;
		for (iDim = 1; iDim < nDim; iDim ++)
			Velocity_Inf[iFluids][iDim] = 0.0;
	}

	Mach_Inf = 	new double [nFluids];
	Mach_Inf[0] = Velocity_Inf[0][0]/sqrt(fabs(Gamma*Pressure_Inf[0]/Density_Inf[0]));
	Mach_Inf[1] = Velocity_Inf[1][0]/sqrt(fabs(Gamma*Pressure_Inf[1]/Density_Inf[1]));
	Mach_Inf[2] = Velocity_Inf[2][0]/sqrt(fabs(Gamma*Pressure_Inf[2]/Density_Inf[2]));

	Energy_Inf = new double [nFluids];
	for ( iFluids = 0; iFluids < nFluids; iFluids ++) {
		Vel2 = 0; 
		for (iDim = 0; iDim < nDim; iDim++) {
			Vel2 += Velocity_Inf[iFluids][iDim]*Velocity_Inf[iFluids][iDim];
		}
		Energy_Inf[iFluids] = Pressure_Inf[iFluids]/(Density_Inf[iFluids]*Gamma_Minus_One)+0.5*Vel2;	
	}

	/*--- Flow at the inlet initialization stuff ---*/
	Density_Inlet = new double [nSpecies];

	for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
		Density_Inlet[iSpecies] = Density_Inf[iSpecies];

	Velocity_Inlet  = new double*[nFluids];
	Mach_Inlet 		= new double [nFluids];
	Pressure_Inlet  = new double [nFluids];
	Energy_Inlet = new double [nFluids];

	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		Mach_Inlet[iFluids]     = Mach_Inf[iFluids];
		Pressure_Inlet[iFluids] = Pressure_Inf[iFluids];
		Energy_Inlet[iFluids]   = Energy_Inf[iFluids];
		Velocity_Inlet[iFluids] = new double [nDim];
		for (iDim = 0; iDim < nDim; iDim ++)
			Velocity_Inlet[iFluids][iDim] = Velocity_Inf[iFluids][iDim];
	}

	cout << " Pressure_Inlet[0] = " << Pressure_Inlet[0] << endl;
	cout << " Pressure_Inlet[1] = " << Pressure_Inlet[1] << endl;
	cout << " Pressure_Inlet[2] = " << Pressure_Inlet[2] << endl;

	cout << " Density_Inlet[0] = " << Density_Inlet[0] << endl;
	cout << " Density_Inlet[1] = " << Density_Inlet[1] << endl;
	cout << " Density_Inlet[2] = " << Density_Inlet[2] << endl;

	cout << " Energy_Inlet[0] = " << Energy_Inlet[0] << endl;
	cout << " Energy_Inlet[1] = " << Energy_Inlet[1] << endl;
	cout << " Energy_Inlet[2] = " << Energy_Inlet[2] << endl;

	cout << " Mach_Inlet[0] = " << Mach_Inlet[0] << endl;
	cout << " Mach_Inlet[1] = " << Mach_Inlet[1] << endl;
	cout << " Mach_Inlet[2] = " << Mach_Inlet[2] << endl;

	/*--- Flow at the Outlet initialization stuff ---*/
	/*
	 Density_Outlet = new double [nSpecies];
	 Density_Outlet[0] =( (Gamma+1)* pow(Mach_Inf[0],2)) /(Gamma_Minus_One * pow(Mach_Inf[0],2) + 2.0) * Density_Inf[0];
	 Density_Outlet[1] =( (Gamma+1)* pow(Mach_Inf[1],2)) /(Gamma_Minus_One * pow(Mach_Inf[1],2) + 2.0) * Density_Inf[1];
	 Density_Outlet[2] = Density_Outlet[1]*M3/M2; 

	 Mach_Outlet = new double [nFluids];
	 Mach_Outlet[0] = sqrt( (pow(Mach_Inf[0],2)*(Gamma_Minus_One) + 2.0)/ ( 2.0*Gamma*pow(Mach_Inf[0],2) - (Gamma_Minus_One)));
	 Mach_Outlet[1] = sqrt( (pow(Mach_Inf[1],2)*(Gamma_Minus_One) + 2.0)/ ( 2.0*Gamma*pow(Mach_Inf[1],2) - (Gamma_Minus_One)));
	 Mach_Outlet[2] = Mach_Inf[2]; 

	 Pressure_Outlet = new double [nFluids];
	 Pressure_Outlet[0] = ((2*Gamma*pow(Mach_Inf[0],2))/(Gamma+1.0) - (Gamma_Minus_One)/(Gamma+1.0) ) * Pressure_Inf[0];
	 Pressure_Outlet[1] = ((2*Gamma*pow(Mach_Inf[1],2))/(Gamma+1.0) - (Gamma_Minus_One)/(Gamma+1.0) ) * Pressure_Inf[1];
	 Pressure_Outlet[2] = Pressure_Inf[2]*Density_Outlet[2]/Density_Inf[2];

	 Velocity_Outlet = new double*[nFluids];
	 for (iFluids = 0;iFluids < nFluids; iFluids ++) {
	 Velocity_Outlet[iFluids] = new double [nDim];
	 Velocity_Outlet[iFluids][0] = 1* Mach_Outlet[iFluids] * sqrt(Gamma*Pressure_Outlet[iFluids]/Density_Outlet[iFluids]);
	 Velocity_Outlet[iFluids][1] = 0 * Mach_Outlet[iFluids] * sqrt(Gamma*Pressure_Outlet[iFluids]/Density_Outlet[iFluids]);
	 }
	 */
	/*--- Flow at the Outlet initialization stuff ---*/
	Density_Outlet = new double [nSpecies];
	Density_Outlet[0] = 0.8745E0;
	Density_Outlet[1] = 1E-3*Density_Outlet[0];
	Density_Outlet[2] = Density_Outlet[1]*M3/M2;

	Pressure_Outlet = new double [nFluids];
	Pressure_Outlet[0] = 0.6167E7;
	Pressure_Outlet[1] = 0.6149E4;
	Pressure_Outlet[2] = 0.5461E2;

	Velocity_Outlet = new double*[nFluids];
	for (iFluids = 0; iFluids < nFluids; iFluids ++) {
		Velocity_Outlet[iFluids] = new double [nDim];
		Velocity_Outlet[iFluids][0] = 188.0;
		for (iDim = 1; iDim < nDim; iDim ++)
			Velocity_Outlet[iFluids][iDim] = 0.0;
	}

	Energy_Outlet = new double [nFluids];
	for ( iFluids = 0; iFluids < nFluids; iFluids ++) {
		Vel2 = 0; 
		for (iDim = 0; iDim < nDim; iDim++) {
			Vel2 += Velocity_Outlet[iFluids][iDim]*Velocity_Outlet[iFluids][iDim];
		}
		Energy_Outlet[iFluids] = Pressure_Outlet[iFluids]/(Density_Outlet[iFluids]*Gamma_Minus_One)+0.5*Vel2;
	}

	Mach_Outlet = new double [nFluids];
	Mach_Outlet[0] = Velocity_Outlet[0][0]/sqrt(Gamma*Pressure_Outlet[0]/Density_Outlet[0]);
	Mach_Outlet[1] = Velocity_Outlet[1][0]/sqrt(Gamma*Pressure_Outlet[1]/Density_Outlet[1]);
	Mach_Outlet[2] = Velocity_Outlet[2][0]/sqrt(Gamma*Pressure_Outlet[2]/Density_Outlet[2]);

	cout << " Pressure_Outlet[0] = " << Pressure_Outlet[0] << endl;
	cout << " Pressure_Outlet[1] = " << Pressure_Outlet[1] << endl;
	cout << " Pressure_Outlet[2] = " << Pressure_Outlet[2] << endl;

	cout << " Density_Outlet[0] = " << Density_Outlet[0] << endl;
	cout << " Density_Outlet[1] = " << Density_Outlet[1] << endl;
	cout << " Density_Outlet[2] = " << Density_Outlet[2] << endl;

	cout << " Energy_Outlet[0] = " << Energy_Outlet[0] << endl;
	cout << " Energy_Outlet[1] = " << Energy_Outlet[1] << endl;
	cout << " Energy_Outlet[2] = " << Energy_Outlet[2] << endl;

	cout << " Velocity_Outlet[0] = " << Velocity_Outlet[0][0] << endl;
	cout << " Velocity_Outlet[1] = " << Velocity_Outlet[1][0] << endl;
	cout << " Velocity_Outlet[2] = " << Velocity_Outlet[2][0] << endl;

	cout << " Mach_Outlet[0] = " << Mach_Outlet[0] << endl;
	cout << " Mach_Outlet[1] = " << Mach_Outlet[1] << endl;
	cout << " Mach_Outlet[2] = " << Mach_Outlet[2] << endl;

	cout << " Temperature_Outlet[0] = " << Pressure_Outlet[0]/(R1*Density_Outlet[0]) << endl;
	cout << " Temperature_Outlet[1] = " << Pressure_Outlet[1]/(R2*Density_Outlet[1]) << endl;
	cout << " Temperature_Outlet[2] = " << Pressure_Outlet[2]/(R3*Density_Outlet[2]) << endl;

	/*--- Restart the solution from file information ---*/
	if (!restart) {
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

			if(geometry->node[iPoint]->GetCoord(0) <= 0.0) {
				node[iPoint] = new CPlasmaMonatomicVariable(Density_Inlet, Velocity_Inlet, Energy_Inlet,nDim, nVar, nSpecies, nFluids, config);
			}
			else {
				node[iPoint] = new CPlasmaMonatomicVariable(Density_Outlet, Velocity_Outlet, Energy_Outlet, nDim, nVar, nSpecies, nFluids, config);
			}
		}
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

		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			point_line >> index;
			for (iVar = 0; iVar < nVar; iVar ++) {
				point_line >> Solution[iVar];
			}
			node[iPoint] = new CPlasmaMonatomicVariable(Solution, nDim, nVar, config);	
		}
		restart_file.close();
	}

}

CPlasmaMonatomicSolution::~CPlasmaMonatomicSolution(void) {
	unsigned short iVar, iDim;

	delete [] Residual;		delete [] Residual_Max;
	delete [] Residual_i;	delete [] Residual_j;
	delete [] Res_Conv;		delete [] Res_Visc;
	delete [] Solution;		delete [] Solution_i; delete [] Solution_j;
	delete [] Vector_i;		delete [] Vector_j;
	delete [] p1_Und_Lapl;  delete [] p2_Und_Lapl;
	delete [] xsol;			delete [] rhs;

	delete [] Energy_Inf;	delete [] Energy_Inlet;		delete [] Energy_Outlet;	
	delete [] Density_Inf;	delete [] Density_Inlet;	delete [] Density_Outlet;
	delete [] Pressure_Inf;	delete [] Pressure_Inlet;	delete [] Pressure_Outlet;
	delete [] Mach_Inf;		delete [] Mach_Inlet;		delete [] Mach_Outlet;	

	for (iVar = 0; iVar < nFluids; iVar ++) {
		delete [] Velocity_Inf[iVar];
		delete [] Velocity_Inlet[iVar];
		delete [] Velocity_Outlet[iVar];
	}
	delete[] Velocity_Inf; 	delete [] Velocity_Inlet; delete [] Velocity_Outlet;

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

}

void CPlasmaMonatomicSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		/*--- Compute squared velocity, sound velocity, pressure, and enthalpy ---*/
		node[iPoint]->SetVelocity2();
		node[iPoint]->SetSoundSpeed(Gamma);
		node[iPoint]->SetPressure(Gamma);
		node[iPoint]->SetEnthalpy();

		/*--- Initialize the convective and viscous residual vector ---*/
		node[iPoint]->Set_ResConv_Zero();
		if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit) 
			node[iPoint]->Set_ResVisc_Zero();
	}
	SetSolution_Gradient_GG(geometry); 

	/*--- Inicialize the jacobian matrices ---*/
	if (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT)
		Jacobian.SetValZero();
}

void CPlasmaMonatomicSolution::SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {
	double *Normal, Area, dV, Mean_SoundSpeed, Mean_ProjVel, Lambda, Local_Delta_Time, Global_Delta_Time = 1E6;
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker, iFluids;
	double u, c, dx;

	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		node[iPoint]->SetMax_Lambda_Inv(0.0);
		node[iPoint]->SetVelocity2();
		node[iPoint]->SetSoundSpeed(Gamma);
	}

	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Point identification, Normal vector and area ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); 
		jPoint = geometry->edge[iEdge]->GetNode(1);
		Normal = geometry->edge[iEdge]->GetNormal();
		Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

		/*--- Mean Values ---*/
		//Lambda = 0.0;
		for ( iFluids = 0; iFluids < nFluids; iFluids ++ ) {
			Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal,iFluids) + node[jPoint]->GetProjVel(Normal,iFluids));
			Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed(iFluids) + node[jPoint]->GetSoundSpeed(iFluids));
			Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed * Area;
			/*--- Inviscid contribution ---*/
			node[iPoint]->AddMax_Lambda_Inv(Lambda,iFluids);
			node[jPoint]->AddMax_Lambda_Inv(Lambda,iFluids);
		}
	}

	/*--- Loop boundary edges ---*/
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) { 
		for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

			/*--- Point identification, Normal vector and area ---*/
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

			//Lambda = 0.0;
			for ( iFluids = 0; iFluids < nFluids; iFluids ++ ) {	
				/*--- Mean Values ---*/
				Mean_ProjVel = node[iPoint]->GetProjVel(Normal,iFluids);
				Mean_SoundSpeed = node[iPoint]->GetSoundSpeed(iFluids);	
				/*--- Inviscid contribution ---*/
				Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed * Area;
				node[iPoint]->AddMax_Lambda_Inv(Lambda, iFluids);
			}
		}
	}
	u = 4800;
	c = 87110;
	c = 732;
	dx = 0.04/810;
	Local_Delta_Time = config->GetCFL(iMesh)*dx/(u+c);
	//cout << " dt = " << Local_Delta_Time << endl;

	/*--- Each element uses their own speed ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		dV = geometry->node[iPoint]->GetVolume();
		/*Lambda_iFluid = max(node[iPoint]->GetMax_Lambda_Inv(0),node[iPoint]->GetMax_Lambda_Inv(1));
		 Lambda_iFluid = max(node[iPoint]->GetMax_Lambda_Inv(2), Lambda_iFluid);
		 Local_Delta_Time = config->GetCFL(iMesh)*dV / Lambda_iFluid;
		 */
		if (config->GetUnsteady_Simulation() == TIME_STEPPING) Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
		else node[iPoint]->SetDelta_Time(Local_Delta_Time);
	}

	/*--- For exact time solution use the minimum delta time of the whole mesh ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
		for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
			node[iPoint]->SetDelta_Time(Global_Delta_Time);
	}
}

void CPlasmaMonatomicSolution::Centred_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
		CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

void CPlasmaMonatomicSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
		CConfig *config, unsigned short iMesh) {

	double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, *Limiter_i = NULL, *Limiter_j = NULL, *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	bool high_order_diss = ((config->GetKind_Upwind() == ROE_2ND) && (iMesh == MESH_0));
	SetSolution_Gradient_GG(geometry); 

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
		U_i = node[iPoint]->GetSolution();
		U_j = node[jPoint]->GetSolution();

		solver->SetConservative(U_i, U_j);

	//	if (geometry->node[iPoint]->GetCoord(0)/.002 > .90) high_order_diss = true;
	//	else high_order_diss = false;

		if (high_order_diss) {
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

		for (iVar = 0; iVar < nVar; iVar ++)
			if (Res_Conv[iVar] != Res_Conv[iVar] )
				cout << "Res_Conv = " << Res_Conv[iVar] << " iPoint = " << iPoint  << " x = " <<  geometry->node[iPoint]->GetCoord(0)/.002<< endl;

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


void CPlasmaMonatomicSolution::SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh) {
	unsigned short iVar;
	unsigned long iPoint;
	double **Gradient;
	double dtime;
	double dx =0.004/81.0;

	dtime = config->GetCFL(0) * dx/(4800+732.0);

	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	for (iVar = 0; iVar < nVar; iVar++)
		Residual[iVar] = 0;

	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 

		Gradient = solution_container[ELEC_SOL]->node[iPoint]->GetGradient(); 

		/*--- Set gradient of phi from electrostatic potential equation ---*/
		solver->SetConsVarGradient(Gradient);

		/*--- Set solution  ---*/
		solver->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

		/*--- Set control volume ---*/
		solver->SetVolume(geometry->node[iPoint]->GetVolume());

		/*--- Compute Residual ---*/
		solver->SetResidual(Residual,Jacobian_i, config);
		solution_container[PLASMA_SOL]->node[iPoint]->SetSource(Residual);

		/*--- Add Residual ---*/
		node[iPoint]->SubtractRes_Conv(Residual);

		if (implicit) {
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
			solution_container[PLASMA_SOL]->node[iPoint]->SetSourceJacobian(Jacobian_i); 

		}
	}
}

void CPlasmaMonatomicSolution::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint, iEdge;
	double Pressure_i = 0, Pressure_j = 0;
	unsigned short iVar;
	bool boundary_i, boundary_j;
	double *Diff = new double[nVar];

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetUnd_LaplZero();

	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		Pressure_i = node[iPoint]->GetPressure();
		Pressure_j = node[jPoint]->GetPressure();

		for (iVar = 0; iVar < nVar; iVar++)
			Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);		
		Diff[nVar-1] = (node[iPoint]->GetSolution(nVar-1)+Pressure_i) - (node[jPoint]->GetSolution(nVar-1)+Pressure_j);

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

	delete [] Diff;
}

void CPlasmaMonatomicSolution::SetSpectral_Radius(CGeometry *geometry, CConfig *config) {
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker;
	double *Normal, Area, ProjVel, Lambda, SoundSpeed_i, SoundSpeed_j, SoundSpeed;

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetLambda(0.0);

	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		SoundSpeed_i = node[iPoint]->GetSoundSpeed();
		SoundSpeed_j = node[jPoint]->GetSoundSpeed();

		Normal = geometry->edge[iEdge]->GetNormal();
		Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

		/*--- Inviscid contribution to the Point i ---*/
		ProjVel = node[iPoint]->GetProjVel(Normal);
		Lambda = 0.5*(fabs(ProjVel) + SoundSpeed_i*Area);
		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);

		/*--- Inviscid contribution to the Point j ---*/
		ProjVel = node[jPoint]->GetProjVel(Normal);
		Lambda = 0.5*(fabs(ProjVel) + SoundSpeed_j*Area);
		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);
	}

	/*--- Loop boundary edges ---*/
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			if (geometry->node[iPoint]->GetDomain()) {
				SoundSpeed = node[iPoint]->GetSoundSpeed();
				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
				/*--- Inviscid contribution to the iPoint ---*/
				ProjVel = node[iPoint]->GetProjVel(Normal);
				Lambda = (fabs(ProjVel) + SoundSpeed*Area);
				node[iPoint]->AddLambda(Lambda);
			}
		}
}

void CPlasmaMonatomicSolution::SetPress_Switch(CGeometry *geometry) {
	unsigned long iEdge, iPoint, jPoint;
	double Pressure_i, Pressure_j;

	/*--- Reset variables ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		p1_Und_Lapl[iPoint] = 0.0;
		p2_Und_Lapl[iPoint] = 0.0;
	}

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
		node[iPoint]->SetSensor(fabs(p1_Und_Lapl[iPoint])/p2_Und_Lapl[iPoint]);
}

void CPlasmaMonatomicSolution::RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, 
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

void CPlasmaMonatomicSolution::ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
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
			node[iPoint]->AddSolution(iVar, -(Residual[iVar])*Delta);
			AddRes_Max( iVar, Residual[iVar]*Residual[iVar]*Vol );
		}
	}
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, sqrt(GetRes_Max(iVar)) );
}

void CPlasmaMonatomicSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
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
	Jacobian.SGSSolution(rhs, xsol, 1e-9, 1000, true);
	//  Jacobian.LU_SGSIteration(rhs, xsol);

	unsigned short interior = 0;
	if (nDim ==2 ) {
		interior = 10;
		/*--- Update solution (system written in terms of increments) ---*/

		for (iPoint = interior; iPoint < geometry->GetnPointDomain(); iPoint++)
			for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);


		for (iPoint = 0; iPoint < interior; iPoint++)
			for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->AddSolution(iVar,xsol[(interior+iPoint)*nVar+iVar]);


		for (iPoint = 0; iPoint < interior; iPoint++)
			for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->SetSolution(iVar,node[iPoint+interior]->GetSolution(iVar));

	}

	if (nDim ==3 ) {
		unsigned short Start = 400;
		unsigned short iJump  = 1;
		unsigned short iStart;
		/*--- Update solution (system written in terms of increments) ---*/

		for (iPoint = interior; iPoint < geometry->GetnPointDomain(); iPoint++)
			for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);


		for (iPoint = 0; iPoint < interior; iPoint++)
			for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->AddSolution(iVar,xsol[(interior+iPoint)*nVar+iVar]);


		for (iJump = 0; iJump <5; iJump ++) {
			iStart = Start + iJump * 405;
			for (iPoint = iStart; iPoint < iStart + 5; iPoint++)
				for (iVar = 0; iVar < nVar; iVar++)
					node[iPoint]->SetSolution(iVar,node[iPoint-5]->GetSolution(iVar));
		}

	}
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, sqrt(GetRes_Max(iVar)));
}

void CPlasmaMonatomicSolution::SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint, iEdge, iVertex;
	unsigned short iDim, iVar, iMarker;
	double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average, 
	Partial_Gradient, Partial_Res, *Normal;

	PrimVar_Vertex = new double [nVar+1];
	PrimVar_i = new double [nVar+1];
	PrimVar_j = new double [nVar+1]; 

	/*--- Set Gradient_Primitive to Zero ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetGradient_PrimitiveZero();

	/*--- Loop interior edges ---*/ 
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

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
				PrimVar_Vertex[0] = node[iPoint]->GetTemperature();
				for (iDim = 0; iDim < nDim; iDim++)
					PrimVar_Vertex[iDim+1] = node[iPoint]->GetVelocity(iDim);
				PrimVar_Vertex[nDim+1] = node[iPoint]->GetPressure();
				PrimVar_Vertex[nDim+2] = node[iPoint]->GetDensity();

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
				//				if (Partial_Gradient > config->GetPrimGrad_Threshold()) Partial_Gradient = config->GetPrimGrad_Threshold();
				//				if (Partial_Gradient < -config->GetPrimGrad_Threshold()) Partial_Gradient = -config->GetPrimGrad_Threshold();
				node[iPoint]->SetGradient_Primitive(iVar, iDim, Partial_Gradient);
			}
		}		
	}

	delete [] PrimVar_Vertex;
	delete [] PrimVar_i;
	delete [] PrimVar_j;
}

void CPlasmaMonatomicSolution::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config) {
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
				node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
			}
		}
	}

	delete [] PrimVar_i;
	delete [] PrimVar_j;
}

void CPlasmaMonatomicSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {

}

void CPlasmaMonatomicSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {

}

void CPlasmaMonatomicSolution::BC_Electrode(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {

}

void CPlasmaMonatomicSolution::BC_Dielectric(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {

}


void CPlasmaMonatomicSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iFluids, loc = 0;
	unsigned short iVar, jVar, iDim;
	double dS, sq_vel,  **P_Matrix, **invP_Matrix, *Face_Normal, 
	*kappa, *U_domain, *U_inlet, *U_update, *W_domain, *W_inlet, *W_update;
	double *vn = NULL, *rho = NULL, *rhoE = NULL, **velocity = NULL, *c = NULL, *pressure = NULL, *enthalpy = NULL;

	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	kappa = new double[nDim]; 
	U_domain = new double[nVar]; U_inlet = new double[nVar]; U_update = new double[nVar];
	W_domain = new double[nVar]; W_inlet = new double[nVar]; W_update = new double[nVar];
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

			/*--- Solution at the inlet ---*/
			for (iFluids = 0; iFluids < nFluids; iFluids ++) {
				loc = iFluids * (nDim+2);
				U_inlet[loc + 0] = GetDensity_Inlet(iFluids);
				U_inlet[loc + 1] = GetDensity_Velocity_Inlet(0,iFluids);
				U_inlet[loc + 2] = GetDensity_Velocity_Inlet(1,iFluids);
				U_inlet[loc + 3] = GetDensity_Energy_Inlet(iFluids);
				if (nDim == 3) {
					U_inlet[loc + 3] = GetDensity_Velocity_Inlet(2,iFluids);
					U_inlet[loc + 4] = GetDensity_Energy_Inlet(iFluids);
				}
			}

			Face_Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			dS = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				dS += Face_Normal[iDim]*Face_Normal[iDim];
			dS = sqrt (dS);

			for (iDim = 0; iDim < nDim; iDim++) {
				kappa[iDim] = -Face_Normal[iDim]/dS;
			}
			/*--- Computation of P and inverse P matrix using values at the infinity ---*/
			vn		 = new double [nFluids];
			rho		 = new double [nSpecies];
			rhoE	 = new double [nFluids];
			c		 = new double [nFluids];
			velocity = new double*[nFluids];
			pressure = new double [nFluids];
			enthalpy = new double [nFluids];

			for (iFluids = 0; iFluids < nFluids; iFluids ++) {
				sq_vel = 0.0;
				vn[iFluids] = 0.0;
				velocity[iFluids] = new double[nDim];
				loc = iFluids * (nDim+2);

				rho[iFluids] = (U_inlet[loc + 0] + U_domain[loc + 0])/2.0;
				rhoE[iFluids] = (U_inlet[loc + nDim + 1] + U_domain[loc + nDim + 1])/2.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					velocity[iFluids][iDim] = (U_inlet[loc + iDim+1]/U_inlet[loc + 0] + U_domain[loc + iDim+1]/U_domain[loc + 0] )/2.0;
					sq_vel += velocity[iFluids][iDim]*velocity[iFluids][iDim];
					vn[iFluids] += velocity[iFluids][iDim]*kappa[iDim]*dS;
				}
				c[iFluids] = sqrt(fabs(Gamma*Gamma_Minus_One*(rhoE[iFluids]/rho[iFluids]-0.5*sq_vel)));
			}

			solver->GetPMatrix_inv(rho, velocity, c, kappa, invP_Matrix);
			solver->GetPMatrix(rho, velocity, c, kappa, P_Matrix);

			/*--- computation of characteristics variables at the infinity ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				W_inlet[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_inlet[iVar] +=invP_Matrix[iVar][jVar]*U_inlet[jVar];
			}

			/*		for (iFluids = 0; iFluids < nFluids; iFluids ++) {
				sq_vel = 0.0;
				vn[iFluids] = 0.0;
				velocity[iFluids] = new double[nDim];
				loc = iFluids * (nDim+2);
				rho[iFluids] = U_domain[loc + 0];
				rhoE[iFluids] = U_domain[loc + nDim + 1];
				for (iDim = 0; iDim < nDim; iDim++) {
					velocity[iFluids][iDim] = U_domain[loc + iDim+1]/rho[iFluids];
					sq_vel += velocity[iFluids][iDim]*velocity[iFluids][iDim];
					vn[iFluids] += velocity[iFluids][iDim]*kappa[iDim]*dS;
				}
				c[iFluids] = sqrt(fabs(Gamma*Gamma_Minus_One*(rhoE[iFluids]/rho[iFluids]-0.5*sq_vel)));
			}

			solver->GetPMatrix_inv(rho, velocity, c, kappa, invP_Matrix);
			solver->GetPMatrix(rho, velocity, c, kappa, P_Matrix);
			 */
			/*--- computation of characteristics variables at the wall ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				W_domain[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_domain[iVar] +=invP_Matrix[iVar][jVar]*U_domain[jVar];
			}

			/*--- fix characteristics value ---*/
			for (iFluids = 0; iFluids < nFluids; iFluids ++) {
				loc = iFluids * (nDim+2);

				if (nDim == 2) {
					if(vn[iFluids] > 0.0) {
						W_update[loc + 0] = W_inlet[loc + 0];
						W_update[loc + 1] = W_inlet[loc + 1];
					}
					else {
						W_update[loc + 0] = W_domain[loc + 0];
						W_update[loc + 1] = W_domain[loc + 1];
					}

					if(vn[iFluids]+c[iFluids]*dS > 0.0) W_update[loc + 2] = W_inlet[loc + 2];
					else W_update[loc + 2] = W_domain[loc + 2];

					if(vn[iFluids]-c[iFluids]*dS > 0.0) W_update[loc + 3] = W_inlet[loc + 3];
					else W_update[loc + 3] = W_domain[loc + 3];
				}

				if (nDim == 3) {
					if(vn[iFluids] > 0.0) {
						W_update[loc + 0] = W_domain[loc + 0];
						W_update[loc + 1] = W_domain[loc + 1];
						W_update[loc + 2] = W_domain[loc + 2];
					}
					else {
						W_update[loc + 0] = W_inlet[loc + 0];
						W_update[loc + 1] = W_inlet[loc + 1];
						W_update[loc + 2] = W_inlet[loc + 2];
					}

					if(vn[iFluids] + c[iFluids]*dS > 0.0) W_update[loc + 3] = W_domain[loc + 3];
					else W_update[loc + 3] = W_inlet[loc + 3];

					if(vn[iFluids]-c[iFluids]*dS > 0.0) W_update[loc + 4] = W_domain[loc + 4];
					else W_update[loc + 4] = W_inlet[loc + 4];
				}
				/*			cout << " iFluids = " << 	iFluids << endl;
				cout << " vn 	  = " << 	vn[iFluids] << endl;
				cout << " vn + c  = " <<	vn[iFluids]+c[iFluids]*dS<< endl;
				cout << " vn - c  = " <<	vn[iFluids]-c[iFluids]*dS<< endl;
				cout << endl;
				cout << endl;
				 */
			}

			/*--- conservative variables using characteristics ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				U_update[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					U_update[iVar] +=P_Matrix[iVar][jVar]*W_update[jVar];
			}

			/*--- Residual computation ---*/
			/*			for (iFluids = 0; iFluids < nFluids; iFluids ++) {
				loc = iFluids * (nDim+2);
				rho[iFluids] = U_update[loc + 0];
				rhoE[iFluids] = U_update[loc + nDim + 1];
				sq_vel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					velocity[iFluids][iDim]  = U_update[loc + iDim+1]/rho[iFluids];
					sq_vel +=velocity[iFluids][iDim]*velocity[iFluids][iDim];
				}
				c[iFluids] = sqrt(fabs(Gamma*Gamma_Minus_One*(rhoE[iFluids]/rho[iFluids] - 0.5*sq_vel)));
				pressure[iFluids] = (c[iFluids] * c[iFluids] * rho[iFluids]) / Gamma;
				enthalpy[iFluids] = (rhoE[iFluids] + pressure[iFluids]) / rho[iFluids];
			}
			 */
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			solver->SetNormal(Vector);

			/*--- Compute the residual using an upwind scheme ---*/
			solver->SetConservative(U_domain, U_update);
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddRes_Conv(Residual);

			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

		}
	}

	for ( iFluids = 0; iFluids < nFluids; iFluids ++ ) {
		delete [] velocity[iFluids];
	}
	delete [] velocity;
	delete [] vn; delete [] rho; delete [] pressure;
	delete [] rhoE; delete [] c; delete [] enthalpy;
	delete [] kappa;
	delete [] U_domain; delete [] U_inlet; delete [] U_update;
	delete [] W_domain; delete [] W_inlet; delete [] W_update;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Matrix[iVar];
		delete [] invP_Matrix[iVar];
	}
	delete [] P_Matrix;
	delete [] invP_Matrix;
}

void CPlasmaMonatomicSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {

	unsigned long iVertex, iPoint;
	unsigned short iFluids, loc = 0;
	unsigned short iVar, jVar, iDim;
	double Area, sq_vel,  **P_Matrix, **invP_Matrix, *Normal,
	*UnitaryNormal, *U_domain, *U_outlet, *U_update, *W_domain, *W_outlet, *W_update;

	double **P_Matrix_domain, **invP_Matrix_domain;
	double *vn = NULL, *rho = NULL, *rhoE = NULL, **velocity = NULL, *c = NULL, *pressure = NULL, *enthalpy = NULL;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	UnitaryNormal = new double[nDim];
	U_domain = new double[nVar]; U_outlet = new double[nVar]; U_update = new double[nVar];
	W_domain = new double[nVar]; W_outlet = new double[nVar]; W_update = new double[nVar];
	P_Matrix = new double* [nVar]; invP_Matrix = new double* [nVar];
	P_Matrix_domain = new double* [nVar]; invP_Matrix_domain = new double* [nVar];

	for (iVar = 0; iVar < nVar; iVar++) {
		P_Matrix[iVar] = new double [nVar];
		invP_Matrix[iVar] = new double [nVar];
		P_Matrix_domain[iVar] = new double [nVar];
		invP_Matrix_domain[iVar] = new double [nVar];
	}

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Solution at the Outlet ---*/
			for (iFluids = 0; iFluids < nFluids; iFluids ++) {
				loc = iFluids * (nDim+2);

				U_outlet[loc + 0] = GetDensity_Outlet(iFluids);
				U_outlet[loc + 1] = GetDensity_Velocity_Outlet(0,iFluids);
				U_outlet[loc + 2] = GetDensity_Velocity_Outlet(1,iFluids);
				U_outlet[loc + 3] = GetDensity_Energy_Outlet(iFluids);
				if (nDim == 3) {
					U_outlet[loc + 3] = GetDensity_Velocity_Outlet(2,iFluids);
					U_outlet[loc + 4] = GetDensity_Energy_Outlet(iFluids);
				}
			}

			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

			Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;

			vn		 = new double [nFluids];
			rho		 = new double [nSpecies];
			rhoE	 = new double [nFluids];
			c		 = new double [nFluids];
			velocity = new double*[nFluids];
			pressure = new double [nFluids];
			enthalpy = new double [nFluids];

			/*--- Computation of P and inverse P matrix using values at the domain ---*/
			for (iFluids = 0; iFluids < nFluids; iFluids ++) {
				sq_vel = 0.0;
				vn[iFluids] = 0.0;
				velocity[iFluids] = new double[nDim];
				loc = iFluids * (nDim + 2);
				rho[iFluids]  = (U_domain[loc + 0] + U_outlet[loc+0]) / 2.0;
				rhoE[iFluids] = (U_domain[loc + nDim + 1] + U_outlet[loc + nDim + 1])/2.0;

				for (iDim = 0; iDim < nDim; iDim++) {
					velocity[iFluids][iDim] = (U_domain[loc + iDim+1]/U_domain[loc + 0] + U_outlet[loc + iDim+1]/U_outlet[loc + 0] );
					sq_vel += velocity[iFluids][iDim]*velocity[iFluids][iDim];
					vn[iFluids] += velocity[iFluids][iDim]*UnitaryNormal[iDim]*Area;
				}
				c[iFluids] = sqrt(fabs(Gamma*Gamma_Minus_One*(rhoE[iFluids]/rho[iFluids]-0.5*sq_vel)));
			}
			solver->GetPMatrix_inv(rho, velocity, c, UnitaryNormal, invP_Matrix);
			solver->GetPMatrix(rho, velocity, c, UnitaryNormal, P_Matrix);

			/*--- computation of characteristics variables at the wall ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				W_domain[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_domain[iVar] +=invP_Matrix[iVar][jVar]*U_domain[jVar];
			}

			/*--- computation of characteristics variables at the infinity ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				W_outlet[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_outlet[iVar] +=invP_Matrix[iVar][jVar]*U_outlet[jVar];
			}

			/*--- fix characteristics value ---*/
			for (iFluids = 0; iFluids < nFluids; iFluids ++) {
				loc = iFluids * (nDim +2);

				if (nDim == 2) {
					if(vn[iFluids] > 0.0) {
						W_update[loc + 0] = W_domain[loc + 0];
						W_update[loc + 1] = W_domain[loc + 1];
					}
					else {
						W_update[loc + 0] = W_outlet[loc + 0];
						W_update[loc + 1] = W_outlet[loc + 1];
					}

					if(vn[iFluids]+c[iFluids]*Area > 0.0) W_update[loc + 2] = W_domain[loc + 2];
					else W_update[loc + 2] = W_outlet[loc + 2];

					if(vn[iFluids]-c[iFluids]*Area > 0.0) W_update[loc + 3] = W_domain[loc + 3];
					else W_update[loc + 3] = W_outlet[loc + 3];

				}

				if (nDim == 3) {
					if(vn[iFluids] > 0.0) {
						W_update[loc + 0] = W_domain[loc + 0];
						W_update[loc + 1] = W_domain[loc + 1];
						W_update[loc + 2] = W_domain[loc + 2];
					}
					else {
						W_update[loc + 0] = W_outlet[loc + 0];
						W_update[loc + 1] = W_outlet[loc + 1];
						W_update[loc + 2] = W_outlet[loc + 2];
					}

					if(vn[iFluids]+c[iFluids]*Area > 0.0) W_update[loc + 3] = W_domain[loc + 3];
					else W_update[loc + 3] = W_outlet[loc + 3];


					if(vn[iFluids]-c[iFluids]*Area > 0.0) W_update[loc + 4] = W_domain[loc + 4];
					else W_update[loc + 4] = W_outlet[loc + 4];
				}
			}

			/*--- conservative variables using characteristics ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				U_update[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					U_update[iVar] +=P_Matrix[iVar][jVar]*W_update[jVar];
			}

			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			solver->SetNormal(Vector);

			/*--- Compute the residual using an upwind scheme ---*/
			solver->SetConservative(U_domain, U_domain);
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddRes_Conv(Residual);


			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

			if (nDim == 2 ) {
				for (short index=0; index<5; index++) {
					node[index]->Set_ResConv_Zero();
					node[index]->Set_ResVisc_Zero();
				}
			}
		}
	}

	for ( iFluids = 0; iFluids < nFluids; iFluids ++ ) {
		delete [] velocity[iFluids];
	}
	delete [] velocity;
	delete [] vn; delete [] rho; delete [] pressure;
	delete [] rhoE; delete [] c; delete [] enthalpy;
	delete [] UnitaryNormal;
	delete [] U_domain; delete [] U_outlet; delete [] U_update;
	delete [] W_domain; delete [] W_outlet; delete [] W_update;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Matrix[iVar];
		delete [] invP_Matrix[iVar];
	}
	delete [] P_Matrix;
	delete [] invP_Matrix;
}

void CPlasmaMonatomicSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iVar, jVar, iDim;
	unsigned short loc = 0, iFluids;
	double Area, sq_vel, *vn, *rho, *rhoE, *c, **P_Matrix, **invP_Matrix, *pressure, *enthalpy, *Normal,
	*UnitaryNormal, **velocity, *U_wall, *U_infty, *U_update, *W_wall, *W_infty, *W_update;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	UnitaryNormal = new double[nDim];
	U_wall = new double[nVar]; U_infty = new double[nVar]; U_update = new double[nVar];
	W_wall = new double[nVar]; W_infty = new double[nVar]; W_update = new double[nVar];
	P_Matrix = new double* [nVar]; invP_Matrix = new double* [nVar];

	rho = new double [nSpecies];
	rhoE = new double [nFluids];
	pressure  = new double [nFluids];
	enthalpy  = new double [nFluids];
	c    = new double [nFluids];
	vn = new double [nFluids];
	velocity = new double * [nFluids];
	for (iFluids =0; iFluids < nFluids; iFluids ++)
		velocity[iFluids] = new double [nDim];

	for (iVar = 0; iVar < nVar; iVar++) {
		P_Matrix[iVar] = new double [nVar];
		invP_Matrix[iVar] = new double [nVar];
	}

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_wall[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Solution at the infinity ---*/
			for(iFluids = 0; iFluids < nFluids; iFluids ++) {
				loc = iFluids * (nDim +2);
				U_infty[loc + 0] = GetDensity_Inf(iFluids);
				U_infty[loc + 1] = GetDensity_Velocity_Inf(0,iFluids);
				U_infty[loc + 2] = GetDensity_Velocity_Inf(1,iFluids);
				U_infty[loc + 3] = GetDensity_Energy_Inf(iFluids);

				if (nDim == 3) {
					U_infty[loc + 3] = GetDensity_Velocity_Inf(2, iFluids);
					U_infty[loc + 4] = GetDensity_Energy_Inf(iFluids);
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
			sq_vel = 0.0;
			for (iFluids = 0; iFluids < nFluids; iFluids ++){
				loc = iFluids*(nDim+2);
				sq_vel = 0.0;
				vn[iFluids] = 0.0;
				rho[iFluids] = U_infty[loc + 0];
				rhoE[iFluids] = U_infty[loc + nDim +1];
				for (iDim = 0; iDim < nDim; iDim++) {
					velocity[iFluids][iDim] = U_infty[loc + iDim+1]/rho[iFluids];
					sq_vel += velocity[iFluids][iDim]*velocity[iFluids][iDim];
					vn[iFluids] += velocity[iFluids][iDim]*UnitaryNormal[iDim]*Area;
				}
				c[iFluids] = sqrt(fabs(Gamma*Gamma_Minus_One*(rhoE[iFluids]/rho[iFluids]-0.5*sq_vel)));
			}

			solver->GetPMatrix_inv(rho, velocity, c, UnitaryNormal, invP_Matrix);
			solver->GetPMatrix(rho, velocity, c, UnitaryNormal, P_Matrix);

			/*--- computation of characteristics variables at the wall ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				W_wall[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_wall[iVar] +=invP_Matrix[iVar][jVar]*U_wall[jVar];
			}

			/*--- computation of characteristics variables at the infinity ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				W_infty[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					W_infty[iVar] +=invP_Matrix[iVar][jVar]*U_infty[jVar];
			}

			/*--- fix characteristics value ---*/
			if (nDim == 2) {
				for (iFluids = 0; iFluids < nFluids; iFluids ++ ) {
					loc = iFluids*(nDim+2);
					if(vn[iFluids] > 0.0) {
						W_update[loc + 0] = W_wall[loc + 0];
						W_update[loc + 1] = W_wall[loc + 1];
					}
					else {
						W_update[loc + 0] = W_infty[loc + 0];
						W_update[loc + 1] = W_infty[loc + 1];
					}

					if(vn[iFluids]+c[iFluids]*Area > 0.0) W_update[loc + 2] = W_wall[loc + 2];
					else W_update[loc + 2] = W_infty[loc + 2];

					if(vn[iFluids]-c[iFluids]*Area > 0.0) W_update[loc + 3] = W_wall[loc + 3];
					else W_update[loc + 3] = W_infty[loc + 3];
				}
			}

			if (nDim == 3) {
				for (iFluids = 0; iFluids < nFluids; iFluids ++ ) {
					loc = iFluids*(nDim+2);
					if(vn[iFluids] > 0.0) {
						W_update[loc + 0] = W_wall[loc + 0];
						W_update[loc + 1] = W_wall[loc + 1];
						W_update[loc + 2] = W_wall[loc + 2];
					}
					else {
						W_update[loc + 0] = W_infty[loc + 0];
						W_update[loc + 1] = W_infty[loc + 1];
						W_update[loc + 2] = W_infty[loc + 2];
					}

					if(vn[iFluids]+c[iFluids]*Area > 0.0) W_update[loc + 3] = W_wall[loc + 3];
					else W_update[loc + 3] = W_infty[loc + 3];

					if(vn[iFluids]-c[iFluids]*Area > 0.0) W_update[loc + 4] = W_wall[loc + 4];
					else W_update[loc + 4] = W_infty[loc + 4];
				}
			}

			/*--- conservative variables using characteristics ---*/
			for (iVar=0; iVar < nVar; iVar++) {
				U_update[iVar] = 0;
				for (jVar=0; jVar < nVar; jVar++)
					U_update[iVar] +=P_Matrix[iVar][jVar]*W_update[jVar];
			}

			/*--- Residual computation ---*/
			for (iFluids = 0; iFluids < nFluids; iFluids ++){
				loc = iFluids * (nDim +2);
				rho[iFluids] = U_update[loc + 0];
				rhoE[iFluids] = U_update[loc + nDim +1];
				sq_vel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					velocity[iFluids][iDim] = U_update[loc + iDim+1]/rho[iFluids];
					sq_vel +=velocity[iFluids][iDim]*velocity[iFluids][iDim];
				}
				c[iFluids] = sqrt(fabs(Gamma*Gamma_Minus_One*(rhoE[iFluids]/rho[iFluids] - 0.5*sq_vel)));
				pressure[iFluids] = (c[iFluids] * c[iFluids] * rho[iFluids]) / Gamma;
				enthalpy[iFluids] = (rhoE[iFluids] + pressure[iFluids]) / rho[iFluids];
			}

			solver->GetInviscidProjFlux(rho, velocity, pressure, enthalpy, UnitaryNormal, Residual);
			for	(iVar = 0; iVar < nVar; iVar++) Residual[iVar] = Residual[iVar]*Area;
			node[iPoint]->AddRes_Conv(Residual);

			/*--- In case we are doing a implicit computation ---*/
			if (implicit) {
				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
				solver->SetNormal(Vector);
				solver->SetConservative(U_wall, U_infty);
				solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			}

		}
	}

	delete [] UnitaryNormal;
	delete [] rho;
	delete [] rhoE;
	delete [] pressure;
	delete [] enthalpy;
	delete [] c;
	delete [] vn;
	for (iFluids = 0; iFluids < nFluids; iFluids ++ ) {
		delete [] velocity[iFluids];
	}
	delete [] velocity;
	delete [] U_wall; delete [] U_infty; delete [] U_update;
	delete [] W_wall; delete [] W_infty; delete [] W_update;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Matrix[iVar];
		delete [] invP_Matrix[iVar];
	}
	delete [] P_Matrix;
	delete [] invP_Matrix;
}

void CPlasmaMonatomicSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iDim, iVar, iFluids, loc;
	double Pressure, *Normal;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			/*--- Compute the projected residual ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

			for (iFluids = 0; iFluids < nFluids; iFluids ++ ) {
				loc = iFluids * (nDim+2);
				Pressure = node[iPoint]->GetPressure(iFluids);
				Residual[loc + 0] = 0; Residual[loc + nDim+1] = 0;
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[loc + iDim+1] = -Pressure*Normal[iDim];
			}
			/*--- Add value to the residual ---*/
			node[iPoint]->AddRes_Conv(Residual);

			/*			for (short index=0; index<15; index++) {
			 node[index]->Set_ResConv_Zero();
			 node[index]->Set_ResVisc_Zero();
			 }*/

			/*--- In case we are doing a implicit computation ---*/
			if (implicit) {
				for (iFluids = 0; iFluids < nFluids; iFluids ++ ) {
					double a2 = Gamma-1.0;
					double phi = 0.5*a2*node[iPoint]->GetVelocity2(iFluids);
					loc = (nDim+2)*iFluids;
					for (iVar = 0; iVar < nVar; iVar++) {
						Jacobian_i[loc + 0][loc + iVar] = 0.0;
						Jacobian_i[loc + nDim+1][loc + iVar] = 0.0;
					}
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_i[loc + iDim+1][loc + 0] = -phi*Normal[iDim];
						for (unsigned short jDim = 0; jDim < nDim; jDim++)
							Jacobian_i[loc + iDim+1][loc + jDim+1] = a2*node[iPoint]->GetVelocity(jDim,iFluids)*Normal[iDim];
						Jacobian_i[loc + iDim+1][loc + nDim+1] = -a2*Normal[iDim];
					}
				}
				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}

		}
	}
}

void CPlasmaMonatomicSolution::BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config,
		unsigned short val_marker, unsigned short val_mesh) {

#ifndef NO_MPI
	unsigned short iVar;
	unsigned long iVertex, iPoint;
	double *Conserv_Var, *Conserv_Undivided_Laplacian = NULL, **Conserv_Grad = NULL, *Grad_Limit = NULL;
	short SendRecv = config->GetMarker_All_SendRecv(val_marker);

	/*--- Send information  ---*/
	if (SendRecv > 0) {

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
			}

			MPI::COMM_WORLD.Bsend(&Buffer_Send_U,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Bsend(&Buffer_Send_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Bsend(&Buffer_Send_Sensor,nBuffer_Scalar,MPI::DOUBLE,send_to, 2);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_Lambda,nBuffer_Scalar,MPI::DOUBLE,send_to, 3);
			MPI::COMM_WORLD.Bsend(&Buffer_Send_Neighbor,nBuffer_Scalar,MPI::UNSIGNED_SHORT,send_to, 4);

		}
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
			unsigned long nBuffer = geometry->nVertex[val_marker]*nVar;
			int receive_from = abs(SendRecv);

			MPI::COMM_WORLD.Recv(&Buffer_Receive_U,nBuffer,MPI::DOUBLE,receive_from, 0);
			if (val_mesh == MESH_0) {
				MPI::COMM_WORLD.Recv(&Buffer_Receive_Ux,nBuffer,MPI::DOUBLE,receive_from, 1);
				MPI::COMM_WORLD.Recv(&Buffer_Receive_Uy,nBuffer,MPI::DOUBLE,receive_from, 2);
				if (nDim == 3) MPI::COMM_WORLD.Recv(&Buffer_Receive_Uz,nBuffer,MPI::DOUBLE,receive_from, 3);
				if (config->GetKind_SlopeLimit_Flow() != NONE) MPI::COMM_WORLD.Recv(&Buffer_Receive_Limit,nBuffer,MPI::DOUBLE,receive_from, 4);
			}

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
			unsigned long nBuffer_Vector = geometry->nVertex[val_marker]*nVar;
			unsigned long nBuffer_Scalar = geometry->nVertex[val_marker];
			int receive_from = abs(SendRecv);

			MPI::COMM_WORLD.Recv(&Buffer_Receive_U,nBuffer_Vector,MPI::DOUBLE,receive_from, 0);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Recv(&Buffer_Receive_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
			if (val_mesh == MESH_0) MPI::COMM_WORLD.Recv(&Buffer_Receive_Sensor,nBuffer_Scalar,MPI::DOUBLE,receive_from, 2);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_Lambda,nBuffer_Scalar,MPI::DOUBLE,receive_from, 3);
			MPI::COMM_WORLD.Recv(&Buffer_Receive_Neighbor,nBuffer_Scalar,MPI::UNSIGNED_SHORT,receive_from, 4);

			for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
				iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
				for (iVar = 0; iVar < nVar; iVar++) {
					node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVertex][iVar]);
					if (val_mesh == MESH_0) node[iPoint]->SetUndivided_Laplacian(iVar, Buffer_Receive_Undivided_Laplacian[iVertex][iVar]);
				}
				if (val_mesh == MESH_0) node[iPoint]->SetSensor(Buffer_Receive_Sensor[iVertex]);
				node[iPoint]->SetLambda(Buffer_Receive_Lambda[iVertex]);
				geometry->node[iPoint]->SetnPoint(Buffer_Receive_Neighbor[iVertex]);
			}
		}

	}
#endif
}

void CPlasmaMonatomicSolution::BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config,
		unsigned short val_marker, unsigned short val_mesh) {

#ifndef NO_MPI
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
			if (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT)
				for (iVar = 0; iVar < nVar; iVar++)
					Jacobian.DeleteValsRowi(iPoint*nVar+iVar);
		}
	}
#endif
}
