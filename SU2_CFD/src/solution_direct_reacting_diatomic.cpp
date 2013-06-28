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

CPlasmaDiatomicSolution::CPlasmaDiatomicSolution(void) : CSolution() { }

CPlasmaDiatomicSolution::CPlasmaDiatomicSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned long iPoint;
	unsigned short iVar, iDim, iSpecies;
	double Vel2;
	string mesh_filename, text_line;
	ifstream restart_file;
	
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	
	/*--- Define geometry constants in the solver structure ---*/
	nDim = geometry->GetnDim();
	nMonatomics = 2; 
	nDiatomics = 0; 
	nSpecies = nMonatomics + nDiatomics; 
	nVar = nMonatomics*(nDim+2) + nDiatomics*(nDim+3);
	node = new CVariable*[geometry->GetnPoint()];
	
	double GammaDiatomic = config->GetGammaDiatomic();
	double GammaMonatomic = config->GetGammaMonatomic();
	

	
	/*--- Define some auxiliar vector related with the residual ---*/
	Residual = new double[nVar];	Residual_Max = new double[nVar];
	Residual_i = new double[nVar];	Residual_j = new double[nVar];
	Res_Conv = new double[nVar];	Res_Visc = new double[nVar];
	
	/*--- Define some auxiliar vector related with the solution ---*/
	Solution = new double[nVar];
	Solution_i = new double[nVar]; Solution_j = new double[nVar];
	
	/*--- Define some auxiliar vector related with the geometry ---*/
	Vector = new double[nDim];
	Vector_i = new double[nDim]; Vector_j = new double[nDim];
		
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
	
	/*--- Flow inifinity initialization stuff ---*/
	Density_Inf = new double [nSpecies];
	Pressure_Inf = new double [nSpecies];
	Velocity_Inf = new double*[nSpecies];
	Mach_Inf = 	new double [nSpecies];
	Energy_Inf = new double [nSpecies];
	Energy_vib_Inf = new double [nDiatomics];
	Energy_Formation = new double [nSpecies];
	Enthalpy_Formation = new double [nSpecies];
	for ( iSpecies = 0; iSpecies < nSpecies; iSpecies++ ){
		Velocity_Inf[iSpecies] = new double [nDim];		
		Density_Inf[iSpecies] = config->GetDensity_FreeStreamND();
		Pressure_Inf[iSpecies] = config->GetPressure_FreeStreamND();
		for (iDim = 0; iDim < nDim; iDim++)
			Velocity_Inf[iSpecies][iDim] = config->GetVelocity_FreeStreamND()[iDim];
		Vel2 = 0.0; 
		for (iDim = 0; iDim < nDim; iDim++) {
			Vel2 += Velocity_Inf[iSpecies][iDim]*Velocity_Inf[iSpecies][iDim];
		}		
		if (iSpecies < nDiatomics) {
			Mach_Inf[iSpecies] = Velocity_Inf[iSpecies][0]/sqrt(fabs(GammaDiatomic*Pressure_Inf[iSpecies]/Density_Inf[iSpecies]));
			Energy_Inf[iSpecies] = Pressure_Inf[iSpecies]/(Density_Inf[iSpecies]*(GammaDiatomic-1.0))+0.5*Vel2;			
			Energy_vib_Inf[iSpecies] = 0.0;	
		}
		else {
			Mach_Inf[iSpecies] = Velocity_Inf[iSpecies][0]/sqrt(fabs(GammaMonatomic*Pressure_Inf[iSpecies]/Density_Inf[iSpecies]));
			Energy_Inf[iSpecies] = Pressure_Inf[iSpecies]/(Density_Inf[iSpecies]*(GammaMonatomic-1.0))+0.5*Vel2;			
		}
		Energy_Formation[iSpecies] = 0.0;
		Enthalpy_Formation[iSpecies] = 0.0;
	}
		
	/*--- Restart the solution from file information ---*/
	if (!restart) {
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
				node[iPoint] = new CPlasmaDiatomicVariable(Density_Inf, Velocity_Inf, Energy_Inf, Energy_vib_Inf, Energy_Formation, Enthalpy_Formation, nDim, nVar, nSpecies, 
																					 nMonatomics, nDiatomics, config);
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
/*		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
			if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
			node[iPoint] = new CPlasmaDiatomicVariable(Solution, nDim, nVar, config);
		}*/
		
		/*--- Close the restart file ---*/
		restart_file.close();

	}
	cout << "Leaving constructor";
}

CPlasmaDiatomicSolution::~CPlasmaDiatomicSolution(void) {
	unsigned short iVar, iDim;
	
	delete [] Residual;		delete [] Residual_Max;
	delete [] Residual_i;	delete [] Residual_j;
	delete [] Res_Conv;		delete [] Res_Visc;
	delete [] Solution;		delete [] Solution_i; delete [] Solution_j;
	delete [] Vector_i;		delete [] Vector_j;
	delete [] xsol;			delete [] rhs;
	
	delete [] Energy_Inf;	
	delete [] Energy_vib_Inf;	
	delete [] Density_Inf;
	delete [] Pressure_Inf;	
	delete [] Mach_Inf;
	
	delete [] Gamma;
	delete [] Energy_Formation;
	delete [] Enthalpy_Formation;
	delete [] Molecular_Mass;
	
	for (iVar = 0; iVar < nFluids; iVar ++) {
		delete [] Velocity_Inf[iVar];
	}
	delete[] Velocity_Inf;
	
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
	
	/*	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
	 delete [] node[iPoint];
	 delete [] node; */
}

void CPlasmaDiatomicSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		/*--- Compute squared velocity, sound velocity, pressure, and enthalpy ---*/
		node[iPoint]->SetVelocity2();
		node[iPoint]->SetSoundSpeed(config->GetGammaMonatomic(), config->GetGammaDiatomic());
		node[iPoint]->SetPressure(config->GetGammaMonatomic(), config->GetGammaDiatomic());
		node[iPoint]->SetEnthalpy();
		
		/*--- Initialize the convective and viscous residual vector ---*/
		node[iPoint]->Set_ResConv_Zero();
		if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit) 
			node[iPoint]->Set_ResVisc_Zero();
	}
		
	/*--- Inicialize the jacobian matrices ---*/
	if (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT)
		Jacobian.SetValZero();
}

void CPlasmaDiatomicSolution::SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {
	double *Normal, Area, dV, Mean_SoundSpeed, Mean_ProjVel, Lambda, Local_Delta_Time, Global_Delta_Time = 1E6;
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker, iSpecies;
	double Lambda_iSpecies;
	
	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for ( iSpecies = 0; iSpecies < nSpecies; iSpecies++ )
			node[iPoint]->SetMax_Lambda_Inv(0.0, iSpecies);
		
		node[iPoint]->SetVelocity2();
		node[iPoint]->SetSoundSpeed(config->GetGammaMonatomic(), config->GetGammaDiatomic());
	}
	
	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		
		/*--- Point identification, Normal vector and area ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); 
		jPoint = geometry->edge[iEdge]->GetNode(1);
		Normal = geometry->edge[iEdge]->GetNormal();
		Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
		
		/*--- Mean Values ---*/
		for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
			Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal,iSpecies, nDiatomics) + node[jPoint]->GetProjVel(Normal,iSpecies,nDiatomics));
			Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed(iSpecies) + node[jPoint]->GetSoundSpeed(iSpecies));
			Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed * Area;
			
			/*--- Inviscid contribution ---*/
			node[iPoint]->AddMax_Lambda_Inv(Lambda, iSpecies);
			node[jPoint]->AddMax_Lambda_Inv(Lambda, iSpecies);
		}
	}
	
	/*--- Loop boundary edges ---*/
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) { 
		for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			
			/*--- Point identification, Normal vector and area ---*/
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
			
			for ( iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {	

				/*--- Mean Values ---*/
				Mean_ProjVel = node[iPoint]->GetProjVel(Normal,iSpecies,nDiatomics);
				Mean_SoundSpeed = node[iPoint]->GetSoundSpeed(iSpecies);	
				
				/*--- Inviscid contribution ---*/
				Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed * Area;
//				cout << "Lambda: " << Lambda << endl;
//				cout << "Boundary ProjVel: " << Mean_ProjVel << endl; cin.get();
				node[iPoint]->AddMax_Lambda_Inv(Lambda, iSpecies);
			}			
		}
	}
	
	/*--- Each element uses their own speed ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		dV = geometry->node[iPoint]->GetVolume();
		
		Lambda_iSpecies = node[iPoint]->GetMax_Lambda_Inv(0);
		for (iSpecies = 1; iSpecies < nSpecies; iSpecies++)
			Lambda_iSpecies = max(Lambda_iSpecies, node[iPoint]->GetMax_Lambda_Inv(iSpecies));
		
		Local_Delta_Time = config->GetCFL(iMesh)*dV / Lambda_iSpecies;
		if (config->GetUnsteady_Simulation() == TIME_STEPPING) Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
		else node[iPoint]->SetDelta_Time(Local_Delta_Time);
	}
	
	/*--- For exact time solution use the minimum delta time of the whole mesh ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
		for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
			node[iPoint]->SetDelta_Time(Global_Delta_Time);
	}
}

void CPlasmaDiatomicSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
																			CConfig *config, unsigned short iMesh) {
	
	double *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
//	bool high_order_diss = ((config->GetKind_Upwind() == ROE_2ND) && (iMesh == MESH_0));
	
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		
		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());
		
		/*--- Set conservative variables w/o reconstruction ---*/
		U_i = node[iPoint]->GetSolution();
		U_j = node[jPoint]->GetSolution();
		solver->SetConservative(U_i, U_j);
		
		
		/*--- Compute the residual ---*/
		solver->SetResidual(Res_Conv, Jacobian_i, Jacobian_j, config);
		
//		if ( (iPoint == 0) || (jPoint == 0)){
//		cout << iEdge << " " << iPoint<< " " << jPoint <<  endl;
//			cout << Res_Conv[0] << " " << Res_Conv[1] << " " <<Res_Conv[2] << " " <<Res_Conv[3] << endl<< endl;
// 			cout << geometry->node[iPoint]->GetCoord(0)<< " " << geometry->node[iPoint]->GetCoord(1) << endl; cin.get();
//		}
/*		cout << U_i[0] << " " << U_i[1] << " " <<U_i[2] << " " <<U_i[3] << endl<< endl;
		cout << U_j[0] << " " << U_j[1] << " " <<U_j[2] << " " <<U_j[3] << endl<< endl;
		
		cout << Jacobian_i[0][0] << " " << Jacobian_i[0][1] << " " <<Jacobian_i[0][2] << " " <<Jacobian_i[0][3] << endl;
		cout << Jacobian_i[1][0] << " " << Jacobian_i[1][1] << " " <<Jacobian_i[1][2] << " " <<Jacobian_i[1][3] << endl;
		cout << Jacobian_i[2][0] << " " << Jacobian_i[2][1] << " " <<Jacobian_i[2][2] << " " <<Jacobian_i[2][3] << endl;
		cout << Jacobian_i[3][0] << " " << Jacobian_i[3][1] << " " <<Jacobian_i[3][2] << " " <<Jacobian_i[3][3] << endl;*/
		
	
		/*--- Update residual value ---*/
		node[iPoint]->AddRes_Conv(Res_Conv);		
		node[jPoint]->SubtractRes_Conv(Res_Conv);
		
//		if ((iPoint== 0) ) {cout << node[iPoint]->GetResConv()[2]; cin.get();}
//		if ((jPoint== 0)) {cout << node[jPoint]->GetResConv()[2]; cin.get();}

		/*--- Set implicit jacobians ---*/
		if (implicit) {		
			Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
			Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
			Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);
		}
	}
}


void CPlasmaDiatomicSolution::SourcePieceWise_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																							 CConfig *config, unsigned short iMesh) {
	unsigned short iVar;
	unsigned long iPoint;
	double **Gradient;
		
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
		solver->SetResidual(Residual, Jacobian_i, config);
				
		/*--- Add Residual ---*/
		node[iPoint]->SubtractRes_Conv(Residual);
		
		if (implicit)
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
	}
}

void CPlasmaDiatomicSolution::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
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

void CPlasmaDiatomicSolution::SetSpectral_Radius(CGeometry *geometry, CConfig *config) {
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

void CPlasmaDiatomicSolution::RungeKutta_Iteration(CGeometry *geometry, CSolution **solution_container, 
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

void CPlasmaDiatomicSolution::ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
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

void CPlasmaDiatomicSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
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
//		cout << iPoint << " " <<Delta << endl;
		
		Jacobian.AddVal2Diag(iPoint,Delta);
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
			Res = local_ResConv[iVar]+local_ResVisc[iVar];
			rhs[total_index] = -(Res+local_TruncationError[iVar]);
//			cout << "RHS: " << iVar <<" "<< local_ResConv[iVar]<< " " <<local_ResVisc[iVar] << endl; cin.get();

			xsol[total_index] = 0.0;
			AddRes_Max( iVar, Res*Res*Vol );
		}
	}
	
//	cin.get();

	/*--- Solve the system ---*/
//	Jacobian.SGSSolution(rhs, xsol, 1e-9, 1000, true);
	Jacobian.LU_SGSIteration(rhs, xsol);
	
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, xsol[iPoint*nVar+iVar]);
		}
	}
	
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, sqrt(GetRes_Max(iVar)));
}

void CPlasmaDiatomicSolution::SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config) {
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

void CPlasmaDiatomicSolution::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config) {
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

void CPlasmaDiatomicSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, 
																						CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iDim, iSpecies, iVar, loc;
	double Pressure, *Normal;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	
	/*--- Bucle over all the vertices ---*/
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {
				
				if (nDim == 2) {
					if ( iSpecies < nDiatomics ) loc = 5*iSpecies;
					else loc = 5*nDiatomics + 4*(iSpecies-nDiatomics);	
					
					Residual[loc+0] = 0.0;
					Residual[loc+nDim+1] = 0.0;
					if ( iSpecies < nDiatomics ) 	
						Residual[loc+nDim+2] = 0.0;
					Pressure = node[iPoint]->GetPressure(iSpecies);
					for (iDim = 0; iDim < nDim; iDim++)
						Residual[loc + iDim+1] = -Pressure*Normal[iDim];
				}	
//				cout << iSpecies << endl;
//				cout << "Boundary Res: " << Residual[1] << " " << Residual[0] << endl; cin.get();
			}
			
			/*--- Add value to the residual ---*/
			node[iPoint]->AddRes_Conv(Residual);
			
//			 cout << Pressure << " " << Residual[0] << " " << Residual[1] << " " <<Residual[2] << " " <<Residual[3] << endl<< endl;
			 
			 
/*			 cout << Jacobian_i[0][0] << " " << Jacobian_i[0][1] << " " <<Jacobian_i[0][2] << " " <<Jacobian_i[0][3] << endl;
			 cout << Jacobian_i[1][0] << " " << Jacobian_i[1][1] << " " <<Jacobian_i[1][2] << " " <<Jacobian_i[1][3] << endl;
			 cout << Jacobian_i[2][0] << " " << Jacobian_i[2][1] << " " <<Jacobian_i[2][2] << " " <<Jacobian_i[2][3] << endl;
			 cout << Jacobian_i[3][0] << " " << Jacobian_i[3][1] << " " <<Jacobian_i[3][2] << " " <<Jacobian_i[3][3] << endl;*/
			
			/*--- In case we are doing a implicit computation ---*/

			if (implicit) {
				double a2;
				double phi;
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
					if ( iSpecies < nDiatomics ) a2 = config->GetGammaDiatomic()-1.0;
					else  a2 = config->GetGammaMonatomic()-1.0;
					phi = 0.5*a2*node[iPoint]->GetVelocity2(iSpecies);
					
					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);	

					for (iVar = 0; iVar < nVar; iVar++) {
						Jacobian_i[loc + 0][loc + iVar] = 0.0;
						Jacobian_i[loc + nDim+1][loc + iVar] = 0.0;
					}
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_i[loc + iDim+1][loc + 0] = -phi*Normal[iDim];
						for (unsigned short jDim = 0; jDim < nDim; jDim++)
							Jacobian_i[loc + iDim+1][loc + jDim+1] = a2*node[iPoint]->GetVelocity(jDim,iSpecies)*Normal[iDim];
						Jacobian_i[loc + iDim+1][loc + nDim+1] = -a2*Normal[iDim];
					}
				}
				
/*					cout << Jacobian_i[0][0] << " " << Jacobian_i[0][1] << " " <<Jacobian_i[0][2] << " " <<Jacobian_i[0][3] << endl;
				 cout << Jacobian_i[1][0] << " " << Jacobian_i[1][1] << " " <<Jacobian_i[1][2] << " " <<Jacobian_i[1][3] << endl;
				 cout << Jacobian_i[2][0] << " " << Jacobian_i[2][1] << " " <<Jacobian_i[2][2] << " " <<Jacobian_i[2][3] << endl;
				 cout << Jacobian_i[3][0] << " " << Jacobian_i[3][1] << " " <<Jacobian_i[3][2] << " " <<Jacobian_i[3][3] << endl;
				cin.get();*/
				
				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}
		}
	}
}

void CPlasmaDiatomicSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
																					 CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iDim, iVar, iSpecies, loc;
	double Pressure, *Normal;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	
	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {
				
				if (nDim == 2) {
					if ( iSpecies < nDiatomics ) loc = 5*iSpecies;
					else loc = 5*nDiatomics + 4*(iSpecies-nDiatomics);	
					
					Residual[loc+0] = 0.0;
					Residual[loc+nDim+1] = 0.0;
					if ( iSpecies < nDiatomics ) 	
						Residual[loc+nDim+2] = 0.0;
					Pressure = node[iPoint]->GetPressure(iSpecies);
					for (iDim = 0; iDim < nDim; iDim++)
						Residual[loc + iDim+1] = -Pressure*Normal[iDim];
				}	
				//				cout << iSpecies << endl;
				//				cout << "Boundary Res: " << Residual[1] << " " << Residual[0] << endl; cin.get();
			}
			
			/*--- Add value to the residual ---*/
			node[iPoint]->AddRes_Conv(Residual);
			
//			cout << "Resid "<< Pressure << " " << Residual[0] << " " << Residual[1] << " " <<Residual[2] << " " <<Residual[3] << endl<< endl;
			
			
			/*			 cout << Jacobian_i[0][0] << " " << Jacobian_i[0][1] << " " <<Jacobian_i[0][2] << " " <<Jacobian_i[0][3] << endl;
			 cout << Jacobian_i[1][0] << " " << Jacobian_i[1][1] << " " <<Jacobian_i[1][2] << " " <<Jacobian_i[1][3] << endl;
			 cout << Jacobian_i[2][0] << " " << Jacobian_i[2][1] << " " <<Jacobian_i[2][2] << " " <<Jacobian_i[2][3] << endl;
			 cout << Jacobian_i[3][0] << " " << Jacobian_i[3][1] << " " <<Jacobian_i[3][2] << " " <<Jacobian_i[3][3] << endl;*/
			
			/*--- In case we are doing a implicit computation ---*/
			
			if (implicit) {
				double a2;
				double phi;
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
					if ( iSpecies < nDiatomics ) a2 = config->GetGammaDiatomic()-1.0;
					else  a2 = config->GetGammaMonatomic()-1.0;
					phi = 0.5*a2*node[iPoint]->GetVelocity2(iSpecies);
					
					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);	
					
					for (iVar = 0; iVar < nVar; iVar++) {
						Jacobian_i[loc + 0][loc + iVar] = 0.0;
						Jacobian_i[loc + nDim+1][loc + iVar] = 0.0;
					}
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_i[loc + iDim+1][loc + 0] = -phi*Normal[iDim];
						for (unsigned short jDim = 0; jDim < nDim; jDim++)
							Jacobian_i[loc + iDim+1][loc + jDim+1] = a2*node[iPoint]->GetVelocity(jDim,iSpecies)*Normal[iDim];
						Jacobian_i[loc + iDim+1][loc + nDim+1] = -a2*Normal[iDim];
					}
				}
				
				/*					cout << Jacobian_i[0][0] << " " << Jacobian_i[0][1] << " " <<Jacobian_i[0][2] << " " <<Jacobian_i[0][3] << endl;
				 cout << Jacobian_i[1][0] << " " << Jacobian_i[1][1] << " " <<Jacobian_i[1][2] << " " <<Jacobian_i[1][3] << endl;
				 cout << Jacobian_i[2][0] << " " << Jacobian_i[2][1] << " " <<Jacobian_i[2][2] << " " <<Jacobian_i[2][3] << endl;
				 cout << Jacobian_i[3][0] << " " << Jacobian_i[3][1] << " " <<Jacobian_i[3][2] << " " <<Jacobian_i[3][3] << endl;
				 cin.get();*/
				
				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}
		}
	}
}

void CPlasmaDiatomicSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim, iSpecies;
	unsigned short loc = 0;
	double *U_wall, *U_infty;
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	
	U_wall = new double[nVar]; U_infty = new double[nVar];
	
	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_wall[iVar] = node[iPoint]->GetSolution(iVar);
			
			/*--- Solution at the infinity ---*/
			for(iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);	
				if (nDim == 2) {
					U_infty[loc + 0] = GetDensity_Inf(iSpecies);
					U_infty[loc + 1] = GetDensity_Velocity_Inf(0,iSpecies);
					U_infty[loc + 2] = GetDensity_Velocity_Inf(1,iSpecies);
					U_infty[loc + 3] = GetDensity_Energy_Inf(iSpecies);
					if ( iSpecies < nDiatomics ) U_infty[loc + 4] = GetDensity_Energy_vib_Inf(iSpecies);
				}
				if (nDim == 3) {	
					U_infty[loc + 0] = GetDensity_Inf(iSpecies);
					U_infty[loc + 1] = GetDensity_Velocity_Inf(0,iSpecies);
					U_infty[loc + 2] = GetDensity_Velocity_Inf(1,iSpecies);
					U_infty[loc + 3] = GetDensity_Velocity_Inf(2,iSpecies);
					U_infty[loc + 4] = GetDensity_Energy_Inf(iSpecies);
					if ( iSpecies < nDiatomics ) U_infty[loc + 5] = GetDensity_Energy_vib_Inf(iSpecies);
				}
			}	
			
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);			
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			
			solver->SetNormal(Vector);
			solver->SetConservative(U_wall, U_infty);
			solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			node[iPoint]->AddRes_Conv(Residual);
			
			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			
		}
	}
	delete [] U_wall; delete [] U_infty;

}

void CPlasmaDiatomicSolution::BC_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
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

void CPlasmaDiatomicSolution::BC_InterProcessor(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
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
