/*!
 * \file solution_adjoint_plasma.cpp
 * \brief Main subrotuines for solving adjoint problems (Euler, Navier-Stokes, etc.).
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

CAdjPlasmaSolution::CAdjPlasmaSolution(void) : CSolution() { }

CAdjPlasmaSolution::CAdjPlasmaSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned long iPoint, iVertex;
	string text_line, mesh_filename;
	unsigned short iDim = 0, iVar, iMarker;
	ifstream restart_file;
	string filename, AdjExt;

	bool restart = config->GetRestart();
	bool axisymmetric = config->GetAxisymmetric();

	/*--- Define geometry constans in the solver structure ---*/
	nDim        = geometry->GetnDim();
	nMonatomics = config->GetnMonatomics();
	nDiatomics  = config->GetnDiatomics();
	nSpecies    = config->GetnSpecies();
	nVar        = nMonatomics*(nDim+2) + nDiatomics*(nDim+3);
	node = new CVariable*[geometry->GetnPoint()];

	/*--- Define some auxiliary vectors related to the residual ---*/
	Residual = new double[nVar];	  Residual_RMS = new double[nVar];
	Residual_i = new double[nVar];	Residual_j = new double[nVar];
	Res_Conv_i = new double[nVar];	Res_Visc_i = new double[nVar];
	Res_Conv_j = new double[nVar];	Res_Visc_j = new double[nVar];
	Res_Sour_i = new double[nVar];	Res_Sour_j = new double[nVar];
  Residual_Max = new double[nVar]; Point_Max = new unsigned long[nVar];

	Residual_Chemistry = new double[nVar]; Residual_MomentumExch = new double[nVar];
	Residual_ElecForce = new double[nVar]; Residual_EnergyExch = new double[nVar];
	Res_Conv = new double[nVar];	Res_Visc = new double[nVar];	Res_Sour = new double[nVar];
	if (axisymmetric) {
		Residual_Axisymmetric = new double[nVar];
	}
	

	/*--- Define some auxiliary vectors related to the solution ---*/
	Solution   = new double[nVar];
	Solution_i = new double[nVar]; Solution_j = new double[nVar];
  
	/*--- Define some auxiliary vectors related to the geometry ---*/
//  Vector = new double[nDim];
	Vector_i = new double[nDim]; Vector_j = new double[nDim];

	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT) {
		Jacobian_ii = new double* [nVar];
		Jacobian_ij = new double* [nVar];
		Jacobian_ji = new double* [nVar];
		Jacobian_jj = new double* [nVar];
		Jacobian_Chemistry    = new double* [nVar];
		Jacobian_ElecForce    = new double* [nVar];
		Jacobian_MomentumExch = new double* [nVar];
		Jacobian_EnergyExch   = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_ii[iVar] = new double [nVar];
			Jacobian_ij[iVar] = new double [nVar];
			Jacobian_ji[iVar] = new double [nVar];
			Jacobian_jj[iVar] = new double [nVar];
			Jacobian_Chemistry[iVar]    = new double [nVar];
			Jacobian_ElecForce[iVar]    = new double [nVar];
			Jacobian_MomentumExch[iVar] = new double [nVar];
			Jacobian_EnergyExch[iVar]   = new double [nVar];
			
		}
		if (axisymmetric) {
			Jacobian_Axisymmetric = new double*[nVar];
			for (iVar = 0; iVar < nVar; iVar++)
				Jacobian_Axisymmetric[iVar] = new double[nVar];
		}
		Initialize_SparseMatrix_Structure(&Jacobian, nVar, nVar, geometry, config);
		xsol = new double [geometry->GetnPoint()*nVar];
		rhs  = new double [geometry->GetnPoint()*nVar];
	}

	/*--- Computation of gradients by least squares ---*/
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim]; 
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nVar]; 
		for (iVar = 0; iVar < nVar; iVar++)
			cvector[iVar] = new double [nDim];
	}
	
	/*--- Sensitivity definition and coefficient in all the markers ---*/
	CSensitivity = new double* [config->GetnMarker_All()];
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		CSensitivity[iMarker] = new double [geometry->nVertex[iMarker]];
	}
	Sens_Geo  = new double[config->GetnMarker_All()];
	Sens_Mach = new double[config->GetnMarker_All()];
	Sens_AoA  = new double[config->GetnMarker_All()];
	Sens_Press = new double[config->GetnMarker_All()];
	Sens_Temp  = new double[config->GetnMarker_All()];
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Sens_Geo[iMarker]  = 0.0;
		Sens_Mach[iMarker] = 0.0;
		Sens_AoA[iMarker]  = 0.0;
		Sens_Press[iMarker] = 0.0;
		Sens_Temp[iMarker]  = 0.0;
		for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
			CSensitivity[iMarker][iVertex] = 0.0;
	}

	/*--- Adjoint flow at the inifinity, initialization stuff ---*/
	PsiRho_Inf = 0.0; PsiE_Inf   = 0.0; PsiEvib_Inf = 0.0;
	Phi_Inf    = new double [nDim];
	Phi_Inf[0] = 0.0; Phi_Inf[1] = 0.0;
	if (nDim == 3) Phi_Inf[2] = 0.0;		
  
	if (!restart || geometry->GetFinestMGLevel() == false) {
		/*--- Restart the solution from infinity ---*/
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CAdjPlasmaVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, PsiEvib_Inf, nDim, nVar, config);
	}
	else {
		
		/*--- Restart the solution from file information ---*/
		mesh_filename = config->GetSolution_AdjFileName();
		
		/*--- Change the name, depending of the objective function ---*/
		filename.assign(mesh_filename);
		filename.erase (filename.end()-4, filename.end());
		switch (config->GetKind_ObjFunc()) {
			case DRAG_COEFFICIENT: AdjExt = "_cd.dat"; break;
			case LIFT_COEFFICIENT: AdjExt = "_cl.dat"; break;
			case SIDEFORCE_COEFFICIENT: AdjExt = "_csf.dat"; break;
			case PRESSURE_COEFFICIENT: AdjExt = "_cp.dat"; break;
			case MOMENT_X_COEFFICIENT: AdjExt = "_cmx.dat"; break;
			case MOMENT_Y_COEFFICIENT: AdjExt = "_cmy.dat"; break;
			case MOMENT_Z_COEFFICIENT: AdjExt = "_cmz.dat"; break;
			case EFFICIENCY: AdjExt = "_eff.dat"; break;
			case EQUIVALENT_AREA: AdjExt = "_ea.dat"; break;
			case NEARFIELD_PRESSURE: AdjExt = "_nfp.dat"; break;
      case FORCE_X_COEFFICIENT: AdjExt = "_cfx.dat"; break;
			case FORCE_Y_COEFFICIENT: AdjExt = "_cfy.dat"; break;
			case FORCE_Z_COEFFICIENT: AdjExt = "_cfz.dat"; break;
      case THRUST_COEFFICIENT: AdjExt = "_ct.dat"; break;
      case TORQUE_COEFFICIENT: AdjExt = "_cq.dat"; break;
      case FIGURE_OF_MERIT: AdjExt = "_merit.dat"; break;
			case FREESURFACE: AdjExt = "_fs.dat"; break;
      case NOISE: AdjExt = "_fwh.dat"; break;
		}
		filename.append(AdjExt);
		restart_file.open(filename.data(), ios::in);
    
		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
			cout << "There is no adjoint restart file!!" << endl;
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
    long iPoint_Local; unsigned long iPoint_Global = 0; unsigned long index;
    
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

        /*--- First value is the point index, then the conservative vars. ---*/
        point_line >> index;
        for (iVar = 0; iVar < nVar; iVar++)
          point_line >> Solution[iVar];
        
        node[iPoint_Local] = new CAdjPlasmaVariable(Solution, nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
      node[iPoint] = new CAdjPlasmaVariable(Solution, nDim, nVar, config);
    }
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
	}
  
  /*--- Define solver parameters needed for execution of destructor ---*/
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTERED) space_centered = true;
	else space_centered = false;
  
}

CAdjPlasmaSolution::~CAdjPlasmaSolution(void) {
	unsigned short iVar, iDim;

	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Jacobian_ii[iVar]; delete [] Jacobian_ij[iVar];
		delete [] Jacobian_ji[iVar]; delete [] Jacobian_jj[iVar];
	}
	delete [] Jacobian_ii; delete [] Jacobian_ij;
	delete [] Jacobian_ji; delete [] Jacobian_jj;
	
	delete [] Residual; delete [] Residual_Max;
	delete [] Residual_i; delete [] Residual_j;
	delete [] Res_Conv_i; delete [] Res_Visc_i;
	delete [] Res_Conv_j; delete [] Res_Visc_j;
	delete [] Res_Sour_i; delete [] Res_Sour_j;
	delete [] Solution; 
	delete [] Solution_i; delete [] Solution_j;
//  delete [] Vector;
	delete [] Vector_i; delete [] Vector_j;
	delete [] xsol; delete [] rhs;
	delete [] Sens_Geo; delete [] Sens_Mach;
	delete [] Sens_AoA; delete [] Sens_Press;
	delete [] Sens_Temp; delete [] Phi_Inf;
	
  if (space_centered) {
		delete [] p1_Und_Lapl;
		delete [] p2_Und_Lapl;
	}
  
	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;
	
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
	
/*	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		delete [] CSensitivity[iMarker];
	 delete [] CSensitivity; */
	
/*	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		delete [] node[iPoint];
	delete [] node; */
}

void CAdjPlasmaSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);  
  
	/*--- Residual initialization ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		
		/*--- Initialize the convective residual vector ---*/
		node[iPoint]->Set_ResConv_Zero();
		
		/*--- Initialize the source residual vector ---*/
		node[iPoint]->Set_ResSour_Zero();
		
		/*--- Initialize the viscous residual vector ---*/
		node[iPoint]->Set_ResVisc_Zero(); 
		
	}
	
	/*--- Implicit solution ---*/
	if (implicit) Jacobian.SetValZero();
	
}

void CAdjPlasmaSolution::SetForceProj_Vector(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	double *ForceProj_Vector, x = 0.0, y = 0.0, z = 0.0, *Normal;
	unsigned short iMarker;
	unsigned long iVertex, iPoint;
	double Alpha      = (config->GetAoA()*PI_NUMBER)/180.0;
	double Beta       = (config->GetAoS()*PI_NUMBER)/180.0;

	ForceProj_Vector = new double [nDim];

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if ((config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE) && 
			(config->GetMarker_All_Monitoring(iMarker) == YES))
			for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				
				x = geometry->node[iPoint]->GetCoord(0); 
				y = geometry->node[iPoint]->GetCoord(1);
				if (nDim == 3) z = geometry->node[iPoint]->GetCoord(2);
				
				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				switch (config->GetKind_ObjFunc()) {	
					case DRAG_COEFFICIENT :
						if (nDim == 2) { ForceProj_Vector[0] = cos(Alpha); ForceProj_Vector[1] = sin(Alpha); }
						if (nDim == 3) { ForceProj_Vector[0] = cos(Alpha)*cos(Beta); ForceProj_Vector[1] = sin(Beta); ForceProj_Vector[2] = sin(Alpha)*cos(Beta); }
						break;
					case LIFT_COEFFICIENT :
						if (nDim == 2) { ForceProj_Vector[0] = -sin(Alpha); ForceProj_Vector[1] = cos(Alpha); }
						if (nDim == 3) { ForceProj_Vector[0] = -sin(Alpha); ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = cos(Alpha); }
						break;
					case SIDEFORCE_COEFFICIENT :
						if (nDim == 2) { cout << "This functional is not possible in 2D!!" << endl;
							cout << "Press any key to exit..." << endl; cin.get(); exit(1);
						}
						if (nDim == 3) { ForceProj_Vector[0] = -sin(Beta) * cos(Alpha); ForceProj_Vector[1] = cos(Beta); ForceProj_Vector[2] = -sin(Beta) * sin(Alpha); }
						break;
          case FORCE_X_COEFFICIENT :
						if (nDim == 2) { ForceProj_Vector[0] = 1.0; ForceProj_Vector[1] = 0.0; }
						if (nDim == 3) { ForceProj_Vector[0] = 1.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 0.0; }
						break;
          case FORCE_Y_COEFFICIENT :
						if (nDim == 2) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 1.0; }
						if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 1.0; ForceProj_Vector[2] = 0.0; }
						break;
          case FORCE_Z_COEFFICIENT :
						if (nDim == 2) {cout << "This functional is not possible in 2D!!" << endl;
							cout << "Press any key to exit..." << endl;
							cin.get(); exit(1);
            }
						if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 1.0; }
						break;
				}
        
        /*--- Store the force projection vector at this node ---*/
				node[iPoint]->SetForceProj_Vector(ForceProj_Vector);
			}
	
	delete [] ForceProj_Vector;
}

void CAdjPlasmaSolution::Centered_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																			 CConfig *config, unsigned short iMesh, unsigned short iRKStep) { 	
	unsigned long iEdge, iPoint, jPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	bool dissipation = ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit);
	unsigned short iSpecies;
	
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		
		/*--- Points in edge, set normal vectors, and number of neighbors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());
		
		solver->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());
		
		/*--- Adjoint variables w/o reconstruction ---*/
		solver->SetAdjointVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
		
		/*--- Set conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(), 
														solution_container[PLASMA_SOL]->node[jPoint]->GetSolution());
		
		
		for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
			solver->SetPressure(solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies), 
													solution_container[PLASMA_SOL]->node[jPoint]->GetPressure(iSpecies), iSpecies);
			solver->SetSoundSpeed(solution_container[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies), 
														solution_container[PLASMA_SOL]->node[jPoint]->GetSoundSpeed(iSpecies),iSpecies);
			solver->SetEnthalpy(solution_container[PLASMA_SOL]->node[iPoint]->GetEnthalpy(iSpecies), 
													solution_container[PLASMA_SOL]->node[jPoint]->GetEnthalpy(iSpecies),iSpecies);
			solver->SetLambda(solution_container[PLASMA_SOL]->node[iPoint]->GetLambda(iSpecies), 
												solution_container[PLASMA_SOL]->node[jPoint]->GetLambda(iSpecies),iSpecies);
		}
				
		/*--- Compute residuals ---*/
		solver->SetResidual(Res_Conv_i, Res_Visc_i, Res_Conv_j, Res_Visc_j, 
												Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
				
//		cout << Res_Visc_i[0] <<" "<< Res_Visc_i[1] <<" "<< Res_Visc_i[2] <<" "<< Res_Visc_i[3] <<" "<< Res_Visc_i[4] <<" "<< Res_Visc_i[5] <<" "<< Res_Visc_i[6] <<" "<< Res_Visc_i[7] <<endl;
//		cout << Res_Visc_j[0] <<" "<< Res_Visc_j[1] <<" "<< Res_Visc_j[2] <<" "<< Res_Visc_j[3] <<" "<< Res_Visc_j[4] <<" "<< Res_Visc_j[5] <<" "<< Res_Visc_j[6] <<" "<< Res_Visc_j[7] <<endl;

		/*--- Update convective and artificial dissipation residuals ---*/
		node[iPoint]->SubtractRes_Conv(Res_Conv_i);
		node[jPoint]->SubtractRes_Conv(Res_Conv_j);
		if (dissipation) {
			node[iPoint]->SubtractRes_Visc(Res_Visc_i);
			node[jPoint]->SubtractRes_Visc(Res_Visc_j);
		}
		
		/*--- Implicit contribution to the residual ---*/
		if (implicit) {
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
			Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
			Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ji);
			Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_jj);
		}		
	}
}

void CAdjPlasmaSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
																				 CConfig *config, unsigned short iMesh) {
	double *Psi_i, *Psi_j, *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
  
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());

		/*--- Adjoint variables w/o reconstruction ---*/
		Psi_i = node[iPoint]->GetSolution(); Psi_j = node[jPoint]->GetSolution();
		solver->SetAdjointVar(Psi_i, Psi_j);
		
		/*--- Conservative variables w/o reconstruction ---*/
		U_i = solution_container[PLASMA_SOL]->node[iPoint]->GetSolution();
		U_j = solution_container[PLASMA_SOL]->node[jPoint]->GetSolution();
		solver->SetConservative(U_i, U_j);
		
		solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
		
		/*--- Add and Subtract Residual ---*/
		node[iPoint]->SubtractRes_Conv(Residual_i);
		node[jPoint]->SubtractRes_Conv(Residual_j);
		
    /*--- Implicit contribution to the residual ---*/
		if (implicit) {
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
			Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
			Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ji);
			Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_jj);
		}
	}
}

void CAdjPlasmaSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CNumerics *second_solver,
																				 CConfig *config, unsigned short iMesh) {
	
	unsigned short iVar, jVar, iSpecies;
	unsigned long iPoint;
	double *Psi_i;
	bool implicit = (config->GetKind_TimeIntScheme_AdjPlasma() == EULER_IMPLICIT);
	bool axisymmetric = config->GetAxisymmetric();
	
	double *Temperature_tr_i, *Temperature_vib_i, *Pressure_i;
	
	Temperature_tr_i = new double [nSpecies];
	Temperature_vib_i = new double [nSpecies];
	Pressure_i = new double[nSpecies];
	
	for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
		Temperature_tr_i[iSpecies] = 0.0;
		Temperature_vib_i[iSpecies] = 0.0;
		Pressure_i[iSpecies] = 0.0;
	}
	
	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		for (iVar = 0; iVar < nVar; iVar++) {
			Residual[iVar] = 0.0;
			Residual_Chemistry[iVar] = 0.0;
			Residual_MomentumExch[iVar] = 0.0;
			Residual_ElecForce[iVar] = 0.0;
			Residual_EnergyExch[iVar] = 0.0;
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_Chemistry[iVar][jVar] = 0.0;
				Jacobian_MomentumExch[iVar][jVar] = 0.0;
				Jacobian_ElecForce[iVar][jVar] = 0.0;
				Jacobian_EnergyExch[iVar][jVar] = 0.0;
			}
		}
		if (axisymmetric) {
			for (iVar = 0; iVar < nVar; iVar++) {
				Residual_Axisymmetric[iVar] = 0.0;
				for (jVar = 0; jVar < nVar; jVar++) 
					Jacobian_Axisymmetric[iVar][jVar] = 0.0;
			}
		}
    
    /*--- Set y coordinate ---*/
    solver->SetCoord(geometry->node[iPoint]->GetCoord(),geometry->node[iPoint]->GetCoord());
    
    /*--- Set solution  ---*/
    solver->SetConservative(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(), solution_container[PLASMA_SOL]->node[iPoint]->GetSolution());
    
    /*--- Set control volume ---*/
    solver->SetVolume(geometry->node[iPoint]->GetVolume());
    
    /*--- Set temperature ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
      Temperature_tr_i[iSpecies] = solution_container[PLASMA_SOL]->node[iPoint]->GetTemperature_tr(iSpecies);
      Temperature_vib_i[iSpecies] = solution_container[PLASMA_SOL]->node[iPoint]->GetTemperature_vib(iSpecies);
      Pressure_i[iSpecies] = solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies);
    }
    solver->SetTemperature_tr(Temperature_tr_i, Temperature_tr_i);
    solver->SetTemperature_vib(Temperature_vib_i, Temperature_vib_i);
    
    /*--- Load auxiliary vector with local adjoint variables ---*/
    Psi_i = node[iPoint]->GetSolution();		
    
    /*--- Set pressure ---*/
    solver->SetPressure(Pressure_i, Pressure_i);
    
    /*--- Axisymmetric source terms ---*/
    if (axisymmetric) {
      solver->SetJacobian_Axisymmetric(Jacobian_Axisymmetric, config);			
      for (iVar = 0; iVar < nVar; iVar ++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          Residual_Axisymmetric[iVar] += Jacobian_Axisymmetric[jVar][iVar]*Psi_i[jVar];
          Jacobian_ii[iVar][jVar] = Jacobian_Axisymmetric[jVar][iVar];
        }
      }
      node[iPoint]->AddRes_Sour(Residual_Axisymmetric);
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);	
    }
    
    /*--- Chemistry source terms ---*/
    solver->SetJacobian_Chemistry(Jacobian_Chemistry, config);
    for (iVar = 0; iVar < nVar; iVar ++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Residual_Chemistry[iVar] += Jacobian_Chemistry[jVar][iVar]*Psi_i[jVar];
        Jacobian_ii[iVar][jVar] = Jacobian_Chemistry[jVar][iVar];
      }
    }
    /*    if (iPoint == 17) {
     cout << "Psi_i: " << endl;
     for (iVar = 0; iVar < nVar; iVar++)
     cout << Psi_i[iVar] << endl;
     cout << endl << endl;
     cout << "Residual Chemistry: " << endl;
     for (iVar = 0; iVar < nVar; iVar++)
     cout << Residual_Chemistry[iVar] << endl;
     cout << endl << endl << "Jacobian Chemistry: " << endl;
     for (iVar = 0; iVar < nVar; iVar++) {
     for (jVar = 0; jVar < nVar; jVar++) {
     cout << Jacobian_Chemistry[iVar][jVar] << "\t";
     }
     cout << endl;
     }
     cin.get();
     }  */
    node[iPoint]->SubtractRes_Sour(Residual_Chemistry);
    if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);	
    
    /*--- Momentum exchange source terms ---*/
    solver->SetJacobian_MomentumExch(Jacobian_MomentumExch, config);			
    for (iVar = 0; iVar < nVar; iVar ++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Residual_MomentumExch[iVar] += Jacobian_MomentumExch[jVar][iVar]*Psi_i[jVar];
        Jacobian_ii[iVar][jVar] = Jacobian_MomentumExch[jVar][iVar];
      }
    }
    node[iPoint]->SubtractRes_Sour(Residual_MomentumExch);
    if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
    /*    if (iPoint == 17) {
     cout << "Psi_i: " << endl;
     for (iVar = 0; iVar < nVar; iVar++)
     cout << Psi_i[iVar] << endl;
     cout << endl << endl;
     cout << "Residual Momentum: " << endl;
     for (iVar = 0; iVar < nVar; iVar++)
     cout << Residual_MomentumExch[iVar] << endl;
     cout << endl << endl << "Jacobian Momentum: " << endl;
     for (iVar = 0; iVar < nVar; iVar++) {
     for (jVar = 0; jVar < nVar; jVar++) {
     cout << Jacobian_MomentumExch[iVar][jVar] << "\t";
     }
     cout << endl;
     }
     cin.get();
     }  */
    
    /*--- Energy exchange source terms ---*/
    solver->SetJacobian_EnergyExch(Jacobian_EnergyExch, config);
    for (iVar = 0; iVar < nVar; iVar ++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Residual_EnergyExch[iVar] += Jacobian_EnergyExch[jVar][iVar]*Psi_i[jVar];
        Jacobian_ii[iVar][jVar] = Jacobian_EnergyExch[jVar][iVar];
      }
    }
    
    /*    if (iPoint == 17) {
     cout << "Psi_i: " << endl;
     for (iVar = 0; iVar < nVar; iVar++)
     cout << Psi_i[iVar] << endl;
     cout << endl << endl;
     cout << "Residual Energyexch: " << endl;
     for (iVar = 0; iVar < nVar; iVar++)
     cout << Residual_Chemistry[iVar] << endl;
     cout << endl << endl << "Jacobian Energyexch: " << endl;
     for (iVar = 0; iVar < nVar; iVar++) {
     for (jVar = 0; jVar < nVar; jVar++) {
     cout << Jacobian_EnergyExch[iVar][jVar] << "\t";
     }
     cout << endl;
     }
     cin.get();
     }  */
    node[iPoint]->SubtractRes_Sour(Residual_EnergyExch);
    if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
  }
	delete[] Temperature_tr_i;
	delete[] Temperature_vib_i;
	delete[] Pressure_i;
}

void CAdjPlasmaSolution::Source_Template(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																				 CConfig *config, unsigned short iMesh) {
}


void CAdjPlasmaSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta, Delta_flow, Res, *local_ResConv, *local_ResVisc, *local_ResSour, *local_Res_TruncError, Vol;
	
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
	
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			local_Res_TruncError = node[iPoint]->GetRes_TruncError();
			local_ResConv = node[iPoint]->GetResConv();
			local_ResVisc = node[iPoint]->GetResVisc();
			local_ResSour = node[iPoint]->GetResSour();
			Vol = geometry->node[iPoint]->GetVolume();
			
			/*--- Modify matrix diagonal to assure diagonal dominance ---*/
			Delta_flow = Vol/(solution_container[PLASMA_SOL]->node[iPoint]->GetDelta_Time());
			Delta = Delta_flow;
			Jacobian.AddVal2Diag(iPoint, Delta);

			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
				Res = local_ResConv[iVar]+local_ResVisc[iVar]+local_ResSour[iVar];
				rhs[total_index] = -(Res+local_Res_TruncError[iVar]);
				AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex());
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
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->AddSolution(iVar, xsol[iPoint*nVar+iVar]);
	
  /*--- MPI solution ---*/
  SetSolution_MPI(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
	
}

void CAdjPlasmaSolution::Inviscid_Sensitivity(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config) {
  unsigned long iVertex, iPoint, Neigh;
	//unsigned short iPos, jPos;
	unsigned short iDim, iMarker, iNeigh;
  double *d = NULL, *Normal = NULL, *Psi = NULL, *U = NULL, Enthalpy, conspsi,
	Area, **PrimVar_Grad = NULL, *ConsPsi_Grad = NULL, ConsPsi, d_press, grad_v,
	v_gradconspsi;
  //double UnitaryNormal[3], *RotVel = NULL, *GridVel = NULL;
  //double Mach_Inf, Beta2;
	//double RefDensity, *RefVelocity = NULL, RefPressure;
  
	//double r, ru, rv, rw, rE, p; // used in farfield sens
	//double dp_dr, dp_dru, dp_drv, dp_drw, dp_drE; // used in farfield sens
	//double dH_dr, dH_dru, dH_drv, dH_drw, dH_drE, H; // used in farfield sens
	//double alpha, beta;
	//double *USens, *U_infty;
  
  /*--- Loop over boundary markers to select those for Euler walls ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Boundary(iMarker) == EULER_WALL)
      
    /*--- Loop over points on the surface to store the auxiliary variable ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Psi = node[iPoint]->GetSolution();
          U = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
          Enthalpy = solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy();
          
          conspsi = U[0]*Psi[0] + U[0]*Enthalpy*Psi[nDim+1];          
          for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];
          
          node[iPoint]->SetAuxVar(conspsi);
          
          /*--- Also load the auxiliary variable for first neighbors ---*/
          for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
            Neigh = geometry->node[iPoint]->GetPoint(iNeigh);
            Psi = node[Neigh]->GetSolution();
            U = solution_container[FLOW_SOL]->node[Neigh]->GetSolution();
            Enthalpy = solution_container[FLOW_SOL]->node[Neigh]->GetEnthalpy();
            conspsi = U[0]*Psi[0] + U[0]*Enthalpy*Psi[nDim+1];
            for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];
            node[Neigh]->SetAuxVar(conspsi);
          }
        }
      }
  
  /*--- Compute surface gradients of the auxiliary variable ---*/
  SetAuxVar_Surface_Gradient(geometry, config);
  
  /*--- Evaluate the shape sensitivity ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Sens_Geo[iMarker] = 0.0;
    
    if (config->GetMarker_All_Boundary(iMarker) == EULER_WALL) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          
          d = node[iPoint]->GetForceProj_Vector();
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
          
          PrimVar_Grad = solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
          ConsPsi_Grad = node[iPoint]->GetAuxVarGradient();
          ConsPsi = node[iPoint]->GetAuxVar();
          
          d_press = 0; grad_v = 0; v_gradconspsi = 0;
          for (iDim = 0; iDim < nDim; iDim++) {
            d_press += d[iDim]*PrimVar_Grad[nDim+1][iDim];
            grad_v += PrimVar_Grad[iDim+1][iDim]*ConsPsi;
//            v_gradconspsi += solution_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim, incompressible) * ConsPsi_Grad[iDim];
          }
          
          /*--- Compute sensitivity for each surface point ---*/
          CSensitivity[iMarker][iVertex] = (d_press + grad_v + v_gradconspsi) * Area;
          Sens_Geo[iMarker] -= CSensitivity[iMarker][iVertex] * Area;
        }
      }
      Total_Sens_Geo += Sens_Geo[iMarker];
    }
  }
}

void CAdjPlasmaSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	
	unsigned long iVertex, iPoint;
	double *d, *Normal, *U, *Psi_Aux, ProjVel = 0.0, bcn, vn = 0.0, Area, *UnitaryNormal, *Coord, Gamma_Minus_One;
  double *Velocity, *Psi, Enthalpy = 0.0, Energy_vib, sq_vel, phin;
  //double phis1, phis2;
	unsigned short iDim, iVar, jVar, jDim, loc, iSpecies;
  double phidotu, phidotn, Energy_el, dPdrho, *dPdrhou, dPdrhoE, dPdrhoEv;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	
	UnitaryNormal = new double[nDim];
	Velocity = new double[nDim];
	Psi      = new double[nVar];
  dPdrhou = new double[nDim];
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		if (geometry->node[iPoint]->GetDomain()) {
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Coord = geometry->node[iPoint]->GetCoord();
			
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {
				
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				
				if ( iSpecies < nDiatomics ) Gamma_Minus_One = config->GetGammaDiatomic()-1.0;
				else  Gamma_Minus_One = config->GetGammaMonatomic()-1.0;
				
				/*--- Create a copy of the adjoint solution ---*/
				Psi_Aux = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];			
				
				/*--- Flow solution ---*/
				U = solution_container[PLASMA_SOL]->node[iPoint]->GetSolution();
				d = node[iPoint]->GetForceProj_Vector();
				
				Area = 0; 
				for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
				Area = sqrt(Area);
				
				for (iDim = 0; iDim < nDim; iDim++)
					UnitaryNormal[iDim]   = -Normal[iDim]/Area;
				
				for (iDim = 0; iDim < nDim; iDim++)
					Velocity[iDim] = U[loc+iDim+1] / U[loc+0];
				
				Enthalpy = solution_container[PLASMA_SOL]->node[iPoint]->GetEnthalpy(iSpecies);
				sq_vel   = 0.5*solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity2(iSpecies);
        
        Energy_vib = 0.0;
        if (iSpecies < nDiatomics)
          Energy_vib = U[loc+nDim+2]/U[loc+0];
				
				/*--- Compute projections ---*/
				ProjVel = 0.0; bcn = 0.0; vn = 0.0, phin = 0.0; sq_vel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
//					ProjVel -= Velocity[iDim]*Normal[iDim];
					ProjVel += Velocity[iDim]*UnitaryNormal[iDim];
					bcn     += d[iDim]*UnitaryNormal[iDim];
					vn      += Velocity[iDim]*UnitaryNormal[iDim];
					phin    += Psi[loc+iDim+1]*UnitaryNormal[iDim];
				}

				/*--- Introduce the boundary condition ---*/
				for (iDim = 0; iDim < nDim; iDim++) 
					Psi[loc+iDim+1] -= ( phin - bcn ) * UnitaryNormal[iDim];

        
        // NEW IMPLEMENTATION (NOT WORKING) ------------------------
        /*--- Pre-compute useful quantities ---*/
        phidotu = 0.0;
        phidotn = 0.0;
        for (iDim = 0; iDim < nDim; iDim ++) {
          phidotu += Velocity[iDim] * Psi[loc+iDim+1];
          phidotn += Psi[loc+iDim+1] * UnitaryNormal[iDim];
        }
        Energy_el = 0;
        dPdrho = Gamma_Minus_One * (sq_vel - config->GetEnthalpy_Formation(iSpecies) - Energy_el);
        for (iDim = 0; iDim < nDim; iDim++)
          dPdrhou[iDim] = -Gamma_Minus_One*Velocity[iDim];
        dPdrhoE = Gamma_Minus_One;
        dPdrhoEv = -Gamma_Minus_One;
        
        /*--- Flux of the Euler wall: Psi^T * (dF/dU dot n) ---*/
        Residual[loc+0] = dPdrho*phidotn - ProjVel*phidotu + ProjVel*(dPdrho-Enthalpy)*Psi[loc+nDim+1];
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[loc+iDim+1] = UnitaryNormal[iDim]*Psi[loc] + ProjVel*Psi[loc+iDim+1] + phidotu*UnitaryNormal[iDim] + dPdrhou[iDim]*phidotn 
          + (dPdrhou[iDim]*ProjVel+Enthalpy*UnitaryNormal[iDim])*Psi[loc+nDim+1];
        Residual[loc+nDim+1] = dPdrhoE*phidotn + ProjVel*(1+dPdrhoE)*Psi[loc+nDim+1];
        
        if (iSpecies < nDiatomics) {
          Residual[loc+0] -= ProjVel*Energy_vib*Psi[loc+nDim+2];
          for (iDim = 0; iDim < nDim; iDim++)
            Residual[loc+iDim+1] += Energy_vib*UnitaryNormal[iDim]*Psi[loc+nDim+2];
          Residual[loc+nDim+1] += 0.0;
          Residual[loc+nDim+2] = dPdrhoEv*phidotn + ProjVel*dPdrhoEv*Psi[loc+nDim+1] + ProjVel*Psi[loc+nDim+2];
        }       
        
        /*--- Calculate Jacobians for implicit time marching ---*/
        if (implicit) {
          
          /*--- Adjoint density ---*/
          Jacobian_ii[loc+0][loc+0] = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_ii[loc+0][loc+iDim+1] = dPdrho*UnitaryNormal[iDim] - ProjVel*Velocity[iDim];
          Jacobian_ii[loc+0][loc+nDim+1] = ProjVel*(dPdrho - Enthalpy);
          
          /*--- Adjoint velocity ---*/
          for (iDim = 0; iDim < nDim; iDim++) {
            Jacobian_ii[loc+iDim+1][0] = UnitaryNormal[iDim];
            for (jDim = 0; jDim < nDim; jDim++)
              Jacobian_ii[loc+iDim+1][loc+jDim+1] = Velocity[jDim]*UnitaryNormal[iDim] + dPdrhou[iDim]*UnitaryNormal[jDim];
            Jacobian_ii[loc+iDim+1][loc+iDim+1] += ProjVel;
            Jacobian_ii[loc+iDim+1][loc+nDim+1] = dPdrhou[iDim]*ProjVel + Enthalpy;                        
          }
          
          /*--- Adjoint energy ---*/
          Jacobian_ii[loc+nDim+1][loc+0] = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) 
            Jacobian_ii[loc+nDim+1][loc+iDim+1] = dPdrhoE*UnitaryNormal[iDim];
          Jacobian_ii[loc+nDim+1][loc+nDim+1] = ProjVel*(1.0+dPdrhoE);
          
          /*--- Adjoint vibrational energy ---*/
          if (iSpecies < nDiatomics) {
            Jacobian_ii[loc+0][loc+nDim+2] = ProjVel*Energy_vib;
            for (iDim = 0; iDim < nDim; iDim++)
              Jacobian_ii[loc+iDim+1][loc+nDim+2] = Energy_vib*UnitaryNormal[iDim];
            Jacobian_ii[loc+nDim+1][loc+nDim+2] = 0.0;
            
            Jacobian_ii[loc+nDim+2][loc+0] = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Jacobian_ii[loc+nDim+2][loc+iDim+1] = dPdrhoEv*UnitaryNormal[iDim];
            Jacobian_ii[loc+nDim+2][loc+nDim+1] = ProjVel*dPdrhoEv;
            Jacobian_ii[loc+nDim+2][loc+nDim+2] = ProjVel;
          }  
        }
              
        /*--- Integrate over the area --*/
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] = Residual[iVar]*Area;
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_ii[iVar][jVar] = Jacobian_ii[iVar][jVar]*Area;          
        }        
			}			
			/*--- Update residual ---*/
			node[iPoint]->SubtractRes_Conv(Residual);
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}	
	delete [] Velocity;
	delete [] UnitaryNormal;
	delete [] Psi;
  delete [] dPdrhou;
}

void CAdjPlasmaSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
																		 CConfig *config, unsigned short val_marker) {
	
	unsigned long iVertex, iPoint;
	double *d, *Normal, *U, *Psi_Aux, ProjVel = 0.0, bcn, vn = 0.0, Area, *UnitaryNormal, *Coord, Gamma_Minus_One;
  double *Velocity, *Psi, Enthalpy = 0.0, Energy_vib, sq_vel, phin;
  //double phis1, phis2;
	unsigned short iDim, iVar, jVar, jDim, loc, iSpecies;
  double phidotu, phidotn, Energy_el, dPdrho, *dPdrhou, dPdrhoE, dPdrhoEv;
  
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	
	UnitaryNormal = new double[nDim];
	Velocity = new double[nDim];
	Psi      = new double[nVar];
  dPdrhou = new double[nDim];
  
  
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		if (geometry->node[iPoint]->GetDomain()) {
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Coord = geometry->node[iPoint]->GetCoord();
			
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ) {
				
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				
				if ( iSpecies < nDiatomics ) Gamma_Minus_One = config->GetGammaDiatomic()-1.0;
				else  Gamma_Minus_One = config->GetGammaMonatomic()-1.0;
				
				/*--- Create a copy of the adjoint solution ---*/
				Psi_Aux = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];			
				
				/*--- Flow solution ---*/
				U = solution_container[PLASMA_SOL]->node[iPoint]->GetSolution();
				d = node[iPoint]->GetForceProj_Vector();
				
				Area = 0; 
				for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
				Area = sqrt(Area);
				
				for (iDim = 0; iDim < nDim; iDim++)
					UnitaryNormal[iDim]   = -Normal[iDim]/Area;
				
				for (iDim = 0; iDim < nDim; iDim++)
					Velocity[iDim] = U[loc+iDim+1] / U[loc+0];
				
				Enthalpy = solution_container[PLASMA_SOL]->node[iPoint]->GetEnthalpy(iSpecies);
				sq_vel   = 0.5*solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity2(iSpecies);
        
        Energy_vib = 0.0;
        if (iSpecies < nDiatomics)
          Energy_vib = U[loc+nDim+2]/U[loc+0];
				
				/*--- Compute projections ---*/
				ProjVel = 0.0; bcn = 0.0; vn = 0.0, phin = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
          ProjVel += Velocity[iDim]*UnitaryNormal[iDim];
					vn      += Velocity[iDim]*UnitaryNormal[iDim];
					phin    += Psi[loc+iDim+1]*UnitaryNormal[iDim];
				}        
        
				/*--- Introduce the boundary condition ---*/
				for (iDim = 0; iDim < nDim; iDim++) 
					Psi[loc+iDim+1] -= phin * UnitaryNormal[iDim];
				
        /*--- Pre-compute useful quantities ---*/
        phidotu = 0.0;
        phidotn = 0.0;
        for (iDim = 0; iDim < nDim; iDim ++) {
          phidotu += Velocity[iDim] * Psi[loc+iDim+1];
          phidotn += Psi[loc+iDim+1] * UnitaryNormal[iDim];
        }
        Energy_el = 0;
        dPdrho = Gamma_Minus_One * (sq_vel - config->GetEnthalpy_Formation(iSpecies) - Energy_el);
        for (iDim = 0; iDim < nDim; iDim++)
          dPdrhou[iDim] = -Gamma_Minus_One*Velocity[iDim];
        dPdrhoE = Gamma_Minus_One;
        dPdrhoEv = -Gamma_Minus_One;
        
        /*--- Flux of the Euler wall: Psi^T * (dF/dU dot n) ---*/
        Residual[loc+0] = dPdrho*phidotn - ProjVel*phidotu + ProjVel*(dPdrho-Enthalpy)*Psi[loc+nDim+1];
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[loc+iDim+1] = UnitaryNormal[iDim]*Psi[loc] + ProjVel*Psi[loc+iDim+1] + phidotu*UnitaryNormal[iDim] + dPdrhou[iDim]*phidotn 
                                + (dPdrhou[iDim]*ProjVel+Enthalpy*UnitaryNormal[iDim])*Psi[loc+nDim+1];
        Residual[loc+nDim+1] = dPdrhoE*phidotn + ProjVel*(1+dPdrhoE)*Psi[loc+nDim+1];
                
        if (iSpecies < nDiatomics) {
          Residual[loc+0] -= ProjVel*Energy_vib*Psi[loc+nDim+2];
          for (iDim = 0; iDim < nDim; iDim++)
            Residual[loc+iDim+1] += Energy_vib*UnitaryNormal[iDim]*Psi[loc+nDim+2];
          Residual[loc+nDim+1] += 0.0;
          Residual[loc+nDim+2] = dPdrhoEv*phidotn + ProjVel*dPdrhoEv*Psi[loc+nDim+1] + ProjVel*Psi[loc+nDim+2];
        }        
        
        /*--- Calculate Jacobians for implicit time marching ---*/
        if (implicit) {
          
          /*--- Adjoint density ---*/
          Jacobian_ii[loc+0][loc+0] = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_ii[loc+0][loc+iDim+1] = dPdrho*UnitaryNormal[iDim] - ProjVel*Velocity[iDim];
          Jacobian_ii[loc+0][loc+nDim+1] = ProjVel*(dPdrho - Enthalpy);
          
          /*--- Adjoint velocity ---*/
          for (iDim = 0; iDim < nDim; iDim++) {
            Jacobian_ii[loc+iDim+1][0] = UnitaryNormal[iDim];
            for (jDim = 0; jDim < nDim; jDim++)
              Jacobian_ii[loc+iDim+1][loc+jDim+1] = Velocity[jDim]*UnitaryNormal[iDim] + dPdrhou[iDim]*UnitaryNormal[jDim];
            Jacobian_ii[loc+iDim+1][loc+iDim+1] += ProjVel;
            Jacobian_ii[loc+iDim+1][loc+nDim+1] = dPdrhou[iDim]*ProjVel + Enthalpy;                        
          }
          
          /*--- Adjoint energy ---*/
          Jacobian_ii[loc+nDim+1][loc+0] = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) 
            Jacobian_ii[loc+nDim+1][loc+iDim+1] = dPdrhoE*UnitaryNormal[iDim];
          Jacobian_ii[loc+nDim+1][loc+nDim+1] = ProjVel*(1.0+dPdrhoE);
                      
          /*--- Adjoint vibrational energy ---*/
          if (iSpecies < nDiatomics) {
            Jacobian_ii[loc+0][loc+nDim+2] = ProjVel*Energy_vib;
            for (iDim = 0; iDim < nDim; iDim++)
              Jacobian_ii[loc+iDim+1][loc+nDim+2] = Energy_vib*UnitaryNormal[iDim];
            Jacobian_ii[loc+nDim+1][loc+nDim+2] = 0.0;
            
            Jacobian_ii[loc+nDim+2][loc+0] = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Jacobian_ii[loc+nDim+2][loc+iDim+1] = dPdrhoEv*UnitaryNormal[iDim];
            Jacobian_ii[loc+nDim+2][loc+nDim+1] = ProjVel*dPdrhoEv;
            Jacobian_ii[loc+nDim+2][loc+nDim+2] = ProjVel;
          }  
        }
        
        /*--- Integrate over the area --*/
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] = Residual[iVar]*Area;
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_ii[iVar][jVar] = Jacobian_ii[iVar][jVar]*Area;          
        }
      }
      
      
      /*--- Update residual ---*/
      node[iPoint]->SubtractRes_Conv(Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
    }
  }
  
  delete [] Velocity;
  delete [] UnitaryNormal;
  delete [] Psi;
  delete [] dPdrhou;
}

void CAdjPlasmaSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver,
                                      CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	double *Normal, *U_domain, *U_infty, *Psi_domain, *Psi_infty;
	unsigned short loc, iSpecies;

	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	Normal = new double[nDim];
	U_domain = new double[nVar]; U_infty = new double[nVar];
	Psi_domain = new double[nVar]; Psi_infty = new double[nVar];
		
	/*--- Loop over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
			conv_solver->SetNormal(Normal);
			
			/*--- Flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(iVar);
			
			/*--- Solution at the infinity ---*/
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++) {
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				U_infty[loc + 0] = solution_container[PLASMA_SOL]->GetDensity_Inf(iSpecies);
        for (iDim = 0; iDim < nDim; iDim++)
          U_infty[loc+iDim+1] = solution_container[PLASMA_SOL]->GetDensity_Velocity_Inf(iDim, iSpecies);
				U_infty[loc+nDim+1] = solution_container[PLASMA_SOL]->GetDensity_Energy_Inf(iSpecies);
        if (iSpecies < nDiatomics)
          U_infty[loc+nDim+2] = solution_container[PLASMA_SOL]->GetDensity_Energy_vib_Inf(iSpecies);          
			}      
			
			conv_solver->SetConservative(U_domain, U_infty);

			/*--- Adjoint flow solution at the farfield ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
				Psi_infty[iVar] = 0.0;
			}
			conv_solver->SetAdjointVar(Psi_domain, Psi_infty);
			
			/*--- Compute the upwind flux ---*/
			conv_solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			/*--- Add and Subtract Residual ---*/
			node[iPoint]->SubtractRes_Conv(Residual_i);
			
			/*--- Implicit contribution to the residual ---*/
			if (implicit) 
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}
	
	delete [] Normal;
	delete [] U_domain; delete [] U_infty;
	delete [] Psi_domain; delete [] Psi_infty;
}

void CAdjPlasmaSolution::MPI_Send_Receive(CGeometry ***geometry, CSolution ****solution_container,
                                          CConfig **config, unsigned short iMGLevel, unsigned short iZone) {
	unsigned short iVar, iMarker;
	unsigned long iVertex, iPoint;
	double *Adjoint_Var; //, *Adjoint_Undivided_Laplacian = NULL;
	
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config[iZone]->GetnMarker_All(); iMarker++)
		if (config[iZone]->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			
			short SendRecv = config[iZone]->GetMarker_All_SendRecv(iMarker);
			unsigned long nVertex = geometry[iZone][iMGLevel]->nVertex[iMarker];
			
			/*--- Send information  ---*/
			if (SendRecv > 0) {
				
				/*--- Sending only occurs with MPI ---*/
#ifndef NO_MPI
				
				double **Adjoint_Grad = NULL;
				unsigned long nBuffer_Vector = geometry[iZone][iMGLevel]->nVertex[iMarker]*nVar;
				unsigned long nBuffer_Scalar = geometry[iZone][iMGLevel]->nVertex[iMarker];
				
				int send_to = SendRecv-1;
				
				/*--- Inviscid part ---*/
				double *Buffer_Send_Psi = new double[nBuffer_Vector];
				
				/*--- Upwind scheme ---*/
				if (config[iZone]->GetKind_ConvNumScheme() == SPACE_UPWIND) {
					
					double *Buffer_Send_Psix = NULL, *Buffer_Send_Psiy = NULL, *Buffer_Send_Psiz = NULL;
					Buffer_Send_Psix = new double[nBuffer_Vector];
					Buffer_Send_Psiy = new double[nBuffer_Vector];
					Buffer_Send_Psiz = new double[nBuffer_Vector];
					
					for (iVertex = 0; iVertex < nVertex; iVertex++) {
						iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();
						Adjoint_Var = node[iPoint]->GetSolution();
						Adjoint_Grad = node[iPoint]->GetGradient();
						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Send_Psi[iVar*nVertex+iVertex] = Adjoint_Var[iVar];
							Buffer_Send_Psix[iVar*nVertex+iVertex] = Adjoint_Grad[iVar][0];					
							Buffer_Send_Psiy[iVar*nVertex+iVertex] = Adjoint_Grad[iVar][1];
							if (nDim == 3) Buffer_Send_Psiz[iVar*nVertex+iVertex] = Adjoint_Grad[iVar][2];
						}
					}
					
					MPI::COMM_WORLD.Bsend(Buffer_Send_Psi,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
					MPI::COMM_WORLD.Bsend(Buffer_Send_Psix,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
					MPI::COMM_WORLD.Bsend(Buffer_Send_Psiy,nBuffer_Vector,MPI::DOUBLE,send_to, 2);
					if (nDim == 3) MPI::COMM_WORLD.Bsend(Buffer_Send_Psiz,nBuffer_Vector,MPI::DOUBLE,send_to, 3);
					
					if (iMGLevel == MESH_0) {
						delete [] Buffer_Send_Psix;
						delete [] Buffer_Send_Psiy;
						delete [] Buffer_Send_Psiz;
					}
				}
				
				/*--- Centered scheme ---*/
				if (config[iZone]->GetKind_ConvNumScheme() == SPACE_CENTERED) {
					
					double *Buffer_Send_Undivided_Laplacian = NULL, *Buffer_Send_Sensor = NULL;
					Buffer_Send_Undivided_Laplacian = new double[nBuffer_Vector];
					Buffer_Send_Sensor = new double [nBuffer_Scalar];
					
					for (iVertex = 0; iVertex < nVertex; iVertex++) {
						iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();
						Adjoint_Var = node[iPoint]->GetSolution();
	//					Adjoint_Undivided_Laplacian = node[iPoint]->GetUnd_Lapl();
						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Send_Psi[iVar*nVertex+iVertex] = Adjoint_Var[iVar];
//							Buffer_Send_Undivided_Laplacian[iVar*nVertex+iVertex] = Adjoint_Undivided_Laplacian[iVar];					
						}
//						Buffer_Send_Sensor[iVertex] = node[iPoint]->GetSensor();
					}
					
					MPI::COMM_WORLD.Bsend(Buffer_Send_Psi,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
//					MPI::COMM_WORLD.Bsend(Buffer_Send_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
//					MPI::COMM_WORLD.Bsend(Buffer_Send_Sensor, nBuffer_Scalar, MPI::DOUBLE, send_to, 2);
					
					delete [] Buffer_Send_Undivided_Laplacian;
					delete [] Buffer_Send_Sensor;
				}
				
				delete [] Buffer_Send_Psi;
				
#endif
			}
			
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
				
				double rotMatrix[3][3], *angles;
				double theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
				unsigned short iPeriodic_Index;
				
				double *newSolution = new double [nVar];
				
				int receive_from = abs(SendRecv)-1;
				unsigned long nBuffer_Vector = geometry[iZone][iMGLevel]->nVertex[iMarker]*nVar;
				unsigned long nBuffer_Scalar = geometry[iZone][iMGLevel]->nVertex[iMarker];
				
				/*--- Inviscid part ---*/
				double *Buffer_Receive_Psi = new double [nBuffer_Vector];
				
				/*--- Upwind scheme ---*/
				if (config[iZone]->GetKind_ConvNumScheme() == SPACE_UPWIND) {
					
					double *Buffer_Receive_Psix = NULL, *Buffer_Receive_Psiy = NULL, *Buffer_Receive_Psiz = NULL;
					if (iMGLevel == MESH_0) {
						Buffer_Receive_Psix = new double [nBuffer_Vector];
						Buffer_Receive_Psiy = new double [nBuffer_Vector];
						Buffer_Receive_Psiz = new double [nBuffer_Vector];
					}
					
#ifdef NO_MPI
					cout << "Upwinding for periodic boundaries in serial not yet supported." << endl;
					cout << "Press any key to exit..." << endl;
					cin.get();
					exit(0);
#else
					MPI::COMM_WORLD.Recv(Buffer_Receive_Psi,nBuffer_Vector,MPI::DOUBLE,receive_from, 0);
					if (iMGLevel == MESH_0) {
						MPI::COMM_WORLD.Recv(Buffer_Receive_Psix,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
						MPI::COMM_WORLD.Recv(Buffer_Receive_Psiy,nBuffer_Vector,MPI::DOUBLE,receive_from, 2);
						if (nDim == 3) MPI::COMM_WORLD.Recv(Buffer_Receive_Psiz,nBuffer_Vector,MPI::DOUBLE,receive_from, 3);
					}
#endif
					
					/*--- Store the received information ---*/
					for (iVertex = 0; iVertex < nVertex; iVertex++) {
						iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();				
						for (iVar = 0; iVar < nVar; iVar++) {
							node[iPoint]->SetSolution(iVar, Buffer_Receive_Psi[iVar*nVertex+iVertex]);
							if (iMGLevel == MESH_0) {
								node[iPoint]->SetGradient(iVar, 0, Buffer_Receive_Psix[iVar*nVertex+iVertex]);
								node[iPoint]->SetGradient(iVar, 1, Buffer_Receive_Psiy[iVar*nVertex+iVertex]);
								if (nDim == 3) node[iPoint]->SetGradient(iVar, 2, Buffer_Receive_Psiz[iVar*nVertex+iVertex]);
							}
						}
					}
					
					if (iMGLevel == MESH_0) {
						delete [] Buffer_Receive_Psix;
						delete [] Buffer_Receive_Psiy;
						delete [] Buffer_Receive_Psiz;
					}
					
				}
				
				/*--- Centered scheme ---*/
				if (config[iZone]->GetKind_ConvNumScheme() == SPACE_CENTERED) {
					
					double *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Receive_Sensor = NULL;
					if (iMGLevel == MESH_0) {
						Buffer_Receive_Undivided_Laplacian = new double [nBuffer_Vector];
						Buffer_Receive_Sensor = new double [nBuffer_Scalar];
					}			
					
#ifdef NO_MPI
					/*--- Allow for periodic boundaries to use SEND_RECEIVE in serial.
					 Serial computations will only use the BC in receive mode, as
					 the proc will always be sending information to itself. ---*/
					
					/*--- Retrieve the donor information from the matching marker ---*/
					unsigned short donor_marker= 1;
					unsigned long donorPoint;
					for (unsigned short iMark = 0; iMark < config[iZone]->GetnMarker_All(); iMark++){
						if (config[iZone]->GetMarker_All_SendRecv(iMark)-1 == receive_from) donor_marker = iMark;
					}
					
					/*--- Get the information from the donor directly. This is a serial
					 computation with access to all nodes. Note that there is an
					 implicit ordering in the list. ---*/
					for (iVertex = 0; iVertex < geometry[iZone][iMGLevel]->nVertex[iMarker]; iVertex++) {
						
						/*--- Get the specific donor point for this vertex ---*/
						donorPoint = geometry[iZone][iMGLevel]->vertex[donor_marker][iVertex]->GetNode();
						
						/*--- Get the solution from the donor and store it in the buffer.
						 This is essentially a dummy send/receive. ---*/
						Adjoint_Var = node[donorPoint]->GetSolution();
//						Adjoint_Undivided_Laplacian = node[donorPoint]->GetUnd_Lapl();
						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Receive_Psi[iVar*nVertex+iVertex] = Adjoint_Var[iVar];
//							Buffer_Receive_Undivided_Laplacian[iVar*nVertex+iVertex] = Adjoint_Undivided_Laplacian[iVar];					
						}
//						Buffer_Receive_Sensor[iVertex] = node[donorPoint]->GetSensor();
					}
#else
					MPI::COMM_WORLD.Recv(Buffer_Receive_Psi,nBuffer_Vector,MPI::DOUBLE,receive_from, 0);
//					MPI::COMM_WORLD.Recv(Buffer_Receive_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
//					MPI::COMM_WORLD.Recv(Buffer_Receive_Sensor,nBuffer_Scalar,MPI::DOUBLE,receive_from, 2);
#endif
					
					for (iVertex = 0; iVertex < geometry[iZone][iMGLevel]->nVertex[iMarker]; iVertex++) {
						
						/*--- Get the current point and its Send/Receive type. ---*/
						iPoint = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetNode();
						iPeriodic_Index = geometry[iZone][iMGLevel]->vertex[iMarker][iVertex]->GetRotation_Type();
						
						/*--- Retrieve the supplied periodic information. ---*/
						angles = config[iZone]->GetPeriodicRotation(iPeriodic_Index);
						
						/*--- Store angles separately for clarity. ---*/
						theta    = angles[0];   phi    = angles[1]; psi    = angles[2];
						cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
						sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
						
						/*--- Compute the rotation matrix. Note that the implicit
						 ordering is rotation about the x-axis, y-axis,
						 then z-axis. Note that this is the transpose of the matrix
						 used during the preprocessing stage. ---*/
						rotMatrix[0][0] = cosPhi*cosPsi; rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi; rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
						rotMatrix[0][1] = cosPhi*sinPsi; rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi; rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
						rotMatrix[0][2] = -sinPhi; rotMatrix[1][2] = sinTheta*cosPhi; rotMatrix[2][2] = cosTheta*cosPhi;
						
						/*--- Copy solution before performing transformation. ---*/
						for (iVar = 0; iVar < nVar; iVar++) newSolution[iVar] = Buffer_Receive_Psi[iVar*nVertex+iVertex];
						
						/*--- Rotate the adjoint velocity components. ---*/
						if (nDim == 2) {
							newSolution[1] = rotMatrix[0][0]*Buffer_Receive_Psi[1*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_Psi[2*nVertex+iVertex];
							newSolution[2] = rotMatrix[1][0]*Buffer_Receive_Psi[1*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_Psi[2*nVertex+iVertex];
						} 
						else {
							newSolution[1] = rotMatrix[0][0]*Buffer_Receive_Psi[1*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_Psi[2*nVertex+iVertex] + rotMatrix[0][2]*Buffer_Receive_Psi[3*nVertex+iVertex];
							newSolution[2] = rotMatrix[1][0]*Buffer_Receive_Psi[1*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_Psi[2*nVertex+iVertex] + rotMatrix[1][2]*Buffer_Receive_Psi[3*nVertex+iVertex];
							newSolution[3] = rotMatrix[2][0]*Buffer_Receive_Psi[1*nVertex+iVertex] + rotMatrix[2][1]*Buffer_Receive_Psi[2*nVertex+iVertex] + rotMatrix[2][2]*Buffer_Receive_Psi[3*nVertex+iVertex];
						}
						
						/*--- Copy transformed solution back into buffer. ---*/
						for (iVar = 0; iVar < nVar; iVar++) Buffer_Receive_Psi[iVar*nVertex+iVertex] = newSolution[iVar];
						
						/*--- Set the solution for this point ---*/
						for (iVar = 0; iVar < nVar; iVar++) {
							node[iPoint]->SetSolution(iVar, Buffer_Receive_Psi[iVar*nVertex+iVertex]);
//							node[iPoint]->SetUndivided_Laplacian(iVar, Buffer_Receive_Undivided_Laplacian[iVar*nVertex+iVertex]);
						}
//						node[iPoint]->SetSensor(Buffer_Receive_Sensor[iVertex]);
						
					}
					delete [] Buffer_Receive_Undivided_Laplacian;
					delete [] Buffer_Receive_Sensor;
				}
				delete [] Buffer_Receive_Psi;
				delete [] newSolution;
			}
		}
}
