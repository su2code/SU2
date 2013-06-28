/*!
 * \file solution_adjoint_plasma.cpp
 * \brief Main subrotuines for solving adjoint problems (Euler, Navier-Stokes, etc.).
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

CAdjPlasmaSolution::CAdjPlasmaSolution(void) : CSolution() { }

CAdjPlasmaSolution::CAdjPlasmaSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned long iPoint, iVertex;
	string text_line, mesh_filename;
	unsigned short iDim, iVar, iMarker;
	ifstream restart_file;
	string filename, AdjExt;

	bool restart = config->GetRestart();

	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	nMonatomics = 7;
	nDiatomics = 0;
	nSpecies = nMonatomics + nDiatomics;
	nVar = nMonatomics*(nDim+2) + nDiatomics*(nDim+3);	
	node = new CVariable*[geometry->GetnPoint()];

	/*--- Define some auxiliary vectors related to the residual ---*/
	Residual = new double[nVar];	  Residual_Max = new double[nVar];
	Residual_i = new double[nVar];	Residual_j = new double[nVar];
	Res_Conv_i = new double[nVar];	Res_Visc_i = new double[nVar];
	Res_Conv_j = new double[nVar];	Res_Visc_j = new double[nVar];
	Res_Sour_i = new double[nVar];	Res_Sour_j = new double[nVar];

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
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_ii[iVar] = new double [nVar];
			Jacobian_ij[iVar] = new double [nVar];
			Jacobian_ji[iVar] = new double [nVar];
			Jacobian_jj[iVar] = new double [nVar];
		}
		Initialize_Jacobian_Structure(geometry, config);
		xsol = new double [geometry->GetnPoint()*nVar];
		rhs  = new double [geometry->GetnPoint()*nVar];
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
	
	/*--- Sensitivity definition and coefficient in all the markers ---*/
	CSensitivity = new double* [config->GetnMarker_All()];
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		CSensitivity[iMarker] = new double [geometry->nVertex[iMarker]];
	}
	CSens_Geo  = new double[config->GetnMarker_All()];
	CSens_Mach = new double[config->GetnMarker_All()];
	CSens_AoA  = new double[config->GetnMarker_All()];
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		CSens_Geo[iMarker]  = 0.0;
		CSens_Mach[iMarker] = 0.0;
		CSens_AoA[iMarker]  = 0.0;
		for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
			CSensitivity[iMarker][iVertex] = 0.0;
	}

	/*--- Adjoint flow at the inifinity, initialization stuff ---*/
	PsiRho_Inf = 0.0; PsiE_Inf   = 0.0;
	Phi_Inf    = new double [nDim];
	Phi_Inf[0] = 0.0; Phi_Inf[1] = 0.0;
	if (nDim == 3) Phi_Inf[2] = 0.0;		
  
	if (!restart) {
		/*--- Restart the solution from infinity ---*/
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CAdjPlasmaVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, nDim, nVar, config);
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
      case FORCE_X_COEFFICIENT: AdjExt = "_cfx.dat"; break;
      case FORCE_Y_COEFFICIENT: AdjExt = "_cfy.dat"; break;
      case FORCE_Z_COEFFICIENT: AdjExt = "_cfz.dat"; break;
		}
		filename.append(AdjExt);
		restart_file.open(filename.data(), ios::in);

		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
			cout << "There is no adjoint restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get(); exit(1);
		}

		/*--- Read the restart file ---*/
/*		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
			if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
			node[iPoint] = new CAdjPlasmaVariable(Solution, nDim, nVar, config);
		}
		restart_file.close();*/
	}
  
  /*--- Define solver parameters needed for execution of destructor ---*/
	if (config->GetKind_ConvNumScheme_Plasma() == SPACE_CENTRED) space_centered = true;
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
	delete [] CSens_Geo; delete [] CSens_Mach;
	delete [] CSens_AoA; delete [] Phi_Inf;
	
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

void CAdjPlasmaSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);  
  
	/*--- Residual initialization ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		
		/*--- Initialize the convective residual vector ---*/
		node[iPoint]->Set_ResConv_Zero();
		
		/*--- Initialize the source residual vector ---*/
		node[iPoint]->Set_ResSour_Zero();
		
		/*--- Initialize the viscous residual vector ---*/
		if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit)
			node[iPoint]->Set_ResVisc_Zero();   
	}
	
	/*--- Implicit solution ---*/
	if (implicit) Jacobian.SetValZero();
	
}

void CAdjPlasmaSolution::SetForceProj_Vector(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	double *ForceProj_Vector, x = 0.0, y = 0.0, z = 0.0, *Normal, C_d, C_l, C_t, C_q;
  double x_origin, y_origin, z_origin, Area;
	unsigned short iMarker, iDim;
	unsigned long iVertex, iPoint;
	double Alpha      = (config->GetAoA()*PI_NUMBER)/180.0;
	double Beta       = (config->GetAoS()*PI_NUMBER)/180.0;
	double RefAreaCoeff    = config->GetRefAreaCoeff();
	double RefLengthMoment  = config->GetRefLengthMoment();
	double *RefOriginMoment = config->GetRefOriginMoment();
  double Scale_GradOF = config->GetScale_GradOF();
  double RefVel2, RefDensity;
  bool rotating_frame = config->GetRotating_Frame();
  
  /*--- If this is a rotating frame problem, use the specified speed
        for computing the force coefficients. Otherwise, use the 
        freestream values which is the standard convention. ---*/
  if (rotating_frame) {
    /*--- Use the rotational origin for computing moments ---*/
    RefOriginMoment = config->GetRotAxisOrigin();
    /*--- Reference area is the rotor disk area ---*/
    RefLengthMoment = config->GetRotRadius();
    RefAreaCoeff = PI_NUMBER*RefLengthMoment*RefLengthMoment;
    /* --- Reference velocity is the rotational speed times rotor radius ---*/
		RefVel2     = (config->GetOmegaMag()*RefLengthMoment)*(config->GetOmegaMag()*RefLengthMoment);
		RefDensity  = config->GetDensity_FreeStreamND();
  } else {
    double *Velocity_Inf = config->GetVelocity_FreeStreamND();
    RefVel2 = 0.0; 
    for (iDim = 0; iDim < nDim; iDim++) 
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    RefDensity  = config->GetDensity_FreeStreamND();
  }
  
/*--- In parallel computations the Cd, and Cl must be recomputed using all the processors ---*/
#ifdef NO_MPI
	C_d = solution_container[PLASMA_SOL]->GetTotal_CDrag();
	C_l = solution_container[PLASMA_SOL]->GetTotal_CLift();
  C_t = solution_container[PLASMA_SOL]->GetTotal_CT();
  C_q = solution_container[PLASMA_SOL]->GetTotal_CQ();
#else
	double *sbuf_force = new double[4];
	double *rbuf_force = new double[4];
	sbuf_force[0] = solution_container[PLASMA_SOL]->GetTotal_CDrag();
	sbuf_force[1] = solution_container[PLASMA_SOL]->GetTotal_CLift();
  sbuf_force[2] = solution_container[PLASMA_SOL]->GetTotal_CT();
	sbuf_force[3] = solution_container[PLASMA_SOL]->GetTotal_CQ();
	MPI::COMM_WORLD.Reduce(sbuf_force, rbuf_force, 4, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
	MPI::COMM_WORLD.Bcast(rbuf_force, 4, MPI::DOUBLE, MASTER_NODE);
	C_d = rbuf_force[0];
	C_l = rbuf_force[1];
  C_t = rbuf_force[2];
  C_q = rbuf_force[3];
	delete [] sbuf_force;
	delete [] rbuf_force;
#endif
	
  /*--- Compute coefficients needed for objective function evaluation. ---*/
	C_d += config->GetCteViscDrag();
	double C_p    = 1.0/(0.5*RefDensity*RefAreaCoeff*RefVel2);
	double invCD  = 1.0 / C_d;
	double CLCD2  = C_l / (C_d*C_d);

	ForceProj_Vector = new double [nDim];
	x_origin = RefOriginMoment[0]; y_origin = RefOriginMoment[1]; z_origin = RefOriginMoment[2];

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
						if (nDim == 2) { ForceProj_Vector[0] = C_p*cos(Alpha); ForceProj_Vector[1] = C_p*sin(Alpha); }
						if (nDim == 3) { ForceProj_Vector[0] = C_p*cos(Alpha)*cos(Beta); ForceProj_Vector[1] = C_p*sin(Beta); ForceProj_Vector[2] = C_p*sin(Alpha)*cos(Beta); }
						break;
					case LIFT_COEFFICIENT :
						if (nDim == 2) { ForceProj_Vector[0] = -C_p*sin(Alpha); ForceProj_Vector[1] = C_p*cos(Alpha); }
						if (nDim == 3) { ForceProj_Vector[0] = -C_p*sin(Alpha); ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = C_p*cos(Alpha); }
						break;
					case SIDEFORCE_COEFFICIENT :
						if (nDim == 2) { cout << "This functional is not possible in 2D!!" << endl;
							cout << "Press any key to exit..." << endl; cin.get(); exit(1);
						}
						if (nDim == 3) { ForceProj_Vector[0] = -C_p*sin(Beta) * cos(Alpha); ForceProj_Vector[1] = C_p*cos(Beta); ForceProj_Vector[2] = -C_p*sin(Beta) * sin(Alpha); }
						break;
					case PRESSURE_COEFFICIENT :
						if (nDim == 2) {
							Area = sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1]);
							ForceProj_Vector[0] = -C_p*Normal[0]/Area; ForceProj_Vector[1] = -C_p*Normal[1]/Area; 
						}
						if (nDim == 3) {
							Area = sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);
							ForceProj_Vector[0] = -C_p*Normal[0]/Area; ForceProj_Vector[1] = -C_p*Normal[1]/Area; ForceProj_Vector[2] = -C_p*Normal[2]/Area; 
						}
						break;
					case MOMENT_X_COEFFICIENT :
						if (nDim == 2) { cout << "This functional is not possible in 2D!!" << endl;
							cout << "Press any key to exit..." << endl; cin.get(); exit(1);
						}
						if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = -C_p*(z - z_origin)/RefLengthMoment; ForceProj_Vector[2] = C_p*(y - y_origin)/RefLengthMoment; }
						break;
					case MOMENT_Y_COEFFICIENT :
						if (nDim == 2) { cout << "This functional is not possible in 2D!!" << endl;
							cout << "Press any key to exit..." << endl;
							cin.get(); exit(1);
						}
						if (nDim == 3) { ForceProj_Vector[0] = -C_p*(z - z_origin)/RefLengthMoment; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = C_p*(x - x_origin)/RefLengthMoment; }
						break;
					case MOMENT_Z_COEFFICIENT :
						if (nDim == 2) { ForceProj_Vector[0] = -C_p*(y - y_origin)/RefLengthMoment; ForceProj_Vector[1] = C_p*(x - x_origin)/RefLengthMoment; }
						if (nDim == 3) { ForceProj_Vector[0] = -C_p*(y - y_origin)/RefLengthMoment; ForceProj_Vector[1] = C_p*(x - x_origin)/RefLengthMoment; ForceProj_Vector[2] = 0; }
						break;
					case EFFICIENCY :
						if (nDim == 2) { ForceProj_Vector[0] = -C_p*(invCD*sin(Alpha)+CLCD2*cos(Alpha)); 
							ForceProj_Vector[1] = C_p*(invCD*cos(Alpha)-CLCD2*sin(Alpha)); }
						if (nDim == 3) { ForceProj_Vector[0] = -C_p*(invCD*sin(Alpha)+CLCD2*cos(Alpha)*cos(Beta)); 
							ForceProj_Vector[1] = -C_p*CLCD2*sin(Beta); 
							ForceProj_Vector[2] = C_p*(invCD*cos(Alpha)-CLCD2*sin(Alpha)*cos(Beta)); }
						break;	
          case FORCE_X_COEFFICIENT :
						if (nDim == 2) { ForceProj_Vector[0] = C_p; ForceProj_Vector[1] = 0.0; }
						if (nDim == 3) { ForceProj_Vector[0] = C_p; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 0.0; }
						break;
          case FORCE_Y_COEFFICIENT :
						if (nDim == 2) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = C_p; }
						if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = C_p; ForceProj_Vector[2] = 0.0; }
						break;
          case FORCE_Z_COEFFICIENT :
						if (nDim == 2) {cout << "This functional is not possible in 2D!!" << endl;
							cout << "Press any key to exit..." << endl;
							cin.get(); exit(1);
            }
						if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = C_p; }
						break;
				}
				
        /*--- Scale the gradient of the objective function, if specified ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          ForceProj_Vector[iDim] *= Scale_GradOF;
        }
        
        /*--- Store the force projection vector at this node ---*/
				node[iPoint]->SetForceProj_Vector(ForceProj_Vector);
			}
	
	delete [] ForceProj_Vector;
}

void CAdjPlasmaSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {

	double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, *Limiter_i = NULL, *Limiter_j = NULL, *Psi_i, *Psi_j, *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	bool high_order_diss = ((config->GetKind_Upwind() == ROE_2ND) && (iMesh == MESH_0));
  
	if (high_order_diss) { 
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry); 
		if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) || 
			(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) SetSolution_Gradient_LS(geometry, config);
    if (config->GetKind_SlopeLimit() != NONE) SetSolution_Limiter(geometry, config);
	}

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
		
		solver->SetSoundSpeed(solution_container[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(), 
													solution_container[PLASMA_SOL]->node[jPoint]->GetSoundSpeed());		
		solver->SetEnthalpy(solution_container[PLASMA_SOL]->node[iPoint]->GetEnthalpy(), 
												solution_container[PLASMA_SOL]->node[jPoint]->GetEnthalpy());
    		
		/*--- High order reconstruction using MUSCL strategy ---*/
		if (high_order_diss) { 
			for (iDim = 0; iDim < nDim; iDim++) {
				Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
				Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}
			
			Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
			if (config->GetKind_SlopeLimit() != NONE) {
				Limiter_i = node[iPoint]->GetLimiter(); Limiter_j = node[jPoint]->GetLimiter();
			}
			
			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
        if (config->GetKind_SlopeLimit() == NONE) {
          Solution_i[iVar] = Psi_i[iVar] + Project_Grad_i;
          Solution_j[iVar] = Psi_j[iVar] + Project_Grad_j;
        } else {
          Solution_i[iVar] = Psi_i[iVar] + Project_Grad_i*Limiter_i[iDim];
          Solution_j[iVar] = Psi_j[iVar] + Project_Grad_j*Limiter_j[iDim];
        }
			}
			/*--- Set conservative variables with reconstruction ---*/
			solver->SetAdjointVar(Solution_i, Solution_j);
		}
		
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

void CAdjPlasmaSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																								 CConfig *config, unsigned short iMesh) {
	
}

void CAdjPlasmaSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta, Delta_flow, Res, *local_ResConv, *local_ResVisc, *local_ResSour, *local_TruncationError, Vol;
	
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);
	
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			local_TruncationError = node[iPoint]->GetTruncationError();
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
				rhs[total_index] = -(Res+local_TruncationError[iVar]);
				xsol[total_index] = 0.0;
				AddRes_Max( iVar, Res*Res*Vol );
			}
		}

	/*--- Solve the system ---*/
//	Jacobian.SGSSolution(rhs, xsol, 1e-9, 100, true, geometry, config);
	Jacobian.LU_SGSIteration(rhs, xsol, geometry, config);
	
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) 
			for (iVar = 0; iVar < nVar; iVar++)
				node[iPoint]->AddSolution(iVar, xsol[iPoint*nVar+iVar]);
	
#ifdef NO_MPI
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, sqrt(GetRes_Max(iVar)));
#endif
	
}

void CAdjPlasmaSolution::Inviscid_Sensitivity(CGeometry *geometry, CSolution **solution_container, CConfig *config) {

}

void CAdjPlasmaSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	
	unsigned long iVertex, iPoint;
	double *d, *Normal, *U, *Psi_Aux, ProjVel = 0.0, bcn, vn = 0.0, Area, *UnitaryNormal, *Coord, Gamma_Minus_One;
  double *Velocity, *Psi, Enthalpy = 0.0, sq_vel, phin, phis1, phis2;
	unsigned short iDim, iVar, jDim, loc, iSpecies;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);
	
	UnitaryNormal = new double[nDim];
	Velocity = new double[nDim];
	Psi      = new double[nVar];
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
					Velocity[iDim] = U[iDim+1] / U[0];
				
				Enthalpy = solution_container[PLASMA_SOL]->node[iPoint]->GetEnthalpy();
				sq_vel   = 0.5*solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity2();
				
				/*--- Compute projections ---*/
				ProjVel = 0.0; bcn = 0.0; vn = 0.0, phin = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					ProjVel -= Velocity[iDim]*Normal[iDim];
					bcn     += d[iDim]*UnitaryNormal[iDim];
					vn      += Velocity[iDim]*UnitaryNormal[iDim];
					phin    += Psi[loc+iDim+1]*UnitaryNormal[iDim];
				}
				
				/*--- Introduce the boundary condition ---*/
				for (iDim = 0; iDim < nDim; iDim++) 
					Psi[loc+iDim+1] -= ( phin - bcn ) * UnitaryNormal[iDim];
				
				/*--- Inner products after introducing BC (Psi has changed) ---*/
				phis1 = 0.0; phis2 = Psi[0] + Enthalpy * Psi[loc+nDim+1];
				for (iDim = 0; iDim < nDim; iDim++) {
					phis1 -= Normal[iDim]*Psi[loc+iDim+1];
					phis2 += Velocity[iDim]*Psi[loc+iDim+1];
				}
				
				/*--- Flux of the Euler wall ---*/
				Residual[loc+0] = ProjVel * Psi[0] - phis2 * ProjVel + phis1 * Gamma_Minus_One * sq_vel;
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[loc+iDim+1] = ProjVel * Psi[iDim+1] - phis2 * Normal[iDim] - phis1 * Gamma_Minus_One * Velocity[iDim];
				Residual[loc+nDim+1] = ProjVel * Psi[nDim+1] + phis1 * Gamma_Minus_One;
				
			}
			
			/*--- Update residual ---*/
			node[iPoint]->SubtractRes_Conv(Residual);
			
			/*--- Implicit stuff ---*/
			if (implicit) {
				
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++ ) {
					
					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
					
					/*--- Adjoint density ---*/
					Jacobian_ii[loc + 0][loc + 0] = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_ii[loc + 0][loc + iDim+1] = -ProjVel * (Velocity[iDim] - UnitaryNormal[iDim] * vn);
					Jacobian_ii[loc + 0][loc + nDim+1] = -ProjVel * Enthalpy;
					
					/*--- Adjoint velocities ---*/
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_ii[loc + iDim+1][loc + 0] = -Normal[iDim];
						for (jDim = 0; jDim < nDim; jDim++)
							Jacobian_ii[loc + iDim+1][loc + jDim+1] = -ProjVel*(UnitaryNormal[jDim]*UnitaryNormal[iDim] - Normal[iDim] * (Velocity[jDim] - UnitaryNormal[jDim] * vn));
						Jacobian_ii[loc + iDim+1][loc + iDim+1] += ProjVel;
						Jacobian_ii[loc + iDim+1][loc + nDim+1] = -Normal[iDim] * Enthalpy;
					}
					
					/*--- Adjoint energy ---*/
					Jacobian_ii[loc + nDim+1][loc + 0] = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_ii[loc + nDim+1][loc + iDim+1] = 0.0;
					Jacobian_ii[loc + nDim+1][loc + nDim+1] = ProjVel;
					
					Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
				}
			}
		}
	}
	
	delete [] Velocity;
	delete [] UnitaryNormal;
	delete [] Psi;
}

void CAdjPlasmaSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
																		 CConfig *config, unsigned short val_marker) {
}

void CAdjPlasmaSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
									 CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	double *Normal, *U_domain, *U_infty, *Psi_domain, *Psi_infty;
	
	bool implicit = (config->GetKind_TimeIntScheme_Plasma() == EULER_IMPLICIT);

	Normal = new double[nDim];
	U_domain = new double[nVar]; U_infty = new double[nVar];
	Psi_domain = new double[nVar]; Psi_infty = new double[nVar];
		
	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
			solver->SetNormal(Normal);
			
			/*--- Flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(iVar);
			
			/*--- Solution at the infinity ---*/
			U_infty[0] = solution_container[PLASMA_SOL]->GetDensity_Inf();
			U_infty[1] = solution_container[PLASMA_SOL]->GetDensity_Velocity_Inf(0);
			U_infty[2] = solution_container[PLASMA_SOL]->GetDensity_Velocity_Inf(1);
			U_infty[3] = solution_container[PLASMA_SOL]->GetDensity_Energy_Inf();
			if (nDim == 3) {
				U_infty[3] = solution_container[PLASMA_SOL]->GetDensity_Velocity_Inf(2);
				U_infty[4] = solution_container[PLASMA_SOL]->GetDensity_Energy_Inf();
			}
			solver->SetConservative(U_domain, U_infty);

			/*--- Adjoint flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
				Psi_infty[iVar] = 0.0;
			}
			solver->SetAdjointVar(Psi_domain, Psi_infty);

				solver->SetSoundSpeed(solution_container[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(), 
															solution_container[PLASMA_SOL]->node[iPoint]->GetSoundSpeed());		
				solver->SetEnthalpy(solution_container[PLASMA_SOL]->node[iPoint]->GetEnthalpy(), 
														solution_container[PLASMA_SOL]->node[iPoint]->GetEnthalpy());
						
			/*--- Compute the upwind flux ---*/
			solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

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

void CAdjPlasmaSolution::MPI_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
																				 unsigned short val_mesh) {
	unsigned short iVar, iMarker;
	unsigned long iVertex, iPoint;
	double *Adjoint_Var, *Adjoint_Undivided_Laplacian = NULL;
	
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			
			short SendRecv = config->GetMarker_All_SendRecv(iMarker);
			unsigned long nVertex = geometry->nVertex[iMarker];
			
			/*--- Send information  ---*/
			if (SendRecv > 0) {
				
				/*--- Sending only occurs with MPI ---*/
#ifndef NO_MPI
				
				double **Adjoint_Grad = NULL;
				unsigned long nBuffer_Vector = geometry->nVertex[iMarker]*nVar;
				unsigned long nBuffer_Scalar = geometry->nVertex[iMarker];
				
				int send_to = SendRecv-1;
				
				/*--- Inviscid part ---*/
				double *Buffer_Send_Psi = new double[nBuffer_Vector];
				
				/*--- Upwind scheme ---*/
				if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {
					
					double *Buffer_Send_Psix = NULL, *Buffer_Send_Psiy = NULL, *Buffer_Send_Psiz = NULL;
					if (val_mesh == MESH_0) {
						Buffer_Send_Psix = new double[nBuffer_Vector];
						Buffer_Send_Psiy = new double[nBuffer_Vector];
						Buffer_Send_Psiz = new double[nBuffer_Vector];
					}
					
					for (iVertex = 0; iVertex < nVertex; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						Adjoint_Var = node[iPoint]->GetSolution();
						if (val_mesh == MESH_0) 
							Adjoint_Grad = node[iPoint]->GetGradient();
						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Send_Psi[iVar*nVertex+iVertex] = Adjoint_Var[iVar];
							if (val_mesh == MESH_0) {
								Buffer_Send_Psix[iVar*nVertex+iVertex] = Adjoint_Grad[iVar][0];					
								Buffer_Send_Psiy[iVar*nVertex+iVertex] = Adjoint_Grad[iVar][1];
								if (nDim == 3) Buffer_Send_Psiz[iVar*nVertex+iVertex] = Adjoint_Grad[iVar][2];
							}
						}
					}
					
					MPI::COMM_WORLD.Bsend(Buffer_Send_Psi,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
					if (val_mesh == MESH_0) {
						MPI::COMM_WORLD.Bsend(Buffer_Send_Psix,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
						MPI::COMM_WORLD.Bsend(Buffer_Send_Psiy,nBuffer_Vector,MPI::DOUBLE,send_to, 2);
						if (nDim == 3) MPI::COMM_WORLD.Bsend(Buffer_Send_Psiz,nBuffer_Vector,MPI::DOUBLE,send_to, 3);
					}
					
					if (val_mesh == MESH_0) {
						delete [] Buffer_Send_Psix;
						delete [] Buffer_Send_Psiy;
						delete [] Buffer_Send_Psiz;
					}
				}
				
				/*--- Centered scheme ---*/
				if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {
					
					double *Buffer_Send_Undivided_Laplacian = NULL, *Buffer_Send_Sensor = NULL;
					if (val_mesh == MESH_0) {
						Buffer_Send_Undivided_Laplacian = new double[nBuffer_Vector];
						Buffer_Send_Sensor = new double [nBuffer_Scalar];
					}
					
					for (iVertex = 0; iVertex < nVertex; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						Adjoint_Var = node[iPoint]->GetSolution();
						if (val_mesh == MESH_0) Adjoint_Undivided_Laplacian = node[iPoint]->GetUnd_Lapl();
						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Send_Psi[iVar*nVertex+iVertex] = Adjoint_Var[iVar];
							if (val_mesh == MESH_0) Buffer_Send_Undivided_Laplacian[iVar*nVertex+iVertex] = Adjoint_Undivided_Laplacian[iVar];					
						}
						if (val_mesh == MESH_0) Buffer_Send_Sensor[iVertex] = node[iPoint]->GetSensor();
					}
					
					MPI::COMM_WORLD.Bsend(Buffer_Send_Psi,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
					if (val_mesh == MESH_0) {
						MPI::COMM_WORLD.Bsend(Buffer_Send_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
						MPI::COMM_WORLD.Bsend(Buffer_Send_Sensor, nBuffer_Scalar, MPI::DOUBLE, send_to, 2);
					}
					
					if (val_mesh == MESH_0) {
						delete [] Buffer_Send_Undivided_Laplacian;
						delete [] Buffer_Send_Sensor;
					}
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
				unsigned long nBuffer_Vector = geometry->nVertex[iMarker]*nVar;
				unsigned long nBuffer_Scalar = geometry->nVertex[iMarker];
				
				/*--- Inviscid part ---*/
				double *Buffer_Receive_Psi = new double [nBuffer_Vector];
				
				/*--- Upwind scheme ---*/
				if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {
					
					double *Buffer_Receive_Psix = NULL, *Buffer_Receive_Psiy = NULL, *Buffer_Receive_Psiz = NULL;
					if (val_mesh == MESH_0) {
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
					if (val_mesh == MESH_0) {
						MPI::COMM_WORLD.Recv(Buffer_Receive_Psix,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
						MPI::COMM_WORLD.Recv(Buffer_Receive_Psiy,nBuffer_Vector,MPI::DOUBLE,receive_from, 2);
						if (nDim == 3) MPI::COMM_WORLD.Recv(Buffer_Receive_Psiz,nBuffer_Vector,MPI::DOUBLE,receive_from, 3);
					}
#endif
					
					/*--- Store the received information ---*/
					for (iVertex = 0; iVertex < nVertex; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();				
						for (iVar = 0; iVar < nVar; iVar++) {
							node[iPoint]->SetSolution(iVar, Buffer_Receive_Psi[iVar*nVertex+iVertex]);
							if (val_mesh == MESH_0) {
								node[iPoint]->SetGradient(iVar, 0, Buffer_Receive_Psix[iVar*nVertex+iVertex]);
								node[iPoint]->SetGradient(iVar, 1, Buffer_Receive_Psiy[iVar*nVertex+iVertex]);
								if (nDim == 3) node[iPoint]->SetGradient(iVar, 2, Buffer_Receive_Psiz[iVar*nVertex+iVertex]);
							}
						}
					}
					
					if (val_mesh == MESH_0) {
						delete [] Buffer_Receive_Psix;
						delete [] Buffer_Receive_Psiy;
						delete [] Buffer_Receive_Psiz;
					}
					
				}
				
				/*--- Centered scheme ---*/
				if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {
					
					double *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Receive_Sensor = NULL;
					if (val_mesh == MESH_0) {
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
					for (unsigned short iMark = 0; iMark < config->GetnMarker_All(); iMark++){
						if (config->GetMarker_All_SendRecv(iMark)-1 == receive_from) donor_marker = iMark;
					}
					
					/*--- Get the information from the donor directly. This is a serial
					 computation with access to all nodes. Note that there is an
					 implicit ordering in the list. ---*/
					for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						
						/*--- Get the specific donor point for this vertex ---*/
						donorPoint = geometry->vertex[donor_marker][iVertex]->GetNode();
						
						/*--- Get the solution from the donor and store it in the buffer.
						 This is essentially a dummy send/receive. ---*/
						Adjoint_Var = node[donorPoint]->GetSolution();
						if (val_mesh == MESH_0) Adjoint_Undivided_Laplacian = node[donorPoint]->GetUnd_Lapl();
						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Receive_Psi[iVar*nVertex+iVertex] = Adjoint_Var[iVar];
							if (val_mesh == MESH_0) Buffer_Receive_Undivided_Laplacian[iVar*nVertex+iVertex] = Adjoint_Undivided_Laplacian[iVar];					
						}
						if (val_mesh == MESH_0) Buffer_Receive_Sensor[iVertex] = node[donorPoint]->GetSensor();
					}
#else
					MPI::COMM_WORLD.Recv(Buffer_Receive_Psi,nBuffer_Vector,MPI::DOUBLE,receive_from, 0);
					if (val_mesh == MESH_0) {
						MPI::COMM_WORLD.Recv(Buffer_Receive_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
						MPI::COMM_WORLD.Recv(Buffer_Receive_Sensor,nBuffer_Scalar,MPI::DOUBLE,receive_from, 2);
					}
#endif
					
					for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						
						/*--- Get the current point and its Send/Receive type. ---*/
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						iPeriodic_Index = geometry->vertex[iMarker][iVertex]->GetRotation_Type();
						
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
							if (val_mesh == MESH_0) node[iPoint]->SetUndivided_Laplacian(iVar, Buffer_Receive_Undivided_Laplacian[iVar*nVertex+iVertex]);
						}
						if (val_mesh == MESH_0) node[iPoint]->SetSensor(Buffer_Receive_Sensor[iVertex]);
						
					}
					if (val_mesh == MESH_0) {
						delete [] Buffer_Receive_Undivided_Laplacian;
						delete [] Buffer_Receive_Sensor;
					}
				}
				delete [] Buffer_Receive_Psi;
				delete [] newSolution;
			}
		}
}
