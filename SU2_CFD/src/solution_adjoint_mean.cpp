/*!
 * \file solution_adjoint_mean.cpp
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

CAdjEulerSolution::CAdjEulerSolution(void) : CSolution() { }

CAdjEulerSolution::CAdjEulerSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned long iPoint, index, iVertex;
	string text_line, mesh_filename;
	unsigned short iDim, iVar, iMarker;
	ifstream restart_file;
	string filename, AdjExt;

	bool restart = config->GetRestart();
	bool incompressible = config->GetIncompressible();

	/*--- Set the gamma value ---*/
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	if (incompressible) nVar = nDim + 1;
	else nVar = nDim + 2;
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

  /*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED) {
		p1_Und_Lapl = new double [geometry->GetnPoint()]; 
		p2_Und_Lapl = new double [geometry->GetnPoint()]; 
  }
  
	/*--- Define some auxiliary vectors related to the geometry ---*/
//  Vector = new double[nDim];
	Vector_i = new double[nDim]; Vector_j = new double[nDim];

	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT) {
		
		/*--- Point to point Jacobians ---*/
		Jacobian_i = new double* [nVar]; 
		Jacobian_j = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new double [nVar]; 
			Jacobian_j[iVar] = new double [nVar];
		}
		
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
			node[iPoint] = new CAdjEulerVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, nDim, nVar, config);
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
			case ELECTRIC_CHARGE: AdjExt = "_ec.dat"; break;
			case EQUIVALENT_AREA: AdjExt = "_ea.dat"; break;
			case NEARFIELD_PRESSURE: AdjExt = "_nfp.dat"; break;
      case FORCE_X_COEFFICIENT: AdjExt = "_cfx.dat"; break;
      case FORCE_Y_COEFFICIENT: AdjExt = "_cfy.dat"; break;
      case FORCE_Z_COEFFICIENT: AdjExt = "_cfz.dat"; break;
      case THRUST_COEFFICIENT: AdjExt = "_ct.dat"; break;
      case TORQUE_COEFFICIENT: AdjExt = "_cq.dat"; break;
      case FIGURE_OF_MERIT: AdjExt = "_merit.dat"; break;
      case FREESURFACE: AdjExt = "_fs.dat"; break;
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
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
			if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
			node[iPoint] = new CAdjEulerVariable(Solution, nDim, nVar, config);
		}
		restart_file.close();
	}
  
  /*--- Define solver parameters needed for execution of destructor ---*/
	if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTRED) space_centered = true;
	else space_centered = false;
  
}

CAdjEulerSolution::~CAdjEulerSolution(void) {
	unsigned short iVar, iDim;

	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Jacobian_ii[iVar]; delete [] Jacobian_ij[iVar];
		delete [] Jacobian_ji[iVar]; delete [] Jacobian_jj[iVar];
	}
	delete [] Jacobian_ii; delete [] Jacobian_ij;
	delete [] Jacobian_ji; delete [] Jacobian_jj;
	
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Jacobian_i[iVar];
		delete [] Jacobian_j[iVar];
	}
	delete [] Jacobian_i;
	delete [] Jacobian_j;

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

void CAdjEulerSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);  
  
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

void CAdjEulerSolution::SetForceProj_Vector(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	double *ForceProj_Vector, x = 0.0, y = 0.0, z = 0.0, *Normal, C_d, C_l, C_t, C_q;
  double x_origin, y_origin, z_origin, WDrag, Area;
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
	C_d = solution_container[FLOW_SOL]->GetTotal_CDrag();
	C_l = solution_container[FLOW_SOL]->GetTotal_CLift();
  C_t = solution_container[FLOW_SOL]->GetTotal_CT();
  C_q = solution_container[FLOW_SOL]->GetTotal_CQ();
#else
	double *sbuf_force = new double[4];
	double *rbuf_force = new double[4];
	sbuf_force[0] = solution_container[FLOW_SOL]->GetTotal_CDrag();
	sbuf_force[1] = solution_container[FLOW_SOL]->GetTotal_CLift();
  sbuf_force[2] = solution_container[FLOW_SOL]->GetTotal_CT();
	sbuf_force[3] = solution_container[FLOW_SOL]->GetTotal_CQ();
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
  double invCQ  = 1.0/C_q;
  double CTRCQ2 = C_t/(RefLengthMoment*C_q*C_q);

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
					case EQUIVALENT_AREA :
						WDrag = config->GetWeightCd();
						if (nDim == 2) { ForceProj_Vector[0] = C_p*cos(Alpha)*WDrag; 
							ForceProj_Vector[1] = C_p*sin(Alpha)*WDrag; }
						if (nDim == 3) { ForceProj_Vector[0] = C_p*cos(Alpha)*cos(Beta)*WDrag; 
							ForceProj_Vector[1] = C_p*sin(Beta)*WDrag; 
							ForceProj_Vector[2] = C_p*sin(Alpha)*cos(Beta)*WDrag; }
						break;	
					case NEARFIELD_PRESSURE :
						WDrag = config->GetWeightCd();
						if (nDim == 2) { ForceProj_Vector[0] = C_p*cos(Alpha)*WDrag; 
							ForceProj_Vector[1] = C_p*sin(Alpha)*WDrag; }
						if (nDim == 3) { ForceProj_Vector[0] = C_p*cos(Alpha)*cos(Beta)*WDrag; 
							ForceProj_Vector[1] = C_p*sin(Beta)*WDrag; 
							ForceProj_Vector[2] = C_p*sin(Alpha)*cos(Beta)*WDrag; }				
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
          case THRUST_COEFFICIENT :
						if (nDim == 2) { ForceProj_Vector[0] = -C_p; ForceProj_Vector[1] = 0.0; }
						if (nDim == 3) { ForceProj_Vector[0] = -C_p; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 0.0; }
						break;
          case TORQUE_COEFFICIENT :
						if (nDim == 2) {cout << "This functional is not possible in 2D!!" << endl;
							cout << "Press any key to exit..." << endl;
							cin.get(); exit(1);
            }						
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = C_p*(z - z_origin)/RefLengthMoment; ForceProj_Vector[2] = -C_p*(y - y_origin)/RefLengthMoment; }						
            break;
          case FIGURE_OF_MERIT :
						if (nDim == 2) {cout << "This functional is not possible in 2D!!" << endl;
							cout << "Press any key to exit..." << endl;
							cin.get(); exit(1);
            }		
            if (nDim == 3) { 
              ForceProj_Vector[0] = -C_p*invCQ; 
							ForceProj_Vector[1] = -C_p*CTRCQ2*(z - z_origin); 
							ForceProj_Vector[2] =  C_p*CTRCQ2*(y - y_origin); 
            }
						break;
					case FREESURFACE :
						if (nDim == 2) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; }
						if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 0.0; }
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

void CAdjEulerSolution::SetIntBoundary_Jump(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iMarker, iVar, jVar, kVar, jc, jrjc, jrjcm1, jrjcp1, jr, jm, jrm1, jrjr, jrp1, jmjm;
	unsigned long iVertex, iPoint, iPointNearField, nPointNearField = 0;
	double aux, *IntBound_Vector, *coord, u, v, w = 0.0, sq_vel, E = 0.0, *U_i, A[5][5], M[5][5], AM[5][5], b[5], WeightSB, sum, MinDist, 
	CoordNF_aux, rho, NearFieldWeight_aux, *CoordNF = NULL, *NearFieldWeight = NULL, Xcoord, Dist, DerivativeOF = 0.0;
	ifstream index_file;
	string text_line;
	double factor = 1.0;
	
	IntBound_Vector = new double [nVar];
	
	/*--- If equivalent area objective function, read the value of 
	 the derivative from a file, this is a preprocess of the direct solution ---*/ 
	
	if (config->GetKind_ObjFunc() == EQUIVALENT_AREA) {
		
		/*--- Read derivative of the objective function at the NearField from file ---*/
		index_file.open("WeightNF.dat", ios::in);
		if (index_file.fail()) {
			cout << "There is no Weight Nearfield Pressure file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}
		
		/*--- Dimensionalization bucle ---*/
		nPointNearField = 0;
		while (!index_file.eof()) {
			getline(index_file, text_line);
			istringstream value_line(text_line);
			value_line >> CoordNF_aux >> NearFieldWeight_aux;
			nPointNearField++;
		}
		index_file.close();
		
		CoordNF = new double [nPointNearField];
		NearFieldWeight = new double [nPointNearField];
		
		/*--- Store the information ---*/
		index_file.open("WeightNF.dat", ios::in);
		nPointNearField = 0;
		while (!index_file.eof()) {
			getline(index_file, text_line);
			istringstream value_line(text_line);
			value_line >> CoordNF_aux >> NearFieldWeight_aux;
			CoordNF[nPointNearField] = CoordNF_aux;
			NearFieldWeight[nPointNearField] = NearFieldWeight_aux*config->GetScale_GradOF();
			nPointNearField++;
		}
		index_file.close();
		
	}
	
	/*--- Compute the jump on the adjoint variables for the upper and the lower side ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				
				if (geometry->node[iPoint]->GetDomain()) {
					
					coord = geometry->node[iPoint]->GetCoord();
					WeightSB = 1.0-config->GetWeightCd(); DerivativeOF = 0.0;
					
					switch (config->GetKind_ObjFunc()) {	
						case EQUIVALENT_AREA : 							
							MinDist = 1E6; 
							for (iPointNearField = 0; iPointNearField < nPointNearField; iPointNearField++) {
								Xcoord = CoordNF[iPointNearField]; Dist = abs(Xcoord-coord[0]);
								if (Dist <= MinDist) {
									MinDist = Dist;
									DerivativeOF = factor*WeightSB*NearFieldWeight[iPointNearField];
								}
							}
							/*--- This is very important!! ---*/
							if ((MinDist > 1E-6) || (coord[nDim-1] > 0.0)) DerivativeOF = 0.0;
							
							break;	
							
						case NEARFIELD_PRESSURE :
							DerivativeOF = factor*WeightSB*(solution_container[FLOW_SOL]->node[iPoint]->GetPressure() 
														- solution_container[FLOW_SOL]->GetPressure_Inf());
							break;
					}
					
					/*--- Compute the jump of the adjoint variables --*/
					U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
					u = U_i[1]/U_i[0]; v = U_i[2]/U_i[0];
					rho = U_i[0];
					if (nDim == 2)	E = U_i[3]/U_i[0];
					if (nDim == 3) { w = U_i[3]/U_i[0]; E = U_i[4]/U_i[0]; }

					if (nDim == 2) {
						sq_vel = u*u+v*v;
						A[0][0] = 0.0;																					A[0][1] = 0.0;									A[0][2] = 1.0;																				A[0][3] = 0.0;
						A[1][0] = -u*v;																					A[1][1] = v;										A[1][2] = u;																					A[1][3] = 0.0;
						A[2][0] = 0.5*(Gamma-3.0)*v*v+0.5*Gamma_Minus_One*u*u;	A[2][1] = -Gamma_Minus_One*u;		A[2][2] = (3.0-Gamma)*v;															A[2][3] = Gamma_Minus_One;
						A[3][0] = -Gamma*v*E+Gamma_Minus_One*v*sq_vel;					A[3][1] = -Gamma_Minus_One*u*v; A[3][2] = Gamma*E-0.5*Gamma_Minus_One*(u*u+3.0*v*v);	A[3][3] = Gamma*v;
							
						M[0][0] = 1.0;				M[0][1] = 0.0;		M[0][2] = 0.0;		M[0][3] = 0.0;
						M[1][0] = u;					M[1][1] = rho;		M[1][2] = 0.0;		M[1][3] = 0.0;
						M[2][0] = v;					M[2][1] = 0.0;		M[2][2] = rho;		M[2][3] = 0.0;
						M[3][0] = 0.5*sq_vel;	M[3][1] = rho*u;	M[3][2] = rho*v;	M[3][3] = 1.0/Gamma_Minus_One;
						
						for (iVar = 0; iVar < nVar; iVar++)
							for (jVar = 0; jVar < nVar; jVar++) {
								aux = 0.0;
								for (kVar = 0; kVar < nVar; kVar++)
									aux += A[iVar][kVar]*M[kVar][jVar];
								AM[iVar][jVar] = aux;
							}
						
						for (iVar = 0; iVar < nVar; iVar++)
							for (jVar = 0; jVar < nVar; jVar++)
								A[iVar][jVar] = AM[jVar][iVar];
												
						b[0] = 0.0;  
						b[1] = 0.0; 
						b[2] = 0.0; 
						b[3] = DerivativeOF; 
					}
					
					if (nDim == 3) {
						sq_vel = u*u+v*v+w*w;
						A[0][0] = 0.0;																	A[0][1] = 0.0;									A[0][2] = 0.0;									A[0][3] = 1.0;																				A[0][4] = 0.0;
						A[1][0] = -u*w;																	A[1][1] = w;										A[1][2] = 0.0;									A[1][3] = u;																					A[1][4] = 0.0;
						A[2][0] = -v*w;																	A[2][1] = 0.0;									A[2][2] = w;										A[2][3] = v;																					A[2][4] = 0.0;
						A[3][0] = -w*w+0.5*(Gamma_Minus_One)*sq_vel;		A[3][1] = -Gamma_Minus_One*u;		A[3][2] = -Gamma_Minus_One*v;		A[3][3] = (3.0-Gamma)*w;															A[3][4] = Gamma_Minus_One;
						A[4][0] = -w*(Gamma*E-Gamma_Minus_One*sq_vel);	A[4][1] = -Gamma_Minus_One*u*w; A[4][2] = -Gamma_Minus_One*v*w; A[4][3] = Gamma*E-0.5*Gamma_Minus_One*(sq_vel+2*w*w); A[4][4] = Gamma*w;
						
						M[0][0] = 1.0;				M[0][1] = 0.0;		M[0][2] = 0.0;		M[0][3] = 0.0;		M[0][4] = 0.0;
						M[1][0] = u;					M[1][1] = rho;		M[1][2] = 0.0;		M[1][3] = 0.0;		M[1][4] = 0.0;
						M[2][0] = v;					M[2][1] = 0.0;		M[2][2] = rho;		M[2][3] = 0.0;		M[2][4] = 0.0;
						M[3][0] = w;					M[3][1] = 0.0;		M[3][2] = 0.0;		M[3][3] = rho;		M[3][4] = 0.0;
						M[4][0] = 0.5*sq_vel;	M[4][1] = rho*u;	M[4][2] = rho*v;	M[4][3] = rho*w;	M[4][4] = 1.0/Gamma_Minus_One;

						for (iVar = 0; iVar < nVar; iVar++)
							for (jVar = 0; jVar < nVar; jVar++) {
								aux = 0.0;
								for (kVar = 0; kVar < nVar; kVar++)
									aux += A[iVar][kVar]*M[kVar][jVar];
								AM[iVar][jVar] = aux;
							}
						
						for (iVar = 0; iVar < nVar; iVar++)
							for (jVar = 0; jVar < nVar; jVar++)
								A[iVar][jVar] = AM[jVar][iVar];
						
						b[0] = 0.0; 
						b[1] = 0.0; 
						b[2] = 0.0;
						b[3] = 0.0; 
						b[4] = DerivativeOF;
					}
					
					/*--- Solve the system using a LU descomposition --*/
					for (jc = 1; jc < nVar; jc++)
						A[0][jc] /= A[0][0];
					
					jrjc = 0;						
					for (;;) {
						jrjc++; jrjcm1 = jrjc-1; jrjcp1 = jrjc+1;
						for (jr = jrjc; jr < nVar; jr++) {
							sum = A[jr][jrjc];
							for (jm = 0; jm <= jrjcm1; jm++)
								sum -= A[jr][jm]*A[jm][jrjc];
							A[jr][jrjc] = sum;
						}
						if ( jrjc == (nVar-1) ) goto stop;
						for (jc = jrjcp1; jc<nVar; jc++) {
							sum = A[jrjc][jc];
							for (jm = 0; jm <= jrjcm1; jm++)
								sum -= A[jrjc][jm]*A[jm][jc];
							A[jrjc][jc] = sum/A[jrjc][jrjc];
						}
					}
					
				stop: 
					
					b[0] = b[0]/A[0][0];
					for (jr = 1; jr<nVar; jr++) {
						jrm1 = jr-1;
						sum = b[jr];
						for (jm = 0; jm<=jrm1; jm++)
							sum -= A[jr][jm]*b[jm];
						b[jr] = sum/A[jr][jr];
					}
					
					for (jrjr = 1; jrjr<nVar; jrjr++) {
						jr = (nVar-1)-jrjr;
						jrp1 = jr+1;
						sum = b[jr];
						for (jmjm = jrp1; jmjm<nVar; jmjm++) {
							jm = (nVar-1)-jmjm+jrp1;
							sum -= A[jr][jm]*b[jm];
						}
						b[jr] = sum;
					}
					
					for (iVar = 0; iVar < nVar; iVar++)
						IntBound_Vector[iVar] = b[iVar];
					
					node[iPoint]->SetIntBoundary_Jump(IntBound_Vector);
				}
			}
	
	delete [] IntBound_Vector;
	
	if (config->GetKind_ObjFunc() == EQUIVALENT_AREA) {
		delete [] CoordNF;
		delete [] NearFieldWeight;
	}
	
}

void CAdjEulerSolution::Centred_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
											  CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	
	unsigned long iEdge, iPoint, jPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool dissipation = ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit);
	bool high_order_diss = ((config->GetKind_Centred() == JST) && (iMesh == MESH_0));		
  bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();

	if (dissipation && high_order_diss) {
		SetDissipation_Switch(geometry, solution_container);
		SetUndivided_Laplacian(geometry, config);
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry); 
		if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) || 
				(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) SetSolution_Gradient_LS(geometry, config);
  }
  
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Points in edge, normal, and neighbors---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());
		solver->SetNeighbor(geometry->node[iPoint]->GetnPoint(), geometry->node[jPoint]->GetnPoint());
	
		/*--- Adjoint variables w/o reconstruction ---*/
		solver->SetAdjointVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
		
		/*--- Conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), 
								solution_container[FLOW_SOL]->node[jPoint]->GetSolution());
		
		if (incompressible) {
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solution_container[FLOW_SOL]->node[jPoint]->GetDensityInc());
			solver->SetBetaInc2(solution_container[FLOW_SOL]->node[iPoint]->GetBetaInc2(), solution_container[FLOW_SOL]->node[jPoint]->GetBetaInc2());
			solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
		}
		else {
			solver->SetSoundSpeed(solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
														solution_container[FLOW_SOL]->node[jPoint]->GetSoundSpeed());
			solver->SetEnthalpy(solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
													solution_container[FLOW_SOL]->node[jPoint]->GetEnthalpy());
		}
		
		solver->SetLambda(solution_container[FLOW_SOL]->node[iPoint]->GetLambda(), 
											solution_container[FLOW_SOL]->node[jPoint]->GetLambda());
		
		if (dissipation && high_order_diss) {
			solver->SetUndivided_Laplacian(node[iPoint]->GetUnd_Lapl(),node[jPoint]->GetUnd_Lapl());
      solver->SetSensor(solution_container[FLOW_SOL]->node[iPoint]->GetSensor(),
                        solution_container[FLOW_SOL]->node[jPoint]->GetSensor());
    }
    
    /*--- Rotational frame - use rotating sensor ---*/
		if (rotating_frame) {
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());
      solver->SetRotFlux(geometry->edge[iEdge]->GetRotFlux());
      solver->SetSensor(node[iPoint]->GetSensor(),node[jPoint]->GetSensor());
    }
    
		/*--- Compute residuals ---*/				
		solver->SetResidual(Res_Conv_i, Res_Visc_i, Res_Conv_j, Res_Visc_j, 
												Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, dissipation, config);
		
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


void CAdjEulerSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {

	double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, *Limiter_i = NULL, *Limiter_j = NULL, *Psi_i, *Psi_j, *U_i, *U_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;
	
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool high_order_diss = ((config->GetKind_Upwind() == ROE_2ND) && (iMesh == MESH_0));
	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
  
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
		U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
		U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
		solver->SetConservative(U_i, U_j);
		
		if (incompressible) {
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solution_container[FLOW_SOL]->node[jPoint]->GetDensityInc());
			solver->SetBetaInc2(solution_container[FLOW_SOL]->node[iPoint]->GetBetaInc2(), solution_container[FLOW_SOL]->node[jPoint]->GetBetaInc2());
			solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
		}
		else {		
			solver->SetSoundSpeed(solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
														solution_container[FLOW_SOL]->node[jPoint]->GetSoundSpeed());		
			solver->SetEnthalpy(solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
													solution_container[FLOW_SOL]->node[jPoint]->GetEnthalpy());
		}
    
    /*--- Rotational frame ---*/
		if (rotating_frame) {
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());
      solver->SetRotFlux(geometry->edge[iEdge]->GetRotFlux());
    }
		
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

void CAdjEulerSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																								 CConfig *config, unsigned short iMesh) {
	  
  unsigned short iVar;
	unsigned long iPoint;

	bool rotating_frame = config->GetRotating_Frame();
	bool axisymmetric = config->GetAxisymmetric();
	bool gravity = (config->GetGravityForce() == YES);
	bool freesurface = (config->GetFreeSurface() == YES);

	for (iVar = 0; iVar < nVar; iVar++)
		Residual[iVar] = 0;

	if (rotating_frame) {
		
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
	
	if (axisymmetric) {

	}

	if (gravity) {
		
	}
	
	if (freesurface) {
		if (config->GetKind_ObjFunc() == FREESURFACE) {	
			for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 
				
				double Volume = geometry->node[iPoint]->GetVolume();
				double **Gradient = solution_container[ADJLEVELSET_SOL]->node[iPoint]->GetGradient();
				double coeff = solution_container[LEVELSET_SOL]->node[iPoint]->GetSolution()[0] 
				/ solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
				
				for (iVar = 1; iVar < nVar; iVar++) {
					Residual[iVar] = coeff*Gradient[0][iVar-1]*Volume;
				}
				node[iPoint]->AddRes_Conv(Residual);
				
			}
		}		
	}
  
}

void CAdjEulerSolution::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
	
  unsigned long iPoint, jPoint, iEdge, Point_Normal = 0, iVertex;
	unsigned short iVar, iMarker;
	double *Diff = new double[nVar];
  double *U_halo = new double[nVar];
  double *Normal;
	
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetUnd_LaplZero();
  
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
    
		for (iVar = 0; iVar < nVar; iVar++)
			Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);
    
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
	}
	
  /*--- Loop over all boundaries and include an extra contribution
	 from a halo node. Find the nearest normal, interior point 
	 for a boundary node and make a linear approximation. ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
				config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
				config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY )
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					
					Point_Normal = geometry->vertex[iMarker][iVertex]->GetClosest_Neighbor();

					/*--- Interpolate & compute difference in the conserved variables ---*/
					for (iVar = 0; iVar < nVar; iVar++) {
						U_halo[iVar] = 2.0*node[iPoint]->GetSolution(iVar) - node[Point_Normal]->GetSolution(iVar);
						Diff[iVar]   = node[iPoint]->GetSolution(iVar) - U_halo[iVar];
					}
					
					/*--- Subtract contribution at the boundary node only ---*/
					if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
				}
			}
  
	delete [] Diff;
  delete [] U_halo;
}

void CAdjEulerSolution::SetDissipation_Switch(CGeometry *geometry, CSolution **solution_container) {
	
	double dx = 0.1;
	double LimK = 0.03;
	double eps2 =  pow((LimK*dx),3);
	
	unsigned long iPoint, jPoint;
	unsigned short iNeigh, nNeigh, iDim;
	double **Gradient_i, *Coord_i, *Coord_j, diff_coord, dist_ij, r_u, r_u_ij, 
	du_max, du_min, u_ij, *Solution_i, *Solution_j, dp, dm;
	
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) 
		
		if (geometry->node[iPoint]->GetDomain()) {
			
			Solution_i = node[iPoint]->GetSolution();
			Gradient_i = node[iPoint]->GetGradient();
			Coord_i = geometry->node[iPoint]->GetCoord();
			nNeigh = geometry->node[iPoint]->GetnPoint();
			
			/*--- Find max and min value of the variable in the control volume around the mesh point ---*/
			du_max = 1.0E-8; du_min = -1.0E-8;
			for (iNeigh = 0; iNeigh < nNeigh; iNeigh++) {
				jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
				Solution_j = node[jPoint]->GetSolution();
				du_max = max(du_max, Solution_j[0] - Solution_i[0]);
				du_min = min(du_min, Solution_j[0] - Solution_i[0]);
			}
			
			r_u = 1.0;
			for (iNeigh = 0; iNeigh < nNeigh; iNeigh++) {
				
				/*--- Unconstrained reconstructed solution ---*/
				jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
				Solution_j = node[jPoint]->GetSolution();
				Coord_j = geometry->node[jPoint]->GetCoord();
				u_ij = Solution_i[0]; dist_ij = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					diff_coord = Coord_j[iDim]-Coord_i[iDim];
					u_ij += 0.5*diff_coord*Gradient_i[0][iDim];
				}
				
				/*--- Venkatakrishnan limiter ---*/
				if ((u_ij - Solution_i[0]) >= 0.0) dp = du_max;
				else	dp = du_min;
				dm = u_ij - Solution_i[0];
				r_u_ij = (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2);
				
				/*--- Take the smallest value of the limiter ---*/
				r_u = min(r_u, r_u_ij);
				
			}
			node[iPoint]->SetSensor(1.0-r_u);
		}
}


void CAdjEulerSolution::ExplicitRK_Iteration(CGeometry *geometry, CSolution **solution_container, 
	CConfig *config, unsigned short iRKStep) {
	double *Residual, *TruncationError, Vol, Delta;
	unsigned short iVar;
	unsigned long iPoint;
	double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);

	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, 0.0 );

	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			Vol = geometry->node[iPoint]->GetVolume();
			Delta = solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() / Vol;
			TruncationError = node[iPoint]->GetTruncationError();
			Residual = node[iPoint]->GetResidual();
			for (iVar = 0; iVar < nVar; iVar++) {
				node[iPoint]->AddSolution(iVar, -(Residual[iVar]+TruncationError[iVar])*Delta*RK_AlphaCoeff);
				AddRes_Max( iVar, Residual[iVar]*Residual[iVar]*Vol );
			}
		}
	
#ifdef NO_MPI
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, sqrt(GetRes_Max(iVar)) );
#endif
	
}

void CAdjEulerSolution::ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	double *Residual, *TruncationError, Vol, Delta;
	unsigned short iVar;
	unsigned long iPoint;
	
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, 0.0 );
	
	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
			Vol = geometry->node[iPoint]->GetVolume();
			Delta = solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() / Vol;
			TruncationError = node[iPoint]->GetTruncationError();
			Residual = node[iPoint]->GetResidual();
			for (iVar = 0; iVar < nVar; iVar++) {
				node[iPoint]->AddSolution(iVar, -(Residual[iVar]+TruncationError[iVar])*Delta);
				AddRes_Max( iVar, Residual[iVar]*Residual[iVar]*Vol );
			}
		}
	
#ifdef NO_MPI
	/*--- Compute the norm-2 of the residual ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max( iVar, sqrt(GetRes_Max(iVar)) );
#endif
	
}

void CAdjEulerSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta, Res, *local_ResConv, *local_ResVisc, *local_ResSour, *local_TruncationError, Vol;
	
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		SetRes_Max(iVar, 0.0);
	
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		/*--- Read the residual ---*/
		local_TruncationError = node[iPoint]->GetTruncationError();
		local_ResConv = node[iPoint]->GetResConv();
		local_ResVisc = node[iPoint]->GetResVisc();
		local_ResSour = node[iPoint]->GetResSour();
		
		/*--- Read the volume ---*/
		Vol = geometry->node[iPoint]->GetVolume();
		
		/*--- Modify matrix diagonal to assure diagonal dominance ---*/
		Delta = Vol/(solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
		
		Jacobian.AddVal2Diag(iPoint, Delta);
		
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			Res = local_ResConv[iVar]+local_ResVisc[iVar]+local_ResSour[iVar];
			rhs[total_index] = -(Res+local_TruncationError[iVar]);
			xsol[total_index] = 0.0;
			AddRes_Max( iVar, Res*Res*Vol );
		}
	}
	
	/*--- Solve the system ---*/
	if (config->GetKind_Linear_Solver() == SYM_GAUSS_SEIDEL) Jacobian.SGSSolution(rhs, xsol, 1e-9, 100, true, geometry, config);
	if (config->GetKind_Linear_Solver() == LU_SGS) Jacobian.LU_SGSIteration(rhs, xsol, geometry, config);
	
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

void CAdjEulerSolution::Inviscid_Sensitivity(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned long iVertex, iPoint, Neigh;
	unsigned short iDim, iMarker, iNeigh;
	double *ForceProj_Vector = NULL, *Normal = NULL, *Psi = NULL, *U = NULL, Enthalpy, conspsi, Mach_Inf,
	ProjPhi, Area, **PrimVar_Grad = NULL, *ConsPsi_Grad = NULL, ConsPsi, d_press, grad_v, Beta2,
	RefVelocity, v_gradconspsi, ProjVel, UnitaryNormal[3], Mach_Inf_3,
	Psi_E, *RotVel = NULL, RefDensity, RefPressure;
	
  bool rotating_frame = config->GetRotating_Frame();
  bool incompressible = config->GetIncompressible();

  /*--- Initialize snesitivities to zero ---*/
	Total_CSens_Geo = 0.0; Total_CSens_Mach = 0.0; Total_CSens_AoA = 0.0;

	/*--- Loop over boundary markers to select those for Euler walls ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == EULER_WALL)
			
		/*--- Loop over points on the surface to store the auxiliar variable ---*/
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					Psi = node[iPoint]->GetSolution();
					U = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
					if (incompressible) {
						Beta2 = solution_container[FLOW_SOL]->node[iPoint]->GetBetaInc2();
						conspsi = Beta2*Psi[0];
					}
					else {
						Enthalpy = solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy();
						conspsi = U[0]*Psi[0] + U[0]*Enthalpy*Psi[nDim+1];
					}
					for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];
					
					node[iPoint]->SetAuxVar(conspsi);
					
					/*--- Also load the auxiliar variable for first neighbors ---*/
					for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
						Neigh = geometry->node[iPoint]->GetPoint(iNeigh);
						Psi = node[Neigh]->GetSolution();
						U = solution_container[FLOW_SOL]->node[Neigh]->GetSolution();
						if (incompressible) {
							Beta2 = solution_container[FLOW_SOL]->node[Neigh]->GetBetaInc2();
							conspsi = Beta2*Psi[0];
						}
						else {							
							Enthalpy = solution_container[FLOW_SOL]->node[Neigh]->GetEnthalpy();
							conspsi = U[0]*Psi[0] + U[0]*Enthalpy*Psi[nDim+1];
						}
						for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];
						node[Neigh]->SetAuxVar(conspsi);
					}
				}
			}

	/*--- Compute surface gradients of the auxiliar variable ---*/
	SetAuxVar_Surface_Gradient(geometry, config);

	/*--- Evaluate the shape sensitivity ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		CSens_Geo[iMarker] = 0.0;
		
		if (config->GetMarker_All_Boundary(iMarker) == EULER_WALL) {
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					
					ForceProj_Vector = node[iPoint]->GetForceProj_Vector();
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					Area = 0; for (iDim = 0; iDim < nDim; iDim++)
						Area += Normal[iDim]*Normal[iDim];
					Area = sqrt(Area);
					
					PrimVar_Grad = solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
					ConsPsi_Grad = node[iPoint]->GetAuxVarGradient();
					ConsPsi = node[iPoint]->GetAuxVar();
					
          /*--- If rotating, make adjustment ---*/
          if (rotating_frame) RotVel = geometry->node[iPoint]->GetRotVel();
          
					d_press = 0; grad_v = 0; v_gradconspsi = 0;
					for (iDim = 0; iDim < nDim; iDim++) {
						d_press += ForceProj_Vector[iDim]*PrimVar_Grad[nDim+1][iDim];
						grad_v += PrimVar_Grad[iDim+1][iDim]*ConsPsi;
						v_gradconspsi += solution_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim) * ConsPsi_Grad[iDim];
						if (rotating_frame) v_gradconspsi -= RotVel[iDim] * ConsPsi_Grad[iDim];
					}
					
					/*--- Compute sensitivity for each surface point ---*/
					CSensitivity[iMarker][iVertex] = (d_press + grad_v + v_gradconspsi) * Area;
					CSens_Geo[iMarker] -= CSensitivity[iMarker][iVertex] * Area;
				}
			}
			Total_CSens_Geo += CSens_Geo[iMarker];
		}
	}
		
	/*--- Mach number sensitivity of the objective function ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == FAR_FIELD) {
			CSens_Mach[iMarker] = 0.0;
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					Psi = node[iPoint]->GetSolution();
					U = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					
          /*--- Use the defined reference values if Mach_Inf = 0.0 ---*/
          Mach_Inf   = config->GetMach_FreeStreamND();
          Mach_Inf_3 = Mach_Inf*Mach_Inf*Mach_Inf;
          if (Mach_Inf == 0.0) {
            if (rotating_frame) {
              RefVelocity = config->GetRotRadius()*config->GetOmegaMag();
              RefPressure = config->GetPressure_FreeStreamND();
              RefDensity  = config->GetDensity_FreeStreamND();
              Mach_Inf    = RefVelocity/sqrt(Gamma*RefPressure/RefDensity);
              Mach_Inf_3  = Mach_Inf*Mach_Inf*Mach_Inf;
            } else {
            RefVelocity = config->GetVelocity_Ref();
            RefPressure = config->GetPressure_Ref();
            RefDensity  = config->GetDensity_Ref();
            Mach_Inf    = RefVelocity/sqrt(Gamma*RefPressure/RefDensity);
            Mach_Inf_3  = Mach_Inf*Mach_Inf*Mach_Inf;
            }
          }
          
					Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) 
						Area += Normal[iDim]*Normal[iDim];
					Area = sqrt(Area);
					
					for (iDim = 0; iDim < nDim; iDim++)
						UnitaryNormal[iDim] = -Normal[iDim]/Area;
					
					ProjPhi = 0; ProjVel = 0;
					for (iDim = 0; iDim < nDim; iDim++) {
						ProjPhi += Psi[iDim+1]*UnitaryNormal[iDim];
						ProjVel += U[iDim+1]*UnitaryNormal[iDim]/U[0];
					}
					Psi_E = Psi[nVar-1];
					CSens_Mach[iMarker] += (2.0/Mach_Inf_3)*Area*((ProjPhi/Gamma)+(ProjVel*Psi_E/Gamma_Minus_One));
				}
			}
			Total_CSens_Mach += CSens_Mach[iMarker];
		}
	}
}

void CAdjEulerSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {

	unsigned long iVertex, iPoint;
	double *d, *Normal, *U, *Psi_Aux, ProjVel = 0.0, bcn, vn = 0.0, Area, *UnitaryNormal, *Coord;
  double *Velocity, *Psi, Enthalpy = 0.0, sq_vel, phin, phis1, phis2, DensityInc = 0.0, BetaInc2 = 0.0;
	unsigned short iDim, iVar, jDim;
	
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();

	UnitaryNormal = new double[nDim];
	Velocity = new double[nDim];
	Psi      = new double[nVar];
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		if (geometry->node[iPoint]->GetDomain()) {
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Coord = geometry->node[iPoint]->GetCoord();
            
			/*--- Create a copy of the adjoint solution ---*/
			Psi_Aux = node[iPoint]->GetSolution();
			for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];			

			/*--- Flow solution ---*/
			U = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
			d = node[iPoint]->GetForceProj_Vector();
			
			Area = 0; 
			for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
			Area = sqrt(Area);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim]   = -Normal[iDim]/Area;

			if (incompressible) {
				
				DensityInc = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
				BetaInc2 = solution_container[FLOW_SOL]->node[iPoint]->GetBetaInc2();
				
				for (iDim = 0; iDim < nDim; iDim++)
					Velocity[iDim] = U[iDim+1] / solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
								
				/*--- Compute projections ---*/
				bcn = 0.0; phin = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					bcn += d[iDim]*UnitaryNormal[iDim];
					phin += Psi[iDim+1]*UnitaryNormal[iDim];
				}
				
				/*--- Introduce the boundary condition ---*/
				for (iDim = 0; iDim < nDim; iDim++) 
					Psi[iDim+1] -= ( phin - bcn ) * UnitaryNormal[iDim];
				
				/*--- Inner products after introducing BC (Psi has changed) ---*/
				phis1 = 0.0; phis2 = Psi[0] * (BetaInc2 / DensityInc);
				for (iDim = 0; iDim < nDim; iDim++) {
					phis1 -= Normal[iDim]*Psi[iDim+1];
					phis2 += Velocity[iDim]*Psi[iDim+1];
				}
				
				/*--- Flux of the Euler wall ---*/
				Residual[0] = phis1;
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[iDim+1] = - phis2 * Normal[iDim];
				
			}
			
			else {
				
				for (iDim = 0; iDim < nDim; iDim++)
					Velocity[iDim] = U[iDim+1] / U[0];
				
				Enthalpy = solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy();
				sq_vel   = 0.5*solution_container[FLOW_SOL]->node[iPoint]->GetVelocity2();
				
				/*--- Compute projections ---*/
				ProjVel = 0.0; bcn = 0.0; vn = 0.0, phin = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					ProjVel -= Velocity[iDim]*Normal[iDim];
					bcn     += d[iDim]*UnitaryNormal[iDim];
					vn      += Velocity[iDim]*UnitaryNormal[iDim];
					phin    += Psi[iDim+1]*UnitaryNormal[iDim];
				}
				
				/*--- Rotating Frame ---*/
				if (rotating_frame) {
					double ProjRotVel = 0.0;
					double *RotVel = geometry->node[iPoint]->GetRotVel();
					for (iDim = 0; iDim < nDim; iDim++) {
						ProjRotVel += RotVel[iDim]*UnitaryNormal[iDim];
					}
					ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux()/Area;
					phin -= Psi[nVar-1]*ProjRotVel;				
				}
				
				/*--- Introduce the boundary condition ---*/
				for (iDim = 0; iDim < nDim; iDim++) 
					Psi[iDim+1] -= ( phin - bcn ) * UnitaryNormal[iDim];
				
				/*--- Inner products after introducing BC (Psi has changed) ---*/
				phis1 = 0.0; phis2 = Psi[0] + Enthalpy * Psi[nVar-1];
				for (iDim = 0; iDim < nDim; iDim++) {
					phis1 -= Normal[iDim]*Psi[iDim+1];
					phis2 += Velocity[iDim]*Psi[iDim+1];
				}
				
				/*--- Flux of the Euler wall ---*/
				Residual[0] = ProjVel * Psi[0] - phis2 * ProjVel + phis1 * Gamma_Minus_One * sq_vel;
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[iDim+1] = ProjVel * Psi[iDim+1] - phis2 * Normal[iDim] - phis1 * Gamma_Minus_One * Velocity[iDim];
				Residual[nVar-1] = ProjVel * Psi[nVar-1] + phis1 * Gamma_Minus_One;
				
				/*--- Rotating Frame ---*/
				if (rotating_frame) {
					double ProjRotVel = 0.0;
					double *RotVel = geometry->node[iPoint]->GetRotVel();
					for (iDim = 0; iDim < nDim; iDim++)
						ProjRotVel -= RotVel[iDim]*Normal[iDim];
					ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux();
					Residual[0] -= ProjRotVel*Psi[0];
					for (iDim = 0; iDim < nDim; iDim++)
						Residual[iDim+1] -= ProjRotVel*Psi[iDim+1];
					Residual[nVar-1] -= ProjRotVel*Psi[nVar-1];
				}
			}
			
			/*--- Update residual ---*/
			node[iPoint]->SubtractRes_Conv(Residual);

			/*--- Implicit stuff ---*/
			if (implicit) {
				if (incompressible) {
					
					/*--- Adjoint density ---*/
					Jacobian_ii[0][0] = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_ii[0][iDim+1] = - Normal[iDim];
					
					/*--- Adjoint velocities ---*/
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_ii[iDim+1][0] = -Normal[iDim] * (BetaInc2 / DensityInc) ;
						for (jDim = 0; jDim < nDim; jDim++)
							Jacobian_ii[iDim+1][jDim+1] = - Normal[iDim] * Velocity[jDim];
					}
					
				}
				else {
					
					/*--- Adjoint density ---*/
					Jacobian_ii[0][0] = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_ii[0][iDim+1] = -ProjVel * (Velocity[iDim] - UnitaryNormal[iDim] * vn);
					Jacobian_ii[0][nVar-1] = -ProjVel * Enthalpy;
					
					/*--- Adjoint velocities ---*/
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_ii[iDim+1][0] = -Normal[iDim];
						for (jDim = 0; jDim < nDim; jDim++)
							Jacobian_ii[iDim+1][jDim+1] = -ProjVel*(UnitaryNormal[jDim]*UnitaryNormal[iDim] - Normal[iDim] * (Velocity[jDim] - UnitaryNormal[jDim] * vn));
						Jacobian_ii[iDim+1][iDim+1] += ProjVel;
						Jacobian_ii[iDim+1][nVar-1] = -Normal[iDim] * Enthalpy;
					}
					
					/*--- Adjoint energy ---*/
					Jacobian_ii[nVar-1][0] = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_ii[nVar-1][iDim+1] = 0.0;
					Jacobian_ii[nVar-1][nVar-1] = ProjVel;
					
					/*--- Contribution from a rotating frame ---*/
					if (rotating_frame) {
						double ProjRotVel = 0.0;
						double *RotVel = geometry->node[iPoint]->GetRotVel();
						for (iDim = 0; iDim < nDim; iDim++)
							ProjRotVel -= RotVel[iDim]*Normal[iDim];
						ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux();
						Jacobian_ii[0][0] -= ProjRotVel;
						for (iDim = 0; iDim < nDim; iDim++)
							Jacobian_ii[iDim+1][iDim+1] -= ProjRotVel;
						Jacobian_ii[nVar-1][nVar-1] -= ProjRotVel;
					}
				}
        
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
			}
		}
	}
	
	delete [] Velocity;
	delete [] UnitaryNormal;
	delete [] Psi;
}

void CAdjEulerSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
																		 CConfig *config, unsigned short val_marker) {
	
	unsigned long iVertex, iPoint;
	double *Normal, *U, *Psi_Aux, ProjVel = 0.0, vn = 0.0, Area, *UnitaryNormal, *Coord;
  double *Velocity, *Psi, Enthalpy = 0.0, sq_vel, phin, phis1, phis2, DensityInc = 0.0, BetaInc2 = 0.0;
	unsigned short iDim, iVar, jDim;
	
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();
	
	UnitaryNormal = new double[nDim];
	Velocity = new double[nDim];
	Psi      = new double[nVar];
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		if (geometry->node[iPoint]->GetDomain()) {
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Coord = geometry->node[iPoint]->GetCoord();
			
			/*--- Create a copy of the adjoint solution ---*/
			Psi_Aux = node[iPoint]->GetSolution();
			for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];			
			
			/*--- Flow solution ---*/
			U = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
			
			Area = 0; 
			for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
			Area = sqrt(Area);
			
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim]   = -Normal[iDim]/Area;
			
			if (incompressible) {
				
				DensityInc = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
				BetaInc2 = solution_container[FLOW_SOL]->node[iPoint]->GetBetaInc2();
				
				for (iDim = 0; iDim < nDim; iDim++)
					Velocity[iDim] = U[iDim+1] / solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
				
				/*--- Compute projections ---*/
				phin = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					phin += Psi[iDim+1]*UnitaryNormal[iDim];
				
				/*--- Introduce the boundary condition ---*/
				for (iDim = 0; iDim < nDim; iDim++) 
					Psi[iDim+1] -= phin * UnitaryNormal[iDim];
				
				/*--- Inner products after introducing BC (Psi has changed) ---*/
				phis1 = 0.0; phis2 = Psi[0] * (BetaInc2 / DensityInc);
				for (iDim = 0; iDim < nDim; iDim++) {
					phis1 -= Normal[iDim]*Psi[iDim+1];
					phis2 += Velocity[iDim]*Psi[iDim+1];
				}
				
				/*--- Flux of the Euler wall ---*/
				Residual[0] = phis1;
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[iDim+1] = - phis2 * Normal[iDim];
				
			}
			
			else {
				
				for (iDim = 0; iDim < nDim; iDim++)
					Velocity[iDim] = U[iDim+1] / U[0];
				
				Enthalpy = solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy();
				sq_vel   = 0.5*solution_container[FLOW_SOL]->node[iPoint]->GetVelocity2();
				
				/*--- Compute projections ---*/
				ProjVel = 0.0; vn = 0.0, phin = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					ProjVel -= Velocity[iDim]*Normal[iDim];
					vn      += Velocity[iDim]*UnitaryNormal[iDim];
					phin    += Psi[iDim+1]*UnitaryNormal[iDim];
				}
				
				/*--- Rotating Frame ---*/
				if (rotating_frame) {
					double ProjRotVel = 0.0;
					double *RotVel = geometry->node[iPoint]->GetRotVel();
					for (iDim = 0; iDim < nDim; iDim++) {
						ProjRotVel += RotVel[iDim]*UnitaryNormal[iDim];
					}
					ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux()/Area;
					phin -= Psi[nVar-1]*ProjRotVel;				
				}
				
				/*--- Introduce the boundary condition ---*/
				for (iDim = 0; iDim < nDim; iDim++) 
					Psi[iDim+1] -= phin * UnitaryNormal[iDim];
				
				/*--- Inner products after introducing BC (Psi has changed) ---*/
				phis1 = 0.0; phis2 = Psi[0] + Enthalpy * Psi[nVar-1];
				for (iDim = 0; iDim < nDim; iDim++) {
					phis1 -= Normal[iDim]*Psi[iDim+1];
					phis2 += Velocity[iDim]*Psi[iDim+1];
				}
				
				/*--- Flux of the Euler wall ---*/
				Residual[0] = ProjVel * Psi[0] - phis2 * ProjVel + phis1 * Gamma_Minus_One * sq_vel;
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[iDim+1] = ProjVel * Psi[iDim+1] - phis2 * Normal[iDim] - phis1 * Gamma_Minus_One * Velocity[iDim];
				Residual[nVar-1] = ProjVel * Psi[nVar-1] + phis1 * Gamma_Minus_One;
				
				/*--- Rotating Frame ---*/
				if (rotating_frame) {
					double ProjRotVel = 0.0;
					double *RotVel = geometry->node[iPoint]->GetRotVel();
					for (iDim = 0; iDim < nDim; iDim++)
						ProjRotVel -= RotVel[iDim]*Normal[iDim];
					ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux();
					Residual[0] -= ProjRotVel*Psi[0];
					for (iDim = 0; iDim < nDim; iDim++)
						Residual[iDim+1] -= ProjRotVel*Psi[iDim+1];
					Residual[nVar-1] -= ProjRotVel*Psi[nVar-1];
				}
			}
			
			/*--- Update residual ---*/
			node[iPoint]->SubtractRes_Conv(Residual);
			
			
			/*--- Implicit stuff ---*/
			if (implicit) {
				
				if (incompressible) {
					
					/*--- Adjoint density ---*/
					Jacobian_ii[0][0] = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_ii[0][iDim+1] = - Normal[iDim];
					
					/*--- Adjoint velocities ---*/
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_ii[iDim+1][0] = -Normal[iDim] * (BetaInc2 / DensityInc) ;
						for (jDim = 0; jDim < nDim; jDim++)
							Jacobian_ii[iDim+1][jDim+1] = - Normal[iDim] * Velocity[jDim];
					}
					
				}
				
				else {
					
					/*--- Adjoint density ---*/
					Jacobian_ii[0][0] = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_ii[0][iDim+1] = -ProjVel * (Velocity[iDim] - UnitaryNormal[iDim] * vn);
					Jacobian_ii[0][nVar-1] = -ProjVel * Enthalpy;
					
					/*--- Adjoint velocities ---*/
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_ii[iDim+1][0] = -Normal[iDim];
						for (jDim = 0; jDim < nDim; jDim++)
							Jacobian_ii[iDim+1][jDim+1] = -ProjVel*(UnitaryNormal[jDim]*UnitaryNormal[iDim] - Normal[iDim] * (Velocity[jDim] - UnitaryNormal[jDim] * vn));
						Jacobian_ii[iDim+1][iDim+1] += ProjVel;
						Jacobian_ii[iDim+1][nVar-1] = -Normal[iDim] * Enthalpy;
					}
					
					/*--- Adjoint energy ---*/
					Jacobian_ii[nVar-1][0] = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_ii[nVar-1][iDim+1] = 0.0;
					Jacobian_ii[nVar-1][nVar-1] = ProjVel;
					
					/*--- Contribution from a rotating frame ---*/
					if (rotating_frame) {
						double ProjRotVel = 0.0;
						double *RotVel = geometry->node[iPoint]->GetRotVel();
						for (iDim = 0; iDim < nDim; iDim++)
							ProjRotVel -= RotVel[iDim]*Normal[iDim];
						ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux();
						Jacobian_ii[0][0] -= ProjRotVel;
						for (iDim = 0; iDim < nDim; iDim++)
							Jacobian_ii[iDim+1][iDim+1] -= ProjRotVel;
						Jacobian_ii[nVar-1][nVar-1] -= ProjRotVel;
					}
				}
        
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
			}
		}
	}
	
	delete [] Velocity;
	delete [] UnitaryNormal;
	delete [] Psi;
}

void CAdjEulerSolution::BC_Interface_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
																							CConfig *config, unsigned short val_marker) {
	
#ifdef NO_MPI
	
	unsigned long iVertex, iPoint, jPoint;
	unsigned short iDim;
	double *Psi_i, *Psi_j, *U_i, *U_j, *Coord;
	
	double  *Normal = new double[nDim];
	
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPoint();
		Coord = geometry->node[iPoint]->GetCoord();
		
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Adjoint variables w/o reconstruction ---*/
			Psi_i = node[iPoint]->GetSolution();
			Psi_j = node[jPoint]->GetSolution();
			
			/*--- Conservative variables w/o reconstruction ---*/
			U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
			U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
			solver->SetConservative(U_i, U_j);
			
			/*--- SoundSpeed enthalpy and lambda variables w/o reconstruction ---*/
			solver->SetSoundSpeed(solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
														solution_container[FLOW_SOL]->node[jPoint]->GetSoundSpeed());
			solver->SetEnthalpy(solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
													solution_container[FLOW_SOL]->node[jPoint]->GetEnthalpy());
			
			/*--- Set face vector, and area ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = - Normal[iDim];
			solver->SetNormal(Normal);

			/*--- Just do a periodic BC ---*/
			solver->SetAdjointVar(Psi_i, Psi_j);
			
			/*--- Compute residual ---*/			
			solver->SetResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			node[iPoint]->SubtractRes_Conv(Res_Conv_i);
			
		}
	}
	
	delete[] Normal;
	
#else
	
	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	unsigned long iVertex, iPoint, jPoint;
	unsigned short iVar, iDim;
	double *Adjoint_Var, Psi_i[5], Psi_j[5], *U_i, *U_j;
	
	double *Normal = new double [nDim]; 
	double *Buffer_Send_Psi = new double[nVar];
	double *Buffer_Receive_Psi = new double[nVar];
	
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
				Adjoint_Var = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Send_Psi[iVar] = Adjoint_Var[iVar];
				MPI::COMM_WORLD.Bsend(Buffer_Send_Psi, nVar, MPI::DOUBLE, jProcessor, iPoint);
			}
		}
	}
	
	
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
			
			/*--- We only receive the information that belong to other boundary ---*/
			if (jProcessor != rank)
				MPI::COMM_WORLD.Recv(Buffer_Receive_Psi, nVar, MPI::DOUBLE, jProcessor, jPoint);
			else {
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Receive_Psi[iVar] = node[jPoint]->GetSolution(iVar); 
			}
			
			/*--- Store the solution for both points ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Psi_i[iVar] = node[iPoint]->GetSolution(iVar); 
				Psi_j[iVar] = Buffer_Receive_Psi[iVar]; 
			}
			
			/*--- Conservative variables w/o reconstruction (the same at both points) ---*/
			U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
			U_j = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
			solver->SetConservative(U_i, U_j);
			
			/*--- SoundSpeed enthalpy and lambda variables w/o reconstruction (the same at both points) ---*/
			solver->SetSoundSpeed(solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
														solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed());
			solver->SetEnthalpy(solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
													solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy());
			
			/*--- Set face vector, and area ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = - Normal[iDim];
			solver->SetNormal(Normal);
			
			/*--- Just do a periodic BC ---*/
			solver->SetAdjointVar(Psi_i, Psi_j);
			
			/*--- Compute residual ---*/			
			solver->SetResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			node[iPoint]->SubtractRes_Conv(Res_Conv_i);
		}
	}
	
	delete[] Buffer_Send_Psi;
	delete[] Buffer_Receive_Psi;
	delete[] Normal;
#endif
	
}

void CAdjEulerSolution::BC_NearField_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
																							CConfig *config, unsigned short val_marker) {
	
#ifdef NO_MPI
	
	unsigned long iVertex, iPoint, jPoint, Pin, Pout;
	unsigned short iVar, iDim;
	double  Psi_out[5], Psi_in[5], Psi_out_ghost[5], Psi_in_ghost[5], 
	MeanPsi[5], *Psi_i, *Psi_j, *U_i, *U_j, *IntBoundary_Jump, *Coord;

	double  *Normal = new double[nDim];
	
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPoint();
		Coord = geometry->node[iPoint]->GetCoord();
		
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Adjoint variables w/o reconstruction ---*/
			Psi_i = node[iPoint]->GetSolution();
			Psi_j = node[jPoint]->GetSolution();
			
			/*--- Conservative variables w/o reconstruction ---*/
			U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
			U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
			solver->SetConservative(U_i, U_j);
			
			/*--- SoundSpeed enthalpy and lambda variables w/o reconstruction ---*/
			solver->SetSoundSpeed(solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
														solution_container[FLOW_SOL]->node[jPoint]->GetSoundSpeed());
			solver->SetEnthalpy(solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
													solution_container[FLOW_SOL]->node[jPoint]->GetEnthalpy());
			
			/*--- Set face vector, and area ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = - Normal[iDim];
			solver->SetNormal(Normal);
			
			/*--- If equivalent area or nearfield pressure condition ---*/
			if ((config->GetKind_ObjFunc() == EQUIVALENT_AREA) || 
					(config->GetKind_ObjFunc() == NEARFIELD_PRESSURE)) {
				
				if (Normal[nDim-1] < 0.0) { Pin = iPoint; Pout = jPoint; }
				else { Pout = iPoint; Pin = jPoint; }
				
				for (iVar = 0; iVar < nVar; iVar++) {
					Psi_out[iVar] = node[Pout]->GetSolution(iVar);
					Psi_in[iVar] = node[Pin]->GetSolution(iVar);	
					MeanPsi[iVar] = 0.5*(Psi_out[iVar] + Psi_in[iVar]);
				}
				
				IntBoundary_Jump = node[iPoint]->GetIntBoundary_Jump();
				
				/*--- Inner point ---*/
				if (iPoint == Pin) {
					for (iVar = 0; iVar < nVar; iVar++)
						Psi_in_ghost[iVar] = 2.0*MeanPsi[iVar] - Psi_in[iVar] - IntBoundary_Jump[iVar];
					solver->SetAdjointVar(Psi_in, Psi_in_ghost);
				}
				
				/*--- Outer point ---*/
				if (iPoint == Pout) {
					for (iVar = 0; iVar < nVar; iVar++)
						Psi_out_ghost[iVar] = 2.0*MeanPsi[iVar] - Psi_out[iVar] + IntBoundary_Jump[iVar];
					solver->SetAdjointVar(Psi_out, Psi_out_ghost);
				}
			}
			else {
				/*--- Just do a periodic BC ---*/
				solver->SetAdjointVar(Psi_i, Psi_j);
			}
			
			/*--- Compute residual ---*/			
			solver->SetResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			node[iPoint]->SubtractRes_Conv(Res_Conv_i);
			
		}
	}
	
	delete[] Normal;

#else
	
	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	unsigned long iVertex, iPoint, jPoint, Pin, Pout;
	unsigned short iVar, iDim;
	double *Adjoint_Var, 
	Psi_out[5], Psi_in[5], Psi_i[5], Psi_j[5], Psi_in_ghost[5], Psi_out_ghost[5], MeanPsi[5], *U_i, *U_j, 
	*IntBoundary_Jump;
	
	double *Normal = new double [nDim]; 
	double *Buffer_Send_Psi = new double[nVar];
	double *Buffer_Receive_Psi = new double[nVar];
	
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
				Adjoint_Var = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Send_Psi[iVar] = Adjoint_Var[iVar];
				MPI::COMM_WORLD.Bsend(Buffer_Send_Psi, nVar, MPI::DOUBLE, jProcessor, iPoint);
			}
		}
	}
	
	
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
			
			/*--- We only receive the information that belong to other boundary ---*/
			if (jProcessor != rank)
				MPI::COMM_WORLD.Recv(Buffer_Receive_Psi, nVar, MPI::DOUBLE, jProcessor, jPoint);
			else {
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Receive_Psi[iVar] = node[jPoint]->GetSolution(iVar); 
			}
			
			/*--- Store the solution for both points ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Psi_i[iVar] = node[iPoint]->GetSolution(iVar); 
				Psi_j[iVar] = Buffer_Receive_Psi[iVar]; 
			}
			
			/*--- Conservative variables w/o reconstruction (the same at both points) ---*/
			U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
			U_j = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
			solver->SetConservative(U_i, U_j);
			
			/*--- SoundSpeed enthalpy and lambda variables w/o reconstruction (the same at both points) ---*/
			solver->SetSoundSpeed(solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
														solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed());
			solver->SetEnthalpy(solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
													solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy());
			
			/*--- Set face vector, and area ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = - Normal[iDim];
			solver->SetNormal(Normal);
			
			/*--- If equivalent area or nearfield pressure condition ---*/
			if ((config->GetKind_ObjFunc() == EQUIVALENT_AREA) || 
					(config->GetKind_ObjFunc() == NEARFIELD_PRESSURE)) {
				
				/*--- Inner nearfield boundary ---*/
				if (Normal[nDim-1] < 0.0)  { 
					Pin = iPoint; Pout = jPoint;
					for (iVar = 0; iVar < nVar; iVar++) {
						Psi_in[iVar] = Psi_i[iVar];
						Psi_out[iVar] = Psi_j[iVar];
						MeanPsi[iVar] = 0.5*(Psi_out[iVar] + Psi_in[iVar]);
					}
				}
				/*--- Outer nearfield boundary ---*/
				else { 
					Pout = iPoint; Pin = jPoint; 
					for (iVar = 0; iVar < nVar; iVar++) {
						Psi_in[iVar] = Psi_j[iVar];
						Psi_out[iVar] = Psi_i[iVar];
						MeanPsi[iVar] = 0.5*(Psi_out[iVar] + Psi_in[iVar]);
					}
				}
				
				IntBoundary_Jump = node[iPoint]->GetIntBoundary_Jump();
				
				/*--- Inner point ---*/
				if (iPoint == Pin) {
					for (iVar = 0; iVar < nVar; iVar++)
						Psi_in_ghost[iVar] = 2.0*MeanPsi[iVar] - Psi_in[iVar] - IntBoundary_Jump[iVar];
					solver->SetAdjointVar(Psi_in, Psi_in_ghost);
				}
				
				/*--- Outer point ---*/
				if (iPoint == Pout) {
					for (iVar = 0; iVar < nVar; iVar++)
						Psi_out_ghost[iVar] = 2.0*MeanPsi[iVar] - Psi_out[iVar] + IntBoundary_Jump[iVar];
					solver->SetAdjointVar(Psi_out, Psi_out_ghost);	
				}
			}
			else {
				/*--- Just do a periodic BC ---*/
				solver->SetAdjointVar(Psi_i, Psi_j);
			}
			
			/*--- Compute residual ---*/			
			solver->SetResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			node[iPoint]->SubtractRes_Conv(Res_Conv_i);
		}
	}
	
	delete[] Buffer_Send_Psi;
	delete[] Buffer_Receive_Psi;
	delete[] Normal;
#endif	
}

void CAdjEulerSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, 
									 CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	double *Normal, *U_domain, *U_infty, *Psi_domain, *Psi_infty, Pressure;
	
  bool rotating_frame = config->GetRotating_Frame();
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
	bool gravity = config->GetGravityForce();

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
				U_domain[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			
			/*--- Solution at the infinity ---*/
			if (incompressible) {
				double yFreeSurface = config->GetFreeSurface_Zero();
				double PressFreeSurface = solution_container[FLOW_SOL]->GetPressure_Inf();
				double Density = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
				double Froude = config->GetFroude();
				double yCoord = geometry->node[iPoint]->GetCoord(1);
				if (gravity) Pressure = PressFreeSurface + Density*(yFreeSurface-yCoord)/(Froude*Froude);
				else Pressure = PressFreeSurface;
				
				U_infty[0] = Pressure;
				U_infty[1] = solution_container[FLOW_SOL]->GetVelocity_Inf(0)*Density;
				U_infty[2] = solution_container[FLOW_SOL]->GetVelocity_Inf(1)*Density;
				if (nDim == 3) U_infty[3] = solution_container[FLOW_SOL]->GetVelocity_Inf(2)*Density;
			}
			else {
				/*--- Flow Solution at the infinity ---*/
				U_infty[0] = solution_container[FLOW_SOL]->GetDensity_Inf();
				U_infty[1] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(0);
				U_infty[2] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(1);
				U_infty[3] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
				if (nDim == 3) {
					U_infty[3] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(2);
					U_infty[4] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
				}
			}
			solver->SetConservative(U_domain, U_infty);

			/*--- Adjoint flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
				Psi_infty[iVar] = 0.0;
			}
			solver->SetAdjointVar(Psi_domain, Psi_infty);

			if (incompressible) {
				solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 
															solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
				solver->SetBetaInc2(solution_container[FLOW_SOL]->node[iPoint]->GetBetaInc2(), 
														solution_container[FLOW_SOL]->node[iPoint]->GetBetaInc2());
				solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
			}
			else {		
				solver->SetSoundSpeed(solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
															solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed());		
				solver->SetEnthalpy(solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
														solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy());
			}
			
			/*--- Rotational frame ---*/
			if (rotating_frame) {
				solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
        solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
			
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

void CAdjEulerSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, Point_Normal;
	unsigned short iVar, iDim;
	double *Normal, *U_domain, *U_inlet, *Psi_domain, *Psi_inlet;
	
  bool rotating_frame = config->GetRotating_Frame();
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
	
	Normal = new double[nDim];
	U_domain = new double[nVar]; U_inlet = new double[nVar];
	Psi_domain = new double[nVar]; Psi_inlet = new double[nVar];
	
	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
			solver->SetNormal(Normal);
			
			/*--- Set the normal point ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor();

			/*--- Flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			
			/*--- Flow Solution at the inlet ---*/
			if (incompressible) {
				U_inlet[0] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(0); 				
				for (iDim = 0; iDim < nDim; iDim++)
					U_inlet[iDim+1] = solution_container[FLOW_SOL]->GetVelocity_Inf(iDim)*
					solution_container[FLOW_SOL]->node[Point_Normal]->GetDensityInc();
				U_inlet[nDim] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(nDim);
			}			
			else {
				U_inlet[0] = solution_container[FLOW_SOL]->GetDensity_Inf();
				U_inlet[1] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(0);
				U_inlet[2] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(1);
				U_inlet[3] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
				if (nDim == 3) {
					U_inlet[3] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(2);
					U_inlet[4] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
				}
			}
			
			solver->SetConservative(U_domain, U_inlet);
			
			/*--- Adjoint flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Adjoint solution at the inlet ---*/
			if (incompressible) {
				Psi_inlet[0] = node[Point_Normal]->GetSolution(0);
				Psi_inlet[1] = 0.0;
				Psi_inlet[2] = 0.0;
			}			
			else {
				for (iVar = 0; iVar < nVar; iVar++)
					Psi_inlet[iVar] = 0.0;
			}
			
			solver->SetAdjointVar(Psi_domain, Psi_inlet);
			
			if (incompressible) {
				solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
				solver->SetBetaInc2(solution_container[FLOW_SOL]->node[iPoint]->GetBetaInc2(), solution_container[FLOW_SOL]->node[iPoint]->GetBetaInc2());
				solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
			}
			else {		
				solver->SetSoundSpeed(solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
															solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed());		
				solver->SetEnthalpy(solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
														solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy());
			}
			
			/*--- Rotational frame ---*/
			if (rotating_frame) {
				solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
        solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
			
			solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			/*--- Add and Subtract Residual ---*/
			node[iPoint]->SubtractRes_Conv(Residual_i);
			
			/*--- Implicit contribution to the residual ---*/
			if (implicit) 
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}
	
	delete [] Normal;
	delete [] U_domain; delete [] U_inlet;
	delete [] Psi_domain; delete [] Psi_inlet;
}

void CAdjEulerSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, Point_Normal;
	unsigned short iVar, iDim;
	double *Normal, *U_domain, *U_outlet, *Psi_domain, *Psi_outlet;
	
  bool rotating_frame = config->GetRotating_Frame();
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
	bool gravity = config->GetGravityForce();

	Normal = new double[nDim];
	U_domain = new double[nVar]; U_outlet = new double[nVar];
	Psi_domain = new double[nVar]; Psi_outlet = new double[nVar];
	
	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
			solver->SetNormal(Normal);
			
			/*--- Set the normal point ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetClosest_Neighbor();

			/*--- Flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
			
			/*--- Flow Solution at the outlet ---*/
			if (incompressible) {
				
				double yFreeSurface = config->GetFreeSurface_Zero();
				double PressFreeSurface = solution_container[FLOW_SOL]->GetPressure_Inf();
				double Density = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
				double Froude = config->GetFroude();
				double yCoord = geometry->node[iPoint]->GetCoord(nDim-1);
				double Pressure;
				if (gravity) Pressure = PressFreeSurface + Density*((yFreeSurface-yCoord)/(Froude*Froude));
				else Pressure = PressFreeSurface;
								
				/*--- The pressure is computed from the infinity ---*/				
				U_outlet[0] = Pressure;
				
				/*--- Velocity computation using internal value ---*/
				for (iDim = 0; iDim < nDim; iDim++)
					U_outlet[iDim+1] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(iDim+1);
			}
			else {
				U_outlet[0] = solution_container[FLOW_SOL]->GetDensity_Inf();
				U_outlet[1] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(0);
				U_outlet[2] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(1);
				U_outlet[3] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
				if (nDim == 3) {
					U_outlet[3] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(2);
					U_outlet[4] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
				}
			}
			
			solver->SetConservative(U_domain, U_outlet);
			
			/*--- Adjoint flow solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Adjoint flow solution at the outlet ---*/
			if (incompressible) {
				Psi_outlet[2] = 0.0;				
//				double coeff = (2.0*U_domain[1])/ solution_container[FLOW_SOL]->node[iPoint]->GetBetaInc2();
				Psi_outlet[1] = 0.0; //node[Point_Normal]->GetSolution(1);
				Psi_outlet[0] = 0.0; //-coeff*Psi_outlet[1];
			}
			else {
				for (iVar = 0; iVar < nVar; iVar++)
					Psi_outlet[iVar] = 0.0;
			}
			
			solver->SetAdjointVar(Psi_domain, Psi_outlet);
			
			if (incompressible) {
				solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
				solver->SetBetaInc2(solution_container[FLOW_SOL]->node[iPoint]->GetBetaInc2(), solution_container[FLOW_SOL]->node[iPoint]->GetBetaInc2());
				solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
			}
			else {		
				solver->SetSoundSpeed(solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed(), 
															solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed());		
				solver->SetEnthalpy(solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy(), 
														solution_container[FLOW_SOL]->node[iPoint]->GetEnthalpy());
			}
			
			/*--- Rotational frame ---*/
			if (rotating_frame) {
				solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
        solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
			
			solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			
			/*--- Add and Subtract Residual ---*/
			node[iPoint]->SubtractRes_Conv(Residual_i);
			
			/*--- Implicit contribution to the residual ---*/
			if (implicit) 
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
		}
	}
	
	delete [] Normal;
	delete [] U_domain; delete [] U_outlet;
	delete [] Psi_domain; delete [] Psi_outlet;
}

void CAdjEulerSolution::MPI_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
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

CAdjNSSolution::CAdjNSSolution(void) : CAdjEulerSolution() { }

CAdjNSSolution::CAdjNSSolution(CGeometry *geometry, CConfig *config) : CAdjEulerSolution() {
	unsigned long iPoint, index, iVertex;
	string text_line, mesh_filename;
	unsigned short iDim, iVar, iMarker;
	ifstream restart_file;
	string filename, AdjExt;

	bool restart = config->GetRestart();
	bool incompressible = config->GetIncompressible();

	/*--- Set the gamma value ---*/
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	
	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	if (incompressible) nVar = nDim + 1;
	else nVar = nDim + 2;
	node = new CVariable*[geometry->GetnPoint()];
	
	/*--- Define some auxiliar vector related with the residual ---*/
	Residual = new double[nVar];	Residual_Max = new double[nVar];
	Residual_i = new double[nVar];	Residual_j = new double[nVar];
	Res_Conv_i = new double[nVar];	Res_Visc_i = new double[nVar];
	Res_Conv_j = new double[nVar];	Res_Visc_j = new double[nVar];
	Res_Sour_j = new double[nVar];	Res_Sour_j = new double[nVar];
	
	/*--- Define some auxiliar vector related with the solution ---*/
	Solution   = new double[nVar];	
	Solution_i = new double[nVar];	Solution_j = new double[nVar];
	
	/*--- Define some auxiliar vector related with the geometry ---*/
	Vector_i = new double[nDim];	Vector_j = new double[nDim];
	
	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT) {
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

	/*--- Sensitivity definition and coefficient in all the markers ---*/
	CSensitivity = new double* [config->GetnMarker_All()];
	for (iMarker=0; iMarker<config->GetnMarker_All(); iMarker++) {
		CSensitivity[iMarker] = new double [geometry->nVertex[iMarker]];
	}
	CSens_Geo = new double[config->GetnMarker_All()];
	CSens_Mach = new double[config->GetnMarker_All()];
	CSens_AoA = new double[config->GetnMarker_All()];
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		CSens_Geo[iMarker] = 0.0;
		CSens_Mach[iMarker] = 0.0;
		CSens_AoA[iMarker] = 0.0;
		for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
			CSensitivity[iMarker][iVertex] = 0.0;
	}
	
	/*--- Adjoint flow infinity inizialization stuff ---*/
	PsiRho_Inf = 0.0; PsiE_Inf = 0.0;
	Phi_Inf = new double [nDim];
	Phi_Inf[0] = 0.0; Phi_Inf[1] = 0.0;
	if (nDim == 3) Phi_Inf[2] = 0.0;
		
	if (!restart) {
		/*--- Restart the solution from infinity ---*/
		for (iPoint = 0 ; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CAdjNSVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, nDim, nVar, config);
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
			case ELECTRIC_CHARGE: AdjExt = "_ec.dat"; break;
			case EQUIVALENT_AREA: AdjExt = "_ea.dat"; break;
			case NEARFIELD_PRESSURE: AdjExt = "_nfp.dat"; break;
      case FORCE_X_COEFFICIENT: AdjExt = "_cfx.dat"; break;
			case FORCE_Y_COEFFICIENT: AdjExt = "_cfy.dat"; break;
			case FORCE_Z_COEFFICIENT: AdjExt = "_cfz.dat"; break;
      case THRUST_COEFFICIENT: AdjExt = "_ct.dat"; break;
      case TORQUE_COEFFICIENT: AdjExt = "_cq.dat"; break;
      case FIGURE_OF_MERIT: AdjExt = "_merit.dat"; break;
			case FREESURFACE: AdjExt = "_fs.dat"; break;
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
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
			if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
			node[iPoint] = new CAdjNSVariable(Solution, nDim, nVar, config);
		}
		restart_file.close();
	}
}

CAdjNSSolution::~CAdjNSSolution(void) {
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
	delete [] Vector_i; delete [] Vector_j;
	delete [] xsol; delete [] rhs;
	delete [] CSens_Geo; delete [] CSens_Mach;
	delete [] CSens_AoA; delete [] Phi_Inf;
	
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


void CAdjNSSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep) {
	unsigned long iPoint;
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);

	/*--- Residual inicialization ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		node[iPoint]->Set_ResConv_Zero();
		node[iPoint]->Set_ResSour_Zero();
		if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit)
			node[iPoint]->Set_ResVisc_Zero();
	}
	
	/*--- Compute gradients for solution reconstruction and viscous term (be careful, if a upwind 
	 strategy is used, then we compute twice the gradient) ---*/
	switch (config->GetKind_Gradient_Method()) {
		case GREEN_GAUSS :
			SetSolution_Gradient_GG(geometry);
			if (config->GetKind_Solver() == ADJ_RANS)
				solution_container[ADJTURB_SOL]->SetSolution_Gradient_GG(geometry);
			break;
		case LEAST_SQUARES : case WEIGHTED_LEAST_SQUARES :
			SetSolution_Gradient_LS(geometry, config);
			if (config->GetKind_Solver() == ADJ_RANS)
				solution_container[ADJTURB_SOL]->SetSolution_Gradient_LS(geometry, config);
			break;
	}
	
	/*--- Implicit solution ---*/
	if (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT)
		Jacobian.SetValZero();
	
}

void CAdjNSSolution::Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
									  CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	
	unsigned long iEdge, iPoint, jPoint;
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);

	if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit) {
		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
			
			/*--- Points in edge, coordinates and normal vector---*/
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);
			solver->SetCoord(geometry->node[iPoint]->GetCoord(),geometry->node[jPoint]->GetCoord());
			solver->SetNormal(geometry->edge[iEdge]->GetNormal());
			
			/*--- Conservative variables w/o reconstruction and adjoint variables w/o reconstruction---*/
			solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), 
															solution_container[FLOW_SOL]->node[jPoint]->GetSolution());
			solver->SetAdjointVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
			
			/*--- Gradient of Adjoint Variables ---*/
			solver->SetAdjointVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
			
			/*--- Viscosity and eddy viscosity---*/
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 
																	solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
			solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(), 
															 solution_container[FLOW_SOL]->node[jPoint]->GetEddyViscosity());
			
			/*--- Compute residual in a non-conservative way, and update ---*/
			solver->SetResidual(Residual_i, Residual_j);
			node[iPoint]->SubtractRes_Visc(Residual_i);
			node[jPoint]->AddRes_Visc(Residual_j);
			
			/*--- Compute residual in a non-conservative way, and update ---*/
			/*		solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			 node[iPoint]->SubtractRes_Visc(Residual_i);
			 node[jPoint]->AddRes_Visc(Residual_j);	
			 Jacobian.SubtractBlock(iPoint, iPoint ,Jacobian_ii);
			 Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
			 Jacobian.AddBlock(jPoint, iPoint, Jacobian_ji);
			 Jacobian.AddBlock(jPoint, jPoint, Jacobian_jj);*/
		}
	}
}

void CAdjNSSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
												  CConfig *config, unsigned short iMesh) {

    unsigned long iPoint;

	if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
		solution_container[FLOW_SOL]->SetPrimVar_Gradient_GG(geometry, config);

	if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) || 
		(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES))
		solution_container[FLOW_SOL]->SetPrimVar_Gradient_LS(geometry, config);

	/*--- Loop over all the points, note that we are supposing that primitive and 
	 adjoint gradients have been computed previously ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		/*--- Conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);
		
		/*--- Gradient of primitive and adjoint variables ---*/
		solver->SetPrimVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
		solver->SetAdjointVarGradient(node[iPoint]->GetGradient(), NULL);
		
		/*--- Laminar viscosity, and eddy viscosity (adjoint with frozen viscosity) ---*/
		solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
		solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(), 0.0);
		
		/*--- Set temperature of the fluid ---*/
		solver->SetTemperature(solution_container[FLOW_SOL]->node[iPoint]->GetTemperature(), 0.0);
		
		/*--- Set volume ---*/
		solver->SetVolume(geometry->node[iPoint]->GetVolume());
		
		/*--- If turbulence computation we must add some coupling terms to the NS adjoint eq. ---*/
		if (config->GetKind_Solver() == ADJ_RANS) {
						
			/*--- Turbulent variables w/o reconstruction ---*/
			solver->SetTurbVar(solution_container[TURB_SOL]->node[iPoint]->GetSolution(), NULL); 

			/*--- Gradient of Turbulent Variables w/o reconstruction ---*/
			solver->SetTurbVarGradient(solution_container[TURB_SOL]->node[iPoint]->GetGradient(), NULL);
			
			/*--- Turbulent adjoint variables w/o reconstruction ---*/
			solver->SetTurbAdjointVar(solution_container[ADJTURB_SOL]->node[iPoint]->GetSolution(), NULL);
			
			/*--- Gradient of Adjoint turbulent variables w/o reconstruction ---*/
			solver->SetTurbAdjointGradient(solution_container[ADJTURB_SOL]->node[iPoint]->GetGradient(), NULL);
			
			/*--- Set distance to the surface ---*/
			solver->SetDistance(geometry->node[iPoint]->GetWallDistance(), 0.0);
		}
		
		/*--- Compute and update residual ---*/
		solver->SetResidual(Residual, config);
		node[iPoint]->AddRes_Conv(Residual);
	}
}

void CAdjNSSolution::SourceConserv_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
												  CConfig *config, unsigned short iMesh) {

	return;

	unsigned long iEdge, iPoint, jPoint;

	if (config->GetKind_Solver() == ADJ_RANS) {

		/*--- Gradient of primitive variables already computed in the previous step ---*/
		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
			
			/*--- Points in edge, and normal vector ---*/
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);
			solver->SetNormal(geometry->edge[iEdge]->GetNormal());

			/*--- Conservative variables w/o reconstruction ---*/
			solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), 
									solution_container[FLOW_SOL]->node[jPoint]->GetSolution());

			/*--- Gradient of primitive variables w/o reconstruction ---*/
			solver->SetPrimVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), 
									   solution_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive());

			/*--- Viscosity ---*/
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 
										solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());

			/*--- Turbulent variables w/o reconstruction ---*/
			solver->SetTurbVar(solution_container[TURB_SOL]->node[iPoint]->GetSolution(), 
							   solution_container[TURB_SOL]->node[jPoint]->GetSolution()); 

			/*--- Turbulent adjoint variables w/o reconstruction ---*/
			solver->SetTurbAdjointVar(solution_container[ADJTURB_SOL]->node[iPoint]->GetSolution(), 
									  solution_container[ADJTURB_SOL]->node[jPoint]->GetSolution());

			/*--- Set distance to the surface ---*/
			solver->SetDistance(geometry->node[iPoint]->GetWallDistance(), geometry->node[jPoint]->GetWallDistance());

			/*--- Add and Subtract Residual ---*/
			solver->SetResidual(Residual, config);
			node[iPoint]->AddRes_Conv(Residual);
			node[jPoint]->SubtractRes_Conv(Residual);
		}
	}
}

void CAdjNSSolution::Viscous_Sensitivity(CGeometry *geometry, CSolution **solution_container, CConfig *config) {

	unsigned long iVertex, iPoint;
	unsigned short iDim, jDim, iMarker;
	double **PsiVar_Grad, **PrimVar_Grad, mue_eff, kappa, div_phi, *Normal, Area, 
	normal_grad_psi5, normal_grad_T, tang_psi_5, sigma_partial, *Psi, 
	Mach_Inf_3, Mach_Inf, ProjPhi, *U, ProjVel, Psi_E;
	double invR = Gamma*config->GetMach_FreeStreamND()*config->GetMach_FreeStreamND();
	double Cp = (Gamma/Gamma_Minus_One)/invR;

	double *UnitaryNormal = new double[nDim];
	double *normal_grad_vel = new double[nDim];
	double *tang_deriv_psi5 = new double[nDim]; 
	double *tang_deriv_T = new double[nDim]; 
	double **Sigma = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++) 
		Sigma[iDim] = new double [nDim];
	
	Total_CSens_Geo = 0.0; Total_CSens_Mach = 0.0; Total_CSens_AoA = 0.0;
	
	/*--- Compute gradient of adjoint variables on the surface ---*/
	SetSurface_Gradient(geometry, config);
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		CSens_Geo[iMarker] = 0.0;
		if (config->GetMarker_All_Boundary(iMarker) == NO_SLIP_WALL) {
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					PsiVar_Grad = node[iPoint]->GetGradient();
					PrimVar_Grad = solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
					mue_eff = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
					kappa = Cp*mue_eff/PRANDTL;
					
					/*--- Components of the effective and adjoint stress tensors ---*/
					div_phi = 0;
					for (iDim = 0; iDim < nDim; iDim++) {
						div_phi += PsiVar_Grad[iDim+1][iDim];
						for (jDim = 0; jDim < nDim; jDim++) 
							Sigma[iDim][jDim] = mue_eff*(PsiVar_Grad[iDim+1][jDim]+PsiVar_Grad[jDim+1][iDim]);
					}
					for (iDim = 0; iDim < nDim; iDim++) 
						Sigma[iDim][iDim] -= TWO3*mue_eff*div_phi;
					
					/*--- Compute face area and the nondimensional normal to the surface ---*/
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					Area = 0; 
					for (iDim = 0; iDim < nDim; iDim++) { Area += Normal[iDim]*Normal[iDim]; } Area = sqrt(Area);
					for (iDim = 0; iDim < nDim; iDim++) { UnitaryNormal[iDim] = Normal[iDim]/Area; } 
					
					/*--- Compute normal derivatives of velocities and tangencial deriv. of psi5 ---*/
					normal_grad_psi5 = 0; normal_grad_T = 0;
					for (iDim = 0; iDim < nDim; iDim++) {
						normal_grad_psi5 += PsiVar_Grad[nVar-1][iDim]*UnitaryNormal[iDim];
						normal_grad_T += PrimVar_Grad[0][iDim]*UnitaryNormal[iDim];
						normal_grad_vel[iDim] = 0;
						for (jDim = 0; jDim < nDim; jDim++)
							normal_grad_vel[iDim] += PrimVar_Grad[iDim+1][jDim]*UnitaryNormal[jDim];
					}
					
					for (iDim = 0; iDim < nDim; iDim++) {
						tang_deriv_psi5[iDim] = PsiVar_Grad[nVar-1][iDim] - normal_grad_psi5*UnitaryNormal[iDim];
						tang_deriv_T[iDim] = PrimVar_Grad[0][iDim] - normal_grad_T*UnitaryNormal[iDim];
					}
					
					/*--- tang_psi_5 = (\partial_tg \psi_5)\cdot (k \partial_tg T) ---*/
					tang_psi_5 = 0;
					for (iDim = 0; iDim < nDim; iDim++) 
						tang_psi_5 += kappa*tang_deriv_psi5[iDim]*tang_deriv_T[iDim]; 
					
					/*--- sigma_partial = \Sigma_{ji} n_i \partial_n v_j ---*/
					sigma_partial = 0;
					for (iDim = 0; iDim < nDim; iDim++)
						for (jDim = 0; jDim < nDim; jDim++) 
							sigma_partial += UnitaryNormal[iDim]*Sigma[iDim][jDim]*normal_grad_vel[jDim];

					/*--- Compute sensitivity for each surface point ---*/
					CSensitivity[iMarker][iVertex] = (sigma_partial - tang_psi_5)*Area;
					CSens_Geo[iMarker] -= CSensitivity[iMarker][iVertex]*Area;
				}
			}
			Total_CSens_Geo += CSens_Geo[iMarker];
		}
	}
	
	/*--- Mach number sensitivity ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == FAR_FIELD) {
			CSens_Mach[iMarker] = 0.0;
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					Psi = node[iPoint]->GetSolution();
					U = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					Mach_Inf = config->GetMach_FreeStreamND();
					Mach_Inf_3 = config->GetMach_FreeStreamND()*config->GetMach_FreeStreamND()*config->GetMach_FreeStreamND();
					
					Area = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) 
						Area += Normal[iDim]*Normal[iDim];
					Area = sqrt(Area);
					
					for (iDim = 0; iDim < nDim; iDim++)
						UnitaryNormal[iDim] = -Normal[iDim]/Area;
					
					ProjPhi = 0; ProjVel = 0;
					for (iDim = 0; iDim < nDim; iDim++) {
						ProjPhi += Psi[iDim+1]*UnitaryNormal[iDim];
						ProjVel += U[iDim+1]*UnitaryNormal[iDim]/U[0];
					}
					Psi_E = Psi[nVar-1];
					CSens_Mach[iMarker] += (2.0/Mach_Inf_3)*Area*((ProjPhi/Gamma)+(ProjVel*Psi_E/Gamma_Minus_One));
				}
			}
			Total_CSens_Mach += CSens_Mach[iMarker];
		}
	}
	
	delete [] UnitaryNormal;
	delete [] normal_grad_vel;
	delete [] tang_deriv_psi5;
	delete [] tang_deriv_T; 
	for (iDim = 0; iDim < nDim; iDim++) 
		delete [] Sigma[iDim];
	delete [] Sigma;
	
}

void CAdjNSSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short val_marker) {

	unsigned long iVertex, iPoint, total_index;
	unsigned short iDim, iVar, jDim;
	double *d, *Normal, *local_res_conv, *local_res_visc, *local_trunc_error, l1psi;
	bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	
	double invR = Gamma*config->GetMach_FreeStreamND()*config->GetMach_FreeStreamND();
	double Cp = (Gamma/Gamma_Minus_One)/invR;
//	double mu2 = 0.404;
	double mu_dyn, Temp, dVisc_T, rho, pressure;
	double div_phi, force_stress, Sigma_5, **PsiVar_Grad;
	double **Tau = new double* [nDim];
	double Gamma = config->GetGamma();
	double Gamma_Minus_One = Gamma - 1.0;
	
	for (iDim = 0; iDim < nDim; iDim++) 
		Tau[iDim] = new double [nDim];
	
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
			
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			d = node[iPoint]->GetForceProj_Vector();
			
			local_res_conv = node[iPoint]->GetResConv();
			local_res_visc = node[iPoint]->GetResVisc();
			local_trunc_error = node[iPoint]->GetTruncationError();
			
			/*--- Strong imposition of psi = ForceProj_Vector ---*/
			l1psi = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				node[iPoint]->SetSolution_Old(iDim+1,d[iDim]); 
				local_res_conv[iDim+1] = 0;
				local_res_visc[iDim+1] = 0;
				local_trunc_error[iDim+1] = 0;
				l1psi += Normal[iDim]*d[iDim];
			}
			
			local_res_conv[nVar-1] += l1psi*Gamma_Minus_One;
			
			/*--- Only change velocity-rows of the Jacobian ---*/
			if (implicit) {
				for (iVar = 1; iVar <= nDim; iVar++) {
					total_index = iPoint*nVar+iVar;
					Jacobian.DeleteValsRowi(total_index);
				}
			}
			
			/*--- Strong imposition of normal_grad_psi5 = g1.d/g2 ---*/
			
			/*--- Components of the effective and adjoint stress tensors ---*/
			PsiVar_Grad = node[iPoint]->GetGradient();
			div_phi = 0;
			for (iDim = 0; iDim < nDim; iDim++) {
				div_phi += PsiVar_Grad[iDim+1][iDim];
				for (jDim = 0; jDim < nDim; jDim++) 
					Tau[iDim][jDim] = (PsiVar_Grad[iDim+1][jDim]+PsiVar_Grad[jDim+1][iDim]);
			}
			for (iDim = 0; iDim < nDim; iDim++) 
				Tau[iDim][iDim] -= TWO3*div_phi;
			
			/*--- force_stress = n_i \Tau_{ij} d_j ---*/
			force_stress = 0;
			for (iDim = 0; iDim < nDim; iDim++)
				for (jDim = 0; jDim < nDim; jDim++) 
					force_stress += Normal[iDim]*Tau[iDim][jDim]*d[jDim];
			
			/*--- \partial \mu_dyn \partial T ---*/
			mu_dyn = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
			Temp = solution_container[FLOW_SOL]->node[iPoint]->GetTemperature();
//			dVisc_T = mu_dyn*(Temp+3.0*mu2)/(2.0*Temp*(Temp+mu2));
			dVisc_T = 0;
			
			/*--- \Sigma_5 ---*/
			Sigma_5 = (Gamma/Cp)*dVisc_T*force_stress;
			
			/*--- Imposition of residuals ---*/
			rho = solution_container[FLOW_SOL]->node[iPoint]->GetDensity();
			pressure = solution_container[FLOW_SOL]->node[iPoint]->GetPressure();
			local_res_visc[0] += pressure*Sigma_5/(Gamma_Minus_One*rho*rho);
			local_res_visc[nVar-1] -= Sigma_5/rho;
		}
	}
	
	
}

void CAdjNSSolution::MPI_Send_Receive(CGeometry *geometry, CSolution **solution_container, CConfig *config, 
										unsigned short val_mesh) {
#ifdef CHECK
#ifndef NO_MPI
	
	unsigned short iVar, iMarker;
	unsigned long iVertex, iPoint;
	double *Adjoint_Var, *Adjoint_Undivided_Laplacian = NULL, **Adjoint_Grad = NULL;
	
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			
			short SendRecv = config->GetMarker_All_SendRecv(iMarker);
			
			/*--- Send information  ---*/
			if (SendRecv > 0) {
				
				/*--- Upwind scheme ---*/
				if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {
					double Buffer_Send_Psi[geometry->nVertex[iMarker]][nVar];
					double Buffer_Send_Psix[geometry->nVertex[iMarker]][nVar];
					double Buffer_Send_Psiy[geometry->nVertex[iMarker]][nVar];
					double Buffer_Send_Psiz[geometry->nVertex[iMarker]][nVar];
					unsigned long nBuffer = geometry->nVertex[iMarker]*nVar;
					int send_to = SendRecv-1;
					
					for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						Adjoint_Var = node[iPoint]->GetSolution();
						if (val_mesh == MESH_0) Adjoint_Grad = node[iPoint]->GetGradient();
						
						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Send_Psi[iVertex][iVar] = Adjoint_Var[iVar];
							if (val_mesh == MESH_0) {
								Buffer_Send_Psix[iVertex][iVar] = Adjoint_Grad[iVar][0];					
								Buffer_Send_Psiy[iVertex][iVar] = Adjoint_Grad[iVar][1];
								if (nDim == 3) Buffer_Send_Psiz[iVertex][iVar] = Adjoint_Grad[iVar][2];
							}
						}
					}
					
					MPI::COMM_WORLD.Bsend(&Buffer_Send_Psi,nBuffer,MPI::DOUBLE,send_to, 0);
					if (val_mesh == MESH_0) {
						MPI::COMM_WORLD.Bsend(&Buffer_Send_Psix,nBuffer,MPI::DOUBLE,send_to, 1);
						MPI::COMM_WORLD.Bsend(&Buffer_Send_Psiy,nBuffer,MPI::DOUBLE,send_to, 2);
						if (nDim == 3) MPI::COMM_WORLD.Bsend(&Buffer_Send_Psiz,nBuffer,MPI::DOUBLE,send_to, 3);
					}
				}
				
				/*--- Centered scheme ---*/
				if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {
					
					double Buffer_Send_Psi[geometry->nVertex[iMarker]][nVar];
					double Buffer_Send_Undivided_Laplacian[geometry->nVertex[iMarker]][nVar];
					unsigned long nBuffer_Vector = geometry->nVertex[iMarker]*nVar;
					int send_to = SendRecv-1;
					
					for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						Adjoint_Var = node[iPoint]->GetSolution();
						if (val_mesh == MESH_0) Adjoint_Undivided_Laplacian = node[iPoint]->GetUnd_Lapl();
						for (iVar = 0; iVar < nVar; iVar++) {
							Buffer_Send_Psi[iVertex][iVar] = Adjoint_Var[iVar];
							if (val_mesh == MESH_0) Buffer_Send_Undivided_Laplacian[iVertex][iVar] = Adjoint_Undivided_Laplacian[iVar];					
						}
					}
					
					MPI::COMM_WORLD.Bsend(&Buffer_Send_Psi,nBuffer_Vector,MPI::DOUBLE,send_to, 0);
					if (val_mesh == MESH_0) MPI::COMM_WORLD.Bsend(&Buffer_Send_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,send_to, 1);
					
				}
			}
			
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
				/*--- Upwind scheme (Not validated)---*/
				if (config->GetKind_ConvNumScheme() == SPACE_UPWIND) {
					double Buffer_Receive_Psi[geometry->nVertex[iMarker]][nVar];
					double Buffer_Receive_Psix[geometry->nVertex[iMarker]][nVar];
					double Buffer_Receive_Psiy[geometry->nVertex[iMarker]][nVar];
					double Buffer_Receive_Psiz[geometry->nVertex[iMarker]][nVar];
					unsigned long nBuffer = geometry->nVertex[iMarker]*nVar;
					int receive_from = abs(SendRecv)-1;
					
					MPI::COMM_WORLD.Recv(&Buffer_Receive_Psi,nBuffer,MPI::DOUBLE,receive_from, 0);
					if (val_mesh == MESH_0) {
						MPI::COMM_WORLD.Recv(&Buffer_Receive_Psix,nBuffer,MPI::DOUBLE,receive_from, 1);
						MPI::COMM_WORLD.Recv(&Buffer_Receive_Psiy,nBuffer,MPI::DOUBLE,receive_from, 2);
						if (nDim == 3) MPI::COMM_WORLD.Recv(&Buffer_Receive_Psiz,nBuffer,MPI::DOUBLE,receive_from, 3);
					}
					
					for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();	
						for (iVar = 0; iVar < nVar; iVar++) {
							node[iPoint]->SetSolution(iVar, Buffer_Receive_Psi[iVertex][iVar]);
							if (val_mesh == MESH_0) {
								node[iPoint]->SetGradient(iVar, 0, Buffer_Receive_Psix[iVertex][iVar]);
								node[iPoint]->SetGradient(iVar, 1, Buffer_Receive_Psiy[iVertex][iVar]);
								if (nDim == 3) node[iPoint]->SetGradient(iVar, 2, Buffer_Receive_Psiz[iVertex][iVar]);
							}
						}
					}
				}
				
				/*--- Centered scheme ---*/
				if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {
					double Buffer_Receive_Psi[geometry->nVertex[iMarker]][nVar];
					double Buffer_Receive_Undivided_Laplacian[geometry->nVertex[iMarker]][nVar];
					unsigned long nBuffer_Vector = geometry->nVertex[iMarker]*nVar;
					int receive_from = abs(SendRecv)-1;
					
					MPI::COMM_WORLD.Recv(&Buffer_Receive_Psi,nBuffer_Vector,MPI::DOUBLE,receive_from, 0);
					if (val_mesh == MESH_0) MPI::COMM_WORLD.Recv(&Buffer_Receive_Undivided_Laplacian,nBuffer_Vector,MPI::DOUBLE,receive_from, 1);
					
					for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						for (iVar = 0; iVar < nVar; iVar++) {
							node[iPoint]->SetSolution(iVar, Buffer_Receive_Psi[iVertex][iVar]);
							if (val_mesh == MESH_0) node[iPoint]->SetUndivided_Laplacian(iVar, Buffer_Receive_Undivided_Laplacian[iVertex][iVar]);
						}
					}
				}			
			}
		}
	
#endif
#endif

}
