/*!
 * \file solution_adjoint_turbulent.cpp
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

CAdjTurbSolution::CAdjTurbSolution(void) : CSolution() {}

CAdjTurbSolution::CAdjTurbSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned long iPoint, nPoint = geometry->GetnPoint();
	unsigned short nMarker, iDim, iVar;//, nNeigh;

	nDim = geometry->GetnDim();
	nMarker = config->GetnMarker_All(); 
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Dimension of the problem --> dependent of the turbulent model ---*/
	switch (config->GetKind_Turb_Model()) {	
		case SA :
			nVar = 1;		
			break;
		case SST :
			nVar = 2;		
			break;
	}
    
    unsigned short nTotalVar = nVar + (nDim + 2); //nTotalVar = # turb vars + # flow vars

	Residual   = new double [nVar]; Residual_RMS = new double[nVar];
	Residual_i = new double [nVar]; Residual_j = new double [nVar];
	Residual_Max = new double [nVar]; Point_Max = new unsigned long[nVar];
  

	Solution   = new double [nVar];
	Solution_i = new double [nVar];
	Solution_j = new double [nVar];

	/*--- Define some auxiliar vector related with the geometry ---*/
	Vector_i = new double[nDim]; Vector_j = new double[nDim];

	/*--- Define some auxiliar vector related with the flow solution ---*/
	FlowSolution_i = new double [nDim+2]; FlowSolution_j = new double [nDim+2];

	/*--- Point to point Jacobians ---*/
	Jacobian_ii = new double* [nVar];
	Jacobian_ij = new double* [nVar];
	Jacobian_ji = new double* [nVar];
	Jacobian_jj = new double* [nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar++) {
		Jacobian_ii[iVar] = new double [nVar];
		Jacobian_ij[iVar] = new double [nVar];
		Jacobian_ji[iVar] = new double [nVar];
		Jacobian_jj[iVar] = new double [nVar];
	}
    
    // Hybrid Jacobians
    Jacobian_i = new double*[nTotalVar];
    Jacobian_j = new double*[nTotalVar];
    for (unsigned short iVar = 0; iVar < nTotalVar; iVar++) {
        Jacobian_i[iVar] = new double [nVar];
        Jacobian_j[iVar] = new double [nVar];
    }
    
    Jacobian_mui = new double [nTotalVar];
    Jacobian_muj = new double [nTotalVar];
    
    Jacobian_gradi = new double**[nTotalVar];
    Jacobian_gradj = new double**[nTotalVar];
    for (unsigned short iVar = 0; iVar < nTotalVar; iVar++) {
        Jacobian_gradi[iVar] = new double*[nDim];
        Jacobian_gradj[iVar] = new double*[nDim];
        for (unsigned short jVar = 0; jVar < nDim; jVar++) {
            Jacobian_gradi[iVar][jVar] = new double[nVar];
            Jacobian_gradj[iVar][jVar] = new double[nVar];
        }
    }
    
    
	/*--- Initialization of the structure of the whole Jacobian ---*/
	Initialize_SparseMatrix_Structure(&Jacobian, nVar, nVar, geometry, config);
    Jacobian.SetValZero();
	xsol = new double [nPoint*nVar];
	rhs  = new double [nPoint*nVar];

	/*--- Initialization of discrete sparse Jacobian for Hybrid ---*/
	// nVar = # turb vars, nTotalVar = # turb vars + # flow vars
	Initialize_SparseMatrix_Structure(&DirectJacobian, nTotalVar, nVar, geometry, config);
	DirectJacobian.SetValZero();

	/*--- Initialization of discrete sparse Jacobian for Hybrid BC ---*/
	// nVar = # turb vars, nTotalVar = # turb vars + # flow vars
	Initialize_SparseMatrix_Structure(&DirectBCJacobian, nTotalVar, nVar, geometry, config);
	DirectBCJacobian.SetValZero();

	/*--- Computation of gradients by least squares ---*/
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim]; 
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nVar+1]; 
		for (iVar = 0; iVar < nVar+1; iVar++)
			cvector[iVar] = new double [nDim];
	}
	
	/*--- Far-Field values and initizalization ---*/
	node = new CVariable* [nPoint];	
	bool restart = config->GetRestart();

	if (!restart || geometry->GetFinestMGLevel() == false) {	
		PsiNu_Inf = 0.0;
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			//nNeigh = geometry->node[iPoint]->GetnNeighbor();
			node[iPoint] = new CAdjTurbVariable(PsiNu_Inf, nDim, nVar, config);
		}
	}
	else {
		unsigned long index;
		double dull_val;
		string filename, AdjExt, text_line;
    ifstream restart_file;

    /*--- Restart the solution from file information ---*/
		string mesh_filename = config->GetSolution_AdjFileName();
		
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
    long iPoint_Local; unsigned long iPoint_Global = 0;
    
    /*--- The first line is the header ---*/
    getline (restart_file, text_line);
    
    while (getline (restart_file,text_line)) {
			istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives 
       on a different processor, the value of iPoint_Local will be -1. 
       Otherwise, the local index for this node on the current processor 
       will be returned and used to instantiate the vars. ---*/
      iPoint_Local = Global2Local[iPoint_Global];
      if (iPoint_Local >= 0) {
        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        //nNeigh = geometry->node[iPoint_Local]->GetnNeighbor();
        node[iPoint_Local] = new CAdjTurbVariable(Solution[0], nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
    	//nNeigh = geometry->node[iPoint_Local]->GetnNeighbor();
        node[iPoint] = new CAdjTurbVariable(Solution[0], nDim, nVar, config);
    }
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
	}
  
}

CAdjTurbSolution::~CAdjTurbSolution(void) {

    unsigned short nTotalVar = nVar + (nDim + 2); //nTotalVar = # turb vars + # flow vars
    
    for (unsigned short iVar = 0; iVar < nTotalVar; iVar++) {
        delete [] Jacobian_i[iVar];
        delete [] Jacobian_j[iVar];
    }
    
    delete [] Jacobian_i;
    delete [] Jacobian_j;
    delete [] Jacobian_mui;
    delete [] Jacobian_muj;
    
    for (unsigned short iVar = 0; iVar < nTotalVar; iVar++) {
        for (unsigned short jVar = 0; jVar < nDim; jVar++) {
            delete [] Jacobian_gradi[iVar][jVar];
            delete [] Jacobian_gradj[iVar][jVar];
        }
        delete [] Jacobian_gradi[iVar];
        delete [] Jacobian_gradj[iVar];
    }
    delete [] Jacobian_gradi;
    delete [] Jacobian_gradj;
    

}

void CAdjTurbSolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) {

	unsigned long Point, iVertex;
    unsigned short iVar;

	for (iVertex = 0; iVertex<geometry->nVertex[val_marker]; iVertex++) {
		Point = geometry->vertex[val_marker][iVertex]->GetNode();
		Solution[0] = 0.0;
		node[Point]->SetSolution_Old(Solution);
        
        if (config->GetExtIter() == 0) {
            for (iVar = 0; iVar < nVar; iVar++)
                Residual[iVar] = EPS;
            //node[Point]->AddResidual(Residual);
        } else
            node[Point]->SetResidualZero();
//		node[Point]->SetRes_TruncErrorZero();
		Jacobian.DeleteValsRowi(Point); // & includes 1 in the diagonal
	}
}

/* void CAdjTurbSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {

	unsigned long Point, iVertex;
	unsigned short iDim;
	double **Grad_i, *Normal, flux, normal2, Vol_i;

	// This tries to impose a zero-flux BC on the far-field by using an approximation of the Green-Gauss theorem

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		Point = geometry->vertex[val_marker][iVertex]->GetNode();

		Grad_i = node[Point]->GetGradient();
		Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
		Vol_i = geometry->node[Point]->GetVolume();

		flux = 0;
		normal2 = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			normal2 += Normal[iDim]*Normal[iDim];
			flux -= Grad_i[0][iDim]*Normal[iDim];
		}

		Solution[0] = node[Point]->GetSolution(0) - 2.0*Vol_i*flux/normal2;
		node[Point]->SetSolution_Old(Solution);
		node[Point]->SetResidualZero();
//		node[Point]->SetRes_TruncErrorZero();
		Jacobian.DeleteValsRowi(Point); // & includes 1 in the diagonal
	}
}*/

void CAdjTurbSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {

	unsigned long Point, iVertex;
    unsigned short iPos, jPos;
    unsigned short nFlowVar = nDim +2;
	
	if ((config->GetKind_Adjoint() != HYBRID) || (config->GetExtIter() == 0)) {
		for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

			Point = geometry->vertex[val_marker][iVertex]->GetNode();
            
            /*--- Set Normal ---*/
            conv_solver->SetNormal(geometry->vertex[val_marker][iVertex]->GetNormal());

            /*--- Set up discrete system of Hybrid Adjoint --*/
			if (config->GetKind_Adjoint() == HYBRID) {
                
                double *Normal;
                
                Normal = new double[nDim];
                
                FlowSolution_i = solution_container[FLOW_SOL]->node[Point]->GetSolution();
                
                /*--- Construct solution state at infinity (farfield) ---*/
                // Compressible
                FlowSolution_j[0] = solution_container[FLOW_SOL]->GetDensity_Inf();
                FlowSolution_j[nDim+1] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
                for (unsigned short iDim = 0; iDim < nDim; iDim++)
                    FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(iDim);
                
                conv_solver->SetConservative(FlowSolution_i, FlowSolution_j); 
                
                /*--- Set Normal (it is necessary to change the sign) ---*/
                geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
                for (unsigned short iDim = 0; iDim < nDim; iDim++)
                    Normal[iDim] = -Normal[iDim];
                conv_solver->SetNormal(Normal);


				double *Turb_i, Turb_j;
				/*--- Turbulence variable w/o reconstruction ---*/
				Turb_i = solution_container[TURB_SOL]->node[Point]->GetSolution();

				/*--- Read farfield conditions from config ---*/
				double Density_Inf   = config->GetDensity_FreeStreamND();
				double Viscosity_Inf = config->GetViscosity_FreeStreamND();

				/*--- Factor_nu_Inf in [3.0, 5.0] ---*/
				double Factor_nu_Inf = 3.0;
				double nu_tilde_Inf  = Factor_nu_Inf*Viscosity_Inf/Density_Inf;

				Turb_j = nu_tilde_Inf;
				conv_solver->SetTurbVar(Turb_i, &Turb_j);

				// BUILD DISCRETE SYSTEM
                /*--- Auto-Diff direct residual ---*/
				conv_solver->SetResidual(Jacobian_i, Jacobian_j, config);

				/*--- Save contribution from explicit U_i, U_j sensitivity ---*/
				DirectBCJacobian.SubtractBlock(Point, Point, Jacobian_i);
                
                for (iPos = 0; iPos < nVar; iPos++)
                    for (jPos = 0; jPos < nVar; jPos++) {
                        Jacobian_ii[iPos][jPos] = Jacobian_i[iPos+nFlowVar][jPos];
                    }
                
                Jacobian.SubtractBlock(Point,Point,Jacobian_ii);

                delete [] Normal;
                
            /*--- Get Continuous Adjoint Residual --*/
			} else {
                
                /*--- Set Conservative variables (for convection) ---*/
                double* U_i = solution_container[FLOW_SOL]->node[Point]->GetSolution();
                conv_solver->SetConservative(U_i, NULL);

				/*--- Turbulent adjoint variables w/o reconstruction ---*/
				double* TurbPsi_i = node[Point]->GetSolution();
				conv_solver->SetTurbAdjointVar(TurbPsi_i, NULL);

				/*--- Add Residuals and Jacobians ---*/
				conv_solver->SetResidual(Residual, Jacobian_ii, NULL, config);
				node[Point]->AddResidual(Residual);
				Jacobian.AddBlock(Point, Point, Jacobian_ii);
			}

		}
	}
}

void CAdjTurbSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

	for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		node[iPoint]->SetResidualZero(); // Initialize the residual vector
    
    if ((config->GetKind_Adjoint() != HYBRID) || ((config->GetKind_Adjoint() == HYBRID) && (config->GetExtIter() == 0)))
        Jacobian.SetValZero();

	/*--- Computation of gradients of the different variables ---*/
	switch (config->GetKind_Gradient_Method()) {
		case GREEN_GAUSS : 
			SetSolution_Gradient_GG(geometry, config);
			solution_container[ADJFLOW_SOL]->SetSolution_Gradient_GG(geometry, config);
			break;
		case WEIGHTED_LEAST_SQUARES : 
			SetSolution_Gradient_LS(geometry, config);
			solution_container[ADJFLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
			break;
	}

}

void CAdjTurbSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {

	unsigned long iEdge, iPoint, jPoint;
	double *U_i, *U_j, *TurbVar_i, *TurbVar_j, *TurbPsi_i, *TurbPsi_j, **TurbVar_Grad_i, **TurbVar_Grad_j, *Limiter_i = NULL,
    *Limiter_j = NULL, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
	unsigned short iDim, iVar;
    unsigned short iPos, jPos;
    unsigned short nFlowVar = nDim +2;
    
	bool high_order_diss = (config->GetKind_Upwind_AdjTurb() == SCALAR_UPWIND_2ND);
	bool limiter = (config->GetKind_SlopeLimit() != NONE);

	if (high_order_diss) { 
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
		if (limiter) SetSolution_Limiter(geometry, config);
	}

	if ((config->GetKind_Adjoint() != HYBRID) || (config->GetExtIter() == 0)) {
		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

            /*--- Initialise flow conditions and geometric info ---*/
			/*--- Points in edge ---*/
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);

			/*--- Conservative variables w/o reconstruction ---*/
			U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
			U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
			solver->SetConservative(U_i, U_j);
            
            /*--- Set normal vectors and length ---*/
            solver->SetNormal(geometry->edge[iEdge]->GetNormal());

			/*--- Set up discrete system of Hybrid Adjoint --*/
			if (config->GetKind_Adjoint() == HYBRID) {

				/*--- Turbulence variable w/o reconstruction ---*/
				TurbVar_i = solution_container[TURB_SOL]->node[iPoint]->GetSolution();
				TurbVar_j = solution_container[TURB_SOL]->node[jPoint]->GetSolution();
				solver->SetTurbVar(TurbVar_i, TurbVar_j);

				// BUILD DISCRETE SYSTEM
                /*--- Auto-Diff direct residual ---*/
				solver->SetResidual(Jacobian_i, Jacobian_j, config);

				/*--- Save contribution from explicit U_i, U_j sensitivity ---*/
				DirectJacobian.AddBlock(iPoint, iPoint, Jacobian_i);
				DirectJacobian.SubtractBlock(iPoint, jPoint, Jacobian_i);
				DirectJacobian.AddBlock(jPoint, iPoint, Jacobian_j);
				DirectJacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
                
                for (iPos = 0; iPos < nVar; iPos++)
                    for (jPos = 0; jPos < nVar; jPos++) {
                        Jacobian_ii[iPos][jPos] = Jacobian_i[iPos+nFlowVar][jPos];
                        Jacobian_jj[iPos][jPos] = Jacobian_j[iPos+nFlowVar][jPos];
                    }
                        
                Jacobian.AddBlock(iPoint,iPoint,Jacobian_ii);
				Jacobian.SubtractBlock(iPoint,jPoint,Jacobian_ii);
				Jacobian.AddBlock(jPoint,iPoint,Jacobian_jj);
				Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_jj);

            /*--- Get Continuous Adjoint Residual --*/
			} else {

				/*--- Turbulent adjoint variables w/o reconstruction ---*/
				TurbPsi_i = node[iPoint]->GetSolution();
				TurbPsi_j = node[jPoint]->GetSolution();
				solver->SetTurbAdjointVar(TurbPsi_i, TurbPsi_j);

				/*--- Gradient of turbulent variables w/o reconstruction ---*/
				TurbVar_Grad_i = solution_container[TURB_SOL]->node[iPoint]->GetGradient();
				TurbVar_Grad_j = solution_container[TURB_SOL]->node[jPoint]->GetGradient();
				solver->SetTurbVarGradient(TurbVar_Grad_i, TurbVar_Grad_j);

				if (high_order_diss) {

					/*--- Conservative solution using gradient reconstruction ---*/
					for (iDim = 0; iDim < nDim; iDim++) {
						Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
						Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
					}
					Gradient_i = solution_container[FLOW_SOL]->node[iPoint]->GetGradient();
					Gradient_j = solution_container[FLOW_SOL]->node[jPoint]->GetGradient();
					for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
						Project_Grad_i = 0; Project_Grad_j = 0;
						for (iDim = 0; iDim < nDim; iDim++) {
							Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
							Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
						}
						FlowSolution_i[iVar] = U_i[iVar] + Project_Grad_i;
						FlowSolution_j[iVar] = U_j[iVar] + Project_Grad_j;
					}
					solver->SetConservative(FlowSolution_i, FlowSolution_j);

					/*--- Adjoint turbulent variables using gradient reconstruction ---*/
					Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
					if (limiter) { Limiter_i = node[iPoint]->GetLimiter(); Limiter_j = node[jPoint]->GetLimiter(); }
					for (iVar = 0; iVar < nVar; iVar++) {
						Project_Grad_i = 0; Project_Grad_j = 0;
						for (iDim = 0; iDim < nDim; iDim++) {
							Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
							Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
						}
						if (limiter) {
							Solution_i[iVar] = TurbPsi_i[iVar] + Project_Grad_i*Limiter_i[iVar];
							Solution_j[iVar] = TurbPsi_j[iVar] + Project_Grad_j*Limiter_j[iVar];
						}
						else {
							Solution_i[iVar] = TurbPsi_i[iVar] + Project_Grad_i;
							Solution_j[iVar] = TurbPsi_j[iVar] + Project_Grad_j;
						}
					}
					solver->SetTurbVar(Solution_i, Solution_j);
				}

				/*--- Set normal vectors and length ---*/
				solver->SetNormal(geometry->edge[iEdge]->GetNormal());

				solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

				/*--- Add and Subtract Residual ---*/
				node[iPoint]->AddResidual(Residual_i);
				node[jPoint]->AddResidual(Residual_j);
				Jacobian.AddBlock(iPoint,iPoint,Jacobian_ii);
				Jacobian.AddBlock(iPoint,jPoint,Jacobian_ij);
				Jacobian.AddBlock(jPoint,iPoint,Jacobian_ji);
				Jacobian.AddBlock(jPoint,jPoint,Jacobian_jj);
			}
		}

	}

}

void CAdjTurbSolution::Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, 
																				unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint, kPoint;
	double *Coord_i, *Coord_j;
	bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);

	if ((config->Get_Beta_RKStep(iRKStep) != 0) || implicit) {
		for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
			
			/*--- Points in edge ---*/
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);
			
            /*--- Get Continuous Adjoint Residual --*/
			if (config->GetKind_Adjoint() != HYBRID) {

				/*--- Points coordinates, and set normal vectors and length ---*/
				Coord_i = geometry->node[iPoint]->GetCoord();
				Coord_j = geometry->node[jPoint]->GetCoord();
				solver->SetCoord(Coord_i, Coord_j);
				solver->SetNormal(geometry->edge[iEdge]->GetNormal());

				/*--- Conservative variables w/o reconstruction, turbulent variables w/o reconstruction,
			 and turbulent adjoint variables w/o reconstruction ---*/
				solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), solution_container[FLOW_SOL]->node[jPoint]->GetSolution());
				solver->SetTurbVar(solution_container[TURB_SOL]->node[iPoint]->GetSolution(), solution_container[TURB_SOL]->node[jPoint]->GetSolution());
				solver->SetTurbAdjointVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());

				/*--- Viscosity ---*/
				solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
						solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());

				/*--- Turbulent adjoint variables w/o reconstruction ---*/
				solver->SetTurbAdjointGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());

				// ATTENTION: CHOOSE ONE OF THE FOLLOWING FORMS TO COMPUTE THE RESIDUAL

				// Add and Subtract Residual (CONSERVATIVE FORM)
				solver->SetResidual(Residual, Jacobian_ii, Jacobian_jj, config);
				node[iPoint]->AddResidual(Residual);
				node[jPoint]->SubtractResidual(Residual);
				Jacobian.AddBlock(iPoint,iPoint,Jacobian_ii);
				Jacobian.AddBlock(iPoint,jPoint,Jacobian_jj);
				Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_ii);
				Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_jj);

				/*		// Add and Subtract Residual (NON-CONSERVATIVE FORM)
			 solver->SetResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
			 node[iPoint]->AddResidual(Residual_i);
			 node[jPoint]->AddResidual(Residual_j);
			 Jacobian.AddBlock(iPoint,iPoint,Jacobian_ii);
			 Jacobian.AddBlock(iPoint,jPoint,Jacobian_ij);
			 Jacobian.AddBlock(jPoint,iPoint,Jacobian_ji);
			 Jacobian.AddBlock(jPoint,jPoint,Jacobian_jj);*/

            /*--- Set up discrete system of Hybrid Adjoint --*/
			} else if ((config->GetKind_Adjoint() == HYBRID) && (config->GetExtIter() == 0)) {

				unsigned short iPos, jPos;
				double Laminar_Viscosity_i, Laminar_Viscosity_j;
				double *U_i, *U_j, *TurbVar_i, *TurbVar_j;
				unsigned short nFlowVar = nDim + 2;
				unsigned short iNode;
				double **Normals;
                double **Jacobian_ik, **Jacobian_jk;

				/*--- Initialise flow conditions and geometric info ---*/
				/*--- Points coordinates, and set normal vectors and length ---*/
				Coord_i = geometry->node[iPoint]->GetCoord();
				Coord_j = geometry->node[jPoint]->GetCoord();
				solver->SetCoord(Coord_i, Coord_j);
				solver->SetNormal(geometry->edge[iEdge]->GetNormal());

				U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
				U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
				solver->SetConservative(U_i, U_j);
                
                // Gradient of primitive variables w/o reconstruction
                solver->SetPrimVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), solution_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive());

				Laminar_Viscosity_i = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
				Laminar_Viscosity_j = solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();
				solver->SetLaminarViscosity(Laminar_Viscosity_i, Laminar_Viscosity_j);
                

				/*--- Eddy Viscosity ---*/
				solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
						solution_container[FLOW_SOL]->node[jPoint]->GetEddyViscosity());

				/*--- Turbulent variables w/o reconstruction, and its gradients ---*/
				TurbVar_i = solution_container[TURB_SOL]->node[iPoint]->GetSolution();
				TurbVar_j = solution_container[TURB_SOL]->node[jPoint]->GetSolution();
				solver->SetTurbVar(TurbVar_i, TurbVar_j);
				solver->SetTurbVarGradient(solution_container[TURB_SOL]->node[iPoint]->GetGradient(),solution_container[TURB_SOL]->node[jPoint]->GetGradient());

				// BUILD DISCRETE SYSTEM
                /*--- Auto-Diff direct residual ---*/
				solver->SetResidual(Jacobian_i, Jacobian_mui, Jacobian_gradi,
						Jacobian_j, Jacobian_muj, Jacobian_gradj, config);

				/*--- Save contribution from explicit U_i, U_j sensitivity ---*/
				DirectJacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
				DirectJacobian.AddBlock(iPoint, jPoint, Jacobian_i);
				DirectJacobian.SubtractBlock(jPoint, iPoint, Jacobian_j);
				DirectJacobian.AddBlock(jPoint, jPoint, Jacobian_j);
                
                for (iPos = 0; iPos < nVar; iPos++)
                    for (jPos = 0; jPos < nVar; jPos++) {
                        Jacobian_ii[iPos][jPos] = Jacobian_i[iPos+nFlowVar][jPos];
                        Jacobian_jj[iPos][jPos] = Jacobian_j[iPos+nFlowVar][jPos];
                    }
                
                Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_ii);
				Jacobian.AddBlock(iPoint,jPoint,Jacobian_ii);
				Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_jj);
				Jacobian.AddBlock(jPoint,jPoint,Jacobian_jj);

				/*--- Extract contribution from implicit U_i sensitivity (from mui sensitivity) ---*/
				this->ConvertSensMu_to_SensU(U_i, Laminar_Viscosity_i,
						Jacobian_i, Jacobian_mui,
						nFlowVar, nVar+nFlowVar, nVar, config);

				/*--- Save contribution from implicit U_i sensitivity (from mui sensitivity) ---*/
				DirectJacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
				DirectJacobian.AddBlock(iPoint, jPoint, Jacobian_i);
                
                for (iPos = 0; iPos < nVar; iPos++)
                    for (jPos = 0; jPos < nVar; jPos++) {
                        Jacobian_ii[iPos][jPos] = Jacobian_i[iPos+nFlowVar][jPos];
                    }
                
                Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_ii);
				Jacobian.AddBlock(iPoint,jPoint,Jacobian_ii);

				/*--- Extract contribution from implicit U_j sensitivity (from muj sensitivity) ---*/
				this->ConvertSensMu_to_SensU(U_j, Laminar_Viscosity_j,
						Jacobian_j, Jacobian_muj,
						nFlowVar, nVar+nFlowVar, nVar, config);

				/*--- Save contribution from implicit U_j sensitivity (from muj sensitivity) ---*/
				DirectJacobian.SubtractBlock(jPoint, iPoint, Jacobian_j);
				DirectJacobian.AddBlock(jPoint, jPoint, Jacobian_j);
                
                for (iPos = 0; iPos < nVar; iPos++)
                    for (jPos = 0; jPos < nVar; jPos++) {
                        Jacobian_jj[iPos][jPos] = Jacobian_j[iPos+nFlowVar][jPos];
                    }
                
				Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_jj);
				Jacobian.AddBlock(jPoint,jPoint,Jacobian_jj);

				/*--- Extract implicit U_i sensitivity (from gradient_i sensitivity) ---*/
                /*--- Set up arrays and store flow and geometry information ---*/
				unsigned short nNeigh;
				nNeigh = geometry->node[iPoint]->GetnPoint();
				double *U_k, *TurbVar_k, *Normal;
				double *Vars_i, **Vars_ks;
				unsigned long iEdge;
                
                Vars_i = new double[nFlowVar+nVar];
                Vars_ks = new double*[nNeigh];
                Normals = new double*[nNeigh];
                for (iPos = 0; iPos<nNeigh; iPos++) {
                    Vars_ks[iPos] = new double[nFlowVar+nVar];
                    Normals[iPos] = new double[nDim];
                }
                
                Jacobian_ik  = new double*[nVar];
                Jacobian_jk  = new double*[nVar];
                for (iPos = 0; iPos<nVar; iPos++) {
                    Jacobian_ik[iPos] = new double[nVar];
                    Jacobian_jk[iPos] = new double[nVar];
                }
                
                for (iPos = 0; iPos<nFlowVar; iPos++) {
                    Vars_i[iPos] = U_i[iPos];
                }
                
                for (iPos = 0; iPos<nVar; iPos++) {
                    Vars_i[iPos+nFlowVar] = TurbVar_i[iPos];
                }
                
                for (iNode = 0; iNode < nNeigh; iNode++){
                    kPoint = geometry->node[iPoint]->GetPoint(iNode);
                    U_k = solution_container[FLOW_SOL]->node[kPoint]->GetSolution();
                    TurbVar_k = solution_container[TURB_SOL]->node[kPoint]->GetSolution();
                    iEdge = geometry->FindEdge(iPoint,kPoint);
                    Normal = geometry->edge[iEdge]->GetNormal();
                    
                    for (iPos = 0; iPos<nFlowVar; iPos++) {
                        Vars_ks[iNode][iPos] = U_k[iPos];
                    }
                    
                    for (iPos = 0; iPos<nVar; iPos++) {
                        Vars_ks[iNode][iPos+nFlowVar] = TurbVar_k[iPos];
                    }
                
                for (iPos = 0; iPos<nDim; iPos++) {
                    Normals[iNode][iPos] = Normal[iPos];
                }
                }

                /*--- For i ---*/
                /*--- Extract Primitive_i sensitivity (from gradient_i sensitivity) ---*/
				ConvertSensGradPrimVar_to_SensPrimVar(Vars_i, Vars_ks,
						Jacobian_i, Jacobian_gradi, iPoint, -1, Normals,
						nVar+nFlowVar, nNeigh, nVar, geometry, config);

				/*--- Extract implicit U_i sensitivity (from gradient_i -> Primitive_i sensitivity) ---*/
				ConvertSensPrimVar_to_SensU(U_i, Jacobian_i, nFlowVar, config);

				/*--- Save contribution from implicit U_i sensitivity (from gradient_i -> Primitive_i sensitivity) ---*/
				DirectJacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
				DirectJacobian.AddBlock(iPoint, jPoint, Jacobian_i);
                
                for (iPos = 0; iPos < nVar; iPos++)
                    for (jPos = 0; jPos < nVar; jPos++) {
                        Jacobian_ii[iPos][jPos] = Jacobian_i[iPos+nFlowVar][jPos];
                    }
                
                Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_ii);
				Jacobian.AddBlock(iPoint,jPoint,Jacobian_ii);

				/*--- For neighbours (k) of i ---*/
				for (iNode = 0; iNode < nNeigh; iNode++){
					kPoint = geometry->node[iPoint]->GetPoint(iNode);

                    /*--- Extract Primitive_k sensitivity (from gradient_i sensitivity) ---*/
						ConvertSensGradPrimVar_to_SensPrimVar(Vars_i, Vars_ks,
								Jacobian_i, Jacobian_gradi, iPoint, iNode, Normals,
								nVar+nFlowVar, nNeigh, nVar, geometry, config);
                        
                        for (iPos = 0; iPos<nFlowVar; iPos++) {
                            U_k[iPos] = Vars_ks[iNode][iPos];
                        }

						/*--- Extract implicit U_k sensitivity (from gradient_i -> Primitive_k sensitivity) ---*/
						ConvertSensPrimVar_to_SensU(U_k, Jacobian_i, nFlowVar, config);

						/*--- Save contribution from implicit U_k sensitivity (from gradient_i -> Primitive_k sensitivity) ---*/
                    if (geometry->CheckEdge(iPoint, kPoint)) {
						DirectJacobian.SubtractBlock(kPoint, iPoint, Jacobian_i);
                        
                        for (iPos = 0; iPos < nVar; iPos++)
                            for (jPos = 0; jPos < nVar; jPos++) {
                                Jacobian_ik[iPos][jPos] = Jacobian_i[iPos+nFlowVar][jPos];
                            }
                        
                        Jacobian.SubtractBlock(kPoint,iPoint,Jacobian_ik);

                    }
                    if (geometry->CheckEdge(jPoint, kPoint)) {
						DirectJacobian.AddBlock(kPoint, jPoint, Jacobian_i);
                    
                        for (iPos = 0; iPos < nVar; iPos++)
                            for (jPos = 0; jPos < nVar; jPos++) {
                                Jacobian_jk[iPos][jPos] = Jacobian_i[iPos+nFlowVar][jPos];
                            }
                    
                        Jacobian.AddBlock(kPoint,jPoint,Jacobian_jk);
                    }
                


				}

                for (iPos = 0; iPos<nNeigh; iPos++) {
                    delete [] Vars_ks[iPos];
                    delete [] Normals[iPos];
                }
                
                delete [] Vars_ks;
                delete [] Normals;

				/*--- Extract implicit U_j sensitivity (from gradient_j sensitivity) ---*/
                /*--- Store flow and geometry information ---*/
				nNeigh = geometry->node[jPoint]->GetnPoint();

                Vars_i = new double[nFlowVar+nVar];
                Vars_ks = new double*[nNeigh];
                Normals = new double*[nNeigh];
                for (iPos = 0; iPos<nNeigh; iPos++) {
                    Vars_ks[iPos] = new double[nFlowVar+nVar];
                    Normals[iPos] = new double[nDim];
                }
                
                
                for (iPos = 0; iPos<nFlowVar; iPos++) {
                    Vars_i[iPos] = U_j[iPos];
                }
                
                for (iPos = 0; iPos<nVar; iPos++) {
                    Vars_i[iPos+nFlowVar] = TurbVar_j[iPos];
                }


                for (iNode = 0; iNode < nNeigh; iNode++){
                    kPoint = geometry->node[jPoint]->GetPoint(iNode);
                    U_k = solution_container[FLOW_SOL]->node[kPoint]->GetSolution();
                    TurbVar_k = solution_container[TURB_SOL]->node[kPoint]->GetSolution();
                    iEdge = geometry->FindEdge(jPoint,kPoint);
                    Normal = geometry->edge[iEdge]->GetNormal();
                    
                    for (iPos = 0; iPos<nFlowVar; iPos++) {
                        Vars_ks[iNode][iPos] = U_k[iPos];
                    }
                    
                    for (iPos = 0; iPos<nVar; iPos++) {
                        Vars_ks[iNode][iPos+nFlowVar] = TurbVar_k[iPos];
                    }
                
                for (iPos = 0; iPos<nDim; iPos++) {
                    Normals[iNode][iPos] = Normal[iPos];
                }
                }

                /*--- For j ---*/
                /*--- Extract Primitive_j sensitivity (from gradient_j sensitivity) ---*/
				ConvertSensGradPrimVar_to_SensPrimVar(Vars_i, Vars_ks,
						Jacobian_j, Jacobian_gradj, jPoint, -1, Normals,
						nVar+nFlowVar, nNeigh, nVar, geometry, config);

				/*--- Extract implicit U_j sensitivity (from gradient_j -> Primitive_j sensitivity) ---*/
				ConvertSensPrimVar_to_SensU(U_j, Jacobian_j, nFlowVar, config);

				/*--- Save contribution from implicit U_j sensitivity (from gradient_j -> Primitive_j sensitivity) ---*/
				DirectJacobian.SubtractBlock(jPoint, iPoint, Jacobian_j);
				DirectJacobian.AddBlock(jPoint, jPoint, Jacobian_j);
                
                for (iPos = 0; iPos < nVar; iPos++)
                    for (jPos = 0; jPos < nVar; jPos++) {
                        Jacobian_jj[iPos][jPos] = Jacobian_j[iPos+nFlowVar][jPos];
                    }
                
				Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_jj);
				Jacobian.AddBlock(jPoint,jPoint,Jacobian_jj);

				/*--- For neighbours (k) of j ---*/
				for (iNode = 0; iNode < nNeigh; iNode++){
					kPoint = geometry->node[jPoint]->GetPoint(iNode);

                    /*--- Extract Primitive_k sensitivity (from gradient_j sensitivity) ---*/
						ConvertSensGradPrimVar_to_SensPrimVar(Vars_i, Vars_ks,
								Jacobian_j, Jacobian_gradj, jPoint, iNode, Normals,
								nVar+nFlowVar, nNeigh, nVar, geometry, config);
                        
                        for (iPos = 0; iPos<nFlowVar; iPos++) {
                            U_k[iPos] = Vars_ks[iNode][iPos];
                        }

						/*--- Extract implicit U_k sensitivity (from gradient_j -> Primitive_k sensitivity) ---*/
						ConvertSensPrimVar_to_SensU(U_k, Jacobian_j, nFlowVar, config);

						/*--- Save contribution from implicit U_k sensitivity (from gradient_j -> Primitive_k sensitivity) ---*/
                    if (geometry->CheckEdge(iPoint, kPoint)) {
						DirectJacobian.SubtractBlock(kPoint, iPoint, Jacobian_i);
                        
                        for (iPos = 0; iPos < nVar; iPos++)
                            for (jPos = 0; jPos < nVar; jPos++) {
                                Jacobian_ik[iPos][jPos] = Jacobian_j[iPos+nFlowVar][jPos];
                            }
                        
                        Jacobian.SubtractBlock(kPoint,iPoint,Jacobian_ik);
                        
                    }
                    if (geometry->CheckEdge(jPoint, kPoint)) {
						DirectJacobian.AddBlock(kPoint, jPoint, Jacobian_i);
                        
                        for (iPos = 0; iPos < nVar; iPos++)
                            for (jPos = 0; jPos < nVar; jPos++) {
                                Jacobian_jk[iPos][jPos] = Jacobian_j[iPos+nFlowVar][jPos];
                            }
                        
                        Jacobian.AddBlock(kPoint,jPoint,Jacobian_jk);
                    }

				}
                
                delete [] Vars_i;
                
                for (iPos = 0; iPos<nNeigh; iPos++) {
                    delete [] Vars_ks[iPos];
                    delete [] Normals[iPos];
                }
                
                delete [] Vars_ks;
                delete [] Normals;
                
                for (iPos = 0; iPos<nVar; iPos++) {
                    delete [] Jacobian_ik[iPos];
                    delete [] Jacobian_jk[iPos];
                }
                
                delete [] Jacobian_ik;
                delete [] Jacobian_jk;


			}

            /*--- Get Hybrid Adjoint Residual --*/
			if (config->GetKind_Adjoint() == HYBRID) {

				double *TurbPsi_i, *TurbPsi_j;
				double **DJ_ij, **DJ_ji;

				unsigned short nFlowVar, nTurbVar, nTotalVar;
				nFlowVar = nDim + 2;
				nTurbVar = nVar;
				nTotalVar = nFlowVar + nTurbVar;

				DJ_ij = new double*[nTotalVar];
				DJ_ji = new double*[nTotalVar];
				for (unsigned short iVar = 0; iVar<nTotalVar; iVar++){
					DJ_ij[iVar] = new double[nTurbVar];
					DJ_ji[iVar] = new double[nTurbVar];
				}
                
                for (unsigned short iVar = 0; iVar<nVar; iVar++){
                    
                    Residual_i[iVar] = 0.0;
                    Residual_j[iVar] = 0.0;
                    
				}


				TurbPsi_i = solution_container[ADJTURB_SOL]->node[iPoint]->GetSolution();
				TurbPsi_j = solution_container[ADJTURB_SOL]->node[jPoint]->GetSolution();

				DirectJacobian.GetBlock(iPoint, jPoint);
				DirectJacobian.ReturnBlock(DJ_ij);
				DirectJacobian.GetBlock(jPoint, iPoint);
				DirectJacobian.ReturnBlock(DJ_ji);

				for (unsigned short iVar = 0; iVar<nVar; iVar++){

					for (unsigned short jVar = 0; jVar<nVar; jVar++) {
						Residual_i[iVar] += DJ_ij[iVar+nFlowVar][jVar]*TurbPsi_j[jVar]; // +?
						Residual_j[iVar] += DJ_ji[iVar+nFlowVar][jVar]*TurbPsi_i[jVar]; // +?
					}

				}
                
                node[iPoint]->AddResidual(Residual_i);
                node[jPoint]->AddResidual(Residual_j);

				for (unsigned short iVar = 0; iVar<nTotalVar; iVar++){
					delete [] DJ_ij[iVar];
					delete [] DJ_ji[iVar];
				}
				delete [] DJ_ij;
				delete [] DJ_ji;

			}
		}
	}

}

void CAdjTurbSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CNumerics *second_solver,
											   CConfig *config, unsigned short iMesh) {
	unsigned long iPoint, jPoint, kPoint, iEdge;
	unsigned short iPos, jPos;
	double *U_i, **GradPrimVar_i, *TurbVar_i;
	double **TurbVar_Grad_i, **TurbVar_Grad_j, *TurbPsi_i, *TurbPsi_j, **PsiVar_Grad_i; // Gradients

    /*--- Piecewise source term ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 

        /*--- Get Continuous Adjoint Residual --*/
		if (config->GetKind_Adjoint() != HYBRID) {

			// Conservative variables w/o reconstruction
			U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
			solver->SetConservative(U_i, NULL);

			// Gradient of primitive variables w/o reconstruction
			GradPrimVar_i = solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
			solver->SetPrimVarGradient(GradPrimVar_i, NULL);

			// Laminar viscosity of the fluid
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);

			// Turbulent variables w/o reconstruction
			TurbVar_i = solution_container[TURB_SOL]->node[iPoint]->GetSolution();
			solver->SetTurbVar(TurbVar_i, NULL);

			// Gradient of Turbulent Variables w/o reconstruction
			TurbVar_Grad_i = solution_container[TURB_SOL]->node[iPoint]->GetGradient();
			solver->SetTurbVarGradient(TurbVar_Grad_i, NULL);

			// Turbulent adjoint variables w/o reconstruction
			TurbPsi_i = node[iPoint]->GetSolution();
			solver->SetTurbAdjointVar(TurbPsi_i, NULL);

			// Gradient of Adjoint flow variables w/o reconstruction
			// (for non-conservative terms depending on gradients of flow adjoint vars.)
			PsiVar_Grad_i = solution_container[ADJFLOW_SOL]->node[iPoint]->GetGradient();
			solver->SetAdjointVarGradient(PsiVar_Grad_i, NULL);

			// Set volume and distances to the surface
			solver->SetVolume(geometry->node[iPoint]->GetVolume());
			solver->SetDistance(geometry->node[iPoint]->GetWallDistance(), 0.0);

			// Add and Subtract Residual
			solver->SetResidual(Residual, Jacobian_ii, NULL, config);
			node[iPoint]->AddResidual(Residual);
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
            
        /*--- Set up discrete system of Hybrid Adjoint --*/
		} else if ((config->GetKind_Adjoint() == HYBRID) && (config->GetExtIter() == 0)) {
			double *EddyViscSens, *TurbVar_i;
			unsigned short nFlowVar = nDim + 2;
			double Laminar_Viscosity_i;
			unsigned short iNode;
            unsigned short nTotalVar = nDim + 2 + nVar;

			EddyViscSens = new double[nVar+nFlowVar];

			/*--- Initialise flow conditions and geometric info ---*/
			/*--- Conservative variables w/o reconstruction ---*/
			U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
			solver->SetConservative(U_i, NULL);

			/*--- Gradient of the primitive and conservative variables ---*/
			GradPrimVar_i = solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
			solver->SetPrimVarGradient(GradPrimVar_i, NULL);


			Laminar_Viscosity_i = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
			solver->SetLaminarViscosity(Laminar_Viscosity_i, 0.0);

			/*--- Turbulent variables w/o reconstruction, and its gradient ---*/
            TurbVar_i = solution_container[TURB_SOL]->node[iPoint]->GetSolution();
			solver->SetTurbVar(TurbVar_i, NULL);
			solver->SetTurbVarGradient(solution_container[TURB_SOL]->node[iPoint]->GetGradient(), NULL);

			/*--- Set volume ---*/
			solver->SetVolume(geometry->node[iPoint]->GetVolume());

			/*--- Set distance to the surface ---*/
			solver->SetDistance(geometry->node[iPoint]->GetWallDistance(), 0.0);

            // BUILD DISCRETE SYSTEM
			/*--- Auto-Diff direct residual ---*/
			solver->SetResidual(Jacobian_i, Jacobian_mui, Jacobian_gradi, config);

			/*--- Save contribution from explicit U_i sensitivity ---*/
			DirectJacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
            
            for (iPos = 0; iPos < nVar; iPos++)
                for (jPos = 0; jPos < nVar; jPos++) {
                    Jacobian_ii[iPos][jPos] = Jacobian_i[iPos+nFlowVar][jPos];
                }
            
            Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_ii);


			/*--- Extract implicit U_i sensitivity (from mui sensitivity) ---*/
			this->ConvertSensMu_to_SensU(U_i, Laminar_Viscosity_i,
					Jacobian_i, Jacobian_mui,
					nFlowVar, nVar+nFlowVar, nVar, config);

			/*--- Save contribution from implicit U_i sensitivity (from mui sensitivity) ---*/
			DirectJacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
            
            for (iPos = 0; iPos < nVar; iPos++)
                for (jPos = 0; jPos < nVar; jPos++) {
                    Jacobian_ii[iPos][jPos] = Jacobian_i[iPos+nFlowVar][jPos];
                }
            
            Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_ii);
            
            /*--- Extract implicit U_i sensitivity (from gradient_i sensitivity) ---*/
            /*--- Set up arrays and store flow and geometry information ---*/
			unsigned short nNeigh;
			nNeigh = geometry->node[iPoint]->GetnPoint();
			double *U_k, *TurbVar_k, *Normal;
			double *Vars_i, **Vars_ks, **Normals;
			unsigned long iEdge;
            
            Vars_i = new double[nFlowVar+nVar];
            Vars_ks = new double*[nNeigh];
            Normals = new double*[nNeigh];
            for (iPos = 0; iPos<nNeigh; iPos++) {
                Vars_ks[iPos] = new double[nFlowVar+nVar];
                Normals[iPos] = new double[nDim];
            }
            
            for (iPos = 0; iPos<nFlowVar; iPos++) {
                Vars_i[iPos] = U_i[iPos];
            }
            
            for (iPos = 0; iPos<nVar; iPos++) {
                Vars_i[iPos+nFlowVar] = TurbVar_i[iPos];
            }

			for (iNode = 0; iNode < nNeigh; iNode++){
				kPoint = geometry->node[iPoint]->GetPoint(iNode);
                U_k = solution_container[FLOW_SOL]->node[kPoint]->GetSolution();
                TurbVar_k = solution_container[TURB_SOL]->node[kPoint]->GetSolution();
                iEdge = geometry->FindEdge(iPoint,kPoint);
                Normal = geometry->edge[iEdge]->GetNormal();
                
                for (iPos = 0; iPos<nFlowVar; iPos++) {
                    Vars_ks[iNode][iPos] = U_k[iPos];
                }
                
                for (iPos = 0; iPos<nVar; iPos++) {
                    Vars_ks[iNode][iPos+nFlowVar] = TurbVar_k[iPos];
                }
                                
                for (iPos = 0; iPos<nDim; iPos++) {
                    Normals[iNode][iPos] = Normal[iPos];
                }
                
			}


            /*--- For i ---*/
            /*--- Extract Primitive_i sensitivity (from gradient_i sensitivity) ---*/
			ConvertSensGradPrimVar_to_SensPrimVar(Vars_i, Vars_ks,
					Jacobian_i, Jacobian_gradi, iPoint, -1, Normals,
					nVar+nFlowVar, nNeigh, nVar, geometry, config);
            
            /*--- Extract implicit U_i sensitivity (from gradient_i -> Primitive_i sensitivity) ---*/
			ConvertSensPrimVar_to_SensU(U_i, Jacobian_i, nFlowVar, config);

			/*--- Save contribution from implicit U_i sensitivity (from gradient_i -> Primitive_i sensitivity) ---*/
			DirectJacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
            
            for (iPos = 0; iPos < nVar; iPos++)
                for (jPos = 0; jPos < nVar; jPos++) {
                    Jacobian_ii[iPos][jPos] = Jacobian_i[iPos+nFlowVar][jPos];
                }
            
            Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_ii);

			/*--- For neighbours (k) of i ---*/
			for (iNode = 0; iNode < nNeigh; iNode++){
				kPoint = geometry->node[iPoint]->GetPoint(iNode);

                /*--- Extract Primitive_k sensitivity (from gradient_i sensitivity) ---*/
                ConvertSensGradPrimVar_to_SensPrimVar(Vars_i, Vars_ks,
							Jacobian_i, Jacobian_gradi, iPoint, iNode, Normals,
							nVar+nFlowVar, nNeigh, nVar, geometry, config);
                    
                for (iPos = 0; iPos<nFlowVar; iPos++) {
                    U_k[iPos] = Vars_ks[iNode][iPos];
                }

                /*--- Extract implicit U_k sensitivity (from gradient_i -> Primitive_k sensitivity) ---*/
                ConvertSensPrimVar_to_SensU(U_k, Jacobian_i, nFlowVar, config);

                /*--- Save contribution from implicit U_k sensitivity (from gradient_i -> Primitive_k sensitivity) ---*/
                DirectJacobian.SubtractBlock(kPoint, iPoint, Jacobian_i);
                
                for (iPos = 0; iPos < nVar; iPos++)
                    for (jPos = 0; jPos < nVar; jPos++) {
                        Jacobian_ij[iPos][jPos] = Jacobian_i[iPos+nFlowVar][jPos];
                    }
                
                Jacobian.SubtractBlock(kPoint,iPoint,Jacobian_ij);

			}


			// GET EDDY VISC SENS
			/*--- Extract eddy viscosity ---*/
			CalcEddyViscSens(U_i, Laminar_Viscosity_i,
					TurbVar_i, EddyViscSens, nFlowVar, config);
            
            /*--- Store eddy viscosity ---*/
            solution_container[ADJTURB_SOL]->node[iPoint]->SetEddyViscSens(EddyViscSens, nTotalVar);

			delete [] EddyViscSens;
            
            delete [] Vars_i;

            for (iPos = 0; iPos<nNeigh; iPos++) {
                delete [] Vars_ks[iPos];
                delete [] Normals[iPos];
            }
            
            delete [] Vars_ks;
            delete [] Normals;
            

		}

		/*--- Get Hybrid Adjoint Residual --*/
		if (config->GetKind_Adjoint() == HYBRID){
			double *EddyViscSens, kappapsi_Volume;
			double **DJ_ii, **DBCJ_ii;

			unsigned short nFlowVar, nTurbVar, nTotalVar;
			nFlowVar = nDim + 2;
			nTurbVar = nVar;
			nTotalVar = nFlowVar + nTurbVar;

			DJ_ii = new double*[nTotalVar];
			DBCJ_ii = new double*[nTotalVar];
			for (unsigned short iVar = 0; iVar<nTotalVar; iVar++){
				DJ_ii[iVar] = new double[nTurbVar];
				DBCJ_ii[iVar] = new double[nTurbVar];
			}
            
            for (unsigned short iVar = 0; iVar<nVar; iVar++){
                
				Residual[iVar] = 0.0;
                
			}

			TurbPsi_i = solution_container[ADJTURB_SOL]->node[iPoint]->GetSolution();

			DirectJacobian.GetBlock(iPoint, iPoint);
			DirectJacobian.ReturnBlock(DJ_ii);

			DirectBCJacobian.GetBlock(iPoint, iPoint);
			DirectBCJacobian.ReturnBlock(DBCJ_ii);

			EddyViscSens = solution_container[ADJTURB_SOL]->node[iPoint]->GetEddyViscSens();
			kappapsi_Volume = solution_container[ADJFLOW_SOL]->node[iPoint]->GetKappaPsiVolume();

			for (unsigned short iVar = 0; iVar<nVar; iVar++){

				for (unsigned short jVar = 0; jVar<nVar; jVar++) {
					Residual[iVar] += DJ_ii[iVar + nFlowVar][jVar]*TurbPsi_i[jVar]; // -?
					Residual[iVar] += DBCJ_ii[iVar + nFlowVar][jVar]*TurbPsi_i[jVar]; // -?
				}

				Residual[iVar] -= EddyViscSens[iVar+nFlowVar]*kappapsi_Volume; // +?

			}
            
            node[iPoint]->AddResidual(Residual);

			for (unsigned short iVar = 0; iVar<nTotalVar; iVar++){
				delete [] DJ_ii[iVar];
				delete [] DBCJ_ii[iVar];
			}
			delete [] DJ_ii;
			delete [] DBCJ_ii;


		}
	}
    
    /*--- Conservative Source Term ---*/
    if ((config->GetKind_Solver() == ADJ_RANS) && (config->GetKind_Adjoint() == CONTINUOUS)) {
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		// Points in edge
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
        
		// Gradient of turbulent variables w/o reconstruction
		TurbVar_Grad_i = solution_container[TURB_SOL]->node[iPoint]->GetGradient();
		TurbVar_Grad_j = solution_container[TURB_SOL]->node[jPoint]->GetGradient();
		solver->SetTurbVarGradient(TurbVar_Grad_i, TurbVar_Grad_j);
        
		// Turbulent adjoint variables w/o reconstruction
		TurbPsi_i = node[iPoint]->GetSolution();
		TurbPsi_j = node[jPoint]->GetSolution();
		solver->SetTurbAdjointVar(TurbPsi_i, TurbPsi_j);
        
		// Set normal vectors and length
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());
        
		// Add and Subtract Residual
		solver->SetResidual(Residual, Jacobian_ii, Jacobian_jj, config);
		node[iPoint]->AddResidual(Residual);
		node[jPoint]->SubtractResidual(Residual);
		Jacobian.AddBlock(iPoint,iPoint,Jacobian_ii);
		Jacobian.AddBlock(iPoint,jPoint,Jacobian_jj);
		Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_ii);
		Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_jj);
        
	}
    }
    
    
}

void CAdjTurbSolution::Source_Template(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
											   CConfig *config, unsigned short iMesh) {
}

void CAdjTurbSolution::SourceConserv_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
											   CConfig *config, unsigned short iMesh) {
    return;

//	unsigned long iEdge, iPoint, jPoint;
//	double *TurbPsi_i, *TurbPsi_j;
//	double **TurbVar_Grad_i, **TurbVar_Grad_j;
//    
//	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
//		// Points in edge
//		iPoint = geometry->edge[iEdge]->GetNode(0);
//		jPoint = geometry->edge[iEdge]->GetNode(1);
//
//		// Gradient of turbulent variables w/o reconstruction
//		TurbVar_Grad_i = solution_container[TURB_SOL]->node[iPoint]->GetGradient();
//		TurbVar_Grad_j = solution_container[TURB_SOL]->node[jPoint]->GetGradient();
//		solver->SetTurbVarGradient(TurbVar_Grad_i, TurbVar_Grad_j);
//
//		// Turbulent adjoint variables w/o reconstruction
//		TurbPsi_i = node[iPoint]->GetSolution();
//		TurbPsi_j = node[jPoint]->GetSolution();
//		solver->SetTurbAdjointVar(TurbPsi_i, TurbPsi_j);
//
//		// Set normal vectors and length
//		solver->SetNormal(geometry->edge[iEdge]->GetNormal());
//
//		// Add and Subtract Residual
//		solver->SetResidual(Residual, Jacobian_ii, Jacobian_jj, config);
//		node[iPoint]->AddResidual(Residual);
//		node[jPoint]->SubtractResidual(Residual);
//		Jacobian.AddBlock(iPoint,iPoint,Jacobian_ii);
//		Jacobian.AddBlock(iPoint,jPoint,Jacobian_jj);
//		Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_ii);
//		Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_jj);
//	}
}

void CAdjTurbSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta, Res, *local_Residual, Vol;

	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		local_Residual = node[iPoint]->GetResidual();
    
    /*--- Read the volume ---*/
		Vol = geometry->node[iPoint]->GetVolume();
    
		/*--- Modify matrix diagonal to assure diagonal dominance ---*/
		Delta = Vol/(solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() + EPS);
    
		Jacobian.AddVal2Diag(iPoint,Delta);
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
      Res = local_Residual[iVar];
			rhs[total_index] = -Res;
			xsol[total_index] = 0.0;
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
  
	/*--- Solve the system ---*/
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
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++)
			node[iPoint]->AddSolution(iVar,xsol[iPoint*nVar+iVar]);
	}
    
  /*--- MPI solution ---*/
  SetSolution_MPI(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CAdjTurbSolution::ConvertSensMu_to_SensU(double *val_Vars, double val_Laminar_Viscosity,
		double **val_Jacobian, double *val_Jacobian_mu,
		unsigned short numVar, unsigned short numTotalVar, unsigned short numEqn, CConfig *config){

	unsigned short iPos, jPos;
	double Laminar_Viscosityd, *Varsd;

	Varsd = new double[numVar];

	// zero val_Jacobian (reusing structure)
	for (iPos = 0; iPos < numTotalVar; iPos++) {
		for (jPos = 0; jPos < numEqn; jPos++) {

			val_Jacobian[iPos][jPos] = 0.0;

		}
	}
    
//    // TEST *****************
//    for (iPos = 0; iPos < numVar; iPos++){
//        cout << "Flow: " << val_Vars[iPos] << endl;
//    }
//    // *****************

	for (iPos = 0; iPos < numVar; iPos++){
		// zero things first
		for (jPos = 0; jPos < numVar; jPos++){
			Varsd[jPos] = 0.0;
		}

		Laminar_Viscosityd = 0.0;

		Varsd[iPos] = 1.0;

		this->CalcLaminarViscosity_ad(val_Vars, Varsd,
				&val_Laminar_Viscosity, &Laminar_Viscosityd, config);

		for (jPos = 0; jPos < nVar; jPos++) {
			val_Jacobian[iPos][jPos] = Laminar_Viscosityd * val_Jacobian_mu[jPos];
		}

//        // TEST AD *****************
//        cout << "AD: " << Laminar_Viscosityd << endl;
//        // AD *****************
        
	}
    
//    // TEST FD *****************
//    double *temp_val_Vars, temp1_val_Laminar_Viscosity, temp2_val_Laminar_Viscosity;
//    double delta = 0.01;
//    temp_val_Vars = new double[numVar];
//    
//    for (iPos = 0; iPos < numVar; iPos++){
//        
//        for (jPos = 0; jPos < numVar; jPos++){
//			temp_val_Vars[jPos] = val_Vars[jPos];
//		}
//        
//        //temp_val_Vars[iPos] = (1.0-delta) * val_Vars[iPos];
//        temp_val_Vars[iPos] = val_Vars[iPos] - delta;
//    
//        this->CalcLaminarViscosity_ad(temp_val_Vars, Varsd,
//                                  &temp1_val_Laminar_Viscosity, &Laminar_Viscosityd, config);
//        
//        //temp_val_Vars[iPos] = (1.0+delta) * val_Vars[iPos];
//        temp_val_Vars[iPos] = val_Vars[iPos] + delta;
//
//        
//        this->CalcLaminarViscosity_ad(temp_val_Vars, Varsd,
//                                      &temp2_val_Laminar_Viscosity, &Laminar_Viscosityd, config);
//        
//        //cout << "FD: " << (temp2_val_Laminar_Viscosity - temp1_val_Laminar_Viscosity)/(2.0*delta*val_Vars[iPos]) << endl;
//              cout << "FD: " << (temp2_val_Laminar_Viscosity - temp1_val_Laminar_Viscosity)/(2.0*delta) << endl;
//    }
//
//    delete [] temp_val_Vars;
//    
//    cin.get();
//    // FD *****************

	delete [] Varsd;

}

void CAdjTurbSolution::ConvertSensGradPrimVar_to_SensPrimVar(double *val_Vars_i, double **val_Vars_ks,
		double **val_Jacobian_i, double ***val_Jacobian_gradi,
		unsigned long iPoint, short kNode, double **Normals,
		unsigned short numTotalVar, unsigned short nNeigh, unsigned short numEqn, CGeometry *geometry, CConfig *config){

	unsigned short iPos, jPos, kPos, lPos;
	unsigned short iNode;
	double *Vars_id, **Vars_ksd;
	double **GradPrimVar_i, **GradPrimVar_id;
    
    Vars_id = new double [numTotalVar];
    Vars_ksd = new double* [nNeigh];
    for (iPos = 0; iPos < nNeigh; iPos++){
        Vars_ksd[iPos] = new double[numTotalVar];
    }
    
    GradPrimVar_i = new double*[numTotalVar];
    GradPrimVar_id = new double*[numTotalVar];
    for (iPos = 0; iPos < numTotalVar; iPos++){
        GradPrimVar_i[iPos] = new double[nDim];
        GradPrimVar_id[iPos] = new double[nDim];
    }

	// zero val_Jacobian (reusing structure)
	for (iPos = 0; iPos < numTotalVar; iPos++) {
		for (jPos = 0; jPos < numEqn; jPos++) {

			val_Jacobian_i[iPos][jPos] = 0.0;

		}
	}
    //
//    for (jPos = 0; jPos < numTotalVar; jPos++){
//        cout << "@i: " << val_Vars_i[jPos] << endl;;
//        for (unsigned short jNode = 0; jNode < nNeigh; jNode++){
//            cout << "@ks: " << val_Vars_ks[jNode][jPos] << endl;
//        }
//    }
//    
//    for (jPos = 0; jPos < nDim; jPos++){
//        for (unsigned short jNode = 0; jNode < nNeigh; jNode++){
//            cout << "@Normals: " << Normals[jNode][jPos] << endl;
//        }
//    }
//
	for (iPos = 0; iPos < numTotalVar; iPos++){
		// zero things first
		for (jPos = 0; jPos < numTotalVar; jPos++){
			Vars_id[jPos] = 0.0;
			for (iNode = 0; iNode < nNeigh; iNode++){
				Vars_ksd[iNode][jPos] = 0.0;
			}
		}

		if (kNode == -1)
			Vars_id[iPos] = 1.0;
		else {
			Vars_ksd[kNode][iPos] = 1.0;
		}

		this->CalcGradient_GG_ad(val_Vars_i, Vars_id, val_Vars_ks, Vars_ksd,
				nNeigh, numTotalVar, Normals, GradPrimVar_i, GradPrimVar_id, config,
				geometry, iPoint);

		for (jPos = 0; jPos < nVar; jPos++) {
			for (kPos = 0; kPos < nDim; kPos++) {
                for (lPos = 0; lPos < numTotalVar; lPos++) {
                    val_Jacobian_i[iPos][jPos] += GradPrimVar_id[lPos][kPos]
                                                    * val_Jacobian_gradi[lPos][kPos][jPos];
                }
			}
		}
        
// TEST AD *****************
//        for (jPos = 0; jPos < numTotalVar; jPos++) {
//            for (kPos = 0; kPos < nDim; kPos++) {
//            //cout << "Direct: " << GradPrimVar_i[jPos] << endl;
//            cout << "AD: " << GradPrimVar_id[jPos][kPos] << endl;
//            }
//        }
//        cout << "--" << iPos << endl;
// AD *****************

	}
    
    // TEST FD *****************
//        double **temp1_val_GradPrimVar_id, **temp2_val_GradPrimVar_id;
//    double *temp_val_Vars_i, **temp_val_Vars_ks;
//    
//        double delta = 0.00001;
//        temp1_val_GradPrimVar_id = new double*[numTotalVar];
//        temp2_val_GradPrimVar_id = new double*[numTotalVar];
//        for (iPos = 0; iPos < numTotalVar; iPos++){
//            temp1_val_GradPrimVar_id[iPos] = new double[nDim];
//            temp2_val_GradPrimVar_id[iPos] = new double[nDim];
//        }
//    
//    temp_val_Vars_i = new double[numTotalVar];
//    temp_val_Vars_ks = new double*[nNeigh];
//    for (iPos = 0; iPos < nNeigh; iPos++){
//        temp_val_Vars_ks[iPos] = new double[numTotalVar];
//    }
//
//    for (iPos = 0; iPos < numTotalVar; iPos++){
//        
//        for (jPos = 0; jPos < numTotalVar; jPos++){
//            temp_val_Vars_i[jPos] = val_Vars_i[jPos];
//            for (unsigned short jNode = 0; jNode < nNeigh; jNode++){
//                temp_val_Vars_ks[jNode][jPos] = val_Vars_ks[jNode][jPos];
//            }
//        }
//        
//		if (kNode == -1)
//			temp_val_Vars_i[iPos] = val_Vars_i[iPos] - delta;
//		else {
//			temp_val_Vars_ks[kNode][iPos] = val_Vars_ks[kNode][iPos] - delta;
//		}
//    
//		this->CalcGradient_GG_ad(temp_val_Vars_i, Vars_id, temp_val_Vars_ks, Vars_ksd,
//                                 nNeigh, numTotalVar, Normals, temp1_val_GradPrimVar_id, GradPrimVar_id, config,
//                                 geometry, iPoint);
//        
//        if (kNode == -1)
//			temp_val_Vars_i[iPos] = val_Vars_i[iPos] + delta;
//		else {
//			temp_val_Vars_ks[kNode][iPos] = val_Vars_ks[kNode][iPos] + delta;
//		}
//        
//		this->CalcGradient_GG_ad(temp_val_Vars_i, Vars_id, temp_val_Vars_ks, Vars_ksd,
//                                 nNeigh, numTotalVar, Normals, temp2_val_GradPrimVar_id, GradPrimVar_id, config,
//                                 geometry, iPoint);
//        
//        for (jPos = 0; jPos < numTotalVar; jPos++) {
//            for (kPos = 0; kPos < nDim; kPos++) {
//                cout << "FD: " << (temp2_val_GradPrimVar_id[jPos][kPos] - temp1_val_GradPrimVar_id[jPos][kPos])/(2.0*delta) << endl;
//            }
//        }
//        cout << "--" << iPos << endl;
//    }
//
//
//    for (iPos = 0; iPos < numTotalVar; iPos++){
//        delete [] temp1_val_GradPrimVar_id[iPos];
//        delete [] temp2_val_GradPrimVar_id[iPos];
//    }
//    
//    delete [] temp1_val_GradPrimVar_id;
//    delete [] temp2_val_GradPrimVar_id;
//    
//    for (iPos = 0; iPos < nNeigh; iPos++){
//        delete [] temp_val_Vars_ks[iPos];
//    }
//    
//    delete [] temp_val_Vars_ks;
//    delete [] temp_val_Vars_i;
//    
//    cin.get();
    // FD *****************
    
    delete [] Vars_id;
    for (iPos = 0; iPos < nNeigh; iPos++){
        delete [] Vars_ksd[iPos];
    }
    
    delete [] Vars_ksd;
    
    for (iPos = 0; iPos < numTotalVar; iPos++){
        delete [] GradPrimVar_i[iPos];
        delete [] GradPrimVar_id[iPos];
    }
    delete [] GradPrimVar_i;
    delete [] GradPrimVar_id;
    
}

void CAdjTurbSolution::ConvertSensPrimVar_to_SensU(double *val_Vars, double **val_Jacobian,
		unsigned short numVar, CConfig *config){

	unsigned short iPos, jPos, lPos;
	double Gas_Constant = config->GetGas_Constant();
	double *Primitive, *Primitived, *val_Varsd;
    double **temp_Jacobian;
    
    val_Varsd = new double[numVar];
    Primitive = new double[numVar];
    Primitived = new double[numVar];
    
    temp_Jacobian = new double*[numVar];
    for (iPos = 0; iPos < numVar; iPos++)
        temp_Jacobian[iPos] = new double[nVar];
    
    for (iPos = 0; iPos < numVar; iPos++)
		for (jPos = 0; jPos < nVar; jPos++)
            temp_Jacobian[iPos][jPos] = 0.0;
    
	for (iPos = 0; iPos < numVar; iPos++){
		// zero things first
		for (jPos = 0; jPos < numVar; jPos++){
			val_Varsd[jPos] = 0.0;
		}

		val_Varsd[iPos] = 1.0;

		this->CalcPrimVar_Compressible_ad(val_Vars, val_Varsd, Gamma, Gas_Constant, numVar,
				0.0, Primitive, Primitived, config);

		for (jPos = 0; jPos < nVar; jPos++) {
            for (lPos = 0; lPos < numVar; lPos++) {
                temp_Jacobian[iPos][jPos] += val_Jacobian[lPos][jPos] * Primitived[lPos];
            }
		}
        
////        // TEST AD *****************
////        for (jPos = 0; jPos < numVar; jPos++) {
////            //cout << "Direct: " << Primitive[jPos] << endl;
////        cout << "AD: " << Primitived[jPos] << endl;
////            }
////        cout << "--" << iPos << endl;
////// AD *****************
	}
    
    for (iPos = 0; iPos < numVar; iPos++)
		for (jPos = 0; jPos < nVar; jPos++)
            val_Jacobian[iPos][jPos] = temp_Jacobian[iPos][jPos];
    
////    // TEST FD *****************
////        double *temp1_val_Primitive, *temp2_val_Primitive, *temp_val_Vars;
////        double delta = 0.000001;
////    temp1_val_Primitive = new double[numVar];
////    temp2_val_Primitive = new double[numVar];
////        temp_val_Vars = new double[numVar];
////    
////    Solution = temp_val_Vars;
////    
////        for (iPos = 0; iPos < numVar; iPos++){
////            
////            for (jPos = 0; jPos < numVar; jPos++) {
////                temp_val_Vars[jPos] = val_Vars[jPos];
////            }
////                            
////            //temp_val_Vars[iPos] = (1.0-delta) * val_Vars[iPos];
////            temp_val_Vars[iPos] = val_Vars[iPos] - delta;
////
////            this->CalcPrimVar_Compressible_ad(Gamma, Gas_Constant, numVar,
////                                              0.0, temp1_val_Primitive, Primitived, config);
////                
////            //temp_val_Vars[iPos] = (1.0-delta) * val_Vars[iPos];
////            temp_val_Vars[iPos] = val_Vars[iPos] + delta;
////
////            this->CalcPrimVar_Compressible_ad(Gamma, Gas_Constant, numVar,
////                                              0.0, temp2_val_Primitive, Primitived, config);
////    
////            //cout << "FD: " << (temp2_val_Laminar_Viscosity - temp1_val_Laminar_Viscosity)/(2.0*delta*val_U_i[iPos]) << endl;
////            for (jPos = 0; jPos < numVar; jPos++) {
////                //cout << "Direct1: " << temp1_val_Primitive[jPos] << endl;
////                //cout << "Direct2: " << temp2_val_Primitive[jPos] << endl;
////                  cout << "FD: " << (temp2_val_Primitive[jPos] - temp1_val_Primitive[jPos])/(2.0*delta) << endl;
////            }
////            cout << "--" << iPos << endl;
////        }
////    
////    delete [] temp1_val_Primitive;
////    delete [] temp2_val_Primitive;
////        delete [] temp_val_Vars;
////    cin.get();
////    // FD *****************
//    
    delete [] val_Varsd;
    delete [] Primitive;
    delete [] Primitived;
    
    for (iPos = 0; iPos < numVar; iPos++)
        delete [] temp_Jacobian[iPos];

    delete [] temp_Jacobian;
}

void CAdjTurbSolution::CalcEddyViscSens(double *val_U_i, double val_Laminar_Viscosity_i,
		double *val_TurbVar_i, double *EddyViscSens, unsigned short nFlowVar, CConfig *config){

	unsigned short iPos, jPos;
	double *U_id, *TurbVar_id, Laminar_Viscosity_id, *Eddy_Viscosity_i, *Eddy_Viscosity_id;
	double **tempEddyViscSens;

	U_id = new double[nFlowVar];
	TurbVar_id = new double[nVar];
	tempEddyViscSens = new double*[nVar+nFlowVar];
	for (iPos = 0; iPos < (nVar+nFlowVar); iPos++)
		tempEddyViscSens[iPos] = new double[1];
	Eddy_Viscosity_i = new double;
	Eddy_Viscosity_id = new double;
    
//// TEST *****************
//        for (iPos = 0; iPos < nFlowVar; iPos++){
//            cout << "Flow: " << val_U_i[iPos] << endl;
//        }
//    for (iPos = 0; iPos < nVar; iPos++){
//        val_TurbVar_i[iPos] = 0.01;
//        cout << "Turb: " << val_TurbVar_i[iPos] << endl;
//    }
//// *****************

	// from U_i
	for (iPos = 0; iPos < nFlowVar; iPos++){
		// zero things first
		for (jPos = 0; jPos < nFlowVar; jPos++){
			U_id[jPos] = 0.0;
		}

		for (jPos = 0; jPos < nVar; jPos++){
			TurbVar_id[jPos] = 0.0;
		}

		Laminar_Viscosity_id = 0.0;


		(*Eddy_Viscosity_id) = 0.0;

		U_id[iPos] = 1.0;

		if (config->GetKind_Turb_Model() == SA)
			this->CalcEddyViscosity_ad(val_U_i, U_id, val_Laminar_Viscosity_i, Laminar_Viscosity_id,
					val_TurbVar_i, TurbVar_id, Eddy_Viscosity_i, Eddy_Viscosity_id, config);

		EddyViscSens[iPos] = (*Eddy_Viscosity_id);
        
//// TEST AD *****************
//        cout << "AD: " << (*Eddy_Viscosity_id) << endl;
//// AD *****************

	}
    
//    // TEST FD *****************
//    double *temp_val_Vars, temp1_val_Eddy_Viscosity, temp2_val_Eddy_Viscosity;
//    double delta = 0.00000001;
//    temp_val_Vars = new double[nFlowVar];
//
//    for (iPos = 0; iPos < nFlowVar; iPos++){
//
//        for (jPos = 0; jPos < nFlowVar; jPos++){
//			temp_val_Vars[jPos] = val_U_i[jPos];
//		}
//
//        //temp_val_Vars[iPos] = (1.0-delta) * val_Vars[iPos];
//        temp_val_Vars[iPos] = val_U_i[iPos] - delta;
//
//        this->CalcEddyViscosity_ad(temp_val_Vars, U_id, val_Laminar_Viscosity_i, Laminar_Viscosity_id,
//                                   val_TurbVar_i, TurbVar_id, &temp1_val_Eddy_Viscosity, Eddy_Viscosity_id, config);
//
//        //temp_val_Vars[iPos] = (1.0+delta) * val_Vars[iPos];
//        temp_val_Vars[iPos] = val_U_i[iPos] + delta;
//
//
//        this->CalcEddyViscosity_ad(temp_val_Vars, U_id, val_Laminar_Viscosity_i, Laminar_Viscosity_id,
//                                   val_TurbVar_i, TurbVar_id, &temp2_val_Eddy_Viscosity, Eddy_Viscosity_id, config);
//
//        //cout << "FD: " << (temp2_val_Laminar_Viscosity - temp1_val_Laminar_Viscosity)/(2.0*delta*val_U_i[iPos]) << endl;
//              cout << "FD: " << (temp2_val_Eddy_Viscosity - temp1_val_Eddy_Viscosity)/(2.0*delta) << endl;
//    }
//
//    delete [] temp_val_Vars;
//    
//    //cin.get();
//    // FD *****************

	// from TurbVar_i

	for (iPos = 0; iPos < nVar; iPos++){
		// zero things first
		for (jPos = 0; jPos < nFlowVar; jPos++){
			U_id[jPos] = 0.0;
		}

		for (jPos = 0; jPos < nVar; jPos++){
			TurbVar_id[jPos] = 0.0;
		}

		Laminar_Viscosity_id = 0.0;

		(*Eddy_Viscosity_id) = 0.0;

		TurbVar_id[iPos] = 1.0;

		if (config->GetKind_Turb_Model() == SA)
			this->CalcEddyViscosity_ad(val_U_i, U_id, val_Laminar_Viscosity_i, Laminar_Viscosity_id,
					val_TurbVar_i, TurbVar_id, Eddy_Viscosity_i, Eddy_Viscosity_id, config);

		EddyViscSens[iPos+nFlowVar] = (*Eddy_Viscosity_id);

//// TEST AD *****************
//        cout << "Direct: " << (*Eddy_Viscosity_i) << endl;
//        cout << "AD: " << (*Eddy_Viscosity_id) << endl;
//// AD *****************
        
	}
    
//    //if (fabs((*Eddy_Viscosity_id)) > 0.0000000000000000001)
//        cin.get(); 
    
//    // TEST FD *****************
//    temp_val_Vars = new double[nVar];
//    
//    for (iPos = 0; iPos < nVar; iPos++){
//        
//        for (jPos = 0; jPos < nVar; jPos++){
//			temp_val_Vars[jPos] = val_TurbVar_i[jPos];
//		}
//        
//        //temp_val_Vars[iPos] = (1.0-delta) * val_Vars[iPos];
//        temp_val_Vars[iPos] = val_TurbVar_i[iPos] - delta;
//        
//        this->CalcEddyViscosity_ad(val_U_i, U_id, val_Laminar_Viscosity_i, Laminar_Viscosity_id,
//                                   temp_val_Vars, TurbVar_id, &temp1_val_Eddy_Viscosity, Eddy_Viscosity_id, config);
//        
//        //temp_val_Vars[iPos] = (1.0+delta) * val_Vars[iPos];
//        temp_val_Vars[iPos] = val_TurbVar_i[iPos] + delta;
//        
//        
//        this->CalcEddyViscosity_ad(val_U_i, U_id, val_Laminar_Viscosity_i, Laminar_Viscosity_id,
//                                   temp_val_Vars, TurbVar_id, &temp2_val_Eddy_Viscosity, Eddy_Viscosity_id, config);
//        
//        //cout << "FD: " << (temp2_val_Laminar_Viscosity - temp1_val_Laminar_Viscosity)/(2.0*delta*val_TurbVar_i[iPos]) << endl;
//        cout << "FD: " << (temp2_val_Eddy_Viscosity - temp1_val_Eddy_Viscosity)/(2.0*delta) << endl;
//    }
//    
//    delete [] temp_val_Vars;
//    
//    //if (fabs((*Eddy_Viscosity_id)) > 0.0000000000000000001)
//     //   cin.get();
//    // FD *****************

	// from mui
	// zero things first
	for (jPos = 0; jPos < nFlowVar; jPos++){
		U_id[jPos] = 0.0;
	}

	for (jPos = 0; jPos < nVar; jPos++){
		TurbVar_id[jPos] = 0.0;
	}

	Laminar_Viscosity_id = 1.0;

	(*Eddy_Viscosity_id) = 0.0;

	if (config->GetKind_Turb_Model() == SA)
		this->CalcEddyViscosity_ad(val_U_i, U_id, val_Laminar_Viscosity_i, Laminar_Viscosity_id,
				val_TurbVar_i, TurbVar_id, Eddy_Viscosity_i, Eddy_Viscosity_id, config);
    
//    // TEST AD *****************
//        cout << "AD: " << (*Eddy_Viscosity_id) << endl;
//    
//    // AD *****************
    
//    // TEST FD *****************
//    double temp_mu;
//    cout << val_Laminar_Viscosity_i << endl;
//    
//    for (iPos = 0; iPos < nVar; iPos++){
//        
//        temp_mu = val_Laminar_Viscosity_i - delta;
//        
//        this->CalcEddyViscosity_ad(val_U_i, U_id, temp_mu, Laminar_Viscosity_id,
//                                   val_TurbVar_i, TurbVar_id, &temp1_val_Eddy_Viscosity, Eddy_Viscosity_id, config);
//        
//        temp_mu = val_Laminar_Viscosity_i + delta;
//        
//        
//        this->CalcEddyViscosity_ad(val_U_i, U_id, temp_mu, Laminar_Viscosity_id,
//                                   val_TurbVar_i, TurbVar_id, &temp2_val_Eddy_Viscosity, Eddy_Viscosity_id, config);
//        
//        //cout << "FD: " << (temp2_val_Laminar_Viscosity - temp1_val_Laminar_Viscosity)/(2.0*delta*val_TurbVar_i[iPos]) << endl;
//        cout << "FD: " << (temp2_val_Eddy_Viscosity - temp1_val_Eddy_Viscosity)/(2.0*delta) << endl;
//    }
//        
//       cin.get();
//    // FD *****************

    

	this->ConvertSensMu_to_SensU(val_U_i, val_Laminar_Viscosity_i,
			tempEddyViscSens, Eddy_Viscosity_id,
			nFlowVar, nVar+nFlowVar, nVar, config);

	for (jPos = 0; jPos < nFlowVar; jPos++){
//        cout << EddyViscSens[jPos] << endl;
//        cout << tempEddyViscSens[jPos][0] << endl;
		EddyViscSens[jPos] += tempEddyViscSens[jPos][0];
        
	}
   // cin.get();

	delete [] U_id;
	delete [] TurbVar_id;
	for (iPos = 0; iPos < (nVar+nFlowVar); iPos++)
		delete [] tempEddyViscSens[iPos];
	delete [] tempEddyViscSens;
	delete [] Eddy_Viscosity_i;
	delete [] Eddy_Viscosity_id;

}
