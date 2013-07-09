/*!
 * \file solution_direct_turbulent.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.5
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

CTurbSolution::CTurbSolution(void) : CSolution() {

  /*--- Array initialization ---*/
	FlowSolution_i = NULL;
	FlowSolution_j = NULL;
	lowerlimit = NULL;
	upperlimit = NULL;

}

CTurbSolution::CTurbSolution(CConfig *config) : CSolution() {

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Array initialization ---*/
	FlowSolution_i = NULL;
	FlowSolution_j = NULL;
	lowerlimit = NULL;
	upperlimit = NULL;
  
}

CTurbSolution::~CTurbSolution(void) {

	if (FlowSolution_i != NULL) delete [] FlowSolution_i;
	if (FlowSolution_j != NULL) delete [] FlowSolution_j;
	if (lowerlimit != NULL) delete [] lowerlimit;
	if (upperlimit != NULL) delete [] upperlimit;
  
}

void CTurbSolution::SetSolution_MPI(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker;
	unsigned long iVertex, iPoint, nVertex, nBuffer_Vector, nBuffer_Scalar;
	double *Buffer_Receive_U = NULL, *Buffer_Receive_muT = NULL;
	short SendRecv;
	int send_to, receive_from;
  
#ifndef NO_MPI
	double *Buffer_Send_U = NULL, *Buffer_Send_muT = NULL;
#endif
  
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
      
			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			nVertex = geometry->nVertex[iMarker];
			nBuffer_Vector = nVertex*nVar;
      nBuffer_Scalar = nVertex;

			send_to = SendRecv-1;
			receive_from = abs(SendRecv)-1;
      
#ifndef NO_MPI
      
			/*--- Send information using MPI  ---*/
			if (SendRecv > 0) {
        Buffer_Send_muT = new double[nBuffer_Scalar];
        Buffer_Send_U = new double[nBuffer_Vector];
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Buffer_Send_muT[iVertex] = node[iPoint]->GetmuT();
          for (iVar = 0; iVar < nVar; iVar++)
            Buffer_Send_U[iVar*nVertex+iVertex] = node[iPoint]->GetSolution(iVar);
				}
        MPI::COMM_WORLD.Bsend(Buffer_Send_U, nBuffer_Vector, MPI::DOUBLE, send_to, 0); delete [] Buffer_Send_U;
        MPI::COMM_WORLD.Bsend(Buffer_Send_muT, nBuffer_Scalar, MPI::DOUBLE, send_to, 1); delete [] Buffer_Send_muT;
			}
      
#endif
      
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
        Buffer_Receive_U = new double [nBuffer_Vector];
        Buffer_Receive_muT = new double [nBuffer_Scalar];
        
#ifdef NO_MPI
				/*--- Get the information from the donor point directly. This is a
				 serial computation with access to all nodes. Note that there is an
				 implicit ordering in the list. ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Buffer_Receive_muT[iVertex] = node[iPoint]->GetmuT();
          for (iVar = 0; iVar < nVar; iVar++)
            Buffer_Receive_U[iVar*nVertex+iVertex] = node[iPoint]->GetSolution(iVar);
        }
#else
        MPI::COMM_WORLD.Recv(Buffer_Receive_U, nBuffer_Vector, MPI::DOUBLE, receive_from, 0);
        MPI::COMM_WORLD.Recv(Buffer_Receive_muT, nBuffer_Scalar, MPI::DOUBLE, receive_from, 1);
#endif
        
				/*--- Do the coordinate transformation ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          node[iPoint]->SetmuT(Buffer_Receive_muT[iVertex]);
          for (iVar = 0; iVar < nVar; iVar++)
            node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertex+iVertex]);
				}
        delete [] Buffer_Receive_U;
        delete [] Buffer_Receive_muT;
			}
		}
	}
}

void CTurbSolution::SetSolution_Limiter_MPI(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker;
	unsigned long iVertex, iPoint, nVertex, nBuffer_Vector;
	double *Buffer_Receive_Limit = NULL;
	short SendRecv;
	int send_to, receive_from;
  
#ifndef NO_MPI
	double *Buffer_Send_Limit;
#endif
    
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
      
			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			nVertex = geometry->nVertex[iMarker];
			nBuffer_Vector			= nVertex*nVar;
      
			send_to = SendRecv-1;
			receive_from = abs(SendRecv)-1;
      
#ifndef NO_MPI
      
			/*--- Send information using MPI  ---*/
			if (SendRecv > 0) {
        
				/*--- Allocate upwind variables ---*/
				Buffer_Send_Limit = new double[nBuffer_Vector];
        
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
          
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
					/*--- Copy upwind data ---*/
          for (iVar = 0; iVar < nVar; iVar++)
            Buffer_Send_Limit[iVar*nVertex+iVertex] = node[iPoint]->GetLimiter(iVar);
          
				}
        
				/*--- Send the buffer, and deallocate information using MPI ---*/
        MPI::COMM_WORLD.Bsend(Buffer_Send_Limit, nBuffer_Vector, MPI::DOUBLE, send_to, 2); delete [] Buffer_Send_Limit;
        
			}
      
#endif
      
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
        
				/*--- Allocate upwind variables ---*/
				Buffer_Receive_Limit = new double [nBuffer_Vector];
        
#ifdef NO_MPI
        
				/*--- Get the information from the donor point directly. This is a
				 serial computation with access to all nodes. Note that there is an
				 implicit ordering in the list. ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          /*--- Copy upwind data ---*/
          for (iVar = 0; iVar < nVar; iVar++)
            Buffer_Receive_Limit[iVar*nVertex+iVertex] = node[iPoint]->GetLimiter(iVar);
          
				}
        
#else
				/*--- Receive the information using MPI---*/
        MPI::COMM_WORLD.Recv(Buffer_Receive_Limit, nBuffer_Vector, MPI::DOUBLE, receive_from, 2);
#endif
        
				/*--- Do the coordinate transformation ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
          
					/*--- Find point and its type of transformation ---*/
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
					/*--- Upwind method. Store the received information ---*/
          for (iVar = 0; iVar < nVar; iVar++)
            node[iPoint]->SetLimiter(iVar, Buffer_Receive_Limit[iVar*nVertex+iVertex]);
          
				}
        
				delete [] Buffer_Receive_Limit;
        
			}
		}
	}
    
}

void CTurbSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	/*--- Convective fluxes across symmetry plane are equal to zero. ---*/
}

void CTurbSolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container,
                                  CNumerics *solver, CConfig *config, unsigned short val_marker) {
	/*--- Convective fluxes across euler wall are equal to zero. ---*/
}

void CTurbSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta, Vol, density_old, density;
    
    bool adjoint = ((config->GetKind_Solver() == ADJ_EULER) || (config->GetKind_Solver() == ADJ_NAVIER_STOKES) ||
                    (config->GetKind_Solver() == ADJ_RANS) || (config->GetKind_Solver() == ADJ_FREE_SURFACE_EULER) ||
                    (config->GetKind_Solver() == ADJ_FREE_SURFACE_NAVIER_STOKES) || (config->GetKind_Solver() == ADJ_FREE_SURFACE_RANS) ||
                    (config->GetKind_Solver() == ADJ_PLASMA_EULER) || (config->GetKind_Solver() == ADJ_PLASMA_NAVIER_STOKES) ||
                    (config->GetKind_Solver() == ADJ_AEROACOUSTIC_EULER));
    
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
    }
    
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
        
        /*--- Read the volume ---*/
		Vol = geometry->node[iPoint]->GetVolume();
        
		/*--- Modify matrix diagonal to assure diagonal dominance ---*/
		Delta = Vol / (config->GetTurb_CFLRedCoeff()*solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
		Jacobian.AddVal2Diag(iPoint,Delta);
        
        /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			xres[total_index] = - xres[total_index];
			xsol[total_index] = 0.0;
			AddRes_RMS(iVar, xres[total_index]*xres[total_index]);
            AddRes_Max(iVar, fabs(xres[total_index]), geometry->node[iPoint]->GetGlobalIndex());
		}
	}
    
    /*--- Initialize residual and solution at the ghost points ---*/
    for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
            total_index = iPoint*nVar + iVar;
            xres[total_index] = 0.0;
            xsol[total_index] = 0.0;
        }
    }
    
	/*--- Solve the linear system (Stationary iterative methods) ---*/
	if (config->GetKind_Linear_Solver() == SYM_GAUSS_SEIDEL)
		Jacobian.SGSSolution(xres, xsol, config->GetLinear_Solver_Error(),
                             config->GetLinear_Solver_Iter(), false, geometry, config);
    
	if (config->GetKind_Linear_Solver() == LU_SGS)
		Jacobian.LU_SGSIteration(xres, xsol, geometry, config);
    
	/*--- Solve the linear system (Krylov subspace methods) ---*/
	if ((config->GetKind_Linear_Solver() == BCGSTAB) ||
        (config->GetKind_Linear_Solver() == GMRES)) {
        
		CSysVector rhs_vec(nPoint, nPointDomain, nVar, xres);
		CSysVector sol_vec(nPoint, nPointDomain, nVar, xsol);
        
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
			system.GMRES(rhs_vec, sol_vec, *mat_vec, *precond, *sol_mpi, config->GetLinear_Solver_Error(),
                         config->GetLinear_Solver_Iter(), false);
        
		sol_vec.CopyToArray(xsol);
		delete mat_vec;
		delete precond;
		delete sol_mpi;
	}
    
	/*--- Update solution (system written in terms of increments) ---*/
	switch (config->GetKind_Turb_Model()){
        case SA:
            if (!adjoint) {
                for (iPoint = 0; iPoint < nPointDomain; iPoint++)
                    for (iVar = 0; iVar < nVar; iVar++)
                        node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*xsol[iPoint*nVar+iVar]);
            }
            break;
        case SST:
            if (!adjoint) {
                for (iPoint = 0; iPoint < nPointDomain; iPoint++){
                    density_old = solution_container[FLOW_SOL]->node[iPoint]->GetSolution_Old(0);
                    density     = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(0);
                    
                    for (iVar = 0; iVar < nVar; iVar++)
                        node[iPoint]->AddConservativeSolution(iVar, config->GetLinear_Solver_Relax()*xsol[iPoint*nVar+iVar],
                                                              density, density_old, lowerlimit[iVar], upperlimit[iVar]);
                }
            }
            break;
	}
    
    /*--- MPI solution ---*/
    SetSolution_MPI(geometry, config);
    
    /*--- Compute the root mean square residual ---*/
    SetResidual_RMS(geometry, config);
    
}

void CTurbSolution::CalcGradient_GG(double *val_U_i, double **val_U_js, unsigned short nNeigh, unsigned short numVar,
		double **Normals, double **grad_U_i, CConfig *config, CGeometry *geometry, unsigned long iPoint) {

	double Volume;
	Volume = geometry->node[iPoint]->GetVolume();

	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//

	//SU2_CPP2C START CTurbSolution::CalcGradient_GG
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *val_U_i **val_U_js
	//SU2_CPP2C OUTVARS **grad_U_i
	//SU2_CPP2C VARS DOUBLE **Normals Volume
	//SU2_CPP2C VARS INT nNeigh numVar
	//SU2_CPP2C CALL_LIST END

	//SU2_CPP2C DEFINE nDim nVar

	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C DECL_LIST END

	unsigned short iVar, iNeigh, jVar, kDim;

	/*--- Set Gradient to Zero ---*/
	for(iVar = 0; iVar < numVar; iVar++)
		for(kDim = 0; kDim < nDim; kDim++)
			grad_U_i[iVar][kDim] = 0.0;

	for(iNeigh = 0; iNeigh < nNeigh; iNeigh++) {

		for(jVar = 0; jVar < numVar; jVar++)
			for(kDim = 0; kDim < nDim; kDim++)
				grad_U_i[jVar][kDim] += 0.5*(val_U_i[jVar] + val_U_js[iNeigh][jVar])
				*Normals[iNeigh][kDim]/Volume;
	}

	//SU2_CPP2C END CTurbSolution::CalcGradient_GG
}

void CTurbSolution::CalcGradient_LS(double *val_U_i, double **val_U_js, unsigned long nNeigh,
		double *coords_i, double **coords_js, double **grad_U_i, CConfig *config) {


	unsigned short iDim, jDim, iNeigh, iVar;
	double r11, r12, r13, r22, r23, r33, r23_a, r23_b;
	double weight, detR2, z11, z12, z13, z22, z23, z33;
	double **Smatrix, **cvector;

	Smatrix = new double*[nDim];
	for (iVar=0; iVar<nDim; iVar++)
		Smatrix[iVar] = new double[nDim];

	cvector = new double*[nVar];
	for (iVar=0; iVar<nVar; iVar++)
		cvector[iVar] = new double[nDim];

	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//

	//SU2_CPP2C START CTurbSolution::CalcGradient_LS
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *val_U_i **val_U_js
	//SU2_CPP2C OUTVARS **grad_U_i
	//SU2_CPP2C VARS DOUBLE *coords_i **coords_js
	//SU2_CPP2C VARS INT nNeigh
	//SU2_CPP2C CALL_LIST END

	//SU2_CPP2C DEFINE nDim nVar

	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS INT SCALAR iDim jDim iNeigh iVar
	//SU2_CPP2C VARS DOUBLE SCALAR r11 r12 r13 r22 r23 r33 r23_a r23_b
	//SU2_CPP2C VARS DOUBLE SCALAR weight detR2 z11 z12 z13 z22 z23 z33
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nDim SIZE=nDim Smatrix
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nVar SIZE=nDim cvector
	//SU2_CPP2C DECL_LIST END

	/*--- Inizialization of variables ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		for (iDim = 0; iDim < nDim; iDim++)
			cvector[iVar][iDim] = 0.0;
	r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0; r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

	for (iNeigh = 0; iNeigh < nNeigh; iNeigh++) {

		weight = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			weight += (coords_js[iDim][iNeigh]-coords_i[iDim])*(coords_js[iDim][iNeigh]-coords_i[iDim]);

		/*--- Sumations for entries of upper triangular matrix R ---*/
		r11 += (coords_js[0][iNeigh]-coords_i[0])*(coords_js[0][iNeigh]-coords_i[0])/(weight);
		r12 += (coords_js[0][iNeigh]-coords_i[0])*(coords_js[1][iNeigh]-coords_i[1])/(weight);
		r22 += (coords_js[1][iNeigh]-coords_i[1])*(coords_js[1][iNeigh]-coords_i[1])/(weight);
		if (nDim == 3) {
			r13 += (coords_js[0][iNeigh]-coords_i[0])*(coords_js[2][iNeigh]-coords_i[2])/(weight);
			r23_a += (coords_js[1][iNeigh]-coords_i[1])*(coords_js[2][iNeigh]-coords_i[2])/(weight);
			r23_b += (coords_js[0][iNeigh]-coords_i[0])*(coords_js[2][iNeigh]-coords_i[2])/(weight);
			r33 += (coords_js[2][iNeigh]-coords_i[2])*(coords_js[2][iNeigh]-coords_i[2])/(weight);
		}

		/*--- Entries of c:= transpose(A)*b ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			for (iDim = 0; iDim < nDim; iDim++)
				cvector[iVar][iDim] += (coords_js[iDim][iNeigh]-coords_i[iDim])*(val_U_js[iVar][iNeigh]-val_U_i[iVar])/(weight);
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
		detR2 = (r11*r22)*(r11*r22);
		Smatrix[0][0] = (r12*r12+r22*r22)/(detR2);
		Smatrix[0][1] = -r11*r12/(detR2);
		Smatrix[1][0] = Smatrix[0][1];
		Smatrix[1][1] = r11*r11/(detR2);
	}
	else {
		detR2 = (r11*r22*r33)*(r11*r22*r33);
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
	for (iVar = 0; iVar < nVar; iVar++) {
		for (iDim = 0; iDim < nDim; iDim++) {
			grad_U_i[iVar][iDim] = 0.0;
			for (jDim = 0; jDim < nDim; jDim++)
				grad_U_i[iVar][iDim] += Smatrix[iDim][jDim]*cvector[iVar][jDim];
		}
	}

	//SU2_CPP2C END CTurbSolution::CalcGradient_LS

	for (iVar=0; iVar<nDim; iVar++)
		delete [] Smatrix[iVar];
	delete [] Smatrix;

	for (iVar=0; iVar<nDim; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;

}

void CTurbSolution::CalcPrimVar_Compressible(double *val_Vars, double Gamma, double Gas_Constant,
		unsigned short numVar, double turb_ke, double* Primitive, CConfig *config) {

	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//

	//SU2_CPP2C START CTurbSolution::CalcPrimVar_Compressible
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *val_Vars
	//SU2_CPP2C OUTVARS *Primitive
	//SU2_CPP2C VARS DOUBLE Gamma Gas_Constant turb_ke
	//SU2_CPP2C VARS INT nNeigh
	//SU2_CPP2C CALL_LIST END

	//SU2_CPP2C DEFINE nDim numVar

	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C DECL_LIST END

	unsigned short iDim;

	double Velocity2;

	//SetVelocity2();								 // Compute the modulus of the velocity.
	Velocity2 = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		Velocity2 += val_Vars[iDim+1]*val_Vars[iDim+1]/(val_Vars[0]*val_Vars[0]);
	//SetPressure(Gamma, turb_ke);   // Requires Velocity2 computation.
	Primitive[nDim+1] = (Gamma-1.0)*val_Vars[0]*(val_Vars[numVar-1]/val_Vars[0]-0.5*Velocity2 - turb_ke); //note: if turb_ke used, need to include its sensitivity...
    //SU2_CPP2C COMMENT START
	//SetSoundSpeed(Gamma);					 // Requires pressure computation.
	Primitive[nDim+4] = sqrt(Gamma*Primitive[nDim+1]/val_Vars[0]);
	//SetEnthalpy();								 // Requires pressure computation.
	Primitive[nDim+3] = (val_Vars[numVar-1] + Primitive[nDim+1]) / val_Vars[0];
    //SU2_CPP2C COMMENT END
	//SetTemperature(Gas_Constant);  // Requires pressure computation.
	Primitive[0] = Primitive[nDim+1] / ( Gas_Constant * val_Vars[0]);
	//SetLaminarViscosity();				 // Requires temperature computation.
	//CalcLaminarViscosity( Solution, LaminarViscosity, config);

	for (iDim = 0; iDim < nDim; iDim++)
		Primitive[iDim+1] = val_Vars[iDim+1] / val_Vars[0];
    //SU2_CPP2C COMMENT START
	Primitive[nDim+2] = val_Vars[0];
//SU2_CPP2C COMMENT END
	//SU2_CPP2C END CTurbSolution::CalcPrimVar_Compressible

}

//void CTurbSolution::CalcPrimVar_Incompressible(double Density_Inf, double Viscosity_Inf, double ArtComp_Factor, double turb_ke, bool freesurface,
//		double* LaminarViscosityInc, double* Primitive) {
//	unsigned short iDim;
//
//	double Velocity2;
//
//	//SetBetaInc2(ArtComp_Factor);	// Set the value of beta.
//	Primitive[nDim+1] = ArtComp_Factor;
//	if (!freesurface) {
//		//SetDensityInc(Density_Inf);							// Set the value of the density
//		Primitive[0] = Density_Inf;
//		//SetLaminarViscosityInc(Viscosity_Inf);	// Set the value of the viscosity
//		LaminarViscosityInc = Viscosity_Inf;
//	}
//	//SetVelocityInc2();						// Requires density computation.
//	Velocity2 = 0.0;
//	for (unsigned short iDim = 0; iDim < nDim; iDim++)
//		Velocity2 += (Solution[iDim+1]/Primitive[0])*(Solution[iDim+1]/Primitive[0]);
//
//	for (iDim = 0; iDim < nDim; iDim++)
//		Primitive[iDim+1] = Solution[iDim+1] / Primitive[0];
//
//}

void CTurbSolution::CalcLaminarViscosity(double *val_U_i, double *val_laminar_viscosity_i, CConfig *config) {

	double Temperature_Ref, Viscosity_Ref, Gas_Constant;

	Temperature_Ref = config->GetTemperature_Ref();
	Viscosity_Ref = config->GetViscosity_Ref();
	Gas_Constant = config->GetGas_ConstantND();

	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//

	//SU2_CPP2C START CTurbSolution::CalcLaminarViscosity
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *val_U_i
	//SU2_CPP2C OUTVARS *val_laminar_viscosity_i
	//SU2_CPP2C VARS DOUBLE Temperature_Ref Viscosity_Ref Gamma_Minus_One Gas_Constant
	//SU2_CPP2C CALL_LIST END

	//SU2_CPP2C DEFINE nDim

	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C DECL_LIST END

	double Temperature, Temperature_Dim;
	double Density, Pressure;

	Density = val_U_i[0]; // Not incompressible
	if (nDim == 2)
		Pressure = Gamma_Minus_One*(val_U_i[3] - (val_U_i[1]*val_U_i[1] + val_U_i[2]*val_U_i[2])/(2.0*val_U_i[0]));
	else
		Pressure = Gamma_Minus_One*(val_U_i[3] - (val_U_i[1]*val_U_i[1] + val_U_i[2]*val_U_i[2] + val_U_i[3]*val_U_i[3])/(2.0*val_U_i[0]));

	Temperature = Pressure/(Gas_Constant*Density);

	/*--- Calculate viscosity from a non-dim. Sutherland's Law ---*/
	Temperature_Dim = Temperature*Temperature_Ref;
	(*val_laminar_viscosity_i) = 1.853E-5*(pow(Temperature_Dim/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_Dim+110.3));
	(*val_laminar_viscosity_i) = (*val_laminar_viscosity_i)/Viscosity_Ref;

	//SU2_CPP2C END CTurbSolution::CalcLaminarViscosity
}

CTurbSASolution::CTurbSASolution(void) : CTurbSolution() { }

CTurbSASolution::CTurbSASolution(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTurbSolution() {
	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double Density_Inf, Viscosity_Inf, Factor_nu_Inf, Density;

	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	bool incompressible = config->GetIncompressible();

	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Dimension of the problem --> dependent of the turbulent model ---*/
	nVar = 1;
	nPoint = geometry->GetnPoint();
	nPointDomain = geometry->GetnPointDomain();
  
  /*--- Define geometry constants in the solver structure ---*/
	nDim = geometry->GetnDim();
	node = new CVariable*[nPoint];
  
	/*--- Single grid simulation ---*/
	if (iMesh == MESH_0) {

		/*--- Define some auxiliar vector related with the residual ---*/
		Residual = new double[nVar]; Residual_RMS = new double[nVar];
		Residual_i = new double[nVar]; Residual_j = new double[nVar];
    Residual_Max = new double[nVar]; Point_Max = new unsigned long[nVar];
    

		/*--- Define some auxiliar vector related with the solution ---*/
		Solution = new double[nVar];
		Solution_i = new double[nVar]; Solution_j = new double[nVar];

		/*--- Define some auxiliar vector related with the geometry ---*/
		Vector_i = new double[nDim]; Vector_j = new double[nDim];

		/*--- Define some auxiliar vector related with the flow solution ---*/
		FlowSolution_i = new double [nDim+2]; FlowSolution_j = new double [nDim+2];

		/*--- Jacobians and vector structures for implicit computations ---*/
    Jacobian_i = new double* [nVar];
    Jacobian_j = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new double [nVar];
      Jacobian_j[iVar] = new double [nVar];
    }
    
    /*--- Initialization of the structure of the whole Jacobian ---*/
    if (rank == MASTER_NODE) cout << "Initialize jacobian structure (SA model)." << endl;
    Initialize_SparseMatrix_Structure(&Jacobian, nVar, nVar, geometry, config);
    xsol = new double [nPoint*nVar];
    xres = new double [nPoint*nVar];

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

	}

	/*--- Read farfield conditions from config ---*/
	Density_Inf   = config->GetDensity_FreeStreamND();
	Viscosity_Inf = config->GetViscosity_FreeStreamND();

	/*--- Factor_nu_Inf in [3.0, 5.0] ---*/
	Factor_nu_Inf = config->GetNuFactor_FreeStream();
	nu_tilde_Inf  = Factor_nu_Inf*Viscosity_Inf/Density_Inf;

	/*--- Eddy viscosity ---*/
	double Ji, Ji_3, fv1, cv1_3 = 7.1*7.1*7.1;
	double muT_Inf;
	Ji = nu_tilde_Inf/Viscosity_Inf*Density_Inf;
	Ji_3 = Ji*Ji*Ji;
	fv1 = Ji_3/(Ji_3+cv1_3);
	muT_Inf = Density_Inf*fv1*nu_tilde_Inf;

	/*--- Restart the solution from file information ---*/
	if (!restart || geometry->GetFinestMGLevel() == false) {
			for (iPoint = 0; iPoint < nPoint; iPoint++)
				node[iPoint] = new CTurbSAVariable(nu_tilde_Inf, muT_Inf, nDim, nVar, config);
	}
	else {

		/*--- Restart the solution from file information ---*/
		ifstream restart_file;
		string filename = config->GetSolution_FlowFileName();
		restart_file.open(filename.data(), ios::in);

		if (restart_file.fail()) {
			cout << "There is no turbulent restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
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
		for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
			Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
		}

		/*--- Read all lines in the restart file ---*/
		long iPoint_Local; unsigned long iPoint_Global = 0; string text_line;
		double U[5], Velocity2, Pressure, Temperature, Temperature_Dim, LaminarViscosity, muT = 0.0;
		double Gas_Constant = config->GetGas_ConstantND();
    
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

				if (incompressible) {
					if (nDim == 2) point_line >> index >> U[0] >> U[1] >> U[2] >> Solution[0];
					if (nDim == 3) point_line >> index >> U[0] >> U[1] >> U[2] >> U[3] >> Solution[0];	
				}
				else {
					if (nDim == 2) point_line >> index >> U[0] >> U[1] >> U[2] >> U[3] >> Solution[0];
					if (nDim == 3) point_line >> index >> U[0] >> U[1] >> U[2] >> U[3] >> U[4] >> Solution[0];
				}

				/*--- Compute the eddy viscosity from the conservative variables and nu ---*/
        if (incompressible) {
          Density = config->GetDensity_FreeStreamND();
          LaminarViscosity  = config->GetViscosity_FreeStreamND();
        }
        else {
          Density = U[0];
          Velocity2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Velocity2 += (U[iDim+1]/Density)*(U[iDim+1]/Density);
          Pressure    = Gamma_Minus_One*Density*(U[nDim+1]/Density-0.5*Velocity2);
          Temperature = Pressure / ( Gas_Constant * Density);
          /*--- Calculate lam. viscosity from a non-dim. Sutherland's Law ---*/
          Temperature_Dim  = Temperature*config->GetTemperature_Ref();
          LaminarViscosity = 1.853E-5*(pow(Temperature_Dim/300.0,3.0/2.0)
                                       *(300.0+110.3)/(Temperature_Dim+110.3));
          LaminarViscosity = LaminarViscosity/config->GetViscosity_Ref();
        }
				Ji   = Solution[0]/(LaminarViscosity/Density);
				Ji_3 = Ji*Ji*Ji;
				fv1  = Ji_3/(Ji_3+cv1_3);
				muT  = Density*fv1*Solution[0];

				/*--- Instantiate the solution at this node ---*/
				node[iPoint_Local] = new CTurbSAVariable(Solution[0], muT, nDim, nVar, config);
			}
			iPoint_Global++;
		}

		/*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CTurbSAVariable(Solution[0], muT, nDim, nVar, config);
		}

		/*--- Close the restart file ---*/
		restart_file.close();

		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
	}

  /*--- MPI solution ---*/
  SetSolution_MPI(geometry, config);
  
}

CTurbSASolution::~CTurbSASolution(void) {

}

void CTurbSASolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the residual vector ---*/
		Set_Residual_Zero(iPoint);
    
    /*--- Basic nu tilde check ---*/
    if (node[iPoint]->GetSolution(0) < 0.0)
      node[iPoint]->SetSolution(0, node[iPoint]->GetSolution_Old(0));
    
  }
  
  /*--- Initialize the jacobian matrices ---*/
	Jacobian.SetValZero();

  /*--- Upwind second order reconstruction ---*/
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
	if (config->GetKind_SlopeLimit() != NONE) SetSolution_Limiter(geometry, config);

}

void CTurbSASolution::Postprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {
	double rho, mu, nu, *nu_hat, muT, Ji, Ji_3, fv1;
	double cv1_3 = 7.1*7.1*7.1;
	unsigned long iPoint;

	bool incompressible = config->GetIncompressible();

	/*--- Compute eddy viscosity ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {

		if (incompressible) {
			rho = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
			mu  = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc();
		}
		else {
			rho = solution_container[FLOW_SOL]->node[iPoint]->GetDensity();
			mu  = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
		}

		nu  = mu/rho;
		nu_hat = node[iPoint]->GetSolution();

		Ji   = nu_hat[0]/nu;
		Ji_3 = Ji*Ji*Ji;
		fv1  = Ji_3/(Ji_3+cv1_3);

		muT = rho*fv1*nu_hat[0];
		node[iPoint]->SetmuT(muT);
    
	}
  
}

void CTurbSASolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {
	double *Turb_i, *Turb_j, *Limiter_i = NULL, *Limiter_j = NULL, *U_i, *U_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;

	bool high_order_diss = (config->GetKind_Upwind_Turb() == SCALAR_UPWIND_2ND);
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement = config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();
	bool limiter = (config->GetKind_SlopeLimit() != NONE);

	if (high_order_diss) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solution_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solution_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
	}

	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());

		/*--- Conservative variables w/o reconstruction ---*/
		U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
		U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
		solver->SetConservative(U_i, U_j);

		/*--- Turbulent variables w/o reconstruction ---*/
		Turb_i = node[iPoint]->GetSolution(); Turb_j = node[jPoint]->GetSolution();
		solver->SetTurbVar(Turb_i,Turb_j);

		/*--- Incompressible density w/o reconstruction ---*/
		if (incompressible)
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), solution_container[FLOW_SOL]->node[jPoint]->GetDensityInc());

		/*--- Rotating Frame ---*/
		if (rotating_frame)
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());

		/*--- Grid Movement ---*/
		if (grid_movement)
			solver->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());

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

			/*--- Turbulent variables using gradient reconstruction ---*/
			Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
			if (limiter) { Limiter_i = node[iPoint]->GetLimiter(); Limiter_j = node[jPoint]->GetLimiter(); }

			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				if (limiter) {
					Solution_i[iVar] = Turb_i[iVar] + Project_Grad_i*Limiter_i[iVar];
					Solution_j[iVar] = Turb_j[iVar] + Project_Grad_j*Limiter_j[iVar];
				}
				else {
					Solution_i[iVar] = Turb_i[iVar] + Project_Grad_i;
					Solution_j[iVar] = Turb_j[iVar] + Project_Grad_j;
				}
			}
			solver->SetTurbVar(Solution_i, Solution_j);
		}

		/*--- Add and subtract Residual ---*/
		solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
		AddResidual(iPoint, Residual);
		SubtractResidual(jPoint, Residual);

		/*--- Implicit part ---*/
		Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
		Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
		Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
	}

}

void CTurbSASolution::Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint;
  
	bool incompressible = config->GetIncompressible();

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Points coordinates, and normal vector ---*/
    solver->SetCoord(geometry->node[iPoint]->GetCoord(),
                     geometry->node[jPoint]->GetCoord());
    solver->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction ---*/
    solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                            solution_container[FLOW_SOL]->node[jPoint]->GetSolution());
    
    /*--- Laminar Viscosity, and density (incompresible solver)  ---*/
    if (incompressible) {
      solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
                                  solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosityInc());
      solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                            solution_container[FLOW_SOL]->node[jPoint]->GetDensityInc());
    } else {
      solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                  solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
    }
    
    /*--- Eddy Viscosity ---*/
    solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
                             solution_container[FLOW_SOL]->node[jPoint]->GetEddyViscosity());
    
    /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
    solver->SetTurbVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    solver->SetTurbVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
    
    /*--- Compute residual, and Jacobians ---*/
    solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    /*--- Add and subtract residual, and update Jacobians ---*/
    SubtractResidual(iPoint, Residual);
    AddResidual(jPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    
  }

}

void CTurbSASolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CNumerics *second_solver,
		CConfig *config, unsigned short iMesh) {
	unsigned long iPoint;
  double LevelSet;
  unsigned short iVar;

	bool incompressible = config->GetIncompressible();
	bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
	bool transition = (config->GetKind_Trans_Model() == LM);
	bool freesurface = ((config->GetKind_Solver() == FREE_SURFACE_RANS) ||
                      (config->GetKind_Solver() == ADJ_FREE_SURFACE_RANS));
  double epsilon          = config->GetFreeSurface_Thickness();
  
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

		/*--- Conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);

		/*--- Gradient of the primitive and conservative variables ---*/
		solver->SetPrimVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);

		/*--- Laminar viscosity and density (incompressible solver) ---*/
		if (incompressible) {
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(), 0.0);
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 0.0);
		} else {
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
		}

		if (transition) {
			solver->SetIntermittency(solution_container[TRANS_SOL]->node[iPoint]->GetIntermittency() );
		}

		/*--- Turbulent variables w/o reconstruction, and its gradient ---*/
		solver->SetTurbVar(node[iPoint]->GetSolution(), NULL);
		solver->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);

		/*--- Set volume ---*/
		solver->SetVolume(geometry->node[iPoint]->GetVolume());

		/*--- Set distance to the surface ---*/
		solver->SetDistance(geometry->node[iPoint]->GetWallDistance(), 0.0);

		/*--- Compute the source term ---*/
		solver->SetResidual(Residual, Jacobian_i, NULL, config);

    /*--- Don't add source term in the interface or air ---*/
    if (freesurface) {
      LevelSet = solution_container[LEVELSET_SOL]->node[iPoint]->GetSolution(0);
      if (LevelSet > -epsilon) for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
    }
    
		/*--- Subtract residual and the jacobian ---*/
		SubtractResidual(iPoint, Residual);

		Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    
	}

	if (time_spectral) {

		double Volume, Source;
		unsigned short nVar_Turb = solution_container[TURB_SOL]->GetnVar();

		/*--- Loop over points ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

			/*--- Get control volume ---*/
			Volume = geometry->node[iPoint]->GetVolume();

			/*--- Access stored time spectral source term ---*/
			for (unsigned short iVar = 0; iVar < nVar_Turb; iVar++) {
				Source = node[iPoint]->GetTimeSpectral_Source(iVar);
				Residual[iVar] = Source*Volume;
			}

			/*--- Add Residual ---*/
			AddResidual(iPoint, Residual);

		}
	}
  
}


void CTurbSASolution::Source_Template(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh) {

}

void CTurbSASolution::BC_HeatFlux_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar;

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Get the velocity vector ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Solution[iVar] = 0.0;

			node[iPoint]->SetSolution_Old(Solution);
			Set_Residual_Zero(iPoint);

			/*--- includes 1 in the diagonal ---*/
			Jacobian.DeleteValsRowi(iPoint);
		}
	}

}

void CTurbSASolution::BC_Isothermal_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config,
                                         unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar;
  
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			/*--- Get the velocity vector ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Solution[iVar] = 0.0;
      
			node[iPoint]->SetSolution_Old(Solution);
			Set_Residual_Zero(iPoint);
      
			/*--- includes 1 in the diagonal ---*/
			Jacobian.DeleteValsRowi(iPoint);
		}
	}
  
}

void CTurbSASolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	double *Normal;

	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement	= config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();
	bool gravity        = config->GetGravityForce();

	Normal = new double[nDim];

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

			/*--- Retrieve solution at the farfield boundary node ---*/
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);

			/*--- Construct solution state at infinity (farfield) ---*/
			if (incompressible) {
				double yFreeSurface = config->GetFreeSurface_Zero();
				double PressFreeSurface = solution_container[FLOW_SOL]->GetPressure_Inf();
				double Density = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
				double Froude = config->GetFroude();
				double yCoord = geometry->node[iPoint]->GetCoord(1);
				double Pressure;
				if (gravity) Pressure = PressFreeSurface + Density*(yFreeSurface-yCoord)/(Froude*Froude);
				else Pressure = PressFreeSurface;

				FlowSolution_j[0] = Pressure;
				FlowSolution_j[1] = solution_container[FLOW_SOL]->GetVelocity_Inf(0)*Density;
				FlowSolution_j[2] = solution_container[FLOW_SOL]->GetVelocity_Inf(1)*Density;
				if (nDim == 3) FlowSolution_j[3] = solution_container[FLOW_SOL]->GetVelocity_Inf(2)*Density;			
			}
			else {
				FlowSolution_j[0] = solution_container[FLOW_SOL]->GetDensity_Inf(); 
				FlowSolution_j[nDim+1] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
				for (iDim = 0; iDim < nDim; iDim++)
					FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(iDim);
			}

			/*--- Rotating Frame ---*/
			if (rotating_frame)
				conv_solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());

			/*--- Grid Movement ---*/
			if (grid_movement)
				conv_solver->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

			/*--- Incompressible --*/
			if (incompressible)
				conv_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 
						solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());

			conv_solver->SetConservative(FlowSolution_i, FlowSolution_j); 

			/*--- Set turbulent variable at the wall, and at infinity ---*/
			for (iVar = 0; iVar < nVar; iVar++) 
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
			Solution_j[0] = nu_tilde_Inf;
			conv_solver->SetTurbVar(Solution_i, Solution_j);

			/*--- Set Normal (it is necessary to change the sign) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = -Normal[iDim]; 
			conv_solver->SetNormal(Normal);

			/*--- Compute residuals and jacobians ---*/
			conv_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);

			/*--- Add residuals and jacobians ---*/
			AddResidual(iPoint, Residual);
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

		}
	}

	delete [] Normal; 

}

void CTurbSASolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned short iVar, iDim, Kind_Inlet = config->GetKind_Inlet();
	unsigned long iVertex, iPoint, Point_Normal;
	double P_Total, T_Total, Velocity[3];
	double Velocity2, H_Total, Temperature, Riemann;
	double Pressure, Density, Energy, *Flow_Dir, Mach2;
	double SoundSpeed2, SoundSpeed_Total2, Vel_Mag;
	double alpha, aa, bb, cc, dd;
	double Two_Gamma_M1 = 2.0/Gamma_Minus_One;
	double Gas_Constant = config->GetGas_ConstantND();
	double *Normal = new double[nDim];
    
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();
    bool freesurface = config->GetFreeSurface();
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
  
	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
			double Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
      
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;
      
			/*--- Current conservative variables at this boundary node (U_domain) ---*/
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
      
			/*--- Build the fictitious intlet state based on characteristics ---*/
			if (incompressible) {
        
				/*--- Pressure computation using the internal value ---*/
				FlowSolution_j[0] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(0);
        
				/*--- The velocity is computed from the interior, and normal
         derivative for the density ---*/
				for (iDim = 0; iDim < nDim; iDim++)
					FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetVelocity_Inf(iDim)*solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
        
				/*--- The y/z velocity is interpolated due to the
         free surface effect on the pressure ---*/
				if (freesurface) FlowSolution_j[nDim] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(nDim);
        
			}
			else {
        
				/*--- Subsonic inflow: there is one outgoing characteristic (u-c),
         therefore we can specify all but one state variable at the inlet.
         The outgoing Riemann invariant provides the final piece of info.
         Adapted from an original implementation in the Stanford University
         multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
         written by Edwin van der Weide, last modified 04-20-2009. ---*/
        
				switch (Kind_Inlet) {
            
            /*--- Total properties have been specified at the inlet. ---*/
          case TOTAL_CONDITIONS:
            
            /*--- Retrieve the specified total conditions for this inlet. ---*/
            P_Total  = config->GetInlet_Ptotal(Marker_Tag);
            T_Total  = config->GetInlet_Ttotal(Marker_Tag);
            Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);
            
            /*--- Non-dim. the inputs if necessary. ---*/
            P_Total /= config->GetPressure_Ref();
            T_Total /= config->GetTemperature_Ref();
            
            /*--- Store primitives and set some variables for clarity. ---*/
            Density = FlowSolution_i[0];
            Velocity2 = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Velocity[iDim] = FlowSolution_i[iDim+1]/Density;
              Velocity2 += Velocity[iDim]*Velocity[iDim];
            }
            Energy      = FlowSolution_i[nVar-1]/Density;
            Pressure    = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
            H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
            SoundSpeed2 = Gamma*Pressure/Density;
            
            /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/
            Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
            for (iDim = 0; iDim < nDim; iDim++)
              Riemann += Velocity[iDim]*UnitaryNormal[iDim];
            
            /*--- Total speed of sound ---*/
            SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy
                                                            + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;
            
            /*--- Dot product of normal and flow direction. This should
             be negative due to outward facing boundary normal convention. ---*/
            alpha = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              alpha += UnitaryNormal[iDim]*Flow_Dir[iDim];
            
            /*--- Coefficients in the quadratic equation for the velocity ---*/
            aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
            bb = -1.0*Gamma_Minus_One*alpha*Riemann;
            cc =  0.5*Gamma_Minus_One*Riemann*Riemann
            -2.0*SoundSpeed_Total2/Gamma_Minus_One;
            
            /*--- Solve quadratic equation for velocity magnitude. Value must
             be positive, so the choice of root is clear. ---*/
            dd = bb*bb - 4.0*aa*cc;
            dd = sqrt(max(0.0,dd));
            Vel_Mag   = (-bb + dd)/(2.0*aa);
            Vel_Mag   = max(0.0,Vel_Mag);
            Velocity2 = Vel_Mag*Vel_Mag;
            
            /*--- Compute speed of sound from total speed of sound eqn. ---*/
            SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
            
            /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
            Mach2 = Velocity2/SoundSpeed2;
            Mach2 = min(1.0,Mach2);
            Velocity2   = Mach2*SoundSpeed2;
            Vel_Mag     = sqrt(Velocity2);
            SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
            
            /*--- Compute new velocity vector at the inlet ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];
            
            /*--- Static temperature from the speed of sound relation ---*/
            Temperature = SoundSpeed2/(Gamma*Gas_Constant);
            
            /*--- Static pressure using isentropic relation at a point ---*/
            Pressure = P_Total*pow((Temperature/T_Total),Gamma/Gamma_Minus_One);
            
            /*--- Density at the inlet from the gas law ---*/
            Density = Pressure/(Gas_Constant*Temperature);
            
            /*--- Using pressure, density, & velocity, compute the energy ---*/
            Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;
            
            /*--- Conservative variables, using the derived quantities ---*/
            FlowSolution_j[0] = Density;
            FlowSolution_j[1] = Velocity[0]*Density;
            FlowSolution_j[2] = Velocity[1]*Density;
            FlowSolution_j[3] = Energy*Density;
            if (nDim == 3) {
              FlowSolution_j[3] = Velocity[2]*Density;
              FlowSolution_j[4] = Energy*Density;
            }
            
            break;
            
            /*--- Mass flow has been specified at the inlet. ---*/
          case MASS_FLOW:
            
            /*--- Retrieve the specified mass flow for the inlet. ---*/
            Density  = config->GetInlet_Ttotal(Marker_Tag);
            Vel_Mag  = config->GetInlet_Ptotal(Marker_Tag);
            Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);
            
            /*--- Non-dim. the inputs if necessary. ---*/
            Density /= config->GetDensity_Ref();
            Vel_Mag /= config->GetVelocity_Ref();
            
            /*--- Get primitives from current inlet state. ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              Velocity[iDim] = node[iPoint]->GetVelocity(iDim, incompressible);
            Pressure    = node[iPoint]->GetPressure(incompressible);
            SoundSpeed2 = Gamma*Pressure/FlowSolution_i[0];
            
            /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/
            Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
            for (iDim = 0; iDim < nDim; iDim++)
              Riemann += Velocity[iDim]*UnitaryNormal[iDim];
            
            /*--- Speed of sound squared for fictitious inlet state ---*/
            SoundSpeed2 = Riemann;
            for (iDim = 0; iDim < nDim; iDim++)
              SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitaryNormal[iDim];
            
            SoundSpeed2 = max(0.0,0.5*Gamma_Minus_One*SoundSpeed2);
            SoundSpeed2 = SoundSpeed2*SoundSpeed2;
            
            /*--- Pressure for the fictitious inlet state ---*/
            Pressure = SoundSpeed2*Density/Gamma;
            
            /*--- Energy for the fictitious inlet state ---*/
            Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Vel_Mag*Vel_Mag;
            
            /*--- Conservative variables, using the derived quantities ---*/
            FlowSolution_j[0] = Density;
            FlowSolution_j[1] = Vel_Mag*Flow_Dir[0]*Density;
            FlowSolution_j[2] = Vel_Mag*Flow_Dir[1]*Density;
            FlowSolution_j[3] = Energy*Density;
            if (nDim == 3) {
              FlowSolution_j[3] = Vel_Mag*Flow_Dir[2]*Density;
              FlowSolution_j[4] = Energy*Density;
            }
            
            break;
				}
			}
      
			/*--- Set the conservative variable states. Note that we really only
       need the density and momentum components so that we can compute
       the velocity for the convective part of the upwind residual. ---*/
			conv_solver->SetConservative(FlowSolution_i, FlowSolution_j);
      
			/*--- Set the turbulent variable states (prescribed for an inflow) ---*/
      Solution_i[0] = node[iPoint]->GetSolution(0);
      Solution_j[0] = nu_tilde_Inf;
      
			conv_solver->SetTurbVar(Solution_i, Solution_j);
      
			/*--- Set various other quantities in the conv_solver class ---*/
			conv_solver->SetNormal(Normal);
      
			if (incompressible)
				conv_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                                   solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			if (rotating_frame) {
				conv_solver->SetRotVel(geometry->node[iPoint]->GetRotVel(),
                               geometry->node[iPoint]->GetRotVel());
				conv_solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
			if (grid_movement)
				conv_solver->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                geometry->node[iPoint]->GetGridVel());
      
			/*--- Compute the residual using an upwind scheme ---*/
      conv_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
      AddResidual(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
			visc_solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
			visc_solver->SetNormal(Normal);
      
			/*--- Conservative variables w/o reconstruction ---*/
			visc_solver->SetConservative(FlowSolution_i, FlowSolution_j);
      
			/*--- Laminar Viscosity, and density (incompresible solver)  ---*/
			if (incompressible) {
				visc_solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
                                         solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc());
				visc_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                                   solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			} else {
				visc_solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                         solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity());
			}
      
			/*--- Eddy Viscosity ---*/
			visc_solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
                                    solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity());
      
			/*--- Turbulent variables w/o reconstruction, and its gradients ---*/
			visc_solver->SetTurbVar(Solution_i, Solution_j);
			visc_solver->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
			/*--- Compute residual, and Jacobians ---*/
			visc_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
      
			/*--- Subtract residual, and update Jacobians ---*/
			SubtractResidual(iPoint, Residual);
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
		}
	}
  
	/*--- Free locally allocated memory ---*/
	delete[] Normal;
  
}

void CTurbSASolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver,
		CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;

	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();

	double *Normal = new double[nDim];

	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {

      /*--- Index of the closest interior node ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
			/*--- Set the conservative variables, density & Velocity same as in
			 the interior, don't need to modify pressure for the convection. ---*/
			conv_solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), 
					solution_container[FLOW_SOL]->node[iPoint]->GetSolution());

			/*--- Set the turbulent variables. Here we use a Neumann BC such
			 that the turbulent variable is copied from the interior of the
			 domain to the outlet before computing the residual.
			 Solution_i --> TurbVar_internal,
			 Solution_j --> TurbVar_outlet ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
				Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
			}
			conv_solver->SetTurbVar(Solution_i, Solution_j);

			/*--- Set Normal (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = -Normal[iDim];
			conv_solver->SetNormal(Normal);

			/*--- Set various quantities in the solver class ---*/
			if (incompressible)
				conv_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 
						solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			if (rotating_frame) {
				conv_solver->SetRotVel(geometry->node[iPoint]->GetRotVel(),
						geometry->node[iPoint]->GetRotVel());
				conv_solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
			if (grid_movement)
				conv_solver->SetGridVel(geometry->node[iPoint]->GetGridVel(),
						geometry->node[iPoint]->GetGridVel());

			/*--- Compute the residual using an upwind scheme ---*/
			conv_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			AddResidual(iPoint, Residual);

			/*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
			visc_solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
			visc_solver->SetNormal(Normal);
      
			/*--- Conservative variables w/o reconstruction ---*/
			visc_solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                                   solution_container[FLOW_SOL]->node[iPoint]->GetSolution());
      
			/*--- Laminar Viscosity, and density (incompresible solver)  ---*/
			if (incompressible) {
				visc_solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
                                         solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc());
				visc_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                                   solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			} else {
				visc_solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                         solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity());
			}
      
			/*--- Eddy Viscosity ---*/
			visc_solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
                                    solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity());
      
			/*--- Turbulent variables w/o reconstruction, and its gradients ---*/
			visc_solver->SetTurbVar(Solution_i, Solution_j);
			visc_solver->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
			/*--- Compute residual, and Jacobians ---*/
			visc_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
      
			/*--- Subtract residual, and update Jacobians ---*/
			SubtractResidual(iPoint, Residual);
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

		}
	}

	/*--- Free locally allocated memory ---*/
	delete[] Normal;

}

void CTurbSASolution::BC_Nacelle_Inflow(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iDim;
  
  bool incompressible = config->GetIncompressible();

	double *Normal = new double[nDim];
  
	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
			/*--- Set the conservative variables, density & Velocity same as in
			 the interior, don't need to modify pressure for the convection. ---*/
			conv_solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                              solution_container[FLOW_SOL]->node[iPoint]->GetSolution());
      
			/*--- Set the turbulent variables. Here we use a Neumann BC such
			 that the turbulent variable is copied from the interior of the
			 domain to the outlet before computing the residual. ---*/
			conv_solver->SetTurbVar(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
      
			/*--- Set Normal (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = -Normal[iDim];
			conv_solver->SetNormal(Normal);
      
			/*--- Compute the residual using an upwind scheme ---*/
			conv_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			AddResidual(iPoint, Residual);
      
			/*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
			visc_solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
			visc_solver->SetNormal(Normal);
      
			/*--- Conservative variables w/o reconstruction ---*/
			visc_solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                                   solution_container[FLOW_SOL]->node[iPoint]->GetSolution());
      
			/*--- Laminar Viscosity, and density (incompresible solver)  ---*/
			if (incompressible) {
				visc_solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
                                         solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc());
				visc_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                                   solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			} else {
				visc_solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                         solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity());
			}
      
			/*--- Eddy Viscosity ---*/
			visc_solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
                               solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity());
      
			/*--- Turbulent variables w/o reconstruction, and its gradients ---*/
			visc_solver->SetTurbVar(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
			visc_solver->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
			/*--- Compute residual, and Jacobians ---*/
			visc_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
      
			/*--- Subtract residual, and update Jacobians ---*/
			SubtractResidual(iPoint, Residual);
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
		}
	}
  
	/*--- Free locally allocated memory ---*/
	delete[] Normal;
  
}

void CTurbSASolution::BC_Nacelle_Exhaust(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned short iVar, iDim;
	unsigned long iVertex, iPoint, Point_Normal;
	double P_Total, T_Total, Velocity[3];
	double Velocity2, H_Total, Temperature, Riemann;
	double Pressure, Density, Energy, Mach2;
	double SoundSpeed2, SoundSpeed_Total2, Vel_Mag;
	double alpha, aa, bb, cc, dd;
	double Gas_Constant = config->GetGas_ConstantND();
	double *Normal = new double[nDim];
	double *Flow_Dir = new double[nDim];
  
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
  bool incompressible = config->GetIncompressible();

	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
			double Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
      
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;
      
			/*--- Current conservative variables at this boundary node (U_domain) ---*/
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
      
			/*--- Retrieve the specified total conditions for this inlet. ---*/
			P_Total  = config->GetNozzle_Ptotal(Marker_Tag);
			T_Total  = config->GetNozzle_Ttotal(Marker_Tag);
      
			/*--- Non-dim. the inputs if necessary. ---*/
			P_Total /= config->GetPressure_Ref();
			T_Total /= config->GetTemperature_Ref();
      
			/*--- Store primitives and set some variables for clarity. ---*/
			Density = FlowSolution_i[0];
			Velocity2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Velocity[iDim] = FlowSolution_i[iDim+1]/Density;
				Velocity2 += Velocity[iDim]*Velocity[iDim];
			}
			Energy      = FlowSolution_i[solution_container[FLOW_SOL]->GetnVar()-1]/Density;
			Pressure    = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
			H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
			SoundSpeed2 = Gamma*Pressure/Density;
      
			/*--- Compute the acoustic Riemann invariant that is extrapolated
			 from the domain interior. ---*/
			Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
			for (iDim = 0; iDim < nDim; iDim++)
				Riemann += Velocity[iDim]*UnitaryNormal[iDim];
      
			/*--- Total speed of sound ---*/
			SoundSpeed_Total2 = Gamma_Minus_One*(H_Total -
                                           (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;
      
			/*--- The flow direction is defined by the surface normal ---*/
			for (iDim = 0; iDim < nDim; iDim++)
				Flow_Dir[iDim] = -UnitaryNormal[iDim];
      
			/*--- Dot product of normal and flow direction. This should
			 be negative due to outward facing boundary normal convention. ---*/
			alpha = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				alpha += UnitaryNormal[iDim]*Flow_Dir[iDim];
      
			/*--- Coefficients in the quadratic equation for the velocity ---*/
			aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
			bb = -1.0*Gamma_Minus_One*alpha*Riemann;
			cc =  0.5*Gamma_Minus_One*Riemann*Riemann - 2.0*SoundSpeed_Total2/Gamma_Minus_One;
      
			/*--- Solve quadratic equation for velocity magnitude. Value must
			 be positive, so the choice of root is clear. ---*/
			dd = bb*bb - 4.0*aa*cc;
			dd = sqrt(max(0.0,dd));
			Vel_Mag   = (-bb + dd)/(2.0*aa);
			Vel_Mag   = max(0.0,Vel_Mag);
			Velocity2 = Vel_Mag*Vel_Mag;
      
			/*--- Compute speed of sound from total speed of sound eqn. ---*/
			SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
      
			/*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
			Mach2 = Velocity2/SoundSpeed2;
			Mach2 = min(1.0,Mach2);
			Velocity2   = Mach2*SoundSpeed2;
			Vel_Mag     = sqrt(Velocity2);
			SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
      
			/*--- Compute new velocity vector at the inlet ---*/
			for (iDim = 0; iDim < nDim; iDim++)
				Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];
      
			/*--- Static temperature from the speed of sound relation ---*/
			Temperature = SoundSpeed2/(Gamma*Gas_Constant);
      
			/*--- Static pressure using isentropic relation at a point ---*/
			Pressure = P_Total*pow((Temperature/T_Total),Gamma/Gamma_Minus_One);
      
			/*--- Density at the inlet from the gas law ---*/
			Density = Pressure/(Gas_Constant*Temperature);
      
			/*--- Using pressure, density, & velocity, compute the energy ---*/
			Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;
      
			/*--- Conservative variables, using the derived quantities ---*/
			FlowSolution_j[0] = Density;
			FlowSolution_j[1] = Velocity[0]*Density;
			FlowSolution_j[2] = Velocity[1]*Density;
			FlowSolution_j[3] = Energy*Density;
			if (nDim == 3) {
				FlowSolution_j[3] = Velocity[2]*Density;
				FlowSolution_j[4] = Energy*Density;
			}
      
			/*--- Set the conservative variable states. Note that we really only
       need the density and momentum components so that we can compute
       the velocity for the convective part of the upwind residual. ---*/
			conv_solver->SetConservative(FlowSolution_i, FlowSolution_j);
      
			/*--- Set the turbulent variable states (prescribed for an inflow) ---*/
      Solution_i[0] = node[iPoint]->GetSolution(0);
      Solution_j[0] = nu_tilde_Inf;
      
			conv_solver->SetTurbVar(Solution_i, Solution_j);
      
			/*--- Set various other quantities in the conv_solver class ---*/
			conv_solver->SetNormal(Normal);
      
			/*--- Compute the residual using an upwind scheme ---*/
      conv_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
      AddResidual(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
			visc_solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
			visc_solver->SetNormal(Normal);
      
			/*--- Conservative variables w/o reconstruction ---*/
			visc_solver->SetConservative(FlowSolution_i, FlowSolution_j);
      
			/*--- Laminar Viscosity, and density (incompresible solver)  ---*/
			if (incompressible) {
				visc_solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
                                         solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc());
				visc_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                                   solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			} else {
				visc_solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                         solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity());
			}
      
			/*--- Eddy Viscosity ---*/
			visc_solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
                                    solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity());
      
			/*--- Turbulent variables w/o reconstruction, and its gradients ---*/
			visc_solver->SetTurbVar(Solution_i, Solution_j);
			visc_solver->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
			/*--- Compute residual, and Jacobians ---*/
			visc_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
      
			/*--- Subtract residual, and update Jacobians ---*/
			SubtractResidual(iPoint, Residual);
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
		}
	}
  
	/*--- Free locally allocated memory ---*/
	delete[] Normal;
    
}

void CTurbSASolution::BC_Interface_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
                                          CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, jPoint;
	unsigned short iVar, iDim;
  
  double *Vector = new double[nDim];

#ifdef NO_MPI
  
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		if (geometry->node[iPoint]->GetDomain()) {
      
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetDonorPoint();
      
			if (iPoint != jPoint) {
        
				/*--- Store the solution for both points ---*/
				for (iVar = 0; iVar < nVar; iVar++) {
					Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
					Solution_j[iVar] = node[jPoint]->GetSolution(iVar);
				}
        
				/*--- Set Conservative Variables ---*/
				solver->SetTurbVar(Solution_i, Solution_j);
        
        /*--- Retrieve flow solution for both points ---*/
        for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
          FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
          FlowSolution_j[iVar] = solution_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
        }
        
				/*--- Set Flow Variables ---*/
        solver->SetConservative(FlowSolution_i, FlowSolution_j);
        
				/*--- Set the normal vector ---*/
				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = -Vector[iDim];
				solver->SetNormal(Vector);
        
				/*--- Add Residuals and Jacobians ---*/
				solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
				AddResidual(iPoint, Residual);
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        
			}
		}
	}
  
#else
  
	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	double *Conserv_Var, *Flow_Var;
	bool compute;
  
  unsigned short Buffer_Size = nVar+solution_container[FLOW_SOL]->GetnVar();
	double *Buffer_Send_U = new double [Buffer_Size];
	double *Buffer_Receive_U = new double [Buffer_Size];
  
	/*--- Do the send process, by the moment we are sending each
	 node individually, this must be changed ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
      
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
      
			if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
			else compute = true;
      
			/*--- We only send the information that belong to other boundary ---*/
			if ((jProcessor != rank) && compute) {
        
				Conserv_Var = node[iPoint]->GetSolution();
        Flow_Var = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
        
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Send_U[iVar] = Conserv_Var[iVar];
        
        for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
					Buffer_Send_U[nVar+iVar] = Flow_Var[iVar];
        
				MPI::COMM_WORLD.Bsend(Buffer_Send_U, Buffer_Size, MPI::DOUBLE, jProcessor, iPoint);
        
			}
		}
	}
  
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		if (geometry->node[iPoint]->GetDomain()) {
      
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
      
			if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
			else compute = true;
      
			if (compute) {
        
				/*--- We only receive the information that belong to other boundary ---*/
				if (jProcessor != rank) {
					MPI::COMM_WORLD.Recv(Buffer_Receive_U, Buffer_Size, MPI::DOUBLE, jProcessor, jPoint);
        }
				else {
          
					for (iVar = 0; iVar < nVar; iVar++)
						Buffer_Receive_U[iVar] = node[jPoint]->GetSolution(iVar);
          
          for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
            Buffer_Send_U[nVar+iVar] = solution_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
          
				}
        
				/*--- Store the solution for both points ---*/
				for (iVar = 0; iVar < nVar; iVar++) {
					Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
					Solution_j[iVar] = Buffer_Receive_U[iVar];
				}
        
        /*--- Set Turbulent Variables ---*/
        solver->SetTurbVar(Solution_i, Solution_j);
        
        /*--- Retrieve flow solution for both points ---*/
        for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
          FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
          FlowSolution_j[iVar] = Buffer_Receive_U[nVar + iVar];
        }
        
				/*--- Set Flow Variables ---*/
        solver->SetConservative(FlowSolution_i, FlowSolution_j);
        
				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = -Vector[iDim];
				solver->SetNormal(Vector);
        
				solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
				AddResidual(iPoint, Residual);
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        
			}
		}
	}
  
	delete[] Buffer_Send_U;
	delete[] Buffer_Receive_U;
  
#endif
  
  delete[] Vector;
  
}

void CTurbSASolution::BC_NearField_Boundary(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
                                            CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, jPoint;
	unsigned short iVar, iDim;
  
  double *Vector = new double[nDim];
  
#ifdef NO_MPI
  
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		if (geometry->node[iPoint]->GetDomain()) {
      
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetDonorPoint();
      
			if (iPoint != jPoint) {
        
				/*--- Store the solution for both points ---*/
				for (iVar = 0; iVar < nVar; iVar++) {
					Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
					Solution_j[iVar] = node[jPoint]->GetSolution(iVar);
				}
        
				/*--- Set Conservative Variables ---*/
				solver->SetTurbVar(Solution_i, Solution_j);
        
        /*--- Retrieve flow solution for both points ---*/
        for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
          FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
          FlowSolution_j[iVar] = solution_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
        }
        
				/*--- Set Flow Variables ---*/
        solver->SetConservative(FlowSolution_i, FlowSolution_j);
        
				/*--- Set the normal vector ---*/
				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = -Vector[iDim];
				solver->SetNormal(Vector);
        
				/*--- Add Residuals and Jacobians ---*/
				solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
				AddResidual(iPoint, Residual);
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        
			}
		}
	}
  
#else
  
	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	double *Conserv_Var, *Flow_Var;
	bool compute;
  
  unsigned short Buffer_Size = nVar+solution_container[FLOW_SOL]->GetnVar();
	double *Buffer_Send_U = new double [Buffer_Size];
	double *Buffer_Receive_U = new double [Buffer_Size];
  
	/*--- Do the send process, by the moment we are sending each
	 node individually, this must be changed ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
      
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
      
			if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
			else compute = true;
      
			/*--- We only send the information that belong to other boundary ---*/
			if ((jProcessor != rank) && compute) {
        
				Conserv_Var = node[iPoint]->GetSolution();
        Flow_Var = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
        
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Send_U[iVar] = Conserv_Var[iVar];
        
        for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
					Buffer_Send_U[nVar+iVar] = Flow_Var[iVar];
        
				MPI::COMM_WORLD.Bsend(Buffer_Send_U, Buffer_Size, MPI::DOUBLE, jProcessor, iPoint);
        
			}
		}
	}
  
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		if (geometry->node[iPoint]->GetDomain()) {
      
			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
			jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];
      
			if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
			else compute = true;
      
			if (compute) {
        
				/*--- We only receive the information that belong to other boundary ---*/
				if (jProcessor != rank) {
					MPI::COMM_WORLD.Recv(Buffer_Receive_U, Buffer_Size, MPI::DOUBLE, jProcessor, jPoint);
        }
				else {
          
					for (iVar = 0; iVar < nVar; iVar++)
						Buffer_Receive_U[iVar] = node[jPoint]->GetSolution(iVar);
          
          for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
            Buffer_Send_U[nVar+iVar] = solution_container[FLOW_SOL]->node[jPoint]->GetSolution(iVar);
          
				}
        
				/*--- Store the solution for both points ---*/
				for (iVar = 0; iVar < nVar; iVar++) {
					Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
					Solution_j[iVar] = Buffer_Receive_U[iVar];
				}
        
        /*--- Set Turbulent Variables ---*/
        solver->SetTurbVar(Solution_i, Solution_j);
        
        /*--- Retrieve flow solution for both points ---*/
        for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++) {
          FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
          FlowSolution_j[iVar] = Buffer_Receive_U[nVar + iVar];
        }
        
				/*--- Set Flow Variables ---*/
        solver->SetConservative(FlowSolution_i, FlowSolution_j);
        
				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = -Vector[iDim];
				solver->SetNormal(Vector);
        
				solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
				AddResidual(iPoint, Residual);
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        
			}
		}
	}
  
	delete[] Buffer_Send_U;
	delete[] Buffer_Receive_U;
  
#endif
  
  delete[] Vector;
  
}

void CTurbSASolution::SetResidual_DualTime(CGeometry *geometry, CSolution **solution_container, CConfig *config,
		unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
	unsigned short iVar, jVar;
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1;
	double Volume_nM1, Volume_n, Volume_nP1, TimeStep;

	bool grid_movement = config->GetGrid_Movement();

	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n   = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();

		/*--- Volume at time n-1 and n ---*/
		if (grid_movement) {
			Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
			Volume_n   = geometry->node[iPoint]->GetVolume_n();
			Volume_nP1 = geometry->node[iPoint]->GetVolume();
		}
		else {
			Volume_nM1 = geometry->node[iPoint]->GetVolume();
			Volume_n   = geometry->node[iPoint]->GetVolume();
			Volume_nP1 = geometry->node[iPoint]->GetVolume();
		}

		/*--- Time Step ---*/
		TimeStep = config->GetDelta_UnstTimeND();

		/*--- Compute Residual ---*/
		for(iVar = 0; iVar < nVar; iVar++) {
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				Residual[iVar] = ( U_time_nP1[iVar]*Volume_nP1 - U_time_n[iVar]*Volume_n ) / TimeStep;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				Residual[iVar] = ( 3.0*U_time_nP1[iVar]*Volume_nP1 - 4.0*U_time_n[iVar]*Volume_n
						+  1.0*U_time_nM1[iVar]*Volume_nM1 ) / (2.0*TimeStep);
		}

		/*--- Add Residual ---*/
		AddResidual(iPoint, Residual);

    /*--- Compute implicit contribution ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++)
        Jacobian_i[iVar][jVar] = 0.0;
      
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
        Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
    }
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    
  }
}

void CTurbSASolution::CalcEddyViscosity(double *val_FlowVars, double val_laminar_viscosity,
		double *val_TurbVar, double *val_eddy_viscosity) {

	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//

	//SU2_CPP2C START CTurbSASolution::CalcEddyViscosity
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *val_FlowVars *val_TurbVar val_laminar_viscosity
	//SU2_CPP2C OUTVARS *val_eddy_viscosity
	//SU2_CPP2C VARS DOUBLE Temperature_Ref Viscosity_Ref Gamma_Minus_One Gas_Constant
	//SU2_CPP2C CALL_LIST END

	//SU2_CPP2C DEFINE nDim

	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C DECL_LIST END

	double Ji, Ji_3, fv1;
	double rho, mu, nu, nu_hat;
	double cv1_3 = 7.1*7.1*7.1;

//	if (incompressible) {
//		rho = ??;
//		mu  = val_laminar_viscosity;
//	}
//	else {
		rho = val_FlowVars[0];
		mu  = val_laminar_viscosity;
//	}

	nu  = mu/rho;
	nu_hat = val_TurbVar[0];

	Ji   = nu_hat/nu;
	Ji_3 = Ji*Ji*Ji;
	fv1  = Ji_3/(Ji_3+cv1_3);

	val_eddy_viscosity[0] = rho*fv1*nu_hat;

	//SU2_CPP2C END CTurbSASolution::CalcEddyViscosity
}

CTurbSSTSolution::CTurbSSTSolution(void) : CTurbSolution() {

  /*--- Array initialization ---*/
  constants = NULL;
  
}

CTurbSSTSolution::CTurbSSTSolution(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTurbSolution() {

	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double dull_val;
	ifstream restart_file;
	string text_line;

	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	bool incompressible = config->GetIncompressible();

	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

  /*--- Array initialization ---*/
  constants = NULL;
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Dimension of the problem --> dependent of the turbulent model ---*/
	nVar = 2;
	nPoint = geometry->GetnPoint();
	nPointDomain = geometry->GetnPointDomain();
  
  /*--- Define geometry constants in the solver structure ---*/
	nDim = geometry->GetnDim();
	node = new CVariable*[nPoint];

	/*--- Single grid simulation ---*/
	if (iMesh == MESH_0) {

		/*--- Define some auxiliary vector related with the residual ---*/
		Residual = new double[nVar]; Residual_RMS = new double[nVar];
		Residual_i = new double[nVar]; Residual_j = new double[nVar];
    Residual_Max = new double[nVar]; Point_Max = new unsigned long[nVar];
    

		/*--- Define some auxiliary vector related with the solution ---*/
		Solution = new double[nVar];
		Solution_i = new double[nVar]; Solution_j = new double[nVar];

		/*--- Define some auxiliary vector related with the geometry ---*/
		Vector_i = new double[nDim]; Vector_j = new double[nDim];

		/*--- Define some auxiliary vector related with the flow solution ---*/
		FlowSolution_i = new double [nDim+2]; FlowSolution_j = new double [nDim+2];

		/*--- Jacobians and vector structures for implicit computations ---*/
    Jacobian_i = new double* [nVar];
    Jacobian_j = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new double [nVar];
      Jacobian_j[iVar] = new double [nVar];
    }
    
    /*--- Initialization of the structure of the whole Jacobian ---*/
    if (rank == MASTER_NODE) cout << "Initialize jacobian structure (SST model)." << endl;
    Initialize_SparseMatrix_Structure(&Jacobian, nVar, nVar, geometry, config);
    xsol = new double [nPoint*nVar];
    xres = new double [nPoint*nVar];

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

	}

	/* --- Initialize value for model constants --- */
	constants = new double[10];
	constants[0] = 0.85;   //sigma_k1
	constants[1] = 1.0;    //sigma_k2
	constants[2] = 0.5;    //sigma_om1
	constants[3] = 0.856;  //sigma_om2
	constants[4] = 0.075;  //beta_1
	constants[5] = 0.0828; //beta_2
	constants[6] = 0.09;   //betaStar
	constants[7] = 0.31;   //a1
	constants[8] = constants[4]/constants[6] - constants[2]*0.41*0.41/sqrt(constants[6]);  //alfa_1
	constants[9] = constants[5]/constants[6] - constants[3]*0.41*0.41/sqrt(constants[6]);  //alfa_2

	/* --- Initialize lower and upper limits--- */
	lowerlimit = new double[nVar];
	upperlimit = new double[nVar];

	lowerlimit[0] = 1.0e-10;
	upperlimit[0] = 1.0e10;

	lowerlimit[1] = 1.0e-4;
	upperlimit[1] = 1.0e15;

	/*--- Flow infinity initialization stuff ---*/
	double rhoInf, *VelInf, muLamInf, Intensity, viscRatio, muT_Inf;

	rhoInf    = config->GetDensity_FreeStreamND();
	VelInf    = config->GetVelocity_FreeStreamND();
	muLamInf  = config->GetViscosity_FreeStreamND();
	Intensity = config->GetTurbulenceIntensity_FreeStream();
	viscRatio = config->GetTurb2LamViscRatio_FreeStream();

	double VelMag = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		VelMag += VelInf[iDim]*VelInf[iDim];
	VelMag = sqrt(VelMag);

	kine_Inf  = 3.0/2.0*(VelMag*VelMag*Intensity*Intensity);
	omega_Inf = rhoInf*kine_Inf/(muLamInf*viscRatio);

	/*--- Eddy viscosity, initialized without stress limiter ---*/
	muT_Inf = rhoInf*kine_Inf/omega_Inf;

	/*--- Restart the solution from file information ---*/
	if (!restart || geometry->GetFinestMGLevel() == false) {
		for (iPoint = 0; iPoint < nPoint; iPoint++)
			node[iPoint] = new CTurbSSTVariable(kine_Inf, omega_Inf, muT_Inf, nDim, nVar, constants, config);
	}
	else {
		/*--- Restart the solution from file information ---*/
		string mesh_filename = config->GetSolution_FlowFileName();
		restart_file.open(mesh_filename.data(), ios::in);

		if (restart_file.fail()) {
			cout << "There is no turbulent restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
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
		for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
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

				if (incompressible) {
					if (nDim == 2) point_line >> index >> rhoInf >> dull_val >> dull_val >> Solution[0] >> Solution[1];
					if (nDim == 3) point_line >> index >> rhoInf >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
				}
				else {
					if (nDim == 2) point_line >> index >> rhoInf >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
					if (nDim == 3) point_line >> index >> rhoInf >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
				}

				/*--- Eddy viscosity, initialized without stress limiter ---*/
				muT_Inf = rhoInf*Solution[0]/Solution[1];

				/*--- Instantiate the solution at this node ---*/
				node[iPoint_Local] = new CTurbSSTVariable(Solution[0], Solution[1], muT_Inf, nDim, nVar, constants, config);
			}
			iPoint_Global++;
		}

		/*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CTurbSSTVariable(Solution[0], Solution[1], muT_Inf, nDim, nVar, constants, config);
		}

		/*--- Close the restart file ---*/
		restart_file.close();

		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
	}
  
  /*--- MPI solution ---*/
  SetSolution_MPI(geometry, config);

}

CTurbSSTSolution::~CTurbSSTSolution(void) {
  
	if (constants != NULL) delete [] constants;
  
}

void CTurbSSTSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < nPoint; iPoint ++){
		Set_Residual_Zero(iPoint);
	}
	Jacobian.SetValZero();

	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

	// Compute mean flow gradients
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
		solution_container[FLOW_SOL]->SetPrimVar_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
		solution_container[FLOW_SOL]->SetPrimVar_Gradient_LS(geometry, config);

	if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
		solution_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
		solution_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);

}

void CTurbSSTSolution::Postprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {
	double rho, mu, dist, omega, kine, vorticity[3], vortMag, strMag, F2, muT, zeta;
	double a1 = constants[7];
	unsigned long iPoint;

	bool incompressible = config->GetIncompressible();

	// Compute mean flow gradients
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
		solution_container[FLOW_SOL]->SetPrimVar_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
		solution_container[FLOW_SOL]->SetPrimVar_Gradient_LS(geometry, config);

	// Compute turbulence variable gradients
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
		SetSolution_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
		SetSolution_Gradient_LS(geometry, config);

	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
		// Compute vorticity and rate of strain magnitude
		solution_container[FLOW_SOL]->node[iPoint]->SetVorticity();
		vorticity[0] = solution_container[FLOW_SOL]->node[iPoint]->GetVorticity(0);
		vorticity[1] = solution_container[FLOW_SOL]->node[iPoint]->GetVorticity(1);
		vorticity[2] = solution_container[FLOW_SOL]->node[iPoint]->GetVorticity(2);
		vortMag = sqrt(vorticity[0]*vorticity[0] + vorticity[1]*vorticity[1] + vorticity[2]*vorticity[2]);
		//vort   = 0.0; //not yet implemented
		solution_container[FLOW_SOL]->node[iPoint]->SetStrainMag();
		strMag = solution_container[FLOW_SOL]->node[iPoint]->GetStrainMag();

		// Compute blending functions and cross diffusion
		if (incompressible) {
			rho  = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
			mu   = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc();
		}
		else {
			rho  = solution_container[FLOW_SOL]->node[iPoint]->GetDensity();
			mu   = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
		}
		dist = geometry->node[iPoint]->GetWallDistance();

		node[iPoint]->SetBlendingFunc(mu,dist,rho);
		F2 = node[iPoint]->GetF2blending();

		// Compute the eddy viscosity
		kine  = node[iPoint]->GetSolution(0);
		omega = node[iPoint]->GetSolution(1);

		//zeta = min(1.0/omega,a1/(vortMag*F2));
		zeta = min(1.0/omega,a1/(strMag*F2));
		//zeta = 1.0/omega;
		muT = min(max(rho*kine*zeta,0.0),1.0);
		node[iPoint]->SetmuT(muT);
	}
}

void CTurbSSTSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {
	double *Turb_i, *Turb_j, *U_i, *U_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;

	bool high_order_diss = (config->GetKind_Upwind_Turb() == SCALAR_UPWIND_2ND);
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement = config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();

	if (high_order_diss) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solution_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solution_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
		if (config->GetKind_SlopeLimit() != NONE) SetSolution_Limiter(geometry, config);
	}

	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());

		/*--- Conservative variables w/o reconstruction ---*/
		U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
		U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
		solver->SetConservative(U_i, U_j);

		/*--- Turbulent variables w/o reconstruction ---*/
		Turb_i = node[iPoint]->GetSolution();
		Turb_j = node[jPoint]->GetSolution();
		solver->SetTurbVar(Turb_i,Turb_j);

		/*--- Incompresible density w/o reconstruction ---*/
		if (incompressible)
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 
					solution_container[FLOW_SOL]->node[jPoint]->GetDensityInc());

		/*--- Rotating Frame ---*/
		if (rotating_frame)
			solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());

		/*--- Grid Movement ---*/
		if (grid_movement)
			solver->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());

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

			/*--- Turbulent variables using gradient reconstruction ---*/
			Gradient_i = node[iPoint]->GetGradient();
			Gradient_j = node[jPoint]->GetGradient();

			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				Solution_i[iVar] = Turb_i[iVar] + Project_Grad_i;
				Solution_j[iVar] = Turb_j[iVar] + Project_Grad_j;

			}
			solver->SetTurbVar(Solution_i, Solution_j);
		}

		/*--- Add and subtract Residual ---*/
		solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
		AddResidual(iPoint, Residual);
		SubtractResidual(jPoint, Residual);

		/*--- Implicit part ---*/
		Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
		Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
		Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    
	}
}

void CTurbSSTSolution::Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint;

	bool incompressible = config->GetIncompressible();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Points coordinates, and normal vector ---*/
    solver->SetCoord(geometry->node[iPoint]->GetCoord(),
                     geometry->node[jPoint]->GetCoord());
    solver->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction ---*/
    solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                            solution_container[FLOW_SOL]->node[jPoint]->GetSolution());
    
    /*--- Laminar Viscosity, and density (incompresible solver)  ---*/
    if (incompressible) {
      solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
                                  solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosityInc());
      solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                            solution_container[FLOW_SOL]->node[jPoint]->GetDensityInc());
    }
    else {
      solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                  solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
    }
    
    /*--- Eddy Viscosity ---*/
    solver->SetEddyViscosity(node[iPoint]->GetmuT(), node[jPoint]->GetmuT());
    
    /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
    solver->SetTurbVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    solver->SetTurbVarGradient(node[iPoint]->GetGradient(),
                               node[jPoint]->GetGradient());
    
    /*--- Menter's first blending function ---*/
    solver->SetF1blending(node[iPoint]->GetF1blending(),node[jPoint]->GetF1blending());
    
    /*--- Compute residual, and Jacobians ---*/
    solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    /*--- Add and subtract residual, and update Jacobians ---*/
    SubtractResidual(iPoint, Residual);
    AddResidual(jPoint, Residual);
    Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
    Jacobian.SubtractBlock(iPoint,jPoint,Jacobian_j);
    Jacobian.AddBlock(jPoint,iPoint,Jacobian_i);
    Jacobian.AddBlock(jPoint,jPoint,Jacobian_j);
  }
  
}

void CTurbSSTSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CNumerics *second_solver,
		CConfig *config, unsigned short iMesh) {
	unsigned long iPoint;
	bool incompressible = config->GetIncompressible();

	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

		/*--- Conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);

		/*--- Gradient of the primitive and conservative variables ---*/
		solver->SetPrimVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);

		/*--- Laminar viscosity, and density (incompresible solver) ---*/
		if (incompressible) {
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(), 0.0);
			solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(), 0.0);
		}
		else {
			solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
		}

		/*--- Eddy Viscosity ---*/
		solver->SetEddyViscosity(node[iPoint]->GetmuT(), 0.0);

		/*--- Turbulent variables w/o reconstruction, and its gradient ---*/
		solver->SetTurbVar(node[iPoint]->GetSolution(), NULL);
		solver->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);

		/*--- Set volume ---*/
		solver->SetVolume(geometry->node[iPoint]->GetVolume());

		/*--- Set distance to the surface ---*/
		solver->SetDistance(geometry->node[iPoint]->GetWallDistance(), 0.0);

		/*--- Menter's first blending function ---*/
		solver->SetF1blending(node[iPoint]->GetF1blending(),0.0);

		/*--- Menter's second blending function ---*/
		solver->SetF2blending(node[iPoint]->GetF2blending(),0.0);

		/*--- Rate of strain magnitude ---*/
		solver->SetStrainMag(solution_container[FLOW_SOL]->node[iPoint]->GetStrainMag(),0.0);

		/*--- Cross diffusion ---*/
		solver->SetCrossDiff(node[iPoint]->GetCrossDiff(),0.0);

		/*--- Compute the source term ---*/
		solver->SetResidual(Residual, Jacobian_i, NULL, config);

		/*--- Subtract residual and the jacobian ---*/
		SubtractResidual(iPoint, Residual);
		Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
	}
}

void CTurbSSTSolution::Source_Template(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
		CConfig *config, unsigned short iMesh) {

}

void CTurbSSTSolution::BC_HeatFlux_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, jPoint, iVertex, total_index;
	unsigned short iDim, iVar;
	double distance, density, laminar_viscosity, beta_1;

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- distance to closest neighbor ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
			distance = 0.0;
			for(iDim = 0; iDim < nDim; iDim++){
				distance += (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim))*
						(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}
			distance = sqrt(distance);

			/*--- Set wall values ---*/
			density = solution_container[FLOW_SOL]->node[jPoint]->GetSolution(0);
			beta_1 = constants[4];
			laminar_viscosity = solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();

			Solution[0] = 0.0;
			Solution[1] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);
			//Solution[1] = 60.0*laminar_viscosity/(beta_1*distance*distance);

			/*--- Set the solution values and zero the residual ---*/
			node[iPoint]->SetSolution_Old(Solution);
			node[iPoint]->SetSolution(Solution);
			Set_Residual_Zero(iPoint);

			/*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				Jacobian.DeleteValsRowi(total_index);
			}

		}
	}
}

void CTurbSSTSolution::BC_Isothermal_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config,
                                          unsigned short val_marker) {
	unsigned long iPoint, jPoint, iVertex, total_index;
	unsigned short iDim, iVar;
	double distance, density, laminar_viscosity, beta_1;
  
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			/*--- distance to closest neighbor ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
			distance = 0.0;
			for(iDim = 0; iDim < nDim; iDim++){
				distance += (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim))*
        (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}
			distance = sqrt(distance);
      
			/*--- Set wall values ---*/
			density = solution_container[FLOW_SOL]->node[jPoint]->GetSolution(0);
			beta_1 = constants[4];
			laminar_viscosity = solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();
      
			Solution[0] = 0.0;
			Solution[1] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);
			//Solution[1] = 60.0*laminar_viscosity/(beta_1*distance*distance);
      
			/*--- Set the solution values and zero the residual ---*/
			node[iPoint]->SetSolution_Old(Solution);
			node[iPoint]->SetSolution(Solution);
			Set_Residual_Zero(iPoint);
      
			/*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				Jacobian.DeleteValsRowi(total_index);
			}
      
		}
	}
}

void CTurbSSTSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	double *Normal;
	unsigned short iVar, iDim;

	Normal = new double[nDim];

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Set conservative variables at the wall, and at infinity ---*/
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);

			FlowSolution_j[0] = solution_container[FLOW_SOL]->GetDensity_Inf();
			FlowSolution_j[nDim+1] = solution_container[FLOW_SOL]->GetDensity_Energy_Inf();
			for (iDim = 0; iDim < nDim; iDim++)
				FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetDensity_Velocity_Inf(iDim);

			conv_solver->SetConservative(FlowSolution_i, FlowSolution_j);

			/*--- Set turbulent variable at the wall, and at infinity ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
			Solution_j[0] = kine_Inf;
			Solution_j[1] = omega_Inf;
			conv_solver->SetTurbVar(Solution_i, Solution_j);

			/*--- Set Normal (it is necessary to change the sign) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = -Normal[iDim];
			conv_solver->SetNormal(Normal);

			/*--- Compute residuals and jacobians ---*/
			conv_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);

			/*--- Add residuals and jacobians ---*/
			AddResidual(iPoint, Residual);
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}

	delete [] Normal;
}

void CTurbSSTSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config,
		unsigned short val_marker) {

	/*--- Local variables and initialization. ---*/
	unsigned short iVar, iDim, Kind_Inlet = config->GetKind_Inlet();

	unsigned long iVertex, iPoint, Point_Normal;

	double P_Total, T_Total, Velocity[3];
	double Velocity2, H_Total, Temperature, Riemann;
	double Pressure, Density, Energy, *Flow_Dir, Mach2;
	double SoundSpeed2, SoundSpeed_Total2, Vel_Mag;
	double alpha, aa, bb, cc, dd;
	double Two_Gamma_M1 = 2.0/Gamma_Minus_One;
	double Gas_Constant = config->GetGas_ConstantND();
	double *Normal = new double[nDim];

	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();

	string Marker_Tag = config->GetMarker_All_Tag(val_marker);

	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

			double Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;

			/*--- Current conservative variables at this boundary node (U_domain) ---*/
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);

			/*--- Build the fictitious intlet state based on characteristics ---*/
			if (incompressible) {

				/*--- Index of the closest interior node ---*/
				Point_Normal = iPoint; //geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

				/*--- Pressure computation using the internal value ---*/
				FlowSolution_j[0] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(0);

				/*--- The velocity is computed from the interior, and normal derivative for the density ---*/
				for (iDim = 0; iDim < nDim; iDim++)
					FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetVelocity_Inf(iDim)*
					solution_container[FLOW_SOL]->node[Point_Normal]->GetDensityInc();

				/*--- Note that the y velocity is recomputed due to the
                free surface effect on the pressure ---*/
				FlowSolution_j[nDim] = solution_container[FLOW_SOL]->node[Point_Normal]->GetSolution(nDim);

			}
			else {

				/*--- Subsonic inflow: there is one outgoing characteristic (u-c),
                therefore we can specify all but one state variable at the inlet.
                The outgoing Riemann invariant provides the final piece of info.
                Adapted from an original implementation in the Stanford University
                multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
                written by Edwin van der Weide, last modified 04-20-2009. ---*/

				switch (Kind_Inlet) {

				/*--- Total properties have been specified at the inlet. ---*/
				case TOTAL_CONDITIONS:

					/*--- Retrieve the specified total conditions for this inlet. ---*/
					P_Total  = config->GetInlet_Ptotal(Marker_Tag);
					T_Total  = config->GetInlet_Ttotal(Marker_Tag);
					Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

					/*--- Non-dim. the inputs if necessary. ---*/
					P_Total /= config->GetPressure_Ref();
					T_Total /= config->GetTemperature_Ref();

					/*--- Store primitives and set some variables for clarity. ---*/
					Density = FlowSolution_i[0];
					Velocity2 = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						Velocity[iDim] = FlowSolution_i[iDim+1]/Density;
						Velocity2 += Velocity[iDim]*Velocity[iDim];
					}
					Energy      = FlowSolution_i[nVar-1]/Density;
					Pressure    = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
					H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
					SoundSpeed2 = Gamma*Pressure/Density;

					/*--- Compute the acoustic Riemann invariant that is extrapolated
                    from the domain interior. ---*/
					Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
					for (iDim = 0; iDim < nDim; iDim++)
						Riemann += Velocity[iDim]*UnitaryNormal[iDim];

					/*--- Total speed of sound ---*/
					SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy
							+ Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

					/*--- Dot product of normal and flow direction. This should
                    be negative due to outward facing boundary normal convention. ---*/
					alpha = 0.0;
					for (iDim = 0; iDim < nDim; iDim++)
						alpha += UnitaryNormal[iDim]*Flow_Dir[iDim];

					/*--- Coefficients in the quadratic equation for the velocity ---*/
					aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
					bb = -1.0*Gamma_Minus_One*alpha*Riemann;
					cc =  0.5*Gamma_Minus_One*Riemann*Riemann
							-2.0*SoundSpeed_Total2/Gamma_Minus_One;

					/*--- Solve quadratic equation for velocity magnitude. Value must
                    be positive, so the choice of root is clear. ---*/
					dd = bb*bb - 4.0*aa*cc;
					dd = sqrt(max(0.0,dd));
					Vel_Mag   = (-bb + dd)/(2.0*aa);
					Vel_Mag   = max(0.0,Vel_Mag);
					Velocity2 = Vel_Mag*Vel_Mag;

					/*--- Compute speed of sound from total speed of sound eqn. ---*/
					SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

					/*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
					Mach2 = Velocity2/SoundSpeed2;
					Mach2 = min(1.0,Mach2);
					Velocity2   = Mach2*SoundSpeed2;
					Vel_Mag     = sqrt(Velocity2);
					SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

					/*--- Compute new velocity vector at the inlet ---*/
					for (iDim = 0; iDim < nDim; iDim++)
						Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

					/*--- Static temperature from the speed of sound relation ---*/
					Temperature = SoundSpeed2/(Gamma*Gas_Constant);

					/*--- Static pressure using isentropic relation at a point ---*/
					Pressure = P_Total*pow((Temperature/T_Total),Gamma/Gamma_Minus_One);

					/*--- Density at the inlet from the gas law ---*/
					Density = Pressure/(Gas_Constant*Temperature);

					/*--- Using pressure, density, & velocity, compute the energy ---*/
					Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;

					/*--- Conservative variables, using the derived quantities ---*/
					FlowSolution_j[0] = Density;
					FlowSolution_j[1] = Velocity[0]*Density;
					FlowSolution_j[2] = Velocity[1]*Density;
					FlowSolution_j[3] = Energy*Density;
					if (nDim == 3) {
						FlowSolution_j[3] = Velocity[2]*Density;
						FlowSolution_j[4] = Energy*Density;
					}

					break;

					/*--- Mass flow has been specified at the inlet. ---*/
				case MASS_FLOW:

					/*--- Retrieve the specified mass flow for the inlet. ---*/
					Density  = config->GetInlet_Ttotal(Marker_Tag);
					Vel_Mag  = config->GetInlet_Ptotal(Marker_Tag);
					Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

					/*--- Non-dim. the inputs if necessary. ---*/
					Density /= config->GetDensity_Ref();
					Vel_Mag /= config->GetVelocity_Ref();

					/*--- Get primitives from current inlet state. ---*/
					for (iDim = 0; iDim < nDim; iDim++)
						Velocity[iDim] = node[iPoint]->GetVelocity(iDim, incompressible);
					Pressure    = node[iPoint]->GetPressure(incompressible);
					SoundSpeed2 = Gamma*Pressure/FlowSolution_i[0];

					/*--- Compute the acoustic Riemann invariant that is extrapolated
                    from the domain interior. ---*/
					Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
					for (iDim = 0; iDim < nDim; iDim++)
						Riemann += Velocity[iDim]*UnitaryNormal[iDim];

					/*--- Speed of sound squared for fictitious inlet state ---*/
					SoundSpeed2 = Riemann;
					for (iDim = 0; iDim < nDim; iDim++)
						SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitaryNormal[iDim];

					SoundSpeed2 = max(0.0,0.5*Gamma_Minus_One*SoundSpeed2);
					SoundSpeed2 = SoundSpeed2*SoundSpeed2;

					/*--- Pressure for the fictitious inlet state ---*/
					Pressure = SoundSpeed2*Density/Gamma;

					/*--- Energy for the fictitious inlet state ---*/
					Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Vel_Mag*Vel_Mag;

					/*--- Conservative variables, using the derived quantities ---*/
					FlowSolution_j[0] = Density;
					FlowSolution_j[1] = Vel_Mag*Flow_Dir[0]*Density;
					FlowSolution_j[2] = Vel_Mag*Flow_Dir[1]*Density;
					FlowSolution_j[3] = Energy*Density;
					if (nDim == 3) {
						FlowSolution_j[3] = Vel_Mag*Flow_Dir[2]*Density;
						FlowSolution_j[4] = Energy*Density;
					}

					break;
				}
			}

			/*--- Set the conservative variable states. Note that we really only
            need the density and momentum components so that we can compute
            the velocity for the convective part of the upwind residual. ---*/
			conv_solver->SetConservative(FlowSolution_i,FlowSolution_j);

			/*--- Set the turbulent variable states. Use free-stream SST
            values for the turbulent state at the inflow. ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      
			Solution_j[0]= kine_Inf;
			Solution_j[1]= omega_Inf;

			conv_solver->SetTurbVar(Solution_i,Solution_j);

			/*--- Set various other quantities in the solver class ---*/
			conv_solver->SetNormal(Normal);
			if (incompressible)
				conv_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
						solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			if (rotating_frame) {
				conv_solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
				conv_solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}

			if (grid_movement)
				conv_solver->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

			/*--- Compute the residual using an upwind scheme ---*/
			conv_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			AddResidual(iPoint, Residual);

			/*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
		}
	}

	/*--- Free locally allocated memory ---*/
	delete[] Normal;
}

void CTurbSSTSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver,
		CConfig *config, unsigned short val_marker) {

	/*--- Local variables and initialization. ---*/
	unsigned long iPoint, iVertex;
	unsigned short iVar, iDim;

	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();

	double *Normal = new double[nDim];

	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Set the conservative variables, density & Velocity same as in
			 the interior, don't need to modify pressure for the convection. ---*/
			conv_solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
					solution_container[FLOW_SOL]->node[iPoint]->GetSolution());

			/*--- Set the turbulent variables. Here we use a Neumann BC such
			 that the turbulent variable is copied from the interior of the
			 domain to the outlet before computing the residual.
			 Solution_i --> TurbVar_internal,
			 Solution_j --> TurbVar_outlet ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
				Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
			}
			conv_solver->SetTurbVar(Solution_i, Solution_j);

			/*--- Set Normal (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = -Normal[iDim];
			conv_solver->SetNormal(Normal);

			/*--- Set various quantities in the solver class ---*/
			if (incompressible)
				conv_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
						solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			if (rotating_frame) {
				conv_solver->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
				conv_solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}

			if (grid_movement)
				conv_solver->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

			/*--- Compute the residual using an upwind scheme ---*/
			conv_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			AddResidual(iPoint, Residual);

			/*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
		}
	}

	/*--- Free locally allocated memory ---*/
	delete[] Normal;
}

void CTurbSSTSolution::SetResidual_DualTime(CGeometry *geometry, CSolution **solution_container, CConfig *config,
		unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {

	unsigned short iVar, jVar;
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1;
	double Volume_nM1, Volume_n, Volume_nP1, TimeStep;

	bool grid_movement = config->GetGrid_Movement();

	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n   = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();

		/*--- Volume at time n-1 and n ---*/
		if (grid_movement) {
			Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
			Volume_n   = geometry->node[iPoint]->GetVolume_n();
			Volume_nP1 = geometry->node[iPoint]->GetVolume();
		}
		else {
			Volume_nM1 = geometry->node[iPoint]->GetVolume();
			Volume_n   = geometry->node[iPoint]->GetVolume();
			Volume_nP1 = geometry->node[iPoint]->GetVolume();
		}

		/*--- Time Step ---*/
		TimeStep = config->GetDelta_UnstTimeND();

		/*--- Compute Residual ---*/
		for(iVar = 0; iVar < nVar; iVar++) {
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				Residual[iVar] = ( U_time_nP1[iVar]*Volume_nP1 - U_time_n[iVar]*Volume_n ) / TimeStep;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				Residual[iVar] = ( 3.0*U_time_nP1[iVar]*Volume_nP1 - 4.0*U_time_n[iVar]*Volume_n
						+  1.0*U_time_nM1[iVar]*Volume_nM1 ) / (2.0*TimeStep);
		}

		/*--- Add Residual ---*/
		AddResidual(iPoint, Residual);

    /*--- Implicit contribution ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++)
        Jacobian_i[iVar][jVar] = 0.0;
      
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
        Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
    }
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    
  }
  
}

double* CTurbSSTSolution::GetConstants() {
	return constants;
}
