/*!
 * \file solution_direct_mean.cpp
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

CEulerSolver::CEulerSolver(void) : CSolver() {

	/*--- Array initialization ---*/
	Velocity_Inlet = NULL;
	Velocity_Outlet = NULL;
	Velocity_Back = NULL;
	CDrag_Inv = NULL;
	CLift_Inv = NULL;
	CSideForce_Inv = NULL;
	CMx_Inv = NULL;
	CMy_Inv = NULL;
	CMz_Inv = NULL;
	CFx_Inv = NULL;
	CFy_Inv = NULL;
	CFz_Inv = NULL;
	CEff_Inv = NULL;
	CMerit_Inv = NULL;
	CT_Inv = NULL;
	CQ_Inv = NULL;
	CEquivArea_Inv = NULL;
	CNearFieldOF_Inv = NULL;
	ForceInviscid = NULL;
	MomentInviscid = NULL;
	FanFace_MassFlow = NULL;
	FanFace_Pressure = NULL;
	FanFace_Mach = NULL;
  FanFace_Area = NULL;
  Exhaust_MassFlow = NULL;
  Exhaust_Area = NULL;
	p1_Und_Lapl = NULL;
	p2_Und_Lapl = NULL;
	PrimVar_i = NULL;
	PrimVar_j = NULL;
	Precon_Mat_inv = NULL;
	CPressure = NULL;
	CHeatTransfer = NULL;
	YPlus = NULL;

}

CEulerSolver::CEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
	unsigned long iPoint, index, counter_local = 0, counter_global = 0;
	unsigned short iVar, iDim, iMarker;
  double Density, Velocity2, Pressure, Temperature;

	unsigned short nZone = geometry->GetnZone();
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	bool incompressible = config->GetIncompressible();
  double Gas_Constant = config->GetGas_ConstantND();
	bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
			(config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
	roe_turkel = false;
  
	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	/*--- Array initialization ---*/
	Velocity_Inlet = NULL;
	Velocity_Outlet = NULL;
	Velocity_Back = NULL;
	CDrag_Inv = NULL;
	CLift_Inv = NULL;
	CSideForce_Inv = NULL;
	CMx_Inv = NULL;
	CMy_Inv = NULL;
	CMz_Inv = NULL;
	CFx_Inv = NULL;
	CFy_Inv = NULL;
	CFz_Inv = NULL;
	CEff_Inv = NULL;
	CMerit_Inv = NULL;
	CT_Inv = NULL;
	CQ_Inv = NULL;
	CEquivArea_Inv = NULL;
	CNearFieldOF_Inv = NULL;
	ForceInviscid = NULL;
	MomentInviscid = NULL;
	FanFace_MassFlow = NULL;
  Exhaust_MassFlow = NULL;
  Exhaust_Area = NULL;
	FanFace_Pressure = NULL;
	FanFace_Mach = NULL;
  FanFace_Area = NULL;
	p1_Und_Lapl = NULL;
	p2_Und_Lapl = NULL;
	PrimVar_i = NULL;
	PrimVar_j = NULL;
	Precon_Mat_inv = NULL;
	CPressure = NULL;
	CHeatTransfer = NULL;
	YPlus = NULL;

	/*--- Set the gamma value ---*/
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Define geometry constants in the solver structure ---*/
	nDim = geometry->GetnDim();	
	if (incompressible) { nVar = nDim + 1; nPrimVar = nDim+2; nPrimVarGrad = nDim+2; }
	else { nVar = nDim + 2; nPrimVar = nDim+5; nPrimVarGrad = nDim+3; }
	nMarker = config->GetnMarker_All();
	nPoint = geometry->GetnPoint();
	nPointDomain = geometry->GetnPointDomain();

	/*--- Allocate the node variables ---*/
	node = new CVariable*[nPoint];

	/*--- Define some auxiliary vectors related to the residual ---*/
	Residual      = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
	Residual_RMS  = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
	Residual_Max  = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
	Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
	Residual_i    = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
	Residual_j    = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
	Res_Conv      = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
	Res_Visc      = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
	Res_Sour      = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;

	/*--- Define some auxiliary vectors related to the solution ---*/
	Solution   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
	Solution_i = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
	Solution_j = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;

	/*--- Define some auxiliary vectors related to the geometry ---*/
	Vector   = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
	Vector_i = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
	Vector_j = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;

	/*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
		p1_Und_Lapl = new double [nPoint]; 
		p2_Und_Lapl = new double [nPoint]; 
	}

	if ((config->GetKind_Upwind_Flow() == ROE_TURKEL_2ND) || (config->GetKind_Upwind_Flow() == ROE_TURKEL_1ST)) {
		Precon_Mat_inv = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar ++)
			Precon_Mat_inv[iVar] = new double[nVar];
		roe_turkel = true;
	}

  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {

		/*--- Block auxiliary Jacobians ---*/
		Jacobian_i = new double* [nVar]; 
		Jacobian_j = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new double [nVar]; 
			Jacobian_j[iVar] = new double [nVar]; 
		}

		/*--- Initialization of the structure for the global Jacobian ---*/
		if (rank == MASTER_NODE) cout << "Initialize jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
    
	} else {
		if (rank == MASTER_NODE)
			cout << "Explicit scheme. No jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
	}

	/*--- Computation of gradients by least squares ---*/
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {

		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim]; 
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];

		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nDim+3];
		for (iVar = 0; iVar < nDim+3; iVar++)
			cvector[iVar] = new double [nDim];
	}

	/*--- Forces definition and coefficient in all the markers ---*/
	CPressure = new double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		CPressure[iMarker] = new double [geometry->nVertex[iMarker]];

	/*--- Non dimensional coefficients ---*/
	ForceInviscid    = new double[nDim];
	MomentInviscid   = new double[3];
	CDrag_Inv        = new double[nMarker];
	CLift_Inv        = new double[nMarker];
	CSideForce_Inv   = new double[nMarker];
	CMx_Inv          = new double[nMarker];
	CMy_Inv          = new double[nMarker];
	CMz_Inv          = new double[nMarker];
	CEff_Inv         = new double[nMarker];
	CFx_Inv          = new double[nMarker];
	CFy_Inv          = new double[nMarker];
	CFz_Inv          = new double[nMarker];

	/*--- Rotational coefficients ---*/
	CMerit_Inv       = new double[nMarker];
	CT_Inv           = new double[nMarker];
	CQ_Inv           = new double[nMarker];

	/*--- Supersonic coefficients ---*/
	CEquivArea_Inv   = new double[nMarker];
	CNearFieldOF_Inv = new double[nMarker];

	/*--- Nacelle simulation ---*/
	FanFace_MassFlow  = new double[nMarker];
  Exhaust_MassFlow  = new double[nMarker];
  Exhaust_Area  = new double[nMarker];
	FanFace_Pressure  = new double[nMarker];
	FanFace_Mach  = new double[nMarker];
  FanFace_Area = new double[nMarker];

	/*--- Init total coefficients ---*/
	Total_CDrag = 0.0;  Total_CLift = 0.0;      Total_CSideForce = 0.0;
	Total_CMx = 0.0;    Total_CMy = 0.0;        Total_CMz = 0.0;
	Total_CEff = 0.0;   Total_CEquivArea = 0.0; Total_CNearFieldOF = 0.0;
	Total_CFx = 0.0;    Total_CFy = 0.0;        Total_CFz = 0.0;
	Total_CT = 0.0;     Total_CQ = 0.0;         Total_CMerit = 0.0;
  Total_Maxq = 0.0;

	/*--- Read farfield conditions ---*/
	Density_Inf  = config->GetDensity_FreeStreamND();
	Pressure_Inf = config->GetPressure_FreeStreamND();
	Velocity_Inf = config->GetVelocity_FreeStreamND();
	Energy_Inf   = config->GetEnergy_FreeStreamND();
	Mach_Inf     = config->GetMach_FreeStreamND();

	if (dual_time && config->GetUnsteady_Farfield()) {

		/*--- Read farfield conditions ---*/
		Density_Inf  = config->GetDensity_FreeStreamND_Time(0);
		Pressure_Inf = config->GetPressure_FreeStreamND_Time(0);
		Velocity_Inf = config->GetVelocity_FreeStreamND_Time(0);
		Energy_Inf   = config->GetEnergy_FreeStreamND_Time(0);
		Mach_Inf     = config->GetMach_FreeStreamND_Time(0);

	}

	/*--- Initializate fan face pressure, fan face mach number, and mass flow rate ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		FanFace_MassFlow[iMarker] = 0.0;
        Exhaust_MassFlow[iMarker] = 0.0;
		FanFace_Mach[iMarker] = Mach_Inf;
		FanFace_Pressure[iMarker] = Pressure_Inf;
        FanFace_Area[iMarker] = 0.0;
        Exhaust_Area[iMarker] = 0.0;
	}

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

	/*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to 
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
	if (!restart || geometry->GetFinestMGLevel() == false || nZone > 1) {

		/*--- Restart the solution from infinity ---*/
		for (iPoint = 0; iPoint < nPoint; iPoint++)
			node[iPoint] = new CEulerVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
	}

	else {

		/*--- Restart the solution from file information ---*/
		ifstream restart_file;
		string filename = config->GetSolution_FlowFileName();

		/*--- Append time step for unsteady restart ---*/
		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
			char buffer[50];
			unsigned long flowIter = config->GetnExtIter() - 1;
			filename.erase (filename.end()-4, filename.end());
			if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
			if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
			if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
			if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
			if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
			string UnstExt = string(buffer);
			filename.append(UnstExt);
		}
		restart_file.open(filename.data(), ios::in);

		/*--- In case there is no restart file ---*/
		if (restart_file.fail()) {
			cout << "There is no flow restart file!! " << filename.data() << "."<< endl;
			cout << "Press any key to exit..." << endl;
			cin.get(); exit(1);
		}    

		/*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
		long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];

		/*--- First, set all indices to a negative value by default ---*/
		for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
			Global2Local[iPoint] = -1;

		/*--- Now fill array with the transform values only for local points ---*/
		for(iPoint = 0; iPoint < nPointDomain; iPoint++)
			Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;

		/*--- Read all lines in the restart file ---*/
		long iPoint_Local; unsigned long iPoint_Global = 0; string text_line;

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
				if (incompressible) {
					if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2];
					if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
				}
				else {
					if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
					if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
				}
				node[iPoint_Local] = new CEulerVariable(Solution, nDim, nVar, config);
			}
			iPoint_Global++;
		}

		/*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++)
			node[iPoint] = new CEulerVariable(Solution, nDim, nVar, config);

		/*--- Close the restart file ---*/
		restart_file.close();

		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
	}

  /*--- Check that the initial solution is physical ---*/
  if (!incompressible) {
    counter_local = 0;
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      Density = node[iPoint]->GetSolution(0);
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);
      Pressure    = Gamma_Minus_One*Density*(node[iPoint]->GetSolution(nDim+1)/Density-0.5*Velocity2);
      Temperature = Pressure / ( Gas_Constant * Density);
      if ((Pressure < 0.0) || (Temperature < 0.0)) {
        Solution[0] = Density_Inf;
        for (iDim = 0; iDim < nDim; iDim++)
          Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
        Solution[nDim+1] = Energy_Inf*Density_Inf;
        node[iPoint]->SetSolution(Solution);
        node[iPoint]->SetSolution_Old(Solution);
        counter_local++;
      }
    }
#ifndef NO_MPI
    MPI::COMM_WORLD.Reduce(&counter_local, &counter_global, 1, MPI::UNSIGNED_LONG, MPI::SUM, MASTER_NODE);
#else
    counter_global = counter_local;
#endif
    if ((rank == MASTER_NODE) && (counter_global != 0))
      cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  }
  
	/*--- Define solver parameters needed for execution of destructor ---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED ) space_centered = true;
	else space_centered = false;

	if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
	else euler_implicit = false;

	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
	else least_squares = false;

	/*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);

}

CEulerSolver::~CEulerSolver(void) {
	unsigned short iVar, iMarker;

	/*--- Array initialization ---*/
	if (Velocity_Inlet != NULL)   delete [] Velocity_Inlet;
	if (Velocity_Outlet != NULL)  delete [] Velocity_Outlet;
	if (Velocity_Back != NULL)    delete [] Velocity_Back;
	if (CDrag_Inv != NULL)        delete [] CDrag_Inv;
	if (CLift_Inv != NULL)        delete [] CLift_Inv;
	if (CSideForce_Inv != NULL)   delete [] CSideForce_Inv;
	if (CMx_Inv != NULL)          delete [] CMx_Inv;
	if (CMy_Inv != NULL)          delete [] CMy_Inv;
	if (CMz_Inv != NULL)          delete [] CMz_Inv;
	if (CFx_Inv != NULL)          delete [] CFx_Inv;
	if (CFy_Inv != NULL)          delete [] CFy_Inv;
	if (CFz_Inv != NULL)          delete [] CFz_Inv;
	if (CEff_Inv != NULL)         delete [] CEff_Inv;
	if (CMerit_Inv != NULL)       delete [] CMerit_Inv;
	if (CT_Inv != NULL)           delete [] CT_Inv;
	if (CQ_Inv != NULL)           delete [] CQ_Inv;
	if (CEquivArea_Inv != NULL)   delete [] CEquivArea_Inv;
	if (CNearFieldOF_Inv != NULL) delete [] CNearFieldOF_Inv;
	if (ForceInviscid != NULL)    delete [] ForceInviscid;
	if (MomentInviscid != NULL)   delete [] MomentInviscid;
	if (FanFace_MassFlow != NULL) delete [] FanFace_MassFlow;
	if (Exhaust_MassFlow != NULL) delete [] Exhaust_MassFlow;
  if (Exhaust_Area != NULL)     delete [] Exhaust_Area;
	if (FanFace_Pressure != NULL) delete [] FanFace_Pressure;
	if (FanFace_Mach != NULL)     delete [] FanFace_Mach;
  if (FanFace_Area != NULL)     delete [] FanFace_Area;
	if (p1_Und_Lapl != NULL)      delete [] p1_Und_Lapl;
	if (p2_Und_Lapl != NULL)      delete [] p2_Und_Lapl;
	if (PrimVar_i != NULL)        delete [] PrimVar_i;
	if (PrimVar_j != NULL)        delete [] PrimVar_j;

	if (Precon_Mat_inv != NULL) {
		for (iVar = 0; iVar < nVar; iVar ++)
			delete Precon_Mat_inv[iVar];
		delete [] Precon_Mat_inv;
	}

	if (CPressure != NULL) {
		for (iMarker = 0; iMarker < nMarker; iMarker++)
			delete CPressure[iMarker];
		delete [] CPressure;
	}

	if (CHeatTransfer != NULL) {
		for (iMarker = 0; iMarker < nMarker; iMarker++) {
			delete CHeatTransfer[iMarker];
		}
		delete [] CHeatTransfer;
	}

	if (YPlus != NULL) {
		for (iMarker = 0; iMarker < nMarker; iMarker++) {
			delete YPlus[iMarker];
		}
		delete [] YPlus;
	}

}

void CEulerSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new double [nBufferR_Vector];
      Buffer_Send_U = new double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
#ifndef NO_MPI
      
//      /*--- Send/Receive using non-blocking communications ---*/
//      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
//      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
//      send_request.Wait(status);
//      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
	}
  
}

void CEulerSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new double [nBufferR_Vector];
      Buffer_Send_U = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
      }
      
#ifndef NO_MPI
      
//      /*--- Send/Receive using non-blocking communications ---*/
//      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
//      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
//      send_request.Wait(status);
//      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Old(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
	}
}

void CEulerSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new double [nBufferR_Vector];
      Buffer_Send_Limit = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter(iVar);
      }
      
#ifndef NO_MPI
      
//      /*--- Send/Receive using non-blocking communications ---*/
//      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Limit, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
//      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Limit, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
//      send_request.Wait(status);
//      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Limit, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetLimiter(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
	}
}

void CEulerSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Send_Undivided_Laplacian = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Undivided_Laplacian = new double [nBufferR_Vector];
      Buffer_Send_Undivided_Laplacian = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Undivided_Laplacian[iVar*nVertexS+iVertex] = node[iPoint]->GetUndivided_Laplacian(iVar);
      }
      
#ifndef NO_MPI
      
//      /*--- Send/Receive using non-blocking communications ---*/
//      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
//      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
//      send_request.Wait(status);
//      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex] = Buffer_Send_Undivided_Laplacian[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Undivided_Laplacian;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetUndivided_Laplacian(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Undivided_Laplacian;
      
    }
    
	}
  
}

void CEulerSolver::Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config) {
	unsigned short iMarker, MarkerS, MarkerR, *Buffer_Receive_Neighbor = NULL, *Buffer_Send_Neighbor = NULL;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new double [nBufferR_Vector];
      Buffer_Send_Lambda = new double[nBufferS_Vector];
      Buffer_Receive_Neighbor = new unsigned short [nBufferR_Vector];
      Buffer_Send_Neighbor = new unsigned short[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetLambda();
        Buffer_Send_Neighbor[iVertex] = geometry->node[iPoint]->GetnPoint();
      }
      
#ifndef NO_MPI
      
//      /*--- Send/Receive using non-blocking communications ---*/
//      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Lambda, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
//      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Lambda, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
//      send_request.Wait(status);
//      recv_request.Wait(status);
//      
//      /*--- Send/Receive using non-blocking communications ---*/
//      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Neighbor, nBufferS_Vector, MPI::UNSIGNED_SHORT, 1, send_to);
//      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Neighbor, nBufferR_Vector, MPI::UNSIGNED_SHORT, 1, receive_from);
//      send_request.Wait(status);
//      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Lambda, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Neighbor, nBufferS_Vector, MPI::UNSIGNED_SHORT, send_to, 1,
                               Buffer_Receive_Neighbor, nBufferR_Vector, MPI::UNSIGNED_SHORT, receive_from, 1);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        Buffer_Receive_Lambda[iVertex] = Buffer_Send_Lambda[iVertex];
        Buffer_Receive_Neighbor[iVertex] = Buffer_Send_Neighbor[iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Lambda;
      delete [] Buffer_Send_Neighbor;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        node[iPoint]->SetLambda(Buffer_Receive_Lambda[iVertex]);
        geometry->node[iPoint]->SetnNeighbor(Buffer_Receive_Neighbor[iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Lambda;
      delete [] Buffer_Receive_Neighbor;
      
    }
    
	}
}

void CEulerSolver::Set_MPI_Dissipation_Switch(CGeometry *geometry, CConfig *config) {
	unsigned short iMarker, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new double [nBufferR_Vector];
      Buffer_Send_Lambda = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetSensor();
      }
      
#ifndef NO_MPI
      
//      /*--- Send/Receive using non-blocking communications ---*/
//      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Lambda, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
//      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Lambda, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
//      send_request.Wait(status);
//      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Lambda, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        Buffer_Receive_Lambda[iVertex] = Buffer_Send_Lambda[iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Lambda;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        node[iPoint]->SetSensor(Buffer_Receive_Lambda[iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Lambda;
      
    }
    
	}
}

void CEulerSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
  double **Gradient = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new double[nDim];
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new double [nBufferR_Vector];
      Buffer_Send_Gradient = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
      }
      
#ifndef NO_MPI
      
//      /*--- Send/Receive using non-blocking communications ---*/
//      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Gradient, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
//      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Gradient, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
//      send_request.Wait(status);
//      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Gradient, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
        
        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
        }
        
        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient(iVar, iDim, Gradient[iVar][iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;
      
    }
    
	}
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}

void CEulerSolver::Set_MPI_PrimVar_Gradient(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
	int send_to, receive_from;
  
#ifndef NO_MPI
  MPI::Status status;
  MPI::Request send_request, recv_request;
#endif
  
  double **Gradient = new double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Gradient[iVar] = new double[nDim];
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nPrimVarGrad*nDim;        nBufferR_Vector = nVertexR*nPrimVarGrad*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new double [nBufferR_Vector];
      Buffer_Send_Gradient = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient_Primitive(iVar, iDim);
      }
      
#ifndef NO_MPI
      
//      /*--- Send/Receive using non-blocking communications ---*/
//      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Gradient, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
//      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Gradient, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
//      send_request.Wait(status);
//      recv_request.Wait(status);
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Gradient, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
        
        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
        }
        
        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient_Primitive(iVar, iDim, Gradient[iVar][iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;
      
    }
    
	}
  
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}

void CEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
	unsigned long iPoint, Point_Fine;
	unsigned short iMesh, iChildren, iVar, iDim;
	double Density, Pressure, yFreeSurface, PressFreeSurface, Froude, yCoord, Velx, Vely, Velz, RhoVelx, RhoVely, RhoVelz, XCoord, YCoord, ZCoord, DensityInc, ViscosityInc, Heaviside, LevelSet, lambda, DensityFreeSurface, Area_Children, Area_Parent, LevelSet_Fine, epsilon, *Solution_Fine, *Solution; //, Factor_nu, nu_tilde;
  
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	unsigned short nDim = geometry[MESH_0]->GetnDim();
	bool freesurface = config->GetFreeSurface();
	bool rans = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS) ||
                 (config->GetKind_Solver() == FREE_SURFACE_RANS) || (config->GetKind_Solver() == ADJ_FREE_SURFACE_RANS));
	bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
    
	if (freesurface) {
        
		for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
            
			for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
                
				/*--- Set initial boundary condition at iter 0 ---*/
				if ((ExtIter == 0) && (!restart)) {
                    
					/*--- Compute the level set value in all the MG levels (basic case, distance to
                     the Y/Z plane, and interpolate the solution to the coarse levels ---*/
					if (iMesh == MESH_0) {
						XCoord = geometry[iMesh]->node[iPoint]->GetCoord(0);
						YCoord = geometry[iMesh]->node[iPoint]->GetCoord(1);
						if (nDim == 2) LevelSet = YCoord - config->GetFreeSurface_Zero();
						else {
							ZCoord = geometry[iMesh]->node[iPoint]->GetCoord(2);
							LevelSet = ZCoord - config->GetFreeSurface_Zero();
						}
						solver_container[iMesh][LEVELSET_SOL]->node[iPoint]->SetSolution(0, LevelSet);
					}
					else {
						Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
						LevelSet = 0.0;
						for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
							Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
							Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
							LevelSet_Fine = solver_container[iMesh-1][LEVELSET_SOL]->node[Point_Fine]->GetSolution(0);
							LevelSet += LevelSet_Fine*Area_Children/Area_Parent;
						}
						solver_container[iMesh][LEVELSET_SOL]->node[iPoint]->SetSolution(0, LevelSet);
					}
                    
					/*--- Compute the flow solution using the level set value. ---*/
					epsilon = config->GetFreeSurface_Thickness();
					Heaviside = 0.0;
					if (LevelSet < -epsilon) Heaviside = 1.0;
					if (fabs(LevelSet) <= epsilon) Heaviside = 1.0 - (0.5*(1.0+(LevelSet/epsilon)+(1.0/PI_NUMBER)*sin(PI_NUMBER*LevelSet/epsilon)));
					if (LevelSet > epsilon) Heaviside = 0.0;
                    
					/*--- Set the value of the incompressible density for free surface flows (density ratio g/l) ---*/
					lambda = config->GetRatioDensity();
					DensityInc = (lambda + (1.0 - lambda)*Heaviside)*config->GetDensity_FreeStreamND();
					solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetDensityInc(DensityInc);
                    
					/*--- Set the value of the incompressible viscosity for free surface flows (viscosity ratio g/l) ---*/
					lambda = config->GetRatioViscosity();
					ViscosityInc = (lambda + (1.0 - lambda)*Heaviside)*config->GetViscosity_FreeStreamND();
					solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetLaminarViscosityInc(ViscosityInc);
                    
					/*--- Set the value of the incompressible viscosity for free surface flows (viscosity ratio g/l) ---*/
					//          if (rans && (iMesh == MESH_0)) {
					//            Factor_nu = config->GetNuFactor_FreeStream();
					//            nu_tilde = Factor_nu*ViscosityInc/DensityInc;
					//            solver_container[iMesh][TURB_SOL]->node[iPoint]->SetSolution(0, nu_tilde);
					//          }
                    
					/*--- Update solution with the new pressure ---*/
					yFreeSurface = config->GetFreeSurface_Zero();
					PressFreeSurface = solver_container[iMesh][FLOW_SOL]->GetPressure_Inf();
					DensityFreeSurface = solver_container[iMesh][FLOW_SOL]->GetDensity_Inf();
					Density = solver_container[iMesh][FLOW_SOL]->node[iPoint]->GetDensityInc();
					Froude = config->GetFroude();
					yCoord = geometry[iMesh]->node[iPoint]->GetCoord(nDim-1);
					Pressure = PressFreeSurface + Density*((yFreeSurface-yCoord)/(Froude*Froude));
					solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(0, Pressure);
                    
					/*--- Update solution with the new velocity ---*/
					Velx = solver_container[iMesh][FLOW_SOL]->GetVelocity_Inf(0);
					Vely = solver_container[iMesh][FLOW_SOL]->GetVelocity_Inf(1);
					RhoVelx = Velx * Density; RhoVely = Vely * Density;
					solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(1, RhoVelx);
					solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(2, RhoVely);
					if (nDim == 3) {
						Velz = solver_container[iMesh][FLOW_SOL]->GetVelocity_Inf(2);
						RhoVelz = Velz * Density;
						solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(3, RhoVelz);
					}
                    
				}
                
			}
            
			/*--- Set the MPI communication ---*/
			solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
			solver_container[iMesh][LEVELSET_SOL]->Set_MPI_Solution(geometry[iMesh], config);
            
            /*--- Set the primitive variables for the level set equation ---*/
            solver_container[iMesh][LEVELSET_SOL]->SetLevelSet_Distance(geometry[iMesh], config, true, false);
            
			/*--- The value of the solution for the first iteration of the dual time ---*/
			for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
				if ((ExtIter == 0) && (dual_time)) {
					solver_container[iMesh][LEVELSET_SOL]->node[iPoint]->Set_Solution_time_n();
					solver_container[iMesh][LEVELSET_SOL]->node[iPoint]->Set_Solution_time_n1();
				}
			}
            
		}
        
	}
    
  if (config->GetEngine_Intake()) {
    
    /*--- Set initial boundary condition at iteration 0 ---*/
    if ((ExtIter == 0) && (!restart)) {
      
      double *Coord = geometry[iMesh]->node[iPoint]->GetCoord();
      double Velocity_FreeStream[3] = {0.0, 0.0, 0.0}, Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStream, Density_FreeStream, Pressure_FreeStream, Density_FreeStreamND, Pressure_FreeStreamND, ModVel_FreeStreamND, Energy_FreeStreamND, ModVel_FreeStream;
      
      double Mach = 0.40;
      double Alpha = config->GetAoA()*PI_NUMBER/180.0;
      double Beta  = config->GetAoS()*PI_NUMBER/180.0;
      double Gamma_Minus_One = Gamma - 1.0;
      double Gas_Constant = config->GetGas_ConstantND();
      double Temperature_FreeStream = config->GetTemperature_FreeStream();
      double Mach2Vel_FreeStream = sqrt(Gamma*Gas_Constant*Temperature_FreeStream);
      
      for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
        
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          
          Velocity_FreeStream[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
          Velocity_FreeStream[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
          Velocity_FreeStream[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
          
          ModVel_FreeStream = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            ModVel_FreeStream += Velocity_FreeStream[iDim]*Velocity_FreeStream[iDim];
          ModVel_FreeStream = sqrt(ModVel_FreeStream);
          
          if (config->GetViscous()) {
            Viscosity_FreeStream = 1.853E-5*(pow(Temperature_FreeStream/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_FreeStream+110.3));
            Density_FreeStream   = config->GetReynolds()*Viscosity_FreeStream/(ModVel_FreeStream*config->GetLength_Reynolds());
            Pressure_FreeStream  = Density_FreeStream*Gas_Constant*Temperature_FreeStream;
          }
          else {
            Pressure_FreeStream  = config->GetPressure_FreeStream();
            Density_FreeStream  = Pressure_FreeStream/(Gas_Constant*Temperature_FreeStream);
          }
          
          Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();
          Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();
          
          for (iDim = 0; iDim < nDim; iDim++)
            Velocity_FreeStreamND[iDim] = Velocity_FreeStream[iDim]/config->GetVelocity_Ref();
          
          ModVel_FreeStreamND = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
          ModVel_FreeStreamND = sqrt(ModVel_FreeStreamND);
          
          Energy_FreeStreamND = Pressure_FreeStreamND/(Density_FreeStreamND*Gamma_Minus_One)+0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;
          
          
          if (((Coord[0] >= 16.0) && (Coord[0] <= 20.0)) &&
              ((Coord[1] >= 0.0) && (Coord[1] <= 0.7)) &&
              ((Coord[2] >= 2.5) && (Coord[2] <= 4.0))) {
            
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(0, Density_FreeStreamND);
            for (iDim = 0; iDim < nDim; iDim++)
              solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(iDim+1, Velocity_FreeStreamND[iDim]);
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(nVar-1, Energy_FreeStreamND);
            
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(0, Density_FreeStreamND);
            for (iDim = 0; iDim < nDim; iDim++)
              solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(iDim+1, Velocity_FreeStreamND[iDim]);
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(nVar-1, Energy_FreeStreamND);
          }
          
				}
        
        /*--- Set the MPI communication ---*/
        solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
        solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution_Old(geometry[iMesh], config);

			}
      
		}
    
	}
  
	/*--- If restart solution, then interpolate the flow solution to
     all the multigrid levels, this is important with the dual time strategy ---*/
	if (restart) {
		Solution = new double[nVar];
		for (iMesh = 1; iMesh <= config->GetMGLevels(); iMesh++) {
			for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
				Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
				for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
				for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
					Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
					Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
					Solution_Fine = solver_container[iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution();
					for (iVar = 0; iVar < nVar; iVar++) {
						Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
					}
				}
				solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(Solution);
                
			}
			solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
		}
		delete [] Solution;
	}
    
    
	/*--- The value of the solution for the first iteration of the dual time ---*/
	for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
		for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
			if ((ExtIter == 0) && (dual_time)) {
				solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
				solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
				if (rans) {
					solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
					solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n1();
				}
			}
		}
	}
    
	if (dual_time && config->GetUnsteady_Farfield()) {
        
		/*--- Read farfield conditions ---*/
		Density_Inf  = config->GetDensity_FreeStreamND_Time(ExtIter);
		Pressure_Inf = config->GetPressure_FreeStreamND_Time(ExtIter);
		Velocity_Inf = config->GetVelocity_FreeStreamND_Time(ExtIter);
		Energy_Inf   = config->GetEnergy_FreeStreamND_Time(ExtIter);
		Mach_Inf     = config->GetMach_FreeStreamND_Time(ExtIter);
        
		/*--- Initializate fan face pressure, fan face mach number, and mass flow rate ---*/
		for (unsigned short iMarker = 0; iMarker < nMarker; iMarker++) {
			FanFace_MassFlow[iMarker] = 0.0;
            Exhaust_MassFlow[iMarker] = 0.0;
            Exhaust_Area[iMarker] = 0.0;
			FanFace_Mach[iMarker] = Mach_Inf;
			FanFace_Pressure[iMarker] = Pressure_Inf;
            FanFace_Area[iMarker] = 0.0;
		}
        
		/*--- Inlet/Outlet boundary conditions, using infinity values ---*/
		Density_Inlet = Density_Inf;		Density_Outlet = Density_Inf;
		Pressure_Inlet = Pressure_Inf;	Pressure_Outlet = Pressure_Inf;
		Energy_Inlet = Energy_Inf;			Energy_Outlet = Energy_Inf;
		Mach_Inlet = Mach_Inf;					Mach_Outlet = Mach_Inf;
		Velocity_Inlet  = new double [nDim]; Velocity_Outlet = new double [nDim];
		for (unsigned short iDim = 0; iDim < nDim; iDim++) {
			Velocity_Inlet[iDim] = Velocity_Inf[iDim];
			Velocity_Outlet[iDim] = Velocity_Inf[iDim];
		}
        
	}
    
    
}

void CEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
	double levelset;
  
	bool freesurface = config->GetFreeSurface();
	bool adjoint = config->GetAdjoint();
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool upwind_2nd = ((config->GetKind_Upwind_Flow() == ROE_2ND) || (config->GetKind_Upwind_Flow() == AUSM_2ND)
                     || (config->GetKind_Upwind_Flow() == HLLC_2ND) || (config->GetKind_Upwind_Flow() == ROE_TURKEL_2ND));
	bool center = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
	bool center_jst = center && (config->GetKind_Centered_Flow() == JST);
	bool low_fidelity = (config->GetLowFidelitySim() && (iMesh == MESH_1));
	bool limiter = ((config->GetKind_SlopeLimit_Flow() != NONE) && !low_fidelity);
	bool incompressible = config->GetIncompressible();
  bool jouleheating = config->GetJouleHeating();
    
    /*--- Compute nacelle inflow and exhaust properties ---*/
    GetNacelle_Properties(geometry, config, iMesh);
    
    /*--- Compute Joule heating ---*/
	if (jouleheating) geometry->SetGeometryPlanes(config);
  
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
		if (freesurface) levelset = solver_container[LEVELSET_SOL]->node[iPoint]->GetPrimVar(0);
		else levelset = 0.0;
    
		/*--- Set the primitive variables incompressible (dens, vx, vy, vz, beta)
     and compressible (temp, vx, vy, vz, press, dens, enthal, sos)---*/
		if (incompressible) node[iPoint]->SetPrimVar_Incompressible(Density_Inf, levelset, config);
		else node[iPoint]->SetPrimVar_Compressible(config);
        
		/*--- Initialize the convective residual vector ---*/
		LinSysRes.SetBlock_Zero(iPoint);
    
	}
  
	/*--- Upwind second order reconstruction ---*/
	if ((upwind_2nd) && ((iMesh == MESH_0) || low_fidelity)) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    
		/*--- Limiter computation ---*/
		if ((limiter) && (iMesh == MESH_0)) SetSolution_Limiter(geometry, config);
	}
  
	/*--- Artificial dissipation ---*/
	if (center) {
		SetMax_Eigenvalue(geometry, config);
		if ((center_jst) && ((iMesh == MESH_0) || low_fidelity)) {
			SetDissipation_Switch(geometry, config);
			SetUndivided_Laplacian(geometry, config);
		}
	}
  
	/*--- Initialize the jacobian matrices ---*/
	if (implicit) Jacobian.SetValZero();
  
}

void CEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
		unsigned short iMesh, unsigned long Iteration) {
	double *Normal, Area, Vol, Mean_SoundSpeed, Mean_ProjVel, Mean_BetaInc2, Lambda, Local_Delta_Time, Mean_DensityInc,
	Global_Delta_Time = 1E6, Global_Delta_UnstTimeND, ProjVel, ProjVel_i, ProjVel_j;
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();
	bool grid_movement = config->GetGrid_Movement();
	bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
			(config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

	Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;

	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		node[iPoint]->SetMax_Lambda_Inv(0.0);

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
			Mean_DensityInc = 0.5 * (node[iPoint]->GetDensityInc() + node[jPoint]->GetDensityInc());
			Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);
		}
		else {
			Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
			Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
		}

		/*--- Adjustment for a rotating frame ---*/
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

		/*--- Adjustment for grid movement ---*/
		if (grid_movement) {
			double *GridVel_i = geometry->node[iPoint]->GetGridVel();
			double *GridVel_j = geometry->node[jPoint]->GetGridVel();
			ProjVel_i = 0.0; ProjVel_j =0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				ProjVel_i += GridVel_i[iDim]*Normal[iDim];
				ProjVel_j += GridVel_j[iDim]*Normal[iDim];
			}
			Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
		}

		/*--- Inviscid contribution ---*/
		Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);

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
				Mean_DensityInc = node[iPoint]->GetDensityInc();
				Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area); 
			}
			else {
				Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
				Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;		
			}

			/*--- Adjustment for a rotating frame ---*/
			if (rotating_frame) {
				double *RotVel = geometry->node[iPoint]->GetRotVel();
				ProjVel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					ProjVel += RotVel[iDim]*Normal[iDim];
				Mean_ProjVel -= ProjVel;
			}

			/*--- Adjustment for grid movement ---*/
			if (grid_movement) {
				double *GridVel = geometry->node[iPoint]->GetGridVel();
				ProjVel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					ProjVel += GridVel[iDim]*Normal[iDim];
				Mean_ProjVel -= ProjVel;
			}

			/*--- Inviscid contribution ---*/
			Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
			if (geometry->node[iPoint]->GetDomain()) {
				node[iPoint]->AddMax_Lambda_Inv(Lambda);
			}

		}
	}

	/*--- Each element uses their own speed, steady state simulation ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
		Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
		Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
		Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
		node[iPoint]->SetDelta_Time(Local_Delta_Time);
	}

	/*--- Check if there is any element with only one neighbor...
   a CV that is inside another CV ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		if (geometry->node[iPoint]->GetnPoint() == 1)
			node[iPoint]->SetDelta_Time(Min_Delta_Time);
	}

	/*--- For exact time solution use the minimum delta time of the whole mesh ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
#ifndef NO_MPI
		double rbuf_time, sbuf_time;
		sbuf_time = Global_Delta_Time;
		MPI::COMM_WORLD.Reduce(&sbuf_time, &rbuf_time, 1, MPI::DOUBLE, MPI::MIN, MASTER_NODE);
		MPI::COMM_WORLD.Bcast(&rbuf_time, 1, MPI::DOUBLE, MASTER_NODE);
		MPI::COMM_WORLD.Barrier();
		Global_Delta_Time = rbuf_time;
#endif
		for(iPoint = 0; iPoint < nPointDomain; iPoint++)
			node[iPoint]->SetDelta_Time(Global_Delta_Time);
	}

	/*--- Recompute the unsteady time step for the dual time strategy 
	 if the unsteady CFL is diferent from 0 ---*/
	if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
		Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);

#ifndef NO_MPI
		double rbuf_time, sbuf_time;
		sbuf_time = Global_Delta_UnstTimeND;
		MPI::COMM_WORLD.Reduce(&sbuf_time, &rbuf_time, 1, MPI::DOUBLE, MPI::MIN, MASTER_NODE);
		MPI::COMM_WORLD.Bcast(&rbuf_time, 1, MPI::DOUBLE, MASTER_NODE);
		MPI::COMM_WORLD.Barrier();
		Global_Delta_UnstTimeND = rbuf_time;
#endif
		config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
	}

	/*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
	if (dual_time)
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
			if (!implicit) {
				Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
				/*--- Check if there is any element with only one neighbor... 
				 a CV that is inside another CV ---*/
				if (geometry->node[iPoint]->GetnPoint() == 1) Local_Delta_Time = 0.0;
				node[iPoint]->SetDelta_Time(Local_Delta_Time);
			}
		}

}

void CEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
		CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool high_order_diss = ((config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0));
	bool low_fidelity = (config->GetLowFidelitySim() && (iMesh == MESH_1));
	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();
	bool grid_movement = config->GetGrid_Movement();

	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Points in edge, set normal vectors, and number of neighbors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
		numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());

		/*--- Set conservative variables w/o reconstruction ---*/
		numerics->SetConservative(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());

		/*--- Set some precomputed flow cuantities ---*/
		if (incompressible) {
			numerics->SetDensityInc(node[iPoint]->GetDensityInc(), node[jPoint]->GetDensityInc());
			numerics->SetBetaInc2(node[iPoint]->GetBetaInc2(), node[jPoint]->GetBetaInc2());
			numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
		}
		else {
			numerics->SetPressure(node[iPoint]->GetPressure(incompressible), node[jPoint]->GetPressure(incompressible));
			numerics->SetSoundSpeed(node[iPoint]->GetSoundSpeed(), node[jPoint]->GetSoundSpeed());
			numerics->SetEnthalpy(node[iPoint]->GetEnthalpy(), node[jPoint]->GetEnthalpy());
		}

		/*--- Set the largest convective eigenvalue ---*/
		numerics->SetLambda(node[iPoint]->GetLambda(), node[jPoint]->GetLambda());

		/*--- Set undivided laplacian an pressure based sensor ---*/
		if ((high_order_diss || low_fidelity)) {
			numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(), node[jPoint]->GetUndivided_Laplacian());
			numerics->SetSensor(node[iPoint]->GetSensor(), node[jPoint]->GetSensor()); 
		}

		/*--- Rotational Frame ---*/
		if (rotating_frame) {
			numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());
			numerics->SetRotFlux(geometry->edge[iEdge]->GetRotFlux());
		}

		/*--- Grid Movement ---*/
		if (grid_movement) {
			numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
		}

		/*--- Compute residuals, and jacobians ---*/
		numerics->ComputeResidual(Res_Conv, Res_Visc, Jacobian_i, Jacobian_j, config);

		/*--- Update convective and artificial dissipation residuals ---*/
		LinSysRes.AddBlock(iPoint, Res_Conv);
		LinSysRes.SubtractBlock(jPoint, Res_Conv);
    LinSysRes.AddBlock(iPoint, Res_Visc);
    LinSysRes.SubtractBlock(jPoint, Res_Visc);

		/*--- Set implicit computation ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
			Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
			Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j); 
		}
	}
}

void CEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
		CConfig *config, unsigned short iMesh) {
	double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, 
	*U_i, *U_j, sqvel, *Limiter_i = NULL, *Limiter_j = NULL;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim, iVar;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool low_fidelity = (config->GetLowFidelitySim() && (iMesh == MESH_1));
	bool high_order_diss = (((config->GetKind_Upwind_Flow() == ROE_2ND) || (config->GetKind_Upwind_Flow() == AUSM_2ND)
			|| (config->GetKind_Upwind_Flow() == HLLC_2ND) || (config->GetKind_Upwind_Flow() == ROE_TURKEL_2ND))
			&& ((iMesh == MESH_0) || low_fidelity));
	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
//	bool gravity = config->GetGravityForce();
	bool grid_movement = config->GetGrid_Movement();
	bool limiter = ((config->GetKind_SlopeLimit_Flow() != NONE) && !low_fidelity);

	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Set a copy of the conservative variables w/o reconstruction, just in case the
     2nd order reconstruction fails ---*/
    U_i = node[iPoint]->GetSolution(); U_j = node[jPoint]->GetSolution();
    numerics->SetConservative_ZeroOrder(U_i, U_j);

		/*--- Roe Turkel preconditioning ---*/
		if (roe_turkel) {
			sqvel = 0.0;
			for (iDim = 0; iDim < nDim; iDim ++)
				sqvel += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
			numerics->SetVelocity2_Inf(sqvel);
		}

		/*--- Rotating frame ---*/
		if (rotating_frame) {
			numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[jPoint]->GetRotVel());
			numerics->SetRotFlux(geometry->edge[iEdge]->GetRotFlux());
		}

		/*--- Grid Movement ---*/
		if (grid_movement)
			numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());

		/*--- High order reconstruction using MUSCL strategy ---*/
		if (high_order_diss) {

			for (iDim = 0; iDim < nDim; iDim++) {
				Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
				Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}

			Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
			if (limiter) { Limiter_j = node[jPoint]->GetLimiter(); Limiter_i = node[iPoint]->GetLimiter(); }

			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				if (limiter) {
					Solution_i[iVar] = U_i[iVar] + Project_Grad_i*Limiter_i[iVar];
					Solution_j[iVar] = U_j[iVar] + Project_Grad_j*Limiter_j[iVar];
				}
				else {
					Solution_i[iVar] = U_i[iVar] + Project_Grad_i;
					Solution_j[iVar] = U_j[iVar] + Project_Grad_j;
				}
			}

			/*--- Set conservative variables with reconstruction ---*/
			numerics->SetConservative(Solution_i, Solution_j);
      
		} else {
      
      /*--- Set conservative variables without reconstruction ---*/
      numerics->SetConservative(U_i, U_j);
      
    }

//		if (gravity) {
//      double YDistance, GradHidrosPress;
//
//			/*--- The zero order reconstruction includes the gradient of the hydrostatic pressure constribution ---*/
//			YDistance = 0.5*(geometry->node[jPoint]->GetCoord(nDim-1)-geometry->node[iPoint]->GetCoord(nDim-1));
//			GradHidrosPress = node[iPoint]->GetDensityInc()/(config->GetFroude()*config->GetFroude());
//			Solution_i[0] = U_i[0] - GradHidrosPress*YDistance;
//			GradHidrosPress = node[jPoint]->GetDensityInc()/(config->GetFroude()*config->GetFroude());
//			Solution_j[0] = U_j[0] + GradHidrosPress*YDistance;
//
//			/*--- No reconstruction for the velocities ---*/
//			for (iVar = 1; iVar < nVar; iVar++) {
//				Solution_i[iVar] = U_i[iVar];
//				Solution_j[iVar] = U_j[iVar];
//			}
//
//			/*--- Set conservative variables with reconstruction only for the pressure ---*/
//			numerics->SetConservative(Solution_i, Solution_j);
//
//			if (high_order_diss) {
//
//				for (iDim = 0; iDim < nDim; iDim++) {
//					Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
//					Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
//				}
//
//				Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
//				if (limiter) { Limiter_j = node[jPoint]->GetLimiter(); Limiter_i = node[iPoint]->GetLimiter(); }
//
//				for (iVar = 0; iVar < nVar; iVar++) {
//					Project_Grad_i = 0; Project_Grad_j = 0;
//					for (iDim = 0; iDim < nDim; iDim++) {
//						Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
//						Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
//					}
//					if (limiter) {
//						/*--- Note that the pressure reconstruction always includes the hydrostatic gradient,
//             and limits the kinematic contribution ---*/
//						if (iVar == 0) {
//							Solution_i[iVar] = Solution_i[iVar] + Limiter_i[iVar]*(U_i[iVar] + Project_Grad_i - Solution_i[iVar]);
//							Solution_j[iVar] = Solution_j[iVar] + Limiter_j[iVar]*(U_j[iVar] + Project_Grad_j - Solution_j[iVar]);
//						}
//						else {
//							Solution_i[iVar] = U_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
//							Solution_j[iVar] = U_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
//						}
//					}
//					else {
//						Solution_i[iVar] = U_i[iVar] + Project_Grad_i;
//						Solution_j[iVar] = U_j[iVar] + Project_Grad_j;
//					}
//				}
//
//				/*--- Set conservative variables with reconstruction ---*/
//				numerics->SetConservative(Solution_i, Solution_j);
//			}
//
//		}

		/*--- Set the density and beta w/o reconstruction (incompressible flows) ---*/
		if (incompressible) {
			numerics->SetBetaInc2(node[iPoint]->GetBetaInc2(), node[jPoint]->GetBetaInc2());
			numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
			numerics->SetDensityInc(node[iPoint]->GetDensityInc(), node[jPoint]->GetDensityInc());
		}

		/*--- Compute the residual ---*/
		numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);

		/*--- Update residual value ---*/
		LinSysRes.AddBlock(iPoint, Res_Conv);
		LinSysRes.SubtractBlock(jPoint, Res_Conv);

		/*--- Set implicit jacobians ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
			Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
			Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
		}

		/*--- Roe Turkel preconditioning, set the value of beta ---*/
		if (roe_turkel) {
			node[iPoint]->SetPreconditioner_Beta(numerics->GetPrecond_Beta());
			node[jPoint]->SetPreconditioner_Beta(numerics->GetPrecond_Beta());
		}

	}

}

void CEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
		CConfig *config, unsigned short iMesh) {
  
	unsigned short iVar;
	unsigned long iPoint;
	bool rotating_frame = config->GetRotating_Frame();
	bool axisymmetric = config->GetAxisymmetric();
	bool incompressible = config->GetIncompressible();
	bool gravity = (config->GetGravityForce() == YES);
	bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
	bool magnet = (config->GetMagnetic_Force() == YES);
	bool jouleheating = config->GetJouleHeating();

	for (iVar = 0; iVar < nVar; iVar++)
		Residual[iVar] = 0;
  
	if (rotating_frame) {

		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) { 

			/*--- Set solution  ---*/
			numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

			/*--- Set control volume ---*/
			numerics->SetVolume(geometry->node[iPoint]->GetVolume());

			/*--- Set rotational velocity ---*/
			numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());

			/*--- Compute Residual ---*/
			numerics->ComputeResidual(Residual, Jacobian_i, config);

			/*--- Add Residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

		}
	}

	if (axisymmetric) {

		bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
		/*--- Zero out Jacobian structure ---*/
		if (implicit) {
			for (iVar = 0; iVar < nVar; iVar ++)
				for (unsigned short jVar = 0; jVar < nVar; jVar ++)
					Jacobian_i[iVar][jVar] = 0.0;
		}

		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) { 

			/*--- Set solution  ---*/
			numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

			if (incompressible) {
				/*--- Set incompressible density  ---*/
				numerics->SetDensityInc(node[iPoint]->GetDensityInc(), node[iPoint]->GetDensityInc());

				/*--- Set beta squared  ---*/
				numerics->SetBetaInc2(node[iPoint]->GetBetaInc2(), node[iPoint]->GetBetaInc2());
			}

			/*--- Set control volume ---*/
			numerics->SetVolume(geometry->node[iPoint]->GetVolume());

			/*--- Set y coordinate ---*/
			numerics->SetCoord(geometry->node[iPoint]->GetCoord(),geometry->node[iPoint]->GetCoord());

			/*--- Compute Source term Residual ---*/
			numerics->ComputeResidual(Residual, Jacobian_i, config);

			/*--- Add Residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- Implicit part ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}

	if (gravity) {

		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) { 

			/*--- Set solution  ---*/
			numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
      
			/*--- Set incompressible density  ---*/
      if (incompressible) numerics->SetDensityInc(node[iPoint]->GetDensityInc(), node[iPoint]->GetDensityInc());

			/*--- Set control volume ---*/
			numerics->SetVolume(geometry->node[iPoint]->GetVolume());

			/*--- Compute Source term Residual ---*/
			numerics->ComputeResidual(Residual, config);
      
			/*--- Add Residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

		}
    
	}

	if (time_spectral) {

		double Volume, Source;

		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

			/*--- Get control volume ---*/
			Volume = geometry->node[iPoint]->GetVolume();

			/*--- Get stored time spectral source term ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Source = node[iPoint]->GetTimeSpectral_Source(iVar);
				Residual[iVar] = Source*Volume;
			}

			/*--- Add Residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

		}
	}

	if (magnet) {
		bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
		for (iVar = 0; iVar < nVar; iVar ++)
			for (unsigned short jVar = 0; jVar < nVar; jVar ++)
				Jacobian_i[iVar][jVar] = 0.0;
		/*--- loop over points ---*/
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

			/*--- Set the coordinate ---*/
			numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());

			/*--- Set solution  ---*/
			numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

			/*--- Set control volume ---*/
			numerics->SetVolume(geometry->node[iPoint]->GetVolume());

			/*--- Compute Residual ---*/
			numerics->ComputeResidual(Residual, Jacobian_i, config);

			LinSysRes.SubtractBlock(iPoint, Residual);
			if (nDim ==3) node[iPoint]->SetMagneticField(numerics->GetMagneticField());
			if (implicit)
				Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

		}
	}

	if (jouleheating) {
		double arc_column_extent = 2.30581;
		bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
		for (iVar = 0; iVar < nVar; iVar ++)
			for (unsigned short jVar = 0; jVar < nVar; jVar ++)
				Jacobian_i[iVar][jVar] = 0.0;

		double integral = 0.0; unsigned long jPoint;
		vector<double> XCoordList = geometry->GetGeometryPlanes();
		vector<vector<double> > Xcoord_plane = geometry->GetXCoord();
		vector<vector<double> > Ycoord_plane = geometry->GetYCoord();
		vector<vector<unsigned long> > Plane_points = geometry->GetPlanarPoints();

		for (unsigned short iPlane = 0; iPlane < XCoordList.size(); iPlane++) {
			integral= 0.0;
			if (Xcoord_plane[iPlane][0] <= arc_column_extent) {
				for (unsigned short iVertex = 1; iVertex < Xcoord_plane[iPlane].size(); iVertex++) {
					iPoint = Plane_points[iPlane][iVertex];
					jPoint = Plane_points[iPlane][iVertex-1];

					/*--- Set solution  ---*/
					numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

					/*--- Set control volume ---*/
					numerics->SetVolume(geometry->node[iPoint]->GetVolume());

					/*--- Set the coordinate ---*/
					numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());

					/*--- Set electrical conductivity  ---*/
					numerics->SetElec_Cond();

					/*--- Sum over points to get the integral  ---*/
					integral +=numerics->GetElec_CondIntegral();
				}

				/*--- Set the integral  ---*/
				numerics->SetElec_CondIntegralsqr(integral*integral);

				for (unsigned short iVertex = 0; iVertex < Xcoord_plane[iPlane].size(); iVertex++) {
					iPoint = Plane_points[iPlane][iVertex];

					/*--- Set solution  ---*/
					numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

					/*--- Set control volume ---*/
					numerics->SetVolume(geometry->node[iPoint]->GetVolume());

					/*--- Set the coordinate ---*/
					numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());

					/*--- Set electrical conductivity  ---*/
					numerics->SetElec_Cond();

					/*--- Compute Residual ---*/
					numerics->ComputeResidual(Residual, Jacobian_i, config);
					//					for (iVar = 0; iVar < nVar; iVar++)
					//						cout << Residual[iVar] << ", ";
					//					cout << endl;
					LinSysRes.SubtractBlock(iPoint, Residual);
					if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
				}
			}
		}
	}
}

void CEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
		CConfig *config, unsigned short iMesh) {

	/* This method should be used to call any new source terms for a particular problem*/
	/* This method calls the new child class in CNumerics, where the new source term should be implemented.  */

	/* Next we describe how to get access to some important quanties for this method */
	/* Access to all points in the current geometric mesh by saying: nPointDomain */
	/* Get the vector of conservative variables at some point iPoint = node[iPoint]->GetSolution() */
	/* Get the volume (or area in 2D) associated with iPoint = node[iPoint]->GetVolume() */
	/* Get the vector of geometric coordinates of point iPoint = node[iPoint]->GetCoord() */

}

void CEulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {
	double *Normal, Area, Mean_SoundSpeed, Mean_ProjVel, Mean_BetaInc2, Lambda, Mean_DensityInc,
	ProjVel, ProjVel_i, ProjVel_j;
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker;

	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();
	bool grid_movement = config->GetGrid_Movement();

	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		node[iPoint]->SetLambda(0.0);
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
			Mean_DensityInc = 0.5 * (node[iPoint]->GetDensityInc() + node[jPoint]->GetDensityInc());
			Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);
		}
		else {
			Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
			Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
		}

		/*--- Adjustment for a rotating frame ---*/
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

		/*--- Adjustment for grid movement ---*/
		if (grid_movement) {
			double *GridVel_i = geometry->node[iPoint]->GetGridVel();
			double *GridVel_j = geometry->node[jPoint]->GetGridVel();
			ProjVel_i = 0.0; ProjVel_j =0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				ProjVel_i += GridVel_i[iDim]*Normal[iDim];
				ProjVel_j += GridVel_j[iDim]*Normal[iDim];
			}
			Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
		}

		/*--- Inviscid contribution ---*/
		Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);

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
				Mean_DensityInc = node[iPoint]->GetDensityInc();
				Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);
			}
			else {
				Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
				Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
			}

			/*--- Adjustment for a rotating frame ---*/
			if (rotating_frame) {
				double *RotVel = geometry->node[iPoint]->GetRotVel();
				ProjVel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					ProjVel += RotVel[iDim]*Normal[iDim];
				Mean_ProjVel -= ProjVel;
			}

			/*--- Adjustment for grid movement ---*/
			if (grid_movement) {
				double *GridVel = geometry->node[iPoint]->GetGridVel();
				ProjVel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					ProjVel += GridVel[iDim]*Normal[iDim];
				Mean_ProjVel -= ProjVel;
			}

			/*--- Inviscid contribution ---*/
			Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
			if (geometry->node[iPoint]->GetDomain()) {
				node[iPoint]->AddLambda(Lambda);
			}

		}
	}

	/*--- MPI parallelization ---*/
	Set_MPI_MaxEigenvalue(geometry, config);

}

void CEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint, iEdge;
	double Pressure_i = 0, Pressure_j = 0, *Diff;
	unsigned short iVar;

	bool incompressible = config->GetIncompressible();

	Diff = new double[nVar];

	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		node[iPoint]->SetUnd_LaplZero();

	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		/*--- Solution differences ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);

		/*--- Correction for compressible flows which use the enthalpy ---*/
		if (!incompressible) {
			Pressure_i = node[iPoint]->GetPressure(incompressible);
			Pressure_j = node[jPoint]->GetPressure(incompressible);
			Diff[nVar-1] = (node[iPoint]->GetSolution(nVar-1) + Pressure_i) - (node[jPoint]->GetSolution(nVar-1) + Pressure_j);
		}

#ifdef STRUCTURED_GRID

		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);

#else

		bool boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
		bool boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();

		/*--- Both points inside the domain, or both in the boundary ---*/
		if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
		}

		/*--- iPoint inside the domain, jPoint on the boundary ---*/
		if (!boundary_i && boundary_j)
			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);

		/*--- jPoint inside the domain, iPoint on the boundary ---*/
		if (boundary_i && !boundary_j)
			if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);

#endif

	}

#ifdef STRUCTURED_GRID

	unsigned long Point_Normal = 0, iVertex;
	double Pressure_mirror = 0, *U_mirror;
	unsigned short iMarker;

	U_mirror = new double[nVar];

	/*--- Loop over all boundaries and include an extra contribution
        from a mirror node. Find the nearest normal, interior point 
        for a boundary node and make a linear approximation. ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {

		if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
				config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
				config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY &&
				config->GetMarker_All_Boundary(iMarker) != PERIODIC_BOUNDARY) {

			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

				if (geometry->node[iPoint]->GetDomain()) {

					Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

					/*--- Interpolate & compute difference in the conserved variables ---*/
					for (iVar = 0; iVar < nVar; iVar++) {
						U_mirror[iVar] = 2.0*node[iPoint]->GetSolution(iVar) - node[Point_Normal]->GetSolution(iVar);
						Diff[iVar]   = node[iPoint]->GetSolution(iVar) - U_mirror[iVar];
					}

					/*--- Correction for compressible flows ---*/
					if (!incompressible) {
						Pressure_mirror = 2.0*node[iPoint]->GetPressure(incompressible) - node[Point_Normal]->GetPressure(incompressible);
						Diff[nVar-1] = (node[iPoint]->GetSolution(nVar-1) + node[iPoint]->GetPressure(incompressible)) - (U_mirror[nVar-1] + Pressure_mirror);
					}

					/*--- Subtract contribution at the boundary node only ---*/
					node[iPoint]->SubtractUnd_Lapl(Diff);
				}
			}
		}
	}

	delete [] U_mirror;

#endif

	delete [] Diff;

	/*--- MPI parallelization ---*/
	Set_MPI_Undivided_Laplacian(geometry, config);

}

void CEulerSolver::SetDissipation_Switch(CGeometry *geometry, CConfig *config) {
	unsigned long iEdge, iPoint, jPoint;
	double Pressure_i, Pressure_j;

	bool incompressible = config->GetIncompressible();

	/*--- Reset variables to store the undivided pressure ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		p1_Und_Lapl[iPoint] = 0.0;
		p2_Und_Lapl[iPoint] = 0.0;
	}

	/*--- Evaluate the pressure sensor ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		/*--- Get the pressure, or density for incompressible solvers ---*/
		if (!incompressible) {
			Pressure_i = node[iPoint]->GetPressure(incompressible);
			Pressure_j = node[jPoint]->GetPressure(incompressible);
		}
		else {
			Pressure_i = node[iPoint]->GetDensityInc();
			Pressure_j = node[jPoint]->GetDensityInc();
		}

#ifdef STRUCTURED_GRID

		/*--- Compute numerator and denominator ---*/
		if (geometry->node[iPoint]->GetDomain()) {
			p1_Und_Lapl[iPoint] += (Pressure_j - Pressure_i);
			p2_Und_Lapl[iPoint] += (Pressure_i + Pressure_j);
		}
		if (geometry->node[jPoint]->GetDomain()) {
			p1_Und_Lapl[jPoint] += (Pressure_i - Pressure_j);
			p2_Und_Lapl[jPoint] += (Pressure_i + Pressure_j);
		}

#else

		bool boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
		bool boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();

		/*--- Both points inside the domain, or both on the boundary ---*/
		if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)){
			if (geometry->node[iPoint]->GetDomain()) { p1_Und_Lapl[iPoint] += (Pressure_j - Pressure_i); p2_Und_Lapl[iPoint] += (Pressure_i + Pressure_j); }
			if (geometry->node[jPoint]->GetDomain()) { p1_Und_Lapl[jPoint] += (Pressure_i - Pressure_j); p2_Und_Lapl[jPoint] += (Pressure_i + Pressure_j); }
		}

		/*--- iPoint inside the domain, jPoint on the boundary ---*/
		if (!boundary_i && boundary_j)
			if (geometry->node[iPoint]->GetDomain()) { p1_Und_Lapl[iPoint] += (Pressure_j - Pressure_i); p2_Und_Lapl[iPoint] += (Pressure_i + Pressure_j); }

		/*--- jPoint inside the domain, iPoint on the boundary ---*/
		if (boundary_i && !boundary_j)
			if (geometry->node[jPoint]->GetDomain()) { p1_Und_Lapl[jPoint] += (Pressure_i - Pressure_j); p2_Und_Lapl[jPoint] += (Pressure_i + Pressure_j); }

#endif

	}

#ifdef STRUCTURED_GRID
	unsigned short iMarker;
	unsigned long iVertex, Point_Normal;
	double Press_mirror;

	/*--- Loop over all boundaries and include an extra contribution
   from a mirror node. Find the nearest normal, interior point
   for a boundary node and make a linear approximation. ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {

		if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
				config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
				config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY &&
				config->GetMarker_All_Boundary(iMarker) != PERIODIC_BOUNDARY) {

			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

				if (geometry->node[iPoint]->GetDomain()) {

					Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

					if (!incompressible) Pressure_i = node[iPoint]->GetPressure(incompressible);
					else Pressure_i = node[iPoint]->GetDensityInc();

					/*--- Interpolate & compute difference in the conserved variables ---*/
					if (!incompressible) {
						Pressure_i = node[iPoint]->GetPressure(incompressible);
						Press_mirror = 2.0*Pressure_i - node[Point_Normal]->GetPressure(incompressible);
					}
					else {
						Pressure_i = node[iPoint]->GetDensityInc();
						Press_mirror = 2.0*Pressure_i - node[Point_Normal]->GetDensityInc();
					}

					/*--- Add contribution at the boundary node only ---*/
					p1_Und_Lapl[iPoint] += (Press_mirror - Pressure_i);
					p2_Und_Lapl[iPoint] += (Pressure_i + Press_mirror);

				}
			}
		}
	}

#endif

	/*--- Set pressure switch for each point ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		node[iPoint]->SetSensor(fabs(p1_Und_Lapl[iPoint]) / p2_Und_Lapl[iPoint]);

	/*--- MPI parallelization ---*/
	Set_MPI_Dissipation_Switch(geometry, config);

}

void CEulerSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) {
	unsigned long iVertex, iPoint;
	unsigned short iDim, iMarker, Boundary, Monitoring;
	double Pressure, *Normal = NULL, dist[3], *Coord, Face_Area, PressInviscid;
	double factor, NFPressOF, RefVel2, RefDensity, RefPressure;

	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();

	double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
	double Beta            = config->GetAoS()*PI_NUMBER/180.0;
	double RefAreaCoeff    = config->GetRefAreaCoeff();
	double RefLengthMoment = config->GetRefLengthMoment();
	double *Origin         = config->GetRefOriginMoment();
	//double Pressure_Inf    = config->GetPressure_FreeStreamND();
	//double *Velocity_Inf   = config->GetVelocity_FreeStreamND();
	bool incompressible    = config->GetIncompressible();

	/*--- If we have a rotating frame problem or an unsteady problem with
   mesh motion, use special reference values for the force coefficients.
   Otherwise, use the freestream values, which is the standard convention. ---*/

	if (rotating_frame) {
    
    /*--- Use the rotational origin for computing moments ---*/
    Origin = config->GetRotAxisOrigin();
    
    /*--- Reference moment length is the rotor radius. Note that
     this assumes that the reference area was set to the rotor disk area. ---*/
    RefLengthMoment = sqrt(RefAreaCoeff/PI_NUMBER);
    
    /*--- Reference velocity is the rotational speed times rotor radius ---*/
    RefVel2 = (config->GetOmegaMag()*RefLengthMoment)*(config->GetOmegaMag()*RefLengthMoment);
  
  } else if (grid_movement) {
		double Gas_Constant = config->GetGas_ConstantND();
		double Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
		double Mach_Motion = config->GetMach_Motion();
		RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
	} else {
		RefVel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
	}

	RefDensity  = Density_Inf;
	RefPressure = Pressure_Inf;

	/*-- Initialization ---*/
	Total_CDrag = 0.0; Total_CLift = 0.0;  Total_CSideForce = 0.0;
	Total_CMx = 0.0;   Total_CMy = 0.0;    Total_CMz = 0.0;
	Total_CFx = 0.0;   Total_CFy = 0.0;    Total_CFz = 0.0;
	Total_CEff = 0.0;  Total_CMerit = 0.0; Total_CNearFieldOF = 0.0;
	Total_CT = 0.0;    Total_CQ = 0.0;     Total_Q = 0.0;
  Total_Maxq = 0.0;
	AllBound_CDrag_Inv = 0.0;        AllBound_CLift_Inv = 0.0;  AllBound_CSideForce_Inv = 0.0;
	AllBound_CMx_Inv = 0.0;          AllBound_CMy_Inv = 0.0;    AllBound_CMz_Inv = 0.0;
	AllBound_CFx_Inv = 0.0;          AllBound_CFy_Inv = 0.0;    AllBound_CFz_Inv = 0.0;
	AllBound_CEff_Inv = 0.0;         AllBound_CMerit_Inv = 0.0; AllBound_CQ_Inv = 0.0;
	AllBound_CNearFieldOF_Inv = 0.0; AllBound_CT_Inv = 0.0;

	/*--- Loop over the Euler and Navier-Stokes markers ---*/
	factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);

	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		Boundary   = config->GetMarker_All_Boundary(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);

		if ((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) ||
				(Boundary == ISOTHERMAL) || (Boundary == NEARFIELD_BOUNDARY)) {

			for (iDim = 0; iDim < nDim; iDim++) ForceInviscid[iDim] = 0.0;
			MomentInviscid[0] = 0.0; MomentInviscid[1] = 0.0; MomentInviscid[2] = 0.0;
			NFPressOF = 0.0; PressInviscid = 0.0;

			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Pressure = node[iPoint]->GetPressure(incompressible);

				CPressure[iMarker][iVertex] = (Pressure - RefPressure)*factor*RefAreaCoeff;

				/*--- Note that the pressure coefficient is computed at the
				 halo cells (for visualization purposes), but not the forces ---*/
				if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					Coord = geometry->node[iPoint]->GetCoord();

					/*--- Quadratic objective function for the near field.
           This uses the infinity pressure regardless of Mach number. ---*/
					NFPressOF += 0.5*(Pressure - Pressure_Inf)*(Pressure - Pressure_Inf)*Normal[nDim-1];

					Face_Area = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						/*--- Total force, distance computation, and face area 
             Note that we have subtracted the Pressure at the infinity, this is important when dealing with
             non-closed surfaces---*/
						ForceInviscid[iDim] -= (Pressure - Pressure_Inf)*Normal[iDim]*factor;
						dist[iDim] = Coord[iDim] - Origin[iDim];
						Face_Area += Normal[iDim]*Normal[iDim];
					}
					Face_Area = sqrt(Face_Area);
					PressInviscid += CPressure[iMarker][iVertex]*Face_Area;

					/*--- Moment with respect to the reference axis ---*/
					if (iDim == 3) {
						MomentInviscid[0] -= (Pressure - Pressure_Inf)*(Normal[2]*dist[1]-Normal[1]*dist[2])*factor/RefLengthMoment;
						MomentInviscid[1] -= (Pressure - Pressure_Inf)*(Normal[0]*dist[2]-Normal[2]*dist[0])*factor/RefLengthMoment;
					}
					MomentInviscid[2]   -= (Pressure - Pressure_Inf)*(Normal[1]*dist[0]-Normal[0]*dist[1])*factor/RefLengthMoment;

				}
			}

			/*--- Transform ForceInviscid and MomentInviscid into non-dimensional coefficient ---*/
			if  (Monitoring == YES) {
				if (nDim == 2) {
					if (Boundary != NEARFIELD_BOUNDARY) {
						CDrag_Inv[iMarker]        =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
						CLift_Inv[iMarker]        = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
						CSideForce_Inv[iMarker]   = 0.0;
						CMx_Inv[iMarker]          = 0.0;
						CMy_Inv[iMarker]          = 0.0;
						CMz_Inv[iMarker]          = MomentInviscid[2];
						CEff_Inv[iMarker]         = CLift_Inv[iMarker]/(CDrag_Inv[iMarker]+config->GetCteViscDrag()+EPS);
						CNearFieldOF_Inv[iMarker] = 0.0;
						CFx_Inv[iMarker]          = ForceInviscid[0];
						CFy_Inv[iMarker]          = ForceInviscid[1];
						CFz_Inv[iMarker]          = 0.0;
						CT_Inv[iMarker]           = -CFx_Inv[iMarker];
						CQ_Inv[iMarker]           = -CMz_Inv[iMarker];
						CMerit_Inv[iMarker]       = 0.0;
					}
					else {
						CDrag_Inv[iMarker] = 0.0; CLift_Inv[iMarker] = 0.0; CSideForce_Inv[iMarker] = 0.0;
						CMx_Inv[iMarker] = 0.0; CMy_Inv[iMarker] = 0.0; CMz_Inv[iMarker] = 0.0;
						CFx_Inv[iMarker] = 0.0; CFy_Inv[iMarker] = 0.0; CFz_Inv[iMarker] = 0.0;
						CEff_Inv[iMarker] = 0.0;
						CNearFieldOF_Inv[iMarker] = NFPressOF;
						CT_Inv[iMarker] = 0.0;
						CQ_Inv[iMarker] = 0.0;
						CMerit_Inv[iMarker] = 0.0;
					}
				}
				if (nDim == 3) {
					if (Boundary != NEARFIELD_BOUNDARY) {
						CDrag_Inv[iMarker] =  ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta);
						CLift_Inv[iMarker] = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
						CSideForce_Inv[iMarker] = -ForceInviscid[0]*sin(Beta)*cos(Alpha) + ForceInviscid[1]*cos(Beta) - ForceInviscid[2]*sin(Beta)*sin(Alpha);
						CMx_Inv[iMarker] = MomentInviscid[0];
						CMy_Inv[iMarker] = MomentInviscid[1];
						CMz_Inv[iMarker] = MomentInviscid[2];
						CEff_Inv[iMarker] = CLift_Inv[iMarker]/(CDrag_Inv[iMarker]+config->GetCteViscDrag()+EPS);
						CNearFieldOF_Inv[iMarker] = 0.0;
						CFx_Inv[iMarker] = ForceInviscid[0];
						CFy_Inv[iMarker] = ForceInviscid[1];
						CFz_Inv[iMarker] = ForceInviscid[2];
						CT_Inv[iMarker] = -CFz_Inv[iMarker];
						CQ_Inv[iMarker] = -CMz_Inv[iMarker];
						/*--- For now, define the FM as a simple rotor efficiency ---*/
						//CMerit_Inv[iMarker] = CT_Inv[iMarker]*sqrt(fabs(CT_Inv[iMarker]))/(sqrt(2.0)*CQ_Inv[iMarker]);
						CMerit_Inv[iMarker] = CT_Inv[iMarker]/CQ_Inv[iMarker];
					}
					else {
						CDrag_Inv[iMarker] = 0.0; CLift_Inv[iMarker] = 0.0; CSideForce_Inv[iMarker] = 0.0;
						CMx_Inv[iMarker] = 0.0; CMy_Inv[iMarker] = 0.0; CMz_Inv[iMarker] = 0.0;
						CFx_Inv[iMarker] = 0.0; CFy_Inv[iMarker] = 0.0; CFz_Inv[iMarker] = 0.0;
						CEff_Inv[iMarker] = 0.0;
						CNearFieldOF_Inv[iMarker] = NFPressOF;
						CT_Inv[iMarker] = 0.0;
						CQ_Inv[iMarker] = 0.0;
						CMerit_Inv[iMarker] = 0.0;
					}
				}

				AllBound_CDrag_Inv += CDrag_Inv[iMarker];
				AllBound_CLift_Inv += CLift_Inv[iMarker];
				AllBound_CSideForce_Inv += CSideForce_Inv[iMarker];
				AllBound_CMx_Inv += CMx_Inv[iMarker];
				AllBound_CMy_Inv += CMy_Inv[iMarker];
				AllBound_CMz_Inv += CMz_Inv[iMarker];
				AllBound_CEff_Inv += CEff_Inv[iMarker];
				AllBound_CNearFieldOF_Inv += CNearFieldOF_Inv[iMarker];
				AllBound_CFx_Inv += CFx_Inv[iMarker];
				AllBound_CFy_Inv += CFy_Inv[iMarker];
				AllBound_CFz_Inv += CFz_Inv[iMarker];
				AllBound_CT_Inv += CT_Inv[iMarker];
				AllBound_CQ_Inv += CQ_Inv[iMarker];
				AllBound_CMerit_Inv += CMerit_Inv[iMarker];
			}
		}
	}

	Total_CDrag += AllBound_CDrag_Inv;
	Total_CLift += AllBound_CLift_Inv;
	Total_CSideForce += AllBound_CSideForce_Inv;
	Total_CMx += AllBound_CMx_Inv;
	Total_CMy += AllBound_CMy_Inv;
	Total_CMz += AllBound_CMz_Inv;
	Total_CEff += AllBound_CEff_Inv;
	Total_CNearFieldOF += AllBound_CNearFieldOF_Inv;
	Total_CFx += AllBound_CFx_Inv;
	Total_CFy += AllBound_CFy_Inv;
	Total_CFz += AllBound_CFz_Inv;
	Total_CT += AllBound_CT_Inv;
	Total_CQ += AllBound_CQ_Inv;
	Total_CMerit += AllBound_CMerit_Inv;
}

void CEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
		CConfig *config, unsigned short iRKStep) {
	double *Residual, *Res_TruncError, Vol, Delta, Res;
	unsigned short iVar;
	unsigned long iPoint;

	double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
	bool adjoint = config->GetAdjoint();

	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
		SetRes_Max(iVar, 0.0, 0);
	}

	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Delta = node[iPoint]->GetDelta_Time() / Vol;

		Res_TruncError = node[iPoint]->GetResTruncError();
		Residual = LinSysRes.GetBlock(iPoint);

		if (!adjoint) {
			for (iVar = 0; iVar < nVar; iVar++) {
				Res = Residual[iVar] + Res_TruncError[iVar];
				node[iPoint]->AddSolution(iVar, -Res*Delta*RK_AlphaCoeff);
				AddRes_RMS(iVar, Res*Res);
				AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex());
			}
		}

	}

	/*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);

	/*--- Compute the root mean square residual ---*/
	SetResidual_RMS(geometry, config);
    
    
}

void CEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
	unsigned short iVar;
	unsigned long iPoint;

	bool adjoint = config->GetAdjoint();

	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
		SetRes_Max(iVar, 0.0, 0);
	}

	/*--- Update the solution ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Delta = node[iPoint]->GetDelta_Time() / Vol;

		local_Res_TruncError = node[iPoint]->GetResTruncError();
		local_Residual = LinSysRes.GetBlock(iPoint);

		if (!adjoint) {
			for (iVar = 0; iVar < nVar; iVar++) {
				Res = local_Residual[iVar] + local_Res_TruncError[iVar];
				node[iPoint]->AddSolution(iVar, -Res*Delta);
				AddRes_RMS(iVar, Res*Res);
				AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex());
			}
		}

	}

	/*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);

	/*--- Compute the root mean square residual ---*/
	SetResidual_RMS(geometry, config);

}

void CEulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	unsigned short iVar, jVar;
	unsigned long iPoint, total_index, IterLinSol = 0;
	double Delta, *local_Res_TruncError, Vol;
    
	bool adjoint = config->GetAdjoint();
    
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
		SetRes_Max(iVar, 0.0, 0);
	}
    
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
        
		/*--- Read the residual ---*/
		local_Res_TruncError = node[iPoint]->GetResTruncError();
        
		/*--- Read the volume ---*/
		Vol = geometry->node[iPoint]->GetVolume();
        
		/*--- Modify matrix diagonal to assure diagonal dominance ---*/
		Delta = Vol / node[iPoint]->GetDelta_Time();
        
		if (roe_turkel) {
			SetPreconditioner(config, iPoint);
			for (iVar = 0; iVar < nVar; iVar ++ )
				for (jVar = 0; jVar < nVar; jVar ++ )
					Precon_Mat_inv[iVar][jVar] = Delta*Precon_Mat_inv[iVar][jVar];
			Jacobian.AddBlock(iPoint, iPoint, Precon_Mat_inv);
		}
		else {
			Jacobian.AddVal2Diag(iPoint, Delta);
		}
        
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar + iVar;
			LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
			AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
			AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex());
		}
	}
    
	/*--- Initialize residual and solution at the ghost points ---*/
	for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar + iVar;
			LinSysRes[total_index] = 0.0;
			LinSysSol[total_index] = 0.0;
		}
	}
    
	/*--- Solve the linear system (Krylov subspace methods) ---*/
  CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
  
  CPreconditioner* precond = NULL;
  if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
    Jacobian.BuildJacobiPreconditioner();
    precond = new CJacobiPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LU_SGS) {
    precond = new CLU_SGSPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    Jacobian.BuildJacobiPreconditioner();
    Jacobian.BuildLineletPreconditioner(geometry, config);
    precond = new CLineletPreconditioner(Jacobian, geometry, config);
  }
  
  CSysSolve system;
  if (config->GetKind_Linear_Solver() == BCGSTAB)
    IterLinSol = system.BCGSTAB(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                                config->GetLinear_Solver_Iter(), false);
  else if (config->GetKind_Linear_Solver() == FGMRES)
    IterLinSol = system.FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                              config->GetLinear_Solver_Iter(), false);
  
  /*--- The the number of iterations of the linear solver ---*/
  SetIterLinSolver(IterLinSol);
  
  /*--- dealocate memory ---*/
  delete mat_vec;
  delete precond;
  
	/*--- Update solution (system written in terms of increments) ---*/
	if (!adjoint) {
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar]);
			}
		}
	}
    
	/*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);

	/*--- Compute the root mean square residual ---*/
	SetResidual_RMS(geometry, config);
    
}

void CEulerSolver::SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint, iEdge, iVertex;
	unsigned short iDim, iVar, iMarker, nPrimVar;
	double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
	Partial_Gradient, Partial_Res, *Normal;

	bool incompressible = config->GetIncompressible();

	if (incompressible) nPrimVar = nDim+2;
	else nPrimVar = nDim+3;

	/*--- Gradient primitive variables compressible (temp, vx, vy, vz, P, rho)
   Gradient primitive variables incompressible (rho, vx, vy, vz, beta) ---*/
	PrimVar_Vertex = new double [nPrimVar];
	PrimVar_i = new double [nPrimVar];
	PrimVar_j = new double [nPrimVar];

	/*--- Set Gradient_Primitive to zero ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		node[iPoint]->SetGradient_PrimitiveZero(nPrimVar);

	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		for (iVar = 0; iVar < nPrimVar; iVar++) {
			PrimVar_i[iVar] = node[iPoint]->GetPrimVar(iVar);
			PrimVar_j[iVar] = node[jPoint]->GetPrimVar(iVar);
		}

		Normal = geometry->edge[iEdge]->GetNormal();
		for (iVar = 0; iVar < nPrimVar; iVar++) {
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

				for (iVar = 0; iVar < nPrimVar; iVar++)
					PrimVar_Vertex[iVar] = node[iPoint]->GetPrimVar(iVar);

				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				for (iVar = 0; iVar < nPrimVar; iVar++)
					for (iDim = 0; iDim < nDim; iDim++) {
						Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
						node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
					}
			}
		}
	}

	/*--- Update gradient value ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		for (iVar = 0; iVar < nPrimVar; iVar++) {
			for (iDim = 0; iDim < nDim; iDim++) {
				Partial_Gradient = node[iPoint]->GetGradient_Primitive(iVar,iDim) / (geometry->node[iPoint]->GetVolume());
				node[iPoint]->SetGradient_Primitive(iVar, iDim, Partial_Gradient);
			}
		}
	}

	delete [] PrimVar_Vertex;
	delete [] PrimVar_i;
	delete [] PrimVar_j;

	Set_MPI_PrimVar_Gradient(geometry, config);

}

void CEulerSolver::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, jDim, iNeigh, nPrimVar;
	unsigned long iPoint, jPoint;
	double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
	r23_b, r33, weight, product;

	bool incompressible = config->GetIncompressible();

	if (incompressible) nPrimVar = nDim+2;
	else nPrimVar = nDim+3;

	/*--- Gradient primitive variables compressible (temp, vx, vy, vz, P, rho)
   Gradient primitive variables incompressible (rho, vx, vy, vz, beta) ---*/
	PrimVar_i = new double [nPrimVar];
	PrimVar_j = new double [nPrimVar];

	/*--- Loop over points of the grid ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		Coord_i = geometry->node[iPoint]->GetCoord();

		for (iVar = 0; iVar < nPrimVar; iVar++)
			PrimVar_i[iVar] = node[iPoint]->GetPrimVar(iVar);

		/*--- Inizialization of variables ---*/
		for (iVar = 0; iVar < nPrimVar; iVar++)
			for (iDim = 0; iDim < nDim; iDim++)
				cvector[iVar][iDim] = 0.0;

		r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0;
		r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

		for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
			Coord_j = geometry->node[jPoint]->GetCoord();

			for (iVar = 0; iVar < nPrimVar; iVar++)
				PrimVar_j[iVar] = node[jPoint]->GetPrimVar(iVar);

			weight = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

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
			for (iVar = 0; iVar < nPrimVar; iVar++)
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
		for (iVar = 0; iVar < nPrimVar; iVar++) {
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

	Set_MPI_PrimVar_Gradient(geometry, config);

}

void CEulerSolver::SetPrimVar_Limiter(CGeometry *geometry, CConfig *config) { }

void CEulerSolver::SetPreconditioner(CConfig *config, unsigned short iPoint) {
	unsigned short iDim, jDim, iVar, jVar;
	double Beta, local_Mach, Beta2, rho, enthalpy, soundspeed, sq_vel;
	double *U_i = NULL;
	double Beta_min = config->GetminTurkelBeta();
	double Beta_max = config->GetmaxTurkelBeta();


	/*--- Variables to calculate the preconditioner parameter Beta ---*/
	local_Mach = sqrt(node[iPoint]->GetVelocity2())/node[iPoint]->GetSoundSpeed();
	Beta 		    = max(Beta_min,min(local_Mach,Beta_max));
	Beta2 		    = Beta*Beta;

	U_i = node[iPoint]->GetSolution();

	rho = U_i[0];
	enthalpy = node[iPoint]->GetEnthalpy();
	soundspeed = node[iPoint]->GetSoundSpeed();
	sq_vel = node[iPoint]->GetVelocity2();

	/*---Calculating the inverse of the preconditioning matrix that multiplies the time derivative  */
	Precon_Mat_inv[0][0] = 0.5*sq_vel;
	Precon_Mat_inv[0][nVar-1] = 1.0;
	for (iDim = 0; iDim < nDim; iDim ++)
		Precon_Mat_inv[0][1+iDim] = -1.0*U_i[iDim+1]/rho;

	for (iDim = 0; iDim < nDim; iDim ++) {
		Precon_Mat_inv[iDim+1][0] = 0.5*sq_vel*U_i[iDim+1]/rho;
		Precon_Mat_inv[iDim+1][nVar-1] = U_i[iDim+1]/rho;
		for (jDim = 0; jDim < nDim; jDim ++) {
			Precon_Mat_inv[iDim+1][1+jDim] = -1.0*U_i[jDim+1]/rho*U_i[iDim+1]/rho;
		}
	}

	Precon_Mat_inv[nVar-1][0] = 0.5*sq_vel*enthalpy;
	Precon_Mat_inv[nVar-1][nVar-1] = enthalpy;
	for (iDim = 0; iDim < nDim; iDim ++)
		Precon_Mat_inv[nVar-1][1+iDim] = -1.0*U_i[iDim+1]/rho*enthalpy;


	for (iVar = 0; iVar < nVar; iVar ++ ) {
		for (jVar = 0; jVar < nVar; jVar ++ ) {
			Precon_Mat_inv[iVar][jVar] = (1.0/(Beta2+EPS) - 1.0) * (Gamma-1.0)/(soundspeed*soundspeed)*Precon_Mat_inv[iVar][jVar];
			if (iVar == jVar)
				Precon_Mat_inv[iVar][iVar] += 1.0;
		}
	}

}

void CEulerSolver::GetNacelle_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh) {
	unsigned short iDim, iMarker, iVar;
	unsigned long iVertex, iPoint;
	double Pressure, Velocity[3], Velocity2, MassFlow, Density, Energy, Area,
  Mach, SoundSpeed, Flow_Dir[3], alpha;
  
  unsigned short nMarker_NacelleInflow = config->GetnMarker_NacelleInflow();
  unsigned short nMarker_NacelleExhaust = config->GetnMarker_NacelleExhaust();
  
  if ((nMarker_NacelleInflow != 0) || (nMarker_NacelleExhaust != 0)) {
    
    int rank = MASTER_NODE;
#ifndef NO_MPI
    rank = MPI::COMM_WORLD.Get_rank();
#endif
    
    /*--- Compute the numerical fan face Mach number, and the total area of the inflow ---*/
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      FanFace_MassFlow[iMarker] = 0.0;
      FanFace_Mach[iMarker] = 0.0;
      FanFace_Pressure[iMarker] = 0.0;
      FanFace_Area[iMarker] = 0.0;
      
      Exhaust_MassFlow[iMarker] = 0.0;
      Exhaust_Area[iMarker] = 0.0;
      
      if (config->GetMarker_All_Boundary(iMarker) == NACELLE_INFLOW) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (geometry->node[iPoint]->GetDomain()) {
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
            
            Density = node[iPoint]->GetSolution(0);
            Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Area += Vector[iDim]*Vector[iDim];
              Velocity[iDim] = node[iPoint]->GetSolution(iDim+1)/Density;
              Velocity2 += Velocity[iDim]*Velocity[iDim];
              MassFlow -= Vector[iDim]*node[iPoint]->GetSolution(iDim+1);
            }
            
            Area       = sqrt (Area);
            Energy     = node[iPoint]->GetSolution(nVar-1)/Density;
            Pressure   = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
            SoundSpeed = sqrt(Gamma*Pressure/Density);
            Mach       = sqrt(Velocity2)/SoundSpeed;
            
            /*--- Compute the FanFace_MassFlow, FanFace_Pressure, FanFace_Mach, and FanFace_Area ---*/
            FanFace_MassFlow[iMarker] += MassFlow;
            FanFace_Pressure[iMarker] += Pressure*Area;
            FanFace_Mach[iMarker] += Mach*Area;
            FanFace_Area[iMarker] += Area;
            
          }
        }
        
      }
      
      if (config->GetMarker_All_Boundary(iMarker) == NACELLE_EXHAUST) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (geometry->node[iPoint]->GetDomain()) {
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
            
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Area += Vector[iDim]*Vector[iDim];
            
            Area = sqrt (Area);
            
            MassFlow = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              MassFlow += Vector[iDim]*node[iPoint]->GetSolution(iDim+1);
            
            /*--- Compute the mass Exhaust_MassFlow ---*/
            Exhaust_MassFlow[iMarker] += MassFlow;
            Exhaust_Area[iMarker] += Area;
            
          }
        }
        
      }
      
    }
    
    /*--- Copy to the appropiate structure ---*/
    unsigned short iMarker_NacelleInflow, iMarker_NacelleExhaust;
    
    double FanFace_MassFlow_Local[nMarker_NacelleInflow];
    double FanFace_Mach_Local[nMarker_NacelleInflow];
    double FanFace_Pressure_Local[nMarker_NacelleInflow];
    double FanFace_Area_Local[nMarker_NacelleInflow];
    
    double FanFace_MassFlow_Total[nMarker_NacelleInflow];
    double FanFace_Mach_Total[nMarker_NacelleInflow];
    double FanFace_Pressure_Total[nMarker_NacelleInflow];
    double FanFace_Area_Total[nMarker_NacelleInflow];
    
    for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
      FanFace_MassFlow_Local[iMarker_NacelleInflow] = 0.0;
      FanFace_Mach_Local[iMarker_NacelleInflow] = 0.0;
      FanFace_Pressure_Local[iMarker_NacelleInflow] = 0.0;
      FanFace_Area_Local[iMarker_NacelleInflow] = 0.0;
      
      FanFace_MassFlow_Total[iMarker_NacelleInflow] = 0.0;
      FanFace_Mach_Total[iMarker_NacelleInflow] = 0.0;
      FanFace_Pressure_Total[iMarker_NacelleInflow] = 0.0;
      FanFace_Area_Total[iMarker_NacelleInflow] = 0.0;
    }
    
    double Exhaust_MassFlow_Local[nMarker_NacelleExhaust];
    double Exhaust_Area_Local[nMarker_NacelleExhaust];
    
    double Exhaust_MassFlow_Total[nMarker_NacelleExhaust];
    double Exhaust_Area_Total[nMarker_NacelleExhaust];
    
    for (iMarker_NacelleExhaust = 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++) {
      Exhaust_MassFlow_Local[iMarker_NacelleExhaust] = 0.0;
      Exhaust_Area_Local[iMarker_NacelleExhaust] = 0.0;
      
      Exhaust_MassFlow_Total[iMarker_NacelleExhaust] = 0.0;
      Exhaust_Area_Total[iMarker_NacelleExhaust] = 0.0;
    }
    
    /*--- Compute the numerical fan face Mach number, and the total area of the inflow ---*/
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_Boundary(iMarker) == NACELLE_INFLOW) {
        
        /*--- Loop over all the boundaries with nacelle inflow bc ---*/
        for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
          
          /*--- Add the FanFace_MassFlow, FanFace_Mach, FanFace_Pressure and FanFace_Area to the particular boundary ---*/
          if (config->GetMarker_All_Tag(iMarker) == config->GetMarker_NacelleInflow(iMarker_NacelleInflow)) {
            FanFace_MassFlow_Local[iMarker_NacelleInflow] += FanFace_MassFlow[iMarker];
            FanFace_Mach_Local[iMarker_NacelleInflow] += FanFace_Mach[iMarker];
            FanFace_Pressure_Local[iMarker_NacelleInflow] += FanFace_Pressure[iMarker];
            FanFace_Area_Local[iMarker_NacelleInflow] += FanFace_Area[iMarker];
          }
          
        }
        
      }
      
      if (config->GetMarker_All_Boundary(iMarker) == NACELLE_EXHAUST) {
        
        /*--- Loop over all the boundaries with nacelle inflow bc ---*/
        for (iMarker_NacelleExhaust= 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++) {
          
          /*--- Add the Exhaust_MassFlow, and Exhaust_Area to the particular boundary ---*/
          if (config->GetMarker_All_Tag(iMarker) == config->GetMarker_NacelleExhaust(iMarker_NacelleExhaust)) {
            Exhaust_MassFlow_Local[iMarker_NacelleExhaust] += Exhaust_MassFlow[iMarker];
            Exhaust_Area_Local[iMarker_NacelleExhaust] += Exhaust_Area[iMarker];
          }
          
        }
        
      }
      
    }
    
#ifndef NO_MPI
    
    MPI::COMM_WORLD.Allreduce(FanFace_MassFlow_Local, FanFace_MassFlow_Total, nMarker_NacelleInflow, MPI::DOUBLE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(FanFace_Mach_Local, FanFace_Mach_Total, nMarker_NacelleInflow, MPI::DOUBLE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(FanFace_Pressure_Local, FanFace_Pressure_Total, nMarker_NacelleInflow, MPI::DOUBLE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(FanFace_Area_Local, FanFace_Area_Total, nMarker_NacelleInflow, MPI::DOUBLE, MPI::SUM);
    
    MPI::COMM_WORLD.Allreduce(Exhaust_MassFlow_Local, Exhaust_MassFlow_Total, nMarker_NacelleExhaust, MPI::DOUBLE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(Exhaust_Area_Local, Exhaust_Area_Total, nMarker_NacelleExhaust, MPI::DOUBLE, MPI::SUM);
    
#else
    
    for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
      FanFace_MassFlow_Total[iMarker_NacelleInflow]   = FanFace_MassFlow_Local[iMarker_NacelleInflow];
      FanFace_Mach_Total[iMarker_NacelleInflow]       = FanFace_Mach_Local[iMarker_NacelleInflow];
      FanFace_Pressure_Total[iMarker_NacelleInflow]   = FanFace_Pressure_Local[iMarker_NacelleInflow];
      FanFace_Area_Total[iMarker_NacelleInflow]       = FanFace_Area_Local[iMarker_NacelleInflow];
    }
    
    for (iMarker_NacelleExhaust = 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++) {
      Exhaust_MassFlow_Total[iMarker_NacelleExhaust]  = Exhaust_MassFlow_Local[iMarker_NacelleExhaust];
      Exhaust_Area_Total[iMarker_NacelleExhaust]      = Exhaust_Area_Local[iMarker_NacelleExhaust];
    }
    
#endif
    
    /*--- Compute the value of FanFace_Area_Total, and FanFace_Pressure_Total, and
     set the value in the config structure for future use ---*/
    for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
      if (FanFace_Area_Total[iMarker_NacelleInflow] != 0.0) FanFace_Mach_Total[iMarker_NacelleInflow] /= FanFace_Area_Total[iMarker_NacelleInflow];
      else FanFace_Mach_Total[iMarker_NacelleInflow] = 0.0;
      if (FanFace_Area_Total[iMarker_NacelleInflow] != 0.0) FanFace_Pressure_Total[iMarker_NacelleInflow] /= FanFace_Area_Total[iMarker_NacelleInflow];
      else FanFace_Pressure_Total[iMarker_NacelleInflow] = 0.0;
      
      if (iMesh == MESH_0) {
        config->SetFanFace_Mach(iMarker_NacelleInflow, FanFace_Mach_Total[iMarker_NacelleInflow]);
        config->SetFanFace_Pressure(iMarker_NacelleInflow, FanFace_Pressure_Total[iMarker_NacelleInflow]);
      }
      
    }
    
    bool write_heads = (((config->GetExtIter() % (config->GetWrt_Con_Freq()*20)) == 0));
    if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads) {
      
      cout.precision(4);
      cout.setf(ios::fixed,ios::floatfield);
      
      cout << endl << "---------------------------- Engine properties --------------------------" << endl;
      for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
        cout << "Nacelle inflow ("<< config->GetMarker_NacelleInflow(iMarker_NacelleInflow)
        << "): MassFlow (kg/s): " << FanFace_MassFlow_Total[iMarker_NacelleInflow] * config->GetDensity_Ref() * config->GetVelocity_Ref()
        << ", Mach: " << FanFace_Mach_Total[iMarker_NacelleInflow]
        << ", Area: " << FanFace_Area_Total[iMarker_NacelleInflow] <<"."<< endl;
      }
      
      for (iMarker_NacelleExhaust = 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++) {
        cout << "Nacelle exhaust ("<< config->GetMarker_NacelleExhaust(iMarker_NacelleExhaust) 
        << "): MassFlow (kg/s): " << Exhaust_MassFlow_Total[iMarker_NacelleExhaust] * config->GetDensity_Ref() * config->GetVelocity_Ref()
        << ", Area: " << Exhaust_Area_Total[iMarker_NacelleExhaust] <<"."<< endl;
      }
      cout << "-------------------------------------------------------------------------" << endl;
      
    }
    
    /*--- Check the flow orientation in the nacelle inflow ---*/
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_Boundary(iMarker) == NACELLE_INFLOW) {
        
        /*--- Loop over all the vertices on this boundary marker ---*/
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          /*--- Normal vector for this vertex (negate for outward convention) ---*/
          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
          
          for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
          
          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Vector[iDim]*Vector[iDim];
          Area = sqrt (Area);
          
          /*--- Compute unitary vector ---*/
          for (iDim = 0; iDim < nDim; iDim++)
            Vector[iDim] /= Area;
          
          /*--- The flow direction is defined by the local velocity on the surface ---*/
          for (iDim = 0; iDim < nDim; iDim++)
            Flow_Dir[iDim] = node[iPoint]->GetSolution(iDim+1) / node[iPoint]->GetSolution(0);
          
          /*--- Dot product of normal and flow direction. ---*/
          alpha = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            alpha += Vector[iDim]*Flow_Dir[iDim];
          
          /*--- Flow in the wrong direction. ---*/
          if (alpha < 0.0) {
            
            /*--- Copy the old solution ---*/
            for (iVar = 0; iVar < nVar; iVar++)
              node[iPoint]->SetSolution(iVar, node[iPoint]->GetSolution_Old(iVar));
            
          }
          
        }
      }
    }
    
  }
  
}

void CEulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
		CNumerics *numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iDim, iVar, jVar, jDim;
	unsigned long iPoint, iVertex;
	double Pressure, *Normal = NULL, *GridVel = NULL, Area, UnitaryNormal[3], Rot_Flux,
			ProjGridVel = 0.0, ProjRotVel = 0.0, Density, a2, phi;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();

	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
			for (iDim = 0; iDim < nDim; iDim++) UnitaryNormal[iDim] = -Normal[iDim]/Area;

			/*--- Set the residual using the pressure ---*/
			if (incompressible) Pressure = node[iPoint]->GetSolution(0);
			else Pressure = node[iPoint]->GetPressure(incompressible);

			Residual[0] = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Residual[iDim+1] = Pressure*UnitaryNormal[iDim]*Area;
			if (!incompressible) Residual[nVar-1] = 0.0;

			/*--- Adjustment to energy equation for a rotating frame ---*/
			if (rotating_frame) {
				ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux();
				Residual[nVar-1] = Pressure*ProjRotVel;
			}

			/*--- Adjustment to energy equation due to grid motion ---*/
			if (grid_movement) {
				ProjGridVel = 0.0;
				GridVel = geometry->node[iPoint]->GetGridVel();
				for (iDim = 0; iDim < nDim; iDim++)
					ProjGridVel += GridVel[iDim]*UnitaryNormal[iDim]*Area;
				Residual[nVar-1] = Pressure*ProjGridVel;
			}

			/*--- Add value to the residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- Form Jacobians for implicit computations ---*/
			if (implicit) {
				if (incompressible)  {
					Density = node[iPoint]->GetDensityInc();
					for (iVar = 0; iVar < nVar; iVar++) {
						for (jVar = 0; jVar < nVar; jVar++)
							Jacobian_i[iVar][jVar] = 0.0;
					}
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_i[iDim+1][0] = -Normal[iDim];
					Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
				}
				else {
					a2 = Gamma-1.0;
					phi = 0.5*a2*node[iPoint]->GetVelocity2();

					for (iVar = 0; iVar < nVar; iVar++) {
						Jacobian_i[0][iVar] = 0.0;
						Jacobian_i[nDim+1][iVar] = 0.0;
					}
					for (iDim = 0; iDim < nDim; iDim++) {
						Jacobian_i[iDim+1][0] = -phi*Normal[iDim];
						for (jDim = 0; jDim < nDim; jDim++)
							Jacobian_i[iDim+1][jDim+1] = a2*node[iPoint]->GetVelocity(jDim, incompressible)*Normal[iDim];
						Jacobian_i[iDim+1][nDim+1] = -a2*Normal[iDim];
					}
					if (grid_movement) {
						ProjGridVel = 0.0;
						GridVel = geometry->node[iPoint]->GetGridVel();
						for (iDim = 0; iDim < nDim; iDim++)
							ProjGridVel += GridVel[iDim]*UnitaryNormal[iDim]*Area;
						for (jDim = 0; jDim < nDim; jDim++)
							Jacobian_i[nDim+1][jDim+1] = a2*node[iPoint]->GetVelocity(jDim, incompressible)*ProjGridVel;
						Jacobian_i[nDim+1][nDim+1] = -a2*ProjGridVel;
					}
					if (rotating_frame) {
						Rot_Flux = -geometry->vertex[val_marker][iVertex]->GetRotFlux();
						for (jDim = 0; jDim < nDim; jDim++)
							Jacobian_i[nDim+1][jDim+1] = a2*node[iPoint]->GetVelocity(jDim, incompressible)*Rot_Flux;
						Jacobian_i[nDim+1][nDim+1] = -a2*Rot_Flux;
					}
					Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
				}
			}
		}
	}

}

void CEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container,
                                  CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iVar, iDim, nPrimVar;
	unsigned long iVertex, iPoint, Point_Normal;
	double Pressure, Density, Height, *Coord;
    
	bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool rotating_frame     = config->GetRotating_Frame();
	bool grid_movement      = config->GetGrid_Movement();
	bool incompressible     = config->GetIncompressible();
	bool gravity            = config->GetGravityForce();
	double FreeSurface_Zero = config->GetFreeSurface_Zero();
	double PressFreeSurface = GetPressure_Inf();
	double Froude           = config->GetFroude();
    double Gas_Constant       = config->GetGas_ConstantND();
	bool viscous              = config->GetViscous();
    
    if (incompressible) nPrimVar = nDim+2;
	else nPrimVar = nDim+3;
    
	double *U_domain = new double[nVar]; double *U_infty = new double[nVar];
    double *V_domain = new double[nPrimVar]; double *V_infty = new double[nPrimVar];
	double *Normal = new double[nDim];
    
	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
            
			/*--- Index of the closest interior node ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
            
			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
			conv_numerics->SetNormal(Normal);
            
			/*--- Retrieve solution at the farfield boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++) U_domain[iVar] = node[iPoint]->GetSolution(iVar);
			for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = node[iPoint]->GetPrimVar(iVar);
            
			/*--- Construct solution state at infinity (farfield) ---*/
			if (incompressible) {
                
				Density = node[iPoint]->GetDensityInc();
				Height = geometry->node[iPoint]->GetCoord(nDim-1);
				Coord = geometry->node[iPoint]->GetCoord();
				if (gravity) Pressure = PressFreeSurface + Density*(FreeSurface_Zero-Height)/(Froude*Froude);
				else Pressure = PressFreeSurface;
                
				U_infty[0] = Pressure;
				for (iDim = 0; iDim < nDim; iDim++)
					U_infty[iDim+1] = GetVelocity_Inf(iDim)*Density;
                
				V_infty[0] = Density;
				for (iDim = 0; iDim < nDim; iDim++)
					V_infty[iDim+1] = GetVelocity_Inf(iDim);
                
			} else {
                
				U_infty[0] = GetDensity_Inf();
				for (iDim = 0; iDim < nDim; iDim++)
					U_infty[iDim+1] = GetDensity_Velocity_Inf(iDim);
				U_infty[nDim+1] = GetDensity_Energy_Inf();
                
                V_infty[0] = GetPressure_Inf() / ( Gas_Constant * GetDensity_Inf());
				for (iDim = 0; iDim < nDim; iDim++)
					V_infty[iDim+1] = GetVelocity_Inf(iDim);
				V_infty[nDim+1] = GetPressure_Inf();
				V_infty[nDim+2] = GetDensity_Inf();
                
			}
            
			/*--- Set various quantities in the solver class ---*/
			conv_numerics->SetConservative(U_domain, U_infty);
            
			if (incompressible) {
				conv_numerics->SetDensityInc(node[iPoint]->GetDensityInc(), node[iPoint]->GetDensityInc());
				conv_numerics->SetBetaInc2(node[iPoint]->GetBetaInc2(), node[iPoint]->GetBetaInc2());
				conv_numerics->SetCoord(Coord, Coord);
			}
            
			if (rotating_frame) {
				conv_numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
				conv_numerics->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
            
			if (grid_movement) {
				conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
			}
            
			/*--- Compute the convective residual residual using an upwind scheme ---*/
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
			LinSysRes.AddBlock(iPoint, Residual);
            
			/*--- Jacobian contribution for implicit integration ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
            
			/*--- Roe Turkel preconditioning, set the value of beta ---*/
			if ((config->GetKind_Upwind() == ROE_TURKEL_2ND) || (config->GetKind_Upwind() == ROE_TURKEL_1ST)) {
				node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
			}
            
			/*--- Viscous contribution ---*/
			if (viscous) {
                
				/*--- Set the normal vector and the coordinates ---*/
				visc_numerics->SetNormal(Normal);
				visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
                
				/*--- Primitive variables, and gradient ---*/
				visc_numerics->SetPrimitive(V_domain, V_infty);
				visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
                
				/*--- Laminar and eddy viscosity ---*/
				if (incompressible) visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosityInc(), node[iPoint]->GetLaminarViscosityInc());
				else visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());
				visc_numerics->SetEddyViscosity(node[iPoint]->GetEddyViscosity(), node[iPoint]->GetEddyViscosity());
                
				/*--- Turbulent kinetic energy ---*/
				if (config->GetKind_Turb_Model() == SST)
					visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
                
				/*--- Compute and update residual ---*/
				visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
				LinSysRes.SubtractBlock(iPoint, Residual);
                
				/*--- Jacobian contribution for implicit integration ---*/
				if (implicit)
					Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
                
			}
            
		}
	}
    
	/*--- Free locally allocated memory ---*/
	delete [] U_domain;
	delete [] U_infty;
    delete [] V_domain;
	delete [] V_infty;
	delete [] Normal;
    
}

void CEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                              CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iVar, iDim, nPrimVar;
	unsigned long iVertex, iPoint, Point_Normal;
	double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
	Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
	alpha, aa, bb, cc, dd, Area, UnitaryNormal[3], Density_Inlet;
    
	bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool rotating_frame       = config->GetRotating_Frame();
	bool grid_movement        = config->GetGrid_Movement();
	bool incompressible       = config->GetIncompressible();
	double Two_Gamma_M1       = 2.0/Gamma_Minus_One;
	double Gas_Constant       = config->GetGas_ConstantND();
	unsigned short Kind_Inlet = config->GetKind_Inlet();
	string Marker_Tag         = config->GetMarker_All_Tag(val_marker);
	bool viscous              = config->GetViscous();
	bool freesurface = config->GetFreeSurface();
    bool gravity = (config->GetGravityForce());
    
    if (incompressible) nPrimVar = nDim+2;
	else nPrimVar = nDim+3;
    
	double *U_domain = new double[nVar];      double *U_inlet = new double[nVar];
	double *V_domain = new double[nPrimVar];  double *V_inlet = new double[nPrimVar];
	double *Normal = new double[nDim];
    
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
			conv_numerics->SetNormal(Normal);
            
			Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
            
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;
            
			/*--- Retrieve solution at this boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++) U_domain[iVar] = node[iPoint]->GetSolution(iVar);
            for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = node[iPoint]->GetPrimVar(iVar);
            
			/*--- Build the fictitious intlet state based on characteristics ---*/
			if (incompressible) {
                
				/*--- Pressure and density using the internal value ---*/
				U_inlet[0] = node[iPoint]->GetSolution(0);
				Density_Inlet = node[iPoint]->GetDensityInc();
                V_inlet[0] = node[iPoint]->GetDensityInc();
                
				/*--- The velocity is computed from the infinity values ---*/
				for (iDim = 0; iDim < nDim; iDim++) {
					U_inlet[iDim+1] = GetVelocity_Inf(iDim)*Density_Inlet;
                    V_inlet[iDim+1] = GetVelocity_Inf(iDim);
                }
                
				/*--- The y/z velocity is interpolated due to the
                 free surface effect on the pressure ---*/
				if (freesurface) {
                    U_inlet[nDim] = node[iPoint]->GetSolution(nDim);
                    V_inlet[nDim] = node[iPoint]->GetPrimVar(nDim);
                }
                
			} else {
                
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
                        if (gravity) P_Total = config->GetInlet_Ptotal(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;
                        else P_Total  = config->GetInlet_Ptotal(Marker_Tag);
                        T_Total  = config->GetInlet_Ttotal(Marker_Tag);
                        Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);
                        
                        /*--- Non-dim. the inputs if necessary. ---*/
                        P_Total /= config->GetPressure_Ref();
                        T_Total /= config->GetTemperature_Ref();
                        
                        /*--- Store primitives and set some variables for clarity. ---*/
                        Density = U_domain[0];
                        Velocity2 = 0.0;
                        for (iDim = 0; iDim < nDim; iDim++) {
                            Velocity[iDim] = U_domain[iDim+1]/Density;
                            Velocity2 += Velocity[iDim]*Velocity[iDim];
                        }
                        Energy      = U_domain[nVar-1]/Density;
                        Pressure    = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
                        H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
                        SoundSpeed2 = Gamma*Pressure/Density;
                        
                        /*--- Compute the acoustic Riemann invariant that is extrapolated
                         from the domain interior. ---*/
                        Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
                        for (iDim = 0; iDim < nDim; iDim++)
                            Riemann += Velocity[iDim]*UnitaryNormal[iDim];
                        
                        /*--- Total speed of sound ---*/
                        SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;
                        
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
                        U_inlet[0] = Density;
                        for (iDim = 0; iDim < nDim; iDim++)
                            U_inlet[iDim+1] = Velocity[iDim]*Density;
                        U_inlet[nDim+1] = Energy*Density;
                        
                        /*--- Primitive variables, using the derived quantities ---*/
                        V_inlet[0] = Temperature;
                        for (iDim = 0; iDim < nDim; iDim++)
                            V_inlet[iDim+1] = Velocity[iDim];
                        V_inlet[nDim+1] = Pressure;
                        V_inlet[nDim+2] = Density;
                        
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
                        SoundSpeed2 = Gamma*Pressure/U_domain[0];
                        
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
                        U_inlet[0] = Density;
                        for (iDim = 0; iDim < nDim; iDim++)
                            U_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim]*Density;
                        U_inlet[nDim+1] = Energy*Density;
                        
                        /*--- Primitive variables, using the derived quantities ---*/
                        V_inlet[0] = Pressure / ( Gas_Constant * Density);
                        for (iDim = 0; iDim < nDim; iDim++)
                            V_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim];
                        V_inlet[nDim+1] = Pressure;
                        V_inlet[nDim+2] = Density;
                        
                        break;
				}
			}
            
			/*--- Set various quantities in the solver class ---*/
			conv_numerics->SetConservative(U_domain, U_inlet);
            
			if (incompressible) {
				conv_numerics->SetDensityInc(node[iPoint]->GetDensityInc(), Density_Inlet);
				conv_numerics->SetBetaInc2(node[iPoint]->GetBetaInc2(), node[iPoint]->GetBetaInc2());
				conv_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
			}
			if (rotating_frame) {
				conv_numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
				conv_numerics->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
			if (grid_movement)
				conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
            
			/*--- Compute the residual using an upwind scheme ---*/
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
			LinSysRes.AddBlock(iPoint, Residual);
            
			/*--- Jacobian contribution for implicit integration ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
            
			/*--- Roe Turkel preconditioning, set the value of beta ---*/
			if ((config->GetKind_Upwind() == ROE_TURKEL_2ND) || (config->GetKind_Upwind() == ROE_TURKEL_1ST)) {
				node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
			}
            
			/*--- Viscous contribution ---*/
			if (viscous) {
                
				/*--- Set the normal vector and the coordinates ---*/
				visc_numerics->SetNormal(Normal);
				visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
                
				/*--- Primitive variables, and gradient ---*/
				visc_numerics->SetPrimitive(V_domain, V_inlet);
				visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
                
				/*--- Laminar and eddy viscosity ---*/
				if (incompressible) visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosityInc(), node[iPoint]->GetLaminarViscosityInc());
				else visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());
				visc_numerics->SetEddyViscosity(node[iPoint]->GetEddyViscosity(), node[iPoint]->GetEddyViscosity());
                
				/*--- Turbulent kinetic energy ---*/
				if (config->GetKind_Turb_Model() == SST)
					visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
                
				/*--- Compute and update residual ---*/
				visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
				LinSysRes.SubtractBlock(iPoint, Residual);
                
				/*--- Jacobian contribution for implicit integration ---*/
				if (implicit)
					Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
                
			}
            
		}
	}
    
	/*--- Free locally allocated memory ---*/
	delete [] U_domain;
	delete [] U_inlet;
    delete [] V_domain;
	delete [] V_inlet;
	delete [] Normal;
    
}

void CEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                               CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iVar, iDim, nPrimVar;
	unsigned long iVertex, iPoint, Point_Normal;
	double LevelSet, Density_Outlet, Pressure, P_Exit, Velocity[3],
	Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
	Area, UnitaryNormal[3], Height;
    
	bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    double Gas_Constant     = config->GetGas_ConstantND();
	bool incompressible     = config->GetIncompressible();
	bool rotating_frame     = config->GetRotating_Frame();
	bool grid_movement      = config->GetGrid_Movement();
	double FreeSurface_Zero = config->GetFreeSurface_Zero();
	double PressFreeSurface = GetPressure_Inf();
	double Froude           = config->GetFroude();
	double epsilon          = config->GetFreeSurface_Thickness();
	double RatioDensity     = config->GetRatioDensity();
	string Marker_Tag       = config->GetMarker_All_Tag(val_marker);
	bool viscous              = config->GetViscous();
    bool freesurface = config->GetFreeSurface();
    bool gravity = (config->GetGravityForce());
    
    if (incompressible) nPrimVar = nDim+2;
	else nPrimVar = nDim+3;
    
	double *U_domain = new double[nVar];      double *U_outlet = new double[nVar];
    double *V_domain = new double[nPrimVar];  double *V_outlet = new double[nPrimVar];
	double *Normal = new double[nDim];
    
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
			conv_numerics->SetNormal(Normal);
            
			Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
            
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;
            
			/*--- Current solution at this boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++) U_domain[iVar] = node[iPoint]->GetSolution(iVar);
            for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = node[iPoint]->GetPrimVar(iVar);
            
			/*--- Build the fictitious intlet state based on characteristics ---*/
			if (incompressible) {
                
                
				if (freesurface) {
                    
					/*--- Density computation at the exit using the level set function ---*/
					Height = geometry->node[iPoint]->GetCoord(nDim-1);
					LevelSet = Height - FreeSurface_Zero;
                    
					/*--- Pressure computation the density at the exit (imposed) ---*/
					if (LevelSet < -epsilon) Density_Outlet = config->GetDensity_FreeStreamND();
					if (LevelSet > epsilon) Density_Outlet = RatioDensity*config->GetDensity_FreeStreamND();
					U_outlet[0] = PressFreeSurface + Density_Outlet*((FreeSurface_Zero-Height)/(Froude*Froude));
                    V_domain[0] = Density_Outlet;
                    
					/*--- Neumman condition in the interface for the pressure and density ---*/
					if (fabs(LevelSet) <= epsilon) {
						U_outlet[0] = node[Point_Normal]->GetSolution(0);
						Density_Outlet = node[Point_Normal]->GetDensityInc();
                        V_domain[0] = Density_Outlet;
					}
                    
				} else {
                    
					/*--- Imposed pressure and density ---*/
					Density_Outlet = GetDensity_Inf();
					U_outlet[0] = GetPressure_Inf();
                    V_domain[0] = Density_Outlet;
                    
				}
                
				/*--- Neumman condition for the velocity ---*/
				for (iDim = 0; iDim < nDim; iDim++) {
					U_outlet[iDim+1] = node[Point_Normal]->GetSolution(iDim+1);
                    V_outlet[iDim+1] = node[Point_Normal]->GetPrimVar(iDim+1);
                }
                
			}
			else {
                
				/*--- Retrieve the specified back pressure for this outlet. ---*/
                if (gravity) P_Exit = config->GetOutlet_Pressure(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;
                else P_Exit = config->GetOutlet_Pressure(Marker_Tag);
                
				/*--- Non-dim. the inputs if necessary. ---*/
				P_Exit = P_Exit/config->GetPressure_Ref();
                
				/*--- Check whether the flow is supersonic at the exit. The type
                 of boundary update depends on this. ---*/
				Density = U_domain[0];
				Velocity2 = 0.0; Vn = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Velocity[iDim] = U_domain[iDim+1]/Density;
					Velocity2 += Velocity[iDim]*Velocity[iDim];
					Vn += Velocity[iDim]*UnitaryNormal[iDim];
				}
				Energy     = U_domain[nVar-1]/Density;
				Pressure   = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
				SoundSpeed = sqrt(Gamma*Pressure/Density);
				Mach_Exit  = sqrt(Velocity2)/SoundSpeed;
                
				if (Mach_Exit >= 1.0) {
                    
					/*--- Supersonic exit flow: there are no incoming characteristics,
                     so no boundary condition is necessary. Set outlet state to current
                     state so that upwinding handles the direction of propagation. ---*/
					for (iVar = 0; iVar < nVar; iVar++) U_outlet[iVar] = U_domain[iVar];
                    for (iVar = 0; iVar < nPrimVar; iVar++) V_outlet[iVar] = V_domain[iVar];
                    
				} else {
                    
					/*--- Subsonic exit flow: there is one incoming characteristic,
                     therefore one variable can be specified (back pressure) and is used
                     to update the conservative variables. Compute the entropy and the
                     acoustic Riemann variable. These invariants, as well as the
                     tangential velocity components, are extrapolated. Adapted from an
                     original implementation in the Stanford University multi-block
                     (SUmb) solver in the routine bcSubsonicOutflow.f90 by Edwin van
                     der Weide, last modified 09-10-2007. ---*/
                    
					Entropy = Pressure*pow(1.0/Density,Gamma);
					Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;
                    
					/*--- Compute the new fictious state at the outlet ---*/
					Density    = pow(P_Exit/Entropy,1.0/Gamma);
					Pressure   = P_Exit;
					SoundSpeed = sqrt(Gamma*P_Exit/Density);
					Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
					Velocity2  = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) {
						Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitaryNormal[iDim];
						Velocity2 += Velocity[iDim]*Velocity[iDim];
					}
					Energy  = P_Exit/(Density*Gamma_Minus_One) + 0.5*Velocity2;
                    
					/*--- Conservative variables, using the derived quantities ---*/
					U_outlet[0] = Density;
					for (iDim = 0; iDim < nDim; iDim++)
						U_outlet[iDim+1] = Velocity[iDim]*Density;
					U_outlet[nDim+1] = Energy*Density;
                    
                    /*--- Conservative variables, using the derived quantities ---*/
                    V_outlet[0] = Pressure / ( Gas_Constant * Density);
					for (iDim = 0; iDim < nDim; iDim++)
						V_outlet[iDim+1] = Velocity[iDim];
					V_outlet[nDim+1] = Pressure;
                    V_outlet[nDim+2] = Density;
                    
				}
			}
            
			/*--- Set various quantities in the solver class ---*/
			conv_numerics->SetConservative(U_domain, U_outlet);
            
			if (incompressible) {
				conv_numerics->SetDensityInc(node[iPoint]->GetDensityInc(), Density_Outlet);
				conv_numerics->SetBetaInc2(node[iPoint]->GetBetaInc2(), node[iPoint]->GetBetaInc2());
				conv_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
			}
			if (rotating_frame) {
				conv_numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(), geometry->node[iPoint]->GetRotVel());
				conv_numerics->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
			if (grid_movement)
				conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
            
			/*--- Compute the residual using an upwind scheme ---*/
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
			LinSysRes.AddBlock(iPoint, Residual);
            
			/*--- Jacobian contribution for implicit integration ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
            
			/*--- Roe Turkel preconditioning, set the value of beta ---*/
			if ((config->GetKind_Upwind() == ROE_TURKEL_2ND) || (config->GetKind_Upwind() == ROE_TURKEL_1ST)) {
				node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
			}
            
			/*--- Viscous contribution ---*/
			if (viscous) {
                
				/*--- Set the normal vector and the coordinates ---*/
				visc_numerics->SetNormal(Normal);
				visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
                
				/*--- Primitive variables, and gradient ---*/
				visc_numerics->SetPrimitive(V_domain, V_outlet);
				visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
                
				/*--- Laminar and eddy viscosity ---*/
				if (incompressible) visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosityInc(), node[iPoint]->GetLaminarViscosityInc());
				else visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());
				visc_numerics->SetEddyViscosity(node[iPoint]->GetEddyViscosity(), node[iPoint]->GetEddyViscosity());
                
				/*--- Turbulent kinetic energy ---*/
				if (config->GetKind_Turb_Model() == SST)
					visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
                
				/*--- Compute and update residual ---*/
				visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
				LinSysRes.SubtractBlock(iPoint, Residual);
                
				/*--- Jacobian contribution for implicit integration ---*/
				if (implicit)
					Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
                
			}
            
		}
	}
    
	/*--- Free locally allocated memory ---*/
	delete [] U_domain;
	delete [] U_outlet;
    delete [] V_domain;
	delete [] V_outlet;
	delete [] Normal;
    
}

void CEulerSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                                         CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iDim, iVar, nPrimVar;
	unsigned long iVertex, iPoint, Point_Normal;
	double Density, Pressure, Temperature, Energy, *Velocity, Velocity2;
	double Gas_Constant = config->GetGas_ConstantND();
    
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();
	bool viscous              = config->GetViscous();
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
    
    if (incompressible) nPrimVar = nDim+2;
	else nPrimVar = nDim+3;
    
    double *U_inlet = new double[nVar]; double *U_domain = new double[nVar];
    double *V_inlet = new double[nPrimVar]; double *V_domain = new double[nPrimVar];
	double *Normal = new double[nDim];
    
	/*--- Supersonic inlet flow: there are no outgoing characteristics,
     so all flow variables can be imposed at the inlet.
     First, retrieve the specified values for the primitive variables. ---*/
	Temperature = config->GetInlet_Temperature(Marker_Tag);
	Pressure    = config->GetInlet_Pressure(Marker_Tag);
	Velocity    = config->GetInlet_Velocity(Marker_Tag);
    
	/*--- Density at the inlet from the gas law ---*/
	Density = Pressure/(Gas_Constant*Temperature);
    
	/*--- Non-dim. the inputs if necessary. ---*/
	Temperature = Temperature/config->GetTemperature_Ref();
	Pressure    = Pressure/config->GetPressure_Ref();
	Density     = Density/config->GetDensity_Ref();
	for (iDim = 0; iDim < nDim; iDim++)
		Velocity[iDim] = Velocity[iDim]/config->GetVelocity_Ref();
    
	/*--- Compute the energy from the specified state ---*/
	Velocity2 = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		Velocity2 += Velocity[iDim]*Velocity[iDim];
	Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;
    
	/*--- Conservative variables, using the derived quantities ---*/
	U_inlet[0] = Density;
    for (iDim = 0; iDim < nDim; iDim++)
        U_inlet[iDim+1] = Velocity[iDim]*Density;
    U_inlet[nDim+1] = Energy*Density;
    
    /*--- Primitive variables, using the derived quantities ---*/
	V_inlet[0] = Temperature;
    for (iDim = 0; iDim < nDim; iDim++)
        V_inlet[iDim+1] = Velocity[iDim];
    V_inlet[nDim+1] = Pressure;
    V_inlet[nDim+2] = Density;
    
	/*--- Loop over all the vertices on this boundary marker ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
            
            /*--- Index of the closest interior node ---*/
            Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
            
			/*--- Current solution at this boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++) U_domain[iVar] = node[iPoint]->GetSolution(iVar);
			for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = node[iPoint]->GetPrimVar(iVar);
            
			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            
			double Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
            
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;
            
			/*--- Set various quantities in the solver class ---*/
			conv_numerics->SetNormal(Normal);
			conv_numerics->SetConservative(U_domain, U_inlet);
			if (incompressible) {
				conv_numerics->SetDensityInc(node[iPoint]->GetDensityInc(),
                                           node[iPoint]->GetDensityInc());
				conv_numerics->SetBetaInc2(node[iPoint]->GetBetaInc2(),
                                         node[iPoint]->GetBetaInc2());
				conv_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                      geometry->node[iPoint]->GetCoord());
			}
			if (rotating_frame) {
				conv_numerics->SetRotVel(geometry->node[iPoint]->GetRotVel(),
                                       geometry->node[iPoint]->GetRotVel());
				conv_numerics->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
			if (grid_movement)
				conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                        geometry->node[iPoint]->GetGridVel());
            
			/*--- Compute the residual using an upwind scheme ---*/
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
			LinSysRes.AddBlock(iPoint, Residual);
            
			/*--- Jacobian contribution for implicit integration ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
            
            /*--- Viscous contribution ---*/
			if (viscous) {
                
				/*--- Set the normal vector and the coordinates ---*/
				visc_numerics->SetNormal(Normal);
				visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
                
				/*--- Primitive variables, and gradient ---*/
				visc_numerics->SetPrimitive(V_domain, V_inlet);
				visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
                
				/*--- Laminar and eddy viscosity ---*/
				visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());
				visc_numerics->SetEddyViscosity(node[iPoint]->GetEddyViscosity(), node[iPoint]->GetEddyViscosity());
                
				/*--- Turbulent kinetic energy ---*/
				if (config->GetKind_Turb_Model() == SST)
					visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
                
				/*--- Compute and update residual ---*/
				visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
				LinSysRes.SubtractBlock(iPoint, Residual);
                
				/*--- Jacobian contribution for implicit integration ---*/
				if (implicit)
					Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
            }
            
		}
	}
    
	/*--- Free locally allocated memory ---*/
	delete [] U_domain;
	delete [] U_inlet;
    delete [] V_domain;
	delete [] V_inlet;
	delete [] Normal;
    
}

void CEulerSolver::BC_Nacelle_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iVar, iDim, nPrimVar;
	unsigned long iVertex, iPoint, Point_Normal;
	double Pressure, P_Fan, Velocity[3], Velocity2, Entropy, Target_FanFace_Mach = 0.0, Density, Energy,
  Riemann, Area, UnitaryNormal[3], Vn, SoundSpeed, Vn_Exit, P_Fan_inc, P_Fan_old, M_Fan_old;
  
	double DampingFactor = config->GetDamp_Nacelle_Inflow();
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool viscous              = config->GetViscous();
  double Gas_Constant = config->GetGas_ConstantND();
  bool incompressible = config->GetIncompressible();
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
  
  if (incompressible) nPrimVar = nDim+2;
	else nPrimVar = nDim+3;
  
  double *U_domain = new double[nVar];      double *U_inflow = new double[nVar];
	double *V_domain = new double[nPrimVar];  double *V_inflow = new double[nPrimVar];
	double *Normal = new double[nDim];
  
	/*--- Retrieve the specified target fan face mach in the nacelle. ---*/
	Target_FanFace_Mach = config->GetFanFace_Mach_Target(Marker_Tag);
  
	/*--- Retrieve the old fan face pressure and mach number in the nacelle (this has been computed in a preprocessing). ---*/
  P_Fan_old = config->GetFanFace_Pressure(Marker_Tag);  // Note that has been computed by the code (non-dimensional).
  M_Fan_old = config->GetFanFace_Mach(Marker_Tag);

  
	/*--- Compute the Pressure increment ---*/
	P_Fan_inc = ((M_Fan_old/Target_FanFace_Mach) - 1.0) * config->GetPressure_FreeStreamND();
  
	/*--- Estimate the new fan face pressure ---*/
	P_Fan = (1.0 - DampingFactor)*P_Fan_old + DampingFactor * (P_Fan_old + P_Fan_inc);
  
	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
			Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
      
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;
      
			/*--- Current solution at this boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++) U_domain[iVar] = node[iPoint]->GetSolution(iVar);
			for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = node[iPoint]->GetPrimVar(iVar);
      
			/*--- Subsonic nacelle inflow: there is one incoming characteristic,
			 therefore one variable can be specified (back pressure) and is used
			 to update the conservative variables.
       
			 Compute the entropy and the acoustic variable. These
			 riemann invariants, as well as the tangential velocity components,
			 are extrapolated. ---*/
      Density = U_domain[0];
			Velocity2 = 0.0; Vn = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Velocity[iDim] = U_domain[iDim+1]/Density;
				Velocity2 += Velocity[iDim]*Velocity[iDim];
				Vn += Velocity[iDim]*UnitaryNormal[iDim];
			}
			Energy     = U_domain[nVar-1]/Density;
			Pressure   = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
			SoundSpeed = sqrt(Gamma*Pressure/Density);
			Entropy = Pressure*pow(1.0/Density,Gamma);
			Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;
      
			/*--- Compute the new fictious state at the outlet ---*/
			Density    = pow(P_Fan/Entropy,1.0/Gamma);
			Pressure   = P_Fan;
			SoundSpeed = sqrt(Gamma*P_Fan/Density);
			Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
			Velocity2  = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitaryNormal[iDim];
				Velocity2 += Velocity[iDim]*Velocity[iDim];
			}
      
			Energy  = P_Fan/(Density*Gamma_Minus_One) + 0.5*Velocity2;
      
			/*--- Conservative variables, using the derived quantities ---*/
			U_inflow[0] = Density;
      for (iDim = 0; iDim < nDim; iDim++)
        U_inflow[iDim+1] = Velocity[iDim]*Density;
			U_inflow[nDim+1] = Energy*Density;
      
      /*--- Conservative variables, using the derived quantities ---*/
			V_inflow[0] = Pressure / ( Gas_Constant * Density);
      for (iDim = 0; iDim < nDim; iDim++)
        V_inflow[iDim+1] = Velocity[iDim];
			V_inflow[nDim+1] = Pressure;
      V_inflow[nDim+2] = Density;
      
			/*--- Set various quantities in the solver class ---*/
			conv_numerics->SetNormal(Normal);
			conv_numerics->SetConservative(U_domain, U_inflow);
      
			/*--- Compute the residual using an upwind scheme ---*/
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
			LinSysRes.AddBlock(iPoint, Residual);
            
			/*--- Jacobian contribution for implicit integration ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
			if (viscous) {
        
				/*--- Set the normal vector and the coordinates ---*/
				visc_numerics->SetNormal(Normal);
				visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
				/*--- Primitive variables, and gradient ---*/
				visc_numerics->SetPrimitive(V_domain, V_inflow);
				visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
				/*--- Laminar and eddy viscosity ---*/
				visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());
				visc_numerics->SetEddyViscosity(node[iPoint]->GetEddyViscosity(), node[iPoint]->GetEddyViscosity());
        
				/*--- Turbulent kinetic energy ---*/
				if (config->GetKind_Turb_Model() == SST)
					visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
				/*--- Compute and update residual ---*/
				visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
				LinSysRes.SubtractBlock(iPoint, Residual);
                
				/*--- Jacobian contribution for implicit integration ---*/
				if (implicit)
					Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
			}
      
		}
	}
  
	delete [] U_domain;
	delete [] U_inflow;
  delete [] V_domain;
	delete [] V_inflow;
	delete [] Normal;
  
}

void CEulerSolver::BC_Nacelle_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iVar, iDim, nPrimVar;
	unsigned long iVertex, iPoint, Point_Normal;
	double P_Exhaust, T_Exhaust, Velocity[3], Velocity2, H_Exhaust, Temperature, Riemann, Area, UnitaryNormal[3], Pressure, Density, Energy, Mach2, SoundSpeed2, SoundSpeed_Exhaust2, Vel_Mag, alpha, aa, bb, cc, dd, Flow_Dir[3];
  
	double Gas_Constant = config->GetGas_ConstantND();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool viscous = config->GetViscous();
  bool incompressible = config->GetIncompressible();
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
  
  if (incompressible) nPrimVar = nDim+2;
	else nPrimVar = nDim+3;
  
  double *U_domain = new double[nVar];      double *U_exhaust = new double[nVar];
  double *V_domain = new double[nPrimVar];  double *V_exhaust = new double[nPrimVar];
	double *Normal = new double[nDim];
  
	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
			Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
      
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;
      
			/*--- Current solution at this boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++) U_domain[iVar] = node[iPoint]->GetSolution(iVar);
      for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = node[iPoint]->GetPrimVar(iVar);
      
			/*--- Subsonic inflow: there is one outgoing characteristic (u-c),
			 therefore we can specify all but one state variable at the inlet.
			 The outgoing Riemann invariant provides the final piece of info. ---*/
      
			/*--- Retrieve the specified total conditions for this inlet. ---*/
			P_Exhaust  = config->GetNozzle_Ptotal(Marker_Tag);
			T_Exhaust  = config->GetNozzle_Ttotal(Marker_Tag);
      
			/*--- Non-dim. the inputs if necessary. ---*/
			P_Exhaust /= config->GetPressure_Ref();
			T_Exhaust /= config->GetTemperature_Ref();
      
			/*--- Store primitives and set some variables for clarity. ---*/
			Density = U_domain[0];
			Velocity2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Velocity[iDim] = U_domain[iDim+1]/Density;
				Velocity2 += Velocity[iDim]*Velocity[iDim];
			}
			Energy      = U_domain[nVar-1]/Density;
			Pressure    = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
			H_Exhaust     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Exhaust;
			SoundSpeed2 = Gamma*Pressure/Density;
      
			/*--- Compute the acoustic Riemann invariant that is extrapolated
			 from the domain interior. ---*/
			Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
			for (iDim = 0; iDim < nDim; iDim++)
				Riemann += Velocity[iDim]*UnitaryNormal[iDim];
      
			/*--- Total speed of sound ---*/
			SoundSpeed_Exhaust2 = Gamma_Minus_One*(H_Exhaust - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;
      
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
			cc =  0.5*Gamma_Minus_One*Riemann*Riemann - 2.0*SoundSpeed_Exhaust2/Gamma_Minus_One;
      
			/*--- Solve quadratic equation for velocity magnitude. Value must
			 be positive, so the choice of root is clear. ---*/
			dd = bb*bb - 4.0*aa*cc;
			dd = sqrt(max(0.0,dd));
			Vel_Mag   = (-bb + dd)/(2.0*aa);
			Vel_Mag   = max(0.0,Vel_Mag);
			Velocity2 = Vel_Mag*Vel_Mag;
      
			/*--- Compute speed of sound from total speed of sound eqn. ---*/
			SoundSpeed2 = SoundSpeed_Exhaust2 - 0.5*Gamma_Minus_One*Velocity2;
      
			/*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
			Mach2 = Velocity2/SoundSpeed2;
			Mach2 = min(1.0,Mach2);
			Velocity2   = Mach2*SoundSpeed2;
			Vel_Mag     = sqrt(Velocity2);
			SoundSpeed2 = SoundSpeed_Exhaust2 - 0.5*Gamma_Minus_One*Velocity2;
      
			/*--- Compute new velocity vector at the inlet ---*/
			for (iDim = 0; iDim < nDim; iDim++)
				Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];
      
			/*--- Static temperature from the speed of sound relation ---*/
			Temperature = SoundSpeed2/(Gamma*Gas_Constant);
      
			/*--- Static pressure using isentropic relation at a point ---*/
			Pressure = P_Exhaust*pow((Temperature/T_Exhaust),Gamma/Gamma_Minus_One);
      
			/*--- Density at the inlet from the gas law ---*/
			Density = Pressure/(Gas_Constant*Temperature);
      
			/*--- Using pressure, density, & velocity, compute the energy ---*/
			Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;
      
			/*--- Conservative variables, using the derived quantities ---*/
			U_exhaust[0] = Density;
      for (iDim = 0; iDim < nDim; iDim++)
        U_exhaust[iDim+1] = Velocity[iDim]*Density;
			U_exhaust[nDim+1] = Energy*Density;
      
      /*--- Primitive variables, using the derived quantities ---*/
			V_exhaust[0] = Temperature;
      for (iDim = 0; iDim < nDim; iDim++)
        V_exhaust[iDim+1] = Velocity[iDim];
			V_exhaust[nDim+1] = Pressure;
			V_exhaust[nDim+2] = Density;
      
			/*--- Set various quantities in the solver class ---*/
			conv_numerics->SetNormal(Normal);
			conv_numerics->SetConservative(U_domain, U_exhaust);
      
			/*--- Compute the residual using an upwind scheme ---*/
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
			LinSysRes.AddBlock(iPoint, Residual);
      
			/*--- Jacobian contribution for implicit integration ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
			if (viscous) {
        
				/*--- Set the normal vector and the coordinates ---*/
				visc_numerics->SetNormal(Normal);
				visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
				/*--- Primitive variables, and gradient ---*/
				visc_numerics->SetPrimitive(V_domain, V_exhaust);
				visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
				/*--- Laminar and eddy viscosity ---*/
				visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());
				visc_numerics->SetEddyViscosity(node[iPoint]->GetEddyViscosity(), node[iPoint]->GetEddyViscosity());
        
				/*--- Turbulent kinetic energy ---*/
				if (config->GetKind_Turb_Model() == SST)
					visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
				/*--- Compute and update residual ---*/
				visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
				LinSysRes.SubtractBlock(iPoint, Residual);
        
				/*--- Jacobian contribution for implicit integration ---*/
				if (implicit)
					Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
			}
      
		}
	}
  
	delete [] U_domain;
	delete [] U_exhaust;
  delete [] V_domain;
	delete [] V_exhaust;
	delete [] Normal;
  
}

void CEulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iDim, iVar, jVar;
	double Pressure, *Normal, Density;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();

	/*--- Loop over all the vertices ---*/
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

			if (incompressible) Pressure = node[iPoint]->GetSolution(0);
			else Pressure = node[iPoint]->GetPressure(incompressible);

			/*--- Set the residual using the pressure ---*/
			Residual[0] = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Residual[iDim+1] = Pressure*UnitaryNormal[iDim]*Area;
			if (!incompressible) Residual[nVar-1] = 0.0;

			/*--- Add value to the residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- In case we are doing a implicit computation ---*/
			if (implicit) {
				if (incompressible)  {
					Density = node[iPoint]->GetDensityInc();
					for (iVar = 0; iVar < nVar; iVar++) {
						for (jVar = 0; jVar < nVar; jVar++)
							Jacobian_i[iVar][jVar] = 0.0;
					}
					for (iDim = 0; iDim < nDim; iDim++)
						Jacobian_i[iDim+1][0] = -Normal[iDim];
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
							Jacobian_i[iDim+1][jDim+1] = a2*node[iPoint]->GetVelocity(jDim, incompressible)*Normal[iDim];
						Jacobian_i[iDim+1][nDim+1] = -a2*Normal[iDim];
					}
					Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
				}
			}
		}
	}
}

void CEulerSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
		CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, jPoint;
	unsigned short iVar, iDim;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

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
				numerics->SetConservative(Solution_i, Solution_j);

				/*--- Set the normal vector ---*/
				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = -Vector[iDim];
				numerics->SetNormal(Vector);

				/*--- Add Residuals and Jacobians ---*/
				numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
				LinSysRes.AddBlock(iPoint, Residual);
				if (implicit) Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);

			}
		}
	}

#else

	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	double *Conserv_Var;
	bool compute;

	double *Buffer_Send_U = new double [nVar];
	double *Buffer_Receive_U = new double [nVar];

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
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Send_U[iVar] = Conserv_Var[iVar];
				MPI::COMM_WORLD.Bsend(Buffer_Send_U, nVar, MPI::DOUBLE, jProcessor, iPoint);
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
				if (jProcessor != rank)
					MPI::COMM_WORLD.Recv(Buffer_Receive_U, nVar, MPI::DOUBLE, jProcessor, jPoint);
				else {
					for (iVar = 0; iVar < nVar; iVar++)
						Buffer_Receive_U[iVar] = node[jPoint]->GetSolution(iVar);
				}

				/*--- Store the solution for both points ---*/
				for (iVar = 0; iVar < nVar; iVar++) {
					Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
					Solution_j[iVar] = Buffer_Receive_U[iVar];
				}

				numerics->SetConservative(Solution_i, Solution_j);

				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = -Vector[iDim];
				numerics->SetNormal(Vector);

				numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
				LinSysRes.AddBlock(iPoint, Residual);
				if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			}
		}
	}

	delete[] Buffer_Send_U;
	delete[] Buffer_Receive_U;

#endif

}

void CEulerSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
		CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, jPoint;
	unsigned short iVar, iDim;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

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
				numerics->SetConservative(Solution_i, Solution_j);

				/*--- Set the normal vector ---*/
				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = -Vector[iDim];
				numerics->SetNormal(Vector);

				/*--- Add Residuals and Jacobians ---*/
				numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
				LinSysRes.AddBlock(iPoint, Residual);
				if (implicit) Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}
		}
	}

#else

	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	double *Conserv_Var;
	bool compute;

	double *Buffer_Send_U = new double [nVar];
	double *Buffer_Receive_U = new double [nVar];

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

			/*--- We only send the information that belong to other boundary, -1 processor
			 means that the boundary condition is not applied ---*/
			if ((jProcessor != rank) && compute) {
				Conserv_Var = node[iPoint]->GetSolution();
				for (iVar = 0; iVar < nVar; iVar++)
					Buffer_Send_U[iVar] = Conserv_Var[iVar];
				MPI::COMM_WORLD.Bsend(Buffer_Send_U, nVar, MPI::DOUBLE, jProcessor, iPoint);
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
				if (jProcessor != rank)
					MPI::COMM_WORLD.Recv(Buffer_Receive_U, nVar, MPI::DOUBLE, jProcessor, jPoint);
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
				numerics->SetConservative(Solution_i, Solution_j);

				/*--- Set Normal ---*/
				geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = -Vector[iDim];
				numerics->SetNormal(Vector);

				/*--- Add Residuals and Jacobians ---*/
				numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
				LinSysRes.AddBlock(iPoint, Residual);
				if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			}
		}
	}

	delete[] Buffer_Send_U;
	delete[] Buffer_Receive_U;

#endif

}

void CEulerSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container,
		CConfig *config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, jPoint, total_index;
	unsigned short iVar;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

#ifdef NO_MPI

	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Find the associate pair to the original node ---*/
			jPoint = geometry->vertex[val_marker][iVertex]->GetDonorPoint();

			/*--- Store the solution for both points ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
				Solution_j[iVar] = node[jPoint]->GetSolution(iVar);
			}

			/*--- Set the residual, truncation error and solution value ---*/
			node[iPoint]->SetSolution_Old(Solution_j);
			node[iPoint]->SetSolution(Solution_j);
			LinSysRes.SetBlock_Zero(iPoint);
			node[iPoint]->SetRes_TruncErrorZero();

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
	double *Conserv_Var;

	double *Buffer_Send_U = new double [nVar];
	double *Buffer_Receive_U = new double [nVar];

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
				MPI::COMM_WORLD.Bsend(Buffer_Send_U, nVar, MPI::DOUBLE, jProcessor, iPoint);
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
				MPI::COMM_WORLD.Recv(Buffer_Receive_U, nVar, MPI::DOUBLE, jProcessor, jPoint);
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
			LinSysRes.SetBlock_Zero(iPoint);
			node[iPoint]->SetRes_TruncErrorZero();

			/*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
			if (implicit)
				for (iVar = 0; iVar < nVar; iVar++) {
					total_index = iPoint*nVar+iVar;
					Jacobian.DeleteValsRowi(total_index);
				}

		}
	}

	delete[] Buffer_Send_U;
	delete[] Buffer_Receive_U;

#endif

}

void CEulerSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {
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
			numerics->SetNormal(Vector);

			/*--- Compute the residual using an upwind scheme ---*/
			numerics->SetConservative(U_domain, Solution);
			numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

		}
	}
}

void CEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
		unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
	unsigned short iVar, jVar;
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1;
	double Volume_nM1, Volume_n, Volume_nP1, TimeStep;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool FlowEq = (RunTime_EqSystem == RUNTIME_FLOW_SYS);
	bool AdjEq = (RunTime_EqSystem == RUNTIME_ADJFLOW_SYS);
	bool incompressible = config->GetIncompressible();
	bool Grid_Movement = config->GetGrid_Movement();

	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n   = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();

		/*--- Volume at time n-1 and n ---*/
		if (Grid_Movement) {
			Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
			Volume_n = geometry->node[iPoint]->GetVolume_n();
			Volume_nP1 = geometry->node[iPoint]->GetVolume();
		}
		else {
			Volume_nM1 = geometry->node[iPoint]->GetVolume();
			Volume_n = geometry->node[iPoint]->GetVolume();
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

		if (incompressible && (FlowEq || AdjEq)) Residual[0] = 0.0;

		/*--- Add Residual ---*/
		LinSysRes.AddBlock(iPoint, Residual);

		if (implicit) {
			for (iVar = 0; iVar < nVar; iVar++) {
				for (jVar = 0; jVar < nVar; jVar++)
					Jacobian_i[iVar][jVar] = 0.0;

				if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
					Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
				if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
					Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
			}
			if (incompressible && (FlowEq || AdjEq)) Jacobian_i[0][0] = 0.0;
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}

}

void CEulerSolver::SetFlow_Displacement(CGeometry **flow_geometry, CVolumetricMovement *flow_grid_movement,
		CConfig *flow_config, CConfig *fea_config, CGeometry **fea_geometry, CSolver ***fea_solution) {
	unsigned short iMarker, iDim;
	unsigned long iVertex, iPoint;
	double *Coord, VarCoord[3];

#ifdef NO_MPI
	unsigned long iPoint_Donor;
	double *CoordDonor, *DisplacementDonor;

	for (iMarker = 0; iMarker < flow_config->GetnMarker_All(); iMarker++) {
		if (flow_config->GetMarker_All_Moving(iMarker) == YES) {
			for(iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
				iPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
				iPoint_Donor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetDonorPoint();
				Coord = flow_geometry[MESH_0]->node[iPoint]->GetCoord();
				CoordDonor = fea_geometry[MESH_0]->node[iPoint_Donor]->GetCoord();
				DisplacementDonor = fea_solution[MESH_0][FEA_SOL]->node[iPoint_Donor]->GetSolution();

				for (iDim = 0; iDim < nDim; iDim++)
					VarCoord[iDim] = (CoordDonor[iDim]+DisplacementDonor[iDim])-Coord[iDim];

				flow_geometry[MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
			}
		}
	}
	flow_grid_movement->SetVolume_Deformation(flow_geometry[MESH_0], flow_config, true);

#else

	int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
	double *Buffer_Send_Coord = new double [nDim];
	double *Buffer_Receive_Coord = new double [nDim];
	unsigned long jPoint;

	/*--- Do the send process, by the moment we are sending each
	 node individually, this must be changed ---*/
	for (iMarker = 0; iMarker < fea_config->GetnMarker_All(); iMarker++) {
		if (fea_config->GetMarker_All_Boundary(iMarker) == LOAD_BOUNDARY) {
			for(iVertex = 0; iVertex < fea_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
				iPoint = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();

				if (fea_geometry[MESH_0]->node[iPoint]->GetDomain()) {

					/*--- Find the associate pair to the original node (index and processor) ---*/
					jPoint = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
					jProcessor = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];

					/*--- We only send the pressure that belong to other boundary ---*/
					if (jProcessor != rank) {
						for (iDim = 0; iDim < nDim; iDim++)
							Buffer_Send_Coord[iDim] = fea_geometry[MESH_0]->node[iPoint]->GetCoord(iDim);

						MPI::COMM_WORLD.Bsend(Buffer_Send_Coord, nDim, MPI::DOUBLE, jProcessor, iPoint);
					}

				}
			}
		}
	}

	/*--- Now the loop is over the fea points ---*/
	for (iMarker = 0; iMarker < flow_config->GetnMarker_All(); iMarker++) {
		if ((flow_config->GetMarker_All_Boundary(iMarker) == EULER_WALL) &&
				(flow_config->GetMarker_All_Moving(iMarker) == YES)) {
			for(iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
				iPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
				if (flow_geometry[MESH_0]->node[iPoint]->GetDomain()) {

					/*--- Find the associate pair to the original node ---*/
					jPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
					jProcessor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];

					/*--- We only receive the information that belong to other boundary ---*/
					if (jProcessor != rank)
						MPI::COMM_WORLD.Recv(Buffer_Receive_Coord, nDim, MPI::DOUBLE, jProcessor, jPoint);
					else {
						for (iDim = 0; iDim < nDim; iDim++)
							Buffer_Send_Coord[iDim] = fea_geometry[MESH_0]->node[jPoint]->GetCoord(iDim);
					}

					/*--- Store the solution for both points ---*/
					Coord = flow_geometry[MESH_0]->node[iPoint]->GetCoord();

					for (iDim = 0; iDim < nDim; iDim++)
						VarCoord[iDim] = Buffer_Send_Coord[iDim]-Coord[iDim];

					flow_geometry[MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);


				}
			}
		}
	}
	delete[] Buffer_Send_Coord;
	delete[] Buffer_Receive_Coord;

#endif

}

void CEulerSolver::GetRestart(CGeometry *geometry, CConfig *config, unsigned short val_iZone) {

	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	/*--- Restart the solution from file information ---*/
	string restart_filename = config->GetSolution_FlowFileName();
	unsigned long iPoint, index, nFlowIter, adjIter, flowIter;
	char buffer[50];
	string UnstExt, text_line;
	ifstream restart_file;
	bool incompressible = config->GetIncompressible();
	bool grid_movement = config->GetGrid_Movement();
	unsigned short nZone = geometry->GetnZone();

	/*--- Multi-zone restart files. ---*/
	if (nZone > 1 && !(config->GetUnsteady_Simulation() == TIME_SPECTRAL)) {
		restart_filename.erase(restart_filename.end()-4, restart_filename.end());
		sprintf (buffer, "_%d.dat", int(val_iZone));
		UnstExt = string(buffer);
		restart_filename.append(UnstExt);
	}

	/*--- For the unsteady adjoint, we integrate backwards through
   physical time, so load in the direct solution files in reverse. ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		flowIter = val_iZone;
		restart_filename.erase(restart_filename.end()-4, restart_filename.end());
		if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.dat", int(val_iZone));
		if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.dat", int(val_iZone));
		UnstExt = string(buffer);
		restart_filename.append(UnstExt);
	} else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		nFlowIter = config->GetnExtIter() - 1;
		adjIter   = config->GetExtIter();
		flowIter  = nFlowIter - adjIter;
		restart_filename.erase (restart_filename.end()-4, restart_filename.end());
		if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
		if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
		if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
		if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
		if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
		UnstExt = string(buffer);
		restart_filename.append(UnstExt);
	} else {
		flowIter = config->GetExtIter();
		restart_filename.erase (restart_filename.end()-4, restart_filename.end());
		if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
		if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
		if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
		if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
		if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
		UnstExt = string(buffer);
		restart_filename.append(UnstExt);
	}

	/*--- Open the flow solution from the restart file ---*/
	if (rank == MASTER_NODE && val_iZone == ZONE_0)
		cout << "Reading in the direct flow solution from iteration " << flowIter << "." << endl;
	restart_file.open(restart_filename.data(), ios::in);

	/*--- In case there is no file ---*/
	if (restart_file.fail()) {
		cout << "There is no flow restart file!! " << restart_filename.data() << "."<< endl;
		cout << "Press any key to exit..." << endl;
		cin.get(); exit(1);
	}

	/*--- Store the previous solution (needed for aeroacoustic adjoint) ---*/
	if (config->GetExtIter() > 0)
		for (iPoint = 0; iPoint < nPoint; iPoint++)
			node[iPoint]->Set_Solution_time_n();

	/*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
	long *Global2Local = NULL;
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
	long iPoint_Local = 0; unsigned long iPoint_Global = 0;

	/*--- The first line is the header ---*/
	getline (restart_file, text_line);

	while (getline (restart_file,text_line)) {
		istringstream point_line(text_line);

		/*--- Retrieve local index. If this node from the restart file lives
     on a different processor, the value of iPoint_Local will be -1, as 
     initialized above. Otherwise, the local index for this node on the 
     current processor will be returned and used to instantiate the vars. ---*/
		iPoint_Local = Global2Local[iPoint_Global];
		if (iPoint_Local >= 0) {

			if (incompressible) {
				if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2];
				if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
			}
			else {
				if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
				if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
			}

			node[iPoint_Local]->SetSolution(Solution);

			/*--- If necessary, read in the grid velocities for the unsteady adjoint ---*/
			if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady() && grid_movement) {
				double Volume, GridVel[3];
				if (nDim == 2) point_line >> Volume >> GridVel[0] >> GridVel[1];
				if (nDim == 3) point_line >> Volume >> GridVel[0] >> GridVel[1] >> GridVel[2];
				if (iPoint_Local >= 0)
					for (unsigned short iDim = 0; iDim < nDim; iDim++)
						geometry->node[iPoint_Local]->SetGridVel(iDim, GridVel[iDim]);
			}

		}
		iPoint_Global++;
	}

	/*--- Set an average grid velocity at any halo nodes for the unsteady adjoint ---*/
	if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady() && grid_movement) {
		unsigned long jPoint;
		unsigned short nNeighbors;
		double AvgVel[3], *GridVel;
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
			AvgVel[0] = 0.0; AvgVel[1] = 0.0; AvgVel[2] = 0.0; nNeighbors = 0;
			/*--- Find & store any neighbor points to the sliding boundary in the donor zone (jZone). ---*/
			for (unsigned short iNeighbor = 0; iNeighbor < geometry->node[iPoint]->GetnPoint(); iNeighbor++) {
				jPoint = geometry->node[iPoint]->GetPoint(iNeighbor);
				if (geometry->node[jPoint]->GetDomain()) {
					GridVel = geometry->node[jPoint]->GetGridVel();
					for (unsigned short iDim = 0; iDim < nDim; iDim++) {
						AvgVel[iDim] += GridVel[iDim];
						nNeighbors++;
					}
				}
			}
			for (unsigned short iDim = 0; iDim < nDim; iDim++)
				geometry->node[iPoint]->SetGridVel(iDim, AvgVel[iDim]/(double)nNeighbors);
		}
	}

	/*--- Close the restart file ---*/
	restart_file.close();

	/*--- Free memory needed for the transformation ---*/
	delete [] Global2Local;

}

CNSSolver::CNSSolver(void) : CEulerSolver() {

	/*--- Array initialization ---*/
	CDrag_Visc = NULL;
	CLift_Visc = NULL;
	CMx_Visc = NULL;
	CMy_Visc = NULL;
	CMz_Visc = NULL;
	CFx_Visc = NULL;
	CFy_Visc = NULL;
	CFz_Visc = NULL;
	CEff_Visc = NULL;
	CMerit_Visc = NULL;
	CT_Visc = NULL;
	CQ_Visc = NULL;
	ForceViscous = NULL;
	MomentViscous = NULL;
	CSkinFriction = NULL;

}

CNSSolver::CNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CEulerSolver() {
	unsigned long iPoint, index, counter_local = 0, counter_global = 0;
	unsigned short iVar, iDim, iMarker;
  double Density, Velocity2, Pressure, Temperature;
  
	unsigned short nZone = geometry->GetnZone();
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	bool incompressible = config->GetIncompressible();
    double Gas_Constant = config->GetGas_ConstantND();

	/*--- Array initialization ---*/
	CDrag_Visc = NULL;
	CLift_Visc = NULL;
	CMx_Visc = NULL;
	CMy_Visc = NULL;
	CMz_Visc = NULL;
	CFx_Visc = NULL;
	CFy_Visc = NULL;
	CFz_Visc = NULL;
	CEff_Visc = NULL;
	CMerit_Visc = NULL;
	CT_Visc = NULL;
	CQ_Visc = NULL;
    Q_Visc = NULL;
    Maxq_Visc = NULL;
	ForceViscous = NULL;
	MomentViscous = NULL;
	CSkinFriction = NULL;
    
	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
    
	/*--- Set the gamma value ---*/
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;

	/*--- Define geometry constants in the solver structure ---*/
	nDim = geometry->GetnDim();
	if (incompressible) { nVar = nDim + 1; nPrimVar = nDim+2; nPrimVarGrad = nDim+2; }
	else { nVar = nDim + 2; nPrimVar = nDim+5; nPrimVarGrad = nDim+3; }
	nMarker = config->GetnMarker_All();
	nPoint = geometry->GetnPoint();
	nPointDomain = geometry->GetnPointDomain();
    
	/*--- Allocate the node variables ---*/
	node = new CVariable*[nPoint];
    
	/*--- Define some auxiliar vector related with the residual ---*/
	Residual      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
	Residual_RMS  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
	Residual_Max  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
	Point_Max  = new unsigned long[nVar]; for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
	Residual_i    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
	Residual_j    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
	Res_Conv      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
	Res_Visc      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
	Res_Sour      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;
    
	/*--- Define some auxiliary vectors related to the solution ---*/
	Solution   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
	Solution_i = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
	Solution_j = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;
    
	/*--- Define some auxiliary vectors related to the geometry ---*/
	Vector   = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
	Vector_i = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
	Vector_j = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
    
	/*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
		p1_Und_Lapl = new double [nPoint];
		p2_Und_Lapl = new double [nPoint];
	}
    
	if ((config->GetKind_Upwind_Flow() == ROE_TURKEL_2ND) || (config->GetKind_Upwind_Flow() == ROE_TURKEL_1ST)) {
		Precon_Mat_inv = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar ++)
			Precon_Mat_inv[iVar] = new double[nVar];
	}
    
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    
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
		if (rank == MASTER_NODE) cout << "Initialize jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
    
	} else {
		if (rank == MASTER_NODE)
			cout << "Explicit scheme. No jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
	}
    
	/*--- Computation of gradients by least squares ---*/
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
        
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
        
		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nDim+3];
		for (iVar = 0; iVar < nDim+3; iVar++)
			cvector[iVar] = new double [nDim];
	}
    
	/*--- Inviscid forces definition and coefficient in all the markers ---*/
	CPressure = new double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		CPressure[iMarker] = new double [geometry->nVertex[iMarker]];
    
	/*--- Heat tranfer in all the markers ---*/
	CHeatTransfer = new double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		CHeatTransfer[iMarker] = new double [geometry->nVertex[iMarker]];
    
	/*--- Y plus in all the markers ---*/
	YPlus = new double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		YPlus[iMarker] = new double [geometry->nVertex[iMarker]];
    
	/*--- Skin friction in all the markers ---*/
	CSkinFriction = new double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		CSkinFriction[iMarker] = new double [geometry->nVertex[iMarker]];
    
	/*--- Non dimensional coefficients ---*/
	ForceInviscid = new double[3];
	MomentInviscid = new double[3];
	CDrag_Inv = new double[nMarker];
	CLift_Inv = new double[nMarker];
	CSideForce_Inv = new double[nMarker];
	CMx_Inv = new double[nMarker];
	CMy_Inv = new double[nMarker];
	CMz_Inv = new double[nMarker];
	CEff_Inv = new double[nMarker];
	CFx_Inv = new double[nMarker];
	CFy_Inv = new double[nMarker];
	CFz_Inv = new double[nMarker];
    
	/*--- Rotational coefficients ---*/
	CMerit_Inv = new double[nMarker];
	CT_Inv = new double[nMarker];
	CQ_Inv = new double[nMarker];
    
	/*--- Supersonic coefficients ---*/
	CEquivArea_Inv = new double[nMarker];
	CNearFieldOF_Inv = new double[nMarker];
    
	/*--- Nacelle simulation ---*/
	FanFace_MassFlow  = new double[nMarker];
	Exhaust_MassFlow  = new double[nMarker];
    Exhaust_Area  = new double[nMarker];
	FanFace_Pressure  = new double[nMarker];
	FanFace_Mach  = new double[nMarker];
    FanFace_Area = new double[nMarker];
    
	/*--- Init total coefficients ---*/
	Total_CDrag = 0.0;	Total_CLift = 0.0;			Total_CSideForce = 0.0;
	Total_CMx = 0.0;		Total_CMy = 0.0;				Total_CMz = 0.0;
	Total_CEff = 0.0;		Total_CEquivArea = 0.0; Total_CNearFieldOF = 0.0;
	Total_CFx = 0.0;		Total_CFy = 0.0;				Total_CFz = 0.0; Total_CMerit = 0.0;
	Total_CT = 0.0;			Total_CQ = 0.0;         Total_Q = 0.0;
    Total_Maxq = 0.0;
    
	ForceViscous = new double[3];
	MomentViscous = new double[3];
	CDrag_Visc = new double[nMarker];
	CLift_Visc = new double[nMarker];
	CMx_Visc = new double[nMarker];
	CMy_Visc = new double[nMarker];
	CMz_Visc = new double[nMarker];
	CEff_Visc = new double[nMarker];
	CFx_Visc = new double[nMarker];
	CFy_Visc = new double[nMarker];
	CFz_Visc = new double[nMarker];
	CMerit_Visc = new double[nMarker];
	CT_Visc = new double[nMarker];
	CQ_Visc = new double[nMarker];
    Q_Visc = new double[nMarker];
    Maxq_Visc = new double[nMarker];
    
	/*--- Read farfield conditions from config ---*/
	Density_Inf   = config->GetDensity_FreeStreamND();
	Pressure_Inf  = config->GetPressure_FreeStreamND();
	Velocity_Inf  = config->GetVelocity_FreeStreamND();
	Energy_Inf    = config->GetEnergy_FreeStreamND();
	Viscosity_Inf = config->GetViscosity_FreeStreamND();
	Mach_Inf      = config->GetMach_FreeStreamND();
	Prandtl_Lam   = config->GetPrandtl_Lam();
	Prandtl_Turb  = config->GetPrandtl_Turb();
    
	/*--- Initializate Fan Face Pressure ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		FanFace_MassFlow[iMarker] = 0.0;
		Exhaust_MassFlow[iMarker] = 0.0;
        Exhaust_Area[iMarker] = 0.0;
		FanFace_Pressure[iMarker] = Pressure_Inf;
		FanFace_Mach[iMarker] = Mach_Inf;
	}
    
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
    
	/*--- Restart the solution from file information ---*/
	if (!restart || geometry->GetFinestMGLevel() == false || nZone > 1) {
        
		/*--- Restart the solution from infinity ---*/
		for (iPoint = 0; iPoint < nPoint; iPoint++)
			node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
	}
    
	else {
        
		/*--- Restart the solution from file information ---*/
		ifstream restart_file;
		string filename = config->GetSolution_FlowFileName();
		restart_file.open(filename.data(), ios::in);
        
		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
			cout << "There is no flow restart file!! " << filename.data() << "."<< endl;
			cout << "Press any key to exit..." << endl;
			cin.get(); exit(1);
		}
        
		/*--- In case this is a parallel simulation, we need to perform the
         Global2Local index transformation first. ---*/
		long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];
        
		/*--- First, set all indices to a negative value by default ---*/
		for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
			Global2Local[iPoint] = -1;
        
		/*--- Now fill array with the transform values only for local points ---*/
		for(iPoint = 0; iPoint < nPointDomain; iPoint++)
			Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
        
		/*--- Read all lines in the restart file ---*/
		long iPoint_Local; unsigned long iPoint_Global = 0; string text_line;
        
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
					if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2];
					if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
				}
				else {
					if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
					if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
				}
				node[iPoint_Local] = new CNSVariable(Solution, nDim, nVar, config);
			}
			iPoint_Global++;
		}
        
		/*--- Instantiate the variable class with an arbitrary solution
         at any halo/periodic nodes. The initial solution can be arbitrary,
         because a send/recv is performed immediately in the solver. ---*/
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++)
			node[iPoint] = new CNSVariable(Solution, nDim, nVar, config);
        
		/*--- Close the restart file ---*/
		restart_file.close();
        
		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
	}
  
  /*--- Check that the initial solution is physical ---*/
  if (!incompressible) {
    counter_local = 0;
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      Density = node[iPoint]->GetSolution(0);
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);
      Pressure    = Gamma_Minus_One*Density*(node[iPoint]->GetSolution(nDim+1)/Density-0.5*Velocity2);
      Temperature = Pressure / ( Gas_Constant * Density);
      if ((Pressure < 0.0) || (Temperature < 0.0)) {
        Solution[0] = Density_Inf;
        for (iDim = 0; iDim < nDim; iDim++)
          Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
        Solution[nDim+1] = Energy_Inf*Density_Inf;
        node[iPoint]->SetSolution(Solution);
        node[iPoint]->SetSolution_Old(Solution);
        counter_local++;
      }
    }
#ifndef NO_MPI
    MPI::COMM_WORLD.Reduce(&counter_local, &counter_global, 1, MPI::UNSIGNED_LONG, MPI::SUM, MASTER_NODE);
#else
    counter_global = counter_local;
#endif
    if ((rank == MASTER_NODE) && (counter_global != 0))
      cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  }
  
	/*--- For incompressible solver set the initial values for the density and viscosity,
	 unless a freesurface problem, this must be constant during the computation ---*/
	if (incompressible) {
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			node[iPoint]->SetDensityInc(Density_Inf);
			node[iPoint]->SetLaminarViscosityInc(Viscosity_Inf);
		}
	}
    
	/*--- Define solver parameters needed for execution of destructor ---*/
	if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
		space_centered = true;
	else space_centered = false;
    
	if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
	else euler_implicit = false;
    
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
	else least_squares = false;
    
	if ((config->GetKind_Upwind_Flow() == ROE_TURKEL_2ND) ||
        (config->GetKind_Upwind_Flow() == ROE_TURKEL_1ST)) roe_turkel = true;
	else roe_turkel = false;
    
	/*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);

}

CNSSolver::~CNSSolver(void) {
	unsigned short iMarker;

	if (CDrag_Visc != NULL) delete [] CDrag_Visc;
	if (CLift_Visc != NULL) delete [] CLift_Visc;
	if (CMx_Visc != NULL) delete [] CMx_Visc;
	if (CMy_Visc != NULL) delete [] CMy_Visc;
	if (CMz_Visc != NULL) delete [] CMz_Visc;
	if (CFx_Visc != NULL) delete [] CFx_Visc;
	if (CFy_Visc != NULL) delete [] CFy_Visc;
	if (CFz_Visc != NULL) delete [] CFz_Visc;
	if (CEff_Visc != NULL) delete [] CEff_Visc;
	if (CMerit_Visc != NULL) delete [] CMerit_Visc;
	if (CT_Visc != NULL) delete [] CT_Visc;
	if (CQ_Visc != NULL) delete [] CQ_Visc;
  if (Q_Visc != NULL) delete [] Q_Visc;
  if (Maxq_Visc != NULL) delete [] Maxq_Visc;
	if (ForceViscous != NULL) delete [] ForceViscous;
	if (MomentViscous != NULL) delete [] MomentViscous;

	if (CSkinFriction != NULL) {
		for (iMarker = 0; iMarker < nMarker; iMarker++) {
			delete CSkinFriction[iMarker];
		}
		delete [] CSkinFriction;
	}

}

void CNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
	double levelset;
    
	bool adjoint = config->GetAdjoint();
	bool freesurface = config->GetFreeSurface();
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool upwind_2nd = ((config->GetKind_Upwind_Flow() == ROE_2ND) || (config->GetKind_Upwind_Flow() == AUSM_2ND)
                       || (config->GetKind_Upwind_Flow() == HLLC_2ND) || (config->GetKind_Upwind_Flow() == ROE_TURKEL_2ND));
	bool center = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
	bool center_jst = center && config->GetKind_Centered_Flow() == JST;
	bool limiter = (config->GetKind_SlopeLimit_Flow() != NONE);
	bool incompressible = config->GetIncompressible();
	unsigned short turb_model = config->GetKind_Turb_Model();
	bool tkeNeeded = (turb_model == SST);
	double turb_ke = 0.0;
  bool jouleheating = config->GetJouleHeating();
    
    /*--- Compute nacelle inflow and exhaust properties ---*/
    GetNacelle_Properties(geometry, config, iMesh);
    
    /*--- Compute Joule heating ---*/
	if (jouleheating) geometry->SetGeometryPlanes(config);
    
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
        
		if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
        
		if (freesurface) levelset = solver_container[LEVELSET_SOL]->node[iPoint]->GetPrimVar(0);
		else levelset = 0.0;
        
		/*--- Set the primitive variables incompressible (dens, vx, vy, vz, beta)
         and compressible (temp, vx, vy, vz, press, dens, enthal, sos)---*/
		if (incompressible) node[iPoint]->SetPrimVar_Incompressible(Density_Inf, Viscosity_Inf, turb_ke, levelset, config);
		else node[iPoint]->SetPrimVar_Compressible(config, turb_ke);
        
		/*--- Set the value of the eddy viscosity ---*/
		if (turb_model != NONE)
			node[iPoint]->SetEddyViscosity(config->GetKind_Turb_Model(), solver_container[TURB_SOL]->node[iPoint]);
		else
			node[iPoint]->SetEddyViscosity(NONE, NULL);
        
		/*--- Initialize the convective, source and viscous residual vector ---*/
		LinSysRes.SetBlock_Zero(iPoint);
        
	}
    
	/*--- Upwind second order reconstruction ---*/
	if ((upwind_2nd) && (iMesh == MESH_0)) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
        
		/*--- Limiter computation ---*/
		if ((limiter) && (iMesh == MESH_0)) SetSolution_Limiter(geometry, config);
	}
    
	/*--- Artificial dissipation ---*/
	if (center) {
		SetMax_Eigenvalue(geometry, config);
		if ((center_jst) && (iMesh == MESH_0)) {
			SetDissipation_Switch(geometry, config);
			SetUndivided_Laplacian(geometry, config);
		}
	}
    
	/*--- Compute gradient of the primitive variables ---*/
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetPrimVar_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetPrimVar_Gradient_LS(geometry, config);
    
	/*--- Initialize the jacobian matrices ---*/
	if (implicit) Jacobian.SetValZero();
    
}

void CNSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) {
	double Mean_BetaInc2, *Normal, Area, Vol, Mean_SoundSpeed, Mean_ProjVel, Lambda, Local_Delta_Time, Local_Delta_Time_Visc, Mean_DensityInc,
	Global_Delta_Time = 1E6, Mean_LaminarVisc, Mean_EddyVisc, Mean_Density, Lambda_1, Lambda_2, K_v = 0.25, Global_Delta_UnstTimeND;
	unsigned long iEdge, iVertex, iPoint = 0, jPoint = 0;
	unsigned short iDim, iMarker;
	double ProjVel, ProjVel_i, ProjVel_j;

	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool rotating_frame = config->GetRotating_Frame();
	bool incompressible = config->GetIncompressible();
	bool grid_movement = config->GetGrid_Movement();
	double turb_model = config->GetKind_Turb_Model();
	bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
			(config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

	Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;

	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed and viscosity ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		node[iPoint]->SetMax_Lambda_Inv(0.0);
		node[iPoint]->SetMax_Lambda_Visc(0.0);

		/*--- Set the value of the eddy viscosity ---*/
		if (turb_model != NONE)
			node[iPoint]->SetEddyViscosity(config->GetKind_Turb_Model(), solver_container[TURB_SOL]->node[iPoint]);
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
			Mean_DensityInc = 0.5 * (node[iPoint]->GetDensityInc() + node[jPoint]->GetDensityInc());
			Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);
		}
		else {
			Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
			Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
		}

		/*--- Contribution due to a rotating frame ---*/
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

		/*--- Adjustment for grid movement ---*/
		if (grid_movement) {
			double *GridVel_i = geometry->node[iPoint]->GetGridVel();
			double *GridVel_j = geometry->node[jPoint]->GetGridVel();
			ProjVel_i = 0.0; ProjVel_j =0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				ProjVel_i += GridVel_i[iDim]*Normal[iDim];
				ProjVel_j += GridVel_j[iDim]*Normal[iDim];
			}
			Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j) ;
		}

		/*--- Inviscid contribution ---*/
		Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed ;
		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);

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

		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Visc(Lambda);

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
				Mean_DensityInc = node[iPoint]->GetDensityInc();
				Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);
			}
			else {
				Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
				Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
			}

			/*--- Contribution due to a rotating frame ---*/
			if (rotating_frame) {
				double *RotVel = geometry->node[iPoint]->GetRotVel();
				ProjVel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					ProjVel += RotVel[iDim]*Normal[iDim];
				Mean_ProjVel -= ProjVel;
			}

			/*--- Adjustment for grid movement ---*/
			if (grid_movement) {
				double *GridVel = geometry->node[iPoint]->GetGridVel();
				ProjVel = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					ProjVel += GridVel[iDim]*Normal[iDim];
				Mean_ProjVel -= ProjVel;
			}

			/*--- Inviscid contribution ---*/
			Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
			if (geometry->node[iPoint]->GetDomain()) {
				node[iPoint]->AddMax_Lambda_Inv(Lambda);
			}

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

			if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);

		}
	}

	/*--- Each element uses their own speed ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();
		Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
		Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Vol*Vol/ node[iPoint]->GetMax_Lambda_Visc();
		Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
		Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
		Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
		Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
		node[iPoint]->SetDelta_Time(Local_Delta_Time);
	}

	/*--- Check if there is any element with only one neighbor...
   a CV that is inside another CV ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		if (geometry->node[iPoint]->GetnPoint() == 1)
			node[iPoint]->SetDelta_Time(Min_Delta_Time);
	}

	/*--- For exact time solution use the minimum delta time of the whole mesh ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
#ifndef NO_MPI
		double rbuf_time, sbuf_time;
		sbuf_time = Global_Delta_Time;
		MPI::COMM_WORLD.Reduce(&sbuf_time, &rbuf_time, 1, MPI::DOUBLE, MPI::MIN, MASTER_NODE);
		MPI::COMM_WORLD.Bcast(&rbuf_time, 1, MPI::DOUBLE, MASTER_NODE);
		MPI::COMM_WORLD.Barrier();
		Global_Delta_Time = rbuf_time;
#endif
		for(iPoint = 0; iPoint < nPointDomain; iPoint++)
			node[iPoint]->SetDelta_Time(Global_Delta_Time);
	}

	/*--- Recompute the unsteady time step for the dual time stratey
	 if the unsteady CFL is diferent from 0 ---*/
	if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
		Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);
#ifndef NO_MPI
		double rbuf_time, sbuf_time;
		sbuf_time = Global_Delta_UnstTimeND;
		MPI::COMM_WORLD.Reduce(&sbuf_time, &rbuf_time, 1, MPI::DOUBLE, MPI::MIN, MASTER_NODE);
		MPI::COMM_WORLD.Bcast(&rbuf_time, 1, MPI::DOUBLE, MASTER_NODE);
		MPI::COMM_WORLD.Barrier();
		Global_Delta_UnstTimeND = rbuf_time;
#endif
		config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
	}

	/*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
	if (dual_time)
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
			if (!implicit) {
				Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
				/*--- Check if there is any element with only one neighbor...
				 a CV that is inside another CV ---*/
				if (geometry->node[iPoint]->GetnPoint() == 1) Local_Delta_Time = 0.0;
				node[iPoint]->SetDelta_Time(Local_Delta_Time);
			}
		}

}

void CNSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iPoint, jPoint, iEdge;
  
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points, coordinates and normal vector in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive variables, and gradient ---*/
    numerics->SetPrimitive(node[iPoint]->GetPrimVar(), node[jPoint]->GetPrimVar());
    numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[jPoint]->GetGradient_Primitive());
    
    /*--- Laminar and eddy viscosity ---*/
    if (incompressible) {
      numerics->SetDensityInc(node[iPoint]->GetDensityInc(), node[jPoint]->GetDensityInc());
      numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosityInc(), node[jPoint]->GetLaminarViscosityInc());
    }
    else
      numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[jPoint]->GetLaminarViscosity());
    
    numerics->SetEddyViscosity(node[iPoint]->GetEddyViscosity(), node[jPoint]->GetEddyViscosity());
    
    /*--- Turbulent kinetic energy ---*/
    if (config->GetKind_Turb_Model() == SST)
      numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0),
                                   solver_container[TURB_SOL]->node[jPoint]->GetSolution(0));
    
    /*--- Compute and update residual ---*/
    numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
    
    LinSysRes.SubtractBlock(iPoint, Res_Visc);
    LinSysRes.AddBlock(jPoint, Res_Visc);
    
    /*--- Implicit part ---*/
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    }
  }
  
}

void CNSSolver::Viscous_Forces(CGeometry *geometry, CConfig *config) {
	unsigned long iVertex, iPoint, iPointNormal;
	unsigned short Boundary, Monitoring, iMarker, iDim, jDim;
	double **Tau, Delta, Viscosity, **Grad_PrimVar, div_vel, *Normal, *TauElem, MomentDist[3], WallDist[3],
	*Coord, *Coord_Normal, *UnitaryNormal, *TauTangent, Area, WallShearStress, TauNormal, factor, RefVel2,
	RefDensity, GradTemperature, Density, Vel[3], VelNormal, VelTangMod, WallDistMod, FrictionVel, VelTang[3], HeatLoad;
    
	double Alpha        = config->GetAoA()*PI_NUMBER/180.0;
	double Beta         = config->GetAoS()*PI_NUMBER/180.0;
	double RefAreaCoeff = config->GetRefAreaCoeff();
	double RefLengthMoment = config->GetRefLengthMoment();
	double *Origin      = config->GetRefOriginMoment();
	double Gas_Constant = config->GetGas_ConstantND();
	double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
    
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();
    
    
	/*--- If we have a rotating frame problem or an unsteady problem with
     mesh motion, use special reference values for the force coefficients.
     Otherwise, use the freestream values, which is the standard convention. ---*/
    
	if (rotating_frame) {
        
        /*--- Use the rotational origin for computing moments ---*/
        Origin = config->GetRotAxisOrigin();
        
        /*--- Reference moment length is the rotor radius. Note that
         this assumes that the reference area was set to the rotor disk area. ---*/
        RefLengthMoment = sqrt(RefAreaCoeff/PI_NUMBER);
        
        /*--- Reference velocity is the rotational speed times rotor radius ---*/
        RefVel2 = (config->GetOmegaMag()*RefLengthMoment)*(config->GetOmegaMag()*RefLengthMoment);
        
    } else if (grid_movement) {
		double Gas_Constant = config->GetGas_ConstantND();
		double Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
		double Mach_Motion = config->GetMach_Motion();
		RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
        
	} else {
		double *Velocity_Inf = config->GetVelocity_FreeStreamND();
		RefVel2 = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
	}
    
	RefDensity  = config->GetDensity_FreeStreamND();
    
	factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
    
	/*-- Initialization --*/
	AllBound_CDrag_Visc = 0.0; AllBound_CLift_Visc = 0.0;
	AllBound_CMx_Visc = 0.0; AllBound_CMy_Visc = 0.0; AllBound_CMz_Visc = 0.0;
	AllBound_CFx_Visc = 0.0; AllBound_CFy_Visc = 0.0; AllBound_CFz_Visc = 0.0;
	AllBound_CEff_Visc = 0.0; AllBound_CMerit_Visc = 0.0;
	AllBound_CT_Visc = 0.0; AllBound_CQ_Visc = 0.0; AllBound_Q_Visc = 0.0;  AllBound_Maxq_Visc = 0.0;
    
	/*--- Vector and variables initialization ---*/
	UnitaryNormal      = new double [nDim];
	TauElem    = new double [nDim];
	TauTangent = new double [nDim];
	Tau        = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Tau[iDim]   = new double [nDim];
    
	/*--- Loop over the Navier-Stokes markers ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);
        
		if ((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL)) {
            
			for (iDim = 0; iDim < nDim; iDim++) ForceViscous[iDim] = 0.0;
			MomentViscous[0] = 0.0; MomentViscous[1] = 0.0; MomentViscous[2] = 0.0;
            HeatLoad = 0.0;
            
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
                
				Coord = geometry->node[iPoint]->GetCoord();
				Coord_Normal = geometry->node[iPointNormal]->GetCoord();
                
				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Grad_PrimVar = node[iPoint]->GetGradient_Primitive();
				if (incompressible) {
					Viscosity = node[iPoint]->GetLaminarViscosityInc();
					Density = node[iPoint]->GetDensityInc();
				}
				else {
					Viscosity = node[iPoint]->GetLaminarViscosity();
					Density = node[iPoint]->GetDensity();
				}
                
				Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
				for (iDim = 0; iDim < nDim; iDim++) {
					UnitaryNormal[iDim] = Normal[iDim]/Area;
					MomentDist[iDim] = Coord[iDim] - Origin[iDim];
				}
                
				div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_PrimVar[iDim+1][iDim];
                
				for (iDim = 0; iDim < nDim; iDim++) {
					for (jDim = 0 ; jDim < nDim; jDim++) {
						Delta = 0.0; if (iDim == jDim) Delta = 1.0;
						Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
                        TWO3*Viscosity*div_vel*Delta;
					}
					TauElem[iDim] = 0.0;
					for (jDim = 0; jDim < nDim; jDim++)
						TauElem[iDim] += Tau[iDim][jDim]*UnitaryNormal[jDim];
				}
                
				/*--- Compute wall shear stress (using the stress tensor) ---*/
				TauNormal = 0.0; for (iDim = 0; iDim < nDim; iDim++) TauNormal += TauElem[iDim] * UnitaryNormal[iDim];
				for (iDim = 0; iDim < nDim; iDim++) TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitaryNormal[iDim];
				WallShearStress = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallShearStress += TauTangent[iDim]*TauTangent[iDim];
				WallShearStress = sqrt(WallShearStress);
                
				/*--- Compute wall shear stress (using mu(delta u/delta y) ---*/
				for (iDim = 0; iDim < nDim; iDim++) Vel[iDim] = node[iPointNormal]->GetVelocity(iDim, incompressible);
				VelNormal = 0.0; for (iDim = 0; iDim < nDim; iDim++) VelNormal += Vel[iDim] * UnitaryNormal[iDim];
				for (iDim = 0; iDim < nDim; iDim++) VelTang[iDim] = Vel[iDim] - VelNormal*UnitaryNormal[iDim];
				VelTangMod = 0.0; for (iDim = 0; iDim < nDim; iDim++) VelTangMod += VelTang[iDim]*VelTang[iDim]; VelTangMod = sqrt(VelTangMod);
				for (iDim = 0; iDim < nDim; iDim++) WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
				WallDistMod = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallDistMod += WallDist[iDim]*WallDist[iDim]; WallDistMod = sqrt(WallDistMod);
				//				WallShearStress = Viscosity*VelTangMod/WallDistMod;
                
				/*--- Compute wall skin friction coefficient, and heat flux on the wall ---*/
				CSkinFriction[iMarker][iVertex] = WallShearStress / (0.5*RefDensity*RefVel2);
                
				/*--- Compute y+ and non-dimensional velocity ---*/
				FrictionVel = sqrt(fabs(WallShearStress)/Density);
				YPlus[iMarker][iVertex] = WallDistMod*FrictionVel/(Viscosity/Density);
                
				/*--- Compute heat flux on the wall ---*/
				GradTemperature = 0.0; for (iDim = 0; iDim < nDim; iDim++) GradTemperature +=  Grad_PrimVar[0][iDim]*(-Normal[iDim]);
				CHeatTransfer[iMarker][iVertex] = (Cp * Viscosity/PRANDTL)*GradTemperature/(0.5*RefDensity*RefVel2);
                HeatLoad += CHeatTransfer[iMarker][iVertex];
                
                if (CHeatTransfer[iMarker][iVertex]/Area > Maxq_Visc[iMarker])
                    Maxq_Visc[iMarker] = CHeatTransfer[iMarker][iVertex]/Area;
                
				/*--- Compute viscous forces, and moment using the stress tensor ---*/
				if ((geometry->node[iPoint]->GetDomain()) && (Monitoring == YES)) {
                    
					for (iDim = 0; iDim < nDim; iDim++) {
						ForceViscous[iDim] += TauElem[iDim]*Area*factor;
						//						ForceViscous[iDim] += WallShearStress*(VelTang[iDim]/VelTangMod)*Area*factor;
					}
                    
					if (iDim == 3) {
                        MomentViscous[0] += (TauElem[2]*MomentDist[1] - TauElem[1]*MomentDist[2])*Area*factor/RefLengthMoment;
                        MomentViscous[1] += (TauElem[0]*MomentDist[2] - TauElem[2]*MomentDist[0])*Area*factor/RefLengthMoment;
                    }
					MomentViscous[2] += (TauElem[1]*MomentDist[0] - TauElem[0]*MomentDist[1])*Area*factor/RefLengthMoment;
                    
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
					CEff_Visc[iMarker] = CLift_Visc[iMarker]/(CDrag_Visc[iMarker]+EPS);
					CFx_Visc[iMarker] = ForceViscous[0];
					CFy_Visc[iMarker] = ForceViscous[1];
					CFz_Visc[iMarker] = 0.0;
					CT_Visc[iMarker] = -CFx_Visc[iMarker];
					CQ_Visc[iMarker] = -CMz_Visc[iMarker];
                    Q_Visc[iMarker]  = HeatLoad;
					CMerit_Visc[iMarker] = 0.0;
				}
				if (nDim == 3) {
					CDrag_Visc[iMarker] =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
					CLift_Visc[iMarker] = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
					CMx_Visc[iMarker] = MomentViscous[0];
					CMy_Visc[iMarker] = MomentViscous[1];
					CMz_Visc[iMarker] = MomentViscous[2];
					CEff_Visc[iMarker] = CLift_Visc[iMarker]/(CDrag_Visc[iMarker]+EPS);
					CFx_Visc[iMarker] = ForceViscous[0];
					CFy_Visc[iMarker] = ForceViscous[1];
					CFz_Visc[iMarker] = ForceViscous[2];
					CT_Visc[iMarker] = -CFz_Visc[iMarker];
					CQ_Visc[iMarker] = -CMz_Visc[iMarker];
                    Q_Visc[iMarker] = HeatLoad;
					CMerit_Visc[iMarker] = CT_Visc[iMarker]/CQ_Visc[iMarker];
				}
                
				AllBound_CDrag_Visc += CDrag_Visc[iMarker];
				AllBound_CLift_Visc += CLift_Visc[iMarker];
				AllBound_CMx_Visc += CMx_Visc[iMarker];
				AllBound_CMy_Visc += CMy_Visc[iMarker];
				AllBound_CMz_Visc += CMz_Visc[iMarker];
				AllBound_CEff_Visc += CEff_Visc[iMarker];
				AllBound_CFx_Visc += CFx_Visc[iMarker];
				AllBound_CFy_Visc += CFy_Visc[iMarker];
				AllBound_CFz_Visc += CFz_Visc[iMarker];
				AllBound_CT_Visc += CT_Visc[iMarker];
				AllBound_CQ_Visc += CQ_Visc[iMarker];
                AllBound_Q_Visc += Q_Visc[iMarker];
                if (Maxq_Visc[iMarker] > AllBound_Maxq_Visc)
                    AllBound_Maxq_Visc = Maxq_Visc[iMarker];
				AllBound_CMerit_Visc += CMerit_Visc[iMarker];
			}
		}
	}
	Total_CDrag += AllBound_CDrag_Visc;
	Total_CLift += AllBound_CLift_Visc;
	Total_CMx += AllBound_CMx_Visc;
	Total_CMy += AllBound_CMy_Visc;
	Total_CMz += AllBound_CMz_Visc;
	Total_CEff = Total_CLift/(Total_CDrag+EPS);
	Total_CFx += AllBound_CFx_Visc;
	Total_CFy += AllBound_CFy_Visc;
	Total_CFz += AllBound_CFz_Visc;
	Total_CT += AllBound_CT_Visc;
	Total_CQ += AllBound_CQ_Visc;
    Total_Q += AllBound_Q_Visc;
    Total_Maxq = AllBound_Maxq_Visc;
	Total_CMerit += AllBound_CMerit_Visc;
    
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Tau[iDim];
	delete [] Tau;
	delete [] UnitaryNormal;
	delete [] TauTangent;
	delete [] TauElem;
}

void CNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
    
	/*--- Local variables ---*/
	unsigned short iDim, iVar;
	unsigned long iVertex, iPoint, total_index;
    
	double Wall_HeatFlux;
	double *Grid_Vel, *Rot_Vel, *Normal, Area;
    
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
	bool grid_movement  = config->GetGrid_Movement();
	bool rotating_frame = config->GetRotating_Frame();
    
	/*--- Identify the boundary by string name ---*/
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
    
	/*--- Get the specified wall heat flux from config ---*/
	Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);
    
	/*--- Loop over all of the vertices on this boundary marker ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
            
			/*--- Compute dual-grid area and boundary normal ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
            
			/*--- Initialize the convective & viscous residuals to zero ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Res_Conv[iVar] = 0.0;
				Res_Visc[iVar] = 0.0;
			}
            
			/*--- Store the corrected velocity at the wall which will
             be zero (v = 0), unless there are moving walls (v = u_wall)---*/
			if (rotating_frame) {
				Rot_Vel = geometry->node[iPoint]->GetRotVel();
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = Rot_Vel[iDim];
			} else if (grid_movement) {
				Grid_Vel = geometry->node[iPoint]->GetGridVel();
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = Grid_Vel[iDim];
			} else {
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
			}
            
			/*--- Set the residual, truncation error, and velocity value ---*/
			node[iPoint]->SetVelocity_Old(Vector, incompressible);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      
			node[iPoint]->SetVel_ResTruncError_Zero();
            
			/*--- Set the residual on the boundary with the specified heat flux ---*/
			Res_Visc[nDim+1] = Wall_HeatFlux * Area;
            
			/*--- Flux contribution due to a rotating frame (energy equation) ---*/
			if (rotating_frame) {
				double ProjRotVel = -geometry->vertex[val_marker][iVertex]->GetRotFlux();
				Res_Conv[nDim+1] = node[iPoint]->GetPressure(incompressible)*ProjRotVel;
			}
            
			/*--- Flux contribution due to grid motion (energy equation) ---*/
			if (grid_movement) {
				double ProjGridVel = 0.0;
				Grid_Vel = geometry->node[iPoint]->GetGridVel();
				for (iDim = 0; iDim < nDim; iDim++)
					ProjGridVel += Grid_Vel[iDim]*(-1.0)*Normal[iDim];
				Res_Conv[nDim+1] = node[iPoint]->GetPressure(incompressible)*ProjGridVel;
			}
            
			/*--- Viscous contribution for moving walls ---*/
			if (rotating_frame || grid_movement) {
                
				unsigned short jDim;
				double total_viscosity, div_vel, Density, val_turb_ke;
				Density = node[iPoint]->GetSolution(0);
				double val_laminar_viscosity = node[iPoint]->GetLaminarViscosity();
				double val_eddy_viscosity = node[iPoint]->GetEddyViscosity();
				double **val_gradprimvar = node[iPoint]->GetGradient_Primitive();
				double tau[3][3], delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
				total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
                
				/*--- Get the appropriate grid velocity at this node ---*/
				if (rotating_frame) {
					Grid_Vel = geometry->node[iPoint]->GetRotVel();
				} else if (grid_movement) {
					Grid_Vel = geometry->node[iPoint]->GetGridVel();
				}
                
				double Flux_Tensor[nVar][nDim];
				for (iVar = 0 ; iVar < nVar; iVar++)
					for (jDim = 0 ; jDim < nDim; jDim++)
						Flux_Tensor[iVar][jDim] = 0.0;
                
				/*--- Turbulent kinetic energy ---*/
				if (config->GetKind_Turb_Model() == SST)
					val_turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
				else
					val_turb_ke = 0.0;
                
				div_vel = 0.0;
				for (iDim = 0 ; iDim < nDim; iDim++)
					div_vel += val_gradprimvar[iDim+1][iDim];
                
				for (iDim = 0 ; iDim < nDim; iDim++)
					for (jDim = 0 ; jDim < nDim; jDim++)
						tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
						- TWO3*total_viscosity*div_vel*delta[iDim][jDim]
                        - TWO3*Density*val_turb_ke*delta[iDim][jDim];
                
				if (nDim == 2) {
					Flux_Tensor[0][0] = 0.0;
					Flux_Tensor[1][0] = 0.0;
					Flux_Tensor[2][0] = 0.0;
					Flux_Tensor[3][0] = tau[0][0]*Grid_Vel[0] + tau[0][1]*Grid_Vel[1];
                    
					Flux_Tensor[0][1] = 0.0;
					Flux_Tensor[1][1] = 0.0;
					Flux_Tensor[2][1] = 0.0;
					Flux_Tensor[3][1] = tau[1][0]*Grid_Vel[0] + tau[1][1]*Grid_Vel[1];
                    
				} else {
					Flux_Tensor[0][0] = 0.0;
					Flux_Tensor[1][0] = 0.0;
					Flux_Tensor[2][0] = 0.0;
					Flux_Tensor[3][0] = 0.0;
					Flux_Tensor[4][0] = tau[0][0]*Grid_Vel[0] + tau[0][1]*Grid_Vel[1] + tau[0][2]*Grid_Vel[2];
                    
					Flux_Tensor[0][1] = 0.0;
					Flux_Tensor[1][1] = 0.0;
					Flux_Tensor[2][1] = 0.0;
					Flux_Tensor[3][1] = 0.0;
					Flux_Tensor[4][1] = tau[1][0]*Grid_Vel[0] + tau[1][1]*Grid_Vel[1] + tau[1][2]*Grid_Vel[2];
                    
					Flux_Tensor[0][2] = 0.0;
					Flux_Tensor[1][2] = 0.0;
					Flux_Tensor[2][2] = 0.0;
					Flux_Tensor[3][2] = 0.0;
					Flux_Tensor[4][2] = tau[2][0]*Grid_Vel[0] + tau[2][1]*Grid_Vel[1] + tau[2][2]*Grid_Vel[2];
				}
                
				for (iVar = 0; iVar < nVar; iVar++) {
					for (iDim = 0; iDim < nDim; iDim++)
						Res_Visc[iVar] += Flux_Tensor[iVar][iDim] * -Normal[iDim];
				}
			}
            
			/*--- Convective contribution to the residual at the wall ---*/
			LinSysRes.AddBlock(iPoint, Res_Conv);
            
			/*--- Viscous contribution to the residual at the wall ---*/
			LinSysRes.SubtractBlock(iPoint, Res_Visc);
            
			/*--- Only change velocity-rows of the Jacobian (includes 1 in the diagonal)/
             Note that we need to add a contribution for moving walls to the Jacobian. ---*/
			if (implicit) {
				/*--- Enforce the no-slip boundary condition in a strong way ---*/
				for (iVar = 1; iVar <= nDim; iVar++) {
					total_index = iPoint*nVar+iVar;
					Jacobian.DeleteValsRowi(total_index);
				}
			}
            
		}
        
	}
}

void CNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
    
	unsigned long iVertex, iPoint, Point_Normal, total_index;
	unsigned short iVar, jVar, iDim;
	double *Normal, *Coord_i, *Coord_j, Area, dist_ij, *Grid_Vel, *Rot_Vel;
	double UnitaryNormal[3];
	double Twall, Temperature, dTdn, dTdrho;
	double Density, Vel2, Energy;
	double Laminar_Viscosity, Eddy_Viscosity, Thermal_Conductivity, Gas_Constant, cp;
	double Theta;
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();
	bool grid_movement = config->GetGrid_Movement();
	bool rotating_frame = config->GetRotating_Frame();
    
	Point_Normal = 0;
    
	/*--- Identify the boundary ---*/
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
    
	/*--- Retrieve the specified wall temperature ---*/
	Twall = config->GetIsothermal_Temperature(Marker_Tag);
	Gas_Constant = config->GetGas_ConstantND();
	cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
    
	/*--- Loop over boundary points ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if (geometry->node[iPoint]->GetDomain()) {
            
			/*--- Compute dual-grid area and boundary normal ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Area = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;
            
			/*--- Calculate useful quantities ---*/
			Theta = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Theta += UnitaryNormal[iDim]*UnitaryNormal[iDim];
            
			/*--- Compute closest normal neighbor ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
            
			/*--- Get cartesian coordinates of i & j and compute inter-node distance ---*/
			Coord_i = geometry->node[iPoint]->GetCoord();
			Coord_j = geometry->node[Point_Normal]->GetCoord();
			dist_ij = 0;
			for (iDim = 0; iDim < nDim; iDim++)
				dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
			dist_ij = sqrt(dist_ij);
            
            
			/*--- Store the corrected velocity at the wall which will
             be zero (v = 0), unless there is grid motion (v = u_wall)---*/
			if (rotating_frame) {
				Rot_Vel = geometry->node[iPoint]->GetRotVel();
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = Rot_Vel[iDim];
			} else if (grid_movement) {
				Grid_Vel = geometry->node[iPoint]->GetGridVel();
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = Grid_Vel[iDim];
			} else {
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
			}
            
			/*--- Initialize viscous residual (and Jacobian if implicit) to zero ---*/
			for (iVar = 0; iVar < nVar; iVar ++)
				Res_Visc[iVar] = 0.0;
			if (implicit) {
				for (iVar = 0; iVar < nVar; iVar ++)
					for (jVar = 0; jVar < nVar; jVar ++)
						Jacobian_i[iVar][jVar] = 0.0;
			}
            
			/*--- Set the residual, truncation error and velocity value on the boundary ---*/
			node[iPoint]->SetVelocity_Old(Vector, incompressible);
			for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
			node[iPoint]->SetVel_ResTruncError_Zero();
            
			/*--- Retrieve temperatures from boundary nearest neighbor --*/
			Temperature = node[Point_Normal]->GetPrimVar(0);
            
			/*--- Calculate temperature gradient normal to the wall using FD ---*/
			dTdn                  = (Twall - Temperature)/dist_ij;
			Laminar_Viscosity     = node[iPoint]->GetLaminarViscosity();
			Eddy_Viscosity        = node[iPoint]->GetEddyViscosity();
			Thermal_Conductivity  = cp * (Laminar_Viscosity/PRANDTL + Eddy_Viscosity/PRANDTL_TURB);
            
			/*--- Set the residual on the boundary, enforcing the no-slip boundary condition ---*/
			Res_Visc[nDim+1] = Thermal_Conductivity * dTdn * Area;
            
			/*--- Calculate Jacobian for implicit time stepping ---*/
			// Note: Momentum equations are enforced in a strong way while the energy equation is enforced weakly
			if (implicit) {
                
				/*--- Calculate useful quantities ---*/
				Density = node[iPoint]->GetPrimVar(nDim+2);
				Energy  = node[iPoint]->GetSolution(nDim+1);
				Temperature = node[iPoint]->GetPrimVar(0);
				Vel2 = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					Vel2 += node[iPoint]->GetPrimVar(iDim+1) * node[iPoint]->GetPrimVar(iDim+1);
				//dTdrho = (Gamma-1.0)/(Density*Gas_Constant) * (-Energy/Density + Vel2);
				dTdrho = 1.0/Density * ( -Twall + (Gamma-1.0)/Gas_Constant*(Vel2/2.0) );
                
				/*--- Enforce the no-slip boundary condition in a strong way ---*/
				for (iVar = 1; iVar <= nDim; iVar++) {
					total_index = iPoint*nVar+iVar;
					Jacobian.DeleteValsRowi(total_index);
				}
                
				/*--- Add contributions to the Jacobian from the weak enforcement of the energy equations ---*/
				Jacobian_i[nDim+1][0]      = -Thermal_Conductivity*Theta/dist_ij * dTdrho * Area;
				Jacobian_i[nDim+1][nDim+1] = -Thermal_Conductivity*Theta/dist_ij * (Gamma-1.0)/(Gas_Constant*Density) * Area;
			}
            
			/*--- Apply calculated residuals and Jacobians to the linear system ---*/
			LinSysRes.SubtractBlock(iPoint, Res_Visc);
			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
		}
	}
}
