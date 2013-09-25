/*!
 * \file solution_direct_tne2.cpp
 * \brief Main subrotuines for solving flows in thermochemical nonequilibrium.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.7
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
#include <math.h>

CTNE2EulerSolver::CTNE2EulerSolver(void) : CSolver() {
  
	/*--- Array initialization ---*/
	Velocity_Inlet  = NULL;
	Velocity_Outlet = NULL;
  Velocity_Inf    = NULL;
	CDrag_Inv       = NULL;
	CLift_Inv       = NULL;
	CSideForce_Inv  = NULL;
	CMx_Inv         = NULL;
	CMy_Inv         = NULL;
	CMz_Inv         = NULL;
	CFx_Inv         = NULL;
	CFy_Inv         = NULL;
	CFz_Inv         = NULL;
	CEff_Inv        = NULL;
	ForceInviscid   = NULL;
	MomentInviscid  = NULL;
	PrimVar_i       = NULL;
	PrimVar_j       = NULL;
	Precon_Mat_inv  = NULL;
	CPressure       = NULL;
	CHeatTransfer   = NULL;
  
}

CTNE2EulerSolver::CTNE2EulerSolver(CGeometry *geometry, CConfig *config,
                                   unsigned short iMesh) : CSolver() {

	unsigned long iPoint, index, counter_local = 0, counter_global = 0;
	unsigned short iVar, iDim, iMarker, iSpecies, nZone;
  double *Mvec_Inf;
  double Alpha, Beta, dull_val;
	bool restart, check_infty, check_node, check_temp, check_press;
  
  /*--- Get MPI rank ---*/
	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
	/*--- Array initialization ---*/
	Velocity_Inlet  = NULL;
	Velocity_Outlet = NULL;
  Velocity_Inf    = NULL;
	CDrag_Inv       = NULL;
	CLift_Inv       = NULL;
	CSideForce_Inv  = NULL;
	CMx_Inv         = NULL;
	CMy_Inv         = NULL;
	CMz_Inv         = NULL;
	CFx_Inv         = NULL;
	CFy_Inv         = NULL;
	CFz_Inv         = NULL;
	CEff_Inv        = NULL;
	ForceInviscid   = NULL;
	MomentInviscid  = NULL;
	PrimVar_i       = NULL;
	PrimVar_j       = NULL;
	Precon_Mat_inv  = NULL;
	CPressure       = NULL;
	CHeatTransfer   = NULL;
  
  /*--- Set booleans for solver settings ---*/
  restart    = (config->GetRestart() || config->GetRestart_Flow());
	roe_turkel = false;
  
  /*--- Define constants in the solver structure ---*/	
	nSpecies     = config->GetnSpecies();
  nMarker      = config->GetnMarker_All();
	nPoint       = geometry->GetnPoint();
	nPointDomain = geometry->GetnPointDomain();
  nZone        = geometry->GetnZone();
  nDim         = geometry->GetnDim();
  
  /*--- Set size of the conserved and primitive vectors ---*/
  //     U: [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
  //     V: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradV: [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  nVar         = nSpecies + nDim + 2;
  nPrimVar     = nSpecies + nDim + 8;
  nPrimVarGrad = nSpecies + nDim + 3;
  
	/*--- Allocate a CVariable array for each node of the mesh ---*/
	node = new CVariable*[nPoint];
  
	/*--- Define some auxiliary vectors related to the residual ---*/
	Residual      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
	Residual_RMS  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
	Residual_Max  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
	Point_Max     = new unsigned long[nVar]; for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
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
	Vector   = new double[nDim];
  for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
	Vector_i = new double[nDim];
  for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
	Vector_j = new double[nDim];
  for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
	if ((config->GetKind_Upwind_TNE2() == ROE_TURKEL_2ND) ||
      (config->GetKind_Upwind_TNE2() == ROE_TURKEL_1ST)   ) {
		Precon_Mat_inv = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar ++)
			Precon_Mat_inv[iVar] = new double[nVar];
		roe_turkel = true;
	}
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
	/*--- Allocate Jacobians for implicit time-stepping ---*/
	if (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT) {
		Jacobian_i = new double* [nVar];
		Jacobian_j = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new double [nVar];
			Jacobian_j[iVar] = new double [nVar];
		}
    
		/*--- Initialization of the structure for the global Jacobian ---*/
		if (rank == MASTER_NODE)
      cout << "Initialize Jacobian structure. MG level: "
           << iMesh <<"." << endl;
		Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
	} else {
		if (rank == MASTER_NODE)
			cout << "Explicit scheme. No jacobian structure (Euler). MG level: "
           << iMesh <<"." << endl;
	}
  
	/*--- Allocate arrays for gradient computation by least squares ---*/
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
    
		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nPrimVarGrad];
		for (iVar = 0; iVar < nPrimVarGrad; iVar++)
			cvector[iVar] = new double [nDim];
	}
  
	/*--- Allocate force & coefficient arrays on boundaries ---*/
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
  
	/*--- Initialize total coefficients ---*/
	Total_CDrag = 0.0;  Total_CLift = 0.0;  Total_CSideForce = 0.0;
	Total_CMx   = 0.0;  Total_CMy   = 0.0;  Total_CMz = 0.0;
	Total_CFx   = 0.0;  Total_CFy   = 0.0;  Total_CFz = 0.0;
  Total_CEff  = 0.0;
  Total_Maxq  = 0.0;
  
	/*--- Read farfield conditions from the config file ---*/
	Pressure_Inf       = config->GetPressure_FreeStream();
  Temperature_Inf    = config->GetTemperature_FreeStream();
  Temperature_ve_Inf = config->GetTemperature_ve_FreeStream();
  MassFrac_Inf       = config->GetMassFrac_FreeStream();
  Mach_Inf           = config->GetMach_FreeStreamND();
  
  /*--- Vectorize free stream Mach number based on AoA & AoS ---*/
  Mvec_Inf = new double[nDim];
  Alpha    = config->GetAoA();
  Beta     = config->GetAoS();
  if (nDim == 2) {
    Mvec_Inf[0] = cos(Alpha)*Mach_Inf;
    Mvec_Inf[1] = sin(Alpha)*Mach_Inf;
  }
  if (nDim == 3) {
    Mvec_Inf[0] = cos(Alpha)*cos(Beta)*Mach_Inf;
    Mvec_Inf[1] = sin(Beta)*Mach_Inf;
    Mvec_Inf[2] = sin(Alpha)*cos(Beta)*Mach_Inf;
  }  
  
  /*--- Create a CVariable that stores the free-stream values ---*/
  node_infty = new CTNE2EulerVariable(Pressure_Inf, MassFrac_Inf,
                                      Mvec_Inf, Temperature_Inf,
                                      Temperature_ve_Inf, nDim, nVar,
                                      nPrimVar, nPrimVarGrad, config);
  check_infty = node_infty->SetPrimVar_Compressible(config);
  
  Velocity_Inf = new double[nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_Inf[iDim] = node_infty->GetVelocity(iDim);
  
	/*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
  
	if (!restart || geometry->GetFinestMGLevel() == false || nZone > 1) {

		/*--- Initialize using freestream values ---*/
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CTNE2EulerVariable(Pressure_Inf, MassFrac_Inf,
                                            Mvec_Inf, Temperature_Inf,
                                            Temperature_ve_Inf, nDim,
                                            nVar, nPrimVar, nPrimVarGrad,
                                            config);
//      node[iPoint]->SetPrimVar_Compressible(config);
    }
	} else {
    
		/*--- Restart the solution from file information ---*/
		ifstream restart_file;
		string filename = config->GetSolution_FlowFileName();
    
		/*--- Append time step for unsteady restart ---*/
		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
			char buffer[50];
			unsigned long flowIter = config->GetnExtIter() - 1;
			filename.erase (filename.end()-4, filename.end());
			if ((int(flowIter) >= 0) && (int(flowIter) < 10))
        sprintf (buffer, "_0000%d.dat", int(flowIter));
			if ((int(flowIter) >= 10) && (int(flowIter) < 100))
        sprintf (buffer, "_000%d.dat", int(flowIter));
			if ((int(flowIter) >= 100) && (int(flowIter) < 1000))
        sprintf (buffer, "_00%d.dat", int(flowIter));
			if ((int(flowIter) >= 1000) && (int(flowIter) < 10000))
        sprintf (buffer, "_0%d.dat", int(flowIter));
			if (int(flowIter) >= 10000)
        sprintf (buffer, "_%d.dat", int(flowIter));
			string UnstExt = string(buffer);
			filename.append(UnstExt);
		}
		restart_file.open(filename.data(), ios::in);
    
		/*--- In case there is no restart file ---*/
		if (restart_file.fail()) {
			cout << "There is no flow restart file!! " << filename.data()
           << "."<< endl;
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
/*        point_line >> index;
        for (iDim = 0; iDim < nDim; iDim++)
          point_line >> dull_val;
        for (iVar = 0; iVar < nVar; iVar++) {
          point_line >> Solution[iVar];
        }*/
        
        point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4] >> Solution[5] >> Solution[6];
        
				node[iPoint_Local] = new CTNE2EulerVariable(Solution, nDim, nVar, nPrimVar,
                                                    nPrimVarGrad, config);
			}
			iPoint_Global++;
		}
    
		/*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    //Solution = node_infty->GetSolution();
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++)
      node[iPoint] = new CTNE2EulerVariable(Pressure_Inf, MassFrac_Inf,
                                            Mvec_Inf, Temperature_Inf,
                                            Temperature_ve_Inf, nDim,
                                            nVar, nPrimVar, nPrimVarGrad,
                                            config);
			//node[iPoint] = new CTNE2EulerVariable(Solution, nDim, nVar, nPrimVar, nPrimVarGrad, config);
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
	}
  
  /*--- Check that the initial solution is physical ---*/
  counter_local = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
 
//    check_node = node[iPoint]->SetPrimVar_Compressible(config);
    
    
    node[iPoint]->SetDensity();
    node[iPoint]->SetVelocity2();
    check_temp = node[iPoint]->SetTemperature(config);
    check_press = node[iPoint]->SetPressure(config);
    
    if (check_temp || check_press) {
      bool ionization;
      unsigned short iEl, nHeavy, nEl, *nElStates;
      double Ru, T, Tve, rhoCvtr, sqvel, rhoE, rhoEve, num, denom, conc;
      double rho, rhos, Ef, Ev, Ee, soundspeed;
      double *xi, *Ms, *thetav, **thetae, **g, *Tref, *hf;
      /*--- Determine the number of heavy species ---*/
      ionization = config->GetIonization();
      if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
      else            { nHeavy = nSpecies;   nEl = 0; }
      
      /*--- Load variables from the config class --*/
      xi        = config->GetRotationModes();      // Rotational modes of energy storage
      Ms        = config->GetMolar_Mass();         // Species molar mass
      thetav    = config->GetCharVibTemp();        // Species characteristic vib. temperature [K]
      thetae    = config->GetCharElTemp();         // Characteristic electron temperature [K]
      g         = config->GetElDegeneracy();       // Degeneracy of electron states
      nElStates = config->GetnElStates();          // Number of electron states
      Tref      = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
      hf        = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]
      
      /*--- Rename & initialize for convenience ---*/
      Ru      = UNIVERSAL_GAS_CONSTANT;         // Universal gas constant [J/(kmol*K)]
      Tve     = Temperature_ve_Inf;             // Vibrational temperature [K]
      T       = Temperature_Inf;                // Translational-rotational temperature [K]
      sqvel   = 0.0;                            // Velocity^2 [m2/s2]
      rhoE    = 0.0;                            // Mixture total energy per mass [J/kg]
      rhoEve  = 0.0;                            // Mixture vib-el energy per mass [J/kg]
      denom   = 0.0;
      conc    = 0.0;
      rhoCvtr = 0.0;
      
      /*--- Calculate mixture density from supplied primitive quantities ---*/
      for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
        denom += MassFrac_Inf[iSpecies] * (Ru/Ms[iSpecies]) * T;
      for (iSpecies = 0; iSpecies < nEl; iSpecies++)
        denom += MassFrac_Inf[nSpecies-1] * (Ru/Ms[nSpecies-1]) * Tve;
      rho = Pressure_Inf / denom;
      
      /*--- Calculate sound speed and extract velocities ---*/
      for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
        conc += MassFrac_Inf[iSpecies]*rho/Ms[iSpecies];
        rhoCvtr += rho*MassFrac_Inf[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
      }
      soundspeed = sqrt((1.0 + Ru/rhoCvtr*conc) * Pressure_Inf/rho);
      for (iDim = 0; iDim < nDim; iDim++)
        sqvel += Mvec_Inf[iDim]*soundspeed * Mvec_Inf[iDim]*soundspeed;
      
      /*--- Calculate energy (RRHO) from supplied primitive quanitites ---*/
      for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
        // Species density
        rhos = MassFrac_Inf[iSpecies]*rho;
        
        // Species formation energy
        Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];
        
        // Species vibrational energy
        if (thetav[iSpecies] != 0.0)
          Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Tve)-1.0);
        else
          Ev = 0.0;
        
        // Species electronic energy
        num = 0.0;
        denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);
        for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
          num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
          denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
        }
        Ee = Ru/Ms[iSpecies] * (num/denom);
        
        // Mixture total energy
        rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (T-Tref[iSpecies])
                        + Ev + Ee + Ef + 0.5*sqvel);
        
        // Mixture vibrational-electronic energy
        rhoEve += rhos * (Ev + Ee);
      }
      for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
        // Species formation energy
        Ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1] * Tref[nSpecies-1];
        
        // Electron t-r mode contributes to mixture vib-el energy
        rhoEve += (3.0/2.0) * Ru/Ms[nSpecies-1] * (Tve - Tref[nSpecies-1]);
      }
      
      /*--- Initialize Solution & Solution_Old vectors ---*/
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        Solution[iSpecies]     = rho*MassFrac_Inf[iSpecies];
      }
      for (iDim = 0; iDim < nDim; iDim++) {
        Solution[nSpecies+iDim]     = rho*Mvec_Inf[iDim]*soundspeed;
      }
      Solution[nSpecies+nDim]       = rhoE;
      Solution[nSpecies+nDim+1]     = rhoEve;
      
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
    cout << "Warning. The original solution contains "<< counter_global
         << " points that are not physical." << endl;

	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
	else least_squares = false;
  
	/*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);
  
  /*--- Deallocate arrays ---*/
  delete [] Mvec_Inf;
}

CTNE2EulerSolver::~CTNE2EulerSolver(void) {
	unsigned short iVar, iMarker;
  
	/*--- Array initialization ---*/
	if (Velocity_Inlet != NULL) delete [] Velocity_Inlet;
	if (Velocity_Outlet != NULL) delete [] Velocity_Outlet;
  if (Velocity_Inf != NULL) delete [] Velocity_Inf;
	if (CDrag_Inv != NULL) delete [] CDrag_Inv;
	if (CLift_Inv != NULL) delete [] CLift_Inv;
	if (CSideForce_Inv != NULL) delete [] CSideForce_Inv;
	if (CMx_Inv != NULL) delete [] CMx_Inv;
	if (CMy_Inv != NULL) delete [] CMy_Inv;
	if (CMz_Inv != NULL) delete [] CMz_Inv;
	if (CFx_Inv != NULL) delete [] CFx_Inv;
	if (CFy_Inv != NULL) delete [] CFy_Inv;
	if (CFz_Inv != NULL) delete [] CFz_Inv;
	if (CEff_Inv != NULL) delete [] CEff_Inv;
	if (ForceInviscid != NULL) delete [] ForceInviscid;
	if (MomentInviscid != NULL) delete [] MomentInviscid;
	if (PrimVar_i != NULL) delete [] PrimVar_i;
	if (PrimVar_j != NULL) delete [] PrimVar_j;
  
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
}

void CTNE2EulerSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
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
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        
        rotMatrix[0][1] = cosPhi*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        
        rotMatrix[0][2] = -sinPhi;
        rotMatrix[1][2] = sinTheta*cosPhi;
        rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex];
        }
        else {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+2] = rotMatrix[2][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
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

void CTNE2EulerSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
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
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex];
        }
        else {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+2] = rotMatrix[2][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
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

void CTNE2EulerSolver::Set_MPI_Primitive(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR, VEL_INDEX;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_V = NULL, *Buffer_Send_V = NULL;
	int send_to, receive_from;
  double *Primitive;
  
  Primitive = new double[nPrimVar];
  VEL_INDEX = node_infty->GetVelIndex();
  
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
			nBufferS_Vector = nVertexS*nPrimVar;    nBufferR_Vector = nVertexR*nPrimVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_V = new double [nBufferR_Vector];
      Buffer_Send_V = new double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVar; iVar++)
          Buffer_Send_V[iVar*nVertexS+iVertex] = node[iPoint]->GetPrimVar(iVar);
      }
      
#ifndef NO_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_V, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_V, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVar; iVar++)
          Buffer_Receive_V[iVar*nVertexR+iVertex] = Buffer_Send_V[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_V;
      
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
        for (iVar = 0; iVar < nPrimVar; iVar++)
          Primitive[iVar] = Buffer_Receive_V[iVar*nVertexR+iVertex];
        
        /*--- Rotate the velocity components. ---*/
        if (nDim == 2) {
          Primitive[VEL_INDEX]   = rotMatrix[0][0]*Buffer_Receive_V[(VEL_INDEX)*nVertexR+iVertex]
                                 + rotMatrix[0][1]*Buffer_Receive_V[(VEL_INDEX+1)*nVertexR+iVertex];
          Primitive[VEL_INDEX+1] = rotMatrix[1][0]*Buffer_Receive_V[(VEL_INDEX)*nVertexR+iVertex]
                                 + rotMatrix[1][1]*Buffer_Receive_V[(VEL_INDEX+1)*nVertexR+iVertex];
        }
        else {
          Primitive[VEL_INDEX]   = rotMatrix[0][0]*Buffer_Receive_V[(VEL_INDEX)*nVertexR+iVertex]
                                 + rotMatrix[0][1]*Buffer_Receive_V[(VEL_INDEX+1)*nVertexR+iVertex]
                                 + rotMatrix[0][2]*Buffer_Receive_V[(VEL_INDEX+2)*nVertexR+iVertex];
          Primitive[VEL_INDEX+1] = rotMatrix[1][0]*Buffer_Receive_V[(VEL_INDEX)*nVertexR+iVertex]
                                 + rotMatrix[1][1]*Buffer_Receive_V[(VEL_INDEX+1)*nVertexR+iVertex]
                                 + rotMatrix[1][2]*Buffer_Receive_V[(VEL_INDEX+2)*nVertexR+iVertex];
          Primitive[VEL_INDEX+2] = rotMatrix[2][0]*Buffer_Receive_V[(VEL_INDEX)*nVertexR+iVertex]
                                 + rotMatrix[2][1]*Buffer_Receive_V[(VEL_INDEX+1)*nVertexR+iVertex]
                                 + rotMatrix[2][2]*Buffer_Receive_V[(VEL_INDEX+2)*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nPrimVar; iVar++)
          node[iPoint]->SetPrimVar(iVar, Primitive[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_V;
    }
	}
  delete [] Primitive;
}


void CTNE2EulerSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
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
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex];
        }
        else {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_Limit[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_Limit[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+2] = rotMatrix[2][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_Limit[(nSpecies+2)*nVertexR+iVertex];
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

void CTNE2EulerSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry,
                                                   CConfig *config) {
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
      MPI::COMM_WORLD.Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector,
                               MPI::DOUBLE, send_to, 0,
                               Buffer_Receive_Undivided_Laplacian,
                               nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
      
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
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex];
        }
        else {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_Undivided_Laplacian[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_Undivided_Laplacian[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+2] = rotMatrix[2][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_Undivided_Laplacian[(nSpecies+2)*nVertexR+iVertex];
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

void CTNE2EulerSolver::Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config) {
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

void CTNE2EulerSolver::Set_MPI_Dissipation_Switch(CGeometry *geometry, CConfig *config) {
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

void CTNE2EulerSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
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
      
      /*--- Copy the solution old that should be sent ---*/
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
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
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

void CTNE2EulerSolver::Set_MPI_PrimVar_Gradient(CGeometry *geometry, CConfig *config) {
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
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
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



void CTNE2EulerSolver::Preprocessing(CGeometry *geometry,
                                     CSolver **solution_container,
                                     CConfig *config, unsigned short iMesh,
                                     unsigned short iRKStep,
                                     unsigned short RunTime_EqSystem) {
#ifdef NO_MPI
	int rank = MASTER_NODE;
#else
	int rank = MPI::COMM_WORLD.Get_rank();
#endif
	
  unsigned long iPoint, ErrorCounter = 0;
	bool implicit   = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool center     = (config->GetKind_ConvNumScheme_TNE2() == SPACE_CENTERED);
	bool upwind_2nd = ((config->GetKind_Upwind_TNE2() == ROE_2ND)  ||
                     (config->GetKind_Upwind_TNE2() == AUSM_2ND) ||
                     (config->GetKind_Upwind_TNE2() == HLLC_2ND) ||
                     (config->GetKind_Upwind_TNE2() == ROE_TURKEL_2ND));
	bool limiter = (config->GetKind_SlopeLimit_TNE2() != NONE);
  bool RightSol;
  
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {

		/*--- Primitive variables [rho1,...,rhoNs,T,Tve,u,v,w,P,rho,h,c] ---*/
		RightSol = node[iPoint]->SetPrimVar_Compressible(config);
    if (!RightSol) ErrorCounter++;
    
    /*--- Initialize the convective residual vector ---*/
		LinSysRes.SetBlock_Zero(iPoint);
    
	}
  
  Set_MPI_Primitive(geometry, config);
  
	/*--- Upwind second order reconstruction ---*/
	if ((upwind_2nd) && (iMesh == MESH_0)) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    
		/*--- Limiter computation ---*/
		if ((limiter) && (iMesh == MESH_0)) SetSolution_Limiter(geometry, config);
	}
  
  /*--- Artificial dissipation ---*/
  if (center)
    SetMax_Eigenvalue(geometry, config);
  
	/*--- Initialize the jacobian matrices ---*/
	if (implicit) Jacobian.SetValZero();
  
  /*--- Error message ---*/
#ifndef NO_MPI
  unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
  MPI::COMM_WORLD.Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI::UNSIGNED_LONG, MPI::SUM);
#endif
  if ((ErrorCounter != 0) && (rank == MASTER_NODE))
    cout <<"The solution contains "<< ErrorCounter << " non-physical points." << endl;
  
}

void CTNE2EulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solution_container, CConfig *config,
                                      unsigned short iMesh, unsigned long Iteration) {
	double *Normal, Area, Vol, Mean_SoundSpeed, Mean_ProjVel, Lambda, Local_Delta_Time,
	Global_Delta_Time = 1E6, Global_Delta_UnstTimeND;
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker;
  
	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  
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
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    
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
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
      
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
}


void CTNE2EulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {
  
	double *Normal, Area, Vol, Mean_SoundSpeed, Mean_ProjVel, Lambda, Local_Delta_Time,
	Global_Delta_Time = 1E6, Global_Delta_UnstTimeND;
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker;
  
	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  
	Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
  
	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		node[iPoint]->SetLambda(0.0);
  
	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Point identification, Normal vector and area ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
    
		Normal = geometry->edge[iEdge]->GetNormal();
		Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
		/*--- Mean Values ---*/
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    
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
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
      
			/*--- Inviscid contribution ---*/
			Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
			if (geometry->node[iPoint]->GetDomain()) {
				node[iPoint]->AddLambda(Lambda);
			}
		}
	}
  
  /*--- Call the MPI routine ---*/
  Set_MPI_MaxEigenvalue(geometry, config);

}



void CTNE2EulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint;
  bool implicit, high_order_diss;

  /*--- Set booleans based on config settings ---*/
	implicit        = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	high_order_diss = ((config->GetKind_Centered_TNE2() == JST) && (iMesh == MESH_0));
  
  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( node[0]->GetPIndex()       );
  numerics->SetTIndex      ( node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( node[0]->GetHIndex()       );
  numerics->SetAIndex      ( node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( node[0]->GetRhoCvveIndex() );
  
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Points in edge, set normal vectors, and number of neighbors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
		numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());
    
		/*--- Pass conservative & primitive variables w/o reconstruction to CNumerics ---*/
		numerics->SetConservative(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    numerics->SetPrimitive(node[iPoint]->GetPrimVar(), node[jPoint]->GetPrimVar());
    
    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdrhos(node[iPoint]->GetdPdrhos(), node[jPoint]->GetdPdrhos());
    
    /*--- Set the largest convective eigenvalue ---*/
		numerics->SetLambda(node[iPoint]->GetLambda(), node[jPoint]->GetLambda());
    
		/*--- Compute residuals, and Jacobians ---*/
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


void CTNE2EulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solution_container, CNumerics *numerics,
                                       CConfig *config, unsigned short iMesh) {
	unsigned long iEdge, iPoint, jPoint;
  bool implicit, high_order_diss, limiter;
  double *U_i, *U_j, *V_i, *V_j;
  
	/*--- Set booleans based on config settings ---*/
	implicit        = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  limiter         = (config->GetKind_SlopeLimit_TNE2() != NONE);
	high_order_diss = (((config->GetKind_Upwind_TNE2() == ROE_2ND)  ||
                      (config->GetKind_Upwind_TNE2() == AUSM_2ND) ||
                      (config->GetKind_Upwind_TNE2() == HLLC_2ND) ||
                      (config->GetKind_Upwind_TNE2() == ROE_TURKEL_2ND))
                     && (iMesh == MESH_0));
  
  if (high_order_diss) cout << "WARNING!!! Upwind_Residual: 2nd order accuracy not in place!" << endl;
  if (limiter) cout << "WARNING!!! Upwind_Residual: Limiter not in place!" << endl;
  
  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( node[0]->GetPIndex()       );
  numerics->SetTIndex      ( node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( node[0]->GetHIndex()       );
  numerics->SetAIndex      ( node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( node[0]->GetRhoCvveIndex() );
  
  /*--- Loop over edges and calculate convective fluxes ---*/
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Retrieve node numbers and pass edge normal to CNumerics ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Pass conserved and primitive variables from CVariable to CNumerics class ---*/
    U_i = node[iPoint]->GetSolution(); U_j = node[jPoint]->GetSolution();
    V_i = node[iPoint]->GetPrimVar();  V_j = node[jPoint]->GetPrimVar();
    numerics->SetPrimitive(V_i, V_j);
    numerics->SetConservative(U_i, U_j);
    
    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdrhos(node[iPoint]->GetdPdrhos(), node[jPoint]->GetdPdrhos());
    
    /*--- Compute the upwind residual ---*/
		numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);    
    
//    unsigned short iVar, jVar;
//    cout << endl << endl << "Upw Jacobian: " << endl;
//    for (iVar =0; iVar < nVar; iVar++) {
//      for (jVar = 0; jVar < nVar; jVar++) {
//        cout << Jacobian_i[iVar][jVar] << "\t";
//      }
//      cout << endl;
//    }
//    cin.get();
		
    /*--- Update the residual values ---*/
		LinSysRes.AddBlock(iPoint, Res_Conv);
		LinSysRes.SubtractBlock(jPoint, Res_Conv);
    
		/*--- Update the implicit Jacobian ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
			Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
			Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
		}
	}
}

void CTNE2EulerSolver::Source_Residual(CGeometry *geometry, CSolver **solution_container, CNumerics *numerics,
                                       CNumerics *second_solver, CConfig *config, unsigned short iMesh) {
  bool implicit;
	unsigned short iVar, jVar;
	unsigned long iPoint;
  
  /*--- Assign booleans ---*/
  implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  
  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( node[0]->GetPIndex()       );
  numerics->SetTIndex      ( node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( node[0]->GetHIndex()       );
  numerics->SetAIndex      ( node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( node[0]->GetRhoCvveIndex() );
  
/*  double delta;
  double *U_i, *V_i;
  double *Res_new, **FDJac;
  Res_new = new double[nVar];
  FDJac = new double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    FDJac[iVar] = new double[nVar];*/
  
  /*--- loop over points ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Initialize Residual & Jacobian arrays to zero ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Residual[iVar] = 0;
    if (implicit)
      for (iVar =0 ; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;
    
    /*--- Set conserved & primitive variables  ---*/
    numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
    numerics->SetPrimitive   (node[iPoint]->GetPrimVar(),  node[iPoint]->GetPrimVar() );

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdrhos(node[iPoint]->GetdPdrhos(), node[iPoint]->GetdPdrhos());
    
    /*--- Set volume of the dual grid cell ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());

    /*--- Compute the non-equilibrium chemistry ---*/
    numerics->ComputeChemistry(Residual, Jacobian_i, config);
    
    /*--- Compute vibrational energy relaxation ---*/
    // NOTE: Jacobians don't account for relaxation time derivatives
    numerics->ComputeVibRelaxation(Residual, Jacobian_i, config);
    
    /*--- Subtract Residual (and Jacobian) ---*/
    LinSysRes.SubtractBlock(iPoint, Residual);
    if (implicit)
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
  }
  // Deallocate
//  delete[] Res_new;
}

void CTNE2EulerSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) {
	unsigned long iVertex, iPoint;
	unsigned short iDim, iMarker, Boundary, Monitoring;
	double Pressure, *Normal = NULL, dist[3], *Coord, Face_Area, PressInviscid;
	double factor, NFPressOF, RefVel2, RefDensity, RefPressure;
  
	bool grid_movement  = config->GetGrid_Movement();
  
	double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
	double Beta            = config->GetAoS()*PI_NUMBER/180.0;
	double RefAreaCoeff    = config->GetRefAreaCoeff();
	double RefLengthMoment = config->GetRefLengthMoment();
	double *Origin         = config->GetRefOriginMoment();
	//double Pressure_Inf    = config->GetPressure_FreeStreamND();
	//double *Velocity_Inf   = config->GetVelocity_FreeStreamND();
  
	/*--- If we have a rotating frame problem or an unsteady problem with
   mesh motion, use special reference values for the force coefficients.
   Otherwise, use the freestream values, which is the standard convention. ---*/
  
  RefVel2     = node_infty->GetVelocity2();
	RefDensity  = node_infty->GetDensity();
	RefPressure = node_infty->GetPressure();
  
	/*-- Initialization ---*/
	Total_CDrag = 0.0; Total_CLift = 0.0;  Total_CSideForce = 0.0;
	Total_CMx = 0.0;   Total_CMy = 0.0;    Total_CMz = 0.0;
	Total_CFx = 0.0;   Total_CFy = 0.0;    Total_CFz = 0.0;
	Total_CEff = 0.0;  Total_Q = 0.0;
  Total_Maxq = 0.0;
	AllBound_CDrag_Inv = 0.0;        AllBound_CLift_Inv = 0.0;  AllBound_CSideForce_Inv = 0.0;
	AllBound_CMx_Inv = 0.0;          AllBound_CMy_Inv = 0.0;    AllBound_CMz_Inv = 0.0;
	AllBound_CFx_Inv = 0.0;          AllBound_CFy_Inv = 0.0;    AllBound_CFz_Inv = 0.0;
	AllBound_CEff_Inv = 0.0;
  
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
				Pressure = node[iPoint]->GetPressure();
        
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
						CFx_Inv[iMarker]          = ForceInviscid[0];
						CFy_Inv[iMarker]          = ForceInviscid[1];
						CFz_Inv[iMarker]          = 0.0;
					}
					else {
						CDrag_Inv[iMarker] = 0.0; CLift_Inv[iMarker] = 0.0; CSideForce_Inv[iMarker] = 0.0;
						CMx_Inv[iMarker] = 0.0; CMy_Inv[iMarker] = 0.0; CMz_Inv[iMarker] = 0.0;
						CFx_Inv[iMarker] = 0.0; CFy_Inv[iMarker] = 0.0; CFz_Inv[iMarker] = 0.0;
						CEff_Inv[iMarker] = 0.0;
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
						CFx_Inv[iMarker] = ForceInviscid[0];
						CFy_Inv[iMarker] = ForceInviscid[1];
						CFz_Inv[iMarker] = ForceInviscid[2];
            
					}
					else {
						CDrag_Inv[iMarker] = 0.0; CLift_Inv[iMarker] = 0.0; CSideForce_Inv[iMarker] = 0.0;
						CMx_Inv[iMarker] = 0.0; CMy_Inv[iMarker] = 0.0; CMz_Inv[iMarker] = 0.0;
						CFx_Inv[iMarker] = 0.0; CFy_Inv[iMarker] = 0.0; CFz_Inv[iMarker] = 0.0;
						CEff_Inv[iMarker] = 0.0;
					}
				}
        
				AllBound_CDrag_Inv += CDrag_Inv[iMarker];
				AllBound_CLift_Inv += CLift_Inv[iMarker];
				AllBound_CSideForce_Inv += CSideForce_Inv[iMarker];
				AllBound_CMx_Inv += CMx_Inv[iMarker];
				AllBound_CMy_Inv += CMy_Inv[iMarker];
				AllBound_CMz_Inv += CMz_Inv[iMarker];
				AllBound_CEff_Inv += CEff_Inv[iMarker];
				AllBound_CFx_Inv += CFx_Inv[iMarker];
				AllBound_CFy_Inv += CFy_Inv[iMarker];
				AllBound_CFz_Inv += CFz_Inv[iMarker];
			}
		}
	}
  
#ifndef NO_MPI
  /*--- Add AllBound information using all the nodes ---*/
  double MyAllBound_CDrag_Inv        = AllBound_CDrag_Inv;        AllBound_CDrag_Inv = 0.0;
	double MyAllBound_CLift_Inv        = AllBound_CLift_Inv;        AllBound_CLift_Inv = 0.0;
	double MyAllBound_CSideForce_Inv   = AllBound_CSideForce_Inv;   AllBound_CSideForce_Inv = 0.0;
	double MyAllBound_CEff_Inv         = AllBound_CEff_Inv;         AllBound_CEff_Inv = 0.0;
	double MyAllBound_CMx_Inv          = AllBound_CMx_Inv;          AllBound_CMx_Inv = 0.0;
	double MyAllBound_CMy_Inv          = AllBound_CMy_Inv;          AllBound_CMy_Inv = 0.0;
	double MyAllBound_CMz_Inv          = AllBound_CMz_Inv;          AllBound_CMz_Inv = 0.0;
	double MyAllBound_CFx_Inv          = AllBound_CFx_Inv;          AllBound_CFx_Inv = 0.0;
	double MyAllBound_CFy_Inv          = AllBound_CFy_Inv;          AllBound_CFy_Inv = 0.0;
	double MyAllBound_CFz_Inv          = AllBound_CFz_Inv;          AllBound_CFz_Inv = 0.0;
  
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CDrag_Inv, &AllBound_CDrag_Inv, 1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CLift_Inv, &AllBound_CLift_Inv, 1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CSideForce_Inv, &AllBound_CSideForce_Inv, 1, MPI::DOUBLE, MPI::SUM);
  AllBound_CEff_Inv = AllBound_CLift_Inv / (AllBound_CDrag_Inv + config->GetCteViscDrag() + EPS);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CMx_Inv, &AllBound_CMx_Inv, 1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CMy_Inv, &AllBound_CMy_Inv, 1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CMz_Inv, &AllBound_CMz_Inv, 1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CFx_Inv, &AllBound_CFx_Inv, 1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CFy_Inv, &AllBound_CFy_Inv, 1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CFz_Inv, &AllBound_CFz_Inv, 1, MPI::DOUBLE, MPI::SUM);  
#endif

  
  
	Total_CDrag += AllBound_CDrag_Inv;
	Total_CLift += AllBound_CLift_Inv;
	Total_CSideForce += AllBound_CSideForce_Inv;
	Total_CMx += AllBound_CMx_Inv;
	Total_CMy += AllBound_CMy_Inv;
	Total_CMz += AllBound_CMz_Inv;
	Total_CEff += AllBound_CEff_Inv;
	Total_CFx += AllBound_CFx_Inv;
	Total_CFy += AllBound_CFy_Inv;
	Total_CFz += AllBound_CFz_Inv;
  
}

void CTNE2EulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solution_container, CConfig *config) {
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

void CTNE2EulerSolver::SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint, iEdge, iVertex;
	unsigned short iDim, iVar, iMarker;
	double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
	Partial_Gradient, Partial_Res, *Normal;
  
	/*--- Gradient primitive variables compressible (temp, vx, vy, vz, P, rho)
   Gradient primitive variables incompressible (rho, vx, vy, vz, beta) ---*/
	PrimVar_Vertex = new double [nPrimVarGrad];
	PrimVar_i = new double [nPrimVarGrad];
	PrimVar_j = new double [nPrimVarGrad];
  
	/*--- Set Gradient_Primitive to zero ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		node[iPoint]->SetGradient_PrimitiveZero(nPrimVarGrad);
  
	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
    
		for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
			PrimVar_i[iVar] = node[iPoint]->GetPrimVar(iVar);
			PrimVar_j[iVar] = node[jPoint]->GetPrimVar(iVar);
		}
    
		Normal = geometry->edge[iEdge]->GetNormal();
		for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
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
        
				for (iVar = 0; iVar < nPrimVarGrad; iVar++)
					PrimVar_Vertex[iVar] = node[iPoint]->GetPrimVar(iVar);
        
				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				for (iVar = 0; iVar < nPrimVarGrad; iVar++)
					for (iDim = 0; iDim < nDim; iDim++) {
						Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
						node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
					}
			}
		}
	}
  
	/*--- Update gradient value ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
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

void CTNE2EulerSolver::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, jDim, iNeigh;
	unsigned long iPoint, jPoint;
	double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
	r23_b, r33, weight, product;
  
	/*--- Gradient primitive variables compressible (temp, vx, vy, vz, P, rho)
   Gradient primitive variables incompressible (rho, vx, vy, vz, beta) ---*/
	PrimVar_i = new double [nPrimVarGrad];
	PrimVar_j = new double [nPrimVarGrad];
  
	/*--- Loop over points of the grid ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		Coord_i = geometry->node[iPoint]->GetCoord();
    
		for (iVar = 0; iVar < nPrimVarGrad; iVar++)
			PrimVar_i[iVar] = node[iPoint]->GetPrimVar(iVar);
    
		/*--- Inizialization of variables ---*/
		for (iVar = 0; iVar < nPrimVarGrad; iVar++)
			for (iDim = 0; iDim < nDim; iDim++)
				cvector[iVar][iDim] = 0.0;
    
		r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0;
		r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;
    
		for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
			Coord_j = geometry->node[jPoint]->GetCoord();
      
			for (iVar = 0; iVar < nPrimVarGrad; iVar++)
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
			for (iVar = 0; iVar < nPrimVarGrad; iVar++)
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
		for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
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


void CTNE2EulerSolver::SetPrimVar_Limiter(CGeometry *geometry, CConfig *config) { }

void CTNE2EulerSolver::SetPreconditioner(CConfig *config, unsigned short iPoint) {
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

void CTNE2EulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solution_container,
                                     CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iSpecies, jSpecies, iVar, jVar;
	unsigned long iPoint, iVertex;
  bool implicit;
  double *Normal, *UnitNormal, *Ms, *dPdrhos;
  double Area, rhoCvtr, rhoCvve, rho_el, Ru;
  double rho, cs, u, v, w, P, rhoE, rhoEve, conc, Beta;
  
  /*--- Allocate arrays ---*/
  UnitNormal = new double[3];
  
  /*--- Set booleans based on configuration options ---*/
	implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  
  /*--- Load parameters from the config class ---*/
  Ms = config->GetMolar_Mass();
  
  /*--- Rename for convenience ---*/
  Ru = UNIVERSAL_GAS_CONSTANT;
  
	/*--- Loop over all the vertices on this boundary (val_marker) ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			/*--- Calculate parameters from the geometry ---*/
      // Note: The vertex normal points out of the geometry by convention,
      //       so to calculate the influence from the boundary condition
      //       to the domain, we negate this vector
      Area   = 0.0;
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
			for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;
      
			/*--- Retrieve the pressure on the vertex ---*/
      P   = node[iPoint]->GetPressure();
      
      /*--- Apply the flow-tangency b.c. to the convective flux ---*/
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        Residual[iSpecies] = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Residual[nSpecies+iDim] = P * UnitNormal[iDim] * Area;
      Residual[nSpecies+nDim]   = 0.0;
			Residual[nSpecies+nDim+1] = 0.0;
      
			/*--- Add value to the residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);
      
			/*--- If using implicit time-stepping, calculate b.c. contribution to Jacobian ---*/
			if (implicit) {
        
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        
        rho     = node[iPoint]->GetDensity();
        u       = node[iPoint]->GetVelocity(0);
        v       = node[iPoint]->GetVelocity(1);
        w       = node[iPoint]->GetVelocity(2);
        rhoCvtr = node[iPoint]->GetRhoCv_tr();
        rhoCvve = node[iPoint]->GetRhoCv_ve();
        rhoE    = node[iPoint]->GetSolution(nSpecies+nDim);
        rhoEve  = node[iPoint]->GetSolution(nSpecies+nDim+1);
        dPdrhos = node[iPoint]->GetdPdrhos();
        
        /*--- If free electrons are present, retrieve the electron gas density ---*/
        if (config->GetIonization()) rho_el = node[iPoint]->GetMassFraction(nSpecies-1) * rho;
        else                         rho_el = 0.0;
        
        conc = 0.0;
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          cs    = node[iPoint]->GetMassFraction(iSpecies);
          conc += cs * rho/Ms[iSpecies];
          
          Jacobian_i[nSpecies][iSpecies]   = dPdrhos[iSpecies] * UnitNormal[0];
          Jacobian_i[nSpecies+1][iSpecies] = dPdrhos[iSpecies] * UnitNormal[1];
          Jacobian_i[nSpecies+2][iSpecies] = dPdrhos[iSpecies] * UnitNormal[2];
          Jacobian_i[nSpecies+3][iSpecies] = 0.0;
          Jacobian_i[nSpecies+4][iSpecies] = 0.0;
          
          Jacobian_i[iSpecies][nSpecies]   = cs * UnitNormal[0];
          Jacobian_i[iSpecies][nSpecies+1] = cs * UnitNormal[1];
          Jacobian_i[iSpecies][nSpecies+2] = cs * UnitNormal[2];
          Jacobian_i[iSpecies][nSpecies+3] = 0.0;
          Jacobian_i[iSpecies][nSpecies+4] = 0.0;
        }
        
        Beta = Ru*conc/rhoCvtr;
        
        Jacobian_i[nSpecies][nSpecies]     = u*UnitNormal[0] - u*Beta*UnitNormal[0];
        Jacobian_i[nSpecies][nSpecies+1]   = u*UnitNormal[1] - v*Beta*UnitNormal[0];
        Jacobian_i[nSpecies][nSpecies+2]   = u*UnitNormal[2] - w*Beta*UnitNormal[0];
        Jacobian_i[nSpecies][nSpecies+3]   = Beta*UnitNormal[0];
        Jacobian_i[nSpecies][nSpecies+4]   = (-Beta + (rho_el/Ms[nSpecies-1])*Ru/rhoCvve)*UnitNormal[0];
        
        Jacobian_i[nSpecies+1][nSpecies]   = v*UnitNormal[0] - u*Beta*UnitNormal[1];
        Jacobian_i[nSpecies+1][nSpecies+1] = v*UnitNormal[1] - v*Beta*UnitNormal[1];
        Jacobian_i[nSpecies+1][nSpecies+2] = v*UnitNormal[2] - w*Beta*UnitNormal[1];
        Jacobian_i[nSpecies+1][nSpecies+3] = Beta*UnitNormal[1];
        Jacobian_i[nSpecies+1][nSpecies+4] = (-Beta + (rho_el/Ms[nSpecies-1])*Ru/rhoCvve)*UnitNormal[1];
        
        Jacobian_i[nSpecies+2][nSpecies]   = w*UnitNormal[0] - u*Beta*UnitNormal[2];
        Jacobian_i[nSpecies+2][nSpecies+1] = w*UnitNormal[1] - v*Beta*UnitNormal[2];
        Jacobian_i[nSpecies+2][nSpecies+2] = w*UnitNormal[2] - w*Beta*UnitNormal[2];
        Jacobian_i[nSpecies+2][nSpecies+3] = Beta*UnitNormal[2];
        Jacobian_i[nSpecies+2][nSpecies+4] = (-Beta + (rho_el/Ms[nSpecies-1])*Ru/rhoCvve)*UnitNormal[2];
        
        Jacobian_i[nSpecies+3][nSpecies]   = (rhoE+P)/rho * UnitNormal[0];
        Jacobian_i[nSpecies+3][nSpecies+1] = (rhoE+P)/rho * UnitNormal[1];
        Jacobian_i[nSpecies+3][nSpecies+2] = (rhoE+P)/rho * UnitNormal[2];
        Jacobian_i[nSpecies+3][nSpecies+3] = 0.0;
        Jacobian_i[nSpecies+3][nSpecies+4] = 0.0;
        
        Jacobian_i[nSpecies+4][nSpecies]   = rhoEve/rho * UnitNormal[0];
        Jacobian_i[nSpecies+4][nSpecies+1] = rhoEve/rho * UnitNormal[1];
        Jacobian_i[nSpecies+4][nSpecies+2] = rhoEve/rho * UnitNormal[2];
        Jacobian_i[nSpecies+4][nSpecies+3] = 0.0;
        Jacobian_i[nSpecies+4][nSpecies+4] = 0.0;
        
        /*--- Integrate over the dual-grid area ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = Jacobian_i[iVar][jVar] * Area;
        
        /*--- Apply the contribution to the system ---*/
        Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
        
			}
		}
	}
}

void CTNE2EulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solution_container,
                                    CNumerics *conv_numerics, CNumerics *visc_numerics,
                                    CConfig *config, unsigned short val_marker) {
	unsigned short iDim;
	unsigned long iVertex, iPoint, Point_Normal;
  bool implicit, viscous;
  double *U_domain, *V_domain, *U_infty, *V_infty, *Normal;
  
  /*--- Set booleans from configuration parameters ---*/
  implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	viscous  = config->GetViscous();
  
  /*--- Allocate arrays ---*/
	Normal = new double[nDim];
  
  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  conv_numerics->SetRhosIndex   ( node[0]->GetRhosIndex()    );
  conv_numerics->SetRhoIndex    ( node[0]->GetRhoIndex()     );
  conv_numerics->SetPIndex      ( node[0]->GetPIndex()       );
  conv_numerics->SetTIndex      ( node[0]->GetTIndex()       );
  conv_numerics->SetTveIndex    ( node[0]->GetTveIndex()     );
  conv_numerics->SetVelIndex    ( node[0]->GetVelIndex()     );
  conv_numerics->SetHIndex      ( node[0]->GetHIndex()       );
  conv_numerics->SetAIndex      ( node[0]->GetAIndex()       );
  conv_numerics->SetRhoCvtrIndex( node[0]->GetRhoCvtrIndex() );
  conv_numerics->SetRhoCvveIndex( node[0]->GetRhoCvveIndex() );
  
	/*--- Loop over all the vertices on this boundary (val_marker) ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			/*--- Retrieve index of the closest interior node ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
			/*--- Pass boundary node normal to CNumerics ---*/
      // Note: The vertex normal points out of the geometry by convention,
      //       so to calculate the influence from the boundary condition
      //       to the domain, we negate this vector
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
			conv_numerics->SetNormal(Normal);
      
			/*--- Retrieve solution at the boundary node & free-stream ---*/
      U_domain = node[iPoint]->GetSolution();
      V_domain = node[iPoint]->GetPrimVar();
      U_infty  = node_infty->GetSolution();
      V_infty  = node_infty->GetPrimVar();
      
      /*--- Pass conserved & primitive variables to CNumerics ---*/
      conv_numerics->SetConservative(U_domain, U_infty);
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Pass supplementary information to CNumerics ---*/
      conv_numerics->SetdPdrhos(node[iPoint]->GetdPdrhos(), node_infty->GetdPdrhos());
      
			/*--- Compute the convective residual (and Jacobian) ---*/
      // Note: This uses the specified boundary num. method specified in definition_structure.cpp
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Apply contribution to the linear system ---*/
      LinSysRes.AddBlock(iPoint, Residual);
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
/*      unsigned short iVar, jVar;
      cout << "TNE2 BC Far Field: " << endl;
      cout << "dPdrhos[0]: " << node[iPoint]->GetdPdrhos()[0] << endl;
      cout << "dPdrhos_inf[0]: " << node_infty->GetdPdrhos()[0] << endl;
      cout << "ConsVarDomain: " << endl;
      for (iVar = 0; iVar < nVar; iVar++)
        cout << U_domain[iVar] << endl;
      cout << endl << endl << "PrimVarDomain: " << endl;
      for (iVar = 0; iVar < nPrimVar; iVar++)
        cout << V_domain[iVar] << endl;
      cout << endl << endl << "ConsVarInfty: " << endl;
      for (iVar = 0; iVar < nVar; iVar++)
        cout << U_infty[iVar] << endl;
      cout << endl << endl << "PrimVarInfty: " << endl;
      for (iVar = 0; iVar < nPrimVar; iVar++)
        cout << V_infty[iVar] << endl;
      
      cout << endl << endl << "Residual: " << endl;
      for (iVar = 0; iVar < nVar; iVar++)
        cout << Residual[iVar] << endl;
      cout << endl << endl << "Jacobian: " << endl;
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          cout << Jacobian_i[iVar][jVar] << "\t";
        }
        cout << endl;
      }
      cin.get();*/
 
      
      
			/*--- Viscous contribution ---*/
			if (viscous) {
        cout << "WARNING!!!  BC_Far_Field: Viscous contribution to boundary not implemented!" << endl;
			}
		}
	}
  
	/*--- Free locally allocated memory ---*/
	delete [] Normal;
}

void CTNE2EulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solution_container,
                                  CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iVar, iDim;
	unsigned long iVertex, iPoint, Point_Normal;
	double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
	Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
	alpha, aa, bb, cc, dd, Area, UnitaryNormal[3];
  
	bool implicit             = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool grid_movement        = config->GetGrid_Movement();
	double Two_Gamma_M1       = 2.0/Gamma_Minus_One;
	double Gas_Constant       = config->GetGas_ConstantND();
	unsigned short Kind_Inlet = config->GetKind_Inlet();
	string Marker_Tag         = config->GetMarker_All_Tag(val_marker);
	bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  
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
            Velocity[iDim] = node[iPoint]->GetVelocity(iDim);
          Pressure    = node[iPoint]->GetPressure(false);
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
      
			/*--- Set various quantities in the solver class ---*/
			conv_numerics->SetConservative(U_domain, U_inlet);
      
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
        
				/*--- Laminar viscosity ---*/
				visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());
        
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

void CTNE2EulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solution_container,
                                   CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iVar, iDim;
	unsigned long iVertex, iPoint, Point_Normal;
	double Pressure, P_Exit, Velocity[3],
	Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
	Area, UnitaryNormal[3];
  
	bool implicit           = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  double Gas_Constant     = config->GetGas_ConstantND();
	bool grid_movement      = config->GetGrid_Movement();
	string Marker_Tag       = config->GetMarker_All_Tag(val_marker);
	bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  
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
      
			/*--- Set various quantities in the solver class ---*/
			conv_numerics->SetConservative(U_domain, U_outlet);
      
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
        
				/*--- Laminar viscosity ---*/
				visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());
        
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

void CTNE2EulerSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solution_container,
                                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iDim, iVar;
	unsigned long iVertex, iPoint, Point_Normal;
	double Density, Pressure, Temperature, Energy, *Velocity, Velocity2;
	double Gas_Constant = config->GetGas_ConstantND();
  
	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool grid_movement  = config->GetGrid_Movement();
	bool viscous              = config->GetViscous();
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
  
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
        
				/*--- Laminar viscosity ---*/
				visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());
        
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

void CTNE2EulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Call the Euler wall routine ---*/
  BC_Euler_Wall(geometry, solver_container, conv_numerics, config, val_marker);
}

void CTNE2EulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solution_container, CConfig *config,
                                              unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
	unsigned short iVar, jVar;
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1;
	double Volume_nM1, Volume_n, Volume_nP1, TimeStep;
  
	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
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
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}
  
}

void CTNE2EulerSolver::GetRestart(CGeometry *geometry, CConfig *config, unsigned short val_iZone) {
  
  unsigned short iVar;
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
      
      /*--- First value is the point index, then the conservative variables ---*/
      point_line >> index;
      
      for (iVar = 0; iVar < nVar; iVar++)
        point_line >> Solution[iVar];
            
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


void CTNE2EulerSolver::SetVolume_Output(CConfig *config, CGeometry *geometry, double **data_container, unsigned short nOutput_Vars) {
  
#ifdef DEBUG_TDE
  
	unsigned short iVar;
	unsigned long iPoint;
  
	/*--- Add up total number of output variables to be written. ---*/
	nOutput_Vars = nVar;
  
	for (iVar = 0; iVar < config->GetnOutput_Vars_Vol(); iVar++ ) {
    
		switch(config->GetOutput_Vars_Vol(iVar)) {
      case PRESSURE:
        nOutput_Vars++;
        break;
      case MACH:
        nOutput_Vars++;
        break;
		}
    
	}
  
	// NEEDS TO BE MAX NUMBER OF POINTS ON ANY PARTITION ?
	data_container = new double*[nOutput_Vars];
	for (iVar = 0; iVar < nOutput_Vars; iVar++ ) {
		data_container[iVar] = new double[nPointDomain];
	}
  
	for (iVar = 0; iVar < config->GetnOutput_Vars_Vol(); iVar++ ) {
    
		switch(config->GetOutput_Vars_Vol(iVar)) {
      case PRESSURE:
        nOutput_Vars++;
        break;
      case MACH:
        nOutput_Vars++;
        break;
		}
    
	}
#endif
}

CTNE2NSSolver::CTNE2NSSolver(void) : CTNE2EulerSolver() {
  
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
  
	ForceViscous = NULL;
	MomentViscous = NULL;
	CSkinFriction = NULL;
  
}

CTNE2NSSolver::CTNE2NSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTNE2EulerSolver() {
	unsigned long iPoint, index, counter_local = 0, counter_global = 0;
	unsigned short iVar, iDim, iMarker;
  double Density, Velocity2, Pressure, Temperature;
  
	unsigned short nZone = geometry->GetnZone();
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
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
	nVar = nDim + 2; nPrimVar = nDim+5; nPrimVarGrad = nDim+3;
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
  
	if ((config->GetKind_Upwind_TNE2() == ROE_TURKEL_2ND) || (config->GetKind_Upwind_TNE2() == ROE_TURKEL_1ST)) {
		Precon_Mat_inv = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar ++)
			Precon_Mat_inv[iVar] = new double[nVar];
	}
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
	/*--- Jacobians and vector structures for implicit computations ---*/
	if (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT) {
    
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
  
	/*--- Init total coefficients ---*/
	Total_CDrag = 0.0;	Total_CLift = 0.0;			Total_CSideForce = 0.0;
	Total_CMx = 0.0;		Total_CMy = 0.0;				Total_CMz = 0.0;
	Total_CEff = 0.0;
	Total_CFx = 0.0;		Total_CFy = 0.0;				Total_CFz = 0.0;
  Total_Q = 0.0;
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
  Q_Visc = new double[nMarker];
  Maxq_Visc = new double[nMarker];
  
	/*--- Read farfield conditions from config ---*/
//	Density_Inf   = config->GetDensity_FreeStreamND();
	Pressure_Inf  = config->GetPressure_FreeStreamND();
	Velocity_Inf  = config->GetVelocity_FreeStreamND();
//	Energy_Inf    = config->GetEnergy_FreeStreamND();
	Viscosity_Inf = config->GetViscosity_FreeStreamND();
	Mach_Inf      = config->GetMach_FreeStreamND();
	Prandtl_Lam   = config->GetPrandtl_Lam();
	Prandtl_Turb  = config->GetPrandtl_Turb();
  
	/*--- Inlet/Outlet boundary conditions, using infinity values ---*/
//	Density_Inlet = Density_Inf;		Density_Outlet = Density_Inf;
	Pressure_Inlet = Pressure_Inf;	Pressure_Outlet = Pressure_Inf;
//	Energy_Inlet = Energy_Inf;			Energy_Outlet = Energy_Inf;
	Mach_Inlet = Mach_Inf;					Mach_Outlet = Mach_Inf;
	Velocity_Inlet  = new double [nDim]; Velocity_Outlet = new double [nDim];
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_Inlet[iDim] = Velocity_Inf[iDim];
		Velocity_Outlet[iDim] = Velocity_Inf[iDim];
	}
  
	/*--- Restart the solution from file information ---*/
	if (!restart || geometry->GetFinestMGLevel() == false || nZone > 1) {
    
		/*--- Restart the solution from infinity ---*/
//		for (iPoint = 0; iPoint < nPoint; iPoint++)
//			node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
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
        if (nDim == 2) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
        if (nDim == 3) point_line >> index >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
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
  counter_local = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    Density = node[iPoint]->GetSolution(0);
    Velocity2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);
    Pressure    = Gamma_Minus_One*Density*(node[iPoint]->GetSolution(nDim+1)/Density-0.5*Velocity2);
    Temperature = Pressure / ( Gas_Constant * Density);
    if ((Pressure < 0.0) || (Temperature < 0.0)) {
//      Solution[0] = Density_Inf;
      for (iDim = 0; iDim < nDim; iDim++)
//        Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
//      Solution[nDim+1] = Energy_Inf*Density_Inf;
      
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
  
  
  if ((rank == MASTER_NODE) && (counter_global != 0)) cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  
	/*--- Define solver parameters needed for execution of destructor ---*/
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
	else least_squares = false;
  
	if ((config->GetKind_Upwind_TNE2() == ROE_TURKEL_2ND) ||
      (config->GetKind_Upwind_TNE2() == ROE_TURKEL_1ST)) roe_turkel = true;
	else roe_turkel = false;
  
	/*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);
  
}

CTNE2NSSolver::~CTNE2NSSolver(void) {
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

void CTNE2NSSolver::Preprocessing(CGeometry *geometry, CSolver **solution_container, CConfig *config,
                                  unsigned short iMesh, unsigned short iRKStep,
                                  unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
  
	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool upwind_2nd = ((config->GetKind_Upwind_TNE2() == ROE_2ND) || (config->GetKind_Upwind_TNE2() == AUSM_2ND)
                     || (config->GetKind_Upwind_TNE2() == HLLC_2ND) || (config->GetKind_Upwind_TNE2() == ROE_TURKEL_2ND));
	bool limiter = (config->GetKind_SlopeLimit_TNE2() != NONE);
	double Gas_Constant = config->GetGas_ConstantND();
  
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
		/*--- Set the primitive variables incompressible (dens, vx, vy, vz, beta)
     and compressible (temp, vx, vy, vz, press, dens, enthal, sos)---*/
		node[iPoint]->SetPrimVar_Compressible(config);
    
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
  
	/*--- Compute gradient of the primitive variables ---*/
	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetPrimVar_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetPrimVar_Gradient_LS(geometry, config);
  
	/*--- Initialize the jacobian matrices ---*/
	if (implicit) Jacobian.SetValZero();
  
}

void CTNE2NSSolver::SetTime_Step(CGeometry *geometry, CSolver **solution_container, CConfig *config,
                                 unsigned short iMesh, unsigned long Iteration) {
	double *Normal, Area, Vol, Mean_SoundSpeed, Mean_ProjVel, Lambda, Local_Delta_Time, Local_Delta_Time_Visc,
	Global_Delta_Time = 1E6, Mean_LaminarVisc, Mean_Density, Lambda_1, Lambda_2, K_v = 0.25, Global_Delta_UnstTimeND;
	unsigned long iEdge, iVertex, iPoint = 0, jPoint = 0;
	unsigned short iDim, iMarker;
	double ProjVel, ProjVel_i, ProjVel_j;
  
	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool grid_movement = config->GetGrid_Movement();
	bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
	Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
  
	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed and viscosity ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		node[iPoint]->SetMax_Lambda_Inv(0.0);
		node[iPoint]->SetMax_Lambda_Visc(0.0);
	}
  
	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Point identification, Normal vector and area ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
    
		Normal = geometry->edge[iEdge]->GetNormal();
		Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
		/*--- Mean Values ---*/
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    
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
    Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
    Mean_Density     = 0.5*(node[iPoint]->GetSolution(0) + node[jPoint]->GetSolution(0));
    
		Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
		Lambda_2 = 0.0;
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
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
      
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
      Mean_LaminarVisc = node[iPoint]->GetLaminarViscosity();
      Mean_Density     = node[iPoint]->GetSolution(0);
      
			Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
			Lambda_2 = 0.0;
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

void CTNE2NSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solution_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iPoint, jPoint, iEdge;
  
	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points, coordinates and normal vector in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive variables, and gradient ---*/
    numerics->SetPrimitive(node[iPoint]->GetPrimVar(), node[jPoint]->GetPrimVar());
    numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[jPoint]->GetGradient_Primitive());
    
    /*--- Laminar viscosity ---*/
    numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[jPoint]->GetLaminarViscosity());
    
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

void CTNE2NSSolver::Viscous_Forces(CGeometry *geometry, CConfig *config) {
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
  
	bool grid_movement  = config->GetGrid_Movement();
  
	/*--- If we have a rotating frame problem or an unsteady problem with
   mesh motion, use special reference values for the force coefficients.
   Otherwise, use the freestream values, which is the standard convention. ---*/
  
	if (grid_movement) {
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
	AllBound_CEff_Visc = 0.0;
	AllBound_Q_Visc = 0.0;  AllBound_Maxq_Visc = 0.0;
  
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
        Viscosity = node[iPoint]->GetLaminarViscosity();
        Density = node[iPoint]->GetDensity();
        
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
				for (iDim = 0; iDim < nDim; iDim++) Vel[iDim] = node[iPointNormal]->GetVelocity(iDim);
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
          Q_Visc[iMarker]  = HeatLoad;
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
          Q_Visc[iMarker] = HeatLoad;
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
        AllBound_Q_Visc += Q_Visc[iMarker];
        if (Maxq_Visc[iMarker] > AllBound_Maxq_Visc)
          AllBound_Maxq_Visc = Maxq_Visc[iMarker];
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
  Total_Q += AllBound_Q_Visc;
  Total_Maxq = AllBound_Maxq_Visc;
  
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Tau[iDim];
	delete [] Tau;
	delete [] UnitaryNormal;
	delete [] TauTangent;
	delete [] TauElem;
}

void CTNE2NSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solution_container,
                                     CNumerics *conv_numerics, CNumerics *visc_numerics,
                                     CConfig *config, unsigned short val_marker) {
  
	/*--- Local variables ---*/
	unsigned short iDim, iVar;
	unsigned long iVertex, iPoint, total_index;
  
	double Wall_HeatFlux;
	double *Grid_Vel, *Normal, Area;
  
	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool grid_movement  = config->GetGrid_Movement();
  
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
			if (grid_movement) {
				Grid_Vel = geometry->node[iPoint]->GetGridVel();
				for (iDim = 0; iDim < nDim; iDim++)
					Vector[iDim] = Grid_Vel[iDim];
			} else {
				for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
			}
      
			/*--- Set the residual, truncation error, and velocity value ---*/
			node[iPoint]->SetVelocity_Old(Vector, false);
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
			node[iPoint]->SetVel_ResTruncError_Zero();
      
			/*--- Set the residual on the boundary with the specified heat flux ---*/
			Res_Visc[nDim+1] = Wall_HeatFlux * Area;
      
			/*--- Flux contribution due to grid motion (energy equation) ---*/
			if (grid_movement) {
				double ProjGridVel = 0.0;
				Grid_Vel = geometry->node[iPoint]->GetGridVel();
				for (iDim = 0; iDim < nDim; iDim++)
					ProjGridVel += Grid_Vel[iDim]*(-1.0)*Normal[iDim];
				Res_Conv[nDim+1] = node[iPoint]->GetPressure(false)*ProjGridVel;
			}
      
			/*--- Viscous contribution for moving walls ---*/
			if (grid_movement) {
        
				unsigned short jDim;
				double total_viscosity, div_vel, Density;
				Density = node[iPoint]->GetSolution(0);
				double val_laminar_viscosity = node[iPoint]->GetLaminarViscosity();
				double **val_gradprimvar = node[iPoint]->GetGradient_Primitive();
				double tau[3][3], delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
				total_viscosity = val_laminar_viscosity;
        
				/*--- Get the appropriate grid velocity at this node ---*/
				if (grid_movement) {
					Grid_Vel = geometry->node[iPoint]->GetGridVel();
				}
        
				double Flux_Tensor[nVar][nDim];
				for (iVar = 0 ; iVar < nVar; iVar++)
					for (jDim = 0 ; jDim < nDim; jDim++)
						Flux_Tensor[iVar][jDim] = 0.0;
        
				div_vel = 0.0;
				for (iDim = 0 ; iDim < nDim; iDim++)
					div_vel += val_gradprimvar[iDim+1][iDim];
        
				for (iDim = 0 ; iDim < nDim; iDim++)
					for (jDim = 0 ; jDim < nDim; jDim++)
						tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
						- TWO3*total_viscosity*div_vel*delta[iDim][jDim];
        
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

void CTNE2NSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solution_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
	unsigned long iVertex, iPoint, Point_Normal, total_index;
	unsigned short iVar, jVar, iDim;
	double *Normal, *Coord_i, *Coord_j, Area, dist_ij, *Grid_Vel;
	double UnitaryNormal[3];
	double Twall, Temperature, dTdn, dTdrho;
	double Density, Vel2, Energy;
	double Laminar_Viscosity, Thermal_Conductivity, Gas_Constant, cp;
	double Theta;
	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool grid_movement = config->GetGrid_Movement();
  
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
			if (grid_movement) {
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
			node[iPoint]->SetVelocity_Old(Vector, false);
			for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
			node[iPoint]->SetVel_ResTruncError_Zero();
      
			/*--- Retrieve temperatures from boundary nearest neighbor --*/
			Temperature = node[Point_Normal]->GetPrimVar(0);
      
			/*--- Calculate temperature gradient normal to the wall using FD ---*/
			dTdn                  = (Twall - Temperature)/dist_ij;
			Laminar_Viscosity     = node[iPoint]->GetLaminarViscosity();
			Thermal_Conductivity  = cp * (Laminar_Viscosity/PRANDTL);
      
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
