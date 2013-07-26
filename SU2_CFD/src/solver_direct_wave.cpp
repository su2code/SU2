/*!
 * \file solution_direct_wave.cpp
 * \brief Main subrotuines for solving the wave equation.
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

CWaveSolver::CWaveSolver(void) : CSolver() { }

CWaveSolver::CWaveSolver(CGeometry *geometry, 
                             CConfig *config) : CSolver() {
	unsigned short nMarker, iVar;
  
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
	nDim    = geometry->GetnDim();
	nMarker = config->GetnMarker_All(); 
	node    = new CVariable*[nPoint];
	nVar    = 2; // solve as a 2 eq. system		
  
	Residual     = new double[nVar]; Residual_RMS = new double[nVar];
	Solution     = new double[nVar];
  Res_Sour     = new double[nVar];
  Residual_Max = new double[nVar]; Point_Max = new unsigned long[nVar];
  

	/*--- Point to point stiffness matrix (only for triangles)---*/	
	StiffMatrix_Elem = new double*[nDim+1];
	for (iVar = 0; iVar < nDim+1; iVar++) {
		StiffMatrix_Elem[iVar] = new double [nDim+1];
	}
  
	StiffMatrix_Node = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		StiffMatrix_Node[iVar] = new double [nVar];
	}
  
	/*--- Initialization of matrix structures ---*/
	StiffMatrixSpace.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
	StiffMatrixTime.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
	Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
  
  /*--- Initialization of linear solver structures ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /* Wave strength coefficient for all of the markers */
  
  CWave = new double[config->GetnMarker_All()];
  Total_CWave = 0.0;
  
  /* Check for a restart (not really used), initialize from zero otherwise */
  
	bool restart = (config->GetRestart());
	if (!restart) {
    
    double *Coord, Radius;
		for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
      
      /*--- Set up the initial condition for the drum problem ---*/
      Coord = geometry->node[iPoint]->GetCoord();
      Radius = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        Radius += Coord[iDim]*Coord[iDim];
      Radius = sqrt(Radius);
      
      /*--- Symmetrically plucked drum ---*/
      //      Solution[0] = 1.0 - Radius; 
      //      Solution[1] = 0.0;
      
      /*--- Off-center strike ---*/
//      if ((Radius > 0.4) && (Radius < 0.6)) {
//        Solution[0] = 10.0; 
//      } else
//          Solution[0] = 0.0;
//      Solution[1] = 0.0;
      
      /*--- Struck drum ---*/
//      Solution[0] = 0.0; 
//      Solution[1] = -1.0;
      
      /*--- Zero initial condition for testing source terms & forcing BCs ---*/
      Solution[0] = 0.0; 
      Solution[1] = 0.0;
      
			node[iPoint] = new CWaveVariable(Solution, nDim, nVar, config);
      
      /* Copy solution to old containers if using dual time */
      
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST) {
        node[iPoint]->Set_Solution_time_n(); 
      } else if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) {
        node[iPoint]->Set_Solution_time_n();
        node[iPoint]->Set_Solution_time_n1();
      }
      
    }
	} else {
    
    cout << "Wave restart file not currently configured!!" << endl;
    cout << "Press any key to exit..." << endl;
    cin.get();
    exit(1);
    
		string mesh_filename = config->GetSolution_FlowFileName();
		ifstream restart_file;
    
		char *cstr; cstr = new char [mesh_filename.size()+1];
		strcpy (cstr, mesh_filename.c_str());
		restart_file.open(cstr, ios::in);
        
		if (restart_file.fail()) {
			cout << "There is no wave restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}
		unsigned long index;
		string text_line;
    
		for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			point_line >> index >> Solution[0] >> Solution[1];
			node[iPoint] = new CWaveVariable(Solution, nDim, nVar, config);
		}
		restart_file.close();
	}
  
  /* Load all of the matrices and Jacobians that are fixed */
  /* during the simulation in order to save effort.        */
  //SetSpace_Matrix(geometry, Solution, config);
  //SetTime_Matrix(geometry, config);
  
}

CWaveSolver::~CWaveSolver(void) {
  
	unsigned short iVar;
  
	delete [] Residual;
	delete [] Residual_Max;
	delete [] Solution;
  
	for (iVar = 0; iVar < nDim+1; iVar++)
		delete [] StiffMatrix_Elem[iVar];
  
  
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] StiffMatrix_Node[iVar];
  
	delete [] StiffMatrix_Elem;
	delete [] StiffMatrix_Node;
  
}

void CWaveSolver::Preprocessing(CGeometry *geometry, 
                                  CSolver **solver_container,
                                  CConfig   *config,
                                  unsigned short iMesh,
                                  unsigned short iRKStep,
                                  unsigned short RunTime_EqSystem) {
  
  /* Set residuals and matrix entries to zero */
  
	for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    LinSysRes.SetBlock_Zero(iPoint);
	}
	
  /* Zero out the entries in the various matrices */
  
	StiffMatrixSpace.SetValZero();
	StiffMatrixTime.SetValZero();
	Jacobian.SetValZero();
  
}

void CWaveSolver::Source_Residual(CGeometry *geometry, 
                                    CSolver **solver_container, 
                                    CNumerics *numerics, CNumerics *second_numerics,
                                    CConfig   *config, 
                                    unsigned short iMesh) {

  /* Local variables and initialization */
  
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double a[3], b[3], Area_Local, Time_Num, Time_Phys;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
	unsigned short iDim;
	
  /*--- Numerical time step. This system is unconditionally stable,
    so a very big step can be used. ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING) 
    Time_Num = config->GetDelta_UnstTimeND();
	else Time_Num = 1E+30;
  
  /*--- Physical timestep for source terms ---*/
  Time_Phys = config->GetDelta_UnstTimeND();
  
	/* Loop through elements to compute contributions from the matrix     */
  /* blocks involving time. These contributions are also added to the   */
  /* Jacobian w/ the time step. Spatial source terms are also computed. */
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
    /* Get node numbers and their coordinate vectors */
    
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		
    Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		
    /* Compute triangle area (2-D) */
		
    for (iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
		Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
		
    /*--- Block contributions to the Jacobian (includes time step) ---*/
    
		StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] =  1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] =  1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
       
    /*--- Noise source term computation for aeroacoustic problems ---*/
    
    if ((config->GetKind_Solver() == AEROACOUSTIC_EULER) || 
        (config->GetKind_Solver() == AEROACOUSTIC_NAVIER_STOKES) ||
        (config->GetKind_Solver() == AEROACOUSTIC_RANS)) { 
      
      double Noise_0, Noise_1, Noise_2, Noise_Total = 0.0;
      
      /* Thickness noise contribution */
      Noise_0 = node[Point_0]->GetThickness_Noise();
      Noise_1 = node[Point_1]->GetThickness_Noise();
      Noise_2 = node[Point_2]->GetThickness_Noise();
      Noise_Total += (1.0/3.0)*(Noise_0+Noise_1+Noise_2);
      
      Noise_0 = node[Point_0]->GetLoading_Noise();
      Noise_1 = node[Point_1]->GetLoading_Noise();
      Noise_2 = node[Point_2]->GetLoading_Noise();
      Noise_Total += (1.0/3.0)*(Noise_0+Noise_1+Noise_2);
      
      Noise_0 = node[Point_0]->GetQuadrupole_Noise();
      Noise_1 = node[Point_1]->GetQuadrupole_Noise();
      Noise_2 = node[Point_2]->GetQuadrupole_Noise();
      Noise_Total += (1.0/3.0)*(Noise_0+Noise_1+Noise_2);
      
      /* Source terms are only applied to the second equation */
      Res_Sour[0] = 0.0;
      Res_Sour[1] = Noise_Total;
      
      /* CHECK SIGN!!! */
      LinSysRes.SubtractBlock(Point_0, Res_Sour);
      LinSysRes.SubtractBlock(Point_1, Res_Sour);
      LinSysRes.SubtractBlock(Point_2, Res_Sour);
      
    }
    
    /*--- Add a sponge condition in the center of the mesh ---*/
    double Sponge_Radius = 1.2; double eps = 0.1;
    double Sponge_0 = 0.0, Sponge_1 = 0.0;
    double Radius_0 =0.0, Radius_1 = 0.0, Radius_2 = 0.0, *Solution;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Radius_0 += Coord_0[iDim]*Coord_0[iDim];
      Radius_1 += Coord_1[iDim]*Coord_1[iDim];
      Radius_2 += Coord_2[iDim]*Coord_2[iDim];
    }
    Radius_0 = sqrt(Radius_0); Radius_1 = sqrt(Radius_1); Radius_2 = sqrt(Radius_2);
    
    /* Check if we're inside sponge radius and compute contributions */
    if (Radius_0 < Sponge_Radius) {
      Solution = node[Point_0]->GetSolution();
      Sponge_0 += (Sponge_Radius-(Radius_0*Radius_0))*Solution[0];
      Sponge_1 += (Sponge_Radius-(Radius_0*Radius_0))*Solution[1];
      
    }
    if (Radius_1 < Sponge_Radius) {
      Solution = node[Point_1]->GetSolution();
      Sponge_0 += (Sponge_Radius-(Radius_1*Radius_1))*Solution[0];
      Sponge_1 += (Sponge_Radius-(Radius_1*Radius_1))*Solution[1];
      
    }
    if (Radius_2 < Sponge_Radius) {
      Solution = node[Point_2]->GetSolution();
      Sponge_0 += (Sponge_Radius-(Radius_2*Radius_2))*Solution[0];
      Sponge_1 += (Sponge_Radius-(Radius_2*Radius_2))*Solution[1];
      
    }
    
    /* Source terms averaged over the elements */
    Res_Sour[0] = eps*(1.0/3.0)*Sponge_0;
    Res_Sour[1] = eps*(1.0/3.0)*Sponge_1;
    
    /* Set sponge source terms */
    LinSysRes.SubtractBlock(Point_0, Res_Sour);
    LinSysRes.SubtractBlock(Point_1, Res_Sour);
    LinSysRes.SubtractBlock(Point_2, Res_Sour);
    
	}
  
}


void CWaveSolver::Source_Template(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CNumerics *numerics,
                                    CConfig   *config,
                                    unsigned short iMesh) {

}

void CWaveSolver::Galerkin_Method(CGeometry *geometry, 
                                    CSolver **solver_container, 
                                    CNumerics *numerics,
                                    CConfig   *config, 
                                    unsigned short iMesh) {
  
  /* Local variables and initialization */
  
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, iPoint, total_index;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, wave_speed_2;
	unsigned short iVar;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
    
		Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
    
		numerics->SetCoord(Coord_0, Coord_1, Coord_2);
    
		numerics->ComputeResidual(StiffMatrix_Elem, config);
    
    /*--- Compute the square of the wave speed ---*/
		wave_speed_2 = config->GetWaveSpeed()*config->GetWaveSpeed();
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[0][0]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_0, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[0][1]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_0, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[0][2]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_0, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][0]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_1, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][1]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_1, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][2]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_1, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[2][0]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_2, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[2][1]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_2, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[2][2]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
    StiffMatrixSpace.AddBlock(Point_2, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
    
	}
	
  /* Prepare solution vector for multiplication */
  
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			LinSysSol[total_index] = node[iPoint]->GetSolution(iVar);
			LinSysRes[total_index] = 0.0;
		}
	
	StiffMatrixSpace.MatrixVectorProduct(LinSysSol, LinSysRes);
	
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			Residual[iVar] = LinSysRes[total_index];
		}
		LinSysRes.SubtractBlock(iPoint, Residual);
	}
  
  
  
}

void CWaveSolver::SetNoise_Source(CSolver ***flow_solution, CGeometry **wave_geometry, CConfig *wave_config) {
	unsigned short iDim, jDim, iMarker;
	unsigned long iVertex, iPoint, iPoint_Donor;
	double Thickness_Noise, Loading_Noise, Quadrupole_Noise;
  double Pressure, Velocity[3], Density, *Solution = NULL, *Solution_Old = NULL;
  double Density_Old, Velocity_Old[3], *Normal = NULL, U_n, U_nM1, L[3][3], Ln[3];
  
  /* Freestream Quantities */
  double deltaT = wave_config->GetDelta_UnstTimeND();
  double rho_0  = wave_config->GetDensity_FreeStreamND();
  double p_0    = wave_config->GetPressure_FreeStreamND();
  double *v_inf = wave_config->GetVelocity_FreeStreamND();
	bool incompressible = wave_config->GetIncompressible();

	for (iMarker = 0; iMarker < wave_config->GetnMarker_All(); iMarker++) {
    
		if (wave_config->GetMarker_All_Boundary(iMarker) == FWH_SURFACE) {
      
			for(iVertex = 0; iVertex < wave_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
        
				iPoint = wave_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
				iPoint_Donor = wave_geometry[MESH_0]->vertex[iMarker][iVertex]->GetDonorPoint();
        
        /* Some geometry information for this boundary node */
        Normal = wave_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNormal();
        
        double Area = 0.0; double UnitaryNormal[3];
        for (iDim = 0; iDim < nDim; iDim++)
          Area += Normal[iDim]*Normal[iDim];
        Area = sqrt (Area);
        
        /* Flipping the normal for now */
        for (iDim = 0; iDim < nDim; iDim++)
          UnitaryNormal[iDim] = -Normal[iDim]/Area;
        
        /* Get primitive variables from the CFD solution */
				Solution = flow_solution[MESH_0][FLOW_SOL]->node[iPoint_Donor]->GetSolution();
        Density  = Solution[0];
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Solution[iDim+1]/Density;
        Pressure = flow_solution[MESH_0][FLOW_SOL]->node[iPoint_Donor]->GetPressure(incompressible);

        /* Get old solution for computing time derivative */
        Solution_Old = flow_solution[MESH_0][FLOW_SOL]->node[iPoint_Donor]->GetSolution_time_n();
        Density_Old = Solution_Old[0];
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity_Old[iDim] = Solution_Old[iDim+1]/Density_Old;
        
        
        /* Compute thickness noise (need to check subtraction of mean flow...) */
        U_n = 0.0; U_nM1 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          // this version subtracts off mean flow
          //U_n   += ( (Velocity[iDim] - v_inf[iDim]) + (Density/rho_0 - 1.0)*(Velocity[iDim] - v_inf[iDim]))*(UnitaryNormal[iDim]*Area);
          //U_nM1 += ( (Velocity_Old[iDim]  - v_inf[iDim]) + (Density_Old/rho_0 - 1.0)*(Velocity_Old[iDim] - v_inf[iDim]))*(UnitaryNormal[iDim]*Area);
          U_n   += ( (Velocity[iDim]) + (Density/rho_0 - 1.0)*(Velocity[iDim] ))*(UnitaryNormal[iDim]*Area);
          U_nM1 += ( (Velocity_Old[iDim]) + (Density_Old/rho_0 - 1.0)*(Velocity_Old[iDim] ))*(UnitaryNormal[iDim]*Area);
        }
        
        /* Approximate the d/dt term with a backward difference */
        Thickness_Noise = rho_0*(U_n-U_nM1)/deltaT;
        
        /* Build loading noise L marix */
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            L[iDim][jDim] = Density*(Velocity[iDim] - v_inf[iDim])*(Velocity[iDim] - v_inf[iDim]);
            if (iDim == jDim)
              L[iDim][jDim] += (Pressure - p_0);
          }
        }
        
        /* Dot with normal vector */
        for (iDim = 0; iDim < nDim; iDim++) {
          Ln[iDim] = 0.0;
          for (jDim = 0; jDim < nDim; jDim++) {
            Ln[iDim] += L[iDim][jDim]*(UnitaryNormal[jDim]*Area);
        }
        }
        
        /* Need to take divergence of Ln */
          
          
        /* Leave as zero for now */

        Loading_Noise   = 0.0;
        Quadrupole_Noise = 0.0;

        /* Set source values at this point */
        
				node[iPoint]->SetThickness_Noise(Thickness_Noise);
        node[iPoint]->SetLoading_Noise(Loading_Noise);
        node[iPoint]->SetQuadrupole_Noise(Quadrupole_Noise);
			}
		}
	}
    
    /* Horribly inefficient - use aux var to take divergence of Ln? need to think about this */

	
}

void CWaveSolver::BC_Euler_Wall(CGeometry *geometry, 
                                  CSolver **solver_container, 
                                  CNumerics *numerics, 
                                  CConfig   *config, 
																	unsigned short val_marker) {
	
  /* Local variables */
  
  unsigned long iPoint, iVertex, total_index, iter;
	double deltaT, omega, time, ampl, wave_sol[2];
  
  /* Set the values needed for periodic forcing */
  
  deltaT = config->GetDelta_UnstTimeND();
  iter   = config->GetExtIter();  
  time   = static_cast<double>(iter)*deltaT;
  omega  = 1.0;
  ampl   = 1.0;
  
    /* Compute sin wave forcing at the boundary */
  
  wave_sol[0] = ampl*sin(omega*time);
  wave_sol[1] = 0.0; //-ampl*cos(omega*time_new)*omega;
  
  /* Set the solution at the boundary nodes and zero the residual */
	
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
			Solution[iVar] = wave_sol[iVar];
			Residual[iVar] = 0.0;
		}
		node[iPoint]->SetSolution(Solution);
		node[iPoint]->SetSolution_Old(Solution);
    LinSysRes.SetBlock(iPoint, Residual);
		LinSysRes.SetBlock(iPoint, Residual);
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      Jacobian.DeleteValsRowi(total_index);
    }
	}
  
}

void CWaveSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, 
                                 CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, 
                                 unsigned short val_marker) {
	
  
//  /* Turned off for now - try to make this a type of sponge */
//  
//  unsigned long iPoint, iVertex, total_index;
//
//	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
//			Solution[iVar] = 0.0;
//			Residual[iVar] = 0.0;
//		}
//    Solution[0] = 0.0;
//    Solution[1] = node[iPoint]->GetSolution(1);
//		node[iPoint]->SetSolution(Solution);
//		node[iPoint]->SetSolution_Old(Solution);
//    SetResidual(iPoint, Residual); 
//		SetResidual(iPoint, Residual);
//    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
//      total_index = iPoint*nVar+iVar;
//      Jacobian.DeleteValsRowi(total_index);
//    }
//	}
  
}

void CWaveSolver::BC_Observer(CGeometry *geometry, 
                                CSolver **solver_container, 
                                CNumerics *numerics, 
                                CConfig   *config, 
                                unsigned short val_marker) {
  
	bool adjoint = config->GetAdjoint();
  
  /*--- This boundary condition only applies to the adjoint problem ---*/
  
  if (adjoint) {
    
    /*--- Local Variables ---*/
    
    double a[3];
    unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
    double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
    double *Solution_0 = NULL, *Solution_1 = NULL, bcn = 0.0;
    double Length_Elem = 0.0;
    unsigned short iDim;

//    This is canceled out in the direct solution at the moment
//    /*--- Compute the constant term for the adjoint BC ---*/
//    double ONE_TWOC2 = 1.0/(2.0*config->GetWaveSpeed()*config->GetWaveSpeed());
    
    /*--- Loop over all surface elements and apply BC as a source term ---*/
    
    for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {		
      Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);	Coord_0 = geometry->node[Point_0]->GetCoord();
      Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);	Coord_1 = geometry->node[Point_1]->GetCoord();
      if (nDim == 3) { Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord();}
      
      
      /*--- Compute area (3D), and length of the surfaces (2D) ---*/
      if (nDim == 2) {
        for (iDim = 0; iDim < nDim; iDim++)
          a[iDim] = Coord_0[iDim]-Coord_1[iDim];
        Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);			
      }
      if (nDim == 3) {
        //        for (iDim = 0; iDim < nDim; iDim++) {
        //          a[iDim] = Coord_0[iDim]-Coord_2[iDim];
        //          b[iDim] = Coord_1[iDim]-Coord_2[iDim];
        //        }
        //        Area_Elem = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
      }	
      
      /*--- Compute adjoint boundary term for this edge from the direct 
       solution and the constant coefficient based on the wave speed 
       (could have a target density signature here as (rho' - rho_targ))---*/
      Solution_0 = node[Point_0]->GetSolution_Direct();
      Solution_1 = node[Point_1]->GetSolution_Direct();
      bcn = -0.5*(Solution_0[0]+Solution_1[0]);
      
      if (nDim == 2) {
        Residual[0] = 0.0; 
        Residual[1] = (1.0/2.0)*bcn*Length_Elem; 
        LinSysRes.AddBlock(Point_0, Residual);
        
        Residual[0] = 0.0;
        Residual[1] = (1.0/2.0)*bcn*Length_Elem;
        LinSysRes.AddBlock(Point_1, Residual);
      }
      else {
        /*--- do nothing for now - only in 2-D ---*/
      }
      
    }
  }
  
}

void CWaveSolver::Wave_Strength(CGeometry *geometry, CConfig *config) {
	
  unsigned long iPoint, iVertex;
  unsigned short iMarker, Boundary, Monitoring;
  double WaveSol, WaveStrength = 0.0, factor;
  
  /* Multiplying rho' by c^2 gives the acoustic pressure, p' */
  factor = config->GetWaveSpeed()*config->GetWaveSpeed();
  
  /* Initialize wave strength to zero */
  
  Total_CWave = 0.0; AllBound_CWave = 0.0;
  
  /*--- Loop over the wave observer markers ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary   = config->GetMarker_All_Boundary(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    if (Boundary == WAVE_OBSERVER) {
      
      WaveStrength = 0.0;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /* Wave solution from the first array position */
        WaveSol = node[iPoint]->GetSolution(0);
        
        /* Add contribution to total strength */
        WaveStrength += factor*WaveSol*WaveSol;
      }
      
      if  (Monitoring == YES) {
        CWave[iMarker] = WaveStrength;
        AllBound_CWave += CWave[iMarker];
      }
      
    }
  }
  
  Total_CWave += AllBound_CWave;
  
}

void CWaveSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep, 
                                        unsigned short iMesh, unsigned short RunTime_EqSystem) {
	
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double a[3], b[3], Area_Local, Time_Num;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
	unsigned short iDim, iVar, jVar;
	double TimeJac = 0.0;
	
	/*--- Numerical time step (this system is uncoditional stable... a very big number can be used) ---*/
	Time_Num = config->GetDelta_UnstTimeND();
	
	/*--- Loop through elements to compute contributions from the matrix
   blocks involving time. These contributions are also added to the 
   Jacobian w/ the time step. Spatial source terms are also computed. ---*/
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
    /* Get node numbers and their coordinate vectors */
    
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		
    Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		
    /* Compute triangle area (2-D) */
		
    for (iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
		Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

		
		/*----------------------------------------------------------------*/
		/*--- Block contributions to the Jacobian (includes time step) ---*/
		/*----------------------------------------------------------------*/
		
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				StiffMatrix_Node[iVar][jVar] = 0.0;	
		
		if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST) TimeJac = 1.0/Time_Num;
		if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) TimeJac = 3.0/(2.0*Time_Num);
		
		/*--- Diagonal value identity matrix ---*/
			StiffMatrix_Node[0][0] = 1.0*TimeJac; 

		
		/*--- Diagonal value ---*/
			StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local*TimeJac);

    Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_1, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_2, StiffMatrix_Node);
		
		/*--- Off Diagonal value ---*/
			StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local*TimeJac);

		Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_1, StiffMatrix_Node);
		Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_2, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_2, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_1, StiffMatrix_Node);

    
	}
	
	unsigned long iPoint, total_index;
	double *U_time_nM1, *U_time_n, *U_time_nP1;
	
	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 
		
		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n   = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();
		
		/*--- Compute Residual ---*/
		for(iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			LinSysRes[total_index] = 0.0;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				LinSysSol[total_index] = (U_time_nP1[iVar] - U_time_n[iVar]);
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				LinSysSol[total_index] = (U_time_nP1[iVar] - (4.0/3.0)*U_time_n[iVar] +  (1.0/3.0)*U_time_nM1[iVar]);
		}
	}
	
	StiffMatrixTime.MatrixVectorProduct(LinSysSol, LinSysRes);
  
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			Residual[iVar] = LinSysRes[total_index];
		}
		LinSysRes.SubtractBlock(iPoint, Residual);
	}
	
}


void CWaveSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	
    unsigned short iVar;
	unsigned long iPoint, total_index;
    
	/*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
    }
	
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
        
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			LinSysSol[total_index] = 0.0;
			AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
            AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex());
		}
	}
    
    /*--- Initialize residual and solution at the ghost points ---*/
    for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
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
    system.BCGSTAB(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                   config->GetLinear_Solver_Iter(), false);
  else if (config->GetKind_Linear_Solver() == FGMRES)
    system.FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                 config->GetLinear_Solver_Iter(), false);
  
  delete mat_vec;
  delete precond;
	
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar]);
		}
	}
	
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}


void CWaveSolver::SetSpace_Matrix(CGeometry *geometry, 
                                    CConfig   *config) {
  
//  /* Local variables and initialization */
//  
//	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, iPoint, total_index;
//	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, wave_speed_2;
//	unsigned short iVar;
//	
//	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
//    
//		Point_0 = geometry->elem[iElem]->GetNode(0);
//		Point_1 = geometry->elem[iElem]->GetNode(1);
//		Point_2 = geometry->elem[iElem]->GetNode(2);
//    
//		Coord_0 = geometry->node[Point_0]->GetCoord();
//		Coord_1 = geometry->node[Point_1]->GetCoord();
//		Coord_2 = geometry->node[Point_2]->GetCoord();
//    
//		numerics->SetCoord(Coord_0, Coord_1, Coord_2);
//    
//		numerics->ComputeResidual(StiffMatrix_Elem, config);
//    
//    /*--- Compute the square of the wave speed ---*/
//		wave_speed_2 = config->GetWaveSpeed()*config->GetWaveSpeed();
//    
//    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
//    StiffMatrix_Node[1][0] = StiffMatrix_Elem[0][0]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
//    StiffMatrixSpace.AddBlock(Point_0, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
//    
//    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
//    StiffMatrix_Node[1][0] = StiffMatrix_Elem[0][1]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
//    StiffMatrixSpace.AddBlock(Point_0, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
//    
//    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
//    StiffMatrix_Node[1][0] = StiffMatrix_Elem[0][2]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
//    StiffMatrixSpace.AddBlock(Point_0, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
//    
//    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
//    StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][0]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
//    StiffMatrixSpace.AddBlock(Point_1, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
//    
//    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
//    StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][1]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
//    StiffMatrixSpace.AddBlock(Point_1, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
//    
//    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
//    StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][2]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
//    StiffMatrixSpace.AddBlock(Point_1, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
//    
//    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
//    StiffMatrix_Node[1][0] = StiffMatrix_Elem[2][0]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
//    StiffMatrixSpace.AddBlock(Point_2, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
//    
//    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
//    StiffMatrix_Node[1][0] = StiffMatrix_Elem[2][1]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
//    StiffMatrixSpace.AddBlock(Point_2, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
//    
//    StiffMatrix_Node[0][0] = 0.0;                               StiffMatrix_Node[0][1] = -1.0;
//    StiffMatrix_Node[1][0] = StiffMatrix_Elem[2][2]*wave_speed_2; StiffMatrix_Node[1][1] = 0.0;
//    StiffMatrixSpace.AddBlock(Point_2, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
//    
//	}
  
}

void CWaveSolver::SetTime_Matrix(CGeometry *geometry, 
                                   CConfig   *config) {
  
  
  /* Local variables and initialization */
  
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double a[3], b[3], Area_Local, Time_Num, Time_Phys;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
	unsigned short iDim;
	
  /*--- Numerical time step. This system is unconditionally stable,
   so a very big step can be used. ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING) 
    Time_Num = config->GetDelta_UnstTimeND();
	else Time_Num = 1E+30;
  
  /*--- Physical timestep for source terms ---*/
  Time_Phys = config->GetDelta_UnstTimeND();
  
	/* Loop through elements to compute contributions from the matrix     */
  /* blocks involving time. These contributions are also added to the   */
  /* Jacobian w/ the time step. Spatial source terms are also computed. */
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
    /* Get node numbers and their coordinate vectors */
    
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		
    Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		
    /* Compute triangle area (2-D) */
		
    for (iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
		Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
		
    /*--- Block contributions to the Jacobian (includes time step) ---*/
    
		StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] =  1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] =  1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
    
    StiffMatrix_Node[0][0] = 1.0/Time_Num; StiffMatrix_Node[0][1] = 0.0;
    StiffMatrix_Node[1][0] = 0.0;            StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local/Time_Num);
    Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
    
  }
  
  //unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	//double a[3], b[3], Area_Local, Time_Num;
	//double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
	unsigned short iVar, jVar;
	double TimeJac = 0.0;
	
	/*--- Numerical time step (this system is uncoditional stable... a very big number can be used) ---*/
  if (config->GetUnsteady_Simulation() == TIME_STEPPING) 
    Time_Num = 1E+30;
	else Time_Num = config->GetDelta_UnstTimeND();
	
	/*--- Loop through elements to compute contributions from the matrix
   blocks involving time. These contributions are also added to the 
   Jacobian w/ the time step. Spatial source terms are also computed. ---*/
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
    /* Get node numbers and their coordinate vectors */
    
		Point_0 = geometry->elem[iElem]->GetNode(0);
		Point_1 = geometry->elem[iElem]->GetNode(1);
		Point_2 = geometry->elem[iElem]->GetNode(2);
		
    Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();
		Coord_2 = geometry->node[Point_2]->GetCoord();
		
    /* Compute triangle area (2-D) */
		
    for (iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
		Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
    
		
		/*----------------------------------------------------------------*/
		/*--- Block contributions to the Jacobian (includes time step) ---*/
		/*----------------------------------------------------------------*/
		
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				StiffMatrix_Node[iVar][jVar] = 0.0;	
		
		if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST) TimeJac = 1.0/Time_Num;
		if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) TimeJac = 3.0/(2.0*Time_Num);
		
		/*--- Diagonal value identity matrix ---*/
    StiffMatrix_Node[0][0] = 1.0*TimeJac; 
    
		
		/*--- Diagonal value ---*/
    StiffMatrix_Node[1][1] = (2.0/12.0)*(Area_Local*TimeJac);
    
    Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_1, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_2, StiffMatrix_Node);
		
		/*--- Off Diagonal value ---*/
    StiffMatrix_Node[1][1] = (1.0/12.0)*(Area_Local*TimeJac);
    
		Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_1, StiffMatrix_Node);
		Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_2, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_2, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_1, StiffMatrix_Node);
    
	}
}

void CWaveSolver::GetRestart(CGeometry *geometry, CConfig *config) {
  
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif
  
  /*--- Restart the solution from file information ---*/
  string restart_filename = config->GetRestart_WaveFileName();
  unsigned long iPoint, index, nFlowIter, adjIter, flowIter;
  char buffer[50];
  string UnstExt, text_line;
  ifstream restart_file;
  
  /*--- For the unsteady adjoint, we integrate backwards through
   physical time, so load in the direct solution files in reverse. ---*/  
  if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    nFlowIter = config->GetnExtIter();
    adjIter   = config->GetExtIter();
    flowIter  = nFlowIter - adjIter - 1;
    restart_filename.erase (restart_filename.end()-4, restart_filename.end());
    if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
    if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
    if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
    if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
    if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
    UnstExt = string(buffer);
    restart_filename.append(UnstExt);
  } else {
    flowIter  =config->GetExtIter();
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
  if (rank == MASTER_NODE)
    cout << "Reading in the direct wave solution from iteration " << flowIter << "." << endl;;
  restart_file.open(restart_filename.data(), ios::in);
  
  /*--- In case there is no file ---*/
  if (restart_file.fail()) {
    cout << "There is no wave restart file!!" << endl;
    cout << "Press any key to exit..." << endl;
    cin.get(); exit(1);
  }
  
  /*--- Read the restart file ---*/
  for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    getline(restart_file,text_line);
    istringstream point_line(text_line);
    point_line >> index >> Solution[0] >> Solution[1];
    node[iPoint]->SetSolution_Direct(Solution);
  }
  
  /*--- Close the restart file ---*/
  restart_file.close();
  
}
