
/*!
 * \file solution_direct_heat.cpp
 * \brief Main subrotuines for solving the heat equation.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.8
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/solver_structure.hpp"

CHeatSolver::CHeatSolver(void) : CSolver() { }

CHeatSolver::CHeatSolver(CGeometry *geometry, CConfig *config) : CSolver() {
  
	unsigned short nMarker, iVar;
  unsigned long iPoint;
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  nPoint =        geometry->GetnPoint();
  nPointDomain =  geometry->GetnPointDomain();
	nDim    =       geometry->GetnDim();
	nMarker =       config->GetnMarker_All();
	node    =       new CVariable*[nPoint];
	nVar    =       1;
  
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
  
	StiffMatrixSpace.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry);
	StiffMatrixTime.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry);
  if (rank == MASTER_NODE) cout << "Initialize jacobian structure (Linear Elasticity)." << endl;
	Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry);
  
  /*--- Initialization of linear solver structures ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Heat coefficient for all of the markers ---*/
  
  CHeat = new double[config->GetnMarker_All()];
  Total_CHeat = 0.0;
  
  /*--- Check for a restart (not really used), initialize from zero otherwise ---*/
  
	bool restart = (config->GetRestart());
	if (!restart) {
    
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
      
      /*--- Zero initial condition for testing source terms & forcing BCs ---*/
      
      Solution[0] = 0.0;
      Solution[1] = 0.0;
      
			node[iPoint] = new CHeatVariable(Solution, nDim, nVar, config);
      
      /*--- Copy solution to old containers if using dual time ---*/
      
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST) {
        node[iPoint]->Set_Solution_time_n();
      } else if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) {
        node[iPoint]->Set_Solution_time_n();
        node[iPoint]->Set_Solution_time_n1();
      }
      
    }
	} else {
    
    cout << "Heat restart file not currently configured!!" << endl;
    exit(1);
    
		string mesh_filename = config->GetSolution_FlowFileName();
		ifstream restart_file;
    
		char *cstr; cstr = new char [mesh_filename.size()+1];
		strcpy (cstr, mesh_filename.c_str());
		restart_file.open(cstr, ios::in);
    
		if (restart_file.fail()) {
			cout << "There is no Heat restart file!!" << endl;
			exit(1);
		}
		unsigned long index;
		string text_line;
    
		for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			point_line >> index >> Solution[0] >> Solution[1];
			node[iPoint] = new CHeatVariable(Solution, nDim, nVar, config);
		}
		restart_file.close();
	}
  
}

CHeatSolver::~CHeatSolver(void) {
  
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

void CHeatSolver::Preprocessing(CGeometry *geometry,
                                CSolver **solver_container,
                                CConfig   *config,
                                unsigned short iMesh,
                                unsigned short iRKStep,
                                unsigned short RunTime_EqSystem) {
  
  unsigned long iPoint;
  
  /*--- Set residuals and matrix entries to zero ---*/
  
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    LinSysSol.SetBlock_Zero(iPoint);
    LinSysAux.SetBlock_Zero(iPoint);
    LinSysRes.SetBlock_Zero(iPoint);
  }
	
  /*--- Zero out the entries in the various matrices ---*/
	StiffMatrixSpace.SetValZero();
	StiffMatrixTime.SetValZero();
	Jacobian.SetValZero();
  
}

void CHeatSolver::Source_Residual(CGeometry *geometry,
                                  CSolver **solver_container,
                                  CNumerics *numerics, CNumerics *second_numerics,
                                  CConfig   *config,
                                  unsigned short iMesh) {

  if (config->GetUnsteady_Simulation() != STEADY) {

    unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
    double a[3], b[3], c[3], d[3], Area_Local = 0.0, Volume_Local = 0.0, Time_Num;
    double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
    unsigned short iDim;

    /*--- Numerical time step (this system is uncoditional stable... a very big number can be used) ---*/
    if (config->GetUnsteady_Simulation() == TIME_STEPPING) Time_Num = config->GetDelta_UnstTimeND();
    else Time_Num = 1E30;
    
    /*--- Loop through elements to compute contributions from the matrix
     blocks involving time. These contributions are also added to the
     Jacobian w/ the time step. Spatial source terms are also computed. ---*/
    
    for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      
      /*--- Get node numbers and their coordinate vectors ---*/
      
      Point_0 = geometry->elem[iElem]->GetNode(0);	Coord_0 = geometry->node[Point_0]->GetCoord();
      Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1 = geometry->node[Point_1]->GetCoord();
      Point_2 = geometry->elem[iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord();
      if (nDim == 3) { Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3 = geometry->node[Point_3]->GetCoord(); }
      
      /*--- Compute area and volume ---*/
      
      if (nDim == 2) {
        for (iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = Coord_0[iDim]-Coord_2[iDim];
          b[iDim] = Coord_1[iDim]-Coord_2[iDim];
        }
        Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
      }
      else {
        for (iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = Coord_0[iDim]-Coord_2[iDim];
          b[iDim] = Coord_1[iDim]-Coord_2[iDim];
          c[iDim] = Coord_3[iDim]-Coord_2[iDim];
        }
        d[0] = a[1]*b[2]-a[2]*b[1];
        d[1] = -(a[0]*b[2]-a[2]*b[0]);
        d[2] = a[0]*b[1]-a[1]*b[0];
        Volume_Local = fabs(c[0]*d[0] + c[1]*d[1] + c[2]*d[2])/6.0;
      }
      
      /*--- Block contributions to the Jacobian (includes time step) ---*/
      
      if (nDim == 2) { StiffMatrix_Node[0][0] = (2.0/12.0)*(Area_Local/Time_Num); }
      else { StiffMatrix_Node[0][0] = (2.0/20.0)*(Volume_Local/Time_Num); }
      Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
      Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
      Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
      if (nDim == 3) Jacobian.AddBlock(Point_3, Point_3, StiffMatrix_Node);
      
      if (nDim == 2) { StiffMatrix_Node[0][0] = (1.0/12.0)*(Area_Local/Time_Num); }
      else { StiffMatrix_Node[0][0] = (1.0/20.0)*(Volume_Local/Time_Num); }
      
      Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
      Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
      Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
      Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
      Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
      Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
      if (nDim == 3) {
        Jacobian.AddBlock(Point_0, Point_3, StiffMatrix_Node);
        Jacobian.AddBlock(Point_1, Point_3, StiffMatrix_Node);
        Jacobian.AddBlock(Point_2, Point_3, StiffMatrix_Node);
        Jacobian.AddBlock(Point_3, Point_0, StiffMatrix_Node);
        Jacobian.AddBlock(Point_3, Point_1, StiffMatrix_Node);
        Jacobian.AddBlock(Point_3, Point_2, StiffMatrix_Node);
      }
      
    }
    
  }
  
}

void CHeatSolver::Galerkin_Method(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                  CConfig *config, unsigned short iMesh) {
  
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0, total_index, iPoint;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3 = NULL, Thermal_Diffusivity;
  
  Thermal_Diffusivity  = -config->GetThermalDiffusivity();

	if (nDim == 2 ) {
    
		for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      
      Point_0 = geometry->elem[iElem]->GetNode(0);  Coord_0 = geometry->node[Point_0]->GetCoord();
      Point_1 = geometry->elem[iElem]->GetNode(1);  Coord_1 = geometry->node[Point_1]->GetCoord();
      Point_2 = geometry->elem[iElem]->GetNode(2);  Coord_2 = geometry->node[Point_2]->GetCoord();
      
      numerics->SetCoord(Coord_0, Coord_1, Coord_2);
      numerics->ComputeResidual(StiffMatrix_Elem, config);
      
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]; StiffMatrixSpace.AddBlock(Point_0, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]; StiffMatrixSpace.AddBlock(Point_0, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]; StiffMatrixSpace.AddBlock(Point_0, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]; StiffMatrixSpace.AddBlock(Point_1, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]; StiffMatrixSpace.AddBlock(Point_1, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]; StiffMatrixSpace.AddBlock(Point_1, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]; StiffMatrixSpace.AddBlock(Point_2, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]; StiffMatrixSpace.AddBlock(Point_2, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]; StiffMatrixSpace.AddBlock(Point_2, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
      
		}
	}
  
	if (nDim == 3 ) {
    
		for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      
      Point_0 = geometry->elem[iElem]->GetNode(0); 	Coord_0 = geometry->node[Point_0]->GetCoord();
      Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1 = geometry->node[Point_1]->GetCoord();
      Point_2 = geometry->elem[iElem]->GetNode(2); 	Coord_2 = geometry->node[Point_2]->GetCoord();
      Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3 = geometry->node[Point_3]->GetCoord();
      
      numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
      numerics->ComputeResidual(StiffMatrix_Elem, config);
      
      Thermal_Diffusivity  = config->GetThermalDiffusivity();
      
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]; StiffMatrixSpace.AddBlock(Point_0, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]; StiffMatrixSpace.AddBlock(Point_0, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]; StiffMatrixSpace.AddBlock(Point_0, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][3]; StiffMatrixSpace.AddBlock(Point_0, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_3, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]; StiffMatrixSpace.AddBlock(Point_1, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]; StiffMatrixSpace.AddBlock(Point_1, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]; StiffMatrixSpace.AddBlock(Point_1, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][3]; StiffMatrixSpace.AddBlock(Point_1, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_3, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]; StiffMatrixSpace.AddBlock(Point_2, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]; StiffMatrixSpace.AddBlock(Point_2, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]; StiffMatrixSpace.AddBlock(Point_2, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][3]; StiffMatrixSpace.AddBlock(Point_2, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_3, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][0]; StiffMatrixSpace.AddBlock(Point_3, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_3, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][1]; StiffMatrixSpace.AddBlock(Point_3, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_3, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][2]; StiffMatrixSpace.AddBlock(Point_3, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_3, Point_2, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][3]; StiffMatrixSpace.AddBlock(Point_3, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_3, Point_3, StiffMatrix_Node);
      
		}
	}
  
  if (config->GetUnsteady_Simulation() != STEADY) {
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      total_index = iPoint*nVar;
      LinSysSol[total_index] = node[iPoint]->GetSolution(0);
      LinSysAux[total_index] = 0.0;
    }
    
    StiffMatrixSpace.MatrixVectorProduct(LinSysSol, LinSysAux);
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      total_index = iPoint*nVar;
      Residual[0] = LinSysAux[total_index];
      LinSysRes.SubtractBlock(iPoint, Residual);
    }
  }
  
}

void CHeatSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

void CHeatSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex, total_index;
	double Twall;
  
  /*--- Identify the boundary ---*/
  
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
  
	/*--- Retrieve the specified wall temperature ---*/
  
	Twall = config->GetIsothermal_Temperature(Marker_Tag);
  
  /*--- Set the solution at the boundary nodes and zero the residual ---*/
	
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    Solution[0] = Twall;
    
    node[iPoint]->SetSolution(Solution);
    node[iPoint]->SetSolution_Old(Solution);

    /*--- Unsteady solution, the equation is solved in terms of increments ---*/
    if (config->GetUnsteady_Simulation() != STEADY) Residual[0] = 0.0;
    
    LinSysRes.SetBlock(iPoint, Residual);
    LinSysSol.SetBlock(iPoint, Residual);

    total_index = iPoint*nVar;
    Jacobian.DeleteValsRowi(total_index);
    
  }
  
}


void CHeatSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep,
                                       unsigned short iMesh, unsigned short RunTime_EqSystem) {
	
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
	double a[3], b[3], c[3], d[3], Area_Local = 0.0, Volume_Local = 0.0, Time_Num;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
	unsigned short iDim, iVar;
	double TimeJac = 0.0;
	
	/*--- Numerical time step (this system is uncoditional stable... a very big number can be used) ---*/
	Time_Num = config->GetDelta_UnstTimeND();
	
	/*--- Loop through elements to compute contributions from the matrix
   blocks involving time. These contributions are also added to the
   Jacobian w/ the time step. Spatial source terms are also computed. ---*/
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
    /*--- Get node numbers and their coordinate vectors ---*/
    
		Point_0 = geometry->elem[iElem]->GetNode(0);	Coord_0 = geometry->node[Point_0]->GetCoord();
		Point_1 = geometry->elem[iElem]->GetNode(1);	Coord_1 = geometry->node[Point_1]->GetCoord();
		Point_2 = geometry->elem[iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord();
		if (nDim == 3) { Point_3 = geometry->elem[iElem]->GetNode(3);	Coord_3 = geometry->node[Point_3]->GetCoord(); }
		
    /*--- Compute area and volume ---*/

		if (nDim == 2) {
			for (iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}
			Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
		}
		else {
			for (iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
				c[iDim] = Coord_3[iDim]-Coord_2[iDim];
			}
			d[0] = a[1]*b[2]-a[2]*b[1]; d[1] = -(a[0]*b[2]-a[2]*b[0]); d[2] = a[0]*b[1]-a[1]*b[0];
			Volume_Local = fabs(c[0]*d[0] + c[1]*d[1] + c[2]*d[2])/6.0;
		}
		
		/*--- Block contributions to the Jacobian (includes time step) ---*/
				
		if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST) TimeJac = 1.0/Time_Num;
		if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) TimeJac = 3.0/(2.0*Time_Num);
		  
		if (nDim == 2) { StiffMatrix_Node[0][0] = (2.0/12.0)*(Area_Local*TimeJac); }
		else { StiffMatrix_Node[0][0] = (2.0/20.0)*(Volume_Local*TimeJac); }
    
    Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_1, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_2, StiffMatrix_Node);
		if (nDim == 3) { Jacobian.AddBlock(Point_3, Point_3, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_2, StiffMatrix_Node); }
		  
		if (nDim == 2) { StiffMatrix_Node[0][0] = (1.0/12.0)*(Area_Local*TimeJac); }
		else { StiffMatrix_Node[0][0] = (1.0/20.0)*(Volume_Local*TimeJac); }
    
		Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_1, StiffMatrix_Node);
		Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_2, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_2, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_0, StiffMatrix_Node);
		Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_1, StiffMatrix_Node);
		if (nDim == 3) {
			Jacobian.AddBlock(Point_0, Point_3, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_3, StiffMatrix_Node);
			Jacobian.AddBlock(Point_1, Point_3, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_3, StiffMatrix_Node);
			Jacobian.AddBlock(Point_2, Point_3, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_3, StiffMatrix_Node);
			Jacobian.AddBlock(Point_3, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_3, Point_0, StiffMatrix_Node);
			Jacobian.AddBlock(Point_3, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_3, Point_1, StiffMatrix_Node);
			Jacobian.AddBlock(Point_3, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_3, Point_2, StiffMatrix_Node);
		}
    
	}
	
	unsigned long iPoint, total_index;
	double *U_time_nM1, *U_time_n, *U_time_nP1;
	
	/*--- loop over points ---*/
  
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		
		/*--- Solution at time n-1, n and n+1 ---*/
    
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n   = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();
		
		/*--- Compute residual ---*/
    
		for(iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				LinSysSol[total_index] = ( U_time_nP1[iVar] - U_time_n[iVar] );
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				LinSysSol[total_index] = ( U_time_nP1[iVar] - (4.0/3.0)*U_time_n[iVar] + (1.0/3.0)*U_time_nM1[iVar] );
		}
	}
	
  /*--- Contribution to the residual ---*/
  
	StiffMatrixTime.MatrixVectorProduct(LinSysSol, LinSysAux);
  
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    total_index = iPoint*nVar;
    Residual[0] = LinSysAux[total_index];
		LinSysRes.SubtractBlock(iPoint, Residual);
	}
	
}

void CHeatSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iVar;
	unsigned long iPoint, total_index;
  bool output;
	
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			LinSysSol[total_index] = 0.0;
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
  
  if (config->GetUnsteady_Simulation() == STEADY) output = true;
  else output = false;
  
  if (config->GetKind_Linear_Solver() == BCGSTAB)
    system.BCGSTAB(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(), config->GetLinear_Solver_Iter(), output);
  else if (config->GetKind_Linear_Solver() == FGMRES)
    system.FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(), config->GetLinear_Solver_Iter(), output);
  
  delete mat_vec;
  delete precond;
  
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			if (config->GetUnsteady_Simulation() == STEADY) node[iPoint]->SetSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      else node[iPoint]->AddSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
		}
	}
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*---  Compute the residual Ax-f ---*/
	Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);
  
  /*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Compute the residual ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			total_index = iPoint*nVar+iVar;
			AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
      AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex());
		}
	}
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}
