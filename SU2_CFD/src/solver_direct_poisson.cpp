/*!
 * \file solution_direct_poisson.cpp
 * \brief Main subrotuines for solving direct problems
 * \author F. Palacios
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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
#include "../include/variables/CPBIncEulerVariable.hpp"

CPoissonSolverFVM::CPoissonSolverFVM(void) : CSolver() { }

CPoissonSolverFVM::CPoissonSolverFVM(CGeometry *geometry, CConfig *config) : CSolver() {
  
  unsigned long  iPoint;
  unsigned short iVar, iDim;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  nDim =          geometry->GetnDim();
  nPoint =        geometry->GetnPoint();
  nPointDomain =  geometry->GetnPointDomain();
  nVar =          1;
  node =          new CVariable*[nPoint];
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  Residual = new su2double[nVar]; Residual_RMS = new su2double[nVar];
  Solution = new su2double[nVar];
  Residual_Max = new su2double[nVar];
  

  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  
  
  /*--- Define some auxiliar vector related with the solution ---*/

  Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];

  /*--- Define some auxiliary vectors related to the geometry ---*/

  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;

  
 
 
  
 if (config->GetKind_TimeIntScheme_Poisson() == EULER_IMPLICIT) { 
	 Jacobian_i = new su2double* [nVar];
	 Jacobian_j = new su2double* [nVar];
	 for (iVar = 0; iVar < nVar; iVar++) {
		 Jacobian_i[iVar] = new su2double [nVar];
		 Jacobian_j[iVar] = new su2double [nVar];
	 }
	 /*--- Initialization of the structure of the whole Jacobian ---*/
	 if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Poisson equation)." << endl;
	 Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
 }
  /*--- Solution and residual vectors ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Computation of gradients by least squares ---*/
  
  Smatrix = new su2double* [nDim]; // S matrix := inv(R)*traspose(inv(R))
  for (iDim = 0; iDim < nDim; iDim++)
    Smatrix[iDim] = new su2double [nDim];
  
  Cvector = new su2double* [nVar]; // c vector := transpose(WA)*(Wb)
  for (iVar = 0; iVar < nVar; iVar++)
    Cvector[iVar] = new su2double [nDim];
  
  /*--- Always instantiate and initialize the variable to a zero value. ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    if (config->GetKind_Incomp_System()==PRESSURE_BASED) 
        node[iPoint] = new CPoissonVariable(0.0, nDim, nVar, config);
    else 
        node[iPoint] = new CPoissonVariable(0.0, nDim, nVar, config);
  
  /*--- Perform the MPI communication of the solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
}

CPoissonSolverFVM::~CPoissonSolverFVM(void) {
  
  unsigned int iVar;
  iVar = 1;
}


void CPoissonSolverFVM::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {
}


void CPoissonSolverFVM::Preprocessing(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Initialize the residual vector ---*/
    LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  /*--- Initialize the Jacobian matrices ---*/

  Jacobian.SetValZero();

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    
}




void CPoissonSolverFVM:: SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {

}


void CPoissonSolverFVM:: Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) {


}


void CPoissonSolverFVM:: Set_MPI_Solution(CGeometry *geometry, CConfig *config) {


}

void CPoissonSolverFVM:: Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {



}

void CPoissonSolverFVM:: Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config){


}

void CPoissonSolverFVM:: LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo){

/*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh, iMeshFine;
  unsigned long iPoint, index, iChildren, Point_Fine;
  
  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine;
  bool grid_movement  = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool steady_restart = config->GetSteadyRestart();
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;

  string UnstExt, text_line;
  ifstream restart_file;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  string restart_filename = config->GetSolution_FlowFileName();

  Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = 0;

  skipVars = geometry[MESH_0]->GetnDim();;

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    restart_filename = config->GetMultizone_FileName(restart_filename, iZone);

  /*--- Modify file name for an unsteady restart ---*/

  if (dual_time || time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      node[iPoint_Local]->SetSolution(Solution);
      iPoint_Global_Local++;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
    if (rank == MASTER_NODE) {
      cout << endl << "The solution file " << restart_filename.data() << " doesn't match with the mesh file!" << endl;
      cout << "It could be empty lines at the end of the file." << endl << endl;
    }
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. We alo call the preprocessing routine
   on the fine level in order to have all necessary quantities updated,
   especially if this is a turbulent simulation (eddy viscosity). ---*/

  //solver[MESH_0][POISSON_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][POISSON_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_POISSON_SYS, false);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][HEAT_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][POISSON_SOL]->node[iPoint]->SetSolution(Solution);
    }
    //solver[iMesh][POISSON_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    solver[iMesh][POISSON_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_POISSON_SYS, false);
  }

  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;




}

void CPoissonSolverFVM::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
										 
su2double Poisson_Coeff_i,Poisson_Coeff_j,**Sol_i_Grad,**Sol_j_Grad,Poissonval_i,Poissonval_j,Normal[3];
su2double Mom_Coeff_i[3],Mom_Coeff_j[3], Vol_i, delT_i, Vol_j, delT_j;
unsigned long iEdge, iPoint, jPoint;
unsigned short iDim;
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		
		/*--- Points coordinates, and normal vector ---*/
		numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
						geometry->node[jPoint]->GetCoord());
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
		
		numerics->SetVolume(geometry->node[iPoint]->GetVolume());
		numerics->SetVolume(geometry->node[jPoint]->GetVolume());
		
		/*--- Primitive variables w/o reconstruction ---*/
		Poissonval_i = node[iPoint]->GetSolution(0);
		Poissonval_j = node[jPoint]->GetSolution(0);
			
		numerics->SetPoissonval(Poissonval_i,Poissonval_j);
			
		Sol_i_Grad = node[iPoint]->GetGradient();
		Sol_j_Grad = node[jPoint]->GetGradient();
    
		numerics->SetConsVarGradient(Sol_i_Grad, Sol_j_Grad);

		if (config->GetKind_Incomp_System()!=PRESSURE_BASED) {
			for (iDim = 0; iDim < nDim; iDim++) {
				Mom_Coeff_i[iDim] = 1.0;
			    Mom_Coeff_j[iDim] = 1.0;
			}
			numerics->SetInvMomCoeff(Mom_Coeff_i,Mom_Coeff_j);
		}
		else {
			Vol_i = geometry->node[iPoint]->GetVolume();
			delT_i = node[iPoint]->GetDelta_Time();
	  
			Vol_j = geometry->node[jPoint]->GetVolume();
			delT_j = node[jPoint]->GetDelta_Time();
			
			
			for (iDim = 0; iDim < nDim; iDim++) {
				Mom_Coeff_i[iDim] = solver_container[FLOW_SOL]->node[iPoint]->Get_Mom_Coeff(iDim) ;
			    Mom_Coeff_j[iDim] = solver_container[FLOW_SOL]->node[jPoint]->Get_Mom_Coeff(iDim) ;
			}
			numerics->SetInvMomCoeff(Mom_Coeff_i,Mom_Coeff_j);
		}

		/*--- Compute residual, and Jacobians ---*/
		numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

		/*--- Add and subtract residual, and update Jacobians ---*/
		LinSysRes.SubtractBlock(iPoint, Residual);
		LinSysRes.AddBlock(jPoint, Residual);
		
       if (config->GetKind_TimeIntScheme_Poisson() == EULER_IMPLICIT) {
		Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
		Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
		Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
		Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
	  }	  
	}
	
}

void CPoissonSolverFVM::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                        unsigned short iMesh){
							
 /*--- Compute gradients so we can use it to find the velocity corrections ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);							
  
  	
}

void CPoissonSolverFVM::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) {

  unsigned short iVar;
  unsigned long iPoint;
  su2double Src_Term;
  //ofstream Coarse_GridFile, Fine_GridFile;

  /*--- Initialize the source residual to zero ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
	  Residual[iVar] = 0.0;
  }
  
  /*if (iMesh == MESH_0 + 1) {
	  Coarse_GridFile.open("coarse_src.txt", ios::out);
  		//cout<<"Writing the _src.txt file."<<endl;
  }*/
  //if (iMesh == MESH_0 ) 
  //Fine_GridFile.open("fine_src.txt", ios::out);
  /*MPI_File Fine_GridFile;
  MPI_Status status;
  MPI_File_open(MPI_COMM_SELF, "Poisson_src.txt",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&Fine_GridFile);*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Load the volume of the dual mesh cell ---*/

    numerics->SetVolume(geometry->node[iPoint]->GetVolume());

    /*--- Compute the source term ---*/
        
    if ((config->GetKind_Incomp_System() == PRESSURE_BASED)) {
		
		Src_Term = solver_container[FLOW_SOL]->node[iPoint]->GetMassFlux() ;
		
		if (Src_Term != Src_Term)Src_Term = 0.0;
		
		if (iMesh == MESH_0 + 1) Src_Term += solver_container[FLOW_SOL]->node[iPoint]->GetMassTruncError();
		
		node[iPoint]->SetSourceTerm(Src_Term);
    }
    numerics->SetSourcePoisson(node[iPoint]->GetSourceTerm());
    
    /*--- Compute the source residual ---*/
   
    numerics->ComputeResidual(Residual, Jacobian_i, config);

    /*--- Add the source residual to the total ---*/

    LinSysRes.AddBlock(iPoint, Residual);
    
   //MPI_File_write(Fine_GridFile, &Src_Term, sizeof(double), MPI_DOUBLE, &status);
   
    //Source term is constant ==> jacobian is zero
    
   /* if (iMesh == MESH_0) {*/
		//Fine_GridFile<<iPoint<<"\t"<<geometry->node[iPoint]->GetCoord(0)<<"\t"<<geometry->node[iPoint]->GetCoord(1)<<"\t";
        //Fine_GridFile<<Src_Term<<"\t"<<solver_container[FLOW_SOL]->node[iPoint]->GetMassFlux()<<endl;
    /*}
    if (iMesh == MESH_0 + 1) {
		Coarse_GridFile<<iPoint<<"\t"<<geometry->node[iPoint]->GetCoord(0)<<"\t"<<geometry->node[iPoint]->GetCoord(1)<<"\t";
        Coarse_GridFile<<Src_Term<<"\t"<<solver_container[FLOW_SOL]->node[iPoint]->GetMassFlux()<<"\t"<<solver_container[FLOW_SOL]->node[iPoint]->GetMassTruncError()<<endl;
    }*/
  }
   //MPI_File_close(&Fine_GridFile);
  /*if (iMesh == MESH_0 + 1) Coarse_GridFile.close();*/
  //if (iMesh == MESH_0) 
  //Fine_GridFile.close();
}

void CPoissonSolverFVM::AssembleCoeffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {



}
void CPoissonSolverFVM::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long iPoint, total_index, IterLinSol = 0;;
  unsigned short iVar;
  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  
	/*--- Build implicit system ---*/

	/*--- Set maximum residual to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
	/*--- Initialize residual and solution at the ghost points ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
	  
	 /*--- Read the residual ---*/
     local_Res_TruncError = node[iPoint]->GetResTruncError();

	/*--- Read the volume ---*/

    Vol = geometry->node[iPoint]->GetVolume();

	/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      //LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar] );
      LinSysRes[total_index] = - (LinSysRes[total_index]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
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
  
  /*--- Solve or smooth the linear system ---*/
  
  IterLinSol = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      node[iPoint]->AddSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
     }
  }


  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);
  
}


void CPoissonSolverFVM::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {



}

void CPoissonSolverFVM::Direct_Solve(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  
}

void CPoissonSolverFVM::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                        unsigned short iMesh, unsigned long Iteration){
							
   
}


su2double CPoissonSolverFVM::GetDirichlet_BC(CGeometry *geometry, CConfig *config, unsigned long Point){
	
	su2double dirichlet_bc;
	
	if (config->GetKind_Incomp_System() == PRESSURE_BASED ) 
	    dirichlet_bc = 0.0;
	else 
	    dirichlet_bc = sin(geometry->node[Point]->GetCoord(0))*cos(geometry->node[Point]->GetCoord(1));
	
	return dirichlet_bc;
	
}


void CPoissonSolverFVM::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container,
                                  CConfig *config, unsigned short val_marker) {
  unsigned long Point, iVertex;
  su2double val_res,pi;
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  unsigned short iVar = 0;
  pi = 4.0*atan(1.0);
  

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    Point = geometry->vertex[val_marker][iVertex]->GetNode();
   
    Solution[0] = GetDirichlet_BC(geometry,config,Point);
    
    
 /*--- Assign the dirichlet BC value to the solution ---*/
    node[Point]->SetSolution(Solution);
    node[Point]->Set_OldSolution();
    
	if (config->GetKind_TimeIntScheme_Poisson()==EULER_IMPLICIT) {
		Jacobian.DeleteValsRowi(Point);
	}
    LinSysRes.SetBlock_Zero(Point, iVar);
    //if (config->GetnMGLevels() > 0) node[Point]->SetVal_ResTruncError_Zero(iVar);
    LinSysSol.SetBlock(Point, Solution);
  }

  
  
  
}


void CPoissonSolverFVM::BC_Neumann(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                unsigned short val_marker) { 
									
									
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double NeumannFlux, Area, *Normal,*Res_Visc;

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  NeumannFlux = 0.0;

  Res_Visc = new su2double[nVar];

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      Res_Visc[0] = NeumannFlux * Area;

      /*--- Viscous contribution to the residual at the wall ---*/

      LinSysRes.SubtractBlock(iPoint, Res_Visc);
    }

   }
}


void CPoissonSolverFVM::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

su2double Poisson_Coeff_i, Poissonval_i;
su2double Velocity_i[3], MassFlux_Part, small=1E-6;
unsigned long iVertex, iPoint, jPoint, Point_Normal;
unsigned short iDim, iVar;
su2double *Normal = new su2double[nDim];
su2double Coeff_Mean;
su2double *MomCoeffxNormal = new su2double[nDim];
su2double dist_ij_2, proj_vector_ij, Edge_Vector[3];

for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      /*--- The farfield boundary is considered as an inlet-outlet boundary, where flow 
       * can either enter or leave. For pressure, it is treated as a fully developed flow
       * and a dirichlet BC is applied. For velocity, based on the sign of massflux, either 
       * a dirichlet or a neumann BC is applied (in Flow_Correction routine). ---*/		
       
       geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
       
       
       for (iVar = 0; iVar < nVar; iVar++) {
		   Residual[iVar] = 0.0;
	   }
       
		/*MassFlux_Part = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
           MassFlux_Part -= solver_container[FLOW_SOL]->node[iPoint]->GetDensity()*(solver_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim))*Normal[iDim];
		} 
		if ((MassFlux_Part < 0.0) && (fabs(MassFlux_Part) > EPS)) {
			 LinSysRes.SubtractBlock(iPoint, Residual);
				 
		    if (config->GetKind_TimeIntScheme_Poisson() == EULER_IMPLICIT) {
			   Jacobian_i[0][0] = 0.0;
			   Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
		    }			    
		}*/
		//else {
			 for (iVar = 0; iVar < nVar; iVar++) {
			   LinSysRes.SetBlock_Zero(iPoint, iVar);
			   Residual[iVar] = 0.0;
		      }
		     
		      node[iPoint]->SetSolution(Residual);
		      node[iPoint]->Set_OldSolution();
		     
		     if (config->GetKind_TimeIntScheme_Poisson()==EULER_IMPLICIT) {
				 Jacobian.DeleteValsRowi(iPoint);
		     }
		//}
   }
  }
   delete [] Normal;
}
                                
                                
                                
void CPoissonSolverFVM::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config, unsigned short val_marker) {
su2double Poisson_Coeff_i,**Sol_i_Grad,Poissonval_i;
su2double Mom_Coeff_i[3],Proj_Mean_GradPoissonVar_Normal[3];
unsigned long iVertex, iPoint, jPoint;
unsigned short iDim, iVar;
su2double *Normal = new su2double[nDim];

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
     /*--- Zero flux (Neumann) BC on pressure ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
      }

	 /*--- Add and subtract residual, and update Jacobians ---*/
		LinSysRes.SubtractBlock(iPoint, Residual);
		
       if (config->GetKind_TimeIntScheme_Poisson() == EULER_IMPLICIT) {
		   Jacobian_i[0][0] = 0.0;
		Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
	  }
     }
	}
	/*--- Free locally allocated memory ---*/
  delete [] Normal;
}
                                 
                                 
void CPoissonSolverFVM::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
								
su2double Poisson_Coeff_i,**Sol_i_Grad,Poissonval_i;
su2double Mom_Coeff_i[3],Proj_Mean_GradPoissonVar_Normal[3];
unsigned long iVertex, iPoint, jPoint, total_index;
unsigned short iDim, iVar;
string Marker_Tag  = config->GetMarker_All_TagBound(val_marker);
su2double *Normal = new su2double[nDim];

/*--- Only fixed velocity inlet is considered ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
     /*--- Zero flux (Neumann) BC on pressure ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
      }
      
      /*--- Add and subtract residual, and update Jacobians ---*/
		LinSysRes.SubtractBlock(iPoint, Residual);
		
       if (config->GetKind_TimeIntScheme_Poisson() == EULER_IMPLICIT) {
		   Jacobian_i[0][0] = 0.0;
		Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
	  }
     }
	}
	/*--- Free locally allocated memory ---*/
  delete [] Normal;
}


void CPoissonSolverFVM::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
								
su2double Poisson_Coeff_i,**Sol_i_Grad,Poissonval_i;
su2double Mom_Coeff_i[3],Proj_Mean_GradPoissonVar_Normal[3];
unsigned long iVertex, iPoint, jPoint;
unsigned short iDim, iVar;
su2double Velocity_i[3], MassFlux_Part;
su2double *Normal = new su2double[nDim];
string Marker_Tag  = config->GetMarker_All_TagBound(val_marker);
unsigned short Kind_Outlet = config->GetKind_Inc_Outlet(Marker_Tag);
/*--- Only fully developed case is considered ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negative for outward convention) ---*/
            
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

     /*--- A Dirichlet BC is applied at the fully developed outlet, pressure is
      * assumed to be uniform and assigned the value of P_ref. ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
      }
      
      switch (Kind_Outlet) {
		  
		  case PRESSURE_OUTLET:
		     for (iVar = 0; iVar < nVar; iVar++)
		       LinSysRes.SetBlock_Zero(iPoint, iVar);
		     
		     node[iPoint]->SetSolution(Residual);
		     node[iPoint]->Set_OldSolution();
		     
		     if (config->GetKind_TimeIntScheme_Poisson()==EULER_IMPLICIT) {
				 Jacobian.DeleteValsRowi(iPoint);
		     }		     
		  break;
		  
		  case OPEN:
		     
		     MassFlux_Part = 0.0;
		     for (iDim = 0; iDim < nDim; iDim++) 
               MassFlux_Part -= solver_container[FLOW_SOL]->node[iPoint]->GetDensity()*(solver_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim))*Normal[iDim];
		     
		     /*if (MassFlux_Part >= 0.0) {
				 for (iVar = 0; iVar < nVar; iVar++)
				   LinSysRes.SetBlock_Zero(iPoint, iVar);
		     
		         node[iPoint]->SetSolution(Residual);
		         node[iPoint]->Set_OldSolution();
		     
		         if (config->GetKind_TimeIntScheme_Poisson()==EULER_IMPLICIT) {
					 Jacobian.DeleteValsRowi(iPoint);
		         }
				 
		     }*/
		    //else {
				 LinSysRes.SubtractBlock(iPoint, Residual);
				 
				 if (config->GetKind_TimeIntScheme_Poisson() == EULER_IMPLICIT) {
					 Jacobian_i[0][0] = 0.0;
					 Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
			    }
			 //}
		  break;
		}
     }
  }
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
}


void CPoissonSolverFVM::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                CConfig *config, unsigned short val_marker) {
									
su2double Poisson_Coeff_i,**Sol_i_Grad,Poissonval_i;
su2double Mom_Coeff_i[3],Proj_Mean_GradPoissonVar_Normal[3];
unsigned long iVertex, iPoint, jPoint, total_index;
unsigned short iDim, iVar;
su2double *Normal = new su2double[nDim];
string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
     /*--- Zero flux (Neumann) BC on pressure ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
      }

	 /*--- Add and subtract residual, and update Jacobians ---*/
		LinSysRes.SubtractBlock(iPoint, Residual);
		
       if (config->GetKind_TimeIntScheme_Poisson() == EULER_IMPLICIT) {
		   Jacobian_i[0][0] = 0.0;
		Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
	  }
     }
	}
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
}


void CPoissonSolverFVM::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
									 
su2double Poisson_Coeff_i,**Sol_i_Grad,Poissonval_i;
su2double Mom_Coeff_i[3],Proj_Mean_GradPoissonVar_Normal[3];
unsigned long iVertex, iPoint, jPoint, total_index;
unsigned short iDim, iVar;
su2double *Normal = new su2double[nDim];
string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
     /*--- Zero flux (Neumann) BC on pressure ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
      }

	 /*--- Add and subtract residual, and update Jacobians ---*/
		LinSysRes.SubtractBlock(iPoint, Residual);
		
       if (config->GetKind_TimeIntScheme_Poisson() == EULER_IMPLICIT) {
		   Jacobian_i[0][0] = 0.0;
		Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
	  }
     }
	}
	/*--- Free locally allocated memory ---*/
  delete [] Normal;
}
