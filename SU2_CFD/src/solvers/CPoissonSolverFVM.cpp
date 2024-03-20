/*!
 * \file solution_direct_poisson.cpp
 * \brief Main subrotuines for solving direct problems
 * \author F. Palacios
 * \version 8.0.0 "Harrier"
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

#include "../../include/solvers/CPoissonSolverFVM.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/gradients/computeGradientsGreenGauss.hpp"
#include "../../include/gradients/computeGradientsLeastSquares.hpp"
#include "../../include/limiters/computeLimiters.hpp"

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
  nVar = 1;
  nPrimVar = 1;
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;

Residual_RMS.resize(nVar,0.0); // TODO: Added PBFlow
  Residual_Max.resize(nVar,0.0);
  
  Residual = new su2double[nVar]; 
//   Residual_RMS = new su2double[nVar];
  Solution = new su2double[nVar];
//   Residual_Max = new su2double[nVar];
  

  /*--- Define some structures for locating max residuals ---*/
//     nVertex.resize(nMarker);
//   for (iMarker = 0; iMarker < nMarker; iMarker++)
//     nVertex[iMarker] = geometry->nVertex[iMarker];
  Point_Max.resize(nVar,0.0);
//   for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord.resize(nVar, nDim) = su2double(0.0);
//   for (iVar = 0; iVar < nVar; iVar++) {
//     Point_Max_Coord[iVar].resize(nDim,0.0);
//     // for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
//   }

  /*--- Define some auxiliar vector related with the solution ---*/

  Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
  
 if (config->GetKind_TimeIntScheme() == EULER_IMPLICIT) { // TODO:PBFlow
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

  /*--- Computation of gradients by least squares ---*/
  
  
  /*--- Always instantiate and initialize the variable to a zero value. ---*/

  nodes = new CPoissonVariable(0.0, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();
  
  /*--- The poisson equation always solved implicitly, so set the
   implicit flag in case we have periodic BCs. Geometry info is 
   communicated in the flow solver. ---*/

  SetImplicitPeriodic(true);

  /*--- Perform the MPI communication of the solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
}

CPoissonSolverFVM::~CPoissonSolverFVM(void) {
  
  if (nodes != nullptr) delete nodes;
}


void CPoissonSolverFVM::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {
}


void CPoissonSolverFVM::Preprocessing(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  unsigned long iPoint;
  /*--- Initialize the residual vector ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint ++)
    LinSysRes.SetBlock_Zero(iPoint);
  
  /*--- Initialize the Jacobian matrices ---*/
  Jacobian.SetValZero();

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
   SetSolution_Gradient_GG(geometry, config,false);

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) 
    SetSolution_Gradient_LS(geometry, config,false);  
}

void CPoissonSolverFVM::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                        unsigned short iMesh){
							
 /*--- Compute gradients so we can use it to find the velocity corrections ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) 
    SetSolution_Gradient_GG(geometry, config,false);

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) 
    SetSolution_Gradient_LS(geometry, config,false);
}

void CPoissonSolverFVM:: LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo){

}

void CPoissonSolverFVM::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  CNumerics* numerics = numerics_container[VISC_TERM];

  CMatrixView<const su2double> Sol_i_Grad,Sol_j_Grad;
  su2double Poisson_Coeff_i,Poisson_Coeff_j,Poissonval_i,Poissonval_j,Normal[MAXNDIM];
  su2double Mom_Coeff_i[MAXNDIM],Mom_Coeff_j[MAXNDIM];
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim;
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    /*--- Points coordinates, and normal vector ---*/
    numerics->SetCoord(geometry->nodes->GetCoord(iPoint),geometry->nodes->GetCoord(jPoint)); //TODO:PB_Flow
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));
    numerics->SetVolume(geometry->nodes->GetVolume(jPoint));

    /*--- Primitive variables w/o reconstruction ---*/
    Poissonval_i = nodes->GetSolution(iPoint, 0);
    Poissonval_j = nodes->GetSolution(jPoint, 0);

    numerics->SetPoissonval(Poissonval_i,Poissonval_j);

    Sol_i_Grad = nodes->GetGradient(iPoint);
    Sol_j_Grad = nodes->GetGradient(jPoint);

	numerics->SetConsVarGradient(Sol_i_Grad, Sol_j_Grad);

    if (config->GetKind_Incomp_System()!=ENUM_INCOMP_SYSTEM::PRESSURE_BASED) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Mom_Coeff_i[iDim] = 1.0;
        Mom_Coeff_j[iDim] = 1.0;
      }
	} else {	
      for (iDim = 0; iDim < nDim; iDim++) {
        Mom_Coeff_i[iDim] = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint)*solver_container[FLOW_SOL]->GetNodes()->Get_Mom_Coeff(iPoint, iDim) ;
        Mom_Coeff_j[iDim] = solver_container[FLOW_SOL]->GetNodes()->GetDensity(jPoint)*solver_container[FLOW_SOL]->GetNodes()->Get_Mom_Coeff(jPoint, iDim) ;
      }
    }
    numerics->SetInvMomCoeff(Mom_Coeff_i,Mom_Coeff_j);

    /*--- Compute and update residual ---*/
    auto residual = numerics->ComputeResidual(config);

    LinSysRes.SubtractBlock(iPoint, residual);
    LinSysRes.AddBlock(jPoint, residual);

    /*--- Implicit part ---*/
    if (implicit)
      Jacobian.UpdateBlocksSub(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
  }
}

void CPoissonSolverFVM::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                   CConfig *config, unsigned short iMesh) {
  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM];

  unsigned short iVar;
  unsigned long iPoint;
  su2double Src_Term;

  /*--- Initialize the source residual to zero ---*/
  for (iVar = 0; iVar < nVar; iVar++)
	  Residual[iVar] = 0.0;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Load the volume of the dual mesh cell ---*/
    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Compute the source term ---*/
    if ((config->GetKind_Incomp_System() == ENUM_INCOMP_SYSTEM::PRESSURE_BASED)) {
		Src_Term = solver_container[FLOW_SOL]->GetNodes()->GetMassFlux(iPoint) ;
		if (Src_Term != Src_Term)Src_Term = 0.0;
		nodes->SetSourceTerm(iPoint, Src_Term);
    }
    numerics->SetSourcePoisson(nodes->GetSourceTerm(iPoint));
    
    /*--- Compute the source residual ---*/
    numerics->ComputeResidual(Residual, Jacobian_i, config);

    /*--- Add the source residual to the total ---*/
    LinSysRes.AddBlock(iPoint, Residual);
    
    /*--- No jacobian contribution from mass flux/constant source term. ---*/
  }
}

void CPoissonSolverFVM::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  /*--- No time integration is done here. The routine is used as a means to solve the jacobian matrix in a way
   * consistent with the rest of the code. Time step is set to zero and no under-relaxation is applied to the
   * jacobian matrix.*/
  
  unsigned long iPoint, total_index, IterLinSol = 0;;
  unsigned short iVar;
  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  
  /*--- Set max linear solver iter which can be different from flow iterations ---*/
//   config->SetLinear_Solver_Iter(config->GetLinear_Solver_Iter_Poisson()); //TODO:PBFlow
//    SetLinear_Solver_Iter(config->GetLinear_Solver_Iter());

  /*--- Set maximum residual to zero ---*/ //TODO:Already done before PBFlow
//   for (iVar = 0; iVar < nVar; iVar++) {
    // Residual_RMS[iVar] = 0.0;
//     SetResidual_RMS(geometry, config);
//   }

    SetResToZero();

  /*--- Initialize residual and solution at the ghost points ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
/*--- Workaround to deal with nodes that are part of multiple boundaries and where
 *    one face might be strong BC and another weak BC (mostly for farfield boundary 
 *    where the boundary face is strong or weak depending on local flux. ---*/
    if (nodes->GetStrongBC(iPoint)) {
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
      LinSysRes.SetBlock_Zero(iPoint);
    }

    /*--- Read the residual ---*/
    local_Res_TruncError = nodes->GetResTruncError(iPoint);

	/*--- Read the volume ---*/
    Vol = geometry->nodes->GetVolume(iPoint);

    /*--- Possible under-relaxation if needed goes here. ---*/
    /*--- Currently, nothing changes. ---*/
    /*su2double *diag = Jacobian.GetBlock(iPoint, iPoint);
    for (iVar = 0; iVar < nVar; iVar++)
      diag[(nVar+1)*iVar] = diag[(nVar+1)*iVar]/1.0; 
    
    Jacobian.SetBlock(iPoint, iPoint, diag);*/

	/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar] );
      LinSysSol[total_index] = 0.0;
      Residual_RMS[iVar] += LinSysRes[total_index]*LinSysRes[total_index];
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint));
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

  /*--- Store the value of the residual. ---*/
  SetResLinSolver(System.GetResidual());

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      nodes->AddSolution(iPoint, iVar, LinSysSol[iPoint*nVar+iVar]);
     }
  }

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  /*--- MPI solution ---*/
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
}
 

void CPoissonSolverFVM::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {



}

void CPoissonSolverFVM::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                        unsigned short iMesh, unsigned long Iteration){

  /*--- Set a value for periodic communication routines. Is not used during any computations. 
   *    Done in the constructor---*/   
}

void CPoissonSolverFVM::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iVertex, iPoint;
  unsigned short iDim, iVar;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {
      /*--- The farfield boundary is considered as an inlet-outlet boundary, where flow 
       * can either enter or leave. For pressure, it is treated as a fully developed flow
       * and a dirichlet BC is applied. For velocity, based on the sign of massflux, either 
       * a dirichlet or a neumann BC is applied (in Flow_Correction routine). ---*/		

      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = 0.0;

      LinSysRes.SetBlock_Zero(iPoint);

      nodes->SetSolution(iPoint, Residual);
      nodes->SetSolution_Old(iPoint,Residual);
      nodes->SetStrongBC(iPoint);
      if (config->GetKind_TimeIntScheme() == EULER_IMPLICIT) { // TODO: PBFlow
        Jacobian.DeleteValsRowi(iPoint);
      }
    }
  }
}
                                
                                
                                
void CPoissonSolverFVM::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
unsigned long iVertex, iPoint, Point_Normal;
unsigned short iVar;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->nodes->GetDomain(iPoint)) {
      /*--- Zero flux (Neumann) BC on pressure ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = 0.0;

	  /*--- Add and subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      if (config->GetKind_TimeIntScheme() == EULER_IMPLICIT) { // TODO: PBFlow
        Jacobian_i[0][0] = 0.0;
        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
	  }
    }
  }
}
                                 
                                 
void CPoissonSolverFVM::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
								
unsigned long iVertex, iPoint;
unsigned short iVar;

  /*--- Only fixed velocity inlet is considered ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {
    /*--- Zero flux (Neumann) BC on pressure ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = 0.0;

      /*--- Add and subtract residual, and update Jacobians ---*/
		LinSysRes.SubtractBlock(iPoint, Residual);

       if (config->GetKind_TimeIntScheme() == EULER_IMPLICIT) { // TODO: PBFlow
         Jacobian_i[0][0] = 0.0;
		 Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
	   }
     }
   }
}


void CPoissonSolverFVM::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iVertex, iPoint;
  unsigned short iDim, iVar;
  su2double MassFlux_Part, Normal[MAXNDIM];
  string Marker_Tag  = config->GetMarker_All_TagBound(val_marker);
  auto Kind_Outlet = config->GetKind_Inc_Outlet(Marker_Tag);

  /*--- Only fully developed case is considered ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {
      
     /*--- A Dirichlet BC is applied at the fully developed outlet, pressure is
      * assumed to be uniform and assigned the value of P_ref. ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = 0.0;

      switch (Kind_Outlet) {
        case INC_OUTLET_TYPE::PRESSURE_OUTLET:
          LinSysRes.SetBlock_Zero(iPoint);
          nodes->SetSolution(iPoint, Residual);
          nodes->SetSolution_Old(iPoint, Residual);
          nodes->SetStrongBC(iPoint);
          if (config->GetKind_TimeIntScheme() == EULER_IMPLICIT)  // TODO: PBFlow
            Jacobian.DeleteValsRowi(iPoint);
        break;
		  
        case INC_OUTLET_TYPE::OPEN:
          /*--- Normal vector for this vertex (negative for outward convention) ---*/
          geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
          MassFlux_Part = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) 
            MassFlux_Part -= solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint)*(solver_container[FLOW_SOL]->GetNodes()->GetVelocity(iPoint, iDim))*Normal[iDim];

          LinSysRes.SubtractBlock(iPoint, Residual);
          if (config->GetKind_TimeIntScheme() == EULER_IMPLICIT) { // TODO:PBFlow
            Jacobian_i[0][0] = 0.0;
            Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
          }
        break;
        }
     }
  }

}


void CPoissonSolverFVM::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                CConfig *config, unsigned short val_marker) {
									
  unsigned long iVertex, iPoint;
  unsigned short iVar;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->nodes->GetDomain(iPoint)) {
      
     /*--- Zero flux (Neumann) BC on pressure ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
      }

	 /*--- Add and subtract residual, and update Jacobians ---*/
		LinSysRes.SubtractBlock(iPoint, Residual);
		
       if (config->GetKind_TimeIntScheme() == EULER_IMPLICIT) { // TODO: PBFlow
		   Jacobian_i[0][0] = 0.0;
		Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
	  }
     }
	}
}


void CPoissonSolverFVM::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
									 
unsigned long iVertex, iPoint, jPoint, total_index;
unsigned short iVar;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->nodes->GetDomain(iPoint)) {
      
     /*--- Zero flux (Neumann) BC on pressure ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
      }

	 /*--- Add and subtract residual, and update Jacobians ---*/
		LinSysRes.SubtractBlock(iPoint, Residual);
		
       if (config->GetKind_TimeIntScheme() == EULER_IMPLICIT) { // TODO: PBFlow
		   Jacobian_i[0][0] = 0.0;
		Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
	  }
     }
	}
}


/*--- Note that velocity indices in residual are hard coded in solver_structure. Need to be careful. ---*/

void CPoissonSolverFVM::BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                               CNumerics *numerics, CConfig *config) {
  
  /*--- Complete residuals for periodic boundary conditions. We loop over
   the periodic BCs in matching pairs so that, in the event that there are
   adjacent periodic markers, the repeated points will have their residuals
   accumulated correctly during the communications. For implicit calculations,
   the Jacobians and linear system are also correctly adjusted here. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
  }
  
}
