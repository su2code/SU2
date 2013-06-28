/*!
 * \file integration_time.cpp
 * \brief Time deppending numerical method.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.
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

#include "../include/integration_structure.hpp"

CMultiGridIntegration::CMultiGridIntegration(CConfig *config) : CIntegration(config) { }

CMultiGridIntegration::~CMultiGridIntegration(void) { }

void CMultiGridIntegration::SetMultiGrid_Solver(CGeometry ***geometry, CSolution ****solution_container, CNumerics *****solver_container, CConfig **config,
																								unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone) {
	unsigned short FinestMesh, iMGLevel;
	double monitor = 0.0;
	bool restart = (config[iZone]->GetRestart() || config[iZone]->GetRestart_Flow());
	bool startup_multigrid = ( config[iZone]->GetRestart_Flow() && (RunTime_EqSystem == RUNTIME_FLOW_SYS) && (Iteration == 0));
	bool direct = ((config[iZone]->GetKind_Solver() == EULER) || (config[iZone]->GetKind_Solver() == NAVIER_STOKES) || 
								 (config[iZone]->GetKind_Solver() == RANS) || (config[iZone]->GetKind_Solver() == FREE_SURFACE_EULER) ||
								 (config[iZone]->GetKind_Solver() == FREE_SURFACE_NAVIER_STOKES) || (config[iZone]->GetKind_Solver() == FREE_SURFACE_RANS) ||
								 (config[iZone]->GetKind_Solver() == FLUID_STRUCTURE_EULER) || (config[iZone]->GetKind_Solver() == FLUID_STRUCTURE_NAVIER_STOKES) || 
								 (config[iZone]->GetKind_Solver() == FLUID_STRUCTURE_RANS) || (config[iZone]->GetKind_Solver() == AEROACOUSTIC_EULER) ||
								 (config[iZone]->GetKind_Solver() == AEROACOUSTIC_NAVIER_STOKES) || (config[iZone]->GetKind_Solver() == AEROACOUSTIC_RANS) ||
								 (config[iZone]->GetKind_Solver() == PLASMA_EULER) || (config[iZone]->GetKind_Solver() == PLASMA_NAVIER_STOKES));
	
	/*--- If restart, update multigrid levels at the first multigrid iteration ---*/
	if ((restart && (Iteration == config[iZone]->GetnStartUpIter())) || startup_multigrid) {
		for (iMGLevel = 0; iMGLevel < config[iZone]->GetMGLevels(); iMGLevel++) {
			SetRestricted_Solution(RunTime_EqSystem, solution_container[iZone][iMGLevel], solution_container[iZone][iMGLevel+1], 
														 geometry[iZone][iMGLevel], geometry[iZone][iMGLevel+1], config[iZone], iMGLevel, true);
		}
	}
	
	/*--- Full multigrid strategy and start up with fine grid only works with the direct problem ---*/
	if (!config[iZone]->GetRestart() && config[iZone]->GetFullMG() &&  direct && ( GetConvergence_FullMG() && (config[iZone]->GetFinestMesh() != MESH_0 ))) {
		SetProlongated_Solution(RunTime_EqSystem, solution_container[iZone][config[iZone]->GetFinestMesh()-1], 
														solution_container[iZone][config[iZone]->GetFinestMesh()], geometry[iZone][config[iZone]->GetFinestMesh()-1],
														geometry[iZone][config[iZone]->GetFinestMesh()], config[iZone]);
		config[iZone]->SubtractFinestMesh();
	}
	
	/*--- Set the current finest grid ---*/
	FinestMesh = config[iZone]->GetFinestMesh();	

	/*--- Perform the Full Approximation Scheme (FAS) multigrid ---*/
	FAS_Multigrid(geometry, solution_container, solver_container, config, FinestMesh,
                config[iZone]->GetMGCycle(), RunTime_EqSystem, Iteration, iZone);

	/*--- Compute non-dimensional parameters and the convergence monitor ---*/
	NonDimensional_Parameters(geometry[iZone], solution_container[iZone], solver_container[iZone], config[iZone],
                                     FinestMesh, RunTime_EqSystem, Iteration, &monitor);

	/*--- Convergence strategy ---*/
	Convergence_Monitoring(geometry[iZone][FinestMesh], config[iZone], Iteration, monitor);

}

void CMultiGridIntegration::FAS_Multigrid(CGeometry ***geometry, CSolution ****solution_container, CNumerics *****solver_container,
                                          CConfig **config, unsigned short iMesh, unsigned short mu, unsigned short RunTime_EqSystem,
                                          unsigned long Iteration, unsigned short iZone) {
  
	unsigned short iPreSmooth, iPostSmooth, iRKStep, iRKLimit = 1;

	bool startup_multigrid = (config[iZone]->GetRestart_Flow() && (RunTime_EqSystem == RUNTIME_FLOW_SYS) && (Iteration == 0));
	unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

	/*--- Do a presmoothing on the grid iMesh to be restricted to the grid iMesh+1 ---*/
	for (iPreSmooth = 0; iPreSmooth < config[iZone]->GetMG_PreSmooth(iMesh); iPreSmooth++) {

		switch (config[iZone]->GetKind_TimeIntScheme()) {
		case RUNGE_KUTTA_EXPLICIT: iRKLimit = config[iZone]->GetnRKStep(); break;
		case EULER_EXPLICIT: case EULER_IMPLICIT: iRKLimit = 1; break; }

		/*--- Time and space integration ---*/
		for (iRKStep = 0; iRKStep < iRKLimit; iRKStep++) {

			/*--- Send-Receive boundary conditions ---*/
			solution_container[iZone][iMesh][SolContainer_Position]->MPI_Send_Receive(geometry, solution_container, config, iMesh, iZone);
			
      /*--- Preprocessing ---*/
			solution_container[iZone][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iMesh], solution_container[iZone][iMesh], solver_container[iZone][iMesh][SolContainer_Position], config[iZone], iRKStep);

			if (iRKStep == 0) {

				/*--- Set the old solution for the Runge-Kutta iteration ---*/
				solution_container[iZone][iMesh][SolContainer_Position]->Set_OldSolution(geometry[iZone][iMesh]);

				/*--- Compute time step, max eigenvalue, and integration scheme (steady and unsteady problems) ---*/
				solution_container[iZone][iMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iMesh], solution_container[iZone][iMesh], config[iZone], iMesh, Iteration);

				/*--- Restrict the solution and gradient for the adjoint problem ---*/
				Adjoint_Setup(geometry, solution_container, config, RunTime_EqSystem, Iteration, iZone);

			}

			/*--- Space integration ---*/
			Space_Integration(geometry[iZone][iMesh], solution_container[iZone][iMesh], solver_container[iZone][iMesh][SolContainer_Position],
					config[iZone], iMesh, iRKStep, RunTime_EqSystem);

			/*--- Time integration ---*/
			Time_Integration(geometry[iZone][iMesh], solution_container[iZone][iMesh], config[iZone], iRKStep,
					RunTime_EqSystem, Iteration);

			/*--- Postprocessing ---*/
			solution_container[iZone][iMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iMesh], solution_container[iZone][iMesh], config[iZone], iMesh);

		}
	}

	/*--- Compute Forcing Term $P_(k+1) = I^(k+1)_k(P_k+F_k(u_k))-F_(k+1)(I^(k+1)_k u_k)$ and update solution for multigrid ---*/
	if ((iMesh < config[iZone]->GetMGLevels() && ((Iteration >= config[iZone]->GetnStartUpIter()) || startup_multigrid)) ) { 

		/*--- 1st step, compute $r_k = P_k+F_k(u_k)$ ---*/
		solution_container[iZone][iMesh][SolContainer_Position]->MPI_Send_Receive(geometry, solution_container, config, iMesh, iZone);
		solution_container[iZone][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iMesh], solution_container[iZone][iMesh], solver_container[iZone][iMesh][SolContainer_Position], config[iZone], NO_RK_ITER);
		Space_Integration(geometry[iZone][iMesh], solution_container[iZone][iMesh], solver_container[iZone][iMesh][SolContainer_Position], config[iZone], iMesh, NO_RK_ITER, RunTime_EqSystem);
		SetResidual_Term(geometry[iZone][iMesh], solution_container[iZone][iMesh][SolContainer_Position]);

		/*--- 2nd step, compute $r_(k+1) = F_(k+1)(I^(k+1)_k u_k)$ ---*/
		SetRestricted_Solution(RunTime_EqSystem, solution_container[iZone][iMesh], solution_container[iZone][iMesh+1], geometry[iZone][iMesh], geometry[iZone][iMesh+1], config[iZone], iMesh, false);
		solution_container[iZone][iMesh+1][SolContainer_Position]->MPI_Send_Receive(geometry, solution_container, config, iMesh+1, iZone);
		solution_container[iZone][iMesh+1][SolContainer_Position]->Preprocessing(geometry[iZone][iMesh+1], solution_container[iZone][iMesh+1], solver_container[iZone][iMesh][SolContainer_Position], config[iZone], NO_RK_ITER);
		Space_Integration(geometry[iZone][iMesh+1], solution_container[iZone][iMesh+1], solver_container[iZone][iMesh+1][SolContainer_Position], config[iZone], iMesh+1, NO_RK_ITER, RunTime_EqSystem);

		/*--- 3rd step, compute $P_(k+1) = I^(k+1)_k(r_k) - r_(k+1) ---*/
		SetForcing_Term(solution_container[iZone][iMesh][SolContainer_Position],solution_container[iZone][iMesh+1][SolContainer_Position],geometry[iZone][iMesh],geometry[iZone][iMesh+1], config[iZone]);

		/*--- Recursive call to FAS_Multigrid ---*/
		for (unsigned short imu = 0; imu <= mu; imu++) {
			if (iMesh == config[iZone]->GetMGLevels()-2) 
				FAS_Multigrid(geometry, solution_container, solver_container, config, iMesh+1, 0, RunTime_EqSystem, Iteration, iZone);
			else 
				FAS_Multigrid(geometry, solution_container, solver_container, config, iMesh+1, mu, RunTime_EqSystem, Iteration, iZone);
		}

		/*--- Compute prolongated solution, and smooth the correction $u^(new)_k = u_k +  Smooth(I^k_(k+1)(u_(k+1)-I^(k+1)_k u_k))$ ---*/
		GetProlongated_Correction(solution_container[iZone][iMesh][SolContainer_Position], solution_container[iZone][iMesh+1][SolContainer_Position],
				geometry[iZone][iMesh],geometry[iZone][iMesh+1], config[iZone]);
		SmoothProlongated_Correction(RunTime_EqSystem, solution_container[iZone][iMesh], geometry[iZone][iMesh],
				config[iZone]->GetMG_CorrecSmooth(iMesh), 1.25, config[iZone]);
		SetProlongated_Correction(solution_container[iZone][iMesh][SolContainer_Position], geometry[iZone][iMesh], config[iZone]);

		/*--- Solution postsmoothing in the prolongated grid ---*/
		for (iPostSmooth = 0; iPostSmooth < config[iZone]->GetMG_PostSmooth(iMesh); iPostSmooth++) {

			switch (config[iZone]->GetKind_TimeIntScheme()) {
			case RUNGE_KUTTA_EXPLICIT: iRKLimit = config[iZone]->GetnRKStep(); break;
			case EULER_EXPLICIT: case EULER_IMPLICIT: iRKLimit = 1; break;
			}

			for (iRKStep = 0; iRKStep < iRKLimit; iRKStep++) {

				solution_container[iZone][iMesh][SolContainer_Position]->MPI_Send_Receive(geometry, solution_container, config, iMesh, iZone);
				solution_container[iZone][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iMesh], solution_container[iZone][iMesh], solver_container[iZone][iMesh][SolContainer_Position], config[iZone], iRKStep);

				if (iRKStep == 0) {
					solution_container[iZone][iMesh][SolContainer_Position]->Set_OldSolution(geometry[iZone][iMesh]);
					solution_container[iZone][iMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iMesh], solution_container[iZone][iMesh], config[iZone], iMesh, Iteration);
				}

				Space_Integration(geometry[iZone][iMesh], solution_container[iZone][iMesh], solver_container[iZone][iMesh][SolContainer_Position], config[iZone], iMesh, iRKStep, RunTime_EqSystem);
				Time_Integration(geometry[iZone][iMesh], solution_container[iZone][iMesh], config[iZone], iRKStep, RunTime_EqSystem, Iteration);

				solution_container[iZone][iMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iMesh], solution_container[iZone][iMesh], config[iZone], iMesh);

			}
		}
	}

}

void CMultiGridIntegration::GetProlongated_Correction(CSolution *sol_fine, CSolution *sol_coarse, CGeometry *geo_fine, 
		CGeometry *geo_coarse, CConfig *config) {

	unsigned long Point_Fine, Point_Coarse, iVertex;
	unsigned short Boundary, iMarker, iChildren, iVar;
	double *Solution_Fine, *Solution_Coarse, *Solution, Area_Parent, Area_Children;
	unsigned short nVar = sol_coarse->GetnVar();
	Solution = new double[nVar];

	for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
		Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
		for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;

		for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
			Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
			Area_Children = geo_fine->node[Point_Fine]->GetVolume();
			Solution_Fine = sol_fine->node[Point_Fine]->GetSolution();
			for (iVar = 0; iVar < nVar; iVar++)
				Solution[iVar] -= Solution_Fine[iVar]*Area_Children/Area_Parent;
		}

		Solution_Coarse = sol_coarse->node[Point_Coarse]->GetSolution();

		for (iVar = 0; iVar < nVar; iVar++)
			Solution[iVar] += Solution_Coarse[iVar];

		for (iVar = 0; iVar < nVar; iVar++)
			sol_coarse->node[Point_Coarse]->SetSolution_Old(Solution);
	}

	/*--- If we are dealing with a dirichlet wall, there is no contribution from the boundaries ---*/
	for (iMarker=0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		if (Boundary == NO_SLIP_WALL)
			for(iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
				Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
				sol_coarse->node[Point_Coarse]->SetVelSolutionOldZero();
			}
	}

	for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++)
		for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
			Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
			sol_fine->node[Point_Fine]->SetResidual(sol_coarse->node[Point_Coarse]->GetSolution_Old());				
		}

	delete [] Solution;
}

void CMultiGridIntegration::SmoothProlongated_Correction (unsigned short RunTime_EqSystem, CSolution **solution, CGeometry *geometry, 
		unsigned short val_nSmooth, double val_smooth_coeff, CConfig *config) {
	double *Residual_Old, *Residual_Sum, *Residual, *Residual_i, *Residual_j;
	unsigned short iVar, iSmooth, iMarker, nneigh;
	unsigned long iEdge, iPoint, jPoint, iVertex;
	unsigned short SolContainer_Position = config->GetContainerPosition(RunTime_EqSystem);
	unsigned short nVar = solution[SolContainer_Position]->GetnVar();

	if (val_nSmooth > 0) {

		Residual = new double [nVar];

		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			Residual_Old = solution[SolContainer_Position]->node[iPoint]->GetResidual();
			solution[SolContainer_Position]->node[iPoint]->SetResidual_Old(Residual_Old);
		}

		/*--- Jacobi iterations ---*/
		for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
				solution[SolContainer_Position]->node[iPoint]->SetResidualSumZero();

			/*--- Loop over Interior edges ---*/
			for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
				iPoint = geometry->edge[iEdge]->GetNode(0);
				jPoint = geometry->edge[iEdge]->GetNode(1);

				Residual_i = solution[SolContainer_Position]->node[iPoint]->GetResidual();			
				Residual_j = solution[SolContainer_Position]->node[jPoint]->GetResidual();

				/*--- Accumulate nearest neighbor Residual to Res_sum for each variable ---*/
				solution[SolContainer_Position]->node[iPoint]->AddResidual_Sum(Residual_j);
				solution[SolContainer_Position]->node[jPoint]->AddResidual_Sum(Residual_i);
			}

			/*--- Loop over all mesh points (Update Residuals with averaged sum) ---*/
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
				nneigh = geometry->node[iPoint]->GetnPoint();
				Residual_Sum = solution[SolContainer_Position]->node[iPoint]->GetResidual_Sum();
				Residual_Old = solution[SolContainer_Position]->node[iPoint]->GetResidual_Old();
				for (iVar = 0; iVar < nVar; iVar++) {
					Residual[iVar] =(Residual_Old[iVar] + val_smooth_coeff*Residual_Sum[iVar])
									/(1.0 + val_smooth_coeff*double(nneigh));
				}
				solution[SolContainer_Position]->node[iPoint]->SetResidual(Residual);
			}

			/*--- Copy boundary values ---*/
			for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
				for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					Residual_Old = solution[SolContainer_Position]->node[iPoint]->GetResidual_Old();
					solution[SolContainer_Position]->node[iPoint]->SetResidual(Residual_Old);
				}
		}

		delete [] Residual;

	}

}

void CMultiGridIntegration::SetProlongated_Correction(CSolution *sol_fine, CGeometry *geo_fine, CConfig *config) {
	unsigned long Point_Fine;
	unsigned short iVar;
	double *Solution_Fine, *Residual_Fine, *Solution;
	unsigned short nVar = sol_fine->GetnVar();
	Solution = new double [nVar];

	for (Point_Fine = 0; Point_Fine < geo_fine->GetnPointDomain(); Point_Fine++) {
		Residual_Fine = sol_fine->node[Point_Fine]->GetResidual();
		Solution_Fine = sol_fine->node[Point_Fine]->GetSolution();
		for (iVar = 0; iVar < nVar; iVar++)
			Solution[iVar] = Solution_Fine[iVar]+config->GetDamp_Correc_Prolong()*Residual_Fine[iVar];
		sol_fine->node[Point_Fine]->SetSolution(Solution);
	}
	delete [] Solution;
}


void CMultiGridIntegration::SetProlongated_Solution(unsigned short RunTime_EqSystem, CSolution **sol_fine, CSolution **sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
	unsigned long Point_Fine, Point_Coarse;
	unsigned short iChildren;
	unsigned short SolContainer_Position = config->GetContainerPosition(RunTime_EqSystem);

	for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++)
		for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
			Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
			sol_fine[SolContainer_Position]->node[Point_Fine]->SetSolution(sol_coarse[SolContainer_Position]->node[Point_Coarse]->GetSolution());
		}
}

void CMultiGridIntegration::SetForcing_Term(CSolution *sol_fine, CSolution *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
	unsigned long Point_Fine, Point_Coarse, iVertex;
	unsigned short iMarker, iVar, iChildren, nVar;
	double *Residual_Fine, *Residual;

	nVar = sol_coarse->GetnVar();
	Residual = new double[nVar];

	for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
		sol_coarse->node[Point_Coarse]->SetRes_TruncErrorZero();

		for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
		for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
			Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
			Residual_Fine = sol_fine->node[Point_Fine]->GetResidual();
			for (iVar = 0; iVar < nVar; iVar++)	
				Residual[iVar] += config->GetDamp_Res_Restric()*Residual_Fine[iVar];
		}
		sol_coarse->node[Point_Coarse]->AddRes_TruncError(Residual);
	}

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == NO_SLIP_WALL)
			for(iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
				Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
				sol_coarse->node[Point_Coarse]->SetVelRes_TruncErrorZero();
			}	
	}

	for(Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
		sol_coarse->node[Point_Coarse]->SubtractRes_TruncError(sol_coarse->node[Point_Coarse]->GetResidual());
	}

	delete [] Residual;
}

void CMultiGridIntegration::SetResidual_Term(CGeometry *geometry, CSolution *solution) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		solution->node[iPoint]->AddResidual(solution->node[iPoint]->GetRes_TruncError());

}

void CMultiGridIntegration::SetRestricted_Residual(CSolution *sol_fine, CSolution *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
	unsigned long iVertex, Point_Fine, Point_Coarse;
	unsigned short iMarker, iVar, iChildren, nVar;
	double *Residual_Fine, *Residual;

	nVar = sol_coarse->GetnVar();
	Residual = new double[nVar];

	for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
		sol_coarse->node[Point_Coarse]->SetRes_TruncErrorZero();

		for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
		for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
			Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
			Residual_Fine = sol_fine->node[Point_Fine]->GetResidual();
			for (iVar = 0; iVar < nVar; iVar++)	
				Residual[iVar] += Residual_Fine[iVar];
		}
		sol_coarse->node[Point_Coarse]->AddRes_TruncError(Residual);
	}

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == NO_SLIP_WALL)
			for(iVertex = 0; iVertex<geo_coarse->nVertex[iMarker]; iVertex++) {
				Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
				sol_coarse->node[Point_Coarse]->SetVelRes_TruncErrorZero();
			}
	}

	delete [] Residual;
}

void CMultiGridIntegration::SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolution **sol_fine, CSolution **sol_coarse, CGeometry *geo_fine, 
		CGeometry *geo_coarse, CConfig *config, unsigned short iMesh, bool InclSharedDomain) {
	unsigned long iVertex, Point_Fine, Point_Coarse;
	unsigned short iMarker, iVar, iChildren;
	double Area_Parent, Area_Children, *Solution_Fine, *Solution;
	unsigned short SolContainer_Position = config->GetContainerPosition(RunTime_EqSystem);
	unsigned short nVar = sol_coarse[SolContainer_Position]->GetnVar();

	Solution = new double[nVar];

	/*--- Compute coarse solution from fine solution ---*/
	for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
		Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();

		for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
		for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {

			Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
			Area_Children = geo_fine->node[Point_Fine]->GetVolume();
			Solution_Fine = sol_fine[SolContainer_Position]->node[Point_Fine]->GetSolution();
			for (iVar = 0; iVar < nVar; iVar++)
				Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;  
		}

		sol_coarse[SolContainer_Position]->node[Point_Coarse]->SetSolution(Solution);

	}

	/*--- Update solution at the boundary ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == NO_SLIP_WALL)
			for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
				Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
				if (SolContainer_Position == FLOW_SOL) sol_coarse[SolContainer_Position]->node[Point_Coarse]->SetVelSolutionZero();
				if (SolContainer_Position == ADJFLOW_SOL) sol_coarse[SolContainer_Position]->node[Point_Coarse]->SetVelSolutionDVector();
				if (SolContainer_Position == TURB_SOL) sol_coarse[SolContainer_Position]->node[Point_Coarse]->SetSolutionZero();
			}

#ifndef NO_MPI
		unsigned long nBuffer;
		double *Conserv_Var;
		short SendRecv;

		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) && (InclSharedDomain)) {
			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			nBuffer = geo_coarse->nVertex[iMarker]*nVar;

			/*--- Send information  ---*/
			if (SendRecv > 0) {
				double *Buffer_Send_U = new double [geo_coarse->nVertex[iMarker]*nVar];
				int send_to = SendRecv-1;
				for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
					Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
					Conserv_Var = sol_coarse[SolContainer_Position]->node[Point_Coarse]->GetSolution();
					for (iVar = 0; iVar < nVar; iVar++)
						Buffer_Send_U[iVertex*nVar + iVar] = Conserv_Var[iVar];
				}
				MPI::COMM_WORLD.Bsend(Buffer_Send_U,nBuffer,MPI::DOUBLE,send_to, 0);

				delete [] Buffer_Send_U;
			}
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
				double *Buffer_Receive_U = new double [geo_coarse->nVertex[iMarker]*nVar];
				int receive_from = abs(SendRecv)-1;
				MPI::COMM_WORLD.Recv(Buffer_Receive_U,nBuffer,MPI::DOUBLE,receive_from, 0);
				for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
					Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();	
					for (iVar = 0; iVar < nVar; iVar++)
						sol_coarse[SolContainer_Position]->node[Point_Coarse]->SetSolution(iVar, Buffer_Receive_U[iVertex*nVar+iVar]);
				}
				delete[] Buffer_Receive_U;
			}
		}	
#endif
	}

	delete [] Solution;

}

void CMultiGridIntegration::SetRestricted_Gradient(unsigned short RunTime_EqSystem, CSolution **sol_fine, CSolution **sol_coarse, CGeometry *geo_fine, 
		CGeometry *geo_coarse, CConfig *config) {
	unsigned long Point_Fine, Point_Coarse;
	unsigned short iVar, iDim, iChildren, nDim, nVar;
	double Area_Parent, Area_Children, **Gradient_fine, **Gradient;
	unsigned short SolContainer_Position = config->GetContainerPosition(RunTime_EqSystem);

	nDim = geo_coarse->GetnDim();
	nVar = sol_coarse[SolContainer_Position]->GetnVar();

	Gradient = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Gradient[iVar] = new double [nDim];

	for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPoint(); Point_Coarse++) {
		Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();

		for (iVar = 0; iVar < nVar; iVar++)
			for (iDim = 0; iDim < nDim; iDim++)
				Gradient[iVar][iDim] = 0.0;

		for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
			Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
			Area_Children = geo_fine->node[Point_Fine]->GetVolume();
			Gradient_fine = sol_fine[SolContainer_Position]->node[Point_Fine]->GetGradient();

			for (iVar = 0; iVar < nVar; iVar++)
				for (iDim = 0; iDim < nDim; iDim++)
					Gradient[iVar][iDim] += Gradient_fine[iVar][iDim]*Area_Children/Area_Parent;
		}
		sol_coarse[SolContainer_Position]->node[Point_Coarse]->SetGradient(Gradient);
	}

	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Gradient[iVar];
	delete [] Gradient;
}

void CMultiGridIntegration::NonDimensional_Parameters(CGeometry **geometry, CSolution ***solution_container, CNumerics ****solver_container, 
																											CConfig *config, unsigned short FinestMesh, unsigned short RunTime_EqSystem, unsigned long Iteration, 
																											double *monitor) {
	
	switch (RunTime_EqSystem) {
			
		case RUNTIME_FLOW_SYS:
			
			/*--- Calculate the inviscid and viscous forces ---*/
			solution_container[FinestMesh][FLOW_SOL]->Inviscid_Forces(geometry[FinestMesh], config);
			if (config->GetKind_ViscNumScheme() != NONE) solution_container[FinestMesh][FLOW_SOL]->Viscous_Forces(geometry[FinestMesh], config);
			
			/*--- Evaluate convergence monitor ---*/
			if (config->GetConvCriteria() == CAUCHY) {
				if (config->GetCauchy_Func_Flow() == DRAG_COEFFICIENT) (*monitor) = solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag();
				if (config->GetCauchy_Func_Flow() == LIFT_COEFFICIENT) (*monitor) = solution_container[FinestMesh][FLOW_SOL]->GetTotal_CLift();
				if (config->GetCauchy_Func_Flow() == NEARFIELD_PRESSURE) (*monitor) = solution_container[FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();
			}

			if (config->GetConvCriteria() == RESIDUAL) {
#ifdef NO_MPI
				(*monitor) = log10(solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0));
#else
				(*monitor) = log10(sqrt(solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0)));
#endif
			}

			break;
			
		case RUNTIME_PLASMA_SYS:

			/*--- Calculate the inviscid and viscous forces ---*/
			solution_container[FinestMesh][PLASMA_SOL]->Inviscid_Forces(geometry[FinestMesh], config);
			if (config->GetKind_ViscNumScheme() != NONE) solution_container[FinestMesh][PLASMA_SOL]->Viscous_Forces(geometry[FinestMesh], config);
			
			/*--- Evaluate convergence monitor ---*/
			if (config->GetConvCriteria() == RESIDUAL) {
#ifdef NO_MPI
				(*monitor) = log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(0));
#else
				(*monitor) = log10(sqrt(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(0)));
#endif
			}
			
			break;
			
		case RUNTIME_ADJFLOW_SYS:
			
			/*--- Calculate the inviscid and viscous sensitivities ---*/
			solution_container[FinestMesh][ADJFLOW_SOL]->Inviscid_Sensitivity(geometry[FinestMesh], solution_container[FinestMesh], solver_container[FinestMesh][ADJFLOW_SOL][BOUND_TERM], config);
			if (config->GetKind_ViscNumScheme() != NONE) solution_container[FinestMesh][ADJFLOW_SOL]->Viscous_Sensitivity(geometry[FinestMesh], solution_container[FinestMesh], solver_container[FinestMesh][ADJFLOW_SOL][BOUND_TERM], config);
			
			/*--- Smooth the inviscid and viscous sensitivities ---*/
			if (config->GetKind_SensSmooth() != NONE) solution_container[FinestMesh][ADJFLOW_SOL]->Smooth_Sensitivity(geometry[FinestMesh], solution_container[FinestMesh], solver_container[FinestMesh][ADJFLOW_SOL][BOUND_TERM], config);

			/*--- Evaluate convergence monitor ---*/
			if (config->GetConvCriteria() == CAUCHY) {
				if (config->GetCauchy_Func_AdjFlow() == SENS_GEOMETRY) (*monitor) = solution_container[FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Geo();
				if (config->GetCauchy_Func_AdjFlow() == SENS_MACH) (*monitor) = solution_container[FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Mach();
			}
			if (config->GetConvCriteria() == RESIDUAL) {
#ifdef NO_MPI
				(*monitor) = log10(solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(0));
#else
				(*monitor) = log10(sqrt(solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(0)));
#endif
			}
			
			break;
			
		case RUNTIME_LINFLOW_SYS:
			
			/*--- Calculate the inviscid and viscous forces ---*/
			solution_container[FinestMesh][LINFLOW_SOL]->Inviscid_DeltaForces(geometry[FinestMesh], solution_container[FinestMesh], config);
			if (config->GetKind_ViscNumScheme() != NONE) solution_container[FinestMesh][LINFLOW_SOL]->Viscous_DeltaForces(geometry[FinestMesh], config);
			
			/*--- Evaluate convergence monitor ---*/
			if (config->GetConvCriteria() == CAUCHY) {
				if (config->GetCauchy_Func_LinFlow() == DELTA_DRAG_COEFFICIENT) (*monitor) = solution_container[FinestMesh][LINFLOW_SOL]->GetTotal_CDeltaDrag();
				if (config->GetCauchy_Func_LinFlow() == DELTA_LIFT_COEFFICIENT) (*monitor) = solution_container[FinestMesh][LINFLOW_SOL]->GetTotal_CDeltaLift();
			}
			if (config->GetConvCriteria() == RESIDUAL) {
#ifdef NO_MPI
				(*monitor) = log10(solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(0));
#else
				(*monitor) = log10(sqrt(solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(0)));
#endif
			}
			
			break;
	}
}

CSingleGridIntegration::CSingleGridIntegration(CConfig *config) : CIntegration(config) { }

CSingleGridIntegration::~CSingleGridIntegration(void) { }

void CSingleGridIntegration::SetSingleGrid_Solver(CGeometry ***geometry, CSolution ****solution_container,
		CNumerics *****solver_container, CConfig **config, unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone) {
	
  double monitor = 0.0;
	unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

	solution_container[iZone][MESH_0][SolContainer_Position]->Set_OldSolution(geometry[iZone][MESH_0]);

	/*--- Send-Receive boundary conditions ---*/
	solution_container[iZone][MESH_0][SolContainer_Position]->MPI_Send_Receive(geometry, solution_container, config, MESH_0, iZone);

	/*--- Preprocessing ---*/
	solution_container[iZone][MESH_0][SolContainer_Position]->Preprocessing(geometry[iZone][MESH_0], solution_container[iZone][MESH_0], solver_container[iZone][MESH_0][SolContainer_Position], config[iZone], 0);
	
	/*--- Time step evaluation ---*/
	solution_container[iZone][MESH_0][SolContainer_Position]->SetTime_Step(geometry[iZone][MESH_0], solution_container[iZone][MESH_0], config[iZone], MESH_0, 0);
	
	/*--- Space integration ---*/
	Space_Integration(geometry[iZone][MESH_0], solution_container[iZone][MESH_0], solver_container[iZone][MESH_0][SolContainer_Position], 
										config[iZone], MESH_0, NO_RK_ITER, RunTime_EqSystem);
	
	/*--- Time integration ---*/
	Time_Integration(geometry[iZone][MESH_0], solution_container[iZone][MESH_0], config[iZone], NO_RK_ITER, 
									 RunTime_EqSystem, Iteration);
	
	/*--- Postprocessing ---*/
	solution_container[iZone][MESH_0][SolContainer_Position]->Postprocessing(geometry[iZone][MESH_0], solution_container[iZone][MESH_0], config[iZone], MESH_0);

	/*--- Compute adimensional parameters and the convergence monitor ---*/
	switch (RunTime_EqSystem) {	
			
		case RUNTIME_WAVE_SYS:
			
			/*--- Calculate the wave strength ---*/
			solution_container[iZone][MESH_0][WAVE_SOL]->Wave_Strength(geometry[iZone][MESH_0], config[iZone]);
			
			/*--- Evaluate convergence monitor ---*/
			if (config[iZone]->GetConvCriteria() == CAUCHY) {
				monitor = solution_container[iZone][MESH_0][WAVE_SOL]->GetTotal_CWave();
			}
			
			if (config[iZone]->GetConvCriteria() == RESIDUAL)
#ifdef NO_MPI
				monitor = log10(solution_container[iZone][MESH_0][WAVE_SOL]->GetRes_Max(0));
#else
			monitor = log10(sqrt(solution_container[iZone][MESH_0][WAVE_SOL]->GetRes_Max(0)));
#endif
			break;
			
		case RUNTIME_FEA_SYS:
#ifdef NO_MPI
			monitor = log10(solution_container[iZone][MESH_0][FEA_SOL]->GetRes_Max(0));
#else
			monitor = log10(sqrt(solution_container[iZone][MESH_0][FEA_SOL]->GetRes_Max(0)));
#endif
			break;
	}
	
	/*--- Convergence strategy ---*/
	Convergence_Monitoring(geometry[iZone][MESH_0], config[iZone], Iteration, monitor);
	
}
