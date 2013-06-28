/*!
 * \file integration_time.cpp
 * \brief Time deppending numerical method.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
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

void CMultiGridIntegration::SetMultiGrid_Solver(CGeometry **geometry, CSolution ***solution_container, CNumerics ****solver_container, CConfig *config,
											 unsigned short RunTime_EqSystem, unsigned long Iteration) {
	unsigned short FinestMesh, iMGLevel;
	double monitor = 0.0;
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	bool startup_multigrid = ( config->GetRestart_Flow() && (RunTime_EqSystem == RUNTIME_FLOW_SYS) && (Iteration == 0));

 	/*--- If restart, update multigrid levels at the first multigrid iteration ---*/
	if ((restart && (Iteration == config->GetnStartUpIter())) || startup_multigrid)
		for (iMGLevel = 0; iMGLevel < config->GetMGLevels(); iMGLevel++) {
			SetRestricted_Solution(RunTime_EqSystem, solution_container[iMGLevel], solution_container[iMGLevel+1], 
								   geometry[iMGLevel], geometry[iMGLevel+1], config, iMGLevel, true);
		}
	
	/*--- Full multigrid strategy and start up with fine grid only works with the direct problem ---*/
	if ( !config->GetRestart() && config->GetFullMG() && 
		( (config->GetKind_Solver() == EULER) || (config->GetKind_Solver() == NAVIER_STOKES) || 
		 (config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == TWO_PHASE_FLOW) ) &&
		( GetConvergence_FullMG() && (config->GetFinestMesh() != MESH_0) ) ) {
		SetProlongated_Solution(RunTime_EqSystem, solution_container[config->GetFinestMesh()-1], 
								solution_container[config->GetFinestMesh()], geometry[config->GetFinestMesh()-1], 
								geometry[config->GetFinestMesh()], config);
		config->SubtractFinestMesh();
	}
	
	/*--- Set the current finest grid ---*/
	FinestMesh = config->GetFinestMesh();	

	/*--- Do the Full Approximation Scheme(FAS) multigrid ---*/
	FAS_Multigrid(geometry, solution_container, solver_container, config, FinestMesh, config->GetMGCycle(), RunTime_EqSystem, Iteration);

	/*--- Compute adimensional parameters and the convergence monitor ---*/
	switch (RunTime_EqSystem) {	
		case RUNTIME_FLOW_SYS:
			
			/*--- Calculates the inviscid forces here ---*/
			solution_container[FinestMesh][FLOW_SOL]->Inviscid_Forces(geometry[FinestMesh], config);
			
			/*--- Calculate the viscous forces here ---*/			
			if (config->GetKind_ViscNumScheme() != NONE) solution_container[FinestMesh][FLOW_SOL]->Viscous_Forces(geometry[FinestMesh], config);
			if (config->GetConvCriteria() == CAUCHY) {
				if (config->GetCauchy_Func_Flow() == DRAG_COEFFICIENT) monitor = solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag();
				if (config->GetCauchy_Func_Flow() == LIFT_COEFFICIENT) monitor = solution_container[FinestMesh][FLOW_SOL]->GetTotal_CLift();
				if (config->GetCauchy_Func_Flow() == NEARFIELD_PRESSURE) monitor = solution_container[FinestMesh][FLOW_SOL]->GetTotal_CNearFieldPress();
			}
			if (config->GetConvCriteria() == RESIDUAL) {
					monitor = log10(solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0));
			}
			break;
		case RUNTIME_PLASMA_SYS:
			
			/*--- Calculates the inviscid forces here ---*/
			solution_container[FinestMesh][PLASMA_SOL]->Inviscid_Forces(geometry[FinestMesh], config);
			
			/*--- Calculate the viscous forces here ---*/			
			if (config->GetKind_ViscNumScheme() != NONE) solution_container[FinestMesh][FLOW_SOL]->Viscous_Forces(geometry[FinestMesh], config);
			if (config->GetConvCriteria() == CAUCHY) {
				if (config->GetCauchy_Func_Flow() == DRAG_COEFFICIENT) monitor = solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CDrag();
				if (config->GetCauchy_Func_Flow() == LIFT_COEFFICIENT) monitor = solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CLift();
				if (config->GetCauchy_Func_Flow() == NEARFIELD_PRESSURE) monitor = solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CNearFieldPress();
			}
			if (config->GetConvCriteria() == RESIDUAL) {
					monitor = log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(0));

			}
			break;
		case RUNTIME_ADJFLOW_SYS: case RUNTIME_FLOW_ADJFLOW_SYS:
			if (config->GetKind_ViscNumScheme() == NONE) solution_container[FinestMesh][ADJFLOW_SOL]->Inviscid_Sensitivity(geometry[FinestMesh], solution_container[FinestMesh], config);
			else solution_container[FinestMesh][ADJFLOW_SOL]->Viscous_Sensitivity(geometry[FinestMesh], solution_container[FinestMesh], config);
			if (config->GetConvCriteria() == CAUCHY) {
				if (config->GetCauchy_Func_AdjFlow() == SENS_GEOMETRY) monitor = solution_container[FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_Geo();
				if (config->GetCauchy_Func_AdjFlow() == SENS_MACH) monitor = solution_container[FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_Mach();
			}
			if (config->GetConvCriteria() == RESIDUAL) {
				monitor = log10(solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(0));
			}
			break;
		case RUNTIME_LINFLOW_SYS:
			solution_container[FinestMesh][LINFLOW_SOL]->Inviscid_DeltaForces(geometry[FinestMesh], solution_container[FinestMesh], config);
			if (config->GetKind_ViscNumScheme() != NONE) solution_container[FinestMesh][LINFLOW_SOL]->Viscous_DeltaForces(geometry[FinestMesh], config);
			if (config->GetConvCriteria() == CAUCHY) {
				if (config->GetCauchy_Func_LinFlow() == DELTA_DRAG_COEFFICIENT) monitor = solution_container[FinestMesh][LINFLOW_SOL]->GetTotal_CDeltaDrag();
				if (config->GetCauchy_Func_LinFlow() == DELTA_LIFT_COEFFICIENT) monitor = solution_container[FinestMesh][LINFLOW_SOL]->GetTotal_CDeltaLift();
			}
			if (config->GetConvCriteria() == RESIDUAL) {
				monitor = log10(solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(0));
			}
			break;
	}
	
	/*--- Convergence strategy ---*/
	Convergence_Monitoring(geometry[FinestMesh], config, Iteration, monitor);
	
}

void CMultiGridIntegration::FAS_Multigrid(CGeometry **geometry, CSolution ***solution_container, 
										  CNumerics ****solver_container, CConfig *config, unsigned short iMesh, 
										  unsigned short mu, unsigned short RunTime_EqSystem, unsigned long Iteration) {
	unsigned short iPreSmooth, iPostSmooth, iRKStep, iRKLimit = 1, iMGLevel;
	bool startup_multigrid = (config->GetRestart_Flow() && (RunTime_EqSystem == RUNTIME_FLOW_SYS) && (Iteration == 0));
	unsigned short SolContainer_Position = config->GetContainerPosition(RunTime_EqSystem);
	/*--- Do a presmoothing on the grid iMesh to be restricted to the grid iMesh+1 ---*/
	for (iPreSmooth = 0; iPreSmooth < config->GetMG_PreSmooth(iMesh); iPreSmooth++) {
		
		/*--- Set the old solution for the Runge-Kutta iteration ---*/
		solution_container[iMesh][SolContainer_Position]->Set_OldSolution(geometry[iMesh]);

		/*--- Compute time step and integration scheme ---*/
		solution_container[iMesh][SolContainer_Position]->SetTime_Step(geometry[iMesh], solution_container[iMesh], config, iMesh);

		switch (config->GetKind_TimeIntScheme()) {
			case RUNGE_KUTTA_EXPLICIT: iRKLimit = config->GetnRKStep(); break;	
			case EULER_EXPLICIT: case EULER_IMPLICIT: iRKLimit = 1; break; }
		
		/*--- Restrict the solution and gradient for the adjoint problem ---*/
		if ( ( ((RunTime_EqSystem == RUNTIME_ADJFLOW_SYS) || (RunTime_EqSystem == RUNTIME_LINFLOW_SYS)) && (Iteration == config->GetnStartUpIter())) 
			|| ((RunTime_EqSystem == RUNTIME_FLOW_ADJFLOW_SYS) || (RunTime_EqSystem == RUNTIME_FLOW_LINFLOW_SYS)) )
			for (iMGLevel = 0; iMGLevel <= config->GetMGLevels(); iMGLevel++) {
				solution_container[iMGLevel][FLOW_SOL]->SetTime_Step(geometry[iMGLevel], solution_container[iMGLevel], config, iMGLevel);
				solution_container[iMGLevel][FLOW_SOL]->SetTotal_CDrag(solution_container[MESH_0][FLOW_SOL]->GetTotal_CDrag());
				solution_container[iMGLevel][FLOW_SOL]->SetTotal_CLift(solution_container[MESH_0][FLOW_SOL]->GetTotal_CLift());
				if (config->GetKind_ConvNumScheme() == SPACE_CENTRED) {
					solution_container[iMGLevel][FLOW_SOL]->SetSpectral_Radius(geometry[iMGLevel], config);
					solution_container[iMGLevel][FLOW_SOL]->SetPress_Switch(geometry[iMGLevel]);
				}
				if (iMGLevel != config->GetMGLevels()) {
					SetRestricted_Solution(RUNTIME_FLOW_SYS, solution_container[iMGLevel], solution_container[iMGLevel+1], 
										   geometry[iMGLevel], geometry[iMGLevel+1], config, iMGLevel, false);
					SetRestricted_Gradient(RUNTIME_FLOW_SYS, solution_container[iMGLevel], solution_container[iMGLevel+1], 
										   geometry[iMGLevel], geometry[iMGLevel+1], config);
				}
			}
		
		/*--- Time and space integration ---*/
		for (iRKStep = 0; iRKStep < iRKLimit; iRKStep++) {

			Space_Integration(geometry[iMesh], solution_container[iMesh], solver_container[iMesh][SolContainer_Position],
												config, iMesh, iRKStep, RunTime_EqSystem);

			Time_Integration(geometry[iMesh], solution_container[iMesh], config, iRKStep,
											 RunTime_EqSystem, Iteration);
		}
	}
	
	/*--- Compute Forcing Term $P_(k+1) = I^(k+1)_k(P_k+F_k(u_k))-F_(k+1)(I^(k+1)_k u_k)$ and update solution for multigrid ---*/
	if ((iMesh < config->GetMGLevels() && ((Iteration >= config->GetnStartUpIter()) || startup_multigrid)) ) { 
		
		/*--- 1st step, compute $r_k = P_k+F_k(u_k)$ ---*/
		Space_Integration(geometry[iMesh], solution_container[iMesh], solver_container[iMesh][SolContainer_Position], config, iMesh, NO_RK_ITER, RunTime_EqSystem);
		SetResidual_Term(geometry[iMesh], solution_container[iMesh][SolContainer_Position]);
		
		/*--- 2nd step, compute $r_(k+1) = F_(k+1)(I^(k+1)_k u_k)$ ---*/
	if (!config->GetReStartMGCycle()) solution_container[iMesh+1][SolContainer_Position]->Set_MultiSolution(geometry[iMesh+1],1);
		SetRestricted_Solution(RunTime_EqSystem, solution_container[iMesh], solution_container[iMesh+1], geometry[iMesh], geometry[iMesh+1], config, iMesh, false);
		Space_Integration(geometry[iMesh+1], solution_container[iMesh+1], solver_container[iMesh+1][SolContainer_Position], config, iMesh+1, NO_RK_ITER, RunTime_EqSystem);
	if (!config->GetReStartMGCycle()) solution_container[iMesh+1][SolContainer_Position]->Set_MultiSolution(geometry[iMesh+1],-1);
		
		/*--- 3rd step, compute $P_(k+1) = I^(k+1)_k(r_k) - r_(k+1) ---*/
		SetForcing_Term(solution_container[iMesh][SolContainer_Position],solution_container[iMesh+1][SolContainer_Position],geometry[iMesh],geometry[iMesh+1], config);
		
		/*--- Recursive call to FAS_Multigrid ---*/
		for (unsigned short imu = 0; imu <= mu; imu++) {
			if (iMesh == config->GetMGLevels()-2) 
				FAS_Multigrid(geometry, solution_container, solver_container, config, iMesh+1, 0, RunTime_EqSystem, Iteration);
			else 
				FAS_Multigrid(geometry, solution_container, solver_container, config, iMesh+1, mu, RunTime_EqSystem, Iteration);
		}
		
		/*--- Compute prolongated solution, and smooth the correction $u^(new)_k = u_k +  Smooth(I^k_(k+1)(u_(k+1)-I^(k+1)_k u_k))$ ---*/
		GetProlongated_Correction(solution_container[iMesh][SolContainer_Position], solution_container[iMesh+1][SolContainer_Position],
															geometry[iMesh],geometry[iMesh+1], config);
		SmoothProlongated_Correction(RunTime_EqSystem, solution_container[iMesh], geometry[iMesh],
																 config->GetMG_CorrecSmooth(iMesh), 1.25, config);
		SetProlongated_Correction(solution_container[iMesh][SolContainer_Position], geometry[iMesh], config);
		
		/*--- Solution postsmoothing in the prolongated grid ---*/
		for (iPostSmooth = 0; iPostSmooth < config->GetMG_PostSmooth(iMesh); iPostSmooth++) {
			solution_container[iMesh][SolContainer_Position]->Set_OldSolution(geometry[iMesh]);
			solution_container[iMesh][SolContainer_Position]->SetTime_Step(geometry[iMesh], solution_container[iMesh], config, iMesh);	
			switch (config->GetKind_TimeIntScheme()) {
				case RUNGE_KUTTA_EXPLICIT: iRKLimit = config->GetnRKStep(); break;	
				case EULER_EXPLICIT: case EULER_IMPLICIT: iRKLimit = 1; break; 
			}
			
			for (iRKStep = 0; iRKStep < config->GetnRKStep(); iRKStep++) {
				Space_Integration(geometry[iMesh], solution_container[iMesh], solver_container[iMesh][SolContainer_Position], config, iMesh, iRKStep, RunTime_EqSystem);
				Time_Integration(geometry[iMesh], solution_container[iMesh], config, iRKStep, RunTime_EqSystem, Iteration);
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
		sol_coarse->node[Point_Coarse]->SetTruncationErrorZero();
		
		for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
		for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
			Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
			Residual_Fine = sol_fine->node[Point_Fine]->GetResidual();
			for (iVar = 0; iVar < nVar; iVar++)	
				Residual[iVar] += config->GetDamp_Res_Restric()*Residual_Fine[iVar];
		}
		sol_coarse->node[Point_Coarse]->AddTruncationError(Residual);
	}
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == NO_SLIP_WALL)
			for(iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
				Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
				sol_coarse->node[Point_Coarse]->SetVelTruncationErrorZero();
			}
	}
	
	for(Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
		sol_coarse->node[Point_Coarse]->SubtractTruncationError(sol_coarse->node[Point_Coarse]->GetResidual());
	}
	
	delete [] Residual;
}

void CMultiGridIntegration::SetResidual_Term(CGeometry *geometry, CSolution *solution) {
	unsigned long iPoint;
	
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		solution->node[iPoint]->AddResidual(solution->node[iPoint]->GetTruncationError());
	
}

void CMultiGridIntegration::SetRestricted_Residual(CSolution *sol_fine, CSolution *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
	unsigned long iVertex, Point_Fine, Point_Coarse;
	unsigned short iMarker, iVar, iChildren, nVar;
	double *Residual_Fine, *Residual;
	
	nVar = sol_coarse->GetnVar();
	Residual = new double[nVar];
	
	for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
		sol_coarse->node[Point_Coarse]->SetTruncationErrorZero();
		
		for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
		for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
			Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
			Residual_Fine = sol_fine->node[Point_Fine]->GetResidual();
			for (iVar = 0; iVar < nVar; iVar++)	
				Residual[iVar] += Residual_Fine[iVar];
		}
		sol_coarse->node[Point_Coarse]->AddTruncationError(Residual);
	}
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == NO_SLIP_WALL)
			for(iVertex = 0; iVertex<geo_coarse->nVertex[iMarker]; iVertex++) {
				Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
				sol_coarse->node[Point_Coarse]->SetVelTruncationErrorZero();
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
				double Buffer_Send_U[geo_coarse->nVertex[iMarker]][nVar];
				int send_to = SendRecv;
				for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
					Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
					Conserv_Var = sol_coarse[SolContainer_Position]->node[Point_Coarse]->GetSolution();
					for (iVar = 0; iVar < nVar; iVar++)
						Buffer_Send_U[iVertex][iVar] = Conserv_Var[iVar];
				}
				MPI::COMM_WORLD.Bsend(&Buffer_Send_U,nBuffer,MPI::DOUBLE,send_to, 0);
			}
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
				double Buffer_Receive_U[geo_coarse->nVertex[iMarker]][nVar];
				int receive_from = abs(SendRecv);
				MPI::COMM_WORLD.Recv(&Buffer_Receive_U,nBuffer,MPI::DOUBLE,receive_from, 0);
				for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
					Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();	
					for (iVar = 0; iVar < nVar; iVar++)
						sol_coarse[SolContainer_Position]->node[Point_Coarse]->SetSolution(iVar, Buffer_Receive_U[iVertex][iVar]);
				}
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

void CMultiGridIntegration::SetSolution_Smoothing (unsigned short RunTime_EqSystem, CSolution **solution, CGeometry *geometry, 
												   unsigned short val_nSmooth, double val_smooth_coeff, CConfig *config) {
	double *Solution_Old, *Solution_Sum, *Solution, *Solution_i, *Solution_j;
	unsigned short iVar, iSmooth, iMarker, nneigh;
	unsigned long iEdge, iPoint, jPoint, iVertex;
	unsigned short SolContainer_Position = config->GetContainerPosition(RunTime_EqSystem);
	unsigned short nVar = solution[SolContainer_Position]->GetnVar();
	
	Solution = new double [nVar];
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		Solution_Old = solution[SolContainer_Position]->node[iPoint]->GetSolution();
		solution[SolContainer_Position]->node[iPoint]->SetSolution_Old(Solution_Old);
	}
	
	/*--- Jacobi iterations ---*/
	for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			solution[SolContainer_Position]->node[iPoint]->SetResidualSumZero();
		
		/*--- Loop over Interior edges ---*/
		for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
			iPoint = geometry->edge[iEdge]->GetNode(0);
			jPoint = geometry->edge[iEdge]->GetNode(1);
			
			Solution_i = solution[SolContainer_Position]->node[iPoint]->GetSolution();			
			Solution_j = solution[SolContainer_Position]->node[jPoint]->GetSolution();
			
			/*--- Accumulate nearest neighbor Solution to Res_sum for each variable ---*/
			solution[SolContainer_Position]->node[iPoint]->AddResidual_Sum(Solution_j);
			solution[SolContainer_Position]->node[jPoint]->AddResidual_Sum(Solution_i);
		}
		
		/*--- Loop over all mesh points (Update Solutions with averaged sum) ---*/
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			nneigh = geometry->node[iPoint]->GetnPoint();
			Solution_Sum = solution[SolContainer_Position]->node[iPoint]->GetResidual_Sum();
			Solution_Old = solution[SolContainer_Position]->node[iPoint]->GetResidual_Old();
			for (iVar = 0; iVar < nVar; iVar++) {
				Solution[iVar] =(Solution_Old[iVar] + val_smooth_coeff*Solution_Sum[iVar])
				/(1.0 + val_smooth_coeff*double(nneigh));
			}
			solution[SolContainer_Position]->node[iPoint]->SetSolution(Solution);
		}
		
		/*--- Copy boundary values ---*/
		for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
			for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Solution_Old = solution[SolContainer_Position]->node[iPoint]->GetSolution_Old();
				solution[SolContainer_Position]->node[iPoint]->SetSolution(Solution_Old);
			}
	}
	
	delete [] Solution;
}

CSingleGridIntegration::CSingleGridIntegration(CConfig *config) : CIntegration(config) { }

CSingleGridIntegration::~CSingleGridIntegration(void) { }

void CSingleGridIntegration::SetSingleGrid_Solver(CGeometry **geometry, CSolution ***solution_container, 
											 CNumerics ****solver_container, CConfig *config, 
											 unsigned short RunTime_EqSystem, unsigned long Iteration) {
	unsigned short SolContainer_Position = config->GetContainerPosition(RunTime_EqSystem);
	
	solution_container[MESH_0][SolContainer_Position]->Set_OldSolution(geometry[MESH_0]);
	
	Space_Integration(geometry[MESH_0], solution_container[MESH_0], solver_container[MESH_0][SolContainer_Position], 
					  config, MESH_0, NO_RK_ITER, RunTime_EqSystem);
	
	Time_Integration(geometry[MESH_0], solution_container[MESH_0], config, NO_RK_ITER, 
					 RunTime_EqSystem, Iteration);
	
}

