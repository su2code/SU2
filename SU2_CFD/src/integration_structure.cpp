/*!
 * \file integration_structure.cpp
 * \brief This subroutine includes the space and time integration structure.
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

CIntegration::CIntegration(CConfig *config) {
	Cauchy_Value = 0;
	Cauchy_Func = 0;
	Old_Func = 0;
	New_Func = 0;
	Cauchy_Counter = 0;
	Convergence = false;
	Convergence_OneShot = false;
	Convergence_FullMG = false;
	Cauchy_Serie = new double [config->GetCauchy_Elems()+1];	
}

CIntegration::~CIntegration(void) {
	delete [] Cauchy_Serie;
}

void CIntegration::Space_Integration(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, 
		CConfig *config, unsigned short iMesh, unsigned short iRKStep,
		unsigned short RunTime_EqSystem) {
	unsigned short iMarker;
	unsigned short MainSolution = config->GetContainerPosition(RunTime_EqSystem);

	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		switch (config->GetMarker_All_Boundary(iMarker)) {
		case SEND_RECEIVE:
			solution_container[MainSolution]->BC_Send_Receive(geometry, NULL, config, iMarker, iMesh);
			break;
		}

	/*--- Do some preprocessing stuff ---*/
	solution_container[MainSolution]->Preprocessing(geometry, solution_container, config, iRKStep);

	/*--- Compute inviscid residuals ---*/
	switch (config->GetKind_ConvNumScheme()) {
	case NONE :
		break;
	case SPACE_CENTRED:
		solution_container[MainSolution]->Centred_Residual(geometry, solution_container, solver[CONV_TERM], config, iMesh, iRKStep);
		break;
	case SPACE_UPWIND:
		solution_container[MainSolution]->Upwind_Residual(geometry, solution_container, solver[CONV_TERM], config, iMesh);
		break;
	}

	/*--- Compute viscous residuals ---*/
	switch (config->GetKind_ViscNumScheme()) {
	case NONE :
		break;
	case AVG_GRAD: case AVG_GRAD_CORRECTED:
		solution_container[MainSolution]->Viscous_Residual(geometry, solution_container, solver[VISC_TERM], config, iMesh, iRKStep);
		break;
	case GALERKIN:
		solution_container[MainSolution]->Galerkin_Method(geometry, solution_container, solver[VISC_TERM], config, iMesh);
		break;
	}

	/*--- Compute source term residuals ---*/
	switch (config->GetKind_SourNumScheme()) {
	case NONE :
		break;
	case PIECEWISE_CONSTANT:
		solution_container[MainSolution]->SourcePieceWise_Residual(geometry, solution_container, solver[SOUR_TERM], config, iMesh);
		solution_container[MainSolution]->SourceConserv_Residual(geometry, solution_container, solver[CONS_SOUR_TERM], config, iMesh);
		break;
	}	

	/*--- Weak boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		switch (config->GetMarker_All_Boundary(iMarker)) {
		case EULER_WALL:
			solution_container[MainSolution]->BC_Euler_Wall(geometry, solution_container, config, iMarker);
			break;
		case INLET_FLOW:
			solution_container[MainSolution]->BC_Inlet(geometry, solution_container, solver[BOUND_TERM], config, iMarker);
			break;
		case OUTLET_FLOW:
			solution_container[MainSolution]->BC_Outlet(geometry, solution_container, solver[BOUND_TERM], config, iMarker);
			break;
		case FAR_FIELD:
			solution_container[MainSolution]->BC_Far_Field(geometry, solution_container, solver[BOUND_TERM], config, iMarker);
			break;
		case SYMMETRY_PLANE:
			solution_container[MainSolution]->BC_Sym_Plane(geometry, solution_container, solver[BOUND_TERM], config, iMarker);
			break;
		case INTERFACE_BOUNDARY:
			solution_container[MainSolution]->BC_Interface_Boundary(geometry, solution_container, solver[BOUND_TERM], config, iMarker);
			break;
		case NEARFIELD_BOUNDARY:
			solution_container[MainSolution]->BC_NearField_Boundary(geometry, solution_container, solver[BOUND_TERM], config, iMarker);
			break;
		case ELECTRODE_BOUNDARY:
			solution_container[MainSolution]->BC_Electrode(geometry, solution_container, solver[BOUND_TERM], config, iMarker);
			break;
		case DIELECTRIC_BOUNDARY:
			solution_container[MainSolution]->BC_Dielectric(geometry, solution_container, solver[BOUND_TERM], config, iMarker);
			break;
		}

	/*--- Strong boundary conditions (Navier-Stokes, and interprocessor) ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		switch (config->GetMarker_All_Boundary(iMarker)) {
		case NO_SLIP_WALL:
			solution_container[MainSolution]->BC_NS_Wall(geometry, solution_container, config, iMarker);
			break;
		case DIRICHLET:
			solution_container[MainSolution]->BC_Dirichlet(geometry, solution_container, config, iMarker);
			break;
		case CUSTOM_BOUNDARY:
			solution_container[MainSolution]->BC_Custom(geometry, solution_container, solver[BOUND_TERM], config, iMarker);
			break;
		case SEND_RECEIVE:
			solution_container[MainSolution]->BC_InterProcessor(geometry, NULL, config, iMarker, iMesh);
			break;
		}

	/*--- Add viscous and convective residuals, and compute the Dual Time Source term ---*/
	if ((RunTime_EqSystem == RUNTIME_FLOW_SYS) || (RunTime_EqSystem == RUNTIME_ADJFLOW_SYS) || (RunTime_EqSystem == RUNTIME_LINFLOW_SYS) 
			|| (RunTime_EqSystem == RUNTIME_FLOW_ADJFLOW_SYS) || (RunTime_EqSystem == RUNTIME_MULTIGRID_SYS) 
			|| (RunTime_EqSystem == RUNTIME_PLASMA_SYS) || (RunTime_EqSystem == RUNTIME_LEVELSET_SYS)) {

		if ((config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) && (RunTime_EqSystem != RUNTIME_LEVELSET_SYS))
			solution_container[MainSolution]->SetResidual_DualTime(geometry, solution_container, config, iRKStep, iMesh);

		solution_container[MainSolution]->SetResidual_Total(geometry, solution_container, config, iRKStep, iMesh);
	}

}

void CIntegration::Time_Integration(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep, 
		unsigned short RunTime_EqSystem, unsigned long Iteration) {
	unsigned short MainSolution = config->GetContainerPosition(RunTime_EqSystem);
	bool restart = (config->GetRestart() || config->GetRestart_Flow());

	
	/*--- In case we are performing a restart, the first iteration doesn't count. This is very important in case
	 we are performing a parallel simulation because at the first iterations gradients, undivided laplacian, etc.
	 are not updated at the send-receive boundaries ---*/
	if ((!restart) || (restart && Iteration != 0)) {

		/*--- Perform a residual smoothing technique ---*/
		if ((RunTime_EqSystem != RUNTIME_TURB_SYS) && (RunTime_EqSystem != RUNTIME_ADJFLOW_SYS) 
				&& (RunTime_EqSystem != RUNTIME_LINFLOW_SYS) && (RunTime_EqSystem != RUNTIME_LEVELSET_SYS) )
			solution_container[MainSolution]->SetResidual_Smoothing (geometry, config->GetnSmooth(), 
					config->GetSmoothCoeff());

		/*--- Perform the time integration ---*/
		switch (config->GetKind_TimeIntScheme()) {
			case (RUNGE_KUTTA_EXPLICIT):
				solution_container[MainSolution]->SetResidual_RKCoeff(geometry, config, iRKStep);
				solution_container[MainSolution]->RungeKutta_Iteration(geometry, solution_container, config, iRKStep);
				break;
			case (EULER_EXPLICIT):
				solution_container[MainSolution]->ExplicitEuler_Iteration(geometry, solution_container, config);
				break;
			case (EULER_IMPLICIT):
				solution_container[MainSolution]->ImplicitEuler_Iteration(geometry, solution_container, config);
				break;
		}
	}
}

void CIntegration::Solving_Linear_System(CGeometry *geometry, CSolution *solution, CSolution **solution_container, CConfig *config, 
		unsigned short iMesh) {

	/*--- Compute the solution of the linear system ---*/
	solution->Solve_LinearSystem(geometry, solution_container, config, iMesh);

	/*--- Compute the residual of the linear system ---*/
	solution->Compute_Residual(geometry, solution_container, config, iMesh);

}

void CIntegration::Convergence_Monitoring(CGeometry *geometry, CConfig *config, unsigned long Iteration, double monitor) {
	unsigned short iCounter;

	bool Already_Converged = Convergence;
	
	if (config->GetConvCriteria() == CAUCHY) {
		if (Iteration  == 0) {
			Cauchy_Value = 0;
			Cauchy_Counter = 0;
			for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
				Cauchy_Serie[iCounter] = 0.0;
		}

		Old_Func = New_Func;
		New_Func = monitor;
		Cauchy_Func = fabs(New_Func - Old_Func);

		Cauchy_Serie[Cauchy_Counter] = Cauchy_Func;
		Cauchy_Counter++;

		if (Cauchy_Counter == config->GetCauchy_Elems()) Cauchy_Counter = 0;

		Cauchy_Value = 1;
		if (Iteration  >= config->GetCauchy_Elems()) {
			Cauchy_Value = 0;
			for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
				Cauchy_Value += Cauchy_Serie[iCounter];
		}

		if (Cauchy_Value >= config->GetCauchy_Eps()) Convergence = false;
		else Convergence = true;

		if (Cauchy_Value >= config->GetCauchy_Eps_OneShot()) Convergence_OneShot = false;
		else Convergence_OneShot = true;

		if (Cauchy_Value >= config->GetCauchy_Eps_FullMG()) Convergence_FullMG = false;
		else Convergence_FullMG = true;
	}

	if (Iteration > config->GetStartConv_Iter()) {
		if (config->GetConvCriteria() == RESIDUAL) {
			if (Iteration == 1 + config->GetStartConv_Iter() ) InitResidual = monitor;

			if (monitor > InitResidual) InitResidual = monitor;

			if (Iteration   > 1 + config->GetStartConv_Iter() ) {
				if (((abs(InitResidual - monitor) >= config->GetOrderMagResidual()) && (monitor < InitResidual))  ||
						(monitor <= config->GetMinLogResidual())) Convergence = true;
				else Convergence = false;

			}
			else Convergence = false;

		}
	}
	else {
		Convergence = false;
		Convergence_OneShot = false;
		Convergence_FullMG = false;

	}

	if (Already_Converged) Convergence = true;
}

void CIntegration::SetDualTime_Solver(CGeometry *geometry, CSolution *flow) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		flow->node[iPoint]->Set_Solution_time_n1();
		flow->node[iPoint]->Set_Solution_time_n();

		geometry->node[iPoint]->SetVolume_n1();
		geometry->node[iPoint]->SetVolume_n();
	}
}

void CIntegration::SetTwoPhase_Solver(CGeometry **geometry, CSolution ***solution_container, CNumerics ****solver_container, CConfig *config,
		unsigned long Iteration) {
	unsigned long iPoint;
	unsigned short iMesh;
	double DensityInc, ViscosityInc, Heaviside, LevelSet, lambda, Density, 
	Pressure, *solution, yFreeSurface, PressFreeSurface, DensityFreeSurface, 
	Froude, yCoord;

	bool gravity = config->GetGravityForce();
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	
	for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
		for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
			LevelSet = solution_container[iMesh][LEVELSET_SOL]->node[iPoint]->GetSolution(0);
			double epsilon = config->GetInterfaseThickness();
			Heaviside = 0.0;
			if (LevelSet < -epsilon) Heaviside = 1.0;
			if (abs(LevelSet) <= epsilon) Heaviside = 1.0 - (0.5*(1.0+(LevelSet/epsilon)+(1.0/PI_NUMBER)*sin(PI_NUMBER*LevelSet/epsilon)));
			if (LevelSet > epsilon) Heaviside = 0.0;

			/*--- Set the value of the incompressible density for two phase flows (density ratio g/l)---*/
			lambda = config->GetRatioDensity();
			DensityInc = (lambda + (1.0 - lambda)*Heaviside)*config->GetDensity_FreeStreamND();
			solution_container[iMesh][FLOW_SOL]->node[iPoint]->SetDensityInc(DensityInc);

			/*--- Set the value of the incompressible viscosity for two phase flows (viscosity ratio g/l)---*/
			lambda = config->GetRatioViscosity();
			ViscosityInc = (lambda + (1.0 - lambda)*Heaviside) / config->GetReynolds();
			solution_container[iMesh][FLOW_SOL]->node[iPoint]->SetLaminarViscosityInc(ViscosityInc);

			/*--- Set initial boundary condition at iter 0 ---*/
			if ((Iteration == 0) && (!restart)) {
				yFreeSurface = config->GetLevelSet_Zero();
				PressFreeSurface = solution_container[iMesh][FLOW_SOL]->GetPressure_Inf();
				DensityFreeSurface = solution_container[iMesh][FLOW_SOL]->GetDensity_Inf();
				Density = solution_container[iMesh][FLOW_SOL]->node[iPoint]->GetDensityInc();
				Froude = config->GetVelocity_Ref() / sqrt( config->GetLength_Ref() * STANDART_GRAVITY);
				yCoord = geometry[iMesh]->node[iPoint]->GetCoord(1);
				if (gravity)
					Pressure = Density*(PressFreeSurface/DensityFreeSurface + (yFreeSurface-yCoord)/(Froude*Froude));
				else
					Pressure = Density*(PressFreeSurface/DensityFreeSurface );

				solution = solution_container[iMesh][FLOW_SOL]->node[iPoint]->GetSolution();
				solution[0] = Pressure/DensityFreeSurface;
//				solution[0] = PressFreeSurface/DensityFreeSurface;
//				solution[0] = PressFreeSurface/Density;

				solution_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(solution);
			}
			
		}
	}
}

void CIntegration::SetInitialCondition(CGeometry **geometry, CSolution ***solution_container, CNumerics ****solver_container, CConfig *config,
																			unsigned long Iteration) {
	unsigned long iPoint;
	unsigned short iMesh;
	double Density, Pressure, *solution, yFreeSurface, PressFreeSurface, DensityFreeSurface, 
	Froude, yCoord;
	
	bool gravity = config->GetGravityForce();
	bool incompressible = config->GetIncompressible();
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	
	if (incompressible) {
		
		for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
			for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
				
				/*--- Set initial boundary condition at the first iteration ---*/
				if ((Iteration == 0) && (!restart)) {
					yFreeSurface = config->GetLevelSet_Zero();
					PressFreeSurface = solution_container[iMesh][FLOW_SOL]->GetPressure_Inf();
					DensityFreeSurface = solution_container[iMesh][FLOW_SOL]->GetDensity_Inf();
					Density = solution_container[iMesh][FLOW_SOL]->GetDensity_Inf();
					Froude = config->GetVelocity_Ref() / sqrt( config->GetLength_Ref() * STANDART_GRAVITY);
					yCoord = geometry[iMesh]->node[iPoint]->GetCoord(1);
					if (gravity)
						Pressure = Density*(PressFreeSurface/DensityFreeSurface + (yFreeSurface-yCoord)/(Froude*Froude));
					else
						Pressure = Density*(PressFreeSurface/DensityFreeSurface);
					
					solution = solution_container[iMesh][FLOW_SOL]->node[iPoint]->GetSolution();
					solution[0] = Pressure/config->GetDensity_FreeStreamND();
					
					solution_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(solution);
				}
			}
		}
	}
}

