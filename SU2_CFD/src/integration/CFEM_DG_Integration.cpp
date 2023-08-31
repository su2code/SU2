/*!
 * \file CFEM_DG_Integration.cpp
 * \brief Definition of time and space integration for the DG solver.
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/integration/CFEM_DG_Integration.hpp"


CFEM_DG_Integration::CFEM_DG_Integration() : CIntegration() { }

void CFEM_DG_Integration::SingleGrid_Iteration(CGeometry ****geometry,
                                               CSolver *****solver_container,
                                               CNumerics ******numerics_container,
                                               CConfig **config,
                                               unsigned short RunTime_EqSystem,
                                               unsigned short iZone,
                                               unsigned short iInst) {

  unsigned short iMesh, iStep, iLimit = 1;
  unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);
  unsigned short FinestMesh = config[iZone]->GetFinestMesh();

  /*--- For now, we assume no geometric multigrid. ---*/
  iMesh = FinestMesh;

  /*--- Check if only the Jacobian of the spatial discretization must
        be computed. If so, call the appropriate function and return. ---*/
  if (config[iZone]->GetJacobian_Spatial_Discretization_Only()) {
    solver_container[iZone][iInst][iMesh][SolContainer_Position]->ComputeSpatialJacobian(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                                                                                         numerics_container[iZone][iInst][iMesh][SolContainer_Position],
                                                                                         config[iZone], iMesh, RunTime_EqSystem);
    return;
  }

  /*--- Determine the number of stages in the time stepping algorithm.
        For the Runge-Kutta schemes this is the number of RK stages,
        while for ADER-DG this information is not used, because a more
        complicated algorithm must be used to facilitate time accurate
        local time stepping.  Note that we are currently hard-coding
        the classical RK4 scheme. ---*/
  bool useADER = false;
  switch (config[iZone]->GetKind_TimeIntScheme()) {
    case RUNGE_KUTTA_EXPLICIT: iLimit = config[iZone]->GetnRKStep(); break;
    case CLASSICAL_RK4_EXPLICIT: iLimit = 4; break;
    case ADER_DG: iLimit = 1; useADER = true; break;
    case EULER_EXPLICIT: case EULER_IMPLICIT: iLimit = 1; break; }

  /*--- In case an unsteady simulation is carried out, it is possible that a
        synchronization time step is specified. If so, set the boolean
        TimeSynSpecified to true, which leads to an outer loop in the
        algorithm below. ---*/
  bool TimeSyncSpecified   = false;
  const su2double TimeSync = config[iZone]->GetTime_Step()/config[iZone]->GetTime_Ref();
  if(config[iZone]->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING &&
     config[iZone]->GetUnst_CFL()            != 0.0           &&
     TimeSync                                != 0.0) TimeSyncSpecified = true;

  /*--- Outer loop, which is only active when a synchronization time has been
        specified for an unsteady simulation. ---*/
  bool syncTimeReached = false;
  su2double timeEvolved = 0.0;
  while( !syncTimeReached ) {

    /* Compute the time step for stability. */
    solver_container[iZone][iInst][iMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iInst][iMesh],
                                                                               solver_container[iZone][iInst][iMesh],
                                                                               config[iZone], iMesh, config[iZone]->GetTimeIter());
    /* Possibly overrule the specified time step when a synchronization time was
       specified and determine whether or not the time loop must be continued.
       When TimeSyncSpecified is false, the loop is always terminated. */
    if( TimeSyncSpecified )
      solver_container[iZone][iInst][iMesh][SolContainer_Position]->CheckTimeSynchronization(config[iZone],
                                                                                             TimeSync, timeEvolved,
                                                                                             syncTimeReached);
    else
      syncTimeReached = true;

    /*--- For ADER in combination with time accurate local time stepping, the
          space and time integration are tightly coupled and cannot be treated
          segregatedly. Therefore a different function is called for ADER to
          carry out the space and time integration. ---*/
    if( useADER ) {
      solver_container[iZone][iInst][iMesh][SolContainer_Position]->ADER_SpaceTimeIntegration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                                                                                              numerics_container[iZone][iInst][iMesh][SolContainer_Position],
                                                                                              config[iZone], iMesh, RunTime_EqSystem);
    }
    else {

      /*--- Time and space integration can be decoupled. ---*/
      for (iStep = 0; iStep < iLimit; iStep++) {

        /*--- Preprocessing ---*/
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                                                                                    config[iZone], iMesh, iStep, RunTime_EqSystem, false);

        /*--- Space integration ---*/
        Space_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                          numerics_container[iZone][iInst][iMesh][SolContainer_Position],
                          config[iZone], iMesh, iStep, RunTime_EqSystem);

        /*--- Time integration, update solution using the old solution plus the solution increment ---*/
        Time_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                         config[iZone], iStep, RunTime_EqSystem);

        /*--- Postprocessing ---*/
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                                                                                     config[iZone], iMesh);
      }
    }
  }

  /*--- Calculate the inviscid and viscous forces ---*/
  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Pressure_Forces(geometry[iZone][iInst][iMesh], config[iZone]);

  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Friction_Forces(geometry[iZone][iInst][iMesh], config[iZone]);

  /*--- Convergence strategy ---*/

  //Convergence_Monitoring(geometry[iZone][iInst][FinestMesh], config[iZone], Iteration, monitor, FinestMesh);
}

void CFEM_DG_Integration::Space_Integration(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CNumerics **numerics,
                                            CConfig *config, unsigned short iMesh,
                                            unsigned short iStep,
                                            unsigned short RunTime_EqSystem) {

  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);

  /*--- Runge-Kutta type of time integration schemes. In the first step, i.e.
        if iStep == 0, set the old solution (working solution for the DG part),
        and if needed, the new solution. ---*/
  if (iStep == 0) {
    solver_container[MainSolver]->Set_OldSolution();

    if (config->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT) {
      solver_container[MainSolver]->Set_NewSolution();
    }
  }

  /*--- Compute the spatial residual by processing the task list. ---*/
  solver_container[MainSolver]->ProcessTaskList_DG(geometry, solver_container, numerics, config, iMesh);
}

void CFEM_DG_Integration::Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iStep,
                                    unsigned short RunTime_EqSystem) {

  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);

  /*--- Perform the time integration ---*/
  switch (config->GetKind_TimeIntScheme()) {
    case (RUNGE_KUTTA_EXPLICIT):
      solver_container[MainSolver]->ExplicitRK_Iteration(geometry, solver_container, config, iStep);
      break;
    case (CLASSICAL_RK4_EXPLICIT):
      solver_container[MainSolver]->ClassicalRK4_Iteration(geometry, solver_container, config, iStep);
      break;
    default:
      SU2_MPI::Error("Time integration scheme not implemented.", CURRENT_FUNCTION);
  }
}
