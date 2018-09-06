/*!
 * \file driver_structure.cpp
 * \brief The main subroutines for driving single-zone problems.
 * \author R. Sanchez
 * \version 6.0.1 "Falcon"
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

#include "../include/driver_structure.hpp"
#include "../include/definition_structure.hpp"


CSinglezoneDriver::CSinglezoneDriver(char* confFile,
                       unsigned short val_nZone,
                       unsigned short val_nDim,
                       bool val_periodic,
                       SU2_Comm MPICommunicator) : CDriver(confFile,
                                                          val_nZone,
                                                          val_nDim,
                                                          val_periodic,
                                                          MPICommunicator) {

  /*--- Initialize the counter for TimeIter ---*/
  TimeIter = 0;

}

CSinglezoneDriver::~CSinglezoneDriver(void) {

}

void CSinglezoneDriver::StartSolver() {

  /*--- Main external loop of the solver. Runs for the number of time steps required. ---*/

  if (rank == MASTER_NODE){
    cout << endl <<"------------------------ Begin Single-zone Solver ------------------------" << endl;
    if (driver_config->GetTime_Domain())
      cout << endl <<"The simulation will run for " << driver_config->GetnTime_Iter() << " time steps." << endl;
  }

  /*--- Set the initial time iteration to the restart iteration. ---*/
  if (config_container[ZONE_0]->GetRestart())
    TimeIter = config_container[ZONE_0]->GetRestart_Iter();

  /*--- Run the problem until the number of time iterations required is reached. ---*/
  while ( TimeIter < config_container[ZONE_0]->GetnTime_Iter() ) {

    /*--- Perform some preprocessing before starting the time-step simulation. ---*/

    Preprocess(TimeIter);

    /*--- Run a time-step iteration of the single-zone problem. ---*/

    Run();

    /*--- Update the solution for dual time stepping strategy ---*/

    Update();

    /*--- Monitor the computations after each iteration. ---*/

    Monitor(TimeIter);

    /*--- Output the solution in files. ---*/

    Output(TimeIter);

    /*--- If the convergence criteria has been met, terminate the simulation. ---*/

    if (StopCalc) break;

    TimeIter++;

  }

}

void CSinglezoneDriver::Preprocess(unsigned long TimeIter) {

  bool unsteady = driver_config->GetTime_Domain();

  /*--- Set the value of the external iteration to TimeIter. -------------------------------------*/
  /*--- TODO: This should be generalised for an homogeneous criteria throughout the code. --------*/
  config_container[ZONE_0]->SetExtIter(TimeIter);

  /*--- Read the target pressure for inverse design. ---------------------------------------------*/
  /*--- TODO: This routine should be taken out of output, and made general for multiple zones. ---*/
  if (config_container[ZONE_0]->GetInvDesign_Cp() == YES)
    output->SetCp_InverseDesign(solver_container[iZone][INST_0][MESH_0][FLOW_SOL],
        geometry_container[ZONE_0][INST_0][MESH_0], config_container[iZone], TimeIter);

  /*--- Read the target heat flux ----------------------------------------------------------------*/
  /*--- TODO: This routine should be taken out of output, and made general for multiple zones. ---*/
  if (config_container[ZONE_0]->GetInvDesign_HeatFlux() == YES)
    output->SetHeatFlux_InverseDesign(solver_container[iZone][INST_0][MESH_0][FLOW_SOL],
        geometry_container[ZONE_0][INST_0][MESH_0], config_container[iZone], TimeIter);

  /*--- Set the initial condition for EULER/N-S/RANS ---------------------------------------------*/
  if ((config_container[ZONE_0]->GetKind_Solver() ==  EULER) ||
      (config_container[ZONE_0]->GetKind_Solver() ==  NAVIER_STOKES) ||
      (config_container[ZONE_0]->GetKind_Solver() ==  RANS) ) {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[ZONE_0][INST_0], solver_container[ZONE_0][INST_0], config_container[ZONE_0], TimeIter);
  }

#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  /*--- Run a predictor step ---*/
  if (config_container[ZONE_0]->GetPredictor())
    iteration_container[ZONE_0][INST_0]->Predictor(output, integration_container, geometry_container, solver_container,
        numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Perform a dynamic mesh update if required. ---*/

  DynamicMeshUpdate(TimeIter);

}

void CSinglezoneDriver::Run() {

  bool UpdateMesh = false;
  unsigned long ExtIter = 0;
  bool Convergence = false;

  unsigned long OuterIter = 0;
  config_container[ZONE_0]->SetOuterIter(OuterIter);

  /*--- Iterate the zone as a block, either to convergence or to a max number of iterations ---*/
  iteration_container[ZONE_0][INST_0]->Iterate_Block(output, integration_container, geometry_container, solver_container,
        numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- A corrector step can help preventing numerical instabilities ---*/
  Corrector();


}

void CSinglezoneDriver::Corrector() {

    if (config_container[ZONE_0]->GetRelaxation())
      iteration_container[ZONE_0][INST_0]->Relaxation(output, integration_container, geometry_container, solver_container,
          numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

}


void CSinglezoneDriver::Update() {

  unsigned short jZone;
  bool UpdateMesh;

  iteration_container[ZONE_0][INST_0]->Update(output, integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

}

void CSinglezoneDriver::DynamicMeshUpdate(unsigned long ExtIter) {

  bool harmonic_balance;

   harmonic_balance = (config_container[iZone]->GetUnsteady_Simulation() == HARMONIC_BALANCE);
    /*--- Dynamic mesh update ---*/
    if (config_container[ZONE_0]->GetGrid_Movement()) {
      iteration_container[ZONE_0][INST_0]->SetGrid_Movement(geometry_container, surface_movement, grid_movement, FFDBox, solver_container, config_container, ZONE_0, INST_0, 0, ExtIter );
    }

}

