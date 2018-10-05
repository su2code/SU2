/*!
 * \file driver_structure_checkpointing.cpp
 * \brief The main subroutines for driving single or multi-zone problems.
 * \author T. Kattmann
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

#include "../include/driver_structure.hpp"
#include "../include/definition_structure.hpp"

//#include "../include/revolve.hpp"
#include "../include/checkpointing.hpp"
#include <stdexcept>

//CDiscAdjFluidDriver::CDiscAdjFluidDriver() { }

//CDiscAdjFluidDriver::~CDiscAdjFluidDriver() { }

void CDiscAdjFluidDriver::PrimalAdvance()
{

  unsigned short iZone, jZone, iPoint;
  bool unsteady = config_container[ZONE_0]->GetUnsteady_Simulation() != STEADY;

  unsigned long IntIter, nIntIter;

  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  unsteady = (config_container[MESH_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[MESH_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  /*--- Zone preprocessing ---*/
  //added: this loads the restart files, in iteration_structure.cpp - CDiscAdjFluidIteration->Preprocess()
  //Here it has to be checked whether it is the beginning of an advance cycle and then for DualTime2nd just 2 steps has to be loaded
  //as it is a simple primal solver restart, maybe just set ExtIter to 0 at the beginning of every advance cycle
  if (ExtIter == 0)
  {
    for (iZone = 0; iZone < nZone; iZone++)
    {
      /*--- loads time steps for adjoint taping ---*/
      //iteration_container[iZone]->Preprocess(output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone);

      /*--- Initializes Primal with freestream values, not explicitly necessary as initialization is done with fressstream valuies anyway ---*/
      //cout << "Initialze Primal by using the values at infinity (freestream)." << endl;
      //solver_container[iZone][MESH_0][FLOW_SOL]->SetFreeStream_Solution(config_container[iZone]);
      //for (iPoint = 0; iPoint < geometry_container[ZONE_0][MESH_0]->GetnPoint(); iPoint++) {
      //solver_container[iZone][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
      //solver_container[iZone][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
      //}
    }
  }

  /*--- applies only to windgust and fsi ---*/
  for (iZone = 0; iZone < nZone; iZone++)
    direct_iteration[iZone]->Preprocess(output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

  /*--- Updating zone interface communication patterns,
   needed only for unsteady simulation since for steady problems
  this is done once in the interpolator_container constructor
   at the beginning of the computation ---*/

  if (unsteady)
  {
    for (iZone = 0; iZone < nZone; iZone++)
    {
      for (jZone = 0; jZone < nZone; jZone++)
        if (jZone != iZone && interpolator_container[iZone][jZone] != NULL)
          interpolator_container[iZone][jZone]->Set_TransferCoeff(config_container);
    }
  }

  nIntIter = config_container[MESH_0]->GetUnst_nIntIter();

  for (IntIter = 0; IntIter < nIntIter; IntIter++)
  {

    ///*--- Print cons[0] for the first 10 points for the timesteps n-1, n, n+1 ---*/
    //for(iZone = 0; iZone < nZone; iZone++) {
    ///*--- Print cons[0] for the first 10 points for the timesteps n-1, n, n+1 ---*/
    //cout << "Print Cons[0] for the first 10 points in the mesh. 1st Solution, 2nd Solution_n, 3rd Solution_n1 " << endl;
    //for (int iPoint = 1000; iPoint < 1010; iPoint++)
    //cout << setprecision(15) << scientific << solver_container[iZone][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution(0) << " ";
    //cout << endl;
    //for (int iPoint = 1000; iPoint < 1010; iPoint++)
    //cout << setprecision(15) << scientific << solver_container[iZone][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n(0) << " ";
    //cout << endl;
    //for (int iPoint = 1000; iPoint < 1010; iPoint++)
    //cout << setprecision(15) << scientific << solver_container[iZone][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1(0) << " ";
    //cout << endl;
    //}

    /*--- Do one iteration of the direct solver  --*/

    /*--- At each pseudo time-step updates transfer data ---*/
    for (iZone = 0; iZone < nZone; iZone++)
      for (jZone = 0; jZone < nZone; jZone++)
        if (jZone != iZone && transfer_container[iZone][jZone] != NULL)
          Transfer_Data(iZone, jZone);

    /*--- For each zone runs one single iteration ---*/

    for (iZone = 0; iZone < nZone; iZone++)
    {

      config_container[iZone]->SetIntIter(IntIter);
      direct_iteration[iZone]->Iterate(output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);
      //cout << "Jacobian set to zero." << endl;
      //See CEuler/NSSolver::Preprocessing, jacobian is not set to zero as it is needed for taping.
      solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->Jacobian.SetValZero();
    }

    ///*--- Print residuals in the first iteration ---*/
    //for (iZone = 0; iZone < nZone; iZone++) {
    //if (rank == MASTER_NODE && ((ExtIter == 0) || (config_container[iZone]->GetUnsteady_Simulation() != STEADY)) ) {
    //cout << "At IntIter: " << IntIter;
    //cout << " Zone " << iZone << ": log10[Conservative 0]: "<< log10(solver_container[iZone][MESH_0][FLOW_SOL]->GetRes_RMS(0)) << endl;
    //if ( config_container[iZone]->GetKind_Turb_Model() != NONE && !config_container[iZone]->GetFrozen_Visc_Disc()) {
    //cout <<"       log10[RMS k]: " << log10(solver_container[iZone][MESH_0][TURB_SOL]->GetRes_RMS(0)) << endl;
    //}
    //}
    //}
  }
}

void CDiscAdjFluidDriver::Update()
{

  /*--- Advance fowrard in time for Primal steps ---*/
  if (rank == MASTER_NODE && false)
    cout << "void CDiscAdjFluidDriver::Update()" << endl;
  if (config_container[ZONE_0]->GetCheckpointing())
  {
    for (iZone = 0; iZone < nZone; iZone++)
      direct_iteration[iZone]->Update(output, integration_container, geometry_container,
                                      solver_container, numerics_container, config_container,
                                      surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }
}

void CDiscAdjFluidDriver::StoreSingleState()
{
  int iPoint, iMesh;
  bool turbulent = (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_RANS);
  int info = 3;

  /*--- Compute 1st timestep at the start of the primal solver to store a full checkpoint (if current state needs to be stored as well, i.e. DT_2nd stores 3 states)  ---*/
  if (r->getcurrent_timestep() == r->gettimestepping_order())
  {
    PreprocessExtIter(r->getcurrent_timestep());
    /*--- Advance one physical primal step by looping over multiple internal iter. But no Update()/solution pushing. ---*/
    PrimalAdvance();
  }

  if (r->getwhere())
  {
    /*--- Save Checkpoint in RAM ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_CheckpointSingleState(r->getcheck()); //TODO getcheck() // allocated checkpoints and ID do not fit together propably only single state checkpoints is the best solution
      if (turbulent)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_CheckpointSingleState(r->getcheck());
      }
    }
  } else {
    /*--- takeshot on DISK, 3 timesteps have to be saved here! Storing wrong Primitives for n and n1. Storage always from sol position (therefore solution_old as intermediate container) ---*/

    /*--- Store Sol, setExtIter necessary for output with correct number  ---*/
    config_container[ZONE_0]->SetExtIter(r->getcurrent_timestep());
    config_container[ZONE_0]->SetKind_Solver(NAVIER_STOKES);
    if (turbulent)
      config_container[ZONE_0]->SetKind_Solver(RANS);
    config_container[ZONE_0]->SetDiscrete_Adjoint(false);

    output->SetResult_Files_Parallel(solver_container, geometry_container, config_container, r->getcurrent_timestep(), nZone, nInst);

    config_container[ZONE_0]->SetKind_Solver(DISC_ADJ_NAVIER_STOKES);
    config_container[ZONE_0]->SetDiscrete_Adjoint(true);
    if (turbulent)
      config_container[ZONE_0]->SetKind_Solver(DISC_ADJ_RANS);
  }

  if (info > 1)
    if (rank == MASTER_NODE)
      cout << "Store single state at " << setw(6) << r->getcurrent_timestep();
  if (r->getwhere())
    if (rank == MASTER_NODE)
      cout << " in RAM" << " in CP " << r->getcheck() << endl;
    else if (rank == MASTER_NODE)
      cout << " in ROM." << endl;
}

/* A restore single state is basically what happens currently during the adjoint run */
/* Only tested for DT2nd */
void CDiscAdjFluidDriver::RestoreSingleState()
{

  int iPoint, iMesh;
  bool turbulent = (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_RANS);
  int info = 3;
  bool rans = ((config_container[ZONE_0]->GetKind_Solver() == RANS) ||
               (config_container[ZONE_0]->GetKind_Solver() == ADJ_RANS) ||
               (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_RANS));

  if (r->getwhere())
  {
    /*--- 1 Load checkpoint (into sol). 2 store sol to sol old. 3 push n to sol. 4 push n1 to n. 5 Put old into n1 ---*/

    /*--- 1 Load checkpoint (into sol) ---*/
    /*--- Restore Checkpoint (conservatives) from RAM, with own routine ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Restore_CheckpointSingleState(r->getcheck());
      if (turbulent)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Restore_CheckpointSingleState(r->getcheck());
      }
    }

    /*--- 2 store sol to sol old ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_OldSolution();
      if (rans)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_OldSolution();
      }
    }
    /*--- 3 push n to sol ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n());

      if (rans)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->SetSolution(solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->GetSolution_time_n());
      }
    }
    /*--- 4 push n1 to n. ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1());
      if (rans)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_Solution_time_n(solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->GetSolution_time_n1());
      }
    }
    /*--- 5 Put old into n1 ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_Old());
      if (rans)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_Solution_time_n1(solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->GetSolution_Old());
      }
    }
  }
  else // get checkpoint from disk
  {
    /*--- 1 Load checkpoint (into sol). 2 store sol to sol old. 3 push n to sol. 4 push n1 to n. 5 Put old into n1 ---*/

    /*--- 1 Load checkpoint (into sol) ---*/
    //last variable updategeo always true except for FSI
    //int val_iter = 1000 + r->getcheck() * 3 + 0;
    int val_iter = r->getcurrent_timestep() - 3; // -3 du to dual time stepping 2nd order
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[ZONE_0][INST_0], solver_container[ZONE_0][INST_0], config_container[ZONE_0], val_iter, true);
    if (turbulent)
      solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->LoadRestart(geometry_container[ZONE_0][INST_0], solver_container[ZONE_0][INST_0], config_container[ZONE_0], val_iter, true);

    /*--- 2 store sol to sol old ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_OldSolution();
      if (rans)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_OldSolution();
      }
    }
    /*--- 3 push n to sol ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n());

      if (rans)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->SetSolution(solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->GetSolution_time_n());
      }
    }
    /*--- 4 push n1 to n. ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1());
      if (rans)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_Solution_time_n(solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->GetSolution_time_n1());
      }
    }
    /*--- 5 Put old into n1 ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_Old());
      if (rans)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_Solution_time_n1(solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->GetSolution_Old());
      }
    }

  } // end of if-else

  if (info > 2 && rank == MASTER_NODE) {
    cout << "Restore single state at " << setw(7) << r->getcurrent_timestep() - 3;
    if (r->getwhere()) {
      cout << "restore in RAM " << "in CP " << r->getcheck() << endl;
    } else {
      cout << "restore on Disk" << endl;
    }
  }
}

void CDiscAdjFluidDriver::StoreFullCheckpoint()
{

  int iPoint, iMesh;
  bool turbulent = (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_RANS);
  int info = 3;

  /*--- Compute 1st timestep at the start of the primal solver to store a full checkpoint (if current state needs to be stored as well, i.e. DT_2nd stores 3 states)  ---*/
  if ( r->getcurrent_timestep() == r->gettimestepping_order() )
  {
    if(rank == MASTER_NODE) std::cout << "Computing first Primal Step before taking full checkpoint." << std::endl;
    PreprocessExtIter(0);
    /*--- Advance one physical primal step by looping over multiple internal iter. But no Update()/solution pushing. ---*/
    PrimalAdvance();
  }

  if (r->getwhere() && false) // Here It needs to stored in single checkpoints only !
  {
    /*--- Save Checkpoint in RAM ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Checkpoint(r->getcheck());
      if (turbulent)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_Checkpoint(r->getcheck());
      }
    }
  } else if (r->getwhere() && true) {
    std::cout << "In new Store Full Checkpoint RAM method" << std::endl;
    /*--- Save Checkpoints in RAM: 1. Store n1 (if DT2nd) in CP(i) 2. Store n in CP(i+1) 3. Store current sol in CP(i+2) ---*/
    
    /* 1. Store n1 (if DT2nd) in CP(i) */
    // if(DT2nd)
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_CheckpointSingleState_time_n1(r->getcheck()-2);
      if (turbulent)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_CheckpointSingleState_time_n1(r->getcheck()-2);
      }
    }
    
    /* 2. Store n in CP(i+1) */
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_CheckpointSingleState_time_n(r->getcheck()-1);
      if (turbulent)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_CheckpointSingleState_time_n(r->getcheck()-1);
      }
    }
    
    /* 3. Store current sol in CP(i+2) */
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_CheckpointSingleState(r->getcheck());
      if (turbulent)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_CheckpointSingleState(r->getcheck());
      }
    }
    
  } else { //Disk
    /*--- takeshot on DISK, 3 timesteps have to be saved here! Storing wrong Primitives for n and n1. Storage always from sol position (therefore solution_old as intermediate container) ---*/

    /*--- Store Sol, setExtIter necessary for output with correct number  ---*/
    config_container[ZONE_0]->SetExtIter(r->getcurrent_timestep() + 0);
    config_container[ZONE_0]->SetKind_Solver(NAVIER_STOKES);
    if (turbulent)
      config_container[ZONE_0]->SetKind_Solver(RANS);
    config_container[ZONE_0]->SetDiscrete_Adjoint(false);
    output->SetResult_Files_Parallel(solver_container, geometry_container, config_container, r->getcurrent_timestep() + 0, nZone, nInst);
    config_container[ZONE_0]->SetKind_Solver(DISC_ADJ_NAVIER_STOKES);
    config_container[ZONE_0]->SetDiscrete_Adjoint(true);
    if (turbulent)
      config_container[ZONE_0]->SetKind_Solver(DISC_ADJ_RANS);

    /*--- store current solution in sol_old (to not kill the advance cycle), pushup n ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_OldSolution();
      if (turbulent)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_OldSolution();
      }
    }

    /*--- Push n to sol ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n());
      if (turbulent)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->SetSolution(solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->GetSolution_time_n());
      }
    }

    /*--- Store Sol (i.e. n) ---*/
    config_container[ZONE_0]->SetExtIter(r->getcurrent_timestep() - 1);
    config_container[ZONE_0]->SetKind_Solver(NAVIER_STOKES);
    config_container[ZONE_0]->SetDiscrete_Adjoint(false);
    if (turbulent)
      config_container[ZONE_0]->SetKind_Solver(RANS);
    output->SetResult_Files_Parallel(solver_container, geometry_container, config_container, r->getcurrent_timestep() - 1, nZone, nInst);
    config_container[ZONE_0]->SetKind_Solver(DISC_ADJ_NAVIER_STOKES);
    config_container[ZONE_0]->SetDiscrete_Adjoint(true);
    if (turbulent)
      config_container[ZONE_0]->SetKind_Solver(DISC_ADJ_RANS);

    /*--- pushup n1 directly to sol, omitting n which stays on the right place ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1());
      if (turbulent)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->SetSolution(solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->GetSolution_time_n1());
      }
    }

    /*--- Store Sol ---*/
    config_container[ZONE_0]->SetExtIter(r->getcurrent_timestep() - 2);
    config_container[ZONE_0]->SetKind_Solver(NAVIER_STOKES);
    config_container[ZONE_0]->SetDiscrete_Adjoint(false);
    if (turbulent)
      config_container[ZONE_0]->SetKind_Solver(RANS);
    output->SetResult_Files_Parallel(solver_container, geometry_container, config_container, r->getcurrent_timestep() - 2, nZone, nInst);
    config_container[ZONE_0]->SetKind_Solver(DISC_ADJ_NAVIER_STOKES);
    config_container[ZONE_0]->SetDiscrete_Adjoint(true);
    if (turbulent)
      config_container[ZONE_0]->SetKind_Solver(DISC_ADJ_RANS);

    /*--- restore Sol from sol_old ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_Old());
      if (turbulent)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->SetSolution(solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->GetSolution_Old());
      }
    }
  }
  if (info > 1 && rank == MASTER_NODE) {
    cout << "Store full checkpoint at " << setw(6) << r->getcurrent_timestep();
      
    if (r->getwhere()) {
      cout << " in RAM " << " in CP " << r->getcheck() << "." << endl;
    } else { 
      cout << " in ROM" << "." << endl;
    }
  }
}

void CDiscAdjFluidDriver::PrimalStep()
{
  int info = 3;
  int val_iter = r->getcurrent_timestep();
  if (rank == MASTER_NODE) cout << "ADVANCE: " << val_iter << endl;
  /*--- Here Flow-solver initial conditions are called, well noe not anymore as first step is taken at takeshot ---*/
  PreprocessExtIter(val_iter);

  /*--- Pushes soltion down for a new timestep. Therefore this must be done after takeshot. Directly before takeshot however update is omitted. ---*/
  //if (j == r->getoldcapo() + 1)
  //  Update();

  /*--- Advance one physical primal step by looping over multiple internal iter ---*/
  PrimalAdvance();

  /*--- No Update as solution needs to at their place for takeshot ---*/
  //if (j != r->getcapo() + 1 - 1)
  //  Update();

  //Monitor(j);

  if (info > 2 && rank == MASTER_NODE)
   cout << "Primal step at " << setw(7) << r->getcurrent_timestep() << endl;
}

void CDiscAdjFluidDriver::RestoreFullCheckpoint()
{
  int info = 3;
  int iPoint;
  bool turbulent = (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_RANS);
  bool rans = ((config_container[ZONE_0]->GetKind_Solver() == RANS) ||
               (config_container[ZONE_0]->GetKind_Solver() == ADJ_RANS) ||
               (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_RANS));
  if (r->getwhere())
  {

    /*--- Restore Checkpoint (conservatives) from RAM, with own routine ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Restore_Checkpoint(r->getcheck());
      if (turbulent)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Restore_Checkpoint(r->getcheck());
      }
    }
  }
  else
  {
    /*--- Load the three restart_files, and put them at legit positions ---*/
    /*--- Load n1, pushback twice, load n, pushback once, load sol  ---*/

    /*--- Load n1 ---*/
    //last variable updategeo always true except for FSI
    int val_iter = r->getcurrent_timestep() - 2;
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[ZONE_0][INST_0], solver_container[ZONE_0][INST_0], config_container[ZONE_0], val_iter, true);
    if (turbulent)
      solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->LoadRestart(geometry_container[ZONE_0][INST_0], solver_container[ZONE_0][INST_0], config_container[ZONE_0], val_iter, true);

    /*--- pushback twice ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
      if (rans)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_Solution_time_n1();
      }
    }

    /*--- Load n ---*/
    val_iter = r->getcurrent_timestep() - 1;
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[ZONE_0][INST_0], solver_container[ZONE_0][INST_0], config_container[ZONE_0], val_iter, true);
    if (turbulent)
      solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->LoadRestart(geometry_container[ZONE_0][INST_0], solver_container[ZONE_0][INST_0], config_container[ZONE_0], val_iter, true);

    /*--- pushback once ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++)
    {
      solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
      if (rans)
      {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
      }
    }

    /*--- Load sol ---*/
    val_iter = r->getcurrent_timestep() - 0;
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[ZONE_0][INST_0], solver_container[ZONE_0][INST_0], config_container[ZONE_0], val_iter, true);
    if (turbulent)
      solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->LoadRestart(geometry_container[ZONE_0][INST_0], solver_container[ZONE_0][INST_0], config_container[ZONE_0], val_iter, true);
  }
  if (info > 2 && rank == MASTER_NODE) {
    cout << " Restore full checkpoint at timestep " << setw(7) << r->getcurrent_timestep() - 0;
      
    if (r->getwhere()) {
      cout << " in RAM" << " in CP " << r->getcheck() << endl;
    } else { 
      cout << " in ROM " << endl;
    }
  }
}
void CDiscAdjFluidDriver::PrimalUpdate()
{
  Update();
  int info = 3;
  if (info > 2 && rank == MASTER_NODE) {
    std::cout << "Primal update at timestep " << r->getcurrent_timestep() << std::endl;
  }
}

void CDiscAdjFluidDriver::AdjointStep()
{
  int info = 3;
  /*--- Do one adjoint step ---*/
  ExtIter = r->getsteps() - r->getcurrent_timestep() - 1;
  PreprocessExtIter(ExtIter);
  DynamicMeshUpdate(ExtIter);
  Run();
  /*--- Needs to be set for correct primal advance ---*/
  solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->Jacobian.SetValZero();
  Monitor(ExtIter);
  Output(ExtIter);
  ExtIter++;

  if (info > 2 && rank == MASTER_NODE)
    cout << "Adjoint Step at timestep " << setw(7) << r->getcurrent_timestep() << endl;
}

/* Equivalent to adjoint step, could probably be completely removed */
void CDiscAdjFluidDriver::Firsturn()
{
  int info = 3;
  /*--- Do first Adjoint Step. As it comes after advance without Update() the correct primals are already set. ---*/
  PreprocessExtIter(ExtIter);
  DynamicMeshUpdate(ExtIter);
  Run();
  /*--- Otherwise wrong primal solution after restore and 1 iteration if Jacobian not set to zero ---*/
  solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->Jacobian.SetValZero();
  /*--- Update(); no update for adjoint solver only for primal, therefore update before PrimalAdvance. Primal Solutions are always restored completely, so no updating necessary---*/
  Monitor(ExtIter);
  Output(ExtIter);
  ExtIter++;

  if (info > 2 && rank == MASTER_NODE)
    cout << " THIS SHOULD NOT RUN at " << setw(6) << r->getcapo() << endl;
}
