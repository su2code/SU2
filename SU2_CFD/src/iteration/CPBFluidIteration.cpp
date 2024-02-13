/*!
 * \file CPBFluidIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 7.0.6 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/iteration/CPBFluidIteration.hpp"
#include "../../include/solvers/CPBIncEulerSolver.hpp"
#include "../../include/solvers/CPBIncNSSolver.hpp"
#include "../../include/output/COutput.hpp"

CPBFluidIteration::CPBFluidIteration(const CConfig *config) : CFluidIteration(config) { }

void CPBFluidIteration::Preprocess(COutput *output,
                                    CIntegration ****integration_container,
                                    CGeometry ****geometry_container,
                                    CSolver *****solver_container,
                                    CNumerics ******numerics_container,
                                    CConfig **config_container,
                                    CSurfaceMovement **surface_movement,
                                    CVolumetricMovement ***grid_movement,
                                    CFreeFormDefBox*** FFDBox,
                                    unsigned short val_iZone,
                                    unsigned short val_iInst) { }

void CPBFluidIteration::Iterate(COutput *output,
                                    CIntegration ****integration,
                                    CGeometry ****geometry,
                                    CSolver *****solver,
                                    CNumerics ******numerics,
                                    CConfig **config,
                                    CSurfaceMovement **surface_movement,
                                    CVolumetricMovement ***grid_movement,
                                    CFreeFormDefBox*** FFDBox,
                                    unsigned short val_iZone,
                                    unsigned short val_iInst) {
  unsigned long InnerIter, TimeIter;
  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  bool unsteady = (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) || (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  bool frozen_visc = (config[val_iZone]->GetContinuous_Adjoint() && config[val_iZone]->GetFrozen_Visc_Cont()) ||
                     (config[val_iZone]->GetDiscrete_Adjoint() && config[val_iZone]->GetFrozen_Visc_Disc());
  bool periodic = (config[val_iZone]->GetnMarker_Periodic() > 0);
  TimeIter = config[val_iZone]->GetTimeIter();
  
  /* --- Setting up iteration values depending on if this is a
   steady or an unsteady simulaiton */
  
  InnerIter = config[val_iZone]->GetInnerIter();

  /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations ---*/
  
  switch( config[val_iZone]->GetKind_Solver() ) {
    case MAIN_SOLVER::INC_EULER: 
      config[val_iZone]->SetGlobalParam(MAIN_SOLVER::INC_EULER, RUNTIME_FLOW_SYS); break;
      
    case MAIN_SOLVER::INC_NAVIER_STOKES:
      config[val_iZone]->SetGlobalParam(MAIN_SOLVER::INC_NAVIER_STOKES, RUNTIME_FLOW_SYS); break;
      
    case MAIN_SOLVER::INC_RANS:
      config[val_iZone]->SetGlobalParam(MAIN_SOLVER::INC_RANS, RUNTIME_FLOW_SYS); break;  
  }

  /*--- Solve the momentum equations. ---*/
  
  integration[val_iZone][val_iInst][FLOW_SOL]->CurrentGridIteration(geometry, solver, numerics, config, FinestMesh, RUNTIME_FLOW_SYS, val_iZone,val_iInst);


  /*--- Set source term for pressure correction equation based on current flow solution ---*/
  solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->SetMomCoeff(geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], periodic, MESH_0);
	
  solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->SetPoissonSourceTerm(geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], MESH_0);
	

  /*--- Solve the poisson equation ---*/
  config[val_iZone]->SetGlobalParam(MAIN_SOLVER::POISSON_EQUATION, RUNTIME_POISSON_SYS);
	
  integration[val_iZone][val_iInst][POISSON_SOL]->CurrentGridIteration(geometry, solver, numerics, config, FinestMesh, RUNTIME_POISSON_SYS, val_iZone,val_iInst);
//   integration[val_iZone][val_iInst][POISSON_SOL]->MultiGrid_Iteration(geometry, solver, numerics, config, RUNTIME_POISSON_SYS, ExtIter, val_iZone,val_iInst);
   
  /*--- Correct pressure and velocities ---*/
  solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->Flow_Correction(geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone]);
    
  /*--- Set the prmitive value based on updated solution ---*/
  solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->Postprocessing(geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], MESH_0);
  
  /*integration[val_iZone][val_iInst][FLOW_SOL]->MultiGrid_CyclePB(geometry, solver, numerics,
                                                                  config, FinestMesh, 0, RUNTIME_FLOW_SYS, val_iZone,val_iInst);*/
                                                                  
  if (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::INC_RANS) {
    
    /*--- Solve the turbulence model ---*/
    
    config[val_iZone]->SetGlobalParam(MAIN_SOLVER::INC_RANS, RUNTIME_TURB_SYS);
    integration[val_iZone][val_iInst][TURB_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config, RUNTIME_TURB_SYS, val_iZone, val_iInst);
    
    /*--- Solve transition model ---*/

    if (config[val_iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      config[val_iZone]->SetGlobalParam(MAIN_SOLVER::INC_RANS, RUNTIME_TRANS_SYS);
      integration[val_iZone][val_iInst][TRANS_SOL]->SingleGrid_Iteration(geometry, solver, numerics,
                                                                        config, RUNTIME_TRANS_SYS, val_iZone, val_iInst);
    }
  }
  

  /*--- Calculate the inviscid and viscous forces ---*/
      
  solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->Pressure_Forces(geometry[val_iZone][val_iInst][MESH_0], config[val_iZone]);
  solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->Momentum_Forces(geometry[val_iZone][val_iInst][MESH_0], config[val_iZone]);
  solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->Friction_Forces(geometry[val_iZone][val_iInst][MESH_0], config[val_iZone]);
    
    
  /*--- Write the convergence history ---*/

  /*if ( unsteady && !config[val_iZone]->GetDiscrete_Adjoint() ) 
	 output->SetConvHistory_Body(NULL, geometry, solver, config, integration, true, 0.0, val_iZone,val_iInst);*/
  
}
