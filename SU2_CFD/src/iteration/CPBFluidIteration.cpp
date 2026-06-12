/*!
 * \file CPBFluidIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../include/output/COutput.hpp"

void CPBFluidIteration::Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                 CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                 CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                 CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) { }

void CPBFluidIteration::Iterate(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                              CSolver***** solver, CNumerics****** numerics, CConfig** config,
                              CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                              CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  SU2_ZONE_SCOPED
  
  const bool unsteady = (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                        (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool frozen_visc = (config[val_iZone]->GetContinuous_Adjoint() && config[val_iZone]->GetFrozen_Visc_Cont()) ||
                           (config[val_iZone]->GetDiscrete_Adjoint() && config[val_iZone]->GetFrozen_Visc_Disc());
  const bool disc_adj = (config[val_iZone]->GetDiscrete_Adjoint());

  const unsigned short nCorrections = config[val_iZone]->GetPISO_corrections();

  /*--- Setting up iteration values depending on if this is a
   steady or an unsteady simulation */

  const auto InnerIter = config[val_iZone]->GetInnerIter();
  const auto TimeIter = config[val_iZone]->GetTimeIter();

  /*--- Solve the Euler, Navier-Stokes, RANS equations. ---*/

  const auto main_solver = config[val_iZone]->GetKind_Solver();
  config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_FLOW_SYS);
  
  /*--- Solve the momentum equations (prediction). ---*/

  integration[val_iZone][val_iInst][FLOW_SOL]->MultiGrid_Iteration(geometry, solver, numerics, config, RUNTIME_FLOW_SYS,
                                                                   val_iZone, val_iInst);

  /*--- Solve the pressure correction (poisson) equation ---*/
  /* The poisson correction is repeated nCorrections number of times to account for the PISO algorithm.
  Note that the pressure- and velocity corrections resulting from the poisson equation 
  are handled in the postprocessing step of the poissonsolver. */

  config[val_iZone]->SetGlobalParam(MAIN_SOLVER::POISSON_EQUATION, RUNTIME_POISSON_SYS);

  for (int i = 0; i < nCorrections; ++i) integration[val_iZone][val_iInst][POISSON_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config, RUNTIME_POISSON_SYS,
                                                                   val_iZone, val_iInst);

  /*--- Pressure-based algorithm finished, now run auxiliary solvers ---*/
  // TODO: this list is the same as the list inside of CFluidIteration and should be made into a separate function 

  /*--- If the flow integration is not fully coupled, run the various single grid integrations. ---*/

  if (config[val_iZone]->GetKind_Turb_Model() != TURB_MODEL::NONE && !frozen_visc) {

    /*--- Solve transition model ---*/

    if (config[val_iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_TRANS_SYS);
      integration[val_iZone][val_iInst][TRANS_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                         RUNTIME_TRANS_SYS, val_iZone, val_iInst);
    }

    /*--- Solve the turbulence model ---*/

    config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_TURB_SYS);
    integration[val_iZone][val_iInst][TURB_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                      RUNTIME_TURB_SYS, val_iZone, val_iInst);
  }

  if (config[val_iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
    config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_SPECIES_SYS);
    integration[val_iZone][val_iInst][SPECIES_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                         RUNTIME_SPECIES_SYS, val_iZone, val_iInst);

    // This only applies if mixture properties are used. But this also doesn't hurt if done w/out mixture properties.
    // In case of turbulence, the Turb-Post computes the correct eddy viscosity based on mixture-density and
    // mixture lam-visc. In order to get the correct mixture properties, based on the just updated mass-fractions, the
    // Flow-Pre has to be called upfront. The updated eddy-visc are copied into the flow-solver Primitive in another
    // Flow-Pre call which is done at the start of the next iteration.
    if (config[val_iZone]->GetKind_Turb_Model() != TURB_MODEL::NONE) {
      solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
      solver[val_iZone][val_iInst][MESH_0][TURB_SOL]->Postprocessing(geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], MESH_0);
    }
  }

  //TODO: this currently does not work, likely because the temperature inside of CPBIncEulerVariable is NOT set like it should.
  if (config[val_iZone]->GetWeakly_Coupled_Heat()) {
    config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_HEAT_SYS);
    integration[val_iZone][val_iInst][HEAT_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                      RUNTIME_HEAT_SYS, val_iZone, val_iInst);
  }

  /*--- Incorporate a weakly-coupled radiation model to the analysis ---*/
  if (config[val_iZone]->AddRadiation()) {
    config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_RADIATION_SYS);
    integration[val_iZone][val_iInst][RAD_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                     RUNTIME_RADIATION_SYS, val_iZone, val_iInst);
  }

  /*--- Adapt the CFL number using an exponential progression with under-relaxation approach. ---*/

  if ((config[val_iZone]->GetCFL_Adapt() == YES) && (!disc_adj)) {
    SU2_OMP_PARALLEL
    solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->AdaptCFLNumber(geometry[val_iZone][val_iInst],
                                                                   solver[val_iZone][val_iInst], config[val_iZone]);
    END_SU2_OMP_PARALLEL
  }

}
