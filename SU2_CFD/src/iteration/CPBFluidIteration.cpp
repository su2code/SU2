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
#include "../../include/integration/CIntegration.hpp"

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
  const bool periodic = (config[val_iZone]->GetnMarker_Periodic() > 0);

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

  /*--- Solve the pressure poisson (correction) equation ---*/
  /* The poisson correction is repeated nCorrections number of times to account for the PISO algorithm.
  Note that the pressure- and velocity corrections resulting from the poisson equation 
  are handled in the postprocessing step of the poissonsolver. */

  config[val_iZone]->SetGlobalParam(MAIN_SOLVER::POISSON_EQUATION, RUNTIME_POISSON_SYS);

  /*--- The momentum coefficients (resulting from the flow solution) are set only once at the start of the corrections ---*/
  solver[val_iZone][val_iInst][MESH_0][POISSON_SOL]->SetMomCoeff(geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], periodic, MESH_0);

  /* TODO: The current piso style iteration works although it has no effect on the stable time step sizes as expected
   * it is thus likely not yet optimal and should be looked at. Additionally, the Poisson solve is very difficult
   * to solve at large Reynolds number cases and would therefore massively benefit from Multigrid. The current 
   * limitation here is that SU2 as of today does not allow multigrid for different solvers and has multigrid 
   * only as an option for the flow solver. */
  for (unsigned short i = 0; i < nCorrections; ++i) {
    integration[val_iZone][val_iInst][POISSON_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config, RUNTIME_POISSON_SYS,
                                                                   val_iZone, val_iInst);
  }

  /*--- Pressure-based algorithm finished, now run auxiliary solvers ---*/

  /*--- If the flow integration is not fully coupled, run the various single grid integrations. ---*/
  CommonSingleGridIterations(output, integration, geometry, solver, numerics, config, surface_movement, 
                             grid_movement, FFDBox, val_iZone, val_iInst, main_solver, frozen_visc);

  /*--- Adapt the CFL number using an exponential progression with under-relaxation approach. ---*/

  if ((config[val_iZone]->GetCFL_Adapt() == YES) && (!disc_adj)) {
    SU2_OMP_PARALLEL
    solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->AdaptCFLNumber(geometry[val_iZone][val_iInst],
                                                                   solver[val_iZone][val_iInst], config[val_iZone]);
    END_SU2_OMP_PARALLEL
  }

}
