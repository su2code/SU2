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

  /*--- Setting up iteration values depending on if this is a
   steady or an unsteady simulation */

  const auto InnerIter = config[val_iZone]->GetInnerIter();
  const auto TimeIter = config[val_iZone]->GetTimeIter();

  /*--- Solve the Euler, Navier-Stokes, RANS equations. ---*/
  
  /*--- Solve the momentum equations. ---*/

  /*--- Solve the pressure correction (poisson) equation ---*/

  /*--- Correct pressure and velocities. ---*/

  /*--- T.B.D. ---*/

}
