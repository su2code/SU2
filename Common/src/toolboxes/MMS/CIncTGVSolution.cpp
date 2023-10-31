/*!
 * \file CIncTGVSolution.cpp
 * \brief Implementations of the member functions of CIncTGVSolution.
 * \author T. Economon, E. van der Weide
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

#include "../../../include/toolboxes/MMS/CIncTGVSolution.hpp"

CIncTGVSolution::CIncTGVSolution() : CVerificationSolution() {}

CIncTGVSolution::CIncTGVSolution(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_iMesh,
                                 CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  /*--- Disable this solution for now, as it has not been tested. ---*/

  SU2_MPI::Error("CIncTGVSolution not yet fully implemented/tested.", CURRENT_FUNCTION);

  /*--- Write a message that the solution is initialized for the
   Taylor-Green vortex test case. ---*/

  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the incompressible Taylor-Green vortex case!!!" << endl;
    cout << endl << flush;
  }

  /*--- Store TGV specific parameters here. ---*/

  tgvLength = 1.0;
  tgvVelocity = 1.0;
  tgvDensity = config->GetDensity_FreeStreamND();
  tgvViscosity = config->GetViscosity_FreeStreamND();

  /*--- We keep a copy of the freestream temperature just to be safe
   when we set the solution, even though this is an isothermal case. ---*/

  Temperature = config->GetTemperature_FreeStreamND();

  /*--- Perform some sanity and error checks for this solution here. ---*/

  if ((config->GetTime_Marching() != TIME_MARCHING::TIME_STEPPING) &&
      (config->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_1ST) &&
      (config->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_2ND))
    SU2_MPI::Error("Unsteady mode must be selected for the incompressible Taylor Green Vortex", CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::INC_EULER && Kind_Solver != MAIN_SOLVER::INC_NAVIER_STOKES &&
      Kind_Solver != MAIN_SOLVER::INC_RANS)
    SU2_MPI::Error("Incompressible flow equations must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::INC_NAVIER_STOKES)
    SU2_MPI::Error("Navier Stokes equations must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);

  if (config->GetKind_FluidModel() != CONSTANT_DENSITY)
    SU2_MPI::Error("Constant density fluid model must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);

  if (config->GetKind_ViscosityModel() != VISCOSITYMODEL::CONSTANT)
    SU2_MPI::Error("Constant viscosity must be selected for the incompressible Taylor Green Vortex", CURRENT_FUNCTION);

  if (config->GetEnergy_Equation())
    SU2_MPI::Error("Energy equation must be disabled (isothermal) for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);

  if (nDim != 2) SU2_MPI::Error("2D calculation required for the incompressible Taylor Green Vortex", CURRENT_FUNCTION);
}

CIncTGVSolution::~CIncTGVSolution() = default;

void CIncTGVSolution::GetBCState(const su2double* val_coords, const su2double val_t, su2double* val_solution) const {
  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CIncTGVSolution::GetSolution(const su2double* val_coords, const su2double val_t, su2double* val_solution) const {
  /* The exact solution is set for the incompressible Taylor-Green
   vortex case. This is the classic solution from the original work
   of Taylor and Green for the specific 2D situation where the
   exact solution can be derived for an incompressible flow. */

  /* Store the termporal term more easily (Taylor expansion). */
  su2double F = 1.0 - 2.0 * (tgvViscosity / tgvDensity) * val_t;

  /* Compute the primitive variables. */
  su2double u = tgvVelocity * F * (sin(val_coords[0] / tgvLength) * cos(val_coords[1] / tgvLength));
  su2double v = -tgvVelocity * F * (cos(val_coords[0] / tgvLength) * sin(val_coords[1] / tgvLength));

  su2double B = (cos(2.0 * val_coords[0] / tgvLength) + cos(2.0 * val_coords[1] / tgvLength));
  su2double p = -(tgvDensity / 4.0) * B * F * F;

  /* Compute the conservative variables. Note that both 2D and 3D
   cases are treated correctly. */
  val_solution[0] = p;
  val_solution[1] = u;
  val_solution[2] = v;
  val_solution[nVar - 1] = Temperature;
}
