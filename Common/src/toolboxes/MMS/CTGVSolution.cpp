/*!
 * \file CTGVSolution.cpp
 * \brief Implementations of the member functions of CTGVSolution.
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

#include "../../../include/toolboxes/MMS/CTGVSolution.hpp"

CTGVSolution::CTGVSolution() : CVerificationSolution() {}

CTGVSolution::CTGVSolution(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_iMesh, CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  /*--- Write a message that the solution is initialized for the
   Taylor-Green vortex test case. ---*/

  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the Taylor-Green vortex case!!!" << endl;
    cout << endl << flush;
  }

  /*--- Store TGV specific parameters here. ---*/

  tgvLength = 1.0;      // Taylor-Green length scale.
  tgvVelocity = 1.0;    // Taylor-Green velocity.
  tgvDensity = 1.0;     // Taylor-Green density.
  tgvPressure = 100.0;  // Taylor-Green pressure.

  /*--- Useful coefficient in which Gamma is present. ---*/

  ovGm1 = 1.0 / (config->GetGamma() - 1.0);

  /*--- Perform some sanity and error checks for this solution here. ---*/

  if ((config->GetTime_Marching() != TIME_MARCHING::TIME_STEPPING) &&
      (config->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_1ST) &&
      (config->GetTime_Marching() != TIME_MARCHING::DT_STEPPING_2ND))
    SU2_MPI::Error("Unsteady mode must be selected for the Taylor Green Vortex", CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::EULER && Kind_Solver != MAIN_SOLVER::NAVIER_STOKES &&
      Kind_Solver != MAIN_SOLVER::RANS && Kind_Solver != MAIN_SOLVER::FEM_EULER &&
      Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::FEM_RANS &&
      Kind_Solver != MAIN_SOLVER::FEM_LES)
    SU2_MPI::Error("Compressible flow equations must be selected for the Taylor Green Vortex", CURRENT_FUNCTION);

  if ((Kind_Solver != MAIN_SOLVER::NAVIER_STOKES) && (Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES))
    SU2_MPI::Error("Navier Stokes equations must be selected for the Taylor Green Vortex", CURRENT_FUNCTION);

  if ((config->GetKind_FluidModel() != STANDARD_AIR) && (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the Taylor Green Vortex", CURRENT_FUNCTION);

  if (config->GetKind_ViscosityModel() != VISCOSITYMODEL::CONSTANT)
    SU2_MPI::Error("Constant viscosity must be selected for the Taylor Green Vortex", CURRENT_FUNCTION);

  if (config->GetKind_ConductivityModel() != CONDUCTIVITYMODEL::CONSTANT_PRANDTL)
    SU2_MPI::Error("Constant Prandtl number must be selected for the Taylor Green Vortex", CURRENT_FUNCTION);
}

CTGVSolution::~CTGVSolution() = default;

void CTGVSolution::GetSolution(const su2double* val_coords, const su2double val_t, su2double* val_solution) const {
  /* The initial conditions are set for the Taylor-Green vortex case, which
   is a DNS case that features vortex breakdown into turbulence. These
   particular settings are for the typical Re = 1600 case (M = 0.08) with
   an initial temperature of 300 K. Note that this condition works in both
   2D and 3D. */

  su2double val_coordsZ = 0.0;
  if (nDim == 3) val_coordsZ = val_coords[2];

  /* Compute the primitive variables. */

  su2double rho = tgvDensity;
  su2double u =
      tgvVelocity * (sin(val_coords[0] / tgvLength) * cos(val_coords[1] / tgvLength) * cos(val_coordsZ / tgvLength));
  su2double v =
      -tgvVelocity * (cos(val_coords[0] / tgvLength) * sin(val_coords[1] / tgvLength) * cos(val_coordsZ / tgvLength));

  su2double factorA = cos(2.0 * val_coordsZ / tgvLength) + 2.0;
  su2double factorB = (cos(2.0 * val_coords[0] / tgvLength) + cos(2.0 * val_coords[1] / tgvLength));

  su2double p = (tgvPressure + tgvDensity * (pow(tgvVelocity, 2.0) / 16.0) * factorA * factorB);

  /* Compute the conservative variables. Note that both 2D and 3D
   cases are treated correctly. */

  val_solution[0] = rho;
  val_solution[1] = rho * u;
  val_solution[2] = rho * v;
  val_solution[3] = 0.0;
  val_solution[nVar - 1] = p * ovGm1 + 0.5 * rho * (u * u + v * v);
}

bool CTGVSolution::ExactSolutionKnown() const { return false; }
