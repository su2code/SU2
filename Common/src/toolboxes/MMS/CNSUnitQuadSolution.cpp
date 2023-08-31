/*!
 * \file CNSUnitQuadSolution.cpp
 * \brief Implementations of the member functions of CNSUnitQuadSolution.
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

#include "../../../include/toolboxes/MMS/CNSUnitQuadSolution.hpp"

CNSUnitQuadSolution::CNSUnitQuadSolution() : CVerificationSolution() {}

CNSUnitQuadSolution::CNSUnitQuadSolution(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_iMesh,
                                         CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  /*--- Write a message that the solution is initialized for the
   Navier-Stokes case on a unit quad. Note that heat conduction
   is neglected for this case. ---*/

  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the Navier Stokes equations " << endl;
    cout << "         on a unit quad. Note that heat conduction " << endl;
    cout << "         is neglected for this case!!!" << endl;
    cout << endl << flush;
  }

  /*--- Coefficients, needed to determine the solution. ---*/
  Gm1 = config->GetGamma() - 1.0;
  flowAngle = config->GetAoA() * PI_NUMBER / 180.0;
  Viscosity = config->GetMu_ConstantND();

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if (config->GetTime_Marching() != TIME_MARCHING::STEADY)
    SU2_MPI::Error("Steady mode must be selected for the NS Unit Quad case", CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::EULER && Kind_Solver != MAIN_SOLVER::NAVIER_STOKES &&
      Kind_Solver != MAIN_SOLVER::RANS && Kind_Solver != MAIN_SOLVER::FEM_EULER &&
      Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::FEM_RANS &&
      Kind_Solver != MAIN_SOLVER::FEM_LES)
    SU2_MPI::Error("Compressible flow equations must be selected for the NS Unit Quad case", CURRENT_FUNCTION);

  if ((Kind_Solver != MAIN_SOLVER::NAVIER_STOKES) && (Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES))
    SU2_MPI::Error("Navier Stokes equations must be selected for the NS Unit Quad case", CURRENT_FUNCTION);

  if (config->GetKind_FluidModel() != IDEAL_GAS)
    SU2_MPI::Error("Ideal gas must be selected for the NS Unit Quad case", CURRENT_FUNCTION);

  if (fabs(Gm1 - 0.5) > 1.e-8) SU2_MPI::Error("Gamma must be 1.5 for the NS Unit Quad case", CURRENT_FUNCTION);

  if (config->GetKind_ViscosityModel() != VISCOSITYMODEL::CONSTANT)
    SU2_MPI::Error("Constant viscosity must be selected for the NS Unit Quad case", CURRENT_FUNCTION);

  if (config->GetKind_ConductivityModel() != CONDUCTIVITYMODEL::CONSTANT_PRANDTL)
    SU2_MPI::Error("Constant Prandtl number must be selected for the NS Unit Quad case", CURRENT_FUNCTION);

  if (config->GetPrandtl_Lam() < 1.e+20)
    SU2_MPI::Error("Laminar Prandtl number must be larger than 1.e+20 for the NS Unit Quad case", CURRENT_FUNCTION);
}

CNSUnitQuadSolution::~CNSUnitQuadSolution() = default;

void CNSUnitQuadSolution::GetBCState(const su2double* val_coords, const su2double val_t,
                                     su2double* val_solution) const {
  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CNSUnitQuadSolution::GetSolution(const su2double* val_coords, const su2double val_t,
                                      su2double* val_solution) const {
  /*--- Compute the flow direction and the coordinates in
        the rotated frame. ---*/
  const su2double cosFlowAngle = cos(flowAngle);
  const su2double sinFlowAngle = sin(flowAngle);

  const su2double xTilde = val_coords[0] * cosFlowAngle - val_coords[1] * sinFlowAngle;
  const su2double yTilde = val_coords[0] * sinFlowAngle + val_coords[1] * cosFlowAngle;

  /*--- Compute the exact solution for this case. Note that it works
        both in 2D and 3D. ---*/
  val_solution[0] = 1.0;
  val_solution[1] = cosFlowAngle * yTilde * yTilde;
  val_solution[2] = -sinFlowAngle * yTilde * yTilde;
  val_solution[3] = 0.0;
  val_solution[nVar - 1] = (2.0 * Viscosity * xTilde + 10.0) / Gm1 + 0.5 * yTilde * yTilde * yTilde * yTilde;
}
