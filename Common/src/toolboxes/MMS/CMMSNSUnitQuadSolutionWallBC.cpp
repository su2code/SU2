/*!
 * \file CMMSNSUnitQuadSolutionWallBC.cpp
 * \brief Implementations of the member functions of CMMSNSUnitQuadSolutionWallBC.
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

#include "../../../include/toolboxes/MMS/CMMSNSUnitQuadSolutionWallBC.hpp"

CMMSNSUnitQuadSolutionWallBC::CMMSNSUnitQuadSolutionWallBC() : CVerificationSolution() {}

CMMSNSUnitQuadSolutionWallBC::CMMSNSUnitQuadSolutionWallBC(unsigned short val_nDim, unsigned short val_nVar,
                                                           unsigned short val_iMesh, CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  /*--- Write a message that the solution is initialized for the manufactured
   solution for the Navier-Stokes equations on a unit quad with no-slip
   wall boundary conditions. ---*/
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the manufactured solution " << endl;
    cout << "         of the Navier-Stokes equations on a unit quad" << endl;
    cout << "         with no-slip wall boundary conditions!!!" << endl;
    cout << endl << flush;
  }

  /*--- Coefficients, needed to determine the solution. ---*/
  const su2double Prandtl = config->GetPrandtl_Lam();

  RGas = config->GetGas_Constant();
  Gamma = config->GetGamma();
  Viscosity = config->GetMu_Constant();
  Conductivity = Viscosity * Gamma * RGas / (Prandtl * (Gamma - 1.0));

  /*--- Initialize TWall to the default value of 300 K (in case the outer wall
        is not modelled as an isothermal wall) and try to retrieve the wall
        temperature from the boundary conditions. ---*/
  TWall = 300.0;
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) {
      const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
      TWall = config->GetIsothermal_Temperature(Marker_Tag);
    }
  }

  /*--- Get the reference values for pressure, density and velocity. ---*/
  Pressure_Ref = config->GetPressure_Ref();
  Density_Ref = config->GetDensity_Ref();
  Velocity_Ref = config->GetVelocity_Ref();

  /*--- The constants for the density and velocities. ---*/
  rho_0 = 1.25;
  u_0 = 135.78;
  v_0 = -67.61;

  /*--- The constants for the temperature solution. ---*/
  a_T1 = 1.05;
  a_T2 = -0.85;

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if (config->GetTime_Marching() != TIME_MARCHING::STEADY)
    SU2_MPI::Error("Steady mode must be selected for the MMS NS Unit Quad case with wall BCs.", CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::EULER && Kind_Solver != MAIN_SOLVER::NAVIER_STOKES &&
      Kind_Solver != MAIN_SOLVER::RANS && Kind_Solver != MAIN_SOLVER::FEM_EULER &&
      Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::FEM_RANS &&
      Kind_Solver != MAIN_SOLVER::FEM_LES)
    SU2_MPI::Error("Compressible flow equations must be selected for the MMS NS Unit Quad case with wall BCs.",
                   CURRENT_FUNCTION);

  if ((Kind_Solver != MAIN_SOLVER::NAVIER_STOKES) && (Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES))
    SU2_MPI::Error("Navier Stokes equations must be selected for the MMS NS Unit Quad case with wall BCs.",
                   CURRENT_FUNCTION);

  if ((config->GetKind_FluidModel() != STANDARD_AIR) && (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the MMS NS Unit Quad case with wall BCs.",
                   CURRENT_FUNCTION);

  if (config->GetKind_ViscosityModel() != VISCOSITYMODEL::CONSTANT)
    SU2_MPI::Error("Sutherland must be selected for viscosity for the MMS NS Unit Quad case with wall BCs.",
                   CURRENT_FUNCTION);

  if (config->GetKind_ConductivityModel() != CONDUCTIVITYMODEL::CONSTANT_PRANDTL)
    SU2_MPI::Error("Constant Prandtl number must be selected for the MMS NS Unit Quad case with wall BCs.",
                   CURRENT_FUNCTION);
}

CMMSNSUnitQuadSolutionWallBC::~CMMSNSUnitQuadSolutionWallBC() = default;

void CMMSNSUnitQuadSolutionWallBC::GetBCState(const su2double* val_coords, const su2double val_t,
                                              su2double* val_solution) const {
  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CMMSNSUnitQuadSolutionWallBC::GetSolution(const su2double* val_coords, const su2double val_t,
                                               su2double* val_solution) const {
  /* Easier storage of the y-coordinate. */
  const su2double y = val_coords[1];

  /* Determine the dimensional solution for the temperature. */
  const su2double Pi = PI_NUMBER;
  const su2double fact = y * y / (a_T1 + a_T2);

  su2double T = 0.25 * TWall * (3.0 + fact * (a_T1 * cos(Pi * (y - 1.0)) + a_T2 * cos(Pi * (y - 1.0) * 2.0)));

  /* Determine the dimensional solution for the velocities. */
  su2double u = u_0 * y * (1.0 - y) * 4.0;
  su2double v = v_0 * y * (1.0 - y) * 4.0;

  /* Compute the pressure from the density and temperature. */
  su2double rho = rho_0;
  su2double p = rho * RGas * T;

  /* Determine the non-dimensional solution. */
  rho /= Density_Ref;
  p /= Pressure_Ref;
  u /= Velocity_Ref;
  v /= Velocity_Ref;

  /* Compute the conservative variables from the primitive ones.
     Note that the implementation below is valid for both 2D and 3D. */
  val_solution[0] = rho;
  val_solution[1] = rho * u;
  val_solution[2] = rho * v;
  val_solution[3] = 0.0;
  val_solution[nDim + 1] = p / (Gamma - 1.0) + 0.5 * rho * (u * u + v * v);
}

void CMMSNSUnitQuadSolutionWallBC::GetMMSSourceTerm(const su2double* val_coords, const su2double val_t,
                                                    su2double* val_source) const {
  /*--- Abbreviate Pi and the y-coordinate. ---*/
  const su2double Pi = PI_NUMBER;
  const su2double y = val_coords[1];

  /*--- The source code for the source terms is generated in Maple.
        See the file CMMSNSUnitQuadSolutionWallBC.mw in the directory
        CreateMMSSourceTerms for the details how to do this. ---*/
  const su2double t1 = (v_0 * rho_0);
  const su2double t2 = 1.0 - y;
  const su2double t6 = rho_0 * u_0;
  const su2double t7 = t2 * t2;
  const su2double t8 = t7 * y;
  const su2double t12 = y * y;
  const su2double t13 = t2 * t12;
  const su2double t20 = (v_0 * v_0);
  const su2double t21 = (t20 * rho_0);
  const su2double t26 = rho_0 * RGas;
  const su2double t27 = a_T1 + a_T2;
  const su2double t28 = 1.0 / t27;
  const su2double t30 = -t2 * Pi;
  const su2double t31 = cos(t30);
  const su2double t32 = t31 * a_T1;
  const su2double t33 = 2.0 * t30;
  const su2double t34 = cos(t33);
  const su2double t36 = t34 * a_T2 + t32;
  const su2double t39 = t28 * t12;
  const su2double t41 = sin(t30);
  const su2double t42 = t41 * Pi * a_T1;
  const su2double t43 = Pi * a_T2;
  const su2double t44 = sin(t33);
  const su2double t50 = (0.2e1 * t36 * t28 * y + (-0.2e1 * t44 * t43 - t42) * t39) * TWall;
  const su2double t52 = t50 * t26 / 0.4e1;
  const su2double t57 = 1.0 / (Gamma - 1.0);
  const su2double t61 = (u_0 * u_0);
  const su2double t62 = y * t61;
  const su2double t64 = t12 * t61;
  const su2double t66 = y * t20;
  const su2double t68 = t12 * t20;
  const su2double t75 = t2 * y;
  const su2double t80 = (t36 * t39 + 0.3e1) * TWall;
  const su2double t92 = v_0 * (t57 * t80 * t26 / 0.4e1 + (8 * (t7 * t64 + t7 * t68) * rho_0) + t80 * t26 / 0.4e1);
  const su2double t97 = t41 * t31;
  const su2double t102 = a_T1 * y;
  const su2double t104 = t31 * t31;
  const su2double t111 = 0.1e1 / t27 / 0.4e1;
  const su2double t115 = Pi * Pi;
  const su2double t116 = t41 * t41;
  const su2double t118 = a_T2 * y;
  const su2double t138 = Viscosity * (-8.0 * y + 4.0);
  const su2double t149 = Viscosity * (2.0 * y - 1.0);
  const su2double t155 =
      0.4e1 * t75 * v_0 *
          (t57 * t50 * t26 / 0.4e1 + (16.0 * (-t2 * t64 - t2 * t68 + t7 * t62 + t7 * t66) * rho_0) + t52) +
      0.4e1 * t2 * t92 - 0.4e1 * y * t92 -
      Conductivity * t111 *
          (-t102 * Pi * t41 - 0.4e1 * y * t43 * t97 + 0.4e1 * t104 * a_T2 - 0.2e1 * a_T2 + 0.2e1 * t32) * TWall -
      Conductivity * t111 *
          (-t102 * t31 * t115 - 0.4e1 * t118 * t115 * t104 + 0.4e1 * t118 * t116 * t115 - 0.12e2 * t43 * t97 -
           0.3e1 * t42) *
          TWall * y +
      (32.0 * t75 * t61 * Viscosity) - (4.0 * t2 * t61 * t138) + (4.0 * t62 * t138) +
      0.128e3 / 0.3e1 * t75 * t20 * Viscosity + 0.64e2 / 0.3e1 * t2 * t20 * t149 - 0.64e2 / 0.3e1 * t66 * t149;

  /*--- Set the source term, which is valid for both 2D and 3D cases.
        Note the scaling for the correct non-dimensionalization. ---*/
  val_source[0] = 4.0 * t2 * t1 - 4.0 * y * t1;
  val_source[1] = (-0.32e2 * v_0 * t13 * t6 + 0.32e2 * v_0 * t8 * t6 + (8.0 * Viscosity * u_0));
  val_source[2] = ((32.0 * t8 * t21) - (32.0 * t13 * t21) + t52 + 0.32e2 / 0.3e1 * Viscosity * v_0);
  val_source[3] = 0.0;
  val_source[nDim + 1] = t155;

  val_source[0] /= Density_Ref * Velocity_Ref;
  val_source[1] /= Pressure_Ref;
  val_source[2] /= Pressure_Ref;
  val_source[nDim + 1] /= Velocity_Ref * Pressure_Ref;
}

bool CMMSNSUnitQuadSolutionWallBC::IsManufacturedSolution() const { return true; }
