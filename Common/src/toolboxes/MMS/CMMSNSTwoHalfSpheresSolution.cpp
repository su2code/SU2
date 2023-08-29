/*!
 * \file CMMSNSTwoHalfSpheresSolution.cpp
 * \brief Implementations of the member functions of CMMSNSTwoHalfSpheresSolution.
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

#include "../../../include/toolboxes/MMS/CMMSNSTwoHalfSpheresSolution.hpp"

CMMSNSTwoHalfSpheresSolution::CMMSNSTwoHalfSpheresSolution() : CVerificationSolution() {}

CMMSNSTwoHalfSpheresSolution::CMMSNSTwoHalfSpheresSolution(unsigned short val_nDim, unsigned short val_nVar,
                                                           unsigned short val_iMesh, CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  /*--- Write a message that the solution is initialized for the manufactured
   solution for the Navier-Stokes equations between two half spheres. ---*/
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the manufactured solution " << endl;
    cout << "         of the Navier-Stokes equations between two half spheres!!!" << endl;
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
  w_0 = 82.75;

  /*--- The constants for the temperature solution. ---*/
  a_T1 = 1.05;
  a_T2 = -0.85;

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if (nDim != 3) SU2_MPI::Error("Grid must be 3D for the MMS NS Two Half Spheres case", CURRENT_FUNCTION);

  if (config->GetTime_Marching() != TIME_MARCHING::STEADY)
    SU2_MPI::Error("Steady mode must be selected for the MMS NS Two Half Spheres case", CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::EULER && Kind_Solver != MAIN_SOLVER::NAVIER_STOKES &&
      Kind_Solver != MAIN_SOLVER::RANS && Kind_Solver != MAIN_SOLVER::FEM_EULER &&
      Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::FEM_RANS &&
      Kind_Solver != MAIN_SOLVER::FEM_LES)
    SU2_MPI::Error("Compressible flow equations must be selected for the MMS NS Two Half Spheres case",
                   CURRENT_FUNCTION);

  if ((Kind_Solver != MAIN_SOLVER::NAVIER_STOKES) && (Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES))
    SU2_MPI::Error("Navier Stokes equations must be selected for the MMS NS Two Half Spheres case", CURRENT_FUNCTION);

  if ((config->GetKind_FluidModel() != STANDARD_AIR) && (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the MMS NS Two Half Spheres case", CURRENT_FUNCTION);

  if (config->GetKind_ViscosityModel() != VISCOSITYMODEL::CONSTANT)
    SU2_MPI::Error("Sutherland must be selected for viscosity for the MMS NS Two Half Spheres case", CURRENT_FUNCTION);

  if (config->GetKind_ConductivityModel() != CONDUCTIVITYMODEL::CONSTANT_PRANDTL)
    SU2_MPI::Error("Constant Prandtl number must be selected for the MMS NS Two Half Spheres case", CURRENT_FUNCTION);
}

CMMSNSTwoHalfSpheresSolution::~CMMSNSTwoHalfSpheresSolution() = default;

void CMMSNSTwoHalfSpheresSolution::GetBCState(const su2double* val_coords, const su2double val_t,
                                              su2double* val_solution) const {
  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CMMSNSTwoHalfSpheresSolution::GetSolution(const su2double* val_coords, const su2double val_t,
                                               su2double* val_solution) const {
  /* Easier storage of the x-, y- and z-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  const su2double z = val_coords[2];

  /* Determine the dimensional solution for the temperature. */
  const su2double Pi = PI_NUMBER;
  const su2double r = sqrt(x * x + y * y + z * z);
  const su2double fact = (r - 1.0) * (r - 1.0) / (a_T1 + a_T2);

  su2double T = 0.25 * TWall * (3.0 + fact * (a_T1 * cos(Pi * (r - 2.0)) + a_T2 * cos(Pi * (r - 2.0) * 2.0)));

  /* Determine the dimensional solution for the velocities. */
  su2double u = u_0 * (r - 1.0) * (2.0 - r) * 4.0;
  su2double v = v_0 * (r - 1.0) * (2.0 - r) * 4.0;
  su2double w = w_0 * (r - 1.0) * (2.0 - r) * 4.0;

  /* Compute the pressure from the density and temperature. */
  su2double rho = rho_0;
  su2double p = rho * RGas * T;

  /* Determine the non-dimensional solution. */
  rho /= Density_Ref;
  p /= Pressure_Ref;
  u /= Velocity_Ref;
  v /= Velocity_Ref;
  w /= Velocity_Ref;

  /* Determine the non-dimensional conserved variables. */
  val_solution[0] = rho;
  val_solution[1] = rho * u;
  val_solution[2] = rho * v;
  val_solution[3] = rho * w;
  val_solution[4] = p / (Gamma - 1.0) + 0.5 * rho * (u * u + v * v + w * w);
}

void CMMSNSTwoHalfSpheresSolution::GetMMSSourceTerm(const su2double* val_coords, const su2double val_t,
                                                    su2double* val_source) const {
  /*--- Abbreviate Pi and the coordinates. ---*/
  const su2double Pi = PI_NUMBER;
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  const su2double z = val_coords[2];

  /*--- The source code for the source terms is generated in Maple.
        See the file CMMSNSTwoHalfSpheresSolution.mw in the directory
        CreateMMSSourceTerms for the details how to do this. ---*/
  const su2double t1 = rho_0 * u_0;
  const su2double t2 = x * x;
  const su2double t3 = y * y;
  const su2double t4 = z * z;
  const su2double t5 = t2 + t3 + t4;
  const su2double t6 = sqrt(t5);
  const su2double t7 = 0.1e1 / t6;
  const su2double t8 = x * t7;
  const su2double t9 = 0.2e1 - t6;
  const su2double t10 = t9 * t8;
  const su2double t12 = t6 - 0.1e1;
  const su2double t13 = t7 * t12;
  const su2double t14 = x * t13;
  const su2double t16 = rho_0 * v_0;
  const su2double t17 = y * t7;
  const su2double t18 = t9 * t17;
  const su2double t20 = y * t13;
  const su2double t22 = rho_0 * w_0;
  const su2double t23 = z * t7;
  const su2double t24 = t9 * t23;
  const su2double t26 = z * t13;
  const su2double t29 = u_0 * u_0;
  const su2double t30 = t29 * rho_0;
  const su2double t32 = t9 * t9;
  const su2double t33 = t7 * t32;
  const su2double t34 = x * t33;
  const su2double t37 = t12 * t12;
  const su2double t41 = rho_0 * RGas;
  const su2double t42 = a_T1 + a_T2;
  const su2double t43 = 0.1e1 / t42;
  const su2double t44 = t43 * t12;
  const su2double t45 = -t9 * Pi;
  const su2double t46 = cos(t45);
  const su2double t47 = t46 * a_T1;
  const su2double t48 = 0.2e1 * t45;
  const su2double t49 = cos(t48);
  const su2double t51 = t49 * a_T2 + t47;
  const su2double t52 = t7 * t51;
  const su2double t56 = t43 * t37;
  const su2double t57 = a_T1 * Pi;
  const su2double t58 = sin(t45);
  const su2double t60 = t58 * t8 * t57;
  const su2double t61 = a_T2 * Pi;
  const su2double t62 = sin(t48);
  const su2double t70 = (x * t52 * t44 / 0.2e1 + (-0.2e1 * t62 * t8 * t61 - t60) * t56 / 0.4e1) * TWall;
  const su2double t71 = t70 * t41;
  const su2double t72 = 0.1e1 / t5;
  const su2double t73 = t72 * Viscosity;
  const su2double t74 = u_0 * x;
  const su2double t76 = v_0 * y;
  const su2double t77 = w_0 * z;
  const su2double t78 = 0.2e1 * t74 - t76 - t77;
  const su2double t84 = (-0.3e1 + 0.2e1 * t6) * Viscosity;
  const su2double t85 = t7 * u_0;
  const su2double t89 = 0.1e1 / t6 / t5;
  const su2double t90 = t89 * t78;
  const su2double t94 = t12 * t1;
  const su2double t95 = v_0 * t32;
  const su2double t99 = t37 * t1;
  const su2double t100 = v_0 * t9;
  const su2double t104 = y * t72;
  const su2double t107 = u_0 * y + v_0 * x;
  const su2double t108 = t107 * Viscosity;
  const su2double t111 = -0.3e1 / 0.2e1 + t6;
  const su2double t112 = Viscosity * t111;
  const su2double t115 = t89 * t107;
  const su2double t119 = w_0 * t32;
  const su2double t120 = t23 * t119;
  const su2double t123 = w_0 * t9;
  const su2double t124 = t23 * t123;
  const su2double t129 = u_0 * z + w_0 * x;
  const su2double t130 = t72 * t129;
  const su2double t131 = z * Viscosity;
  const su2double t134 = t111 * t129;
  const su2double t135 = t89 * Viscosity;
  const su2double t136 = z * t135;
  const su2double t139 = 0.32e2 * t34 * t12 * t30 - 0.32e2 * t10 * t37 * t30 + t71 + 0.16e2 / 0.3e1 * t78 * x * t73 +
                         0.16e2 / 0.3e1 * t85 * t84 - 0.8e1 / 0.3e1 * x * t90 * t84 + 0.32e2 * t17 * t95 * t94 -
                         0.32e2 * t17 * t100 * t99 + 0.8e1 * t108 * t104 + 0.16e2 * t85 * t112 -
                         0.8e1 * y * t115 * t112 + 0.32e2 * t120 * t94 - 0.32e2 * t124 * t99 + 0.8e1 * t131 * t130 -
                         0.8e1 * t136 * t134;
  const su2double t146 = x * t72;
  const su2double t155 = v_0 * v_0;
  const su2double t156 = t155 * rho_0;
  const su2double t158 = y * t33;
  const su2double t168 = t58 * t17 * t57;
  const su2double t176 = (y * t52 * t44 / 0.2e1 + (-0.2e1 * t62 * t17 * t61 - t168) * t56 / 0.4e1) * TWall;
  const su2double t177 = t176 * t41;
  const su2double t179 = t74 - 0.2e1 * t76 + t77;
  const su2double t180 = t111 * t179;
  const su2double t187 = t12 * t16;
  const su2double t190 = t37 * t16;
  const su2double t195 = v_0 * z + w_0 * y;
  const su2double t196 = t72 * t195;
  const su2double t199 = t111 * t195;
  const su2double t202 = 0.32e2 * t8 * t95 * t94 - 0.32e2 * t8 * t100 * t99 + 0.8e1 * t108 * t146 +
                         0.80e2 / 0.3e1 * t7 * v_0 * t112 - 0.8e1 * x * t115 * t112 + 0.32e2 * t158 * t12 * t156 -
                         0.32e2 * t18 * t37 * t156 + t177 + 0.16e2 / 0.3e1 * y * t180 * t135 -
                         0.16e2 / 0.3e1 * y * t179 * t73 + 0.32e2 * t120 * t187 - 0.32e2 * t124 * t190 +
                         0.8e1 * t131 * t196 - 0.8e1 * t136 * t199;
  const su2double t209 = t111 * w_0;
  const su2double t231 = w_0 * w_0;
  const su2double t232 = t231 * rho_0;
  const su2double t234 = z * t33;
  const su2double t244 = t58 * t23 * t57;
  const su2double t252 = (z * t52 * t44 / 0.2e1 + (-0.2e1 * t62 * t23 * t61 - t244) * t56 / 0.4e1) * TWall;
  const su2double t253 = t252 * t41;
  const su2double t255 = t74 + t76 - 0.2e1 * t77;
  const su2double t256 = t111 * t255;
  const su2double t263 = 0.32e2 * t8 * t119 * t94 - 0.32e2 * t8 * t123 * t99 + 0.80e2 / 0.3e1 * t7 * Viscosity * t209 +
                         0.8e1 * x * Viscosity * t130 - 0.8e1 * x * t135 * t134 + 0.32e2 * t17 * t119 * t187 -
                         0.32e2 * t17 * t123 * t190 + 0.8e1 * y * Viscosity * t196 - 0.8e1 * y * t135 * t199 +
                         0.32e2 * t234 * t12 * t232 - 0.32e2 * t24 * t37 * t232 + t253 +
                         0.16e2 / 0.3e1 * z * t256 * t135 - 0.16e2 / 0.3e1 * z * t255 * t73;
  const su2double t265 = t12 * u_0;
  const su2double t266 = x * t9;
  const su2double t270 = t115 * t112;
  const su2double t271 = t12 * v_0;
  const su2double t275 = t135 * t134;
  const su2double t276 = t12 * w_0;
  const su2double t280 = y * t9;
  const su2double t285 = t280 * t271;
  const su2double t288 = t135 * t199;
  const su2double t292 = z * t9;
  const su2double t300 = t292 * t276;
  const su2double t303 = Viscosity * t146;
  const su2double t305 = t9 * t12;
  const su2double t309 = t107 * t112;
  const su2double t310 = v_0 * t72;
  const su2double t316 = t305 * Viscosity * w_0;
  const su2double t319 = Viscosity * t134;
  const su2double t320 = w_0 * t72;
  const su2double t329 =
      -0.32e2 / 0.3e1 * t266 * t265 * t90 * t84 - 0.32e2 * t266 * t271 * t270 - 0.32e2 * t266 * t276 * t275 -
      0.32e2 * t280 * t265 * t270 + 0.64e2 / 0.3e1 * t285 * t180 * t135 - 0.32e2 * t280 * t276 * t288 -
      0.32e2 * t292 * t265 * t275 - 0.32e2 * t292 * t271 * t288 + 0.64e2 / 0.3e1 * t300 * t256 * t135 +
      0.32e2 * t305 * v_0 * t107 * t303 + 0.32e2 * t266 * t310 * t309 + 0.32e2 * t316 * x * t130 +
      0.32e2 * t266 * t320 * t319 + 0.32e2 * t305 * u_0 * t107 * Viscosity * t104;
  const su2double t330 = u_0 * t72;
  const su2double t334 = t179 * t73;
  const su2double t337 = v_0 * t111;
  const su2double t344 = Viscosity * t199;
  const su2double t364 = t255 * t73;
  const su2double t370 = Pi * Pi;
  const su2double t371 = t7 * t370;
  const su2double t375 = a_T2 * t46 + a_T1 / 0.4e1;
  const su2double t377 = t12 * t375 * t46;
  const su2double t379 = t58 * t58;
  const su2double t381 = t7 * t370 * t379;
  const su2double t385 = t375 * t58;
  const su2double t389 = t7 * Pi * t46;
  const su2double t399 = 0.1e1 / t42 / 0.4e1;
  const su2double t400 = Conductivity * t399;
  const su2double t406 = t46 * t46;
  const su2double t410 = Pi * t12 * t385 + a_T2 * (-t406 + 0.1e1 / 0.2e1) - t47 / 0.2e1;
  const su2double t412 = t12 * t89 * t410;
  const su2double t421 = t78 * t84;
  const su2double t422 = x * t12;
  const su2double t426 =
      0.32e2 * t280 * t330 * t309 - 0.64e2 / 0.3e1 * t285 * t334 - 0.64e2 / 0.3e1 * t280 * t337 * t334 +
      0.32e2 * t316 * y * t196 + 0.32e2 * t280 * t320 * t344 + 0.32e2 * t305 * Viscosity * u_0 * z * t130 +
      0.32e2 * t292 * t330 * t319 + 0.32e2 * t305 * Viscosity * v_0 * z * t196 + 0.32e2 * t292 * t310 * t344 -
      0.64e2 / 0.3e1 * t300 * t364 - 0.64e2 / 0.3e1 * t292 * t209 * t364 +
      0.4e1 * t400 * TWall * x * t12 * t7 *
          (t377 * x * t371 - t12 * a_T2 * x * t381 + Pi * t8 * t385 + 0.2e1 * a_T2 * t58 * x * t389 + t60 / 0.2e1) -
      0.4e1 * t400 * t2 * TWall * t412 + 0.64e2 / 0.3e1 * t305 * u_0 * t78 * t303 - 0.32e2 / 0.3e1 * t422 * t330 * t421;
  const su2double t457 = y * t12;
  const su2double t483 = t400 * TWall * t12;
  const su2double t486 = t410 * t4;
  const su2double t490 = z * t12;
  const su2double t506 = (0.3e1 / 0.4e1 + t51 * t56 / 0.4e1) * TWall;
  const su2double t508 = 1.0 / (Gamma - 1.0);
  const su2double t511 = t37 * t29;
  const su2double t513 = t37 * t155;
  const su2double t515 = t37 * t231;
  const su2double t521 = t508 * t506 * t41 + 0.8e1 * (t32 * t511 + t32 * t513 + t32 * t515) * rho_0 + t506 * t41;
  const su2double t522 = u_0 * t521;
  const su2double t527 =
      -0.32e2 * t422 * t310 * t309 - 0.32e2 * t422 * t320 * t319 +
      0.4e1 * t400 * y * TWall * t12 * t7 *
          (t377 * y * t371 - t12 * a_T2 * y * t381 + Pi * t17 * t385 + 0.2e1 * a_T2 * t58 * y * t389 + t168 / 0.2e1) -
      0.4e1 * t400 * TWall * t3 * t412 - 0.32e2 * t457 * t330 * t309 + 0.64e2 / 0.3e1 * t457 * t337 * t334 -
      0.32e2 * t457 * t320 * t344 +
      0.4e1 * t483 * t7 *
          (t377 * z * t371 - t12 * a_T2 * z * t381 + Pi * t23 * t385 + 0.2e1 * a_T2 * t58 * z * t389 + t244 / 0.2e1) *
          z -
      0.4e1 * t483 * t89 * t486 - 0.32e2 * t490 * t330 * t319 - 0.32e2 * t490 * t310 * t344 +
      0.64e2 / 0.3e1 * t490 * t209 * t364 + 0.32e2 / 0.3e1 * t266 * t330 * t421 + 0.4e1 * t10 * t522 -
      0.4e1 * t14 * t522;
  const su2double t528 = v_0 * t521;
  const su2double t533 = w_0 * t521;
  const su2double t541 = Conductivity * t399 * TWall;
  const su2double t545 = t9 * t13;
  const su2double t558 = t72 * t410;
  const su2double t570 = t12 * t29;
  const su2double t573 = t12 * t155;
  const su2double t576 = t12 * t231;
  const su2double t616 =
      0.4e1 * t18 * t528 - 0.4e1 * t20 * t528 + 0.4e1 * t24 * t533 - 0.4e1 * t26 * t533 +
      0.12e2 * t541 * t12 * t7 * t410 + 0.320e3 / 0.3e1 * t545 * t155 * t112 +
      0.320e3 / 0.3e1 * t545 * Viscosity * t111 * t231 + 0.64e2 * t545 * t29 * t112 +
      0.64e2 / 0.3e1 * t545 * t29 * t84 + 0.4e1 * t541 * t2 * t558 + 0.4e1 * t541 * t3 * t558 +
      0.4e1 * t541 * t72 * t486 +
      0.4e1 * t305 * u_0 *
          (t508 * t70 * t41 +
           0.16e2 * (-t10 * t511 - t10 * t513 - t10 * t515 + t34 * t570 + t34 * t573 + t34 * t576) * rho_0 + t71) +
      0.4e1 * t305 * v_0 *
          (t508 * t176 * t41 +
           0.16e2 * (t158 * t570 + t158 * t573 + t158 * t576 - t18 * t511 - t18 * t513 - t18 * t515) * rho_0 + t177) +
      0.4e1 * t305 * w_0 *
          (t508 * t252 * t41 +
           0.16e2 * (t234 * t570 + t234 * t573 + t234 * t576 - t24 * t511 - t24 * t513 - t24 * t515) * rho_0 + t253);

  /*--- Set the source term. Note the scaling for the correct non-dimensionalization. ---*/
  val_source[0] = 0.4e1 * t10 * t1 - 0.4e1 * t14 * t1 + 0.4e1 * t18 * t16 - 0.4e1 * t20 * t16 + 0.4e1 * t24 * t22 -
                  0.4e1 * t26 * t22;
  val_source[1] = t139;
  val_source[2] = t202;
  val_source[3] = t263;
  val_source[4] = t329 + t426 + t527 + t616;

  val_source[0] /= Density_Ref * Velocity_Ref;
  val_source[1] /= Pressure_Ref;
  val_source[2] /= Pressure_Ref;
  val_source[3] /= Pressure_Ref;
  val_source[4] /= Velocity_Ref * Pressure_Ref;
}

bool CMMSNSTwoHalfSpheresSolution::IsManufacturedSolution() const { return true; }
