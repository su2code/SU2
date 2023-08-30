/*!
 * \file CMMSNSTwoHalfCirclesSolution.cpp
 * \brief Implementations of the member functions of CMMSNSTwoHalfCirclesSolution.
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

#include "../../../include/toolboxes/MMS/CMMSNSTwoHalfCirclesSolution.hpp"

CMMSNSTwoHalfCirclesSolution::CMMSNSTwoHalfCirclesSolution() : CVerificationSolution() {}

CMMSNSTwoHalfCirclesSolution::CMMSNSTwoHalfCirclesSolution(unsigned short val_nDim, unsigned short val_nVar,
                                                           unsigned short val_iMesh, CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  /*--- Write a message that the solution is initialized for the manufactured
   solution for the Navier-Stokes equations between two half circles. ---*/
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the manufactured solution " << endl;
    cout << "         of the Navier-Stokes equations between two half circles!!!" << endl;
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
    SU2_MPI::Error("Steady mode must be selected for the MMS NS Two Half Circles case", CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::EULER && Kind_Solver != MAIN_SOLVER::NAVIER_STOKES &&
      Kind_Solver != MAIN_SOLVER::RANS && Kind_Solver != MAIN_SOLVER::FEM_EULER &&
      Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::FEM_RANS &&
      Kind_Solver != MAIN_SOLVER::FEM_LES)
    SU2_MPI::Error("Compressible flow equations must be selected for the MMS NS Two Half Circles case",
                   CURRENT_FUNCTION);

  if ((Kind_Solver != MAIN_SOLVER::NAVIER_STOKES) && (Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES))
    SU2_MPI::Error("Navier Stokes equations must be selected for the MMS NS Two Half Circles case", CURRENT_FUNCTION);

  if ((config->GetKind_FluidModel() != STANDARD_AIR) && (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the MMS NS Two Half Circles case", CURRENT_FUNCTION);

  if (config->GetKind_ViscosityModel() != VISCOSITYMODEL::CONSTANT)
    SU2_MPI::Error("Sutherland must be selected for viscosity for the MMS NS Two Half Circles case", CURRENT_FUNCTION);

  if (config->GetKind_ConductivityModel() != CONDUCTIVITYMODEL::CONSTANT_PRANDTL)
    SU2_MPI::Error("Constant Prandtl number must be selected for the MMS NS Two Half Circles case", CURRENT_FUNCTION);
}

CMMSNSTwoHalfCirclesSolution::~CMMSNSTwoHalfCirclesSolution() = default;

void CMMSNSTwoHalfCirclesSolution::GetBCState(const su2double* val_coords, const su2double val_t,
                                              su2double* val_solution) const {
  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CMMSNSTwoHalfCirclesSolution::GetSolution(const su2double* val_coords, const su2double val_t,
                                               su2double* val_solution) const {
  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  /* Determine the dimensional solution for the temperature. */
  const su2double Pi = PI_NUMBER;
  const su2double r = sqrt(x * x + y * y);
  const su2double fact = (r - 1.0) * (r - 1.0) / (a_T1 + a_T2);

  su2double T = 0.25 * TWall * (3.0 + fact * (a_T1 * cos(Pi * (r - 2.0)) + a_T2 * cos(Pi * (r - 2.0) * 2.0)));

  /* Determine the dimensional solution for the velocities. */
  su2double u = u_0 * (r - 1.0) * (2.0 - r) * 4.0;
  su2double v = v_0 * (r - 1.0) * (2.0 - r) * 4.0;

  /* Compute the pressure from the density and temperature. */
  su2double rho = rho_0;
  su2double p = rho * RGas * T;

  /* Determine the non-dimensional solution. */
  rho /= Density_Ref;
  p /= Pressure_Ref;
  u /= Velocity_Ref;
  v /= Velocity_Ref;

  /* Determine the non-dimensional conserved variables.
     Note that the implementation below is valid for both 2D and 3D. */
  val_solution[0] = rho;
  val_solution[1] = rho * u;
  val_solution[2] = rho * v;
  val_solution[3] = 0.0;
  val_solution[nDim + 1] = p / (Gamma - 1.0) + 0.5 * rho * (u * u + v * v);
}

void CMMSNSTwoHalfCirclesSolution::GetMMSSourceTerm(const su2double* val_coords, const su2double val_t,
                                                    su2double* val_source) const {
  /*--- Abbreviate Pi and the coordinates. ---*/
  const su2double Pi = PI_NUMBER;
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  /*--- The source code for the source terms is generated in Maple.
        See the file CMMSNSTwoHalfCirclesSolution.mw in the directory
        CreateMMSSourceTerms for the details how to do this. ---*/
  const su2double t1 = rho_0 * u_0;
  const su2double t2 = x * x;
  const su2double t3 = y * y;
  const su2double t4 = t2 + t3;
  const su2double t5 = sqrt(t4);
  const su2double t6 = 0.1e1 / t5;
  const su2double t7 = x * t6;
  const su2double t8 = 0.2e1 - t5;
  const su2double t9 = t8 * t7;
  const su2double t11 = t5 - 0.1e1;
  const su2double t12 = t6 * t11;
  const su2double t13 = x * t12;
  const su2double t15 = rho_0 * v_0;
  const su2double t16 = y * t6;
  const su2double t17 = t8 * t16;
  const su2double t19 = y * t12;
  const su2double t22 = u_0 * u_0;
  const su2double t23 = t22 * rho_0;
  const su2double t25 = t8 * t8;
  const su2double t26 = t6 * t25;
  const su2double t27 = x * t26;
  const su2double t30 = t11 * t11;
  const su2double t34 = rho_0 * RGas;
  const su2double t35 = a_T1 + a_T2;
  const su2double t36 = 0.1e1 / t35;
  const su2double t37 = t36 * t11;
  const su2double t38 = -t8 * Pi;
  const su2double t39 = cos(t38);
  const su2double t40 = t39 * a_T1;
  const su2double t41 = 0.2e1 * t38;
  const su2double t42 = cos(t41);
  const su2double t44 = t42 * a_T2 + t40;
  const su2double t45 = t6 * t44;
  const su2double t49 = t36 * t30;
  const su2double t50 = a_T1 * Pi;
  const su2double t51 = sin(t38);
  const su2double t53 = t51 * t7 * t50;
  const su2double t54 = a_T2 * Pi;
  const su2double t55 = sin(t41);
  const su2double t62 = (0.2e1 * x * t45 * t37 + (-0.2e1 * t55 * t7 * t54 - t53) * t49) * TWall;
  const su2double t64 = t62 * t34 / 0.4e1;
  const su2double t65 = 0.1e1 / t4;
  const su2double t66 = t65 * Viscosity;
  const su2double t67 = u_0 * x;
  const su2double t69 = v_0 * y;
  const su2double t70 = 0.2e1 * t67 - t69;
  const su2double t75 = -0.3e1 + 0.2e1 * t5;
  const su2double t76 = t75 * Viscosity;
  const su2double t81 = 0.1e1 / t5 / t4;
  const su2double t82 = t81 * t70;
  const su2double t86 = t11 * t1;
  const su2double t87 = v_0 * t25;
  const su2double t91 = t30 * t1;
  const su2double t92 = v_0 * t8;
  const su2double t99 = t3 * u_0;
  const su2double t102 = t81 * t75;
  const su2double t105 = v_0 * x;
  const su2double t113 = (0.4e1 * y * t102 * t105 - 0.8e1 * y * t65 * t105 - 0.4e1 * t6 * t75 * u_0 +
                          0.4e1 * t102 * t99 - 0.8e1 * t65 * t99) *
                         Viscosity;
  const su2double t121 = u_0 * y;
  const su2double t131 = t2 * v_0;
  const su2double t137 = (0.4e1 * x * t102 * t121 - 0.8e1 * x * t65 * t121 - 0.4e1 * t6 * t75 * v_0 +
                          0.4e1 * t102 * t131 - 0.8e1 * t65 * t131) *
                         Viscosity;
  const su2double t138 = v_0 * v_0;
  const su2double t139 = t138 * rho_0;
  const su2double t141 = y * t26;
  const su2double t151 = t51 * t16 * t50;
  const su2double t158 = (0.2e1 * y * t45 * t37 + (-0.2e1 * t55 * t16 * t54 - t151) * t49) * TWall;
  const su2double t160 = t158 * t34 / 0.4e1;
  const su2double t161 = -0.3e1 / 0.2e1 + t5;
  const su2double t167 = t67 - 0.2e1 * t69;
  const su2double t168 = t65 * t167;
  const su2double t172 = t161 * t167;
  const su2double t173 = Viscosity * t81;
  const su2double t182 = a_T2 * t39 + a_T1 / 0.4e1;
  const su2double t183 = t182 * t51;
  const su2double t186 = t39 * t39;
  const su2double t190 = -Pi * t11 * t183 + a_T2 * (t186 - 0.1e1 / 0.2e1) + t40 / 0.2e1;
  const su2double t191 = 0.1e1 / t35 / 0.4e1;
  const su2double t193 = Conductivity * t191 * t190;
  const su2double t201 = t11 * TWall * t6;
  const su2double t208 = t6 * t75;
  const su2double t212 = 0.4e1 * (-t208 * t105 - t208 * t121) * Viscosity;
  const su2double t213 = v_0 * t212;
  const su2double t218 = u_0 * t212;
  const su2double t230 = t11 * v_0;
  const su2double t231 = y * t8;
  const su2double t236 = t11 * u_0;
  const su2double t237 = x * t8;
  const su2double t246 = (t44 * t49 + 0.3e1) * TWall;
  const su2double t248 = 1.0 / (Gamma - 1.0);
  const su2double t252 = t30 * t22;
  const su2double t254 = t30 * t138;
  const su2double t261 = t248 * t246 * t34 / 0.4e1 + 0.8e1 * (t25 * t252 + t25 * t254) * rho_0 + t246 * t34 / 0.4e1;
  const su2double t262 = v_0 * t261;
  const su2double t267 = -0.4e1 * t193 * t2 * TWall * t65 - 0.4e1 * t193 * TWall * t65 * t3 - 0.8e1 * t193 * t201 +
                         0.64e2 / 0.3e1 * t8 * t12 * t22 * t76 - 0.4e1 * t9 * t213 + 0.4e1 * t13 * t213 -
                         0.4e1 * t17 * t218 + 0.4e1 * t19 * t218 +
                         0.128e3 / 0.3e1 * t8 * t11 * Viscosity * t6 * t161 * t138 +
                         0.64e2 / 0.3e1 * t231 * t230 * t173 * t172 - 0.32e2 / 0.3e1 * t237 * t236 * t82 * t76 -
                         0.4e1 * t8 * t236 * t113 - 0.4e1 * t19 * t262 + 0.4e1 * t17 * t262;
  const su2double t271 = u_0 * t261;
  const su2double t277 = Pi * Pi;
  const su2double t278 = t6 * t277;
  const su2double t281 = t11 * t182 * t39;
  const su2double t283 = t51 * t51;
  const su2double t285 = t6 * t277 * t283;
  const su2double t292 = t6 * Pi * t39;
  const su2double t300 = Conductivity * t191;
  const su2double t310 = t70 * t76;
  const su2double t311 = u_0 * t65;
  const su2double t316 = Viscosity * v_0;
  const su2double t317 = t8 * t11;
  const su2double t321 = t65 * t172;
  const su2double t364 = t11 * t22;
  const su2double t367 = t11 * t138;
  const su2double t391 =
      -0.4e1 * t8 * t230 * t137 + 0.4e1 * t9 * t271 - 0.4e1 * t13 * t271 -
      0.4e1 * t300 *
          (-t281 * y * t278 + t11 * a_T2 * y * t285 - Pi * t16 * t183 - 0.2e1 * a_T2 * t51 * y * t292 - t151 / 0.2e1) *
          t11 * TWall * t16 +
      0.4e1 * t300 * t190 * t11 * TWall * t81 * t3 + 0.32e2 / 0.3e1 * t237 * t311 * t310 -
      0.64e2 / 0.3e1 * t317 * t316 * y * t168 - 0.64e2 / 0.3e1 * t231 * t316 * t321 +
      0.64e2 / 0.3e1 * y * t11 * t316 * t321 -
      0.4e1 * t300 * x *
          (-t281 * x * t278 + t11 * a_T2 * x * t285 - Pi * t7 * t183 - 0.2e1 * a_T2 * t51 * x * t292 - t53 / 0.2e1) *
          t201 +
      0.4e1 * t300 * t2 * t190 * t11 * TWall * t81 + 0.64e2 / 0.3e1 * t317 * u_0 * t70 * x * t66 -
      0.32e2 / 0.3e1 * x * t11 * t311 * t310 +
      0.4e1 * t317 * u_0 *
          (t248 * t62 * t34 / 0.4e1 + 0.16e2 * (-t9 * t252 - t9 * t254 + t27 * t364 + t27 * t367) * rho_0 + t64) +
      0.4e1 * t317 * v_0 *
          (t248 * t158 * t34 / 0.4e1 + 0.16e2 * (t141 * t364 + t141 * t367 - t17 * t252 - t17 * t254) * rho_0 + t160);

  /*--- Set the source term, which is valid for both 2D and 3D cases.
        Note the scaling for the correct non-dimensionalization. ---*/
  val_source[0] = -0.4e1 * t13 * t1 + 0.4e1 * t9 * t1 + 0.4e1 * t17 * t15 - 0.4e1 * t19 * t15;
  val_source[1] = 0.32e2 * t27 * t11 * t23 - 0.32e2 * t9 * t30 * t23 + t64 + 0.16e2 / 0.3e1 * t70 * x * t66 +
                  0.16e2 / 0.3e1 * t6 * u_0 * t76 - 0.8e1 / 0.3e1 * x * t82 * t76 + 0.32e2 * t16 * t87 * t86 -
                  0.32e2 * t16 * t92 * t91 - t113;
  val_source[2] = 0.32e2 * t7 * t87 * t86 - 0.32e2 * t7 * t92 * t91 - t137 + 0.32e2 * t141 * t11 * t139 -
                  0.32e2 * t17 * t30 * t139 + t160 + 0.32e2 / 0.3e1 * Viscosity * t6 * t161 * v_0 -
                  0.16e2 / 0.3e1 * y * Viscosity * t168 + 0.16e2 / 0.3e1 * y * t173 * t172;
  val_source[3] = 0.0;
  val_source[nDim + 1] = t267 + t391;

  val_source[0] /= Density_Ref * Velocity_Ref;
  val_source[1] /= Pressure_Ref;
  val_source[2] /= Pressure_Ref;
  val_source[nDim + 1] /= Velocity_Ref * Pressure_Ref;
}

bool CMMSNSTwoHalfCirclesSolution::IsManufacturedSolution() const { return true; }
