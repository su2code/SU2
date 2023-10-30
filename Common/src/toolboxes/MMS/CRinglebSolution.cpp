/*!
 * \file CRinglebSolution.cpp
 * \brief Implementations of the member functions of CRinglebSolution.
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

#include "../../../include/toolboxes/MMS/CRinglebSolution.hpp"

CRinglebSolution::CRinglebSolution() : CVerificationSolution() {}

CRinglebSolution::CRinglebSolution(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_iMesh,
                                   CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  /*--- Write a message that the solution is initialized for the
   Ringleb test case. ---*/
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the Ringleb case!!!" << endl;
    cout << endl << flush;
  }

  /*--- Useful coefficients in which Gamma is present. ---*/
  Gamma = config->GetGamma();
  Gm1 = Gamma - 1.0;
  tovGm1 = 2.0 / Gm1;
  tGamOvGm1 = Gamma * tovGm1;

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if (config->GetTime_Marching() != TIME_MARCHING::STEADY)
    SU2_MPI::Error("Steady mode must be selected for the Ringleb case", CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::EULER && Kind_Solver != MAIN_SOLVER::NAVIER_STOKES &&
      Kind_Solver != MAIN_SOLVER::RANS && Kind_Solver != MAIN_SOLVER::FEM_EULER &&
      Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::FEM_RANS &&
      Kind_Solver != MAIN_SOLVER::FEM_LES)
    SU2_MPI::Error("Compressible flow equations must be selected for the Ringleb case", CURRENT_FUNCTION);

  if ((Kind_Solver != MAIN_SOLVER::EULER) && (Kind_Solver != MAIN_SOLVER::FEM_EULER))
    SU2_MPI::Error("Euler equations must be selected for the Ringleb case", CURRENT_FUNCTION);

  if ((config->GetKind_FluidModel() != STANDARD_AIR) && (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the Ringleb case", CURRENT_FUNCTION);
}

CRinglebSolution::~CRinglebSolution() = default;

void CRinglebSolution::GetBCState(const su2double* val_coords, const su2double val_t, su2double* val_solution) const {
  /*--- The exact solution is prescribed on the boundaries for the
        Ringleb flow. Note that a (much) more difficult test case is to
        use inviscid wall boundary conditions for the inner and outer
        walls of the channel. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CRinglebSolution::GetSolution(const su2double* val_coords, const su2double val_t, su2double* val_solution) const {
  /* Easier storage of the coordinates and abbreviate y*y. */
  const su2double x = val_coords[0], y = val_coords[1];
  const su2double y2 = y * y;

  /* Initial guess for q (velocity magnitude) and k (streamline parameter). */
  su2double k = 1.2;
  su2double q = 1.0;

  /* Newton algorithm to solve for the variables q and k for the given x and y. */
  const int iterMax = 500;
  su2double duMaxPrev = 10.0;

  int iter;
  for (iter = 0; iter < iterMax; ++iter) {
    /* Compute the speed of sound, the density, the parameter JJ
       and its derivatives w.r.t. q. */
    const su2double a = sqrt(1.0 - 0.5 * Gm1 * q * q);
    const su2double rho = pow(a, tovGm1);
    const su2double JJ =
        1.0 / a + 1.0 / (3.0 * a * a * a) + 1.0 / (5.0 * a * a * a * a * a) - 0.5 * log((1.0 + a) / (1.0 - a));

    const su2double dadq = -0.5 * Gm1 * q / a;
    const su2double drhodq = 2.0 * rho * dadq / (Gm1 * a);
    const su2double dJJdq = dadq / (pow(a, 6) * (a * a - 1.0));

    /* Determine the values of the nonlinear equations to solve
       and its corresponding Jacobian matrix. */
    const su2double y2c = (k * k - q * q) / (k * k * k * k * rho * rho * q * q);
    const su2double f[] = {(2.0 / (k * k) - 1.0 / (q * q)) / (2.0 * rho) - 0.5 * JJ - x, y2c - y2};
    su2double Jac[2][2];
    Jac[0][0] = -(1.0 / (k * k) - 0.50 / (q * q)) * drhodq / (rho * rho) + 1.0 / (rho * q * q * q) - 0.5 * dJJdq;
    Jac[0][1] = -2.0 / (rho * k * k * k);
    Jac[1][0] = -2.0 / (k * k * rho * rho * q * q * q) - 2.0 * y2c * drhodq / rho;
    Jac[1][1] = (4.0 * q * q - 2.0 * k * k) / (k * k * k * k * k * rho * rho * q * q);

    /* Determine the update dU. */
    const su2double det = Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0];
    const su2double dU[] = {(f[0] * Jac[1][1] - f[1] * Jac[0][1]) / det, (f[1] * Jac[0][0] - f[0] * Jac[1][0]) / det};

    /* Determine the underrelaxation coefficient alp. */
    const su2double dUMax = max(fabs(dU[0]), fabs(dU[1]));
    su2double alp = 1.0;
    if (dUMax > 1.0)
      alp = 0.04;
    else if (dUMax > 0.1)
      alp = 0.2;

    /* Update q and k. */
    q -= alp * dU[0];
    k -= alp * dU[1];

    /* Convergence check, which is independent of the precision used. */
    if ((dUMax < 1.e-3) && (dUMax >= duMaxPrev)) break;
    duMaxPrev = dUMax;
  }

  /* Check if the Newton algorithm actually converged. */
  if (iter == iterMax) SU2_MPI::Error("Newton algorithm did not converge", CURRENT_FUNCTION);

  /* Compute the speed of sound, density and pressure. */
  const su2double a = sqrt(1.0 - 0.5 * Gm1 * q * q);
  const su2double rho = pow(a, tovGm1);
  const su2double p = pow(a, tGamOvGm1) / Gamma;

  /* Determine the derivative of x w.r.t. q and ydxdq. */
  const su2double dadq = -0.5 * Gm1 * q / a;
  const su2double drhodq = 2.0 * rho * dadq / (Gm1 * a);
  const su2double dJJdq = dadq / (pow(a, 6) * (a * a - 1.0));

  const su2double dxdq =
      -(1.0 / (k * k) - 0.5 / (q * q)) * drhodq / (rho * rho) + 1.0 / (rho * q * q * q) - 0.5 * dJJdq;
  const su2double ydxdq = y * dxdq;

  /* Determine the derivative of 1/2 y2 w.r.t. q, which is ydydq. The reason is
     that ydydq is always well defined, while dydq is singular for y = 0. */
  const su2double ydydq =
      drhodq * (q * q - k * k) / (k * k * k * k * rho * rho * rho * q * q) - 1.0 / (k * k * rho * rho * q * q * q);

  /* Determine the direction of the streamline. */
  const su2double vecLen = sqrt(ydxdq * ydxdq + ydydq * ydydq);

  su2double velDir[] = {ydxdq / vecLen, ydydq / vecLen};
  if (velDir[1] > 0.0) {
    velDir[0] = -velDir[0];
    velDir[1] = -velDir[1];
  }

  /* Compute the conservative variables. Note that both 2D and 3D
     cases are treated correctly. */
  val_solution[0] = rho;
  val_solution[1] = rho * q * velDir[0];
  val_solution[2] = rho * q * velDir[1];
  val_solution[3] = 0.0;
  val_solution[nVar - 1] = p / Gm1 + 0.5 * rho * q * q;
}
