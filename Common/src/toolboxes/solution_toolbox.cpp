/*!
 * \file solution_toolbox.cpp
 * \brief Implementations of the exact solution classes.
 * \author T. Economon, E. van der Weide
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../include/toolboxes/solution_toolbox.hpp"

CVerificationSolution::CVerificationSolution(void) { }

CVerificationSolution::CVerificationSolution(unsigned short val_nDim,
                                             unsigned short val_nVar,
                                             CConfig*       config) {
  
  /*--- Store the rank and size for the calculation. ---*/
  
  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();
  
  /*--- Store the dimension and number of variables. ---*/
  
  nDim = val_nDim;
  nVar = val_nVar;
  
}

CVerificationSolution::~CVerificationSolution(void) { }

void CVerificationSolution::GetSolution(const unsigned short val_nParams,
                                        const su2double      *val_params,
                                        const su2double      *val_coords,
                                        const su2double      val_t,
                                        su2double            *val_solution) {

  SU2_MPI::Error("Function must be overwritten by the derived class", CURRENT_FUNCTION);
}

void CVerificationSolution::GetInitialCondition(const unsigned short val_nParams,
                                                const su2double      *val_params,
                                                const su2double      *val_coords,
                                                su2double            *val_solution) {
  
  /*--- Initial conditions call the GetSolution() method at t = 0. ---*/
  GetSolution(val_nParams, val_params, val_coords, 0.0, val_solution);
}

void CVerificationSolution::GetBCState(const unsigned short val_nParams,
                                       const su2double      *val_params,
                                       const su2double      *val_coords,
                                       const su2double      val_t,
                                       su2double            *val_solution) {

  SU2_MPI::Error("Function must be overwritten by the derived class", CURRENT_FUNCTION);
}

void CVerificationSolution::GetMMSSourceTerm(const unsigned short val_nParams,
                                             const su2double      *val_params,
                                             const su2double      *val_coords,
                                             const su2double      val_t,
                                             su2double            *val_source) {

  /* Default implementation of the source terms for the method of manufactured
     solutions. Simply set them to zero. */
  for(unsigned short iVar=0; iVar<nVar; ++iVar)
    val_source[iVar] = 0.0;
}

bool CVerificationSolution::IsManufacturedSolution(void) {return false;}

bool CVerificationSolution::ExactSolutionKnown(void) {return true;}

void CVerificationSolution::GetLocalError(const unsigned short val_nParams,
                                          const su2double      *val_params,
                                          const su2double      *val_coords,
                                          const su2double      val_t,
                                          const su2double      *val_solution,
                                          su2double            *val_error) {
  
  /*--- Get the value of the verification solution first.
        Use val_error to store this solution. ---*/
  
  GetSolution(val_nParams, val_params, val_coords, val_t, val_error);
  
  /*--- Compute the local error as the difference between the current
   numerical solution and the verification solution. ---*/
  
  for (unsigned short iVar=0; iVar<nVar; ++iVar)
    val_error[iVar] = val_solution[iVar] - val_error[iVar];
  
}

CInviscidVortexSolution::CInviscidVortexSolution(void) : CVerificationSolution() { }

CInviscidVortexSolution::CInviscidVortexSolution(unsigned short val_nDim,
                                                 unsigned short val_nVar,
                                                 CConfig*       config)
  : CVerificationSolution(val_nDim, val_nVar, config) {

  /*--- Write a message that the solution is initialized for the
   inviscid vortex test case. ---*/
  if (rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the inviscid vortex case!!!" << endl;
    cout << endl << flush;
  }

  /*--- Store the inviscid vortex specific parameters here. ---*/
  x0Vortex    = -0.5;     // Initial x-coordinate of the vortex center.
  y0Vortex    =  0.0;     // Initial y-coordinate of the vortex center.
  RVortex     =  0.1;     // Radius of the vortex.
  epsVortex   =  1.0;     // Strength of the vortex.

  /* Get the Mach number and advection angle (in degrees). */
  MachVortex  = config->GetMach();
  thetaVortex = config->GetAoA();

  /*--- Useful coefficients in which Gamma is present. ---*/
  Gamma    = config->GetGamma();
  Gm1      = Gamma - 1.0;
  ovGm1    = 1.0/Gm1;
  gamOvGm1 = ovGm1*Gamma;

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if((config->GetUnsteady_Simulation() != TIME_STEPPING) &&
     (config->GetUnsteady_Simulation() != DT_STEPPING_1ST) &&
     (config->GetUnsteady_Simulation() != DT_STEPPING_2ND))
    SU2_MPI::Error("Unsteady mode must be selected for the inviscid vortex",
                   CURRENT_FUNCTION);

  if(config->GetKind_Regime() != COMPRESSIBLE)
    SU2_MPI::Error("Compressible flow equations must be selected for the inviscid vortex",
                   CURRENT_FUNCTION);

  if((config->GetKind_Solver() != EULER) &&
     (config->GetKind_Solver() != FEM_EULER))
    SU2_MPI::Error("Euler equations must be selected for the inviscid vortex",
                   CURRENT_FUNCTION);

  if((config->GetKind_FluidModel() != STANDARD_AIR) &&
     (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the inviscid vortex",
                   CURRENT_FUNCTION);

  if(fabs(config->GetPressure_FreeStreamND() - 1.0) > 1.e-8)
    SU2_MPI::Error("Free-stream pressure must be 1.0 for the inviscid vortex",
                   CURRENT_FUNCTION);

  if(fabs(config->GetDensity_FreeStreamND() - 1.0) > 1.e-8)
    SU2_MPI::Error("Free-stream density must be 1.0 for the inviscid vortex",
                   CURRENT_FUNCTION);
}

CInviscidVortexSolution::~CInviscidVortexSolution(void) { }

void CInviscidVortexSolution::GetBCState(const unsigned short val_nParams,
                                         const su2double      *val_params,
                                         const su2double      *val_coords,
                                         const su2double      val_t,
                                         su2double            *val_solution) {

  /*--- For the case that the inviscid vortex is run with boundary
        conditions (other possibility is with periodic conditions),
        the exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_nParams, val_params, val_coords, val_t, val_solution);
}

void CInviscidVortexSolution::GetSolution(const unsigned short val_nParams,
                                          const su2double      *val_params,
                                          const su2double      *val_coords,
                                          const su2double      val_t,
                                          su2double            *val_solution) {

  /* Compute the free stream velocities in x- and y-direction. */
  const su2double VelInf = MachVortex*sqrt(Gamma);
  const su2double uInf   = VelInf*cos(thetaVortex*PI_NUMBER/180.0);
  const su2double vInf   = VelInf*sin(thetaVortex*PI_NUMBER/180.0);

  /* Compute the coordinates relative to the center of the vortex. */
  const su2double dx = val_coords[0] - (x0Vortex + val_t*uInf);
  const su2double dy = val_coords[1] - (y0Vortex + val_t*vInf);

  /* Compute the components of the velocity. */
  su2double f  = 1.0 - (dx*dx + dy*dy)/(RVortex*RVortex);
  su2double t1 = epsVortex*dy*exp(0.5*f)/(2.0*PI_NUMBER*RVortex);
  su2double u  = uInf - VelInf*t1;

  t1          = epsVortex*dx*exp(0.5*f)/(2.0*PI_NUMBER*RVortex);
  su2double v = vInf + VelInf*t1;

  /* Compute the density and the pressure. */
  t1 = 1.0 - epsVortex*epsVortex*Gm1
     *       MachVortex*MachVortex*exp(f)/(8.0*PI_NUMBER*PI_NUMBER);

  su2double rho = pow(t1,ovGm1);
  su2double p   = pow(t1,gamOvGm1);

  /* Compute the conservative variables. Note that both 2D and 3D
     cases are treated correctly. */
  val_solution[0]      = rho;
  val_solution[1]      = rho*u;
  val_solution[2]      = rho*v;
  val_solution[3]      = 0.0;
  val_solution[nVar-1] = p*ovGm1 + 0.5*rho*(u*u + v*v);
}

CRinglebSolution::CRinglebSolution(void) : CVerificationSolution() { }

CRinglebSolution::CRinglebSolution(unsigned short val_nDim,
                                   unsigned short val_nVar,
                                   CConfig*       config)
  : CVerificationSolution(val_nDim, val_nVar, config) {

  /*--- Write a message that the solution is initialized for the
   Ringleb test case. ---*/
  if (rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the Ringleb case!!!" << endl;
    cout << endl << flush;
  }

  /*--- Useful coefficients in which Gamma is present. ---*/
  Gamma     = config->GetGamma();
  Gm1       = Gamma - 1.0;
  tovGm1    = 2.0/Gm1;
  tGamOvGm1 = Gamma*tovGm1;

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if(config->GetUnsteady_Simulation() != STEADY)
    SU2_MPI::Error("Steady mode must be selected for the Ringleb case",
                   CURRENT_FUNCTION);

  if(config->GetKind_Regime() != COMPRESSIBLE)
    SU2_MPI::Error("Compressible flow equations must be selected for the Ringleb case",
                   CURRENT_FUNCTION);

  if((config->GetKind_Solver() != EULER) &&
     (config->GetKind_Solver() != FEM_EULER))
    SU2_MPI::Error("Euler equations must be selected for the Ringleb case",
                   CURRENT_FUNCTION);

  if((config->GetKind_FluidModel() != STANDARD_AIR) &&
     (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the Ringleb case",
                   CURRENT_FUNCTION);
}

CRinglebSolution::~CRinglebSolution(void) { }

void CRinglebSolution::GetBCState(const unsigned short val_nParams,
                                  const su2double      *val_params,
                                  const su2double      *val_coords,
                                  const su2double      val_t,
                                  su2double            *val_solution) {

  /*--- The exact solution is prescribed on the boundaries for the
        Ringleb flow. Note that a (much) more difficult test case is to
        use inviscid wall boundary conditions for the inner and outer
        walls of the channel. ---*/
  GetSolution(val_nParams, val_params, val_coords, val_t, val_solution);
}

void CRinglebSolution::GetSolution(const unsigned short val_nParams,
                                   const su2double      *val_params,
                                   const su2double      *val_coords,
                                   const su2double      val_t,
                                   su2double            *val_solution) {

  /* Easier storage of the coordinates and abbreviate y*y. */
  const su2double x  = val_coords[0], y = val_coords[1];
  const su2double y2 = y*y;

  /* Initial guess for q (velocity magnitude) and k (streamline parameter). */
  su2double k = 1.2;
  su2double q = 1.0;

  /* Newton algorithm to solve for the variables q and k for the given x and y. */
  const int iterMax = 500;
  su2double duMaxPrev = 10.0;

  int iter;
  for(iter=0; iter<iterMax; ++iter) {

    /* Compute the speed of sound, the density, the parameter JJ
       and its derivatives w.r.t. q. */
    const su2double a   = sqrt(1.0 - 0.5*Gm1*q*q);
    const su2double rho = pow(a,tovGm1);
    const su2double JJ  = 1.0/a + 1.0/(3.0*a*a*a) + 1.0/(5.0*a*a*a*a*a)
                        - 0.5*log((1.0+a)/(1.0-a));

    const su2double dadq   = -0.5*Gm1*q/a;
    const su2double drhodq =  2.0*rho*dadq/(Gm1*a);
    const su2double dJJdq  =  dadq/(pow(a,6)*(a*a-1.0));

    /* Determine the values of the nonlinear equations to solve
       and its corresponding Jacobian matrix. */
    const su2double y2c = (k*k - q*q)/(k*k*k*k*rho*rho*q*q);
    const su2double f[] = {(2.0/(k*k) - 1.0/(q*q))/(2.0*rho) - 0.5*JJ - x,
                           y2c - y2};
    su2double Jac[2][2];
    Jac[0][0] = -(1.0/(k*k) - 0.50/(q*q))*drhodq/(rho*rho)
              + 1.0/(rho*q*q*q) - 0.5*dJJdq;
    Jac[0][1] = -2.0/(rho*k*k*k);
    Jac[1][0] = -2.0/(k*k*rho*rho*q*q*q) - 2.0*y2c*drhodq/rho;
    Jac[1][1] = (4.0*q*q - 2.0*k*k)/(k*k*k*k*k*rho*rho*q*q);

    /* Determine the update dU. */
    const su2double det  = Jac[0][0]*Jac[1][1] - Jac[0][1]*Jac[1][0];
    const su2double dU[] = {(f[0]*Jac[1][1] - f[1]*Jac[0][1])/det,
                            (f[1]*Jac[0][0] - f[0]*Jac[1][0])/det};

    /* Determine the underrelaxation coefficient alp. */
    const su2double dUMax = max(fabs(dU[0]), fabs(dU[1]));
    su2double alp = 1.0;
    if(     dUMax > 1.0) alp = 0.04;
    else if(dUMax > 0.1) alp = 0.2;

    /* Update q and k. */
    q -= alp*dU[0];
    k -= alp*dU[1];

    /* Convergence check, which is independent of the precision used. */
    if((dUMax < 1.e-3) && (dUMax >= duMaxPrev)) break;
    duMaxPrev = dUMax;
  }

  /* Check if the Newton algorithm actually converged. */
  if(iter == iterMax)
    SU2_MPI::Error("Newton algorithm did not converge", CURRENT_FUNCTION);
  
  /* Compute the speed of sound, density and pressure. */
  const su2double a   = sqrt(1.0 - 0.5*Gm1*q*q);
  const su2double rho = pow(a,tovGm1);
  const su2double p   = pow(a,tGamOvGm1)/Gamma;

  /* Determine the derivative of x w.r.t. q and ydxdq. */
  const su2double dadq   = -0.5*Gm1*q/a;
  const su2double drhodq =  2.0*rho*dadq/(Gm1*a);
  const su2double dJJdq  =  dadq/(pow(a,6)*(a*a-1.0));

  const su2double dxdq  = -(1.0/(k*k) - 0.5/(q*q))*drhodq/(rho*rho)
                        +   1.0/(rho*q*q*q) - 0.5*dJJdq;
  const su2double ydxdq = y*dxdq;

  /* Determine the derivative of 1/2 y2 w.r.t. q, which is ydydq. The reason is
     that ydydq is always well defined, while dydq is singular for y = 0. */
  const su2double ydydq = drhodq*(q*q-k*k)/(k*k*k*k*rho*rho*rho*q*q)
                        - 1.0/(k*k*rho*rho*q*q*q);

  /* Determine the direction of the streamline. */
  const su2double vecLen = sqrt(ydxdq*ydxdq + ydydq*ydydq);

  su2double velDir[] = {ydxdq/vecLen, ydydq/vecLen};
  if(velDir[1] > 0.0){velDir[0] = -velDir[0]; velDir[1] = -velDir[1];}

  /* Compute the conservative variables. Note that both 2D and 3D
     cases are treated correctly. */
  val_solution[0]      = rho;
  val_solution[1]      = rho*q*velDir[0];
  val_solution[2]      = rho*q*velDir[1];
  val_solution[3]      = 0.0;
  val_solution[nVar-1] = p/Gm1 + 0.5*rho*q*q;
}

CNSUnitQuadSolution::CNSUnitQuadSolution(void) : CVerificationSolution() { }

CNSUnitQuadSolution::CNSUnitQuadSolution(unsigned short val_nDim,
                                         unsigned short val_nVar,
                                         CConfig*       config)
  : CVerificationSolution(val_nDim, val_nVar, config) {

  /*--- Write a message that the solution is initialized for the
   Navier-Stokes case on a unit quad. Note that heat conduction
   is neglected for this case. ---*/
  
  if (rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the Navier Stokes equations " << endl;
    cout << "         on a unit quad. Note that heat conduction " << endl;
    cout << "         is neglected for this case!!!" << endl;
    cout << endl << flush;
  }

  /*--- Coefficients, needed to determine the solution. ---*/
  Gm1       = config->GetGamma() - 1.0;
  flowAngle = config->GetAoA()*PI_NUMBER/180.0;
  Viscosity = config->GetMu_ConstantND();

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if(config->GetUnsteady_Simulation() != STEADY)
    SU2_MPI::Error("Steady mode must be selected for the NS Unit Quad case",
                   CURRENT_FUNCTION);

  if(config->GetKind_Regime() != COMPRESSIBLE)
    SU2_MPI::Error("Compressible flow equations must be selected for the NS Unit Quad case",
                   CURRENT_FUNCTION);

  if((config->GetKind_Solver() != NAVIER_STOKES) &&
     (config->GetKind_Solver() != FEM_NAVIER_STOKES))
    SU2_MPI::Error("Navier Stokes equations must be selected for the NS Unit Quad case",
                   CURRENT_FUNCTION);

  if(config->GetKind_FluidModel() != IDEAL_GAS)
    SU2_MPI::Error("Ideal gas must be selected for the NS Unit Quad case",
                   CURRENT_FUNCTION);

  if(fabs(Gm1-0.5) > 1.e-8)
    SU2_MPI::Error("Gamma must be 1.5 for the NS Unit Quad case",
                   CURRENT_FUNCTION);

  if(config->GetKind_ViscosityModel() != CONSTANT_VISCOSITY)
    SU2_MPI::Error("Constant viscosity must be selected for the NS Unit Quad case",
                   CURRENT_FUNCTION);

  if(config->GetKind_ConductivityModel() != CONSTANT_PRANDTL)
    SU2_MPI::Error("Constant Prandtl number must be selected for the NS Unit Quad case",
                   CURRENT_FUNCTION);

  if(config->GetPrandtl_Lam() < 1.e+20)
     SU2_MPI::Error("Laminar Prandtl number must be larger than 1.e+20 for the NS Unit Quad case",
                   CURRENT_FUNCTION);
}

CNSUnitQuadSolution::~CNSUnitQuadSolution(void) { }

void CNSUnitQuadSolution::GetBCState(const unsigned short val_nParams,
                                     const su2double      *val_params,
                                     const su2double      *val_coords,
                                     const su2double      val_t,
                                     su2double            *val_solution) {

  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_nParams, val_params, val_coords, val_t, val_solution);
}

void CNSUnitQuadSolution::GetSolution(const unsigned short val_nParams,
                                      const su2double      *val_params,
                                      const su2double      *val_coords,
                                      const su2double      val_t,
                                      su2double            *val_solution) {

  /*--- Compute the flow direction and the coordinates in
        the rotated frame. ---*/
  const su2double cosFlowAngle = cos(flowAngle);
  const su2double sinFlowAngle = sin(flowAngle);

  const su2double xTilde = val_coords[0]*cosFlowAngle - val_coords[1]*sinFlowAngle;
  const su2double yTilde = val_coords[0]*sinFlowAngle + val_coords[1]*cosFlowAngle;

  /*--- Compute the exact solution for this case. Note that it works
        both in 2D and 3D. ---*/
  val_solution[0]      =  1.0;
  val_solution[1]      =  cosFlowAngle*yTilde*yTilde;
  val_solution[2]      = -sinFlowAngle*yTilde*yTilde;
  val_solution[3]      =  0.0;
  val_solution[nVar-1] =  (2.0*Viscosity*xTilde + 10.0)/Gm1
                       +  0.5*yTilde*yTilde*yTilde*yTilde;
}

CTGVSolution::CTGVSolution(void) : CVerificationSolution() { }

CTGVSolution::CTGVSolution(unsigned short val_nDim,
                           unsigned short val_nVar,
                           CConfig*       config)
  : CVerificationSolution(val_nDim, val_nVar, config) {
  
  /*--- Write a message that the solution is initialized for the
   Taylor-Green vortex test case. ---*/
  
  if (rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the Taylor-Green vortex case!!!" << endl;
    cout << endl << flush;
  }
  
  /*--- Store TGV specific parameters here. ---*/
  
  tgvLength   = 1.0;     // Taylor-Green length scale.
  tgvVelocity = 1.0;     // Taylor-Green velocity.
  tgvDensity  = 1.0;     // Taylor-Green density.
  tgvPressure = 100.0;   // Taylor-Green pressure.
  
  /*--- Useful coefficient in which Gamma is present. ---*/
  
  ovGm1 = 1.0/(config->GetGamma() - 1.0);
  
  /*--- Perform some sanity and error checks for this solution here. ---*/

  if((config->GetUnsteady_Simulation() != TIME_STEPPING) &&
     (config->GetUnsteady_Simulation() != DT_STEPPING_1ST) &&
     (config->GetUnsteady_Simulation() != DT_STEPPING_2ND))
    SU2_MPI::Error("Unsteady mode must be selected for the Taylor Green Vortex",
                   CURRENT_FUNCTION); 

  if(config->GetKind_Regime() != COMPRESSIBLE)
    SU2_MPI::Error("Compressible flow equations must be selected for the Taylor Green Vortex",
                   CURRENT_FUNCTION);

  if((config->GetKind_Solver() != NAVIER_STOKES) &&
     (config->GetKind_Solver() != FEM_NAVIER_STOKES))
    SU2_MPI::Error("Navier Stokes equations must be selected for the Taylor Green Vortex",
                   CURRENT_FUNCTION);

  if((config->GetKind_FluidModel() != STANDARD_AIR) &&
     (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the Taylor Green Vortex",
                   CURRENT_FUNCTION);

  if(config->GetKind_ViscosityModel() != CONSTANT_VISCOSITY)
    SU2_MPI::Error("Constant viscosity must be selected for the Taylor Green Vortex",
                   CURRENT_FUNCTION);

  if(config->GetKind_ConductivityModel() != CONSTANT_PRANDTL)
    SU2_MPI::Error("Constant Prandtl number must be selected for the Taylor Green Vortex",
                   CURRENT_FUNCTION);
}

CTGVSolution::~CTGVSolution(void) { }

void CTGVSolution::GetSolution(const unsigned short val_nParams,
                               const su2double      *val_params,
                               const su2double      *val_coords,
                               const su2double      val_t,
                               su2double            *val_solution) {
  
  /* The initial conditions are set for the Taylor-Green vortex case, which
   is a DNS case that features vortex breakdown into turbulence. These
   particular settings are for the typical Re = 1600 case (M = 0.08) with
   an initial temperature of 300 K. Note that this condition works in both
   2D and 3D. */
  
  su2double val_coordsZ      = 0.0;
  if (nDim == 3) val_coordsZ = val_coords[2];
  
  /* Compute the primitive variables. */
  
  su2double rho =  tgvDensity;
  su2double u   =  tgvVelocity * (sin(val_coords[0]/tgvLength)*
                                  cos(val_coords[1]/tgvLength)*
                                  cos(val_coordsZ  /tgvLength));
  su2double v   = -tgvVelocity * (cos(val_coords[0]/tgvLength)*
                                  sin(val_coords[1]/tgvLength)*
                                  cos(val_coordsZ  /tgvLength));
  
  su2double factorA = cos(2.0*val_coordsZ/tgvLength) + 2.0;
  su2double factorB = (cos(2.0*val_coords[0]/tgvLength) +
                       cos(2.0*val_coords[1]/tgvLength));
  
  su2double p = (tgvPressure +
                 tgvDensity*(pow(tgvVelocity,2.0)/16.0)*factorA*factorB);
  
  /* Compute the conservative variables. Note that both 2D and 3D
   cases are treated correctly. */
  
  val_solution[0]      = rho;
  val_solution[1]      = rho*u;
  val_solution[2]      = rho*v;
  val_solution[3]      = 0.0;
  val_solution[nVar-1] = p*ovGm1 + 0.5*rho*(u*u + v*v);
}

bool CTGVSolution::ExactSolutionKnown(void) {return false;}

CIncTGVSolution::CIncTGVSolution(void) : CVerificationSolution() { }

CIncTGVSolution::CIncTGVSolution(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 CConfig*       config)
: CVerificationSolution(val_nDim, val_nVar, config) {
  
  /*--- Write a message that the solution is initialized for the
   Taylor-Green vortex test case. ---*/
  
  if (rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the incompressible Taylor-Green vortex case!!!" << endl;
    cout << endl << flush;
  }
  
  /*--- Store TGV specific parameters here. ---*/
  
  tgvLength    = 1.0;
  tgvVelocity  = 1.0;
  tgvDensity   = config->GetDensity_FreeStreamND();
  tgvViscosity = config->GetViscosity_FreeStreamND();
  
  /*--- We keep a copy of the freestream temperature just to be safe
   when we set the solution, even though this is an isothermal case. ---*/
  
  Temperature = config->GetTemperature_FreeStreamND();
  
  /*--- Perform some sanity and error checks for this solution here. ---*/
  
  if((config->GetUnsteady_Simulation() != TIME_STEPPING) &&
     (config->GetUnsteady_Simulation() != DT_STEPPING_1ST) &&
     (config->GetUnsteady_Simulation() != DT_STEPPING_2ND))
    SU2_MPI::Error("Unsteady mode must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_Regime() != INCOMPRESSIBLE)
    SU2_MPI::Error("Incompressible flow equations must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_Solver() != NAVIER_STOKES)
    SU2_MPI::Error("Navier Stokes equations must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_FluidModel() != CONSTANT_DENSITY)
    SU2_MPI::Error("Constant density fluid model must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_ViscosityModel() != CONSTANT_VISCOSITY)
    SU2_MPI::Error("Constant viscosity must be selected for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
  if(config->GetEnergy_Equation())
    SU2_MPI::Error("Energy equation must be disabled (isothermal) for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
  if(nDim != 2)
    SU2_MPI::Error("2D calculation required for the incompressible Taylor Green Vortex",
                   CURRENT_FUNCTION);
}

CIncTGVSolution::~CIncTGVSolution(void) { }

void CIncTGVSolution::GetBCState(const unsigned short val_nParams,
                                 const su2double      *val_params,
                                 const su2double      *val_coords,
                                 const su2double      val_t,
                                 su2double            *val_solution) {
  
  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_nParams, val_params, val_coords, val_t, val_solution);
}

void CIncTGVSolution::GetSolution(const unsigned short val_nParams,
                                  const su2double      *val_params,
                                  const su2double      *val_coords,
                                  const su2double      val_t,
                                  su2double            *val_solution) {
  
  /* The exact solution is set for the incompressible Taylor-Green
   vortex case. This is the classic solution from the original work
   of Taylor and Green for the specific 2D situation where the
   exact solution can be derived for an incompressible flow. */
  
  /* Store the termporal term more easily (Taylor expansion). */
  su2double F = 1.0 - 2.0*(tgvViscosity/tgvDensity)*val_t;
  
  /* Compute the primitive variables. */
  su2double u =  tgvVelocity * F * (sin(val_coords[0]/tgvLength)*
                                    cos(val_coords[1]/tgvLength));
  su2double v = -tgvVelocity * F * (cos(val_coords[0]/tgvLength)*
                                    sin(val_coords[1]/tgvLength));
  
  su2double B = (cos(2.0*val_coords[0]/tgvLength) +
                 cos(2.0*val_coords[1]/tgvLength));
  su2double p = -(tgvDensity/4.0)*B*F*F;
  
  /* Compute the conservative variables. Note that both 2D and 3D
   cases are treated correctly. */
  val_solution[0]      = p;
  val_solution[1]      = u;
  val_solution[2]      = v;
  val_solution[nVar-1] = Temperature;
}

CMMSNSUnitQuadSolution::CMMSNSUnitQuadSolution(void) : CVerificationSolution() { }

CMMSNSUnitQuadSolution::CMMSNSUnitQuadSolution(unsigned short val_nDim,
                                               unsigned short val_nVar,
                                               CConfig*       config)
  : CVerificationSolution(val_nDim, val_nVar, config) {
  
  /*--- Write a message that the solution is initialized for the manufactured
   solution for the Navier-Stokes equations on a unit quad. ---*/
  if (rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the manufactured solution " << endl;
    cout << "         of the Navier-Stokes equations on a unit quad!!!" << endl;
    cout << endl << flush;
  }

  /*--- Coefficients, needed to determine the solution. ---*/
  const su2double Prandtl = config->GetPrandtl_Lam();

  RGas         = config->GetGas_ConstantND();
  Gamma        = config->GetGamma();
  Viscosity    = config->GetMu_ConstantND();
  Conductivity = Viscosity*Gamma*RGas/(Prandtl*(Gamma-1.0));

  /*--- Constants, which describe this manufactured solution. This is a viscous
        solution on the unit quad, where the primitive variables vary as a
        combination of sine and cosine functions. The unit quad is probably not
        necessary, and an arbitrary domain should work as well. ---*/
  L       =      1.0;
  a_Px    =      1.0;
  a_Pxy   =      0.75;
  a_Py    =      1.25;
  a_rhox  =      0.75;
  a_rhoxy =      1.25;
  a_rhoy  =      1.0;
  a_ux    =      1.6666666667;
  a_uxy   =      0.6;
  a_uy    =      1.5;
  a_vx    =      1.5;
  a_vxy   =      0.9;
  a_vy    =      1.0;
  P_0     = 100000.0;
  P_x     = -30000.0;
  P_xy    = -25000.0;
  P_y     =  20000.0;
  rho_0   =      1.0;
  rho_x   =      0.1;
  rho_xy  =      0.08;
  rho_y   =      0.15;
  u_0     =     70.0;
  u_x     =      4.0;
  u_xy    =      7.0;
  u_y     =    -12.0;
  v_0     =     90.0;
  v_x     =    -20.0;
  v_xy    =    -11.0;
  v_y     =      4.0;

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if(config->GetUnsteady_Simulation() != STEADY)
    SU2_MPI::Error("Steady mode must be selected for the MMS NS Unit Quad case",
                   CURRENT_FUNCTION);

  if(config->GetKind_Regime() != COMPRESSIBLE)
    SU2_MPI::Error("Compressible flow equations must be selected for the MMS NS Unit Quad case",
                   CURRENT_FUNCTION);

  if((config->GetKind_Solver() != NAVIER_STOKES) &&
     (config->GetKind_Solver() != FEM_NAVIER_STOKES))
    SU2_MPI::Error("Navier Stokes equations must be selected for the MMS NS Unit Quad case",
                   CURRENT_FUNCTION);

  if((config->GetKind_FluidModel() != STANDARD_AIR) &&
     (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the MMS NS Unit Quad case",
                   CURRENT_FUNCTION);

  if(config->GetKind_ViscosityModel() != CONSTANT_VISCOSITY)
    SU2_MPI::Error("Constant viscosity must be selected for the MMS NS Unit Quad case",
                   CURRENT_FUNCTION);

  if(config->GetKind_ConductivityModel() != CONSTANT_PRANDTL)
    SU2_MPI::Error("Constant Prandtl number must be selected for the MMS NS Unit Quad case",
                   CURRENT_FUNCTION);
}

CMMSNSUnitQuadSolution::~CMMSNSUnitQuadSolution(void) { }

void CMMSNSUnitQuadSolution::GetBCState(const unsigned short val_nParams,
                                        const su2double      *val_params,
                                        const su2double      *val_coords,
                                        const su2double      val_t,
                                        su2double            *val_solution) {

  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_nParams, val_params, val_coords, val_t, val_solution);
}

void CMMSNSUnitQuadSolution::GetSolution(const unsigned short val_nParams,
                                         const su2double      *val_params,
                                         const su2double      *val_coords,
                                         const su2double      val_t,
                                         su2double            *val_solution) {

  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  /* Determine the solution for the density, velocity
     components and pressure. */
  const su2double LInv    = 1.0/L;
  const su2double PiLInv  = PI_NUMBER*LInv;
  const su2double PiL2Inv = PiLInv*LInv;

  const su2double rho = rho_0 + rho_x *sin(a_rhox *PiLInv*x)
                      +         rho_y *cos(a_rhoy *PiLInv*y)
                      +         rho_xy*cos(a_rhoxy*PiL2Inv*x*y);

  const su2double u = u_0 + u_x *sin(a_ux *PiLInv*x)
                    +       u_y *cos(a_uy *PiLInv*y)
                    +       u_xy*cos(a_uxy*PiL2Inv*x*y);

  const su2double v = v_0 + v_x *cos(a_vx *PiLInv*x)
                    +       v_y *sin(a_vy *PiLInv*y)
                    +       v_xy*cos(a_vxy*PiL2Inv*x*y);

  const su2double p = P_0 + P_x *cos(a_Px *PiLInv*x)
                    +       P_y *sin(a_Py *PiLInv*y)
                    +       P_xy*sin(a_Pxy*PiL2Inv*x*y);

  /* Compute the conservative variables from the primitive ones.
     Note that the implementation below is valid for both 2D and 3D. */
  val_solution[0]      = rho;
  val_solution[1]      = rho*u;
  val_solution[2]      = rho*v;
  val_solution[3]      = 0.0;
  val_solution[nDim+1] = p/(Gamma-1.0) + 0.5*rho*(u*u + v*v);
}

void CMMSNSUnitQuadSolution::GetMMSSourceTerm(const unsigned short val_nParams,
                                              const su2double      *val_params,
                                              const su2double      *val_coords,
                                              const su2double      val_t,
                                              su2double            *val_source) {

  /*--- The source code for the source terms is generated in Maple. ---*/
  const su2double Pi = PI_NUMBER;
  const su2double fourThird = 4.0/3.0;
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  const su2double t1 = rho_x * a_rhox;
  const su2double t2 = 0.1e1 / L;
  const su2double t3 = t2 * Pi;
  const su2double t4 = a_rhox * Pi;
  const su2double t5 = x * t2;
  const su2double t6 = t5 * t4;
  const su2double t7 = cos(t6);
  const su2double t10 = a_rhoxy * rho_xy;
  const su2double t11 = Pi * t10;
  const su2double t12 = L * L;
  const su2double t13 = 0.1e1 / t12;
  const su2double t14 = y * t13;
  const su2double t16 = x * t13;
  const su2double t17 = y * t16;
  const su2double t18 = t17 * a_rhoxy * Pi;
  const su2double t19 = sin(t18);
  const su2double t22 = t7 * t3 * t1 - t19 * t14 * t11;
  const su2double t24 = t5 * a_ux * Pi;
  const su2double t25 = sin(t24);
  const su2double t26 = t25 * u_x;
  const su2double t28 = y * t2;
  const su2double t29 = t28 * a_uy * Pi;
  const su2double t30 = cos(t29);
  const su2double t31 = t30 * u_y;
  const su2double t33 = t17 * a_uxy * Pi;
  const su2double t34 = cos(t33);
  const su2double t36 = t34 * u_xy + t26 + t31 + u_0;
  const su2double t37 = t36 * t22;
  const su2double t38 = sin(t6);
  const su2double t39 = t38 * rho_x;
  const su2double t40 = a_rhoy * Pi;
  const su2double t41 = t28 * t40;
  const su2double t42 = cos(t41);
  const su2double t43 = t42 * rho_y;
  const su2double t44 = cos(t18);
  const su2double t46 = t44 * rho_xy + rho_0 + t39 + t43;
  const su2double t48 = cos(t24);
  const su2double t52 = u_xy * a_uxy * Pi;
  const su2double t53 = sin(t33);
  const su2double t56 = t48 * t3 * u_x * a_ux - t53 * t14 * t52;
  const su2double t57 = t56 * t46;
  const su2double t58 = rho_y * a_rhoy;
  const su2double t59 = sin(t41);
  const su2double t60 = t59 * t3;
  const su2double t64 = -t19 * t16 * t11 - t60 * t58;
  const su2double t66 = t5 * a_vx * Pi;
  const su2double t67 = cos(t66);
  const su2double t68 = t67 * v_x;
  const su2double t70 = t28 * a_vy * Pi;
  const su2double t71 = sin(t70);
  const su2double t72 = t71 * v_y;
  const su2double t74 = t17 * a_vxy * Pi;
  const su2double t75 = cos(t74);
  const su2double t77 = t75 * v_xy + t68 + t72 + v_0;
  const su2double t80 = cos(t70);
  const su2double t84 = v_xy * a_vxy * Pi;
  const su2double t85 = sin(t74);
  const su2double t88 = t80 * t3 * v_y * a_vy - t85 * t16 * t84;
  const su2double t91 = t36 * t36;
  const su2double t93 = t36 * t46;
  const su2double t98 = t5 * a_Px * Pi;
  const su2double t99 = sin(t98);
  const su2double t101 = t99 * t3 * P_x * a_Px;
  const su2double t103 = P_xy * a_Pxy * Pi;
  const su2double t104 = a_Pxy * Pi;
  const su2double t105 = t17 * t104;
  const su2double t106 = cos(t105);
  const su2double t107 = t106 * t14;
  const su2double t108 = t107 * t103;
  const su2double t109 = a_uxy * a_uxy;
  const su2double t111 = t13 * Pi * t109;
  const su2double t112 = y * y;
  const su2double t116 = a_vxy * a_vxy;
  const su2double t118 = t13 * Pi * t116;
  const su2double t120 = v_xy * x;
  const su2double t122 = t120 * t75 * y * t118;
  const su2double t124 = a_vxy * t85;
  const su2double t125 = v_xy * t124;
  const su2double t127 = a_ux * a_ux;
  const su2double t131 = (-u_xy * t34 * t112 * t111 + t122 / 0.2e1 + t125 / 0.2e1 - t26 * Pi * t127) * Pi;
  const su2double t132 = t13 * Viscosity;
  const su2double t137 = u_y * a_uy;
  const su2double t138 = sin(t29);
  const su2double t143 = -t138 * t3 * t137 - t53 * t16 * t52;
  const su2double t147 = a_uy * a_uy;
  const su2double t150 = x * x;
  const su2double t161 = (-t13 * (u_xy * t34 * t150 * t111 + t31 * Pi * t147) * Pi - t13 * (t122 + t125) * Pi) * Viscosity;
  const su2double t165 = v_x * a_vx;
  const su2double t166 = sin(t66);
  const su2double t171 = -t85 * t14 * t84 - t166 * t3 * t165;
  const su2double t174 = u_xy * x;
  const su2double t176 = t174 * t34 * y * t111;
  const su2double t177 = a_uxy * t53;
  const su2double t178 = u_xy * t177;
  const su2double t182 = a_vx * a_vx;
  const su2double t192 = (-t13 * (t176 + t178) * Pi - t13 * (v_xy * t75 * t112 * t118 + t68 * Pi * t182) * Pi) * Viscosity;
  const su2double t193 = t77 * t77;
  const su2double t200 = t28 * a_Py * Pi;
  const su2double t201 = cos(t200);
  const su2double t203 = t201 * t3 * P_y * a_Py;
  const su2double t204 = t106 * t16;
  const su2double t205 = t204 * t103;
  const su2double t210 = a_vy * a_vy;
  const su2double t215 = (0.2e1 * v_xy * t75 * t150 * t118 + 0.2e1 * t72 * Pi * t210 - t176 - t178) * Pi;
  const su2double t219 = -t101 + t108;
  const su2double t221 = 1 / (Gamma - 1);
  const su2double t223 = t91 + t193;
  const su2double t233 = cos(t98);
  const su2double t234 = t233 * P_x;
  const su2double t235 = sin(t200);
  const su2double t236 = t235 * P_y;
  const su2double t237 = sin(t105);
  const su2double t238 = t237 * P_xy;
  const su2double t239 = P_0 + t234 + t236 + t238;
  const su2double t243 =  t221 * t239 + t223 * t46 / 0.2e1 + P_0 + t234 + t236 + t238;
  const su2double t245 = a_Px * a_Px;
  const su2double t246 = Pi * t245;
  const su2double t248 = a_Pxy * a_Pxy;
  const su2double t250 = P_xy * t112 * t248;
  const su2double t251 = t13 * Pi;
  const su2double t252 = t237 * t251;
  const su2double t258 = P_x * t99;
  const su2double t260 = a_Pxy * y;
  const su2double t261 = t106 * P_xy;
  const su2double t264 = (t258 * L * a_Px - t261 * t260) * rho_xy;
  const su2double t267 = t19 * y * t251;
  const su2double t272 = a_rhoxy * a_rhoxy;
  const su2double t273 = t272 * rho_xy;
  const su2double t277 = t44 * t13 * Pi * t239;
  const su2double t285 = t39 + t43 + rho_0;
  const su2double t288 = t237 * t13 * Pi * t285;
  const su2double t290 = a_rhox * a_rhox;
  const su2double t291 = t290 * rho_x;
  const su2double t296 = P_xy * t7;
  const su2double t306 = t234 + t236 + P_0;
  const su2double t308 = t38 * t3;
  const su2double t311 = t39 + rho_0;
  const su2double t318 = t44 * t44;
  const su2double t319 = rho_xy * rho_xy;
  const su2double t321 = t285 * rho_xy;
  const su2double t324 = t42 * t42;
  const su2double t325 = rho_y * rho_y;
  const su2double t327 = t311 * rho_y;
  const su2double t330 = t7 * t7;
  const su2double t331 = rho_x * rho_x;
  const su2double t336 = rho_0 * rho_0;
  const su2double t337 = 0.2e1 * rho_x * rho_0 * t38 + t319 * t318 + 0.2e1 * t44 * t321 + t325 * t324 + 0.2e1 * t42 * t327 - t331 * t330 + t331 + t336;
  const su2double t338 = 0.1e1 / t337;
  const su2double t340 = 0.1e1 / RGas;
  const su2double t341 = t340 * t13;
  const su2double t349 = t106 * t285 * P_xy;
  const su2double t365 = t337 * t337;
  const su2double t366 = 0.1e1 / t365;
  const su2double t369 = a_rhoxy * t319 * t44;
  const su2double t376 = a_rhoxy * t321;
  const su2double t398 = u_xy * y * t177;
  const su2double t399 = t120 * t124;
  const su2double t402 = a_ux * u_x * t48;
  const su2double t404 = a_vy * v_y * t80;
  const su2double t428 = (-t13 * (t137 * t138 * L + t174 * t177) * Pi - t13 * (t165 * t166 * L + v_xy * y * t124) * Pi) * Viscosity;
  const su2double t430 = t203 + t205;
  const su2double t442 = a_Py * a_Py;
  const su2double t446 = P_xy * t150 * t248;
  const su2double t452 = P_y * t201;
  const su2double t454 = a_Pxy * x;
  const su2double t457 = (t452 * L * a_Py + t261 * t454) * rho_xy;
  const su2double t460 = t19 * x * t251;
  const su2double t474 = a_rhoy * a_rhoy;
  const su2double t475 = t474 * rho_y;
  const su2double t480 = P_xy * t59;
  const su2double t486 = t235 * t2;
  const su2double t564 = t36 * ( t221 * t219 + t223 * t22 / 0.2e1 + (t171 * t77 + t56 * t36) * t46 - t101 + t108) + t56 * t243
                       + Conductivity * t341 * t338 * (t44 * (t234 * t246 + t252 * t250) * rho_xy - t267 * a_rhoxy * t264
                       - t19 * t219 * y * t10 - t277 * t112 * t273 - t106 * t7 * t2 * t4 * P_xy * rho_x * t260 + t288 * t250
                       + (t42 * P_x * t233 * t2 * Pi * rho_y * t245 - t238 * t38 * t2 * Pi * t291
                       + t107 * t104 * t296 * t1 + t311 * t234 * t2 * t246 - t308 * t306 * t291) * L) * Pi
                       - 0.2e1 * (t44 * t7 * t3 * rho_xy * rho_x * a_rhox + t42 * t7 * t3 * rho_y * rho_x * a_rhox
                       + rho_x * rho_0 * t7 * t2 * t4 + t308 * a_rhox * t331 * t7 - t267 * t369 - t267 * t376) * Conductivity * t341 * t366
                       * (t44 * t264 - t19 * t239 * y * t10 - t349 * t260 + (t311 * P_x * t99 * a_Px + t42 * t258 * a_Px * rho_y
                       + t237 * t296 * t1 + t7 * t306 * t1) * L) * Pi - fourThird * t36 * t132 * t131
                       - fourThird * t56 * t132 * (-t398 + t399 / 0.2e1 + (t402 - t404 / 0.2e1) * L) * Pi
                       - t77 * t192 - t171 * t428 + t77 * ( t221 * t430 + t223 * t64 / 0.2e1 + (t143 * t36 + t88 * t77) * t46 + t203 + t205)
                       + t88 * t243 - Conductivity * t340 * t338 * t13 * (t44 * (-t236 * Pi * t442 - t252 * t446) * rho_xy
                       - t460 * a_rhoxy * t457 + t19 * t430 * x * t10 + t277 * t150 * t273 - t106 * t59 * t2 * t40 * P_xy * rho_y * t454
                       - t288 * t446 + (-P_y * t38 * t486 * Pi * rho_x * t442 - t42 * P_y * t486 * Pi * rho_y * t442
                       - P_y * t486 * Pi * rho_0 * t442 + t238 * t42 * t2 * Pi * t475 + t204 * t104 * t480 * t58
                       + t42 * t3 * t306 * t475) * L) * Pi + 0.2e1 * (-t44 * t59 * t3 * rho_xy * rho_y * a_rhoy
                       - t60 * a_rhoy * t325 * t42 - t60 * a_rhoy * t327 - t460 * t369 - t460 * t376) * Conductivity * t340 * t366 * t13
                       * (t44 * t457 + t19 * t239 * x * t10 + t349 * t454 + (P_y * t38 * t201 * a_Py * rho_x + t42 * t452 * a_Py * rho_y
                       + t452 * a_Py * rho_0 + t237 * t480 * t58 + t59 * t306 * t58) * L) * Pi - t36 * t161
                       - t143 * t428 + 0.2e1 / 0.3e1 * t77 * t132 * t215 + 0.2e1 / 0.3e1 * t88 * t132 * (-t398 + 0.2e1 * t399
                       + (t402 - 0.2e1 * t404) * L) * Pi;

  val_source[0]      = t88 * t46 + t77 * t64 + t37 + t57;
  val_source[1]      = t91 * t22 + 0.2e1 * t56 * t93 - t101 + t108 - fourThird * t132 * t131 + t77 * t36 * t64 + t77 * t143 * t46 + t88 * t93 - t161;
  val_source[2]      = t77 * t37 + t77 * t57 + t171 * t93 - t192 + t193 * t64 + 0.2e1 * t88 * t77 * t46 + t203 + t205 + 0.2e1 / 0.3e1 * t132 * t215;
  val_source[3]      = 0.0;
  val_source[nDim+1] = t564;    
}

bool CMMSNSUnitQuadSolution::IsManufacturedSolution(void) {
  return true;
}

CMMSNSTwoHalfSpheresSolution::CMMSNSTwoHalfSpheresSolution(void) : CVerificationSolution() { }

CMMSNSTwoHalfSpheresSolution::CMMSNSTwoHalfSpheresSolution(unsigned short val_nDim,
                                                           unsigned short val_nVar,
                                                           CConfig*       config)
  : CVerificationSolution(val_nDim, val_nVar, config) {
  
  /*--- Write a message that the solution is initialized for the manufactured
   solution for the Navier-Stokes equations between two half spheres. ---*/
  if (rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the manufactured solution " << endl;
    cout << "         of the Navier-Stokes equations between two half spheres!!!" << endl;
    cout << endl << flush;
  }

  /*--- Coefficients, needed to determine the solution. ---*/
  const su2double Prandtl = config->GetPrandtl_Lam();

  RGas         = config->GetGas_Constant();
  Gamma        = config->GetGamma();
  Viscosity    = config->GetMu_Constant();
  Conductivity = Viscosity*Gamma*RGas/(Prandtl*(Gamma-1.0));

  /*--- Initialize TWall to the default value of 300 K (in case the outer wall
        is not modelled as an isothermal wall) and try to retrieve the wall
        temperature from the boundary conditions. ---*/
  TWall = 300.0;
  for(unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if(config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) {
      const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
      TWall = config->GetIsothermal_Temperature(Marker_Tag);
    }
  }

  /*--- Get the reference values for pressure, density and velocity. ---*/
  Pressure_Ref = config->GetPressure_Ref();
  Density_Ref  = config->GetDensity_Ref();
  Velocity_Ref = config->GetVelocity_Ref();

  /*--- The constants for the density and velocities. ---*/
  rho_0 =   1.25;
  u_0   = 135.78;
  v_0   = -67.61;
  w_0   =  82.75;

  /*--- The constants for the temperature solution. ---*/
  a_T1 =  1.05;
  a_T2 = -0.85;

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if(nDim != 3)
    SU2_MPI::Error("Grid must be 3D for the MMS NS Two Half Spheres case",
                   CURRENT_FUNCTION);

  if(config->GetUnsteady_Simulation() != STEADY)
    SU2_MPI::Error("Steady mode must be selected for the MMS NS Two Half Spheres case",
                   CURRENT_FUNCTION);

  if(config->GetKind_Regime() != COMPRESSIBLE)
    SU2_MPI::Error("Compressible flow equations must be selected for the MMS NS Two Half Spheres case",
                   CURRENT_FUNCTION);

  if((config->GetKind_Solver() != NAVIER_STOKES) &&
     (config->GetKind_Solver() != FEM_NAVIER_STOKES))
    SU2_MPI::Error("Navier Stokes equations must be selected for the MMS NS Two Half Spheres case",
                   CURRENT_FUNCTION);

  if((config->GetKind_FluidModel() != STANDARD_AIR) &&
     (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the MMS NS Two Half Spheres case",
                   CURRENT_FUNCTION);

  if(config->GetKind_ViscosityModel() != CONSTANT_VISCOSITY)
    SU2_MPI::Error("Sutherland must be selected for viscosity for the MMS NS Two Half Spheres case",
                   CURRENT_FUNCTION);

  if(config->GetKind_ConductivityModel() != CONSTANT_PRANDTL)
    SU2_MPI::Error("Constant Prandtl number must be selected for the MMS NS Two Half Spheres case",
                   CURRENT_FUNCTION);
};

CMMSNSTwoHalfSpheresSolution::~CMMSNSTwoHalfSpheresSolution(void) { }

void CMMSNSTwoHalfSpheresSolution::GetBCState(const unsigned short val_nParams,
                                              const su2double      *val_params,
                                              const su2double      *val_coords,
                                              const su2double      val_t,
                                              su2double            *val_solution) {

  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_nParams, val_params, val_coords, val_t, val_solution);
}

void CMMSNSTwoHalfSpheresSolution::GetSolution(const unsigned short val_nParams,
                                               const su2double      *val_params,
                                               const su2double      *val_coords,
                                               const su2double      val_t,
                                               su2double            *val_solution) {

  /* Easier storage of the x-, y- and z-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  const su2double z = val_coords[2];

  /* Determine the dimensional solution for the temperature. */
  const su2double Pi   = PI_NUMBER;
  const su2double r    = sqrt(x*x + y*y + z*z);
  const su2double fact = (r-1.0)*(r-1.0)/(a_T1 + a_T2);

  su2double T = 0.25*TWall*(3.0 + fact*(a_T1*cos(Pi*(r-2.0))
              +                         a_T2*cos(Pi*(r-2.0)*2.0)));

  /* Determine the dimensional solution for the velocities. */
  su2double u = u_0*(r-1.0)*(2.0-r)*4.0;
  su2double v = v_0*(r-1.0)*(2.0-r)*4.0;
  su2double w = w_0*(r-1.0)*(2.0-r)*4.0;

  /* Compute the pressure from the density and temperature. */
  su2double rho = rho_0;
  su2double p   = rho*RGas*T;

  /* Determine the non-dimensional solution. */
  rho /= Density_Ref;
  p   /= Pressure_Ref;
  u   /= Velocity_Ref;
  v   /= Velocity_Ref;
  w   /= Velocity_Ref;

  /* Determine the non-dimensional conserved variables. */
  val_solution[0] = rho;
  val_solution[1] = rho*u;
  val_solution[2] = rho*v;
  val_solution[3] = rho*w;
  val_solution[4] = p/(Gamma-1.0) + 0.5*rho*(u*u + v*v + w*w);
}

void CMMSNSTwoHalfSpheresSolution::GetMMSSourceTerm(const unsigned short val_nParams,
                                                    const su2double      *val_params,
                                                    const su2double      *val_coords,
                                                    const su2double      val_t,
                                                    su2double            *val_source) {

  /*--- Abbreviate Pi and the coordinates. ---*/
  const su2double Pi = PI_NUMBER;
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  const su2double z = val_coords[2];

  /*--- The source code for the source terms is generated in Maple. ---*/
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
  const su2double t139 = 0.32e2 * t34 * t12 * t30 - 0.32e2 * t10 * t37 * t30 + t71 + 0.16e2 / 0.3e1 * t78 * x * t73
                       + 0.16e2 / 0.3e1 * t85 * t84 - 0.8e1 / 0.3e1 * x * t90 * t84 + 0.32e2 * t17 * t95 * t94
                       - 0.32e2 * t17 * t100 * t99 + 0.8e1 * t108 * t104 + 0.16e2 * t85 * t112 - 0.8e1 * y * t115 * t112
                       + 0.32e2 * t120 * t94 - 0.32e2 * t124 * t99 + 0.8e1 * t131 * t130 - 0.8e1 * t136 * t134;
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
  const su2double t202 = 0.32e2 * t8 * t95 * t94 - 0.32e2 * t8 * t100 * t99 + 0.8e1 * t108 * t146
                       + 0.80e2 / 0.3e1 * t7 * v_0 * t112 - 0.8e1 * x * t115 * t112 + 0.32e2 * t158 * t12 * t156
                       - 0.32e2 * t18 * t37 * t156 + t177 + 0.16e2 / 0.3e1 * y * t180 * t135
                       - 0.16e2 / 0.3e1 * y * t179 * t73 + 0.32e2 * t120 * t187 - 0.32e2 * t124 * t190
                       + 0.8e1 * t131 * t196 - 0.8e1 * t136 * t199;
  const su2double t209 = t111 * w_0;
  const su2double t231 = w_0 * w_0;
  const su2double t232 = t231 * rho_0;
  const su2double t234 = z * t33;
  const su2double t244 = t58 * t23 * t57;
  const su2double t252 = (z * t52 * t44 / 0.2e1 + (-0.2e1 * t62 * t23 * t61 - t244) * t56 / 0.4e1) * TWall;
  const su2double t253 = t252 * t41;
  const su2double t255 = t74 + t76 - 0.2e1 * t77;
  const su2double t256 = t111 * t255;
  const su2double t263 = 0.32e2 * t8 * t119 * t94 - 0.32e2 * t8 * t123 * t99 + 0.80e2 / 0.3e1 * t7 * Viscosity * t209
                       + 0.8e1 * x * Viscosity * t130 - 0.8e1 * x * t135 * t134 + 0.32e2 * t17 * t119 * t187
                       - 0.32e2 * t17 * t123 * t190 + 0.8e1 * y * Viscosity * t196 - 0.8e1 * y * t135 * t199
                       + 0.32e2 * t234 * t12 * t232 - 0.32e2 * t24 * t37 * t232 + t253
                       + 0.16e2 / 0.3e1 * z * t256 * t135 - 0.16e2 / 0.3e1 * z * t255 * t73;
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
  const su2double t329 = -0.32e2 / 0.3e1 * t266 * t265 * t90 * t84 - 0.32e2 * t266 * t271 * t270
                       - 0.32e2 * t266 * t276 * t275 - 0.32e2 * t280 * t265 * t270 + 0.64e2 / 0.3e1 * t285 * t180 * t135
                       - 0.32e2 * t280 * t276 * t288 - 0.32e2 * t292 * t265 * t275 - 0.32e2 * t292 * t271 * t288
                       + 0.64e2 / 0.3e1 * t300 * t256 * t135 + 0.32e2 * t305 * v_0 * t107 * t303 + 0.32e2 * t266 * t310 * t309
                       + 0.32e2 * t316 * x * t130 + 0.32e2 * t266 * t320 * t319 + 0.32e2 * t305 * u_0 * t107 * Viscosity * t104;
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
  const su2double t426 = 0.32e2 * t280 * t330 * t309 - 0.64e2 / 0.3e1 * t285 * t334 - 0.64e2 / 0.3e1 * t280 * t337 * t334
                       + 0.32e2 * t316 * y * t196 + 0.32e2 * t280 * t320 * t344 + 0.32e2 * t305 * Viscosity * u_0 * z * t130
                       + 0.32e2 * t292 * t330 * t319 + 0.32e2 * t305 * Viscosity * v_0 * z * t196 + 0.32e2 * t292 * t310 * t344
                       - 0.64e2 / 0.3e1 * t300 * t364 - 0.64e2 / 0.3e1 * t292 * t209 * t364
                       + 0.4e1 * t400 * TWall * x * t12 * t7 * (t377 * x * t371 - t12 * a_T2 * x * t381 + Pi * t8 * t385
                       + 0.2e1 * a_T2 * t58 * x * t389 + t60 / 0.2e1) - 0.4e1 * t400 * t2 * TWall
                       * t412 + 0.64e2 / 0.3e1 * t305 * u_0 * t78 * t303 - 0.32e2 / 0.3e1 * t422 * t330 * t421;
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
  const su2double t527 = -0.32e2 * t422 * t310 * t309 - 0.32e2 * t422 * t320 * t319
                       + 0.4e1 * t400 * y * TWall * t12 * t7 * (t377 * y * t371 - t12 * a_T2 * y * t381 + Pi * t17 * t385
                       + 0.2e1 * a_T2 * t58 * y * t389 + t168 / 0.2e1) - 0.4e1 * t400 * TWall * t3 * t412
                       - 0.32e2 * t457 * t330 * t309 + 0.64e2 / 0.3e1 * t457 * t337 * t334 - 0.32e2 * t457 * t320 * t344
                       + 0.4e1 * t483 * t7 * (t377 * z * t371 - t12 * a_T2 * z * t381 + Pi * t23 * t385
                       + 0.2e1 * a_T2 * t58 * z * t389 + t244 / 0.2e1) * z - 0.4e1 * t483 * t89 * t486
                       - 0.32e2 * t490 * t330 * t319 - 0.32e2 * t490 * t310 * t344 + 0.64e2 / 0.3e1 * t490 * t209 * t364
                       + 0.32e2 / 0.3e1 * t266 * t330 * t421 + 0.4e1 * t10 * t522 - 0.4e1 * t14 * t522;
  const su2double t528 = v_0 * t521;
  const su2double t533 = w_0 * t521;
  const su2double t541 = Conductivity * t399 * TWall;
  const su2double t545 = t9 * t13;
  const su2double t558 = t72 * t410;
  const su2double t570 = t12 * t29;
  const su2double t573 = t12 * t155;
  const su2double t576 = t12 * t231;
  const su2double t616 = 0.4e1 * t18 * t528 - 0.4e1 * t20 * t528 + 0.4e1 * t24 * t533 - 0.4e1 * t26 * t533
                       + 0.12e2 * t541 * t12 * t7 * t410 + 0.320e3 / 0.3e1 * t545 * t155 * t112
                       + 0.320e3 / 0.3e1 * t545 * Viscosity * t111 * t231 + 0.64e2 * t545 * t29 * t112
                       + 0.64e2 / 0.3e1 * t545 * t29 * t84 + 0.4e1 * t541 * t2 * t558 + 0.4e1 * t541 * t3 * t558 
                       + 0.4e1 * t541 * t72 * t486 + 0.4e1 * t305 * u_0 * (t508 * t70 * t41
                       + 0.16e2 * (-t10 * t511 - t10 * t513 - t10 * t515 + t34 * t570 + t34 * t573 + t34 * t576) * rho_0 + t71)
                       + 0.4e1 * t305 * v_0 * (t508 * t176 * t41 + 0.16e2 * (t158 * t570 + t158 * t573 + t158 * t576
                       - t18 * t511 - t18 * t513 - t18 * t515) * rho_0 + t177) + 0.4e1 * t305 * w_0 * (t508 * t252 * t41
                       + 0.16e2 * (t234 * t570 + t234 * t573 + t234 * t576 - t24 * t511 - t24 * t513 - t24 * t515) * rho_0 + t253);

  /*--- Set the source term. Note the scaling for the correct non-dimensionalization. ---*/
  val_source[0] = 0.4e1 * t10 * t1 - 0.4e1 * t14 * t1 + 0.4e1 * t18 * t16 - 0.4e1 * t20 * t16 + 0.4e1 * t24 * t22 - 0.4e1 * t26 * t22;
  val_source[1] = t139;
  val_source[2] = t202;
  val_source[3] = t263;
  val_source[4] = t329 + t426 + t527 + t616;

  val_source[0] /= Density_Ref*Velocity_Ref;
  val_source[1] /= Pressure_Ref;
  val_source[2] /= Pressure_Ref;
  val_source[3] /= Pressure_Ref;
  val_source[4] /= Velocity_Ref*Pressure_Ref;
}

bool CMMSNSTwoHalfSpheresSolution::IsManufacturedSolution(void) {
  return true;
}

CMMSIncEulerSolution::CMMSIncEulerSolution(void) : CVerificationSolution() { }

CMMSIncEulerSolution::CMMSIncEulerSolution(unsigned short val_nDim,
                                           unsigned short val_nVar,
                                           CConfig*       config)
: CVerificationSolution(val_nDim, val_nVar, config) {
  
  /*--- Write a message that the solution is initialized for the manufactured
   solution for the incompressible Navier-Stokes equations. ---*/
  if (rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the manufactured solution " << endl;
    cout << "         of the incompressible Euler equations!!!" << endl;
    cout << endl << flush;
  }
  
  /*--- Coefficients, needed to determine the solution. ---*/
  Density     = config->GetDensity_FreeStreamND();
  Temperature = config->GetTemperature_FreeStreamND();
  
  /*--- Constants, which describe this manufactured solution. This is a
   solution where the primitive variables vary as a combination
   of sine and cosine functions. The solution is from Salari K, and
   Knupp P, "Code verification by the method of manufactured solutions,"
   SAND 2000-1444, Sandia National Laboratories, Albuquerque, NM, 2000. ---*/
  
  P_0     =   1.0;
  u_0     =   1.0;
  v_0     =   1.0;
  epsilon = 0.001;
  
  /*--- Perform some sanity and error checks for this solution here. ---*/
  if(config->GetUnsteady_Simulation() != STEADY)
    SU2_MPI::Error("Steady mode must be selected for the MMS incompressible Euler case",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_Regime() != INCOMPRESSIBLE)
    SU2_MPI::Error("Incompressible flow equations must be selected for the MMS incompressible Euler case",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_Solver() != EULER)
    SU2_MPI::Error("Euler equations must be selected for the MMS incompressible Euler case",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_FluidModel() != CONSTANT_DENSITY)
    SU2_MPI::Error("Constant density fluid model must be selected for the MMS incompressible Euler case",
                   CURRENT_FUNCTION);
  
  if(config->GetEnergy_Equation())
    SU2_MPI::Error("Energy equation must be disabled (isothermal) for the MMS incompressible Euler case",
                   CURRENT_FUNCTION);
}

CMMSIncEulerSolution::~CMMSIncEulerSolution(void) { }

void CMMSIncEulerSolution::GetBCState(const unsigned short val_nParams,
                                      const su2double      *val_params,
                                      const su2double      *val_coords,
                                      const su2double      val_t,
                                      su2double            *val_solution) {
  
  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_nParams, val_params, val_coords, val_t, val_solution);
}

void CMMSIncEulerSolution::GetSolution(const unsigned short val_nParams,
                                       const su2double      *val_params,
                                       const su2double      *val_coords,
                                       const su2double      val_t,
                                       su2double            *val_solution) {
  
  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  
  /* Compute the primitives from the defined solution. */
  const su2double u = u_0*(sin(x*x + y*y) + epsilon);
  const su2double v = v_0*(cos(x*x + y*y) + epsilon);
  const su2double p = P_0*(sin(x*x + y*y) +     2.0);
  
  /* For the incompressible solver, we return the primitive variables
   directly, as they are used for the working variables in the solver.
   Note that the implementation below is valid for both 2D and 3D. */
  val_solution[0]      = p;
  val_solution[1]      = u;
  val_solution[2]      = v;
  val_solution[3]      = 0.0;
  val_solution[nVar-1] = Temperature;
  
}

void CMMSIncEulerSolution::GetMMSSourceTerm(const unsigned short val_nParams,
                                            const su2double      *val_params,
                                            const su2double      *val_coords,
                                            const su2double      val_t,
                                            su2double            *val_source) {
  
  /*--- Easier storage of the x- and y-coordinates. ---*/
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  
  /*--- The expressions for the source terms are generated
   automatically by the sympy package in python.---*/
  val_source[0] = 2*Density*(u_0*x*cos(pow(x, 2) + pow(y, 2)) - v_0*y*sin(pow(x, 2) + pow(y, 2)));
  val_source[1] = 4*Density*pow(u_0, 2)*x*(epsilon + sin(pow(x, 2) + pow(y, 2)))*cos(pow(x, 2) + pow(y, 2)) - 2*Density*u_0*v_0*y*(epsilon + sin(pow(x, 2) + pow(y, 2)))*sin(pow(x, 2) + pow(y, 2)) + 2*Density*u_0*v_0*y*(epsilon + cos(pow(x, 2) + pow(y, 2)))*cos(pow(x, 2) + pow(y, 2)) + 2*P_0*x*cos(pow(x, 2) + pow(y, 2));
  val_source[2] = -2*Density*u_0*v_0*x*(epsilon + sin(pow(x, 2) + pow(y, 2)))*sin(pow(x, 2) + pow(y, 2)) + 2*Density*u_0*v_0*x*(epsilon + cos(pow(x, 2) + pow(y, 2)))*cos(pow(x, 2) + pow(y, 2)) - 4*Density*pow(v_0, 2)*y*(epsilon + cos(pow(x, 2) + pow(y, 2)))*sin(pow(x, 2) + pow(y, 2)) + 2*P_0*y*cos(pow(x, 2) + pow(y, 2));
  val_source[3]      = 0.0;
  val_source[nVar-1] = 0.0;
  
}

bool CMMSIncEulerSolution::IsManufacturedSolution(void) {
  return true;
}

CMMSIncNSSolution::CMMSIncNSSolution(void) : CVerificationSolution() { }

CMMSIncNSSolution::CMMSIncNSSolution(unsigned short val_nDim,
                                     unsigned short val_nVar,
                                     CConfig*       config)
: CVerificationSolution(val_nDim, val_nVar, config) {
  
  /*--- Write a message that the solution is initialized for the manufactured
   solution for the incompressible Navier-Stokes equations. ---*/
  if (rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the manufactured solution " << endl;
    cout << "         of the incompressible Navier-Stokes equations!!!" << endl;
    cout << endl << flush;
  }
  
  /*--- Coefficients, needed to determine the solution. ---*/
  Viscosity   = config->GetViscosity_FreeStreamND();
  Density     = config->GetDensity_FreeStreamND();
  Temperature = config->GetTemperature_FreeStreamND();
  
  /*--- Constants, which describe this manufactured solution. This is a
   viscous solution where the primitive variables vary as a combination
   of sine and cosine functions. The solution is from Salari K, and
   Knupp P, "Code verification by the method of manufactured solutions,"
   SAND 2000-1444, Sandia National Laboratories, Albuquerque, NM, 2000. ---*/
  
  P_0     =   1.0;
  u_0     =   1.0;
  v_0     =   1.0;
  epsilon = 0.001;
  
  /*--- Perform some sanity and error checks for this solution here. ---*/
  if(config->GetUnsteady_Simulation() != STEADY)
    SU2_MPI::Error("Steady mode must be selected for the MMS incompressible NS case",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_Regime() != INCOMPRESSIBLE)
    SU2_MPI::Error("Incompressible flow equations must be selected for the MMS incompressible NS case",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_Solver() != NAVIER_STOKES)
    SU2_MPI::Error("Navier Stokes equations must be selected for the MMS incompressible NS case",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_FluidModel() != CONSTANT_DENSITY)
    SU2_MPI::Error("Constant density fluid model must be selected for the MMS incompressible NS case",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_ViscosityModel() != CONSTANT_VISCOSITY)
    SU2_MPI::Error("Constant viscosity must be selected for the MMS incompressible NS case",
                   CURRENT_FUNCTION);
  
  if(config->GetEnergy_Equation())
    SU2_MPI::Error("Energy equation must be disabled (isothermal) for the MMS incompressible NS case",
                   CURRENT_FUNCTION);
}

CMMSIncNSSolution::~CMMSIncNSSolution(void) { }

void CMMSIncNSSolution::GetBCState(const unsigned short val_nParams,
                                   const su2double      *val_params,
                                   const su2double      *val_coords,
                                   const su2double      val_t,
                                   su2double            *val_solution) {
  
  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_nParams, val_params, val_coords, val_t, val_solution);
}

void CMMSIncNSSolution::GetSolution(const unsigned short val_nParams,
                                    const su2double      *val_params,
                                    const su2double      *val_coords,
                                    const su2double      val_t,
                                    su2double            *val_solution) {
  
  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  
  /* Compute the primitives from the defined solution. */
  const su2double u = u_0*(sin(x*x + y*y) + epsilon);
  const su2double v = v_0*(cos(x*x + y*y) + epsilon);
  const su2double p = P_0*(sin(x*x + y*y) +     2.0);
  
  /* For the incompressible solver, we return the primitive variables
   directly, as they are used for the working variables in the solver.
   Note that the implementation below is valid for both 2D and 3D. */
  val_solution[0]      = p;
  val_solution[1]      = u;
  val_solution[2]      = v;
  val_solution[3]      = 0.0;
  val_solution[nVar-1] = Temperature;
  
}

void CMMSIncNSSolution::GetMMSSourceTerm(const unsigned short val_nParams,
                                         const su2double      *val_params,
                                         const su2double      *val_coords,
                                         const su2double      val_t,
                                         su2double            *val_source) {
  
  /*--- Easier storage of the x- and y-coordinates. ---*/
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  
  /*--- The expressions for the source terms are generated
   automatically by the sympy package in python.---*/
  val_source[0] = 2*Density*(u_0*x*cos(pow(x, 2) + pow(y, 2)) - v_0*y*sin(pow(x, 2) + pow(y, 2)));
  val_source[1] = 4*Density*pow(u_0, 2)*x*(epsilon + sin(pow(x, 2) + pow(y, 2)))*cos(pow(x, 2) + pow(y, 2)) - 2*Density*u_0*v_0*y*(epsilon + sin(pow(x, 2) + pow(y, 2)))*sin(pow(x, 2) + pow(y, 2)) + 2*Density*u_0*v_0*y*(epsilon + cos(pow(x, 2) + pow(y, 2)))*cos(pow(x, 2) + pow(y, 2)) + 2*P_0*x*cos(pow(x, 2) + pow(y, 2)) - 0.666666666666667*Viscosity*(-8.0*u_0*pow(x, 2)*sin(pow(x, 2) + pow(y, 2)) + 4.0*u_0*cos(pow(x, 2) + pow(y, 2)) + 4*v_0*x*y*cos(pow(x, 2) + pow(y, 2))) + 2*Viscosity*(2*u_0*pow(y, 2)*sin(pow(x, 2) + pow(y, 2)) - u_0*cos(pow(x, 2) + pow(y, 2)) + 2*v_0*x*y*cos(pow(x, 2) + pow(y, 2)));
  val_source[2] = -2*Density*u_0*v_0*x*(epsilon + sin(pow(x, 2) + pow(y, 2)))*sin(pow(x, 2) + pow(y, 2)) + 2*Density*u_0*v_0*x*(epsilon + cos(pow(x, 2) + pow(y, 2)))*cos(pow(x, 2) + pow(y, 2)) - 4*Density*pow(v_0, 2)*y*(epsilon + cos(pow(x, 2) + pow(y, 2)))*sin(pow(x, 2) + pow(y, 2)) + 2*P_0*y*cos(pow(x, 2) + pow(y, 2)) + 0.666666666666667*Viscosity*(-4*u_0*x*y*sin(pow(x, 2) + pow(y, 2)) + 8.0*v_0*pow(y, 2)*cos(pow(x, 2) + pow(y, 2)) + 4.0*v_0*sin(pow(x, 2) + pow(y, 2))) + 2*Viscosity*(2*u_0*x*y*sin(pow(x, 2) + pow(y, 2)) + 2*v_0*pow(x, 2)*cos(pow(x, 2) + pow(y, 2)) + v_0*sin(pow(x, 2) + pow(y, 2)));
  val_source[3]      = 0.0;
  val_source[nVar-1] = 0.0;
  
}

bool CMMSIncNSSolution::IsManufacturedSolution(void) {
  return true;
}

CUserDefinedSolution::CUserDefinedSolution(void) : CVerificationSolution() { }

CUserDefinedSolution::CUserDefinedSolution(unsigned short val_nDim,
                                           unsigned short val_nVar,
                                           CConfig*       config)
  : CVerificationSolution(val_nDim, val_nVar, config) {
  
  /*--- Write a message that the solution is initialized for a
   user-defined verification case. ---*/
  
  if (rank == MASTER_NODE) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for a user-defined verification case!!!" << endl;
    cout << endl << flush;
  }

  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
}

CUserDefinedSolution::~CUserDefinedSolution(void) { }

void CUserDefinedSolution::GetBCState(const unsigned short val_nParams,
                                      const su2double      *val_params,
                                      const su2double      *val_coords,
                                      const su2double      val_t,
                                      su2double            *val_solution) {

  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
}

void CUserDefinedSolution::GetSolution(const unsigned short val_nParams,
                                       const su2double      *val_params,
                                       const su2double      *val_coords,
                                       const su2double      val_t,
                                       su2double            *val_solution) {

  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
}

void CUserDefinedSolution::GetMMSSourceTerm(const unsigned short val_nParams,
                                            const su2double      *val_params,
                                            const su2double      *val_coords,
                                            const su2double      val_t,
                                            su2double            *val_source) {

  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
}

bool CUserDefinedSolution::IsManufacturedSolution(void) {
  SU2_MPI::Error("User must implement this function", CURRENT_FUNCTION);
  return false;  /* True if manufactured. */
}
