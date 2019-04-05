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
  Viscosity = config->GetViscosity_FreeStreamND();

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
  Viscosity    = config->GetViscosity_FreeStreamND();
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

  /*--- Coefficients, needed to determine the solution and MMS source term.
        Note that the dimensional values must be got. ---*/
  Gamma   = config->GetGamma();
  RGas    = config->GetGas_Constant();
  muRef   = config->GetMu_Ref();
  TRef    = config->GetMu_Temperature_Ref();
  S       = config->GetMu_S();
  Prandtl = config->GetPrandtl_Lam();

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

  /*--- Constants, which describe this manufactured solution. The temperature only
        varies radially in order to have the zero normal gradient as the inner wall
        of r == 1. The other variables vary in a more complicated way.
        First the constants for the pressure solution. ---*/
  a_Px  =      1.40;
  a_Py  =      1.25;
  a_Pz  =      1.90;
  a_Pxy =      0.75;
  a_Pxz =      0.83;
  a_Pyz =      1.17;
  P_0   = 100000.0;
  P_x   = -15000.0;
  P_y   =  10000.0;
  P_z   =  12500.0;
  P_xy  = -17500.0;
  P_xz  =   7500.0;
  P_yz  =  -5000.0;

  /*--- The constants for the temperature solution. ---*/
  a_T1 =  2.1;
  a_T2 = -1.9;
  a_T3 = -2.0;
  a_T4 =  2.5;

  /*--- The constants for the x-velocity solution. ---*/
  a_ux  =   1.6666666667;
  a_uy  =   1.5;
  a_uz  =   2.3333333333;
  a_uxy =   0.6;
  a_uxz =   0.9777777777;
  a_uyz =   1.0833333333;
  u_0   =  70.0;
  u_x   =   4.0;
  u_y   = -12.0;
  u_z   = -17.0;
  u_xy  =   7.0;
  u_xz  = -10.5;
  u_yz  =   8.3333333333;

  /*--- The constants for the y-velocity solution. ---*/
  a_vx  =   1.4;
  a_vy  =   2.2;
  a_vz  =   1.8;
  a_vxy =   0.9;
  a_vxz =   1.111111;
  a_vyz =   0.666667;
  v_0   =  60.0;
  v_x   = -20.0;
  v_y   =   9.0;
  v_z   =  17.333333;
  v_xy  = -11.0;
  v_xz  =  -5.666667;
  v_yz  =   7.133333;

  /*--- The constants for the z-velocity solution. ---*/
  a_wx  =   2.4;
  a_wy  =   0.8;
  a_wz  =   1.222222;
  a_wxy =   0.3;
  a_wxz =   1.7;
  a_wyz =   0.75;
  w_0   =  30.0;
  w_x   =  13.0;
  w_y   = -24.0;
  w_z   =   8.5;
  w_xy  =   9.8;
  w_xz  = -17.0;
  w_yz  =  12.3333333; 

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

  if(config->GetKind_ViscosityModel() != SUTHERLAND)
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

  /* Determine the dimensional solution for the pressure. */
  const su2double Pi = PI_NUMBER;

  su2double p = P_0 + P_x *cos(a_Px *Pi*x)
              +       P_y *sin(a_Py *Pi*y)
              +       P_z *cos(a_Pz *Pi*z)
              +       P_xy*sin(a_Pxy*Pi*x*y)
              +       P_xz*sin(a_Pxz*Pi*x*z)
              +       P_yz*cos(a_Pyz*Pi*y*z);

  /* Determine the dimensional solution for the temperature. */
  const su2double r    = sqrt(x*x + y*y + z*z);
  const su2double fact = (r-1.0)*(r-1.0)/(a_T1 + a_T2 + a_T3 + a_T4);

  su2double T = 0.5*TWall*(1.0 + fact*(a_T1*cos(Pi*(r-2.0))
              +                        a_T2*cos(Pi*(r-2.0)*2.0)
              +                        a_T3*cos(Pi*(r-2.0)*3.0)
              +                        a_T4*cos(Pi*(r-2.0)*4.0)));

  /* Determine the dimensional solution for the x-velocity. */
  su2double u = (u_0 + u_x *sin(a_ux *Pi*x)
              +        u_y *cos(a_uy *Pi*y)
              +        u_z *cos(a_uz *Pi*z)
              +        u_xy*cos(a_uxy*Pi*x*y)
              +        u_xz*sin(a_uxz*Pi*x*z)
              +        u_yz*sin(a_uyz*Pi*y*z))*(r-1.0)*(2.0-r)*4.0;

  /* Determine the dimensional solution for the y-velocity. */
  su2double v = (v_0 + v_x *cos(a_vx *Pi*x)
              +        v_y *sin(a_vy *Pi*y)
              +        v_z *sin(a_vz *Pi*z)
              +        v_xy*cos(a_vxy*Pi*x*y)
              +        v_xz*sin(a_vxz*Pi*x*z)
              +        v_yz*sin(a_vyz*Pi*y*z))*(r-1.0)*(2.0-r)*4.0;

  /* Determine the dimensional solution for the y-velocity. */
  su2double w = (w_0 + w_x *cos(a_wx *Pi*x)
              +        w_y *cos(a_wy *Pi*y)
              +        w_z *cos(a_wz *Pi*z)
              +        w_xy*sin(a_wxy*Pi*x*y)
              +        w_xz*sin(a_wxz*Pi*x*z)
              +        w_yz*sin(a_wyz*Pi*y*z))*(r-1.0)*(2.0-r)*4.0;

  /* Compute the density and determine the non-dimensional solution. */
  su2double rho = p/(RGas*T);

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
  const su2double t3 = a_Px * Pi * x;
  const su2double t4 = sin(t3);
  const su2double t6 = t4 * Pi * P_x * a_Px;
  const su2double t7 = P_xy * a_Pxy;
  const su2double t8 = y * Pi;
  const su2double t10 = x * y;
  const su2double t11 = t10 * a_Pxy * Pi;
  const su2double t12 = cos(t11);
  const su2double t14 = t12 * t8 * t7;
  const su2double t15 = P_xz * a_Pxz;
  const su2double t16 = z * Pi;
  const su2double t18 = x * z;
  const su2double t19 = t18 * a_Pxz * Pi;
  const su2double t20 = cos(t19);
  const su2double t22 = t20 * t16 * t15;
  const su2double t23 = -t6 + t14 + t22;
  const su2double t24 = 0.1e1 / RGas;
  const su2double t25 = t24 * t23;
  const su2double t26 = 0.1e1 / TWall;
  const su2double t27 = t26 * t25;
  const su2double t28 = x * x;
  const su2double t29 = y * y;
  const su2double t30 = z * z;
  const su2double t31 = t28 + t29 + t30;
  const su2double t32 = sqrt(t31);
  const su2double t33 = t32 - 0.1e1;
  const su2double t34 = t33 * t33;
  const su2double t35 = a_T1 + a_T2 + a_T3 + a_T4;
  const su2double t36 = 0.1e1 / t35;
  const su2double t37 = t36 * t34;
  const su2double t38 = t32 - 0.2e1;
  const su2double t39 = t38 * Pi;
  const su2double t40 = cos(t39);
  const su2double t42 = 0.2e1 * t39;
  const su2double t43 = cos(t42);
  const su2double t45 = 0.3e1 * t39;
  const su2double t46 = cos(t45);
  const su2double t48 = 0.4e1 * t39;
  const su2double t49 = cos(t48);
  const su2double t51 = t40 * a_T1 + t43 * a_T2 + t46 * a_T3 + t49 * a_T4;
  const su2double t53 = t51 * t37 + 0.1e1;
  const su2double t54 = 0.2e1 / t53;
  const su2double t55 = a_ux * Pi;
  const su2double t56 = x * t55;
  const su2double t57 = sin(t56);
  const su2double t58 = t57 * u_x;
  const su2double t59 = a_uy * Pi;
  const su2double t60 = y * t59;
  const su2double t61 = cos(t60);
  const su2double t62 = t61 * u_y;
  const su2double t63 = a_uz * Pi;
  const su2double t64 = z * t63;
  const su2double t65 = cos(t64);
  const su2double t66 = t65 * u_z;
  const su2double t67 = a_uxy * Pi;
  const su2double t68 = t10 * t67;
  const su2double t69 = cos(t68);
  const su2double t70 = t69 * u_xy;
  const su2double t71 = a_uxz * Pi;
  const su2double t72 = t18 * t71;
  const su2double t73 = sin(t72);
  const su2double t74 = t73 * u_xz;
  const su2double t76 = y * z;
  const su2double t77 = t76 * a_uyz * Pi;
  const su2double t78 = sin(t77);
  const su2double t79 = t78 * u_yz;
  const su2double t80 = u_0 + t58 + t62 + t66 + t70 + t74 + t79;
  const su2double t82 = -t38 * t33;
  const su2double t85 = cos(t3);
  const su2double t86 = t85 * P_x;
  const su2double t88 = a_Py * Pi * y;
  const su2double t89 = sin(t88);
  const su2double t90 = t89 * P_y;
  const su2double t92 = a_Pz * Pi * z;
  const su2double t93 = cos(t92);
  const su2double t94 = t93 * P_z;
  const su2double t95 = sin(t11);
  const su2double t96 = t95 * P_xy;
  const su2double t97 = sin(t19);
  const su2double t98 = t97 * P_xz;
  const su2double t100 = t76 * a_Pyz * Pi;
  const su2double t101 = cos(t100);
  const su2double t102 = t101 * P_yz;
  const su2double t103 = P_0 + t86 + t90 + t94 + t96 + t98 + t102;
  const su2double t104 = t24 * t103;
  const su2double t106 = 0.4e1 * pow(t53, -0.2e1);
  const su2double t108 = t106 * t26 * t104;
  const su2double t109 = t33 * t80;
  const su2double t110 = t36 * t33;
  const su2double t111 = 0.1e1 / t32;
  const su2double t112 = t111 * t51;
  const su2double t115 = a_T1 * Pi;
  const su2double t116 = x * t111;
  const su2double t117 = sin(t39);
  const su2double t118 = t117 * t116;
  const su2double t120 = a_T2 * Pi;
  const su2double t121 = sin(t42);
  const su2double t125 = a_T3 * Pi;
  const su2double t126 = sin(t45);
  const su2double t130 = a_T4 * Pi;
  const su2double t131 = sin(t48);
  const su2double t138 = x * t112 * t110 + (-0.2e1 * t121 * t116 * t120 - 0.3e1 * t126 * t116 * t125 - 0.4e1 * t131 * t116 * t130 - t118 * t115) * t37 / 0.2e1;
  const su2double t142 = t26 * t104;
  const su2double t144 = cos(t56);
  const su2double t147 = u_xy * a_uxy;
  const su2double t148 = sin(t68);
  const su2double t151 = u_xz * a_uxz;
  const su2double t152 = cos(t72);
  const su2double t155 = t144 * Pi * u_x * a_ux - t148 * t8 * t147 + t152 * t16 * t151;
  const su2double t159 = t54 * t26;
  const su2double t160 = t159 * t104;
  const su2double t161 = t111 * t80;
  const su2double t162 = -t38 * x;
  const su2double t163 = t162 * t161;
  const su2double t165 = t116 * t109;
  const su2double t168 = cos(t88);
  const su2double t170 = t168 * Pi * P_y * a_Py;
  const su2double t171 = x * Pi;
  const su2double t173 = t12 * t171 * t7;
  const su2double t174 = P_yz * a_Pyz;
  const su2double t175 = sin(t100);
  const su2double t177 = t175 * t16 * t174;
  const su2double t178 = t170 + t173 - t177;
  const su2double t179 = t24 * t178;
  const su2double t180 = t26 * t179;
  const su2double t181 = a_vx * Pi;
  const su2double t182 = x * t181;
  const su2double t183 = cos(t182);
  const su2double t184 = t183 * v_x;
  const su2double t185 = a_vy * Pi;
  const su2double t186 = y * t185;
  const su2double t187 = sin(t186);
  const su2double t188 = t187 * v_y;
  const su2double t189 = a_vz * Pi;
  const su2double t190 = z * t189;
  const su2double t191 = sin(t190);
  const su2double t192 = t191 * v_z;
  const su2double t193 = a_vxy * Pi;
  const su2double t194 = t10 * t193;
  const su2double t195 = cos(t194);
  const su2double t196 = t195 * v_xy;
  const su2double t197 = a_vxz * Pi;
  const su2double t198 = t18 * t197;
  const su2double t199 = sin(t198);
  const su2double t200 = t199 * v_xz;
  const su2double t201 = a_vyz * Pi;
  const su2double t202 = t76 * t201;
  const su2double t203 = sin(t202);
  const su2double t204 = t203 * v_yz;
  const su2double t205 = v_0 + t184 + t188 + t192 + t196 + t200 + t204;
  const su2double t209 = t33 * t205;
  const su2double t212 = y * t111;
  const su2double t213 = t117 * t212;
  const su2double t227 = y * t112 * t110 + (-0.2e1 * t121 * t212 * t120 - 0.3e1 * t126 * t212 * t125
                       -                     0.4e1 * t131 * t212 * t130 - t213 * t115) * t37 / 0.2e1;
  const su2double t232 = cos(t186);
  const su2double t235 = v_xy * a_vxy;
  const su2double t236 = sin(t194);
  const su2double t237 = t236 * t171;
  const su2double t239 = v_yz * a_vyz;
  const su2double t240 = cos(t202);
  const su2double t243 = t232 * Pi * v_y * a_vy + t240 * t16 * t239 - t237 * t235;
  const su2double t247 = t111 * t205;
  const su2double t248 = -t38 * y;
  const su2double t249 = t248 * t247;
  const su2double t251 = t212 * t209;
  const su2double t254 = sin(t92);
  const su2double t256 = t254 * Pi * P_z * a_Pz;
  const su2double t258 = t20 * t171 * t15;
  const su2double t260 = t175 * t8 * t174;
  const su2double t261 = -t256 + t258 - t260;
  const su2double t262 = t24 * t261;
  const su2double t263 = t26 * t262;
  const su2double t264 = a_wx * Pi;
  const su2double t265 = x * t264;
  const su2double t266 = cos(t265);
  const su2double t267 = t266 * w_x;
  const su2double t268 = a_wy * Pi;
  const su2double t269 = y * t268;
  const su2double t270 = cos(t269);
  const su2double t271 = t270 * w_y;
  const su2double t272 = a_wz * Pi;
  const su2double t273 = z * t272;
  const su2double t274 = cos(t273);
  const su2double t275 = t274 * w_z;
  const su2double t276 = a_wxy * Pi;
  const su2double t277 = t10 * t276;
  const su2double t278 = sin(t277);
  const su2double t279 = t278 * w_xy;
  const su2double t280 = a_wxz * Pi;
  const su2double t281 = t18 * t280;
  const su2double t282 = sin(t281);
  const su2double t283 = t282 * w_xz;
  const su2double t284 = a_wyz * Pi;
  const su2double t285 = t76 * t284;
  const su2double t286 = sin(t285);
  const su2double t287 = t286 * w_yz;
  const su2double t288 = w_0 + t267 + t271 + t275 + t279 + t283 + t287;
  const su2double t292 = t33 * t288;
  const su2double t295 = z * t111;
  const su2double t296 = t117 * t295;
  const su2double t310 = z * t112 * t110 + (-0.2e1 * t121 * t295 * t120 - 0.3e1 * t126 * t295 * t125
                       -                     0.4e1 * t131 * t295 * t130 - t296 * t115) * t37 / 0.2e1;
  const su2double t315 = sin(t273);
  const su2double t318 = w_xz * a_wxz;
  const su2double t319 = cos(t281);
  const su2double t320 = t319 * t171;
  const su2double t322 = w_yz * a_wyz;
  const su2double t323 = cos(t285);
  const su2double t324 = t323 * t8;
  const su2double t326 = -t315 * Pi * w_z * a_wz + t320 * t318 + t324 * t322;
  const su2double t330 = t111 * t288;
  const su2double t331 = -t38 * z;
  const su2double t332 = t331 * t330;
  const su2double t334 = t295 * t292;
  const su2double t336 = t138 * t38 * t109 * t108 + t227 * t38 * t209 * t108
                       + t310 * t38 * t292 * t108 + t82 * t155 * t54 * t142 + t82 * t243 * t54 * t142
                       + t82 * t326 * t54 * t142 + t82 * t205 * t54 * t180 + t82 * t288 * t54 * t263
                       + t82 * t80 * t54 * t27 + t163 * t160 - t165 * t160 + t249 * t160 - t251 * t160
                       + t332 * t160 - t334 * t160;
  const su2double t337 = sqrt(0.8e1);
  const su2double t339 = t28 + t29 + t30 - 0.2e1 * t32 + 0.1e1;
  const su2double t340 = t339 * a_T4;
  const su2double t341 = t40 * t40;
  const su2double t342 = t341 * t341;
  const su2double t343 = t342 * t340;
  const su2double t344 = t339 * a_T3;
  const su2double t345 = t341 * t40;
  const su2double t346 = t345 * t344;
  const su2double t349 = a_T2 - 0.4e1 * a_T4;
  const su2double t350 = t339 * t349;
  const su2double t351 = t341 * t350;
  const su2double t354 = a_T1 - 0.3e1 * a_T3;
  const su2double t355 = t339 * t354;
  const su2double t356 = t40 * t355;
  const su2double t358 = -a_T4 + a_T2;
  const su2double t366 = t343 + t346 / 0.2e1 + t351 / 0.4e1 + t356 / 0.8e1 + t32 * t358 / 0.4e1
                       - t28 * t358 / 0.8e1 - t29 * t358 / 0.8e1 - t30 * t358 / 0.8e1
                       + a_T1 / 0.8e1 + a_T3 / 0.8e1 + a_T4 / 0.4e1;
  const su2double t368 = 0.1e1 / TRef;
  const su2double t371 = sqrt(t36 * t368 * t366 * TWall);
  const su2double t372 = 0.1e1 / t371;
  const su2double t373 = t372 * t337;
  const su2double t374 = TRef + S;
  const su2double t375 = t374 * t373;
  const su2double t376 = TWall * TWall;
  const su2double t377 = sqrt(0.2e1);
  const su2double t378 = t377 * t376;
  const su2double t379 = v_0 * y;
  const su2double t380 = 0.3e1 / 0.2e1 * t379;
  const su2double t381 = w_0 * z;
  const su2double t382 = u_0 * x;
  const su2double t384 = y * t184;
  const su2double t385 = y * t188;
  const su2double t386 = y * t192;
  const su2double t387 = y * t196;
  const su2double t388 = y * t200;
  const su2double t389 = y * t204;
  const su2double t390 = x * t58;
  const su2double t392 = x * t62;
  const su2double t394 = x * t66;
  const su2double t396 = x * t70;
  const su2double t398 = x * t74;
  const su2double t400 = x * t79;
  const su2double t402 = t379 + t381 - 0.2e1 * t382 + t384 + t385 + t386 + t387 + t388 + t389
                       - 0.2e1 * t390 - 0.2e1 * t392 - 0.2e1 * t394 - 0.2e1 * t396
                       - 0.2e1 * t398 - 0.2e1 * t400;
  const su2double t403 = z * t267;
  const su2double t404 = z * t271;
  const su2double t405 = z * t275;
  const su2double t406 = z * t279;
  const su2double t407 = z * t283;
  const su2double t408 = z * t287;
  const su2double t409 = z * t151;
  const su2double t410 = t28 + t29 + t30 + 0.2e1;
  const su2double t411 = t410 * Pi;
  const su2double t412 = t152 * t411;
  const su2double t413 = t412 * t409;
  const su2double t414 = y * t147;
  const su2double t415 = t148 * t411;
  const su2double t416 = t415 * t414;
  const su2double t417 = w_xz * x;
  const su2double t418 = a_wxz * t417;
  const su2double t419 = t319 * t411;
  const su2double t420 = t419 * t418;
  const su2double t422 = y * t322;
  const su2double t423 = t323 * t411;
  const su2double t424 = t423 * t422;
  const su2double t426 = z * t239;
  const su2double t427 = t240 * t411;
  const su2double t429 = t427 * t426 / 0.2e1;
  const su2double t430 = v_xy * x;
  const su2double t431 = a_vxy * t430;
  const su2double t432 = t236 * t411;
  const su2double t434 = t432 * t431 / 0.2e1;
  const su2double t435 = t410 * w_z;
  const su2double t437 = t315 * t435 * t272;
  const su2double t442 = t144 * t410 * a_ux * Pi * u_x;
  const su2double t443 = t410 * v_y;
  const su2double t446 = t232 * t443 * t185 / 0.2e1;
  const su2double t447 = t403 + t404 + t405 + t406 + t407 + t408 - t413 + t416 + t420 / 0.2e1
                       + t424 / 0.2e1 + t429 - t434 - t437 / 0.2e1 - t442 + t446;
  const su2double t448 = t447 + t402;
  const su2double t452 = 0.3e1 / 0.2e1 * t384;
  const su2double t453 = 0.3e1 / 0.2e1 * t385;
  const su2double t454 = 0.3e1 / 0.2e1 * t386;
  const su2double t455 = 0.3e1 / 0.2e1 * t387;
  const su2double t456 = 0.3e1 / 0.2e1 * t388;
  const su2double t457 = 0.3e1 / 0.2e1 * t389;
  const su2double t463 = -t380 + t32 * t448 - 0.3e1 / 0.2e1 * t381 + 0.3e1 * t382 - t452 - t453
                       - t454 - t455 - t456 - t457 + 0.3e1 * t390 + 0.3e1 * t392 + 0.3e1 * t394
                       + 0.3e1 * t396 + 0.3e1 * t398;
  const su2double t471 = t319 * t31;
  const su2double t472 = Pi * t471;
  const su2double t473 = t418 * t472;
  const su2double t475 = t323 * t31;
  const su2double t476 = Pi * t475;
  const su2double t477 = t422 * t476;
  const su2double t479 = t31 * t148;
  const su2double t480 = Pi * t479;
  const su2double t481 = t414 * t480;
  const su2double t484 = Pi * t31 * t152;
  const su2double t485 = t409 * t484;
  const su2double t487 = t240 * t31;
  const su2double t488 = Pi * t487;
  const su2double t490 = 0.3e1 / 0.2e1 * t426 * t488;
  const su2double t491 = t236 * t31;
  const su2double t492 = Pi * t491;
  const su2double t494 = 0.3e1 / 0.2e1 * t431 * t492;
  const su2double t496 = u_x * t55;
  const su2double t497 = t496 * t31 * t144;
  const su2double t501 = w_z * t272 * t315 * t31;
  const su2double t506 = 0.3e1 / 0.2e1 * v_y * t185 * t232 * t31;
  const su2double t507 = 0.3e1 * t400 - 0.3e1 / 0.2e1 * t403 - 0.3e1 / 0.2e1 * t404
                       - 0.3e1 / 0.2e1 * t405 - 0.3e1 / 0.2e1 * t406 - 0.3e1 / 0.2e1 * t407
                       - 0.3e1 / 0.2e1 * t408 - 0.3e1 / 0.2e1 * t473 - 0.3e1 / 0.2e1 * t477
                       - 0.3e1 * t481 + 0.3e1 * t485 - t490 + t494 + 0.3e1 * t497
                       + 0.3e1 / 0.2e1 * t501 - t506;
  const su2double t508 = t507 + t463;
  const su2double t519 = -t28 * t358 - t29 * t358 - t30 * t358 + 0.2e1 * t32 * t358 + a_T1 + a_T3
                       + 0.2e1 * a_T4 + 0.8e1 * t343 + 0.4e1 * t346 + 0.2e1 * t351 + t356;
  const su2double t521 = t28 + t29 + t30 + 0.1e1;
  const su2double t526 = t32 * t521 - 0.2e1 * t28 - 0.2e1 * t29 - 0.2e1 * t30;
  const su2double t527 = t526 * TWall;
  const su2double t528 = t342 * a_T4;
  const su2double t530 = 0.4e1 * t528 * t527;
  const su2double t531 = TWall * a_T3;
  const su2double t534 = 0.2e1 * t345 * t526 * t531;
  const su2double t535 = t341 * t349;
  const su2double t536 = t535 * t527;
  const su2double t537 = t40 * t354;
  const su2double t539 = t537 * t527 / 0.2e1;
  const su2double t540 = t358 * TWall;
  const su2double t547 = a_T1 / 0.2e1;
  const su2double t548 = a_T3 / 0.2e1;
  const su2double t551 = t35 * S;
  const su2double t552 = -t28 * t540 / 0.2e1 - t29 * t540 / 0.2e1 - t30 * t540 / 0.2e1
                       + TWall * (t547 + t548 + a_T4) + t551;
  const su2double t555 = t358 * t31 * TWall;
  const su2double t556 = t32 * t552 + t530 + t534 + t536 + t539 + t555;
  const su2double t557 = 0.1e1 / t556;
  const su2double t558 = t557 * t519 * muRef;
  const su2double t559 = TRef * TRef;
  const su2double t560 = 0.1e1 / t559;
  const su2double t561 = x - t116;
  const su2double t563 = 0.2e1 * t342 * t561 * a_T4;
  const su2double t564 = t345 * t340;
  const su2double t565 = t111 * Pi;
  const su2double t567 = t117 * x * t565;
  const su2double t568 = t567 * t564;
  const su2double t571 = 0.2e1 * t345 * t561 * a_T3;
  const su2double t573 = t341 * t344;
  const su2double t574 = t567 * t573;
  const su2double t577 = 0.2e1 * t341 * t561 * t349;
  const su2double t579 = t40 * t350;
  const su2double t580 = t567 * t579;
  const su2double t583 = 0.2e1 * t40 * t561 * t354;
  const su2double t585 = Pi * t355;
  const su2double t586 = t118 * t585;
  const su2double t588 = t111 * t358 / 0.4e1;
  const su2double t592 = t563 - 0.4e1 * t568 + t571 / 0.2e1 - 0.3e1 / 0.2e1 * t574 + t577 / 0.4e1
                       - t580 / 0.2e1 + t583 / 0.8e1 - t586 / 0.8e1 + x * t588 - x * t358 / 0.4e1;
  const su2double t598 = a_uyz * u_yz;
  const su2double t599 = z * t598;
  const su2double t600 = cos(t77);
  const su2double t601 = t600 * t411;
  const su2double t604 = a_vxz * v_xz;
  const su2double t605 = z * t604;
  const su2double t606 = cos(t198);
  const su2double t607 = t606 * t411;
  const su2double t610 = u_xy * x;
  const su2double t611 = a_uxy * t610;
  const su2double t614 = y * t235;
  const su2double t617 = y * t70;
  const su2double t618 = x * t196;
  const su2double t619 = y * t74;
  const su2double t620 = y * t79;
  const su2double t621 = x * t200;
  const su2double t622 = x * t204;
  const su2double t623 = t410 * u_y;
  const su2double t624 = sin(t60);
  const su2double t628 = Pi * v_x;
  const su2double t630 = sin(t182);
  const su2double t634 = x * t184;
  const su2double t635 = x * t188;
  const su2double t636 = x * t192;
  const su2double t637 = y * t62;
  const su2double t638 = y * t66;
  const su2double t639 = y * t58;
  const su2double t640 = v_0 * x;
  const su2double t641 = u_0 * y;
  const su2double t642 = t601 * t599 / 0.2e1 + t607 * t605 / 0.2e1 - t415 * t611 / 0.2e1
                       - t432 * t614 / 0.2e1 + t617 + t618 + t619 + t620 + t621 + t622
                       - t624 * t623 * t59 / 0.2e1 - t630 * t410 * a_vx * t628 / 0.2e1
                       + t634 + t635 + t636 + t637 + t638 + t639 + t640 + t641;
  const su2double t645 = Pi * t600 * t31;
  const su2double t649 = Pi * t606 * t31;
  const su2double t668 = v_x * t181 * t630 * t31 + u_y * t59 * t624 * t31 - t622 - t634 - t635
                       - t636 - t637 - t638 - t639 - t640 - t641;
  const su2double t670 = (0.3e1 / 0.2e1 * t668 + t32 * t642 - 0.3e1 / 0.2e1 * t599 * t645
                       - 0.3e1 / 0.2e1 * t605 * t649 + 0.3e1 / 0.2e1 * t611 * t480
                       + 0.3e1 / 0.2e1 * t614 * t492 - 0.3e1 / 0.2e1 * t617 - 0.3e1 / 0.2e1 * t618
                       - 0.3e1 / 0.2e1 * t619 - 0.3e1 / 0.2e1 * t620 - 0.3e1 / 0.2e1 * t621) * t374;
  const su2double t671 = t519 * t670;
  const su2double t672 = t337 * t377;
  const su2double t673 = t372 * t672;
  const su2double t674 = t673 * t671;
  const su2double t676 = t560 * t376 * muRef;
  const su2double t682 = TWall * (-t28 * t358 / 0.2e1 - t29 * t358 / 0.2e1 - t30 * t358 / 0.2e1
                       + t547 + t548 + a_T4) + t551;
  const su2double t684 = t32 * t682 + t530 + t534 + t536 + t539 + t555;
  const su2double t685 = 0.1e1 / t684;
  const su2double t686 = y - t212;
  const su2double t688 = 0.2e1 * t342 * t686 * a_T4;
  const su2double t690 = t117 * y * t565;
  const su2double t691 = t690 * t564;
  const su2double t694 = 0.2e1 * t345 * t686 * a_T3;
  const su2double t696 = t690 * t573;
  const su2double t699 = 0.2e1 * t341 * t686 * t349;
  const su2double t701 = t690 * t579;
  const su2double t704 = 0.2e1 * t40 * t686 * t354;
  const su2double t706 = t213 * t585;
  const su2double t711 = t688 - 0.4e1 * t691 + t694 / 0.2e1 - 0.3e1 / 0.2e1 * t696 + t699 / 0.4e1
                       - t701 / 0.2e1 + t704 / 0.8e1 - t706 / 0.8e1 + y * t588 - y * t358 / 0.4e1;
  const su2double t713 = t36 * t711 * t685;
  const su2double t716 = u_xz * x;
  const su2double t717 = a_uxz * t716;
  const su2double t720 = y * t598;
  const su2double t723 = a_wxy * w_xy;
  const su2double t724 = y * t723;
  const su2double t725 = cos(t277);
  const su2double t726 = t725 * t411;
  const su2double t729 = z * t318;
  const su2double t732 = z * t70;
  const su2double t733 = z * t74;
  const su2double t734 = z * t79;
  const su2double t735 = x * t279;
  const su2double t736 = x * t283;
  const su2double t737 = x * t287;
  const su2double t738 = t410 * u_z;
  const su2double t739 = sin(t64);
  const su2double t743 = Pi * w_x;
  const su2double t745 = sin(t265);
  const su2double t749 = x * t267;
  const su2double t750 = x * t271;
  const su2double t751 = x * t275;
  const su2double t752 = z * t62;
  const su2double t753 = z * t66;
  const su2double t754 = z * t58;
  const su2double t755 = w_0 * x;
  const su2double t756 = u_0 * z;
  const su2double t757 = t412 * t717 / 0.2e1 + t601 * t720 / 0.2e1 + t726 * t724 / 0.2e1
                       + t419 * t729 / 0.2e1 + t732 + t733 + t734 + t735 + t736 + t737
                       - t739 * t738 * t63 / 0.2e1 - t745 * t410 * a_wx * t743 / 0.2e1
                       + t749 + t750 + t751 + t752 + t753 + t754 + t755 + t756;
  const su2double t764 = Pi * t725 * t31;
  const su2double t775 = t31 * u_z;
  const su2double t781 = w_x * t264 * t745 * t31 + t739 * t775 * t63 - t737 - t749 - t750
                       - t751 - t752 - t753 - t754 - t755 - t756;
  const su2double t783 = (0.3e1 / 0.2e1 * t781 + t32 * t757 - 0.3e1 / 0.2e1 * t484 * t717
                       - 0.3e1 / 0.2e1 * t645 * t720 - 0.3e1 / 0.2e1 * t724 * t764
                       - 0.3e1 / 0.2e1 * t729 * t472 - 0.3e1 / 0.2e1 * t732 - 0.3e1 / 0.2e1 * t733
                       - 0.3e1 / 0.2e1 * t734 - 0.3e1 / 0.2e1 * t735 - 0.3e1 / 0.2e1 * t736) * t374;
  const su2double t784 = t519 * t783;
  const su2double t785 = t673 * t784;
  const su2double t786 = z - t295;
  const su2double t788 = 0.2e1 * t342 * t786 * a_T4;
  const su2double t790 = t117 * z * t565;
  const su2double t791 = t790 * t564;
  const su2double t794 = 0.2e1 * t345 * t786 * a_T3;
  const su2double t796 = t790 * t573;
  const su2double t799 = 0.2e1 * t341 * t786 * t349;
  const su2double t801 = t790 * t579;
  const su2double t804 = 0.2e1 * t40 * t786 * t354;
  const su2double t806 = t296 * t585;
  const su2double t811 = t788 - 0.4e1 * t791 + t794 / 0.2e1 - 0.3e1 / 0.2e1 * t796 + t799 / 0.4e1
                       - t801 / 0.2e1 + t804 / 0.8e1 - t806 / 0.8e1 + z * t588 - z * t358 / 0.4e1;
  const su2double t813 = t36 * t811 * t685;
  const su2double t816 = t371 * t337;
  const su2double t819 = t377 * TWall * t374 * t816;
  const su2double t820 = muRef * t508;
  const su2double t822 = t556 * t556;
  const su2double t823 = 0.1e1 / t822;
  const su2double t824 = t368 * t823;
  const su2double t825 = t32 * x;
  const su2double t827 = t111 * t521;
  const su2double t830 = x * t827 + 0.2e1 * t825 - 0.4e1 * x;
  const su2double t831 = t830 * TWall;
  const su2double t833 = 0.4e1 * t528 * t831;
  const su2double t834 = t345 * a_T4;
  const su2double t835 = t834 * t527;
  const su2double t837 = 0.16e2 * t567 * t835;
  const su2double t840 = 0.2e1 * t345 * t830 * t531;
  const su2double t842 = t341 * t526 * t531;
  const su2double t844 = 0.6e1 * t567 * t842;
  const su2double t845 = t535 * t831;
  const su2double t847 = t40 * t349 * t527;
  const su2double t849 = 0.2e1 * t567 * t847;
  const su2double t851 = t537 * t831 / 0.2e1;
  const su2double t852 = t354 * t527;
  const su2double t854 = t567 * t852 / 0.2e1;
  const su2double t856 = t111 * t552;
  const su2double t859 = 0.2e1 * x * t540;
  const su2double t860 = -t825 * t540 + x * t856 + t833 - t837 + t840 - t844 + t845 - t849 + t851 - t854 + t859;
  const su2double t866 = t337 * t377 * t519;
  const su2double t867 = t866 * t670;
  const su2double t868 = muRef * t371;
  const su2double t869 = TWall * t868;
  const su2double t870 = t684 * t684;
  const su2double t871 = 0.1e1 / t870;
  const su2double t872 = t871 * t368;
  const su2double t877 = 0.2e1 * t32 * y + y * t827 - 0.4e1 * y;
  const su2double t878 = t877 * TWall;
  const su2double t880 = 0.4e1 * t528 * t878;
  const su2double t882 = 0.16e2 * t690 * t835;
  const su2double t885 = 0.2e1 * t345 * t877 * t531;
  const su2double t887 = 0.6e1 * t690 * t842;
  const su2double t888 = t535 * t878;
  const su2double t890 = 0.2e1 * t690 * t847;
  const su2double t892 = t537 * t878 / 0.2e1;
  const su2double t894 = t690 * t852 / 0.2e1;
  const su2double t896 = t32 * TWall;
  const su2double t899 = t111 * t682;
  const su2double t901 = y * TWall;
  const su2double t903 = 0.2e1 * t358 * t901;
  const su2double t904 = -t896 * y * t358 + y * t899 + t880 - t882 + t885 - t887 + t888 - t890 + t892 - t894 + t903;
  const su2double t905 = t904 * t872;
  const su2double t909 = t866 * t783;
  const su2double t914 = 0.2e1 * t32 * z + z * t827 - 0.4e1 * z;
  const su2double t915 = t914 * TWall;
  const su2double t917 = 0.4e1 * t528 * t915;
  const su2double t919 = 0.16e2 * t790 * t835;
  const su2double t922 = 0.2e1 * t345 * t914 * t531;
  const su2double t924 = 0.6e1 * t790 * t842;
  const su2double t925 = t535 * t915;
  const su2double t927 = 0.2e1 * t790 * t847;
  const su2double t929 = t537 * t915 / 0.2e1;
  const su2double t931 = t790 * t852 / 0.2e1;
  const su2double t936 = z * TWall;
  const su2double t938 = 0.2e1 * t358 * t936;
  const su2double t939 = -t896 * z * t358 + z * t899 + t917 - t919 + t922 - t924 + t925 - t927 + t929 - t931 + t938;
  const su2double t940 = t939 * t872;
  const su2double t944 = t34 * t80;
  const su2double t945 = t38 * t38;
  const su2double t946 = t243 * t945;
  const su2double t950 = t159 * t262;
  const su2double t951 = t288 * t945;
  const su2double t952 = t951 * t944;
  const su2double t958 = t152 * t171;
  const su2double t960 = t600 * t8;
  const su2double t962 = -t739 * Pi * u_z * a_uz + t958 * t151 + t960 * t598;
  const su2double t967 = t326 * t945;
  const su2double t971 = t80 * t80;
  const su2double t972 = t34 * t971;
  const su2double t978 = t155 * t945 * t944;
  const su2double t981 = t159 * t179;
  const su2double t982 = t205 * t945;
  const su2double t983 = t982 * t944;
  const su2double t989 = t148 * t171;
  const su2double t993 = -t624 * Pi * u_y * a_uy + t600 * t16 * t598 - t989 * t147;
  const su2double t999 = t945 * t34;
  const su2double t1003 = -0.2e1 / 0.3e1 * t36 * t592 * t560 * t558 * t508 * t378 * t375
                        + t713 * t676 * t674 + t813 * t676 * t785
                        + 0.4e1 / 0.3e1 * t860 * t824 * t519 * t820 * t819
                        - 0.2e1 * t905 * t869 * t867 - 0.2e1 * t940 * t869 * t909
                        + 0.16e2 * t946 * t944 * t160 + 0.16e2 * t952 * t950 
                        + 0.16e2 * t951 * t34 * t962 * t160 + 0.16e2 * t967 * t944 * t160
                        - 0.16e2 * t138 * t945 * t972 * t108 + 0.32e2 * t978 * t160
                        + 0.16e2 * t983 * t981 + 0.16e2 * t982 * t34 * t993 * t160
                        + 0.16e2 * t999 * t971 * t54 * t27 + t14;
  const su2double t1004 = t33 * t971;
  const su2double t1005 = t111 * t945;
  const su2double t1006 = x * t1005;
  const su2double t1007 = t1006 * t1004;
  const su2double t1010 = -t111 * t38;
  const su2double t1011 = x * t1010;
  const su2double t1012 = t1011 * t972;
  const su2double t1015 = t310 * t951;
  const su2double t1023 = a_uz * a_uz;
  const su2double t1024 = Pi * Pi;
  const su2double t1025 = t1024 * t1023;
  const su2double t1029 = a_wxz * a_wxz;
  const su2double t1031 = t1024 * t1029 * t31;
  const su2double t1032 = t282 * z;
  const su2double t1035 = 0.3e1 / 0.2e1 * t417 * t1032 * t1031;
  const su2double t1038 = t739 * z * u_z * t63;
  const su2double t1042 = t315 * x * w_z * t272;
  const su2double t1044 = w_xz * t28;
  const su2double t1046 = t319 * t280 * t1044;
  const su2double t1050 = z * w_x * t745 * t264;
  const su2double t1052 = u_xz * t28;
  const su2double t1053 = a_uxz * a_uxz;
  const su2double t1054 = t1053 * t1052;
  const su2double t1055 = t31 * t1024;
  const su2double t1059 = a_uyz * a_uyz;
  const su2double t1060 = t1059 * u_yz;
  const su2double t1061 = t29 * t1060;
  const su2double t1065 = t958 * t409;
  const su2double t1068 = t323 * t171 * t422;
  const su2double t1072 = z * w_xy * t725;
  const su2double t1073 = t1072 * y * t276;
  const su2double t1075 = 0.3e1 / 0.2e1 * t65 * t775 * t1025 + t1035 + 0.9e1 / 0.2e1 * t1038
                        + 0.3e1 / 0.2e1 * t1042 - 0.3e1 / 0.2e1 * t1046 + 0.3e1 * t1050
                        + 0.3e1 / 0.2e1 * t73 * t1055 * t1054 + 0.3e1 / 0.2e1 * t78 * t1055 * t1061
                        - 0.9e1 / 0.2e1 * t1065 - 0.3e1 / 0.2e1 * t1068 - 0.3e1 * t1073;
  const su2double t1078 = w_xz * t319 * t30 * t280;
  const su2double t1082 = 0.3e1 / 0.2e1 * w_xz * t280 * t471;
  const su2double t1083 = t960 * t599;
  const su2double t1084 = 0.9e1 / 0.2e1 * t1083;
  const su2double t1085 = 0.3e1 / 0.2e1 * u_0;
  const su2double t1086 = 0.3e1 / 0.2e1 * t58;
  const su2double t1087 = 0.3e1 / 0.2e1 * t62;
  const su2double t1088 = 0.3e1 / 0.2e1 * t66;
  const su2double t1089 = 0.3e1 / 0.2e1 * t70;
  const su2double t1090 = 0.3e1 / 0.2e1 * t74;
  const su2double t1091 = 0.3e1 / 0.2e1 * t79;
  const su2double t1093 = t410 * t1024;
  const su2double t1094 = t73 * t1093;
  const su2double t1097 = 0.2e1 * t1083;
  const su2double t1098 = t78 * t1093;
  const su2double t1102 = t419 * t318 / 0.2e1;
  const su2double t1106 = t1032 * t1093 * t1029 * t417 / 0.2e1;
  const su2double t1112 = t79 + t1046 + t1068 - 0.2e1 * t1038 - t65 * t738 * t1025 / 0.2e1
                        - t1050 - t1042 + t62 + t66 + t58 + u_0;
  const su2double t1115 = t111 * t757;
  const su2double t1117 = -0.3e1 * t1078 - t1082 - t1084 - t1085 - t1086 - t1087 - t1088
                        - t1089 - t1090 - t1091 + t32 * (t1112 + 0.2e1 * t1065
                        - t1094 * t1054 / 0.2e1 + t1097 - t1098 * t1061 / 0.2e1 + t1073
                        + t1102 + t1078 - t1106 + t70 + t74) + z * t1115;
  const su2double t1119 = (t1117 + t1075) * t374;
  const su2double t1121 = t368 * TWall;
  const su2double t1122 = t685 * t1121;
  const su2double t1123 = t1122 * t868;
  const su2double t1132 = 0.2e1 * t111 * t358;
  const su2double t1136 = z * t1132 - 0.2e1 * z * t358 + 0.8e1 * t788 - 0.32e2 * t791
                        + 0.4e1 * t794 - 0.12e2 * t796 + 0.2e1 * t799 - 0.4e1 * t801 + t804 - t806;
  const su2double t1151 = x * t1132 - 0.2e1 * x * t358 + 0.8e1 * t563 - 0.32e2 * t568
                        + 0.4e1 * t571 - 0.12e2 * t574 + 0.2e1 * t577 - 0.4e1 * t580 + t583 - t586;
  const su2double t1158 = t80 * t159 * t104;
  const su2double t1159 = t945 * t33;
  const su2double t1160 = y * t247;
  const su2double t1164 = -t38 * t34;
  const su2double t1168 = a_uy * a_uy;
  const su2double t1176 = t624 * y * u_y * t59;
  const su2double t1180 = v_xy * t236 * t29 * t193;
  const su2double t1184 = t232 * x * v_y * t185;
  const su2double t1186 = v_xy * t28;
  const su2double t1188 = t236 * t193 * t1186;
  const su2double t1192 = y * v_x * t630 * t181;
  const su2double t1197 = y * v_xz * t606 * z * t197;
  const su2double t1199 = t989 * t414;
  const su2double t1202 = t240 * t171 * t426;
  const su2double t1206 = 0.3e1 / 0.2e1 * v_xy * t193 * t491;
  const su2double t1207 = 0.3e1 / 0.2e1 * u_y * t61 * t1024 * t1168 * t31 + 0.9e1 / 0.2e1 * t1176
                        + 0.3e1 * t1180 - 0.3e1 / 0.2e1 * t1184 + 0.3e1 / 0.2e1 * t1188 + 0.3e1 * t1192
                        - 0.3e1 * t1197 + 0.9e1 / 0.2e1 * t1199 - 0.3e1 / 0.2e1 * t1202 + t1206 - t1084;
  const su2double t1218 = u_xy * t28;
  const su2double t1219 = a_uxy * a_uxy;
  const su2double t1221 = t69 * t1093;
  const su2double t1225 = t432 * t235 / 0.2e1;
  const su2double t1226 = a_vxy * a_vxy;
  const su2double t1228 = t195 * y;
  const su2double t1231 = t1228 * t1093 * t1226 * t430 / 0.2e1;
  const su2double t1238 = t74 + t79 + t1202 - 0.2e1 * t1176 - t61 * t623 * t1024 * t1168 / 0.2e1
                        - t1192 + t1184 + t62 + t66 + t58 + u_0;
  const su2double t1241 = t111 * t642;
  const su2double t1244 = t1024 * t1226 * t31;
  const su2double t1247 = 0.3e1 / 0.2e1 * t430 * t1228 * t1244;
  const su2double t1248 = t1024 * t1219;
  const su2double t1250 = t31 * t69;
  const su2double t1251 = u_xy * t1250;
  const su2double t1254 = 0.3e1 / 0.2e1 * u_yz * t78 * t30 * t1024 * t1059 * t31 - t1085 - t1086 - t1087
                        - t1088 - t1089 - t1090 - t1091 + t32 * (t1238 + t1097 - t1098 * t30 * t1060 / 0.2e1
                        + t1197 - 0.2e1 * t1199 - t1221 * t1219 * t1218 / 0.2e1 - t1225 - t1180 - t1231
                        + t70 - t1188) + y * t1241 + t1247 + 0.3e1 / 0.2e1 * t1251 * t28 * t1248;
  const su2double t1256 = (t1254 + t1207) * t374;
  const su2double t1269 = y * t1132 - 0.2e1 * y * t358 + 0.8e1 * t688 - 0.32e2 * t691 + 0.4e1 * t694
                        - 0.12e2 * t696 + 0.2e1 * t699 - 0.4e1 * t701 + t704 - t706;
  const su2double t1275 = z * t330;
  const su2double t1276 = t1275 * t1159;
  const su2double t1279 = t1275 * t1164;
  const su2double t1312 = t1094 * t30 * t1053 * u_xz + t1221 * t29 * t1219 * u_xy - 0.4e1 * t1065
                        + t1068 + t1073 - t1106 - t1231 - 0.2e1 * t58 - 0.2e1 * t62 - 0.2e1 * t66
                        - 0.2e1 * t70 - 0.2e1 * t74 - 0.2e1 * t79 - 0.2e1 * u_0;
  const su2double t1315 = t496 * x * t144;
  const su2double t1318 = a_ux * a_ux;
  const su2double t1322 = t57 * t410 * t1318 * u_x * t1024 - t1042 + t1046 - t1050 + t1078 + t1102
                        - t1180 + t1184 - t1188 - t1192 + t1197 + 0.4e1 * t1199 + t1202 - t1225
                        - 0.4e1 * t1315;
  const su2double t1329 = -0.3e1 * u_xz * t31 * t73 * t30 * t1024 * t1053 - 0.3e1 * t1251 * t29 * t1248
                        + 0.3e1 * u_0 + 0.3e1 * t58 + 0.3e1 * t62 + 0.3e1 * t66 + 0.3e1 * t70
                        + 0.3e1 * t74 + 0.3e1 * t79 + t1035 + t1247 + t32 * (t1322 + t1312)
                        + x * t111 * t448 + 0.9e1 * t1065 - 0.3e1 * t1068;
  const su2double t1348 = -0.3e1 / 0.2e1 * t1073 - 0.3e1 / 0.2e1 * t1197 - 0.9e1 * t1199 - 0.3e1 * t1202
                        + 0.3e1 * t1042 - 0.3e1 * t1046 + 0.3e1 / 0.2e1 * t1050 - 0.3e1 / 0.2e1 * t1078
                        - t1082 + 0.3e1 / 0.2e1 * t1180 - 0.3e1 * t1184 + 0.3e1 * t1188 + 0.3e1 / 0.2e1 * t1192
                        + t1206 - 0.3e1 * u_x * t31 * t57 * t1024 * t1318 + 0.9e1 * t1315;
  const su2double t1349 = t1348 + t1329;
  const su2double t1351 = t557 * t519;
  const su2double t1352 = t368 * t1351;
  const su2double t1356 = t22 + 0.32e2 * t1007 * t160 - 0.32e2 * t1012 * t160 - 0.16e2 * t1015 * t944 * t108
                        - 0.16e2 * t227 * t982 * t944 * t108 + 0.2e1 * t1123 * t866 * t1119
                        + 0.2e1 * t1123 * t337 * t377 * t1136 * t783 - 0.4e1 / 0.3e1 * t368 * t557 * t1151 * t820 * t819
                        + 0.32e2 * t1160 * t1159 * t1158 - 0.32e2 * t1160 * t1164 * t1158 + 0.2e1 * t1123 * t866 * t1256
                        + 0.2e1 * t1123 * t337 * t377 * t1269 * t670 + 0.32e2 * t1276 * t1158 - 0.32e2 * t1279 * t1158
                        - 0.4e1 / 0.3e1 * t1352 * muRef * t1349 * t819 - t6;
  const su2double t1360 = t36 * t592 * t685 * t676;
  const su2double t1372 = t379 - t381 / 0.2e1 - t382 / 0.2e1 + t384 + t385 + t386 + t387 + t388 + t389 - t390 / 0.2e1
                        - t392 / 0.2e1 - t394 / 0.2e1 - t396 / 0.2e1 - t398 / 0.2e1 - t400 / 0.2e1;
  const su2double t1385 = -t403 / 0.2e1 - t404 / 0.2e1 - t405 / 0.2e1 - t406 / 0.2e1 - t407 / 0.2e1 - t408 / 0.2e1
                        - t413 / 0.4e1 + t416 / 0.4e1 - t420 / 0.4e1 - t424 / 0.4e1 + t429 - t434
                        + t437 / 0.4e1 - t442 / 0.4e1 + t446;
  const su2double t1386 = t1385 + t1372;
  const su2double t1393 = -t380 + 0.3e1 / 0.4e1 * t381 + 0.3e1 / 0.4e1 * t382 + t32 * t1386 - t452 - t453 - t454
                        - t455 - t456 - t457 + 0.3e1 / 0.4e1 * t390 + 0.3e1 / 0.4e1 * t392 + 0.3e1 / 0.4e1 * t394
                        + 0.3e1 / 0.4e1 * t396 + 0.3e1 / 0.4e1 * t398;
  const su2double t1407 = 0.3e1 / 0.4e1 * t400 + 0.3e1 / 0.4e1 * t403 + 0.3e1 / 0.4e1 * t404 + 0.3e1 / 0.4e1 * t405
                        + 0.3e1 / 0.4e1 * t406 + 0.3e1 / 0.4e1 * t407 + 0.3e1 / 0.4e1 * t408 + 0.3e1 / 0.4e1 * t473
                        + 0.3e1 / 0.4e1 * t477 - 0.3e1 / 0.4e1 * t481 + 0.3e1 / 0.4e1 * t485 - t490 + t494
                        + 0.3e1 / 0.4e1 * t497 - 0.3e1 / 0.4e1 * t501 - t506;
  const su2double t1408 = t1407 + t1393;
  const su2double t1416 = t519 * t374;
  const su2double t1417 = t377 * t1416;
  const su2double t1419 = muRef * t373 * t1417;
  const su2double t1421 = x * v_xz * a_vxz;
  const su2double t1424 = y * t239;
  const su2double t1428 = x * w_xy * a_wxy;
  const su2double t1431 = z * t322;
  const su2double t1434 = z * t196;
  const su2double t1435 = z * t200;
  const su2double t1436 = z * t204;
  const su2double t1437 = y * t279;
  const su2double t1438 = y * t283;
  const su2double t1439 = y * t287;
  const su2double t1440 = t410 * v_z;
  const su2double t1441 = cos(t190);
  const su2double t1445 = t410 * w_y;
  const su2double t1446 = sin(t269);
  const su2double t1450 = y * t267;
  const su2double t1451 = y * t271;
  const su2double t1452 = y * t275;
  const su2double t1453 = z * t184;
  const su2double t1454 = z * t188;
  const su2double t1455 = z * t192;
  const su2double t1456 = w_0 * y;
  const su2double t1457 = v_0 * z;
  const su2double t1458 = t607 * t1421 / 0.2e1 + t427 * t1424 / 0.2e1 + t726 * t1428 / 0.2e1
                        + t423 * t1431 / 0.2e1 + t1434 + t1435 + t1436 + t1437 + t1438 + t1439
                        + t1441 * t1440 * t189 / 0.2e1 - t1446 * t1445 * t268 / 0.2e1 + t1450
                        + t1451 + t1452 + t1453 + t1454 + t1455 + t1456 + t1457;
  const su2double t1477 = t31 * w_y;
  const su2double t1480 = -v_z * t189 * t1441 * t31 + t1446 * t1477 * t268 - t1439 - t1450 - t1451
                        - t1452 - t1453 - t1454 - t1455 - t1456 - t1457;
  const su2double t1481 = 0.3e1 / 0.2e1 * t1480 + t32 * t1458 - 0.3e1 / 0.2e1 * t1421 * t649
                        - 0.3e1 / 0.2e1 * t1424 * t488 - 0.3e1 / 0.2e1 * t764 * t1428
                        - 0.3e1 / 0.2e1 * t476 * t1431 - 0.3e1 / 0.2e1 * t1434 - 0.3e1 / 0.2e1 * t1435
                        - 0.3e1 / 0.2e1 * t1436 - 0.3e1 / 0.2e1 * t1437 - 0.3e1 / 0.2e1 * t1438;
  const su2double t1483 = t560 * t376 * t1481;
  const su2double t1490 = -t896 * x * t358 + x * t899 + t833 - t837 + t840 - t844 + t845 - t849 + t851 - t854 + t859;
  const su2double t1492 = t1490 * t872 * t869;
  const su2double t1495 = muRef * t1408;
  const su2double t1497 = t32 * t358;
  const su2double t1500 = -t1497 * t901 + y * t856 + t880 - t882 + t885 - t887 + t888 - t890 + t892 - t894 + t903;
  const su2double t1505 = t371 * t672;
  const su2double t1506 = t1505 * t1416;
  const su2double t1507 = t1481 * muRef;
  const su2double t1508 = TWall * t1507;
  const su2double t1518 = t205 * t159 * t104;
  const su2double t1523 = t1136 * t374;
  const su2double t1525 = t1122 * t1507;
  const su2double t1528 = a_wyz * a_wyz;
  const su2double t1530 = t1024 * t1528 * t31;
  const su2double t1531 = t286 * z;
  const su2double t1534 = w_yz * y * t1531 * t1530;
  const su2double t1538 = t1441 * z * v_z * t189;
  const su2double t1540 = Pi * t29;
  const su2double t1542 = t323 * t1540 * t322;
  const su2double t1543 = 0.3e1 / 0.2e1 * t1542;
  const su2double t1546 = z * w_y * t1446 * t268;
  const su2double t1550 = t315 * y * w_z * t272;
  const su2double t1551 = 0.3e1 / 0.2e1 * t1550;
  const su2double t1552 = t111 * t1458;
  const su2double t1554 = t606 * t171;
  const su2double t1555 = t1554 * t605;
  const su2double t1556 = 0.9e1 / 0.2e1 * t1555;
  const su2double t1560 = y * v_yz * t240 * z * t201;
  const su2double t1561 = 0.9e1 / 0.2e1 * t1560;
  const su2double t1563 = t319 * t8 * t418;
  const su2double t1564 = 0.3e1 / 0.2e1 * t1563;
  const su2double t1566 = t1072 * x * t276;
  const su2double t1570 = w_yz * t323 * t30 * t284;
  const su2double t1572 = 0.3e1 / 0.2e1 * t1534 - 0.9e1 / 0.2e1 * t1538 - t1543 + 0.3e1 * t1546 + t1551
                        + z * t1552 - t1556 - t1561 - t1564 - 0.3e1 * t1566 - 0.3e1 * t1570;
  const su2double t1574 = w_yz * t284 * t475;
  const su2double t1576 = a_vxz * a_vxz;
  const su2double t1578 = t1024 * t1576 * t31;
  const su2double t1583 = a_vyz * a_vyz;
  const su2double t1585 = t1024 * t1583 * t31;
  const su2double t1590 = 0.3e1 / 0.2e1 * v_0;
  const su2double t1591 = 0.3e1 / 0.2e1 * t184;
  const su2double t1592 = 0.3e1 / 0.2e1 * t188;
  const su2double t1593 = 0.3e1 / 0.2e1 * t192;
  const su2double t1594 = 0.3e1 / 0.2e1 * t196;
  const su2double t1595 = 0.3e1 / 0.2e1 * t200;
  const su2double t1596 = 0.3e1 / 0.2e1 * t204;
  const su2double t1597 = 0.2e1 * t1555;
  const su2double t1600 = t199 * t1093;
  const su2double t1603 = 0.2e1 * t1560;
  const su2double t1604 = t1583 * v_yz;
  const su2double t1606 = t203 * t1093;
  const su2double t1609 = t423 * t322;
  const su2double t1611 = t1528 * w_yz;
  const su2double t1614 = t1531 * t1093 * y * t1611;
  const su2double t1618 = a_vz * a_vz;
  const su2double t1623 = t204 + t1563 + t1542 + 0.2e1 * t1538 - t191 * t1440 * t1024 * t1618 / 0.2e1
                        - t1546 - t1550 + t184 + t188 + t192 + v_0;
  const su2double t1631 = -0.3e1 / 0.2e1 * t1574 + 0.3e1 / 0.2e1 * v_xz * t199 * t28 * t1578
                        + 0.3e1 / 0.2e1 * v_yz * t203 * t29 * t1585 - t1590 - t1591 - t1592
                        - t1593 - t1594 - t1595 - t1596 + t32 * (t1623 + t1597
                        - t1600 * t1576 * v_xz * t28 / 0.2e1 + t1603 - t1606 * t29 * t1604 / 0.2e1
                        + t1566 + t1609 / 0.2e1 + t1570 - t1614 / 0.2e1 + t196 + t200)
                        + 0.3e1 / 0.2e1 * v_z * t191 * t1024 * t1618 * t31;
  const su2double t1632 = t1631 + t1572;
  const su2double t1637 = x * t247;
  const su2double t1646 = t610 * t1250 * y * t1248;
  const su2double t1652 = a_vx * a_vx;
  const su2double t1659 = t148 * t67 * t1218;
  const su2double t1663 = u_xy * t148 * t29 * t67;
  const su2double t1664 = 0.3e1 / 0.2e1 * t1663;
  const su2double t1667 = t624 * x * u_y * t59;
  const su2double t1671 = t630 * x * a_vx * t628;
  const su2double t1673 = u_x * t144;
  const su2double t1675 = y * t1673 * t55;
  const su2double t1676 = 0.3e1 / 0.2e1 * t1675;
  const su2double t1678 = t600 * t171;
  const su2double t1679 = t1678 * t599;
  const su2double t1681 = 0.3e1 / 0.2e1 * t1646 + 0.3e1 / 0.2e1 * v_xy * t195 * t29 * t1244
                        + 0.3e1 / 0.2e1 * v_x * t183 * t1024 * t1652 * t31 + 0.3e1 * t1659 + t1664
                        + 0.3e1 * t1667 + 0.9e1 / 0.2e1 * t1671 - t1676 + x * t1241 - 0.3e1 * t1679 - t1556;
  const su2double t1682 = t237 * t614;
  const su2double t1683 = 0.9e1 / 0.2e1 * t1682;
  const su2double t1687 = y * u_xz * t152 * z * t71;
  const su2double t1688 = 0.3e1 / 0.2e1 * t1687;
  const su2double t1694 = u_xy * t67 * t479;
  const su2double t1700 = t415 * t147;
  const su2double t1705 = t69 * y * t1093 * t1219 * t610;
  const su2double t1707 = 0.2e1 * t1682;
  const su2double t1710 = t195 * t1093;
  const su2double t1720 = t1687 + t200 + t204 - t1667 - 0.2e1 * t1671 - t183 * t410 * t1652 * v_x * t1024 / 0.2e1
                        + t184 + t188 + t192 + t1675 + v_0;
  const su2double t1723 = t1683 - t1688 + 0.3e1 / 0.2e1 * v_xz * t199 * t30 * t1578 + 0.3e1 / 0.2e1 * t1694
                        - t1590 - t1591 - t1592 - t1593 - t1594 - t1595 - t1596 + t32 * (t1720 + t1679
                        + t1597 - t1600 * t30 * t1576 * v_xz / 0.2e1 - t1700 / 0.2e1 - t1659 - t1705 / 0.2e1
                        - t1707 - t1710 * t29 * t1226 * v_xy / 0.2e1 - t1663 + t196);
  const su2double t1725 = (t1723 + t1681) * t374;
  const su2double t1730 = t337 * t377 * t1151;
  const su2double t1734 = t159 * t25;
  const su2double t1737 = t1360 * t674 + 0.4e1 / 0.3e1 * t36 * t711 * t560 * t558 * t1408 * t378 * t375
                        + t813 * t1483 * t1419 - 0.2e1 * t1492 * t867 - 0.8e1 / 0.3e1 * t1500 * t824 * t519 * t1495 * t819
                        - 0.2e1 * t940 * t1508 * t1506 + 0.8e1 / 0.3e1 * t368 * t557 * t1269 * t1495 * t819
                        + 0.32e2 * t1276 * t1518 - 0.32e2 * t1279 * t1518 + 0.2e1 * t1525 * t1505 * t1523
                        + 0.2e1 * t1122 * t1632 * muRef * t1506 + 0.32e2 * t1637 * t1159 * t1158
                        - 0.32e2 * t1637 * t1164 * t1158 + 0.2e1 * t1123 * t866 * t1725
                        + 0.2e1 * t1123 * t1730 * t670 + 0.16e2 * t983 * t1734;
  const su2double t1738 = t34 * t155;
  const su2double t1749 = -t630 * Pi * v_x * a_vx + t606 * t16 * t604 - t236 * t8 * t235;
  const su2double t1750 = t1749 * t945;
  const su2double t1754 = t205 * t205;
  const su2double t1755 = t34 * t1754;
  const su2double t1760 = t34 * t205;
  const su2double t1761 = t946 * t1760;
  const su2double t1764 = t951 * t1760;
  const su2double t1773 = t1441 * Pi * v_z * a_vz + t240 * t8 * t239 + t1554 * t604;
  const su2double t1792 = t33 * t1754;
  const su2double t1793 = y * t1005;
  const su2double t1794 = t1793 * t1792;
  const su2double t1797 = y * t1010;
  const su2double t1798 = t1797 * t1755;
  const su2double t1811 = 0.3e1 / 0.2e1 * v_xy * t195 * t28 * t1244 + 0.3e1 / 0.4e1 * t1679 + t1683
                        + t1688 - t1561 + t1564 + 0.3e1 / 0.4e1 * t1566 + 0.3e1 / 0.2e1 * v_yz * t203 * t30 * t1585
                        - t1590 - t1591 - t1592 - t1593 - t1594 - t1595 - t1596;
  const su2double t1822 = -t1679 / 0.2e1 - t1707 - t1687 / 0.2e1 + t1603 - t1563 / 0.2e1 - t1566 / 0.2e1
                        - t1710 * t1226 * t1186 / 0.2e1 - t1606 * t30 * t1604 / 0.2e1 + v_0 + t184
                        + t188 + t192 + t196 + t200;
  const su2double t1825 = v_y * t232;
  const su2double t1827 = y * t1825 * t185;
  const su2double t1837 = a_vy * a_vy;
  const su2double t1844 = t204 + t1705 / 0.4e1 + t1614 / 0.4e1 + 0.2e1 * t1827 - t1542 / 0.2e1 + t1546 / 0.2e1
                        + t1659 / 0.2e1 + t1663 / 0.2e1 + t1667 / 0.2e1 - t1675 / 0.2e1 + t1550 / 0.2e1
                        - t1570 / 0.2e1 - t187 * t443 * t1024 * t1837 / 0.2e1 + t1700 / 0.4e1 - t1609 / 0.4e1;
  const su2double t1863 = t32 * (t1844 + t1822) + y * t111 * t1386 - 0.3e1 / 0.4e1 * t1646 - 0.3e1 / 0.4e1 * t1534
                        + 0.3e1 / 0.2e1 * v_y * t187 * t1024 * t1837 * t31 - 0.9e1 / 0.2e1 * t1827 + t1543
                        - 0.3e1 / 0.4e1 * t1546 - 0.3e1 / 0.4e1 * t1659 - t1664 - 0.3e1 / 0.4e1 * t1667
                        + t1676 - t1551 - 0.3e1 / 0.4e1 * t1694 + 0.3e1 / 0.4e1 * t1570 + 0.3e1 / 0.4e1 * t1574;
  const su2double t1864 = t1863 + t1811;
  const su2double t1869 = 0.16e2 * t982 * t1738 * t160 + 0.16e2 * t1750 * t944 * t160
                        - 0.16e2 * t227 * t945 * t1755 * t108 + 0.32e2 * t1761 * t160 + 0.16e2 * t1764 * t950
                        + 0.16e2 * t951 * t34 * t1773 * t160 + 0.16e2 * t967 * t1760 * t160
                        + 0.16e2 * t999 * t1754 * t54 * t180 + t173 - t177 - 0.16e2 * t1015 * t1760 * t108
                        - 0.16e2 * t138 * t982 * t944 * t108 + 0.32e2 * t1794 * t160 - 0.32e2 * t1798 * t160
                        + 0.8e1 / 0.3e1 * t1352 * muRef * t1864 * t819 + t170;
  const su2double t1877 = t379 - 0.2e1 * t381 + t382 + t384 + t385 + t386 + t387 + t388 + t389 + t390
                        + t392 + t394 + t396 + t398 + t400;
  const su2double t1887 = -0.2e1 * t403 - 0.2e1 * t404 - 0.2e1 * t405 - 0.2e1 * t406 - 0.2e1 * t407
                        - 0.2e1 * t408 + t413 / 0.2e1 - t416 / 0.2e1 - t420 - t424 + t429
                        - t434 + t437 + t442 / 0.2e1 + t446;
  const su2double t1888 = t1887 + t1877;
  const su2double t1895 = -t380 + 0.3e1 * t381 - 0.3e1 / 0.2e1 * t382 + t32 * t1888 - t452 - t453 - t454
                        - t455 - t456 - t457 - 0.3e1 / 0.2e1 * t390 - 0.3e1 / 0.2e1 * t392
                        - 0.3e1 / 0.2e1 * t394 - 0.3e1 / 0.2e1 * t396 - 0.3e1 / 0.2e1 * t398;
  const su2double t1909 = -0.3e1 / 0.2e1 * t400 + 0.3e1 * t403 + 0.3e1 * t404 + 0.3e1 * t405 + 0.3e1 * t406
                        + 0.3e1 * t407 + 0.3e1 * t408 + 0.3e1 * t473 + 0.3e1 * t477 + 0.3e1 / 0.2e1 * t481
                        - 0.3e1 / 0.2e1 * t485 - t490 + t494 - 0.3e1 / 0.2e1 * t497 - 0.3e1 * t501 - t506;
  const su2double t1910 = t1909 + t1895;
  const su2double t1923 = muRef * t1910;
  const su2double t1927 = -t1497 * t936 + z * t856 + t917 - t919 + t922 - t924 + t925 - t927 + t929 - t931 + t938;
  const su2double t1932 = x * t330;
  const su2double t1940 = a_wx * a_wx;
  const su2double t1950 = t1678 * t720;
  const su2double t1952 = t725 * t171;
  const su2double t1953 = t1952 * t724;
  const su2double t1954 = 0.9e1 / 0.2e1 * t1953;
  const su2double t1955 = t320 * t729;
  const su2double t1960 = z * u_xy * t148 * y * t67;
  const su2double t1964 = u_xz * t152 * t30 * t71;
  const su2double t1967 = 0.3e1 / 0.2e1 * t484 * t151;
  const su2double t1968 = a_wxy * a_wxy;
  const su2double t1975 = t1053 * t716;
  const su2double t1976 = t73 * z;
  const su2double t1979 = 0.3e1 / 0.2e1 * t1976 * t1055 * t1975;
  const su2double t1980 = x * t1115 + 0.3e1 / 0.2e1 * w_x * t266 * t1024 * t1940 * t31
                        + 0.3e1 / 0.2e1 * w_xz * t282 * t30 * t1031 - 0.3e1 * t1950 - t1954
                        - 0.9e1 / 0.2e1 * t1955 + 0.3e1 / 0.2e1 * t1960 - 0.3e1 / 0.2e1 * t1964
                        - t1967 + 0.3e1 / 0.2e1 * w_xy * t278 * t29 * t1024 * t1968 * t31 + t1979;
  const su2double t1981 = 0.3e1 / 0.2e1 * w_0;
  const su2double t1983 = t412 * t151 / 0.2e1;
  const su2double t1985 = t152 * t71 * t1052;
  const su2double t1988 = t1976 * t1093 * t1975 / 0.2e1;
  const su2double t1989 = 0.2e1 * t1953;
  const su2double t1992 = t278 * t1093;
  const su2double t1998 = t282 * t1093;
  const su2double t2004 = t739 * x * u_z * t63;
  const su2double t2007 = t745 * x * a_wx * t743;
  const su2double t2015 = z * t1673 * t55;
  const su2double t2016 = t279 + t283 + t287 - t2004 - 0.2e1 * t2007 - t266 * t410 * t1940 * w_x * t1024 / 0.2e1
                        + t267 + t271 + t275 + t2015 + w_0;
  const su2double t2023 = 0.3e1 / 0.2e1 * t267;
  const su2double t2024 = 0.3e1 / 0.2e1 * t271;
  const su2double t2025 = 0.3e1 / 0.2e1 * t275;
  const su2double t2026 = 0.3e1 / 0.2e1 * t279;
  const su2double t2027 = 0.3e1 / 0.2e1 * t283;
  const su2double t2028 = 0.3e1 / 0.2e1 * t287;
  const su2double t2029 = -t1981 + t32 * (t2016 + t1983 + t1985 - t1988 + t1950 + t1989
                        - t1992 * t29 * t1968 * w_xy / 0.2e1 + 0.2e1 * t1955 - t1998 * t30 * t1029 * w_xz / 0.2e1
                        - t1960 + t1964) - 0.3e1 * t1985 + 0.3e1 * t2004 + 0.9e1 / 0.2e1 * t2007 - 0.3e1 / 0.2e1 * t2015
                        - t2023 - t2024 - t2025 - t2026 - t2027 - t2028;
  const su2double t2031 = (t2029 + t1980) * t374;
  const su2double t2038 = y * t330;
  const su2double t2045 = t1269 * t374;
  const su2double t2050 = t203 * z;
  const su2double t2054 = 0.3e1 / 0.2e1 * v_yz * y * t2050 * t1585;
  const su2double t2055 = a_wy * a_wy;
  const su2double t2056 = t1024 * t2055;
  const su2double t2061 = t606 * t8 * t1421;
  const su2double t2063 = t324 * t1431;
  const su2double t2068 = z * v_xy * t236 * x * t193;
  const su2double t2072 = v_yz * t240 * t30 * t201;
  const su2double t2076 = 0.3e1 / 0.2e1 * v_yz * t201 * t487;
  const su2double t2078 = t1968 * w_xy * t28;
  const su2double t2082 = t30 * t1611;
  const su2double t2086 = y * t1552 + t2054 + 0.3e1 / 0.2e1 * t270 * t1477 * t2056 - 0.3e1 * t2061
                        - t1954 - 0.9e1 / 0.2e1 * t2063 + 0.3e1 / 0.2e1 * t2068 - 0.3e1 / 0.2e1 * t2072
                        - t2076 + 0.3e1 / 0.2e1 * t278 * t1055 * t2078 + 0.3e1 / 0.2e1 * t286 * t1055 * t2082;
  const su2double t2088 = t427 * t239 / 0.2e1;
  const su2double t2090 = t240 * t1540 * t239;
  const su2double t2094 = t2050 * t1093 * y * t1604 / 0.2e1;
  const su2double t2098 = t286 * t1093;
  const su2double t2104 = t1441 * y * v_z * t189;
  const su2double t2107 = t1446 * y * w_y * t268;
  const su2double t2113 = z * t1825 * t185;
  const su2double t2114 = t279 + t283 + t287 + t2104 - 0.2e1 * t2107 - t270 * t1445 * t2056 / 0.2e1
                        + t267 + t271 + t275 + t2113 + w_0;
  const su2double t2121 = -t1981 + t32 * (t2114 + t2061 + t2088 + t2090 - t2094 + t1989
                        - t1992 * t2078 / 0.2e1 + 0.2e1 * t2063 - t2098 * t2082 / 0.2e1 - t2068 + t2072)
                        - 0.3e1 * t2090 - 0.3e1 * t2104 + 0.9e1 / 0.2e1 * t2107 - 0.3e1 / 0.2e1 * t2113
                        - t2023 - t2024 - t2025 - t2026 - t2027 - t2028;
  const su2double t2122 = t2121 + t2086;
  const su2double t2148 = -0.3e1 / 0.2e1 * t1950 + 0.9e1 * t1955 + 0.3e1 * t1960 - 0.3e1 / 0.2e1 * t2061
                        + 0.9e1 * t2063 + 0.3e1 * t2068 - 0.3e1 * w_xz * t282 * t28 * t1031
                        - 0.3e1 * w_yz * t286 * t29 * t1530 + 0.3e1 * w_0 + 0.3e1 * t267 + 0.3e1 * t271
                        + 0.3e1 * t275 + 0.3e1 * t279 + 0.3e1 * t283 + 0.3e1 * t287;
  const su2double t2161 = t1998 * t1029 * t1044 + t2098 * t29 * t1611 + t1950 - 0.4e1 * t1955 - t1960
                        + t2061 - 0.4e1 * t2063 - t2068 - 0.2e1 * t267 - 0.2e1 * t271 - 0.2e1 * t275
                        - 0.2e1 * t279 - 0.2e1 * t283 - 0.2e1 * w_0;
  const su2double t2163 = a_wz * a_wz;
  const su2double t2169 = t315 * z * w_z * t272;
  const su2double t2171 = t274 * t435 * t1024 * t2163 + t1964 + t1983 + t1985 - t1988 - t2004 + t2015
                        + t2072 + t2088 + t2090 - t2094 + t2104 + t2113 + 0.4e1 * t2169 - 0.2e1 * t287;
  const su2double t2190 = t32 * (t2171 + t2161) + z * t111 * t1888 - 0.3e1 * t1964 - t1967
                        - 0.3e1 / 0.2e1 * t1985 + 0.3e1 / 0.2e1 * t2004 - 0.3e1 * t2015 - 0.3e1 * t2072
                        - t2076 - 0.3e1 / 0.2e1 * t2090 - 0.3e1 / 0.2e1 * t2104 + t1979 + t2054
                        - 0.3e1 * t2113 - 0.3e1 * w_z * t274 * t1024 * t2163 * t31 - 0.9e1 * t2169;
  const su2double t2191 = t2190 + t2148;
  const su2double t2201 = t1360 * t785 + t713 * t1483 * t1419 - 0.2e1 / 0.3e1 * t36 * t811 * t560 * t558 * t1910 * t378 * t375
                        - 0.2e1 * t1492 * t909 - 0.2e1 * t905 * t1508 * t1506 + 0.4e1 / 0.3e1 * t1927 * t824 * t519 * t1923 * t819
                        + 0.32e2 * t1932 * t1159 * t1158 - 0.32e2 * t1932 * t1164 * t1158 + 0.2e1 * t1123 * t866 * t2031
                        + 0.2e1 * t1123 * t1730 * t783 + 0.32e2 * t2038 * t1159 * t1518 - 0.32e2 * t2038 * t1164 * t1518
                        + 0.2e1 * t1525 * t1505 * t2045 + 0.2e1 * t1122 * t2122 * muRef * t1506
                        - 0.4e1 / 0.3e1 * t1352 * muRef * t2191 * t819 - 0.4e1 / 0.3e1 * t368 * t557 * t1136 * t1923 * t819;
  const su2double t2214 = -t745 * Pi * w_x * a_wx + t319 * t16 * t318 + t725 * t8 * t723;
  const su2double t2215 = t2214 * t945;
  const su2double t2231 = -t1446 * Pi * w_y * a_wy + t323 * t16 * t322 + t1952 * t723;
  const su2double t2232 = t2231 * t945;
  const su2double t2236 = t288 * t288;
  const su2double t2237 = t34 * t2236;
  const su2double t2242 = t34 * t288;
  const su2double t2243 = t967 * t2242;
  const su2double t2258 = t33 * t2236;
  const su2double t2259 = z * t1005;
  const su2double t2260 = t2259 * t2258;
  const su2double t2263 = z * t1010;
  const su2double t2264 = t2263 * t2237;
  const su2double t2267 = -0.16e2 * t138 * t951 * t944 * t108 - 0.16e2 * t227 * t951 * t1760 * t108
                        - 0.16e2 * t310 * t945 * t2237 * t108 + 0.16e2 * t951 * t34 * t243 * t160
                        + 0.16e2 * t999 * t2236 * t54 * t263 + 0.16e2 * t951 * t1738 * t160
                        + 0.16e2 * t2232 * t1760 * t160 + 0.16e2 * t2215 * t944 * t160 + 0.32e2 * t2243 * t160
                        + 0.32e2 * t2260 * t160 - 0.32e2 * t2264 * t160 + 0.16e2 * t952 * t1734
                        + 0.16e2 * t1764 * t981 - t256 + t258 - t260;
  const su2double t2269 = t374 * t816;
  const su2double t2270 = t377 * TWall;
  const su2double t2274 = t82 * t288 * t368;
  const su2double t2279 = t1910 * t2270 * t2269;
  const su2double t2293 = t82 * t80 * t368;
  const su2double t2298 = t508 * t2270 * t2269;
  const su2double t2312 = t368 * TWall * muRef;
  const su2double t2313 = t205 * t685;
  const su2double t2314 = t82 * t2313;
  const su2double t2315 = t2314 * t2312;
  const su2double t2322 = t1505 * t671;
  const su2double t2330 = t288 * t685;
  const su2double t2331 = t82 * t2330;
  const su2double t2332 = t2331 * t2312;
  const su2double t2339 = -0.16e2 / 0.3e1 * t2274 * t558 * t2191 * t2270 * t2269 - 0.16e2 / 0.3e1 * t2274 * t557 * t1136 * muRef * t2279
                        - 0.16e2 / 0.3e1 * t82 * t326 * t368 * t558 * t2279 - 0.16e2 / 0.3e1 * t2293 * t558 * t1349 * t2270 * t2269
                        - 0.16e2 / 0.3e1 * t2293 * t557 * t1151 * muRef * t2298 - 0.16e2 / 0.3e1 * t82 * t155 * t368 * t558 * t2298
                        + 0.8e1 * t2315 * t1505 * t519 * t1725 + 0.8e1 * t2315 * t1505 * t1151 * t670
                        + 0.8e1 * t82 * t1749 * t685 * t2312 * t2322 + 0.8e1 * t2332 * t1505 * t519 * t2031
                        + 0.8e1 * t2332 * t1505 * t1151 * t783;
  const su2double t2340 = t1505 * t784;
  const su2double t2348 = t80 * t685;
  const su2double t2350 = t82 * t2348 * t2312;
  const su2double t2365 = t82 * t205 * t368;
  const su2double t2370 = t1408 * t2270 * t2269;
  const su2double t2382 = t117 * t33;
  const su2double t2385 = t2382 * t130 - a_T3 / 0.4e1;
  const su2double t2391 = 0.3e1 / 0.8e1 * t2382 * t125 - a_T2 / 0.8e1 + a_T4 / 0.2e1;
  const su2double t2393 = t349 * Pi;
  const su2double t2398 = t2382 * t2393 / 0.8e1 - a_T1 / 0.16e2 + 0.3e1 / 0.16e2 * a_T3;
  const su2double t2400 = t354 * Pi;
  const su2double t2405 = -t528 / 0.2e1 + t345 * t2385 + t341 * t2391 + t40 * t2398 + t2382 * t2400 / 0.32e2 + a_T2 / 0.16e2 - a_T4 / 0.16e2;
  const su2double t2406 = t2405 * t111;
  const su2double t2408 = 0.1e1 / t35 / 0.2e1;
  const su2double t2409 = t2408 * t33;
  const su2double t2412 = t53 * TWall / 0.2e1;
  const su2double t2413 = t368 * t2412;
  const su2double t2414 = sqrt(t2413);
  const su2double t2416 = t374 * t2414 * t2413;
  const su2double t2417 = t2412 + S;
  const su2double t2418 = 0.1e1 / t2417;
  const su2double t2420 = Gamma * RGas;
  const su2double t2422 = 1 / (Gamma - 1);
  const su2double t2423 = 1 / Prandtl;
  const su2double t2424 = t2423 * t2422;
  const su2double t2426 = t2424 * t2420 * t2418 * t2416;
  const su2double t2429 = 0.1e1 / t31;
  const su2double t2440 = muRef * t2408 * TWall;
  const su2double t2449 = 0.8e1 * t82 * t2214 * t685 * t2312 * t2340 + 0.8e1 * t2350 * t1505 * t519 * t1256
                        + 0.8e1 * t2350 * t1505 * t1269 * t670 + 0.8e1 * t82 * t993 * t685 * t2312 * t2322
                        + 0.32e2 / 0.3e1 * t2365 * t558 * t1864 * t2270 * t2269 + 0.32e2 / 0.3e1 * t2365 * t557 * t1269 * muRef * t2370
                        + 0.32e2 / 0.3e1 * t82 * t243 * t368 * t558 * t2370 + 0.96e2 * t2426 * muRef * t2409 * TWall * t2406
                        + 0.32e2 * t2426 * muRef * t2408 * t29 * TWall * t2405 * t2429 + 0.32e2 * t2426 * t2440 * t2405 * t2429 * t28
                        + 0.8e1 * t82 * t962 * t685 * t2312 * t2340;
  const su2double t2452 = muRef * t816;
  const su2double t2455 = t368 * TWall * t1481;
  const su2double t2459 = t2452 * t1417;
  const su2double t2498 = t376 * t2405;
  const su2double t2501 = t2414 * muRef * t2409;
  const su2double t2504 = t2420 * t2418 * t374;
  const su2double t2518 = 0.8e1 * t2314 * t2455 * t2452 * t377 * t1523 + 0.8e1 * t2314 * t368 * TWall * t1632 * t2459
                        + 0.8e1 * t82 * t1773 * t685 * t2455 * t2459 + 0.32e2 * t2426 * t2440 * t2405 * t2429 * t30
                        + 0.8e1 * t2331 * t2455 * t2452 * t377 * t2045 + 0.8e1 * t2331 * t368 * TWall * t2122 * t2459
                        + 0.8e1 * t82 * t2231 * t685 * t2455 * t2459 + 0.8e1 * t2350 * t1505 * t519 * t1119
                        + 0.8e1 * t2350 * t1505 * t1136 * t783 + 0.48e2 * t368 * t138 * t2424 * t2504 * t2501 * t2498 * t116
                        + 0.48e2 * t368 * t227 * t2424 * t2504 * t2501 * y * t376 * t2406;
  const su2double t2526 = Pi * t834;
  const su2double t2531 = t33 * a_T4 * t1024;
  const su2double t2532 = t40 * t116;
  const su2double t2537 = Pi * t341 * t2385;
  const su2double t2542 = t33 * a_T3 * t1024;
  const su2double t2547 = Pi * t40 * t2391;
  const su2double t2552 = t33 * t349 * t1024;
  const su2double t2556 = Pi * t2398;
  const su2double t2561 = t33 * t354 * t1024;
  const su2double t2567 = muRef * t2408;
  const su2double t2568 = t2567 * t33 * TWall;
  const su2double t2573 = t1923 * t2270 * t2269;
  const su2double t2578 = 0.1e1 / t32 / t31;
  const su2double t2585 = t820 * t2270 * t2269;
  const su2double t2587 = t368 * t823 * t519;
  const su2double t2599 = t868 * t672;
  const su2double t2600 = t2599 * t671;
  const su2double t2601 = t871 * t1121;
  const su2double t2602 = -t1490 * t38;
  const su2double t2615 = t2599 * t784;
  const su2double t2624 = 0.48e2 * t368 * t310 * t2424 * t2504 * t2501 * t2498 * t295
                        + 0.32e2 * t2426 * t2568 * (0.2e1 * t118 * t2526 + t345 * (t118 * t130 + t2532 * t2531)
                        - 0.3e1 * t118 * t2537 + 0.3e1 / 0.8e1 * t341 * (t118 * t125 + t2532 * t2542) - 0.2e1 * t118 * t2547
                        + t40 * (t118 * t2393 + t2532 * t2552) / 0.8e1 - t118 * t2556 + t118 * t2400 / 0.32e2
                        + t2532 * t2561 / 0.32e2) * t116 + 0.16e2 / 0.3e1 * t334 * t1352 * t2573
                        - 0.32e2 * t2426 * t2568 * t2405 * t2578 * t28 - 0.16e2 / 0.3e1 * t860 * t38 * t109 * t2587 * t2585
                        - 0.16e2 / 0.3e1 * t163 * t1352 * t2585 + 0.16e2 / 0.3e1 * t165 * t1352 * t2585
                        - 0.8e1 * t2602 * t209 * t2601 * t2600 + 0.8e1 * t162 * t247 * t1122 * t2600 - 0.8e1 * t116 * t209 * t1122 * t2600
                        - 0.8e1 * t2602 * t292 * t2601 * t2615 + 0.8e1 * t162 * t330 * t1122 * t2615;
  const su2double t2634 = t40 * t212;
  const su2double t2659 = t2567 * t33 * y;
  const su2double t2670 = -t904 * t38;
  const su2double t2684 = t1495 * t2270 * t2269;
  const su2double t2697 = t1507 * t816 * t1417;
  const su2double t2706 = -0.8e1 * t116 * t292 * t1122 * t2615 + 0.32e2 * t2426 * t2659 * TWall * (0.2e1 * t213 * t2526
                        + t345 * (t213 * t130 + t2634 * t2531) - 0.3e1 * t213 * t2537
                        + 0.3e1 / 0.8e1 * t341 * (t213 * t125 + t2634 * t2542) - 0.2e1 * t213 * t2547 + t40 * (t213 * t2393
                        + t2634 * t2552) / 0.8e1 - t213 * t2556 + t213 * t2400 / 0.32e2 + t2634 * t2561 / 0.32e2) * t111
                        - 0.32e2 * t2426 * t2567 * t33 * t29 * TWall * t2405 * t2578 - 0.8e1 * t2670 * t109 * t2601 * t2600
                        + 0.8e1 * t248 * t161 * t1122 * t2600 - 0.8e1 * t212 * t109 * t1122 * t2600
                        + 0.32e2 / 0.3e1 * t1500 * t38 * t209 * t2587 * t2684 + 0.32e2 / 0.3e1 * t249 * t1352 * t2684
                        - 0.32e2 / 0.3e1 * t251 * t1352 * t2684 - 0.8e1 * t2670 * t292 * t2601 * t2697
                        + 0.8e1 * t248 * t330 * t1122 * t2697;
  const su2double t2714 = t40 * t295;
  const su2double t2745 = -t939 * t38;
  const su2double t2780 = t2567 * t33 * t376;
  const su2double t2782 = t2417 * t2417;
  const su2double t2785 = Gamma / t2782 * t2416;
  const su2double t2786 = t2422 * RGas;
  const su2double t2792 = -0.8e1 * t212 * t292 * t1122 * t2697 + 0.32e2 * t2426 * t2568 * (0.2e1 * t296 * t2526
                        + t345 * (t296 * t130 + t2714 * t2531) - 0.3e1 * t296 * t2537 + 0.3e1 / 0.8e1 * t341 * (t296 * t125
                        + t2714 * t2542) - 0.2e1 * t296 * t2547 + t40 * (t296 * t2393 + t2714 * t2552) / 0.8e1 - t296 * t2556
                        + t296 * t2400 / 0.32e2 + t2714 * t2561 / 0.32e2) * t295 - 0.32e2 * t2426 * t2568 * t2405 * t2578 * t30
                        - 0.8e1 * t2745 * t109 * t2601 * t2615 + 0.8e1 * t331 * t161 * t1122 * t2615
                        - 0.8e1 * t295 * t109 * t1122 * t2615 - 0.8e1 * t2745 * t209 * t2601 * t2697
                        + 0.8e1 * t331 * t247 * t1122 * t2697 - 0.8e1 * t295 * t209 * t1122 * t2697
                        - 0.16e2 / 0.3e1 * t1927 * t38 * t292 * t2587 * t2573 - 0.16e2 / 0.3e1 * t332 * t1352 * t2573
                        - 0.32e2 * t138 * t2423 * t2786 * t2785 * t2780 * t2405 * t116;
  const su2double t2806 = t36 * t592 * t82;
  const su2double t2811 = muRef * t372 * t672;
  const su2double t2812 = t2811 * t671;
  const su2double t2813 = t560 * t376;
  const su2double t2814 = t2313 * t2813;
  const su2double t2818 = t2811 * t784;
  const su2double t2819 = t2330 * t2813;
  const su2double t2823 = t2348 * t2813;
  const su2double t2825 = t36 * t711 * t82;
  const su2double t2837 = t1507 * t373 * t1417;
  const su2double t2842 = t36 * t811 * t82;
  const su2double t2863 = -0.32e2 * t227 * t2423 * t2786 * t2785 * t2659 * t376 * t2406
                        - 0.8e1 / 0.3e1 * t2806 * t80 * t560 * t1351 * t820 * t378 * t375 + 0.4e1 * t2806 * t2814 * t2812
                        + 0.4e1 * t2806 * t2819 * t2818 + 0.4e1 * t2825 * t2823 * t2812
                        + 0.16e2 / 0.3e1 * t2825 * t205 * t560 * t1351 * t1495 * t378 * t375 + 0.4e1 * t2825 * t2819 * t2837
                        + 0.4e1 * t2842 * t2823 * t2818 + 0.4e1 * t2842 * t2814 * t2837
                        - 0.8e1 / 0.3e1 * t2842 * t288 * t560 * t1351 * t1923 * t378 * t375
                        - 0.32e2 * t310 * t2423 * t2786 * t2785 * t2780 * t2405 * t295;
  const su2double t2868 = t945 * t1755 + t945 * t2237 + t945 * t972;
  const su2double t2869 = 0.16e2 * t2868 * t159;
  const su2double t2872 = t2422 * t103 + t2869 * t104 / 0.2e1 + P_0 + t86 + t90 + t94 + t96 + t98 + t102;
  const su2double t2873 = t288 * t2872;
  const su2double t2875 = t111 * t33;
  const su2double t2878 = t205 * t2872;
  const su2double t2882 = t80 * t2872;
  const su2double t2889 = 0.16e2 * t2868 * t106;
  const su2double t2953 = t2263 * t2873 - z * t2875 * t2873 + t1797 * t2878 - y * t2875 * t2878 + t1011 * t2882
                        - x * t2875 * t2882 + t82 * t288 * (t2422 * t261 + t2869 * t262 / 0.2e1 - t310 * t2889 * t142 / 0.2e1
                        + 0.16e2 * (t1773 * t945 * t1760 + t962 * t945 * t944 + t2259 * t1004 - t2263 * t1755
                        + t2259 * t1792 - t2263 * t972 + t2243 + t2260 - t2264) * t159 * t104 - t256 + t258 - t260)
                        + t82 * t326 * t2872 + t82 * t155 * t2872 + t82 * t205 * (t2422 * t178 + t2869 * t179 / 0.2e1
                        - t227 * t2889 * t142 / 0.2e1 + 0.16e2 * (t993 * t945 * t944 + t1793 * t1004 + t1793 * t2258
                        - t1797 * t2237 - t1797 * t972 + t2232 * t2242 + t1761 + t1794 - t1798) * t159 * t104 + t170
                        + t173 - t177) + t82 * t243 * t2872 + t82 * t80 * (t2422 * t23 + t2869 * t25 / 0.2e1 - t138 * t2889 * t142 / 0.2e1
                        + 0.16e2 * (t1006 * t1792 + t1006 * t2258 - t1011 * t1755 - t1011 * t2237 + t1750 * t1760
                        + t2215 * t2242 + t1007 - t1012 + t978) * t159 * t104 - t6 + t14 + t22);

  /*--- Set the source term. Note the scaling for the correct non-dimensionalization. ---*/
  val_source[0] = (0.4e1 * t336)/(Density_Ref*Velocity_Ref);
  val_source[1] = (t1356 + t1003)/Pressure_Ref;
  val_source[2] = (t1869 + t1737)/Pressure_Ref;
  val_source[3] = (t2267 + t2201)/Pressure_Ref;
  val_source[4] = (0.4e1 * t2953 + t2863 + t2792 + t2706 + t2624 + t2518 + t2449 + t2339)/(Velocity_Ref*Pressure_Ref);
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
  val_source[1] = 4*Density*pow(u_0, 2)*x*(epsilon + sin(pow(x, 2) + pow(y, 2)))*cos(pow(x, 2) + pow(y, 2)) - 2*Density*u_0*v_0*y*(epsilon + sin(pow(x, 2) + pow(y, 2)))*sin(pow(x, 2) + pow(y, 2)) + 2*Density*u_0*v_0*y*(epsilon + cos(pow(x, 2) + pow(y, 2)))*cos(pow(x, 2) + pow(y, 2)) + 2*P_0*x*cos(pow(x, 2) + pow(y, 2)) - 0.666666666666667*Viscosity*(-8.0*u_0*pow(x, 2)*sin(pow(x, 2) + pow(y, 2)) + 4.0*u_0*cos(pow(x, 2) + pow(y, 2)) + 4*v_0*x*y*cos(pow(x, 2) + pow(y, 2))) + 4*u_0*pow(y, 2)*sin(pow(x, 2) + pow(y, 2)) - 2*u_0*cos(pow(x, 2) + pow(y, 2)) + 4*v_0*x*y*cos(pow(x, 2) + pow(y, 2));
  val_source[2] = -2*Density*u_0*v_0*x*(epsilon + sin(pow(x, 2) + pow(y, 2)))*sin(pow(x, 2) + pow(y, 2)) + 2*Density*u_0*v_0*x*(epsilon + cos(pow(x, 2) + pow(y, 2)))*cos(pow(x, 2) + pow(y, 2)) - 4*Density*pow(v_0, 2)*y*(epsilon + cos(pow(x, 2) + pow(y, 2)))*sin(pow(x, 2) + pow(y, 2)) + 2*P_0*y*cos(pow(x, 2) + pow(y, 2)) + 0.666666666666667*Viscosity*(-4*u_0*x*y*sin(pow(x, 2) + pow(y, 2)) + 8.0*v_0*pow(y, 2)*cos(pow(x, 2) + pow(y, 2)) + 4.0*v_0*sin(pow(x, 2) + pow(y, 2))) + 4*u_0*x*y*sin(pow(x, 2) + pow(y, 2)) + 4*v_0*pow(x, 2)*cos(pow(x, 2) + pow(y, 2)) + 2*v_0*sin(pow(x, 2) + pow(y, 2));
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
