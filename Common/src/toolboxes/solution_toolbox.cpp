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
  
  /*--- Some repeated terms for shorthand. ---*/
  const su2double sxy = sin(x*x + y*y);
  const su2double cxy = cos(x*x + y*y);
  
  /*--- The expressions for the source terms are given in
   the paper by Salari & Knupp. Note that the leading 1.0/rho
   term does not appear here, because our formulation allows
   for variable density (source should be multiplied by density). ---*/
  val_source[0]      = 2.0*u_0*x*cxy - 2.0*v_0*y*sxy;
  val_source[1]      = 2.0*((P_0*x + Density*u_0*(2.0*epsilon*u_0*x - 2.0*Viscosity + epsilon*v_0))*cxy + Density*u_0*(v_0*y*cos(2.0*(x*x + y*y)) + (2.0*x*x*Viscosity - epsilon*v_0*y + 2.0*Viscosity*y*y + 2.0*u_0*x*cxy)*sxy));
  
  val_source[2]      = 2.0*((epsilon*Density*u_0*v_0*x + 2.0*Density*v_0*x*x*Viscosity + P_0*y + 2.0*Density*v_0*Viscosity*y*y)*cxy + (-1.0*Density) *v_0*(-1.0*u_0*x*cos(2.0*(x*x + y*y)) + (epsilon*u_0*x - 2.0*Viscosity + 2.0*epsilon*v_0*y + 2.0*v_0*y*cxy)*sxy));
  
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
