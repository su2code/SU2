/*!
 * \file CInviscidShockVortexSolution.cpp
 * \brief Implementations of the member functions of CInviscidShockVortexSolution.
 * \author E. Molina, based on CInviscidVortexSolution
 * \version 7.0.3 "Blackbird"
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

#include "../../../include/toolboxes/MMS/CInviscidShockVortexSolution.hpp"

CInviscidShockVortexSolution::CInviscidShockVortexSolution(void) : CVerificationSolution() { }

CInviscidShockVortexSolution::CInviscidShockVortexSolution(unsigned short val_nDim,
                                                 unsigned short val_nVar,
                                                 unsigned short val_iMesh,
                                                 CConfig*       config)
  : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {

  /*--- Write a message that the solution is initialized for the
   inviscid vortex test case. ---*/
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the inviscid shock vortex interaction case!!!" << endl;
    cout << endl << flush;
  }

  /*--- Useful coefficients in which Gamma is present. ---*/
  Gamma    = config->GetGamma();
  Gm1      = Gamma - 1.0;
  ovGm1    = 1.0/Gm1;
  gamOvGm1 = ovGm1*Gamma;
  RGas     = config->GetGas_ConstantND();

  /*--- Vortex Location ---*/
  x_c   =  0.25;
  y_c   =  0.5;

  /*--- Vortex Size ---*/
  a = 0.075;
  b = 0.175;

  /*--- Vortex Strength ---*/
  M_v = 1.1;
  v_m = M_v * sqrt(Gamma);

  /*--- Shock conditions. ---*/
  M_s   = 1.1;

  /*--- Shock Initialization. ---*/
  /*--- Upstream conditions. ---*/
  rho_u = 1.0;
  u_u   = M_s * sqrt(Gamma);
  v_u   = 1e-20;
  p_u   = 1.0;
  t_u = p_u / (rho_u * RGas);

  /*--- Downstream conditions. ---*/
//  rho_d = rho_u * (Gamma + 1.0) * M_s * M_s / (2.0 + (Gamma - 1.0) * M_s * M_s);
//  u_d   = u_u * (2.0 + (Gamma - 1.0) * M_s * M_s) / ((Gamma + 1.0) * M_s * M_s);
//  v_d   = 1e-20;
//  p_d   = p_u * (1.0 + (2.0 * Gamma / (Gamma + 1.0)) * (M_s * M_s - 1.0));

  rho_d = rho_u*((Gamma+1)*M_s*M_s)/((Gamma-1)*M_s*M_s+2);
  u_d = (rho_u*u_u)/rho_d; // Conservation of mass
  v_d = 0.0;
  p_d = p_u*(2*Gamma*M_s*M_s-(Gamma-1))/(Gamma+1);

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if((config->GetTime_Marching() != TIME_STEPPING) &&
     (config->GetTime_Marching() != DT_STEPPING_1ST) &&
     (config->GetTime_Marching() != DT_STEPPING_2ND))
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

CInviscidShockVortexSolution::~CInviscidShockVortexSolution(void) { }

void CInviscidShockVortexSolution::GetBCState(const su2double *val_coords,
                                         const su2double val_t,
                                         su2double       *val_solution) const {

  /*--- For the case that the inviscid vortex is run with boundary
        conditions (other possibility is with periodic conditions),
        the exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

bool CInviscidShockVortexSolution::ExactSolutionKnown(void) const {return false;}

void CInviscidShockVortexSolution::GetSolution(const su2double *val_coords,
                                          const su2double val_t,
                                          su2double       *val_solution) const {
  su2double p, u, v, rho;
  su2double epsilon = 0.3;
  su2double rc = 0.05;
  su2double alph = 0.204;

  /*--- Distance from Vortex ---*/
  su2double dx = (val_coords[0] - x_c);
  su2double dy = (val_coords[1] - y_c);
  su2double r  = sqrt(dx * dx + dy * dy);
  su2double tau = r/rc;
  su2double theta = atan2(dy, dx);

  /*--- Shock conditions ---*/
  if (val_coords[0] <= 0.5){

    /*--- SUPERIMPOSE VORTEX ---*/
    rho = pow(1-(Gamma-1)*epsilon*epsilon*exp(2*alph*(1-tau*tau))/(4*alph*Gamma),(1/(Gamma-1)));
    p = pow(rho,Gamma);
    u = u_u + epsilon*tau*exp(alph*(1.-tau*tau))*sin(theta);
    v = v_u - epsilon*tau*exp(alph*(1.-tau*tau))*cos(theta);
  }
  else{
    rho = rho_d;
    p = p_d;
    u = u_d;
    v = v_d;
  }

  /* Compute the conservative variables. Note that both 2D and 3D
     cases are treated correctly. */
  val_solution[0]      = rho;
  val_solution[1]      = rho*u;
  val_solution[2]      = rho*v;
  val_solution[3]      = 0.0;
  val_solution[nVar-1] = ovGm1*p + 0.5*rho*(u*u + v*v);

}
