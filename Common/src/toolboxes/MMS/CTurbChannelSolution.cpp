/*!
 * \file CTurbChannelSolution.cpp
 * \brief Implementations of the member functions of CTurbChannelSolution
 * \author E.Molina
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

#include "../../../include/toolboxes/MMS/CTurbChannelSolution.hpp"
#include <random>

CTurbChannelSolution::CTurbChannelSolution(void) : CVerificationSolution() { }

CTurbChannelSolution::CTurbChannelSolution(unsigned short val_nDim,
                           unsigned short val_nVar,
                           unsigned short val_iMesh,
                           CConfig*       config)
  : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
    
  /*--- Store specific parameters here. ---*/
    
  ReynoldsFriction = 550.0; // Friction Reynolds Number.

  /*--- Turbulent flow. Use the relation of Malaspinas and Sagaut,JCP 275 2014,
  to compute the Reynolds number. ---*/
    
  ReynoldsMeanVelocity = pow((8./ 0.073),(4./7.))
                        * pow(ReynoldsFriction,(8./7.));
    
  /*--- Useful coefficients  ---*/
  RGas            = config->GetGas_ConstantND();
  su2double Gamma = config->GetGamma();
  su2double Cv    = RGas/(Gamma - 1.0);
  su2double Prandtl_Lam = config->GetPrandtl_Lam();
  su2double muLam = config->GetViscosity_FreeStreamND();
  su2double rhoDim = config->GetDensity_FreeStreamND();
    
  pRef  = config->GetPressure_FreeStreamND();
  ovGm1 = 1.0/(Gamma - 1.0);

  // Friction velocity
  su2double uTau = (2.0*ReynoldsFriction*muLam)/(rhoDim*config->GetLength_Reynolds());

  // Compute the wall shear stress and the body force.
  su2double tauWall = rhoDim*uTau*uTau;

  su2double fBodyX = 2.0*tauWall/config->GetLength_Reynolds();

  // The velocity profile is approximated with the following fit
  // u = a0*(1-|2*y/h|^alpha), where a0 = uMean*(alpha+1)/alpha.
  // The coefficient alpha can be computed from the friction and mean
  // velocity Reynolds numbers and uMean from the Mach number.
  // Compute these coefficients below.
  su2double uMean = config->GetModVel_FreeStreamND();
  alpha = 2.0*ReynoldsFriction*ReynoldsFriction/ReynoldsMeanVelocity - 1.0;
  a0    = uMean*(alpha+1.0)/alpha;

  // From a simple analysis of the fully developed equations it follows
  // that the reduced temperature distribution is theta = 1-|2*y/h|^(alpha+2).
  // When the removed heat equals the work of the body force vector, the
  // temperature in the middle of the channel can be estimated. This is
  // done below.
  halfChan = 0.5*config->GetLength_Reynolds();
  su2double kCondLam = muLam*Gamma*Cv/Prandtl_Lam;
  TWall    = 273.15; // Need to fix this
  TMiddle  = TWall + halfChan*halfChan*uMean*fBodyX/(kCondLam*(alpha+2.0));
  cout << muLam << " " << pRef << " " << kCondLam << " " << TMiddle <<endl;
    
  /*--- Write a message that the solution is initialized for the
   Turbulent Channel test case. ---*/
  
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the Turbulent Channel case!!!" << endl;
    cout << setprecision(9) << "Friction Reynolds: " << ReynoldsFriction << " Reynolds Number based on velocity: " << ReynoldsMeanVelocity << endl;
    cout << setprecision(9) << "Body Force: " << fBodyX << endl;
    cout << endl << flush;
  }

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

//  if(config->GetKind_ViscosityModel() != CONSTANT_VISCOSITY)
//    SU2_MPI::Error("Constant viscosity must be selected for the Taylor Green Vortex",
//                   CURRENT_FUNCTION);

  if(config->GetKind_ConductivityModel() != CONSTANT_PRANDTL)
    SU2_MPI::Error("Constant Prandtl number must be selected for the Taylor Green Vortex",
                   CURRENT_FUNCTION);
  
}

CTurbChannelSolution::~CTurbChannelSolution(void) { }

void CTurbChannelSolution::GetSolution(const su2double *val_coords,
                               const su2double val_t,
                               su2double       *val_solution) {
  
  /* The initial conditions are set for the turbulent channel flow test case.*/

  su2double val_coordsZ      = 0.0;
  if (nDim == 3) val_coordsZ = val_coords[2];

  /* Compute the primitive variables. */

  // Determine the non-dimensional y-coordinate and compute the velocity.
  su2double y, u, t;
  if (val_coords[1] < halfChan){
    y = fabs(val_coords[1]/halfChan - halfChan);
  }
  else{
    y = val_coords[1]/halfChan - halfChan;
  }
  
  u = a0*(1.0 - pow(y,alpha));
  t = TWall + (TMiddle-TWall)*(1.0 - pow(y,(alpha+2.0)));
  
  // Determine the seed of the random generator based on the coordinates.
  // This means that the random generator will always produce the same
  // number for the same coordinates.
  const int seed = (int) (36*val_coords[0] + 821*val_coords[1] + 18955*val_coords[2]);
  std::default_random_engine rand_gen(seed);
  std::uniform_real_distribution<su2double> u01(-0.1,0.1);
  
  // Compute the dimensional temperature and the density.
  su2double rho = pRef/(RGas*t);

  /* Compute the conservative variables. Note that both 2D and 3D
   cases are treated correctly. */

  val_solution[0]      = rho;
  val_solution[1]      = rho*(u + u*u01(rand_gen));
  val_solution[2]      = rho*(u*u01(rand_gen));
  val_solution[3]      = rho*(u*u01(rand_gen));
  val_solution[nVar-1] = pRef*ovGm1 + 0.5*rho*(u*u);
}

bool CTurbChannelSolution::ExactSolutionKnown(void) {return false;}

