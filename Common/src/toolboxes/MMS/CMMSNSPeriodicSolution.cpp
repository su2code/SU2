/*!
 * \file CMMSNSPeriodicSolution.cpp
 * \brief Implementations of the member functions of CMMSNSPeriodicSolution.
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

#include "../../../include/toolboxes/MMS/CMMSNSPeriodicSolution.hpp"

CMMSNSPeriodicSolution::CMMSNSPeriodicSolution(void) : CVerificationSolution() { }

CMMSNSPeriodicSolution::CMMSNSPeriodicSolution(unsigned short val_nDim,
                                               unsigned short val_nVar,
                                               unsigned short val_iMesh,
                                               CConfig*       config)
: CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
  
  /*--- Write a message that the solution is initialized for the manufactured
   solution for the Navier-Stokes equations on a periodic quad. ---*/
  
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the manufactured solution " << endl;
    cout << "         of the Navier-Stokes equations on a periodic quad!!!" << endl;
    cout << endl << flush;
  }
  
  /*--- Coefficients, needed to determine the solution. ---*/
  
  const su2double Prandtl = config->GetPrandtl_Lam();
  
  RGas         = config->GetGas_ConstantND();
  Gamma        = config->GetGamma();
  Viscosity    = config->GetMu_ConstantND();
  Conductivity = Viscosity*Gamma*RGas/(Prandtl*(Gamma-1.0));
  
  /*--- Constants, which describe this manufactured solution. This is a viscous
   solution on a [-1,1]x[-1,1] quad, where the conservative variables vary as a
   combination of sine functions from Lopez-Morales et al., "Verification and
   Validation of HiFiLES: a High-Order LES unstructured solver on
   multi-GPU platforms," AIAA Paper 2014-3168. The solution is periodic and
   can be used to test the periodic BCs in the solver. ---*/
  
  k = PI_NUMBER;
  a = 3.0;
  
  /*--- Perform some sanity and error checks for this solution here. ---*/
  
  if(config->GetUnsteady_Simulation() != STEADY)
    SU2_MPI::Error("Steady mode must be selected for the MMS NS Periodic Quad case",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_Regime() != COMPRESSIBLE)
    SU2_MPI::Error("Compressible flow equations must be selected for the MMS NS Periodic Quad case",
                   CURRENT_FUNCTION);
  
  if((config->GetKind_Solver() != NAVIER_STOKES) &&
     (config->GetKind_Solver() != FEM_NAVIER_STOKES))
    SU2_MPI::Error("Navier Stokes equations must be selected for the MMS NS Periodic Quad case",
                   CURRENT_FUNCTION);
  
  if((config->GetKind_FluidModel() != STANDARD_AIR) &&
     (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the MMS NS Periodic Quad case",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_ViscosityModel() != CONSTANT_VISCOSITY)
    SU2_MPI::Error("Constant viscosity must be selected for the MMS NS Periodic Quad case",
                   CURRENT_FUNCTION);
  
  if(config->GetKind_ConductivityModel() != CONSTANT_PRANDTL)
    SU2_MPI::Error("Constant Prandtl number must be selected for the MMS NS Periodic Quad case",
                   CURRENT_FUNCTION);
}

CMMSNSPeriodicSolution::~CMMSNSPeriodicSolution(void) { }

void CMMSNSPeriodicSolution::GetBCState(const su2double *val_coords,
                                        const su2double val_t,
                                        su2double       *val_solution) {
  
  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CMMSNSPeriodicSolution::GetSolution(const su2double *val_coords,
                                         const su2double val_t,
                                         su2double       *val_solution) {
  
  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  
  /*--- Set the conservative variables. This is a viscous manufactured
   solution on a [-1,1]x[-1,1] quad, where the conservative variables vary as a
   combination of sine functions from Lopez-Morales et al., "Verification and
   Validation of HiFiLES: a High-Order LES unstructured solver on
   multi-GPU platforms," AIAA Paper 2014-3168. The solution is periodic and
   can be used to test the periodic BCs in the solver. ---*/
  val_solution[0]      = sin(k*(x+y)) + a;
  val_solution[1]      = sin(k*(x+y)) + a;
  val_solution[2]      = sin(k*(x+y)) + a;
  val_solution[3]      = 0.0;
  val_solution[nDim+1] = pow(sin(k*(x+y)) + a,2.0);
  
}

void CMMSNSPeriodicSolution::GetPrimitiveGradient(const su2double *val_coords,
                                                  const su2double val_t,
                                                  su2double       **val_gradient) {
  
  /* Analytic values for the primitive variable gradients.
   The source code for the gradient values is generated in Python.
   See the file CMMSNSPeriodicSolution.py in the directory
   CreateMMSSourceTerms for the details on how to do this.*/
  
  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  
  /* Gradient of temperature in the x- and y-directions. */
  val_gradient[0][0] = -k*(Gamma - 1.0)*(-a + pow(a + sin(k*(x + y)), 2.0) - sin(k*(x + y)))*cos(k*(x + y))/(RGas*pow(a + sin(k*(x + y)), 2)) + (Gamma - 1.0)*(2.0*k*pow(a + sin(k*(x + y)), 1.0)*cos(k*(x + y)) - k*cos(k*(x + y)))/(RGas*(a + sin(k*(x + y))));
  val_gradient[0][1] = -k*(Gamma - 1.0)*(-a + pow(a + sin(k*(x + y)), 2.0) - sin(k*(x + y)))*cos(k*(x + y))/(RGas*pow(a + sin(k*(x + y)), 2)) + (Gamma - 1.0)*(2.0*k*pow(a + sin(k*(x + y)), 1.0)*cos(k*(x + y)) - k*cos(k*(x + y)))/(RGas*(a + sin(k*(x + y))));
  if (nDim == 3) val_gradient[0][2] = 0.0;
  
  /* Gradient of x-velocity in the x- and y-directions. */
  val_gradient[1][0] = 0.0;
  val_gradient[1][1] = 0.0;
  if (nDim == 3) val_gradient[1][2] = 0.0;
  
  /* Gradient of y-velocity in the x- and y-directions. */
  val_gradient[2][0] = 0.0;
  val_gradient[2][1] = 0.0;
  if (nDim == 3) val_gradient[2][2] = 0.0;
  
  /* Gradient of z-velocity in the x- and y-directions. */
  val_gradient[3][0] = 0.0;
  val_gradient[3][1] = 0.0;
  if (nDim == 3) val_gradient[3][2] = 0.0;
  
  /* Gradient of pressure in the x- and y-directions. */
  val_gradient[nDim+1][0] = (Gamma - 1.0)*(2.0*k*pow(a + sin(k*(x + y)), 1.0)*cos(k*(x + y)) - k*cos(k*(x + y)));
  val_gradient[nDim+1][1] = (Gamma - 1.0)*(2.0*k*pow(a + sin(k*(x + y)), 1.0)*cos(k*(x + y)) - k*cos(k*(x + y)));
  if (nDim == 3) val_gradient[nDim+1][2] = 0.0;
  
  
}

void CMMSNSPeriodicSolution::GetMMSSourceTerm(const su2double *val_coords,
                                              const su2double val_t,
                                              su2double       *val_source) {
  
  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];
  
  /* The source code for the source terms is generated in Python.
   See the file CMMSNSPeriodicSolution.py in the directory
   CreateMMSSourceTerms for the details how to do this. */
  val_source[0]      = 2*k*cos(k*(x + y));
  val_source[1]      = k*((Gamma - 1.0)*(2.0*pow(a + sin(k*(x + y)), 1.0) - 1) + 2)*cos(k*(x + y));
  val_source[2]      = k*((Gamma - 1.0)*(2.0*pow(a + sin(k*(x + y)), 1.0) - 1) + 2)*cos(k*(x + y));
  val_source[3]      = 0.0;
  val_source[nDim+1] = k*(2*Conductivity*k*(Gamma - 1.0)*(-pow(a + sin(k*(x + y)), 2)*(-2.0*a*sin(k*(x + y)) + 1.0*sin(k*(x + y)) + 2.0*cos(2*k*(x + y))) + (a + sin(k*(x + y)))*(2*(2.0*pow(a + sin(k*(x + y)), 1.0) - 1)*pow(cos(k*(x + y)), 2) + (a - pow(a + sin(k*(x + y)), 2.0) + sin(k*(x + y)))*sin(k*(x + y))) + 2*(a - pow(a + sin(k*(x + y)), 2.0) + sin(k*(x + y)))*pow(cos(k*(x + y)), 2)) + RGas*pow(a + sin(k*(x + y)), 3)*(4.0*Gamma*pow(a + sin(k*(x + y)), 1.0) - 2.0*Gamma + 2.0)*cos(k*(x + y)))/(RGas*pow(a + sin(k*(x + y)), 3));
  
}

bool CMMSNSPeriodicSolution::IsManufacturedSolution(void) {
  return true;
}

bool CMMSNSPeriodicSolution::ExactPrimitiveGradientKnown(void) {
  return true;
}

