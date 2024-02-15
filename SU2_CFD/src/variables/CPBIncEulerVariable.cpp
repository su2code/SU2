/*!
 * \file variable_direct_mean_inc.cpp
 * \brief Definition of the variable classes for incompressible flow.
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include "../../include/variables/CPBIncEulerVariable.hpp"
#include <limits>

CPBIncEulerVariable::CPBIncEulerVariable(su2double density, su2double pressure, const su2double *velocity, unsigned long nPoint,
                                     unsigned long ndim, unsigned long nvar, CConfig *config) : CFlowVariable(nPoint, ndim, nvar, nDim+4, nDim+2, config) {
                                    //  Gradient_Reconstruction(config->GetReconstructionGradientRequired() ? Gradient_Aux : Gradient_Primitive) {

  bool dual_time    = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                      (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  bool viscous      = config->GetViscous();
  bool axisymmetric = config->GetAxisymmetric();

  /*--- Allocate and initialize the primitive variables and gradients ---*/

  nPrimVar = nDim+4; nPrimVarGrad = nDim+2;

  /*--- Allocate residual structures ---*/

  Res_TruncError.resize(nPoint,nVar) = su2double(0.0);

  /*--- Only for residual smoothing (multigrid) ---*/

  for (unsigned long iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    if (config->GetMG_CorrecSmooth(iMesh) > 0) {
      Residual_Sum.resize(nPoint,nVar);
      Residual_Old.resize(nPoint,nVar);
      break;
    }
  }

  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
    Undivided_Laplacian.resize(nPoint,nVar);

  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/

  Limiter_Primitive.resize(nPoint,nPrimVarGrad) = su2double(0.0);

  Limiter.resize(nPoint,nVar) = su2double(0.0);

  Solution_Max.resize(nPoint,nPrimVarGrad) = su2double(0.0);
  Solution_Min.resize(nPoint,nPrimVarGrad) = su2double(0.0);

  /*--- Solution initialization ---*/
  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = density*velocity[iVar];

  Solution_Old = Solution;

  /*--- Allocate and initialize solution for dual time strategy ---*/

  if (dual_time) {
    Solution_time_n = Solution;
    Solution_time_n1 = Solution;
  }

  /*--- Incompressible flow, primitive variables nDim+4, (P, vx, vy, vz, rho, lamMu, EddyMu) ---*/

  Primitive.resize(nPoint,nPrimVar) = su2double(0.0);

  /*--- Incompressible flow, gradients primitive variables nDim+2, (P, vx, vy, vz, rho),
        Might need P and rho for running the adjoint problem in future.---*/

  Gradient_Primitive.resize(nPoint,nPrimVarGrad,nDim,0.0);

  if (config->GetReconstructionGradientRequired()) {
    Gradient_Aux.resize(nPoint,nPrimVarGrad,nDim,0.0);
  }
  
  if (config->GetLeastSquaresRequired()) {
    Rmatrix.resize(nPoint,nDim,nDim,0.0);
  }
  
  /*--- If axisymmetric and viscous, we need an auxiliary gradient. ---*/
  
  if (axisymmetric && viscous) Grad_AuxVar.resize(nPoint,nAuxVar,nDim, 0.0);

  Velocity2.resize(nPoint) = su2double(0.0);
  Max_Lambda_Inv.resize(nPoint) = su2double(0.0);
  Delta_Time.resize(nPoint) = su2double(0.0);

  /* Under-relaxation parameter. */
  UnderRelaxation.resize(nPoint) = su2double(1.0);
  LocalCFL.resize(nPoint) = su2double(0.0);
  
  /* Non-physical point (first-order) initialization. */
  Non_Physical.resize(nPoint) = false;
  Non_Physical_Counter.resize(nPoint) = 0;
  
 /*--- Store coefficients of momentum equation ---*/
  Mom_Coeff.resize(nPoint, nDim) = su2double(0.0);
  Mom_Coeff_nb.resize(nPoint, nDim) = su2double(0.0);
  
 /*--- Strong BC flag ---*/
 strong_bc.resize(nPoint) = false;
 
 MassFlux.resize(nPoint) = su2double(0.0);
}

bool CPBIncEulerVariable::SetPrimVar(unsigned long iPoint, su2double Density_Inf,  CConfig *config) {

  unsigned long iVar;
  bool check_dens = false, physical = true;

  /*--- Set the value of the density ---*/
  
  check_dens = SetDensity(iPoint, Density_Inf);

  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/

  SetVelocity(iPoint);

  return physical;

}

CPoissonVariable::CPoissonVariable(su2double val_SourceTerm, unsigned long nPoint,
                                   unsigned short val_nDim, unsigned short val_nVar,
                                       CConfig *config) : CVariable(nPoint, val_nDim,
                                                                    val_nVar,
                                                                    config) {
  unsigned short iVar,iMesh ,nMGSmooth = 0;

  /*--- Allocate residual structures ---*/

  Res_TruncError.resize(nPoint,nVar) = su2double(0.0);
    
  /*--- Only for residual smoothing (multigrid) ---*/

  for (unsigned long iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    if (config->GetMG_CorrecSmooth(iMesh) > 0) {
      Residual_Sum.resize(nPoint,nVar);
      Residual_Old.resize(nPoint,nVar);
      break;
    }
  }
  
 /*--- Gradient related fields ---*/

  Gradient.resize(nPoint,val_nVar,val_nDim,0.0);
  
  if (config->GetLeastSquaresRequired()) {
    Rmatrix.resize(nPoint,nDim,nDim,0.0);
  }
 
 /*--- Intitialize the source term of Poisson equation. ---*/
 
  SourceTerm.resize(nPoint) = su2double(0.0);
  Delta_Time.resize(nPoint) = su2double(0.0);
  
  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = su2double(0.0);

  Solution_Old = Solution;
  
  /*--- Strong BC flag ---*/
 strong_bc.resize(nPoint) = false;
}
