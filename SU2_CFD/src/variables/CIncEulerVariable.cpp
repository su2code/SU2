/*!
 * \file CIncEulerVariable.cpp
 * \brief Definition of the variable classes for incompressible flow.
 * \author F. Palacios, T. Economon
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

#include "../../include/variables/CIncEulerVariable.hpp"


CIncEulerVariable::CIncEulerVariable(su2double pressure, const su2double *velocity, su2double temperature, unsigned long npoint,
                                     unsigned long ndim, unsigned long nvar, CConfig *config) : CVariable(npoint, ndim, nvar, config) {

  bool dual_time    = (config->GetTime_Marching() == DT_STEPPING_1ST) ||
                      (config->GetTime_Marching() == DT_STEPPING_2ND);
  bool viscous      = config->GetViscous();
  bool axisymmetric = config->GetAxisymmetric();

  /*--- Allocate and initialize the primitive variables and gradients ---*/

  nPrimVar = nDim+9; nPrimVarGrad = nDim+4;

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

  su2double val_solution[5] = {pressure, velocity[0], velocity[1], temperature, temperature};
  if(nDim==3) val_solution[3] = velocity[2];

  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = val_solution[iVar];

  Solution_Old = Solution;

  /*--- Allocate and initialize solution for dual time strategy ---*/

  if (dual_time) {
    Solution_time_n = Solution;
    Solution_time_n1 = Solution;
  }

  /*--- Incompressible flow, primitive variables nDim+9, (P, vx, vy, vz, T, rho, beta, lamMu, EddyMu, Kt_eff, Cp, Cv) ---*/

  Primitive.resize(nPoint,nPrimVar) = su2double(0.0);

  /*--- Incompressible flow, gradients primitive variables nDim+4, (P, vx, vy, vz, T, rho, beta),
        We need P, and rho for running the adjoint problem ---*/

  Gradient_Primitive.resize(nPoint,nPrimVarGrad,nDim,0.0);

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    Rmatrix.resize(nPoint,nDim,nDim,0.0);
  }

  /*--- If axisymmetric and viscous, we need an auxiliary gradient. ---*/

  if (axisymmetric && viscous) Grad_AuxVar.resize(nPoint,nDim);

  if (config->GetMultizone_Problem())
    Set_BGSSolution_k();

  Density_Old.resize(nPoint) = su2double(0.0);
  Velocity2.resize(nPoint) = su2double(0.0);
  Max_Lambda_Inv.resize(nPoint) = su2double(0.0);
  Delta_Time.resize(nPoint) = su2double(0.0);
  Lambda.resize(nPoint) = su2double(0.0);
  Sensor.resize(nPoint) = su2double(0.0);

}

void CIncEulerVariable::SetGradient_PrimitiveZero() {
  Gradient_Primitive.storage.setConstant(0.0);
}

bool CIncEulerVariable::SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) {

  unsigned long iVar;
  bool check_dens = false, check_temp = false, physical = true;

  /*--- Store the density from the previous iteration. ---*/

  Density_Old(iPoint) = GetDensity(iPoint);

  /*--- Set the value of the pressure ---*/

  SetPressure(iPoint);

  /*--- Set the value of the temperature directly ---*/

  su2double Temperature = Solution(iPoint,nDim+1);
  check_temp = SetTemperature(iPoint,Temperature);

  /*--- Use the fluid model to compute the new value of density.
  Note that the thermodynamic pressure is constant and decoupled
  from the dynamic pressure being iterated. ---*/

  /*--- Use the fluid model to compute the new value of density. ---*/

  FluidModel->SetTDState_T(Temperature);

  /*--- Set the value of the density ---*/

  check_dens = SetDensity(iPoint, FluidModel->GetDensity());

  /*--- Non-physical solution found. Revert to old values. ---*/

  if (check_dens || check_temp) {

    /*--- Copy the old solution ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint, iVar) = Solution_Old(iPoint, iVar);

    /*--- Recompute the primitive variables ---*/

    Temperature = Solution(iPoint, nDim+1);
    SetTemperature(iPoint, Temperature);
    FluidModel->SetTDState_T(Temperature);
    SetDensity(iPoint, FluidModel->GetDensity());

    /*--- Flag this point as non-physical. ---*/

    physical = false;

  }

  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/

  SetVelocity(iPoint);

  /*--- Set specific heats (only necessary for consistency with preconditioning). ---*/

  SetSpecificHeatCp(iPoint, FluidModel->GetCp());
  SetSpecificHeatCv(iPoint, FluidModel->GetCv());

  return physical;

}
