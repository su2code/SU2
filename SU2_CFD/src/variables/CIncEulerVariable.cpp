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

CIncEulerVariable::CIncEulerVariable(void) : CVariable() {

  /*--- Array initialization ---*/

  Primitive          = NULL;
  Gradient_Primitive = NULL;
  Limiter_Primitive  = NULL;

  Grad_AuxVar = NULL;

  nPrimVar     = 0;
  nPrimVarGrad = 0;

  nSecondaryVar     = 0;
  nSecondaryVarGrad = 0;

  Solution_BGS_k = NULL;

}

CIncEulerVariable::CIncEulerVariable(su2double val_pressure, su2double *val_velocity, su2double val_temperature,
                                     unsigned short val_nDim, unsigned short val_nvar, CConfig *config) :
                                     CVariable(val_nDim, val_nvar, config) {

  unsigned short iVar, iDim, iMesh, nMGSmooth = 0;

  bool dual_time    = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                       (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool viscous      = config->GetViscous();
  bool axisymmetric = config->GetAxisymmetric();
  bool fsi          = config->GetFSI_Simulation();
  bool multizone = config->GetMultizone_Problem();

  /*--- Array initialization ---*/

  Primitive          = NULL;
  Gradient_Primitive = NULL;
  Limiter_Primitive  = NULL;

  Grad_AuxVar = NULL;

  nPrimVar     = 0;
  nPrimVarGrad = 0;

  nSecondaryVar     = 0;
  nSecondaryVarGrad = 0;

  /*--- Allocate and initialize the primitive variables and gradients ---*/

  nPrimVar = nDim+9; nPrimVarGrad = nDim+4;

  /*--- Allocate residual structures ---*/

  Res_TruncError = new su2double [nVar];

  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }

  /*--- Only for residual smoothing (multigrid) ---*/

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);

  if (nMGSmooth > 0) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }

  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    Undivided_Laplacian = new su2double [nVar];
  }

  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/

  Limiter_Primitive = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Limiter_Primitive[iVar] = 0.0;

  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;

  Solution_Max = new su2double [nPrimVarGrad];
  Solution_Min = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }

  /*--- Solution and old solution initialization ---*/

  Solution[0] = val_pressure;
  Solution_Old[0] = val_pressure;
  for (iDim = 0; iDim < nDim; iDim++) {
    Solution[iDim+1] = val_velocity[iDim];
    Solution_Old[iDim+1] = val_velocity[iDim];
  }
  Solution[nDim+1] = val_temperature;
  Solution_Old[nDim+1] = val_temperature;

  /*--- Allocate and initialize solution for dual time strategy ---*/

  if (dual_time) {
    Solution_time_n[0]  =  val_pressure;
    Solution_time_n1[0] =  val_pressure;
    for (iDim = 0; iDim < nDim; iDim++) {
      Solution_time_n[iDim+1] = val_velocity[iDim];
      Solution_time_n1[iDim+1] = val_velocity[iDim];
    }
    Solution[nDim+1] = val_temperature;
    Solution_Old[nDim+1] = val_temperature;
  }

  /*--- Incompressible flow, primitive variables nDim+9, (P, vx, vy, vz, T, rho, beta, lamMu, EddyMu, Kt_eff, Cp, Cv) ---*/

  Primitive = new su2double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;

  /*--- Incompressible flow, gradients primitive variables nDim+4, (P, vx, vy, vz, T, rho, beta)
   * We need P, and rho for running the adjoint problem ---*/

  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }

  /*--- If axisymmetric and viscous, we need an auxiliary gradient. ---*/

  if (axisymmetric && viscous)
    Grad_AuxVar = new su2double[nDim];

  Solution_BGS_k = NULL;
  if (fsi || multizone){
    Solution_BGS_k  = new su2double [nVar];
    Solution_BGS_k[0] = val_pressure;
    for (iDim = 0; iDim < nDim; iDim++) {
      Solution_BGS_k[iDim+1] = val_velocity[iDim]*config->GetDensity_FreeStreamND();
    }
    Solution_BGS_k[nDim+1] = val_temperature;
  }

}

CIncEulerVariable::CIncEulerVariable(su2double *val_solution, unsigned short val_nDim, unsigned short val_nvar,
                                     CConfig *config) : CVariable(val_nDim, val_nvar, config) {

  unsigned short iVar, iDim, iMesh, nMGSmooth = 0;

  bool dual_time    = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool viscous      = config->GetViscous();
  bool axisymmetric = config->GetAxisymmetric();
  bool fsi = config->GetFSI_Simulation();
  bool multizone = config->GetMultizone_Problem();

  /*--- Array initialization ---*/

  Primitive          = NULL;
  Gradient_Primitive = NULL;
  Limiter_Primitive  = NULL;

  Grad_AuxVar = NULL;

  nPrimVar     = 0;
  nPrimVarGrad = 0;

  nSecondaryVar     = 0;
  nSecondaryVarGrad = 0;

  /*--- Allocate and initialize the primitive variables and gradients ---*/

  nPrimVar = nDim+9; nPrimVarGrad = nDim+4;

  /*--- Allocate residual structures ---*/

  Res_TruncError = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }

  /*--- Only for residual smoothing (multigrid) ---*/

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);

  if (nMGSmooth > 0) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }

  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
    Undivided_Laplacian = new su2double [nVar];

  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/

  Limiter_Primitive = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Limiter_Primitive[iVar] = 0.0;

  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;

  Solution_Max = new su2double [nPrimVarGrad];
  Solution_Min = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }

  /*--- Solution initialization ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Solution[iVar] = val_solution[iVar];
    Solution_Old[iVar] = val_solution[iVar];
  }

  /*--- Allocate and initialize solution for dual time strategy ---*/

  if (dual_time) {
    Solution_time_n = new su2double [nVar];
    Solution_time_n1 = new su2double [nVar];

    for (iVar = 0; iVar < nVar; iVar++) {
      Solution_time_n[iVar] = val_solution[iVar];
      Solution_time_n1[iVar] = val_solution[iVar];
    }
  }

  /*--- Incompressible flow, primitive variables nDim+9, (P, vx, vy, vz, T, rho, beta, lamMu, EddyMu, Kt_eff, Cp, Cv) ---*/

  Primitive = new su2double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;

  /*--- Incompressible flow, gradients primitive variables nDim+4, (P, vx, vy, vz, T, rho, beta),
        We need P, and rho for running the adjoint problem ---*/

  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }

  /*--- If axisymmetric and viscous, we need an auxiliary gradient. ---*/

  if (axisymmetric && viscous)
    Grad_AuxVar = new su2double[nDim];

  Solution_BGS_k = NULL;
  if (fsi || multizone){
      Solution_BGS_k  = new su2double [nVar];
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_BGS_k[iVar] = val_solution[iVar];
      }
  }

}

CIncEulerVariable::~CIncEulerVariable(void) {
  unsigned short iVar;

  if (Primitive         != NULL) delete [] Primitive;
  if (Limiter_Primitive != NULL) delete [] Limiter_Primitive;

  if (Gradient_Primitive != NULL) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      if (Gradient_Primitive!=NULL) delete [] Gradient_Primitive[iVar];
    delete [] Gradient_Primitive;
  }

  if (Solution_BGS_k  != NULL) delete [] Solution_BGS_k;

}

void CIncEulerVariable::SetGradient_PrimitiveZero(unsigned short val_primvar) {
  unsigned short iVar, iDim;

  for (iVar = 0; iVar < val_primvar; iVar++)
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
}


su2double CIncEulerVariable::GetProjVel(su2double *val_vector) {
  su2double ProjVel;
  unsigned short iDim;

  ProjVel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ProjVel += Primitive[iDim+1]*val_vector[iDim];

  return ProjVel;
}

bool CIncEulerVariable::SetPrimVar(CFluidModel *FluidModel) {

  unsigned short iVar;
  bool check_dens = false, check_temp = false, physical = true;

  /*--- Store the density from the previous iteration. ---*/

  Density_Old = GetDensity();

  /*--- Set the value of the pressure ---*/

  SetPressure();

  /*--- Set the value of the temperature directly ---*/

  su2double Temperature = Solution[nDim+1];
  check_temp = SetTemperature(Temperature);

  /*--- Use the fluid model to compute the new value of density.
  Note that the thermodynamic pressure is constant and decoupled
  from the dynamic pressure being iterated. ---*/

  /*--- Use the fluid model to compute the new value of density. ---*/

  FluidModel->SetTDState_T(Temperature);

  /*--- Set the value of the density ---*/

  check_dens = SetDensity(FluidModel->GetDensity());

  /*--- Non-physical solution found. Revert to old values. ---*/

  if (check_dens || check_temp) {

    /*--- Copy the old solution ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];

    /*--- Recompute the primitive variables ---*/

    Temperature = Solution[nDim+1];
    SetTemperature(Temperature);
    FluidModel->SetTDState_T(Temperature);
    SetDensity(FluidModel->GetDensity());

    /*--- Flag this point as non-physical. ---*/

    physical = false;

  }

  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/

  SetVelocity();

  /*--- Set specific heats (only necessary for consistency with preconditioning). ---*/

  SetSpecificHeatCp(FluidModel->GetCp());
  SetSpecificHeatCv(FluidModel->GetCv());

  return physical;

}
