/*!
 * \file CAdjEulerVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, T. Economon
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


#include "../../include/variables/CAdjEulerVariable.hpp"

CAdjEulerVariable::CAdjEulerVariable(su2double psirho, const su2double *phi, su2double psie, unsigned long npoint, unsigned long ndim,
                                     unsigned long nvar, CConfig *config) : CVariable(npoint, ndim, nvar, config),
                                     Gradient_Reconstruction(config->GetReconstructionGradientRequired() ? Gradient_Aux : Gradient) {

  bool dual_time = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                   (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);

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

  Gradient.resize(nPoint,nVar,nDim,0.0);

  if (config->GetReconstructionGradientRequired()) {
    Gradient_Aux.resize(nPoint,nVar,nDim,0.0);
  }

  if (config->GetLeastSquaresRequired()) {
    Rmatrix.resize(nPoint,nDim,nDim,0.0);
  }

  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
  if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED)
    Undivided_Laplacian.resize(nPoint,nVar);

  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/
  Limiter.resize(nPoint,nVar) = su2double(0.0);
  Solution_Max.resize(nPoint,nVar) = su2double(0.0);
  Solution_Min.resize(nPoint,nVar) = su2double(0.0);

  /*--- Solution initialization ---*/
  su2double val_solution[5] = {psirho, phi[0], phi[1], psie, psie};
  if(nDim==3) val_solution[3] = phi[2];

  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = val_solution[iVar];

  Solution_Old = Solution;

  /*--- Allocate and initializate solution for dual time strategy ---*/
  if (dual_time) {
    Solution_time_n = Solution;
    Solution_time_n1 = Solution;
  }

  /*--- Allocate auxiliar vector for sensitivity computation ---*/
  nAuxVar = 1;
  AuxVar.resize(nPoint, nAuxVar) = su2double(0.0);
  Grad_AuxVar.resize(nPoint, nAuxVar, nDim);

  /*--- Allocate and initializate projection vector for wall boundary condition ---*/
  ForceProj_Vector.resize(nPoint,nDim) = su2double(0.0);

  /*--- Allocate space for the harmonic balance source terms ---*/
  if (config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE)
    HB_Source.resize(nPoint,nVar) = su2double(0.0);

  if (config->GetMultizone_Problem())
    Set_BGSSolution_k();

  Sensor.resize(nPoint);

  /* Non-physical point (first-order) initialization. */
  Non_Physical.resize(nPoint) = false;
  Non_Physical_Counter.resize(nPoint) = 0;

}

bool CAdjEulerVariable::SetPrimVar(unsigned long iPoint, su2double SharpEdge_Distance, bool check, CConfig *config) {

  bool RightVol = true;

  su2double adj_limit = config->GetAdjointLimit();

  bool check_dens = (fabs(Solution(iPoint,0)) > adj_limit);

  /*--- Check that the adjoint solution is bounded ---*/

  if (check_dens) {

    /*--- Copy the old solution ---*/

    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = Solution_Old(iPoint,iVar);

    RightVol = false;

  }

  return RightVol;
}
