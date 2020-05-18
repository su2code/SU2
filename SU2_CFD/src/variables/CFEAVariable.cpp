/*!
 * \file CFEAVariable.cpp
 * \brief Definition of the variables for FEM elastic structural problems.
 * \author R. Sanchez
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/variables/CFEAVariable.hpp"


CFEAVariable::CFEAVariable(const su2double *val_fea, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
  : CVariable(npoint, ndim, nvar, config) {

  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);
  bool body_forces        = config->GetDeadLoad();
  bool incremental_load   = config->GetIncrementalLoad();
  bool prestretch_fem     = config->GetPrestretch();  // Structure is prestretched
  bool discrete_adjoint   = config->GetDiscrete_Adjoint();
  bool refgeom            = config->GetRefGeom(); // Reference geometry needs to be stored
  bool dynamic_analysis   = config->GetTime_Domain();
  bool fsi_analysis       = config->GetFSI_Simulation();

  VonMises_Stress.resize(nPoint) = su2double(0.0);

  if (nDim==2) Stress.resize(nPoint,3);
  else         Stress.resize(nPoint,6);

  /*--- Initialization of variables ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = val_fea[iVar];

  if (dynamic_analysis) {
    Solution_Vel.resize(nPoint,nVar);
    Solution_Accel.resize(nPoint,nVar);

    for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
      for (unsigned long iVar = 0; iVar < nVar; iVar++) {
        Solution_Vel(iPoint,iVar) = val_fea[iVar+nVar];
        Solution_Accel(iPoint,iVar) = val_fea[iVar+2*nVar];
      }
    }
    Solution_Vel_time_n = Solution_Vel;
    Solution_Accel_time_n = Solution_Accel;
  }

  if (fsi_analysis) {
    Solution_Pred = Solution;
    Solution_Pred_Old = Solution;
  }

  /*--- If we are going to use incremental analysis, we need a way to store the old solution ---*/

  if (incremental_load && nonlinear_analysis) Solution_Old.resize(nPoint,nVar) = su2double(0.0);

  /*--- If we are running a discrete adjoint iteration, we need this vector for cross-dependencies ---*/

  else if (discrete_adjoint && fsi_analysis) Solution_Old = Solution;

  /*--- Body residual ---*/
  if (body_forces) Residual_Ext_Body.resize(nPoint,nVar) = su2double(0.0);

  if (refgeom) Reference_Geometry.resize(nPoint,nVar);

  if (prestretch_fem) Prestretch.resize(nPoint,nVar);

  if (config->GetMultizone_Problem())
    Set_BGSSolution_k();
}

void CFEAVariable::SetSolution_Vel_time_n() { Solution_Vel_time_n = Solution_Vel; }

void CFEAVariable::SetSolution_Accel_time_n() { Solution_Accel_time_n = Solution_Accel; }

void CFEAVariable::Register_femSolution_time_n() {
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      AD::RegisterInput(Solution_time_n(iPoint,iVar));
}

void CFEAVariable::RegisterSolution_Vel(bool input) {
  if (input) {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
      for (unsigned long iVar = 0; iVar < nVar; iVar++)
        AD::RegisterInput(Solution_Vel(iPoint,iVar));
  }
  else {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
      for (unsigned long iVar = 0; iVar < nVar; iVar++)
        AD::RegisterOutput(Solution_Vel(iPoint,iVar));
  }
}

void CFEAVariable::RegisterSolution_Vel_time_n() {
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      AD::RegisterInput(Solution_Vel_time_n(iPoint,iVar));
}

void CFEAVariable::RegisterSolution_Accel(bool input) {
  if (input) {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
      for (unsigned long iVar = 0; iVar < nVar; iVar++)
        AD::RegisterInput(Solution_Accel(iPoint,iVar));
  }
  else {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
      for (unsigned long iVar = 0; iVar < nVar; iVar++)
        AD::RegisterOutput(Solution_Accel(iPoint,iVar));
  }
}

void CFEAVariable::RegisterSolution_Accel_time_n() {
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      AD::RegisterInput(Solution_Accel_time_n(iPoint,iVar));
}
