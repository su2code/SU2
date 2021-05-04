/*!
 * \file CFEAVariable.cpp
 * \brief Definition of the variables for FEM elastic structural problems.
 * \author R. Sanchez
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

CFEAVariable::CFEAVariable(const su2double *val_fea, unsigned long npoint, unsigned long ndim,
                           unsigned long nvar, CConfig *config) :
  CVariable(npoint, ndim, config->GetTime_Domain()? 3*nvar : nvar, config),
  dynamic_analysis(config->GetTime_Domain()) {

  /*--- In time domain CVariable::nVar is mult. by 3 ^^^ (for vel. and accel.) and the original value then restored. ---*/
  nVar = nvar;

  const bool nonlinear_analysis = (config->GetGeometricConditions() == STRUCT_DEFORMATION::LARGE);
  const bool body_forces        = config->GetDeadLoad();
  const bool incremental_load   = config->GetIncrementalLoad();
  const bool prestretch_fem     = config->GetPrestretch();  // Structure is prestretched
  const bool discrete_adjoint   = config->GetDiscrete_Adjoint();
  const bool refgeom            = config->GetRefGeom(); // Reference geometry needs to be stored
  const bool multizone          = config->GetMultizone_Problem();
  const bool fsi_analysis       = config->GetFSI_Simulation() || multizone;

  VonMises_Stress.resize(nPoint) = su2double(0.0);

  if (nDim==2) Stress.resize(nPoint,3);
  else         Stress.resize(nPoint,6);

  /*--- Initialization of variables ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = val_fea[iVar];

  if (dynamic_analysis) {
    for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
      for (unsigned long iVar = 0; iVar < nVar; iVar++) {
        Solution_Vel_time_n(iPoint,iVar) = Solution_Vel(iPoint,iVar) = val_fea[iVar+nVar];
        Solution_Accel_time_n(iPoint,iVar) = Solution_Accel(iPoint,iVar) = val_fea[iVar+2*nVar];
        if (fsi_analysis) Solution_Vel_Pred(iPoint,iVar) = val_fea[iVar+nVar];
      }
    }
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

  if (multizone) Set_BGSSolution_k();

  if (config->GetTopology_Optimization()) {
    nAuxVar = 1;
    AuxVar.resize(nPoint);
  }
}
