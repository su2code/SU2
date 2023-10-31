/*!
 * \file CFEAVariable.cpp
 * \brief Definition of the variables for FEM elastic structural problems.
 * \author R. Sanchez
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

#include "../../include/variables/CFEAVariable.hpp"

CFEAVariable::CFEAVariable(const su2double *val_fea, unsigned long npoint, unsigned long ndim,
                           unsigned long nvar, CConfig *config) :
  CVariable(npoint, ndim, config->GetTime_Domain()? 3*nvar : nvar, config) {

  /*--- In time domain CVariable::nVar is mult. by 3 ^^^ (for vel. and accel.)
   * and the original value then restored (below). ---*/
  nVar = nvar;
  /*--- This simplifies the discrete adjoint of this solver, as it allows an abstract
   * treatment of the "state" (disp. vel. accel.) whose details are only important
   * for the primal solver. In time domain the primal "believes" it has nVar variables,
   * which it uses for linear solvers, and then handles velocity and acceleration
   * explicitly (for time integration). Whereas the discrete adjoint "thinks" the
   * primal solution has 3*nVar variables. This is a little different from simply
   * giving names to parts of the solution, it requires the methods of CVariable that
   * deal with adjoints to deduce "nVar" from the container, rather than relying on
   * the nVar member (which is manipulated above, so that CVariable::SetSolution, etc.
   * still work as expected for the primal solver). ---*/

  const bool dynamic_analysis   = config->GetTime_Domain();
  const bool body_forces        = config->GetDeadLoad();
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

  if (fsi_analysis) {
    Solution_Pred.resize(nPoint,nVar);
    Solution_Pred_Old.resize(nPoint,nVar);

    for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
      for (unsigned long iVar = 0; iVar < nVar; iVar++)
        Solution_Pred(iPoint,iVar) = Solution_Pred_Old(iPoint,iVar) = val_fea[iVar];
  }

  if (dynamic_analysis) {
    if (fsi_analysis) Solution_Vel_Pred.resize(nPoint,nVar);

    for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
      for (unsigned long iVar = 0; iVar < nVar; iVar++) {
        Solution_Vel_time_n(iPoint,iVar) = Solution_Vel(iPoint,iVar) = val_fea[iVar+nVar];
        Solution_Accel_time_n(iPoint,iVar) = Solution_Accel(iPoint,iVar) = val_fea[iVar+2*nVar];
        if (fsi_analysis) Solution_Vel_Pred(iPoint,iVar) = val_fea[iVar+nVar];
      }
    }
  }

  if (discrete_adjoint && fsi_analysis) Solution_Old = Solution;

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
