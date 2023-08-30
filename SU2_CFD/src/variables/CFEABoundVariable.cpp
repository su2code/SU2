/*!
 * \file CFEABoundVariable.cpp
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


#include "../../include/variables/CFEABoundVariable.hpp"


CFEABoundVariable::CFEABoundVariable(const su2double *val_fea, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
  : CFEAVariable(val_fea, npoint, ndim, nvar, config) {

  VertexMap.Reset(nPoint);
}

void CFEABoundVariable::AllocateBoundaryVariables(CConfig *config) {

  if (VertexMap.GetIsValid()) return; // nothing to do

  /*--- Count number of vertices and build map ---*/

  unsigned long nBoundPt = VertexMap.Build();

  /*--- Allocate ---*/

  fsi_analysis = config->GetFSI_Simulation();

  /*--- Surface residual ---*/
  Residual_Ext_Surf.resize(nBoundPt,nVar) = su2double(0.0);

  /*--- Flow traction ---*/
  if (fsi_analysis) FlowTraction.resize(nBoundPt,nVar) = su2double(0.0);

  /*--- Generalized alpha integration method requires storing the old residuals ---*/
  if (config->GetKind_TimeIntScheme_FEA() == STRUCT_TIME_INT::GENERALIZED_ALPHA) {
    Residual_Ext_Surf_n.resize(nBoundPt,nVar) = su2double(0.0);

    if (fsi_analysis) FlowTraction_n.resize(nBoundPt,nVar) = su2double(0.0);
  }
}

void CFEABoundVariable::Set_FlowTraction_n() { FlowTraction_n = FlowTraction; }

void CFEABoundVariable::Set_SurfaceLoad_Res_n() { Residual_Ext_Surf_n = Residual_Ext_Surf; }

void CFEABoundVariable::Clear_SurfaceLoad_Res() { Residual_Ext_Surf.setConstant(0.0); }

void CFEABoundVariable::RegisterFlowTraction(bool reset) {
  if (!fsi_analysis) return;
  for (unsigned long iVertex = 0; iVertex < FlowTraction.rows(); iVertex++)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      if (reset) AD::ResetInput(FlowTraction(iVertex,iVar));
      else AD::RegisterInput(FlowTraction(iVertex,iVar));
}
