/*!
 * \file CFEABoundVariable.cpp
 * \brief Definition of the variables for FEM elastic structural problems.
 * \author R. Sanchez
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

#include "../../include/variables/CFEABoundVariable.hpp"


CFEABoundVariable::CFEABoundVariable(void) : CFEAVariable() {

  FlowTraction          = NULL;    // Nodal traction due to the fluid (fsi)
  Residual_Ext_Surf     = NULL;    // Residual component due to external surface forces

  FlowTraction_n        = NULL;    // Nodal traction due to the fluid (fsi) at time n (for gen-alpha methods)
  Residual_Ext_Surf_n   = NULL;    // Residual component due to external surface forces at time n (for gen-alpha methods)

}

CFEABoundVariable::CFEABoundVariable(su2double *val_fea, unsigned short val_nDim, unsigned short val_nvar,
                                     CConfig *config) : CFEAVariable(val_fea, val_nDim, val_nvar, config) {

  unsigned short iVar;
  bool gen_alpha = (config->GetKind_TimeIntScheme_FEA() == GENERALIZED_ALPHA);
  bool fsi_analysis = config->GetFSI_Simulation();

  FlowTraction          = NULL;
  Residual_Ext_Surf     = NULL;
  FlowTraction_n        = NULL;
  Residual_Ext_Surf_n   = NULL;

  /*--- Surface residual ---*/
  Residual_Ext_Surf = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) Residual_Ext_Surf[iVar] = 0.0;

  /*--- Flow traction ---*/
  if (fsi_analysis){
    FlowTraction = new su2double [nVar];
    for (iVar = 0; iVar < nVar; iVar++) FlowTraction[iVar] = 0.0;
  }

  /*--- Generalized alpha integration method requires storing the old residuals ---*/
  if (gen_alpha) {
    Residual_Ext_Surf_n = new su2double [nVar];
    for (iVar = 0; iVar < nVar; iVar++) Residual_Ext_Surf_n[iVar] = 0.0;

    if (fsi_analysis) {
      FlowTraction_n = new su2double [nVar];
      for (iVar = 0; iVar < nVar; iVar++) FlowTraction_n[iVar] = 0.0;
    }
  }

}

CFEABoundVariable::~CFEABoundVariable(void) {

  if (FlowTraction          != NULL) delete [] FlowTraction;
  if (Residual_Ext_Surf     != NULL) delete [] Residual_Ext_Surf;

  if (FlowTraction_n         != NULL) delete [] FlowTraction_n;
  if (Residual_Ext_Surf_n    != NULL) delete [] Residual_Ext_Surf_n;

}
