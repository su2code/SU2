/*!
 * \file CPassiveScalarVariable.cpp
 * \brief Definition of the variable fields for the passive scalar class.
 * \author T. Economon
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

#include "../../include/variables/CPassiveScalarVariable.hpp"

CPassiveScalarVariable::CPassiveScalarVariable(su2double     *val_scalar_inf,
                                               unsigned long npoint,
                                               unsigned long ndim,
                                               unsigned long nvar,
                                               CConfig       *config)
: CScalarVariable(npoint, ndim, nvar, config) {
  
  for (unsigned long iPoint=0; iPoint<nPoint; ++iPoint) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      Solution(iPoint,iVar) = val_scalar_inf[iVar];
    }
  }
  
  Solution_Old = Solution;
  
  /*--- Allocate and initialize solution for the dual time strategy ---*/
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  
  if (dual_time) {
    Solution_time_n  = Solution;
    Solution_time_n1 = Solution;
  }
  
  /*--- Allocate space for the mass diffusivity. ---*/
  
  Diffusivity.resize(nPoint,nVar) = su2double(0.0);
  
  /*--- If axisymmetric and viscous, we need an auxiliary gradient. ---*/
  
  if (config->GetAxisymmetric() && config->GetViscous()) {
    AuxVar.resize(nPoint);
    Grad_AuxVar.resize(nPoint,nDim);
  }
  
}
