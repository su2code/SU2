/*!
 * \file CTransLMVariable.cpp
 * \brief Definition of the solution fields.
 * \author A. Aranake
 * \version 7.0.2 "Blackbird"
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


#include "../../include/variables/CTransLMVariable.hpp"

/* develop version
CTransLMVariable::CTransLMVariable(su2double intermittency, su2double REth, unsigned long npoint, unsigned long ndim,
                                   unsigned long nvar, CConfig *config) : CTurbVariable(npoint, ndim, nvar, config) {

  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
  {
    Solution_Old(iPoint,0) = Solution(iPoint,0) = intermittency;
    Solution_Old(iPoint,1) = Solution(iPoint,1) = REth;
  }

  if (config->GetMultizone_Problem())
    Set_BGSSolution_k();

  gamma_sep.resize(nPoint);
}*/

// LM branch modified
//CTransLMVariable::CTransLMVariable(void) : CTurbVariable() {}

CTransLMVariable::CTransLMVariable(const su2double      val_intermittency,
                                   const su2double      val_REth,
                                   const unsigned long val_nPoint,
                                   const unsigned long val_nDim,
                                   const unsigned long val_nvar,
                                   CConfig              *config)
: CTurbVariable(val_nPoint, val_nDim, val_nvar, config) {

  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  
  gamma_sep.resize(nPoint) = val_intermittency;
  /* Initialization of variables */
  for (unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
  {
    Solution(iPoint,0) = val_intermittency; Solution_Old(iPoint,0) = val_intermittency;
    Solution(iPoint,1) = val_REth;          Solution_Old(iPoint,1) = val_REth;
  }
  /* Initialize gamma_sep to the intermittency. */
  //gamma_sep.resize(nPoint);//= val_intermittency;
  
  muT.resize(nPoint) = su2double(0.0);

  /*--- Allocate and initialize solution for the dual time strategy ---*/
  if (dual_time) {
    for (unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
    {
      Solution_time_n(iPoint,0)  = val_intermittency; Solution_time_n(iPoint,0)  = val_REth;
      Solution_time_n1(iPoint,1) = val_intermittency; Solution_time_n1(iPoint,1) = val_REth;
    }
  }
}
