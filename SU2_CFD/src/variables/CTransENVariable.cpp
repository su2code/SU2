/*!
 * \file CTransENVariable.cpp
 * \brief Definition of the solution fields.
 * \author R. Roos
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/variables/CTransENVariable.hpp"

CTransENVariable::CTransENVariable(su2double AmplificationFactor, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
  : CTurbVariable(npoint, ndim, nvar, config) {

  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
  {
    Solution(iPoint,0) = AmplificationFactor;
  }

  Solution_Old = Solution;
  
}

void CTransENVariable::SetAmplificationFactor(unsigned long iPoint, su2double val_AmplificationFactor) {
  AD::StartPreacc();
  AmplificationFactor(iPoint) = val_AmplificationFactor;
  AD::SetPreaccOut(AmplificationFactor(iPoint)); //Still need to be fixed?
  AD::EndPreacc();
}
