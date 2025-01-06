/*!
 * \file CTransAFTVariable.cpp
 * \brief Definition of the solution fields.
 * \author S. Kang
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../include/variables/CTransAFTVariable.hpp"

CTransAFTVariable::CTransAFTVariable(su2double AF, su2double LnIntermittency, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
  : CTurbVariable(npoint, ndim, nvar, config) {

  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
  {
    Solution(iPoint,0) = AF;
    Solution(iPoint,1) = LnIntermittency;
  }

  Solution_Old = Solution;

  /*--- Setting CTransLMVariable of intermittency---*/
  TempVar1.resize(nPoint) = 0.0;
  TempVar2.resize(nPoint) = 0.0;
  TempVar3.resize(nPoint) = 0.0;
  TempVar4.resize(nPoint) = 0.0;
  TempVar5.resize(nPoint) = 0.0;
  TempVar6.resize(nPoint) = 0.0;
  TempVar7.resize(nPoint) = 0.0;
  TempVar8.resize(nPoint) = 0.0;
  TempVar9.resize(nPoint) = 0.0;
  TempVar10.resize(nPoint) = 0.0;
  TempVar11.resize(nPoint) = 0.0;
  TempVar12.resize(nPoint) = 0.0;
  TempVar13.resize(nPoint) = 0.0;
  TempVar14.resize(nPoint) = 0.0;
  TempVar15.resize(nPoint) = 0.0;

  nAuxVar = 2;
  Grad_AuxVar.resize(nPoint,nAuxVar,nDim,0.0);
  AuxVar.resize(nPoint,nAuxVar) = su2double(0.0);

}


void CTransAFTVariable::SetAFT_Wonder_Func(unsigned long iPoint, su2double var1, su2double var2 
          , su2double var3, su2double var4, su2double var5, su2double var6, su2double var7
          , su2double var8, su2double var9, su2double var10, su2double var11, su2double var12
          , su2double var13, su2double var14, su2double var15) {

  TempVar1(iPoint) = var1;
  TempVar2(iPoint) = var2;
  TempVar3(iPoint) = var3;
  TempVar4(iPoint) = var4;
  TempVar5(iPoint) = var5;
  TempVar6(iPoint) = var6;
  TempVar7(iPoint) = var7;
  TempVar8(iPoint) = var8;
  TempVar9(iPoint) = var9;
  TempVar10(iPoint) = var10;
  TempVar11(iPoint) = var11;
  TempVar12(iPoint) = var12;
  TempVar13(iPoint) = var13;
  TempVar14(iPoint) = var14;
  TempVar15(iPoint) = var15;
  
}
