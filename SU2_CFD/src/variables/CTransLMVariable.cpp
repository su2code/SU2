/*!
 * \file CTransLMVariable.cpp
 * \brief Definition of the solution fields.
 * \author A. Aranake, S. Kang
 * \version 7.5.0 "Blackbird"
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


#include "../../include/variables/CTransLMVariable.hpp"

CTransLMVariable::CTransLMVariable(su2double Intermittency, su2double ReThetaT, su2double gammaSep, su2double gammaEff, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config) 
  : CTurbVariable(npoint, ndim, nvar, config) {

  LM_ParsedOptions options = config->GetLMParsedOptions();

  if (!options.SLM) {
    for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
    {
      Solution(iPoint,0) = Intermittency;
      Solution(iPoint,1) = ReThetaT;
    }
  } else {
    for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
    {
      Solution(iPoint,0) = Intermittency;
    }
  }

  Solution_Old = Solution;

  /*--- Setting CTransLMVariable of intermittency_Eff---*/
  Intermittency_Eff.resize(nPoint) = gammaEff;
  Intermittency_Sep.resize(nPoint) = gammaSep;

  if (options.SLM) {
    Re_t.resize(nPoint) = ReThetaT;

    nAuxVar = 1;
    Grad_AuxVar.resize(nPoint, nAuxVar, nDim, 0.0);
    AuxVar.resize(nPoint, nAuxVar) = su2double(0.0);
    Tu.resize(nPoint) = 0.0;

    normal_x.resize(nPoint) = 0.0;
    normal_y.resize(nPoint) = 0.0;
  }


  Re_v.resize(nPoint) = 0.0;
  Corr_Rec.resize(nPoint) = 0.0;
  Prod.resize(nPoint) = 0.0;
  Destr.resize(nPoint) = 0.0;
  F_onset1.resize(nPoint) = 0.0;
  F_onset2.resize(nPoint) = 0.0;
  F_onset3.resize(nPoint) = 0.0;
  F_onset.resize(nPoint) = 0.0;
  Lambda_theta.resize(nPoint) = 0.0;
  duds.resize(nPoint) = 0.0;
  
}


void CTransLMVariable::SetIntermittencyEff(unsigned long iPoint, su2double val_Intermittency_sep) {

  /*--- Effective intermittency ---*/
  Intermittency_Eff(iPoint) = max(Solution(iPoint,0), val_Intermittency_sep);

}

void CTransLMVariable::SetIntermittencySep(unsigned long iPoint, su2double val_Intermittency_sep) {
  Intermittency_Sep(iPoint) = val_Intermittency_sep;
}

void CTransLMVariable::SetCorr_Rec(unsigned long iPoint, su2double val_Corr_Rec) {
  Corr_Rec(iPoint) = val_Corr_Rec;
}

void CTransLMVariable::SetRe_t(unsigned long iPoint, su2double val_Re_t) {
  Re_t(iPoint) = val_Re_t;
}
void CTransLMVariable::SetTu(unsigned long iPoint, su2double val_Tu) {
  Tu(iPoint) = val_Tu;
}
void CTransLMVariable::SetLambda_theta(unsigned long iPoint, su2double val_Lambda_theta) {
  Lambda_theta(iPoint) = val_Lambda_theta;
}
void CTransLMVariable::Setduds(unsigned long iPoint, su2double val_duds) {
  duds(iPoint) = val_duds;
}
void CTransLMVariable::SetRe_v(unsigned long iPoint, su2double val_Re_v) {
  Re_v(iPoint) = val_Re_v;
}
void CTransLMVariable::SetProd(unsigned long iPoint, su2double val_Prod) {
  Prod(iPoint) = val_Prod;
}
void CTransLMVariable::SetDestr(unsigned long iPoint, su2double val_Destr) {
  Destr(iPoint) = val_Destr;
}
void CTransLMVariable::SetF_onset1(unsigned long iPoint, su2double val_F_onset1) {
  F_onset1(iPoint) = val_F_onset1;
}
void CTransLMVariable::SetF_onset2(unsigned long iPoint, su2double val_F_onset2) {
  F_onset2(iPoint) = val_F_onset2;
}
void CTransLMVariable::SetF_onset3(unsigned long iPoint, su2double val_F_onset3) {
  F_onset3(iPoint) = val_F_onset3;
}
void CTransLMVariable::SetF_onset(unsigned long iPoint, su2double val_F_onset) {
  F_onset(iPoint) = val_F_onset;
}
void CTransLMVariable::SetNormal(unsigned long iPoint, su2double val_normal_x, su2double val_normal_y) {
  normal_x(iPoint) = val_normal_x;
  normal_y(iPoint) = val_normal_y;
}