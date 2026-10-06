/*!
 * \file CTransLMVariable.cpp
 * \brief Definition of the solution fields.
 * \author A. Aranake, S. Kang
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

  TransitionData.resize(nPoint);
  if (options.SLM) {
    for (auto& data : TransitionData) data.momentumThicknessReynolds = ReThetaT;

    /*--- Wall-normal velocity, and with cross-flow for SA the three components of the vorticity direction. ---*/
    const bool crossFlowSA = options.CrossFlow && TurbModelFamily(config->GetKind_Turb_Model()) == TURB_FAMILY::SA;
    nAuxVar = crossFlowSA ? 4 : 1;
    Grad_AuxVar.resize(nPoint, nAuxVar, nDim, su2double(0.0));
    AuxVar.resize(nPoint, nAuxVar) = su2double(0.0);
    WallNormal.resize(nPoint, 3) = su2double(0.0);
  }

}

void CTransLMVariable::SetIntermittencyEff(unsigned long iPoint, su2double val_Intermittency_sep) {

  /*--- Effective intermittency ---*/
  Intermittency_Eff(iPoint) = max(Solution(iPoint,0), val_Intermittency_sep);

}

void CTransLMVariable::SetIntermittencySep(unsigned long iPoint, su2double val_Intermittency_sep) {
  Intermittency_Sep(iPoint) = val_Intermittency_sep;
}
