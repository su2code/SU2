/*!
 * \file CAdjFlowOutput.cpp
 * \brief Main subroutines for flow discrete adjoint output
 * \author T. Kattmann
 * \version 7.2.1 "Blackbird"
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

#include "../../include/output/CAdjFlowOutput.hpp"

#include "../../../Common/include/CConfig.hpp"

CAdjFlowOutput::CAdjFlowOutput(CConfig* config, unsigned short nDim) : COutput(config, nDim, false) {
  turb_model = config->GetKind_Turb_Model();
  cont_adj = config->GetContinuous_Adjoint();
  frozen_visc = (config->GetFrozen_Visc_Disc() && !cont_adj) || (config->GetFrozen_Visc_Cont() && cont_adj);
}