/*!
 * \file COutputFactory.cpp
 * \brief Main subroutines for output solver information
 * \author T. Albring
 * \version 7.0.4 "Blackbird"
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

#include "../../include/output/COutputFactory.hpp"
#include "../../include/output/COutputLegacy.hpp"
#include "../../include/output/COutput.hpp"
#include "../../include/output/CMultizoneOutput.hpp"
#include "../../include/output/CElasticityOutput.hpp"
#include "../../include/output/CAdjElasticityOutput.hpp"
#include "../../include/output/CFlowCompOutput.hpp"
#include "../../include/output/CAdjFlowOutput.hpp"
#include "../../include/output/CFlowCompFEMOutput.hpp"
#include "../../include/output/CFlowIncOutput.hpp"
#include "../../include/output/CAdjFlowIncOutput.hpp"
#include "../../include/output/CHeatOutput.hpp"
#include "../../include/output/CAdjHeatOutput.hpp"

COutput* COutputFactory::createOutput(ENUM_MAIN_SOLVER kindSolver, CConfig* config, int nDim){

  COutput* output = nullptr;

  switch(kindSolver){
    case EULER: case NAVIER_STOKES: case RANS:
      output = new CFlowCompOutput(config, nDim);
      break;
    case INC_EULER: case INC_NAVIER_STOKES: case INC_RANS:
      output = new CFlowIncOutput(config, nDim);
      break;
    case HEAT_EQUATION:
      output = new CHeatOutput(config, nDim);
      break;
    case FEM_ELASTICITY:
      output = new CElasticityOutput(config, nDim);
      break;
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
      output = new CAdjFlowCompOutput(config, nDim);
      break;
    case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
      output = new CAdjFlowIncOutput(config, nDim);
      break;
    case DISC_ADJ_FEM:
      output = new CAdjElasticityOutput(config, nDim);
      break;
    case DISC_ADJ_HEAT:
      output = new CAdjHeatOutput(config, nDim);
      break;
    case FEM_EULER: case FEM_LES: case FEM_RANS: case FEM_NAVIER_STOKES:
      output = new CFlowCompFEMOutput(config, nDim);
      break;
    default:
      output = new COutput(config, nDim, false);
      break;
  }

  return output;
}

COutput* COutputFactory::createMultizoneOutput(CConfig *driverConfig, CConfig** config_container, int nDim){

  COutput* output = new CMultizoneOutput(driverConfig, config_container, nDim);

  return output;
}

COutputLegacy* COutputFactory::createLegacyOutput(CConfig *config){

  COutputLegacy* output = new COutputLegacy(config);

  return output;

}

