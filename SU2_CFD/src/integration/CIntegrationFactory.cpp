/*!
 * \file CIntegrationFactory.cpp
 * \brief Main subroutines for CIntegrationFactory .
 * \author T. Albring
 * \version 7.0.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/integration/CIntegrationFactory.hpp"
#include "../../include/integration/CSingleGridIntegration.hpp"
#include "../../include/integration/CMultiGridIntegration.hpp"
#include "../../include/integration/CStructuralIntegration.hpp"
#include "../../include/integration/CFEM_DG_Integration.hpp"

CIntegration** CIntegrationFactory::createIntegrationContainer(ENUM_MAIN_SOLVER kindMainSolver, CConfig *config){

  CIntegration **integration = new CIntegration* [MAX_SOLS]();

  switch (kindMainSolver) {
    case TEMPLATE_SOLVER:
      integration[FLOW_SOL] = createIntegration(SUB_SOLVER::TEMPLATE, config);
      break;
    case INC_EULER:
      integration[FLOW_SOL] = createIntegration(SUB_SOLVER::INC_EULER, config);
      break;
    case EULER:
      integration[FLOW_SOL] = createIntegration(SUB_SOLVER::EULER, config);
      break;
    case INC_NAVIER_STOKES:
      integration[FLOW_SOL] = createIntegration(SUB_SOLVER::INC_NAVIER_STOKES, config);
      integration[HEAT_SOL] = createIntegration(SUB_SOLVER::HEAT, config);
      break;
    case NAVIER_STOKES:
      integration[FLOW_SOL] = createIntegration(SUB_SOLVER::NAVIER_STOKES, config);
      break;
    case RANS:
      integration[FLOW_SOL] = createIntegration(SUB_SOLVER::NAVIER_STOKES, config);
      integration[TURB_SOL] = createIntegration(SUB_SOLVER::TURB, config);
      break;
    case INC_RANS:
      integration[FLOW_SOL] = createIntegration(SUB_SOLVER::INC_NAVIER_STOKES, config);
      integration[HEAT_SOL] = createIntegration(SUB_SOLVER::HEAT, config);
      integration[TURB_SOL] = createIntegration(SUB_SOLVER::TURB, config);
      break;
    case HEAT_EQUATION:
      integration[HEAT_SOL] = createIntegration(SUB_SOLVER::HEAT, config);
      break;
    case ADJ_EULER:
      integration[FLOW_SOL]    = createIntegration(SUB_SOLVER::EULER, config);
      integration[ADJFLOW_SOL] = createIntegration(SUB_SOLVER::CONT_ADJ_EULER, config);
      break;
    case ADJ_NAVIER_STOKES:
      integration[FLOW_SOL]    = createIntegration(SUB_SOLVER::NAVIER_STOKES, config);
      integration[ADJFLOW_SOL] = createIntegration(SUB_SOLVER::CONT_ADJ_NAVIER_STOKES, config);
      break;
    case ADJ_RANS:
      integration[FLOW_SOL]    = createIntegration(SUB_SOLVER::NAVIER_STOKES, config);
      integration[ADJFLOW_SOL] = createIntegration(SUB_SOLVER::CONT_ADJ_NAVIER_STOKES, config);
      integration[TURB_SOL]    = createIntegration(SUB_SOLVER::TURB, config);
      integration[ADJTURB_SOL] = createIntegration(SUB_SOLVER::CONT_ADJ_TURB, config);
      break;
    case DISC_ADJ_EULER:
      integration[FLOW_SOL]    = createIntegration(SUB_SOLVER::EULER, config);
      integration[ADJFLOW_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_FLOW, config);
      break;
    case DISC_ADJ_NAVIER_STOKES:
      integration[FLOW_SOL]    = createIntegration(SUB_SOLVER::NAVIER_STOKES, config);
      integration[ADJFLOW_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_FLOW, config);
      break;
    case DISC_ADJ_RANS:
      integration[FLOW_SOL]    = createIntegration(SUB_SOLVER::NAVIER_STOKES, config);
      integration[ADJFLOW_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_FLOW, config);
      integration[TURB_SOL]    = createIntegration(SUB_SOLVER::TURB, config);
      integration[ADJTURB_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_TURB, config);
      break;
    case DISC_ADJ_INC_EULER:
      integration[FLOW_SOL]    = createIntegration(SUB_SOLVER::INC_EULER, config);
      integration[ADJFLOW_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_FLOW, config);
      break;
    case DISC_ADJ_INC_NAVIER_STOKES:
      integration[FLOW_SOL]    = createIntegration(SUB_SOLVER::INC_NAVIER_STOKES, config);
      integration[ADJFLOW_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_FLOW, config);
      integration[HEAT_SOL]    = createIntegration(SUB_SOLVER::HEAT, config);
      integration[ADJHEAT_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_HEAT, config);
      break;
    case DISC_ADJ_INC_RANS:
      integration[FLOW_SOL]    = createIntegration(SUB_SOLVER::INC_NAVIER_STOKES, config);
      integration[ADJFLOW_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_FLOW, config);
      integration[HEAT_SOL]    = createIntegration(SUB_SOLVER::HEAT, config);
      integration[ADJHEAT_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_HEAT, config);
      integration[TURB_SOL]    = createIntegration(SUB_SOLVER::TURB, config);
      integration[ADJTURB_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_TURB, config);
      break;
    case DISC_ADJ_HEAT:
      integration[HEAT_SOL]    = createIntegration(SUB_SOLVER::HEAT, config);
      integration[ADJHEAT_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_HEAT, config);
      break;
    case FEM_ELASTICITY:
      integration[FEA_SOL] = createIntegration(SUB_SOLVER::FEA, config);
      break;
    case DISC_ADJ_FEM:
      integration[FEA_SOL]    = createIntegration(SUB_SOLVER::FEA, config);
      integration[ADJFEA_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_FEA, config);
      break;
    case FEM_EULER:
      integration[FLOW_SOL] = createIntegration(SUB_SOLVER::DG_EULER, config);
      break;
    case FEM_NAVIER_STOKES: case FEM_LES:
      integration[FLOW_SOL] = createIntegration(SUB_SOLVER::DG_NAVIER_STOKES, config);
      break;
    case FEM_RANS:
      SU2_MPI::Error("FEM RANS not available", CURRENT_FUNCTION);
      break;
    case DISC_ADJ_FEM_EULER:
      integration[FLOW_SOL]    = createIntegration(SUB_SOLVER::DG_EULER, config);
      integration[ADJFLOW_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_FLOW, config);
      break;
    case DISC_ADJ_FEM_NS:
      integration[FLOW_SOL]    = createIntegration(SUB_SOLVER::DG_NAVIER_STOKES, config);
      integration[ADJFLOW_SOL] = createIntegration(SUB_SOLVER::DISC_ADJ_FLOW, config);
      break;
    case DISC_ADJ_FEM_RANS:
      SU2_MPI::Error("Adjoint FEM RANS not available", CURRENT_FUNCTION);
      break;
     default:
      integration = nullptr;
  }

  return integration;
}

CIntegration* CIntegrationFactory::createIntegration(SUB_SOLVER kindSubSolver, CConfig *config){

  CIntegration *integration = nullptr;

  switch(kindSubSolver){
    case SUB_SOLVER::EULER: case SUB_SOLVER::INC_EULER:
    case SUB_SOLVER::NAVIER_STOKES: case SUB_SOLVER::INC_NAVIER_STOKES:
    case SUB_SOLVER::CONT_ADJ_EULER: case SUB_SOLVER::CONT_ADJ_NAVIER_STOKES:
      integration = new CMultiGridIntegration(config);
      break;
    case SUB_SOLVER::TURB: case SUB_SOLVER::TURB_SA: case SUB_SOLVER::TURB_SST:
    case SUB_SOLVER::TRANSITION: case SUB_SOLVER::HEAT: case SUB_SOLVER::CONT_ADJ_TURB:
    case SUB_SOLVER::TEMPLATE:
      integration = new CSingleGridIntegration(config);
      break;
    case SUB_SOLVER::DISC_ADJ_FLOW: case SUB_SOLVER::DISC_ADJ_HEAT: case SUB_SOLVER::DISC_ADJ_TURB:
    case SUB_SOLVER::DISC_ADJ_FEA:
      integration = new CIntegration(config);
      break;
    case SUB_SOLVER::DG_EULER: case SUB_SOLVER::DG_NAVIER_STOKES:
      integration = new CFEM_DG_Integration(config);
      break;
    case SUB_SOLVER::FEA:
      integration = new CStructuralIntegration(config);
      break;
    default:
      SU2_MPI::Error("No proper integration class found for requested sub solver", CURRENT_FUNCTION);
      break;
  }

  return integration;

}
