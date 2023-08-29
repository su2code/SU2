/*!
 * \file CIntegrationFactory.cpp
 * \brief Main subroutines for CIntegrationFactory .
 * \author T. Albring
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../include/integration/CNewtonIntegration.hpp"
#include "../../include/integration/CStructuralIntegration.hpp"
#include "../../include/integration/CFEM_DG_Integration.hpp"

CIntegration** CIntegrationFactory::CreateIntegrationContainer(MAIN_SOLVER kindMainSolver,
                                                               const CSolver* const* solver_container){

  auto **integration = new CIntegration* [MAX_SOLS]();

  for (unsigned int iSol = 0; iSol < MAX_SOLS; iSol++){
    if (solver_container[iSol] != nullptr){
      const SolverMetaData &solverInfo = CSolverFactory::GetSolverMeta(solver_container[iSol]);
      integration[iSol] = CreateIntegration(solverInfo.integrationType);
    }
  }

  return integration;
}

CIntegration* CIntegrationFactory::CreateIntegration(INTEGRATION_TYPE integrationType){

  CIntegration *integration = nullptr;

  switch(integrationType){
    case INTEGRATION_TYPE::DEFAULT:
      integration = new CIntegration();
      break;
    case INTEGRATION_TYPE::SINGLEGRID:
      integration = new CSingleGridIntegration();
      break;
    case INTEGRATION_TYPE::MULTIGRID:
      integration = new CMultiGridIntegration();
      break;
    case INTEGRATION_TYPE::NEWTON:
      integration = new CNewtonIntegration();
      break;
    case INTEGRATION_TYPE::STRUCTURAL:
      integration = new CStructuralIntegration();
      break;
    case INTEGRATION_TYPE::FEM_DG:
      integration = new CFEM_DG_Integration();
      break;
    case INTEGRATION_TYPE::NONE:
      integration = nullptr;
      break;
  }

  return integration;

}
