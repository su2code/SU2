/*!
 * \file CAdjFluidIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 7.1.0 "Blackbird"
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

#include "../../include/iteration/CIterationFactory.hpp"
#include "../../include/iteration/CIteration.hpp"
#include "../../include/iteration/CAdjFluidIteration.hpp"
#include "../../include/iteration/CDiscAdjFEAIteration.hpp"
#include "../../include/iteration/CDiscAdjFluidIteration.hpp"
#include "../../include/iteration/CDiscAdjHeatIteration.hpp"
#include "../../include/iteration/CFluidIteration.hpp"
#include "../../include/iteration/CFEMFluidIteration.hpp"
#include "../../include/iteration/CTurboIteration.hpp"
#include "../../include/iteration/CHeatIteration.hpp"
#include "../../include/iteration/CFEAIteration.hpp"

CIteration* CIterationFactory::CreateIteration(ENUM_MAIN_SOLVER kindSolver, const CConfig* config){

  CIteration *iteration = nullptr;

  const auto rank = SU2_MPI::GetRank();

  /*--- Loop over all zones and instantiate the physics iteration. ---*/

  switch (kindSolver) {

    case EULER: case NAVIER_STOKES: case RANS:
    case INC_EULER: case INC_NAVIER_STOKES: case INC_RANS:
    case NEMO_EULER: case NEMO_NAVIER_STOKES:
      if(config->GetBoolTurbomachinery()){
        if (rank == MASTER_NODE)
          cout << "Euler/Navier-Stokes/RANS turbomachinery fluid iteration." << endl;
        iteration = new CTurboIteration(config);

      }
      else{
        if (rank == MASTER_NODE)
          cout << "Euler/Navier-Stokes/RANS fluid iteration." << endl;
        iteration = new CFluidIteration(config);
      }
      break;

    case FEM_EULER: case FEM_NAVIER_STOKES: case FEM_RANS: case FEM_LES:
      if (rank == MASTER_NODE)
        cout << "Finite element Euler/Navier-Stokes/RANS/LES flow iteration." << endl;
      iteration = new CFEMFluidIteration(config);
      break;

    case HEAT_EQUATION:
      if (rank == MASTER_NODE)
        cout << "Heat iteration (finite volume method)." << endl;
      iteration = new CHeatIteration(config);
      break;

    case FEM_ELASTICITY:
      if (rank == MASTER_NODE)
        cout << "FEM iteration." << endl;
      iteration = new CFEAIteration(config);
      break;

    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
      if (rank == MASTER_NODE)
        cout << "Adjoint Euler/Navier-Stokes/RANS fluid iteration." << endl;
      iteration = new CAdjFluidIteration(config);
      break;

    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
      if (rank == MASTER_NODE)
        cout << "Discrete adjoint Euler/Navier-Stokes/RANS fluid iteration." << endl;
      iteration = new CDiscAdjFluidIteration(config);
      break;

    case DISC_ADJ_FEM_EULER : case DISC_ADJ_FEM_NS : case DISC_ADJ_FEM_RANS :
      if (rank == MASTER_NODE)
        cout << "Discrete adjoint finite element Euler/Navier-Stokes/RANS fluid iteration." << endl;
      iteration = new CDiscAdjFluidIteration(config);
      break;

    case DISC_ADJ_FEM:
      if (rank == MASTER_NODE)
        cout << "Discrete adjoint FEM structural iteration." << endl;
      iteration = new CDiscAdjFEAIteration(config);
      break;

    case DISC_ADJ_HEAT:
      if (rank == MASTER_NODE)
        cout << "Discrete adjoint heat iteration." << endl;
      iteration = new CDiscAdjHeatIteration(config);
      break;

    case NO_SOLVER: case TEMPLATE_SOLVER: case MULTIPHYSICS:
      SU2_MPI::Error("No iteration found for specified solver.", CURRENT_FUNCTION);
      break;
  }

  return iteration;
}
