/*!
 * \file CSingleGridIntegration.cpp
 * \brief Single (fine) grid integration class implementation.
 * \author F. Palacios, T. Economon
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

#include "../../include/integration/CSingleGridIntegration.hpp"


CSingleGridIntegration::CSingleGridIntegration(CConfig *config) : CIntegration(config) { }

void CSingleGridIntegration::SingleGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                                                  CNumerics ******numerics_container, CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) {
  unsigned short iMesh;

  unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

  unsigned short FinestMesh = config[iZone]->GetFinestMesh();

  /*--- Preprocessing ---*/

  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone], FinestMesh, 0, RunTime_EqSystem, false);

  /*--- Set the old solution ---*/

  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Set_OldSolution(geometry[iZone][iInst][FinestMesh]);

  /*--- Time step evaluation ---*/

  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone], FinestMesh, config[iZone]->GetTimeIter());

  /*--- Space integration ---*/

  Space_Integration(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], numerics_container[iZone][iInst][FinestMesh][SolContainer_Position],
                    config[iZone], FinestMesh, NO_RK_ITER, RunTime_EqSystem);

  /*--- Time integration ---*/

  Time_Integration(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone], NO_RK_ITER,
                   RunTime_EqSystem);

  /*--- Postprocessing ---*/

  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone], FinestMesh);

  if (RunTime_EqSystem == RUNTIME_HEAT_SYS) {
    solver_container[iZone][iInst][FinestMesh][HEAT_SOL]->Heat_Fluxes(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone]);
  }

  /*--- If turbulence model, copy the turbulence variables to the coarse levels ---*/

  if (RunTime_EqSystem == RUNTIME_TURB_SYS) {
    for (iMesh = FinestMesh; iMesh < config[iZone]->GetnMGLevels(); iMesh++) {
      SetRestricted_Solution(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position], geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
      SetRestricted_EddyVisc(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position], geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
    }
  }
}

void CSingleGridIntegration::SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse;
  unsigned short iVar, iChildren;
  su2double Area_Parent, Area_Children, *Solution_Fine, *Solution;

  unsigned short nVar = sol_coarse->GetnVar();

  Solution = new su2double[nVar];

  /*--- Compute coarse solution from fine solution ---*/

  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();

    for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;

    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {

      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
    }

    sol_coarse->GetNodes()->SetSolution(Point_Coarse,Solution);

  }

  /*--- MPI the new interpolated solution ---*/

  sol_coarse->InitiateComms(geo_coarse, config, SOLUTION);
  sol_coarse->CompleteComms(geo_coarse, config, SOLUTION);

  delete [] Solution;

}

void CSingleGridIntegration::SetRestricted_EddyVisc(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {

  unsigned long iVertex, Point_Fine, Point_Coarse;
  unsigned short iMarker, iChildren;
  su2double Area_Parent, Area_Children, EddyVisc_Fine, EddyVisc;

  /*--- Compute coarse Eddy Viscosity from fine solution ---*/

  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();

    EddyVisc = 0.0;

    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      EddyVisc_Fine = sol_fine->GetNodes()->GetmuT(Point_Fine);
      EddyVisc += EddyVisc_Fine*Area_Children/Area_Parent;
    }

    sol_coarse->GetNodes()->SetmuT(Point_Coarse,EddyVisc);

  }

  /*--- Update solution at the no slip wall boundary, only the first
   variable (nu_tilde -in SA and SA_NEG- and k -in SST-), to guarantee that the eddy viscoisty
   is zero on the surface ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE     )) {
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->GetNodes()->SetmuT(Point_Coarse,0.0);
      }
    }
  }

  /*--- MPI the new interpolated solution (this also includes the eddy viscosity) ---*/

  sol_coarse->InitiateComms(geo_coarse, config, SOLUTION_EDDY);
  sol_coarse->CompleteComms(geo_coarse, config, SOLUTION_EDDY);

}

