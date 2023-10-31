/*!
 * \file CSingleGridIntegration.cpp
 * \brief Single (fine) grid integration class implementation.
 * \author F. Palacios, T. Economon
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

#include "../../include/integration/CSingleGridIntegration.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"


CSingleGridIntegration::CSingleGridIntegration() : CIntegration() { }

void CSingleGridIntegration::SingleGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                                                  CNumerics ******numerics_container, CConfig **config,
                                                  unsigned short RunTime_EqSystem, unsigned short iZone,
                                                  unsigned short iInst) {

  const unsigned short Solver_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

  /*--- Start an OpenMP parallel region covering the entire iteration. ---*/

  SU2_OMP_PARALLEL_(if(solver_container[iZone][iInst][MESH_0][Solver_Position]->GetHasHybridParallel()))
  {

  unsigned short FinestMesh = config[iZone]->GetFinestMesh();

  CGeometry* geometry_fine = geometry[iZone][iInst][FinestMesh];
  CSolver** solvers_fine = solver_container[iZone][iInst][FinestMesh];

  /*--- Preprocessing ---*/

  solvers_fine[Solver_Position]->Preprocessing(geometry_fine, solvers_fine, config[iZone],
                                               FinestMesh, 0, RunTime_EqSystem, false);

  /*--- Set the old solution ---*/

  solvers_fine[Solver_Position]->Set_OldSolution();

  /*--- Time step evaluation ---*/

  solvers_fine[Solver_Position]->SetTime_Step(geometry_fine, solvers_fine, config[iZone],
                                              FinestMesh, config[iZone]->GetTimeIter());

  /*--- Space integration ---*/

  Space_Integration(geometry_fine, solvers_fine,
                    numerics_container[iZone][iInst][FinestMesh][Solver_Position],
                    config[iZone], FinestMesh, NO_RK_ITER, RunTime_EqSystem);

  /*--- Time integration ---*/

  Time_Integration(geometry_fine, solvers_fine, config[iZone], NO_RK_ITER, RunTime_EqSystem);

  /*--- Postprocessing ---*/

  solvers_fine[Solver_Position]->Postprocessing(geometry_fine, solvers_fine, config[iZone], FinestMesh);

  if (RunTime_EqSystem == RUNTIME_HEAT_SYS) {
    SU2_OMP_SAFE_GLOBAL_ACCESS(solvers_fine[HEAT_SOL]->Heat_Fluxes(geometry_fine, solvers_fine, config[iZone]);)
  }

  /*--- If turbulence model, copy the turbulence variables to the coarse levels ---*/

  if (RunTime_EqSystem == RUNTIME_TURB_SYS) {

    for (unsigned short iMesh = FinestMesh; iMesh < config[iZone]->GetnMGLevels(); iMesh++) {

      SetRestricted_Solution(RunTime_EqSystem,
                             solver_container[iZone][iInst][iMesh][Solver_Position],
                             solver_container[iZone][iInst][iMesh+1][Solver_Position],
                             geometry[iZone][iInst][iMesh],
                             geometry[iZone][iInst][iMesh+1],
                             config[iZone]);

      SetRestricted_EddyVisc(RunTime_EqSystem,
                             solver_container[iZone][iInst][iMesh][Solver_Position],
                             solver_container[iZone][iInst][iMesh+1][Solver_Position],
                             geometry[iZone][iInst][iMesh],
                             geometry[iZone][iInst][iMesh+1],
                             config[iZone]);
    }

  }

  }
  END_SU2_OMP_PARALLEL
}

void CSingleGridIntegration::SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                    CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  CSolver::MultigridRestriction(*geo_fine, sol_fine->GetNodes()->GetSolution(),
                                *geo_coarse, sol_coarse->GetNodes()->GetSolution());
  sol_coarse->InitiateComms(geo_coarse, config, SOLUTION);
  sol_coarse->CompleteComms(geo_coarse, config, SOLUTION);
}

void CSingleGridIntegration::SetRestricted_EddyVisc(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                    CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {

  unsigned long iVertex, Point_Fine, Point_Coarse;
  unsigned short iMarker, iChildren;
  su2double Area_Parent, Area_Children, EddyVisc_Fine, EddyVisc;

  /*--- Compute coarse Eddy Viscosity from fine solution ---*/

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

    Area_Parent = geo_coarse->nodes->GetVolume(Point_Coarse);

    EddyVisc = 0.0;

    for (iChildren = 0; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      Area_Children = geo_fine->nodes->GetVolume(Point_Fine);
      EddyVisc_Fine = sol_fine->GetNodes()->GetmuT(Point_Fine);
      EddyVisc += EddyVisc_Fine*Area_Children/Area_Parent;
    }

    sol_coarse->GetNodes()->SetmuT(Point_Coarse,EddyVisc);

  }
  END_SU2_OMP_FOR

  /*--- Update solution at the no slip wall boundary, only the first
   variable (nu_tilde -in SA and SA_NEG- and k -in SST-), to guarantee that the eddy viscoisty
   is zero on the surface ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetViscous_Wall(iMarker)) {
      SU2_OMP_FOR_STAT(32)
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->GetNodes()->SetmuT(Point_Coarse,0.0);
      }
      END_SU2_OMP_FOR
    }
  }

  /*--- MPI the new interpolated solution (this also includes the eddy viscosity) ---*/

  sol_coarse->InitiateComms(geo_coarse, config, SOLUTION_EDDY);
  sol_coarse->CompleteComms(geo_coarse, config, SOLUTION_EDDY);

}
