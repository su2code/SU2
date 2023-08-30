/*!
 * \file CStructuralIntegration.cpp
 * \brief Space and time integration for structural problems.
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

#include "../../include/integration/CStructuralIntegration.hpp"


CStructuralIntegration::CStructuralIntegration() : CIntegration() { }

void CStructuralIntegration::Structural_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                                                  CNumerics ******numerics_container, CConfig **config,
                                                  unsigned short RunTime_EqSystem, unsigned short iZone,
                                                  unsigned short iInst) {

  unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

  CSolver* solver = solver_container[iZone][iInst][MESH_0][SolContainer_Position];


  solver->Preprocessing(geometry[iZone][iInst][MESH_0],
                        solver_container[iZone][iInst][MESH_0],
                        config[iZone],
                        numerics_container[iZone][iInst][MESH_0][SolContainer_Position],
                        MESH_0,
                        NO_RK_ITER,
                        RunTime_EqSystem,
                        false);

  Space_Integration_FEM(geometry[iZone][iInst][MESH_0],
                        solver_container[iZone][iInst][MESH_0],
                        numerics_container[iZone][iInst][MESH_0][SolContainer_Position],
                        config[iZone],
                        RunTime_EqSystem);

  Time_Integration_FEM(geometry[iZone][iInst][MESH_0],
                       solver_container[iZone][iInst][MESH_0],
                       numerics_container[iZone][iInst][MESH_0][SolContainer_Position],
                       config[iZone],
                       RunTime_EqSystem);

  solver->Postprocessing(geometry[iZone][iInst][MESH_0],
                         config[iZone],
                         numerics_container[iZone][iInst][MESH_0][SolContainer_Position]);
}

void CStructuralIntegration::Space_Integration_FEM(CGeometry *geometry, CSolver **solver_container,
                                                   CNumerics **numerics, CConfig *config,
                                                   unsigned short RunTime_EqSystem) {
  const bool first_iter = (config->GetInnerIter() == 0);
  const bool linear_analysis = (config->GetGeometricConditions() == STRUCT_DEFORMATION::SMALL);
  const auto IterativeScheme = config->GetKind_SpaceIteScheme_FEA();

  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);
  CSolver* solver = solver_container[MainSolver];

  /*--- Mass Matrix was computed during preprocessing, see notes therein. ---*/

  if (linear_analysis) {
    /*--- If the analysis is linear, only a the constitutive term of the stiffness matrix has to be computed. ---*/
    solver->Compute_StiffMatrix(geometry, numerics, config);
  }
  else {
    /*--- If the analysis is nonlinear the stress terms also need to be computed. ---*/
    /*--- For full Newton-Raphson the stiffness matrix and the nodal term are updated every time. ---*/
    if (IterativeScheme == STRUCT_SPACE_ITE::NEWTON) {
      solver->Compute_StiffMatrix_NodalStressRes(geometry, numerics, config);
    }

    /*--- If the method is modified Newton-Raphson, the stiffness matrix is only computed once at the beginning
     * of the time step, then only the Nodal Stress Term has to be computed on each iteration. ---*/
    if (IterativeScheme == STRUCT_SPACE_ITE::MOD_NEWTON) {
      if (first_iter)
        solver->Compute_StiffMatrix_NodalStressRes(geometry, numerics, config);
      else
        solver->Compute_NodalStressRes(geometry, numerics, config);
    }
  }

  /*--- Apply the NATURAL BOUNDARY CONDITIONS (loads). ---*/
  /*--- If there are FSI loads, they have to be previously applied at other level involving both zones. ---*/

  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    switch (config->GetMarker_All_KindBC(iMarker)) {

      case LOAD_DIR_BOUNDARY:
        solver->BC_Dir_Load(geometry, config, iMarker);
        break;

      case LOAD_SINE_BOUNDARY:
        solver->BC_Sine_Load(geometry, config, iMarker);
        break;

      case LOAD_BOUNDARY:
        solver->BC_Normal_Load(geometry, config, iMarker);
        break;

      case DAMPER_BOUNDARY:
        solver->BC_Damper(geometry, config, iMarker);
        break;
    }
  }

}

void CStructuralIntegration::Time_Integration_FEM(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics,
                                                  CConfig *config, unsigned short RunTime_EqSystem) {
  unsigned short iMarker;

  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);

  /*--- Set the Jacobian according to the different time integration methods ---*/

  switch (config->GetKind_TimeIntScheme_FEA()) {
    case (STRUCT_TIME_INT::NEWMARK_IMPLICIT):
      solver_container[MainSolver]->ImplicitNewmark_Iteration(geometry, numerics, config);
      break;
    case (STRUCT_TIME_INT::GENERALIZED_ALPHA):
      solver_container[MainSolver]->GeneralizedAlpha_Iteration(geometry, numerics, config);
      break;
  }

  /*--- Apply ESSENTIAL BOUNDARY CONDITIONS ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case CLAMPED_BOUNDARY:
        solver_container[MainSolver]->BC_Clamped(geometry, config, iMarker);
        break;
      case SYMMETRY_PLANE:
        solver_container[MainSolver]->BC_Sym_Plane(geometry, config, iMarker);
        break;
      case DISP_DIR_BOUNDARY:
        solver_container[MainSolver]->BC_DispDir(geometry, config, iMarker);
        break;
    }
  }

  /*--- Solve linear system ---*/

  solver_container[MainSolver]->Solve_System(geometry, config);

  /*--- Update solution ---*/

  switch (config->GetKind_TimeIntScheme_FEA()) {
    case (STRUCT_TIME_INT::NEWMARK_IMPLICIT):
      solver_container[MainSolver]->ImplicitNewmark_Update(geometry, config);
      break;
    case (STRUCT_TIME_INT::GENERALIZED_ALPHA):
      solver_container[MainSolver]->GeneralizedAlpha_UpdateDisp(geometry, config);
      break;
  }

  /*--- Reinforce ESSENTIAL BOUNDARY CONDITIONS: avoids accumulation of numerical error ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case CLAMPED_BOUNDARY:
        solver_container[MainSolver]->BC_Clamped_Post(geometry, config, iMarker);
        break;
    }
  }

}

void CStructuralIntegration::SetDualTime_Solver(const CGeometry *geometry, CSolver *solver, const CConfig *config, unsigned short iMesh) {

  const bool fsi = config->GetFSI_Simulation();

  /*--- Update the solution according to the integration scheme used ---*/

  switch (config->GetKind_TimeIntScheme_FEA()) {
    case (STRUCT_TIME_INT::NEWMARK_IMPLICIT):
      if (fsi && config->GetRelaxation()) solver->ImplicitNewmark_Relaxation(geometry, config);
      break;
    case (STRUCT_TIME_INT::GENERALIZED_ALPHA):
      solver->GeneralizedAlpha_UpdateSolution(geometry, config);
      solver->GeneralizedAlpha_UpdateLoads(geometry, config);
      break;
  }

  /*--- Store the solution at t+1 as solution at t, both for the local points and for the halo points ---*/

  solver->GetNodes()->Set_Solution_time_n();

  /*--- If FSI problem, save the last Aitken relaxation parameter of the previous time step ---*/

  if (fsi) solver->SetWAitken_Dyn_tn1();
}
