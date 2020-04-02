/*!
 * \file CStructuralIntegration.cpp
 * \brief Space and time integration for structural problems.
 * \author F. Palacios, T. Economon
 * \version 7.0.3 "Blackbird"
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
                         solver_container[iZone][iInst][MESH_0],
                         config[iZone],
                         numerics_container[iZone][iInst][MESH_0][SolContainer_Position],
                         MESH_0);
}

void CStructuralIntegration::Space_Integration_FEM(CGeometry *geometry, CSolver **solver_container,
                                                   CNumerics **numerics, CConfig *config,
                                                   unsigned short RunTime_EqSystem) {
  bool dynamic = config->GetTime_Domain();
  bool first_iter = (config->GetInnerIter() == 0);
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);
  unsigned short IterativeScheme = config->GetKind_SpaceIteScheme_FEA();

  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);
  CSolver* solver = solver_container[MainSolver];

  /*--- Initial calculation, different logic for restarted simulations. ---*/
  bool initial_calc = false;
  if (config->GetRestart())
    initial_calc = (config->GetTimeIter() == config->GetRestart_Iter()) && first_iter;
  else
    initial_calc = (config->GetTimeIter() == 0) && first_iter;

  /*--- Mass Matrix computed during preprocessing, see notes therein. ---*/

  /*--- If the analysis is linear, only a the constitutive term of the stiffness matrix has to be computed. ---*/
  /*--- This is done only once, at the beginning of the calculation. From then on, K is constant. ---*/
  /*--- For correct differentiation of dynamic cases the matrix needs to be computed every time. ---*/
  if (linear_analysis && (dynamic || initial_calc))
    solver->Compute_StiffMatrix(geometry, numerics, config);

  if (nonlinear_analysis) {

    /*--- If the analysis is nonlinear the stress terms also need to be computed. ---*/
    /*--- For full Newton-Raphson the stiffness matrix and the nodal term are updated every time. ---*/
    if (IterativeScheme == NEWTON_RAPHSON) {
      solver->Compute_StiffMatrix_NodalStressRes(geometry, numerics, config);
    }

    /*--- If the method is modified Newton-Raphson, the stiffness matrix is only computed once at the beginning
     * of the time step, then only the Nodal Stress Term has to be computed on each iteration. ---*/
    if (IterativeScheme == MODIFIED_NEWTON_RAPHSON) {
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
        solver->BC_Dir_Load(geometry, numerics[FEA_TERM], config, iMarker);
        break;

      case LOAD_SINE_BOUNDARY:
        solver->BC_Sine_Load(geometry, numerics[FEA_TERM], config, iMarker);
        break;

      case LOAD_BOUNDARY:
        solver->BC_Normal_Load(geometry, numerics[FEA_TERM], config, iMarker);
        break;

      case DAMPER_BOUNDARY:
        solver->BC_Damper(geometry, numerics[FEA_TERM], config, iMarker);
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
    case (CD_EXPLICIT):
      solver_container[MainSolver]->ImplicitNewmark_Iteration(geometry, solver_container, config);
      break;
    case (NEWMARK_IMPLICIT):
      solver_container[MainSolver]->ImplicitNewmark_Iteration(geometry, solver_container, config);
      break;
    case (GENERALIZED_ALPHA):
      solver_container[MainSolver]->GeneralizedAlpha_Iteration(geometry, solver_container, config);
      break;
  }

  /*--- Apply ESSENTIAL BOUNDARY CONDITIONS ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case CLAMPED_BOUNDARY:
        solver_container[MainSolver]->BC_Clamped(geometry, numerics[FEA_TERM], config, iMarker);
        break;
      case DISP_DIR_BOUNDARY:
        solver_container[MainSolver]->BC_DispDir(geometry, numerics[FEA_TERM], config, iMarker);
        break;
      case DISPLACEMENT_BOUNDARY:
        solver_container[MainSolver]->BC_Normal_Displacement(geometry, numerics[CONV_BOUND_TERM], config, iMarker);
        break;
    }
  }

  /*--- Solver linearized system ---*/

  solver_container[MainSolver]->Solve_System(geometry, config);

  /*--- Update solution ---*/

  switch (config->GetKind_TimeIntScheme_FEA()) {
    case (CD_EXPLICIT):
      solver_container[MainSolver]->ImplicitNewmark_Update(geometry, solver_container, config);
      break;
    case (NEWMARK_IMPLICIT):
      solver_container[MainSolver]->ImplicitNewmark_Update(geometry, solver_container, config);
      break;
    case (GENERALIZED_ALPHA):
      solver_container[MainSolver]->GeneralizedAlpha_UpdateDisp(geometry, solver_container, config);
      break;
  }

  /*--- Reinforce ESSENTIAL BOUNDARY CONDITIONS: avoids accumulation of numerical error ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case CLAMPED_BOUNDARY:
        solver_container[MainSolver]->BC_Clamped_Post(geometry, numerics[FEA_TERM], config, iMarker);
        break;
    }
  }

  /*--- Perform the MPI communication of the solution ---*/

  solver_container[MainSolver]->InitiateComms(geometry, config, SOLUTION_FEA);
  solver_container[MainSolver]->CompleteComms(geometry, config, SOLUTION_FEA);

}
