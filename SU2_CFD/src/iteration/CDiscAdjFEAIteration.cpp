/*!
 * \file CDiscAdjFEAIteration.cpp
 * \brief Main subroutines used by SU2_CFD
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

#include "../../include/iteration/CDiscAdjFEAIteration.hpp"
#include "../../include/iteration/CFEAIteration.hpp"
#include "../../include/solvers/CFEASolver.hpp"
#include "../../include/output/COutput.hpp"

CDiscAdjFEAIteration::CDiscAdjFEAIteration(const CConfig *config) : CIteration(config), CurrentRecording(NONE) {
  fem_iteration = new CFEAIteration(config);
}

CDiscAdjFEAIteration::~CDiscAdjFEAIteration() = default;

void CDiscAdjFEAIteration::Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                      CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                      CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                      CFreeFormDefBox*** FFDBox, unsigned short iZone, unsigned short iInst) {
  auto solvers0 = solver[iZone][iInst][MESH_0];
  auto geometry0 = geometry[iZone][iInst][MESH_0];
  auto dirNodes = solvers0[FEA_SOL]->GetNodes();
  auto adjNodes = solvers0[ADJFEA_SOL]->GetNodes();

  auto StoreDirectSolution = [&]() {
    for (auto iPoint = 0ul; iPoint < geometry0->GetnPoint(); iPoint++) {
      adjNodes->SetSolution_Direct(iPoint, dirNodes->GetSolution(iPoint));
    }
  };

  /*--- For the dynamic adjoint, load direct solutions from restart files.
   * For steady, store the direct solution to be able to reset it later. ---*/

  if (config[iZone]->GetTime_Domain()) {
    const int TimeIter = config[iZone]->GetTimeIter();
    const int Direct_Iter = SU2_TYPE::Int(config[iZone]->GetUnst_AdjointIter()) - TimeIter - 1;

    /*--- We need the already converged solution at timesteps n and n-1. On the first
     * adjoint time step we load both, then n-1 becomes "n" and we only load n-2. ---*/

    if (TimeIter > 0) {
      /*--- Save n-1 to become n. ---*/
      for (auto iPoint = 0ul; iPoint < geometry0->GetnPoint(); iPoint++) {
        adjNodes->SetSolution_Direct(iPoint, dirNodes->GetSolution_time_n(iPoint));
      }
    }

    /*--- Load solution at timestep n-1 and push back to correct array. ---*/

    LoadDynamic_Solution(geometry, solver, config, iZone, iInst, Direct_Iter - 1);
    dirNodes->Set_Solution_time_n();

    /*--- Load or set solution at timestep n. ---*/

    if (TimeIter == 0) {
      LoadDynamic_Solution(geometry, solver, config, iZone, iInst, Direct_Iter);
      StoreDirectSolution();
    } else {
      /*--- Set n-1 as n. ---*/
      for (auto iPoint = 0ul; iPoint < geometry0->GetnPoint(); iPoint++)
        for (auto iVar = 0u; iVar < solvers0[ADJFEA_SOL]->GetnVar(); iVar++)
          dirNodes->SetSolution(iPoint, iVar, adjNodes->GetSolution_Direct(iPoint)[iVar]);
    }

  } else {
    /*--- Steady. ---*/
    StoreDirectSolution();
  }

  solvers0[ADJFEA_SOL]->Preprocessing(geometry0, solvers0, config[iZone], MESH_0, 0, RUNTIME_ADJFEA_SYS, false);

}

void CDiscAdjFEAIteration::LoadDynamic_Solution(CGeometry**** geometry, CSolver***** solver, CConfig** config,
                                                unsigned short iZone, unsigned short iInst,
                                                int val_DirectIter) {
  /*--- Set to false to prevent updating Solution_time_n when loading primal solutions of unsteady cases. ---*/
  const bool update_geo = false;
  auto*** solvers = solver[iZone][iInst];

  if (val_DirectIter >= 0) {
    if (rank == MASTER_NODE && iZone == ZONE_0)
      cout << " Loading FEA solution from direct iteration " << val_DirectIter << ".\n";
    solvers[MESH_0][FEA_SOL]->LoadRestart(geometry[iZone][iInst], solvers, config[iZone], val_DirectIter, update_geo);
  } else {
    /*--- If there is no solution file we set the initial conditions. ---*/
    if (rank == MASTER_NODE && iZone == ZONE_0)
      cout << " Setting static conditions at direct iteration " << val_DirectIter << ".\n";
    solvers[MESH_0][FEA_SOL]->SetInitialCondition(geometry[iZone][iInst], solvers, config[iZone], val_DirectIter);
  }
}

void CDiscAdjFEAIteration::IterateDiscAdj(CGeometry**** geometry, CSolver***** solver, CConfig** config,
                                          unsigned short iZone, unsigned short iInst, bool CrossTerm) {

  /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

  solver[iZone][iInst][MESH_0][ADJFEA_SOL]->ExtractAdjoint_Solution(geometry[iZone][iInst][MESH_0], config[iZone],
                                                                    CrossTerm);

  solver[iZone][iInst][MESH_0][ADJFEA_SOL]->ExtractAdjoint_Variables(geometry[iZone][iInst][MESH_0], config[iZone]);
}

void CDiscAdjFEAIteration::RegisterInput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                         unsigned short iZone, unsigned short iInst, RECORDING kind_recording) {
  if (kind_recording != RECORDING::MESH_COORDS) {
    /*--- Register structural displacements as input ---*/

    solver[iZone][iInst][MESH_0][ADJFEA_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

    /*--- Register variables as input ---*/

    solver[iZone][iInst][MESH_0][ADJFEA_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
  } else {
    /*--- Register topology optimization densities (note direct solver) ---*/

    solver[iZone][iInst][MESH_0][FEA_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);

    /*--- Register mesh coordinates for geometric sensitivities ---*/

    geometry[iZone][iInst][MESH_0]->RegisterCoordinates();
  }
}

void CDiscAdjFEAIteration::SetDependencies(CSolver***** solver, CGeometry**** geometry, CNumerics****** numerics,
                                           CConfig** config, unsigned short iZone, unsigned short iInst,
                                           RECORDING kind_recording) {
  auto dir_solver = solver[iZone][iInst][MESH_0][FEA_SOL];
  auto adj_solver = solver[iZone][iInst][MESH_0][ADJFEA_SOL];
  auto structural_geometry = geometry[iZone][iInst][MESH_0];
  auto structural_numerics = numerics[iZone][iInst][MESH_0][FEA_SOL];

  /*--- Some numerics are only instanciated under these conditions ---*/
  const bool fsi = config[iZone]->GetFSI_Simulation() || config[iZone]->GetMultizone_Problem();
  const bool nonlinear = config[iZone]->GetGeometricConditions() == STRUCT_DEFORMATION::LARGE;
  const bool de_effects = config[iZone]->GetDE_Effects() && nonlinear;
  const bool element_based = dir_solver->IsElementBased() && nonlinear;

  SU2_OMP_PARALLEL
  {

  const int thread = omp_get_thread_num();
  const int offset = thread*MAX_TERMS;
  const int fea_term = FEA_TERM+offset;
  const int mat_nhcomp = MAT_NHCOMP+offset;
  const int mat_idealde = MAT_IDEALDE+offset;
  const int mat_knowles = MAT_KNOWLES+offset;
  const int de_term = DE_TERM+offset;

  for (unsigned short iProp = 0; iProp < config[iZone]->GetnElasticityMat(); iProp++) {
    su2double E = adj_solver->GetVal_Young(iProp);
    su2double nu = adj_solver->GetVal_Poisson(iProp);
    su2double rho = adj_solver->GetVal_Rho(iProp);
    su2double rhoDL = adj_solver->GetVal_Rho_DL(iProp);

    /*--- Add dependencies for E and Nu ---*/

    structural_numerics[fea_term]->SetMaterial_Properties(iProp, E, nu);

    /*--- Add dependencies for Rho and Rho_DL ---*/

    structural_numerics[fea_term]->SetMaterial_Density(iProp, rho, rhoDL);

    /*--- Add dependencies for element-based simulations. ---*/

    if (element_based) {
      /*--- Neo Hookean Compressible ---*/
      structural_numerics[mat_nhcomp]->SetMaterial_Properties(iProp, E, nu);
      structural_numerics[mat_nhcomp]->SetMaterial_Density(iProp, rho, rhoDL);

      /*--- Ideal DE ---*/
      structural_numerics[mat_idealde]->SetMaterial_Properties(iProp, E, nu);
      structural_numerics[mat_idealde]->SetMaterial_Density(iProp, rho, rhoDL);

      /*--- Knowles ---*/
      structural_numerics[mat_knowles]->SetMaterial_Properties(iProp, E, nu);
      structural_numerics[mat_knowles]->SetMaterial_Density(iProp, rho, rhoDL);
    }
  }

  if (de_effects) {
    for (unsigned short iEField = 0; iEField < adj_solver->GetnEField(); iEField++) {
      structural_numerics[fea_term]->Set_ElectricField(iEField, adj_solver->GetVal_EField(iEField));
      structural_numerics[de_term]->Set_ElectricField(iEField, adj_solver->GetVal_EField(iEField));
    }
  }

  /*--- Add dependencies for element-based simulations. ---*/

  switch (config[iZone]->GetDV_FEA()) {
    case YOUNG_MODULUS:
    case POISSON_RATIO:
    case DENSITY_VAL:
    case DEAD_WEIGHT:
    case ELECTRIC_FIELD:

      for (unsigned short iDV = 0; iDV < adj_solver->GetnDVFEA(); iDV++) {
        su2double dvfea = adj_solver->GetVal_DVFEA(iDV);

        structural_numerics[fea_term]->Set_DV_Val(iDV, dvfea);

        if (de_effects) structural_numerics[de_term]->Set_DV_Val(iDV, dvfea);

        if (element_based) {
          structural_numerics[mat_nhcomp]->Set_DV_Val(iDV, dvfea);
          structural_numerics[mat_idealde]->Set_DV_Val(iDV, dvfea);
          structural_numerics[mat_knowles]->Set_DV_Val(iDV, dvfea);
        }
      }
      break;
  }

  /*--- MPI dependencies. ---*/

  dir_solver->InitiateComms(structural_geometry, config[iZone], SOLUTION_FEA);
  dir_solver->CompleteComms(structural_geometry, config[iZone], SOLUTION_FEA);

  structural_geometry->InitiateComms(structural_geometry, config[iZone], COORDINATES);
  structural_geometry->CompleteComms(structural_geometry, config[iZone], COORDINATES);

  }
  END_SU2_OMP_PARALLEL

  /*--- FSI specific dependencies. ---*/
  if (fsi) {
    /*--- Set relation between solution and predicted displacements, which are the transferred ones. ---*/
    dir_solver->PredictStruct_Displacement(structural_geometry, config[iZone]);
  }

  /*--- Topology optimization dependencies. ---*/

  /*--- We only differentiate wrt to this variable in the adjoint secondary recording. ---*/
  if (config[iZone]->GetTopology_Optimization() && (kind_recording == RECORDING::MESH_COORDS)) {
    /*--- The filter may require the volumes of the elements. ---*/
    structural_geometry->SetElemVolume();
    /// TODO: Ideally there would be a way to capture this dependency without the `static_cast`, but
    ///       making it a virtual method of CSolver does not feel "right" as its purpose could be confused.
    static_cast<CFEASolver*>(dir_solver)->FilterElementDensities(structural_geometry, config[iZone]);
  }

}

void CDiscAdjFEAIteration::RegisterOutput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                          unsigned short iZone, unsigned short iInst) {
  /*--- Register conservative variables as output of the iteration ---*/

  solver[iZone][iInst][MESH_0][ADJFEA_SOL]->RegisterOutput(geometry[iZone][iInst][MESH_0], config[iZone]);
}

void CDiscAdjFEAIteration::InitializeAdjoint(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                             unsigned short iZone, unsigned short iInst) {
  /*--- Initialize the adjoints the conservative variables ---*/

  AD::ResizeAdjoints();
  solver[iZone][iInst][MESH_0][ADJFEA_SOL]->SetAdjoint_Output(geometry[iZone][iInst][MESH_0], config[iZone]);
}

bool CDiscAdjFEAIteration::Monitor(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                   CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                   CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                   CFreeFormDefBox*** FFDBox, unsigned short iZone, unsigned short iInst) {
  /*--- Write the convergence history (only screen output) ---*/

  output->SetHistoryOutput(geometry[iZone][INST_0][MESH_0], solver[iZone][INST_0][MESH_0], config[iZone],
                            config[iZone]->GetTimeIter(), config[iZone]->GetOuterIter(),
                            config[iZone]->GetInnerIter());

  return output->GetConvergence();
}

void CDiscAdjFEAIteration::Postprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                       CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                       CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                       CFreeFormDefBox*** FFDBox, unsigned short iZone, unsigned short iInst) {
  auto solvers0 = solver[iZone][iInst][MESH_0];

  // TEST: for implementation of python framework in standalone structural problems
  if (config[iZone]->GetAdvanced_FEAElementBased() && (rank == MASTER_NODE)) {
    /*--- Header of the temporary output file ---*/
    ofstream myfile_res;
    bool outputDVFEA = true;

    switch (config[iZone]->GetDV_FEA()) {
      case YOUNG_MODULUS:
        myfile_res.open("grad_young.opt");
        break;
      case POISSON_RATIO:
        myfile_res.open("grad_poisson.opt");
        break;
      case DENSITY_VAL:
      case DEAD_WEIGHT:
        myfile_res.open("grad_density.opt");
        break;
      case ELECTRIC_FIELD:
        myfile_res.open("grad_efield.opt");
        break;
      default:
        outputDVFEA = false;
        break;
    }

    if (outputDVFEA) {
      const auto nDV = solvers0[ADJFEA_SOL]->GetnDVFEA();
      myfile_res << "INDEX\tGRAD\n";
      myfile_res.precision(15);

      for (auto iDV = 0u; iDV < nDV; iDV++) {
        myfile_res << iDV << "\t";
        myfile_res << scientific << solvers0[ADJFEA_SOL]->GetTotal_Sens_DVFEA(iDV) << "\n";
      }
    }
  }
}
