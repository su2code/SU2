/*!
 * \file driver_oneshot_singlezone.cpp
 * \brief The main subroutines for driving adjoint single-zone problems.
 * \author E. C. Bunschoten
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

#include "../../include/drivers/COneShotSinglezoneDriver.hpp"
#include "../../include/output/tools/CWindowingTools.hpp"
#include "../../include/output/COutputFactory.hpp"
#include "../../include/output/COutput.hpp"
#include "../../include/iteration/CIterationFactory.hpp"
#include "../../include/iteration/CTurboIteration.hpp"
#include "../../../Common/include/toolboxes/CQuasiNewtonInvLeastSquares.hpp"
#include "../../include/solvers/CMeshSolver.hpp"
#include "../../include/solvers/CDiscAdjMeshSolver.hpp"

COneShotSinglezoneDriver::COneShotSinglezoneDriver(char* confFile,
             unsigned short val_nZone,
             SU2_Comm MPICommunicator) : CDiscAdjSinglezoneDriver(confFile,val_nZone,MPICommunicator){

  solver[MESH_SOL] = new CMeshSolver(geometry, config);
  solver[ADJMESH_SOL] = new CDiscAdjMeshSolver(geometry, config, solver[MESH_SOL]);
  surface_movement[ZONE_0] = new CSurfaceMovement();

  unsigned short nDV = config_container[ZONE_0]->GetnDV();
  Gradient = new su2double*[nDV];

  for (auto iDV = 0u; iDV < nDV; iDV++) {
    /*--- Initialize to zero ---*/
    unsigned short nValue = config_container[ZONE_0]->GetnDV_Value(iDV);

    Gradient[iDV] = new su2double[nValue];
  }

  cout << static_cast<int>(config->GetKind_Solver()) << endl;

}

COneShotSinglezoneDriver::~COneShotSinglezoneDriver() {

  delete direct_iteration;
  delete direct_output;

}

void COneShotSinglezoneDriver::Run() {

  CQuasiNewtonInvLeastSquares<passivedouble> fixPtCorrector;
  if (config->GetnQuasiNewtonSamples() > 1) {
    fixPtCorrector.resize(config->GetnQuasiNewtonSamples(),
                          geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(),
                          GetTotalNumberOfVariables(ZONE_0,true),
                          geometry_container[ZONE_0][INST_0][MESH_0]->GetnPointDomain());

    if (TimeIter != 0) GetAllSolutions(ZONE_0, true, fixPtCorrector);
  }

  for (auto Adjoint_Iter = 0ul; Adjoint_Iter < nAdjoint_Iter; Adjoint_Iter++) {
    DeformMesh();

    MainRecording();
    /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
     *--- of the previous iteration. The values are passed to the AD tool.
     *--- Issues with iteration number should be dealt with once the output structure is in place. ---*/

    config->SetInnerIter(Adjoint_Iter);

    iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

    /*--- Initialize the adjoint of the objective function with 1.0. ---*/

    SetAdjObjFunction();

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/

    iteration->IterateDiscAdj(geometry_container, solver_container,
                              config_container, ZONE_0, INST_0, false);

    /*--- Monitor the pseudo-time ---*/

    StopCalc = iteration->Monitor(output_container[ZONE_0], integration_container, geometry_container,
                                  solver_container, numerics_container, config_container,
                                  surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();

    SecondaryRecording();

    ComputeDesignVarGradients();

    /*--- Output files for steady state simulations. ---*/

    if (!config->GetTime_Domain()) {
      iteration->Output(output_container[ZONE_0], geometry_container, solver_container,
                        config_container, Adjoint_Iter, false, ZONE_0, INST_0);
      direct_iteration->Output(output_container[ZONE_0], geometry_container, solver_container,
                        config_container, Adjoint_Iter, false, ZONE_0, INST_0);
    }

    UpdateDesignVars();

    if (StopCalc) break;

    /*--- Correct the solution with the quasi-Newton approach. ---*/

    if (fixPtCorrector.size()) {
      GetAllSolutions(ZONE_0, true, fixPtCorrector.FPresult());
      SetAllSolutions(ZONE_0, true, fixPtCorrector.compute());
    }

  }

}

void COneShotSinglezoneDriver::ComputeDesignVarGradients() {

    su2double *VarCoord, *Normal, Area, Sensitivity, my_Gradient, localGradient;

    for (auto iDV = 0u; iDV < config->GetnDV(); iDV++) {
        for (auto iDV_Value = 0u; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++) Gradient[iDV][iDV_Value] = 0;
    }
    
    /*--- Start recording of operations. ---*/

    AD::StartRecording();

    /*--- Register design variables as input and set them to zero
    * (since we want to have the derivative at alpha = 0, i.e. for the current design). ---*/

    for (auto iDV = 0u; iDV < config->GetnDV(); iDV++) {
        for (auto iDV_Value = 0u; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++) {

        AD::RegisterInput(config->GetDV_Value(iDV, iDV_Value));
        }
    }

    /*--- Call the surface deformation routine. ---*/

    surface_movement[ZONE_0]->SetSurface_Deformation(geometry, config);

    /*--- Stop the recording. --- */

    AD::StopRecording();

    /*--- Create a structure to identify points that have been already visited.
   * We need that to make sure to set the sensitivity of surface points only once.
   * Markers share points, so we would visit them more than once in the loop over the markers below). ---*/

  vector<bool> visited(geometry->GetnPoint(), false);

  /*--- Initialize the derivatives of the output of the surface deformation routine
   * with the discrete adjoints from the CFD solution. ---*/

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      auto nVertex = geometry->nVertex[iMarker];
      for (auto iVertex = 0u; iVertex < nVertex; iVertex++) {
        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (!visited[iPoint]) {
          VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          Area = 0.0;
          for (auto iDim = 0u; iDim < nDim; iDim++) {
            Area += Normal[iDim] * Normal[iDim];
          }
          Area = sqrt(Area);

          for (auto iDim = 0u; iDim < nDim; iDim++) {
            if (config->GetDiscrete_Adjoint()) {
              Sensitivity = geometry->GetSensitivity(iPoint, iDim);
            } else {
              Sensitivity = -Normal[iDim] * geometry->vertex[iMarker][iVertex]->GetAuxVar() / Area;
            }
            SU2_TYPE::SetDerivative(VarCoord[iDim], SU2_TYPE::GetValue(Sensitivity));
          }
          visited[iPoint] = true;
        }
      }
    }
  }

  /*--- Compute derivatives and extract gradient. ---*/

  AD::ComputeAdjoint();

  for (auto iDV = 0u; iDV < config->GetnDV(); iDV++) {
    for (auto iDV_Value = 0u; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++) {
      auto my_Gradient = SU2_TYPE::GetDerivative(config->GetDV_Value(iDV, iDV_Value));

      SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

      Gradient[iDV][iDV_Value] += localGradient;
    }
  }

  AD::Reset();

}

void COneShotSinglezoneDriver::DeformMesh() {

    /*--- Set the stiffness of each element mesh into the mesh numerics. ---*/

    solver[MESH_SOL]->SetMesh_Stiffness(numerics[MESH_SOL], config);

    /*--- Deform the volume grid around the new boundary locations. ---*/
    /*--- Force the number of levels to be 0 because in this driver we do not build MG levels. ---*/
    const auto nMGLevels = config->GetnMGLevels();
    config->SetMGLevels(0);
    solver[MESH_SOL]->DeformMesh(geometry_container[ZONE_0][INST_0],numerics[MESH_SOL],config);
    config->SetMGLevels(nMGLevels);
}

void COneShotSinglezoneDriver::UpdateDesignVars() {
    su2double step_size=1e-9;

    for (auto iDV=0u; iDV < config->GetnDV(); iDV++) {
        for (auto iDV_Value = 0u; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++)
            config->SetDV_Value(iDV, iDV_Value, config->GetDV_Value(iDV, iDV_Value) - step_size*Gradient[iDV][iDV_Value]);
    }
}