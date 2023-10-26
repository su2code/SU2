/*!
 * \file CRadP1Solver.cpp
 * \brief Main subroutines for solving generic radiation problems (P1, M1, discrete ordinates...)
 * \author Ruben Sanchez
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

#include "../../include/solvers/CRadSolver.hpp"
#include "../../include/variables/CRadVariable.hpp"

CRadSolver::CRadSolver() : CSolver() {

}

CRadSolver::CRadSolver(CGeometry* geometry, CConfig *config) : CSolver() {

  Absorption_Coeff = config->GetAbsorption_Coeff();
  Scattering_Coeff = config->GetScattering_Coeff();

  Absorption_Coeff = max(Absorption_Coeff,0.01);

}

void CRadSolver::SetVolumetricHeatSource(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;
  unsigned short iDim;

  su2double CP[3]={0.0,0.0,0.0};
  su2double alpha = config->GetHeatSource_Rot_Z() * PI_NUMBER/180.0;
  su2double OP_rot[3]={0.0,0.0,0.0};
  const su2double *OP;
  const su2double *OC = config->GetHeatSource_Center();
  const su2double *Axes = config->GetHeatSource_Axes();
  su2double check;
  // Reset the boolean for all points
  nodes->ResetVol_HeatSource();
  // Loop over all points and determine whether they are inside
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    check = 0;
    OP = geometry->nodes->GetCoord(iPoint);
    // Reference point with respect to center of the ellipse
    for (iDim = 0; iDim < nDim; iDim++) CP[iDim] = OP[iDim]-OC[iDim];
    // Rotate point with respect to Z axis
    OP_rot[0] = OC[0] + CP[0]*cos(alpha) + CP[1]*sin(alpha);
    OP_rot[1] = OC[1] - CP[0]*sin(alpha) + CP[1]*cos(alpha);
    // Check if rotated point is inside the ellipse
    for (iDim = 0; iDim < nDim; iDim++) check += pow(OP_rot[iDim]-OC[iDim],2.0)/pow(Axes[iDim], 2.0);
    if (check <=1) nodes->SetVol_HeatSource(iPoint);
  }

}

void CRadSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/

  unsigned short iVar;
  unsigned long index;

  bool rans = ((config->GetKind_Solver()== MAIN_SOLVER::INC_RANS) ||
               (config->GetKind_Solver()== MAIN_SOLVER::DISC_ADJ_INC_RANS));

  string UnstExt, text_line;
  ifstream restart_file;

  string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);
  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;

  /*--- Skip flow variables ---*/

  unsigned short skipVars = 0;

  if (nDim == 2) skipVars += 6;
  if (nDim == 3) skipVars += 8;

  bool incompressible       = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
  bool energy               = config->GetEnergy_Equation();
  bool weakly_coupled_heat  = config->GetWeakly_Coupled_Heat();

  if (incompressible && ((!energy) && (!weakly_coupled_heat))) skipVars--;

  /*--- Skip turbulent variables ---*/

  if (rans) skipVars += solver[MESH_0][TURB_SOL]->GetnVar();

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {


    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      nodes->SetSolution(iPoint_Local, Solution);
      iPoint_Global_Local++;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) {
    SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- MPI communication ---*/
  solver[MESH_0][RAD_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][RAD_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);

  /*--- Preprocess the fluid solver to compute the primitive variables ---*/
  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);

  /*--- Postprocess the radiation solver to compute the source term that goes into the fluid equations ---*/
  solver[MESH_0][RAD_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);

  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars;
  delete [] Restart_Data;
  Restart_Vars = nullptr; Restart_Data = nullptr;

}
