/*!
 * \file CFluidIteration.cpp
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

#include "../../include/iteration/CFluidIteration.hpp"
#include "../../include/output/COutput.hpp"

void CFluidIteration::Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                 CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                 CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                 CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  unsigned long TimeIter = config[val_iZone]->GetTimeIter();

  bool fsi = config[val_iZone]->GetFSI_Simulation();
  unsigned long OuterIter = config[val_iZone]->GetOuterIter();

  /*--- Set the initial condition for FSI problems with subiterations ---*/
  /*--- This is done only in the first block subiteration.---*/
  /*--- From then on, the solver reuses the partially converged solution obtained in the previous subiteration ---*/
  if (fsi && !config[val_iZone]->GetDiscrete_Adjoint() && (OuterIter == 0)) {
    solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->SetInitialCondition(
        geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst], config[val_iZone], TimeIter);
  }

  /*--- Apply a Wind Gust ---*/

  if (config[val_iZone]->GetWind_Gust()) {
    SetWind_GustField(config[val_iZone], geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst]);
  }
}

void CFluidIteration::Iterate(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                              CSolver***** solver, CNumerics****** numerics, CConfig** config,
                              CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                              CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {

  const bool unsteady = (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                        (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool frozen_visc = (config[val_iZone]->GetContinuous_Adjoint() && config[val_iZone]->GetFrozen_Visc_Cont()) ||
                           (config[val_iZone]->GetDiscrete_Adjoint() && config[val_iZone]->GetFrozen_Visc_Disc());
  const bool disc_adj = (config[val_iZone]->GetDiscrete_Adjoint());

  /*--- Setting up iteration values depending on if this is a
   steady or an unsteady simulation */

  const auto InnerIter = config[val_iZone]->GetInnerIter();
  const auto TimeIter = config[val_iZone]->GetTimeIter();

  /*--- Update global parameters ---*/

  const auto main_solver = config[val_iZone]->GetKind_Solver();
  config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_FLOW_SYS);

  /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/

  integration[val_iZone][val_iInst][FLOW_SOL]->MultiGrid_Iteration(geometry, solver, numerics, config, RUNTIME_FLOW_SYS,
                                                                   val_iZone, val_iInst);

  /*--- If the flow integration is not fully coupled, run the various single grid integrations. ---*/

  if (config[val_iZone]->GetKind_Turb_Model() != TURB_MODEL::NONE && !frozen_visc) {

    /*--- Solve transition model ---*/

    if (config[val_iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_TRANS_SYS);
      integration[val_iZone][val_iInst][TRANS_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                         RUNTIME_TRANS_SYS, val_iZone, val_iInst);
    }

    /*--- Solve the turbulence model ---*/

    config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_TURB_SYS);
    integration[val_iZone][val_iInst][TURB_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                      RUNTIME_TURB_SYS, val_iZone, val_iInst);
  }

  if (config[val_iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
    config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_SPECIES_SYS);
    integration[val_iZone][val_iInst][SPECIES_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                         RUNTIME_SPECIES_SYS, val_iZone, val_iInst);

    // This only applies if mixture properties are used. But this also doesn't hurt if done w/out mixture properties.
    // In case of turbulence, the Turb-Post computes the correct eddy viscosity based on mixture-density and
    // mixture lam-visc. In order to get the correct mixture properties, based on the just updated mass-fractions, the
    // Flow-Pre has to be called upfront. The updated eddy-visc are copied into the flow-solver Primitive in another
    // Flow-Pre call which is done at the start of the next iteration.
    if (config[val_iZone]->GetKind_Turb_Model() != TURB_MODEL::NONE) {
      solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
      solver[val_iZone][val_iInst][MESH_0][TURB_SOL]->Postprocessing(geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0], config[val_iZone], MESH_0);
    }
  }

  if (config[val_iZone]->GetWeakly_Coupled_Heat()) {
    config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_HEAT_SYS);
    integration[val_iZone][val_iInst][HEAT_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                      RUNTIME_HEAT_SYS, val_iZone, val_iInst);
  }

  /*--- Incorporate a weakly-coupled radiation model to the analysis ---*/
  if (config[val_iZone]->AddRadiation()) {
    config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_RADIATION_SYS);
    integration[val_iZone][val_iInst][RAD_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                     RUNTIME_RADIATION_SYS, val_iZone, val_iInst);
  }

  /*--- Adapt the CFL number using an exponential progression with under-relaxation approach. ---*/

  if ((config[val_iZone]->GetCFL_Adapt() == YES) && (!disc_adj)) {
    SU2_OMP_PARALLEL
    solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->AdaptCFLNumber(geometry[val_iZone][val_iInst],
                                                                   solver[val_iZone][val_iInst], config[val_iZone]);
    END_SU2_OMP_PARALLEL
  }

  /*--- Call Dynamic mesh update if AEROELASTIC motion was specified ---*/

  if ((config[val_iZone]->GetGrid_Movement()) && (config[val_iZone]->GetAeroelastic_Simulation()) && unsteady) {
    SetGrid_Movement(geometry[val_iZone][val_iInst], surface_movement[val_iZone], grid_movement[val_iZone][val_iInst],
                     solver[val_iZone][val_iInst], config[val_iZone], InnerIter, TimeIter);

    /*--- Apply a Wind Gust ---*/

    if (config[val_iZone]->GetWind_Gust()) {
      if (InnerIter % config[val_iZone]->GetAeroelasticIter() == 0 && InnerIter != 0)
        SetWind_GustField(config[val_iZone], geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst]);
    }
  }
}

void CFluidIteration::Update(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                             CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                             CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                             unsigned short val_iInst) {
  unsigned short iMesh;

  /*--- Dual time stepping strategy ---*/

  if ((config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
      (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)) {
    /*--- Update dual time solver on all mesh levels ---*/

    for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
      integration[val_iZone][val_iInst][FLOW_SOL]->SetDualTime_Solver(geometry[val_iZone][val_iInst][iMesh],
                                                                      solver[val_iZone][val_iInst][iMesh][FLOW_SOL],
                                                                      config[val_iZone], iMesh);

      integration[val_iZone][val_iInst][FLOW_SOL]->SetDualTime_Geometry(geometry[val_iZone][val_iInst][iMesh],
                                                                        solver[val_iZone][val_iInst][iMesh][MESH_SOL],
                                                                        config[val_iZone], iMesh);
    }

    SetDualTime_Aeroelastic(config[val_iZone]);

    /*--- Update dual time solver for the turbulence model ---*/

    if ((config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::RANS) || (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_RANS) ||
        (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::INC_RANS) ||
        (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_INC_RANS)) {
      integration[val_iZone][val_iInst][TURB_SOL]->SetDualTime_Solver(geometry[val_iZone][val_iInst][MESH_0],
                                                                      solver[val_iZone][val_iInst][MESH_0][TURB_SOL],
                                                                      config[val_iZone], MESH_0);
    }

    /*--- Update dual time solver for the transition model ---*/

    if (config[val_iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      integration[val_iZone][val_iInst][TRANS_SOL]->SetDualTime_Solver(geometry[val_iZone][val_iInst][MESH_0],
                                                                       solver[val_iZone][val_iInst][MESH_0][TRANS_SOL],
                                                                       config[val_iZone], MESH_0);
    }

    /*--- Update dual time solver for the weakly coupled energy equation ---*/

    if (config[val_iZone]->GetWeakly_Coupled_Heat()) {
      integration[val_iZone][val_iInst][HEAT_SOL]->SetDualTime_Solver(geometry[val_iZone][val_iInst][MESH_0],
                                                                      solver[val_iZone][val_iInst][MESH_0][HEAT_SOL],
                                                                      config[val_iZone], MESH_0);
    }
  }
}

bool CFluidIteration::Monitor(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                              CSolver***** solver, CNumerics****** numerics, CConfig** config,
                              CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                              CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  bool StopCalc = false;

  StopTime = SU2_MPI::Wtime();

  UsedTime = StopTime - StartTime;

  output->SetHistoryOutput(geometry[val_iZone][val_iInst][MESH_0], solver[val_iZone][val_iInst][MESH_0],
                           config[val_iZone], config[val_iZone]->GetTimeIter(), config[val_iZone]->GetOuterIter(),
                           config[val_iZone]->GetInnerIter());

  /*--- If convergence was reached --*/
  StopCalc = output->GetConvergence();

  /* --- Checking convergence of Fixed CL mode to target CL, and perform finite differencing if needed  --*/

  if (config[val_iZone]->GetFixed_CL_Mode()) {
    StopCalc = MonitorFixed_CL(output, geometry[val_iZone][INST_0][MESH_0], solver[val_iZone][INST_0][MESH_0],
                               config[val_iZone]);
  }

  return StopCalc;
}

void CFluidIteration::Postprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                  CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                  CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                  CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {

  /*--- Temporary: enable only for single-zone driver. This should be removed eventually when generalized. ---*/
  if (!config[val_iZone]->GetMultizone_Problem()) {

    /*--- Compute the tractions at the vertices ---*/
    solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->ComputeVertexTractions(geometry[val_iZone][val_iInst][MESH_0],
                                                                           config[val_iZone]);
  }
}

void CFluidIteration::Solve(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                            CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                            CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                            unsigned short val_iInst) {
  /*--- Boolean to determine if we are running a static or dynamic case ---*/
  bool steady = !config[val_iZone]->GetTime_Domain();

  unsigned long Inner_Iter, nInner_Iter = config[val_iZone]->GetnInner_Iter();
  bool StopCalc = false;

  /*--- Synchronization point before a single solver iteration.
        Compute the wall clock time required. ---*/

  StartTime = SU2_MPI::Wtime();

  /*--- Preprocess the solver ---*/
  Preprocess(output, integration, geometry, solver, numerics, config, surface_movement, grid_movement, FFDBox,
             val_iZone, INST_0);

  /*--- For steady-state flow simulations, we need to loop over ExtIter for the number of time steps ---*/
  /*--- However, ExtIter is the number of FSI iterations, so nIntIter is used in this case ---*/

  for (Inner_Iter = 0; Inner_Iter < nInner_Iter; Inner_Iter++) {
    config[val_iZone]->SetInnerIter(Inner_Iter);

    /*--- Run a single iteration of the solver ---*/
    Iterate(output, integration, geometry, solver, numerics, config, surface_movement, grid_movement, FFDBox, val_iZone,
            INST_0);

    /*--- Monitor the pseudo-time ---*/
    StopCalc = Monitor(output, integration, geometry, solver, numerics, config, surface_movement, grid_movement, FFDBox,
                       val_iZone, INST_0);

    /*--- Output files at intermediate iterations if the problem is single zone ---*/

    if (singlezone && steady) {
      Output(output, geometry, solver, config, Inner_Iter, StopCalc, val_iZone, val_iInst);
    }

    /*--- If the iteration has converged, break the loop ---*/
    if (StopCalc) break;
  }

  if (multizone && steady) {
    Output(output, geometry, solver, config, config[val_iZone]->GetOuterIter(), StopCalc, val_iZone, val_iInst);
  }
}

void CFluidIteration::SetWind_GustField(CConfig* config, CGeometry** geometry, CSolver*** solver) {
  // The gust is imposed on the flow field via the grid velocities. This method called the Field Velocity Method is
  // described in the NASA TM–2012-217771 - Development, Verification and Use of Gust Modeling in the NASA Computational
  // Fluid Dynamics Code FUN3D the desired gust is prescribed as the negative of the grid velocity.

  unsigned short iDim, nDim = geometry[MESH_0]->GetnDim();

  /*--- Gust Parameters from config ---*/
  unsigned short Gust_Type = config->GetGust_Type();
  su2double xbegin = config->GetGust_Begin_Loc();   // Location at which the gust begins.
  su2double L = config->GetGust_WaveLength();       // Gust size
  su2double tbegin = config->GetGust_Begin_Time();  // Physical time at which the gust begins.
  su2double gust_amp = config->GetGust_Ampl();      // Gust amplitude
  su2double n = config->GetGust_Periods();          // Number of gust periods
  unsigned short GustDir = config->GetGust_Dir();   // Gust direction

  /*--- Variables needed to compute the gust ---*/
  unsigned short Kind_Grid_Movement = config->GetKind_GridMovement();
  unsigned long iPoint;
  unsigned short iMGlevel, nMGlevel = config->GetnMGLevels();

  su2double x, y, x_gust, Gust[3] = {0.0}, NewGridVel[3] = {0.0};
  const su2double* GridVel = nullptr;

  su2double Physical_dt = config->GetDelta_UnstTime();
  unsigned long TimeIter = config->GetTimeIter();
  if (config->GetDiscrete_Adjoint()) TimeIter = config->GetUnst_AdjointIter() - TimeIter - 1;

  su2double Physical_t = TimeIter * Physical_dt;

  su2double Uinf = solver[MESH_0][FLOW_SOL]->GetVelocity_Inf(0);  // Assumption gust moves at infinity velocity

  // Print some information to check that we are doing the right thing. Not sure how to convert the index back to a string...
  if (rank == MASTER_NODE) {
	  cout << endl << "Setting up a wind gust type " << Gust_Type << " with amplitude of " << gust_amp << " in direction " << GustDir << endl;
	  cout << " U_inf      = " << Uinf << endl;
	  cout << " Physical_t = " << Physical_t << endl;
	  su2double loc_x = (xbegin + L + Uinf * (Physical_t - tbegin));
	  cout << " Location_x = " << loc_x << endl;
  }

  // Vortex variables
  unsigned long nVortex = 0;
  vector<su2double> x0, y0, vort_strenth, r_core;  // vortex is positive in clockwise direction.
  if (Gust_Type == VORTEX) {
    InitializeVortexDistribution(nVortex, x0, y0, vort_strenth, r_core);
  }

  /*--- Check to make sure gust lenght is not zero or negative (vortex gust doesn't use this). ---*/
  if (L <= 0.0 && Gust_Type != VORTEX) {
    SU2_MPI::Error("The gust length needs to be positive", CURRENT_FUNCTION);
  }

  /*--- Loop over all multigrid levels ---*/

  for (iMGlevel = 0; iMGlevel <= nMGlevel; iMGlevel++) {
    /*--- Loop over each node in the volume mesh ---*/

    for (iPoint = 0; iPoint < geometry[iMGlevel]->GetnPoint(); iPoint++) {
      /*--- Reset the Grid Velocity to zero if there is no grid movement ---*/
      if (Kind_Grid_Movement == GUST && !(config->GetFSI_Simulation()) && !(config->GetDeform_Mesh())) {
        for (iDim = 0; iDim < nDim; iDim++) geometry[iMGlevel]->nodes->SetGridVel(iPoint, iDim, 0.0);
      }

      /*--- initialize the gust and derivatives to zero everywhere ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        Gust[iDim] = 0.0;
      }

      /*--- Begin applying the gust ---*/

      if (Physical_t >= tbegin) {
        x = geometry[iMGlevel]->nodes->GetCoord(iPoint)[0];  // x-location of the node.
        y = geometry[iMGlevel]->nodes->GetCoord(iPoint)[1];  // y-location of the node.

        // Gust coordinate
        x_gust = (x - xbegin - Uinf * (Physical_t - tbegin)) / L;

        /*--- Calculate the specified gust ---*/
        switch (Gust_Type) {
          case TOP_HAT:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = gust_amp;
              // Still need to put the gust derivatives. Think about this.
            }
            break;

          case SINE:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = gust_amp * (sin(2 * PI_NUMBER * x_gust));
            }
            break;

          case ONE_M_COSINE:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = gust_amp * 0.5 * (1 - cos(2 * PI_NUMBER * x_gust));
            }
            break;

          case EOG:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = -0.37 * gust_amp * sin(3 * PI_NUMBER * x_gust) * (1 - cos(2 * PI_NUMBER * x_gust));
            }
            break;

          case VORTEX:

            /*--- Use vortex distribution ---*/
            // Algebraic vortex equation.
            for (unsigned long i = 0; i < nVortex; i++) {
              su2double r2 = pow(x - (x0[i] + Uinf * (Physical_t - tbegin)), 2) + pow(y - y0[i], 2);
              su2double r = sqrt(r2);
              su2double v_theta = vort_strenth[i] / (2 * PI_NUMBER) * r / (r2 + pow(r_core[i], 2));
              Gust[0] = Gust[0] + v_theta * (y - y0[i]) / r;
              Gust[1] = Gust[1] - v_theta * (x - (x0[i] + Uinf * (Physical_t - tbegin))) / r;
            }
            break;

          case NONE:
          default:

            /*--- There is no wind gust specified. ---*/
            if (rank == MASTER_NODE) {
              cout << "No wind gust specified." << endl;
            }
            break;
        }
      }

      GridVel = geometry[iMGlevel]->nodes->GetGridVel(iPoint);

      /*--- Store new grid velocity ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        NewGridVel[iDim] = GridVel[iDim] - Gust[iDim];
        geometry[iMGlevel]->nodes->SetGridVel(iPoint, iDim, NewGridVel[iDim]);
      }
    }
  }
}

void CFluidIteration::InitializeVortexDistribution(unsigned long& nVortex, vector<su2double>& x0, vector<su2double>& y0,
                                                   vector<su2double>& vort_strength, vector<su2double>& r_core) {
  /*--- Read in Vortex Distribution ---*/
  std::string line;
  std::ifstream file;
  su2double x_temp, y_temp, vort_strength_temp, r_core_temp;
  file.open("vortex_distribution.txt");
  /*--- In case there is no vortex file ---*/
  if (file.fail()) {
    SU2_MPI::Error("There is no vortex data file!!", CURRENT_FUNCTION);
  }

  // Ignore line containing the header
  getline(file, line);
  // Read in the information of the vortices (xloc, yloc, lambda(strength), eta(size, gradient))
  while (file.good()) {
    getline(file, line);
    std::stringstream ss(line);
    if (!line.empty()) {  // ignore blank lines if they exist.
      ss >> x_temp;
      ss >> y_temp;
      ss >> vort_strength_temp;
      ss >> r_core_temp;
      x0.push_back(x_temp);
      y0.push_back(y_temp);
      vort_strength.push_back(vort_strength_temp);
      r_core.push_back(r_core_temp);
    }
  }
  file.close();
  // number of vortices
  nVortex = x0.size();
}

bool CFluidIteration::MonitorFixed_CL(COutput *output, CGeometry *geometry, CSolver **solver, CConfig *config) {

  CSolver* flow_solver= solver[FLOW_SOL];

  bool fixed_cl_convergence = flow_solver->FixedCL_Convergence(config, output->GetConvergence());

  /* --- If Fixed CL mode has ended and Finite Differencing has started: --- */

  if (flow_solver->GetStart_AoA_FD() && flow_solver->GetIter_Update_AoA() == config->GetInnerIter()){

    /* --- Print convergence history and volume files since fixed CL mode has converged--- */
    if (rank == MASTER_NODE) output->PrintConvergenceSummary();

    output->SetResultFiles(geometry, config, solver,
                            config->GetInnerIter(), true);

    /* --- Set finite difference mode in config (disables output) --- */
    config->SetFinite_Difference_Mode(true);
  }

  /* --- Set convergence based on fixed CL convergence  --- */
  return fixed_cl_convergence;
}

void CFluidIteration::SetDualTime_Aeroelastic(CConfig* config) const {

  /*--- Store old aeroelastic solutions ---*/

  if (config->GetGrid_Movement() && config->GetAeroelastic_Simulation()) {

    config->SetAeroelastic_n1();
    config->SetAeroelastic_n();

    /*--- Also communicate plunge and pitch to the master node. Needed for output in case of parallel run ---*/

#ifdef HAVE_MPI
    su2double plunge, pitch, *plunge_all = nullptr, *pitch_all = nullptr;
    unsigned short iMarker, iMarker_Monitoring;
    unsigned long iProcessor, owner, *owner_all = nullptr;

    string Marker_Tag, Monitoring_Tag;
    int nProcessor = size;

    /*--- Only if master node allocate memory ---*/

    if (rank == MASTER_NODE) {
      plunge_all = new su2double[nProcessor];
      pitch_all  = new su2double[nProcessor];
      owner_all  = new unsigned long[nProcessor];
    }

    /*--- Find marker and give it's plunge and pitch coordinate to the master node ---*/

    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {

      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag) { owner = 1; break;
        }           owner = 0;
       

      }
      plunge = config->GetAeroelastic_plunge(iMarker_Monitoring);
      pitch  = config->GetAeroelastic_pitch(iMarker_Monitoring);

      /*--- Gather the data on the master node. ---*/

      SU2_MPI::Gather(&plunge, 1, MPI_DOUBLE, plunge_all, 1, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
      SU2_MPI::Gather(&pitch, 1, MPI_DOUBLE, pitch_all, 1, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
      SU2_MPI::Gather(&owner, 1, MPI_UNSIGNED_LONG, owner_all, 1, MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());

      /*--- Set plunge and pitch on the master node ---*/

      if (rank == MASTER_NODE) {
        for (iProcessor = 0; iProcessor < (unsigned long)nProcessor; iProcessor++) {
          if (owner_all[iProcessor] == 1) {
            config->SetAeroelastic_plunge(iMarker_Monitoring, plunge_all[iProcessor]);
            config->SetAeroelastic_pitch(iMarker_Monitoring, pitch_all[iProcessor]);
            break;
          }
        }
      }
    }

    delete [] plunge_all;
    delete [] pitch_all;
    delete [] owner_all;
#endif
  }

}
