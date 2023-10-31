/*!
 * \file CHeatSolver.cpp
 * \brief Main subroutines for solving the heat equation
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

#include "../../include/solvers/CHeatSolver.hpp"
#include <cstddef>
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/solvers/CScalarSolver.inl"

/*--- Explicit instantiation of the parent class of CHeatSolver. ---*/
template class CScalarSolver<CHeatVariable>;

CHeatSolver::CHeatSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
  : CScalarSolver<CHeatVariable>(geometry, config, false),
    flow(config->GetFluidProblem()), heat_equation(config->GetHeatProblem()) {

  /*--- Dimension of the problem --> temperature is the only conservative variable ---*/

  nVar = 1;
  nPrimVar = 1;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Define some structures for locating max residuals ---*/

  Residual_RMS.resize(nVar,0.0);
  Residual_Max.resize(nVar,0.0);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

  /*--- Initialization of the structure of the whole Jacobian ---*/

  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (heat equation) MG level: " << iMesh << "." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  if (ReducerStrategy) EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);

  if (config->GetExtraOutput()) {
    if (nDim == 2) { nOutputVariables = 13; }
    else if (nDim == 3) { nOutputVariables = 19; }
    OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
    OutputHeadingNames = new string[nOutputVariables];
  }

  HeatFlux_per_Marker.resize(nMarker, 0.0);
  AverageT_per_Marker.resize(nMarker, 0.0);
  Surface_Areas.resize(config->GetnMarker_HeatFlux(), 0.0);

  Set_Heatflux_Areas(geometry, config);

  /*--- Set the reference values for temperature ---*/

  su2double Temperature_FreeStream = config->GetTemperature_FreeStream();
  su2double Temperature_Ref = 0.0;

  if (config->GetRef_Inc_NonDim() == DIMENSIONAL) {
    Temperature_Ref = 1.0;
  }
  else if (config->GetRef_Inc_NonDim() == INITIAL_VALUES) {
    Temperature_Ref = Temperature_FreeStream;
  }
  else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
    Temperature_Ref = config->GetInc_Temperature_Ref();
  }
  config->SetTemperature_Ref(Temperature_Ref);
  config->SetTemperature_FreeStreamND(Temperature_FreeStream / Temperature_Ref);
  Solution_Inf[0] = config->GetTemperature_FreeStreamND();

  /*--- Set the reference values for heat fluxes. If the heat solver runs stand-alone,
   *    thermal conductivity is read directly from config file ---*/

  if (heat_equation) {
    su2double rho_cp = config->GetMaterialDensity(0)*config->GetSpecific_Heat_Cp();
    config->SetThermalDiffusivity(config->GetThermal_Conductivity_Constant() / rho_cp);

    /*--- Fluxes are computed via thermal diffusivity (not conductivity), so we have to divide by rho*cp ---*/
    config->SetHeat_Flux_Ref(rho_cp*Temperature_Ref);
  }
  else if (flow) {
    config->SetHeat_Flux_Ref(config->GetViscosity_Ref()*config->GetSpecific_Heat_Cp());
  }

  config->SetDelta_UnstTimeND(config->GetDelta_UnstTime() / config->GetTime_Ref());

  /*--- Store the value of the temperature and the heat flux density at the boundaries,
   used for communications with donor cells ---*/

  const unsigned short nConjVariables = 4;
  AllocVectorOfMatrices(nVertex, nConjVariables, ConjugateVar, Solution_Inf[0]);

  /*--- Heat flux in all the markers ---*/

  AllocVectorOfVectors(nVertex, HeatFlux);

  if (config->GetMultizone_Problem()){
    /*--- Initialize the BGS residuals. ---*/
    Residual_BGS.resize(nVar,1.0);
    Residual_Max_BGS.resize(nVar,1.0);
    Point_Max_BGS.resize(nVar,0);
    Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
  }

  /*--- Initialize the nodes vector. ---*/

  nodes = new CHeatVariable(Solution_Inf[0], nPoint, nDim, nVar, config);

  SetBaseClassPointerToNodes();

  /*--- Communicate and store volume and the number of neighbors for any dual CVs that lie on on periodic markers. ---*/
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic() / 2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_VOLUME);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_VOLUME);
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_NEIGHBORS);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_NEIGHBORS);
  }
  /*--- Store if implicit scheme is used. This has implications on the Residual and Jacobian handling for periodic
   * boundaries  ---*/
  const bool euler_implicit = (config->GetKind_TimeIntScheme_Heat() == EULER_IMPLICIT);
  SetImplicitPeriodic(euler_implicit);

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Store the initial CFL number for all grid points. ---*/

  const su2double CFL = config->GetCFL(MGLevel) * config->GetCFLRedCoeff_Turb();
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  END_SU2_OMP_FOR
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name. ---*/

  SolverName = "HEAT";
}

void CHeatSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)
  CommonPreprocessing(geometry, config, Output);

  /*--- Need to clear EdgeFluxes and Jacobian when only the viscous part is called for solid heat transfer,
   * for the weakly coupled energy equation the convection part does this by setting instead of incrementing. ---*/
  if (!Output && !flow && ReducerStrategy) {
    EdgeFluxes.SetValZero();
    if (config->GetKind_TimeIntScheme() == EULER_IMPLICIT) {
      Jacobian.SetValZero();
    } else {
      SU2_OMP_BARRIER
    }
  }
}

void CHeatSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter,
                              bool val_update_geo) {
  const string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
  /*--- Skip coordinates ---*/

  unsigned short skipVars = nDim;

  if (flow) {
    // P, vx, vy (,vz)
    skipVars += 1 + nDim;
  }

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/
  unsigned long counter = 0;
  for (auto iPoint_Global = 0ul; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    const auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {
      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      const auto index = counter*Restart_Vars[1] + skipVars;
      for (auto iVar = 0u; iVar < nVar; iVar++) nodes->SetSolution(iPoint_Local, iVar, Restart_Data[index + iVar]);

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }
  }

  /*--- Detect a wrong solution file ---*/

  if (counter != nPointDomain) {
    SU2_MPI::Error(string("The solution file ") + restart_filename + string(" does not match with the mesh file!\n") +
                   string("This can be caused by empty lines at the end of the file."), CURRENT_FUNCTION);
  }
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. We alo call the preprocessing routine
   on the fine level in order to have all necessary quantities updated,
   especially if this is a turbulent simulation (eddy viscosity). ---*/

  solver[MESH_0][HEAT_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][HEAT_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);

  solver[MESH_0][HEAT_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_HEAT_SYS, false);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (auto iMesh = 1u; iMesh <= config->GetnMGLevels(); iMesh++) {
    MultigridRestriction(*geometry[iMesh - 1], solver[iMesh - 1][HEAT_SOL]->GetNodes()->GetSolution(),
                         *geometry[iMesh], solver[iMesh][HEAT_SOL]->GetNodes()->GetSolution());
    solver[iMesh][HEAT_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][HEAT_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);

    solver[iMesh][HEAT_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_HEAT_SYS, false);
  }

  /*--- Delete the class memory that is used to load the restart. ---*/

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    delete[] Restart_Vars;
    Restart_Vars = nullptr;
    delete[] Restart_Data;
    Restart_Data = nullptr;
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

void CHeatSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                                  CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {
  /*--- For solid heat transfer there is no convection. ---*/
  if (!flow) return;
  CScalarSolver<CHeatVariable>::Upwind_Residual(geometry, solver_container, numerics_container, config, iMesh);
}

void CHeatSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  /*--- For fluid problems the viscous residual is included in the convective residual. ---*/
  if (flow) return;
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  CNumerics* numerics = numerics_container[VISC_TERM + omp_get_thread_num() * MAX_TERMS];

  bool pausePreacc = false;
  if (ReducerStrategy)
    pausePreacc = AD::PausePreaccumulation();
  else
    AD::StartNoSharedReading();

  for (auto color : EdgeColoring) {
    SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
    for (auto k = 0ul; k < color.size; ++k) {
      auto iEdge = color.indices[k];
      Viscous_Residual(iEdge, geometry, solver_container, numerics, config);
    }
    END_SU2_OMP_FOR
  }

  /*--- Restore preaccumulation and adjoint evaluation state. ---*/
  AD::ResumePreaccumulation(pausePreacc);
  if (!ReducerStrategy) AD::EndNoSharedReading();

  if (ReducerStrategy) {
    SumEdgeFluxes(geometry);
    if (implicit) Jacobian.SetDiagonalAsColumnSum();
  }
}

void CHeatSolver::Set_Heatflux_Areas(CGeometry *geometry, CConfig *config) {

  string HeatFlux_Tag, Marker_Tag;

  su2double *Local_Surface_Areas, Local_HeatFlux_Areas_Monitor, Area, *Normal;
  Local_Surface_Areas = new su2double[config->GetnMarker_HeatFlux()];

  for (auto  iMarker_HeatFlux = 0u; iMarker_HeatFlux < config->GetnMarker_HeatFlux(); iMarker_HeatFlux++ ) {
    Local_Surface_Areas[iMarker_HeatFlux] = 0.0;
  }
  Local_HeatFlux_Areas_Monitor = 0.0;

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {

    const auto Monitoring = config->GetMarker_All_Monitoring(iMarker);

    for (auto iMarker_HeatFlux = 0u; iMarker_HeatFlux < config->GetnMarker_HeatFlux(); iMarker_HeatFlux++ ) {

      HeatFlux_Tag = config->GetMarker_HeatFlux_TagBound(iMarker_HeatFlux);
      Marker_Tag = config->GetMarker_All_TagBound(iMarker);

      if (Marker_Tag == HeatFlux_Tag) {

        Local_Surface_Areas[iMarker_HeatFlux] = 0.0;

        for(auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++ ) {

          const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if(!geometry->nodes->GetDomain(iPoint)) continue;

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = GeometryToolbox::Norm(nDim, Normal);

          Local_Surface_Areas[iMarker_HeatFlux] += Area;

          if(Monitoring == YES)
            Local_HeatFlux_Areas_Monitor += Area;

        }
      }
    }
  }

  SU2_MPI::Allreduce(Local_Surface_Areas, Surface_Areas.data(), Surface_Areas.size(), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&Local_HeatFlux_Areas_Monitor, &Total_HeatFlux_Areas_Monitor, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

  Total_HeatFlux_Areas = 0.0;
  for(auto iMarker_HeatFlux = 0u; iMarker_HeatFlux < config->GetnMarker_HeatFlux(); iMarker_HeatFlux++ ) {
    Total_HeatFlux_Areas += Surface_Areas[iMarker_HeatFlux];
  }

  delete[] Local_Surface_Areas;
}

void CHeatSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  su2double Twall = config->GetIsothermal_Temperature(Marker_Tag) / config->GetTemperature_Ref();
  const bool IsPyCustom = config->GetMarker_All_PyCustom(val_marker);

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    if (IsPyCustom) {
      Twall = geometry->GetCustomBoundaryTemperature(val_marker, iVertex) / config->GetTemperature_Ref();
    }
    IsothermalBoundaryCondition(geometry, solver_container[FLOW_SOL], config, val_marker, iVertex, Twall);
  }
  END_SU2_OMP_FOR
}

void CHeatSolver::BC_HeatFlux_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                   CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  const bool IsPyCustom = config->GetMarker_All_PyCustom(val_marker);

  su2double Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag) / config->GetHeat_Flux_Ref();
  if (config->GetIntegrated_HeatFlux()) {
    Wall_HeatFlux /= geometry->GetSurfaceArea(config, val_marker);
  }

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    if (!geometry->nodes->GetDomain(iPoint)) continue;

    if (IsPyCustom) {
      Wall_HeatFlux = geometry->GetCustomBoundaryHeatFlux(val_marker, iVertex) / config->GetHeat_Flux_Ref();
    }
    /*--- Viscous contribution to the residual at the wall. ---*/
    const auto* Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
    const su2double Area = GeometryToolbox::Norm(nDim, Normal);
    const su2double flux = Wall_HeatFlux * Area;
    LinSysRes(iPoint, 0) -= flux;
  }
  END_SU2_OMP_FOR
}

void CHeatSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                           CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  if (!flow) return;

  const bool viscous = config->GetViscous();
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    su2double Normal[MAXNDIM];
    geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
    for (auto iDim = 0u; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

    /*--- Normal vector for this vertex (negate for outward convention) ---*/

    conv_numerics->SetNormal(Normal);

    /*--- Retrieve solution at this boundary node ---*/

    const auto* V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

    /*--- Retrieve the specified velocity for the inlet. ---*/

    const auto* V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

    conv_numerics->SetPrimitive(V_domain, V_inlet);

    if (dynamic_grid)
      conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

    const su2double Temp_i = nodes->GetTemperature(iPoint);
    const su2double Temp_j = V_inlet[prim_idx.Temperature()];
    conv_numerics->SetScalarVar(&Temp_i, &Temp_j);

    /*--- Compute the residual using an upwind scheme ---*/

    auto residual = conv_numerics->ComputeResidual(config);

    /*--- Update residual value ---*/

    LinSysRes.AddBlock(iPoint, residual);

    /*--- Jacobian contribution for implicit integration ---*/

    if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    /*--- Viscous contribution ---*/

    if (viscous) {
      IsothermalBoundaryCondition(geometry, solver_container[FLOW_SOL], config, val_marker, iVertex, Temp_j);
    }
  }
  END_SU2_OMP_FOR
}

void CHeatSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  if (!flow) return;

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    const auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

    /*--- Normal vector for this vertex (negate for outward convention) ---*/

    su2double Normal[MAXNDIM];
    geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
    for (auto iDim = 0u; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

    conv_numerics->SetNormal(Normal);

    /*--- Retrieve solution at this boundary node ---*/

    const auto* V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

    /*--- Retrieve the specified velocity for the inlet. ---*/

    const auto* V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

    conv_numerics->SetPrimitive(V_domain, V_outlet);

    if (dynamic_grid)
      conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

    const su2double Temp_i = nodes->GetTemperature(iPoint);
    const su2double Temp_j = nodes->GetTemperature(Point_Normal);
    conv_numerics->SetScalarVar(&Temp_i, &Temp_j);

    /*--- Compute the residual using an upwind scheme ---*/

    auto residual = conv_numerics->ComputeResidual(config);

    /*--- Update residual value ---*/

    LinSysRes.AddBlock(iPoint, residual);

    /*--- Jacobian contribution for implicit integration ---*/

    if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
  }
  END_SU2_OMP_FOR
}

void CHeatSolver::BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const su2double Temperature_Ref = config->GetTemperature_Ref();
  const su2double rho_cp_solid = config->GetMaterialDensity(0) * config->GetSpecific_Heat_Cp();

  if (flow) {
    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

      const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

      if (geometry->nodes->GetDomain(iPoint)) {
        const su2double T_Conjugate = GetConjugateHeatVariable(val_marker, iVertex, 0) / Temperature_Ref;

        nodes->SetSolution_Old(iPoint, &T_Conjugate);
        LinSysRes(iPoint, 0) = 0.0;
        nodes->SetRes_TruncErrorZero(iPoint);

        if (implicit) Jacobian.DeleteValsRowi(iPoint);
      }
    }
    END_SU2_OMP_FOR
  }
  else if (heat_equation) {
    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

      const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

      if (geometry->nodes->GetDomain(iPoint)) {

        su2double const* Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
        const su2double Area = GeometryToolbox::Norm(nDim, Normal);

        const su2double thermal_diffusivity = GetConjugateHeatVariable(val_marker, iVertex, 2) / rho_cp_solid;
        su2double HeatFlux = 0;

        if ((config->GetKind_CHT_Coupling() == CHT_COUPLING::DIRECT_TEMPERATURE_ROBIN_HEATFLUX) ||
            (config->GetKind_CHT_Coupling() == CHT_COUPLING::AVERAGED_TEMPERATURE_ROBIN_HEATFLUX)) {

          const su2double Tinterface = nodes->GetTemperature(iPoint);
          const su2double Tnormal_Conjugate = GetConjugateHeatVariable(val_marker, iVertex, 3) / Temperature_Ref;

          const su2double HeatFluxDensity = thermal_diffusivity * (Tinterface - Tnormal_Conjugate);
          HeatFlux = HeatFluxDensity * Area;

          if (implicit) {
            su2double Jacobian_i[] = {-thermal_diffusivity*Area};
            Jacobian.SubtractBlock2Diag(iPoint, &Jacobian_i);
          }
        }
        else {
          const su2double HeatFluxDensity =
              GetConjugateHeatVariable(val_marker, iVertex, 1) / config->GetHeat_Flux_Ref();
          HeatFlux = HeatFluxDensity*Area;
        }

        LinSysRes(iPoint, 0) += HeatFlux;
      }
    }
    END_SU2_OMP_FOR
  }
}

void CHeatSolver::Heat_Fluxes(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned long iPointNormal;
  unsigned short Boundary, Monitoring;
  su2double *Coord, *Coord_Normal, *Normal, Area, dist, Twall, dTdn;
  string Marker_Tag, HeatFlux_Tag;

  const su2double thermal_diffusivity = flow ? config->GetViscosity_FreeStreamND()/config->GetPrandtl_Lam() :
                                               config->GetThermalDiffusivity();

  AllBound_HeatFlux = 0.0;
  AllBound_AverageT = 0.0;

  for (auto iMarker = 0u; iMarker < nMarker; iMarker++ ) {

    AverageT_per_Marker[iMarker]  = 0.0;
    HeatFlux_per_Marker[iMarker]  = 0.0;

    Boundary = config->GetMarker_All_KindBC(iMarker);
    Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);

    if ( Boundary == ISOTHERMAL ) {

      Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

      for(auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++ ) {

        const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if(geometry->nodes->GetDomain(iPoint)) {

          iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

          Coord = geometry->nodes->GetCoord(iPoint);
          Coord_Normal = geometry->nodes->GetCoord(iPointNormal);

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = GeometryToolbox::Norm(nDim, Normal);

          dist = 0.0;
          for (auto iDim = 0u; iDim < nDim; iDim++) dist += (Coord_Normal[iDim]-Coord[iDim])*(Coord_Normal[iDim]-Coord[iDim]);
          dist = sqrt(dist);

          dTdn = (Twall - nodes->GetTemperature(iPointNormal))/dist;

          HeatFlux[iMarker][iVertex] = thermal_diffusivity*dTdn*config->GetHeat_Flux_Ref();

          HeatFlux_per_Marker[iMarker] += HeatFlux[iMarker][iVertex]*Area;

        }
      }
    }
    else if ( Boundary == CHT_WALL_INTERFACE || Boundary == HEAT_FLUX ) {

      for(auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++ ) {

        const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if(geometry->nodes->GetDomain(iPoint)) {

          iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

          Twall = nodes->GetTemperature(iPoint);

          Coord = geometry->nodes->GetCoord(iPoint);
          Coord_Normal = geometry->nodes->GetCoord(iPointNormal);

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = GeometryToolbox::Norm(nDim, Normal);

          dist = 0.0;
          for (auto iDim = 0u; iDim < nDim; iDim++) dist += (Coord_Normal[iDim]-Coord[iDim])*(Coord_Normal[iDim]-Coord[iDim]);
          dist = sqrt(dist);

          dTdn = (Twall - nodes->GetTemperature(iPointNormal))/dist;

          HeatFlux[iMarker][iVertex] = thermal_diffusivity*dTdn*config->GetHeat_Flux_Ref();

          HeatFlux_per_Marker[iMarker] += HeatFlux[iMarker][iVertex]*Area;

          /*--- We do only aim to compute averaged temperatures on the (interesting) heat flux walls ---*/

          if ( Boundary == HEAT_FLUX ) {

            AverageT_per_Marker[iMarker] += Twall*config->GetTemperature_Ref()*Area;
          }
        }
      }
    }

    if (Monitoring == YES) {

      AllBound_HeatFlux += HeatFlux_per_Marker[iMarker];
      AllBound_AverageT += AverageT_per_Marker[iMarker];
    }
  }

  su2double Send_Bound_HeatFlux = AllBound_HeatFlux;
  su2double Send_Bound_AverageT = AllBound_AverageT;
  SU2_MPI::Allreduce(&Send_Bound_HeatFlux, &AllBound_HeatFlux, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&Send_Bound_AverageT, &AllBound_AverageT, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

  if (Total_HeatFlux_Areas_Monitor != 0.0) {
    Total_AverageT = AllBound_AverageT/Total_HeatFlux_Areas_Monitor;
  }
  else {
    Total_AverageT = 0.0;
  }

  Total_HeatFlux = AllBound_HeatFlux;
}

void CHeatSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                               unsigned short iMesh, unsigned long Iteration) {

  const bool implicit      = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool time_stepping = (config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING);
  const bool dual_time     = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                             (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const su2double K_v = 0.25;

  const su2double prandtl_lam = config->GetPrandtl_Lam();
  const su2double prandtl_turb = config->GetPrandtl_Turb();
  const su2double constant_thermal_diffusivity = config->GetThermalDiffusivity();

  const CVariable* flow_nodes = solver_container[FLOW_SOL] ? solver_container[FLOW_SOL]->GetNodes() : nullptr;

  /*--- Init thread-shared variables to compute min/max values.
   *    Critical sections are used for this instead of reduction
   *    clauses for compatibility with OpenMP 2.0 (Windows...). ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    Min_Delta_Time = 1e30;
    Max_Delta_Time = 0.0;
    Global_Delta_UnstTimeND = 1e30;
  } END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- Loop domain points. ---*/

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {

    /*--- Set maximum eigenvalues to zero. ---*/

    nodes->SetMax_Lambda_Visc(iPoint, 0.0);

    /*--- Loop over the neighbors of point i. ---*/

    for (unsigned short iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); ++iNeigh) {
      const auto jPoint = geometry->nodes->GetPoint(iPoint,iNeigh);
      const auto iEdge = geometry->nodes->GetEdge(iPoint,iNeigh);
      const auto* Normal = geometry->edges->GetNormal(iEdge);
      const su2double Area2 = GeometryToolbox::SquaredNorm(nDim, Normal);

      if (flow) {
        const su2double laminar_viscosity = 0.5 * (flow_nodes->GetLaminarViscosity(iPoint) +
                                                   flow_nodes->GetLaminarViscosity(jPoint));
        const su2double eddy_viscosity = 0.5 * (flow_nodes->GetEddyViscosity(iPoint) +
                                                flow_nodes->GetEddyViscosity(jPoint));

        const su2double thermal_diffusivity = laminar_viscosity / prandtl_lam + eddy_viscosity / prandtl_turb;
        nodes->AddMax_Lambda_Visc(iPoint, thermal_diffusivity * Area2);
      } else {
        nodes->AddMax_Lambda_Visc(iPoint, constant_thermal_diffusivity * Area2);
      }
    }
  }
  END_SU2_OMP_FOR

  /*--- Loop boundary edges ---*/

  for (unsigned short iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {

      SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
      for (auto iVertex = 0ul; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        /*--- Point identification, Normal vector and area ---*/

        const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (!geometry->nodes->GetDomain(iPoint)) continue;

        const auto* Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        const su2double Area2 = GeometryToolbox::SquaredNorm(nDim, Normal);

        if (flow) {
          const su2double thermal_diffusivity = flow_nodes->GetLaminarViscosity(iPoint) / prandtl_lam +
                                                flow_nodes->GetEddyViscosity(iPoint) / prandtl_turb;
          nodes->AddMax_Lambda_Visc(iPoint, thermal_diffusivity * Area2);
        } else {
          nodes->AddMax_Lambda_Visc(iPoint, constant_thermal_diffusivity * Area2);
        }
      }
      END_SU2_OMP_FOR
    }
  }

  /*--- Each element uses their own speed, steady state simulation. ---*/
  {
    /*--- Thread-local variables for min/max reduction. ---*/
    su2double minDt = 1e30, maxDt = 0.0;

    SU2_OMP_FOR_(schedule(static,omp_chunk_size) SU2_NOWAIT)
    for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {

      su2double local_delta_time = 0.0;

      const su2double Vol = geometry->nodes->GetVolume(iPoint);
      if (Vol != 0.0) {
        local_delta_time = nodes->GetLocalCFL(iPoint) * K_v * pow(Vol, 2) / nodes->GetMax_Lambda_Visc(iPoint);
        if (flow) {
          switch (config->GetKind_TimeStep_Heat()) {
            case BYFLOW:
            case CONVECTIVE:
              local_delta_time = flow_nodes->GetDelta_Time(iPoint);
              break;
            case MINIMUM:
              local_delta_time = min(local_delta_time, flow_nodes->GetDelta_Time(iPoint));
              break;
            case VISCOUS:
              break;
            default:
              break;
          }
        }
        minDt = min(minDt, local_delta_time);
        maxDt = max(maxDt, local_delta_time);
      }
      nodes->SetDelta_Time(iPoint, local_delta_time);
    }
    END_SU2_OMP_FOR
    /*--- Min/max over threads. ---*/
    SU2_OMP_CRITICAL
    {
      Min_Delta_Time = min(Min_Delta_Time, minDt);
      Max_Delta_Time = max(Max_Delta_Time, maxDt);
      Global_Delta_Time = Min_Delta_Time;
    }
    END_SU2_OMP_CRITICAL
  }

  /*--- Compute the min/max dt (in parallel, now over mpi ranks). ---*/

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    if (config->GetComm_Level() == COMM_FULL) {
      su2double rbuf_time;
      SU2_MPI::Allreduce(&Min_Delta_Time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      Min_Delta_Time = rbuf_time;

      SU2_MPI::Allreduce(&Max_Delta_Time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
      Max_Delta_Time = rbuf_time;
    }
  } END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- For exact time solution use the minimum delta time of the whole mesh. ---*/
  if (time_stepping) {

    /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step,
     *    otherwise it computes the time step based on the unsteady CFL. ---*/

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    {
      if (config->GetUnst_CFL() == 0.0) {
        Global_Delta_Time = config->GetDelta_UnstTime();
      }
      else {
        Global_Delta_Time = Min_Delta_Time;
      }
      Max_Delta_Time = Global_Delta_Time;

      config->SetDelta_UnstTimeND(Global_Delta_Time);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    /*--- Sets the regular CFL equal to the unsteady CFL. ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {
      nodes->SetLocalCFL(iPoint, config->GetUnst_CFL());
      nodes->SetDelta_Time(iPoint, Global_Delta_Time);
    }
    END_SU2_OMP_FOR

  }

  /*--- Recompute the unsteady time step for the dual time strategy if the unsteady CFL is diferent from 0.
    * This is only done once because in dual time the time step cannot be variable. ---*/

  if (dual_time && (Iteration == config->GetRestart_Iter()) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {

    /*--- Thread-local variable for reduction. ---*/
    su2double glbDtND = 1e30;

    SU2_OMP_FOR_(schedule(static,omp_chunk_size) SU2_NOWAIT)
    for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {
      glbDtND = min(glbDtND, config->GetUnst_CFL()*Global_Delta_Time / nodes->GetLocalCFL(iPoint));
    }
    END_SU2_OMP_FOR
    SU2_OMP_CRITICAL
    Global_Delta_UnstTimeND = min(Global_Delta_UnstTimeND, glbDtND);
    END_SU2_OMP_CRITICAL

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    {
      SU2_MPI::Allreduce(&Global_Delta_UnstTimeND, &glbDtND, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      Global_Delta_UnstTimeND = glbDtND;

      config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/

  if (dual_time && !implicit) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {
      su2double dt = min((2.0/3.0)*config->GetDelta_UnstTimeND(), nodes->GetDelta_Time(iPoint));
      nodes->SetDelta_Time(iPoint, dt);
    }
    END_SU2_OMP_FOR
  }
}

void CHeatSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  const bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  const bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                          (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));

  /*--- If restart solution, then interpolate the flow solution to
   all the multigrid levels, this is important with the dual time strategy ---*/

  if (restart && (TimeIter == 0)) {

    for (auto iMesh = 1u; iMesh <= config->GetnMGLevels(); iMesh++) {
      MultigridRestriction(*geometry[iMesh - 1], solver_container[iMesh - 1][HEAT_SOL]->GetNodes()->GetSolution(),
                           *geometry[iMesh], solver_container[iMesh][HEAT_SOL]->GetNodes()->GetSolution());
      solver_container[iMesh][HEAT_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
      solver_container[iMesh][HEAT_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    }
  }

  /*--- The value of the solution for the first iteration of the dual time ---*/

  if (dual_time && (TimeIter == 0 || (restart && (long)TimeIter == (long)config->GetRestart_Iter()))) {

    /*--- Push back the initial condition to previous solution containers
     for a 1st-order restart or when simply intitializing to freestream. ---*/

    for (auto iMesh = 0u; iMesh <= config->GetnMGLevels(); iMesh++) {
      solver_container[iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n();
      solver_container[iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n1();
    }

    if ((restart && (long)TimeIter == (long)config->GetRestart_Iter()) &&
        (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)) {

      /*--- Load an additional restart file for a 2nd-order restart ---*/

      solver_container[MESH_0][HEAT_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetRestart_Iter()-1), true);

      /*--- Push back this new solution to time level N. ---*/

      for (auto iMesh = 0u; iMesh <= config->GetnMGLevels(); iMesh++) {
        solver_container[iMesh][HEAT_SOL]->GetNodes()->Set_Solution_time_n();
      }
    }
  }
}

void CHeatSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                       unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
  if (flow) {
    CScalarSolver<CHeatVariable>::SetResidual_DualTime(geometry, solver_container, config, iRKStep, iMesh, RunTime_EqSystem);
    return;
  }
  const bool first_order  = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST);
  const bool second_order = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Store the physical time step ---*/

  const su2double TimeStep = config->GetDelta_UnstTimeND();

  /*--- Compute the dual time-stepping source term ---*/
  /*--- Loop over all nodes (excluding halos) ---*/

  AD::StartNoSharedReading();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {

    /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
      we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
      previous solutions that are stored in memory. ---*/

    const auto* U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
    const auto* U_time_n = nodes->GetSolution_time_n(iPoint);
    const auto* U_time_nP1 = nodes->GetSolution(iPoint);

    /*--- CV volume at time n+1. As we are on a static mesh, the volume
      of the CV will remained fixed for all time steps. ---*/

    const su2double Volume_nP1 = geometry->nodes->GetVolume(iPoint);

    /*--- Compute the dual time-stepping source term based on the chosen
      time discretization scheme (1st- or 2nd-order).---*/

    for (auto iVar = 0u; iVar < nVar; iVar++) {
      if (first_order)
        LinSysRes(iPoint, iVar) += (U_time_nP1[iVar] - U_time_n[iVar]) * Volume_nP1 / TimeStep;
      if (second_order)
        LinSysRes(iPoint, iVar) += (3.0 * U_time_nP1[iVar] - 4.0 * U_time_n[iVar] +
                                    1.0 * U_time_nM1[iVar]) * Volume_nP1 / (2.0 * TimeStep);
    }

    /*--- Compute the Jacobian contribution due to the dual time source term. ---*/
    if (implicit) {
      if (first_order) Jacobian.AddVal2Diag(iPoint, Volume_nP1 / TimeStep);
      if (second_order) Jacobian.AddVal2Diag(iPoint, (Volume_nP1 * 3.0) / (2.0 * TimeStep));
    }
  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();
}
