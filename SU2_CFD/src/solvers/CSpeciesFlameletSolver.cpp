/*!
 * \file CSpeciesFlameletSolver.cpp
 * \brief Main subroutines of CSpeciesFlameletSolver class
 * \author D. Mayer, T. Economon, N. Beishuizen, E. Bunschoten
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/solvers/CSpeciesFlameletSolver.hpp"

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/fluid/CFluidFlamelet.hpp"
#include "../../include/solvers/CSpeciesSolver.hpp"
#include "../../include/variables/CFlowVariable.hpp"
#include "../../include/variables/CSpeciesFlameletVariable.hpp"

CSpeciesFlameletSolver::CSpeciesFlameletSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh)
    : CSpeciesSolver(geometry, config, true) {

  /*--- Retrieve options from config. ---*/
  flamelet_config_options = config->GetFlameletParsedOptions();
  
  /*--- Dimension of the problem. ---*/
  nVar = flamelet_config_options.n_scalars;
  include_mixture_fraction = (flamelet_config_options.n_control_vars == 3);

  Initialize(geometry, config, iMesh, nVar);

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CSpeciesFlameletVariable(Solution_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- Store the initial CFL number for all grid points. ---*/

  const su2double CFL = config->GetCFL(MGLevel) * config->GetCFLRedCoeff_Species();
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0u; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  END_SU2_OMP_FOR
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name. ---*/
  SolverName = "FLAMELET";
}

void CSpeciesFlameletSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                           unsigned short iMesh, unsigned short iRKStep,
                                           unsigned short RunTime_EqSystem, bool Output) {
  unsigned long n_not_in_domain_local = 0, n_not_in_domain_global = 0;
  vector<su2double> scalars_vector(nVar);
  unsigned long spark_iter_start, spark_duration;
  bool ignition = false;
  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

  /*--- Retrieve spark ignition parameters for spark-type ignition. ---*/
  if ((flamelet_config_options.ignition_method == FLAMELET_INIT_TYPE::SPARK)) {
    auto spark_init = flamelet_config_options.spark_init;
    spark_iter_start = ceil(spark_init[4]);
    spark_duration = ceil(spark_init[5]);
    unsigned long iter = config->GetMultizone_Problem() ? config->GetOuterIter() : config->GetInnerIter();
    ignition = ((iter >= spark_iter_start) && (iter <= (spark_iter_start + spark_duration)));
  }

  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto i_point = 0u; i_point < nPoint; i_point++) {
    CFluidModel* fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();
    su2double* scalars = nodes->GetSolution(i_point);
    for (auto iVar = 0u; iVar < nVar; iVar++) scalars_vector[iVar] = scalars[iVar];

    /*--- Compute total source terms from the production and consumption. ---*/
    unsigned long misses = SetScalarSources(config, fluid_model_local, i_point, scalars_vector);

    if (ignition) {
      /*--- Apply source terms within spark radius. ---*/
      su2double dist_from_center = 0,
                spark_radius = flamelet_config_options.spark_init[3];
      dist_from_center = GeometryToolbox::SquaredDistance(nDim, geometry->nodes->GetCoord(i_point), flamelet_config_options.spark_init.data());
      if (dist_from_center < pow(spark_radius,2)) {
        for (auto iVar = 0u; iVar < nVar; iVar++)
          nodes->SetScalarSource(i_point, iVar, nodes->GetScalarSources(i_point)[iVar] + flamelet_config_options.spark_reaction_rates[iVar]);
      }
    }

    nodes->SetTableMisses(i_point, misses);
    n_not_in_domain_local += misses;
    /*--- Obtain passive look-up scalars. ---*/
    SetScalarLookUps(fluid_model_local, i_point, scalars_vector);

    /*--- Set mass diffusivity based on thermodynamic state. ---*/
    auto T = flowNodes->GetTemperature(i_point);
    fluid_model_local->SetTDState_T(T, scalars);
    /*--- set the diffusivity in the fluid model to the diffusivity obtained from the lookup table ---*/
    for (auto i_scalar = 0u; i_scalar < nVar; ++i_scalar) {
      nodes->SetDiffusivity(i_point, fluid_model_local->GetMassDiffusivity(i_scalar), i_scalar);
    }

    /*--- Obtain preferential diffusion scalar values. ---*/
    if (flamelet_config_options.preferential_diffusion)
      SetPreferentialDiffusionScalars(fluid_model_local, i_point, scalars_vector);

    if (!Output) LinSysRes.SetBlock_Zero(i_point);
  }
  END_SU2_OMP_FOR
  /* --- Sum up some global counters over processes. --- */
  SU2_MPI::Reduce(&n_not_in_domain_local, &n_not_in_domain_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE,
                  SU2_MPI::GetComm());
  if ((rank == MASTER_NODE) && (n_not_in_domain_global > 0))
    cout << "Number of points outside manifold domain: " << n_not_in_domain_global << endl;

  /*--- Compute preferential diffusion scalar gradients. ---*/
  if (flamelet_config_options.preferential_diffusion) {
    switch (config->GetKind_Gradient_Method()) {
      case GREEN_GAUSS:
        SetAuxVar_Gradient_GG(geometry, config);
        break;
      case WEIGHTED_LEAST_SQUARES:
        SetAuxVar_Gradient_LS(geometry, config);
        break;
      default:
        break;
    }
  }
  /*--- Clear Residual and Jacobian. Upwind second order reconstruction and gradients ---*/
  CommonPreprocessing(geometry, config, Output);
}

void CSpeciesFlameletSolver::SetInitialCondition(CGeometry** geometry, CSolver*** solver_container, CConfig* config,
                                                 unsigned long ExtIter) {
  const bool Restart = (config->GetRestart() || config->GetRestart_Flow());

  if ((!Restart) && ExtIter == 0) {
    if (rank == MASTER_NODE) {
      cout << "Initializing progress variable and total enthalpy (using temperature)" << endl;
    }

    su2double flame_offset[3] = {0, 0, 0}, flame_normal[3] = {0, 0, 0}, flame_thickness = 0, flame_burnt_thickness = 0,
              flamenorm = 0;
    bool flame_front_ignition = (flamelet_config_options.ignition_method == FLAMELET_INIT_TYPE::FLAME_FRONT);

    if (flame_front_ignition) {
      /*--- Collect flame front ignition parameters. ---*/
      auto flame_init = flamelet_config_options.flame_init;
      for (auto iDim = 0u; iDim < 3; ++iDim) {
        flame_offset[iDim] = flame_init[iDim];
        flame_normal[iDim] = flame_init[3 + iDim];
      }
      flame_thickness = flame_init[6];
      flame_burnt_thickness = flame_init[7];
      flamenorm = GeometryToolbox::Norm(nDim, flame_normal);
    }

    const su2double temp_inlet = config->GetInc_Temperature_Init();
    su2double enth_inlet = config->GetSpecies_Init()[I_ENTH];

    su2double prog_burnt = 0, prog_unburnt, point_loc;
    su2double scalar_init[MAXNVAR];

    if (rank == MASTER_NODE) {
      cout << "initial condition: T = " << temp_inlet << endl;
      for (auto iCV = 0u; iCV < flamelet_config_options.n_control_vars; iCV++) {
        const auto& cv_name = flamelet_config_options.controlling_variable_names[iCV];
        cout << "initial condition: " << cv_name << " = " << config->GetSpecies_Init()[iCV] << endl;
      }
      switch (flamelet_config_options.ignition_method) {
        case FLAMELET_INIT_TYPE::FLAME_FRONT:
          cout << "Ignition with a straight flame front" << endl;
          break;
        case FLAMELET_INIT_TYPE::SPARK:
          cout << "Ignition with an artificial spark" << endl;
          break;
        case FLAMELET_INIT_TYPE::NONE:
          cout << "No solution ignition (cold flow)" << endl;
          break;
        default:
          break;
      }
    }

    CFluidModel* fluid_model_local;

    unsigned long n_not_iterated_local = 0, n_not_in_domain_local = 0, n_points_unburnt_local = 0,
                  n_points_burnt_local = 0, n_points_flame_local = 0, n_not_iterated_global, n_not_in_domain_global,
                  n_points_burnt_global, n_points_flame_global, n_points_unburnt_global;

    for (unsigned long i_mesh = 0; i_mesh <= config->GetnMGLevels(); i_mesh++) {
      fluid_model_local = solver_container[i_mesh][FLOW_SOL]->GetFluidModel();

      for (auto iVar = 0u; iVar < nVar; iVar++) scalar_init[iVar] = config->GetSpecies_Init()[iVar];

      /*--- Set enthalpy based on initial temperature and scalars. ---*/
      n_not_iterated_local += GetEnthFromTemp(fluid_model_local, temp_inlet, config->GetSpecies_Init(), &enth_inlet);
      scalar_init[I_ENTH] = enth_inlet;

      prog_unburnt = config->GetSpecies_Init()[I_PROGVAR];
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (unsigned long i_point = 0; i_point < nPoint; i_point++) {
        auto coords = geometry[i_mesh]->nodes->GetCoord(i_point);

        if (flame_front_ignition) {

          prog_burnt = GetBurntProgressVariable(fluid_model_local, scalar_init);
          
          /*--- Determine if point is above or below the plane, assuming the normal
            is pointing towards the burned region. ---*/
          point_loc = 0.0;
          for (unsigned short i_dim = 0; i_dim < geometry[i_mesh]->GetnDim(); i_dim++) {
            point_loc += flame_normal[i_dim] * (coords[i_dim] - flame_offset[i_dim]);
          }

          /*--- Compute the exact distance from point to plane. ---*/
          point_loc = point_loc / flamenorm;

          /* --- Unburnt region upstream of the flame. --- */
          if (point_loc <= 0) {
            scalar_init[I_PROGVAR] = prog_unburnt;
            n_points_unburnt_local++;

            /* --- Flame zone where we lineary increase from unburnt to burnt conditions. --- */
          } else if ((point_loc > 0) && (point_loc <= flame_thickness)) {
            scalar_init[I_PROGVAR] = prog_unburnt + point_loc * (prog_burnt - prog_unburnt) / flame_thickness;
            n_points_flame_local++;

            /* --- Burnt region behind the flame zone. --- */
          } else if ((point_loc > flame_thickness) && (point_loc <= flame_thickness + flame_burnt_thickness)) {
            scalar_init[I_PROGVAR] = prog_burnt;
            n_points_burnt_local++;

            /* --- Unburnt region downstream of the flame and burnt region. --- */
          } else {
            scalar_init[I_PROGVAR] = prog_unburnt;
            n_points_unburnt_local++;
          }
        } else {
          scalar_init[I_PROGVAR] = prog_unburnt;
        }
        /* --- Perform manifold evaluation to check whether initial scalar is out-of-bounds. --- */
        fluid_model_local->SetTDState_T(temp_inlet, scalar_init);
        n_not_in_domain_local += fluid_model_local->GetExtrapolation();

        /* --- Initialize the auxiliary transported scalars  (not controlling variables). --- */
        for (int i_scalar = flamelet_config_options.n_control_vars; i_scalar < flamelet_config_options.n_scalars; ++i_scalar) {
          scalar_init[i_scalar] = config->GetSpecies_Init()[i_scalar];
        }

        solver_container[i_mesh][SPECIES_SOL]->GetNodes()->SetSolution(i_point, scalar_init);
      }

      solver_container[i_mesh][SPECIES_SOL]->InitiateComms(geometry[i_mesh], config, MPI_QUANTITIES::SOLUTION);
      solver_container[i_mesh][SPECIES_SOL]->CompleteComms(geometry[i_mesh], config, MPI_QUANTITIES::SOLUTION);

      solver_container[i_mesh][FLOW_SOL]->InitiateComms(geometry[i_mesh], config, MPI_QUANTITIES::SOLUTION);
      solver_container[i_mesh][FLOW_SOL]->CompleteComms(geometry[i_mesh], config, MPI_QUANTITIES::SOLUTION);

      solver_container[i_mesh][FLOW_SOL]->Preprocessing(geometry[i_mesh], solver_container[i_mesh], config, i_mesh,
                                                        NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      END_SU2_OMP_FOR
    }

    /* --- Sum up some global counters over processes. --- */
    SU2_MPI::Reduce(&n_not_in_domain_local, &n_not_in_domain_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE,
                    SU2_MPI::GetComm());
    SU2_MPI::Reduce(&n_not_iterated_local, &n_not_iterated_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE,
                    SU2_MPI::GetComm());
    SU2_MPI::Reduce(&n_points_unburnt_local, &n_points_unburnt_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE,
                    SU2_MPI::GetComm());
    SU2_MPI::Reduce(&n_points_burnt_local, &n_points_burnt_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE,
                    SU2_MPI::GetComm());
    SU2_MPI::Reduce(&n_points_flame_local, &n_points_flame_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE,
                    SU2_MPI::GetComm());

    if (rank == MASTER_NODE) {
      cout << endl;
      if (flame_front_ignition) {
        cout << " Number of points in unburnt region: " << n_points_unburnt_global << "." << endl;
        cout << " Number of points in burnt region  : " << n_points_burnt_global << "." << endl;
        cout << " Number of points in flame zone    : " << n_points_flame_global << "." << endl;
      }

      if (n_not_in_domain_global > 0)
        cout << " Initial condition: Number of points outside of table domain: " << n_not_in_domain_global << " !!!"
             << endl;

      if (n_not_iterated_global > 0)
        cout << " Initial condition: Number of points in which enthalpy could not be iterated: "
             << n_not_iterated_global << " !!!" << endl;
    }
  }
}

void CSpeciesFlameletSolver::SetPreconditioner(CGeometry* geometry, CSolver** solver_container, CConfig* config) {
  const bool variable_density = (config->GetVariable_Density_Model());
  const bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0u; iPoint < nPointDomain; iPoint++) {
    /*--- Access the primitive variables at this node. ---*/

    su2double Density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    su2double BetaInc2 = solver_container[FLOW_SOL]->GetNodes()->GetBetaInc2(iPoint);
    su2double Temperature = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);

    su2double SolP = solver_container[FLOW_SOL]->LinSysSol(iPoint, prim_idx.Pressure());
    su2double SolT = solver_container[FLOW_SOL]->LinSysSol(iPoint, prim_idx.Temperature());

    /*--- We need the derivative of the equation of state to build the
     preconditioning matrix. For now, the only option is the ideal gas
     law, but in the future, dRhodT should be in the fluid model. ---*/

    su2double dRhodT = 0.0;
    if (variable_density) {
      dRhodT = -Density / Temperature;
    }

    /*--- Passive scalars have no impact on the density. ---*/

    su2double dRhodC = 0.0;

    /*--- Modify matrix diagonal with term including volume and time step. ---*/

    su2double Vol = geometry->nodes->GetVolume(iPoint);
    su2double Delta =
        Vol / (config->GetCFLRedCoeff_Species() * solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint));

    /*--- Calculating the inverse of the preconditioning matrix
     that multiplies the time derivative during time integration. ---*/

    if (implicit) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        su2double scalar = nodes->GetSolution(iPoint, iVar);

        /*--- Compute the lag terms for the decoupled linear system from
         the mean flow equations and add to the residual for the scalar.
         In short, we are effectively making these terms explicit. ---*/

        su2double artcompc1 = SolP * scalar / (Density * BetaInc2);
        su2double artcompc2 = SolT * dRhodT * scalar / (Density);

        LinSysRes(iPoint, iVar) += artcompc1 + artcompc2;

        /*--- Add the extra Jacobian term to the scalar system. ---*/

        su2double Jaccomp = scalar * dRhodC + Density;
        su2double JacTerm = Jaccomp * Delta;

        Jacobian.AddVal2Diag(iPoint, iVar, JacTerm);
      }
    }
  }
  END_SU2_OMP_FOR
}

void CSpeciesFlameletSolver::Source_Residual(CGeometry* geometry, CSolver** solver_container,
                                             CNumerics** numerics_container, CConfig* config, unsigned short iMesh) {
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto i_point = 0u; i_point < nPointDomain; i_point++) {
    /*--- Add source terms from the lookup table directly to the residual. ---*/
    for (auto i_var = 0; i_var < nVar; i_var++) {
      LinSysRes(i_point, i_var) -= nodes->GetScalarSources(i_point)[i_var] * geometry->nodes->GetVolume(i_point);
    }
  }
  END_SU2_OMP_FOR

  /*--- call the species solver for the shared sources (axisymmetric) ---*/
  CSpeciesSolver::Source_Residual(geometry, solver_container, numerics_container, config, iMesh);
}

void CSpeciesFlameletSolver::BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                      CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double temp_inlet = config->GetInletTtotal(Marker_Tag);

  /*--- We compute inlet enthalpy from the temperature and progress variable. ---*/
  su2double enth_inlet;
  GetEnthFromTemp(solver_container[FLOW_SOL]->GetFluidModel(), temp_inlet, config->GetInlet_SpeciesVal(Marker_Tag),
                  &enth_inlet);

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    Inlet_SpeciesVars[val_marker][iVertex][I_ENTH] = enth_inlet;
  }
  END_SU2_OMP_FOR

  /*--- Call the general inlet boundary condition implementation. ---*/
  CSpeciesSolver::BC_Inlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CSpeciesFlameletSolver::BC_Isothermal_Wall_Generic(CGeometry* geometry, CSolver** solver_container,
                                                        CNumerics* conv_numerics, CNumerics* visc_numerics,
                                                        CConfig* config, unsigned short val_marker, bool cht_mode) {
  const bool implicit = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  const string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  CFluidModel* fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();
  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());
  unsigned long n_not_iterated = 0;

  /*--- Loop over all the vertices on this boundary marker. ---*/
  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (unsigned long iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    unsigned long iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    su2double temp_wall, enth_wall;
    if (cht_mode)
      temp_wall = solver_container[FLOW_SOL]->GetConjugateHeatVariable(val_marker, iVertex, 0);
    else
      temp_wall = config->GetIsothermal_Temperature(Marker_Tag);

    /*--- Check if the node belongs to the domain (i.e., not a halo node). ---*/

    if (geometry->nodes->GetDomain(iPoint)) {
      if (config->GetMarker_StrongBC(Marker_Tag) == true) {
        /*--- Initial guess for enthalpy value. ---*/
        enth_wall = nodes->GetSolution(iPoint, I_ENTH);

        /*--- Set enthalpy on the wall. ---*/
        n_not_iterated += GetEnthFromTemp(fluid_model_local, temp_wall, nodes->GetSolution(iPoint), &enth_wall);

        /*--- Impose the value of the enthalpy as a strong boundary
        condition (Dirichlet) and remove any
        contribution to the residual at this node. ---*/

        nodes->SetSolution(iPoint, I_ENTH, enth_wall);
        nodes->SetSolution_Old(iPoint, I_ENTH, enth_wall);

        LinSysRes(iPoint, I_ENTH) = 0.0;

        nodes->SetVal_ResTruncError_Zero(iPoint, I_ENTH);

        if (implicit) {
          unsigned long total_index = iPoint * nVar + I_ENTH;
          Jacobian.DeleteValsRowi(total_index);
        }
      } else {
        /*--- Weak BC formulation. ---*/
        const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

        const su2double Area = GeometryToolbox::Norm(nDim, Normal);

        const auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        /*--- Get coordinates of i & nearest normal and compute distance. ---*/

        const auto Coord_i = geometry->nodes->GetCoord(iPoint);
        const auto Coord_j = geometry->nodes->GetCoord(Point_Normal);
        su2double Edge_Vector[MAXNDIM];
        GeometryToolbox::Distance(nDim, Coord_j, Coord_i, Edge_Vector);
        su2double dist_ij_2 = GeometryToolbox::SquaredNorm(nDim, Edge_Vector);
        su2double dist_ij = sqrt(dist_ij_2);

        /*--- Compute the normal gradient in temperature using Twall. ---*/
        ///TODO: Account for preferential diffusion in computation of the heat flux
        su2double dTdn = -(flowNodes->GetTemperature(Point_Normal) - temp_wall) / dist_ij;

        /*--- Get thermal conductivity. ---*/

        su2double thermal_conductivity = flowNodes->GetThermalConductivity(iPoint);

        /*--- Apply a weak boundary condition for the energy equation.
        Compute the residual due to the prescribed heat flux. ---*/

        LinSysRes(iPoint, I_ENTH) -= thermal_conductivity * dTdn * Area;
      }
    }
  }
  END_SU2_OMP_FOR

  if (rank == MASTER_NODE && n_not_iterated > 0) {
    cout << " !!! Wall bc (" << Marker_Tag
         << "): Number of points in which enthalpy could not be iterated: " << n_not_iterated << " !!!" << endl;
  }
}

void CSpeciesFlameletSolver::BC_Isothermal_Wall(CGeometry* geometry, CSolver** solver_container,
                                                CNumerics* conv_numerics, CNumerics* visc_numerics, CConfig* config,
                                                unsigned short val_marker) {
  BC_Isothermal_Wall_Generic(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CSpeciesFlameletSolver::BC_ConjugateHeat_Interface(CGeometry* geometry, CSolver** solver_container,
                                                        CNumerics* conv_numerics, CConfig* config,
                                                        unsigned short val_marker) {
  BC_Isothermal_Wall_Generic(geometry, solver_container, conv_numerics, nullptr, config, val_marker, true);
}

unsigned long CSpeciesFlameletSolver::SetScalarSources(const CConfig* config, CFluidModel* fluid_model_local,
                                                       unsigned long iPoint, const vector<su2double>& scalars) {
  /*--- Compute total source terms from the production and consumption. ---*/

  vector<su2double> table_sources(flamelet_config_options.n_control_vars + 2 * flamelet_config_options.n_user_scalars);
  unsigned long misses = fluid_model_local->EvaluateDataSet(scalars, FLAMELET_LOOKUP_OPS::SOURCES, table_sources);
  table_sources[I_PROGVAR] = fmax(0, table_sources[I_PROGVAR]);
  nodes->SetTableMisses(iPoint, misses);

  /*--- The source term for progress variable is always positive, we clip from below to makes sure. --- */

  vector<su2double> source_scalar(flamelet_config_options.n_scalars);
  for (auto iCV = 0u; iCV < flamelet_config_options.n_control_vars; iCV++) source_scalar[iCV] = table_sources[iCV];

  /*--- Source term for the auxiliary species transport equations. ---*/
  for (size_t i_aux = 0; i_aux < flamelet_config_options.n_user_scalars; i_aux++) {
    /*--- The source term for the auxiliary equations consists of a production term and a consumption term:
          S_TOT = S_PROD + S_CONS * Y ---*/
    su2double y_aux = scalars[flamelet_config_options.n_control_vars + i_aux];
    su2double source_prod = table_sources[flamelet_config_options.n_control_vars + 2 * i_aux];
    su2double source_cons = table_sources[flamelet_config_options.n_control_vars + 2 * i_aux + 1];
    source_scalar[flamelet_config_options.n_control_vars + i_aux] = source_prod + source_cons * y_aux;
  }
  for (auto i_scalar = 0u; i_scalar < nVar; i_scalar++)
    nodes->SetScalarSource(iPoint, i_scalar, source_scalar[i_scalar]);
  return misses;
}

unsigned long CSpeciesFlameletSolver::SetScalarLookUps(CFluidModel* fluid_model_local,
                                                       unsigned long iPoint, const vector<su2double>& scalars) {
  /*--- Retrieve the passive look-up variables from the manifold. ---*/
  unsigned long misses{0};
  /*--- Skip if no passive look-ups are listed ---*/
  if (flamelet_config_options.n_lookups > 0) {
    vector<su2double> lookup_scalar(flamelet_config_options.n_lookups);
    misses = fluid_model_local->EvaluateDataSet(scalars, FLAMELET_LOOKUP_OPS::LOOKUP, lookup_scalar);

    for (auto i_lookup = 0u; i_lookup < flamelet_config_options.n_lookups; i_lookup++) {
      nodes->SetLookupScalar(iPoint, lookup_scalar[i_lookup], i_lookup);
    }
  }

  return misses;
}

unsigned long CSpeciesFlameletSolver::SetPreferentialDiffusionScalars(CFluidModel* fluid_model_local,
                                                                      unsigned long iPoint,
                                                                      const vector<su2double>& scalars) {
  /*--- Retrieve the preferential diffusion scalar values from the manifold. ---*/

  vector<su2double> beta_scalar(FLAMELET_PREF_DIFF_SCALARS::N_BETA_TERMS);
  unsigned long misses = fluid_model_local->EvaluateDataSet(scalars, FLAMELET_LOOKUP_OPS::PREFDIF, beta_scalar);

  for (auto i_beta = 0u; i_beta < FLAMELET_PREF_DIFF_SCALARS::N_BETA_TERMS; i_beta++) {
    nodes->SetAuxVar(iPoint, i_beta, beta_scalar[i_beta]);
  }
  return misses;
}

void CSpeciesFlameletSolver::Viscous_Residual(const unsigned long iEdge, const CGeometry* geometry, CSolver** solver_container,
                                              CNumerics* numerics, const CConfig* config) {
  /*--- Overloaded viscous residual method which accounts for preferential diffusion.  ---*/
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT),
             PreferentialDiffusion = flamelet_config_options.preferential_diffusion;

  /*--- Points in edge ---*/
  auto iPoint = geometry->edges->GetNode(iEdge, 0);
  auto jPoint = geometry->edges->GetNode(iEdge, 1);

  auto SolverSpecificNumerics = [&](unsigned long iPoint, unsigned long jPoint) {
    /*--- Mass diffusivity coefficients. ---*/

    numerics->SetDiffusionCoeff(nodes->GetDiffusivity(iPoint), nodes->GetDiffusivity(jPoint));
  };

  /*--- Regular viscous scalar residual computation. ---*/

  Viscous_Residual_impl(SolverSpecificNumerics, iEdge, geometry, solver_container, numerics, config);

  /*--- Viscous residual due to preferential diffusion ---*/
  if (PreferentialDiffusion) {
    CFlowVariable* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

    su2double scalar_i[MAXNVAR] = {0},
              scalar_j[MAXNVAR] = {0},
              diff_coeff_beta_i[MAXNVAR] = {0},
              diff_coeff_beta_j[MAXNVAR] = {0};

    // Number of active transport scalars
    const auto n_CV = flamelet_config_options.n_control_vars;

    su2activematrix scalar_grad_i(MAXNVAR, MAXNDIM), scalar_grad_j(MAXNVAR, MAXNDIM);
    /*--- Looping over spatial dimensions to fill in the diffusion scalar gradients. ---*/
    /*--- The scalar gradient is subtracted to account for regular viscous diffusion. ---*/
    for (auto iScalar = 0u; iScalar < n_CV; ++iScalar) {
      for (auto iDim = 0u; iDim < nDim; ++iDim) {
        switch (iScalar) {
          case I_PROGVAR:
            scalar_grad_i[iScalar][iDim] =
                nodes->GetAuxVarGradient(iPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_PROGVAR, iDim) -
                nodes->GetGradient(iPoint, I_PROGVAR, iDim);
            scalar_grad_j[iScalar][iDim] =
                nodes->GetAuxVarGradient(jPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_PROGVAR, iDim) -
                nodes->GetGradient(jPoint, I_PROGVAR, iDim);
            break;
          case I_ENTH:
            scalar_grad_i[iScalar][iDim] =
                nodes->GetAuxVarGradient(iPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH, iDim) -
                nodes->GetGradient(iPoint, I_ENTH, iDim);
            scalar_grad_j[iScalar][iDim] =
                nodes->GetAuxVarGradient(jPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH, iDim) -
                nodes->GetGradient(jPoint, I_ENTH, iDim);
            break;
          case I_MIXFRAC:
            scalar_grad_i[iScalar][iDim] =
                nodes->GetAuxVarGradient(iPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_MIXFRAC, iDim) -
                nodes->GetGradient(iPoint, I_MIXFRAC, iDim);
            scalar_grad_j[iScalar][iDim] =
                nodes->GetAuxVarGradient(jPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_MIXFRAC, iDim) -
                nodes->GetGradient(jPoint, I_MIXFRAC, iDim);
            break;
          default:
            break;
        }
      }
    }
    /*--- No preferential diffusion modification for passive species. ---*/
    for (auto iScalar = n_CV; iScalar < nVar; ++iScalar) {
      for (auto iDim = 0u; iDim < nDim; ++iDim) {
        scalar_grad_i[iScalar][iDim] = 0;
        scalar_grad_j[iScalar][iDim] = 0;
      }
    }

    for (auto iScalar = 0u; iScalar < n_CV; ++iScalar) {
      /*--- Filling in the preferential diffusion scalars (beta_pv, beta_h2, beta_Z). ---*/
      switch (iScalar) {
        case I_PROGVAR:
          scalar_i[iScalar] = nodes->GetAuxVar(iPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_PROGVAR) -
                              nodes->GetSolution(iPoint, iScalar);
          scalar_j[iScalar] = nodes->GetAuxVar(jPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_PROGVAR) -
                              nodes->GetSolution(jPoint, iScalar);
          break;
        case I_ENTH:
          scalar_i[iScalar] =
              nodes->GetAuxVar(iPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH) - nodes->GetSolution(iPoint, iScalar);
          scalar_j[iScalar] =
              nodes->GetAuxVar(jPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH) - nodes->GetSolution(jPoint, iScalar);
          break;
        case I_MIXFRAC:
          scalar_i[iScalar] = nodes->GetAuxVar(iPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_MIXFRAC) -
                              nodes->GetSolution(iPoint, iScalar);
          scalar_j[iScalar] = nodes->GetAuxVar(jPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_MIXFRAC) -
                              nodes->GetSolution(jPoint, iScalar);
          break;
        default:
          break;
      }
      diff_coeff_beta_i[iScalar] = nodes->GetDiffusivity(iPoint, iScalar);
      diff_coeff_beta_j[iScalar] = nodes->GetDiffusivity(jPoint, iScalar);
    }

    for (auto iScalar = n_CV; iScalar < nVar; ++iScalar) {
      scalar_i[iScalar] = 0;
      scalar_j[iScalar] = 0;
      diff_coeff_beta_i[iScalar] = 0;
      diff_coeff_beta_j[iScalar] = 0;
    }

    numerics->SetScalarVar(scalar_i, scalar_j);

    numerics->SetScalarVarGradient(CMatrixView<su2double>(scalar_grad_i), CMatrixView<su2double>(scalar_grad_j));

    numerics->SetDiffusionCoeff(diff_coeff_beta_i, diff_coeff_beta_j);

    /*--- Computing first preferential residual component. ---*/
    auto residual_PD = numerics->ComputeResidual(config);

    if (ReducerStrategy) {
      EdgeFluxes.SubtractBlock(iEdge, residual_PD);

      if (implicit) Jacobian.UpdateBlocksSub(iEdge, residual_PD.jacobian_i, residual_PD.jacobian_j);
    } else {
      LinSysRes.SubtractBlock(iPoint, residual_PD);
      LinSysRes.AddBlock(jPoint, residual_PD);
      /*--- Set implicit computation ---*/
      if (implicit) Jacobian.UpdateBlocksSub(iEdge, iPoint, jPoint, residual_PD.jacobian_i, residual_PD.jacobian_j);
    }

    /* Computing the second preferential diffusion terms due to heat flux */
    for (auto iScalar = 0u; iScalar < nVar; ++iScalar) {
      for (auto iDim = 0u; iDim < nDim; ++iDim) {
        if (iScalar == I_ENTH) {
          /* Setting the temperature gradient */
          scalar_grad_i[iScalar][iDim] = flowNodes->GetGradient_Primitive(iPoint, prim_idx.Temperature(), iDim);
          scalar_grad_j[iScalar][iDim] = flowNodes->GetGradient_Primitive(jPoint, prim_idx.Temperature(), iDim);
        } else {
          scalar_grad_i[iScalar][iDim] = 0;
          scalar_grad_j[iScalar][iDim] = 0;
        }
      }

      if (iScalar == I_ENTH) {
        scalar_i[iScalar] = flowNodes->GetTemperature(iPoint);
        scalar_j[iScalar] = flowNodes->GetTemperature(jPoint);
        diff_coeff_beta_i[iScalar] = nodes->GetAuxVar(iPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH_THERMAL) *
                                     nodes->GetDiffusivity(iPoint, iScalar);
        diff_coeff_beta_j[iScalar] = nodes->GetAuxVar(jPoint, FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH_THERMAL) *
                                     nodes->GetDiffusivity(jPoint, iScalar);
      } else {
        scalar_i[iScalar] = 0;
        scalar_j[iScalar] = 0;
        diff_coeff_beta_i[iScalar] = 0;
        diff_coeff_beta_j[iScalar] = 0;
      }
    }

    numerics->SetScalarVar(scalar_i, scalar_j);

    numerics->SetScalarVarGradient(CMatrixView<su2double>(scalar_grad_i), CMatrixView<su2double>(scalar_grad_j));

    numerics->SetDiffusionCoeff(diff_coeff_beta_i, diff_coeff_beta_j);

    auto residual_thermal = numerics->ComputeResidual(config);

    if (ReducerStrategy) {
      EdgeFluxes.SubtractBlock(iEdge, residual_thermal);
    } else {
      LinSysRes.SubtractBlock(iPoint, residual_thermal);
      LinSysRes.AddBlock(jPoint, residual_thermal);
      /* No implicit part for the preferential diffusion of heat */
    }
  }
}

unsigned long CSpeciesFlameletSolver::GetEnthFromTemp(CFluidModel* fluid_model, su2double const val_temp,
                                                      const su2double* scalar_solution, su2double* val_enth) {
  /*--- convergence criterion for temperature in [K], high accuracy needed for restarts. ---*/
  su2double delta_temp_final = 0.001;
  su2double enth_iter = scalar_solution[I_ENTH];
  su2double delta_enth;
  su2double delta_temp_iter = 1e10;
  unsigned long exit_code = 0;
  const int counter_limit = 1000;

  int counter = 0;

  su2double val_scalars[MAXNVAR];
  for (auto iVar = 0u; iVar < nVar; iVar++) val_scalars[iVar] = scalar_solution[iVar];

  while ((abs(delta_temp_iter) > delta_temp_final) && (counter++ < counter_limit)) {
    /*--- Add all quantities and their names to the look up vectors. ---*/
    val_scalars[I_ENTH] = enth_iter;
    fluid_model->SetTDState_T(val_temp, val_scalars);

    su2double Temperature = fluid_model->GetTemperature();
    su2double Cp = fluid_model->GetCp();

    delta_temp_iter = val_temp - Temperature;

    delta_enth = Cp * delta_temp_iter;

    enth_iter += delta_enth;
  }

  *val_enth = enth_iter;

  if (counter >= counter_limit) {
    exit_code = 1;
  }

  return exit_code;
}

su2double CSpeciesFlameletSolver::GetBurntProgressVariable(CFluidModel* fluid_model, const su2double* scalar_solution) {
  su2double scalars[MAXNVAR], delta = 1e-3;
  for (auto iVar = 0u; iVar < nVar; iVar++) scalars[iVar] = scalar_solution[iVar];

  bool outside = false;
  while (!outside) {
    fluid_model->SetTDState_T(300, scalars);
    if (fluid_model->GetExtrapolation() == 1) outside = true;
    scalars[I_PROGVAR] += delta;
  }
  su2double pv_burnt = scalars[I_PROGVAR] - delta;
  return pv_burnt;
}
