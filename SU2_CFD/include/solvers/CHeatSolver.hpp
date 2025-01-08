/*!
 * \file CHeatSolver.hpp
 * \brief Headers of the CHeatSolver class
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

#pragma once

#include "CScalarSolver.hpp"
#include "../variables/CHeatVariable.hpp"

/*!
 * \class CHeatSolver
 * \brief Main class for defining the finite-volume heat solver.
 * \author O. Burghardt
 * \version 8.0.0 "Harrier"
 */
class CHeatSolver final : public CScalarSolver<CHeatVariable> {
protected:
  static constexpr size_t MAXNDIM = 3; /*!< \brief Max number of space dimensions, used in some static arrays. */
  static constexpr size_t MAXNVAR = 1; /*!< \brief Max number of variables, for static arrays. */

  const bool flow; /*!< \brief Use solver as a scalar transport equation of Temperature for the inc solver. */
  const bool heat_equation; /*!< \brief use solver for heat conduction in solids. */

  su2double Global_Delta_Time = 0.0, Global_Delta_UnstTimeND = 0.0;

  vector<vector<su2double>> HeatFlux;
  vector<su2double> HeatFlux_per_Marker;
  su2double Total_HeatFlux;
  su2double AllBound_HeatFlux;
  vector<su2double> AverageT_per_Marker;
  su2double Total_AverageT;
  su2double AllBound_AverageT;
  vector<su2double> Surface_Areas;
  su2double Total_HeatFlux_Areas;
  su2double Total_HeatFlux_Areas_Monitor;
  vector<su2activematrix> ConjugateVar;

  /*!
   * \brief Applies an isothermal condition to a vertex of a marker.
   */
  void IsothermalBoundaryCondition(CGeometry *geometry, CSolver *flow_solver, const CConfig *config,
                                   unsigned short iMarker, unsigned long iVertex, const su2double& temperature) {

    const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
    if (!geometry->nodes->GetDomain(iPoint)) return;

    const bool implicit = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;
    const su2double prandtl_lam = config->GetPrandtl_Lam();
    const su2double const_diffusivity = config->GetThermalDiffusivity();

    const auto Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

    const auto* Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
    const su2double Area = GeometryToolbox::Norm(nDim, Normal);

    const auto* Coord_i = geometry->nodes->GetCoord(iPoint);
    const auto* Coord_j = geometry->nodes->GetCoord(Point_Normal);
    const su2double dist_ij = GeometryToolbox::Distance(nDim, Coord_i, Coord_j);

    const su2double dTdn = -(nodes->GetTemperature(Point_Normal) - temperature) / dist_ij;

    su2double thermal_diffusivity = const_diffusivity;
    if (flow) {
      thermal_diffusivity = flow_solver->GetNodes()->GetLaminarViscosity(iPoint) / prandtl_lam;
    }
    LinSysRes(iPoint, 0) -= thermal_diffusivity * dTdn * Area;

    if (implicit) {
      su2double Jacobian_i[] = {-thermal_diffusivity / dist_ij * Area};
      Jacobian.SubtractBlock2Diag(iPoint, &Jacobian_i);
    }
  }

  /*!
   * \brief Compute the viscous flux for the scalar equation at a particular edge.
   * \param[in] iEdge - Edge for which we want to compute the flux
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \note Calls a generic implementation after defining a SolverSpecificNumerics object.
   */
  inline void Viscous_Residual(unsigned long iEdge, CGeometry* geometry, CSolver** solver_container,
                                       CNumerics* numerics, CConfig* config) override {
    const CVariable* flow_nodes = flow ? solver_container[FLOW_SOL]->GetNodes() : nullptr;

    const su2double const_diffusivity = config->GetThermalDiffusivity();
    const su2double pr_lam = config->GetPrandtl_Lam();
    const su2double pr_turb = config->GetPrandtl_Turb();

    su2double thermal_diffusivity_i{}, thermal_diffusivity_j{};

    /*--- Computes the thermal diffusivity to use in the viscous numerics. ---*/
    auto compute_thermal_diffusivity = [&](unsigned long iPoint, unsigned long jPoint) {
      if (flow) {
        thermal_diffusivity_i = flow_nodes->GetLaminarViscosity(iPoint) / pr_lam +
                                flow_nodes->GetEddyViscosity(iPoint) / pr_turb;
        thermal_diffusivity_j = flow_nodes->GetLaminarViscosity(jPoint) / pr_lam +
                                flow_nodes->GetEddyViscosity(jPoint) / pr_turb;
        numerics->SetDiffusionCoeff(&thermal_diffusivity_i, &thermal_diffusivity_j);
      }
      else {
        numerics->SetDiffusionCoeff(&const_diffusivity, &const_diffusivity);
      }
    };
    /*--- Compute residual and Jacobians. ---*/
    Viscous_Residual_impl(compute_thermal_diffusivity, iEdge, geometry, solver_container, numerics, config);
  }

public:

  /*!
   * \brief Constructor of the class.
   */
  CHeatSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Restart residual and compute gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry,
                    CSolver **solver_container,
                    CConfig *config,
                    unsigned short iMesh,
                    unsigned short iRKStep,
                    unsigned short RunTime_EqSystem,
                    bool Output) override;

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry,
                   CSolver ***solver,
                   CConfig *config,
                   int val_iter,
                   bool val_update_geo) override;

  /*!
   * \brief Compute the spatial integration using a upwind scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Upwind_Residual(CGeometry *geometry,
                      CSolver **solver_container,
                      CNumerics **numerics_container,
                      CConfig *config,
                      unsigned short iMesh) override;

  /*!
   * \brief Compute the viscous residuals for the turbulent equation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Viscous_Residual(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics **numerics_container,
                        CConfig *config,
                        unsigned short iMesh,
                        unsigned short iRKStep) override;


  void Set_Heatflux_Areas(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Impose the Navier-Stokes boundary condition (strong).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry *geometry,
                          CSolver **solver_container,
                          CNumerics *conv_numerics,
                          CNumerics *visc_numerics,
                          CConfig *config,
                          unsigned short val_marker) override;

  /*!
   * \brief Impose a constant heat-flux condition at the wall.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFlux_Wall(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) override;

  /*!
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet(CGeometry *geometry,
                CSolver **solver_container,
                CNumerics *conv_numerics,
                CNumerics *visc_numerics,
                CConfig *config,
                unsigned short val_marker) override;
  /*!
   * \brief Impose the outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Outlet(CGeometry *geometry,
                 CSolver **solver_container,
                 CNumerics *conv_numerics,
                 CNumerics *visc_numerics,
                 CConfig *config,
                 unsigned short val_marker) override;

  /*!
   * \brief Impose the (received) conjugate heat variables.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ConjugateHeat_Interface(CGeometry *geometry,
                                  CSolver **solver_container,
                                  CNumerics *numerics,
                                  CConfig *config,
                                  unsigned short val_marker) override;

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker        - marker index
   * \param[in] val_vertex        - vertex index
   * \param[in] pos_var           - variable position (in vector of all conjugate heat variables)
   */
  inline su2double GetConjugateHeatVariable(unsigned short val_marker,
                                            unsigned long val_vertex,
                                            unsigned short pos_var) const override {
    return ConjugateVar[val_marker][val_vertex][pos_var];
  }

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker        - marker index
   * \param[in] val_vertex        - vertex index
   * \param[in] pos_var           - variable position (in vector of all conjugate heat variables)
   * \param[in] relaxation factor - relaxation factor for the change of the variables
   * \param[in] val_var           - value of the variable
   */
  inline void SetConjugateHeatVariable(unsigned short val_marker,
                                       unsigned long val_vertex,
                                       unsigned short pos_var,
                                       su2double relaxation_factor,
                                       su2double val_var) override {
    ConjugateVar[val_marker][val_vertex][pos_var] = relaxation_factor*val_var + (1.0-relaxation_factor)*ConjugateVar[val_marker][val_vertex][pos_var];
  }


  /*!
   * \brief Evaluate heat-flux related objectives.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Heat_Fluxes(CGeometry *geometry,
                   CSolver **solver_container,
                   CConfig *config) override;

  /*!
   * \brief Get value of the heat load (integrated heat flux).
   * \return Value of the heat load (integrated heat flux).
   */
  inline su2double GetTotal_HeatFlux() const override { return Total_HeatFlux; }

  /*!
   * \brief Get value of the integral-averaged temperature.
   * \return Value of the integral-averaged temperature.
   */
  inline su2double GetTotal_AvgTemperature() const override { return Total_AverageT; }

  /*!
   * \brief Compute objective output
   * \param[in] config - Definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   */
  void Evaluate_ObjFunc(const CConfig *config, CSolver**) override {
    const auto weight = config->GetWeight_ObjFunc(0);

    switch (config->GetKind_ObjFunc()) {
      case TOTAL_HEATFLUX:
        Total_ComboObj = weight * Total_HeatFlux;
        break;
      case AVG_TEMPERATURE:
        Total_ComboObj = weight * Total_AverageT;
        break;
      case CUSTOM_OBJFUNC:
        Total_ComboObj = weight * Total_Custom_ObjFunc;
        break;
      default:
        Total_ComboObj = 0.0;
    }
  }

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry,
                    CSolver **solver_container,
                    CConfig *config,
                    unsigned short iMesh,
                    unsigned long Iteration) override;

  /*!
   * \brief Set the initial condition for the FEM structural problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry **geometry,
                           CSolver ***solver_container,
                           CConfig *config,
                           unsigned long TimeIter) override;

  /*!
   * \brief Set the total residual adding the term that comes from the Dual Time-Stepping Strategy.
   * \param[in] geometry - Geometric definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void SetResidual_DualTime(CGeometry *geometry,
                            CSolver **solver_container,
                            CConfig *config,
                            unsigned short iRKStep,
                            unsigned short iMesh,
                            unsigned short RunTime_EqSystem) override;

  /*!
   * \brief Get the heat flux.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the heat flux.
   */
  inline su2double GetHeatFlux(unsigned short val_marker, unsigned long val_vertex) const override {
    return HeatFlux[val_marker][val_vertex];
  }
};
