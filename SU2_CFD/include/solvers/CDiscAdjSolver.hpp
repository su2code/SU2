/*!
 * \file CDiscAdjSolver.hpp
 * \brief Headers of the CDiscAdjSolver class
 * \author T. Albring
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

#include "CSolver.hpp"
#include "../variables/CDiscAdjVariable.hpp"

/*!
 * \class CDiscAdjSolver
 * \ingroup DiscAdj
 * \brief Main class for defining the discrete adjoint solver.
 * \author T. Albring
 */
class CDiscAdjSolver final : public CSolver {
protected:
  static constexpr size_t MAXNDIM = 3;  /*!< \brief Max number of space dimensions, used in some static arrays. */
  static constexpr size_t MAXNVAR = 32; /*!< \brief Max number of variables, for static arrays. */

  static constexpr size_t OMP_MAX_SIZE = 1024; /*!< \brief Max chunk size for light point loops. */

  unsigned long omp_chunk_size; /*!< \brief Chunk size used in light point loops. */

  unsigned short KindDirect_Solver;
  CSolver *direct_solver;
  vector<vector<su2double> > CSensitivity; /*!< \brief Shape sensitivity coefficient for each boundary and vertex. */
  vector<su2double> Sens_Geo;    /*!< \brief Total shape sensitivity for each monitored boundary. */
  su2double Total_Sens_Mach;     /*!< \brief Total mach sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_AoA;      /*!< \brief Total angle of attack sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_Geo;      /*!< \brief Total shape sensitivity coefficient for all the boundaries. */
  su2double Total_Sens_Press;    /*!< \brief Total farfield sensitivity to pressure. */
  su2double Total_Sens_Temp;     /*!< \brief Total farfield sensitivity to temperature. */
  su2double Total_Sens_BPress;   /*!< \brief Total sensitivity to outlet pressure. */
  su2double Total_Sens_Density;  /*!< \brief Total sensitivity to initial density (incompressible). */
  su2double Total_Sens_ModVel;   /*!< \brief Total sensitivity to inlet velocity (incompressible). */
  su2double Mach, Alpha, Beta, Pressure, Temperature, BPressure, ModVel;
  su2double TemperatureRad, Total_Sens_Temp_Rad;

  CDiscAdjVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

public:

  /*!
   * \brief Constructor of the class.
   */
  CDiscAdjSolver() = default;

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Initialize the discrete adjoint solver with the corresponding direct solver.
   * \param[in] Kind_Solver - The kind of direct solver.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CDiscAdjSolver(CGeometry *geometry, CConfig *config, CSolver* solver, unsigned short Kind_Solver, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjSolver() override;

  /*!
   * \brief Performs the preprocessing of the adjoint AD-based solver.
   *        Registers all necessary variables on the tape. Called while tape is active.
   * \param[in] geometry_container - The geometry container holding all grid levels.
   * \param[in] config_container - The particular config.
   */
  void RegisterSolution(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Performs the preprocessing of the adjoint AD-based solver.
   *        Registers all necessary variables that are output variables on the tape.
   *        Called while tape is active.
   * \param[in] geometry_container - The geometry container holding all grid levels.
   * \param[in] config_container - The particular config.
   */
  void RegisterOutput(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Sets the adjoint values of the output of the flow (+turb.) iteration
   *        before evaluation of the tape.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] config - The particular config.
   */
  void SetAdjoint_Output(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Sets the adjoint values of the input variables of the flow (+turb.) iteration
   *        after tape has been evaluated.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] config - The particular config.
   * \param[in] CrossTerm - Boolean for CrossTerm.
   */
  void ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config, bool CrossTerm) override;

  /*!
   * \brief Set the surface sensitivity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSurface_Sensitivity(CGeometry *geometry, CConfig* config) override;

  /*!
   * \brief Extract and set the geometrical sensitivity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetSensitivity(CGeometry *geometry, CConfig *config, CSolver*) override;

  /*!
   * \brief Provide the total shape sensitivity coefficient.
   * \return Value of the geometrical sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Sens_Geo() const override { return Total_Sens_Geo; }

  /*!
   * \brief Set the total Mach number sensitivity coefficient.
   * \return Value of the Mach sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Sens_Mach() const override { return Total_Sens_Mach; }

  /*!
   * \brief Set the total angle of attack sensitivity coefficient.
   * \return Value of the angle of attack sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Sens_AoA() const override { return Total_Sens_AoA; }

  /*!
   * \brief Set the total farfield pressure sensitivity coefficient.
   * \return Value of the farfield pressure sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Sens_Press() const override { return Total_Sens_Press; }

  /*!
   * \brief Set the total farfield temperature sensitivity coefficient.
   * \return Value of the farfield temperature sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Sens_Temp() const override { return Total_Sens_Temp; }

  /*!
   * \author H. Kline
   * \brief Get the total Back pressure number sensitivity coefficient.
   * \return Value of the Back sensitivity coefficient
   *         (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Sens_BPress() const override { return Total_Sens_BPress; }

  /*!
   * \brief Get the total density sensitivity coefficient.
   * \return Value of the density sensitivity.
   */
  inline su2double GetTotal_Sens_Density() const override { return Total_Sens_Density; }

  /*!
   * \brief Get the total velocity magnitude sensitivity coefficient.
   * \return Value of the velocity magnitude sensitivity.
   */
  inline su2double GetTotal_Sens_ModVel() const override { return Total_Sens_ModVel; }

  /*!
   * \brief Get the shape sensitivity coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the sensitivity coefficient.
   */
  inline su2double GetCSensitivity(unsigned short val_marker,
                                   unsigned long val_vertex) const override {
    return CSensitivity[val_marker][val_vertex];
  }

  /*!
   * \brief Prepare the solver for a new recording.
   * \param[in] kind_recording - Kind of AD recording.
   */
  void SetRecording(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reset - If true reset variables to their initial values.
   */
  void RegisterVariables(CGeometry *geometry,
                         CConfig *config,
                         bool reset = false) override;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Update the dual-time derivatives.
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
   * \brief Depends on the direct solver.
   */
  inline bool GetHasHybridParallel() const override {
    if (direct_solver) return direct_solver->GetHasHybridParallel();
    return false;
  }

};
