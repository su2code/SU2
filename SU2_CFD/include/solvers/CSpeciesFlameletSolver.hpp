/*!
 * \file CSpeciesFlameletSolver.hpp
 * \brief Headers of the CSpeciesFlameletSolver class
 * \author D. Mayer, N. Beishuizen, T. Economon
 * \version 7.5.1 "Blackbird"
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

#include "CSpeciesSolver.hpp"

/*!
 * \class CSpeciesFlameletSolver
 * \brief Main class for defining the flamelet model solver.
 * \author N. Beishuizen
 * \ingroup Scalar_Transport
 */
class CSpeciesFlameletSolver final : public CSpeciesSolver {
 private:
  vector<su2activematrix> conjugate_var; /*!< \brief CHT variables for each boundary and vertex. */
  bool include_mixture_fraction = false; /*!< \brief include mixture fraction as a controlling variable. */
  /*!
   * \brief Compute the preconditioner for low-Mach flows.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPreconditioner(CGeometry* geometry, CSolver** solver_container, CConfig* config);

  /*!
   * \brief Compute the primitive variables (diffusivities).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - Boolean to determine whether to print output.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver** solver_container, CConfig* config, bool Output);

  /*!
   * \brief Generic implementation of the isothermal wall also covering CHT cases,
   * for which the wall temperature is given by GetConjugateHeatVariable.
   */
  void BC_Isothermal_Wall_Generic(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                  CNumerics* visc_numerics, CConfig* config, unsigned short val_marker,
                                  bool cht_mode = false);

  /*!
   * \brief Reverse look-up to retrieve enthalpy value based on temperature and other controlling variables.
   * \param[in] fluid_model - pointer to flamelet fluid model.
   * \param[in] val_temp - temperature value used for reverse look-up.
   * \param[in] val_enth - pointer to enthalpy value to be retrieved.
   * \param[in] scalar_solution - local scalar solution.
   * \param[out] Converged - 0 if Newton solver converged, 1 if not.
   */
  unsigned long GetEnthFromTemp(CFluidModel * fluid_model, su2double const val_temp, su2double * val_enth, const su2double * scalar_solution);

  /*!
   * \brief Find maximum progress variable value within the manifold for the current solution.
   * \param[in] fluid_model - pointer to flamelet fluid model.
   * \param[in] scalars - local scalar solution.
   * \return - maximum progress variable value within manifold bounds.
   */
  su2double GetBurntProgressVariable(CFluidModel * fluid_model, const su2double * scalars);

  /*!
   * \brief Retrieve scalar source terms from manifold.
   * \param[in] config - definition of particular problem.
   * \param[in] fluid_model_local - pointer to flamelet fluid model.
   * \param[in] iPoint - node ID.
   * \param[in] scalars - local scalar solution.
   * \param[in] table_source_names - variable names of scalar source terms.
   * \return - within manifold bounds (0) or outside manifold bounds (1).
   */
  unsigned long SetScalarSources(CConfig *config, CFluidModel *fluid_model_local, unsigned long iPoint, vector<su2double> &scalars);

  /*!
   * \brief Retrieve passive look-up data from manifold.
   * \param[in] config - definition of particular problem.
   * \param[in] fluid_model_local - pointer to flamelet fluid model.
   * \param[in] iPoint - node ID.
   * \param[in] scalars - local scalar solution.
   * \param[in] table_lookup_names - variable names of scalar source terms.
   * \return - within manifold bounds (0) or outside manifold bounds (1).
   */
  unsigned long SetScalarLookUps(CConfig *config, CFluidModel *fluid_model_local, unsigned long iPoint, vector<su2double> &scalars);

 public:
  /*!
   * \brief Constructor.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] FluidModel
   */
  CSpeciesFlameletSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh);

  /*!
   * \brief Restart residual and compute gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - Boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh,
                     unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) override;

  /*!
   * \brief Set the initial condition for the scalar transport problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry** geometry, CSolver*** solver_container, CConfig* config,
                           unsigned long ExtIter) override;

  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics** numerics_container, CConfig* config,
                       unsigned short iMesh) override;

  /*!
   * \brief Impose the Navier-Stokes wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                          CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Impose the (received) conjugate heat variables.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ConjugateHeat_Interface(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics, CConfig* config,
                                  unsigned short val_marker) override;

  /*!
   * \brief Get the conjugate heat variables.
   * \param[in] val_marker - The marker index.
   * \param[in] val_vertex - The vertex index.
   * \param[in] pos_var - The variable position (in vector of all conjugate heat variables).
   */
  inline su2double GetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex,
                                            unsigned short pos_var) const override {
    return conjugate_var[val_marker][val_vertex][pos_var];
  }

  /*!
   * \brief Set the conjugate heat variables.
   * \param[in] val_marker - The marker index.
   * \param[in] val_vertex - The vertex index.
   * \param[in] pos_var - The variable position (in vector of all conjugate heat variables).
   * \param[in] relaxation_factor - The relaxation factor for the change of the variables.
   * \param[in] val_var - The value of the variable.
   */
  inline void SetConjugateHeatVariable(unsigned short val_marker, unsigned long val_vertex, unsigned short pos_var,
                                       su2double relaxation_factor, su2double val_var) override {
    conjugate_var[val_marker][val_vertex][pos_var] =
        relaxation_factor * val_var + (1.0 - relaxation_factor) * conjugate_var[val_marker][val_vertex][pos_var];
  }
};
