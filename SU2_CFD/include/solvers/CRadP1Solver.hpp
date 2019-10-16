/*!
 * \file CRadP1Solver.hpp
 * \brief Declaration and inlines of the class to compute a P1 radiation problem.
 * \author Ruben Sanchez
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "CRadSolver.hpp"

class CRadP1Solver final: public CRadSolver {
private:

protected:

  su2double Temperature_Inf;      /*!< \brief Temperature at the infinity. */

public:
  /*!
   * \brief Constructor of the class.
   */
  CRadP1Solver(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] FluidModel
   */
  CRadP1Solver(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CRadP1Solver(void);

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
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                      unsigned short iMesh);


  /*!
   * \brief Update the solution using an implicit solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config);

  /*!
   * \brief Compute the viscous residuals for the turbulent equation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                        CConfig *config, unsigned short iMesh, unsigned short iRKStep);

  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] second_numerics - Description of the second numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                       CConfig *config, unsigned short iMesh);

  /*!
   * \brief Impose the Marshak BC on Heatflux wall boundaries.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief Impose an isothermal wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                          CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief Impose the Far Field boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief Impose the Marshak BC on inlet boundaries.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                       CNumerics *visc_numerics, CConfig *config, unsigned short val_marker){
    BC_Marshak(geometry, solver_container, config, val_marker);
  }


  /*!
   * \brief Impose the Marshak BC on outlet boundaries.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker){
    BC_Marshak(geometry, solver_container, config, val_marker);
  }

  /*!
   * \brief Impose via the Marshak BC on Euler walls.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                            CNumerics *visc_numerics, CConfig *config, unsigned short val_marker){
    BC_Marshak(geometry, solver_container, config, val_marker);
  }

  /*!
  * \brief Impose the Marshak BC on sliding interfaces.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] solver_container - Container vector with all the solutions.
  * \param[in] conv_numerics - Description of the numerical method.
  * \param[in] visc_numerics - Description of the numerical method.
  * \param[in] config - Definition of the particular problem.
  */
  inline void BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config){
    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
        BC_Marshak(geometry, solver_container, config, iMarker);
      }
    }
  }

  /*!
   * \brief Impose the Marshak BC on conjugate heat interfaces.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline void BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container,
                                         CNumerics *numerics, CConfig *config, unsigned short val_marker){
    BC_Marshak(geometry, solver_container, config, val_marker);
  }

  /*!
   * \brief Impose the Marshak boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Marshak(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_marker);

  /*!
   * \brief Compute the time step for solving the Radiation P1 equation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                    unsigned short iMesh, unsigned long Iteration);

  /*!
   * \brief Set the freestream temperature.
   * \param[in] Value of freestream temperature.
   */
  inline void SetTemperature_Inf(su2double t_inf) {Temperature_Inf = t_inf;}


  /*!
   * \brief Get the temperature value at infinity.
   * \return Value of the temperature at infinity.
   */
  inline su2double GetTemperature_Inf(void) { return Temperature_Inf; }

};
