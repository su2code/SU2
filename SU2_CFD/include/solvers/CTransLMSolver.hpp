/*!
 * \file CTransLMSolver.hpp
 * \brief Headers of the CTransLMSolver class
 * \author A. Aranake, M. Cerabona
 * \version 7.3.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/solvers/CTransSolver.hpp"
//#include "../../include/solvers/CTurbSolver.hpp"
#include "../../include/variables/CTransLMVariable.hpp"

/*!
 * \class CTransLMSolver
 * \brief Main class for defining the transition model solver.
 * \ingroup Transition_Model
 * \author A. Aranake, M. Cerabona
 */
class CTransLMSolver final : public CTransSolver {
private:

   su2double constants[9] = {0.0}; /*!< \brief Constants for the model. */

    // non credo che la parte SetTurbVars_WF serva -> tolta anche dal .cpp

public:
    /*!
   * \brief Constructor.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
    CTransLMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

    /*!
   * \brief Destructor of the class.
   */
    ~CTransLMSolver() = default;

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
   * \brief Computes the eddy viscosity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
    // non so se è necessario rifarlo dopo il SST
    void Postprocessing(CGeometry *geometry,
                        CSolver **solver_container,
                        CConfig *config,
                        unsigned short iMesh) override;

    /*!
    * \brief Compute the viscous flux for the transition equation at a particular edge.
    * \param[in] iEdge - Edge for which we want to compute the flux
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \note Calls a generic implementation after defining a SolverSpecificNumerics object.
    */
    void Viscous_Residual(unsigned long iEdge, CGeometry *geometry, CSolver **solver_container,
                          CNumerics *numerics, CConfig *config) override;

    /*!
     * \brief Source term computation.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] numerics_container - Description of the numerical method.
     * \param[in] config - Definition of the particular problem.
     * \param[in] iMesh - Index of the mesh in multigrid computations.
     */
    void Source_Residual(CGeometry *geometry,
                         CSolver **solver_container,
                         CNumerics **numerics_container,
                         CConfig *config,
                         unsigned short iMesh) override;

    /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
    void Source_Template(CGeometry *geometry,
                         CSolver **solver_container,
                         CNumerics *numerics,
                         CConfig *config,
                         unsigned short iMesh) override;

    /*!
   * \brief Impose the Langtry Menter transition wall boundary condition.
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
   * \brief Impose the Navier-Stokes wall boundary condition.
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
   * \brief Impose the Far Field boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
     */
    void BC_Far_Field(CGeometry *geometry,
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
    * \brief Get the constants for the LM model.
    * \return A pointer to an array containing a set of constants
    */
    inline const su2double *GetConstants() const override { return constants; }


    // in queste due sotto attenzione al contenitore Solution_Inf ->
    // va già bene così?
    /*!
    * \brief Get the value of the intermittency.
    * \return Value of the intermittency.
    */
    inline su2double GetIntermittency_Inf(void) const override { return Solution_Inf[0]; }

    /*!
     * \brief Get the value of the Re_theta_t.
     * \return Value of the Re_theta_t.
     */
    inline su2double GetRethetat_Inf(void) const override { return Solution_Inf[1]; }



};
