/*!
 * \file CTransLMSolver.hpp
 * \brief Headers of the CTransLMSolver class
 * \author A. Aranake
 * \version 7.4.0 "Blackbird"
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

#include "CTurbSolver.hpp"

/*!
 * \class CTransLMSolver
 * \brief Main class for defining the transition model solver.
 * \ingroup Turbulence_Model
 * \author A. Aranake, S. Kang.
 */

class CTransLMSolver final : public CTurbSolver {
private:

  TURB_TRANS_CORRELATION TransCorrelation;
  TURB_FAMILY TurbFamily;

  su2double ReThetaC_Correlations(const su2double Tu, const su2double Re_Theta_t){

    su2double rethetac = 0.0;

      switch (TransCorrelation) {
      case TURB_TRANS_CORRELATION::MALAN: {
        rethetac = min(0.615 * Re_Theta_t + 61.5, Re_Theta_t);
        break;
      }

      case TURB_TRANS_CORRELATION::SULUKSNA: {
        rethetac = min(0.1 * exp(-0.0022 * Re_Theta_t + 12), 300.0);
        break;
      }

      case TURB_TRANS_CORRELATION::KRAUSE: {
        rethetac = 0.91 * Re_Theta_t + 5.32;
        break;
      }

      case TURB_TRANS_CORRELATION::KRAUSE_HYPER: {
        const su2double FirstTerm = -0.042 * pow(Tu, 3);
        const su2double SecondTerm = 0.4233 * pow(Tu, 2);
        const su2double ThirdTerm = 0.0118 * pow(Tu, 1);
        rethetac = Re_Theta_t / (FirstTerm + SecondTerm + ThirdTerm + 1.0744);
        break;
      }

      case TURB_TRANS_CORRELATION::MEDIDA_BAEDER: {
        const su2double FirstTerm = 4.45 * pow(Tu, 3);
        const su2double SecondTerm = 5.7 * pow(Tu, 2);
        const su2double ThirdTerm = 1.37 * pow(Tu, 1);
        rethetac = (FirstTerm - SecondTerm + ThirdTerm + 0.585) * Re_Theta_t;
        break;
      }

      case TURB_TRANS_CORRELATION::MEDIDA: {
        rethetac = 0.62 * Re_Theta_t;
        break;
      }

      case TURB_TRANS_CORRELATION::MENTER_LANGTRY: {

        if (Re_Theta_t <= 1870) {
          const su2double FirstTerm = (-396.035 * pow(10, -2));
          const su2double SecondTerm = (10120.656 * pow(10, -4)) * Re_Theta_t;
          const su2double ThirdTerm = (-868.230 * pow(10, -6)) * pow(Re_Theta_t, 2);
          const su2double ForthTerm = (696.506 * pow(10, -9)) * pow(Re_Theta_t, 3);
          const su2double FifthTerm = (-174.105 * pow(10, -12)) * pow(Re_Theta_t, 4);
          rethetac = FirstTerm + SecondTerm + ThirdTerm + ForthTerm + FifthTerm;
        } else {
          rethetac = Re_Theta_t - (593.11 + 0.482 * (Re_Theta_t - 1870.0));
        }

        break;
      }
      case TURB_TRANS_CORRELATION::DEFAULT: SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.", CURRENT_FUNCTION); break;

    }

    return rethetac;
  }

public:
  /*!
   * \overload
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
   * \brief Computes the effective intermtittency.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Postprocessing(CGeometry *geometry,
                      CSolver **solver_container,
                      CConfig *config,
                      unsigned short iMesh) override;

  /*!
   * \brief Compute the viscous flux for the LM equation at a particular edge.
   * \param[in] iEdge - Edge for which we want to compute the flux
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \note Calls a generic implementation after defining a SolverSpecificNumerics object.
   */
  void Viscous_Residual(unsigned long iEdge, CGeometry* geometry, CSolver** solver_container,
                        CNumerics* numerics, CConfig* config) override;

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
   * \brief Get the value of the intermittency.
   * \return Value of the turbulent kinetic energy.
   */
  inline su2double GetIntermittency_Inf(void) const override { return Solution_Inf[0]; }

  /*!
   * \brief Get the value of the intermittency.
   * \return Value of the turbulent kinetic energy.
   */
  inline su2double GetReThetaT_Inf(void) const override { return Solution_Inf[1]; }

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry** geometry, CSolver*** solver, CConfig* config, int val_iter, bool val_update_geo) final;

};
