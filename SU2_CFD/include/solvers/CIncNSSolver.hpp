/*!
 * \file CIncNSSolver.hpp
 * \brief Headers of the CIncNSSolver class
 * \author F. Palacios, T. Economon, T. Albring
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

#include "CIncEulerSolver.hpp"

/*!
 * \class CIncNSSolver
 * \brief Main class for defining the incompressible Navier-Stokes flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios, T. Economon, T. Albring
 */
class CIncNSSolver final : public CIncEulerSolver {

  /*!
   * \brief Generic implementation of the isothermal, heatflux and heat-transfer/convection walls.
   */
  void BC_Wall_Generic(const CGeometry *geometry,
                       const CConfig *config,
                       unsigned short val_marker,
                       unsigned short kind_boundary);

  /*!
   * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver **solver_container,
                                       const CConfig *config) override;

  /*!
   * \brief Compute the viscous contribution for a particular edge.
   * \param[in] iEdge - Edge for which the flux and Jacobians are to be computed.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void Viscous_Residual(unsigned long iEdge, CGeometry *geometry, CSolver **solver_container,
                        CNumerics *numerics, CConfig *config) override;

/*!
   * \brief Computes the wall shear stress (Tau_Wall) on the surface using a wall function.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetTau_Wall_WF(CGeometry *geometry,
                     CSolver** solver_container,
                     const CConfig* config);

  /*!
   * \brief Compute necessary quantities (massflow, integrated heatflux, avg density)
   *        for streamwise periodic cases. Also sets new delta P for prescribed massflow.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   */
  void GetStreamwise_Periodic_Properties(const CGeometry *geometry,
                                         CConfig *config,
                                         const unsigned short iMesh);

  /*!
   * \brief Compute recovered pressure/temperature for streamwise periodic flow and store in CVariable.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   */
  void Compute_Streamwise_Periodic_Recovered_Values(CConfig *config, const CGeometry *geometry,
                                                    const unsigned short iMesh);

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CIncNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Restart residual and compute gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
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
   * \brief Impose a no-slip condition.
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
   * \brief Impose an isothermal temperature condition at the wall.
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
   * \brief Impose a heat flux by prescribing a heat transfer coefficient and a temperature at infinity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatTransfer_Wall(const CGeometry *geometry,
                            const CConfig *config,
                            const unsigned short val_marker) override;

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

};
