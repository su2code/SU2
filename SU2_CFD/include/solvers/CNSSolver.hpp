/*!
 * \file CNSSolver.hpp
 * \brief Headers of the CNSSolver class
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

#include "CEulerSolver.hpp"

/*!
 * \class CNSSolver
 * \brief Main class for defining the Navier-Stokes flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios
 */
class CNSSolver final : public CEulerSolver {
private:

  vector<su2double> Surface_Buffet_Metric;  /*!< \brief Integrated separation sensor for each monitoring surface. */
  vector<su2double> Buffet_Metric;          /*!< \brief Integrated separation sensor for each boundary. */
  vector<vector<su2double> > Buffet_Sensor; /*!< \brief Separation sensor for each boundary and vertex. */
  su2double Total_Buffet_Metric = 0.0;      /*!< \brief Integrated separation sensor for all the boundaries. */

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   */
  void SetRoe_Dissipation(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver **solver_container,
                                       const CConfig *config) override;

  /*!
   * \brief Common code for wall boundaries, add the residual and Jacobian
   * contributions due to grid motion associated with a particular boundary point.
   */
  void AddDynamicGridResidualContribution(unsigned long iPoint,
                                          unsigned long Point_Normal,
                                          const CGeometry* geometry,
                                          const su2double* UnitNormal,
                                          su2double Area,
                                          const su2double* GridVel,
                                          su2double** Jacobian_i,
                                          su2double& Res_Conv,
                                          su2double& Res_Visc) const;

  /*!
   * \brief Get the wall temperature at a given vertex of a given marker for CHT problems.
   */
  su2double GetCHTWallTemperature(const CConfig* config,
                                  unsigned short val_marker,
                                  unsigned long iVertex,
                                  su2double thermal_conductivity,
                                  su2double dist_ij,
                                  su2double There,
                                  su2double Temperature_Ref) const;

  /*!
   * \brief Generic implementation of the isothermal wall also covering CHT cases,
   * for which the wall temperature is given by GetCHTWallTemperature.
   */
  void BC_Isothermal_Wall_Generic(CGeometry *geometry,
                                  CSolver **solver_container,
                                  CNumerics *conv_numerics,
                                  CNumerics *visc_numerics,
                                  CConfig *config,
                                  unsigned short val_marker,
                                  bool cht_mode = false);

  /*!
   * \brief Generic implementation of the heatflux and heat-transfer/convection walls.
   */
  void BC_HeatFlux_Wall_Generic(const CGeometry *geometry,
                                const CConfig *config,
                                unsigned short val_marker,
                                unsigned short kind_boundary);

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
  void SetTau_Wall_WF(CGeometry *geometry, CSolver** solver_container, const CConfig* config);

public:
  /*!
   * \brief Constructor of the class.
   */
  CNSSolver() = default;

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CNSSolver() = default;

  /*!
   * \brief Provide the buffet metric.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the buffet metric on the surface <i>val_marker</i>.
   */
  inline su2double GetSurface_Buffet_Metric(unsigned short val_marker) const override { return Surface_Buffet_Metric[val_marker]; }

  /*!
   * \brief Get the buffet metric.
   * \return Value of the buffet metric.
   */
  inline su2double GetTotal_Buffet_Metric() const override { return Total_Buffet_Metric; }

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
   * \brief Compute weighted-sum "combo" objective output
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Container vector with all the solutions.
   */
  void Evaluate_ObjFunc(const CConfig *config, CSolver **solver) override;

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
   * \brief Impose a heat flux by prescribing a heat transfer coefficient and a temperature at infinity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatTransfer_Wall(const CGeometry *geometry,
                            const CConfig *config,
                            const unsigned short val_marker) override;

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
   * \brief Impose the Navier-Stokes boundary condition (strong) with values from a CHT coupling.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ConjugateHeat_Interface(CGeometry *geometry,
                                  CSolver **solver_container,
                                  CNumerics *numerics,
                                  CConfig *config,
                                  unsigned short val_marker) override;

  /*!
   * \brief Compute the buffet sensor.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Buffet_Monitoring(const CGeometry *geometry, const CConfig *config) override;

  /*!
   * \brief Get the value of the buffet sensor
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the buffet sensor.
   */
  inline su2double GetBuffetSensor(unsigned short val_marker, unsigned long val_vertex) const override {
    return Buffet_Sensor[val_marker][val_vertex];
  }
};
