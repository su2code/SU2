/*!
 * \file CEulerSolver.hpp
 * \brief Headers of the CEulerSolver class
 * \author F. Palacios, T. Economon
 * \version 7.0.3 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "CSolver.hpp"
#include "CNEMOEulerSolver.hpp"
#include "../variables/CNEMONSVariable.hpp"

/*!
 * \class CNEMONSSolver
 * \brief Main class for defining the NEMO Navier-Stokes flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author S. R. Copeland, F. Palacios, W. Maier.
 * \version 6.1
 */
class CNEMONSSolver : public CNEMOEulerSolver {
private:

  su2double Viscosity_Inf; /*!< \brief Viscosity at the infinity. */
  su2double Tke_Inf;       /*!< \brief Turbulent kinetic energy at infinity. */
  su2double Prandtl_Lam,   /*!< \brief Laminar Prandtl number. */
  Prandtl_Turb;            /*!< \brief Turbulent Prandtl number. */
  su2double *CD_Visc,	   /*!< \brief Drag coefficient (viscous contribution) for each boundary. */
  *CSF_Visc,               /*!< \brief Side force coefficient (viscous contribution) for each boundary. */
  *CL_Visc,		           /*!< \brief Lift coefficient (viscous contribution) for each boundary. */
  *CMx_Visc,			   /*!< \brief Moment x coefficient (viscous contribution) for each boundary. */
  *CMy_Visc,			   /*!< \brief Moment y coefficient (viscous contribution) for each boundary. */
  *CMz_Visc,			   /*!< \brief Moment z coefficient (viscous contribution) for each boundary. */
  *CoPx_Visc,              /*!< \brief Moment x coefficient (viscous contribution) for each boundary. */
  *CoPy_Visc,              /*!< \brief Moment y coefficient (viscous contribution) for each boundary. */
  *CoPz_Visc,              /*!< \brief Moment z coefficient (viscous contribution) for each boundary. */
  *CFx_Visc,			   /*!< \brief Force x coefficient (viscous contribution) for each boundary. */
  *CFy_Visc,			   /*!< \brief Force y coefficient (viscous contribution) for each boundary. */
  *CFz_Visc,			   /*!< \brief Force z coefficient (viscous contribution) for each boundary. */
  *Surface_CL_Visc,        /*!< \brief Lift coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CD_Visc,        /*!< \brief Drag coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CSF_Visc,       /*!< \brief Side-force coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CEff_Visc,      /*!< \brief Side-force coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFx_Visc,       /*!< \brief Force x coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFy_Visc,       /*!< \brief Force y coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CFz_Visc,       /*!< \brief Force z coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMx_Visc,       /*!< \brief Moment x coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMy_Visc,       /*!< \brief Moment y coefficient (viscous contribution) for each monitoring surface. */
  *Surface_CMz_Visc,       /*!< \brief Moment z coefficient (viscous contribution) for each monitoring surface. */
  *Surface_Buffet_Metric,  /*!< \brief Integrated separation sensor for each monitoring surface. */
  *CEff_Visc,              /*!< \brief Efficiency (Cl/Cd) (Viscous contribution) for each boundary. */
  *CMerit_Visc,            /*!< \brief Rotor Figure of Merit (Viscous contribution) for each boundary. */
  *Buffet_Metric,          /*!< \brief Integrated separation sensor for each boundary. */
  *CT_Visc,                /*!< \brief Thrust coefficient (viscous contribution) for each boundary. */
  *CQ_Visc,                /*!< \brief Torque coefficient (viscous contribution) for each boundary. */
  *HF_Visc,                /*!< \brief Heat load (viscous contribution) for each boundary. */
  *MaxHF_Visc,             /*!< \brief Maximum heat flux (viscous contribution) for each boundary. */
  ***HeatConjugateVar,     /*!< \brief Conjugate heat transfer variables for each boundary and vertex. */
  ***CSkinFriction,        /*!< \brief Skin friction coefficient for each boundary and vertex. */
  **Buffet_Sensor;         /*!< \brief Separation sensor for each boundary and vertex. */
  su2double Total_Buffet_Metric;  /*!< \brief Integrated separation sensor for all the boundaries. */
  su2double *ForceViscous,        /*!< \brief Viscous force for each boundary. */
  *MomentViscous;                 /*!< \brief Inviscid moment for each boundary. */
  su2double AllBound_CD_Visc,     /*!< \brief Drag coefficient (viscous contribution) for all the boundaries. */
  AllBound_CL_Visc,              /*!< \brief Lift coefficient (viscous contribution) for all the boundaries. */
  AllBound_CSF_Visc,      /*!< \brief Sideforce coefficient (viscous contribution) for all the boundaries. */
  AllBound_CMx_Visc,      /*!< \brief Moment x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMy_Visc,      /*!< \brief Moment y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMz_Visc,      /*!< \brief Moment z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPx_Visc,     /*!< \brief Moment x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPy_Visc,     /*!< \brief Moment y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CoPz_Visc,     /*!< \brief Moment z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CEff_Visc,     /*!< \brief Efficient coefficient (Viscous contribution) for all the boundaries. */
  AllBound_CFx_Visc,      /*!< \brief Force x coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFy_Visc,      /*!< \brief Force y coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CFz_Visc,      /*!< \brief Force z coefficient (inviscid contribution) for all the boundaries. */
  AllBound_CMerit_Visc,   /*!< \brief Rotor Figure of Merit coefficient (Viscous contribution) for all the boundaries. */
  AllBound_CT_Visc,       /*!< \brief Thrust coefficient (viscous contribution) for all the boundaries. */
  AllBound_CQ_Visc,       /*!< \brief Torque coefficient (viscous contribution) for all the boundaries. */
  AllBound_HF_Visc,       /*!< \brief Heat load (viscous contribution) for all the boundaries. */
  AllBound_MaxHF_Visc;    /*!< \brief Maximum heat flux (viscous contribution) for all boundaries. */
  su2double StrainMag_Max, Omega_Max; /*!< \brief Maximum Strain Rate magnitude and Omega. */
  su2double *primitives_aux; /*!< \brief Primitive auxiliary variables (Y_s, T, Tve, ...) in compressible flows. */

 

public:

  /*!
   * \brief Constructor of the class.
   */
  CNEMONSSolver(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CNEMONSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CNEMONSSolver(void);

  /*!
   * \brief Compute the viscosity at the infinity.
   * \return Value of the viscosity at the infinity.
   */
  inline su2double GetViscosity_Inf(void) const override { return Viscosity_Inf; }

  /*!
   * \brief Restart residual and compute gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);

  /*!
   * \brief Compute the gradient of the primitive variables using Green-Gauss method,
   *        and stores the result in the <i>Gradient_Primitive</i> variable.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reconstruction - indicator that the gradient being computed is for upwind reconstruction.
   */
  void SetPrimitive_Gradient_GG(CGeometry *geometry,
                                CConfig *config,
                                bool reconstruction = false) final;

  /*!
   * \brief Compute the gradient of the primitive variables using a Least-Squares method,
   *        and stores the result in the <i>Gradient_Primitive</i> variable.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reconstruction - indicator that the gradient being computed is for upwind reconstruction.
   */
  void SetPrimitive_Gradient_LS(CGeometry *geometry,
                                CConfig *config,
                                bool reconstruction = false) final;

  /*!
   * \brief Impose the symmetry boundary condition using the residual.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

//   /*!
//   * \brief A virtual member.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] config - Definition of the particular problem.
//   */
//  void SetPrimitive_Gradient_GG(CGeometry *geometry,
//                                CConfig *config);
//
//  /*!
//   * \brief A virtual member.
//   * \param[in] geometry - Geometrical definition of the problem.
//   * \param[in] config - Definition of the particular problem.
//   */
//  void SetPrimitive_Gradient_LS(CGeometry *geometry,
//                                CConfig *config);

  /*!
   * \brief Impose a constant heat-flux condition at the wall.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief Impose a constant heat-flux condition at the wall.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFluxCatalytic_Wall(CGeometry *geometry,
                                 CSolver **solver_container,
                                 CNumerics *conv_numerics,
                                 CNumerics *visc_numerics,
                                 CConfig *config, unsigned short val_marker);

  /*!
   * \brief Impose a constant heat-flux condition at the wall.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFluxNonCatalytic_Wall(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CNumerics *conv_numerics,
                                    CNumerics *visc_numerics,
                                    CConfig *config, unsigned short val_marker);

  /*!
   * \brief Impose the Navier-Stokes boundary condition (strong).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                          CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  /*!
   * \brief Impose the Navier-Stokes boundary condition (strong).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_IsothermalCatalytic_Wall(CGeometry *geometry,
                                   CSolver **solver_container,
                                   CNumerics *conv_numerics,
                                   CNumerics *visc_numerics,
                                   CConfig *config, unsigned short val_marker);

  /*!
   * \brief Impose the Navier-Stokes boundary condition (strong).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method for convective terms.
   * \param[in] visc_numerics - Description of the numerical method for viscous terms.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_IsothermalNonCatalytic_Wall(CGeometry *geometry,
                                      CSolver **solver_container,
                                      CNumerics *conv_numerics,
                                      CNumerics *visc_numerics,
                                      CConfig *config,
                                      unsigned short val_marker);

  /*!
   * \brief Compute the viscous forces and all the addimensional coefficients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Friction_Forces(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Get the non dimensional lift coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the lift coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCL_Visc(unsigned short val_marker) const override { return CL_Visc[val_marker]; }
  
  /*!
   * \brief Get the non dimensional drag coefficient (viscous contribution).
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \return Value of the drag coefficient (viscous contribution) on the surface <i>val_marker</i>.
   */
  inline su2double GetCD_Visc(unsigned short val_marker) const override { return CD_Visc[val_marker]; }

  /*!
   * \brief Get the inviscid contribution to the lift coefficient.
   * \return Value of the lift coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CL_Visc() const override { return AllBound_CL_Visc; }

  /*!
   * \brief Get the inviscid contribution to the drag coefficient.
   * \return Value of the drag coefficient (inviscid contribution).
   */
  inline su2double GetAllBound_CD_Visc() const override { return AllBound_CD_Visc; }


  /*!
   * \brief Compute the viscous residuals.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                        CConfig *config, unsigned short iMesh, unsigned short iRKStep);

  /*!
   * \brief Compute the viscous residuals.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
 //void Source_Residual(CGeometry *geometry, CSolver **solution_container, CNumerics **numerics_container,
 //                     CNumerics *second_solver, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Get the skin friction coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the skin friction coefficient.
   */
    inline su2double GetCSkinFriction(unsigned short val_marker,
                                    unsigned long val_vertex,
                                    unsigned short val_dim) const override {
    return CSkinFriction[val_marker][val_dim][val_vertex];
  }

  /*!
   * \brief Get the skin friction coefficient.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the heat transfer coefficient.
   */
  inline su2double GetHeatFlux(unsigned short val_marker, unsigned long val_vertex) const override {
    return HeatFlux[val_marker][val_vertex];
  }


  /*!
   * \brief Get the y plus.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the y plus.
   */
  inline su2double GetYPlus(unsigned short val_marker, unsigned long val_vertex) const override {
    return YPlus[val_marker][val_vertex];
  }

 

};