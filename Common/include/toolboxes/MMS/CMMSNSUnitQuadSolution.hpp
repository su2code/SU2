/*!
 * \file CMMSNSUnitQuadSolution.hpp
 * \brief Header file for the class CMMSNSUnitQuadSolution.
 *        The implementations are in the <i>CMMSNSUnitQuadSolution.cpp</i> file.
 * \author T. Economon, E. van der Weide
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

#include <cmath>
#include "CVerificationSolution.hpp"

/*!
 * \class CMMSNSUnitQuadSolution
 * \brief Class to define the required data for the manufactured solution of the
          laminar Navier-Stokes equations on a unit quad.
 * \author E. van der Weide, T. Economon
 */
class CMMSNSUnitQuadSolution final : public CVerificationSolution {
 protected:
  /*--- Variables that define the solution and MMS source term. ---*/
  su2double Gamma;        /*!< \brief Specific heat ratio. */
  su2double RGas;         /*!< \brief Gas constant. */
  su2double Viscosity;    /*!< \brief Viscosity, must be constant. */
  su2double Conductivity; /*!< \brief Thermal conductivity, must be constant. */

  /*--- Constants, which describe this manufactured solution. This is a viscous
        solution on the unit quad, where the primitive variables vary as a
        combination of sine and cosine functions. The unit quad is probably not
        necessary, and an arbitrary domain should work as well. ---*/

  su2double L;       /*!< \brief Length scale. */
  su2double a_Px;    /*!< \brief Parameter for the pressure solution. */
  su2double a_Pxy;   /*!< \brief Parameter for the pressure solution. */
  su2double a_Py;    /*!< \brief Parameter for the pressure solution. */
  su2double a_rhox;  /*!< \brief Parameter for the density solution. */
  su2double a_rhoxy; /*!< \brief Parameter for the density solution. */
  su2double a_rhoy;  /*!< \brief Parameter for the density solution. */
  su2double a_ux;    /*!< \brief Parameter for the x-velocity solution. */
  su2double a_uxy;   /*!< \brief Parameter for the x-velocity solution. */
  su2double a_uy;    /*!< \brief Parameter for the x-velocity solution. */
  su2double a_vx;    /*!< \brief Parameter for the y-velocity solution. */
  su2double a_vxy;   /*!< \brief Parameter for the y-velocity solution. */
  su2double a_vy;    /*!< \brief Parameter for the y-velocity solution. */
  su2double P_0;     /*!< \brief Parameter for the pressure solution. */
  su2double P_x;     /*!< \brief Parameter for the pressure solution. */
  su2double P_xy;    /*!< \brief Parameter for the pressure solution. */
  su2double P_y;     /*!< \brief Parameter for the pressure solution. */
  su2double rho_0;   /*!< \brief Parameter for the density solution. */
  su2double rho_x;   /*!< \brief Parameter for the density solution. */
  su2double rho_xy;  /*!< \brief Parameter for the density solution. */
  su2double rho_y;   /*!< \brief Parameter for the density solution. */
  su2double u_0;     /*!< \brief Parameter for the x-velocity solution. */
  su2double u_x;     /*!< \brief Parameter for the x-velocity solution. */
  su2double u_xy;    /*!< \brief Parameter for the x-velocity solution. */
  su2double u_y;     /*!< \brief Parameter for the x-velocity solution. */
  su2double v_0;     /*!< \brief Parameter for the y-velocity solution. */
  su2double v_x;     /*!< \brief Parameter for the y-velocity solution. */
  su2double v_xy;    /*!< \brief Parameter for the y-velocity solution. */
  su2double v_y;     /*!< \brief Parameter for the y-velocity solution. */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CMMSNSUnitQuadSolution(void);

  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CMMSNSUnitQuadSolution(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_iMesh, CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CMMSNSUnitQuadSolution(void) override;

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const su2double* val_coords, const su2double val_t, su2double* val_solution) const override;

  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const su2double* val_coords, const su2double val_t, su2double* val_solution) const override;

  /*!
   * \brief Get the source term for the manufactured solution (MMS).
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetMMSSourceTerm(const su2double* val_coords, const su2double val_t, su2double* val_source) const override;

  /*!
   * \brief Whether or not this verification solution is a manufactured solution.
   * \return  - True, because this is a manufactured solution.
   */
  bool IsManufacturedSolution(void) const override;
};
