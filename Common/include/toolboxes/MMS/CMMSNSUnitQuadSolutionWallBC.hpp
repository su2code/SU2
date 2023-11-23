/*!
 * \file CMMSNSUnitQuadSolutionWallBC.hpp
 * \brief Header file for the class CMMSNSUnitQuadSolutionWallBC.
 *        The implementations are in the <i>CMMSNSUnitQuadSolutionWallBC.cpp</i> file.
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
 * \class CMMSNSUnitQuadSolutionWallBC
 * \brief Class to define the required data for the manufactured solution of the
          laminar Navier-Stokes equations on a unit quad with wall boundary conditions.
 * \author E. van der Weide, T. Economon
 */
class CMMSNSUnitQuadSolutionWallBC final : public CVerificationSolution {
 protected:
  /*--- Variables that define the solution and MMS source term. ---*/
  su2double Gamma;        /*!< \brief Specific heat ratio. */
  su2double RGas;         /*!< \brief Gas constant. */
  su2double Viscosity;    /*!< \brief Viscosity, must be constant. */
  su2double Conductivity; /*!< \brief Thermal conductivity, must be constant. */
  su2double TWall;        /*!< \brief Prescribed wall temperature at the upper wall. */

  su2double Pressure_Ref; /*!< \brief Reference pressure for non-dimensionalization. */
  su2double Density_Ref;  /*!< \brief Reference density for non-dimensionalization.  */
  su2double Velocity_Ref; /*!< \brief Reference velocity for non-dimensionalization. */

  /*--- Constants, which describe this manufactured solution. The primitive variables
        rho, T, u and v are described by analytical functions in such a way that the
        lower wall is an adiabatic no-slip wall and the upper wall is an isothermal
        no-slip wall. ---*/
  su2double rho_0; /*!< \brief Constant density. */
  su2double u_0;   /*!< \brief Maximum x-velocity in the domain. */
  su2double v_0;   /*!< \brief Maximum y-velocity in the domain. */

  su2double a_T1; /*!< \brief Parameter for the temperature solution. */
  su2double a_T2; /*!< \brief Parameter for the temperature solution. */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CMMSNSUnitQuadSolutionWallBC(void);

  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CMMSNSUnitQuadSolutionWallBC(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_iMesh,
                               CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CMMSNSUnitQuadSolutionWallBC(void) override;

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
