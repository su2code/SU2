/*!
 * \file CNSUnitQuadSolution.hpp
 * \brief Header file for the class CNSUnitQuadSolution.hpp.
 *        The implementations are in the <i>CNSUnitQuadSolution.cpp</i> file.
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
 * \class CNSUnitQuadSolution
 * \brief Class to define the required data for the Navier-Stokes solution
          on a unit quad, heat conduction is neglected.
 * \author E. van der Weide, T. Economon
 */
class CNSUnitQuadSolution final : public CVerificationSolution {
 protected:
  /*--- Variables that define the soltion. ---*/
  su2double Gm1;       /*!< \brief Gamma minus 1 */
  su2double flowAngle; /*!< \brief Angle of the velocity vector in radians. */
  su2double Viscosity; /*!< \brief Viscosity, must be constant. */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CNSUnitQuadSolution(void);

  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CNSUnitQuadSolution(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_iMesh, CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CNSUnitQuadSolution(void) override;

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
};
