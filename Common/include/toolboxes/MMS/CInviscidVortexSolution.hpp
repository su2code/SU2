/*!
 * \file CInviscidVortexSolution.hpp
 * \brief Header file for the class CInviscidVortexSolution.
 *        The implementations are in the <i>CInviscidVortexSolution.cpp</i> file.
 * \author T. Economon, E. van der Weide
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

#include <cmath>
#include "CVerificationSolution.hpp"

/*!
 * \class CInviscidVortexSolution
 * \brief Class to define the required data for the Inviscid Vortex.
 * \author E. van der Weide, T. Economon
 */
class CInviscidVortexSolution: public CVerificationSolution {

protected:

  /*--- Specific conditions for the inviscid vortex. ---*/
  su2double MachVortex;     /*!< \brief Mach number of the undisturbed flow. */
  su2double x0Vortex;       /*!< \brief Initial x-coordinate of the vortex center. */
  su2double y0Vortex;       /*!< \brief Initial y-coordinate of the vortex center. */
  su2double RVortex;        /*!< \brief Radius of the vortex. */
  su2double epsVortex;      /*!< \brief Strength of the vortex. */
  su2double thetaVortex;    /*!< \brief Advection angle (in degrees) of the vortex. */

  /*--- Variables involving gamma. */
  su2double Gamma;        /*!< \brief Gamma */
  su2double Gm1;          /*!< \brief Gamma minus 1 */
  su2double ovGm1;        /*!< \brief 1 over Gamma minus 1 */
  su2double gamOvGm1;     /*!< \brief Gamma over Gamma minus 1 */

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CInviscidVortexSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CInviscidVortexSolution(unsigned short val_nDim,
                          unsigned short val_nvar,
                          unsigned short val_iMesh,
                          CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CInviscidVortexSolution(void);

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetSolution(const su2double *val_coords,
                   const su2double val_t,
                   su2double       *val_solution);

  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const su2double *val_coords,
                  const su2double val_t,
                  su2double       *val_solution);
};
