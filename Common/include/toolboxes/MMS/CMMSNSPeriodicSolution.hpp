/*!
 * \file CMMSNSPeriodicSolution.hpp
 * \brief Header file for the class CMMSNSPeriodicSolution.
 *        The implementations are in the <i>CMMSNSPeriodicSolution</i> file.
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
 * \class CMMSNSPeriodicSolution
 * \brief Class to define the required data for the manufactured solution of the
 laminar Navier-Stokes equations on a periodic quad [-1,1]x[-1,1].
 * \author T. Economon, E. van der Weide
 */
class CMMSNSPeriodicSolution: public CVerificationSolution {
  
protected:
  
  /*--- Variables that define the solution and MMS source term. ---*/
  
  su2double Gamma;        /*!< \brief Specific heat ratio. */
  su2double RGas;         /*!< \brief Gas constant. */
  su2double Viscosity;    /*!< \brief Viscosity, must be constant. */
  su2double Conductivity; /*!< \brief Thermal conductivity, must be constant. */
  
  /* Constants, which describe this manufactured solution. This is a viscous
   solution on a [-1,1]x[-1,1] quad, where the conservative variables vary as a
   combination of sine functions. The solution is periodic and can be
   used to test the periodic BCs in the solver. */
  
  su2double k;  /*!< \brief Factor in the sine function. */
  su2double a;  /*!< \brief Parameter added as an offset. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CMMSNSPeriodicSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CMMSNSPeriodicSolution(unsigned short val_nDim,
                         unsigned short val_nvar,
                         unsigned short val_iMesh,
                         CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CMMSNSPeriodicSolution(void);
  
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
   * \brief Get the exact primitive variable gradients at the current position and time.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_gradient - Array where the exact solution is stored.
   */
  void GetPrimitiveGradient(const su2double *val_coords,
                            const su2double val_t,
                            su2double       **val_gradient);
  
  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetBCState(const su2double *val_coords,
                  const su2double val_t,
                  su2double       *val_solution);
  
  /*!
   * \brief Get the source term for the manufactured solution (MMS).
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetMMSSourceTerm(const su2double *val_coords,
                        const su2double val_t,
                        su2double       *val_source);
  
  /*!
   * \brief Whether or not this verification solution is a manufactured solution.
   * \return  - True, because this is a manufactured solution.
   */
  bool IsManufacturedSolution(void);
  
  /*!
   * \brief Whether or not we have the exact primitive gradient for this solution.
   * \return  - True since we have the exact gradient from the manufactured solution.
   */
  bool ExactPrimitiveGradientKnown(void);
  
};
