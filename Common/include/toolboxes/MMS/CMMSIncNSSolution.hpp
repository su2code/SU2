/*!
 * \file CMMSIncNSSolution.hpp
 * \brief Header file for the class CMMSIncNSSolution.
 *        The implementations are in the <i>CMMSIncNSSolution.cpp</i> file.
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
 * \class CMMSIncNSSolution
 * \brief Class to define the required data for the manufactured solution of the
 *        laminar incompressible Navier-Stokes equations.
 * \author T. Economon, E. van der Weide
 */
class CMMSIncNSSolution: public CVerificationSolution {
  
protected:
  
  /*--- Variables that define the solution and MMS source term. ---*/
  su2double Viscosity;    /*!< \brief Viscosity, must be constant. */
  su2double Density;      /*!< \brief Density, must be constant. */
  su2double Temperature;  /*!< \brief Temperature, just to be safe. */
  
  /*--- Constants, which describe this manufactured solution. This is a
   viscous solution where the primitive variables vary as a combination
   of sine and cosine functions. The solution is from Salari K, and
   Knupp P, "Code verification by the method of manufactured solutions,"
   SAND 2000-1444, Sandia National Laboratories, Albuquerque, NM, 2000. ---*/
  
  su2double P_0;      /*!< \brief Parameter for the pressure solution. */
  su2double u_0;      /*!< \brief Parameter for the x-velocity solution. */
  su2double v_0;      /*!< \brief Parameter for the y-velocity solution. */
  su2double epsilon;  /*!< \brief Parameter for the velocity solutions. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CMMSIncNSSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Configuration of the particular problem.
   */
  CMMSIncNSSolution(unsigned short val_nDim,
                    unsigned short val_nvar,
                    unsigned short val_iMesh,
                    CConfig*       config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CMMSIncNSSolution(void);
  
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
};
