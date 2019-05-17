/*!
 * \file CVerificationSolution.hpp
 * \brief Header file for the base class CVerificationSolution.
 *        The implementations are in the <i>CVerificationSolution.cpp</i> file.
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
#include "../../config_structure.hpp"

/*!
 * \class CVerificationSolution
 * \brief Class for holding verification PDE solutions, e.g., phi = phi(x,y,z,t),
 *        used for initial conditions, analytic solutions, manufactured solutions.
 * \author T. Economon, E. van der Weide
 */
class CVerificationSolution {
  
protected:
  
  int rank;  /*!< \brief MPI Rank. */
  int size;  /*!< \brief MPI Size. */
  
  unsigned short nDim;  /*!< \brief Number of dimension of the problem. */
  unsigned short nVar;  /*!< \brief Number of variables of the problem  */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CVerificationSolution(void);
  
  /*!
   * \overload
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nvar  - Number of variables of the problem.
   * \param[in] val_iMesh - Multigrid level of the solver.
   * \param[in] config    - Definition of the particular problem.
   */
  CVerificationSolution(unsigned short val_nDim,
                        unsigned short val_nvar,
                        unsigned short val_iMesh,
                        CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CVerificationSolution(void);

  /*!
   * \brief Get the exact solution at the current position and time.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  virtual void GetSolution(const su2double *val_coords,
                           const su2double val_t,
                           su2double       *val_solution);
  
  /*!
   * \brief Get the exact solution at the current position and t = 0.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  void GetInitialCondition(const su2double *val_coords,
                           su2double       *val_solution);
  
  /*!
   * \brief Get the boundary conditions state for an exact solution.
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  virtual void GetBCState(const su2double *val_coords,
                          const su2double val_t,
                          su2double       *val_solution);
  
  /*!
   * \brief Get the source term for the manufactured solution (MMS).
   * \param[in] val_coords   - Cartesian coordinates of the current position.
   * \param[in] val_t        - Current physical time.
   * \param[in] val_solution - Array where the exact solution is stored.
   */
  virtual void GetMMSSourceTerm(const su2double *val_coords,
                                const su2double val_t,
                                su2double       *val_source);

  /*!
   * \brief Whether or not this verification solution is a manufactured solution.
   * \return  - False as default value. Overwrite this function for a
                manufactured solution.
   */
  virtual bool IsManufacturedSolution(void);

  /*!
   * \brief Whether or not the exact solution is known for this verification solution.
   * \return  - True as default value. Overwrite this function if the exacti
                solution is not known.
   */
  virtual bool ExactSolutionKnown(void);
  
  /*!
   * \brief Get the local error defined as the local solution minus the verification solution.
   * \param[in]  val_coords   - Cartesian coordinates of the current position.
   * \param[in]  val_solution - Array where the exact solution is stored.
   * \param[out] val_error    - Array where the local error is stored.
   */
  void GetLocalError(const su2double *val_coords,
                     const su2double val_t,
                     const su2double *GetLocalErrorval_solution,
                     su2double       *val_error);
};
