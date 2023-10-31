/*!
 * \file task_definition.hpp
 * \brief Header of the task definition class for the SU2 solvers.
 * \author E. van der Weide, T. Economon
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

#include "../../Common/include/parallelization/mpi_structure.hpp"

#include <iostream>

using namespace std;

/*!
 * \class CTaskDefinition
 * \brief Class for defining a task to be carried out
 * \author: E. van der Weide, T. Economon
 * \version 8.0.0 "Harrier"
 */
class CTaskDefinition {

public:

  /*!
   * \brief Enumerated type, which defines the tasks to be carried out.
   */
  enum SOLVER_TASK {
    NO_TASK                                           =  0,   /*!< \brief Default value used for checking. */
    ADER_PREDICTOR_STEP_COMM_ELEMENTS                 =  1,   /*!< \brief ADER predictor step for elements whose solution must be communicated. */
    ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS             =  2,   /*!< \brief ADER predictor step for internal elements. */
    INITIATE_MPI_COMMUNICATION                        =  3,   /*!< \brief Start the communication of the conserved variables. */
    COMPLETE_MPI_COMMUNICATION                        =  4,   /*!< \brief Complete the communication of the conserved variables. */
    INITIATE_REVERSE_MPI_COMMUNICATION                =  5,   /*!< \brief Start the communication of the residuals. */
    COMPLETE_REVERSE_MPI_COMMUNICATION                =  6,   /*!< \brief Complete the communication of the residuals. */
    ADER_TIME_INTERPOLATE_OWNED_ELEMENTS              =  7,   /*!< \brief Carry out a time interpolation to an integration point for the owned elements. */
    ADER_TIME_INTERPOLATE_HALO_ELEMENTS               =  8,   /*!< \brief Carry out a time interpolation to an integration point for the halo elements. */
    SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS          =  9,   /*!< \brief Compute the shock capturing artificial viscosity for the owned elements, if needed. */
    SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS           = 10,   /*!< \brief Compute the shock capturing artificial viscosity for the halo elements, if needed. */
    VOLUME_RESIDUAL                                   = 11,   /*!< \brief Compute the contribution to the residual from the volume integral. */
    SURFACE_RESIDUAL_OWNED_ELEMENTS                   = 12,   /*!< \brief Compute the contribution to the residual from the interior surface integral between owned elements. */
    SURFACE_RESIDUAL_HALO_ELEMENTS                    = 13,   /*!< \brief Compute the contribution to the residual from the interior surface integral between an owned and halo element. */
    BOUNDARY_CONDITIONS_DEPEND_ON_OWNED               = 14,   /*!< \brief Compute the contribution to the residual from the boundary conditions that only depend on owned elements. */
    BOUNDARY_CONDITIONS_DEPEND_ON_HALO                = 15,   /*!< \brief Compute the contribution to the residual from the boundary conditions that depend on halo elements. */
    SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS      = 16,   /*!< \brief Sum up all contributions to the residual for the owned elements. */
    SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS       = 17,   /*!< \brief Sum up all contributions to the residual for the halo elements. */
    ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS = 18,   /*!< \brief Accumlate the ADER space time residual for the owned elements. */
    ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS  = 19,   /*!< \brief Accumlate the ADER space time residual for the halo elements. */
    MULTIPLY_INVERSE_MASS_MATRIX                      = 20,   /*!< \brief Multiply the accumulated residual with the inverse of the mass matrix. */
    ADER_UPDATE_SOLUTION                              = 21    /*!< \brief Update the solution for the ADER scheme. */
  };

  SOLVER_TASK    task;                  /*!< \brief Task to be carried out. */
  unsigned short timeLevel;             /*!< \brief Time level of the task to be carried out. */
  unsigned short intPointADER;          /*!< \brief Time integration point for ADER, if relevant for the task. */
  bool           secondPartTimeIntADER; /*!< \brief Whether or not this is the second part of the time interval for elements
                                                    adjacent to a lower time level. */
  unsigned short nIndMustBeCompleted;   /*!< \brief Number of relevant indices in indMustBeCompleted. */
  int            indMustBeCompleted[5]; /*!< \brief Indices in the list of tasks that must be completed before this task can be carried out. */

  /*!
   * \brief Constructor of the class.
   */
  CTaskDefinition(void);

  /*!
   * \brief Alternative constructor of the class, which initializes the member
            variables to the given arguments.
   * \param[in] val_task                - Task to be set.
   * \param[in] val_timeLevel           - Time level to be set.
   * \param[in] val_ind0MustBeCompleted - Completed index to be set, defaulted to -1.
   * \param[in] val_ind1MustBeCompleted - Completed index to be set, defaulted to -1.
   * \param[in] val_ind2MustBeCompleted - Completed index to be set, defaulted to -1.
   * \param[in] val_ind3MustBeCompleted - Completed index to be set, defaulted to -1.
   * \param[in] val_ind4MustBeCompleted - Completed index to be set, defaulted to -1.
   */
  CTaskDefinition(SOLVER_TASK    val_task,
                  unsigned short val_timeLevel,
                  int            val_ind0MustBeCompleted = -1,
                  int            val_ind1MustBeCompleted = -1,
                  int            val_ind2MustBeCompleted = -1,
                  int            val_ind3MustBeCompleted = -1,
                  int            val_ind4MustBeCompleted = -1);

  /*!
   * \brief Destructor of the class.
   */
  ~CTaskDefinition(void);

  /*!
   * \brief Copy constructor of the class.
   */
  CTaskDefinition(const CTaskDefinition &other);

  /*!
   * \brief Assignment operator.
   */
  CTaskDefinition& operator=(const CTaskDefinition &other);

private:
  /*!
   * \brief Function that copies the data from other into the member variables.
   * \param[in] other - Object from which the data must be copied.
   */
  void Copy(const CTaskDefinition &other);
};

#include "task_definition.inl"
