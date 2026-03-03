/*!
 * \file CReorderElements.hpp
 * \brief Header file for the class CReorderElements.
 *        The implementations are in the <i>CReorderElements.cpp</i> file.
 * \author E. van der Weide
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include <iostream>

using namespace std;

/*!
 * \class CReorderElements
 * \brief Class, used to reorder the owned elements after the partitioning.
 * \author E. van der Weide
 * \version 8.2.0 "Harrier"
 */
class CReorderElements {
 private:
  unsigned long globalElemID; /*!< \brief Global element ID of the element. */
  unsigned short timeLevel;   /*!< \brief Time level of the element. Only relevant
                                          for time accurate local time stepping. */
  bool commSolution;          /*!< \brief Whether or not the solution must be
                                          communicated to other ranks. */
  unsigned short elemType;    /*!< \brief Short hand for the element type, Which
                                          stored info of the VTK_Type, polynomial
                                          degree of the solution and whether or
                                          not the Jacobian is constant. */
 public:
  /*!
   * \brief Constructor of the class, set the member variables to the arguments.
   */
  CReorderElements(const unsigned long val_GlobalElemID, const unsigned short val_TimeLevel,
                   const bool val_CommSolution, const unsigned short val_VTK_Type, const unsigned short val_nPolySol,
                   const bool val_JacConstant);

  /*!
   * \brief Default constructor of the class. Disabled.
   */
  CReorderElements() = delete;

  /*!
   * \brief Less than operator of the class. Needed for the sorting.
   */
  bool operator<(const CReorderElements& other) const;

  /*!
   * \brief Function to make available the variable commSolution.
   * \return Whether or not the solution of the element must be communicated.
   */
  inline bool GetCommSolution() const { return commSolution; }

  /*!
   * \brief Function to make available the element type of the element.
   * \return The value of elemType, which stores the VTK type, polynomial degree
             and whether or not the Jacobian is constant.
   */
  inline unsigned short GetElemType() const { return elemType; }

  /*!
   * \brief Function to make available the global element ID.
   * \return The global element ID of the element.
   */
  inline unsigned long GetGlobalElemID() const { return globalElemID; }

  /*!
   * \brief Function to make available the time level.
   * \return The time level of the element.
   */
  inline unsigned short GetTimeLevel() const { return timeLevel; }

  /*!
   * \brief Function, which sets the value of commSolution.
   * \param[in] val_CommSolution  - value to which commSolution must be set.
   */
  inline void SetCommSolution(const bool val_CommSolution) { commSolution = val_CommSolution; }
};
