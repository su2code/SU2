/*!
 * \file CSobolevSmoothingVariable.cpp
 * \brief Definition of the solution for gradient smoothing.
 * \author T.Dick
 * \version 7.0.0 "Blackbird"
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

#include "CVariable.hpp"

class CSobolevSmoothingVariable : public CVariable {
public:

    bool* boundary_vertex;  /*!< \brief Stores if a point belongs to the boundary of a boundary. */
    unsigned long nBoundPoints;

  /*!
   * \brief Constructor of the class.
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSobolevSmoothingVariable(unsigned long npoint, unsigned long ndim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSobolevSmoothingVariable();

  /*!
   * \brief Mark a point as boundary of a boundary
   */
  void MarkAsBoundaryPoint(unsigned long iPoint);

  /*!
   * \brief return wether a point is a boundary of a boundary
   */
  bool IsBoundaryPoint(unsigned long iPoint);

  /*!
   * \brief return the number of marked points
   */
  unsigned int GetNBoundPoints();

};
