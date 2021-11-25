/*!
 * \file CGradSmoothing.hpp
 * \brief Declarations and inlines of the numerics class for gradient smoothing.
 * \author T.Dick
 * \version 7.2.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "CNumerics.hpp"
#include "../../../Common/include/geometry/elements/CElement.hpp"

/*!
 * \class CGradSmoothing
 * \brief Class for computing the stiffness matrix of the sobolev problem
 * \ingroup Grad_Smooth
 * \author T. Dick
 */
class CGradSmoothing : public CNumerics {

    su2double **val_DHiDHj;
    su2double *Ni_Vec;

public:

  /*!
   * \brief Constructor of the class.
   */
  CGradSmoothing(void);

  /*!
   * \brief Constructor of the class (overload).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGradSmoothing(unsigned short val_nDim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CGradSmoothing(void);

  void Compute_Tangent_Matrix(CElement *element_container, const CConfig *config);

};
