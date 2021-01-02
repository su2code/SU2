/*!
 * \file CGemmBase.hpp
 * \brief Base class for carrying out a GEMM multiplication.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
 * \class CGemmBase
 * \brief Base class to carry out a GEMM multiplication.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 */
class CGemmBase {

public:
  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Default constructor of the class.
   */
  CGemmBase() = default;

  /*!
   * \brief Destructor.
   */
  virtual ~CGemmBase() = default;

  /*-----------------------------------------------------------------------------------*/
  /*---                        Public enumerated types.                             ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Types of GEMMs that are present. Needed for the tensor products.
   */
  enum ENUM_GEMM_TYPE {
    DOFS_TO_INT = 0,    /*!< \brief Create data in the integration points from the DOFs. */
    INT_TO_DOFS = 1     /*!< \brief Create data in the DOFs from the integration points. */
  };
};
