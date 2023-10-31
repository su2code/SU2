/*!
 * \file CGradSmoothing.hpp
 * \brief Declarations and inlines of the numerics class for gradient smoothing.
 * \author T.Dick
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

#include "../../../Common/include/geometry/elements/CElement.hpp"
#include "../../../SU2_CFD/include/numerics/CNumerics.hpp"

/*!
 * \class CGradSmoothing
 * \brief Class for computing the stiffness matrix of the Sobolev problem
 * \ingroup GradSmooth
 * \author T. Dick
 */
class CGradSmoothing final : public CNumerics {
  su2activematrix val_DHiDHj;
  su2activevector Ni_Vec;

 public:
  /*!
   * \brief Default constructor
   */
  CGradSmoothing() = delete;

  /*!
   * \brief Constructor of the class (overload).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGradSmoothing(unsigned short val_nDim, const CConfig* config);

  /*!
   * \brief Build the tangent stiffness matrix of an element.
   * \param[in,out] element_container - Element whose tangent matrix is being built.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Tangent_Matrix(CElement *element_container, const CConfig *config) override;

};
