/*!
 * \file CGaussVariable.hpp
 * \brief Light-weight class to store Gaussian point information.
 * \author R. Sanchez
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

#include "../../containers/C2DContainer.hpp"

/*!
 * \class CGaussVariable
 * \ingroup FemAlgos
 * \brief Main class for defining the gaussian points.
 * \version 8.0.0 "Harrier"
 */
class CGaussVariable {
 protected:
  su2activematrix GradNi_Xj; /*!< \brief Gradient of the shape functions N[i] wrt the reference configuration. */
  su2activematrix GradNi_xj; /*!< \brief Gradient of the shape functions N[i] wrt the current configuration. */
  su2activevector Ni;        /*!< \brief Shape functions N[i] at the gaussian point. */
  su2double J_X = 0.0; /*!< \brief Element Jacobian evaluated at this Gauss Point wrt the reference configuration. */
  su2double J_x = 0.0; /*!< \brief Element Jacobian evaluated at this Gauss Point wrt the current configuration. */
  unsigned short iGaussPoint = 0; /*!< \brief Identifier of the Gauss point considered. */

 public:
  /*!
   * \brief Deleted default constructor as this class does not allow resizing once created.
   */
  CGaussVariable() = delete;

  /*!
   * \brief Class constructor
   * \param[in] val_iGauss - ID of the Gaussian Point
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGaussVariable(unsigned short val_iGauss, unsigned short val_nDim, unsigned short val_nNodes)
      : J_X(0.0), J_x(0.0), iGaussPoint(val_iGauss) {
    /* --- For the structural mechanics solver the dimensions (nNodes x nDim) are sufficient.
     * For the Sobolev smoothing solver dimensions (nNodes x (nDim+1)) are necessary
     * when working on a curved design surface embedded in 3D. ---*/
    GradNi_Xj.resize(val_nNodes, val_nDim + 1) = su2double(0.0);
    GradNi_xj = GradNi_Xj;

    Ni.resize(val_nNodes) = su2double(0.0);
  }

  /*!
   * \brief Destructor of the class.
   */
  ~CGaussVariable(void) = default;

  inline void SetGradNi_Xj(su2double val_GradNi_Xj, unsigned short val_iDim, unsigned short val_Ni) {
    GradNi_Xj(val_Ni, val_iDim) = val_GradNi_Xj;
  }

  inline void SetGradNi_xj(su2double val_GradNi_xj, unsigned short val_iDim, unsigned short val_Ni) {
    GradNi_xj(val_Ni, val_iDim) = val_GradNi_xj;
  }

  inline void SetNi(su2double val_ShapeNi, unsigned short val_Ni) { Ni(val_Ni) = val_ShapeNi; }

  inline void SetJ_X(su2double valJ_X) { J_X = valJ_X; }

  inline void SetJ_x(su2double valJ_x) { J_x = valJ_x; }

  inline su2double GetGradNi_Xj(unsigned short val_Ni, unsigned short val_iDim) const {
    return GradNi_Xj(val_Ni, val_iDim);
  }

  inline su2double GetGradNi_xj(unsigned short val_Ni, unsigned short val_iDim) const {
    return GradNi_xj(val_Ni, val_iDim);
  }

  inline su2double GetNi(unsigned short val_Ni) const { return Ni(val_Ni); }

  inline su2double GetJ_X(void) const { return J_X; }

  inline su2double GetJ_x(void) const { return J_x; }

  inline unsigned short Get_iGauss(void) const { return iGaussPoint; }
};
