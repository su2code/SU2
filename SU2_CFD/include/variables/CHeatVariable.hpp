/*!
 * \file CHeatVariable.hpp
 * \brief Class for defining the variables of the finite-volume heat equation solver.
 * \author F. Palacios, T. Economon
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

/*!
 * \class CHeatVariable
 * \brief Class for defining the variables of the finite-volume heat equation solver.
 * \author O. Burghardt
 * \version 7.3.1 "Blackbird"
 */
class CHeatVariable final : public CVariable {
protected:
  CVectorOfMatrix& Gradient_Reconstruction;  /*!< \brief Reference to the gradient of the primitive variables for MUSCL reconstruction for the convective term */
  CVectorOfMatrix Gradient_Aux;              /*!< \brief Auxiliary structure to store a second gradient for reconstruction, if required. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] heat - Values of the Heat solution (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CHeatVariable(su2double heat, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CHeatVariable() override = default;

  /*!
   * \brief Get the array of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \return Array of the reconstruction variables gradient at a node.
   */
  inline CMatrixView<su2double> GetGradient_Reconstruction(unsigned long iPoint) final { return Gradient_Reconstruction[iPoint]; }

  /*!
   * \brief Get the reconstruction gradient for primitive variable at all points.
   * \return Reference to variable reconstruction gradient.
   */
  inline CVectorOfMatrix& GetGradient_Reconstruction() final { return Gradient_Reconstruction; }
  inline const CVectorOfMatrix& GetGradient_Reconstruction() const final { return Gradient_Reconstruction; }

  /*!
   * \brief Get the temperature of the point.
   * \return Value of the temperature of the point.
   */
  inline su2double GetTemperature(unsigned long iPoint) const final { return Solution(iPoint,0); }

};
