/*!
 * \file CFEALinearElasticity.hpp
 * \brief Declaration and inlines of the linear elasticity FE numerics class.
 * \author Ruben Sanchez
 * \version 7.0.1 "Blackbird"
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

#include "CFEAElasticity.hpp"


/*!
 * \class CFEALinearElasticity
 * \brief Class for computing the stiffness matrix of a linear, elastic problem.
 * \ingroup FEM_Discr
 * \author R.Sanchez
 * \version 7.0.1 "Blackbird"
 */
class CFEALinearElasticity : public CFEAElasticity {
protected:
  su2activematrix nodalDisplacement;  /*!< \brief Nodal displacements, used to compute nodal stresses. */

  /*!
   * \brief Default constructor, protected to avoid instantiation without arguments.
   */
  CFEALinearElasticity() : CFEAElasticity() {}

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEALinearElasticity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CFEALinearElasticity(void) = default;

  /*!
   * \brief Build the tangent stiffness matrix of an element.
   * \param[in,out] element_container - Element whose tangent matrix is being built.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Tangent_Matrix(CElement *element_container, CConfig *config) final;

  /*!
   * \brief Compute averaged nodal stresses (for post processing).
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Averaged_NodalStress(CElement *element_container, CConfig *config) final;

private:
  /*!
   * \brief Compute the constitutive matrix.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Constitutive_Matrix(CElement *element_container, CConfig *config) final;

};
