/*!
 * \file CFEALinearElasticity.hpp
 * \brief Declaration and inlines of the linear elasticity FE numerics class.
 * \author Ruben Sanchez
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

#include "CFEAElasticity.hpp"


/*!
 * \class CFEALinearElasticity
 * \brief Class for computing the stiffness matrix of a linear, elastic problem.
 * \ingroup Elasticity_Equations
 * \author R.Sanchez
 * \version 8.0.0 "Harrier"
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
  CFEALinearElasticity(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEALinearElasticity(void) override = default;

  /*!
   * \brief Build the tangent stiffness matrix of an element.
   * \param[in,out] element_container - Element whose tangent matrix is being built.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Tangent_Matrix(CElement *element_container, const CConfig *config) final;

  /*!
   * \brief Compute averaged nodal stresses (for post processing).
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  su2double Compute_Averaged_NodalStress(CElement *element_container, const CConfig *config) final;

private:
  /*!
   * \brief Compute the constitutive matrix.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Constitutive_Matrix(CElement *element_container, const CConfig *config) final;

};


/*!
 * \class CFEAMeshElasticity
 * \brief Particular case of linear elasticity used for mesh deformation.
 * \ingroup Elasticity_Equations
 * \author R.Sanchez
 * \version 8.0.0 "Harrier"
 */
class CFEAMeshElasticity final : public CFEALinearElasticity {

  bool element_based;

public:
  /*!
   * \brief Default constructor deleted as instantiation with no argument would not allocate fields.
   */
  CFEAMeshElasticity() = delete;

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEAMeshElasticity(unsigned short val_nDim, unsigned short val_nVar, unsigned long val_nElem, const CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEAMeshElasticity(void) override = default;

  /*!
   * \brief Set the element-based local Young's modulus in mesh problems
   * \param[in] iElem - Element index.
   * \param[in] val_E - Value of elasticity modulus.
   */
  inline void SetMeshElasticProperties(unsigned long iElem, su2double val_E) override {
    if (element_based) E_i[iElem] = val_E;
  }

private:
  /*!
   * \brief Set element material properties.
   * \param[in] element_container - Element defining the properties.
   * \param[in] config - Definition of the problem.
   */
  inline void SetElement_Properties(const CElement *element, const CConfig *config) override {
    if (element_based) {
      E = E_i[element->Get_iProp()];
      Compute_Lame_Parameters();
    }
  }

};
