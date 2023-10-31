/*!
 * \file CFEANonlinearElasticity.hpp
 * \brief Declaration and inlines of the nonlinear elasticity FE numerics class.
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
 * \class CFEANonlinearElasticity
 * \brief Abstract class for computing the stiffness matrix of a nonlinear elasticity problem.
 *        This class does not implement a particular model, that will be done by its children.
 * \note  In addition to Compute_Constitutive_Matrix, derived classes MUST further implement
 *        Compute_Plane_Stress_Term and Compute_Stress_Tensor.
 * \ingroup Elasticity_Equations
 * \author R.Sanchez
 * \version 8.0.0 "Harrier"
 */
class CFEANonlinearElasticity : public CFEAElasticity {

protected:

  su2double **F_Mat;             /*!< \brief Deformation gradient. */
  su2double **b_Mat;             /*!< \brief Left Cauchy-Green Tensor. */
  su2double **currentCoord;      /*!< \brief Current coordinates. */
  su2double **Stress_Tensor;     /*!< \brief Cauchy stress tensor */

  su2double **FmT_Mat;           /*!< \brief Deformation gradient inverse and transpose. */

  su2double **KAux_P_ab;         /*!< \brief Auxiliar matrix for the pressure term */
  su2double *KAux_t_a;           /*!< \brief Auxiliar matrix for the pressure term */

  su2double J_F;                 /*!< \brief Jacobian of the transformation (determinant of F) */

  su2double f33;                 /*!< \brief Plane stress term for non-linear 2D plane stress analysis */

  bool nearly_incompressible;    /*!< \brief Boolean to consider nearly_incompressible effects */

  su2double **F_Mat_Iso;         /*!< \brief Isocoric component of the deformation gradient. */
  su2double **b_Mat_Iso;         /*!< \brief Isocoric component of the left Cauchy-Green tensor. */

  su2double C10, D1;             /*!< \brief C10 = Mu/2. D1 = Kappa/2. */
  su2double J_F_Iso;             /*!< \brief J_F_Iso: det(F)^-1/3. */

  su2double cijkl[3][3][3][3];   /*!< \brief Constitutive tensor i,j,k,l (defined only for incompressibility - near inc.). */

  bool maxwell_stress;           /*!< \brief Consider the effects of the dielectric loads */

  su2double *EField_Ref_Unit,    /*!< \brief Electric Field, unitary, in the reference configuration. */
  *EField_Ref_Mod;               /*!< \brief Electric Field, modulus, in the reference configuration. */
  su2double *EField_Curr_Unit;   /*!< \brief Auxiliary vector for the unitary Electric Field in the current configuration. */
  unsigned short nElectric_Field,
  nDim_Electric_Field;

  su2double *ke_DE_i;            /*!< \brief Electric Constant for Dielectric Elastomers. */

  su2double ke_DE;               /*!< \brief Electric Constant for Dielectric Elastomers. */
  su2double EFieldMod_Ref;       /*!< \brief Modulus of the electric field in the reference configuration. */

public:
  /*!
   * \brief Default constructor deleted to avoid instantiation of derived classes without sizes.
   */
  CFEANonlinearElasticity() = delete;

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEANonlinearElasticity(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEANonlinearElasticity(void) override;

  /*!
   * \brief Set element electric field.
   * \param[in] i_DV - Index of the variable.
   * \param[in] val_EField - Value of the field.
   */
  inline void Set_ElectricField(unsigned short i_DV, su2double val_EField) final {
    EField_Ref_Mod[i_DV] = val_EField;
  }

  /*!
   * \brief Build the tangent stiffness matrix of an element.
   * \param[in,out] element_container - Element whose tangent matrix is being built.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Tangent_Matrix(CElement *element_container, const CConfig *config) final;

  /*!
   * \brief Compute the nodal stress terms for an element.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_NodalStress_Term(CElement *element_container, const CConfig *config) final;

  /*!
   * \brief Compute averaged nodal stresses (for post processing).
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  su2double Compute_Averaged_NodalStress(CElement *element_container, const CConfig *config) final;

protected:
  /*!
   * \brief Compute the plane stress term.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  virtual void Compute_Plane_Stress_Term(CElement *element_container, const CConfig *config) = 0;

  /*!
   * \brief Compute the stress tensor.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  virtual void Compute_Stress_Tensor(CElement *element_container, const CConfig *config) = 0;

  /*!
   * \brief Update an element with Maxwell's stress.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Add_MaxwellStress(CElement *element_container, const CConfig *config);

  /*!
   * \brief Set element electric properties.
   * \param[in] element_container - Element defining the properties.
   * \param[in] config - Definition of the problem.
   */
  void SetElectric_Properties(const CElement *element_container, const CConfig *config);

  /*!
   * \brief TODO: Describe what this does.
   */
  void Compute_FmT_Mat(void);

  /*!
   * \brief TODO: Describe what this does.
   */
  void Compute_Isochoric_F_b(void);

  /*!
   * \brief Assign elements of constitutive tensor to matrix D.
   */
  void Assign_cijkl_D_Mat(void);

};
