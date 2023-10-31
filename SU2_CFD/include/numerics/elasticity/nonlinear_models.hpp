/*!
 * \file nonlinear_models.hpp
 * \brief Declarations of nonlinear constitutive models.
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

#include "CFEANonlinearElasticity.hpp"


/*!
 * \class CFEM_NeoHookean_Comp
 * \brief Class for computing the constitutive and stress tensors for a neo-Hookean material model, compressible.
 * \ingroup Elasticity_Equations
 * \author R.Sanchez
 * \version 8.0.0 "Harrier"
 */
class CFEM_NeoHookean_Comp final : public CFEANonlinearElasticity {

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_NeoHookean_Comp(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_NeoHookean_Comp(void) override = default;

private:
  /*!
   * \brief Compute the plane stress term.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Plane_Stress_Term(CElement *element_container, const CConfig *config) override;

  /*!
   * \brief Compute the constitutive matrix.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Constitutive_Matrix(CElement *element_container, const CConfig *config) override;

  /*!
   * \brief Compute the stress tensor.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Stress_Tensor(CElement *element_container, const CConfig *config) override;

};


/*!
 * \class CFEM_NeoHookean_Comp
 * \brief Constitutive and stress tensors for a Knowles stored-energy function, nearly incompressible.
 * \ingroup Elasticity_Equations
 * \author R.Sanchez
 * \version 8.0.0 "Harrier"
 */
class CFEM_Knowles_NearInc final : public CFEANonlinearElasticity {

  su2double trbbar, term1, term2, Ek, Pr; /*!< \brief Variables of the model calculation. */
  su2double Bk, Nk;                       /*!< \brief Parameters b and n of the model. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_Knowles_NearInc(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_Knowles_NearInc(void) override = default;

private:
  /*!
   * \brief Compute the plane stress term.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Plane_Stress_Term(CElement *element_container, const CConfig *config) override;

  /*!
   * \brief Compute the constitutive matrix.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Constitutive_Matrix(CElement *element_container, const CConfig *config) override;

  /*!
   * \brief Compute the stress tensor.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Stress_Tensor(CElement *element_container, const CConfig *config) override;

};


/*!
 * \class CFEM_DielectricElastomer
 * \brief Class for computing the constitutive and stress tensors for a dielectric elastomer.
 * \ingroup Elasticity_Equations
 * \author R.Sanchez
 * \version 8.0.0 "Harrier"
 */
class CFEM_DielectricElastomer final : public CFEANonlinearElasticity {

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_DielectricElastomer(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_DielectricElastomer(void) override = default;

private:
  /*!
   * \brief Compute the plane stress term.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  inline void Compute_Plane_Stress_Term(CElement *element_container, const CConfig *config) override { };

  /*!
   * \brief Compute the constitutive matrix.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Constitutive_Matrix(CElement *element_container, const CConfig *config) override;

  /*!
   * \brief Compute the stress tensor.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Stress_Tensor(CElement *element_container, const CConfig *config) override;

};


/*!
 * \class CFEM_IdealDE
 * \brief Class for computing the constitutive and stress tensors for a nearly-incompressible ideal DE.
 * \ingroup Elasticity_Equations
 * \author R.Sanchez
 * \version 8.0.0 "Harrier"
 */
class CFEM_IdealDE final : public CFEANonlinearElasticity {

  su2double trbbar, Eg, Eg23, Ek, Pr; /*!< \brief Variables of the model calculation. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_IdealDE(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_IdealDE(void) override = default;

private:
  /*!
   * \brief Compute the plane stress term.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Plane_Stress_Term(CElement *element_container, const CConfig *config) override;

  /*!
   * \brief Compute the constitutive matrix.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Constitutive_Matrix(CElement *element_container, const CConfig *config) override;

  /*!
   * \brief Compute the stress tensor.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Stress_Tensor(CElement *element_container, const CConfig *config) override;

};
