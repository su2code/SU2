/*!
 * \file CFEAElasticity.hpp
 * \brief Declaration and inlines of the base class for elasticity problems.
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

#include "../CNumerics.hpp"
#include "../../../../Common/include/geometry/elements/CElement.hpp"

/*!
 * \class CFEAElasticity
 * \ingroup Elasticity_Equations
 * \brief Abstract class for computing the tangent matrix and the residual for structural problems.
 * \note  At the next level of abstraction (linear or not) a class must define the constitutive term.
 *        The methods we override in this class with an empty implementation are here just to better
 *        document the public interface of this class hierarchy.
 * \author R.Sanchez
 * \version 8.0.0 "Harrier"
 */
class CFEAElasticity : public CNumerics {

protected:

  enum : unsigned short {DIM_STRAIN_2D = 3,   /*!< \brief Exx, Eyy, Gxy. */
                         DIM_STRAIN_3D = 6};  /*!< \brief Exx, Eyy, Ezz, Gxy, Gxz, Gyz. */

  enum : unsigned short {NNODES_2D = 4,   /*!< \brief Maximum number of nodes for 2D problems. */
                         NNODES_3D = 8};  /*!< \brief Maximum number of nodes for 3D problems. */

  su2double E         = 1.0;              /*!< \brief Aux. variable, Young's modulus of elasticity. */
  su2double Nu        = 0.0;              /*!< \brief Aux. variable, Poisson's ratio. */
  su2double Rho_s     = 0.0;              /*!< \brief Aux. variable, Structural density. */
  su2double Rho_s_DL  = 0.0;              /*!< \brief Aux. variable, Structural density (for dead loads). */

  su2double Mu        = 0.0;              /*!< \brief Aux. variable, Lame's coeficient. */
  su2double Lambda    = 0.0;              /*!< \brief Aux. variable, Lame's coeficient. */
  su2double Kappa     = 0.0;              /*!< \brief Aux. variable, Compressibility constant. */

  su2double *E_i      = nullptr;          /*!< \brief Young's modulus of elasticity. */
  su2double *Nu_i     = nullptr;          /*!< \brief Poisson's ratio. */
  su2double *Rho_s_i  = nullptr;          /*!< \brief Structural density. */
  su2double *Rho_s_DL_i = nullptr;        /*!< \brief Structural density (for dead loads). */

  su2double **Ba_Mat = nullptr;           /*!< \brief Matrix B for node a - Auxiliary. */
  su2double **Bb_Mat = nullptr;           /*!< \brief Matrix B for node b - Auxiliary. */
  su2double *Ni_Vec  = nullptr;           /*!< \brief Vector of shape functions - Auxiliary. */
  su2double **D_Mat  = nullptr;           /*!< \brief Constitutive matrix - Auxiliary. */
  su2double **KAux_ab = nullptr;          /*!< \brief Node ab stiffness matrix - Auxiliary. */
  su2double **GradNi_Ref_Mat = nullptr;   /*!< \brief Gradients of Ni - Auxiliary. */
  su2double **GradNi_Curr_Mat = nullptr;  /*!< \brief Gradients of Ni - Auxiliary. */

  su2double *FAux_Dead_Load = nullptr;    /*!< \brief Auxiliar vector for the dead loads */

  su2double *DV_Val = nullptr;            /*!< \brief For optimization cases, value of the design variables. */
  unsigned short n_DV = 0;                /*!< \brief For optimization cases, number of design variables. */

  bool plane_stress = false;              /*!< \brief Checks if we are solving a plane stress case */

public:
  /*!
   * \brief Default constructor
   */
  CFEAElasticity() = default;

  /*!
   * \brief Constructor of the class (overload).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEAElasticity(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEAElasticity(void) override;

  /*!
   * \brief Set elasticity modulus and Poisson ratio.
   * \param[in] iVal - Index of the property.
   * \param[in] val_E - Value of the elasticity (Young) modulus.
   * \param[in] val_Nu - Value of the Poisson ratio.
   */
  inline void SetMaterial_Properties(unsigned short iVal, su2double val_E, su2double val_Nu) final {
    E_i[iVal] = val_E;
    Nu_i[iVal] = val_Nu;
  }

  /*!
   * \brief Set densities, real and for gravity loading purposes.
   * \param[in] iVal - Index of the property.
   * \param[in] val_Rho - Material density.
   * \param[in] val_Rho_DL - Density for gravity (dead) loads.
   */
  inline void SetMaterial_Density(unsigned short iVal, su2double val_Rho, su2double val_Rho_DL) final {
    Rho_s_i[iVal] = val_Rho;
    Rho_s_DL_i[iVal] = val_Rho_DL;
  }

  /*!
   * \brief Set element electric field.
   * \param[in] i_DV - Index of the variable.
   * \param[in] val_EField - Value of the field.
   */
  inline void Set_ElectricField(unsigned short i_DV, su2double val_EField) override { }

  /*!
   * \brief Set the element-based local Young's modulus in mesh problems
   * \param[in] iElem - Element index.
   * \param[in] val_E - Value of elasticity modulus.
   */
  inline void SetMeshElasticProperties(unsigned long iElem, su2double val_E) override { }

  /*!
   * \brief Set the value of a design variable.
   * \param[in] i_DV - Index of the variable.
   * \param[in] val_DV - Value of the variable.
   */
  inline void Set_DV_Val(unsigned short i_DV, su2double val_DV) final { DV_Val[i_DV] = val_DV; }

  /*!
   * \brief Get the value of a design variable.
   * \param[in] i_DV - Index of the variable.
   * \return Value of the variable.
   */
  inline su2double Get_DV_Val(unsigned short i_DV) const final { return DV_Val[i_DV]; }

  /*!
   * \brief Build the mass matrix of an element.
   * \param[in,out] element_container - Element whose mass matrix is being built.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Mass_Matrix(CElement *element_container, const CConfig *config) final;

  /*!
   * \brief Compute the nodal gravity loads for an element.
   * \param[in,out] element_container - The element for which the dead loads are computed.
   * \param[in] config - Definition of the problem.
   */
  void Compute_Dead_Load(CElement *element_container, const CConfig *config) final;

  /*!
   * \brief Build the tangent stiffness matrix of an element.
   * \param[in,out] element_container - Element whose tangent matrix is being built.
   * \param[in] config - Definition of the problem.
   */
  inline void Compute_Tangent_Matrix(CElement *element_container, const CConfig *config) override { };

  /*!
   * \brief Compute averaged nodal stresses (for post processing).
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  inline su2double Compute_Averaged_NodalStress(CElement *element_container, const CConfig *config) override { return 0; };

  /*!
   * \brief Compute VonMises stress from components Sxx Syy Sxy Szz Sxz Syz.
   */
  template<class T>
  static su2double VonMisesStress(unsigned short nDim, const T& stress) {
    if (nDim == 2) {
      su2double Sxx = stress[0], Syy = stress[1], Sxy = stress[2];

      su2double S1, S2; S1 = S2 = (Sxx+Syy)/2;
      su2double tauMax = sqrt(pow((Sxx-Syy)/2, 2) + pow(Sxy,2));
      S1 += tauMax;
      S2 -= tauMax;

      return sqrt(S1*S1+S2*S2-2*S1*S2);
    }
    else {
      su2double Sxx = stress[0], Syy = stress[1], Szz = stress[3];
      su2double Sxy = stress[2], Sxz = stress[4], Syz = stress[5];

      return sqrt(0.5*(pow(Sxx - Syy, 2) +
                       pow(Syy - Szz, 2) +
                       pow(Szz - Sxx, 2) +
                       6.0*(Sxy*Sxy+Sxz*Sxz+Syz*Syz)));
    }
  }

protected:
  /*!
   * \brief Compute the constitutive matrix, must be implemented by derived classes.
   * \param[in,out] element_container - The finite element.
   * \param[in] config - Definition of the problem.
   */
  virtual void Compute_Constitutive_Matrix(CElement *element_container, const CConfig *config) = 0;

  /*!
   * \brief Set element material properties.
   * \param[in] element_container - Element defining the properties.
   * \param[in] config - Definition of the problem.
   */
  virtual void SetElement_Properties(const CElement *element_container, const CConfig *config);

  /*!
   * \brief Read design variables from file.
   * \param[in] config - Definition of the problem.
   */
  void ReadDV(const CConfig *config);

  /*!
   * \brief Update the Lame parameters (required in AD to account for all dependencies).
   */
  inline void Compute_Lame_Parameters(void) {
    Mu     = E / (2.0*(1.0 + Nu));
    Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));
    Kappa  = Lambda + (2/3)*Mu;
  }

  /*!
   * \brief Kronecker delta.
   * \param[in] iVar - Index i.
   * \param[in] jVar - Index j.
   * \return 1 if i=j, 0 otherwise.
   */
  inline static su2double deltaij(unsigned short iVar, unsigned short jVar) {
    return su2double(iVar==jVar);
  }

};
