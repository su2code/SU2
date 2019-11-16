/*!
 * \file gauss_structure.hpp
 * \brief Headers of the Finite Element structure (gaussian points)
 *        The subroutines and functions are in the <i>gauss_structure.cpp</i> file.
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "toolboxes/C2DContainer.hpp"

/*!
 * \class CGaussVariable
 * \brief Main class for defining the gaussian points.
 * \author R. Sanchez
 */
class CGaussVariable {
protected:

  su2activematrix GradNi_Xj;      /*!< \brief Gradient of the shape functions N[i] wrt the reference configuration. */
  su2activematrix GradNi_xj;      /*!< \brief Gradient of the shape functions N[i] wrt the current configuration. */
  su2activevector Ni;             /*!< \brief Shape functions N[i] at the gaussian point. */
  su2double J_X = 0.0;            /*!< \brief Element Jacobian evaluated at this Gauss Point wrt the reference configuration. */
  su2double J_x = 0.0;            /*!< \brief Element Jacobian evaluated at this Gauss Point wrt the current configuration. */
  unsigned short iGaussPoint = 0; /*!< \brief Identifier of the Gauss point considered. */

public:
  /*!
   * \brief Class constructor
   * \param[in] val_iGauss - ID of the Gaussian Point
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGaussVariable(unsigned short val_iGauss, unsigned short val_nDim, unsigned short val_nNodes)
    : J_X(0.0), J_x(0.0), iGaussPoint(val_iGauss)
  {
    GradNi_Xj.resize(val_nNodes,val_nDim) = su2double(0.0);
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

/*!
 * \class CElementProperty
 * \brief Main class for defining the element properties.
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 */
class CProperty {
protected:

  unsigned long iMat_Prop = 0;  /*!< \brief Index of the properties (E, Nu) for the structural model used. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] valMat_Prop - Index of the physical properties (E,nu,rho,rho_dead_load) assigned to the element.
   */
  CProperty(unsigned long valMat_Prop) : iMat_Prop(valMat_Prop) {}

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CProperty(void) = default;

  /*!
   * \brief Get the material model to use for the element.
   */
  inline virtual unsigned long GetMat_Mod(void) const { return 0; } 

  /*!
   * \brief Get index of the physical properties.
   */
  inline unsigned long GetMat_Prop(void) const { return iMat_Prop; }

  /*!
   * \brief Get index of the electric properties.
   */
  inline virtual unsigned long GetElectric_Prop(void) const { return 0; }

  /*!
   * \brief Get index of the design variable.
   */
  inline virtual unsigned long GetDV(void) const { return 0; }

  /*!
   * \brief Set the Design density (topology optimization variable).
   */
  inline virtual void SetDesignDensity(su2double valDensity) {}

  /*!
   * \brief Get the value of the Design density.
   */
  inline virtual su2double GetDesignDensity(void) const { return 0.0; }

  /*!
   * \brief Set the Physical density (used to penalize element stiffness by the FEM solver).
   */
  inline virtual void SetPhysicalDensity(su2double valDensity) {}

  /*!
   * \brief Get the value of the Physical density.
   */
  inline virtual su2double GetPhysicalDensity(void) const { return 0.0; }

  /*!
   * \brief Extract the derivative of the Design density.
   */
  inline virtual su2double GetAdjointDensity(void) const { return 0.0; }

  /*!
   * \brief Register the Design density as an AD input variable.
   */
  inline virtual void RegisterDensity(void) {}
};

/*!
 * \class CElementProperty
 * \brief Main class for defining the element properties.
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 */
class CElementProperty final : public CProperty {
protected:

  unsigned long iMat_Mod = 0;        /*!< \brief Index of the material model used. */
  unsigned long iElectric_Prop = 0;  /*!< \brief Index of the electric properties (Em) for the structural model used. */
  unsigned long iDV = 0;             /*!< \brief Index of the group of design variables to which the element belongs. */
  su2double design_rho = 1.0;        /*!< \brief Value of the design density for material-based topology optimization. */
  su2double physical_rho = 1.0;      /*!< \brief Value of the physical density for material-based topology optimization. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] valMat_Model - Type of material model (i.e. numerics) for the element, see FEA_TERM etc. in option_structure.hpp.
   * \param[in] valMat_Prop - Index of the physical properties (E,nu,rho,rho_dead_load) assigned to the element.
   * \param[in] valElectric_Prop - Index of the electric properties.
   * \param[in] valDV - Index of the design variable assigned to the element (bound to a material property by "DESIGN_VARIABLE_FEA").
   * \param[in] valDensity - Value for Design and Physical densities (topology optimization variables).
   */
  CElementProperty(unsigned long valMat_Model, unsigned long valMat_Prop,
                   unsigned long valElectric_Prop, unsigned long valDV,
                   su2double valDensity = 1.0) : CProperty(valMat_Prop),
    iMat_Mod(valMat_Model), iElectric_Prop(valElectric_Prop),
    iDV(valDV), design_rho(valDensity), physical_rho(valDensity) {}

  /*!
   * \brief Destructor of the class.
   */
  ~CElementProperty(void) = default;

  /*!
   * \brief Get the material model to use for the element.
   */
  inline unsigned long GetMat_Mod(void) const override { return iMat_Mod; }

  /*!
   * \brief Get index of the electric properties.
   */
  inline unsigned long GetElectric_Prop(void) const override { return iElectric_Prop; }

  /*!
   * \brief Get index of the design variable.
   */
  inline unsigned long GetDV(void) const override { return iDV; }
  
  /*!
   * \brief Set the Design density (topology optimization variable).
   */
  inline void SetDesignDensity(su2double valDensity) override { design_rho = valDensity; }
  
  /*!
   * \brief Get the value of the Design density.
   */
  inline su2double GetDesignDensity(void) const override { return design_rho; }
  
  /*!
   * \brief Set the Physical density (used to penalize element stiffness by the FEM solver).
   */
  inline void SetPhysicalDensity(su2double valDensity) override { physical_rho = valDensity; }
  
  /*!
   * \brief Get the value of the Physical density.
   */
  inline su2double GetPhysicalDensity(void) const override { return physical_rho; }
  
  /*!
   * \brief Extract the derivative of the Design density.
   */
  inline su2double GetAdjointDensity(void) const override { return SU2_TYPE::GetDerivative(design_rho); }
  
  /*!
   * \brief Register the Design density as an AD input variable.
   */
  inline void RegisterDensity(void) override { AD::RegisterInput(design_rho); }
};
