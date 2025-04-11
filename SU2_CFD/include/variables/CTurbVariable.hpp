/*!
 * \file CTurbVariable.hpp
 * \brief Base class for defining the variables of the turbulence model.
 * \author F. Palacios, T. Economon
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

#include "CScalarVariable.hpp"

/*!
 * \class CTurbVariable
 * \brief Base class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */
class CTurbVariable : public CScalarVariable {
protected:
  VectorType muT; /*!< \brief Eddy viscosity. */

public:
  static constexpr size_t MAXNVAR = 2;
  VectorType turb_index;
  VectorType intermittency;         /*!< \brief Value of the intermittency for the trans. model. */

  VectorType Pk;
  VectorType Pw;
  VectorType Dk;
  VectorType Dw;
  VectorType PkLim;

  VectorType Res_k;
  VectorType Res_w;
  
  VectorType Jac_00;
  VectorType Jac_01;
  VectorType Jac_10;
  VectorType Jac_11;
  VectorType Jac_add;

  /*!
   * \brief Constructor of the class.
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTurbVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbVariable() override = default;

  /*!
   * \brief Get the value of the eddy viscosity.
   * \param[in] iPoint - Point index.
   * \return the value of the eddy viscosity.
   */
  inline su2double GetmuT(unsigned long iPoint) const final { return muT(iPoint); }

  /*!
   * \brief Set the value of the eddy viscosity.
   * \param[in] iPoint - Point index.
   * \param[in] val_muT - Value of the eddy viscosity.
   */
  inline void SetmuT(unsigned long iPoint, su2double val_muT) final { muT(iPoint) = val_muT; }

  /*!
    * \brief Set the value of the turbulence index.
   * \param[in] iPoint - Point index.
   * \param[in] val_turb_index - Value of the turbulence index.
   */
  inline void SetTurbIndex(unsigned long iPoint, su2double val_turb_index) final { turb_index(iPoint) = val_turb_index; }

  /*!
   * \brief Get the value of the turbulence index.
   * \param[in] iPoint - Point index.
   * \return Value of the intermittency of the turbulence index.
   */
  inline su2double GetTurbIndex(unsigned long iPoint) const final { return turb_index(iPoint); }

  /*!
   * \brief Get the intermittency of the transition model.
   * \param[in] iPoint - Point index.
   * \return Value of the intermittency of the transition model.
   */
  inline su2double GetIntermittency(unsigned long iPoint) const final { return intermittency(iPoint); }

  /*!
   * \brief Set the intermittency of the transition model.
   * \param[in] iPoint - Point index.
   * \param[in] val_intermittency - New value of the intermittency.
   */
  inline void SetIntermittency(unsigned long iPoint, su2double val_intermittency) final { intermittency(iPoint) = val_intermittency; }
  inline void SetProdDestr(unsigned long iPoint, su2double* val_ProdDestr) final { 
    Pk(iPoint) = val_ProdDestr[0];
    Dk(iPoint) = val_ProdDestr[1];
    Pw(iPoint) = val_ProdDestr[2];
    Dw(iPoint) = val_ProdDestr[3]; 
    PkLim(iPoint) = val_ProdDestr[4]; 
  }

  inline void SetResidualHere(unsigned long iPoint, su2double* val_ResidualHere) final { 
    Res_k(iPoint) = val_ResidualHere[0];
    Res_w(iPoint) = val_ResidualHere[1];
  }

  inline void SetJacobianHere(unsigned long iPoint, su2double* val_JacobianHere) final { 
    Jac_00(iPoint) = val_JacobianHere[0];
    Jac_01(iPoint) = val_JacobianHere[1];
    Jac_10(iPoint) = val_JacobianHere[2];
    Jac_11(iPoint) = val_JacobianHere[3]; 
    Jac_add(iPoint) = val_JacobianHere[4]; 
  }


  inline su2double GetProdTKE(unsigned long iPoint) const final { 
    return Pk(iPoint); 
  }
  inline su2double GetDestrTKE(unsigned long iPoint) const final { 
    return Dk(iPoint); 
  }
  inline su2double GetProdW(unsigned long iPoint) const final { 
    return Pw(iPoint); 
  }
  inline su2double GetDestrW(unsigned long iPoint) const final { 
    return Dw(iPoint); 
  }
  inline su2double GetPkLim(unsigned long iPoint) const final { 
    return PkLim(iPoint); 
  }

  inline su2double GetResTKE(unsigned long iPoint) const final { 
    return Res_k(iPoint); 
  }
  inline su2double GetResW(unsigned long iPoint) const final { 
    return Res_w(iPoint); 
  }
  
  inline su2double GetJac_00(unsigned long iPoint) const final { 
    return Jac_00(iPoint); 
  }
  inline su2double GetJac_01(unsigned long iPoint) const final { 
    return Jac_01(iPoint); 
  }
  inline su2double GetJac_10(unsigned long iPoint) const final { 
    return Jac_10(iPoint); 
  }
  inline su2double GetJac_11(unsigned long iPoint) const final { 
    return Jac_11(iPoint); 
  }
  inline su2double GetJac_add(unsigned long iPoint) const final { 
    return Jac_add(iPoint); 
  }
};

