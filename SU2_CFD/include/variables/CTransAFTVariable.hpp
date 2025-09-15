/*!
 * \file CTransAFTVariable.hpp
 * \brief Declaration of the variables of the transition model.
 * \author S. Kang.
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

#include "CTurbVariable.hpp"

/*!
 * \class CTransAFTVariable
 * \brief Transition model variables.
 * \ingroup Turbulence_Model
 * \author S. Kang.
 */

class CTransAFTVariable final : public CTurbVariable {
protected:
  VectorType TempVar1, TempVar2, TempVar3, TempVar4, TempVar5, TempVar6;
  VectorType TempVar7, TempVar8, TempVar9, TempVar10, TempVar11;
  VectorType TempVar12, TempVar13, TempVar14, TempVar15;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] AF - Amplification Factor (initialization value).
   * \param[in] LnIntermittency - Natural logarithm intermittency (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTransAFTVariable(su2double AF, su2double LnIntermittency, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTransAFTVariable() override = default;

  void SetAFT_Wonder_Func(unsigned long iPoint, su2double tempVar1, su2double tempVar2, su2double tempVar3, su2double tempVar4, su2double tempVar5, 
                                    su2double tempVar6, su2double tempVar7, su2double tempVar8, su2double tempVar9, su2double tempVar10, 
                                    su2double tempVar11, su2double tempVar12, su2double tempVar13, su2double tempVar14, su2double tempVar15) override;

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var1(unsigned long iPoint) const override { return TempVar1(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var2(unsigned long iPoint) const override { return TempVar2(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var3(unsigned long iPoint) const override { return TempVar3(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var4(unsigned long iPoint) const override { return TempVar4(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var5(unsigned long iPoint) const override { return TempVar5(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var6(unsigned long iPoint) const override { return TempVar6(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var7(unsigned long iPoint) const override { return TempVar7(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var8(unsigned long iPoint) const override { return TempVar8(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var9(unsigned long iPoint) const override { return TempVar9(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var10(unsigned long iPoint) const override { return TempVar10(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var11(unsigned long iPoint) const override { return TempVar11(iPoint); }


  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var12(unsigned long iPoint) const override { return TempVar12(iPoint); }


  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var13(unsigned long iPoint) const override { return TempVar13(iPoint); }


  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var14(unsigned long iPoint) const override { return TempVar14(iPoint); }


  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetAFT_Wonder_Func_var15(unsigned long iPoint) const override { return TempVar15(iPoint); }


};
