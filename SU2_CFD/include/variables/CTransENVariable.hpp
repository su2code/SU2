/*!
 * \file CTransENVariable.hpp
 * \brief Declaration of the variables of the e^N transition model.
 * \author F. Palacios, T. Economon
 * \version 7.4.0 "Blackbird"
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

#include "CTurbVariable.hpp"

/*!
 * \class CTransENVariable
 * \brief e^N Transition model variables.
 * \ingroup Turbulence_Model
 * \author R. Roos
 */

class CTransENVariable final : public CTurbVariable {
protected:
  VectorType AmplificationFactor;
  VectorType ModifiedIntermittency;

  /*--- Debug terms ---*/
  VectorType Prod;
  VectorType Uedge;
  VectorType HL;
  VectorType H12;
  VectorType FG;
  VectorType FC;
  VectorType REY;
  VectorType REY0;
  VectorType Dist;
  VectorType Strain;

  VectorType normal_x;
  VectorType normal_y;
  VectorType normal_z;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] AmplificationFactor - Amplification factor(n) (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTransENVariable(su2double AmplificationFactor,su2double ModifiedIntermittency, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTransENVariable() override = default;

  /*!
   * \brief Set Amplification Factor.
   */
  void SetAmplificationFactor(unsigned long iPoint, su2double val_AmplificationFactor) override ;

  /*!
   * \brief Set Amplification Factor.
   */
  void SetModifiedIntermittency(unsigned long iPoint, su2double val_Gamma) override ;

  void SetNormal(unsigned long iPoint, su2double val_normal_x, su2double val_normal_y, su2double val_normal_z) override;

  void SetProd(unsigned long iPoint, su2double val_Prod) override;
  void SetUedge(unsigned long iPoint, su2double val_Uedge) override;
  void SetHL(unsigned long iPoint, su2double val_HL) override;
  void SetH12(unsigned long iPoint, su2double val_H12) override;
  void SetFG(unsigned long iPoint, su2double val_FG) override;
  void SetFC(unsigned long iPoint, su2double val_FC) override;
  void SetREY(unsigned long iPoint, su2double val_REY) override;
  void SetREY0(unsigned long iPoint, su2double val_REY0) override;
  void SetDist(unsigned long iPoint, su2double val_Dist) override;
  void SetStrain(unsigned long iPoint, su2double val_Strain) override;

  /*!
   * \brief Value of AmplificationFactor.
   */
  inline su2double GetAmplificationFactor(unsigned long iPoint) const override { return AmplificationFactor(iPoint); }

  /*!
   * \brief Value of AmplificationFactor.
   */
  inline su2double GetModifiedIntermittency(unsigned long iPoint) const override { return ModifiedIntermittency(iPoint); }


  inline su2double GetNormal_x(unsigned long iPoint) const override {return normal_x(iPoint);};
  inline su2double GetNormal_y(unsigned long iPoint) const override {return normal_y(iPoint);};
  inline su2double GetNormal_z(unsigned long iPoint) const override {return normal_z(iPoint);};

  inline su2double GetProd(unsigned long iPoint) const override { return Prod(iPoint); }
  inline su2double GetUedge(unsigned long iPoint) const override { return Uedge(iPoint); }
  inline su2double GetHL(unsigned long iPoint) const override { return HL(iPoint); }
  inline su2double GetH12(unsigned long iPoint) const override { return H12(iPoint); }
  inline su2double GetFG(unsigned long iPoint) const override { return FG(iPoint); }
  inline su2double GetFC(unsigned long iPoint) const override { return FC(iPoint); }
  inline su2double GetREY(unsigned long iPoint) const override { return REY(iPoint); }
  inline su2double GetREY0(unsigned long iPoint) const override { return REY0(iPoint); }
  inline su2double GetDist(unsigned long iPoint) const override { return Dist(iPoint); }
  inline su2double GetStrain(unsigned long iPoint) const override { return Strain(iPoint); }

};
