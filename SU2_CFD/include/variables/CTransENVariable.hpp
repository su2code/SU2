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
  VectorType Prod_n;
  VectorType Prod_g;
  VectorType Dest_g;
  VectorType GammaN;
  VectorType HL;
  VectorType H12;
  VectorType FG;
  VectorType FC;
  VectorType REV;
  VectorType REV0;
  VectorType Dist;
  VectorType Strain;
  VectorType Fonset1;
  VectorType Fonset;
  VectorType Fturb;

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

  void SetProdN(unsigned long iPoint, su2double val_ProdN) override;
  void SetProdG(unsigned long iPoint, su2double val_ProdG) override;
  void SetDestG(unsigned long iPoint, su2double val_ProdG) override;
  void SetGammaN(unsigned long iPoint, su2double val_GammaN) override;
  void SetHL(unsigned long iPoint, su2double val_HL) override;
  void SetH12(unsigned long iPoint, su2double val_H12) override;
  void SetFG(unsigned long iPoint, su2double val_FG) override;
  void SetFC(unsigned long iPoint, su2double val_FC) override;
  void SetREV(unsigned long iPoint, su2double val_REV) override;
  void SetREV0(unsigned long iPoint, su2double val_REV0) override;
  void SetDist(unsigned long iPoint, su2double val_Dist) override;
  void SetStrain(unsigned long iPoint, su2double val_Strain) override;
  void SetFonset1(unsigned long iPoint, su2double val_Fonset1) override;
  void SetFonset(unsigned long iPoint, su2double val_Fonset) override;
  void SetFturb(unsigned long iPoint, su2double val_Fturb) override;

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

  inline su2double GetProdN(unsigned long iPoint) const override { return Prod_n(iPoint); }
  inline su2double GetProdG(unsigned long iPoint) const override { return Prod_g(iPoint); }
  inline su2double GetDestG(unsigned long iPoint) const override { return Dest_g(iPoint); }
  inline su2double GetGammaN(unsigned long iPoint) const override { return GammaN(iPoint); }
  inline su2double GetHL(unsigned long iPoint) const override { return HL(iPoint); }
  inline su2double GetH12(unsigned long iPoint) const override { return H12(iPoint); }
  inline su2double GetFG(unsigned long iPoint) const override { return FG(iPoint); }
  inline su2double GetFC(unsigned long iPoint) const override { return FC(iPoint); }
  inline su2double GetREV(unsigned long iPoint) const override { return REV(iPoint); }
  inline su2double GetREV0(unsigned long iPoint) const override { return REV0(iPoint); }
  inline su2double GetDist(unsigned long iPoint) const override { return Dist(iPoint); }
  inline su2double GetStrain(unsigned long iPoint) const override { return Strain(iPoint); }
  inline su2double GetFonset1(unsigned long iPoint) const override { return Fonset1(iPoint); }
  inline su2double GetFonset(unsigned long iPoint) const override { return Fonset(iPoint); }
  inline su2double GetFturb(unsigned long iPoint) const override { return Fturb(iPoint); }

};
