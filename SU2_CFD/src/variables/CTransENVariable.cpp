/*!
 * \file CTransENVariable.cpp
 * \brief Definition of the solution fields.
 * \author R. Roos
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


#include "../../include/variables/CTransENVariable.hpp"

CTransENVariable::CTransENVariable(su2double AmplificationFactor,su2double ModifiedIntermittency, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
  : CTurbVariable(npoint, ndim, nvar, config) {

  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
  {
    Solution(iPoint,0) = AmplificationFactor;
    Solution(iPoint,1) = ModifiedIntermittency;
  }

  Solution_Old = Solution;
  
  /* Normal values and velocity used for Boundary Pressure Gradient H_L */
  nAuxVar = 1;
  Grad_AuxVar.resize(nPoint, nAuxVar, nDim, 0.0);
  AuxVar.resize(nPoint, nAuxVar) = su2double(0.0);

  normal_x.resize(nPoint) = 0.0;
  normal_y.resize(nPoint) = 0.0;
  normal_z.resize(nPoint) = 0.0;

  /* List of debug variables */
  Prod_n.resize(nPoint) = 0.0;
  Prod_g.resize(nPoint) = 0.0;
  Dest_g.resize(nPoint) = 0.0;
  GammaN.resize(nPoint) = 0.0;
  HL.resize(nPoint) = 0.0;
  H12.resize(nPoint) = 0.0;
  FG.resize(nPoint) = 0.0;
  FC.resize(nPoint) = 0.0;
  REV.resize(nPoint) = 0.0;
  REV0.resize(nPoint) = 0.0;
  Dist.resize(nPoint) = 0.0;
  Strain.resize(nPoint) = 0.0;
  Fonset1.resize(nPoint) = 0.0;
  Fonset.resize(nPoint) = 0.0;
  Fturb.resize(nPoint) = 0.0;

}

void CTransENVariable::SetAmplificationFactor(unsigned long iPoint, su2double val_AmplificationFactor) {
  AmplificationFactor(iPoint) = val_AmplificationFactor;
}

void CTransENVariable::SetModifiedIntermittency(unsigned long iPoint, su2double val_Gamma) {
  ModifiedIntermittency(iPoint) = val_Gamma;
}

void CTransENVariable::SetNormal(unsigned long iPoint, su2double val_normal_x, su2double val_normal_y, su2double val_normal_z) {
  normal_x(iPoint) = val_normal_x;
  normal_y(iPoint) = val_normal_y;
  normal_z(iPoint) = val_normal_z;
}

void CTransENVariable::SetProdN(unsigned long iPoint, su2double val_ProdN) {
  Prod_n(iPoint) = val_ProdN;
}

void CTransENVariable::SetProdG(unsigned long iPoint, su2double val_ProdG) {
  Prod_g(iPoint) = val_ProdG;
}

void CTransENVariable::SetDestG(unsigned long iPoint, su2double val_DestG) {
  Dest_g(iPoint) = val_DestG;
}

void CTransENVariable::SetGammaN(unsigned long iPoint, su2double val_GammaN) {
  GammaN(iPoint) = val_GammaN;
}

void CTransENVariable::SetHL(unsigned long iPoint, su2double val_HL) {
  HL(iPoint) = val_HL;
}

void CTransENVariable::SetH12(unsigned long iPoint, su2double val_H12) {
  H12(iPoint) = val_H12;
}

void CTransENVariable::SetFG(unsigned long iPoint, su2double val_FG) {
  FG(iPoint) = val_FG;
}

void CTransENVariable::SetFC(unsigned long iPoint, su2double val_FC) {
  FC(iPoint) = val_FC;
}

void CTransENVariable::SetREV(unsigned long iPoint, su2double val_REV) {
  REV(iPoint) = val_REV;
}

void CTransENVariable::SetREV0(unsigned long iPoint, su2double val_REV0) {
  REV0(iPoint) = val_REV0;
}

void CTransENVariable::SetDist(unsigned long iPoint, su2double val_Dist) {
  Dist(iPoint) = val_Dist;
}

void CTransENVariable::SetStrain(unsigned long iPoint, su2double val_Strain) {
  Strain(iPoint) = val_Strain;
}

void CTransENVariable::SetFonset1(unsigned long iPoint, su2double val_Fonset1) {
  Fonset1(iPoint) = val_Fonset1;
}

void CTransENVariable::SetFonset(unsigned long iPoint, su2double val_Fonset) {
  Fonset(iPoint) = val_Fonset;
}

void CTransENVariable::SetFturb(unsigned long iPoint, su2double val_Fturb) {
  Fturb(iPoint) = val_Fturb;
}
