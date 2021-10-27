/*!
 * \file scalarLegacy_sources.cpp
 * \brief Implementation of numerics classes for integration of
 *        turbulence source-terms.
 * \author F. Palacios, T. Economon
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/flamelet/scalarLegacy_sources.hpp"

CSourcePieceWise_transportedScalar_general::CSourcePieceWise_transportedScalar_general(unsigned short val_nDim,
                                                   unsigned short val_nVar,
                                                   const CConfig* config) :
                          CNumerics(val_nDim, val_nVar, config) {

  incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
  axisymmetric = config->GetAxisymmetric();
  viscous = config->GetViscous();
  implicit = (config->GetKind_TimeIntScheme_Scalar() == EULER_IMPLICIT);
  flame  = (config->GetKind_Scalar_Model() == PROGRESS_VARIABLE);
  inc_rans = (config->GetKind_Solver() == INC_RANS) || (config->GetKind_Solver() == DISC_ADJ_INC_RANS);

  Sc_t = config->GetSchmidt_Turb();

  Residual       = new su2double [nVar];
  scalar_sources = new su2double [nVar];
  Jacobian_i     = new su2double* [nVar];

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar] ();
  }
}

CSourcePieceWise_transportedScalar_general::~CSourcePieceWise_transportedScalar_general(void){
  delete [] Residual;
   if (Jacobian_i != nullptr) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete [] Jacobian_i[iVar];
    }
    delete [] Jacobian_i;
  }
}

CNumerics::ResidualType<> CSourcePieceWise_transportedScalar_general::ComputeResidual(const CConfig* config) {

  AD::StartPreacc();
  AD::SetPreaccIn(ScalarVar_i, nVar);
  AD::SetPreaccIn(scalar_sources, nVar);
  AD::SetPreaccIn(Volume); 

  // FIXME dan: the next two lines crashes when I run SU2_CFD_AD on the asym probe case
  // AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
  // AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);

  //if (config->GetKind_Scalar_Model() == PROGRESS_VARIABLE)
  //  AD::SetPreaccIn(sourcepv_i);

  //unsigned short iDim;

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+6);

    //Density_i = V_i[nDim+2];
    //Laminar_Viscosity_i = V_i[nDim+4];
    // we do not know if this exists
    //Eddy_Viscosity_i = V_i[nDim+5];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7);

    //Density_i = V_i[nDim+2];
    //Laminar_Viscosity_i = V_i[nDim+5];
    // we do not know if this exists
    //Eddy_Viscosity_i = V_i[nDim+6];
  }

  /*--- Implicit part for production term (to do). ---*/
  for (auto i_var = 0; i_var < nVar; i_var++) {
    Residual[i_var] = scalar_sources[i_var] * Volume;
    for (auto j_var = 0; j_var < nVar; j_var++) {
      Jacobian_i[i_var][j_var] = 0.0;
    }
  }
  // FIXME dan: add source term derivatives to jacobian
  //Jacobian[i][j] = dSource_i / dScalar_j

   /*--- Add the production terms to the residuals. ---*/

   /*--- Contribution due to 2D axisymmetric formulation ---*/
   
   if (axisymmetric) ResidualAxisymmetric();

   /*--- Implicit part ---*/

   //Jacobian_i[0][0] =0.0;// -beta_star*ScalarVar_i[1]*Volume;
   //Jacobian_i[0][1] = 0.0;//-beta_star*ScalarVar_i[0]*Volume;
   //Jacobian_i[1][0] = 0.0;
   //Jacobian_i[1][1] = 0.0;//-2.0*beta_blended*ScalarVar_i[1]*Volume;

  AD::SetPreaccOut(Residual, nVar);
  AD::EndPreacc();

  return ResidualType<>(Residual, Jacobian_i, nullptr);

}
