/*!
 * \file turb_sources.cpp
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

#include "../../../include/numerics/flamelet/scalar_sources.hpp"

CSourcePieceWise_transportedScalar_general::CSourcePieceWise_transportedScalar_general(unsigned short val_nDim,
                                                   unsigned short val_nVar,
                                                   const CConfig* config) :
                          CNumerics(val_nDim, val_nVar, config) {

  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  axisymmetric = config->GetAxisymmetric();
  viscous = config->GetViscous();
  implicit = (config->GetKind_TimeIntScheme_Scalar() == EULER_IMPLICIT);

  /*--- "Allocate" the Jacobian using the static buffer. ---*/
  Jacobian_i[0] = Jacobian_Buffer;
  Jacobian_i[1] = Jacobian_Buffer+2;

}

CNumerics::ResidualType<> CSourcePieceWise_transportedScalar_general::ComputeResidual(const CConfig* config) {

  AD::StartPreacc();
  //AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(ScalarVar_i, nVar);
  AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);
  //AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F2_i); AD::SetPreaccIn(CDkw_i);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);
  //AD::SetPreaccIn(Vorticity_i, 3);

  unsigned short iDim;

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+6);

    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+4];
    Eddy_Viscosity_i = V_i[nDim+5];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7);

    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];
  }

  Residual[0] = 0.0;       Residual[1] = 0.0;
  Jacobian_i[0][0] = 0.0;  Jacobian_i[0][1] = 0.0;
  Jacobian_i[1][0] = 0.0;  Jacobian_i[1][1] = 0.0;

   /*--- Add the production terms to the residuals. ---*/

   //Residual[0] += pk*Volume;
   //Residual[1] += pw*Volume;

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
