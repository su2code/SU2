/*!
 * \file hllc.cpp
 * \brief Implementations of HLLC schemes.
 * \author F. Palacios, T. Economon
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

#include "../../../../include/numerics/flow/convection/hllc.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwHLLC_Flow::CUpwHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  kappa = config->GetRoe_Kappa();
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  Gamma = config->GetGamma();

  Gamma_Minus_One = Gamma - 1.0;

  IntermediateState = new su2double [nVar];
  dSm_dU            = new su2double [nVar];
  dPI_dU            = new su2double [nVar];
  drhoStar_dU       = new su2double [nVar];
  dpStar_dU         = new su2double [nVar];
  dEStar_dU         = new su2double [nVar];

  Velocity_i        = new su2double [nDim];
  Velocity_j        = new su2double [nDim];
  RoeVelocity       = new su2double [nDim];

  Flux = new su2double [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CUpwHLLC_Flow::~CUpwHLLC_Flow() {

  delete [] IntermediateState;
  delete [] dSm_dU;
  delete [] dPI_dU;
  delete [] drhoStar_dU;
  delete [] dpStar_dU;
  delete [] dEStar_dU;

  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Flux;
  delete [] Jacobian_i;
  delete [] Jacobian_j;

}

CNumerics::ResidualType<> CUpwHLLC_Flow::ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Face area (norm or the normal vector) ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

  /*-- Unit Normal ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim] / Area;

  /*-- Fluid velocity at node i,j ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim]  = V_i[iDim+1];
    Velocity_j[iDim]  = V_j[iDim+1];
  }

  /*--- Primitive variables at point i ---*/

  Pressure_i = V_i[nDim+1];
  Density_i  = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];

  /*--- Primitive variables at point j ---*/

  Pressure_j = V_j[nDim+1];
  Density_j  = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];


  sq_vel_i = 0.0;
  sq_vel_j = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    sq_vel_i += Velocity_i[iDim] * Velocity_i[iDim];
    sq_vel_j += Velocity_j[iDim] * Velocity_j[iDim];
  }

  Energy_i = Enthalpy_i - Pressure_i / Density_i;
  Energy_j = Enthalpy_j - Pressure_j / Density_j;

  SoundSpeed_i = sqrt( (Enthalpy_i - 0.5 * sq_vel_i) * Gamma_Minus_One );
  SoundSpeed_j = sqrt( (Enthalpy_j - 0.5 * sq_vel_j) * Gamma_Minus_One );

  /*--- Projected velocities ---*/

  ProjVelocity_i = 0;
  ProjVelocity_j = 0;

  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim] * UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim] * UnitNormal[iDim];
  }

  /*--- Projected Grid Velocity ---*/

  ProjInterfaceVel = 0;

  if (dynamic_grid) {

    for (iDim = 0; iDim < nDim; iDim++)
      ProjInterfaceVel += 0.5 * ( GridVel_i[iDim] + GridVel_j[iDim] )*UnitNormal[iDim];

    SoundSpeed_i -= ProjInterfaceVel;
    SoundSpeed_j += ProjInterfaceVel;

    ProjVelocity_i -= ProjInterfaceVel;
    ProjVelocity_j -= ProjInterfaceVel;
  }

  /*--- Roe's averaging ---*/

  Rrho = ( sqrt(Density_i) + sqrt(Density_j) );

  sq_velRoe        = 0.0;
  RoeProjVelocity  = - ProjInterfaceVel;

  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = ( Velocity_i[iDim] * sqrt(Density_i) + Velocity_j[iDim] * sqrt(Density_j) ) / Rrho;
    sq_velRoe        +=  RoeVelocity[iDim] * RoeVelocity[iDim];
    RoeProjVelocity  +=  RoeVelocity[iDim] * UnitNormal[iDim];
  }

  /*--- Mean Roe variables iPoint and jPoint ---*/

  RoeDensity = sqrt( Density_i * Density_j );
  RoeEnthalpy = ( sqrt(Density_j) * Enthalpy_j + sqrt(Density_i) * Enthalpy_i) / Rrho;

  /*--- Roe-averaged speed of sound ---*/

  //RoeSoundSpeed2 = Gamma_Minus_One * ( RoeEnthalpy - 0.5 * sq_velRoe );
  RoeSoundSpeed  = sqrt( Gamma_Minus_One * ( RoeEnthalpy - 0.5 * sq_velRoe  ) ) - ProjInterfaceVel;


  /*--- Speed of sound at L and R ---*/

  sL = min( RoeProjVelocity - RoeSoundSpeed, ProjVelocity_i - SoundSpeed_i);
  sR = max( RoeProjVelocity + RoeSoundSpeed, ProjVelocity_j + SoundSpeed_j);

  /*--- speed of contact surface ---*/

  RHO = Density_j * (sR - ProjVelocity_j) - Density_i * (sL - ProjVelocity_i);
  sM = ( Pressure_i - Pressure_j - Density_i * ProjVelocity_i * ( sL - ProjVelocity_i ) + Density_j * ProjVelocity_j * ( sR - ProjVelocity_j ) ) / RHO;

  /*--- Pressure at right and left (Pressure_j=Pressure_i) side of contact surface ---*/

  pStar = Density_j * ( ProjVelocity_j - sR ) * ( ProjVelocity_j - sM ) + Pressure_j;


  if (sM > 0.0) {

    if (sL > 0.0) {

      /*--- Compute Left Flux ---*/

      Flux[0] = Density_i * ProjVelocity_i;
      for (iDim = 0; iDim < nDim; iDim++)
        Flux[iDim+1] = Density_i * Velocity_i[iDim] * ProjVelocity_i + Pressure_i * UnitNormal[iDim];
      Flux[nVar-1] = Enthalpy_i * Density_i * ProjVelocity_i;

    }
    else {

      /*--- Compute Flux Left Star from Left Star State ---*/

      rhoSL = ( sL - ProjVelocity_i ) / ( sL - sM );

      IntermediateState[0] = rhoSL * Density_i;
      for (iDim = 0; iDim < nDim; iDim++)
        IntermediateState[iDim+1] = rhoSL * ( Density_i * Velocity_i[iDim] + ( pStar - Pressure_i ) / ( sL - ProjVelocity_i ) * UnitNormal[iDim] ) ;
      IntermediateState[nVar-1] = rhoSL * ( Density_i * Energy_i - ( Pressure_i * ProjVelocity_i - pStar * sM) / ( sL - ProjVelocity_i ) );


      Flux[0] = sM * IntermediateState[0];
      for (iDim = 0; iDim < nDim; iDim++)
        Flux[iDim+1] = sM * IntermediateState[iDim+1] + pStar * UnitNormal[iDim];
      Flux[nVar-1] = sM * ( IntermediateState[nVar-1] + pStar ) + pStar * ProjInterfaceVel;
    }
  }
  else {

    if (sR < 0.0) {

      /*--- Compute Right Flux ---*/

      Flux[0] = Density_j * ProjVelocity_j;
      for (iDim = 0; iDim < nDim; iDim++)
        Flux[iDim+1] = Density_j * Velocity_j[iDim] * ProjVelocity_j + Pressure_j * UnitNormal[iDim];
      Flux[nVar-1] = Enthalpy_j * Density_j * ProjVelocity_j;
    }
    else {

      /*--- Compute Flux Right Star from Right Star State ---*/

      rhoSR = ( sR - ProjVelocity_j ) / ( sR - sM );

      IntermediateState[0] = rhoSR * Density_j;
      for (iDim = 0; iDim < nDim; iDim++)
        IntermediateState[iDim+1] = rhoSR * ( Density_j * Velocity_j[iDim] + ( pStar - Pressure_j ) / ( sR - ProjVelocity_j ) * UnitNormal[iDim] ) ;
      IntermediateState[nVar-1] = rhoSR * ( Density_j * Energy_j - ( Pressure_j * ProjVelocity_j - pStar * sM ) / ( sR - ProjVelocity_j ) );


      Flux[0] = sM * IntermediateState[0];
      for (iDim = 0; iDim < nDim; iDim++)
        Flux[iDim+1] = sM * IntermediateState[iDim+1] + pStar * UnitNormal[iDim];
      Flux[nVar-1] = sM * (IntermediateState[nVar-1] + pStar ) + pStar * ProjInterfaceVel;
    }
  }


  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] *= Area;

  /*--- Return early if the Jacobians do not need to be computed. ---*/

  if (implicit)
  {
  if (sM > 0.0) {

    if (sL > 0.0) {

      /*--- Compute Jacobian based on Left State ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          Jacobian_j[iVar][jVar] = 0;

      GetInviscidProjJac(Velocity_i, &Energy_i, UnitNormal, 1.0, Jacobian_i);

    }
    else {
      /*--- Compute Jacobian based on Left Star State ---*/

      EStar = IntermediateState[nVar-1];
      Omega = 1/(sL-sM);
      OmegaSM = Omega * sM;


      /*--------- Left Jacobian ---------*/


      /*--- Computing pressure derivatives d/dU_L (PI) ---*/

      dPI_dU[0] = 0.5 * Gamma_Minus_One * sq_vel_i;
      for (iDim = 0; iDim < nDim; iDim++)
        dPI_dU[iDim+1] = - Gamma_Minus_One * Velocity_i[iDim];
      dPI_dU[nVar-1] = Gamma_Minus_One;


      /*--- Computing d/dU_L (Sm) ---*/

      dSm_dU[0] = ( - ProjVelocity_i * ProjVelocity_i + sM * sL + dPI_dU[0] ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = ( UnitNormal[iDim] * ( 2 * ProjVelocity_i - sL - sM ) + dPI_dU[iDim+1] ) / RHO;
      dSm_dU[nVar-1] = dPI_dU[nVar-1] / RHO;


      /*--- Computing d/dU_L (rhoStar) ---*/

      drhoStar_dU[0] = Omega * ( sL + IntermediateState[0] * dSm_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        drhoStar_dU[iDim+1] = Omega * ( - UnitNormal[iDim] + IntermediateState[0] * dSm_dU[iDim+1] );
      drhoStar_dU[nVar-1] = Omega * IntermediateState[0] * dSm_dU[nVar-1];


      /*--- Computing d/dU_L (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_i * (sR - ProjVelocity_j) * dSm_dU[iVar];


      /*--- Computing d/dU_L (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );

      dEStar_dU[0] += Omega * ProjVelocity_i * ( Enthalpy_i - dPI_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        dEStar_dU[iDim+1] += Omega * ( - UnitNormal[iDim] * Enthalpy_i - ProjVelocity_i * dPI_dU[iDim+1] );
      dEStar_dU[nVar-1] += Omega * ( sL - ProjVelocity_i - ProjVelocity_i * dPI_dU[nVar-1] );



      /*--- Jacobian First Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_i[0][iVar] = sM * drhoStar_dU[iVar] + IntermediateState[0] * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (jDim = 0; jDim < nDim; jDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Jacobian_i[jDim+1][iVar] = ( OmegaSM + 1 ) * ( UnitNormal[jDim] * dpStar_dU[iVar] + IntermediateState[jDim+1] * dSm_dU[iVar] );

        Jacobian_i[jDim+1][0] += OmegaSM * Velocity_i[jDim] * ProjVelocity_i;

        Jacobian_i[jDim+1][jDim+1] += OmegaSM * (sL - ProjVelocity_i);

        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_i[jDim+1][iDim+1] -= OmegaSM * Velocity_i[jDim] * UnitNormal[iDim];

        for (iVar = 0; iVar < nVar; iVar++)
          Jacobian_i[jDim+1][iVar] -= OmegaSM * dPI_dU[iVar] * UnitNormal[jDim];
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_i[nVar-1][iVar] = sM * ( dEStar_dU[iVar] + dpStar_dU[iVar] ) + ( EStar + pStar ) * dSm_dU[iVar];


      /*--------- Right Jacobian ---------*/


      /*--- Computing d/dU_R (Sm) ---*/

      dSm_dU[0] = ( ProjVelocity_j * ProjVelocity_j - sM * sR - 0.5 * Gamma_Minus_One * sq_vel_j ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = - ( UnitNormal[iDim] * ( 2 * ProjVelocity_j - sR - sM) - Gamma_Minus_One * Velocity_j[iDim] ) / RHO;
      dSm_dU[nVar-1]  = - Gamma_Minus_One / RHO;


      /*--- Computing d/dU_R (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_j * (sL - ProjVelocity_i) * dSm_dU[iVar];


      /*--- Computing d/dU_R (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );



      /*--- Jacobian First Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_j[0][iVar] = IntermediateState[0] * ( OmegaSM + 1 ) * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Jacobian_j[iDim+1][iVar] = ( OmegaSM + 1 ) * ( IntermediateState[iDim+1] * dSm_dU[iVar] + UnitNormal[iDim] * dpStar_dU[iVar] );
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_j[nVar-1][iVar] = sM * (dEStar_dU[iVar] + dpStar_dU[iVar]) + (EStar + pStar) * dSm_dU[iVar];
    }
  }
  else {
    if (sR < 0.0) {

      /*--- Compute Jacobian based on Right State ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0;

      GetInviscidProjJac(Velocity_j, &Energy_j, UnitNormal, 1.0, Jacobian_j);

    }
    else {
      /*--- Compute Jacobian based on Right Star State ---*/

      EStar = IntermediateState[nVar-1];
      Omega = 1/(sR-sM);
      OmegaSM = Omega * sM;


      /*--------- Left Jacobian ---------*/


      /*--- Computing d/dU_L (Sm) ---*/

      dSm_dU[0] = ( - ProjVelocity_i * ProjVelocity_i + sM * sL + 0.5 * Gamma_Minus_One * sq_vel_i ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = ( UnitNormal[iDim] * ( 2 * ProjVelocity_i - sL - sM ) - Gamma_Minus_One * Velocity_i[iDim] ) / RHO;
      dSm_dU[nVar-1] = Gamma_Minus_One / RHO;


      /*--- Computing d/dU_L (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_i * (sR - ProjVelocity_j) * dSm_dU[iVar];


      /*--- Computing d/dU_L (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );



      /*--- Jacobian First Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_i[0][iVar] = IntermediateState[0] * ( OmegaSM + 1 ) * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Jacobian_i[iDim+1][iVar] = (OmegaSM + 1) * ( IntermediateState[iDim+1] * dSm_dU[iVar] + UnitNormal[iDim] * dpStar_dU[iVar] );
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_i[nVar-1][iVar] = sM * (dEStar_dU[iVar] + dpStar_dU[iVar]) + (EStar + pStar) * dSm_dU[iVar];



      /*--------- Right Jacobian ---------*/


      /*--- Computing pressure derivatives d/dU_R (PI) ---*/

      dPI_dU[0] = 0.5 * Gamma_Minus_One * sq_vel_j;
      for (iDim = 0; iDim < nDim; iDim++)
        dPI_dU[iDim+1] = - Gamma_Minus_One * Velocity_j[iDim];
      dPI_dU[nVar-1] = Gamma_Minus_One;



      /*--- Computing d/dU_R (Sm) ---*/

      dSm_dU[0] = - ( - ProjVelocity_j * ProjVelocity_j + sM * sR + dPI_dU[0] ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = - ( UnitNormal[iDim] * ( 2 * ProjVelocity_j - sR - sM) + dPI_dU[iDim+1] ) / RHO;
      dSm_dU[nVar-1]  = - dPI_dU[nVar-1] / RHO;


      /*--- Computing d/dU_R (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_j * (sL - ProjVelocity_i) * dSm_dU[iVar];


      /*--- Computing d/dU_R (rhoStar) ---*/

      drhoStar_dU[0] = Omega * ( sR + IntermediateState[0] * dSm_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        drhoStar_dU[iDim+1] = Omega * ( - UnitNormal[iDim] + IntermediateState[0] * dSm_dU[iDim+1] );
      drhoStar_dU[nVar-1] = Omega * IntermediateState[0] * dSm_dU[nVar-1];


      /*--- Computing d/dU_R (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );

      dEStar_dU[0] += Omega * ProjVelocity_j * ( Enthalpy_j - dPI_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        dEStar_dU[iDim+1] += Omega * ( - UnitNormal[iDim] * Enthalpy_j - ProjVelocity_j * dPI_dU[iDim+1] );
      dEStar_dU[nVar-1] += Omega * ( sR - ProjVelocity_j - ProjVelocity_j * dPI_dU[nVar-1] );



      /*--- Jacobian First Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_j[0][iVar] = sM * drhoStar_dU[iVar] + IntermediateState[0] * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (jDim = 0; jDim < nDim; jDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Jacobian_j[jDim+1][iVar] = ( OmegaSM + 1 ) * ( UnitNormal[jDim] * dpStar_dU[iVar] + IntermediateState[jDim+1] * dSm_dU[iVar] );

        Jacobian_j[jDim+1][0] += OmegaSM * Velocity_j[jDim] * ProjVelocity_j;

        Jacobian_j[jDim+1][jDim+1] += OmegaSM * (sR - ProjVelocity_j);

        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_j[jDim+1][iDim+1] -= OmegaSM * Velocity_j[jDim] * UnitNormal[iDim];

        for (iVar = 0; iVar < nVar; iVar++)
          Jacobian_j[jDim+1][iVar] -= OmegaSM * dPI_dU[iVar] * UnitNormal[jDim];
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_j[nVar-1][iVar] = sM * ( dEStar_dU[iVar] + dpStar_dU[iVar] ) + ( EStar + pStar ) * dSm_dU[iVar];

    }
  }


  /*--- Jacobians of the inviscid flux, scale = k because Flux ~ 0.5*(fc_i+fc_j)*Normal ---*/

  Area *= kappa;

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      Jacobian_i[iVar][jVar] *=   Area;
      Jacobian_j[iVar][jVar] *=   Area;
    }
  }
  } // end if implicit

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

CUpwGeneralHLLC_Flow::CUpwGeneralHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  kappa = config->GetRoe_Kappa();
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  Gamma = config->GetGamma();

  IntermediateState = new su2double [nVar];
  dSm_dU            = new su2double [nVar];
  dPI_dU            = new su2double [nVar];
  drhoStar_dU       = new su2double [nVar];
  dpStar_dU         = new su2double [nVar];
  dEStar_dU         = new su2double [nVar];

  Velocity_i        = new su2double [nDim];
  Velocity_j        = new su2double [nDim];
  RoeVelocity       = new su2double [nDim];

  Flux = new su2double [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CUpwGeneralHLLC_Flow::~CUpwGeneralHLLC_Flow() {

  delete [] IntermediateState;
  delete [] dSm_dU;
  delete [] dPI_dU;
  delete [] drhoStar_dU;
  delete [] dpStar_dU;
  delete [] dEStar_dU;

  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Flux;
  delete [] Jacobian_i;
  delete [] Jacobian_j;

}

CNumerics::ResidualType<> CUpwGeneralHLLC_Flow::ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Face area (norm or the normal vector) ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

  /*-- Unit Normal ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim] / Area;

  /*-- Fluid velocity at node i,j ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim]  = V_i[iDim+1];
    Velocity_j[iDim]  = V_j[iDim+1];
  }

  /*--- Primitive variables at point i ---*/

  Pressure_i = V_i[nDim+1];
  Density_i  = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];

  /*--- Primitive variables at point j ---*/

  Pressure_j = V_j[nDim+1];
  Density_j  = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];


  sq_vel_i = 0.0;
  sq_vel_j = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    sq_vel_i += Velocity_i[iDim] * Velocity_i[iDim];
    sq_vel_j += Velocity_j[iDim] * Velocity_j[iDim];
  }

  Energy_i         = Enthalpy_i - Pressure_i / Density_i;
  StaticEnthalpy_i = Enthalpy_i - 0.5 * sq_vel_i;
  StaticEnergy_i   = Energy_i - 0.5 * sq_vel_i;

  Kappa_i = S_i[1] / Density_i;
  Chi_i   = S_i[0] - Kappa_i * StaticEnergy_i;
  SoundSpeed_i = sqrt(Chi_i + StaticEnthalpy_i * Kappa_i);


  Energy_j         = Enthalpy_j - Pressure_j / Density_j;
  StaticEnthalpy_j = Enthalpy_j - 0.5 * sq_vel_j;
  StaticEnergy_j   = Energy_j - 0.5 * sq_vel_j;

  Kappa_j = S_j[1] / Density_j;
  Chi_j   = S_j[0] - Kappa_j * StaticEnergy_j;
  SoundSpeed_j = sqrt(Chi_j + StaticEnthalpy_j * Kappa_j);

  /*--- Projected velocities ---*/

  ProjVelocity_i = 0.0;
  ProjVelocity_j = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim] * UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim] * UnitNormal[iDim];
  }


  /*--- Projected Grid Velocity ---*/

  ProjInterfaceVel = 0;

  if (dynamic_grid) {

    for (iDim = 0; iDim < nDim; iDim++)
      ProjInterfaceVel += 0.5 * ( GridVel_i[iDim] + GridVel_j[iDim] )*UnitNormal[iDim];

    SoundSpeed_i -= ProjInterfaceVel;
    SoundSpeed_j += ProjInterfaceVel;

    ProjVelocity_i -= ProjInterfaceVel;
    ProjVelocity_j -= ProjInterfaceVel;
  }

  /*--- Roe's averaging ---*/

  Rrho = ( sqrt(Density_i) + sqrt(Density_j) );

  sq_velRoe        = 0.0;
  RoeProjVelocity  = - ProjInterfaceVel;

  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = ( Velocity_i[iDim] * sqrt(Density_i) + Velocity_j[iDim] * sqrt(Density_j) ) / Rrho;
    sq_velRoe        +=  RoeVelocity[iDim] * RoeVelocity[iDim];
    RoeProjVelocity  +=  RoeVelocity[iDim] * UnitNormal[iDim];
  }

  /*--- Mean Roe variables iPoint and jPoint ---*/

  RoeKappa = 0.5 * ( Kappa_i + Kappa_j );
  RoeChi   = 0.5 * ( Chi_i + Chi_j );
  RoeDensity = sqrt( Density_i * Density_j );
  RoeEnthalpy = ( sqrt(Density_j) * Enthalpy_j + sqrt(Density_i) * Enthalpy_i) / Rrho;

  VinokurMontagne();

  /*--- Roe-averaged speed of sound ---*/

  //RoeSoundSpeed2 = RoeChi + RoeKappa * ( RoeEnthalpy - 0.5 * sq_velRoe );
  RoeSoundSpeed  = sqrt( RoeChi + RoeKappa * ( RoeEnthalpy - 0.5 * sq_velRoe ) ) - ProjInterfaceVel;

  /*--- Speed of sound at L and R ---*/

  sL = min( RoeProjVelocity - RoeSoundSpeed, ProjVelocity_i - SoundSpeed_i );
  sR = max( RoeProjVelocity + RoeSoundSpeed, ProjVelocity_j + SoundSpeed_j );

  /*--- speed of contact surface ---*/

  RHO = Density_j * (sR - ProjVelocity_j) - Density_i * (sL - ProjVelocity_i);
  sM = ( Pressure_i - Pressure_j - Density_i * ProjVelocity_i * ( sL - ProjVelocity_i ) + Density_j * ProjVelocity_j * ( sR - ProjVelocity_j ) ) / RHO;

  /*--- Pressure at right and left (Pressure_j=Pressure_i) side of contact surface ---*/

  pStar = Density_j * ( ProjVelocity_j - sR ) * ( ProjVelocity_j - sM ) + Pressure_j;


  if (sM > 0.0) {

    if (sL > 0.0) {

      /*--- Compute Left Flux ---*/

      Flux[0] = Density_i * ProjVelocity_i;
      for (iDim = 0; iDim < nDim; iDim++)
        Flux[iDim+1] = Density_i * Velocity_i[iDim] * ProjVelocity_i + Pressure_i * UnitNormal[iDim];
      Flux[nVar-1] = Enthalpy_i * Density_i * ProjVelocity_i;
    }
    else {

      /*--- Compute Flux Left Star from Left Star State ---*/

      rhoSL = ( sL - ProjVelocity_i ) / ( sL - sM );

      IntermediateState[0] = rhoSL * Density_i;
      for (iDim = 0; iDim < nDim; iDim++)
        IntermediateState[iDim+1] = rhoSL * ( Density_i * Velocity_i[iDim] + ( pStar - Pressure_i ) / ( sL - ProjVelocity_i ) * UnitNormal[iDim] ) ;
      IntermediateState[nVar-1] = rhoSL * ( Density_i * Energy_i - ( Pressure_i * ProjVelocity_i - pStar * sM) / ( sL - ProjVelocity_i ) );


      Flux[0] = sM * IntermediateState[0];
      for (iDim = 0; iDim < nDim; iDim++)
        Flux[iDim+1] = sM * IntermediateState[iDim+1] + pStar * UnitNormal[iDim];
      Flux[nVar-1] = sM * ( IntermediateState[nVar-1] + pStar )  + pStar * ProjInterfaceVel;
    }
  }
  else {

    if (sR < 0.0) {

      /*--- Compute Right Flux ---*/

      Flux[0] = Density_j * ProjVelocity_j;
      for (iDim = 0; iDim < nDim; iDim++)
        Flux[iDim+1] = Density_j * Velocity_j[iDim] * ProjVelocity_j + Pressure_j * UnitNormal[iDim];
      Flux[nVar-1] = Enthalpy_j * Density_j * ProjVelocity_j;
    }
    else {

      /*--- Compute Flux Right Star from Right Star State ---*/

      rhoSR = ( sR - ProjVelocity_j ) / ( sR - sM );

      IntermediateState[0] = rhoSR * Density_j;
      for (iDim = 0; iDim < nDim; iDim++)
        IntermediateState[iDim+1] = rhoSR * ( Density_j * Velocity_j[iDim] + ( pStar - Pressure_j ) / ( sR - ProjVelocity_j ) * UnitNormal[iDim] ) ;
      IntermediateState[nVar-1] = rhoSR * ( Density_j * Energy_j - ( Pressure_j * ProjVelocity_j - pStar * sM ) / ( sR - ProjVelocity_j ) );


      Flux[0] = sM * IntermediateState[0];
      for (iDim = 0; iDim < nDim; iDim++)
        Flux[iDim+1] = sM * IntermediateState[iDim+1] + pStar * UnitNormal[iDim];
      Flux[nVar-1] = sM * (IntermediateState[nVar-1] + pStar )  + pStar * ProjInterfaceVel;
    }
  }

  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] *= Area;

  /*--- Return early if the Jacobians do not need to be computed. ---*/

  if (implicit)
  {
  if (sM > 0.0) {

    if (sL > 0.0) {

      /*--- Compute Jacobian based on Left State ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          Jacobian_j[iVar][jVar] = 0;


      GetInviscidProjJac(Velocity_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, 1.0, Jacobian_i);

    }
    else {
      /*--- Compute Jacobian based on Left Star State ---*/

      EStar = IntermediateState[nVar-1];
      Omega = 1/(sL-sM);
      OmegaSM = Omega * sM;


      /*--------- Left Jacobian ---------*/


      /*--- Computing pressure derivatives d/dU_L (PI) ---*/

      dPI_dU[0] = Chi_i - 0.5 * Kappa_i * sq_vel_i;
      for (iDim = 0; iDim < nDim; iDim++)
        dPI_dU[iDim+1] = - Kappa_i * Velocity_i[iDim];
      dPI_dU[nVar-1] = Kappa_i;


      /*--- Computing d/dU_L (Sm) ---*/

      dSm_dU[0] = ( - ProjVelocity_i * ProjVelocity_i + sM * sL + dPI_dU[0] ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = ( UnitNormal[iDim] * ( 2 * ProjVelocity_i - sL - sM ) + dPI_dU[iDim+1] ) / RHO;
      dSm_dU[nVar-1] = dPI_dU[nVar-1] / RHO;


      /*--- Computing d/dU_L (rhoStar) ---*/

      drhoStar_dU[0] = Omega * ( sL + IntermediateState[0] * dSm_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        drhoStar_dU[iDim+1] = Omega * ( - UnitNormal[iDim] + IntermediateState[0] * dSm_dU[iDim+1] );
      drhoStar_dU[nVar-1] = Omega * IntermediateState[0] * dSm_dU[nVar-1];


      /*--- Computing d/dU_L (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_i * (sR - ProjVelocity_j) * dSm_dU[iVar];


      /*--- Computing d/dU_L (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );

      dEStar_dU[0] += Omega * ProjVelocity_i * ( Enthalpy_i - dPI_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        dEStar_dU[iDim+1] += Omega * ( - UnitNormal[iDim] * Enthalpy_i - ProjVelocity_i * dPI_dU[iDim+1] );
      dEStar_dU[nVar-1] += Omega * ( sL - ProjVelocity_i - ProjVelocity_i * dPI_dU[nVar-1] );



      /*--- Jacobian First Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_i[0][iVar] = sM * drhoStar_dU[iVar] + IntermediateState[0] * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (jDim = 0; jDim < nDim; jDim++) {

        for (iVar = 0; iVar < nVar; iVar++)
          Jacobian_i[jDim+1][iVar] = ( OmegaSM + 1 ) * ( UnitNormal[jDim] * dpStar_dU[iVar] + IntermediateState[jDim+1] * dSm_dU[iVar] );

        Jacobian_i[jDim+1][0] += OmegaSM * Velocity_i[jDim] * ProjVelocity_i;

        Jacobian_i[jDim+1][jDim+1] += OmegaSM * (sL - ProjVelocity_i);

        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_i[jDim+1][iDim+1] -= OmegaSM * Velocity_i[jDim] * UnitNormal[iDim];

        for (iVar = 0; iVar < nVar; iVar++)
          Jacobian_i[jDim+1][iVar] -= OmegaSM * dPI_dU[iVar] * UnitNormal[jDim];
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_i[nVar-1][iVar] = sM * ( dEStar_dU[iVar] + dpStar_dU[iVar] ) + ( EStar + pStar ) * dSm_dU[iVar];


      /*--------- Right Jacobian ---------*/


      /*--- Computing pressure derivatives d/dU_R (PI) ---*/

      dPI_dU[0] = Chi_j - 0.5 * Kappa_j * sq_vel_j;
      for (iDim = 0; iDim < nDim; iDim++)
        dPI_dU[iDim+1] = - Kappa_j * Velocity_j[iDim];
      dPI_dU[nVar-1] = Kappa_j;


      /*--- Computing d/dU_R (Sm) ---*/

      dSm_dU[0] = ( ProjVelocity_j * ProjVelocity_j - sM * sR - dPI_dU[0] ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = - ( UnitNormal[iDim] * ( 2 * ProjVelocity_j - sR - sM) + dPI_dU[iDim+1] ) / RHO;
      dSm_dU[nVar-1]  = - dPI_dU[nVar-1] / RHO;


      /*--- Computing d/dU_R (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_j * (sL - ProjVelocity_i) * dSm_dU[iVar];


      /*--- Computing d/dU_R (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );



      /*--- Jacobian First Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_j[0][iVar] = IntermediateState[0] * ( OmegaSM + 1 ) * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Jacobian_j[iDim+1][iVar] = ( OmegaSM + 1 ) * ( IntermediateState[iDim+1] * dSm_dU[iVar] + UnitNormal[iDim] * dpStar_dU[iVar] );
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_j[nVar-1][iVar] = sM * (dEStar_dU[iVar] + dpStar_dU[iVar]) + (EStar + pStar) * dSm_dU[iVar];
    }
  }
  else {
    if (sR < 0.0) {

      /*--- Compute Jacobian based on Right State ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0;

      GetInviscidProjJac(Velocity_j, &Enthalpy_j, &Chi_j, &Kappa_j, UnitNormal, 1.0, Jacobian_j);

    }
    else {
      /*--- Compute Jacobian based on Right Star State ---*/

      EStar = IntermediateState[nVar-1];
      Omega = 1/(sR-sM);
      OmegaSM = Omega * sM;


      /*--------- Left Jacobian ---------*/


      /*--- Computing pressure derivatives d/dU_L (PI) ---*/

      dPI_dU[0] = Chi_i - 0.5 * Kappa_i * sq_vel_i;
      for (iDim = 0; iDim < nDim; iDim++)
        dPI_dU[iDim+1] = - Kappa_i * Velocity_i[iDim];
      dPI_dU[nVar-1] = Kappa_i;


      /*--- Computing d/dU_L (Sm) ---*/

      dSm_dU[0] = ( - ProjVelocity_i * ProjVelocity_i + sM * sL + dPI_dU[0] ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = ( UnitNormal[iDim] * ( 2 * ProjVelocity_i - sL - sM ) + dPI_dU[iDim+1] ) / RHO;
      dSm_dU[nVar-1] = dPI_dU[nVar-1] / RHO;


      /*--- Computing d/dU_L (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_i * (sR - ProjVelocity_j) * dSm_dU[iVar];


      /*--- Computing d/dU_L (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );



      /*--- Jacobian First Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_i[0][iVar] = IntermediateState[0] * ( OmegaSM + 1 ) * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Jacobian_i[iDim+1][iVar] = (OmegaSM + 1) * ( IntermediateState[iDim+1] * dSm_dU[iVar] + UnitNormal[iDim] * dpStar_dU[iVar] );
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_i[nVar-1][iVar] = sM * (dEStar_dU[iVar] + dpStar_dU[iVar]) + (EStar + pStar) * dSm_dU[iVar];



      /*--------- Right Jacobian ---------*/


      /*--- Computing pressure derivatives d/dU_R (PI) ---*/

      dPI_dU[0] = Chi_j - 0.5 * Kappa_j * sq_vel_j;
      for (iDim = 0; iDim < nDim; iDim++)
        dPI_dU[iDim+1] = - Kappa_j * Velocity_j[iDim];
      dPI_dU[nVar-1] = Kappa_j;


      /*--- Computing d/dU_R (Sm) ---*/

      dSm_dU[0] = - ( - ProjVelocity_j * ProjVelocity_j + sM * sR + dPI_dU[0] ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = - ( UnitNormal[iDim] * ( 2 * ProjVelocity_j - sR - sM) + dPI_dU[iDim+1] ) / RHO;
      dSm_dU[nVar-1]  = - dPI_dU[nVar-1] / RHO;


      /*--- Computing d/dU_R (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_j * (sL - ProjVelocity_i) * dSm_dU[iVar];


      /*--- Computing d/dU_R (rhoStar) ---*/

      drhoStar_dU[0] = Omega * ( sR + IntermediateState[0] * dSm_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        drhoStar_dU[iDim+1] = Omega * ( - UnitNormal[iDim] + IntermediateState[0] * dSm_dU[iDim+1] );
      drhoStar_dU[nVar-1] = Omega * IntermediateState[0] * dSm_dU[nVar-1];


      /*--- Computing d/dU_R (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );

      dEStar_dU[0] += Omega * ProjVelocity_j * ( Enthalpy_j - dPI_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        dEStar_dU[iDim+1] += Omega * ( - UnitNormal[iDim] * Enthalpy_j - ProjVelocity_j * dPI_dU[iDim+1] );
      dEStar_dU[nVar-1] += Omega * ( sR - ProjVelocity_j - ProjVelocity_j * dPI_dU[nVar-1] );



      /*--- Jacobian First Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_j[0][iVar] = sM * drhoStar_dU[iVar] + IntermediateState[0] * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (jDim = 0; jDim < nDim; jDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Jacobian_j[jDim+1][iVar] = ( OmegaSM + 1 ) * ( UnitNormal[jDim] * dpStar_dU[iVar] + IntermediateState[jDim+1] * dSm_dU[iVar] );

        Jacobian_j[jDim+1][0] += OmegaSM * Velocity_j[jDim] * ProjVelocity_j;

        Jacobian_j[jDim+1][jDim+1] += OmegaSM * (sR - ProjVelocity_j);

        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_j[jDim+1][iDim+1] -= OmegaSM * Velocity_j[jDim] * UnitNormal[iDim];

        for (iVar = 0; iVar < nVar; iVar++)
          Jacobian_j[jDim+1][iVar] -= OmegaSM * dPI_dU[iVar] * UnitNormal[jDim];
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_j[nVar-1][iVar] = sM * ( dEStar_dU[iVar] + dpStar_dU[iVar] ) + ( EStar + pStar ) * dSm_dU[iVar];
    }
  }


  /*--- Jacobians of the inviscid flux, scale = kappa because Flux ~ 0.5*(fc_i+fc_j)*Normal ---*/

  Area *= kappa;

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      Jacobian_i[iVar][jVar] *= Area;
      Jacobian_j[iVar][jVar] *= Area;
    }
  }
  } // end if implicit

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

void CUpwGeneralHLLC_Flow::VinokurMontagne() {

  su2double delta_rhoStaticEnergy, delta_rho, delta_p, err_P, s, D;

  delta_rho = Density_j - Density_i;
  delta_p   = Pressure_j - Pressure_i;

  RoeKappaStaticEnthalpy = 0.5 * ( StaticEnthalpy_i * Kappa_i + StaticEnthalpy_j * Kappa_j );

  s = RoeChi + RoeKappaStaticEnthalpy;

  D = s*s * delta_rho * delta_rho + delta_p * delta_p;

  delta_rhoStaticEnergy = Density_j * StaticEnergy_j - Density_i * StaticEnergy_i;

  err_P = delta_p - RoeChi * delta_rho - RoeKappa * delta_rhoStaticEnergy;

  if (abs((D - delta_p*err_P)/Density_i) > 1e-3 && abs(delta_rho/Density_i) > 1e-3 && s/Density_i > 1e-3) {

    RoeKappa = ( D * RoeKappa ) / ( D - delta_p * err_P );
    RoeChi   = ( D * RoeChi+ s*s * delta_rho * err_P ) / ( D - delta_p * err_P );

  }
}
