/*!
 * \file CUpwGeneralHLLC_Flow.cpp
 * \brief Implementation of numerics class CUpwGeneralHLLC_Flow.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/flow/convection_upwind/CUpwGeneralHLLC_Flow.hpp"

CUpwGeneralHLLC_Flow::CUpwGeneralHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
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

}

CUpwGeneralHLLC_Flow::~CUpwGeneralHLLC_Flow(void) {
  
  delete [] IntermediateState;
  delete [] dSm_dU;
  delete [] dPI_dU;
  delete [] drhoStar_dU;
  delete [] dpStar_dU;
  delete [] dEStar_dU;

  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;

}

void CUpwGeneralHLLC_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  /*--- Face area (norm or the normal vector) ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim] * Normal[iDim];

  Area = sqrt(Area);
  
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

    val_residual[0] = Density_i * ProjVelocity_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = Density_i * Velocity_i[iDim] * ProjVelocity_i + Pressure_i * UnitNormal[iDim];
    val_residual[nVar-1] = Enthalpy_i * Density_i * ProjVelocity_i;
  }
  else {

    /*--- Compute Flux Left Star from Left Star State ---*/

                rhoSL = ( sL - ProjVelocity_i ) / ( sL - sM );

    IntermediateState[0] = rhoSL * Density_i;
    for (iDim = 0; iDim < nDim; iDim++)
      IntermediateState[iDim+1] = rhoSL * ( Density_i * Velocity_i[iDim] + ( pStar - Pressure_i ) / ( sL - ProjVelocity_i ) * UnitNormal[iDim] ) ;
    IntermediateState[nVar-1] = rhoSL * ( Density_i * Energy_i - ( Pressure_i * ProjVelocity_i - pStar * sM) / ( sL - ProjVelocity_i ) );


    val_residual[0] = sM * IntermediateState[0];
    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = sM * IntermediateState[iDim+1] + pStar * UnitNormal[iDim];
    val_residual[nVar-1] = sM * ( IntermediateState[nVar-1] + pStar )  + pStar * ProjInterfaceVel;
  }
  }
  else {

  if (sR < 0.0) {

    /*--- Compute Right Flux ---*/

    val_residual[0] = Density_j * ProjVelocity_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = Density_j * Velocity_j[iDim] * ProjVelocity_j + Pressure_j * UnitNormal[iDim];
    val_residual[nVar-1] = Enthalpy_j * Density_j * ProjVelocity_j;
  }
  else {

    /*--- Compute Flux Right Star from Right Star State ---*/

                rhoSR = ( sR - ProjVelocity_j ) / ( sR - sM );

    IntermediateState[0] = rhoSR * Density_j;
    for (iDim = 0; iDim < nDim; iDim++)
      IntermediateState[iDim+1] = rhoSR * ( Density_j * Velocity_j[iDim] + ( pStar - Pressure_j ) / ( sR - ProjVelocity_j ) * UnitNormal[iDim] ) ;
    IntermediateState[nVar-1] = rhoSR * ( Density_j * Energy_j - ( Pressure_j * ProjVelocity_j - pStar * sM ) / ( sR - ProjVelocity_j ) );


    val_residual[0] = sM * IntermediateState[0];
    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = sM * IntermediateState[iDim+1] + pStar * UnitNormal[iDim];
    val_residual[nVar-1] = sM * (IntermediateState[nVar-1] + pStar )  + pStar * ProjInterfaceVel;
  }
  }

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] *= Area;


  if (implicit) {

  if (sM > 0.0) {

    if (sL > 0.0) {

      /*--- Compute Jacobian based on Left State ---*/
  
      for (iVar = 0; iVar < nVar; iVar++) 
        for (jVar = 0; jVar < nVar; jVar++) 
          val_Jacobian_j[iVar][jVar] = 0;


      GetInviscidProjJac(Velocity_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, 1.0, val_Jacobian_i);

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
        val_Jacobian_i[0][iVar] = sM * drhoStar_dU[iVar] + IntermediateState[0] * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (jDim = 0; jDim < nDim; jDim++) {

        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_i[jDim+1][iVar] = ( OmegaSM + 1 ) * ( UnitNormal[jDim] * dpStar_dU[iVar] + IntermediateState[jDim+1] * dSm_dU[iVar] );

        val_Jacobian_i[jDim+1][0] += OmegaSM * Velocity_i[jDim] * ProjVelocity_i;

        val_Jacobian_i[jDim+1][jDim+1] += OmegaSM * (sL - ProjVelocity_i);
        
        for (iDim = 0; iDim < nDim; iDim++)
          val_Jacobian_i[jDim+1][iDim+1] -= OmegaSM * Velocity_i[jDim] * UnitNormal[iDim];

        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_i[jDim+1][iVar] -= OmegaSM * dPI_dU[iVar] * UnitNormal[jDim];
      }

      /*--- Jacobian Last Row ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_i[nVar-1][iVar] = sM * ( dEStar_dU[iVar] + dpStar_dU[iVar] ) + ( EStar + pStar ) * dSm_dU[iVar];




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
        val_Jacobian_j[0][iVar] = IntermediateState[0] * ( OmegaSM + 1 ) * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_j[iDim+1][iVar] = ( OmegaSM + 1 ) * ( IntermediateState[iDim+1] * dSm_dU[iVar] + UnitNormal[iDim] * dpStar_dU[iVar] );
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_j[nVar-1][iVar] = sM * (dEStar_dU[iVar] + dpStar_dU[iVar]) + (EStar + pStar) * dSm_dU[iVar];
    }
  }
  else {
    if (sR < 0.0) {

      /*--- Compute Jacobian based on Right State ---*/
  
      for (iVar = 0; iVar < nVar; iVar++) 
        for (jVar = 0; jVar < nVar; jVar++) 
          val_Jacobian_i[iVar][jVar] = 0;

      GetInviscidProjJac(Velocity_j, &Enthalpy_j, &Chi_j, &Kappa_j, UnitNormal, 1.0, val_Jacobian_j);
    
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
        val_Jacobian_i[0][iVar] = IntermediateState[0] * ( OmegaSM + 1 ) * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_i[iDim+1][iVar] = (OmegaSM + 1) * ( IntermediateState[iDim+1] * dSm_dU[iVar] + UnitNormal[iDim] * dpStar_dU[iVar] );
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_i[nVar-1][iVar] = sM * (dEStar_dU[iVar] + dpStar_dU[iVar]) + (EStar + pStar) * dSm_dU[iVar];



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
        val_Jacobian_j[0][iVar] = sM * drhoStar_dU[iVar] + IntermediateState[0] * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (jDim = 0; jDim < nDim; jDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_j[jDim+1][iVar] = ( OmegaSM + 1 ) * ( UnitNormal[jDim] * dpStar_dU[iVar] + IntermediateState[jDim+1] * dSm_dU[iVar] );

        val_Jacobian_j[jDim+1][0] += OmegaSM * Velocity_j[jDim] * ProjVelocity_j;

        val_Jacobian_j[jDim+1][jDim+1] += OmegaSM * (sR - ProjVelocity_j);
        
        for (iDim = 0; iDim < nDim; iDim++)
          val_Jacobian_j[jDim+1][iDim+1] -= OmegaSM * Velocity_j[jDim] * UnitNormal[iDim];

        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_j[jDim+1][iVar] -= OmegaSM * dPI_dU[iVar] * UnitNormal[jDim];
      }
      
      /*--- Jacobian Last Row ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_j[nVar-1][iVar] = sM * ( dEStar_dU[iVar] + dpStar_dU[iVar] ) + ( EStar + pStar ) * dSm_dU[iVar];  
    }
  }


  /*--- Jacobians of the inviscid flux, scale = kappa because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/

  Area *= kappa;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_Jacobian_i[iVar][jVar] *= Area;
      val_Jacobian_j[iVar][jVar] *= Area;
    }
  }

  }

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

#ifdef CHECK

int UgpWithCvCompFlow::calcEulerFluxMatrices_HLLC(su2double (*val_Jacobian_i)[5], su2double (*val_Jacobian_j)[5], su2double (*val_Jacobian_i_Scal)[6], su2double (*val_Jacobian_j_Scal)[6],
                                                  const su2double Density_i, const su2double *uL, const su2double pL, const su2double TL, const su2double h0, const su2double RL, const su2double gammaL, const su2double *scalL, const su2double kL,
                                                  const su2double Density_j, const su2double *uR, const su2double pR, const su2double TR, const su2double h1, const su2double RR, const su2double gammaR, const su2double *scalR, const su2double kR,
                                                  const su2double area, const su2double *nVec, const int nScal, const su2double surfVeloc)
{

  su2double unL  = vecDotVec3d(uL, nVec);
  su2double uLuL = vecDotVec3d(uL, uL);
  su2double cL   = sqrt(gammaL*pL/Density_i);
  su2double hL   = gammaL/(gammaL-1.0)*pL/Density_i + 0.5*uLuL + kL;
  //  su2double hL   = h0 + 0.5*uLuL + kL;
  su2double eL   = hL*Density_i-pL;
  
  su2double unR  = vecDotVec3d(uR, nVec);
  su2double uRuR = vecDotVec3d(uR, uR);
  su2double cR   = sqrt(gammaR*pR/Density_j);
  su2double hR   = gammaR/(gammaR-1.0)*pR/Density_j + 0.5*uRuR + kR;
  //  su2double hR   = h1 + 0.5*uRuR + kR;
  su2double eR   = hR*Density_j-pR;
  
  
  // Roe's aveaging
  su2double Rrho = sqrt(Density_j/Density_i);
  su2double tmp = 1.0/(1.0+Rrho);
  su2double velRoe[3];
  for (int i=0; i<3; i++)
    velRoe[i] = tmp*(uL[i] + uR[i]*Rrho);
  su2double uRoe  = vecDotVec3d(velRoe, nVec);
  su2double hRoe = tmp*(hL + hR*Rrho);
  
  //  su2double cRoe  = sqrt((gammaL-1.0)*(hRoe- 0.5*vecDotVec3d(velRoe, velRoe)));
  su2double gamPdivRho = tmp*((gammaL*pL/Density_i+0.5*(gammaL-1.0)*uLuL) + (gammaR*pR/Density_j+0.5*(gammaR-1.0)*uRuR)*Rrho);
  su2double cRoe  = sqrt(gamPdivRho - ((gammaL+gammaR)*0.5-1.0)*0.5*vecDotVec3d(velRoe, velRoe));
  
  // speed of sound at L and R
  su2double sL = min(uRoe-cRoe, unL-cL);
  su2double sR = max(uRoe+cRoe, unR+cR);
  
  // speed of contact surface
  su2double sM = (pL-pR-Density_i*unL*(sL-unL)+Density_j*unR*(sR-unR))/(Density_j*(sR-unR)-Density_i*(sL-unL));
  
  // pressure at right and left (pR=pL) side of contact surface
  su2double pStar = Density_j*(unR-sR)*(unR-sM)+pR;
  
  if (sM >= 0.0) {
    
    if (sL > 0.0) {
      
      su2double nVecArea[3];
      for (int i=0; i<3; i++) nVecArea[i] = nVec[i]*area;
      
      calcJacobianA(val_Jacobian_i, uL, pL, Density_i, nVecArea, 0.5*(gammaL+gammaR), 0.0);
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_j[i][j] = 0.0;
      
    }
    else {
      
      su2double invSLmSs = 1.0/(sL-sM);
      su2double sLmuL = sL-unL;
      su2double rhoSL = Density_i*sLmuL*invSLmSs;
      su2double rhouSL[3];
      
      for (int i=0; i<3; i++)
        rhouSL[i] = (Density_i*uL[i]*sLmuL+(pStar-pL)*nVec[i])*invSLmSs;
      
      su2double eSL = (sLmuL*eL-pL*unL+pStar*sM)*invSLmSs;
      su2double gammaLM1 = (gammaL-1.0);
      su2double gammaRM1 = (gammaR-1.0);
      su2double invrhotld = 1.0/(Density_j*(sR-unR)-Density_i*(sL-unL));
      
      su2double dSMdUL[5], dSMdUR[5];
      su2double dpsdUL[5], dpsdUR[5];
      
      dSMdUL[0] = -unL*unL + uLuL*gammaLM1/2.0 + sM*sL;
      dSMdUL[1] =  nVec[0]*(2.0*unL-sL-sM) - gammaLM1*uL[0];
      dSMdUL[2] =  nVec[1]*(2.0*unL-sL-sM) - gammaLM1*uL[1];
      dSMdUL[3] =  nVec[2]*(2.0*unL-sL-sM) - gammaLM1*uL[2];
      dSMdUL[4] =  gammaLM1;
      
      for (iVar = 0; iVar < nVar; iVar++)
      {
        dSMdUL[i] *= invrhotld;
        dpsdUL[i] = Density_j*(sR-unR)*dSMdUL[i];
      }
      
      dSMdUR[0] =  unR*unR - uRuR*gammaRM1/2.0 - sM*sR;
      dSMdUR[1] = -nVec[0]*(2.0*unR-sR-sM) + gammaRM1*uR[0];
      dSMdUR[2] = -nVec[1]*(2.0*unR-sR-sM) + gammaRM1*uR[1];
      dSMdUR[3] = -nVec[2]*(2.0*unR-sR-sM) + gammaRM1*uR[2];
      dSMdUR[4] = -gammaRM1;
      
      for (iVar = 0; iVar < nVar; iVar++)
      {
        dSMdUR[i] *= invrhotld;
        dpsdUR[i] = Density_i*(sL-unL)*dSMdUR[i];
      }
      
      calcSubSonicJacobeanHLLC(val_Jacobian_i, val_Jacobian_j,
                               Density_i, uL, pL, eL, unL, uLuL, sL,
                               rhoSL, rhouSL, eSL, dSMdUL,
                               dSMdUR, dpsdUL, dpsdUR, sM, pStar, 0.5*(gammaL+gammaR), nVec);
      
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[0][i] =  val_Jacobian_i[0][i]*sM + dSMdUL[i]*rhoSL;
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[1][i] =  val_Jacobian_i[1][i]*sM + dSMdUL[i]*rhouSL[0] + dpsdUL[i]*nVec[0];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[2][i] =  val_Jacobian_i[2][i]*sM + dSMdUL[i]*rhouSL[1] + dpsdUL[i]*nVec[1];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[3][i] =  val_Jacobian_i[3][i]*sM + dSMdUL[i]*rhouSL[2] + dpsdUL[i]*nVec[2];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[4][i] = (val_Jacobian_i[4][i]+dpsdUL[i])*sM + (eSL+pStar)*dSMdUL[i];
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_i[i][j] *= area;
      
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[0][i] =  val_Jacobian_j[0][i]*sM + dSMdUR[i]*rhoSL;
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[1][i] =  val_Jacobian_j[1][i]*sM + dSMdUR[i]*rhouSL[0] + dpsdUR[i]*nVec[0];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[2][i] =  val_Jacobian_j[2][i]*sM + dSMdUR[i]*rhouSL[1] + dpsdUR[i]*nVec[1];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[3][i] =  val_Jacobian_j[3][i]*sM + dSMdUR[i]*rhouSL[2] + dpsdUR[i]*nVec[2];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[4][i] = (val_Jacobian_j[4][i]+dpsdUR[i])*sM + (eSL+pStar)*dSMdUR[i];
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_j[i][j] *= area;
      
    }
  }
  
  else {
    
    if (sR >= 0.0) {
      
      su2double invSRmSs = 1.0/(sR-sM);
      su2double sRmuR = sR-unR;
      su2double rhoSR = Density_j*sRmuR*invSRmSs;
      su2double rhouSR[3];
      for (int i=0; i<3; i++)
        rhouSR[i] = (Density_j*uR[i]*sRmuR+(pStar-pR)*nVec[i])*invSRmSs;
      su2double eSR = (sRmuR*eR-pR*unR+pStar*sM)*invSRmSs;
      su2double gammaLM1 = (gammaL-1.0);
      su2double gammaRM1 = (gammaR-1.0);
      su2double invrhotld = 1.0/(Density_j*(sR-unR)-Density_i*(sL-unL));
      
      su2double dSMdUL[5], dSMdUR[5];
      su2double dpsdUL[5], dpsdUR[5];
      
      dSMdUL[0] = -unL*unL + uLuL*gammaLM1/2.0 + sM*sL;
      dSMdUL[1] =  nVec[0]*(2.0*unL-sL-sM) - gammaLM1*uL[0];
      dSMdUL[2] =  nVec[1]*(2.0*unL-sL-sM) - gammaLM1*uL[1];
      dSMdUL[3] =  nVec[2]*(2.0*unL-sL-sM) - gammaLM1*uL[2];
      dSMdUL[4] =  gammaLM1;
      
      for (iVar = 0; iVar < nVar; iVar++) {
        dSMdUL[i] *= invrhotld;
        dpsdUL[i] = Density_j*(sR-unR)*dSMdUL[i];
      }
      
      dSMdUR[0] =  unR*unR - uRuR*gammaRM1/2.0 - sM*sR;
      dSMdUR[1] = -nVec[0]*(2.0*unR-sR-sM) + gammaRM1*uR[0];
      dSMdUR[2] = -nVec[1]*(2.0*unR-sR-sM) + gammaRM1*uR[1];
      dSMdUR[3] = -nVec[2]*(2.0*unR-sR-sM) + gammaRM1*uR[2];
      dSMdUR[4] = -gammaRM1;
      
      for (iVar = 0; iVar < nVar; iVar++) {
        dSMdUR[i] *= invrhotld;
        dpsdUR[i] = Density_i*(sL-unL)*dSMdUR[i];
      }
      
      calcSubSonicJacobeanHLLC(val_Jacobian_j, val_Jacobian_i,
                               Density_j, uR, pR, eR, unR, uRuR, sR,
                               rhoSR, rhouSR, eSR,
                               dSMdUR, dSMdUL, dpsdUR, dpsdUL, sM, pStar, 0.5*(gammaL+gammaR), nVec);
      
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[0][i] =  val_Jacobian_i[0][i]*sM + dSMdUL[i]*rhoSR;
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[1][i] =  val_Jacobian_i[1][i]*sM + dSMdUL[i]*rhouSR[0] + dpsdUL[i]*nVec[0];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[2][i] =  val_Jacobian_i[2][i]*sM + dSMdUL[i]*rhouSR[1] + dpsdUL[i]*nVec[1];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[3][i] =  val_Jacobian_i[3][i]*sM + dSMdUL[i]*rhouSR[2] + dpsdUL[i]*nVec[2];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[4][i] = (val_Jacobian_i[4][i]+dpsdUL[i])*sM + (eSR+pStar)*dSMdUL[i];
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_i[i][j] *= area;
      
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[0][i] =  val_Jacobian_j[0][i]*sM + dSMdUR[i]*rhoSR;
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[1][i] =  val_Jacobian_j[1][i]*sM + dSMdUR[i]*rhouSR[0] + dpsdUR[i]*nVec[0];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[2][i] =  val_Jacobian_j[2][i]*sM + dSMdUR[i]*rhouSR[1] + dpsdUR[i]*nVec[1];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[3][i] =  val_Jacobian_j[3][i]*sM + dSMdUR[i]*rhouSR[2] + dpsdUR[i]*nVec[2];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[4][i] = (val_Jacobian_j[4][i]+dpsdUR[i])*sM + (eSR+pStar)*dSMdUR[i];
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_j[i][j] *= area;
      
    }
    
    else {
      
      su2double nVecArea[3];
      for (int i=0; i<3; i++)        nVecArea[i] = nVec[i]*area;
      calcJacobianA(val_Jacobian_j, uR, pR, Density_j, nVecArea, 0.5*(gammaL+gammaR), 0.0);
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_i[i][j] = 0.0;
      
    }
    
  }
  
}

void UgpWithCvCompFlow::calcSubSonicJacobeanHLLC(su2double (*AL)[5], su2double (*AR)[5],
                                                 su2double Density_i, const su2double *uL, su2double pL, su2double eL, su2double qL, su2double psiL, su2double SL,
                                                 su2double rhoSL, su2double *rhouSL, su2double eSL,
                                                 su2double *dSMdUL, su2double *dSMdUR, su2double *dpsdUL, su2double *dpsdUR, su2double SM, su2double pS,
                                                 su2double gamma, const su2double *nV) // nV is not normalized
{
  
  su2double gammaMinus1 = (gamma-1.0);
  su2double omL = 1.0/(SL-SM);
  
  AL[0][0] =  SL    + rhoSL*dSMdUL[0];
  AL[0][1] = -nV[0] + rhoSL*dSMdUL[1];
  AL[0][2] = -nV[1] + rhoSL*dSMdUL[2];
  AL[0][3] = -nV[2] + rhoSL*dSMdUL[3];
  AL[0][4] =        + rhoSL*dSMdUL[4];
  
  AL[1][0] =    qL*uL[0]       - nV[0]*psiL*gammaMinus1/2.0   + nV[0]*dpsdUL[0] + rhouSL[0]*dSMdUL[0];
  AL[1][1] =  SL - qL          + nV[0]*(gamma-2.0)*uL[0]      + nV[0]*dpsdUL[1] + rhouSL[0]*dSMdUL[1];
  AL[1][2] =     - uL[0]*nV[1] + nV[0]*gammaMinus1*uL[1]      + nV[0]*dpsdUL[2] + rhouSL[0]*dSMdUL[2];
  AL[1][3] =     - uL[0]*nV[2] + nV[0]*gammaMinus1*uL[2]      + nV[0]*dpsdUL[3] + rhouSL[0]*dSMdUL[3];
  AL[1][4] = -gammaMinus1*nV[0]                               + nV[0]*dpsdUL[4] + rhouSL[0]*dSMdUL[4];
  
  AL[2][0] =    qL*uL[1]       - nV[1]*psiL*gammaMinus1/2.0   + nV[1]*dpsdUL[0] + rhouSL[1]*dSMdUL[0];
  AL[2][1] =     - uL[1]*nV[0] + nV[1]*gammaMinus1*uL[0]      + nV[1]*dpsdUL[1] + rhouSL[1]*dSMdUL[1];
  AL[2][2] =  SL - qL          + nV[1]*(gamma-2.0)*uL[1]      + nV[1]*dpsdUL[2] + rhouSL[1]*dSMdUL[2];
  AL[2][3] =     - uL[1]*nV[2] + nV[1]*gammaMinus1*uL[2]      + nV[1]*dpsdUL[3] + rhouSL[1]*dSMdUL[3];
  AL[2][4] = -gammaMinus1*nV[1]                               + nV[1]*dpsdUL[4] + rhouSL[1]*dSMdUL[4];
  
  AL[3][0] =    qL*uL[2]       - nV[2]*psiL*gammaMinus1/2.0   + nV[2]*dpsdUL[0] + rhouSL[2]*dSMdUL[0];
  AL[3][1] =     - uL[2]*nV[0] + nV[2]*gammaMinus1*uL[0]      + nV[2]*dpsdUL[1] + rhouSL[2]*dSMdUL[1];
  AL[3][2] =     - uL[2]*nV[1] + nV[2]*gammaMinus1*uL[1]      + nV[2]*dpsdUL[2] + rhouSL[2]*dSMdUL[2];
  AL[3][3] =  SL - qL          + nV[2]*(gamma-2.0)*uL[2]      + nV[2]*dpsdUL[3] + rhouSL[2]*dSMdUL[3];
  AL[3][4] = -gammaMinus1*nV[2]                               + nV[2]*dpsdUL[4] + rhouSL[2]*dSMdUL[4];
  
  AL[4][0] =      qL*(eL+pL)/Density_i - qL*psiL*(gamma-1.0)/2.0   + SM*dpsdUL[0] + (pS+eSL)*dSMdUL[0];
  AL[4][1] = - nV[0]*(eL+pL)/Density_i + gammaMinus1*uL[0]*qL      + SM*dpsdUL[1] + (pS+eSL)*dSMdUL[1];
  AL[4][2] = - nV[1]*(eL+pL)/Density_i + gammaMinus1*uL[1]*qL      + SM*dpsdUL[2] + (pS+eSL)*dSMdUL[2];
  AL[4][3] = - nV[2]*(eL+pL)/Density_i + gammaMinus1*uL[2]*qL      + SM*dpsdUL[3] + (pS+eSL)*dSMdUL[3];
  AL[4][4] =   SL-qL*gamma                                    + SM*dpsdUL[4] + (pS+eSL)*dSMdUL[4];
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      AL[i][j] *= omL;
  
  for (iVar = 0; iVar < nVar; iVar++)    AR[0][i] = omL*rhoSL*dSMdUR[i];
  for (iVar = 0; iVar < nVar; iVar++)    AR[1][i] = omL*(nV[0]*dpsdUR[i]+rhouSL[0]*dSMdUR[i]);
  for (iVar = 0; iVar < nVar; iVar++)    AR[2][i] = omL*(nV[1]*dpsdUR[i]+rhouSL[1]*dSMdUR[i]);
  for (iVar = 0; iVar < nVar; iVar++)    AR[3][i] = omL*(nV[2]*dpsdUR[i]+rhouSL[2]*dSMdUR[i]);
  for (iVar = 0; iVar < nVar; iVar++)    AR[4][i] = omL*(dpsdUR[i]*SM+(pS+eSL)*dSMdUR[i]);
  
}

void UgpWithCvCompFlow::calcJacobianA(su2double (*A)[5], const su2double *vel, su2double pp, su2double rrho, const su2double *nV, su2double gamma, su2double surfVeloc) // nV is not normalized
{
 
  su2double kapm1 = (gamma - 1.0);
  
  su2double nVel[3];
  nVel[0] = vel[0]*nV[0];
  nVel[1] = vel[1]*nV[1];
  nVel[2] = vel[2]*nV[2];
  su2double U_k = nVel[0]+nVel[1]+nVel[2];
  su2double vSquHlf = 0.5*vecDotVec3d(vel, vel);
  su2double c = sqrt(gamma*pp/rrho);
  su2double inv_kap_m1 = 1.0/kapm1;
  
  A[0][0] =-surfVeloc;
  A[0][1] = nV[0];
  A[0][2] = nV[1];
  A[0][3] = nV[2];
  A[0][4] = 0.0;
  
  A[1][0] = -vel[0]*(nVel[1]+nVel[2])+nV[0]*(kapm1*vSquHlf-vel[0]*vel[0]);
  A[1][1] = (2.-gamma)*nVel[0]+U_k-surfVeloc;
  A[1][2] = vel[0]*nV[1]-kapm1*vel[1]*nV[0];
  A[1][3] = vel[0]*nV[2]-kapm1*vel[2]*nV[0];
  A[1][4] = kapm1*nV[0];
  
  A[2][0] = -vel[1]*(nVel[0]+nVel[2])+nV[1]*(kapm1*vSquHlf-vel[1]*vel[1]);
  A[2][1] = -kapm1*vel[0]*nV[1]+ vel[1]*nV[0];
  A[2][2] = (2.-gamma)*nVel[1]+U_k-surfVeloc;
  A[2][3] = vel[1]*nV[2]-kapm1*vel[2]*nV[1];
  A[2][4] = kapm1*nV[1];
  
  A[3][0] = -vel[2]*(nVel[0]+nVel[1])+nV[2]*(kapm1*vSquHlf-vel[2]*vel[2]);
  A[3][1] = -kapm1*vel[0]*nV[2]+vel[2]*nV[0];
  A[3][2] = -kapm1*vel[1]*nV[2]+vel[2]*nV[1];
  A[3][3] = (2.-gamma)*nVel[2]+U_k-surfVeloc;
  A[3][4] = kapm1*nV[2];
  
  A[4][0] = U_k*((gamma-2.)*vSquHlf-c*c*inv_kap_m1);
  A[4][1] = c*c*inv_kap_m1*nV[0]-kapm1*vel[0]*(nVel[1]+nVel[2])-(kapm1*vel[0]*vel[0]-vSquHlf)*nV[0];
  A[4][2] = c*c*inv_kap_m1*nV[1]-kapm1*vel[1]*(nVel[0]+nVel[2])-(kapm1*vel[1]*vel[1]-vSquHlf)*nV[1];
  A[4][3] = c*c*inv_kap_m1*nV[2]-kapm1*vel[2]*(nVel[0]+nVel[1])-(kapm1*vel[2]*vel[2]-vSquHlf)*nV[2];
  A[4][4] = gamma*U_k-surfVeloc;
  
}


#endif

