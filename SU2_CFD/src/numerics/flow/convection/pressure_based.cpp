/*!
 * \file centered.cpp
 * \brief Implementations of centered schemes.
 * \author F. Palacios, T. Economon
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/flow/convection/pressure_based.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CCentLinearPB_Flow::CCentLinearPB_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  dynamic_grid = config->GetDynamic_Grid();
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  MeanMassFlux = new su2double [nDim];
  Flux = new su2double [nDim];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CNumerics::ResidualType<> CCentLinearPB_Flow::ComputeResidual(const CConfig *config) {

  su2double MeanDensity, MassFlux, ProjGridMassFlux = 0.0;

  /*--- Primitive variables at point i and j ---*/

  DensityInc_i = V_i[nDim+2];  DensityInc_j = V_j[nDim+2];
  MeanDensity  = 0.5*(DensityInc_i + DensityInc_j);

  /*--- Get mass flux ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim]   = V_i[iDim+1];
    Velocity_j[iDim]   = V_j[iDim+1];

    MeanMassFlux[iDim] = 0.5 * (DensityInc_i * Velocity_i[iDim] + DensityInc_j * Velocity_j[iDim]);
    MassFlux += MeanMassFlux[iDim] * Normal[iDim];
  }

  if (dynamic_grid) {
    ProjGridMassFlux = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridMassFlux += 0.5*(DensityInc_i*GridVel_i[iDim]+DensityInc_j*GridVel_j[iDim])*Normal[iDim];
    MassFlux -= ProjGridMassFlux;
  }

  /*--- Pure central flux of momentum: (rho*u)_avg * u_avg = (rho*u)_avg * (rho*u)_avg / rho_avg  ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] = MassFlux * MeanMassFlux[iVar] / MeanDensity; 

  /*--- Compute symmetric jacobian for implicit formulations if required ---*/
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        Jacobian_i[iVar][jVar] = Jacobian_j[iVar][jVar] = 0.0;

    for (iVar = 0; iVar < nDim; iVar++) {
      for (jVar = 0; jVar < nDim; jVar++) {

        su2double delta = (iVar == jVar) ? 1.0 : 0.0;

        su2double dF_dmBar = 1/MeanDensity * (Normal[jVar] * MeanMassFlux[iVar] + MassFlux * delta);

        Jacobian_i[iVar][jVar] = 0.5 * dF_dmBar;
        Jacobian_j[iVar][jVar] = 0.5 * dF_dmBar;
      }
    }
  }

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}

CCentLinearPB_Flow::~CCentLinearPB_Flow(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanMassFlux;
  delete [] Flux;
    
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Jacobian_i;
  delete [] Jacobian_j;
  
}

// TODO: This class was directly copied from the old version by Akshay and it has not been thorougly checked
CUpwPB_Flow::CUpwPB_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  gravity = config->GetGravityForce();
  Froude = config->GetFroude();
  dynamic_grid = config->GetDynamic_Grid();
  
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  Velocity_upw = new su2double [nDim];
  MeanVelocity = new su2double [nDim];
  Flux = new su2double [nDim];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  Jacobian_upw = new su2double* [nVar];
  
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
    Jacobian_upw[iVar] = new su2double [nVar];
  }
  
}

CUpwPB_Flow::~CUpwPB_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] Velocity_upw;
  delete [] MeanVelocity;
  delete [] Flux;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
    
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
    delete [] Jacobian_upw[iVar];
  }
  delete [] Jacobian_i;
  delete [] Jacobian_j;
  delete [] Jacobian_upw;
  
}

CNumerics::ResidualType<> CUpwPB_Flow::ComputeResidual(const CConfig *config) {
	
	
  su2double MeanDensity, Flux0, Flux1, MeanPressure, Area, FF, Vel0, Vel1, ProjGridVelFlux = 0.0;
   
   
  /*--- Primitive variables at point i and j ---*/
  Pressure_i =    V_i[0];       Pressure_j = V_j[0];
  DensityInc_i =  V_i[nDim+2];  DensityInc_j = V_j[nDim+2];
  MeanDensity = 0.5*(DensityInc_i + DensityInc_j);
  MeanPressure = 0.5*(Pressure_i + Pressure_j);
  
  Area = 0.0;
  for(iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*--- (rho*u_i) ---*/
  Face_Flux = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim] + Velocity_j[iDim]);
    Face_Flux += MeanDensity*MeanVelocity[iDim]*Normal[iDim];
  }

  if (dynamic_grid) {
	  ProjGridVelFlux = 0.0;
	  for (iDim = 0; iDim < nDim; iDim++) {
		  ProjGridVelFlux += 0.5*MeanDensity*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
	  }
	  Face_Flux -= ProjGridVelFlux;
  }

  /*--- Find upwind direction. ---*/
  Flux0 = 0.5*(Face_Flux + fabs(Face_Flux)) ;
  Flux1 = 0.5*(Face_Flux - fabs(Face_Flux)) ;
  
  Upw_i = round(fabs(Flux0/(fabs(Face_Flux)+EPS)));
  Upw_j = round(fabs(Flux1/(fabs(Face_Flux)+EPS)));

  /*--- Find flux. ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Flux[iVar] = Flux0*V_i[iVar+1] + Flux1*V_j[iVar+1];
    Velocity_upw[iVar] = Upw_i*V_i[iVar+1] + Upw_j*V_j[iVar+1]; 
    if (dynamic_grid) Velocity_upw[iVar] -= (Upw_i*GridVel_i[iVar] + Upw_j*GridVel_j[iVar]); 
  }

  if (implicit) {
	  for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++) {
          Jacobian_j[iVar][jVar] = 0.0;
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_upw[iVar][jVar] = 0.0;
	    }

  unsigned short iDim;
  su2double proj_vel;

  proj_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    proj_vel += Velocity_upw[iDim]*Normal[iDim];

  if (nDim == 2) {
    Jacobian_upw[0][0] = 0.5*MeanDensity*(Velocity_upw[0]*Normal[0] + proj_vel);
    Jacobian_upw[0][1] = 0.5*MeanDensity*Velocity_upw[0]*Normal[1];

    Jacobian_upw[1][0] = 0.5*MeanDensity*Velocity_upw[1]*Normal[0];
    Jacobian_upw[1][1] = 0.5*MeanDensity*(Velocity_upw[1]*Normal[1] + proj_vel);
  }
  else {
    Jacobian_upw[0][0] = 0.5*MeanDensity*(proj_vel+Velocity_upw[0]*Normal[0]);
    Jacobian_upw[0][1] = 0.5*MeanDensity*(Velocity_upw[0]*Normal[1]);
    Jacobian_upw[0][2] = 0.5*MeanDensity*(Velocity_upw[0]*Normal[2]);

    Jacobian_upw[1][0] = 0.5*MeanDensity*(Velocity_upw[1]*Normal[0]);
    Jacobian_upw[1][1] = 0.5*MeanDensity*(proj_vel+Velocity_upw[1]*Normal[1]);
    Jacobian_upw[1][2] = 0.5*MeanDensity*(Velocity_upw[1]*Normal[2]);

    Jacobian_upw[2][0] = 0.5*MeanDensity*(Velocity_upw[2]*Normal[0]);
    Jacobian_upw[2][1] = 0.5*MeanDensity*(Velocity_upw[2]*Normal[1]);
    Jacobian_upw[2][2] = 0.5*MeanDensity*(proj_vel+Velocity_upw[2]*Normal[2]);

  }

	 for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] = Upw_i*Jacobian_upw[iVar][jVar];
        Jacobian_j[iVar][jVar] = Upw_j*Jacobian_upw[iVar][jVar];
	}
  } 
  
  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}