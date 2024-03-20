/*!
 * \file pbflow.cpp
 * \brief This file contains the numerical methods for incompressible flow.
 * \author F. Palacios, T. Economon
 * \version 6.0.1 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../include/numerics/pbflow.hpp"

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
  DensityInc_i =  V_i[nDim+1];  DensityInc_j = V_j[nDim+1];
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

	GetInviscidPBProjJac(MeanDensity, Velocity_upw,  Normal, 0.5, Jacobian_upw);

	 for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] = Upw_i*Jacobian_upw[iVar][jVar];
        Jacobian_j[iVar][jVar] = Upw_j*Jacobian_upw[iVar][jVar];
	}
  } 
  
  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}

CCentPB_Flow::CCentPB_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
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

CCentPB_Flow::~CCentPB_Flow(void) {
  
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

CNumerics::ResidualType<> CCentPB_Flow::ComputeResidual(const CConfig *config) {
  /*--- To do list
   * 1. Find upwind direction, upwind node and downwind node
   * 2. Compute normalised variable for the upwind node using the gradient of the upwind node
   * 3. Find face velocity using central difference scheme
   * 4. Use the ULTIMATE limiter to find the adjusted face velocity
   * 5. Find residual as FaceFlux * AdjustedFaceVelocity
   * */
  
  
  su2double MeanDensity, Flux0, Flux1, MeanPressure, Area, FF, Vel0, Vel1, ProjGridVelFlux = 0.0;
  su2double dissipation, kappa=0.15, ProjVel_i = 0.0, ProjVel_j = 0.0;

  /*--- Conservative variables at point i and j ---*/
  
  Pressure_i =    V_i[0];       Pressure_j = V_j[0];
  DensityInc_i =  V_i[nDim+1];  DensityInc_j = V_j[nDim+1];

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
  }
  
  Area = 0.0;
  for(iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  Face_Flux = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim] + Velocity_j[iDim]);
    Face_Flux += MeanDensity*MeanVelocity[iDim]*Normal[iDim];
    ProjVel_i += Velocity_i[iDim]*Normal[iDim];
    ProjVel_j += Velocity_j[iDim]*Normal[iDim];
  }
  
  if (dynamic_grid) {
	  ProjGridVelFlux = 0.0; 
	  for (iDim = 0; iDim < nDim; iDim++) { 
		  ProjGridVelFlux   += 0.5*MeanDensity*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim]; 
	  }
	  Face_Flux -= ProjGridVelFlux;
  } 

  Flux0 = 0.5*(Face_Flux + fabs(Face_Flux)) ;
  Flux1 = 0.5*(Face_Flux - fabs(Face_Flux)) ;

  Upw_i = round(fabs(Flux0/(fabs(Face_Flux)+EPS)));
  Upw_j = round(fabs(Flux1/(fabs(Face_Flux)+EPS)));

  for (iVar = 0; iVar < nVar; iVar++) {
	  Velocity_upw[iVar] = Upw_i*V_i[iVar+1] + Upw_j*V_j[iVar+1]; 
	  if (dynamic_grid) Velocity_upw[iVar] -= (Upw_i*GridVel_i[iVar] + Upw_j*GridVel_j[iVar]); 
	  Flux[iVar] = Face_Flux*MeanVelocity[iVar];
  }
  
  //Dissipation
  su2double lambda_i = 0.0, lambda_j = 0.0;
  lambda_i = 2.0*abs(ProjVel_i);
  lambda_j = 2.0*abs(ProjVel_j);

  su2double lambda_mean = 0.5*(lambda_i + lambda_j);
  if (lambda_mean < EPS) {
    lambda_i = 2.0*abs(config->GetVelocity_Ref()*Area);
    lambda_j = 2.0*abs(config->GetVelocity_Ref()*Area);
    lambda_mean = abs(config->GetVelocity_Ref()*Area);
  }

  su2double Param_p = 0.3, SF=0.0;
  su2double Phi_i = pow((lambda_i/(4.0*lambda_mean)),Param_p);
  su2double Phi_j = pow((lambda_j/(4.0*lambda_mean)),Param_p);

  if ((Phi_i + Phi_j) != 0.0)
    SF = 4.0*Phi_i*Phi_j/(Phi_i + Phi_j);

  su2double sc0 = 3.0*(Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j);
  su2double E_0 = kappa*sc0*2.0/3.0;

  su2double diss[MAXNDIM];
  for (iDim = 0; iDim < nDim; iDim++) 
    diss[iDim] = E_0*(DensityInc_i*Velocity_i[iDim] - DensityInc_j*Velocity_j[iDim])*SF*lambda_mean;

  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] += diss[iVar];

  /*--- For implicit schemes, compute jacobians based on the upwind scheme for stability issues. (See ANSYS user guide) ---*/
  if (implicit) {
	  for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_j[iVar][jVar] = 0.0;
        Jacobian_i[iVar][jVar] = 0.0;
        Jacobian_upw[iVar][jVar] = 0.0;
	}

	GetInviscidPBProjJac(MeanDensity, MeanVelocity,  Normal, 0.5, Jacobian_upw);
	//GetInviscidPBProjJac(&DensityInc_i, Velocity_upw,  Normal, 0.5, val_Jacobian_upw);
	/*GetInviscidPBProjJac(&DensityInc_i, Velocity_i,  Normal, 0.5, Jacobian_i);
	GetInviscidPBProjJac(&DensityInc_i, Velocity_j,  Normal, 0.5, Jacobian_j);*/

	 for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] = Upw_i*Jacobian_upw[iVar][jVar];
        Jacobian_j[iVar][jVar] = Upw_j*Jacobian_upw[iVar][jVar];
	}
	for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar][iVar] += E_0*DensityInc_i*SF*lambda_mean;
      Jacobian_i[iVar][iVar] -= E_0*DensityInc_j*SF*lambda_mean;
    }
	
  }
  
  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

CAvgGradPBInc_Flow::CAvgGradPBInc_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  /*--- Primitive flow variables. ---*/

  PrimVar_i    = new su2double [nDim+4];
  PrimVar_j    = new su2double [nDim+4];
  Mean_PrimVar = new su2double [nDim+4];

  /*--- Incompressible flow, primitive variables nDim+2, (P, vx, vy, vz, rho) ---*/
  
  Mean_GradPrimVar = new su2double*[nVar];
  
  /*--- Incompressible flow, gradient primitive variables nDim+2, (P, vx, vy, vz, rho) ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradPrimVar[iVar] = new su2double[nDim];
  
}

CAvgGradPBInc_Flow::~CAvgGradPBInc_Flow(void) {

  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Mean_PrimVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Jacobian_i;
  delete [] Jacobian_j;
    
}

CNumerics::ResidualType<> CAvgGradPBInc_Flow::ComputeResidual(const CConfig *config) {
  
  /*--- Normalized normal vector ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  for (iVar = 0; iVar < nDim+4; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }

  /*--- Density and transport properties ---*/
  
  Laminar_Viscosity_i    = V_i[nDim+2];  Laminar_Viscosity_j    = V_j[nDim+2];
  Eddy_Viscosity_i       = V_i[nDim+3];  Eddy_Viscosity_j       = V_j[nDim+3];

  /*--- Mean transport properties ---*/
  
  Mean_Laminar_Viscosity    = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity       = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke              = 0.0;//0.5*(turb_ke_i + turb_ke_j);

  /*--- Mean gradient approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (iDim = 0; iDim < nDim; iDim++)
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar+1][iDim] + PrimVar_Grad_j[iVar+1][iDim]);
      
  /*--- Get projected flux tensor ---*/
  
  GetViscousPBIncProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, Mean_turb_ke);
  

  /*--- Implicit part ---*/
  
  if (implicit) {
    
    dist_ij = 0.0; proj_vector_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      dist_ij        += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      proj_vector_ij += (Coord_j[iDim]-Coord_i[iDim])*Normal[iDim];
    }
    proj_vector_ij = proj_vector_ij/dist_ij;
    dist_ij = sqrt(dist_ij);

    if (dist_ij == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    else {
      GetViscousPBIncProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, dist_ij, UnitNormal,
                                Area, Jacobian_i, Jacobian_j);
    }
    
  }
  
  return ResidualType<>(Proj_Flux_Tensor, Jacobian_i, Jacobian_j);
}

CAvgGradCorrectedPBInc_Flow::CAvgGradCorrectedPBInc_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  /*--- Primitive flow variables. ---*/

  Mean_PrimVar = new su2double [nDim+4];
  PrimVar_i    = new su2double [nDim+4];
  PrimVar_j    = new su2double [nDim+4];
  Proj_Mean_GradPrimVar_Edge = new su2double [nVar];
  Edge_Vector = new su2double [nDim];
  
  Mean_GradPrimVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradPrimVar[iVar] = new su2double [nDim];
  
}

CAvgGradCorrectedPBInc_Flow::~CAvgGradCorrectedPBInc_Flow(void) {

  delete [] Mean_PrimVar;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Proj_Mean_GradPrimVar_Edge;
  delete [] Edge_Vector;
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Jacobian_i;
  delete [] Jacobian_j;
  
}

CNumerics::ResidualType<> CAvgGradCorrectedPBInc_Flow::ComputeResidual(const CConfig *config) {

 unsigned short nPrimVarGrad = nDim+2, nPrimVar = nDim+4;

 /*--- Normalized normal vector ---*/

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Conversion to Primitive Variables (P, u, v, w, rho, mu, muT) ---*/
  
  for (iVar = 0; iVar < nPrimVar; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }
  
  /*--- Density and transport properties ---*/
  
  DensityInc_i           = V_i[nDim+1];  DensityInc_j           = V_j[nDim+1];
  Laminar_Viscosity_i    = V_i[nDim+2];  Laminar_Viscosity_j    = V_j[nDim+2];
  Eddy_Viscosity_i       = V_i[nDim+3];  Eddy_Viscosity_j       = V_j[nDim+3];
  
  /*--- Mean transport properties ---*/
  
  Mean_Laminar_Viscosity    = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity       = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke              = 0.0;//0.5*(turb_ke_i + turb_ke_j);
 
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }

  /*--- Projection of the mean gradient in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar+1][iDim] + PrimVar_Grad_j[iVar+1][iDim]);
      Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
    }
    if (dist_ij_2 != 0.0) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                         (PrimVar_j[iVar+1]-PrimVar_i[iVar+1]))*Edge_Vector[iDim] / dist_ij_2;
      }
    }
  }
  
  /*--- Get projected flux tensor ---*/
  
  GetViscousPBIncProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, Mean_turb_ke);

  /*--- Implicit part ---*/

  if (implicit) {

    proj_vector_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      proj_vector_ij += (Coord_j[iDim]-Coord_i[iDim])*Normal[iDim];
    }
    proj_vector_ij = proj_vector_ij/dist_ij_2;

    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    else {
      GetViscousPBIncProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, sqrt(dist_ij_2), UnitNormal,
                                Area, Jacobian_i, Jacobian_j);
    }
    
  }

 return ResidualType<>(Proj_Flux_Tensor, Jacobian_i, Jacobian_j);
}
