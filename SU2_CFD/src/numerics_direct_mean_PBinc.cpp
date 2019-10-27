/*!
 * \file numerics_direct_mean_inc.cpp
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

#include "../include/numerics_structure.hpp"
#include <limits>

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
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  val_Jacobian_upw = new su2double* [nVar];
  
  
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
    val_Jacobian_upw[iVar] = new su2double [nVar];
  }
  
}

CUpwPB_Flow::~CUpwPB_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] Velocity_upw;
  delete [] MeanVelocity;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
    delete [] val_Jacobian_upw[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  delete [] val_Jacobian_upw;
  
}

void CUpwPB_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
	
	
  su2double MeanDensity, Flux0, Flux1, MeanPressure, Area, FF, Vel0, Vel1, ProjGridVelFlux = 0.0;
   
   
  /*--- Primitive variables at point i and j ---*/
  Pressure_i =    V_i[0];       Pressure_j = V_j[0];
  DensityInc_i =  V_i[nDim+1];  DensityInc_j = V_j[nDim+1];
  MeanDensity = 0.5*(DensityInc_i + DensityInc_j);
  MeanPressure = 0.5*(Pressure_i + Pressure_j);
  
  Area = 0.0;
  for(iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
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
		  ProjGridVelFlux   += 0.5*MeanDensity*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim]; 
	  }
	  Face_Flux -= ProjGridVelFlux;
  } 
  
 
  Flux0 = 0.5*(Face_Flux + fabs(Face_Flux)) ;
  Flux1 = 0.5*(Face_Flux - fabs(Face_Flux)) ;
  
  Upw_i = round(fabs(Flux0/(fabs(Face_Flux)+EPS)));
  Upw_j = round(fabs(Flux1/(fabs(Face_Flux)+EPS)));
     
  for (iVar = 0; iVar < nVar; iVar++) {
	  val_residual[iVar] = Flux0*V_i[iVar+1] + Flux1*V_j[iVar+1];
	  Velocity_upw[iVar] = Upw_i*V_i[iVar+1] + Upw_j*V_j[iVar+1]; 
	  if (dynamic_grid) Velocity_upw[iVar] -= (Upw_i*GridVel_i[iVar] + Upw_j*GridVel_j[iVar]); 
  }
    
    	
  if (implicit) {
	  
	  for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_j[iVar][jVar] = 0.0;
        val_Jacobian_i[iVar][jVar] = 0.0;
        val_Jacobian_upw[iVar][jVar] = 0.0;
	}
	  
	GetInviscidPBProjJac(&DensityInc_i, Velocity_upw,  Normal, 0.5, val_Jacobian_upw);
	
	
	/*if (Face_Flux > 0) 
         GetInviscidPBProjJac(&DensityInc_i, Velocity_i,  Normal, 0.5, val_Jacobian_i);
    else
         GetInviscidPBProjJac(&DensityInc_j, Velocity_j,  Normal, 0.5, val_Jacobian_j);*/
	
	 for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_i[iVar][jVar] = Upw_i*val_Jacobian_upw[iVar][jVar];
        val_Jacobian_j[iVar][jVar] = Upw_j*val_Jacobian_upw[iVar][jVar];
	}
  } 
}

CCentPB_Flow::CCentPB_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  gravity = config->GetGravityForce();
  Froude = config->GetFroude();
  
  /*--- Artificial dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_Flow();
  
  /*--- Allocate some structures ---*/
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  MeanVelocity = new su2double [nDim];
  ProjFlux = new su2double [nVar];
  
}

CCentPB_Flow::~CCentPB_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
  
}

void CCentPB_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                                           CConfig *config) {
  
  su2double U_i[3] = {0.0,0.0,0.0}, U_j[3] = {0.0,0.0,0.0};

  /*--- Conservative variables at point i and j ---*/
  
  Pressure_i =    V_i[0];       Pressure_j = V_j[0];
  DensityInc_i =  V_i[nDim+1];  DensityInc_j = V_j[nDim+1];
  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }
  
  /*--- Recompute conservative variables ---*/

  //U_i[0] = Pressure_i; U_j[0] = Pressure_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim] = DensityInc_i*Velocity_i[iDim]; U_j[iDim] = DensityInc_j*Velocity_j[iDim];
  }
  
  /*--- Compute mean values of the variables ---*/
  
  MeanDensity = 0.5*(DensityInc_i+DensityInc_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  
  /*--- Get projected flux tensor ---*/
  
  GetInviscidPBProjFlux(&MeanDensity, MeanVelocity, &MeanPressure,  Normal, ProjFlux);
  
  /*--- Compute inviscid residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];
  
  /*--- Jacobians of the inviscid flux ---*/
  
  if (implicit) {
    GetInviscidPBProjJac(&MeanDensity, MeanVelocity,  Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }
  
  /*--- Computes differences btw. conservative variables ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  
  /*--- Compute the local espectral radius and the stretching factor ---*/
  
  ProjVelocity_i = 0; ProjVelocity_j = 0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  //SoundSpeed_i = sqrt(ProjVelocity_i*ProjVelocity_i + (BetaInc2_i/DensityInc_i)*Area*Area);
  //SoundSpeed_j = sqrt(ProjVelocity_j*ProjVelocity_j + (BetaInc2_j/DensityInc_j)*Area*Area);

  Local_Lambda_i = fabs(2.0*ProjVelocity_i);
  Local_Lambda_j = fabs(2.0*ProjVelocity_j);
  MeanLambda = 0.5*(Local_Lambda_i + Local_Lambda_j);
  
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc0 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc0*su2double(nDim)/3.0;
  
  /*--- Compute viscous part of the residual ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;
  
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++) {
      val_Jacobian_i[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
      val_Jacobian_j[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
    }
  }
  
}

CPressureSource::CPressureSource(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CPressureSource::~CPressureSource(void) { }

void CPressureSource::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {

    unsigned short iDim;

    Pressure_i = V_i[0]; Pressure_j = V_j[0];
    MeanPressure = 0.5*(Pressure_i + Pressure_j);

    for (iDim =0;iDim<nDim;iDim++)
       val_residual[iDim] = MeanPressure*Normal[iDim];

}



CAvgGradPBInc_Flow::CAvgGradPBInc_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  energy   = config->GetEnergy_Equation();

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
  
}

void CAvgGradPBInc_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
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
  
  /*--- Update viscous residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  }

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
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    else {
      GetViscousPBIncProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, dist_ij, UnitNormal,
                                Area, val_Jacobian_i, val_Jacobian_j);
    }
    
  }
}

CAvgGradCorrectedPBInc_Flow::CAvgGradCorrectedPBInc_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  energy   = config->GetEnergy_Equation();

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
  
}

void CAvgGradCorrectedPBInc_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

 su2double  Mean_GradVar_Edge[5], GradEdge[5];
 su2double  Mean_GradVar_Face[5][3], Mean_GradVar[5][3];
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
  
  /*--- Correct the face gradient for odd-even decoupling ---*/
  /*--- Steps are 
   * 1. Interpolate the gradient at the face -> du/ds|_in
   * 2. Find the projection of the interpolated gradient on the edge vector -> (du/ds|_in) . e
   * 3. Find the gradient as the difference between the neighboring nodes -> (u_j - u_i)/ds
   * 4. Correct the gradient at the face using du/ds = du/ds|_in + ( (u_j - u_i)/ds - (du/ds|_in . e) ) . e ---*/
  
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
	  Mean_GradVar_Edge[iVar] = 0.0;	  
	  
	  GradEdge[iVar] = (PrimVar_j[iVar] - PrimVar_i[iVar])/sqrt(dist_ij_2);
	  
	  for (iDim = 0; iDim < nDim; iDim++) {
		  Mean_GradVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);  
		  
		  Mean_GradVar_Edge[iVar] += Mean_GradVar[iVar][iDim]*Edge_Vector[iDim]/sqrt(dist_ij_2);
       }
   }
   
   for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
	   for (iDim = 0; iDim < nDim; iDim++) {
		   Mean_GradVar_Face[iVar][iDim] = Mean_GradVar[iVar][iDim] + 
                                        (GradEdge[iVar] - Mean_GradVar_Edge[iVar])*Edge_Vector[iDim]/sqrt(dist_ij_2); 
       }
   }
  
  
  /*--- Projection of the mean gradient in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar+1][iDim] + PrimVar_Grad_j[iVar+1][iDim]);
      //Mean_GradPrimVar[iVar][iDim] = Mean_GradVar_Face[iVar+1][iDim];
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
  
  /*--- Update viscous residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  }

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
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    else {
      GetViscousPBIncProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, sqrt(dist_ij_2), UnitNormal,
                                Area, val_Jacobian_i, val_Jacobian_j);
    }
    
  }

}
