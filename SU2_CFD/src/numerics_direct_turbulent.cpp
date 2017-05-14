/*!
 * \file numerics_direct_turbulent.cpp
 * \brief This file contains all the convective term discretization.
 * \author F. Palacios, A. Bueno
 * \version 4.3.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

CUpwSca_TurbSA::CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                               CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  
}

CUpwSca_TurbSA::~CUpwSca_TurbSA(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  
}

void CUpwSca_TurbSA::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  q_ij = 0.0;
  
  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+1); AD::SetPreaccIn(V_j, nDim+1);
  AD::SetPreaccIn(TurbVar_i[0]); AD::SetPreaccIn(TurbVar_j[0]);
  AD::SetPreaccIn(Normal, nDim);
  if (grid_movement){
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  } else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  val_residual[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_j[0][0] = a1;
  }
  
  AD::SetPreaccOut(val_residual[0]);
  AD::EndPreacc();
}

CAvgGrad_TurbSA::CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGrad_TurbSA::~CAvgGrad_TurbSA(void) {
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGrad_TurbSA::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {    
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+nu_e*proj_vector_ij)/sigma;
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

CAvgGrad_TurbSA_Neg::CAvgGrad_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  cn1   = 16.0;
  fn    = 0.0;

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGrad_TurbSA_Neg::~CAvgGrad_TurbSA_Neg(void) {
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGrad_TurbSA_Neg::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5);   AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7);   AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  
  nu_ij = 0.5*(nu_i+nu_j);
  nu_tilde_ij = 0.5*(TurbVar_i[0]+TurbVar_j[0]);

  Xi = nu_tilde_ij/nu_ij;
  
  if (nu_tilde_ij > 0.0) {
    nu_e = nu_ij + nu_tilde_ij;
  }
  else {
    fn = (cn1 + Xi*Xi*Xi)/(cn1 - Xi*Xi*Xi);
    nu_e = nu_ij + fn*nu_tilde_ij;
  }
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+nu_e*proj_vector_ij)/sigma;
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

CAvgGradCorrected_TurbSA::CAvgGradCorrected_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGradCorrected_TurbSA::~CAvgGradCorrected_TurbSA(void) {
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGradCorrected_TurbSA::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_i, nVar); AD::SetPreaccIn(TurbVar_j, nVar);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5);   AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7);   AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}

CAvgGradCorrected_TurbSA_Neg::CAvgGradCorrected_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  cn1   = 16.0;
  fn    = 0.0;

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGradCorrected_TurbSA_Neg::~CAvgGradCorrected_TurbSA_Neg(void) {
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGradCorrected_TurbSA_Neg::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim);  AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_i, nVar); AD::SetPreaccIn(TurbVar_j, nVar);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  
  nu_ij = 0.5*(nu_i+nu_j);
  nu_tilde_ij = 0.5*(TurbVar_i[0]+TurbVar_j[0]);
  
  Xi = nu_tilde_ij/nu_ij;
  
  if (nu_tilde_ij > 0.0) {
    nu_e = nu_ij + nu_tilde_ij;
  }
  else {
    fn = (cn1 + Xi*Xi*Xi)/(cn1 - Xi*Xi*Xi);
    nu_e = nu_ij + fn*nu_tilde_ij;
  }
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

CSourcePieceWise_TurbSA::CSourcePieceWise_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  rotating_frame = config->GetRotating_Frame();
  
  /*--- Spalart-Allmaras closure constants ---*/
  
  cv1_3 = pow(7.1, 3.0);
  k2    = pow(0.41, 2.0);
  cb1   = 0.1355;
  cw2   = 0.3;
  ct3   = 1.2;
  ct4   = 0.5;
  cw3_6 = pow(2.0, 6.0);
  sigma = 2./3.;
  cb2   = 0.622;
  cb2_sigma = cb2/sigma;
  cw1 = cb1/k2+(1.0+cb2)/sigma;
  
}

CSourcePieceWise_TurbSA::~CSourcePieceWise_TurbSA(void) { }

void CSourcePieceWise_TurbSA::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
//  AD::StartPreacc();
//  AD::SetPreaccIn(V_i, nDim+6);
//  AD::SetPreaccIn(Vorticity_i, nDim);
//  AD::SetPreaccIn(StrainMag_i);
//  AD::SetPreaccIn(TurbVar_i[0]);
//  AD::SetPreaccIn(TurbVar_Grad_i[0], nDim);
//  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);

  if (incompressible) {
    Density_i = V_i[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
  }
  
  val_residual[0] = 0.0;
  Production      = 0.0;
  Destruction     = 0.0;
  CrossProduction = 0.0;
  val_Jacobian_i[0][0] = 0.0;
  
  /*--- Evaluate Omega ---*/
  
  Omega = sqrt(Vorticity_i[0]*Vorticity_i[0] + Vorticity_i[1]*Vorticity_i[1] + Vorticity_i[2]*Vorticity_i[2]);
  
  /*--- Rotational correction term ---*/
  
  if (rotating_frame) { Omega += 2.0*min(0.0, StrainMag_i-Omega); }
  
  if (dist_i > 1e-10) {
    
    /*--- Production term ---*/
    
    dist_i_2 = dist_i*dist_i;
    nu = Laminar_Viscosity_i/Density_i;
    Ji = TurbVar_i[0]/nu;
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    ft2 = ct3*exp(-ct4*Ji_2);
    S = Omega;
    inv_k2_d2 = 1.0/(k2*dist_i_2);
    
    Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
    Shat = max(Shat, 1.0e-10);
    inv_Shat = 1.0/Shat;
    
    /*--- Production term ---*/;

//    Original SA model
//    Production = cb1*(1.0-ft2)*Shat*TurbVar_i[0]*Volume;
    
    Production = cb1*Shat*TurbVar_i[0]*Volume;

    /*--- Destruction term ---*/
    
    r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
    g = r + cw2*(pow(r,6.0)-r);
    g_6 = pow(g,6.0);
    glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
    fw = g*glim;
    
//    Original SA model
//    Destruction = (cw1*fw-cb1*ft2/k2)*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
    
    Destruction = cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;

    /*--- Diffusion term ---*/
    
    norm2_Grad = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];
    
    CrossProduction = cb2_sigma*norm2_Grad*Volume;
    
    val_residual[0] = Production - Destruction + CrossProduction;
    
    /*--- Implicit part, production term ---*/
    
    dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
    dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
    if ( Shat <= 1.0e-10 ) dShat = 0.0;
    else dShat = (fv2+TurbVar_i[0]*dfv2)*inv_k2_d2;
    val_Jacobian_i[0][0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;
    
    /*--- Implicit part, destruction term ---*/
    
    dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
    if (r == 10.0) dr = 0.0;
    dg = dr*(1.+cw2*(6.0*pow(r,5.0)-1.0));
    dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
    val_Jacobian_i[0][0] -= cw1*(dfw*TurbVar_i[0] + 2.0*fw)*TurbVar_i[0]/dist_i_2*Volume;
    
  }

//  AD::SetPreaccOut(val_residual[0]);
//  AD::EndPreacc();
  
}

CSourcePieceWise_TurbSA_Neg::CSourcePieceWise_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar,
                                                         CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  rotating_frame = config->GetRotating_Frame();
  
  /*--- Negative Spalart-Allmaras closure constants ---*/
  
  cv1_3 = pow(7.1, 3.0);
  k2    = pow(0.41, 2.0);
  cb1   = 0.1355;
  cw2   = 0.3;
  ct3   = 1.2;
  ct4   = 0.5;
  cw3_6 = pow(2.0, 6.0);
  sigma = 2./3.;
  cb2   = 0.622;
  cb2_sigma = cb2/sigma;
  cw1 = cb1/k2+(1.0+cb2)/sigma;
  
}

CSourcePieceWise_TurbSA_Neg::~CSourcePieceWise_TurbSA_Neg(void) {
  
}

void CSourcePieceWise_TurbSA_Neg::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
//  AD::StartPreacc();
//  AD::SetPreaccIn(V_i, nDim+6);
//  AD::SetPreaccIn(Vorticity_i, nDim);
//  AD::SetPreaccIn(StrainMag_i);
//  AD::SetPreaccIn(TurbVar_i[0]);
//  AD::SetPreaccIn(TurbVar_Grad_i[0], nDim);
//  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);

  if (incompressible) {
    Density_i = V_i[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
  }
  
  val_residual[0] = 0.0;
  Production      = 0.0;
  Destruction     = 0.0;
  CrossProduction = 0.0;
  val_Jacobian_i[0][0] = 0.0;
  
  /*--- Evaluate Omega ---*/
  
  Omega = sqrt(Vorticity_i[0]*Vorticity_i[0] + Vorticity_i[1]*Vorticity_i[1] + Vorticity_i[2]*Vorticity_i[2]);

  /*--- Rotational correction term ---*/
  
  if (rotating_frame) { Omega += 2.0*min(0.0, StrainMag_i-Omega); }
  
  if (dist_i > 1e-10) {
    
    if (TurbVar_i[0] > 0.0) {
      
      /*--- Production term ---*/
      
      dist_i_2 = dist_i*dist_i;
      nu = Laminar_Viscosity_i/Density_i;
      Ji = TurbVar_i[0]/nu;
      Ji_2 = Ji*Ji;
      Ji_3 = Ji_2*Ji;
      fv1 = Ji_3/(Ji_3+cv1_3);
      fv2 = 1.0 - Ji/(1.0+Ji*fv1);
      ft2 = ct3*exp(-ct4*Ji_2);
      S = Omega;
      inv_k2_d2 = 1.0/(k2*dist_i_2);
      
      Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
      Shat = max(Shat, 1.0e-10);
      inv_Shat = 1.0/Shat;
      
      /*--- Production term ---*/;
      
      //    Original SA model
      //    Production = cb1*(1.0-ft2)*Shat*TurbVar_i[0]*Volume;
      
      Production = cb1*Shat*TurbVar_i[0]*Volume;
      
      /*--- Destruction term ---*/
      
      r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
      g = r + cw2*(pow(r,6.0)-r);
      g_6 =	pow(g,6.0);
      glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
      fw = g*glim;
      
      //    Original SA model
      //    Destruction = (cw1*fw-cb1*ft2/k2)*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
      
      Destruction = cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
      
      /*--- Diffusion term ---*/
      
      norm2_Grad = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];
      
      CrossProduction = cb2_sigma*norm2_Grad*Volume;
      
      val_residual[0] = Production - Destruction + CrossProduction;
      
      /*--- Implicit part, production term ---*/
      
      dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
      dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
      if ( Shat <= 1.0e-10 ) dShat = 0.0;
      else dShat = (fv2+TurbVar_i[0]*dfv2)*inv_k2_d2;
      val_Jacobian_i[0][0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;
      
      /*--- Implicit part, destruction term ---*/
      
      dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
      if (r == 10.0) dr = 0.0;
      dg = dr*(1.+cw2*(6.0*pow(r,5.0)-1.0));
      dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
      val_Jacobian_i[0][0] -= cw1*(dfw*TurbVar_i[0] +	2.0*fw)*TurbVar_i[0]/dist_i_2*Volume;
      
    }
    
    else {
      
      /*--- Production term ---*/
      
      dist_i_2 = dist_i*dist_i;
      
      /*--- Production term ---*/;
      
      Production = cb1*(1.0-ct3)*Omega*TurbVar_i[0]*Volume;
      
      /*--- Destruction term ---*/
      
      Destruction = cw1*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
      
      /*--- Diffusion term ---*/
      
      norm2_Grad = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];
      
      CrossProduction = cb2_sigma*norm2_Grad*Volume;
      
      val_residual[0] = Production + Destruction + CrossProduction;
      
      /*--- Implicit part, production term ---*/
      
      val_Jacobian_i[0][0] += cb1*(1.0-ct3)*Omega*Volume;
      
      /*--- Implicit part, destruction term ---*/
      
      val_Jacobian_i[0][0] += 2.0*cw1*TurbVar_i[0]/dist_i_2*Volume;
      
    }
    
  }

//  AD::SetPreaccOut(val_residual, nVar);
//  AD::EndPreacc();
}

CUpwSca_TurbSST::CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  
}

CUpwSca_TurbSST::~CUpwSca_TurbSST(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  
}

void CUpwSca_TurbSST::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(TurbVar_i,2);
  AD::SetPreaccIn(TurbVar_j,2);
  AD::SetPreaccIn(Normal, nDim);
  if (grid_movement){
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+2);
    AD::SetPreaccIn(V_j, nDim+2);

    Density_i = V_i[nDim+1];
    Density_j = V_j[nDim+1];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+3);
    AD::SetPreaccIn(V_j, nDim+3);

    Density_i = V_i[nDim+2];
    Density_j = V_j[nDim+2];
  }
  
  q_ij = 0.0;
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  
  val_residual[0] = a0*Density_i*TurbVar_i[0]+a1*Density_j*TurbVar_j[0];
  val_residual[1] = a0*Density_i*TurbVar_i[1]+a1*Density_j*TurbVar_j[1];
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;		val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[1][0] = 0.0;		val_Jacobian_i[1][1] = a0;
    
    val_Jacobian_j[0][0] = a1;		val_Jacobian_j[0][1] = 0.0;
    val_Jacobian_j[1][0] = 0.0;		val_Jacobian_j[1][1] = a1;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}

CAvgGrad_TurbSST::CAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar, su2double *constants, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma_k1  = constants[0];
  sigma_om1 = constants[2];
  sigma_k2  = constants[1];
  sigma_om2 = constants[3];
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Normal = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGrad_TurbSST::~CAvgGrad_TurbSST(void) {
  
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Normal;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGrad_TurbSST::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  su2double sigma_kine_i, sigma_kine_j, sigma_omega_i, sigma_omega_j;
  su2double diff_i_kine, diff_i_omega, diff_j_kine, diff_j_omega;
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F1_j);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute the blended constant for the viscous terms ---*/
  sigma_kine_i  = F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2;
  sigma_kine_j  = F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2;
  sigma_omega_i = F1_i*sigma_om1 + (1.0 - F1_i)*sigma_om2;
  sigma_omega_j = F1_j*sigma_om1 + (1.0 - F1_j)*sigma_om2;
  
  /*--- Compute mean effective viscosity ---*/
  diff_i_kine  = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine  = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_omega = Laminar_Viscosity_i + sigma_omega_i*Eddy_Viscosity_i;
  diff_j_omega = Laminar_Viscosity_j + sigma_omega_j*Eddy_Viscosity_j;
  
  diff_kine  = 0.5*(diff_i_kine + diff_j_kine);    // Could instead use weighted average!
  diff_omega = 0.5*(diff_i_omega + diff_j_omega);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
  }
  
  val_residual[0] = diff_kine*Proj_Mean_GradTurbVar_Corrected[0];
  val_residual[1] = diff_omega*Proj_Mean_GradTurbVar_Corrected[1];
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {
    Jacobian_i[0][0] = -diff_kine*proj_vector_ij/Density_i;		Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;									    Jacobian_i[1][1] = -diff_omega*proj_vector_ij/Density_i;
    
    Jacobian_j[0][0] = diff_kine*proj_vector_ij/Density_j; 		Jacobian_j[0][1] = 0.0;
    Jacobian_j[1][0] = 0.0;									    Jacobian_j[1][1] = diff_omega*proj_vector_ij/Density_j;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}


CAvgGradCorrected_TurbSST::CAvgGradCorrected_TurbSST(unsigned short val_nDim, unsigned short val_nVar, su2double *constants, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma_k1  = constants[0];
  sigma_om1 = constants[2];
  sigma_k2  = constants[1];
  sigma_om2 = constants[3];
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Normal = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGradCorrected_TurbSST::~CAvgGradCorrected_TurbSST(void) {
  
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Normal;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGradCorrected_TurbSST::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  su2double sigma_kine_i, sigma_kine_j, sigma_omega_i, sigma_omega_j;
  su2double diff_i_kine, diff_i_omega, diff_j_kine, diff_j_omega;
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(TurbVar_i, nVar); AD::SetPreaccIn(TurbVar_j ,nVar);
  AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F1_j);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute the blended constant for the viscous terms ---*/
  sigma_kine_i  = F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2;
  sigma_kine_j  = F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2;
  sigma_omega_i = F1_i*sigma_om1 + (1.0 - F1_i)*sigma_om2;
  sigma_omega_j = F1_j*sigma_om1 + (1.0 - F1_j)*sigma_om2;
  
  /*--- Compute mean effective viscosity ---*/
  diff_i_kine  = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine  = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_omega = Laminar_Viscosity_i + sigma_omega_i*Eddy_Viscosity_i;
  diff_j_omega = Laminar_Viscosity_j + sigma_omega_j*Eddy_Viscosity_j;
  
  diff_kine  = 0.5*(diff_i_kine + diff_j_kine);    // Could instead use weighted average!
  diff_omega = 0.5*(diff_i_omega + diff_j_omega);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  dist_ij_2 = 0.0; proj_vector_ij = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  }
  
  val_residual[0] = diff_kine*Proj_Mean_GradTurbVar_Corrected[0];
  val_residual[1] = diff_omega*Proj_Mean_GradTurbVar_Corrected[1];
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {
    Jacobian_i[0][0] = -diff_kine*proj_vector_ij/Density_i;		Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;									    Jacobian_i[1][1] = -diff_omega*proj_vector_ij/Density_i;
    
    Jacobian_j[0][0] = diff_kine*proj_vector_ij/Density_j; 		Jacobian_j[0][1] = 0.0;
    Jacobian_j[1][0] = 0.0;									    Jacobian_j[1][1] = diff_omega*proj_vector_ij/Density_j;
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

CSourcePieceWise_TurbSST::CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar, su2double *constants,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  /*--- Closure constants ---*/
  beta_star     = constants[6];
  sigma_omega_1 = constants[2];
  sigma_omega_2 = constants[3];
  beta_1        = constants[4];
  beta_2        = constants[5];
  alfa_1        = constants[8];
  alfa_2        = constants[9];
  a1            = constants[7];
}

CSourcePieceWise_TurbSST::~CSourcePieceWise_TurbSST(void) { }

void CSourcePieceWise_TurbSST::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(TurbVar_i, nVar);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);
  AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F2_i); AD::SetPreaccIn(CDkw_i);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);

  unsigned short iDim;
  su2double alfa_blended, beta_blended;
  su2double diverg, pk, pw, zeta;
  
  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5);

    Density_i = V_i[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7);

    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];
  }
  
  val_residual[0] = 0.0;        val_residual[1] = 0.0;
  val_Jacobian_i[0][0] = 0.0;	val_Jacobian_i[0][1] = 0.0;
  val_Jacobian_i[1][0] = 0.0;	val_Jacobian_i[1][1] = 0.0;
  
  //  cout<<" F1_i: "<<F1_i<<"\n";

  /*--- Computation of blended constants for the source terms---*/
  
  alfa_blended = F1_i*alfa_1 + (1.0 - F1_i)*alfa_2;
  beta_blended = F1_i*beta_1 + (1.0 - F1_i)*beta_2;
  
  if (dist_i > 1e-10) {
    
    /*--- Production ---*/
    diverg = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      diverg += PrimVar_Grad_i[iDim+1][iDim];
    
    pk = Eddy_Viscosity_i*StrainMag_i*StrainMag_i - 2.0/3.0*Density_i*TurbVar_i[0]*diverg;
    pk = min(pk,20.0*beta_star*Density_i*TurbVar_i[1]*TurbVar_i[0]);
    pk = max(pk,0.0);

    zeta = max(TurbVar_i[1], StrainMag_i*F2_i/a1);
    pw = StrainMag_i*StrainMag_i - 2.0/3.0*zeta*diverg;
    pw = max(pw,0.0);
    
    val_residual[0] += pk*Volume;
    val_residual[1] += alfa_blended*Density_i*pw*Volume;
    
    /*--- Dissipation ---*/
    val_residual[0] -= beta_star*Density_i*TurbVar_i[1]*TurbVar_i[0]*Volume;
    val_residual[1] -= beta_blended*Density_i*TurbVar_i[1]*TurbVar_i[1]*Volume;
    
    /*--- Cross diffusion ---*/
    val_residual[1] += (1.0 - F1_i)*CDkw_i*Volume;
    
    /*--- Implicit part ---*/

    // original:
    /*
    val_Jacobian_i[0][0] = -beta_star*TurbVar_i[1]*Volume;    val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[1][0] = 0.0;                               val_Jacobian_i[1][1] = -2.0*beta_blended*TurbVar_i[1]*Volume;
    */


    // swh:
    val_Jacobian_i[0][0] = -beta_star*Density_i*TurbVar_i[1]*Volume;
    val_Jacobian_i[0][1] = -beta_star*Density_i*TurbVar_i[0]*Volume;
    val_Jacobian_i[1][0] = 0.0;
    val_Jacobian_i[1][1] = -2.0*beta_blended*Density_i*TurbVar_i[1]*Volume;

  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}



//swh
CUpwSca_TurbKE::CUpwSca_TurbKE(unsigned short val_nDim, unsigned short val_nVar,
                               CConfig *config)
  :
  CNumerics(val_nDim, val_nVar, config) {

  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();

  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
}

CUpwSca_TurbKE::~CUpwSca_TurbKE(void) {

  delete [] Velocity_i;
  delete [] Velocity_j;

}

void CUpwSca_TurbKE::ComputeResidual(su2double *val_residual,
                                     su2double **val_Jacobian_i,
                                     su2double **val_Jacobian_j,
                                     CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+3);
  AD::SetPreaccIn(V_j, nDim+3);
  AD::SetPreaccIn(TurbVar_i,nVar);
  AD::SetPreaccIn(TurbVar_j,nVar);
  AD::SetPreaccIn(Normal, nDim);

  if (incompressible) {
    Density_i = V_i[nDim+1];
    Density_j = V_j[nDim+1];
  }
  else {
    Density_i = V_i[nDim+2];
    Density_j = V_j[nDim+2];
  }

  q_ij = 0.0;
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }

  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));

  val_residual[0] = a0*Density_i*TurbVar_i[0]+a1*Density_j*TurbVar_j[0];
  val_residual[1] = a0*Density_i*TurbVar_i[1]+a1*Density_j*TurbVar_j[1];
  val_residual[2] = a0*Density_i*TurbVar_i[2]+a1*Density_j*TurbVar_j[2];
  val_residual[3] = 0.0; // no convection in f scalar


  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[0][2] = 0.0;
    val_Jacobian_i[0][3] = 0.0;

    val_Jacobian_i[1][0] = 0.0;
    val_Jacobian_i[1][1] = a0;
    val_Jacobian_i[1][2] = 0.0;
    val_Jacobian_i[1][3] = 0.0;

    val_Jacobian_i[2][0] = 0.0;
    val_Jacobian_i[2][1] = 0.0;
    val_Jacobian_i[2][2] = a0;
    val_Jacobian_i[2][3] = 0.0;

    val_Jacobian_i[3][0] = 0.0;
    val_Jacobian_i[3][1] = 0.0;
    val_Jacobian_i[3][2] = 0.0;
    val_Jacobian_i[3][3] = 0.0;


    val_Jacobian_j[0][0] = a1;
    val_Jacobian_j[0][1] = 0.0;
    val_Jacobian_j[0][2] = 0.0;
    val_Jacobian_j[0][3] = 0.0;

    val_Jacobian_j[1][0] = 0.0;
    val_Jacobian_j[1][1] = a1;
    val_Jacobian_j[1][2] = 0.0;
    val_Jacobian_j[1][3] = 0.0;

    val_Jacobian_j[2][0] = 0.0;
    val_Jacobian_j[2][1] = 0.0;
    val_Jacobian_j[2][2] = a1;
    val_Jacobian_j[2][3] = 0.0;

    val_Jacobian_j[3][0] = 0.0;
    val_Jacobian_j[3][1] = 0.0;
    val_Jacobian_j[3][2] = 0.0;
    val_Jacobian_j[3][3] = 0.0;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

CAvgGrad_TurbKE::CAvgGrad_TurbKE(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 su2double *constants,
                                 CConfig *config)
  :
  CNumerics(val_nDim, val_nVar, config) {

  unsigned short iVar;

  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  sigma_k = constants[1];
  sigma_e = constants[2];
  sigma_z = constants[3];

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Normal = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];

}

CAvgGrad_TurbKE::~CAvgGrad_TurbKE(void) {

  unsigned short iVar;

  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Normal;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
  delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;

}

void CAvgGrad_TurbKE::ComputeResidual(su2double *val_residual,
                                      su2double **Jacobian_i,
                                      su2double **Jacobian_j,
                                      CConfig *config) {

  su2double sigma_kine_i, sigma_kine_j, sigma_epsi_i, sigma_epsi_j;
  su2double sigma_zeta_i, sigma_zeta_j;
  su2double diff_i_kine, diff_j_kine;
  su2double diff_i_epsi, diff_j_epsi;
  su2double diff_i_zeta, diff_j_zeta;
  su2double diff_i_f, diff_j_f;

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(Lm_i); AD::SetPreaccIn(Lm_j);
  AD::SetPreaccIn(Volume);

  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }

  /*--- Compute the blended constant for the viscous terms ---*/
  // there are already stored as inverses
  sigma_kine_i = sigma_k;
  sigma_kine_j = sigma_k;
  sigma_epsi_i = sigma_e;
  sigma_epsi_j = sigma_e;
  sigma_zeta_i = sigma_z;
  sigma_zeta_j = sigma_z;

  /*--- Compute mean effective viscosity ---*/
  diff_i_kine = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_epsi = Laminar_Viscosity_i + sigma_epsi_i*Eddy_Viscosity_i;
  diff_j_epsi = Laminar_Viscosity_j + sigma_epsi_j*Eddy_Viscosity_j;
  diff_i_zeta = Laminar_Viscosity_i + sigma_zeta_i*Eddy_Viscosity_i;
  diff_j_zeta = Laminar_Viscosity_j + sigma_zeta_j*Eddy_Viscosity_j;

  diff_kine = 0.5*(diff_i_kine + diff_j_kine);    // Could instead use weighted average!
  diff_epsi = 0.5*(diff_i_epsi + diff_j_epsi);
  diff_zeta = 0.5*(diff_i_zeta + diff_j_zeta);
  //  diff_f = Lm_i*Lm_i; //here
  diff_f = 1.0; //Density_i;
  //  cout << "lm_i:" << Lm_i << "\n";

  /*--- Compute vector going from iPoint to jPoint ---*/
  su2double n_mag=0.0;
  su2double s_mag=0.0;
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
    n_mag += Normal[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;

  s_mag = sqrt(dist_ij_2);
  n_mag = sqrt(n_mag);

  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
  }

  val_residual[0] = diff_kine*Proj_Mean_GradTurbVar_Corrected[0];
  val_residual[1] = diff_epsi*Proj_Mean_GradTurbVar_Corrected[1];
  val_residual[2] = diff_zeta*Proj_Mean_GradTurbVar_Corrected[2];
  val_residual[3] = diff_f*Proj_Mean_GradTurbVar_Corrected[3];

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/ //here
  if (implicit) {
    Jacobian_i[0][0] = -diff_kine*proj_vector_ij/Density_i;
    Jacobian_i[0][1] = 0.0;
    Jacobian_i[0][2] = 0.0;
    Jacobian_i[0][3] = 0.0;

    Jacobian_i[1][0] = 0.0;
    Jacobian_i[1][1] = -diff_epsi*proj_vector_ij/Density_i;
    Jacobian_i[1][2] = 0.0;
    Jacobian_i[1][3] = 0.0;

    Jacobian_i[2][0] = 0.0;
    Jacobian_i[2][1] = 0.0;
    Jacobian_i[2][2] = -diff_zeta*proj_vector_ij/Density_i;
    Jacobian_i[2][3] = 0.0;

    Jacobian_i[3][0] = 0.0;
    Jacobian_i[3][1] = 0.0;
    Jacobian_i[3][2] = 0.0;
    Jacobian_i[3][3] = -diff_f*proj_vector_ij;

    Jacobian_j[0][0] = diff_kine*proj_vector_ij/Density_j;
    Jacobian_j[0][1] = 0.0;
    Jacobian_j[0][2] = 0.0;
    Jacobian_j[0][3] = 0.0;

    Jacobian_j[1][0] = 0.0;
    Jacobian_j[1][1] = diff_epsi*proj_vector_ij/Density_j;
    Jacobian_j[1][2] = 0.0;
    Jacobian_j[1][3] = 0.0;

    Jacobian_j[2][0] = 0.0;
    Jacobian_j[2][1] = 0.0;
    Jacobian_j[2][2] = diff_zeta*proj_vector_ij/Density_j;
    Jacobian_j[2][3] = 0.0;

    Jacobian_j[3][0] = 0.0;
    Jacobian_j[3][1] = 0.0;
    Jacobian_j[3][2] = 0.0;
    Jacobian_j[3][3] = diff_f*proj_vector_ij;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}


CAvgGradCorrected_TurbKE::CAvgGradCorrected_TurbKE(unsigned short val_nDim,
                                                   unsigned short val_nVar,
                                                   su2double *constants,
                                                   CConfig *config)
  :
  CNumerics(val_nDim, val_nVar, config) {

  unsigned short iVar;

  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  sigma_k = constants[1];
  sigma_e = constants[2];
  sigma_z = constants[3];

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Normal = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];

}

CAvgGradCorrected_TurbKE::~CAvgGradCorrected_TurbKE(void) {

  unsigned short iVar;

  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Normal;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;

}

void CAvgGradCorrected_TurbKE::ComputeResidual(su2double *val_residual,
                                               su2double **Jacobian_i,
                                               su2double **Jacobian_j,
                                               CConfig *config) {

  su2double sigma_kine_i, sigma_kine_j, sigma_epsi_i, sigma_epsi_j, sigma_zeta_i, sigma_zeta_j;
  su2double diff_i_kine, diff_i_epsi, diff_j_kine, diff_j_epsi;
  su2double diff_i_zeta, diff_j_zeta, diff_i_f, diff_j_f;

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(TurbVar_i, nVar); AD::SetPreaccIn(TurbVar_j ,nVar);
  AD::SetPreaccIn(Lm_i); AD::SetPreaccIn(Lm_j);
  AD::SetPreaccIn(Volume);

  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }

  /*--- Compute the blended constant for the viscous terms ---*/
  sigma_kine_i = sigma_k;
  sigma_kine_j = sigma_k;
  sigma_epsi_i = sigma_e;
  sigma_epsi_j = sigma_e;
  sigma_zeta_i = sigma_z;
  sigma_zeta_j = sigma_z;

  /*--- Compute mean effective viscosity ---*/
  diff_i_kine = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_epsi = Laminar_Viscosity_i + sigma_epsi_i*Eddy_Viscosity_i;
  diff_j_epsi = Laminar_Viscosity_j + sigma_epsi_j*Eddy_Viscosity_j;
  diff_i_zeta = Laminar_Viscosity_i + sigma_zeta_i*Eddy_Viscosity_i;
  diff_j_zeta = Laminar_Viscosity_j + sigma_zeta_j*Eddy_Viscosity_j;

  diff_kine = 0.5*(diff_i_kine + diff_j_kine);    // Could instead use weighted average!
  diff_epsi = 0.5*(diff_i_epsi + diff_j_epsi);
  diff_zeta = 0.5*(diff_i_zeta + diff_j_zeta);
  //  diff_f = Lm_i*Lm_i; //here
  diff_f = 1.0; //Density_i;
  //  cout << "lm_i:" << Lm_i << "\n";

  /*--- Compute vector going from iPoint to jPoint ---*/
  su2double n_mag=0.0;
  su2double s_mag=0.0;
  dist_ij_2 = 0.0; proj_vector_ij = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
    n_mag += Normal[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;

  s_mag = sqrt(dist_ij_2);
  n_mag = sqrt(n_mag);

  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
      (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;

    // NB: Jacobian corresponds to just this part
    //Proj_Mean_GradTurbVar_Corrected[iVar] = (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;

  }

  val_residual[0] = diff_kine*Proj_Mean_GradTurbVar_Corrected[0];
  val_residual[1] = diff_epsi*Proj_Mean_GradTurbVar_Corrected[1];
  val_residual[2] = diff_zeta*Proj_Mean_GradTurbVar_Corrected[2];
  val_residual[3] = diff_f*Proj_Mean_GradTurbVar_Corrected[3];

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {

    Jacobian_i[0][0] = -diff_kine*proj_vector_ij/Density_i;
    Jacobian_i[0][1] = 0.0;
    Jacobian_i[0][2] = 0.0;
    Jacobian_i[0][3] = 0.0;

    Jacobian_i[1][0] = 0.0;
    Jacobian_i[1][1] = -diff_epsi*proj_vector_ij/Density_i;
    Jacobian_i[1][2] = 0.0;
    Jacobian_i[1][3] = 0.0;

    Jacobian_i[2][0] = 0.0;
    Jacobian_i[2][1] = 0.0;
    Jacobian_i[2][2] = -diff_zeta*proj_vector_ij/Density_i;
    Jacobian_i[2][3] = 0.0;

    Jacobian_i[3][0] = 0.0;
    Jacobian_i[3][1] = 0.0;
    Jacobian_i[3][2] = 0.0;
    Jacobian_i[3][3] = -diff_f*proj_vector_ij;
    Jacobian_i[3][3] -= 0.5 * diff_f * n_mag/Volume;
    Jacobian_i[3][3] += 0.5 * diff_f * proj_vector_ij*s_mag/Volume;


    Jacobian_j[0][0] = diff_kine*proj_vector_ij/Density_j;
    Jacobian_j[0][1] = 0.0;
    Jacobian_j[0][2] = 0.0;
    Jacobian_j[0][3] = 0.0;

    Jacobian_j[1][0] = 0.0;
    Jacobian_j[1][1] = diff_epsi*proj_vector_ij/Density_j;
    Jacobian_j[1][2] = 0.0;
    Jacobian_j[1][3] = 0.0;

    Jacobian_j[2][0] = 0.0;
    Jacobian_j[2][1] = 0.0;
    Jacobian_j[2][2] = diff_zeta*proj_vector_ij/Density_j;
    Jacobian_j[2][3] = 0.0;

    Jacobian_j[3][0] = 0.0;
    Jacobian_j[3][1] = 0.0;
    Jacobian_j[3][2] = 0.0;
    Jacobian_j[3][3] = diff_f*proj_vector_ij;
    Jacobian_j[3][3] += 0.5 * diff_f * n_mag/Volume;
    Jacobian_j[3][3] -= 0.5 * diff_f * proj_vector_ij*s_mag/Volume;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

CSourcePieceWise_TurbKE::CSourcePieceWise_TurbKE(unsigned short val_nDim,
                                                 unsigned short val_nVar,
                                                 su2double *constants,
                                                 CConfig *config)
  :
  CNumerics(val_nDim, val_nVar, config) {

  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  /*--- Closure constants ---*/
  /*sigma_k = constants[0];
  sigma_e = constants[1];
  C_mu    = constants[2];
  C_e1    = constants[3];
  C_e2    = constants[4];*/

  /* RNG
  C_mu    = constants[0];
  sigma_k = constants[1];
  sigma_e = constants[2];
  C_e1    = constants[4];
  C_e2    = constants[5];
  eta_o   = constants[6];
  beta    = constants[7];
  */

  /* zeta f */
  C_mu    = constants[0];
  sigma_k = constants[1];
  sigma_e = constants[2];
  sigma_z = constants[3];
  C_e1o   = constants[4];
  C_e2    = constants[5];
  C_1     = constants[6];
  C_2p    = constants[7];
  C_T     = constants[8];
  C_L     = constants[9];
  C_eta   = constants[10];

}

CSourcePieceWise_TurbKE::~CSourcePieceWise_TurbKE(void) { }

void CSourcePieceWise_TurbKE::ComputeResidual(su2double *val_residual,
                                              su2double **val_Jacobian_i,
                                              su2double **val_Jacobian_j,
                                              CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+7);
  AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(TurbVar_i, nVar);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);

  // Pull variables out of V_i
  if (incompressible) {
    Density_i = V_i[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];
  }

  // for readability...
  const su2double tke = TurbVar_i[0];
  const su2double tdr = TurbVar_i[1];
  const su2double v20 = TurbVar_i[2];
  const su2double f   = TurbVar_i[3];

  const su2double tke_lim = max(tke, 1e-8); // 1e-8 arbitrary
  const su2double tdr_lim = max(tdr, 36.0*Laminar_Viscosity_i/Density_i); // 36*lam visc ensures T3 <= 1

  // if (tke < 0.0) {
  //   std::cout << "tke = " << tke
  //             << " at x = " << Coord_i[0] << ", " << Coord_i[1]
  //             << std::endl;
  // }

  // if (tdr < 0.0) {
  //   std::cout << "tdr = " << tdr
  //             << " at x = " << Coord_i[0] << ", " << Coord_i[1]
  //             << std::endl;
  // }

  // if (v20 < -1e-15) {
  //   std::cout << "v2  = " << v20
  //             << " at x = " << Coord_i[0] << ", " << Coord_i[1]
  //             << std::endl;
  // }

  // if (f < -1e-15) {
  //   std::cout << "f   = " << f
  //             << " at x = " << Coord_i[0] << ", " << Coord_i[1]
  //             << std::endl;
  // }

  // make sure v2 is well-behaved
  const su2double scale = 1.0e-8;
  //su2double zeta = max(v20/tke, scale);
  su2double zeta = max(v20/tke_lim, scale);
  zeta = min(zeta,2.0/3.0);

  const su2double v2 = max(v20, zeta*tke);

  const su2double rho = Density_i;

  const su2double mu  = Laminar_Viscosity_i;
  const su2double muT = Eddy_Viscosity_i;

  const su2double nu  = mu/rho;
  const su2double nuT = muT/rho;

  const su2double S   = StrainMag_i; //*sqrt(2.0) already included
  const su2double Vol = Volume;

  su2double diverg = 0.0;
  for (unsigned int iDim = 0; iDim < nDim; iDim++)
    diverg += PrimVar_Grad_i[iDim+1][iDim];

  //----------------------------------------------------------------------------
  // Scream if tke, tdr, etc are out of bounds?
  //
  // if(tke < 0.0) std::cout << "WTF!?! k is negative!!!" << std::endl;
  // if(tdr < 0.0) std::cout << "WTF!?! epsilon is negative!!!" << std::endl;
  //  if (tke>=20.0) {cout << "TKE: " << tke << "\n";}
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Impose some bounds on turb variables?
  //
  // // denominator floors
  // su2double VelMag, *VelInf, L_Inf, scalar_min;
  // su2double tke_raw, tdr_raw, v2_raw;
  // tke_raw = tke;
  // tdr_raw = tdr;
  // v2_raw = v2;
  // VelInf = config->GetVelocity_FreeStreamND();
  // L_Inf = config->GetLength_Reynolds();
  // for (iDim = 0; iDim < nDim; iDim++)
  // VelMag += VelInf[iDim]*VelInf[iDim];
  // VelMag = sqrt(VelMag);

  // su2double solve_tol = config->GetLinear_Solver_Error();
  // su2double Re = config->GetReynolds();
  // su2double iRe = 1.0/Re;
  // su2double scale;
  // scale = 1.0e-8;
  // zeta = max(v2_raw/tke_raw, scale);
  // zeta = min(zeta,2.0/3.0);
  // scalar_min = scale/(VelMag*VelMag); // setting based on tke min being 1e-8
  // tke = max(tke, scalar_min*VelMag*VelMag);
  // tdr = max(tdr, scalar_min*VelMag*VelMag*VelMag/L_Inf);
  // //v2 = max(v2, 2.0/3.0*scalar_min*VelMag*VelMag);
  // f = max(f, scalar_min*VelMag/L_Inf);

  // v2 = max(v2, zeta*tke);
  // //  S = max(S,scalar_min*VelMag/L_Inf); // no checked...
  // S = max(S,1.0E-14);
  //----------------------------------------------------------------------------

  // NB: We determine time and length scales here due to Jacobian branching
  // TODO: Could replace max and min with differentiable approximations
  // TODO: Write a function that evaluates T and L along with derivatives

  //--- Model time scale ---//
  // const su2double T1     = tke/tdr;
  // const su2double T1_rk  =  1.0/(rho*tdr);
  // const su2double T1_re  = - T1/(rho*tdr);
  // const su2double T1_rv2 = 0.0;

  //const su2double T1     = tke_lim/tdr;
  const su2double T1     = tke_lim/tdr_lim;
  const su2double T1_rk  =  1.0/(rho*tdr);
  const su2double T1_re  = - T1/(rho*tdr);
  const su2double T1_rv2 = 0.0;


  const su2double T2     = 1.0E14; //0.6/(sqrt(6.0)*C_mu*S*zeta);
  const su2double T2_rk  = 0.0;
  const su2double T2_re  = 0.0;
  const su2double T2_rv2 = 0.0;

  //const su2double T3     = C_T*sqrt(nu/tdr);
  const su2double T3     = C_T*sqrt(nu/tdr_lim);
  const su2double T3_rk  = 0.0;
  const su2double T3_re  = -0.5*T3/(rho*tdr);
  const su2double T3_rv2 = 0.0;

  // T = max(min(T1,T2),T3)
  su2double T     = T1;
  su2double T_rk  = T1_rk;
  su2double T_re  = T1_re;
  su2double T_rv2 = T1_rv2;

  // Use smooth version of maximum?
  // // T = smooth_max(T1,T3)
  // const su2double del  = T1 - T3;
  // const su2double sabs = sqrt( del*del + 1.0 );

  // T = 0.5*(T1 + T3 + sabs );

  // T_rk  = 0.5*(T1_rk  + T3_rk  + (T1_rk  - T3_rk )*del/sabs );
  // T_re  = 0.5*(T1_re  + T3_re  + (T1_re  - T3_re )*del/sabs );
  // T_rv2 = 0.5*(T1_rv2 + T3_rv2 + (T1_rv2 - T3_rv2)*del/sabs );

  // Use maximum?
  if (T>T2) {
    T = T2;
    T_rk = T2_rk; T_re = T2_re; T_rv2 = T2_rv2;
  }

  if (T<T3) {
    T = T3;
    T_rk = T3_rk; T_re = T3_re; T_rv2 = T3_rv2;
  }

  const su2double Tsq = T*T;


  //--- Model length scale ---//
  //const su2double L1     = pow(tke,1.5)/tdr;
  const su2double L1     = pow(tke_lim,1.5)/tdr_lim;
  const su2double L1_rk  =    -L1/(rho*tke);
  const su2double L1_re  = 1.5*L1/(rho*tdr);
  const su2double L1_rv2 = 0.0;

  const su2double L2     = 1.0E14; //sqrt(tke)/(sqrt(6.0)*C_mu*S*zeta);
  const su2double L2_rk  = 0.0;
  const su2double L2_re  = 0.0;
  const su2double L2_rv2 = 0.0;

  //const su2double L3     = C_eta*pow(pow(nu,3.0)/tdr,0.25);
  const su2double L3     = C_eta*pow(pow(nu,3.0)/tdr_lim,0.25);
  const su2double L3_rk  = 0.0;
  const su2double L3_re  = -0.25*L3/(rho*tdr);
  const su2double L3_rv2 = 0.0;

  // L = max(min(L1,L2),L3)... mult by C_L below
  su2double L     = L1;
  su2double L_rk  = L1_rk;
  su2double L_re  = L1_re;
  su2double L_rv2 = L1_rv2;

  if (L>L2) {
    L = L2;
    L_rk = L2_rk; L_re = L2_re; L_rv2 = L2_rv2;
  }

  if (L<L3) {
    L = L3;
    L_rk = L3_rk; L_re = L3_re; L_rv2 = L3_rv2;
  }

  L *= C_L;
  L_rk *= C_L; L_re *= C_L; L_rv2 *= C_L;

  const su2double Lsq = L*L;


  //--- v2-f ---//

  // 4 equations.  For each equation, we identify production and
  // dissipation terms.  This is somewhat artificial for f.  The
  // Jacobian is abused to help keep k and epsilon positive.


  // TKE equation...
  su2double Pk, Pk_rk, Pk_re, Pk_rv2;
  su2double Dk, Dk_rk, Dk_re, Dk_rv2;

  //... production
  // NB: we ignore the jacobian of production here

  const su2double xx = Coord_i[0];
  const su2double xfac = 0.5*(1.0 + tanh(20.0*(xx-0.3)));
  //Pk     = (muT*S*S - 2.0/3.0*rho*tke*diverg)*xfac;
  Pk     = muT*S*S - 2.0/3.0*rho*tke*diverg;
  //Pk     = muT*S*S;
  //Pk     = 0.0; // Should be very robust with production off

  Pk_rk  = 0.0;
  //Pk_rk  = min(C_mu*zeta*T*S*S, 0.9/(Vol*TimeStep));
  //Pk_rk  = C_mu*zeta*T*S*S;
  Pk_re  = 0.0;
  Pk_rv2 = 0.0;

  // //... dissipation
  // Dk     = rho*tdr;

  // Dk_rk  = 0.0;
  // Dk_re  = 1.0;
  // Dk_rv2 = 0.0;

  //... dissipation
  Dk     = rho*tke/T1;

  Dk_rk  = 1.0/T1;
  Dk_re  = 0.0;
  Dk_rv2 = 0.0;


  // Dissipation equation...
  su2double Pe, Pe_rk, Pe_re, Pe_rv2;
  su2double De, De_rk, De_re, De_rv2;

  // NB: C_e1 depends on tke and v2 in v2-f
  //const su2double C_e1 = C_e1o*(1.0+0.045*sqrt(tke/v2));
  const su2double C_e1 = C_e1o*(1.0+0.045*sqrt(1.0/zeta));
  //const su2double C_e1 = C_e1o;

  // ... production
  //Pe = C_e1*Pk/T;
  Pe = C_e1*C_mu*rho*v2*S*S;

  Pe_rk  = 0.0;
  Pe_re  = 0.0;
  Pe_rv2 = 0.0;

  // ... dissipation
  // De = C_e2*rho*tdr/T;

  // De_rk  =        - C_e2*rho*tdr*T_rk /Tsq;
  // De_re  = C_e2/T - C_e2*rho*tdr*T_re /Tsq;
  // De_rv2 =        - C_e2*rho*tdr*T_rv2/Tsq;

  // De = C_e2*rho*tdr/T1;

  // De_rk  =         - C_e2*rho*tdr*T1_rk /Tsq;
  // De_re  = C_e2/T1 - C_e2*rho*tdr*T1_re /Tsq;
  // De_rv2 =         - C_e2*rho*tdr*T1_rv2/Tsq;

  De = C_e2*rho*tdr/T;

  De_rk  = 0.0;
  De_re  = C_e2/T;
  De_rv2 = 0.0;


  // v2 equation...
  su2double Pv2, Pv2_rk, Pv2_re, Pv2_rv2, Pv2_f;
  su2double Dv2, Dv2_rk, Dv2_re, Dv2_rv2, Dv2_f;

  // ... production
  Pv2 = rho*tke*f;

  Pv2_rk  = 0.0; //f;
  Pv2_rk  = 0.0;
  Pv2_re  = 0.0;
  Pv2_rv2 = 0.0; //min(f/zeta, 1.0/(Vol*TimeStep)); //0.0;
  Pv2_f   = 0.0; //rho*tke; // keep this?

  // // ... dissipation
  // Dv2     =  6.0*(v2/tke)*rho*tdr;

  // Dv2_rk  = -6.0*(v2/tke)*(tdr/tke);
  // Dv2_re  =  6.0*(v2/tke);
  // Dv2_rv2 =  6.0*(tdr/tke);
  // Dv2_f   =  0.0;

  // ... dissipation
  Dv2     =  6.0*rho*v2/T1;

  Dv2_rk  = 0.0;
  Dv2_re  = 0.0;
  Dv2_rv2 = 6.0/T1;
  Dv2_f   = 0.0;


  // f equation...
  su2double Pf;
  su2double Df, Df_f;

  //... production
  const su2double C1m6 = C_1 - 6.0;
  const su2double ttC1m1 = (2.0/3.0)*(C_1 - 1.0);
  //const su2double C_2f = C_2p + 0.5*(2.0/3.0-C_2p)*(1.0+tanh(50.0*(v2/tke-0.55)));
  const su2double C_2f = C_2p + 0.5*(2.0/3.0-C_2p)*(1.0+tanh(50.0*(zeta-0.55)));

  //Pf = (C_2f*Pk/tke - (C1m6*v2/tke - ttC1m1)*rho/T) / Lsq;
  //Pf = (C_2f*Pk/(rho*tke) - (C1m6*v2/tke - ttC1m1)/T) / Lsq;
  //Pf = (C_2f*Pk/(rho*tke) - (C1m6*(2.0/3.0) - ttC1m1)/T) / Lsq;
  //Pf = ( - (C1m6*(2.0/3.0) - ttC1m1)/T) / Lsq;
  //Pf = (C_2f*Pk/(rho*tke_lim) - (C1m6*(2.0/3.0) - ttC1m1)/T) / Lsq;
  //Pf = (C_2f*Pk/(rho*tke_lim) - (C1m6*zeta - ttC1m1)/T) / Lsq;
  //Pf = (C_2f*Pk/(rho*tke_lim) - (C1m6*(2.0/3.0) - ttC1m1)/T) / Lsq;
  Pf = (C_2f*C_mu*zeta*T*S*S - (C1m6*zeta - ttC1m1)/T) / Lsq;
  //Pf = 0.0;

  // not keeping any derivatives of Pf

  //... dissipation
  Df = f/Lsq;

  Df_f = 1.0/Lsq;


  // check for nans
  bool found_nan = (std::isnan(Pk)  || std::isnan(Dk)  ||
                    std::isnan(Pe)  || std::isnan(De)  ||
                    std::isnan(Pv2) || std::isnan(Dv2) ||
                    std::isnan(Pf)  || std::isnan(Df)  ||
                    std::isnan(Pk_rk)  || std::isnan(Pk_re)  || std::isnan(Pk_rv2)  ||
                    std::isnan(Pe_rk)  || std::isnan(Pe_re)  || std::isnan(Pe_rv2)  ||
                    std::isnan(Pv2_rk) || std::isnan(Pv2_re) || std::isnan(Pv2_rv2) ||
                    std::isnan(Dk_rk)  || std::isnan(Dk_re)  || std::isnan(Dk_rv2)  ||
                    std::isnan(De_rk)  || std::isnan(De_re)  || std::isnan(De_rv2)  ||
                    std::isnan(Dv2_rk) || std::isnan(Dv2_re) || std::isnan(Dv2_rv2) );

  if (found_nan && Coord_i[1] > 1e-14) {
    std::cout << "WTF!?! Found a nan at x = " << Coord_i[0] << ", " << Coord_i[1] << std::endl;
    std::cout << "turb state = " << tke << ", " << tdr << ", " << v2 << ", " << f << std::endl;
    std::cout << "T1         = " << T1 << ", T3 " << T3 << std::endl;
    std::cout << "T          = " << T  << ", C_e1 = " << C_e1 << std::endl;
    std::cout << "TKE eqn    = " << Pk << " - " << Dk << std::endl;
    std::cout << "TDR eqn    = " << Pe << " - " << De << std::endl;
    std::cout << "v2  eqn    = " << Pv2 << " - " << Dv2 << std::endl;
  }


  // form source term and Jacobian...

  // TKE
  val_residual[0]      = (Pk      - Dk     ) * Vol;

  val_Jacobian_i[0][0] = (Pk_rk   - Dk_rk  ) * Vol;
  val_Jacobian_i[0][1] = (Pk_re   - Dk_re  ) * Vol;
  val_Jacobian_i[0][2] = (Pk_rv2  - Dk_rv2 ) * Vol;
  val_Jacobian_i[0][3] = 0.0;

  // Dissipation
  val_residual[1]      = (Pe      - De     ) * Vol;

  val_Jacobian_i[1][0] = (Pe_rk   - De_rk  ) * Vol;
  val_Jacobian_i[1][1] = (Pe_re   - De_re  ) * Vol;
  val_Jacobian_i[1][2] = (Pe_rv2  - De_rv2 ) * Vol;
  val_Jacobian_i[1][3] = 0.0;

  // v2
  val_residual[2]      = (Pv2     - Dv2    ) * Vol;

  val_Jacobian_i[2][0] = (Pv2_rk  - Dv2_rk ) * Vol;
  val_Jacobian_i[2][1] = (Pv2_re  - Dv2_re ) * Vol;
  val_Jacobian_i[2][2] = (Pv2_rv2 - Dv2_rv2) * Vol;
  val_Jacobian_i[2][3] = (Pv2_f   - Dv2_f  ) * Vol;

  // f
  val_residual[3]      = (Pf      - Df     ) * Vol;

  val_Jacobian_i[3][0] = 0.0;
  val_Jacobian_i[3][1] = 0.0;
  val_Jacobian_i[3][2] = 0.0;
  val_Jacobian_i[3][3] = (        - Df_f   ) * Vol;


  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}



CUpwSca_TurbML::CUpwSca_TurbML(unsigned short val_nDim, unsigned short val_nVar,
                               CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  
}

CUpwSca_TurbML::~CUpwSca_TurbML(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  
}

void CUpwSca_TurbML::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  
  q_ij = 0.0;
  
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  } else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  val_residual[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_j[0][0] = a1;
  }
  
  
}

CAvgGrad_TurbML::CAvgGrad_TurbML(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGrad_TurbML::~CAvgGrad_TurbML(void) {
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGrad_TurbML::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+nu_e*proj_vector_ij)/sigma;
  }
  
}

CAvgGradCorrected_TurbML::CAvgGradCorrected_TurbML(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iVar;
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGradCorrected_TurbML::~CAvgGradCorrected_TurbML(void) {
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGradCorrected_TurbML::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
  }
  
}
