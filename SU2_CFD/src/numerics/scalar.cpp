/*!
 * \file scalar.cpp
 * \brief This file contains the numerical methods for scalar transport eqns.
 * \author T. Economon, D. Mayer
 * \version 6.1.0 "Falcon"
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

#include "../../include/numerics/scalar.hpp"

CUpwtransportedScalar::CUpwtransportedScalar(unsigned short val_nDim,
                       unsigned short val_nVar,
                       const CConfig* config) :
  CNumerics(val_nDim, val_nVar, config),
  // FIXME dan: this has to work for turbulence scalars and others. 
  //implicit(config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT)
  implicit(config->GetKind_TimeIntScheme_Scalar() == EULER_IMPLICIT),
  incompressible(config->GetKind_Regime() == INCOMPRESSIBLE),
  dynamic_grid(config->GetDynamic_Grid())
{
  Flux = new su2double [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CUpwtransportedScalar::~CUpwtransportedScalar(void) {

  delete [] Flux;
  if (Jacobian_i != nullptr) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete [] Jacobian_i[iVar];
      delete [] Jacobian_j[iVar];
    }
    delete [] Jacobian_i;
    delete [] Jacobian_j;
  }
}

CNumerics::ResidualType<> CUpwtransportedScalar::ComputeResidual(const CConfig* config) {

  unsigned short iDim;

  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);

  // FIXME daniel: this has to work for TurbVar and for scalars
  // AD::SetPreaccIn(TurbVar_i, nVar);  AD::SetPreaccIn(TurbVar_j, nVar);

  AD::SetPreaccIn(scalar_i, nVar);  AD::SetPreaccIn(scalar_j, nVar);

  if (dynamic_grid) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }

  ExtraADPreaccIn();

  Density_i = V_i[nDim+2];
  Density_j = V_j[nDim+2];

  q_ij = 0.0;
  if (dynamic_grid) {
    for (iDim = 0; iDim < nDim; iDim++) {
      su2double Velocity_i = V_i[iDim+1] - GridVel_i[iDim];
      su2double Velocity_j = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i+Velocity_j)*Normal[iDim];
    }
  }
  else {
    for (iDim = 0; iDim < nDim; iDim++) {
      q_ij += 0.5*(V_i[iDim+1]+V_j[iDim+1])*Normal[iDim];
    }
  }

  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));

  FinishResidualCalc(config);

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

CUpwtransportedScalar_General::CUpwtransportedScalar_General(unsigned short val_nDim,
                                       unsigned short val_nVar,
                                       CConfig *config)
: CUpwtransportedScalar(val_nDim, val_nVar, config) { }

CUpwtransportedScalar_General::~CUpwtransportedScalar_General(void) { }

void CUpwtransportedScalar_General::ExtraADPreaccIn() {
  AD::SetPreaccIn(scalar_i, nVar);  AD::SetPreaccIn(scalar_j, nVar);
  AD::SetPreaccIn(V_i, nDim+3); AD::SetPreaccIn(V_j, nDim+3);
}

void CUpwtransportedScalar_General::FinishResidualCalc(su2double *val_residual,
                                            su2double **val_Jacobian_i,
                                            su2double **val_Jacobian_j,
                                            CConfig *config) {
  
  unsigned short iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = (a0*Density_i*scalar_i[iVar] +
                          a1*Density_j*scalar_j[iVar]);
    if (implicit) {
      for (jVar = 0; jVar < nVar; jVar++) {
        if (iVar == jVar) {
          val_Jacobian_i[iVar][jVar] = a0*Density_i;
          val_Jacobian_j[iVar][jVar] = a1*Density_j;
        } else {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
  }
  
}

void CUpwtransportedScalar_General::FinishResidualCalc(const CConfig* config) {

  unsigned short iVar, jVar;

  for (iVar = 0; iVar < nVar; iVar++) {

    Flux[iVar] = a0*Density_i*scalar_i[iVar] + a1*Density_j*scalar_j[iVar];

    if (implicit) {
      for (jVar = 0; jVar < nVar; jVar++) {
        if (iVar == jVar) {
          Jacobian_i[iVar][jVar] = a0*Density_i;
          Jacobian_j[iVar][jVar] = a1*Density_j;
        } else {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
  }

}

CAvgGradtransportedScalar::CAvgGradtransportedScalar(unsigned short val_nDim,
                               unsigned short val_nVar,
                               bool correct_grad,
                               CConfig *config)
: CNumerics(val_nDim, val_nVar, config), correct_gradient(correct_grad) {
  
  implicit       = (config->GetKind_TimeIntScheme_Scalar() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  Edge_Vector = new su2double[nDim];
  
  Proj_Mean_GradScalarVar_Normal = new su2double[nVar];
  Proj_Mean_GradScalarVar_Edge   = new su2double[nVar];
  Proj_Mean_GradScalarVar        = new su2double[nVar];
  
  Mean_GradScalarVar = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradScalarVar[iVar] = new su2double[nDim];
  
}

CAvgGradtransportedScalar::~CAvgGradtransportedScalar(void) {
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradScalarVar_Normal;
  delete [] Proj_Mean_GradScalarVar_Edge;
  delete [] Proj_Mean_GradScalarVar;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradScalarVar[iVar];
  delete [] Mean_GradScalarVar;
  
}

void CAvgGradtransportedScalar::ComputeResidual(su2double *val_residual,
                                     su2double **Jacobian_i,
                                     su2double **Jacobian_j,
                                     CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(scalar_grad_i, nVar, nDim);
  AD::SetPreaccIn(scalar_grad_j, nVar, nDim);
  AD::SetPreaccIn(Diffusion_Coeff_i, nVar);
  AD::SetPreaccIn(Diffusion_Coeff_j, nVar);
  if (correct_gradient) {
    AD::SetPreaccIn(scalar_i, nVar); AD::SetPreaccIn(scalar_j, nVar);
  }
  ExtraADPreaccIn();
  
  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+6); AD::SetPreaccIn(V_j, nDim+6);
    
    Density_i           = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+4];  Laminar_Viscosity_j = V_j[nDim+4];
    Eddy_Viscosity_i    = V_i[nDim+5];     Eddy_Viscosity_j = V_j[nDim+5];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);
    
    Density_i           = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i    = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
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
    Proj_Mean_GradScalarVar_Normal[iVar] = 0.0;
    Proj_Mean_GradScalarVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradScalarVar[iVar][iDim] = 0.5*(scalar_grad_i[iVar][iDim] +
                                            scalar_grad_j[iVar][iDim]);
      Proj_Mean_GradScalarVar_Normal[iVar] += Mean_GradScalarVar[iVar][iDim] *
      Normal[iDim];
      if (correct_gradient)
        Proj_Mean_GradScalarVar_Edge[iVar] += Mean_GradScalarVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradScalarVar[iVar] = Proj_Mean_GradScalarVar_Normal[iVar];
    if (correct_gradient) {
      Proj_Mean_GradScalarVar[iVar] -= Proj_Mean_GradScalarVar_Edge[iVar]*proj_vector_ij -
      (scalar_j[iVar]-scalar_i[iVar])*proj_vector_ij;
    }
  }
  
  FinishResidualCalc(val_residual, Jacobian_i, Jacobian_j, config);
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}

CAvgGradtransportedScalar_General::CAvgGradtransportedScalar_General(unsigned short val_nDim,
                                               unsigned short val_nVar, bool correct_grad,
                                               CConfig *config)
: CAvgGradtransportedScalar(val_nDim, val_nVar, correct_grad, config) {
  
  Mean_Diffusivity = new su2double[nVar];
  
}

CAvgGradtransportedScalar_General::~CAvgGradtransportedScalar_General(void) {
  
  if (Mean_Diffusivity != NULL) delete [] Mean_Diffusivity;
  
}

void CAvgGradtransportedScalar_General::ExtraADPreaccIn() { }

void CAvgGradtransportedScalar_General::FinishResidualCalc(su2double *val_residual,
                                                su2double **Jacobian_i,
                                                su2double **Jacobian_j,
                                                CConfig *config) {
  
  unsigned short iVar, jVar;

  for (iVar = 0; iVar < nVar; iVar++) {
    
    /*--- Get the diffusion coefficient(s). ---*/

    Mean_Diffusivity[iVar] = 0.5*(Diffusion_Coeff_i[iVar] + Diffusion_Coeff_j[iVar]);

    /*--- Compute the viscous residual. ---*/
    
    val_residual[iVar] = Mean_Diffusivity[iVar]*Proj_Mean_GradScalarVar[iVar];
    
    /*--- Use TSL approx. to compute derivatives of the gradients. ---*/

    if (implicit) {
      for (jVar = 0; jVar < nVar; jVar++) {
        if (iVar == jVar) {
          Jacobian_i[iVar][jVar] = -Mean_Diffusivity[iVar]*proj_vector_ij;
          Jacobian_j[iVar][jVar] =  Mean_Diffusivity[iVar]*proj_vector_ij;
        } else {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    
  }
  
}

CSourcePieceWise_Scalar::CSourcePieceWise_Scalar(unsigned short val_nDim,
                                                 unsigned short val_nVar,
                                                 CConfig *config) :
CNumerics(val_nDim, val_nVar, config) {
  
  implicit       = (config->GetKind_TimeIntScheme_Scalar() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
}

CSourcePieceWise_Scalar::~CSourcePieceWise_Scalar(void) { }

void CSourcePieceWise_Scalar::ComputeResidual(su2double *val_residual,
                                              su2double **val_Jacobian_i,
                                              su2double **val_Jacobian_j,
                                              CConfig *config) {
  
  unsigned short iVar, jVar;
  
  Density_i = V_i[nDim+2];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.0;
    if (implicit) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_i[iVar][jVar] = 0.0;
      }
    }
  }
  
}

CSourceAxisymmetric_Scalar::CSourceAxisymmetric_Scalar(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  energy   = config->GetEnergy_Equation();
  viscous  = config->GetViscous();
  
}

CSourceAxisymmetric_Scalar::~CSourceAxisymmetric_Scalar(void) { }

void CSourceAxisymmetric_Scalar::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config) {
  
  su2double yinv, Velocity_i[3];
  unsigned short iDim, iVar, jVar;
  
  if (Coord_i[1] > EPS) {
    yinv          = 1.0/Coord_i[1];
    Density_i     = V_i[nDim+2];
    
    /*--- Set primitive variables at points iPoint. ---*/
    
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity_i[iDim] = V_i[iDim+1];
    
    /*--- Inviscid component of the source term. ---*/
    
    for (iVar=0; iVar < nVar; iVar++)
      val_residual[iVar] = yinv*Volume*Density_i*scalar_i[iVar]*Velocity_i[1];
    
    if (implicit) {
      
      for (iVar=0; iVar < nVar; iVar++) {
        for (jVar=0; jVar < nVar; jVar++) {
          if (iVar == jVar) Jacobian_i[iVar][jVar] = Velocity_i[1];
          Jacobian_i[iVar][jVar] *= yinv*Volume*Density_i;
        }
      }
      
    }
    
    /*--- Add the viscous terms if necessary. ---*/
    
    if (viscous) {
      
      for (iVar=0; iVar < nVar; iVar++)
        val_residual[iVar] -= Volume*yinv*Diffusion_Coeff_i[iVar]*scalar_grad_i[iVar][1];
      
    }
    
  } else {
    
    for (iVar=0; iVar < nVar; iVar++)
      val_residual[iVar] = 0.0;
    
    if (implicit) {
      for (iVar=0; iVar < nVar; iVar++) {
        for (jVar=0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;
      }
    }
    
  }
  
}

