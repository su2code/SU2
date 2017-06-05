/*!
 * \file numerics_direct_hybrid.cpp
 * \brief This file contains all the discretizations for hybrid parameter(s).
 * \author C. Pederson
 * \version 5.0.0 "Raven"
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
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

CUpwSca_HybridConv::CUpwSca_HybridConv(unsigned short val_nDim,
                                           unsigned short val_nVar,
                                           CConfig *config)
: CNumerics(val_nDim, val_nVar, config) {

  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();

  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];

}

CUpwSca_HybridConv::~CUpwSca_HybridConv(void) {

  delete [] Velocity_i;
  delete [] Velocity_j;

}

void CUpwSca_HybridConv::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  q_ij = 0.0;

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+1); AD::SetPreaccIn(V_j, nDim+1);
  AD::SetPreaccIn(HybridParameter_i[0]); AD::SetPreaccIn(HybridParameter_j[0]);
  AD::SetPreaccIn(Normal, nDim);
  if (grid_movement) {
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
  val_residual[0] = a0*HybridParameter_i[0]+a1*HybridParameter_j[0];

  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_j[0][0] = a1;
  }

  AD::SetPreaccOut(val_residual[0]);
  AD::EndPreacc();
}

CSourcePieceWise_HybridConv::CSourcePieceWise_HybridConv(unsigned short val_nDim, unsigned short val_nVar,
                                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  rotating_frame = config->GetRotating_Frame();
  transition = (config->GetKind_Trans_Model() == BC);

}

CSourcePieceWise_HybridConv::~CSourcePieceWise_HybridConv(void) { }

void CSourcePieceWise_HybridConv::ComputeResidual(su2double *val_residual,
                                                    su2double **val_Jacobian_i,
                                                    su2double **val_Jacobian_j,
                                                    CConfig *config) {
  su2double udMdx; // Dot product of resolved velocity and resolution tensor gradient
  su2double max_udMdx; // Max value of the convection of the resolution tensor
  su2double S_r, // Gentle switch between raising and lowering hybrid parameter
            S_c, // Source for removing turbulent scales
            T_c; // Timescale for removal of turbulent scales

  if (Resolution_Adequacy >=1.0)
    S_r = tanh(Resolution_Adequacy - 1.0);
  else
    S_r = tanh(1.0 - 1.0/Resolution_Adequacy);


  for (unsigned int iDim = 0; iDim < nDim; iDim++) {
    for (unsigned int jDim = 0; jDim < nDim; jDim++) {
      udMdx = 0.0;
      for (unsigned int kDim = 0; kDim < nDim; nDim++) {
        // XXX: This is NOT upwinded.  It could be improved to upwind the product.
        udMdx += V_i[kDim+1] * Resolution_Tensor_Gradient[kDim][iDim][jDim];
      }
      if (iDim == 0 and jDim == 0) max_udMdx = udMdx;
      else if (udMdx > max_udMdx) max_udMdx = udMdx;
    }
  }

  T_c = TurbL/max_udMdx;

  if (S_r >= 0.0 && T_c >= 0.0)
    S_c = 1.0;
  else
    S_c = 0.0;

  val_residual[0] = (S_r/TurbT + S_c/T_c) * Volume;

  // Jacobian of \alpha with respect to \alpha
  val_Jacobian_i[0][0] = 0.0; val_Jacobian_j[0][0] = 0.0;
}

