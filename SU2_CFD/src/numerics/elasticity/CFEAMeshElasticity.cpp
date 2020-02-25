/*!
 * \file CFEAMeshElasticity.cpp
 * \brief This file contains the routines for setting the mesh pseudo-elastic problem.
 * \author Ruben Sanchez
 * \version 7.0.1 "Blackbird"
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

#include "../../../include/numerics/elasticity/CFEAMeshElasticity.hpp"


CFEAMeshElasticity::CFEAMeshElasticity(unsigned short val_nDim, unsigned short val_nVar,
                                       unsigned long val_nElem, CConfig *config) :
                                       CFEALinearElasticity() {
  DV_Val         = NULL;
  FAux_Dead_Load = NULL;
  Rho_s_i        = NULL;
  Rho_s_DL_i     = NULL;
  Nu_i           = NULL;

  nDim = val_nDim;
  nVar = val_nVar;

  unsigned long iVar;

  E = config->GetDeform_ElasticityMod();
  Nu = config->GetDeform_PoissonRatio();
  Compute_Lame_Parameters();

  switch (config->GetDeform_Stiffness_Type()) {
  case INVERSE_VOLUME:
  case SOLID_WALL_DISTANCE:
    element_based = true;
    Nu = config->GetDeform_Coeff();
    break;
  case CONSTANT_STIFFNESS:
    element_based = false;
    break;
  }

  E_i  = NULL;
  if (element_based){
    E_i = new su2double[val_nElem];
    for (iVar = 0; iVar < val_nElem; iVar++){
      E_i[iVar] = E;
    }
  }

  KAux_ab = new su2double* [nDim];
  for (iVar = 0; iVar < nDim; iVar++) {
    KAux_ab[iVar] = new su2double[nDim];
  }

  unsigned short nStrain = (nDim==2) ? DIM_STRAIN_2D : DIM_STRAIN_3D;
  unsigned short nNodes = (nDim==2) ? NNODES_2D : NNODES_3D;

  Ba_Mat = new su2double* [nStrain];
  Bb_Mat = new su2double* [nStrain];
  D_Mat  = new su2double* [nStrain];
  Ni_Vec  = new su2double [nNodes];
  GradNi_Ref_Mat = new su2double* [nNodes];
  GradNi_Curr_Mat = new su2double* [nNodes];
  for (iVar = 0; iVar < nStrain; iVar++) {
    Ba_Mat[iVar] = new su2double[nDim];
    Bb_Mat[iVar] = new su2double[nDim];
    D_Mat[iVar] = new su2double[nStrain];
  }
  for (iVar = 0; iVar < nNodes; iVar++) {
    GradNi_Ref_Mat[iVar] = new su2double[nDim];
    GradNi_Curr_Mat[iVar] = new su2double[nDim];
  }

}
