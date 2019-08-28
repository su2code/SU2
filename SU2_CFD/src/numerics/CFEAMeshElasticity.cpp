/*!
 * \file CFEAMeshElasticity.cpp
 * \brief This file contains the routines for setting the mesh pseudo-elastic problem.
 * \author Ruben Sanchez
 * \version 6.2.0 "Falcon"
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

#include "../../include/numerics/CFEAMeshElasticity.hpp"
#include <limits>

CFEAMeshElasticity::CFEAMeshElasticity(unsigned short val_nDim, unsigned short val_nVar, unsigned long val_nElem, CConfig *config) : CFEALinearElasticity() {

  DV_Val         = NULL;
  FAux_Dead_Load = NULL;
  Rho_s_i        = NULL;
  Rho_s_DL_i     = NULL;
  Nu_i           = NULL;

  nDim = val_nDim;
  nVar = val_nVar;

  const unsigned short DIM_STRAIN_2D = 3; //Exx, Eyy, Gxy
  const unsigned short DIM_STRAIN_3D = 6; //Exx, Eyy, Ezz, Gxy, Gxz, Gyz

  const unsigned short NNODES_2D = 4;     // Maximum number of nodes for 2D problems
  const unsigned short NNODES_3D = 8;     // Maximum number of nodes for 3D problems

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
    E_i         = new su2double[val_nElem];
    for (iVar = 0; iVar < val_nElem; iVar++){
      E_i[iVar] = E;
    }
  }

  KAux_ab = new su2double* [nDim];
  for (iVar = 0; iVar < nDim; iVar++) {
    KAux_ab[iVar] = new su2double[nDim];
  }

  if (nDim == 2) {
    Ba_Mat = new su2double* [DIM_STRAIN_2D];
    Bb_Mat = new su2double* [DIM_STRAIN_2D];
    D_Mat  = new su2double* [DIM_STRAIN_2D];
    Ni_Vec  = new su2double [NNODES_2D];
    GradNi_Ref_Mat = new su2double* [NNODES_2D];
    GradNi_Curr_Mat = new su2double* [NNODES_2D];
    for (iVar = 0; iVar < DIM_STRAIN_2D; iVar++) {
      Ba_Mat[iVar] = new su2double[nDim];
      Bb_Mat[iVar] = new su2double[nDim];
      D_Mat[iVar]  = new su2double[DIM_STRAIN_2D];
    }
    for (iVar = 0; iVar < NNODES_2D; iVar++) {
      GradNi_Ref_Mat[iVar]   = new su2double[nDim];
      GradNi_Curr_Mat[iVar]   = new su2double[nDim];
    }
  }
  else if (nDim == 3) {
    Ba_Mat = new su2double* [DIM_STRAIN_3D];
    Bb_Mat = new su2double* [DIM_STRAIN_3D];
    D_Mat  = new su2double* [DIM_STRAIN_3D];
    Ni_Vec  = new su2double [NNODES_3D];
    GradNi_Ref_Mat = new su2double* [NNODES_3D];
    GradNi_Curr_Mat = new su2double* [NNODES_3D];
    for (iVar = 0; iVar < DIM_STRAIN_3D; iVar++) {
      Ba_Mat[iVar] = new su2double[nDim];
      Bb_Mat[iVar] = new su2double[nDim];
      D_Mat[iVar]  = new su2double[DIM_STRAIN_3D];
    }
    for (iVar = 0; iVar < NNODES_3D; iVar++) {
      GradNi_Ref_Mat[iVar]   = new su2double[nDim];
      GradNi_Curr_Mat[iVar]   = new su2double[nDim];
    }
  }

}

CFEAMeshElasticity::~CFEAMeshElasticity(void) {

}

