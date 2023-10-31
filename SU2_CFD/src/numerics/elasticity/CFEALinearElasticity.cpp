/*!
 * \file CFEALinearElasticity.cpp
 * \brief Classes for linear elasticity problems.
 * \author R. Sanchez
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

#include "../../../include/numerics/elasticity/CFEALinearElasticity.hpp"


CFEALinearElasticity::CFEALinearElasticity(unsigned short val_nDim, unsigned short val_nVar,
                                           const CConfig *config) : CFEAElasticity(val_nDim, val_nVar, config) {
  if (nDim == 2)
    nodalDisplacement.resize(NNODES_2D,nDim);
  else
    nodalDisplacement.resize(NNODES_3D,nDim);
}

void CFEALinearElasticity::Compute_Tangent_Matrix(CElement *element, const CConfig *config) {

  unsigned short iVar, jVar, kVar;
  unsigned short iGauss, nGauss;
  unsigned short iNode, jNode, nNode;
  unsigned short iDim;
  unsigned short bDim;

  su2double Weight, Jac_X;

  su2double AuxMatrix[3][6], *res_aux = new su2double[nVar];

  /*--- Set element properties and recompute the constitutive matrix, this is needed
        for multiple material cases and for correct differentiation ---*/
  SetElement_Properties(element, config);

  /*--- Register pre-accumulation inputs, material props and nodal coords ---*/
  AD::StartPreacc();
  AD::SetPreaccIn(E);
  AD::SetPreaccIn(Nu);
  element->SetPreaccIn_Coords();
  /*--- Recompute Lame parameters as they depend on the material properties ---*/
  Compute_Lame_Parameters();

  Compute_Constitutive_Matrix(element, config);

  /*--- Initialize auxiliary matrices ---*/

  bDim = (nDim == 2) ? DIM_STRAIN_2D : DIM_STRAIN_3D;

  for (iVar = 0; iVar < bDim; iVar++) {
    for (jVar = 0; jVar < nDim; jVar++) {
      Ba_Mat[iVar][jVar] = 0.0;
      Bb_Mat[iVar][jVar] = 0.0;
    }
  }

  for (iVar = 0; iVar < 3; iVar++) {
    for (jVar = 0; jVar < 6; jVar++) {
      AuxMatrix[iVar][jVar] = 0.0;
    }
  }

  element->ClearElement();       /*--- Restarts the element: avoids adding over previous results in other elements --*/
  element->ComputeGrad_Linear();
  nNode = element->GetnNodes();
  nGauss = element->GetnGaussPoints();

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    Weight = element->GetWeight(iGauss);
    Jac_X = element->GetJ_X(iGauss);

    /*--- Retrieve the values of the gradients of the shape functions for each node ---*/
    /*--- This avoids repeated operations ---*/
    for (iNode = 0; iNode < nNode; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
      }
    }

    for (iNode = 0; iNode < nNode; iNode++) {

      if (nDim == 2) {
        Ba_Mat[0][0] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[1][1] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[2][0] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[2][1] = GradNi_Ref_Mat[iNode][0];
      }
      else {
        Ba_Mat[0][0] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[1][1] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[2][2] = GradNi_Ref_Mat[iNode][2];
        Ba_Mat[3][0] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[3][1] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[4][0] = GradNi_Ref_Mat[iNode][2];
        Ba_Mat[4][2] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[5][1] = GradNi_Ref_Mat[iNode][2];
        Ba_Mat[5][2] = GradNi_Ref_Mat[iNode][1];
      }

      /*--- Compute the BT.D Matrix ---*/

      for (iVar = 0; iVar < nDim; iVar++) {
        for (jVar = 0; jVar < bDim; jVar++) {
          AuxMatrix[iVar][jVar] = 0.0;
          for (kVar = 0; kVar < bDim; kVar++) {
            AuxMatrix[iVar][jVar] += Ba_Mat[kVar][iVar]*D_Mat[kVar][jVar];
          }
        }
      }

      /*--- Assumming symmetry ---*/
      for (jNode = iNode; jNode < nNode; jNode++) {
        if (nDim == 2) {
          Bb_Mat[0][0] = GradNi_Ref_Mat[jNode][0];
          Bb_Mat[1][1] = GradNi_Ref_Mat[jNode][1];
          Bb_Mat[2][0] = GradNi_Ref_Mat[jNode][1];
          Bb_Mat[2][1] = GradNi_Ref_Mat[jNode][0];
        }
        else {
          Bb_Mat[0][0] = GradNi_Ref_Mat[jNode][0];
          Bb_Mat[1][1] = GradNi_Ref_Mat[jNode][1];
          Bb_Mat[2][2] = GradNi_Ref_Mat[jNode][2];
          Bb_Mat[3][0] = GradNi_Ref_Mat[jNode][1];
          Bb_Mat[3][1] = GradNi_Ref_Mat[jNode][0];
          Bb_Mat[4][0] = GradNi_Ref_Mat[jNode][2];
          Bb_Mat[4][2] = GradNi_Ref_Mat[jNode][0];
          Bb_Mat[5][1] = GradNi_Ref_Mat[jNode][2];
          Bb_Mat[5][2] = GradNi_Ref_Mat[jNode][1];
        }

        for (iVar = 0; iVar < nDim; iVar++) {
          for (jVar = 0; jVar < nDim; jVar++) {
            KAux_ab[iVar][jVar] = 0.0;
            for (kVar = 0; kVar < bDim; kVar++) {
              KAux_ab[iVar][jVar] += Weight * AuxMatrix[iVar][kVar] * Bb_Mat[kVar][jVar] * Jac_X;
            }
          }
        }

        element->Add_Kab(iNode, jNode, KAux_ab);
        /*--- Symmetric terms --*/
        if (iNode != jNode) {
          element->Add_Kab_T(jNode, iNode, KAux_ab);
        }

      }

    }

  }

  /*--- Compute residual ---*/
  for(iNode = 0; iNode<nNode; ++iNode)
  {
    for(jNode = 0; jNode<nNode; ++jNode)
    {
      const su2double *Kab = element->Get_Kab(iNode,jNode);

      for (iVar = 0; iVar < nVar; iVar++) {
        res_aux[iVar] = 0.0;
        for (jVar = 0; jVar < nVar; jVar++)
          res_aux[iVar] += Kab[iVar*nVar+jVar]*
            (element->GetCurr_Coord(jNode,jVar)-element->GetRef_Coord(jNode,jVar));
      }
      element->Add_Kt_a(iNode, res_aux);
    }
  }

  /*--- Register the stress residual as preaccumulation output ---*/
  element->SetPreaccOut_Kt_a();
  AD::EndPreacc();

  delete[] res_aux;
}


void CFEALinearElasticity::Compute_Constitutive_Matrix(CElement *element_container, const CConfig *config) {

  /*--- Compute the D Matrix (for plane stress and 2-D)---*/


  if (nDim == 2) {
    if (plane_stress) {
      /*--- We enable plane stress cases ---*/

      D_Mat[0][0] = E/(1-Nu*Nu);        D_Mat[0][1] = (E*Nu)/(1-Nu*Nu);  D_Mat[0][2] = 0.0;
      D_Mat[1][0] = (E*Nu)/(1-Nu*Nu);   D_Mat[1][1] = E/(1-Nu*Nu);       D_Mat[1][2] = 0.0;
      D_Mat[2][0] = 0.0;                D_Mat[2][1] = 0.0;               D_Mat[2][2] = ((1-Nu)*E)/(2*(1-Nu*Nu));
    }
    else {
      /*--- Assuming plane strain as a general case ---*/

      D_Mat[0][0] = Lambda + 2.0*Mu;  D_Mat[0][1] = Lambda;           D_Mat[0][2] = 0.0;
      D_Mat[1][0] = Lambda;           D_Mat[1][1] = Lambda + 2.0*Mu;  D_Mat[1][2] = 0.0;
      D_Mat[2][0] = 0.0;              D_Mat[2][1] = 0.0;              D_Mat[2][2] = Mu;
    }

  }
  else {

    su2double Lbda_2Mu = Lambda + 2.0*Mu;

    D_Mat[0][0] = Lbda_2Mu;  D_Mat[0][1] = Lambda;    D_Mat[0][2] = Lambda;    D_Mat[0][3] = 0.0;  D_Mat[0][4] = 0.0;  D_Mat[0][5] = 0.0;
    D_Mat[1][0] = Lambda;    D_Mat[1][1] = Lbda_2Mu;  D_Mat[1][2] = Lambda;    D_Mat[1][3] = 0.0;  D_Mat[1][4] = 0.0;  D_Mat[1][5] = 0.0;
    D_Mat[2][0] = Lambda;    D_Mat[2][1] = Lambda;    D_Mat[2][2] = Lbda_2Mu;  D_Mat[2][3] = 0.0;  D_Mat[2][4] = 0.0;  D_Mat[2][5] = 0.0;
    D_Mat[3][0] = 0.0;       D_Mat[3][1] = 0.0;       D_Mat[3][2] = 0.0;       D_Mat[3][3] = Mu;   D_Mat[3][4] = 0.0;  D_Mat[3][5] = 0.0;
    D_Mat[4][0] = 0.0;       D_Mat[4][1] = 0.0;       D_Mat[4][2] = 0.0;       D_Mat[4][3] = 0.0;  D_Mat[4][4] = Mu;   D_Mat[4][5] = 0.0;
    D_Mat[5][0] = 0.0;       D_Mat[5][1] = 0.0;       D_Mat[5][2] = 0.0;       D_Mat[5][3] = 0.0;  D_Mat[5][4] = 0.0;  D_Mat[5][5] = Mu;

  }

}


su2double CFEALinearElasticity::Compute_Averaged_NodalStress(CElement *element, const CConfig *config) {

  unsigned short iVar, jVar;
  unsigned short iGauss, nGauss;
  unsigned short iNode, nNode;
  unsigned short iDim, bDim;

  su2double avgStress[DIM_STRAIN_3D] = {0.0};

  /*--- Set element properties and recompute the constitutive matrix, this is needed
        for multiple material cases and for correct differentiation ---*/
  SetElement_Properties(element, config);

  /*--- Register pre-accumulation inputs ---*/
  AD::StartPreacc();
  AD::SetPreaccIn(E);
  AD::SetPreaccIn(Nu);
  element->SetPreaccIn_Coords();
  /*--- Recompute Lame parameters as they depend on the material properties ---*/
  Compute_Lame_Parameters();

  Compute_Constitutive_Matrix(element, config);

  /*--- Initialize auxiliary matrices ---*/

  bDim = (nDim == 2) ? DIM_STRAIN_2D : DIM_STRAIN_3D;

  for (iVar = 0; iVar < bDim; iVar++) {
    for (jVar = 0; jVar < nDim; jVar++) {
      Ba_Mat[iVar][jVar] = 0.0;
    }
  }

  element->ClearStress(); /*--- Clears the stress in the element to avoid adding over previous results. --*/
  element->ComputeGrad_Linear();
  nNode = element->GetnNodes();
  nGauss = element->GetnGaussPoints();

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    /*--- Retrieve the values of the gradients of the shape functions for each node ---*/
    /*--- This avoids repeated operations ---*/
    for (iNode = 0; iNode < nNode; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
        nodalDisplacement[iNode][iDim] = element->GetCurr_Coord(iNode, iDim) - element->GetRef_Coord(iNode, iDim);
      }
    }

    su2double Strain[DIM_STRAIN_3D] = {0.0};

    for (iNode = 0; iNode < nNode; iNode++) {

      /*--- Set matrix B ---*/
      if (nDim == 2) {
        Ba_Mat[0][0] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[1][1] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[2][0] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[2][1] = GradNi_Ref_Mat[iNode][0];
      }
      else {
        Ba_Mat[0][0] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[1][1] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[2][2] = GradNi_Ref_Mat[iNode][2];
        Ba_Mat[3][0] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[3][1] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[4][0] = GradNi_Ref_Mat[iNode][2];
        Ba_Mat[4][2] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[5][1] = GradNi_Ref_Mat[iNode][2];
        Ba_Mat[5][2] = GradNi_Ref_Mat[iNode][1];
      }

      /*--- Compute the Strain Vector as B*u ---*/

      for (iVar = 0; iVar < bDim; iVar++) {
        for (jVar = 0; jVar < nDim; jVar++) {
          Strain[iVar] += Ba_Mat[iVar][jVar]*nodalDisplacement[iNode][jVar];
        }
      }

    }

    /*--- Compute the Stress Vector as D*epsilon ---*/

    su2double Stress[DIM_STRAIN_3D] = {0.0};

    for (iVar = 0; iVar < bDim; iVar++) {
      for (jVar = 0; jVar < bDim; jVar++) {
        Stress[iVar] += D_Mat[iVar][jVar]*Strain[jVar];
      }
      avgStress[iVar] += Stress[iVar] / nGauss;
    }

    for (iNode = 0; iNode < nNode; iNode++) {

      su2double Ni_Extrap = element->GetNi_Extrap(iNode, iGauss);

      if (nDim == 2) {
        for(iVar = 0; iVar < 3; ++iVar)
          element->Add_NodalStress(iNode, iVar, Stress[iVar] * Ni_Extrap);
      }
      else {
        /*--- If nDim is 3 and we compute it this way, the 3rd component is the Szz,
         *    while in the output it is the 4th component for practical reasons. ---*/
        element->Add_NodalStress(iNode, 0, Stress[0] * Ni_Extrap);
        element->Add_NodalStress(iNode, 1, Stress[1] * Ni_Extrap);
        element->Add_NodalStress(iNode, 2, Stress[3] * Ni_Extrap);
        element->Add_NodalStress(iNode, 3, Stress[2] * Ni_Extrap);
        element->Add_NodalStress(iNode, 4, Stress[4] * Ni_Extrap);
        element->Add_NodalStress(iNode, 5, Stress[5] * Ni_Extrap);
      }
    }

  }

  if (nDim == 3) std::swap(avgStress[2], avgStress[3]);
  auto elStress = VonMisesStress(nDim, avgStress);

  /*--- We only differentiate w.r.t. an avg VM stress for the element as
   * considering all nodal stresses would use too much memory. ---*/
  AD::SetPreaccOut(elStress);
  AD::EndPreacc();

  return elStress;
}


CFEAMeshElasticity::CFEAMeshElasticity(unsigned short val_nDim, unsigned short val_nVar,
                                       unsigned long val_nElem, const CConfig *config) :
                                       CFEALinearElasticity() {
  DV_Val         = nullptr;
  FAux_Dead_Load = nullptr;
  Rho_s_i        = nullptr;
  Rho_s_DL_i     = nullptr;
  Nu_i           = nullptr;

  nDim = val_nDim;
  nVar = val_nVar;

  unsigned long iVar;

  E = 1.0;
  Nu = config->GetDeform_PoissonRatio();
  Compute_Lame_Parameters();

  switch (config->GetDeform_Stiffness_Type()) {
  case INVERSE_VOLUME:
  case SOLID_WALL_DISTANCE:
    element_based = true;
    break;
  case CONSTANT_STIFFNESS:
    element_based = false;
    break;
  }

  E_i  = nullptr;
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
