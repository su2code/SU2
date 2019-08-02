/*!
 * \file numerics_direct_elasticity_linear.cpp
 * \brief This file contains the routines for setting the FEM elastic structural problem.
 * \author R. Sanchez
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

#include "../include/numerics_structure.hpp"
#include <limits>

CFEALinearElasticity::CFEALinearElasticity(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEAElasticity(val_nDim, val_nVar, config) {

  unsigned short iVar;

  if (nDim == 2) {
    nodalDisplacement = new su2double* [4];  /*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
    for (iVar = 0; iVar < 4; iVar++) nodalDisplacement[iVar] = new su2double[nDim];
  }
  else if (nDim == 3) {
    nodalDisplacement = new su2double* [8];  /*--- As of now, 8 is the maximum number of nodes for 3D problems ---*/
    for (iVar = 0; iVar < 8; iVar++) nodalDisplacement[iVar] = new su2double[nDim];
  }

}

CFEALinearElasticity::~CFEALinearElasticity(void) {

}

void CFEALinearElasticity::Compute_Tangent_Matrix(CElement *element, CConfig *config) {

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
  AD::SetPreaccIn(Rho_s);
  AD::SetPreaccIn(Rho_s_DL);
  element->SetPreaccIn_Coords();
  /*--- Recompute Lame parameters as they depend on the material properties ---*/
  Compute_Lame_Parameters();
  
  Compute_Constitutive_Matrix(element, config);

  /*--- Initialize auxiliary matrices ---*/

  if (nDim == 2) bDim = 3;
  else bDim = 6;

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

  element->clearElement();       /*--- Restarts the element: avoids adding over previous results in other elements --*/
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
      else if (nDim == 3) {
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
        else if (nDim ==3) {
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

        element->Add_Kab(KAux_ab,iNode, jNode);
        /*--- Symmetric terms --*/
        if (iNode != jNode) {
          element->Add_Kab_T(KAux_ab, jNode, iNode);
        }

      }

    }

  }
  
  // compute residual
  for(iNode = 0; iNode<nNode; ++iNode)
  {
    for(jNode = 0; jNode<nNode; ++jNode)
    {
      su2double *Kab = element->Get_Kab(iNode,jNode);
      
      for (iVar = 0; iVar < nVar; iVar++) {
          res_aux[iVar] = 0.0;
          for (jVar = 0; jVar < nVar; jVar++)
            res_aux[iVar] += Kab[iVar*nVar+jVar]*
              (element->GetCurr_Coord(jNode,jVar)-element->GetRef_Coord(jNode,jVar));
      }
      element->Add_Kt_a(res_aux,iNode);
    }
  }
  
  /*--- Register the stress residual as preaccumulation output ---*/
  element->SetPreaccOut_Kt_a();
  AD::EndPreacc();
  
  delete[] res_aux;
}


void CFEALinearElasticity::Compute_Constitutive_Matrix(CElement *element_container, CConfig *config) {

     /*--- Compute the D Matrix (for plane stress and 2-D)---*/


  if (nDim == 2) {
    if (plane_stress) {

      /*--- We enable plane stress cases ---*/

      D_Mat[0][0] = E/(1-Nu*Nu);        D_Mat[0][1] = (E*Nu)/(1-Nu*Nu);  D_Mat[0][2] = 0.0;
      D_Mat[1][0] = (E*Nu)/(1-Nu*Nu);      D_Mat[1][1] = E/(1-Nu*Nu);      D_Mat[1][2] = 0.0;
      D_Mat[2][0] = 0.0;                 D_Mat[2][1] = 0.0;               D_Mat[2][2] = ((1-Nu)*E)/(2*(1-Nu*Nu));
    }
    else {
      /*--- Assuming plane strain as a general case ---*/

      D_Mat[0][0] = Lambda + 2.0*Mu;  D_Mat[0][1] = Lambda;            D_Mat[0][2] = 0.0;
      D_Mat[1][0] = Lambda;           D_Mat[1][1] = Lambda + 2.0*Mu;   D_Mat[1][2] = 0.0;
      D_Mat[2][0] = 0.0;              D_Mat[2][1] = 0.0;               D_Mat[2][2] = Mu;
    }

  }
  else if (nDim == 3) {

    D_Mat[0][0] = Lambda + 2.0*Mu;  D_Mat[0][1] = Lambda;      D_Mat[0][2] = Lambda;      D_Mat[0][3] = 0.0;  D_Mat[0][4] = 0.0;  D_Mat[0][5] = 0.0;
    D_Mat[1][0] = Lambda;      D_Mat[1][1] = Lambda + 2.0*Mu;  D_Mat[1][2] = Lambda;      D_Mat[1][3] = 0.0;  D_Mat[1][4] = 0.0;  D_Mat[1][5] = 0.0;
    D_Mat[2][0] = Lambda;      D_Mat[2][1] = Lambda;      D_Mat[2][2] = Lambda + 2.0*Mu;  D_Mat[2][3] = 0.0;  D_Mat[2][4] = 0.0;  D_Mat[2][5] = 0.0;
    D_Mat[3][0] = 0.0;        D_Mat[3][1] = 0.0;        D_Mat[3][2] = 0.0;        D_Mat[3][3] = Mu;  D_Mat[3][4] = 0.0;  D_Mat[3][5] = 0.0;
    D_Mat[4][0] = 0.0;        D_Mat[4][1] = 0.0;        D_Mat[4][2] = 0.0;        D_Mat[4][3] = 0.0;  D_Mat[4][4] = Mu;  D_Mat[4][5] = 0.0;
    D_Mat[5][0] = 0.0;        D_Mat[5][1] = 0.0;        D_Mat[5][2] = 0.0;        D_Mat[5][3] = 0.0;  D_Mat[5][4] = 0.0;  D_Mat[5][5] = Mu;

  }

}

void CFEALinearElasticity::Compute_Averaged_NodalStress(CElement *element, CConfig *config) {

  unsigned short iVar, jVar;
  unsigned short iGauss, nGauss;
  unsigned short iNode, nNode;
  unsigned short iDim, bDim;

  /*--- Auxiliary vector ---*/
  su2double Strain[6], Stress[6];
  
  /*--- Set element properties and recompute the constitutive matrix, this is needed
        for multiple material cases and for correct differentiation ---*/
  SetElement_Properties(element, config);
  Compute_Constitutive_Matrix(element, config);

  /*--- Initialize auxiliary matrices ---*/

  if (nDim == 2) bDim = 3;
  else bDim = 6;

  for (iVar = 0; iVar < bDim; iVar++) {
    for (jVar = 0; jVar < nDim; jVar++) {
      Ba_Mat[iVar][jVar] = 0.0;
    }
  }

  element->clearStress();       /*--- Clears the stress in the element: avoids adding over previous results in other elements --*/
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

    for (iVar = 0; iVar < bDim; iVar++) {
      Strain[iVar] = 0.0;
    }

    for (iNode = 0; iNode < nNode; iNode++) {

      /*--- Set matrix B ---*/
      if (nDim == 2) {
        Ba_Mat[0][0] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[1][1] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[2][0] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[2][1] = GradNi_Ref_Mat[iNode][0];
      }
      else if (nDim ==3) {
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

    for (iVar = 0; iVar < bDim; iVar++) {
      Stress[iVar] = 0.0;
      for (jVar = 0; jVar < bDim; jVar++) {
        Stress[iVar] += D_Mat[iVar][jVar]*Strain[jVar];
      }
    }

    for (iNode = 0; iNode < nNode; iNode++) {
      /*--- If nDim is 3 and we compute it this way, the 3rd component is the Szz, while in the ---*/
      /*--- output it is the 4th component for practical reasons ---*/
      if (nDim == 2) {
        element->Add_NodalStress(Stress[0] * element->GetNi_Extrap(iNode, iGauss), iNode, 0);
        element->Add_NodalStress(Stress[1] * element->GetNi_Extrap(iNode, iGauss), iNode, 1);
        element->Add_NodalStress(Stress[2] * element->GetNi_Extrap(iNode, iGauss), iNode, 2);
      }
      else if (nDim == 3) {
        element->Add_NodalStress(Stress[0] * element->GetNi_Extrap(iNode, iGauss), iNode, 0);
        element->Add_NodalStress(Stress[1] * element->GetNi_Extrap(iNode, iGauss), iNode, 1);
        element->Add_NodalStress(Stress[3] * element->GetNi_Extrap(iNode, iGauss), iNode, 2);
        element->Add_NodalStress(Stress[2] * element->GetNi_Extrap(iNode, iGauss), iNode, 3);
        element->Add_NodalStress(Stress[4] * element->GetNi_Extrap(iNode, iGauss), iNode, 4);
        element->Add_NodalStress(Stress[5] * element->GetNi_Extrap(iNode, iGauss), iNode, 5);
      }
    }



  }


}
