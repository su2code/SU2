/*!
 * \file CFEALinearElasticity.cpp
 * \brief Classes for linear elasticity problems.
 * \author R. Sanchez
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
  unsigned short iGauss;
  unsigned short iNode, jNode;
  unsigned short iDim;

  su2double AuxMatrix[MAXNDIM][DIM_STRAIN_3D] = {};
  su2double Ba_Mat[DIM_STRAIN_3D][MAXNDIM] = {}, Bb_Mat[DIM_STRAIN_3D][MAXNDIM] = {};

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

  /*--- Restart the element: avoids adding over previous results in other elements --*/
  element->ClearElement();
  element->ComputeGrad_Linear();
  const auto nNode = element->GetnNodes();
  const auto nGauss = element->GetnGaussPoints();
  const auto bDim = (nDim == 2) ? DIM_STRAIN_2D : DIM_STRAIN_3D;

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    const su2double Weight = element->GetWeight(iGauss);
    const su2double Jac_X = element->GetJ_X(iGauss);

    /*--- Retrieve the values of the gradients of the shape functions for each node ---*/
    /*--- This avoids repeated operations ---*/
    su2double GradNi_Ref_Mat[NNODES_3D][MAXNDIM] = {};
    for (iNode = 0; iNode < nNode; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
      }
    }

    const su2double thermalStress = ThermalStressTerm * (element->GetTemperature(iGauss) - ReferenceTemperature);

    for (iNode = 0; iNode < nNode; iNode++) {

      su2double KAux_t_a[MAXNDIM] = {};
      for (iVar = 0; iVar < nDim; iVar++) {
        KAux_t_a[iVar] += Weight * thermalStress * GradNi_Ref_Mat[iNode][iVar] * Jac_X;
      }
      element->Add_Kt_a(iNode, KAux_t_a);

      FillBMat(iNode, GradNi_Ref_Mat, Ba_Mat);

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
        FillBMat(jNode, GradNi_Ref_Mat, Bb_Mat);

        su2double KAux_ab[MAXNDIM][MAXNDIM] = {};
        for (iVar = 0; iVar < nDim; iVar++) {
          for (jVar = 0; jVar < nDim; jVar++) {
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
  for (iNode = 0; iNode<nNode; ++iNode) {
    for (jNode = 0; jNode<nNode; ++jNode) {
      const su2double *Kab = element->Get_Kab(iNode,jNode);

      su2double res_aux[MAXNDIM] = {};
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          res_aux[iVar] += Kab[iVar*nVar+jVar] *
              (element->GetCurr_Coord(jNode,jVar) - element->GetRef_Coord(jNode,jVar));
        }
      }
      element->Add_Kt_a(iNode, res_aux);
    }
  }

  /*--- Register the stress residual as preaccumulation output ---*/
  element->SetPreaccOut_Kt_a();
  AD::EndPreacc();
}


void CFEALinearElasticity::Compute_Constitutive_Matrix(CElement *element_container, const CConfig *config) {

  /*--- Compute the D Matrix (for plane stress and 2-D)---*/

  if (nDim == 2) {
    if (plane_stress) {
      su2double D = E / (1 - Nu * Nu);
      D_Mat[0][0] = D;        D_Mat[0][1] = Nu * D;  D_Mat[0][2] = 0.0;
      D_Mat[1][0] = Nu * D;   D_Mat[1][1] = D;       D_Mat[1][2] = 0.0;
      D_Mat[2][0] = 0.0;      D_Mat[2][1] = 0.0;     D_Mat[2][2] = (1 - Nu) * D / 2;
    } else {
      D_Mat[0][0] = Lambda + 2.0*Mu;  D_Mat[0][1] = Lambda;           D_Mat[0][2] = 0.0;
      D_Mat[1][0] = Lambda;           D_Mat[1][1] = Lambda + 2.0*Mu;  D_Mat[1][2] = 0.0;
      D_Mat[2][0] = 0.0;              D_Mat[2][1] = 0.0;              D_Mat[2][2] = Mu;
    }
  } else {
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
  unsigned short iGauss;
  unsigned short iNode;
  unsigned short iDim;

  su2double avgStress[DIM_STRAIN_3D] = {};
  su2double Ba_Mat[DIM_STRAIN_3D][MAXNDIM] = {};

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

  /*--- Clears the stress in the element to avoid adding over previous results. --*/
  element->ClearStress();
  element->ComputeGrad_Linear();

  const auto nNode = element->GetnNodes();
  const auto nGauss = element->GetnGaussPoints();
  const auto bDim = (nDim == 2) ? DIM_STRAIN_2D : DIM_STRAIN_3D;

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    /*--- Retrieve the values of the gradients of the shape functions for each node ---*/
    /*--- This avoids repeated operations ---*/
    su2double GradNi_Ref_Mat[NNODES_3D][MAXNDIM] = {};
    for (iNode = 0; iNode < nNode; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
        nodalDisplacement[iNode][iDim] = element->GetCurr_Coord(iNode, iDim) - element->GetRef_Coord(iNode, iDim);
      }
    }

    /*--- Compute the Strain Vector as B*u ---*/
    su2double Strain[DIM_STRAIN_3D] = {};

    for (iNode = 0; iNode < nNode; iNode++) {
      FillBMat(iNode, GradNi_Ref_Mat, Ba_Mat);

      for (iVar = 0; iVar < bDim; iVar++) {
        for (jVar = 0; jVar < nDim; jVar++) {
          Strain[iVar] += Ba_Mat[iVar][jVar]*nodalDisplacement[iNode][jVar];
        }
      }
    }

    /*--- Compute the Stress Vector as D*epsilon + thermal stress ---*/

    su2double Stress[DIM_STRAIN_3D] = {0.0};

    const su2double thermalStress = ThermalStressTerm * (element->GetTemperature(iGauss) - ReferenceTemperature);

    for (iVar = 0; iVar < bDim; iVar++) {
      for (jVar = 0; jVar < bDim; jVar++) {
        Stress[iVar] += D_Mat[iVar][jVar]*Strain[jVar];
      }
      if (iVar < nDim) Stress[iVar] += thermalStress;
      avgStress[iVar] += Stress[iVar] / nGauss;
    }

    if (nDim == 2 && config->GetElas2D_Formulation() == STRUCT_2DFORM::PLANE_STRAIN) {
      Stress[3] = Nu * (Stress[0] + Stress[1]) + (1 - 2 * Nu) * thermalStress;
      avgStress[3] += Stress[3] / nGauss;
    }

    /*--- If nDim is 3 and we compute it this way, the 3rd component is the Szz,
     *    while in the output it is the 4th component for practical reasons. ---*/
    if (nDim == 3) std::swap(Stress[2], Stress[3]);

    for (iNode = 0; iNode < nNode; iNode++) {
      su2double Ni_Extrap = element->GetNi_Extrap(iNode, iGauss);
      for (iVar = 0; iVar < DIM_STRAIN_3D; ++iVar) {
        element->Add_NodalStress(iNode, iVar, Stress[iVar] * Ni_Extrap);
      }
    }

  }

  /*--- See note regarding output order above. ---*/
  if (nDim == 3) std::swap(avgStress[2], avgStress[3]);
  auto elStress = CFEAElasticity::VonMisesStress(nDim, avgStress);

  /*--- We only differentiate w.r.t. an avg VM stress for the element as
   * considering all nodal stresses would use too much memory. ---*/
  AD::SetPreaccOut(elStress);
  AD::EndPreacc();

  return elStress;
}


CFEAMeshElasticity::CFEAMeshElasticity(unsigned short val_nDim, unsigned short val_nVar,
                                       unsigned long val_nElem, const CConfig *config) :
                                       CFEALinearElasticity() {
  nDim = val_nDim;
  nVar = val_nVar;

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

  if (element_based){
    E_i.reset(new su2double[val_nElem]);
    for (unsigned long iVar = 0; iVar < val_nElem; iVar++){
      E_i[iVar] = E;
    }
  }

}
