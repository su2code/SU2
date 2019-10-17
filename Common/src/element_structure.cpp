/*!
 * \file element_structure.cpp
 * \brief Definition of the Finite Element structure (elements)
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

#include "../include/element_structure.hpp"

unsigned short CElement::nDim = 0;

CElement::CElement(void) {
  
  CurrentCoord = NULL;
  RefCoord = NULL;
  GaussWeight = NULL;
  GaussCoord = NULL;
  
  GaussPoint = NULL;
  
  NodalStress = NULL;
  NodalExtrap = NULL;
  
  nNodes = 0;
  nGaussPoints = 0;
  
  el_Pressure = 0.0;
  
  Mab = NULL;
  Kab = NULL;
  Ks_ab = NULL;
  Kt_a = NULL;
  dNiXj = NULL;
  
  FDL_a = NULL;
  
  iDe = 0;
  iDV = 0;
  iProp = 0;
}


CElement::CElement(unsigned short val_nDim, CConfig *config) {
  
  /*--- Initializate the number of dimension and some structures we need ---*/
  nDim = val_nDim;
  
  CurrentCoord = NULL;
  RefCoord = NULL;
  GaussWeight = NULL;
  GaussCoord = NULL;
  
  GaussPoint = NULL;
  
  NodalStress = NULL;
  NodalExtrap = NULL;
  
  nNodes = 0;
  nGaussPoints = 0;
  
  el_Pressure = 0.0;
  
  Mab = NULL;
  Kab = NULL;
  Ks_ab = NULL;
  Kt_a = NULL;
  dNiXj = NULL;
  
  FDL_a = NULL;
  
  iDe = 0;
  iDV = 0;
  iProp = 0;

}

CElement::~CElement(void) {
  
  unsigned short iNode, iGauss, jNode;
  
  if (GaussPoint != NULL) {
    for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
      delete GaussPoint[iGauss];
    }
    delete [] GaussPoint;
  }
  
  if (NodalExtrap != NULL) {
    for (iNode = 0; iNode < nNodes; iNode++) {
      delete [] NodalExtrap[iNode];
    }
    delete [] NodalExtrap;
  }
  
  if (NodalStress != NULL) {
    for (iNode = 0; iNode < nNodes; iNode++) {
      delete [] NodalStress[iNode];
    }
    delete [] NodalStress;
  }
  
  if (CurrentCoord != NULL) {
    for (iNode = 0; iNode < nNodes; iNode++) {
      delete [] CurrentCoord [iNode];
    }
    delete [] CurrentCoord;
  }
  
  if (RefCoord != NULL) {
    for (iNode = 0; iNode < nNodes; iNode++) {
      delete [] RefCoord [iNode];
    }
    delete [] RefCoord;
  }
  
  if (GaussWeight != NULL) {
    delete [] GaussWeight;
  }
  
  if (GaussCoord != NULL) {
    for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
      delete [] GaussCoord[iGauss];
    }
    delete [] GaussCoord;
  }
  
  if (Mab != NULL) {
    for (iNode = 0; iNode < nNodes; iNode++) {
      delete [] Mab[iNode];
    }
    delete [] Mab;
  }

  if (Kab != NULL) {
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (jNode = 0; jNode < nNodes; jNode++) {
        delete [] Kab [iNode][jNode];
      }
      delete [] Kab[iNode];
    }
    delete [] Kab;
  }

  if (Ks_ab != NULL) {
    for (iNode = 0; iNode < nNodes; iNode++) {
      delete [] Ks_ab[iNode];
    }
    delete [] Ks_ab;
  }
  
  if (Kt_a != NULL) {
    for (iNode = 0; iNode < nNodes; iNode++) {
      delete [] Kt_a[iNode];
    }
    delete [] Kt_a;
  }

  if (dNiXj != NULL) {
    for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
      for (iNode = 0; iNode < nNodes; iNode++) {
        delete [] dNiXj [iGauss][iNode];
      }
      delete [] dNiXj[iGauss];
    }
    delete [] dNiXj;
  }

  if (FDL_a != NULL) {
    for (iNode = 0; iNode < nNodes; iNode++) {
      delete [] FDL_a[iNode];
    }
    delete [] FDL_a;
  }

}

void CElement::AllocateStructures(const bool body_forces) {

  /*--- Derived classes should call this method after setting nGauss and nNodes. ---*/

  unsigned short iNode, jNode, iGauss, nDimSq = nDim*nDim;

  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }

  CurrentCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    CurrentCoord [iNode] = new su2double[nDim];
  }

  RefCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    RefCoord [iNode] = new su2double[nDim];
  }

  GaussWeight = new su2double [nGaussPoints];

  GaussCoord = new su2double*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussCoord [iGauss] = new su2double[nDim];
  }

  dNiXj = new su2double **[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    dNiXj [iGauss] = new su2double*[nNodes];
    for (iNode = 0; iNode < nNodes; iNode++) {
      dNiXj [iGauss][iNode] = new su2double[nDim];
    }
  }

  NodalExtrap = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalExtrap[iNode] = new su2double[nGaussPoints];
  }

  NodalStress = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalStress[iNode] = new su2double[6];
  }

  Mab = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Mab[iNode] = new su2double [nNodes];
  }

  Kab = new su2double **[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kab [iNode] = new su2double*[nNodes];
    for (jNode = 0; jNode < nNodes; jNode++) {
      Kab [iNode][jNode] = new su2double[nDimSq];
    }
  }

  Ks_ab = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Ks_ab[iNode] = new su2double [nNodes];
  }

  Kt_a = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kt_a[iNode] = new su2double [nDim];
  }

  /*--- Body forces ---*/

  if (body_forces) {

    FDL_a = new su2double *[nNodes];
    for (iNode = 0; iNode < nNodes; iNode++) {
      FDL_a[iNode] = new su2double [nDim];
    }

  }

}

void CElement::Add_Kab(su2double **val_Kab, unsigned short nodeA, unsigned short nodeB) {
  
  unsigned short iDim, jDim;
  
  for(iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Kab[nodeA][nodeB][iDim*nDim+jDim] += val_Kab[iDim][jDim];
    }
  }
}

void CElement::Add_Kab_T(su2double **val_Kab, unsigned short nodeA, unsigned short nodeB) {
  
  unsigned short iDim, jDim;
  
  for(iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Kab[nodeA][nodeB][iDim*nDim+jDim] += val_Kab[jDim][iDim];
    }
  }
}

void CElement::Add_Kt_a(su2double *val_Kt_a, unsigned short nodeA) {
  
  unsigned short iDim;
  
  for(iDim = 0; iDim < nDim; iDim++) {
    Kt_a[nodeA][iDim] += val_Kt_a[iDim];
  }
  
}

void CElement::Add_FDL_a(su2double *val_FDL_a, unsigned short nodeA) {
  
  unsigned short iDim;
  
  for(iDim = 0; iDim < nDim; iDim++) {
    FDL_a[nodeA][iDim] += val_FDL_a[iDim];
  }
  
}


void CElement::clearElement(void) {
  
  unsigned short iNode, jNode, iDim, nDimSq;
  
  nDimSq = nDim*nDim;
  
  for(iNode = 0; iNode < nNodes; iNode++) {
    for(iDim = 0; iDim < nDim; iDim++) {
      if (Kt_a != NULL) Kt_a[iNode][iDim] = 0.0;
      if (FDL_a != NULL) FDL_a[iNode][iDim] = 0.0;
    }
    for (jNode = 0; jNode < nNodes; jNode++) {
      if (Ks_ab != NULL) Ks_ab[iNode][jNode] = 0.0;
      if (Mab != NULL) Mab[iNode][jNode] = 0.0;
      for(iDim = 0; iDim < nDimSq; iDim++) {
        if (Kab != NULL) Kab[iNode][jNode][iDim] = 0.0;
      }
    }
  }
}

void CElement::clearStress(void) {
  
  unsigned short iNode, iStress, nStress;
  
  if (nDim == 2) nStress = 3;
  else nStress = 6;
  
  for(iNode = 0; iNode < nNodes; iNode++) {
    for (iStress = 0; iStress < nStress; iStress++) {
      NodalStress[iNode][iStress] = 0.0;
    }
  }
  
}

void CElement::Set_ElProperties(CProperty *input_element) {

  /*--- Set the properties local to the element ---*/

  iDV = input_element->GetDV();
  iProp = input_element->GetMat_Prop();
  iDe = input_element->GetElectric_Prop();

}

void CElement::ComputeGrad_Linear(void) {

  if (nDim==2)
    ComputeGrad_2D(REFERENCE);
  else
    ComputeGrad_3D(REFERENCE);

}

void CElement::ComputeGrad_NonLinear(void) {

  if (nDim==2) {
    ComputeGrad_2D(REFERENCE);
    ComputeGrad_2D(CURRENT);
  }
  else {
    ComputeGrad_3D(REFERENCE);
    ComputeGrad_3D(CURRENT);
  }

}

void CElement::ComputeGrad_2D(const FrameType mode) {

  su2double Jacobian[2][2], ad[2][2];
  su2double detJac, GradNi_Xj;
  unsigned short iNode, iDim, jDim, iGauss;

  /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
        for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed) ---*/
  su2double **Coord = (mode==REFERENCE) ? RefCoord : CurrentCoord;

  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {

    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/

    for (iDim = 0; iDim < 2; iDim++)
      for (jDim = 0; jDim < 2; jDim++)
        Jacobian[iDim][jDim] = 0.0;

    for (iNode = 0; iNode < nNodes; iNode++)
      for (iDim = 0; iDim < 2; iDim++)
        for (jDim = 0; jDim < 2; jDim++)
          Jacobian[iDim][jDim] += Coord[iNode][jDim] * dNiXj[iGauss][iNode][iDim];

    /*--- Adjoint to Jacobian ---*/

    ad[0][0] =  Jacobian[1][1];  ad[0][1] = -Jacobian[0][1];
    ad[1][0] = -Jacobian[1][0];  ad[1][1] =  Jacobian[0][0];

    /*--- Determinant of Jacobian ---*/

    detJac = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];

    if (mode==REFERENCE)
      GaussPoint[iGauss]->SetJ_X(detJac);
    else
      GaussPoint[iGauss]->SetJ_x(detJac);

    /*--- Jacobian inverse (it was already computed as transpose) ---*/

    for (iDim = 0; iDim < 2; iDim++)
      for (jDim = 0; jDim < 2; jDim++)
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;

    /*--- Derivatives with respect to global coordinates ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < 2; iDim++) {
        GradNi_Xj = 0.0;
        for (jDim = 0; jDim < 2; jDim++)
          GradNi_Xj += Jacobian[iDim][jDim] * dNiXj[iGauss][iNode][jDim];

        if (mode==REFERENCE)
          GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
        else
          GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj, iDim, iNode);
      }
    }

  }

}

void CElement::ComputeGrad_3D(const FrameType mode) {

  su2double Jacobian[3][3], ad[3][3];
  su2double detJac, GradNi_Xj;
  unsigned short iNode, iDim, jDim, iGauss;

  /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
        for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed) ---*/
  su2double **Coord = (mode==REFERENCE) ? RefCoord : CurrentCoord;

  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {

    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/

    for (iDim = 0; iDim < 3; iDim++)
      for (jDim = 0; jDim < 3; jDim++)
        Jacobian[iDim][jDim] = 0.0;

    for (iNode = 0; iNode < nNodes; iNode++)
      for (iDim = 0; iDim < 3; iDim++)
        for (jDim = 0; jDim < 3; jDim++)
          Jacobian[iDim][jDim] += Coord[iNode][jDim] * dNiXj[iGauss][iNode][iDim];

    /*--- Adjoint to Jacobian ---*/

    ad[0][0] = Jacobian[1][1]*Jacobian[2][2]-Jacobian[1][2]*Jacobian[2][1];
    ad[0][1] = Jacobian[0][2]*Jacobian[2][1]-Jacobian[0][1]*Jacobian[2][2];
    ad[0][2] = Jacobian[0][1]*Jacobian[1][2]-Jacobian[0][2]*Jacobian[1][1];
    ad[1][0] = Jacobian[1][2]*Jacobian[2][0]-Jacobian[1][0]*Jacobian[2][2];
    ad[1][1] = Jacobian[0][0]*Jacobian[2][2]-Jacobian[0][2]*Jacobian[2][0];
    ad[1][2] = Jacobian[0][2]*Jacobian[1][0]-Jacobian[0][0]*Jacobian[1][2];
    ad[2][0] = Jacobian[1][0]*Jacobian[2][1]-Jacobian[1][1]*Jacobian[2][0];
    ad[2][1] = Jacobian[0][1]*Jacobian[2][0]-Jacobian[0][0]*Jacobian[2][1];
    ad[2][2] = Jacobian[0][0]*Jacobian[1][1]-Jacobian[0][1]*Jacobian[1][0];

    /*--- Determinant of Jacobian ---*/

    detJac = Jacobian[0][0]*ad[0][0]+Jacobian[0][1]*ad[1][0]+Jacobian[0][2]*ad[2][0];

    if (mode==REFERENCE)
      GaussPoint[iGauss]->SetJ_X(detJac);
    else
      GaussPoint[iGauss]->SetJ_x(detJac);

    /*--- Jacobian inverse (it was already computed as transpose) ---*/

    for (iDim = 0; iDim < 3; iDim++)
      for (jDim = 0; jDim < 3; jDim++)
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;

    /*--- Derivatives with respect to global coordinates ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < 3; iDim++) {
        GradNi_Xj = 0.0;
        for (jDim = 0; jDim < 3; jDim++)
          GradNi_Xj += Jacobian[iDim][jDim] * dNiXj[iGauss][iNode][jDim];

        if (mode==REFERENCE)
          GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
        else
          GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj, iDim, iNode);
      }
    }

  }

}
