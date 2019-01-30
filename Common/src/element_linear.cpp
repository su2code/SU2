/*!
 * \file element_linear.cpp
 * \brief Definition of the linear element structure for structural applications
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

CTRIA1::CTRIA1(void) : CElement() {
  
}

CTRIA1::CTRIA1(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {
  
  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;
  
  bool body_forces = config->GetDeadLoad();	// Body forces (dead loads).
  
  nNodes = 3;
  nGaussPoints = 1;
  
  nDimSq = nDim*nDim;
  
  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }
  
  NodalExtrap = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalExtrap[iNode] = new su2double[nGaussPoints];
  }
  
  NodalStress = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalStress[iNode] = new su2double[3];
  }
  
  /*--- Initialize structure for current and reference configuration ---*/
  
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
  
  GaussCoord[0][0] = 0.333333333333333;  GaussCoord[0][1] = 0.333333333333333;  GaussWeight[0] = 0.5;
  
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
  
  if (body_forces) {
    FDL_a = new su2double *[nNodes];
    for (iNode = 0; iNode < nNodes; iNode++) {
      FDL_a[iNode] = new su2double [nDim];
    }
  }
  else {
    FDL_a = NULL;
  }
  
  su2double Xi, Eta, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    
    val_Ni = Xi;
    GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = Eta;
    GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 1-Xi-Eta;
    GaussPoint[iGauss]->SetNi(val_Ni,2);
  }
  
  /*--- Shape functions evaluated at the nodes for extrapolation of the stresses at the Gaussian Points ---*/
  /*--- The stress is constant at a TRIA element ---*/
  NodalExtrap[0][0] = 1.0;
  NodalExtrap[1][0] = 1.0;
  NodalExtrap[2][0] = 1.0;
  
}

CTRIA1::~CTRIA1(void) {
  
}

void CTRIA1::ComputeGrad_Linear(void) {
  
  su2double Jacobian[2][2], dNiXj[3][2];
  su2double detJac, GradNi_Xj;
  su2double ad[2][2];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = 1.0; 	dNiXj[0][1] = 0.0;
    dNiXj[1][0] = 0.0; 	dNiXj[1][1] = 1.0;
    dNiXj[2][0] = -1.0; 	dNiXj[2][1] = -1.0;
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad[0][0] = Jacobian[1][1];
    ad[0][1] = -Jacobian[0][1];
    ad[1][0] = -Jacobian[1][0];
    ad[1][1] = Jacobian[0][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
      }
    }
  }
  
}

void CTRIA1::ComputeGrad_NonLinear(void) {
  
  su2double Jac_Ref[2][2], Jac_Curr[2][2], dNiXj[3][2];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[2][2], ad_Curr[2][2];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = 1.0; 	dNiXj[0][1] = 0.0;
    dNiXj[1][0] = 0.0; 	dNiXj[1][1] = 1.0;
    dNiXj[2][0] = -1.0; 	dNiXj[2][1] = -1.0;
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad_Ref[0][0] = Jac_Ref[1][1];
    ad_Ref[0][1] = -Jac_Ref[0][1];
    ad_Ref[1][0] = -Jac_Ref[1][0];
    ad_Ref[1][1] = Jac_Ref[0][0];
    
    ad_Curr[0][0] = Jac_Curr[1][1];
    ad_Curr[0][1] = -Jac_Curr[0][1];
    ad_Curr[1][0] = -Jac_Curr[1][0];
    ad_Curr[1][1] = Jac_Curr[0][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac_Ref = ad_Ref[0][0]*ad_Ref[1][1]-ad_Ref[0][1]*ad_Ref[1][0];
    detJac_Curr = ad_Curr[0][0]*ad_Curr[1][1]-ad_Curr[0][1]*ad_Curr[1][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac_Ref);
    GaussPoint[iGauss]->SetJ_x(detJac_Curr);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }
  
  
}
  
su2double CTRIA1::ComputeArea(void){
  
  unsigned short iDim;
  su2double a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0};
  su2double Area = 0.0;
  
  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = RefCoord[0][iDim]-RefCoord[2][iDim];
    b[iDim] = RefCoord[1][iDim]-RefCoord[2][iDim];
  }
  
  Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  
  return Area;
  
}

su2double CTRIA1::ComputeCurrentArea(void){

  unsigned short iDim;
  su2double a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0};
  su2double Area = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = CurrentCoord[0][iDim]-CurrentCoord[2][iDim];
    b[iDim] = CurrentCoord[1][iDim]-CurrentCoord[2][iDim];
  }

  Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

  return Area;

}
  
CQUAD4::CQUAD4(void) : CElement() {
  
}

CQUAD4::CQUAD4(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {
  
  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;
  
  bool body_forces = config->GetDeadLoad();	// Body forces (dead loads).
  
  nNodes = 4;
  nGaussPoints = 4;
  
  nDimSq = nDim*nDim;
  
  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }
  
  NodalExtrap = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalExtrap[iNode] = new su2double[nGaussPoints];
  }
  
  NodalStress = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalStress[iNode] = new su2double[3];
  }
  
  /*--- Initialize structure for current and reference configuration ---*/
  
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
  
  GaussCoord[0][0] = -0.577350269189626;  GaussCoord[0][1] = -0.577350269189626;  GaussWeight[0] = 1.0;
  GaussCoord[1][0] = 0.577350269189626;   GaussCoord[1][1] = -0.577350269189626;  GaussWeight[1] = 1.0;
  GaussCoord[2][0] = 0.577350269189626;   GaussCoord[2][1] = 0.577350269189626;   GaussWeight[2] = 1.0;
  GaussCoord[3][0] = -0.577350269189626;  GaussCoord[3][1] = 0.577350269189626;   GaussWeight[3] = 1.0;
  
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
  
  if (body_forces) {
    FDL_a = new su2double *[nNodes];
    for (iNode = 0; iNode < nNodes; iNode++) {
      FDL_a[iNode] = new su2double [nDim];
    }
  }
  else {
    FDL_a = NULL;
  }
  
  /*--- Store the shape functions (they only need to be computed once) ---*/
  su2double Xi, Eta, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    
    val_Ni = 0.25*(1.0-Xi)*(1.0-Eta);		GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = 0.25*(1.0+Xi)*(1.0-Eta);		GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 0.25*(1.0+Xi)*(1.0+Eta);		GaussPoint[iGauss]->SetNi(val_Ni,2);
    val_Ni = 0.25*(1.0-Xi)*(1.0+Eta);		GaussPoint[iGauss]->SetNi(val_Ni,3);
  }
  
  su2double ExtrapCoord[4][2];
  
  ExtrapCoord[0][0] = -1.732050807568877;  ExtrapCoord[0][1] = -1.732050807568877;
  ExtrapCoord[1][0] = 1.732050807568877;   ExtrapCoord[1][1] = -1.732050807568877;
  ExtrapCoord[2][0] = 1.732050807568877;   ExtrapCoord[2][1] = 1.732050807568877;
  ExtrapCoord[3][0] = -1.732050807568877;  ExtrapCoord[3][1] = 1.732050807568877;
  
  /*--- Store the shape functions (they only need to be computed once) ---*/
  for (iNode = 0; iNode < nNodes; iNode++) {
    Xi = ExtrapCoord[iNode][0];
    Eta = ExtrapCoord[iNode][1];
    
    NodalExtrap[iNode][0] = 0.25*(1.0-Xi)*(1.0-Eta);
    NodalExtrap[iNode][1] = 0.25*(1.0+Xi)*(1.0-Eta);
    NodalExtrap[iNode][2] = 0.25*(1.0+Xi)*(1.0+Eta);
    NodalExtrap[iNode][3] = 0.25*(1.0-Xi)*(1.0+Eta);
    
  }
  
}

CQUAD4::~CQUAD4(void) {
  
}

void CQUAD4::ComputeGrad_Linear(void) {
  
  su2double Xi, Eta;
  su2double Jacobian[2][2], dNiXj[4][2];
  su2double detJac, GradNi_Xj;
  su2double ad[2][2];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = -0.25*(1.0-Eta); dNiXj[0][1] = -0.25*(1.0-Xi);
    dNiXj[1][0] =  0.25*(1.0-Eta); dNiXj[1][1] = -0.25*(1.0+Xi);
    dNiXj[2][0] =  0.25*(1.0+Eta); dNiXj[2][1] =  0.25*(1.0+Xi);
    dNiXj[3][0] = -0.25*(1.0+Eta); dNiXj[3][1] =  0.25*(1.0-Xi);
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad[0][0] = Jacobian[1][1];
    ad[0][1] = -Jacobian[0][1];
    ad[1][0] = -Jacobian[1][0];
    ad[1][1] = Jacobian[0][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
      }
    }
  }
  
}

void CQUAD4::ComputeGrad_NonLinear(void) {
  
  su2double Xi, Eta;
  su2double Jac_Ref[2][2], Jac_Curr[2][2], dNiXj[4][2];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[2][2], ad_Curr[2][2];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = -0.25*(1.0-Eta); dNiXj[0][1] = -0.25*(1.0-Xi);
    dNiXj[1][0] =  0.25*(1.0-Eta); dNiXj[1][1] = -0.25*(1.0+Xi);
    dNiXj[2][0] =  0.25*(1.0+Eta); dNiXj[2][1] =  0.25*(1.0+Xi);
    dNiXj[3][0] = -0.25*(1.0+Eta); dNiXj[3][1] =  0.25*(1.0-Xi);
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad_Ref[0][0] = Jac_Ref[1][1];
    ad_Ref[0][1] = -Jac_Ref[0][1];
    ad_Ref[1][0] = -Jac_Ref[1][0];
    ad_Ref[1][1] = Jac_Ref[0][0];
    
    ad_Curr[0][0] = Jac_Curr[1][1];
    ad_Curr[0][1] = -Jac_Curr[0][1];
    ad_Curr[1][0] = -Jac_Curr[1][0];
    ad_Curr[1][1] = Jac_Curr[0][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac_Ref = ad_Ref[0][0]*ad_Ref[1][1]-ad_Ref[0][1]*ad_Ref[1][0];
    detJac_Curr = ad_Curr[0][0]*ad_Curr[1][1]-ad_Curr[0][1]*ad_Curr[1][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac_Ref);
    GaussPoint[iGauss]->SetJ_x(detJac_Curr);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }
  
}
  
su2double CQUAD4::ComputeArea(void){
  
  unsigned short iDim;
  su2double a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0};
  su2double Area = 0.0;
  
  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = RefCoord[0][iDim]-RefCoord[2][iDim];
    b[iDim] = RefCoord[1][iDim]-RefCoord[2][iDim];
  }
  
  Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  
  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = RefCoord[0][iDim]-RefCoord[3][iDim];
    b[iDim] = RefCoord[2][iDim]-RefCoord[3][iDim];
  }
  
  Area += 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  
  return Area;
  
}

su2double CQUAD4::ComputeCurrentArea(void){

  unsigned short iDim;
  su2double a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0};
  su2double Area = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = CurrentCoord[0][iDim]-CurrentCoord[2][iDim];
    b[iDim] = CurrentCoord[1][iDim]-CurrentCoord[2][iDim];
  }

  Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = CurrentCoord[0][iDim]-CurrentCoord[3][iDim];
    b[iDim] = CurrentCoord[2][iDim]-CurrentCoord[3][iDim];
  }

  Area += 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

  return Area;

}
  
  
CQUAD1::CQUAD1(void) : CElement() {

}

CQUAD1::CQUAD1(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {

	unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;

  nNodes = 4;
  nGaussPoints = 1;

  nDimSq = nDim*nDim;

  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }

  NodalExtrap = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalExtrap[iNode] = new su2double[nGaussPoints];
  }

  NodalStress = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalStress[iNode] = new su2double[3];
  }

  /*--- Initialize structure for current and reference configuration ---*/

  CurrentCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++){
    CurrentCoord [iNode] = new su2double[nDim];
  }

  RefCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++){
    RefCoord [iNode] = new su2double[nDim];
  }

  GaussWeight = new su2double [nGaussPoints];

  GaussCoord = new su2double*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++){
    GaussCoord [iGauss] = new su2double[nDim];
  }

  GaussCoord[0][0] = 0.0;  GaussCoord[0][1] = 0.0;  GaussWeight[0] = 4.0;

  Kk_ab = new su2double **[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kk_ab [iNode] = new su2double*[nNodes];
    for (jNode = 0; jNode < nNodes; jNode++) {
      Kk_ab [iNode][jNode] = new su2double[nDimSq];
    }
  }

  /*--- Store the shape functions (they only need to be computed once) ---*/
  su2double Xi, Eta, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++){
      Xi = GaussCoord[iGauss][0];
      Eta = GaussCoord[iGauss][1];

      val_Ni = 0.25*(1.0-Xi)*(1.0-Eta);    GaussPoint[iGauss]->SetNi(val_Ni,0);
      val_Ni = 0.25*(1.0+Xi)*(1.0-Eta);    GaussPoint[iGauss]->SetNi(val_Ni,1);
      val_Ni = 0.25*(1.0+Xi)*(1.0+Eta);    GaussPoint[iGauss]->SetNi(val_Ni,2);
      val_Ni = 0.25*(1.0-Xi)*(1.0+Eta);    GaussPoint[iGauss]->SetNi(val_Ni,3);
  }

  /*--- Shape functions evaluated at the nodes for extrapolation of the stresses at the Gaussian Points ---*/
  /*--- The stress is constant at a QUAD1 element ---*/
  NodalExtrap[0][0] = 1.0;
  NodalExtrap[1][0] = 1.0;
  NodalExtrap[2][0] = 1.0;
  NodalExtrap[3][0] = 1.0;

}

CQUAD1::~CQUAD1(void) {

  unsigned short iVar, jVar;

  for (iVar = 0; iVar < nGaussPoints; iVar++){
    delete [] GaussCoord[iVar];
    delete [] GaussPoint[iVar];
  }

  for (iVar = 0; iVar < nNodes; iVar++){
    for (jVar = 0; jVar < nNodes; jVar++){
      delete [] Kab[iVar][jVar];
    }
    delete [] CurrentCoord[iVar];
    delete [] RefCoord[iVar];
    delete [] Mab[iVar];
    delete [] Kab[iVar];
    delete [] Ks_ab[iVar];
    delete [] Kt_a[iVar];
    if (FDL_a != NULL) delete [] FDL_a[iVar];
    delete [] NodalExtrap[iVar];
  }

  delete [] GaussCoord;
  delete [] GaussPoint;
  delete [] CurrentCoord;
  delete [] RefCoord;
  delete [] Mab;
  delete [] Kab;
  delete [] Ks_ab;
  delete [] Kt_a;
  delete [] GaussWeight;
  delete [] NodalExtrap;

  if (FDL_a != NULL) delete [] FDL_a;

}

void CQUAD1::ComputeGrad_Linear(void){

    su2double Xi, Eta;
    su2double Jacobian[2][2], dNiXj[4][2];
    su2double detJac, GradNi_Xj;
    su2double ad[2][2];
    unsigned short iNode, iDim, jDim, iGauss;

    for (iGauss = 0; iGauss < nGaussPoints; iGauss++){

      Xi = GaussCoord[iGauss][0];
      Eta = GaussCoord[iGauss][1];

      /*--- dN/d xi, dN/d eta ---*/

      dNiXj[0][0] = -0.25*(1.0-Eta); dNiXj[0][1] = -0.25*(1.0-Xi);
      dNiXj[1][0] =  0.25*(1.0-Eta); dNiXj[1][1] = -0.25*(1.0+Xi);
      dNiXj[2][0] =  0.25*(1.0+Eta); dNiXj[2][1] =  0.25*(1.0+Xi);
      dNiXj[3][0] = -0.25*(1.0+Eta); dNiXj[3][1] =  0.25*(1.0-Xi);

      /*--- Jacobian transformation ---*/
      /*--- This does dX/dXi transpose ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
      }

      /*--- Adjoint to Jacobian ---*/

      ad[0][0] = Jacobian[1][1];
      ad[0][1] = -Jacobian[0][1];
      ad[1][0] = -Jacobian[1][0];
      ad[1][1] = Jacobian[0][0];

      /*--- Determinant of Jacobian ---*/

      detJac = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];

      GaussPoint[iGauss]->SetJ_X(detJac);

      /*--- Jacobian inverse (it was already computed as transpose) ---*/

      for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
      }

      /*--- Derivatives with respect to global coordinates ---*/

      for (iNode = 0; iNode < nNodes; iNode++) {
        for (iDim = 0; iDim < nDim; iDim++){
          GradNi_Xj = 0.0;
          for (jDim = 0; jDim < nDim; jDim++){
            GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
          }
          GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
        }
      }
    }

}

void CQUAD1::ComputeGrad_NonLinear(void){

    su2double Xi, Eta;
    su2double Jac_Ref[2][2], Jac_Curr[2][2], dNiXj[4][2];
    su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
    su2double ad_Ref[2][2], ad_Curr[2][2];
    unsigned short iNode, iDim, jDim, iGauss;

    for (iGauss = 0; iGauss < nGaussPoints; iGauss++){

      Xi = GaussCoord[iGauss][0];
      Eta = GaussCoord[iGauss][1];

      /*--- dN/d xi, dN/d eta ---*/

      dNiXj[0][0] = -0.25*(1.0-Eta); dNiXj[0][1] = -0.25*(1.0-Xi);
      dNiXj[1][0] =  0.25*(1.0-Eta); dNiXj[1][1] = -0.25*(1.0+Xi);
      dNiXj[2][0] =  0.25*(1.0+Eta); dNiXj[2][1] =  0.25*(1.0+Xi);
      dNiXj[3][0] = -0.25*(1.0+Eta); dNiXj[3][1] =  0.25*(1.0-Xi);

      /*--- Jacobian transformation ---*/
      /*--- This does dX/dXi transpose ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
      }

      /*--- Adjoint to Jacobian ---*/

      ad_Ref[0][0] = Jac_Ref[1][1];
      ad_Ref[0][1] = -Jac_Ref[0][1];
      ad_Ref[1][0] = -Jac_Ref[1][0];
      ad_Ref[1][1] = Jac_Ref[0][0];

      ad_Curr[0][0] = Jac_Curr[1][1];
      ad_Curr[0][1] = -Jac_Curr[0][1];
      ad_Curr[1][0] = -Jac_Curr[1][0];
      ad_Curr[1][1] = Jac_Curr[0][0];

      /*--- Determinant of Jacobian ---*/

      detJac_Ref = ad_Ref[0][0]*ad_Ref[1][1]-ad_Ref[0][1]*ad_Ref[1][0];
      detJac_Curr = ad_Curr[0][0]*ad_Curr[1][1]-ad_Curr[0][1]*ad_Curr[1][0];

      GaussPoint[iGauss]->SetJ_X(detJac_Ref);
      GaussPoint[iGauss]->SetJ_x(detJac_Curr);

      /*--- Jacobian inverse (it was already computed as transpose) ---*/

      for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
      }

      /*--- Derivatives with respect to global coordinates ---*/

      for (iNode = 0; iNode < nNodes; iNode++) {
        for (iDim = 0; iDim < nDim; iDim++){
          GradNi_Xj_Ref = 0.0;
          GradNi_Xj_Curr = 0.0;
          for (jDim = 0; jDim < nDim; jDim++){
            GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
            GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
          }
          GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
          GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
        }
      }
    }

}

void CQUAD1::ComputeGrad_Pressure(void){

    su2double Xi, Eta;
    su2double Jac_Ref[2][2], Jac_Curr[2][2], dNiXj[4][2];
    su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
    su2double ad_Ref[2][2], ad_Curr[2][2];
    unsigned short iNode, iDim, jDim, iGauss;

    for (iGauss = 0; iGauss < nGaussPoints; iGauss++){

      Xi = GaussCoord[iGauss][0];
      Eta = GaussCoord[iGauss][1];

      /*--- dN/d xi, dN/d eta ---*/

      dNiXj[0][0] = -0.25*(1.0-Eta); dNiXj[0][1] = -0.25*(1.0-Xi);
      dNiXj[1][0] =  0.25*(1.0-Eta); dNiXj[1][1] = -0.25*(1.0+Xi);
      dNiXj[2][0] =  0.25*(1.0+Eta); dNiXj[2][1] =  0.25*(1.0+Xi);
      dNiXj[3][0] = -0.25*(1.0+Eta); dNiXj[3][1] =  0.25*(1.0-Xi);

      /*--- Jacobian transformation ---*/
      /*--- This does dX/dXi transpose ---*/

      for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < 4; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
      }

      /*--- Adjoint to Jacobian ---*/

      ad_Ref[0][0] = Jac_Ref[1][1];
      ad_Ref[0][1] = -Jac_Ref[0][1];
      ad_Ref[1][0] = -Jac_Ref[1][0];
      ad_Ref[1][1] = Jac_Ref[0][0];

      ad_Curr[0][0] = Jac_Curr[1][1];
      ad_Curr[0][1] = -Jac_Curr[0][1];
      ad_Curr[1][0] = -Jac_Curr[1][0];
      ad_Curr[1][1] = Jac_Curr[0][0];

      /*--- Determinant of Jacobian ---*/

      detJac_Ref = ad_Ref[0][0]*ad_Ref[1][1]-ad_Ref[0][1]*ad_Ref[1][0];
      detJac_Curr = ad_Curr[0][0]*ad_Curr[1][1]-ad_Curr[0][1]*ad_Curr[1][0];

      GaussPoint[iGauss]->SetJ_X(detJac_Ref);
      GaussPoint[iGauss]->SetJ_x(detJac_Curr);

      /*--- Jacobian inverse (it was already computed as transpose) ---*/

      for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
      }

      /*--- Derivatives with respect to global coordinates ---*/

      for (iNode = 0; iNode < nNodes; iNode++) {
        for (iDim = 0; iDim < nDim; iDim++){
          GradNi_Xj_Ref = 0.0;
          GradNi_Xj_Curr = 0.0;
          for (jDim = 0; jDim < nDim; jDim++){
            GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
            GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
          }
          GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
          GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
        }
      }
    }

}

CTETRA1::CTETRA1(void) : CElement() {
  
}

CTETRA1::CTETRA1(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {
  
  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;
  
  bool body_forces = config->GetDeadLoad();	// Body forces (dead loads).
  
  nNodes = 4;
  nGaussPoints = 1;
  
  nDimSq = nDim*nDim;
  
  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }
  
  NodalExtrap = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalExtrap[iNode] = new su2double[nGaussPoints];
  }
  
  NodalStress = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalStress[iNode] = new su2double[6];
  }
  
  /*--- Initialize structure for current and reference configuration ---*/
  
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
  
  GaussCoord[0][0] = 0.25;  GaussCoord[0][1] = 0.25; GaussCoord[0][2] = 0.25;  GaussWeight[0] = 0.166666666666666;
  
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
  
  if (body_forces) {
    FDL_a = new su2double *[nNodes];
    for (iNode = 0; iNode < nNodes; iNode++) {
      FDL_a[iNode] = new su2double [nDim];
    }
  }
  else {
    FDL_a = NULL;
  }
  
  /*--- Store the shape functions (they only need to be computed once) ---*/
  su2double Xi, Eta, Zeta, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];
    
    val_Ni = Xi;						GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = Eta;						GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 1.0 - Xi - Eta - Zeta;	GaussPoint[iGauss]->SetNi(val_Ni,2);
    val_Ni = Zeta;					GaussPoint[iGauss]->SetNi(val_Ni,3);
  }
  
  /*--- Shape functions evaluated at the nodes for extrapolation of the stresses at the Gaussian Points ---*/
  /*--- The stress is constant at a TETRA element ---*/
  NodalExtrap[0][0] = 1.0;
  NodalExtrap[1][0] = 1.0;
  NodalExtrap[2][0] = 1.0;
  NodalExtrap[3][0] = 1.0;
  
}

CTETRA1::~CTETRA1(void) {
  
}

void CTETRA1::ComputeGrad_Linear(void) {
  
  su2double Jacobian[3][3], dNiXj[4][3];
  su2double detJac, GradNi_Xj;
  su2double ad[3][3];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = 1.0;   dNiXj[0][1] = 0.0;   dNiXj[0][2] = 0.0;
    dNiXj[1][0] = 0.0;   dNiXj[1][1] = 1.0;   dNiXj[1][2] = 0.0;
    dNiXj[2][0] = -1.0;  dNiXj[2][1] = -1.0;  dNiXj[2][2] = -1.0;
    dNiXj[3][0] = 0.0;   dNiXj[3][1] = 0.0;   dNiXj[3][2] = 1.0;
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
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
    
    GaussPoint[iGauss]->SetJ_X(detJac);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
      }
    }
  }
  
}

void CTETRA1::ComputeGrad_NonLinear(void) {
  
  su2double Jac_Ref[3][3], Jac_Curr[3][3], dNiXj[4][3];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[3][3], ad_Curr[3][3];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = 1.0;   dNiXj[0][1] = 0.0;   dNiXj[0][2] = 0.0;
    dNiXj[1][0] = 0.0;   dNiXj[1][1] = 1.0;   dNiXj[1][2] = 0.0;
    dNiXj[2][0] = -1.0;  dNiXj[2][1] = -1.0;  dNiXj[2][2] = -1.0;
    dNiXj[3][0] = 0.0;   dNiXj[3][1] = 0.0;   dNiXj[3][2] = 1.0;
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad_Ref[0][0] = Jac_Ref[1][1]*Jac_Ref[2][2]-Jac_Ref[1][2]*Jac_Ref[2][1];
    ad_Ref[0][1] = Jac_Ref[0][2]*Jac_Ref[2][1]-Jac_Ref[0][1]*Jac_Ref[2][2];
    ad_Ref[0][2] = Jac_Ref[0][1]*Jac_Ref[1][2]-Jac_Ref[0][2]*Jac_Ref[1][1];
    ad_Ref[1][0] = Jac_Ref[1][2]*Jac_Ref[2][0]-Jac_Ref[1][0]*Jac_Ref[2][2];
    ad_Ref[1][1] = Jac_Ref[0][0]*Jac_Ref[2][2]-Jac_Ref[0][2]*Jac_Ref[2][0];
    ad_Ref[1][2] = Jac_Ref[0][2]*Jac_Ref[1][0]-Jac_Ref[0][0]*Jac_Ref[1][2];
    ad_Ref[2][0] = Jac_Ref[1][0]*Jac_Ref[2][1]-Jac_Ref[1][1]*Jac_Ref[2][0];
    ad_Ref[2][1] = Jac_Ref[0][1]*Jac_Ref[2][0]-Jac_Ref[0][0]*Jac_Ref[2][1];
    ad_Ref[2][2] = Jac_Ref[0][0]*Jac_Ref[1][1]-Jac_Ref[0][1]*Jac_Ref[1][0];
    
    ad_Curr[0][0] = Jac_Curr[1][1]*Jac_Curr[2][2]-Jac_Curr[1][2]*Jac_Curr[2][1];
    ad_Curr[0][1] = Jac_Curr[0][2]*Jac_Curr[2][1]-Jac_Curr[0][1]*Jac_Curr[2][2];
    ad_Curr[0][2] = Jac_Curr[0][1]*Jac_Curr[1][2]-Jac_Curr[0][2]*Jac_Curr[1][1];
    ad_Curr[1][0] = Jac_Curr[1][2]*Jac_Curr[2][0]-Jac_Curr[1][0]*Jac_Curr[2][2];
    ad_Curr[1][1] = Jac_Curr[0][0]*Jac_Curr[2][2]-Jac_Curr[0][2]*Jac_Curr[2][0];
    ad_Curr[1][2] = Jac_Curr[0][2]*Jac_Curr[1][0]-Jac_Curr[0][0]*Jac_Curr[1][2];
    ad_Curr[2][0] = Jac_Curr[1][0]*Jac_Curr[2][1]-Jac_Curr[1][1]*Jac_Curr[2][0];
    ad_Curr[2][1] = Jac_Curr[0][1]*Jac_Curr[2][0]-Jac_Curr[0][0]*Jac_Curr[2][1];
    ad_Curr[2][2] = Jac_Curr[0][0]*Jac_Curr[1][1]-Jac_Curr[0][1]*Jac_Curr[1][0];
    
    
    /*--- Determinant of Jacobian ---*/
    
    detJac_Ref = Jac_Ref[0][0]*ad_Ref[0][0]+Jac_Ref[0][1]*ad_Ref[1][0]+Jac_Ref[0][2]*ad_Ref[2][0];
    detJac_Curr = Jac_Curr[0][0]*ad_Curr[0][0]+Jac_Curr[0][1]*ad_Curr[1][0]+Jac_Curr[0][2]*ad_Curr[2][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac_Ref);
    GaussPoint[iGauss]->SetJ_x(detJac_Curr);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }
  
}
  
su2double CTETRA1::ComputeVolume(void){

  unsigned short iDim;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0};
  su2double Volume = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[1][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[2][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[3][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  return Volume;

}

su2double CTETRA1::ComputeCurrentVolume(void){

  unsigned short iDim;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0};
  su2double Volume = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = CurrentCoord[1][iDim] - CurrentCoord[0][iDim];
    r2[iDim] = CurrentCoord[2][iDim] - CurrentCoord[0][iDim];
    r3[iDim] = CurrentCoord[3][iDim] - CurrentCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  return Volume;

}

CHEXA8::CHEXA8(void) : CElement() {
  
}

CHEXA8::CHEXA8(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {
  
  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;
  
  bool body_forces = config->GetDeadLoad();	// Body forces (dead loads).
  
  nNodes = 8;
  nGaussPoints = 8;
  
  nDimSq = nDim*nDim;
  
  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }
  
  NodalExtrap = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalExtrap[iNode] = new su2double[nGaussPoints];
  }
  
  NodalStress = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalStress[iNode] = new su2double[6];
  }
  
  /*--- Initialize structure for current and reference configuration ---*/
  
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
  
  GaussCoord[0][0] = -0.577350269189626;  GaussCoord[0][1] = -0.577350269189626;  GaussCoord[0][2] = -0.577350269189626;	GaussWeight[0] = 1.0;
  GaussCoord[1][0] = 0.577350269189626;   GaussCoord[1][1] = -0.577350269189626;  GaussCoord[1][2] = -0.577350269189626;  GaussWeight[1] = 1.0;
  GaussCoord[2][0] = 0.577350269189626;   GaussCoord[2][1] = 0.577350269189626;  	GaussCoord[2][2] = -0.577350269189626;  GaussWeight[2] = 1.0;
  GaussCoord[3][0] = -0.577350269189626;  GaussCoord[3][1] = 0.577350269189626;  	GaussCoord[3][2] = -0.577350269189626;  GaussWeight[3] = 1.0;
  GaussCoord[4][0] = -0.577350269189626;  GaussCoord[4][1] = -0.577350269189626;  GaussCoord[4][2] = 0.577350269189626;  	GaussWeight[4] = 1.0;
  GaussCoord[5][0] = 0.577350269189626;   GaussCoord[5][1] = -0.577350269189626;  GaussCoord[5][2] = 0.577350269189626;  	GaussWeight[5] = 1.0;
  GaussCoord[6][0] = 0.577350269189626;   GaussCoord[6][1] = 0.577350269189626;  	GaussCoord[6][2] = 0.577350269189626;  	GaussWeight[6] = 1.0;
  GaussCoord[7][0] = -0.577350269189626;  GaussCoord[7][1] = 0.577350269189626;  	GaussCoord[7][2] = 0.577350269189626;  	GaussWeight[7] = 1.0;
  
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
  
  if (body_forces) {
    FDL_a = new su2double *[nNodes];
    for (iNode = 0; iNode < nNodes; iNode++) {
      FDL_a[iNode] = new su2double [nDim];
    }
  }
  else {
    FDL_a = NULL;
  }
  
  
  /*--- Store the shape functions (they only need to be computed once) ---*/
  su2double Xi, Eta, Zeta, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];
    
    val_Ni = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0-Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0-Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0-Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,2);
    val_Ni = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0-Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,3);
    val_Ni = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0+Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,4);
    val_Ni = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0+Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,5);
    val_Ni = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0+Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,6);
    val_Ni = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0+Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,7);
  }
  
  
  su2double ExtrapCoord[8][3];
  
  ExtrapCoord[0][0] = -1.732050807568877;  ExtrapCoord[0][1] = -1.732050807568877;  	ExtrapCoord[0][2] = -1.732050807568877;
  ExtrapCoord[1][0] = 1.732050807568877;   ExtrapCoord[1][1] = -1.732050807568877;  	ExtrapCoord[1][2] = -1.732050807568877;
  ExtrapCoord[2][0] = 1.732050807568877;   ExtrapCoord[2][1] = 1.732050807568877;  	ExtrapCoord[2][2] = -1.732050807568877;
  ExtrapCoord[3][0] = -1.732050807568877;  ExtrapCoord[3][1] = 1.732050807568877;  	ExtrapCoord[3][2] = -1.732050807568877;
  ExtrapCoord[4][0] = -1.732050807568877;  ExtrapCoord[4][1] = -1.732050807568877;  	ExtrapCoord[4][2] = 1.732050807568877;
  ExtrapCoord[5][0] = 1.732050807568877;   ExtrapCoord[5][1] = -1.732050807568877;  	ExtrapCoord[5][2] = 1.732050807568877;
  ExtrapCoord[6][0] = 1.732050807568877;   ExtrapCoord[6][1] = 1.732050807568877;  	ExtrapCoord[6][2] = 1.732050807568877;
  ExtrapCoord[7][0] = -1.732050807568877;  ExtrapCoord[7][1] = 1.732050807568877;  	ExtrapCoord[7][2] = 1.732050807568877;
  
  
  /*--- Store the shape functions (they only need to be computed once) ---*/
  for (iNode = 0; iNode < nNodes; iNode++) {
    Xi = ExtrapCoord[iNode][0];
    Eta = ExtrapCoord[iNode][1];
    Zeta = ExtrapCoord[iNode][2];
    
    NodalExtrap[iNode][0] = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0-Zeta);
    NodalExtrap[iNode][1] = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0-Zeta);
    NodalExtrap[iNode][2] = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0-Zeta);
    NodalExtrap[iNode][3] = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0-Zeta);
    NodalExtrap[iNode][4] = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0+Zeta);
    NodalExtrap[iNode][5] = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0+Zeta);
    NodalExtrap[iNode][6] = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0+Zeta);
    NodalExtrap[iNode][7] = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0+Zeta);
  }
  
}

CHEXA8::~CHEXA8(void) {
  
}

void CHEXA8::ComputeGrad_Linear(void) {
  
  su2double Xi, Eta, Zeta;
  su2double Jacobian[3][3], dNiXj[8][3];
  su2double detJac, GradNi_Xj;
  su2double ad[3][3];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];
    
    /*--- dN/d xi ---*/
    
    dNiXj[0][0] = -0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[1][0] = 0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[2][0] = 0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[3][0] = -0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[4][0] = -0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[5][0] = 0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[6][0] = 0.125*(1.0+Eta)*(1.0+Zeta);
    dNiXj[7][0] = -0.125*(1.0+Eta)*(1.0+Zeta);
    
    /*--- dN/d eta ---*/
    
    dNiXj[0][1] = -0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[1][1] = -0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[2][1] = 0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[3][1] = 0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[4][1] = -0.125*(1.0-Xi)*(1.0+Zeta);
    dNiXj[5][1] = -0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[6][1] = 0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[7][1] = 0.125*(1.0-Xi)*(1.0+Zeta);
    
    /*--- dN/d mu ---*/
    
    dNiXj[0][2] = -0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[1][2] = -0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[2][2] = -0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[3][2] = -0.125*(1.0-Xi)*(1.0+Eta);
    dNiXj[4][2] = 0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[5][2] = 0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[6][2] = 0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[7][2] = 0.125*(1.0-Xi)*(1.0+Eta);
    
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
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
    
    GaussPoint[iGauss]->SetJ_X(detJac);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
      }
    }
  }
  
  
}

void CHEXA8::ComputeGrad_NonLinear(void) {
  
  su2double Xi, Eta, Zeta;
  su2double Jac_Ref[3][3], Jac_Curr[3][3], dNiXj[8][3];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[3][3], ad_Curr[3][3];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];
    
    /*--- dN/d xi, dN/d eta ---*/
    
    /*--- dN/d xi ---*/
    
    dNiXj[0][0] = -0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[1][0] = 0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[2][0] = 0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[3][0] = -0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[4][0] = -0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[5][0] = 0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[6][0] = 0.125*(1.0+Eta)*(1.0+Zeta);
    dNiXj[7][0] = -0.125*(1.0+Eta)*(1.0+Zeta);
    
    /*--- dN/d eta ---*/
    
    dNiXj[0][1] = -0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[1][1] = -0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[2][1] = 0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[3][1] = 0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[4][1] = -0.125*(1.0-Xi)*(1.0+Zeta);
    dNiXj[5][1] = -0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[6][1] = 0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[7][1] = 0.125*(1.0-Xi)*(1.0+Zeta);
    
    /*--- dN/d mu ---*/
    
    dNiXj[0][2] = -0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[1][2] = -0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[2][2] = -0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[3][2] = -0.125*(1.0-Xi)*(1.0+Eta);
    dNiXj[4][2] = 0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[5][2] = 0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[6][2] = 0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[7][2] = 0.125*(1.0-Xi)*(1.0+Eta);
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad_Ref[0][0] = Jac_Ref[1][1]*Jac_Ref[2][2]-Jac_Ref[1][2]*Jac_Ref[2][1];
    ad_Ref[0][1] = Jac_Ref[0][2]*Jac_Ref[2][1]-Jac_Ref[0][1]*Jac_Ref[2][2];
    ad_Ref[0][2] = Jac_Ref[0][1]*Jac_Ref[1][2]-Jac_Ref[0][2]*Jac_Ref[1][1];
    ad_Ref[1][0] = Jac_Ref[1][2]*Jac_Ref[2][0]-Jac_Ref[1][0]*Jac_Ref[2][2];
    ad_Ref[1][1] = Jac_Ref[0][0]*Jac_Ref[2][2]-Jac_Ref[0][2]*Jac_Ref[2][0];
    ad_Ref[1][2] = Jac_Ref[0][2]*Jac_Ref[1][0]-Jac_Ref[0][0]*Jac_Ref[1][2];
    ad_Ref[2][0] = Jac_Ref[1][0]*Jac_Ref[2][1]-Jac_Ref[1][1]*Jac_Ref[2][0];
    ad_Ref[2][1] = Jac_Ref[0][1]*Jac_Ref[2][0]-Jac_Ref[0][0]*Jac_Ref[2][1];
    ad_Ref[2][2] = Jac_Ref[0][0]*Jac_Ref[1][1]-Jac_Ref[0][1]*Jac_Ref[1][0];
    
    ad_Curr[0][0] = Jac_Curr[1][1]*Jac_Curr[2][2]-Jac_Curr[1][2]*Jac_Curr[2][1];
    ad_Curr[0][1] = Jac_Curr[0][2]*Jac_Curr[2][1]-Jac_Curr[0][1]*Jac_Curr[2][2];
    ad_Curr[0][2] = Jac_Curr[0][1]*Jac_Curr[1][2]-Jac_Curr[0][2]*Jac_Curr[1][1];
    ad_Curr[1][0] = Jac_Curr[1][2]*Jac_Curr[2][0]-Jac_Curr[1][0]*Jac_Curr[2][2];
    ad_Curr[1][1] = Jac_Curr[0][0]*Jac_Curr[2][2]-Jac_Curr[0][2]*Jac_Curr[2][0];
    ad_Curr[1][2] = Jac_Curr[0][2]*Jac_Curr[1][0]-Jac_Curr[0][0]*Jac_Curr[1][2];
    ad_Curr[2][0] = Jac_Curr[1][0]*Jac_Curr[2][1]-Jac_Curr[1][1]*Jac_Curr[2][0];
    ad_Curr[2][1] = Jac_Curr[0][1]*Jac_Curr[2][0]-Jac_Curr[0][0]*Jac_Curr[2][1];
    ad_Curr[2][2] = Jac_Curr[0][0]*Jac_Curr[1][1]-Jac_Curr[0][1]*Jac_Curr[1][0];
    
    
    /*--- Determinant of Jacobian ---*/
    
    detJac_Ref = Jac_Ref[0][0]*ad_Ref[0][0]+Jac_Ref[0][1]*ad_Ref[1][0]+Jac_Ref[0][2]*ad_Ref[2][0];
    detJac_Curr = Jac_Curr[0][0]*ad_Curr[0][0]+Jac_Curr[0][1]*ad_Curr[1][0]+Jac_Curr[0][2]*ad_Curr[2][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac_Ref);
    GaussPoint[iGauss]->SetJ_x(detJac_Curr);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }
  
  
}

su2double CHEXA8::ComputeVolume(void){

  unsigned short iDim;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0};
  su2double Volume = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[1][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[2][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[5][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[2][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[7][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[5][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[2][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[3][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[7][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[5][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[7][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[4][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[7][iDim] - RefCoord[2][iDim];
    r2[iDim] = RefCoord[5][iDim] - RefCoord[2][iDim];
    r3[iDim] = RefCoord[6][iDim] - RefCoord[2][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  return Volume;

}

su2double CHEXA8::ComputeCurrentVolume(void){

  unsigned short iDim;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0};
  su2double Volume = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = CurrentCoord[1][iDim] - CurrentCoord[0][iDim];
    r2[iDim] = CurrentCoord[2][iDim] - CurrentCoord[0][iDim];
    r3[iDim] = CurrentCoord[5][iDim] - CurrentCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = CurrentCoord[2][iDim] - CurrentCoord[0][iDim];
    r2[iDim] = CurrentCoord[7][iDim] - CurrentCoord[0][iDim];
    r3[iDim] = CurrentCoord[5][iDim] - CurrentCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = CurrentCoord[2][iDim] - CurrentCoord[0][iDim];
    r2[iDim] = CurrentCoord[3][iDim] - CurrentCoord[0][iDim];
    r3[iDim] = CurrentCoord[7][iDim] - CurrentCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = CurrentCoord[5][iDim] - CurrentCoord[0][iDim];
    r2[iDim] = CurrentCoord[7][iDim] - CurrentCoord[0][iDim];
    r3[iDim] = CurrentCoord[4][iDim] - CurrentCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = CurrentCoord[7][iDim] - CurrentCoord[2][iDim];
    r2[iDim] = CurrentCoord[5][iDim] - CurrentCoord[2][iDim];
    r3[iDim] = CurrentCoord[6][iDim] - CurrentCoord[2][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  return Volume;

}

CHEXA1::CHEXA1(void) : CElement() {

}

CHEXA1::CHEXA1(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {

  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;
  nNodes = 8;
  nGaussPoints = 1;

  nDimSq = nDim*nDim;

  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }

  NodalExtrap = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalExtrap[iNode] = new su2double[nGaussPoints];
  }

  NodalStress = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalStress[iNode] = new su2double[6];
  }

  /*--- Initialize structure for current and reference configuration ---*/

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

  GaussCoordP[0][0] = 0.0;  GaussCoordP[0][1] = 0.0;  GaussCoordP[0][1] = 0.0;  GaussWeightP[0] = 8.0;

  Kk_ab = new su2double **[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kk_ab [iNode] = new su2double*[nNodes];
    for (jNode = 0; jNode < nNodes; jNode++) {
      Kk_ab [iNode][jNode] = new su2double[nDimSq];
    }
  }

  /*--- Store the shape functions (they only need to be computed once) ---*/
  su2double Xi, Eta, Zeta, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];

    val_Ni = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0-Zeta);   GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0-Zeta);   GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0-Zeta);   GaussPoint[iGauss]->SetNi(val_Ni,2);
    val_Ni = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0-Zeta);   GaussPoint[iGauss]->SetNi(val_Ni,3);
    val_Ni = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0+Zeta);   GaussPoint[iGauss]->SetNi(val_Ni,4);
    val_Ni = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0+Zeta);   GaussPoint[iGauss]->SetNi(val_Ni,5);
    val_Ni = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0+Zeta);   GaussPoint[iGauss]->SetNi(val_Ni,6);
    val_Ni = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0+Zeta);   GaussPoint[iGauss]->SetNi(val_Ni,7);
  }

  /*--- Shape functions evaluated at the nodes for extrapolation of the stresses at the Gaussian Points ---*/
  /*--- The stress is constant at a HEXA1 element ---*/
  NodalExtrap[0][0] = 1.0;
  NodalExtrap[1][0] = 1.0;
  NodalExtrap[2][0] = 1.0;
  NodalExtrap[3][0] = 1.0;
  NodalExtrap[4][0] = 1.0;
  NodalExtrap[5][0] = 1.0;
  NodalExtrap[6][0] = 1.0;
  NodalExtrap[7][0] = 1.0;

}

CHEXA1::~CHEXA1(void) {

}

void CHEXA1::ComputeGrad_Linear(void){

  su2double Xi, Eta, Zeta;
  su2double Jacobian[3][3], dNiXj[8][3];
  su2double detJac, GradNi_Xj;
  su2double ad[3][3];
  unsigned short iNode, iDim, jDim, iGauss;

  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {

    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];

    /*--- dN/d xi ---*/

    dNiXj[0][0] = -0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[1][0] = 0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[2][0] = 0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[3][0] = -0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[4][0] = -0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[5][0] = 0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[6][0] = 0.125*(1.0+Eta)*(1.0+Zeta);
    dNiXj[7][0] = -0.125*(1.0+Eta)*(1.0+Zeta);

    /*--- dN/d eta ---*/

    dNiXj[0][1] = -0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[1][1] = -0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[2][1] = 0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[3][1] = 0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[4][1] = -0.125*(1.0-Xi)*(1.0+Zeta);
    dNiXj[5][1] = -0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[6][1] = 0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[7][1] = 0.125*(1.0-Xi)*(1.0+Zeta);

    /*--- dN/d mu ---*/

    dNiXj[0][2] = -0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[1][2] = -0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[2][2] = -0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[3][2] = -0.125*(1.0-Xi)*(1.0+Eta);
    dNiXj[4][2] = 0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[5][2] = 0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[6][2] = 0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[7][2] = 0.125*(1.0-Xi)*(1.0+Eta);


    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }

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

    GaussPoint[iGauss]->SetJ_X(detJac);

    /*--- Jacobian inverse (it was already computed as transpose) ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
    }

    /*--- Derivatives with respect to global coordinates ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
      }
    }
  }

}

void CHEXA1::ComputeGrad_NonLinear(void) {

  su2double Xi, Eta, Zeta;
  su2double Jac_Ref[3][3], Jac_Curr[3][3], dNiXj[8][3];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[3][3], ad_Curr[3][3];
  unsigned short iNode, iDim, jDim, iGauss;

  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {

    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];

    /*--- dN/d xi, dN/d eta ---*/
    
    /*--- dN/d xi ---*/
    
    dNiXj[0][0] = -0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[1][0] = 0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[2][0] = 0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[3][0] = -0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[4][0] = -0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[5][0] = 0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[6][0] = 0.125*(1.0+Eta)*(1.0+Zeta);
    dNiXj[7][0] = -0.125*(1.0+Eta)*(1.0+Zeta);

    /*--- dN/d eta ---*/

    dNiXj[0][1] = -0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[1][1] = -0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[2][1] = 0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[3][1] = 0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[4][1] = -0.125*(1.0-Xi)*(1.0+Zeta);
    dNiXj[5][1] = -0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[6][1] = 0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[7][1] = 0.125*(1.0-Xi)*(1.0+Zeta);

    /*--- dN/d mu ---*/

    dNiXj[0][2] = -0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[1][2] = -0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[2][2] = -0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[3][2] = -0.125*(1.0-Xi)*(1.0+Eta);
    dNiXj[4][2] = 0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[5][2] = 0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[6][2] = 0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[7][2] = 0.125*(1.0-Xi)*(1.0+Eta);

    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }

    /*--- Adjoint to Jacobian ---*/

    ad_Ref[0][0] = Jac_Ref[1][1]*Jac_Ref[2][2]-Jac_Ref[1][2]*Jac_Ref[2][1];
    ad_Ref[0][1] = Jac_Ref[0][2]*Jac_Ref[2][1]-Jac_Ref[0][1]*Jac_Ref[2][2];
    ad_Ref[0][2] = Jac_Ref[0][1]*Jac_Ref[1][2]-Jac_Ref[0][2]*Jac_Ref[1][1];
    ad_Ref[1][0] = Jac_Ref[1][2]*Jac_Ref[2][0]-Jac_Ref[1][0]*Jac_Ref[2][2];
    ad_Ref[1][1] = Jac_Ref[0][0]*Jac_Ref[2][2]-Jac_Ref[0][2]*Jac_Ref[2][0];
    ad_Ref[1][2] = Jac_Ref[0][2]*Jac_Ref[1][0]-Jac_Ref[0][0]*Jac_Ref[1][2];
    ad_Ref[2][0] = Jac_Ref[1][0]*Jac_Ref[2][1]-Jac_Ref[1][1]*Jac_Ref[2][0];
    ad_Ref[2][1] = Jac_Ref[0][1]*Jac_Ref[2][0]-Jac_Ref[0][0]*Jac_Ref[2][1];
    ad_Ref[2][2] = Jac_Ref[0][0]*Jac_Ref[1][1]-Jac_Ref[0][1]*Jac_Ref[1][0];

    ad_Curr[0][0] = Jac_Curr[1][1]*Jac_Curr[2][2]-Jac_Curr[1][2]*Jac_Curr[2][1];
    ad_Curr[0][1] = Jac_Curr[0][2]*Jac_Curr[2][1]-Jac_Curr[0][1]*Jac_Curr[2][2];
    ad_Curr[0][2] = Jac_Curr[0][1]*Jac_Curr[1][2]-Jac_Curr[0][2]*Jac_Curr[1][1];
    ad_Curr[1][0] = Jac_Curr[1][2]*Jac_Curr[2][0]-Jac_Curr[1][0]*Jac_Curr[2][2];
    ad_Curr[1][1] = Jac_Curr[0][0]*Jac_Curr[2][2]-Jac_Curr[0][2]*Jac_Curr[2][0];
    ad_Curr[1][2] = Jac_Curr[0][2]*Jac_Curr[1][0]-Jac_Curr[0][0]*Jac_Curr[1][2];
    ad_Curr[2][0] = Jac_Curr[1][0]*Jac_Curr[2][1]-Jac_Curr[1][1]*Jac_Curr[2][0];
    ad_Curr[2][1] = Jac_Curr[0][1]*Jac_Curr[2][0]-Jac_Curr[0][0]*Jac_Curr[2][1];
    ad_Curr[2][2] = Jac_Curr[0][0]*Jac_Curr[1][1]-Jac_Curr[0][1]*Jac_Curr[1][0];


    /*--- Determinant of Jacobian ---*/

    detJac_Ref = Jac_Ref[0][0]*ad_Ref[0][0]+Jac_Ref[0][1]*ad_Ref[1][0]+Jac_Ref[0][2]*ad_Ref[2][0];
    detJac_Curr = Jac_Curr[0][0]*ad_Curr[0][0]+Jac_Curr[0][1]*ad_Curr[1][0]+Jac_Curr[0][2]*ad_Curr[2][0];

    GaussPoint[iGauss]->SetJ_X(detJac_Ref);
    GaussPoint[iGauss]->SetJ_x(detJac_Curr);

    /*--- Jacobian inverse (it was already computed as transpose) ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }

    /*--- Derivatives with respect to global coordinates ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }

}

su2double CHEXA1::ComputeVolume(void){

  unsigned short iDim;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0};
  su2double Volume = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[1][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[2][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[5][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[2][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[7][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[5][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[2][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[3][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[7][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[5][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[7][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[4][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[7][iDim] - RefCoord[2][iDim];
    r2[iDim] = RefCoord[5][iDim] - RefCoord[2][iDim];
    r3[iDim] = RefCoord[6][iDim] - RefCoord[2][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  return Volume;

}

void CHEXA1::ComputeGrad_Pressure(void) {

  su2double Xi, Eta, Zeta;
  su2double Jac_Ref[3][3], Jac_Curr[3][3], dNiXj[8][3];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[3][3], ad_Curr[3][3];
  unsigned short iNode, iDim, jDim, iGauss;

  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {

    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];

    /*--- dN/d xi, dN/d eta ---*/

    /*--- dN/d xi ---*/

    dNiXj[0][0] = -0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[1][0] = 0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[2][0] = 0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[3][0] = -0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[4][0] = -0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[5][0] = 0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[6][0] = 0.125*(1.0+Eta)*(1.0+Zeta);
    dNiXj[7][0] = -0.125*(1.0+Eta)*(1.0+Zeta);

    /*--- dN/d eta ---*/

    dNiXj[0][1] = -0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[1][1] = -0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[2][1] = 0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[3][1] = 0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[4][1] = -0.125*(1.0-Xi)*(1.0+Zeta);
    dNiXj[5][1] = -0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[6][1] = 0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[7][1] = 0.125*(1.0-Xi)*(1.0+Zeta);

    /*--- dN/d mu ---*/

    dNiXj[0][2] = -0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[1][2] = -0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[2][2] = -0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[3][2] = -0.125*(1.0-Xi)*(1.0+Eta);
    dNiXj[4][2] = 0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[5][2] = 0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[6][2] = 0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[7][2] = 0.125*(1.0-Xi)*(1.0+Eta);

    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }

    /*--- Adjoint to Jacobian ---*/

    ad_Ref[0][0] = Jac_Ref[1][1]*Jac_Ref[2][2]-Jac_Ref[1][2]*Jac_Ref[2][1];
    ad_Ref[0][1] = Jac_Ref[0][2]*Jac_Ref[2][1]-Jac_Ref[0][1]*Jac_Ref[2][2];
    ad_Ref[0][2] = Jac_Ref[0][1]*Jac_Ref[1][2]-Jac_Ref[0][2]*Jac_Ref[1][1];
    ad_Ref[1][0] = Jac_Ref[1][2]*Jac_Ref[2][0]-Jac_Ref[1][0]*Jac_Ref[2][2];
    ad_Ref[1][1] = Jac_Ref[0][0]*Jac_Ref[2][2]-Jac_Ref[0][2]*Jac_Ref[2][0];
    ad_Ref[1][2] = Jac_Ref[0][2]*Jac_Ref[1][0]-Jac_Ref[0][0]*Jac_Ref[1][2];
    ad_Ref[2][0] = Jac_Ref[1][0]*Jac_Ref[2][1]-Jac_Ref[1][1]*Jac_Ref[2][0];
    ad_Ref[2][1] = Jac_Ref[0][1]*Jac_Ref[2][0]-Jac_Ref[0][0]*Jac_Ref[2][1];
    ad_Ref[2][2] = Jac_Ref[0][0]*Jac_Ref[1][1]-Jac_Ref[0][1]*Jac_Ref[1][0];

    ad_Curr[0][0] = Jac_Curr[1][1]*Jac_Curr[2][2]-Jac_Curr[1][2]*Jac_Curr[2][1];
    ad_Curr[0][1] = Jac_Curr[0][2]*Jac_Curr[2][1]-Jac_Curr[0][1]*Jac_Curr[2][2];
    ad_Curr[0][2] = Jac_Curr[0][1]*Jac_Curr[1][2]-Jac_Curr[0][2]*Jac_Curr[1][1];
    ad_Curr[1][0] = Jac_Curr[1][2]*Jac_Curr[2][0]-Jac_Curr[1][0]*Jac_Curr[2][2];
    ad_Curr[1][1] = Jac_Curr[0][0]*Jac_Curr[2][2]-Jac_Curr[0][2]*Jac_Curr[2][0];
    ad_Curr[1][2] = Jac_Curr[0][2]*Jac_Curr[1][0]-Jac_Curr[0][0]*Jac_Curr[1][2];
    ad_Curr[2][0] = Jac_Curr[1][0]*Jac_Curr[2][1]-Jac_Curr[1][1]*Jac_Curr[2][0];
    ad_Curr[2][1] = Jac_Curr[0][1]*Jac_Curr[2][0]-Jac_Curr[0][0]*Jac_Curr[2][1];
    ad_Curr[2][2] = Jac_Curr[0][0]*Jac_Curr[1][1]-Jac_Curr[0][1]*Jac_Curr[1][0];


    /*--- Determinant of Jacobian ---*/

    detJac_Ref = Jac_Ref[0][0]*ad_Ref[0][0]+Jac_Ref[0][1]*ad_Ref[1][0]+Jac_Ref[0][2]*ad_Ref[2][0];
    detJac_Curr = Jac_Curr[0][0]*ad_Curr[0][0]+Jac_Curr[0][1]*ad_Curr[1][0]+Jac_Curr[0][2]*ad_Curr[2][0];

    GaussPoint[iGauss]->SetJ_X(detJac_Ref);
    GaussPoint[iGauss]->SetJ_x(detJac_Curr);

    /*--- Jacobian inverse (it was already computed as transpose) ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }

    /*--- Derivatives with respect to global coordinates ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }

}


CPYRAM5::CPYRAM5(void) : CElement() {

}

CPYRAM5::CPYRAM5(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {

  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;

  nNodes = 5;
  nGaussPoints = 5;

  nDimSq = nDim*nDim;

  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }

  /* These structures are only needed for non-linear elasticity in structural analysiss
   * These elements are only used for mesh deformation
   */

  NodalExtrap  = NULL;
  NodalStress  = NULL;
  CurrentCoord = NULL;
  Mab          = NULL;
  Ks_ab        = NULL;
  Kt_a         = NULL;
  FDL_a        = NULL;

  /* Initialize structure only for reference configuration
   * These elements are only used for mesh deformation
   */

  RefCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++){
    RefCoord [iNode] = new su2double[nDim];
  }

  GaussWeight = new su2double [nGaussPoints];

  GaussCoord = new su2double*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++){
    GaussCoord [iGauss] = new su2double[nDim];
  }

  GaussCoord[0][0] = 0.5;   GaussCoord[0][1] = 0.0;   GaussCoord[0][2] = 0.1531754163448146;  GaussWeight[0] = 0.133333333333333;
  GaussCoord[1][0] = 0.0;   GaussCoord[1][1] = 0.5;   GaussCoord[1][2] = 0.1531754163448146;  GaussWeight[1] = 0.133333333333333;
  GaussCoord[2][0] = -0.5;  GaussCoord[2][1] = 0.0;   GaussCoord[2][2] = 0.1531754163448146;  GaussWeight[2] = 0.133333333333333;
  GaussCoord[3][0] = 0.0;   GaussCoord[3][1] = -0.5;  GaussCoord[3][2] = 0.1531754163448146;  GaussWeight[3] = 0.133333333333333;
  GaussCoord[4][0] = 0.0;   GaussCoord[4][1] = 0.0;   GaussCoord[4][2] = 0.6372983346207416;  GaussWeight[4] = 0.133333333333333;

  Kab = new su2double **[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++){
    Kab [iNode] = new su2double*[nNodes];
    for (jNode = 0; jNode < nNodes; jNode++){
      Kab [iNode][jNode] = new su2double[nDimSq];
    }
  }

  /*--- Store the shape functions (they only need to be computed once) ---*/
  su2double Xi, Eta, Zeta, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++){
      Xi = GaussCoord[iGauss][0];
      Eta = GaussCoord[iGauss][1];
      Zeta = GaussCoord[iGauss][2];

      val_Ni = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0-Zeta); GaussPoint[iGauss]->SetNi(val_Ni,0);
      val_Ni = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0-Zeta); GaussPoint[iGauss]->SetNi(val_Ni,1);
      val_Ni = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0-Zeta); GaussPoint[iGauss]->SetNi(val_Ni,2);
      val_Ni = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0-Zeta); GaussPoint[iGauss]->SetNi(val_Ni,3);
      val_Ni = 0.5*(1.0+Zeta);                      GaussPoint[iGauss]->SetNi(val_Ni,4);

  }

}

CPYRAM5::~CPYRAM5(void) {

}

void CPYRAM5::ComputeGrad_Linear(void){

    su2double Xi, Eta, Zeta;
    su2double Jacobian[3][3], dNiXj[8][3];
    su2double detJac, GradNi_Xj;
    su2double ad[3][3];
    unsigned short iNode, iDim, jDim, iGauss;

    for (iGauss = 0; iGauss < nGaussPoints; iGauss++){

      Xi = GaussCoord[iGauss][0];
      Eta = GaussCoord[iGauss][1];
      Zeta = GaussCoord[iGauss][2];

      /*--- dN/d xi ---*/

      dNiXj[0][0] = -0.125*(1.0-Eta)*(1.0-Zeta);
      dNiXj[1][0] = 0.125*(1.0-Eta)*(1.0-Zeta);
      dNiXj[2][0] = 0.125*(1.0+Eta)*(1.0-Zeta);
      dNiXj[3][0] = -0.125*(1.0+Eta)*(1.0-Zeta);
      dNiXj[4][0] = 0.0;

      /*--- dN/d eta ---*/

      dNiXj[0][1] = -0.125*(1.0-Xi)*(1.0-Zeta);
      dNiXj[1][1] = -0.125*(1.0+Xi)*(1.0-Zeta);
      dNiXj[2][1] = 0.125*(1.0+Xi)*(1.0-Zeta);
      dNiXj[3][1] = 0.125*(1.0-Xi)*(1.0-Zeta);
      dNiXj[4][1] = 0.0;

      /*--- dN/d mu ---*/

      dNiXj[0][2] = -0.125*(1.0-Xi)*(1.0-Eta);
      dNiXj[1][2] = -0.125*(1.0+Xi)*(1.0-Eta);
      dNiXj[2][2] = -0.125*(1.0+Xi)*(1.0+Eta);
      dNiXj[3][2] = -0.125*(1.0-Xi)*(1.0+Eta);
      dNiXj[4][2] = 0.5;

      /*--- Jacobian transformation ---*/
      /*--- This does dX/dXi transpose ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
      }

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

      GaussPoint[iGauss]->SetJ_X(detJac);

      /*--- Jacobian inverse (it was already computed as transpose) ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
      }

      /*--- Derivatives with respect to global coordinates ---*/

      for (iNode = 0; iNode < nNodes; iNode++) {
        for (iDim = 0; iDim < nDim; iDim++){
          GradNi_Xj = 0.0;
          for (jDim = 0; jDim < nDim; jDim++){
            GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
          }
          GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
        }
      }
    }

}

su2double CPYRAM5::ComputeVolume(void){

  unsigned short iDim;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0};
  su2double Volume = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[1][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[2][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[4][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[2][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[3][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[4][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  return Volume;

}


CPRISM6::CPRISM6(void) : CElement() {

}

CPRISM6::CPRISM6(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {

  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;

  nNodes = 6;
  nGaussPoints = 6;

  nDimSq = nDim*nDim;

  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }

  /* These structures are only needed for non-linear elasticity in structural analysiss
   * These elements are only used for mesh deformation
   */

  NodalExtrap  = NULL;
  NodalStress  = NULL;
  CurrentCoord = NULL;
  Mab          = NULL;
  Ks_ab        = NULL;
  Kt_a         = NULL;
  FDL_a        = NULL;

  /* Initialize structure only for reference configuration
   * These elements are only used for mesh deformation
   */

  RefCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++){
    RefCoord [iNode] = new su2double[nDim];
  }

  GaussWeight = new su2double [nGaussPoints];

  GaussCoord = new su2double*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++){
    GaussCoord [iGauss] = new su2double[nDim];
  }

  GaussCoord[0][0] = 0.5;                 GaussCoord[0][1] = 0.5;                 GaussCoord[0][2] = -0.577350269189626;  GaussWeight[0] = 0.166666666666666;
  GaussCoord[1][0] = -0.577350269189626;  GaussCoord[1][1] = 0.0;                 GaussCoord[1][2] = 0.5;                 GaussWeight[1] = 0.166666666666666;
  GaussCoord[2][0] = 0.5;                 GaussCoord[2][1] = -0.577350269189626;  GaussCoord[2][2] = 0.0;                 GaussWeight[2] = 0.166666666666666;
  GaussCoord[3][0] = 0.5;                 GaussCoord[3][1] = 0.5;                 GaussCoord[3][2] = 0.577350269189626;   GaussWeight[3] = 0.166666666666666;
  GaussCoord[4][0] = 0.577350269189626;   GaussCoord[4][1] = 0.0;                 GaussCoord[4][2] = 0.5;                 GaussWeight[4] = 0.166666666666666;
  GaussCoord[5][0] = 0.5;                 GaussCoord[5][1] = 0.577350269189626;   GaussCoord[5][2] = 0.0;                 GaussWeight[5] = 0.166666666666666;

  Kab = new su2double **[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++){
    Kab [iNode] = new su2double*[nNodes];
    for (jNode = 0; jNode < nNodes; jNode++){
      Kab [iNode][jNode] = new su2double[nDimSq];
    }
  }

  /*--- Store the shape functions (they only need to be computed once) ---*/
  su2double Xi, Eta, Zeta, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++){
      Xi = GaussCoord[iGauss][0];
      Eta = GaussCoord[iGauss][1];
      Zeta = GaussCoord[iGauss][2];

      val_Ni = 0.5*Eta*(1.0-Xi);              GaussPoint[iGauss]->SetNi(val_Ni,0);
      val_Ni = 0.5*Zeta*(1.0-Xi);             GaussPoint[iGauss]->SetNi(val_Ni,1);
      val_Ni = 0.5*(1.0-Eta-Zeta)*(1.0-Xi);   GaussPoint[iGauss]->SetNi(val_Ni,2);
      val_Ni = 0.5*Eta*(Xi+1.0);              GaussPoint[iGauss]->SetNi(val_Ni,3);
      val_Ni = 0.5*Zeta*(Xi+1.0);             GaussPoint[iGauss]->SetNi(val_Ni,4);
      val_Ni = 0.5*(1.0-Eta-Zeta)*(Xi+1.0);   GaussPoint[iGauss]->SetNi(val_Ni,5);

  }

}

CPRISM6::~CPRISM6(void) {

}

void CPRISM6::ComputeGrad_Linear(void){

    su2double Xi, Eta, Zeta;
    su2double Jacobian[3][3], dNiXj[8][3];
    su2double detJac, GradNi_Xj;
    su2double ad[3][3];
    unsigned short iNode, iDim, jDim, iGauss;

    for (iGauss = 0; iGauss < nGaussPoints; iGauss++){

      Xi = GaussCoord[iGauss][0];
      Eta = GaussCoord[iGauss][1];
      Zeta = GaussCoord[iGauss][2];

      /*--- dN/d xi ---*/

      dNiXj[0][0] = -0.5*Eta;
      dNiXj[1][0] = -0.5*Zeta;
      dNiXj[2][0] = -0.5*(1.0-Eta-Zeta);
      dNiXj[3][0] = 0.5*Eta;
      dNiXj[4][0] = 0.5*Zeta;
      dNiXj[5][0] = 0.5*(1.0-Eta-Zeta);

      /*--- dN/d eta ---*/

      dNiXj[0][1] = 0.5*(1.0-Xi);
      dNiXj[1][1] = 0.0;
      dNiXj[2][1] = -0.5*(1.0-Xi);
      dNiXj[3][1] = 0.5*(Xi+1.0);
      dNiXj[4][1] = 0.0;
      dNiXj[5][1] = -0.5*(Xi+1.0);

      /*--- dN/d mu ---*/

      dNiXj[0][2] = 0.0;
      dNiXj[1][2] = 0.5*(1.0-Xi);
      dNiXj[2][2] = -0.5*(1.0-Xi);
      dNiXj[3][2] = 0.0;
      dNiXj[4][2] = 0.5*(Xi+1.0);
      dNiXj[5][2] = -0.5*(Xi+1.0);

      /*--- Jacobian transformation ---*/
      /*--- This does dX/dXi transpose ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
      }

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

      GaussPoint[iGauss]->SetJ_X(detJac);

      /*--- Jacobian inverse (it was already computed as transpose) ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
      }

      /*--- Derivatives with respect to global coordinates ---*/

      for (iNode = 0; iNode < nNodes; iNode++) {
        for (iDim = 0; iDim < nDim; iDim++){
          GradNi_Xj = 0.0;
          for (jDim = 0; jDim < nDim; jDim++){
            GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
          }
          GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
        }
      }
    }

}

su2double CPRISM6::ComputeVolume(void){

  unsigned short iDim;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0};
  su2double Volume = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[2][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[1][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[5][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[5][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[1][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[4][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = RefCoord[5][iDim] - RefCoord[0][iDim];
    r2[iDim] = RefCoord[4][iDim] - RefCoord[0][iDim];
    r3[iDim] = RefCoord[3][iDim] - RefCoord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  return Volume;

}



