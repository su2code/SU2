/*!
 * \file numerics_direct_elasticity.cpp
 * \brief This file contains all the convective term discretization.
 * \author F. Palacios
 * \version 3.2.8.1 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (fpalacios@stanford.edu).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
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

CGalerkin_FEA::CGalerkin_FEA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	E = config->GetElasticyMod();
	Nu = config->GetPoissonRatio();
	Mu = E / (2.0*(1.0 + Nu));
	Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));
  
}

CGalerkin_FEA::~CGalerkin_FEA(void) { }

double CGalerkin_FEA::ShapeFunc_Triangle(double Xi, double Eta, double CoordCorners[8][3], double DShapeFunction[8][4]) {
  
  int i, j, k;
  double c0, c1, xsj;
  double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = 1-Xi-Eta;
  DShapeFunction[1][3] = Xi;
  DShapeFunction[2][3] = Eta;
  
  /*--- dN/d xi, dN/d eta, dN/d mu ---*/
  
  DShapeFunction[0][0] = -1.0;  DShapeFunction[0][1] = -1.0;
  DShapeFunction[1][0] = 1;     DShapeFunction[1][1] = 0.0;
  DShapeFunction[2][0] = 0;     DShapeFunction[2][1] = 1;
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 3; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to Jacobian ---*/
  
  ad[0][0] = xs[1][1];
  ad[0][1] = -xs[0][1];
  ad[1][0] = -xs[1][0];
  ad[1][1] = xs[0][0];
  
  /*--- Determinant of Jacobian ---*/
  
  xsj = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 3; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]; // dN/dy
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
  }
  
  return xsj;
  
}

double CGalerkin_FEA::ShapeFunc_Rectangle(double Xi, double Eta, double CoordCorners[8][3], double DShapeFunction[8][4]) {
  
  int i, j, k;
  double c0, c1, xsj;
  double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = 0.25*(1.0-Xi)*(1.0-Eta);
  DShapeFunction[1][3] = 0.25*(1.0+Xi)*(1.0-Eta);
  DShapeFunction[2][3] = 0.25*(1.0+Xi)*(1.0+Eta);
  DShapeFunction[3][3] = 0.25*(1.0-Xi)*(1.0+Eta);
  
  /*--- dN/d xi, dN/d eta, dN/d mu ---*/
  
  DShapeFunction[0][0] = -0.25*(1.0-Eta); DShapeFunction[0][1] = -0.25*(1.0-Xi);
  DShapeFunction[1][0] =  0.25*(1.0-Eta); DShapeFunction[1][1] = -0.25*(1.0+Xi);
  DShapeFunction[2][0] =  0.25*(1.0+Eta); DShapeFunction[2][1] =  0.25*(1.0+Xi);
  DShapeFunction[3][0] = -0.25*(1.0+Eta); DShapeFunction[3][1] =  0.25*(1.0-Xi);
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 4; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to Jacobian ---*/
  
  ad[0][0] = xs[1][1];
  ad[0][1] = -xs[0][1];
  ad[1][0] = -xs[1][0];
  ad[1][1] = xs[0][0];
  
  /*--- Determinant of Jacobian ---*/
  
  xsj = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 4; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]; // dN/dy
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
  }
  
  return xsj;
  
}

double CGalerkin_FEA::ShapeFunc_Tetra(double Xi, double Eta, double Zeta, double CoordCorners[8][3], double DShapeFunction[8][4]) {
  
  int i, j, k;
  double c0, c1, c2, xsj;
  double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = Xi;
  DShapeFunction[1][3] = Eta;
  DShapeFunction[2][3] = Mu;
  DShapeFunction[3][3] = 1.0 - Xi - Eta - Mu;
  
  /*--- dN/d xi, dN/d eta, dN/d mu ---*/
  
  DShapeFunction[0][0] = 1.0;   DShapeFunction[0][1] = 0.0;   DShapeFunction[0][2] = 0.0;
  DShapeFunction[1][0] = 0.0;   DShapeFunction[1][1] = 1.0;   DShapeFunction[1][2] = 0.0;
  DShapeFunction[2][0] = 0.0;   DShapeFunction[2][1] = 0.0;   DShapeFunction[2][2] = 1.0;
  DShapeFunction[3][0] = -1.0;  DShapeFunction[3][1] = -1.0;  DShapeFunction[3][2] = -1.0;
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 4; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to Jacobian ---*/
  
  ad[0][0] = xs[1][1]*xs[2][2]-xs[1][2]*xs[2][1];
  ad[0][1] = xs[0][2]*xs[2][1]-xs[0][1]*xs[2][2];
  ad[0][2] = xs[0][1]*xs[1][2]-xs[0][2]*xs[1][1];
  ad[1][0] = xs[1][2]*xs[2][0]-xs[1][0]*xs[2][2];
  ad[1][1] = xs[0][0]*xs[2][2]-xs[0][2]*xs[2][0];
  ad[1][2] = xs[0][2]*xs[1][0]-xs[0][0]*xs[1][2];
  ad[2][0] = xs[1][0]*xs[2][1]-xs[1][1]*xs[2][0];
  ad[2][1] = xs[0][1]*xs[2][0]-xs[0][0]*xs[2][1];
  ad[2][2] = xs[0][0]*xs[1][1]-xs[0][1]*xs[1][0];
  
  /*--- Determinant of Jacobian ---*/
  
  xsj = xs[0][0]*ad[0][0]+xs[0][1]*ad[1][0]+xs[0][2]*ad[2][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 4; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]+xs[0][2]*DShapeFunction[k][2]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]+xs[1][2]*DShapeFunction[k][2]; // dN/dy
    c2 = xs[2][0]*DShapeFunction[k][0]+xs[2][1]*DShapeFunction[k][1]+xs[2][2]*DShapeFunction[k][2]; // dN/dz
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d mu
  }
  
  return xsj;
  
}

double CGalerkin_FEA::ShapeFunc_Pyram(double Xi, double Eta, double Zeta, double CoordCorners[8][3], double DShapeFunction[8][4]) {
  
  int i, j, k;
  double c0, c1, c2, xsj;
  double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  double Den = 4.0*(1.0 - Mu);
  
  DShapeFunction[0][3] = (-Xi+Eta+Mu-1.0)*(-Xi-Eta+Mu-1.0)/Den;
  DShapeFunction[1][3] = (-Xi-Eta+Mu-1.0)*(Xi-Eta+Mu-1.0)/Den;
  DShapeFunction[2][3] = (Xi+Eta+Mu-1.0)*(Xi-Eta+Mu-1.0)/Den;
  DShapeFunction[3][3] = (Xi+Eta+Mu-1.0)*(-Xi+Eta+Mu-1.0)/Den;
  DShapeFunction[4][3] = Mu;
  
  /*--- dN/d xi, dN/d eta, dN/d mu ---*/
  
  DShapeFunction[0][0] = 0.5 + (0.5*Xi)/(1.0 - Mu);
  DShapeFunction[0][1] = (0.5*Eta)/(-1.0 + Mu);
  DShapeFunction[0][2] = (-0.25 - 0.25*Eta*Eta + (0.5 - 0.25*Mu)*Mu + 0.25*Xi*Xi)/((-1.0 + Mu)*(-1.0 + Mu));
  
  DShapeFunction[1][0] = (0.5*Xi)/(-1.0 + Mu);
  DShapeFunction[1][1] = (-0.5 - 0.5*Eta + 0.5*Mu)/(-1.0 + Mu);
  DShapeFunction[1][2] = (-0.25 + 0.25*Eta*Eta + (0.5 - 0.25*Mu)*Mu - 0.25*Xi*Xi)/((-1.0 + Mu)*(-1.0 + Mu));
  
  DShapeFunction[2][0] = -0.5 + (0.5*Xi)/(1.0 - 1.0*Mu);
  DShapeFunction[2][1] = (0.5*Eta)/(-1.0 + Mu);
  DShapeFunction[2][2] = (-0.25 - 0.25*Eta*Eta + (0.5 - 0.25*Mu)*Mu + 0.25*Xi*Xi)/((-1.0 + Mu)*(-1.0 + Mu));
  
  DShapeFunction[3][0] = (0.5*Xi)/(-1.0 + Mu);
  DShapeFunction[3][1] = (0.5 - 0.5*Eta - 0.5*Mu)/(-1.0 + Mu);
  DShapeFunction[3][2] = (-0.25 + 0.25*Eta*Eta + (0.5 - 0.25*Mu)*Mu - 0.25*Xi*Xi)/((-1.0 + Mu)*(-1.0 + Mu));
  
  DShapeFunction[4][0] = 0.0;
  DShapeFunction[4][1] = 0.0;
  DShapeFunction[4][2] = 1.0;
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 5; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to Jacobian ---*/
  
  ad[0][0] = xs[1][1]*xs[2][2]-xs[1][2]*xs[2][1];
  ad[0][1] = xs[0][2]*xs[2][1]-xs[0][1]*xs[2][2];
  ad[0][2] = xs[0][1]*xs[1][2]-xs[0][2]*xs[1][1];
  ad[1][0] = xs[1][2]*xs[2][0]-xs[1][0]*xs[2][2];
  ad[1][1] = xs[0][0]*xs[2][2]-xs[0][2]*xs[2][0];
  ad[1][2] = xs[0][2]*xs[1][0]-xs[0][0]*xs[1][2];
  ad[2][0] = xs[1][0]*xs[2][1]-xs[1][1]*xs[2][0];
  ad[2][1] = xs[0][1]*xs[2][0]-xs[0][0]*xs[2][1];
  ad[2][2] = xs[0][0]*xs[1][1]-xs[0][1]*xs[1][0];
  
  /*--- Determinant of Jacobian ---*/
  
  xsj = xs[0][0]*ad[0][0]+xs[0][1]*ad[1][0]+xs[0][2]*ad[2][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 5; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]+xs[0][2]*DShapeFunction[k][2]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]+xs[1][2]*DShapeFunction[k][2]; // dN/dy
    c2 = xs[2][0]*DShapeFunction[k][0]+xs[2][1]*DShapeFunction[k][1]+xs[2][2]*DShapeFunction[k][2]; // dN/dz
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d mu
  }
  
  return xsj;
  
}

double CGalerkin_FEA::ShapeFunc_Wedge(double Xi, double Eta, double Zeta, double CoordCorners[8][3], double DShapeFunction[8][4]) {
  
  int i, j, k;
  double c0, c1, c2, xsj;
  double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = 0.5*Eta*(1.0-Xi);
  DShapeFunction[1][3] = 0.5*Mu*(1.0-Xi);;
  DShapeFunction[2][3] = 0.5*(1.0-Eta-Mu)*(1.0-Xi);
  DShapeFunction[3][3] = 0.5*Eta*(Xi+1.0);
  DShapeFunction[4][3] = 0.5*Mu*(Xi+1.0);
  DShapeFunction[5][3] = 0.5*(1.0-Eta-Mu)*(Xi+1.0);
  
  /*--- dN/d Xi, dN/d Eta, dN/d Mu ---*/
  
  DShapeFunction[0][0] = -0.5*Eta;            DShapeFunction[0][1] = 0.5*(1.0-Xi);      DShapeFunction[0][2] = 0.0;
  DShapeFunction[1][0] = -0.5*Mu;             DShapeFunction[1][1] = 0.0;               DShapeFunction[1][2] = 0.5*(1.0-Xi);
  DShapeFunction[2][0] = -0.5*(1.0-Eta-Mu);   DShapeFunction[2][1] = -0.5*(1.0-Xi);     DShapeFunction[2][2] = -0.5*(1.0-Xi);
  DShapeFunction[3][0] = 0.5*Eta;             DShapeFunction[3][1] = 0.5*(Xi+1.0);      DShapeFunction[3][2] = 0.0;
  DShapeFunction[4][0] = 0.5*Mu;              DShapeFunction[4][1] = 0.0;               DShapeFunction[4][2] = 0.5*(Xi+1.0);
  DShapeFunction[5][0] = 0.5*(1.0-Eta-Mu);    DShapeFunction[5][1] = -0.5*(Xi+1.0);     DShapeFunction[5][2] = -0.5*(Xi+1.0);
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 6; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to Jacobian ---*/
  
  ad[0][0] = xs[1][1]*xs[2][2]-xs[1][2]*xs[2][1];
  ad[0][1] = xs[0][2]*xs[2][1]-xs[0][1]*xs[2][2];
  ad[0][2] = xs[0][1]*xs[1][2]-xs[0][2]*xs[1][1];
  ad[1][0] = xs[1][2]*xs[2][0]-xs[1][0]*xs[2][2];
  ad[1][1] = xs[0][0]*xs[2][2]-xs[0][2]*xs[2][0];
  ad[1][2] = xs[0][2]*xs[1][0]-xs[0][0]*xs[1][2];
  ad[2][0] = xs[1][0]*xs[2][1]-xs[1][1]*xs[2][0];
  ad[2][1] = xs[0][1]*xs[2][0]-xs[0][0]*xs[2][1];
  ad[2][2] = xs[0][0]*xs[1][1]-xs[0][1]*xs[1][0];
  
  /*--- Determinant of Jacobian ---*/
  
  xsj = xs[0][0]*ad[0][0]+xs[0][1]*ad[1][0]+xs[0][2]*ad[2][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 6; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]+xs[0][2]*DShapeFunction[k][2]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]+xs[1][2]*DShapeFunction[k][2]; // dN/dy
    c2 = xs[2][0]*DShapeFunction[k][0]+xs[2][1]*DShapeFunction[k][1]+xs[2][2]*DShapeFunction[k][2]; // dN/dz
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d mu
  }
  
  return xsj;
  
}

double CGalerkin_FEA::ShapeFunc_Hexa(double Xi, double Eta, double Zeta, double CoordCorners[8][3], double DShapeFunction[8][4]) {
  
  int i, j, k;
  double a0, a1, a2, c0, c1, c2, xsj;
  double ss[3], xs[3][3], ad[3][3];
  double s0[8] = {-0.5, 0.5, 0.5,-0.5,-0.5, 0.5,0.5,-0.5};
  double s1[8] = {-0.5,-0.5, 0.5, 0.5,-0.5,-0.5,0.5, 0.5};
  double s2[8] = {-0.5,-0.5,-0.5,-0.5, 0.5, 0.5,0.5, 0.5};
  
  ss[0] = Xi;
  ss[1] = Eta;
  ss[2] = Mu;
  
  /*--- Shape functions ---*/
  
  for (i = 0; i < 8; i++) {
    a0 = 0.5+s0[i]*ss[0]; // shape function in xi-direction
    a1 = 0.5+s1[i]*ss[1]; // shape function in eta-direction
    a2 = 0.5+s2[i]*ss[2]; // shape function in mu-direction
    DShapeFunction[i][0] = s0[i]*a1*a2; // dN/d xi
    DShapeFunction[i][1] = s1[i]*a0*a2; // dN/d eta
    DShapeFunction[i][2] = s2[i]*a0*a1; // dN/d mu
    DShapeFunction[i][3] = a0*a1*a2; // actual shape function N
  }
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 8; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to Jacobian ---*/
  
  ad[0][0] = xs[1][1]*xs[2][2]-xs[1][2]*xs[2][1];
  ad[0][1] = xs[0][2]*xs[2][1]-xs[0][1]*xs[2][2];
  ad[0][2] = xs[0][1]*xs[1][2]-xs[0][2]*xs[1][1];
  ad[1][0] = xs[1][2]*xs[2][0]-xs[1][0]*xs[2][2];
  ad[1][1] = xs[0][0]*xs[2][2]-xs[0][2]*xs[2][0];
  ad[1][2] = xs[0][2]*xs[1][0]-xs[0][0]*xs[1][2];
  ad[2][0] = xs[1][0]*xs[2][1]-xs[1][1]*xs[2][0];
  ad[2][1] = xs[0][1]*xs[2][0]-xs[0][0]*xs[2][1];
  ad[2][2] = xs[0][0]*xs[1][1]-xs[0][1]*xs[1][0];
  
  /*--- Determinant of Jacobian ---*/
  
  xsj = xs[0][0]*ad[0][0]+xs[0][1]*ad[1][0]+xs[0][2]*ad[2][0];
  
  /*--- Jacobian inverse ---*/
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 8; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]+xs[0][2]*DShapeFunction[k][2]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]+xs[1][2]*DShapeFunction[k][2]; // dN/dy
    c2 = xs[2][0]*DShapeFunction[k][0]+xs[2][1]*DShapeFunction[k][1]+xs[2][2]*DShapeFunction[k][2]; // dN/dz
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d mu
  }
  
  return xsj;
  
}

void CGalerkin_FEA::SetFEA_StiffMatrix2D(double **StiffMatrix_Elem, double CoordCorners[8][3], unsigned short nNodes) {
  
  double B_Matrix[3][8], D_Matrix[3][3], Aux_Matrix[8][3];
  double Xi = 0.0, Eta = 0.0, Det = 0.0;
  unsigned short iNode, iVar, jVar, kVar, iGauss, nGauss = 0;
  double DShapeFunction[8][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
  double Location[4][3], Weight[4];
  unsigned short nVar = 2;
  
  /*--- Integration formulae from "Shape functions and points of
   integration of the Résumé" by Josselin DELMAS (2013) ---*/

  for (iVar = 0; iVar < nNodes*nVar; iVar++) {
    for (jVar = 0; jVar < nNodes*nVar; jVar++) {
      StiffMatrix_Elem[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Triangle. Nodes of numerical integration at 1 point (order 1). ---*/
  
  if (nNodes == 3) {
    nGauss = 1;
    Location[0][0] = 0.333333333333333;  Location[0][1] = 0.333333333333333;  Weight[0] = 0.5;
  }
  
  /*--- Rectangle. Nodes of numerical integration at 4 points (order 2). ---*/
  
  if (nNodes == 4) {
    nGauss = 4;
    Location[0][0] = -0.577350269189626;  Location[0][1] = -0.577350269189626;  Weight[0] = 1.0;
    Location[1][0] = 0.577350269189626;   Location[1][1] = -0.577350269189626;  Weight[1] = 1.0;
    Location[2][0] = 0.577350269189626;   Location[2][1] = 0.577350269189626;   Weight[2] = 1.0;
    Location[3][0] = -0.577350269189626;  Location[3][1] = 0.577350269189626;   Weight[3] = 1.0;
  }
  
  for (iGauss = 0; iGauss < nGauss; iGauss++) {
    
    Xi = Location[iGauss][0]; Eta = Location[iGauss][1];
    
    if (nNodes == 3) Det = ShapeFunc_Triangle(Xi, Eta, CoordCorners, DShapeFunction);
    if (nNodes == 4) Det = ShapeFunc_Rectangle(Xi, Eta, CoordCorners, DShapeFunction);
    
    /*--- Compute the B Matrix ---*/
    
    for (iVar = 0; iVar < 3; iVar++)
      for (jVar = 0; jVar < nNodes*nVar; jVar++)
        B_Matrix[iVar][jVar] = 0.0;
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      B_Matrix[0][0+iNode*nVar] = DShapeFunction[iNode][0];
      B_Matrix[1][1+iNode*nVar] = DShapeFunction[iNode][1];
      
      B_Matrix[2][0+iNode*nVar] = DShapeFunction[iNode][1];
      B_Matrix[2][1+iNode*nVar] = DShapeFunction[iNode][0];
    }
    
    /*--- Compute the D Matrix (for plane strain and 3-D)---*/
    
    D_Matrix[0][0] = Lambda + 2.0*Mu;		D_Matrix[0][1] = Lambda;            D_Matrix[0][2] = 0.0;
    D_Matrix[1][0] = Lambda;            D_Matrix[1][1] = Lambda + 2.0*Mu;   D_Matrix[1][2] = 0.0;
    D_Matrix[2][0] = 0.0;               D_Matrix[2][1] = 0.0;               D_Matrix[2][2] = Mu;
    
    /*--- Compute the BT.D Matrix ---*/
    
    for (iVar = 0; iVar < nNodes*nVar; iVar++) {
      for (jVar = 0; jVar < 3; jVar++) {
        Aux_Matrix[iVar][jVar] = 0.0;
        for (kVar = 0; kVar < 3; kVar++)
          Aux_Matrix[iVar][jVar] += B_Matrix[kVar][iVar]*D_Matrix[kVar][jVar];
      }
    }
    
    /*--- Compute the BT.D.B Matrix (stiffness matrix), and add to the original
     matrix using Gauss integration ---*/
    
    for (iVar = 0; iVar < nNodes*nVar; iVar++) {
      for (jVar = 0; jVar < nNodes*nVar; jVar++) {
        for (kVar = 0; kVar < 3; kVar++) {
          StiffMatrix_Elem[iVar][jVar] += Weight[iGauss] * Aux_Matrix[iVar][kVar]*B_Matrix[kVar][jVar] * Det;
        }
      }
    }
    
  }
  
}

void CGalerkin_FEA::SetFEA_StiffMatrix3D(double **StiffMatrix_Elem, double CoordCorners[8][3], unsigned short nNodes) {
  
  double B_Matrix[6][24], D_Matrix[6][6], Aux_Matrix[24][6];
  double Xi = 0.0, Eta = 0.0, Det = 0.0;
  unsigned short iNode, iVar, jVar, kVar, iGauss, nGauss = 0;
  double DShapeFunction[8][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
  double Location[8][3], Weight[8];
  
  unsigned short nVar = 3;
  
  /*--- Integration formulae from "Shape functions and points of
   integration of the Résumé" by Josselin Delmas (2013) ---*/
  
  for (iVar = 0; iVar < nNodes*nVar; iVar++) {
    for (jVar = 0; jVar < nNodes*nVar; jVar++) {
      StiffMatrix_Elem[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Tetrahedrons. Nodes of numerical integration at 1 point (order 1). ---*/
  
  if (nNodes == 4) {
    nGauss = 1;
    Location[0][0] = 0.25;  Location[0][1] = 0.25;  Location[0][2] = 0.25;  Weight[0] = 0.166666666666666;
  }
  
  /*--- Pyramids. Nodes numerical integration at 5 points. ---*/
  
  if (nNodes == 5) {
    nGauss = 5;
    Location[0][0] = 0.5;   Location[0][1] = 0.0;   Location[0][2] = 0.1531754163448146;  Weight[0] = 0.133333333333333;
    Location[1][0] = 0.0;   Location[1][1] = 0.5;   Location[1][2] = 0.1531754163448146;  Weight[1] = 0.133333333333333;
    Location[2][0] = -0.5;  Location[2][1] = 0.0;   Location[2][2] = 0.1531754163448146;  Weight[2] = 0.133333333333333;
    Location[3][0] = 0.0;   Location[3][1] = -0.5;  Location[3][2] = 0.1531754163448146;  Weight[3] = 0.133333333333333;
    Location[4][0] = 0.0;   Location[4][1] = 0.0;   Location[4][2] = 0.6372983346207416;  Weight[4] = 0.133333333333333;
  }
  
  /*--- Wedge. Nodes of numerical integration at 6 points (order 3 in Xi, order 2 in Eta and Mu ). ---*/
  
  if (nNodes == 6) {
    nGauss = 6;
    Location[0][0] = 0.5;                 Location[0][1] = 0.5;                 Location[0][2] = -0.577350269189626;  Weight[0] = 0.166666666666666;
    Location[1][0] = -0.577350269189626;  Location[1][1] = 0.0;                 Location[1][2] = 0.5;                 Weight[1] = 0.166666666666666;
    Location[2][0] = 0.5;                 Location[2][1] = -0.577350269189626;  Location[2][2] = 0.0;                 Weight[2] = 0.166666666666666;
    Location[3][0] = 0.5;                 Location[3][1] = 0.5;                 Location[3][2] = 0.577350269189626;   Weight[3] = 0.166666666666666;
    Location[4][0] = 0.577350269189626;   Location[4][1] = 0.0;                 Location[4][2] = 0.5;                 Weight[4] = 0.166666666666666;
    Location[5][0] = 0.5;                 Location[5][1] = 0.577350269189626;   Location[5][2] = 0.0;                 Weight[5] = 0.166666666666666;
  }
  
  /*--- Hexahedrons. Nodes of numerical integration at 6 points (order 3). ---*/
  
  if (nNodes == 8) {
    nGauss = 8;
    Location[0][0] = -0.577350269189626;  Location[0][1] = -0.577350269189626;  Location[0][2] = -0.577350269189626;  Weight[0] = 1.0;
    Location[1][0] = -0.577350269189626;  Location[1][1] = -0.577350269189626;  Location[1][2] = 0.577350269189626;   Weight[1] = 1.0;
    Location[2][0] = -0.577350269189626;  Location[2][1] = 0.577350269189626;   Location[2][2] = -0.577350269189626;  Weight[2] = 1.0;
    Location[3][0] = -0.577350269189626;  Location[3][1] = 0.577350269189626;   Location[3][2] = 0.577350269189626;   Weight[3] = 1.0;
    Location[4][0] = 0.577350269189626;   Location[4][1] = -0.577350269189626;  Location[4][2] = -0.577350269189626;  Weight[4] = 1.0;
    Location[5][0] = 0.577350269189626;   Location[5][1] = -0.577350269189626;  Location[5][2] = 0.577350269189626;   Weight[5] = 1.0;
    Location[6][0] = 0.577350269189626;   Location[6][1] = 0.577350269189626;   Location[6][2] = -0.577350269189626;  Weight[6] = 1.0;
    Location[7][0] = 0.577350269189626;   Location[7][1] = 0.577350269189626;   Location[7][2] = 0.577350269189626;   Weight[7] = 1.0;
  }
  
  for (iGauss = 0; iGauss < nGauss; iGauss++) {
    
    Xi = Location[iGauss][0]; Eta = Location[iGauss][1];  Mu = Location[iGauss][2];
    
    if (nNodes == 4) Det = ShapeFunc_Tetra(Xi, Eta, Mu, CoordCorners, DShapeFunction);
    if (nNodes == 5) Det = ShapeFunc_Pyram(Xi, Eta, Mu, CoordCorners, DShapeFunction);
    if (nNodes == 6) Det = ShapeFunc_Wedge(Xi, Eta, Mu, CoordCorners, DShapeFunction);
    if (nNodes == 8) Det = ShapeFunc_Hexa(Xi, Eta, Mu, CoordCorners, DShapeFunction);
    
    /*--- Compute the B Matrix ---*/
    
    for (iVar = 0; iVar < 6; iVar++)
      for (jVar = 0; jVar < nNodes*nVar; jVar++)
        B_Matrix[iVar][jVar] = 0.0;
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      B_Matrix[0][0+iNode*nVar] = DShapeFunction[iNode][0];
      B_Matrix[1][1+iNode*nVar] = DShapeFunction[iNode][1];
      B_Matrix[2][2+iNode*nVar] = DShapeFunction[iNode][2];
      
      B_Matrix[3][0+iNode*nVar] = DShapeFunction[iNode][1];
      B_Matrix[3][1+iNode*nVar] = DShapeFunction[iNode][0];
      
      B_Matrix[4][1+iNode*nVar] = DShapeFunction[iNode][2];
      B_Matrix[4][2+iNode*nVar] = DShapeFunction[iNode][1];
      
      B_Matrix[5][0+iNode*nVar] = DShapeFunction[iNode][2];
      B_Matrix[5][2+iNode*nVar] = DShapeFunction[iNode][0];
    }
    
    /*--- Compute the D Matrix (for plane strain and 3-D)---*/
    
    D_Matrix[0][0] = Lambda + 2.0*Mu;	D_Matrix[0][1] = Lambda;					D_Matrix[0][2] = Lambda;					D_Matrix[0][3] = 0.0;	D_Matrix[0][4] = 0.0;	D_Matrix[0][5] = 0.0;
    D_Matrix[1][0] = Lambda;					D_Matrix[1][1] = Lambda + 2.0*Mu;	D_Matrix[1][2] = Lambda;					D_Matrix[1][3] = 0.0;	D_Matrix[1][4] = 0.0;	D_Matrix[1][5] = 0.0;
    D_Matrix[2][0] = Lambda;					D_Matrix[2][1] = Lambda;					D_Matrix[2][2] = Lambda + 2.0*Mu;	D_Matrix[2][3] = 0.0;	D_Matrix[2][4] = 0.0;	D_Matrix[2][5] = 0.0;
    D_Matrix[3][0] = 0.0;							D_Matrix[3][1] = 0.0;							D_Matrix[3][2] = 0.0;							D_Matrix[3][3] = Mu;	D_Matrix[3][4] = 0.0;	D_Matrix[3][5] = 0.0;
    D_Matrix[4][0] = 0.0;							D_Matrix[4][1] = 0.0;							D_Matrix[4][2] = 0.0;							D_Matrix[4][3] = 0.0;	D_Matrix[4][4] = Mu;	D_Matrix[4][5] = 0.0;
    D_Matrix[5][0] = 0.0;							D_Matrix[5][1] = 0.0;							D_Matrix[5][2] = 0.0;							D_Matrix[5][3] = 0.0;	D_Matrix[5][4] = 0.0;	D_Matrix[5][5] = Mu;
    
    
    /*--- Compute the BT.D Matrix ---*/
    
    for (iVar = 0; iVar < nNodes*nVar; iVar++) {
      for (jVar = 0; jVar < 6; jVar++) {
        Aux_Matrix[iVar][jVar] = 0.0;
        for (kVar = 0; kVar < 6; kVar++)
          Aux_Matrix[iVar][jVar] += B_Matrix[kVar][iVar]*D_Matrix[kVar][jVar];
      }
    }
    
    /*--- Compute the BT.D.B Matrix (stiffness matrix), and add to the original
     matrix using Gauss integration ---*/
    
    for (iVar = 0; iVar < nNodes*nVar; iVar++) {
      for (jVar = 0; jVar < nNodes*nVar; jVar++) {
        for (kVar = 0; kVar < 6; kVar++) {
          StiffMatrix_Elem[iVar][jVar] += Weight[iGauss] * Aux_Matrix[iVar][kVar]*B_Matrix[kVar][jVar] * Det;
        }
      }
    }
    
  }
  
}
