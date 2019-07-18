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

  /*--- Allocate internal structures ---*/

  nNodes = 3;
  nGaussPoints = 1;
  AllocateStructures(config->GetDeadLoad());

  /*--- Gauss coordinates and weights ---*/

  GaussCoord[0][0] = 1.0/3.0;  GaussCoord[0][1] = 1.0/3.0;  GaussWeight[0] = 0.5;

  /*--- Store the values of the shape functions and their derivatives ---*/

  unsigned short iGauss;
  su2double Xi, Eta, val_Ni;

  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {

    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];

    val_Ni = Xi;        GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = Eta;       GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 1-Xi-Eta;  GaussPoint[iGauss]->SetNi(val_Ni,2);

    /*--- dN/d xi, dN/d eta ---*/

    dNiXj[iGauss][0][0] =  1.0;  dNiXj[iGauss][0][1] =  0.0;
    dNiXj[iGauss][1][0] =  0.0;  dNiXj[iGauss][1][1] =  1.0;
    dNiXj[iGauss][2][0] = -1.0;  dNiXj[iGauss][2][1] = -1.0;

  }

  /*--- Shape functions evaluated at the nodes for extrapolation of the stresses at the Gaussian Points ---*/
  /*--- The stress is constant over a TRIA1 element ---*/

  NodalExtrap[0][0] = 1.0;
  NodalExtrap[1][0] = 1.0;
  NodalExtrap[2][0] = 1.0;

}

su2double CTRIA1::ComputeArea(const FrameType mode){
  
  unsigned short iDim;
  su2double a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0};
  su2double Area = 0.0;
  
  /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
        for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed) ---*/
  su2double **Coord = (mode==REFERENCE) ? RefCoord : CurrentCoord;
  
  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = Coord[0][iDim]-Coord[2][iDim];
    b[iDim] = Coord[1][iDim]-Coord[2][iDim];
  }
  
  Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  
  return Area;
  
}

CTRIA1::~CTRIA1(void) {

}

CQUAD4::CQUAD4(void) : CElement() {

}

CQUAD4::CQUAD4(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {

  /*--- Allocate internal structures ---*/

  nNodes = 4;
  nGaussPoints = 4;
  AllocateStructures(config->GetDeadLoad());

  /*--- Gauss coordinates and weights ---*/

  su2double oneOnSqrt3 = 0.577350269189626;

  GaussCoord[0][0] = -oneOnSqrt3;  GaussCoord[0][1] = -oneOnSqrt3;  GaussWeight[0] = 1.0;
  GaussCoord[1][0] =  oneOnSqrt3;  GaussCoord[1][1] = -oneOnSqrt3;  GaussWeight[1] = 1.0;
  GaussCoord[2][0] =  oneOnSqrt3;  GaussCoord[2][1] =  oneOnSqrt3;  GaussWeight[2] = 1.0;
  GaussCoord[3][0] = -oneOnSqrt3;  GaussCoord[3][1] =  oneOnSqrt3;  GaussWeight[3] = 1.0;

  /*--- Store the values of the shape functions and their derivatives ---*/

  unsigned short iNode, iGauss;
  su2double Xi, Eta, val_Ni;

  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {

    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];

    val_Ni = 0.25*(1.0-Xi)*(1.0-Eta);		GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = 0.25*(1.0+Xi)*(1.0-Eta);		GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 0.25*(1.0+Xi)*(1.0+Eta);		GaussPoint[iGauss]->SetNi(val_Ni,2);
    val_Ni = 0.25*(1.0-Xi)*(1.0+Eta);		GaussPoint[iGauss]->SetNi(val_Ni,3);

    /*--- dN/d xi, dN/d eta ---*/

    dNiXj[iGauss][0][0] = -0.25*(1.0-Eta);  dNiXj[iGauss][0][1] = -0.25*(1.0-Xi);
    dNiXj[iGauss][1][0] =  0.25*(1.0-Eta);  dNiXj[iGauss][1][1] = -0.25*(1.0+Xi);
    dNiXj[iGauss][2][0] =  0.25*(1.0+Eta);  dNiXj[iGauss][2][1] =  0.25*(1.0+Xi);
    dNiXj[iGauss][3][0] = -0.25*(1.0+Eta);  dNiXj[iGauss][3][1] =  0.25*(1.0-Xi);

  }

  /*--- Store the extrapolation functions (used to compute nodal stresses) ---*/

  su2double ExtrapCoord[4][2], sqrt3 = 1.732050807568877;;

  ExtrapCoord[0][0] = -sqrt3;  ExtrapCoord[0][1] = -sqrt3;
  ExtrapCoord[1][0] =  sqrt3;  ExtrapCoord[1][1] = -sqrt3;
  ExtrapCoord[2][0] =  sqrt3;  ExtrapCoord[2][1] =  sqrt3;
  ExtrapCoord[3][0] = -sqrt3;  ExtrapCoord[3][1] =  sqrt3;

  for (iNode = 0; iNode < nNodes; iNode++) {

    Xi = ExtrapCoord[iNode][0];
    Eta = ExtrapCoord[iNode][1];

    NodalExtrap[iNode][0] = 0.25*(1.0-Xi)*(1.0-Eta);
    NodalExtrap[iNode][1] = 0.25*(1.0+Xi)*(1.0-Eta);
    NodalExtrap[iNode][2] = 0.25*(1.0+Xi)*(1.0+Eta);
    NodalExtrap[iNode][3] = 0.25*(1.0-Xi)*(1.0+Eta);

  }

}

su2double CQUAD4::ComputeArea(const FrameType mode){
  
  unsigned short iDim;
  su2double a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0};
  su2double Area = 0.0;
  
  /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
        for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed)---*/
  su2double **Coord = (mode==REFERENCE) ? RefCoord : CurrentCoord;
  
  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = Coord[0][iDim]-Coord[2][iDim];
    b[iDim] = Coord[1][iDim]-Coord[2][iDim];
  }
  
  Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  
  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = Coord[0][iDim]-Coord[3][iDim];
    b[iDim] = Coord[2][iDim]-Coord[3][iDim];
  }
  
  Area += 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  
  return Area;
  
}

CQUAD4::~CQUAD4(void) {

}

CTETRA1::CTETRA1(void) : CElement() {

}

CTETRA1::CTETRA1(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {

  /*--- Allocate internal structures ---*/

  nNodes = 4;
  nGaussPoints = 1;
  AllocateStructures(config->GetDeadLoad());

  /*--- Gauss coordinates and weights ---*/

  GaussCoord[0][0] = 0.25;  GaussCoord[0][1] = 0.25; GaussCoord[0][2] = 0.25;  GaussWeight[0] = 1.0/6.0;

  /*--- Store the values of the shape functions and their derivatives ---*/

  unsigned short iGauss;
  su2double Xi, Eta, Zeta, val_Ni;

  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {

    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];

    val_Ni = Xi;						  GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = Eta;						  GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 1.0-Xi-Eta-Zeta;	GaussPoint[iGauss]->SetNi(val_Ni,2);
    val_Ni = Zeta;					  GaussPoint[iGauss]->SetNi(val_Ni,3);

    /*--- dN/d xi, dN/d eta, dN/d zeta ---*/

    dNiXj[iGauss][0][0] =  1.0;  dNiXj[iGauss][0][1] =  0.0;  dNiXj[iGauss][0][2] =  0.0;
    dNiXj[iGauss][1][0] =  0.0;  dNiXj[iGauss][1][1] =  1.0;  dNiXj[iGauss][1][2] =  0.0;
    dNiXj[iGauss][2][0] = -1.0;  dNiXj[iGauss][2][1] = -1.0;  dNiXj[iGauss][2][2] = -1.0;
    dNiXj[iGauss][3][0] =  0.0;  dNiXj[iGauss][3][1] =  0.0;  dNiXj[iGauss][3][2] =  1.0;

  }

  /*--- Shape functions evaluated at the nodes for extrapolation of the stresses at the Gaussian Points ---*/
  /*--- The stress is constant at a TETRA1 element ---*/

  NodalExtrap[0][0] = 1.0;
  NodalExtrap[1][0] = 1.0;
  NodalExtrap[2][0] = 1.0;
  NodalExtrap[3][0] = 1.0;

}

su2double CTETRA1::ComputeVolume(const FrameType mode){

  unsigned short iDim;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0};
  su2double Volume = 0.0;
  
  /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
        for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed)---*/
  su2double **Coord = (mode==REFERENCE) ? RefCoord : CurrentCoord;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord[1][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[2][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[3][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  return Volume;

}

CTETRA1::~CTETRA1(void) {

}

CHEXA8::CHEXA8(void) : CElement() {

}

CHEXA8::CHEXA8(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {

  /*--- Allocate internal structures ---*/

  nNodes = 8;
  nGaussPoints = 8;
  AllocateStructures(config->GetDeadLoad());

  /*--- Gauss coordinates and weights ---*/

  su2double oneOnSqrt3 = 0.577350269189626;

  GaussCoord[0][0] = -oneOnSqrt3;  GaussCoord[0][1] = -oneOnSqrt3;  GaussCoord[0][2] = -oneOnSqrt3;	 GaussWeight[0] = 1.0;
  GaussCoord[1][0] =  oneOnSqrt3;  GaussCoord[1][1] = -oneOnSqrt3;  GaussCoord[1][2] = -oneOnSqrt3;  GaussWeight[1] = 1.0;
  GaussCoord[2][0] =  oneOnSqrt3;  GaussCoord[2][1] =  oneOnSqrt3;  GaussCoord[2][2] = -oneOnSqrt3;  GaussWeight[2] = 1.0;
  GaussCoord[3][0] = -oneOnSqrt3;  GaussCoord[3][1] =  oneOnSqrt3;  GaussCoord[3][2] = -oneOnSqrt3;  GaussWeight[3] = 1.0;
  GaussCoord[4][0] = -oneOnSqrt3;  GaussCoord[4][1] = -oneOnSqrt3;  GaussCoord[4][2] =  oneOnSqrt3;  GaussWeight[4] = 1.0;
  GaussCoord[5][0] =  oneOnSqrt3;  GaussCoord[5][1] = -oneOnSqrt3;  GaussCoord[5][2] =  oneOnSqrt3;  GaussWeight[5] = 1.0;
  GaussCoord[6][0] =  oneOnSqrt3;  GaussCoord[6][1] =  oneOnSqrt3;  GaussCoord[6][2] =  oneOnSqrt3;  GaussWeight[6] = 1.0;
  GaussCoord[7][0] = -oneOnSqrt3;  GaussCoord[7][1] =  oneOnSqrt3;  GaussCoord[7][2] =  oneOnSqrt3;  GaussWeight[7] = 1.0;

  /*--- Store the values of the shape functions and their derivatives ---*/

  unsigned short iNode, iGauss;
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

    /*--- dN/d xi ---*/

    dNiXj[iGauss][0][0] = -0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[iGauss][1][0] =  0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[iGauss][2][0] =  0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[iGauss][3][0] = -0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[iGauss][4][0] = -0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[iGauss][5][0] =  0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[iGauss][6][0] =  0.125*(1.0+Eta)*(1.0+Zeta);
    dNiXj[iGauss][7][0] = -0.125*(1.0+Eta)*(1.0+Zeta);

    /*--- dN/d eta ---*/

    dNiXj[iGauss][0][1] = -0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[iGauss][1][1] = -0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[iGauss][2][1] =  0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[iGauss][3][1] =  0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[iGauss][4][1] = -0.125*(1.0-Xi)*(1.0+Zeta);
    dNiXj[iGauss][5][1] = -0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[iGauss][6][1] =  0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[iGauss][7][1] =  0.125*(1.0-Xi)*(1.0+Zeta);

    /*--- dN/d zeta ---*/

    dNiXj[iGauss][0][2] = -0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[iGauss][1][2] = -0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[iGauss][2][2] = -0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[iGauss][3][2] = -0.125*(1.0-Xi)*(1.0+Eta);
    dNiXj[iGauss][4][2] =  0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[iGauss][5][2] =  0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[iGauss][6][2] =  0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[iGauss][7][2] =  0.125*(1.0-Xi)*(1.0+Eta);

  }

  /*--- Store the extrapolation functions ---*/

  su2double ExtrapCoord[8][3], sqrt3 = 1.732050807568877;

  ExtrapCoord[0][0] = -sqrt3;  ExtrapCoord[0][1] = -sqrt3;  ExtrapCoord[0][2] = -sqrt3;
  ExtrapCoord[1][0] =  sqrt3;  ExtrapCoord[1][1] = -sqrt3;  ExtrapCoord[1][2] = -sqrt3;
  ExtrapCoord[2][0] =  sqrt3;  ExtrapCoord[2][1] =  sqrt3;  ExtrapCoord[2][2] = -sqrt3;
  ExtrapCoord[3][0] = -sqrt3;  ExtrapCoord[3][1] =  sqrt3;  ExtrapCoord[3][2] = -sqrt3;
  ExtrapCoord[4][0] = -sqrt3;  ExtrapCoord[4][1] = -sqrt3;  ExtrapCoord[4][2] =  sqrt3;
  ExtrapCoord[5][0] =  sqrt3;  ExtrapCoord[5][1] = -sqrt3;  ExtrapCoord[5][2] =  sqrt3;
  ExtrapCoord[6][0] =  sqrt3;  ExtrapCoord[6][1] =  sqrt3;  ExtrapCoord[6][2] =  sqrt3;
  ExtrapCoord[7][0] = -sqrt3;  ExtrapCoord[7][1] =  sqrt3;  ExtrapCoord[7][2] =  sqrt3;

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

su2double CHEXA8::ComputeVolume(const FrameType mode){

  unsigned short iDim;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0};
  su2double Volume = 0.0;
  
  /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
        for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed)---*/
  su2double **Coord = (mode==REFERENCE) ? RefCoord : CurrentCoord;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord[1][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[2][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[5][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord[2][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[7][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[5][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord[2][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[3][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[7][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord[5][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[7][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[4][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord[7][iDim] - Coord[2][iDim];
    r2[iDim] = Coord[5][iDim] - Coord[2][iDim];
    r3[iDim] = Coord[6][iDim] - Coord[2][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  return Volume;

}

CHEXA8::~CHEXA8(void) {

}

CPYRAM5::CPYRAM5(void) : CElement() {

}

CPYRAM5::CPYRAM5(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {

  /*--- Allocate internal structures ---*/

  nNodes = 5;
  nGaussPoints = 5;
  AllocateStructures(config->GetDeadLoad());

  /*--- Gauss coordinates and weights ---*/

  GaussCoord[0][0] = 0.5;  GaussCoord[0][1] = 0.0;  GaussCoord[0][2] = 0.1531754163448146;  GaussWeight[0] = 2.0/15.0;
  GaussCoord[1][0] = 0.0;  GaussCoord[1][1] = 0.5;  GaussCoord[1][2] = 0.1531754163448146;  GaussWeight[1] = 2.0/15.0;
  GaussCoord[2][0] =-0.5;  GaussCoord[2][1] = 0.0;  GaussCoord[2][2] = 0.1531754163448146;  GaussWeight[2] = 2.0/15.0;
  GaussCoord[3][0] = 0.0;  GaussCoord[3][1] =-0.5;  GaussCoord[3][2] = 0.1531754163448146;  GaussWeight[3] = 2.0/15.0;
  GaussCoord[4][0] = 0.0;  GaussCoord[4][1] = 0.0;  GaussCoord[4][2] = 0.6372983346207416;  GaussWeight[4] = 2.0/15.0;

  /*--- Store the values of the shape functions and their derivatives ---*/

  unsigned short iNode, iGauss;
  su2double Xi, Eta, Zeta, val_Ni;

  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {

    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];

    val_Ni = 0.25*(-Xi+Eta+Zeta-1.0)*(-Xi-Eta+Zeta-1.0)/(1.0-Zeta);  GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = 0.25*(-Xi-Eta+Zeta-1.0)*( Xi-Eta+Zeta-1.0)/(1.0-Zeta);  GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 0.25*( Xi+Eta+Zeta-1.0)*( Xi-Eta+Zeta-1.0)/(1.0-Zeta);  GaussPoint[iGauss]->SetNi(val_Ni,2);
    val_Ni = 0.25*( Xi+Eta+Zeta-1.0)*(-Xi+Eta+Zeta-1.0)/(1.0-Zeta);  GaussPoint[iGauss]->SetNi(val_Ni,3);
    val_Ni = Zeta;                                                   GaussPoint[iGauss]->SetNi(val_Ni,4);

    /*--- dN/d xi ---*/

    dNiXj[iGauss][0][0] = 0.5*(Zeta-Xi-1.0)/(Zeta-1.0);
    dNiXj[iGauss][1][0] = 0.5*Xi/(Zeta-1.0);
    dNiXj[iGauss][2][0] = 0.5*(1.0-Zeta-Xi)/(Zeta-1.0);
    dNiXj[iGauss][3][0] = dNiXj[iGauss][1][0];
    dNiXj[iGauss][4][0] = 0.0;

    /*--- dN/d eta ---*/

    dNiXj[iGauss][0][1] = 0.5*Eta/(Zeta-1.0);
    dNiXj[iGauss][1][1] = 0.5*(Zeta-Eta-1.0)/(Zeta-1.0);
    dNiXj[iGauss][2][1] = dNiXj[iGauss][0][1];
    dNiXj[iGauss][3][1] = 0.5*(1.0-Zeta-Eta)/(Zeta-1.0);
    dNiXj[iGauss][4][1] = 0.0;

    /*--- dN/d zeta ---*/

    dNiXj[iGauss][0][2] = 0.25*(-1.0 + 2.0*Zeta - Zeta*Zeta - Eta*Eta + Xi*Xi)/((1.0-Zeta)*(1.0-Zeta));
    dNiXj[iGauss][1][2] = 0.25*(-1.0 + 2.0*Zeta - Zeta*Zeta + Eta*Eta - Xi*Xi)/((1.0-Zeta)*(1.0-Zeta));
    dNiXj[iGauss][2][2] = dNiXj[iGauss][0][2];
    dNiXj[iGauss][3][2] = dNiXj[iGauss][1][2];
    dNiXj[iGauss][4][2] = 1.0;

  }

  /*--- Store the extrapolation functions ---*/

  su2double ExtrapCoord[5][3];

  ExtrapCoord[0][0] =  2.0;  ExtrapCoord[0][1] =  0.0;  ExtrapCoord[0][2] = -0.316397779494322;
  ExtrapCoord[1][0] =  0.0;  ExtrapCoord[1][1] =  2.0;  ExtrapCoord[1][2] = -0.316397779494322;
  ExtrapCoord[2][0] = -2.0;  ExtrapCoord[2][1] =  0.0;  ExtrapCoord[2][2] = -0.316397779494322;
  ExtrapCoord[3][0] =  0.0;  ExtrapCoord[3][1] = -2.0;  ExtrapCoord[3][2] = -0.316397779494322;
  ExtrapCoord[4][0] =  0.0;  ExtrapCoord[4][1] =  0.0;  ExtrapCoord[4][2] =  1.749193338482970;

  for (iNode = 0; iNode < nNodes; iNode++) {

    Xi = ExtrapCoord[iNode][0];
    Eta = ExtrapCoord[iNode][1];
    Zeta = ExtrapCoord[iNode][2];

    NodalExtrap[iNode][0] = 0.25*(-Xi+Eta+Zeta-1.0)*(-Xi-Eta+Zeta-1.0)/(1.0-Zeta);
    NodalExtrap[iNode][1] = 0.25*(-Xi-Eta+Zeta-1.0)*( Xi-Eta+Zeta-1.0)/(1.0-Zeta);
    NodalExtrap[iNode][2] = 0.25*( Xi+Eta+Zeta-1.0)*( Xi-Eta+Zeta-1.0)/(1.0-Zeta);
    NodalExtrap[iNode][3] = 0.25*( Xi+Eta+Zeta-1.0)*(-Xi+Eta+Zeta-1.0)/(1.0-Zeta);
    NodalExtrap[iNode][4] = Zeta;

  }

}

su2double CPYRAM5::ComputeVolume(const FrameType mode){

  unsigned short iDim;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0};
  su2double Volume = 0.0;
  
  /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
        for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed)---*/
  su2double **Coord = (mode==REFERENCE) ? RefCoord : CurrentCoord;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord[1][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[2][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[4][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord[2][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[3][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[4][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  return Volume;

}

CPYRAM5::~CPYRAM5(void) {

}

CPRISM6::CPRISM6(void) : CElement() {

}

CPRISM6::CPRISM6(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {

  /*--- Allocate internal structures ---*/

  nNodes = 6;
  nGaussPoints = 6;
  AllocateStructures(config->GetDeadLoad());

  /*--- Gauss coordinates and weights ---*/

  /*--- There is some inconsistency between the shape functions and the order of the nodes
        that causes "negative" stiffness, the remedy is to use negative weights. ---*/

  su2double oneOnSqrt3 = 0.577350269189626;
  GaussCoord[0][0] = -oneOnSqrt3;  GaussCoord[0][1] = 1.0/6.0;  GaussCoord[0][2] = 1.0/6.0;  GaussWeight[0] = -1.0/6.0;
  GaussCoord[1][0] = -oneOnSqrt3;  GaussCoord[1][1] = 2.0/3.0;  GaussCoord[1][2] = 1.0/6.0;  GaussWeight[1] = -1.0/6.0;
  GaussCoord[2][0] = -oneOnSqrt3;  GaussCoord[2][1] = 1.0/6.0;  GaussCoord[2][2] = 2.0/3.0;  GaussWeight[2] = -1.0/6.0;
  GaussCoord[3][0] =  oneOnSqrt3;  GaussCoord[3][1] = 1.0/6.0;  GaussCoord[3][2] = 1.0/6.0;  GaussWeight[3] = -1.0/6.0;
  GaussCoord[4][0] =  oneOnSqrt3;  GaussCoord[4][1] = 2.0/3.0;  GaussCoord[4][2] = 1.0/6.0;  GaussWeight[4] = -1.0/6.0;
  GaussCoord[5][0] =  oneOnSqrt3;  GaussCoord[5][1] = 1.0/6.0;  GaussCoord[5][2] = 2.0/3.0;  GaussWeight[5] = -1.0/6.0;

  /*--- Store the values of the shape functions and their derivatives ---*/

  unsigned short iNode, iGauss;
  su2double Xi, Eta, Zeta, val_Ni;

  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {

      Xi = GaussCoord[iGauss][0];
      Eta = GaussCoord[iGauss][1];
      Zeta = GaussCoord[iGauss][2];

      val_Ni = 0.5*Eta*(1.0-Xi);              GaussPoint[iGauss]->SetNi(val_Ni,0);
      val_Ni = 0.5*Zeta*(1.0-Xi);             GaussPoint[iGauss]->SetNi(val_Ni,1);
      val_Ni = 0.5*(1.0-Eta-Zeta)*(1.0-Xi);   GaussPoint[iGauss]->SetNi(val_Ni,2);
      val_Ni = 0.5*Eta*(Xi+1.0);              GaussPoint[iGauss]->SetNi(val_Ni,3);
      val_Ni = 0.5*Zeta*(Xi+1.0);             GaussPoint[iGauss]->SetNi(val_Ni,4);
      val_Ni = 0.5*(1.0-Eta-Zeta)*(Xi+1.0);   GaussPoint[iGauss]->SetNi(val_Ni,5);

      /*--- dN/d xi ---*/

      dNiXj[iGauss][0][0] = -0.5*Eta;
      dNiXj[iGauss][1][0] = -0.5*Zeta;
      dNiXj[iGauss][2][0] = -0.5*(1.0-Eta-Zeta);
      dNiXj[iGauss][3][0] =  0.5*Eta;
      dNiXj[iGauss][4][0] =  0.5*Zeta;
      dNiXj[iGauss][5][0] =  0.5*(1.0-Eta-Zeta);

      /*--- dN/d eta ---*/

      dNiXj[iGauss][0][1] =  0.5*(1.0-Xi);
      dNiXj[iGauss][1][1] =  0.0;
      dNiXj[iGauss][2][1] = -0.5*(1.0-Xi);
      dNiXj[iGauss][3][1] =  0.5*(Xi+1.0);
      dNiXj[iGauss][4][1] =  0.0;
      dNiXj[iGauss][5][1] = -0.5*(Xi+1.0);

      /*--- dN/d mu ---*/

      dNiXj[iGauss][0][2] =  0.0;
      dNiXj[iGauss][1][2] =  0.5*(1.0-Xi);
      dNiXj[iGauss][2][2] = -0.5*(1.0-Xi);
      dNiXj[iGauss][3][2] =  0.0;
      dNiXj[iGauss][4][2] =  0.5*(Xi+1.0);
      dNiXj[iGauss][5][2] = -0.5*(Xi+1.0);

  }

  /*--- Store the extrapolation functions ---*/

  su2double ExtrapCoord[6][3], sqrt3 = 1.732050807568877;

  ExtrapCoord[0][0] = -sqrt3;  ExtrapCoord[0][1] = -1.0/3.0;  ExtrapCoord[0][2] = -1.0/3.0;
  ExtrapCoord[1][0] = -sqrt3;  ExtrapCoord[1][1] =  5.0/3.0;  ExtrapCoord[1][2] = -1.0/3.0;
  ExtrapCoord[2][0] = -sqrt3;  ExtrapCoord[2][1] = -1.0/3.0;  ExtrapCoord[2][2] =  5.0/3.0;
  ExtrapCoord[3][0] =  sqrt3;  ExtrapCoord[3][1] = -1.0/3.0;  ExtrapCoord[3][2] = -1.0/3.0;
  ExtrapCoord[4][0] =  sqrt3;  ExtrapCoord[4][1] =  5.0/3.0;  ExtrapCoord[4][2] = -1.0/3.0;
  ExtrapCoord[5][0] =  sqrt3;  ExtrapCoord[5][1] = -1.0/3.0;  ExtrapCoord[5][2] =  5.0/3.0;

  for (iNode = 0; iNode < nNodes; iNode++) {

    Xi = ExtrapCoord[iNode][0];
    Eta = ExtrapCoord[iNode][1];
    Zeta = ExtrapCoord[iNode][2];

    NodalExtrap[iNode][0] = 0.5*Eta*(1.0-Xi);
    NodalExtrap[iNode][1] = 0.5*Zeta*(1.0-Xi);
    NodalExtrap[iNode][2] = 0.5*(1.0-Eta-Zeta)*(1.0-Xi);
    NodalExtrap[iNode][3] = 0.5*Eta*(Xi+1.0);
    NodalExtrap[iNode][4] = 0.5*Zeta*(Xi+1.0);
    NodalExtrap[iNode][5] = 0.5*(1.0-Eta-Zeta)*(Xi+1.0);

  }

}

su2double CPRISM6::ComputeVolume(const FrameType mode){

  unsigned short iDim;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0};
  su2double Volume = 0.0;

  /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
        for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed)---*/
  su2double **Coord = (mode==REFERENCE) ? RefCoord : CurrentCoord;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord[2][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[1][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[5][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord[5][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[1][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[4][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord[5][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[4][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[3][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

  return Volume;

}

CPRISM6::~CPRISM6(void) {

}
