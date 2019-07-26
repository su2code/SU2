/*!
 * \file numerics_direct_elasticity.cpp
 * \brief This file contains the routines for setting the tangent matrix and residual of a sobolev gradient problem
 * \author T.Dick
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

CGradSmoothing::CGradSmoothing(void) : CFEAElasticity () {

  unsigned short iDim;
  Aux_Mat = NULL;
  if (nDim == 2) {
    Aux_Mat = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Aux_Mat[iDim]      = new su2double[3];
    }
  }
  else if (nDim == 3) {
    Aux_Mat = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Aux_Mat[iDim]      = new su2double[6];
    }
  }

}

CGradSmoothing::CGradSmoothing(unsigned short val_nDim, CConfig *config)
  : CFEAElasticity(val_nDim, val_nDim, config) {

  unsigned short iDim;
  Aux_Mat = NULL;
  if (nDim == 2) {
    Aux_Mat = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Aux_Mat[iDim]      = new su2double[3];
    }
  }
  else if (nDim == 3) {
    Aux_Mat = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Aux_Mat[iDim]      = new su2double[6];
    }
  }

}

CGradSmoothing::~CGradSmoothing(void) {

  unsigned short iDim;
  for (iDim = 0; iDim < nDim; iDim++) {
    delete [] Aux_Mat[iDim];
  }
  delete [] Aux_Mat;

}

void CGradSmoothing::Compute_Tangent_Matrix(CElement *element, CConfig *config) {

  unsigned short iDim, jDim, kDim;
  unsigned short iGauss, nGauss;
  unsigned short iNode, jNode, nNode;
  unsigned short bDim;

  su2double Weight, Jac_X, GradNiXGradNj;

  su2double val_HiHj;
  su2double *val_DHiHj;
  val_DHiHj = new su2double[nDim];

  /*--- Set element properties  if they need to be set ---*/
  SetElement_Properties(element, config);

  Compute_Constitutive_Matrix(element, config);

  element->clearElement(true);       /*--- Restarts the element: avoids adding over previous results in other elements --*/
  element->ComputeGrad_Linear();
  nNode = element->GetnNodes();
  nGauss = element->GetnGaussPoints();

  /*--- contribution from the gradients of the shape functions ---*/

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    Weight = element->GetWeight(iGauss);
    Jac_X = element->GetJ_X(iGauss);

    for (iNode = 0; iNode < nNode; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
        // std::cout << "GradNi: " << iNode << ", " << iDim << ", " << GradNi_Ref_Mat[iNode][iDim] << ", " << element->GetRef_Coord(iNode,iDim) << std::endl;
      }
    }

    for (iNode = 0; iNode < nNode; iNode++) {

      /*--- Assumming symmetry ---*/
      for (jNode = iNode; jNode < nNode; jNode++) {

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            GradNiXGradNj += GradNi_Ref_Mat[iNode][iDim]* GradNi_Ref_Mat[jNode][jDim];
          }
        }

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            if (iDim == jDim) {
              KAux_ab[iDim][jDim] = Weight * Jac_X * GradNiXGradNj;
              // std::cout << "Grad Term: " << iGauss << ", " << iNode << ", " << jNode << ", " << Weight << ", " << Jac_X << ", " << GradNiXGradNj << std::endl;
            } else {
              KAux_ab[iDim][jDim] = 0;
            }
          }
        }

        GradNiXGradNj=0;

        element->Add_DHiDHj(KAux_ab,iNode, jNode);
        /*--- Symmetric terms --*/
        if (iNode != jNode) {
          element->Add_DHiDHj_T(KAux_ab, jNode, iNode);
        }

      }

    }

  }

  /*--- contribution from the shape functions themselfes --*/

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    Weight = element->GetWeight(iGauss);
    // do we need the determinant of the Jacobian in this case?
    // Jac_X = element->GetJ_X(iGauss);

    for (iNode = 0; iNode < nNode; iNode++) {
      Ni_Vec[iNode] = element->GetNi(iNode,iGauss);
    }

    for (iNode = 0; iNode < nNode; iNode++) {
      for (jNode = 0; jNode < nNode; jNode++) {
        val_HiHj = Weight * Ni_Vec[iNode] * Ni_Vec[jNode];
        // std::cout << "Function Term: " << iGauss << ", " << iNode << ", " << jNode << ", " << Weight << ", " << Jac_X << ", " << val_HiHj << std::endl;
        element->Add_HiHj(val_HiHj, iNode, jNode);
      }
    }

  }

  /*--- contribution from Neumann boundary terms ---*/
  /*
  for (iGauss = 0; iGauss < nGauss; iGauss++) {


    Weight = element->GetWeight(iGauss);

    for (iNode = 0; iNode < nNode; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
      }
      Ni_Vec[iNode] = element->GetNi(iNode,iGauss);
    }

    for (iNode = 0; iNode < nNode; iNode++) {
      for (jNode = 0; jNode < nNode; jNode++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          val_DHiHj[iDim] =  Weight * GradNi_Ref_Mat[iNode][iDim] * Ni_Vec[jNode];
        }
        element->Add_DHiHj(valDHiHj, iNode, jNode);
        if (iNode != jNode) {
          element->Add_DHiHj(val_DHiHj, jNode, iNode);
        }
      }

    }

  }
  */

  delete [] val_DHiHj;

}

void CGradSmoothing::SetElement_Properties(CElement *element, CConfig *config) {

 // no special element properties needed so far ???

}

void CGradSmoothing::Compute_Constitutive_Matrix(CElement *element_container, CConfig *config) {

  /*--- Compute the D Matrix, so far identity for test purposes ---*/

  if (nDim == 2) {

    D_Mat[0][0] = 1.0;  D_Mat[0][1] = 0.0;  D_Mat[0][2] = 0.0;
    D_Mat[1][0] = 0.0;  D_Mat[1][1] = 1.0;  D_Mat[1][2] = 0.0;
    D_Mat[2][0] = 0.0;  D_Mat[2][1] = 0.0;  D_Mat[2][2] = 1.0;
  }
  else if (nDim == 3) {

    D_Mat[0][0] = 1.0;  D_Mat[0][1] = 0.0;  D_Mat[0][2] = 0.0;  D_Mat[0][3] = 0.0;  D_Mat[0][4] = 0.0;  D_Mat[0][5] = 0.0;
    D_Mat[1][0] = 0.0;  D_Mat[1][1] = 1.0;  D_Mat[1][2] = 0.0;  D_Mat[1][3] = 0.0;  D_Mat[1][4] = 0.0;  D_Mat[1][5] = 0.0;
    D_Mat[2][0] = 0.0;  D_Mat[2][1] = 0.0;  D_Mat[2][2] = 1.0;  D_Mat[2][3] = 0.0;  D_Mat[2][4] = 0.0;  D_Mat[2][5] = 0.0;
    D_Mat[3][0] = 0.0;  D_Mat[3][1] = 0.0;  D_Mat[3][2] = 0.0;  D_Mat[3][3] = 1.0;  D_Mat[3][4] = 0.0;  D_Mat[3][5] = 0.0;
    D_Mat[4][0] = 0.0;  D_Mat[4][1] = 0.0;  D_Mat[4][2] = 0.0;  D_Mat[4][3] = 0.0;  D_Mat[4][4] = 1.0;  D_Mat[4][5] = 0.0;
    D_Mat[5][0] = 0.0;  D_Mat[5][1] = 0.0;  D_Mat[5][2] = 0.0;  D_Mat[5][3] = 0.0;  D_Mat[5][4] = 0.0;  D_Mat[5][5] = 1.0;

  }

}
