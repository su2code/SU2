/*!
 * \file numerics_direct_elasticity_nonlinear.cpp
 * \brief This file contains the routines for setting the tangent matrix and residual of a FEM nonlinear elastic structural problem.
 * \author R. Sanchez
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

CFEM_NonlinearElasticity::CFEM_NonlinearElasticity(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_Elasticity(val_nDim, val_nVar, config) {

  unsigned short iVar;

  F_Mat = new su2double *[3];
  b_Mat = new su2double *[3];
  Stress_Tensor = new su2double *[3];
  for (iVar = 0; iVar < 3; iVar++) {
    F_Mat[iVar] = new su2double [3];
    b_Mat[iVar] = new su2double [3];
    Stress_Tensor[iVar] = new su2double [3];
  }

  KAux_t_a = new su2double [nDim];

  KAux_P_ab = new su2double* [nDim];
  for (iVar = 0; iVar < nDim; iVar++) {
    KAux_P_ab[iVar] = new su2double[nDim];
  }

  if (nDim == 2) {
    currentCoord = new su2double* [4];  /*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
    for (iVar = 0; iVar < 4; iVar++) currentCoord[iVar] = new su2double[nDim];
  }
  else if (nDim == 3) {
    currentCoord = new su2double* [8];  /*--- As of now, 8 is the maximum number of nodes for 3D problems ---*/
    for (iVar = 0; iVar < 8; iVar++) currentCoord[iVar] = new su2double[nDim];
  }

  J_F = 0.0;
  f33 = 1.0;


}

CFEM_NonlinearElasticity::~CFEM_NonlinearElasticity(void) {

  unsigned short iVar;

  for (iVar = 0; iVar < 3; iVar++) {
    delete [] F_Mat[iVar];
    delete [] b_Mat[iVar];
    delete [] Stress_Tensor[iVar];
  }

  for (iVar = 0; iVar < nDim; iVar++) {
    delete [] KAux_P_ab[iVar];
  }

  if (nDim == 2) {
    for (iVar = 0; iVar < 4; iVar++) {
      delete [] currentCoord[iVar];
    }
  }
  else if (nDim == 3) {
    for (iVar = 0; iVar < 8; iVar++) {
      delete [] currentCoord[iVar];
    }
  }

  delete [] F_Mat;
  delete [] b_Mat;
  delete [] Stress_Tensor;
  delete [] KAux_t_a;
  delete [] KAux_P_ab;
  delete [] currentCoord;

}


void CFEM_NonlinearElasticity::Compute_Tangent_Matrix(CElement *element, CConfig *config) {

  unsigned short iVar, jVar, kVar;
  unsigned short iGauss, nGauss;
  unsigned short iNode, jNode, nNode;
  unsigned short iDim, bDim;

  su2double Ks_Aux_ab;

  su2double Weight, Jac_x;

  su2double AuxMatrixKc[3][6];
  su2double AuxMatrixKs[3];

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
      AuxMatrixKc[iVar][jVar] = 0.0;
    }
  }

  for (iVar = 0; iVar < 3; iVar++) {
    AuxMatrixKs[iVar] = 0.0;
  }

  element->clearElement();       /*--- Restarts the element: avoids adding over previous results in other elements --*/
  element->ComputeGrad_NonLinear();

  nNode = element->GetnNodes();
  nGauss = element->GetnGaussPoints();

  /*--- Full integration of the constitutive and stress term ---*/

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    Weight = element->GetWeight(iGauss);
    Jac_x = element->GetJ_x(iGauss);

    /*--- Initialize the deformation gradient for each Gauss Point ---*/

    for (iVar = 0; iVar < 3; iVar++) {
      for (jVar = 0; jVar < 3; jVar++) {
        F_Mat[iVar][jVar] = 0.0;
        b_Mat[iVar][jVar] = 0.0;
      }
    }

    /*--- Retrieve the values of the gradients of the shape functions for each node ---*/
    /*--- This avoids repeated operations ---*/

    for (iNode = 0; iNode < nNode; iNode++) {

      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
        GradNi_Curr_Mat[iNode][iDim] = element->GetGradNi_x(iNode,iGauss,iDim);
        currentCoord[iNode][iDim] = element->GetCurr_Coord(iNode, iDim);
      }

      /*--- Compute the deformation gradient ---*/

      for (iVar = 0; iVar < nDim; iVar++) {
        for (jVar = 0; jVar < nDim; jVar++) {
          F_Mat[iVar][jVar] += currentCoord[iNode][iVar]*GradNi_Ref_Mat[iNode][jVar];
        }
      }

      /*--- This implies plane strain --> Consider the possible implementation for plane stress --*/
      if (nDim == 2) {
        F_Mat[2][2] = 1.0;
      }

    }

    if (nDim == 2) {
      if (plane_stress) {
        // Compute the value of the term 33 for the deformation gradient
        Compute_Plane_Stress_Term(element, config);
        F_Mat[2][2] = f33;
      }
      else {
        F_Mat[2][2] = 1.0;
      }
    }

    /*--- Determinant of F --> Jacobian of the transformation ---*/

    J_F =   F_Mat[0][0]*F_Mat[1][1]*F_Mat[2][2]+
        F_Mat[0][1]*F_Mat[1][2]*F_Mat[2][0]+
        F_Mat[0][2]*F_Mat[1][0]*F_Mat[2][1]-
        F_Mat[0][2]*F_Mat[1][1]*F_Mat[2][0]-
        F_Mat[1][2]*F_Mat[2][1]*F_Mat[0][0]-
        F_Mat[2][2]*F_Mat[0][1]*F_Mat[1][0];

    /*--- Compute the left Cauchy deformation tensor ---*/

    for (iVar = 0; iVar < 3; iVar++) {
      for (jVar = 0; jVar < 3; jVar++) {
        for (kVar = 0; kVar < 3; kVar++) {
          b_Mat[iVar][jVar] += F_Mat[iVar][kVar]*F_Mat[jVar][kVar];
        }
      }
    }

    /*--- Compute the constitutive matrix ---*/

    Compute_Constitutive_Matrix(element, config);
    Compute_Stress_Tensor(element, config);


    for (iNode = 0; iNode < nNode; iNode++) {

      /*--------------------------------------------------------------------------------*/
      /*---------------------------- NODAL STRESS TERM ---------------------------------*/
      /*--------------------------------------------------------------------------------*/
        /*--- Compute the nodal stress term for each gaussian point and for each node, ---*/
        /*--- and add it to the element structure to be retrieved from the solver      ---*/

      for (iVar = 0; iVar < nDim; iVar++) {
        KAux_t_a[iVar] = 0.0;
        for (jVar = 0; jVar < nDim; jVar++) {
          KAux_t_a[iVar] += Weight * Stress_Tensor[iVar][jVar] * GradNi_Curr_Mat[iNode][jVar] * Jac_x;
        }
      }

      element->Add_Kt_a(KAux_t_a, iNode);

      /*--------------------------------------------------------------------------------*/
      /*----------------------- CONSTITUTIVE AND STRESS TERM ---------------------------*/
      /*--------------------------------------------------------------------------------*/

      if (nDim == 2) {
        Ba_Mat[0][0] = GradNi_Curr_Mat[iNode][0];
        Ba_Mat[1][1] = GradNi_Curr_Mat[iNode][1];
        Ba_Mat[2][0] = GradNi_Curr_Mat[iNode][1];
        Ba_Mat[2][1] = GradNi_Curr_Mat[iNode][0];
      }
      else if (nDim ==3) {
        Ba_Mat[0][0] = GradNi_Curr_Mat[iNode][0];
        Ba_Mat[1][1] = GradNi_Curr_Mat[iNode][1];
        Ba_Mat[2][2] = GradNi_Curr_Mat[iNode][2];
        Ba_Mat[3][0] = GradNi_Curr_Mat[iNode][1];
        Ba_Mat[3][1] = GradNi_Curr_Mat[iNode][0];
        Ba_Mat[4][0] = GradNi_Curr_Mat[iNode][2];
        Ba_Mat[4][2] = GradNi_Curr_Mat[iNode][0];
        Ba_Mat[5][1] = GradNi_Curr_Mat[iNode][2];
        Ba_Mat[5][2] = GradNi_Curr_Mat[iNode][1];
      }

        /*--- Compute the BT.D Matrix ---*/

      for (iVar = 0; iVar < nDim; iVar++) {
        for (jVar = 0; jVar < bDim; jVar++) {
          AuxMatrixKc[iVar][jVar] = 0.0;
          for (kVar = 0; kVar < bDim; kVar++) {
            AuxMatrixKc[iVar][jVar] += Ba_Mat[kVar][iVar]*D_Mat[kVar][jVar];
          }
        }
      }

        /*--- Compute the BT.D Matrix ---*/

      for (iVar = 0; iVar < nDim; iVar++) {
        AuxMatrixKs[iVar] = 0.0;
        for (jVar = 0; jVar < nDim; jVar++) {
          AuxMatrixKs[iVar] += GradNi_Curr_Mat[iNode][jVar]*Stress_Tensor[jVar][iVar]; // DOUBLE CHECK
        }
      }

      /*--- Assumming symmetry ---*/
      for (jNode = iNode; jNode < nNode; jNode++) {
        if (nDim == 2) {
          Bb_Mat[0][0] = GradNi_Curr_Mat[jNode][0];
          Bb_Mat[1][1] = GradNi_Curr_Mat[jNode][1];
          Bb_Mat[2][0] = GradNi_Curr_Mat[jNode][1];
          Bb_Mat[2][1] = GradNi_Curr_Mat[jNode][0];
        }
        else if (nDim ==3) {
          Bb_Mat[0][0] = GradNi_Curr_Mat[jNode][0];
          Bb_Mat[1][1] = GradNi_Curr_Mat[jNode][1];
          Bb_Mat[2][2] = GradNi_Curr_Mat[jNode][2];
          Bb_Mat[3][0] = GradNi_Curr_Mat[jNode][1];
          Bb_Mat[3][1] = GradNi_Curr_Mat[jNode][0];
          Bb_Mat[4][0] = GradNi_Curr_Mat[jNode][2];
          Bb_Mat[4][2] = GradNi_Curr_Mat[jNode][0];
          Bb_Mat[5][1] = GradNi_Curr_Mat[jNode][2];
          Bb_Mat[5][2] = GradNi_Curr_Mat[jNode][1];
        }

        /*--- KAux_ab is the term for the constitutive part of the tangent matrix ---*/
        for (iVar = 0; iVar < nDim; iVar++) {
          for (jVar = 0; jVar < nDim; jVar++) {
            KAux_ab[iVar][jVar] = 0.0;
            for (kVar = 0; kVar < bDim; kVar++) {
              KAux_ab[iVar][jVar] += Weight * AuxMatrixKc[iVar][kVar] * Bb_Mat[kVar][jVar] * Jac_x;
            }
          }
        }

        /*--- Ks_Aux_ab is the term for the constitutive part of the tangent matrix ---*/
        Ks_Aux_ab = 0.0;
        for (iVar = 0; iVar < nDim; iVar++) {
          Ks_Aux_ab += Weight * AuxMatrixKs[iVar] * GradNi_Curr_Mat[jNode][iVar] * Jac_x;
        }

        element->Add_Kab(KAux_ab,iNode, jNode);
        element->Add_Ks_ab(Ks_Aux_ab,iNode, jNode);
        /*--- Symmetric terms --*/
        if (iNode != jNode) {
          element->Add_Kab_T(KAux_ab, jNode, iNode);
          element->Add_Ks_ab(Ks_Aux_ab,jNode, iNode);
        }

      }

    }

  }

}

void CFEM_NonlinearElasticity::Compute_MeanDilatation_Term(CElement *element, CConfig *config) {

  unsigned short iVar, jVar;
  unsigned short iGauss, nGauss;
  unsigned short iNode, jNode, nNode;
  su2double Weight, Jac_X, Jac_x;
  unsigned short iDim ;

  su2double GradNi_Mat_Term;
  su2double Vol_current, Vol_reference;
  su2double Avg_kappa;
  su2double el_Pressure;


  /*--- Under integration of the pressure term, if the calculations assume incompressibility or near incompressibility ---*/

  element->ComputeGrad_Pressure(); // Check if we can take this out!

  /*--- nGauss is here the number of Gaussian Points for the pressure term ---*/
  nGauss = element->GetnGaussPointsP();
  nNode = element->GetnNodes();

  /*--- Initialize the Gradient auxiliary Matrix ---*/
  for (iNode = 0; iNode < nNode; iNode++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      GradNi_Curr_Mat[iNode][iDim] = 0.0;
    }
  }

  Vol_current = 0.0;
  Vol_reference = 0.0;

  /*--------------------------------------------------------------------------------*/
  /*-------------------------- INCOMPRESSIBLE TERM ---------------------------------*/
  /*--------------------------------------------------------------------------------*/

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    Weight = element->GetWeight_P(iGauss);
    Jac_X = element->GetJ_X_P(iGauss);
    Jac_x = element->GetJ_x_P(iGauss);

    /*--- Retrieve the values of the gradients of the shape functions for each node ---*/
    /*--- This avoids repeated operations ---*/

    /*--- We compute the average gradient ---*/
    for (iNode = 0; iNode < nNode; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Mat_Term = element->GetGradNi_x_P(iNode,iGauss,iDim);
        GradNi_Curr_Mat[iNode][iDim] += Weight * GradNi_Mat_Term * Jac_x;
      }
    }

    Vol_reference += Weight * Jac_X;
    Vol_current += Weight * Jac_x;

  }

  if ((Vol_current > 0.0) && (Vol_reference > 0.0)) {

    /*--- It is necessary to divide over the current volume to obtain the averaged gradients ---*/
    for (iNode = 0; iNode < nNode; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Curr_Mat[iNode][iDim] = GradNi_Curr_Mat[iNode][iDim] / Vol_current;
      }
    }

    Avg_kappa = Kappa * Vol_current / Vol_reference;

    el_Pressure = Kappa * ((Vol_current / Vol_reference) - 1);

    element->SetElement_Pressure(el_Pressure);

  }
  else {
    cout << "Warning: Negative volume computed during FE structural analysis. Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  for (iNode = 0; iNode < nNode; iNode++) {

    for (jNode = 0; jNode < nNode; jNode++) {

      /*--- KAux_P_ab is the term for the incompressibility part of the tangent matrix ---*/
      for (iVar = 0; iVar < nDim; iVar++) {
        for (jVar = 0; jVar < nDim; jVar++) {
          KAux_P_ab[iVar][jVar] = Avg_kappa * Vol_current * GradNi_Curr_Mat[iNode][iVar] * GradNi_Curr_Mat[jNode][jVar];
        }
      }

      element->Set_Kk_ab(KAux_P_ab,iNode, jNode);

    }

  }

}


void CFEM_NonlinearElasticity::Compute_NodalStress_Term(CElement *element, CConfig *config) {

  unsigned short iVar, jVar, kVar;
  unsigned short iGauss, nGauss;
  unsigned short iNode, nNode;
  unsigned short iDim;

  su2double Weight, Jac_x;

  element->clearElement();       /*--- Restarts the element: avoids adding over previous results in other elements --*/
  element->ComputeGrad_NonLinear();  /*--- Check if we can take this out... so we don't have to do it twice ---*/

  nNode = element->GetnNodes();
  nGauss = element->GetnGaussPoints();

  /*--- Full integration of the nodal stress ---*/

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    Weight = element->GetWeight(iGauss);
    Jac_x = element->GetJ_x(iGauss);

    /*--- Initialize the deformation gradient for each Gauss Point ---*/

    for (iVar = 0; iVar < 3; iVar++) {
      for (jVar = 0; jVar < 3; jVar++) {
        F_Mat[iVar][jVar] = 0.0;
        b_Mat[iVar][jVar] = 0.0;
      }
    }

    /*--- Retrieve the values of the gradients of the shape functions for each node ---*/
    /*--- This avoids repeated operations ---*/

    for (iNode = 0; iNode < nNode; iNode++) {

      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
        GradNi_Curr_Mat[iNode][iDim] = element->GetGradNi_x(iNode,iGauss,iDim);
        currentCoord[iNode][iDim] = element->GetCurr_Coord(iNode, iDim);
      }

      /*--- Compute the deformation gradient ---*/

      for (iVar = 0; iVar < nDim; iVar++) {
        for (jVar = 0; jVar < nDim; jVar++) {
          F_Mat[iVar][jVar] += currentCoord[iNode][iVar]*GradNi_Ref_Mat[iNode][jVar];
        }
      }

      /*--- This implies plane strain --> Consider the possible implementation for plane stress --*/
      if (nDim == 2) {
        F_Mat[2][2] = 1.0;
      }

    }

    if (nDim == 2) {
      if (plane_stress) {
        // Compute the value of the term 33 for the deformation gradient
        Compute_Plane_Stress_Term(element, config);
        F_Mat[2][2] = f33;
      }
      else {
        F_Mat[2][2] = 1.0;
      }
    }

    /*--- Determinant of F --> Jacobian of the transformation ---*/

    J_F =   F_Mat[0][0]*F_Mat[1][1]*F_Mat[2][2]+
        F_Mat[0][1]*F_Mat[1][2]*F_Mat[2][0]+
        F_Mat[0][2]*F_Mat[1][0]*F_Mat[2][1]-
        F_Mat[0][2]*F_Mat[1][1]*F_Mat[2][0]-
        F_Mat[1][2]*F_Mat[2][1]*F_Mat[0][0]-
        F_Mat[2][2]*F_Mat[0][1]*F_Mat[1][0];

    /*--- Compute the left Cauchy deformation tensor ---*/

    for (iVar = 0; iVar < 3; iVar++) {
      for (jVar = 0; jVar < 3; jVar++) {
        for (kVar = 0; kVar < 3; kVar++) {
          b_Mat[iVar][jVar] += F_Mat[iVar][kVar]*F_Mat[jVar][kVar];
        }
      }
    }

    /*--- Compute the stress tensor ---*/

    Compute_Stress_Tensor(element, config);


    for (iNode = 0; iNode < nNode; iNode++) {

        /*--- Compute the nodal stress term for each gaussian point and for each node, ---*/
        /*--- and add it to the element structure to be retrieved from the solver      ---*/

      for (iVar = 0; iVar < nDim; iVar++) {
        KAux_t_a[iVar] = 0.0;
        for (jVar = 0; jVar < nDim; jVar++) {
          KAux_t_a[iVar] += Weight * Stress_Tensor[iVar][jVar] * GradNi_Curr_Mat[iNode][jVar] * Jac_x;
        }
      }

      element->Add_Kt_a(KAux_t_a, iNode);

    }

  }

}

void CFEM_NonlinearElasticity::Compute_Averaged_NodalStress(CElement *element, CConfig *config) {

  unsigned short iVar, jVar, kVar;
  unsigned short iGauss, nGauss;
  unsigned short iDim, iNode, nNode;

  su2double Weight, Jac_x;

  element->clearStress();
  element->clearElement();       /*--- Restarts the element: avoids adding over previous results in other elements --*/
  element->ComputeGrad_NonLinear();

  nNode = element->GetnNodes();
  nGauss = element->GetnGaussPoints();

  /*--- Computation of the deformation gradient ---*/

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    Weight = element->GetWeight(iGauss);
    Jac_x = element->GetJ_x(iGauss);

    /*--- Initialize the deformation gradient for each Gauss Point ---*/

    for (iVar = 0; iVar < 3; iVar++) {
      for (jVar = 0; jVar < 3; jVar++) {
        F_Mat[iVar][jVar] = 0.0;
        b_Mat[iVar][jVar] = 0.0;
      }
    }

    /*--- Retrieve the values of the gradients of the shape functions for each node ---*/
    /*--- This avoids repeated operations ---*/

    for (iNode = 0; iNode < nNode; iNode++) {

      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
        GradNi_Curr_Mat[iNode][iDim] = element->GetGradNi_x(iNode,iGauss,iDim);
        currentCoord[iNode][iDim] = element->GetCurr_Coord(iNode, iDim);
      }

      /*--- Compute the deformation gradient ---*/

      for (iVar = 0; iVar < nDim; iVar++) {
        for (jVar = 0; jVar < nDim; jVar++) {
          F_Mat[iVar][jVar] += currentCoord[iNode][iVar]*GradNi_Ref_Mat[iNode][jVar];
        }
      }

      /*--- This implies plane strain --> Consider the possible implementation for plane stress --*/
      if (nDim == 2) {
        F_Mat[2][2] = 1.0;
      }

    }

    if (nDim == 2) {
      if (plane_stress) {
        // Compute the value of the term 33 for the deformation gradient
        Compute_Plane_Stress_Term(element, config);
        F_Mat[2][2] = f33;
      }
      else {
        F_Mat[2][2] = 1.0;
      }
    }

    /*--- Determinant of F --> Jacobian of the transformation ---*/

    J_F =   F_Mat[0][0]*F_Mat[1][1]*F_Mat[2][2]+
        F_Mat[0][1]*F_Mat[1][2]*F_Mat[2][0]+
        F_Mat[0][2]*F_Mat[1][0]*F_Mat[2][1]-
        F_Mat[0][2]*F_Mat[1][1]*F_Mat[2][0]-
        F_Mat[1][2]*F_Mat[2][1]*F_Mat[0][0]-
        F_Mat[2][2]*F_Mat[0][1]*F_Mat[1][0];

    /*--- Compute the left Cauchy deformation tensor ---*/

    for (iVar = 0; iVar < 3; iVar++) {
      for (jVar = 0; jVar < 3; jVar++) {
        for (kVar = 0; kVar < 3; kVar++) {
          b_Mat[iVar][jVar] += F_Mat[iVar][kVar]*F_Mat[jVar][kVar];
        }
      }
    }

    /*--- Compute the stress tensor ---*/

    Compute_Stress_Tensor(element, config);

    for (iNode = 0; iNode < nNode; iNode++) {


        /*--- Compute the nodal stress term for each gaussian point and for each node, ---*/
        /*--- and add it to the element structure to be retrieved from the solver      ---*/

      for (iVar = 0; iVar < nDim; iVar++) {
        KAux_t_a[iVar] = 0.0;
        for (jVar = 0; jVar < nDim; jVar++) {
          KAux_t_a[iVar] += Weight * Stress_Tensor[iVar][jVar] * GradNi_Curr_Mat[iNode][jVar] * Jac_x;
        }
      }

      element->Add_Kt_a(KAux_t_a, iNode);

        /*--- Compute the average nodal stresses for each node ---*/

      element->Add_NodalStress(Stress_Tensor[0][0] * element->GetNi_Extrap(iNode, iGauss), iNode, 0);
      element->Add_NodalStress(Stress_Tensor[1][1] * element->GetNi_Extrap(iNode, iGauss), iNode, 1);
      element->Add_NodalStress(Stress_Tensor[0][1] * element->GetNi_Extrap(iNode, iGauss), iNode, 2);
      if (nDim == 3) {
        element->Add_NodalStress(Stress_Tensor[2][2] * element->GetNi_Extrap(iNode, iGauss), iNode, 3);
        element->Add_NodalStress(Stress_Tensor[0][2] * element->GetNi_Extrap(iNode, iGauss), iNode, 4);
        element->Add_NodalStress(Stress_Tensor[1][2] * element->GetNi_Extrap(iNode, iGauss), iNode, 5);
      }

    }

  }


}


CFEM_NeoHookean_Comp::CFEM_NeoHookean_Comp(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_NonlinearElasticity(val_nDim, val_nVar, config) {


}

CFEM_NeoHookean_Comp::~CFEM_NeoHookean_Comp(void) {

}

void CFEM_NeoHookean_Comp::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

  su2double j_red = 1.0;
  su2double fx = 0.0, fpx = 1.0;
  su2double xkm1 = 1.0, xk = 1.0;
  su2double cte = 0.0;

  unsigned short iNR, nNR;
  su2double NRTOL;

  // Maximum number of iterations and tolerance (relative)
  nNR = 10;
  NRTOL = 1E-25;

  // j_red: reduced jacobian, for the 2x2 submatrix of F
  j_red = F_Mat[0][0] * F_Mat[1][1] - F_Mat[1][0] * F_Mat[0][1];
  // cte: constant term in the NR method
  cte = Lambda*log(j_red) - Mu;

  // f(f33)  = mu*f33^2 + lambda*ln(f33) + (lambda*ln(j_red)-mu) = 0
  // f'(f33) = 2*mu*f33 + lambda/f33

  for (iNR = 0; iNR < nNR; iNR++) {
    fx  = Mu*pow(xk,2.0) + Lambda*log(xk) + cte;
    fpx = 2*Mu*xk + (Lambda / xk);
    xkm1 = xk - fx / fpx;
    if (((xkm1 - xk) / xk) < NRTOL) break;
    xk = xkm1;
  }

  f33 = xkm1;

}

void CFEM_NeoHookean_Comp::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

  su2double Mu_p = 0.0, Lambda_p = 0.0;

  /*--- This can be done in a better way ---*/
  if (J_F != 0.0) {
    Mu_p = (Mu - Lambda*log(J_F))/J_F;
    Lambda_p = Lambda/J_F;
  }

  /*--- Assuming plane strain ---*/

  if (nDim == 2) {
      D_Mat[0][0] = Lambda_p + 2.0 * Mu_p;  D_Mat[0][1] = Lambda_p;          D_Mat[0][2] = 0.0;
      D_Mat[1][0] = Lambda_p;          D_Mat[1][1] = Lambda_p + 2.0 * Mu_p;  D_Mat[1][2] = 0.0;
      D_Mat[2][0] = 0.0;            D_Mat[2][1] = 0.0;            D_Mat[2][2] = Mu_p;
  }
  else if (nDim == 3) {
      D_Mat[0][0] = Lambda_p + 2.0 * Mu_p;  D_Mat[0][1] = Lambda_p;          D_Mat[0][2] = Lambda_p;        D_Mat[0][3] = 0.0;  D_Mat[0][4] = 0.0;  D_Mat[0][5] = 0.0;
      D_Mat[1][0] = Lambda_p;          D_Mat[1][1] = Lambda_p + 2.0 * Mu_p;  D_Mat[1][2] = Lambda_p;        D_Mat[1][3] = 0.0;  D_Mat[1][4] = 0.0;  D_Mat[1][5] = 0.0;
      D_Mat[2][0] = Lambda_p;          D_Mat[2][1] = Lambda_p;          D_Mat[2][2] = Lambda_p + 2.0 * Mu_p;  D_Mat[2][3] = 0.0;  D_Mat[2][4] = 0.0;  D_Mat[2][5] = 0.0;
      D_Mat[3][0] = 0.0;            D_Mat[3][1] = 0.0;            D_Mat[3][2] = 0.0;          D_Mat[3][3] = Mu_p;  D_Mat[3][4] = 0.0;  D_Mat[3][5] = 0.0;
      D_Mat[4][0] = 0.0;            D_Mat[4][1] = 0.0;            D_Mat[4][2] = 0.0;          D_Mat[4][3] = 0.0;  D_Mat[4][4] = Mu_p;  D_Mat[4][5] = 0.0;
      D_Mat[5][0] = 0.0;            D_Mat[5][1] = 0.0;            D_Mat[5][2] = 0.0;          D_Mat[5][3] = 0.0;  D_Mat[5][4] = 0.0;  D_Mat[5][5] = Mu_p;
  }

}

void CFEM_NeoHookean_Comp::Compute_Stress_Tensor(CElement *element, CConfig *config) {

  unsigned short iVar,jVar;
  su2double Mu_J = 0.0, Lambda_J = 0.0;
  su2double dij = 0.0;

  /*--- This can be done in a better way ---*/
  if (J_F != 0.0) {
    Mu_J = Mu/J_F;
    Lambda_J = Lambda/J_F;
  }

  for (iVar = 0; iVar < 3; iVar++) {
    for (jVar = 0; jVar < 3; jVar++) {
      if (iVar == jVar) dij = 1.0;
      else if (iVar != jVar) dij = 0.0;
      Stress_Tensor[iVar][jVar] = Mu_J * (b_Mat[iVar][jVar] - dij) + Lambda_J * log(J_F) * dij;
    }
  }

///  if (plane_stress) {
//  cout << "Deformation gradient (F): " << endl;
//  cout << F_Mat[0][0] << " " << F_Mat[0][1] << " " << F_Mat[0][2] << endl;
//  cout << F_Mat[1][0] << " " << F_Mat[1][1] << " " << F_Mat[1][2] << endl;
//  cout << F_Mat[2][0] << " " << F_Mat[2][1] << " " << F_Mat[2][2] << endl;
//  cout << endl;
//  cout << "Left Cauchy-Green tensor (b): " << endl;
//  cout << b_Mat[0][0] << " " << b_Mat[0][1] << " " << b_Mat[0][2] << endl;
//  cout << b_Mat[1][0] << " " << b_Mat[1][1] << " " << b_Mat[1][2] << endl;
//  cout << b_Mat[2][0] << " " << b_Mat[2][1] << " " << b_Mat[2][2] << endl;
//  cout << endl;
//  cout << "Stress Tensor (sigma): " << endl;
//  cout << Stress_Tensor[0][0] << " " << Stress_Tensor[0][1] << " " << Stress_Tensor[0][2] << endl;
//  cout << Stress_Tensor[1][0] << " " << Stress_Tensor[1][1] << " " << Stress_Tensor[1][2] << endl;
//  cout << Stress_Tensor[2][0] << " " << Stress_Tensor[2][1] << " " << Stress_Tensor[2][2] << endl;
//  cout << endl;
//  }

}

CFEM_NeoHookean_Incomp::CFEM_NeoHookean_Incomp(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_NonlinearElasticity(val_nDim, val_nVar, config) {


}

CFEM_NeoHookean_Incomp::~CFEM_NeoHookean_Incomp(void) {

}

void CFEM_NeoHookean_Incomp::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

}

void CFEM_NeoHookean_Incomp::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

  unsigned short iVar;
  su2double el_P;
  su2double Ib = 0.0, Jft;

  /*--- First invariant of b -> Ib = tr(b) ---*/
  for (iVar = 0; iVar < 3; iVar++) {
    Ib += b_Mat[iVar][iVar];
  }

  /*--- Retrieve element pressure ---*/

  el_P = element->GetElement_Pressure();

  /*--- J^(-5/3) ---*/
  Jft = pow(J_F, -1.666666666666667);

  if (nDim == 2) {

    /*--- Diagonal terms ---*/
    D_Mat[0][0] = 2.0 * Mu * Jft * ((4.0 / 9.0) * Ib - (2.0 / 3.0) * b_Mat[0][0]) - el_P;
    D_Mat[1][1] = 2.0 * Mu * Jft * ((4.0 / 9.0) * Ib - (2.0 / 3.0) * b_Mat[1][1]) - el_P;

    D_Mat[2][2] = (1.0 / 3.0) * Mu * Jft * Ib - el_P;

    /*--- Q2 off diagonal terms (and symmetric) ---*/

    D_Mat[0][1] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];
    D_Mat[1][0] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];

    D_Mat[0][2] = 0.0;
    D_Mat[2][0] = 0.0;

  }
  else if (nDim == 3) {

    /*--- Diagonal terms ---*/
    D_Mat[0][0] = 2.0 * Mu * Jft * ((4.0 / 9.0) * Ib - (2.0 / 3.0) * b_Mat[0][0]) - el_P;
    D_Mat[1][1] = 2.0 * Mu * Jft * ((4.0 / 9.0) * Ib - (2.0 / 3.0) * b_Mat[1][1]) - el_P;
    D_Mat[2][2] = 2.0 * Mu * Jft * ((4.0 / 9.0) * Ib - (2.0 / 3.0) * b_Mat[2][2]) - el_P;

    D_Mat[3][3] = (1.0 / 3.0) * Mu * Jft * Ib - el_P;
    D_Mat[4][4] = (1.0 / 3.0) * Mu * Jft * Ib - el_P;
    D_Mat[5][5] = (1.0 / 3.0) * Mu * Jft * Ib - el_P;

    /*--- Q1 off diagonal terms (and symmetric) ---*/

    D_Mat[0][1] = 2.0 * Mu * Jft * ((1.0 / 9.0) * Ib - (1.0 / 3.0) * b_Mat[0][0] - (1.0 / 3.0) * b_Mat[1][1]) + el_P;
    D_Mat[0][2] = 2.0 * Mu * Jft * ((1.0 / 9.0) * Ib - (1.0 / 3.0) * b_Mat[0][0] - (1.0 / 3.0) * b_Mat[2][2]) + el_P;
    D_Mat[1][2] = 2.0 * Mu * Jft * ((1.0 / 9.0) * Ib - (1.0 / 3.0) * b_Mat[1][1] - (1.0 / 3.0) * b_Mat[2][2]) + el_P;

    D_Mat[1][0] = 2.0 * Mu * Jft * ((1.0 / 9.0) * Ib - (1.0 / 3.0) * b_Mat[0][0] - (1.0 / 3.0) * b_Mat[1][1]) + el_P;
    D_Mat[2][0] = 2.0 * Mu * Jft * ((1.0 / 9.0) * Ib - (1.0 / 3.0) * b_Mat[0][0] - (1.0 / 3.0) * b_Mat[2][2]) + el_P;
    D_Mat[2][1] = 2.0 * Mu * Jft * ((1.0 / 9.0) * Ib - (1.0 / 3.0) * b_Mat[1][1] - (1.0 / 3.0) * b_Mat[2][2]) + el_P;

    /*--- Q2 off diagonal terms (and symmetric) ---*/

    D_Mat[0][3] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];
    D_Mat[1][3] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];
    D_Mat[2][3] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];

    D_Mat[0][4] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][2];
    D_Mat[1][4] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][2];
    D_Mat[2][4] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][2];

    D_Mat[0][5] = (-2.0 / 3.0) * Mu * Jft * b_Mat[1][2];
    D_Mat[1][5] = (-2.0 / 3.0) * Mu * Jft * b_Mat[1][2];
    D_Mat[2][5] = (-2.0 / 3.0) * Mu * Jft * b_Mat[1][2];


    D_Mat[3][0] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];
    D_Mat[3][1] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];
    D_Mat[3][2] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];

    D_Mat[4][0] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][2];
    D_Mat[4][1] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][2];
    D_Mat[4][2] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][2];

    D_Mat[5][0] = (-2.0 / 3.0) * Mu * Jft * b_Mat[1][2];
    D_Mat[5][1] = (-2.0 / 3.0) * Mu * Jft * b_Mat[1][2];
    D_Mat[5][2] = (-2.0 / 3.0) * Mu * Jft * b_Mat[1][2];

    /*--- Q3 off diagonal terms (and symmetric) are all zeros ---*/

    D_Mat[3][4] = 0.0;
    D_Mat[3][5] = 0.0;
    D_Mat[4][5] = 0.0;

    D_Mat[4][3] = 0.0;
    D_Mat[5][3] = 0.0;
    D_Mat[5][4] = 0.0;

  }

}

void CFEM_NeoHookean_Incomp::Compute_Stress_Tensor(CElement *element, CConfig *config) {

  unsigned short iDim,jDim;
  su2double dij = 0.0, el_P;
  su2double Ib = 0.0, Jft;

  /*--- First invariant of b -> Ib = tr(b) ---*/
  for (iDim = 0; iDim < 3; iDim++) {
    Ib += b_Mat[iDim][iDim];
  }

  /*--- Retrieve element pressure ---*/

  el_P = element->GetElement_Pressure();

  /*--- J^(-5/3) ---*/
  Jft = pow(J_F, -1.666666666666667);

  for (iDim = 0; iDim < 3; iDim++) {
    for (jDim = 0; jDim < 3; jDim++) {
      if (iDim == jDim) dij = 1.0;
      else if (iDim != jDim) dij = 0.0;
      Stress_Tensor[iDim][jDim] = Mu *  Jft * (b_Mat[iDim][jDim] - ((1.0 / 3.0) * Ib * dij )) + el_P * dij;
    }
  }


}



