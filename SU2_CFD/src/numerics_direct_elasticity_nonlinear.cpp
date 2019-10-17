/*!
 * \file numerics_direct_elasticity_nonlinear.cpp
 * \brief This file contains the routines for setting the tangent matrix and residual of a FEM nonlinear elastic structural problem.
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

CFEANonlinearElasticity::CFEANonlinearElasticity(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEAElasticity(val_nDim, val_nVar, config) {

  nearly_incompressible = (config->GetMaterialCompressibility() == NEARLY_INCOMPRESSIBLE_MAT);

  unsigned short iVar;

  F_Mat = new su2double *[3];
  b_Mat = new su2double *[3];
  FmT_Mat = new su2double *[3];
  Stress_Tensor = new su2double *[3];
  for (iVar = 0; iVar < 3; iVar++) {
    F_Mat[iVar] = new su2double [3];
    b_Mat[iVar] = new su2double [3];
    FmT_Mat[iVar] = new su2double [3];
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

  J_F = 1.0; J_F_Iso = 1.0;
  f33 = 1.0;

  C10 = Mu/2.0;
  D1  = 2.0/Kappa;

  F_Mat_Iso = NULL;
  b_Mat_Iso = NULL;
  cijkl     = NULL;

  F_Mat_Iso = new su2double *[3];
  b_Mat_Iso = new su2double *[3];
  for (iVar = 0; iVar < 3; iVar++){
    F_Mat_Iso[iVar] = new su2double [3];
    b_Mat_Iso[iVar] = new su2double [3];
  }

  unsigned short jVar, kVar;
  cijkl = new su2double ***[3];
  for (iVar = 0; iVar < 3; iVar++){
    cijkl[iVar] = new su2double **[3];
    for (jVar = 0; jVar < 3; jVar++){
      cijkl[iVar][jVar] = new su2double *[3];
      for (kVar = 0; kVar < 3; kVar++){
        cijkl[iVar][jVar][kVar] = new su2double [3];
      }
    }
  }

  maxwell_stress = config->GetDE_Effects();

  ke_DE         = 0.0;
  nElectric_Field   = 0;
  nDim_Electric_Field = 0;

  EField_Ref_Unit   = NULL;
  EField_Ref_Mod    = NULL;
  EField_Curr_Unit  = NULL;
  if (maxwell_stress == true){

    su2double *Electric_Field_Dir = config->Get_Electric_Field_Dir();
    unsigned short iVar, iDim;
    su2double ref_Efield_mod;

    // Initialize value of the electric field in the reference configuration
    EFieldMod_Ref = 0.0;

    ke_DE = config->GetElectric_Constant(0);

    nDim_Electric_Field = config->GetnDim_Electric_Field();

    if (nDim != nDim_Electric_Field) cout << "DIMENSIONS DON'T AGREE (Fix this)" << endl;

     /*--- DV_Val: Vector to store the value of the design variable. ---*/

    nElectric_Field = config->GetnElectric_Field();

    /*--- We initialize the modulus ---*/
    ref_Efield_mod = 0.0;
    /*--- Normalize the electric field vector ---*/
    for (iDim = 0; iDim < nDim_Electric_Field; iDim++) {
      ref_Efield_mod += Electric_Field_Dir[iDim]*Electric_Field_Dir[iDim];
    }
    ref_Efield_mod = sqrt(ref_Efield_mod);

    if (ref_Efield_mod == 0){
      SU2_MPI::Error("The electric field has not been defined!!!!!", CURRENT_FUNCTION);
    }

    /*--- Initialize pointer for the electric field ---*/
    EField_Ref_Unit = new su2double[nDim_Electric_Field];
    /*--- Assign values to the auxiliary Electric_Field structure ---*/
    for (iDim = 0; iDim < nDim_Electric_Field; iDim++) {
      EField_Ref_Unit[iDim] = Electric_Field_Dir[iDim]/ref_Efield_mod;
    }

    /*--- Auxiliary vector for hosting the electric field modulus in the reference configuration ---*/
    EField_Ref_Mod = new su2double[nElectric_Field];
    for (iVar = 0; iVar < nElectric_Field; iVar++)
      EField_Ref_Mod[iVar] = config->Get_Electric_Field_Mod(iVar);

    /*--- Auxiliary vector for computing the electric field in the current configuration ---*/
    EField_Curr_Unit = new su2double[nDim_Electric_Field];
    for (iDim = 0; iDim < nDim_Electric_Field; iDim++) {
      EField_Curr_Unit[iDim] = 0.0;
    }

    /*--- Auxiliary vector for storing the electric field constant ---*/
    unsigned short nElectric_Constant = config->GetnElectric_Constant();
    ke_DE_i = new su2double[nElectric_Field];

    if (nElectric_Constant == nElectric_Field)
      for (iVar = 0; iVar < nElectric_Field; iVar++) ke_DE_i[iVar] = config->GetElectric_Constant(iVar);
    else
      for (iVar = 0; iVar < nElectric_Field; iVar++) ke_DE_i[iVar] = config->GetElectric_Constant(0);

    switch (config->GetDV_FEA()) {
      case YOUNG_MODULUS:
      case POISSON_RATIO:
      case DENSITY_VAL:
      case DEAD_WEIGHT:
      case ELECTRIC_FIELD:
        break;
      default:
        switch (config->GetDirectDiff()){
          case D_EFIELD:
            SU2_TYPE::SetDerivative(EField_Ref_Mod[config->GetnID_DV()],1.0);
            break;
          default:
            break;
        }

    }

  }

}

CFEANonlinearElasticity::~CFEANonlinearElasticity(void) {

  unsigned short iVar, jVar, kVar;

  for (iVar = 0; iVar < 3; iVar++) {
    delete [] F_Mat[iVar];
    delete [] b_Mat[iVar];
    delete [] FmT_Mat[iVar];
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
  delete [] FmT_Mat;
  delete [] Stress_Tensor;
  delete [] KAux_t_a;
  delete [] KAux_P_ab;
  delete [] currentCoord;

  if (F_Mat_Iso != NULL) {
    for (iVar = 0; iVar < 3; iVar++){
      if (F_Mat_Iso[iVar] != NULL) delete [] F_Mat_Iso[iVar];
    }
    delete [] F_Mat_Iso;
  }
  if (b_Mat_Iso != NULL){
    for (iVar = 0; iVar < 3; iVar++){
      if (b_Mat_Iso[iVar] != NULL) delete [] b_Mat_Iso[iVar];
    }
    delete [] b_Mat_Iso;
  }

  if (cijkl != NULL){
    for (iVar = 0; iVar < 3; iVar++){
      for (jVar = 0; jVar < 3; jVar++){
        for (kVar = 0; kVar < 3;kVar++){
          delete [] cijkl[iVar][jVar][kVar];
        }
        delete [] cijkl[iVar][jVar];
      }
      delete [] cijkl[iVar];
    }
    delete [] cijkl;
  }

  if (EField_Ref_Unit  != NULL)   delete [] EField_Ref_Unit;
  if (EField_Ref_Mod   != NULL)   delete [] EField_Ref_Mod;
  if (EField_Curr_Unit != NULL)   delete [] EField_Curr_Unit;

}


void CFEANonlinearElasticity::Compute_Tangent_Matrix(CElement *element, CConfig *config) {

  unsigned short iVar, jVar, kVar;
  unsigned short iGauss, nGauss;
  unsigned short iNode, jNode, nNode;
  unsigned short iDim, bDim;

  su2double Ks_Aux_ab;

  su2double Weight, Jac_x;

  su2double AuxMatrixKc[3][6];
  su2double AuxMatrixKs[3];

  /*--- TODO: Initialize values for the material model considered ---*/
//  cout << "PROPERTY: " << element->Get_iProp() << " and DV " << element->Get_iDV() << endl;
  SetElement_Properties(element, config);
  if (maxwell_stress) SetElectric_Properties(element, config);
  
  /*--- Register pre-accumulation inputs, material props and nodal coords ---*/
  AD::StartPreacc();
  AD::SetPreaccIn(E);
  AD::SetPreaccIn(Nu);
  AD::SetPreaccIn(Rho_s);
  AD::SetPreaccIn(Rho_s_DL);
  if (maxwell_stress) {
    AD::SetPreaccIn(EFieldMod_Ref);
    AD::SetPreaccIn(ke_DE);
  }
  element->SetPreaccIn_Coords();
  /*--- Recompute Lame parameters as they depend on the material properties ---*/
  Compute_Lame_Parameters();

  /*-----------------------------------------------------------*/

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
    }

    if (nDim == 2) {
      if (plane_stress) {
        // Compute the value of the term 33 for the deformation gradient
        Compute_Plane_Stress_Term(element, config);
        F_Mat[2][2] = f33;
      }
      else { // plane strain
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

    Compute_Stress_Tensor(element, config);
//    if (maxwell_stress) Add_MaxwellStress(element, config);
    Compute_Constitutive_Matrix(element, config);


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

  /*--- Register the stress residual as preaccumulation output ---*/
  element->SetPreaccOut_Kt_a();
  AD::EndPreacc();

}

void CFEANonlinearElasticity::Compute_NodalStress_Term(CElement *element, CConfig *config) {

  unsigned short iVar, jVar, kVar;
  unsigned short iGauss, nGauss;
  unsigned short iNode, nNode;
  unsigned short iDim;

  /*--- TODO: Initialize values for the material model considered ---*/
  SetElement_Properties(element, config);
  if (maxwell_stress) SetElectric_Properties(element, config);
  
  /*--- Register pre-accumulation inputs, material props and nodal coords ---*/
  AD::StartPreacc();
  AD::SetPreaccIn(E);
  AD::SetPreaccIn(Nu);
  AD::SetPreaccIn(Rho_s);
  AD::SetPreaccIn(Rho_s_DL);
  if (maxwell_stress) {
    AD::SetPreaccIn(EFieldMod_Ref);
    AD::SetPreaccIn(ke_DE);
  }
  element->SetPreaccIn_Coords();
  /*--- Recompute Lame parameters as they depend on the material properties ---*/
  Compute_Lame_Parameters();
  
  /*-----------------------------------------------------------*/

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
    }

    if (nDim == 2) {
      if (plane_stress) {
        // Compute the value of the term 33 for the deformation gradient
        Compute_Plane_Stress_Term(element, config);
        F_Mat[2][2] = f33;
      }
      else { // plane strain
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
//    if (maxwell_stress) Add_MaxwellStress(element, config);

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

  /*--- Register the stress residual as preaccumulation output ---*/
  element->SetPreaccOut_Kt_a();
  AD::EndPreacc();

}

void CFEANonlinearElasticity::Add_MaxwellStress(CElement *element, CConfig *config){

//  Adds the Maxwell stress to the output of the stress Sxx, Syy, Szz, SVM...

  unsigned short iDim, jDim;

  su2double E0 = 0.0, E1 = 0.0, E2 = 0.0;
  su2double E0_2 = 0.0, E1_2 = 0.0, E2_2 = 0.0;
  su2double E_2 = 0.0;

  Compute_FmT_Mat();

  for (iDim = 0; iDim < nDim; iDim++){
    EField_Curr_Unit[iDim] = 0.0;
    for (jDim = 0; jDim < nDim; jDim++){
      EField_Curr_Unit[iDim] += FmT_Mat[iDim][jDim] * EField_Ref_Unit[jDim];
    }
  }

  E0 = EFieldMod_Ref*EField_Curr_Unit[0];         E0_2 = pow(E0,2);
  E1 = EFieldMod_Ref*EField_Curr_Unit[1];         E1_2 = pow(E1,2);
  if (nDim == 3) {E2 = EFieldMod_Ref*EField_Curr_Unit[2]; E2_2 = pow(E2,2);}

  E_2 = E0_2+E1_2+E2_2;

  Stress_Tensor[0][0] += ke_DE*(E0_2-0.5*E_2);  Stress_Tensor[0][1] += ke_DE*E0*E1;        Stress_Tensor[0][2] += ke_DE*E0*E2;
  Stress_Tensor[1][0] += ke_DE*E1*E0;        Stress_Tensor[1][1] += ke_DE*(E1_2-0.5*E_2);  Stress_Tensor[1][2] += ke_DE*E1*E2;
  Stress_Tensor[2][0] += ke_DE*E2*E0;        Stress_Tensor[2][1] += ke_DE*E2*E1;        Stress_Tensor[2][2] += ke_DE*(E2_2-0.5*E_2);

}

void CFEANonlinearElasticity::SetElectric_Properties(CElement *element, CConfig *config){

  // Set the modulus of the electric field in the current element
  /*--- These variables are set as preaccumulation inputs in Compute_Tangent_Matrix and
  Compute_NodalStress_Term, if you add variables here be sure to register them in those routines too. ---*/
  EFieldMod_Ref = EField_Ref_Mod[element->Get_iDe()];
  ke_DE = ke_DE_i[element->Get_iDe()];

  switch (config->GetDV_FEA()) {
    case ELECTRIC_FIELD:
      EFieldMod_Ref   = DV_Val[element->Get_iDV()] * EFieldMod_Ref;
      break;
  }

}

void CFEANonlinearElasticity::Compute_FmT_Mat(void) {

  FmT_Mat[0][0] = (F_Mat[1][1]*F_Mat[2][2] - F_Mat[1][2]*F_Mat[2][1]) / J_F;
  FmT_Mat[0][1] = (F_Mat[1][2]*F_Mat[2][0] - F_Mat[2][2]*F_Mat[1][0]) / J_F;
  FmT_Mat[0][2] = (F_Mat[1][0]*F_Mat[2][1] - F_Mat[1][1]*F_Mat[2][0]) / J_F;

  FmT_Mat[1][0] = (F_Mat[0][2]*F_Mat[2][1] - F_Mat[0][1]*F_Mat[2][2]) / J_F;
  FmT_Mat[1][1] = (F_Mat[0][0]*F_Mat[2][2] - F_Mat[2][0]*F_Mat[0][2]) / J_F;
  FmT_Mat[1][2] = (F_Mat[0][1]*F_Mat[2][1] - F_Mat[0][0]*F_Mat[2][0]) / J_F;

  FmT_Mat[2][0] = (F_Mat[0][1]*F_Mat[1][2] - F_Mat[0][2]*F_Mat[1][1]) / J_F;
  FmT_Mat[2][1] = (F_Mat[0][2]*F_Mat[1][0] - F_Mat[0][0]*F_Mat[1][2]) / J_F;
  FmT_Mat[2][2] = (F_Mat[0][0]*F_Mat[1][1] - F_Mat[0][1]*F_Mat[1][0]) / J_F;

}

void CFEANonlinearElasticity::Compute_Isochoric_F_b(void) {

  unsigned short iVar, jVar, kVar;

  J_F_Iso = pow(J_F,-0.333333333333333);

  // Isochoric deformation tensor
  for (iVar = 0; iVar < 3; iVar++){
    for (jVar = 0; jVar < 3; jVar++){
      F_Mat_Iso[iVar][jVar] = F_Mat[iVar][jVar] * J_F_Iso;
    }
  }

  // Isochoric left Cauchy-Green tensor

  for (iVar = 0; iVar < 3; iVar++){
    for (jVar = 0; jVar < 3; jVar++){
      b_Mat_Iso[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < 3; kVar++){
        b_Mat_Iso[iVar][jVar] += F_Mat_Iso[iVar][kVar]*F_Mat_Iso[jVar][kVar];
      }
    }
  }

}

void CFEANonlinearElasticity::Assign_cijkl_D_Mat(void) {

  unsigned short iVar, jVar;

  if (nDim == 2){
    D_Mat[0][0] = cijkl[0][0][0][0];
    D_Mat[1][1] = cijkl[1][1][1][1];

    D_Mat[0][1] = cijkl[0][0][1][1];
    D_Mat[1][0] = cijkl[1][1][0][0];

    D_Mat[0][2] = cijkl[0][0][0][1];
    D_Mat[2][0] = cijkl[1][0][0][0];

    D_Mat[1][2] = cijkl[1][1][0][1];
    D_Mat[2][1] = cijkl[1][0][1][1];

    D_Mat[2][2] = cijkl[0][1][0][1];
  }
  else{
    D_Mat[0][0] = cijkl[0][0][0][0];
    D_Mat[1][1] = cijkl[1][1][1][1];
    D_Mat[2][2] = cijkl[2][2][2][2];
    D_Mat[3][3] = cijkl[0][1][0][1];
    D_Mat[4][4] = cijkl[0][2][0][2];
    D_Mat[5][5] = cijkl[1][2][1][2];

    D_Mat[0][1] = cijkl[0][0][1][1];
    D_Mat[0][2] = cijkl[0][0][2][2];
    D_Mat[0][3] = cijkl[0][0][0][1];
    D_Mat[0][4] = cijkl[0][0][0][2];
    D_Mat[0][5] = cijkl[0][0][1][2];

    D_Mat[1][2] = cijkl[1][1][2][2];
    D_Mat[1][3] = cijkl[1][1][0][1];
    D_Mat[1][4] = cijkl[1][1][0][2];
    D_Mat[1][5] = cijkl[1][1][1][2];

    D_Mat[2][3] = cijkl[2][2][0][1];
    D_Mat[2][4] = cijkl[2][2][0][2];
    D_Mat[2][5] = cijkl[2][2][1][2];

    D_Mat[3][4] = cijkl[0][1][0][2];
    D_Mat[3][5] = cijkl[0][1][1][2];

    D_Mat[4][5] = cijkl[0][2][1][2];

    for (jVar = 0; jVar < 6; jVar++){
      for (iVar = 0; iVar < jVar; iVar++){
        D_Mat[jVar][iVar] = D_Mat[iVar][jVar];
      }
    }

  }

}


void CFEANonlinearElasticity::Compute_Averaged_NodalStress(CElement *element, CConfig *config) {

  unsigned short iVar, jVar, kVar;
  unsigned short iGauss, nGauss;
  unsigned short iDim, iNode, nNode;

  /*--- TODO: Initialize values for the material model considered ---*/
  SetElement_Properties(element, config);
  if (maxwell_stress) SetElectric_Properties(element, config);
  /*-----------------------------------------------------------*/

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
    }

    if (nDim == 2) {
      if (plane_stress) {
        // Compute the value of the term 33 for the deformation gradient
        Compute_Plane_Stress_Term(element, config);
        F_Mat[2][2] = f33;
      }
      else { // plane strain
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
    if (maxwell_stress) Add_MaxwellStress(element, config);

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
                                   CConfig *config) : CFEANonlinearElasticity(val_nDim, val_nVar, config) {


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

}


CFEM_Knowles_NearInc::CFEM_Knowles_NearInc(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEANonlinearElasticity(val_nDim, val_nVar, config) {


  /* -- The formulation adopted for this material model has been described by:
   * --
   * -- Suchocki, C., A Finite Element Implementation of Knowles stored-energy function:
   * -- theory, coding and applications, Archive of Mechanical Engineering, Vol. 58, pp. 319-346 (2011).
   * --
   * -- DOI: 10.2478/v10180-011-0021-7
   */

  Bk = config->GetKnowles_B();
  Nk = config->GetKnowles_N();

  trbbar   = 0.0;
  Ek    = 0.0;
  Pr    = 0.0;


}

CFEM_Knowles_NearInc::~CFEM_Knowles_NearInc(void) {

}

void CFEM_Knowles_NearInc::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

  SU2_MPI::Error("This material model cannot (yet) be used for plane stress.",CURRENT_FUNCTION);

}

void CFEM_Knowles_NearInc::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

  /* -- Suchocki (2011) (full reference in class constructor). ---*/

  /*--- Computation of the tensor cijkl ---*/

  unsigned short iVar, jVar, kVar, lVar;

  for (iVar = 0; iVar < 3; iVar++){
    for (jVar = 0; jVar < 3; jVar++){
      for (kVar = 0; kVar < 3; kVar++){
        for (lVar = 0; lVar < 3; lVar++){
          cijkl[iVar][jVar][kVar][lVar] =
            term1 * ((1.0/2.0)*( deltaij(iVar,kVar)*b_Mat_Iso[jVar][lVar]
                      +deltaij(jVar,lVar)*b_Mat_Iso[iVar][kVar]
                      +deltaij(iVar,lVar)*b_Mat_Iso[jVar][kVar]
                      +deltaij(jVar,kVar)*b_Mat_Iso[iVar][lVar])
                 +(2.0/3.0)*(trbbar*deltaij(iVar,jVar)*deltaij(kVar,lVar)
                      -b_Mat_Iso[iVar][jVar]*deltaij(kVar,lVar)
                      -deltaij(iVar,jVar)*b_Mat_Iso[kVar][lVar]))
             +term2 * ( b_Mat_Iso[iVar][jVar]*b_Mat_Iso[kVar][lVar]
                - trbbar*(b_Mat_Iso[iVar][jVar]*deltaij(kVar,lVar)
                     +deltaij(iVar,jVar)*b_Mat_Iso[kVar][lVar])
                + pow(trbbar,2.0) * deltaij(iVar,jVar) * deltaij(kVar,lVar))
             +Kappa * (2.0 * J_F - 1.0) * deltaij(iVar,jVar) * deltaij(kVar,lVar);

        }
      }
    }
  }

  /*--- Reorganizing the tensor into the matrix D ---*/

  Assign_cijkl_D_Mat();


}

void CFEM_Knowles_NearInc::Compute_Stress_Tensor(CElement *element, CConfig *config) {

  /* -- Suchocki (2011) (full reference in class constructor). ---*/

  unsigned short iVar, jVar;

  /*--- Compute the isochoric deformation gradient Fbar and left Cauchy-Green tensor bbar ---*/
  Compute_Isochoric_F_b();

  trbbar = (b_Mat_Iso[0][0] + b_Mat_Iso[1][1] + b_Mat_Iso[2][2]) / 3.0;
  term1 = (Mu / J_F) * pow((1 + (Bk / Nk) * (3.0 * trbbar - 3.0)), (Nk-1.0));
  term2 = 2.0 * (Mu / J_F) * (Bk * (Nk - 1.0) / Nk) *
      pow((1.0 + (Bk / Nk) * (3.0 * trbbar - 3.0)), (Nk-2.0));

  Ek = Kappa * (2.0 * J_F - 1.0);
  Pr = Kappa * (J_F - 1.0);

  for (iVar = 0; iVar < 3; iVar++){
    for (jVar = 0; jVar < 3; jVar++){
      Stress_Tensor[iVar][jVar] = term1 * (b_Mat_Iso[iVar][jVar] -
                        (deltaij(iVar,jVar) * trbbar))
                    + deltaij(iVar,jVar) * Pr;
    }
  }

}


CFEM_IdealDE::CFEM_IdealDE(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEANonlinearElasticity(val_nDim, val_nVar, config) {

  /* -- The formulation adopted for this material model has been described by:
   * --
   * -- Zhao, X. and Suo, Z., Method to analyze programmable deformation of
   * -- dielectric elastomer layers, Applied Physics Letters 93, 251902 (2008).
   * --
   * -- http://imechanica.org/node/4234
   */

  trbbar   = 0.0;
  Eg    = 0.0;
  Eg23  = 0.0;
  Ek    = 0.0;
  Pr    = 0.0;

}

CFEM_IdealDE::~CFEM_IdealDE(void) {

}

void CFEM_IdealDE::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

  SU2_MPI::Error("This material model cannot (yet) be used for plane stress.",CURRENT_FUNCTION);

}

void CFEM_IdealDE::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

  /* -- Zhao, X. and Suo, Z. (2008) (full reference in class constructor). ---*/

  if (nDim == 2){

    D_Mat[0][0] = Eg23*(b_Mat_Iso[0][0]+trbbar)+Ek;
    D_Mat[1][1] = Eg23*(b_Mat_Iso[1][1]+trbbar)+Ek;

    D_Mat[0][1] = -Eg23*(b_Mat_Iso[0][0]+b_Mat_Iso[1][1]-trbbar)+Ek;
    D_Mat[1][0] = -Eg23*(b_Mat_Iso[0][0]+b_Mat_Iso[1][1]-trbbar)+Ek;

    D_Mat[0][2] = Eg23*b_Mat_Iso[0][1]/2.0;
    D_Mat[2][0] = Eg23*b_Mat_Iso[0][1]/2.0;

    D_Mat[1][2] = Eg23*b_Mat_Iso[0][1]/2.0;
    D_Mat[2][1] = Eg23*b_Mat_Iso[0][1]/2.0;

    D_Mat[2][2] = Eg*(b_Mat_Iso[0][0]+b_Mat_Iso[1][1])/2.0;

  }
  else if (nDim == 3){

  }

}

void CFEM_IdealDE::Compute_Stress_Tensor(CElement *element, CConfig *config) {

  /* -- Zhao, X. and Suo, Z. (2008) (full reference in class constructor). ---*/

  unsigned short iVar, jVar;
  su2double dij = 0.0;

  /*--- Compute the isochoric deformation gradient Fbar and left Cauchy-Green tensor bbar ---*/
  Compute_Isochoric_F_b();

  cout.precision(15);
  // Stress terms

  trbbar = (b_Mat_Iso[0][0] + b_Mat_Iso[1][1] + b_Mat_Iso[2][2]) / 3.0;
  Eg = Mu / J_F;
  Ek = Kappa * (2.0 * J_F - 1.0);
  Pr = Kappa * (J_F - 1.0);
  Eg23 = 2.0 * Eg / 3.0;

  // Stress tensor

  for (iVar = 0; iVar < 3; iVar++){
    for (jVar = 0; jVar < 3; jVar++){
      if (iVar == jVar) dij = 1.0;
      else if (iVar != jVar) dij = 0.0;
      Stress_Tensor[iVar][jVar] = Eg * ( b_Mat_Iso[iVar][jVar] - dij * trbbar) + dij * Pr ;
    }
  }

}

CFEM_DielectricElastomer::CFEM_DielectricElastomer(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEANonlinearElasticity(val_nDim, val_nVar, config) {


}

CFEM_DielectricElastomer::~CFEM_DielectricElastomer(void) {

}

void CFEM_DielectricElastomer::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

}

void CFEM_DielectricElastomer::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

  /*--- This reduces performance by now, but it is temporal ---*/

  if (nDim == 2){
      D_Mat[0][0] = 0.0;  D_Mat[0][1] = 0.0;  D_Mat[0][2] = 0.0;
      D_Mat[1][0] = 0.0;  D_Mat[1][1] = 0.0;  D_Mat[1][2] = 0.0;
      D_Mat[2][0] = 0.0;  D_Mat[2][1] = 0.0;  D_Mat[2][2] = 0.0;
  }
  else if (nDim == 3){
      D_Mat[0][0] = 0.0;  D_Mat[0][1] = 0.0;  D_Mat[0][2] = 0.0;  D_Mat[0][3] = 0.0;  D_Mat[0][4] = 0.0;  D_Mat[0][5] = 0.0;
      D_Mat[1][0] = 0.0;  D_Mat[1][1] = 0.0;  D_Mat[1][2] = 0.0;  D_Mat[1][3] = 0.0;  D_Mat[1][4] = 0.0;  D_Mat[1][5] = 0.0;
      D_Mat[2][0] = 0.0;  D_Mat[2][1] = 0.0;  D_Mat[2][2] = 0.0;  D_Mat[2][3] = 0.0;  D_Mat[2][4] = 0.0;  D_Mat[2][5] = 0.0;
      D_Mat[3][0] = 0.0;  D_Mat[3][1] = 0.0;  D_Mat[3][2] = 0.0;  D_Mat[3][3] = 0.0;  D_Mat[3][4] = 0.0;  D_Mat[3][5] = 0.0;
      D_Mat[4][0] = 0.0;  D_Mat[4][1] = 0.0;  D_Mat[4][2] = 0.0;  D_Mat[4][3] = 0.0;  D_Mat[4][4] = 0.0;  D_Mat[4][5] = 0.0;
      D_Mat[5][0] = 0.0;  D_Mat[5][1] = 0.0;  D_Mat[5][2] = 0.0;  D_Mat[5][3] = 0.0;  D_Mat[5][4] = 0.0;  D_Mat[5][5] = 0.0;
  }


}

void CFEM_DielectricElastomer::Compute_Stress_Tensor(CElement *element, CConfig *config) {

  unsigned short iDim, jDim;

  su2double E0 = 0.0, E1 = 0.0, E2 = 0.0;
  su2double E0_2 = 0.0, E1_2 = 0.0, E2_2 = 0.0;
  su2double E_2 = 0.0;

  Compute_FmT_Mat();

  for (iDim = 0; iDim < nDim; iDim++){
    EField_Curr_Unit[iDim] = 0.0;
    for (jDim = 0; jDim < nDim; jDim++){
      EField_Curr_Unit[iDim] += FmT_Mat[iDim][jDim] * EField_Ref_Unit[jDim];
    }
  }

  E0 = EFieldMod_Ref*EField_Curr_Unit[0];          E0_2 = pow(E0,2);
  E1 = EFieldMod_Ref*EField_Curr_Unit[1];          E1_2 = pow(E1,2);
  if (nDim == 3) {E2 = EFieldMod_Ref*EField_Curr_Unit[2];  E2_2 = pow(E2,2);}

  E_2 = E0_2+E1_2+E2_2;

  Stress_Tensor[0][0] = ke_DE*(E0_2-0.5*E_2);  Stress_Tensor[0][1] = ke_DE*E0*E1;      Stress_Tensor[0][2] = ke_DE*E0*E2;
  Stress_Tensor[1][0] = ke_DE*E1*E0;      Stress_Tensor[1][1] = ke_DE*(E1_2-0.5*E_2);  Stress_Tensor[1][2] = ke_DE*E1*E2;
  Stress_Tensor[2][0] = ke_DE*E2*E0;      Stress_Tensor[2][1] = ke_DE*E2*E1;      Stress_Tensor[2][2] = ke_DE*(E2_2-0.5*E_2);


}

