/*!
 * \file CFEANonlinearElasticity.cpp
 * \brief This file contains the routines for setting the tangent matrix and
 *        residual of a FEM nonlinear elastic structural problem.
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

#include "../../../include/numerics/elasticity/CFEANonlinearElasticity.hpp"


CFEANonlinearElasticity::CFEANonlinearElasticity(unsigned short val_nDim, unsigned short val_nVar,
                                   const CConfig *config) : CFEAElasticity(val_nDim, val_nVar, config) {

  nearly_incompressible = (config->GetMaterialCompressibility() == STRUCT_COMPRESS::NEARLY_INCOMP);

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

  unsigned short nNodes = (nDim==2) ? NNODES_2D : NNODES_3D;

  currentCoord = new su2double* [nNodes];
  for (iVar = 0; iVar < nNodes; iVar++)
    currentCoord[iVar] = new su2double[nDim];

  J_F = 1.0; J_F_Iso = 1.0;
  f33 = 1.0;

  C10 = Mu/2.0;
  D1  = 2.0/Kappa;

  F_Mat_Iso = nullptr;
  b_Mat_Iso = nullptr;

  F_Mat_Iso = new su2double *[3];
  b_Mat_Iso = new su2double *[3];
  for (iVar = 0; iVar < 3; iVar++){
    F_Mat_Iso[iVar] = new su2double [3];
    b_Mat_Iso[iVar] = new su2double [3];
  }

  maxwell_stress = config->GetDE_Effects();

  ke_DE         = 0.0;
  nElectric_Field   = 0;
  nDim_Electric_Field = 0;

  EField_Ref_Unit   = nullptr;
  EField_Ref_Mod    = nullptr;
  EField_Curr_Unit  = nullptr;

  if (maxwell_stress) {

    const su2double *Electric_Field_Dir = config->Get_Electric_Field_Dir();
    unsigned short iVar, iDim;
    su2double ref_Efield_mod;

    // Initialize value of the electric field in the reference configuration
    EFieldMod_Ref = 0.0;

    ke_DE = config->GetElectric_Constant(0);

    nDim_Electric_Field = config->GetnDim_Electric_Field();

    if (nDim != nDim_Electric_Field) {
      SU2_MPI::Error("The electric field dimensions do not agree with the geometry.", CURRENT_FUNCTION);
    }

    /*--- DV_Val: Vector to store the value of the design variable. ---*/

    nElectric_Field = config->GetnElectric_Field();

    /*--- We initialize the modulus ---*/
    ref_Efield_mod = 0.0;
    /*--- Normalize the electric field vector ---*/
    for (iDim = 0; iDim < nDim_Electric_Field; iDim++) {
      ref_Efield_mod += Electric_Field_Dir[iDim]*Electric_Field_Dir[iDim];
    }
    ref_Efield_mod = sqrt(ref_Efield_mod);

    if (ref_Efield_mod == 0) {
      SU2_MPI::Error("The electric field has not been defined.", CURRENT_FUNCTION);
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
        break;
    }

  }

}

CFEANonlinearElasticity::~CFEANonlinearElasticity() {

  unsigned short iVar;

  for (iVar = 0; iVar < 3; iVar++) {
    delete [] F_Mat[iVar];
    delete [] b_Mat[iVar];
    delete [] FmT_Mat[iVar];
    delete [] Stress_Tensor[iVar];
  }

  for (iVar = 0; iVar < nDim; iVar++) {
    delete [] KAux_P_ab[iVar];
  }

  unsigned short nNodes = (nDim==2) ? NNODES_2D : NNODES_3D;

  for (iVar = 0; iVar < nNodes; iVar++) {
    delete [] currentCoord[iVar];
  }

  delete [] F_Mat;
  delete [] b_Mat;
  delete [] FmT_Mat;
  delete [] Stress_Tensor;
  delete [] KAux_t_a;
  delete [] KAux_P_ab;
  delete [] currentCoord;

  if (F_Mat_Iso != nullptr) {
    for (iVar = 0; iVar < 3; iVar++){
      if (F_Mat_Iso[iVar] != nullptr) delete [] F_Mat_Iso[iVar];
    }
    delete [] F_Mat_Iso;
  }
  if (b_Mat_Iso != nullptr){
    for (iVar = 0; iVar < 3; iVar++){
      if (b_Mat_Iso[iVar] != nullptr) delete [] b_Mat_Iso[iVar];
    }
    delete [] b_Mat_Iso;
  }

  delete [] EField_Ref_Unit;
  delete [] EField_Ref_Mod;
  delete [] EField_Curr_Unit;

}


void CFEANonlinearElasticity::Compute_Tangent_Matrix(CElement *element, const CConfig *config) {

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
  if (maxwell_stress) {
    AD::SetPreaccIn(EFieldMod_Ref);
    AD::SetPreaccIn(ke_DE);
  }
  element->SetPreaccIn_Coords();
  /*--- Recompute Lame parameters as they depend on the material properties ---*/
  Compute_Lame_Parameters();

  /*-----------------------------------------------------------*/

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
      AuxMatrixKc[iVar][jVar] = 0.0;
    }
  }

  for (iVar = 0; iVar < 3; iVar++) {
    AuxMatrixKs[iVar] = 0.0;
  }

  element->ClearElement();  /*--- Restarts the element to avoid adding over previous results. --*/
  element->ComputeGrad_Linear();
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

      element->Add_Kt_a(iNode, KAux_t_a);

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

        element->Add_Kab(iNode, jNode, KAux_ab);
        element->Add_Ks_ab(iNode, jNode, Ks_Aux_ab);
        /*--- Symmetric terms --*/
        if (iNode != jNode) {
          element->Add_Kab_T(jNode, iNode, KAux_ab);
          element->Add_Ks_ab(jNode, iNode, Ks_Aux_ab);
        }

      }

    }

  }

  /*--- Register the stress residual as preaccumulation output ---*/
  element->SetPreaccOut_Kt_a();
  AD::EndPreacc();

}

void CFEANonlinearElasticity::Compute_NodalStress_Term(CElement *element, const CConfig *config) {

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
  if (maxwell_stress) {
    AD::SetPreaccIn(EFieldMod_Ref);
    AD::SetPreaccIn(ke_DE);
  }
  element->SetPreaccIn_Coords();
  /*--- Recompute Lame parameters as they depend on the material properties ---*/
  Compute_Lame_Parameters();

  /*-----------------------------------------------------------*/

  su2double Weight, Jac_x;

  element->ClearElement();       /*--- Restarts the element to avoid adding over previous results. --*/
  element->ComputeGrad_Linear(); /*--- TODO: Check if we can take this out so we don't have to do it twice. ---*/
  element->ComputeGrad_NonLinear();

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

      element->Add_Kt_a(iNode, KAux_t_a);

    }

  }

  /*--- Register the stress residual as preaccumulation output ---*/
  element->SetPreaccOut_Kt_a();
  AD::EndPreacc();

}

void CFEANonlinearElasticity::Add_MaxwellStress(CElement *element, const CConfig *config) {

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

  Stress_Tensor[0][0] += ke_DE*(E0_2-0.5*E_2);  Stress_Tensor[0][1] += ke_DE*E0*E1;           Stress_Tensor[0][2] += ke_DE*E0*E2;
  Stress_Tensor[1][0] += ke_DE*E1*E0;           Stress_Tensor[1][1] += ke_DE*(E1_2-0.5*E_2);  Stress_Tensor[1][2] += ke_DE*E1*E2;
  Stress_Tensor[2][0] += ke_DE*E2*E0;           Stress_Tensor[2][1] += ke_DE*E2*E1;           Stress_Tensor[2][2] += ke_DE*(E2_2-0.5*E_2);

}

void CFEANonlinearElasticity::SetElectric_Properties(const CElement *element, const CConfig *config) {

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

void CFEANonlinearElasticity::Compute_FmT_Mat() {

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

void CFEANonlinearElasticity::Compute_Isochoric_F_b() {

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

void CFEANonlinearElasticity::Assign_cijkl_D_Mat() {

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

su2double CFEANonlinearElasticity::Compute_Averaged_NodalStress(CElement *element, const CConfig *config) {

  unsigned short iVar, jVar, kVar;
  unsigned short iGauss, nGauss;
  unsigned short iDim, iNode, nNode;

  su2double avgStress[DIM_STRAIN_3D] = {0.0};

  /*--- TODO: Initialize values for the material model considered ---*/
  SetElement_Properties(element, config);
  if (maxwell_stress) SetElectric_Properties(element, config);
  /*-----------------------------------------------------------*/

  /*--- Register pre-accumulation inputs ---*/
  /*--- WARNING: Outputs must be registered outside of this method, this allows more
   * flexibility in selecting what is captured by AD, capturing the entire stress
   * tensor would use more memory than that used by the stress residuals. ---*/
  AD::StartPreacc();
  AD::SetPreaccIn(E);
  AD::SetPreaccIn(Nu);
  if (maxwell_stress) {
    AD::SetPreaccIn(EFieldMod_Ref);
    AD::SetPreaccIn(ke_DE);
  }
  element->SetPreaccIn_Coords();
  /*--- Recompute Lame parameters as they depend on the material properties ---*/
  Compute_Lame_Parameters();

  su2double Weight, Jac_x;

  element->ClearStress();
  element->ClearElement(); /*--- Restarts the element to avoid adding over previous results. ---*/
  element->ComputeGrad_Linear();
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
        currentCoord[iNode][iDim] = element->GetCurr_Coord(iNode,iDim);
        GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
        GradNi_Curr_Mat[iNode][iDim] = element->GetGradNi_x(iNode,iGauss,iDim);
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

    J_F = F_Mat[0][0]*F_Mat[1][1]*F_Mat[2][2]+
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

    avgStress[0] += Stress_Tensor[0][0] / nGauss;
    avgStress[1] += Stress_Tensor[1][1] / nGauss;
    avgStress[2] += Stress_Tensor[0][1] / nGauss;
    if (nDim == 3) {
      avgStress[3] += Stress_Tensor[2][2] / nGauss;
      avgStress[4] += Stress_Tensor[0][2] / nGauss;
      avgStress[5] += Stress_Tensor[1][2] / nGauss;
    }

    for (iNode = 0; iNode < nNode; iNode++) {

      /*--- Compute the nodal stress term for each gaussian point and for each node, ---*/
      /*--- and add it to the element structure to be retrieved from the solver      ---*/

      for (iVar = 0; iVar < nDim; iVar++) {
        KAux_t_a[iVar] = 0.0;
        for (jVar = 0; jVar < nDim; jVar++) {
          KAux_t_a[iVar] += Weight * Stress_Tensor[iVar][jVar] * GradNi_Curr_Mat[iNode][jVar] * Jac_x;
        }
      }

      element->Add_Kt_a(iNode, KAux_t_a);

      /*--- Compute the average nodal stresses for each node ---*/

      su2double Ni_Extrap = element->GetNi_Extrap(iNode, iGauss);

      element->Add_NodalStress(iNode, 0, Stress_Tensor[0][0] * Ni_Extrap);
      element->Add_NodalStress(iNode, 1, Stress_Tensor[1][1] * Ni_Extrap);
      element->Add_NodalStress(iNode, 2, Stress_Tensor[0][1] * Ni_Extrap);
      if (nDim == 3) {
        element->Add_NodalStress(iNode, 3, Stress_Tensor[2][2] * Ni_Extrap);
        element->Add_NodalStress(iNode, 4, Stress_Tensor[0][2] * Ni_Extrap);
        element->Add_NodalStress(iNode, 5, Stress_Tensor[1][2] * Ni_Extrap);
      }

    }

  }

  auto elStress = VonMisesStress(nDim, avgStress);

  /*--- We only differentiate w.r.t. an avg VM stress for the element as
   * considering all nodal stresses would use too much memory. ---*/
  AD::SetPreaccOut(elStress);
  AD::EndPreacc();

  return elStress;
}
