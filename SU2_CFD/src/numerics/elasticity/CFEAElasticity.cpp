/*!
 * \file CFEAElasticity.cpp
 * \brief Base class for all elasticity problems.
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

#include "../../../include/numerics/elasticity/CFEAElasticity.hpp"
#include "../../../../Common/include/parallelization/omp_structure.hpp"


CFEAElasticity::CFEAElasticity(unsigned short val_nDim, unsigned short val_nVar,
                               const CConfig *config) : CNumerics() {

  nDim = val_nDim;
  nVar = val_nVar;

  bool body_forces = config->GetDeadLoad();  // Body forces (dead loads).
  bool pseudo_static = config->GetPseudoStatic();

  unsigned short iVar;

  /*--- Initialize vector structures for multiple material definition ---*/
  const auto nProp = config->GetnElasticityMat();

  E_i = new su2double[nProp];
  for (iVar = 0; iVar < nProp; iVar++)
    E_i[iVar] = config->GetElasticyMod(iVar);

  Nu_i = new su2double[nProp];
  for (iVar = 0; iVar < nProp; iVar++)
    Nu_i[iVar] = config->GetPoissonRatio(iVar);

  Rho_s_i = new su2double[nProp];     // For inertial effects
  Rho_s_DL_i = new su2double[nProp];  // For dead loads

  for (iVar = 0; iVar < nProp; iVar++) {
    Rho_s_DL_i[iVar] = config->GetMaterialDensity(iVar);
    Rho_s_i[iVar] = pseudo_static ? 0.0 : config->GetMaterialDensity(iVar);
  }

  // Initialization
  E   = E_i[0];
  Nu  = Nu_i[0];
  Rho_s = Rho_s_i[0];
  Rho_s_DL = Rho_s_DL_i[0];

  Compute_Lame_Parameters();

  // Auxiliary vector for body forces (dead load)
  FAux_Dead_Load = nullptr;
  if (body_forces) FAux_Dead_Load = new su2double [nDim];

  plane_stress = (config->GetElas2D_Formulation() == STRUCT_2DFORM::PLANE_STRESS);

  KAux_ab = new su2double* [nDim];
  for (iVar = 0; iVar < nDim; iVar++) {
    KAux_ab[iVar] = new su2double[nDim];
  }

  unsigned short nStrain = (nDim==2) ? DIM_STRAIN_2D : DIM_STRAIN_3D;
  unsigned short nNodes = (nDim==2) ? NNODES_2D : NNODES_3D;

  Ba_Mat = new su2double* [nStrain];
  Bb_Mat = new su2double* [nStrain];
  D_Mat  = new su2double* [nStrain];
  Ni_Vec = new su2double [nNodes];
  GradNi_Ref_Mat = new su2double* [nNodes];
  GradNi_Curr_Mat = new su2double* [nNodes];
  for (iVar = 0; iVar < nStrain; iVar++) {
    Ba_Mat[iVar] = new su2double[nDim];
    Bb_Mat[iVar] = new su2double[nDim];
    D_Mat[iVar] = new su2double[nStrain];
  }
  for (iVar = 0; iVar < nNodes; iVar++) {
    GradNi_Ref_Mat[iVar] = new su2double[nDim];
    GradNi_Curr_Mat[iVar] = new su2double[nDim];
  }

  DV_Val      = nullptr;
  n_DV        = 0;
  switch (config->GetDV_FEA()) {
    case YOUNG_MODULUS:
    case POISSON_RATIO:
    case DENSITY_VAL:
    case DEAD_WEIGHT:
    case ELECTRIC_FIELD:
      ReadDV(config);
      for (unsigned short iDV = 0; iDV < n_DV; iDV++)
      if ((config->GetDirectDiff() == D_YOUNG) ||
          (config->GetDirectDiff() == D_POISSON) ||
          (config->GetDirectDiff() == D_RHO) ||
          (config->GetDirectDiff() == D_RHO_DL) ||
          (config->GetDirectDiff() == D_EFIELD)){
            SU2_TYPE::SetDerivative(DV_Val[config->GetnID_DV()],1.0);
          }
      break;

    default:
      switch (config->GetDirectDiff()){
        case D_YOUNG:
          SU2_TYPE::SetDerivative(E_i[config->GetnID_DV()],1.0);
          break;
        case D_POISSON:
          SU2_TYPE::SetDerivative(Nu_i[config->GetnID_DV()],1.0);
          break;
        case D_RHO:
          SU2_TYPE::SetDerivative(Rho_s_i[config->GetnID_DV()],1.0);
          break;
        case D_RHO_DL:
          SU2_TYPE::SetDerivative(Rho_s_DL_i[config->GetnID_DV()],1.0);
          break;
      }

  }
}

CFEAElasticity::~CFEAElasticity() {

  unsigned short iVar;
  unsigned short nStrain = (nDim==2) ? DIM_STRAIN_2D : DIM_STRAIN_3D;
  unsigned short nNodes = (nDim==2) ? NNODES_2D : NNODES_3D;

  for (iVar = 0; iVar < nDim; iVar++) {
    delete [] KAux_ab[iVar];
  }

  for (iVar = 0; iVar < nStrain; iVar++) {
    delete [] Ba_Mat[iVar];
    delete [] Bb_Mat[iVar];
    delete [] D_Mat[iVar];
  }
  for (iVar = 0; iVar < nNodes; iVar++) {
    delete [] GradNi_Ref_Mat[iVar];
    delete [] GradNi_Curr_Mat[iVar];
  }

  delete [] KAux_ab;
  delete [] Ba_Mat;
  delete [] Bb_Mat;
  delete [] D_Mat;
  delete [] GradNi_Ref_Mat;
  delete [] GradNi_Curr_Mat;

  delete[] DV_Val;

  delete [] FAux_Dead_Load;

  delete [] E_i;
  delete [] Nu_i;
  delete [] Rho_s_i;
  delete [] Rho_s_DL_i;
  delete [] Ni_Vec;
}


void CFEAElasticity::Compute_Mass_Matrix(CElement *element, const CConfig *config) {

  /*--- Initialize values for the material model considered ---*/
  SetElement_Properties(element, config);
  /*-----------------------------------------------------------*/

  unsigned short iGauss, nGauss;
  unsigned short iNode, jNode, nNode;

  su2double Weight, Jac_X, val_Mab;

  /*--- Register pre-accumulation inputs, density and reference coords. ---*/
  AD::StartPreacc();
  AD::SetPreaccIn(Rho_s);
  element->SetPreaccIn_Coords(false);

  element->ClearElement();       /*--- Restarts the element: avoids adding over previous results in other elements --*/
  element->ComputeGrad_Linear(); /*--- Need to compute the gradients to obtain the Jacobian ---*/

  nNode = element->GetnNodes();
  nGauss = element->GetnGaussPoints();

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    Weight = element->GetWeight(iGauss);
    Jac_X = element->GetJ_X(iGauss);      /*--- The mass matrix is computed in the reference configuration ---*/

    /*--- Retrieve the values of the shape functions for each node ---*/
    /*--- This avoids repeated operations ---*/
    for (iNode = 0; iNode < nNode; iNode++) {
      Ni_Vec[iNode] = element->GetNi(iNode,iGauss);
    }

    for (iNode = 0; iNode < nNode; iNode++) {

      /*--- Assumming symmetry ---*/
      for (jNode = iNode; jNode < nNode; jNode++) {

        val_Mab = Weight * Ni_Vec[iNode] * Ni_Vec[jNode] * Jac_X * Rho_s;

        element->Add_Mab(iNode, jNode, val_Mab);
        /*--- Symmetric terms --*/
        if (iNode != jNode) {
          element->Add_Mab(jNode, iNode, val_Mab);
        }

      }

    }

  }

  /*--- Register the mass matrix as preaccumulation output. ---*/
  element->SetPreaccOut_Mab();
  AD::EndPreacc();

}


void CFEAElasticity::Compute_Dead_Load(CElement *element, const CConfig *config) {

  /*--- Initialize values for the material model considered ---*/
  SetElement_Properties(element, config);
  /*-----------------------------------------------------------*/

  /*--- Register pre-accumulation inputs, density and reference coords. ---*/
  AD::StartPreacc();
  AD::SetPreaccIn(Rho_s_DL);
  element->SetPreaccIn_Coords(false);

  unsigned short iGauss, nGauss;
  unsigned short iNode, iDim, nNode;

  su2double Weight, Jac_X;

  /* -- Gravity directionality:
   * -- For 2D problems, we assume the direction for gravity is -y
   * -- For 3D problems, we assume the direction for gravity is -z
   */
  su2double g_force[3] = {0.0,0.0,0.0};

  if (nDim == 2) g_force[1] = -1*STANDARD_GRAVITY;
  else if (nDim == 3) g_force[2] = -1*STANDARD_GRAVITY;

  element->ClearElement();       /*--- Restart the element to avoid adding over previous results. --*/
  element->ComputeGrad_Linear(); /*--- Need to compute the gradients to obtain the Jacobian. ---*/

  nNode = element->GetnNodes();
  nGauss = element->GetnGaussPoints();

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    Weight = element->GetWeight(iGauss);
    Jac_X = element->GetJ_X(iGauss);      /*--- The dead load is computed in the reference configuration ---*/

    /*--- Retrieve the values of the shape functions for each node ---*/
    /*--- This avoids repeated operations ---*/
    for (iNode = 0; iNode < nNode; iNode++) {
      Ni_Vec[iNode] = element->GetNi(iNode,iGauss);
    }

    for (iNode = 0; iNode < nNode; iNode++) {

      for (iDim = 0; iDim < nDim; iDim++) {
        FAux_Dead_Load[iDim] = Weight * Ni_Vec[iNode] * Jac_X * Rho_s_DL * g_force[iDim];
      }

      element->Add_FDL_a(iNode, FAux_Dead_Load);

    }

  }

  /*--- Register the dead load as preaccumulation output. ---*/
  element->SetPreaccOut_FDL_a();
  AD::EndPreacc();

}


void CFEAElasticity::SetElement_Properties(const CElement *element, const CConfig *config) {

  /*--- These variables are set as preaccumulation inputs in Compute_Tangent_Matrix and
  Compute_NodalStress_Term, if you add variables here be sure to register them in those routines too. ---*/
  E   = E_i[element->Get_iProp()];
  Nu  = Nu_i[element->Get_iProp()];
  Rho_s = Rho_s_i[element->Get_iProp()];
  Rho_s_DL = Rho_s_DL_i[element->Get_iProp()];

  switch (config->GetDV_FEA()) {
    case YOUNG_MODULUS:
      E   = DV_Val[element->Get_iDV()] * E;
      break;
    case POISSON_RATIO:
      Nu  = DV_Val[element->Get_iDV()] * Nu;
      break;
    case DENSITY_VAL:
      Rho_s = DV_Val[element->Get_iDV()] * Rho_s;
      break;
    case DEAD_WEIGHT:
      Rho_s_DL = DV_Val[element->Get_iDV()] * Rho_s_DL;
      break;
  }

  Compute_Lame_Parameters();

}


void CFEAElasticity::ReadDV(const CConfig *config) {

  int rank = SU2_MPI::GetRank();
  bool master_node = false;
  SU2_OMP_MASTER
  master_node = (rank == MASTER_NODE);
  END_SU2_OMP_MASTER

  unsigned long index;

  string filename;
  ifstream properties_file;

  /*--- Choose the filename of the design variable ---*/

  string input_name;

  switch (config->GetDV_FEA()) {
    case YOUNG_MODULUS:
      input_name = "dv_young.opt";
      break;
    case POISSON_RATIO:
      input_name = "dv_poisson.opt";
      break;
    case DENSITY_VAL:
    case DEAD_WEIGHT:
      input_name = "dv_density.opt";
      break;
    case ELECTRIC_FIELD:
      input_name = "dv_efield.opt";
      break;
    default:
      input_name = "dv.opt";
      break;
  }

  filename = input_name;

  if (master_node) cout << "Filename: " << filename << "." << endl;

  properties_file.open(filename.data(), ios::in);

  /*--- In case there is no file, all elements get the same property (0) ---*/

  if (properties_file.fail()) {

    if (master_node)
      cout << "There is no design variable file." << endl;

    n_DV   = 1;
    DV_Val = new su2double[n_DV];
    for (unsigned short iDV = 0; iDV < n_DV; iDV++)
      DV_Val[iDV] = 1.0;

  }
  else{

    string text_line;

    /*--- First pass: determine number of design variables ---*/

    unsigned short iDV = 0;

    /*--- Skip the first line: it is the header ---*/

    getline (properties_file, text_line);

    while (getline (properties_file, text_line)) iDV++;

    /*--- Close the restart file ---*/

    properties_file.close();

    n_DV = iDV;
    DV_Val = new su2double[n_DV];

    /*--- Reopen the file (TODO: improve this) ---*/

    properties_file.open(filename.data(), ios::in);

    /*--- Skip the first line: it is the header ---*/

    getline (properties_file, text_line);

    iDV = 0;
    while (getline (properties_file, text_line)) {

      istringstream point_line(text_line);

      point_line >> index >> DV_Val[iDV];

      iDV++;

    }

    /*--- Close the restart file ---*/

    properties_file.close();

  }

}
