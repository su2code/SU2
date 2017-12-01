/*!
 * \file numerics_direct_elasticity.cpp
 * \brief This file contains the routines for setting the tangent matrix and residual of a FEM linear elastic structural problem.
 * \author R. Sanchez
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
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

CFEM_Elasticity::CFEM_Elasticity(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  bool body_forces = config->GetDeadLoad();  // Body forces (dead loads).
  bool pseudo_static = config->GetPseudoStatic();

  unsigned short iVar;  // TODO: This needs to be changed; however, the config option is only limited to short...

  // Initialize values
  E = 0.0;  Nu = 0.0;     Rho_s = 0.0; Rho_s_DL = 0.0;
  Mu = 0.0; Lambda = 0.0; Kappa = 0.0;

  /*--- Initialize vector structures for multiple material definition ---*/
  /*--- TODO: This needs to be changed, to adapt for cases in which the numbers in the config are not the same ---*/
  E_i         = new su2double[config->GetnElasticityMod()];
  for (iVar = 0; iVar < config->GetnElasticityMod(); iVar++)
    E_i[iVar]        = config->GetElasticyMod(iVar);

  Nu_i        = new su2double[config->GetnPoissonRatio()];
  for (iVar = 0; iVar < config->GetnPoissonRatio(); iVar++)
    Nu_i[iVar]        = config->GetPoissonRatio(iVar);

  Rho_s_i     = new su2double[config->GetnMaterialDensity()];     // For inertial effects
  Rho_s_DL_i  = new su2double[config->GetnMaterialDensity()];     // For dead loads
  for (iVar = 0; iVar < config->GetnMaterialDensity(); iVar++){
    Rho_s_DL_i[iVar]        = config->GetMaterialDensity(iVar);
    if (pseudo_static) Rho_s_i[iVar] = 0.0;
    else               Rho_s_i[iVar] = config->GetMaterialDensity(iVar);
  }

  if (pseudo_static) Rho_s_i[0] = 0.0;             // Pseudo-static: no inertial effects considered

  // Initialization
  E   = E_i[0];
  Nu  = Nu_i[0];
  Rho_s = Rho_s_i[0];
  Rho_s_DL = Rho_s_DL_i[0];

  Mu     = E / (2.0*(1.0 + Nu));
  Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));
  Kappa  = Lambda + (2/3)*Mu;

  // Auxiliary vector for body forces (dead load)
  if (body_forces) FAux_Dead_Load = new su2double [nDim]; else FAux_Dead_Load = NULL;

  plane_stress = (config->GetElas2D_Formulation() == PLANE_STRESS);

  KAux_ab = new su2double* [nDim];
  for (iVar = 0; iVar < nDim; iVar++) {
    KAux_ab[iVar] = new su2double[nDim];
  }


  if (nDim == 2) {
    Ba_Mat = new su2double* [3];
    Bb_Mat = new su2double* [3];
    D_Mat  = new su2double* [3];
    Ni_Vec  = new su2double [4];      /*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
    GradNi_Ref_Mat = new su2double* [4];  /*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
    GradNi_Curr_Mat = new su2double* [4];  /*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
    for (iVar = 0; iVar < 3; iVar++) {
      Ba_Mat[iVar]      = new su2double[nDim];
      Bb_Mat[iVar]      = new su2double[nDim];
      D_Mat[iVar]       = new su2double[3];
    }
    for (iVar = 0; iVar < 4; iVar++) {
      GradNi_Ref_Mat[iVar]   = new su2double[nDim];
      GradNi_Curr_Mat[iVar]   = new su2double[nDim];
    }
  }
  else if (nDim == 3) {
    Ba_Mat = new su2double* [6];
    Bb_Mat = new su2double* [6];
    D_Mat  = new su2double* [6];
    Ni_Vec  = new su2double [8];      /*--- As of now, 8 is the maximum number of nodes for 3D problems ---*/
    GradNi_Ref_Mat = new su2double* [8];  /*--- As of now, 8 is the maximum number of nodes for 3D problems ---*/
    GradNi_Curr_Mat = new su2double* [8];  /*--- As of now, 8 is the maximum number of nodes for 3D problems ---*/
    for (iVar = 0; iVar < 6; iVar++) {
      Ba_Mat[iVar]      = new su2double[nDim];
      Bb_Mat[iVar]      = new su2double[nDim];
      D_Mat[iVar]        = new su2double[6];
    }
    for (iVar = 0; iVar < 8; iVar++) {
      GradNi_Ref_Mat[iVar]   = new su2double[nDim];
      GradNi_Curr_Mat[iVar]   = new su2double[nDim];
    }
  }

  DV_Val      = NULL;
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

CFEM_Elasticity::~CFEM_Elasticity(void) {

  unsigned short iVar;

  for (iVar = 0; iVar < nDim; iVar++) {
    delete [] KAux_ab[iVar];
  }

  if (nDim == 2) {
    for (iVar = 0; iVar < 3; iVar++) {
      delete [] Ba_Mat[iVar];
      delete [] Bb_Mat[iVar];
      delete [] D_Mat[iVar];
    }
    for (iVar = 0; iVar < 4; iVar++) {
      delete [] GradNi_Ref_Mat[iVar];
      delete [] GradNi_Curr_Mat[iVar];
    }
  }
  else if (nDim == 3) {
    for (iVar = 0; iVar < 6; iVar++) {
      delete [] Ba_Mat[iVar];
      delete [] Bb_Mat[iVar];
      delete [] D_Mat[iVar];
    }
    for (iVar = 0; iVar < 8; iVar++) {
      delete [] GradNi_Ref_Mat[iVar];
      delete [] GradNi_Curr_Mat[iVar];
    }
  }

  delete [] KAux_ab;
  delete [] Ba_Mat;
  delete [] Bb_Mat;
  delete [] D_Mat;
  delete [] GradNi_Ref_Mat;
  delete [] GradNi_Curr_Mat;

  if (DV_Val != NULL) delete[] DV_Val;

  if (FAux_Dead_Load     != NULL) delete [] FAux_Dead_Load;

}

void CFEM_Elasticity::Compute_Mass_Matrix(CElement *element, CConfig *config) {

  /*--- Initialize values for the material model considered ---*/
  SetElement_Properties(element, config);
  /*-----------------------------------------------------------*/

  unsigned short iGauss, nGauss;
  unsigned short iNode, jNode, nNode;

  su2double Weight, Jac_X;

  su2double val_Mab;

  element->clearElement();       /*--- Restarts the element: avoids adding over previous results in other elements --*/
  element->ComputeGrad_Linear();    /*--- Need to compute the gradients to obtain the Jacobian ---*/

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

        element->Add_Mab(val_Mab,iNode, jNode);
        /*--- Symmetric terms --*/
        if (iNode != jNode) {
          element->Add_Mab(val_Mab, jNode, iNode);
        }

      }

    }

  }

}

void CFEM_Elasticity::Compute_Dead_Load(CElement *element, CConfig *config) {

  /*--- Initialize values for the material model considered ---*/
  SetElement_Properties(element, config);
  /*-----------------------------------------------------------*/

  unsigned short iGauss, nGauss;
  unsigned short iNode, iDim, nNode;

  su2double Weight, Jac_X;

  /* -- Gravity directionality:
   * -- For 2D problems, we assume the direction for gravity is -y
   * -- For 3D problems, we assume the direction for gravity is -z
   */
  su2double g_force[3] = {0.0,0.0,0.0};

  if (nDim == 2) g_force[1] = -1*STANDART_GRAVITY;
  else if (nDim == 3) g_force[2] = -1*STANDART_GRAVITY;

  element->clearElement();       /*--- Restarts the element: avoids adding over previous results in other elements and sets initial values to 0--*/
  element->ComputeGrad_Linear();    /*--- Need to compute the gradients to obtain the Jacobian ---*/

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

      element->Add_FDL_a(FAux_Dead_Load,iNode);

    }

  }

}

void CFEM_Elasticity::SetElement_Properties(CElement *element, CConfig *config) {

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

  Mu     = E / (2.0*(1.0 + Nu));
  Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));
  Kappa  = Lambda + (2/3)*Mu;

//  cout << "YOUNG: " << E << " and Mu " << Mu << endl;

}

void CFEM_Elasticity::ReadDV(CConfig *config) {

  unsigned long index;

  string filename;
  ifstream properties_file;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

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

  if (rank == MASTER_NODE) cout << "Filename: " << filename << "." << endl;

  properties_file.open(filename.data(), ios::in);

  /*--- In case there is no file, all elements get the same property (0) ---*/

  if (properties_file.fail()) {

    if (rank == MASTER_NODE)
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
