/*!
 * \file numerics_direct_elasticity.cpp
 * \brief This file contains the routines for setting the tangent matrix and residual of a FEM linear elastic structural problem.
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

CFEAElasticity::CFEAElasticity(void) : CNumerics () {

  E        = 1.0;          /*!< \brief Aux. variable, Young's modulus of elasticity. */
  Nu       = 0.0;          /*!< \brief Aux. variable, Poisson's ratio. */
  Rho_s    = 0.0;          /*!< \brief Aux. variable, Structural density. */
  Rho_s_DL = 0.0;          /*!< \brief Aux. variable, Structural density (for dead loads). */
  Mu       = 0.0;          /*!< \brief Aux. variable, Lame's coeficient. */
  Lambda   = 0.0;          /*!< \brief Aux. variable, Lame's coeficient. */
  Kappa    = 0.0;          /*!< \brief Aux. variable, Compressibility constant. */

  E_i        = NULL;       /*!< \brief Young's modulus of elasticity (list). */
  Nu_i       = NULL;       /*!< \brief Poisson's ratio (list). */
  Rho_s_i    = NULL;       /*!< \brief Structural density (list). */
  Rho_s_DL_i = NULL;       /*!< \brief Structural density (for dead loads) (list). */

  plane_stress = false;    /*!< \brief Checks if we are solving a plane stress case */

  Ba_Mat          = NULL;  /*!< \brief Matrix B for node a - Auxiliary. */
  Bb_Mat          = NULL;  /*!< \brief Matrix B for node b - Auxiliary. */
  Ni_Vec          = NULL;  /*!< \brief Vector of shape functions - Auxiliary. */
  D_Mat           = NULL;  /*!< \brief Constitutive matrix - Auxiliary. */
  KAux_ab         = NULL;  /*!< \brief Node ab stiffness matrix - Auxiliary. */
  GradNi_Ref_Mat  = NULL;  /*!< \brief Gradients of Ni - Auxiliary. */
  GradNi_Curr_Mat = NULL;  /*!< \brief Gradients of Ni - Auxiliary. */

  FAux_Dead_Load  = NULL;  /*!< \brief Auxiliar vector for the dead loads */

  DV_Val = NULL;           /*!< \brief For optimization cases, value of the design variables. */
  n_DV = 0.0;              /*!< \brief For optimization cases, number of design variables. */

}

CFEAElasticity::CFEAElasticity(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  bool body_forces = config->GetDeadLoad();  // Body forces (dead loads).
  bool pseudo_static = config->GetPseudoStatic();

  unsigned short iVar;  // TODO: This needs to be changed; however, the config option is only limited to short...

  // Initialize values
  E = 0.0;  Nu = 0.0;     Rho_s = 0.0; Rho_s_DL = 0.0;
  Mu = 0.0; Lambda = 0.0; Kappa = 0.0;

  /*--- Initialize vector structures for multiple material definition ---*/
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

  Compute_Lame_Parameters();

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

CFEAElasticity::~CFEAElasticity(void) {

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

  if (E_i != NULL) delete [] E_i;
  if (Nu_i != NULL) delete [] Nu_i;
  if (Rho_s_i != NULL) delete [] Rho_s_i;
  if (Rho_s_DL_i != NULL) delete [] Rho_s_DL_i;
  if (Ni_Vec != NULL) delete [] Ni_Vec;
}

void CFEAElasticity::Compute_Mass_Matrix(CElement *element, CConfig *config) {

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

void CFEAElasticity::Compute_Dead_Load(CElement *element, CConfig *config) {

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

  if (nDim == 2) g_force[1] = -1*STANDARD_GRAVITY;
  else if (nDim == 3) g_force[2] = -1*STANDARD_GRAVITY;

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

void CFEAElasticity::SetElement_Properties(CElement *element, CConfig *config) {

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

void CFEAElasticity::ReadDV(CConfig *config) {

  int rank = SU2_MPI::GetRank();
  
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
