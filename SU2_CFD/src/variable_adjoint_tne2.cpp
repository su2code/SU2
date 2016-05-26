/*!
 * \file variable_adjoint_tne2.cpp
 * \brief Definition of the solution fields.
 * \author S. Copeland
 * \version 4.0.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/variable_structure.hpp"

CAdjTNE2EulerVariable::CAdjTNE2EulerVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
	Psi = NULL;
	ForceProj_Vector = NULL;
  
}

CAdjTNE2EulerVariable::CAdjTNE2EulerVariable(su2double *val_psirho,
                                             su2double *val_phi,
                                             su2double val_psie,
                                             su2double val_psieve,
                                             unsigned short val_nDim,
                                             unsigned short val_nvar,
                                             CConfig *config) : CVariable(val_nDim,
                                                                          val_nvar,
                                                                          config) {

	unsigned short iDim, iMesh, iSpecies, iVar, nMGSmooth = 0;
  
  /*--- Get properties ---*/
  nSpecies = config->GetnSpecies();
  nDim     = val_nDim;
  
  /*--- Array initialization ---*/
	Psi              = NULL;
	ForceProj_Vector = NULL;
  
	/*--- Allocate residual structures ---*/
  Res_TruncError = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
		Res_TruncError[iVar] = 0.0;
	}
  
  /*--- Only for residual smoothing (multigrid) ---*/
	for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
		nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
  
  if (nMGSmooth > 0) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }
	
	/*--- Allocate undivided laplacian ---*/
	if (config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_CENTERED)
		Undivided_Laplacian = new su2double [nVar];
  
  /*--- Allocate limiter ---*/
	if (config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_UPWIND) {
		Limiter      = new su2double [nVar];
		Solution_Max = new su2double [nVar];
		Solution_Min = new su2double [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Limiter[iVar]      = 0.0;
			Solution_Max[iVar] = 0.0;
			Solution_Min[iVar] = 0.0;
		}
	}
  
  /*--- Allocate and initialize solution ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Solution[iSpecies]     = val_psirho[iSpecies];
    Solution_Old[iSpecies] = val_psirho[iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    Solution[nSpecies+iDim]     = val_phi[iDim];
    Solution_Old[nSpecies+iDim] = val_phi[iDim];
  }
  Solution[nSpecies+nDim]       = val_psie;
  Solution_Old[nSpecies+nDim]   = val_psie;
  Solution[nSpecies+nDim+1]     = val_psieve;
  Solution_Old[nSpecies+nDim+1] = val_psieve;
  
  /*--- Allocate auxiliar vector for sensitivity computation ---*/
	Grad_AuxVar = new su2double [nDim];
	
	/*--- Allocate and initialize projection vector for wall boundary condition ---*/
	ForceProj_Vector = new su2double [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		ForceProj_Vector[iDim] = 0.0;
	
}

CAdjTNE2EulerVariable::CAdjTNE2EulerVariable(su2double *val_solution,
                                             unsigned short val_nDim,
                                             unsigned short val_nvar,
                                             CConfig *config) : CVariable(val_nDim, val_nvar, config) {
	unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Array initialization ---*/
	Psi              = NULL;
	ForceProj_Vector = NULL;
  
	/*--- Allocate residual structures ---*/
  Res_TruncError = new su2double [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
		Res_TruncError[iVar] = 0.0;
	}
  
  /*--- Only for residual smoothing (multigrid) ---*/
	for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
		nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
  
  if (nMGSmooth > 0) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }
  
	/*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
	if (config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_CENTERED)
		Undivided_Laplacian = new su2double [nVar];
	if (config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_UPWIND) {
		Limiter      = new su2double [nVar];
		Solution_Max = new su2double [nVar];
		Solution_Min = new su2double [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Limiter[iVar]      = 0.0;
			Solution_Max[iVar] = 0.0;
			Solution_Min[iVar] = 0.0;
		}
	}
  
	/*--- Solution initialization ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_solution[iVar];
		Solution_Old[iVar] = val_solution[iVar];
	}
  
	/*--- Allocate and initializate solution for dual time strategy ---*/
	if (dual_time) {
		Solution_time_n = new su2double [nVar];
		Solution_time_n1 = new su2double [nVar];
    
		for (iVar = 0; iVar < nVar; iVar++) {
			Solution_time_n[iVar] = val_solution[iVar];
			Solution_time_n1[iVar] = val_solution[iVar];
		}
	}
  
  /*--- Allocate auxiliar vector for sensitivity computation ---*/
	Grad_AuxVar = new su2double [nDim];
	
	/*--- Allocate and initializate projection vector for wall boundary condition ---*/
	ForceProj_Vector = new su2double [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		ForceProj_Vector[iDim] = 0.0;
  
}

CAdjTNE2EulerVariable::~CAdjTNE2EulerVariable(void) {
	
	if (Psi               != NULL) delete [] Psi;
	if (ForceProj_Vector  != NULL) delete [] ForceProj_Vector;
  
}

bool CAdjTNE2EulerVariable::SetPrimVar_Compressible(su2double SharpEdge_Distance,
                                                    bool check,
                                                    CConfig *config) {
	unsigned short iVar;
  bool check_dens = false, RightVol = true;
  
  su2double adj_limit = config->GetAdjointLimit();
  
  check_dens = (fabs(Solution[0]) > adj_limit);
  
  /*--- Check that the adjoint solution is bounded ---*/
  
  if (check_dens) {
    
    /*--- Copy the old solution ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    
    RightVol = false;
    
  }
  
  return RightVol;
  
}


CAdjTNE2NSVariable::CAdjTNE2NSVariable(void) : CAdjTNE2EulerVariable() { }

CAdjTNE2NSVariable::CAdjTNE2NSVariable(su2double *val_psirho,
                                       su2double *val_phi,
                                       su2double val_psie,
                                       su2double val_psieve,
                                       unsigned short val_nDim,
                                       unsigned short val_nvar,
                                       CConfig *config) : CAdjTNE2EulerVariable(val_psirho,
                                                                                val_phi,
                                                                                val_psie,
                                                                                val_psieve,
                                                                                val_nDim,
                                                                                val_nvar,
                                                                                config) {

}

CAdjTNE2NSVariable::CAdjTNE2NSVariable(su2double *val_solution,
                                       unsigned short val_nDim,
                                       unsigned short val_nvar,
                                       CConfig *config) : CAdjTNE2EulerVariable(val_solution,
                                                                                val_nDim,
                                                                                val_nvar,
                                                                                config) {
  
}

CAdjTNE2NSVariable::~CAdjTNE2NSVariable(void) { }
