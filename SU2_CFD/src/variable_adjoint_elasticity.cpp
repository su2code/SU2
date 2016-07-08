/*!
 * \file variable_adjoint_elasticity.cpp
 * \brief Definition of the variables for FEM adjoint elastic structural problems.
 * \author R. Sanchez
 * \version 4.1.0 "Cardinal"
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
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

CFEM_ElasVariable_Adj::CFEM_ElasVariable_Adj(void) : CVariable() {

	Reference_Geometry		= NULL;		// Reference geometry for optimization purposes
	Gradient_Adj			= NULL;		// Adjoint gradient dS/dv for structural problems (temporary)

}

CFEM_ElasVariable_Adj::CFEM_ElasVariable_Adj(su2double *val_fea, unsigned short val_nDim, unsigned short val_nvar, CConfig *config) : CVariable(val_nDim, val_nvar, config) {

	unsigned short iVar;
	bool refgeom = config->GetRefGeom();				// Reference geometry needs to be stored

	/*--- Initialization of variables ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_fea[iVar];
	}

	Gradient_Adj    = new su2double [nVar];

    switch (config->GetKind_ObjFunc()) {
    	case REFERENCE_GEOMETRY:
    		Reference_Geometry = new su2double [nVar];
    		break;
    	default:
    		Reference_Geometry = NULL;
    		break;
    }


}

CFEM_ElasVariable_Adj::~CFEM_ElasVariable_Adj(void) {

	if (Reference_Geometry 		!= NULL) delete [] Reference_Geometry;
	if (Gradient_Adj 			!= NULL) delete [] Gradient_Adj;

}


CDiscAdjFEAVariable::CDiscAdjFEAVariable() : CVariable(){

  Sensitivity           = NULL;
  Solution_Direct       = NULL;
  DualTime_Derivative   = NULL;
  DualTime_Derivative_n = NULL;

  Solution_Direct_Vel   = NULL;
  Solution_Direct_Accel = NULL;

  Solution_Vel   = NULL;
  Solution_Accel = NULL;

  Solution_Old_Vel      = NULL;
  Solution_Old_Accel    = NULL;

  Solution_time_n   = NULL;
  Solution_Vel_time_n   = NULL;
  Solution_Accel_time_n = NULL;

}

CDiscAdjFEAVariable::CDiscAdjFEAVariable(su2double* val_solution, unsigned short val_ndim,
                               unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config){

  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

  DualTime_Derivative   = NULL;
  DualTime_Derivative_n = NULL;

  Solution_Direct_Vel   = NULL;
  Solution_Direct_Accel = NULL;

  Solution_Vel          = NULL;
  Solution_Accel        = NULL;

  Solution_Old_Vel      = NULL;
  Solution_Old_Accel    = NULL;

  Solution_time_n   = NULL;
  Solution_Vel_time_n   = NULL;
  Solution_Accel_time_n = NULL;
//  if (dynamic){
//    DualTime_Derivative     = new su2double[nVar];
//    DualTime_Derivative_n   = new su2double[nVar];
//
//    Solution_Direct_Vel     = new su2double[nVar];
//    Solution_Direct_Accel   = new su2double[nVar];
//
//    Solution_Vel            = new su2double[nVar];
//    Solution_Accel          = new su2double[nVar];
//
//    Solution_Old_Vel        = new su2double[nVar];
//    Solution_Old_Accel      = new su2double[nVar];
//
//    Solution_Vel_time_n     = new su2double[nVar];
//    Solution_Accel_time_n   = new su2double[nVar];
//  }

  Solution_Direct = new su2double[nVar];

  Sensitivity = new su2double[nDim];

  unsigned short iVar,iDim;

  for (iDim = 0; iDim < nDim; iDim++){
    Sensitivity[iDim] = 0.0;
  }

  for (iVar = 0; iVar < nVar; iVar++){
    Solution[iVar] = val_solution[iVar];
  }


//  if (dynamic){
//    for (iVar = 0; iVar < nVar; iVar++){
//      Solution_time_n[iVar]         = 0.0;
//
//      DualTime_Derivative[iVar]     = 0.0;
//      DualTime_Derivative_n[iVar]   = 0.0;
//
//      Solution_Direct_Vel[iVar]     = 0.0;
//      Solution_Direct_Accel[iVar]   = 0.0;
//
//      Solution_Vel[iVar]            = 0.0;
//      Solution_Accel[iVar]          = 0.0;
//
//      Solution_Vel_time_n[iVar]     = 0.0;
//      Solution_Accel_time_n[iVar]   = 0.0;
//
//      Solution_Old_Vel[iVar]        = 0.0;
//      Solution_Old_Accel[iVar]      = 0.0;
//
//    }
//  }

}

CDiscAdjFEAVariable::CDiscAdjFEAVariable(su2double* val_solution, su2double* val_solution_accel, su2double* val_solution_vel, unsigned short val_ndim,
                               unsigned short val_nvar, CConfig *config) : CVariable(val_ndim, val_nvar, config){

  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

  DualTime_Derivative     = new su2double[nVar];
  DualTime_Derivative_n   = new su2double[nVar];

  Solution_Direct_Vel     = new su2double[nVar];
  Solution_Direct_Accel   = new su2double[nVar];

  Solution_Vel            = new su2double[nVar];
  Solution_Accel          = new su2double[nVar];

  Solution_Old_Vel        = new su2double[nVar];
  Solution_Old_Accel      = new su2double[nVar];

  Solution_time_n         = new su2double[nVar];
  Solution_Vel_time_n     = new su2double[nVar];
  Solution_Accel_time_n   = new su2double[nVar];

  Solution_Direct = new su2double[nVar];

  Sensitivity = new su2double[nDim];

  unsigned short iVar,iDim;

  for (iDim = 0; iDim < nDim; iDim++){
    Sensitivity[iDim] = 0.0;
  }

  for (iVar = 0; iVar < nVar; iVar++){
    Solution[iVar] = val_solution[iVar];
  }

  for (iVar = 0; iVar < nVar; iVar++){
    Solution_Accel[iVar] = val_solution_accel[iVar];
  }

  for (iVar = 0; iVar < nVar; iVar++){
    Solution_Vel[iVar] = val_solution_vel[iVar];
  }

  /*--- Initialize the rest to 0 ---*/

  for (iVar = 0; iVar < nVar; iVar++){
    Solution_time_n[iVar]         = 0.0;

    DualTime_Derivative[iVar]     = 0.0;
    DualTime_Derivative_n[iVar]   = 0.0;

    Solution_Direct_Vel[iVar]     = 0.0;
    Solution_Direct_Accel[iVar]   = 0.0;

    Solution_Vel_time_n[iVar]     = 0.0;
    Solution_Accel_time_n[iVar]   = 0.0;

    Solution_Old_Vel[iVar]        = 0.0;
    Solution_Old_Accel[iVar]      = 0.0;

  }

}


CDiscAdjFEAVariable::~CDiscAdjFEAVariable(){

  if (Sensitivity           != NULL) delete [] Sensitivity;
  if (Solution_Direct       != NULL) delete [] Solution_Direct;

  if (DualTime_Derivative   != NULL) delete [] DualTime_Derivative;
  if (DualTime_Derivative_n != NULL) delete [] DualTime_Derivative_n;

  if (Solution_Direct_Vel   != NULL) delete [] Solution_Direct_Vel;
  if (Solution_Direct_Accel != NULL) delete [] Solution_Direct_Accel;

  if (Solution_Vel          != NULL) delete [] Solution_Vel;
  if (Solution_Accel        != NULL) delete [] Solution_Accel;

  if (Solution_time_n       != NULL) delete [] Solution_time_n;
  if (Solution_Vel_time_n   != NULL) delete [] Solution_Vel_time_n;
  if (Solution_Accel_time_n != NULL) delete [] Solution_Accel_time_n;

  if (Solution_Old_Vel      != NULL) delete [] Solution_Old_Vel;
  if (Solution_Old_Accel    != NULL) delete [] Solution_Old_Accel;

}

