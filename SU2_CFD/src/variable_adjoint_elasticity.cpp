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
