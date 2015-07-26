/*!
 * \file variable_fem_elasticity.cpp
 * \brief Definition of the variables for FEM elastic structural problems.
 * \author R. Sanchez
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

CFEM_ElasVariable::CFEM_ElasVariable(void) : CVariable() {

	dynamicFEA = false;
	VonMises_Stress = 0.0;

	nConnectedElements = 0;

	Stress 				= NULL;		// Nodal stress (for output purposes)
	FlowTraction 		= NULL;		// Nodal traction due to the fluid (fsi)
//	Residual_Int 		= NULL;		// Internal component of the residual
	Residual_Ext_Surf 	= NULL;		// Residual component due to external surface forces
	Residual_Ext_Body 	= NULL;		// Residual component due to body forces

}

CFEM_ElasVariable::CFEM_ElasVariable(double *val_fea, unsigned short val_nDim, unsigned short val_nvar, CConfig *config) : CVariable(val_nDim, val_nvar, config) {

	unsigned short iVar, iDim, jDim;
	bool fsi = config->GetFSI_Simulation();
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Nonlinear analysis.
	bool body_forces = false;		// Bool for adding body forces in the future.

	nConnectedElements = 0;
	VonMises_Stress = 0.0;

	dynamicFEA = (config->GetDynamic_Analysis() == DYNAMIC);

	if (nDim == 2) Stress = new double [3];
	else if (nDim == 3) Stress = new double [6];

	if (fsi) FlowTraction = new double [nVar]; else FlowTraction = NULL;

//	if (nonlinear_analysis) Residual_Int = new double [nVar];	else Residual_Int = NULL;
	if (body_forces) Residual_Ext_Body = new double [nVar];	else Residual_Ext_Body = NULL;
	Residual_Ext_Surf = new double [nVar];


}

CFEM_ElasVariable::~CFEM_ElasVariable(void) {

	delete [] Stress;

	if (FlowTraction 		!= NULL) delete [] FlowTraction;
//	if (Residual_Int 		!= NULL) delete [] Residual_Int;
	if (Residual_Ext_Body 	!= NULL) delete [] Residual_Ext_Body;
	if (Residual_Ext_Surf 	!= NULL) delete [] FlowTraction;

}
