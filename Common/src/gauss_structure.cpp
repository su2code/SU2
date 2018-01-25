/*!
 * \file gauss_structure.cpp
 * \brief Definition of the Gaussian Points structure for Finite Element applications
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

#include "../include/gauss_structure.hpp"

CGaussVariable::CGaussVariable(void) {

	GradNi_Xj = NULL;
	GradNi_xj = NULL;
	J_X = 0.0;
	J_x = 0.0;
	iGaussPoint = 0;
	Ni = NULL;

}

CGaussVariable::CGaussVariable(unsigned short val_iGauss, unsigned short val_nDim, unsigned short val_nNodes) {

	 GradNi_Xj = new su2double* [val_nNodes];
	 for (unsigned short iNode = 0; iNode < val_nNodes; iNode++)
		 GradNi_Xj[iNode] = new su2double [val_nDim];

	 GradNi_xj = new su2double* [val_nNodes];
	 for (unsigned short iNode = 0; iNode < val_nNodes; iNode++)
		 GradNi_xj[iNode] = new su2double [val_nDim];

	 J_X = 0.0;
	 J_x = 0.0;

	 iGaussPoint = val_iGauss;

	 Ni = new su2double [val_nNodes];

}

CGaussVariable::~CGaussVariable(void) {

  if (GradNi_Xj            != NULL) delete [] GradNi_Xj;
  if (GradNi_xj            != NULL) delete [] GradNi_xj;
  if (Ni            	   != NULL) delete [] Ni;

}

CElementProperty::CElementProperty(void) {

  iMat_Mod = 0;
  iMat_Prop = 0;
  iElectric_Prop = 0;
  iDV = 0;

}

CElementProperty::CElementProperty(unsigned long valMat_Model, unsigned long valMat_Prop, unsigned long valElectric_Prop, unsigned long valDV) {

  iMat_Mod = valMat_Model;
  iMat_Prop = valMat_Prop;
  iElectric_Prop = valElectric_Prop;
  iDV = valDV;

}

CElementProperty::~CElementProperty(void) {

}
