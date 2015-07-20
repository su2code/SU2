/*!
 * \file element_structure.cpp
 * \brief Definition of the finite element structure (elements)
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


#include "../include/element_structure.hpp"

unsigned short CElement::nDim = 0;

CElement::CElement(void) {

	CurrentCoord = NULL;
	RefCoord = NULL;
	GaussWeight = NULL;
	GaussCoord = NULL;

	GaussPoint = NULL;

	nNodes = 0;
	nGaussPoints = 0;

}


CElement::CElement(unsigned short val_nDim, CConfig *config) {

	/*--- Initializate the number of dimension and some structures we need ---*/
	nDim = val_nDim;

	CurrentCoord = NULL;
	RefCoord = NULL;
	GaussWeight = NULL;
	GaussCoord = NULL;

}

CElement::~CElement(void) {
	unsigned short iVar;

}

double CElement::GetGradNi_X(unsigned short iNode, unsigned short iGauss, unsigned short iDim){

	return GaussPoint[iGauss]->GetGradNi_Xj(iNode,iDim);

}

void CElement::Add_Kab(double **val_Kab, unsigned short nodeA, unsigned short nodeB){

	unsigned short iDim, jDim;

	for(iDim = 0; iDim < nDim; iDim++) {
		for (jDim = 0; jDim < nDim; jDim++) {
			Kab[nodeA][nodeB][iDim*nDim+jDim] += val_Kab[iDim][jDim];
		}
	}
}

void CElement::Add_Kab_T(double **val_Kab, unsigned short nodeA, unsigned short nodeB){

	unsigned short iDim, jDim;

	for(iDim = 0; iDim < nDim; iDim++) {
		for (jDim = 0; jDim < nDim; jDim++) {
			Kab[nodeA][nodeB][iDim*nDim+jDim] += val_Kab[jDim][iDim];
		}
	}
}

void CElement::clearElement(void){

	unsigned short iNode, jNode, iDim, nDimSq;

	nDimSq = nDim*nDim;

	for(iNode = 0; iNode < nNodes; iNode++) {
		for (jNode = 0; jNode < nNodes; jNode++) {
			for(iDim = 0; iDim < nDimSq; iDim++){
				Kab[iNode][jNode][iDim] = 0.0;
			}
		}
	}
}

