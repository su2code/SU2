/*!
 * \file element_structure.inl
 * \brief In-Line subroutines of the <i>element_structure.hpp</i> file.
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
 
#pragma once

inline void CElement::ComputeGrad_Linear(void) { }

inline void CElement::ComputeGrad_NonLinear(void) { }

inline void CElement::OutputGradN_X(CGeometry *geometry, CConfig *config) {  }

inline unsigned short CElement::GetnNodes(void) { return nNodes;}

inline unsigned short CElement::GetnGaussPoints(void) { return nGaussPoints;}

inline void CElement::SetRef_Coord(double val_CoordRef, unsigned short iNode, unsigned short iDim) { RefCoord[iNode][iDim] = val_CoordRef;}

inline void CElement::SetCurr_Coord(double val_CoordCurr, unsigned short iNode, unsigned short iDim) { CurrentCoord[iNode][iDim] = val_CoordCurr;}

inline double CElement::GetRef_Coord(unsigned short iNode, unsigned short iDim) { return RefCoord[iNode][iDim];}

inline double CElement::GetCurr_Coord(unsigned short iNode, unsigned short iDim) { return CurrentCoord[iNode][iDim];}

inline double CElement::GetWeight(unsigned short iGauss) { return GaussWeight[iGauss];}

inline double CElement::GetJ_X(unsigned short iGauss) {return GaussPoint[iGauss]->GetJ_X();}

inline double *CElement::Get_Kab(unsigned short nodeA, unsigned short nodeB){ return Kab[nodeA][nodeB];}
