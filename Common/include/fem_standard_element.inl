/*!
 * \file fem_standard_element.inl
 * \brief In-Line subroutines of the <i>fem_standard_element.hpp</i> file.
 * \author E. van der Weide
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

inline FEMStandardElementBaseClass::FEMStandardElementBaseClass(){}

inline FEMStandardElementBaseClass::~FEMStandardElementBaseClass(){}

inline unsigned short FEMStandardElementBaseClass::GetNIntegration(void){return nIntegration;}

inline unsigned short FEMStandardElementBaseClass::GetOrderExact(void){return orderExact;}

inline FEMStandardElementClass::FEMStandardElementClass(){}

inline FEMStandardElementClass::~FEMStandardElementClass(){}

inline FEMStandardElementClass::FEMStandardElementClass(const FEMStandardElementClass &other){Copy(other);}

inline FEMStandardElementClass& FEMStandardElementClass::operator=(const FEMStandardElementClass &other){Copy(other); return (*this);}

inline su2double* FEMStandardElementClass::GetDrBasisFunctionsIntegration(void){return drLagBasisIntegration.data();}

inline su2double* FEMStandardElementClass::GetDsBasisFunctionsIntegration(void){return dsLagBasisIntegration.data();}

inline su2double* FEMStandardElementClass::GetDtBasisFunctionsIntegration(void){return dtLagBasisIntegration.data();}

inline unsigned short* FEMStandardElementClass::GetConnFace0(void){return connFace0.data();}

inline unsigned short* FEMStandardElementClass::GetConnFace1(void){return connFace1.data();}

inline unsigned short* FEMStandardElementClass::GetConnFace2(void){return connFace2.data();}

inline unsigned short* FEMStandardElementClass::GetConnFace3(void){return connFace3.data();}

inline unsigned short* FEMStandardElementClass::GetConnFace4(void){return connFace4.data();}

inline unsigned short* FEMStandardElementClass::GetConnFace5(void){return connFace5.data();}

inline unsigned short FEMStandardElementClass::GetNDOFs(void){return nDOFs;}
