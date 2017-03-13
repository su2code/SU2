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

inline unsigned short FEMStandardElementBaseClass::GetVTK_Type(void) const {return VTK_Type;}

inline const su2double* FEMStandardElementBaseClass::GetWeightsIntegration(void) const {return wIntegration.data();}

inline unsigned short FEMStandardElementBaseClass::GetNIntegration(void) const {return nIntegration;}

inline unsigned short FEMStandardElementBaseClass::GetOrderExact(void){return orderExact;}

inline FEMStandardElementClass::FEMStandardElementClass(){}

inline FEMStandardElementClass::~FEMStandardElementClass(){}

inline FEMStandardElementClass::FEMStandardElementClass(const FEMStandardElementClass &other) : FEMStandardElementBaseClass(other) {Copy(other);}

inline FEMStandardElementClass& FEMStandardElementClass::operator=(const FEMStandardElementClass &other){Copy(other); return (*this);}

inline su2double* FEMStandardElementClass::GetBasisFunctionsIntegration(void){return lagBasisIntegration.data();}

inline const su2double* FEMStandardElementClass::GetBasisFunctionsIntegrationTrans(void) const {return lagBasisIntegrationTrans.data();}

inline const su2double* FEMStandardElementClass::GetBasisFunctionsSolDOFs(void) const {return lagBasisSolDOFs.data();}

inline su2double* FEMStandardElementClass::GetDrBasisFunctionsIntegration(void){return drLagBasisIntegration.data();}

inline su2double* FEMStandardElementClass::GetDsBasisFunctionsIntegration(void){return dsLagBasisIntegration.data();}

inline su2double* FEMStandardElementClass::GetDtBasisFunctionsIntegration(void){return dtLagBasisIntegration.data();}

inline const su2double* FEMStandardElementClass::GetMatVandermondeInv(void) const {return matVandermondeInv.data();}

inline const su2double* FEMStandardElementClass::GetMatBasisFunctionsIntegration(void) const {return matBasisIntegration.data();}

inline const su2double* FEMStandardElementClass::GetDerMatBasisFunctionsIntTrans(void) const {return matDerBasisIntTrans.data();}

inline const su2double* FEMStandardElementClass::GetMatDerBasisFunctionsOwnDOFs(void) const {return matDerBasisOwnDOFs.data();}

inline const su2double* FEMStandardElementClass::GetMatDerBasisFunctionsSolDOFs(void) const {return matDerBasisSolDOFs.data();}

inline unsigned short* FEMStandardElementClass::GetConnFace0(void){return connFace0.data();}

inline unsigned short* FEMStandardElementClass::GetConnFace1(void){return connFace1.data();}

inline unsigned short* FEMStandardElementClass::GetConnFace2(void){return connFace2.data();}

inline unsigned short* FEMStandardElementClass::GetConnFace3(void){return connFace3.data();}

inline unsigned short* FEMStandardElementClass::GetConnFace4(void){return connFace4.data();}

inline unsigned short* FEMStandardElementClass::GetConnFace5(void){return connFace5.data();}

inline unsigned short FEMStandardElementClass::GetNDOFs(void) const {return nDOFs;}

inline unsigned short FEMStandardElementClass::GetNPoly(void) const {return nPoly;}

inline unsigned short FEMStandardElementClass::GetVTK_Type1(void) const {return VTK_Type1;}

inline unsigned short FEMStandardElementClass::GetNSubElemsType1(void) const {return subConn1ForPlotting.size()/GetNDOFsPerSubElem(GetVTK_Type1());}

inline const unsigned short* FEMStandardElementClass::GetSubConnType1(void) const {return subConn1ForPlotting.data();}

inline unsigned short FEMStandardElementClass::GetVTK_Type2(void) const {return VTK_Type2;}

inline unsigned short FEMStandardElementClass::GetNSubElemsType2(void) const {return subConn2ForPlotting.size()/GetNDOFsPerSubElem(GetVTK_Type2());}

inline const unsigned short* FEMStandardElementClass::GetSubConnType2(void) const {return subConn2ForPlotting.data();}

inline const vector<su2double>* FEMStandardElementClass::GetRDOFs(void) const {return &rDOFs;}

inline const vector<su2double>* FEMStandardElementClass::GetSDOFs(void) const {return &sDOFs;}

inline const vector<su2double>* FEMStandardElementClass::GetTDOFs(void) const {return &tDOFs;}

inline FEMStandardInternalFaceClass::FEMStandardInternalFaceClass(){}

inline FEMStandardInternalFaceClass::~FEMStandardInternalFaceClass(){}

inline FEMStandardInternalFaceClass::FEMStandardInternalFaceClass(const FEMStandardInternalFaceClass &other) : FEMStandardElementBaseClass(other) {Copy(other);}

inline FEMStandardInternalFaceClass& FEMStandardInternalFaceClass::operator=(const FEMStandardInternalFaceClass &other){Copy(other); return (*this);}

inline su2double* FEMStandardInternalFaceClass::GetDrBasisElemIntegrationSide0(void) {return drLagBasisElemIntegrationSide0.data();}

inline su2double* FEMStandardInternalFaceClass::GetDrBasisElemIntegrationSide1(void) {return drLagBasisElemIntegrationSide1.data();}

inline su2double* FEMStandardInternalFaceClass::GetDsBasisElemIntegrationSide0(void) {return dsLagBasisElemIntegrationSide0.data();}

inline su2double* FEMStandardInternalFaceClass::GetDsBasisElemIntegrationSide1(void) {return dsLagBasisElemIntegrationSide1.data();}

inline su2double* FEMStandardInternalFaceClass::GetDtBasisElemIntegrationSide0(void) {return dtLagBasisElemIntegrationSide0.data();}

inline su2double* FEMStandardInternalFaceClass::GetDtBasisElemIntegrationSide1(void) {return dtLagBasisElemIntegrationSide1.data();}

inline const su2double* FEMStandardInternalFaceClass::GetMatDerBasisElemIntegrationSide0(void) const {return matDerBasisElemIntegrationSide0.data();}

inline const su2double* FEMStandardInternalFaceClass::GetMatDerBasisElemIntegrationSide1(void) const {return matDerBasisElemIntegrationSide1.data();}

inline const su2double* FEMStandardInternalFaceClass::GetMatDerBasisElemIntegrationTransposeSide0(void) const {return matDerBasisElemIntegrationTransposeSide0.data();}

inline const su2double* FEMStandardInternalFaceClass::GetMatDerBasisElemIntegrationTransposeSide1(void) const {return matDerBasisElemIntegrationTransposeSide1.data();}

inline const su2double* FEMStandardInternalFaceClass::GetBasisFaceIntegrationSide0(void) const {return lagBasisFaceIntegrationSide0.data();}

inline const su2double* FEMStandardInternalFaceClass::GetBasisFaceIntegrationSide1(void) const {return lagBasisFaceIntegrationSide1.data();}

inline const su2double* FEMStandardInternalFaceClass::GetBasisFaceIntegrationTransposeSide0(void) const {return lagBasisFaceIntegrationTransposeSide0.data();}

inline const su2double* FEMStandardInternalFaceClass::GetBasisFaceIntegrationTransposeSide1(void) const {return lagBasisFaceIntegrationTransposeSide1.data();}

inline su2double* FEMStandardInternalFaceClass::GetDrBasisFaceIntegrationSide0(void) {return drLagBasisFaceIntegrationSide0.data();}

inline su2double* FEMStandardInternalFaceClass::GetDrBasisFaceIntegrationSide1(void) {return drLagBasisFaceIntegrationSide1.data();}

inline su2double* FEMStandardInternalFaceClass::GetDsBasisFaceIntegrationSide0(void) {return dsLagBasisFaceIntegrationSide0.data();}

inline su2double* FEMStandardInternalFaceClass::GetDsBasisFaceIntegrationSide1(void) {return dsLagBasisFaceIntegrationSide1.data();}

inline unsigned short FEMStandardInternalFaceClass::GetNDOFsElemSide0(void) const {return nDOFsElemSide0;}

inline unsigned short FEMStandardInternalFaceClass::GetNDOFsElemSide1(void) const {return nDOFsElemSide1;}

inline unsigned short FEMStandardInternalFaceClass::GetNDOFsFaceSide0(void) const {return nDOFsFaceSide0;}

inline unsigned short FEMStandardInternalFaceClass::GetNDOFsFaceSide1(void) const {return nDOFsFaceSide1;}

inline su2double FEMStandardInternalFaceClass::GetPenaltyConstant(void) const {return penaltyConstantFace;}

inline FEMStandardBoundaryFaceClass::FEMStandardBoundaryFaceClass(){}

inline FEMStandardBoundaryFaceClass::~FEMStandardBoundaryFaceClass(){}

inline FEMStandardBoundaryFaceClass::FEMStandardBoundaryFaceClass(const FEMStandardBoundaryFaceClass &other) : FEMStandardElementBaseClass(other) {Copy(other);}

inline FEMStandardBoundaryFaceClass& FEMStandardBoundaryFaceClass::operator=(const FEMStandardBoundaryFaceClass &other){Copy(other); return (*this);}

inline const su2double* FEMStandardBoundaryFaceClass::GetDrBasisElemIntegration(void) const {return drLagBasisElemIntegration.data();}

inline const su2double* FEMStandardBoundaryFaceClass::GetDsBasisElemIntegration(void) const {return dsLagBasisElemIntegration.data();}

inline const su2double* FEMStandardBoundaryFaceClass::GetDtBasisElemIntegration(void) const {return dtLagBasisElemIntegration.data();}

inline const su2double* FEMStandardBoundaryFaceClass::GetMatDerBasisElemIntegration(void) const {return matDerBasisElemIntegration.data();}

inline const su2double* FEMStandardBoundaryFaceClass::GetMatDerBasisElemIntegrationTranspose(void) const {return matDerBasisElemIntegrationTranspose.data();}

inline const su2double* FEMStandardBoundaryFaceClass::GetBasisFaceIntegration(void) const {return lagBasisFaceIntegration.data();}

inline const su2double* FEMStandardBoundaryFaceClass::GetBasisFaceIntegrationTranspose(void) const {return lagBasisFaceIntegrationTranspose.data();}

inline const su2double* FEMStandardBoundaryFaceClass::GetDrBasisFaceIntegration(void) const {return drLagBasisFaceIntegration.data();}

inline const su2double* FEMStandardBoundaryFaceClass::GetDsBasisFaceIntegration(void) const {return dsLagBasisFaceIntegration.data();}

inline unsigned short FEMStandardBoundaryFaceClass::GetNDOFsElem(void) const {return nDOFsElem;}

inline unsigned short FEMStandardBoundaryFaceClass::GetNDOFsFace(void) const {return nDOFsFace;}

inline unsigned short FEMStandardBoundaryFaceClass::GetNSubFaces(void) const {return subConnForPlotting.size()/GetNDOFsPerSubFace();}

inline su2double FEMStandardBoundaryFaceClass::GetPenaltyConstant(void) const {return penaltyConstantFace;}

inline const unsigned short* FEMStandardBoundaryFaceClass::GetSubFaceConn(void) const {return subConnForPlotting.data();}
