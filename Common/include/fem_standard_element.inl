/*!
 * \file fem_standard_element.inl
 * \brief In-Line subroutines of the <i>fem_standard_element.hpp</i> file.
 * \author E. van der Weide
 * \version 6.2.0 "Falcon"
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

inline CFEMStandardElementBase::CFEMStandardElementBase(){}

inline CFEMStandardElementBase::~CFEMStandardElementBase(){}

inline unsigned short CFEMStandardElementBase::GetVTK_Type(void) const {return VTK_Type;}

inline const su2double* CFEMStandardElementBase::GetWeightsIntegration(void) const {return wIntegration.data();}

inline unsigned short CFEMStandardElementBase::GetNIntegration(void) const {return nIntegration;}

inline unsigned short CFEMStandardElementBase::GetOrderExact(void){return orderExact;}

inline CFEMStandardElement::CFEMStandardElement(){}

inline CFEMStandardElement::~CFEMStandardElement(){}

inline CFEMStandardElement::CFEMStandardElement(const CFEMStandardElement &other) : CFEMStandardElementBase(other) {Copy(other);}

inline CFEMStandardElement& CFEMStandardElement::operator=(const CFEMStandardElement &other){Copy(other); return (*this);}

inline su2double* CFEMStandardElement::GetBasisFunctionsIntegration(void){return lagBasisIntegration.data();}

inline const su2double* CFEMStandardElement::GetBasisFunctionsIntegrationTrans(void) const {return lagBasisIntegrationTrans.data();}

inline const su2double* CFEMStandardElement::GetBasisFunctionsSolDOFs(void) const {return lagBasisSolDOFs.data();}

inline su2double* CFEMStandardElement::GetDrBasisFunctionsIntegration(void){return drLagBasisIntegration.data();}

inline su2double* CFEMStandardElement::GetDsBasisFunctionsIntegration(void){return dsLagBasisIntegration.data();}

inline su2double* CFEMStandardElement::GetDtBasisFunctionsIntegration(void){return dtLagBasisIntegration.data();}

inline const su2double* CFEMStandardElement::GetMatVandermondeInv(void) const {return matVandermondeInv.data();}

inline const su2double* CFEMStandardElement::GetMatBasisFunctionsIntegration(void) const {return matBasisIntegration.data();}

inline const su2double* CFEMStandardElement::GetDerMatBasisFunctionsIntTrans(void) const {return matDerBasisIntTrans.data();}

inline const su2double* CFEMStandardElement::GetMatDerBasisFunctionsOwnDOFs(void) const {return matDerBasisOwnDOFs.data();}

inline const su2double* CFEMStandardElement::GetMatDerBasisFunctionsSolDOFs(void) const {return matDerBasisSolDOFs.data();}

inline const su2double* CFEMStandardElement::GetMat2ndDerBasisFunctionsInt(void) const {return mat2ndDerBasisInt.data();}

inline unsigned short* CFEMStandardElement::GetConnFace0(void){return connFace0.data();}

inline unsigned short* CFEMStandardElement::GetConnFace1(void){return connFace1.data();}

inline unsigned short* CFEMStandardElement::GetConnFace2(void){return connFace2.data();}

inline unsigned short* CFEMStandardElement::GetConnFace3(void){return connFace3.data();}

inline unsigned short* CFEMStandardElement::GetConnFace4(void){return connFace4.data();}

inline unsigned short* CFEMStandardElement::GetConnFace5(void){return connFace5.data();}

inline unsigned short CFEMStandardElement::GetNDOFs(void) const {return nDOFs;}

inline unsigned short CFEMStandardElement::GetNPoly(void) const {return nPoly;}

inline unsigned short CFEMStandardElement::GetVTK_Type1(void) const {return VTK_Type1;}

inline unsigned short CFEMStandardElement::GetNSubElemsType1(void) const {return subConn1ForPlotting.size()/GetNDOFsPerSubElem(GetVTK_Type1());}

inline const unsigned short* CFEMStandardElement::GetSubConnType1(void) const {return subConn1ForPlotting.data();}

inline unsigned short CFEMStandardElement::GetVTK_Type2(void) const {return VTK_Type2;}

inline unsigned short CFEMStandardElement::GetNSubElemsType2(void) const {return subConn2ForPlotting.size()/GetNDOFsPerSubElem(GetVTK_Type2());}

inline const unsigned short* CFEMStandardElement::GetSubConnType2(void) const {return subConn2ForPlotting.data();}

inline const vector<su2double>* CFEMStandardElement::GetRDOFs(void) const {return &rDOFs;}

inline const vector<su2double>* CFEMStandardElement::GetSDOFs(void) const {return &sDOFs;}

inline const vector<su2double>* CFEMStandardElement::GetTDOFs(void) const {return &tDOFs;}

inline CFEMStandardInternalFace::CFEMStandardInternalFace(){}

inline CFEMStandardInternalFace::~CFEMStandardInternalFace(){}

inline CFEMStandardInternalFace::CFEMStandardInternalFace(const CFEMStandardInternalFace &other) : CFEMStandardElementBase(other) {Copy(other);}

inline CFEMStandardInternalFace& CFEMStandardInternalFace::operator=(const CFEMStandardInternalFace &other){Copy(other); return (*this);}

inline su2double* CFEMStandardInternalFace::GetDrBasisElemIntegrationSide0(void) {return drLagBasisElemIntegrationSide0.data();}

inline su2double* CFEMStandardInternalFace::GetDrBasisElemIntegrationSide1(void) {return drLagBasisElemIntegrationSide1.data();}

inline su2double* CFEMStandardInternalFace::GetDsBasisElemIntegrationSide0(void) {return dsLagBasisElemIntegrationSide0.data();}

inline su2double* CFEMStandardInternalFace::GetDsBasisElemIntegrationSide1(void) {return dsLagBasisElemIntegrationSide1.data();}

inline su2double* CFEMStandardInternalFace::GetDtBasisElemIntegrationSide0(void) {return dtLagBasisElemIntegrationSide0.data();}

inline su2double* CFEMStandardInternalFace::GetDtBasisElemIntegrationSide1(void) {return dtLagBasisElemIntegrationSide1.data();}

inline const su2double* CFEMStandardInternalFace::GetMatDerBasisElemIntegrationSide0(void) const {return matDerBasisElemIntegrationSide0.data();}

inline const su2double* CFEMStandardInternalFace::GetMatDerBasisElemIntegrationSide1(void) const {return matDerBasisElemIntegrationSide1.data();}

inline const su2double* CFEMStandardInternalFace::GetMatDerBasisElemIntegrationTransposeSide0(void) const {return matDerBasisElemIntegrationTransposeSide0.data();}

inline const su2double* CFEMStandardInternalFace::GetMatDerBasisElemIntegrationTransposeSide1(void) const {return matDerBasisElemIntegrationTransposeSide1.data();}

inline const su2double* CFEMStandardInternalFace::GetBasisFaceIntegrationSide0(void) const {return lagBasisFaceIntegrationSide0.data();}

inline const su2double* CFEMStandardInternalFace::GetBasisFaceIntegrationSide1(void) const {return lagBasisFaceIntegrationSide1.data();}

inline const su2double* CFEMStandardInternalFace::GetBasisFaceIntegrationTransposeSide0(void) const {return lagBasisFaceIntegrationTransposeSide0.data();}

inline const su2double* CFEMStandardInternalFace::GetBasisFaceIntegrationTransposeSide1(void) const {return lagBasisFaceIntegrationTransposeSide1.data();}

inline su2double* CFEMStandardInternalFace::GetDrBasisFaceIntegrationSide0(void) {return drLagBasisFaceIntegrationSide0.data();}

inline su2double* CFEMStandardInternalFace::GetDrBasisFaceIntegrationSide1(void) {return drLagBasisFaceIntegrationSide1.data();}

inline su2double* CFEMStandardInternalFace::GetDsBasisFaceIntegrationSide0(void) {return dsLagBasisFaceIntegrationSide0.data();}

inline su2double* CFEMStandardInternalFace::GetDsBasisFaceIntegrationSide1(void) {return dsLagBasisFaceIntegrationSide1.data();}

inline unsigned short CFEMStandardInternalFace::GetNDOFsElemSide0(void) const {return nDOFsElemSide0;}

inline unsigned short CFEMStandardInternalFace::GetNDOFsElemSide1(void) const {return nDOFsElemSide1;}

inline unsigned short CFEMStandardInternalFace::GetNDOFsFaceSide0(void) const {return nDOFsFaceSide0;}

inline unsigned short CFEMStandardInternalFace::GetNDOFsFaceSide1(void) const {return nDOFsFaceSide1;}

inline su2double CFEMStandardInternalFace::GetPenaltyConstant(void) const {return penaltyConstantFace;}

inline CFEMStandardBoundaryFace::CFEMStandardBoundaryFace(){}

inline CFEMStandardBoundaryFace::~CFEMStandardBoundaryFace(){}

inline CFEMStandardBoundaryFace::CFEMStandardBoundaryFace(const CFEMStandardBoundaryFace &other) : CFEMStandardElementBase(other) {Copy(other);}

inline CFEMStandardBoundaryFace& CFEMStandardBoundaryFace::operator=(const CFEMStandardBoundaryFace &other){Copy(other); return (*this);}

inline const su2double* CFEMStandardBoundaryFace::GetDrBasisElemIntegration(void) const {return drLagBasisElemIntegration.data();}

inline const su2double* CFEMStandardBoundaryFace::GetDsBasisElemIntegration(void) const {return dsLagBasisElemIntegration.data();}

inline const su2double* CFEMStandardBoundaryFace::GetDtBasisElemIntegration(void) const {return dtLagBasisElemIntegration.data();}

inline const su2double* CFEMStandardBoundaryFace::GetMatDerBasisElemIntegration(void) const {return matDerBasisElemIntegration.data();}

inline const su2double* CFEMStandardBoundaryFace::GetMatDerBasisElemIntegrationTranspose(void) const {return matDerBasisElemIntegrationTranspose.data();}

inline const su2double* CFEMStandardBoundaryFace::GetBasisFaceIntegration(void) const {return lagBasisFaceIntegration.data();}

inline const su2double* CFEMStandardBoundaryFace::GetBasisFaceIntegrationTranspose(void) const {return lagBasisFaceIntegrationTranspose.data();}

inline const su2double* CFEMStandardBoundaryFace::GetDrBasisFaceIntegration(void) const {return drLagBasisFaceIntegration.data();}

inline const su2double* CFEMStandardBoundaryFace::GetDsBasisFaceIntegration(void) const {return dsLagBasisFaceIntegration.data();}

inline unsigned short CFEMStandardBoundaryFace::GetNDOFsElem(void) const {return nDOFsElem;}

inline unsigned short CFEMStandardBoundaryFace::GetNDOFsFace(void) const {return nDOFsFace;}

inline unsigned short CFEMStandardBoundaryFace::GetNSubFaces(void) const {return subConnForPlotting.size()/GetNDOFsPerSubFace();}

inline su2double CFEMStandardBoundaryFace::GetPenaltyConstant(void) const {return penaltyConstantFace;}

inline const unsigned short* CFEMStandardBoundaryFace::GetSubFaceConn(void) const {return subConnForPlotting.data();}
