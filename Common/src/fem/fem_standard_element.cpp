/*!
 * \file fem_standard_element.cpp
 * \brief Functions for the FEM standard elements.
 * \author E. van der Weide
 * \version 7.0.5 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fem/fem_standard_element.hpp"
#include "../../include/fem/fem_gauss_jacobi_quadrature.hpp"
#include "../../include/blas_structure.hpp"

/*----------------------------------------------------------------------------------*/
/*          Public member functions of CFEMStandardElementBase.                     */
/*----------------------------------------------------------------------------------*/

unsigned short CFEMStandardElementBase::GetNDOFsStatic(unsigned short VTK_Type,
                                                       unsigned short nPoly,
                                                       unsigned long  typeErrorMessage) {
  unsigned short nDOFsEdge = nPoly + 1;
  unsigned short nDOFs = 0;    // To avoid a compiler warning.

  switch(VTK_Type) {

    case LINE:
      nDOFs = nDOFsEdge;
      break;

    case TRIANGLE:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)/2;
      break;

    case QUADRILATERAL:
      nDOFs = nDOFsEdge*nDOFsEdge;
      break;

    case TETRAHEDRON:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)*(nDOFsEdge+2)/6;
      break;

    case HEXAHEDRON:
      nDOFs = nDOFsEdge*nDOFsEdge*nDOFsEdge;
      break;

    case PRISM:
      nDOFs = nDOFsEdge*nDOFsEdge*(nDOFsEdge+1)/2;
      break;

    case PYRAMID:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)*(2*nDOFsEdge+1)/6;
      break;

    default:
      ostringstream message;
      message << "Unknown FEM element type, " << typeErrorMessage
              << ", encountered.";
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  return nDOFs;
}


void CFEMStandardElementBase::InverseMatrix(unsigned short    n,
                                            vector<su2double> &A) {

}

void CFEMStandardElementBase::Vandermonde1D(unsigned short          nDOFs,
                                            const vector<su2double> &r,
                                            vector<su2double>       &V) {

}

void CFEMStandardElementBase::GradVandermonde1D(unsigned short          nDOFs,
                                                const vector<su2double> &r,
                                                vector<su2double>       &VDr) {

}

/*----------------------------------------------------------------------------------*/
/*         Protected member functions of CFEMStandardElementBase.                   */
/*----------------------------------------------------------------------------------*/

CFEMStandardElementBase::CFEMStandardElementBase(unsigned short val_VTK_Type,
                                                 unsigned short val_nPoly,
                                                 bool           val_constJac,
                                                 CConfig        *config,
                                                 unsigned short val_orderExact) {

}

void CFEMStandardElementBase::CheckSumDerivativesLagrangianBasisFunctions(
                                          const unsigned short    nPoints,
                                          const unsigned short    nDOFs,
                                          const vector<su2double> &dLagBasisPoints) {

}

void CFEMStandardElementBase::CheckSumLagrangianBasisFunctions(
                                          const unsigned short nPoints,
                                          const unsigned short nDOFs,
                                          vector<su2double>    &lagBasisPoints) {

}

void CFEMStandardElementBase::Copy(const CFEMStandardElementBase &other) {

}

void CFEMStandardElementBase::DerivativesBasisFunctionsAdjacentElement(
                                        unsigned short    VTK_TypeElem,
                                        unsigned short    nPolyElem,
                                        const bool        swapFaceInElement,
                                        unsigned short    &nDOFsElem,
                                        vector<su2double> &drLagBasisIntegration,
                                        vector<su2double> &dsLagBasisIntegration,
                                        vector<su2double> &dtLagBasisIntegration) {

}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesLine(
                                       const unsigned short    nPoly,
                                       const vector<su2double> &rPoints,
                                       unsigned short          &nDOFs,
                                       vector<su2double>       &rDOFs,
                                       vector<su2double>       &matVandermondeInv,
                                       vector<su2double>       &lagBasisPoints,
                                       vector<su2double>       &drLagBasisPoints) {

}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesTriangle(
                                       const unsigned short    nPoly,
                                       const vector<su2double> &rPoints,
                                       const vector<su2double> &sPoints,
                                       unsigned short          &nDOFs,
                                       vector<su2double>       &rDOFs,
                                       vector<su2double>       &sDOFs,
                                       vector<su2double>       &matVandermondeInv,
                                       vector<su2double>       &lagBasisPoints,
                                       vector<su2double>       &drLagBasisPoints,
                                       vector<su2double>       &dsLagBasisPoints) {

}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesQuadrilateral(
                                       const unsigned short    nPoly,
                                       const vector<su2double> &rPoints,
                                       const vector<su2double> &sPoints,
                                       unsigned short          &nDOFs,
                                       vector<su2double>       &rDOFs,
                                       vector<su2double>       &sDOFs,
                                       vector<su2double>       &matVandermondeInv,
                                       vector<su2double>       &lagBasisPoints,
                                       vector<su2double>       &drLagBasisPoints,
                                       vector<su2double>       &dsLagBasisPoints) {

}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesTetrahedron(
                                       const unsigned short    nPoly,
                                       const vector<su2double> &rPoints,
                                       const vector<su2double> &sPoints,
                                       const vector<su2double> &tPoints,
                                       unsigned short          &nDOFs,
                                       vector<su2double>       &rDOFs,
                                       vector<su2double>       &sDOFs,
                                       vector<su2double>       &tDOFs,
                                       vector<su2double>       &matVandermondeInv,
                                       vector<su2double>       &lagBasisPoints,
                                       vector<su2double>       &drLagBasisPoints,
                                       vector<su2double>       &dsLagBasisPoints,
                                       vector<su2double>       &dtLagBasisPoints) {
}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesPyramid(
                                       const unsigned short    nPoly,
                                       const vector<su2double> &rPoints,
                                       const vector<su2double> &sPoints,
                                       const vector<su2double> &tPoints,
                                       unsigned short          &nDOFs,
                                       vector<su2double>       &rDOFs,
                                       vector<su2double>       &sDOFs,
                                       vector<su2double>       &tDOFs,
                                       vector<su2double>       &matVandermondeInv,
                                       vector<su2double>       &lagBasisPoints,
                                       vector<su2double>       &drLagBasisPoints,
                                       vector<su2double>       &dsLagBasisPoints,
                                       vector<su2double>       &dtLagBasisPoints) {
}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesPrism(
                                       const unsigned short    nPoly,
                                       const vector<su2double> &rPoints,
                                       const vector<su2double> &sPoints,
                                       const vector<su2double> &tPoints,
                                       unsigned short          &nDOFs,
                                       vector<su2double>       &rDOFs,
                                       vector<su2double>       &sDOFs,
                                       vector<su2double>       &tDOFs,
                                       vector<su2double>       &matVandermondeInv,
                                       vector<su2double>       &lagBasisPoints,
                                       vector<su2double>       &drLagBasisPoints,
                                       vector<su2double>       &dsLagBasisPoints,
                                       vector<su2double>       &dtLagBasisPoints) {
}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesHexahedron(
                                       const unsigned short    nPoly,
                                       const vector<su2double> &rPoints,
                                       const vector<su2double> &sPoints,
                                       const vector<su2double> &tPoints,
                                       unsigned short          &nDOFs,
                                       vector<su2double>       &rDOFs,
                                       vector<su2double>       &sDOFs,
                                       vector<su2double>       &tDOFs,
                                       vector<su2double>       &matVandermondeInv,
                                       vector<su2double>       &lagBasisPoints,
                                       vector<su2double>       &drLagBasisPoints,
                                       vector<su2double>       &dsLagBasisPoints,
                                       vector<su2double>       &dtLagBasisPoints) {
}

void CFEMStandardElementBase::MatMulRowMajor(const unsigned short nDOFs,
                                             const unsigned short nPoints,
                                             const vector<su2double> &A,
                                             const vector<su2double> &B,
                                             vector<su2double>       &C) {

}

void CFEMStandardElementBase::SubConnForPlottingLine(
                                         const unsigned short   nPoly,
                                         vector<unsigned short> &subConn) {

}

void CFEMStandardElementBase::SubConnForPlottingQuadrilateral(
                                         const unsigned short   nPoly,
                                         vector<unsigned short> &subConn) {

}

void CFEMStandardElementBase::SubConnForPlottingTriangle(
                                         const unsigned short   nPoly,
                                         vector<unsigned short> &subConn) {

}


void CFEMStandardElementBase::Vandermonde2D_Triangle(unsigned short          nPoly,
                                                     unsigned short          nDOFs,
                                                     const vector<su2double> &r,
                                                     const vector<su2double> &s,
                                                     vector<su2double>       &V) {

}

void CFEMStandardElementBase::GradVandermonde2D_Triangle(unsigned short          nPoly,
                                                         unsigned short          nDOFs,
                                                         const vector<su2double> &r,
                                                         const vector<su2double> &s,
                                                         vector<su2double>       &VDr,
                                                         vector<su2double>       &VDs) {

}

void CFEMStandardElementBase::Vandermonde2D_Quadrilateral(unsigned short          nPoly,
                                                          unsigned short          nDOFs,
                                                          const vector<su2double> &r,
                                                          const vector<su2double> &s,
                                                          vector<su2double>       &V) {

}

void CFEMStandardElementBase::GradVandermonde2D_Quadrilateral(unsigned short          nPoly,
                                                              unsigned short          nDOFs,
                                                              const vector<su2double> &r,
                                                              const vector<su2double> &s,
                                                              vector<su2double>       &VDr,
                                                              vector<su2double>       &VDs) {

}

void CFEMStandardElementBase::Vandermonde3D_Tetrahedron(unsigned short          nPoly,
                                                        unsigned short          nDOFs,
                                                        const vector<su2double> &r,
                                                        const vector<su2double> &s,
                                                        const vector<su2double> &t,
                                                        vector<su2double>       &V) {

}

void CFEMStandardElementBase::GradVandermonde3D_Tetrahedron(unsigned short          nPoly,
                                                            unsigned short          nDOFs,
                                                            const vector<su2double> &r,
                                                            const vector<su2double> &s,
                                                            const vector<su2double> &t,
                                                            vector<su2double>       &VDr,
                                                            vector<su2double>       &VDs,
                                                            vector<su2double>       &VDt) {

}

void CFEMStandardElementBase::Vandermonde3D_Pyramid(unsigned short          nPoly,
                                                    unsigned short          nDOFs,
                                                    const vector<su2double> &r,
                                                    const vector<su2double> &s,
                                                    const vector<su2double> &t,
                                                    vector<su2double>       &V) {

}

void CFEMStandardElementBase::GradVandermonde3D_Pyramid(unsigned short          nPoly,
                                                        unsigned short          nDOFs,
                                                        const vector<su2double> &r,
                                                        const vector<su2double> &s,
                                                        const vector<su2double> &t,
                                                        vector<su2double>       &VDr,
                                                        vector<su2double>       &VDs,
                                                        vector<su2double>       &VDt) {

}

void CFEMStandardElementBase::Vandermonde3D_Prism(unsigned short          nPoly,
                                                  unsigned short          nDOFs,
                                                  const vector<su2double> &r,
                                                  const vector<su2double> &s,
                                                  const vector<su2double> &t,
                                                  vector<su2double>       &V) {

}

void CFEMStandardElementBase::GradVandermonde3D_Prism(unsigned short          nPoly,
                                                      unsigned short          nDOFs,
                                                      const vector<su2double> &r,
                                                      const vector<su2double> &s,
                                                      const vector<su2double> &t,
                                                      vector<su2double>       &VDr,
                                                      vector<su2double>       &VDs,
                                                      vector<su2double>       &VDt) {

}

void CFEMStandardElementBase::Vandermonde3D_Hexahedron(unsigned short          nPoly,
                                                       unsigned short          nDOFs,
                                                       const vector<su2double> &r,
                                                       const vector<su2double> &s,
                                                       const vector<su2double> &t,
                                                       vector<su2double>       &V) {

}

void CFEMStandardElementBase::GradVandermonde3D_Hexahedron(unsigned short          nPoly,
                                                           unsigned short          nDOFs,
                                                           const vector<su2double> &r,
                                                           const vector<su2double> &s,
                                                           const vector<su2double> &t,
                                                           vector<su2double>       &VDr,
                                                           vector<su2double>       &VDs,
                                                           vector<su2double>       &VDt) {

}

su2double CFEMStandardElementBase::ViscousPenaltyParameter(
                                       const unsigned short VTK_TypeElem,
                                       const unsigned short nPolyElem) const {

  return 0.0;
}

/*----------------------------------------------------------------------------------*/
/*          Private member functions of CFEMStandardElementBase.                    */
/*----------------------------------------------------------------------------------*/

void CFEMStandardElementBase::GaussLegendrePoints1D(vector<su2double> &GLPoints,
                                                    vector<su2double> &GLWeights) {

}

su2double CFEMStandardElementBase::NormJacobi(unsigned short n,
                                              unsigned short alpha,
                                              unsigned short beta,
                                              su2double      x) {
  return 0.0;
}

su2double CFEMStandardElementBase::GradNormJacobi(unsigned short n,
                                                  unsigned short alpha,
                                                  unsigned short beta,
                                                  su2double      x) {

  return 0.0;
}

/*----------------------------------------------------------------------------------*/
/*           Public member functions of CFEMStandardElement.                        */
/*----------------------------------------------------------------------------------*/

CFEMStandardElement::CFEMStandardElement(unsigned short          val_VTK_Type,
                                         unsigned short          val_nPoly,
                                         bool                    val_constJac,
                                         CConfig                 *config,
                                         unsigned short          val_orderExact,
                                         const vector<su2double> *rLocSolDOFs,
                                         const vector<su2double> *sLocSolDOFs,
                                         const vector<su2double> *tLocSolDOFs)

  : CFEMStandardElementBase(val_VTK_Type, val_nPoly, val_constJac,
                            config, val_orderExact) {

}

void CFEMStandardElement::BasisFunctionsInPoint(const su2double   *parCoor,
                                                vector<su2double> &lagBasis) {

}

void CFEMStandardElement::BasisFunctionsAndDerivativesInPoint(
                                           const su2double            *parCoor,
                                           vector<su2double>          &lagBasis,
                                           vector<vector<su2double> > &dLagBasis) {

}

bool CFEMStandardElement::SameStandardElement(unsigned short val_VTK_Type,
                                              unsigned short val_nPoly,
                                              bool           val_constJac) {
  return false;
}

/*----------------------------------------------------------------------------------*/
/*           Private member functions of CFEMStandardElement.                       */
/*----------------------------------------------------------------------------------*/

void CFEMStandardElement::Copy(const CFEMStandardElement &other) {

}

void CFEMStandardElement::CreateBasisFunctionsAndMatrixDerivatives(
                                       const vector<su2double> &rLoc,
                                       const vector<su2double> &sLoc,
                                       const vector<su2double> &tLoc,
                                             vector<su2double> &matVandermondeInv,
                                             vector<su2double> &lagBasis,
                                             vector<su2double> &matDerBasis) {

}

void CFEMStandardElement::DataStandardLine(void) {

}

void CFEMStandardElement::DataStandardTriangle(void) {

}

void CFEMStandardElement::DataStandardQuadrilateral(void) {

}

void CFEMStandardElement::DataStandardTetrahedron(void) {

}

void CFEMStandardElement::DataStandardPyramid(void) {

}

void CFEMStandardElement::DataStandardPrism(void) {

}

void CFEMStandardElement::DataStandardHexahedron(void) {

}

void CFEMStandardElement::SubConnTetrahedron(void) {

}

void CFEMStandardElement::SubConnPyramid(void) {

}

void CFEMStandardElement::SubConnPrism(void) {

}

void CFEMStandardElement::SubConnHexahedron(void) {

}

unsigned short CFEMStandardElement::GetNDOFsPerSubElem(unsigned short val_VTK_Type) const {

  return 0;
}

void CFEMStandardElement::ChangeDirectionQuadConn(vector<unsigned short> &connQuad,
                                                  unsigned short         vert0,
                                                  unsigned short         vert1,
                                                  unsigned short         vert2,
                                                  unsigned short         vert3) const {

}

void CFEMStandardElement::ChangeDirectionTriangleConn(vector<unsigned short> &connTriangle,
                                                      unsigned short         vert0,
                                                      unsigned short         vert1,
                                                      unsigned short         vert2) const {

}

/*----------------------------------------------------------------------------------*/
/*         Public member functions of CFEMStandardInternalFace.                     */
/*----------------------------------------------------------------------------------*/

CFEMStandardInternalFace::CFEMStandardInternalFace(unsigned short val_VTK_TypeFace,
                                                   unsigned short val_VTK_TypeSide0,
                                                   unsigned short val_nPolySide0,
                                                   unsigned short val_VTK_TypeSide1,
                                                   unsigned short val_nPolySide1,
                                                   bool           val_constJac,
                                                   bool           val_swapFaceInElementSide0,
                                                   bool           val_swapFaceInElementSide1,
                                                   CConfig        *config,
                                                   unsigned short val_orderExact)

  : CFEMStandardElementBase(val_VTK_TypeFace, max(val_nPolySide0, val_nPolySide1),
                            val_constJac, config, val_orderExact) {

}

bool CFEMStandardInternalFace::SameStandardMatchingFace(unsigned short val_VTK_TypeFace,
                                                        bool           val_constJac,
                                                        unsigned short val_VTK_TypeSide0,
                                                        unsigned short val_nPolySide0,
                                                        unsigned short val_VTK_TypeSide1,
                                                        unsigned short val_nPolySide1,
                                                        bool           val_swapFaceInElementSide0,
                                                        bool           val_swapFaceInElementSide1) {
  return false;
}

/*----------------------------------------------------------------------------------*/
/*        Private member functions of CFEMStandardInternalFace.                     */
/*----------------------------------------------------------------------------------*/

void CFEMStandardInternalFace::Copy(const CFEMStandardInternalFace &other) {

}

/*----------------------------------------------------------------------------------*/
/*         Public member functions of CFEMStandardBoundaryFace.                     */
/*----------------------------------------------------------------------------------*/

CFEMStandardBoundaryFace::CFEMStandardBoundaryFace(unsigned short val_VTK_TypeFace,
                                                   unsigned short val_VTK_TypeElem,
                                                   unsigned short val_nPolyElem,
                                                   bool           val_constJac,
                                                   bool           val_swapFaceInElement,
                                                   CConfig        *config,
                                                   unsigned short val_orderExact)

  : CFEMStandardElementBase(val_VTK_TypeFace, val_nPolyElem,
                            val_constJac, config, val_orderExact) {

}

bool CFEMStandardBoundaryFace::SameStandardBoundaryFace(unsigned short val_VTK_TypeFace,
                                                        bool           val_constJac,
                                                        unsigned short val_VTK_TypeElem,
                                                        unsigned short val_nPolyElem,
                                                        bool           val_swapFaceInElem) {
  return false;
}

/*----------------------------------------------------------------------------------*/
/*        Private member functions of CFEMStandardBoundaryFace.                     */
/*----------------------------------------------------------------------------------*/

void CFEMStandardBoundaryFace::Copy(const CFEMStandardBoundaryFace &other) {

}
