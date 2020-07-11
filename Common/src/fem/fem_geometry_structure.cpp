/*!
 * \file fem_geometry_structure.cpp
 * \brief Functions for creating the primal grid for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
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

#include "../../include/fem/fem_geometry_structure.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridFEM.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridBoundFEM.hpp"
#include "../../include/adt_structure.hpp"

CMeshFEM::CMeshFEM(CGeometry *geometry, CConfig *config) {

}

void CMeshFEM::ComputeGradientsCoordinatesFace(const unsigned short nIntegration,
                                               const unsigned short nDOFs,
                                               const su2double      *matDerBasisInt,
                                               const unsigned long  *DOFs,
                                               su2double            *derivCoor,
                                               CConfig              *config) {

}

void CMeshFEM::ComputeGradientsCoorWRTParam(const unsigned short nIntegration,
                                            const unsigned short nDOFs,
                                            const su2double      *matDerBasisInt,
                                            const unsigned long  *DOFs,
                                            su2double            *derivCoor,
                                            CConfig              *config) {

}

void CMeshFEM::MetricTermsBoundaryFaces(CBoundaryFEM *boundary,
                                        CConfig      *config) {

}

void CMeshFEM::SetPositive_ZArea(CConfig *config) {

}

CMeshFEM_DG::CMeshFEM_DG(CGeometry *geometry, CConfig *config)
  : CMeshFEM(geometry, config) {
}

void CMeshFEM_DG::SetGlobal_to_Local_Point(void) {
}

void CMeshFEM_DG::CoordinatesIntegrationPoints(void) {

}

void CMeshFEM_DG::CoordinatesSolDOFs(void) {

}

void CMeshFEM_DG::CreateFaces(CConfig *config) {

}

void CMeshFEM_DG::CreateStandardVolumeElements(CConfig *config) {

}

void CMeshFEM_DG::SetSendReceive(CConfig *config) {

}

void CMeshFEM_DG::CreateConnectivitiesFace(
                                const unsigned short        VTK_TypeFace,
                                const unsigned long         *cornerPointsFace,
                                const unsigned short        VTK_TypeElem,
                                const unsigned short        nPolyGrid,
                                const vector<unsigned long> &elemNodeIDsGrid,
                                const unsigned short        nPolyConn,
                                const unsigned long         *connElem,
                                bool                        &swapFaceInElement,
                                unsigned long               *modConnFace,
                                unsigned long               *modConnElem) {

}

void CMeshFEM_DG::CreateConnectivitiesLineAdjacentQuadrilateral(
                                const unsigned long         *cornerPointsLine,
                                const unsigned short        nPolyGrid,
                                const vector<unsigned long> &quadNodeIDsGrid,
                                const unsigned short        nPolyConn,
                                const unsigned long         *connQuad,
                                unsigned long               *modConnLine,
                                unsigned long               *modConnQuad) {

}

void CMeshFEM_DG::CreateConnectivitiesLineAdjacentTriangle(
                                const unsigned long         *cornerPointsLine,
                                const unsigned short        nPolyGrid,
                                const vector<unsigned long> &triaNodeIDsGrid,
                                const unsigned short        nPolyConn,
                                const unsigned long         *connTria,
                                unsigned long               *modConnLine,
                                unsigned long               *modConnTria) {

}

void CMeshFEM_DG::CreateConnectivitiesQuadrilateralAdjacentHexahedron(
                                const unsigned long         *cornerPointsQuad,
                                const unsigned short        nPolyGrid,
                                const vector<unsigned long> &hexaNodeIDsGrid,
                                const unsigned short        nPolyConn,
                                const unsigned long         *connHexa,
                                unsigned long               *modConnQuad,
                                unsigned long               *modConnHexa) {

}

void CMeshFEM_DG::CreateConnectivitiesQuadrilateralAdjacentPrism(
                                const unsigned long         *cornerPointsQuad,
                                const unsigned short        nPolyGrid,
                                const vector<unsigned long> &prismNodeIDsGrid,
                                const unsigned short        nPolyConn,
                                const unsigned long         *connPrism,
                                bool                        &swapFaceInElement,
                                unsigned long               *modConnQuad,
                                unsigned long               *modConnPrism) {

}

void CMeshFEM_DG::CreateConnectivitiesQuadrilateralAdjacentPyramid(
                                const unsigned long         *cornerPointsQuad,
                                const unsigned short        nPolyGrid,
                                const vector<unsigned long> &pyraNodeIDsGrid,
                                const unsigned short        nPolyConn,
                                const unsigned long         *connPyra,
                                unsigned long               *modConnQuad,
                                unsigned long               *modConnPyra) {

}

void CMeshFEM_DG::CreateConnectivitiesTriangleAdjacentPrism(
                                const unsigned long         *cornerPointsTria,
                                const unsigned short        nPolyGrid,
                                const vector<unsigned long> &prismNodeIDsGrid,
                                const unsigned short        nPolyConn,
                                const unsigned long         *connPrism,
                                unsigned long               *modConnTria,
                                unsigned long               *modConnPrism) {

}

void CMeshFEM_DG::CreateConnectivitiesTriangleAdjacentPyramid(
                                const unsigned long         *cornerPointsTria,
                                const unsigned short        nPolyGrid,
                                const vector<unsigned long> &pyraNodeIDsGrid,
                                const unsigned short        nPolyConn,
                                const unsigned long         *connPyra,
                                bool                        &swapFaceInElement,
                                unsigned long               *modConnTria,
                                unsigned long               *modConnPyra) {

}

void CMeshFEM_DG::CreateConnectivitiesTriangleAdjacentTetrahedron(
                                const unsigned long         *cornerPointsTria,
                                const unsigned short        nPolyGrid,
                                const vector<unsigned long> &tetNodeIDsGrid,
                                const unsigned short        nPolyConn,
                                const unsigned long         *connTet,
                                unsigned long               *modConnTria,
                                unsigned long               *modConnTet) {

}

void CMeshFEM_DG::MetricTermsMatchingFaces(CConfig *config) {

}

void CMeshFEM_DG::LengthScaleVolumeElements(void) {

}

void CMeshFEM_DG::MetricTermsSurfaceElements(CConfig *config) {

}

void CMeshFEM_DG::MetricTermsVolumeElements(CConfig *config) {

}

void CMeshFEM_DG::TimeCoefficientsPredictorADER_DG(CConfig *config) {

}

void CMeshFEM_DG::VolumeMetricTermsFromCoorGradients(
                                       const unsigned short nEntities,
                                       const su2double      *gradCoor,
                                       vector<su2double>    &metricTerms) {

}

void CMeshFEM_DG::WallFunctionPreprocessing(CConfig *config) {

}

void CMeshFEM_DG::HighOrderContainmentSearch(const su2double      *coor,
                                             const unsigned long  parElem,
                                             const unsigned short subElem,
                                             const su2double      *weightsSubElem,
                                             su2double            *parCoor) {

}

void CMeshFEM_DG::InitStaticMeshMovement(CConfig              *config,
                                         const unsigned short Kind_Grid_Movement,
                                         const unsigned short iZone) {

}

CDummyMeshFEM_DG::CDummyMeshFEM_DG(CConfig *config): CMeshFEM_DG() {

}
