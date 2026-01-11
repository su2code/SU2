/*!
 * \file CLinearElasticity.hpp
 * \brief Headers of the CLinearElasticity class.
 * \author F. Palacios, A. Bueno, T. Economon, S. Padron.
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
#include "CVolumetricMovement.hpp"
#include "../linear_algebra/CSysMatrix.hpp"
#include "../linear_algebra/CSysVector.hpp"
#include "../linear_algebra/CSysSolve.hpp"

/*!
 * \class CLinearElasticity
 * \brief Class for moving the volumetric numerical grid using the linear elasticity analogy.
 * \author F. Palacios, A. Bueno, T. Economon, S. Padron.
 */
class CLinearElasticity final : public CVolumetricMovement {
 protected:
  unsigned short nVar; /*!< \brief Number of variables. */

  unsigned long nPoint;       /*!< \brief Number of points. */
  unsigned long nPointDomain; /*!< \brief Number of points in the domain. */

  unsigned long nIterMesh; /*!< \brief Number of iterations in the mesh update. +*/

#ifndef CODI_FORWARD_TYPE
  CSysMatrix<su2mixedfloat> StiffMatrix; /*!< \brief Stiffness matrix of the elasticity problem. */
  CSysSolve<su2mixedfloat> System;       /*!< \brief Linear solver/smoother. */
#else
  CSysMatrix<su2double> StiffMatrix;
  CSysSolve<su2double> System;
#endif
  CSysVector<su2double> LinSysSol;
  CSysVector<su2double> LinSysRes;

 public:
  /*!
   * \brief Constructor of the class.
   */
  CLinearElasticity(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CLinearElasticity() override;

  /*!
   * \brief Grid deformation using the spring analogy method.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] UpdateGeo - Update geometry.
   * \param[in] Derivative - Compute the derivative (disabled by default). Does not actually deform the grid if enabled.
   */
  void SetVolume_Deformation(CGeometry* geometry, CConfig* config, bool UpdateGeo, bool Derivative,
                             bool ForwardProjectionDerivative) override;

 private:
  /*!
   * \brief Update the value of the coordinates after the grid movement.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void UpdateGridCoord(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Update the derivatives of the coordinates after the grid movement.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void UpdateGridCoord_Derivatives(CGeometry* geometry, CConfig* config, bool ForwardProjectionDerivative);

  /*!
   * \brief Compute the minimum distance to the nearest solid surface.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeSolid_Wall_Distance(CGeometry* geometry, CConfig* config, su2double& MinDistance,
                                  su2double& MaxDistance) const;

  /*!
   * \brief Compute the stiffness matrix for grid deformation using spring analogy.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \return Value of the length of the smallest edge of the grid.
   */
  su2double SetFEAMethodContributions_Elem(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Build the stiffness matrix for a 3-D hexahedron element. The result will be placed in StiffMatrix_Elem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] StiffMatrix_Elem - Element stiffness matrix to be filled.
   * \param[in] CoordCorners - Index value for Node 1 of the current hexahedron.
   * \param[in] PointCorners - Index values for element corners
   * \param[in] nNodes - Number of nodes defining the element.
   * \param[in] scale
   */
  void SetFEA_StiffMatrix2D(CGeometry* geometry, CConfig* config, su2double** StiffMatrix_Elem,
                            unsigned long PointCorners[8], su2double CoordCorners[8][3], unsigned short nNodes,
                            su2double ElemVolume, su2double ElemDistance);

  /*!
   * \brief Build the stiffness matrix for a 3-D hexahedron element. The result will be placed in StiffMatrix_Elem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] StiffMatrix_Elem - Element stiffness matrix to be filled.
   * \param[in] CoordCorners - Index value for Node 1 of the current hexahedron.
   * \param[in] PointCorners - Index values for element corners
   * \param[in] nNodes - Number of nodes defining the element.
   * \param[in] scale
   */
  void SetFEA_StiffMatrix3D(CGeometry* geometry, CConfig* config, su2double** StiffMatrix_Elem,
                            unsigned long PointCorners[8], su2double CoordCorners[8][3], unsigned short nNodes,
                            su2double ElemVolume, su2double ElemDistance);

  /*!
   * \brief Add the stiffness matrix for a 2-D triangular element to the global stiffness matrix for the entire mesh
   * (node-based). \param[in] geometry - Geometrical definition of the problem. \param[in] StiffMatrix_Elem - Element
   * stiffness matrix to be filled. \param[in] PointCorners - Index values for element corners \param[in] nNodes -
   * Number of nodes defining the element.
   */
  void AddFEA_StiffMatrix(CGeometry* geometry, su2double** StiffMatrix_Elem, unsigned long PointCorners[8],
                          unsigned short nNodes);

  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Zeta - Local coordinates.
   * \param[in] CoordCorners - Coordiantes of the corners.
   * \param[in] DShapeFunction - Shape function information
   */
  su2double ShapeFunc_Hexa(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3],
                           su2double DShapeFunction[8][4]);

  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Zeta - Local coordinates.
   * \param[in] CoordCorners - Coordiantes of the corners.
   * \param[in] DShapeFunction - Shape function information
   */
  su2double ShapeFunc_Tetra(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3],
                            su2double DShapeFunction[8][4]);

  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Zeta - Local coordinates.
   * \param[in] CoordCorners - Coordiantes of the corners.
   * \param[in] DShapeFunction - Shape function information
   */
  su2double ShapeFunc_Pyram(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3],
                            su2double DShapeFunction[8][4]);

  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Zeta - Local coordinates.
   * \param[in] CoordCorners - Coordiantes of the corners.
   * \param[in] DShapeFunction - Shape function information
   */
  su2double ShapeFunc_Prism(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3],
                            su2double DShapeFunction[8][4]);

  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] CoordCorners - Coordiantes of the corners.
   * \param[in] DShapeFunction - Shape function information
   */
  su2double ShapeFunc_Triangle(su2double Xi, su2double Eta, su2double CoordCorners[8][3],
                               su2double DShapeFunction[8][4]);

  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] CoordCorners - Coordiantes of the corners.
   * \param[in] DShapeFunction - Shape function information
   */
  su2double ShapeFunc_Quadrilateral(su2double Xi, su2double Eta, su2double CoordCorners[8][3],
                                    su2double DShapeFunction[8][4]);

  /*!
   * \brief Check the domain points vertex that are going to be moved.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetDomainDisplacements(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Check the boundary vertex that are going to be moved.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetBoundaryDisplacements(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Set the derivatives of the boundary nodes.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetBoundaryDerivatives(CGeometry* geometry, CConfig* config, bool ForwardProjectionDerivative);
};
