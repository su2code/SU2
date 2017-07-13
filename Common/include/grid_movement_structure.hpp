/*!
 * \file grid_movement_structure.hpp
 * \brief Headers of the main subroutines for doing the numerical grid 
 *        movement (including volumetric movement, surface movement and Free From 
 *        technique definition). The subroutines and functions are in 
 *        the <i>grid_movement_structure.cpp</i> file.
 * \author F. Palacios, T. Economon, S. Padron
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

#pragma once

#include "surface_deformation.hpp"
#include "matrix_structure.hpp"
#include "vector_structure.hpp"
#include "linear_solvers_structure.hpp"

using namespace std;

/*!
 * \class CGridMovement
 * \brief Class for moving the surface and volumetric 
 *        numerical grid (2D and 3D problems).
 * \author F. Palacios
 * \version 5.0.0 "Raven"
 */
class CGridMovement {
public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CGridMovement(void);

	/*! 
	 * \brief Destructor of the class. 
	 */
         virtual ~CGridMovement(void);
};


/*! 
 * \class CVolumetricMovement
 * \brief Class for moving the volumetric numerical grid.
 * \author F. Palacios, A. Bueno, T. Economon, S. Padron.
 * \version 5.0.0 "Raven"
 */
class CVolumetricMovement : public CGridMovement {
protected:

	unsigned short nDim;		/*!< \brief Number of dimensions. */
	unsigned short nVar;		/*!< \brief Number of variables. */
  
	unsigned long nPoint;		/*!< \brief Number of points. */
	unsigned long nPointDomain;		/*!< \brief Number of points in the domain. */

	unsigned long nIterMesh;	/*!< \brief Number of iterations in the mesh update. +*/

  CSysMatrix StiffMatrix; /*!< \brief Matrix to store the point-to-point stiffness. */
  CSysVector LinSysSol;
  CSysVector LinSysRes;

public:

	/*! 
	 * \brief Constructor of the class.
	 */
	CVolumetricMovement(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CVolumetricMovement(void);
  
	/*!
	 * \brief Update the value of the coordinates after the grid movement.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void UpdateGridCoord(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Update the dual grid after the grid movement (edges and control volumes).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void UpdateDualGrid(CGeometry *geometry, CConfig *config);
  
	/*! 
	 * \brief Update the coarse multigrid levels after the grid movement.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void UpdateMultiGrid(CGeometry **geometry, CConfig *config);
  
  /*!
	 * \brief Compute the stiffness matrix for grid deformation using spring analogy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \return Value of the length of the smallest edge of the grid.
	 */
	su2double SetFEAMethodContributions_Elem(CGeometry *geometry, CConfig *config);
  
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
  void SetFEA_StiffMatrix3D(CGeometry *geometry, CConfig *config, su2double **StiffMatrix_Elem, unsigned long PointCorners[8], su2double CoordCorners[8][3],
                            unsigned short nNodes, su2double ElemVolume, su2double ElemDistance);
  
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
  void SetFEA_StiffMatrix2D(CGeometry *geometry, CConfig *config, su2double **StiffMatrix_Elem, unsigned long PointCorners[8], su2double CoordCorners[8][3],
                            unsigned short nNodes, su2double ElemVolume, su2double ElemDistance);
    
  /*!
	 * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Zeta - Local coordinates.
	 * \param[in] CoordCorners - Coordiantes of the corners.
   * \param[in] DShapeFunction - Shape function information
	 */
  su2double ShapeFunc_Hexa(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]);
  
  /*!
	 * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Zeta - Local coordinates.
	 * \param[in] CoordCorners - Coordiantes of the corners.
   * \param[in] DShapeFunction - Shape function information
	 */
  su2double ShapeFunc_Tetra(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]);
  
  /*!
	 * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Zeta - Local coordinates.
	 * \param[in] CoordCorners - Coordiantes of the corners.
   * \param[in] DShapeFunction - Shape function information
	 */
  su2double ShapeFunc_Pyram(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]);
  
  /*!
	 * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Zeta - Local coordinates.
	 * \param[in] CoordCorners - Coordiantes of the corners.
   * \param[in] DShapeFunction - Shape function information
	 */
  su2double ShapeFunc_Prism(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]);
  
  /*!
	 * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
	 * \param[in] CoordCorners - Coordiantes of the corners.
   * \param[in] DShapeFunction - Shape function information
	 */
  su2double ShapeFunc_Triangle(su2double Xi, su2double Eta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]);
  
  /*!
	 * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
	 * \param[in] CoordCorners - Coordiantes of the corners.
   * \param[in] DShapeFunction - Shape function information
	 */
  su2double ShapeFunc_Quadrilateral(su2double Xi, su2double Eta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]);
  
  /*!
	 * \brief Compute the shape functions for hexahedron
	 * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
	 */
  su2double GetHexa_Volume(su2double CoordCorners[8][3]);
  
  /*!
	 * \brief Compute the shape functions for hexahedron
	 * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
	 */
  su2double GetTetra_Volume(su2double CoordCorners[8][3]);
  
  /*!
	 * \brief Compute the shape functions for hexahedron
	 * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
	 */
  su2double GetPrism_Volume(su2double CoordCorners[8][3]);
  
  /*!
	 * \brief Compute the shape functions for hexahedron
	 * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
	 */
  su2double GetPyram_Volume(su2double CoordCorners[8][3]);
  
  /*!
	 * \brief Compute the shape functions for hexahedron
	 * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
	 */
  su2double GetTriangle_Area(su2double CoordCorners[8][3]);
  
  /*!
	 * \brief Compute the shape functions for hexahedron
	 * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
	 */
  su2double GetQuadrilateral_Area(su2double CoordCorners[8][3]);
    
  /*!
	 * \brief Add the stiffness matrix for a 2-D triangular element to the global stiffness matrix for the entire mesh (node-based).
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] StiffMatrix_Elem - Element stiffness matrix to be filled.
   * \param[in] PointCorners - Index values for element corners
   * \param[in] nNodes - Number of nodes defining the element.
	 */
  void AddFEA_StiffMatrix(CGeometry *geometry, su2double **StiffMatrix_Elem, unsigned long PointCorners[8], unsigned short nNodes);
  
  /*!
	 * \brief Check for negative volumes (all elements) after performing grid deformation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
  void ComputeDeforming_Element_Volume(CGeometry *geometry, su2double &MinVolume, su2double &MaxVolume);
  
  
  /*!
	 * \brief Compute the minimum distance to the nearest deforming surface.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
	 */
  void ComputeDeforming_Wall_Distance(CGeometry *geometry, CConfig *config, su2double &MinDistance, su2double &MaxDistance);
    
	/*!
	 * \brief Check the boundary vertex that are going to be moved.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetBoundaryDisplacements(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Check the domain points vertex that are going to be moved.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetDomainDisplacements(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Unsteady grid movement using rigid mesh rotation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
	 */
	void Rigid_Rotation(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter);
	
  /*!
	 * \brief Unsteady pitching grid movement using rigid mesh motion.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
	 */
	void Rigid_Pitching(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter);
  
  /*!
	 * \brief Unsteady plunging grid movement using rigid mesh motion.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
	 */
	void Rigid_Plunging(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter);
  
  /*!
	 * \brief Unsteady translational grid movement using rigid mesh motion.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Zone number in the mesh.
   * \param[in] iter - Physical time iteration number.
	 */
	void Rigid_Translation(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter);
  
  /*!
   * \brief Scale the volume grid by a multiplicative factor.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] UpdateGeo - Update geometry.
   */
  void SetVolume_Scaling(CGeometry *geometry, CConfig *config, bool UpdateGeo);
  
  /*!
   * \brief Translate the volume grid by a specified displacement vector.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] UpdateGeo - Update geometry.
   */
  void SetVolume_Translation(CGeometry *geometry, CConfig *config, bool UpdateGeo);
  
  /*!
   * \brief Rotate the volume grid around a specified axis and angle.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] UpdateGeo - Update geometry.
   */
  void SetVolume_Rotation(CGeometry *geometry, CConfig *config, bool UpdateGeo);
  
  /*!
	 * \brief Grid deformation using the spring analogy method.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] UpdateGeo - Update geometry.
   * \param[in] Derivative - Compute the derivative (disabled by default). Does not actually deform the grid if enabled.
	 */
  void SetVolume_Deformation(CGeometry *geometry, CConfig *config, bool UpdateGeo, bool Derivative = false);

  /*!
   * \brief Set the derivatives of the boundary nodes.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetBoundaryDerivatives(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Update the derivatives of the coordinates after the grid movement.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void UpdateGridCoord_Derivatives(CGeometry *geometry, CConfig *config);

	/*!
	 * \brief Compute the determinant of a 3 by 3 matrix.
	 * 3 by 3 matrix elements
	 * \param[in] A00
	 * \param[in] A01
	 * \param[in] A02
	 * \param[in] A10
	 * \param[in] A11
	 * \param[in] A12
	 * \param[in] A20
	 * \param[in] A21
	 * \param[in] A22
	 * \result Determinant of the matrix
	 */
	su2double Determinant_3x3(su2double A00, su2double A01, su2double A02, su2double A10, su2double A11, su2double A12, su2double A20, su2double A21, su2double A22);


	/*!
	 * \brief Store the number of iterations when moving the mesh.
	 * \param[in] val_nIterMesh - Number of iterations.
	 */
	void Set_nIterMesh(unsigned long val_nIterMesh);

	/*!
	 * \brief Retrieve the number of iterations when moving the mesh.
	 * \param[out] Number of iterations.
	 */
	unsigned long Get_nIterMesh(void);
};


#include "grid_movement_structure.inl"
