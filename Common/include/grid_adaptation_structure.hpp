/*!
 * \file grid_adaptation_structure.hpp
 * \brief Headers of the main subroutines for doing the numerical grid
 *        adaptation.
 * \author F. Palacios
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "./mpi_structure.hpp"

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>

#include "geometry_structure.hpp"
#include "config_structure.hpp"

using namespace std;

/*!
 * \class CGridAdaptation
 * \brief Parent class for defining the grid adaptation.
 * \author F. Palacios
 * \version 2.0.7
 */
class CGridAdaptation {
protected:
  int rank, 	/*!< \brief MPI Rank. */
  size;       	/*!< \brief MPI Size. */
	unsigned long nPoint_new,	/*!< \brief Number of new points. */
	nElem_new;					/*!< \brief Number of new elements. */
	unsigned short nDim,	/*!< \brief Number of dimensions of the problem. */
	nVar;					/*!< \brief Number of variables in the problem. */
	su2double **ConsVar_Sol,	/*!< \brief Conservative variables (original solution). */
	**ConsVar_Res,			/*!< \brief Conservative variables (residual). */
	**ConsVar_Adapt;		/*!< \brief Conservative variables (adapted solution). */
	su2double **AdjVar_Sol,	/*!< \brief Adjoint variables (original solution). */
	**AdjVar_Res,			/*!< \brief Adjoint variables (residual). */
	**AdjVar_Adapt;			/*!< \brief Adjoint variables (adapted solution). */
	su2double **LinVar_Sol,	/*!< \brief Linear variables (original solution). */
	**LinVar_Res,			/*!< \brief Linear variables (residual). */
	**LinVar_Adapt;			/*!< \brief Linear variables (adapted solution). */
	su2double **Gradient,		/*!< \brief Gradient value. */
	**Gradient_Flow,		/*!< \brief Gradient of the flow variables. */
	**Gradient_Adj;			/*!< \brief Fradient of the adjoint variables. */
	su2double *Index;			/*!< \brief Adaptation index (indicates the value of the adaptation). */
	
public:

	/*! 
	 * \brief Constructor of the class.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CGridAdaptation(CGeometry *geometry, CConfig *config);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CGridAdaptation(void);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void GetFlowSolution(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void GetFlowResidual(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void GetAdjSolution(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void GetAdjResidual(CGeometry *geometry, CConfig *config);
	
	/*!
	 * \brief Do a complete adaptation of the computational grid.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] strength - Adaptation Strength.	 
	 */		
	void SetComplete_Refinement(CGeometry *geometry, unsigned short strength);
	
	/*! 
	 * \brief Do not do any kind of adaptation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] strength - Adaptation Strength.	 
	 */		
	void SetNo_Refinement(CGeometry *geometry, unsigned short strength);
	
	/*! 
	 * \brief Do an adaptation of the computational grid on the wake.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] strength - Adaptation Strength.	 
	 */		
	void SetWake_Refinement(CGeometry *geometry, unsigned short strength);
	
	/*! 
	 * \brief Do an adaptation of the computational grid on the supersonic shock region.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */		
	void SetSupShock_Refinement(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Do an adaptation of the computational grid on a near field boundary.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */		
	void SetNearField_Refinement(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Do a complete adaptation of the computational grid using a homothetic technique (2D).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] geo_adapt - Geometrical definition of the adapted grid.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void SetHomothetic_Adaptation2D(CGeometry *geometry, CPhysicalGeometry *geo_adapt, CConfig *config);
	
	/*! 
	 * \brief Do a complete adaptation of the computational grid using a homothetic technique (3D).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] geo_adapt - Geometrical definition of the adapted grid.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void SetHomothetic_Adaptation3D(CGeometry *geometry, CPhysicalGeometry *geo_adapt, CConfig *config);
		
	/*! 
	 * \brief Find the adaptation code for each element in the fine grid.
	 * \param[in] AdaptCode - Edge combination to stablish the right elemeent division.
	 * \return Adaptation code for the element.
	 */	
	long CheckTriangleCode(bool *AdaptCode);
	
	/*! 
	 * \brief Find the adaptation code for each element in the fine grid.
	 * \param[in] AdaptCode - Edge combination to stablish the right elemeent division.
	 * \return Adaptation code for the element.
	 */	
	long CheckRectCode(bool *AdaptCode);
	
	/*! 
	 * \brief Find the adaptation code for each element in the fine grid.
	 * \param[in] AdaptCode - Edge combination to stablish the right elemeent division.
	 * \return Adaptation code for the element.
	 */	
	long CheckRectExtCode(bool *AdaptCode);
	
	/*! 
	 * \brief Find the adaptation code for each element in the fine grid.
	 * \param[in] AdaptCode - Edge combination to stablish the right elemeent division.
	 * \return Adaptation code for the element.
	 */	
	long CheckTetraCode(bool *AdaptCode);
	
	/*! 
	 * \brief Find the adaptation code for each element in the fine grid.
	 * \param[in] AdaptCode - Edge combination to stablish the right elemeent division.
	 * \return Adaptation code for the element.
	 */	
	long CheckHexaCode(bool *AdaptCode);
	
	/*! 
	 * \brief Find the adaptation code for each element in the fine grid.
	 * \param[in] AdaptCode - Edge combination to stablish the right elemeent division.
	 * \return Adaptation code for the element.
	 */	
	long CheckPyramCode(bool *AdaptCode);
	
	/*! 
	 * \brief Division pattern of the element.
	 * \param[in] code - number that identify the division.
	 * \param[in] nodes - Nodes that compose the element, including new nodes.
	 * \param[in] edges - Edges that compose the element.
	 * \param[out] Division - Division pattern. 
	 * \param[out] nPart - Number of new elements after the division. 
	 */	
	void TriangleDivision(long code, long *nodes, long *edges, long **Division, long *nPart);
	
	/*! 
	 * \brief Division pattern of the element.
	 * \param[in] code - number that identify the division.
	 * \param[in] nodes - Nodes that compose the element, including new nodes.
	 * \param[in] edges - Edges that compose the element.
	 * \param[out] Division - Division pattern. 
	 * \param[out] nPart - Number of new elements after the division. 
	 */	
	void RectDivision(long code, long *nodes, long **Division, long *nPart);
	
	/*! 
	 * \brief Division pattern of the element.
	 * \param[in] code - number that identify the division.
	 * \param[in] nodes - Nodes that compose the element, including new nodes.
	 * \param[in] edges - Edges that compose the element.
	 * \param[out] Division - Division pattern. 
	 * \param[out] nPart - Number of new elements after the division. 
	 */	
	void RectExtDivision(long code, long *nodes, long **Division, long *nPart);
	
	/*! 
	 * \brief Division pattern of the element.
	 * \param[in] code - number that identify the division.
	 * \param[in] nodes - Nodes that compose the element, including new nodes.
	 * \param[in] edges - Edges that compose the element.
	 * \param[out] Division - Division pattern. 
	 * \param[out] nPart - Number of new elements after the division. 
	 */	
	void TetraDivision(long code, long *nodes, long *edges, long **Division, long *nPart);
	
	/*! 
	 * \brief Division pattern of the element.
	 * \param[in] code - number that identify the division.
	 * \param[in] nodes - Nodes that compose the element, including new nodes.
	 * \param[in] edges - Edges that compose the element.
	 * \param[out] Division - Division pattern. 
	 * \param[out] nPart - Number of new elements after the division. 
	 */	
	void HexaDivision(long code, long *nodes, long **Division, long *nPart);
	
	/*! 
	 * \brief Division pattern of the element.
	 * \param[in] code - number that identify the division.
	 * \param[in] nodes - Nodes that compose the element, including new nodes.
	 * \param[in] edges - Edges that compose the element.
	 * \param[out] Division - Division pattern. 
	 * \param[out] nPart - Number of new elements after the division. 
	 */	
	void PyramDivision(long code, long *nodes, long **Division, long *nPart);
	
	/*! 
	 * \brief Do a complete adaptation of the computational grid.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] strength _________________________.
	 */		
	void SetIndicator_Flow(CGeometry *geometry, CConfig *config, unsigned short strength);
	
	/*! 
	 * \brief Do a complete adaptation of the computational grid.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] strength _________________________.
	 */		
	void SetIndicator_Adj(CGeometry *geometry, CConfig *config, unsigned short strength);
	
	/*! 
	 * \brief Do a complete adaptation of the computational grid.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */		
	void SetIndicator_FlowAdj(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void SetIndicator_Robust(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void SetIndicator_Computable(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void SetIndicator_Computable_Robust(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Write the restart file with the adapted grid.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] mesh_flowfilename - _________________________.
	 */		
	void SetRestart_FlowSolution(CConfig *config, CPhysicalGeometry *geo_adapt, string mesh_flowfilename);
	
	/*! 
	 * \brief Write the restart file with the adapted grid.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] mesh_adjfilename - _________________________.
	 */		
	void SetRestart_AdjSolution(CConfig *config, CPhysicalGeometry *geo_adapt, string mesh_adjfilename);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] mesh_linfilename - _________________________.
	 */	
	void SetRestart_LinSolution(CConfig *config, CPhysicalGeometry *geo_adapt, string mesh_linfilename);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] max_elem - _________________________.
	 */	
	void SetSensorElem(CGeometry *geometry, CConfig *config, unsigned long max_elem);

};

#include "grid_adaptation_structure.inl"


