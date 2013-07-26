/*!
 * \file grid_adaptation_structure.hpp
 * \brief Header file for the adaptation subroutines.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
 */

#pragma once

#include "geometry_structure.hpp"
#include "config_structure.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

class CGridAdaptation {
protected:
	unsigned long nPoint_new,	/*!< \brief Number of new points. */
	nElem_new;					/*!< \brief Number of new elements. */
	unsigned short nDim,	/*!< \brief Number of dimensions of the problem. */
	nVar;					/*!< \brief Number of variables in the problem. */
	double **ConsVar_Sol,	/*!< \brief Conservative variables (original solution). */
	**ConsVar_Res,			/*!< \brief Conservative variables (residual). */
	**ConsVar_Adapt;		/*!< \brief Conservative variables (adapted solution). */
	double **AdjVar_Sol,	/*!< \brief Adjoint variables (original solution). */
	**AdjVar_Res,			/*!< \brief Adjoint variables (residual). */
	**AdjVar_Adapt;			/*!< \brief Adjoint variables (adapted solution). */
	double **LinVar_Sol,	/*!< \brief Linear variables (original solution). */
	**LinVar_Res,			/*!< \brief Linear variables (residual). */
	**LinVar_Adapt;			/*!< \brief Linear variables (adapted solution). */
	double **Gradient,		/*!< \brief Gradient value. */
	**Gradient_Flow,		/*!< \brief Gradient of the flow variables. */
	**Gradient_Adj;			/*!< \brief Fradient of the adjoint variables. */
	double *Index;			/*!< \brief Adaptation index (indicates the value of the adaptation). */
	
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
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void GetLinSolution(CGeometry *geometry, CConfig *config);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void GetLinResidual(CGeometry *geometry, CConfig *config);
	
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
	 * \brief Do an adaptation of the computational grid on the two phase problem interphase.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] strength - Adaptation Strength.	 
	 */		
	void SetTwoPhase_Refinement(CGeometry *geometry, unsigned short strength);
	
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
	int CheckTriangleCode(bool *AdaptCode);
	
	/*! 
	 * \brief Find the adaptation code for each element in the fine grid.
	 * \param[in] AdaptCode - Edge combination to stablish the right elemeent division.
	 * \return Adaptation code for the element.
	 */	
	int CheckRectCode(bool *AdaptCode);
	
	/*! 
	 * \brief Find the adaptation code for each element in the fine grid.
	 * \param[in] AdaptCode - Edge combination to stablish the right elemeent division.
	 * \return Adaptation code for the element.
	 */	
	int CheckRectExtCode(bool *AdaptCode);
	
	/*! 
	 * \brief Find the adaptation code for each element in the fine grid.
	 * \param[in] AdaptCode - Edge combination to stablish the right elemeent division.
	 * \return Adaptation code for the element.
	 */	
	int CheckTetraCode(bool *AdaptCode);
	
	/*! 
	 * \brief Find the adaptation code for each element in the fine grid.
	 * \param[in] AdaptCode - Edge combination to stablish the right elemeent division.
	 * \return Adaptation code for the element.
	 */	
	int CheckHexaCode(bool *AdaptCode);
	
	/*! 
	 * \brief Find the adaptation code for each element in the fine grid.
	 * \param[in] AdaptCode - Edge combination to stablish the right elemeent division.
	 * \return Adaptation code for the element.
	 */	
	int CheckPyramCode(bool *AdaptCode);
	
	/*! 
	 * \brief Division pattern of the element.
	 * \param[in] code - number that identify the division.
	 * \param[in] nodes - Nodes that compose the element, including new nodes.
	 * \param[in] edges - Edges that compose the element.
	 * \param[out] Division - Division pattern. 
	 * \param[out] nPart - Number of new elements after the division. 
	 */	
	void TriangleDivision(int code, int *nodes, int *edges, int **Division, int *nPart);
	
	/*! 
	 * \brief Division pattern of the element.
	 * \param[in] code - number that identify the division.
	 * \param[in] nodes - Nodes that compose the element, including new nodes.
	 * \param[in] edges - Edges that compose the element.
	 * \param[out] Division - Division pattern. 
	 * \param[out] nPart - Number of new elements after the division. 
	 */	
	void RectDivision(int code, int *nodes, int **Division, int *nPart);
	
	/*! 
	 * \brief Division pattern of the element.
	 * \param[in] code - number that identify the division.
	 * \param[in] nodes - Nodes that compose the element, including new nodes.
	 * \param[in] edges - Edges that compose the element.
	 * \param[out] Division - Division pattern. 
	 * \param[out] nPart - Number of new elements after the division. 
	 */	
	void RectExtDivision(int code, int *nodes, int **Division, int *nPart);
	
	/*! 
	 * \brief Division pattern of the element.
	 * \param[in] code - number that identify the division.
	 * \param[in] nodes - Nodes that compose the element, including new nodes.
	 * \param[in] edges - Edges that compose the element.
	 * \param[out] Division - Division pattern. 
	 * \param[out] nPart - Number of new elements after the division. 
	 */	
	void TetraDivision(int code, int *nodes, int *edges, int **Division, int *nPart);
	
	/*! 
	 * \brief Division pattern of the element.
	 * \param[in] code - number that identify the division.
	 * \param[in] nodes - Nodes that compose the element, including new nodes.
	 * \param[in] edges - Edges that compose the element.
	 * \param[out] Division - Division pattern. 
	 * \param[out] nPart - Number of new elements after the division. 
	 */	
	void HexaDivision(int code, int *nodes, int **Division, int *nPart);
	
	/*! 
	 * \brief Division pattern of the element.
	 * \param[in] code - number that identify the division.
	 * \param[in] nodes - Nodes that compose the element, including new nodes.
	 * \param[in] edges - Edges that compose the element.
	 * \param[out] Division - Division pattern. 
	 * \param[out] nPart - Number of new elements after the division. 
	 */	
	void PyramDivision(int code, int *nodes, int **Division, int *nPart);
	
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
	void SetRestart_FlowSolution(CConfig *config, string mesh_flowfilename);
	
	/*! 
	 * \brief Write the restart file with the adapted grid.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] mesh_adjfilename - _________________________.
	 */		
	void SetRestart_AdjSolution(CConfig *config, string mesh_adjfilename);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] mesh_linfilename - _________________________.
	 */	
	void SetRestart_LinSolution(CConfig *config, string mesh_linfilename);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] max_elem - _________________________.
	 */	
	void SetSensorElem(CGeometry *geometry, CConfig *config, unsigned long max_elem);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] mesh_filename - _________________________.
	 */	
	void WriteAdaptSensor(CGeometry *geometry, char mesh_filename[200]);
};

#include "grid_adaptation_structure.inl"


