/*!
 * \file grid_adaptation_structure.hpp
 * \brief Header file for the adaptation subroutines.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.1.
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
	unsigned long nPoint_new,	/*!< \brief _______________. */
	nElem_new;					/*!< \brief _______________. */
	unsigned short nDim,	/*!< \brief _______________. */
	nVar;					/*!< \brief _______________. */
	double **ConsVar_Sol,	/*!< \brief _______________. */
	**ConsVar_Res,			/*!< \brief _______________. */
	**ConsVar_Adapt;		/*!< \brief _______________. */
	double **AdjVar_Sol,	/*!< \brief _______________. */
	**AdjVar_Res,			/*!< \brief _______________. */
	**AdjVar_Adapt;			/*!< \brief _______________. */
	double **LinVar_Sol,	/*!< \brief _______________. */
	**LinVar_Res,			/*!< \brief _______________. */
	**LinVar_Adapt;			/*!< \brief _______________. */
	double **Gradient,		/*!< \brief _______________. */
	**Gradient_Flow,		/*!< \brief _______________. */
	**Gradient_Adj;			/*!< \brief _______________. */
	double *Index;			/*!< \brief _______________. */
	bool change_geometry;	/*!< \brief _______________. */
	
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
	 * \brief Do a complete adaptation of the computational grid using a homothetic technique.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] geo_adapt - Geometrical definition of the adapted grid.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void SetHomothetic_Adaptation(CGeometry *geometry, CPhysicalGeometry *geo_adapt, CConfig *config);
	
	/*! 
	 * \brief Create a domain interface.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] geo_adapt - Geometrical definition of the adapted grid.
	 * \param[in] config - Definition of the particular problem.
	 */	
	void SetDomain_Interface(CGeometry *geometry, CPhysicalGeometry *geo_adapt, CConfig *config);
	
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
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] mesh_flowfilename - _________________________.
	 */		
	void SetReStart_FlowSolution(CGeometry *geometry, CConfig *config, string mesh_flowfilename);
	
	/*! 
	 * \brief Write the restart file with the adapted grid.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] mesh_adjfilename - _________________________.
	 */		
	void SetReStart_AdjSolution(CGeometry *geometry, CConfig *config, string mesh_adjfilename);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] mesh_linfilename - _________________________.
	 */	
	void SetReStart_LinSolution(CGeometry *geometry, CConfig *config, string mesh_linfilename);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] edge_code - _________________________.
	 * \param[in] new_tetra - _________________________.
	 * \param[in] nTetra - _________________________.
	 */		
	void SetTetraPattern (unsigned long  edge_code[6], bool new_tetra[8][10], unsigned short &nTetra);
	
	/*! 
	 * \brief Read the flow solution from the restart file.
	 * \param[in] edge_code - _________________________.
	 * \param[in] new_triangle - _________________________.
	 * \param[in] nTriangle - _________________________.
	 */	
	void SetTrianglePattern (unsigned long  edge_code[3], bool new_triangle[4][6], unsigned short &nTriangle);
	
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


