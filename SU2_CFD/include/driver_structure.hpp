/*!
 * \file driver_structure.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 4.2.0 "Cardinal"
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
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../../Common/include/mpi_structure.hpp"
#include "iteration_structure.hpp"
#include "solver_structure.hpp"
#include "integration_structure.hpp"
#include "output_structure.hpp"
#include "numerics_structure.hpp"
#include "transfer_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/interpolation_structure.hpp"

using namespace std;

/*! 
 * \class CDriver
 * \brief Parent class for driving an iteration of a single or multi-zone problem.
 * \author T. Economon
 * \version 4.2.0 "Cardinal"
 */
class CDriver {
protected:
  unsigned short nZone;	/*!< \brief Total number of zones in the problem. */

public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - Total number of zones.
	 */
  CDriver(CIteration **iteration_container,
          CSolver ****solver_container,
          CGeometry ***geometry_container,
          CIntegration ***integration_container,
          CNumerics *****numerics_container,
          CInterpolator ***interpolator_container,
          CTransfer ***transfer_container,
          CConfig **config,
          unsigned short val_nZone,
          unsigned short val_nDim);
	
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CDriver(void);

	/*!
	 * \brief A virtual member.
	 * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
	 */
  
  virtual void Run(CIteration **iteration_container,
                   COutput *output,
                   CIntegration ***integration_container,
                   CGeometry ***geometry_container,
                   CSolver ****solver_container,
                   CNumerics *****numerics_container,
                   CConfig **config_container,
                   CSurfaceMovement **surface_movement,
                   CVolumetricMovement **grid_movement,
                   CFreeFormDefBox*** FFDBox,
                   CInterpolator ***interpolator_container,
                   CTransfer ***transfer_container){
  };
  /*!
   * \brief Definition of the physics iteration class or within a single zone.
   * \param[in] iteration_container - Pointer to the iteration container to be instantiated.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   */
  void Iteration_Preprocessing(CIteration **iteration_container, CConfig **config, unsigned short iZone);

  /*!
   * \brief Definition and allocation of all solution classes.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Preprocessing(CSolver ***solver_container, CGeometry **geometry, CConfig *config);

  /*!
   * \brief Definition and allocation of all solution classes.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Postprocessing(CSolver ***solver_container, CGeometry **geometry, CConfig *config);

  /*!
   * \brief Definition and allocation of all integration classes.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Integration_Preprocessing(CIntegration **integration_container, CGeometry **geometry, CConfig *config);

  /*!
   * \brief Definition and allocation of all integration classes.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Integration_Postprocessing(CIntegration **integration_container, CGeometry **geometry, CConfig *config);

  /*!
   * \brief Definition and allocation of all interface classes.
   * \param[in] transfer_container - Definition of the transfer of information and the physics involved in the interface.
   * \param[in] interpolator_container - Definition of the interpolation method between non-matching discretizations of the interface.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] nZone - Total number of zones.
   * \param[in] nDim - Total number of dimensions.
   */
  void Interface_Preprocessing(CTransfer ***transfer_container, CInterpolator ***interpolator_container,
		  	  	  	  	  	  	   CGeometry ***geometry_container, CConfig **config_container, CSolver ****solver_container,
		  	  	  	  	  	  	   unsigned short nZone, unsigned short nDim);


  /*!
   * \brief Definition and allocation of all solver classes.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Numerics_Preprocessing(CNumerics ****numerics_container, CSolver ***solver_container, CGeometry **geometry, CConfig *config);


  /*!
   * \brief Definition and allocation of all solver classes.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Numerics_Postprocessing(CNumerics ****numerics_container, CSolver ***solver_container, CGeometry **geometry, CConfig *config);


  /*!
   * \brief Deallocation routine
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - Total number of zones.
   */
  void Postprocessing(CIteration **iteration_container,
                   CSolver ****solver_container,
                   CGeometry ***geometry_container,
                   CIntegration ***integration_container,
                   CNumerics *****numerics_container,
                   CInterpolator ***interpolator_container,
                   CTransfer ***transfer_container,
                   CConfig **config_container,
                   unsigned short val_nZone);

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone in which the displacements will be predicted.
   * \param[in] targetZone - zone which receives the predicted displacements.
   */
  virtual void Predict_Displacements(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 unsigned short donorZone, unsigned short targetZone){};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone in which the tractions will be predicted.
   * \param[in] targetZone - zone which receives the predicted traction.
   */
  virtual void Predict_Tractions(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 unsigned short donorZone, unsigned short targetZone){};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone in which the displacements will be transferred.
   * \param[in] targetZone - zone which receives the tractions transferred.
   */
  virtual void Transfer_Displacements(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 CTransfer ***transfer_container, unsigned short donorZone, unsigned short targetZone){};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone from which the tractions will be transferred.
   * \param[in] targetZone - zone which receives the tractions transferred.
   */
  virtual void Transfer_Tractions(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 CTransfer ***transfer_container, unsigned short donorZone, unsigned short targetZone){};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iFSIIter - Fluid-Structure Interaction subiteration.
   */
  virtual void Relaxation_Displacements(COutput *output, CGeometry ***geometry_container, CSolver ****solver_container,
			CConfig **config_container, unsigned short donorZone, unsigned short targetZone, unsigned long iFSIIter){};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iFSIIter - Fluid-Structure Interaction subiteration.
   */
  virtual void Relaxation_Tractions(COutput *output, CGeometry ***geometry_container, CSolver ****solver_container,
			CConfig **config_container, unsigned short donorZone, unsigned short targetZone, unsigned long iFSIIter){};

  /*!
   * \brief A virtual member.
   * \param[in] zoneFlow - zone of the flow equations.
   * \param[in] zoneStruct - zone of the structural equations.
   */
  virtual void Update(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 CTransfer ***transfer_container, unsigned short zoneFlow, unsigned short zoneStruct){};

};
/*!
 * \class CSingleZoneDriver
 * \brief Class for driving an iteration of the physics within a single zone.
 * \author T. Economon
 * \version 4.2.0 "Cardinal"
 */
class CSingleZoneDriver : public CDriver {
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] iteration_container - Container vector with all the iteration methods.
	 * \param[in] solver_container - Container vector with all the solutions.
	 * \param[in] geometry_container - Geometrical definition of the problem.
	 * \param[in] integration_container - Container vector with all the integration methods.
	 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] val_nZone - Total number of zones.
	 */
  CSingleZoneDriver(CIteration **iteration_container,
                    CSolver ****solver_container,
                    CGeometry ***geometry_container,
                    CIntegration ***integration_container,
                    CNumerics *****numerics_container,
                    CInterpolator ***interpolator_container,
                    CTransfer ***transfer_container,
                    CConfig **config,
                    unsigned short val_nZone,
                    unsigned short val_nDim);
	
	/*!
	 * \brief Destructor of the class.
	 */
	~CSingleZoneDriver(void);
	
	/*! 
	 * \brief Run a single iteration of the physics within a single zone.
	 * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
	 */
  
  void Run(CIteration **iteration_container,
           COutput *output,
           CIntegration ***integration_container,
           CGeometry ***geometry_container,
           CSolver ****solver_container,
           CNumerics *****numerics_container,
           CConfig **config_container,
           CSurfaceMovement **surface_movement,
           CVolumetricMovement **grid_movement,
           CFreeFormDefBox*** FFDBox,
           CInterpolator ***interpolator_container,
           CTransfer ***transfer_container);


};


/*!
 * \class CMultiZoneDriver
 * \brief Class for driving an iteration of the physics within multiple zones.
 * \author T. Economon
 * \version 4.2.0 "Cardinal"
 */
class CMultiZoneDriver : public CDriver {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - Total number of zones.
   */
  CMultiZoneDriver(CIteration **iteration_container,
                   CSolver ****solver_container,
                   CGeometry ***geometry_container,
                   CIntegration ***integration_container,
                   CNumerics *****numerics_container,
                   CInterpolator ***interpolator_container,
                   CTransfer ***transfer_container,
                   CConfig **config,
                   unsigned short val_nZone,
                   unsigned short val_nDim);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CMultiZoneDriver(void);
  
  /*!
   * \brief Run a single iteration of the physics within multiple zones.
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   */
  
  void Run(CIteration **iteration_container,
           COutput *output,
           CIntegration ***integration_container,
           CGeometry ***geometry_container,
           CSolver ****solver_container,
           CNumerics *****numerics_container,
           CConfig **config_container,
           CSurfaceMovement **surface_movement,
           CVolumetricMovement **grid_movement,
           CFreeFormDefBox*** FFDBox,
           CInterpolator ***interpolator_container,
           CTransfer ***transfer_container);

};


/*!
 * \class CSpectralDriver
 * \brief Class for driving an iteration of a spectral method problem using multiple zones.
 * \author T. Economon
 * \version 4.2.0 "Cardinal"
 */
class CSpectralDriver : public CDriver {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - Total number of zones.
   */
  CSpectralDriver(CIteration **iteration_container,
                   CSolver ****solver_container,
                   CGeometry ***geometry_container,
                   CIntegration ***integration_container,
                   CNumerics *****numerics_container,
                   CInterpolator ***interpolator_container,
                   CTransfer ***transfer_container,
                   CConfig **config,
                   unsigned short val_nZone,
                   unsigned short val_nDim);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSpectralDriver(void);
  
  /*!
   * \brief Run a single iteration of a spectral method problem.
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   */
  
  void Run(CIteration **iteration_container,
          COutput *output,
          CIntegration ***integration_container,
          CGeometry ***geometry_container,
          CSolver ****solver_container,
          CNumerics *****numerics_container,
          CConfig **config_container,
          CSurfaceMovement **surface_movement,
          CVolumetricMovement **grid_movement,
          CFreeFormDefBox*** FFDBox,
          CInterpolator ***interpolator_container,
          CTransfer ***transfer_container);
  
  /*!
   * \brief Computation and storage of the time spectral source terms.
   * \author T. Economon, K. Naik
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] nZone - Total number of zones (periodic instances).
   * \param[in] iZone - Current zone number.
   */
  void SetTimeSpectral(CGeometry ***geometry_container, CSolver ****solver_container,
                       CConfig **config_container, unsigned short nZone, unsigned short iZone);
  
  /*!
   * \brief Computation of the Time-Spectral operator matrix.
   * \author K. Naik
   * \param[in] D - su2double pointer to the operator matrix.
   * \param[in] nZone - Total number of zones (periodic instances).
   */
  void ComputeTimeSpectral_Operator(su2double **D, su2double period, unsigned short nZone);
  
  /*!
   * \brief Computation and storage of the time-spectral mesh velocities.
   * \author K. Naik, T. Economon
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] nZone - Total number of zones (periodic instances).
   */
  void SetTimeSpectral_Velocities(CGeometry ***geometry_container,
                                  CConfig **config_container, unsigned short nZone);
  
};


/*!
 * \class CFSIDriver
 * \brief Class for driving a BGS iteration for a fluid-structure interaction problem in multiple zones.
 * \author R. Sanchez.
 * \version 4.2.0 "Cardinal"
 */
class CFSIDriver : public CDriver {
public:

	/*!
	 * \brief Constructor of the class.
	 * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - Total number of zones.
	 */
  CFSIDriver(CIteration **iteration_container,
             CSolver ****solver_container,
             CGeometry ***geometry_container,
             CIntegration ***integration_container,
             CNumerics *****numerics_container,
             CInterpolator ***interpolator_container,
             CTransfer ***transfer_container,
             CConfig **config,
             unsigned short val_nZone,
             unsigned short val_nDim);

	/*!
	 * \brief Destructor of the class.
	 */
	~CFSIDriver(void);

	/*!
	 * \brief Run a Block Gauss-Seidel iteration of the FSI problem.
	 * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
	 */
  
  void Run(CIteration **iteration_container,
           COutput *output,
           CIntegration ***integration_container,
           CGeometry ***geometry_container,
           CSolver ****solver_container,
           CNumerics *****numerics_container,
           CConfig **config_container,
           CSurfaceMovement **surface_movement,
           CVolumetricMovement **grid_movement,
           CFreeFormDefBox*** FFDBox,
           CInterpolator ***interpolator_container,
           CTransfer ***transfer_container);

  /*!
   * \brief Predict the structural displacements to pass them into the fluid solver on a BGS implementation.
   * \param[in] donorZone - zone in which the displacements will be predicted.
   * \param[in] targetZone - zone which receives the predicted displacements.
   */
  void Predict_Displacements(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 unsigned short donorZone, unsigned short targetZone);

  /*!
   * \brief Predict the fluid tractions to pass them into the structural solver on a BGS implementation.
   * \param[in] donorZone - zone in which the tractions will be predicted.
   * \param[in] targetZone - zone which receives the predicted traction.
   */
  void Predict_Tractions(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 unsigned short donorZone, unsigned short targetZone);

  /*!
   * \brief Transfer the displacements computed on the structural solver into the fluid solver.
   * \param[in] donorZone - zone in which the displacements will be transferred.
   * \param[in] targetZone - zone which receives the tractions transferred.
   */
  void Transfer_Displacements(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 CTransfer ***transfer_container, unsigned short donorZone, unsigned short targetZone);

  /*!
   * \brief Transfer the tractions computed on the fluid solver into the structural solver.
   * \param[in] donorZone - zone from which the tractions will be transferred.
   * \param[in] targetZone - zone which receives the tractions transferred.
   */
  void Transfer_Tractions(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 CTransfer ***transfer_container, unsigned short donorZone, unsigned short targetZone);

  /*!
   * \brief Apply a relaxation method into the computed displacements.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iFSIIter - Fluid-Structure Interaction subiteration.
   */
  void Relaxation_Displacements(COutput *output, CGeometry ***geometry_container, CSolver ****solver_container,
			CConfig **config_container, unsigned short donorZone, unsigned short targetZone, unsigned long iFSIIter);

  /*!
   * \brief Apply a relaxation method into the computed tractions.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iFSIIter - Fluid-Structure Interaction subiteration.
   */
  void Relaxation_Tractions(COutput *output, CGeometry ***geometry_container, CSolver ****solver_container,
			CConfig **config_container, unsigned short donorZone, unsigned short targetZone, unsigned long iFSIIter);

  /*!
   * \brief Enforce the coupling condition at the end of the time step
   * \param[in] zoneFlow - zone of the flow equations.
   * \param[in] zoneStruct - zone of the structural equations.
   */
  void Update(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 CTransfer ***transfer_container, unsigned short zoneFlow, unsigned short zoneStruct);

};
