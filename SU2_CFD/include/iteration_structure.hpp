/*!
 * \file iteration_structure.hpp
 * \brief Headers of the main subroutines used by SU2_CFD.
 *        The subroutines and functions are in the <i>definition_structure.cpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 4.0.2 "Cardinal"
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

#include "../../Common/include/mpi_structure.hpp"

#include <ctime>

#include "solver_structure.hpp"
#include "integration_structure.hpp"
#include "output_structure.hpp"
#include "numerics_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../Common/include/config_structure.hpp"

using namespace std;

/*!
 * \class CIteration
 * \brief Parent class for defining a single iteration of a physics problem.
 * \author T. Economon
 * \version 4.0.2 "Cardinal"
 */
class CIteration {
protected:
  unsigned short nZone;	/*!< \brief Total number of zones in the problem. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CIteration(CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CIteration(void);
  
  /*!
   * \brief A virtual member.
   * \param[in] ??? - Description here.
   */
  virtual void Preprocess(COutput *output,
                          CIntegration ***integration_container,
                          CGeometry ***geometry_container,
                          CSolver ****solver_container,
                          CNumerics *****numerics_container,
                          CConfig **config_container,
                          CSurfaceMovement **surface_movement,
                          CVolumetricMovement **grid_movement,
                          CFreeFormDefBox*** FFDBox,
                          unsigned short val_iZone);
  
  /*!
   * \brief A virtual member.
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
  virtual void Iterate(COutput *output,
                       CIntegration ***integration_container,
                       CGeometry ***geometry_container,
                       CSolver ****solver_container,
                       CNumerics *****numerics_container,
                       CConfig **config_container,
                       CSurfaceMovement **surface_movement,
                       CVolumetricMovement **grid_movement,
                       CFreeFormDefBox*** FFDBox,
                       unsigned short val_iZone);
  
  /*!
   * \brief A virtual member.
   * \param[in] ??? - Description here.
   */
  virtual void Update(COutput *output,
                      CIntegration ***integration_container,
                      CGeometry ***geometry_container,
                      CSolver ****solver_container,
                      CNumerics *****numerics_container,
                      CConfig **config_container,
                      CSurfaceMovement **surface_movement,
                      CVolumetricMovement **grid_movement,
                      CFreeFormDefBox*** FFDBox,
                      unsigned short val_iZone);
  
  /*!
   * \brief A virtual member.
   * \param[in] ??? - Description here.
   */
  virtual void Monitor();
  
  /*!
   * \brief A virtual member.
   * \param[in] ??? - Description here.
   */
  virtual void Output();
  
  /*!
   * \brief A virtual member.
   * \param[in] ??? - Description here.
   */
  virtual void Postprocess();

};


/*!
 * \class CMeanFlowIteration
 * \brief Class for driving an iteration of the mean flow system.
 * \author T. Economon
 * \version 4.0.2 "Cardinal"
 */
class CMeanFlowIteration : public CIteration {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CMeanFlowIteration(CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CMeanFlowIteration(void);
  
  /*!
   * \brief Preprocessing to prepare for an iteration of the physics.
   * \param[in] ??? - Description here.
   */
  void Preprocess(COutput *output,
                  CIntegration ***integration_container,
                  CGeometry ***geometry_container,
                  CSolver ****solver_container,
                  CNumerics *****numerics_container,
                  CConfig **config_container,
                  CSurfaceMovement **surface_movement,
                  CVolumetricMovement **grid_movement,
                  CFreeFormDefBox*** FFDBox,
                  unsigned short val_iZone);
  
  /*!
   * \brief Perform a single iteration of the mean flow system.
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
  void Iterate(COutput *output,
               CIntegration ***integration_container,
               CGeometry ***geometry_container,
               CSolver ****solver_container,
               CNumerics *****numerics_container,
               CConfig **config_container,
               CSurfaceMovement **surface_movement,
               CVolumetricMovement **grid_movement,
               CFreeFormDefBox*** FFDBox,
               unsigned short val_iZone);
  
  /*!
   * \brief Updates the containers for the mean flow system.
   * \param[in] ??? - Description here.
   */
  void Update(COutput *output,
              CIntegration ***integration_container,
              CGeometry ***geometry_container,
              CSolver ****solver_container,
              CNumerics *****numerics_container,
              CConfig **config_container,
              CSurfaceMovement **surface_movement,
              CVolumetricMovement **grid_movement,
              CFreeFormDefBox*** FFDBox,
              unsigned short val_iZone);
  
  /*!
   * \brief Monitors the convergence and other metrics for the mean flow system.
   * \param[in] ??? - Description here.
   */
  void Monitor();
  
  /*!
   * \brief Outputs desired files and quantities for the mean flow system.
   * \param[in] ??? - Description here.
   */
  void Output();
  
  /*!
   * \brief Postprocesses the mean flow system before heading to another physics system or the next iteration.
   * \param[in] ??? - Description here.
   */
  void Postprocess();
  
  /*!
   * \brief Imposes a gust via the grid velocities.
   * \author S. Padron
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   */
  void SetWind_GustField(CConfig *config_container, CGeometry **geometry_container, CSolver ***solver_container);
  
  /*!
   * \brief Reads and initializes the vortex positions, strengths and gradient.
   * \author S. Padron
   * \param[in] nVortex - number of vortices.
   * \param[in] x0 - Vector of x-loc of the vortices.
   * \param[in] y0 - Vector of y-loc of the vortices.
   * \param[in] vort_strength - Vector of vortex strengths.
   * \param[in] r_core - Vector of vortex core size.
   */
  void InitializeVortexDistribution(unsigned long &nVortex, vector<su2double>& x0, vector<su2double>& y0, vector<su2double>& vort_strength, vector<su2double>& r_core);
  
  /*!
   * \brief compute and set mixing-plane quantities at the interface between two adjacent zone.
   * \author S. Vitale
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - zone of the problem.
   */
  void SetMixingPlane(CGeometry ***geometry_container, CSolver ****solver_container, CConfig **config_container, unsigned short iZone);

  /*!
   * \brief compute and set mixing-plane quantities at the interface between two adjacent zone.
   * \author S. Vitale
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] output - Definition of the output for the particular problem.
   * \param[in] iZone - zone of the problem.
   */
  void SetTurboPerformance(CGeometry ***geometry_container, CSolver ****solver_container, CConfig **config_container, COutput *output, unsigned short iZone);

  /*!
   * \brief compute and set mixing-plane quantities at the interface between two adjacent zone.
   * \author S. Vitale
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] output - Definition of the output for the particular problem.
   * \param[in] iZone - zone of the problem.
   */
  void SetMPITurboPerformance(CGeometry ***geometry_container, CSolver ****solver_container, CConfig **config_container, COutput *output, unsigned short iZone);


};

/*!
 * \class CWaveIteration
 * \brief Class for driving an iteration of the wave system.
 * \author T. Economon
 * \version 4.0.2 "Cardinal"
 */
class CWaveIteration : public CIteration {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CWaveIteration(CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CWaveIteration(void);
  
  /*!
   * \brief Preprocessing to prepare for an iteration of the physics.
   * \param[in] ??? - Description here.
   */
  void Preprocess(COutput *output,
                  CIntegration ***integration_container,
                  CGeometry ***geometry_container,
                  CSolver ****solver_container,
                  CNumerics *****numerics_container,
                  CConfig **config_container,
                  CSurfaceMovement **surface_movement,
                  CVolumetricMovement **grid_movement,
                  CFreeFormDefBox*** FFDBox,
                  unsigned short val_iZone);
  
  /*!
   * \brief Perform a single iteration of the wave system.
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
  void Iterate(COutput *output,
               CIntegration ***integration_container,
               CGeometry ***geometry_container,
               CSolver ****solver_container,
               CNumerics *****numerics_container,
               CConfig **config_container,
               CSurfaceMovement **surface_movement,
               CVolumetricMovement **grid_movement,
               CFreeFormDefBox*** FFDBox,
               unsigned short val_iZone);
  
  /*!
   * \brief Updates the containers for the wave system.
   * \param[in] ??? - Description here.
   */
  void Update(COutput *output,
              CIntegration ***integration_container,
              CGeometry ***geometry_container,
              CSolver ****solver_container,
              CNumerics *****numerics_container,
              CConfig **config_container,
              CSurfaceMovement **surface_movement,
              CVolumetricMovement **grid_movement,
              CFreeFormDefBox*** FFDBox,
              unsigned short val_iZone);
  
  /*!
   * \brief Monitors the convergence and other metrics for the wave system.
   * \param[in] ??? - Description here.
   */
  void Monitor();
  
  /*!
   * \brief Outputs desired files and quantities for the wave system.
   * \param[in] ??? - Description here.
   */
  void Output();
  
  /*!
   * \brief Postprocesses the wave system before heading to another physics system or the next iteration.
   * \param[in] ??? - Description here.
   */
  void Postprocess();
  
};

/*!
 * \class CHeatIteration
 * \brief Class for driving an iteration of the heat system.
 * \author T. Economon
 * \version 4.0.2 "Cardinal"
 */
class CHeatIteration : public CIteration {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CHeatIteration(CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CHeatIteration(void);
  
  /*!
   * \brief Preprocessing to prepare for an iteration of the physics.
   * \param[in] ??? - Description here.
   */
  void Preprocess(COutput *output,
                  CIntegration ***integration_container,
                  CGeometry ***geometry_container,
                  CSolver ****solver_container,
                  CNumerics *****numerics_container,
                  CConfig **config_container,
                  CSurfaceMovement **surface_movement,
                  CVolumetricMovement **grid_movement,
                  CFreeFormDefBox*** FFDBox,
                  unsigned short val_iZone);
  
  /*!
   * \brief Perform a single iteration of the heat system.
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
  void Iterate(COutput *output,
               CIntegration ***integration_container,
               CGeometry ***geometry_container,
               CSolver ****solver_container,
               CNumerics *****numerics_container,
               CConfig **config_container,
               CSurfaceMovement **surface_movement,
               CVolumetricMovement **grid_movement,
               CFreeFormDefBox*** FFDBox,
               unsigned short val_iZone);
  
  /*!
   * \brief Updates the containers for the heat system.
   * \param[in] ??? - Description here.
   */
  void Update(COutput *output,
              CIntegration ***integration_container,
              CGeometry ***geometry_container,
              CSolver ****solver_container,
              CNumerics *****numerics_container,
              CConfig **config_container,
              CSurfaceMovement **surface_movement,
              CVolumetricMovement **grid_movement,
              CFreeFormDefBox*** FFDBox,
              unsigned short val_iZone);
  
  /*!
   * \brief Monitors the convergence and other metrics for the heat system.
   * \param[in] ??? - Description here.
   */
  void Monitor();
  
  /*!
   * \brief Outputs desired files and quantities for the heat system.
   * \param[in] ??? - Description here.
   */
  void Output();
  
  /*!
   * \brief Postprocesses the heat system before heading to another physics system or the next iteration.
   * \param[in] ??? - Description here.
   */
  void Postprocess();
  
};

/*!
 * \class CPoissonIteration
 * \brief Class for driving an iteration of the poisson system.
 * \author T. Economon
 * \version 4.0.2 "Cardinal"
 */
class CPoissonIteration : public CIteration {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CPoissonIteration(CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CPoissonIteration(void);
  
  /*!
   * \brief Preprocessing to prepare for an iteration of the physics.
   * \param[in] ??? - Description here.
   */
  void Preprocess(COutput *output,
                  CIntegration ***integration_container,
                  CGeometry ***geometry_container,
                  CSolver ****solver_container,
                  CNumerics *****numerics_container,
                  CConfig **config_container,
                  CSurfaceMovement **surface_movement,
                  CVolumetricMovement **grid_movement,
                  CFreeFormDefBox*** FFDBox,
                  unsigned short val_iZone);
  
  /*!
   * \brief Perform a single iteration of the poisson system.
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
  void Iterate(COutput *output,
               CIntegration ***integration_container,
               CGeometry ***geometry_container,
               CSolver ****solver_container,
               CNumerics *****numerics_container,
               CConfig **config_container,
               CSurfaceMovement **surface_movement,
               CVolumetricMovement **grid_movement,
               CFreeFormDefBox*** FFDBox,
               unsigned short val_iZone);
  
  /*!
   * \brief Updates the containers for the poisson system.
   * \param[in] ??? - Description here.
   */
  void Update(COutput *output,
              CIntegration ***integration_container,
              CGeometry ***geometry_container,
              CSolver ****solver_container,
              CNumerics *****numerics_container,
              CConfig **config_container,
              CSurfaceMovement **surface_movement,
              CVolumetricMovement **grid_movement,
              CFreeFormDefBox*** FFDBox,
              unsigned short val_iZone);
  
  /*!
   * \brief Monitors the convergence and other metrics for the poisson system.
   * \param[in] ??? - Description here.
   */
  void Monitor();
  
  /*!
   * \brief Outputs desired files and quantities for the poisson system.
   * \param[in] ??? - Description here.
   */
  void Output();
  
  /*!
   * \brief Postprocesses the poisson system before heading to another physics system or the next iteration.
   * \param[in] ??? - Description here.
   */
  void Postprocess();
  
};

/*!
 * \class CFEAIteration
 * \brief Class for driving an iteration of the FEA system.
 * \author T. Economon
 * \version 4.0.2 "Cardinal"
 */
class CFEAIteration : public CIteration {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CFEAIteration(CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CFEAIteration(void);
  
  /*!
   * \brief Preprocessing to prepare for an iteration of the physics.
   * \param[in] ??? - Description here.
   */
  void Preprocess(COutput *output,
                  CIntegration ***integration_container,
                  CGeometry ***geometry_container,
                  CSolver ****solver_container,
                  CNumerics *****numerics_container,
                  CConfig **config_container,
                  CSurfaceMovement **surface_movement,
                  CVolumetricMovement **grid_movement,
                  CFreeFormDefBox*** FFDBox,
                  unsigned short val_iZone);
  
  /*!
   * \brief Perform a single iteration of the FEA system.
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
  void Iterate(COutput *output,
               CIntegration ***integration_container,
               CGeometry ***geometry_container,
               CSolver ****solver_container,
               CNumerics *****numerics_container,
               CConfig **config_container,
               CSurfaceMovement **surface_movement,
               CVolumetricMovement **grid_movement,
               CFreeFormDefBox*** FFDBox,
               unsigned short val_iZone);
  
  /*!
   * \brief Updates the containers for the FEA system.
   * \param[in] ??? - Description here.
   */
  void Update(COutput *output,
              CIntegration ***integration_container,
              CGeometry ***geometry_container,
              CSolver ****solver_container,
              CNumerics *****numerics_container,
              CConfig **config_container,
              CSurfaceMovement **surface_movement,
              CVolumetricMovement **grid_movement,
              CFreeFormDefBox*** FFDBox,
              unsigned short val_iZone);
  
  /*!
   * \brief Monitors the convergence and other metrics for the FEA system.
   * \param[in] ??? - Description here.
   */
  void Monitor();
  
  /*!
   * \brief Outputs desired files and quantities for the FEA system.
   * \param[in] ??? - Description here.
   */
  void Output();
  
  /*!
   * \brief Postprocesses the FEA system before heading to another physics system or the next iteration.
   * \param[in] ??? - Description here.
   */
  void Postprocess();
  
};

/*!
 * \class CAdjMeanFlowIteration
 * \brief Class for driving an iteration of the adjoint mean flow system.
 * \author T. Economon
 * \version 4.0.2 "Cardinal"
 */
class CAdjMeanFlowIteration : public CIteration {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CAdjMeanFlowIteration(CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAdjMeanFlowIteration(void);
  
  /*!
   * \brief Preprocessing to prepare for an iteration of the physics.
   * \param[in] ??? - Description here.
   */
  void Preprocess(COutput *output,
                  CIntegration ***integration_container,
                  CGeometry ***geometry_container,
                  CSolver ****solver_container,
                  CNumerics *****numerics_container,
                  CConfig **config_container,
                  CSurfaceMovement **surface_movement,
                  CVolumetricMovement **grid_movement,
                  CFreeFormDefBox*** FFDBox,
                  unsigned short val_iZone);
  
  /*!
   * \brief Perform a single iteration of the adjoint mean flow system.
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
  void Iterate(COutput *output,
               CIntegration ***integration_container,
               CGeometry ***geometry_container,
               CSolver ****solver_container,
               CNumerics *****numerics_container,
               CConfig **config_container,
               CSurfaceMovement **surface_movement,
               CVolumetricMovement **grid_movement,
               CFreeFormDefBox*** FFDBox,
               unsigned short val_iZone);
  
  /*!
   * \brief Updates the containers for the adjoint mean flow system.
   * \param[in] ??? - Description here.
   */
  void Update(COutput *output,
              CIntegration ***integration_container,
              CGeometry ***geometry_container,
              CSolver ****solver_container,
              CNumerics *****numerics_container,
              CConfig **config_container,
              CSurfaceMovement **surface_movement,
              CVolumetricMovement **grid_movement,
              CFreeFormDefBox*** FFDBox,
              unsigned short val_iZone);
  
  /*!
   * \brief Monitors the convergence and other metrics for the adjoint mean flow system.
   * \param[in] ??? - Description here.
   */
  void Monitor();
  
  /*!
   * \brief Outputs desired files and quantities for the adjoint mean flow system.
   * \param[in] ??? - Description here.
   */
  void Output();
  
  /*!
   * \brief Postprocesses the adjoint mean flow system before heading to another physics system or the next iteration.
   * \param[in] ??? - Description here.
   */
  void Postprocess();
  
};

/*!
 * \class CDiscAdjMeanFlowIteration
 * \brief Class for driving an iteration of the discrete adjoint mean flow system.
 * \author T. Economon
 * \version 4.0.2 "Cardinal"
 */
class CDiscAdjMeanFlowIteration : public CIteration {

private:

  CIteration* meanflow_iteration; /*!< \brief Pointer to the mean flow iteration class. */
  unsigned short CurrentRecording; /*!< \brief Stores the current status of the recording. */
  bool turbulent;       /*!< \brief Stores the turbulent flag. */

  enum RECORDING{
    NONE = 0,      /*!< \brief Indicates that nothing is recorded. */
    FLOW_VARIABLES = 1, /*!< \brief Indicates that the current recording
                                    can be used to compute the gradients with respect
                                    to the conservative flow variables. */
    GEOMETRY_VARIABLES = 2 /*!< \brief Indicates that the current recording
                                       can be used to compute the gradients with respect
                                       to the geometry variables. */
  };


public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjMeanFlowIteration(CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjMeanFlowIteration(void);
  
  /*!
   * \brief Preprocessing to prepare for an iteration of the physics.
   * \param[in] ??? - Description here.
   */
  void Preprocess(COutput *output,
                  CIntegration ***integration_container,
                  CGeometry ***geometry_container,
                  CSolver ****solver_container,
                  CNumerics *****numerics_container,
                  CConfig **config_container,
                  CSurfaceMovement **surface_movement,
                  CVolumetricMovement **grid_movement,
                  CFreeFormDefBox*** FFDBox,
                  unsigned short val_iZone);
  
  /*!
   * \brief Perform a single iteration of the adjoint mean flow system.
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
  void Iterate(COutput *output,
               CIntegration ***integration_container,
               CGeometry ***geometry_container,
               CSolver ****solver_container,
               CNumerics *****numerics_container,
               CConfig **config_container,
               CSurfaceMovement **surface_movement,
               CVolumetricMovement **grid_movement,
               CFreeFormDefBox*** FFDBox,
               unsigned short val_iZone);
  
  /*!
   * \brief Updates the containers for the discrete adjoint mean flow system.
   * \param[in] ??? - Description here.
   */
  void Update(COutput *output,
              CIntegration ***integration_container,
              CGeometry ***geometry_container,
              CSolver ****solver_container,
              CNumerics *****numerics_container,
              CConfig **config_container,
              CSurfaceMovement **surface_movement,
              CVolumetricMovement **grid_movement,
              CFreeFormDefBox*** FFDBox,
              unsigned short val_iZone);
  
  /*!
   * \brief Monitors the convergence and other metrics for the discrete adjoint mean flow system.
   * \param[in] ??? - Description here.
   */
  void Monitor();
  
  /*!
   * \brief Outputs desired files and quantities for the discrete adjoint mean flow system.
   * \param[in] ??? - Description here.
   */
  void Output();
  
  /*!
   * \brief Postprocesses the discrete adjoint mean flow system before heading to another physics system or the next iteration.
   * \param[in] ??? - Description here.
   */
  void Postprocess(); 

  /*!
   * \brief Registers all input variables of the mean flow iteration.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] kind_recording - Kind of recording, either FLOW_VARIABLES or GEOMETRY_VARIABLES
   */
  void RegisterInput(CSolver ****solver_container, CGeometry*** geometry_container, CConfig** config_container, unsigned short iZone, unsigned short kind_recording);

  /*!
   * \brief Registers all output variables of the mean flow iteration.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   */
  void RegisterOutput(CSolver ****solver_container, CGeometry*** geometry_container, CConfig** config_container, unsigned short iZone);

  /*!
   * \brief Initializes the adjoints of the output variables of the meanflow iteration.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   */
  void InitializeAdjoint(CSolver ****solver_container, CGeometry*** geometry_container, CConfig** config_container, unsigned short iZone);


  /*!
   * \brief Record a single iteration of the direct mean flow system.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] kind_recording - The kind of recording (geometry or flow).
   */
  void SetRecording(COutput *output,
                      CIntegration ***integration_container,
                      CGeometry ***geometry_container,
                      CSolver ****solver_container,
                      CNumerics *****numerics_container,
                      CConfig **config_container,
                      CSurfaceMovement **surface_movement,
                      CVolumetricMovement **grid_movement,
                      CFreeFormDefBox*** FFDBox,
                      unsigned short val_iZone,
                      unsigned short kind_recording);
};


/*!
 * \brief Iteration function for Fluid-Structure Interaction applications.
 * \author F. Palacios, R. Sanchez.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxes of the problem.
 * \param[in] ExtIter - Current physical time iteration number.
 * \param[in] nFluidIt - Number of fluid iterations within a fixed time step.
 */
void FluidStructureIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
														 CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container, 
														 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
														 unsigned long iFluidIt, unsigned long nFluidIt);

/*!
 * \brief Updates the positions and grid velocities for dynamic meshes between physical time steps.
 * \author T. Economon
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] surface_movement - Surface movement classes of the problem.
 * \param[in] grid_movement - Volume grid movement classes of the problem.
 * \param[in] FFDBox - FFD FFDBoxes of the problem.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iZone - Index of the zone.
 * \param[in] IntIter - Current sudo time iteration number.
 * \param[in] ExtIter - Current physical time iteration number.
 */
void SetGrid_Movement(CGeometry **geometry_container, CSurfaceMovement *surface_movement, 
                      CVolumetricMovement *grid_movement, CFreeFormDefBox **FFDBox,
                      CSolver ***solver_container, CConfig *config_container, unsigned short iZone, unsigned long IntIter, unsigned long ExtIter);
