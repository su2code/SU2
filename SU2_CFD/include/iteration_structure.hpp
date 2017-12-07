/*!
 * \file iteration_structure.hpp
 * \brief Headers of the main subroutines used by SU2_CFD.
 *        The subroutines and functions are in the <i>definition_structure.cpp</i> file.
 * \author F. Palacios, T. Economon
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

#include "../../Common/include/mpi_structure.hpp"

#include <ctime>

#include "solver_structure.hpp"
#include "integration_structure.hpp"
#include "output_structure.hpp"
#include "numerics_structure.hpp"
#include "transfer_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../Common/include/config_structure.hpp"

using namespace std;

/*!
 * \class CIteration
 * \brief Parent class for defining a single iteration of a physics problem.
 * \author T. Economon
 * \version 5.0.0 "Raven"
 */
class CIteration {
protected:
  unsigned short nZone;  /*!< \brief Total number of zones in the problem. */
  
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
  virtual void SetGrid_Movement(CGeometry ***geometry_container, CSurfaceMovement **surface_movement,
                      CVolumetricMovement **grid_movement, CFreeFormDefBox ***FFDBox,
                      CSolver ****solver_container, CConfig **config_container,
                      unsigned short val_iZone, unsigned long IntIter, unsigned long ExtIter);
  
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
  virtual void Postprocess(COutput *output,
                            CIntegration ***integration_container,
                            CGeometry ***geometry_container,
                            CSolver ****solver_container,
                            CNumerics *****numerics_container,
                            CConfig **config_container,
                            CSurfaceMovement **surface_movement,
                            CVolumetricMovement **grid_movement,
                            CFreeFormDefBox*** FFDBox,
                            unsigned short val_iZone);

  virtual void InitializeAdjoint(CSolver ****solver_container,
                                 CGeometry ***geometry_container,
                                 CConfig **config_container,
                                 unsigned short iZone){}

  virtual void InitializeAdjoint_CrossTerm(CSolver ****solver_container,
                                 CGeometry ***geometry_container,
                                 CConfig **config_container,
                                 unsigned short iZone){}

  virtual void RegisterInput(CSolver ****solver_container,
                             CGeometry*** geometry_container,
                             CConfig** config_container,
                             unsigned short iZone,
                             unsigned short kind_recording){}

  virtual void SetDependencies(CSolver ****solver_container,
                               CGeometry ***geometry_container,
                               CConfig **config_container,
                               unsigned short iZone,
                               unsigned short kind_recording){}

  virtual void SetDependencies(CSolver ****solver_container,
                               CGeometry ***geometry_container,
                               CNumerics *****numerics_container,
                               CConfig **config_container,
                               unsigned short iZone,
                               unsigned short kind_recording){}

  virtual void RegisterOutput(CSolver ****solver_container,
                              CGeometry*** geometry_container,
                              CConfig** config_container,
                              COutput* output,
                              unsigned short iZone){}

  virtual void LoadUnsteady_Solution(CGeometry ***geometry_container,
                                         CSolver ****solver_container,
                                         CConfig **config_container,
                                         unsigned short val_iZone,
                                         int val_DirectIter){}

  virtual void LoadDynamic_Solution(CGeometry ***geometry_container,
                                        CSolver ****solver_container,
                                        CConfig **config_container,
                                        unsigned short val_iZone,
                                        int val_DirectIter){}

};


/*!
 * \class CFluidIteration
 * \brief Class for driving an iteration of the fluid system.
 * \author T. Economon
 * \version 5.0.0 "Raven"
 */
class CFluidIteration : public CIteration {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CFluidIteration(CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CFluidIteration(void);
  
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
   * \brief Perform a single iteration of the fluid system.
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
   * \brief Updates the containers for the fluid system.
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
   * \brief Monitors the convergence and other metrics for the fluid system.
   * \param[in] ??? - Description here.
   */
  void Monitor();
  
  /*!
   * \brief Outputs desired files and quantities for the fluid system.
   */
  void Output();
  
  /*!
   * \brief Postprocesses the fluid system before heading to another physics system or the next iteration.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   */
  void Postprocess(COutput *output,
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
  
};



/*!
 * \class CTurboIteration
 * \brief Class for driving an iteration for turbomachinery simulation.
 * \author T. Economon
 * \version 5.0.0 "Raven"
 */
class CTurboIteration : public CFluidIteration {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CTurboIteration(CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurboIteration(void);

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
   * \brief Postprocesses the fluid system before heading to another physics system or the next iteration.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   */
  void Postprocess(COutput *output,
                   CIntegration ***integration_container,
                   CGeometry ***geometry_container,
                   CSolver ****solver_container,
                   CNumerics *****numerics_container,
                   CConfig **config_container,
                   CSurfaceMovement **surface_movement,
                   CVolumetricMovement **grid_movement,
                   CFreeFormDefBox*** FFDBox,
                   unsigned short val_iZone);


};

/*!
 * \class CWaveIteration
 * \brief Class for driving an iteration of the wave system.
 * \author T. Economon
 * \version 5.0.0 "Raven"
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
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   * \param[in] val_iZone - zone of the problem.
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
   * \param[in] val_iZone - zone of the problem.
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
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   * \param[in] val_iZone - zone of the problem.
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
   */
  void Monitor();
  
  /*!
   * \brief Outputs desired files and quantities for the wave system.
   */
  void Output();
  
  /*!
   * \brief Postprocess ???.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   */
  void Postprocess(COutput *output,
                   CIntegration ***integration_container,
                   CGeometry ***geometry_container,
                   CSolver ****solver_container,
                   CNumerics *****numerics_container,
                   CConfig **config_container,
                   CSurfaceMovement **surface_movement,
                   CVolumetricMovement **grid_movement,
                   CFreeFormDefBox*** FFDBox,
                   unsigned short val_iZone);
  
};

/*!
 * \class CHeatIteration
 * \brief Class for driving an iteration of the heat system.
 * \author T. Economon
 * \version 5.0.0 "Raven"
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
   * \brief Postprocess ???.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   */
  void Postprocess(COutput *output,
                   CIntegration ***integration_container,
                   CGeometry ***geometry_container,
                   CSolver ****solver_container,
                   CNumerics *****numerics_container,
                   CConfig **config_container,
                   CSurfaceMovement **surface_movement,
                   CVolumetricMovement **grid_movement,
                   CFreeFormDefBox*** FFDBox,
                   unsigned short val_iZone);
  
};

/*!
 * \class CPoissonIteration
 * \brief Class for driving an iteration of the poisson system.
 * \author T. Economon
 * \version 5.0.0 "Raven"
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
   * \brief Postprocess ???.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   */
  void Postprocess(COutput *output,
                   CIntegration ***integration_container,
                   CGeometry ***geometry_container,
                   CSolver ****solver_container,
                   CNumerics *****numerics_container,
                   CConfig **config_container,
                   CSurfaceMovement **surface_movement,
                   CVolumetricMovement **grid_movement,
                   CFreeFormDefBox*** FFDBox,
                   unsigned short val_iZone);
  
};

/*!
 * \class CFEM_StructuralAnalysis
 * \brief Class for driving an iteration of structural analysis.
 * \author R. Sanchez
 * \version 4.0.1 "Cardinal"
 */
class CFEM_StructuralAnalysis : public CIteration {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CFEM_StructuralAnalysis(CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_StructuralAnalysis(void);

  /*!
   * \brief Preprocessing to prepare for an iteration of the physics.
   * \param[in] ??? - Description here.
   */
  void Preprocess();
  using CIteration::Preprocess;


  /*!
   * \brief Perform a single iteration for structural analysis using the Finite Element Method.
   * \author R. Sanchez.
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
   * \brief Updates the containers for the FEM system.
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
   * \brief Monitors the convergence and other metrics for the FEM system.
   * \param[in] ??? - Description here.
   */
  void Monitor();

  /*!
   * \brief Outputs desired files and quantities for the FEM system.
   * \param[in] ??? - Description here.
   */
  void Output();

  /*!
   * \brief Postprocess ???.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   */
  void Postprocess(COutput *output,
                   CIntegration ***integration_container,
                   CGeometry ***geometry_container,
                   CSolver ****solver_container,
                   CNumerics *****numerics_container,
                   CConfig **config_container,
                   CSurfaceMovement **surface_movement,
                   CVolumetricMovement **grid_movement,
                   CFreeFormDefBox*** FFDBox,
                   unsigned short val_iZone);


};

/*!
 * \class CAdjFluidIteration
 * \brief Class for driving an iteration of the adjoint fluid system.
 * \author T. Economon
 * \version 5.0.0 "Raven"
 */
class CAdjFluidIteration : public CIteration {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CAdjFluidIteration(CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAdjFluidIteration(void);
  
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
   * \brief Perform a single iteration of the adjoint fluid system.
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
   * \brief Updates the containers for the adjoint fluid system.
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
   * \brief Monitors the convergence and other metrics for the adjoint fluid system.
   */
  void Monitor();
  
  /*!
   * \brief Outputs desired files and quantities for the adjoint fluid system.
   */
  void Output();
  
  /*!
   * \brief Postprocess ???.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   */
  void Postprocess(COutput *output,
                   CIntegration ***integration_container,
                   CGeometry ***geometry_container,
                   CSolver ****solver_container,
                   CNumerics *****numerics_container,
                   CConfig **config_container,
                   CSurfaceMovement **surface_movement,
                   CVolumetricMovement **grid_movement,
                   CFreeFormDefBox*** FFDBox,
                   unsigned short val_iZone);

  
};

/*!
 * \class CDiscAdjFluidIteration
 * \brief Class for driving an iteration of the discrete adjoint fluid system.
 * \author T. Economon
 * \version 5.0.0 "Raven"
 */
class CDiscAdjFluidIteration : public CIteration {

private:

  CFluidIteration* meanflow_iteration; /*!< \brief Pointer to the mean flow iteration class. */
  unsigned short CurrentRecording; /*!< \brief Stores the current status of the recording. */
  bool turbulent;       /*!< \brief Stores the turbulent flag. */

public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjFluidIteration(CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFluidIteration(void);
  
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
   * \brief Perform a single iteration of the adjoint fluid system.
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
   * \brief Updates the containers for the discrete adjoint fluid system.
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
   * \brief Monitors the convergence and other metrics for the discrete adjoint fluid system.
   * \param[in] ??? - Description here.
   */
  void Monitor();
  
  /*!
   * \brief Outputs desired files and quantities for the discrete adjoint fluid system.
   * \param[in] ??? - Description here.
   */
  void Output();
  
  /*!
   * \brief Postprocess ???.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   */
  void Postprocess(COutput *output,
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
   * \brief Registers all input variables of the fluid iteration.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   */
  void InitializeAdjoint(CSolver ****solver_container,
                         CGeometry*** geometry_container,
                         CConfig** config_container,
                         unsigned short iZone);

  /*!
   * \brief Registers all output variables of the fluid iteration.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] kind_recording - Kind of recording, either FLOW_CONS_VARS or MESH_COORDS
   */
  void RegisterInput(CSolver ****solver_container,
                     CGeometry*** geometry_container,
                     CConfig** config_container,
                     unsigned short iZone,
                     unsigned short kind_recording);

  /*!
   * \brief Initializes the adjoints of the output variables of the fluid iteration.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   */
  void RegisterOutput(CSolver ****solver_container,
                      CGeometry*** geometry_container,
                      CConfig** config_container,
											                      COutput* output,
                      unsigned short iZone);

  /*!
   * \brief Initializes the adjoints of the output variables of the meanflow iteration - without the contribution of the objective function
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   */
  void InitializeAdjoint_CrossTerm(CSolver ****solver_container, CGeometry ***geometry_container, CConfig **config_container, unsigned short iZone);

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

  /*!
   * \brief Compute necessary variables that depend on the conservative variables or the mesh node positions
   * (e.g. turbulent variables, normals, volumes).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] kind_recording - The kind of recording (geometry or flow).
   */
  void SetDependencies(CSolver ****solver_container,
                       CGeometry ***geometry_container,
                       CConfig **config_container,
                       unsigned short iZone,
                       unsigned short kind_recording);

  /*!
   * \brief load unsteady solution for unsteady problems
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_DirectIter - Direct iteration to load.
   */
  void LoadUnsteady_Solution(CGeometry ***geometry_container,
                      CSolver ****solver_container,
                      CConfig **config_container,
                      unsigned short val_iZone,
                      int val_DirectIter);



};

/*!
 * \class CDiscAdjFEAIteration
 * \brief Class for driving an iteration of the discrete adjoint FEM system.
 * \author R. Sanchez
 * \version 5.0.0 "Raven"
 */
class CDiscAdjFEAIteration : public CIteration {

private:

  CFEM_StructuralAnalysis* fem_iteration; /*!< \brief Pointer to the mean flow iteration class. */
  unsigned short CurrentRecording;        /*!< \brief Stores the current status of the recording. */


public:

  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjFEAIteration(CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFEAIteration(void);

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
  void Postprocess(COutput *output,
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
   * \brief Registers all input variables of the FEM iteration.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] kind_recording - Kind of recording, either FEM_VARIABLES or MESH_COORDS
   */
  void RegisterInput(CSolver ****solver_container, CGeometry*** geometry_container, CConfig** config_container, unsigned short iZone, unsigned short kind_recording);

  /*!
   * \brief Registers all output variables of the FEM iteration.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   */
  void RegisterOutput(CSolver ****solver_container, CGeometry*** geometry_container, CConfig** config_container, unsigned short iZone);
  using CIteration::RegisterOutput;
  
  /*!
   * \brief Initializes the adjoints of the output variables of the FEM iteration.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   */
  void InitializeAdjoint(CSolver ****solver_container, CGeometry*** geometry_container, CConfig** config_container, unsigned short iZone);

  /*!
   * \brief Initializes the adjoints of the output variables of the FEM iteration - without the contribution of the objective function
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   */
  void InitializeAdjoint_CrossTerm(CSolver ****solver_container, CGeometry ***geometry_container, CConfig **config_container, unsigned short iZone);

  /*!
   * \brief Record a single iteration of the direct FEM system.
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

  /*!
   * \brief Compute necessary variables that depend on the variables in the numerics (E, Nu...)
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] kind_recording - The kind of recording (geometry or flow).
   */
  void SetDependencies(CSolver ****solver_container,
                       CGeometry ***geometry_container,
                       CNumerics *****numerics_container,
                       CConfig **config_container,
                       unsigned short iZone,
                       unsigned short kind_recording);

  /*!
   * \brief load solution for dynamic problems
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_DirectIter - Direct iteration to load.
   */
  void LoadDynamic_Solution(CGeometry ***geometry_container,
                      CSolver ****solver_container,
                      CConfig **config_container,
                      unsigned short val_iZone,
                      int val_DirectIter);

};


/*!
 * \brief Iteration function for structural analysis using the Finite Element Method.
 * \author R. Sanchez.
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
 */
void FEM_StructuralIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                  CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                  CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox);
