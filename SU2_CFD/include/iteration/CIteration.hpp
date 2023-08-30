/*!
 * \file CIteration.hpp
 * \brief Headers of the iteration classes used by SU2_CFD.
 *        Each CIteration class represents an available physics package.
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/grid_movement/CSurfaceMovement.hpp"
#include "../../../Common/include/grid_movement/CVolumetricMovement.hpp"
#include "../../../Common/include/grid_movement/CFreeFormDefBox.hpp"
#include "../../../Common/include/parallelization/mpi_structure.hpp"
#include "../integration/CIntegration.hpp"

using namespace std;

class COutput;

/*!
 * \class CIteration
 * \ingroup Drivers
 * \brief Parent class for defining a single iteration of a physics problem.
 * \author T. Economon
 */
class CIteration {
 protected:
  const int rank,             /*!< \brief MPI Rank. */
      size;                   /*!< \brief MPI Size. */
  const unsigned short nZone; /*!< \brief Total number of zones in the problem. */
  const unsigned short nInst; /*!< \brief Total number of instances in the problem. */

  const bool multizone, /*!< \brief Flag for multizone problems. */
      singlezone;       /*!< \brief Flag for singlezone problems. */

  su2double StartTime{0.0}, /*!< \brief Tracking wall time. */
      StopTime{0.0}, UsedTime{0.0};

 public:
  /*!
   * \brief Constructor of the class.
   */
  explicit CIteration(const CConfig* config)
      : rank(SU2_MPI::GetRank()),
        size(SU2_MPI::GetSize()),
        nZone(config->GetnZone()),
        nInst(config->GetnTimeInstances()),
        multizone(config->GetMultizone_Problem()),
        singlezone(!multizone) {}

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CIteration(void) = default;

  /*!
   * \brief Updates the positions and grid velocities for dynamic meshes between physical time steps.
   * \author T. Economon
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] IntIter - Current psuedo time iteration number.
   * \param[in] TimeIter - Current physical time iteration number.
   */
  virtual void SetGrid_Movement(CGeometry** geometry, CSurfaceMovement* surface_movement,
                                CVolumetricMovement* grid_movement, CSolver*** solver, CConfig* config,
                                unsigned long IntIter, unsigned long TimeIter);
  /*!
   * \brief Run the mesh deformation algorithms.
   * \author R. Sanchez
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] kind_recording - Current kind of recording.
   */
  void SetMesh_Deformation(CGeometry** geometry, CSolver** solver_container, CNumerics*** numerics_container,
                           CConfig* config_container, RECORDING kind_recording);

  /*!
   * \brief A virtual member.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instantiation.
   */
  virtual void Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                          CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                          CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                          unsigned short val_iInst){}

  /*!
   * \brief A virtual member.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instantiation.
   */
  virtual void Iterate(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                       CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                       CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                       unsigned short val_iInst){}

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instantiation.
   * \param[in] CrossTerm - Boolean for CrossTerm.
   */
  virtual void IterateDiscAdj(CGeometry**** geometry, CSolver***** solver, CConfig** config,
                              unsigned short val_iZone, unsigned short val_iInst, bool CrossTerm) {}

  /*!
   * \brief A virtual member.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instantiation.
   */
  virtual void Solve(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                     CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                     CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                     unsigned short val_iInst){}

  /*!
   * \brief A virtual member.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instantiation.
   */
  virtual void Update(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                      CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                      CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                      unsigned short val_iInst){}

  /*!
   * \brief A virtual member.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instantiation.
   */
  virtual void Predictor(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                         CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                         CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                         unsigned short val_iInst){}

  /*!
   * \brief A virtual member.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instantiation.
   */
  virtual void Relaxation(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                          CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                          CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                          unsigned short val_iInst){}

  /*!
   * \brief A virtual member.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instantiation.
   */
  virtual bool Monitor(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                       CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                       CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                       unsigned short val_iInst) {
    return false;
  }

  /*!
   * \brief A virtual member.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] InnerIter - Current psuedo time iteration number.
   * \param[in] StopCalc - Boolean storing whether to stop the simulation
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instantiation.
   */
  void Output(COutput* output, CGeometry**** geometry, CSolver***** solver, CConfig** config, unsigned long InnerIter,
              bool StopCalc, unsigned short val_iZone, unsigned short val_iInst);

  /*!
   * \brief A virtual member.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instantiation.
   */
  virtual void Postprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                           CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                           CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                           unsigned short val_iInst) {}

  virtual void InitializeAdjoint(CSolver***** solver, CGeometry**** geometry, CConfig** config, unsigned short iZone,
                                 unsigned short iInst) {}

  virtual void RegisterInput(CSolver***** solver, CGeometry**** geometry, CConfig** config, unsigned short iZone,
                             unsigned short iInst, RECORDING kind_recording) {}

  virtual void SetDependencies(CSolver***** solver, CGeometry**** geometry, CNumerics****** numerics, CConfig** config,
                               unsigned short iZone, unsigned short iInst, RECORDING kind_recording) {}

  virtual void RegisterOutput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                              unsigned short iZone, unsigned short iInst) {}
};
