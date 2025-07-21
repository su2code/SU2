/*!
 * \file CFEAIteration.hpp
 * \brief Headers of the iteration classes used by SU2_CFD.
 *        Each CIteration class represents an available physics package.
 * \author F. Palacios, T. Economon
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "CIteration.hpp"

class CFEAIteration;

/*!
 * \class CDiscAdjFEAIteration
 * \ingroup DiscAdj
 * \brief Class for driving an iteration of the discrete adjoint FEM system.
 * \author R. Sanchez
 */
class CDiscAdjFEAIteration final : public CIteration {
 private:
  unsigned short CurrentRecording; /*!< \brief Stores the current status of the recording. */

  /*!
   * \brief load solution for dynamic problems
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instance.
   * \param[in] val_DirectIter - Direct iteration to load.
   */
  void LoadDynamic_Solution(CGeometry**** geometry, CSolver***** solver, CConfig** config, unsigned short val_iZone,
                            unsigned short val_iInst, int val_DirectIter);

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  explicit CDiscAdjFEAIteration(const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFEAIteration(void) override;

  /*!
   * \brief Preprocessing to prepare for an iteration of the physics.
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
   * \param[in] val_iInst - Index of the instance.
   */
  void Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                  CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                  CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                  unsigned short val_iInst) override;

  /*!
   * \brief Perform a single iteration of the adjoint FEA problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instance.
   * \param[in] CrossTerm - Boolean for CrossTerm.
   */
  void IterateDiscAdj(CGeometry**** geometry, CSolver***** solver, CConfig** config,
                      unsigned short val_iZone, unsigned short val_iInst, bool CrossTerm) override;

  /*!
   * \brief Monitors the convergence and other metrics for the discrete adjoint FEA problem.
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
   * \param[in] val_iInst - Index of the instance.
   */
  bool Monitor(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
               CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
               CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
               unsigned short val_iInst) override;

  /*!
   * \brief Postprocesses the discrete adjoint FEA problem before heading to another physics system or the next
   * iteration.
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
   * \param[in] val_iInst - Index of the instance.
   */
  void Postprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                   CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                   CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                   unsigned short val_iInst) override;

  /*!
   * \brief Registers all input variables of the FEM iteration.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] kind_recording - Kind of recording, either FEM_VARIABLES or MESH_COORDS
   */
  void RegisterInput(CSolver***** solver, CGeometry**** geometry, CConfig** config, unsigned short iZone,
                     unsigned short iInst, RECORDING kind_recording) override;

  /*!
   * \brief Registers all output variables of the FEM iteration.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] iInst - Index of the zone.
   */
  void RegisterOutput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                      unsigned short iZone, unsigned short iInst) override;

  /*!
   * \brief Initializes the adjoints of the output variables of the FEM iteration.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] iInst - Index of the zone.
   */
  void InitializeAdjoint(CSolver***** solver, CGeometry**** geometry, CConfig** config, unsigned short iZone,
                         unsigned short iInst) override;

  /*!
   * \brief Compute necessary variables that depend on the variables in the numerics (E, Nu...)
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] iInst - Index of the zone.
   * \param[in] kind_recording - The kind of recording (geometry or flow).
   */
  void SetDependencies(CSolver***** solver, CGeometry**** geometry, CNumerics****** numerics, CConfig** config,
                       unsigned short iZone, unsigned short iInst, RECORDING kind_recording) override;

};
