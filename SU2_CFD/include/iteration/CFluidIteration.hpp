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

#include "CIteration.hpp"

/*!
 * \class CFluidIteration
 * \ingroup Drivers
 * \brief Class for driving an iteration of the fluid system.
 * \author T. Economon
 */
class CFluidIteration : public CIteration {
 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  explicit CFluidIteration(const CConfig* config) : CIteration(config) {}

  /*!
   * \brief Preprocessing to prepare for an iteration of the physics.
   * \param[in] ??? - Description here.
   */
  void Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                  CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                  CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                  unsigned short val_iInst) override;

  /*!
   * \brief Perform a single iteration of the fluid system.
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
   * \param[in] val_iInst - Index of the instance layer.
   */
  void Iterate(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
               CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
               CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
               unsigned short val_iInst) override;

  /*!
   * \brief Iterate the fluid system for a number of Inner_Iter iterations.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   */
  void Solve(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
             CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
             CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
             unsigned short val_iInst) override;

  /*!
   * \brief Updates the containers for the fluid system.
   * \param[in] ??? - Description here.
   */
  void Update(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
              CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
              CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
              unsigned short val_iInst) override;

  /*!
   * \brief Monitors the convergence and other metrics for the fluid system.
   * \param[in] ??? - Description here.
   */
  bool Monitor(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
               CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
               CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
               unsigned short val_iInst) override;

  /*!
   * \brief Postprocesses the fluid system before heading to another physics system or the next iteration.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Postprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                   CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                   CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                   unsigned short val_iInst) override;

 private:

  /*!
   * \brief Imposes a gust via the grid velocities.
   * \author S. Padron
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   */
  void SetWind_GustField(CConfig* config, CGeometry** geometry, CSolver*** solver);

  /*!
   * \brief Reads and initializes the vortex positions, strengths and gradient.
   * \author S. Padron
   * \param[in] nVortex - number of vortices.
   * \param[in] x0 - Vector of x-loc of the vortices.
   * \param[in] y0 - Vector of y-loc of the vortices.
   * \param[in] vort_strength - Vector of vortex strengths.
   * \param[in] r_core - Vector of vortex core size.
   */
  void InitializeVortexDistribution(unsigned long& nVortex, vector<su2double>& x0, vector<su2double>& y0,
                                    vector<su2double>& vort_strength, vector<su2double>& r_core);

  /*!
   * \brief Fixed CL monitoring function
   * \author J. Mukhopadhaya
   * \param[in] output - Pointer to the COutput class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Pointer to the flow solver
   * \param[in] config - Definition of the particular problem.
   * \return Boolean indicating weather calculation should be stopped
   */
  bool MonitorFixed_CL(COutput* output, CGeometry* geometry, CSolver** solver, CConfig* config);

  /*!
   * \brief Store old aeroelastic solutions
   * \param[in,out] config - Definition of the particular problem.
   */
  void SetDualTime_Aeroelastic(CConfig* config) const;

};
