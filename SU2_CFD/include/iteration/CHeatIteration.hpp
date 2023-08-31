/*!
 * \file CHeatIteration.hpp
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

#include "CFluidIteration.hpp"

/*!
 * \class CHeatIteration
 * \ingroup Drivers
 * \brief Class for driving an iteration of the heat system.
 * \author T. Economon
 */
class CHeatIteration : public CFluidIteration {
 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  explicit CHeatIteration(const CConfig* config) : CFluidIteration(config) {}

  /*!
   * \brief Perform a single iteration of the heat system.
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
   * \brief Updates the containers for the heat system.
   */
  void Update(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
              CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
              CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
              unsigned short val_iInst) override;

  /*!
   * \brief Override the preprocessing of CFluidIteration.
   */
  inline void Preprocess(COutput*, CIntegration****, CGeometry****, CSolver*****, CNumerics******, CConfig**,
                         CSurfaceMovement**, CVolumetricMovement***, CFreeFormDefBox***, unsigned short,
                         unsigned short) override {}

  /*!
   * \brief Override the postprocessing of CFluidIteration.
   */
  inline void Postprocess(COutput*, CIntegration****, CGeometry****, CSolver*****, CNumerics******, CConfig**,
                          CSurfaceMovement**, CVolumetricMovement***, CFreeFormDefBox***, unsigned short,
                          unsigned short) override {}
};
