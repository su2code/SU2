/*!
 * \file CTurboIteration.hpp
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
 * \class CTurboIteration
 * \ingroup Drivers
 * \brief Class for driving an iteration for turbomachinery simulation.
 * \author T. Economon
 */
class CTurboIteration : public CFluidIteration {
 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  explicit CTurboIteration(const CConfig *config) : CFluidIteration(config) {}

  /*!
   * \brief Preprocessing to prepare for an iteration of the physics.
   * \param[in] ??? - Description here.
   */
  void Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
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
};
