/*!
 * \file CDiscAdjFlowTractionInterface.hpp
 * \brief Declaration and inlines of the class to transfer flow tractions
 *        from a fluid zone into a structural zone in a discrete adjoint simulation.
 * \author Ruben Sanchez
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

#include "CFlowTractionInterface.hpp"

/*!
 * \brief Fluid-structure interface (forces or tractions) for discrete adjoint problems.
 * \note This is a small specialization of CFlowTractionInterface that ommits load ramping.
 * \ingroup Interfaces
 */
class CDiscAdjFlowTractionInterface : public CFlowTractionInterface {
public:
  /*!
   * \overload
   * \param[in] val_nVar - Number of variables that need to be transferred.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjFlowTractionInterface(unsigned short val_nVar, unsigned short val_nConst,
                                const CConfig *config, bool conservative_);

  /*!
   * \brief Retrieve some constants needed for the calculations.
   * \param[in] donor_solution - Solution from the donor mesh.
   * \param[in] target_solution - Solution from the target mesh.
   * \param[in] donor_geometry - Geometry of the donor mesh.
   * \param[in] target_geometry - Geometry of the target mesh.
   * \param[in] donor_config - Definition of the problem at the donor mesh.
   * \param[in] target_config - Definition of the problem at the target mesh.
   */
  void GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                             CGeometry *donor_geometry, CGeometry *target_geometry,
                             const CConfig *donor_config, const CConfig *target_config) override;

};
