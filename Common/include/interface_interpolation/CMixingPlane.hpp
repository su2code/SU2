/*!
 * \file CMixingPlane.hpp
 * \brief Header of mixing plane interpolation methods.
 * \author J. Kelly
 * \version 8.3.0 "Harrier"
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
 #include "CInterpolator.hpp"

/*!
 * \brief Mixing plane interpolation.
 * \note This contains several interpolation methods used in the mixing plane interpolation
 * and enables the mixing state class structure for proper recording in AD mode
 * \ingroup Interfaces
 */
class CMixingPlane final : public CInterpolator {
    public:
    CMixingPlane(CGeometry**** geometry_container, const CConfig* const* config, unsigned int iZone,
                                   unsigned int jZone);

    void SetTransferCoeff(const CConfig* const* config) override;

    /*!
     * \brief Write interpolation details to file.
     * \param[in] filename - Name of output file.
     * \param[in] config - Configuration for all zones.
     */
    void WriteInterpolationDetails(const string& filename, const CConfig* const* config) override;
};