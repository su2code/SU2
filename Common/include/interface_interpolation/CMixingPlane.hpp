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
#include "../option_structure.hpp"

/*!
 * \brief Mixing plane interpolation.
 * \note This contains several interpolation methods used in the mixing plane interpolation
 * and enables the mixing state class structure for proper recording in AD mode
 * \ingroup Interfaces
 */
class CMixingPlane final : public CInterpolator {
 public:
  CMixingPlane(CGeometry**** geometry_container, const CConfig* const* config, unsigned int iZone, unsigned int jZone);

  void SetTransferCoeff(CGeometry**** geometry, const CConfig* const* config) override;

  inline CSpanDonorInfo MapMatchingSpan(unsigned short iSpanTarget) { return {iSpanTarget, 0.0}; }

  inline CSpanDonorInfo MapNearestSpan(const su2double iSpanTargetValue, const su2double* spanValuesDonor, unsigned long nSpanDonor) {
    unsigned short tSpan = 0;         // Nearest donor span index
    auto minDist = std::numeric_limits<su2double>::max();

    for (auto iSpanDonor = 0u; iSpanDonor < nSpanDonor - 1; iSpanDonor++) {
      const auto dist = abs(iSpanTargetValue - spanValuesDonor[iSpanDonor]);
      if (dist < minDist) {
        minDist = dist;
        tSpan = iSpanDonor;
      }
    }
    return {tSpan, 0.0};
  };

  inline CSpanDonorInfo MapLinearInterpolationSpan(const su2double iSpanTargetValue, const su2double* spanValuesDonor, unsigned long nSpanDonor, int rank) {
    unsigned short kSpan = 0;         // Lower bound donor span for interpolation
    auto minDist = std::numeric_limits<su2double>::max();
    su2double coeff = 0.0;  // Interpolation coefficient

    if (iSpanTargetValue < spanValuesDonor[0]) {
      PrintClampingWarning(rank, true);
      return {0, 0.0};
    } 
    
    if (iSpanTargetValue > spanValuesDonor[nSpanDonor - 1]) {
      PrintClampingWarning(rank, false);
      return {nSpanDonor - 1, 0.0};
    }

    for (auto iSpanDonor = 0u; iSpanDonor < nSpanDonor - 1; iSpanDonor++) {
      const auto dist = abs(iSpanTargetValue - spanValuesDonor[iSpanDonor]);
      if (dist < minDist && iSpanTargetValue > spanValuesDonor[iSpanDonor]) {
        kSpan = iSpanDonor;
        minDist = dist;
        break;
      }
    }
    coeff = (iSpanTargetValue - spanValuesDonor[kSpan]) /
                (spanValuesDonor[kSpan + 1] - spanValuesDonor[kSpan]);
    return {kSpan, coeff};
  };
  
  inline void PrintClampingWarning(int rank, bool atHub) {
    if (rank != MASTER_NODE) return;
    cout << "Warning! Target spans exist outside the bounds of donor spans! Clamping interpolator..." << endl;
    cout << (atHub ? "This is an issue at the hub." : "This is an issue at the shroud.") << endl;
    cout << "Setting coeff = 0.0 and transferring endwall value!" << endl;
  };

  /*!
   * \brief Write interpolation details to file.
   * \param[in] filename - Name of output file.
   * \param[in] config - Configuration for all zones.
   */
  void WriteInterpolationDetails(const string& filename, const CConfig* const* config) override;
};
