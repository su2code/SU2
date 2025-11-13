/*!
 * \file CMixingPlane.cpp
 * \brief Implementation of mixing plane interpolation methods.
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

#include "../../include/interface_interpolation/CMixingPlane.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"

CMixingPlane::CMixingPlane(CGeometry**** geometry_container, const CConfig* const* config, unsigned int iZone,
                                   unsigned int jZone)
    : CInterpolator(geometry_container, config, iZone, jZone) {
  SetTransferCoeff(config);
}

void CMixingPlane::SetTransferCoeff(const CConfig* const* config) {
    const auto nMarkerInt = config[donorZone]->GetnMarker_MixingPlaneInterface() / 2;
    const auto nDim = donor_geometry->GetnDim();

    const auto donor_config = config[donorZone];
    const auto target_config = config[targetZone];

    const auto nMarkerDonor   = donor_geometry->GetnMarker(); // Number of markers on the interfacce
    const auto nMarkerTarget  = target_geometry->GetnMarker();
    //TODO turbo this approach only works if all the turboamchinery marker
    //    of all zones have the same amount of span wise sections.
    //TODO turbo initialization needed for the MPI routine should be place somewhere else.
    auto nSpanDonor     = donor_config->GetnSpanWiseSections();
    auto nSpanTarget    = target_config->GetnSpanWiseSections();

    targetSpans.resize(nMarkerInt);

    /*--- On the donor side ---*/
    for (auto iMarkerInt = 0u; iMarkerInt < nMarkerInt; iMarkerInt++){
        const auto markDonor = donor_config->FindInterfaceMarker(iMarkerInt);
        const auto markTarget = target_config->FindInterfaceMarker(iMarkerInt);

        // Spans are defined on one processor only, check to see if other processors have interface
#ifdef HAVE_MPI
    auto buffMarkDonor = new int[size];
    auto buffMarkTarget = new int[size];
    for (int iSize = 0; iSize<size; iSize++){
        buffMarkDonor[iSize] = -1;
        buffMarkTarget[iSize] = -1;
    }

    SU2_MPI::Allgather(&markDonor, 1, MPI_INT, buffMarkDonor, 1, MPI_INT, SU2_MPI::GetComm());
    SU2_MPI::Allgather(&markTarget, 1, MPI_INT, buffMarkTarget, 1, MPI_INT, SU2_MPI::GetComm());

    delete [] buffMarkDonor;
    delete [] buffMarkTarget;
#endif

        if (!CheckInterfaceBoundary(markDonor, markTarget)) continue;

        if (markDonor != -1) nSpanDonor = donor_config->GetnSpanWiseSections();
        if (markTarget != -1) nSpanTarget = target_config->GetnSpanWiseSections();

        if (nSpanTarget) targetSpans[markTarget].resize(nSpanTarget);

        const auto spanValuesDonor = donor_geometry->GetSpanWiseValue(markDonor);
        const auto spanValuesTarget = target_geometry->GetSpanWiseValue(markTarget);

        /*--- Interpolation of values at hub & shroud ---*/
        targetSpans[iMarkerInt].globalSpan[0] = 0;
        targetSpans[iMarkerInt].coefficient[0] = 0.0;
        targetSpans[iMarkerInt].globalSpan[nSpanTarget] = nSpanTarget;
        targetSpans[iMarkerInt].coefficient[nSpanTarget] = 0.0;

        for(auto iSpanTarget = 1; iSpanTarget < nSpanTarget - 2; iSpanTarget++){
            auto tSpan = 0; // Nearest donor span index
            auto kSpan = 0; // Lower bound donor span for interpolation
            su2double coeff = 0.0; // Interpolation coefficient
            su2double minDist = 10E+06;
            
            switch(donor_config->GetKind_MixingPlaneInterface()){
                case MATCHING:
                    targetSpans[iMarkerInt].globalSpan[iSpanTarget] = iSpanTarget;
                    targetSpans[iMarkerInt].coefficient[iSpanTarget] = 0.0;
                    break;
                    
                case NEAREST_SPAN:
                    // Find the nearest donor span
                    for (auto iSpanDonor = 0; iSpanDonor < nSpanDonor - 1; iSpanDonor++) {
                        const auto dist = abs(spanValuesTarget[iSpanTarget] - spanValuesDonor[iSpanDonor]);
                        if(dist < minDist){
                            minDist = dist;
                            tSpan = iSpanDonor;
                        }
                    }
                    targetSpans[iMarkerInt].globalSpan[iSpanTarget] = tSpan;
                    targetSpans[iMarkerInt].coefficient[iSpanTarget] = 0.0;
                    break;
                    
                case LINEAR_INTERPOLATION:
                    // Find the donor span interval that brackets the target span
                    for (auto iSpanDonor = 1; iSpanDonor < nSpanDonor - 2; iSpanDonor++) {
                        if(spanValuesTarget[iSpanTarget] >= spanValuesDonor[iSpanDonor] && 
                        spanValuesTarget[iSpanTarget] <= spanValuesDonor[iSpanDonor + 1]){
                            kSpan = iSpanDonor;
                            break;
                        }
                    }
                    // Calculate interpolation coefficient
                    coeff = (spanValuesTarget[iSpanTarget] - spanValuesDonor[kSpan]) / 
                            (spanValuesDonor[kSpan + 1] - spanValuesDonor[kSpan]);
                    targetSpans[iMarkerInt].globalSpan[iSpanTarget] = kSpan;
                    targetSpans[iMarkerInt].coefficient[iSpanTarget] = coeff;
                    break;
                    
                default:
                    SU2_MPI::Error("MixingPlane interface option not implemented yet", CURRENT_FUNCTION);
                    break;
            }


            // for (auto iSpanDonor = 1; iSpanDonor < nSpanDonor - 2; iSpanDonor++) {
            //     const auto test = abs(SpanValuesTarget[iSpanTarget] - SpanValuesDonor[iSpanDonor]);
            //     const auto test2 = abs(SpanValuesTarget[iSpanTarget] - SpanValuesDonor[iSpanDonor]);
            //     if(test < dist && SpanValuesTarget[iSpanTarget] > SpanValuesDonor[iSpanDonor]){
            //         dist = test;
            //         kSpan = iSpanDonor;
            //     }
            //     if(test2 < dist2){
            //         dist2 = test2;
            //         tSpan = iSpanDonor;
            //     }
            // }
            // switch(donor_config->GetKind_MixingPlaneInterface()){
            //     case MATCHING:
            //         targetSpans.globalSpan[iSpanTarget]  = iSpanTarget;
            //         targetSpans.coefficent[iSpanTarget]  = 0.0;
            //         break;
            //     case NEAREST_SPAN:
            //         targetSpans.globalSpan[iSpanTarget]  = tSpan;
            //         targetSpans.coefficent[iSpanTarget]  = 0.0;
            //         break;
            //     case LINEAR_INTERPOLATION:
            //         targetSpans.globalSpan[iSpanTarget]  = kSpan;
            //         targetSpans.coefficent[iSpanTarget]  = (SpanValuesTarget[iSpanTarget] - SpanValuesDonor[kSpan])
            //                                         /(SpanValuesDonor[kSpan + 1] - SpanValuesDonor[kSpan]);
            //         break;
            //     default:
            //         SU2_MPI::Error("MixingPlane interface option not implemented yet", CURRENT_FUNCTION);
            //         break;
            // }
    }
}
}