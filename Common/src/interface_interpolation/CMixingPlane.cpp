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

    targetSpans.resize(config[donorZone]->GetnMarker_MixingPlaneInterface());

    /*--- On the donor side ---*/
    for (auto iMarkerInt = 0u; iMarkerInt < nMarkerInt; iMarkerInt++){
        int markDonor = -1, markTarget = -1;
        unsigned short donorFlag = 0, targetFlag = 0;

        markDonor = donor_config->FindMixingPlaneInterfaceMarker(iMarkerInt, donor_config->GetnMarker_All());
        donorFlag = (markDonor != -1) ? donor_config->GetMarker_All_MixingPlaneInterface(markDonor) : -1;

#ifdef HAVE_MPI
    auto buffMarkerDonor = new int[size];
    auto buffDonorFlag = new int[size];
    for (int iSize=0; iSize<size; iSize++){
        buffMarkerDonor[iSize] = -1;
        buffDonorFlag[iSize] = -1;
    }

    SU2_MPI::Allgather(&markDonor, 1 , MPI_INT, buffMarkerDonor, 1, MPI_INT, SU2_MPI::GetComm());
    SU2_MPI::Allgather(&donorFlag, 1 , MPI_INT, buffDonorFlag, 1, MPI_INT, SU2_MPI::GetComm());

    markDonor= -1;
    donorFlag= -1;

    for (int iSize=0; iSize<size; iSize++) {
        if(buffMarkerDonor[iSize] >= 0.0) {
            markDonor = buffMarkerDonor[iSize];
            donorFlag = buffDonorFlag[iSize];
            break;
        }
    }
    delete [] buffMarkerDonor;
    delete [] buffDonorFlag;
#endif

        markTarget = target_config->FindMixingPlaneInterfaceMarker(iMarkerInt, target_config->GetnMarker_All());
        targetFlag = (markTarget != -1) ? target_config->GetMarker_All_MixingPlaneInterface(markTarget) : -1;

        if (markTarget == -1 || markDonor == -1) continue;

        nSpanDonor = donor_config->GetnSpanWiseSections();
        nSpanTarget = target_config->GetnSpanWiseSections();

        targetSpans[iMarkerInt].resize(nSpanTarget+1);

        const auto spanValuesDonor = donor_geometry->GetSpanWiseValue(donorFlag);
        const auto spanValuesTarget = target_geometry->GetSpanWiseValue(targetFlag);

        /*--- Interpolation at hub, shroud & 1D values ---*/
        targetSpans[iMarkerInt][0].donorSpan = 0;
        targetSpans[iMarkerInt][0].coefficient = 0.0;
        if (nDim > 2) {
            targetSpans[iMarkerInt][nSpanTarget-2].donorSpan = nSpanDonor-2;
            targetSpans[iMarkerInt][nSpanTarget-2].coefficient = 0.0;
            targetSpans[iMarkerInt][nSpanTarget-1].donorSpan = nSpanDonor-1;
            targetSpans[iMarkerInt][nSpanTarget-1].coefficient = 0.0;
        }

        for(auto iSpanTarget = 1; iSpanTarget < nSpanTarget - 2; iSpanTarget++){
            auto &targetSpan = targetSpans[iMarkerInt][iSpanTarget];

            auto tSpan = 0; // Nearest donor span index
            auto kSpan = 0; // Lower bound donor span for interpolation
            su2double coeff = 0.0; // Interpolation coefficient
            su2double minDist = 10E+06;
            
            switch(donor_config->GetKind_MixingPlaneInterface()){
                case MATCHING:
                    targetSpan.donorSpan = iSpanTarget;
                    targetSpan.coefficient = 0.0;
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
                    targetSpan.donorSpan = tSpan;
                    targetSpan.coefficient = 0.0;
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
                    targetSpan.donorSpan = kSpan;
                    targetSpan.coefficient = coeff;
                    break;
                    
                default:
                    SU2_MPI::Error("MixingPlane interface option not implemented yet", CURRENT_FUNCTION);
                    break;
            }
        }
    }
}