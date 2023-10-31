/*!
 * \file CMixingPlaneInterface.cpp
 * \brief Declaration and inlines of the class to transfer average variables
 *        needed for MixingPlane computation from a generic zone into another one.
 * \author S. Vitale
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

#include "../../../include/interfaces/cfd/CMixingPlaneInterface.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../include/solvers/CSolver.hpp"

CMixingPlaneInterface::CMixingPlaneInterface(unsigned short val_nVar, unsigned short val_nConst){
  nVar = val_nVar;
  Donor_Variable     = new su2double[nVar + 5]();
  Target_Variable    = new su2double[nVar + 5]();
}

void CMixingPlaneInterface::SetSpanWiseLevels(const CConfig *donor_config, const CConfig *target_config){

  unsigned short iSpan;
  nSpanMaxAllZones = donor_config->GetnSpanMaxAllZones();


  SpanValueCoeffTarget = new su2double[target_config->GetnSpanWiseSections() + 1];
  SpanLevelDonor       = new unsigned short[target_config->GetnSpanWiseSections() + 1];


  for (iSpan = 0; iSpan < target_config->GetnSpanWiseSections() + 1;iSpan++){
    SpanValueCoeffTarget[iSpan] = 0.0;
    SpanLevelDonor[iSpan]       = 1;
  }
}

void CMixingPlaneInterface::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry,
                                              const CConfig *donor_config, unsigned long Marker_Donor,
                                              unsigned long iSpan, unsigned long rank) {

  unsigned short nDim = nVar - 2;
  bool turbulent = (donor_config->GetKind_Turb_Model() != TURB_MODEL::NONE);



  Donor_Variable[0] = donor_solution->GetAverageDensity(Marker_Donor, iSpan);
  Donor_Variable[1] = donor_solution->GetAveragePressure(Marker_Donor, iSpan);
  Donor_Variable[2] = donor_solution->GetAverageTurboVelocity(Marker_Donor, iSpan)[0];
  Donor_Variable[3] = donor_solution->GetAverageTurboVelocity(Marker_Donor, iSpan)[1];

  if(nDim == 3){
    Donor_Variable[4] = donor_solution->GetAverageTurboVelocity(Marker_Donor, iSpan)[2];
  }
  else{
    Donor_Variable[4] = -1.0;
  }

  if(turbulent){
    Donor_Variable[5] = donor_solution->GetAverageNu(Marker_Donor, iSpan);
    Donor_Variable[6] = donor_solution->GetAverageKine(Marker_Donor, iSpan);
    Donor_Variable[7] = donor_solution->GetAverageOmega(Marker_Donor, iSpan);
  }
  else{
    Donor_Variable[5] = -1.0;
    Donor_Variable[6] = -1.0;
    Donor_Variable[7] = -1.0;
  }

}


void CMixingPlaneInterface::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                                               const CConfig *target_config, unsigned long Marker_Target,
                                               unsigned long iSpan, unsigned long rank) {

  unsigned short nDim = nVar - 2;
  bool turbulent = (target_config->GetKind_Turb_Model() != TURB_MODEL::NONE);


  target_solution->SetExtAverageDensity(Marker_Target, iSpan, Target_Variable[0]);
  target_solution->SetExtAveragePressure(Marker_Target, iSpan, Target_Variable[1]);
  target_solution->SetExtAverageTurboVelocity(Marker_Target, iSpan, 0, Target_Variable[2]);
  target_solution->SetExtAverageTurboVelocity(Marker_Target, iSpan, 1, Target_Variable[3]);

  if(nDim == 3){
    target_solution->SetExtAverageTurboVelocity(Marker_Target, iSpan, 2, Target_Variable[4]);
  }

  if(turbulent){
    target_solution->SetExtAverageNu(Marker_Target, iSpan, Target_Variable[5]);
    target_solution->SetExtAverageKine(Marker_Target, iSpan, Target_Variable[6]);
    target_solution->SetExtAverageOmega(Marker_Target, iSpan,  Target_Variable[7]);
  }

}

void CMixingPlaneInterface::SetAverageValues(CSolver *donor_solution, CSolver *target_solution,
                                             unsigned short donorZone){
  unsigned short iSpan;

  for(iSpan = 0; iSpan<nSpanMaxAllZones +1; iSpan++){
    /*--- trnasfer inviscid quantities ---*/
    target_solution->SetDensityIn(donor_solution->GetDensityIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetPressureIn(donor_solution->GetPressureIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetTurboVelocityIn(donor_solution->GetTurboVelocityIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetDensityOut(donor_solution->GetDensityOut(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetPressureOut(donor_solution->GetPressureOut(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetTurboVelocityOut(donor_solution->GetTurboVelocityOut(donorZone, iSpan), donorZone, iSpan);

    /*--- transfer turbulent quantities ---*/
    target_solution->SetKineIn(donor_solution->GetKineIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetOmegaIn(donor_solution->GetOmegaIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetNuIn(donor_solution->GetNuIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetKineOut(donor_solution->GetKineOut(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetOmegaOut(donor_solution->GetOmegaOut(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetNuOut(donor_solution->GetNuOut(donorZone, iSpan), donorZone, iSpan);

  }
}

void CMixingPlaneInterface::SetAverageTurboGeoValues(CGeometry *donor_geometry, CGeometry *target_geometry,
                                                     unsigned short donorZone){
  unsigned short iSpan;

  for(iSpan = 0; iSpan<nSpanMaxAllZones+1; iSpan++){
    target_geometry->SetTurboRadiusIn(donor_geometry->GetTurboRadiusIn(donorZone, iSpan), donorZone, iSpan);
    target_geometry->SetSpanAreaIn(donor_geometry->GetSpanAreaIn(donorZone, iSpan), donorZone, iSpan);
    target_geometry->SetTangGridVelIn(donor_geometry->GetTangGridVelIn(donorZone, iSpan), donorZone, iSpan);
    target_geometry->SetTurboRadiusOut(donor_geometry->GetTurboRadiusOut(donorZone, iSpan), donorZone, iSpan);
    target_geometry->SetSpanAreaOut(donor_geometry->GetSpanAreaOut(donorZone, iSpan), donorZone, iSpan);
    target_geometry->SetTangGridVelOut(donor_geometry->GetTangGridVelOut(donorZone, iSpan), donorZone, iSpan);
  }

}
