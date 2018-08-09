/*!
 * \file output_structure.inl
 * \brief In-Line subroutines of the <i>output_structure.hpp</i> file.
 * \author J. Smith
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

inline su2double COutput::GetEntropyGen(unsigned short iMarkerTP, unsigned short iSpan) { return EntropyGen[iMarkerTP][iSpan]; }

inline su2double COutput::GetEntropyGenAvg_HB() { return EntropyGenAverage_HB; }

inline su2double COutput::GetPower_HB() { return Power_HB; }

inline su2double COutput::GetTotalWorkDone_HB() { return TotalWorkDone_Surface_HB; }

inline su2double COutput::GetFlowAngleOut(unsigned short iMarkerTP, unsigned short iSpan) { return FlowAngleOut[iMarkerTP][iSpan]*180.0/PI_NUMBER; }

inline su2double COutput::GetMassFlowIn(unsigned short iMarkerTP, unsigned short iSpan) { return MassFlowIn[iMarkerTP][iSpan]; }

inline su2double COutput::GetTotalPressureLoss(unsigned short iMarkerTP, unsigned short iSpan) { return TotalPressureLoss[iMarkerTP][iSpan]; }

inline su2double COutput::GetKineticEnergyLoss(unsigned short iMarkerTP, unsigned short iSpan) { return KineticEnergyLoss[iMarkerTP][iSpan]; }
