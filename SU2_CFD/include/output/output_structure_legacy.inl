/*!
 * \file output_structure.inl
 * \brief In-Line subroutines of the <i>output.hpp</i> file.
 * \author J. Smith
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

inline su2double COutputLegacy::GetEntropyGen(unsigned short iMarkerTP, unsigned short iSpan) { return EntropyGen[iMarkerTP][iSpan]; }

inline su2double COutputLegacy::GetFlowAngleOut(unsigned short iMarkerTP, unsigned short iSpan) { return FlowAngleOut[iMarkerTP][iSpan]*180.0/PI_NUMBER; }

inline su2double COutputLegacy::GetMassFlowIn(unsigned short iMarkerTP, unsigned short iSpan) { return MassFlowIn[iMarkerTP][iSpan]; }

inline bool COutputLegacy::PrintOutput(unsigned long iIter, unsigned long iFreq) { return (iIter % iFreq == 0); }
