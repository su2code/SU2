/*!
 * \file mpi_structure.cpp
 * \brief Main subroutines for the mpi structures.
 * \author T. Albring
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

#include "../include/mpi_structure.hpp"

int CBaseMPIWrapper::Rank = 0;
int CBaseMPIWrapper::Size = 1;
CBaseMPIWrapper::Comm CBaseMPIWrapper::currentComm = MPI_COMM_WORLD;

#ifdef HAVE_MPI
int  CBaseMPIWrapper::MinRankError;
bool CBaseMPIWrapper::winMinRankErrorInUse = false;
CBaseMPIWrapper::Win CBaseMPIWrapper::winMinRankError;
#endif

#ifdef HAVE_MPI
#if defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE
//AMPI_ADOUBLE_TYPE* AMPI_ADOUBLE;
#include <medi/medi.cpp>
#endif // defined CODI_REVERSE_TYPE || defined CODI_FORWARD_TYPE

#endif// HAVE_MPI
