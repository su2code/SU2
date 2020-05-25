/*!
 * \file codi_reverse_structure.hpp
 * \brief Header for codi reverse type definition.
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

#pragma once

#include "codi.hpp"
#include "codi/tools/dataStore.hpp"

#ifndef CODI_INDEX_TAPE
#  define CODI_INDEX_TAPE 0
#endif

#ifndef CODI_PRIMAL_TAPE
#  define CODI_PRIMAL_TAPE 0
#endif

#ifndef CODI_PRIMAL_INDEX_TAPE
#  define CODI_PRIMAL_INDEX_TAPE 0
#endif

#if CODI_INDEX_TAPE
  typedef codi::RealReverseIndex su2double;
#elif CODI_PRIMAL_TAPE
  typedef codi::RealReversePrimal su2double;
#elif CODI_PRIMAL_INDEX_TAPE
  typedef codi::RealReversePrimalIndex su2double;
#else
  typedef codi::RealReverse su2double;
#endif
