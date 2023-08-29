/*!
 * \file omp_structure.cpp
 * \brief Source file counterpart for omp_structure.hpp.
 * \note Contains OpDiLib initialization, finalization and includes the OpDiLib source file.
 * \author J. Blühdorn
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

#include "omp_structure.hpp"

void omp_initialize() {
#ifdef HAVE_OPDI
#if !defined(HAVE_OMPT)
  opdi::backend = new opdi::MacroBackend;
  opdi::backend->init();
#endif
  opdi::logic = new opdi::OmpLogic;
  opdi::logic->init();
  opdi::tool = new CoDiOpDiLibTool<su2double>;
  opdi::tool->init();
#endif
}

void omp_finalize() {
#ifdef HAVE_OPDI
  opdi::tool->finalize();
  opdi::logic->finalize();
  opdi::backend->finalize();
  delete opdi::tool;
  delete opdi::logic;
#if !defined(HAVE_OMPT)
  delete opdi::backend;
#endif
#endif
}

#ifdef HAVE_OPDI
#include "opdi.cpp"
#endif
