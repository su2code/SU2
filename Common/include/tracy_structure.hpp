/*!
 * \file tracy_structure.hpp
 * \brief Tracy profiler interface header.
 *        Provides compatibility macros if the code is built without
 *        tracy support. Profiling macros are defined here so that
 *        they can be completely "disabled" when compiling without tracy.
 * \note Do not include tracy headers explicitly anywhere, use this header instead.
 * \note To enable tracy, define the TRACY_ENABLE macro during compilation.
 * \author Divyaprakash
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#ifdef HAVE_TRACY
#include "tracy/Tracy.hpp"
#define SU2_ZONE_SCOPED ZoneScoped
#define SU2_ZONE_SCOPED_N(name) ZoneScopedN(name)
#else
#define SU2_ZONE_SCOPED
#define SU2_ZONE_SCOPED_N(name)
#endif
