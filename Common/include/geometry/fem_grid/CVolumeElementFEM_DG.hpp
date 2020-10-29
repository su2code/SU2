/*!
 * \file CVolumeElementFEM_DG.hpp
 * \brief Class for a volume element for the DG-FEM solver.
 *        The implementations are in the <i>CVolumeElementFEM_DG.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
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

#include "CVolumeElementFEM_Base.hpp"

using namespace std;

/*!
 * \class CVolumeElementFEM_DG
 * \brief Class to store a volume element for the DG-FEM solver.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
 */
class CVolumeElementFEM_DG final: public CVolumeElementFEM_Base {
public:

  short periodIndexToDonor;    /*!< \brief The index of the periodic transformation to the donor
                                           element. Only for halo elements. A -1 indicates no
                                           periodic transformation. */
  unsigned short timeLevel;    /*!< \brief Time level of the element when time accurate local
                                           time stepping is employed. */

  unsigned int factTimeLevel;  /*!< \brief Number of local time steps for this element
                                           compared to the largest time step when time
                                           accurate local time stepping is employed. */

  unsigned long offsetDOFsSolGlobal; /*!< \brief Global offset of the solution DOFs of this element. */

  CFEMStandardElementBase *standardElemSol; /*!< \brief Pointer to the standard element for the solution. */
};
