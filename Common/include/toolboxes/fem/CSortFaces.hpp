/*!
 * \file CSortFaces.hpp
 * \brief Header file for the class CSortFaces.
 *        The implementations are in the <i>CSortFaces.cpp</i> file.
 * \author E. van der Weide
 * \version 8.2.0 "Harrier"
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

#pragma once

#include "CFaceOfElement.hpp"

using namespace std;

/*!
 * \class CSortFaces
 * \brief Functor, used for a different sorting of the faces than the < operator
 *        of CFaceOfElement.
 * \author E. van der Weide
 * \version 8.2.0 "Harrier"
 */
class CVolumeElementFEM;  // Forward declaration to avoid problems.
class CSortFaces {
 private:
  unsigned long nVolElemOwned; /*!< \brief Number of locally owned volume elements. */
  unsigned long nVolElemTot;   /*!< \brief Total number of local volume elements . */

  const CVolumeElementFEM* volElem; /*!< \brief The locally stored volume elements. */

 public:
  /*!
   * \brief Constructor of the class. Set the values of the member variables.
   */
  CSortFaces(unsigned long val_nVolElemOwned, unsigned long val_nVolElemTot, const CVolumeElementFEM* val_volElem) {
    nVolElemOwned = val_nVolElemOwned;
    nVolElemTot = val_nVolElemTot;
    volElem = val_volElem;
  }

  /*!
   * \brief Default constructor of the class. Disabled.
   */
  CSortFaces() = delete;

  /*!
   * \brief Operator used for the comparison.
   * \param[in] f0 - First face in the comparison.
   * \param[in] f1 - Second face in the comparison.
   */
  bool operator()(const CFaceOfElement& f0, const CFaceOfElement& f1);
};
