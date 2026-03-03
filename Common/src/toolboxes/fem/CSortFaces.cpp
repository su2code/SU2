/*!
 * \file CSortFaces.cpp
 * \brief Functor, used for a different sorting of the surface elements
 *        than the < operator of CFaceOfElement.
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

#include "../../../include/toolboxes/fem/CSortFaces.hpp"
#include "../../../include/fem/fem_geometry_structure.hpp"

bool CSortFaces::operator()(const CFaceOfElement& f0, const CFaceOfElement& f1) {
  /*--- Comparison in case both faces are boundary faces. ---*/
  if (f0.faceIndicator >= 0 && f1.faceIndicator >= 0) {
    /* Both faces are boundary faces. The first comparison is the boundary
       marker, which is stored in faceIndicator. */
    if (f0.faceIndicator != f1.faceIndicator) return f0.faceIndicator < f1.faceIndicator;

    /* Both faces belong to the same boundary marker. The second comparison is
       based on the on the local volume ID's of the adjacent elements. As the
       volumes are sorted according to the time levels for time accurate
       local stepping, there is no need to do this check seperately here. */
    unsigned long ind0 = f0.elemID0 < nVolElemTot ? f0.elemID0 : f0.elemID1;
    unsigned long ind1 = f1.elemID0 < nVolElemTot ? f1.elemID0 : f1.elemID1;

    return ind0 < ind1;
  }

  /*--- Comparison in case both faces are internal faces. ---*/
  if (f0.faceIndicator == -1 && f1.faceIndicator == -1) {
    /* Both faces are internal faces. First determine the minimum and maximum
       ID of its adjacent elements.  */
    unsigned long elemIDMin0 = min(f0.elemID0, f0.elemID1);
    unsigned long elemIDMax0 = max(f0.elemID0, f0.elemID1);

    unsigned long elemIDMin1 = min(f1.elemID0, f1.elemID1);
    unsigned long elemIDMax1 = max(f1.elemID0, f1.elemID1);

    /* Determine the situation. */
    if (elemIDMax0 < nVolElemTot && elemIDMax1 < nVolElemTot) {
      /* Both faces are matching internal faces. Determine whether or not these
         faces are local faces, i.e. faces between locally owned elements. */
      const bool face0IsLocal = elemIDMax0 < nVolElemOwned;
      const bool face1IsLocal = elemIDMax1 < nVolElemOwned;

      /* Check if both faces have the same status, i.e. either local or
         not local. */
      if (face0IsLocal == face1IsLocal) {
        /* Both faces are either local or not local. Determine the time level
           of the faces, which is the minimum value of the adjacent volume
           elements. */
        const unsigned short timeLevel0 = min(volElem[elemIDMin0].timeLevel, volElem[elemIDMax0].timeLevel);
        const unsigned short timeLevel1 = min(volElem[elemIDMin1].timeLevel, volElem[elemIDMax1].timeLevel);

        /* Internal faces with the same status are first sorted according to
           their time level. Faces with the smallest time level are numbered
           first. Note this is only relevant for time accurate local time
           stepping. */
        if (timeLevel0 != timeLevel1) return timeLevel0 < timeLevel1;

        /* The faces belong to the same time level. They are sorted according
           to their element ID's in order to increase cache performance. */
        if (elemIDMin0 != elemIDMin1) return elemIDMin0 < elemIDMin1;
        return elemIDMax0 < elemIDMax1;
      } /* One face is a local face and the other is not. Make sure that
           the local faces are numbered first. */
      return face0IsLocal;

    } else if (elemIDMax0 >= nVolElemTot && elemIDMax1 >= nVolElemTot) {
      /* Both faces are non-matching internal faces. Sort them according to
         their relevant element ID. The time level is not taken into account
         yet, because non-matching faces are not possible at the moment with
         time accurate local time stepping. */
      return elemIDMin0 < elemIDMin1;
    } else {
      /* One face is a matching internal face and the other face is a
         non-matching internal face. Make sure that the non-matching face
         is numbered after the matching face. This is accomplished by comparing
         the maximum element ID's. */
      return elemIDMax0 < elemIDMax1;
    }
  }

  /*--- One face is a boundary face and the other face is an internal face.
        Make sure that the boundary face is numbered first. This can be
        accomplished by using the > operator for faceIndicator. ---*/
  return f0.faceIndicator > f1.faceIndicator;
}
