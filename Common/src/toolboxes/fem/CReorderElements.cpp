/*!
 * \file CReorderElements.cpp
 * \brief Helper class used to reorder the owned elements
 *        after the partitioning
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

#include "../../../include/toolboxes/fem/CReorderElements.hpp" 

CReorderElements::CReorderElements(const unsigned long  val_GlobalElemID,
                                   const unsigned short val_TimeLevel,
                                   const bool           val_CommSolution,
                                   const unsigned short val_VTK_Type,
                                   const unsigned short val_nPolySol,
                                   const bool           val_JacConstant) {

  /* Copy the global elment ID, time level and whether or not this element
     must be communicated. */
  globalElemID = val_GlobalElemID;
  timeLevel    = val_TimeLevel;
  commSolution = val_CommSolution;

  /* Create the element type used in this class, which stores information of
     the VTK type, the polynomial degree of the solution and whether or not the
     Jacobian of the transformation is constant. As it is possible that the
     polynomial degree of the solution is zero, this convention is different
     from the convention used in the SU2 grid file. */
  elemType = val_VTK_Type + 100*val_nPolySol;
  if( !val_JacConstant ) elemType += 50;
}

bool CReorderElements::operator< (const CReorderElements &other) const {

  /* Elements with the lowest time level are stored first. */
  if(timeLevel != other.timeLevel) return timeLevel < other.timeLevel;

  /* Next comparison is whether or not the element must communicate its
     solution data to other ranks. Elements which do not need to do this
     are stored first. */
  if(commSolution != other.commSolution) return other.commSolution;

  /* The comparisons below, element type and global ID are there to
     give a unique distribution. However, these are not necessary
     for a correct functioning of the code, in contrast to the
     two first criteria. */
  if(elemType != other.elemType) return elemType < other.elemType;
  return globalElemID < other.globalElemID;
}
