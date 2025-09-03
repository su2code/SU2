/*!
 * \file SwapBytes.cpp
 * \brief Function to swap bytes of primitive data types
 * \author P. Gomes
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

#include "../../include/toolboxes/SwapBytes.hpp"

/*--- Function to swap bytes, in case we need to convert between
  big and little endian storage. ---*/
void SwapBytes(char* buffer, size_t nBytes, unsigned long nVar) {
  /*--- Store half the number of bytes in kk. ---*/
  const int kk = (int)nBytes / 2;

  /*--- Loop over the number of variables in the buffer. ---*/
  for (unsigned long j = 0; j < nVar; j++) {
    /*--- Initialize ii and jj, which are used to store the
     indices of the bytes to be swapped. ---*/
    unsigned long ii = j * nBytes;
    unsigned long jj = ii + nBytes - 1;

    /*--- Swap the bytes. ---*/
    for (int i = 0; i < kk; i++) {
      char tmp = buffer[jj];
      buffer[jj] = buffer[ii];
      buffer[ii] = tmp;

      ii++;
      jj--;
    }
  }
}
